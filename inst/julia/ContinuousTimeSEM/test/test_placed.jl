using ForwardDiff, LinearAlgebra, Statistics

# The placed objective (`laplace_placed.jl`): the hybrid with each flagged
# unit's product rule placed at its mode and curvature at every evaluation, and
# its exact gradient through the placement.
#
# What has to hold, each checked on its own:
#
#  1. At the continuation's centre the placed value is the continuation's
#     placed value, unit by unit: the same rule at the same point.
#  2. The gradient is the derivative of the value: central differences of the
#     placed value, at the centre and away from it, on a nonlinear one-level
#     fixture, a nonlinear two-level one (a leaf's conditional mode chained
#     through the root's nodes) and a linear two-level one.
#  3. On a Gaussian integrand the rule is exact wherever it sits, so the
#     placement terms vanish and the gradient is the frozen rule's.
#  4. The Hessian, localised to each column's readers, is the plain central
#     difference of the exact gradient.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

# Subjects nested in two studies on the nonlinear objective: a T0MEANS effect
# a subject, whose state-dependent drift makes the integrand non-Gaussian, and
# a MANIFESTMEANS effect a study. Raw layout: 1-4 the model parameters, 5 the
# subject scale, 6 the study scale.
_fresh_placed_nested() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 4], sd_index=[5, 6], cor_index=Int[], sd_scale=[1.0, 1.0],
    level_nre=[1, 1], group=vcat(1:5, [1, 1, 2, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.1])

# Every unit flagged (a negative tolerance), every rule a product rule.
_placed_from(laplace, theta) = ctsem_laplace_placed(ctsem_laplace_continuation(
    laplace, theta; tolerance=-1.0, product_maxdim=5))

@testset "the placed value is the continuation's placed value at its centre" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("nested", _fresh_placed_nested),
                           ("two levels, linear", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, product_maxdim=5)
        P = ctsem_laplace_placed(o)
        @test (label, length(P.placed)) == (label, length(o.flagged))
        r = ctsem_laplace_placed_evaluate(P, theta; gradient=false)
        @test (label, isapprox(r.value, o.centre_value; rtol=1e-10)) == (label, true)
        @test (label, maximum(abs, r.unit_loglik[o.flagged] .-
            o.centre_quadrature[o.flagged]) < 1e-9) == (label, true)
    end
end

@testset "the placed gradient is the derivative of the placed value" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("nested", _fresh_placed_nested),
                           ("two levels, linear", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        P = _placed_from(laplace, theta)
        n = length(theta)
        for shift in (zeros(n), fill(0.05, n), [0.1 * (-1)^j for j in 1:n])
            r = ctsem_laplace_placed_check(P, theta .+ shift; step=1e-6)
            @test (label, all(isfinite, r.gradient)) == (label, true)
            @test (label, all(isapprox.(r.gradient, r.differences; rtol=1e-5,
                atol=1e-6))) == (label, true)
        end
    end
end

@testset "on a Gaussian integrand the placement terms vanish" begin
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, product_maxdim=5)
    P = ctsem_laplace_placed(o)
    placed = ctsem_laplace_placed_evaluate(P, theta; gradient=true).gradient
    frozen = ctsem_laplace_continuation_evaluate(o, theta; gradient=true).gradient
    @test all(isapprox.(placed, frozen; rtol=1e-7, atol=1e-8))
end

@testset "the localised Hessian is the difference of the exact gradient" begin
    laplace, values = _fresh_placed_nested()
    theta = collect(Float64, values)
    P = _placed_from(laplace, theta)
    H = ctsem_laplace_placed_hessian(P, theta; step=1e-4)
    @test all(isfinite, H)
    @test H == transpose(H)
    home = deepcopy(P.anchors)
    n = length(theta)
    D = zeros(n, n)
    for t in 1:n
        h = 1e-4 * max(1.0, abs(theta[t]))
        xp = copy(theta); xp[t] += h
        xm = copy(theta); xm[t] -= h
        P.anchors = deepcopy(home); P.cache = nothing; P.value_cache = nothing
        gp = ctsem_laplace_placed_evaluate(P, xp; gradient=true).gradient
        P.anchors = deepcopy(home); P.cache = nothing; P.value_cache = nothing
        gm = ctsem_laplace_placed_evaluate(P, xm; gradient=true).gradient
        D[:, t] .= (gp .- gm) ./ (2h)
    end
    D = (D .+ transpose(D)) ./ 2
    @test all(isapprox.(H, D; rtol=1e-6, atol=1e-8))
end
