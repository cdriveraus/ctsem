using ForwardDiff, LinearAlgebra, Statistics

# The placed objective (`laplace_placed.jl`): the hybrid with each flagged
# unit's rule -- product Gauss-Hermite or lattice-t (`point_sets.jl`) -- placed
# at its mode and curvature at every evaluation, and its exact gradient
# through the placement.
#
# What has to hold, each checked on its own:
#
#  1. At the continuation's centre the placed value is the continuation's
#     placed value, unit by unit: the same rule at the same point.
#  2. The gradient is the derivative of the value: central differences of the
#     placed value, at the centre and away from it, on a nonlinear one-level
#     fixture, a nonlinear two-level one (a leaf's conditional mode chained
#     through the root's nodes), a linear two-level one, and on lattice
#     rules: four effects a subject, and a four-effect study over subjects.
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

# Four effects a subject on the nonlinear objective: T0MEANS, the
# state-dependent drift, DIFFUSION and MANIFESTMEANS. Raw layout: 1-4 the model
# parameters, 5-8 the scales, 9-14 the correlations.
_fresh_placed_wide() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 2, 3, 4], sd_index=[5, 6, 7, 8], cor_index=collect(9:14),
    sd_scale=[1.0, 1.0, 1.0, 1.0]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.6, -0.6, -0.4, 0.1, 0.0, -0.1, 0.05, 0.0, 0.1])

# Every subject's T0MEANS under a study's four effects: each study a root of
# four effects, on the lattice, over leaves of one, on product rules. Raw
# layout: 1-4 the model parameters, 5 the subject scale, 6-9 the study scales,
# 10-15 their correlations.
_fresh_placed_wide_nested() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 1, 2, 3, 4], sd_index=[5, 6, 7, 8, 9], cor_index=collect(10:15),
    sd_scale=[1.0, 1.0, 1.0, 1.0, 1.0], level_nre=[1, 4],
    group=vcat(1:5, [1, 1, 2, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.5, -0.6, -0.6, -0.4, 0.1, 0.0, -0.1, 0.05,
     0.0, 0.1])

# Every unit flagged (a negative tolerance), each on its screened rule
# (`draws_se` out of reach, so the refine keeps it).
_placed_from(laplace, theta) = ctsem_laplace_placed(ctsem_laplace_continuation(
    laplace, theta; tolerance=-1.0, draws_se=1e9))

@testset "the placed value is the continuation's placed value at its centre" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("nested", _fresh_placed_nested),
                           ("two levels, linear", _fresh_twolevel),
                           ("lattice", _fresh_placed_wide),
                           ("lattice root", _fresh_placed_wide_nested))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws_se=1e9)
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
                           ("two levels, linear", _fresh_twolevel),
                           ("lattice", _fresh_placed_wide),
                           ("lattice root", _fresh_placed_wide_nested))
        laplace, values = fresh()
        theta = collect(Float64, values)
        P = _placed_from(laplace, theta)
        lattice = label in ("lattice", "lattice root")
        @test (label, lattice == any(U -> any(ps -> ps.family === :lattice,
            P.sets[U]), P.placed)) == (label, true)
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
    o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0)
    P = ctsem_laplace_placed(o)
    placed = ctsem_laplace_placed_evaluate(P, theta; gradient=true).gradient
    frozen = ctsem_laplace_continuation_evaluate(o, theta; gradient=true).gradient
    @test all(isapprox.(placed, frozen; rtol=1e-7, atol=1e-8))
end

@testset "a lattice rule past its target is raised, and the error is read where it sits" begin
    C = ContinuousTimeSEM
    laplace, values = _fresh_placed_wide()
    theta = collect(Float64, values)
    P = _placed_from(laplace, theta)
    ctsem_laplace_placed_evaluate(P, theta; gradient=false)
    U = P.placed[1]
    @test P.sets[U][1].family === :lattice && P.sets[U][1].count == C._LATTICE_START
    before = P.se[U]
    @test isfinite(before) && before > 0
    # Under a target the start's spread exceeds, the unit is rough and its
    # points double until the spread is within it, or the lattice's most.
    P.cont.draws_se = before / 3
    @test U in C._placed_rough(P, 1.0)
    @test C._placed_raise!(P, theta, [U]) == 1
    @test P.sets[U][1].count > C._LATTICE_START
    @test P.cont.unit_rules[U].points == P.sets[U][1].count
    @test P.se[U] < before
    @test P.se[U] <= P.cont.draws_se || P.sets[U][1].count == C._LATTICE_MAX
    # The error at the last evaluation: the lattice units' spread there in
    # quadrature.
    ctsem_laplace_placed_evaluate(P, theta .+ 0.05; gradient=false)
    e = ctsem_laplace_placed_error(P)
    lat = [V for V in P.placed if P.sets[V][1].family === :lattice]
    @test isapprox(e.error_draws, sqrt(sum(abs2, P.se[lat])); rtol=1e-12)
    @test e.centre_error[lat] == P.se[lat]
    @test isapprox(e.error, e.error_draws + e.error_nodes + e.error_residual; rtol=1e-12)
    # The slope check: the lattice units' gradients on two independent sets
    # of shifts, through the curvature, in standard errors. Every unit a
    # bar of zero leaves carrying something, none under an infinite one, and
    # the other sets are lent and taken back.
    H = ctsem_laplace_placed_hessian(P, theta)
    sets = copy(P.sets)
    loose = C._placed_slope_check(P, theta, H; bar=Inf)
    @test isfinite(loose.sd) && loose.sd > 0
    @test isempty(loose.units)
    tight = C._placed_slope_check(P, theta, H; bar=0.0)
    @test tight.sd == loose.sd
    @test !isempty(tight.units) && issubset(tight.units, lat)
    @test loose.doublings == 0 && 1 <= tight.doublings <= 5
    @test all(V -> P.sets[V] === sets[V], P.placed)
    # Without a curvature there is nothing to carry the error through.
    @test isempty(C._placed_slope_check(P, theta, zeros(1, 1)).units)
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
