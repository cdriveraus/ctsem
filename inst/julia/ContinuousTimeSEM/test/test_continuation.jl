using ForwardDiff, LinearAlgebra

# The quadrature continuation (`laplace_continuation.jl`): a hybrid objective
# whose flagged units carry fixed-node quadrature rules and whose other units
# keep their Laplace terms, with the Fisher-identity gradient.
#
# What has to hold, each checked on its own:
#
#  1. The gradient is the gradient of the value returned -- central
#     differences of the fixed-node value, elementwise, at points away from the
#     centre, on a nonlinear one-level fixture and on units that nest.
#  2. At its centre the fixed rule IS the adaptive quadrature: the same nodes,
#     so the same value as `ctsem_laplace_quadrature`, unit by unit.
#  3. Where nothing is flagged the hybrid IS the Laplace objective, value and
#     gradient, bit for bit; and on a Gaussian integrand the soft-direction
#     rule is exact, as every rule of this family must be.
#  4. Re-placing the nodes and putting them back restores the objective.
#  5. A round never leaves its trust region.
#  6. At its centre, on a Gaussian integrand, the soft rule's fixed-node
#     gradient and Hessian are the marginal's -- which a stiff complement
#     taken at one node gets wrong, although its value is exact.
#  7. The stiff complement's rule has a standard normal's moments through
#     degree four, cross moments included.
#  8. A unit wider than `maxdim` is left on Laplace, unscored.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

"""Central differences of the continuation's value, against its gradient."""
function _continuation_fd_check(o, x; h=1e-5)
    g = ctsem_laplace_continuation_evaluate(o, x; gradient=true).gradient
    fd = similar(g)
    for j in eachindex(x)
        hj = h * max(1.0, abs(x[j]))
        xp = copy(x); xp[j] += hj
        xm = copy(x); xm[j] -= hj
        vp = ctsem_laplace_continuation_evaluate(o, xp; gradient=false).value
        vm = ctsem_laplace_continuation_evaluate(o, xm; gradient=false).value
        fd[j] = (vp - vm) / (2hj)
    end
    return (g=g, fd=fd)
end

@testset "the fixed-node gradient is the derivative of the fixed-node value" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("two levels", _fresh_twolevel),
                           ("three levels", _fresh_threelevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        # Tolerance zero flags every unit with a gap, so the rule is what is
        # being differentiated wherever there is one to differentiate.
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0)
        n = length(theta)
        # Away from the centre, where the nodes are no longer the adaptive
        # rule's and only the fixed-node derivative is exact.
        for shift in (zeros(n), fill(0.05, n), [0.1 * (-1)^j for j in 1:n])
            x = theta .+ shift
            r = _continuation_fd_check(o, x)
            @test (label, all(isfinite, r.g)) == (label, true)
            for j in eachindex(x)
                @test (label, j, isapprox(r.g[j], r.fd[j]; atol=1e-6, rtol=1e-5)) ==
                    (label, j, true)
            end
        end
    end
end

@testset "at its centre the fixed rule is the adaptive quadrature" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("three levels", _fresh_threelevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0)
        info = ctsem_laplace_continuation_info(o)
        reference = ctsem_laplace_quadrature(laplace, theta; nodes=5,
            contributions=true)
        @test (label, isapprox(info.quadrature, reference.value; rtol=1e-12)) ==
            (label, true)
        @test (label, maximum(abs, info.quadrature_units .- reference.subject_loglik)) <
            (label, 1e-9)
        # Every unit with a gap is flagged, and the hybrid at the centre is then
        # the quadrature value itself.
        @test (label, isapprox(ctsem_laplace_continuation_evaluate(o, theta;
            gradient=false).value, reference.value; rtol=1e-12)) == (label, true)
    end
end

@testset "with nothing flagged the hybrid is the Laplace objective" begin
    laplace, values = _fresh_linear()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta)
    info = ctsem_laplace_continuation_info(o)
    # A Gaussian integrand: the rule and Laplace agree to rounding, the screen
    # passes, and nothing is flagged.
    @test info.nflagged == 0
    @test info.screen < 1e-8
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=true)
    hyb = ctsem_laplace_continuation_evaluate(o, theta; gradient=true)
    @test hyb.value == lap.value
    # Gradients are summed per worker, so a threaded run may reassociate.
    @test isapprox(hyb.gradient, lap.gradient; rtol=1e-12, atol=1e-14)
    # And the nonlinear fixture with the tolerance out of reach: nothing is
    # flagged there either, and the same identity holds away from the centre.
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=Inf)
    @test ctsem_laplace_continuation_info(o).nflagged == 0
    x = theta .+ 0.05
    lap = ctsem_laplace_evaluate(laplace, x; gradient=true)
    hyb = ctsem_laplace_continuation_evaluate(o, x; gradient=true)
    @test hyb.value == lap.value
    @test isapprox(hyb.gradient, lap.gradient; rtol=1e-12, atol=1e-14)
end

@testset "the soft-direction rule is exact on a Gaussian integrand" begin
    # `product_maxdim = 0` puts every block on the soft rule, including the
    # outer ones of the nested fixture; each integrates a Gaussian exactly, so
    # the rule's value is the Laplace value, which is exact there.
    for (label, fresh) in (("one level", _fresh_linear),
                           ("two levels", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        reference = ctsem_laplace_evaluate(laplace, theta; gradient=false)
        for maxdirs in (1, 2)
            o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0,
                product_maxdim=0, soft_maxdirs=maxdirs)
            info = ctsem_laplace_continuation_info(o)
            @test (label, maxdirs, info.soft_blocks > 0 || info.nflagged == 0) ==
                (label, maxdirs, true)
            @test (label, maxdirs, isapprox(info.quadrature, reference.value;
                atol=1e-7)) == (label, maxdirs, true)
        end
    end
    # On the nonlinear fixture the soft rule is a genuine approximation, and
    # its gradient must still be the derivative of its value.
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, product_maxdim=0)
    r = _continuation_fd_check(o, theta .+ 0.05)
    @test isapprox(r.g, r.fd; atol=1e-6, rtol=1e-5)
end

@testset "re-placing the nodes and putting them back restores the objective" begin
    laplace, values = _fresh_threelevel()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0)
    x = theta .+ 0.03
    before = ctsem_laplace_continuation_evaluate(o, x; gradient=true)
    moved = ctsem_laplace_continuation_recentre!(o, x)
    @test moved.centre == x
    after = ctsem_laplace_continuation_evaluate(o, x; gradient=false)
    # The nodes moved, so the value at `x` is now the adaptive rule's there.
    @test isapprox(after.value, ctsem_laplace_quadrature(laplace, x; nodes=5).value;
        rtol=1e-12)
    back = ctsem_laplace_continuation_revert!(o)
    @test back.centre == theta
    again = ctsem_laplace_continuation_evaluate(o, x; gradient=true)
    @test again.value == before.value
    @test isapprox(again.gradient, before.gradient; rtol=1e-12, atol=1e-14)
end

@testset "a round stays inside its trust region and does not lose objective" begin
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    fit = ctsem_laplace_optimize(laplace, theta; maxiter=200, tune_chunks=false)
    est = fit.minimizer
    H = ctsem_laplace_hessian(laplace, est)
    E = eigen(Symmetric(-(H .+ transpose(H)) ./ 2))
    keep = abs.(E.values) .> 1e-8 * maximum(abs, E.values)
    B = E.vectors[:, keep] * Diagonal(1 ./ sqrt.(abs.(E.values[keep])))
    o = ctsem_laplace_continuation(laplace, est; tolerance=0.0)
    start = ctsem_laplace_continuation_evaluate(o, est; gradient=false).value
    for radius in (0.05, 0.5)
        r = ctsem_laplace_continuation_optimize(o, est, B, radius; maxiter=50)
        if r.stationary
            @test r.minimizer == est
            continue
        end
        y = B \ (r.minimizer .- est)
        @test norm(y) <= radius * (1 + 1e-9)
        @test r.value >= start - 1e-10
        @test r.value == ctsem_laplace_continuation_evaluate(o, r.minimizer;
            gradient=false).value
    end
    # `stationary_only` evaluates and reports, and moves nothing.
    s = ctsem_laplace_continuation_optimize(o, est, B, 0.5; stationary_only=true)
    @test s.minimizer == est
    @test s.iterations == 0
    @test isfinite(s.start_gain)
end

@testset "at its centre the soft rule's gradient and curvature are the marginal's" begin
    # With one soft direction, the two-effect subject blocks keep a stiff
    # complement of one. Its value at the centre is exact on a Gaussian
    # integrand whatever rule it gets; its derivatives are not. Held at a
    # single node, the complement's log determinant drops out of the gradient
    # and its covariance out of the curvature -- on the gated-gaps A14 config
    # that gave the continuation standard errors 0.07 to 0.45 of Laplace's.
    # Laplace is exact here, so its gradient and Hessian are the reference.
    for (label, fresh) in (("one level", _fresh_linear),
                           ("two levels", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0,
            product_maxdim=0, soft_maxdirs=1)
        # A unit whose rule matches Laplace to the last bit has no gap and is
        # not flagged; the rest carry the rule.
        info = ctsem_laplace_continuation_info(o)
        @test (label, info.nflagged > 0, info.soft_blocks > 0) == (label, true, true)
        lap = ctsem_laplace_evaluate(laplace, theta; gradient=true)
        hyb = ctsem_laplace_continuation_evaluate(o, theta; gradient=true)
        @test (label, isapprox(hyb.value, lap.value; atol=1e-7)) == (label, true)
        for j in eachindex(theta)
            @test (label, j, isapprox(hyb.gradient[j], lap.gradient[j];
                rtol=1e-6, atol=1e-8)) == (label, j, true)
        end
        # Central differences, as the Laplace Hessian is taken, match it to the
        # finite-difference error. The default forward scheme matches those to
        # its own truncation, which is O(step) in absolute terms -- the third
        # derivative times the step, whatever the size of the entry -- so it is
        # held to a share of the matrix's scale rather than of each entry's:
        # a small cross term can differ by more than 0.5% of itself, as one
        # here does, where nothing a covariance reads moves.
        H = ctsem_laplace_continuation_hessian(o, theta; scheme=:central)
        Hl = ctsem_laplace_hessian(laplace, theta)
        for i in eachindex(theta), j in eachindex(theta)
            @test (label, i, j, isapprox(H[i, j], Hl[i, j]; rtol=1e-4, atol=1e-6)) ==
                (label, i, j, true)
        end
        Hf = ctsem_laplace_continuation_hessian(o, theta)
        scale = maximum(abs, H)
        for i in eachindex(theta), j in eachindex(theta)
            @test isapprox(Hf[i, j], H[i, j]; atol=1e-3 * scale)
        end
        @test_throws ArgumentError ctsem_laplace_continuation_hessian(o, theta;
            scheme=:backward)
    end
end

@testset "the stiff complement's rule integrates the moments a derivative needs" begin
    # The unscented transform's points are exact to degree three, which is
    # every Gaussian integrand, but not for the cross moments of a
    # two-dimensional complement; on the gated-gaps A14 config that left the
    # fixed-node derivative 100 times the re-placed rule's. The product rule
    # that replaced it: weights summing to one, and the moments through degree
    # four of a standard normal, cross ones included.
    rule = ContinuousTimeSEM._continuation_stiff_rule
    @test rule(0, 5).logweights == [0.0]
    for nh in (1, 2, 3)
        r = rule(nh, 5)
        w = exp.(r.logweights)
        X = r.points
        @test (nh, size(X)) == (nh, (nh, 5^nh))
        @test (nh, isapprox(sum(w), 1.0; atol=1e-12)) == (nh, true)
        moment(f) = sum(w[c] * f(X[:, c]) for c in eachindex(w))
        for i in 1:nh
            @test (nh, i, isapprox(moment(x -> x[i]), 0.0; atol=1e-12)) == (nh, i, true)
            @test (nh, i, isapprox(moment(x -> x[i]^2), 1.0; atol=1e-12)) == (nh, i, true)
            @test (nh, i, isapprox(moment(x -> x[i]^4), 3.0; atol=1e-10)) == (nh, i, true)
            for j in (i + 1):nh
                @test (nh, i, j, isapprox(moment(x -> x[i] * x[j]), 0.0; atol=1e-12)) ==
                    (nh, i, j, true)
                @test (nh, i, j, isapprox(moment(x -> x[i]^2 * x[j]^2), 1.0;
                    atol=1e-10)) == (nh, i, j, true)
            end
        end
    end
end

@testset "a unit wider than maxdim is not scored and keeps its Laplace term" begin
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    # Each unit is a study (one effect) over two subjects (two effects each):
    # three effects on every root-to-leaf path of its block tree.
    @test ContinuousTimeSEM._continuation_unit_width(laplace, 1) == 3
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, maxdim=2)
    info = ctsem_laplace_continuation_info(o)
    @test (info.nwide, info.nflagged, info.rule_failures) == (info.nunits, 0, 0)
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=true)
    hyb = ctsem_laplace_continuation_evaluate(o, theta; gradient=true)
    @test hyb.value == lap.value
    @test isapprox(hyb.gradient, lap.gradient; rtol=1e-12, atol=1e-14)
    @test ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace, theta;
        tolerance=0.0, maxdim=3)).nwide == 0
end
