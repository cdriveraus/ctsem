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
        # finite-difference error. The exact scheme, the default, matches those
        # entry by entry to their own error. The forward scheme matches them to
        # its truncation, which is O(step) in absolute terms -- the third
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
        He = ctsem_laplace_continuation_hessian(o, theta)
        scale = maximum(abs, H)
        for i in eachindex(theta), j in eachindex(theta)
            @test (label, i, j, isapprox(He[i, j], H[i, j]; rtol=1e-6,
                atol=1e-7 * scale)) == (label, i, j, true)
        end
        Hf = ctsem_laplace_continuation_hessian(o, theta; scheme=:forward)
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

# --- two states, for the exact Hessian ------------------------------------------
#
# The exact Hessian differentiates the reverse pass, and one latent hides whole
# classes of reverse-pass error: a covariance cotangent that is a scalar is
# symmetric whatever the code assumes. So two latents with coupled dynamics, a
# random effect on a transformed DRIFT -- the integrand is not Gaussian, so
# every unit carries a gap and is flagged -- and random CINT and T0MEANS on
# different latents.
#
# Raw layout: 1-6 the model parameters, then the population scales, then the
# correlations.
if !isdefined(@__MODULE__, :_CONTINUATION_2D_OBJECTIVE)
const _CONTINUATION_2D_OBJECTIVE = let
    df = _laplace_test_dataframe(
        drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3],
        cint=[0.0; 0.0;;], diffusion=[0.3 0.0; 0.05 0.25],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.2 0.0; 0.0 0.2],
        t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DRIFT, 1, 2) => (2, "param[2]"),
            (:JAx, 1, 2) => (2, "param[2]"),
            (:CINT, 1, 1) => (3, "param[3]"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:T0MEANS, 2, 1) => (5, "param[5]"),
            (:MANIFESTMEANS, 2, 1) => (6, "param[6]"),
        ),
    )
    starts, times, data = _laplace_test_data(6, 2; seed=7)
    ctsem_objective(ekf_from_data_frame(df), starts, times, data)
end
end

# Three effects a subject (DRIFT[1,1], CINT[1], T0MEANS[2]), so wider than the
# product rule: scales 7-9, correlations 10-12.
_fresh_2d3() = (ctsem_laplace_objective(_CONTINUATION_2D_OBJECTIVE, [1, 3, 5],
    [7, 8, 9], [10, 11, 12], [1.0, 1.0, 1.0]),
    [0.3, 0.2, 0.1, -0.5, 0.2, -0.1, 0.4, 0.3, 0.5, 0.2, -0.1, 0.15])
# Two (DRIFT[1,1], T0MEANS[2]), on the product rule: scales 7-8, correlation 9.
_fresh_2d2() = (ctsem_laplace_objective(_CONTINUATION_2D_OBJECTIVE, [1, 5],
    [7, 8], [9], [1.0, 1.0]), [0.3, 0.2, 0.1, -0.5, 0.2, -0.1, 0.4, 0.5, 0.2])

@testset "a member's exact Hessian at a node is forward mode's, cross terms included" begin
    # The assembly `J' Hs J` plus the second derivative of `L u`, against
    # forward-over-forward through the primal filter and the Cholesky factor,
    # which shares none of that code. Every entry: the blocks between the
    # population parameters and the rest are where an assembly error would
    # sit, and they are not small here.
    CT = ContinuousTimeSEM
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    npar = length(theta)
    spec = laplace.spec
    Ls = CT._laplace_popchols(theta, spec)
    dL = CT._laplace_level_chol_derivatives(theta, spec)
    positions = Vector{Int}[CT._laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    layout = CT._continuation_hessian_layout(spec, npar, positions, dL)
    d2L = CT._continuation_level_chol_second_derivatives(theta, spec, positions)
    @test layout.active == collect(1:6)
    @test layout.pop == collect(7:12)
    CT._laplace_ensure_pool!(laplace)
    U, m = 2, 1
    u = [0.7, -1.1, 0.4]
    offsets = laplace.units.offsets[U][m]
    subject = laplace.objective.subject_objectives[laplace.units.members[U][m]]
    f = y -> subject(CT._laplace_member_values(y, spec, CT._laplace_popchols(y, spec),
        u, offsets))
    reference = ForwardDiff.hessian(f, theta,
        ForwardDiff.HessianConfig(f, theta, ForwardDiff.Chunk{4}()))
    gradient = ForwardDiff.gradient(f, theta)
    scale = maximum(abs, reference)
    pop, act = layout.pop, layout.active
    @test maximum(abs, reference[pop, act]) > 1e-2 * scale
    @test maximum(abs, reference[pop, pop]) > 1e-2 * scale
    for width in (1, 2, 6)
        S = ForwardDiff.Dual{CT._LaplaceSeedInner,Float64,width}
        aws = CT._laplace_workspace!(laplace, S, npar)
        sc = CT._continuation_hessian_scratch(laplace, S, npar, layout, 1)
        g = zeros(npar)
        H = zeros(npar, npar)
        ll, status = CT._continuation_member_hessian!(g, H, laplace, U, m, theta, Ls,
            dL, d2L, layout, u, aws, sc)
        @test (width, status) == (width, 0)
        @test (width, isapprox(ll, f(theta); rtol=1e-12)) == (width, true)
        for i in 1:npar
            @test (width, i, isapprox(g[i], gradient[i]; rtol=1e-10,
                atol=1e-12 * maximum(abs, gradient))) == (width, i, true)
        end
        for i in 1:npar, j in 1:npar
            @test (width, i, j, isapprox(H[i, j], reference[i, j]; rtol=1e-9,
                atol=1e-11 * scale)) == (width, i, j, true)
        end
    end
end

@testset "the exact Hessian is the derivative of the hybrid's gradient" begin
    # Two states, on the soft and the product rule; and units that nest, where
    # an outer block's node carries its children's rules and the identity is
    # applied at each level of the tree.
    for (label, fresh) in (("three effects, soft rule", _fresh_2d3),
                           ("two effects, product rule", _fresh_2d2),
                           ("two levels", _fresh_twolevel),
                           ("three levels", _fresh_threelevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0)
        info = ctsem_laplace_continuation_info(o)
        # Every unit on its rule, so what is differenced is the prior alone and
        # the reference below measures the exact part.
        @test (label, info.nflagged) == (label, info.nunits)
        x = theta .+ [0.02 * (-1)^j for j in eachindex(theta)]
        He = ctsem_laplace_continuation_hessian(o, x)
        @test (label, all(isfinite, He)) == (label, true)
        # Its byproducts are the gradient path's.
        flagged = ContinuousTimeSEM._continuation_flagged_hessian(o, x)
        plain = ContinuousTimeSEM._continuation_flagged(o, x, true)
        @test (label, isapprox(flagged.values, plain.values; rtol=1e-12)) == (label, true)
        @test (label, isapprox(flagged.gradient, plain.gradient; rtol=1e-11)) ==
            (label, true)
        # Central differences at two steps, extrapolated: O(step^4).
        Hc = ctsem_laplace_continuation_hessian(o, x; scheme=:central, step=1e-4)
        Hr = (4 .* Hc .- ctsem_laplace_continuation_hessian(o, x; scheme=:central,
            step=2e-4)) ./ 3
        scale = maximum(abs, He)
        for i in eachindex(x), j in eachindex(x)
            @test (label, i, j, isapprox(He[i, j], Hr[i, j]; rtol=1e-7,
                atol=1e-9 * scale)) == (label, i, j, true)
        end
        # A sweep's width changes nothing but rounding: a direction a sweep,
        # two (several sweeps a node, the last one short), all at once.
        for width in (1, 2, 5, 6)
            Hw = ctsem_laplace_continuation_hessian(o, x; width=width)
            @test (label, width, maximum(abs, Hw .- He) <= 1e-11 * scale) ==
                (label, width, true)
        end
    end
end
