using ForwardDiff, LinearAlgebra, Statistics

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
#     gradient, bit for bit.
#  4. The point sets: the product rule integrates a normal's moments exactly,
#     the lattice-t rule integrates one within its stated error, and a
#     lattice rule placed at a unit integrates it within that error.
#  5. Every unit is scored, whatever its width, on the rule its cost chooses.
#  6. Under `draws_all`, one block or nested, a unit is integrated by frozen
#     draws.

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

@testset "the point sets integrate a normal as their family should" begin
    C = ContinuousTimeSEM
    # The product rule: exact through degree four, cross moments included.
    for k in (1, 2, 3)
        ps = C._point_set_product(k, 5)
        w = exp.(ps.logweights) .* exp.(-vec(sum(abs2, ps.points; dims=1)) ./ 2) ./
            (2 * pi)^(k / 2)
        X = ps.points
        moment(f) = sum(w[c] * f(X[:, c]) for c in eachindex(w))
        @test (k, isapprox(sum(w), 1.0; atol=1e-12)) == (k, true)
        for i in 1:k
            @test (k, i, isapprox(moment(x -> x[i]^2), 1.0; atol=1e-12)) == (k, i, true)
            @test (k, i, isapprox(moment(x -> x[i]^4), 3.0; atol=1e-10)) == (k, i, true)
            for j in (i + 1):k
                @test (k, i, j, isapprox(moment(x -> x[i]^2 * x[j]^2), 1.0;
                    atol=1e-10)) == (k, i, j, true)
            end
        end
    end
    # One too large to evaluate is refused, not built.
    @test_throws ArgumentError C._point_set_product(13, 5)
    # The t quantile in closed form is the general one, and its density
    # normalises.
    for p in (1e-9, 0.01, 0.3, 0.5, 0.77, 1 - 1e-6)
        x = ContinuousTimeSEM.SpecialFunctions.beta_inc_inv(2.0, 0.5,
            2 * min(p, 1 - p))[1]
        t = sqrt(4 * (1 - x) / x) * (p < 0.5 ? -1 : 1)
        @test (p, isapprox(C._t4_quantile(p), t; rtol=1e-8, atol=1e-12)) == (p, true)
    end
    @test isapprox(sum(exp(C._t4_logdensity(y)) for y in -400:0.01:400) * 0.01, 1.0;
        atol=1e-4)
    # The generating vector: odd, and the same each time it is built.
    z = C._lattice_cbc(5, 64)
    @test all(isodd, z) && z[1] == 1 && z == C._lattice_cbc(5, 64)
    # The lattice-t rule: the integral of a standard normal within a few of
    # its stated errors, at every size, the error falling with the points;
    # and a key fixes the points.
    for k in (2, 5, 9)
        last = Inf
        for n in (16, 128, 1024)
            ps = C._point_set_lattice(k, n, (1, k))
            @test size(ps.points) == (k, n * C._LATTICE_SHIFTS)
            t = ps.logweights .- vec(sum(abs2, ps.points; dims=1)) ./ 2 .-
                k * log(2 * pi) / 2
            peak = maximum(t)
            sums = zeros(ps.nrep)
            for j in eachindex(t)
                sums[ps.replicate[j]] += exp(t[j] - peak)
            end
            value = peak + log(sum(sums))
            se = C._point_set_error(sums, ps.nrep)
            @test (k, n, abs(value) < 5 * se + 1e-12) == (k, n, true)
            @test (k, n, se < last) == (k, n, true)
            last = se
        end
        @test C._point_set_lattice(k, 16, (1, k)).points ==
            C._point_set_lattice(k, 16, (1, k)).points
        @test C._point_set_lattice(k, 16, (1, k)).points !=
            C._point_set_lattice(k, 16, (2, k)).points
    end
end

@testset "a lattice rule placed at a unit integrates it within its error" begin
    C = ContinuousTimeSEM
    # Laplace is exact on the linear fixture, and the product rule at 9 nodes
    # is on the nonlinear three-effect one.
    for (label, fresh, nodes) in (("linear", _fresh_linear, 0), ("2d3", _fresh_2d3, 9))
        laplace, values = fresh()
        theta = C._laplace_extend(laplace, collect(Float64, values))
        Ls = C._laplace_popchols(theta, laplace.spec)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=1e6)
        opts = C._continuation_opts(o)
        ws = C._laplace_workspace!(laplace, Float64, length(theta))
        reference = nodes == 0 ?
            ctsem_laplace_evaluate(laplace, theta; gradient=false).unit_loglik :
            ctsem_laplace_quadrature(laplace, theta; nodes=nodes,
                contributions=true).subject_loglik
        for U in 1:length(laplace.units.members), n in (16, 256)
            r = C._continuation_unit_rule(laplace, U, theta, Ls, ws,
                merge(opts, (choice=C.CTSEMRuleChoice(:lattice, 5, n),)))
            @test (label, U, n, r.se > 0 && abs(r.value - reference[U]) <
                5 * r.se + 1e-9) == (label, U, n, true)
        end
    end
end

@testset "every unit is scored, whatever its width, on the rule its cost chooses" begin
    C = ContinuousTimeSEM
    # Each unit is a study (one effect) over two subjects (two effects each):
    # three effects on every root-to-leaf path, each block on a product rule.
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0)
    info = ctsem_laplace_continuation_info(o)
    @test info.rule_failures == 0
    @test all(isfinite, info.quadrature_units)
    @test all(U -> o.unit_rules[U].family === :mixed, 1:info.nunits)
    @test all(b -> C._rule_block_product(o.unit_rules[1], b.size),
        laplace.units.blocks[1])
    # Four effects a subject: screened on the lattice, whose 128 points are
    # fewer than the product rule's 625, and each refined to the cheapest
    # rule that reaches the target.
    laplace, values = _fresh_2d4()
    theta = collect(Float64, values)
    @test C._rule_default(laplace, 1, 5) == C.CTSEMRuleChoice(:lattice, 5, 16)
    @test C._rule_cost(laplace, 1, C._rule_default(laplace, 1, 5)) == 128
    o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0)
    info = ctsem_laplace_continuation_info(o)
    @test (info.rule_failures, info.nflagged) == (0, info.nunits)
    for U in o.flagged
        @test (U, o.importance[U] || info.centre_error[U] <= o.draws_se) == (U, true)
    end
    families = C._continuation_steps(laplace, 1, C._continuation_opts(o))
    @test [f[1].family for f in families] == [:product, :lattice]
    @test [c.nodes for c in families[1]] == [3, 5, 7, 9]
    @test [c.points for c in families[2]] == [16 * 2^i for i in 0:7]
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

# Three effects a subject (DRIFT[1,1], CINT[1], T0MEANS[2]): scales 7-9,
# correlations 10-12.
_fresh_2d3() = (ctsem_laplace_objective(_CONTINUATION_2D_OBJECTIVE, [1, 3, 5],
    [7, 8, 9], [10, 11, 12], [1.0, 1.0, 1.0]),
    [0.3, 0.2, 0.1, -0.5, 0.2, -0.1, 0.4, 0.3, 0.5, 0.2, -0.1, 0.15])

# Four (DRIFT[1,1], DRIFT[1,2], CINT[1], T0MEANS[2]): scales 7-10, correlations
# 11-16. Past three effects the lattice is the cheaper first rule.
_fresh_2d4() = (ctsem_laplace_objective(_CONTINUATION_2D_OBJECTIVE, [1, 2, 3, 5],
    [7, 8, 9, 10], collect(11:16), [1.0, 1.0, 1.0, 1.0]),
    [0.3, 0.2, 0.1, -0.5, 0.2, -0.1, 0.4, 0.2, 0.3, 0.5, 0.2, -0.1, 0.15, 0.0,
     0.1, -0.05])

# --- frozen importance draws -------------------------------------------------
#
# Under `draws_all` (the random-effects check), and where no rule converges, a
# unit gets draws held fixed in its standardised effects in place of nodes
# (`_continuation_importance_rule`).
# What has to hold: its value is the unit's integral to within its own error,
# its gradient is the derivative of its value, a gap within that error flags
# nothing, a unit the caller names is flagged whatever its gap, and with
# `screen = false` only that unit is scored.

@testset "a unit on frozen draws is valued as quadrature values it" begin
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    # Three effects a subject: the product rule's 125 nodes are the reference.
    reference = ctsem_laplace_quadrature(laplace, theta; nodes=5,
        contributions=true)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
        draws=256, draws_se=0.02, force=1:length(laplace.units.members))
    info = ctsem_laplace_continuation_info(o)
    @test (info.nimportance, info.nflagged) == (info.nunits, info.nunits)
    @test all(se -> 0 < se <= 0.02, info.centre_se)
    for U in 1:info.nunits
        @test abs(info.quadrature_units[U] - reference.subject_loglik[U]) <
            4 * info.centre_se[U] + 1e-3
    end
    # Away from the centre the draws are held, and the gradient is the
    # derivative of the value they give.
    r = _continuation_fd_check(o, theta .+ 0.05)
    @test isapprox(r.g, r.fd; atol=1e-6, rtol=1e-5)
    # Without `draws_all` a unit takes draws only where no rule converges:
    # here none.
    off = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace,
        theta; tolerance=0.0, draws=256))
    @test (off.nimportance, off.rule_failures) == (0, 0)
    @test_throws ArgumentError ctsem_laplace_continuation(laplace, theta;
        draws_all=true)
end

@testset "an importance rule on a Gaussian integrand is the Laplace value" begin
    laplace, values = _fresh_linear()
    theta = collect(Float64, values)
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=false)
    # On t components (the default) and on Gaussian ones alike. Laplace is
    # exact here, so every gap is the draws' error. That error's standard
    # error is the spread of eight replicate lattices, so a gap over it is a
    # t on 7 degrees of freedom: over 40 seeds of six units the mean squared
    # ratio was 1.3-1.9 (t7's 1.4), and over 60 the largest 4.8-6.2, under
    # Julia 1.12.7 and 1.13.3 alike (dev3, 2026-10-29). The bars sit past
    # t7's 1e-4 tail for each unit and five combined errors for their signed
    # sum, whose largest over 240 runs was 3.7. A bar of four for each unit
    # fails about one run in sixteen by t7, and failed 1.13.3's streams; the
    # screen at tolerance 0, which flags any unit past two errors, flagged in
    # 15-21 runs of 40. Its rules are tested exactly below.
    for df in (4, Inf)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
            draws=256, draws_df=df)
        info = ctsem_laplace_continuation_info(o)
        @test info.nimportance == info.nunits
        gap = info.quadrature_units .- lap.unit_loglik
        for U in 1:info.nunits
            @test abs(gap[U]) < 8 * info.centre_se[U] + 1e-8
        end
        @test abs(sum(gap)) < 5 * sqrt(sum(abs2, info.centre_se)) + 1e-8
    end
end

@testset "a t component's density is normalised and its draws stretch by a chi-square" begin
    C = ContinuousTimeSEM.CTSEMImportanceComponent
    dens = ContinuousTimeSEM._continuation_component_logdensity
    work = zeros(2)
    vecs = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    h = 0.1
    grid = -30:h:30
    for df in (4.0, Inf)
        c = C([0.2, -0.1], vecs, [2.0, 0.5], df)
        total = sum(exp(dens(c, [x, y], work)) for x in grid, y in grid) * h^2
        @test isapprox(total, 1.0; atol=1e-3)
    end
    # The stretch at the median normal is the median chi-square's on 4 degrees
    # of freedom (3.35669), and a Gaussian component is not stretched.
    @test isapprox(ContinuousTimeSEM._continuation_t_scale(4.0, 0.0),
        sqrt(4 / 3.35669); rtol=1e-4)
    @test ContinuousTimeSEM._continuation_t_scale(Inf, 1.7) == 1.0
    # A t component's draws are its Gaussian draws stretched: the same points
    # in direction, heavier in the tails. On the three-effect fixture both
    # integrate to the product rule's value within their error.
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    reference = ctsem_laplace_quadrature(laplace, theta; nodes=5,
        contributions=true)
    for df in (4, Inf)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
            draws=256, draws_se=0.02, draws_df=df,
            force=1:length(laplace.units.members))
        info = ctsem_laplace_continuation_info(o)
        for U in 1:info.nunits
            @test (df, U, abs(info.quadrature_units[U] - reference.subject_loglik[U]) <
                4 * info.centre_se[U] + 1e-3) == (df, U, true)
        end
    end
    @test_throws ArgumentError ctsem_laplace_continuation(laplace, theta;
        draws_all=true, draws=256, draws_df=0)
end

@testset "a unit named is flagged whatever its gap, and screen = false scores it alone" begin
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=false)
    o = ctsem_laplace_continuation(laplace, theta; draws_all=true, draws=256,
        force=[2], screen=false)
    info = ctsem_laplace_continuation_info(o)
    @test (info.nflagged, info.flagged) == (1, [2])
    others = setdiff(1:info.nunits, 2)
    @test info.quadrature_units[others] == info.laplace_units[others]
    hyb = ctsem_laplace_continuation_evaluate(o, theta; gradient=false)
    @test isapprox(hyb.value, lap.value - lap.unit_loglik[2] +
        info.quadrature_units[2]; rtol=1e-12)
    # A placement again keeps it, and a point given for it joins its proposal.
    C = ContinuousTimeSEM
    d = laplace.units.dims[2]
    C._continuation_add_points!(o, [2], reshape(fill(0.3, d), d, 1))
    C._continuation_place!(o, C._laplace_extend(laplace, theta .+ 0.01))
    @test ctsem_laplace_continuation_info(o).flagged == [2]
    @test length(o.points[2]) == 1
    @test_throws ArgumentError ctsem_laplace_continuation(laplace, theta;
        draws=256, force=[info.nunits + 1])
end

@testset "the screen sums small gaps of one sign, and each unit's own against its error" begin
    C = ContinuousTimeSEM
    # Ten gaps each under twice its error, sharing a sign: their sum, 0.04,
    # is resolved many times over against a combined error of 0.0095, so the
    # screen flags until what is left is within the tolerance plus twice that
    # error. The old per-unit rule flagged none (the validation's acnonlin).
    gaps = fill(0.004, 10)
    se = fill(0.003, 10)
    @test C._continuation_screen(gaps, se) ≈ 0.04 - 2 * sqrt(10) * 0.003
    # Four: with seven left the sum beyond twice its error is 0.012, with six
    # 0.009.
    flagged = C._continuation_flags(gaps, 0.01; se=se)
    @test length(flagged) == 4
    @test C._continuation_screen(gaps[setdiff(1:10, flagged)], se[setdiff(1:10, flagged)]) <= 0.01
    # Gaps that cancel are noise, whatever their size against their own error
    # says, unless one alone exceeds the tolerance beyond its error.
    @test isempty(C._continuation_flags([0.004, -0.004], 0.01; se=[0.003, 0.003]))
    @test C._continuation_flags([0.05, -0.05], 0.01; se=[0.003, 0.003]) == [1, 2]
    # With no error the rule is the product rules': the plain absolute sum.
    @test C._continuation_flags([0.3, 0.2, 0.005, 0.004], 0.01) == [1, 2]
    @test C._continuation_screen([0.3, -0.2], [0.0, 0.0]) ≈ 0.5
end

@testset "the importance rule's standard error is the spread of its replicates" begin
    # On the nonlinear three-effect fixture, the reported se against the spread
    # of the rule's value over seeds, unit by unit: the independent-draws
    # formula overstated it 3- to 10-fold at this width (the validation note);
    # the replicate spread should be within a factor of two of the truth, with
    # six seeds' own uncertainty allowed for.
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    seeds = 1:6
    vals = Float64[]; ses = Float64[]
    nunits = length(laplace.units.members)
    V = zeros(nunits, length(seeds)); S = zeros(nunits, length(seeds))
    for (i, seed) in enumerate(seeds)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
            draws=256, draws_max=256, draws_seed=seed,
            force=1:nunits)
        info = ctsem_laplace_continuation_info(o)
        V[:, i] = info.quadrature_units
        S[:, i] = info.centre_se
    end
    spread = [std(V[U, :]) for U in 1:nunits]
    reported = [sum(S[U, :]) / length(seeds) for U in 1:nunits]
    ratio = reported ./ spread
    @test all(isfinite, ratio)
    @test 0.4 < median(ratio) < 2.5
end

@testset "no unit is named a jump on a smooth objective" begin
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    g = ctsem_laplace_evaluate(laplace, theta; gradient=true).gradient
    s = -g ./ sqrt(sum(abs2, g))
    jump = ContinuousTimeSEM.CTSEMJump(true, theta, s, 1e-3, 1e-5, 0.3)
    named = ContinuousTimeSEM._ctsem_jump_units(laplace, jump)
    @test isempty(named.units)
    @test named.point_units == [0]
end

# `_fresh_2d3` with a standard normal prior on every raw parameter. Without
# one the scales run to zero, past raw 1e7, where Laplace is exact and every
# rule agrees, so a run has nothing to converge on.
function _fresh_2d3_prior()
    starts, times, data = _laplace_test_data(6, 2; seed=7)
    prior = ctsem_objective(_CONTINUATION_2D_OBJECTIVE.params, starts, times, data;
        prior_index=collect(1:12), prior_scale=ones(12))
    return (ctsem_laplace_objective(prior, [1, 3, 5], [7, 8, 9], [10, 11, 12],
        [1.0, 1.0, 1.0]), last(_fresh_2d3()))
end

@testset "a run on the placed objective starts from the pairs the last one ended with" begin
    # A run asked to (`lbfgs_pairs = :objective`, as the refined fit's climb
    # asks, `.ctLaplacePlacedClimb()`) leaves its curvature on the objective,
    # and the next run that asks starts from it. From a point off the first
    # run's optimum the warm run needs fewer iterations than a cold one.
    C = ContinuousTimeSEM
    laplace, values = _fresh_2d3_prior()
    theta = collect(Float64, values)
    n = length(theta)
    P = ctsem_laplace_placed(ctsem_laplace_continuation(laplace, theta;
        force=[1, 2], screen=false))
    @test length(P.placed) == 2
    stage(from; kw...) = C.ctsem_optimize(P, from; maxiter=500, g_tol=1e-6,
        newton=false, certify=false, progress=false, verbose=false,
        overshoot_probe=:off, kw...)
    @test P.lbfgs === nothing
    before = stage(theta; lbfgs_pairs=:objective)
    held = P.lbfgs
    @test held isa C.CTSEMLBFGSMemory
    @test !isempty(held.S) && all(v -> length(v) == n, held.S)
    x1 = collect(Float64, before.minimizer)
    @test maximum(abs, x1) < 10
    x2 = x1 .+ 0.1
    cold = stage(x2)
    # A run that does not ask leaves the held pairs alone.
    @test P.lbfgs === held
    warm = stage(x2; lbfgs_pairs="objective")
    @test P.lbfgs !== held
    @test cold.iterations > 2
    @test warm.iterations < cold.iterations
    @test isapprox(warm.maximum_loglik, cold.maximum_loglik; atol=1e-8)
    @test isapprox(warm.minimizer, cold.minimizer; atol=1e-4)
    # A memory for another number of parameters is refused, given or held.
    bad = C.CTSEMLBFGSMemory(5)
    @test C._ctsem_lbfgs_push!(bad, ones(n + 1), fill(2.0, n + 1))
    @test_throws DimensionMismatch stage(x1; lbfgs_pairs=bad)
    P.lbfgs = bad
    @test_throws DimensionMismatch stage(x1; lbfgs_pairs=:objective)
    @test_throws ArgumentError stage(x1; lbfgs_pairs=:elsewhere)
    # And an objective that keeps none starts cold.
    @test C._ctsem_lbfgs_held(laplace) === nothing
end

# --- nested units by frozen draws ---------------------------------------------
#
# Under `draws_all` a unit whose block tree nests gets draws per block, each
# block's conditional on a draw of the blocks above
# (`_continuation_importance_rule`). What has to hold, as for one block: its
# value is the unit's integral to within its own error, at two levels and at
# three; its gradient is the derivative of its value; and its standard error
# read again at the centre is the one it was placed to.

# Two levels on a nonlinear integrand: five subjects in two studies, each
# subject's DRIFT (through `-log1p_exp`) and T0MEANS varying, each study's
# MANIFESTMEANS. Raw layout: 1-4 the model parameters, 5-6 the subject scales,
# 7 their correlation, 8 the study scale. The drift's scale is wide enough
# (raw 1.5) that Laplace is 0.30 and 0.23 below the integral, which the
# product rule at 7 nodes has to 3e-4 (against 9).
_fresh_nested_nonlinear() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 2, 4], sd_index=[5, 6, 8], cor_index=[7],
    sd_scale=[1.0, 1.0, 1.0], level_nre=[2, 1],
    group=vcat(1:5, [1, 1, 1, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, 1.5, 0.2, -0.3])

@testset "each flagged unit reports its rule's error, and the hybrid its total" begin
    # A product rule's error is its value against the rule two nodes fewer,
    # placed at the same centre; a frozen rule's is its standard error; an
    # unflagged unit has none; and the hybrid's error is the draws' in
    # quadrature plus the node errors plus what the screen leaves. With the
    # draws' target out of reach no rule is refined (below).
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    five = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace, theta;
        tolerance=-1.0, draws_se=1e9))
    three = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace, theta;
        tolerance=-1.0, nodes=3, draws_se=1e9))
    @test five.nflagged == five.nunits
    for U in five.flagged
        @test isapprox(five.centre_error[U],
            abs(five.quadrature_units[U] - three.quadrature_units[U]); atol=1e-8)
    end
    @test five.error_draws == 0
    @test five.error_nodes ≈ sum(five.centre_error[U] for U in five.flagged)
    @test five.error ≈ five.error_draws + five.error_nodes + five.error_residual
    # Nothing flagged: no unit has an error, and the hybrid's is the residual.
    none = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace, theta;
        tolerance=1e6))
    @test none.nflagged == 0
    @test all(isnan, none.centre_error)
    @test none.error == none.error_residual
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    nunits = length(laplace.units.members)
    drawn = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace, theta;
        tolerance=0.0, draws_all=true, draws=256, force=1:nunits))
    @test all(drawn.centre_error[U] == drawn.centre_se[U] for U in drawn.flagged)
    @test drawn.error_draws ≈ sqrt(sum(abs2, drawn.centre_se[drawn.flagged]))
    @test drawn.error_nodes == 0
end

@testset "a rule that has not converged is refined, and takes draws only where no family converges" begin
    # On the nonlinear fixture, at a point where some unit's 5-node rule
    # differs from its 3-node rule by more than the draws' target: with the
    # target out of reach the errors are those differences; at the target
    # every flagged unit either reads within it, carries a finer rule, or is
    # on draws, and the hybrid's error is no larger. See
    # `_continuation_refine!`.
    C = ContinuousTimeSEM
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    coarse = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws_se=1e9)
    c = ctsem_laplace_continuation_info(coarse)
    worst = maximum(c.centre_error[U] for U in c.flagged if isfinite(c.centre_error[U]))
    @test worst > 0
    # A target the worst unit's cheap check exceeds: that unit is refined or
    # drawn, the rest keep their rule and their cheap measure.
    target = 0.5 * worst
    fine = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws=256,
        draws_se=target)
    f = ctsem_laplace_continuation_info(fine)
    @test f.nflagged == c.nflagged
    refined = 0
    for (k, U) in enumerate(fine.flagged)
        nodes = size(fine.rules[k].roots[1].points, 2)
        if c.centre_error[U] <= target
            @test f.centre_error[U] == c.centre_error[U]
            @test nodes == 5 && !fine.importance[U]
        else
            refined += 1
            if fine.importance[U]
                @test f.centre_error[U] == f.centre_se[U]
            else
                @test f.centre_error[U] <= target || nodes > 5
                @test nodes <= C._CONTINUATION_NODES_MAX
            end
        end
    end
    @test refined >= 1
    @test f.error <= c.error + 1e-12
    # A target only an exact rule reaches, on the three-effect nonlinear
    # fixture, whose product rules differ at every node count. Past nine
    # nodes a unit whose product rules converge -- the last node-step
    # difference below the one before -- keeps its finest product rule with
    # the geometric tail as its error, never less than that difference, or
    # the lattice's finest where that one's spread is smaller; it takes draws
    # only where they do not converge. With draws off it keeps the rule with
    # the smallest error either way. (The one-effect fixture above reads its
    # 7-node rule within 1e-12 of its 5-node rule, so nothing there
    # escalates.) See `_continuation_refine!`.
    laplace, values = _fresh_2d3()
    theta = collect(Float64, values)
    tiny = 1e-12
    at(m) = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace,
        theta; tolerance=-1.0, nodes=m, draws_se=1e9))
    i5, i7, i9 = at(5), at(7), at(C._CONTINUATION_NODES_MAX)
    function tail(U)
        d9 = abs(i9.quadrature_units[U] - i7.quadrature_units[U])
        d7 = abs(i7.quadrature_units[U] - i5.quadrature_units[U])
        converging = d7 > 0 && d9 < d7
        r = converging ? d9 / d7 : NaN
        return (d9=d9, converging=converging,
            error=converging ? d9 * max(1.0, r / (1 - r)) : d9)
    end
    # The lattice's finest rule, and its spread, for the comparison.
    full = C._laplace_extend(laplace, theta)
    Ls = C._laplace_popchols(full, laplace.spec)
    ws = C._laplace_workspace!(laplace, Float64, length(full))
    finest(U, o) = C._continuation_unit_rule(laplace, U, full, Ls, ws,
        merge(C._continuation_opts(o),
            (choice=C.CTSEMRuleChoice(:lattice, 5, C._LATTICE_MAX),)))
    drawn = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws=256,
        draws_se=tiny)
    d = ctsem_laplace_continuation_info(drawn)
    checked = 0
    for U in drawn.flagged
        t = tail(U)
        t.d9 <= tiny && continue
        checked += 1
        @test drawn.importance[U] == !t.converging
        if drawn.importance[U]
            @test d.centre_error[U] == d.centre_se[U]
        else
            f = finest(U, drawn)
            if f.se < t.error
                @test drawn.unit_rules[U] == C.CTSEMRuleChoice(:lattice, 5, C._LATTICE_MAX)
                @test isapprox(d.centre_error[U], f.se; rtol=1e-8)
            else
                @test isapprox(d.quadrature_units[U], i9.quadrature_units[U]; rtol=1e-10)
                @test isapprox(d.centre_error[U], t.error; rtol=1e-8)
                @test d.centre_error[U] >= t.d9
            end
        end
    end
    @test checked >= 1
    nodrawn = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws=0,
        draws_se=tiny)
    n = ctsem_laplace_continuation_info(nodrawn)
    @test !any(nodrawn.importance)
    for U in nodrawn.flagged
        t = tail(U)
        t.d9 <= tiny && continue
        f = finest(U, nodrawn)
        @test isapprox(n.centre_error[U], min(t.error, f.se); rtol=1e-8)
    end
end

@testset "a nested unit under draws_all gets frozen draws per block" begin
    # At the default `draws_df`: a t at the top block and Gaussians beneath.
    # A t beneath read 0.03-0.12 nats high here (`_continuation_importance_rule`).
    for (label, fresh, draws, draws_max) in (
            ("two levels", _fresh_twolevel, 64, 4096),
            # Draws at three levels multiply: one batch of 32 at each.
            ("three levels", _fresh_threelevel, 32, 32))
        laplace, values = fresh()
        theta = collect(Float64, values)
        # Linear Gaussian, so Laplace is the integral.
        lap = ctsem_laplace_evaluate(laplace, theta; gradient=false)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
            draws=draws, draws_max=draws_max, force=1:length(laplace.units.members))
        info = ctsem_laplace_continuation_info(o)
        @test (label, info.nimportance, info.nflagged) ==
            (label, info.nunits, info.nunits)
        # Each draw of a root carries one rule per block beneath it.
        root = o.rules[1].roots[1]
        @test (label, !isempty(root.children) && length(root.children) ==
            size(root.points, 2)) == (label, true)
        @test (label, info.importance_draws > sum(size(r.roots[1].points, 2)
            for r in o.rules)) == (label, true)
        for U in 1:info.nunits
            @test (label, U, abs(info.quadrature_units[U] - lap.unit_loglik[U]) <
                4 * info.centre_se[U] + 1e-3) == (label, U, true)
        end
        @test (label, isapprox(ctsem_laplace_continuation_se(o, theta),
            info.centre_se[info.flagged]; rtol=1e-8)) == (label, true)
        r = _continuation_fd_check(o, theta .+ 0.05)
        for j in eachindex(theta)
            @test (label, j, isapprox(r.g[j], r.fd[j]; atol=1e-6, rtol=1e-5)) ==
                (label, j, true)
        end
    end
end

@testset "nested draws on a nonlinear integrand integrate as the product rule does" begin
    laplace, values = _fresh_nested_nonlinear()
    theta = collect(Float64, values)
    # Seven nodes on each of a path's three effects: the reference.
    reference = ctsem_laplace_continuation_info(ctsem_laplace_continuation(laplace,
        theta; nodes=7, tolerance=-1.0, draws_se=1e9)).quadrature_units
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=false).unit_loglik
    o = ctsem_laplace_continuation(laplace, theta; tolerance=0.0, draws_all=true,
        draws=128, draws_se=0.02, force=1:length(laplace.units.members))
    info = ctsem_laplace_continuation_info(o)
    @test info.nimportance == info.nunits
    for U in 1:info.nunits
        @test abs(info.quadrature_units[U] - reference[U]) <
            4 * info.centre_se[U] + 1e-3
        # And it is a bar Laplace fails.
        @test abs(lap[U] - reference[U]) > 8 * info.centre_se[U]
    end
    r = _continuation_fd_check(o, theta .+ 0.05)
    @test isapprox(r.g, r.fd; atol=1e-6, rtol=1e-5)
end
