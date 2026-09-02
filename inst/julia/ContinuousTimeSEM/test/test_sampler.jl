# The Hamiltonian sampler: the joint target, the trajectories, and the metric.
#
# Two kinds of test here, and the second kind is the one that earns its runtime.
#
# The gradient is checked against ForwardDiff over the whole `(theta, u)` vector,
# which is straightforward and catches a wrong chain rule immediately.
#
# The *sampler* is checked against targets whose posterior is known in closed
# form, because a wrong sampler does not error -- it returns a plausible-looking
# posterior that is quietly the wrong one. Two real bugs were found this way and
# neither would have been found any other way:
#
#   * a single candidate buffer shared across recursion depths, so a child
#     overwrote the candidate its parent was holding. On a standard normal that
#     gave posterior sd 1.03-1.10 and E[x^4] 3.7-4.9, against 1 and 3.
#   * the subtree no-U-turn check reaching one state back *outside* the subtree,
#     which breaks the reversibility argument the algorithm rests on. That gave
#     sd 0.974 on an isotropic 4-D normal.
#
# Neither shrank with ten times the draws, which is the signature to test for:
# these assertions are on moments that converge, at tolerances tight enough that
# a few percent of bias fails them.

using ForwardDiff, LinearAlgebra, Random, Statistics

# `_fresh_linear`, `_fresh_nonlinear`, `_fresh_twolevel` and `_fresh_threelevel`
# come from `laplace_fixtures.jl`, which `test_laplace.jl` also includes. Shared
# rather than rebuilt because `ekf_from_columns` parses transform strings through
# `eval`, so a rebuilt model is uniquely typed and recompiles the whole filter.
isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

################################################################################
# The joint density and its gradient
################################################################################

"""The joint log density written the obvious way, for ForwardDiff to chew on."""
function _sampler_reference(sampler, x::AbstractVector{T}) where {T}
    lp = sampler.laplace
    spec = lp.spec
    theta = x[1:sampler.npar]
    Ls = ContinuousTimeSEM._laplace_popchols(theta, spec)
    total = zero(T)
    for U in 1:sampler.nunits
        range = (sampler.uoffsets[U] + 1):(sampler.uoffsets[U] + sampler.udims[U])
        u = x[range]
        for (m, i) in enumerate(lp.units.members[U])
            shifted = ContinuousTimeSEM._laplace_member_values(theta, spec, Ls, u,
                lp.units.offsets[U][m])
            total += lp.objective.subject_objectives[i](shifted)
        end
        total -= dot(u, u) / 2
    end
    return total + ContinuousTimeSEM._ctsem_log_prior(lp.objective, theta)
end

@testset "the joint gradient matches ForwardDiff at every level" begin
    for (label, fresh) in (("one level", _fresh_linear),
        ("nonlinear", _fresh_nonlinear), ("two levels", _fresh_twolevel),
        ("three levels", _fresh_threelevel))
        laplace, values = fresh()
        sampler = ContinuousTimeSEM.ctsem_sampler(laplace, length(values))
        x = ContinuousTimeSEM.ctsem_sample_start(sampler, values)
        # Non-zero effects, or the population scale's contribution through
        # `dL * u` is multiplied by zero and the part most likely to be wrong is
        # never exercised at all.
        for a in (sampler.npar + 1):sampler.ndim
            x[a] = 0.35 * sin(1.7a)
        end
        got = ContinuousTimeSEM.ctsem_sample_density(sampler, x)
        @test got.value ≈ _sampler_reference(sampler, x) rtol = 1e-12
        @test got.gradient ≈ ForwardDiff.gradient(z -> _sampler_reference(sampler, z), x) rtol = 1e-8
    end
end

@testset "the sampled dimension is the population vector plus every effect" begin
    laplace, values = _fresh_twolevel()
    sampler = ContinuousTimeSEM.ctsem_sampler(laplace, length(values))
    @test sampler.npar == length(values)
    @test ContinuousTimeSEM.ctsem_sample_dimension(sampler) ==
        length(values) + sum(laplace.units.dims)
    # The effects start at zero, which is their prior mean and, after a fit,
    # close to their posterior mode.
    x = ContinuousTimeSEM.ctsem_sample_start(sampler, values)
    @test x[1:sampler.npar] == collect(Float64, values)
    @test all(iszero, x[(sampler.npar + 1):end])
end

@testset "a point that cannot be evaluated is a rejection, not an error" begin
    laplace, values = _fresh_linear()
    sampler = ContinuousTimeSEM.ctsem_sampler(laplace, length(values))
    x = ContinuousTimeSEM.ctsem_sample_start(sampler, values)

    # A leapfrog trajectory routinely steps somewhere the model cannot be
    # evaluated. Throwing there would end the chain rather than the trajectory,
    # so every such point must come back as -Inf with a finite gradient the
    # caller will discard.
    poisoned = copy(x)
    poisoned[1] = NaN
    got = ContinuousTimeSEM.ctsem_sample_density(sampler, poisoned)
    @test got.value == -Inf
    @test all(isfinite, got.gradient)

    # A population scale so large the covariance will not factorize takes the
    # same route -- outside the support, not a bug.
    scaled = copy(x)
    for level in laplace.spec.levels, j in level.sd_index
        scaled[j] = 1e6
    end
    @test isfinite(ContinuousTimeSEM.ctsem_sample_density(sampler, scaled).value) ||
        ContinuousTimeSEM.ctsem_sample_density(sampler, scaled).value == -Inf

    # Far from the data the filter degrades rather than failing, and that is
    # worth pinning: a merely terrible point must stay usable, because the
    # sampler needs its gradient to climb back out.
    distant = copy(x)
    distant[1] = 1e6
    far = ContinuousTimeSEM.ctsem_sample_density(sampler, distant)
    @test isfinite(far.value)
    @test far.value < -1e6
    @test all(isfinite, far.gradient)
end

################################################################################
# The metric
################################################################################

@testset "the metric blocks follow the model's own sparsity" begin
    laplace, values = _fresh_twolevel()
    sampler = ContinuousTimeSEM.ctsem_sampler(laplace, length(values))
    metric = ctsem_sample_metric(sampler, values)
    # Every coordinate covered exactly once: a coordinate the metric omits gets
    # zero momentum and never moves, which is a silent failure.
    covered = zeros(Int, sampler.ndim)
    for r in metric.ranges, i in r
        covered[i] += 1
    end
    @test all(==(1), covered)
    @test first(metric.ranges) == 1:sampler.npar
    @test all(issuccess(cholesky(Symmetric(C * transpose(C)); check=false))
              for C in metric.factors)
end

@testset "momentum drawn from the metric has the metric's inverse covariance" begin
    # Sigma is the *inverse* mass matrix, so momentum has covariance Sigma^-1.
    # Getting this backwards is an easy slip that leaves the sampler correct but
    # catastrophically slow, so it is asserted rather than assumed.
    S = [4.0 1.2; 1.2 1.0]
    metric = ContinuousTimeSEM._metric_from_covariances([1:2], [copy(S)])
    rng = Random.Xoshiro(3)
    p = zeros(2)
    draws = Matrix{Float64}(undef, 2, 200_000)
    for t in axes(draws, 2)
        ContinuousTimeSEM._metric_momentum!(p, metric, rng)
        draws[:, t] = p
    end
    @test cov(draws; dims=2) ≈ inv(S) rtol = 0.03

    # And the velocity is Sigma * p, with the kinetic energy sharing its work.
    v = zeros(2); scratch = zeros(2); q = [0.7, -0.3]
    kinetic = ContinuousTimeSEM._metric_velocity!(v, metric, q, scratch)
    @test v ≈ S * q
    @test kinetic ≈ dot(q, S * q) / 2
end

################################################################################
# The sampler, against posteriors that are known exactly
################################################################################

"""Run the sampler over a plain log density, outside the ctsem model."""
function _sample_target(logdensity!, ndim, metric; ndraws=20_000, nwarmup=500,
    seed=1, maxdepth=10)
    rng = Random.Xoshiro(seed)
    ws = ContinuousTimeSEM._NUTSWorkspace(ndim, maxdepth)
    x = zeros(ndim); g = zeros(ndim)
    logp = logdensity!(g, x)
    eps = ContinuousTimeSEM._init_stepsize(logdensity!, metric, rng, x, g, logp, ws)
    da = ContinuousTimeSEM._DualAverage(eps, 0.8)
    for _ in 1:nwarmup
        step = ContinuousTimeSEM._nuts_transition!(ws, logdensity!, metric, rng,
            x, g, logp, eps, maxdepth, 1000.0)
        logp = step.logp
        eps = ContinuousTimeSEM._dual_update!(da, step.accept)
    end
    eps = ContinuousTimeSEM._dual_final(da)
    draws = Matrix{Float64}(undef, ndim, ndraws)
    divergent = 0
    for t in 1:ndraws
        step = ContinuousTimeSEM._nuts_transition!(ws, logdensity!, metric, rng,
            x, g, logp, eps, maxdepth, 1000.0)
        logp = step.logp
        divergent += step.divergent
        draws[:, t] = x
    end
    return draws, divergent
end

@testset "a standard normal is sampled with the right moments" begin
    ld! = function (g, x)
        g[1] = -x[1]
        return -x[1]^2 / 2
    end
    draws, divergent = _sample_target(ld!, 1, ContinuousTimeSEM.ctsem_identity_metric(1);
        ndraws=60_000, seed=4)
    x = vec(draws)
    @test divergent == 0
    @test mean(x) ≈ 0 atol = 0.02
    # The buffer-aliasing bug gave 1.03-1.10 here and did not shrink with more
    # draws, so 1.5% is a real assertion and not a formality.
    @test std(x) ≈ 1 rtol = 0.015
    @test mean(x .^ 4) ≈ 3 rtol = 0.05
    @test quantile(x, 0.975) ≈ 1.96 atol = 0.05
end

@testset "an isotropic normal is sampled with the right scale in four dimensions" begin
    # Four dimensions specifically: the subtree U-turn bug was worst here,
    # giving sd 0.974 where one dimension gave 0.999 and sixteen gave 1.000.
    ld! = function (g, x)
        g .= .-x
        return -dot(x, x) / 2
    end
    draws, divergent = _sample_target(ld!, 4, ContinuousTimeSEM.ctsem_identity_metric(4);
        ndraws=60_000, seed=6)
    @test divergent == 0
    for j in 1:4
        @test std(view(draws, j, :)) ≈ 1 rtol = 0.02
        @test mean(view(draws, j, :)) ≈ 0 atol = 0.03
    end
end

@testset "a correlated normal is sampled with the right covariance" begin
    S = [1.0 0.8 0.3; 0.8 1.5 -0.2; 0.3 -0.2 0.6]
    P = inv(S)
    ld! = function (g, x)
        mul!(g, P, x)
        g .= .-g
        return -dot(x, P * x) / 2
    end
    # One dense block: the block metric must reproduce the target exactly when
    # it is handed the answer.
    metric = ContinuousTimeSEM._metric_from_covariances([1:3], [copy(S)])
    draws, divergent = _sample_target(ld!, 3, metric; ndraws=40_000, seed=8)
    @test divergent == 0
    @test cov(draws; dims=2) ≈ S rtol = 0.06
end

@testset "closed-form targets stay this close at four times the draws" begin
    # The three tests above run each target at one draw count, so a small,
    # persistent bias close to their tolerance would pass unnoticed -- this
    # file's own header comment names the diagnostic ("neither shrank with ten
    # times the draws") without automating it. A consistent estimator's error
    # falls as 1/sqrt(N); a bias does not. Rather than compare two random
    # draws against each other, which is itself noisy enough to flip sign
    # target to target (measured, on this exact setup: the correlated case's
    # covariance error went *up* from 40k to 160k draws by chance, on a error
    # an order of magnitude under either tolerance), each target here is run
    # once more at four times its draw count and checked against half the
    # tolerance the smaller-draw-count test above uses -- half being what
    # 1/sqrt(4) predicts for noise, and tight enough that a bias sitting near
    # the original tolerance has no room left to hide in.
    ld1! = function (g, x)
        g[1] = -x[1]
        return -x[1]^2 / 2
    end
    draws1, divergent1 = _sample_target(ld1!, 1,
        ContinuousTimeSEM.ctsem_identity_metric(1); ndraws=240_000, seed=4)
    x = vec(draws1)
    @test divergent1 == 0
    @test mean(x) ≈ 0 atol = 0.01
    @test std(x) ≈ 1 rtol = 0.0075
    @test mean(x .^ 4) ≈ 3 rtol = 0.025
    @test quantile(x, 0.975) ≈ 1.96 atol = 0.025

    ld4! = function (g, x)
        g .= .-x
        return -dot(x, x) / 2
    end
    draws4, divergent4 = _sample_target(ld4!, 4,
        ContinuousTimeSEM.ctsem_identity_metric(4); ndraws=240_000, seed=6)
    @test divergent4 == 0
    for j in 1:4
        @test std(view(draws4, j, :)) ≈ 1 rtol = 0.01
        @test mean(view(draws4, j, :)) ≈ 0 atol = 0.015
    end

    S = [1.0 0.8 0.3; 0.8 1.5 -0.2; 0.3 -0.2 0.6]
    P = inv(S)
    ld3! = function (g, x)
        mul!(g, P, x)
        g .= .-g
        return -dot(x, P * x) / 2
    end
    metric3 = ContinuousTimeSEM._metric_from_covariances([1:3], [copy(S)])
    draws3, divergent3 = _sample_target(ld3!, 3, metric3; ndraws=160_000, seed=8)
    @test divergent3 == 0
    @test cov(draws3; dims=2) ≈ S rtol = 0.03
end

################################################################################
# Diagnostics
################################################################################

@testset "R-hat and ESS say what they should on draws with known behaviour" begin
    rng = Random.Xoshiro(12)
    # Independent draws from one distribution: R-hat at 1, ESS near the count.
    agreeing = randn(rng, 1, 4000)
    diagnostics = ctsem_sample_diagnostics(agreeing, 4)
    @test diagnostics.rhat[1] < 1.02
    @test diagnostics.ess[1] > 1000

    # Chains centred somewhere different each: R-hat must notice.
    disagreeing = randn(rng, 1, 4000)
    for c in 1:4, t in 1:1000
        disagreeing[1, (c - 1) * 1000 + t] += 3.0 * c
    end
    @test ctsem_sample_diagnostics(disagreeing, 4).rhat[1] > 1.5

    # A slowly mixing chain has fewer effective draws than draws.
    correlated = zeros(1, 4000)
    for c in 1:4
        value = 0.0
        for t in 1:1000
            value = 0.95 * value + sqrt(1 - 0.95^2) * randn(rng)
            correlated[1, (c - 1) * 1000 + t] = value
        end
    end
    @test ctsem_sample_diagnostics(correlated, 4).ess[1] < 1000
end

@testset "the adaptation schedule brackets its windows the way Stan does" begin
    # Too short for a window: the Laplace metric is used as it stands, which is
    # a reasonable mode here rather than a degenerate one.
    @test isempty(ContinuousTimeSEM._adapt_windows(100))
    windows = ContinuousTimeSEM._adapt_windows(1000)
    @test !isempty(windows)
    @test first(windows) > 75            # after the opening buffer
    @test last(windows) == 1000 - 50     # and before the closing one
    @test issorted(windows)
end

################################################################################
# End to end
################################################################################

@testset "a model samples end to end and reports its diagnostics" begin
    laplace, values = _fresh_linear()
    out = ctsem_sample(laplace, values; nchains=2, nwarmup=120, ndraws=120,
        seed=5, maxdepth=7)
    @test size(out.draws) == (length(values), 2 * 120)
    @test out.npar == length(values)
    @test out.ndim > out.npar
    @test !out.saved_effects
    @test length(out.effect_mean) == out.ndim - out.npar
    @test length(out.effect_sd) == out.ndim - out.npar
    @test all(isfinite, out.draws)
    @test length(out.rhat) == out.npar
    @test length(out.stepsize) == 2
    @test all(>(0), out.stepsize)
    @test out.ndivergent >= 0

    # The effects are kept when asked for, and then the draws carry them.
    with = ctsem_sample(laplace, values; nchains=1, nwarmup=60, ndraws=60,
        seed=5, maxdepth=7, save_effects=true)
    @test size(with.draws, 1) == with.ndim
    @test with.saved_effects
end
