# The marginal-filter categorical generator, `_generate_binary!`
# (kalman_filters.jl:352), against the round-trip identity its route is
# documented to satisfy.
#
# `ctsem_generate` (kalman_trace.jl:728) drives `_generate_binary!` for a
# categorical row and is checked against the defining identity -- a dataset
# the filter generated must get back, from itself, the same log likelihood it
# reported while generating it -- for Gaussian and state-dependent-nonlinear
# models at test_kalman_trace.jl's "generated data is a draw from the model
# the filter conditions on". Grep across this suite finds `ctsem_generate`
# never called with a categorical `manifesttype`, so that identity, and the
# values `_generate_binary!` actually produces, are unverified for the branch
# that runs it.
#
# This is deliberately not `ctsem_generate_states` (state_sampling.jl), the
# *other* generation route, already checked against closed-form conditional
# distributions at test_state_sampling.jl's "conditional draws have the
# distribution the model gives them". Per R/ctGenerate.R, `intoverstates =
# 'auto'` resolves to FALSE for every categorical model, so an R-level
# `ctGenerate()` test (test-julia-binary.R and friends) exercises only that
# other route, never this one -- see CLAUDE.md, "Generation changed under
# you". This file is this route's own check, at the Julia level, against
# `ctsem_generate` directly.

using Random, Statistics

# One latent, one categorical manifest, with real dynamics (nonzero DIFFUSION
# and T0VAR) so the filter's categorical update has a moving state to
# condition on between rows, not just a fixed marginal -- the case the
# round-trip identity actually has work to do on.
function _generate_categorical_objective(kind::Int; nrows=8, ncategories=0,
    manifestmeans=0.3, thresholds=nothing, censormin=Float64[],
    censormax=Float64[], manifestvar=0.0, t0var=1.0)
    mats = Pair[:DRIFT => [-0.5;;], :JAx => [-0.5;;],
        :CINT => reshape([0.0], :, 1), :DIFFUSION => [0.6;;],
        :LAMBDA => [1.0;;], :Jy => [1.0;;],
        :MANIFESTMEANS => reshape([manifestmeans], :, 1),
        :MANIFESTVAR => [manifestvar;;],
        :T0VAR => [t0var;;], :T0MEANS => reshape([0.0], :, 1),
        :PARS => reshape([0.0], :, 1)]
    thresholds === nothing ||
        push!(mats, :THRESHOLDS => reshape(thresholds, 1, :))
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    values = Float64[]
    for (name, mat) in mats
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name)
            push!(rows, i)
            push!(cols, j)
            push!(values, mat[i, j])
        end
    end
    axis = ContinuousTimeSEM.retrieve_axes(matrices, rows, cols)
    n = length(values)
    sp = ContinuousTimeSEM.EKFParameters(falses(n), fill(false, n),
        fill(false, n), fill(false, n), fill(false, n), Function[], Function[],
        Function[], Function[], Int[], axis, fill(true, n),
        AbstractFloat[values...], Int[], Int[], Int[], [1],
        # Positional, so every trailing argument has to be counted: the
        # `nasymptotes` slot sits between `ncategories` and `censormin` and is
        # named here rather than left out. Leaving it out is silent, because
        # `Vector{Int}([0.0])` succeeds -- a censormin of 0.0 became an
        # asymptote count, the censormax became the censormin, and every
        # censored draw came back clamped at a lower limit of 5.0.
        true, [kind], [ncategories], Int[], censormin, censormax)
    nrows_ = nrows
    times = collect(0.0:1.0:(nrows_ - 1))
    data = zeros(1, nrows_)
    return (sp=sp,
        objective=ContinuousTimeSEM.ctsem_objective(sp, [1], times, data),
        times=times)
end

@testset "categorical generation is a draw from the filter's own predictive" begin
    nrows = 8
    cases = [
        # label, kind, valuecheck(y -> Bool), extra kwargs...
        ("binary", ContinuousTimeSEM.CTSEM_OBS_BINARY,
            y -> all(v -> v == 0.0 || v == 1.0, y), 1001,
            (manifestmeans=0.3,)),
        ("ordinal", ContinuousTimeSEM.CTSEM_OBS_ORDINAL,
            y -> all(v -> v in (1.0, 2.0, 3.0, 4.0), y), 1002,
            (manifestmeans=0.0, ncategories=4, thresholds=[-1.0, 1.0, 1.0])),
        ("count", ContinuousTimeSEM.CTSEM_OBS_COUNT,
            y -> all(v -> v >= 0 && v == round(v), y), 1003,
            (manifestmeans=0.8,)),
        ("censored", ContinuousTimeSEM.CTSEM_OBS_CENSORED,
            y -> all(v -> v >= 0.0 - 1e-8 && v <= 5.0 + 1e-8, y), 1004,
            (manifestmeans=2.5, censormin=[0.0], censormax=[5.0],
                manifestvar=0.5)),
    ]

    for (label, kind, valuecheck, seed, kwargs) in cases
        @testset "$label" begin
            setup = _generate_categorical_objective(kind; nrows=nrows, kwargs...)
            base = randn(MersenneTwister(seed), 1, nrows)
            g = ContinuousTimeSEM.ctsem_generate(setup.objective, Float64[], base)

            @test size(g.Y) == (1, nrows)
            @test all(isfinite, g.Y)
            @test valuecheck(vec(g.Y))
            # Not every draw identical, or the model's dynamics are not being
            # exercised and the round-trip below would be trivial.
            @test length(unique(vec(g.Y))) > 1

            # The defining identity: the filter's own likelihood for the
            # generated data must be the likelihood it reported while
            # generating it. Exact, per test_kalman_trace.jl's Gaussian
            # counterpart -- nothing here is stochastic once `g.Y` is fixed.
            regenerated = ContinuousTimeSEM.ctsem_objective(setup.sp, [1],
                setup.times, g.Y)
            @test ContinuousTimeSEM.ctsem_evaluate(regenerated, Float64[];
                gradient=false).value ≈ sum(g.subject_loglik) rtol = 1e-12
            @test sum(g.llrow) ≈ sum(g.subject_loglik) rtol = 1e-12

            # Different base normals must give different data: the draws
            # drive the result rather than the marginal alone.
            other_base = randn(MersenneTwister(seed + 1), 1, nrows)
            other = ContinuousTimeSEM.ctsem_generate(setup.objective, Float64[],
                other_base)
            @test !isapprox(other.Y, g.Y; rtol=1e-6)
        end
    end
end



# The count branch, which is the one categorical kind whose support is
# unbounded and so the one that cannot be drawn by inverting its marginal
# the way binary and ordinal are. It is sampled in the two stages that define
# the marginal instead -- predictor, then count given it -- so what has to be
# checked is that the composition has the marginal it should, and that the
# second stage's randomness is still governed by `set.seed()` now that it does
# not come from `base`.
@testset "count generation" begin
    # One row per dataset, so every draw is from the same t0 predictive and the
    # target is closed form: the linear predictor is MANIFESTMEANS + LAMBDA
    # T0MEANS and its variance LAMBDA T0VAR LAMBDA'. Rows after the first are
    # not iid -- the filter conditions on what it drew -- so replicates have to
    # come from repeated datasets rather than from one long one.
    #
    # `predictor_var` is the variance of that predictor, which is what the
    # closed forms below are written in. T0VAR is supplied to the fixture as its
    # *factor* rather than as a variance -- an entry of 0.5 measures as 0.249 on
    # a Gaussian indicator -- so the square root is what goes in.
    function count_draws(n; manifestmeans, predictor_var, seedbase = 0)
        setup = _generate_categorical_objective(
            ContinuousTimeSEM.CTSEM_OBS_COUNT; nrows = 1,
            manifestmeans = manifestmeans, t0var = sqrt(predictor_var))
        ys = Vector{Float64}(undef, n)
        for i in 1:n
            base = randn(MersenneTwister(seedbase + i), 1, 1)
            g = ContinuousTimeSEM.ctsem_generate(setup.objective, Float64[],
                base; seed = seedbase + i)
            ys[i] = g.Y[1, 1]
        end
        return ys
    end

    @testset "the draw has the Poisson-lognormal marginal" begin
        # Written out here rather than taken from the code under test. The
        # `E[y]` term in the variance is the Poisson's own share: a draw that
        # carried only the state's would be too narrow, which is precisely the
        # kind of error a posterior predictive check exists to detect and would
        # instead be reporting as a property of the model.
        for (mm, pvar, n) in ((0.8, 0.25, 40000), (2.0, 0.25, 40000),
                (0.8, 1.0, 40000))
            ys = count_draws(n; manifestmeans = mm, predictor_var = pvar)
            @test all(v -> v >= 0 && v == round(v), ys)
            expected_mean = exp(mm + pvar / 2)
            expected_var = expected_mean + expected_mean^2 * expm1(pvar)
            @test mean(ys) ≈ expected_mean rtol = 0.05
            @test sqrt(var(ys)) ≈ sqrt(expected_var) rtol = 0.1
        end
    end

    @testset "a rate past Int64 range is a number, not an InexactError" begin
        # The reported bug: the walk that used to serve this branch sized
        # itself with `Int(ceil(mean_rate + ...))`, and Julia evaluates
        # arguments before the `min` that was meant to cap it, so a rate of
        # 3.3e22 threw from inside generation. Nothing here may throw, and
        # nothing may come back `Inf` either -- an infinite draw takes the rest
        # of the row's likelihood with it, which is worse than the report.
        for mm in (30.0, 52.0, 200.0, 400.0)
            ys = count_draws(40; manifestmeans = mm, predictor_var = 1.0)
            @test all(isfinite, ys)
            @test all(v -> v >= 0, ys)
            # Far above the handover the Poisson's relative spread is
            # negligible against the state's, so the draw is essentially
            # `exp(η)` and its median sits at `exp(manifestmeans)`. Clamped at
            # `_CTSEM_COUNT_MAX_LOG_RATE`, which is what keeps it finite rather
            # than `Inf` once the predictor passes 200.
            @test log(median(ys)) ≈ min(mm, 200.0) atol = 1.5
        end
    end

    @testset "set.seed() still governs the second stage" begin
        # `base` no longer carries all the randomness a count needs, so the
        # property that a generated dataset is reproducible from the seed the
        # user set is no longer automatic -- it holds because R draws `seed`
        # with its own generator and the engine derives its streams from it.
        setup = _generate_categorical_objective(
            ContinuousTimeSEM.CTSEM_OBS_COUNT; nrows = 6, manifestmeans = 2.0)
        base = randn(MersenneTwister(77), 1, 6)
        gen(seed) = ContinuousTimeSEM.ctsem_generate(setup.objective,
            Float64[], base; seed = seed).Y
        @test gen(11) == gen(11)
        @test gen(11) != gen(12)
        # And `base` still moves the draw on its own, so the state's deviate
        # has not quietly stopped mattering.
        other = randn(MersenneTwister(78), 1, 6)
        @test ContinuousTimeSEM.ctsem_generate(setup.objective, Float64[],
            other; seed = 11).Y != gen(11)
    end

    @testset "streams are per subject, not per dataset" begin
        # `_ctsem_generate_rng` keys on the subject index as well as the seed,
        # so no two subjects share draws and the order they are visited in
        # cannot change a result. Distinct subject indices must give distinct
        # streams for the same seed.
        a = ContinuousTimeSEM._ctsem_generate_rng(4321, 1)
        b = ContinuousTimeSEM._ctsem_generate_rng(4321, 2)
        c = ContinuousTimeSEM._ctsem_generate_rng(4321, 1)
        @test rand(a, 8) != rand(b, 8)
        @test rand(c, 8) == rand(ContinuousTimeSEM._ctsem_generate_rng(4321, 1), 8)
    end

    @testset "log(y!) survives a count past Int64 range" begin
        # The same conversion on the likelihood side: `_log_factorial` counted
        # in `Int`, so evaluating the likelihood of a saturated draw threw
        # rather than evaluating Stirling on the float it already had. Also
        # reachable from data, not only from generation.
        @test ContinuousTimeSEM._log_factorial(0.0) == 0.0
        @test ContinuousTimeSEM._log_factorial(5.0) ≈ log(120.0)
        @test ContinuousTimeSEM._log_factorial(1e6) ≈ 1.2815518384658169e7 rtol = 1e-12
        @test isfinite(ContinuousTimeSEM._log_factorial(3.8e90))
        for y in (400.0, 1e6, 5e21, 3.8e90)
            @test ContinuousTimeSEM._category_loglikelihood(50.0, y, (),
                ContinuousTimeSEM.CTSEM_OBS_COUNT) <= 0
        end
    end
end
