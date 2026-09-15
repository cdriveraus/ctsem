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
    censormax=Float64[], manifestvar=0.0)
    mats = Pair[:DRIFT => [-0.5;;], :JAx => [-0.5;;],
        :CINT => reshape([0.0], :, 1), :DIFFUSION => [0.6;;],
        :LAMBDA => [1.0;;], :Jy => [1.0;;],
        :MANIFESTMEANS => reshape([manifestmeans], :, 1),
        :MANIFESTVAR => [manifestvar;;],
        :T0VAR => [1.0;;], :T0MEANS => reshape([0.0], :, 1),
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
        true, [kind], [ncategories], censormin, censormax)
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


# A count is the one categorical kind with no upper bound, so the marginal walk
# that serves the others cannot always serve it -- the walk costs one
# quadrature per value, and past a rate of a few hundred every individual
# marginal probability underflows, so it accumulates no mass at all. These are
# the two regimes `_generate_count_marginal` splits into and the properties
# that have to hold across both of them.
@testset "count generation across both regimes" begin
    nodes, weights = ContinuousTimeSEM._gauss_hermite(
        ContinuousTimeSEM._CTSEM_BINARY_NODES[])
    draw(etabar, s, z) = ContinuousTimeSEM._generate_count_marginal(etabar, s,
        ContinuousTimeSEM._standard_normal_cdf(z), z, nodes, weights)

    @testset "a rate past Int64 range is a number, not an InexactError" begin
        # The reported bug: `Int(ceil(mean_rate + ...))` was evaluated before
        # the `min` meant to cap it, so a rate of 3.3e22 threw from inside
        # generation. Nothing here may throw, and nothing may come back `Inf`
        # either -- an infinite draw takes the rest of the row's likelihood
        # with it, which is a worse failure than the one reported.
        for etabar in (30.0, 52.0, 200.0, 400.0), s in (0.0, 1.0, 20.0, 400.0)
            for z in (-3.0, 0.0, 2.5)
                y = draw(etabar, s, z)
                @test isfinite(y)
                @test y >= 0
                @test y == round(y)
            end
        end
    end

    @testset "the walk has the marginal distribution" begin
        # Where the walk runs it inverts the marginal exactly, so sweeping `z`
        # over the standard normal and weighting by its density must reproduce
        # each value's marginal probability -- which the quadrature reports
        # independently as that value's own log marginal likelihood.
        etabar, s = 1.5, 0.8
        zs = range(-4.0, 4.0; length = 8001)
        ys = [draw(etabar, s, z) for z in zs]
        w = [exp(-z * z / 2) for z in zs]
        w ./= sum(w)
        for k in 0:8
            empirical = sum(w[i] for i in eachindex(ys) if ys[i] == k; init = 0.0)
            logZ, _, _ = ContinuousTimeSEM._binary_moments(etabar, s, float(k),
                nodes, weights, (), ContinuousTimeSEM.CTSEM_OBS_COUNT)
            @test empirical ≈ exp(logZ) atol = 5e-3
        end
    end

    @testset "the closed form is the Poisson-lognormal quantile" begin
        # Past the walk, the draw is the marginal's lognormal quantile matched
        # to `E[y] = exp(etabar + s^2/2)` and `Var[y]/E[y]^2 = expm1(s^2) +
        # 1/E[y]`. Those moments are written here rather than taken from the
        # code under test, and the `1/E[y]` term is the Poisson's own share of
        # the spread -- drop it and a draw at a small `s` comes out far too
        # narrow, which is exactly the kind of error a posterior predictive
        # check exists to detect and would instead be reporting.
        for (etabar, s) in ((12.0, 1.0), (20.0, 2.0), (9.0, 0.05), (52.0, 3.0))
            mean_y = exp(etabar + s^2 / 2)
            relvar = expm1(s^2) + 1 / mean_y
            sigma2 = log1p(relvar)
            mu = log(mean_y) - sigma2 / 2
            for z in (-2.0, -0.5, 0.0, 1.0, 2.5)
                @test draw(etabar, s, z) ≈ round(exp(mu + sqrt(sigma2) * z)) rtol = 1e-10
            end
        end
    end

    @testset "the draw increases with the deviate, to within one count" begin
        # Both regimes map `z` to `y` the same way -- increasing, through the
        # marginal -- so this holds across the handover as well as within
        # either side of it. It is what inverting a CDF means, and it is the
        # property that says the two regimes agree where they meet: an earlier
        # version that sampled the mixture in two stages rather than inverting
        # it dropped by half the value at the seam, which no summary of the
        # generated data would have shown as anything but a wrong model.
        #
        # One count of slack, not zero. The walk inverts the marginal exactly
        # and the closed form inverts a lognormal matched to its first two
        # moments, so at the single `z` where one hands over to the other they
        # can disagree by a count -- measured at 82 against 81, and only in the
        # one of these six cases whose marginal is heavy-tailed enough to reach
        # the handover at all.
        for (etabar, s) in ((1.5, 0.8), (2.2, 1.0), (6.0, 0.5), (6.2, 0.3),
                (12.0, 1.0), (30.0, 2.0))
            ys = [draw(etabar, s, z) for z in range(-4.0, 4.0; length = 2001)]
            @test all(isfinite, ys)
            @test maximum(ys[i] - ys[i + 1] for i in 1:(length(ys) - 1)) <= 1.0
        end
    end

    @testset "log(y!) survives a count past Int64 range" begin
        # The second conversion on the same road: `_log_factorial` counted in
        # `Int`, so the likelihood of a saturated draw threw rather than
        # evaluating Stirling on the float it already had.
        @test ContinuousTimeSEM._log_factorial(0.0) == 0.0
        @test ContinuousTimeSEM._log_factorial(5.0) ≈ log(120.0)
        @test ContinuousTimeSEM._log_factorial(1e6) ≈ 1.2815518384658169e7 rtol = 1e-12
        @test isfinite(ContinuousTimeSEM._log_factorial(3.8e90))
        # And the count likelihood it feeds stays a log probability out there.
        for y in (400.0, 1e6, 5e21, 3.8e90)
            @test ContinuousTimeSEM._category_loglikelihood(50.0, y, (),
                ContinuousTimeSEM.CTSEM_OBS_COUNT) <= 0
        end
    end
end
