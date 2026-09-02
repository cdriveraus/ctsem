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
