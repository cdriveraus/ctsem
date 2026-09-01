# The state-explicit path: sampled states, a conditional observation model, and
# the joint density of the two.
#
# The load-bearing test here is the marginalisation identity. A joint density
# over states and data is easy to get plausibly wrong -- a missing constant, a
# factor used where its transpose belongs, an innovation applied before the
# transition instead of after -- and every one of those still returns a smooth,
# finite number that moves in the right direction when you change a parameter.
# What none of them survives is integrating the states back out and landing on
# the marginal likelihood the Kalman filter reports for the same data, which is
# exact for a linear Gaussian model and so is a reference rather than a second
# opinion.
#
# The rest check the pieces that identity cannot see: that the conditional
# draws have the distributions they claim (a wrong sampler returns a plausible
# wrong dataset), and that the deterministic skeleton is the model's own
# trajectory.

using LinearAlgebra, Random, Statistics, ForwardDiff

# Matrices in, `EKFParameters` out, everything fixed. Written as a vector of
# pairs rather than a Dict so the flattened layout is deterministic.
function _state_params(mats::Vector{<:Pair}; manifesttype=Int[], ncategories=Int[],
    diffusion_state_indices=[1])
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
    return ContinuousTimeSEM.EKFParameters(falses(n), fill(false, n),
        fill(false, n), fill(false, n), fill(false, n), Function[], Function[],
        Function[], Function[], Int[], axis, fill(true, n),
        AbstractFloat[values...], Int[], Int[], Int[], diffusion_state_indices,
        true, manifesttype, ncategories)
end

# One latent, one manifest, and every matrix the filter reads. `sd` values, not
# variances: MANIFESTVAR, DIFFUSION and T0VAR go through `sdcovsqrt2cov`.
function _scalar_gaussian_params(; drift=-0.5, cint=0.2, diffusion=0.2,
    t0means=0.3, t0var=0.5, lambda=1.0, manifestmeans=1.2, manifestvar=0.4)
    return _state_params([
        :DRIFT => [drift;;], :JAx => [drift;;],
        :CINT => reshape([cint], :, 1),
        :DIFFUSION => [diffusion;;],
        :LAMBDA => [lambda;;], :Jy => [lambda;;],
        :MANIFESTMEANS => reshape([manifestmeans], :, 1),
        :MANIFESTVAR => [manifestvar;;],
        :T0VAR => [t0var;;], :T0MEANS => reshape([t0means], :, 1),
        :PARS => reshape([0.0], :, 1)])
end

# The covariance a scalar standard deviation becomes.
#
# Not `sd^2`: `sdcovsqrt2cov` regularises the correlation it builds with an
# epsilon of 1e-5, so a variance of 0.25 arrives as 0.2500025. That is a
# property of the parameterisation, tested on its own in
# `test_constrain_cor_sqrt.jl`, and taking it as given here is what keeps this
# file about the transition and the density rather than about that transform.
# Squaring by hand instead put a relative error of 1e-5 into every expected
# value below, which is far too small to look like a bug and far too large for
# an equality.
function _sdcov(sd)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 1)
    ContinuousTimeSEM.sdcovsqrt2cov!(buffer, [sd;;], 0, Val(1))
    return buffer.out[1, 1]
end

# The discrete-time quantities for the scalar model, from the closed forms
# rather than from the engine: a test that discretised through the same code it
# is checking would agree with itself.
function _scalar_discrete(dt; drift=-0.5, cint=0.2, diffusion=0.2)
    A = exp(drift * dt)
    b = (A - 1) / drift * cint
    asymptotic = -_sdcov(diffusion) / (2 * drift)
    Q = asymptotic * (1 - A^2)
    return (A=A, b=b, Q=Q)
end

@testset "innovation count follows the design, not the parameters" begin
    sp = _scalar_gaussian_params()
    times = [0.0, 1.0, 2.5]
    data = reshape([0.4, 0.9, 1.1], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
    # One latent innovation for the first row, one per row after it.
    @test ContinuousTimeSEM.ctsem_state_dimension(objective) == 3

    # Two subjects on the same schedule doubles it.
    two = ContinuousTimeSEM.ctsem_objective(sp, [1, 4],
        vcat(times, times), hcat(data, data))
    @test ContinuousTimeSEM.ctsem_state_dimension(two) == 6

    # A bounded maxtimestep splits each interval, and each substep gets its own
    # innovation: the first interval is one unit and the second one and a half,
    # so at a bound of 0.5 that is two substeps and three.
    bounded = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data,
        zeros(0, 3), zeros(1, 0), 0.5)
    @test ContinuousTimeSEM.ctsem_state_dimension(bounded) == 1 + 2 + 3
end

@testset "zero innovations give the deterministic trajectory" begin
    sp = _scalar_gaussian_params()
    times = [0.0, 1.0, 2.5]
    data = reshape([0.4, 0.9, 1.1], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
    ndim = ContinuousTimeSEM.ctsem_state_dimension(objective)
    generated = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[],
        zeros(ndim), zeros(1, 3))

    expected = 0.3
    @test generated.states[1, 1] ≈ expected
    for (index, dt) in enumerate((1.0, 1.5))
        step = _scalar_discrete(dt)
        expected = step.A * expected + step.b
        @test generated.states[1, index + 1] ≈ expected
    end
    # And with no measurement noise left to draw, the observation is the state
    # through LAMBDA and MANIFESTMEANS.
    @test generated.Y[1, 3] ≈ expected + 1.2
end

@testset "the state path builds the transition the model specifies" begin
    sp = _scalar_gaussian_params()
    times = [0.0, 1.0]
    data = reshape([0.4, 0.9], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)

    # Linear in the innovations, so one basis vector at a time recovers each
    # column of the map exactly.
    base = zeros(1, 2)
    baseline = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[],
        zeros(2), base).states
    first = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[],
        [1.0, 0.0], base).states
    second = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[],
        [0.0, 1.0], base).states

    step = _scalar_discrete(1.0)
    t0sd = sqrt(_sdcov(0.5))
    # T0VAR is given as a standard deviation of 0.5.
    @test first[1, 1] - baseline[1, 1] ≈ t0sd
    # ... and propagates through the transition to the next row.
    @test first[1, 2] - baseline[1, 2] ≈ step.A * t0sd
    # The second innovation is the process noise over that interval, and none
    # of it reaches the row before it.
    @test second[1, 1] - baseline[1, 1] ≈ 0.0
    @test second[1, 2] - baseline[1, 2] ≈ sqrt(step.Q)
end

@testset "the joint density is the conditional likelihood plus the innovations" begin
    sp = _scalar_gaussian_params()
    times = [0.0, 1.0, 2.5]
    y = [0.4, 0.9, 1.1]
    data = reshape(y, 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
    z = [0.7, -1.3, 0.25]

    # The trajectory, written out here rather than taken from the engine.
    states = zeros(3)
    states[1] = 0.3 + sqrt(_sdcov(0.5)) * z[1]
    for (index, dt) in enumerate((1.0, 1.5))
        step = _scalar_discrete(dt)
        states[index + 1] = step.A * states[index] + step.b + sqrt(step.Q) * z[index + 1]
    end

    sd = sqrt(_sdcov(0.4))
    expected = sum(-((y[t] - states[t] - 1.2) / sd)^2 / 2 - log(sd) -
        log(2pi) / 2 for t in 1:3)
    expected += sum(-value^2 / 2 - log(2pi) / 2 for value in z)

    @test ContinuousTimeSEM.ctsem_joint_loglikelihood(objective, Float64[], z) ≈ expected

    parts = ContinuousTimeSEM.ctsem_joint_evaluate(objective, Float64[], z;
        gradient=false)
    @test parts.value ≈ expected
    @test parts.observation + parts.state_prior + parts.parameter_prior ≈ parts.value
    @test parts.ndim == 3
end

@testset "integrating the states back out gives the filter's marginal" begin
    # The identity the whole path rests on. Exact for a linear Gaussian model,
    # where the EKF's marginal likelihood is the true one, so any disagreement
    # is this file's fault rather than an approximation showing.
    sp = _scalar_gaussian_params()
    times = [0.0, 1.0]
    data = reshape([0.4, 0.9], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
    marginal = objective(Float64[])

    nodes, weights = ContinuousTimeSEM._gauss_hermite(60)
    # Gauss-Hermite is written for weight exp(-x^2), so a standard normal
    # expectation is the rule at `sqrt(2) x` divided by `sqrt(pi)` per
    # dimension.
    total = 0.0
    for i in eachindex(nodes), j in eachindex(nodes)
        z = [sqrt(2) * nodes[i], sqrt(2) * nodes[j]]
        joint = ContinuousTimeSEM.ctsem_joint_loglikelihood(objective, Float64[], z)
        conditional = joint - sum(-value^2 / 2 - log(2pi) / 2 for value in z)
        total += weights[i] * weights[j] * exp(conditional)
    end
    @test log(total / pi) ≈ marginal atol = 1e-8
end

@testset "the joint density differentiates in both of its arguments" begin
    # One free parameter, so the gradient has a parameter block and a state
    # block and the packing of the two can be checked against a difference.
    matrices = [:DRIFT => [-0.5;;], :JAx => [-0.5;;],
        :CINT => reshape([0.0], :, 1), :DIFFUSION => [0.2;;],
        :LAMBDA => [1.0;;], :Jy => [1.0;;],
        :MANIFESTMEANS => reshape([0.0], :, 1), :MANIFESTVAR => [0.4;;],
        :T0VAR => [0.5;;], :T0MEANS => reshape([0.0], :, 1),
        :PARS => reshape([0.0], :, 1)]
    sp = _state_params(matrices)
    times = [0.0, 1.0]
    data = reshape([0.4, 0.9], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)

    z = [0.3, -0.8]
    result = ContinuousTimeSEM.ctsem_joint_evaluate(objective, Float64[], z)
    @test length(result.gradient) == 2
    @test all(isfinite, result.gradient)
    step = 1e-6
    for i in 1:2
        up = copy(z); up[i] += step
        down = copy(z); down[i] -= step
        numeric = (ContinuousTimeSEM.ctsem_joint_loglikelihood(objective, Float64[], up) -
            ContinuousTimeSEM.ctsem_joint_loglikelihood(objective, Float64[], down)) / (2 * step)
        @test result.gradient[i] ≈ numeric atol = 1e-5
    end
end

# Conditional draws.
#
# Two sample sizes, because a sampler with the wrong scale can still produce a
# mean that passes at one of them; the error of a correct one falls as the
# square root of the count and the error of a wrong one does not.
function _categorical_objective(kind::Int; ncategories=0, manifestmeans=0.0,
    nrows=4000)
    matrices = Pair[:DRIFT => [-0.5;;], :JAx => [-0.5;;],
        :CINT => reshape([0.0], :, 1), :DIFFUSION => [0.0;;],
        :LAMBDA => [1.0;;], :Jy => [1.0;;],
        :MANIFESTMEANS => reshape([manifestmeans], :, 1),
        :MANIFESTVAR => [0.0;;],
        :T0VAR => [0.0;;], :T0MEANS => reshape([0.0], :, 1),
        :PARS => reshape([0.0], :, 1)]
    if kind == ContinuousTimeSEM.CTSEM_OBS_ORDINAL
        # Cumulated by `_ordinal_thresholds!` from gaps, so these are the first
        # threshold and the two increments after it.
        push!(matrices, :THRESHOLDS => [-1.0 1.0 1.0])
    end
    sp = _state_params(matrices; manifesttype=[kind], ncategories=[ncategories])
    # No diffusion and no T0 variance: every state is exactly zero, so each row
    # is an independent draw at the same linear predictor and the sample is the
    # conditional distribution itself.
    times = collect(0.0:1.0:(nrows - 1))
    data = zeros(1, nrows)
    return ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
end

@testset "conditional draws have the distribution the model gives them" begin
    function draws(objective, nrows, seed)
        ndim = ContinuousTimeSEM.ctsem_state_dimension(objective)
        # Seeded here, because the engine has no RNG of its own: every standard
        # normal it uses arrives from the caller, which is what makes a
        # generated dataset reproducible from R's `set.seed()` and a failure
        # here reproducible from this line.
        base = randn(MersenneTwister(seed), 1, nrows)
        generated = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[],
            zeros(ndim), base)
        return vec(generated.Y)
    end

    # Bernoulli at a linear predictor of 0.7.
    for nrows in (1000, 4000)
        objective = _categorical_objective(ContinuousTimeSEM.CTSEM_OBS_BINARY;
            manifestmeans=0.7, nrows=nrows)
        y = draws(objective, nrows, 20260901)
        @test all(v -> v == 0.0 || v == 1.0, y)
        @test mean(y) ≈ 1 / (1 + exp(-0.7)) atol = 3 / sqrt(nrows)
    end

    # Poisson at a rate of exp(1.1).
    for nrows in (1000, 4000)
        objective = _categorical_objective(ContinuousTimeSEM.CTSEM_OBS_COUNT;
            manifestmeans=1.1, nrows=nrows)
        y = draws(objective, nrows, 20260902)
        rate = exp(1.1)
        @test all(v -> v >= 0 && v == round(v), y)
        @test mean(y) ≈ rate atol = 4 * sqrt(rate / nrows)
        # A Poisson has variance equal to its mean, which is what separates it
        # from every other count-shaped answer a wrong draw could give.
        @test var(y) ≈ rate atol = 0.35 * rate
    end

    # Four ordinal categories at thresholds -1, 0, 1.
    for nrows in (1000, 4000)
        objective = _categorical_objective(ContinuousTimeSEM.CTSEM_OBS_ORDINAL;
            ncategories=4, nrows=nrows)
        y = draws(objective, nrows, 20260903)
        @test all(v -> v in (1.0, 2.0, 3.0, 4.0), y)
        cumulative = [1 / (1 + exp(-(-1.0))), 0.5, 1 / (1 + exp(-1.0))]
        expected = [cumulative[1], cumulative[2] - cumulative[1],
            cumulative[3] - cumulative[2], 1 - cumulative[3]]
        for k in 1:4
            @test mean(y .== k) ≈ expected[k] atol = 3 / sqrt(nrows)
        end
    end
end

@testset "a count far into its own tail does not move the state" begin
    # The failure this path exists to make impossible. The filter's categorical
    # update is an assumed-density projection, so an improbable count moves the
    # state it conditions on and the next row is drawn from a rate that has
    # already moved; here no observation touches a state at all, so a large
    # draw is a large draw and nothing else.
    sp = _state_params([:DRIFT => [-0.4;;], :JAx => [-0.4;;],
            :CINT => reshape([0.0], :, 1), :DIFFUSION => [1.0;;],
            :LAMBDA => [1.0;;], :Jy => [1.0;;],
            :MANIFESTMEANS => reshape([0.8], :, 1), :MANIFESTVAR => [0.0;;],
            :T0VAR => [1.0;;], :T0MEANS => reshape([0.0], :, 1),
            :PARS => reshape([0.0], :, 1)];
        manifesttype=[ContinuousTimeSEM.CTSEM_OBS_COUNT])
    nrows = 12
    times = collect(0.0:1.0:(nrows - 1))
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], times, zeros(1, nrows))
    ndim = ContinuousTimeSEM.ctsem_state_dimension(objective)

    # A five-sigma first innovation, which is what put a count of 252 into the
    # report, and zero afterwards: the state must decay back rather than run.
    z = zeros(ndim)
    z[1] = 5.0
    generated = ContinuousTimeSEM.ctsem_generate_states(objective, Float64[], z,
        fill(2.5, 1, nrows))
    @test all(isfinite, generated.states)
    @test all(isfinite, generated.Y)
    # exp(0.8 + 5) is about three hundred, and every later row has a state
    # decaying towards zero, so nothing may exceed the first row's rate.
    @test generated.Y[1, 1] > 100
    @test all(generated.Y[1, 2:end] .< generated.Y[1, 1])
    # The state relaxes at the rate DRIFT gives it and at no other: over the
    # eleven remaining intervals that is `exp(-0.4 * 11)`, about 1.2% of where
    # it started. The filter path, conditioning on the count it drew, walked
    # the other way.
    @test generated.states[1, end] ≈ generated.states[1, 1] * exp(-0.4 * 11)
    @test abs(generated.states[1, end]) < abs(generated.states[1, 1]) / 50
end


# Fitting over the joint density.
#
# Two implementations of the same arithmetic are compared throughout: the
# per-subject one the engine uses, which exploits the fact that a subject's
# innovations enter no other subject's likelihood, and the naive dense one over
# the whole vector. They must agree exactly, because the blocks the fast one
# skips are structurally zero rather than small -- so a disagreement is a
# scatter written to the wrong offset, which is the failure mode of every
# blocked implementation and is invisible in the value alone.

# Free parameters, so there is a population block to differentiate. DRIFT and
# the manifest mean are estimated; everything else is fixed.
function _joint_params(; manifesttype=Int[], ncategories=Int[])
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    parnumber = Int[]
    values = Float64[]
    transforms = String[]
    function add!(name::Symbol, mat::AbstractMatrix, pars=nothing, tf="")
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name); push!(rows, i); push!(cols, j)
            push!(parnumber, pars === nothing ? 0 : pars[i, j])
            push!(values, pars === nothing || pars[i, j] == 0 ? mat[i, j] : NaN)
            push!(transforms, pars === nothing || pars[i, j] == 0 ? "" : tf)
        end
    end
    # A transform on the drift, as ctsem's own models have. Without one a raw
    # value of zero *is* a drift of zero, which makes JAx singular and the
    # discretisation NaN -- so the optimiser would start at an invalid point
    # and the test would be about that rather than about the optimiser.
    add!(:DRIFT, [-0.5;;], [1;;], "-log1p_exp(param[1])")
    add!(:JAx, [-0.5;;], [1;;], "-log1p_exp(param[1])")
    add!(:CINT, reshape([0.0], :, 1))
    add!(:DIFFUSION, [0.4;;])
    add!(:LAMBDA, [1.0;;])
    add!(:Jy, [1.0;;])
    add!(:MANIFESTMEANS, reshape([0.0], :, 1), reshape([2], :, 1))
    add!(:MANIFESTVAR, [0.3;;])
    add!(:T0VAR, [0.6;;])
    add!(:T0MEANS, reshape([0.0], :, 1))
    add!(:PARS, reshape([0.0], :, 1))
    n = length(matrices)
    kwargs = Dict{Symbol,Any}()
    if !isempty(manifesttype)
        kwargs[:manifesttype] = manifesttype
        kwargs[:ncategories] = ncategories
    end
    return ContinuousTimeSEM.ekf_from_columns(matrices, rows, cols, parnumber,
        values, transforms, fill("", n), fill("", n), fill("", n);
        diffusion_state_indices=[1], kwargs...)
end

function _joint_setup(; nsubjects=3, nobs=4, manifesttype=Int[], ncategories=Int[],
    seed=99)
    sp = _joint_params(manifesttype=manifesttype, ncategories=ncategories)
    times = repeat(collect(0.0:1.0:(nobs - 1)), nsubjects)
    starts = [1 + (i - 1) * nobs for i in 1:nsubjects]
    rng = MersenneTwister(seed)
    data = if isempty(manifesttype)
        reshape(randn(rng, nsubjects * nobs) .+ 1.0, 1, :)
    else
        reshape(Float64.(rand(rng, 0:3, nsubjects * nobs)), 1, :)
    end
    objective = ContinuousTimeSEM.ctsem_objective(sp, starts, times, data)
    return ContinuousTimeSEM.ctsem_joint_objective(objective, 2)
end

@testset "the joint objective is the joint density" begin
    joint = _joint_setup()
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    @test ndim == 2 + 3 * 4
    x = [(-0.3 + 0.1 * i) for i in 1:ndim]
    theta = x[1:2]
    z = x[3:end]
    @test joint(x) ≈ ContinuousTimeSEM.ctsem_joint_loglikelihood(
        joint.objective, theta, z)
end

@testset "the blocked gradient is the dense one" begin
    joint = _joint_setup()
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    x = [sin(Float64(i)) * 0.5 for i in 1:ndim]
    result = ContinuousTimeSEM.ctsem_evaluate(joint, x; contributions=true)
    @test result.value ≈ joint(x)
    dense = ForwardDiff.gradient(joint, x)
    @test result.gradient ≈ dense
    # And against a difference, so that a shared error in both AD paths would
    # still be caught.
    step = 1e-6
    for i in (1, 2, 3, ndim)
        up = copy(x); up[i] += step
        down = copy(x); down[i] -= step
        @test (joint(up) - joint(down)) / (2 * step) ≈ result.gradient[i] atol = 1e-5
    end
    # Subject contributions sum to the density less the parameter prior, which
    # is zero here.
    @test sum(result.subject_loglik) ≈ result.value
end

@testset "the blocked Hessian is the dense one" begin
    joint = _joint_setup(nsubjects=2, nobs=3)
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    x = [cos(Float64(i)) * 0.4 for i in 1:ndim]
    full = ContinuousTimeSEM.ctsem_joint_hessian(joint, x; profile=false)
    dense = ForwardDiff.hessian(joint, x)
    dense = (dense .+ transpose(dense)) ./ 2
    @test full ≈ dense atol = 1e-8

    # The profiled form is the Schur complement of that same matrix, formed
    # here from the dense one rather than block by block.
    npar = 2
    Hpp = dense[1:npar, 1:npar]
    Hpz = dense[1:npar, (npar + 1):end]
    Hzz = dense[(npar + 1):end, (npar + 1):end]
    expected = Hpp - Hpz * (Hzz \ transpose(Hpz))
    profiled = ContinuousTimeSEM.ctsem_joint_hessian(joint, x; profile=true)
    @test profiled ≈ expected atol = 1e-7
    # And it is genuinely different from the corner of the joint matrix, which
    # is the mistake it exists to avoid: with the states held fixed the
    # curvature is larger, so the interval it would give is narrower.
    @test !isapprox(profiled, Hpp; atol=1e-3)
end

@testset "the optimiser drives the joint target" begin
    joint = _joint_setup(nsubjects=4, nobs=5)
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(joint)
    start = zeros(ndim)
    result = ContinuousTimeSEM.ctsem_optimize(joint, start; maxiter=200,
        progress=false, verbose=false)
    @test isfinite(result.maximum_loglik)
    @test result.maximum_loglik > joint(start)
    @test length(result.minimizer) == ndim
    @test result.gradient_norm < 1e-4
    # Saturation is judged on the population block alone: a large innovation is
    # an unusual trajectory, not a transform pinned at its floating-point
    # limit, and reading it as one would report a good fit as failed.
    @test !result.saturated
end

@testset "a count model fits over the joint density" begin
    # The mixed case the marginal filter approximates and this one does not,
    # on data generated from the model itself so that the parameters are
    # identified and the answer is checkable. Four subjects of five counts --
    # the first size tried here -- left the drift unidentified: the optimiser
    # ran its transform to a raw value of -29, where a drift of zero is flat to
    # machine precision, and the profiled curvature in that direction came back
    # at 7e-4. The engine reported that correctly as not converged; the test
    # was simply asking a question the data could not answer.
    generator = _joint_setup(nsubjects=12, nobs=8,
        manifesttype=[ContinuousTimeSEM.CTSEM_OBS_COUNT], ncategories=[0])
    ndim = ContinuousTimeSEM.ctsem_joint_dimension(generator)
    # -log1p_exp(-0.4328) is -0.5, and the manifest mean is a log rate.
    truth = [-0.4328, 1.0]
    rng = MersenneTwister(4242)
    drawn = ContinuousTimeSEM.ctsem_generate_states(generator.objective, truth,
        randn(rng, ndim - 2), randn(rng, 1, 96))
    @test all(v -> v >= 0 && v == round(v), drawn.Y)

    sp = _joint_params(manifesttype=[ContinuousTimeSEM.CTSEM_OBS_COUNT],
        ncategories=[0])
    times = repeat(collect(0.0:1.0:7.0), 12)
    starts = [1 + (i - 1) * 8 for i in 1:12]
    objective = ContinuousTimeSEM.ctsem_objective(sp, starts, times, drawn.Y)
    joint = ContinuousTimeSEM.ctsem_joint_objective(objective, 2)

    result = ContinuousTimeSEM.ctsem_optimize(joint, zeros(ndim); maxiter=500,
        progress=false, verbose=false)
    @test isfinite(result.maximum_loglik)
    @test all(isfinite, result.minimizer)
    @test result.converged

    # The manifest mean is a log rate the counts see directly, and the joint
    # mode recovers it.
    @test result.minimizer[2] ≈ truth[2] atol = 0.35

    # The drift is not recovered, and that is the estimator rather than a bug:
    # with the trajectory free to be chosen, less mean reversion means smaller
    # innovations, so the innovation prior pushes the drift toward zero and the
    # joint mode follows it there. Measured here: a raw value of -18.5, where
    # `-log1p_exp` is 9e-9 -- a random walk, against a generating drift of -0.5.
    #
    # This is the bias ctFit warns about at `intoverstates=FALSE` with
    # `optimize=TRUE`, and it is asserted rather than avoided so that a change
    # which quietly removed it would be noticed. Integrating the states out
    # (`intoverstates=TRUE`) or sampling them (`optimize=FALSE`) is what
    # estimates a drift.
    @test result.minimizer[1] < truth[1] - 5

    profiled = ContinuousTimeSEM.ctsem_joint_hessian(joint, result.minimizer)
    @test profiled !== nothing
    @test all(isfinite, profiled)
    # Negative semi-definite at a maximum, which is what makes it usable as an
    # observed information. Semi- and not strictly: the drift direction is flat
    # here, for the reason above, and reporting a curvature of 4e-8 for it is
    # the honest answer rather than a failure -- the identifiability report is
    # what turns that into a statement the user sees.
    @test maximum(eigvals(Symmetric(-profiled))) > 1
    @test minimum(eigvals(Symmetric(-profiled))) > -1e-6
end
