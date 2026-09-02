# The EKF categorical quadrature kernel: `_binary_moments`, `_binary_mode` and
# `_ekf_binary_update!`, binary_measurement.jl.
#
# These are the functions that actually run during a fit with a categorical
# indicator -- they compute the adaptive Gauss-Hermite moments of the linear
# predictor and lift them back to the state -- and none of the three was
# called from any file in this suite. test_quadrature.jl exercises a different
# quadrature user (the Laplace random-effects integral), and
# test_state_sampling.jl exercises `_category_loglikelihood` at a known,
# sampled state, where eta is not uncertain and no integration runs.
#
# The checks here are against a dense reference computed independently in this
# file: plain trapezoidal integration over a wide, fine grid, with no adaptive
# node placement and nothing shared with the adaptive Gauss-Hermite rule under
# test. All four categorical kinds go through the same shared kernel (`kind`
# is an explicit integer, see binary_measurement.jl:135-144), so all four are
# covered here.

using LinearAlgebra

# E[1], E[eta], E[eta^2] under N(eta; etabar, s^2) * P(y | eta), by trapezoidal
# integration in log space (so a likelihood that is a large negative number
# stays representable, exactly the concern `_category_loglikelihood`'s own
# docstring raises about the probability form).
function _dense_categorical_moments(etabar, s, y, thresholds, kind;
    halfwidth=40.0, n=200001)
    lo = etabar - halfwidth * s
    hi = etabar + halfwidth * s
    h = (hi - lo) / (n - 1)
    logprior(eta) = -(eta - etabar)^2 / (2 * s^2) - log(sqrt(2 * pi) * s)
    logvals = Vector{Float64}(undef, n)
    m = -Inf
    @inbounds for i in 1:n
        eta = lo + (i - 1) * h
        lv = logprior(eta) +
            ContinuousTimeSEM._category_loglikelihood(eta, y, thresholds, kind)
        logvals[i] = lv
        lv > m && (m = lv)
    end
    w = exp.(logvals .- m)
    etas = [lo + (i - 1) * h for i in 1:n]
    # Trapezoidal rule: the sum minus half the two endpoints, times the step.
    Z = h * (sum(w) - 0.5 * (w[1] + w[end]))
    M1 = h * (sum(w .* etas) - 0.5 * (w[1] * etas[1] + w[end] * etas[end]))
    M2 = h * (sum(w .* etas .^ 2) -
        0.5 * (w[1] * etas[1]^2 + w[end] * etas[end]^2))
    mean = M1 / Z
    return (logZ=log(Z) + m, mean=mean, variance=M2 / Z - mean^2)
end

# Same fixed one-latent, one-manifest model `test_workspace_and_ekf.jl` uses,
# with `manifesttype` added so a workspace built from it carries the field
# `_ekf_binary_update!` and its callers read. LAMBDA = 1, MANIFESTMEANS = 0
# means the linear predictor is just the state itself, which keeps the
# dense-reference comparison a one-line correspondence.
function _one_dim_categorical_params(kind::Int)
    mats = Pair[:DRIFT => [-0.5;;], :JAx => [-0.5;;],
        :CINT => reshape([0.0], :, 1), :DIFFUSION => [0.2;;],
        :LAMBDA => [1.0;;], :Jy => [1.0;;],
        :MANIFESTMEANS => reshape([0.0], :, 1), :MANIFESTVAR => [0.0;;],
        :T0VAR => [0.3;;], :T0MEANS => reshape([0.0], :, 1),
        :PARS => reshape([0.0], :, 1)]
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
        AbstractFloat[values...], Int[], Int[], Int[], [1],
        true, [kind], [0])
end

const _BINARY = ContinuousTimeSEM.CTSEM_OBS_BINARY
const _ORDINAL = ContinuousTimeSEM.CTSEM_OBS_ORDINAL
const _COUNT = ContinuousTimeSEM.CTSEM_OBS_COUNT
const _CENSORED = ContinuousTimeSEM.CTSEM_OBS_CENSORED

# (label, etabar, s, y, thresholds, kind). Censored carries (lower, upper, sd)
# in the thresholds slot, matching `_censor_limits`; one case sits at each
# limit and one strictly inside, since the kernel branches on that.
const _MOMENT_CASES = [
    ("binary y=1", 0.4, 1.3, 1.0, (), _BINARY),
    ("binary y=0", -0.6, 0.7, 0.0, (), _BINARY),
    ("ordinal k=2 of 4", 0.2, 0.9, 2, [-1.0, 0.4, 1.9], _ORDINAL),
    ("ordinal k=1 of 4 (end category)", 0.2, 0.9, 1, [-1.0, 0.4, 1.9], _ORDINAL),
    ("count y=3", 1.1, 0.7, 3, (), _COUNT),
    ("count y=0", -0.3, 1.0, 0, (), _COUNT),
    ("censored, interior", 0.5, 0.8, 2.5, [0.0, 5.0, 0.6], _CENSORED),
    ("censored, at lower limit", 0.5, 0.8, 0.0, [0.0, 5.0, 0.6], _CENSORED),
    ("censored, at upper limit", 4.5, 0.9, 5.0, [0.0, 5.0, 0.6], _CENSORED),
]

@testset "_binary_moments matches a dense reference quadrature" begin
    nodes, weights = ContinuousTimeSEM._gauss_hermite(
        ContinuousTimeSEM._CTSEM_BINARY_NODES[])
    for (label, etabar, s, y, thresholds, kind) in _MOMENT_CASES
        dense = _dense_categorical_moments(etabar, s, y, thresholds, kind)
        logZ, offset, variance = ContinuousTimeSEM._binary_moments(etabar, s,
            y, nodes, weights, thresholds, kind)
        @testset "$label" begin
            @test logZ ≈ dense.logZ atol = 1e-6
            @test (etabar + offset) ≈ dense.mean atol = 1e-6
            @test variance ≈ dense.variance atol = 1e-6
        end
    end
end

@testset "_binary_mode solves the scalar posterior's first-order condition" begin
    for (label, etabar, s, y, thresholds, kind) in _MOMENT_CASES
        s2 = s^2
        offset, curvature = ContinuousTimeSEM._binary_mode(etabar, s2, y,
            thresholds, kind)
        score, information = ContinuousTimeSEM._category_score(etabar + offset,
            y, thresholds, kind)
        @testset "$label" begin
            # d/deta [ -(eta-etabar)^2/(2s^2) + log P(y|eta) ] = 0 at the mode:
            # the prior's score, -(eta-etabar)/s^2 = -offset/s^2, plus the
            # observation's score, must cancel.
            @test (-offset / s2 + score) ≈ 0.0 atol = 1e-8
            # Curvature is the prior's precision plus the observed information,
            # by construction (`_binary_mode`'s own docstring): checked as an
            # identity independent of the dense reference above.
            @test curvature ≈ (1 / s2 + information) atol = 1e-10
        end
    end
end

@testset "_ekf_binary_update! moves a one-state workspace to the posterior moments" begin
    # Not the scalar moments in isolation, as above, but the full state-space
    # mechanics: c = P*lambda, s^2 = lambda'*P*lambda, the mean shift and the
    # covariance shrink `1 - V[eta|y]/s^2`, written into the workspace in
    # place. One latent state with LAMBDA = 1 collapses that exactly onto the
    # scalar dense reference, so this checks the matrix code, not a different
    # answer. Ordinal is not repeated here: `_ordinal_thresholds!` is what
    # supplies its thresholds from the workspace in the real filter, and that
    # plumbing belongs with the row-application tests, not this scalar check.
    for (label, etabar, s, y, thresholds, kind) in (
        ("binary", 0.4, 1.3, 1.0, (), _BINARY),
        ("count", 1.1, 0.7, 3.0, (), _COUNT),
        ("censored, interior", 0.5, 0.8, 2.5, [0.0, 5.0, 0.6], _CENSORED),
        ("censored, at limit", 0.5, 0.8, 0.0, [0.0, 5.0, 0.6], _CENSORED),
    )
        sp = _one_dim_categorical_params(kind)
        ws = ContinuousTimeSEM._init_continuous_ekf_workspace(Float64, sp)
        ws.state[1] = etabar
        ws.P_predict.data[1, 1] = s^2
        logZ = ContinuousTimeSEM._ekf_binary_update!(ws, [1.0], 0.0, y, 1,
            thresholds, kind)
        dense = _dense_categorical_moments(etabar, s, y, thresholds, kind)
        @testset "$label" begin
            @test logZ ≈ dense.logZ atol = 1e-6
            @test ws.state[1] ≈ dense.mean atol = 1e-6
            @test ws.P_predict.data[1, 1] ≈ dense.variance atol = 1e-6
            # The update mirrors the lower triangle to the upper one; trivial
            # at n=1, checked anyway since it is part of what the function does.
            @test ws.P_predict.data[1, 1] == ws.P_predict[1, 1]
        end
    end
end
