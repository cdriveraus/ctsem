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

# A censored row's partials are written by hand (`_censored_partials!`); forward
# mode through `_binary_moment_dual!`, which they replaced, is the oracle. The
# extras are (lower, upper, sd) and, with a response scale, (half-width,
# power): the cases cover each limit, a value recorded exactly and as an
# interval, an end held at a limit, and a power at one (read, so the general
# transform), near zero (the series branch) and negative.
@testset "censored moments' partials by hand match forward mode" begin
    CT = ContinuousTimeSEM
    cases = [
        ("interior", 0.5, 0.64, 2.5, [0.0, 5.0, 0.6]),
        ("at lower limit", 0.5, 0.64, 0.0, [0.0, 5.0, 0.6]),
        ("at upper limit", 4.5, 0.81, 5.0, [0.0, 5.0, 0.6]),
        ("one-sided, deep", -3.0, 0.5, 0.0, [0.0, Inf, 0.4]),
        ("interval", 2.1, 0.7, 3.0, [0.0, 10.0, 0.8, 0.5, 0.7]),
        ("interval, power one", 2.1, 0.7, 3.0, [0.0, 10.0, 0.8, 0.5, 1.0]),
        ("interval, power near zero", 1.0, 0.7, 3.0, [0.0, 10.0, 0.8, 0.5, 1e-9]),
        ("interval, negative power", 1.0, 0.4, 2.0, [0.0, 10.0, 0.5, 0.5, -0.3]),
        ("interval, end held at the floor", 0.3, 0.6, 0.2, [0.0, 10.0, 0.8, 0.5, 0.7]),
        ("at the floor, with a scale", 0.3, 0.6, 0.0, [0.0, 10.0, 0.8, 0.5, 0.7]),
        ("at the ceiling, with a scale", 8.5, 2.0, 10.0, [0.0, 10.0, 0.8, 0.5, 0.7]),
        ("narrow interval far out", -2.0, 0.3, 6.0, [0.0, 10.0, 0.3, 0.05, 0.9]),
        ("exact value, with a power", 2.1, 0.7, 3.0, [0.0, 10.0, 0.8, 0.0, 0.7]),
    ]
    for (label, a, b, y, extras) in cases
        n = length(extras)
        oracle = zeros(3, 2 + n)
        lz, m, v, ok = CT._binary_moment_dual!(oracle, a, b, y, Float64[], Float64[],
            extras, _CENSORED)
        hlz, hm, hv, J = CT._binary_moment_jacobian(a, b, y, Float64[], Float64[],
            extras, _CENSORED)
        @testset "$label" begin
            @test ok
            @test (hlz, hm, hv) == (lz, m, v) ||
                isapprox([hlz, hm, hv], [lz, m, v]; rtol=1e-13, atol=1e-14)
            @test J ≈ oracle rtol = 1e-9 atol = 1e-12
        end
    end
    # Under an outer dual, as inside the Laplace curvature: the hand partials'
    # own derivatives in the outer seed match forward mode's.
    outer(v, p=0.0) = ForwardDiff.Dual{ForwardDiff.Tag{typeof(sin),Float64}}(v, p)
    D = typeof(outer(0.0))
    a, b = outer(2.1, 1.0), outer(0.7, 0.3)
    extras = D[outer(0.0), outer(10.0), outer(0.8, 0.2), outer(0.5), outer(0.7, -0.1)]
    oracle = zeros(D, 3, 7)
    CT._binary_moment_dual!(oracle, a, b, 3.0, Float64[], Float64[], extras, _CENSORED)
    _, _, _, J = CT._binary_moment_jacobian(a, b, 3.0, Float64[], Float64[], extras,
        _CENSORED)
    @test ForwardDiff.value.(J) ≈ ForwardDiff.value.(oracle) rtol = 1e-9 atol = 1e-12
    @test ForwardDiff.partials.(J, 1) ≈ ForwardDiff.partials.(oracle, 1) rtol = 1e-8 atol = 1e-11
end

# An asymptote item's partials are carried by hand through its mixture
# (`_asymptote_partials!`); forward mode, which they replaced, is the oracle:
# each response, asymptotes near their bounds and crossed, a wide and a
# narrow prior, a prior far below and above the item.
@testset "asymptote moments' partials by hand match forward mode" begin
    CT = ContinuousTimeSEM
    nodes, weights = CT._gauss_hermite(CT._CTSEM_BINARY_NODES[])
    cases = [
        ("guess, correct", 0.3, 0.8, 1.0, [0.2, 1.0]),
        ("guess, incorrect", 0.3, 0.8, 0.0, [0.2, 1.0]),
        ("both, correct", -0.4, 1.5, 1.0, [0.15, 0.9]),
        ("both, incorrect", -0.4, 1.5, 0.0, [0.15, 0.9]),
        ("near the bounds", 1.2, 0.5, 1.0, [1e-6, 1 - 1e-6]),
        ("crossed", 0.5, 0.7, 1.0, [0.8, 0.3]),
        ("crossed, incorrect", 0.5, 0.7, 0.0, [0.8, 0.3]),
        ("prior far below", -5.0, 2.0, 1.0, [0.2, 0.95]),
        ("prior far above", 4.0, 0.3, 0.0, [0.25, 0.85]),
        ("narrow prior", 0.1, 1e-3, 1.0, [0.2, 0.9]),
    ]
    for (label, a, b, y, extras) in cases
        oracle = zeros(3, 4)
        lz, m, v, ok = CT._binary_moment_dual!(oracle, a, b, y, nodes, weights,
            extras, CT.CTSEM_OBS_BINARY)
        hlz, hm, hv, J = CT._binary_moment_jacobian(a, b, y, nodes, weights, extras,
            CT.CTSEM_OBS_BINARY)
        @testset "$label" begin
            @test ok
            @test isapprox([hlz, hm, hv], [lz, m, v]; rtol=1e-13, atol=1e-14)
            @test J ≈ oracle rtol = 1e-8 atol = 1e-11
        end
    end
    # Under an outer dual, as inside the Laplace curvature.
    outer(v, p=0.0) = ForwardDiff.Dual{ForwardDiff.Tag{typeof(sin),Float64}}(v, p)
    D = typeof(outer(0.0))
    a, b = outer(-0.4, 1.0), outer(1.5, 0.3)
    extras = D[outer(0.15, 0.2), outer(0.9, -0.1)]
    oracle = zeros(D, 3, 4)
    CT._binary_moment_dual!(oracle, a, b, 1.0, nodes, weights, extras,
        CT.CTSEM_OBS_BINARY)
    _, _, _, J = CT._binary_moment_jacobian(a, b, 1.0, nodes, weights, extras,
        CT.CTSEM_OBS_BINARY)
    @test ForwardDiff.value.(J) ≈ ForwardDiff.value.(oracle) rtol = 1e-8 atol = 1e-11
    @test ForwardDiff.partials.(J, 1) ≈ ForwardDiff.partials.(oracle, 1) rtol = 1e-7 atol = 1e-10
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
            #
            # Asserted the way the solve stops -- on the objective one more
            # Newton step is predicted to gain -- and not on the raw gradient,
            # which is what this line used to do. A gradient bound asks for a
            # number of digits that depends on the curvature, so the same mode
            # passes it at one prior sd and fails it at another: at `atol =
            # 1e-8` three of the seven cases here failed on gradients of
            # around 5e-7, at which the mode is placed to better than 1e-6 and
            # the quadrature it positions is nowhere near that accurate.
            gradient = -offset / s2 + score
            @test ContinuousTimeSEM._newton_gain(gradient,
                gradient / curvature) <=
                ContinuousTimeSEM._mode_tolerance(
                    ContinuousTimeSEM._mode_objective(etabar, 1 / s2, offset,
                        y, thresholds, kind))
            # And the same thing said as a distance, which is what a reader
            # wants to know and which the rule above does not state: one more
            # Newton step would move the mode by less than this.
            @test abs(gradient / curvature) < 1e-5
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
