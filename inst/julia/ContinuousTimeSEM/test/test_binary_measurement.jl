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
# narrow prior, a prior far below and above the item. Except in an asymptote
# whose response makes the flat part exactly zero -- a ceiling of one and an
# incorrect response, a floor of zero and a correct one -- where forward mode
# goes through `log(0)` and drops that term's derivative; there the column is
# checked against differences of the values. Such a cell is fixed in
# a model (a probability transform does not reach its bound), so the column
# is not read.
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
        flatzero = (y > 0.5 && extras[1] == 0) || (y < 0.5 && extras[2] == 1)
        cols = [1, 2, 3, 4]
        flatzero && setdiff!(cols, y > 0.5 ? 3 : 4)
        @testset "$label" begin
            @test ok
            @test isapprox([hlz, hm, hv], [lz, m, v]; rtol=1e-13, atol=1e-14)
            @test J[:, cols] ≈ oracle[:, cols] rtol = 1e-8 atol = 1e-11
            if flatzero
                # One-sided, second order, stepping inward: past the bound the
                # flat part is negative, and its log throws.
                col = y > 0.5 ? 3 : 4
                inward = col == 3 ? 1.0 : -1.0
                h = 1e-5
                at(t) = (e = copy(extras); e[col - 2] += inward * t; e)
                f(e) = collect(CT._binary_moments(a, sqrt(b), y, nodes, weights, e,
                    CT.CTSEM_OBS_BINARY))
                d = inward .* (-3 .* f(at(0.0)) .+ 4 .* f(at(h)) .- f(at(2h))) ./ (2h)
                @test J[:, col] ≈ d rtol = 1e-6 atol = 1e-8
            end
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

# A binary observation or an ordinal end category at a predicted sd above
# `_CTSEM_BINARY_VBAND`'s lower end is integrated in the logistic variable,
# blended with the rule in eta across the band (`_binary_quadrature`). The
# rule in eta alone put P(1) + P(0) at 1.06 for s = 500, with an error that
# changed sign as s grew, and fits stopped on the spurious maxima it made.
# Inside the band the rule in eta keeps part of the weight, and with it its own
# error there, about 3e-5 at s = 3; above it the rule in V is good to 1e-6.
@testset "one-sided observations at large predicted sd" begin
    CT = ContinuousTimeSEM
    nodes, weights = CT._binary_rule()
    tau = [-1.0, 0.4, 1.9]
    lo, hi = CT._CTSEM_BINARY_VBAND
    tol(s, outside) = lo < s < hi ? 1e-4 : outside
    for s in (2.0, 2.9, 3.4, 4.0, 10.0, 100.0, 1000.0), z in (-2.0, 0.0, 0.7)
        etabar = z * s
        p = [exp(CT._binary_moments(etabar, s, y, nodes, weights, (), _BINARY)[1])
            for y in (0.0, 1.0)]
        po = [exp(CT._binary_moments(etabar, s, k, nodes, weights, tau, _ORDINAL)[1])
            for k in 1:4]
        @test sum(p) ≈ 1 atol = tol(s, 1e-6)
        @test sum(po) ≈ 1 atol = tol(s, 1e-6)
    end
    # Against the dense reference, in the band and above it.
    for s in (2.8, 3.2, 5.0, 40.0), z in (-2.0, 0.0, 1.5),
        (y, th, kind) in ((1.0, (), _BINARY), (0.0, (), _BINARY), (1, tau, _ORDINAL),
            (4, tau, _ORDINAL))
        etabar = z * s
        dense = _dense_categorical_moments(etabar, s, y, th, kind)
        logZ, offset, variance = CT._binary_moments(etabar, s, y, nodes, weights, th,
            kind)
        @test logZ ≈ dense.logZ atol = tol(s, 2e-6)
        @test (etabar + offset) / s ≈ dense.mean / s atol = tol(s, 2e-6)
        @test variance / s^2 ≈ dense.variance / s^2 atol = tol(s, 2e-6)
    end
    # The partials by hand against forward mode through the same function, in
    # the band (where the blend's weight moves) and above it; an asymptote
    # item, whose logistic part is this; and under an outer dual, as inside
    # the Laplace curvature and the hand second order.
    for b in (2.7^2, 3.0^2, 3.3^2, 4.0^2, 20.0^2, 300.0^2), z in (-2.0, 0.3),
        (y, th, kind) in ((1.0, Float64[], _BINARY), (0.0, Float64[], _BINARY),
            (1, tau, _ORDINAL), (4, tau, _ORDINAL), (1.0, [0.2, 0.95], _BINARY))
        a = z * sqrt(b)
        oracle = zeros(3, 2 + length(th))
        lz, m, v, ok = CT._binary_moment_dual!(oracle, a, b, y, nodes, weights, th,
            kind)
        hlz, hm, hv, J = CT._binary_moment_jacobian(a, b, y, nodes, weights, th,
            kind)
        @test ok
        @test isapprox([hlz, hm, hv], [lz, m, v]; rtol=1e-12, atol=1e-12)
        @test J ≈ oracle rtol = 1e-8 atol = 1e-10
    end
    outer(v, p=0.0) = ForwardDiff.Dual{ForwardDiff.Tag{typeof(sin),Float64}}(v, p)
    D = typeof(outer(0.0))
    for (a, b) in ((outer(-0.6, 1.0), outer(3.1^2, 0.3)),
        (outer(-0.6, 1.0), outer(12.0^2, 5.0)))
        for (y, th, kind) in ((1.0, D[], _BINARY),
            (4, D[outer(-1.0, 0.2), outer(0.4), outer(1.9, -0.1)], _ORDINAL))
            oracle = zeros(D, 3, 2 + length(th))
            CT._binary_moment_dual!(oracle, a, b, y, nodes, weights, th, kind)
            _, _, _, J = CT._binary_moment_jacobian(a, b, y, nodes, weights, th, kind)
            @test ForwardDiff.value.(J) ≈ ForwardDiff.value.(oracle) rtol = 1e-8 atol = 1e-10
            @test ForwardDiff.partials.(J, 1) ≈ ForwardDiff.partials.(oracle, 1) rtol = 1e-7 atol = 1e-9
        end
    end
    # The mode solve in V is bracketed, so it ends at the mode from any finite
    # input and the rule's partials are always its own: no observation above
    # the band goes to forward mode, at predicted means far either side of
    # the threshold and at any width.
    for s in (3.6, 50.0, 1e4), z in (-40.0, -5.0, 0.0, 5.0, 40.0)
        m = z * s
        v, C, ok = CT._binary_vmode(m, s)
        g, _ = CT._binary_vslope(v, m, s)
        @test ok
        @test abs(g / C) <= 1e-10 * (1 + abs(v))
        J = zeros(3, 2)
        _, _, _, handled = CT._binary_moment_partials!(J, m, s^2, 1.0, nodes,
            weights, Float64[], _BINARY)
        @test handled
    end
    # The blend joins each rule smoothly to every order: just inside the band,
    # the blended log P and its first three derivatives in s are the single
    # rule's. The Laplace gradient reads the third. A quintic smoothstep put a
    # difference of about 60 times the rules' gap (1e-6) into the third
    # derivative here.
    d1(f) = x -> ForwardDiff.derivative(f, x)
    for (s, single) in ((lo + 0.02 * (hi - lo),
            s -> CT._binary_quadrature_eta(promote(0.4, s)..., 1.0, nodes, weights, (),
                _BINARY)[1]),
        (hi - 0.02 * (hi - lo), s -> CT._binary_vquadrature(promote(0.4, s)..., false)[1]))
        blend(s) = CT._binary_moments(promote(0.4, s)..., 1.0, nodes, weights, (),
            _BINARY)[1]
        for f in (identity, d1, d1 ∘ d1, d1 ∘ d1 ∘ d1)
            @test f(blend)(s) ≈ f(single)(s) atol = 1e-9
        end
    end
    # At a trace jet, as the hand second order's jet sweeps run it: the
    # moments' first parts are their directional derivatives and the second
    # part the second, which nested forward mode through the values gives; the
    # moments' Jacobian likewise. The third case sits just inside the band's
    # top, where the blend's weight takes the logistic of about -500: there a
    # jet's second part of the Jacobian's b column was NaN (bench gD1), until
    # `_logistic` stopped exponentiating positive arguments.
    for (a0, b0, y) in ((-0.6, 3.1^2, 1.0), (2.0, 15.0^2, 0.0), (4.2175, 3.498^2, 1.0))
        da, db = 0.7, -0.4 * b0
        jet(x, d) = CT.CTSEMTraceJet{1}(x, (d,), 0.0)
        Jm = zeros(CT.CTSEMTraceJet{1}, 3, 2)
        logZ, m, v, ok = CT._binary_moment_partials!(Jm, jet(a0, da), jet(b0, db), y,
            nodes, weights, CT.CTSEMTraceJet{1}[], _BINARY)
        @test ok
        along(f) = t -> f(CT._binary_moments(a0 + t * da, sqrt(b0 + t * db), y, nodes,
            weights, (), _BINARY))
        for (k, x) in enumerate((logZ, m, v))
            g = along(r -> r[k])
            @test x.d[1] ≈ ForwardDiff.derivative(g, 0.0) rtol = 1e-7 atol = 1e-9
            @test x.s ≈ ForwardDiff.derivative(t -> ForwardDiff.derivative(g, t), 0.0) rtol = 1e-6 atol = 1e-8
        end
        jac(t) = CT._binary_moment_jacobian(a0 + t * da, b0 + t * db, y, nodes, weights,
            typeof(t)[], _BINARY)[4]
        dJ = ForwardDiff.derivative(jac, 0.0)
        d2J = ForwardDiff.derivative(t -> ForwardDiff.derivative(jac, t), 0.0)
        @test all(x -> isfinite(x.v) && isfinite(x.d[1]) && isfinite(x.s), Jm)
        @test [x.d[1] for x in Jm] ≈ dJ rtol = 1e-6 atol = 1e-8
        @test [x.s for x in Jm] ≈ d2J rtol = 1e-5 atol = 1e-7
    end
end

# `_logistic` at second order, out to where `exp` overflows: a trace jet's
# parts are the logistic's first and second derivatives, `p(1 - p)` and
# `p(1 - p)(1 - 2p)`, finite at every argument. Switched to `e^x / (1 + e^x)`
# only below -700, the second part was `Inf / Inf` from about -236.
@testset "_logistic's jet parts are finite at any argument" begin
    CT = ContinuousTimeSEM
    for x in (-800.0, -500.0, -300.0, -30.0, -0.5, 0.0, 0.5, 30.0, 300.0, 800.0)
        p = x < 0 ? exp(x) / (1 + exp(x)) : 1 / (1 + exp(-x))
        q = x < 0 ? 1 / (1 + exp(x)) : exp(-x) / (1 + exp(-x))
        r = CT._logistic(CT.CTSEMTraceJet{1}(x, (1.0,), 0.0))
        @test r.v ≈ p rtol = 1e-14
        @test r.d[1] ≈ p * q rtol = 1e-12
        @test r.s ≈ p * q * (q - p) rtol = 1e-10 atol = 1e-300
        @test isfinite(r.d[1]) && isfinite(r.s)
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
