"""
A choice and its response time as one observation of a latent drift: the
Wiener first-passage density of the drift diffusion model (manifesttype 5,
`CTSEM_OBS_WIENER`).

# The observation

The data value is the signed response time: `y = RT` for a response at the
upper boundary, `y = -RT` at the lower. With diffusion coefficient 1, boundary
separation `a`, relative start `w` and nondecision time `tau`, the decision
time is `t = RT - tau` and, for drift `v` (Navarro & Fuss, 2009),

    f(t, boundary | v) = a^-2 exp(v x - v^2 t / 2) g(t / a^2, w*)

with `x = a (up - w)` the end point relative to the start, `w*` the start's
distance from the boundary reached (`1 - w` at the upper, `w` at the lower)
and `g` the density for drift 0 and boundary 1. As a function of `v` this is

    f = N(x / t; v, 1 / t) exp(h),    h = log a - 2 log t + log S(t / a^2, w*)

where `S` is `g` with its Gaussian factor taken out (`_wiener_logS`), and the
`x^2 / 2t` terms of the kernel and of `g` cancel. So a trial is an exact linear
Gaussian observation of the drift: pseudo-observation `z = x / t` with variance
`q = 1 / t`, plus a term `h` that does not involve the drift.

The drift is the row's linear predictor `eta = lambda'x + mu`. Trial to trial
drift variability `sv` is a Gaussian term added to it, so it joins the
predictor's variance `b` exactly as a count's dispersion does
(`_count_dispersion`), and the moments here never read it. Under the predicted
`eta ~ N(ebar, b)` everything is closed form, as for a censored row
(`_censored_moments`):

    logZ = h + log N(z; ebar, q + b)
    E[eta | y] - ebar = (z - ebar) b / (q + b)
    V[eta | y] = b q / (q + b)

and the update is the exact Kalman update.

# The extras

`(sv, a, tau, w)`: `sv` is MANIFESTVAR's diagonal entry, read as a standard
deviation as a count's dispersion is, and `a`, `tau`, `w` are the row's WIENER
cells. The moments read the last three (`_extras_read`), which keeps the
partials within the five the sweeps' dual carries (`_TanCatDual`).

# The edge

Where `tau` passes a response time the density is zero. A unit's inner Newton
starts at the population values, which lie past the edge for anyone whose
fastest response is below the population `tau`, so the likelihood needs a
finite continuation there. The decision time is clamped smoothly at
`_WIENER_T0` (exact to 1e-14 s for decision times above about 30 ms), and a
trial past the edge is charged `slope d + quad d^2` nats for `d` seconds past
it. The charge is at least the Jacobian the density of the response time needs
there, so the density integrates to at most one for any `tau`. Three other
shapes failed on real data: no Jacobian (an unbounded maximum as `a` goes to
zero, found at `tau = 1e41`), a doubly exponential wall (Newton could not get
back), and the Jacobian alone (a second maximum that gave up a few fast
trials). The quadratic term lands a step from past the edge within
`slope / (2 quad)` = 5 ms of it.

# The series

`S(u, w)` is the small-time series where `u` is small and the large-time series
where it is large, blended on `log u` by a logistic weight; with six terms each
both are accurate across the blend (within 2e-5 nats of RWiener's density for
decision times 0.03-8 s and boundaries 0.6-2.5). The term counts are fixed: a
count that changed with the parameters would make the objective jump.
"""

const _WIENER_T0 = 0.005
const _WIENER_EDGE = 1000.0
const _WIENER_SLOPE = 1e5
const _WIENER_QUAD = 1e7
const _WIENER_KSMALL = 6
const _WIENER_KLARGE = 6
const _WIENER_BLEND = 10.0
const _WIENER_LARGE_C = 0.5 * log(2 * pi) + log(pi)
const _WIENER_HALFLOG2PI = 0.5 * log(2 * pi)

"""`log(1 + e^x)` without overflow, at any number type the engine runs."""
@inline _wiener_softplus(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

"""A series sum held at `floor` where it is not positive: there the series is
outside its range and the blend gives it no weight."""
@inline _wiener_floorpos(S, floor::Float64) =
    _primal(S) > floor ? (S, true) : (zero(S) + floor, false)

"""
    _wiener_logS_small(u, ws)

`log S` from the small-time series, `S = sum_k (ws + 2k) exp(-(2k^2 + 2k ws) / u)`,
`k = -K..K`, and its partials in `u` and `ws`. Every exponent is at most zero.
"""
@inline function _wiener_logS_small(u, ws)
    S = zero(u * ws)
    Su = zero(S)
    Sw = zero(S)
    for k in -_WIENER_KSMALL:_WIENER_KSMALL
        if k == 0
            S += ws
            Sw += one(S)
            continue
        end
        c = 2 * k * k + 2 * k * ws
        e = exp(-c / u)
        term = (ws + 2 * k) * e
        S += term
        Su += term * c / (u * u)
        Sw += e - term * (2 * k) / u
    end
    S, ok = _wiener_floorpos(S, 1e-12)
    ok || return (log(S), zero(S), zero(S))
    return (log(S), Su / S, Sw / S)
end

"""
    _wiener_logS_large(u, ws)

`log S` from the large-time series, `log g + log(2 pi u^3) / 2 + ws^2 / 2u` with
`g = pi sum_k k exp(-k^2 pi^2 u / 2) sin(k pi ws)`, `k = 1..K`, and its partials
in `u` and `ws`.
"""
@inline function _wiener_logS_large(u, ws)
    L = zero(u * ws)
    Lu = zero(L)
    Lw = zero(L)
    for k in 1:_WIENER_KLARGE
        r = (k * k - 1) * pi * pi / 2
        e = exp(-r * u)
        s = sin(k * pi * ws)
        c = cos(k * pi * ws)
        L += k * e * s
        Lu -= k * e * s * r
        Lw += k * e * c * k * pi
    end
    L, ok = _wiener_floorpos(L, 1e-300)
    value = _WIENER_LARGE_C + 1.5 * log(u) + ws * ws / (2 * u) - pi * pi / 2 * u +
        log(L)
    du = 1.5 / u - ws * ws / (2 * u * u) - pi * pi / 2
    dw = ws / u
    ok && (du += Lu / L; dw += Lw / L)
    return (value, du, dw)
end

"""
    _wiener_logS(u, ws)

`log S(u, ws)` blended between the two series on `log u`, and its partials in
`u` and `ws`.
"""
@inline function _wiener_logS(u, ws)
    ls, lsu, lsw = _wiener_logS_small(u, ws)
    ll, llu, llw = _wiener_logS_large(u, ws)
    ω = _logistic(_WIENER_BLEND * log(u))
    dω = _WIENER_BLEND * ω * (1 - ω) / u
    value = (1 - ω) * ls + ω * ll
    du = (1 - ω) * lsu + ω * llu + dω * (ll - ls)
    dw = (1 - ω) * lsw + ω * llw
    return (value, du, dw)
end

"""
    _wiener_kernel(y, a, tau, w)

A trial's Gaussian kernel in the drift and the term free of it, with their
partials in the row's parameters: `(z, q, h, z_a, z_tau, z_w, q_tau, h_a,
h_tau, h_w)` for pseudo-observation `z = x / t`, kernel variance `q = 1 / t`
and `h` (the module's header). `y` is the signed response time, data, so which
boundary was reached is decided on its value.
"""
@inline function _wiener_kernel(y::Real, a, tau, w)
    up = y > 0
    sgn = up ? 1 : -1
    gap = abs(y) - tau - _WIENER_T0
    kg = _WIENER_EDGE * gap
    t = _WIENER_T0 + _wiener_softplus(kg) / _WIENER_EDGE
    t_tau = -_logistic(kg)
    d = _wiener_softplus(-kg) / _WIENER_EDGE
    charge = _WIENER_SLOPE * d + _WIENER_QUAD * d * d
    charge_tau = (_WIENER_SLOPE + 2 * _WIENER_QUAD * d) * _logistic(-kg)
    ws = up ? 1 - w : w
    ws_w = up ? -1 : 1
    u = t / (a * a)
    lS, lS_u, lS_w = _wiener_logS(u, ws)
    z = sgn * a * ws / t
    q = 1 / t
    h = log(a) - 2 * log(t) + lS - charge
    z_a = sgn * ws / t
    z_w = sgn * a * ws_w / t
    z_tau = -z * t_tau / t
    q_tau = -t_tau / (t * t)
    h_a = 1 / a - 2 * lS_u * u / a
    h_w = lS_w * ws_w
    h_tau = -2 * t_tau / t + lS_u * t_tau / (a * a) - charge_tau
    return (z, q, h, z_a, z_tau, z_w, q_tau, h_a, h_tau, h_w)
end

"""The `(sv, a, tau, w)` a Wiener row carries (`_ordinal_thresholds!`), as
they are, for the reason `_censor_limits` gives."""
@inline _wiener_extras(extras) = (extras[1], extras[2], extras[3], extras[4])

"""
    _wiener_loglikelihood(eta, y, extras)

A trial's log density given the drift's predictor `eta`, its variability `sv`
integrated: `h + log N(z; eta, q + sv^2)`.
"""
@inline function _wiener_loglikelihood(eta, y::Real, extras)
    sv, a, tau, w = _wiener_extras(extras)
    z, q, h = _wiener_kernel(y, a, tau, w)
    V = q + sv * sv
    dev = z - eta
    return h - dev * dev / (2 * V) - log(V) / 2 - _WIENER_HALFLOG2PI
end

"""`(score, information)` of `_wiener_loglikelihood` in `eta`."""
@inline function _wiener_score(eta, y::Real, extras)
    sv, a, tau, w = _wiener_extras(extras)
    z, q, _ = _wiener_kernel(y, a, tau, w)
    V = q + sv * sv
    return ((z - eta) / V, inv(V))
end

"""
    _wiener_moments(ebar, s, y, extras)

`(logZ, mean - ebar, variance)` for a Wiener row under `eta ~ N(ebar, s^2)`, in
closed form; `s^2` already holds `sv^2` (`_ekf_binary_update!`).
"""
@inline function _wiener_moments(ebar::T, s::T, y::Real, extras) where {T}
    _, a, tau, w = _wiener_extras(extras)
    z, q, h = _wiener_kernel(y, a, tau, w)
    b = s * s
    V = q + b
    dev = z - ebar
    logZ = h - dev * dev / (2 * V) - log(V) / 2 - _WIENER_HALFLOG2PI
    return (logZ, dev * b / V, b * q / V)
end

"""
    _wiener_partials!(J, ebar, s, y, extras)

`_wiener_moments(ebar, s, y, extras)` with its partials written into `J` by
hand, as `_binary_moment_jacobian` lays them out: rows `logZ`, mean offset and
variance; columns `ebar`, `b = s^2`, then the extras `sv` (not read: zero),
`a`, `tau`, `w`. Generic in `T`, so the hand-coded sweeps take its second order
by running it at their duals and jets. Returns `(logZ, offset, variance, true)`.

Through the kernel's `z`, `q` and `h`, with `V = q + b` and `dev = z - ebar`:
`d logZ = dh - (dev/V) dz + L dq`, `L = (dev^2 - V) / 2V^2`;
`d m = (b/V) dz - (dev b / V^2) dq`; `d v = (b/V)^2 dq`.
"""
function _wiener_partials!(J::Matrix{T}, ebar::T, s::T, y::Real, extras) where {T}
    _, a, tau, w = _wiener_extras(extras)
    z, q, h, z_a, z_tau, z_w, q_tau, h_a, h_tau, h_w = _wiener_kernel(y, a, tau, w)
    b = s * s
    V = q + b
    dev = z - ebar
    logZ = h - dev * dev / (2 * V) - log(V) / 2 - _WIENER_HALFLOG2PI
    L = (dev * dev - V) / (2 * V * V)
    @inbounds begin
        J[1, 1] = dev / V
        J[1, 2] = L
        J[2, 1] = -b / V
        J[2, 2] = dev * q / (V * V)
        J[3, 2] = (q / V)^2
        if size(J, 2) >= 6
            r = b / V
            mq = -dev * b / (V * V)
            J[1, 4] = h_a - (dev / V) * z_a
            J[2, 4] = r * z_a
            J[1, 5] = h_tau - (dev / V) * z_tau + L * q_tau
            J[2, 5] = r * z_tau + mq * q_tau
            J[3, 5] = r * r * q_tau
            J[1, 6] = h_w - (dev / V) * z_w
            J[2, 6] = r * z_w
        end
    end
    return (logZ, dev * b / V, b * q / V, true)
end

"""
    _wiener_upper_probability(v, a, w)

The probability that a diffusion with drift `v` from `w a` reaches the upper
boundary of `[0, a]` first, `(1 - e^{-2 v a w}) / (1 - e^{-2 v a})`, written to
stay finite for drifts of either sign and continuous through zero, where it
is `w`.
"""
@inline function _wiener_upper_probability(v::Real, a::Real, w::Real)
    k = 2 * v * a
    abs(k) < 1e-8 && return w + k * w * (1 - w) / 2
    # For a negative drift, numerator and denominator times e^k, so that no
    # exponential taken is above one.
    return k > 0 ? expm1(-k * w) / expm1(-k) :
        (exp(k * (1 - w)) - exp(k)) / -expm1(k)
end

"""
    _wiener_survival(t, v, a, w)

`P(T > t, lower boundary)` for a diffusion with drift `v` from `w a`: the
large-time series of the density integrated from `t` to infinity,
`(pi / a^2) e^{-v a w} sum_k k sin(k pi w) e^{-l_k t} / l_k`,
`l_k = v^2/2 + k^2 pi^2 / 2a^2`, to 200 terms. Generation only.
"""
function _wiener_survival(t::Float64, v::Float64, a::Float64, w::Float64)
    total = 0.0
    for k in 1:200
        l = v * v / 2 + k * k * pi * pi / (2 * a * a)
        e = -v * a * w - l * t
        term = k * sin(k * pi * w) * exp(min(e, 700.0)) / l
        total += term
        k > 10 && abs(term) < 1e-16 * abs(total) && break
    end
    return pi / (a * a) * total
end

"""
    _wiener_draw(v, a, w, u1, u2)

One first passage, as `(up, t)`: the boundary from `u1` against
`_wiener_upper_probability`, and the decision time from `u2` by inverting that
boundary's survival function (`_wiener_survival`, the upper boundary as the
lower with the drift and start reflected), by bisection to 1e-12 s.
"""
function _wiener_draw(v::Float64, a::Float64, w::Float64, u1::Float64, u2::Float64)
    pup = _wiener_upper_probability(v, a, w)
    up = u1 < pup
    vv, ww = up ? (-v, 1 - w) : (v, w)
    p = up ? pup : 1 - pup
    target = (1 - u2) * p
    lo, hi = 0.0, a * a
    while _wiener_survival(hi, vv, a, ww) > target && hi < 1e6
        lo = hi
        hi *= 2
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        _wiener_survival(mid, vv, a, ww) > target ? (lo = mid) : (hi = mid)
        hi - lo < 1e-12 && break
    end
    return (up, (lo + hi) / 2)
end

"""
    _wiener_choice_moments(ebar, s2, extras)

The predicted choice of a Wiener row whose drift predictor is `N(ebar, s2)`
(`s2` without `sv^2`, which is added here): the probability of the upper
boundary averaged over the drift, by 41 Gauss-Hermite nodes, and the Bernoulli
variance of that. What ctKalman and ctPredict report for the item.
"""
function _wiener_choice_moments(ebar::T, s2::T, extras) where {T}
    sv, a, _, w = _wiener_extras(extras)
    spread = sqrt(max(s2 + sv * sv, zero(T)))
    nodes, weights = _gauss_hermite(41)
    m = zero(T)
    @inbounds for k in eachindex(nodes)
        v = ebar + sqrt(T(2)) * spread * T(nodes[k])
        m += T(weights[k]) / sqrt(T(pi)) * _wiener_upper_probability(v, a, w)
    end
    return (m, m * (1 - m))
end
