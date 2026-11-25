"""
The joint posterior sampled through a learnt map.

# The target

The joint density over the population parameters and the effects
(`ctsem_sample_density!`), on whichever objective it is given. On the hybrid
objective the effects the filter integrates exactly -- an effect on
MANIFESTMEANS, T0MEANS, CINT or TDPREDEFFECT through the affine transform, in a
linear Gaussian model -- are integrated (`_HYBRID_FULL`), so the target is the
exact joint over theta and the effects that remain: fewer dimensions, and none
of those effects' dependence on the others for the sampler to cross. On the
layer it is the joint `ctsem_sample` draws.

# The map

NUTS runs in latent coordinates `(w, z_1..z_N)`, mapped to `(theta, u_1..u_N)`
triangularly:

    theta = theta0 + C w
    u_U   = m_U(theta) + R_U(theta)^-1 z_U          each unit, given theta

A member `s` of unit `U` (one subject) sees its unit's coordinates only through
`b_s = theta[P] + J_s(theta) u_U`, its raw parameters at the positions `P` that
vary at some level, `J_s` holding each level's root at the member's offsets as
`_laplace_member_values!` shifts them. Each member's likelihood is summarised
as a Gaussian in `b_s` (precision `Lam_s`, linear term `h_s`), and for that
summary the unit's conditional given theta is

    Q_U = I + sum_s J_s' Lam_s J_s,
    m_U = Q_U^-1 sum_s J_s' (h_s - Lam_s theta[P] - G_s (theta[R] - thref)),
    R_U = chol(Q_U):

partial non-centring (Papaspiliopoulos, Roberts and Skold 2007), centred where a
unit's data dominate and non-centred where its prior does, continuous in theta,
for any number of levels, members a unit and ranks. `G_s` couples a member's
summary to the free parameters outside `P` its likelihood also reads (`cross`):
those every subject reads and its own sampled TI predictor values, never
another subject's (`_jmap_layout`).
The log Jacobian is `log det C - sum_U log det R_U`, exactly, so the draws are
the target's whatever the summaries are; they decide only how near normal the
latent posterior is.

# Learning

The summaries start from the Laplace fit -- each member's curvature in `b_s` at
its unit's mode, by central differences of its score, which also give `G_s` --
and theta's whitening from the fit's covariance. In warm-up windows, from the
chains' draws pooled over the chains that moved: theta's centre and `C` from
their covariance shrunk toward the last (`_smc_shrunk_cov`), each member's
summary by least squares of its scores on `b_s` in the directions its draws
span, the last summary kept in the others. Score matching is exact for a
Gaussian conditional from a few draws. A window's draws are warm-up, so a
chain that stalls within one finds its step size again where it stands
(`_ChainRetune`); a chain that stopped over a window restarts from a draw of
those that moved; one check window under the final map does the same
without a refit; the kept draws all come after the last refit, under one map,
and stop by the rule every sampler stops by (`_sample_until_target`), judged
on theta.

# Reporting the effects

The effects are reported as natural deviations (`_laplace_deviations!`). Where
the target integrated some, each kept draw takes them from their exact
conditional given the draw -- Gaussian, its mean and precision the hybrid
filter's closing mode and curvature (`_hybrid_member_modes`), mapped to raw
deviations by `_hybrid_dev_map` -- so a sampled fit reports every effect.

Measured on the posterior race's fixtures (`review/SAMPLING-general-joint-map-
2026-10-25.md`): on the collapsed target, with the map, 1.1-17 effective draws
per 1,000 gradients against the layer's NUTS's 0.0-0.9, accurate against long
references throughout.
"""

export ctsem_sample_map


################################################################################
# Theta's triangular map
################################################################################
#
# The whitening `u = Ct^-1 (theta - theta0)` is linear, and theta's posterior is
# not a linear transform of a normal: locations bend with other coordinates,
# scales move with them, marginals are skewed and heavy-tailed. Between the
# latent `w` and `u` sits a monotone triangular map (Marzouk, Moselhy, Parno &
# Spantini 2016; Parno & Marzouk 2018), learnt from the pooled draws of the
# learning windows. Component k, from `u` to `w`:
#
#     v_k = (u_k - f_k(u_<k)) exp(-a_k(u_<k)),   w_k = sinh(delta_k asinh(v_k) - eps_k)
#
# `f_k` quadratic in each earlier coordinate and in the products of pairs of
# them, `a_k` linear in them and capped by tanh, and a sinh-arcsinh marginal (Jones & Pewsey 2009; skew `eps`, tail
# weight `delta`). Both directions are in closed form and each component is
# increasing in its own coordinate, so it is a bijection everywhere and the
# sampler exact whatever was learnt
# (`review/SAMPLING-general-joint-map-2026-10-25.md`).
#
# The windows give a few hundred draws, and the map must not extrapolate past
# them. The earlier coordinates enter `f_k` and `a_k` through `c tanh(u / c)`
# (`_jmap_tri_sat`): the same map near the centre, a location and scale that
# stop bending beyond about `c` sd. Read as `u`, quadratics fitted from 480
# draws moved later locations by up to 22 sd at a draw 14 sd out on gN1, and
# its worst held-out draw lost 970 nats to the normal; saturated, 9, and the
# mean gain grew. `c` is 5: against 3, on nine race fixtures at twelve seeds,
# the map was about as efficient (1.08 times the effective draws a gradient)
# and one run of 108 ended under min ESS 60 to five. And
# `delta` is held at most 1, tails no lighter than the normal's: a small sample
# misses the tails, a fit from it makes them light, and the latent target is
# then heavy-tailed, where NUTS mixes slowly
# (`scripts/marginal-default-2026-10/tri_saturate.R`).
#
# A component conditions on the earlier coordinates its draws depend on most
# (`_jmap_tri_parents`), as many as its held-out likelihood says the draws
# support (`_jmap_tri_fit`): fitted at 0, 1, 2, 4, 8, ... parents on three
# quarters of the draws and scored on the last quarter, which for pooled chains
# is mostly the last chain's, it keeps the best size and stops climbing after
# two sizes that do not improve on it (adaptive transport maps: Baptista,
# Bigoni and Marzouk 2023). So the map is as large as its draws can learn,
# whether they are a stage's 200 particles or thousands of warm-up draws, and
# its coefficients do not grow with the parameters unless the draws carry
# them. On the time series' reference draws the held-out log density of a map
# so chosen was 1.2 nats above four parents without products and 0.3 above the
# best single size with them (`scripts/dais-sample-2026-11/map_ml.R`). Until
# 2026-11-03 every component read four parents; until 2026-11-01 every earlier
# coordinate, and above 30 parameters, where the draws could not fit that,
# there was no map (`review/SAMPLING-sparse-map-2026-11-01.md`).

const _JMAP_TRI_CAP = 3.0
"""Beyond about this many whitened sd an earlier coordinate stops bending a
later one's location and scale."""
const _JMAP_TRI_SAT = 5.0
_jmap_tri_sat(u) = _JMAP_TRI_SAT * tanh(u / _JMAP_TRI_SAT)
_jmap_tri_dsat(u) = 1 - tanh(u / _JMAP_TRI_SAT)^2

"""
Theta's triangular map. Component `k` conditions on the earlier coordinates
`parents[k]`; `f1[k]`, `f2[k]` and `a1[k]` are its location's linear and
quadratic coefficients and its log scale's linear ones, one for each parent,
and `f3[k]` its location's coefficients on the products of pairs of parents,
`(q, r)` for `q < r` in order.
"""
struct _JMapTri
    n::Int
    parents::Vector{Vector{Int}}
    f0::Vector{Float64}
    f1::Vector{Vector{Float64}}
    f2::Vector{Vector{Float64}}
    f3::Vector{Vector{Float64}}
    a0::Vector{Float64}
    a1::Vector{Vector{Float64}}
    eps::Vector{Float64}
    ld::Vector{Float64}
end

# Component `k`'s location and raw log scale at whitened `u`, of which it reads
# only its parents.
@inline function _jmap_tri_terms(tri::_JMapTri, k::Int, u::AbstractVector{T}) where {T}
    f = T(tri.f0[k])
    b = T(tri.a0[k])
    P = tri.parents[k]
    f1 = tri.f1[k]
    f2 = tri.f2[k]
    f3 = tri.f3[k]
    a1 = tri.a1[k]
    @inbounds for q in eachindex(P)
        s = _jmap_tri_sat(u[P[q]])
        f += f1[q] * s + f2[q] * s^2
        b += a1[q] * s
    end
    j = 0
    @inbounds for q in eachindex(P), r in (q + 1):length(P)
        j += 1
        f += f3[j] * _jmap_tri_sat(u[P[q]]) * _jmap_tri_sat(u[P[r]])
    end
    return f, b
end

"""The whitened coordinates at latent `w` and the log Jacobian of that map."""
function _jmap_tri_forward(tri::_JMapTri, w::AbstractVector{T}) where {T<:Real}
    n = tri.n
    A = _JMAP_TRI_CAP
    u = Vector{T}(undef, n)
    lj = zero(T)
    @inbounds for k in 1:n
        f, b = _jmap_tri_terms(tri, k, u)
        a = A * tanh(b / A)
        delta = exp(tri.ld[k])
        r = (asinh(w[k]) + tri.eps[k]) / delta
        u[k] = f + exp(a) * sinh(r)
        lj += a + log(cosh(r)) - log(delta) - log1p(w[k]^2) / 2
    end
    return u, lj
end

"""The latent at whitened `u`."""
function _jmap_tri_inverse(tri::_JMapTri, u::AbstractVector)
    n = tri.n
    A = _JMAP_TRI_CAP
    w = zeros(n)
    @inbounds for k in 1:n
        f, b = _jmap_tri_terms(tri, k, u)
        a = A * tanh(b / A)
        v = (u[k] - f) * exp(-a)
        w[k] = sinh(exp(tri.ld[k]) * asinh(v) - tri.eps[k])
    end
    return w
end

"""The latent gradient from the whitened one, `gu`, the log Jacobian's own
gradient included: the forward pass, then a reverse pass from the last
component, each handing its location's and scale's adjoints to its parents."""
function _jmap_tri_pullback(tri::_JMapTri, w::AbstractVector, gu::AbstractVector)
    n = tri.n
    A = _JMAP_TRI_CAP
    u = zeros(n)
    a = zeros(n)
    da = zeros(n)
    r = zeros(n)
    @inbounds for k in 1:n
        f, b = _jmap_tri_terms(tri, k, u)
        t = tanh(b / A)
        a[k] = A * t
        da[k] = 1 - t^2
        r[k] = (asinh(w[k]) + tri.eps[k]) / exp(tri.ld[k])
        u[k] = f + exp(a[k]) * sinh(r[k])
    end
    ub = collect(Float64, gu)
    gw = zeros(n)
    @inbounds for k in n:-1:1
        e = exp(a[k])
        abar = ub[k] * e * sinh(r[k]) + 1
        rbar = ub[k] * e * cosh(r[k]) + tanh(r[k])
        wk = w[k]
        gw[k] = rbar / (exp(tri.ld[k]) * sqrt(1 + wk^2)) - wk / (1 + wk^2)
        fbar = ub[k]
        bbar = abar * da[k]
        P = tri.parents[k]
        f1 = tri.f1[k]
        f2 = tri.f2[k]
        f3 = tri.f3[k]
        a1 = tri.a1[k]
        # The location's slope in each parent's saturated value: its own terms
        # and, through each product, the other parent's value.
        slope = [f1[q] + 2 * f2[q] * _jmap_tri_sat(u[P[q]]) for q in eachindex(P)]
        jj = 0
        for q in eachindex(P), r in (q + 1):length(P)
            jj += 1
            slope[q] += f3[jj] * _jmap_tri_sat(u[P[r]])
            slope[r] += f3[jj] * _jmap_tri_sat(u[P[q]])
        end
        for q in eachindex(P)
            j = P[q]
            ub[j] += (fbar * slope[q] + bbar * a1[q]) * _jmap_tri_dsat(u[j])
        end
    end
    return gw
end

# A component's negative log likelihood over the whitened draws, each weighted
# by `w` (mean one), in its coefficients `p` = f0, f1, f2, f3, a0, a1 (over its
# parents and their pairs), eps, log delta, with a unit ridge on all but the
# level and the scale. The coefficients enter the location and the raw log
# scale linearly, through the designs `Phi` = [1 s s^2 s_q s_r] and `Psi` =
# [1 s] (a row a draw, `s` the parents' saturated values,
# `_jmap_tri_designs`), so the gradient and Hessian are those designs weighted
# by each draw's first and second derivatives in the location, the raw log
# scale, eps and log delta. `pinned`: `p` stops at eps and log delta is 0.
struct _JMapTriComponent
    x::Vector{Float64}
    Phi::Matrix{Float64}
    Psi::Matrix{Float64}
    pinned::Bool
    w::Vector{Float64}
end
_JMapTriComponent(x, Phi, Psi, pinned) = _JMapTriComponent(x, Phi, Psi, pinned, ones(length(x)))

_jmap_tri_component_value(C::_JMapTriComponent, p::AbstractVector) =
    _jmap_tri_component_nll(C, p) + _jmap_tri_component_ridge(C, p)

function _jmap_tri_component_ridge(C::_JMapTriComponent, p::AbstractVector{T}) where {T}
    nf = size(C.Phi, 2)
    reg = zero(T)
    for j in 2:length(p)
        j == nf + 1 && continue
        reg += p[j]^2
    end
    return reg / 2
end

function _jmap_tri_component_nll(C::_JMapTriComponent, p::AbstractVector{T}) where {T}
    nf = size(C.Phi, 2)
    nb = size(C.Psi, 2)
    f = C.Phi * p[1:nf]
    b = C.Psi * p[(nf + 1):(nf + nb)]
    epsk = p[nf + nb + 1]
    lam = C.pinned ? zero(T) : p[nf + nb + 2]
    delta = exp(lam)
    A = _JMAP_TRI_CAP
    s = zero(T)
    @inbounds for i in eachindex(C.x)
        a = A * tanh(b[i] / A)
        v = (C.x[i] - f[i]) * exp(-a)
        r = delta * asinh(v) - epsk
        s += C.w[i] * (sinh(r)^2 / 2 - (lam + log(cosh(r)) - log1p(v^2) / 2) + a)
    end
    return s
end

function _jmap_tri_component_derivatives(C::_JMapTriComponent, p::Vector{Float64})
    Phi = C.Phi
    Psi = C.Psi
    N = length(C.x)
    nf = size(Phi, 2)
    nb = size(Psi, 2)
    ie = nf + nb + 1
    il = ie + 1
    np = length(p)
    f = Phi * p[1:nf]
    b = Psi * p[(nf + 1):(nf + nb)]
    epsk = p[ie]
    delta = C.pinned ? 1.0 : exp(p[il])
    A = _JMAP_TRI_CAP
    Lf = zeros(N); Lb = zeros(N); Lff = zeros(N); Lfb = zeros(N); Lbb = zeros(N)
    Lfe = zeros(N); Lfl = zeros(N); Lbe = zeros(N); Lbl = zeros(N)
    ge = 0.0; gl = 0.0; Hee = 0.0; Hel = 0.0; Hll = 0.0
    @inbounds for i in 1:N
        t = tanh(b[i] / A)
        a = A * t
        ab = 1 - t^2
        abb = -2 * t * ab / A
        ea = exp(-a)
        v = (C.x[i] - f[i]) * ea
        q = asinh(v)
        sq = sqrt(1 + v^2)
        q1 = 1 / sq
        q2 = -v / sq^3
        r = delta * q - epsk
        sr = sinh(r)
        cr = cosh(r)
        g1 = sr * cr - sr / cr
        g2 = cr^2 + sr^2 - 1 / cr^2
        h1 = v / (1 + v^2)
        h2 = (1 - v^2) / (1 + v^2)^2
        dq = delta * q
        dq1 = delta * q1
        # In the standardised value v, eps and log delta.
        lv = g1 * dq1 + h1
        lvv = g2 * dq1^2 + g1 * delta * q2 + h2
        lve = -g2 * dq1
        lvl = g2 * dq * dq1 + g1 * dq1
        wi = C.w[i]
        ge -= wi * g1
        gl += wi * (g1 * dq - 1)
        Hee += wi * g2
        Hel -= wi * g2 * dq
        Hll += wi * (g2 * dq^2 + g1 * dq)
        # In the location and the capped log scale, v = (x - f) exp(-a), and
        # the log Jacobian's own `+ a`.
        La = -lv * v + 1
        Laa = lvv * v^2 + lv * v
        Lfa = lvv * ea * v + lv * ea
        # The capped log scale through the raw one.
        Lf[i] = -wi * lv * ea
        Lff[i] = wi * lvv * ea^2
        Lb[i] = wi * La * ab
        Lbb[i] = wi * (Laa * ab^2 + La * abb)
        Lfb[i] = wi * Lfa * ab
        Lfe[i] = -wi * lve * ea
        Lfl[i] = -wi * lvl * ea
        Lbe[i] = -wi * lve * v * ab
        Lbl[i] = -wi * lvl * v * ab
    end
    g = zeros(np)
    H = zeros(np, np)
    fi = 1:nf
    bi = (nf + 1):(nf + nb)
    g[fi] .= transpose(Phi) * Lf
    g[bi] .= transpose(Psi) * Lb
    g[ie] = ge
    H[fi, fi] .= transpose(Phi) * (Lff .* Phi)
    H[fi, bi] .= transpose(Phi) * (Lfb .* Psi)
    H[bi, bi] .= transpose(Psi) * (Lbb .* Psi)
    H[fi, ie] .= transpose(Phi) * Lfe
    H[bi, ie] .= transpose(Psi) * Lbe
    H[ie, ie] = Hee
    if !C.pinned
        g[il] = gl
        H[fi, il] .= transpose(Phi) * Lfl
        H[bi, il] .= transpose(Psi) * Lbl
        H[ie, il] = Hel
        H[il, il] = Hll
    end
    for j in 1:np, i in (j + 1):np
        H[i, j] = H[j, i]
    end
    for j in 2:np
        j == nf + 1 && continue
        g[j] += p[j]
        H[j, j] += 1
    end
    return g, H
end

# Damped Newton on the component from `p`, on the Hessian's eigenvalues floored.
function _jmap_tri_newton(C::_JMapTriComponent, p::Vector{Float64}, iters::Int)
    val = _jmap_tri_component_value(C, p)
    for _ in 1:iters
        g, H = _jmap_tri_component_derivatives(C, p)
        (all(isfinite, g) && all(isfinite, H)) || break
        ev = eigen(Symmetric(H))
        lam = max.(ev.values, 1e-6 * maximum(abs, ev.values) + 1e-8)
        step = -(ev.vectors * ((transpose(ev.vectors) * g) ./ lam))
        t = 1.0
        improved = false
        for _ in 1:30
            pn = p .+ t .* step
            vn = _jmap_tri_component_value(C, pn)
            if isfinite(vn) && vn < val - 1e-10 * abs(val)
                improved = abs(val - vn) > 1e-9 * max(1.0, abs(val))
                p = pn
                val = vn
                break
            end
            t /= 2
        end
        improved || break
    end
    return p
end

# The squared correlation of `x` and `y` over the draws.
function _jmap_cor2(x::AbstractVector, y::AbstractVector)
    mx = sum(x) / length(x)
    my = sum(y) / length(y)
    sxy = 0.0; sxx = 0.0; syy = 0.0
    @inbounds for i in eachindex(x, y)
        dx = x[i] - mx
        dy = y[i] - my
        sxy += dx * dy
        sxx += dx^2
        syy += dy^2
    end
    return sxx > 0 && syy > 0 ? sxy^2 / (sxx * syy) : 0.0
end

# How much component `k`'s draws depend on each earlier coordinate -- its
# location linearly or bending, or its spread: the squared correlations of
# `u_k` with `s_j` and with `s_j^2`, and of `(u_k - mean)^2` with `s_j`, summed.
# `X` holds the draws (a row a coordinate); `Sd` and `S2d` the saturated values
# and their squares a column a coordinate, so each is read contiguously.
function _jmap_tri_scores(X::AbstractMatrix, Sd::AbstractMatrix, S2d::AbstractMatrix,
    k::Int)
    x = X[k, :]
    xc = x .- sum(x) / length(x)
    x2 = xc .^ 2
    score = zeros(k - 1)
    for j in 1:(k - 1)
        s = view(Sd, :, j)
        score[j] = _jmap_cor2(xc, s) + _jmap_cor2(xc, view(S2d, :, j)) + _jmap_cor2(x2, s)
    end
    return score
end

# The earlier coordinates component `k` conditions on: every one where there are
# at most `m`, otherwise the `m` its draws depend on most (`_jmap_tri_scores`).
function _jmap_tri_parents(score::Vector{Float64}, k::Int, m::Int)
    k - 1 <= m && return collect(1:(k - 1))
    m <= 0 && return Int[]
    return sort(partialsortperm(score, 1:m; rev=true))
end

"""
A component's designs over parents `P` at draws `S` (saturated values, a column
a draw): the location's [1, s, s^2, s_q s_r for q < r] and the raw log scale's
[1, s], a row a draw.
"""
function _jmap_tri_designs(S::AbstractMatrix, P::Vector{Int})
    N = size(S, 2)
    m = length(P)
    Phi = Matrix{Float64}(undef, N, 1 + 2m + m * (m - 1) ÷ 2)
    Psi = Matrix{Float64}(undef, N, 1 + m)
    @inbounds for i in 1:N
        Phi[i, 1] = 1.0
        Psi[i, 1] = 1.0
        for q in 1:m
            s = S[P[q], i]
            Phi[i, 1 + q] = s
            Phi[i, 1 + m + q] = s^2
            Psi[i, 1 + q] = s
        end
        j = 1 + 2m
        for q in 1:m, r in (q + 1):m
            j += 1
            Phi[i, j] = S[P[q], i] * S[P[r], i]
        end
    end
    return Phi, Psi
end

# A component fitted by damped Newton from `start` (the identity where none is
# given); one whose tails come out lighter than the normal's is fitted again
# with its tail weight held at the normal's. Its coefficients, log delta last,
# and whether it was held.
function _jmap_tri_component_fit(x, Phi, Psi, w, iters; start=nothing)
    np = size(Phi, 2) + size(Psi, 2) + 2
    p = _jmap_tri_newton(_JMapTriComponent(x, Phi, Psi, false, w),
        start === nothing ? zeros(np) : start, iters)
    p[end] > 0 || return p, false
    return vcat(_jmap_tri_newton(_JMapTriComponent(x, Phi, Psi, true, w), p[1:(end - 1)],
        iters), 0.0), true
end

# A component's coefficients `p` on parents `P`, laid out for parents `Q` that
# hold them (both in order): each parent's and each pair's where it now sits,
# zero for the new; `nothing` where `Q` does not hold `P`. A larger size starts
# from the smaller one's fit, the same model with the new terms at zero.
function _jmap_tri_embed(p::Vector{Float64}, P::Vector{Int}, Q::Vector{Int})
    pos = [findfirst(==(j), Q) for j in P]
    any(isnothing, pos) && return nothing
    m = length(P)
    M = length(Q)
    nf = 1 + 2m + m * (m - 1) ÷ 2
    nF = 1 + 2M + M * (M - 1) ÷ 2
    out = zeros(nF + 1 + M + 2)
    out[1] = p[1]
    for q in 1:m
        out[1 + pos[q]] = p[1 + q]
        out[1 + M + pos[q]] = p[1 + m + q]
        out[nF + 1 + pos[q]] = p[nf + 1 + q]
    end
    j = 0
    for q in 1:m, r in (q + 1):m
        j += 1
        lo, hi = pos[q], pos[r]
        out[1 + 2M + (lo - 1) * M - (lo - 1) * lo ÷ 2 + (hi - lo)] = p[1 + 2m + j]
    end
    out[nF + 1] = p[nf + 1]
    out[end - 1] = p[end - 1]
    out[end] = p[end]
    return out
end

"""
    _jmap_tri_fit(U; weights, iters, parallel, progress)

The triangular map fitted to whitened draws `U` (a column each, weighted by
`weights` where given) by maximum likelihood, a component at a time. Each
component's size is its held-out likelihood's: fitted at 0, 1, 2, 4, 8, ...
parents (`_jmap_tri_parents`, at most every earlier coordinate) on the first
three quarters of the draws and scored on the last quarter, the best size kept
and the ladder left after two sizes that do not improve on it, or where a size
has more coefficients than the draws it is fitted on carry (their effective
size, where they are weighted); then refitted on every draw. Each size's
Newton starts from the last size's fit, and the refit from the size kept: the
same fits, in fewer iterations. `nothing` where too few draws fit even the
marginal.

The components depend on the draws alone, not on each other's fits, so with
`parallel` (a caller that has sized the pool, `_laplace_ensure_pool!`) they
run on it, OpenBLAS held to one thread meanwhile so that the workers' products
do not queue on its pool; a component's fit is the same wherever it runs.
`progress(done, n)` is called from the calling task as components finish
(`_ctsem_call_hook`): on the reanalysis theta had 2,868 coordinates (2,803 of
them sampled TI predictor values), fitted one after another in silence after
every window.
"""
function _jmap_tri_fit(U::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing,
    iters::Int=60, parallel::Bool=false, progress=nothing)
    n, N = size(U)
    nv = cld(N, 4)
    nt = N - nv
    nt > 4 || return nothing
    X = Matrix{Float64}(U)
    S = _jmap_tri_sat.(X)
    w = weights === nothing ? ones(N) : collect(Float64, weights)
    wt = w[1:nt] ./ (sum(w[1:nt]) / nt)
    wv = w[(nt + 1):N] ./ (sum(w[(nt + 1):N]) / nv)
    wall = w ./ (sum(w) / N)
    # What the training draws carry: their effective size, which weighted
    # draws (a tempering stage's cloud) hold fewer of than their count.
    neff = sum(wt)^2 / sum(abs2, wt)
    tr = 1:nt
    va = (nt + 1):N
    Xt = X[:, tr]
    St = S[:, tr]
    Sv = S[:, va]
    # The training draws' saturated values and their squares, a column a
    # coordinate, for the parents' scores.
    Sd = permutedims(St)
    S2d = Sd .^ 2
    fits = Vector{Tuple{Vector{Int},Vector{Float64},Int}}(undef, n)
    done = Threads.Atomic{Int}(0)
    caller = current_task()
    body = function (k)
        fits[k] = _jmap_tri_component_ladder(X, S, Xt, St, Sv, Sd, S2d, wt, wv, wall,
            neff, va, k, iters)
        Threads.atomic_add!(done, 1)
        current_task() === caller && _ctsem_call_hook(progress, done[], n)
        return true
    end
    if parallel && n > 1
        nb = LinearAlgebra.BLAS.get_num_threads()
        LinearAlgebra.BLAS.set_num_threads(1)
        try
            _laplace_parallel(body, 1:n)
        finally
            LinearAlgebra.BLAS.set_num_threads(nb)
        end
    else
        foreach(body, 1:n)
    end
    _ctsem_call_hook(progress, n, n)
    par = Vector{Vector{Int}}(undef, n)
    f1 = Vector{Vector{Float64}}(undef, n)
    f2 = Vector{Vector{Float64}}(undef, n)
    f3 = Vector{Vector{Float64}}(undef, n)
    a1 = Vector{Vector{Float64}}(undef, n)
    f0 = zeros(n); a0 = zeros(n); epsv = zeros(n); ld = zeros(n)
    for k in 1:n
        P, p, nf = fits[k]
        mm = length(P)
        par[k] = P
        f0[k] = p[1]
        f1[k] = p[2:(mm + 1)]
        f2[k] = p[(mm + 2):(2mm + 1)]
        f3[k] = p[(2mm + 2):nf]
        a0[k] = p[nf + 1]
        a1[k] = p[(nf + 2):(nf + mm + 1)]
        epsv[k] = p[nf + mm + 2]
        ld[k] = p[nf + mm + 3]
    end
    return _JMapTri(n, par, f0, f1, f2, f3, a0, a1, epsv, ld)
end

# Component `k` of `_jmap_tri_fit`: its parents' ladder on the training draws,
# scored on the held-out ones, then the size kept refitted on every draw. Its
# parents, coefficients and the location's design width.
function _jmap_tri_component_ladder(X, S, Xt, St, Sv, Sd, S2d, wt, wv, wall, neff, va,
    k::Int, iters::Int)
    score = _jmap_tri_scores(Xt, Sd, S2d, k)
    xk = Xt[k, :]
    xv = X[k, va]
    best = (nll=Inf, P=Int[], p=Float64[])
    prev = nothing
    worse = 0
    m = 0
    while true
        P = _jmap_tri_parents(score, k, m)
        Phi, Psi = _jmap_tri_designs(St, P)
        size(Phi, 2) + size(Psi, 2) + 2 <= neff || break
        p0 = prev === nothing ? nothing : _jmap_tri_embed(prev.p, prev.P, P)
        p, held = _jmap_tri_component_fit(xk, Phi, Psi, wt, iters; start=p0)
        prev = (P=P, p=p)
        Phv, Psv = _jmap_tri_designs(Sv, P)
        nll = _jmap_tri_component_nll(_JMapTriComponent(xv, Phv, Psv, held, wv), p)
        if nll < best.nll
            best = (nll=nll, P=P, p=p)
            worse = 0
        else
            worse += 1
        end
        (worse >= 2 || length(P) >= k - 1) && break
        m = m == 0 ? 1 : 2m
    end
    P = best.P
    Phi, Psi = _jmap_tri_designs(S, P)
    p, _ = _jmap_tri_component_fit(X[k, :], Phi, Psi, wall, iters;
        start=isempty(best.p) ? nothing : best.p)
    return (P, p, size(Phi, 2))
end

"""
    CTSEMJointMap

The map: theta's centre and whitening factor, and per member of each unit its
Gaussian summary in its raw effects (`Lam`, `h`) and coupling `G` to the free
parameters outside its effects that its subject reads: `theta[R]`, those every
subject reads, then `theta[own[s]]`, its subject's own sampled TI predictor
values (`_jmap_layout`). `P` are the raw positions varying at some level,
`row` a raw position's row of `P`, `thref` the parameters (extended) the
couplings are centred at.
"""
struct CTSEMJointMap{S}
    sampler::S
    npar::Int
    P::Vector{Int}
    row::Vector{Int}
    R::Vector{Int}
    own::Vector{Vector{Int}}
    members::Vector{Vector{Int}}
    Lam::Vector{Matrix{Float64}}
    h::Vector{Vector{Float64}}
    G::Vector{Matrix{Float64}}
    thref::Vector{Float64}
    theta0::Vector{Float64}
    Ct::Matrix{Float64}
    tri::Union{Nothing,_JMapTri}
end

_jmap_with(m::CTSEMJointMap; Lam=m.Lam, h=m.h, G=m.G, theta0=m.theta0, Ct=m.Ct,
    tri=m.tri) =
    CTSEMJointMap(m.sampler, m.npar, m.P, m.row, m.R, m.own, m.members, Lam, h, G,
        m.thref, theta0, Ct, tri)

"""Member `s`'s coupled positions, `R` then its own: the columns of `G[s]`."""
_jmap_coupled(m::CTSEMJointMap, s::Int) = vcat(m.R, m.own[s])

"""
    _JMapTheta

Theta's part of the map alone, `theta = theta0 + Ct T(w)`, for a sampler that
moves theta and nothing else through it: tempered SMC's moves and final
chains (`_smc_geometry`, sample_smc.jl), refitted from the cloud at every
stage. The functions of theta's part below take either.
"""
struct _JMapTheta
    npar::Int
    theta0::Vector{Float64}
    Ct::Matrix{Float64}
    tri::Union{Nothing,_JMapTri}
end

_jmap_with(m::_JMapTheta; tri=m.tri) = _JMapTheta(m.npar, m.theta0, m.Ct, tri)

const _JMapThetaLike = Union{CTSEMJointMap,_JMapTheta}

"""
    _jmap_learn_tri(m, pooled; weights)

The map with theta's triangular map fitted from `pooled` (draws of theta, a
column each, weighted where `weights` are given: a tempering stage's cloud),
whitened by `m`'s own centre and factor, each component as large as the draws
support (`_jmap_tri_fit`); `m` as it is where they support none. At most 4000
draws are read, spread through the pool in its order, so the held-out quarter
stays the pool's last. `parallel` and `progress` are `_jmap_tri_fit`'s.
"""
function _jmap_learn_tri(m::_JMapThetaLike, pooled::AbstractMatrix;
    weights::Union{Nothing,AbstractVector}=nothing, parallel::Bool=false,
    progress=nothing)
    N = size(pooled, 2)
    pick = N > 4000 ? round.(Int, range(1, N; length=4000)) : collect(1:N)
    U = LowerTriangular(m.Ct) \ (pooled[:, pick] .- m.theta0)
    tri = _jmap_tri_fit(U; weights=weights === nothing ? nothing : weights[pick],
        parallel=parallel, progress=progress)
    return tri === nothing ? m : _jmap_with(m; tri=tri)
end

# Theta's whitening and its inverse go through the engine's own products
# (small_linalg.jl), never BLAS: tempered SMC moves its particles through them
# on tasks, and its draws are the same whether the particles run one after
# another or on copies (test_smc.jl) only where every product's sums are
# taken in one order.

"""Theta at latent `w`, and the log Jacobian of theta's triangular map (zero
without one; the whitening's, `log det Ct`, is a constant and left out)."""
function _jmap_theta(m::_JMapThetaLike, w::AbstractVector{<:Real})
    u, lj = m.tri === nothing ? (w, 0.0) : _jmap_tri_forward(m.tri, w)
    theta = Vector{promote_type(eltype(u), Float64)}(m.theta0)
    _ctsem_mulvec!(theta, m.Ct, u, true, true)
    return theta, lj
end

"""Theta at each column of latent `W`."""
function _jmap_thetas(m::_JMapThetaLike, W::AbstractMatrix)
    out = zeros(m.npar, size(W, 2))
    for t in axes(W, 2)
        out[:, t] .= _jmap_theta(m, view(W, :, t))[1]
    end
    return out
end

"""Theta's latent coordinates."""
function _jmap_theta_inverse(m::_JMapThetaLike, theta::AbstractVector)
    L = m.Ct
    u = collect(Float64, theta) .- m.theta0
    @inbounds for i in eachindex(u)
        acc = u[i]
        for k in 1:(i - 1)
            acc -= L[i, k] * u[k]
        end
        u[i] = acc / L[i, i]
    end
    m.tri === nothing && return u
    return _jmap_tri_inverse(m.tri, u)
end

"""The latent gradient in theta's coordinates from the joint's gradient in
theta, `gtheta`, through the whitening and the triangular map, its log
Jacobian's gradient included."""
function _jmap_theta_pullback(m::_JMapThetaLike, w::AbstractVector,
    gtheta::AbstractVector)
    gu = _ctsem_mulTvec!(zeros(m.npar), m.Ct, gtheta)
    m.tri === nothing && return gu
    return _jmap_tri_pullback(m.tri, w, gu)
end


# The positions `P` that vary, each one's row, and with coupling (`cross`) what
# each member's summary is coupled to: `R`, the free parameters outside `P`
# that every subject reads, and per member (flat index) its subject's own
# sampled TI predictor values (`TIMissingRecipe`). A subject's score is exactly
# zero at another subject's values, so a coupling to them is a zero column: on
# the reanalysis (8,446 members, 2,803 sampled values) every member carried
# 2,860 columns where it read 62 and its own, about 1.5 GB a matrix of `G` and
# of the refit's sums, read at every density evaluation.
function _jmap_layout(sampler::CTSEMSampler, cross::Bool)
    laplace = sampler.laplace
    spec = laplace.spec
    nv = length(_laplace_extend(laplace, zeros(sampler.npar)))
    P = sort(unique(reduce(vcat, [Vector{Int}(level.re_index) for level in spec.levels];
        init=Int[])))
    row = zeros(Int, nv)
    for (j, p) in enumerate(P)
        row[p] = j
    end
    nflat = length(sampler.flat_unit)
    own = [Int[] for _ in 1:nflat]
    R = Int[]
    if cross
        sampled = falses(nv)
        for j in laplace.objective.ti_missing_parameter
            sampled[j] = true
        end
        R = [j for j in 1:sampler.npar if row[j] == 0 && !sampled[j]]
        for s in 1:nflat
            U = sampler.flat_unit[s]
            i = laplace.units.members[U][sampler.flat_member[s]]
            tip = laplace.objective.subject_objectives[i].tipreds
            tip isa TIMissingRecipe || continue
            own[s] = [j for j in tip.parameter_index if j <= sampler.npar && row[j] == 0]
        end
    end
    members = [Int[] for _ in 1:sampler.nunits]
    for s in eachindex(sampler.flat_unit)
        sampler.udims[sampler.flat_unit[s]] > 0 && push!(members[sampler.flat_unit[s]], s)
    end
    return P, row, R, own, members
end

# Member s's J_s, |P| by its unit's dimension, from the levels' roots `Ls`.
function _jmap_jacobian!(J::AbstractMatrix{T}, m::CTSEMJointMap, Ls, s::Int) where {T}
    sampler = m.sampler
    laplace = sampler.laplace
    spec = laplace.spec
    U = sampler.flat_unit[s]
    offsets = laplace.units.offsets[U][sampler.flat_member[s]]
    fill!(J, zero(T))
    @inbounds for l in eachindex(spec.levels)
        level = spec.levels[l]
        k = nrandomeffects(level)
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        base = offsets[l]
        L = Ls[l]
        for p in 1:k
            i = m.row[level.re_index[p]]
            for q in 1:r
                J[i, base + q] += L[p, q]
            end
        end
    end
    return J
end

# Member s's linear term `h_s - Lam_s mu - G_s shift_s` at row `a`, where its
# shift is `shift` (theta[R] less its reference) then its own values less
# theirs, read from the extended theta `th`.
@inline function _jmap_lin(m::CTSEMJointMap, s::Int, a::Int, mu, shift, th)
    Lam = m.Lam[s]
    acc = m.h[s][a]
    nP = length(m.P)
    @inbounds for b in 1:nP
        acc -= Lam[a, b] * mu[b]
    end
    G = m.G[s]
    nR = length(shift)
    @inbounds for b in 1:nR
        acc -= G[a, b] * shift[b]
    end
    own = m.own[s]
    @inbounds for j in eachindex(own)
        acc -= G[a, nR + j] * (th[own[j]] - m.thref[own[j]])
    end
    return acc
end

# Unit U's conditional given theta for the members' summaries: Q_U factored
# (`_JMapUnitFactor`) and the mean m_U. `Ls`, `mu = theta[P]`, `shift =
# theta[R] - thref[R]` and the extended theta `th` at theta. Lock-free kernels
# throughout: it runs on every leapfrog step of chains that run at once. A unit
# of nested blocks -- persons in a study -- is built and factored on its blocks
# (`_jmap_unit_blocked`); a unit of one block in place (`_ctsem_cholesky`,
# upper).
function _jmap_unit(m::CTSEMJointMap, Ls, mu::AbstractVector{T}, shift::AbstractVector,
    th::AbstractVector, U::Int) where {T}
    blocks = m.sampler.laplace.units.blocks[U]
    _JMAP_BLOCKS[] && length(blocks) > 1 &&
        return _jmap_unit_blocked(m, Ls, mu, shift, th, U, blocks)
    d = m.sampler.udims[U]
    nP = length(m.P)
    Q = Matrix{T}(I, d, d)
    c = zeros(T, d)
    J = zeros(T, nP, d)
    LJ = zeros(T, nP, d)
    lin = zeros(T, nP)
    for s in m.members[U]
        _jmap_jacobian!(J, m, Ls, s)
        Lam = m.Lam[s]
        _ctsem_mul!(LJ, Lam, J)
        _ctsem_mulTN!(Q, J, LJ, true, true)
        @inbounds for a in 1:nP
            lin[a] = _jmap_lin(m, s, a, mu, shift, th)
        end
        @inbounds for q in 1:d, a in 1:nP
            c[q] += J[a, q] * lin[a]
        end
    end
    F = _ctsem_cholesky(Q, d)
    F.ok || throw(DomainError(U, "a unit's conditional precision did not factor"))
    return _JMapUnitFactor(F), F \ c
end

"""Whether the map builds a nested unit's `Q_U` from its blocks, or densely."""
const _JMAP_BLOCKS = Ref(true)

"""
    ctsem_set_map_blocks!(on)

Build a nested unit's conditional precision in the learnt map from the unit's
blocks (`true`, the default), or densely (`false`, the oracle): the same
precision and mean, its root differing by a rotation, at the unit's width
squared in storage and cubed a leapfrog step. Returns the previous setting.
"""
function ctsem_set_map_blocks!(on::Bool)
    previous = _JMAP_BLOCKS[]
    _JMAP_BLOCKS[] = on
    return previous
end
export ctsem_set_map_blocks!

"""
A unit's conditional precision `Q_U` factored as the map's root: densely, `Q_U
= U'U` and the root `U` (`dense`), or for a unit of nested blocks from its
block elimination, `Q_U = K K'` and the root `K'` (`factors`, `coupling`,
`blocks`; `_laplace_block_draw`). Either root makes the change of variables
exact; they differ by a rotation of the unit's latent coordinates.
"""
struct _JMapUnitFactor{T}
    dense::CTSEMCholesky{T,Matrix{T}}
    factors::Vector{CTSEMCholesky{T,Matrix{T}}}
    coupling::Vector{Vector{Matrix{T}}}
    blocks::Vector{CTSEMLaplaceBlock}
end
_JMapUnitFactor(F::CTSEMCholesky{T,Matrix{T}}) where {T} = _JMapUnitFactor{T}(F,
    CTSEMCholesky{T,Matrix{T}}[], Vector{Matrix{T}}[], CTSEMLaplaceBlock[])
_JMapUnitFactor(factors::Vector{CTSEMCholesky{T,Matrix{T}}},
    coupling::Vector{Vector{Matrix{T}}}, blocks::Vector{CTSEMLaplaceBlock}) where {T} =
    _JMapUnitFactor{T}(CTSEMCholesky(zeros(T, 0, 0), 0, true), factors, coupling, blocks)
_jmap_blocked(F::_JMapUnitFactor) = !isempty(F.factors)

# The root's algebra: `R^-1 b`, `R^-T b`, `R b` and `log det R`, for either root.
_jmap_usolve(F::_JMapUnitFactor, b::AbstractVector) = _jmap_blocked(F) ?
    _laplace_block_draw(F.factors, F.coupling, F.blocks, b) : _jmap_usolve(F.dense, b)
_jmap_utsolve(F::_JMapUnitFactor, b::AbstractVector) = _jmap_blocked(F) ?
    _laplace_block_lsolve(F.factors, F.coupling, F.blocks, b) : _jmap_utsolve(F.dense, b)
_jmap_umul(F::_JMapUnitFactor, b::AbstractVector) = _jmap_blocked(F) ?
    _laplace_block_rmul(F.factors, F.coupling, F.blocks, b) : _jmap_umul(F.dense, b)
function _jmap_halflogdet(F::_JMapUnitFactor{T}) where {T}
    _jmap_blocked(F) && return _laplace_block_halflogdet(F.factors, F.blocks)
    acc = zero(T)
    @inbounds for i in 1:F.dense.d
        acc += log(F.dense.U[i, i])
    end
    return acc
end

"""A unit's blocks by their offset in its coordinates."""
_jmap_blockat(blocks::Vector{CTSEMLaplaceBlock}) =
    Dict{Int,Int}(b.offset => t for (t, b) in enumerate(blocks))

"""
Member `s`'s chain of blocks in unit `U`, innermost first: `(level, block,
column)` for each level it has coordinates at, `column` where that block starts
in the chain's own columns; and the chain's width. Each block after the first
is an ancestor of every block before it.
"""
function _jmap_member_chain(m::CTSEMJointMap, U::Int, s::Int,
    blocks::Vector{CTSEMLaplaceBlock}, blockat::Dict{Int,Int})
    laplace = m.sampler.laplace
    offsets = laplace.units.offsets[U][m.sampler.flat_member[s]]
    chain = Tuple{Int,Int,Int}[]
    col = 0
    for (l, level) in enumerate(laplace.spec.levels)
        (nrandomeffects(level) == 0 || nlatent(level) == 0) && continue
        b = blockat[offsets[l]]
        push!(chain, (l, b, col))
        col += blocks[b].size
    end
    return chain, col
end

"""Member `s`'s `J_s` on its chain's columns, as `_jmap_jacobian!` lays it on
the unit's."""
function _jmap_chain_jacobian!(J::AbstractMatrix{T}, m::CTSEMJointMap, Ls,
    chain::Vector{Tuple{Int,Int,Int}}) where {T}
    spec = m.sampler.laplace.spec
    fill!(J, zero(T))
    @inbounds for (l, _, col) in chain
        level = spec.levels[l]
        L = Ls[l]
        for p in 1:nrandomeffects(level)
            i = m.row[level.re_index[p]]
            for q in 1:nlatent(level)
                J[i, col + q] += L[p, q]
            end
        end
    end
    return J
end

"""Where block `e` sits in block `c`'s chain, `c` then its ancestors: 1 for `c`
itself, `nothing` where `e` is neither."""
function _jmap_chain_slot(blocks::Vector{CTSEMLaplaceBlock}, c::Int, e::Int)
    c == e && return 1
    t = findfirst(==(e), blocks[c].ancestors)
    return t === nothing ? nothing : t + 1
end

# Unit U's conditional for a unit of nested blocks: `Q_U` in the unit's block
# storage, each member adding to its own chain's diagonal blocks and couplings,
# factored by the block elimination, its mean by a block solve. Storage and
# work are linear in the members, where the dense build was the width squared
# and its factor the width cubed, at every leapfrog step.
function _jmap_unit_blocked(m::CTSEMJointMap, Ls, mu::AbstractVector{T},
    shift::AbstractVector, th::AbstractVector, U::Int,
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    nP = length(m.P)
    Q = CTSEMBlockMatrix(T, blocks)
    for D in Q.diag
        for i in axes(D, 1)
            D[i, i] = one(T)
        end
    end
    c = zeros(T, m.sampler.udims[U])
    blockat = _jmap_blockat(blocks)
    lin = zeros(T, nP)
    for s in m.members[U]
        chain, nc = _jmap_member_chain(m, U, s, blocks, blockat)
        Jc = zeros(T, nP, nc)
        _jmap_chain_jacobian!(Jc, m, Ls, chain)
        LJ = zeros(T, nP, nc)
        _ctsem_mul!(LJ, m.Lam[s], Jc)
        @inbounds for a in 1:nP
            lin[a] = _jmap_lin(m, s, a, mu, shift, th)
        end
        for k1 in eachindex(chain)
            _, b1, c1 = chain[k1]
            r1 = blocks[b1].size
            o1 = blocks[b1].offset
            J1 = view(Jc, :, (c1 + 1):(c1 + r1))
            @inbounds for q in 1:r1, a in 1:nP
                c[o1 + q] += J1[a, q] * lin[a]
            end
            _ctsem_mulTN!(Q.diag[b1], J1, view(LJ, :, (c1 + 1):(c1 + r1)), true, true)
            for k2 in (k1 + 1):length(chain)
                _, b2, c2 = chain[k2]
                t = findfirst(==(b2), blocks[b1].ancestors)::Int
                _ctsem_mulTN!(Q.coupling[b1][t], J1,
                    view(LJ, :, (c2 + 1):(c2 + blocks[b2].size)), true, true)
            end
        end
    end
    ok, _, factors, coupling = _laplace_block_factor(Q, blocks)
    ok || throw(DomainError(U, "a unit's conditional precision did not factor"))
    return _JMapUnitFactor(factors, coupling, blocks),
        _laplace_block_solve(factors, coupling, blocks, c)
end

"""
Murray's `P` (`_jmap_unit_adjoint!`) on blocks `a` by `a2`, in the block
elimination's order -- block, then within it -- which is the order of the
lower factor `K`: off the diagonal the later coordinate's `z` times the
earlier's `w`, halved and negated; on it `-(z_i w_i + 1) / 2`.
"""
function _jmap_murray_block(z::Vector{Float64}, w::Vector{Float64},
    blocks::Vector{CTSEMLaplaceBlock}, a::Int, a2::Int)
    oa, ka = blocks[a].offset, blocks[a].size
    ob, kb = blocks[a2].offset, blocks[a2].size
    P = Matrix{Float64}(undef, ka, kb)
    @inbounds for j in 1:kb, i in 1:ka
        ii, jj = oa + i, ob + j
        P[i, j] = if a == a2 && i == j
            -(z[ii] * w[ii] + 1) / 2
        elseif a > a2 || (a == a2 && i > j)
            -z[ii] * w[jj] / 2
        else
            -z[jj] * w[ii] / 2
        end
    end
    return P
end

# Unit U's reverse pass on its blocks: `_jmap_unit_adjoint!`'s algebra with the
# root `K'` of the block elimination, which is the Cholesky factor of `Q` in the
# elimination's order, so Murray's adjoint holds there. `X = K^-T` has a block
# at `(c, a)` only for `a` in `c`'s chain, and `Q`, so its gradient `S`, is read
# only on the pattern -- each block with itself and its ancestors -- whose
# entries need `X` and `P` on chains alone. Each member then takes `S` on its
# own chain, as the dense pass takes it on the unit.
function _jmap_unit_adjoint_blocked!(work, m::CTSEMJointMap, Ls, mu,
    shift, th, U::Int, F::_JMapUnitFactor{Float64}, u::Vector{Float64},
    g::Vector{Float64}, z::AbstractVector)
    sampler = m.sampler
    spec = sampler.laplace.spec
    blocks, factors, coupling = F.blocks, F.factors, F.coupling
    nb = length(blocks)
    nP = length(m.P)
    nR = length(m.R)
    coupled = !isempty(work.Rbar)
    zv = collect(Float64, z)
    y = _laplace_block_draw(factors, coupling, blocks, zv)
    w = _laplace_block_lsolve(factors, coupling, blocks, g)
    v = _laplace_block_draw(factors, coupling, blocks, w)
    mean = u .- y
    part(x, c) = view(x, (blocks[c].offset + 1):(blocks[c].offset + blocks[c].size))
    # X's blocks along each chain, outermost block first: `X[c][1] = U_c^-1`,
    # and at the chain's k-th ancestor `-D_c^-1` of the couplings times the
    # ancestors' own blocks there.
    X = Vector{Vector{Matrix{Float64}}}(undef, nb)
    for c in nb:-1:1
        anc = blocks[c].ancestors
        Xc = Vector{Matrix{Float64}}(undef, 1 + length(anc))
        Xc[1] = _ctsem_cholesky_uinv(factors[c])
        for (k, a) in enumerate(anc)
            acc = zeros(blocks[c].size, blocks[a].size)
            for (t, at) in enumerate(anc)
                slot = _jmap_chain_slot(blocks, at, a)
                slot === nothing && continue
                acc .+= coupling[c][t] * X[at][slot]
            end
            Xc[k + 1] = -(factors[c] \ acc)
        end
        X[c] = Xc
    end
    # S on the pattern: `S[c][k]` at block c and its chain's k-th block.
    S = Vector{Vector{Matrix{Float64}}}(undef, nb)
    for c in 1:nb
        chc = vcat(c, blocks[c].ancestors)
        Sc = Vector{Matrix{Float64}}(undef, length(chc))
        for (k, e) in enumerate(chc)
            che = vcat(e, blocks[e].ancestors)
            acc = zeros(blocks[c].size, blocks[e].size)
            for (i1, a) in enumerate(chc), (i2, a2) in enumerate(che)
                acc .+= X[c][i1] * _jmap_murray_block(zv, w, blocks, a, a2) *
                    transpose(X[e][i2])
            end
            vc, ve, mc, me = part(v, c), part(v, e), part(mean, c), part(mean, e)
            @inbounds for j in axes(acc, 2), i in axes(acc, 1)
                acc[i, j] -= (vc[i] * me[j] + mc[i] * ve[j]) / 2
            end
            Sc[k] = acc
        end
        S[c] = Sc
    end
    lin = work.lin
    linb = work.linb
    blockat = _jmap_blockat(blocks)
    for s in m.members[U]
        chain, nc = _jmap_member_chain(m, U, s, blocks, blockat)
        Jc = zeros(nP, nc)
        _jmap_chain_jacobian!(Jc, m, Ls, chain)
        Lam = m.Lam[s]
        LJ = zeros(nP, nc)
        _ctsem_mul!(LJ, Lam, Jc)
        Sch = zeros(nc, nc)
        vch = zeros(nc)
        for (k1, (_, b1, c1)) in enumerate(chain)
            r1 = blocks[b1].size
            vch[(c1 + 1):(c1 + r1)] .= part(v, b1)
            for (k2, (_, b2, c2)) in enumerate(chain)
                cols = (c2 + 1):(c2 + blocks[b2].size)
                if k2 >= k1
                    Sch[(c1 + 1):(c1 + r1), cols] .= S[b1][_jmap_chain_slot(blocks, b1, b2)::Int]
                else
                    Sch[(c1 + 1):(c1 + r1), cols] .=
                        transpose(S[b2][_jmap_chain_slot(blocks, b2, b1)::Int])
                end
            end
        end
        @inbounds for a in 1:nP
            lin[a] = _jmap_lin(m, s, a, mu, shift, th)
            t = 0.0
            for q in 1:nc
                t += Jc[a, q] * vch[q]
            end
            linb[a] = t
        end
        Jb = zeros(nP, nc)
        _ctsem_mul!(Jb, LJ, Sch, 2.0, false)
        @inbounds for q in 1:nc, a in 1:nP
            Jb[a, q] += lin[a] * vch[q]
        end
        @inbounds for a in 1:nP
            acc = 0.0
            for b in 1:nP
                acc += Lam[a, b] * linb[b]
            end
            work.mubar[a] -= acc
        end
        if coupled
            G = m.G[s]
            own = m.own[s]
            @inbounds for b in 1:(nR + length(own))
                acc = 0.0
                for a in 1:nP
                    acc += G[a, b] * linb[a]
                end
                work.Rbar[b <= nR ? m.R[b] : own[b - nR]] -= acc
            end
        end
        # The roots' cotangents, where `_jmap_chain_jacobian!` placed each entry.
        @inbounds for (l, _, col) in chain
            level = spec.levels[l]
            Lb = work.Lbar[l]
            for p in 1:nrandomeffects(level)
                i = m.row[level.re_index[p]]
                for q in 1:nlatent(level)
                    Lb[p, q] += Jb[i, col + q]
                end
            end
        end
    end
    return work
end

# `U \ b` and `U' \ b` for a factor's upper triangle.
function _jmap_usolve(F::CTSEMCholesky, b::AbstractVector)
    U = F.U
    d = F.d
    x = Vector{promote_type(eltype(U), eltype(b))}(undef, d)
    @inbounds for i in d:-1:1
        t = b[i]
        for k in (i + 1):d
            t -= U[i, k] * x[k]
        end
        x[i] = t / U[i, i]
    end
    return x
end

function _jmap_utsolve(F::CTSEMCholesky, b::AbstractVector)
    U = F.U
    d = F.d
    x = Vector{promote_type(eltype(U), eltype(b))}(undef, d)
    @inbounds for i in 1:d
        t = b[i]
        for k in 1:(i - 1)
            t -= U[k, i] * x[k]
        end
        x[i] = t / U[i, i]
    end
    return x
end

function _jmap_umul(F::CTSEMCholesky, b::AbstractVector)
    U = F.U
    d = F.d
    y = Vector{promote_type(eltype(U), eltype(b))}(undef, d)
    @inbounds for i in 1:d
        t = zero(eltype(y))
        for k in i:d
            t += U[i, k] * b[k]
        end
        y[i] = t
    end
    return y
end

"""
    _jmap_effects(m, theta, xt)

Every unit's effects at theta for latent `xt`, the map's log Jacobian in the
units, and each unit's factor. Units are independent, so they run on the pool
(`_laplace_partition`), each writing its own entries; the log Jacobian is
summed over units in their order after the join, so it does not depend on how
the units were split. A unit whose conditional precision does not factor is a
`DomainError`, as it was before the split.
"""
function _jmap_effects(m::CTSEMJointMap, theta::AbstractVector{T},
    xt::AbstractVector) where {T}
    sampler = m.sampler
    laplace = sampler.laplace
    th = _laplace_extend(laplace, theta)
    Ls = _laplace_popchols(th, laplace.spec)
    mu = th[m.P]
    shift = th[m.R] .- m.thref[m.R]
    nunits = sampler.nunits
    us = Vector{Vector{T}}(undef, nunits)
    Fs = Vector{_JMapUnitFactor{T}}(undef, nunits)
    ljs = zeros(T, nunits)
    ok = _laplace_partition(laplace, nunits) do mine, _w
        for U in mine
            d = sampler.udims[U]
            if d == 0
                us[U] = T[]
                Fs[U] = _JMapUnitFactor(CTSEMCholesky(zeros(T, 0, 0), 0, true))
                continue
            end
            unit = try
                _jmap_unit(m, Ls, mu, shift, th, U)
            catch err
                _ctsem_must_propagate(err) && rethrow()
                nothing
            end
            unit === nothing && return false
            F, mean = unit
            us[U] = mean .+ _jmap_usolve(F, view(xt, _sample_urange(sampler, U)))
            Fs[U] = F
            ljs[U] = -_jmap_halflogdet(F)
        end
        return true
    end
    ok || throw(DomainError(theta, "a unit's conditional precision did not factor"))
    logjac = zero(T)
    for k in 1:nunits
        logjac += ljs[k]
    end
    return us, logjac, Fs
end

"""
    CTSEMJointMapDensity(map, count)

The latent density for NUTS: the joint density at the mapped point plus the
map's log Jacobian, with its gradient through the map. One a chain -- `x` and
`gj` are its scratch -- all counting gradients into `count`.
"""
struct CTSEMJointMapDensity{S}
    map::CTSEMJointMap{S}
    x::Vector{Float64}
    gj::Vector{Float64}
    count::Threads.Atomic{Int}
end
CTSEMJointMapDensity(m::CTSEMJointMap, count::Threads.Atomic{Int}) =
    CTSEMJointMapDensity(m, zeros(m.sampler.ndim), zeros(m.sampler.ndim), count)

# Whether latent `xt` maps theta past the random-effect correlation bound
# (sample_nuts.jl): there the wall (`_ctsem_corbound`) is negative. Asked at
# every leapfrog step, so it is a few small factorizations beside the density's
# subject passes. On the race's fixtures nine in ten of the divergent steps
# were past it themselves, deep past, where trajectories that started near a
# correlation the posterior reaches, or anywhere, had crossed it.
function _divergence_at_bound(d::CTSEMJointMapDensity, xt::AbstractVector)
    m = d.map
    laplace = m.sampler.laplace
    isempty(laplace.objective.corbound) && return false
    return try
        theta, _ = _jmap_theta(m, view(xt, 1:m.npar))
        _ctsem_corbound(laplace.objective, _laplace_extend(laplace,
            collect(Float64, theta))) < 0
    catch err
        _ctsem_must_propagate(err) && rethrow()
        false
    end
end

function (d::CTSEMJointMapDensity)(g::AbstractVector{Float64}, xt::AbstractVector)
    Threads.atomic_add!(d.count, 1)
    m = d.map
    sampler = m.sampler
    npar = m.npar
    theta, ltheta = _jmap_theta(m, view(xt, 1:npar))
    mapped = try
        _jmap_effects(m, theta, xt)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    if mapped === nothing
        fill!(g, 0.0)
        return -Inf
    end
    us, logjac, Fs = mapped
    x = d.x
    x[1:npar] .= theta
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) || (x[r] .= us[U])
    end
    v = ctsem_sample_density!(d.gj, sampler, x)
    if !isfinite(v)
        fill!(g, 0.0)
        return -Inf
    end
    gj = d.gj
    gus = [gj[_sample_urange(sampler, U)] for U in 1:sampler.nunits]
    gtheta = gj[1:npar]
    _jmap_theta_adjoint!(gtheta, m, theta, xt, gus, us, Fs)
    fill!(g, 0.0)
    g[1:npar] .= _jmap_theta_pullback(m, view(xt, 1:npar), gtheta)
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) || (g[r] .= _jmap_utsolve(Fs[U], gus[U]))
    end
    return v + logjac + ltheta
end

# The effects' and the log Jacobian's dependence on theta at fixed latent z,
# contracted with the joint's gradient in the effects `gus`, added into
# `gtheta`: one reverse pass through each unit's factor (`_jmap_unit_adjoint!`),
# the units on the pool, each worker accumulating the cotangents of the levels'
# roots, the effects' means `theta[P]` and the coupling's `theta[R]`, which are
# then carried to theta once: the roots through `_laplace_level_chol_derivatives`
# (the joint density's own derivatives of them), the rest directly. This was 37
# forward-mode passes of every unit on the affect pilot (435 positions a dual of
# 12 lanes), 172 s a gradient against the joint density's 1.1 s
# (`review/PILOT-sampling-and-updating-2026-11-06.md`).
function _jmap_theta_adjoint!(gtheta::Vector{Float64}, m::CTSEMJointMap,
    theta::Vector{Float64}, xt::AbstractVector, gus::Vector{Vector{Float64}},
    us::Vector{Vector{Float64}}, Fs)
    sampler = m.sampler
    laplace = sampler.laplace
    spec = laplace.spec
    npar = m.npar
    th = _laplace_extend(laplace, theta)
    Ls = _laplace_popchols(th, spec)
    mu = th[m.P]
    shift = th[m.R] .- m.thref[m.R]
    nP = length(m.P)
    maximum(sampler.udims; init=0) == 0 && return gtheta
    # Dense scratch for the units of one block; a nested unit's pass allocates
    # on its blocks (`_jmap_unit_adjoint_blocked!`).
    dmax = maximum((sampler.udims[U] for U in 1:sampler.nunits
        if !_jmap_blocked(Fs[U])); init=0)
    Lbar = [zeros(size(L)) for L in Ls]
    mubar = zeros(nP)
    # The couplings' cotangent by raw position: `theta[R]` for every member,
    # and each member's own values.
    cbar = zeros(npar)
    coupled = !isempty(m.R) || any(!isempty, m.own)
    reduction = ReentrantLock()
    _laplace_partition(laplace, sampler.nunits) do mine, _w
        work = _JMapAdjointWork(nP, coupled ? npar : 0, dmax, Ls)
        for U in mine
            sampler.udims[U] == 0 && continue
            _jmap_unit_adjoint!(work, m, Ls, mu, shift, th, U, Fs[U], us[U], gus[U],
                view(xt, _sample_urange(sampler, U)))
        end
        Base.@lock reduction begin
            for l in eachindex(Lbar)
                Lbar[l] .+= work.Lbar[l]
            end
            mubar .+= work.mubar
            coupled && (cbar .+= work.Rbar)
        end
        return true
    end
    @inbounds for a in 1:nP
        p = m.P[a]
        p <= npar && (gtheta[p] += mubar[a])
    end
    coupled && (gtheta .+= cbar)
    dL = _laplace_level_chol_derivatives(th, spec)
    @inbounds for l in eachindex(spec.levels)
        level = spec.levels[l]
        k = nrandomeffects(level)
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        Lb = Lbar[l]
        for (t, j) in enumerate(sampler.positions[l])
            j <= npar || continue
            D = dL[l][t]
            acc = 0.0
            for q in 1:r, p in 1:k
                acc += D[p, q] * Lb[p, q]
            end
            gtheta[j] += acc
        end
    end
    return gtheta
end

# A worker's scratch for `_jmap_unit_adjoint!`, sized for the largest unit, and
# its accumulators: `Rbar` the couplings' cotangent by raw position, as long as
# theta where any member is coupled and empty otherwise.
struct _JMapAdjointWork
    X::Matrix{Float64}
    S::Matrix{Float64}
    XP::Matrix{Float64}
    J::Matrix{Float64}
    LJ::Matrix{Float64}
    Jb::Matrix{Float64}
    y::Vector{Float64}
    w::Vector{Float64}
    v::Vector{Float64}
    mean::Vector{Float64}
    lin::Vector{Float64}
    linb::Vector{Float64}
    Lbar::Vector{Matrix{Float64}}
    mubar::Vector{Float64}
    Rbar::Vector{Float64}
end
_JMapAdjointWork(nP::Int, nR::Int, d::Int, Ls) = _JMapAdjointWork(zeros(d, d), zeros(d, d),
    zeros(d, d), zeros(nP, d), zeros(nP, d), zeros(nP, d), zeros(d), zeros(d), zeros(d),
    zeros(d), zeros(nP), zeros(nP), [zeros(size(L)) for L in Ls], zeros(nP), zeros(nR))

# Unit U's reverse pass. Forward (`_jmap_unit`, `_jmap_effects`):
#
#     Q = I + sum_s J_s' Lam_s J_s,   c = sum_s J_s' lin_s,
#     lin_s = h_s - Lam_s mu - G_s shift,   Q = U'U,
#     u = Q^-1 c + U^-1 z,   logjac = -sum_i log U_ii.
#
# With g the joint's gradient in u, `y = U^-1 z`, `w = U^-T g`, `v = Q^-1 g`,
# `mean = u - y`, and `X = U^-1`, the symmetric gradient in Q is
#
#     S = X P X' - (v mean' + mean v') / 2,
#     P_ij = -(z_max(i,j) w_min(i,j)) / 2 off the diagonal, -(z_i w_i + 1) / 2 on it,
#
# the first term the Cholesky factor's adjoint (Murray 2016) for `g'y` and the
# log Jacobian -- with `L = U'` and `Lbar` lower, `L'Lbar` restricted to its
# lower triangle is `-z w'` there less the identity, since `U y = z` -- and the
# second `g'Q^-1 c`'s. The gradient in c is v. Then per member
# `Jbar_s = 2 Lam_s J_s S + lin_s v'`, `linbar_s = J_s v`, so
# `mubar -= Lam_s linbar_s`, `shiftbar_s -= G_s' linbar_s` (at `R` and the
# member's own positions), and the levels' roots take `Jbar_s` back at the
# positions `_jmap_jacobian!` put them.
function _jmap_unit_adjoint!(work::_JMapAdjointWork, m::CTSEMJointMap, Ls, mu, shift,
    th, U::Int, Fu::_JMapUnitFactor{Float64}, u::Vector{Float64}, g::Vector{Float64},
    z::AbstractVector)
    _jmap_blocked(Fu) &&
        return _jmap_unit_adjoint_blocked!(work, m, Ls, mu, shift, th, U, Fu, u, g, z)
    F = Fu.dense
    sampler = m.sampler
    laplace = sampler.laplace
    spec = laplace.spec
    d = sampler.udims[U]
    nP = length(m.P)
    nR = length(m.R)
    coupled = !isempty(work.Rbar)
    X = view(work.X, 1:d, 1:d)
    S = view(work.S, 1:d, 1:d)
    XP = view(work.XP, 1:d, 1:d)
    J = view(work.J, :, 1:d)
    LJ = view(work.LJ, :, 1:d)
    Jb = view(work.Jb, :, 1:d)
    y = view(work.y, 1:d)
    w = view(work.w, 1:d)
    v = view(work.v, 1:d)
    mean = view(work.mean, 1:d)
    lin = work.lin
    linb = work.linb
    _ctsem_cholesky_uinv!(X, F)
    @inbounds for i in 1:d
        a = 0.0
        b = 0.0
        for k in i:d
            a += X[i, k] * z[k]
        end
        for k in 1:i
            b += X[k, i] * g[k]
        end
        y[i] = a
        w[i] = b
    end
    @inbounds for i in 1:d
        a = 0.0
        for k in i:d
            a += X[i, k] * w[k]
        end
        v[i] = a
        mean[i] = u[i] - y[i]
    end
    # P into S, then XP = X P, then S = XP X' less the mean's term.
    @inbounds for j in 1:d, i in 1:d
        S[i, j] = i == j ? -(z[i] * w[i] + 1) / 2 :
            (i > j ? -z[i] * w[j] / 2 : -z[j] * w[i] / 2)
    end
    _ctsem_mul!(XP, X, S)
    _ctsem_mulNT!(S, XP, X)
    @inbounds for j in 1:d, i in 1:d
        S[i, j] -= (v[i] * mean[j] + mean[i] * v[j]) / 2
    end
    for s in m.members[U]
        _jmap_jacobian!(J, m, Ls, s)
        Lam = m.Lam[s]
        _ctsem_mul!(LJ, Lam, J)
        @inbounds for a in 1:nP
            lin[a] = _jmap_lin(m, s, a, mu, shift, th)
            t = 0.0
            for q in 1:d
                t += J[a, q] * v[q]
            end
            linb[a] = t
        end
        _ctsem_mul!(Jb, LJ, S, 2.0, false)
        @inbounds for q in 1:d, a in 1:nP
            Jb[a, q] += lin[a] * v[q]
        end
        @inbounds for a in 1:nP
            acc = 0.0
            for b in 1:nP
                acc += Lam[a, b] * linb[b]
            end
            work.mubar[a] -= acc
        end
        if coupled
            G = m.G[s]
            own = m.own[s]
            @inbounds for b in 1:(nR + length(own))
                acc = 0.0
                for a in 1:nP
                    acc += G[a, b] * linb[a]
                end
                work.Rbar[b <= nR ? m.R[b] : own[b - nR]] -= acc
            end
        end
        # The roots' cotangents, where `_jmap_jacobian!` placed each entry.
        offsets = laplace.units.offsets[U][sampler.flat_member[s]]
        @inbounds for l in eachindex(spec.levels)
            level = spec.levels[l]
            k = nrandomeffects(level)
            r = nlatent(level)
            (k == 0 || r == 0) && continue
            base = offsets[l]
            Lb = work.Lbar[l]
            for p in 1:k
                i = m.row[level.re_index[p]]
                for q in 1:r
                    Lb[p, q] += Jb[i, base + q]
                end
            end
        end
    end
    return work
end

"""Latent coordinates of `(theta, effects x)` under `m`: the map's inverse."""
function _jmap_inverse(m::CTSEMJointMap, theta::AbstractVector, x::AbstractVector)
    sampler = m.sampler
    laplace = sampler.laplace
    xt = zeros(sampler.ndim)
    xt[1:m.npar] .= _jmap_theta_inverse(m, theta)
    th = _laplace_extend(laplace, collect(Float64, theta))
    Ls = _laplace_popchols(th, laplace.spec)
    mu = th[m.P]
    shift = th[m.R] .- m.thref[m.R]
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) && continue
        F, mean = _jmap_unit(m, Ls, mu, shift, th, U)
        xt[r] .= _jmap_umul(F, x[r] .- mean)
    end
    return xt
end

"""`(theta, x)`: the parameters and the joint sampler's vector at latent `xt`."""
function _jmap_point(m::CTSEMJointMap, xt::AbstractVector)
    sampler = m.sampler
    theta, _ = _jmap_theta(m, view(xt, 1:m.npar))
    us, _, _ = _jmap_effects(m, theta, xt)
    x = zeros(sampler.ndim)
    x[1:m.npar] .= theta
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) || (x[r] .= us[U])
    end
    return theta, x
end

_jmap_carry(m::CTSEMJointMap, m2::CTSEMJointMap, xt::AbstractVector) =
    _jmap_inverse(m2, _jmap_point(m, xt)...)

################################################################################
# The members' scores and summaries
################################################################################

_jmap_score!(gsub, laplace, i::Int, aws, shifted, ekf) =
    _laplace_subject_value_gradient!(gsub, laplace.objective.subject_objectives[i],
        aws, shifted, ekf, _HYBRID_FULL)

"""A member's training rows as sums: `b_s`, its score in `b_s` and at its
coupled positions (`_jmap_coupled`), and the cross products the regressions
read."""
mutable struct _JMapStats
    n::Int
    sb::Vector{Float64}
    ss::Vector{Float64}
    st::Vector{Float64}
    bb::Matrix{Float64}
    bs::Matrix{Float64}
    bt::Matrix{Float64}
end
_JMapStats(nP::Int, nR::Int) = _JMapStats(0, zeros(nP), zeros(nP), zeros(nR),
    zeros(nP, nP), zeros(nP, nP), zeros(nP, nR))

"""Every member's empty sums, sized by its coupling."""
_jmap_stats(m::CTSEMJointMap) = [_JMapStats(length(m.P), size(m.G[s], 2))
                                 for s in eachindex(m.G)]

# Members `mine` at the effects `x` (theta extended `th`, roots `Ls`): one row
# each into `stats`, none where the subject is not finite. Into the sums in
# place: a row is a pass over every subject, and copies of the score at the
# coupled positions, and of their outer products, were a gigabyte a row on
# the reanalysis.
function _jmap_scores_range!(stats, m::CTSEMJointMap, flat::Vector{Int}, th, Ls, x, mine)
    sampler = m.sampler
    laplace = sampler.laplace
    aws = _laplace_workspace!(laplace, Float64, length(th))
    ekf = _laplace_ekf_workspace!(laplace, Float64)
    shifted = similar(th)
    gsub = similar(th)
    P = m.P
    R = m.R
    nP = length(P)
    nR = length(R)
    for k in mine
        s = flat[k]
        U = sampler.flat_unit[s]
        mm = sampler.flat_member[s]
        _laplace_member_values!(shifted, th, laplace.spec, Ls,
            view(x, _sample_urange(sampler, U)), laplace.units.offsets[U][mm])
        v = _jmap_score!(gsub, laplace, laplace.units.members[U][mm], aws, shifted, ekf)
        own = m.own[s]
        isfinite(v) || continue
        (all(a -> isfinite(gsub[a]), P) && all(a -> isfinite(gsub[a]), R) &&
            all(a -> isfinite(gsub[a]), own)) || continue
        st = stats[s]
        st.n += 1
        @inbounds for a in 1:nP
            ba = shifted[P[a]]
            st.sb[a] += ba
            st.ss[a] += gsub[P[a]]
            for c in 1:nP
                st.bb[a, c] += ba * shifted[P[c]]
                st.bs[a, c] += ba * gsub[P[c]]
            end
        end
        @inbounds for c in 1:(nR + length(own))
            t = gsub[c <= nR ? R[c] : own[c - nR]]
            st.st[c] += t
            for a in 1:nP
                st.bt[a, c] += shifted[P[a]] * t
            end
        end
    end
    return true
end

function _jmap_scores!(stats, m::CTSEMJointMap, flat::Vector{Int}, x::Vector{Float64})
    laplace = m.sampler.laplace
    th = _laplace_extend(laplace, x[1:m.npar])
    Ls = _laplace_popchols(th, laplace.spec)
    return _laplace_partition((mine, _w) -> _jmap_scores_range!(stats, m, flat, th, Ls, x, mine),
        laplace, length(flat))
end

# Members `mine` of `flat`: each one's curvature in b_s at its unit's mode
# `modes[U]` by central differences of its score, its coupling to its coupled
# positions (`R`, then its own `own[s]`) from the same differences, and the
# linear term that puts the summary's maximum at the mode.
function _jmap_curvature_range!(Lam, h, G, sampler, P, R, own, flat, th, Ls, modes, mine)
    laplace = sampler.laplace
    aws = _laplace_workspace!(laplace, Float64, length(th))
    ekf = _laplace_ekf_workspace!(laplace, Float64)
    shifted = similar(th)
    work = similar(th)
    gplus = similar(th)
    gminus = similar(th)
    g0 = similar(th)
    nP = length(P)
    for k in mine
        s = flat[k]
        U = sampler.flat_unit[s]
        mm = sampler.flat_member[s]
        i = laplace.units.members[U][mm]
        _laplace_member_values!(shifted, th, laplace.spec, Ls, modes[U],
            laplace.units.offsets[U][mm])
        isfinite(_jmap_score!(g0, laplace, i, aws, shifted, ekf)) || continue
        Rs = vcat(R, own[s])
        nR = length(Rs)
        H = zeros(nP, nP)
        Gt = zeros(nR, nP)
        for j in 1:nP
            p = P[j]
            step = 1e-4 * max(1.0, abs(shifted[p]))
            copyto!(work, shifted)
            work[p] += step
            _jmap_score!(gplus, laplace, i, aws, work, ekf)
            work[p] -= 2step
            _jmap_score!(gminus, laplace, i, aws, work, ekf)
            for a in 1:nP
                H[a, j] = -(gplus[P[a]] - gminus[P[a]]) / (2step)
            end
            for a in 1:nR
                Gt[a, j] = -(gplus[Rs[a]] - gminus[Rs[a]]) / (2step)
            end
        end
        (all(isfinite, H) && all(isfinite, Gt)) || continue
        E = eigen(Symmetric((H .+ transpose(H)) ./ 2))
        Lam[s] = E.vectors * Diagonal(max.(E.values, 0.0)) * transpose(E.vectors)
        h[s] = g0[P] .+ Lam[s] * shifted[P]
        G[s] = Matrix(transpose(Gt))
    end
    return true
end

"""
    _jmap_member_update(Lp, Gp, st, thbar, thref)

A member's summary refitted from its sums: least squares of its score in `b`
on `b` in the directions the rows' `b` span, the previous summary `Lp` kept in
the directions they do not; with coupling, its score at its coupled positions
on `b` the same way, for `G`. `thbar` is those positions' mean over the rows,
`thref` their reference. `nothing` where the rows say too little.
"""
function _jmap_member_update(Lp::Matrix{Float64}, Gp::Matrix{Float64}, st::_JMapStats,
    thbar::Vector{Float64}, thref::Vector{Float64}; tol::Real=1e-8)
    n = st.n
    n >= 8 || return nothing
    bbar = st.sb ./ n
    sbar = st.ss ./ n
    Cbb = st.bb ./ n .- bbar * transpose(bbar)
    Cbs = st.bs ./ n .- bbar * transpose(sbar)
    E = eigen(Symmetric((Cbb .+ transpose(Cbb)) ./ 2))
    top = maximum(E.values)
    top > 0 || return nothing
    V = E.vectors[:, E.values .> tol * top]
    g = size(V, 2)
    n >= 4 * (g + 1) || return nothing
    Caa = Symmetric(transpose(V) * Cbb * V)
    Y = -(transpose(transpose(V) * Cbs) / Caa)
    Mg = transpose(V) * Y
    Mg = (Mg .+ transpose(Mg)) ./ 2
    Cg = Y .- V * (transpose(V) * Y)
    Pc = I - V * transpose(V)
    Lam = V * Mg * transpose(V) .+ Cg * transpose(V) .+ V * transpose(Cg) .+ Pc * Lp * Pc
    F = eigen(Symmetric((Lam .+ transpose(Lam)) ./ 2))
    Lam = F.vectors * Diagonal(max.(F.values, 0.0)) * transpose(F.vectors)
    G = Gp
    if size(Gp, 2) > 0
        tbar = st.st ./ n
        Cbt = st.bt ./ n .- bbar * transpose(tbar)
        YT = -(transpose(transpose(V) * Cbt) / Caa)
        G = Pc * Gp .+ V * transpose(YT)
    end
    hnew = sbar .+ Lam * bbar
    size(G, 2) > 0 && (hnew = hnew .+ G * (thbar .- thref))
    return Lam, hnew, G
end

"""The flat indices of the members whose units have coordinates."""
_jmap_flat(m::CTSEMJointMap) = reduce(vcat, m.members; init=Int[])

"""
    _jmap_build(sampler, theta0, cov; cross)

The map from the Laplace fit: every unit's mode at `theta0`, each member's
curvature there, theta whitened by `cov`. Returns the map and its cost in
full passes over the subjects.
"""
function _jmap_build(sampler::CTSEMSampler, theta0::Vector{Float64},
    cov::AbstractMatrix{Float64}; cross::Bool=false)
    laplace = sampler.laplace
    npar = sampler.npar
    _laplace_ensure_pool!(laplace)
    P, row, R, own, members = _jmap_layout(sampler, cross)
    nP = length(P)
    nR = length(R)
    th = _laplace_extend(laplace, theta0)
    Ls = _laplace_popchols(th, laplace.spec)
    modes = [zeros(d) for d in sampler.udims]
    # Units on the pool, as the Laplace evaluation solves them: one thread took
    # 123 s of the affect pilot's 1,324 units from their held modes.
    _laplace_parallel(laplace, 1:sampler.nunits) do U
        sampler.udims[U] == 0 && return true
        _laplace_solve_unit_mode!(laplace, U, th, Ls)
        modes[U] = copy(laplace.modes[U])
        return true
    end
    nflat = length(sampler.flat_unit)
    Lam = [zeros(nP, nP) for _ in 1:nflat]
    h = [zeros(nP) for _ in 1:nflat]
    G = [zeros(nP, nR + length(own[s])) for s in 1:nflat]
    flat = reduce(vcat, members; init=Int[])
    _laplace_partition((mine, _w) -> _jmap_curvature_range!(Lam, h, G, sampler, P, R, own,
        flat, th, Ls, modes, mine), laplace, length(flat))
    S = Symmetric((Matrix{Float64}(cov) .+ transpose(Matrix{Float64}(cov))) ./ 2)
    Ct = Matrix(cholesky(Symmetric(Matrix(S) + 1e-12 * max(tr(S), 1.0) / npar * I)).L)
    m = CTSEMJointMap(sampler, npar, P, row, R, own, members, Lam, h, G, copy(th),
        copy(theta0), Ct, nothing)
    # The differences' passes, and about five for the units' mode solves.
    return m, isempty(flat) ? 0 : 2 * nP + 6
end

"""
    _jmap_refit(m, Wl, rng; maxrows)

The map refitted from latent draws `Wl` (columns, under `m`), at most `maxrows`
of them: theta's centre and whitening, each member's summary. Returns the new
map, the members refitted, theta's shrinkage, and the rows' cost in full passes.
`progress(done, rows)` is called after each row (`_ctsem_call_hook`): a row is
a pass over every subject, and on the reanalysis 400 of them ran for hours
with nothing said.
"""
function _jmap_refit(m::CTSEMJointMap, Wl::AbstractMatrix, rng; maxrows::Int=400,
    progress=nothing)
    sampler = m.sampler
    npar = m.npar
    n = size(Wl, 2)
    pick = n > maxrows ? sort(randperm(rng, n)[1:maxrows]) : collect(1:n)
    W = Wl[:, pick]
    nr = length(pick)
    thetas = _jmap_thetas(m, W[1:npar, :])
    shrunk = _smc_shrunk_cov(Matrix{Float64}(thetas), ones(nr), m.Ct * transpose(m.Ct))
    theta0 = vec(mean(thetas; dims=2))
    S = shrunk.cov
    Ct = Matrix(cholesky(Symmetric(S + 1e-12 * max(tr(S), 1.0) / npar * I)).L)
    nP = length(m.P)
    flat = _jmap_flat(m)
    stats = _jmap_stats(m)
    rows = 0
    if nP > 0 && !isempty(flat)
        for t in 1:nr
            _, x = _jmap_point(m, view(W, :, t))
            _jmap_scores!(stats, m, flat, x)
            _ctsem_call_hook(progress, t, nr)
        end
        rows = nr
    end
    Lam = copy(m.Lam)
    h = copy(m.h)
    G = copy(m.G)
    refit = 0
    for s in flat
        cpos = _jmap_coupled(m, s)
        upd = _jmap_member_update(m.Lam[s], m.G[s], stats[s], theta0[cpos], m.thref[cpos])
        upd === nothing && continue
        Lam[s], h[s], G[s] = upd
        refit += 1
    end
    return _jmap_with(m; Lam=Lam, h=h, G=G, theta0=theta0, Ct=Ct), refit, shrunk.lambda, rows
end

"""Transitions in a row that moved nothing (divergent, or accepting under
`_STALL_ACCEPT` on average) before a learning window's chain finds its step
size again (`_ChainRetune`). At a step size tuned to accept 0.8 on average,
five in a row are rare enough that a retune of a moving chain costs little
more than its few leapfrogs and settle; each stalled transition on the
reanalysis cost about 4.5 minutes."""
const _JMAP_STALL = 5

"""Which chains moved over their last window: one whose draws spread, in
theta's latent coordinates, under a tenth of the pooled spread in most of them
has stopped. All stopped reads as all moving: there is then none to restart
from."""
function _jmap_moved(results::Vector{_ChainResult}, npar::Int)
    W = reduce(hcat, [r.draws[1:npar, :] for r in results])
    pooled = max.(vec(std(W; dims=2)), 1e-300)
    moved = [median(vec(std(r.draws[1:npar, :]; dims=2)) ./ pooled) >= 0.1 for r in results]
    any(moved) || fill!(moved, true)
    return moved
end

"""The chains' step sizes over a window, smallest and largest, and the deepest
of their mean tree depths: a window lasts as long as its slowest chain, and on
a simulated analogue of the hc-wellbeing reanalysis three chains finished a
window 15 minutes after the first."""
function _jmap_chain_spread(results::Vector{_ChainResult})
    steps = [r.stepsize for r in results]
    depth = maximum(r -> isempty(r.depth) ? 0.0 : sum(r.depth) / length(r.depth), results)
    return (step_min=minimum(steps), step_max=maximum(steps), depth_max=Float64(depth))
end

"""Each chain's start under `m2`: a chain that moved goes on from its last
state, a stopped one from a random draw of those that moved."""
function _jmap_restarts(m::CTSEMJointMap, m2::CTSEMJointMap, results, moved, Wl, rng)
    return [moved[c] ? _jmap_carry(m, m2, results[c].x) :
            _jmap_carry(m, m2, Wl[:, rand(rng, 1:size(Wl, 2))]) for c in eachindex(results)]
end

################################################################################
# Reporting the effects
################################################################################

# For each unit, each of the report layout's blocks (`report.units.blocks[U]`)
# matched to the sampled objective's block of the same level and members, 0
# where the sampled objective has none -- a unit whose effects are all
# integrated keeps no block there. The two objectives are built from one fit,
# so their units and members agree; refused if they do not.
function _jmap_report_match(report::CTSEMLaplaceObjective, laplace::CTSEMLaplaceObjective)
    length(report.units.members) == length(laplace.units.members) ||
        throw(ArgumentError("the report objective's units are not the sampled one's"))
    match = Vector{Vector{Int}}(undef, length(report.units.members))
    for U in eachindex(report.units.members)
        report.units.members[U] == laplace.units.members[U] ||
            throw(ArgumentError("the report objective's members are not the sampled one's"))
        mine = laplace.units.blocks[U]
        match[U] = [something(findfirst(b -> b.level == blk.level && b.members == blk.members,
            mine), 0) for blk in report.units.blocks[U]]
    end
    return match
end

# One draw's natural deviations in the report layout's coordinates: each
# block's coordinate part from the sampled objective's matching block, through
# that objective's roots `Ls` -- zero where all of the block's effects are
# integrated, whose part is added from their conditional afterwards.
function _jmap_layer_deviations!(dev, report, laplace, sampler, Ls, x, devstarts, match)
    fill!(dev, 0.0)
    @inbounds for U in eachindex(report.units.members)
        ubase = sampler.uoffsets[U]
        for (b, hb) in enumerate(match[U])
            hb == 0 && continue
            blk = laplace.units.blocks[U][hb]
            L = Ls[blk.level]
            s = devstarts[U][b]
            for p in axes(L, 1)
                acc = 0.0
                for q in axes(L, 2)
                    acc += L[p, q] * x[ubase + blk.offset + q]
                end
                dev[s + p] = acc
            end
        end
    end
    return dev
end

"""
    _sample_report_effects(laplace, report, sampler, results, nchains, ndraws, npar,
        save_effects, seed, verbose, point)

A joint sampler's kept draws as theta and natural deviations, the one report
of the joint target whichever sampler drew it. `point(c, t)` gives chain `c`'s
draw `t` as `(theta, x)`, `x` the joint sampler's vector on `laplace`. The
deviations are laid out on `report`, the layer of the same fit (every effect a
coordinate; `laplace` itself where `report` is `nothing`): each block's
coordinate part from `laplace`'s matching block, and where `laplace`'s filter
integrates effects, those drawn from their exact conditional given the draw
(`_jmap_hybrid_range!`), so a sampled fit reports every effect. Returns
`(pooled, effects, hybrid)`: the pooled draws with theta (and the deviations,
under `save_effects`), the effects' summary, and whether any were integrated.
"""
function _sample_report_effects(laplace::CTSEMLaplaceObjective, report, sampler,
    results, nchains::Int, ndraws::Int, npar::Int, save_effects::Bool, seed::Integer,
    verbose::Bool, point)
    layout = report === nothing ? laplace : report
    match = _jmap_report_match(layout, laplace)
    devstarts, ndev = _laplace_deviation_offsets(layout)
    kept = save_effects ? npar + ndev : npar
    pooled = _sample_pool(results, nchains, ndraws, npar; rows=kept)
    draws = pooled.draws
    total = nchains * ndraws
    effect_sum = zeros(Float64, ndev)
    effect_sq = zeros(Float64, ndev)
    blocks = _laplace_deviation_blocks(layout)
    effect_cross = zeros(Float64, _effect_cross_length(blocks))
    hyb = laplace.objective.params.hybrid
    hybrid = _hybrid_on(hyb)
    lv = laplace.spec.levels[1]
    mblocks = hybrid ? _jmap_member_blocks(layout, devstarts) : nothing
    dev = zeros(Float64, ndev)
    failed = 0
    for c in 1:nchains, t in 1:ndraws
        column = (c - 1) * ndraws + t
        theta, x = point(c, t)
        th = _laplace_extend(laplace, theta)
        Ls = _laplace_popchols(th, laplace.spec)
        _jmap_layer_deviations!(dev, layout, laplace, sampler, Ls, x, devstarts, match)
        if hybrid && ndev > 0
            B = _hybrid_dev_map(hyb, lv, th)
            ok = _laplace_partition((mine, _w) -> _jmap_hybrid_range!(dev, laplace, th,
                Ls, x, sampler, B, mblocks, seed, column, mine), laplace, sampler.nunits)
            ok || (failed += 1)
        end
        @inbounds for j in 1:npar
            draws[j, column] = theta[j]
        end
        @inbounds for j in 1:ndev
            effect_sum[j] += dev[j]
            effect_sq[j] += dev[j] * dev[j]
            save_effects && (draws[npar + j, column] = dev[j])
        end
        _effect_cross_add!(effect_cross, dev, blocks)
    end
    failed > 0 && verbose && println(_console(), "  ", failed,
        " draw(s) had a subject whose integrated effects could not be evaluated")
    return (pooled=pooled,
        effects=_effect_summary(effect_sum, effect_sq, effect_cross, total, blocks),
        hybrid=hybrid)
end

# Each member's subject-level block start in the natural deviations, -1 where
# it has none.
function _jmap_member_blocks(laplace::CTSEMLaplaceObjective, devstarts)
    out = [fill(-1, length(laplace.units.members[U])) for U in eachindex(laplace.units.members)]
    for U in eachindex(laplace.units.members), (b, blk) in enumerate(laplace.units.blocks[U])
        blk.level == 1 && length(blk.members) == 1 && (out[U][blk.members[1]] = devstarts[U][b])
    end
    return out
end

# Units `mine`: the filter-integrated deviations of every member, drawn from
# their conditional given the draw (`_hybrid_member_modes`), added through `B`
# into `dev`. Items a member's pass did not integrate take their prior. False
# where a member could not be evaluated.
function _jmap_hybrid_range!(dev, laplace, th, Ls, x, sampler, B, mblocks, seed, column, mine)
    aws = _laplace_workspace!(laplace, Float64, length(th))
    nh = size(B, 2)
    re = laplace.spec.levels[1].re_index
    ok = true
    for U in mine
        hook = _hybrid_member_modes(laplace, U, th, Ls,
            Vector{Float64}(view(x, _sample_urange(sampler, U))), aws)
        if hook === nothing
            ok = false
            continue
        end
        rng = Random.Xoshiro(hash((seed, column, U)))
        for mm in eachindex(laplace.units.members[U])
            start = mblocks[U][mm]
            start < 0 && continue
            d = randn(rng, nh)
            q = hook.q[mm]
            if q > 0
                F = _ctsem_cholesky(copy(hook.A[mm]), q)
                if F.ok
                    e = _jmap_usolve(F, d[1:q])
                    full = zeros(nh)
                    seen = falses(nh)
                    for (c, item) in enumerate(hook.items[mm])
                        full[item] = hook.dhat[mm][c] + e[c]
                        seen[item] = true
                    end
                    for j in 1:nh
                        seen[j] || (full[j] = d[j])
                    end
                    d = full
                end
            end
            for p in axes(B, 1)
                acc = 0.0
                for j in 1:nh
                    acc += B[p, j] * d[j]
                end
                dev[start + p] += acc
            end
        end
    end
    return ok
end

################################################################################
# The sampler
################################################################################

"""
    _JMapLearnReport(reporter, callback, offset, nwarmup, total, logp, retuned)

The learning windows' progress, reported as one phase counted over the whole
schedule (`total` iterations a chain). It is the callback every chain of each
window calls (`_run_chains`, rate-limited there), and maps that window's own
counts -- warm-up from one, then draws from one -- onto the schedule:
`offset` is the iterations a chain before this window, `nwarmup` this
window's warm-up. It prints on `reporter` and hands the caller's callback the
same five fields under the phase `"learning"`, which is all a worker process
reports through (`.ctBackendProgressFileWriter()`). Reported as they came,
the counts would restart eight times; reported not at all, as until
2026-11-18, the line was blank from the sampler's start to its kept draws --
hours on a hard model -- and a run in worker processes showed nothing then.

Every chain reports (`_JMapChainWatch`, `_run_chains`' `chain_callback`), one
at a time under `lock`, and the count is the chains' mean position: a window
ends when its last chain does, and chain 1 alone, as it reported until
2026-11-22, fell silent whenever it finished first -- 15 minutes on a
simulated analogue of the hc-wellbeing reanalysis while three chains finished
window 2, and on the reanalysis itself hours that read as a silent refit.
Beside it the line gives the log density of the chain furthest behind, the
window's divergences over every chain, each chain's iterations into the
window where they differ, and the step sizes retuned where a chain stalled
(`retuned`, `_ChainRetune`), which no log density shows of a chain that is
not the one printed.
"""
mutable struct _JMapLearnReport
    reporter::CTSEMProgress
    callback::CTSEMCallback
    offset::Int
    nwarmup::Int
    total::Int
    logp::Float64
    retuned::Threads.Atomic{Int}
    at::Vector{Int}
    lps::Vector{Float64}
    divs::Vector{Int}
    lock::ReentrantLock
end

_JMapLearnReport(reporter::CTSEMProgress, callback::CTSEMCallback, total::Int,
    nchains::Int) = _JMapLearnReport(reporter, callback, 0, 0, total, NaN,
        Threads.Atomic{Int}(0), zeros(Int, nchains), fill(NaN, nchains),
        zeros(Int, nchains), ReentrantLock())

"""A window's start: every chain at the schedule's `offset`, none diverged."""
function _jmap_learn_window!(r::_JMapLearnReport, nwarmup::Int)
    r.nwarmup = nwarmup
    fill!(r.at, r.offset)
    fill!(r.divs, 0)
    return r
end

function (r::_JMapLearnReport)(c::Int, phase, iteration, total, logp, divergent)
    Base.@lock r.lock begin
        r.at[c] = r.offset + Int(iteration) + (phase == "sampling" ? r.nwarmup : 0)
        r.lps[c] = Float64(logp)
        r.divs[c] = Int(divergent)
        at = sum(r.at) ÷ length(r.at)
        k = argmin(r.at)
        r.logp = r.lps[k]
        div = sum(r.divs)
        n = r.retuned[]
        lag = maximum(r.at) > minimum(r.at)
        # Printed once the chain shown has finished a transition.
        !isnan(r.logp) && _due(r.reporter) && _progress_line(r.reporter, at, r.total,
            @sprintf("logp %11.2f", r.logp), @sprintf("div %d", div),
            lag ? "chains " * join(r.at .- r.offset, " ") : "",
            n > 0 ? @sprintf("retuned %d", n) : "")
        _invoke_callback(r.callback, "learning", at, r.total, r.logp, div)
    end
    return nothing
end

"""Chain `c`'s callback into the learning report."""
struct _JMapChainWatch
    report::_JMapLearnReport
    chain::Int
end

(w::_JMapChainWatch)(phase, iteration, total, logp, divergent) =
    w.report(w.chain, phase, iteration, total, logp, divergent)

"""
    _JMapStageReport(phase, label, verbose, overwrite, sink, callback, logp)

A window-end refit stage's progress, `(done, total)` at a time: its rows
(`"scoring"`, `_jmap_refit`) or theta's triangular fit (`"mapping"`,
`_jmap_learn_tri`), printed as its own line and handed to the caller's
callback under that phase with the learning line's last log density (the
furthest-behind chain's), which is what the
line and a worker process report until the chains run again. After window 2
on the reanalysis the two ran 13 hours with nothing printed.
"""
mutable struct _JMapStageReport
    reporter::CTSEMProgress
    callback::CTSEMCallback
    phase::String
    logp::Float64
end

_JMapStageReport(phase::String, label::String, verbose::Bool, overwrite::Bool, sink,
    callback, logp::Float64) =
    _JMapStageReport(CTSEMProgress(verbose; label=label, overwrite=overwrite, sink=sink),
        CTSEMCallback(callback; every=1.0), phase, logp)

function (r::_JMapStageReport)(done, total)
    _due(r.reporter) && _progress_line(r.reporter, done, total)
    _invoke_callback(r.callback, r.phase, done, total, r.logp, 0)
    return nothing
end

"""
    _jmap_final(m, counter, nchains, parallel, seed, centre, metric, nwarmup,
        ndraws, ...)

The kept draws under the final map `m`: a re-tuning warm-up and `ndraws` a
chain from `starts`, then extended by the one stopping rule
(`_sample_until_target`) judged on theta rather than on its latent
coordinates. `_sample_to_target` with the map's pooling, in a function of its
own so that its continuation captures nothing that is assigned again.
"""
function _jmap_final(m::CTSEMJointMap, counter::Threads.Atomic{Int}, nchains::Int,
    parallel::Bool, seed::Integer, centre::Vector{Float64}, metric::CTSEMMetric,
    nwarmup::Int, ndraws::Int, maxdepth::Int, target_accept::Float64,
    maxdelta::Float64, init_scale::Float64, adapt_metric::Bool, adapt, eps::Float64,
    starts, min_ess::Float64, mean_ess::Float64, max_draws::Int, rhat_target::Float64,
    verbose::Bool, overwrite::Bool, progress_callback, progress_sink, t0::Float64,
    coordinate)
    npar = m.npar
    held = Ref(_sample_chains(nchains, parallel, seed, centre, metric, nwarmup, ndraws,
        maxdepth, target_accept, maxdelta, init_scale, adapt_metric, adapt,
        c -> CTSEMJointMapDensity(m, counter); init_eps=eps, starts=starts,
        progress=verbose, overwrite=overwrite, progress_callback=progress_callback,
        progress_sink=progress_sink))
    extend! = function (wanted, attempt)
        held[] = _merge_chains(held[],
            _continue_chains(nchains, parallel, seed + 1000 * attempt, wanted, maxdepth,
                maxdelta, c -> CTSEMJointMapDensity(m, counter), held[]; progress=verbose,
                overwrite=overwrite, progress_callback=progress_callback,
                progress_sink=progress_sink))
        return nothing
    end
    run = _sample_until_target(() -> _jmap_thetas(m, _pool_draws(held[], npar)),
        extend!, nchains, ndraws, min_ess, mean_ess, max_draws, rhat_target, verbose;
        t0=t0, coordinate=coordinate)
    return (results=held[], ndraws=run.total, trace=run.trace)
end

"""
The map sampler after its first learning window (`ctsem_sample_map` with
`first_only`): the window's chains, for the same call to go on from
(`resume`), and where each chain's draws are centred. `distance` is a chain
mean's squared distance from the estimate, in theta's latent coordinates under
the starting map, whitened by the Hessian the caller handed in, and in the
wider of that Hessian's normal and the window's own spread
(`_jmap_window_distances`); NaN for a chain that stopped, and with no Hessian,
whose whitening is not the estimate's marginal covariance. A mean past `bar`
lies where neither reaches: the posterior's mass is centred away from the
estimate, which no curvature at the estimate can show (R's check for a higher
maximum, `.ctBackendMaximumCheck()`, climbs the fit's objective from
`means`).
`map` is the starting map and `extra` its cost (`_jmap_build`), kept so that
going on neither builds it again nor differs by what a second build's mode
solves would.
"""
struct CTSEMMapStart{M}
    map::M
    extra::Int
    results::Vector{_ChainResult}
    gradients::Int
    secs::Float64
    bar::Float64
    distance::Vector{Float64}
    means::Matrix{Float64}
    flagged::Vector{Bool}
end

function _jmap_start(m::CTSEMJointMap, extra::Int, results::Vector{_ChainResult}, moved,
    gradients::Int, secs::Float64, bar::Float64, metric::Bool)
    npar = m.npar
    means = zeros(npar, length(results))
    distance = fill(NaN, length(results))
    windows = Matrix{Float64}[Matrix{Float64}(view(r.draws, 1:npar, :)) for r in results]
    for c in eachindex(results)
        means[:, c] .= vec(sum(_jmap_thetas(m, windows[c]); dims=2)) ./ size(windows[c], 2)
    end
    if metric
        within = _jmap_window_distances(windows, collect(Bool, moved))
        for c in eachindex(results)
            moved[c] && (distance[c] = within[c])
        end
    end
    flagged = [bar > 0 && isfinite(d) && d > bar for d in distance]
    return CTSEMMapStart(m, extra, results, gradients, secs, bar, distance, means, flagged)
end

"""
Each chain's squared distance from the estimate over a window, from its draws
in latent coordinates whitened at the estimate (`windows`, a matrix a chain):
its mean against the wider, direction by direction, of the Hessian's normal
there (the identity) and the window's spread within the chains that moved
(`spread`), i.e. `m' (V max(L, 1) V')^-1 m` for that spread's eigenvectors `V`
and eigenvalues `L`. The bar on it (R's `.ctMaximumBar()`, the 0.999 point of a
draw of the Hessian's normal) holds for the mean of any number of a chain's
draws while the posterior is no wider than that normal; where it is wider
along a direction the chains cross slowly, the mean is near one draw of it and
the Hessian's metric alone put it past the bar often. On test-julia-sample.R's
random-drift fixture (centred at the estimate, sds up to 4.6 times the
Hessian's) that metric flagged 23 runs of 80 and this one none, its largest
26.1 of 27.9; on the affect pilot's window (a higher maximum 41.6 nats up) the
two chains read about 2,400 and 3,000 of 950 here against 4,500 and 5,400
(review/OPEN.md, 2026-11-22). The spread between chains is left out, so that
one chain that finds another basin while the others stay counts as found.
"""
function _jmap_window_distances(windows::Vector{Matrix{Float64}}, spread::AbstractVector{Bool})
    means = Vector{Float64}[vec(sum(W; dims=2)) ./ size(W, 2) for W in windows]
    distance = Float64[sum(abs2, w) for w in means]
    blocks = Matrix{Float64}[]
    for c in eachindex(windows)
        n = size(windows[c], 2)
        (spread[c] && n > 1) || continue
        push!(blocks, (windows[c] .- means[c]) ./ sqrt(n - 1))
    end
    isempty(blocks) && return distance
    X = reduce(hcat, blocks) ./ sqrt(length(blocks))
    all(isfinite, X) || return distance
    F = svd(X)
    for c in eachindex(means)
        z = transpose(F.U) * means[c]
        inside = 0.0
        scaled = 0.0
        for j in eachindex(z)
            inside += z[j]^2
            scaled += z[j]^2 / max(F.S[j]^2, 1.0)
        end
        distance[c] = max(distance[c] - inside, 0.0) + scaled
    end
    return distance
end

"""
    ctsem_map_centre(start)

What R reads of a `CTSEMMapStart`: the bar, each chain's squared distance and
whether it is past the bar, the chains' mean theta (a column a chain), and
what the window cost.
"""
ctsem_map_centre(s::CTSEMMapStart) = (bar=s.bar, distance=s.distance, means=s.means,
    flagged=s.flagged, gradients=s.gradients, secs=s.secs,
    draws=size(s.results[1].draws, 2), nchains=length(s.results))

export ctsem_map_centre

"""
    ctsem_sample_map(laplace, values; kwargs...)

Sample the joint posterior through a learnt map (this file's header). `laplace`
is the hybrid objective where the model has effects the filter integrates, the
layer otherwise; the keywords are `ctsem_sample`'s, plus the map's: `windows`
(draws a chain in each learning window, after `nwarmup` warm-up iterations on
the map from the fit), `nsettle` (the re-tuning warm-up after each refit and
before the kept draws), `ncheck` (a check window under the final map),
`maxrows` (the most draws a refit reads), `cross` (each member's summary
coupled to the free parameters outside its effects), `report` (the layer's
objective, whose layout the effects are reported in: a unit whose effects are
all integrated has no block on the hybrid objective), `theta_map` (theta's
triangular map learnt; false keeps its whitening alone).

`first_only` returns after the first window's draws, before anything is learnt
from them, as a `CTSEMMapStart`: where each chain is centred against
`centre_bar` (`ctsem_map_centre`). The same call with that as `resume` goes on
from those draws, to the draws an uninterrupted call returns. So a caller can
act on the first window before the rest is paid for, and a run it lets go on
loses nothing.

Returns `ctsem_sample`'s record, with the draws of theta and, with
`save_effects`, of the natural deviations, plus the map's: `gradients` (every
density evaluation and score pass, in gradients), the learning windows, and
the kept draws' latent mean and sd per coordinate (`latent_mean`,
`latent_sd`, theta's first) and the most parents a component of theta's map
conditions on (`theta_parents`, -1 without a map).
"""
function ctsem_sample_map(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    npar::Integer=length(values), nchains::Integer=4, nwarmup::Integer=75,
    ndraws::Integer=500, maxdepth::Integer=10, target_accept::Real=0.8,
    maxdelta::Real=1000.0, seed::Integer=20260828, init_scale::Real=1.0,
    adapt_metric::Bool=false, save_effects::Bool=false,
    hessian::Union{Nothing,AbstractMatrix}=nothing, verbose::Bool=false,
    min_ess::Real=0.0, mean_ess::Real=0.0, max_draws::Integer=0,
    rhat_target::Real=1.01, progress_overwrite::Bool=true, progress_callback=nothing,
    progress_sink=nothing, stepsize::Real=0.0, coordinate=nothing,
    windows::Union{Nothing,AbstractVector{<:Integer}}=nothing, nsettle::Integer=25,
    ncheck::Integer=25, maxrows::Integer=400, cross::Bool=true,
    report::Union{Nothing,CTSEMLaplaceObjective}=nothing, theta_map::Bool=true,
    first_only::Bool=false, centre_bar::Real=0.0,
    resume::Union{Nothing,CTSEMMapStart}=nothing)

    t0 = time()
    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    _check_sample_arguments(nchains, ndraws, target_accept; nwarmup=nwarmup)
    sampler = ctsem_sampler(laplace, npar)
    npar = sampler.npar
    D = sampler.ndim
    theta0 = collect(Float64, values)[1:npar]
    width = _laplace_ensure_pool!(laplace)
    parallel = nchains > 1 && Threads.nthreads() > 1

    # Theta's whitening: the marginal covariance, from the Hessian the caller
    # hands in (the fit's, rescaled by its check where it has one), its flat
    # directions bounded by the prior; the conditional block where none is given.
    prec = _ctsem_prior_precision(laplace.objective, npar)
    # Names of their own here: `th` and `Ls` again would make the draw loop's
    # below, which a closure there captures, function locals assigned twice,
    # each a `Core.Box`.
    # A call going on from its first window (`resume`) takes that window's map.
    built = if resume !== nothing
        resume.map.sampler.laplace === laplace ||
            throw(ArgumentError("resume is another objective's first window"))
        (resume.map::CTSEMJointMap{typeof(sampler)}, resume.extra)
    elseif hessian === nothing
        th0 = _laplace_extend(laplace, theta0)
        Ls0 = _laplace_popchols(th0, laplace.spec)
        _laplace_parallel(laplace, 1:sampler.nunits) do U
            _laplace_solve_unit_mode!(laplace, U, th0, Ls0)
            return true
        end
        _jmap_build(sampler, theta0,
            _conditional_population_covariance(sampler, th0, Ls0, nothing); cross=cross)
    else
        H = Matrix{Float64}(hessian)[1:npar, 1:npar]
        _jmap_build(sampler, theta0,
            _prior_bounded_covariance(Symmetric(-(H .+ transpose(H)) ./ 2), prec); cross=cross)
    end
    m, extra = built
    counter = Threads.Atomic{Int}(0)

    ranges = UnitRange{Int}[1:npar]
    covs = Matrix{Float64}[Matrix{Float64}(I, npar, npar)]
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) && continue
        push!(ranges, r)
        push!(covs, Matrix{Float64}(I, length(r), length(r)))
    end
    metric = _metric_from_covariances(ranges, covs)
    adapt = [b == 1 for b in eachindex(ranges)]
    (verbose && resume === nothing) && println(_console(), "Sampling through a learnt map: ",
        nchains, " chain(s), ",
        D, " dimensions (", npar, " population + ", D - npar, " effects), ",
        width, " worker(s) shared")

    rng = Random.Xoshiro(UInt64(seed) + 0x6a6d6170)
    centre = zeros(D)
    starts = nothing
    eps = Float64(stepsize)
    windowlog = NamedTuple[]
    # Theta's triangular map is learnt from the windows' pooled draws, over a
    # third window.
    learn = theta_map
    windows = windows === nothing ? (learn ? [40, 80, 160] : [40, 80]) : windows
    checking = ncheck > 0 && !isempty(windows)
    # The windows and the check window, reported as one phase; handed to the
    # chains only where something reads it. A loop, not a generator, which
    # would capture `nwarmup`, assigned twice here, in a box.
    nlearn = checking ? Int(nsettle) + Int(ncheck) : 0
    for (w, nw) in enumerate(windows)
        nlearn += (w == 1 ? nwarmup : Int(nsettle)) + Int(nw)
    end
    learning = _JMapLearnReport(CTSEMProgress(verbose; label="learning the map",
            overwrite=progress_overwrite, sink=progress_sink),
        CTSEMCallback(progress_callback; every=0.0), nlearn, nchains)
    watch = (verbose || progress_callback !== nothing) ?
        (c -> _JMapChainWatch(learning, c)) : nothing
    # Every window's draws are warm-up, so a chain that stalls in one finds its
    # step size again where it stands (`_ChainRetune`).
    retune = _ChainRetune(_JMAP_STALL, Int(nsettle), Float64(target_accept),
        learning.retuned)
    pooled = zeros(npar, 0)
    # The first window, where the call is asked for that alone or goes on from
    # one (`CTSEMMapStart`).
    (first_only && isempty(windows)) && throw(ArgumentError(
        "first_only returns after the first learning window, and windows is empty"))
    if resume !== nothing
        (!isempty(windows) && length(resume.results) == nchains &&
            size(resume.results[1].draws) == (D, Int(windows[1]))) ||
            throw(ArgumentError("resume is not this call's first window"))
        Threads.atomic_add!(counter, resume.gradients - extra)
    end
    for (w, nw) in enumerate(windows)
        mw = m
        _jmap_learn_window!(learning, w == 1 ? nwarmup : Int(nsettle))
        before = learning.retuned[]
        results = (w == 1 && resume !== nothing) ? resume.results :
            _sample_chains(nchains, parallel, seed + 101 * w, centre, metric,
            w == 1 ? nwarmup : Int(nsettle), Int(nw), Int(maxdepth), Float64(target_accept),
            Float64(maxdelta), Float64(init_scale), adapt_metric, adapt,
            c -> CTSEMJointMapDensity(mw, counter); init_eps=eps, starts=starts,
            chain_callback=watch, retune=retune)
        learning.offset += learning.nwarmup + Int(nw)
        moved = _jmap_moved(results, npar)
        if w == 1 && first_only
            _progress_break(learning.reporter)
            return _jmap_start(m, extra, results, moved, counter[] + extra, time() - t0,
                Float64(centre_bar), hessian !== nothing)
        end
        _progress_break(learning.reporter)
        Wl = _pool_draws(results[moved], D)
        scoring = _JMapStageReport("scoring", "refitting the map", verbose,
            progress_overwrite, progress_sink, progress_callback, learning.logp)
        at = time()
        m2, refit, lambda, rows = _jmap_refit(m, Wl, rng; maxrows=Int(maxrows),
            progress=scoring)
        _progress_done(scoring.reporter, @sprintf("%d rows", rows))
        refit_secs = time() - at
        map_secs = 0.0
        if learn
            pooled = hcat(pooled, _jmap_thetas(m, Wl[1:npar, :]))
            mapping = _JMapStageReport("mapping", "learning theta's map", verbose,
                progress_overwrite, progress_sink, progress_callback, learning.logp)
            at = time()
            m2 = _jmap_learn_tri(m2, pooled; parallel=true, progress=mapping)
            _progress_done(mapping.reporter, @sprintf("%d coordinates", npar))
            map_secs = time() - at
        end
        extra += rows
        starts = _jmap_restarts(m, m2, results, moved, Wl, rng)
        eps = median([r.stepsize for r in results[moved]])
        retuned = learning.retuned[] - before
        spread = _jmap_chain_spread(results)
        push!(windowlog, (window=Int(nw), stopped=count(!, moved), refit=refit,
            shrink=lambda, stepsize=eps, retuned=retuned, refit_secs=refit_secs,
            map_secs=map_secs, spread...))
        verbose && println(_console(), "  map window ", w, " of ", length(windows), ": ",
            nw, " draws a chain, ", count(!, moved), " chain(s) restarted, ",
            retuned > 0 ? string(retuned, " step size(s) retuned, ") : "", refit,
            " summaries refitted, step ", round(eps; sigdigits=3), " (chains ",
            round(spread.step_min; sigdigits=3), " to ", round(spread.step_max; sigdigits=3),
            ", tree depth up to ", round(spread.depth_max; digits=1), " on average)")
        m = m2
    end
    if checking
        mc = m
        _jmap_learn_window!(learning, Int(nsettle))
        before = learning.retuned[]
        results = _sample_chains(nchains, parallel, seed + 997, centre, metric, Int(nsettle),
            Int(ncheck), Int(maxdepth), Float64(target_accept), Float64(maxdelta), 0.0,
            adapt_metric, adapt, c -> CTSEMJointMapDensity(mc, counter); init_eps=eps,
            starts=starts, chain_callback=watch, retune=retune)
        moved = _jmap_moved(results, npar)
        Wl = _pool_draws(results[moved], D)
        starts = _jmap_restarts(m, m, results, moved, Wl, rng)
        eps = median([r.stepsize for r in results[moved]])
        push!(windowlog, (window=0, stopped=count(!, moved), refit=0, shrink=NaN,
            stepsize=eps, retuned=learning.retuned[] - before, refit_secs=0.0,
            map_secs=0.0, _jmap_chain_spread(results)...))
    end
    _progress_done(learning.reporter, @sprintf("%d iterations a chain", learning.total))
    learned = time() - t0 + (resume === nothing ? 0.0 : resume.secs)
    learn_gradients = counter[] + extra

    # The kept draws, under the final map, to the targets, judged on theta.
    mf = m
    final = _jmap_final(mf, counter, nchains, parallel, seed, centre, metric,
        isempty(windows) ? nwarmup : Int(nsettle), ndraws, Int(maxdepth),
        Float64(target_accept), Float64(maxdelta), Float64(init_scale), adapt_metric,
        adapt, eps, starts, Float64(min_ess), Float64(mean_ess),
        max(Int(max_draws), ndraws), Float64(rhat_target), verbose, progress_overwrite,
        progress_callback, progress_sink, t0, coordinate)
    results = final.results
    ndraws = final.ndraws

    # Theta and the natural deviations at every kept draw, in the report
    # objective's layout.
    reported = _sample_report_effects(laplace, report, sampler, results, nchains,
        ndraws, npar, save_effects, seed, verbose,
        (c, t) -> _jmap_point(mf, view(results[c].draws, :, t)))
    pooled = reported.pooled
    hybrid = reported.hybrid
    ndev = length(reported.effects.mean)
    total = nchains * ndraws
    # The kept draws' latent mean and sd per coordinate: 0 and 1 where the map
    # is right, so they say where it is not.
    lsum = zeros(D)
    lsq = zeros(D)
    for c in 1:nchains
        X = results[c].draws
        @inbounds for t in 1:ndraws, j in 1:D
            lsum[j] += X[j, t]
            lsq[j] += X[j, t]^2
        end
    end
    lmean = lsum ./ total
    lsd = sqrt.(max.(lsq ./ total .- lmean .^ 2, 0.0))
    return merge(_sample_result(pooled; nchains=nchains, ndraws=ndraws, npar=npar,
            ndim=npar + ndev, maxdepth=Int(maxdepth),
            stepsize=[r.stepsize for r in results],
            warmup_divergent=[r.warmup_divergent for r in results], trace=final.trace,
            saved_effects=save_effects, effects=reported.effects),
        (sampler="map", placement=_saem_placement_summary(nothing),
         map=(gradients=counter[] + extra + (hybrid ? total : 0),
              learn_gradients=learn_gradients, learn_secs=learned,
              # Iterations a chain before its kept draws: the learning
              # schedule and the final phase's warm-up.
              warmup=learning.total + (isempty(windows) ? nwarmup : Int(nsettle)),
              windows=windowlog, coordinates=D - npar, latent_mean=lmean,
              latent_sd=lsd, integrated=hybrid,
              # The most free parameters one member's summary is coupled to:
              # those every subject reads, and its own sampled TI values.
              coupled=length(mf.R) + maximum(length, mf.own; init=0),
              theta_parents=mf.tri === nothing ? -1 :
                maximum(length, mf.tri.parents; init=0))))
end

"""
    ctsem_sample_map(objective::CTSEMObjective, values; kwargs...)

The map on a filter objective: a model without random effects, or one whose
effects the filter carries as states (`intoverpop = 'augmented'`). Its
Laplace view has no effects (`ctsem_laplace_objective` with none: each subject
a unit of none), so the joint posterior is the marginal one and the map is
theta's alone, through the same code as every other target.
"""
ctsem_sample_map(objective::CTSEMObjective, values::AbstractVector; kwargs...) =
    ctsem_sample_map(ctsem_laplace_objective(objective), values; kwargs...)

"""
    ctsem_sample_map_check(laplace, values; hessian, cross, seed, nu, h)

The map's exactness at a random latent point, with the map built as
`ctsem_sample_map` builds it from `values`: the latent density's gradient
against central differences (in every population coordinate and `nu` effect
coordinates), and the map's log Jacobian against the log absolute determinant
of its Jacobian by differences, each as the largest error over the largest
entry. Both near zero say the change of variables is exact; for tests.
"""
function ctsem_sample_map_check(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    hessian::Union{Nothing,AbstractMatrix}=nothing, cross::Bool=true, seed::Integer=1,
    nu::Integer=6, h::Real=1e-5, tri::Bool=false)
    sampler = ctsem_sampler(laplace, length(values))
    npar = sampler.npar
    theta0 = collect(Float64, values)[1:npar]
    prec = _ctsem_prior_precision(laplace.objective, npar)
    cov = if hessian === nothing
        Matrix{Float64}(I, npar, npar)
    else
        H = Matrix{Float64}(hessian)[1:npar, 1:npar]
        _prior_bounded_covariance(Symmetric(-(H .+ transpose(H)) ./ 2), prec)
    end
    m, _ = _jmap_build(sampler, theta0, cov; cross=cross)
    D = sampler.ndim
    rng = Random.Xoshiro(UInt64(seed))
    if tri
        n = npar
        par = [[j for j in 1:(k - 1) if rand(rng) < 0.6] for k in 1:n]
        coef(sc) = [sc .* randn(rng, length(par[k])) for k in 1:n]
        pairs(sc) = [sc .* randn(rng, length(par[k]) * (length(par[k]) - 1) ÷ 2) for k in 1:n]
        m = _jmap_with(m; tri=_JMapTri(n, par, 0.1 .* randn(rng, n), coef(0.2), coef(0.1),
            pairs(0.1), 0.1 .* randn(rng, n), coef(0.2), 0.2 .* randn(rng, n),
            0.2 .* randn(rng, n)))
    end
    xt = 0.3 .* randn(rng, D)
    d = CTSEMJointMapDensity(m, Threads.Atomic{Int}(0))
    g = zeros(D)
    d(g, xt)
    idx = vcat(1:npar, (npar + 1):min(D, npar + Int(nu)))
    gp = zeros(D)
    gerr = 0.0
    for a in idx
        xp = copy(xt); xp[a] += h
        xm = copy(xt); xm[a] -= h
        gerr = max(gerr, abs((d(gp, xp) - d(gp, xm)) / (2h) - g[a]))
    end
    J = zeros(D, D)
    for a in 1:D
        xp = copy(xt); xp[a] += h
        xm = copy(xt); xm[a] -= h
        theta_p, x_p = _jmap_point(m, xp)
        theta_m, x_m = _jmap_point(m, xm)
        J[:, a] .= (x_p .- x_m) ./ (2h)
    end
    theta, ltheta = _jmap_theta(m, view(xt, 1:npar))
    _, logjac, _ = _jmap_effects(m, theta, xt)
    expected = sum(log, diag(m.Ct)) + ltheta + logjac
    numeric = logabsdet(J)[1]
    return (gradient=gerr / max(maximum(abs, g[idx]), 1e-300),
        jacobian=abs(numeric - expected) / max(abs(expected), 1.0), ndim=D,
        coordinates=D - npar)
end

export ctsem_sample_map_check
