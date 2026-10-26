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
summary to the free parameters outside `P` its likelihood also reads (`cross`).
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
Gaussian conditional from a few draws. A chain that stopped restarts from a
draw of those that moved; one check window under the final map does the same
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

"""Width of the seeded duals the map's dependence on theta is taken in: one
width for every model, so no model compiles a type of its own."""
const _JMAP_CHUNK = 8

"""Tag for the map's seeded duals."""
struct _JMapTag end

const _JMapDual = ForwardDiff.Dual{_JMapTag,Float64,_JMAP_CHUNK}

"""
    CTSEMJointMap

The map: theta's centre and whitening factor, and per member of each unit its
Gaussian summary in its raw effects (`Lam`, `h`) and coupling `G` to
`theta[R]`. `P` are the raw positions varying at some level, `row` a raw
position's row of `P`, `dep` the free positions the effects' map reads.
"""
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
# `f_k` quadratic in each earlier coordinate, `a_k` linear in them and capped by
# tanh, and a sinh-arcsinh marginal (Jones & Pewsey 2009; skew `eps`, tail
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
# mean gain grew. And
# `delta` is held at most 1, tails no lighter than the normal's: a small sample
# misses the tails, a fit from it makes them light, and the latent target is
# then heavy-tailed, where NUTS mixes slowly
# (`scripts/marginal-default-2026-10/tri_saturate.R`).

const _JMAP_TRI_CAP = 3.0
"""Beyond about this many whitened sd an earlier coordinate stops bending a
later one's location and scale."""
const _JMAP_TRI_SAT = 3.0
_jmap_tri_sat(u) = _JMAP_TRI_SAT * tanh(u / _JMAP_TRI_SAT)
"""Theta's triangular map is learnt for up to this many free parameters; a
component's fit has `3k + 1` coefficients and a Newton step over them."""
const _JMAP_TRI_MAXPAR = 30

struct _JMapTri
    n::Int
    f0::Vector{Float64}
    f1::Matrix{Float64}
    f2::Matrix{Float64}
    a0::Vector{Float64}
    a1::Matrix{Float64}
    eps::Vector{Float64}
    ld::Vector{Float64}
end

"""The whitened coordinates at latent `w` and the log Jacobian of that map."""
function _jmap_tri_forward(tri::_JMapTri, w::AbstractVector{T}) where {T<:Real}
    n = tri.n
    A = _JMAP_TRI_CAP
    u = Vector{T}(undef, n)
    lj = zero(T)
    @inbounds for k in 1:n
        f = T(tri.f0[k])
        a = T(tri.a0[k])
        for j in 1:(k - 1)
            s = _jmap_tri_sat(u[j])
            f += tri.f1[k, j] * s + tri.f2[k, j] * s^2
            a += tri.a1[k, j] * s
        end
        a = A * tanh(a / A)
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
        f = tri.f0[k]
        a = tri.a0[k]
        for j in 1:(k - 1)
            s = _jmap_tri_sat(u[j])
            f += tri.f1[k, j] * s + tri.f2[k, j] * s^2
            a += tri.a1[k, j] * s
        end
        a = A * tanh(a / A)
        v = (u[k] - f) * exp(-a)
        w[k] = sinh(exp(tri.ld[k]) * asinh(v) - tri.eps[k])
    end
    return w
end

struct _JMapTriForward
    tri::_JMapTri
end
function (F::_JMapTriForward)(w::AbstractVector{T}) where {T}
    u, lj = _jmap_tri_forward(F.tri, w)
    return vcat(u, lj)
end

"""The latent gradient from the whitened one, `gu`, the log Jacobian's own
gradient included."""
function _jmap_tri_pullback(tri::_JMapTri, w::AbstractVector, gu::AbstractVector)
    n = tri.n
    J = ForwardDiff.jacobian(_JMapTriForward(tri), collect(Float64, w))
    return transpose(view(J, 1:n, :)) * gu .+ view(J, n + 1, :)
end

# A component's negative log likelihood over whitened draws `U` (a column each),
# coefficients `p` = f0, f1[1:k-1], f2[1:k-1], a0, a1[1:k-1], eps, log delta,
# with a unit ridge on all but the level and the scale. `pinned`: `p` stops at
# eps and log delta is 0.
struct _JMapTriNLL{M}
    U::M
    k::Int
    pinned::Bool
end
function (F::_JMapTriNLL)(p::AbstractVector{T}) where {T}
    U = F.U
    k = F.k
    A = _JMAP_TRI_CAP
    delta = F.pinned ? one(T) : exp(p[3k + 1])
    epsk = p[3k]
    s = zero(T)
    @inbounds for i in axes(U, 2)
        f = p[1]
        a = p[2k]
        for j in 1:(k - 1)
            sj = _jmap_tri_sat(U[j, i])
            f += p[1 + j] * sj + p[k + j] * sj^2
            a += p[2k + j] * sj
        end
        a = A * tanh(a / A)
        v = (U[k, i] - f) * exp(-a)
        r = delta * asinh(v) - epsk
        s += sinh(r)^2 / 2 - (log(delta) + log(cosh(r)) - log1p(v^2) / 2) + a
    end
    reg = zero(T)
    for j in 2:length(p)
        j == 2k && continue
        reg += p[j]^2
    end
    return s + reg / 2
end

# Damped Newton on `F` from `p`, on the Hessian's eigenvalues floored.
function _jmap_tri_newton(F::_JMapTriNLL, p::Vector{Float64}, iters::Int)
    val = F(p)
    for _ in 1:iters
        g = ForwardDiff.gradient(F, p)
        H = ForwardDiff.hessian(F, p)
        ev = eigen(Symmetric((H .+ transpose(H)) ./ 2))
        lam = max.(ev.values, 1e-6 * maximum(abs, ev.values) + 1e-8)
        step = -(ev.vectors * ((transpose(ev.vectors) * g) ./ lam))
        t = 1.0
        improved = false
        for _ in 1:30
            pn = p .+ t .* step
            vn = F(pn)
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

"""
    _jmap_tri_fit(U; iters)

The triangular map fitted to whitened draws `U` (a column each) by maximum
likelihood, a component at a time, by damped Newton from the identity; a
component whose tails come out lighter than the normal's is fitted again with
its tail weight held at the normal's.
"""
function _jmap_tri_fit(U::AbstractMatrix; iters::Int=60)
    n = size(U, 1)
    f0 = zeros(n); f1 = zeros(n, n); f2 = zeros(n, n)
    a0 = zeros(n); a1 = zeros(n, n); epsv = zeros(n); ld = zeros(n)
    X = Matrix{Float64}(U)
    for k in 1:n
        p = _jmap_tri_newton(_JMapTriNLL(X, k, false), zeros(3k + 1), iters)
        if p[3k + 1] > 0
            p = vcat(_jmap_tri_newton(_JMapTriNLL(X, k, true), p[1:(3k)], iters), 0.0)
        end
        f0[k] = p[1]
        for j in 1:(k - 1)
            f1[k, j] = p[1 + j]
            f2[k, j] = p[k + j]
            a1[k, j] = p[2k + j]
        end
        a0[k] = p[2k]
        epsv[k] = p[3k]
        ld[k] = p[3k + 1]
    end
    return _JMapTri(n, f0, f1, f2, a0, a1, epsv, ld)
end

struct CTSEMJointMap{S}
    sampler::S
    npar::Int
    P::Vector{Int}
    row::Vector{Int}
    dep::Vector{Int}
    R::Vector{Int}
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
    CTSEMJointMap(m.sampler, m.npar, m.P, m.row, m.dep, m.R, m.members, Lam, h, G,
        m.thref, theta0, Ct, tri)

"""
    _jmap_learn_tri(m, pooled)

The map with theta's triangular map fitted from `pooled`, theta's draws over
every learning window so far, whitened by `m`'s own centre and factor; `m` as
it is while there are fewer than ten draws a coefficient of the largest
component. At most 4000 draws are read, spread through the pool.
"""
function _jmap_learn_tri(m::CTSEMJointMap, pooled::AbstractMatrix)
    n = m.npar
    N = size(pooled, 2)
    N >= 10 * (3n + 1) || return m
    pick = N > 4000 ? round.(Int, range(1, N; length=4000)) : collect(1:N)
    U = LowerTriangular(m.Ct) \ (pooled[:, pick] .- m.theta0)
    return _jmap_with(m; tri=_jmap_tri_fit(U))
end

"""Theta at latent `w`, and the log Jacobian of theta's triangular map (zero
without one; the whitening's, `log det Ct`, is a constant and left out)."""
function _jmap_theta(m::CTSEMJointMap, w::AbstractVector{<:Real})
    m.tri === nothing && return m.theta0 .+ m.Ct * w, 0.0
    u, lj = _jmap_tri_forward(m.tri, w)
    return m.theta0 .+ m.Ct * u, lj
end

"""Theta at each column of latent `W`."""
function _jmap_thetas(m::CTSEMJointMap, W::AbstractMatrix)
    m.tri === nothing && return m.theta0 .+ m.Ct * W
    out = zeros(m.npar, size(W, 2))
    for t in axes(W, 2)
        out[:, t] .= _jmap_theta(m, view(W, :, t))[1]
    end
    return out
end

"""Theta's latent coordinates."""
function _jmap_theta_inverse(m::CTSEMJointMap, theta::AbstractVector)
    u = LowerTriangular(m.Ct) \ (theta .- m.theta0)
    m.tri === nothing && return u
    return _jmap_tri_inverse(m.tri, u)
end

"""The latent gradient in theta's coordinates from the joint's gradient in
theta, `gtheta`, through the whitening and the triangular map, its log
Jacobian's gradient included."""
function _jmap_theta_pullback(m::CTSEMJointMap, w::AbstractVector, gtheta::AbstractVector)
    gu = transpose(m.Ct) * gtheta
    m.tri === nothing && return gu
    return _jmap_tri_pullback(m.tri, w, gu)
end


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
    positions = reduce(vcat, sampler.positions; init=Int[])
    R = cross ? setdiff(1:sampler.npar, P) : Int[]
    # The duals take the map's dependence on the effects' positions and the
    # level scales; the coupling's on theta[R] is linear and taken in closed
    # form (`_jmap_coupling_gradient!`), so `dep` holds no more of a large
    # model's parameters than its random effects reach.
    dep = sort(unique(filter(i -> 1 <= i <= sampler.npar, vcat(P, positions))))
    members = [Int[] for _ in 1:sampler.nunits]
    for s in eachindex(sampler.flat_unit)
        sampler.udims[sampler.flat_unit[s]] > 0 && push!(members[sampler.flat_unit[s]], s)
    end
    return P, row, dep, R, members
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

# Unit U's conditional given theta for the members' summaries: Q_U factored in
# place (`_ctsem_cholesky`, upper) and the mean m_U. `Ls`, `mu = theta[P]` and
# `thR = theta[R]` at theta. Lock-free kernels throughout: it runs on every
# leapfrog step of chains that run at once.
function _jmap_unit(m::CTSEMJointMap, Ls, mu::AbstractVector{T}, thR::AbstractVector,
    U::Int) where {T}
    d = m.sampler.udims[U]
    nP = length(m.P)
    Q = Matrix{T}(I, d, d)
    c = zeros(T, d)
    J = zeros(T, nP, d)
    LJ = zeros(T, nP, d)
    lin = zeros(T, nP)
    coupled = !isempty(m.R)
    shift = coupled ? thR .- m.thref : thR
    for s in m.members[U]
        _jmap_jacobian!(J, m, Ls, s)
        Lam = m.Lam[s]
        _ctsem_mul!(LJ, Lam, J)
        _ctsem_mulTN!(Q, J, LJ, true, true)
        @inbounds for a in 1:nP
            acc = m.h[s][a]
            for b in 1:nP
                acc -= Lam[a, b] * mu[b]
            end
            if coupled
                G = m.G[s]
                for b in eachindex(shift)
                    acc -= G[a, b] * shift[b]
                end
            end
            lin[a] = acc
        end
        @inbounds for q in 1:d, a in 1:nP
            c[q] += J[a, q] * lin[a]
        end
    end
    F = _ctsem_cholesky(Q, d)
    F.ok || throw(DomainError(U, "a unit's conditional precision did not factor"))
    return F, F \ c
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
    _jmap_effects(m, theta, xt; thR)

Every unit's effects at theta for latent `xt`, the map's log Jacobian in the
units, and each unit's factor, for any number type (the seeded duals take the
map's dependence on theta through it). `thR` holds the coupling's theta[R] at
given values -- the duals' pass hands them in plain, since that dependence is
taken in closed form.
"""
function _jmap_effects(m::CTSEMJointMap, theta::AbstractVector{T}, xt::AbstractVector;
    thR::Union{Nothing,AbstractVector}=nothing) where {T}
    sampler = m.sampler
    laplace = sampler.laplace
    th = _laplace_extend(laplace, theta)
    Ls = _laplace_popchols(th, laplace.spec)
    mu = th[m.P]
    thR = thR === nothing ? th[m.R] : thR
    us = Vector{Vector{T}}(undef, sampler.nunits)
    Fs = Vector{CTSEMCholesky{T,Matrix{T}}}(undef, sampler.nunits)
    logjac = zero(T)
    for U in 1:sampler.nunits
        d = sampler.udims[U]
        if d == 0
            us[U] = T[]
            Fs[U] = CTSEMCholesky(zeros(T, 0, 0), 0, true)
            continue
        end
        F, mean = _jmap_unit(m, Ls, mu, thR, U)
        us[U] = mean .+ _jmap_usolve(F, view(xt, _sample_urange(sampler, U)))
        Fs[U] = F
        @inbounds for i in 1:d
            logjac -= log(F.U[i, i])
        end
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
    _jmap_theta_gradient!(gtheta, m, theta, xt, gus)
    isempty(m.R) || _jmap_coupling_gradient!(gtheta, m, theta, gus, Fs)
    fill!(g, 0.0)
    g[1:npar] .= _jmap_theta_pullback(m, view(xt, 1:npar), gtheta)
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) || (g[r] .= _jmap_utsolve(Fs[U], gus[U]))
    end
    return v + logjac + ltheta
end

# The effects' and the Jacobian's dependence on theta at fixed latent z,
# contracted with the joint's gradient in the effects, added into `gtheta`:
# `_JMAP_CHUNK` positions of `dep` a pass, seeded as duals.
function _jmap_theta_gradient!(gtheta::Vector{Float64}, m::CTSEMJointMap,
    theta::Vector{Float64}, xt::AbstractVector, gus::Vector{Vector{Float64}})
    dep = m.dep
    nd = length(dep)
    nd == 0 && return gtheta
    th = Vector{_JMapDual}(undef, length(theta))
    for start in 1:_JMAP_CHUNK:nd
        @inbounds for j in eachindex(theta)
            th[j] = _JMapDual(theta[j])
        end
        stop = min(start + _JMAP_CHUNK - 1, nd)
        @inbounds for a in start:stop
            th[dep[a]] = _JMapDual(theta[dep[a]],
                ForwardDiff.Partials(ntuple(b -> b == a - start + 1 ? 1.0 : 0.0, _JMAP_CHUNK)))
        end
        uu, lj, _ = _jmap_effects(m, th, xt; thR=theta[m.R])
        acc = lj
        for U in eachindex(uu)
            u = uu[U]
            g = gus[U]
            @inbounds for q in eachindex(u)
                acc += g[q] * u[q]
            end
        end
        p = ForwardDiff.partials(acc)
        @inbounds for a in start:stop
            gtheta[dep[a]] += p[a - start + 1]
        end
    end
    return gtheta
end

# The coupling's share of theta's gradient: each unit's mean moves with
# theta[R] by `-Q_U^-1 sum_s J_s' G_s`, and the Jacobian not at all, so the
# contraction with the effects' gradient `g_U` is `-sum_s G_s' J_s Q_U^-1 g_U`,
# a solve a unit.
function _jmap_coupling_gradient!(gtheta::Vector{Float64}, m::CTSEMJointMap,
    theta::Vector{Float64}, gus::Vector{Vector{Float64}}, Fs)
    sampler = m.sampler
    laplace = sampler.laplace
    th = _laplace_extend(laplace, theta)
    Ls = _laplace_popchols(th, laplace.spec)
    nP = length(m.P)
    nR = length(m.R)
    for U in 1:sampler.nunits
        d = sampler.udims[U]
        d == 0 && continue
        v = Fs[U] \ gus[U]
        J = zeros(nP, d)
        w = zeros(nP)
        for s in m.members[U]
            _jmap_jacobian!(J, m, Ls, s)
            @inbounds for a in 1:nP
                acc = 0.0
                for q in 1:d
                    acc += J[a, q] * v[q]
                end
                w[a] = acc
            end
            G = m.G[s]
            @inbounds for b in 1:nR
                acc = 0.0
                for a in 1:nP
                    acc += G[a, b] * w[a]
                end
                gtheta[m.R[b]] -= acc
            end
        end
    end
    return gtheta
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
    thR = th[m.R]
    for U in 1:sampler.nunits
        r = _sample_urange(sampler, U)
        isempty(r) && continue
        F, mean = _jmap_unit(m, Ls, mu, thR, U)
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

"""A member's training rows as sums: `b_s`, its score in `b_s` and in
`theta[R]`, and the cross products the regressions read."""
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

# Members `mine` at the effects `x` (theta extended `th`, roots `Ls`): one row
# each into `stats`, none where the subject is not finite.
function _jmap_scores_range!(stats, m::CTSEMJointMap, flat::Vector{Int}, th, Ls, x, mine)
    sampler = m.sampler
    laplace = sampler.laplace
    aws = _laplace_workspace!(laplace, Float64, length(th))
    ekf = _laplace_ekf_workspace!(laplace, Float64)
    shifted = similar(th)
    gsub = similar(th)
    for k in mine
        s = flat[k]
        U = sampler.flat_unit[s]
        mm = sampler.flat_member[s]
        _laplace_member_values!(shifted, th, laplace.spec, Ls,
            view(x, _sample_urange(sampler, U)), laplace.units.offsets[U][mm])
        v = _jmap_score!(gsub, laplace, laplace.units.members[U][mm], aws, shifted, ekf)
        b = shifted[m.P]
        sb = gsub[m.P]
        sr = gsub[m.R]
        (isfinite(v) && all(isfinite, sb) && all(isfinite, sr)) || continue
        st = stats[s]
        st.n += 1
        st.sb .+= b
        st.ss .+= sb
        st.st .+= sr
        st.bb .+= b * transpose(b)
        st.bs .+= b * transpose(sb)
        st.bt .+= b * transpose(sr)
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
# `modes[U]` by central differences of its score, its coupling to theta[R]
# from the same differences, and the linear term that puts the summary's
# maximum at the mode.
function _jmap_curvature_range!(Lam, h, G, sampler, P, R, flat, th, Ls, modes, mine)
    laplace = sampler.laplace
    aws = _laplace_workspace!(laplace, Float64, length(th))
    ekf = _laplace_ekf_workspace!(laplace, Float64)
    shifted = similar(th)
    work = similar(th)
    gplus = similar(th)
    gminus = similar(th)
    g0 = similar(th)
    nP = length(P)
    nR = length(R)
    for k in mine
        s = flat[k]
        U = sampler.flat_unit[s]
        mm = sampler.flat_member[s]
        i = laplace.units.members[U][mm]
        _laplace_member_values!(shifted, th, laplace.spec, Ls, modes[U],
            laplace.units.offsets[U][mm])
        isfinite(_jmap_score!(g0, laplace, i, aws, shifted, ekf)) || continue
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
                Gt[a, j] = -(gplus[R[a]] - gminus[R[a]]) / (2step)
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
the directions they do not; with coupling, its score in `theta[R]` on `b` the
same way, for `G`. `thbar` is `theta[R]`'s mean over the rows. `nothing` where
the rows say too little.
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
    P, row, dep, R, members = _jmap_layout(sampler, cross)
    nP = length(P)
    nR = length(R)
    th = _laplace_extend(laplace, theta0)
    Ls = _laplace_popchols(th, laplace.spec)
    modes = [zeros(d) for d in sampler.udims]
    for U in 1:sampler.nunits
        sampler.udims[U] == 0 && continue
        _laplace_solve_unit_mode!(laplace, U, th, Ls)
        modes[U] = copy(laplace.modes[U])
    end
    nflat = length(sampler.flat_unit)
    Lam = [zeros(nP, nP) for _ in 1:nflat]
    h = [zeros(nP) for _ in 1:nflat]
    G = [zeros(nP, nR) for _ in 1:nflat]
    flat = reduce(vcat, members; init=Int[])
    _laplace_partition((mine, _w) -> _jmap_curvature_range!(Lam, h, G, sampler, P, R, flat,
        th, Ls, modes, mine), laplace, length(flat))
    S = Symmetric((Matrix{Float64}(cov) .+ transpose(Matrix{Float64}(cov))) ./ 2)
    Ct = Matrix(cholesky(Symmetric(Matrix(S) + 1e-12 * max(tr(S), 1.0) / npar * I)).L)
    m = CTSEMJointMap(sampler, npar, P, row, dep, R, members, Lam, h, G, th[R],
        copy(theta0), Ct, nothing)
    # The differences' passes, and about five for the units' mode solves.
    return m, isempty(flat) ? 0 : 2 * nP + 6
end

"""
    _jmap_refit(m, Wl, rng; maxrows)

The map refitted from latent draws `Wl` (columns, under `m`), at most `maxrows`
of them: theta's centre and whitening, each member's summary. Returns the new
map, the members refitted, theta's shrinkage, and the rows' cost in full passes.
"""
function _jmap_refit(m::CTSEMJointMap, Wl::AbstractMatrix, rng; maxrows::Int=400)
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
    nflat = length(sampler.flat_unit)
    stats = [_JMapStats(nP, length(m.R)) for _ in 1:nflat]
    rows = 0
    if nP > 0 && !isempty(flat)
        for t in 1:nr
            _, x = _jmap_point(m, view(W, :, t))
            _jmap_scores!(stats, m, flat, x)
        end
        rows = nr
    end
    thbar = isempty(m.R) ? Float64[] : vec(mean(thetas[m.R, :]; dims=2))
    Lam = copy(m.Lam)
    h = copy(m.h)
    G = copy(m.G)
    refit = 0
    for s in flat
        upd = _jmap_member_update(m.Lam[s], m.G[s], stats[s], thbar, m.thref)
        upd === nothing && continue
        Lam[s], h[s], G[s] = upd
        refit += 1
    end
    return _jmap_with(m; Lam=Lam, h=h, G=G, theta0=theta0, Ct=Ct), refit, shrunk.lambda, rows
end

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
    ctsem_sample_map(laplace, values; kwargs...)

Sample the joint posterior through a learnt map (this file's header). `laplace`
is the hybrid objective where the model has effects the filter integrates, the
layer otherwise; the keywords are `ctsem_sample`'s, plus the map's: `windows`
(draws a chain in each learning window, after `nwarmup` warm-up iterations on
the map from the fit), `nsettle` (the re-tuning warm-up after each refit and
before the kept draws), `ncheck` (a check window under the final map),
`maxrows` (the most draws a refit reads), `cross` (each member's summary
coupled to the free parameters outside its effects) and `report` (the layer's
objective, whose layout the effects are reported in: a unit whose effects are
all integrated has no block on the hybrid objective).

Returns `ctsem_sample`'s record, with the draws of theta and, with
`save_effects`, of the natural deviations, plus the map's: `gradients` (every
density evaluation and score pass, in gradients) and the learning windows.
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
    report::Union{Nothing,CTSEMLaplaceObjective}=nothing, theta_map::Bool=true)

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
    cov = if hessian === nothing
        th = _laplace_extend(laplace, theta0)
        Ls = _laplace_popchols(th, laplace.spec)
        for U in 1:sampler.nunits
            _laplace_solve_unit_mode!(laplace, U, th, Ls)
        end
        _conditional_population_covariance(sampler, th, Ls, nothing)
    else
        H = Matrix{Float64}(hessian)[1:npar, 1:npar]
        _prior_bounded_covariance(Symmetric(-(H .+ transpose(H)) ./ 2), prec)
    end
    m, extra = _jmap_build(sampler, theta0, cov; cross=cross)
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
    verbose && println(_console(), "Sampling through a learnt map: ", nchains, " chain(s), ",
        D, " dimensions (", npar, " population + ", D - npar, " effects), ",
        width, " worker(s) shared")

    rng = Random.Xoshiro(UInt64(seed) + 0x6a6d6170)
    centre = zeros(D)
    starts = nothing
    eps = Float64(stepsize)
    windowlog = NamedTuple[]
    # Theta's triangular map is learnt from the windows' pooled draws, over a
    # third window, up to `_JMAP_TRI_MAXPAR` parameters.
    learn = theta_map && npar <= _JMAP_TRI_MAXPAR
    windows = windows === nothing ? (learn ? [40, 80, 160] : [40, 80]) : windows
    pooled = zeros(npar, 0)
    for (w, nw) in enumerate(windows)
        mw = m
        results = _sample_chains(nchains, parallel, seed + 101 * w, centre, metric,
            w == 1 ? nwarmup : Int(nsettle), Int(nw), Int(maxdepth), Float64(target_accept),
            Float64(maxdelta), Float64(init_scale), adapt_metric, adapt,
            c -> CTSEMJointMapDensity(mw, counter); init_eps=eps, starts=starts)
        moved = _jmap_moved(results, npar)
        Wl = _pool_draws(results[moved], D)
        m2, refit, lambda, rows = _jmap_refit(m, Wl, rng; maxrows=Int(maxrows))
        if learn
            pooled = hcat(pooled, _jmap_thetas(m, Wl[1:npar, :]))
            m2 = _jmap_learn_tri(m2, pooled)
        end
        extra += rows
        starts = _jmap_restarts(m, m2, results, moved, Wl, rng)
        eps = median([r.stepsize for r in results[moved]])
        push!(windowlog, (window=Int(nw), stopped=count(!, moved), refit=refit,
            shrink=lambda, stepsize=eps))
        verbose && println(_console(), "  map window ", w, " of ", length(windows), ": ",
            nw, " draws a chain, ", count(!, moved), " chain(s) restarted, ", refit,
            " summaries refitted, step ", round(eps; sigdigits=3))
        m = m2
    end
    if ncheck > 0 && !isempty(windows)
        mc = m
        results = _sample_chains(nchains, parallel, seed + 997, centre, metric, Int(nsettle),
            Int(ncheck), Int(maxdepth), Float64(target_accept), Float64(maxdelta), 0.0,
            adapt_metric, adapt, c -> CTSEMJointMapDensity(mc, counter); init_eps=eps,
            starts=starts)
        moved = _jmap_moved(results, npar)
        Wl = _pool_draws(results[moved], D)
        starts = _jmap_restarts(m, m, results, moved, Wl, rng)
        eps = median([r.stepsize for r in results[moved]])
        push!(windowlog, (window=0, stopped=count(!, moved), refit=0, shrink=NaN,
            stepsize=eps))
    end
    learned = time() - t0
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
    for c in 1:nchains
        r = results[c]
        for t in 1:ndraws
            column = (c - 1) * ndraws + t
            theta, x = _jmap_point(mf, view(r.draws, :, t))
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
    end
    failed > 0 && verbose && println(_console(), "  ", failed,
        " draw(s) had a subject whose integrated effects could not be evaluated")
    return merge(_sample_result(pooled; nchains=nchains, ndraws=ndraws, npar=npar,
            ndim=npar + ndev, maxdepth=Int(maxdepth),
            stepsize=[r.stepsize for r in results],
            warmup_divergent=[r.warmup_divergent for r in results], trace=final.trace,
            saved_effects=save_effects,
            effects=_effect_summary(effect_sum, effect_sq, effect_cross, total, blocks)),
        (sampler="map", placement=_saem_placement_summary(nothing),
         map=(gradients=counter[] + extra + (hybrid ? total : 0),
              learn_gradients=learn_gradients, learn_secs=learned,
              windows=windowlog, coordinates=D - npar,
              integrated=hybrid, coupled=length(mf.R))))
end

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
        lower(sc) = [j < k ? sc * randn(rng) : 0.0 for k in 1:n, j in 1:n]
        m = _jmap_with(m; tri=_JMapTri(n, 0.1 .* randn(rng, n), lower(0.2), lower(0.1),
            0.1 .* randn(rng, n), lower(0.2), 0.2 .* randn(rng, n), 0.2 .* randn(rng, n)))
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
