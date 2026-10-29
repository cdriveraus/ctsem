# The gated floor's exact gradient by a reverse pass: laplace.jl's
# `_laplace_gated_term`, differentiated by hand, in passes whose count does not
# depend on the number of parameters.
#
# Forward mode over the parameters (`_laplace_gated_unit_gradient`'s fallback)
# repeats the whole rule -- the inner mode, the curvature at the mode and at
# two quadrature nodes, the conditional Newton points -- once per chunk of
# parameters, in nested duals. On the SNSF pilot (346 parameters) that made a
# gated unit cost a few hundred ordinary ones.
#
# The rule touches the parameters only through the unit objective `g` at five
# points -- the mode `z0`, the nodes `ub_j` and the Newton points `zz_j` -- and
# there only through `g`, `grad g` and `M = -hess g`. So the reverse runs in
# two layers:
#
#   * the rule's own small dense algebra, reversed in Float64: the softmax over
#     the nodes, the smoothstep weight, the smallest eigenpair (first-order
#     perturbation, as the forward dual code takes it), the complement basis
#     (which moves only with `v1`), the clip (Daleckii-Krein), the Newton solve
#     and the logdets. That leaves adjoints `(alpha, beta, Gamma)` of
#     `alpha g + beta' grad g + tr(Gamma M)` at each point, and the adjoint of
#     each point's position;
#   * each point's term back to the parameters by seeded passes in `z`
#     (`_laplace_gated_vjp!`): the value pass for `alpha`, one first-order pass
#     along `beta`, one second-order pass per eigendirection of `Gamma` -- each
#     returning the derivative with respect to every parameter at once, the
#     same forward-over-reverse the seeded gradient uses, with its seed types.
#
# The mode's own dependence closes it: `dz0 = M0^-1 d(grad g)/dtheta`, one more
# first-order pass along `M0^-1 z0bar`.
#
# The unit's coordinates are its own, or for an objective whose filter
# integrates deviations (`CTSEMHybrid`) the joint ones `(u, d_1, ...)`, where
# the gated rule runs (`_hybrid_joint_gated`); the two differ only in the three
# functions handed in.

"""
Whether `_laplace_gated_unit_gradient` takes the reverse pass, or forward mode
over the parameters (the oracle the reverse pass is tested against).
"""
const _LAPLACE_GATED_ADJOINT = Ref(true)

"""
    ctsem_set_gated_adjoint!(on)

Choose the gated floor's gradient route: the reverse pass (`true`, the default)
or forward mode over the parameters. Returns the previous setting.
"""
function ctsem_set_gated_adjoint!(on::Bool)
    previous = _LAPLACE_GATED_ADJOINT[]
    _LAPLACE_GATED_ADJOINT[] = on
    return previous
end
export ctsem_set_gated_adjoint!

"""
    _laplace_unit_theta_gradient(laplace, U, theta, Ls, dL, positions, z, layout)

At fixed coordinates `z`, the unit objective's gradient with respect to the raw
vector -- through each member's filter, and for the level parameters through
`L z` (`_laplace_chol_chain!`) -- and with respect to `z`, as
`(theta, z)`, or `nothing` if a member cannot be evaluated. Generic in
`eltype(z)`, so a seeded `z` returns both gradients' directional derivatives.
`layout` is `nothing` for the unit's own coordinates, or the joint layout for
`z = (u, d_1, ...)`, the filter at the fixed deviations `d_m`.
"""
function _laplace_unit_theta_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, dL, positions,
    z::AbstractVector{S}, layout) where {S}
    spec = laplace.spec
    units = laplace.units
    members = units.members[U]
    npar = length(theta)
    D = length(z)
    d = layout === nothing ? D : layout.d
    aws = _laplace_workspace!(laplace, S, npar)
    grad = _laplace_scratch_vector!(laplace, S, npar, :gated_theta_grad)
    shift = _laplace_scratch_vector!(laplace, S, npar, :gated_theta_shift)
    u = _laplace_scratch_vector!(laplace, S, d, :gated_theta_u)
    @inbounds for a in 1:d
        u[a] = z[a]
    end
    Gtheta = zeros(S, npar)
    Gz = zeros(S, D)
    mode = layout === nothing ? _HYBRID_PROFILE : _HYBRID_FIXED
    for m in eachindex(members)
        offsets = units.offsets[U][m]
        shifted = _laplace_member_values!(shift, theta, spec, Ls, u, offsets)
        if layout !== nothing
            hb = aws.ekf_ws.hyb
            q = layout.q[m]
            length(hb.dfix) >= q || (hb.dfix = zeros(S, q))
            @inbounds for c in 1:q
                hb.dfix[c] = z[layout.offs[m] + c]
            end
        end
        loglik = _laplace_subject_value_gradient!(grad,
            laplace.objective.subject_objectives[members[m]], aws, shifted, nothing,
            mode)
        _laplace_finite(loglik) || return nothing
        @inbounds for t in 1:npar
            Gtheta[t] += grad[t]
        end
        _laplace_chol_chain!(Gtheta, grad, spec, dL, positions, u, offsets)
        _laplace_loading_chain!(Gz, grad, spec, Ls, offsets)
        if layout !== nothing
            hb = aws.ekf_ws.hyb
            q = layout.q[m]
            _hybrid_ndev(hb) == q || error(string("member ", m, " integrates ",
                _hybrid_ndev(hb), " deviations in joint coordinates laid out for ", q))
            _hybrid_fixed_gradient!(Gz, layout.offs[m], hb; add = true)
        end
    end
    @inbounds for a in 1:D
        Gz[a] -= z[a]
    end
    return (theta=Gtheta, z=Gz)
end

"""
The gated rule's unit, its point's parameters and layout, which its passes
hand the hand-coded sweeps (`_tan_laplace_theta_gradient!`,
`_tan_laplace_theta_jet!`). Called, `_laplace_unit_theta_gradient` at `zz`:
what the tests' oracle differentiates (test/dual_oracles.jl).
"""
struct _GatedThetaGradient{L,DL,P,Y}
    laplace::L
    U::Int
    theta::Vector{Float64}
    Ls::Vector{Matrix{Float64}}
    dL::DL
    positions::P
    layout::Y
end
(g::_GatedThetaGradient)(zz) = _laplace_unit_theta_gradient(g.laplace, g.U, g.theta,
    g.Ls, g.dL, g.positions, zz, g.layout)

"""
`G`'s gradients and their derivatives along `beta`, `(theta, z)` as
`G(z + eps beta)`'s partials, by the hand-coded sweeps. `nothing` if a member
cannot be evaluated.
"""
function _gated_directional(G::_GatedThetaGradient, z::Vector{Float64},
    beta::Vector{Float64})
    npar, D = length(G.theta), length(z)
    Gt = zeros(Float64, npar, 2)
    Gz = zeros(Float64, D, 2)
    _tan_laplace_theta_gradient!(Gt, Gz, G.laplace, G.U, G.theta, G.Ls, G.dL,
        G.positions, z, reshape(beta, D, 1), 1, G.layout) || return nothing
    return (theta=Gt[:, 2], z=Gz[:, 2])
end

"""
    _laplace_gated_vjp!(thetabar, G, z, alpha, beta, Gamma)

Add `d/dtheta [alpha g + beta' grad g + tr(Gamma M)]` at fixed `z` to
`thetabar` and return its `z`-gradient, or `nothing` if a pass fails. `G` is
`_laplace_unit_theta_gradient` at the unit's coordinates; `beta` and `Gamma`
may be `nothing`. `M = -hess g`, so `tr(Gamma M) = -sum_i sigma_i gamma_i' hess
g gamma_i` over `Gamma`'s eigenpairs: a sum of second directional derivatives,
along `sqrt(|sigma_i|) gamma_i`, which a jet (`CTSEMTraceJet`) carries
`_CTSEM_JET_WIDTH` directions at a time (at `_ctsem_jet_width` of a pass's
count), the positive and the negative `sigma_i` in passes of their own.
"""
function _laplace_gated_vjp!(thetabar::Vector{Float64}, G, z::Vector{Float64},
    alpha::Float64, beta, Gamma)
    D = length(z)
    npar = length(thetabar)
    zbar = zeros(Float64, D)
    if alpha != 0
        r = G(z)
        r === nothing && return nothing
        @inbounds for t in 1:npar
            thetabar[t] += alpha * r.theta[t]
        end
        @inbounds for a in 1:D
            zbar[a] += alpha * r.z[a]
        end
    end
    if beta !== nothing && any(!iszero, beta)
        r = _gated_directional(G, z, collect(Float64, beta))
        r === nothing && return nothing
        @inbounds for t in 1:npar
            thetabar[t] += r.theta[t]
        end
        @inbounds for a in 1:D
            zbar[a] += r.z[a]
        end
    end
    if Gamma !== nothing
        E = _ctsem_symeig(Gamma)
        for sgn in (1.0, -1.0)
            dirs = [i for i in 1:D if sign(E.values[i]) == sgn]
            for first in 1:_CTSEM_JET_WIDTH:length(dirs)
                chunk = dirs[first:min(first + _CTSEM_JET_WIDTH - 1, length(dirs))]
                ok = _ctsem_jet_width(length(chunk)) == _CTSEM_JET_NARROW ?
                    _laplace_gated_jet_pass!(thetabar, zbar, G, z, E, chunk, sgn,
                        Val(_CTSEM_JET_NARROW)) :
                    _laplace_gated_jet_pass!(thetabar, zbar, G, z, E, chunk, sgn,
                        Val(_CTSEM_JET_WIDTH))
                ok || return nothing
            end
        end
    end
    return zbar
end

# One jet pass of `_laplace_gated_vjp!`'s trace term along `chunk`'s
# eigendirections, at width `W`: `false` where the pass fails.
function _laplace_gated_jet_pass!(thetabar, zbar, G::_GatedThetaGradient, z, E, chunk,
    sgn::Float64, ::Val{W}) where {W}
    D = length(z)
    # By the hand-coded jet sweeps (`_tan_laplace_theta_jet!`).
    B = zeros(D, length(chunk))
    @inbounds for (j, c) in enumerate(chunk), a in 1:D
        B[a, j] = sqrt(abs(E.values[c])) * E.vectors[a, c]
    end
    Gts, Gzs = zeros(length(thetabar)), zeros(D)
    _tan_laplace_theta_jet!(Gts, Gzs, G.laplace, G.U, G.theta, G.Ls, G.dL,
        G.positions, z, B, length(chunk), G.layout, CTSEMTraceJet{W}) || return false
    @inbounds for t in eachindex(thetabar)
        thetabar[t] -= sgn * Gts[t]
    end
    @inbounds for a in 1:D
        zbar[a] -= sgn * Gzs[a]
    end
    return true
end

# Dense helpers for the rule's own algebra: loops, never BLAS (no lock).
function _gated_mul(A::AbstractMatrix, B::AbstractMatrix)
    C = zeros(Float64, size(A, 1), size(B, 2))
    @inbounds for j in axes(B, 2), k in axes(A, 2)
        b = B[k, j]
        iszero(b) && continue
        for i in axes(A, 1)
            C[i, j] += A[i, k] * b
        end
    end
    return C
end
_gated_t(A::AbstractMatrix) = Matrix{Float64}(transpose(A))
function _gated_mulvec(A::AbstractMatrix, x::AbstractVector)
    y = zeros(Float64, size(A, 1))
    @inbounds for k in axes(A, 2)
        xk = x[k]
        for i in axes(A, 1)
            y[i] += A[i, k] * xk
        end
    end
    return y
end
function _gated_tmulvec(A::AbstractMatrix, x::AbstractVector)
    y = zeros(Float64, size(A, 2))
    @inbounds for j in axes(A, 2)
        acc = 0.0
        for i in axes(A, 1)
            acc += A[i, j] * x[i]
        end
        y[j] = acc
    end
    return y
end
_gated_sym(A) = (A .+ transpose(A)) ./ 2

"""`Q f(mu) Q'` and the clip's pieces, `f = max(., 1)` as `_laplace_clip_first_order`."""
function _gated_clip(B::Matrix{Float64})
    E = _ctsem_symeig(B)
    Q, mu = E.vectors, E.values
    f = [max(x, 1.0) for x in mu]
    return (Q=Q, mu=mu, f=f)
end

"""Adjoint of `clip(B)` given `Cbar`, by Daleckii-Krein with the forward's
coincidence rule."""
function _gated_clip_adjoint(Q::Matrix{Float64}, mu::Vector{Float64},
    Cbar::Matrix{Float64})
    n = length(mu)
    Ct = _gated_mul(_gated_t(Q), _gated_mul(Cbar, Q))
    fp(x) = x > 1.0 ? 1.0 : 0.0
    @inbounds for j in 1:n, i in 1:n
        coef = (i == j || abs(mu[i] - mu[j]) <= 1e-12 * max(1.0, abs(mu[i]))) ?
            fp(mu[i]) : (max(mu[i], 1.0) - max(mu[j], 1.0)) / (mu[i] - mu[j])
        Ct[i, j] *= coef
    end
    return _gated_mul(Q, _gated_mul(Ct, _gated_t(Q)))
end

"""`Q diag(v) Q'`."""
function _gated_spectral(Q::Matrix{Float64}, v::Vector{Float64})
    n = length(v)
    S = zeros(Float64, size(Q, 1), size(Q, 1))
    @inbounds for k in 1:n
        vk = v[k]
        iszero(vk) && continue
        for j in axes(Q, 1), i in axes(Q, 1)
            S[i, j] += Q[i, k] * vk * Q[j, k]
        end
    end
    return S
end

"""
    _laplace_gated_reverse(laplace, z0, objf, curvf, G, dcheck, npar)

The gated term's gradient at the unit's mode `z0`, by the reverse pass at the
top of this file, or `nothing` where the forward rule would not take its
eigenvalue path (the caller then differentiates it the other way). `objf(z)`
is the unit objective's `(value, gradient)`, `curvf(z)` its curvature as a
dense matrix, `G` `_laplace_unit_theta_gradient` in the same coordinates, and
`dcheck` the dimension the eigenvalue cap compares, as the forward rule's.
"""
function _laplace_gated_reverse(laplace::CTSEMLaplaceObjective, z0::Vector{Float64},
    objf, curvf, G, dcheck::Integer, npar::Integer)
    D = length(z0)
    lo, hi = laplace.gate_lo, laplace.gate_hi
    dcheck > _LAPLACE_EIGEN_MAXDIM[] && return nothing
    r0 = objf(z0)
    g0 = Float64(r0.value)
    isfinite(g0) || return nothing
    M0 = Matrix{Float64}(curvf(z0))
    M0 = _gated_sym(M0)
    Mw = copy(M0)
    F0 = _ctsem_cholesky(Mw, D)
    indefinite = !issuccess(F0)
    Ld = indefinite ? NaN : logdet(F0)
    E0 = _ctsem_symeig(M0)
    lam, V = E0.values, E0.vectors
    if D > 1 && (lam[2] - lam[1]) <= 1e-3 * max(abs(lam[2]), 1.0)
        return nothing
    end
    lam1 = lam[1]
    v1 = V[:, 1]
    # The smoothstep weight and its derivative.
    w, wp = if lam1 <= lo
        (1.0, 0.0)
    elseif lam1 >= hi
        (0.0, 0.0)
    else
        s = (lam1 - lo) / (hi - lo)
        (1 - (3 * s^2 - 2 * s^3), -(6 * s - 6 * s^2) / (hi - lo))
    end
    lap = indefinite ? NaN : g0 - Ld / 2
    thetabar = zeros(Float64, npar)
    gbar0 = 0.0
    Mbar0 = zeros(Float64, D, D)
    zbar0 = zeros(Float64, D)
    M0inv = _gated_spectral(V, [1 / l for l in lam])
    if iszero(w) && !indefinite
        # The forward rule returns the plain term.
        gbar0 = 1.0
        Mbar0 .-= M0inv ./ 2
        return _laplace_gated_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
    end
    lam1c = max(lam1, 1.0)
    Vh = V[:, 2:D]
    xs, ws = _gauss_hermite(3)
    nn = length(xs)
    terms = zeros(Float64, nn)
    # What each node keeps for the reverse.
    nodes = Vector{Any}(undef, nn)
    for j in 1:nn
        x = xs[j]
        if abs(x) < 1e-12 || D == 1
            if D == 1
                a = sqrt(2) * x / sqrt(lam1c)
                ub = z0 .+ a .* v1
                if abs(x) < 1e-12
                    h = g0
                    nodes[j] = (kind=:centre1,)
                else
                    rb = objf(ub)
                    h = Float64(rb.value)
                    nodes[j] = (kind=:node1, a=a, ub=ub)
                end
            else
                B0 = _gated_mul(_gated_t(Vh), _gated_mul(M0, Vh))
                C0 = _gated_clip(_gated_sym(B0))
                h = g0 - sum(log, C0.f) / 2
                nodes[j] = (kind=:centre, C=C0)
            end
        else
            a = sqrt(2) * x / sqrt(lam1c)
            ub = z0 .+ a .* v1
            rb = objf(ub)
            rgrad = Vector{Float64}(rb.gradient)
            Mb = _gated_sym(Matrix{Float64}(curvf(ub)))
            B = _gated_sym(_gated_mul(_gated_t(Vh), _gated_mul(Mb, Vh)))
            C = _gated_clip(B)
            Bcinv = _gated_spectral(C.Q, [1 / f for f in C.f])
            y = _gated_tmulvec(Vh, rgrad)
            zeta = _gated_mulvec(Bcinv, y)
            zz = ub .+ _gated_mulvec(Vh, zeta)
            rz = objf(zz)
            isfinite(rz.value) || return nothing
            h = Float64(rz.value) - sum(log, C.f) / 2
            nodes[j] = (kind=:node, a=a, ub=ub, rgrad=rgrad, Mb=Mb, C=C, Bcinv=Bcinv,
                y=y, zeta=zeta, zz=zz)
        end
        terms[j] = h + log(ws[j]) + x^2
    end
    peak = maximum(terms)
    isfinite(peak) || return nothing
    ex = [exp(t - peak) for t in terms]
    total = sum(ex)
    soft = peak + log(total) - log(lam1c) / 2 - log(pi) / 2
    # Reverse.
    if indefinite
        lapbar, softbar, wbar = 0.0, 1.0, 0.0
    else
        lapbar, softbar, wbar = 1 - w, w, soft - lap
    end
    lam1bar = wbar * wp
    lam1 > 1.0 && (lam1bar -= softbar / (2 * lam1c))
    Vhbar = zeros(Float64, D, max(D - 1, 0))
    v1bar = zeros(Float64, D)
    for j in 1:nn
        hbar = softbar * ex[j] / total
        node = nodes[j]
        if node.kind === :centre1
            gbar0 += hbar
        elseif node.kind === :centre
            gbar0 += hbar
            C = node.C
            Bbar = _gated_spectral(C.Q, [(-hbar / 2) * (m > 1.0 ? 1 / m : 0.0) for m in C.mu])
            Mbar0 .+= _gated_mul(Vh, _gated_mul(Bbar, _gated_t(Vh)))
            Vhbar .+= 2 .* _gated_mul(M0, _gated_mul(Vh, Bbar))
        elseif node.kind === :node1
            ubbar = _laplace_gated_vjp!(thetabar, G, node.ub, hbar, nothing, nothing)
            ubbar === nothing && return nothing
            zbar0 .+= ubbar
            v1bar .+= node.a .* ubbar
            lam1 > 1.0 && (lam1bar += dot(v1, ubbar) * (-node.a / (2 * lam1c)))
        else
            zzbar = _laplace_gated_vjp!(thetabar, G, node.zz, hbar, nothing, nothing)
            zzbar === nothing && return nothing
            C = node.C
            Bcbar = .-(hbar / 2) .* node.Bcinv
            ubbar = copy(zzbar)
            zeta = node.zeta
            @inbounds for k in 1:(D - 1), i in 1:D
                Vhbar[i, k] += zzbar[i] * zeta[k]
            end
            zetabar = _gated_tmulvec(Vh, zzbar)
            ybar = _gated_mulvec(node.Bcinv, zetabar)
            @inbounds for jj in 1:(D - 1), ii in 1:(D - 1)
                Bcbar[ii, jj] -= (ybar[ii] * zeta[jj] + zeta[ii] * ybar[jj]) / 2
            end
            Bbar = _gated_sym(_gated_clip_adjoint(C.Q, C.mu, _gated_sym(Bcbar)))
            Mbbar = _gated_mul(Vh, _gated_mul(Bbar, _gated_t(Vh)))
            Vhbar .+= 2 .* _gated_mul(node.Mb, _gated_mul(Vh, Bbar))
            rbar = _gated_mulvec(Vh, ybar)
            @inbounds for k in 1:(D - 1), i in 1:D
                Vhbar[i, k] += node.rgrad[i] * ybar[k]
            end
            ub2 = _laplace_gated_vjp!(thetabar, G, node.ub, 0.0, rbar, _gated_sym(Mbbar))
            ub2 === nothing && return nothing
            ubbar .+= ub2
            zbar0 .+= ubbar
            v1bar .+= node.a .* ubbar
            lam1 > 1.0 && (lam1bar += dot(v1, ubbar) * (-node.a / (2 * lam1c)))
        end
    end
    # The complement basis moves with v1 alone: dVh = -v1 dv1' Vh.
    if D > 1
        t = _gated_tmulvec(Vhbar, v1)
        v1bar .-= _gated_mulvec(Vh, t)
    end
    # The smallest eigenpair, as the forward's first-order perturbation.
    @inbounds for jj in 1:D, ii in 1:D
        Mbar0[ii, jj] += lam1bar * v1[ii] * v1[jj]
    end
    for k in 2:D
        c = dot(v1bar, view(V, :, k)) / (lam1 - lam[k])
        iszero(c) && continue
        @inbounds for jj in 1:D, ii in 1:D
            Mbar0[ii, jj] += c * (V[ii, k] * v1[jj] + v1[ii] * V[jj, k]) / 2
        end
    end
    if !indefinite
        gbar0 += lapbar
        Mbar0 .-= (lapbar / 2) .* M0inv
    end
    return _laplace_gated_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
end

"""The mode's own terms, then its dependence: `dz0 = M0^-1 d(grad g)/dtheta`."""
function _laplace_gated_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
    zb = _laplace_gated_vjp!(thetabar, G, z0, gbar0, nothing, _gated_sym(Mbar0))
    zb === nothing && return nothing
    zbar0 .+= zb
    beta = _gated_mulvec(M0inv, zbar0)
    _laplace_gated_vjp!(thetabar, G, z0, 0.0, beta, nothing) === nothing && return nothing
    return thetabar
end

"""
    _laplace_gated_unit_adjoint(laplace, U, theta)

Unit `U`'s gated term's gradient by `_laplace_gated_reverse`, in the unit's own
coordinates or, for a filter-integrating objective, the joint ones; `nothing`
where the reverse pass does not apply.
"""
function _laplace_gated_unit_adjoint(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64})
    st = _laplace_gated_setup(laplace, U, theta, laplace.modes[U])
    st === nothing && return nothing
    return _laplace_gated_reverse(laplace, st.z0, st.objf, st.curvf, st.G, st.dcheck,
        length(theta))
end

"""
    _laplace_gated_setup(laplace, U, theta, uhat)

What the gated rule's reverse pass reads at unit `U`'s mode `uhat`: the
coordinates `z0` -- the unit's own, or for an objective whose filter integrates
deviations the joint ones, the deviations at their modes -- the unit objective
`objf` and its curvature `curvf` there, the passes' `G` (`_GatedThetaGradient`)
and the dimension `dcheck` the eigenvalue cap compares; `nothing` where the
unit has no coordinate or a member cannot be evaluated. Shared by the gated
term's gradient (`_laplace_gated_unit_adjoint`) and a mode's fade in its
unit's sum (`_laplace_mode_fade`).
"""
function _laplace_gated_setup(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64})
    spec = laplace.spec
    npar = length(theta)
    isempty(uhat) && return nothing
    members = laplace.units.members[U]
    Ls = _laplace_popchols(theta, spec)
    dL = _laplace_level_chol_derivatives(theta, spec)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    aws = _laplace_workspace!(laplace, Float64, npar)
    if _hybrid_on(laplace.objective.params.hybrid)
        modes = _HybridModeHook{Float64}(length(members))
        r = _laplace_unit_loglik_gradient(laplace, U, theta, Ls, uhat, aws,
            eachindex(members), nothing, _HYBRID_PROFILE, modes)
        isfinite(r.value) || return nothing
        layout = _hybrid_layout(length(uhat), modes.q)
        z0 = zeros(Float64, layout.n)
        z0[1:layout.d] .= uhat
        for m in eachindex(members), c in 1:layout.q[m]
            z0[layout.offs[m] + c] = modes.dhat[m][c]
        end
        objf = let layout = layout
            zz -> _hybrid_joint_objective_gradient(laplace, U, theta, Ls, zz, aws, layout)
        end
        curvf = let layout = layout
            zz -> _hybrid_joint_curvature(laplace, U, theta, Ls, zz, layout).diag[1]
        end
        dcheck = length(uhat) + length(members) *
            _hybrid_items(laplace.objective.params.hybrid)
    else
        layout = nothing
        z0 = copy(uhat)
        blocks = laplace.units.blocks[U]
        objf = zz -> _laplace_unit_objective_gradient(laplace, U, theta, Ls, zz, aws)
        curvf = zz -> _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta,
            Ls, zz), blocks, length(zz))
        dcheck = length(uhat)
    end
    G = _GatedThetaGradient(laplace, Int(U), theta, Ls, dL, positions, layout)
    return (z0=z0, objf=objf, curvf=curvf, G=G, dcheck=dcheck)
end

"""
    _laplace_mode_fade(laplace, U, theta, uhat; gradient=false)

The smallest eigenvalue `lam` of unit `U`'s curvature at its mode `uhat`, in
the coordinates the gated rule decides on (`_laplace_gated_setup`), and with
`gradient` its derivative with respect to the raw vector, the mode's own
movement included: `v' dM v` along its eigenvector, the gated rule's closing
pass seeded with `v v'` (`_laplace_gated_close!`). What fades a mode out of its
unit's sum at a fold (`_laplace_mode_sum`). `nothing` where the unit cannot be
set up or a pass fails.
"""
function _laplace_mode_fade(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64}; gradient::Bool=false)
    geo = _laplace_modes_geometry(laplace, U, theta, [uhat])
    geo === nothing && return nothing
    E = geo.geos[1].E
    lam = E.values[1]
    gradient || return (lam=lam, grad=nothing)
    npar = length(theta)
    lam > 0 || return (lam=lam, grad=zeros(Float64, npar))
    v1 = E.vectors[:, 1]
    grad = _laplace_mode_close(laplace, U, theta, uhat, E, 0.0, v1 .* transpose(v1),
        zeros(Float64, length(v1)))
    grad === nothing && return nothing
    return (lam=lam, grad=grad)
end

"""
    _laplace_modes_geometry(laplace, U, theta, modes)

Unit `U`'s maxima `modes` as their combination reads them
(`_laplace_mode_sum`), in the gated rule's coordinates
(`_laplace_gated_setup`): for each, its point `z0`, the unit objective `g`
there and its curvature's eigendecomposition `E` (`geos`); for each pair, the
valley between them, the lowest `g` at the seven interior points
`_laplace_modes_separated` reads (`valley[k, j]`), and where on the segment it
lies, as the fraction of the way from `k` (`vt[k, j]`). Numbers only, since the
setup's closures hold the calling task's workspace. `nothing` where a mode's
quadratic cannot be had.
"""
function _laplace_modes_geometry(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, modes::AbstractVector)
    K = length(modes)
    geos = Vector{Any}(undef, K)
    objf = nothing
    for k in 1:K
        st = _laplace_gated_setup(laplace, U, theta, modes[k])
        st === nothing && return nothing
        k == 1 && (objf = st.objf)
        g = Float64(st.objf(st.z0).value)
        M0 = _gated_sym(Matrix{Float64}(st.curvf(st.z0)))
        (isfinite(g) && all(isfinite, M0)) || return nothing
        geos[k] = (z0=copy(st.z0), g=g, E=_ctsem_symeig(M0))
    end
    valley = fill(Inf, K, K)
    vt = zeros(Float64, K, K)
    for k in 1:K, j in (k + 1):K
        zk, zj = geos[k].z0, geos[j].z0
        lo, at = Inf, 0.5
        for i in 1:7
            t = i / 8
            v = Float64(objf(zk .+ t .* (zj .- zk)).value)
            # Not finite on the way: separate, as `_laplace_modes_separated` reads it.
            isfinite(v) || (v = -Inf)
            if v < lo
                lo, at = v, t
            end
        end
        valley[k, j] = valley[j, k] = lo
        vt[k, j] = at
        vt[j, k] = 1 - at
    end
    return (geos=geos, valley=valley, vt=vt)
end

"""
    _laplace_mode_point_gradient(laplace, U, theta, uhat, z)

Unit `U`'s objective's gradient with respect to the raw vector and to the
coordinates at the point `z` of the gated rule's coordinates, set up at its
maximum `uhat` (`_laplace_gated_setup`): `(theta, z)`, or `nothing` where it
cannot be had. What a valley between two maxima moves by (`_laplace_mode_sum`).
"""
function _laplace_mode_point_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64}, z::Vector{Float64})
    st = _laplace_gated_setup(laplace, U, theta, uhat)
    st === nothing && return nothing
    r = st.G(z)
    (r === nothing || !all(isfinite, r.theta) || !all(isfinite, r.z)) && return nothing
    return r
end

"""
    _laplace_mode_close(laplace, U, theta, uhat, E, gbar, Mbar, zbar)

The derivative with respect to the raw vector of `gbar g + tr(Mbar M) +
zbar' z0` at unit `U`'s mode `uhat`, in the gated rule's coordinates
(`_laplace_gated_setup`), the mode's own movement included: the gated rule's
closing pass (`_laplace_gated_close!`), `E` the curvature's eigendecomposition
there (`_laplace_modes_geometry`). `nothing` where the curvature is not
positive definite, a pass fails or the result is not finite.
"""
function _laplace_mode_close(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64}, E, gbar::Float64,
    Mbar::AbstractMatrix{Float64}, zbar::AbstractVector{Float64})
    all(>(0), E.values) || return nothing
    st = _laplace_gated_setup(laplace, U, theta, uhat)
    st === nothing && return nothing
    M0inv = _gated_spectral(E.vectors, [1 / l for l in E.values])
    grad = _laplace_gated_close!(zeros(Float64, length(theta)), st.G, st.z0, gbar,
        Matrix{Float64}(Mbar), collect(Float64, zbar), M0inv)
    (grad === nothing || !all(isfinite, grad)) && return nothing
    return grad
end
