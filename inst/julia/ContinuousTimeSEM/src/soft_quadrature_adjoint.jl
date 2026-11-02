# The soft-direction quadrature's exact gradient by a reverse pass: laplace.jl's
# `_laplace_soft_quadrature_term`, differentiated by hand, in passes whose count does not
# depend on the number of parameters.
#
# Forward mode over the parameters (`_laplace_soft_quadrature_unit_gradient`'s fallback)
# repeats the whole rule -- the inner mode, the curvature at the mode and at
# two quadrature nodes, the conditional Newton points -- once per chunk of
# parameters, in nested duals. On the SNSF pilot (346 parameters) that made a
# soft-quadrature unit cost a few hundred ordinary ones.
#
# The rule touches the parameters only through the unit objective `g` at five
# points -- the mode `z0`, the nodes `ub_j` and the Newton points `zz_j` -- and
# there only through `g`, `grad g` and `M = -hess g`. So the reverse runs in
# two layers:
#
#   * the rule's own small dense algebra, reversed in Float64: the softmax over
#     the nodes, the smoothstep weight, the smallest eigenpair (first-order
#     perturbation, as the forward dual code takes it), the complement basis
#     (which moves only with `v1`), each node's curvature floor
#     (Daleckii-Krein), the Newton solve
#     and the logdets. That leaves adjoints `(alpha, beta, Gamma)` of
#     `alpha g + beta' grad g + tr(Gamma M)` at each point, and the adjoint of
#     each point's position;
#   * each point's term back to the parameters by seeded passes in `z`
#     (`_laplace_soft_quadrature_vjp!`): the value pass for `alpha`, one first-order pass
#     along `beta`, one second-order pass per eigendirection of `Gamma` -- each
#     returning the derivative with respect to every parameter at once, the
#     same forward-over-reverse the seeded gradient uses, with its seed types.
#
# The mode's own dependence closes it: `dz0 = M0^-1 d(grad g)/dtheta`, one more
# first-order pass along `M0^-1 z0bar`.
#
# The unit's coordinates are its own, or for an objective whose filter
# integrates deviations (`CTSEMHybrid`) the joint ones `(u, d_1, ...)`, where
# the soft-direction quadrature runs (`_hybrid_joint_soft_quadrature`); the two differ only in the three
# functions handed in.

"""
Whether `_laplace_soft_quadrature_unit_gradient` takes the reverse pass, or forward mode
over the parameters (the oracle the reverse pass is tested against).
"""
const _LAPLACE_SOFT_QUADRATURE_ADJOINT = Ref(true)

"""
    ctsem_set_soft_quadrature_adjoint!(on)

Choose the soft-direction quadrature's gradient route: the reverse pass (`true`, the default)
or forward mode over the parameters. Returns the previous setting.
"""
function ctsem_set_soft_quadrature_adjoint!(on::Bool)
    previous = _LAPLACE_SOFT_QUADRATURE_ADJOINT[]
    _LAPLACE_SOFT_QUADRATURE_ADJOINT[] = on
    return previous
end
export ctsem_set_soft_quadrature_adjoint!

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
    grad = _laplace_scratch_vector!(laplace, S, npar, :soft_theta_grad)
    shift = _laplace_scratch_vector!(laplace, S, npar, :soft_theta_shift)
    u = _laplace_scratch_vector!(laplace, S, d, :soft_theta_u)
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
The soft-direction quadrature's unit, its point's parameters and layout, which its passes
hand the hand-coded sweeps (`_tan_laplace_theta_gradient!`,
`_tan_laplace_theta_jet!`), and the coordinates each member reads
(`_soft_quadrature_member_coordinates`). Called, `_laplace_unit_theta_gradient`
at `zz`: what the tests' oracle differentiates (test/dual_oracles.jl).
"""
struct _SoftQuadratureThetaGradient{L,DL,P,Y}
    laplace::L
    U::Int
    theta::Vector{Float64}
    Ls::Vector{Matrix{Float64}}
    dL::DL
    positions::P
    layout::Y
    coords::Vector{Vector{Int}}
end
_SoftQuadratureThetaGradient(laplace, U, theta, Ls, dL, positions, layout) =
    _SoftQuadratureThetaGradient(laplace, Int(U), theta, Ls, dL, positions, layout,
        _soft_quadrature_member_coordinates(laplace, U, layout))
(g::_SoftQuadratureThetaGradient)(zz) = _laplace_unit_theta_gradient(g.laplace, g.U, g.theta,
    g.Ls, g.dL, g.positions, zz, g.layout)

"""
    _soft_quadrature_member_coordinates(laplace, U, layout)

The coordinates of the unit's `z` each member's filter reads: those of every
block holding it -- its own and its ancestors' -- and in joint coordinates
(`layout`) its own fixed deviations. A member's curvature is zero outside
them, which is what lets a trace term be taken member by member
(`_soft_quadrature_trace!`).
"""
function _soft_quadrature_member_coordinates(laplace::CTSEMLaplaceObjective, U::Integer,
    layout)
    nm = length(laplace.units.members[U])
    coords = [Int[] for _ in 1:nm]
    for block in laplace.units.blocks[U], m in block.members
        append!(coords[m], (block.offset + 1):(block.offset + block.size))
    end
    if layout !== nothing
        for m in 1:nm
            append!(coords[m], layout.offs[m] .+ (1:layout.q[m]))
        end
    end
    return coords
end

"""
`G`'s gradients and their derivatives along `beta`, `(theta, z)` as
`G(z + eps beta)`'s partials, by the hand-coded sweeps. `nothing` if a member
cannot be evaluated.
"""
function _soft_quadrature_directional(G::_SoftQuadratureThetaGradient, z::Vector{Float64},
    beta::Vector{Float64})
    npar, D = length(G.theta), length(z)
    Gt = zeros(Float64, npar, 2)
    Gz = zeros(Float64, D, 2)
    _tan_laplace_theta_gradient!(Gt, Gz, G.laplace, G.U, G.theta, G.Ls, G.dL,
        G.positions, z, reshape(beta, D, 1), 1, G.layout) || return nothing
    return (theta=Gt[:, 2], z=Gz[:, 2])
end

"""
    _laplace_soft_quadrature_vjp!(thetabar, G, z, alpha, beta, Gamma)

Add `d/dtheta [alpha g + beta' grad g + tr(Gamma M)]` at fixed `z` to
`thetabar` and return its `z`-gradient, or `nothing` if a pass fails. `G` is
`_laplace_unit_theta_gradient` at the unit's coordinates; `beta` and `Gamma`
may be `nothing`. `Gamma` is a dense matrix or a `_SoftGamma`; its trace term
is `_soft_quadrature_trace!`'s.
"""
function _laplace_soft_quadrature_vjp!(thetabar::Vector{Float64}, G, z::Vector{Float64},
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
        r = _soft_quadrature_directional(G, z, collect(Float64, beta))
        r === nothing && return nothing
        @inbounds for t in 1:npar
            thetabar[t] += r.theta[t]
        end
        @inbounds for a in 1:D
            zbar[a] += r.z[a]
        end
    end
    if Gamma !== nothing
        _soft_quadrature_trace!(thetabar, zbar, G, z, Gamma) || return nothing
    end
    return zbar
end

"""
    _soft_quadrature_trace!(thetabar, zbar, G, z, Gamma)

Add the derivatives of `tr(Gamma M)` at fixed `z`, with respect to the raw
vector and to `z`, to `thetabar` and `zbar`; `false` if a pass fails. `M = I -
sum_m hess g_m`, member `m`'s Hessian zero outside the coordinates it reads
(`G.coords[m]`), so the trace is the members' `-tr(Gamma_m hess g_m)`,
`Gamma_m` `Gamma`'s block on those coordinates (`_soft_member_gamma`): over
`Gamma_m`'s eigenpairs a sum of second directional derivatives along
`sqrt(|sigma_i|) gamma_i`, which a jet (`CTSEMTraceJet`) carries
`_CTSEM_JET_WIDTH` at a time (at `_ctsem_jet_width` of the widest member's
count), each member along its own, the positive and the negative `sigma_i` in
passes of their own. A nested unit then costs two sweeps of its members where
`Gamma`'s own eigendirections, each reaching every member, cost about `d / 8`.
"""
function _soft_quadrature_trace!(thetabar::Vector{Float64}, zbar::Vector{Float64},
    G::_SoftQuadratureThetaGradient, z::Vector{Float64}, Gamma)
    coords = G.coords
    nm = length(coords)
    pos = Vector{Matrix{Float64}}(undef, nm)
    neg = Vector{Matrix{Float64}}(undef, nm)
    for m in 1:nm
        Gm = _soft_member_gamma(Gamma, G, m)
        n = size(Gm, 1)
        if n == 0
            pos[m] = neg[m] = zeros(Float64, 0, 0)
            continue
        end
        E = _ctsem_symeig(Gm)
        ip = findall(>(0), E.values)
        ineg = findall(<(0), E.values)
        pos[m] = [sqrt(E.values[c]) * E.vectors[i, c] for i in 1:n, c in ip]
        neg[m] = [sqrt(-E.values[c]) * E.vectors[i, c] for i in 1:n, c in ineg]
    end
    W = _CTSEM_JET_WIDTH
    Gts, Gzs = zeros(length(thetabar)), zeros(length(z))
    for (dirs, sgn) in ((pos, 1.0), (neg, -1.0))
        most = maximum(d -> size(d, 2), dirs; init=0)
        for first in 1:W:most
            nds = [clamp(size(dirs[m], 2) - first + 1, 0, W) for m in 1:nm]
            chunk = [dirs[m][:, first:(first + nds[m] - 1)] for m in 1:nm]
            ok = _ctsem_jet_width(maximum(nds)) == _CTSEM_JET_NARROW ?
                _tan_laplace_theta_jet!(Gts, Gzs, G.laplace, G.U, G.theta, G.Ls, G.dL,
                    G.positions, z, coords, chunk, nds, G.layout,
                    CTSEMTraceJet{_CTSEM_JET_NARROW}) :
                _tan_laplace_theta_jet!(Gts, Gzs, G.laplace, G.U, G.theta, G.Ls, G.dL,
                    G.positions, z, coords, chunk, nds, G.layout,
                    CTSEMTraceJet{_CTSEM_JET_WIDTH})
            ok || return false
            @inbounds for t in eachindex(thetabar)
                thetabar[t] -= sgn * Gts[t]
            end
            @inbounds for a in eachindex(zbar)
                zbar[a] -= sgn * Gzs[a]
            end
        end
    end
    return true
end

"""`Gamma`'s block on member `m`'s coordinates (`G.coords[m]`), symmetrised."""
function _soft_member_gamma(Gamma::AbstractMatrix, G::_SoftQuadratureThetaGradient,
    m::Integer)
    S = G.coords[m]
    return [(Gamma[S[i], S[j]] + Gamma[S[j], S[i]]) / 2 for i in eachindex(S),
        j in eachindex(S)]
end

"""
A trace term's adjoint as the soft-direction quadrature's factorized reverse
forms it (`_laplace_soft_quadrature_reverse_factored`): `c M^-1 + sum_k w_k
sym(a_k b_k')`, `M` a block curvature held as the entries of its inverse in its
own sparsity pattern (`_laplace_selected_inverse`) -- all a member's block of
it needs, since a member's blocks form a chain. Nothing `d x d` is formed.
"""
struct _SoftGamma
    c::Float64
    Cdiag::Vector{Matrix{Float64}}
    Ccoup::Vector{Vector{Matrix{Float64}}}
    blocks::Vector{CTSEMLaplaceBlock}
    # Each coordinate's block and its place there.
    block_of::Vector{Int}
    place::Vector{Int}
    w::Vector{Float64}
    a::Vector{Vector{Float64}}
    b::Vector{Vector{Float64}}
end

"""`_SoftGamma` from `c`, the factorization `F` of `M` (`_soft_shifted_factor`)
and its blocks, and the low-rank terms."""
function _SoftGamma(c::Float64, F, blocks::Vector{CTSEMLaplaceBlock},
    w::Vector{Float64}=Float64[], a::Vector{Vector{Float64}}=Vector{Float64}[],
    b::Vector{Vector{Float64}}=Vector{Float64}[])
    Cdiag, Ccoup = iszero(c) ?
        ([zeros(Float64, k.size, k.size) for k in blocks],
            [[zeros(Float64, k.size, blocks[s].size) for s in k.ancestors] for k in blocks]) :
        _laplace_selected_inverse(F.factors, F.coupling, blocks)
    D = sum(k.size for k in blocks; init=0)
    block_of, place = zeros(Int, D), zeros(Int, D)
    for (bi, k) in enumerate(blocks), i in 1:k.size
        block_of[k.offset + i] = bi
        place[k.offset + i] = i
    end
    return _SoftGamma(c, Cdiag, Ccoup, blocks, block_of, place, w, a, b)
end

function _soft_member_gamma(Gamma::_SoftGamma, G::_SoftQuadratureThetaGradient,
    m::Integer)
    S = G.coords[m]
    n = length(S)
    out = zeros(Float64, n, n)
    if !iszero(Gamma.c)
        for j in 1:n, i in 1:n
            out[i, j] = Gamma.c * _soft_inverse_entry(Gamma, S[i], S[j])
        end
    end
    for k in eachindex(Gamma.w)
        wk, ak, bk = Gamma.w[k], Gamma.a[k], Gamma.b[k]
        iszero(wk) && continue
        @inbounds for j in 1:n, i in 1:n
            out[i, j] += wk * (ak[S[i]] * bk[S[j]] + bk[S[i]] * ak[S[j]]) / 2
        end
    end
    return out
end

# `M^-1[p, q]` for two coordinates whose blocks are the same or one an
# ancestor of the other, as every pair of one member's are.
function _soft_inverse_entry(Gamma::_SoftGamma, p::Int, q::Int)
    bp, bq = Gamma.block_of[p], Gamma.block_of[q]
    ip, iq = Gamma.place[p], Gamma.place[q]
    bp == bq && return Gamma.Cdiag[bp][ip, iq]
    t = findfirst(==(bq), Gamma.blocks[bp].ancestors)
    t === nothing || return Gamma.Ccoup[bp][t][ip, iq]
    t = findfirst(==(bp), Gamma.blocks[bq].ancestors)
    t === nothing || return Gamma.Ccoup[bq][t][iq, ip]
    error("coordinates ", p, " and ", q, " are in blocks with no ancestor relation")
end

# Dense helpers for the rule's own algebra: loops, never BLAS (no lock).
function _dense_mul(A::AbstractMatrix, B::AbstractMatrix)
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
_dense_t(A::AbstractMatrix) = Matrix{Float64}(transpose(A))
function _dense_mulvec(A::AbstractMatrix, x::AbstractVector)
    y = zeros(Float64, size(A, 1))
    @inbounds for k in axes(A, 2)
        xk = x[k]
        for i in axes(A, 1)
            y[i] += A[i, k] * xk
        end
    end
    return y
end
function _dense_tmulvec(A::AbstractMatrix, x::AbstractVector)
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
_dense_sym(A) = (A .+ transpose(A)) ./ 2

"""A soft-direction node's curvature `B` with each eigenvalue floored, `Q f(mu)
Q'`, as its pieces: `f` the floor `_laplace_soft_node_curvature_floor`, as the
forward rule takes it (`_laplace_soft_node_floored_curvature`)."""
function _soft_node_floored_curvature_pieces(B::Matrix{Float64})
    E = _ctsem_symeig(B)
    Q, mu = E.vectors, E.values
    f = [_laplace_soft_node_curvature_floor(x) for x in mu]
    return (Q=Q, mu=mu, f=f)
end

"""The adjoint of a soft-direction node's floored curvature
(`_soft_node_floored_curvature_pieces`) given `Cbar`, by Daleckii-Krein with the
forward's coincidence rule."""
function _soft_node_floored_curvature_adjoint(Q::Matrix{Float64},
    mu::Vector{Float64}, Cbar::Matrix{Float64})
    n = length(mu)
    Ct = _dense_mul(_dense_t(Q), _dense_mul(Cbar, Q))
    f = _laplace_soft_node_curvature_floor
    fp = _laplace_soft_node_curvature_floor_slope
    @inbounds for j in 1:n, i in 1:n
        coef = (i == j || abs(mu[i] - mu[j]) <= 1e-12 * max(1.0, abs(mu[i]))) ?
            fp(mu[i]) : (f(mu[i]) - f(mu[j])) / (mu[i] - mu[j])
        Ct[i, j] *= coef
    end
    return _dense_mul(Q, _dense_mul(Ct, _dense_t(Q)))
end

"""`Q diag(v) Q'`."""
function _dense_spectral(Q::Matrix{Float64}, v::Vector{Float64})
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
    _laplace_soft_quadrature_reverse(laplace, z0, objf, curvf, G, npar)

The soft-quadrature term's gradient at the unit's mode `z0`, by the reverse pass at the
top of this file in dense algebra, or `nothing` where the forward rule would
not take its eigenvalue path (the caller then differentiates it the other
way). `objf(z)` is the unit objective's `(value, gradient)`, `curvf(z)` its
curvature as a dense matrix, `G` `_laplace_unit_theta_gradient` in the same
coordinates. Where the floor at the nodes is the identity
`_laplace_soft_quadrature_reverse_factored` takes the same gradient from block
factorizations, and this is its oracle.
"""
function _laplace_soft_quadrature_reverse(laplace::CTSEMLaplaceObjective, z0::Vector{Float64},
    objf, curvf, G, npar::Integer)
    D = length(z0)
    lo, hi = laplace.soft_blend_lo, laplace.soft_blend_hi
    r0 = objf(z0)
    g0 = Float64(r0.value)
    isfinite(g0) || return nothing
    M0 = Matrix{Float64}(curvf(z0))
    M0 = _dense_sym(M0)
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
    M0inv = _dense_spectral(V, [1 / l for l in lam])
    if iszero(w) && !indefinite
        # The forward rule returns the plain term.
        gbar0 = 1.0
        Mbar0 .-= M0inv ./ 2
        return _laplace_soft_quadrature_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
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
                B0 = _dense_mul(_dense_t(Vh), _dense_mul(M0, Vh))
                C0 = _soft_node_floored_curvature_pieces(_dense_sym(B0))
                h = g0 - sum(log, C0.f) / 2
                nodes[j] = (kind=:centre, C=C0)
            end
        else
            a = sqrt(2) * x / sqrt(lam1c)
            ub = z0 .+ a .* v1
            rb = objf(ub)
            rgrad = Vector{Float64}(rb.gradient)
            Mb = _dense_sym(Matrix{Float64}(curvf(ub)))
            B = _dense_sym(_dense_mul(_dense_t(Vh), _dense_mul(Mb, Vh)))
            C = _soft_node_floored_curvature_pieces(B)
            Bcinv = _dense_spectral(C.Q, [1 / f for f in C.f])
            y = _dense_tmulvec(Vh, rgrad)
            zeta = _dense_mulvec(Bcinv, y)
            zz = ub .+ _dense_mulvec(Vh, zeta)
            hb, hs = Float64(rb.value), Float64(objf(zz).value)
            h = _laplace_soft_node_stepped(hb, hs) - sum(log, C.f) / 2
            nodes[j] = (kind=:node, a=a, ub=ub, rgrad=rgrad, Mb=Mb, C=C, Bcinv=Bcinv,
                y=y, zeta=zeta, zz=zz, ps=_laplace_soft_node_stepped_slope(hb, hs))
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
            # -sum log f(mu) / 2 is a spectral function of B: its gradient is
            # Q diag(-f'(mu) / (2 f(mu))) Q', f the node floor.
            Bbar = _dense_spectral(C.Q, [(-hbar / 2) *
                _laplace_soft_node_curvature_floor_slope(C.mu[i]) / C.f[i]
                for i in eachindex(C.mu)])
            Mbar0 .+= _dense_mul(Vh, _dense_mul(Bbar, _dense_t(Vh)))
            Vhbar .+= 2 .* _dense_mul(M0, _dense_mul(Vh, Bbar))
        elseif node.kind === :node1
            ubbar = _laplace_soft_quadrature_vjp!(thetabar, G, node.ub, hbar, nothing, nothing)
            ubbar === nothing && return nothing
            zbar0 .+= ubbar
            v1bar .+= node.a .* ubbar
            lam1 > 1.0 && (lam1bar += dot(v1, ubbar) * (-node.a / (2 * lam1c)))
        else
            # The stepped point's share of the node's value and the node's own
            # point the rest (`_laplace_soft_node_stepped`): a step taken back
            # entirely is not read at all.
            hz = hbar * node.ps
            zzbar = iszero(hz) ? zeros(Float64, D) :
                _laplace_soft_quadrature_vjp!(thetabar, G, node.zz, hz, nothing, nothing)
            zzbar === nothing && return nothing
            C = node.C
            Bcbar = .-(hbar / 2) .* node.Bcinv
            ubbar = copy(zzbar)
            zeta = node.zeta
            @inbounds for k in 1:(D - 1), i in 1:D
                Vhbar[i, k] += zzbar[i] * zeta[k]
            end
            zetabar = _dense_tmulvec(Vh, zzbar)
            ybar = _dense_mulvec(node.Bcinv, zetabar)
            @inbounds for jj in 1:(D - 1), ii in 1:(D - 1)
                Bcbar[ii, jj] -= (ybar[ii] * zeta[jj] + zeta[ii] * ybar[jj]) / 2
            end
            Bbar = _dense_sym(_soft_node_floored_curvature_adjoint(C.Q, C.mu, _dense_sym(Bcbar)))
            Mbbar = _dense_mul(Vh, _dense_mul(Bbar, _dense_t(Vh)))
            Vhbar .+= 2 .* _dense_mul(node.Mb, _dense_mul(Vh, Bbar))
            rbar = _dense_mulvec(Vh, ybar)
            @inbounds for k in 1:(D - 1), i in 1:D
                Vhbar[i, k] += node.rgrad[i] * ybar[k]
            end
            ub2 = _laplace_soft_quadrature_vjp!(thetabar, G, node.ub, hbar - hz, rbar,
                _dense_sym(Mbbar))
            ub2 === nothing && return nothing
            ubbar .+= ub2
            zbar0 .+= ubbar
            v1bar .+= node.a .* ubbar
            lam1 > 1.0 && (lam1bar += dot(v1, ubbar) * (-node.a / (2 * lam1c)))
        end
    end
    # The complement basis moves with v1 alone: dVh = -v1 dv1' Vh.
    if D > 1
        t = _dense_tmulvec(Vhbar, v1)
        v1bar .-= _dense_mulvec(Vh, t)
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
    return _laplace_soft_quadrature_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
end

"""The mode's own terms, then its dependence: `dz0 = M0^-1 d(grad g)/dtheta`."""
function _laplace_soft_quadrature_close!(thetabar, G, z0, gbar0, Mbar0, zbar0, M0inv)
    zb = _laplace_soft_quadrature_vjp!(thetabar, G, z0, gbar0, nothing, _dense_sym(Mbar0))
    zb === nothing && return nothing
    zbar0 .+= zb
    beta = _dense_mulvec(M0inv, zbar0)
    _laplace_soft_quadrature_vjp!(thetabar, G, z0, 0.0, beta, nothing) === nothing && return nothing
    return thetabar
end

"""
    _laplace_soft_quadrature_reverse_factored(laplace, z0, objf, curvb, blocks, G, npar)

`_laplace_soft_quadrature_reverse` by block factorizations, where the forward
rule takes them (`_laplace_soft_quadrature_factored`): `nothing` wherever that
declines, and the caller then takes the dense reverse. `curvb(z)` is the
unit's curvature as a block matrix on `blocks`.

The same reverse in the forward's terms. With `q = Mb^-1 v1`, `beta = v1' q`
and `P = Mb^-1 - q q' / beta`, a node's complement logdet is `logdet Mb + log
beta` and its step `s = P grad`; `P` is `V (V' Mb V)^-1 V'`, so `dP = -P dMb
P`. Each point's trace adjoint is then a multiple of its curvature's inverse
plus outer products (`_SoftGamma`), and the eigenvector's own derivative,
`-(M0 - lambda1)^+ dM0 v1`, is a solve in `v1`'s complement
(`_soft_deflated_solve`).
"""
function _laplace_soft_quadrature_reverse_factored(laplace::CTSEMLaplaceObjective,
    z0::Vector{Float64}, objf, curvb, blocks::Vector{CTSEMLaplaceBlock}, G,
    npar::Integer)
    D = length(z0)
    D >= 2 || return nothing
    lo, hi = laplace.soft_blend_lo, laplace.soft_blend_hi
    r0 = objf(z0)
    g0 = Float64(r0.value)
    isfinite(g0) || return nothing
    M0 = curvb(z0)
    F0 = _soft_shifted_factor(M0, blocks, 0.0)
    F0 === nothing && return nothing
    thetabar = zeros(Float64, npar)
    plain() = _soft_factored_close!(thetabar, G, z0, 1.0, _SoftGamma(-0.5, F0, blocks),
        zeros(Float64, D), F0, blocks)
    # Above the gate the forward rule returns the plain term.
    _soft_shifted_factor(M0, blocks, hi) === nothing || return plain()
    eig = _soft_smallest_eigenpair(M0, blocks, F0)
    eig === nothing && return nothing
    lam1, v1 = eig.lambda, eig.v
    _soft_second_above(M0, blocks, v1, lam1, _soft_degenerate_bar(lam1)) ||
        return nothing
    w, wp = if lam1 <= lo
        (1.0, 0.0)
    elseif lam1 >= hi
        (0.0, 0.0)
    else
        s = (lam1 - lo) / (hi - lo)
        (1 - (3 * s^2 - 2 * s^3), -(6 * s - 6 * s^2) / (hi - lo))
    end
    iszero(w) && return plain()
    top = _laplace_soft_node_curvature_floor_top()
    _soft_second_above(M0, blocks, v1, lam1, top) || return nothing
    lap = g0 - F0.logdet / 2
    lam1c = max(lam1, 1.0)
    xs, ws = _gauss_hermite(3)
    nn = length(xs)
    terms = zeros(Float64, nn)
    nodes = Vector{Any}(undef, nn)
    for j in 1:nn
        x = xs[j]
        if abs(x) < 1e-12
            h = g0 - (F0.logdet - log(lam1)) / 2
            nodes[j] = (kind=:centre,)
        else
            a = sqrt(2) * x / sqrt(lam1c)
            ub = z0 .+ a .* v1
            rb = objf(ub)
            rgrad = Vector{Float64}(rb.gradient)
            Mb = curvb(ub)
            _soft_shifted_factor(Mb, blocks, top) === nothing && return nothing
            Fb = _soft_shifted_factor(Mb, blocks, 0.0)
            Fb === nothing && return nothing
            q = _soft_solve(Fb, blocks, v1)
            beta = dot(v1, q)
            y = _soft_solve(Fb, blocks, rgrad)
            s = y .- q .* (dot(v1, y) / beta)
            zz = ub .+ s
            hb, hs = Float64(rb.value), Float64(objf(zz).value)
            h = _laplace_soft_node_stepped(hb, hs) - (Fb.logdet + log(beta)) / 2
            nodes[j] = (kind=:node, a=a, ub=ub, Fb=Fb, q=q, beta=beta, y=y, s=s, zz=zz,
                ps=_laplace_soft_node_stepped_slope(hb, hs))
        end
        terms[j] = h + log(ws[j]) + x^2
    end
    peak = maximum(terms)
    isfinite(peak) || return nothing
    ex = [exp(t - peak) for t in terms]
    total = sum(ex)
    soft = peak + log(total) - log(lam1c) / 2 - log(pi) / 2
    # Reverse.
    lapbar, softbar, wbar = 1 - w, w, soft - lap
    lam1bar = wbar * wp
    lam1 > 1.0 && (lam1bar -= softbar / (2 * lam1c))
    gbar0, c0 = 0.0, 0.0
    v1bar = zeros(Float64, D)
    zbar0 = zeros(Float64, D)
    for j in 1:nn
        hbar = softbar * ex[j] / total
        node = nodes[j]
        if node.kind === :centre
            # -(logdet M0 - log lambda1) / 2.
            gbar0 += hbar
            c0 -= hbar / 2
            lam1bar += hbar / (2 * lam1)
            continue
        end
        # The stepped point's share, as in the dense reverse.
        hz = hbar * node.ps
        zzbar = iszero(hz) ? zeros(Float64, D) :
            _laplace_soft_quadrature_vjp!(thetabar, G, node.zz, hz, nothing, nothing)
        zzbar === nothing && return nothing
        q, beta, y, s = node.q, node.beta, node.y, node.s
        ys = _soft_solve(node.Fb, blocks, zzbar)
        t = ys .- q .* (dot(v1, ys) / beta)
        # -(hbar / 2) P from the logdet, -sym(s t') from the step.
        Gb = _SoftGamma(-hbar / 2, node.Fb, blocks, [hbar / (2 * beta), -1.0], [q, s],
            [q, t])
        ub2 = _laplace_soft_quadrature_vjp!(thetabar, G, node.ub, hbar - hz, t, Gb)
        ub2 === nothing && return nothing
        ubbar = zzbar .+ ub2
        # v1 through `log beta` and through P.
        v1y, v1ys = dot(v1, y), dot(v1, ys)
        @inbounds for i in 1:D
            v1bar[i] += -hbar * q[i] / beta - (v1y * ys[i] + v1ys * y[i]) / beta +
                2 * v1ys * v1y * q[i] / beta^2 + node.a * ubbar[i]
            zbar0[i] += ubbar[i]
        end
        lam1 > 1.0 && (lam1bar += dot(v1, ubbar) * (-node.a / (2 * lam1c)))
    end
    # The eigenvector's derivative, -(M0 - lambda1)^+ dM0 v1.
    vperp = v1bar .- v1 .* dot(v1, v1bar)
    y1 = _soft_deflated_solve(eig, blocks, vperp)
    y1 === nothing && return nothing
    gbar0 += lapbar
    c0 -= lapbar / 2
    G0 = _SoftGamma(c0, F0, blocks, [lam1bar, -1.0], [v1, y1], [v1, v1])
    return _soft_factored_close!(thetabar, G, z0, gbar0, G0, zbar0, F0, blocks)
end

"""`_laplace_soft_quadrature_close!` on a factorization: the mode's own terms,
then `dz0 = M0^-1 d(grad g)/dtheta` by a block solve."""
function _soft_factored_close!(thetabar, G, z0, gbar0, Gamma0::_SoftGamma, zbar0, F0,
    blocks)
    zb = _laplace_soft_quadrature_vjp!(thetabar, G, z0, gbar0, nothing, Gamma0)
    zb === nothing && return nothing
    zbar0 .+= zb
    beta = _soft_solve(F0, blocks, zbar0)
    _laplace_soft_quadrature_vjp!(thetabar, G, z0, 0.0, beta, nothing) === nothing &&
        return nothing
    return thetabar
end

"""
    _laplace_soft_quadrature_unit_adjoint(laplace, U, theta)

Unit `U`'s soft-quadrature term's gradient by `_laplace_soft_quadrature_reverse`, in the unit's own
coordinates or, for a filter-integrating objective, the joint ones; `nothing`
where the reverse pass does not apply.
"""
function _laplace_soft_quadrature_unit_adjoint(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64})
    st = _laplace_soft_quadrature_setup(laplace, U, theta, laplace.modes[U])
    st === nothing && return nothing
    if _LAPLACE_SOFT_QUADRATURE_FACTORED[]
        g = _laplace_soft_quadrature_reverse_factored(laplace, st.z0, st.objf, st.curvb,
            st.blocks, st.G, length(theta))
        g === nothing || return g
    end
    return _laplace_soft_quadrature_reverse(laplace, st.z0, st.objf, st.curvf, st.G,
        length(theta))
end

"""
    _laplace_soft_quadrature_setup(laplace, U, theta, uhat)

What the soft-direction quadrature's reverse pass reads at unit `U`'s mode `uhat`: the
coordinates `z0` -- the unit's own, or for an objective whose filter integrates
deviations the joint ones, the deviations at their modes -- the unit objective
`objf` and its curvature there, as a dense matrix (`curvf`) and as a block
matrix (`curvb`) on `blocks`, and the passes' `G` (`_SoftQuadratureThetaGradient`);
`nothing` where the unit has no coordinate or a member cannot be evaluated.
Shared by the soft-quadrature term's gradient
(`_laplace_soft_quadrature_unit_adjoint`) and a mode's fade in its unit's sum
(`_laplace_mode_fade`).
"""
function _laplace_soft_quadrature_setup(laplace::CTSEMLaplaceObjective, U::Integer,
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
        # `let` wherever a closure captures a name both branches assign, or the
        # name is boxed.
        curvb = let layout = layout
            zz -> _hybrid_joint_curvature(laplace, U, theta, Ls, zz, layout)
        end
        curvf = let curvb = curvb
            zz -> curvb(zz).diag[1]
        end
        blocks = [CTSEMLaplaceBlock(0, layout.n, collect(eachindex(members)), 1, Int[])]
    else
        layout = nothing
        z0 = copy(uhat)
        blocks = laplace.units.blocks[U]
        objf = zz -> _laplace_unit_objective_gradient(laplace, U, theta, Ls, zz, aws)
        curvb = zz -> _laplace_unit_curvature(laplace, U, theta, Ls, zz)
        curvf = let curvb = curvb, blocks = blocks
            zz -> _laplace_block_dense(curvb(zz), blocks, length(zz))
        end
    end
    G = _SoftQuadratureThetaGradient(laplace, Int(U), theta, Ls, dL, positions, layout)
    return (z0=z0, objf=objf, curvf=curvf, curvb=curvb, blocks=blocks, G=G)
end

"""
    _laplace_mode_fade(laplace, U, theta, uhat; gradient=false)

The smallest eigenvalue `lam` of unit `U`'s curvature at its mode `uhat`, in
the coordinates the soft-direction quadrature decides on (`_laplace_soft_quadrature_setup`), and with
`gradient` its derivative with respect to the raw vector, the mode's own
movement included: `v' dM v` along its eigenvector, the soft-direction quadrature's closing
pass seeded with `v v'` (`_laplace_soft_quadrature_close!`). What fades a mode out of its
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
(`_laplace_mode_sum`), in the soft-direction quadrature's coordinates
(`_laplace_soft_quadrature_setup`): for each, its point `z0`, the unit objective `g`
there and its curvature's eigendecomposition `E` (`geos`); for each pair, the
valley between them, the lowest `g` on the segment between them
(`valley[k, j]`), and where on the segment it lies, as the fraction of the way
from `k` (`vt[k, j]`). Numbers only, since the setup's closures hold the
calling task's workspace. `nothing` where a mode's quadratic cannot be had.

The valley is the segment's minimum (`_laplace_segment_valley`), and at it the
valley's derivative is the one its point gives with the point held, which is
what `_laplace_mode_sum_gradient!` reads.
"""
function _laplace_modes_geometry(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, modes::AbstractVector)
    K = length(modes)
    geos = Vector{Any}(undef, K)
    objf = nothing
    for k in 1:K
        st = _laplace_soft_quadrature_setup(laplace, U, theta, modes[k])
        st === nothing && return nothing
        k == 1 && (objf = st.objf)
        g = Float64(st.objf(st.z0).value)
        M0 = _dense_sym(Matrix{Float64}(st.curvf(st.z0)))
        (isfinite(g) && all(isfinite, M0)) || return nothing
        geos[k] = (z0=copy(st.z0), g=g, E=_ctsem_symeig(M0))
    end
    valley = fill(Inf, K, K)
    vt = zeros(Float64, K, K)
    for k in 1:K, j in (k + 1):K
        h = let objf = objf, zk = geos[k].z0, zj = geos[j].z0
            t -> _laplace_segment_value(objf, zk, zj, t)
        end
        v = _laplace_segment_valley(h, geos[k].g, geos[j].g)
        valley[k, j] = valley[j, k] = v.value
        vt[k, j] = v.t
        vt[j, k] = 1 - v.t
    end
    return (geos=geos, valley=valley, vt=vt)
end

"""
    _laplace_segment_valley(h, gk, gj)

The valley between two maxima of heights `gk` and `gj` at the ends of a
segment, `h(t)` the objective at the fraction `t` of the way from the first:
the lowest of the seven interior points `_laplace_modes_separated` reads, and
where it may fade either maximum (less than `_LAPLACE_MODE_DEPTH` below the
lower), the lowest of the segment's minima at every dip those points show,
each refined between its point's neighbours (`_laplace_segment_min`).
`(value, t)`.

Not the lowest point alone: seven points are a minimum of seven smooth
functions of the parameters, kinked wherever two of them cross. On the affect
pilot a unit's valley moved from the point at 3/4 to the one at 7/8 just where
a run's certification stalled; its lower maximum's fade and share bent there,
and the Newton step was refused at every radius. Nor the dip of the lowest
point alone, refined: where two dips lie on the segment the lowest point can
pass from one to the other while the deeper stays deeper, and the valley would
jump between their minima. Refined at every dip it moves smoothly, except
where two dips' own minima cross, a bend of the minimum's kind. A dip narrower
than the points' spacing can lie between two of them, unseen until it
reaches one; every dip an eighth of the segment wide is seen.
"""
function _laplace_segment_valley(h, gk::Float64, gj::Float64)
    hs = Float64[h(i / 8) for i in 1:7]
    lo, i0 = findmin(hs)
    at = i0 / 8
    isfinite(lo) || return (value=lo, t=at)
    # A dip deeper than this fades neither maximum, and refining it only
    # deepens it; where the lowest point is one, so is the valley.
    reach = min(gk, gj) - _LAPLACE_MODE_DEPTH
    lo > reach || return (value=lo, t=at)
    for i in 1:7
        hs[i] > reach || continue
        left = i == 1 ? gk : hs[i - 1]
        right = i == 7 ? gj : hs[i + 1]
        (hs[i] <= left && hs[i] <= right) || continue
        m = _laplace_segment_min(h, (i - 1) / 8, i / 8, (i + 1) / 8, hs[i])
        if m.value <= lo
            lo, at = m.value, m.t
        end
    end
    return (value=lo, t=at)
end

# `g` at the fraction `t` of the way from `zk` to `zj`, `-Inf` where it is not
# finite: separate, as `_laplace_modes_separated` reads it.
function _laplace_segment_value(objf, zk::Vector{Float64}, zj::Vector{Float64}, t::Float64)
    v = Float64(objf(zk .+ t .* (zj .- zk)).value)
    return isfinite(v) ? v : -Inf
end

"""
    _laplace_segment_min(h, a, x, b, hx)

The minimum of `h` on `[a, b]` from `x` inside it, `h(x) = hx` no higher than
at either end: Brent's method, golden sections and parabolic steps, to
`sqrt(eps)` in `t`, the closest a minimum's place can be read from its values.
`(t, value)`; `value` is `-Inf` where `h` was not finite at a trial.
"""
function _laplace_segment_min(h, a::Float64, x::Float64, b::Float64, hx::Float64)
    golden = 0.3819660112501051
    tol = sqrt(eps(Float64))
    w = v = x
    hw = hv = hx
    d = e = 0.0
    for _ in 1:100
        mid = (a + b) / 2
        tol1 = tol * abs(x) + 1e-12
        abs(x - mid) <= 2 * tol1 - (b - a) / 2 && break
        parabolic = false
        if abs(e) > tol1
            r = (x - w) * (hx - hv)
            q = (x - v) * (hx - hw)
            p = (x - v) * q - (x - w) * r
            q = 2 * (q - r)
            q > 0 && (p = -p)
            q = abs(q)
            if abs(p) < abs(q * e / 2) && p > q * (a - x) && p < q * (b - x)
                e, d = d, p / q
                (x + d - a < 2 * tol1 || b - x - d < 2 * tol1) &&
                    (d = copysign(tol1, mid - x))
                parabolic = true
            end
        end
        if !parabolic
            e = x >= mid ? a - x : b - x
            d = golden * e
        end
        u = abs(d) >= tol1 ? x + d : x + copysign(tol1, d)
        hu = h(u)
        isfinite(hu) || return (t=u, value=-Inf)
        if hu <= hx
            u >= x ? (a = x) : (b = x)
            v, hv, w, hw, x, hx = w, hw, x, hx, u, hu
        else
            u < x ? (a = u) : (b = u)
            if hu <= hw || w == x
                v, hv, w, hw = w, hw, u, hu
            elseif hu <= hv || v == x || v == w
                v, hv = u, hu
            end
        end
    end
    return (t=x, value=hx)
end

"""
    _laplace_mode_point_gradient(laplace, U, theta, uhat, z)

Unit `U`'s objective's gradient with respect to the raw vector and to the
coordinates at the point `z` of the soft-direction quadrature's coordinates, set up at its
maximum `uhat` (`_laplace_soft_quadrature_setup`): `(theta, z)`, or `nothing` where it
cannot be had. What a valley between two maxima moves by (`_laplace_mode_sum`).
"""
function _laplace_mode_point_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64}, z::Vector{Float64})
    st = _laplace_soft_quadrature_setup(laplace, U, theta, uhat)
    st === nothing && return nothing
    r = st.G(z)
    (r === nothing || !all(isfinite, r.theta) || !all(isfinite, r.z)) && return nothing
    return r
end

"""
    _laplace_mode_close(laplace, U, theta, uhat, E, gbar, Mbar, zbar)

The derivative with respect to the raw vector of `gbar g + tr(Mbar M) +
zbar' z0` at unit `U`'s mode `uhat`, in the soft-direction quadrature's coordinates
(`_laplace_soft_quadrature_setup`), the mode's own movement included: the soft-direction quadrature's
closing pass (`_laplace_soft_quadrature_close!`), `E` the curvature's eigendecomposition
there (`_laplace_modes_geometry`). `nothing` where the curvature is not
positive definite, a pass fails or the result is not finite.
"""
function _laplace_mode_close(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, uhat::Vector{Float64}, E, gbar::Float64,
    Mbar::AbstractMatrix{Float64}, zbar::AbstractVector{Float64})
    all(>(0), E.values) || return nothing
    st = _laplace_soft_quadrature_setup(laplace, U, theta, uhat)
    st === nothing && return nothing
    M0inv = _dense_spectral(E.vectors, [1 / l for l in E.values])
    grad = _laplace_soft_quadrature_close!(zeros(Float64, length(theta)), st.G, st.z0, gbar,
        Matrix{Float64}(Mbar), collect(Float64, zbar), M0inv)
    (grad === nothing || !all(isfinite, grad)) && return nothing
    return grad
end
