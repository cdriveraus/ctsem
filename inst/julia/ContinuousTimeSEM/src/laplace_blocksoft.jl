# A nested unit's leaves integrate their own soft directions.
#
# The soft-direction quadrature (`_laplace_soft_quadrature_term`) corrects
# Laplace where a unit's curvature has a direction softer than the prior's:
# it integrates the softest direction at three nodes and the rest by Laplace at
# each. On a nested unit, persons in a study, that direction is nearly always
# one person's: on the affect pilot's study units 90-97% of it lay in a single
# person's block and at most 3.5% in the study's. Taken over the whole unit,
# that one person put the whole study on the rule -- its nodes and its reverse
# pass over every member, three trace passes a gradient -- where with the
# study's effects as parameters only that person paid it; and where two
# persons were soft, it integrated one mixed direction and left the rest on
# Laplace, or floored the pair.
#
# Laplace's term on a nested unit factors by the blocks: with the leaves
# (blocks no other lies under) eliminated first,
#
#     logdet M = logdet E + sum_b logdet M_bb,
#
# `M_bb` leaf `b`'s curvature given the rest of the unit at its mode and `E`
# the rest's curvature with the leaves eliminated (its Schur complement). Each
# piece is a Laplace approximation of its own: the leaf's integral given the
# rest, and the rest's. So each is corrected by its own soft rule, decided on
# its own eigenvalue:
#
#     T = (1 - w) (T_lap + sum_b S_b) + w T_whole,
#
# `S_b` the soft rule's correction on leaf `b`'s objective given the rest at
# the unit's mode (its blended term less its Laplace term, zero where `M_bb`'s
# eigenvalues are above the band), `w` the band's weight at `E`'s smallest
# eigenvalue, and `T_whole` the whole unit's rule, which takes over as the
# rest turns soft. On the pilot `E`'s smallest eigenvalue was 2.4-10.5 and 0-3
# persons a study were soft given the study, so each study unit takes plain
# Laplace and the rule on those few persons alone.
#
# A leaf's correction depends on the rest only through its members' shifted
# values, `v = theta + L_b y + sum_a L_a z_a`, so its derivative in an
# ancestor's coordinates is `L_a'` its derivative in that level's population
# means, read off the leaf's own reverse pass; what that leaves on the rest is
# carried back through the unit's mode by one solve and one pass, as the block
# sum's is (`_laplace_unit_combo_gradient!`).

const _LAPLACE_BLOCK_SOFT = Ref(true)

"""
    ctsem_set_block_soft!(on)

Whether a nested unit's soft-direction quadrature goes leaf by leaf
(`_laplace_block_soft_term`; `true`, the default) or over the whole unit, the
oracle. Returns the previous setting.
"""
function ctsem_set_block_soft!(on::Bool)
    previous = _LAPLACE_BLOCK_SOFT[]
    _LAPLACE_BLOCK_SOFT[] = on
    return previous
end
export ctsem_set_block_soft!

"""Whether unit `U`'s soft-direction quadrature goes leaf by leaf: a nested
unit, the quadrature on, on the unit's own coordinates (the one-step hybrid
decides on joint ones)."""
_laplace_block_soft_on(laplace::CTSEMLaplaceObjective, U::Integer) =
    _LAPLACE_BLOCK_SOFT[] && _LAPLACE_PRIOR_FLOOR[] &&
    laplace.soft_directions === :quadrature && length(laplace.units.blocks[U]) > 1 &&
    !_laplace_joint_hybrid(laplace)

"""A nested unit's leaves: its blocks with coordinates that no block lies under."""
function _laplace_leaves(blocks::Vector{CTSEMLaplaceBlock})
    under = fill(false, length(blocks))
    for block in blocks, a in block.ancestors
        under[a] = true
    end
    return Int[b for b in eachindex(blocks) if !under[b] && blocks[b].size > 0]
end

"""Whether leaf `b`'s own curvature, `M`'s diagonal block, has an eigenvalue
under the band's top: `M_bb - hi I` does not factor."""
function _laplace_leaf_soft_gate(laplace::CTSEMLaplaceObjective, M::CTSEMBlockMatrix,
    b::Int, hi::Real)
    D = M.diag[b]
    n = size(D, 1)
    A = _laplace_scratch_matrix!(laplace, Float64, n, n, :leaf_gate)
    @inbounds for j in 1:n, i in 1:n
        A[i, j] = (_laplace_deepvalue(D[i, j]) + _laplace_deepvalue(D[j, i])) / 2 -
            (i == j ? hi : 0.0)
    end
    return !_ctsem_cholesky!(A, n)
end

"""
    _laplace_leaf_closures(laplace, U, values, Ls, aws, base, st)

Leaf `st.b`'s objective and curvature in its own coordinates, the rest of the
unit held at `base` (`_laplace_block_sum`'s `h_b`): `objf(y)` its members' log
likelihood less its coordinates' prior, with its gradient and the members'
`hlogdet`, as `_laplace_unit_objective_gradient` returns a unit's; `curvf(y)`
its curvature on the subtree's blocks; `densef(y)` that dense.
"""
function _laplace_leaf_closures(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls, aws, base::AbstractVector{T}, st) where {T}
    point = let base = base, st = st
        y -> begin
            u = Vector{eltype(y)}(undef, length(base))
            @inbounds for i in eachindex(base)
                u[i] = base[i]
            end
            @inbounds for (k, i) in enumerate(st.idx)
                u[i] = y[k]
            end
            u
        end
    end
    objf = let laplace = laplace, U = U, values = values, Ls = Ls, aws = aws, st = st,
            point = point
        y -> begin
            u = point(y)
            r = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws, st.members)
            n = length(st.idx)
            g = Vector{eltype(u)}(undef, n)
            v = r.value
            @inbounds for (k, i) in enumerate(st.idx)
                g[k] = r.gradient[i] - u[i]
                v -= u[i]^2 / 2
            end
            (value=v, gradient=g, hlogdet=r.hlogdet)
        end
    end
    curvf = let laplace = laplace, U = U, values = values, Ls = Ls, st = st, point = point
        y -> _laplace_subtree_curvature(laplace, U, values, Ls, point(y), st)
    end
    densef = let curvf = curvf, st = st
        y -> _laplace_block_dense(curvf(y), st.sub, length(st.idx))
    end
    return (objf=objf, curvf=curvf, densef=densef, point=point)
end

"""
    _laplace_leaf_soft(laplace, U, values, Ls, aws, base, st, Mb)

Leaf `st.b`'s correction `S_b` at `base` (the file's header), `Mb` its
curvature there on the subtree's blocks: the soft rule's term on its objective
given the rest (`_laplace_soft_quadrature_term` with the leaf's closures) less
its Laplace term. Also whether the rule scored it (`flagged`; a near-tied
pair is handed off to the floored Laplace term, still flagged) or floored a
degenerate pair where the curvature does not factor (`floor`), and its Laplace
logdet. `nothing` where its value
is not finite or `Mb` does not factor.
"""
function _laplace_leaf_soft(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls, aws, base::AbstractVector{T}, st,
    Mb::CTSEMBlockMatrix{T}) where {T}
    cl = _laplace_leaf_closures(laplace, U, values, Ls, aws, base, st)
    y = T[base[i] for i in st.idx]
    inner = cl.objf(y)
    isfinite(_laplace_deepvalue(inner.value)) || return nothing
    ok, ld, _, _ = _laplace_block_factor(Mb, st.sub)
    ok || return nothing
    logdetM = ld + 2 * inner.hlogdet
    lap = inner.value - logdetM / 2
    r = _laplace_soft_quadrature_term(laplace, U, values, Ls, y, aws; M=Mb,
        logdetM=logdetM, inner=inner, objective=cl.objf, curvature=cl.curvf,
        blocks=st.sub)
    S = (r.flagged || r.logdet_floor) ? r.value - lap : zero(lap)
    return (S=S, flagged=r.flagged, floor=r.logdet_floor, logdet=logdetM, cl=cl, y=y)
end

"""
    _laplace_rest_lowest(M, blocks, leaves)

The smallest eigenvalue of the rest's curvature with the leaves eliminated,
`E = M_RR - sum_b M_Rb M_bb^-1 M_bR` (`R` every coordinate outside a leaf), and
its eigenvector lifted to the unit, `v = (w on R, -M_bb^-1 M_bR w on each
leaf)`, along which `v' dM v` is the eigenvalue's derivative. `nothing` where
the unit has no rest or a leaf does not factor. Exact to first order in a dual
`M` (`_laplace_eig_first_order`).
"""
function _laplace_rest_lowest(M::CTSEMBlockMatrix{T}, blocks::Vector{CTSEMLaplaceBlock},
    leaves::Vector{Int}) where {T}
    isleaf = fill(false, length(blocks))
    for b in leaves
        isleaf[b] = true
    end
    rest = Int[b for b in eachindex(blocks) if !isleaf[b] && blocks[b].size > 0]
    isempty(rest) && return nothing
    roff = fill(-1, length(blocks))
    r = 0
    for b in rest
        roff[b] = r
        r += blocks[b].size
    end
    E = zeros(T, r, r)
    for b in rest
        k, o = blocks[b].size, roff[b]
        D = M.diag[b]
        @inbounds for j in 1:k, i in 1:k
            E[o + i, o + j] = (D[i, j] + D[j, i]) / 2
        end
        for (t, a) in enumerate(blocks[b].ancestors)
            roff[a] < 0 && continue
            C = M.coupling[b][t]
            oa = roff[a]
            @inbounds for j in 1:blocks[a].size, i in 1:k
                E[o + i, oa + j] = C[i, j]
                E[oa + j, o + i] = C[i, j]
            end
        end
    end
    # Each leaf's `M_bb^-1 M_bR`, kept for the lift.
    X = Vector{Matrix{T}}(undef, length(blocks))
    for b in leaves
        nb = blocks[b].size
        Gb = zeros(T, nb, r)
        for (t, a) in enumerate(blocks[b].ancestors)
            roff[a] < 0 && continue
            C = M.coupling[b][t]
            @inbounds for j in 1:blocks[a].size, i in 1:nb
                Gb[i, roff[a] + j] = C[i, j]
            end
        end
        Xb = _laplace_spd_solve(M.diag[b], Gb)
        Xb === nothing && return nothing
        X[b] = Xb
        @inbounds for j in 1:r, i in 1:r
            acc = zero(T)
            for p in 1:nb
                acc += Gb[p, i] * Xb[p, j]
            end
            E[i, j] -= acc
        end
    end
    @inbounds for j in 1:r, i in (j + 1):r
        s = (E[i, j] + E[j, i]) / 2
        E[i, j] = s
        E[j, i] = s
    end
    eig = _laplace_eig_first_order(E)
    lambda = eig.values[1]
    w = eig.primal.vectors[:, 1]
    d = sum(k.size for k in blocks; init=0)
    v = zeros(Float64, d)
    for b in rest, i in 1:blocks[b].size
        v[blocks[b].offset + i] = w[roff[b] + i]
    end
    for b in leaves
        Xb = X[b]
        @inbounds for i in 1:blocks[b].size
            acc = 0.0
            for j in 1:r
                acc += _laplace_deepvalue(Xb[i, j]) * w[j]
            end
            v[blocks[b].offset + i] = -acc
        end
    end
    return (lambda=lambda, v=v)
end

"""`D \\ B` for a symmetric positive definite `D` by the engine's own Cholesky
(no LAPACK lock: `small_linalg.jl`), `nothing` where `D` does not factor."""
function _laplace_spd_solve(D::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    n = size(D, 1)
    R = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        R[i, j] = (D[i, j] + D[j, i]) / 2
    end
    _ctsem_cholesky!(R, n) || return nothing
    X = Matrix{T}(B)
    m = size(B, 2)
    # R'R X = B: forward with R', back with R (R upper).
    @inbounds for c in 1:m
        for i in 1:n
            acc = X[i, c]
            for k in 1:(i - 1)
                acc -= R[k, i] * X[k, c]
            end
            X[i, c] = acc / R[i, i]
        end
        for i in n:-1:1
            acc = X[i, c]
            for k in (i + 1):n
                acc -= R[i, k] * X[k, c]
            end
            X[i, c] = acc / R[i, i]
        end
    end
    return X
end

"""The band's weight at `lambda` and its slope there (`_laplace_smoothstep_weight`)."""
function _laplace_smoothstep_weight_slope(lambda::Real, lo::Real, hi::Real)
    lambda <= lo && return (1.0, 0.0)
    lambda >= hi && return (0.0, 0.0)
    s = (lambda - lo) / (hi - lo)
    return (1 - (3 * s^2 - 2 * s^3), -(6 * s - 6 * s^2) / (hi - lo))
end

"""
    _laplace_block_soft_term(laplace, U, values, Ls, u, aws, M, logdetM, inner)

Nested unit `U`'s term at its maximum `u` by leaf (the file's header), `M` its
curvature there, `logdetM` its joint logdet and `inner` its objective. Returns
the term, whether any leaf or the rest moved it off plain Laplace (`active`,
the route its gradient takes), and its pieces; `nothing` where a piece cannot
be scored, and the whole unit's rule then decides.
"""
function _laplace_block_soft_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls, u::AbstractVector{T}, aws, M::CTSEMBlockMatrix{T},
    logdetM, inner) where {T}
    blocks = laplace.units.blocks[U]
    lo, hi = laplace.soft_blend_lo, laplace.soft_blend_hi
    lap = inner.value - logdetM / 2
    leaves = _laplace_leaves(blocks)
    S = zero(lap)
    soft = Int[]
    for b in leaves
        _laplace_leaf_soft_gate(laplace, M, b, hi) || continue
        st = _laplace_subtree(blocks, b)
        r = _laplace_leaf_soft(laplace, U, values, Ls, aws, u, st,
            _laplace_subtree_matrix(M, blocks, st))
        r === nothing && return nothing
        (r.flagged || r.floor) || continue
        S += r.S
        push!(soft, b)
    end
    rest = _laplace_rest_lowest(M, blocks, leaves)
    w = rest === nothing ? zero(lap) : _laplace_smoothstep_weight(rest.lambda, lo, hi)
    term = lap + S
    if !iszero(_laplace_deepvalue(w))
        g = _laplace_soft_quadrature_term(laplace, U, values, Ls, u, aws; M=M,
            logdetM=logdetM, inner=inner)
        term = (1 - w) * (lap + S) + w * g.value
    end
    return (term=term, active=!isempty(soft) || !iszero(_laplace_deepvalue(w)),
        soft=soft, w=w)
end

################################################################################
# Gradient
################################################################################

"""
A `_SoftQuadratureThetaGradient` read in a subtree's own coordinates, the rest
of the unit held at `base`: what the soft rule's reverse pass
(`_laplace_soft_quadrature_reverse_factored`) asks of a unit, asked of a leaf.
`G` is restricted to the subtree's members and prior (its `group`, `prior`);
its passes' gradients in the rest's coordinates are dropped, the caller
reading them off the means instead (`_laplace_leaf_soft_gradient!`).
"""
struct _SubtreeThetaGradient{G}
    G::G
    base::Vector{Float64}
    idx::Vector{Int}
    loc::Vector{Int}
end

function _subtree_point(Gs::_SubtreeThetaGradient, y::AbstractVector)
    u = copy(Gs.base)
    @inbounds for (k, i) in enumerate(Gs.idx)
        u[i] = y[k]
    end
    return u
end

function (Gs::_SubtreeThetaGradient)(y)
    g = Gs.G(_subtree_point(Gs, y))
    g === nothing && return nothing
    return (theta=g.theta, z=g.z[Gs.idx])
end

function _soft_quadrature_directional(Gs::_SubtreeThetaGradient, y::Vector{Float64},
    beta::Vector{Float64})
    B = zeros(Float64, length(Gs.base))
    @inbounds for (k, i) in enumerate(Gs.idx)
        B[i] = beta[k]
    end
    r = _soft_quadrature_directional(Gs.G, _subtree_point(Gs, y), B)
    r === nothing && return nothing
    return (theta=r.theta, z=r.z[Gs.idx])
end

function _soft_quadrature_trace!(thetabar::Vector{Float64}, zbar::Vector{Float64},
    Gs::_SubtreeThetaGradient, y::Vector{Float64}, Gamma)
    zf = zeros(Float64, length(Gs.base))
    _soft_quadrature_trace!(thetabar, zf, Gs.G, _subtree_point(Gs, y),
        _SubtreeMemberGamma(Gamma, Gs.loc)) || return false
    @inbounds for (k, i) in enumerate(Gs.idx)
        zbar[k] += zf[i]
    end
    return true
end

"""
A trace term's adjoint `Gamma` on a subtree's coordinates (a dense matrix or a
`_SoftGamma` on the subtree's blocks), read by members whose coordinates are
the unit's: `loc` places a unit coordinate in the subtree, zero outside it,
where the adjoint is zero.
"""
struct _SubtreeMemberGamma{T}
    inner::T
    loc::Vector{Int}
end

function _soft_member_gamma(Gamma::_SubtreeMemberGamma, G::_SoftQuadratureThetaGradient,
    m::Integer)
    S = G.coords[m]
    n = length(S)
    out = zeros(Float64, n, n)
    @inbounds for j in 1:n, i in 1:n
        p, q = Gamma.loc[S[i]], Gamma.loc[S[j]]
        (p == 0 || q == 0) && continue
        out[i, j] = _soft_gamma_entry(Gamma.inner, p, q)
    end
    return out
end

"""Entry `(p, q)` of a trace adjoint, symmetrised: a dense matrix's, or a
`_SoftGamma`'s `c M^-1 + sum_k w_k sym(a_k b_k')`."""
_soft_gamma_entry(Gamma::AbstractMatrix, p::Int, q::Int) = (Gamma[p, q] + Gamma[q, p]) / 2
function _soft_gamma_entry(Gamma::_SoftGamma, p::Int, q::Int)
    acc = iszero(Gamma.c) ? 0.0 : Gamma.c * _soft_inverse_entry(Gamma, p, q)
    for k in eachindex(Gamma.w)
        acc += Gamma.w[k] * (Gamma.a[k][p] * Gamma.b[k][q] + Gamma.b[k][p] * Gamma.a[k][q]) / 2
    end
    return acc
end

"""
    _laplace_leaf_soft_gradient!(thetabar, zbar, laplace, U, theta, Ls, dL,
        positions, aws, base, st, Mb, c)

`c dS_b` (`_laplace_leaf_soft`) into `thetabar` given the rest at `base`, its
movement in the leaf's own maximum included (the soft reverse's close, and
Laplace's), and into `zbar` what it leaves on the rest's coordinates: `c L_a'`
its derivative in level `a`'s population means, for each ancestor `a` (the
file's header). Returns `S_b`, or `nothing` where a piece fails.
"""
function _laplace_leaf_soft_gradient!(thetabar::Vector{Float64}, zbar::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dL, positions, aws, base::Vector{Float64}, st,
    Mb::CTSEMBlockMatrix{Float64}, c::Float64)
    spec = laplace.spec
    blocks = laplace.units.blocks[U]
    npar = length(theta)
    r = _laplace_leaf_soft(laplace, U, theta, Ls, aws, base, st, Mb)
    r === nothing && return nothing
    (r.flagged || r.floor) || return 0.0
    c == 0 && return r.S
    G = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dL, positions, nothing;
        group=collect(Int, st.members), prior=collect(Int, st.idx))
    Gs = _SubtreeThetaGradient(G, base, collect(Int, st.idx), st.loc)
    y = r.y
    tb = zeros(Float64, npar)
    if r.flagged
        g = _laplace_soft_quadrature_reverse_factored(laplace, y, r.cl.objf, r.cl.curvf,
            st.sub, Gs, npar)
        g === nothing && (g = _laplace_soft_quadrature_reverse(laplace, y, r.cl.objf,
            r.cl.densef, Gs, npar))
        g === nothing && return nothing
        tb .+= g
    elseif r.logdet < 0
        # A floored pair: `h - max(logdet, 0) / 2` is `h`, flat in `y` at the
        # leaf's maximum.
        gh = Gs(y)
        gh === nothing && return nothing
        tb .+= gh.theta
    else
        return r.S
    end
    # Less Laplace's `h - logdet / 2`, through the leaf's maximum.
    F = _soft_shifted_factor(Mb, st.sub, 0.0)
    F === nothing && return nothing
    tl = zeros(Float64, npar)
    _soft_factored_close!(tl, Gs, y, 1.0, _SoftGamma(-0.5, F, st.sub),
        zeros(Float64, length(y)), F, st.sub) === nothing && return nothing
    tb .-= tl
    all(isfinite, tb) || return nothing
    @inbounds for t in 1:npar
        thetabar[t] += c * tb[t]
    end
    for a in blocks[st.b].ancestors
        l = blocks[a].level
        rho = spec.levels[l].re_index
        L = Ls[l]
        @inbounds for q in 1:blocks[a].size
            acc = 0.0
            for p in eachindex(rho)
                acc += L[p, q] * tb[rho[p]]
            end
            zbar[blocks[a].offset + q] += c * acc
        end
    end
    return r.S
end

"""
    _laplace_block_soft_gradient!(row, laplace, U, theta, Ls, dLlevels,
        dLsolved, M, curvature)

The gradient of `_laplace_block_soft_term` at the unit's maximum
(`laplace.modes[U]`), `M` its curvature and `curvature` its factorization:
`(1 - w) (dT_lap + sum_b dS_b) + w dT_whole + (T_whole - T_lap - sum_b S_b)
dw`, `T_lap` by the seeded route, each leaf's `S_b` by its own reverse pass,
the rest's coordinates they move carried through the unit's maximum, `T_whole`
by the whole unit's soft reverse, and `dw` through `v' dM v`
(`_laplace_rest_lowest`). `false` where a piece fails.
"""
function _laplace_block_soft_gradient!(row::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dLlevels, dLsolved, M::CTSEMBlockMatrix{Float64},
    curvature)
    spec = laplace.spec
    blocks = laplace.units.blocks[U]
    z = laplace.modes[U]
    npar, d = length(theta), length(z)
    lo, hi = laplace.soft_blend_lo, laplace.soft_blend_hi
    aws = _laplace_workspace!(laplace, Float64, npar)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    factors, elim = curvature
    leaves = _laplace_leaves(blocks)
    rest = _laplace_rest_lowest(M, blocks, leaves)
    w, wp = rest === nothing ? (0.0, 0.0) :
        _laplace_smoothstep_weight_slope(rest.lambda, lo, hi)
    fill!(row, 0.0)
    S = 0.0
    if w < 1
        gS = zeros(Float64, npar)
        zbar = zeros(Float64, d)
        for b in leaves
            _laplace_leaf_soft_gate(laplace, M, b, hi) || continue
            st = _laplace_subtree(blocks, b)
            s = _laplace_leaf_soft_gradient!(gS, zbar, laplace, U, theta, Ls, dLlevels,
                positions, aws, z, st, _laplace_subtree_matrix(M, blocks, st), 1.0)
            s === nothing && return false
            S += s
        end
        if any(!iszero, zbar)
            beta = _laplace_block_solve(factors, elim, blocks, zbar)
            Gw = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dLlevels, positions,
                nothing)
            g = _soft_quadrature_directional(Gw, z, beta)
            g === nothing && return false
            gS .+= g.theta
        end
        _laplace_seeded_unit_gradient!(row, laplace, U, theta, Ls, dLlevels, dLsolved,
            M, factors, elim) || return false
        row .+= gS
        row .*= (1 - w)
    end
    w > 0 || return all(isfinite, row)
    # The whole unit's term as its rule scored it: integrated, or handed off to
    # the floored Laplace term where its softest pair is near a tie (alike
    # members usually are; both flagged), floored, or plain
    # Laplace above the band; each differentiated as such.
    inner = _laplace_unit_objective_gradient(laplace, U, theta, Ls, z, aws)
    ok, ld, _, _ = _laplace_block_factor(M, blocks)
    ok || return false
    logdetM = ld + 2 * inner.hlogdet
    whole = _laplace_soft_quadrature_term(laplace, U, theta, Ls, z, aws; M=M,
        logdetM=logdetM, inner=inner)
    gw = zeros(Float64, npar)
    if whole.flagged
        g = _laplace_soft_quadrature_unit_adjoint(laplace, U, theta)
        g === nothing && return false
        gw .= g
    elseif whole.logdet_floor && logdetM < 0
        # The floor binding: `g - max(logdet, 0) / 2` is `g`.
        _laplace_floored_unit_gradient!(gw, laplace, U, theta, Ls, dLlevels) ||
            return false
    else
        _laplace_seeded_unit_gradient!(gw, laplace, U, theta, Ls, dLlevels, dLsolved,
            M, factors, elim) || return false
    end
    row .+= w .* gw
    if wp != 0
        coef = (whole.value - (inner.value - logdetM / 2) - S) * wp
        glam = _laplace_rest_lambda_gradient(laplace, U, theta, Ls, dLlevels, positions,
            z, M, rest.v)
        glam === nothing && return false
        row .+= coef .* glam
    end
    return all(isfinite, row)
end

"""The derivative of `_laplace_rest_lowest`'s eigenvalue in the raw vector,
`v' dM v` along its lifted eigenvector with the unit's maximum moving: the soft
rule's closing pass seeded with `v v'`."""
function _laplace_rest_lambda_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls, dLlevels, positions, z::Vector{Float64},
    M::CTSEMBlockMatrix{Float64}, v::Vector{Float64})
    blocks = laplace.units.blocks[U]
    F = _soft_shifted_factor(M, blocks, 0.0)
    F === nothing && return nothing
    Gw = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dLlevels, positions, nothing)
    out = zeros(Float64, length(theta))
    _soft_factored_close!(out, Gw, z, 0.0, _SoftGamma(0.0, F, blocks, [1.0], [v], [v]),
        zeros(Float64, length(z)), F, blocks) === nothing && return nothing
    return out
end

################################################################################
# A nested unit whose curvature does not factor
################################################################################
#
# Where the inner solve ends at a point that is not a strict maximum -- in a
# fit, a line search's long trial where some members' dynamics grow a
# hundredfold over a gap -- a unit's curvature `M` does not factor: it is
# indefinite, or a pivot's condition is past what `_laplace_block_factor`
# accepts. The whole unit's soft rule took such a unit densely: four Jacobi
# decompositions of `d x d` and a Gram-Schmidt basis. On the affect pilot's
# study units, `d` over a thousand, one evaluation took one to two hours on one
# thread while the session's other threads waited on it
# (CT-SEM/review/PILOT-stalls-read-by-gdb-2026-11-25.md).
#
# The line search rejects such points, so their term must be finite, cheap and
# continuous, not accurate. A nested unit takes the Laplace term with each pivot
# of its block elimination floored as the soft rule's nodes floor eigenvalues
# (`_laplace_soft_node_curvature_floor`, `f`):
#
#     T = g - (1/2) sum_b logdet f(P_b),
#
# the blocks eliminated innermost first, as `_laplace_block_factor` takes them,
# each pivot `P_b` replaced by `f(P_b)` in the elimination as well
# (`_laplace_pivot_floor_elimination`). Where every pivot's eigenvalues are at
# or above the floor's top this is Laplace's own `g - logdet M / 2`, however
# badly conditioned `M` is; elsewhere a soft or negative pivot direction is
# credited what the floor credits a node's. Its cost is the elimination's,
# linear in the members, and it is smooth in `M` (`f` is C2). Its gradient is
# the elimination's reverse (`_laplace_pivot_floor_adjoint`) through the soft
# rule's trace passes and closing solve (`_soft_factored_close!`).

"""
One pivot of `_laplace_pivot_floor_elimination`: `logdet f(P)` and `f(P)^-1`,
exact to first order in `P`'s dual parts (Daleckii-Krein, as
`_laplace_soft_node_floored_curvature`), and the primal eigenpairs `Q, mu` and
floored values `f` its reverse reads.
"""
struct _PivotFloor{T}
    logdet::T
    inverse::Matrix{T}
    Q::Matrix{Float64}
    mu::Vector{Float64}
    f::Vector{Float64}
end

# `h = 1 / f`'s divided difference between two eigenvalues of a pivot, its
# slope where they coincide (the floor's own coincidence rule,
# `_soft_node_floored_curvature_adjoint`).
function _pivot_floor_inverse_divided(a::Float64, b::Float64)
    f = _laplace_soft_node_curvature_floor
    if abs(a - b) <= 1e-12 * max(1.0, abs(a))
        return -_laplace_soft_node_curvature_floor_slope(a) / f(a)^2
    end
    return (1 / f(a) - 1 / f(b)) / (a - b)
end

"""
    _laplace_pivot_floor(P)

Pivot `P` with each eigenvalue floored (`_laplace_soft_node_curvature_floor`), as
a `_PivotFloor`: read as symmetric through both triangles averaged.
"""
function _laplace_pivot_floor(P::AbstractMatrix{T}) where {T}
    n = size(P, 1)
    P0 = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        P0[i, j] = (Float64(_laplace_deepvalue(P[i, j])) +
            Float64(_laplace_deepvalue(P[j, i]))) / 2
    end
    E = _ctsem_symeig(P0)
    Q, mu = E.vectors, E.values
    f = Float64[_laplace_soft_node_curvature_floor(x) for x in mu]
    # The middle factor `C` of `f(P)^-1 = Q C Q'`: `diag(1 / f)`, plus at a
    # dual `P` its first-order part.
    C = zeros(T, n, n)
    logdet = zero(T)
    @inbounds for i in 1:n
        C[i, i] = 1 / f[i]
        logdet += log(f[i])
    end
    if T !== Float64
        D = Matrix{T}(undef, n, n)
        @inbounds for j in 1:n, i in 1:n
            D[i, j] = (P[i, j] + P[j, i]) / 2 - P0[i, j]
        end
        Dt = transpose(Q) * D * Q
        @inbounds for j in 1:n, i in 1:n
            C[i, j] += _pivot_floor_inverse_divided(mu[i], mu[j]) * Dt[i, j]
        end
        @inbounds for i in 1:n
            logdet += _laplace_soft_node_curvature_floor_slope(mu[i]) / f[i] * Dt[i, i]
        end
    end
    inverse = zeros(T, n, n)
    @inbounds for j in 1:n, k in 1:n, l in 1:n
        c = C[l, k]
        iszero(c) && continue
        qjk = Q[j, k]
        for i in 1:n
            inverse[i, j] += Q[i, l] * c * qjk
        end
    end
    return _PivotFloor{T}(logdet, inverse, Q, mu, f)
end

# `A * B` and `A' * B` by loops, generic in the element types: never BLAS,
# which takes a process-wide lock (small_linalg.jl).
function _pivot_floor_mul(A::AbstractMatrix, B::AbstractMatrix)
    C = zeros(promote_type(eltype(A), eltype(B)), size(A, 1), size(B, 2))
    @inbounds for j in axes(B, 2), k in axes(A, 2)
        b = B[k, j]
        iszero(b) && continue
        for i in axes(A, 1)
            C[i, j] += A[i, k] * b
        end
    end
    return C
end
function _pivot_floor_mulTN(A::AbstractMatrix, B::AbstractMatrix)
    C = zeros(promote_type(eltype(A), eltype(B)), size(A, 2), size(B, 2))
    @inbounds for j in axes(B, 2), i in axes(A, 2)
        acc = zero(eltype(C))
        for k in axes(A, 1)
            acc += A[k, i] * B[k, j]
        end
        C[i, j] = acc
    end
    return C
end

"""
    _laplace_pivot_floor_elimination(M, blocks)

`M` eliminated as `_laplace_block_factor` eliminates it, innermost block first,
with each pivot replaced by its floored form (`_laplace_pivot_floor`) both in
the log determinant and in the Schur updates it makes to its ancestors:
`(logdet, pivots, coupling)`, `logdet` the sum of `logdet f(P_b)`, `pivots`
the `_PivotFloor`s and `coupling` the couplings as the elimination used them.
Generic in the element type. `nothing` where a pivot's eigenvalues are not
finite.
"""
function _laplace_pivot_floor_elimination(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    nb = length(blocks)
    diag = [Matrix{T}(d) for d in M.diag]
    coupling = [[Matrix{T}(c) for c in row] for row in M.coupling]
    pivots = Vector{_PivotFloor{T}}(undef, nb)
    logdet = zero(T)
    for b in 1:nb
        p = _laplace_pivot_floor(diag[b])
        all(isfinite, p.mu) || return nothing
        pivots[b] = p
        logdet += p.logdet
        ancestors = blocks[b].ancestors
        isempty(ancestors) && continue
        W = [_pivot_floor_mul(p.inverse, coupling[b][t]) for t in eachindex(ancestors)]
        for t in eachindex(ancestors)
            a = ancestors[t]
            diag[a] .-= _laplace_symmetrise(_pivot_floor_mulTN(coupling[b][t], W[t]))
            for s in eachindex(ancestors)
                s == t && continue
                slot = findfirst(==(ancestors[s]), blocks[a].ancestors)
                slot === nothing && continue
                coupling[a][slot] .-= _pivot_floor_mulTN(coupling[b][t], W[s])
            end
        end
    end
    return (logdet=logdet, pivots=pivots, coupling=coupling)
end

"""
    _laplace_pivot_floor_adjoint(F, blocks)

The gradient of `_laplace_pivot_floor_elimination`'s `logdet` with respect to
the curvature it eliminated, `F` its Float64 result: `(diag, coupling)` in the
curvature's own block pattern, the derivative in each stored block, a coupling
block standing for both its entries and their transposes. The elimination
reversed, outermost block first: a pivot's share is `Q diag(f'/f) Q'` and,
through its inverse `h(P) = f(P)^-1` in the Schur updates, the Daleckii-Krein
pullback of what those updates' adjoints ask of it; the couplings' shares are
the updates' own.
"""
function _laplace_pivot_floor_adjoint(F, blocks::Vector{CTSEMLaplaceBlock})
    nb = length(blocks)
    Dbar = [zeros(Float64, k.size, k.size) for k in blocks]
    Cbar = [[zeros(Float64, k.size, blocks[a].size) for a in k.ancestors] for k in blocks]
    for b in nb:-1:1
        p = F.pivots[b]
        Q, mu, f = p.Q, p.mu, p.f
        n = length(mu)
        # `logdet f(P)`'s own share.
        ct = zeros(Float64, n, n)
        @inbounds for i in 1:n
            ct[i, i] = _laplace_soft_node_curvature_floor_slope(mu[i]) / f[i]
        end
        ancestors = blocks[b].ancestors
        if !isempty(ancestors)
            C = F.coupling[b]
            G = p.inverse
            Gbar = zeros(Float64, n, n)
            for t in eachindex(ancestors)
                a = ancestors[t]
                # diag[a] -= C_t' G C_t
                Gbar .-= _pivot_floor_mul(_pivot_floor_mul(C[t], Dbar[a]), transpose(C[t]))
                Cbar[b][t] .-= 2 .* _pivot_floor_mul(_pivot_floor_mul(G, C[t]), Dbar[a])
                for s in eachindex(ancestors)
                    s == t && continue
                    slot = findfirst(==(ancestors[s]), blocks[a].ancestors)
                    slot === nothing && continue
                    # coupling[a][slot] -= C_t' G C_s
                    X = Cbar[a][slot]
                    Gbar .-= _pivot_floor_mul(_pivot_floor_mul(C[t], X), transpose(C[s]))
                    Cbar[b][t] .-= _pivot_floor_mul(_pivot_floor_mul(G, C[s]), transpose(X))
                    Cbar[b][s] .-= _pivot_floor_mul(_pivot_floor_mul(G, C[t]), X)
                end
            end
            # Through `G = h(P)`: Daleckii-Krein on the symmetric part.
            Gt = _pivot_floor_mulTN(Q, _pivot_floor_mul(_laplace_symmetrise(Gbar), Q))
            @inbounds for j in 1:n, i in 1:n
                ct[i, j] += _pivot_floor_inverse_divided(mu[i], mu[j]) * Gt[i, j]
            end
        end
        Dbar[b] = _laplace_symmetrise(_pivot_floor_mul(_pivot_floor_mul(Q, ct), transpose(Q)))
    end
    return (diag=Dbar, coupling=Cbar)
end

"""
    _laplace_pivot_floor_term(laplace, U, values, Ls, u, aws, inner)

Nested unit `U`'s term at `u`, where its curvature does not factor: `g - (1/2)
sum_b logdet f(P_b)` (the section's header), from the curvature at `u` formed
afresh, `inner` the unit objective there. `NaN` where the elimination cannot
be made.
"""
function _laplace_pivot_floor_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls, u::AbstractVector{T}, aws, inner) where {T}
    M = _laplace_unit_curvature(laplace, U, values, Ls, u)::CTSEMBlockMatrix{T}
    F = _laplace_pivot_floor_elimination(M, laplace.units.blocks[U])
    F === nothing && return _ctsem_nan(T)
    return inner.value - inner.hlogdet - F.logdet / 2
end

"""
    _laplace_pivot_floor_gradient!(row, laplace, U, theta, Ls, dLlevels)

`_laplace_pivot_floor_term`'s gradient at unit `U`'s mode into `row`: `g`'s own,
and `logdet`'s through `M`'s trace passes, `-Gamma / 2` with `Gamma` the
elimination's adjoint (`_laplace_pivot_floor_adjoint`), closed through the
mode's movement by the implicit function theorem (`_soft_factored_close!`), on
`M`'s own elimination where it is not singular (`_soft_node_factor`) and on
the floored one's otherwise. `false` where a pass fails.
"""
function _laplace_pivot_floor_gradient!(row::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dLlevels)
    spec = laplace.spec
    blocks = laplace.units.blocks[U]
    z = laplace.modes[U]
    M = _laplace_copy_blockmatrix(_laplace_unit_curvature(laplace, U, theta, Ls, z))
    F = _laplace_pivot_floor_elimination(M, blocks)
    F === nothing && return false
    bar = _laplace_pivot_floor_adjoint(F, blocks)
    Gamma = _soft_gamma_entries(-0.5, bar.diag,
        [[c ./ 2 for c in row_] for row_ in bar.coupling], blocks)
    F0 = _soft_node_factor(M, blocks)
    F0 === nothing && (F0 = (logdet=F.logdet,
        factors=[_LDLPivot(p.f, p.Q) for p in F.pivots], coupling=F.coupling))
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    G = _SoftQuadratureThetaGradient(laplace, U, theta, Ls, dLlevels, positions, nothing)
    fill!(row, 0.0)
    _soft_factored_close!(row, G, z, 1.0, Gamma, zeros(Float64, length(z)), F0,
        blocks) === nothing && return false
    return all(isfinite, row)
end
