# The gated floor for deviations the subject filter integrates (`CTSEMHybrid`,
# hybrid_filter.jl): after laplace.jl, whose types it reads.

################################################################################
# Joint coordinates for the gated floor
################################################################################
#
# The profile reproduces Laplace over the random effects and the deviations
# jointly (see hybrid_filter.jl) wherever the plain Laplace term applies.
# The gated floor also asks for the joint curvature's smallest eigenvalue, and
# the profile's can be larger: it is the Schur complement's. So a unit is
# gated on the joint curvature, assembled from the profile curvature `S`, each
# member's `M` and the cross curvature `H = -d(b - M dhat)/du` this file's
# buffers leave on the profile pass:
#
#     joint - c I  is positive definite  iff
#     S - c I - sum_m H_m' [(M_m - c I)^-1 - M_m^-1] H_m  is,
#
# and a unit it flags is scored by the gated rule itself in joint
# coordinates `z = (u, d_1, ..., d_n)`, `d_m` member `m`'s observed columns,
# with the filter at fixed deviations (`_HYBRID_FIXED`). That is the layer's
# own computation, so the result is the layer's.
#
# The gate runs for every unit on every evaluation, so its record and its test
# live in per-slot scratch and allocate nothing after a slot's first unit of a
# given shape. The joint coordinates are for flagged units only, a few percent
# of them, and are written for clarity over allocation.

"""
    CTSEMHybridCross{T}

One unit's cross curvature with its members' filter-integrated deviations, read
by `_hybrid_after!` from the curvature pass: per member, the observed column
count `q`, `H` (`-d(b - M dhat)/du` over the unit's `d` coordinates), `A`
(`M`) and `dhat` (the conditional mode), all at the outer element type. Each
array is capacity sized and read over its leading `q` rows; `_hybrid_cross!`
hands out a slot's own, reset.
"""
mutable struct CTSEMHybridCross{T}
    d::Int
    q::Vector{Int}
    H::Vector{Matrix{T}}
    A::Vector{Matrix{T}}
    dhat::Vector{Vector{T}}
end
CTSEMHybridCross{T}(d::Integer, nmembers::Integer) where {T} = CTSEMHybridCross{T}(
    Int(d), zeros(Int, nmembers), [zeros(T, 0, d) for _ in 1:nmembers],
    [zeros(T, 0, 0) for _ in 1:nmembers], [T[] for _ in 1:nmembers])

"""The calling slot's `CTSEMHybridCross` for a unit of `d` coordinates and
`nmembers` members, with every count zeroed and every `H` cleared: a sweep
fills only the columns it seeded, and the rest must read zero."""
function _hybrid_cross!(laplace::CTSEMLaplaceObjective, ::Type{T}, d::Integer,
    nmembers::Integer, tag::Symbol=:hybrid_cross) where {T}
    slot = _laplace_slot()
    _laplace_check_slot(laplace, slot)
    store = laplace.workspaces[slot]
    key = (tag, T, Int(d), Int(nmembers))
    cross = get(store, key, nothing)
    if cross === nothing
        cross = CTSEMHybridCross{T}(d, nmembers)
        store[key] = cross
    end
    cross = cross::CTSEMHybridCross{T}
    fill!(cross.q, 0)
    for H in cross.H
        fill!(H, zero(T))
    end
    return cross
end

# Hooks around one member's pass in `_laplace_unit_loglik_gradient`; nothing
# for no hook.
_hybrid_before!(::Nothing, m, ws) = nothing
_hybrid_after!(::Nothing, m, ws, u) = nothing
_hybrid_before!(::CTSEMHybridCross, m, ws) = nothing

function _hybrid_after!(cross::CTSEMHybridCross{T}, m, ws, u) where {T}
    hb = ws.ekf_ws.hyb
    q = hb.active ? hb.q : 0
    cross.q[m] = q
    q == 0 && return nothing
    nh = length(hb.gv)
    if size(cross.H[m], 1) < q
        cross.H[m] = zeros(T, nh, cross.d)
        cross.A[m] = zeros(T, nh, nh)
        cross.dhat[m] = zeros(T, nh)
    end
    A, dhat, H = cross.A[m], cross.dhat[m], cross.H[m]
    @inbounds for b in 1:q
        dhat[b] = _hybrid_strip(hb.dhat[b])
        for a in 1:q
            A[a, b] = _hybrid_strip(hb.Mh[a, b])
        end
    end
    eltype(u) <: ForwardDiff.Dual || return nothing
    # Each coordinate's seed, wherever this sweep put one: a chunked Jacobian
    # or one block's sweep seeds some coordinates and not others, and each
    # sweep fills the columns it seeded.
    @inbounds for a in 1:cross.d
        p = ForwardDiff.partials(u[a])
        for j in eachindex(p)
            iszero(_laplace_deepvalue(p[j])) && continue
            for i in 1:q
                H[i, a] = -ForwardDiff.partials(hb.r[i])[j]
            end
            break
        end
    end
    return nothing
end

"""
    _hybrid_joint_exceeds(laplace, M, blocks, cross, shift)

Whether the joint curvature less `shift I` is positive definite, from the
profile curvature `M` and `cross` (the identity at the top of this file). A
decision, so on the deep values, in the calling slot's scratch.
"""
function _hybrid_joint_exceeds(laplace::CTSEMLaplaceObjective, M::CTSEMBlockMatrix,
    blocks::Vector{CTSEMLaplaceBlock}, cross::CTSEMHybridCross, shift::Real)
    d = cross.d
    D = _laplace_scratch_matrix!(laplace, Float64, d, d, :hybrid_joint_D)
    fill!(D, 0.0)
    @inbounds for (b, block) in enumerate(blocks)
        ob = block.offset
        Db = M.diag[b]
        for j in 1:block.size, i in 1:block.size
            D[ob + i, ob + j] = _laplace_deepvalue(Db[i, j])
        end
        for (t, a) in enumerate(block.ancestors)
            oa = blocks[a].offset
            Cb = M.coupling[b][t]
            for j in 1:blocks[a].size, i in 1:block.size
                v = _laplace_deepvalue(Cb[i, j])
                D[ob + i, oa + j] = v
                D[oa + j, ob + i] = v
            end
        end
    end
    @inbounds for i in 1:d
        D[i, i] -= shift
    end
    for m in eachindex(cross.q)
        q = cross.q[m]
        q == 0 && continue
        cap = size(cross.A[m], 1)
        A = _laplace_scratch_matrix!(laplace, Float64, cap, cap, :hybrid_joint_A)
        As = _laplace_scratch_matrix!(laplace, Float64, cap, cap, :hybrid_joint_As)
        KH = _laplace_scratch_matrix!(laplace, Float64, cap, d, :hybrid_joint_KH)
        col = _laplace_scratch_vector!(laplace, Float64, cap, :hybrid_joint_col)
        a1 = _laplace_scratch_vector!(laplace, Float64, cap, :hybrid_joint_a1)
        a2 = _laplace_scratch_vector!(laplace, Float64, cap, :hybrid_joint_a2)
        Am, H = cross.A[m], cross.H[m]
        @inbounds for j in 1:q, i in 1:q
            v = _laplace_deepvalue(Am[i, j])
            A[i, j] = v
            As[i, j] = i == j ? v - shift : v
        end
        Fa = _ctsem_cholesky(A, q)
        Fs = _ctsem_cholesky(As, q)
        (issuccess(Fa) && issuccess(Fs)) || return false
        # K H with K = (A - shift I)^-1 - A^-1, a column of H at a time.
        @inbounds for j in 1:d
            for i in 1:q
                col[i] = _laplace_deepvalue(H[i, j])
            end
            ldiv!(a1, Fs, col)
            ldiv!(a2, Fa, col)
            for i in 1:q
                KH[i, j] = a1[i] - a2[i]
            end
        end
        @inbounds for j in 1:d, i in 1:d
            acc = 0.0
            for k in 1:q
                acc += _laplace_deepvalue(H[k, i]) * KH[k, j]
            end
            D[i, j] -= acc
        end
    end
    return issuccess(_ctsem_cholesky(D, d))
end

# Collects each member's conditional mode, at the pass's own element type,
# from a pass that needs nothing else.
struct _HybridModeHook{T}
    q::Vector{Int}
    dhat::Vector{Vector{T}}
end
_hybrid_before!(::_HybridModeHook, m, ws) = nothing
function _hybrid_after!(h::_HybridModeHook{T}, m, ws, u) where {T}
    hb = ws.ekf_ws.hyb
    q = hb.active ? hb.q : 0
    h.q[m] = q
    h.dhat[m] = T[hb.dhat[c] for c in 1:q]
    return nothing
end

# Hands each member its fixed deviations from `z`, and collects the deviations'
# own gradient `b - (M - I) d` from its pass.
struct _HybridJointHook{T,Z,L}
    z::Z
    layout::L
    grad::Vector{T}
end

function _hybrid_before!(h::_HybridJointHook, m, ws)
    hb = ws.ekf_ws.hyb
    q = h.layout.q[m]
    length(hb.dfix) >= q || (hb.dfix = zeros(eltype(hb.dfix), q))
    @inbounds for c in 1:q
        hb.dfix[c] = h.z[h.layout.offs[m] + c]
    end
    return nothing
end

function _hybrid_after!(h::_HybridJointHook, m, ws, u)
    hb = ws.ekf_ws.hyb
    q = h.layout.q[m]
    (hb.active ? hb.q : 0) == q || error(string("member ", m, " integrates ",
        hb.q, " deviations in joint coordinates laid out for ", q))
    base = h.layout.offs[m] - h.layout.d
    @inbounds for c in 1:q
        acc = hb.gv[c] * hb.s[c]
        for e in 1:q
            acc -= hb.gv[c] * hb.S[c, e] * hb.gv[e] * hb.dfix[e]
        end
        h.grad[base + c] = acc
    end
    return nothing
end

"""
    _hybrid_joint_objective_gradient(laplace, U, values, Ls, z, aws, layout)

The unit's joint objective at `z = (u, d_1, ...)` -- its members' log
likelihoods at fixed deviations, less `z'z / 2` -- and its gradient, as
`_laplace_unit_objective_gradient` returns them (`hlogdet` zero: nothing is
profiled).
"""
function _hybrid_joint_objective_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, z::AbstractVector{T},
    aws, layout) where {T}
    d = layout.d
    u = T[z[a] for a in 1:d]
    hook = _HybridJointHook{T,typeof(z),typeof(layout)}(z, layout, zeros(T, layout.n - d))
    r = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws,
        eachindex(laplace.units.members[U]), nothing, _HYBRID_FIXED, hook)
    isfinite(r.value) ||
        return (value=r.value, gradient=fill(T(NaN), length(z)), hlogdet=zero(T))
    gradient = vcat(r.gradient, hook.grad)
    @inbounds for a in eachindex(z)
        gradient[a] -= z[a]
    end
    return (value=r.value - dot(z, z) / 2, gradient=gradient, hlogdet=zero(T))
end

# Fixed deviations handed in, the cross curvature at them read out: the two
# hooks one curvature pass in joint coordinates needs.
struct _HybridFixedCross{J,C}
    joint::J
    cross::C
end
_hybrid_before!(h::_HybridFixedCross, m, ws) = _hybrid_before!(h.joint, m, ws)
_hybrid_after!(h::_HybridFixedCross, m, ws, u) = _hybrid_after!(h.cross, m, ws, u)

"""
    _hybrid_joint_curvature(laplace, U, values, Ls, z, layout)

`I - d2(log likelihood)/dz dz` in joint coordinates, as one dense block, from a
curvature pass over the unit's own `d` coordinates alone. The deviations enter
the log likelihood exactly quadratically for given `u`, so their block is `M`
and their cross block `-d(b - M d)/du`, both read from that same pass with the
filter at the fixed deviations: one sweep of the unit's width rather than one
of the joint width, which is most of the gated rule's cost here.
"""
function _hybrid_joint_curvature(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, z::AbstractVector{T},
    layout) where {T}
    d = layout.d
    members = laplace.units.members[U]
    blocks = laplace.units.blocks[U]
    u = T[z[a] for a in 1:d]
    cross = _hybrid_cross!(laplace, T, d, length(members), :hybrid_joint_cross)
    hook = _HybridFixedCross(_HybridJointHook{T,typeof(z),typeof(layout)}(z, layout,
        T[]), cross)
    Muu = _laplace_unit_curvature(laplace, U, values, Ls, u; hook=hook,
        mode=_HYBRID_FIXED)
    D = layout.n
    dense = zeros(T, D, D)
    @inbounds for (b, block) in enumerate(blocks)
        ob = block.offset
        for j in 1:block.size, i in 1:block.size
            dense[ob + i, ob + j] = Muu.diag[b][i, j]
        end
        for (t, a) in enumerate(block.ancestors)
            oa = blocks[a].offset
            C = Muu.coupling[b][t]
            for j in 1:blocks[a].size, i in 1:block.size
                dense[ob + i, oa + j] = C[i, j]
                dense[oa + j, ob + i] = C[i, j]
            end
        end
    end
    @inbounds for m in eachindex(members)
        q = layout.q[m]
        q == 0 && continue
        cross.q[m] == q || error(string("member ", m, " integrates ", cross.q[m],
            " deviations in joint coordinates laid out for ", q))
        o = layout.offs[m]
        H, A = cross.H[m], cross.A[m]
        for i in 1:q
            for a in 1:d
                dense[o + i, a] = H[i, a]
                dense[a, o + i] = H[i, a]
            end
            for j in 1:q
                dense[o + i, o + j] = A[i, j]
            end
        end
    end
    return CTSEMBlockMatrix{T}([dense], [Matrix{T}[]])
end

"""
    _hybrid_gated_term(laplace, U, values, Ls, u, aws; M, logdetM, inner, cross,
                       flagged)

`_laplace_gated_term` for an objective whose filter integrates deviations: the
plain term with the joint logdet where the joint curvature clears the gate,
and the gated rule in joint coordinates where it does not -- what the layer
holding the deviations as coordinates computes, unit for unit. `flagged` says
the unit is already known to be gated (the gated gradient, differentiating a
unit its primal flagged), which needs the conditional modes and not the
curvature the decision reads.
"""
function _hybrid_gated_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws; M=nothing, logdetM=nothing, inner=nothing, cross=nothing,
    flagged::Bool=false) where {T}
    blocks = laplace.units.blocks[U]
    members = laplace.units.members[U]
    if flagged
        modes = _HybridModeHook{T}(zeros(Int, length(members)),
            [T[] for _ in eachindex(members)])
        r = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws,
            eachindex(members), nothing, _HYBRID_PROFILE, modes)
        isfinite(r.value) || return (value=r.value, flagged=true, floor=false)
        return _hybrid_joint_gated(laplace, U, values, Ls, u, aws, length(u),
            modes.q, modes.dhat)
    end
    inner === nothing &&
        (inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws))
    g = inner.value
    (isfinite(g) && !isempty(u)) || return (value=g, flagged=false, floor=true)
    if M === nothing || cross === nothing
        cross = _hybrid_cross!(laplace, T, length(u), length(members))
        M = _laplace_unit_curvature(laplace, U, values, Ls, u; hook=cross)
        logdetM = nothing
    end
    ok = true
    if logdetM === nothing
        ok, logdetM, _, _ = _laplace_block_factor(M, blocks)
        logdetM += 2 * inner.hlogdet
    end
    hi = laplace.gate_hi
    ok && _hybrid_joint_exceeds(laplace, M, blocks, cross, hi) &&
        return (value=g - logdetM / 2, flagged=false, floor=false)
    return _hybrid_joint_gated(laplace, U, values, Ls, u, aws, cross.d, cross.q,
        cross.dhat)
end

"""
    _hybrid_joint_gated(laplace, U, values, Ls, u, aws, d, q, dhat)

The gated rule in joint coordinates at the joint mode `(u, dhat_1, ...)`, `q`
each member's observed columns.
"""
function _hybrid_joint_gated(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws, d::Int, q::Vector{Int}, dhat) where {T}
    members = laplace.units.members[U]
    offs = zeros(Int, length(q))
    at = d
    for m in eachindex(q)
        offs[m] = at
        at += q[m]
    end
    layout = (d=d, q=copy(q), offs=offs, n=at)
    z = Vector{T}(undef, layout.n)
    @inbounds for a in 1:layout.d
        z[a] = u[a]
    end
    @inbounds for m in eachindex(members), c in 1:layout.q[m]
        z[layout.offs[m] + c] = dhat[m][c]
    end
    objective = zz -> _hybrid_joint_objective_gradient(laplace, U, values, Ls, zz,
        aws, layout)
    curvature = zz -> _hybrid_joint_curvature(laplace, U, values, Ls, zz, layout)
    jblocks = [CTSEMLaplaceBlock(0, layout.n, collect(eachindex(members)), 1, Int[])]
    # The layer's unit has a coordinate for every own variance of every member,
    # observed or not, and its eigenvalue cap counts them all.
    dlayer = length(u) + length(members) * length(laplace.objective.params.hybrid.manifests)
    return _laplace_gated_term(laplace, U, values, Ls, z, aws; M=curvature(z),
        inner=objective(z), objective=objective, curvature=curvature,
        blocks=jblocks, dcheck=dlayer)
end
