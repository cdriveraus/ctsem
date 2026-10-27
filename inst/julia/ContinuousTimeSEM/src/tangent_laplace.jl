# A Laplace unit's curvature (`_laplace_columns_jacobian`) by the hand-coded
# second order: each member's sweep carries the unit's coordinates as tangents
# (`tangent_sweep.jl`, `tangent_rows.jl`), through the member's shifted values
# at the curvature's dual (`_laplace_curvature_dual`), each coordinate of a
# chunk a lane, where the other route runs the member's adjoint gradient at
# that dual. A coordinate moves only the cells its effects reach, so a random
# intercept's sweep carries tangents on the mean alone, which the dual carries
# as zeros through every covariance. The gated floor's cross curvature with the
# deviations the filter integrates (`CTSEMHybridCross`) is read off the same
# sweeps. The dual route serves the other hooks, an element type other than
# Float64, and what the sweep does not cover (`_tan_refusal`); it is the oracle
# (`ctsem_set_tangent_hessian!(false)`).

"""A slot's hand-coded curvature: the member's layer at the curvature's dual,
its sweep, and their scratch."""
struct TanLaplaceContext{W,DC,A,R}
    awsW::A
    mod::TanModel{W}
    acc::TanAcc{W}
    rows::R
    sw::TanSweep{W}
    Hb::Matrix{Float64}         # the member's raw columns, nvalues x W
    xw::Vector{DC}              # its shifted values
    uw::Vector{DC}              # the unit's coordinates, seeded
    Lb::Matrix{Float64}         # a hybrid block's factor
end

"""
The slot's hand-coded curvature for `laplace` in hybrid `mode`, built at its
first use from `values`; `nothing` where the sweep does not cover the model.
"""
function _tan_laplace_context(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    mode::Int, slot::Int=_laplace_slot())
    store = laplace.workspaces[slot]
    key = (:tan_curvature, length(values), mode)
    haskey(store, key) && return store[key]
    objective = laplace.objective
    DC = _laplace_curvature_dual(laplace)
    ctx = _tan_refusal(objective, mode) === nothing ?
        _tan_laplace_build(objective, values, mode, DC,
            Val(ForwardDiff.npartials(DC))) : nothing
    store[key] = ctx
    return ctx
end

function _tan_laplace_build(objective::CTSEMObjective, values::AbstractVector, mode::Int,
    ::Type{DC}, ::Val{W}) where {DC,W}
    sp = objective.params
    nv = length(values)
    subjects = objective.subject_objectives
    aws0 = CTSEMAdjointWorkspace(Float64, sp, nv, 0)
    m0 = _tan_layer_mats(aws0, sp, collect(Float64, values), subjects[1].tipreds)
    n, m, ntd = size(m0.JAx, 1), size(m0.LAM, 1), size(m0.TDE, 2)
    na, dyn = aws0.affine_dim, collect(aws0.diffusion_state_indices)
    rows = _tan_rows_model(sp)
    h = sp.hybrid
    nh = _hybrid_on(h) && mode != _HYBRID_OFF ? _hybrid_items(h) : 0
    gmask = any(t -> t != 0, sp.manifesttype) ?
        Bool[k > length(sp.manifesttype) || sp.manifesttype[k] == 0 for k in 1:m] : Bool[]
    mod = TanModel(n, m, W; ct=sp.continuous_time, stat=sp.stationary, rows=rows,
        gmask=gmask, na=na, dyn=dyn, ntd=ntd, hyb=h, nh=nh,
        wlog=mode == _HYBRID_FULL ? 1.0 : 0.0)
    awsW = CTSEMAdjointWorkspace(DC, sp, nv, ntd)
    return TanLaplaceContext(awsW, mod, TanAcc(n, m, ntd, nh, W),
        rows ? TanRows(awsW, sp, n, m, ntd, Val(W), DC) : nothing,
        TanSweep(n, m, ntd, nh, W, na, length(dyn)), zeros(nv, W), zeros(DC, nv),
        DC[], h.block && nh > 0 ? zeros(length(h.perm), length(h.perm)) : zeros(0, 0))
end

"""
`Jw`'s first `length(columns)` columns, the Jacobian of `members`' log
likelihood gradient in the unit's coordinates `columns` at `base`, by the
hand-coded sweeps, the rest zero; `false` where they do not serve, and the
caller takes the dual. A member whose filter fails makes it `NaN`, as the
dual's gradient does. `cross`, where given, is filled as `_hybrid_after!` fills
it from the dual's passes.
"""
function _tan_laplace_columns_jacobian!(Jw::Matrix{Float64}, laplace::CTSEMLaplaceObjective,
    U::Integer, values::AbstractVector{Float64}, Ls::Vector{<:AbstractMatrix},
    base::AbstractVector{Float64}, columns::Vector{Int}, members::Vector{Int}, mode::Int,
    cross::Union{Nothing,CTSEMHybridCross{Float64}}=nothing)
    _CTSEM_TANGENT_HESSIAN[] || return false
    ctx = _tan_laplace_context(laplace, values, mode)
    ctx === nothing && return false
    # A run-time context, so one dynamic call; everything under it is typed.
    return _tan_laplace_columns!(Jw, ctx, laplace, Int(U), values, Ls, base, columns,
        members, mode, cross)::Bool
end

function _tan_laplace_columns!(Jw::Matrix{Float64}, ctx::TanLaplaceContext{W},
    laplace::CTSEMLaplaceObjective, U::Int, values::AbstractVector{Float64},
    Ls::Vector{<:AbstractMatrix}, base::AbstractVector{Float64}, columns::Vector{Int},
    members::Vector{Int}, mode::Int, cross) where {W}
    fill!(Jw, 0.0)
    npos = length(members)
    # The hybrid's factor is the population's, the same for every member and
    # unmoved by the coordinates: values only.
    _tan_laplace_hybrid!(ctx, values) || (fill!(Jw, NaN); return true)
    if _laplace_pool_width() <= 1 || npos < 2 || isempty(_LAPLACE_FREE)
        st = _tan_laplace_members!(Jw, ctx, laplace, U, values, Ls, base, columns,
            members, 1:npos, cross)
        st == _TAN_REFUSED && return false
        st == _TAN_OK || fill!(Jw, NaN)
        return true
    end
    # As `_laplace_unit_loglik_gradient`'s reduction: each worker its own
    # context and columns, summed after the join.
    width = _laplace_pool_width()
    used = fill(false, width)
    status = fill(_TAN_OK, width)
    d, nx = size(Jw)
    _laplace_partition(laplace, npos) do mine, _w
        local slot, c, J
        slot = _laplace_slot()
        c = _tan_laplace_context(laplace, values, mode, slot)::typeof(ctx)
        _tan_laplace_hybrid!(c, values)
        J = _laplace_scratch_matrix!(laplace, Float64, d, nx, :tan_curv_J, slot)
        fill!(J, 0.0)
        used[slot] = true
        status[slot] = _tan_laplace_members!(J, c, laplace, U, values, Ls, base, columns,
            members, mine, cross)
        return status[slot] == _TAN_OK
    end
    any(==(_TAN_REFUSED), status) && return false
    any(!=(_TAN_OK), status) && (fill!(Jw, NaN); return true)
    @inbounds for w in 1:width
        used[w] || continue
        Jw .+= _laplace_scratch_matrix!(laplace, Float64, d, nx, :tan_curv_J, w)
    end
    return true
end

"""The hybrid's factor into the context's model, at `values`; `false` where
its block has none."""
function _tan_laplace_hybrid!(ctx::TanLaplaceContext, values::AbstractVector{Float64})
    mod = ctx.mod
    nh = mod.nh
    nh == 0 && return true
    h = mod.hyb
    G = mod.G
    fill!(G.v, 0.0)
    _tan_setnz!(G, false)
    if !h.block
        for k in 1:nh
            G.v[k, k] = _tan_hybrid_own(h, values, k)[1]
        end
    else
        _hybrid_block_factor!(ctx.Lb, h, values) || return false
        for c in 1:nh, f in 1:c
            G.v[c, f] = _tan_hybrid_block_G(h, ctx.Lb, c, f)
        end
    end
    return true
end

"""
The members at positions `members[mine]` swept along each chunk of `columns`,
their curvature columns pushed through the loadings into `J`, and each one's
cross curvature into `cross` where given. A status code.
"""
function _tan_laplace_members!(J::Matrix{Float64}, ctx::TanLaplaceContext{W,DC},
    laplace::CTSEMLaplaceObjective, U::Int, values::AbstractVector{Float64},
    Ls::Vector{<:AbstractMatrix}, base::AbstractVector{Float64}, columns::Vector{Int},
    members::Vector{Int}, mine, cross) where {W,DC}
    spec, units = laplace.spec, laplace.units
    objective = laplace.objective
    d = length(base)
    uw = ctx.uw
    length(uw) == d || resize!(uw, d)
    for start in 1:W:length(columns)
        cols = start:min(start + W - 1, length(columns))
        nl = length(cols)
        @inbounds for c in 1:d
            uw[c] = DC(base[c])
        end
        @inbounds for (a, t) in enumerate(cols)
            uw[columns[t]] = DC(base[columns[t]],
                ForwardDiff.Partials(ntuple(w -> w == a ? 1.0 : 0.0, Val(W))))
        end
        for p in mine
            m = members[p]
            offsets = units.offsets[U][m]
            xw = _laplace_member_values!(ctx.xw, values, spec, Ls, uw, offsets)
            # A member none of the chunk's coordinates reach has no columns;
            # its cross curvature's `M` and mode are still read.
            cross === nothing && !_tan_any_lanes(xw) && continue
            i = units.members[U][m]
            st = _tan_laplace_member!(ctx, objective, i, xw, nl)
            st == _TAN_OK || return st
            cross === nothing || _tan_laplace_cross!(cross, m, ctx.sw.hs, ctx.mod,
                objective.subject_objectives[i].data, columns, cols)
            @inbounds for (a, t) in enumerate(cols)
                _laplace_loading_chain!(view(J, :, t), view(ctx.Hb, :, a), spec, Ls,
                    offsets)
            end
        end
    end
    return _TAN_OK
end

"""
Member `m`'s cross curvature with the deviations its filter integrates, into
`cross` as `_hybrid_after!` reads it off the dual's pass: `M` and the
conditional mode `dhat` at their values, and `H = -d(b - M dhat)/du` along the
chunk's coordinates, `dhat` held fixed. In the filter's coordinates
(`_hybrid_factor!`): an own variance's items in order, those a manifest
intercept's subject never observes left out, where the sweep keeps every
item's (`tangent_hybrid.jl`); a block's columns, every one.
"""
function _tan_laplace_cross!(cross::CTSEMHybridCross{Float64}, m::Int, hs::TanHybState{W},
    mod::TanModel{W}, data::AbstractMatrix, columns::Vector{Int}, cols) where {W}
    h = mod.hyb
    nh = mod.nh
    seen(k) = h.kind[k] != 1 ||
        any(c -> _ctsem_observed(data[h.target[k], c]), axes(data, 2))
    nobs = nh == 0 ? 0 : count(seen, 1:nh)
    q = nobs == 0 ? 0 : (h.block ? nh : nobs)
    cross.q[m] = q
    q == 0 && return nothing
    if size(cross.H[m], 1) < q
        cross.H[m] = zeros(nh, cross.d)
        cross.A[m] = zeros(nh, nh)
        cross.dhat[m] = zeros(nh)
    end
    A, dhat, H = cross.A[m], cross.dhat[m], cross.H[m]
    M, b, beta = hs.M, hs.b, hs.beta
    keep(k) = h.block || seen(k)
    ci = 0
    @inbounds for k in 1:nh
        keep(k) || continue
        ci += 1
        dhat[ci] = beta.v[k]
        cj = 0
        for k2 in 1:nh
            keep(k2) || continue
            cj += 1
            A[ci, cj] = M.v[k, k2]
        end
        for (a, t) in enumerate(cols)
            r = b.d[a, k]
            for k2 in 1:nh
                r -= M.d[a, k, k2] * beta.v[k2]
            end
            H[ci, columns[t]] = -r
        end
    end
    return nothing
end

"""Whether any of `xw`'s partials is nonzero."""
function _tan_any_lanes(xw::Vector{DC}) where {DC}
    @inbounds for x in xw
        iszero(ForwardDiff.partials(x)) || return true
    end
    return false
end

"""
Subject `i`'s curvature columns at its shifted values `xw` into `ctx.Hb`, the
raw values' gradient's tangents along the chunk's `nl` lanes, unweighted as
the Laplace layer's members are, and the gradient itself into `g0` where
given. The hybrid factor's cotangent is left in `ctx.acc.Gbar`. A status code.
"""
function _tan_laplace_member!(ctx::TanLaplaceContext{W,DC}, objective::CTSEMObjective,
    i::Int, xw::Vector, nl::Int; g0=nothing) where {W,DC}
    sp = objective.params
    so = objective.subject_objectives[i]
    mod, sw, Hb = ctx.mod, ctx.sw, ctx.Hb
    _tan_layer_fill!(mod, ctx.awsW, sp, xw, nl, so.tipreds) || return _TAN_REFUSED
    fill!(Hb, 0.0)
    ctx.rows === nothing ||
        return _tan_sweep_rows!(ctx.acc, Hb, sw, mod, ctx.rows, objective, i, 1.0, xw, nl;
            g0=g0)
    _tan_sweep_begin!(sw, mod)
    sa = _tan_acc_zero!(ctx.acc)
    st = _tan_sweep!(sa, sw, mod, so.timesteps, so.data, so.tdpreds, so.max_timestep, 1.0)
    st == _TAN_OK || return st
    _tan_flush_frechet!(sa, sw.fw, mod, sw.sT)
    # The layer's reverse (`_tan_layer_columns!`), its values kept too.
    C = NamedTuple{_TAN_LAYER_NAMES}(map(k -> _tan_duals(_tan_layer_field(sa, k), DC),
        _TAN_LAYER_NAMES))
    g = _tan_layer_reverse(ctx.awsW, sp, xw, C, so.tipreds)
    @inbounds for w in 1:nl, j in eachindex(g)
        Hb[j, w] = ForwardDiff.partials(g[j], w)
    end
    if g0 !== nothing
        @inbounds for j in eachindex(g)
            g0[j] = ForwardDiff.value(g[j])
        end
    end
    return _TAN_OK
end

"""
`_laplace_unit_seeded_gradient`'s sweeps by the hand-coded second order: each
member's gradient in its shifted values (`d0`) and its derivative along the
raw direction `d1` (`d1c`), one lane of the curvature's context; the members
on the pool as the dual's are. `nothing` where the sweep does not serve, and
the caller takes the dual.
"""
function _tan_laplace_seeded(laplace::CTSEMLaplaceObjective, U::Integer,
    values::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64}, members,
    d1::Vector{Float64})
    _CTSEM_TANGENT_HESSIAN[] || return nothing
    ctx = _tan_laplace_context(laplace, values, _HYBRID_PROFILE)
    ctx === nothing && return nothing
    return _tan_laplace_seeded!(ctx, laplace, Int(U), values, Ls, u, members, d1)
end

function _tan_laplace_seeded!(ctx::TanLaplaceContext{W,DC}, laplace::CTSEMLaplaceObjective,
    U::Int, values::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    members, d1::Vector{Float64}) where {W,DC}
    spec, units = laplace.spec, laplace.units
    objective = laplace.objective
    unitmembers = units.members[U]
    npar = length(values)
    nm = length(members)
    h = objective.params.hybrid
    hyb = ctx.mod.nh > 0
    # A block's factor's derivatives in its level's parameters, once a call.
    dG = hyb && h.block ? _tan_hybrid_block_dG(h, values,
        ForwardDiff.Dual{_tan_tag(_ctsem_hessian_dual(objective)),Float64,1}) : nothing
    hyb && h.block && dG === nothing &&
        return (ok=false, d0=_LAPLACE_NO_MATRIX, d1c=_LAPLACE_NO_MATRIX)
    d0 = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d0)
    d1c = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d1c)
    refused = Ref(false)
    one_member! = function (c::Int, m)
        local cx, shifted, xw, st, g0, g1
        cx = _tan_laplace_context(laplace, values, _HYBRID_PROFILE)::typeof(ctx)
        _tan_laplace_hybrid!(cx, values) || return false
        shifted = _laplace_member_values!(
            _laplace_scratch_vector!(laplace, Float64, npar, :sweep_shift),
            values, spec, Ls, u, units.offsets[U][m])
        xw = cx.xw
        @inbounds for t in 1:npar
            xw[t] = DC(shifted[t], ForwardDiff.Partials(ntuple(w -> w == 1 ? d1[t] : 0.0,
                Val(W))))
        end
        g0 = view(d0, :, c)
        st = _tan_laplace_member!(cx, objective, unitmembers[m], xw, 1; g0=g0)
        st == _TAN_REFUSED && (refused[] = true; return false)
        st == _TAN_OK || return false
        g1 = view(d1c, :, c)
        @inbounds for t in 1:npar
            g1[t] = cx.Hb[t, 1]
        end
        hyb && _tan_laplace_gbar!(g0, g1, cx.acc.Gbar, h, values, dG)
        return true
    end
    ok = _laplace_run_members(one_member!, laplace, members)
    refused[] && return nothing
    ok || return (ok=false, d0=_LAPLACE_NO_MATRIX, d1c=_LAPLACE_NO_MATRIX)
    return (ok=true, d0=d0, d1c=d1c)
end

"""
The hybrid factor's cotangent `Gbar` (with its tangent along lane 1) pushed
onto its raw positions, as `_hybrid_scatter!` does: the gradient into `g0` and
its derivative along a direction that does not move the factor into `g1`, so
by G's first derivative alone -- an own sd's in closed form, a block's from
`dG` (`_tan_hybrid_block_dG`).
"""
function _tan_laplace_gbar!(g0::AbstractVector{Float64}, g1::AbstractVector{Float64},
    Gbar::TanMat, h::CTSEMHybrid, values::AbstractVector{Float64}, dG)
    nh = _hybrid_items(h)
    if !h.block
        @inbounds for k in 1:nh
            idx = h.psi_index[k]
            _, dg, _ = _tan_hybrid_own(h, values, k)
            g0[idx] += Gbar.v[k, k] * dg
            g1[idx] += Gbar.d[1, k, k] * dg
        end
        return nothing
    end
    @inbounds for j in eachindex(h.lparams)
        a0 = 0.0
        a1 = 0.0
        for c in 1:nh, f in 1:c
            a0 += Gbar.v[c, f] * dG[c, f, j]
            a1 += Gbar.d[1, c, f] * dG[c, f, j]
        end
        g0[h.lparams[j]] += a0
        g1[h.lparams[j]] += a1
    end
    return nothing
end

"""
A block's factor's first derivatives in its level's parameters at `values`,
`dG[:, :, j]` along `h.lparams[j]`, from `_hybrid_block_factor!` at the
width-1 dual `D1`; `nothing` where the block does not factor.
"""
function _tan_hybrid_block_dG(h::CTSEMHybrid, values::AbstractVector{Float64},
    ::Type{D1}) where {D1}
    nh = _hybrid_items(h)
    k = length(h.perm)
    dG = zeros(nh, nh, length(h.lparams))
    vd = Vector{D1}(undef, length(values))
    Lb = zeros(D1, k, k)
    for (j, p) in enumerate(h.lparams)
        @inbounds for i in eachindex(values)
            vd[i] = D1(values[i], ForwardDiff.Partials((i == p ? 1.0 : 0.0,)))
        end
        _hybrid_block_factor!(Lb, h, vd) || return nothing
        for c in 1:nh, f in 1:c
            dG[c, f, j] = ForwardDiff.partials(_tan_hybrid_block_G(h, Lb, c, f), 1)
        end
    end
    return dG
end
