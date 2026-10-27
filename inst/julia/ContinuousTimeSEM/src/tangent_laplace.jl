# A Laplace unit's curvature (`_laplace_columns_jacobian`) by the hand-coded
# second order: each member's sweep carries the unit's coordinates as tangents
# (`tangent_sweep.jl`, `tangent_rows.jl`), through the member's shifted values
# at the curvature's dual (`_laplace_curvature_dual`), each coordinate of a
# chunk a lane, where the other route runs the member's adjoint gradient at
# that dual. A coordinate moves only the cells its effects reach, so a random
# intercept's sweep carries tangents on the mean alone, which the dual carries
# as zeros through every covariance. The gated floor's cross curvature with the
# deviations the filter integrates (`CTSEMHybridCross`) is read off the same
# sweeps, at their profile or at the fixed values of its joint coordinates
# (`_HybridFixedCross`). The dual route serves the other hooks, an element type
# other than Float64, and what the sweep does not cover (`_tan_refusal`); it is
# the oracle (`ctsem_set_tangent_hessian!(false)`).
#
# A member whose coordinates reach only the means -- the initial state, the
# intercepts, the TD effects -- of a linear model needs no reverse sweep: its
# innovations are affine in the coordinates and their covariances do not move
# with them, so its curvature is exactly the Gauss-Newton sum over rows of the
# whitened innovations' tangents (`_tan_laplace_gn_member!`), from the forward
# sweep alone.

"""A slot's hand-coded curvature: the member's layer at the curvature's dual,
its sweep, and their scratch."""
struct TanLaplaceContext{W,DC,A,R,B}
    awsW::A
    mod::TanModel{W}
    acc::TanAcc{W}
    rows::R
    sw::TanSweep{W}
    Hb::Matrix{Float64}         # the member's raw columns, nvalues x W
    xw::Vector{DC}              # its shifted values
    uw::Vector{DC}              # the unit's coordinates, seeded
    Lb::Matrix{Float64}         # a hybrid block's factor
    # The Gauss-Newton route (`_tan_laplace_gn_member!`): whether the model
    # admits it at all -- linear, no integrated deviations, every state in the
    # dynamics -- each chunk's whitened innovation tangents, a row's W at a
    # time, the columns' positions by coordinate, a member's coordinates, and
    # the members left to the full sweeps.
    gn::Bool
    gnE::Vector{Vector{Float64}}
    colpos::Vector{Int}
    coords::Vector{Int}
    rest::Vector{Int}
    buf::B                      # the member layer's (`_tan_layer_buffers`)
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
    # Fixed deviations are the joint curvature's, which this route serves where
    # it serves their profile; the Hessian's refusal of them is its own.
    ctx = _tan_refusal(objective, mode == _HYBRID_FIXED ? _HYBRID_PROFILE : mode) === nothing ?
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
        wlog=mode == _HYBRID_FULL ? 1.0 : 0.0, fixed=mode == _HYBRID_FIXED && nh > 0)
    awsW = CTSEMAdjointWorkspace(DC, sp, nv, ntd)
    return TanLaplaceContext(awsW, mod, TanAcc(n, m, ntd, nh, W),
        rows ? TanRows(awsW, sp, n, m, ntd, Val(W), DC) : nothing,
        TanSweep(n, m, ntd, nh, W, na, length(dyn)), zeros(nv, W), zeros(DC, nv),
        DC[], h.block && nh > 0 ? zeros(length(h.perm), length(h.perm)) : zeros(0, 0),
        !rows && nh == 0 && !mod.aug, Vector{Float64}[], Int[], Int[], Int[],
        _tan_layer_buffers(DC, awsW, sp, nv, ntd))
end

"""
`Jw`'s first `length(columns)` columns, the Jacobian of `members`' log
likelihood gradient in the unit's coordinates `columns` at `base`, by the
hand-coded sweeps, the rest zero; `false` where they do not serve, and the
caller takes the dual. A member whose filter fails makes it `NaN`, as the
dual's gradient does. `cross`, a hook where given (`CTSEMHybridCross`, or
`_HybridFixedCross` at fixed deviations), is read and filled as the dual's
passes read and fill it (`_hybrid_before!`, `_hybrid_after!`).
"""
function _tan_laplace_columns_jacobian!(Jw::Matrix{Float64}, laplace::CTSEMLaplaceObjective,
    U::Integer, values::AbstractVector{Float64}, Ls::Vector{<:AbstractMatrix},
    base::AbstractVector{Float64}, columns::Vector{Int}, members::Vector{Int}, mode::Int,
    cross=nothing)
    (mode == _HYBRID_FIXED) == (cross isa _HybridFixedCross) || return false
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
    # The members the Gauss-Newton route serves first; the rest swept below.
    if cross === nothing && ctx.gn
        colpos = ctx.colpos
        length(colpos) < d && append!(colpos, zeros(Int, d - length(colpos)))
        @inbounds for (t, c) in enumerate(columns)
            colpos[c] = t
        end
        rest = empty!(ctx.rest)
        st = _TAN_OK
        for p in mine
            st = _tan_laplace_gn_member!(J, ctx, laplace, U, values, Ls, base, members[p])
            st == _TAN_FULL && (push!(rest, p); st = _TAN_OK)
            st == _TAN_OK || break
        end
        @inbounds for c in columns
            colpos[c] = 0
        end
        st == _TAN_OK || return st
        mine = rest
        isempty(mine) && return _TAN_OK
    end
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
            data = objective.subject_objectives[i].data
            cross isa _HybridFixedCross && _tan_laplace_dfix!(ctx.mod, cross.joint, m, data)
            st = _tan_laplace_member!(ctx, objective, i, xw, nl)
            st == _TAN_OK || return st
            cross === nothing || _tan_laplace_cross!(cross isa _HybridFixedCross ?
                cross.cross : cross, m, ctx.sw.hs, ctx.mod, data, columns, cols)
            @inbounds for (a, t) in enumerate(cols)
                _laplace_loading_chain!(view(J, :, t), view(ctx.Hb, :, a), spec, Ls,
                    offsets)
            end
        end
    end
    return _TAN_OK
end

"""
Whether a subject observes hybrid item `k` in the filter's sense
(`_hybrid_start!`): any item but a manifest intercept, and that where any row
observes its manifest.
"""
_tan_hybrid_seen(h::CTSEMHybrid, k::Int, data::AbstractMatrix) = h.kind[k] != 1 ||
    any(c -> _ctsem_observed(data[h.target[k], c]), axes(data, 2))

"""
Member `m`'s fixed deviations from the joint coordinates (`_HybridJointHook`,
as `_hybrid_before!` hands them to the filter) into the model's `dfix`, one
per item: the filter's coordinates mapped as `_tan_laplace_cross!` maps them,
an item the subject does not observe at zero.
"""
function _tan_laplace_dfix!(mod::TanModel, joint, m::Int, data::AbstractMatrix)
    h = mod.hyb
    dfix = mod.dfix
    fill!(dfix.v, 0.0)
    _tan_setnz!(dfix, false)
    q = joint.layout.q[m]
    q == 0 && return nothing
    off = joint.layout.offs[m]
    c = 0
    @inbounds for k in 1:mod.nh
        (h.block || _tan_hybrid_seen(h, k, data)) || continue
        c += 1
        c <= q && (dfix.v[k] = joint.z[off + c])
    end
    return nothing
end

"""
Member `m`'s cross curvature with the deviations its filter integrates, into
`cross` as `_hybrid_after!` reads it off the dual's pass: `M` and the
conditional mode `dhat` at their values, and `H = -d(b - M d)/du` along the
chunk's coordinates at the deviations `d` the pass is at -- the conditional
mode, or the fixed ones -- held fixed. In the filter's coordinates
(`_hybrid_factor!`): an own variance's items in order, those a manifest
intercept's subject never observes left out, where the sweep keeps every
item's (`tangent_hybrid.jl`); a block's columns, every one.
"""
function _tan_laplace_cross!(cross::CTSEMHybridCross{Float64}, m::Int, hs::TanHybState{W},
    mod::TanModel{W}, data::AbstractMatrix, columns::Vector{Int}, cols) where {W}
    h = mod.hyb
    nh = mod.nh
    seen(k) = _tan_hybrid_seen(h, k, data)
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
    dref = mod.fixed ? mod.dfix : beta
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
                r -= M.d[a, k, k2] * dref.v[k2]
            end
            H[ci, columns[t]] = -r
        end
    end
    return nothing
end

"""
What the Kalman trace (`CTSEMKalmanTrace`) keeps of a hand-coded forward sweep
(`_tan_sweep!`, `_tan_sweep_rows!` with `rec`): at each row's prediction (kind
1) and update (2) the state and its covariance as duals of the lanes, recorded
by the engine's own (`_record_row_prior!`, `_record_row_update!`) at the cells
the row's filter holds -- the layer's, or by row `TanRows.allW` -- so the trace
is the filter's at the dual, response scale and all, without the filter there.
"""
struct TanTraceRecorder{TR,WS,V,SO,SP}
    trace::TR
    ws::WS
    cells::V
    so::SO
    sp::SP
    i::Int
    tp::Base.RefValue{AbstractVector}   # the TI predictors the cells read
end

function _tan_record!(rec::TanTraceRecorder, kind::Int, t::Int, x::TanVec, P::TanMat,
    sw::TanSweep{W}) where {W}
    ws, sp = rec.ws, rec.sp
    n = length(x.v)
    _tan_to_duals!(view(ws.state, 1:n), x)
    _tan_to_duals!(kind == 1 ? ws.P_predict.data : ws.P_update.data, P)
    pars = CTSEMPars(rec.cells, sp.parameter_axis)
    if kind == 1
        sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, ws.sp.covmatcode, ws.manifest_dim)
        _record_row_prior!(rec.trace, ws, pars, t)
        return nothing
    end
    DC = eltype(ws.state)
    o = (t - 1) * (1 + W)
    ll = DC(sw.llrow[o + 1], ForwardDiff.Partials(ntuple(a -> sw.llrow[o + 1 + a], Val(W))))
    times = rec.so.timesteps
    ctx = CTSEMRowContext(ws.state, pars, view(rec.so.tdpreds, :, t), rec.tp[], times[t],
        t == 1 ? 0.0 : times[t] - times[t - 1], rec.i, t)
    _record_row_update!(rec.trace, ws, pars, sp, ctx, t, ll)
    return nothing
end

"""
Subject `i`'s plug-in Kalman trace at the duals `xw` (the curvature's, their
first `nl` lanes seeded), by the hand-coded forward sweep, into `trace`, as
`_extended_kalman_filter_continuous!` at that dual writes it with the
deviations the filter integrates off: `true` where the pass came out finite,
`false` where not, `nothing` where the sweep does not serve and the caller
filters at the dual.
"""
function _tan_laplace_trace!(trace::CTSEMKalmanTrace, laplace::CTSEMLaplaceObjective,
    values::AbstractVector{Float64}, i::Int, xw::Vector, nl::Int)
    _CTSEM_TANGENT_HESSIAN[] || return nothing
    ctx = _tan_laplace_context(laplace, values, _HYBRID_OFF)
    ctx === nothing && return nothing
    st = _tan_laplace_trace_sweep!(trace, ctx, laplace.objective, i, xw, nl)::Int
    st == _TAN_REFUSED && return nothing
    return st == _TAN_OK
end

function _tan_laplace_trace_sweep!(trace::CTSEMKalmanTrace, ctx::TanLaplaceContext{W,DC},
    objective::CTSEMObjective, i::Int, xw::Vector, nl::Int) where {W,DC}
    eltype(xw) === DC || return _TAN_REFUSED
    sp = objective.params
    so = objective.subject_objectives[i]
    mod, sw = ctx.mod, ctx.sw
    _tan_layer_fill!(mod, ctx.awsW, sp, xw, nl, so.tipreds) || return _TAN_REFUSED
    ws = ctx.awsW.ekf_ws
    if ctx.rows === nothing
        rec = TanTraceRecorder(trace, ws, ws.all_params, so, sp, i,
            Ref{AbstractVector}(_ctsem_tipred_vector(so.tipreds, xw)))
        _tan_sweep_begin!(sw, mod)
        return _tan_sweep!(ctx.acc, sw, mod, so.timesteps, so.data, so.tdpreds,
            so.max_timestep, 1.0; forward=true, rec=rec)
    end
    R = ctx.rows
    rec = TanTraceRecorder(trace, ws, R.allW, so, sp, i, R.tp)
    return _tan_sweep_rows!(ctx.acc, zeros(0, 0), sw, mod, R, objective, i, 1.0, xw, nl;
        forward=true, rec=rec)
end

# A member the Gauss-Newton route leaves to the full sweeps.
const _TAN_FULL = -1

"""
Whether the layer just filled into `mod` moves only the means -- the initial
state, the intercepts, the TD effects -- along its lanes, in a model whose
every other matrix it leaves fixed: the innovations then are affine in the
lanes and their covariances are not moved by them.
"""
_tan_gn_quadratic(mod::TanModel) = !(_tan_nz(mod.JAx) || _tan_nz(mod.DRIFT) ||
    _tan_nz(mod.Qc) || _tan_nz(mod.LAMBDA) || _tan_nz(mod.Jy) || _tan_nz(mod.R) ||
    _tan_nz(mod.T0V) || _tan_nz(mod.Jtd))

"""Each updated row's whitened innovation tangents, `W` a row, into `E`."""
struct TanGNRecorder
    E::Vector{Float64}
end

function _tan_record!(rec::TanGNRecorder, kind::Int, t::Int, x::TanVec, P::TanMat,
    sw::TanSweep{W}) where {W}
    kind == 2 || return nothing
    mo = length(sw.upds[t].o)
    mo == 0 && return nothing
    e = _tan_row(sw.fw, mo).e
    @inbounds for i in 1:mo, a in 1:W
        push!(rec.E, e.d[a, i])
    end
    return nothing
end

"""
Unit member `m`'s curvature over the columns at `ctx.colpos` into `J`, by
Gauss-Newton where its coordinates move only the means (`_tan_gn_quadratic`):
a row's log likelihood is `-1/2 (e'e + logdet S)` with `e = Ui' (y - Lambda x
- mu)` affine in them and `S` fixed, so its second derivative is `-de' de`
exactly, from the forward sweeps alone, a chunk of the member's coordinates
each. `_TAN_FULL` where they move anything else, and the caller sweeps the
member in full. A status code otherwise.
"""
function _tan_laplace_gn_member!(J::Matrix{Float64}, ctx::TanLaplaceContext{W,DC},
    laplace::CTSEMLaplaceObjective, U::Int, values::AbstractVector{Float64},
    Ls::Vector{<:AbstractMatrix}, base::AbstractVector{Float64}, m::Int) where {W,DC}
    spec, units = laplace.spec, laplace.units
    objective = laplace.objective
    offsets = units.offsets[U][m]
    i = units.members[U][m]
    so = objective.subject_objectives[i]
    sp = objective.params
    coords = empty!(ctx.coords)
    for (l, level) in enumerate(spec.levels)
        (nrandomeffects(level) == 0 || nlatent(level) == 0) && continue
        for q in 1:nlatent(level)
            push!(coords, offsets[l] + q)
        end
    end
    # Nothing it can add to these columns.
    any(c -> ctx.colpos[c] != 0, coords) || return _TAN_OK
    nc = length(coords)
    nch = cld(nc, W)
    while length(ctx.gnE) < nch
        push!(ctx.gnE, Float64[])
    end
    uw, mod, sw = ctx.uw, ctx.mod, ctx.sw
    for ch in 1:nch
        lanes = ((ch - 1) * W + 1):min(ch * W, nc)
        @inbounds for c in coords
            uw[c] = DC(base[c])
        end
        @inbounds for (a, k) in enumerate(lanes)
            c = coords[k]
            uw[c] = DC(base[c], ForwardDiff.Partials(ntuple(w -> w == a ? 1.0 : 0.0, Val(W))))
        end
        xw = _laplace_member_values!(ctx.xw, values, spec, Ls, uw, offsets)
        _tan_layer_fill!(mod, ctx.awsW, sp, xw, length(lanes), so.tipreds, ctx.buf) ||
            return _TAN_REFUSED
        _tan_gn_quadratic(mod) || return _TAN_FULL
        rec = TanGNRecorder(empty!(ctx.gnE[ch]))
        _tan_sweep_begin!(sw, mod)
        st = _tan_sweep!(ctx.acc, sw, mod, so.timesteps, so.data, so.tdpreds,
            so.max_timestep, 1.0; forward=true, rec=rec)
        st == _TAN_OK || return st
    end
    nr = length(ctx.gnE[1]) ÷ W
    @inbounds for b in 1:nc
        t = ctx.colpos[coords[b]]
        t == 0 && continue
        Eb = ctx.gnE[(b - 1) ÷ W + 1]
        lb = (b - 1) % W + 1
        for a in 1:nc
            Ea = ctx.gnE[(a - 1) ÷ W + 1]
            la = (a - 1) % W + 1
            s = 0.0
            for r in 0:(nr - 1)
                s += Ea[r * W + la] * Eb[r * W + lb]
            end
            J[coords[a], t] -= s
        end
    end
    return _TAN_OK
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
    rows = ctx.rows !== nothing
    _tan_layer_fill!(mod, ctx.awsW, sp, xw, nl, so.tipreds, rows ? nothing : ctx.buf) ||
        return _TAN_REFUSED
    fill!(Hb, 0.0)
    rows && return _tan_sweep_rows!(ctx.acc, Hb, sw, mod, ctx.rows, objective, i, 1.0, xw,
        nl; g0=g0)
    _tan_sweep_begin!(sw, mod)
    sa = _tan_acc_zero!(ctx.acc)
    st = _tan_sweep!(sa, sw, mod, so.timesteps, so.data, so.tdpreds, so.max_timestep, 1.0)
    st == _TAN_OK || return st
    _tan_flush_frechet!(sa, sw.fw, mod, sw.sT)
    # The layer's reverse (`_tan_layer_columns!`), its values kept too, at
    # the cells the fill left and in the context's buffers.
    buf = ctx.buf
    for k in _TAN_LAYER_NAMES
        _tan_to_duals!(getproperty(buf.C, k), _tan_layer_field(sa, k))
    end
    g = _tan_layer_reverse(ctx.awsW, sp, xw, buf.C, so.tipreds, buf)
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
        g1 = view(d1c, :, c:c)
        @inbounds for t in 1:npar
            g1[t, 1] = cx.Hb[t, 1]
        end
        hyb && _tan_laplace_gbar!(g0, g1, cx.acc.Gbar, h, values, dG, 1)
        return true
    end
    ok = _laplace_run_members(one_member!, laplace, members)
    refused[] && return nothing
    ok || return (ok=false, d0=_LAPLACE_NO_MATRIX, d1c=_LAPLACE_NO_MATRIX)
    return (ok=true, d0=d0, d1c=d1c)
end

"""
The hybrid factor's cotangent `Gbar` (with its tangents along its first `nl`
lanes) pushed onto its raw positions, as `_hybrid_scatter!` does: the gradient
into `g0` and its derivatives along directions that do not move the factor
into `H`'s columns, so by G's first derivative alone -- an own sd's in closed
form, a block's from `dG` (`_tan_hybrid_block_dG`).
"""
function _tan_laplace_gbar!(g0::AbstractVector{Float64}, H::AbstractMatrix{Float64},
    Gbar::TanMat, h::CTSEMHybrid, values::AbstractVector{Float64}, dG, nl::Int)
    nh = _hybrid_items(h)
    if !h.block
        @inbounds for k in 1:nh
            idx = h.psi_index[k]
            _, dg, _ = _tan_hybrid_own(h, values, k)
            g0[idx] += Gbar.v[k, k] * dg
            for a in 1:nl
                H[idx, a] += Gbar.d[a, k, k] * dg
            end
        end
        return nothing
    end
    @inbounds for j in eachindex(h.lparams)
        a0 = 0.0
        for c in 1:nh, f in 1:c
            a0 += Gbar.v[c, f] * dG[c, f, j]
        end
        g0[h.lparams[j]] += a0
        for a in 1:nl
            a1 = 0.0
            for c in 1:nh, f in 1:c
                a1 += Gbar.d[a, c, f] * dG[c, f, j]
            end
            H[h.lparams[j], a] += a1
        end
    end
    return nothing
end

"""
`_laplace_unit_theta_gradient` at `z` and its derivatives along the first `nl`
columns of `B` (directions in `z`), by the hand-coded sweeps, into `Gt`
(`npar x (1 + nl)`) and `Gz` (`length(z) x (1 + nl)`), the value first: what
the gated rule's first-order pass and the mode's cross derivatives seeded a
dual in `z` for. In joint coordinates (`layout`) each member's filter is at the
deviations `z` fixes, their directions among the lanes. `nothing` where the
sweep does not serve, `false` where a member cannot be evaluated.
"""
function _tan_laplace_theta_gradient!(Gt::Matrix{Float64}, Gz::Matrix{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dL, positions, z::Vector{Float64},
    B::AbstractMatrix{Float64}, nl::Int, layout)
    _CTSEM_TANGENT_HESSIAN[] || return nothing
    mode = layout === nothing ? _HYBRID_PROFILE : _HYBRID_FIXED
    ctx = _tan_laplace_context(laplace, theta, mode)
    ctx === nothing && return nothing
    st = _tan_laplace_theta_gradient_!(Gt, Gz, ctx, laplace, Int(U), theta, Ls, dL,
        positions, z, B, nl, layout)::Int
    st == _TAN_REFUSED && return nothing
    return st == _TAN_OK
end

function _tan_laplace_theta_gradient_!(Gt::Matrix{Float64}, Gz::Matrix{Float64},
    ctx::TanLaplaceContext{W,DC}, laplace::CTSEMLaplaceObjective, U::Int,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, dL, positions,
    z::Vector{Float64}, B::AbstractMatrix{Float64}, nl::Int, layout) where {W,DC}
    spec, units = laplace.spec, laplace.units
    objective = laplace.objective
    members = units.members[U]
    npar = length(theta)
    D = length(z)
    d = layout === nothing ? D : layout.d
    fill!(Gt, 0.0)
    fill!(Gz, 0.0)
    _tan_laplace_hybrid!(ctx, theta) || return _TAN_INVALID
    mod = ctx.mod
    h = objective.params.hybrid
    hyb = mod.nh > 0
    dG = hyb && h.block ? _tan_hybrid_block_dG(h, theta,
        ForwardDiff.Dual{_tan_tag(_ctsem_hessian_dual(objective)),Float64,1}) : nothing
    hyb && h.block && dG === nothing && return _TAN_INVALID
    uw = ctx.uw
    length(uw) == d || resize!(uw, d)
    @inbounds for a in 1:d
        uw[a] = DC(z[a], ForwardDiff.Partials(ntuple(w -> w <= nl ? B[a, w] : 0.0, Val(W))))
    end
    u = view(z, 1:d)
    g0 = view(Gt, :, 1)
    gm = zeros(npar)
    for m in eachindex(members)
        offsets = units.offsets[U][m]
        i = members[m]
        xw = _laplace_member_values!(ctx.xw, theta, spec, Ls, uw, offsets)
        layout === nothing || _tan_laplace_dfix_lanes!(mod, layout, m, z, B, nl,
            objective.subject_objectives[i].data)
        st = _tan_laplace_member!(ctx, objective, i, xw, nl; g0=gm)
        st == _TAN_OK || return st
        Hb = ctx.Hb
        hyb && _tan_laplace_gbar!(gm, Hb, ctx.acc.Gbar, h, theta, dG, nl)
        @inbounds for t in 1:npar
            g0[t] += gm[t]
        end
        @inbounds for a in 1:nl, t in 1:npar
            Gt[t, 1 + a] += Hb[t, a]
        end
        # Through `L z`: the population parameters' share, and `z`'s.
        _laplace_chol_chain!(g0, gm, spec, dL, positions, u, offsets)
        _laplace_loading_chain!(view(Gz, :, 1), gm, spec, Ls, offsets)
        for a in 1:nl
            _laplace_chol_chain!(view(Gt, :, 1 + a), view(Hb, :, a), spec, dL, positions,
                u, offsets)
            _laplace_chol_chain!(view(Gt, :, 1 + a), gm, spec, dL, positions,
                view(B, 1:d, a), offsets)
            _laplace_loading_chain!(view(Gz, :, 1 + a), view(Hb, :, a), spec, Ls, offsets)
        end
        layout === nothing || _tan_laplace_fixed_gradient!(Gz, ctx.sw.hs, mod, layout, m,
            nl, objective.subject_objectives[i].data)
    end
    @inbounds for a in 1:D
        Gz[a, 1] -= z[a]
        for w in 1:nl
            Gz[a, 1 + w] -= B[a, w]
        end
    end
    return _TAN_OK
end

"""
Subject `i`'s log likelihood at `shifted`, its gradient `g` and that gradient's
derivatives along the raw positions `active` (`Hs`, rows and columns in
`active`'s order), by the hand-coded sweeps a chunk of the curvature's lanes a
pass: what `_continuation_member_hessian!` seeds a dual for, as
`(loglik, status)`, status 1 where the filter fails. `nothing` where the sweep
does not serve, and where the filter integrates deviations, whose factor the
population parameters among `active` would move.
"""
function _tan_laplace_member_hessian!(Hs::AbstractMatrix{Float64}, g::Vector{Float64},
    laplace::CTSEMLaplaceObjective, i::Integer, theta::Vector{Float64},
    shifted::AbstractVector{Float64}, active::Vector{Int})
    _CTSEM_TANGENT_HESSIAN[] || return nothing
    ctx = _tan_laplace_context(laplace, theta, _HYBRID_FULL)
    (ctx === nothing || ctx.mod.nh > 0) && return nothing
    return _tan_laplace_member_hessian_!(Hs, g, ctx, laplace.objective, Int(i), shifted,
        active)
end

"""Whether `_tan_laplace_member_hessian!` serves `laplace`'s members at `theta`."""
function _tan_laplace_serves_member_hessian(laplace::CTSEMLaplaceObjective,
    theta::AbstractVector)
    _CTSEM_TANGENT_HESSIAN[] || return false
    ctx = _tan_laplace_context(laplace, theta, _HYBRID_FULL)
    return ctx !== nothing && ctx.mod.nh == 0
end

function _tan_laplace_member_hessian_!(Hs::AbstractMatrix{Float64}, g::Vector{Float64},
    ctx::TanLaplaceContext{W,DC}, objective::CTSEMObjective, i::Int,
    shifted::AbstractVector{Float64}, active::Vector{Int}) where {W,DC}
    na = length(active)
    xw = ctx.xw
    loglik = NaN
    for lo in 0:W:(na - 1)
        ncol = min(W, na - lo)
        @inbounds for t in eachindex(shifted)
            xw[t] = DC(shifted[t])
        end
        @inbounds for c in 1:ncol
            t = active[lo + c]
            xw[t] = DC(shifted[t], ForwardDiff.Partials(ntuple(q -> q == c ? 1.0 : 0.0,
                Val(W))))
        end
        st = _tan_laplace_member!(ctx, objective, i, xw, ncol; g0=lo == 0 ? g : nothing)
        st == _TAN_REFUSED && return nothing
        st == _TAN_OK || return (NaN, 1)
        lo == 0 && (loglik = ctx.sw.ll[1])
        @inbounds for c in 1:ncol, r in 1:na
            Hs[r, lo + c] = ctx.Hb[active[r], c]
        end
    end
    return (loglik, 0)
end

"""
Member `m`'s fixed deviations from joint coordinates `z` into the model's
`dfix`, as `_tan_laplace_dfix!`, with their tangents from the rows of `B`
that lay them out.
"""
function _tan_laplace_dfix_lanes!(mod::TanModel{W}, layout, m::Int, z::Vector{Float64},
    B::AbstractMatrix{Float64}, nl::Int, data::AbstractMatrix) where {W}
    h = mod.hyb
    dfix = mod.dfix
    fill!(dfix.v, 0.0)
    fill!(dfix.d, 0.0)
    q = layout.q[m]
    off = layout.offs[m]
    c = 0
    @inbounds for k in 1:mod.nh
        (h.block || _tan_hybrid_seen(h, k, data)) || continue
        c += 1
        c <= q || continue
        dfix.v[k] = z[off + c]
        for a in 1:nl
            dfix.d[a, k] = B[off + c, a]
        end
    end
    _tan_setnz_from!(dfix)
    return nothing
end

"""
Member `m`'s gradient in its fixed deviations, `b - (M - I) d`
(`_hybrid_fixed_gradient!`), and its tangents, added into `Gz`'s rows that
lay them out, in the filter's coordinates as `_tan_laplace_cross!` maps them.
"""
function _tan_laplace_fixed_gradient!(Gz::Matrix{Float64}, hs::TanHybState{W},
    mod::TanModel{W}, layout, m::Int, nl::Int, data::AbstractMatrix) where {W}
    h = mod.hyb
    nh = mod.nh
    q = layout.q[m]
    q == 0 && return nothing
    off = layout.offs[m]
    M, b, dfix = hs.M, hs.b, mod.dfix
    c = 0
    @inbounds for k in 1:nh
        (h.block || _tan_hybrid_seen(h, k, data)) || continue
        c += 1
        c <= q || continue
        v = b.v[k]
        for k2 in 1:nh
            v -= (M.v[k, k2] - (k == k2 ? 1.0 : 0.0)) * dfix.v[k2]
        end
        Gz[off + c, 1] += v
        for a in 1:nl
            t = b.d[a, k]
            for k2 in 1:nh
                t -= M.d[a, k, k2] * dfix.v[k2] +
                    (M.v[k, k2] - (k == k2 ? 1.0 : 0.0)) * dfix.d[a, k2]
            end
            Gz[off + c, 1 + a] += t
        end
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
