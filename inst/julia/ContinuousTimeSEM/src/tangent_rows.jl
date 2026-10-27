# Models whose cells read the state, a TD predictor or the time, which the
# filter evaluates by row in three groups -- before each prediction substep,
# before a row's TD impulse, before its update (`_record_group!`). The sweep
# keeps the matrices' algebra by hand; the cells are the engine's own code at
# the Hessian's dual, each lane a column: forward by
# `apply_complex_transforms_at_indices!` on the subject's cells, back by
# `_reverse_group!` on their cotangents. A group's evaluation refreshes the
# model's matrices from the cells and keeps a copy, which each reverse step that
# read them restores; a group's reverse first moves the matrices' cotangents
# onto the cells, through the covariance constructions at the cells of the
# stretch it began, then pulls them through its transforms.

"""The model's matrices as a group left them, and its covariances' cells."""
struct TanSnap{W,DW}
    JAx::TanMat{W}; DRIFT::TanMat{W}; DmJ::TanMat{W}; CINT::TanVec{W}; Qc::TanMat{W}
    LAMBDA::TanMat{W}; Jy::TanMat{W}; MM::TanVec{W}; R::TanMat{W}
    TDE::TanMat{W}; Jtd::TanMat{W}
    rawD::Matrix{DW}; rawR::Matrix{DW}
end
TanSnap{W,DW}(n::Int, m::Int, ntd::Int) where {W,DW} = TanSnap{W,DW}(TanMat{W}(n, n),
    TanMat{W}(n, n), TanMat{W}(n, n), TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(m, n),
    TanMat{W}(m, n), TanVec{W}(m), TanMat{W}(m, m), TanMat{W}(n, ntd), TanMat{W}(n, n),
    zeros(DW, n, n), zeros(DW, m, m))

function _tan_snap_save!(s::TanSnap, mod::TanModel, allW::AbstractVector, sp)
    _tan_copy!(s.JAx, mod.JAx); _tan_copy!(s.DRIFT, mod.DRIFT); _tan_copy!(s.DmJ, mod.DmJ)
    _tan_copy!(s.CINT, mod.CINT); _tan_copy!(s.Qc, mod.Qc)
    _tan_copy!(s.LAMBDA, mod.LAMBDA); _tan_copy!(s.Jy, mod.Jy); _tan_copy!(s.MM, mod.MM)
    _tan_copy!(s.R, mod.R); _tan_copy!(s.TDE, mod.TDE); _tan_copy!(s.Jtd, mod.Jtd)
    pars = CTSEMPars(allW, sp.parameter_axis)
    n, m = size(s.rawD, 1), size(s.rawR, 1)
    s.rawD .= view(pars.DIFFUSION, 1:n, 1:n)
    s.rawR .= view(pars.MANIFESTVAR, 1:m, 1:m)
    return s
end

function _tan_snap_load!(mod::TanModel, s::TanSnap)
    _tan_copy!(mod.JAx, s.JAx); _tan_copy!(mod.DRIFT, s.DRIFT); _tan_copy!(mod.DmJ, s.DmJ)
    _tan_copy!(mod.CINT, s.CINT); _tan_copy!(mod.Qc, s.Qc)
    _tan_copy!(mod.LAMBDA, s.LAMBDA); _tan_copy!(mod.Jy, s.Jy); _tan_copy!(mod.MM, s.MM)
    _tan_copy!(mod.R, s.R); _tan_copy!(mod.TDE, s.TDE); _tan_copy!(mod.Jtd, s.Jtd)
    return mod
end

"""
A worker's state for subjects whose cells change by row: its workspace at the
Hessian's dual (whose tape holds the group records), the subject's cells and
their cotangents, the state and its cotangent as duals, the snapshots, the
layer's own (`base`), and the forward's program: `(kind, index, snapshot)`,
kind 1 a group (its tape slot), 2 a prediction substep, 3 a TD impulse, 4 an
update (the row), each with the snapshot it read.
"""
mutable struct TanRows{W,DW,A}
    awsW::A
    allW::Vector{DW}
    allbar::Vector{DW}
    xW::Vector{DW}
    xbarW::Vector{DW}
    snaps::Vector{TanSnap{W,DW}}
    base::TanSnap{W,DW}
    prog::Vector{NTuple{3,Int}}
    ng::Int
end
function TanRows(awsW, sp, n::Int, m::Int, ntd::Int, ::Val{W}, ::Type{DW}) where {W,DW}
    nall = length(sp.mutables)
    return TanRows{W,DW,typeof(awsW)}(awsW, zeros(DW, nall), zeros(DW, nall),
        zeros(DW, n), zeros(DW, n), TanSnap{W,DW}[], TanSnap{W,DW}(n, m, ntd),
        NTuple{3,Int}[], 0)
end

"""`out` the duals of `X`'s values and tangents, lane `w` its `w`-th."""
function _tan_to_duals!(out::AbstractVector{DW}, X::Union{TanMat{W},TanVec{W}}) where {W,DW}
    @inbounds for idx in eachindex(X.v)
        o = W * (idx - 1)
        out[idx] = DW(X.v[idx], ForwardDiff.Partials(ntuple(w -> X.d[o + w], Val(W))))
    end
    return out
end

"""`X += D`, the duals' values and their first `nl` lanes as tangents."""
function _tan_add_lanes!(X::TanVec{W}, D::AbstractVector, nl::Int) where {W}
    @inbounds for i in eachindex(D)
        X.v[i] += ForwardDiff.value(D[i])
        p = ForwardDiff.partials(D[i])
        for w in 1:nl
            iszero(p[w]) && continue
            # Zero tangents are stored zero, so marking them live needs no fill.
            X.nz[] = true
            X.d[w + W * (i - 1)] += p[w]
        end
    end
    return X
end

"""
The model's matrices from the subject's cells `R.allW` as the filter builds
them, their lanes the tangents, and `DmJ`.
"""
function _tan_rows_refresh!(mod::TanModel{W}, R::TanRows, sp, nl::Int) where {W}
    mw = _tan_mats_of(R.awsW, sp, R.allW)
    for k in _TAN_LAYER_NAMES
        _tan_lanes!(_tan_layer_field(mod, k), getproperty(mw, k), nl)
    end
    _tan_comb!(mod.DmJ, ((1.0, mod.DRIFT), (-1.0, mod.JAx)))
    return mod
end

"""
Group `g`'s cells at the state `x` (a context as the filter's for row `t`,
`time`, `dt`): recorded on the workspace's tape, applied to the cells, the
model refreshed and kept. Returns the snapshot the model now holds, `snap`
where the group is empty.
"""
function _tan_rows_group!(R::TanRows{W,DW}, g::Int, x::TanVec{W}, mod::TanModel{W},
    sw::TanSweep{W}, sp, so, subject::Int, t::Int, time::Float64, dt::Float64,
    snap::Int, nl::Int) where {W,DW}
    awsW = R.awsW
    indices = g == 1 ? awsW.predict_indices : g == 2 ? awsW.td_indices : awsW.update_indices
    isempty(indices) && return snap
    _tan_to_duals!(R.xW, x)
    pars = CTSEMPars(R.allW, sp.parameter_axis)
    ctx = CTSEMRowContext(R.xW, pars, view(so.tdpreds, :, t), so.tipreds, time, dt,
        subject, t)
    slot = _record_group_slot!(awsW.tape, g, indices, R.allW, ctx)
    # One call per group: the three transform collections are different types.
    if g == 1
        apply_complex_transforms_at_indices!(R.allW, indices, sp.predict_transforms, ctx)
    elseif g == 2
        apply_complex_transforms_at_indices!(R.allW, indices, sp.td_transforms, ctx)
    else
        apply_complex_transforms_at_indices!(R.allW, indices, sp.update_transforms, ctx)
    end
    _tan_rows_refresh!(mod, R, sp, nl)
    # The drift's factors for the predictions that read it.
    g == 1 && _tan_sweep_begin!(sw, mod)
    R.ng += 1
    R.ng > length(R.snaps) && push!(R.snaps, TanSnap{W,DW}(sw.n, sw.m, mod.ntd))
    _tan_snap_save!(R.snaps[R.ng], mod, R.allW, sp)
    push!(R.prog, (1, slot, R.ng))
    return R.ng
end

"""
The matrices' cotangents in `acc` moved onto the cells' (`R.allbar`), the
covariances through their constructions at `s`'s cells (T0VAR's at the
subject's own, which no group writes), and `acc` emptied.
"""
function _tan_rows_flush!(R::TanRows{W,DW}, acc::TanAcc{W}, sp, s::TanSnap) where {W,DW}
    C = NamedTuple{_TAN_LAYER_NAMES}(map(k -> _tan_duals(_tan_layer_field(acc, k), DW),
        _TAN_LAYER_NAMES))
    pars = CTSEMPars(R.allW, sp.parameter_axis)
    _tan_bar_of!(R.allbar, R.awsW, sp, C, s.rawD, s.rawR, pars.T0VAR)
    _tan_acc_zero!(acc)
    return R
end

"""
One subject, `i` of `objective`, weight `w`, whose cells change by row: the
sweep in the filter's order, its groups' cells evaluated as it goes; the
reverse walk restoring each step's matrices, and each group's reverse moving
the cotangents onto the cells and through it; then the cells' cotangents back
to the raw parameters, lane `w` into `Hb[:, w]`. `R.allW` holds the subject's
cells from the layer at `xw` (`_tan_layer_fill!`), `mod` its matrices.
"""
function _tan_sweep_rows!(acc::TanAcc{W}, Hb::Matrix{Float64}, sw::TanSweep{W},
    mod::TanModel{W}, R::TanRows{W,DW}, objective::CTSEMObjective, i::Int, w::Float64,
    xw::Vector, nl::Int) where {W,DW}
    sp = objective.params
    so = objective.subject_objectives[i]
    awsW = R.awsW
    n, fw, rv = sw.n, sw.fw, sw.rv
    times, Y, TD, rule = so.timesteps, so.data, so.tdpreds, so.max_timestep
    T = length(times)
    copyto!(R.allW, awsW.ekf_ws.all_params)
    fill!(R.allbar, zero(DW))
    _ctsem_tipred_fill!(awsW.tipreds, so.tipreds, xw)
    _ctsem_tipred_bar_reset!(awsW, so.tipreds, _ctsem_tipred_length(so.tipreds))
    awsW.tape.ngroups = 0
    empty!(R.prog)
    R.ng = 0
    _tan_snap_save!(R.base, mod, R.allW, sp)
    _tan_acc_zero!(acc)
    x, P = sw.x, sw.P
    _tan_copy!(x, mod.T0M)
    _tan_copy!(P, mod.T0V)
    fill!(sw.ll, 0.0)
    hs = sw.hs
    _tan_sweep_begin!(sw, mod)
    snap = 0
    npred = 0
    for t in 1:T
        if t == 1
            # The predict group once at the first row, with no prediction.
            snap = _tan_rows_group!(R, 1, x, mod, sw, sp, so, i, 1, times[1], 0.0, snap, nl)
        else
            dtrow = times[t] - times[t - 1]
            k = _ctsem_substeps(dtrow, rule, t)
            h = dtrow / k
            for s in 1:k
                snap = _tan_rows_group!(R, 1, x, mod, sw, sp, so, i, t,
                    times[t - 1] + s * h, h, snap, nl)
                npred += 1
                npred > length(sw.preds) && push!(sw.preds, TanPredRec(n, sw.nh, W))
                _tan_predict!(x, P, sw.preds[npred], fw, mod, h, hs) || return _TAN_REFUSED
                push!(R.prog, (2, npred, snap))
            end
        end
        dt = t == 1 ? 0.0 : times[t] - times[t - 1]
        snap = _tan_rows_group!(R, 2, x, mod, sw, sp, so, i, t, times[t], dt, snap, nl)
        if mod.ntd > 0
            _tan_td!(x, P, sw, mod, view(TD, :, t), t)
            push!(R.prog, (3, t, snap))
        end
        snap = _tan_rows_group!(R, 3, x, mod, sw, sp, so, i, t, times[t], dt, snap, nl)
        t > length(sw.upds) && push!(sw.upds, TanUpdRec(n, sw.m, sw.nh, W))
        _tan_update!(x, P, sw.upds[t], fw, mod, view(Y, :, t), sw.ll, hs) ||
            return _TAN_INVALID
        push!(R.prog, (4, t, snap))
    end
    isfinite(sw.ll[1]) || return _TAN_INVALID

    # The reverse walk.
    xbar, Pbar = _tan_zero!(sw.xbar), _tan_zero!(sw.Pbar)
    held = -1     # the snapshot the model holds
    factored = -1 # the snapshot the drift's factors are of
    for e in length(R.prog):-1:1
        kind, idx, s = R.prog[e]
        if kind == 1
            # The stretch this group began is over: its cotangents onto the
            # cells, then through the group.
            _tan_rows_flush!(R, acc, sp, R.snaps[s])
            fill!(R.xbarW, zero(DW))
            _reverse_group!(R.allbar, R.xbarW, awsW.tape.groups[idx], sp, awsW)
            _tan_add_lanes!(xbar, R.xbarW, nl)
            continue
        end
        if s != held
            _tan_snap_load!(mod, s == 0 ? R.base : R.snaps[s])
            held = s
        end
        if kind == 4
            _tan_reverse_row_update!(xbar, Pbar, sw.upds[idx], fw, mod, acc, 1.0, hs) ||
                return _TAN_INVALID
        elseif kind == 3
            _tan_reverse_td!(xbar, Pbar, sw, mod, acc, 1.0, view(TD, :, idx), idx)
        else
            if s != factored
                _tan_sweep_begin!(sw, mod)
                factored = s
            end
            _tan_reverse_predict!(xbar, Pbar, sw.preds[idx], rv, fw, mod, acc, 1.0, hs)
            # The exponential's direction belongs to this step's drift.
            _tan_flush_frechet!(acc, fw, mod, sw.sT)
        end
    end
    _tan_lincomb!(acc.T0M, true, acc.T0M, 1.0, xbar)
    _tan_sym!(rv.Ps, Pbar)
    _tan_comb!(acc.T0V, ((1.0, rv.Ps),), true)
    _tan_rows_flush!(R, acc, sp, R.base)
    # The cells back to the raw parameters, and the subject's weight: the walk
    # runs unweighted, since a group's pullback sends the cells' cotangents
    # into the state's, which the kernels carry unweighted.
    vb = _tan_values_reverse(awsW, sp, xw, R.allbar, so.tipreds)
    @inbounds for l in 1:nl, j in eachindex(vb)
        Hb[j, l] += w * ForwardDiff.partials(vb[j], l)
    end
    return _TAN_OK
end
