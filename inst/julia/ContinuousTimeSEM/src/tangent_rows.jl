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
    mc::TanVec{W}; momon::Base.RefValue{Bool}
    rawD::Matrix{DW}; rawR::Matrix{DW}
end
TanSnap{W,DW}(n::Int, m::Int, ntd::Int) where {W,DW} = TanSnap{W,DW}(TanMat{W}(n, n),
    TanMat{W}(n, n), TanMat{W}(n, n), TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(m, n),
    TanMat{W}(m, n), TanVec{W}(m), TanMat{W}(m, m), TanMat{W}(n, ntd), TanMat{W}(n, n),
    TanVec{W}(n), Ref(false), zeros(DW, n, n), zeros(DW, m, m))

function _tan_snap_save!(s::TanSnap, mod::TanModel, allW::AbstractVector, sp)
    _tan_copy!(s.JAx, mod.JAx); _tan_copy!(s.DRIFT, mod.DRIFT); _tan_copy!(s.DmJ, mod.DmJ)
    _tan_copy!(s.CINT, mod.CINT); _tan_copy!(s.Qc, mod.Qc)
    _tan_copy!(s.LAMBDA, mod.LAMBDA); _tan_copy!(s.Jy, mod.Jy); _tan_copy!(s.MM, mod.MM)
    _tan_copy!(s.R, mod.R); _tan_copy!(s.TDE, mod.TDE); _tan_copy!(s.Jtd, mod.Jtd)
    _tan_copy!(s.mc, mod.mc); s.momon[] = mod.momon[]
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
    _tan_copy!(mod.mc, s.mc); mod.momon[] = s.momon[]
    return mod
end

"""
One step of an iterated update (`_ekf_iterated_rows!`), with tangents: from
the prior mean `xp` and covariance `P`, the measurement linearised at the
iterate `xl` -- loading `H` (Jy's observed rows), intercept `c`, `R` the
observed block of the covariance the prior evaluated -- and its update's mean.
"""
mutable struct TanRelinRec{W}
    o::Vector{Int}
    xp::TanVec{W}; P::TanMat{W}
    H::TanMat{W}; Lam::TanMat{W}; c::TanVec{W}; xl::TanVec{W}; R::TanMat{W}
    y::TanVec{W}
    first::Bool
    alpha::Float64
end
TanRelinRec(n::Int, mo::Int, ::Val{W}) where {W} = TanRelinRec{W}(Int[], TanVec{W}(n),
    TanMat{W}(n, n), TanMat{W}(mo, n), TanMat{W}(mo, n), TanVec{W}(mo), TanVec{W}(n),
    TanMat{W}(mo, mo), TanVec{W}(mo), true, 1.0)

"""
A worker's state for subjects whose cells change by row: its workspace at the
Hessian's dual (whose tape holds the group records), the subject's cells and
their cotangents, the state and its cotangent as duals, the snapshots, the
layer's own (`base`), and the forward's program: `(kind, index, snapshot)`,
kind 1 a group (its tape slot), 2 a prediction substep, 3 a TD impulse, 4 an
update (the row), 5 a step of an iterated update, each with the snapshot it
read.
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
    # An iterated update's steps (kind 5), the row's measurement covariance as
    # its prior evaluated it, which its steps keep, and the prior mean's
    # cotangent the reverse carries beside the iterates'.
    relins::Vector{TanRelinRec{W}}
    nrl::Int
    Rrow::TanSnap{W,DW}
    xpbar::TanVec{W}
    xn::TanVec{W}
    # A moment-matched step's covariance and cotangents as duals, and its mean
    # correction's cotangent, from its prediction's reverse (`affbar`).
    PW::Matrix{DW}
    PbarW::Matrix{DW}
    cbar::TanVec{W}
end
function TanRows(awsW, sp, n::Int, m::Int, ntd::Int, ::Val{W}, ::Type{DW}) where {W,DW}
    nall = length(sp.mutables)
    return TanRows{W,DW,typeof(awsW)}(awsW, zeros(DW, nall), zeros(DW, nall),
        zeros(DW, n), zeros(DW, n), TanSnap{W,DW}[], TanSnap{W,DW}(n, m, ntd),
        NTuple{3,Int}[], 0, TanRelinRec{W}[], 0, TanSnap{W,DW}(n, m, ntd), TanVec{W}(n),
        TanVec{W}(n), zeros(DW, n, n), zeros(DW, n, n), TanVec{W}(n))
end

"""`out` the duals of `X`'s values and tangents, lane `w` its `w`-th."""
function _tan_to_duals!(out::AbstractArray{DW}, X::Union{TanMat{W},TanVec{W}}) where {W,DW}
    @inbounds for idx in eachindex(X.v)
        o = W * (idx - 1)
        out[idx] = DW(X.v[idx], ForwardDiff.Partials(ntuple(w -> X.d[o + w], Val(W))))
    end
    return out
end

"""`X += D`, the duals' values and their first `nl` lanes as tangents."""
function _tan_add_lanes!(X::Union{TanMat{W},TanVec{W}}, D::AbstractArray, nl::Int) where {W}
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
    snap::Int, nl::Int, keepR::Bool=false) where {W,DW}
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
    mod.momon[] = false
    # An iterated update's steps keep the measurement covariance its prior
    # evaluated (`_ekf_iterated_rows!`).
    keepR && _tan_copy!(mod.R, R.Rrow.R)
    # The drift's factors for the predictions that read it.
    g == 1 && _tan_sweep_begin!(sw, mod)
    R.ng += 1
    R.ng > length(R.snaps) && push!(R.snaps, TanSnap{W,DW}(sw.n, sw.m, mod.ntd))
    _tan_snap_save!(R.snaps[R.ng], mod, R.allW, sp)
    keepR && (R.snaps[R.ng].rawR .= R.Rrow.rawR)
    push!(R.prog, (1, slot, R.ng))
    return R.ng
end

"""
The prediction substep's dynamics moment-matched (`_ekf_moment_cells!`), run as
the engine's own at the Hessian's dual on the subject's cells, after the
predict group at the mean: the cubature's points from `x` and the incoming
covariance `P`, each point's group recorded on the workspace's tape, the
averaged cells written and the mean correction into `mod.mc`; the model
refreshed and kept. Returns the snapshot, `snap` where the covariance does
not factor and the step stays the EKF's.
"""
function _tan_rows_moment!(R::TanRows{W,DW}, x::TanVec{W}, P::TanMat{W}, mod::TanModel{W},
    sw::TanSweep{W}, sp, so, subject::Int, t::Int, time::Float64, dt::Float64,
    snap::Int, nl::Int) where {W,DW}
    awsW = R.awsW
    ws = awsW.ekf_ws
    _tan_to_duals!(ws.state, x)
    _tan_to_duals!(R.PW, P)
    pars = CTSEMPars(R.allW, sp.parameter_axis)
    ctx = CTSEMRowContext(ws.state, pars, view(so.tdpreds, :, t), so.tipreds, time, dt,
        subject, t)
    _ekf_moment_cells!(ws, pars, sp, R.allW, ctx, awsW.tape, R.PW) || return snap
    _tan_rows_refresh!(mod, R, sp, nl)
    _tan_lanes!(mod.mc, view(ws.moment_c, 1:sw.n), nl)
    mod.momon[] = true
    _tan_sweep_begin!(sw, mod)
    R.ng += 1
    R.ng > length(R.snaps) && push!(R.snaps, TanSnap{W,DW}(sw.n, sw.m, mod.ntd))
    _tan_snap_save!(R.snaps[R.ng], mod, R.allW, sp)
    push!(R.prog, (6, awsW.tape.nmoments, R.ng))
    return R.ng
end

"""
Row `t`'s categorical observations (`_ekf_binary_rows!`), the engine's own at
the Hessian's dual on the worker's workspace: the state and covariance in as
duals and back, the row's log likelihood's share into `sw.ll`, the record on
the workspace's tape (program kind 7). Returns the snapshot, or -1 where an
observation has no support.
"""
function _tan_rows_categorical!(x::TanVec{W}, P::TanMat{W}, R::TanRows{W,DW},
    mod::TanModel{W}, sw::TanSweep{W}, sp, Y::AbstractMatrix, t::Int, snap::Int,
    nl::Int) where {W,DW}
    any(k -> isfinite(Y[k, t]) && !mod.gmask[k], axes(Y, 1)) || return snap
    awsW = R.awsW
    ws = awsW.ekf_ws
    n = sw.n
    observed = Int[k for k in axes(Y, 1) if isfinite(Y[k, t])]
    _tan_to_duals!(ws.state, x)
    _tan_to_duals!(ws.P_predict.data, P)
    pars = CTSEMPars(R.allW, sp.parameter_axis)
    record = _record_binary!(awsW.tape, ws, pars, Y, t, _ekf_binary_subset(ws, observed), n)
    ll = _ekf_binary_rows!(ws, pars, Y, t, observed, nothing, record)
    ll === nothing && return -1
    _tan_from_lanes!(x, view(ws.state, 1:n), nl)
    _tan_from_lanes!(P, ws.P_predict.data, nl)
    sw.ll[1] += ForwardDiff.value(ll)
    p = ForwardDiff.partials(ll)
    @inbounds for w in 1:nl
        sw.ll[1 + w] += p[w]
    end
    push!(R.prog, (7, awsW.tape.nbinaries, snap))
    return snap
end

"""`X`'s values and tangents set from the duals `D`, lanes past `nl` zero."""
function _tan_from_lanes!(X::Union{TanMat{W},TanVec{W}}, D::AbstractArray, nl::Int) where {W}
    return _tan_lanes!(X, D, nl)
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
    awsW.tape.nmoments = 0
    awsW.tape.nbinaries = 0
    mod.momon[] = false
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
    R.nrl = 0
    relin = awsW.ekf_ws.relinearise && _CTSEM_RELINEARISE[]
    moment = awsW.ekf_ws.moment && _CTSEM_MOMENT_PREDICT[]
    cat = !isempty(mod.gmask)
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
                moment && (snap = _tan_rows_moment!(R, x, P, mod, sw, sp, so, i, t,
                    times[t - 1] + s * h, h, snap, nl))
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
        upd = sw.upds[t]
        upd.relin = 0
        # A row's categorical observations first, each an exact scalar
        # conditioning, then its Gaussian ones from there.
        if cat && _tan_rows_categorical!(x, P, R, mod, sw, sp, Y, t, snap, nl) < 0
            return _TAN_INVALID
        end
        if relin && any(k -> isfinite(Y[k, t]) && (isempty(mod.gmask) || mod.gmask[k]),
            axes(Y, 1))
            snap = _tan_relin_row!(x, P, upd, R, mod, sw, sp, so, i, t, times[t], dt, snap,
                nl, view(Y, :, t))
            snap < 0 && return _TAN_INVALID
        end
        _tan_update!(x, P, upd, fw, mod, view(Y, :, t), sw.ll, hs) || return _TAN_INVALID
        push!(R.prog, (4, t, snap))
    end
    isfinite(sw.ll[1]) || return _TAN_INVALID

    # The reverse walk.
    xbar, Pbar = _tan_zero!(sw.xbar), _tan_zero!(sw.Pbar)
    _tan_zero!(R.xpbar)
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
        if kind == 7
            # A row's categorical observations, by the engine's own reverse at
            # the Hessian's dual (`_reverse_binary!`).
            _tan_to_duals!(R.xbarW, xbar)
            _tan_to_duals!(R.PbarW, Pbar)
            _reverse_binary!(R.xbarW, R.PbarW, CTSEMPars(R.allbar, sp.parameter_axis),
                awsW.tape.binaries[idx], n)
            _tan_from_lanes!(xbar, R.xbarW, nl)
            _tan_from_lanes!(Pbar, R.PbarW, nl)
            continue
        elseif kind == 6
            # A moment-matched step's cubature: the averaged cells' cotangents
            # onto the cells, its correction's from the prediction's reverse,
            # then through the points to the mean and the covariance.
            _tan_rows_flush!(R, acc, sp, R.snaps[s])
            ms = awsW.moment
            _tan_to_duals!(view(ms.cbar, 1:n), R.cbar)
            fill!(R.xbarW, zero(DW))
            fill!(R.PbarW, zero(DW))
            _reverse_moment!(R.xbarW, R.PbarW, R.allbar,
                CTSEMPars(R.allbar, sp.parameter_axis), awsW.tape.moments[idx], awsW.tape,
                sp, awsW, n, sw.m)
            _tan_add_lanes!(xbar, R.xbarW, nl)
            _tan_add_lanes!(Pbar, R.PbarW, nl)
            continue
        elseif kind == 5
            _tan_reverse_relin!(xbar, R.xpbar, Pbar, R.relins[idx], acc, fw) ||
                return _TAN_INVALID
        elseif kind == 4
            upd = sw.upds[idx]
            _tan_reverse_row_update!(xbar, Pbar, upd, fw, mod, acc, 1.0, hs) ||
                return _TAN_INVALID
            # An iterated update's prior mean and linearisation point.
            if upd.relin > 0
                u = _tan_row(fw, length(upd.o)).u
                _tan_relin_route!(xbar, R.xpbar, upd.Lam, u.H, u.ybarbar, upd.relin == 1,
                    _tan_row(fw, length(upd.o)).H)
            end
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
            # A moment-matched step's correction's cotangent, its offset's.
            mod.momon[] && _tan_copy!(R.cbar, rv.affbar)
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

"""
The measurement linearised at the state `x` on the observed rows `o`, into
`rr`: Jy's and LAMBDA's rows, the intercept `MANIFESTMEANS + (LAMBDA - Jy) x`
(`_ekf_relin_offset!`), the linearisation point, and `R`'s block.
"""
function _tan_relin_linearise!(rr::TanRelinRec{W}, mod::TanModel{W}, x::TanVec{W},
    o::Vector{Int}) where {W}
    n, mo = length(x.v), length(o)
    if size(rr.H.v, 1) != mo
        rr.H, rr.Lam, rr.c, rr.R, rr.y = TanMat{W}(mo, n), TanMat{W}(mo, n), TanVec{W}(mo),
            TanMat{W}(mo, mo), TanVec{W}(mo)
    end
    resize!(rr.o, mo)
    copyto!(rr.o, o)
    _tan_gather!(rr.H, mod.Jy, o, 1:n)
    _tan_gather!(rr.Lam, mod.LAMBDA, o, 1:n)
    _tan_gather!(rr.c, mod.MM, o)
    _tan_mulvec!(rr.c, rr.Lam, x, true, true)
    _tan_mulvec!(rr.c, rr.H, x, -1.0, true)
    _tan_copy!(rr.xl, x)
    _tan_gather!(rr.R, mod.R, o, o)
    return rr
end

"""
The mean of the update linearised as `rr` says, from its prior
(`_ekf_relin_mean!`): `xp + Pr H' S^-1 (y - c - H xp)`, `Pr = P + eps I`,
`S = sym(H Pr H' + R) + eps I`, into `xn`. `false` where S does not factor or
the mean is not finite.
"""
function _tan_relin_mean!(xn::TanVec{W}, rr::TanRelinRec{W}, rs::TanRowScratch{W}) where {W}
    yt = rs.yt
    _tan_copy!(yt, rr.c)
    _tan_mulvec!(yt, rr.H, rr.xp, true, true)
    _tan_lincomb!(yt, true, rr.y, -1.0, yt)
    _tan_copy!(rs.Pr, rr.P)
    _tan_diag!(rs.Pr, 1e-10)
    _tan_mulNT!(rs.PHt, rs.Pr, rr.H)
    _tan_mul!(rs.S, rr.H, rs.PHt)
    _tan_comb!(rs.S, ((1.0, rr.R),), true)
    _tan_sym!(rs.Xc, rs.S)
    _tan_copy!(rs.S, rs.Xc)
    _tan_diag!(rs.S, 1e-10)
    _tan_cholinv!(rs.Ui, rs.S, rs.Xc, rs.Tc, rs.F) || return false
    _tan_mulTvec!(rs.e, rs.Ui, yt, true, false, Val(:AupperT))
    _tan_mulvec!(rs.a, rs.Ui, rs.e, true, false, Val(:Aupper))
    _tan_copy!(xn, rr.xp)
    _tan_mulvec!(xn, rs.PHt, rs.a, true, true)
    return all(isfinite, xn.v)
end

"""
Row `t`'s iterated update before its final update (`_ekf_iterated_rows!`):
the frozen count of steps where the subject's plan gives one, else to the
mode, halving the steps where they stop shrinking; each step recorded
(program kind 5) and the update group evaluated at the new iterate, the
measurement covariance kept as the prior evaluated it. Leaves `x` at the prior
mean and `upd` set for the final update, linearised at the last iterate.
Returns the snapshot the model holds, or -1 where a step fails.
"""
function _tan_relin_row!(x::TanVec{W}, P::TanMat{W}, upd::TanUpdRec{W},
    R::TanRows{W,DW}, mod::TanModel{W}, sw::TanSweep{W}, sp, so, i::Int, t::Int,
    time::Float64, dt::Float64, snap::Int, nl::Int,
    y::AbstractVector{Float64}) where {W,DW}
    n = sw.n
    o = Int[k for k in eachindex(y) if isfinite(y[k]) &&
        (isempty(mod.gmask) || mod.gmask[k])]
    mo = length(o)
    rs = _tan_row(sw.fw, mo)
    # The covariance the prior evaluated, and the prior mean.
    _tan_snap_save!(R.Rrow, mod, R.allW, sp)
    xp = R.xn
    _tan_copy!(xp, x)
    first = true
    steps = _ctsem_relin_steps(so.relin, t)
    alpha = 1.0
    last = Inf
    for it in 1:(steps >= 0 ? steps : _CTSEM_RELIN_MAXIT)
        R.nrl += 1
        R.nrl > length(R.relins) && push!(R.relins, TanRelinRec(n, mo, Val(W)))
        rr = R.relins[R.nrl]
        _tan_relin_linearise!(rr, mod, x, o)
        _tan_copy!(rr.xp, xp)
        _tan_copy!(rr.P, P)
        @inbounds for k in 1:mo
            rr.y.v[k] = y[o[k]]
        end
        _tan_setnz!(rr.y, false)
        xn = rs.u.xbar   # scratch: the reverse's, unused in the forward
        _tan_relin_mean!(xn, rr, rs) || return -1
        if steps < 0
            step = 0.0
            @inbounds for j in 1:n
                step = max(step, abs(xn.v[j] - x.v[j]) / sqrt(abs(P.v[j, j]) + 1e-10))
            end
            if step <= _CTSEM_RELIN_TOL
                R.nrl -= 1
                break
            end
            step >= last && (alpha = max(alpha / 2, 1 / 1024))
            last = step
        end
        rr.first = first
        rr.alpha = alpha
        push!(R.prog, (5, R.nrl, snap))
        # x <- x + alpha (xn - x).
        _tan_lincomb!(x, 1.0 - alpha, x, alpha, xn)
        first = false
        snap = _tan_rows_group!(R, 3, x, mod, sw, sp, so, i, t, time, dt, snap, nl, true)
    end
    # The final update, linearised at the last iterate, from the prior.
    upd.relin = first ? 1 : 2
    mo == size(upd.Lam.v, 1) || (upd.c = TanVec{W}(mo); upd.Lam = TanMat{W}(mo, n))
    _tan_gather!(upd.Lam, mod.LAMBDA, o, 1:n)
    _tan_gather!(upd.c, mod.MM, o)
    _tan_mulvec!(upd.c, upd.Lam, x, true, true)
    _tan_gather!(rs.H, mod.Jy, o, 1:n)
    _tan_mulvec!(upd.c, rs.H, x, -1.0, true)
    _tan_copy!(upd.xl, x)
    _tan_copy!(x, xp)
    return snap
end

"""
The end of an iterated update's reverse (`_relin_route!`): `xbar`, the prior
mean's cotangent, joins `xpbar` and leaves as the linearisation point's,
through the intercept `c = MANIFESTMEANS + (Lam - H) xl`, `-(Lam - H)' ybar`;
at the first step the point is the prior mean, and `xpbar` is spent. `D` is
scratch the shape of `H`.
"""
function _tan_relin_route!(xbar::TanVec{W}, xpbar::TanVec{W}, Lam::TanMat{W}, H::TanMat{W},
    ybar::TanVec{W}, first::Bool, D::TanMat{W}) where {W}
    _tan_lincomb!(xpbar, true, xpbar, true, xbar)
    _tan_comb!(D, ((1.0, Lam), (-1.0, H)))
    _tan_mulTvec!(xbar, D, ybar, -1.0, false)
    if first
        _tan_lincomb!(xbar, true, xbar, true, xpbar)
        _tan_zero!(xpbar)
    end
    return xbar
end

"""
Undo one step of an iterated update (`_reverse_relin!`): `xbar` the step's
cotangent, `alpha` of it the mean's; the prior mean's share into `xpbar`, the
prior covariance's into `Pbar`, the measurement's and its covariance's into
`acc`, and `xbar` left as the linearisation point's. `false` where S does not
factor.
"""
function _tan_reverse_relin!(xbar::TanVec{W}, xpbar::TanVec{W}, Pbar::TanMat{W},
    rr::TanRelinRec{W}, acc::TanAcc{W}, fw::TanForward{W}) where {W}
    o = rr.o
    mo = length(o)
    n = fw.n
    rs = _tan_row(fw, mo)
    u = rs.u
    H, xp, xl = rr.H, rr.xp, rr.xl
    # The forward's intermediates, with its ridges and its factorization.
    Pr, PHt, S, Ui = rs.Pr, rs.PHt, rs.S, rs.Ui
    _tan_copy!(Pr, rr.P)
    _tan_diag!(Pr, 1e-10)
    _tan_mulNT!(PHt, Pr, H)
    _tan_mul!(S, H, PHt)
    _tan_comb!(S, ((1.0, rr.R),), true)
    _tan_sym!(rs.Xc, S)
    _tan_copy!(S, rs.Xc)
    _tan_diag!(S, 1e-10)
    _tan_cholinv!(Ui, S, rs.Xc, rs.Tc, rs.F) || return false
    yt, e, a = rs.yt, rs.e, rs.a
    _tan_copy!(yt, rr.c)
    _tan_mulvec!(yt, H, xp, true, true)
    _tan_lincomb!(yt, true, rr.y, -1.0, yt)
    _tan_mulTvec!(e, Ui, yt, true, false, Val(:AupperT))
    _tan_mulvec!(a, Ui, e, true, false, Val(:Aupper))
    # The step: alpha xbar reaches the mean, the rest the point it left.
    keep = u.We
    _tan_lincomb!(keep, 1.0 - rr.alpha, xbar, false, xbar)
    _tan_scale!(xbar, rr.alpha)
    v = u.HKx
    _tan_mulvec!(v, Pr, xbar)                              # Pr xbar
    wv, ybar = u.Kx, u.yt
    _tan_mulvec!(wv, H, v)                                 # abar = H Pr xbar
    _tan_mulTvec!(u.a, Ui, wv, true, false, Val(:AupperT))
    _tan_mulvec!(ybar, Ui, u.a, true, false, Val(:Aupper))   # ybar = S^-1 abar
    Ha = u.Ky
    _tan_mulTvec!(Ha, H, a)                                # H' a
    Sbar = u.Y
    _tan_ger!(Sbar, -0.5, ybar, a, false)                  # Sbar = -sym(ybar a')
    _tan_ger!(Sbar, -0.5, a, ybar, true)
    _tan_ger!(Pbar, 0.5, xbar, Ha, true)                   # Pbar += sym(xbar (H'a)')
    _tan_ger!(Pbar, 0.5, Ha, xbar, true)
    SH, Hbar = rs.L, u.Hbar
    _tan_mul!(SH, Sbar, H)
    _tan_mulTN!(Pbar, H, SH, true, true)                   # Pbar += H' Sbar H
    _tan_mul!(Hbar, SH, Pr, 2.0, false)                    # Hbar = 2 Sbar H Pr
    _tan_ger!(Hbar, true, a, v, true)                      # + a v'
    dxv = u.Ppx
    _tan_lincomb!(dxv, true, xl, -1.0, xp)
    _tan_ger!(Hbar, true, ybar, dxv, true)                 # + ybar (xl - xp)'
    allstates = 1:n
    _tan_scatter!(acc.MM, -1.0, ybar, o)
    G = rs.H
    _tan_ger!(G, true, ybar, xl, false)
    _tan_scatter!(acc.LAMBDA, -1.0, G, o, allstates)
    _tan_scatter!(acc.Jy, 1.0, Hbar, o, allstates)
    _tan_scatter!(acc.R, 1.0, Sbar, o, o)
    _tan_mulTvec!(xbar, H, ybar, -1.0, true)               # xp's: xbar - H' ybar
    _tan_relin_route!(xbar, xpbar, rr.Lam, H, ybar, rr.first, rs.H)
    _tan_lincomb!(xbar, true, xbar, true, keep)
    return true
end
