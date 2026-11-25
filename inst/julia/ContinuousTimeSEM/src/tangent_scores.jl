# The row scores (`ctsem_row_scores`) by the hand-coded forward sweep: each
# data row's log predictive density with its tangents along the scores'
# directions, a chunk of directions a pass as the Hessian's sweeps carry, the
# directions the lanes of the parameter layer at the Hessian's dual. The
# sweep's forward (`_tan_sweep!`), stopped before its reverse, with each row's
# share of the log likelihood kept (`_tan_llrow!`).
# The filter at the Hessian's dual is the tests' oracle (test/dual_oracles.jl).
# A forward sweep, so subjects whose reverse is unstable (`_CTSEM_ADJOINT_GROWTH`)
# take it too.

"""
The row scores into `out` (rows of every subject in order, a column per
direction) by the hand-coded forward sweep, the filter's deviations off; a
subject whose filter fails leaves its rows as they were (`NaN`). An error for
what the sweep does not take (`_tan_check`, `_TAN_NOT_JACOBIAN`).
"""
function _tan_row_scores!(out::Matrix{Float64}, objective::CTSEMObjective,
    theta::Vector{Float64}, directions::AbstractMatrix, rows::Vector{Int})
    _tan_check(objective, _HYBRID_OFF)
    _tan_row_scores!(out, objective, theta, directions, rows, Val(_TAN_WIDTH)) ||
        throw(ArgumentError(_TAN_NOT_JACOBIAN))
    return out
end

function _tan_row_scores!(out::Matrix{Float64}, objective::CTSEMObjective,
    theta::Vector{Float64}, directions::AbstractMatrix, rows::Vector{Int},
    ::Val{W})::Bool where {W}
    sp = objective.params
    npar = length(theta)
    DW = _ctsem_hessian_dual(objective)
    ForwardDiff.npartials(DW) == W || return false
    subjects = objective.subject_objectives
    isempty(subjects) && return true
    rowsmode = _tan_rows_cells(sp)
    shared = !rowsmode && isempty(sp.ti_parameter_indices) && isempty(sp.ti_tail)
    nsub = length(subjects)
    nworkers = _ctsem_nchunks(nsub)
    nblocks = _ctsem_nblocks(nsub, nworkers)
    ranges = _ctsem_chunk_ranges(nsub, nblocks)
    aws0 = CTSEMAdjointWorkspace(Float64, sp, npar, 0)
    m0 = _tan_layer_mats(aws0, sp, theta, subjects[1].tipreds)
    n, m, ntd = size(m0.JAx, 1), size(m0.LAM, 1), size(m0.TDE, 2)
    na, dyn = aws0.affine_dim, collect(aws0.diffusion_state_indices)
    gmask, cat = _tan_cat_masks(sp, m)
    mod = TanModel(n, m, W; ct=sp.continuous_time, rows=rowsmode,
        drows=rowsmode && _tan_drift_by_row(sp), gmask=gmask, na=na, dyn=dyn, ntd=ntd, cat=cat)
    awsW = CTSEMAdjointWorkspace(DW, sp, npar, ntd)
    layers = map(1:nworkers) do _
        awsWk = CTSEMAdjointWorkspace(DW, sp, npar, ntd)
        (awsW=awsWk, mod=shared ? mod : _tan_model_like(mod),
            acc=TanAcc(n, m, ntd, 0, W, cat.nxc), rows=TanRows(awsWk, sp, n, m, ntd, Val(W), DW),
            buf=_tan_layer_buffers(DW, awsWk, sp, npar, ntd))
    end
    offsets = cumsum(vcat(0, rows))
    status = zeros(Int, nblocks)
    K = size(directions, 2)
    for start in 1:W:K
        cols = start:min(start + W - 1, K)
        nl = length(cols)
        xw = [DW(theta[j], ForwardDiff.Partials(ntuple(a -> a <= nl ?
            Float64(directions[j, cols[a]]) : 0.0, Val(W)))) for j in 1:npar]
        shared && (_tan_layer_fill!(mod, awsW, sp, xw, nl, subjects[1].tipreds) ||
            return false)
        fill!(status, _TAN_OK)
        _ctsem_pull_blocks(nblocks, nworkers) do b, wk
            status[b] = _tan_score_block!(out, offsets, rows, cols,
                _tan_sweep_scratch(n, m, ntd, 0, Val(W), na, length(dyn)), layers[wk],
                objective, ranges[b], xw, nl, shared)
            return nothing
        end
        any(==(_TAN_REFUSED), status) && return false
    end
    return true
end

"""One block's subjects' row scores along the chunk's directions `cols`, on the
worker's state `L`: the layer and any stationary moments once for the block
where it is shared (`shared`, `_tan_block_start!`), else each subject's."""
function _tan_score_block!(out::Matrix{Float64}, offsets::Vector{Int}, rows::Vector{Int},
    cols::UnitRange{Int}, sw::TanSweep{W}, L, objective::CTSEMObjective, range,
    xw::Vector, nl::Int, shared::Bool) where {W}
    subjects, sp = objective.subject_objectives, objective.params
    R = L.rows
    noH = zeros(0, 0)
    if shared
        _tan_layer_values!(L.awsW, sp, xw, subjects[first(range)].tipreds)
        copyto!(R.allW, L.awsW.ekf_ws.all_params)
        _tan_sweep_begin!(sw, L.mod)
        # Where the moments do not exist every subject's sweep says so.
        _tan_block_start!(R, sw, L.mod, sp, nl)
    end
    for i in range
        so = subjects[i]
        shared || _tan_layer_fill!(L.mod, L.awsW, sp, xw, nl, so.tipreds, L.buf) ||
            return _TAN_REFUSED
        st = _tan_sweep!(L.acc, noH, sw, L.mod, R, objective, i, 1.0, xw, nl;
            forward=true, shared=shared)
        st == _TAN_REFUSED && return _TAN_REFUSED
        # A subject whose filter fails keeps its rows NaN, as the filter's.
        st == _TAN_OK || continue
        @inbounds for r in 1:rows[i], (a, col) in enumerate(cols)
            out[offsets[i] + r, col] = sw.llrow[(r - 1) * (1 + W) + 1 + a]
        end
    end
    return _TAN_OK
end
