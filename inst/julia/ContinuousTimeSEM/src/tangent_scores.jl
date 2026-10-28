# The row scores (`ctsem_row_scores`) by the hand-coded forward sweep: each
# data row's log predictive density with its tangents along the scores'
# directions, a chunk of directions a pass as the Hessian's sweeps carry, the
# directions the lanes of the parameter layer at the Hessian's dual. The
# forward of `tangent_sweep.jl` and `tangent_rows.jl`, stopped before their
# reverse, with each row's share of the log likelihood kept (`_tan_llrow!`).
# The filter at the Hessian's dual (`_row_scores!`) is the oracle.

"""
The row scores into `out` (rows of every subject in order, a column per
direction) by the hand-coded forward sweep, the filter's deviations off as
`_row_scores!` has them. `false` where the sweep does not cover the model or
this point, and the caller takes the dual route; a subject whose filter fails
leaves its rows as they were (`NaN`).
"""
function _tan_row_scores!(out::Matrix{Float64}, objective::CTSEMObjective,
    theta::Vector{Float64}, directions::AbstractMatrix, rows::Vector{Int})::Bool
    _CTSEM_TANGENT_HESSIAN[] || return false
    _tan_refusal(objective, _HYBRID_OFF) === nothing || return false
    return _tan_row_scores!(out, objective, theta, directions, rows, Val(_TAN_WIDTH))
end

function _tan_row_scores!(out::Matrix{Float64}, objective::CTSEMObjective,
    theta::Vector{Float64}, directions::AbstractMatrix, rows::Vector{Int},
    ::Val{W})::Bool where {W}
    sp = objective.params
    npar = length(theta)
    DW = _ctsem_hessian_dual(objective)
    ForwardDiff.npartials(DW) == W || return false
    subjects = objective.subject_objectives
    rowsmode = _tan_rows_model(sp)
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
    mod = TanModel(n, m, W; ct=sp.continuous_time, stat=sp.stationary, rows=rowsmode,
        gmask=gmask, na=na, dyn=dyn, ntd=ntd, cat=cat)
    awsW = CTSEMAdjointWorkspace(DW, sp, npar, ntd)
    layers = shared ? nothing : map(1:nworkers) do _
        awsWk = CTSEMAdjointWorkspace(DW, sp, npar, ntd)
        (awsW=awsWk, mod=_tan_model_like(mod), acc=TanAcc(n, m, ntd, 0, W, cat.nxc),
            rows=rowsmode ? TanRows(awsWk, sp, n, m, ntd, Val(W), DW) : nothing,
            buf=_tan_layer_buffers(DW, awsWk, sp, npar, ntd))
    end
    accs = [TanAcc(n, m, ntd, 0, W, cat.nxc) for _ in 1:nworkers]
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
                _tan_sweep_scratch(n, m, ntd, 0, Val(W), na, length(dyn)), mod, accs[wk],
                layers === nothing ? nothing : layers[wk], objective, ranges[b], xw, nl)
            return nothing
        end
        any(==(_TAN_REFUSED), status) && return false
    end
    return true
end

"""One block's subjects' row scores along the chunk's directions `cols`."""
function _tan_score_block!(out::Matrix{Float64}, offsets::Vector{Int}, rows::Vector{Int},
    cols::UnitRange{Int}, sw::TanSweep{W}, mod::TanModel{W}, acc::TanAcc{W}, L,
    objective::CTSEMObjective, range, xw::Vector, nl::Int) where {W}
    subjects, sp = objective.subject_objectives, objective.params
    L === nothing && _tan_sweep_begin!(sw, mod)
    for i in range
        so = subjects[i]
        if L === nothing
            st = _tan_sweep!(acc, sw, mod, so.timesteps, so.data, so.tdpreds,
                so.max_timestep, 1.0; forward=true)
        else
            _tan_layer_fill!(L.mod, L.awsW, sp, xw, nl, so.tipreds, L.buf) ||
                return _TAN_REFUSED
            st = if L.rows !== nothing
                _tan_sweep_rows!(L.acc, zeros(0, 0), sw, L.mod, L.rows, objective, i, 1.0,
                    xw, nl; forward=true)
            else
                _tan_sweep_begin!(sw, L.mod)
                _tan_sweep!(L.acc, sw, L.mod, so.timesteps, so.data, so.tdpreds,
                    so.max_timestep, 1.0; forward=true)
            end
        end
        st == _TAN_REFUSED && return _TAN_REFUSED
        # A subject whose filter fails keeps its rows NaN, as the filter's.
        st == _TAN_OK || continue
        @inbounds for r in 1:rows[i], (a, col) in enumerate(cols)
            out[offsets[i] + r, col] = sw.llrow[(r - 1) * (1 + W) + 1 + a]
        end
    end
    return _TAN_OK
end
