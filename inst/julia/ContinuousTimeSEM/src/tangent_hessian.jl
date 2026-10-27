# The marginal route's Hessian by the hand-coded second order: each subject's
# filter sweep carries W = 12 directions as Float64 tangents
# (`tangent_sweep.jl`), where `ctsem_hessian`'s other route runs the adjoint
# gradient at a 12-wide dual. The parameter layer -- subject values, the
# regular transforms, the covariance constructions, their pullbacks and the
# prior -- is the engine's own code, run once per column at a width-1 dual:
# forward along the column for the parameter matrices' tangents, and back with
# their cotangents carrying the column's tangent, so the result's partial is
# the column. That layer is shared by every subject where no TI effect reaches
# it, which is the only case taken here (`_tan_refusal`); everything else goes
# to the dual, as does a point where an interval needs the series. The dual
# route is the oracle (`ctsem_set_tangent_hessian!(false)`,
# test_tangent_hessian.jl).

const _CTSEM_TANGENT_HESSIAN = Ref(true)

"""
    ctsem_set_tangent_hessian!(on)

Whether `ctsem_hessian` and `ctsem_hessian_columns` take the hand-coded second
order where it covers the model (`true`, the default) or forward mode over the
adjoint gradient at its dual (`false`): the oracle the hand-coded route is
checked against.
"""
function ctsem_set_tangent_hessian!(on::Bool)
    _CTSEM_TANGENT_HESSIAN[] = on
    return on
end
export ctsem_set_tangent_hessian!

const _TAN_WIDTH = 12

"""
Why the hand-coded Hessian does not serve `objective`, or `nothing`: what its
sweep covers is a linear Gaussian continuous-time model with every state
diffusing, its parameter layer shared by every subject, no stationary or
population block, no integrated deviations, and no forward-mode subjects.
"""
function _tan_refusal(objective::CTSEMObjective)
    subjects = objective.subject_objectives
    isempty(subjects) && return "no subjects"
    sp = objective.params
    sp.continuous_time || return "discrete time"
    any(sp.predict_transforms_indices) && return "per-row predict transforms"
    any(sp.update_transforms_indices) && return "per-row update transforms"
    any(sp.td_transforms_indices) && return "per-row TD transforms"
    any(t -> t != 0, sp.manifesttype) && return "categorical indicators"
    (sp.stationary || _stationary_elements_on(sp.stationary_elements)) &&
        return "stationary moments"
    isempty(sp.population_indices) || return "a population block"
    isempty(sp.ti_parameter_indices) || return "TI predictor effects"
    isempty(sp.ti_tail) || return "a TI predictor tail"
    isempty(objective.ti_missing_parameter) || return "sampled TI predictor values"
    _hybrid_on(sp.hybrid) && return "integrated deviations"
    isfinite(_CTSEM_ADJOINT_GROWTH[]) && return "forward-mode subjects"
    for so in subjects
        so.tipreds isa AbstractVector || return "sampled TI predictor values"
        size(so.tdpreds, 1) == 0 || return "TD predictors"
    end
    return nothing
end

"""The parameter matrices at `xv` (any element type), through the filter's layer."""
function _tan_layer_mats(aws, sp, xv::AbstractVector, tipreds)
    ws = aws.ekf_ws
    _materialize_subject_values!(ws.subject_values, xv, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = CTSEMPars(ws.all_params, sp.parameter_axis)
    n, m = ws.state_dim, ws.manifest_dim
    # As the filter builds them: MANIFESTVAR's and DIFFUSION's covariances as
    # the construction leaves them, T0VAR's with its lower triangle copied up.
    cov(mat, d, buf) = (sdcovsqrt2cov!(buf, mat, sp.covmatcode, d);
        Matrix(view(buf.out, 1:d, 1:d)))
    T0V = cov(pars.T0VAR, n, ws.bufferQ)
    _copy_lower_to_upper!(T0V, n)
    return (JAx=Matrix(view(pars.JAx, 1:n, 1:n)), DRIFT=Matrix(view(pars.DRIFT, 1:n, 1:n)),
        CINT=collect(view(pars.CINT, 1:n)), Qc=cov(pars.DIFFUSION, n, ws.bufferQ),
        LAM=Matrix(view(pars.LAMBDA, 1:m, 1:n)), Jy=Matrix(view(pars.Jy, 1:m, 1:n)),
        MM=collect(view(pars.MANIFESTMEANS, 1:m)), R=cov(pars.MANIFESTVAR, m, ws.bufferΘ),
        T0M=collect(view(pars.T0MEANS, 1:n)), T0V=T0V)
end

"""
The raw gradient of `<C, matrices(xv)>`, the matrices' cotangents `C` pushed
back through the layer: the adjoint's own pullbacks (`_sdcovsqrt2cov_pullback!`,
`_ctsem_regular_pullback!`, `_ctsem_ti_pullback!`), at `xv`'s element type.
"""
function _tan_layer_reverse(aws, sp, xv::AbstractVector{S}, C, tipreds) where {S}
    ws = aws.ekf_ws
    n, m = ws.state_dim, ws.manifest_dim
    _materialize_subject_values!(ws.subject_values, xv, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = CTSEMPars(ws.all_params, sp.parameter_axis)
    tb = zeros(S, length(sp.mutables))
    tc = CTSEMPars(tb, sp.parameter_axis)
    view(tc.JAx, 1:n, 1:n) .+= C.JAx
    view(tc.CINT, 1:n) .+= C.CINT
    view(tc.LAMBDA, 1:m, 1:n) .+= C.LAM
    view(tc.MANIFESTMEANS, 1:m) .+= C.MM
    view(tc.T0MEANS, 1:n) .+= C.T0M
    pullback!(target, mat, cbar, d) = _sdcovsqrt2cov_pullback!(target, mat, cbar, d;
        scratch=aws.covsqrt_scratch, covmatcode=sp.covmatcode)
    pullback!(tc.DIFFUSION, pars.DIFFUSION, C.Qc, n)
    pullback!(tc.MANIFESTVAR, pars.MANIFESTVAR, C.R, m)
    t0 = zeros(S, n, n)
    pullback!(t0, pars.T0VAR, C.T0V, n)
    view(tc.T0VAR, 1:n, 1:n) .+= t0
    svb = zeros(S, length(ws.subject_values))
    _ctsem_regular_pullback!(svb, tb, ws.subject_values, sp, aws.regular_supports,
        aws.regular_dual_scratch)
    vb = zeros(S, length(xv))
    _ctsem_ti_pullback!(vb, svb, sp, tipreds, xv)
    return vb
end

const _TAN_LAYER_NAMES = (:JAx, :CINT, :Qc, :LAM, :MM, :R, :T0M, :T0V)
_tan_tag(::Type{<:ForwardDiff.Dual{T}}) where {T} = T
@inline _tan_layer_field(X, k::Symbol) = getproperty(X, k === :LAM ? :LAMBDA : k)

"""A worker's sweep, kept in its scratch store between calls."""
function _tan_sweep_scratch(n::Int, m::Int, ::Val{W}) where {W}
    store = _ctsem_scratch()
    key = (:tangent_sweep, W, n, m)
    sw = get(store, key, nothing)
    if sw === nothing
        sw = TanSweep(n, m, W)
        store[key] = sw
    end
    return sw::TanSweep{W}
end

"""
The Hessian's columns `cols` into `H` (`length(x)` by `length(cols)`) by the
hand-coded second order. `false` where it does not serve this point -- the
drift is not the Jacobian, or an interval needs the series -- and the caller
takes the dual route; `_tan_refusal(objective)` must be `nothing`.
"""
_tan_hessian_columns!(H::Matrix{Float64}, objective::CTSEMObjective,
    x::Vector{Float64}, cols::Vector{Int}, progress)::Bool =
    _tan_hessian_columns!(H, objective, x, cols, progress, Val(_TAN_WIDTH))

function _tan_hessian_columns!(H::Matrix{Float64}, objective::CTSEMObjective,
    x::Vector{Float64}, cols::Vector{Int}, progress, ::Val{W})::Bool where {W}
    sp = objective.params
    npar = length(x)
    # Width 1 at the Hessian's own tag: the layer's own seeded duals nest inside
    # it as they nest inside the 12-wide one.
    D1 = ForwardDiff.Dual{_tan_tag(_ctsem_hessian_dual(objective)),Float64,1}
    subjects = objective.subject_objectives
    # The layer is every subject's, there being no TI effect: any subject's
    # predictor row gives it.
    tipreds = subjects[1].tipreds
    aws0 = CTSEMAdjointWorkspace(Float64, sp, npar, 0)
    aws1 = CTSEMAdjointWorkspace(D1, sp, npar, 0)
    m0 = _tan_layer_mats(aws0, sp, x, tipreds)
    n, m = size(m0.JAx, 1), size(m0.LAM, 1)
    # Every state diffusing, none static, and the drift its own Jacobian.
    (aws0.diffusion_state_indices == 1:n && aws0.affine_dim == n) || return false
    (m0.JAx == m0.DRIFT && m0.Jy == m0.LAM) || return false
    seed(c) = [D1(x[i], ForwardDiff.Partials((i == c ? 1.0 : 0.0,))) for i in 1:npar]
    mod = TanModel(n, m, W)
    nsub = length(subjects)
    nworkers = _ctsem_nchunks(nsub)
    nblocks = _ctsem_nblocks(nsub, nworkers)
    ranges = _ctsem_chunk_ranges(nsub, nblocks)
    accs = [TanAcc(n, m, W) for _ in 1:nblocks]
    status = zeros(Int, nblocks)
    weights = objective.subject_weights
    for start in 1:W:length(cols)
        chunk = start:min(start + W - 1, length(cols))
        # The matrices' tangents along the chunk's columns.
        for k in _TAN_LAYER_NAMES
            X = _tan_layer_field(mod, k)
            X.v .= getproperty(m0, k)
            fill!(X.d, 0.0)
        end
        for (w, p) in enumerate(chunk)
            mw = _tan_layer_mats(aws1, sp, seed(cols[p]), tipreds)
            (mw.JAx == mw.DRIFT && mw.Jy == mw.LAM) || return false
            for k in _TAN_LAYER_NAMES
                X = _tan_layer_field(mod, k)
                P = ForwardDiff.partials.(getproperty(mw, k), 1)
                X.d isa Matrix ? (X.d[w, :] .= P) : (X.d[w, :, :] .= P)
            end
        end
        foreach(k -> _tan_setnz_from!(_tan_layer_field(mod, k)), _TAN_LAYER_NAMES)
        # The subjects, in the gradient's blocks on its workers; each block's
        # exponential directions flushed once for all its subjects.
        fill!(status, _TAN_OK)
        _ctsem_pull_blocks(nblocks, nworkers) do b, _
            sw = _tan_sweep_scratch(n, m, Val(W))
            acc = _tan_acc_zero!(accs[b])
            _tan_sweep_begin!(sw, mod)
            for i in ranges[b]
                so = subjects[i]
                wt = isempty(weights) ? 1.0 : weights[i]
                wt == 0 && continue
                st = _tan_sweep!(acc, sw, mod, so.timesteps, so.data, so.max_timestep, wt)
                if st != _TAN_OK
                    status[b] = st
                    return nothing
                end
            end
            _tan_flush_frechet!(acc, sw.fw, mod, sw.sT)
            return nothing
        end
        any(==(_TAN_REFUSED), status) && return false
        if any(==(_TAN_INVALID), status)
            # The gradient there is an invalid point's, NaN; so is its Hessian.
            fill!(H, NaN)
            return true
        end
        tot = accs[1]
        for b in 2:nblocks
            _tan_acc_add!(tot, accs[b])
        end
        # Each column: the layer's reverse at a width-1 dual seeded along it,
        # the cotangents carrying their tangent along it, and the prior.
        for (w, p) in enumerate(chunk)
            C = NamedTuple{_TAN_LAYER_NAMES}(map(_TAN_LAYER_NAMES) do k
                X = _tan_layer_field(tot, k)
                dv = X.d isa Matrix ? X.d[w, :] : X.d[w, :, :]
                D1.(X.v, ForwardDiff.Partials.(tuple.(dv)))
            end)
            xs = seed(cols[p])
            g = _tan_layer_reverse(aws1, sp, xs, C, tipreds)
            _ctsem_log_prior_gradient!(g, objective, xs)
            @inbounds for i in 1:npar
                H[i, p] = ForwardDiff.partials(g[i], 1)
            end
        end
        progress === nothing || progress(last(chunk), length(cols), "parameters")
    end
    return true
end
