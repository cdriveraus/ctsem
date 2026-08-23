using ForwardDiff
using ChainRulesCore
using DiffResults
using LinearAlgebra

export extended_kalman_log_likelihood_continuous, grad_log_likelihood_ekf_continuous
export res_and_grad_likelihood_ekf_continuous
export ContinuousEKFGradientWorkspace, grad_log_likelihood_ekf_continuous!
export res_and_grad_likelihood_ekf_continuous!

################################################################################
# Continuous-time Extended Kalman filter
################################################################################

mutable struct ContinuousEKFObjective{P,D,TS,TD,TI,MT}
    value_ws::Any
    dual_ws::Any
    params::P
    data::D
    timesteps::TS
    tdpreds::TD
    tipreds::TI
    subject::Int
    max_timestep::MT
end

function _validate_continuous_ekf_inputs(timesteps::Vector, data::Matrix)
    isempty(timesteps) && throw(ArgumentError("timesteps must contain at least one observation"))
    length(timesteps) == size(data, 2) || throw(DimensionMismatch(
        "timesteps length must equal the number of data columns",
    ))
    @inbounds for idx in 2:length(timesteps)
        timesteps[idx] >= timesteps[idx - 1] || throw(ArgumentError(
            "timesteps must be nondecreasing within one subject; evaluate separate subjects in separate EKF workspaces",
        ))
    end
    return nothing
end

function ContinuousEKFObjective(params::EKFParameters, data::Matrix, timesteps::Vector;
    tdpreds::AbstractMatrix=zeros(eltype(data), 0, size(data, 2)),
    tipreds::AbstractVector=eltype(data)[], subject::Integer=1,
    max_timestep::Real=Inf)
    _validate_continuous_ekf_inputs(timesteps, data)
    size(tdpreds, 2) == size(data, 2) || throw(DimensionMismatch("TD predictor columns must match observations"))
    max_timestep > 0 || throw(ArgumentError("max_timestep must be positive"))
    return ContinuousEKFObjective{typeof(params),typeof(data),typeof(timesteps),typeof(tdpreds),typeof(tipreds),typeof(max_timestep)}(
        nothing, nothing, params, data, timesteps, tdpreds, tipreds, Int(subject), max_timestep)
end

function _get_or_init_objective_workspace!(objective::ContinuousEKFObjective, ::Type{T}) where {T}
    if T <: ForwardDiff.Dual
        ws = objective.dual_ws
        if ws === nothing || !(ws isa ContinuousEKFWorkspace) || eltype(ws.all_params) != T
            ws = _init_continuous_ekf_workspace(T, objective.params)
            objective.dual_ws = ws
        end
        return ws::ContinuousEKFWorkspace{T}
    end

    ws = objective.value_ws
    if ws === nothing || !(ws isa ContinuousEKFWorkspace) || eltype(ws.all_params) != T
        ws = _init_continuous_ekf_workspace(T, objective.params)
        objective.value_ws = ws
    end
    return ws::ContinuousEKFWorkspace{T}
end

function (objective::ContinuousEKFObjective)(p::AbstractVector{T}) where {T}
    ws = _get_or_init_objective_workspace!(objective, eltype(p))
    return _extended_kalman_filter_continuous!(ws, p, objective.data, objective.timesteps,
        objective.params, objective.tdpreds, objective.tipreds, objective.subject,
        objective.max_timestep)::T
end

"""
    ContinuousEKFGradientWorkspace(values, params, timesteps, data)

Create reusable ForwardDiff and EKF workspaces for repeated continuous-time EKF
gradient evaluations.

The workspace is valid for parameter vectors with the same length and scalar
type as `values`, and for the same `params`, `timesteps`, and `data`.
"""
struct ContinuousEKFGradientWorkspace{O,C,R,DR}
    objective::O
    cfg::C
    result::R
    diff_result::DR
end

function ContinuousEKFGradientWorkspace(values::AbstractVector, params::EKFParameters, timesteps::Vector, data::Matrix)
    objective = ContinuousEKFObjective(params, data, timesteps)
    cfg = ForwardDiff.GradientConfig(objective, values)
    result = similar(values)
    diff_result = DiffResults.GradientResult(values)
    return ContinuousEKFGradientWorkspace(objective, cfg, result, diff_result)
end

function _check_gradient_workspace_args(result::AbstractVector, values::AbstractVector)
    @boundscheck length(result) == length(values) || throw(DimensionMismatch("result length must match values length"))
    return nothing
end

"""
    _ekf_predict_step!(ws, pars, Δt)

Advance the EKF state and covariance by one continuous-time prediction step.

The continuous dynamics in `pars` are discretized for the interval `Δt`, and the
prediction is written back into `ws.state` and `ws.P_predict`.
"""
@inline function _ekf_predict_step!(ws::ContinuousEKFWorkspace, pars, Δt)
    # Refresh process noise and map continuous-time dynamics to discrete time.
    # Discretization provides:
    #   A_d = dDRIFT, b_d = dINT, Q_d = dDIFFUSION
    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.DIFFUSION, 0, ws.state_dim)
    _compute_discrete_time_form!(ws.discrete_ca, ws.bufferQ, ws.bufferQ.out, pars, Δt,
        ws.exp_buffer, ws.lyap_buffer, ws.state, ws.diffusion_state_indices,
        ws.diffusion_buffer, ws.discretization_buffer, ws.state_dim,
        ws.discretization_cache)

    # Predict mean and covariance.
    #   x_{t|t-1} = A_d * x_{t-1|t-1} + b_d
    #   P_{t|t-1} = A_d * P_{t-1|t-1} * A_d' + Q_d
    _matvec_mul!(ws.bufferQ.r, ws.discrete_ca.dDRIFT, ws.state, ws.state_dim, ws.state_dim)
    @. ws.state = ws.bufferQ.r + ws.discrete_ca.dINT

    # Stan ridges the incoming posterior covariance before every propagation
    # (`etacov = quad_form_sym(makesym(etacov,verbose,1), eJAx')` in
    # ctModelWriter.R): the +1e-10 diagonal ridge is baked into the transition
    # itself, not just applied defensively before a Cholesky call. ws.P_update
    # is fully overwritten again before this predict step is next entered
    # (either by the measurement update or by a copy from ws.P_predict for the
    # next substep), so ridging it in place here is safe and matches Stan row
    # for row rather than only at the final innovation-covariance Cholesky.
    n = _val(ws.state_dim)
    @inbounds for i in 1:n
        ws.P_update.data[i, i] += 1e-10
    end

    mul!(ws.bufferQ.intermediate, ws.discrete_ca.dDRIFT, ws.P_update)
    _mul_right_transpose!(ws.P_predict.data, ws.bufferQ.intermediate, ws.discrete_ca.dDRIFT, ws.state_dim, ws.state_dim, ws.state_dim)
    ws.P_predict.data .+= ws.discrete_ca.dDIFFUSION
    _copy_lower_to_upper!(ws.P_predict.data, ws.state_dim)
    return nothing
end

@inline function _copy_lower_to_upper!(matrix, state_dim::Val{n}) where {n}
    @inbounds for col in 1:n
        for row in (col + 1):n
            matrix[col, row] = matrix[row, col]
        end
    end
    return matrix
end

"""
    _ridge_diagonal!(mat, d, amount)

Add `amount` to each of the first `d` diagonal entries of `mat` in place.

Shared by the predict step and `_ekf_masked_update_step!` so the Stan-matching
`makesym()` ridge is applied identically everywhere instead of being
duplicated as separate loops.
"""
@inline function _ridge_diagonal!(mat, d::Int, amount)
    @inbounds for i in 1:d
        mat[i, i] += amount
    end
    return mat
end

"""
    _symmetrize_and_ridge!(mat, d)

Average `mat` with its transpose over its first `d`×`d` block (guarding
against asymmetry from floating-point roundoff), then apply the same
`+1e-10` ridge as `_ridge_diagonal!`.
"""
@inline function _symmetrize_and_ridge!(mat, d::Int)
    @inbounds for j in 1:d
        for i in (j + 1):d
            symmetric_value = (mat[i, j] + mat[j, i]) / 2
            mat[i, j] = symmetric_value
            mat[j, i] = symmetric_value
        end
    end
    return _ridge_diagonal!(mat, d, 1e-10)
end

@inline _ctsem_observed(x) = !ismissing(x) && isfinite(x)

"""
    _ekf_update_observed!(ws, pars, data, obs_col, log2π_const)

Update the filter using only observed manifest entries in one row.

A fully missing row leaves the predicted state/covariance in place and
contributes zero to the likelihood -- matching Stan's own
`if(si==0 || nobs_y[rowi] > 0 || dosmoother)` gate around the whole
measurement block. Otherwise this calls `_ekf_masked_update_step!` with the
observed row indices (`1:manifest_dim` when every row is observed, a
`Vector{Int}` otherwise): there is one Kalman-update implementation, not a
separate fast path for the fully-observed case, since Julia compiles a
distinct, `UnitRange`-specialised method for the former that is exactly as
allocation-free and BLAS-eligible as a hand-written fast path would be.
"""
function _ekf_update_observed!(ws::ContinuousEKFWorkspace, pars,
    data::AbstractMatrix, obs_col::Int, log2π_const, trace=nothing, generate=nothing)
    m_full = _val(ws.manifest_dim)
    n_observed = 0
    @inbounds for i in 1:m_full
        _ctsem_observed(data[i, obs_col]) && (n_observed += 1)
    end
    if n_observed == 0
        copyto!(ws.P_update.data, ws.P_predict.data)
        return zero(eltype(ws.state))
    end

    if n_observed == m_full
        observed = 1:m_full
    else
        # `observed` (which manifest rows are present) necessarily allocates
        # here -- its length varies row to row, so it cannot be a workspace
        # field sized once at construction -- but everything downstream of it
        # uses views into existing workspace buffers instead of allocating
        # fresh matrices (see `_ekf_masked_update_step!`).
        observed = Vector{Int}(undef, n_observed)
        idx = 0
        @inbounds for i in 1:m_full
            if _ctsem_observed(data[i, obs_col])
                idx += 1
                observed[idx] = i
            end
        end
    end

    # Recorded before the update runs: `_ekf_masked_update_step!` overwrites
    # ws.state in place, and the reverse pass needs the prior mean.
    _record_update!(trace, ws, pars, data, obs_col, observed,
        ws.state, ws.P_predict.data, _val(ws.state_dim))

    factor = _ekf_masked_update_step!(ws, pars, data, obs_col, observed, generate)
    factor === nothing && return nothing
    return _kalman_loglikelihood_cholesky!(view(ws.ll_buffer, 1:n_observed),
        factor, view(ws.ỹ, 1:n_observed), log2π_const)
end

"""
    _ekf_masked_update_step!(ws, pars, data, obs_col, observed)

Apply one EKF measurement update using only the manifest rows listed in
`observed` (an arbitrary subset, or `1:manifest_dim` for full observation).

Every quantity that is conceptually manifest_dim-by-manifest_dim (or
manifest_dim-long) is a `view` into a `ContinuousEKFWorkspace` buffer here,
restricted to `observed`, so no new matrices are allocated regardless of how
many rows call this or how many manifest variables are observed in them.

Matches Stan's own distinction between the manifest function and its
Jacobian: `pars.LAMBDA` (evaluated at the current state) predicts the
manifest mean for the innovation, while `pars.Jy` (its Jacobian) propagates
covariance -- these coincide for models with a fixed/linear LAMBDA but can
differ for state-dependent measurement models.
"""
function _ekf_masked_update_step!(ws::ContinuousEKFWorkspace, pars,
    data::AbstractMatrix, obs_col::Int, observed::AbstractVector{Int}, generate=nothing)
    m = length(observed)
    n = _val(ws.state_dim)

    Lv = view(pars.LAMBDA, observed, :)
    μv = view(pars.MANIFESTMEANS, observed)
    predview = view(ws.bufferΘ.r, 1:m)
    _matvec_mul!(predview, Lv, ws.state)
    yv = view(ws.ỹ, 1:m)
    @inbounds for i in 1:m
        yv[i] = data[observed[i], obs_col] - (predview[i] + μv[i])
    end

    # Stan's generated model never Cholesky-factorizes a covariance matrix
    # directly: every use passes through `makesym()`, which unconditionally
    # adds 1e-10 to the diagonal (see ctModelWriter.R, `out[coli,coli] =
    # mat[coli,coli] + 1e-10`). Two separate makesym() calls feed into the
    # innovation covariance: one ridges the predicted state covariance etacov
    # itself before it is projected through Jy (`ycov[o,o] =
    # quad_form_sym(makesym(etacov,verbose,1), Jy[o,]') + MANIFESTcov[o,o]`),
    # and a second ridges ycov again right before the Cholesky/division used
    # for the Kalman gain. Only replicating the second one (as a single ridge
    # on the innovation covariance) understates Stan's regularisation whenever
    # Jy isn't the identity, and for models with MANIFESTVAR fixed to exactly
    # zero and LAMBDA the identity, the innovation covariance is numerically
    # identical to the predicted state covariance P, which then sits exactly
    # on the PD/non-PD boundary whenever a trial DIFFUSION or T0VAR estimate
    # is small. Without both ridges, this Julia port hits that boundary far
    # more readily than Stan does, which both corrupts L-BFGS's line search
    # (trial points get rejected that Stan would have accepted) and forces
    # the optimizer to search a narrower, more boundary-constrained trust
    # region overall. The ridge on P_predict here is transient (mirroring
    # makesym()'s non-mutating copy): it is removed again below before
    # ws.P_predict is reused, unridged, for the Kalman gain and
    # covariance-update terms that Stan's `etacov * Jy[od,]'` also computes
    # from the unridged etacov.
    _copy_lower_to_upper!(ws.P_predict.data, ws.state_dim)
    _ridge_diagonal!(ws.P_predict.data, n, 1e-10)
    Hv = view(pars.Jy, observed, :)
    # PHt = P_{t|t-1} * Jy[observed,:]' (n x m); since P is symmetric,
    # transpose(PHt) is simultaneously Jy[observed,:] * P_{t|t-1}, so this one
    # buffer covers both orientations of that product.
    PHt = view(ws.K, :, 1:m)
    mul!(PHt, ws.P_predict.data, transpose(Hv))
    Sv = view(ws.S.UL.data, 1:m, 1:m)
    mul!(Sv, Hv, PHt)
    Rv = view(ws.bufferΘ.out, observed, observed)
    Sv .+= Rv
    _ridge_diagonal!(ws.P_predict.data, n, -1e-10)
    _symmetrize_and_ridge!(Sv, m)
    factor = cholesky!(Sv, check=false)
    issuccess(factor) || return nothing

    # Data generation, if asked for: draw this row's observation from its own
    # prior predictive and carry on as though it had been read. Taken here
    # rather than earlier because S is the innovation covariance the update is
    # about to use, already factorized -- so the drawn innovation is L z by
    # construction and cannot drift from the covariance the filter then
    # conditions on.
    generate === nothing ||
        _generate_row!(generate, ws, factor, yv, predview, μv, observed, obs_col)

    # State update: x_{t|t} = x_{t|t-1} + PHt * (S^{-1} * ỹ)
    row_sq = view(ws.bufferΘ.row_sq, 1:m)
    ldiv!(row_sq, factor, yv)
    mul!(ws.state, PHt, row_sq, one(eltype(ws.state)), one(eltype(ws.state)))

    # Covariance update (Joseph form):
    #   G_t = PHt * S^{-1}; A_t = I - G_t * Jy[observed,:]
    #   P_{t|t} = A_t * P_{t|t-1} * A_t' + G_t * Theta[observed,observed] * G_t'
    KRv = view(ws.KR, :, 1:m)
    copyto!(KRv, PHt)
    rdiv!(KRv, factor)
    mul!(ws.bufferQ.intermediate, KRv, Hv)
    oneT = one(eltype(ws.K))
    zeroT = zero(oneT)
    @inbounds for j in 1:n, i in 1:n
        ws.bufferQ.intermediate[i, j] = (i == j ? oneT : zeroT) - ws.bufferQ.intermediate[i, j]
    end
    mul!(ws.bufferQ.out, ws.bufferQ.intermediate, ws.P_predict.data)
    mul!(ws.P_update.data, ws.bufferQ.out, transpose(ws.bufferQ.intermediate))

    GR = view(ws.K, :, 1:m)  # PHt is no longer needed; reuse the same scratch for G_t * Theta
    mul!(GR, KRv, Rv)
    mul!(ws.bufferQ.out, GR, transpose(KRv))
    ws.P_update.data .+= ws.bufferQ.out
    _copy_lower_to_upper!(ws.P_update.data, ws.state_dim)
    return factor
end

"""
    _ekf_update_step!(ws, pars, data, obs_col)

Apply one EKF measurement update assuming every manifest row is observed.

A thin convenience wrapper around `_ekf_masked_update_step!` with
`observed = 1:manifest_dim`, kept for callers (e.g. `benchmark_setup.jl`)
that already know their data has no missingness and so have no need to
compute an observed-row count first.
"""
@inline function _ekf_update_step!(ws::ContinuousEKFWorkspace, pars, data::AbstractMatrix, obs_col::Int)
    return _ekf_masked_update_step!(ws, pars, data, obs_col, 1:_val(ws.manifest_dim))
end

@inline function _apply_td_impulse!(ws::ContinuousEKFWorkspace, pars, tdpreds)
    isempty(tdpreds) && return nothing
    _matvec_mul!(ws.bufferQ.r, pars.TDPREDEFFECT, tdpreds, ws.state_dim, Val(length(tdpreds)))
    ws.state .+= ws.bufferQ.r
    mul!(ws.bufferQ.intermediate, pars.Jtd, ws.P_predict.data)
    _mul_right_transpose!(ws.bufferQ.out, ws.bufferQ.intermediate, pars.Jtd,
        ws.state_dim, ws.state_dim, ws.state_dim)
    copyto!(ws.P_predict.data, ws.bufferQ.out)
    return nothing
end

# A finite-but-huge sentinel (e.g. -1e100) is invisible to the isfinite()
# checks in ctsem_optimize's fg!, and under ForwardDiff it silently
# contributes an exact-zero gradient for whatever subject/row triggered it
# (0 * finite == 0), rather than flagging the evaluation as invalid. That let
# Optim's L-BFGS accept trial points whose reported gradient had quietly
# dropped one subject's true partials, corrupting the line search/curvature
# updates without ever producing a NaN or Inf anywhere else. NaN poisons both
# the primal value and, via IEEE's 0 * NaN == NaN, every ForwardDiff partial
# once this subject's contribution is summed into the total objective, so the
# existing isfinite(value)/isfinite(gradient) guards in fg! correctly reject
# the whole trial point instead of silently accepting a corrupted gradient.
@inline _invalid_ekf_loglikelihood(ws::ContinuousEKFWorkspace) =
    -one(eltype(ws.state)) * NaN

"""
    _extended_kalman_filter_continuous!(ws, params, data, timesteps, sp)

Evaluate the continuous-time EKF log-likelihood using a preallocated workspace.

`params` contains the free parameter values, `sp` describes how to materialize
them, `data` is arranged as variables by time, and `timesteps` gives the time
for each column.
"""
function _extended_kalman_filter_continuous!(
    ws::ContinuousEKFWorkspace{T},
    params::AbstractVector{T},
    data::Matrix,
    timesteps::Vector,
    sp,
    tdpreds::AbstractMatrix=zeros(eltype(params), 0, size(data, 2)),
    tipreds::AbstractVector=eltype(params)[],
    subject::Integer=1,
    max_timestep::Real=Inf,
    trace=nothing,
    generate=nothing,
)::T where {T}
    # Materialize transformed/free/fixed values into the full parameter vector.
    _materialize_subject_values!(ws.subject_values, params, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    # `trace` is the adjoint tape (see adjoint_ekf.jl) or `nothing`. Every
    # `_record_*!`/`_begin_predict!` call below has a `::Nothing` method, so the
    # ordinary primal and ForwardDiff paths compile to exactly the code they
    # did before this hook existed. Keeping the recording inline in the one
    # real loop -- rather than maintaining a parallel traced copy of it -- is
    # what stops the reverse pass silently drifting away from the forward one.
    _record_subject_values!(trace, ws.subject_values)

    # Populate ComponentArray storage from transformed parameters.
    pars = ws.pars
    all_params = getdata(pars)
    all_params .= ws.all_params

    # Initial prior:
    #   x_{1|0} = T0MEANS
    #   P_{1|0} = T0VAR
    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.T0VAR, 0, ws.state_dim)
    copyto!(ws.P_predict.data, ws.bufferQ.out)
    _copy_lower_to_upper!(ws.P_predict.data, ws.state_dim)
    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0, ws.manifest_dim)
    copyto!(ws.state, pars.T0MEANS)
    _record_init!(trace, pars, _val(ws.state_dim))
    _record_theta!(trace, pars, _val(ws.manifest_dim))

    # Each row follows one contract: prediction, TD impulse, measurement.
    # The first row has no prediction interval but can still contain an impulse.
    first_context = CTSEMRowContext(ws.state, pars, view(tdpreds, :, 1), tipreds,
        timesteps[1], zero(eltype(params)), Int(subject), 1)
    _record_group!(trace, 2, ws.td_param_indices, all_params, first_context)
    apply_complex_transforms_at_indices!(all_params, ws.td_param_indices, sp.td_transforms, first_context)
    _record_td!(trace, ws, pars, first_context.tdpreds, _val(ws.state_dim))
    _apply_td_impulse!(ws, pars, first_context.tdpreds)
    update_context = CTSEMRowContext(ws.state, pars, first_context.tdpreds, tipreds,
        timesteps[1], zero(eltype(params)), Int(subject), 1)
    _record_group!(trace, 3, ws.update_param_indices, all_params, update_context)
    apply_complex_transforms_at_indices!(all_params, ws.update_param_indices, sp.update_transforms, update_context)
    log2π_const = log(2π)
    _record_row_prior!(trace, ws, pars, 1)
    ll = _ekf_update_observed!(ws, pars, data, 1, log2π_const, trace, generate)
    ll === nothing && return _invalid_ekf_loglikelihood(ws)
    generate === nothing || (generate.llrow[generate.offset+1] = ll)
    _record_row_update!(trace, ws, pars, sp, update_context, 1, ll)

    # Main EKF loop for t >= 2:
    #   (1) predict from t-1 to t using Δt
    #   (2) update with observation y_t
    #   (3) accumulate log p(y_t | y_{1:t-1})
    prev_timestep = timesteps[1]
    nsteps = min(length(timesteps), size(data, 2))
    @inbounds for t_idx in 2:nsteps
        curr_timestep = timesteps[t_idx]
        Δt = curr_timestep - prev_timestep

        # Match Stan's nonlinear integration contract: re-materialize the
        # local affine model at each bounded substep before prediction.
        n_substeps = max(1, ceil(Int, Δt / max_timestep))
        substep_dt = Δt / n_substeps
        @inbounds for substep in 1:n_substeps
            substep_time = prev_timestep + substep * substep_dt
            predict_context = CTSEMRowContext(ws.state, pars, view(tdpreds, :, t_idx), tipreds,
                substep_time, substep_dt, Int(subject), t_idx)
            _record_group!(trace, 1, ws.predict_param_indices, all_params, predict_context)
            apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices, sp.predict_transforms, predict_context)
            predict_snapshot = _begin_predict!(trace, ws, _val(ws.state_dim))
            _ekf_predict_step!(ws, pars, substep_dt)
            _record_predict!(trace, ws, pars, predict_snapshot, substep_dt, _val(ws.state_dim))
            _record_transition!(trace, ws, pars, t_idx, substep, n_substeps, Δt)
            # The next bounded step begins from this step's predicted covariance.
            copyto!(ws.P_update.data, ws.P_predict.data)
        end

        td_context = CTSEMRowContext(ws.state, pars, view(tdpreds, :, t_idx), tipreds,
            curr_timestep, Δt, Int(subject), t_idx)
        _record_group!(trace, 2, ws.td_param_indices, all_params, td_context)
        apply_complex_transforms_at_indices!(all_params, ws.td_param_indices, sp.td_transforms, td_context)
        _record_td!(trace, ws, pars, td_context.tdpreds, _val(ws.state_dim))
        _apply_td_impulse!(ws, pars, td_context.tdpreds)
        _record_td_transition!(trace, pars, t_idx, size(tdpreds, 1))

        measurement_context = CTSEMRowContext(ws.state, pars, td_context.tdpreds, tipreds,
            curr_timestep, Δt, Int(subject), t_idx)
        _record_group!(trace, 3, ws.update_param_indices, all_params, measurement_context)
        apply_complex_transforms_at_indices!(all_params, ws.update_param_indices, sp.update_transforms, measurement_context)
        ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0, ws.manifest_dim)
        _record_theta!(trace, pars, _val(ws.manifest_dim))
        _record_row_prior!(trace, ws, pars, t_idx)
        row_ll = _ekf_update_observed!(ws, pars, data, t_idx, log2π_const, trace, generate)
        row_ll === nothing && return _invalid_ekf_loglikelihood(ws)
        generate === nothing || (generate.llrow[generate.offset+t_idx] = row_ll)
        _record_row_update!(trace, ws, pars, sp, measurement_context, t_idx, row_ll)
        ll += row_ll

        prev_timestep = curr_timestep
    end
    return ll
end

"""
    _extended_kalman_filter_continuous(params, data, timesteps, sp)

Evaluate the continuous-time EKF log-likelihood with a freshly allocated
workspace.
"""
function _extended_kalman_filter_continuous(params, data::Matrix, timesteps::Vector, sp)
    _validate_continuous_ekf_inputs(timesteps, data)
    ws = _init_continuous_ekf_workspace(eltype(params), sp)
    return _extended_kalman_filter_continuous!(ws, params, data, timesteps, sp)
end

"""
    extended_kalman_log_likelihood_continuous(values, params, timesteps, data)

Return the continuous-time extended Kalman filter log-likelihood.

`values` contains free parameter values and `params` contains the parameter
metadata produced by `ekf_from_data_frame`.
"""
function extended_kalman_log_likelihood_continuous(values::Vector, params::EKFParameters, timesteps::Vector, data::Matrix)
    return _extended_kalman_filter_continuous(values, data, timesteps, params)
end

"""
    grad_log_likelihood_ekf_continuous(values, params, timesteps, data)

Return the ForwardDiff gradient of the continuous-time EKF log-likelihood with
respect to `values`.
"""
function grad_log_likelihood_ekf_continuous(values::Vector, params::EKFParameters, timesteps::Vector, data::Matrix)
    workspace = ContinuousEKFGradientWorkspace(values, params, timesteps, data)
    return grad_log_likelihood_ekf_continuous!(workspace.result, values, workspace)
end

"""
    grad_log_likelihood_ekf_continuous!(result, values, workspace)

Write the ForwardDiff gradient of the continuous-time EKF log-likelihood into
`result` using a reusable `ContinuousEKFGradientWorkspace`.
"""
function grad_log_likelihood_ekf_continuous!(
    result::AbstractVector,
    values::AbstractVector,
    workspace::ContinuousEKFGradientWorkspace,
)
    _check_gradient_workspace_args(result, values)
    ForwardDiff.gradient!(result, workspace.objective, values, workspace.cfg)
    return result
end

"""
    grad_log_likelihood_ekf_continuous!(values, workspace)

Write the gradient into the reusable result buffer stored by `workspace`.
"""
function grad_log_likelihood_ekf_continuous!(
    values::AbstractVector,
    workspace::ContinuousEKFGradientWorkspace,
)
    return grad_log_likelihood_ekf_continuous!(workspace.result, values, workspace)
end

"""
    res_and_grad_likelihood_ekf_continuous(values, params, timesteps, data)

Return the continuous-time EKF log-likelihood and its gradient.

The result is a named tuple `(value = ..., gradient = ...)`.
"""
function res_and_grad_likelihood_ekf_continuous(values::Vector, params::EKFParameters, timesteps::Vector, data::Matrix)
    workspace = ContinuousEKFGradientWorkspace(values, params, timesteps, data)
    return res_and_grad_likelihood_ekf_continuous!(workspace, values)
end

"""
    res_and_grad_likelihood_ekf_continuous!(workspace, values)

Return the continuous-time EKF log-likelihood and gradient using reusable
ForwardDiff and EKF workspaces.
"""
function res_and_grad_likelihood_ekf_continuous!(
    workspace::ContinuousEKFGradientWorkspace,
    values::AbstractVector,
)
    _check_gradient_workspace_args(workspace.result, values)
    ForwardDiff.gradient!(
        workspace.diff_result,
        workspace.objective,
        values,
        workspace.cfg,
    )
    return (
        value = DiffResults.value(workspace.diff_result),
        gradient = DiffResults.gradient(workspace.diff_result),
    )
end
