using LinearAlgebra

################################################################################
# Per-row Kalman output
################################################################################
#
# What `ctKalman()`/`ctPredict()` need: prior, filtered and smoothed states and
# observations for every row, and -- for a model with random effects -- each
# subject's own parameters.
#
# All of it is a byproduct of the pass the likelihood already makes, so this
# rides on that pass rather than adding a second one. The mechanism is the
# `trace` argument the filter already carries for the adjoint: a trace is
# whatever implements the `_record_*!` hooks, the adjoint tape implements the
# ones it needs, and `CTSEMKalmanTrace` implements the ones it needs. Dispatch
# rather than branching means the primal and adjoint paths compile to exactly
# the code they did before this file existed -- a `nothing` trace, and now a
# tape, resolve every hook to a no-op inlined away.
#
# Three hooks are new (`_record_row_prior!`, `_record_row_update!`,
# `_record_transition!`), because the adjoint has no use for the quantities they
# capture and so never needed a call site there. They sit in the filter loop
# rather than inside `_ekf_update_observed!` so that they see the row index and
# fire even on a fully missing row, which the update returns early from.
#
# Two details follow Stan rather than improving on it, because the point is that
# the same model predicts the same way whichever backend ran it:
#
#   * the manifest mean and covariance are recorded over *all* manifest
#     variables, observed or not, and are not re-evaluated at the updated state
#     (`LAMBDA`/`Jy` stay as the predicted state left them);
#   * with bounded substeps, the interval transition handed to the smoother is
#     `exp(JAx * Δt)` from the last substep rather than the product of the
#     substep transitions. The product is the exact Jacobian of what the filter
#     actually did; Stan calls its version an approximation, and switching is a
#     one-line change if Stan does.

export CTSEMKalmanTrace, ctsem_kalman

const _CTSEM_KALMAN_PRIOR = 1
const _CTSEM_KALMAN_UPD = 2
const _CTSEM_KALMAN_SMOOTH = 3

"""
    CTSEMKalmanTrace(nlatent, nmanifest, nrows, nsubjects, nall)

Sink for per-row Kalman output, indexed by *global* data row so that one
allocation covers the whole dataset and each subject writes its own slice.

The leading dimension of `eta`, `etacov`, `y` and `ycov` is 1 = prior, 2 =
updated, 3 = smoothed, matching Stan's `etaa`/`ycova` layout so the R side needs
no second unpacking path.
"""
mutable struct CTSEMKalmanTrace{T}
    nlatent::Int
    nmanifest::Int
    nrows::Int
    eta::Array{T,3}
    etacov::Array{T,4}
    y::Array{T,3}
    ycov::Array{T,4}
    llrow::Vector{T}
    transition::Vector{Matrix{T}}
    Jy::Vector{Matrix{T}}
    subject::Vector{Int}
    subject_params::Vector{Vector{T}}
    subject_t0::Vector{Vector{T}}
    # Set by the driver before each subject's forward pass.
    offset::Int
    current_subject::Int
end

function CTSEMKalmanTrace(::Type{T}, nlatent::Int, nmanifest::Int, nrows::Int,
    nsubjects::Int, nall::Int) where {T}
    return CTSEMKalmanTrace{T}(nlatent, nmanifest, nrows,
        zeros(T, 3, nrows, nlatent), zeros(T, 3, nrows, nlatent, nlatent),
        zeros(T, 3, nrows, nmanifest), zeros(T, 3, nrows, nmanifest, nmanifest),
        zeros(T, nrows), [Matrix{T}(I, nlatent, nlatent) for _ in 1:nrows],
        [zeros(T, nmanifest, nlatent) for _ in 1:nrows], zeros(Int, nrows),
        [zeros(T, nall) for _ in 1:nsubjects], [zeros(T, nlatent) for _ in 1:nsubjects],
        0, 0)
end

# No-ops for the trace types that do not want these, so the primal and adjoint
# paths are untouched by their presence.
@inline _record_row_prior!(::Nothing, args...) = nothing
@inline _record_row_update!(::Nothing, args...) = nothing
@inline _record_transition!(::Nothing, args...) = nothing
@inline _record_row_prior!(::CTSEMAdjointTape, args...) = nothing
@inline _record_row_update!(::CTSEMAdjointTape, args...) = nothing
@inline _record_transition!(::CTSEMAdjointTape, args...) = nothing

# ... and no-ops for every adjoint hook, so the same filter runs with a Kalman
# trace in place of a tape.
@inline _record_subject_values!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_init!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_theta!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_group!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_td!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_update!(::CTSEMKalmanTrace, args...) = nothing
@inline _begin_predict!(::CTSEMKalmanTrace, args...) = nothing
@inline _record_predict!(::CTSEMKalmanTrace, args...) = nothing

"""
    _kalman_store!(trace, kind, row, ws, pars, state, statecov)

Write one (state, covariance) pair and the manifest quantities implied by it.
"""
function _kalman_store!(trace::CTSEMKalmanTrace, kind::Int, row::Int, ws, pars,
    state, statecov)
    n = trace.nlatent
    m = trace.nmanifest
    Θ = ws.bufferΘ.out
    @inbounds for i in 1:n
        trace.eta[kind, row, i] = state[i]
        for j in 1:n
            trace.etacov[kind, row, i, j] = statecov[i, j]
        end
    end
    @inbounds for i in 1:m
        value = pars.MANIFESTMEANS[i, 1]
        for j in 1:n
            value += pars.LAMBDA[i, j] * state[j]
        end
        trace.y[kind, row, i] = value
    end
    @inbounds for i in 1:m
        for j in 1:m
            value = Θ[i, j]
            for a in 1:n
                inner = zero(eltype(Θ))
                for b in 1:n
                    inner += statecov[a, b] * pars.Jy[j, b]
                end
                value += pars.Jy[i, a] * inner
            end
            trace.ycov[kind, row, i, j] = value
        end
    end
    return nothing
end

function _record_row_prior!(trace::CTSEMKalmanTrace, ws, pars, row::Int)
    r = trace.offset + row
    _kalman_store!(trace, _CTSEM_KALMAN_PRIOR, r, ws, pars, ws.state, ws.P_predict.data)
    @inbounds copyto!(trace.Jy[r], pars.Jy)
    trace.subject[r] = trace.current_subject
    return nothing
end

function _record_row_update!(trace::CTSEMKalmanTrace, ws, pars, row::Int, row_ll)
    r = trace.offset + row
    _kalman_store!(trace, _CTSEM_KALMAN_UPD, r, ws, pars, ws.state, ws.P_update.data)
    trace.llrow[r] = row_ll === nothing ? NaN : row_ll
    return nothing
end

function _record_transition!(trace::CTSEMKalmanTrace, ws, pars, row::Int, substep::Int,
    n_substeps::Int, Δt)
    r = trace.offset + row
    if n_substeps == 1
        @inbounds copyto!(trace.transition[r], ws.discrete_ca.dDRIFT)
    elseif substep == n_substeps
        # See the header note: Stan recomputes the whole-interval exponential
        # from the last substep's Jacobian rather than composing the substeps.
        @inbounds copyto!(trace.transition[r], exp(Matrix(pars.JAx) .* Δt))
    end
    return nothing
end

"""
    _kalman_smooth!(trace, first, nobs)

Backward RTS pass over one subject's rows, already filtered.

Fixed-interval smoothing needs the last row before it can produce the first, so
it cannot ride inside the forward loop; it is a genuinely separate pass and is
written as one rather than disguised as part of the first.
"""
function _kalman_smooth!(trace::CTSEMKalmanTrace{T}, first::Int, nobs::Int) where {T}
    n = trace.nlatent
    m = trace.nmanifest
    last = first + nobs - 1
    P = Matrix{T}(undef, n, n)
    for r in last:-1:first
        if r == last
            @inbounds trace.eta[_CTSEM_KALMAN_SMOOTH, r, :] .= trace.eta[_CTSEM_KALMAN_UPD, r, :]
            @inbounds trace.etacov[_CTSEM_KALMAN_SMOOTH, r, :, :] .=
                trace.etacov[_CTSEM_KALMAN_UPD, r, :, :]
            @inbounds trace.y[_CTSEM_KALMAN_SMOOTH, r, :] .= trace.y[_CTSEM_KALMAN_UPD, r, :]
            @inbounds trace.ycov[_CTSEM_KALMAN_SMOOTH, r, :, :] .=
                trace.ycov[_CTSEM_KALMAN_UPD, r, :, :]
            continue
        end
        # gain = P_upd[r] A[r+1]' inv(P_prior[r+1]), with Stan's makesym() ridge
        # on the matrix being inverted, solved rather than inverted.
        @inbounds P .= trace.etacov[_CTSEM_KALMAN_PRIOR, r+1, :, :]
        P .= (P .+ P') ./ 2
        @inbounds for i in 1:n
            P[i, i] += 1e-10
        end
        Pupd = @inbounds trace.etacov[_CTSEM_KALMAN_UPD, r, :, :]
        cross = Pupd * trace.transition[r+1]'
        gain = (Symmetric(P) \ cross')'

        etanext_s = @inbounds trace.eta[_CTSEM_KALMAN_SMOOTH, r+1, :]
        etanext_p = @inbounds trace.eta[_CTSEM_KALMAN_PRIOR, r+1, :]
        etaupd = @inbounds trace.eta[_CTSEM_KALMAN_UPD, r, :]
        etasm = etaupd .+ gain * (etanext_s .- etanext_p)

        covnext_s = @inbounds trace.etacov[_CTSEM_KALMAN_SMOOTH, r+1, :, :]
        covnext_p = @inbounds trace.etacov[_CTSEM_KALMAN_PRIOR, r+1, :, :]
        covsm = Pupd .+ gain * (covnext_s .- covnext_p) * gain'

        Jyr = trace.Jy[r]
        yupd = @inbounds trace.y[_CTSEM_KALMAN_UPD, r, :]
        ycovupd = @inbounds trace.ycov[_CTSEM_KALMAN_UPD, r, :, :]
        ysm = yupd .+ Jyr * (etasm .- etaupd)
        ycovsm = ycovupd .+ Jyr * (covsm .- Pupd) * Jyr'

        @inbounds trace.eta[_CTSEM_KALMAN_SMOOTH, r, :] .= etasm
        @inbounds trace.etacov[_CTSEM_KALMAN_SMOOTH, r, :, :] .= covsm
        @inbounds trace.y[_CTSEM_KALMAN_SMOOTH, r, :] .= ysm
        @inbounds trace.ycov[_CTSEM_KALMAN_SMOOTH, r, :, :] .= ycovsm
    end
    return nothing
end

"""
    ctsem_kalman(objective, values; subject_matrices=true)

Prior, filtered and smoothed states and observations for every data row, plus
each subject's own model matrices.

`values` is one raw parameter vector. The returned named tuple carries `eta`,
`etacov`, `y`, `ycov` with leading dimension 1 = prior, 2 = updated, 3 =
smoothed; `llrow`; `subject`; `subject_loglik`; and, unless
`subject_matrices=false`, a `size`-by-`nsubjects` matrix packed in the layout
`ctsem_parameter_layout` describes.

A subject's matrices are its parameter vector as of its last row with T0MEANS
replaced by its smoothed initial state. That is where individual differences
come out: ctsem carries an individually varying parameter as an augmented
latent state with no drift and no diffusion, so the smoothed t0 estimate of that
state is the subject's value for it.
"""
function ctsem_kalman(objective::CTSEMObjective, values::AbstractVector;
    subject_matrices::Bool=true)

    sp = objective.params
    raw = Vector{Float64}(values)
    ws = _init_continuous_ekf_workspace(Float64, sp)
    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)
    nall = length(sp.mutables)
    subjects = objective.subject_objectives
    nrows = sum(size(sub.data, 2) for sub in subjects)

    trace = CTSEMKalmanTrace(Float64, n, m, nrows, length(subjects), nall)
    loglik = zeros(Float64, length(subjects))

    offset = 0
    for (i, sub) in enumerate(subjects)
        nobs = size(sub.data, 2)
        trace.offset = offset
        trace.current_subject = i
        value = _extended_kalman_filter_continuous!(ws, raw, sub.data,
            collect(sub.timesteps), sp, sub.tdpreds, sub.tipreds, i, sub.max_timestep,
            trace)
        loglik[i] = value
        if isfinite(value)
            _kalman_smooth!(trace, offset + 1, nobs)
            # getdata(ws.pars), not ws.all_params: the filter copies the
            # materialized vector into the ComponentArray once and every
            # state-dependent transform writes there afterwards, so
            # ws.all_params still holds the pre-transform values and its
            # state-dependent cells (JAx, Jy, ...) were never written at all.
            copyto!(trace.subject_params[i], getdata(ws.pars))
            @inbounds trace.subject_t0[i] .= trace.eta[_CTSEM_KALMAN_SMOOTH, offset + 1, :]
        end
        offset += nobs
    end

    result = (eta=trace.eta, etacov=trace.etacov, y=trace.y, ycov=trace.ycov,
        llrow=trace.llrow, subject=trace.subject, subject_loglik=loglik)
    subject_matrices || return result
    return merge(result, (subject_matrices=_ctsem_subject_matrices(objective, trace),))
end

"""
    _ctsem_subject_matrices(objective, trace)

Pack every subject's matrices into the flat layout `ctsem_parameter_layout`
describes, one column per subject.
"""
function _ctsem_subject_matrices(objective::CTSEMObjective, trace::CTSEMKalmanTrace)
    sp = objective.params
    layout = ctsem_parameter_layout(objective)
    nsubjects = length(trace.subject_params)
    out = zeros(Float64, layout.size, nsubjects)
    pars = ComponentVector(zeros(Float64, length(sp.mutables)), sp.parameter_axis)
    for i in 1:nsubjects
        copyto!(getdata(pars), trace.subject_params[i])
        @inbounds for k in 1:trace.nlatent
            pars.T0MEANS[k, 1] = trace.subject_t0[i][k]
        end
        _ctsem_pack_matrices!(view(out, :, i), pars, sp, layout)
    end
    return out
end
