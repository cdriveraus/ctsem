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
# Manifest quantities are recorded over *all* manifest variables, observed or
# not, as Stan does. Three other things Stan does approximately are done exactly
# here; none touches the likelihood or the gradient, because all three affect
# only what is reported:
#
#   * the measurement model is re-evaluated at the *updated* state before the
#     filtered observation estimate is recorded. Stan applies LAMBDA and
#     MANIFESTMEANS as the prior state left them, and its own comment on that
#     block says it "could be improved by recomputing all state dependent pars,
#     error covariances etc. at each step". This matters in ordinary models, not
#     exotic ones: ctsem represents an individually varying parameter as an
#     augmented latent state and MANIFESTMEANS is individually varying by
#     default, so the measurement intercept *is* a state whose pre-update value
#     Stan reports;
#   * with bounded substeps the interval transition is the product of the
#     substep transitions -- the derivative of what the filter actually
#     computed -- rather than the whole-interval exponential recomputed from the
#     last substep;
#   * a TD impulse maps the covariance through Jtd and so belongs in the
#     interval transition the smoother uses. Stan saves only the exponential.

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
    # Scratch for re-evaluating the measurement model at the updated state
    # without leaking anything back into the filter.
    param_snapshot::Vector{T}
    theta_snapshot::Matrix{T}
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
        zeros(T, nall), zeros(T, nmanifest, nmanifest), 0, 0)
end

# No-ops for the trace types that do not want these, so the primal and adjoint
# paths are untouched by their presence.
@inline _record_row_prior!(::Nothing, args...) = nothing
@inline _record_row_update!(::Nothing, args...) = nothing
@inline _record_transition!(::Nothing, args...) = nothing
@inline _record_td_transition!(::Nothing, args...) = nothing
@inline _record_row_prior!(::CTSEMAdjointTape, args...) = nothing
@inline _record_row_update!(::CTSEMAdjointTape, args...) = nothing
@inline _record_transition!(::CTSEMAdjointTape, args...) = nothing
@inline _record_td_transition!(::CTSEMAdjointTape, args...) = nothing

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
    trace.subject[r] = trace.current_subject
    return nothing
end

"""
    _record_row_update!(trace, ws, pars, sp, context, row, row_ll)

Record the filtered estimates, with the measurement model re-evaluated at the
updated state.

`context` already refers to `ws.state`, which the measurement update mutated in
place, so re-running the transform groups against it evaluates them where the
filter now is rather than where it was before the observation arrived.

The parameter vector and the manifest covariance are restored afterwards, so
this cannot reach the likelihood: the filter continues from exactly where it
was.
"""
function _record_row_update!(trace::CTSEMKalmanTrace, ws, pars, sp, context,
    row::Int, row_ll)
    r = trace.offset + row
    all_params = getdata(pars)
    copyto!(trace.param_snapshot, all_params)
    copyto!(trace.theta_snapshot, ws.bufferΘ.out)

    apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
        sp.predict_transforms, context)   # PARS the measurement cells may read
    apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
        sp.update_transforms, context)
    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0, ws.manifest_dim)

    _kalman_store!(trace, _CTSEM_KALMAN_UPD, r, ws, pars, ws.state, ws.P_update.data)
    # Jy about the updated state: that is where the smoother takes its
    # measurement correction.
    @inbounds copyto!(trace.Jy[r], pars.Jy)
    trace.llrow[r] = row_ll === nothing ? NaN : row_ll

    copyto!(all_params, trace.param_snapshot)
    copyto!(ws.bufferΘ.out, trace.theta_snapshot)
    return nothing
end

# The interval transition, composed as the filter composed it: the substeps
# propagate x <- A_s x + b_s, so the interval Jacobian is A_S ... A_1. For a
# state-independent JAx this equals the whole-interval exponential exactly; for
# a state-dependent one it is the derivative of what was actually computed.
function _record_transition!(trace::CTSEMKalmanTrace, ws, pars, row::Int, substep::Int,
    n_substeps::Int, Δt)
    r = trace.offset + row
    @inbounds if substep == 1
        copyto!(trace.transition[r], ws.discrete_ca.dDRIFT)
    else
        trace.transition[r] = ws.discrete_ca.dDRIFT * trace.transition[r]
    end
    return nothing
end

# A TD impulse sits between the previous row's posterior and this row's prior and
# maps the covariance through Jtd, so it is part of the interval transition.
function _record_td_transition!(trace::CTSEMKalmanTrace, pars, row::Int, ntdpred::Int)
    ntdpred == 0 && return nothing
    r = trace.offset + row
    @inbounds trace.transition[r] = Matrix(pars.Jtd) * trace.transition[r]
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

        # The measurement correction is taken about the updated state, where Jy
        # was recorded and y_upd evaluated. For a measurement equation linear in
        # the augmented state -- every intoverpop model's -- the two together are
        # exact: y_upd = Jy x_upd, so y_sm = Jy x_sm.
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

`values` is one raw parameter vector, used for every subject -- or a matrix
whose row `i` is subject `i`'s own. The second form is what a Laplace fit needs,
where each subject is filtered at its own realized parameters rather than at a
shared population vector. The returned named tuple carries `eta`,
`etacov`, `y`, `ycov` with leading dimension 1 = prior, 2 = updated, 3 =
smoothed; `llrow`; `subject`; `subject_loglik`; `transition` (the interval
Jacobian that reached each row, `nrows` by `n` by `n`); and, unless
`subject_matrices=false`, a `size`-by-`nsubjects` matrix packed in the layout
`ctsem_parameter_layout` describes.

A subject's matrices are its parameter vector as of its last row with T0MEANS
replaced by its smoothed initial state. That is where individual differences
come out: ctsem carries an individually varying parameter as an augmented
latent state with no drift and no diffusion, so the smoothed t0 estimate of that
state is the subject's value for it.
"""
function ctsem_kalman(objective::CTSEMObjective, values::AbstractVecOrMat;
    subject_matrices::Bool=true)

    sp = objective.params
    persubject = values isa AbstractMatrix
    persubject && size(values, 1) == length(objective.subject_objectives) ||
        persubject && throw(DimensionMismatch(
            "one row of parameters per subject is required"))
    raw = persubject ? Vector{Float64}(view(values, 1, :)) : Vector{Float64}(values)
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
        persubject && copyto!(raw, view(values, i, :))
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

    # The interval transitions come back as an array rather than a vector of
    # matrices, both because the R bridge marshals a list element by element and
    # because this is the quantity the smoother is most easily got wrong in --
    # exposing it makes it checkable from outside.
    transitions = zeros(Float64, nrows, n, n)
    for r in 1:nrows, i in 1:n, j in 1:n
        @inbounds transitions[r, i, j] = trace.transition[r][i, j]
    end
    result = (eta=trace.eta, etacov=trace.etacov, y=trace.y, ycov=trace.ycov,
        llrow=trace.llrow, subject=trace.subject, subject_loglik=loglik,
        transition=transitions)
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


################################################################################
# Posterior-predictive data generation
################################################################################
#
# The same forward pass, drawing each row's observation from its own prior
# predictive instead of reading it. Unlike the recorders above this is not a
# passive observer -- it changes what the filter consumes, so the state carried
# forward is conditioned on the drawn data rather than the real data, which is
# exactly what makes the result a draw from the model.
#
# It rides on the filter for the same reason the recorders do: a separate
# simulator would have to re-derive the prior predictive at every row and would
# be free to disagree with the filter about it. Here the draw is taken from the
# already-factorized innovation covariance the update is about to use, so the
# two cannot come apart.
#
# The standard normal draws come from the caller, not from an engine RNG, so
# that a seed set in R gives the same data whichever engine ran it.

export CTSEMGenerateSpec, ctsem_generate

"""
    CTSEMGenerateSpec(base, out, llrow)

`base` and `out` are `nmanifest` by `nrows`, in the same layout as the observed
data. `offset` is set by the driver before each subject's pass.

Missing entries are never generated: the generated dataset keeps the original
missingness, because the point of it is comparison against the observations that
are actually there.
"""
mutable struct CTSEMGenerateSpec
    base::Matrix{Float64}
    out::Matrix{Float64}
    llrow::Vector{Float64}
    offset::Int
end

CTSEMGenerateSpec(base, out, llrow) = CTSEMGenerateSpec(base, out, llrow, 0)

function _generate_row!(gen::CTSEMGenerateSpec, ws, factor, yv, predview, μv,
    observed::AbstractVector{Int}, obs_col::Int)
    r = gen.offset + obs_col
    m = length(observed)
    z = view(ws.bufferΘ.row_sq, 1:m)
    @inbounds for i in 1:m
        z[i] = gen.base[observed[i], r]
    end
    # The innovation of a draw from N(mean, S) is exactly L z.
    mul!(yv, factor.L, z)
    @inbounds for i in 1:m
        gen.out[observed[i], r] = predview[i] + μv[i] + yv[i]
    end
    return nothing
end

"""
    ctsem_generate(objective, values, base)

One posterior-predictive dataset.

`base` is `nmanifest` by `nrows` standard normals. Returns `Y` (the generated
observations, `NaN` wherever the original data was missing), `llrow` (each row's
log likelihood *of the generated data*) and `subject_loglik`.
"""
function ctsem_generate(objective::CTSEMObjective, values::AbstractVector,
    base::AbstractMatrix)

    sp = objective.params
    raw = Vector{Float64}(values)
    ws = _init_continuous_ekf_workspace(Float64, sp)
    m = _val(ws.manifest_dim)
    subjects = objective.subject_objectives
    nrows = sum(size(sub.data, 2) for sub in subjects)
    size(base) == (m, nrows) ||
        throw(DimensionMismatch("base must be $(m) by $(nrows)"))

    generate = CTSEMGenerateSpec(Matrix{Float64}(base), fill(NaN, m, nrows),
        zeros(Float64, nrows))
    loglik = zeros(Float64, length(subjects))

    offset = 0
    for (i, sub) in enumerate(subjects)
        generate.offset = offset
        loglik[i] = _extended_kalman_filter_continuous!(ws, raw, sub.data,
            collect(sub.timesteps), sp, sub.tdpreds, sub.tipreds, i, sub.max_timestep,
            nothing, generate)
        offset += size(sub.data, 2)
    end
    return (Y=generate.out, llrow=generate.llrow, subject_loglik=loglik)
end
