using LinearAlgebra
using ComponentArrays

################################################################################
# Reverse-mode EKF
################################################################################
#
# The forward filter (`_extended_kalman_filter_continuous!`) optionally records
# a *tape* while it runs; this file replays that tape backwards.
#
# The tape is deliberately an ordered program of small typed records rather
# than a fixed-shape array-of-checkpoints. The forward pass is not a fixed
# sequence: the number of prediction substeps depends on `max_timestep` and the
# gap between observations, TD impulses only happen when there are TD
# predictors, state-dependent transform groups only exist for nonlinear models,
# and a fully missing row skips the measurement update entirely. Recording the
# order the forward pass actually took, rather than assuming a canonical one,
# is what lets one reverse pass cover every case the roadmap requires --
# linear, state/PARS-expression, partially missing, fully missing, and
# unequal-length multi-subject models -- without special-casing any of them.
#
# Conventions used throughout the reverse pass:
#
#   * `x̄` and `P̄` are the running cotangents of "the state and covariance at
#     this point in the forward pass". Every reverse step consumes them as the
#     cotangent of its outputs and leaves them as the cotangent of its inputs.
#     Identity copies in the primal (`P_update .= P_predict` between substeps,
#     and on a fully missing row) therefore need no reverse code at all.
#   * `θ̄` is the cotangent of `all_params`, in exactly the layout `ws.pars`
#     views. Matrix cotangents are accumulated into it through a
#     `ComponentVector` view, so `θ̄ca.DRIFT` and friends name the right cells.
#   * `Θ̄` is the running cotangent of the *manifest covariance* Θ, which is
#     accumulated by measurement updates and consumed by the `:theta` tape
#     entry at the point where the forward pass built Θ from MANIFESTVAR. It is
#     kept separate from `θ̄` precisely because that build happens at a
#     different point in the row for the first row than for later rows.
#   * Covariance cotangents are symmetrised before use. The primal maintains
#     symmetry explicitly (`_copy_lower_to_upper!`, `_symmetrize_and_ridge!`),
#     and those operations are the identity on symmetric perturbations, so
#     keeping the cotangents symmetric is what makes it valid to treat them as
#     such.
#   * The `+1e-10` ridges Stan applies (and this port matches) are affine, so
#     they are invisible to the reverse pass except through the *value* at
#     which downstream derivatives are evaluated -- which is why the reverse
#     recomputes `Pr = P + 1e-10 I` rather than reusing the unridged `P`.

const _CTSEM_RIDGE = 1e-10

################################################################################
# Tape records
################################################################################

"""One prediction substep: everything the reverse pass needs to undo it."""
struct CTSEMPredictRecord{T}
    state_in::Vector{T}      # state entering the substep
    P_in::Matrix{T}          # posterior covariance entering the substep (unridged)
    A::Matrix{T}             # eJAx = exp(JAx * dt), also the discrete drift
    JAx::Matrix{T}
    DRIFT::Matrix{T}
    DIFFUSION::Matrix{T}     # raw SD/correlation-sqrt parameters
    Xlyap::Matrix{T}         # Lyapunov solution on the dynamic block
    affine::Vector{T}        # `r` in the derivation: the local affine correction
    dINT_dynamic::Vector{T}  # solved discrete intercept on the dynamic block
    dt::T
end

"""One TD-predictor impulse."""
struct CTSEMTDRecord{T}
    P_in::Matrix{T}
    Jtd::Matrix{T}
    tdpreds::Vector{T}
end

"""
One measurement update.

Stores only inputs, not intermediates: the innovation, gain and Cholesky
factor are cheap to recompute from these and doing so keeps the primal's
`_ekf_masked_update_step!` free of tracing code (it reuses `ws.K` as scratch
partway through, so its intermediates are not all live at the end anyway).
"""
struct CTSEMUpdateRecord{T}
    observed::Vector{Int}
    state_in::Vector{T}
    P_in::Matrix{T}          # prior covariance, unridged
    Lambda::Matrix{T}        # LAMBDA[observed, :]
    manifestmeans::Vector{T} # MANIFESTMEANS[observed]
    H::Matrix{T}             # Jy[observed, :]
    R::Matrix{T}             # Θ[observed, observed]
    y::Vector{T}             # observed data for this row
end

"""
One application of a group of state-dependent transforms.

`params_before` holds the parameter values from *before* the group ran, and
both halves of that matter.

**Before, not after.** Within a group the transforms are applied in ascending
flattened-index order into the same buffer they read from, so a transform can
read a cell that a later transform in the same group is about to overwrite --
and when it does, it sees the value from the *previous* row, not this one.
(`PARS` sorts last in the parameter axis, so a state-dependent `DRIFT`
expression referencing a state-dependent `PARS` cell is exactly this case.)
Recording the pre-group state lets the reverse pass reconstruct, for each
transform, the values that transform actually saw; recording the post-group
state would silently evaluate some derivatives at the wrong point.

**Only the relevant cells, not the whole vector.** A group's transforms read
only the indices in their recorded read sets and write only their own
`indices`; every other entry of `all_params` is dead as far as this group is
concerned. Storing the lot cost ~23 KB per group per row on a 20-latent model
and made state-dependent models allocate ~50 MB per gradient. `params_before`
therefore holds just the values at `tape.group_relevant[group]`, in that
order.
"""
struct CTSEMGroupRecord{T}
    group::Int               # 1 = predict, 2 = td, 3 = update
    params_before::Vector{T} # values at `group_relevant[group]`, pre-group
    state::Vector{T}
    tdpreds::Vector{T}
    time::T
    dt::T
    row::Int
end

"""Construction of Θ from MANIFESTVAR."""
struct CTSEMThetaRecord{T}
    MANIFESTVAR::Matrix{T}
end

"""The `t = 1` prior: state from T0MEANS, covariance from T0VAR."""
struct CTSEMInitRecord{T}
    T0VAR::Matrix{T}
end

"""
    CTSEMAdjointTape

Ordered record of one subject's forward pass.

`program` lists `(kind, index)` pairs in forward order; the reverse pass walks
it backwards and dispatches on `kind` into the matching record vector.
"""
mutable struct CTSEMAdjointTape{T}
    program::Vector{Tuple{Symbol,Int}}
    predicts::Vector{CTSEMPredictRecord{T}}
    tds::Vector{CTSEMTDRecord{T}}
    updates::Vector{CTSEMUpdateRecord{T}}
    groups::Vector{CTSEMGroupRecord{T}}
    thetas::Vector{CTSEMThetaRecord{T}}
    inits::Vector{CTSEMInitRecord{T}}
    subject_values::Vector{T}
    # For each transform group (1 = predict, 2 = td, 3 = update), the
    # `all_params` indices that group can read or write: the union of its
    # transforms' recorded read sets and their written indices. This is what
    # a group record snapshots instead of the whole parameter vector.
    group_relevant::Vector{Vector{Int}}
end

CTSEMAdjointTape(::Type{T}, group_relevant=[Int[], Int[], Int[]]) where {T} =
    CTSEMAdjointTape{T}(
        Tuple{Symbol,Int}[], CTSEMPredictRecord{T}[], CTSEMTDRecord{T}[],
        CTSEMUpdateRecord{T}[], CTSEMGroupRecord{T}[], CTSEMThetaRecord{T}[],
        CTSEMInitRecord{T}[], T[], group_relevant)

function _tape_reset!(tape::CTSEMAdjointTape)
    empty!(tape.program); empty!(tape.predicts); empty!(tape.tds)
    empty!(tape.updates); empty!(tape.groups); empty!(tape.thetas)
    empty!(tape.inits)
    return tape
end

@inline _tape_push!(tape::CTSEMAdjointTape, kind::Symbol, index::Int) =
    push!(tape.program, (kind, index))

################################################################################
# Recording hooks
################################################################################
#
# Every hook has a `::Nothing` method so the untraced primal pays nothing: with
# `trace === nothing` the calls have no arguments the compiler cannot see
# through and disappear entirely. This is why the forward filter has exactly
# one loop rather than a primal loop and a parallel traced copy -- the reverse
# pass cannot drift out of step with a forward pass it is literally recorded
# from.

@inline _record_subject_values!(::Nothing, args...) = nothing
@inline _record_init!(::Nothing, args...) = nothing
@inline _record_theta!(::Nothing, args...) = nothing
@inline _record_group!(::Nothing, args...) = nothing
@inline _record_td!(::Nothing, args...) = nothing
@inline _record_update!(::Nothing, args...) = nothing
@inline _begin_predict!(::Nothing, args...) = nothing
@inline _record_predict!(::Nothing, args...) = nothing

function _record_subject_values!(tape::CTSEMAdjointTape, subject_values)
    resize!(tape.subject_values, length(subject_values))
    copyto!(tape.subject_values, subject_values)
    return nothing
end

function _record_init!(tape::CTSEMAdjointTape, pars, n::Int)
    push!(tape.inits, CTSEMInitRecord(Matrix(pars.T0VAR[1:n, 1:n])))
    _tape_push!(tape, :init, length(tape.inits))
    return nothing
end

function _record_theta!(tape::CTSEMAdjointTape, pars, m::Int)
    push!(tape.thetas, CTSEMThetaRecord(Matrix(pars.MANIFESTVAR[1:m, 1:m])))
    _tape_push!(tape, :theta, length(tape.thetas))
    return nothing
end

"""
Record a transform group. **Must be called before the group is applied** --
see `CTSEMGroupRecord` for why the pre-group parameter values are the ones
the reverse pass needs.
"""
function _record_group!(tape::CTSEMAdjointTape{T}, group::Int,
    indices::AbstractVector{Int}, all_params, ctx::CTSEMRowContext) where {T}
    # Linear models have no state-dependent cells at all; skipping the
    # snapshot there keeps the tape out of the common case entirely.
    isempty(indices) && return nothing
    relevant = tape.group_relevant[group]
    compact = Vector{T}(undef, length(relevant))
    @inbounds for j in eachindex(relevant)
        compact[j] = all_params[relevant[j]]
    end
    push!(tape.groups, CTSEMGroupRecord(group, compact, collect(vec(ctx.state)),
        collect(ctx.tdpreds), float(ctx.time), float(ctx.dt), ctx.row))
    _tape_push!(tape, :group, length(tape.groups))
    return nothing
end

function _record_td!(tape::CTSEMAdjointTape, ws, pars, tdpreds, n::Int)
    # Mirrors `_apply_td_impulse!`'s own early return: with no TD predictors
    # the impulse is a no-op and contributes nothing to reverse.
    isempty(tdpreds) && return nothing
    push!(tape.tds, CTSEMTDRecord(Matrix(ws.P_predict.data[1:n, 1:n]),
        Matrix(pars.Jtd[1:n, 1:n]), collect(tdpreds)))
    _tape_push!(tape, :td, length(tape.tds))
    return nothing
end

function _record_update!(tape::CTSEMAdjointTape, ws, pars, data, obs_col::Int,
    observed::AbstractVector{Int}, state_in, P_in, n::Int)
    o = collect(observed)
    push!(tape.updates, CTSEMUpdateRecord(
        o, collect(vec(state_in)), Matrix(P_in),
        Matrix(pars.LAMBDA[o, 1:n]), collect(pars.MANIFESTMEANS[o]),
        Matrix(pars.Jy[o, 1:n]), Matrix(ws.bufferΘ.out[o, o]),
        [float(data[i, obs_col]) for i in o]))
    _tape_push!(tape, :update, length(tape.updates))
    return nothing
end

"""Snapshot the substep inputs, which `_ekf_predict_step!` overwrites."""
function _begin_predict!(tape::CTSEMAdjointTape{T}, ws, n::Int) where {T}
    return (collect(vec(ws.state))::Vector{T}, Matrix(ws.P_update.data[1:n, 1:n])::Matrix{T})
end

function _record_predict!(tape::CTSEMAdjointTape, ws, pars, snapshot, Δt, n::Int)
    state_in, P_in = snapshot
    dyn = ws.diffusion_state_indices
    k = length(dyn)
    push!(tape.predicts, CTSEMPredictRecord(
        state_in, P_in,
        Matrix(ws.discrete_ca.eJAx[1:n, 1:n]),
        Matrix(pars.JAx[1:n, 1:n]), Matrix(pars.DRIFT[1:n, 1:n]),
        Matrix(pars.DIFFUSION[1:n, 1:n]),
        Matrix(ws.diffusion_buffer.out[1:k, 1:k]),
        collect(ws.diffusion_buffer.r[1:k]),
        collect(ws.discrete_ca.dINT[dyn]),
        float(Δt)))
    _tape_push!(tape, :predict, length(tape.predicts))
    return nothing
end

################################################################################
# Reverse pass
################################################################################

"""
    _flush_frechet!(aws)

Push any deferred matrix-exponential Fréchet derivative into the JAx cotangent.

`A = exp(JAx * dt)` is the single most expensive step in the reverse pass: its
adjoint is a matrix exponential of *twice* the state dimension. Profiling the
equivalent C++ port put it at ~66% of the whole reverse pass on a 20-latent
model -- ~540 µs per prediction substep against ~30 µs for everything else in
that substep combined.

But the Fréchet derivative is **linear in its direction argument**, and a
balanced panel design hands it the same `JAx * dt` at every substep of every
subject. So instead of `Σ_steps dt * L(A', Ā_step)`, accumulate the directions
and evaluate `dt * L(A', Σ_steps Ā_step)` once. That is exact, not an
approximation, and it turns an O(rows) count of block exponentials into
O(distinct `(JAx, dt)` pairs) -- one, for a linear model on an equally spaced
panel.

The guard, in `_reverse_predict!`, is an exact comparison against the actual
`JAx * dt` that produced the pending directions, so nothing here assumes
linearity: a state-dependent model changes `JAx` per row, misses every time and
pays only the comparison, while irregular observation times batch within each
distinct interval.

Correctness depends on nothing consuming the JAx cotangent while a flush is
outstanding. Two things can: a state-dependent transform group that *writes* a
JAx cell (whose reverse zeroes that cell's cotangent), and the parameter layer.
Both flush first.
"""
function _flush_frechet!(aws)
    aws.frechet_pending || return nothing
    aws.frechet_pending = false
    contribution = aws.frechet_dt .*
        _ctsem_exp_frechet_adjoint(aws.frechet_A, aws.frechet_accum)
    if aws.defer_frechet
        aws.jax_bar_deferred .+= contribution
        return nothing
    end
    θ̄ca = ComponentVector(aws.theta_bar, aws.sp.parameter_axis)
    n = aws.n
    @inbounds for j in 1:n, i in 1:n
        θ̄ca.JAx[i, j] += contribution[i, j]
    end
    return nothing
end

"""
    _reverse_predict!(x̄, P̄, θ̄ca, record, dyn, n)

Undo one prediction substep.

Forward, on the dynamic index subset `D = dyn` of size `k`:

    A            = exp(JAx * dt)
    affine[i]    = CINT[Dᵢ] + Σⱼ (DRIFT[Dᵢ,j] - JAx[Dᵢ,j]) x[j]
    s[i]         = -affine[i] + Σ_q A[Dᵢ,D_q] affine[q]
    dINT[D]      = JAx[D,D] \\ s
    X            = lyap(JAx[D,D], Qc[D,D])
    dDIFF[D,D]   = X - A[D,D] X A[D,D]'
    x⁺           = A x + dINT
    P⁺           = A (P + εI) A' + dDIFF

Note `A` is the matrix exponential itself: this port sets `dDRIFT = eJAx`, so
the discrete drift and the exponential are the same object and share one
cotangent. `D` is `diffusion_state_indices` -- the states with their own
diffusion. For an ordinary model that is all of them, but `intoverpop`-style
augmentation appends static "carrier" states with structurally zero diffusion,
and the Lyapunov solve, the discrete-intercept solve and the `dDIFF` term all
live on that sub-block while `A` still spans the whole augmented state.
"""
function _reverse_predict!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMPredictRecord{T}, dyn::AbstractVector{Int}, n::Int,
    lyap_buffer, aws) where {T}
    A = record.A
    JAx = record.JAx
    x = record.state_in
    k = length(dyn)
    Ad = A[dyn, dyn]
    JAxd = JAx[dyn, dyn]

    Ā = zeros(T, n, n)
    JAx_bar = zeros(T, n, n)

    # --- mean: x⁺ = A x + dINT
    mul!(Ā, x̄, transpose(x), one(T), one(T))
    dINT_bar = copy(x̄)
    x̄_new = transpose(A) * x̄

    # --- covariance: P⁺ = A (P + εI) A' + dDIFF
    Ps = _symmetrized(P̄)
    Ptilde = copy(record.P_in)
    _ridge_diagonal!(Ptilde, n, _CTSEM_RIDGE)
    # d(A P̃ A')/dA contracted with P̄ is P̄ A P̃' + P̄' A P̃; both P̄ (symmetrised
    # just above) and P̃ are symmetric, so that collapses to twice one term.
    Ā .+= 2 .* (Ps * A * Ptilde)
    P̄_new = transpose(A) * Ps * A
    dDIFF_bar = Ps

    # --- dDIFF[D,D] = X - Ad X Ad'
    X = record.Xlyap
    Qb = _symmetrized(dDIFF_bar[dyn, dyn])
    X̄ = Qb .- transpose(Ad) * Qb * Ad
    Ād = -(Qb * Ad * transpose(X) .+ transpose(Qb) * Ad * X)

    # --- X = lyap(JAx[D,D], Qc[D,D])
    JAxd_bar, Qcd_bar = _ctsem_lyap_pullback(JAxd, X, X̄, lyap_buffer)

    # --- dINT[D] = JAx[D,D] \ s   (a linear solve: s̄ = JAxd⁻ᵀ dINT_bar, M̄ = -s̄ dINT')
    s̄ = transpose(JAxd) \ dINT_bar[dyn]
    JAxd_bar .-= s̄ * transpose(record.dINT_dynamic)

    # --- s = -affine + Ad affine
    affine_bar = -s̄ .+ transpose(Ad) * s̄
    Ād .+= s̄ * transpose(record.affine)

    # --- affine[i] = CINT[Dᵢ] + Σⱼ (DRIFT[Dᵢ,j] - JAx[Dᵢ,j]) x[j]
    @inbounds for i in 1:k
        θ̄ca.CINT[dyn[i]] += affine_bar[i]
    end
    @inbounds for j in 1:n, i in 1:k
        contribution = affine_bar[i] * x[j]
        θ̄ca.DRIFT[dyn[i], j] += contribution
        JAx_bar[dyn[i], j] -= contribution
        x̄_new[j] += (record.DRIFT[dyn[i], j] - JAx[dyn[i], j]) * affine_bar[i]
    end

    # Scatter the dynamic-block contributions back into the full matrices.
    @inbounds for j in 1:k, i in 1:k
        Ā[dyn[i], dyn[j]] += Ād[i, j]
        JAx_bar[dyn[i], dyn[j]] += JAxd_bar[i, j]
    end

    @inbounds for j in 1:n, i in 1:n
        θ̄ca.JAx[i, j] += JAx_bar[i, j]
    end

    # --- A = exp(JAx * dt). The direction is *queued* rather than pushed
    # through immediately, so substeps sharing a `JAx * dt` cost one block
    # exponential between them instead of one each -- see `_flush_frechet!`.
    # The exponential itself still goes through `Base.exp` rather than the
    # package's `my_exp!`, deliberately: buffering it and routing through
    # `my_exp!` was tried and made the adjoint *slower*, up to 2x on small
    # models, because `my_exp!` always runs the degree-13 Padé approximant
    # while `Base.exp` picks a lower degree adaptively from the matrix norm,
    # and `JAx * dt` is typically small-norm.
    scaled = record.dt .* JAx
    # `==` on two `Matrix{T}` is an elementwise comparison and allocates
    # nothing; a `Val(n)`-dispatched helper would, because `n` is a runtime
    # value here and constructing the `Val` costs a dynamic dispatch per row.
    if aws.frechet_pending && aws.frechet_dt == record.dt && aws.frechet_A == scaled
        aws.frechet_accum .+= Ā
    else
        _flush_frechet!(aws)
        copyto!(aws.frechet_A, scaled)
        aws.frechet_dt = record.dt
        copyto!(aws.frechet_accum, Ā)
        aws.frechet_pending = true
    end

    # --- Qc = sdcovsqrt2cov(DIFFUSION); only the dynamic block was consumed.
    Qc_bar = zeros(T, n, n)
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Qcd_bar[i, j]
    end
    diffusion_bar = zeros(T, n, n)
    _sdcovsqrt2cov_pullback!(diffusion_bar, record.DIFFUSION, Qc_bar, n)
    @inbounds for j in 1:n, i in 1:n
        θ̄ca.DIFFUSION[i, j] += diffusion_bar[i, j]
    end

    copyto!(x̄, x̄_new)
    copyto!(P̄, P̄_new)
    return nothing
end

"""Undo one TD-predictor impulse: `x⁺ = x + TDPREDEFFECT td`, `P⁺ = Jtd P Jtd'`."""
function _reverse_td!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMTDRecord{T}, n::Int) where {T}
    td = record.tdpreds
    @inbounds for j in eachindex(td), i in 1:n
        θ̄ca.TDPREDEFFECT[i, j] += x̄[i] * td[j]
    end
    Ps = _symmetrized(P̄)
    Jtd = record.Jtd
    Jtd_bar = Ps * Jtd * transpose(record.P_in) .+ transpose(Ps) * Jtd * record.P_in
    @inbounds for j in 1:n, i in 1:n
        θ̄ca.Jtd[i, j] += Jtd_bar[i, j]
    end
    copyto!(P̄, transpose(Jtd) * Ps * Jtd)
    # The mean update is a pure translation, so x̄ passes through unchanged.
    return nothing
end

"""
    _reverse_update!(x̄, P̄, Θ̄, θ̄ca, record, n)

Undo one measurement update, including its log-likelihood contribution.

Forward (on the observed subset, with `ε` the Stan-matching ridge):

    Pr  = P + εI
    PHt = Pr H'
    S   = sym(H PHt + R) + εI
    ỹ   = y - (Λ x + μ)
    α   = S⁻¹ ỹ
    x⁺  = x + PHt α
    G   = PHt S⁻¹;   M = I - G H
    P⁺  = M P M' + G R G'                (Joseph form, on the *unridged* P)
    ll  = -½ (m log 2π + logdet S + ỹ'α)

The seed for `ll` is 1: the reverse pass differentiates the summed
log-likelihood, so every row contributes with unit weight.
"""
function _reverse_update!(x̄::Vector{T}, P̄::Matrix{T}, Θ̄::Matrix{T}, θ̄ca,
    record::CTSEMUpdateRecord{T}, n::Int) where {T}
    o = record.observed
    m = length(o)
    H = record.H
    Λ = record.Lambda
    R = record.R
    x = record.state_in

    # Recompute the forward intermediates, in the same order and with the same
    # ridges as `_ekf_masked_update_step!`, so the derivative is taken at the
    # values the primal actually used.
    Pr = copy(record.P_in)
    _ridge_diagonal!(Pr, n, _CTSEM_RIDGE)
    PHt = Pr * transpose(H)
    S = _symmetrized(H * PHt .+ R)
    _ridge_diagonal!(S, m, _CTSEM_RIDGE)
    Sinv = inv(S)                     # m is the manifest count: small and dense
    ỹ = record.y .- (Λ * x .+ record.manifestmeans)
    α = Sinv * ỹ
    G = PHt * Sinv
    M = Matrix{T}(I, n, n) .- G * H

    # --- log-likelihood contribution (unit seed)
    S̄ = -0.5 .* (Sinv .- α * transpose(α))
    ỹ̄ = -α

    # --- x⁺ = x + PHt α
    x̄_new = copy(x̄)
    PHt_bar = x̄ * transpose(α)
    ᾱ = transpose(PHt) * x̄

    # --- α = S⁻¹ ỹ
    β = Sinv * ᾱ
    ỹ̄ = ỹ̄ .+ β
    S̄ .-= β * transpose(α)

    # --- P⁺ = M P M' + G R G'
    Ps = _symmetrized(P̄)
    P_in = record.P_in
    M̄ = Ps * M * transpose(P_in) .+ transpose(Ps) * M * P_in
    P̄_new = transpose(M) * Ps * M
    Ḡ = Ps * G * transpose(R) .+ transpose(Ps) * G * R
    R̄ = transpose(G) * Ps * G

    # --- M = I - G H
    Ḡ .-= M̄ * transpose(H)
    H̄ = -transpose(G) * M̄

    # --- G = PHt S⁻¹
    PHt_bar .+= Ḡ * Sinv
    S̄ .-= transpose(G) * Ḡ * Sinv

    # --- S = sym(H PHt + R) + εI
    S̄0 = _symmetrized(S̄)
    H̄ .+= S̄0 * transpose(PHt)
    PHt_bar .+= transpose(H) * S̄0
    R̄ .+= S̄0

    # --- PHt = Pr H'
    P̄_new .+= PHt_bar * H
    H̄ .+= transpose(PHt_bar) * Pr

    # --- ỹ = y - (Λ x + μ)
    Λ̄ = -ỹ̄ * transpose(x)
    x̄_new .-= transpose(Λ) * ỹ̄

    # --- scatter back into the full-size parameter and Θ cotangents
    @inbounds for i in 1:m
        θ̄ca.MANIFESTMEANS[o[i]] -= ỹ̄[i]
        for j in 1:n
            θ̄ca.LAMBDA[o[i], j] += Λ̄[i, j]
            θ̄ca.Jy[o[i], j] += H̄[i, j]
        end
        for j in 1:m
            Θ̄[o[i], o[j]] += R̄[i, j]
        end
    end

    copyto!(x̄, x̄_new)
    copyto!(P̄, P̄_new)
    return nothing
end

@inline function _symmetrized(A::AbstractMatrix)
    return (A .+ transpose(A)) ./ 2
end

"""
    _ctsem_reverse_tape!(tape, sp, aws, n, m)

Walk one subject's tape backwards, *accumulating* into `aws.theta_bar`, the
cotangent of the materialized parameter matrices.

Deliberately stops there rather than going on to the free parameters: whether
`theta_bar` can be shared across subjects (and the parameter layer therefore
unwound once instead of once per subject) is a property of the model, not of
one tape. See `_ctsem_parameter_layer_shareable`.

`theta_bar` is *not* zeroed here. The caller owns that, because sharing is the
whole point.
"""
function _ctsem_reverse_tape!(tape::CTSEMAdjointTape{T},
    sp::EKFParameters, aws, n::Int, m::Int) where {T}
    x̄ = aws.x_bar
    P̄ = aws.P_bar
    Θ̄ = aws.manifest_cov_bar
    θ̄ = aws.theta_bar
    fill!(x̄, zero(T)); fill!(P̄, zero(T)); fill!(Θ̄, zero(T))
    # Constructed here rather than read from the workspace on purpose: a
    # `ComponentVector` view is a free wrapper around `θ̄`, but a workspace
    # field holding one has to be `Any`-typed (its axis type is model
    # specific), and reading it would make every `θ̄ca.DRIFT[i,j] +=` in the
    # loop below a dynamic dispatch. Measured at 45% of the whole reverse
    # pass when this was done the other way.
    θ̄ca = ComponentVector(θ̄, sp.parameter_axis)
    dyn = aws.diffusion_state_indices

    for entry_index in length(tape.program):-1:1
        kind, index = tape.program[entry_index]
        if kind === :update
            _reverse_update!(x̄, P̄, Θ̄, θ̄ca, tape.updates[index], n)
        elseif kind === :td
            _reverse_td!(x̄, P̄, θ̄ca, tape.tds[index], n)
        elseif kind === :predict
            _reverse_predict!(x̄, P̄, θ̄ca, tape.predicts[index], dyn, n, aws.lyap_buffer, aws)
        elseif kind === :group
            # Only a group that writes a JAx cell can consume a queued Fréchet
            # contribution (its reverse zeroes that cell's cotangent).
            # Flushing before *every* group instead would break the batch on
            # ctsem's default model, whose individually varying MANIFESTMEANS
            # puts a group between every pair of prediction substeps in an
            # otherwise entirely linear model.
            aws.groups_write_jax && _flush_frechet!(aws)
            _reverse_group!(θ̄, x̄, tape.groups[index], sp, aws)
        elseif kind === :theta
            manifestvar_bar = zeros(T, m, m)
            _sdcovsqrt2cov_pullback!(manifestvar_bar, tape.thetas[index].MANIFESTVAR, Θ̄, m)
            @inbounds for j in 1:m, i in 1:m
                θ̄ca.MANIFESTVAR[i, j] += manifestvar_bar[i, j]
            end
            fill!(Θ̄, zero(T))
        elseif kind === :init
            @inbounds for i in 1:n
                θ̄ca.T0MEANS[i] += x̄[i]
            end
            fill!(x̄, zero(T))
            t0var_bar = zeros(T, n, n)
            _sdcovsqrt2cov_pullback!(t0var_bar, tape.inits[index].T0VAR, _symmetrized(P̄), n)
            @inbounds for j in 1:n, i in 1:n
                θ̄ca.T0VAR[i, j] += t0var_bar[i, j]
            end
            fill!(P̄, zero(T))
        else
            throw(ArgumentError("adjoint: unknown tape entry $(kind)"))
        end
    end

    return θ̄
end

"""
    _ctsem_parameter_layer!(values_bar, θ̄, subject_values, sp, aws, tipreds)

Unwind the parameter layer: whatever cotangent is left on `all_params` was put
there by the regular (non-state-dependent) transforms, so push it through them
and then through the TI-predictor effects into `values_bar`.
"""
function _ctsem_parameter_layer!(values_bar::AbstractVector{T}, θ̄::AbstractVector{T},
    subject_values::AbstractVector, sp::EKFParameters, aws,
    tipreds::AbstractVector) where {T}
    aws.defer_frechet || _flush_frechet!(aws)
    subject_values_bar = aws.subject_values_bar
    resize!(subject_values_bar, length(subject_values))
    fill!(subject_values_bar, zero(T))
    _ctsem_regular_pullback!(subject_values_bar, θ̄, subject_values, sp,
        aws.regular_supports, aws.regular_dual_scratch)
    _ctsem_ti_pullback!(values_bar, subject_values_bar, sp, tipreds)
    return values_bar
end

"""Reverse one recorded state-dependent transform group."""
function _reverse_group!(θ̄::Vector{T}, x̄::Vector{T}, record::CTSEMGroupRecord{T},
    sp::EKFParameters, aws) where {T}
    transforms, indices, supports = if record.group == 1
        sp.predict_transforms, aws.predict_indices, aws.predict_supports
    elseif record.group == 2
        sp.td_transforms, aws.td_indices, aws.td_supports
    else
        sp.update_transforms, aws.update_indices, aws.update_supports
    end
    isempty(indices) && return nothing
    # `group_scratch` is a persistent full-length buffer whose entries outside
    # this group's relevant set are never written and stay `NaN` for the life
    # of the workspace (see `CTSEMAdjointWorkspace`). Scatter the compact
    # record into the relevant positions; anything a transform reads that we
    # failed to record therefore reads `NaN` and poisons the gradient loudly
    # rather than silently using a stale number.
    relevant = aws.tape.group_relevant[record.group]
    working = aws.group_scratch
    @inbounds for j in eachindex(relevant)
        working[relevant[j]] = record.params_before[j]
    end
    pars = ComponentVector(working, sp.parameter_axis)
    ctx = CTSEMRowContext(record.state, pars, record.tdpreds, aws.tipreds,
        record.time, record.dt, 1, record.row)
    _ctsem_complex_group_pullback!(θ̄, x̄, transforms, indices, supports, ctx,
        relevant, aws.dual_context)
    return nothing
end
