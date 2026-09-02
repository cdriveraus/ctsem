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
mutable struct CTSEMPredictRecord{T}
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
mutable struct CTSEMTDRecord{T}
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
mutable struct CTSEMUpdateRecord{T}
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
mutable struct CTSEMGroupRecord{T}
    group::Int               # 1 = predict, 2 = td, 3 = update
    params_before::Vector{T} # values at `group_relevant[group]`, pre-group
    state::Vector{T}
    tdpreds::Vector{T}
    time::T
    dt::T
    row::Int
end

"""Construction of Θ from MANIFESTVAR."""
mutable struct CTSEMThetaRecord{T}
    MANIFESTVAR::Matrix{T}
end

"""The `t = 1` prior: state from T0MEANS, covariance from T0VAR."""
mutable struct CTSEMInitRecord{T}
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
    binaries::Vector{CTSEMBinaryRecord{T}}
    # How many of each vector the *current* pass has written. The vectors are
    # never emptied, so a subject after the first writes into records that
    # already exist and already have arrays of the right shape -- see
    # `_tape_reset!`. `program` is the exception: `empty!` keeps a Vector's
    # capacity, so pushing into it again allocates nothing.
    npredicts::Int
    ntds::Int
    nupdates::Int
    ngroups::Int
    nthetas::Int
    ninits::Int
    nbinaries::Int
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
        CTSEMInitRecord{T}[], CTSEMBinaryRecord{T}[], 0, 0, 0, 0, 0, 0, 0,
        T[], group_relevant)

"""
    _tape_reset!(tape)

Begin a new pass, keeping the records the last one built.

Emptying the record vectors is correct and was what this did; it is also why a
24-row subject of a one-latent model allocated 1,263 arrays in its forward
sweep alone, every one of them a 1x1 `Matrix{Float64}` that the previous
subject had already allocated at exactly that shape. The tape is per chunk and
a chunk's subjects share a model, so after the first subject every shape is
already there.

Resetting counts instead means the records are *reused*: `_record_predict!` and
friends copy into `tape.predicts[i]` when `i` is within reach and only build
when the tape has to grow. Nothing outside the recording hooks reads past the
count, and the reverse pass reads records by the index `program` recorded, so a
stale record beyond the count is unreachable rather than merely unread.
"""
function _tape_reset!(tape::CTSEMAdjointTape)
    empty!(tape.program)
    tape.npredicts = 0; tape.ntds = 0; tape.nupdates = 0
    tape.ngroups = 0; tape.nthetas = 0; tape.ninits = 0
    tape.nbinaries = 0
    return tape
end

@inline _tape_push!(tape::CTSEMAdjointTape, kind::Symbol, index::Int) =
    push!(tape.program, (kind, index))

"""
    _tape_fill!(destination, source)

Overwrite a pooled record field, reallocating only when the shape changed.

A vector is `resize!`d, which keeps its capacity across resets and so is free
after the first pass. A matrix cannot be reshaped in place, so a genuine shape
change -- which for these records means a row with a different set of observed
manifests -- builds a new one and the caller stores it back. That is why the
records are mutable.
"""
@inline function _tape_fill!(destination::Vector{T}, source) where {T}
    resize!(destination, length(source))
    copyto!(destination, source)
    return destination
end

@inline function _tape_fill!(destination::Matrix{T}, source::AbstractMatrix) where {T}
    size(destination) == size(source) || return Matrix{T}(source)
    copyto!(destination, source)
    return destination
end

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

function _record_init!(tape::CTSEMAdjointTape{T}, pars, n::Int) where {T}
    index = (tape.ninits += 1)
    source = view(pars.T0VAR, 1:n, 1:n)
    if index <= length(tape.inits)
        record = tape.inits[index]
        record.T0VAR = _tape_fill!(record.T0VAR, source)
    else
        push!(tape.inits, CTSEMInitRecord{T}(Matrix{T}(source)))
    end
    _tape_push!(tape, :init, index)
    return nothing
end

function _record_theta!(tape::CTSEMAdjointTape{T}, pars, m::Int) where {T}
    index = (tape.nthetas += 1)
    source = view(pars.MANIFESTVAR, 1:m, 1:m)
    if index <= length(tape.thetas)
        record = tape.thetas[index]
        record.MANIFESTVAR = _tape_fill!(record.MANIFESTVAR, source)
    else
        push!(tape.thetas, CTSEMThetaRecord{T}(Matrix{T}(source)))
    end
    _tape_push!(tape, :theta, index)
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
    index = (tape.ngroups += 1)
    # Explicit `{T}` for the same reason as in `_record_update!`: the TD
    # predictors and the time/interval are data, so they arrive as Float64
    # whatever the tape's element type is, and the record requires every field
    # to share it.
    if index <= length(tape.groups)
        record = tape.groups[index]
        record.group = group
        # `params_before` is as long as the group's relevant set, and the three
        # groups have different sets -- so this is the one pooled field whose
        # length genuinely varies between neighbouring records. `resize!` keeps
        # the capacity, so it still costs nothing after the first pass.
        resize!(record.params_before, length(relevant))
        @inbounds for j in eachindex(relevant)
            record.params_before[j] = all_params[relevant[j]]
        end
        _tape_fill!(record.state, vec(ctx.state))
        _tape_fill!(record.tdpreds, ctx.tdpreds)
        record.time = T(ctx.time)
        record.dt = T(ctx.dt)
        record.row = ctx.row
    else
        compact = Vector{T}(undef, length(relevant))
        @inbounds for j in eachindex(relevant)
            compact[j] = all_params[relevant[j]]
        end
        push!(tape.groups, CTSEMGroupRecord{T}(group, compact, collect(vec(ctx.state)),
            T[x for x in ctx.tdpreds], T(ctx.time), T(ctx.dt), ctx.row))
    end
    _tape_push!(tape, :group, index)
    return nothing
end

function _record_td!(tape::CTSEMAdjointTape{T}, ws, pars, tdpreds, n::Int) where {T}
    # Mirrors `_apply_td_impulse!`'s own early return: with no TD predictors
    # the impulse is a no-op and contributes nothing to reverse.
    isempty(tdpreds) && return nothing
    index = (tape.ntds += 1)
    P_in = view(ws.P_predict.data, 1:n, 1:n)
    Jtd = view(pars.Jtd, 1:n, 1:n)
    if index <= length(tape.tds)
        record = tape.tds[index]
        record.P_in = _tape_fill!(record.P_in, P_in)
        record.Jtd = _tape_fill!(record.Jtd, Jtd)
        _tape_fill!(record.tdpreds, tdpreds)
    else
        push!(tape.tds, CTSEMTDRecord{T}(Matrix{T}(P_in), Matrix{T}(Jtd),
            T[x for x in tdpreds]))
    end
    _tape_push!(tape, :td, index)
    return nothing
end

function _record_update!(tape::CTSEMAdjointTape{T}, ws, pars, data, obs_col::Int,
    observed::AbstractVector{Int}, state_in, P_in, n::Int) where {T}
    index = (tape.nupdates += 1)
    # `T[...]` rather than `[float(...)]`: the observed data is a constant with
    # respect to the parameters, so it is Float64 whatever the tape's element
    # type is. Every other field is `T`, and `CTSEMUpdateRecord{T}` requires
    # them all to agree -- which they did as long as `T` was `Float64` too.
    # Differentiating this gradient (see `ctsem_hessian`) makes `T` a dual, and
    # the record then has to carry the data as a dual with zero partials.
    if index > length(tape.updates)
        o = collect(observed)
        push!(tape.updates, CTSEMUpdateRecord{T}(
            o, collect(vec(state_in)), Matrix(P_in),
            Matrix(pars.LAMBDA[o, 1:n]), collect(pars.MANIFESTMEANS[o]),
            Matrix(pars.Jy[o, 1:n]), Matrix(ws.bufferΘ.out[o, o]),
            T[data[i, obs_col] for i in o]))
        _tape_push!(tape, :update, index)
        return nothing
    end
    record = tape.updates[index]
    o = _tape_fill!(record.observed, observed)
    # Views rather than slices throughout: `pars.LAMBDA[o, 1:n]` materialises a
    # matrix only to copy it and throw it away, and this runs once per row.
    _tape_fill!(record.state_in, state_in)
    record.P_in = _tape_fill!(record.P_in, P_in)
    record.Lambda = _tape_fill!(record.Lambda, view(pars.LAMBDA, o, 1:n))
    _tape_fill!(record.manifestmeans, view(pars.MANIFESTMEANS, o))
    record.H = _tape_fill!(record.H, view(pars.Jy, o, 1:n))
    record.R = _tape_fill!(record.R, view(ws.bufferΘ.out, o, o))
    resize!(record.y, length(o))
    @inbounds for j in eachindex(o)
        record.y[j] = T(data[o[j], obs_col])
    end
    _tape_push!(tape, :update, index)
    return nothing
end

"""An all-zero predict record shaped for `n` states and `k` dynamic states."""
_empty_predict_record(::Type{T}, n::Int, k::Int) where {T} =
    CTSEMPredictRecord{T}(zeros(T, n), zeros(T, n, n), zeros(T, n, n),
        zeros(T, n, n), zeros(T, n, n), zeros(T, n, n), zeros(T, k, k),
        zeros(T, k), zeros(T, k), zero(T))

"""
Snapshot the substep inputs, which `_ekf_predict_step!` overwrites.

Returns the pooled record the matching `_record_predict!` will finish, rather
than a fresh tuple -- the two are called in lockstep around one
`_ekf_predict_step!`, so the slot is known here and the snapshot can go straight
into it. Growing the tape here rather than there keeps the return type a plain
`CTSEMPredictRecord`, so the filter's local stays concrete.
"""
function _begin_predict!(tape::CTSEMAdjointTape{T}, ws, n::Int) where {T}
    index = tape.npredicts + 1
    if index > length(tape.predicts)
        push!(tape.predicts,
            _empty_predict_record(T, n, length(ws.diffusion_state_indices)))
    end
    record = tape.predicts[index]
    _tape_fill!(record.state_in, ws.state)
    record.P_in = _tape_fill!(record.P_in, view(ws.P_update.data, 1:n, 1:n))
    return record
end

function _record_predict!(tape::CTSEMAdjointTape{T}, ws, pars,
    record::CTSEMPredictRecord{T}, Δt, n::Int) where {T}
    dyn = ws.diffusion_state_indices
    k = length(dyn)
    index = (tape.npredicts += 1)
    record.A = _tape_fill!(record.A, view(ws.discrete_ca.eJAx, 1:n, 1:n))
    record.JAx = _tape_fill!(record.JAx, view(pars.JAx, 1:n, 1:n))
    record.DRIFT = _tape_fill!(record.DRIFT, view(pars.DRIFT, 1:n, 1:n))
    record.DIFFUSION = _tape_fill!(record.DIFFUSION, view(pars.DIFFUSION, 1:n, 1:n))
    record.Xlyap = _tape_fill!(record.Xlyap, view(ws.diffusion_buffer.out, 1:k, 1:k))
    _tape_fill!(record.affine, view(ws.diffusion_buffer.r, 1:k))
    _tape_fill!(record.dINT_dynamic, view(ws.discrete_ca.dINT, dyn))
    record.dt = T(Δt)
    _tape_push!(tape, :predict, index)
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
    _reverse_predict!(x̄, P̄, θ̄ca, record, dyn, n, lyap_buffer, aws)

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

    # The discrete-time reverse pass is the continuous one with its three hard
    # pieces removed rather than a second implementation of it: A is JAx (no
    # Frechet derivative of an exponential), dDIFFUSION is the diffusion
    # covariance itself (no Lyapunov pullback), and dINT is the affine offset
    # (no linear solve). Written out separately rather than branched inside the
    # shared recursion, which would put a test in every step of it.
    if !aws.sp.continuous_time
        return _reverse_predict_discrete!(x̄, P̄, θ̄ca, record, dyn, n)
    end

    # Every temporary is a view of a scratch buffer; see `CTSEMReverseScratch`.
    # The derivation below is unchanged -- what changed is only where the memory
    # comes from. Buffers are shared with `_reverse_update!`, which is safe
    # because the two are called at different points of the tape walk and
    # neither is ever live while the other runs.
    sc = aws.reverse_scratch
    Ā          = _rs(sc.Abar, n, n)
    JAx_bar    = _rs(sc.nn3, n, n)
    P̄_new      = _rs(sc.Pbar_new, n, n)
    Ps         = _rs(sc.Ps, n, n)
    Ptilde     = _rs(sc.Pr, n, n)
    nn1        = _rs(sc.nn1, n, n)
    nn2        = _rs(sc.nn2, n, n)
    scaled     = _rs(sc.nn4, n, n)
    Qc_bar     = _rs(sc.nn5, n, n)
    diffusion_bar = _rs(sc.nn6, n, n)
    Ad         = _rs(sc.kk1, k, k)
    JAxd       = _rs(sc.kk2, k, k)
    Qb         = _rs(sc.kk3, k, k)
    X̄          = _rs(sc.kk4, k, k)
    Ād         = _rs(sc.kk5, k, k)
    kk6        = _rs(sc.kk6, k, k)
    x̄_new      = _rs(sc.xbar_new, n)
    dINT_bar   = _rs(sc.nv1, n)
    s̄          = _rs(sc.kv1, k)
    affine_bar = _rs(sc.kv2, k)

    @inbounds for j in 1:k, i in 1:k
        Ad[i, j] = A[dyn[i], dyn[j]]
        JAxd[i, j] = JAx[dyn[i], dyn[j]]
    end
    fill!(Ā, zero(T))
    fill!(JAx_bar, zero(T))

    # --- mean: x⁺ = A x + dINT
    _ctsem_outer!(Ā, x̄, x, one(T), one(T))
    copyto!(dINT_bar, x̄)
    _ctsem_mulTvec!(x̄_new, A, x̄)

    # --- covariance: P⁺ = A (P + εI) A' + dDIFF
    _symmetrize_into!(Ps, P̄)
    copyto!(Ptilde, record.P_in)
    _ridge_diagonal!(Ptilde, n, _CTSEM_RIDGE)
    # d(A P̃ A')/dA contracted with P̄ is P̄ A P̃' + P̄' A P̃; both P̄ (symmetrised
    # just above) and P̃ are symmetric, so that collapses to twice one term.
    _ctsem_mul!(nn1, Ps, A)
    _ctsem_mul!(Ā, nn1, Ptilde, T(2), one(T))
    _ctsem_mulTN!(nn2, A, Ps)
    _ctsem_mul!(P̄_new, nn2, A)
    # `dDIFF_bar` is `Ps`, so `Ps` has to survive until `Qb` is taken from it.

    # --- dDIFF[D,D] = X - Ad X Ad'
    X = record.Xlyap
    @inbounds for j in 1:k, i in 1:k
        Qb[i, j] = (Ps[dyn[i], dyn[j]] + Ps[dyn[j], dyn[i]]) / 2
    end
    _ctsem_mulTN!(kk6, Ad, Qb)
    _ctsem_mul!(X̄, kk6, Ad)
    X̄ .= Qb .- X̄                                   # X̄ = Qb - Ad' Qb Ad
    _ctsem_mul!(kk6, Qb, Ad)
    _ctsem_mulNT!(Ād, kk6, X)
    _ctsem_mulTN!(kk6, Qb, Ad)
    _ctsem_mul!(Ād, kk6, X, one(T), one(T))
    Ād .= .-Ād                                     # Ād = -(Qb Ad X' + Qb' Ad X)

    # --- X = lyap(JAx[D,D], Qc[D,D])
    JAxd_bar, Qcd_bar = _ctsem_lyap_pullback(JAxd, X, X̄, lyap_buffer)

    # --- dINT[D] = JAx[D,D] \ s   (a linear solve: s̄ = JAxd⁻ᵀ dINT_bar, M̄ = -s̄ dINT')
    @inbounds for i in 1:k; s̄[i] = dINT_bar[dyn[i]]; end
    # `transpose(JAxd) \ s̄` would be a LAPACK `getrf!`/`getrs!` pair, once per
    # prediction substep, on a matrix of the dynamic-state size. At that size
    # the call is mostly OpenBLAS's process-global buffer lock -- see
    # `small_linalg.jl` -- so it goes through the engine's own LU instead.
    @inbounds for j in 1:k, i in 1:k; kk6[i, j] = JAxd[j, i]; end
    _solve_square_system_generic!(kk6, s̄, sc.piv, Val(k))
    _ctsem_outer!(JAxd_bar, s̄, record.dINT_dynamic, -one(T), one(T))

    # --- s = -affine + Ad affine
    _ctsem_mulTvec!(affine_bar, Ad, s̄)
    affine_bar .-= s̄
    _ctsem_outer!(Ād, s̄, record.affine, one(T), one(T))

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
    scaled .= record.dt .* JAx
    # `==` on two matrices is an elementwise comparison and allocates nothing;
    # a `Val(n)`-dispatched helper would, because `n` is a runtime value here
    # and constructing the `Val` costs a dynamic dispatch per row.
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
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Qcd_bar[i, j]
    end
    fill!(diffusion_bar, zero(T))
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
    _reverse_update!(x̄, P̄, Θ̄, θ̄ca, record, n, sc)

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
    record::CTSEMUpdateRecord{T}, n::Int, sc::CTSEMReverseScratch{T}) where {T}
    o = record.observed
    m = length(o)
    H = record.H
    Λ = record.Lambda
    R = record.R
    x = record.state_in

    # Every temporary is a view of a scratch buffer; see `CTSEMReverseScratch`.
    # The names and the order are the derivation's, unchanged -- what changed is
    # only where the memory comes from.
    Pr     = _rs(sc.Pr, n, n)
    PHt    = _rs(sc.PHt, n, m)
    S      = _rs(sc.S, m, m)
    Sinv   = _rs(sc.Sinv, m, m)
    G      = _rs(sc.G, n, m)
    M      = _rs(sc.M, n, n)
    S̄      = _rs(sc.Sbar, m, m)
    S̄0     = _rs(sc.Sbar0, m, m)
    R̄      = _rs(sc.Rbar, m, m)
    M̄      = _rs(sc.Mbar, n, n)
    P̄_new  = _rs(sc.Pnew, n, n)
    Ps     = _rs(sc.Ps, n, n)
    Ḡ      = _rs(sc.Gbar, n, m)
    H̄      = _rs(sc.Hbar, m, n)
    Λ̄      = _rs(sc.Lbar, m, n)
    PHt_bar = _rs(sc.PHtbar, n, m)
    ỹ      = _rs(sc.ytilde, m)
    ỹ̄      = _rs(sc.ybar, m)
    α      = _rs(sc.alpha, m)
    ᾱ      = _rs(sc.alphabar, m)
    β      = _rs(sc.beta, m)
    x̄_new  = _rs(sc.xnew, n)
    nn1    = _rs(sc.nn1, n, n)
    nn2    = _rs(sc.nn2, n, n)
    nm1    = _rs(sc.nm1, n, m)
    mm1    = _rs(sc.mm1, m, m)
    mm2    = _rs(sc.mm2, m, m)

    # Recompute the forward intermediates, in the same order and with the same
    # ridges as `_ekf_masked_update_step!`, so the derivative is taken at the
    # values the primal actually used.
    copyto!(Pr, record.P_in)
    _ridge_diagonal!(Pr, n, _CTSEM_RIDGE)
    _ctsem_mulNT!(PHt, Pr, H)                       # PHt = Pr H'
    _ctsem_mul!(mm1, H, PHt)                                 # mm1 = H PHt
    mm1 .+= R
    _symmetrize_into!(S, mm1)                         # S = sym(H PHt + R)
    _ridge_diagonal!(S, m, _CTSEM_RIDGE)
    # `S` is the innovation covariance: symmetric, positive definite, and the
    # size of the manifest vector. `inv` on it is a LAPACK `getrf!`/`getri!`
    # pair, which at this size is mostly OpenBLAS's process-global buffer lock
    # -- it was 6.4% of the reverse pass's samples and it does not thread. The
    # engine's own Cholesky plus `m` triangular solves is the same inverse
    # without the lock. See `small_linalg.jl`.
    copyto!(mm2, S)
    if _ctsem_cholesky!(mm2, m)
        F = CTSEMCholesky(mm2, m, true)
        @inbounds for j in 1:m
            column = view(Sinv, :, j)
            fill!(column, zero(T))
            column[j] = one(T)
            ldiv!(column, F, column)
        end
    else
        # Not `inv(S)`. The reverse recomputes `S` with its own kernels rather
        # than reusing the forward's factorization, so it can fail here on a
        # marginally definite `S` the forward accepted -- and `inv` of a
        # symmetric matrix that just failed to factorize returns something
        # large but finite, so the whole reverse pass would run on to a finite,
        # wrong gradient that `fg!` accepts and feeds to the L-BFGS curvature
        # update. Poison instead, which is how the forward already treats this
        # state: an invalid trial point.
        #
        # `θ̄ca` and not only `x̄`/`P̄`/`Θ̄`, because `_ctsem_regular_pullback!`
        # drops cotangents at non-mutable positions -- a poison confined to the
        # state cotangents could be filtered out of a model whose measurement
        # matrices are all fixed.
        fill!(x̄, T(NaN))
        fill!(P̄, T(NaN))
        fill!(Θ̄, T(NaN))
        fill!(θ̄ca, T(NaN))
        return nothing
    end
    copyto!(ỹ, record.manifestmeans)
    _ctsem_mulvec!(ỹ, Λ, x, one(T), one(T))                     # ỹ = Λ x + μ
    ỹ .= record.y .- ỹ
    _ctsem_mulvec!(α, Sinv, ỹ)
    _ctsem_mul!(G, PHt, Sinv)
    _ctsem_mul!(M, G, H)                                     # M = I - G H
    M .= .-M
    @inbounds for i in 1:n; M[i, i] += one(T); end

    # --- log-likelihood contribution (unit seed)
    _ctsem_outer!(mm1, α, α)
    S̄ .= -0.5 .* (Sinv .- mm1)
    ỹ̄ .= .-α

    # --- x⁺ = x + PHt α
    copyto!(x̄_new, x̄)
    _ctsem_outer!(PHt_bar, x̄, α)
    _ctsem_mulTvec!(ᾱ, PHt, x̄)

    # --- α = S⁻¹ ỹ
    _ctsem_mulvec!(β, Sinv, ᾱ)
    ỹ̄ .+= β
    _ctsem_outer!(mm1, β, α)
    S̄ .-= mm1

    # --- P⁺ = M P M' + G R G'
    _symmetrize_into!(Ps, P̄)
    P_in = record.P_in
    _ctsem_mul!(nn1, Ps, M)                           # nn1 = Ps M
    _ctsem_mulNT!(M̄, nn1, P_in)
    _ctsem_mulTN!(nn2, Ps, M)
    _ctsem_mul!(M̄, nn2, P_in, one(T), one(T))                # M̄ = Ps M P' + Ps' M P
    _ctsem_mulTN!(nn1, M, Ps)
    _ctsem_mul!(P̄_new, nn1, M)                               # P̄_new = M' Ps M
    _ctsem_mul!(nm1, Ps, G)
    _ctsem_mulNT!(Ḡ, nm1, R)
    _ctsem_mulTN!(nm1, Ps, G)
    _ctsem_mul!(Ḡ, nm1, R, one(T), one(T))                   # Ḡ = Ps G R' + Ps' G R
    _ctsem_mul!(nm1, Ps, G)
    _ctsem_mulTN!(R̄, G, nm1)                        # R̄ = G' Ps G

    # --- M = I - G H
    _ctsem_mulNT!(Ḡ, M̄, H, -one(T), one(T))
    _ctsem_mulTN!(H̄, G, M̄)
    H̄ .= .-H̄

    # --- G = PHt S⁻¹
    _ctsem_mul!(PHt_bar, Ḡ, Sinv, one(T), one(T))
    _ctsem_mulTN!(mm1, G, Ḡ)
    _ctsem_mul!(S̄, mm1, Sinv, -one(T), one(T))

    # --- S = sym(H PHt + R) + εI
    _symmetrize_into!(S̄0, S̄)
    _ctsem_mulNT!(H̄, S̄0, PHt, one(T), one(T))
    _ctsem_mulTN!(PHt_bar, H, S̄0, one(T), one(T))
    R̄ .+= S̄0

    # --- PHt = Pr H'
    _ctsem_mul!(P̄_new, PHt_bar, H, one(T), one(T))
    _ctsem_mulTN!(H̄, PHt_bar, Pr, one(T), one(T))

    # --- ỹ = y - (Λ x + μ)
    _ctsem_outer!(Λ̄, ỹ̄, x)
    Λ̄ .= .-Λ̄
    _ctsem_mulTvec!(x̄_new, Λ, ỹ̄, -one(T), one(T))

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
            _reverse_update!(x̄, P̄, Θ̄, θ̄ca, tape.updates[index], n, aws.reverse_scratch)
        elseif kind === :binary
            # After the Gaussian block of the same row, because the forward
            # applied binary observations before it.
            _reverse_binary!(x̄, P̄, θ̄ca, tape.binaries[index], n)
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

"""
    _reverse_predict_discrete!(XBAR, PBAR, THETA, record, dyn, n)

Undo one prediction step of a discrete-time model.

Forward:

    A          = JAx
    dINT[i]    = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j]
    dDIFF[D,D] = Qc[D,D]
    x_next     = A x + dINT
    P_next     = A (P + eps I) A' + dDIFF

which is the continuous form with the exponential, the Lyapunov solve and the
intercept solve all collapsed -- a discrete model's DRIFT, CINT and DIFFUSION
are already the one-step quantities. The affine offset runs over every state
rather than only the diffusing ones, matching the forward pass: with no solve to
keep away from the singular augmented block there is no reason to restrict it.
"""
function _reverse_predict_discrete!(XBAR::Vector{T}, PBAR::Matrix{T}, THETA,
    record::CTSEMPredictRecord{T}, dyn::AbstractVector{Int}, n::Int) where {T}
    A = record.A
    JAx = record.JAx
    x = record.state_in
    k = length(dyn)

    JAx_bar = zeros(T, n, n)

    # --- mean: x_next = A x + dINT
    Abar = XBAR * transpose(x)
    dINT_bar = copy(XBAR)
    xbar_new = transpose(A) * XBAR

    # --- covariance: P_next = A (P + eps I) A' + dDIFF
    Ps = _symmetrized(PBAR)
    Ptilde = copy(record.P_in)
    _ridge_diagonal!(Ptilde, n, _CTSEM_RIDGE)
    Abar .+= 2 .* (Ps * A * Ptilde)
    Pbar_new = transpose(A) * Ps * A

    # --- dDIFF[D,D] = Qc[D,D]: the cotangent passes straight through.
    Qcd_bar = _symmetrized(Ps[dyn, dyn])

    # --- dINT[i] = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j]
    @inbounds for i in 1:n
        THETA.CINT[i] += dINT_bar[i]
    end
    @inbounds for j in 1:n, i in 1:n
        contribution = dINT_bar[i] * x[j]
        THETA.DRIFT[i, j] += contribution
        JAx_bar[i, j] -= contribution
        xbar_new[j] += (record.DRIFT[i, j] - JAx[i, j]) * dINT_bar[i]
    end

    # --- A is JAx itself, so its cotangent simply adds.
    JAx_bar .+= Abar
    @inbounds for j in 1:n, i in 1:n
        THETA.JAx[i, j] += JAx_bar[i, j]
    end

    # --- Qc = sdcovsqrt2cov(DIFFUSION); only the dynamic block was consumed.
    Qc_bar = zeros(T, n, n)
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Qcd_bar[i, j]
    end
    diffusion_bar = zeros(T, n, n)
    _sdcovsqrt2cov_pullback!(diffusion_bar, record.DIFFUSION, Qc_bar, n)
    @inbounds for j in 1:n, i in 1:n
        THETA.DIFFUSION[i, j] += diffusion_bar[i, j]
    end

    copyto!(XBAR, xbar_new)
    copyto!(PBAR, Pbar_new)
    return nothing
end
