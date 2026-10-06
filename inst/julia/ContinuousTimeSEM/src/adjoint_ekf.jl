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

# Distinct `JAx * dt` keys the reverse pass holds Fréchet directions for at
# once. A shared wave schedule with nineteen distinct intervals batched to one
# Fréchet evaluation per substep under a last-value rule -- 3800 per gradient
# on 200 subjects of 20 rows -- and batches to nineteen with a table. Fully
# irregular times fill the table every 32 substeps and flush it, which costs
# what the last-value rule did plus a scan of 32 scalars per substep.
const _CTSEM_FRECHET_TABLE = 32

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
    # Which pieces took the series route (`series_discretization.jl`), and the
    # diffusion block the noise series integrated, which its pullback needs
    # where the closed form needed only `Xlyap`. Empty until a substep uses
    # it, so a model that never takes the series carries no extra storage.
    series_intercept::Bool
    series_noise::Bool
    Qd::Matrix{T}
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
    # The raw population covariance matrix, at the point the forward pass
    # constructed it. `0` by `0` for a model whose population covariance is not
    # a matrix of its own, which is every model that does not use intoverpop.
    POPRAW::Matrix{T}
end

"""
The stationary initial moments, `_ctsem_stationary!`: the drift over the
leading block of genuine dynamics, the mean it gave, the covariance over the
diffusion states, and the raw DIFFUSION, which keys the deferred diffusion
pullback exactly as a prediction's does.
"""
mutable struct CTSEMStationaryRecord{T}
    DRIFT::Matrix{T}
    mean::Vector{T}
    X::Matrix{T}
    DIFFUSION::Matrix{T}
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
    stationaries::Vector{CTSEMStationaryRecord{T}}
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
    nstationaries::Int
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
        CTSEMInitRecord{T}[], CTSEMBinaryRecord{T}[],
        CTSEMStationaryRecord{T}[], 0, 0, 0, 0, 0, 0, 0, 0,
        T[], group_relevant)

"""
Transition growth above which a subject's gradient is taken by forward mode
rather than by the reverse pass. The growth is that of the product of every
transition since the last observed row (`CTSEMGrowth`).

The reverse pass is exact in exact arithmetic and not stable in double
precision once a prediction's covariance is very large in some direction: the
covariance cotangent there comes from an explicit `S^-1` whose rounding
(about 1e-6 on an innovation covariance of 1e10) is far larger than its true
value (about 1e-10), and the reverse through `e^{JAx dt}` then multiplies that
direction by up to the square of the growth. Forward mode never forms that
cotangent. Found on the SNSF pilot: a person's drift with eigenvalues -1.05
and +0.28 over a 43-day gap, growth 2e5 -- the adjoint gradient off by up to
1.7e4 (relative) against BigFloat, the same adjoint in BigFloat exact, Float64
forward mode good to 1.6e-6. A stable drift's transition shrinks, so this
fires only on explosive (or strongly non-normal) intervals.

`Inf`, the reverse pass always, by default: forward mode needs a dual type the
fit has not otherwise used, so the first subject that takes it compiles the
filter again -- inside the Laplace route's seeded sweeps at a nested type. On a
one-latent model, 55 s on the plain route and 140 s on the Laplace one, once
per session (dev2); on the SNSF pilot's 32-effect model a single iteration
then ran 35 minutes on one core with the memory climbing past 25 GB. Opted
into per fit by optimcontrol's `explosive_forward`, at
`_CTSEM_EXPLOSIVE_GROWTH`. `ctsem_set_adjoint_growth!` moves it.
"""
const _CTSEM_ADJOINT_GROWTH = Ref(Inf)

"""
Transition growth past which a subject's filter pass counts as explosive: what
the progress line counts, what the fit's warning at the estimate and the
sampler's check of its draws report. Measured and reported whatever
`_CTSEM_ADJOINT_GROWTH` says, since the measurement costs nothing.
"""
const _CTSEM_EXPLOSIVE_GROWTH = Ref(100.0)

# Filter passes that crossed `_CTSEM_EXPLOSIVE_GROWTH`, since the session
# started; read as differences (`ctsem_explosive_passes`).
const _CTSEM_EXPLOSIVE_PASSES = Threads.Atomic{Int}(0)

"""Filter passes so far that crossed `_CTSEM_EXPLOSIVE_GROWTH`."""
ctsem_explosive_passes() = _CTSEM_EXPLOSIVE_PASSES[]

"""Set the transition growth above which subjects take the forward-mode gradient."""
function ctsem_set_adjoint_growth!(x::Real)
    x > 0 || throw(ArgumentError("threshold must be positive"))
    _CTSEM_ADJOINT_GROWTH[] = Float64(x)
    return _CTSEM_ADJOINT_GROWTH[]
end

"""
    _track_growth!(growth, E, n, subject)

Fold the transition `E[1:n, 1:n]` just applied into `growth` (`CTSEMGrowth`),
count the pass when it first exceeds `_CTSEM_EXPLOSIVE_GROWTH`, and note
`subject` then while `ctsem_explosive_subjects` is recording. A hand loop rather than
`_ctsem_mul!`: `E` may hold duals and only its primal values are wanted.
"""
function _track_growth!(growth::CTSEMGrowth, E::AbstractMatrix, n::Int, subject)
    phi, tmp = growth.phi, growth.tmp
    @inbounds if growth.steps == 0
        for j in 1:n, i in 1:n
            phi[i, j] = Float64(_primal(E[i, j]))
        end
    else
        for j in 1:n, i in 1:n
            acc = 0.0
            for k in 1:n
                acc += Float64(_primal(E[i, k])) * phi[k, j]
            end
            tmp[i, j] = acc
        end
        growth.phi, growth.tmp = tmp, phi
        phi = tmp
    end
    largest = 0.0
    @inbounds for j in 1:n
        column = 0.0
        for i in 1:n
            column += abs(phi[i, j])
        end
        largest = max(largest, column)
    end
    growth.steps += 1
    if largest > growth.max
        threshold = _CTSEM_EXPLOSIVE_GROWTH[]
        if growth.max <= threshold < largest
            Threads.atomic_add!(_CTSEM_EXPLOSIVE_PASSES, 1)
            _CTSEM_EXPLOSIVE_RECORDING[] && _note_explosive(subject)
        end
        growth.max = largest
    end
    return nothing
end

"""Whether the subject just filtered on `ws` needs the forward-mode gradient."""
@inline _growth_unstable(ws) = ws.growth.max > _CTSEM_ADJOINT_GROWTH[]

# Forward-mode subject gradients taken since the last reset, for the progress
# line (`_ctsem_forward_progress`) and the fit's result.
const _CTSEM_FORWARD_GRADIENTS = Threads.Atomic{Int}(0)

"""Forward-mode subject gradients taken since the last reset; `reset` zeroes the count."""
function ctsem_forward_gradients(; reset::Bool=false)
    return reset ? Threads.atomic_xchg!(_CTSEM_FORWARD_GRADIENTS, 0) :
        _CTSEM_FORWARD_GRADIENTS[]
end

const _CTSEM_EXPLOSIVE_RECORDING = Ref(false)
const _CTSEM_EXPLOSIVE_SUBJECTS = Set{Int}()
const _CTSEM_EXPLOSIVE_LOCK = ReentrantLock()

_note_explosive(subject) =
    lock(() -> push!(_CTSEM_EXPLOSIVE_SUBJECTS, Int(subject)), _CTSEM_EXPLOSIVE_LOCK)

"""
    ctsem_explosive_subjects(f)

Run `f()` and return its value with the sorted subject indices whose filter,
during it, exceeded `_CTSEM_EXPLOSIVE_GROWTH` -- the subjects whose predictions
grow that much between observations at the points `f` evaluates. Not
reentrant.
"""
function ctsem_explosive_subjects(f)
    lock(() -> empty!(_CTSEM_EXPLOSIVE_SUBJECTS), _CTSEM_EXPLOSIVE_LOCK)
    _CTSEM_EXPLOSIVE_RECORDING[] = true
    value = try
        f()
    finally
        _CTSEM_EXPLOSIVE_RECORDING[] = false
    end
    return value, lock(() -> sort!(collect(_CTSEM_EXPLOSIVE_SUBJECTS)),
        _CTSEM_EXPLOSIVE_LOCK)
end

"""
    _with_forward_count(f)

`f()`'s named tuple with `forward_gradients`, the forward-mode subject gradients
taken during it, and `explosive_passes`, the filter passes that crossed
`_CTSEM_EXPLOSIVE_GROWTH`, added: how a sampler reports them.
"""
function _with_forward_count(f)
    n0 = ctsem_forward_gradients()
    e0 = ctsem_explosive_passes()
    result = f()
    return merge(result, (forward_gradients=ctsem_forward_gradients() - n0,
        explosive_passes=ctsem_explosive_passes() - e0))
end

"""
    ctsem_explosive_draws(objective, draws)

Evaluate `objective` at each column of `draws` and report how many columns had
a subject whose filter exceeded `_CTSEM_EXPLOSIVE_GROWTH`, and which subjects
(`[0]` for none). A Laplace objective filters each subject at its effects' mode
given the column.
"""
function ctsem_explosive_draws(objective, draws::AbstractMatrix)
    hits = 0
    subjects = Set{Int}()
    for j in axes(draws, 2)
        _, s = ctsem_explosive_subjects(() ->
            ctsem_evaluate(objective, collect(Float64, view(draws, :, j)); gradient=false))
        isempty(s) && continue
        hits += 1
        union!(subjects, s)
    end
    return (checked=size(draws, 2), explosive=hits,
        subjects=isempty(subjects) ? [0] : sort!(collect(subjects)))
end

export ctsem_explosive_draws

"""
The progress line's counts since `forward0` and `explosive0`, each when any:
filter passes past `_CTSEM_EXPLOSIVE_GROWTH`, and forward-mode subject
gradients (optimcontrol's `explosive_forward`).
"""
function _ctsem_forward_progress(forward0::Int, explosive0::Int)
    e = ctsem_explosive_passes() - explosive0
    n = ctsem_forward_gradients() - forward0
    e > 0 || n > 0 || return ()
    n > 0 || return (@sprintf("explosive %d", e),)
    return (@sprintf("explosive %d", e), @sprintf("fwd-grad %d", n))
end

"""
    _ctsem_forward_subject_gradient!(g, subject_objective, values)

`g` set to the gradient of one subject's log likelihood at `values`, by forward
mode, on a workspace of its own (the subject's cached one can be shared by
sampler chains filtering the same subject). Behind a barrier, so the filter is
compiled for the new dual type only when a subject actually needs this. Works
at any element type, the nested duals of the Hessian and the Laplace sweeps
included: ForwardDiff's tag for this call is created after theirs.
"""
function _ctsem_forward_subject_gradient!(g::AbstractVector, so, values::AbstractVector)
    Threads.atomic_add!(_CTSEM_FORWARD_GRADIENTS, 1)
    _ctsem_barrier(_ctsem_forward_subject_gradient_impl!, g, so, values)
    return g
end

function _ctsem_forward_subject_gradient_impl!(g, so, values)
    f = x -> _ekf_run(_init_continuous_ekf_workspace(eltype(x), so.params), so, x,
        _ctsem_tipred_vector(so.tipreds, x), nothing)
    ForwardDiff.gradient!(g, f, values)
    return nothing
end

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
    tape.nbinaries = 0; tape.nstationaries = 0
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

# Gathering copies, for the rows a record keeps.
#
# `view(pars.LAMBDA, o, 1:n)` with `o` a *vector* of indices builds a SubArray
# that is not strided, and unlike a range-indexed view it does not stay on the
# stack: line-level tracking put 96 KB per 300 subject evaluations on one of
# these, and there are four per row. The view existed only to avoid
# materialising `pars.LAMBDA[o, 1:n]`, which is the right instinct -- copying
# straight out of the parent with the index achieves it without the view.
@inline function _tape_gather!(destination::Vector{T}, source, rows) where {T}
    resize!(destination, length(rows))
    @inbounds for k in eachindex(rows)
        destination[k] = source[rows[k]]
    end
    return destination
end

@inline function _tape_gather!(destination::Matrix{T}, source, rows, ncol::Int) where {T}
    size(destination) == (length(rows), ncol) ||
        return T[source[rows[k], j] for k in eachindex(rows), j in 1:ncol]
    @inbounds for j in 1:ncol, k in eachindex(rows)
        destination[k, j] = source[rows[k], j]
    end
    return destination
end

@inline function _tape_gather_square!(destination::Matrix{T}, source, rows) where {T}
    m = length(rows)
    size(destination) == (m, m) ||
        return T[source[rows[i], rows[j]] for i in 1:m, j in 1:m]
    @inbounds for j in 1:m, i in 1:m
        destination[i, j] = source[rows[i], rows[j]]
    end
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
@inline _record_stationary!(::Nothing, args...) = nothing

function _record_subject_values!(tape::CTSEMAdjointTape, subject_values)
    resize!(tape.subject_values, length(subject_values))
    copyto!(tape.subject_values, subject_values)
    return nothing
end

function _record_init!(tape::CTSEMAdjointTape{T}, pars, n::Int,
        popraw=nothing) where {T}
    index = (tape.ninits += 1)
    source = view(pars.T0VAR, 1:n, 1:n)
    popsource = popraw === nothing ? view(zeros(T, 0, 0), 1:0, 1:0) : popraw
    if index <= length(tape.inits)
        record = tape.inits[index]
        record.T0VAR = _tape_fill!(record.T0VAR, source)
        record.POPRAW = _tape_fill!(record.POPRAW, popsource)
    else
        push!(tape.inits, CTSEMInitRecord{T}(Matrix{T}(source),
            Matrix{T}(popsource)))
    end
    _tape_push!(tape, :init, index)
    return nothing
end

"""
    _record_stationary!(tape, ws, pars)

Record `_ctsem_stationary!` just after it ran: the mean is in `ws.state` and the
covariance still in `ws.diffusion_buffer.out`, where it left them.
"""
function _record_stationary!(tape::CTSEMAdjointTape{T}, ws, pars) where {T}
    index = (tape.nstationaries += 1)
    n = _val(ws.state_dim)
    naff = _val(ws.affine_buffer.dim)
    k = length(ws.diffusion_state_indices)
    drift = view(pars.DRIFT, 1:naff, 1:naff)
    mean = view(ws.state, 1:naff)
    X = view(ws.diffusion_buffer.out, 1:k, 1:k)
    diffusion = view(pars.DIFFUSION, 1:n, 1:n)
    if index <= length(tape.stationaries)
        record = tape.stationaries[index]
        record.DRIFT = _tape_fill!(record.DRIFT, drift)
        record.mean = _tape_fill!(record.mean, mean)
        record.X = _tape_fill!(record.X, X)
        record.DIFFUSION = _tape_fill!(record.DIFFUSION, diffusion)
    else
        push!(tape.stationaries, CTSEMStationaryRecord{T}(Matrix{T}(drift),
            Vector{T}(mean), Matrix{T}(X), Matrix{T}(diffusion)))
    end
    _tape_push!(tape, :stationary, index)
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
    record.Lambda = _tape_gather!(record.Lambda, pars.LAMBDA, o, n)
    _tape_gather!(record.manifestmeans, pars.MANIFESTMEANS, o)
    record.H = _tape_gather!(record.H, pars.Jy, o, n)
    record.R = _tape_gather_square!(record.R, ws.bufferΘ.out, o)
    resize!(record.y, length(o))
    @inbounds for j in eachindex(o)
        record.y[j] = T(data[o[j], obs_col])
    end
    _tape_push!(tape, :update, index)
    return nothing
end

"""An all-zero predict record shaped for `n` states and `k` dynamic states."""
# `k` sizes the Lyapunov block (the diffusing states), `naff` the affine one
# (the leading genuine dynamics). They are different sets; see `affine_dim`.
_empty_predict_record(::Type{T}, n::Int, k::Int, naff::Int=k) where {T} =
    CTSEMPredictRecord{T}(zeros(T, n), zeros(T, n, n), zeros(T, n, n),
        zeros(T, n, n), zeros(T, n, n), zeros(T, n, n), zeros(T, k, k),
        zeros(T, naff), zeros(T, naff), zero(T), false, false, zeros(T, 0, 0))

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
            _empty_predict_record(T, n, length(ws.diffusion_state_indices),
                _val(ws.affine_buffer.dim)))
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
    # The affine offset and the solved intercept live on the affine block, not
    # the diffusion one, so they come from `affine_buffer` and from the leading
    # entries of `dINT`. (Unused by the discrete-time reverse pass, which has
    # neither a solve nor a Lyapunov term.)
    naff = _val(ws.affine_buffer.dim)
    _tape_fill!(record.affine, view(ws.affine_buffer.r, 1:naff))
    _tape_fill!(record.dINT_dynamic, view(ws.discrete_ca.dINT, 1:naff))
    record.dt = T(Δt)
    series = ws.discretization_buffer.series
    record.series_intercept = series.intercept
    record.series_noise = series.noise
    if series.noise
        record.Qd = _tape_fill!(record.Qd, view(series.Q, 1:k, 1:k))
    end
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

The key, in `_frechet_enqueue!`, is an exact comparison against the actual
`JAx * dt` that produced each pending direction, so nothing here assumes
linearity: a state-dependent model changes `JAx` per row, matches nothing and
pays only the comparison, while a schedule with few distinct intervals
batches to one evaluation per distinct interval however its rows are ordered.

Correctness depends on nothing consuming the JAx cotangent while a flush is
outstanding. Two things can: a state-dependent transform group that *writes* a
JAx cell (whose reverse zeroes that cell's cotangent), and the parameter layer.
Both flush first.
"""
function _flush_frechet!(aws)
    count = aws.frechet_count
    count == 0 && return nothing
    aws.frechet_count = 0
    n = aws.n
    L = aws.frechet_L
    @inbounds for e in 1:count
        # `L(A', Ā)` by the Al-Mohy-Higham recurrence, into the workspace's
        # own buffers; the `_ctsem_exp_frechet_adjoint` wrapper does the same
        # computation but allocates, which at one flush per substep is what
        # the GC share of the reverse pass was made of.
        _CTSEM_OPCOUNT.frechet[] += 1
        my_exp_frechet!(aws.frechet_Y, L, transpose(aws.frechet_As[e]),
            aws.frechet_accums[e], aws.frechet_buffer)
        dt = aws.frechet_dts[e]
        if aws.defer_frechet
            _frechet_add!(aws.jax_bar_deferred, L, dt, n)
        else
            _frechet_add!(ComponentVector(aws.theta_bar, aws.sp.parameter_axis).JAx, L, dt, n)
        end
    end
    return nothing
end

"""`target .+= dt .* L` over the leading `n × n` block, without a temporary."""
@inline function _frechet_add!(target, L, dt, n::Int)
    @inbounds for j in 1:n, i in 1:n
        target[i, j] += dt * L[i, j]
    end
    return nothing
end

"""
Queue the direction `Ā` for the exponential of `scaled = JAx * dt`.

A direction whose key is already in the table is added to that entry; a new
key takes the next free slot, and a full table is flushed first. `dt` is
compared before the matrix, so a schedule of irregular intervals pays one
scalar comparison per live entry and no matrix comparison at all.
"""
function _frechet_enqueue!(aws, dt, scaled::AbstractMatrix, Ā::AbstractMatrix)
    @inbounds for e in 1:aws.frechet_count
        if aws.frechet_dts[e] == dt && aws.frechet_As[e] == scaled
            aws.frechet_accums[e] .+= Ā
            return nothing
        end
    end
    aws.frechet_count == length(aws.frechet_dts) && _flush_frechet!(aws)
    e = aws.frechet_count + 1
    aws.frechet_dts[e] = dt
    copyto!(aws.frechet_As[e], scaled)
    copyto!(aws.frechet_accums[e], Ā)
    aws.frechet_count = e
    return nothing
end

"""
    _reverse_predict!(x̄, P̄, θ̄ca, record, dyn, n, lyap_buffer, aws)

Undo one prediction substep.

Forward, on the dynamic index subset `D = dyn` of size `k`:

    A            = exp(JAx * dt)
    affine[i]    = CINT[i] + Σⱼ (DRIFT[i,j] - JAx[i,j]) x[j]      (i ≤ naff)
    s[i]         = -affine[i] + Σ_q A[i,q] affine[q]              (q ≤ naff)
    dINT[1:naff] = JAx[1:naff,1:naff] \\ s
    X            = lyap(JAx[D,D], Qc[D,D])
    dDIFF[D,D]   = X - A[D,D] X A[D,D]'
    x⁺           = A x + dINT
    P⁺           = A (P + εI) A' + dDIFF

Note `A` is the matrix exponential itself: this port sets `dDRIFT = eJAx`, so
the discrete drift and the exponential are the same object and share one
cotangent.

Where a pivot was small against the interval the forward pass took either of
two lines as the integral it is instead -- `dINT[1:naff] = Φ affine` with
`Φ = ∫₀^dt e^{JAx s} ds`, and `dDIFF[D,D] = ∫₀^dt e^{JAx s} Qc e^{JAx' s} ds`
-- and the record says which (`series_intercept`, `series_noise`). Neither
reads `A`, so on those lines nothing reaches `Ā`; see
`series_discretization.jl`.

Two different sub-blocks appear, and conflating them was a bug. `D = dyn` is
`diffusion_state_indices`, the states with their own diffusion; the Lyapunov
solve and the `dDIFF` term live there. `1:naff` is the leading block of genuine
dynamics -- everything that is not a static random-effect carrier -- and the
intercept solve and the affine offset live there. `dyn` is a subset of
`1:naff`, equal to it for any model with a full DIFFUSION, and strictly smaller
when a latent has no diffusion of its own and no coupling to one that has some.
Running the affine offset over `dyn` dropped such a state's CINT from the
likelihood entirely. `A` spans the whole augmented state throughout.
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
        return _reverse_predict_discrete!(x̄, P̄, θ̄ca, record, dyn, n;
            covmatcode=aws.sp.covmatcode)
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
    Ad         = _rs(sc.kk1, k, k)
    JAxd       = _rs(sc.kk2, k, k)
    Qb         = _rs(sc.kk3, k, k)
    X̄          = _rs(sc.kk4, k, k)
    Ād         = _rs(sc.kk5, k, k)
    kk6        = _rs(sc.kk6, k, k)
    x̄_new      = _rs(sc.xbar_new, n)
    dINT_bar   = _rs(sc.nv1, n)
    # The intercept solve and the affine offset run over the leading `naff`
    # states, a different and larger block than `dyn`; both are live at once,
    # so these cannot share the `kk`/`kv` slots.
    naff       = aws.affine_dim
    s̄          = _rs(sc.av1, naff)
    affine_bar = _rs(sc.av2, naff)
    aa1        = _rs(sc.aa1, naff, naff)

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
    @inbounds for j in 1:k, i in 1:k
        Qb[i, j] = (Ps[dyn[i], dyn[j]] + Ps[dyn[j], dyn[i]]) / 2
    end
    JAxd_bar = sc.kk7
    Qcd_bar = sc.kk8
    if record.series_noise
        # --- dDIFF[D,D] = V(JAx[D,D], Qc[D,D], dt), the integral itself; see
        # `series_discretization.jl`. It reads no exponential, so Ād is zero.
        fill!(Ād, zero(T))
        fill!(JAxd_bar, zero(T))
        fill!(Qcd_bar, zero(T))
        _series_noise!(sc.series, JAxd, record.Qd, record.dt, k)
        _series_noise_pullback!(JAxd_bar, Qcd_bar, sc.series, k, Qb)
    else
        X = record.Xlyap
        _ctsem_mulTN!(kk6, Ad, Qb)
        _ctsem_mul!(X̄, kk6, Ad)
        X̄ .= Qb .- X̄                                   # X̄ = Qb - Ad' Qb Ad
        _ctsem_mul!(kk6, Qb, Ad)
        _ctsem_mulNT!(Ād, kk6, X)
        _ctsem_mulTN!(kk6, Qb, Ad)
        _ctsem_mul!(Ād, kk6, X, one(T), one(T))
        Ād .= .-Ād                                     # Ād = -(Qb Ad X' + Qb' Ad X)

        # --- X = lyap(JAx[D,D], Qc[D,D])
        _ctsem_lyap_pullback!(JAxd_bar, Qcd_bar, sc.kk9, sc.kk10, JAxd, X, X̄, lyap_buffer)
    end

    if record.series_intercept
        # --- dINT[1:naff] = Phi affine, Phi = int_0^dt e^{JAx s} ds over the
        # leading block, taken as a series; see `series_discretization.jl`.
        # Like the noise series it reads no exponential, so nothing reaches Ā.
        Phi = _series_intercept!(sc.series, JAx, record.dt, naff)
        @inbounds for i in 1:naff
            acc = zero(T)
            for q in 1:naff
                acc += Phi[q, i] * dINT_bar[q]
            end
            affine_bar[i] = acc
        end
        @inbounds for j in 1:naff, i in 1:naff
            aa1[i, j] = dINT_bar[i] * record.affine[j]
        end
        _series_intercept_pullback!(JAx_bar, sc.series, naff, aa1)
    else
        # --- dINT[1:naff] = JAx[1:naff,1:naff] \ s   (s̄ = JAx⁻ᵀ dINT_bar, M̄ = -s̄ dINT')
        # The affine block is a leading one, so everything below indexes directly
        # with no gather, and its two outer products go straight into the full
        # cotangents rather than through a sub-block that then has to be scattered.
        @inbounds for i in 1:naff; s̄[i] = dINT_bar[i]; end
        # `transpose(JAx) \ s̄` would be a LAPACK `getrf!`/`getrs!` pair, once per
        # prediction substep, on a matrix of the dynamic-state size. At that size
        # the call is mostly OpenBLAS's process-global buffer lock -- see
        # `small_linalg.jl` -- so it goes through the engine's own LU instead.
        @inbounds for j in 1:naff, i in 1:naff; aa1[i, j] = JAx[j, i]; end
        _solve_square_system_generic!(aa1, s̄, sc.piv, naff)
        @inbounds for j in 1:naff, i in 1:naff
            JAx_bar[i, j] -= s̄[i] * record.dINT_dynamic[j]
        end

        # --- s = -affine + A affine, on the same leading block
        @inbounds for i in 1:naff
            acc = zero(T)
            for q in 1:naff
                acc += A[q, i] * s̄[q]
            end
            affine_bar[i] = acc - s̄[i]
        end
        @inbounds for j in 1:naff, i in 1:naff
            Ā[i, j] += s̄[i] * record.affine[j]
        end
    end

    # --- affine[i] = CINT[i] + Σⱼ (DRIFT[i,j] - JAx[i,j]) x[j]
    @inbounds for i in 1:naff
        θ̄ca.CINT[i] += affine_bar[i]
    end
    @inbounds for j in 1:n, i in 1:naff
        contribution = affine_bar[i] * x[j]
        θ̄ca.DRIFT[i, j] += contribution
        JAx_bar[i, j] -= contribution
        x̄_new[j] += (record.DRIFT[i, j] - JAx[i, j]) * affine_bar[i]
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
    _frechet_enqueue!(aws, record.dt, scaled, Ā)

    # --- Qc = sdcovsqrt2cov(DIFFUSION); only the dynamic block was consumed.
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Qcd_bar[i, j]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, θ̄ca, n)

    copyto!(x̄, x̄_new)
    copyto!(P̄, P̄_new)
    return nothing
end

"""
    _reverse_stationary!(x̄, P̄, θ̄ca, record, dyn, n, aws)

Undo `_ctsem_stationary!`. The forward pass overwrote the mean over the leading
`naff` states and every entry of those rows and columns of the covariance, so
their cotangents end here; the `:init` entry before it sees only what the
overwrite left, which in those rows is nothing.

  * mean, `m = -A \\ c` with `A = DRIFT[1:naff, 1:naff]`: `w = A' \\ m̄`, then
    `Ā -= w m'` and `c̄ -= w`.
  * covariance, `A_D X + X A_D' + Q_D = 0` over the diffusion states: the
    prediction's own Lyapunov pullback, and `Q̄` deferred through DIFFUSION's
    construction exactly as a prediction's is, keyed on the same raw matrix.
"""
function _reverse_stationary!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMStationaryRecord{T}, dyn::AbstractVector{Int}, n::Int,
    aws) where {T}
    sc = aws.reverse_scratch
    naff = aws.affine_dim
    k = length(dyn)

    w = _rs(sc.av1, naff)
    At = _rs(sc.aa1, naff, naff)
    @inbounds for i in 1:naff
        w[i] = x̄[i]
    end
    @inbounds for j in 1:naff, i in 1:naff
        At[i, j] = record.DRIFT[j, i]
    end
    _solve_square_system_generic!(At, w, sc.piv, naff)
    @inbounds for i in 1:naff
        θ̄ca.CINT[i] -= w[i]
    end
    @inbounds for j in 1:naff, i in 1:naff
        θ̄ca.DRIFT[i, j] -= w[i] * record.mean[j]
    end

    Ad = _rs(sc.kk1, k, k)
    X̄ = _rs(sc.kk4, k, k)
    @inbounds for j in 1:k, i in 1:k
        Ad[i, j] = record.DRIFT[dyn[i], dyn[j]]
        X̄[i, j] = P̄[dyn[i], dyn[j]]
    end
    Ād = sc.kk7
    Q̄d = sc.kk8
    _ctsem_lyap_pullback!(Ād, Q̄d, sc.kk9, sc.kk10, Ad, record.X, X̄,
        aws.lyap_buffer)
    @inbounds for j in 1:k, i in 1:k
        θ̄ca.DRIFT[dyn[i], dyn[j]] += Ād[i, j]
    end
    Qc_bar = _rs(sc.nn5, n, n)
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Q̄d[i, j]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, θ̄ca, n)

    @inbounds for i in 1:naff
        x̄[i] = zero(T)
    end
    @inbounds for j in 1:n, i in 1:naff
        P̄[i, j] = zero(T)
        P̄[j, i] = zero(T)
    end
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
    # The same size ternary the forward uses (`kalman_filters.jl`): above
    # `_CTSEM_SMALL_CHOLESKY` LAPACK's blocking earns its lock back, and using a
    # different factorization from the forward's is itself a way for the reverse
    # to decline an `S` the forward accepted.
    # One return type from both branches -- see the forward's copy of this in
    # `kalman_filters.jl` for why the `Union` cost what it did.
    F = if m <= _CTSEM_SMALL_CHOLESKY[]
        _ctsem_cholesky(mm2, m)
    else
        CTSEMCholesky(mm2, m, issuccess(cholesky!(mm2, check=false)))
    end
    if issuccess(F)
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
    # The seven broadcasts in this function were rewritten as @inbounds loops
    # once, on the evidence of a profile that put 9% of the gradient in
    # `_setindex!` and 4% in a `==` whose callers were broadcast's aliasing
    # check and its CartesianIndices iterator. Measured on dev1 over three runs
    # at one and ten threads: 1.5% on a one-latent model and nothing at all on
    # five or twelve latents. Most of that self time is not broadcast. Reverted,
    # and recorded here so it is not rediscovered.
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
            # applied binary observations before it. By the workspace's type,
            # so a model with no categorical indicator, whose tape never holds
            # one, does not compile the categorical reverse; see
            # `_ekf_categorical_call`.
            _ekf_categorical_call(aws.ekf_ws, _reverse_binary!, x̄, P̄, θ̄ca,
                tape.binaries[index], n)
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
            aws.groups_write_manifestvar && _flush_manifestvar!(aws, θ̄ca, m)
            aws.groups_write_diffusion && _flush_diffusion!(aws, θ̄ca, n)
            _reverse_group!(θ̄, x̄, tape.groups[index], sp, aws)
        elseif kind === :theta
            # Deferred rather than pushed back here; see `_defer_manifestvar!`.
            _defer_manifestvar!(aws, tape.thetas[index].MANIFESTVAR, Θ̄, θ̄ca, m)
            fill!(Θ̄, zero(T))
        elseif kind === :stationary
            _reverse_stationary!(x̄, P̄, θ̄ca, tape.stationaries[index], dyn, n, aws)
        elseif kind === :init
            @inbounds for i in 1:n
                θ̄ca.T0MEANS[i] += x̄[i]
            end
            fill!(x̄, zero(T))
            # The population block first, while P̄ still holds its
            # cotangent, then mask those rows and columns out so the T0VAR
            # pullback sees only what survived the overwrite. See the note on
            # `_apply_population_block!` in kalman_filters.jl for the forward
            # order this mirrors.
            popidx = aws.sp.population_indices
            if !isempty(popidx)
                k = length(popidx)
                # The forward pass wrote `k_a k_b block[a,b]` into P, so the
                # cotangent picks up the same factors on its way back. Missing
                # them leaves the likelihood right and the gradient wrong.
                popscale = aws.sp.population_scale
                popbar = zeros(T, k, k)
                @inbounds for b in 1:k, a in 1:k
                    popbar[a, b] = popscale[a] * popscale[b] *
                        P̄[popidx[a], popidx[b]]
                end
                popraw_bar = zeros(T, k, k)
                # The POPULATION code, not the model's: the forward pass
                # built this block with it, and a reverse pass differentiating
                # the other construction gives a gradient for a different
                # model than the likelihood.
                _sdcovsqrt2cov_pullback!(popraw_bar, tape.inits[index].POPRAW,
                    _symmetrized(popbar), k; scratch=aws.covsqrt_scratch,
                    covmatcode=aws.sp.population_covmatcode)
                # Written through the flat cotangent by range rather than by
                # name: `θ̄ca.RAWPOPVAR` would have to compile for models
                # whose axis has no such block, and cannot.
                rng = aws.sp.population_range
                popslot = reshape(view(θ̄, rng), k, k)
                @inbounds for b in 1:k, a in 1:k
                    popslot[a, b] += popraw_bar[a, b]
                end
                @inbounds for sidx in popidx
                    for j in 1:n
                        P̄[sidx, j] = zero(T)
                    end
                    for i in 1:n
                        P̄[i, sidx] = zero(T)
                    end
                end
            end
            # Its own workspace slot rather than `aws.covsqrt_scratch`, which
            # the pullback on the next line consumes.
            t0var_bar = _rs(aws.reverse_scratch.t0var_bar, n, n)
            fill!(t0var_bar, zero(T))
            sym = _symmetrize_into!(_rs(aws.reverse_scratch.sym, n, n), P̄)
            _sdcovsqrt2cov_pullback!(t0var_bar, tape.inits[index].T0VAR, sym, n;
                scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
            @inbounds for j in 1:n, i in 1:n
                θ̄ca.T0VAR[i, j] += t0var_bar[i, j]
            end
            fill!(P̄, zero(T))
        else
            throw(ArgumentError("adjoint: unknown tape entry $(kind)"))
        end
    end
    _flush_manifestvar!(aws, θ̄ca, m)
    _flush_diffusion!(aws, θ̄ca, n)

    return θ̄
end

"""
    _defer_manifestvar!(aws, mat, cov_bar, theta_bar_ca, m)

Hold this row's measurement-covariance cotangent instead of pushing it back.

`sdcovsqrt2cov`'s reverse is linear in the cotangent, so rows sharing a
MANIFESTVAR can sum their cotangents and pay for one pullback rather than one
each. Whether they share it is asked of the matrix itself rather than assumed
from the model: a state-dependent MANIFESTVAR changes between rows, and the
tape stores each row's own copy, so the comparison is exact and costs O(m^2)
against a pullback's O(m^3) plus a correlation square root.

The batch is released by `_flush_manifestvar!`, which the tape walk calls when
the matrix changes, before a transform group that writes MANIFESTVAR cells
(that group's reverse zeroes those cotangents, so a contribution arriving after
it would be attributed to the raw parameter instead of through the transform),
and at the end of the subject.
"""
function _defer_manifestvar!(aws, mat, cov_bar, θ̄ca, m::Int)
    if aws.mvar_pending && _ctsem_same_matrix(aws.mvar_mat, mat, m)
        @inbounds for j in 1:m, i in 1:m
            aws.mvar_bar[i, j] += cov_bar[i, j]
        end
        return nothing
    end
    _flush_manifestvar!(aws, θ̄ca, m)
    @inbounds for j in 1:m, i in 1:m
        aws.mvar_mat[i, j] = mat[i, j]
        aws.mvar_bar[i, j] = cov_bar[i, j]
    end
    aws.mvar_pending = true
    return nothing
end

"""
    _defer_diffusion!(aws, mat, cov_bar, theta_bar_ca, n)

`_defer_manifestvar!` for DIFFUSION, and the same argument applies with more
force: the prediction reverse runs once per *substep*, so a model with a
substep mesh pushes the same matrix back several times per row.
"""
function _defer_diffusion!(aws, mat, cov_bar, θ̄ca, n::Int)
    if aws.dvar_pending && _ctsem_same_matrix(aws.dvar_mat, mat, n)
        @inbounds for j in 1:n, i in 1:n
            aws.dvar_bar[i, j] += cov_bar[i, j]
        end
        return nothing
    end
    _flush_diffusion!(aws, θ̄ca, n)
    @inbounds for j in 1:n, i in 1:n
        aws.dvar_mat[i, j] = mat[i, j]
        aws.dvar_bar[i, j] = cov_bar[i, j]
    end
    aws.dvar_pending = true
    return nothing
end

"""Push back whatever `_defer_diffusion!` is holding, if anything."""
function _flush_diffusion!(aws, θ̄ca, n::Int)
    aws.dvar_pending || return nothing
    aws.dvar_pending = false
    _sdcovsqrt2cov_pullback!(θ̄ca.DIFFUSION, aws.dvar_mat, aws.dvar_bar, n;
        scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
    return nothing
end

"""Push back whatever `_defer_manifestvar!` is holding, if anything."""
function _flush_manifestvar!(aws, θ̄ca, m::Int)
    aws.mvar_pending || return nothing
    aws.mvar_pending = false
    # Straight into the cotangent: the pullback accumulates rather than
    # overwrites, and writes only the lower triangle and the diagonal, which is
    # exactly what the intermediate buffer used to be copied for.
    _sdcovsqrt2cov_pullback!(θ̄ca.MANIFESTVAR, aws.mvar_mat, aws.mvar_bar, m;
        scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
    return nothing
end

"""Elementwise equality over the leading `d` by `d` block."""
@inline function _ctsem_same_matrix(a::AbstractMatrix, b::AbstractMatrix, d::Int)
    @inbounds for j in 1:d, i in 1:d
        a[i, j] == b[i, j] || return false
    end
    return true
end

"""
    _ctsem_parameter_layer!(values_bar, θ̄, subject_values, sp, aws, tipreds, values)

Unwind the parameter layer: whatever cotangent is left on `all_params` was put
there by the regular (non-state-dependent) transforms, so push it through them
and then through the TI-predictor effects into `values_bar`.

`tipreds` is untyped rather than `::AbstractVector` because a subject with a
sampled (missing) TI predictor cell passes its `TIMissingRecipe` here instead
-- see `_ctsem_ti_pullback!`, which dispatches on it. `values` is the raw
trial point the forward pass ran at, needed only by that recipe method (the
plain-vector one ignores it).
"""
function _ctsem_parameter_layer!(values_bar::AbstractVector{T}, θ̄::AbstractVector{T},
    subject_values::AbstractVector, sp::EKFParameters, aws,
    tipreds, values::AbstractVector{T}) where {T}
    aws.defer_frechet || _flush_frechet!(aws)
    subject_values_bar = aws.subject_values_bar
    resize!(subject_values_bar, length(subject_values))
    fill!(subject_values_bar, zero(T))
    _ctsem_regular_pullback!(subject_values_bar, θ̄, subject_values, sp,
        aws.regular_supports, aws.regular_dual_scratch)
    _ctsem_ti_pullback!(values_bar, subject_values_bar, sp, tipreds, values)
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
    _reverse_predict_discrete!(XBAR, PBAR, THETA, record, dyn, n;
        covmatcode=sp.covmatcode)

Undo one prediction step of a discrete-time model.

Forward:

    A          = JAx
    dINT[i]    = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j]
    dDIFF[D,D] = Qc[D,D]
    x_next     = A x + dINT
    P_next     = A (P + eps I) A' + dDIFF

which is the continuous form with the exponential, the Lyapunov solve and the
intercept solve all collapsed -- a discrete model's DRIFT, CINT and DIFFUSION
are already the one-step quantities. The affine offset runs over every state,
matching the forward pass; `_compute_one_step_form!`'s docstring says why that
requires DRIFT's augmented diagonal to be 1 rather than 0.
"""
function _reverse_predict_discrete!(XBAR::Vector{T}, PBAR::Matrix{T}, THETA,
    record::CTSEMPredictRecord{T}, dyn::AbstractVector{Int}, n::Int;
    covmatcode::Int=0) where {T}
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
    _sdcovsqrt2cov_pullback!(diffusion_bar, record.DIFFUSION, Qc_bar, n;
        covmatcode=covmatcode)
    @inbounds for j in 1:n, i in 1:n
        THETA.DIFFUSION[i, j] += diffusion_bar[i, j]
    end

    copyto!(XBAR, xbar_new)
    copyto!(PBAR, Pbar_new)
    return nothing
end
