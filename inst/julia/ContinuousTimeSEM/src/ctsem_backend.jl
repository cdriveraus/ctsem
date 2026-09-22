"""
    CTSEMOptimisable

Anything `ctsem_optimize` can maximise: a prepared likelihood over a vector,
answering `ctsem_evaluate(objective, x; gradient, gradient_method)`.

There are three. `CTSEMObjective` is the marginal one, over the parameters
alone with the latent states integrated out by the filter. `CTSEMJointObjective`
(`state_sampling.jl`) is the state-explicit one, over the parameters *and* the
innovations that build the states. `CTSEMLaplaceObjective` (`laplace.jl`) is
the multilevel one, over the population parameters with each unit's random
effects at their own mode. The optimiser needs to know nothing about the
difference, and does not: it asks for a value and a gradient at a vector, and
all three answer.

The third answers that contract but does not yet go through
`ctsem_optimize`: it has its own copy of the driver, which is why the
convergence verdict below is shared rather than duplicated a fourth time. What
remains to merge is the `fg!` closure -- whose only semantic difference is that
a trial point where a unit's inner Newton did not converge is invalid, since
the objective is not defined away from the mode -- the trace's extra `inner`
column, the progress line's extra field, the call accounting, and five extra
result fields.

An abstract type rather than a `Union` because the second is defined several
files later -- a signature naming it here could not be parsed.
"""
abstract type CTSEMOptimisable end

"""
    CTSEMObjective(params, subject_starts, timesteps, data)

Prepared continuous-time likelihood for the R `ctsem` backend. `data` is
stored as manifest variables by observation. Subject starts are prepared by R,
which owns data ordering and predictor preprocessing. This type retains only
the numerical per-subject EKF workspaces needed for repeated evaluation.
"""
mutable struct CTSEMObjective{P,O} <: CTSEMOptimisable
    params::P
    subject_objectives::O
    # Lazily built reverse-mode workspace (see adjoint.jl). Cached here so the
    # one-off model inspection it does -- discovering each transform's read set
    # -- happens once per objective rather than once per gradient evaluation.
    adjoint_ws::Any
    # Prior specification, as index/scale pairs over the raw parameter vector.
    # Deliberately *data*, not logic: which raw parameters carry which prior is
    # ctsem semantics that the R side already knows (it reads `laplaceprior`,
    # `tipredeffectscale`, `rawpopsdbase` and the parameter ordering out of
    # `standata`), while all the engine needs to do is evaluate a normal
    # log-density and its derivative. Keeping the decision on the R side also
    # keeps one copy of it rather than three.
    prior_index::Vector{Int}
    prior_scale::Vector{Float64}
    prior_weight::Float64
    # The conditional-imputation term for sampled (missing) TI predictor
    # values: one Gaussian per missing cell, `values[ti_missing_parameter[k]]
    # ~ Normal(ti_missing_mu[k], ti_missing_sigma[k])`. `mu`/`sigma` are
    # pre-computed once on the R side from complete cases (see
    # SPEC-tipred-sampling.md, "Decisions") and held fixed here -- this is
    # data, exactly like `prior_index`/`prior_scale` above, and for the same
    # reason. Empty for every model with no missing TI predictor cells, which
    # is what keeps `_ctsem_ti_missing_loglik` free for them.
    ti_missing_parameter::Vector{Int}
    ti_missing_mu::Vector{Float64}
    ti_missing_sigma::Vector{Float64}
end

CTSEMObjective(params, subject_objectives) =
    CTSEMObjective(params, subject_objectives, nothing, Int[], Float64[], 1.0,
        Int[], Float64[], Float64[])

"""
    _ctsem_log_prior(objective, values)

The prior contribution to the log posterior, and its gradient.

Matches the generated Stan model's `model` block term for term:
`normal_lpdf(x / scale | 0, 1)`, i.e. `-x^2/(2 scale^2) - log(2pi)/2`, summed
over the parameters the R side flagged. Note the missing `-log(scale)`: Stan
takes the density *of the scaled quantity*, so that Jacobian term is absent
there too. It is a constant and so irrelevant to optimisation, but including it
would put this engine's log probability a constant away from Stan's, which is
exactly what the parity tests would then report as a mismatch.
"""
function _ctsem_log_prior(objective::CTSEMObjective, values::AbstractVector{T}) where {T}
    isempty(objective.prior_index) && return zero(T)
    total = zero(T)
    weight = objective.prior_weight
    @inbounds for k in eachindex(objective.prior_index)
        scaled = values[objective.prior_index[k]] / objective.prior_scale[k]
        total += -0.5 * scaled * scaled - 0.5 * log(2 * pi)
    end
    return weight * total
end

"""Accumulate the prior's gradient contribution into `gradient`."""
function _ctsem_log_prior_gradient!(gradient::AbstractVector{T},
    objective::CTSEMObjective, values::AbstractVector{T}, weight::Real=1.0) where {T}
    isempty(objective.prior_index) && return gradient
    scale_weight = objective.prior_weight * weight
    @inbounds for k in eachindex(objective.prior_index)
        idx = objective.prior_index[k]
        scale = objective.prior_scale[k]
        gradient[idx] -= scale_weight * values[idx] / (scale * scale)
    end
    return gradient
end

"""
    _ctsem_ti_missing_loglik(objective, values)

The conditional-imputation contribution to the log posterior: a proper
Gaussian log-density (the `-log(sigma)` term is *not* dropped, unlike
`_ctsem_log_prior` above -- these are genuine data values on their own scale,
not a standardised quantity, and the pre-computed `sigma` differs per
predictor, so the normalising constant is part of the answer here).

`_ctsem_ti_missing_loglik_gradient!`, right below, is its analytic derivative
-- both the adjoint path (`ctsem_adjoint_gradient`, which now supports a
sampled TI predictor value in full) and `:forward` (ForwardDiff
differentiating straight through this loop) reach this term; the adjoint calls
the gradient function once per evaluation rather than differentiating it.
Empty and free for every model with no missing TI predictor cells.
"""
function _ctsem_ti_missing_loglik(objective::CTSEMObjective, values::AbstractVector{T}) where {T}
    isempty(objective.ti_missing_parameter) && return zero(T)
    total = zero(T)
    @inbounds for k in eachindex(objective.ti_missing_parameter)
        idx = objective.ti_missing_parameter[k]
        sigma = objective.ti_missing_sigma[k]
        z = (values[idx] - objective.ti_missing_mu[k]) / sigma
        total += -0.5 * z * z - log(sigma) - 0.5 * log(2 * pi)
    end
    return total
end

"""
    _ctsem_ti_missing_loglik_gradient!(gradient, objective, values)

Accumulate `_ctsem_ti_missing_loglik`'s gradient contribution into `gradient`,
analytically: `d/d(values[idx]) [-0.5 z^2 - log(sigma) - 0.5 log(2 pi)]` with
`z = (values[idx] - mu) / sigma` is `-z / sigma`, i.e.
`-(values[idx] - mu) / sigma^2`. Each sampled cell gets its own raw-parameter
index (`.ctJuliaBackend.R` reserves one past every other block), so no two
`k` ever write the same `idx` -- but `+=` regardless, matching
`_ctsem_log_prior_gradient!`'s style, since summing is always safe and never
assumes that non-collision.

Mirrors `_ctsem_log_prior_gradient!` exactly: a closed-form function of the
raw parameters alone, so it is added once per evaluation rather than folded
into any per-subject pass. Empty and free for every model with no missing TI
predictor cells.
"""
function _ctsem_ti_missing_loglik_gradient!(gradient::AbstractVector{T},
    objective::CTSEMObjective, values::AbstractVector{T}) where {T}
    isempty(objective.ti_missing_parameter) && return gradient
    @inbounds for k in eachindex(objective.ti_missing_parameter)
        idx = objective.ti_missing_parameter[k]
        sigma = objective.ti_missing_sigma[k]
        gradient[idx] -= (values[idx] - objective.ti_missing_mu[k]) / (sigma * sigma)
    end
    return gradient
end

export CTSEMObjective, ctsem_objective, ctsem_evaluate, ctsem_optimize

"""The substep policy for one subject's rows `r`: a rule is shared, a mesh is sliced."""
_ctsem_subject_substeps(rule::Real, r) = rule
_ctsem_subject_substeps(mesh::AbstractVector, r) = Int[Int(m) for m in view(mesh, r)]

function _ctsem_row_loglikelihood(objective::CTSEMObjective, values::AbstractVector)
    contributions = Vector{eltype(values)}()
    for subject_objective in objective.subject_objectives
        previous = zero(eltype(values))
        nrows = size(subject_objective.data, 2)
        for row in 1:nrows
            prefix = ContinuousEKFObjective(objective.params,
                Matrix(subject_objective.data[:, 1:row]),
                collect(subject_objective.timesteps[1:row]);
                tdpreds=Matrix(subject_objective.tdpreds[:, 1:row]),
                tipreds=subject_objective.tipreds, subject=subject_objective.subject,
                max_timestep=_ctsem_subject_substeps(subject_objective.max_timestep, 1:row))
            current = prefix(values)
            push!(contributions, current - previous)
            previous = current
        end
    end
    contributions
end

function _ctsem_subject_ranges(subject_starts::AbstractVector, timesteps::AbstractVector,
    data::AbstractMatrix)
    n = length(timesteps)
    n > 0 || throw(ArgumentError("data must contain at least one observation"))
    n == size(data, 2) || throw(DimensionMismatch("timesteps length must equal data columns"))
    starts = Int.(subject_starts)
    !isempty(starts) || throw(ArgumentError("subject_starts must contain the first observation"))
    starts[1] == 1 || throw(ArgumentError("subject_starts must begin at one"))
    all((1 .<= starts) .& (starts .<= n)) || throw(ArgumentError("subject_starts are outside the observation range"))
    all(diff(starts) .> 0) || throw(ArgumentError("subject_starts must be strictly increasing"))
    stops = vcat(starts[2:end] .- 1, n)
    out = UnitRange{Int}[start:stop for (start, stop) in zip(starts, stops)]
    for r in out
        _validate_continuous_ekf_inputs(collect(view(timesteps, r)), Matrix(view(data, :, r)))
    end
    out
end

function CTSEMObjective(params::EKFParameters, subject_starts::AbstractVector,
    timesteps::AbstractVector, data::AbstractMatrix,
    tdpred_data::AbstractMatrix=zeros(eltype(data), 0, size(data, 2)),
    tipred_data::AbstractMatrix=zeros(eltype(data), length(subject_starts), 0),
    max_timestep=Inf; prior_index=Int[], prior_scale=Float64[],
    prior_weight::Real=1.0,
    # One entry per missing/sampled TI predictor cell, parallel arrays, empty
    # for every model with none (the overwhelmingly common case). `subject`
    # and `predictor` are consumed only here, to sort each cell into its
    # subject's `TIMissingRecipe`; `parameter`/`mu`/`sigma` are also kept flat
    # on the returned objective for `_ctsem_ti_missing_loglik`.
    ti_missing_subject=Int[], ti_missing_predictor=Int[], ti_missing_parameter=Int[],
    ti_missing_mu=Float64[], ti_missing_sigma=Float64[])
    ranges = _ctsem_subject_ranges(subject_starts, timesteps, data)
    # A mesh (one substep count per row of the whole dataset) is sliced per
    # subject exactly as `timesteps` is; a `maxtimestep` rule is shared.
    max_timestep isa AbstractVector && length(max_timestep) != length(timesteps) &&
        throw(DimensionMismatch("a substep mesh needs one entry per row of the data"))
    size(tdpred_data, 2) == size(data, 2) || throw(DimensionMismatch("TD predictor columns must match observations"))
    size(tipred_data, 1) == length(ranges) || throw(DimensionMismatch("TI predictor rows must match subjects"))
    nmissing = length(ti_missing_subject)
    (length(ti_missing_predictor) == nmissing && length(ti_missing_parameter) == nmissing &&
        length(ti_missing_mu) == nmissing && length(ti_missing_sigma) == nmissing) ||
        throw(DimensionMismatch("ti_missing_* vectors must all have equal length"))
    # Group missing cells by subject once, rather than scanning all of them
    # for every subject -- irrelevant at fit sizes seen so far, but O(n) not
    # O(n^2) costs nothing to keep that way.
    by_subject = Dict{Int,Vector{Int}}()
    for k in 1:nmissing
        push!(get!(() -> Int[], by_subject, Int(ti_missing_subject[k])), k)
    end
    # Copy each subject once. This avoids R proxy/view lifetime issues and makes
    # the objective safe to retain in a Julia session.
    #
    # A plain comprehension, not `Any[...]`. Every per-subject call in the
    # package goes through this vector, and an `Any` element type makes each of
    # them a dynamic dispatch -- which boxes the log likelihood coming back and
    # the `EKFParameters` going in, on every subject of every evaluation. The
    # comprehension narrows to the concrete subject type whenever the subjects
    # agree, which is every model whose subjects share a substep rule and none
    # of whose TI predictor cells are sampled; where they do not agree it
    # widens on its own and nothing is worse than it was.
    objects = [
        ContinuousEKFObjective(params, Matrix(view(data, :, r)), collect(view(timesteps, r));
            tdpreds=Matrix(view(tdpred_data, :, r)),
            # A subject absent from `by_subject` (i.e. every model with no
            # missing TI predictor cells) takes exactly the path it always
            # has: a plain `Vector{Float64}`, nothing wrapped around it.
            tipreds=if haskey(by_subject, i)
                ks = by_subject[i]
                TIMissingRecipe(vec(tipred_data[i, :]), Int.(ti_missing_predictor[ks]),
                    Int.(ti_missing_parameter[ks]))
            else
                vec(tipred_data[i, :])
            end,
            subject=i, max_timestep=_ctsem_subject_substeps(max_timestep, r))
        for (i, r) in enumerate(ranges)
    ]
    length(prior_index) == length(prior_scale) ||
        throw(DimensionMismatch("prior index and scale vectors must have equal length"))
    return CTSEMObjective(params, objects, nothing, Int.(prior_index),
        Float64.(prior_scale), Float64(prior_weight), Int.(ti_missing_parameter),
        Float64.(ti_missing_mu), Float64.(ti_missing_sigma))
end

ctsem_objective(params::EKFParameters, subject_starts, timesteps, data,
    tdpred_data=zeros(eltype(data), 0, size(data, 2)),
    tipred_data=zeros(eltype(data), length(subject_starts), 0), max_timestep=Inf;
    prior_index=Int[], prior_scale=Float64[], prior_weight::Real=1.0,
    ti_missing_subject=Int[], ti_missing_predictor=Int[], ti_missing_parameter=Int[],
    ti_missing_mu=Float64[], ti_missing_sigma=Float64[]) =
    CTSEMObjective(params, subject_starts, timesteps, data, tdpred_data, tipred_data,
        max_timestep; prior_index=prior_index, prior_scale=prior_scale,
        prior_weight=prior_weight, ti_missing_subject=ti_missing_subject,
        ti_missing_predictor=ti_missing_predictor, ti_missing_parameter=ti_missing_parameter,
        ti_missing_mu=ti_missing_mu, ti_missing_sigma=ti_missing_sigma)

################################################################################
# Threading
################################################################################
#
# The subject loop is the only parallelism here, and it is the natural one: the
# log-likelihood is a sum over subjects, each subject's `ContinuousEKFObjective`
# already owns its own primal and dual workspaces, and nothing is shared between
# them but the read-only `EKFParameters`.
#
# Work is split into contiguous *chunks* handed to `Threads.@spawn`, with each
# chunk owning the workspace it uses, rather than indexing workspaces by
# `threadid()`. That distinction matters: a task can migrate between threads at
# any yield point, so `threadid()` is not stable for the duration of a task and
# indexing mutable scratch by it is a data race waiting to happen.
#
# Threading changes the summation order, so a threaded result differs from a
# serial one at the last bits. `test_threading.jl` asserts agreement to 1e-12
# rather than bitwise, which is the correct gate for a floating-point reduction.

"""Cap on the number of chunks; 0 means "use `Threads.nthreads()`"."""
const _CTSEM_MAX_CHUNKS = Ref(0)

export ctsem_set_max_chunks!, ctsem_max_chunks
"""
    ctsem_set_max_chunks!(n)

Limit how many chunks the subject loop is split into. `0` (the default) means
use `Threads.nthreads()`, i.e. whatever the Julia process was started with.
Setting `1` forces the serial path, which is what the threading tests compare
against; the process-level thread count cannot be changed after startup.
"""
function ctsem_set_max_chunks!(n::Integer)
    n >= 0 || throw(ArgumentError("max chunks must be non-negative"))
    _CTSEM_MAX_CHUNKS[] = Int(n)
    return Int(n)
end

"""Current chunk cap, and the thread count it is resolved against."""
ctsem_max_chunks() = (max_chunks=_CTSEM_MAX_CHUNKS[], nthreads=Threads.nthreads())

@inline function _ctsem_nchunks(nsubjects::Int)
    requested = _CTSEM_MAX_CHUNKS[]
    available = requested == 0 ? Threads.nthreads() : min(requested, Threads.nthreads())
    return max(1, min(available, nsubjects))
end

"""
    ctsem_tune_chunks!(evaluate; ceiling=0, verbose=false)

Pick the chunk count that is actually fastest for this model, and pin it.

`cores` is a *ceiling*, not an instruction. Splitting the subject loop wider is
not free and it is not even monotone, and how far from monotone depends on the
model: a one-subject, one-latent evaluation once measured **3.7x slower on 23
threads than on one**, where a twenty-latent one ran 2.7x faster. What differs
is how much arithmetic sits between two allocations, so the ratio moves
whenever either side of that does.

It has moved. The reverse sweep used to allocate on the order of ten thousand
small arrays per subject and the collector was what the threads queued for;
after the allocation work it is about one per subject evaluation. Measured on
dev1 (23 cores, idle), 200 subjects of 24 rows, mean seconds per summed
gradient over thirty repeats, before and after that work:

                     1 thread        10 threads       MB/gradient   GC at 10
  1 latent    0.0186 -> 0.0172   0.0092 -> 0.0075    1.8 -> 0.1    14% -> 0%
  5 latent    0.0730 -> 0.0684   0.0125 -> 0.0099    4.8 -> 0.2    10% -> 0%
  12 latent   0.3628 -> 0.3235   0.0426 -> 0.0342   19.1 -> 0.3    29% -> 0%

So roughly a fifth off a ten-thread gradient, a tenth off a serial one, and
1-to-10 scaling from 8.5x to 9.5x on the twelve-latent model. The wider point
for this function is the last column and the spread behind it: the *minimum*
time barely moved, because a fast iteration was one with no collection in it,
while the mean and its run-to-run spread moved a lot. A tuner timing candidates
was measuring the collector as much as the division.

So the count cannot be chosen from `cores` alone, and it cannot be derived from
the model shape without a constant nobody has measured for the machine in
front of them. It can be measured, and cheaply: a handful of evaluations
against a fit that will do thousands. The ladder doubles from one, stops as
soon as two successive candidates fail to beat the best by a clear margin, and
leaves `_CTSEM_MAX_CHUNKS` at the winner.

`evaluate` is a zero-argument closure doing one representative evaluation --
whatever the caller's inner loop will actually be doing.
"""
function ctsem_tune_chunks!(evaluate; ceiling::Integer=0, verbose::Bool=false)
    limit = ceiling <= 0 ? _CTSEM_MAX_CHUNKS[] : Int(ceiling)
    limit = limit <= 0 ? Threads.nthreads() : min(limit, Threads.nthreads())
    if limit <= 1
        _CTSEM_MAX_CHUNKS[] = 1
        return (chunks=1, timings=[(1, NaN)])
    end
    previous = _CTSEM_MAX_CHUNKS[]
    # Each candidate spends the whole budget: `n` chunks of units, each dividing
    # a unit's members `limit / n` ways, so `n * width <= limit` always and the
    # ceiling the caller asked for is what runs. `n = limit` is the old
    # behaviour, chunks only; `n = 1` is members only.
    #
    # Both ends matter and neither dominates. A unit's innermost level is one
    # block per subject, which only chunks can spread; its outermost is a single
    # block over every member, which only the member axis can. And where units
    # are few and unequal -- thirteen studies, the largest 18.4% of the rows --
    # chunks alone cap the speedup at 5.4x however many cores are given.
    timings = Tuple{Int,Float64}[]
    best_chunks = 1
    best_time = Inf
    misses = 0
    try
        # One untimed pass at the serial setting so compilation, workspace
        # construction and the first mode solve are not charged to chunk 1.
        _CTSEM_MAX_CHUNKS[] = 1
        evaluate()
        serial = @elapsed evaluate()
        # A short evaluation cannot be timed once. Measured on a 100-subject
        # model whose evaluation is ~12 ms, three consecutive ladders gave chunk
        # 2 as 0.076 s, 0.0125 s and 0.0129 s -- a sixfold swing on the same
        # work, enough to pick a different winner each run. The minimum of a few
        # repeats is the right statistic under contention (a scheduler can only
        # ever make a timing longer), and three repeats of a 12 ms evaluation
        # cost nothing against the fit that follows.
        repeats = serial < 0.005 ? 5 : serial < 0.05 ? 3 : 1
        # The largest relative spread seen between repeats of one candidate,
        # which is this machine's noise floor for this evaluation, measured
        # rather than assumed. A candidate has to beat the incumbent by more
        # than this to be believed. Starts at 5%, the fixed margin this
        # replaces, so it can only ever become more demanding.
        noise = 0.05
        # The ladder is doubled from one, but only while an evaluation is cheap
        # enough that trying six settings is free against the fit. A model whose
        # evaluation takes a second is also a model with enough arithmetic
        # between allocations to thread well, so the ladder collapses to its
        # ends and the tuning costs four evaluations rather than a dozen.
        candidates = if serial > 1.0
            limit > 3 ? [1, max(2, limit ÷ 2), limit] : [1, limit]
        else
            built = Int[]; c = 1
            while c < limit; push!(built, c); c *= 2; end
            push!(built, limit); built
        end
        for n in candidates
            _CTSEM_MAX_CHUNKS[] = n
            evaluate()                       # warm this setting's workspaces
            elapsed = Inf
            slowest = 0.0
            for _ in 1:repeats
                one = @elapsed evaluate()
                elapsed = min(elapsed, one)
                slowest = max(slowest, one)
            end
            if repeats > 1 && elapsed > 0
                noise = max(noise, (slowest - elapsed) / elapsed)
            end
            push!(timings, (n, elapsed))
            verbose && println(_console(), "Chunk tuning: ", n, " worker(s) ",
                round(elapsed; digits=4), " s")
            # `1 - noise`, not a fixed 0.95. Splitting a small model's subject
            # loop finely is close to free either way, so the ladder was picking
            # between candidates that differ by less than the timing varies:
            # three identical 12-subject fits at `cores = 12` chose 12, 12 and
            # then 4. Whichever it lands on is then pinned for the whole fit, and
            # the wrong end of that is not harmless -- one chunk per subject is
            # where the allocator contention this tuner exists to avoid begins.
            #
            # Requiring a candidate to beat the incumbent by more than the
            # observed run-to-run spread means an unmeasurable difference cannot
            # decide the answer. Where nothing clears the bar the earliest
            # candidate stands, and the ladder climbs from one, so ties resolve
            # toward fewer chunks -- the conservative direction.
            if elapsed < best_time * (1 - noise)
                best_time = min(best_time, elapsed)
                best_chunks = n
                misses = 0
            else
                best_time = min(best_time, elapsed)
                misses += 1
                misses >= 2 && break
            end
        end
    catch err
        err isa InterruptException && rethrow()
        _CTSEM_MAX_CHUNKS[] = previous
        rethrow()
    end
    _CTSEM_MAX_CHUNKS[] = best_chunks
    verbose && println(_console(), "Chunk tuning: using ", best_chunks,
        " worker(s) of at most ", limit)
    return (chunks=best_chunks, timings=timings)
end

export ctsem_tune_chunks!

"""
    _ctsem_chunk_assignment(weights, nchunks)

Assign `1:length(weights)` to `nchunks` chunks, heaviest first into whichever
chunk is currently lightest.

Contiguous equal-count ranges are the right split when every element costs the
same, which is true of subjects and false of *units*: with subjects nested in
studies a unit is a study, and three studies of 200, 20 and 20 subjects split
three ways leave one chunk doing five times the work of the others while they
wait at the barrier. Sorting by cost first removes that, and for equal weights
it reproduces the same balance the contiguous split gives.

The chunks are returned as index vectors rather than ranges, since they are no
longer contiguous. Summation order across chunks changes with the assignment,
which is why `test_threading.jl` compares to a tolerance rather than bitwise.
"""
function _ctsem_chunk_assignment(weights::AbstractVector{<:Real}, nchunks::Int)
    n = length(weights)
    nchunks = max(1, min(nchunks, n))
    chunks = [Int[] for _ in 1:nchunks]
    nchunks == 1 && (append!(chunks[1], 1:n); return chunks)
    load = zeros(Float64, nchunks)
    for i in sortperm(weights; rev=true)
        c = argmin(load)
        push!(chunks[c], i)
        load[c] += max(0.0, Float64(weights[i]))
    end
    # Ascending within a chunk, so that a chunk still walks its units in the
    # order the serial loop would -- warm starts and diagnostics read better,
    # and nothing downstream depends on the order.
    for c in 1:nchunks; sort!(chunks[c]); end
    return chunks
end

"""Contiguous, near-equal partition of `1:n` into `nchunks` ranges."""
function _ctsem_chunk_ranges(n::Int, nchunks::Int)
    nchunks = max(1, min(nchunks, n))
    base, extra = divrem(n, nchunks)
    ranges = Vector{UnitRange{Int}}(undef, nchunks)
    start = 1
    @inbounds for c in 1:nchunks
        len = base + (c <= extra ? 1 : 0)
        ranges[c] = start:(start + len - 1)
        start += len
    end
    return ranges
end

function (objective::CTSEMObjective)(values::AbstractVector)
    subjects = objective.subject_objectives
    nsubjects = length(subjects)
    T = eltype(values)
    nchunks = _ctsem_nchunks(nsubjects)
    if nchunks <= 1
        total = zero(T)
        @inbounds for subject_objective in subjects
            total += subject_objective(values)
        end
        return total + _ctsem_log_prior(objective, values) +
               _ctsem_ti_missing_loglik(objective, values)
    end

    ranges = _ctsem_chunk_ranges(nsubjects, nchunks)
    partials = Vector{T}(undef, nchunks)
    Threads.@sync for c in 1:nchunks
        Threads.@spawn begin
            accumulator = zero(T)
            @inbounds for i in ranges[c]
                accumulator += subjects[i](values)
            end
            partials[c] = accumulator
        end
    end
    total = zero(T)
    @inbounds for c in 1:nchunks
        total += partials[c]
    end
    return total + _ctsem_log_prior(objective, values) +
           _ctsem_ti_missing_loglik(objective, values)
end

"""
Evaluate a prepared objective and optionally return its gradient.

`gradient_method` selects how the gradient is computed:

  * `:adjoint` -- the reverse-mode pass in `adjoint.jl`, and the default. One
    traced forward sweep plus one reverse sweep per subject, independent of
    the parameter count.
  * `:forward` -- ForwardDiff. Cost scales with the number of free parameters,
    since ForwardDiff needs one dual pass per chunk of them. Kept for
    cross-checking, and marginally quicker on very small nonlinear models.

The two are checked against each other (and against finite differences) by
`test/test_adjoint_gradient_validation.jl`. `:adjoint` never silently falls
back to `:forward`: an unsupported model or an invalid trial point produces a
non-finite gradient or a thrown error, not a quietly different answer.
"""
function ctsem_evaluate(objective::CTSEMObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    # Accepts a String as well as a Symbol: R passes this across the
    # JuliaConnectoR boundary, which marshals character vectors to `String`.
    method = Symbol(gradient_method)
    method in (:forward, :adjoint) ||
        throw(ArgumentError("gradient_method must be :forward or :adjoint, got :$(method)"))
    # A sampled TI predictor value is fully covered by the reverse pass now:
    # `_ctsem_ti_pullback!`'s `TIMissingRecipe` method (adjoint_parameters.jl)
    # supplies the product-rule term through `_materialize_subject_values!`'s
    # TI effect, and `ctsem_adjoint_gradient` adds
    # `_ctsem_ti_missing_loglik_gradient!`'s term for the imputation
    # log-density itself. See test_ti_missing_predictor.jl for the
    # cross-checks against ForwardDiff and FiniteDiff this replaced a refusal
    # with. `gradient_method=:forward` remains available, and identical.
    if gradient && method === :adjoint
        result = ctsem_adjoint_gradient(objective, collect(values))
        value = result.value
        grad = result.gradient
    else
        value = objective(values)
        grad = gradient ? ForwardDiff.gradient(objective, values) : nothing
    end
    if contributions
        subject = [obj(values) for obj in objective.subject_objectives]
        return (value=value, gradient=grad, subject_loglik=subject,
            row_loglik=_ctsem_row_loglikelihood(objective, values))
    end
    return (value=value, gradient=grad)
end

"""
How many curvature pairs L-BFGS keeps.

Optim's default is 10; ctsem's Stan path (`mize`) uses 100. Measured twice on a
100-subject, 50-wave model, once through the engine from zero starts and once
end to end through `ctFit`:

    m     engine, s   iters   |g| at stop   met g_tol      ctFit, s   iters
    10       4.08      42       3.5e-07        no            5.59      66
    20       3.73      38       4.7e-10        yes           5.35      57
    50       3.77      37       6.7e-10        yes           5.14      54
    100      3.82      37       6.7e-10        yes            --       --

Both runs reach the same optimum (log likelihood identical to four decimals)
and both say the same thing about direction: more memory, fewer iterations,
less time, and the gains flatten by twenty. They disagree about ten. The first
run stopped there with a gradient of 3.5e-07 against a `g_tol` of 1e-8 --
converged on the `f` test after the line search ran out, which is precisely the
false "converged" this route used to report -- and the second run met the
criterion at ten perfectly well. So "ten silently fails" is not a property of
the setting; it is a thing that *can* happen at ten and did not happen twice.

Twenty is the knee on both. The `lbfgs_memory` keyword overrides it, and R
reaches that through `optimcontrol`.
"""
const _CTSEM_LBFGS_MEMORY = 20

"""
Which coordinates of a minimizer carry ctsem parameter transforms, and so can
saturate. Every one of them, for the marginal objective; `state_sampling.jl`
gives its own answer.
"""
_ctsem_saturation_range(::CTSEMOptimisable, minimizer) = eachindex(minimizer)

"""
The `EKFParameters` behind an optimisable objective -- what
`_ctsem_saturated_parameters` (`parameter_transforms.jl`) needs to look up
each raw coordinate's materialising transform. `state_sampling.jl` gives its
own answer for `CTSEMJointObjective`, unwrapping to the same `CTSEMObjective`.
"""
_ctsem_params(o::CTSEMObjective) = o.params

"""
Which coordinates the overshoot probe pulls back.

* `:magnitude` -- prefixes of the coordinates ordered by `|raw|`, pulled back
  together. The default; see `_ctsem_overshot`.
* `:saturation` -- each coordinate the saturation detector flagged, one at a
  time. What this used to do, kept so the two can be compared.
* `:off` -- no probe. `overshot` is then always false, and a fit that walked
  into a flat region reports success, which is what the probe exists to stop.
"""
const _CTSEM_OVERSHOOT_PROBE = Ref(:magnitude)

"""
    ctsem_set_overshoot_probe!(mode)

Set the overshoot probe to `:magnitude`, `:saturation` or `:off`. Returns the
previous setting.
"""
function ctsem_set_overshoot_probe!(mode)
    previous = _CTSEM_OVERSHOOT_PROBE[]
    _CTSEM_OVERSHOOT_PROBE[] = _ctsem_overshoot_mode(mode)
    return previous
end

"""
    _ctsem_overshoot_mode(mode)

`mode` as one of the three symbols, or an error naming them.

Takes a string as well, because the R side passes one: a bare symbol does not
survive the bridge and spelling it as `Symbol(...)` at each call site is one
more place for `:magnitde` to be accepted silently.
"""
function _ctsem_overshoot_mode(mode)
    symbol = mode isa Symbol ? mode : Symbol(mode)
    symbol in (:magnitude, :saturation, :off) || throw(ArgumentError(
        "overshoot probe must be :magnitude, :saturation or :off, got $(repr(mode))"))
    return symbol
end

export ctsem_set_overshoot_probe!

"""
    _ctsem_probe_value(objective, x)

The objective at a probe point, or `-Inf` where the route cannot use it.

Untyped, because the fallback's contract is exactly "whatever answers
`ctsem_evaluate`" -- which is what the probe called directly before, and what
the mock objectives in `test_ctsem_backend.jl` supply. Anything that does not
answer raises inside the `try` and is refused, `nothing` included.

The laplace route overrides this: its value is only the objective at a mode, so
a probe point whose inner solve did not converge has to be refused rather than
compared. Without that the probe can be handed a number computed away from any
mode, find it larger, and report a maximum as an overshoot.
"""
function _ctsem_probe_value(objective, x)
    value = try
        ctsem_evaluate(objective, x; gradient=false).value
    catch
        -Inf
    end
    isfinite(value) ? value : -Inf
end

"""
What multiples of its current value each coordinate in a pulled-back set is
tried at.

Contractions and reflections. Contraction is what a saturated coordinate needs
-- raw 12.1 at 0.25 is raw 3.0, which is where that transform came back to life
and gained 142 nats -- but it can only ever reach zero, and zero is not a
neutral place to reach: it is correlation 0, drift -0.693 and sd 0.693,
depending entirely on the transform. So the ladder goes past it.

The negative rungs are there because a contraction cannot change a sign, and on
one measured fit the escape needed a population correlation to go from -2.358
to +0.327. `-1.0` is a pure reflection, which is the move a correlation pinned
at +1 needs if it belongs at -1.

Cheap either way: these are value-only evaluations, and more rungs buy a better
ranking at linear cost. See `_ctsem_overshot` for what the ranking is for.
"""
const _CTSEM_PULLBACK_FRACTIONS = (0.5, 0.25, 0.1, 0.0, -0.1, -0.25, -0.5, -1.0)

"""
    _ctsem_pullback_sets(minimizer)

The coordinate sets the magnitude probe pulls back: every prefix of the
coordinates sorted by `|raw|` descending.

There is no threshold in that, which is the point. A *cutoff* on raw magnitude
is the heuristic `ctsem_optimize`'s saturation note explains was abandoned,
because an identity transform never goes flat however large its coordinate. An
*ordering* by magnitude claims nothing about any coordinate -- it only decides
what to try first, and the objective decides whether it was right. A coordinate
wrongly included costs one evaluation.

Every prefix rather than a geometric ladder, and that is measured rather than
cautious. Two local optima of one 50-subject model, each with a population
scale collapsed and its correlations following it out:

    dataset 13   escaping set {6,9,10,11}   +62.2   prefix of size 4
    dataset 12   escaping set {6,9,10}      +53.7   prefix of size 3

A ladder of 1, 2, 4, 8 finds the first and jumps over the second. In both, an
exhaustive search over every subset of the six population coordinates found
nothing better than the prefix -- so the ordering is doing the work and the
only question was how finely to sample it.

`n` sets of `length(fractions)` evaluations, worst case, and it stops at the
first improvement.
"""
function _ctsem_pullback_sets(minimizer)
    n = length(minimizer)
    n == 0 && return Vector{Int}[]
    order = sortperm(collect(minimizer); by=abs, rev=true)
    return [order[1:k] for k in 1:n]
end

"""
    _ctsem_overshot(objective, minimizer, saturated_parameters, value, tolerance)

Whether the optimizer walked *past* an optimum into a flat region, rather than
stopping at one.

A transform gone flat, or a population scale collapsed to zero, covers two
outcomes that look identical from the gradient, because a flat transform
reports zero gradient in both:

  1. The optimizer overstepped. Measured on a binary model: one L-BFGS
     iteration to raw 20.9, log likelihood -730.7 where the profile peak is
     -714.5. That fit did not converge, and the point it stopped at is not a
     maximum.
  2. The data do not identify that coordinate, so the optimizer correctly ran
     it to the edge while every other parameter converged. A population
     standard deviation with no individual differences behind it is the common
     case, and it is a *finding*, not a failure -- one such fit stopped with a
     largest gradient of 5e-10 and matched Stan's log likelihood to the digit.

Treating both as "not converged" is what made the flag useless: over 64
optimisation replications of a benchmark whose log likelihoods matched Stan's,
45 reported `converged=false`, every one of them for a collapsed population
scale and the correlation that goes with it.

The two are told apart by asking the only question that separates them -- is
this a maximum? At a maximum no move improves the objective. So coordinates are
pulled back toward zero, into the region where their transforms still respond,
and an improvement larger than `tolerance` means the reported point is not a
maximum, which is case 1; no improvement is case 2. The probe is
self-validating in a way the saturation *flag* is not: a coordinate wrongly
included costs one evaluation, because zeroing a well-estimated parameter makes
the objective worse and nothing is reported.

## Why not one coordinate at a time

Because the move out of a degenerate corner is not along a coordinate. Measured
on a 50-subject model with three correlated random effects, at an optimum
149 nats below the one two other starts reached, with the T0MEANS population
scale collapsed to raw -15.6 and its correlation on the cap:

    pulled back alone     popsd_T0m_eta1            -0.0001 to -771
    pulled back alone     rawcor_drift__T0m_eta1      -18.7 to -657
    the two flagged ones, together                    -18.7 to -2065
    the four largest |raw|, together, to zero                 +62.2

Every single-coordinate probe says maximum. It is not one: the scale and all
three of its correlations have to come back at once, because a dimension with
no variance leaves its correlations unconstrained, and any one of them alone
still describes a degenerate covariance. Two of those four are what the
saturation detector flags -- the other two have transforms that are merely
unresponsive rather than flat -- so selecting the set from the flag cannot work
either. See `_ctsem_pullback_sets` for what is selected instead.

At most `length(fractions) * npar` value-only evaluations, once per fit, and it
stops at the first set that improves -- so `gain` is a lower bound on what is
left rather than the best pullback available. The question it answers is
whether the estimate is a maximum, and any improvement settles that.
"""
function _ctsem_overshot(objective, minimizer, saturated_parameters, value,
        tolerance; fractions=_CTSEM_PULLBACK_FRACTIONS,
        mode=_CTSEM_OVERSHOOT_PROBE[])
    gain = 0.0
    empty_result = (overshot=false, gain=gain, coordinates=Int[], point=Float64[])
    (!isfinite(value) || isempty(minimizer)) && return empty_result
    mode = _ctsem_overshoot_mode(mode)
    mode === :off && return empty_result
    sets = if mode === :saturation
        [[p] for p in saturated_parameters if 1 <= p <= length(minimizer)]
    else
        _ctsem_pullback_sets(minimizer)
    end
    isempty(sets) && return empty_result
    keep = collect(minimizer)
    probe = collect(minimizer)
    # No ranking of the points that did not improve. Ordering them by how
    # little they lost is a proxy for nothing: the least bad candidate is the
    # one that disturbed the estimate least, which makes it the likeliest to
    # fall straight back into the basin a caller is trying to leave. What to do
    # when nothing improves is decided by which coordinate is on a boundary,
    # not by which near miss was nearest -- see `.ctBackendStallEscape()`.
    for set in sets
        for f in fractions
            for p in set; probe[p] = f * keep[p]; end
            trial = _ctsem_probe_value(objective, probe)
            isfinite(trial) && (gain = max(gain, trial - value))
            if gain > tolerance
                # The point as well as the verdict. It costs a copy and it is
                # what a caller needs to *act* on an overshoot rather than only
                # report it -- see `ctsem_pullback`.
                return (overshot=true, gain=gain, coordinates=sort(set),
                    point=copy(probe))
            end
        end
        for p in set; probe[p] = keep[p]; end
    end
    return (overshot=false, gain=gain, coordinates=Int[], point=Float64[])
end

"""
    ctsem_pullback(objective, values; tolerance, probe)

The pulled-back point that improves on `values`, or nothing to report.

`_ctsem_overshot` answers "is this a maximum" and finds an improving point on
the way to saying no. This hands that point back, so a fit that stopped in a
degenerate corner can be resumed from somewhere better instead of only being
told it is stuck. The improvement is measured, not predicted, so a resume from
here cannot start worse than it stopped.

Returned as a named tuple rather than `nothing` for the empty case, because a
zero-length vector deadlocks the R bridge: `found` is the flag to read.
"""
function ctsem_pullback(objective, values::AbstractVector; tolerance::Real=1e-6,
        probe=_CTSEM_OVERSHOOT_PROBE[])
    x = collect(Float64, values)
    nothing_found = (found=false, gain=0.0, point=x, coordinates=[0])
    current = _ctsem_probe_value(objective, x)
    isfinite(current) || return nothing_found
    out = _ctsem_overshot(objective, x, Int[], current, tolerance; mode=probe)
    out.overshot && !isempty(out.point) || return nothing_found
    return (found=true, gain=out.gain, point=collect(Float64, out.point),
        coordinates=isempty(out.coordinates) ? [0] : out.coordinates)
end

export ctsem_pullback

"""
A line search that records the directional derivative it is handed.

`LineSearches` receives `dphi0 = g'p` as its last positional argument, every
iteration, because it needs it to test the Wolfe conditions. For a quasi-Newton
direction `-dphi0 / 2` is the improvement that step was predicted to make, in
objective units -- so observing it here gives a stopping rule on the quantity
that matters, for no arithmetic at all.

Deliberately an observer: it forwards every call unchanged, so which line
search actually runs is unaffected and `linesearch` still reports what it
always did.
"""
mutable struct CTSEMDirectional{LS}
    inner::LS
    dphi0::Float64
end

CTSEMDirectional(inner) = CTSEMDirectional(inner, NaN)

function (ls::CTSEMDirectional)(args...)
    # Last positional argument in both of LineSearches' call forms.
    last = args[end]
    ls.dphi0 = last isa Real ? Float64(last) : NaN
    ls.inner(args...)
end

"""The improvement the last accepted step was predicted to make, or `Inf`."""
function _ctsem_predicted_gain(ls::CTSEMDirectional)
    isfinite(ls.dphi0) ? abs(ls.dphi0) / 2 : Inf
end

"""
The per-objective pieces of the shared optimiser.

`ctsem_optimize` drives every route. What the routes differ by is here, as one
small generic function each, so that anything added to the optimiser is added
once. Three things were not, and each cost a bug: the convergence verdict was
fixed on one route and left on the other, the saturation guard had to be copied
across afterwards, and `gap_tol` went into one and every laplace fit died with a
MethodError.

The defaults are the marginal and joint routes' behaviour, so those two need no
methods at all; `laplace.jl` supplies the differences.

* `_ctsem_optimise_label` names the route in its messages.
* `_ctsem_optimise_setup!` resets whatever per-run state the route keeps.
* `_ctsem_optimise_log` returns a mutable accumulator for the route's own call
  accounting, or `nothing` when it keeps none. It is told `verbose`, because a
  route that narrates its accounting as it goes has to know whether anyone
  asked.
* `_ctsem_optimise_trace_keys` and `_ctsem_optimise_trace_values` are the
  columns the per-iteration trace carries: keys are fixed at construction so the
  vectors stay type-stable, values are read each iteration.
* `_ctsem_optimise_progress_extra` is anything further the progress line should
  carry, already formatted.
* `_ctsem_optimise_trial` evaluates a trial point and says whether it may be
  used. This is the one hook with real semantics rather than plumbing: on the
  laplace route a point where a unit's inner Newton did not reach its tolerance
  is *invalid* even though its value is finite, because the objective is only
  defined at the mode -- accepting it lets the outer optimiser follow a function
  of something other than theta.
* `_ctsem_saturated_for` detects a flat transform, which the laplace route does
  differently: a level of correlations can saturate together, which the per-cell
  check cannot see.
* `_ctsem_optimise_result_extra` adds the route's own result fields.
* `_ctsem_optimise_verbose_shape` and `_ctsem_optimise_verbose_report` are the
  route's own `verbose` lines, before the run and after it: what is about to be
  fitted, and what the run did. The second is where the call accounting in its
  log gets read.
"""
_ctsem_optimise_label(::CTSEMOptimisable) = "ctsem_optimize"
_ctsem_optimise_setup!(::CTSEMOptimisable) = nothing
_ctsem_optimise_log(::CTSEMOptimisable, verbose::Bool) = nothing
_ctsem_optimise_trace_keys(::CTSEMOptimisable) = (:objective, :gradient_norm)
_ctsem_optimise_trace_values(::CTSEMOptimisable) = ()
_ctsem_optimise_progress_extra(::CTSEMOptimisable) = ()
_ctsem_optimise_verbose_shape(::CTSEMOptimisable) = nothing
_ctsem_optimise_verbose_report(::CTSEMOptimisable, log) = nothing

_ctsem_saturated_for(o::CTSEMOptimisable, minimizer) =
    _ctsem_saturated_parameters(_ctsem_params(o), minimizer,
        _ctsem_saturation_range(o, minimizer))

"""The row-wise contributions only the filter routes can decompose."""
_ctsem_optimise_result_extra(::CTSEMOptimisable, final, log) =
    (row_loglik=final.row_loglik,)

"""
    _ctsem_optimise_trial(objective, x, want_gradient, gradient_method, limit, log)

One trial point: what it evaluated to, and whether the optimiser may use it.

`nothing` for the evaluation, or `valid == false`, both mean the same thing to
the caller -- hand back the sentinel objective and a zero gradient so the line
search shrinks. Deliberately not the last valid gradient: that feeds L-BFGS a
secant pair whose gradient never belonged to the point, which corrupts the
curvature history.
"""
function _ctsem_optimise_trial(o::CTSEMOptimisable, x, want_gradient::Bool,
        gradient_method, limit::Real, log)
    evaluated = try
        ctsem_evaluate(o, x; gradient=want_gradient,
            gradient_method=gradient_method)
    catch
        nothing
    end
    valid = evaluated !== nothing && isfinite(evaluated.value)
    if valid && want_gradient
        valid = all(isfinite, evaluated.gradient) &&
            all(abs(value) < limit for value in evaluated.gradient)
    end
    return (evaluated=evaluated, valid=valid)
end

"""
    CTSEMPinnedObjective(objective, index, value)

`objective` with the coordinates in `index` held at `value`.

A profile point and an escape attempt are the same operation, which is why this
exists once and has two callers. Both fix a coordinate somewhere other than
where the optimiser left it and re-optimise everything else; the difference is
only what the answer is used for. Profiling reads the constrained maximum as a
statement about identification; escaping releases the coordinate afterwards and
keeps the result if it beats where the fit was.

Pinning is what separates this from resuming a fit from a displaced point,
which the escape loop already does. An unpinned resume can slide straight back
down the direction it was pushed along -- measured on this package at -3721.05
against -2950.59 -- because nothing stops the coordinate returning to the basin
while the rest of the model stays put. Holding it still forces the other
parameters to accommodate the displaced value first, and only then is it let
go.

The parameter vector keeps its full length and the pinned entries are
overwritten on the way in, rather than optimising a shorter vector. Every index
in the engine -- `matsetup` rows, transform lookups, saturation and pullback
coordinate sets, the preconditioner -- is positional in the full vector, so a
reduced vector would need all of them remapped, and a missed one would not
error. It would return a number.

The gradient comes back zero in the pinned entries, which is what actually
holds them: L-BFGS builds its direction from gradients and secant pairs, and a
coordinate contributing zero to both keeps whatever it started with. The caller
should still overwrite the pinned entries of the minimizer, because "cannot
move in exact arithmetic" is not the same claim as "did not move".
"""
struct CTSEMPinnedObjective{O} <: CTSEMOptimisable
    objective::O
    index::Vector{Int}
    value::Vector{Float64}
end

"""
    ctsem_pin(objective, index, value)

`objective` with `index` pinned at `value`; see `CTSEMPinnedObjective`.

Validated here rather than at the first evaluation: a misspelled index should
cost nothing, not a whole optimisation that silently pinned the wrong
coordinate.
"""
function ctsem_pin(objective::CTSEMOptimisable, index::AbstractVector,
        value::AbstractVector)
    idx = collect(Int, index)
    val = collect(Float64, value)
    length(idx) == length(val) ||
        throw(ArgumentError("ctsem_pin: index and value must be the same length"))
    allunique(idx) ||
        throw(ArgumentError("ctsem_pin: each coordinate may be pinned only once"))
    all(isfinite, val) ||
        throw(ArgumentError("ctsem_pin: a pinned value must be finite"))
    all(>=(1), idx) ||
        throw(ArgumentError("ctsem_pin: coordinates are 1-based"))
    return CTSEMPinnedObjective(objective, idx, val)
end

export ctsem_pin

"""The trial point as the inner objective sees it: pinned entries restored."""
function _ctsem_pin_expand(p::CTSEMPinnedObjective, x::AbstractVector)
    y = collect(x)
    @inbounds for (position, i) in enumerate(p.index)
        1 <= i <= length(y) || continue
        y[i] = p.value[position]
    end
    return y
end

"""An evaluation with the pinned coordinates' gradient entries removed."""
function _ctsem_pin_project(p::CTSEMPinnedObjective, evaluated)
    evaluated === nothing && return nothing
    hasproperty(evaluated, :gradient) || return evaluated
    gradient = evaluated.gradient
    gradient === nothing && return evaluated
    g = collect(gradient)
    @inbounds for i in p.index
        1 <= i <= length(g) || continue
        g[i] = zero(eltype(g))
    end
    return merge(evaluated, (gradient=g,))
end

# Everything else is the inner objective's. The laplace route overrides most of
# this protocol -- its own trial validity, its own trace columns, its own
# verbose report -- and wrapping it must not quietly return any of that to the
# generic default. In particular `_ctsem_optimise_trial` is forwarded rather
# than reimplemented, so a pinned laplace fit still refuses a point whose inner
# Newton did not converge.
_ctsem_optimise_label(p::CTSEMPinnedObjective) =
    string(_ctsem_optimise_label(p.objective), " (pinned)")
_ctsem_optimise_setup!(p::CTSEMPinnedObjective) = _ctsem_optimise_setup!(p.objective)
_ctsem_optimise_log(p::CTSEMPinnedObjective, verbose::Bool) =
    _ctsem_optimise_log(p.objective, verbose)
_ctsem_optimise_trace_keys(p::CTSEMPinnedObjective) =
    _ctsem_optimise_trace_keys(p.objective)
_ctsem_optimise_trace_values(p::CTSEMPinnedObjective) =
    _ctsem_optimise_trace_values(p.objective)
_ctsem_optimise_progress_extra(p::CTSEMPinnedObjective) =
    _ctsem_optimise_progress_extra(p.objective)
_ctsem_optimise_verbose_shape(p::CTSEMPinnedObjective) =
    _ctsem_optimise_verbose_shape(p.objective)
_ctsem_optimise_verbose_report(p::CTSEMPinnedObjective, log) =
    _ctsem_optimise_verbose_report(p.objective, log)
_ctsem_optimise_result_extra(p::CTSEMPinnedObjective, final, log) =
    _ctsem_optimise_result_extra(p.objective, final, log)
_ctsem_params(p::CTSEMPinnedObjective) = _ctsem_params(p.objective)

"""
A pinned coordinate cannot saturate, overshoot or stall, because it cannot
move. Leaving it in the range would let the pullback probe pick it up and
report a gain from moving something this objective is holding still.
"""
_ctsem_saturation_range(p::CTSEMPinnedObjective, minimizer) =
    [i for i in _ctsem_saturation_range(p.objective, minimizer) if !(i in p.index)]

"""
Forwarded rather than left to the generic default, which would rebuild the
answer from `_ctsem_params` and the range. The laplace route overrides this and
its override is not reconstructible from those: it is the only thing that knows
a population correlation has reached its cap, which is exactly the state that
stalled the fit the stall check was written for. Falling back here would have
lost that silently on every pinned laplace stage.
"""
_ctsem_saturated_for(p::CTSEMPinnedObjective, minimizer) =
    [i for i in _ctsem_saturated_for(p.objective, minimizer) if !(i in p.index)]

function _ctsem_optimise_trial(p::CTSEMPinnedObjective, x, want_gradient::Bool,
        gradient_method, limit::Real, log)
    trial = _ctsem_optimise_trial(p.objective, _ctsem_pin_expand(p, x),
        want_gradient, gradient_method, limit, log)
    return (evaluated=_ctsem_pin_project(p, trial.evaluated), valid=trial.valid)
end

_ctsem_probe_value(p::CTSEMPinnedObjective, x) =
    _ctsem_probe_value(p.objective, _ctsem_pin_expand(p, x))

ctsem_evaluate(p::CTSEMPinnedObjective, x::AbstractVector; kwargs...) =
    _ctsem_pin_project(p, ctsem_evaluate(p.objective, _ctsem_pin_expand(p, x);
        kwargs...))

"""
    _ctsem_metric(precondition, n)

A diagonal metric for L-BFGS, or `nothing` to leave it alone.

L-BFGS starts with a *scalar* initial inverse Hessian and takes a first step of
unit length in the raw coordinates, so it has no per-coordinate scaling until
secant pairs accumulate. When one coordinate maps to model quantities ten times
faster than the others -- which a `meanscale` of 10 against 1 is exactly -- the
step that suits the rest is ten times too long for it, and the line search has
to shrink the *whole* step to accommodate the worst one. Measured on a
40-subject count model's laplace route: gradient components spanning 210:1
across coordinates, 30 function evaluations spent on 4 iterations, stopping with
303 nats still available and not one trial point rejected. At 26:1 the same fit
converges.

`precondition[i]` is `|d value / d raw|` for coordinate `i`, so `P` is its
square: `P` stands in for the Hessian, and a transform with factor `s`
multiplies second derivatives by `s^2`. A step is then the same amount of model
in every coordinate, which is the metric L-BFGS would have built for itself
after enough iterations to be worth having.

Non-finite or non-positive entries become 1 -- no preconditioning for that
coordinate -- rather than disqualifying the whole metric, since a single
un-differentiable transform should not cost every other parameter its scaling.
"""
function _ctsem_metric(precondition, n::Integer)
    precondition === nothing && return nothing
    scale = collect(Float64, precondition)
    length(scale) == n || return nothing
    @inbounds for i in eachindex(scale)
        (isfinite(scale[i]) && scale[i] > 0) || (scale[i] = 1.0)
    end
    all(isequal(1.0), scale) && return nothing
    Diagonal(scale .^ 2)
end

"""What the whole run has gained, and what its last `window` iterations did.

Only for reporting: `_ctsem_stalled` computes both itself rather than calling
these, so the rule and the message cannot drift apart on a rounding.
"""
function _ctsem_trace_progress(trace::CTSEMTrace)
    values = get(trace.values, :objective, Float64[])
    length(values) < 2 ? 0.0 : values[end] - values[1]
end

function _ctsem_window_gain(trace::CTSEMTrace, window::Integer)
    values = get(trace.values, :objective, Float64[])
    length(values) > window ? values[end] - values[end - window] : 0.0
end

"""
What the final iteration gained, in objective units.

`Inf` when there is no pair of iterations to compare, so a run that recorded one
row or none cannot satisfy a convergence test with it. Read from the trace
rather than tracked separately because the trace already records every
iteration's objective and a second accumulator of the same numbers is a second
thing to keep in step.
"""
function _ctsem_last_gain(trace::CTSEMTrace)
    values = get(trace.values, :objective, Float64[])
    length(values) < 2 && return Inf
    abs(values[end] - values[end - 1])
end

"""
    _ctsem_stalled(trace, window, fraction)

Whether the last `window` iterations gained a negligible share of the progress
the fit has made.

Half of a conjunction, and deliberately the weak half. On its own it cannot
tell a fit that is stuck from one that is converging slowly, and those want
opposite treatment -- a slow fit should be left alone, because the estimate has
to settle properly before its curvature is worth anything. So this only decides
*when to look*, and `_ctsem_flat_coordinates` decides whether there is anything
to find. A false positive here costs one cheap derivative pass.

The share of progress made, rather than a number of nats or a share of what is
predicted to remain. Nats are not comparable across models. What remains is
estimated in flight by `1/2 g'Bg`, and that is the one quantity that goes wrong
exactly here: measured on the fit this was written for it decayed nine orders
and then bounced back two, because a direction going flat sends the gradient to
zero and `B` to infinity together. Progress already made is neither, and it
carries its own safety -- a fit that started near its optimum has little of it,
so the bar is small and this cannot fire.

`false` before there is a window to look at.
"""
function _ctsem_stalled(trace::CTSEMTrace, window::Integer, fraction::Real)
    values = get(trace.values, :objective, Float64[])
    window >= 1 || return false
    length(values) > window || return false
    progress = values[end] - values[1]
    isfinite(progress) && progress > 0 || return false
    gained = values[end] - values[end - window]
    isfinite(gained) || return false
    return gained <= fraction * progress
end

"""
The conjunction that decides a fit has stopped for a reason, and the hysteresis
that keeps it from asking the same question every iteration.

`_ctsem_stalled` says the fit is not getting anywhere. `_ctsem_flat_coordinates`
says whether a transform has gone flat, which is the structural reason a fit
stops getting anywhere. Either alone is a mistake this codebase has already
made: saturation alone reported 45 of 64 good fits as failures, because a
population scale with no individual differences behind it saturates early and
legitimately while the fit goes on to a perfectly good optimum; and a progress
test alone cannot tell stuck from slow.

When the progress test fires and nothing is flat, the fit is given more rope
rather than asked again immediately: `cooldown` iterations of quiet, and the bar
tightened by `tighten`, at most `tightenings` times so it cannot drift somewhere
unprincipled. A fit that is merely slow therefore costs a handful of derivative
passes over its whole run and is then left alone, which is what a long careful
optimisation needs.

## Stalled and flat is not enough

It was, and `test_state_sampling.jl`'s count model over the joint density is
why it is not. That fit *correctly* ends saturated: its drift arrives at raw
-18.5 where `-log1p_exp` is flat, the objective is at its supremum there, and
pulling the coordinate back finds nothing better. Stalled-and-flat is exactly
what a fit looks like while it converges *into* a flat region, so the two
halves alone stopped it before it arrived and turned a converged fit into a
failed one.

What separates them is already being computed: whether a pullback finds
anything. The runaway drift this was built for gains 142 nats on the spot; the
count model gains nothing. So the third condition is that there is somewhere
better to go, and the fit is only stopped when all three hold -- which is also
what makes stopping safe, because the point to resume from is in hand.

A fit that is stalled and flat with nothing better nearby is treated as the
slow case: cooldown and tighten. If it is nonetheless in the wrong place, the
zero-and-refit escape after the fit is what finds out, because only a refit
can.
"""
mutable struct CTSEMStallWatch
    window::Int
    fraction::Float64
    cooldown::Int
    tighten::Float64
    tightenings_left::Int
    quiet_until::Int
    triggers::Int
    flat::Vector{Int}
    point::Vector{Float64}
    gain::Float64
end

CTSEMStallWatch(; window::Integer=80, fraction::Real=1e-2, cooldown::Integer=30,
    tighten::Real=0.1, tightenings::Integer=2) =
    CTSEMStallWatch(Int(window), Float64(fraction), Int(cooldown),
        Float64(tighten), Int(tightenings), 0, 0, Int[], Float64[], 0.0)

"""
    _ctsem_stall_verdict!(watch, trace, iteration, params, values, range, ratio)

Whether to stop: the progress test fired *and* something is flat.

Mutates `watch` with the hysteresis, and records which coordinates were flat so
the caller can say what ended the run. `params === nothing`, or a range with
nothing in it, means there is no transform layer to ask -- the conjunction can
then never complete, which is the right answer rather than half of one.
"""
function _ctsem_stall_verdict!(watch::CTSEMStallWatch, trace::CTSEMTrace,
        iteration::Integer, objective, params, values, range, ratio::Real;
        tolerance::Real=1e-6)
    watch.window >= 1 || return false
    params === nothing && return false
    iteration >= watch.quiet_until || return false
    _ctsem_stalled(trace, watch.window, watch.fraction) || return false
    watch.triggers += 1
    # Two detectors, because neither sees what the other does.
    # `_ctsem_flat_coordinates` measures a transform against its own live value
    # and so is scale free, but it walks `regular_transforms` and the
    # population scales and correlations are not in it. `_ctsem_saturated_for`
    # is the route's own, and on the laplace route it is the only thing that
    # knows a correlation has reached its cap -- which is exactly the state
    # that stalled the fit this was measured on. Missing it meant the
    # conjunction never fired on that fit at all.
    # Guarded separately, not as one expression. A detector that cannot answer
    # for this objective must not take the other one down with it: wrapping
    # both in a single `try` meant one `MethodError` reported nothing flat
    # anywhere, which is the failure mode that reads as "no problem found".
    relative = try
        _ctsem_flat_coordinates(params, values, range; ratio=ratio)
    catch
        Int[]
    end
    route = try
        _ctsem_saturated_for(objective, values)
    catch
        Int[]
    end
    flat = sort!(unique(vcat(relative, route)))
    if !isempty(flat)
        # And the third condition: somewhere better to go. Without it a fit
        # converging into a flat region is stopped before it arrives.
        #
        # The value is taken at `values` rather than from the trace's last row,
        # even though the trace has one and this costs an evaluation. They are
        # not the same point: the trace records completed iterations, and the
        # point handed in here is the last one the objective accepted, which
        # may be a line-search trial taken after that row was written. A gain
        # measured against the wrong baseline is not a gain.
        value = _ctsem_probe_value(objective, values)
        out = if isfinite(value)
            try
                _ctsem_overshot(objective, values, flat, value, tolerance)
            catch
                nothing
            end
        else
            nothing
        end
        if out !== nothing && out.overshot && !isempty(out.point)
            watch.flat = flat
            watch.point = collect(Float64, out.point)
            watch.gain = Float64(out.gain)
            return true
        end
    end
    # Stalled with nothing flat, or flat with nothing better nearby: slow
    # rather than stuck. Wait, and ask less readily next time.
    watch.quiet_until = Int(iteration) + watch.cooldown
    if watch.tightenings_left > 0
        watch.tightenings_left -= 1
        watch.fraction *= watch.tighten
    end
    return false
end

"""
    _ctsem_optimise_verdict(objective, minimizer, start_values, value, gradient_norm,
                            predicted_gain, last_gain, saturated_parameters,
                            g_tol, converge_tol; label, verbose)

Whether an optimiser arrived at a maximum, in one place for every route.

Optim's own `converged` is the disjunction of three criteria, and a line search
that fails on its first try satisfies the `f` one trivially: the objective did
not change because nothing was accepted. Measured on a 100-subject model, that
route stopped with `g_converged=false` and a largest gradient of 3.5e-07
against `g_tol=1e-8` and reported success. So convergence is judged on what is
still available, and on three separate ways of not being at a maximum.

`stalled` is an optimiser that never left its starting values while the
gradient there is not zero -- a fit that has not fitted anything, whatever its
flags say. `overshot` is an optimiser that stepped into the flat region of a
parameter's transform, where the gradient underflows to zero and every
tolerance passes: `_ctsem_overshot` separates that from a coordinate the data
simply do not identify by pulling the coordinate back and asking whether the
objective improves. Saturation on its own is a finding and not a failure, and
only the overshoot disqualifies the fit.

Optim's `g_tol` is an *absolute* gradient bound, and a log likelihood of order
1e3 puts 1e-8 out of reach however good the fit is -- L-BFGS runs out of line
search first and reports nothing converged -- so `Optim.g_converged` is kept
only as a sufficient condition, never a necessary one. `finite_gradient` is
tested explicitly because it can be true at a point whose gradient is NaN, set
on an earlier iterate, and `NaN <= tolerance` is false: one draw in ten reported
convergence with a NaN gradient. The failure that shaped all of this is still
the one to beat -- a 12-subject model stopping after a single iteration with a
gradient of 1.2e9 and a log likelihood of -7.7e6, reporting success -- and a
point like that has an enormous predicted gain, so it fails this criterion the
way it failed the last one.

## The criterion is not a gradient

`converge_tol` is in *nats*, and what is compared against it is
`predicted_gain` -- `1/2 g' B g` from the L-BFGS metric, read off the line
search by `CTSEMDirectional`. Both sides are objective differences, which is the
only comparison here that survives its own units.

A gradient does not. Under a reparameterisation `theta -> A theta` it becomes
`A^-T g`, so `|g| < c` says something different in every parameterisation, and
there is no divisor that repairs it: `max(1, |value|)`, the worst gradient the
run saw, and `sum(abs, subject_loglik)` were each tried here, and the best of
them is only the least arbitrary. An objective *difference*, by contrast, is a
log likelihood ratio: the likelihood's additive constants cancel, the
reparameterisation cancels, and 0.005 nats is negligible with six subjects or
six thousand while 10 nats matters with either. That is why `gaptol` is an
absolute number and correctly so, and this is the same number.

`1/2 g' B g` is the same quantity `.ctBackendOptimGap()` certifies with,
`1/2 g' H^-1 g`, under a weaker metric: `B` is L-BFGS's approximation and
depends on the retained secant pairs and the initial scaling, so it is affine
invariant only as far as that approximation is good. It can therefore stop a fit
and not certify one -- which is the division of labour already in place. When a
Hessian is computed, `.ctBackendCertifiedVerdict()` replaces this verdict with
the exact one; what is left here is the fits that certify nothing (`estonly`,
`certify = FALSE`).

`last_gain` -- the objective change over the final iteration -- is the second
sufficient condition, and it is here because `1/2 g'Bg` alone is not enough.
Measured on five fits against the exact gap:

    fit                    1/2 g'Bg    last step    exact gap
    saturated, converged      1.002            0     2.1e-15
    iteration cap 2           0.061        0.093       0.312
    iteration cap 8         1.5e-14            0     1.4e-15
    clean, converged        4.4e-15            0     3.4e-15
    clean, cap 3            3.0e-15      2.8e-14     2.3e-21

The metric fails twice. On the saturated fit it reports a whole nat where the
truth is 2e-15: a transform gone flat has a gradient underflowing to zero and a
`B` blowing up in the same direction, and `0 * Inf` is not a measurement.
`.ctBackendOptimGap()` handles that case by excluding flat directions from the
trusted subspace and *probing* them; nothing available here can. And on the
capped fit it reports a fifth of the truth, so it does not even err in the safe
direction.

What `last_gain` claims is weaker -- that no further progress was achievable,
not that nothing is left -- and it can be satisfied by a fit that stopped for a
reason other than being at a maximum. `stalled` and `overshot` are what stand
against the two ways that happens here. But it is in objective units and it
cannot be corrupted by the metric, which is where the other one breaks, so the
two together cover what neither does alone.

The overshoot bar is the same `converge_tol`, for the same reason: what
`_ctsem_overshot` measures is an objective gain, so it belongs against an
objective tolerance. One number, one unit, two uses it is dimensionally entitled
to.

`saturated_parameters` arrives as an argument because the routes detect it
differently: the marginal and joint ones ask `_ctsem_saturated_parameters` over
`_ctsem_saturation_range`, and the laplace one adds the levels whose
correlations saturate together (`_laplace_saturated_parameters`).
"""
function _ctsem_optimise_verdict(objective, minimizer, start_values, value,
        gradient_norm, predicted_gain, last_gain, saturated_parameters,
        g_tol, converge_tol;
        label::AbstractString="ctsem_optimize", verbose::Bool=false,
        overshoot_probe=_CTSEM_OVERSHOOT_PROBE[])
    moved = isempty(minimizer) ? 0.0 : maximum(abs, minimizer .- start_values)
    saturated = !isempty(saturated_parameters)
    stalled = moved == 0 && (!isfinite(value) || gradient_norm > max(g_tol, 1e-6))
    finite_gradient = isfinite(gradient_norm)
    # Both `Inf` before there is anything to report, so a fit that never took a
    # step cannot pass here on an uninitialised number.
    converged_enough = isfinite(value) && finite_gradient &&
        (predicted_gain <= converge_tol || last_gain <= converge_tol)
    overshoot = _ctsem_overshot(objective, minimizer, saturated_parameters,
        value, converge_tol; mode=overshoot_probe)
    overshot = overshoot.overshot
    verbose && stalled && println(_console(), label, ": the optimizer made no ",
        "progress from its starting values; reporting this as not converged")
    verbose && overshot && println(_console(), label, ": pulling raw ",
        "parameter(s) ", overshoot.coordinates, " back toward zero improves ",
        "the objective by ", overshoot.gain, ", so the estimate is not a ",
        "maximum; reporting this as not converged")
    verbose && saturated && !overshot && println(_console(), label,
        ": raw parameter(s) ", saturated_parameters, " have a materialising ",
        "transform that is flat to machine precision at the estimate, but no ",
        "pullback improves the objective, so this is a maximum with those ",
        "coordinates unidentified rather than a failed fit")
    (moved=moved, saturated=saturated, stalled=stalled,
        finite_gradient=finite_gradient,
        converged_enough=converged_enough, overshoot=overshoot,
        overshot=overshot)
end

"""Optimize a prepared likelihood entirely within Julia using L-BFGS."""
function ctsem_optimize(objective::CTSEMOptimisable, start::AbstractVector;
    maxiter::Integer=1000, g_tol::Real=1e-8, f_tol::Real=0.0,
    x_tol::Real=0.0, verbose::Bool=false, gradient_method=:adjoint,
    tune_chunks::Bool=true, lbfgs_memory::Integer=_CTSEM_LBFGS_MEMORY,
    progress_overwrite::Bool=true, progress_sink=nothing,
    progress_callback=nothing,
    progress::Bool=verbose, progress_label::AbstractString="optimise",
    progress_budget::Bool=false, progress_every::Real=0.0,
    gap_tol::Real=0.0, converge_tol::Real=1e-6,
    precondition=nothing, initial_alpha::Real=0.1,
    overshoot_probe=_CTSEM_OVERSHOOT_PROBE[],
    stall_window::Integer=80, stall_fraction::Real=1e-2,
    stall_cooldown::Integer=30, stall_tighten::Real=0.1,
    stall_tightenings::Integer=2, stall_ratio::Real=1e-3)
    start_values = collect(start)
    # Validated here rather than at the probe, which runs after the fit: a
    # misspelled mode should cost nothing, not a whole optimisation.
    overshoot_probe = _ctsem_overshoot_mode(overshoot_probe)
    invalid_objective = floatmax(eltype(start_values)) / 1e8
    gradient_limit = sqrt(floatmax(eltype(start_values)))
    label = _ctsem_optimise_label(objective)
    _ctsem_optimise_setup!(objective)
    # The route's own call accounting, or `nothing` where it keeps none. Passed
    # to every trial so the counting happens where the validity is decided.
    call_log = _ctsem_optimise_log(objective, verbose)
    # `evaluated`, not `result`: see `ctsem_laplace_optimize`. A closure's
    # assignment binds to the enclosing local of the same name, and the outer
    # Optim result below is called `result`.
    fg! = function (F, G, x)
        # The route decides what a usable trial point is: see
        # `_ctsem_optimise_trial`. On the laplace route a finite value at a
        # point whose inner Newton did not converge is *not* usable, because the
        # objective is only defined at the mode.
        trial = _ctsem_optimise_trial(objective, x, G !== nothing,
            gradient_method, gradient_limit, call_log)
        evaluated = trial.evaluated
        valid = trial.valid
        if !valid
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        if G !== nothing
            G .= -evaluated.gradient
        end
        # The point the callback will be asked about; see `current_x`.
        copyto!(current_x, x)
        return F === nothing ? nothing : -evaluated.value
    end
    # A callback rather than Optim's `show_trace`, which prints one dense line
    # per iteration whatever the model costs -- thousands on a fast one, and on
    # a slow one nothing for minutes. The objective and the gradient norm are
    # what say whether this is going anywhere: a log posterior that has stopped
    # moving while the gradient is still large is a fit in trouble, and that is
    # visible here long before the convergence flags are set.
    # Progress is not verbosity. Someone watching a fit wants to know it is
    # going somewhere; they do not thereby want the chunk-tuning timings and
    # the model-shape summary that `verbose` also turns on. Separating them is
    # what lets progress be the default without making the default noisy.
    reporter = CTSEMProgress(progress; label=progress_label,
        overwrite=progress_overwrite, every=progress_every, sink=progress_sink)
    # The trace records every iteration whatever `verbose` says: it costs a
    # push onto a vector, and a fit that turns out to have gone somewhere odd
    # is exactly the one nobody thought to turn reporting on for.
    # `predicted_gain` sits between the two fixed columns and the route's own,
    # so `_ctsem_optimise_trace_keys` keeps naming only what a route adds.
    trace = CTSEMTrace(:objective, :gradient_norm, :predicted_gain,
        _ctsem_optimise_trace_keys(objective)[3:end]...)
    watcher = CTSEMCallback(progress_callback)
    # How far this run has come toward the nearest of its stopping rules; see
    # `CTSEMConvergence`. Fed every iteration rather than every printed line,
    # because the scales it interpolates on are the worst gradient and the
    # largest objective change the fit ever had, and the printed lines are a
    # time-sampled subset.
    convergence = CTSEMConvergence(g_tol, Int(maxiter))
    # An observer around the line search that runs; it changes nothing about
    # which one that is.
    directional = CTSEMDirectional(Optim.LineSearches.BackTracking())
    # Optim's own `iterations`, `f_calls` and `g_calls` do not survive a
    # callback stop: measured on a two-latent model, four runs that stopped
    # after 4, 8 and 10 iterations all reported 1 iteration and 2 gradient
    # calls, while the trace -- one row per iteration, recorded here -- showed
    # 5, 9 and 11 rows and log likelihoods five orders apart. So the iteration
    # count is taken from what this callback saw, which is the same number on a
    # run that ends any other way.
    seen_iterations = Ref(0)
    stopped_by_gap = Ref(false)
    stopped_by_stall = Ref(false)
    # Where the fit currently is. Optim's callback is handed convergence
    # numbers, not the point they describe -- `extended_trace` would carry it
    # but costs a copy of every iterate -- and the stall conjunction has to ask
    # the transforms what they are doing *at this point*. `fg!` sees every
    # trial, so the last one it accepted is recorded there and read here.
    current_x = copy(start_values)
    stall = CTSEMStallWatch(window=stall_window, fraction=stall_fraction,
        cooldown=stall_cooldown, tighten=stall_tighten,
        tightenings=stall_tightenings)
    watch = function (state)
        latest = state isa AbstractVector ? last(state) : state
        _record!(trace, latest.iteration, -latest.value, latest.g_norm,
            _ctsem_predicted_gain(directional),
            _ctsem_optimise_trace_values(objective)...)
        seen_iterations[] = max(seen_iterations[], Int(latest.iteration))
        percent = _convergence_percent!(convergence, latest.g_norm,
            latest.value, Int(latest.iteration))
        if _due(reporter)
            # Never on a budget stage: its own fraction is exact, and an
            # estimate would replace a correct denominator with a guess.
            _progress_optimise(reporter, latest.iteration, Int(maxiter),
                @sprintf("logpost %11.2f", -latest.value),
                @sprintf("|g| %9.2e", latest.g_norm),
                _ctsem_optimise_progress_extra(objective)...;
                budget=progress_budget,
                percent=progress_budget ? NaN : percent)
        end
        # Its own cadence, so passing a callback with `verbose = 0` -- the
        # obvious combination for a front end that draws rather than prints --
        # still reports.
        # The point comes with the numbers. Without it a caller can watch a
        # multi-hour fit and still have nothing to restart from when it is
        # interrupted, because nothing is written until the fit returns --
        # which has cost a run here. `current_x` is the last accepted trial,
        # the same point the reported objective and gradient describe.
        _invoke_callback(watcher, latest.iteration, Int(maxiter),
            -latest.value, latest.g_norm, current_x)
        # Stop when the step just taken was predicted to gain less objective
        # than asked for. The predicted gain of the *next* step is not knowable
        # here, and the last one is the standard stand-in: a quasi-Newton
        # direction that has stopped promising anything is not about to start.
        #
        # A proxy, and treated as one -- `B` is limited memory and carries its
        # own scaling, so this is not the invariant decrement. The exact check
        # after the fit is what certifies, and what resumes with a tightened
        # rule when this stopped too early.
        if gap_tol > 0 && _ctsem_predicted_gain(directional) < gap_tol
            stopped_by_gap[] = true
            return true
        end
        # And stop when the fit has stopped getting anywhere *and* a transform
        # has gone flat, which is the structural reason it stopped. Neither
        # half is a stopping rule on its own -- see `_ctsem_stall_verdict!` for
        # why, and for the cooldown that keeps a merely slow fit from being
        # asked over and over.
        if _ctsem_stall_verdict!(stall, trace, latest.iteration, objective,
                _ctsem_params(objective), current_x,
                _ctsem_saturation_range(objective, current_x), stall_ratio;
                tolerance=converge_tol)
            stopped_by_stall[] = true
            return true
        end
        return false
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol,
        f_reltol=f_tol, x_abstol=x_tol, show_trace=false, store_trace=false,
        callback=watch, extended_trace=false)
    # See `ctsem_tune_chunks!`: the subject loop is not monotone in the chunk
    # count, so the count is measured on this model rather than taken from
    # `cores`.
    tuning = tune_chunks ? ctsem_tune_chunks!(
        () -> ctsem_evaluate(objective, start_values; gradient=true,
            gradient_method=gradient_method); verbose=verbose) : nothing
    # After the tuning, which is what decided the chunk count it reports.
    verbose && _ctsem_optimise_verbose_shape(objective)
    # A first step of unit *length*, not of unit alpha.
    #
    # L-BFGS has no curvature history on its first iteration, so it takes the
    # steepest-descent direction with whatever the initial step guess gives.
    # Optim's default is `InitialStatic()`, an unscaled alpha of one -- which
    # means the first step is as long as the gradient, and ctsem's gradients
    # are routinely of magnitude tens. Watched on a binary model: one step from
    # raw 0 to raw 20.9. That step *improved* the objective, so the line search
    # was right to take it, but it overshot the optimum at raw ~1 and landed
    # where every ctsem transform is flat to machine precision, and a zero
    # gradient ends the optimisation.
    #
    # `scaled=true` divides alpha by the gradient norm, so the first step has
    # length one in parameter space regardless of how steep the objective is.
    # That is the standard remedy and it is what the saturation guard below
    # would otherwise spend its life reporting.
    #
    # Backtracking, not Hager-Zhang. A Wolfe line search needs the directional
    # derivative at every trial point, so every trial costs a full adjoint;
    # Armijo needs only the objective, and on this engine that is genuinely
    # cheaper -- `ctsem_evaluate(..., gradient=false)` takes the forward-only
    # branch, with no reverse pass and no Frechet blocks, those existing solely
    # for derivatives. Armijo imposes no curvature condition and so stops
    # polishing sooner, which was the reason it was not the default: measured
    # here, 47 iterations against 36, finishing near 1e-5 where Hager-Zhang
    # reaches 1e-9. Both halves of that objection have gone. It counted
    # iterations rather than their cost, and the precision it was protecting
    # was protecting a *gradient* convergence test -- where the estimate stands
    # relative to the optimum is now measured exactly after the fit, in
    # objective units, and the fit continues when it matters.
    linesearch = "backtracking"
    result = Optim.optimize(Optim.only_fg!(fg!), start_values,
        Optim.LBFGS(m=Int(lbfgs_memory),
            alphaguess=Optim.LineSearches.InitialStatic(
                alpha=Float64(initial_alpha), scaled=true),
            linesearch=directional,
            P=_ctsem_metric(precondition, length(start_values))), options)
    # No rescue stage. The one that stood here existed because Hager-Zhang can
    # run out of line search and return the iterate it had reached while Optim
    # reports a finished optimisation -- measured on a binary model: two
    # iterations, 68 objective evaluations, a final gradient of 475, and the
    # last step uphill. It answered by restarting with backtracking, which is
    # now what runs in the first place, and over the eight fits that reproduced
    # that stall backtracking converged every one.
    #
    # A stop that is short for any other reason is not the line search's
    # problem to solve twice: the certification measures what the estimate
    # still has to gain, and `.ctBackendCorrectResult()` continues the fit from
    # a damped Newton step with a tightened stopping rule. Two mechanisms for
    # one job, where the second can only act in cases the first did not fix, is
    # a way to be surprised rather than a safety net.
    # Before the final evaluation, so what it reports is the run rather than
    # the extra call: see `_ctsem_optimise_verbose_report`.
    verbose && _ctsem_optimise_verbose_report(objective, call_log)
    minimizer = collect(Optim.minimizer(result))
    final = ctsem_evaluate(objective, minimizer; gradient=true,
        contributions=true, gradient_method=gradient_method)

    # Optim's own `converged` is the disjunction of three criteria, and a line
    # search that fails on its first try satisfies the `f` one trivially: the
    # objective did not change because nothing was accepted. Measured on a
    # 100-subject model, this route stopped with `g_converged=false` and a
    # largest gradient of 3.5e-7 against `g_tol=1e-8` and reported success.
    #
    gradient_norm = isempty(final.gradient) ? 0.0 : maximum(abs, final.gradient)
    # `max(shown, ...)`: when the backtracking fallback ran and was kept,
    # `Optim.iterations` describes that second run alone, which can be fewer
    # than the user already watched go past. The closing line closes what was
    # on screen.
    iterations = max(reporter.shown, Optim.iterations(result))
    # Running out of iterations is a different outcome from converging, and the
    # closing line used to report both as a bare count. On a stage whose cap is
    # the plan (`budget`) reaching it is not news; anywhere else it is the one
    # thing about the fit the user most needs to know.
    capped = !progress_budget && iterations >= Int(maxiter)
    progress && _progress_done(reporter,
        @sprintf("%d iterations%s", iterations,
            capped ? " -- ITERATION CAP REACHED, not converged" : ""),
        @sprintf("logpost %.4f", final.value),
        @sprintf("|g| %.2e", gradient_norm))
    # Forced, whatever the cadence says: a rate-limited callback on a fit that
    # finishes inside one interval would otherwise never fire at all, and the
    # final state is the one a live plot most needs.
    _invoke_callback(watcher, Optim.iterations(result), Int(maxiter),
        final.value, gradient_norm, minimizer; force=true)
    # A saturated transform reports a zero gradient, and a zero gradient is
    # indistinguishable from an optimum.
    #
    # This used to be judged on the raw coordinate's *magnitude* -- flagged
    # once |raw| >= 20, on the reasoning that every transform ctsem writes is
    # flat to machine precision by then. That reasoning is only half right: it
    # is a statement about how far a transform's *derivative* has collapsed,
    # and raw magnitude is a proxy for that which fails in both directions. An
    # identity transform (`MANIFESTMEANS`, `T0MEANS` with no scale) has
    # derivative exactly 1 at any raw magnitude and never saturates, so a
    # model fit to data with a mean of 25 or 100 -- an unremarkable fit --
    # reported not converged purely because its mean landed past 20 on the raw
    # scale it is never transformed away from. And different transforms go
    # flat at different raw magnitudes in the first place: a drift diagonal's
    # `-(1e-6 + 2log1p_exp(-2x))` is already down to derivative 8e-9 by raw 10,
    # while a correlation's `2/(1+exp(-x))-1` is still at 9e-5 there. No single
    # cutoff on the raw value fits both, and precomputing one per parameter
    # (per transform, per state-dependent case) is the tedium the derivative
    # check below avoids entirely.
    #
    # So saturation is judged directly on the transform's own derivative at
    # the estimate -- see `_ctsem_saturated_parameters`
    # (`parameter_transforms.jl`) -- which is the one quantity that actually
    # says whether a raw coordinate still does anything: an optimiser that
    # oversteps into the flat region finds `gradient_norm` underflowing to
    # zero there, satisfies any tolerance, and reports success. Observed on a
    # binary model: one L-BFGS iteration to raw 20.9, declared converged, log
    # likelihood -730.7 where the profile peak is -714.5.
    #
    # `stalled` does not catch it, because the optimiser did move -- it moved
    # too far. So saturation is its own verdict: at a flagged coordinate the
    # parameter is not identified by the data but by the transform's
    # floating-point limit, and calling that converged is the wrong answer
    # confidently delivered.
    # Over the *transformed* coordinates only. Saturation is a statement about
    # ctsem's parameter transforms going flat, and the joint target's vector
    # also carries state innovations, which have no transform and no flat
    # region -- a trajectory five standard deviations out is unusual data, not
    # an unidentified parameter, and reading it as saturation would report
    # every such fit as failed.
    saturated_parameters = _ctsem_saturated_for(objective, minimizer)
    # Before the verdict, because the verdict probes. A route's own result
    # fields can be live per-run state -- the laplace route's `inner_converged`,
    # `inner_iterations` and the rest are -- and every probe point overwrites
    # them, so read after the probe they describe wherever it went last rather
    # than the estimate.
    result_extra = _ctsem_optimise_result_extra(objective, final, call_log)
    # And the verdict itself, which every route reaches the same way and in one
    # place: see `_ctsem_optimise_verdict`.
    verdict = _ctsem_optimise_verdict(objective, minimizer, start_values,
        final.value, gradient_norm, _ctsem_predicted_gain(directional),
        _ctsem_last_gain(trace), saturated_parameters, g_tol, converge_tol;
        label=label, verbose=verbose, overshoot_probe=overshoot_probe)
    saturated = verdict.saturated
    stalled = verdict.stalled
    finite_gradient = verdict.finite_gradient
    converged_enough = verdict.converged_enough
    overshoot = verdict.overshoot
    overshot = verdict.overshot
    verbose && stopped_by_stall[] && println(_console(), label, ": the last ",
        stall.window, " iterations gained ",
        _ctsem_window_gain(trace, stall.window), " against ",
        _ctsem_trace_progress(trace), " gained over the run, and raw ",
        "parameter(s) ", stall.flat, " have lost all but ", stall_ratio,
        " of what their transform does when live, and pulling them back gains ",
        stall.gain, " -- so the fit has stopped getting anywhere, this is why, ",
        "and there is somewhere better to go. Stopping here rather than ",
        "running to the iteration cap")
    verbose && !stalled && !(finite_gradient &&
        (Optim.g_converged(result) || converged_enough)) &&
        println(_console(), label, ": the optimizer stopped with an estimated ",
            _ctsem_predicted_gain(directional), " log likelihood still ",
            "available and ", _ctsem_last_gain(trace),
            " gained on its last iteration, against a tolerance of ",
            converge_tol, " (largest gradient ", gradient_norm,
            "); reporting this as not converged")
    return (
        minimizer=minimizer,
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        result_extra...,
        # The larger of the two: they agree unless the callback stopped the
        # run, in which case Optim's is the one that stopped being updated.
        iterations=max(Optim.iterations(result), seen_iterations[]),
        # Left as Optim reports them, and undercounted for the same reason when
        # the run was stopped by the callback -- there is no second source for
        # these, and inventing one would be worse than a number whose limit is
        # written down. `stopped_by_gap` is what says the run is such a case.
        f_calls=Optim.f_calls(result),
        g_calls=Optim.g_calls(result),
        stopped_by_gap=stopped_by_gap[],
        # Whether the run ended because it stopped making progress rather than
        # because it arrived. `f_calls` and `g_calls` undercount here for the
        # same reason they do under `stopped_by_gap`: Optim stops updating them
        # when a callback ends the run.
        stopped_by_stall=stopped_by_stall[],
        stall_window=Int(stall.window),
        # Which coordinates were flat when it stopped, and how many times the
        # progress test fired without finding one. `0` for none, for the reason
        # `saturated_parameters` gives.
        stall_parameters=isempty(stall.flat) ? [0] : stall.flat,
        stall_triggers=stall.triggers,
        # The point the in-flight probe found, so the caller resuming from it
        # does not pay for the same ladder twice. Empty unless it stopped here.
        stall_point=isempty(stall.point) ? Float64[0.0] : stall.point,
        stall_gain=stall.gain,
        stalled=stalled,
        chunks=ctsem_max_chunks().max_chunks,
        # Reported because the answer depends on it. Dividing a unit's members
        # reassociates a sum, so two fits at different widths differ in the last
        # few digits -- deterministic at a fixed width, and the thing to hold
        # constant for a before-and-after.
        chunk_timings=tuning === nothing ? Tuple{Int,Int,Float64}[] : tuning.timings,
        gradient_norm=gradient_norm,
        converge_tol=Float64(converge_tol),
        last_gain=_ctsem_last_gain(trace),
        # What the last step was predicted to gain, and the rule it was judged
        # against. `Inf` when no line search ran, and `0` when the rule was off.
        predicted_gain=_ctsem_predicted_gain(directional),
        gap_tol=Float64(gap_tol),
        linesearch=linesearch,
        # See `ctsem_laplace_optimize`: `Optim.converged` includes the x and
        # f criteria, which a line search that stops making progress satisfies
        # trivially, so convergence is judged on the gradient alone.
        #
        # `overshot`, not `saturated`. `converged` answers one question -- did
        # the optimizer arrive at a maximum -- and a coordinate the data do not
        # identify is a separate finding, reported separately in `saturated`
        # and `saturated_parameters`. Keying convergence on saturation made the
        # flag false on two thirds of good fits; see `_ctsem_overshot`.
        converged=!stalled && !overshot && finite_gradient &&
            (Optim.g_converged(result) || converged_enough),
        saturated=saturated,
        # The pullback verdict and its margin. `overshot` is the half of
        # saturation that is a convergence failure; `overshoot_gain` is how
        # much the objective improved when the flagged coordinate was pulled
        # back, so a user can see whether it was 16 log units or 1e-13.
        overshot=overshot,
        overshoot_gain=overshoot.gain,
        # Which coordinates the pullback moved. Not the same as
        # `saturated_parameters` any more and no longer derivable from it: the
        # probe selects by magnitude order, so the set that proved the estimate
        # is not a maximum can include coordinates whose transform is merely
        # unresponsive rather than flat. `0` means none, for the reason the
        # saturated list gives.
        overshoot_parameters=isempty(overshoot.coordinates) ? [0] :
            overshoot.coordinates,
        # And where it went to prove it. The probe has already paid for this
        # point -- it is the trial that beat the estimate -- so handing it back
        # lets a caller resume from somewhere measurably better for no further
        # evaluations, rather than only being told the estimate is not a
        # maximum. Used by `.ctBackendStallEscape()`.
        #
        # `[0.0]` rather than `Float64[]` for the reason `saturated_parameters`
        # gives: a zero-length vector deadlocks the R bridge. A real point is
        # `length(minimizer)` long, which is how the caller tells them apart.
        overshoot_point=isempty(overshoot.point) ? Float64[0.0] : overshoot.point,
        # Which raw parameters, not just whether one did -- most of the
        # diagnostic value, and free once the derivatives are computed.
        #
        # NEVER an empty vector. A zero-length vector deadlocks the
        # JuliaConnectoR bridge in both directions, and "nothing saturated" is
        # the normal case, so returning Int[] here hung every healthy fit the
        # moment its result crossed back to R. 0 means none; any other entry is
        # a raw parameter index.
        saturated_parameters=isempty(saturated_parameters) ? [0] : saturated_parameters,
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
        trace=_trace_result(trace),
    )
end
