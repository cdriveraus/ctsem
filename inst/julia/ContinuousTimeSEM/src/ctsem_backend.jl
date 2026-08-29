"""
    CTSEMObjective(params, subject_starts, timesteps, data)

Prepared continuous-time likelihood for the R `ctsem` backend. `data` is
stored as manifest variables by observation. Subject starts are prepared by R,
which owns data ordering and predictor preprocessing. This type retains only
the numerical per-subject EKF workspaces needed for repeated evaluation.
"""
mutable struct CTSEMObjective{P,O}
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
end

CTSEMObjective(params, subject_objectives) =
    CTSEMObjective(params, subject_objectives, nothing, Int[], Float64[], 1.0)

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

export CTSEMObjective, ctsem_objective, ctsem_evaluate, ctsem_optimize

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
                max_timestep=subject_objective.max_timestep)
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
    max_timestep::Real=Inf; prior_index=Int[], prior_scale=Float64[],
    prior_weight::Real=1.0)
    ranges = _ctsem_subject_ranges(subject_starts, timesteps, data)
    size(tdpred_data, 2) == size(data, 2) || throw(DimensionMismatch("TD predictor columns must match observations"))
    size(tipred_data, 1) == length(ranges) || throw(DimensionMismatch("TI predictor rows must match subjects"))
    # Copy each subject once. This avoids R proxy/view lifetime issues and makes
    # the objective safe to retain in a Julia session.
    objects = Any[
        ContinuousEKFObjective(params, Matrix(view(data, :, r)), collect(view(timesteps, r));
            tdpreds=Matrix(view(tdpred_data, :, r)), tipreds=vec(tipred_data[i, :]), subject=i,
            max_timestep=max_timestep)
        for (i, r) in enumerate(ranges)
    ]
    length(prior_index) == length(prior_scale) ||
        throw(DimensionMismatch("prior index and scale vectors must have equal length"))
    return CTSEMObjective(params, objects, nothing, Int.(prior_index),
        Float64.(prior_scale), Float64(prior_weight))
end

ctsem_objective(params::EKFParameters, subject_starts, timesteps, data,
    tdpred_data=zeros(eltype(data), 0, size(data, 2)),
    tipred_data=zeros(eltype(data), length(subject_starts), 0), max_timestep::Real=Inf;
    prior_index=Int[], prior_scale=Float64[], prior_weight::Real=1.0) =
    CTSEMObjective(params, subject_starts, timesteps, data, tdpred_data, tipred_data,
        max_timestep; prior_index=prior_index, prior_scale=prior_scale,
        prior_weight=prior_weight)

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
not free and it is not even monotone: the reverse sweep allocates on the order
of ten thousand small arrays per subject, and above roughly ten million
allocations a second the allocator, not the arithmetic, is what the threads are
queueing for. Measured on a 23-core machine, one 24-row subject of a
one-latent, one-indicator model runs **3.7x slower on 23 threads than on one**,
while the same code on a twenty-latent model runs 2.7x *faster*. Both are the
same subject loop; what differs is how much arithmetic sits between two
allocations.

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
        repeats = serial < 0.05 ? 3 : 1
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
            evaluate()                       # warm this chunk count's workspaces
            elapsed = Inf
            for _ in 1:repeats
                elapsed = min(elapsed, @elapsed evaluate())
            end
            push!(timings, (n, elapsed))
            verbose && println("Chunk tuning: ", n, " chunk(s) ",
                round(elapsed; digits=4), " s")
            if elapsed < best_time * 0.95
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
    verbose && println("Chunk tuning: using ", best_chunks, " chunk(s) of at most ", limit)
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
        return total + _ctsem_log_prior(objective, values)
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
    return total + _ctsem_log_prior(objective, values)
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
reaches that through `backendcontrol`.
"""
const _CTSEM_LBFGS_MEMORY = 20

"""Optimize a prepared likelihood entirely within Julia using L-BFGS."""
function ctsem_optimize(objective::CTSEMObjective, start::AbstractVector;
    maxiter::Integer=1000, g_tol::Real=1e-8, f_tol::Real=0.0,
    x_tol::Real=0.0, verbose::Bool=false, gradient_method=:adjoint,
    tune_chunks::Bool=true, lbfgs_memory::Integer=_CTSEM_LBFGS_MEMORY,
    progress_overwrite::Bool=true, progress_callback=nothing,
    progress::Bool=verbose)
    start_values = collect(start)
    invalid_objective = floatmax(eltype(start_values)) / 1e8
    gradient_limit = sqrt(floatmax(eltype(start_values)))
    fg! = function (F, G, x)
        result = try
            ctsem_evaluate(objective, x; gradient=G !== nothing,
                gradient_method=gradient_method)
        catch
            nothing
        end
        valid = result !== nothing && isfinite(result.value)
        if valid && G !== nothing
            valid = all(isfinite, result.gradient) && all(abs(value) < gradient_limit for value in result.gradient)
        end
        if !valid
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        if G !== nothing
            G .= -result.gradient
        end
        return F === nothing ? nothing : -result.value
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
    reporter = CTSEMProgress(progress; label="optimise",
        overwrite=progress_overwrite)
    # The trace records every iteration whatever `verbose` says: it costs a
    # push onto a vector, and a fit that turns out to have gone somewhere odd
    # is exactly the one nobody thought to turn reporting on for.
    trace = CTSEMTrace(:objective, :gradient_norm)
    watcher = CTSEMCallback(progress_callback)
    watch = function (state)
        latest = state isa AbstractVector ? last(state) : state
        _record!(trace, latest.iteration, -latest.value, latest.g_norm)
        if _due(reporter)
            _progress_line(reporter, latest.iteration, Int(maxiter),
                @sprintf("logpost %11.2f", -latest.value),
                @sprintf("|g| %9.2e", latest.g_norm))
        end
        # Its own cadence, so passing a callback with `verbose = 0` -- the
        # obvious combination for a front end that draws rather than prints --
        # still reports.
        _invoke_callback(watcher, latest.iteration, Int(maxiter),
            -latest.value, latest.g_norm)
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
    result = Optim.optimize(Optim.only_fg!(fg!), start_values,
        Optim.LBFGS(m=Int(lbfgs_memory),
            alphaguess=Optim.LineSearches.InitialStatic(scaled=true)), options)
    minimizer = collect(Optim.minimizer(result))
    final = ctsem_evaluate(objective, minimizer; gradient=true,
        contributions=true, gradient_method=gradient_method)

    # The same convergence verdict `ctsem_laplace_optimize` reaches, and for the
    # same reasons -- this route was simply left behind when that one was fixed,
    # which is worse than it sounds because this is ctsem's *default* route.
    #
    # Optim's own `converged` is the disjunction of three criteria, and a line
    # search that fails on its first try satisfies the `f` one trivially: the
    # objective did not change because nothing was accepted. Measured on a
    # 100-subject model, this route stopped with `g_converged=false` and a
    # largest gradient of 3.5e-7 against `g_tol=1e-8` and reported success.
    #
    # `g_tol` is an *absolute* bound, so a log likelihood of order 1e3 puts 1e-8
    # out of reach however good the fit is. `scaled_tolerance` is the criterion
    # that scales with the problem; it is an addition to the strict test, never
    # a loosening of it, and `stalled` is what stops a fit that never moved from
    # passing either.
    moved = isempty(minimizer) ? 0.0 : maximum(abs, minimizer .- start_values)
    gradient_norm = isempty(final.gradient) ? 0.0 : maximum(abs, final.gradient)
    progress && _progress_done(reporter,
        @sprintf("%d iterations", Optim.iterations(result)),
        @sprintf("logpost %.4f", final.value),
        @sprintf("|g| %.2e", gradient_norm))
    # Forced, whatever the cadence says: a rate-limited callback on a fit that
    # finishes inside one interval would otherwise never fire at all, and the
    # final state is the one a live plot most needs.
    _invoke_callback(watcher, Optim.iterations(result), Int(maxiter),
        final.value, gradient_norm; force=true)
    # A saturated transform reports a zero gradient, and a zero gradient is
    # indistinguishable from an optimum.
    #
    # Every transform ctsem writes -- `log1p_exp(2x)`, `-(1e-6 + 2log1p_exp(-2x))`,
    # `2/(1+exp(-x))-1` -- is flat to machine precision by |x| ~ 20: the
    # exponential underflows and the derivative is *exactly* zero, not merely
    # small. An optimiser that oversteps into that region then finds
    # `gradient_norm = 4e-16`, satisfies any tolerance, and reports success.
    # Observed on a binary model: one L-BFGS iteration to raw 20.9, declared
    # converged, log likelihood -730.7 where the profile peak is -714.5.
    #
    # `stalled` does not catch it, because the optimiser did move -- it moved
    # too far. So saturation is its own verdict: past this point the parameter
    # is not identified by the data but by the transform's floating-point
    # limit, and calling that converged is the wrong answer confidently
    # delivered.
    saturated = isempty(minimizer) ? false : maximum(abs, minimizer) >= _CTSEM_SATURATION[]
    stalled = moved == 0 && (!isfinite(final.value) || gradient_norm > max(g_tol, 1e-6))
    scaled_tolerance = max(g_tol, 1e-6 * max(one(gradient_norm), abs(final.value)))
    converged_enough = isfinite(final.value) && gradient_norm <= scaled_tolerance
    verbose && stalled && println("ctsem_optimize: the optimizer made no progress ",
        "from its starting values; reporting this as not converged")
    verbose && saturated && println("ctsem_optimize: a parameter reached ",
        maximum(abs, minimizer), " on the unconstrained scale, where its ",
        "transform is flat to machine precision; reporting this as not converged")

    return (
        minimizer=minimizer,
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        row_loglik=final.row_loglik,
        iterations=Optim.iterations(result),
        f_calls=Optim.f_calls(result),
        g_calls=Optim.g_calls(result),
        stalled=stalled,
        chunks=ctsem_max_chunks().max_chunks,
        chunk_timings=tuning === nothing ? Tuple{Int,Float64}[] : tuning.timings,
        gradient_norm=gradient_norm,
        scaled_tolerance=scaled_tolerance,
        converged=!stalled && !saturated && (Optim.converged(result) || converged_enough),
        saturated=saturated,
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
        trace=_trace_result(trace),
    )
end
