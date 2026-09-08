"""
Running the sampler: chains, warmup, and what comes back.

# Two parallel axes, and which gets the threads

The gradient parallelises over units -- the log likelihood is a sum over them,
so each chunk accumulates a partial value and partial gradient and they are
summed. Chains parallelise too, and share nothing at all.

Threads go to chains first, and the unit loop takes what is left.

The unit axis, measured on a fit rather than on a warm re-evaluation: a
5-latent, 200-subject, 12-occasion, 110-parameter model, `inits` fixed so every
run starts from the same point, `estonly` so the uncertainty phase is not
included, iterations capped at 15.

    cores      1        2        4        8       16
    fit      160.7 s  63.4 s   43.0 s   32.2 s   29.4 s
    speedup    --      1.72x    2.53x    3.38x    3.70x

**Read that as a shape, not as a scaling curve.** It came off an i9-12900KS,
8 P-cores with SMT plus 8 E-cores without, so "8 cores" is eight heterogeneous
threads and the sixteenth is nothing like the first. On 23 homogeneous cores the
same axis reaches 5.28x at eight chunks. See the end of this docstring; every
timing here taken on that machine has the same flaw.

Three repeats of `cores = 1` agreed to 7%, and the sweep was run *descending*
so that a speedup could not be confused with having run later in the session.
Both matter, and neither is a formality:

  - Run ascending on a contended machine, three identical `cores = 1` fits gave
    271 s, 301 s and 125 s -- a 140% spread, wider than any effect being looked
    for. Every fit-level number taken while another job shared the machine is
    worthless, and several were.
  - Times also fall with position in the session even when the work is
    identical, so an ascending sweep confounds "more cores" with "ran later".
    Ascending gave 4.18x at sixteen where descending gives 3.70x; the direction
    survives the reversal, which is what makes it real.

`logposterior` came back as -11644.64 from every one of those fits, at every
chunk count. The split is over a sum, so it should be exact, and it is.

**Chains do not gain linearly, and the reason for serving them first is not the
one that used to be written here.** Measured on the same model, one converged
fit reused so that every run starts from an identical estimate and Hessian,
`maxdepth = 3` so both configurations do equal work, order balanced A B B A A B
B A, and each cell repeated: two chains at once on T/2 chunks against the same
two chains in turn on T chunks.

    threads    8       16
    gain     1.12x   1.20x

against a baseline spread of 0.9-4.3%. Nowhere near the 2x that "linearly"
implies.

The ordering survives anyway, because the *other* axis has already stopped
paying. Per-gradient cost from the same runs: one chain on 8 chunks 0.0122 s,
one chain on 16 chunks 0.0134 s -- **10% slower on twice the threads**. So
threads nine to sixteen are worth 1.20x given to a second chain and 0.91x given
to the subject loop. That is the argument for the policy; the crossover sits
where the unit loop saturates, between 8 and 16 chunks for this 200-subject
target, and it will sit elsewhere for a different one.

Two things temper it. The advantage is contingent on the chains being balanced,
and **it inverts when they are not**: on the chain pair NUTS produced at the
first seed tried, one chain adapted to a degenerate step size of 3.7e-9 while
the other saturated `maxdepth`, a 2.08:1 imbalance, and running them at once was
12% *slower* at 8 threads because the light chain finished early and left half
the threads idle. One chain in 13 probed failed to adapt a usable step size, so
roughly 15% of two-chain runs and 27% of four-chain runs are lopsided like that,
and milder imbalance is universal. And capping `maxdepth` at 3 to equalise the
work also removes most of the natural per-draw variation, so at the default of
10 the imbalance -- and the loss to it -- is larger than measured here.

**Two claims that stood here have been withdrawn, and why matters to whoever
measures next.**

This recorded eight processes on disjoint subject groups reaching 2.60x against
3.38x for eight threads, concluding the ceiling was not the shared allocator.
Both figures came from a development machine that is an i9-12900KS: 8 P-cores
with SMT plus 8 E-cores without. "Eight threads" there is eight *heterogeneous*
threads, an E-core being roughly half a P-core, so no curve measured on it can
be read as N times one core. Every local timing in this file's history carries
that flaw.

Re-measured on 23 homogeneous cores, eight independent processes reach **8.19x**
throughput: 0.14888 s per gradient alone, 0.14157-0.14548 s each with eight
running. They do not interfere at all, though each allocates 544 MB per
gradient.

Nor is the subject loop bandwidth bound, which the 2.60x had been taken to
suggest. Against two references measured in the same session -- a dependent FMA
chain with no memory traffic, and a deliberately bandwidth-bound array sum --
the real work lands near the arithmetic ceiling and beats the bandwidth one:

    chunks    gradient   scores   Hessian   arithmetic   bandwidth
    8           5.28x     5.00x    4.66x      7.76x        3.19x

A bandwidth-bound workload cannot outrun the bandwidth-bound reference. By
volume the loop uses about 3% of what the machine has, 3.3 GB/s against 110
GB/s. What remains is garbage collection, 22-24% of gradient wall time.

So threads versus processes on the unit axis is **open**, not settled against
processes. What would settle it is the scatter/gather barrier a process split
needs on every gradient -- not a scaling ceiling that turned out to be the
measuring machine.

Nothing needs deciding by the caller either way: `workspace_slot` on
`ctsem_sample_density!` is present exactly when chains are running concurrently,
and a chain given a block of slots chunks its units inside it.

Splitting units across *processes* rather than threads is a different question
and the answer differs by phase. Sampling does on the order of `draws x
(2^depth - 1)` gradients, each needing a scatter of `theta` and a gather of the
gradient, so a per-gradient barrier is paid ~10^5 times; chains, which need no
barrier at all until the end, are the better process axis there. Optimisation
does hundreds to low thousands of gradients at seconds apiece, where a round
trip carrying `2 x npar` doubles is lost in the noise -- so if the unit axis
ever does hit a ceiling on threads, that is the phase where moving it to
processes would pay, and it would cut per-worker memory rather than multiply it.

# What crosses the bridge

`save_effects` defaults to false, and that is a bandwidth decision rather than a
statistical one. A hundred subjects with two effects each and four chains of a
thousand draws is 800,000 numbers -- 6.4 MB, and the bridge is slow enough
inbound that returning them costs longer than many fits do. The per-effect
posterior mean and standard deviation come back always, because they are `ndim`
numbers rather than `ndim * ndraws` and answer most of what the draws would be
used for.

The bridge, measured directly rather than inferred from one transfer:

  - **~0.8 ms per message on Linux once the socket is tuned**, and ~41 or ~82 ms
    before it, depending on whether one direction split its writes or both did.
    Windows was always ~0.3 ms.
  - **outbound ~870 MB/s** -- sending 32 MB costs 37 ms, so arguments are
    effectively free.
  - **inbound ~2.0 MB/s** -- the only per-byte term that matters.

The fixed cost used to be the part that surprised, and it was not a latency at
all: JuliaConnectoR assembles a message from many small writes, and on Linux
that write-write-read pattern met Nagle on the sender and a delayed ACK on the
receiver, stalling ~40 ms per message whatever it carried. `ctJuliaSetup` now
turns both halves off; see `bridge_tuning.jl`. **The figures that motivated
call-count batching were taken before that fix**, so read them as history: one
gradient cost 0.417 s from R against 0.005 s in the engine, and the post-fit
uncertainty phase measured ~99% bridge on Linux -- of which almost all was the
stall. Batching narrow calls cut that phase 1.48-1.90x and is still worth having,
since fewer messages is still fewer messages, but it is now worth milliseconds
rather than seconds, and nothing here justifies contorting an interface to save
a call.
"""

using LinearAlgebra
using Random
using Printf

export ctsem_sample

"""
    _sample_initial_point(sampler, values, metric, rng, scale, logdensity!)

A starting position drawn from the Laplace approximation itself.

Chains must not all start at the same point -- R-hat compares between-chain to
within-chain variance, and identical starts make it meaningless until the chains
have forgotten them. Drawing from the fitted approximation gives genuinely
different starts that are nonetheless all in the typical set, which is what the
approximation is *for*. Falling back toward the mode on a non-finite draw keeps
a bad approximation from preventing a start altogether.
"""
function _sample_initial_point(centre::Vector{Float64},
    metric::CTSEMMetric, rng::AbstractRNG, scale::Float64, logdensity!,
    gradient::Vector{Float64})
    x = copy(centre)
    logp = logdensity!(gradient, x)
    isfinite(logp) || throw(ArgumentError(
        "the log density is not finite at the supplied parameter values, so no " *
        "chain can start there; check the estimate the sample was asked to start from"))
    scale <= 0 && return (x, logp)
    draw = zeros(length(centre))
    for attempt in 0:6
        shrink = scale * 0.5^attempt
        for b in eachindex(metric.ranges)
            r = metric.ranges[b]
            C = metric.factors[b]
            k = length(r)
            for i in 1:k
                draw[r[i]] = randn(rng)
            end
            for i in k:-1:1
                acc = 0.0
                for j in 1:i
                    acc += C[i, j] * draw[r[j]]
                end
                x[r[i]] = centre[r[i]] + shrink * acc
            end
        end
        candidate = logdensity!(gradient, x)
        isfinite(candidate) && return (x, candidate)
    end
    copyto!(x, centre)
    return (x, logdensity!(gradient, x))
end

"""The largest or smallest finite entry, or NaN when there is none."""
function _finite_extremum(values::AbstractVector{Float64}, reduce)
    finite = filter(isfinite, values)
    return isempty(finite) ? NaN : reduce(finite)
end

"""
One chain's output, and enough of its state to continue it later.

`x`, `stepsize` and `metric` are the whole of what a chain needs to resume: the
sampler is Markov in the position, and the tuning is the step size and the
metric. The random stream is deliberately *not* preserved -- continuing on a
fresh stream is still a valid chain, and reconstructing a `Xoshiro` across the R
boundary would be a lot of fragile machinery for no statistical gain.
"""
struct _ChainResult
    draws::Matrix{Float64}      # ndim x ndraws
    accept::Vector{Float64}
    divergent::Vector{Bool}
    depth::Vector{Int}
    energy::Vector{Float64}
    stepsize::Float64
    warmup_divergent::Int
    warmup_used::Int
    x::Vector{Float64}
    metric::CTSEMMetric
end

# How often warmup re-examines its divergence rate, what rate counts as
# persisting, and how far the acceptance target may be pushed. The rate is
# deliberately well above zero: an isolated divergence early in warmup, at a
# step size dual averaging is still moving, says nothing worth acting on.
const _ACCEPT_CHECK_STRIDE = 50
const _ACCEPT_DIVERGENCE_RATE = 0.05
const _ACCEPT_STEP = 0.05
const _ACCEPT_MAX = 0.95
const _ACCEPT_MAX_RAISES = 3

function _run_chain(logdensity!, centre::Vector{Float64},
    metric::CTSEMMetric, rng::AbstractRNG,
    nwarmup::Int, ndraws::Int, maxdepth::Int, target_accept::Float64,
    maxdelta::Float64, init_scale::Float64, adapt_metric::Bool,
    adapt::Union{Nothing,Vector{Bool}}; settle_tol::Float64=0.0,
    resume::Union{Nothing,_ChainResult}=nothing,
    progress::CTSEMProgress=CTSEMProgress(false),
    callback::CTSEMCallback=CTSEMCallback(nothing))

    ndim = length(centre)
    ws = _NUTSWorkspace(ndim, maxdepth)
    g = zeros(ndim)
    # Continuing a chain: its position, step size and metric are the whole of
    # its state, so warmup is skipped entirely rather than repeated.
    if resume !== nothing
        x = copy(resume.x)
        logp = logdensity!(g, x)
        return _continue_chain(logdensity!, ws, rng, x, g, logp,
            resume.stepsize, resume.metric, ndraws, maxdepth, maxdelta, progress,
            callback)
    end
    x, logp = _sample_initial_point(centre, metric, rng, init_scale,
        logdensity!, g)

    current = metric
    eps = _init_stepsize(logdensity!, current, rng, x, g, logp, ws)
    da = _DualAverage(eps, target_accept)
    # Divergences since the last checkpoint, for the automatic raise below.
    check_divergent = 0
    check_start = 0
    raises = 0
    windows = adapt_metric ? _adapt_windows(nwarmup) : Int[]
    window_draws = Vector{Vector{Float64}}()
    warmup_divergent = 0
    warmup_used = nwarmup
    settled = 0

    iteration = 0
    depths = 0
    while iteration < nwarmup
        iteration += 1
        step = _nuts_transition!(ws, logdensity!, current, rng, x, g, logp, eps,
            maxdepth, maxdelta)
        logp = step.logp
        depths += step.depth
        step.divergent && (warmup_divergent += 1)
        step.divergent && (check_divergent += 1)
        eps = _dual_update!(da, step.accept)
        # Raise `target_accept` rather than only reporting divergences.
        #
        # A divergence means the integrator could not follow the geometry at the
        # step size it was using, and the standard remedy -- a higher acceptance
        # target, hence shorter steps -- was until now something the user had to
        # apply by hand, after the run, having read a warning about it. Warmup
        # already holds what is needed to apply it *during* the run, which is
        # when it is worth something: the draws a too-large step size ruined are
        # warmup draws, discarded either way.
        #
        # Checked on a fixed stride rather than at the metric windows, so it
        # still works when warmup is too short to hold a window and the Laplace
        # metric is used as it stands. Bounded at 0.95 and at three raises: past
        # that the steps are short enough that trajectories lengthen to
        # compensate, and geometry surviving 0.95 wants a reparameterisation
        # rather than a smaller step. `_dual_restart!` because the averaging is
        # chasing a new target from here, and its accumulated `hbar` is evidence
        # about the old one.
        if iteration - check_start >= _ACCEPT_CHECK_STRIDE
            rate = check_divergent / (iteration - check_start)
            if rate > _ACCEPT_DIVERGENCE_RATE && da.target < _ACCEPT_MAX &&
                    raises < _ACCEPT_MAX_RAISES
                da.target = min(_ACCEPT_MAX, da.target + _ACCEPT_STEP)
                raises += 1
                eps = _dual_restart!(da, eps)
            end
            check_divergent = 0
            check_start = iteration
        end
        if _due(progress)
            _progress_line(progress, iteration, nwarmup,
                @sprintf("logp %11.2f", logp), @sprintf("eps %8.2e", eps),
                @sprintf("depth %4.1f", depths / iteration),
                @sprintf("div %d", warmup_divergent),
                # Shown only once it has moved, so the ordinary run reads as it
                # always did.
                raises > 0 ? @sprintf("accept %.2f", da.target) : "")
        end
        # Its own cadence; see `ctsem_optimize`. Unconditional -- the callback
        # rate-limits itself -- so a GUI watching a chain sees warmup progress
        # even when nothing is being printed.
        _invoke_callback(callback, "warmup", iteration, nwarmup, logp,
            warmup_divergent)
        isempty(windows) && continue
        push!(window_draws, copy(x))
        iteration in windows || continue
        # Re-estimate from this window only. Earlier draws were taken under a
        # different metric and a different step size, and pooling them biases
        # the estimate toward wherever the chain used to be rather than where
        # it is.
        previous = current
        current = _estimate_metric(window_draws, current; adapt=adapt)
        empty!(window_draws)
        eps = _dual_restart!(da, _init_stepsize(logdensity!, current, rng, x, g,
            logp, ws))
        # Stop warming up once the metric has stopped moving.
        #
        # Stan's schedule runs a fixed number of iterations whatever happens,
        # which is the safe default and often twice what a well-started chain
        # needs -- and here the chain *is* well started, from a fit, with a
        # metric read off its curvature. Two consecutive windows that agree to
        # within `settle_tol` mean the remaining windows would re-estimate the
        # same matrix, so the terminal buffer is run and warmup ends.
        #
        # **Measured, and it loses badly.** On the N=200 augmented marginal
        # route: adaptive warmup took 1795 s to reach min ESS 142.8 (hitting an
        # 8000-draw budget without meeting its targets), against 559 s for min
        # ESS 246.3 on the fixed schedule -- 0.080 against 0.441 ESS/s, a factor
        # of 5.5 the wrong way. The reason the reasoning above is wrong: a
        # metric that has stopped *moving* is not the same as a metric that is
        # *good*. Two agreeing windows early mean the estimate has converged to
        # what a few hundred draws can say, and that estimate is then used for
        # every remaining draw, so the sampling phase pays for the whole run
        # what the warmup saved once. Hence `settle_tol = 0` -- off -- as the
        # default, and this comment rather than a removal, so the idea is not
        # re-invented and re-measured.
        if settle_tol > 0 && _metric_settled(previous, current, settle_tol)
            settled += 1
            if settled >= 2
                remaining = min(nwarmup - iteration, _ADAPT_TERM_BUFFER)
                for _ in 1:remaining
                    iteration += 1
                    step = _nuts_transition!(ws, logdensity!, current, rng, x, g,
                        logp, eps, maxdepth, maxdelta)
                    logp = step.logp
                    step.divergent && (warmup_divergent += 1)
                    eps = _dual_update!(da, step.accept)
                end
                warmup_used = iteration
                break
            end
        else
            settled = 0
        end
    end
    eps = _dual_final(da)
    # Forced: warmup may finish inside a callback interval (a short run, or
    # the settle_tol early-stop above), and the transition to sampling is
    # exactly the state a live watcher should not miss.
    _invoke_callback(callback, "warmup", warmup_used, nwarmup, logp,
        warmup_divergent; force=true)

    draws = Matrix{Float64}(undef, ndim, ndraws)
    accept = Vector{Float64}(undef, ndraws)
    divergent = Vector{Bool}(undef, ndraws)
    depth = Vector{Int}(undef, ndraws)
    energy = Vector{Float64}(undef, ndraws)
    progress.label = "sampling"
    progress.started = time()
    for iteration in 1:ndraws
        step = _nuts_transition!(ws, logdensity!, current, rng, x, g, logp, eps,
            maxdepth, maxdelta)
        logp = step.logp
        @inbounds for j in 1:ndim
            draws[j, iteration] = x[j]
        end
        accept[iteration] = step.accept
        divergent[iteration] = step.divergent
        depth[iteration] = step.depth
        energy[iteration] = step.energy
        if _due(progress)
            _progress_line(progress, iteration, ndraws,
                @sprintf("logp %11.2f", logp),
                @sprintf("depth %4.1f", sum(view(depth, 1:iteration)) / iteration),
                @sprintf("div %d", count(view(divergent, 1:iteration))))
        end
        _invoke_callback(callback, "sampling", iteration, ndraws, logp,
            count(view(divergent, 1:iteration)))
    end
    # Forced for the same reason as the optimiser's final call: a rate-limited
    # callback on a chain that finishes inside one interval would otherwise
    # never report the finished state at all.
    _invoke_callback(callback, "sampling", ndraws, ndraws, logp,
        count(divergent); force=true)
    return _ChainResult(draws, accept, divergent, depth, energy, eps,
        warmup_divergent, warmup_used, copy(x), current)
end

"""
    _continue_chain(...)

Draw another batch from a chain that has already warmed up.

No warmup and no adaptation: the step size and metric come in fixed, which is
what makes this a continuation of the same chain rather than a new one that
happens to start nearby. Adapting again would break detailed balance across the
join, because the transition kernel would then depend on the draws it had
already produced.
"""
function _continue_chain(logdensity!, ws::_NUTSWorkspace, rng::AbstractRNG,
    x::Vector{Float64}, g::Vector{Float64}, logp::Float64, eps::Float64,
    metric::CTSEMMetric, ndraws::Int, maxdepth::Int, maxdelta::Float64,
    progress::CTSEMProgress=CTSEMProgress(false),
    callback::CTSEMCallback=CTSEMCallback(nothing))
    progress.label = "sampling"
    progress.started = time()
    ndim = length(x)
    draws = Matrix{Float64}(undef, ndim, ndraws)
    accept = Vector{Float64}(undef, ndraws)
    divergent = Vector{Bool}(undef, ndraws)
    depth = Vector{Int}(undef, ndraws)
    energy = Vector{Float64}(undef, ndraws)
    for iteration in 1:ndraws
        step = _nuts_transition!(ws, logdensity!, metric, rng, x, g, logp, eps,
            maxdepth, maxdelta)
        logp = step.logp
        @inbounds for j in 1:ndim
            draws[j, iteration] = x[j]
        end
        accept[iteration] = step.accept
        divergent[iteration] = step.divergent
        depth[iteration] = step.depth
        energy[iteration] = step.energy
        if _due(progress)
            _progress_line(progress, iteration, ndraws,
                @sprintf("logp %11.2f", logp),
                @sprintf("depth %4.1f", sum(view(depth, 1:iteration)) / iteration),
                @sprintf("div %d", count(view(divergent, 1:iteration))))
        end
        _invoke_callback(callback, "sampling", iteration, ndraws, logp,
            count(view(divergent, 1:iteration)))
    end
    _invoke_callback(callback, "sampling", ndraws, ndraws, logp,
        count(divergent); force=true)
    return _ChainResult(draws, accept, divergent, depth, energy, eps, 0, 0,
        copy(x), metric)
end

"""
    _sample_to_target(...)

Sample until the effective sample size targets are met, or the draw budget runs
out, whichever comes first.

Draws arrive in batches, and after each one the diagnostics are recomputed over
*everything so far*. Effective sample size and R-hat are properties of the whole
run rather than of the latest batch, and stopping on a batch-local figure would
be stopping on noise. Between batches every chain continues from where it
stopped with its own step size and metric held fixed, so the result is one run
of `total` draws rather than a concatenation of short ones.

The per-chain state comes back with the draws, which is the same mechanism a
caller uses to ask for more after looking at the answer.
"""
function _sample_to_target(density_for, centre, metric, nchains::Int,
    parallel::Bool, seed::Integer, nwarmup::Int, ndraws::Int, maxdepth::Int,
    target_accept::Float64, maxdelta::Float64, init_scale::Float64,
    adapt_metric::Bool, adapt, settle_tol::Float64, min_ess::Float64,
    mean_ess::Float64, max_draws::Int, rhat_target::Float64, npar::Int,
    resume, verbose::Bool, overwrite::Bool=true; progress_callback=nothing,
    progress_sink=nothing)

    results = if resume === nothing
        _sample_chains(nchains, parallel, seed, centre, metric, nwarmup, ndraws,
            maxdepth, target_accept, maxdelta, init_scale, adapt_metric, adapt,
            density_for; settle_tol=settle_tol, progress=verbose,
            overwrite=overwrite, progress_callback=progress_callback,
            progress_sink=progress_sink)
    else
        _continue_chains(nchains, parallel, seed, ndraws, maxdepth, maxdelta,
            density_for, resume; progress=verbose, overwrite=overwrite,
            progress_callback=progress_callback, progress_sink=progress_sink)
    end
    total = ndraws
    attempt = 0
    was_met = false
    while min_ess > 0 || mean_ess > 0
        pooled = _pool_draws(results, npar)
        diagnostics = ctsem_sample_diagnostics(pooled, nchains)
        finite_ess = filter(isfinite, diagnostics.ess)
        finite_rhat = filter(isfinite, diagnostics.rhat)
        worst = isempty(finite_ess) ? 0.0 : minimum(finite_ess)
        average = isempty(finite_ess) ? 0.0 : sum(finite_ess) / length(finite_ess)
        rhat = isempty(finite_rhat) ? Inf : maximum(finite_rhat)
        met = worst >= min_ess && average >= mean_ess && rhat <= rhat_target
        # Confirmed once before stopping. Stopping the moment a target is first
        # met is a rule correlated with the quantity it tests: effective size is
        # estimated with error, so a first crossing is more often a favourable
        # error than a real one, and the realised size settles below target. One
        # extra batch removes most of that, and costs one batch.
        confirmed = met && was_met
        if verbose
            println(_console(), "  ", total, " draws per chain: min ESS ",
                round(worst; digits=1), ", mean ESS ", round(average; digits=1),
                ", worst R-hat ", round(rhat; digits=3),
                confirmed ? " -- targets met" :
                met ? " -- targets met, confirming" : "")
        end
        was_met = met
        confirmed && break
        if total >= max_draws
            verbose && println(_console(), "  draw budget of ", max_draws,
                " per chain reached before the targets were met")
            break
        end
        # Ask for as many more as the shortfall suggests, bounded below by a
        # quarter of the last batch -- a tiny follow-up costs a round of
        # scheduling and moves the estimate hardly at all -- and above by four
        # times it, so one badly mixing coordinate cannot demand an enormous
        # single batch on the strength of an early, noisy estimate.
        shortfall = max(min_ess / max(worst, 1.0), mean_ess / max(average, 1.0))
        wanted = clamp(ceil(Int, ndraws * (shortfall - 1)), fld(ndraws, 4), ndraws * 4)
        wanted = min(wanted, max_draws - total)
        wanted <= 0 && break
        attempt += 1
        results = _merge_chains(results,
            _continue_chains(nchains, parallel, seed + 1000 * attempt, wanted,
                maxdepth, maxdelta, density_for, results; progress=verbose,
                overwrite=overwrite, progress_callback=progress_callback,
                progress_sink=progress_sink))
        total += wanted
    end
    return (results=results, ndraws=total)
end

"""Draws from every chain, population part only, chain-major."""
function _pool_draws(results::Vector{_ChainResult}, npar::Int)
    ndraws = size(first(results).draws, 2)
    pooled = Matrix{Float64}(undef, npar, length(results) * ndraws)
    for (c, r) in enumerate(results)
        @inbounds for t in 1:ndraws, j in 1:npar
            pooled[j, (c - 1) * ndraws + t] = r.draws[j, t]
        end
    end
    return pooled
end

"""Continue every chain from its own state."""
function _continue_chains(nchains::Int, parallel::Bool, seed::Integer,
    ndraws::Int, maxdepth::Int, maxdelta::Float64, density_for,
    previous::Vector{_ChainResult}; progress::Bool=false,
    overwrite::Bool=true, progress_callback=nothing, progress_sink=nothing)
    results = Vector{_ChainResult}(undef, nchains)
    runner = function (c)
        reporter = CTSEMProgress(progress && c == 1; label="sampling",
            overwrite=overwrite, sink=progress_sink)
        # Only chain 1 gets a live callback too, and for the same reason as
        # the printed line: several threads calling back into R at once is
        # not merely unreadable, it is unsafe. See `_sample_chains`.
        watcher = CTSEMCallback(c == 1 ? progress_callback : nothing)
        results[c] = _run_chain(density_for(c), previous[c].x,
            previous[c].metric, Random.Xoshiro(UInt64(seed) + UInt64(c)),
            0, ndraws, maxdepth, 0.8, maxdelta, 0.0, false, nothing;
            resume=previous[c], progress=reporter, callback=watcher)
        # Close the line. Nothing did, so the last in-place update was left open
        # with no newline on it and whatever R printed next landed inside it --
        # "div 0Laplace fit: trajectories are conditional...". The optimiser
        # routes have closed theirs for a while; this one never has.
        _progress_done(reporter, @sprintf("%d draws", ndraws))
        return nothing
    end
    if parallel
        Threads.@sync for c in 1:nchains
            Threads.@spawn runner(c)
        end
    else
        for c in 1:nchains
            runner(c)
        end
    end
    return results
end

"""Glue two batches of the same chains into one."""
function _merge_chains(a::Vector{_ChainResult}, b::Vector{_ChainResult})
    return [_ChainResult(hcat(a[c].draws, b[c].draws),
        vcat(a[c].accept, b[c].accept), vcat(a[c].divergent, b[c].divergent),
        vcat(a[c].depth, b[c].depth), vcat(a[c].energy, b[c].energy),
        b[c].stepsize, a[c].warmup_divergent, a[c].warmup_used,
        b[c].x, b[c].metric) for c in eachindex(a)]
end

"""Have two successive metric estimates stopped disagreeing?"""

function _metric_settled(a::CTSEMMetric, b::CTSEMMetric, tol::Real)
    length(a.factors) == length(b.factors) || return false
    for k in eachindex(a.factors)
        A = a.factors[k]
        B = b.factors[k]
        size(A) == size(B) || return false
        # Every entry, not just the diagonal. The diagonals are the scale the
        # step size is tuned against, which is true and not sufficient: the
        # off-diagonal structure is the whole reason for a block metric, and it
        # goes on being learned after the scales have settled.
        for j in axes(A, 2), i in axes(A, 1)
            reference = max(abs(A[j, j]), abs(B[j, j]), eps(Float64))
            abs(A[i, j] - B[i, j]) / reference > tol && return false
        end
    end
    return true
end

"""Iterations of step-size-only adaptation after the metric has settled."""
const _ADAPT_TERM_BUFFER = 50

"""
    _sample_chains(nchains, parallel, seed, centre, metric, ..., density_for)

Run the chains, concurrently when there are threads for them.

`density_for(c)` builds chain `c`'s log-density closure. It is a function of the
chain rather than one shared closure because a target may need per-chain
scratch: the joint sampler hands each chain its own adjoint workspace, since
every chain filters every subject and they would otherwise share buffers.
"""
function _sample_chains(nchains::Int, parallel::Bool, seed::Integer,
    centre::Vector{Float64}, metric::CTSEMMetric, nwarmup::Int, ndraws::Int,
    maxdepth::Int, target_accept::Float64, maxdelta::Float64,
    init_scale::Float64, adapt_metric::Bool, adapt::Union{Nothing,Vector{Bool}},
    density_for; settle_tol::Float64=0.0, progress::Bool=false,
    overwrite::Bool=true, progress_callback=nothing, progress_sink=nothing)
    results = Vector{_ChainResult}(undef, nchains)
    runner = function (c)
        # Only the first chain reports. Four threads writing lines interleave
        # into something unreadable, and a lock to prevent that would serialise
        # the work being reported on; one chain is representative when they are
        # all doing the same thing.
        reporter = CTSEMProgress(progress && c == 1; label="warmup",
            overwrite=overwrite, sink=progress_sink)
        # Same restriction on the callback, and for a sharper reason than
        # readability: several `Threads.@spawn`ed chains calling back into R
        # at once is a concurrency hazard, not just noise. One representative
        # chain is what a GUI gets, exactly as one representative chain is
        # what the console gets.
        watcher = CTSEMCallback(c == 1 ? progress_callback : nothing)
        results[c] = _run_chain(density_for(c), centre, metric,
            Random.Xoshiro(UInt64(seed) + UInt64(c)), nwarmup, ndraws, maxdepth,
            target_accept, maxdelta, init_scale, adapt_metric, adapt;
            settle_tol=settle_tol, progress=reporter, callback=watcher)
        # See `_continue_chains`. One reporter spans both phases -- `_run_chain`
        # relabels it from "warmup" to "sampling" partway -- so the closing line
        # names both rather than whichever phase it ended in.
        _progress_done(reporter,
            @sprintf("%d warmup + %d draws", nwarmup, ndraws))
        return nothing
    end
    if parallel
        Threads.@sync for c in 1:nchains
            Threads.@spawn runner(c)
        end
    else
        for c in 1:nchains
            runner(c)
        end
    end
    return results
end

"""
    ctsem_sample(laplace, values; kwargs...)

Sample the joint posterior over population parameters and random effects.

`values` is where to start -- a Laplace estimate, normally, which is also what
the initial metric is read from. Everything else has a Stan-compatible default:
four chains, 500 warmup and 500 sampling iterations each, target acceptance 0.8,
maximum tree depth 10.

Returns population draws as `npar x (nchains * ndraws)`, chain-major, with the
random effects summarised rather than returned unless `save_effects` is set.
Diagnostics come back alongside: split R-hat and effective sample size per
parameter, divergences, tree depths, step sizes and E-BFMI.

`adapt_metric=false` keeps the Laplace metric throughout, which is worth trying
when warmup is short: it is already a good metric, and re-estimating it from a
few hundred draws can be worse than leaving it alone.
"""
function ctsem_sample(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    npar::Integer=length(values), nchains::Integer=4, nwarmup::Integer=500,
    ndraws::Integer=500, maxdepth::Integer=10, target_accept::Real=0.8,
    maxdelta::Real=1000.0, seed::Integer=20260828, init_scale::Real=1.0,
    adapt_metric::Bool=true, adapt_effects::Bool=false, save_effects::Bool=false,
    hessian::Union{Nothing,AbstractMatrix}=nothing, verbose::Bool=false,
    min_ess::Real=0.0, mean_ess::Real=0.0, max_draws::Integer=0,
    rhat_target::Real=1.01, settle_tol::Real=0.0, resume=nothing,
    progress_overwrite::Bool=true, progress_callback=nothing,
    progress_sink=nothing)

    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    nchains >= 1 || throw(ArgumentError("nchains must be positive"))
    ndraws >= 1 || throw(ArgumentError("ndraws must be positive"))
    nwarmup >= 0 || throw(ArgumentError("nwarmup must be non-negative"))
    0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0, 1)"))

    sampler = ctsem_sampler(laplace, npar)
    start = collect(Float64, values)
    # Builds the metric *and* solves each unit's conditional mode at `start`,
    # leaving it on the Laplace object, which `ctsem_sample_start` then uses to
    # place the effects. So the order of these two lines matters: reversed, the
    # effects are placed at whatever the objective's last call left behind --
    # zero on a fresh one -- and `ctsem_sample_metric` records what that was
    # measured to do to a chain.
    metric = ctsem_sample_metric(sampler, start; hessian=hessian)

    # Grown before anything is spawned: a concurrent push! onto the shared
    # workspace vector is a race, and one chain per slot is the whole reason
    # chains can run at all.
    parallel = nchains > 1 && Threads.nthreads() > 1
    # Chains first, then the unit loop with whatever threads remain. Two chains
    # on ten threads previously used two; they now take five apiece. Bounded by
    # `ctsem_set_max_chunks!` as well, so `cores` still caps the total.
    per_chain = parallel ?
        max(1, min(ctsem_max_chunks().max_chunks,
                   Threads.nthreads() ÷ nchains)) : 1
    if parallel
        while length(laplace.workspaces) < nchains * per_chain
            push!(laplace.workspaces, Dict{Any,Any}())
        end
    end
    verbose && println(_console(), "Sampling: ", nchains, " chain(s), ", ctsem_sample_dimension(sampler),
        " dimensions (", sampler.npar, " population + ",
        ctsem_sample_dimension(sampler) - sampler.npar, " effects), ",
        parallel ? (per_chain > 1 ?
            string(per_chain, " threads each") : "one thread each") :
            "unit-parallel", ", metric in ",
        length(metric.ranges), " block(s)")

    # Which metric blocks warmup may re-estimate. The population block always
    # benefits: its Laplace value is a local quadratic fit and the posterior is
    # not quadratic. The effect blocks come from a conditional covariance that
    # is *exact* for a linear model, and replacing one with a k x k estimate
    # from a few hundred draws is as likely to add noise as to remove bias --
    # measured at 356 effective draws keeping them against 261 re-estimating.
    adapt = adapt_effects ? nothing : [b == 1 for b in eachindex(metric.ranges)]
    centre = ctsem_sample_start(sampler, start)
    run = _sample_to_target(
        c -> ((g, x) -> ctsem_sample_density!(g, sampler, x;
            workspace_slot=parallel ? (c - 1) * per_chain + 1 : nothing,
            workspace_chunks=per_chain)),
        centre, metric, nchains, parallel, seed, nwarmup, ndraws, Int(maxdepth),
        Float64(target_accept), Float64(maxdelta), Float64(init_scale),
        adapt_metric, adapt, Float64(settle_tol), Float64(min_ess),
        Float64(mean_ess), max(Int(max_draws), ndraws), Float64(rhat_target),
        sampler.npar, resume, verbose, progress_overwrite;
        progress_callback=progress_callback, progress_sink=progress_sink)
    results = run.results
    ndraws = run.ndraws

    ndim = sampler.ndim
    total = nchains * ndraws
    kept = save_effects ? ndim : sampler.npar
    draws = Matrix{Float64}(undef, kept, total)
    accept = Vector{Float64}(undef, total)
    divergent = Vector{Bool}(undef, total)
    depth = Vector{Int}(undef, total)
    energy = Vector{Float64}(undef, total)
    # Running mean and sum of squares for the effects, so their summary costs
    # nothing whether or not the draws themselves are kept.
    effect_sum = zeros(Float64, ndim - sampler.npar)
    effect_sq = zeros(Float64, ndim - sampler.npar)
    for c in 1:nchains
        r = results[c]
        for t in 1:ndraws
            column = (c - 1) * ndraws + t
            @inbounds for j in 1:kept
                draws[j, column] = r.draws[j, t]
            end
            @inbounds for j in (sampler.npar + 1):ndim
                v = r.draws[j, t]
                effect_sum[j - sampler.npar] += v
                effect_sq[j - sampler.npar] += v * v
            end
            accept[column] = r.accept[t]
            divergent[column] = r.divergent[t]
            depth[column] = r.depth[t]
            energy[column] = r.energy[t]
        end
    end
    effect_mean = effect_sum ./ total
    effect_sd = sqrt.(max.(effect_sq ./ total .- effect_mean .^ 2, 0.0) .*
        (total / max(total - 1, 1)))

    diagnostics = ctsem_sample_diagnostics(draws, nchains)
    return (
        draws=draws,
        npar=sampler.npar,
        ndim=ndim,
        nchains=nchains,
        ndraws=ndraws,
        saved_effects=save_effects,
        effect_mean=effect_mean,
        effect_sd=effect_sd,
        rhat=diagnostics.rhat,
        ess=diagnostics.ess,
        accept=accept,
        divergent=divergent,
        depth=depth,
        energy=energy,
        stepsize=[r.stepsize for r in results],
        warmup_divergent=[r.warmup_divergent for r in results],
        ndivergent=count(divergent),
        max_depth=Int(maxdepth),
        nsaturated=count(==(Int(maxdepth)), depth),
        ebfmi=_ebfmi(energy, nchains),
        # `init=NaN` would be wrong rather than defensive: Julia's `min` and
        # `max` propagate NaN, so it poisons every result instead of only the
        # empty one.
        worst_rhat=_finite_extremum(diagnostics.rhat, maximum),
        min_ess=_finite_extremum(diagnostics.ess, minimum),
    )
end

export ctsem_sample_marginal

"""
    ctsem_sample_marginal(objective, values; kwargs...)

Sample the population parameters with the random effects already integrated out.

The counterpart to `ctsem_sample`, and usually the faster one. Where that
samples `theta` and every subject's effects jointly -- `npar + sum_U dim(u_U)`
coordinates, so 207 for two hundred subjects with one effect each -- this
samples `theta` alone, because the objective it is given has already done the
integral. The dimension therefore does not grow with the subject count at all.

Two objectives qualify and both work here unchanged, because `ctsem_evaluate`
and `ctsem_hessian` are defined for each:

  * a `CTSEMObjective` built from an augmented model, where the filter carries
    the effects as states and integrates them analytically;
  * a `CTSEMLaplaceObjective`, where the Laplace approximation does it.

The geometry is much kinder. The joint target has a `L(theta) u` product in it,
which is the funnel that leaves E-BFMI at 0.16-0.49 and costs most of the
sampling efficiency; a `theta`-only target has no such product. Measured against
Stan on identical data, Stan samples the augmented marginal in 674s where the
joint sampler takes 1114s at two hundred subjects -- and the two agree on every
posterior mean and standard deviation to three decimals, which is what says the
difference is dimension rather than correctness.

# Chains run one after another here

The joint sampler gives each chain its own adjoint workspace, because every
chain filters every subject and they would otherwise share buffers. `ctsem_
evaluate` has no such per-chain slot: its workspaces are per *chunk*, and chunks
partition subjects, which is safe for one caller and a race for several. So
chains are sequential and each gradient uses the subject-loop parallelism
instead. That is the smaller of the two axes -- about twofold against linear --
and it is a real cost, but the dimension collapse is worth far more than the
axis is.

# The Laplace objective is stateful, and that matters more here

`CTSEMLaplaceObjective` retains each unit's inner mode between calls and
warm-starts the next solve from it. Across an optimizer's trajectory that is
free accuracy; across a sampler's it makes the density a function of where the
chain has *been* as well as where it is. The inner problem is concave in `u` for
these models, so the warm start converges to the same mode either way and the
dependence does not bite -- but it is an assumption rather than a guarantee, and
`test_sampler.jl`'s "the inner mode does not depend on how the chain got there"
pins it: the same `theta`, evaluated after three different histories (a nearby
warm start, a warm start from far away, and a cold start with no prior call at
all), gives the same value, gradient and retained modes to numerical tolerance.
"""
function ctsem_sample_marginal(objective, values::AbstractVector;
    nchains::Integer=4, nwarmup::Integer=500, ndraws::Integer=500,
    maxdepth::Integer=10, target_accept::Real=0.8, maxdelta::Real=1000.0,
    seed::Integer=20260828, init_scale::Real=1.0, adapt_metric::Bool=true,
    hessian::Union{Nothing,AbstractMatrix}=nothing, gradient_method=:adjoint,
    verbose::Bool=false, min_ess::Real=0.0, mean_ess::Real=0.0,
    max_draws::Integer=0, rhat_target::Real=1.01, settle_tol::Real=0.0,
    resume=nothing, progress_overwrite::Bool=true, progress_callback=nothing,
    progress_sink=nothing)

    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    nchains >= 1 || throw(ArgumentError("nchains must be positive"))
    ndraws >= 1 || throw(ArgumentError("ndraws must be positive"))
    nwarmup >= 0 || throw(ArgumentError("nwarmup must be non-negative"))
    0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0, 1)"))

    centre = collect(Float64, values)
    npar = length(centre)
    logdensity! = function (g, x)
        result = try
            ctsem_evaluate(objective, x; gradient=true,
                gradient_method=gradient_method)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        if result === nothing || !isfinite(result.value) ||
            result.gradient === nothing || !all(isfinite, result.gradient)
            fill!(g, 0.0)
            return -Inf
        end
        copyto!(g, result.gradient)
        return result.value
    end

    # One dense block: at this dimension the whole covariance is affordable to
    # estimate and to factorize, and there is no sparsity to exploit -- the
    # population parameters are all coupled.
    H = hessian === nothing ? ctsem_hessian(objective, centre) : Matrix(hessian)
    information = Symmetric((-(H .+ transpose(H))) ./ 2)
    metric = _metric_from_covariances([1:npar], [_bounded_inverse(information)])

    verbose && println(_console(), "Sampling: ", nchains, " chain(s), ", npar,
        " dimensions (effects integrated out), chains sequential, ",
        "metric in 1 block")

    run = _sample_to_target(_ -> logdensity!, centre, metric, nchains, false,
        seed, nwarmup, ndraws, Int(maxdepth), Float64(target_accept),
        Float64(maxdelta), Float64(init_scale), adapt_metric, nothing,
        Float64(settle_tol), Float64(min_ess), Float64(mean_ess),
        max(Int(max_draws), ndraws), Float64(rhat_target), npar, resume,
        verbose, progress_overwrite; progress_callback=progress_callback,
        progress_sink=progress_sink)
    results = run.results
    ndraws = run.ndraws

    total = nchains * ndraws
    draws = Matrix{Float64}(undef, npar, total)
    accept = Vector{Float64}(undef, total)
    divergent = Vector{Bool}(undef, total)
    depth = Vector{Int}(undef, total)
    energy = Vector{Float64}(undef, total)
    for c in 1:nchains
        r = results[c]
        for t in 1:ndraws
            column = (c - 1) * ndraws + t
            @inbounds for j in 1:npar
                draws[j, column] = r.draws[j, t]
            end
            accept[column] = r.accept[t]
            divergent[column] = r.divergent[t]
            depth[column] = r.depth[t]
            energy[column] = r.energy[t]
        end
    end
    diagnostics = ctsem_sample_diagnostics(draws, nchains)
    return (
        draws=draws, npar=npar, ndim=npar, nchains=nchains, ndraws=ndraws,
        saved_effects=false, effect_mean=Float64[], effect_sd=Float64[],
        rhat=diagnostics.rhat, ess=diagnostics.ess,
        accept=accept, divergent=divergent, depth=depth, energy=energy,
        stepsize=[r.stepsize for r in results],
        warmup_divergent=[r.warmup_divergent for r in results],
        ndivergent=count(divergent), max_depth=Int(maxdepth),
        nsaturated=count(==(Int(maxdepth)), depth),
        ebfmi=_ebfmi(energy, nchains),
        worst_rhat=_finite_extremum(diagnostics.rhat, maximum),
        min_ess=_finite_extremum(diagnostics.ess, minimum),
    )
end
