"""
Running the sampler: chains, warmup, and what comes back.

# Chains are the parallel axis

The gradient parallelises over units and gains about twofold; chains share
nothing at all and gain linearly. So with several chains asked for, each takes a
thread and evaluates its own gradient serially in its own adjoint workspace, and
the unit-level chunking is switched off inside it. One chain gets the unit-level
parallelism instead. Nothing needs deciding by the caller: `workspace_slot` on
`ctsem_sample_density!` is present exactly when chains are running concurrently.

# What crosses the bridge

`save_effects` defaults to false, and that is a bandwidth decision rather than a
statistical one. A hundred subjects with two effects each and four chains of a
thousand draws is 800,000 numbers -- 6.4 MB, which JuliaConnectoR moves at
roughly 0.85 MB/s on Windows, so returning them costs longer than many fits do.
The per-effect posterior mean and standard deviation come back always, because
they are `ndim` numbers rather than `ndim * ndraws` and answer most of what the
draws would be used for.
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

function _run_chain(logdensity!, centre::Vector{Float64},
    metric::CTSEMMetric, rng::AbstractRNG,
    nwarmup::Int, ndraws::Int, maxdepth::Int, target_accept::Float64,
    maxdelta::Float64, init_scale::Float64, adapt_metric::Bool,
    adapt::Union{Nothing,Vector{Bool}}; settle_tol::Float64=0.0,
    resume::Union{Nothing,_ChainResult}=nothing,
    progress::CTSEMProgress=CTSEMProgress(false))

    ndim = length(centre)
    ws = _NUTSWorkspace(ndim, maxdepth)
    g = zeros(ndim)
    # Continuing a chain: its position, step size and metric are the whole of
    # its state, so warmup is skipped entirely rather than repeated.
    if resume !== nothing
        x = copy(resume.x)
        logp = logdensity!(g, x)
        return _continue_chain(logdensity!, ws, rng, x, g, logp,
            resume.stepsize, resume.metric, ndraws, maxdepth, maxdelta, progress)
    end
    x, logp = _sample_initial_point(centre, metric, rng, init_scale,
        logdensity!, g)

    current = metric
    eps = _init_stepsize(logdensity!, current, rng, x, g, logp, ws)
    da = _DualAverage(eps, target_accept)
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
        eps = _dual_update!(da, step.accept)
        if _due(progress)
            _progress_line(progress, iteration, nwarmup,
                @sprintf("logp %11.2f", logp), @sprintf("eps %8.2e", eps),
                @sprintf("depth %4.1f", depths / iteration),
                @sprintf("div %d", warmup_divergent))
        end
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
    end
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
    progress::CTSEMProgress=CTSEMProgress(false))
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
    end
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
    resume, verbose::Bool, overwrite::Bool=true)

    results = if resume === nothing
        _sample_chains(nchains, parallel, seed, centre, metric, nwarmup, ndraws,
            maxdepth, target_accept, maxdelta, init_scale, adapt_metric, adapt,
            density_for; settle_tol=settle_tol, progress=verbose,
            overwrite=overwrite)
    else
        _continue_chains(nchains, parallel, seed, ndraws, maxdepth, maxdelta,
            density_for, resume; progress=verbose, overwrite=overwrite)
    end
    total = ndraws
    attempt = 0
    while min_ess > 0 || mean_ess > 0
        pooled = _pool_draws(results, npar)
        diagnostics = ctsem_sample_diagnostics(pooled, nchains)
        finite_ess = filter(isfinite, diagnostics.ess)
        finite_rhat = filter(isfinite, diagnostics.rhat)
        worst = isempty(finite_ess) ? 0.0 : minimum(finite_ess)
        average = isempty(finite_ess) ? 0.0 : sum(finite_ess) / length(finite_ess)
        rhat = isempty(finite_rhat) ? Inf : maximum(finite_rhat)
        met = worst >= min_ess && average >= mean_ess && rhat <= rhat_target
        if verbose
            println("  ", total, " draws per chain: min ESS ",
                round(worst; digits=1), ", mean ESS ", round(average; digits=1),
                ", worst R-hat ", round(rhat; digits=3),
                met ? " -- targets met" : "")
        end
        met && break
        if total >= max_draws
            verbose && println("  draw budget of ", max_draws,
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
                overwrite=overwrite))
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
    overwrite::Bool=true)
    results = Vector{_ChainResult}(undef, nchains)
    runner = function (c)
        reporter = CTSEMProgress(progress && c == 1; label="sampling",
            overwrite=overwrite)
        results[c] = _run_chain(density_for(c), previous[c].x,
            previous[c].metric, Random.Xoshiro(UInt64(seed) + UInt64(c)),
            0, ndraws, maxdepth, 0.8, maxdelta, 0.0, false, nothing;
            resume=previous[c], progress=reporter)
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
    overwrite::Bool=true)
    results = Vector{_ChainResult}(undef, nchains)
    runner = function (c)
        # Only the first chain reports. Four threads writing lines interleave
        # into something unreadable, and a lock to prevent that would serialise
        # the work being reported on; one chain is representative when they are
        # all doing the same thing.
        reporter = CTSEMProgress(progress && c == 1; label="warmup",
            overwrite=overwrite)
        results[c] = _run_chain(density_for(c), centre, metric,
            Random.Xoshiro(UInt64(seed) + UInt64(c)), nwarmup, ndraws, maxdepth,
            target_accept, maxdelta, init_scale, adapt_metric, adapt;
            settle_tol=settle_tol, progress=reporter)
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
    progress_overwrite::Bool=true)

    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    nchains >= 1 || throw(ArgumentError("nchains must be positive"))
    ndraws >= 1 || throw(ArgumentError("ndraws must be positive"))
    nwarmup >= 0 || throw(ArgumentError("nwarmup must be non-negative"))
    0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0, 1)"))

    sampler = ctsem_sampler(laplace, npar)
    start = collect(Float64, values)
    # Builds the metric *and* leaves each unit's conditional mode on the Laplace
    # object, which `ctsem_sample_start` then uses to place the effects. So the
    # order of these two lines matters.
    metric = ctsem_sample_metric(sampler, start; hessian=hessian)

    # Grown before anything is spawned: a concurrent push! onto the shared
    # workspace vector is a race, and one chain per slot is the whole reason
    # chains can run at all.
    parallel = nchains > 1 && Threads.nthreads() > 1
    if parallel
        while length(laplace.workspaces) < nchains
            push!(laplace.workspaces, Dict{Any,Any}())
        end
    end
    verbose && println("Sampling: ", nchains, " chain(s), ", ctsem_sample_dimension(sampler),
        " dimensions (", sampler.npar, " population + ",
        ctsem_sample_dimension(sampler) - sampler.npar, " effects), ",
        parallel ? "one thread each" : "unit-parallel", ", metric in ",
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
            workspace_slot=parallel ? c : nothing)),
        centre, metric, nchains, parallel, seed, nwarmup, ndraws, Int(maxdepth),
        Float64(target_accept), Float64(maxdelta), Float64(init_scale),
        adapt_metric, adapt, Float64(settle_tol), Float64(min_ess),
        Float64(mean_ess), max(Int(max_draws), ndraws), Float64(rhat_target),
        sampler.npar, resume, verbose, progress_overwrite)
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
`test_sampler.jl` pins it by evaluating the same `theta` from two different
histories and requiring the same answer.
"""
function ctsem_sample_marginal(objective, values::AbstractVector;
    nchains::Integer=4, nwarmup::Integer=500, ndraws::Integer=500,
    maxdepth::Integer=10, target_accept::Real=0.8, maxdelta::Real=1000.0,
    seed::Integer=20260828, init_scale::Real=1.0, adapt_metric::Bool=true,
    hessian::Union{Nothing,AbstractMatrix}=nothing, gradient_method=:adjoint,
    verbose::Bool=false, min_ess::Real=0.0, mean_ess::Real=0.0,
    max_draws::Integer=0, rhat_target::Real=1.01, settle_tol::Real=0.0,
    resume=nothing, progress_overwrite::Bool=true)

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

    verbose && println("Sampling: ", nchains, " chain(s), ", npar,
        " dimensions (effects integrated out), chains sequential, ",
        "metric in 1 block")

    run = _sample_to_target(_ -> logdensity!, centre, metric, nchains, false,
        seed, nwarmup, ndraws, Int(maxdepth), Float64(target_accept),
        Float64(maxdelta), Float64(init_scale), adapt_metric, nothing,
        Float64(settle_tol), Float64(min_ess), Float64(mean_ess),
        max(Int(max_draws), ndraws), Float64(rhat_target), npar, resume,
        verbose, progress_overwrite)
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
