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
function _sample_initial_point(sampler::CTSEMSampler, values::AbstractVector,
    metric::CTSEMMetric, rng::AbstractRNG, scale::Float64, logdensity!,
    gradient::Vector{Float64})
    centre = ctsem_sample_start(sampler, values)
    x = copy(centre)
    logp = logdensity!(gradient, x)
    isfinite(logp) || throw(ArgumentError(
        "the log density is not finite at the supplied parameter values, so no " *
        "chain can start there; check the estimate the sample was asked to start from"))
    scale <= 0 && return (x, logp)
    draw = zeros(sampler.ndim)
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

"""One chain's output, before the chains are stitched together."""
struct _ChainResult
    draws::Matrix{Float64}      # ndim x ndraws
    accept::Vector{Float64}
    divergent::Vector{Bool}
    depth::Vector{Int}
    energy::Vector{Float64}
    stepsize::Float64
    warmup_divergent::Int
end

function _run_chain(sampler::CTSEMSampler, values::AbstractVector,
    metric::CTSEMMetric, rng::AbstractRNG, slot::Union{Nothing,Int},
    nwarmup::Int, ndraws::Int, maxdepth::Int, target_accept::Float64,
    maxdelta::Float64, init_scale::Float64, adapt_metric::Bool)

    ndim = sampler.ndim
    logdensity! = (g, x) -> ctsem_sample_density!(g, sampler, x; workspace_slot=slot)
    ws = _NUTSWorkspace(ndim, maxdepth)
    g = zeros(ndim)
    x, logp = _sample_initial_point(sampler, values, metric, rng, init_scale,
        logdensity!, g)

    current = metric
    eps = _init_stepsize(logdensity!, current, rng, x, g, logp, ws)
    da = _DualAverage(eps, target_accept)
    windows = adapt_metric ? _adapt_windows(nwarmup) : Int[]
    window_draws = Vector{Vector{Float64}}()
    warmup_divergent = 0

    for iteration in 1:nwarmup
        step = _nuts_transition!(ws, logdensity!, current, rng, x, g, logp, eps,
            maxdepth, maxdelta)
        logp = step.logp
        step.divergent && (warmup_divergent += 1)
        eps = _dual_update!(da, step.accept)
        if !isempty(windows)
            push!(window_draws, copy(x))
            if iteration in windows
                # Re-estimate from this window only. Earlier draws were taken
                # under a different metric and a different step size, and
                # pooling them biases the estimate toward wherever the chain
                # used to be rather than where it is.
                current = _estimate_metric(window_draws, metric.ranges)
                empty!(window_draws)
                eps = _dual_restart!(da, _init_stepsize(logdensity!, current,
                    rng, x, g, logp, ws))
            end
        end
    end
    eps = _dual_final(da)

    draws = Matrix{Float64}(undef, ndim, ndraws)
    accept = Vector{Float64}(undef, ndraws)
    divergent = Vector{Bool}(undef, ndraws)
    depth = Vector{Int}(undef, ndraws)
    energy = Vector{Float64}(undef, ndraws)
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
    end
    return _ChainResult(draws, accept, divergent, depth, energy, eps, warmup_divergent)
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
    adapt_metric::Bool=true, save_effects::Bool=false,
    hessian::Union{Nothing,AbstractMatrix}=nothing, verbose::Bool=false)

    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    nchains >= 1 || throw(ArgumentError("nchains must be positive"))
    ndraws >= 1 || throw(ArgumentError("ndraws must be positive"))
    nwarmup >= 0 || throw(ArgumentError("nwarmup must be non-negative"))
    0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0, 1)"))

    sampler = ctsem_sampler(laplace, npar)
    start = collect(Float64, values)
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

    results = Vector{_ChainResult}(undef, nchains)
    runner = function (c)
        results[c] = _run_chain(sampler, start, metric,
            Random.Xoshiro(UInt64(seed) + UInt64(c)), parallel ? c : nothing,
            nwarmup, ndraws, Int(maxdepth), Float64(target_accept),
            Float64(maxdelta), Float64(init_scale), adapt_metric)
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
