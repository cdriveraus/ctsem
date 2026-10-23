"""
Tempered sequential Monte Carlo on the marginal routes, with the NUTS
transition as its move.

`ctsem_sample_marginal` runs four chains of NUTS from a normal around the
estimate. Where the posterior's bulk sits far from that normal -- a time series
whose drift and diffusion parameters lie along curved ridges, with 41 of 75
directions three times wider than the Hessian says and its mass 95 normal sds
from the mode -- the chains take their whole budget to reach it and end at
R-hat 1.07-1.10; on a loading shell they end at 1.27. This samples the same
posterior by tempering instead: `nparticles` draws from the start normal
`q0 = N(values, widen^2 S0)` are carried to the posterior along
`(1 - beta) log q0 + beta log pi`, each stage's `beta` set so the incremental
weights keep an effective sample of `ess_frac` of the particles, then
systematic resampling, then `nmoves` NUTS transitions of every particle at the
stage's density. At `beta = 1` every particle is a chain at the posterior, and
the final phase draws `ndraws` from each and extends them under the one
stopping rule every sampler here stops by (`_sample_until_target`).

What it reuses, which is the point: the transition (`_nuts_transition!`
through `_continue_chain`), the pilot that sets a step size (`_run_chain`
with dual averaging), the metric (`_metric_from_covariances`), the chain
record (`_ChainResult`), the diagnostics (`ctsem_sample_diagnostics`) and the
stopping rule. It adds the tempering, the particles and three rules each of
which cost a run to learn (review/SAMPLING-smc-2026-10-23.md):

  * **The metric is the weighted cloud's covariance, shrunk toward the last
    stage's.** A raw covariance of 200 particles is singular below 200
    dimensions and noisy above it, and the particles then lag every stage's
    density: on a 100-dimensional normal with 150 particles the evidence came
    out 17 nats low. Shrunk toward its diagonal it threw away the correlations
    NUTS needs (trees from depth 4.9 to 6.4, twice the gradients on the time
    series). Whitened by the previous metric and shrunk toward the identity
    there by the Schaefer-Strimmer intensity, it keeps what was learned and
    costs what the first version cost.
  * **The step size is carried across stages**, a pilot of `pilot`
    dual-averaging iterations on one particle restarted from the last stage's
    step; found again from scratch each stage it jittered 0.05-0.08 at
    `beta` near 1.
  * **Particle `i`'s moves at stage `t` draw from `Xoshiro(hash((seed, t,
    i)))`**, so the result does not depend on how the work is laid out.

Moving a stage's particles until the cloud stopped drifting was measured and
does not pay once the metric is right (two blocks added in 24 stages on the
time series, for nothing), so a stage is one block of `nmoves`.

Particles run one after another on the shared objective, each gradient
parallel over the subject loop, as the marginal sampler's chains do and for
the same reason (`ctsem_evaluate`'s workspaces are per chunk). Spreading
particles over objective copies gave 3.4x on six copies in the prototype and
is the next lever; worker processes, as the chains have, the one after.

The final phase's chains are many and short (200 of 20 or so), and split
R-hat over them measures their autocorrelation rather than any disagreement:
for converged chains of `n` draws it reads about `sqrt(1 + (tau - 1) / n)`,
1.12-1.27 at ten draws on the time series where four long chains read 1.008.
The particles' independence after resampling and moves is the between-chain
guarantee here, so the caller passes `rhat_target = Inf` and the effective
size decides; R-hat is still returned, for the record.

No evidence is reported. The stages' log mean weights (`stage_logz`) come
back for the engine tests, whose closed-form targets have a known constant,
because a lag shows there where the final moments hide it; nothing in ctsem
reads them, and nothing should.
"""

using LinearAlgebra
using Random

# The engine does not depend on Statistics; a mean is a sum.
_smc_mean(v) = sum(v) / length(v)

export ctsem_sample_smc

"""Systematic resampling: indices drawn with probabilities `w`."""
function _smc_systematic(w::Vector{Float64}, rng::AbstractRNG)
    n = length(w)
    c = cumsum(w ./ sum(w))
    u = (rand(rng) .+ (0:(n - 1))) ./ n
    idx = Vector{Int}(undef, n)
    j = 1
    @inbounds for i in 1:n
        while j < n && c[j] < u[i]
            j += 1
        end
        idx[i] = j
    end
    return idx
end

"""The effective sample size of log weights."""
function _smc_ess(logw::Vector{Float64})
    m = maximum(logw)
    isfinite(m) || return 0.0
    w = exp.(logw .- m)
    return sum(w)^2 / sum(abs2, w)
end

"""
    _smc_shrunk_cov(X, w, Cprev)

The weighted covariance of X's columns (weights `w`, normalised here), shrunk
toward `Cprev`: the data are whitened by `Cprev`'s factor, their covariance
there is shrunk toward the identity scaled to its trace by the
Schaefer-Strimmer intensity -- the entries' estimated sampling variances over
their squared distances from that target, clamped to [0, 1] -- and the result
is mapped back. Zero where the cloud is large beside the dimension, near one
where it is small or has not moved, so whatever `Cprev` knew of the
correlations is kept wherever the particles say nothing new. Returns the
intensity too.
"""
function _smc_shrunk_cov(X::AbstractMatrix{Float64}, w::Vector{Float64},
    Cprev::AbstractMatrix{Float64})
    d, n = size(X)
    w = w ./ sum(w)
    neff = 1 / sum(abs2, w)
    m = X * w
    Lt = cholesky(Symmetric(Matrix{Float64}((Cprev .+ transpose(Cprev)) ./ 2))).L
    Z = Lt \ (X .- m)
    S = (Z .* transpose(w)) * transpose(Z)
    S = (S .+ transpose(S)) ./ 2
    mu = tr(S) / d
    scale = neff > 1 ? 1 / (neff - 1) : 1.0
    num = 0.0
    den = 0.0
    @inbounds for i in 1:d, j in i:d
        s = S[i, j]
        v = 0.0
        for k in 1:n
            v += w[k] * (Z[i, k] * Z[j, k] - s)^2
        end
        f = i == j ? 1.0 : 2.0
        num += f * v * scale
        den += f * (s - (i == j ? mu : 0.0))^2
    end
    lambda = den > 0 ? clamp(num / den, 0.0, 1.0) : 1.0
    C = (1 - lambda) .* S
    @inbounds for i in 1:d
        C[i, i] += lambda * mu
    end
    C = Lt * C * transpose(Lt)
    return (cov=(C .+ transpose(C)) ./ 2, lambda=lambda)
end

"""
    ctsem_sample_smc(objective, values; kwargs...)

Sample the marginal posterior of `objective` by tempered SMC from the normal
`N(values, widen^2 S0)`, `S0` the bounded inverse of `-hessian` (the metric
the marginal sampler would use; the exact Hessian at `values` when none is
given). See the file comment for the method. Returns what
`ctsem_sample_marginal` returns, the final phase's `nparticles` chains of
`ndraws` (extended to the targets) as the draws, with `sampler = "smc"`, the
stage record and the tempering's seconds beside.
"""
function ctsem_sample_smc(objective, values::AbstractVector;
    hessian::Union{Nothing,AbstractMatrix}=nothing, widen::Real=2.0,
    nparticles::Integer=200, nmoves::Integer=3, maxdepth::Integer=7,
    target_accept::Real=0.8, maxdelta::Real=1000.0, seed::Integer=20260828,
    ndraws::Integer=20, min_ess::Real=0.0, mean_ess::Real=0.0, max_draws::Integer=0,
    rhat_target::Real=Inf, ess_frac::Real=0.5, maxstages::Integer=500,
    pilot::Integer=30, gradient_method=:adjoint, verbose::Bool=false,
    progress_overwrite::Bool=true, progress_callback=nothing, progress_sink=nothing)

    t0 = time()
    forward0 = ctsem_forward_gradients()
    explosive0 = ctsem_explosive_passes()
    np = Int(nparticles)
    nmoves = Int(nmoves)
    maxdepth = Int(maxdepth)
    ndraws = Int(ndraws)
    np >= 2 || throw(ArgumentError("nparticles must be at least 2"))
    nmoves >= 1 || throw(ArgumentError("nmoves must be positive"))
    ndraws >= 1 || throw(ArgumentError("ndraws must be positive"))
    0 < target_accept < 1 || throw(ArgumentError("target_accept must be in (0, 1)"))
    0 < ess_frac < 1 || throw(ArgumentError("ess_frac must be in (0, 1)"))
    widen > 0 || throw(ArgumentError("widen must be positive"))

    m0 = collect(Float64, values)
    ndim = length(m0)
    npar = ndim
    maxdelta = Float64(maxdelta)
    target_accept = Float64(target_accept)
    callback = CTSEMCallback(progress_callback)

    # The log density and its gradient, refusing a point the route cannot use
    # as the marginal sampler does.
    logdensity! = function (g, x)
        result = try
            ctsem_evaluate(objective, x; gradient=true, gradient_method=gradient_method)
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
    values_at = X -> ctsem_evaluate_batch(objective, X)

    # The start: the metric's covariance, widened.
    H = hessian === nothing ? ctsem_hessian(objective, m0) : Matrix(hessian)
    information = Symmetric((-(H .+ transpose(H))) ./ 2)
    S0 = _bounded_inverse(information) .* Float64(widen)^2
    S0 = (S0 .+ transpose(S0)) ./ 2
    F0 = cholesky(Symmetric(S0))
    logdet0 = 2 * sum(log, diag(F0.L))
    P0 = inv(F0)
    logq = x -> begin
        z = F0.L \ (x .- m0)
        -0.5 * dot(z, z) - 0.5 * logdet0 - 0.5 * ndim * log(2pi)
    end
    gradq = x -> -(P0 * (x .- m0))
    tempered_at = function (b)
        return function (g, x)
            v = logdensity!(g, x)
            isfinite(v) || return -Inf
            qx = logq(x)
            g .= b .* g .+ (1 - b) .* gradq(x)
            return b * v + (1 - b) * qx
        end
    end

    verbose && println(_console(), "Sampling by tempered SMC: ", np, " particles, ",
        npar, " dimensions (effects integrated out), ", nmoves,
        " NUTS moves a stage, from the metric's normal at ", widen, "x")

    rng = Random.Xoshiro(UInt64(seed))
    X = m0 .+ F0.L * randn(rng, ndim, np)
    lp = values_at(X)
    lq = [logq(view(X, :, i)) for i in 1:np]
    beta = 0.0
    step = 0.0
    Cprev = copy(S0)
    metric = _metric_from_covariances([1:ndim], [copy(S0)])
    stages = NamedTuple[]
    ws = _NUTSWorkspace(ndim, maxdepth)
    g = zeros(ndim)
    for stage in 1:Int(maxstages)
        _ctsem_interrupt_check()
        u = lp .- lq
        u[.!isfinite.(u)] .= -Inf
        target = Float64(ess_frac) * np
        delta = 1.0 - beta
        if _smc_ess(delta .* u) < target
            lo, hi = 0.0, delta
            for _ in 1:60
                mid = (lo + hi) / 2
                _smc_ess(mid .* u) >= target ? (lo = mid) : (hi = mid)
            end
            delta = max(lo, 1e-12)
        end
        lw = delta .* u
        mx = maximum(lw)
        w = exp.(lw .- mx)
        inc = mx + log(_smc_mean(w))
        beta = min(1.0, beta + delta)
        shrunk = _smc_shrunk_cov(X, w, Cprev)
        Cprev = shrunk.cov
        metric = _metric_from_covariances([1:ndim], [shrunk.cov])
        idx = _smc_systematic(w, rng)
        X = X[:, idx]; lp = lp[idx]; lq = lq[idx]
        td! = tempered_at(beta)
        # The step: one particle's dual averaging from the carried step, twice
        # as long on the first stage, which has none. A failed pilot keeps the
        # carried step.
        npilot = step > 0 ? Int(pilot) : 2 * Int(pilot)
        piloted = try
            _run_chain(td!, X[:, 1], metric, Random.Xoshiro(hash((UInt64(seed), stage, 0))),
                npilot, 1, maxdepth, target_accept, maxdelta, 0.0, false, nothing;
                init_eps=step)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        if piloted !== nothing && isfinite(piloted.stepsize) && piloted.stepsize > 0
            X[:, 1] .= piloted.x
            step = piloted.stepsize
        elseif step <= 0
            step = 0.1
        end
        acc = 0.0
        depth = 0.0
        divergent = 0
        moved = 0
        for i in 1:np
            _ctsem_interrupt_check()
            x = X[:, i]
            logp = td!(g, x)
            isfinite(logp) || continue
            r = _continue_chain(td!, ws, Random.Xoshiro(hash((UInt64(seed), stage, i))),
                x, g, logp, step, metric, nmoves, maxdepth, maxdelta)
            moved += sum(abs2, r.x .- X[:, i]) > 0
            X[:, i] .= r.x
            acc += _smc_mean(r.accept)
            depth += _smc_mean(r.depth)
            divergent += count(r.divergent)
        end
        lp = values_at(X)
        lq = [logq(view(X, :, i)) for i in 1:np]
        push!(stages, (beta=beta, step=step, accept=acc / np, depth=depth / np,
            moved=moved / np, divergent=divergent, shrink=shrunk.lambda, logz=inc,
            secs=time() - t0))
        verbose && println(_console(), "  stage ", stage, ": beta ", round(beta; sigdigits=4),
            ", step ", round(step; sigdigits=3), ", accept ", round(acc / np; digits=2),
            ", depth ", round(depth / np; digits=1), ", ", round(time() - t0; digits=0), " s")
        _invoke_callback(callback, "tempering", round(Int, 100 * beta), 100,
            _smc_mean(filter(isfinite, lp)), divergent; force=beta >= 1.0)
        beta >= 1.0 && break
    end
    beta >= 1.0 || @warn "tempered SMC did not reach the posterior in $(maxstages) stages; the draws are from beta = $(beta)"
    smc_secs = time() - t0

    # The final phase: every particle a chain at the posterior, the step and
    # metric the last stage left, drawn and extended by the one stopping rule.
    results = Vector{_ChainResult}(undef, np)
    for i in 1:np
        _ctsem_interrupt_check()
        x = X[:, i]
        logp = logdensity!(g, x)
        if !isfinite(logp)
            # A particle whose density is not finite here cannot have moved to
            # it; it sits where the last finite move left it, which is X.
            logp = logdensity!(g, X[:, mod1(i + 1, np)])
            x = X[:, mod1(i + 1, np)]
        end
        results[i] = _continue_chain(logdensity!, ws, Random.Xoshiro(hash((UInt64(seed), 0, i))),
            x, g, logp, step, metric, ndraws, maxdepth, maxdelta)
    end
    held = Ref(results)
    extend! = function (wanted, attempt)
        held[] = _merge_chains(held[],
            _continue_chains(np, false, seed + 1000 * attempt, wanted, maxdepth, maxdelta,
                _ -> logdensity!, held[]; progress=false))
        return nothing
    end
    run = _sample_until_target(() -> _pool_draws(held[], npar), extend!, np, ndraws,
        Float64(min_ess), Float64(mean_ess), max(Int(max_draws), ndraws),
        Float64(rhat_target), verbose; t0=t0)
    results = held[]
    ndraws = run.total

    total = np * ndraws
    draws = Matrix{Float64}(undef, npar, total)
    accept = Vector{Float64}(undef, total)
    divergent = Vector{Bool}(undef, total)
    depth = Vector{Int}(undef, total)
    energy = Vector{Float64}(undef, total)
    for c in 1:np
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
    diagnostics = ctsem_sample_diagnostics(draws, np)
    return (
        draws=draws, npar=npar, ndim=npar, nchains=np, ndraws=ndraws,
        saved_effects=false, effect_mean=Float64[], effect_sd=Float64[],
        effect_cov=Float64[], effect_cov_rows=Int[], effect_cov_cols=Int[],
        rhat=diagnostics.rhat, ess=diagnostics.ess, ess_tail=diagnostics.ess_tail,
        accept=accept, divergent=divergent, depth=depth, energy=energy,
        stepsize=fill(step, np), warmup_divergent=zeros(Int, np),
        ndivergent=count(divergent), max_depth=maxdepth,
        nsaturated=count(==(maxdepth), depth),
        ebfmi=_ebfmi(energy, np),
        worst_rhat=_finite_extremum(diagnostics.rhat, maximum),
        min_ess=_finite_extremum(vcat(diagnostics.ess, diagnostics.ess_tail), minimum),
        target_trace=run.trace,
        forward_gradients=ctsem_forward_gradients() - forward0,
        explosive_passes=ctsem_explosive_passes() - explosive0,
        sampler="smc", nparticles=np, nmoves=nmoves, widen=Float64(widen),
        nstages=length(stages),
        stage_beta=[s.beta for s in stages], stage_step=[s.step for s in stages],
        stage_accept=[s.accept for s in stages], stage_depth=[s.depth for s in stages],
        stage_moved=[s.moved for s in stages], stage_divergent=[s.divergent for s in stages],
        stage_shrink=[s.shrink for s in stages], stage_logz=[s.logz for s in stages],
        stage_secs=[s.secs for s in stages], smc_secs=smc_secs,
    )
end
