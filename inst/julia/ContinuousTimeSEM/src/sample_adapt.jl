"""
Warmup adaptation and convergence diagnostics for the sampler.

Separate from `sample_nuts.jl` because none of it is about Hamiltonian dynamics:
these are the schedules and estimators that decide *how* the sampler is tuned
and *whether* its output can be believed, and both are worth reading without the
trajectory machinery in the way.

Stan's schemes throughout, and deliberately so. They are the ones the applied
literature's rules of thumb are calibrated against -- "R-hat below 1.01", "ESS
above 400", "no divergences" mean what people expect them to mean only if the
estimators behind them are the familiar ones.
"""

using LinearAlgebra
using Random

################################################################################
# Step size
################################################################################

"""
    _init_stepsize(logdensity!, metric, rng, x, g, logp, ws)

A step size in the right order of magnitude, by doubling or halving from 1 until
one step's acceptance probability crosses one half.

Crude on purpose: dual averaging converges from anywhere reasonable, and its
cost is entirely in how far away "anywhere reasonable" was. Starting from a
fixed 1 on a model whose scales differ by orders of magnitude spends most of the
warmup walking back.
"""
function _init_stepsize(logdensity!, metric::CTSEMMetric, rng::AbstractRNG,
    x::Vector{Float64}, g::Vector{Float64}, logp::Float64, ws::_NUTSWorkspace)
    point = ws.working
    trial = function (e)
        copyto!(point.x, x); copyto!(point.g, g); point.logp = logp
        _metric_momentum!(point.p, metric, rng)
        h0 = -logp + _metric_velocity!(ws.velocity, metric, point.p, ws.scratch)
        kinetic = _leapfrog!(point, logdensity!, e, metric, ws.velocity, ws.scratch)
        (!isfinite(kinetic) || !isfinite(point.logp)) && return -Inf
        return h0 - (-point.logp + kinetic)
    end
    eps = 1.0
    direction = trial(eps) > log(0.5) ? 1 : -1
    for _ in 1:100
        ratio = trial(eps)
        if direction == 1
            ratio > log(0.5) || break
            eps *= 2
        else
            ratio < log(0.5) || break
            eps /= 2
        end
        (eps > 1e7 || eps < 1e-10) && break
    end
    return clamp(eps, 1e-10, 1e7)
end

"""Nesterov dual averaging on `log(eps)`, as Stan runs it."""
mutable struct _DualAverage
    mu::Float64
    target::Float64
    gamma::Float64
    t0::Float64
    kappa::Float64
    counter::Int
    hbar::Float64
    logeps::Float64
    logepsbar::Float64
end

# `logepsbar` starts at `log(eps)` rather than at zero, and that is a fix
# rather than a preference. `_dual_final` returns `exp(logepsbar)`, so a zero
# start means "step size 1" for as long as no update has been taken -- and one
# update is enough to erase it, since the first update's weight is
# `1^(-kappa) = 1`. The only way to reach sampling with zero updates is to
# restart on the last warmup iteration or to ask for no warmup at all, and both
# happen: see `_dual_restart!`.
_DualAverage(eps::Float64, target::Float64) =
    _DualAverage(log(10 * eps), target, 0.05, 10.0, 0.75, 0, 0.0, log(eps),
        log(eps))

"""Take one acceptance statistic and return the next step size to use."""
function _dual_update!(da::_DualAverage, accept::Float64)
    da.counter += 1
    eta = 1 / (da.counter + da.t0)
    da.hbar = (1 - eta) * da.hbar + eta * (da.target - accept)
    da.logeps = da.mu - sqrt(da.counter) / da.gamma * da.hbar
    weight = da.counter^(-da.kappa)
    da.logepsbar = weight * da.logeps + (1 - weight) * da.logepsbar
    return exp(da.logeps)
end

"""
Restart the averaging around a new centre, after the metric has changed.

`logepsbar = log(eps)`, not zero, for the reason given at the constructor: the
averaged step size is what sampling uses, and a restart that is never followed
by an update would otherwise hand sampling `exp(0) = 1`.

That is not hypothetical. The divergence-triggered acceptance raise in
`_run_chain` checks on a fixed stride of 50, so on a 500-iteration warmup it
fires *on iteration 500* -- and a raise there restarts the averaging with no
iteration left to update it. The chain then samples at step size 1 on a target whose adapted step
size is around 0.03, diverges on every transition, and returns 500 draws of a
stuck chain: one such chain in a four-chain run took R-hat to 4.11 and the
minimum effective size to 2.2. The same zero also explains why `nwarmup = 0`
appeared to be useless -- it was not that the step size went unadapted, it was
that `_init_stepsize`'s answer was discarded and replaced by 1.
"""
function _dual_restart!(da::_DualAverage, eps::Float64)
    da.mu = log(10 * eps)
    da.counter = 0
    da.hbar = 0.0
    da.logeps = log(eps)
    da.logepsbar = log(eps)
    return eps
end

"""The averaged step size, which is what sampling should use."""
_dual_final(da::_DualAverage) = exp(da.logepsbar)

################################################################################
# Metric adaptation schedule
################################################################################

"""
Warmup iterations excluded from the metric estimate, and from the divergence
rate that raises `target_accept`. Stan's `init_buffer`, and a constant here
because `_run_chain` has to agree with `_adapt_windows` about where the first
window starts -- it did not, see the note on `window_draws` there.
"""
const _ADAPT_INIT_BUFFER = 75

"""
    _adapt_windows(nwarmup; init_buffer, term_buffer, base_window)

The warmup iterations at which the metric is re-estimated.

Stan's schedule: an opening buffer where only the step size moves, then doubling
windows each ending in a metric update, then a closing buffer where the metric
is fixed and the step size settles against it. Both buffers earn their keep -- a
metric estimated under a badly wrong step size is itself wrong, and a step size
tuned against a metric that is about to change is thrown away.

Returns an empty schedule when warmup is too short to hold one window, in which
case the Laplace metric is simply used as it stands. That is a perfectly
reasonable mode of operation here and not a degenerate one: it is already a good
metric, which is the whole point of starting from a fit.
"""

function _adapt_windows(nwarmup::Int; init_buffer::Int=_ADAPT_INIT_BUFFER,
    term_buffer::Int=50, base_window::Int=25)
    nwarmup < init_buffer + term_buffer + base_window && return Int[]
    ends = Int[]
    start = init_buffer + 1
    window = base_window
    limit = nwarmup - term_buffer
    while start <= limit
        stop = start + window - 1
        # Absorb a next window that would not fit into this one, rather than
        # ending on a stub too short to estimate anything from.
        if stop + 2 * window > limit
            stop = limit
        end
        stop = min(stop, limit)
        push!(ends, stop)
        start = stop + 1
        window *= 2
    end
    return ends
end

"""
    _estimate_metric(draws, ranges)

Block covariances from warmup draws, shrunk toward the identity.

`(n / (n + 5)) S + 1e-3 (5 / (n + 5)) I` is Stan's regularisation, and it is
what lets a `k x k` block be estimated from fewer than `k` draws: the shrinkage
holds it positive definite until there are enough draws to say otherwise. That
matters more here than in a diagonal sampler, because the blocks are the point.
"""
function _estimate_metric(draws::Vector{Vector{Float64}}, base::CTSEMMetric;
    adapt::Union{Nothing,Vector{Bool}}=nothing)
    ranges = base.ranges
    n = length(draws)
    covariances = Vector{Matrix{Float64}}(undef, length(ranges))
    factors = Vector{Union{Nothing,Matrix{Float64}}}(nothing, length(ranges))
    for b in eachindex(ranges)
        r = ranges[b]
        k = length(r)
        # A block the caller declined to adapt keeps the factor it came in with.
        if adapt !== nothing && !adapt[b]
            factors[b] = base.factors[b]
            covariances[b] = Matrix(1.0I, k, k)
            continue
        end
        if n < 2
            covariances[b] = Matrix(1.0I, k, k)
            continue
        end
        mean = zeros(k)
        for d in draws
            for (j, i) in enumerate(r)
                mean[j] += d[i]
            end
        end
        mean ./= n
        S = zeros(k, k)
        centred = zeros(k)
        for d in draws
            for (j, i) in enumerate(r)
                centred[j] = d[i] - mean[j]
            end
            for a in 1:k, c in 1:k
                S[a, c] += centred[a] * centred[c]
            end
        end
        S ./= (n - 1)
        S .*= n / (n + 5.0)
        for i in 1:k
            S[i, i] += 1e-3 * (5.0 / (n + 5.0))
        end
        covariances[b] = S
    end
    estimated = _metric_from_covariances(ranges, covariances)
    adapt === nothing && return estimated
    kept = [factors[b] === nothing ? estimated.factors[b] : factors[b]
            for b in eachindex(ranges)]
    return CTSEMMetric(ranges, kept)
end

################################################################################
# Convergence diagnostics
################################################################################

"""Split R-hat: each chain is halved first, so within-chain drift is caught."""
function _split_rhat(chains::AbstractMatrix{Float64})
    ndraws, nchains = size(chains)
    ndraws < 4 && return NaN
    half = div(ndraws, 2)
    pieces = Vector{Vector{Float64}}()
    for c in 1:nchains
        push!(pieces, collect(view(chains, 1:half, c)))
        push!(pieces, collect(view(chains, (ndraws - half + 1):ndraws, c)))
    end
    m = length(pieces)
    means = [sum(p) / half for p in pieces]
    vars = [sum(abs2, p .- means[j]) / (half - 1) for (j, p) in enumerate(pieces)]
    W = sum(vars) / m
    W <= 0 && return NaN
    grand = sum(means) / m
    B = half * sum(abs2, means .- grand) / (m - 1)
    varplus = ((half - 1) / half) * W + B / half
    return sqrt(varplus / W)
end

"""
    _ess(chains)

Effective sample size by Geyer's initial monotone positive sequence: sum the
autocorrelations in pairs, stop at the first non-positive pair, then enforce
that the pair sums are non-increasing.

The pairing is not an optimisation. An autocorrelation estimate for a reversible
chain is positive in consecutive pairs even when individual lags are not, so
truncating on single lags stops early and overstates the sample size -- exactly
the direction of error a diagnostic must not have.

Each lag's autocorrelation is computed only when the sequence reaches it: the
sum stops at the first non-positive pair, usually within a few dozen lags, and
computing all thousand lags up front cost the same result many times over.
"""
function _ess(chains::AbstractMatrix{Float64})
    ndraws, nchains = size(chains)
    ndraws < 4 && return NaN
    means = [sum(view(chains, :, c)) / ndraws for c in 1:nchains]
    vars = [sum(abs2, view(chains, :, c) .- means[c]) / (ndraws - 1) for c in 1:nchains]
    W = sum(vars) / nchains
    W <= 0 && return NaN
    grand = sum(means) / nchains
    B = nchains > 1 ? ndraws * sum(abs2, means .- grand) / (nchains - 1) : 0.0
    varplus = ((ndraws - 1) / ndraws) * W + (nchains > 1 ? B / ndraws : 0.0)
    varplus <= 0 && return NaN

    maxlag = min(ndraws - 2, 1000)
    corr = function (t)
        acov = 0.0
        for c in 1:nchains
            s = 0.0
            m = means[c]
            @inbounds for i in 1:(ndraws - t)
                s += (chains[i, c] - m) * (chains[i + t, c] - m)
            end
            acov += s / ndraws
        end
        acov /= nchains
        return 1 - (W - acov) / varplus
    end

    # Initial positive sequence, then made monotone.
    pairs = Float64[]
    t = 1
    while t + 1 <= maxlag
        pair = corr(t) + corr(t + 1)
        pair <= 0 && break
        push!(pairs, pair)
        t += 2
    end
    for i in 2:length(pairs)
        pairs[i] = min(pairs[i], pairs[i - 1])
    end
    tau = max(1.0, -1.0 + 2 * sum(pairs))
    return nchains * ndraws / tau
end

"""
    _rank_normalise(x)

The draws replaced by normal scores of their ranks over all chains,
`Phi^-1((r - 3/8) / (S + 1/4))`, ties at their average rank (Vehtari, Gelman,
Simpson, Carpenter and Buerkner 2021). R-hat and effective size of these are
the bulk diagnostics: defined whatever the tails do, and the same under any
monotone transform of the parameter -- which matters here, where a raw scale
parameter and its constrained value are both reported.
"""
function _rank_normalise(x::AbstractMatrix{Float64})
    v = vec(x)
    S = length(v)
    order = sortperm(v)
    ranks = Vector{Float64}(undef, S)
    i = 1
    while i <= S
        j = i
        while j < S && v[order[j + 1]] == v[order[i]]
            j += 1
        end
        r = (i + j) / 2
        for k in i:j
            ranks[order[k]] = r
        end
        i = j + 1
    end
    z = similar(x)
    for k in 1:S
        z[k] = sqrt(2.0) * erfinv(2 * (ranks[k] - 3 / 8) / (S + 1 / 4) - 1)
    end
    return z
end

"""
    _tail_ess(x)

The smaller of the effective sizes of `x <= q05` and `x >= q95` as indicators,
the quantiles over all chains: how many effective draws the 5% and 95% points
rest on. A chain that has not visited a tail lately leaves bulk effective size
untouched and this one low, which is the failure it is for: a study-level
standard deviation reached its lower tail only in bursts thousands of draws
apart, and bulk size read about 160 where the tail held about 40
(review/POSTERIOR-race-2026-10-02.md, gN3).
"""
function _tail_ess(x::AbstractMatrix{Float64})
    v = sort(vec(x))
    S = length(v)
    lo = v[clamp(ceil(Int, 0.05 * S), 1, S)]
    hi = v[clamp(ceil(Int, 0.95 * S), 1, S)]
    a = _ess(Float64.(x .<= lo))
    b = _ess(Float64.(x .>= hi))
    return _finite_extremum([a, b], minimum)
end

export ctsem_sample_diagnostics
"""
    ctsem_sample_diagnostics(draws, nchains)

R-hat, bulk and tail effective sample size for every sampled coordinate, as
Stan reports them now (Vehtari et al. 2021). R-hat is the larger of the split
R-hats of the rank-normalised draws and of their rank-normalised distances from
the median: the first sees chains that disagree on location, the second chains
that agree on location but not on spread -- chains that have missed a tail.
`ess` is the bulk effective size (rank-normalised draws) and `ess_tail` the
tails' (`_tail_ess`).

`draws` is `ndim x (nchains * ndraws)`, chain-major -- the layout `ctsem_sample`
returns, and the one that crosses the R bridge as a single matrix.
"""
function ctsem_sample_diagnostics(draws::AbstractMatrix{Float64}, nchains::Integer)
    ndim, total = size(draws)
    nchains = Int(nchains)
    nchains >= 1 || throw(ArgumentError("nchains must be positive"))
    ndraws = div(total, nchains)
    ndraws * nchains == total ||
        throw(DimensionMismatch("draw count is not a multiple of the chain count"))
    rhat = fill(NaN, ndim)
    ess = fill(NaN, ndim)
    ess_tail = fill(NaN, ndim)
    buffer = Matrix{Float64}(undef, ndraws, nchains)
    for j in 1:ndim
        @inbounds for c in 1:nchains, t in 1:ndraws
            buffer[t, c] = draws[j, (c - 1) * ndraws + t]
        end
        all(isfinite, buffer) || continue
        z = _rank_normalise(buffer)
        sorted = sort(vec(buffer))
        n = length(sorted)
        med = isodd(n) ? sorted[(n + 1) ÷ 2] : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2
        folded = _rank_normalise(abs.(buffer .- med))
        rhat[j] = _finite_extremum([_split_rhat(z), _split_rhat(folded)], maximum)
        ess[j] = _ess(z)
        ess_tail[j] = _tail_ess(buffer)
    end
    return (rhat=rhat, ess=ess, ess_tail=ess_tail)
end

"""
    _ebfmi(energy, nchains)

Energy-Bayesian fraction of missing information, per chain.

Low values (below about 0.3) say the sampler is not exploring the energy
distribution -- the marginal energy has heavier tails than the transitions can
cross, which on a hierarchical model usually means a funnel the metric cannot
straighten. It is the one diagnostic that catches that failure when R-hat and
divergences do not.
"""
function _ebfmi(energy::Vector{Float64}, nchains::Integer)
    nchains = Int(nchains)
    ndraws = div(length(energy), nchains)
    out = fill(NaN, nchains)
    for c in 1:nchains
        e = view(energy, ((c - 1) * ndraws + 1):(c * ndraws))
        length(e) < 3 && continue
        mean = sum(e) / length(e)
        variance = sum(abs2, e .- mean) / (length(e) - 1)
        variance <= 0 && continue
        delta = 0.0
        for i in 2:length(e)
            delta += (e[i] - e[i - 1])^2
        end
        out[c] = (delta / (length(e) - 1)) / variance
    end
    return out
end
