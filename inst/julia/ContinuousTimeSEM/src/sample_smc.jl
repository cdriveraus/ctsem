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
stopping rule; on the Laplace route the coreset (`CTSEMCoreset`, its
non-negative weights by `_ctsem_nnls`, its evaluation on copies by
`_laplace_coreset_each`'s pattern) and the pool (`_laplace_parallel`). It adds
the tempering, the particles and the rules below, each of which cost a run to
learn (review/SAMPLING-smc-2026-10-23.md):

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
    i)))`**, so the result does not depend on how the work is laid out: one
    particle after another, or spread over objective copies.

Moving a stage's particles until the cloud stopped drifting was measured and
does not pay once the metric is right (two blocks added in 24 stages on the
time series, for nothing), so a stage is one block of `nmoves`.

**On the Laplace route the particles move on a coreset** (`coreset = true`,
the default): a few units with non-negative weights whose weighted terms
follow the full sum over the current cloud (`_smc_coreset_fit`, fitted from
every unit's term at every particle, which the stage's weights need anyway).
A move on the coreset density costs the coreset's share of a gradient. The
moves are invariant for the tempered coreset density `q_t = q0^(1 - beta)
c^beta`, not for the stage's `pi_t = q0^(1 - beta) pi^beta`, and the
particles enter them distributed as `pi_t`, so after finitely many moves
they are neither. The correction is the SMC sampler's backward kernel (Del
Moral, Doucet and Jasra 2006) taken as the moves' own reversal under `q_t`:
a particle moved from `x` to `y` carries `r(y) / r(x)`, `r = (pi / c)^beta`,
whose log is `beta ((lp(y) - lc(y)) - (lp(x) - lc(x)))` -- the coreset's
error where the particle ended less where it began; one for a particle that
did not move (`_smc_move_correction`). The full values at `y` are the next
stage's weights' anyway, and the coreset's at `x` are one coreset batch
before the moves. Its expectation is one under the stage's density, so the
stages' log mean weights still sum to the log constant with a coreset as
without, and the engine tests check that on closed forms. An earlier version
took `exp(beta (lp(y) - lc(y)))` alone, which is right only once the moves
have mixed to `q_t`; with a few moves it tilted the cloud by the coreset's
error where the particles began (found by the chart session, 2026-10-23).
The correction's own effective size (`stage_correction_ess`) is the check: a
coreset whose error alone drops it under the stage's target is refitted at
half the tolerance, after a stage that moves `beta` by nothing and only
corrects; two such stages in a row and the run goes on without a coreset.
The final phase's chains run on a coreset fitted to the posterior's cloud;
the same reversal makes draw `k` of a chain that started at `x0` carry
`(pi / c)(x_k) / (pi / c)(x0)`, times whatever the last stage's moves left
owing at `x0`, and the draws are resampled by those weights
(`final_correction_ess` is the fraction kept); the chain diagnostics are the
chains' own, before the resampling. A coreset that follows the cloud poorly
makes the run slower, not its posterior different, which is what makes it
safe to run by default: `coreset = false` is the oracle.

**Particles run in parallel** (`parallel = true`), laid out as any set of
independent points is (`_each_point!`, parallel_points.jl): on the Laplace
pool with a copy of the objective a worker, or over deep copies of a filter
objective whose subjects are fewer than the pool (a time series is one), so
that a model with nothing to parallelise over inside a gradient gets the
whole ceiling across particles; one after another otherwise.

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

"""
    _importance_ess(logw)

The effective sample size of log weights, `(sum w)^2 / sum w^2`: 0 where none
is finite, and at most their number, which equal weights reach exactly but
the arithmetic can pass by a rounding error. SMC's stages and corrections and
DAIS's damping read it.
"""
function _importance_ess(logw::AbstractVector{Float64})
    m = maximum(logw)
    isfinite(m) || return 0.0
    w = exp.(logw .- m)
    return min(Float64(length(w)), sum(w)^2 / sum(abs2, w))
end

"""`f` of each column of `X`. A function of its own rather than a
comprehension where `X` is rebound: a comprehension is a closure, and one
over a local assigned twice boxes it (CLAUDE.md, "Failures that pass every
check short of loading")."""
_smc_columns(f, X::AbstractMatrix) = [f(view(X, :, i)) for i in axes(X, 2)]

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

# ------------------------------------------------------------------ the coreset

"""The largest share of a full evaluation a stage's coreset may cost. Above
it the moves save too little for the fit and the correction to be worth
having; the full values at the stage's end are paid either way."""
const _SMC_CORESET_BUDGET = 0.5

"""Folds the particles are split into to hold each out once when a coreset
is chosen (`_smc_coreset_fit`)."""
const _SMC_CORESET_FOLDS = 4

"""
    _smc_coreset_fit(terms, total, cost, tol, budget, rng)

A coreset for the cloud whose units' terms are the columns of `terms`
(`nunits x nparticles`) and whose full weighted sum is `total`: the units
chosen by forward selection on the centred terms (`_coreset_select`, refitted
by NNLS each step) -- each step the unit whose term, orthogonalised against
those in, best explains what is left of the centred total with a positive
weight -- their weights by non-negative least squares (`_ctsem_nnls`), at the
smallest size on the path whose held-out error is within `tol`. Held out means: the particles are split into
`_SMC_CORESET_FOLDS` folds at random (`rng`), the weights are fitted on the
other folds and the error read on the one left out, as a standard deviation
over its particles, since a constant offset of the coreset's value cancels in
the normalised weights that read it. The sizes tried are
`_ctsem_coreset_sizes`', 1.25 apart, within `budget` of a full evaluation
(`cost`: each unit's share). `nothing` where no size within the budget passes,
or the usable particles are too few to hold any out.

Returns `(units, weights, error, share)`; the weights replace the units'
sampling weights, since the total they were fitted to carries those.
"""
function _smc_coreset_fit(terms::AbstractMatrix{Float64}, total::AbstractVector{Float64},
    cost::AbstractVector{Float64}, tol::Real, budget::Real, rng::AbstractRNG)
    N, n = size(terms)
    cols = Int[j for j in 1:n if isfinite(total[j]) && all(isfinite, view(terms, :, j))]
    m = length(cols)
    m >= 2 * _SMC_CORESET_FOLDS || return nothing
    A = Matrix{Float64}(transpose(view(terms, :, cols)))
    b = Float64[total[j] for j in cols]
    # Centred: the fit is to the variation over the cloud, the level is free.
    @inbounds for U in 1:N
        mu = 0.0
        for i in 1:m
            mu += A[i, U]
        end
        mu /= m
        for i in 1:m
            A[i, U] -= mu
        end
    end
    b .-= _smc_mean(b)
    perm = randperm(rng, m)
    nf = _SMC_CORESET_FOLDS
    edges = round.(Int, range(0, m; length=nf + 1))
    # The first size whose held-out error passes, with a positive weight:
    # scored one there, after which no larger can score more.
    value = (units, w) -> begin
        worst = 0.0
        for f in 1:nf
            held = perm[(edges[f] + 1):edges[f + 1]]
            train = vcat(perm[1:edges[f]], perm[(edges[f + 1] + 1):m])
            At = A[train, units]
            wf = _ctsem_nnls(transpose(At) * At, transpose(At) * b[train])
            e = A[held, units] * wf .- b[held]
            e .-= _smc_mean(e)
            worst = max(worst, sqrt(sum(abs2, e) / length(e)))
        end
        passes = isfinite(worst) && worst <= tol && any(>(0), w)
        (passes ? 1.0 : 0.0, passes ? 1.0 : 2.0, worst)
    end
    found = _coreset_select(A, b, cost, budget, _ctsem_coreset_sizes(N), value;
        refit=true, resid_tol=1e-10)
    found === nothing && return nothing
    keep = findall(>(0), found.weights)
    kept = found.units[keep]
    return (units=kept, weights=found.weights[keep], error=found.extra,
        share=sum(cost[U] for U in kept))
end

# What the sampler asks of an objective whose log density is a sum of units'
# terms, so that its particles can move on a coreset of them. The Laplace
# route answers; any other objective is moved in full. A test can answer for
# a closed-form target of its own.

"""Whether `objective` is a sum of at least two units' terms a coreset can
stand in for."""
_smc_coreset_capable(objective) = false
_smc_coreset_capable(laplace::CTSEMLaplaceObjective) = length(laplace.units.members) >= 2

"""Each unit's share of a full evaluation."""
function _smc_unit_cost(laplace::CTSEMLaplaceObjective)
    w = _laplace_unit_weights(laplace)
    return w ./ sum(w)
end

"""Every column of `X` in full: the log density (`value`, `-Inf` where the
point cannot be used), each unit's term (`terms`, `nunits x n`), and their
weighted sum without the prior (`total`)."""
function _smc_unit_terms(laplace::CTSEMLaplaceObjective, X::AbstractMatrix{Float64})
    r = ctsem_coreset_terms(laplace, X)
    return (value=r.value, total=r.total, terms=r.unit_loglik)
end

"""The coreset's log density (prior included) at every column of `X`, `-Inf`
where unusable; the draws several at once on copies."""
function _smc_coreset_values(laplace::CTSEMLaplaceObjective, X::AbstractMatrix{Float64},
    coreset::CTSEMCoreset)
    results = _laplace_coreset_each(laplace, X, coreset, false)
    return Float64[r === nothing ? -Inf : r.value for r in results]
end

"""The coreset's log density and gradient at `x`, or `nothing` where the
point cannot be used."""
function _smc_coreset_eval(laplace::CTSEMLaplaceObjective, coreset::CTSEMCoreset,
    x::AbstractVector)
    r = _laplace_coreset_try(laplace, Vector{Float64}(x), true, coreset)
    return r === nothing ? nothing : (value=r.value, gradient=r.gradient)
end

"""
    _smc_terms(objective, X, want_terms)

Every column of `X` evaluated in full: the log density (`value`, `-Inf` where
the point cannot be used), and, where the objective is a sum of units' terms
and `want_terms`, each unit's term and their weighted sum without the prior,
for fitting a coreset (`_smc_unit_terms`).
"""
function _smc_terms(objective, X::AbstractMatrix{Float64}, want_terms::Bool)
    want_terms && _smc_coreset_capable(objective) && return _smc_unit_terms(objective, X)
    return (value=ctsem_evaluate_batch(objective, X), total=nothing, terms=nothing)
end

"""
    _smc_move_correction(beta, lp_y, lc_y, lp_x, lc_x, moved)

The log weight a move made on the coreset density at `beta` leaves owing, for
a particle that moved from `x` to `y`: `beta ((lp(y) - lc(y)) - (lp(x) -
lc(x)))`, `lp` the full log density and `lc` the coreset's -- the backward
kernel taken as the moves' reversal under the tempered coreset density (see
the file comment); zero for one that did not move; `-Inf` where a value at a
moved particle is not finite.
"""
function _smc_move_correction(beta::Real, lp_y::Real, lc_y::Real, lp_x::Real, lc_x::Real,
    moved::Bool)
    moved || return 0.0
    v = beta * ((lp_y - lc_y) - (lp_x - lc_x))
    return isfinite(v) ? Float64(v) : -Inf
end

"""
    _smc_logdensity(objective, coreset, gradient_method)

The log density and its gradient as a `(g, x)` closure, writing the gradient
into `g` and returning the value, `-Inf` where the route cannot use the
point: the full objective, or the coreset's density where one is given
(`_smc_coreset_eval`).
"""
function _smc_logdensity(objective, coreset, gradient_method)
    if coreset !== nothing
        return function (g, x)
            r = _smc_coreset_eval(objective, coreset, x)
            if r === nothing || !isfinite(r.value) || !all(isfinite, r.gradient)
                fill!(g, 0.0)
                return -Inf
            end
            copyto!(g, r.gradient)
            return r.value
        end
    end
    return function (g, x)
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
end

"""
    ctsem_sample_smc(objective, values; kwargs...)

Sample the marginal posterior of `objective` by tempered SMC from the normal
`N(values, widen^2 S0)`, `S0` the bounded inverse of `-hessian` (the metric
the marginal sampler would use; the exact Hessian at `values` when none is
given). See the file comment for the method, the coreset (`coreset`,
`coreset_tol`: the held-out error in log units a stage's coreset may have)
and the particles' parallelism (`parallel`). Returns what
`ctsem_sample_marginal` returns, the final phase's `nparticles` chains of
`ndraws` (extended to the targets) as the draws, with `sampler = "smc"`, the
stage record, the coreset record and the tempering's seconds beside.
"""
function ctsem_sample_smc(objective, values::AbstractVector;
    hessian::Union{Nothing,AbstractMatrix}=nothing, widen::Real=2.0,
    nparticles::Integer=200, nmoves::Integer=3, maxdepth::Integer=7,
    target_accept::Real=0.8, maxdelta::Real=1000.0, seed::Integer=20260828,
    ndraws::Integer=20, min_ess::Real=0.0, mean_ess::Real=0.0, max_draws::Integer=0,
    rhat_target::Real=Inf, ess_frac::Real=0.5, maxstages::Integer=500,
    pilot::Integer=30, gradient_method=:adjoint, coreset::Bool=true,
    coreset_tol::Real=0.3, parallel::Bool=true, verbose::Bool=false,
    progress_overwrite::Bool=true, progress_callback=nothing, progress_sink=nothing)

    t0 = time()
    forward0 = ctsem_forward_gradients()
    explosive0 = ctsem_explosive_passes()
    # Converted once into names of their own: a closure that captures a
    # variable assigned twice boxes it (CLAUDE.md, "No Core.Box").
    np = Int(nparticles)
    nmv = Int(nmoves)
    mdepth = Int(maxdepth)
    nfirst = Int(ndraws)
    np >= 2 || throw(ArgumentError("nparticles must be at least 2"))
    nmv >= 1 || throw(ArgumentError("nmoves must be positive"))
    _check_sample_arguments(np, nfirst, target_accept)
    0 < ess_frac < 1 || throw(ArgumentError("ess_frac must be in (0, 1)"))
    widen > 0 || throw(ArgumentError("widen must be positive"))
    coreset_tol > 0 || throw(ArgumentError("coreset_tol must be positive"))

    m0 = collect(Float64, values)
    ndim = length(m0)
    npar = ndim
    mdelta = Float64(maxdelta)
    tacc = Float64(target_accept)
    callback = CTSEMCallback(progress_callback)

    # The coreset applies where the objective is a sum of units' terms: the
    # Laplace route with at least two units (`_smc_coreset_capable`).
    capable = _smc_coreset_capable(objective)
    unit_cost = capable ? _smc_unit_cost(objective) : Float64[]
    nunits = capable ? length(unit_cost) : 1
    use_coreset = coreset && capable
    tol_now = Float64(coreset_tol)
    corrections = 0
    mode = _points_mode(objective, parallel)
    # Each worker's NUTS workspace and gradient buffer (`_each_point!`).
    nuts_scratch = () -> (ws=_NUTSWorkspace(ndim, mdepth), g=zeros(ndim))
    fit_rng = Random.Xoshiro(hash((UInt64(seed), :coreset)))

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
    # The tempered density at `b` on an objective (a copy, or the shared
    # one), on the coreset where one is given.
    tempered_for = function (obj, cs, b)
        ld = _smc_logdensity(obj, cs, gradient_method)
        return function (g, x)
            v = ld(g, x)
            isfinite(v) || return -Inf
            qx = logq(x)
            g .= b .* g .+ (1 - b) .* gradq(x)
            return b * v + (1 - b) * qx
        end
    end

    # The stages report as every other long call does (progress.jl): through
    # the caller's sink where it gave one (a front end, a knitr run), in place
    # on the console otherwise. Printed straight to the console they reached
    # neither.
    reporter = CTSEMProgress(verbose; label="tempering", overwrite=progress_overwrite,
        sink=progress_sink)
    if _due(reporter, true)
        _emit(reporter, string("  tempered SMC: ", np, " particles, ", npar,
            " dimensions, ", nmv, " NUTS moves a stage, from the metric's normal at ",
            widen, "x", use_coreset ? ", moves on a coreset of the $nunits units" : "",
            mode == :serial ? "" : ", particles in parallel on copies"))
    end

    rng = Random.Xoshiro(UInt64(seed))
    X = m0 .+ F0.L * randn(rng, ndim, np)
    full = _smc_terms(objective, X, use_coreset)
    lp = full.value
    lq = _smc_columns(logq, X)
    beta = 0.0
    step = 0.0
    Cprev = copy(S0)
    metric = _metric_from_covariances([1:ndim], [copy(S0)])
    stages = NamedTuple[]
    # What the last moves left owing, a log weight per particle
    # (`_smc_move_correction`): zero after moves on the full objective.
    owed = zeros(np)
    for stage in 1:Int(maxstages)
        _ctsem_interrupt_check()
        u = lp .- lq
        u[.!isfinite.(u)] .= -Inf
        target = Float64(ess_frac) * np
        corr_ess = all(iszero, owed) ? Float64(np) : _importance_ess(owed)
        delta = 1.0 - beta
        corrected = false
        if corr_ess < target
            # The coreset's error alone spends the stage's effective size: this
            # stage corrects it and moves `beta` by nothing, and the next fit
            # is held to half the tolerance.
            delta = 0.0
            corrected = true
            corrections += 1
            tol_now /= 2
        elseif _importance_ess(delta .* u .+ owed) < target
            lo, hi = 0.0, delta
            for _ in 1:60
                mid = (lo + hi) / 2
                _importance_ess(mid .* u .+ owed) >= target ? (lo = mid) : (hi = mid)
            end
            delta = max(lo, 1e-12)
        end
        if !corrected
            corrections = 0
        elseif corrections >= 2
            # Two correction stages in a row: the coreset is not following
            # this cloud, and the full objective moves the particles from here.
            use_coreset = false
        end
        lw = delta .* u .+ owed
        mx = maximum(lw)
        w = exp.(lw .- mx)
        inc = mx + log(_smc_mean(w))
        beta = min(1.0, beta + delta)
        shrunk = _smc_shrunk_cov(X, w, Cprev)
        Cprev = shrunk.cov
        metric = _metric_from_covariances([1:ndim], [shrunk.cov])
        # The coreset for this stage's moves, from every unit's term at the
        # cloud before it is resampled (the duplicates after it say nothing
        # new); none where no size within the budget holds out.
        fit = nothing
        if use_coreset && full.terms !== nothing
            fit = _smc_coreset_fit(full.terms, full.total, unit_cost, tol_now,
                _SMC_CORESET_BUDGET, fit_rng)
        end
        cs_stage = fit === nothing ? nothing : CTSEMCoreset(fit.units, fit.weights)
        idx = _smc_systematic(w, rng)
        Xpre = X[:, idx]
        lp_pre = lp[idx]
        lq = lq[idx]
        # Where each particle begins its moves, pilot included, and the
        # coreset's value there: the backward kernel's weight reads the
        # coreset's error at both ends of the move.
        lc_pre = cs_stage === nothing ? lp_pre : _smc_coreset_values(objective, Xpre, cs_stage)
        td! = tempered_for(objective, cs_stage, beta)
        # The step: one particle's dual averaging from the carried step, twice
        # as long on the first stage, which has none. A failed pilot keeps the
        # carried step.
        npilot = step > 0 ? Int(pilot) : 2 * Int(pilot)
        piloted = try
            _run_chain(td!, Xpre[:, 1], metric, Random.Xoshiro(hash((UInt64(seed), stage, 0))),
                npilot, 1, mdepth, tacc, mdelta, 0.0, false, nothing;
                init_eps=step)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        Xstart = copy(Xpre)
        if piloted !== nothing && isfinite(piloted.stepsize) && piloted.stepsize > 0
            Xstart[:, 1] .= piloted.x
            step = piloted.stepsize
        elseif step <= 0
            step = 0.1
        end
        Xnew = copy(Xstart)
        acc = zeros(np)
        dep = zeros(np)
        ndiv = zeros(Int, np)
        mv = zeros(Int, np)
        # A Vector{Bool}, not a BitVector: particles on other tasks write their
        # own entries at once, and a BitVector's entries share words.
        okp = zeros(Bool, np)
        stage_beta = beta
        stage_step = step
        stage_metric = metric
        _each_point!(objective, np, mode; coreset=cs_stage, scratch=nuts_scratch) do i, obj, s
            tdi! = tempered_for(obj, cs_stage, stage_beta)
            # The chain advances its start in place; the start is kept apart.
            xi = Xstart[:, i]
            logp = tdi!(s.g, xi)
            isfinite(logp) || return nothing
            r = _continue_chain(tdi!, s.ws, Random.Xoshiro(hash((UInt64(seed), stage, i))),
                copy(xi), s.g, logp, stage_step, stage_metric, nmv, mdepth, mdelta)
            Xnew[:, i] .= r.x
            acc[i] = _smc_mean(r.accept)
            dep[i] = _smc_mean(r.depth)
            ndiv[i] = count(r.divergent)
            mv[i] = r.x != xi ? 1 : 0
            okp[i] = true
            return nothing
        end
        X = Xnew
        nok = max(1, count(okp))
        full = _smc_terms(objective, X, use_coreset)
        lp = full.value
        if cs_stage === nothing
            fill!(owed, 0.0)
        else
            lc = _smc_coreset_values(objective, X, cs_stage)
            @inbounds for i in 1:np
                owed[i] = _smc_move_correction(beta, lp[i], lc[i], lp_pre[i], lc_pre[i],
                    view(X, :, i) != view(Xpre, :, i))
            end
        end
        lq = _smc_columns(logq, X)
        push!(stages, (beta=beta, step=step, accept=sum(acc) / nok, depth=sum(dep) / nok,
            moved=sum(mv) / np, divergent=sum(ndiv), shrink=shrunk.lambda, logz=inc,
            coreset_size=cs_stage === nothing ? 0 : length(cs_stage.units),
            coreset_share=fit === nothing ? 1.0 : fit.share,
            coreset_error=fit === nothing ? NaN : fit.error,
            correction_ess=corr_ess / np, corrected=corrected,
            secs=time() - t0))
        if _due(reporter, true)
            _emit(reporter, join(filter(!isempty, [
                @sprintf("  tempering stage %d", stage),
                @sprintf("beta %.4g", beta), @sprintf("step %.3g", step),
                @sprintf("accept %.2f", sum(acc) / nok), @sprintf("depth %.1f", sum(dep) / nok),
                cs_stage === nothing ? "" : "coreset $(length(cs_stage.units)) of $nunits units",
                corrected ? "correction only" : "",
                _duration(time() - t0)]), " | "))
        end
        _invoke_callback(callback, "tempering", round(Int, 100 * beta), 100,
            _smc_mean(filter(isfinite, lp)), sum(ndiv); force=beta >= 1.0)
        beta >= 1.0 && break
    end
    _progress_done(reporter, "$(length(stages)) stages")
    beta >= 1.0 || @warn "tempered SMC did not reach the posterior in $(maxstages) stages; the draws are from beta = $(beta)"
    smc_secs = time() - t0

    # The final phase: every particle a chain at the posterior, the step and
    # metric the last stage left, drawn and extended by the one stopping rule;
    # on the coreset fitted to the posterior's cloud where one holds out, its
    # draws corrected to the full density afterwards, with what the last
    # stage's moves left owing.
    final_fit = nothing
    if use_coreset && full.terms !== nothing
        final_fit = _smc_coreset_fit(full.terms, full.total, unit_cost, tol_now,
            _SMC_CORESET_BUDGET, fit_rng)
    end
    cs_final = final_fit === nothing ? nothing : CTSEMCoreset(final_fit.units, final_fit.weights)
    Xf = X
    lp_f = lp
    owed_f = copy(owed)
    first_results = Vector{_ChainResult}(undef, np)
    start = collect(1:np)
    final_step = step
    final_metric = metric
    _each_point!(objective, np, mode; coreset=cs_final, scratch=nuts_scratch) do i, obj, s
        ld = _smc_logdensity(obj, cs_final, gradient_method)
        k = i
        logp = ld(s.g, Xf[:, k])
        if !isfinite(logp)
            # A particle whose density is not finite here cannot have moved to
            # it; its chain starts from the next particle instead.
            k = mod1(i + 1, np)
            logp = ld(s.g, Xf[:, k])
        end
        start[i] = k
        first_results[i] = _continue_chain(ld, s.ws, Random.Xoshiro(hash((UInt64(seed), 0, i))),
            Xf[:, k], s.g, logp, final_step, final_metric, nfirst, mdepth, mdelta)
        return nothing
    end
    held = Ref(first_results)
    extend! = function (wanted, attempt)
        previous = held[]
        more = Vector{_ChainResult}(undef, np)
        _each_point!(objective, np, mode; coreset=cs_final, scratch=nuts_scratch) do i, obj, s
            ld = _smc_logdensity(obj, cs_final, gradient_method)
            more[i] = _run_chain(ld, previous[i].x, previous[i].metric,
                Random.Xoshiro(UInt64(seed + 1000 * attempt) + UInt64(i)),
                0, wanted, mdepth, 0.8, mdelta, 0.0, false, nothing;
                resume=previous[i])
            return nothing
        end
        held[] = _merge_chains(previous, more)
        return nothing
    end
    run = _sample_until_target(() -> _pool_draws(held[], npar), extend!, np, nfirst,
        Float64(min_ess), Float64(mean_ess), max(Int(max_draws), nfirst),
        Float64(rhat_target), verbose; t0=t0)
    chains = held[]
    nper = run.total

    total = np * nper
    pooled = _sample_pool(chains, np, nper, npar)
    draws, accept, divergent = pooled.draws, pooled.accept, pooled.divergent
    depth, energy = pooled.depth, pooled.energy
    # The chains' own diagnostics, before any correction.
    diagnostics = ctsem_sample_diagnostics(draws, np)
    final_correction = 1.0
    if cs_final !== nothing || any(!iszero, owed_f)
        # Each chain's base: what the last stage's moves left owing at its
        # start, less the final coreset's error there (the reversal of a chain
        # on the coreset reads its error at both ends); each draw adds the
        # error where it sits. A systematic resample by these weights returns
        # the draws unweighted.
        base = zeros(np)
        lc0 = cs_final === nothing ? nothing : _smc_coreset_values(objective, Xf, cs_final)
        for ch in 1:np
            k = start[ch]
            base[ch] = cs_final === nothing ? owed_f[k] : owed_f[k] - (lp_f[k] - lc0[k])
        end
        lw = Vector{Float64}(undef, total)
        if cs_final === nothing
            for j in 1:total
                lw[j] = base[div(j - 1, nper) + 1]
            end
        else
            lpf = ctsem_evaluate_batch(objective, draws)
            lpk = _smc_coreset_values(objective, draws, cs_final)
            for j in 1:total
                v = base[div(j - 1, nper) + 1] + (lpf[j] - lpk[j])
                lw[j] = isfinite(v) ? v : -Inf
            end
        end
        nfin = count(isfinite, lw)
        final_correction = nfin > 0 ? _importance_ess(lw) / nfin : 0.0
        if nfin > 0
            mx = maximum(lw)
            wk = exp.(lw .- mx)
            idx = _smc_systematic(wk, Random.Xoshiro(hash((UInt64(seed), :final))))
            draws = draws[:, idx]
            accept = accept[idx]
            divergent = divergent[idx]
            depth = depth[idx]
            energy = energy[idx]
        end
        final_correction < 0.5 && @warn "the final chains' coreset kept $(round(100 * final_correction; digits=0))% of the draws' effective size after the correction to the full density; a larger coreset_tol budget or coreset = false would give more"
    end
    corrected_pool = (draws=draws, accept=accept, divergent=divergent, depth=depth,
        energy=energy)
    return merge(_sample_result(corrected_pool; nchains=np, ndraws=nper, npar=npar,
        ndim=npar, maxdepth=mdepth, stepsize=fill(step, np),
        warmup_divergent=zeros(Int, np), trace=run.trace, diagnostics=diagnostics), (
        forward_gradients=ctsem_forward_gradients() - forward0,
        explosive_passes=ctsem_explosive_passes() - explosive0,
        sampler="smc", nparticles=np, nmoves=nmv, widen=Float64(widen),
        nstages=length(stages),
        stage_beta=[s.beta for s in stages], stage_step=[s.step for s in stages],
        stage_accept=[s.accept for s in stages], stage_depth=[s.depth for s in stages],
        stage_moved=[s.moved for s in stages], stage_divergent=[s.divergent for s in stages],
        stage_shrink=[s.shrink for s in stages], stage_logz=[s.logz for s in stages],
        stage_secs=[s.secs for s in stages], smc_secs=smc_secs,
        coreset=use_coreset || any(s.coreset_size > 0 for s in stages),
        nunits=nunits,
        stage_coreset_size=[s.coreset_size for s in stages],
        stage_coreset_share=[s.coreset_share for s in stages],
        stage_coreset_error=[s.coreset_error for s in stages],
        stage_correction_ess=[s.correction_ess for s in stages],
        stage_corrected=[s.corrected for s in stages],
        final_coreset_size=cs_final === nothing ? 0 : length(cs_final.units),
        final_correction_ess=final_correction,
        parallel=String(mode),
    ))
end
