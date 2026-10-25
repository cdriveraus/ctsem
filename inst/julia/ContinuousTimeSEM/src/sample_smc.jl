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
the final phase runs as many of them as its effective-size target needs,
draws `ndraws` from each and extends them under the one stopping rule every
sampler here stops by (`_sample_until_target`).

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
  * **The step size is carried across stages**, set by a pilot of
    `2 pilot` dual-averaging iterations on one particle at the first stage
    and by every particle's moves after it (`_smc_next_step`): the moves'
    mean acceptance at one stage sets the step of the next. Found again from
    scratch each stage it jittered 0.05-0.08 at `beta` near 1; a pilot each
    stage ran one particle while the others waited.
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
Since the correction is exact, a coreset's error costs only efficiency, and a
stage's coreset is the size of least expected cost per unit of tempering
(`_smc_stage_price`): its share of the moves, the correction's values, and
the stage's increment, which shrinks as the correction spends the stage's
effective size. The correction's own effective size (`stage_correction_ess`)
is the check: the error a stage's correction shows, against what its fit
expected, sizes the next coreset (`inflate`); a coreset whose error alone
drops it under the stage's target costs a stage that moves `beta` by nothing
and only corrects; two such stages in a row and the run goes on without a
coreset. The final phase's chains run on a coreset fitted to the posterior's
cloud, priced per effective draw (`_smc_draw_price`), and run until the
targets divided by the share of effective size the correction is expected to
keep;
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

**The final phase runs as many chains as its target needs**, at most one a
particle. It starts as many as the last stage's moves suggest
(`_smc_final_chains`): across the cloud, each parameter's correlation between
where its moves began and where they ended is its autocorrelation at that
lag, extended to one draw as an AR(1) decay (`_smc_move_tau`). Where their
effective size falls short it starts more, from the particles not yet used,
and only once every particle runs does it extend the chains. An unused
particle is an independent start, worth many draws of a chain that mixes
slowly, and the moves' short lag does not see a slow tail: on the time
series, where long NUTS chains take about 44 draws an effective draw, the
lag-3 estimate said 3.6, and 45 chains extended to 89 draws reached an
effective size of 87 where 200 chains of 20 had reached 536 (dev1). The
chains are taken in an order spread over the cloud (`_smc_spread_order`) and
share the particles' draw budget. Run from every particle on the loading
shells, the final phase took 4,000-5,000 draws for a target of 200 and
reached a minimum effective size of 1,188-2,088, most of a run's time.

The final phase's chains are short (20 draws or so), and split
R-hat over them measures their autocorrelation rather than any disagreement:
for converged chains of `n` draws it reads about `sqrt(1 + (tau - 1) / n)`,
1.12-1.27 at ten draws on the time series where four long chains read 1.008.
The particles' independence after resampling and moves is the between-chain
guarantee here, so the caller passes `rhat_target = Inf` and the effective
size decides. R-hat is returned as NaN: chains that start from particles
already drawn from the posterior agree whether or not they mixed, so the
number measures nothing else, and beside other samplers' R-hat in a summary
its 1.07-1.53 on runs that matched long NUTS read as a failure to converge.

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

"""Folds the particles are split into to hold each out once when a coreset
is chosen (`_smc_coreset_fit`)."""
const _SMC_CORESET_FOLDS = 4

"""A gradient's cost in full values, as the engine counts evaluations: on the
Laplace route both are mostly the units' inner solves."""
const _SMC_GRADIENT_VALUES = 2.0

"""
    _smc_leapfrogs(depths)

The leapfrog steps NUTS transitions of tree depths `depths` take, on
average: `2^depth - 1` a transition, which a tree stopped by a U-turn inside
its last doubling takes fewer than. Before any transition has run, a tree of
depth three.
"""
_smc_leapfrogs(depths) = isempty(depths) ? 7.0 : _smc_mean(2.0 .^ depths .- 1)

"""
    _smc_stage_price(share, s, work, beta, budget)

A tempering stage's cost a particle, in full evaluations, per unit of the
tempering it does, on a coreset whose units cost `share` of a full
evaluation and whose error differs across a particle's moves by `s` (an sd,
log units): the moves' `work` at the coreset's share (`nmoves` transitions of
their leapfrog steps, each a gradient), one full value at the moves' end and
the coreset's at both ends for the backward kernel's weight. That weight is
`beta` times the error's difference, and the stage's effective size allows
log weights a variance `budget`, `ln(1 / ess_frac)` for log-normal weights:
the correction spends `beta^2 s^2` of it and the tempering increment, which
the stages needed go as one over, the rest. `Inf` where the error alone
spends it. Without a coreset a stage costs `work + 1`.
"""
function _smc_stage_price(share::Real, s::Real, work::Real, beta::Real, budget::Real)
    left = 1 - (beta * s)^2 / budget
    left > 0 || return Inf
    return (work * share + 1 + 2 * share) / sqrt(left)
end

"""
    _smc_draw_price(share, s, work)

The final phase's cost an effective draw, in full evaluations, on a coreset
whose units cost `share` of a full evaluation and whose error at a draw
differs from its chain's start by `s`: the transition's `work` at the
coreset's share, a full value and a coreset value for the draw's weight, over
the share of effective size the weights keep, `exp(-s^2)` for log-normal
weights. Without a coreset a draw costs `work`.
"""
_smc_draw_price(share::Real, s::Real, work::Real) = (work * share + 1 + share) * exp(s^2)

"""
What a coreset costs a tempering stage or the final phase
(`_smc_stage_price`, `_smc_draw_price`), as a function of its share and the
spread of its error's difference, and what the same work costs without one
(`_smc_full_price`). `inflate` is how much more variance the last stage's
correction showed than its coreset's fit expected: the error's spread over
the cloud is measured where the particles are, and moves on the coreset take
them where it can be larger."""
struct _SMCPrice
    work::Float64
    beta::Float64
    budget::Float64
    inflate::Float64
    final::Bool
end

(p::_SMCPrice)(share::Real, s::Real) = p.final ?
    _smc_draw_price(share, s * sqrt(p.inflate), p.work) :
    _smc_stage_price(share, s * sqrt(p.inflate), p.work, p.beta, p.budget)
_smc_full_price(p::_SMCPrice) = p.final ? p.work : p.work + 1

"""
    _smc_coreset_fit(terms, total, cost, price, full, rng)

A coreset for the cloud whose units' terms are the columns of `terms`
(`nunits x nparticles`) and whose full weighted sum is `total`, the size of
least `price(share, s)` below `full`, what the work costs without one.
`share` is the units' cost (`cost`: each unit's share of a full evaluation)
and `s` the spread of the coreset's error's difference between two of the
cloud's points, `sqrt(2)` times its held-out spread over the cloud: the
weights that correct a coreset read its error at both ends of a move or a
chain, and a sampler's moves decorrelate the error between them nearly
fully (lag-three correlations of 0.1-0.5 on the race's fixtures). Held out
means: the particles are split into `_SMC_CORESET_FOLDS` folds at random
(`rng`), and for each a selection path on the others (`_coreset_cv`, on the
centred terms, refitted by NNLS each step) is read on the fold, as a
standard deviation over its particles, since a constant offset of the
coreset's value cancels in the normalised weights that read it. The units
are chosen out of sample as their weights are; the coreset is then fitted
on every particle at the size chosen. `nothing` where no size is cheaper
than none, or the usable particles are too few to hold any out.

Returns `(units, weights, error, price, share)`: `error` the held-out spread
over the cloud, `price` the size's; the weights replace the units' sampling
weights, since the total they were fitted to carries those.
"""
function _smc_coreset_fit(terms::AbstractMatrix{Float64}, total::AbstractVector{Float64},
    cost::AbstractVector{Float64}, price, full::Real, rng::AbstractRNG)
    N, n = size(terms)
    cols = Int[j for j in 1:n if isfinite(total[j]) && all(isfinite, view(terms, :, j))]
    m = length(cols)
    m >= 2 * _SMC_CORESET_FOLDS || return nothing
    A = Matrix{Float64}(transpose(view(terms, :, cols)))
    b = Float64[total[j] for j in cols]
    splits = _coreset_folds(m, _SMC_CORESET_FOLDS, rng)
    fullcost = Float64(full)
    judge = (paths, held) -> begin
        sq = 0.0
        count_ = 0
        share = Inf
        for k in eachindex(paths)
            P = paths[k]
            units = P.chosen
            e = A[held[k], units] * _coreset_weights(P) .- b[held[k]]
            e .-= _smc_mean(e)
            sq += sum(abs2, e)
            count_ += length(e)
            share = min(share, sum(cost[U] for U in units))
        end
        sigma = sqrt(sq / count_)
        c = price(share, sqrt(2.0) * sigma)
        # Any larger coreset costs at least its share without error.
        bound = fullcost - price(share, 0.0)
        ok = isfinite(c) && c < fullcost
        (ok ? fullcost - c : 0.0, bound, (error=sigma, price=c))
    end
    found = _coreset_cv(A, b, cost, Inf, _ctsem_coreset_sizes(N), splits, judge;
        refit=true, resid_tol=1e-10, centre=true)
    found === nothing && return nothing
    # Centred: the fit is to the variation over the cloud, the level is free.
    Ac = A .- sum(A; dims=1) ./ m
    bc = b .- _smc_mean(b)
    final = _coreset_select(Ac, bc, cost, Inf, [found.size], (u, w) -> (1.0, 0.0, nothing);
        refit=true, resid_tol=1e-10)
    final === nothing && return nothing
    keep = findall(>(0), final.weights)
    isempty(keep) && return nothing
    kept = final.units[keep]
    return (units=kept, weights=final.weights[keep], error=found.extra.error,
        price=found.extra.price, share=sum(cost[U] for U in kept))
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
    _smc_inflation(owed, X, Xpre, beta, error, previous)

How much more variance the stage's correction `owed` showed, over the
particles that moved, than its coreset's fit expected, `beta^2` times twice
its held-out spread `error` squared: what the next stage's coreset is priced
with. `previous` where too few moved to say.
"""
function _smc_inflation(owed::AbstractVector{Float64}, X::AbstractMatrix{Float64},
    Xpre::AbstractMatrix{Float64}, beta::Real, error::Real, previous::Real)
    v = Float64[owed[i] for i in eachindex(owed)
        if isfinite(owed[i]) && view(X, :, i) != view(Xpre, :, i)]
    expected = 2 * (beta * error)^2
    (length(v) >= 20 && expected > 0 && isfinite(expected)) || return Float64(previous)
    mu = _smc_mean(v)
    return sum(abs2, v .- mu) / (length(v) - 1) / expected
end

"""
    _smc_move_tau(Xfrom, Xto, nmoves)

The draws an effective draw takes for a chain moved by the last stage's
kernel, read off the cloud. Where the particles sit at the posterior, each
parameter's correlation across them between where they began the stage's
`nmoves` moves (`Xfrom`) and where they ended (`Xto`) is its autocorrelation
at that lag. An AR(1) decay from it, `r = rho^(1 / nmoves)` a draw, gives the
integrated autocorrelation time `(1 + r) / (1 - r)`. Returns the largest over
the parameters, and at least 1: a negative correlation is not credited.
"""
function _smc_move_tau(Xfrom::AbstractMatrix{Float64}, Xto::AbstractMatrix{Float64},
    nmoves::Int)
    tau = 1.0
    for j in axes(Xto, 1)
        a = view(Xfrom, j, :)
        b = view(Xto, j, :)
        (std(a) > 0 && std(b) > 0) || continue
        rho = cor(a, b)
        (isfinite(rho) && rho > 0) || continue
        r = min(rho^(1 / nmoves), 0.999)
        tau = max(tau, (1 + r) / (1 - r))
    end
    return tau
end

"""The fewest chains the final phase runs, where it has more particles: the
stopping rule reads the effective size across chains, and a chain of 20 draws
says little about its own autocorrelation."""
const _SMC_FINAL_MIN_CHAINS = 16

"""How many of `np` particles the final phase starts as chains: enough that
their first `ndraws` draws each reach 1.25 times the effective-size `target`
at `tau` draws an effective draw, and at least `_SMC_FINAL_MIN_CHAINS`, at most
one a particle. Every particle where there is no target."""
function _smc_final_chains(tau::Float64, target::Float64, ndraws::Int, np::Int)
    (target > 0 && isfinite(tau)) || return np
    n = ceil(Int, min(1.25 * tau * target / ndraws, Float64(np)))
    return clamp(n, min(np, _SMC_FINAL_MIN_CHAINS), np)
end

"""The particle indices `1:np` in an order whose every prefix is spread
evenly over the cloud: the van der Corput sequence's points, each taken to
the index it falls in, first come. A systematic resample leaves an
ancestor's copies side by side, so a spread prefix takes fewer copies of one
ancestor than the first indices would, and chains started later from the
order's next indices fill the gaps."""
function _smc_spread_order(np::Int)
    order = Int[]
    seen = falses(np)
    k = 0
    while length(order) < np
        v, f, n = 0.0, 0.5, k
        while n > 0
            v += f * (n & 1)
            n >>= 1
            f /= 2
        end
        i = 1 + floor(Int, v * np)
        if !seen[i]
            seen[i] = true
            push!(order, i)
        end
        k += 1
    end
    return order
end

"""
    _smc_next_step(step, accept, target)

The step for the next stage's moves from this stage's: their mean acceptance
`accept` at `step`. Where a leapfrog trajectory's energy error is normal with
sd `sigma`, its mean acceptance is `2 Phi(-sigma / 2)`, and `sigma` grows as
the step's square (Beskos, Pillai, Roberts, Sanz-Serna and Stuart 2013); so
the step for `target` is `step sqrt(sigma(target) / sigma(accept))`, held
within a factor of two a stage. Every particle's moves inform it, where a
pilot would read one particle's.
"""
function _smc_next_step(step::Float64, accept::Float64, target::Float64)
    isfinite(accept) || return step
    accept >= 1 && return 2 * step
    accept <= 0 && return step / 2
    sig = a -> -2 * sqrt(2) * erfinv(a - 1)
    return step * clamp(sqrt(sig(target) / sig(accept)), 0.5, 2.0)
end

"""
    _smc_logdensity(objective, coreset, gradient_method)

The log density and its gradient as a `(g, x)` closure, writing the gradient
into `g` and returning the value, `-Inf` where the route cannot use the
point: the full objective, or the coreset's density where one is given
(`_smc_coreset_eval`).
"""
function _smc_logdensity(objective, coreset, gradient_method)
    coreset === nothing &&
        return (g, x) -> _marginal_density!(objective, gradient_method, g, x)
    return (g, x) -> _coreset_density!(objective, coreset, g, x)
end

"""
    _smc_final_weights(objective, cs, owed, lp, Xf, start, draws, nper)

The final chains' log weights to the full density, a draw each (`draws`,
chain-major, `nper` a chain). Each chain's base is what the last stage's
moves left owing at its start (`owed`, by particle), less the final
coreset's error there, since the reversal of a chain on the coreset reads
its error at both ends; each draw adds the error where it sits. Without a
final coreset, the base alone. A systematic resample by these weights
returns the draws unweighted.
"""
function _smc_final_weights(objective, cs, owed::AbstractVector{Float64},
    lp::AbstractVector{Float64}, Xf::AbstractMatrix{Float64}, start::AbstractVector{Int},
    draws::AbstractMatrix{Float64}, nper::Int)
    nch = length(start)
    total = size(draws, 2)
    lc0 = cs === nothing ? nothing : _smc_coreset_values(objective, Xf[:, start], cs)
    base = [cs === nothing ? owed[start[ch]] :
        owed[start[ch]] - (lp[start[ch]] - lc0[ch]) for ch in 1:nch]
    lw = Vector{Float64}(undef, total)
    if cs === nothing
        for j in 1:total
            lw[j] = base[div(j - 1, nper) + 1]
        end
        return lw
    end
    lpf = ctsem_evaluate_batch(objective, draws)
    lpk = _smc_coreset_values(objective, draws, cs)
    for j in 1:total
        v = base[div(j - 1, nper) + 1] + (lpf[j] - lpk[j])
        lw[j] = isfinite(v) ? v : -Inf
    end
    return lw
end

"""The share of the draws' effective size that log weights `lw` keep."""
function _smc_kept(lw::AbstractVector{Float64})
    nfin = count(isfinite, lw)
    return nfin > 0 ? _importance_ess(lw) / nfin : 0.0
end

"""
    ctsem_sample_smc(objective, values; kwargs...)

Sample the marginal posterior of `objective` by tempered SMC from the normal
`N(values, widen^2 S0)`, `S0` the bounded inverse of `-hessian` (the metric
the marginal sampler would use; the exact Hessian at `values` when none is
given). See the file comment for the method, the coreset (`coreset`) and the
particles' parallelism (`parallel`). Returns what
`ctsem_sample_marginal` returns, the final phase's chains of `ndraws`
(extended to the targets) as the draws, with `sampler = "smc"`, the stage
record, the coreset record and the tempering's seconds beside. The final
phase runs as many chains as `min_ess` and `mean_ess` need (`final_tau`, the
draws an effective draw takes, says how many), at most `nparticles`; they
share the particles' budget, `max_draws` a particle.
"""
function ctsem_sample_smc(objective, values::AbstractVector;
    hessian::Union{Nothing,AbstractMatrix}=nothing, widen::Real=2.0,
    nparticles::Integer=200, nmoves::Integer=3, maxdepth::Integer=7,
    target_accept::Real=0.8, maxdelta::Real=1000.0, seed::Integer=20260828,
    ndraws::Integer=20, min_ess::Real=0.0, mean_ess::Real=0.0, max_draws::Integer=0,
    rhat_target::Real=Inf, ess_frac::Real=0.5, maxstages::Integer=500,
    pilot::Integer=30, gradient_method=:adjoint, coreset::Bool=true,
    parallel::Bool=true, verbose::Bool=false,
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
    # The log-weight variance a stage's effective size allows, and how much
    # more the last stage's correction showed than its coreset's fit expected.
    vbudget = log(1 / Float64(ess_frac))
    inflate = 1.0
    leap = _smc_leapfrogs(Int[])
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
    # Where the last stage's moves began, which the final phase reads its
    # chain count from (`_smc_move_tau`).
    Xfrom = X
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
            # is priced with the error it showed (`inflate`).
            delta = 0.0
            corrected = true
            corrections += 1
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
        # new); none where no size is cheaper than moving in full.
        fit = nothing
        if use_coreset && full.terms !== nothing
            price = _SMCPrice(nmv * leap * _SMC_GRADIENT_VALUES, beta, vbudget, inflate, false)
            fit = _smc_coreset_fit(full.terms, full.total, unit_cost, price,
                _smc_full_price(price), fit_rng)
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
        # The step: on the first stage, which has none, one particle's dual
        # averaging; after it, what the last stage's moves set
        # (`_smc_next_step`). A failed pilot starts from 0.1.
        Xstart = copy(Xpre)
        tpilot = time()
        if step <= 0
            piloted = try
                _run_chain(td!, Xpre[:, 1], metric,
                    Random.Xoshiro(hash((UInt64(seed), stage, 0))), 2 * Int(pilot), 1,
                    mdepth, tacc, mdelta, 0.0, false, nothing; init_eps=step)
            catch err
                err isa InterruptException && rethrow()
                nothing
            end
            if piloted !== nothing && isfinite(piloted.stepsize) && piloted.stepsize > 0
                Xstart[:, 1] .= piloted.x
                step = piloted.stepsize
            else
                step = 0.1
            end
        end
        pilot_secs = time() - tpilot
        Xnew = copy(Xstart)
        acc = zeros(np)
        dep = zeros(np)
        lf = zeros(np)
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
            lf[i] = _smc_leapfrogs(r.depth)
            ndiv[i] = count(r.divergent)
            mv[i] = r.x != xi ? 1 : 0
            okp[i] = true
            return nothing
        end
        X = Xnew
        Xfrom = Xstart
        nok = max(1, count(okp))
        count(okp) > 0 && (leap = sum(lf) / nok)
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
            inflate = _smc_inflation(owed, X, Xpre, beta, fit.error, inflate)
        end
        lq = _smc_columns(logq, X)
        push!(stages, (beta=beta, step=step, accept=sum(acc) / nok, depth=sum(dep) / nok,
            moved=sum(mv) / np, divergent=sum(ndiv), shrink=shrunk.lambda, logz=inc,
            coreset_size=cs_stage === nothing ? 0 : length(cs_stage.units),
            coreset_share=fit === nothing ? 1.0 : fit.share,
            coreset_error=fit === nothing ? NaN : fit.error,
            correction_ess=corr_ess / np, corrected=corrected,
            pilot_secs=pilot_secs, secs=time() - t0))
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
        # The next stage's step, and the final phase's after the last.
        step = _smc_next_step(step, sum(acc) / nok, tacc)
        beta >= 1.0 && break
    end
    _progress_done(reporter, "$(length(stages)) stages")
    beta >= 1.0 || @warn "tempered SMC did not reach the posterior in $(maxstages) stages; the draws are from beta = $(beta)"
    smc_secs = time() - t0

    # The final phase: particles run on as chains at the posterior, with the
    # step and metric the last stage left, until the one stopping rule is met;
    # on the coreset fitted to the posterior's cloud where one holds out, their
    # draws corrected to the full density afterwards, with what the last
    # stage's moves left owing. As many chains as the target needs, more
    # chains before longer ones (see the file comment), sharing the particles'
    # draw budget.
    tau = _smc_move_tau(Xfrom, X, nmv)
    order = _smc_spread_order(np)
    final_fit = nothing
    if use_coreset && full.terms !== nothing
        price = _SMCPrice(leap * _SMC_GRADIENT_VALUES, 1.0, vbudget, inflate, true)
        final_fit = _smc_coreset_fit(full.terms, full.total, unit_cost, price,
            _smc_full_price(price), fit_rng)
    end
    cs_final = final_fit === nothing ? nothing : CTSEMCoreset(final_fit.units, final_fit.weights)
    # The share of the chains' effective size the correction to the full
    # density is expected to keep, which the targets are divided by: the
    # chains' diagnostics are their own, before it.
    expected_kept = final_fit === nothing ? 1.0 : exp(-2 * inflate * final_fit.error^2)
    fmin_ess = Float64(min_ess) / expected_kept
    fmean_ess = Float64(mean_ess) / expected_kept
    target_ess = max(fmin_ess, fmean_ess)
    Xf = X
    lp_f = lp
    owed_f = copy(owed)
    # Each chain's particle, and the chains.
    start = Int[]
    held = Ref(_ChainResult[])
    final_step = step
    final_metric = metric
    start_chains! = function (ks)
        off = length(start)
        res = Vector{_ChainResult}(undef, length(ks))
        from = copy(ks)
        _each_point!(objective, length(ks), mode; coreset=cs_final,
            scratch=nuts_scratch) do i, obj, s
            ld = _smc_logdensity(obj, cs_final, gradient_method)
            k = ks[i]
            logp = ld(s.g, Xf[:, k])
            if !isfinite(logp)
                # A particle whose density is not finite here cannot have moved
                # to it; its chain starts from the next particle instead.
                k = mod1(k + 1, np)
                logp = ld(s.g, Xf[:, k])
            end
            from[i] = k
            res[i] = _continue_chain(ld, s.ws,
                Random.Xoshiro(hash((UInt64(seed), 0, off + i))), Xf[:, k], s.g, logp,
                final_step, final_metric, nfirst, mdepth, mdelta)
            return nothing
        end
        append!(start, from)
        held[] = vcat(held[], res)
        return nothing
    end
    start_chains!(order[1:_smc_final_chains(tau, target_ess, nfirst, np)])
    while target_ess > 0 && length(start) < np
        n = length(start)
        v = _sample_verdict(_pool_draws(held[], npar), fill(nfirst, n), nfirst, nfirst,
            fmin_ess, fmean_ess, nfirst, Float64(rhat_target), false)
        v.met && break
        shortfall = max(fmin_ess / max(v.worst, 1.0), fmean_ess / max(v.average, 1.0))
        start_chains!(order[(n + 1):min(np, max(n + 1, ceil(Int, 1.25 * shortfall * n)))])
    end
    nch = length(start)
    max_final = cld(max(Int(max_draws), nfirst) * np, nch)
    extend! = function (wanted, attempt)
        previous = held[]
        more = Vector{_ChainResult}(undef, length(previous))
        _each_point!(objective, length(previous), mode; coreset=cs_final,
            scratch=nuts_scratch) do i, obj, s
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
    run = _sample_until_target(() -> _pool_draws(held[], npar), extend!, nch, nfirst,
        fmin_ess, fmean_ess, max_final, Float64(rhat_target), verbose; t0=t0)
    nper = run.total
    trace = run.trace
    reweigh = cs_final !== nothing || any(!iszero, owed_f)
    pooled = _sample_pool(held[], nch, nper, npar)
    lw = reweigh ? _smc_final_weights(objective, cs_final, owed_f, lp_f, Xf, start,
        pooled.draws, nper) : Float64[]
    final_correction = reweigh ? _smc_kept(lw) : 1.0
    # The draws the correction keeps are what the targets ask for. Where it
    # kept less than the chains were run for -- the coreset's error is
    # measured over the cloud, and the chains on it go where it can be larger
    # -- they run on, once, to the targets over the share it kept, and the
    # correction is taken again over every draw.
    if cs_final !== nothing && !isempty(trace.ess) &&
        trace.ess[end] * final_correction < max(Float64(min_ess), Float64(mean_ess)) &&
        final_correction < expected_kept && nper < max_final
        keep_share = max(final_correction, 0.01)
        extra = _sample_until_target(() -> _pool_draws(held[], npar), extend!, nch, nper,
            Float64(min_ess) / keep_share, Float64(mean_ess) / keep_share, max_final, Float64(rhat_target),
            verbose; t0=t0)
        if extra.total > nper
            nper = extra.total
            trace = (draws=vcat(trace.draws, extra.trace.draws),
                secs=vcat(trace.secs, extra.trace.secs), ess=vcat(trace.ess, extra.trace.ess),
                rhat=vcat(trace.rhat, extra.trace.rhat))
            pooled = _sample_pool(held[], nch, nper, npar)
            lw = _smc_final_weights(objective, cs_final, owed_f, lp_f, Xf, start,
                pooled.draws, nper)
            final_correction = _smc_kept(lw)
        end
    end

    draws, accept, divergent = pooled.draws, pooled.accept, pooled.divergent
    depth, energy = pooled.depth, pooled.energy
    # The chains' own diagnostics, before any correction; no R-hat (see the
    # file comment).
    diagnostics = ctsem_sample_diagnostics(draws, nch)
    diagnostics = merge(diagnostics, (rhat=fill(NaN, length(diagnostics.rhat)),))
    if reweigh
        nfin = count(isfinite, lw)
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
        final_correction < 0.5 * expected_kept && @warn "the final chains' coreset kept $(round(100 * final_correction; digits=0))% of the draws' effective size after the correction to the full density, against $(round(100 * expected_kept; digits=0))% expected; coreset = false would give more"
    end
    corrected_pool = (draws=draws, accept=accept, divergent=divergent, depth=depth,
        energy=energy)
    return merge(_sample_result(corrected_pool; nchains=nch, ndraws=nper, npar=npar,
        ndim=npar, maxdepth=mdepth, stepsize=fill(step, nch),
        warmup_divergent=zeros(Int, nch), trace=trace, diagnostics=diagnostics), (
        forward_gradients=ctsem_forward_gradients() - forward0,
        explosive_passes=ctsem_explosive_passes() - explosive0,
        sampler="smc", nparticles=np, nmoves=nmv, widen=Float64(widen),
        nstages=length(stages),
        stage_beta=[s.beta for s in stages], stage_step=[s.step for s in stages],
        stage_accept=[s.accept for s in stages], stage_depth=[s.depth for s in stages],
        stage_moved=[s.moved for s in stages], stage_divergent=[s.divergent for s in stages],
        stage_shrink=[s.shrink for s in stages], stage_logz=[s.logz for s in stages],
        stage_secs=[s.secs for s in stages], stage_pilot_secs=[s.pilot_secs for s in stages],
        smc_secs=smc_secs,
        coreset=use_coreset || any(s.coreset_size > 0 for s in stages),
        nunits=nunits,
        stage_coreset_size=[s.coreset_size for s in stages],
        stage_coreset_share=[s.coreset_share for s in stages],
        stage_coreset_error=[s.coreset_error for s in stages],
        stage_correction_ess=[s.correction_ess for s in stages],
        stage_corrected=[s.corrected for s in stages],
        final_tau=tau,
        final_coreset_size=cs_final === nothing ? 0 : length(cs_final.units),
        final_coreset_error=final_fit === nothing ? NaN : final_fit.error,
        final_kept_expected=expected_kept,
        final_correction_ess=final_correction,
        parallel=String(mode),
    ))
end
