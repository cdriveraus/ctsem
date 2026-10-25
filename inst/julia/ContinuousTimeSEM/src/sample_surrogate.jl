# NUTS on a coreset, each draw's endpoint checked against the full density.
#
# On the Laplace route a gradient costs every unit's inner solve, and a NUTS
# transition takes a gradient at every leapfrog step. Where a few units with
# weights follow the sum of all of them over the posterior -- the coreset SMC
# moves its particles on (`_smc_coreset_fit`) -- the trajectory can run on
# that instead, and the full density is needed only where the trajectory
# ends. The transition is then delayed acceptance (Christen and Fox 2005;
# Liu's surrogate transition, Monte Carlo Strategies in Scientific Computing
# 2001, section 9.4.3): from `x`, a NUTS transition targeting the coreset's
# density `c` proposes `y`, accepted with probability
#
#     min(1, pi(y) c(x) / (pi(x) c(y))),
#
# `pi` the full density, and otherwise the chain stays at `x`. A NUTS
# transition at a fixed step size and metric is reversible with respect to
# the density it targets, and the composite is then reversible with respect
# to `pi`, whatever the coreset: a poor one costs acceptance, not accuracy.
# `coreset = false` is the oracle.
#
# A chain fits its coreset twice, from the full density's unit terms at
# `_SURROGATE_TRAIN` states: before warm-up, at draws of the normal it starts
# in -- the fit's estimate and the metric's covariance (`_surrogate_seed`) --
# and a third of the way through it (`_surrogate_train_at`), at its own
# states so far, which by then cover the posterior. Warm-up runs on the
# coreset from the start, with the endpoint test, so the step size is tuned to
# the transitions that will be sampled with. Run in full, warm-up was most of
# a short run's cost: on a 400-subject model whose coreset was a fortieth of
# the units, its 200 transitions took 12 of the run's 15 s, and its first
# third still most of the 11 s once the rest ran on the coreset. Warm-up
# draws are discarded, so the test there only keeps the chain where the
# posterior is, and a normal that fits the posterior poorly costs acceptance
# until the second fit. Warm-up's end and each batch of draws whose endpoints
# were accepted too rarely to pay for the coreset give it up for the rest of
# the run: one change, after which the kernel is fixed again. A chain resumed
# in a new call has no warm-up to fit from and samples in full.
#
# What one draw costs: the trajectory's leapfrog steps at the coreset's share
# of a gradient, and one full value, no gradient, at a moved endpoint; a
# refused endpoint wastes the draw. So the coreset is the size of least cost
# per accepted transition (`_surrogate_price`), and is given up where the
# acceptance measured makes that more than a transition in full.

"""Warm-up states a chain's coreset is fitted from: its last this many."""
const _SURROGATE_TRAIN = 64

"""Moved endpoints a chain's acceptance is read over, to price its coreset
again or give it up."""
const _SURROGATE_MIN_PROPOSED = 20

"""
    _surrogate_price(share, s, work)

A transition's cost, in full evaluations, per accepted transition, on a
coreset whose units cost `share` of a full evaluation and whose error differs
between a transition's ends by `s` (an sd, log units): the trajectory's
`work` (its leapfrog steps, a gradient each) at the coreset's share and one
full value at its end, over the endpoint test's acceptance, about
`2 Phi(-s / 2)` for a normal error. Without a coreset a transition costs
`work`.
"""
_surrogate_price(share::Real, s::Real, work::Real) =
    (work * share + 1) / (2 * _ctsem_normal_cdf(-s / 2))

"""What a coreset costs a chain's transition (`_surrogate_price`), with the
error's spread scaled by `inflate`, the variance the last coreset's measured
acceptance implies over what its fit expected."""
struct _SurrogatePrice
    work::Float64
    inflate::Float64
end

(p::_SurrogatePrice)(share::Real, s::Real) = _surrogate_price(share, s * sqrt(p.inflate), p.work)

"""
    _surrogate_work(d)

The leapfrog steps a transition of chain `d` takes on its coreset, each a
gradient in full values (`_SMC_GRADIENT_VALUES`), on average over the
transitions so far; before any, a tree of depth three's.
"""
function _surrogate_work(d)
    leap = d.transitions >= 10 ? d.leapfrogs / d.transitions : _smc_leapfrogs(Int[])
    return _SMC_GRADIENT_VALUES * leap
end

"""
    _surrogate_inflation(d)

How much more variance in its error's difference chain `d`'s coreset showed
than its fit expected: the spread an acceptance `a` implies, `-2
Phi^-1(a / 2)`, squared, over twice its held-out spread squared. One where
fewer than `_SURROGATE_MIN_PROPOSED` endpoints were tested, or every one was
accepted.
"""
function _surrogate_inflation(d)
    (d.coreset === nothing || d.proposed < _SURROGATE_MIN_PROPOSED ||
        !(d.error > 0)) && return 1.0
    a = d.accepted / d.proposed
    0 < a < 1 || return a <= 0 ? 100.0 : 1.0
    implied = 2 * sqrt(2.0) * erfcinv(a)
    return implied^2 / (2 * d.error^2)
end

"""
    _SurrogateDensity(objective, gradient_method, ndim, seed)

One chain's log density on the marginal route: the full objective until it
has a coreset (`_surrogate_fit!`), the coreset after, with the state the
endpoint test needs. Called as the chain's `logdensity!(g, x)`.
"""
mutable struct _SurrogateDensity{O}
    objective::O
    gradient_method::Symbol
    coreset::Union{Nothing,CTSEMCoreset}
    # The full and the coreset log density at the chain's current point.
    lp::Float64
    lc::Float64
    xold::Vector{Float64}
    gold::Vector{Float64}
    proposed::Int
    accepted::Int
    batch_proposed::Int
    batch_accepted::Int
    # Coreset gradients taken, and transitions run on the coreset.
    leapfrogs::Int
    transitions::Int
    size::Int
    share::Float64
    error::Float64
    fits::Int
    abandoned::Bool
    seed::UInt64
end

_SurrogateDensity(objective, gradient_method, ndim::Integer, seed::Integer) =
    _SurrogateDensity(objective, Symbol(gradient_method), nothing, NaN, NaN,
        zeros(Int(ndim)), zeros(Int(ndim)), 0, 0, 0, 0, 0, 0, 0, 0.0, NaN, 0, false,
        UInt64(seed))

"""The marginal samplers' full log density and gradient at `x`, into `g`;
`-Inf` and a zero gradient where the point cannot be used. NUTS
(`ctsem_sample_marginal`), SMC (`_smc_logdensity`) and a chain's density
before its coreset all evaluate through it."""
function _marginal_density!(objective, gradient_method, g::AbstractVector{Float64},
    x::AbstractVector)
    result = try
        ctsem_evaluate(objective, x; gradient=true, gradient_method=gradient_method)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    if result === nothing || !isfinite(result.value) ||
        result.gradient === nothing || !all(isfinite, result.gradient)
        fill!(g, 0.0)
        return -Inf
    end
    copyto!(g, result.gradient)
    return Float64(result.value)
end

"""The coreset's log density and gradient at `x`, into `g`; `-Inf` and a
zero gradient where the point cannot be used (`_smc_coreset_eval`)."""
function _coreset_density!(objective, coreset::CTSEMCoreset, g::AbstractVector{Float64},
    x::AbstractVector)
    r = _smc_coreset_eval(objective, coreset, x)
    if r === nothing || !isfinite(r.value) || !all(isfinite, r.gradient)
        fill!(g, 0.0)
        return -Inf
    end
    copyto!(g, r.gradient)
    return Float64(r.value)
end

"""The full log density at `x` without its gradient, `-Inf` where unusable."""
function _marginal_value(objective, gradient_method, x::AbstractVector)
    result = try
        ctsem_evaluate(objective, x; gradient=false, gradient_method=gradient_method)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    (result === nothing || !isfinite(result.value)) && return -Inf
    return Float64(result.value)
end

function (d::_SurrogateDensity)(g::AbstractVector{Float64}, x::AbstractVector)
    cs = d.coreset
    cs === nothing && return _marginal_density!(d.objective, d.gradient_method, g, x)
    d.leapfrogs += 1
    return _coreset_density!(d.objective, cs, g, x)
end

# The hooks `_run_chain` and `_sample_draws!` call. Any other density answers
# with nothing to do, so a plain closure samples exactly as before.

"""The warm-up iteration a chain fits its coreset at: a third of the way in
where that leaves `_SURROGATE_TRAIN` states before it, the end otherwise."""
_surrogate_train_at(nwarmup::Integer) =
    nwarmup >= 3 * _SURROGATE_TRAIN ? Int(nwarmup) ÷ 3 : Int(nwarmup)

"""Whether a chain's density still fits a coreset: before warm-up and once
in it (`_SURROGATE_FITS`), and not after giving one up."""
_surrogate_training(logdensity!) = false
_surrogate_training(d::_SurrogateDensity) = d.fits < _SURROGATE_FITS && !d.abandoned &&
    _smc_coreset_capable(d.objective)

"""The fits a chain makes: one from draws of the normal it starts in, one
from its own warm-up states."""
const _SURROGATE_FITS = 2

"""
    _surrogate_seed(logdensity!, centre, metric, n)

`n` draws from the normal a chain starts in, `centre` and the metric's
covariance, to fit a first coreset from before warm-up. Drawn on the
density's own stream, so the chain's draws use the stream they would
without a coreset.
"""
_surrogate_seed(logdensity!, centre, metric, n) = Vector{Float64}[]
function _surrogate_seed(d::_SurrogateDensity, centre::Vector{Float64},
    metric::CTSEMMetric, n::Int)
    rng = Random.Xoshiro(d.seed ⊻ 0x9e3779b97f4a7c15)
    states = Vector{Vector{Float64}}(undef, n)
    for s in 1:n
        x = copy(centre)
        for b in eachindex(metric.ranges)
            r = metric.ranges[b]
            C = metric.factors[b]
            z = randn(rng, length(r))
            for i in eachindex(r)
                acc = 0.0
                for j in 1:i
                    acc += C[i, j] * z[j]
                end
                x[r[i]] = centre[r[i]] + acc
            end
        end
        states[s] = x
    end
    return states
end

"""Keep `x` among the last `_SURROGATE_TRAIN` warm-up states."""
function _surrogate_keep!(states::Vector{Vector{Float64}}, x::Vector{Float64})
    length(states) >= _SURROGATE_TRAIN && popfirst!(states)
    push!(states, copy(x))
    return states
end

"""
    _surrogate_fit!(logdensity!, states, x, g, logp)

Fit a chain's coreset from the full density's unit terms at `states`, and
switch the density to it if one passes; a chain that has one keeps it where
none does. `logp` is the density in use at `x` (the full one before any
coreset); the full value there is kept for the endpoint test. Returns the
log density in use at `x` afterwards, its gradient in `g`.
"""
_surrogate_fit!(logdensity!, states, x, g, logp) = logp
function _surrogate_fit!(d::_SurrogateDensity, states::Vector{Vector{Float64}},
    x::Vector{Float64}, g::Vector{Float64}, logp::Float64)
    d.fits += 1
    (_smc_coreset_capable(d.objective) && isfinite(logp) &&
        length(states) >= 2 * _SMC_CORESET_FOLDS) || return logp
    full = try
        _smc_unit_terms(d.objective, reduce(hcat, states))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    full === nothing && return logp
    price = _SurrogatePrice(_surrogate_work(d), _surrogate_inflation(d))
    fit = _smc_coreset_fit(full.terms, full.total, _smc_unit_cost(d.objective), price,
        price.work, Random.Xoshiro(d.seed + UInt64(d.fits)))
    fit === nothing && return logp
    previous = d.coreset
    d.coreset = CTSEMCoreset(fit.units, fit.weights)
    gold = copy(g)
    lc = d(g, x)
    if !isfinite(lc)
        d.coreset = previous
        copyto!(g, gold)
        return logp
    end
    previous === nothing && (d.lp = logp)
    d.lc = lc
    d.size = length(fit.units)
    d.share = fit.share
    d.error = fit.error
    # Acceptance is read for the coreset in use: warm-up's end gives up only
    # one whose own moves were refused too often.
    d.batch_proposed = 0
    d.batch_accepted = 0
    return lc
end

"""Before a transition: remember where it starts."""
_surrogate_begin!(logdensity!, x, g) = nothing
function _surrogate_begin!(d::_SurrogateDensity, x::Vector{Float64}, g::Vector{Float64})
    d.coreset === nothing && return nothing
    copyto!(d.xold, x)
    copyto!(d.gold, g)
    return nothing
end

"""
    _surrogate_correct!(logdensity!, x, g, logp, rng)

After a transition on the coreset from `d.xold` to `x` (log density `logp`
there): accept it with probability `min(1, pi(x) c(xold) / (pi(xold) c(x)))`,
or put the chain back at `xold`. Returns the log density at the chain's point,
its gradient in `g`. A transition that did not move owes nothing.
"""
_surrogate_correct!(logdensity!, x, g, logp, rng) = logp
function _surrogate_correct!(d::_SurrogateDensity, x::Vector{Float64},
    g::Vector{Float64}, logp::Float64, rng::AbstractRNG)
    d.coreset === nothing && return logp
    d.transitions += 1
    x == d.xold && return logp
    d.proposed += 1
    d.batch_proposed += 1
    lp = _marginal_value(d.objective, d.gradient_method, x)
    logratio = (lp - logp) - (d.lp - d.lc)
    if isfinite(lp) && log(rand(rng)) < logratio
        d.accepted += 1
        d.batch_accepted += 1
        d.lp = lp
        d.lc = logp
        return logp
    end
    copyto!(x, d.xold)
    copyto!(g, d.gold)
    return d.lc
end

"""
    _surrogate_batch!(logdensity!, x, g, logp)

At the end of a batch of draws: a chain whose moved endpoints were accepted
too rarely for its coreset to cost less per accepted transition than a
transition in full (`_surrogate_price` at the acceptance measured) gives it
up, and samples in full from here. Returns the log density at `x`, its
gradient in `g`, under whichever density is now in use.
"""
_surrogate_batch!(logdensity!, x, g, logp) = logp
function _surrogate_batch!(d::_SurrogateDensity, x::Vector{Float64},
    g::Vector{Float64}, logp::Float64)
    d.coreset === nothing && return logp
    work = _surrogate_work(d)
    poor = d.batch_proposed >= _SURROGATE_MIN_PROPOSED &&
        d.batch_accepted * work < (work * d.share + 1) * d.batch_proposed
    d.batch_proposed = 0
    d.batch_accepted = 0
    poor || return logp
    d.coreset = nothing
    d.abandoned = true
    return _marginal_density!(d.objective, d.gradient_method, g, x)
end

"""What each chain's density did, for the sampler's record."""
function _surrogate_record(densities)
    densities === nothing && return (coreset=false, coreset_size=Int[],
        coreset_share=Float64[], coreset_error=Float64[], coreset_proposed=Int[],
        coreset_accepted=Int[], coreset_abandoned=Bool[])
    return (coreset=any(d -> d.size > 0, densities),
        coreset_size=[d.size for d in densities],
        coreset_share=[d.share for d in densities],
        coreset_error=[d.error for d in densities],
        coreset_proposed=[d.proposed for d in densities],
        coreset_accepted=[d.accepted for d in densities],
        coreset_abandoned=[d.abandoned for d in densities])
end
