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
# The coreset is fitted once per chain, at the end of its warm-up, from the
# full density's unit terms at the chain's last warm-up states
# (`_SURROGATE_TRAIN`), which by then are the posterior's; warm-up itself runs
# on the full density, so the step size and metric are tuned to it. A batch
# of draws whose endpoint acceptance falls below `_SURROGATE_MIN_ACCEPT`
# gives the coreset up for the rest of the run: one change, at a batch's end,
# after which the kernel is fixed again. A chain resumed in a new call has no
# warm-up to fit from and samples in full.
#
# What one draw costs: the trajectory's leapfrog steps at the coreset's share
# of a gradient (`_SMC_CORESET_BUDGET` at most), and one full value, no
# gradient, at a moved endpoint.

"""Warm-up states a chain's coreset is fitted from: its last this many."""
const _SURROGATE_TRAIN = 64

"""The held-out error, in log units over the training states, a chain's
coreset may have. The endpoint test reads the difference of that error
between two states, so about 0.3 keeps acceptance near NUTS's own."""
const _SURROGATE_TOL = 0.3

"""Below this share of moved endpoints accepted over a batch, the chain gives
its coreset up; read only over at least `_SURROGATE_MIN_PROPOSED`."""
const _SURROGATE_MIN_ACCEPT = 0.5
const _SURROGATE_MIN_PROPOSED = 20

"""
    _SurrogateDensity(objective, gradient_method, ndim, seed)

One chain's log density on the marginal route: the full objective until its
warm-up ends, its coreset after (`_surrogate_warmed!`), with the state the
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
    size::Int
    share::Float64
    error::Float64
    fitted::Bool
    abandoned::Bool
    seed::UInt64
end

_SurrogateDensity(objective, gradient_method, ndim::Integer, seed::Integer) =
    _SurrogateDensity(objective, Symbol(gradient_method), nothing, NaN, NaN,
        zeros(Int(ndim)), zeros(Int(ndim)), 0, 0, 0, 0, 0, 0.0, NaN, false, false,
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
    return _coreset_density!(d.objective, cs, g, x)
end

# The hooks `_run_chain` and `_sample_draws!` call. Any other density answers
# with nothing to do, so a plain closure samples exactly as before.

"""Whether a chain's density wants its warm-up states kept, to fit from."""
_surrogate_training(logdensity!) = false
_surrogate_training(d::_SurrogateDensity) = !d.fitted && _smc_coreset_capable(d.objective)

"""Keep `x` among the last `_SURROGATE_TRAIN` warm-up states."""
function _surrogate_keep!(states::Vector{Vector{Float64}}, x::Vector{Float64})
    length(states) >= _SURROGATE_TRAIN && popfirst!(states)
    push!(states, copy(x))
    return states
end

"""
    _surrogate_warmed!(logdensity!, states, x, g, logp)

At the end of a chain's warm-up: fit its coreset from the full density's
unit terms at `states`, and switch the density to it if one passes. Returns
the log density the sampling phase starts from at `x`, its gradient in `g`:
the coreset's where one was chosen, `logp` (the full) otherwise.
"""
_surrogate_warmed!(logdensity!, states, x, g, logp) = logp
function _surrogate_warmed!(d::_SurrogateDensity, states::Vector{Vector{Float64}},
    x::Vector{Float64}, g::Vector{Float64}, logp::Float64)
    d.fitted && return logp
    d.fitted = true
    (_smc_coreset_capable(d.objective) && isfinite(logp) &&
        length(states) >= 2 * _SMC_CORESET_FOLDS) || return logp
    full = try
        _smc_unit_terms(d.objective, reduce(hcat, states))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    full === nothing && return logp
    fit = _smc_coreset_fit(full.terms, full.total, _smc_unit_cost(d.objective),
        _SURROGATE_TOL, _SMC_CORESET_BUDGET, Random.Xoshiro(d.seed))
    fit === nothing && return logp
    d.coreset = CTSEMCoreset(fit.units, fit.weights)
    d.size = length(fit.units)
    d.share = fit.share
    d.error = fit.error
    gfull = copy(g)
    lc = d(g, x)
    if !isfinite(lc)
        d.coreset = nothing
        d.size = 0
        copyto!(g, gfull)
        return logp
    end
    d.lp = logp
    d.lc = lc
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
too rarely gives its coreset up, and samples in full from here. Returns the
log density at `x`, its gradient in `g`, under whichever density is now in
use.
"""
_surrogate_batch!(logdensity!, x, g, logp) = logp
function _surrogate_batch!(d::_SurrogateDensity, x::Vector{Float64},
    g::Vector{Float64}, logp::Float64)
    d.coreset === nothing && return logp
    poor = d.batch_proposed >= _SURROGATE_MIN_PROPOSED &&
        d.batch_accepted < _SURROGATE_MIN_ACCEPT * d.batch_proposed
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
