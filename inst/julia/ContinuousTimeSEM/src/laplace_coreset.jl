# Coreset evaluation -----------------------------------------------------------
#
# A log density that is a sum of per-unit terms -- the Laplace route's units,
# each an inner mode solve, or the marginal route's subjects, each a filter
# pass -- varies over the region a sampler or an optimiser's stage visits in
# about as many dimensions as it moves, however many units there are, so a few
# of them with fitted weights (`CTSEMCoreset`) follow the sum. These entry
# points give R's `.ctCoresetLpg()` and the engine's samplers what they need:
# every unit's term at a set of draws, to fit and check the weights, and the
# coreset's value and gradient at many draws, evaluating only its units, on the
# pool a full evaluation uses. Each takes either route's objective, through its
# view as a sum over units (`_ctsem_unit_sum`).
#
# Each refuses a draw the way `ctsem_evaluate_batch` does (`_ctsem_probe_value`):
# a point the model cannot evaluate, or whose inner solve did not converge,
# comes back `-Inf` and is never raised, since a sampler's tails are full of them.

export ctsem_coreset_terms, ctsem_coreset_terms_gradient, ctsem_coreset_batch,
    ctsem_coreset_batch_gradient, ctsem_set_coreset_draw_parallel!

"""
Whether a coreset's draws are evaluated several at once (`_laplace_coreset_each`).
On; off evaluates them one after another, each over its units in parallel,
which gives the same numbers and is there to time against and to check that it
does. `ctsem_set_coreset_draw_parallel!` sets it.
"""
const _LAPLACE_CORESET_DRAWS = Ref(true)

"""Set whether a coreset's draws are evaluated several at once."""
ctsem_set_coreset_draw_parallel!(on::Bool) = (_LAPLACE_CORESET_DRAWS[] = on)

_coreset_units(laplace::CTSEMLaplaceObjective) = length(laplace.units.members)

"""
    _ctsem_unit_sum(objective)

`objective` as a sum over units, which is what every coreset stands in for --
the optimiser's coreset stages, the importance sampler's, DAIS's, tempered
SMC's and NUTS's endpoint test alike -- or `nothing` for an objective that is
no such sum. A Laplace objective is its own. The marginal route's objective is
a sum over its subjects, and a Laplace objective with no level over it is that
sum: a unit a subject, nothing to solve (`_laplace_filter_only`), the subjects
summed through the marginal route's own adjoint chunk with their sampling
weights or a coreset's. So a coreset is one thing on both routes, and the
marginal route has no coreset code of its own. The view is built when first
asked for and kept on the objective (`unit_sum`); it reads the objective's
prior, weights and subjects as they are at each evaluation. None for a
CTSEMObjective whose filter integrates deviations (`params.hybrid`): that is
the inner objective of a Laplace one, whose own evaluation profiles them,
where the view would take their integral, another function.
"""
_ctsem_unit_sum(objective) = nothing
_ctsem_unit_sum(laplace::CTSEMLaplaceObjective) = laplace
function _ctsem_unit_sum(objective::CTSEMObjective)
    _hybrid_on(objective.params.hybrid) && return nothing
    built = objective.unit_sum
    built === nothing || return built::CTSEMLaplaceObjective
    # Under a lock: a sampler's chains or particles can ask at once.
    return lock(_CTSEM_UNIT_SUM_LOCK) do
        objective.unit_sum === nothing && (objective.unit_sum =
            CTSEMLaplaceObjective(objective, CTSEMLaplaceSpec(CTSEMLaplaceLevel[])))
        objective.unit_sum::CTSEMLaplaceObjective
    end
end
const _CTSEM_UNIT_SUM_LOCK = ReentrantLock()

"""`_ctsem_unit_sum(objective)`, refusing an objective that is no sum over units."""
function _ctsem_units_of(objective)
    L = _ctsem_unit_sum(objective)
    L === nothing && throw(ArgumentError(string("a coreset needs an objective ",
        "that is a sum over units, not a ", nameof(typeof(objective)))))
    return L
end

# The entry points on any objective that is a sum over units.
ctsem_coreset_terms(objective, values::AbstractMatrix) =
    ctsem_coreset_terms(_ctsem_units_of(objective), values)
ctsem_coreset_terms_gradient(objective, values::AbstractMatrix,
    scored::Union{Integer,AbstractVector}) =
    ctsem_coreset_terms_gradient(_ctsem_units_of(objective), values, scored)
ctsem_coreset_batch(objective, values::AbstractMatrix, units::AbstractVector,
    weights::AbstractVector) =
    ctsem_coreset_batch(_ctsem_units_of(objective), values, units, weights)
ctsem_coreset_batch_gradient(objective, values::AbstractMatrix, units::AbstractVector,
    weights::AbstractVector) =
    ctsem_coreset_batch_gradient(_ctsem_units_of(objective), values, units, weights)

"""
    ctsem_coreset_terms(laplace, values)

Every unit's term at each column of `values` (`npar x ndraws`), for fitting a
coreset's weights, with

- `value`: the log posterior, `-Inf` where the draw cannot be used,
- `total`: the units' terms summed with their sampling weights (the log
  posterior less the log prior), what a coreset's weighted sum stands in for,
- `unit_loglik`: `nunits x ndraws`, each unit's term, NaN where unusable.

`ctsem_unit_terms` is the leave-one-unit-out version, which raises on what a
proposal's tails produce; this is for draws from a sampler.
"""
function ctsem_coreset_terms(laplace::_CTSEMUnits, values::AbstractMatrix)
    ncol = size(values, 2)
    ncol >= 1 || throw(ArgumentError("values must have at least one column"))
    nunits = _coreset_units(laplace)
    value = fill(-Inf, ncol)
    total = fill(NaN, ncol)
    terms = fill(NaN, nunits, ncol)
    for j in 1:ncol
        r = _coreset_column(laplace, Vector{Float64}(view(values, :, j)), false, nothing)
        r === nothing && continue
        value[j] = r.value
        total[j] = r.total
        terms[:, j] .= r.terms
    end
    return (value=value, total=total, unit_loglik=terms)
end

"""
    _coreset_column(objective, x, gradient, scores)

The full density at `x` as a coreset reads it: `value`, `total` (the value less
the log prior, what the units' weighted terms stand in for), each unit's term
(`terms`) and, with `gradient`, the gradient at `x`'s length, with each unit's
score in `scores`' rows where that is given. `nothing` where the draw cannot
be used.
"""
function _coreset_column(laplace::CTSEMLaplaceObjective, x::Vector{Float64},
    gradient::Bool, scores::Union{Nothing,Matrix{Float64}})
    r = try
        ctsem_laplace_evaluate(laplace, x; gradient=gradient, scores=scores)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    (r !== nothing && r.converged && isfinite(r.value) && (!gradient ||
        (r.gradient !== nothing && all(isfinite, r.gradient)))) || return nothing
    return (value=r.value, total=_laplace_weighted_total(laplace, r.unit_loglik),
        terms=r.unit_loglik, gradient=r.gradient)
end

"""
    ctsem_coreset_terms_gradient(laplace, values, scored)

A full gradient batch that also teaches a coreset: per column of `values`
(`npar x ndraws`) `value`, `total` and `unit_loglik` as `ctsem_coreset_terms`
gives them, `gradient` (`npar x ndraws`) as `ctsem_evaluate_batch_gradient`
does, and for the columns `scored` names each unit's score, `scores`
(`nunits x npar x length(scored)`): its term's gradient, with the units'
sampling weights and without the prior, so that the units' scores sum to the
gradient of `total`. Given a count instead, the scores are the `scored`
usable columns' of highest `total`, where a sampler's weight sits; every
column's are formed as its gradient is, so choosing costs nothing, and only
that many are held. `scored` in the result says which columns they are. One
gradient evaluation a column, what the batch gradient costs; `-Inf`, NaN
where a draw is unusable, by `_ctsem_probe_value_gradient`'s rule.
"""
function ctsem_coreset_terms_gradient(laplace::_CTSEMUnits,
    values::AbstractMatrix, scored::Union{Integer,AbstractVector})
    n, ncol = size(values)
    ncol >= 1 || throw(ArgumentError("values must have at least one column"))
    nunits = _coreset_units(laplace)
    value = fill(-Inf, ncol)
    total = fill(NaN, ncol)
    terms = fill(NaN, nunits, ncol)
    gradient = fill(NaN, n, ncol)
    listed = scored isa AbstractVector
    k = listed ? length(scored) : clamp(Int(scored), 0, ncol)
    slots = [fill(NaN, nunits, n) for _ in 1:k]
    slotcol = listed ? collect(Int, scored) : zeros(Int, k)
    slottotal = fill(-Inf, k)
    S = zeros(Float64, nunits, n)
    for j in 1:ncol
        want = listed ? findfirst(==(j), slotcol) : (k > 0 ? argmin(slottotal) : nothing)
        want === nothing || fill!(S, 0.0)
        r = _coreset_column(laplace, Vector{Float64}(view(values, :, j)), true,
            want === nothing ? nothing : S)
        r === nothing && continue
        value[j] = r.value
        total[j] = r.total
        terms[:, j] .= r.terms
        gradient[:, j] .= r.gradient
        want === nothing && continue
        if listed
            copyto!(slots[want], S)
        elseif total[j] > slottotal[want]
            copyto!(slots[want], S)
            slotcol[want] = j
            slottotal[want] = total[j]
        end
    end
    keep = listed ? collect(1:k) : sort(findall(>(0), slotcol); by=i -> slotcol[i])
    scores = Array{Float64}(undef, nunits, n, length(keep))
    for (a, i) in enumerate(keep)
        scores[:, :, a] .= slots[i]
    end
    return (value=value, total=total, unit_loglik=terms, gradient=gradient,
        scores=scores, scored=slotcol[keep])
end

"""
    ctsem_coreset_batch(laplace, values, units, weights)

The log prior plus `sum(weights .* terms)` over `units` at each column of
`values`, solving only those units (`CTSEMCoreset`); `-Inf` where unusable.
"""
function ctsem_coreset_batch(laplace::_CTSEMUnits, values::AbstractMatrix,
    units::AbstractVector, weights::AbstractVector)
    ncol = size(values, 2)
    ncol >= 1 || throw(ArgumentError("values must have at least one column"))
    coreset = _laplace_coreset(laplace, units, weights)
    out = fill(-Inf, ncol)
    results = _laplace_coreset_each(laplace, values, coreset, false)
    for j in 1:ncol
        r = results[j]
        r === nothing || (out[j] = r.value)
    end
    return out
end

"""
    ctsem_coreset_batch_gradient(laplace, values, units, weights)

`ctsem_coreset_batch` with the gradient: an `(n + 1) x ncol` matrix whose first
row holds the values and the rest the gradients, as
`ctsem_evaluate_batch_gradient` returns them; value `-Inf` and gradient NaN
where unusable.
"""
function ctsem_coreset_batch_gradient(laplace::_CTSEMUnits,
    values::AbstractMatrix, units::AbstractVector, weights::AbstractVector)
    n, ncol = size(values)
    ncol >= 1 || throw(ArgumentError("values must have at least one column"))
    coreset = _laplace_coreset(laplace, units, weights)
    out = fill(NaN, n + 1, ncol)
    results = _laplace_coreset_each(laplace, values, coreset, true)
    for j in 1:ncol
        r = results[j]
        if r === nothing
            out[1, j] = -Inf
        else
            out[1, j] = r.value
            out[2:end, j] .= r.gradient
        end
    end
    return out
end

"""
    _laplace_coreset_each(laplace, values, coreset, gradient)

The coreset evaluated at every column of `values`, several draws at once: the
pool works through the draws, each on its own copy of the objective's per-unit
state (`_laplace_draw_copy`), and a draw's units take whatever slots the draws
leave free. Over a coreset's units alone most of the pool sat idle whenever it
had more workers than the coreset has units. One result per column, `nothing`
where the draw cannot be used (`_laplace_coreset_try`).
"""
function _laplace_coreset_each(laplace::CTSEMLaplaceObjective, values::AbstractMatrix,
    coreset::CTSEMCoreset, gradient::Bool)
    ncol = size(values, 2)
    results = Vector{Any}(nothing, ncol)
    if !_LAPLACE_CORESET_DRAWS[] || ncol == 1
        for j in 1:ncol
            results[j] = _laplace_coreset_try(laplace, Vector{Float64}(view(values, :, j)),
                gradient, coreset)
        end
        return results
    end
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, 1:ncol) do j
        results[j] = _laplace_coreset_try(_laplace_draw_copy(laplace, coreset),
            Vector{Float64}(view(values, :, j)), gradient, coreset)
        return true
    end
    return results
end

"""
    _laplace_draw_copy(laplace, coreset)

`laplace` for one draw evaluated while others are: fresh modes for the
coreset's units and its own per-unit flags, the rest shared. An evaluation
writes each unit's mode and flags on the objective, so draws evaluated at once
on one objective would overwrite each other's; the subjects and the
specification are only read, and the workspaces are per slot, which the pool
gives one task at a time. Modes are rebound per unit
(`_laplace_solve_unit_mode!`), never written in place, so a shallow copy of the
others is enough: every field is the original's unless replaced here, so a
field added to the objective is shared by a draw's copy without an edit.
"""
function _laplace_draw_copy(laplace::CTSEMLaplaceObjective, coreset::CTSEMCoreset)
    draw = _ctsem_shallow_copy(laplace)
    modes = copy(laplace.modes)
    for U in coreset.units
        modes[U] = zeros(Float64, length(laplace.modes[U]))
    end
    draw.modes = modes
    alt = copy(laplace.alt_modes)
    for U in coreset.units
        alt[U] = Vector{Float64}[]
    end
    draw.alt_modes = alt
    blk = copy(laplace.block_modes)
    for U in coreset.units
        blk[U] = CTSEMBlockMode[]
    end
    draw.block_modes = blk
    draw.inner_iterations = copy(laplace.inner_iterations)
    draw.inner_gradient = copy(laplace.inner_gradient)
    draw.inner_converged = copy(laplace.inner_converged)
    draw.hessian_repaired = copy(laplace.hessian_repaired)
    draw.mode_repaired = copy(laplace.mode_repaired)
    draw.logdet_floored = copy(laplace.logdet_floored)
    draw.last_values = Float64[]
    draw.soft_quadrature_units = 0
    return draw
end

"""A new object of `x`'s type holding `x`'s fields, by its default
constructor: for a copy that replaces a few of them, with the list of fields
in the struct alone."""
@generated function _ctsem_shallow_copy(x::T) where {T}
    return :(T($((:(getfield(x, $i)) for i in 1:fieldcount(T))...)))
end

"""One evaluation for the entry points above, or `nothing` where the draw
cannot be used: an error that is not a bug (`_ctsem_must_propagate`), an inner
solve that did not converge, or a non-finite result. With `scores`, each
unit's score is written to its row, and they must be finite too."""
function _laplace_coreset_try(laplace::CTSEMLaplaceObjective, x::Vector{Float64},
    gradient::Bool, coreset::Union{Nothing,CTSEMCoreset};
    scores::Union{Nothing,Matrix{Float64}}=nothing)
    r = try
        ctsem_laplace_evaluate(laplace, x; gradient=gradient, coreset=coreset,
            scores=scores)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    r === nothing && return nothing
    ok = r.converged && isfinite(r.value) &&
        (!gradient || (r.gradient !== nothing && all(isfinite, r.gradient))) &&
        (scores === nothing || all(isfinite, scores))
    return ok ? r : nothing
end

# The overshoot probe's screen -------------------------------------------------

export ctsem_set_probe_screen!

"""
Whether the overshoot probe judges its trials on samples of units before
evaluating them in full (`CTSEMProbeScreen`). On; off evaluates every trial in
full, which is the oracle. `ctsem_set_probe_screen!` sets it.
"""
const _CTSEM_PROBE_SCREEN = Ref(true)

"""Set whether the overshoot probe judges its trials on samples of units."""
ctsem_set_probe_screen!(on::Bool) = (_CTSEM_PROBE_SCREEN[] = on; nothing)

"""The screen's first sample is this share of the units a trial moves, at least
`_CTSEM_SCREEN_FIRST` of them; each later one is twice the last."""
const _CTSEM_SCREEN_SHARE = 16

"""The fewest units the screen's first sample takes."""
const _CTSEM_SCREEN_FIRST = 64

"""The fewest units a stratum is sampled by; a stratum no larger is taken whole."""
const _CTSEM_SCREEN_STRATUM = 16

"""The standard errors the screen's bound adds to its estimate."""
const _CTSEM_SCREEN_Z = 4.0

"""
    CTSEMProbeScreen

What the overshoot probe (`_ctsem_overshot`) needs to judge a trial point from
some of an objective's units: the objective as a sum over units
(`_ctsem_unit_sum`), the estimate, and, formed at the first trial it judges,
each unit's term at the estimate and which units read each parameter.

The probe asks of a trial only whether it beats the estimate by more than its
tolerance, and nearly every trial loses by hundreds of nats or more. Each was a
full evaluation: 42 s on a 7,141-unit ordinal Laplace fit at eight threads, and
up to two hundred of them a search. A trial's change from the estimate is a sum
over units, and a random sample of the units estimates it without bias wherever
the trial lies. A coreset fitted to draws near the estimate does not: its
weights say nothing about a point tenths of a coordinate away. So each trial is
first evaluated on a sample (`_ctsem_probe_trial`) of a sixteenth of the units
it moves (at least 64), then twice that, and so on, and is set aside once its
gain's upper bound is under the tolerance. A trial the samples have not
settled by half the units is evaluated in full, and so is every trial whose
bound allows a gain, so an overshoot is always found on the full objective.
The share rather than a fixed count, because a gain carried by one unit in a
hundred is missed by a sample of 64 half the time, and by a sixteenth of 7,141
units about once in ninety.

The units are stratified by which of the trial's moved coordinates they read.
A unit that reads none of them is unchanged and not sampled. A unit that reads
a coordinate few others read (a group's predictor effect, an item few
answered, a spline's local basis) falls in a stratum of its own, which is taken
whole where it is small, so a change in those units is never left to the
chance of drawing them. Readers are the units with a nonzero score at the
estimate or beside it, as the localised Hessian reads them
(`_laplace_hessian_readers`). A parameter no unit reads there is taken as read
by every unit, since a transform flat past the arithmetic scores zero too.

Each unit's first-order change, its score at the estimate times the trial's
move, is known for every unit, so within a stratum only the rest is sampled: a
regression estimator, the sampled changes regressed on their first-order
predictions. Near the estimate the predictions are large and of either sign
while their sum, the gradient's, is near zero, and the total is the small
second-order rest: sampled plainly, the first trials of a 3,000-unit fit's
search lost 29 to 290 nats under a spread no half of the units could settle.
Far out the predictions say little, the regression gives them little weight,
and the estimate is the plain one.

The bound is the stratified estimate, four standard errors, and a tail term:
each stratum's unsampled units may include, at the rule of three's rate `3/m`,
units whose residual exceeds the stratum's mean by as much as the largest
sampled one does. What the screen can miss is a unit further out than every
sampled one, at a rate under `3/m`. Each trial is sampled afresh, seeded by
the point, and the gain the probe looks for (a scale back from a collapse, a
coordinate off a flat transform) moves every unit that reads it.

A sampled unit whose inner solve fails refuses the trial, as the full
evaluation would. Where every stratum is taken whole the change is exact, and
only a gain is evaluated again in full.

Measured on a simulated 3,000-unit ordinal item-level Laplace fit (seven
effects a unit, 54 parameters, one item answered by one unit in twenty; dev4
at eight threads; `scripts/probescreen-2026-11`, review/PROBE-screen-2026-11-20.md).
The search over every coordinate at its estimate, 162 trials and no gain, took
486,000 unit evaluations and 162 s unscreened, 105,140 with the plain
estimator and 56,184 (18 s) with the regression, the same verdict. Every trial
of the search judged both ways at the estimate and at three collapsed scales
(item 2's, and item 1's where its effect varied in 3 and in 0.3 per cent of
the units): no trial that gained was set aside, and the nearest that was
gained 8e-5. Judged on samples drawn without the nine units whose effect
varied, a trial gaining 9.2 nats still estimated a gain at every sample size:
their first-order predictions, known whether drawn or not, are in the
regression's total. Its cost is the gradients at the estimate, about three values: where the
first trials gain, the probe costs that much more than it did.
"""
mutable struct CTSEMProbeScreen{L}
    units::L
    minimizer::Vector{Float64}
    ready::Union{Nothing,Bool}
    terms::Vector{Float64}
    weights::Vector{Float64}
    value::Float64
    prior::Float64
    readers::BitMatrix
    scores::Matrix{Float64}
end

"""
    _ctsem_probe_screen(objective, minimizer)

A screen for the overshoot probe at `minimizer` (`CTSEMProbeScreen`), or
`nothing` where none can serve: screening off, an objective that is no sum over
units, or fewer units than two first samples, where a sample is most of a full
evaluation. It evaluates nothing until its first trial.
"""
function _ctsem_probe_screen(objective, minimizer::AbstractVector)
    _CTSEM_PROBE_SCREEN[] || return nothing
    L = _ctsem_probe_units(objective)
    L === nothing && return nothing
    _coreset_units(L) >= 2 * _CTSEM_SCREEN_FIRST || return nothing
    return CTSEMProbeScreen(L, collect(Float64, minimizer), nothing, Float64[],
        Float64[], NaN, NaN, falses(0, 0), zeros(0, 0))
end

"""
    _ctsem_probe_units(objective)

The objective as the sum over units the probe's screen samples: its
`_ctsem_unit_sum`, and a variational objective itself, whose units' bounds it
evaluates a few at a time without touching the trial (variational.jl). The
screen asks of it `_coreset_units`, `_laplace_unit_weight`,
`_laplace_coreset`, `_laplace_coreset_try` and `_ctsem_screen_prior`.
"""
_ctsem_probe_units(objective) = _ctsem_unit_sum(objective)

"""The prior term at `x` of a sum over units the screen samples."""
_ctsem_screen_prior(L::CTSEMLaplaceObjective, x::AbstractVector) =
    _ctsem_log_prior(L.objective, _laplace_extend(L, x))

"""Forms `s`'s units' terms, scores and readers at its estimate, once, from one
gradient with the units' scores (two where some parameter has readers and
non-readers); whether it can serve."""
function _ctsem_screen_ready!(s::CTSEMProbeScreen)
    s.ready === nothing || return s.ready::Bool
    s.ready = false
    L = s.units
    N = _coreset_units(L)
    n = length(s.minimizer)
    weights = Float64[_laplace_unit_weight(L, U) for U in 1:N]
    scores = zeros(Float64, N, n)
    _CTSEM_OPCOUNT.overshoot_gradients[] += 1
    r = _laplace_coreset_try(L, s.minimizer, true, nothing; scores=scores)
    r === nothing && return false
    # A unit of no weight adds nothing, and the marginal route leaves it unevaluated.
    all(U -> weights[U] == 0 || isfinite(r.unit_loglik[U]), 1:N) || return false
    readers = scores .!= 0
    if any(j -> 0 < count(view(readers, :, j)) < N, 1:n)
        # At the localised Hessian's difference step.
        moved = zeros(Float64, N, n)
        _CTSEM_OPCOUNT.overshoot_gradients[] += 1
        _laplace_coreset_try(L, _laplace_reader_point(s.minimizer, 1e-4), true,
            nothing; scores=moved) === nothing && return false
        readers .|= moved .!= 0
    end
    for j in 1:n
        any(view(readers, :, j)) || (readers[:, j] .= true)
    end
    for U in 1:N
        weights[U] == 0 && (readers[U, :] .= false)
    end
    s.terms = Float64[weights[U] == 0 ? 0.0 : weights[U] * r.unit_loglik[U] for U in 1:N]
    s.weights = weights
    s.value = r.value
    s.prior = _ctsem_screen_prior(L, s.minimizer)
    s.readers = readers
    s.scores = scores
    s.ready = true
    return true
end

"""The units that read any of `moved`, grouped by which of them each reads:
each group in unit order, the groups by their first unit."""
function _ctsem_screen_strata(readers::BitMatrix, moved::Vector{Int})
    isempty(moved) && return Vector{Int}[]
    groups = Dict{BitVector,Vector{Int}}()
    for U in axes(readers, 1)
        reads = readers[U, moved]
        any(reads) && push!(get!(Vector{Int}, groups, reads), U)
    end
    return sort!(collect(values(groups)); by=first)
end

"""How many units of a stratum of `size` a sample of `target` of `affected`
takes: its share, at least `_CTSEM_SCREEN_STRATUM`, and all of a small one."""
_ctsem_screen_allot(size::Int, target::Int, affected::Int) =
    size <= _CTSEM_SCREEN_STRATUM ? size :
    min(size, max(_CTSEM_SCREEN_STRATUM, cld(target * size, affected)))

"""
    _ctsem_screen_bound(changes, predicted, totals, sizes)

The strata's change estimated from their sampled `changes`, and its upper
bound: four standard errors and the tail term over (`CTSEMProbeScreen`). A
stratum taken whole adds its changes. A sampled one adds the regression
estimate on the sampled units' first-order `predicted` changes, whose sum over
the whole stratum is its entry of `totals`: with `b` the sampled changes'
slope on the predictions, `b * total + size * mean(change - b * predicted)`,
its variance the residuals'.
"""
function _ctsem_screen_bound(changes::Vector{Vector{Float64}},
        predicted::Vector{Vector{Float64}}, totals::Vector{Float64},
        sizes::Vector{Int})
    estimate = 0.0
    variance = 0.0
    tail = 0.0
    for h in eachindex(changes)
        y = changes[h]
        m = length(y)
        n = sizes[h]
        if m == n
            estimate += sum(y)
            continue
        end
        p = predicted[h]
        pm = mean(p)
        ym = mean(y)
        spread = sum(q -> (q - pm)^2, p)
        slope = spread > 0 ? sum((y[k] - ym) * (p[k] - pm) for k in 1:m) / spread : 0.0
        r = y .- slope .* p
        rm = mean(r)
        estimate += slope * totals[h] + n * rm
        # Two degrees of freedom spent on the residuals: their mean and slope.
        variance += n^2 * (1 - m / n) * sum(e -> (e - rm)^2, r) / ((m - 2) * m)
        tail += (n - m) * min(1.0, 3 / m) * (maximum(r) - rm)
    end
    return (estimate=estimate,
        upper=estimate + _CTSEM_SCREEN_Z * sqrt(variance) + tail)
end

"""Each of `units`' terms at `x`, NaN elsewhere, or `nothing` where one of
them cannot be evaluated."""
function _ctsem_screen_terms(L, x::Vector{Float64}, units::Vector{Int})
    _CTSEM_OPCOUNT.overshoot_samples[] += 1
    _CTSEM_OPCOUNT.overshoot_units[] += length(units)
    r = _laplace_coreset_try(L, x, false,
        _laplace_coreset(L, units, ones(length(units))))
    (r === nothing || !all(U -> isfinite(r.unit_loglik[U]), units)) && return nothing
    return r.unit_loglik
end

"""
    _ctsem_probe_trial(screen, objective, x, value, tolerance)

The overshoot probe's value at trial point `x`, judged against the estimate's
`value`. With no screen, or one that cannot serve, it is `_ctsem_probe_value`.
Otherwise it is what the screen's samples settle (`CTSEMProbeScreen`): `-Inf`
for a trial they show cannot gain `tolerance` or that a unit refuses, the exact
value where they took every unit that moves, and the full value where they
settle nothing or allow a gain.
"""
_ctsem_probe_trial(::Nothing, objective, x, value, tolerance) =
    _ctsem_probe_full(objective, x)

function _ctsem_probe_trial(s::CTSEMProbeScreen, objective, x::Vector{Float64},
        value::Real, tolerance::Real)
    _ctsem_screen_ready!(s) || return _ctsem_probe_full(objective, x)
    L = s.units
    N = length(s.terms)
    moved = Int[p for p in eachindex(x) if x[p] != s.minimizer[p]]
    strata = _ctsem_screen_strata(s.readers, moved)
    sizes = Int[length(h) for h in strata]
    affected = sum(sizes; init=0)
    # What no unit carries: the prior's change, and how far the estimate's own
    # value lies from the one the probe compares against.
    shift = _ctsem_screen_prior(L, x) - s.prior + (s.value - value)
    rng = Random.Xoshiro(hash(x))
    for h in strata
        shuffle!(rng, h)
    end
    # Each unit's first-order change, and their sum over each stratum.
    step = Float64[x[p] - s.minimizer[p] for p in moved]
    first_order(U) = sum((s.scores[U, moved[k]] * step[k] for k in eachindex(moved));
        init=0.0)
    totals = Float64[sum(first_order, h; init=0.0) for h in strata]
    taken = zeros(Int, length(strata))
    changes = [Float64[] for _ in strata]
    predicted = [Float64[] for _ in strata]
    target = max(_CTSEM_SCREEN_FIRST, cld(affected, _CTSEM_SCREEN_SHARE))
    while true
        want = Int[_ctsem_screen_allot(sizes[h], target, affected) for h in eachindex(strata)]
        # Past half the units a sample costs most of what it would spare.
        sum(want; init=0) > N ÷ 2 && break
        units = Int[]
        for h in eachindex(strata), k in (taken[h] + 1):want[h]
            push!(units, strata[h][k])
        end
        if !isempty(units)
            terms = _ctsem_screen_terms(L, x, units)
            if terms === nothing
                _CTSEM_OPCOUNT.overshoot_screened[] += 1
                return -Inf
            end
            for h in eachindex(strata), k in (taken[h] + 1):want[h]
                U = strata[h][k]
                push!(changes[h], s.weights[U] * terms[U] - s.terms[U])
                push!(predicted[h], first_order(U))
            end
            copyto!(taken, want)
        end
        if taken == sizes
            gain = shift + sum(sum, changes; init=0.0)
            gain > tolerance && return _ctsem_probe_full(objective, x)
            _CTSEM_OPCOUNT.overshoot_screened[] += 1
            return Float64(value) + gain
        end
        if shift + _ctsem_screen_bound(changes, predicted, totals, sizes).upper < tolerance
            _CTSEM_OPCOUNT.overshoot_screened[] += 1
            return -Inf
        end
        target *= 2
    end
    return _ctsem_probe_full(objective, x)
end
