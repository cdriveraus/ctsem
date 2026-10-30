# Coreset evaluation for importance sampling ----------------------------------
#
# A Laplace log density is a sum of per-unit terms, each an inner mode solve.
# Over the region an importance sampler visits those terms vary in about as many
# dimensions as the sampler moves, however many units there are, so a few of
# them with fitted weights (`CTSEMCoreset`) follow the sum. These entry points
# give R's `.ctCoresetLpg()` what it needs: every unit's term at a set of draws,
# to fit and check the weights, and the coreset's value and gradient at many
# draws, solving only its units, on the pool a full evaluation uses.
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
function ctsem_coreset_terms(laplace::CTSEMLaplaceObjective, values::AbstractMatrix)
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
function ctsem_coreset_terms_gradient(laplace::CTSEMLaplaceObjective,
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
function ctsem_coreset_batch(laplace::CTSEMLaplaceObjective, values::AbstractMatrix,
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
function ctsem_coreset_batch_gradient(laplace::CTSEMLaplaceObjective,
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
    draw.inner_iterations = copy(laplace.inner_iterations)
    draw.inner_gradient = copy(laplace.inner_gradient)
    draw.inner_converged = copy(laplace.inner_converged)
    draw.hessian_repaired = copy(laplace.hessian_repaired)
    draw.mode_repaired = copy(laplace.mode_repaired)
    draw.logdet_floored = copy(laplace.logdet_floored)
    draw.last_values = Float64[]
    draw.gated_units = 0
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
solve that did not converge, or a non-finite result."""
function _laplace_coreset_try(laplace::CTSEMLaplaceObjective, x::Vector{Float64},
    gradient::Bool, coreset::Union{Nothing,CTSEMCoreset})
    r = try
        ctsem_laplace_evaluate(laplace, x; gradient=gradient, coreset=coreset)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    r === nothing && return nothing
    ok = r.converged && isfinite(r.value) &&
        (!gradient || (r.gradient !== nothing && all(isfinite, r.gradient)))
    return ok ? r : nothing
end
