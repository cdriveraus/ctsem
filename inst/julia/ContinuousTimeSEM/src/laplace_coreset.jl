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

export ctsem_coreset_terms, ctsem_coreset_batch, ctsem_coreset_batch_gradient

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
    nunits = length(laplace.units.members)
    value = fill(-Inf, ncol)
    total = fill(NaN, ncol)
    terms = fill(NaN, nunits, ncol)
    for j in 1:ncol
        r = _laplace_coreset_try(laplace, Vector{Float64}(view(values, :, j)), false, nothing)
        r === nothing && continue
        value[j] = r.value
        total[j] = _laplace_weighted_total(laplace, r.unit_loglik)
        terms[:, j] .= r.unit_loglik
    end
    return (value=value, total=total, unit_loglik=terms)
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
    for j in 1:ncol
        r = _laplace_coreset_try(laplace, Vector{Float64}(view(values, :, j)), false, coreset)
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
    for j in 1:ncol
        r = _laplace_coreset_try(laplace, Vector{Float64}(view(values, :, j)), true, coreset)
        if r === nothing
            out[1, j] = -Inf
        else
            out[1, j] = r.value
            out[2:end, j] .= r.gradient
        end
    end
    return out
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
