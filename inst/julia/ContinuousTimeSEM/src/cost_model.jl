# What an evaluation costs: one model every pricer reads.
#
# The importance sampler (coreset_sampler.jl), tempered SMC (sample_smc.jl),
# NUTS's endpoint test (sample_surrogate.jl) and the optimiser's coreset
# stages and Newton finish (optimiser.jl) each weigh work against the
# objective's full value, the unit. They price it here, never with constants
# of their own:
#
# - `_ctsem_gradient_values(objective)`: a gradient, in values;
# - `_ctsem_unit_shares(objective)`: each unit's share of a value, and
#   `_ctsem_coreset_share`, a coreset's;
# - `_ctsem_hessian_cost(objective, npar)`: a Hessian, in gradients.
#
# Fixed by route, never timed: a price that moved with the machine's load
# would make a seeded sampler's choices, and so its draws, differ between runs.
# The figures were measured at each bench fit's estimate, one Julia thread, in
# the engine (`dev/optimbench`'s `costs` grid, juliaFit f979c331, dev3,
# 2026-10-31). Every pricer once took a gradient as two values, which was
# right on neither route.

"""A gradient's cost in values on the marginal route and on a Laplace objective
with no coordinate, whose gradient is that route's adjoint: the taped forward
pass and its reverse. Measured 1.8-4.3 over the bench's 5 to 66 parameters,
2.7-4.3 where every indicator is Gaussian and 1.8-2.0 where one is not, whose
forward pass integrates the categories; 2.05 on a 1,500-subject model."""
const _CTSEM_GRADIENT_VALUES_FILTER = 3.0

"""A gradient's cost in values on a Laplace objective with coordinates,
evaluated as a sampler evaluates it, outside an optimisation: the value is
mostly the units' inner solves, which the gradient pays too, and the gradient
adds the seeded sweeps of the trace terms and the reverse pass. Measured
1.21-1.49 on the bench's nineteen Laplace models with coordinates. Within an
optimisation, the solves anchored at the last accepted point's modes, a value
is 1.3-2.5 times cheaper and a gradient 1.4-2.0 values, but no pricer reads
the price there: the optimiser's coreset stages weigh gradients against
gradients."""
const _CTSEM_GRADIENT_VALUES_LAPLACE = 1.4

"""
    _ctsem_gradient_values(objective)

What a gradient of `objective` costs in values of it: the price every sampler's
move, leapfrog or gradient batch is weighed by against a full value.
"""
_ctsem_gradient_values(::Any) = _CTSEM_GRADIENT_VALUES_FILTER
_ctsem_gradient_values(laplace::CTSEMLaplaceObjective) =
    _laplace_filter_only(laplace) ? _CTSEM_GRADIENT_VALUES_FILTER :
    _CTSEM_GRADIENT_VALUES_LAPLACE

"""`_ctsem_gradient_values` for R, which counts what its own samplers spend in
the same unit (`imis_is()`'s record and budget)."""
ctsem_gradient_values(objective) = Float64(_ctsem_gradient_values(objective))
export ctsem_gradient_values

"""
    _ctsem_unit_shares(laplace)

Each unit's share of a full evaluation of `laplace`, a sum over units
(`_ctsem_unit_sum`), summing to one: its rows times its members
(`_laplace_unit_weights`), which is what a coreset of those units costs
against all of them. On AnomAuth, whose units' rows vary by a fifth, it
explained 37-80 per cent of the variance in a unit's time; where units are
alike in it, their times still varied by a few per cent to twofold with the
work each inner solve took, which no count before the solve can see.
"""
function _ctsem_unit_shares(laplace::CTSEMLaplaceObjective)
    w = _laplace_unit_weights(laplace)
    s = sum(w)
    return s > 0 ? w ./ s : fill(1 / max(length(w), 1), length(w))
end

"""A coreset's share of a full evaluation: its units' shares summed."""
_ctsem_coreset_share(shares::AbstractVector{<:Real}, units) =
    sum(shares[U] for U in units; init=0.0)

"""
    _ctsem_hessian_cost(objective, npar)

What one Hessian costs, in gradients: `2 npar` on the Laplace route (central
differences of the gradient), and about `npar` on the marginal route, whose
forward-over-reverse sweeps carry one lane per parameter at about a gradient
each (measured: 750 gradients' time for 715 parameters), as on a Laplace
objective with no coordinate, whose Hessian is that route's
(`_laplace_filter_hessian`). The unit the finish weighs a Hessian against its
steps in, each of which costs a gradient. The refined fit's
(`CTSEMLaplacePlaced`) is in laplace_placed.jl, beside its type.

On the bench the Laplace figure held to within 5% of the finish's own
(anchored) gradients. The marginal one is high below about a hundred
parameters: 0.3-0.5 `npar` gradients over 5 to 66, where a sweep's twelve
lanes share most of its cost, so a small fit's finish weighs its Hessians at
up to three times what they cost. Left there: those Hessians are a few
gradients either way.
"""
_ctsem_hessian_cost(laplace::CTSEMLaplaceObjective, npar::Integer) =
    _laplace_filter_only(laplace) ? Int(npar) : 2 * Int(npar)
_ctsem_hessian_cost(::Any, npar::Integer) = Int(npar)
