using LinearAlgebra
using ForwardDiff
using ComponentArrays

################################################################################
# Parameter transforms
################################################################################

"""
    CTSEMRowContext(state, pars, tdpreds, tipreds, time, dt, subject, row)

All row-dependent expressions receive this stable forward-pass context. Keeping
the expression boundary here lets a future adjoint share exactly the same
primitive operations as the ForwardDiff path.
"""
struct CTSEMRowContext{S,P,TD,TI,TM,DT}
    state::S
    pars::P
    tdpreds::TD
    tipreds::TI
    time::TM
    dt::DT
    subject::Int
    row::Int
end

"""
    apply_complex_transforms_at_indices!(all_params, indices, transforms, context)

Apply state-dependent parameter transforms at selected flattened positions.

Each transform is called as `transform(context)` and its result is written
to the corresponding entry of `all_params`.
"""
function apply_complex_transforms_at_indices!(
    all_params::AbstractVector,
    indices::AbstractVector{Int},
    transforms,
    context::CTSEMRowContext,
)
    @boundscheck length(indices) == length(transforms) || throw(DimensionMismatch("Transform indices/functions mismatch"))
    @inbounds for k in eachindex(indices)
        all_params[indices[k]] = transforms[k](context)
    end
    return all_params
end

function apply_complex_transforms_at_indices!(all_params, indices, transforms, state, pars)
    @boundscheck length(indices) == length(transforms) || throw(DimensionMismatch("Transform indices/functions mismatch"))
    @inbounds for k in eachindex(indices)
        all_params[indices[k]] = transforms[k](state, pars)
    end
    return all_params
end

"""Copy population raw parameters and add TI effects into an isolated subject buffer."""
function _materialize_subject_values!(subject_values::AbstractVector, values::AbstractVector,
    sp::EKFParameters, tipreds::AbstractVector)
    subject_values isa Vector || throw(ArgumentError("subject parameter buffer must be resizable"))
    resize!(subject_values, length(values))
    copyto!(subject_values, values)
    @boundscheck length(sp.ti_parameter_indices) == length(sp.ti_predictor_indices) == length(sp.ti_coefficient_indices) ||
        throw(DimensionMismatch("TI effect metadata has inconsistent lengths"))
    @inbounds for i in eachindex(sp.ti_parameter_indices)
        parameter = sp.ti_parameter_indices[i]
        predictor = sp.ti_predictor_indices[i]
        coefficient = sp.ti_coefficient_indices[i]
        parameter <= length(subject_values) || throw(BoundsError(subject_values, parameter))
        predictor <= length(tipreds) || throw(BoundsError(tipreds, predictor))
        coefficient <= length(values) || throw(BoundsError(values, coefficient))
        subject_values[parameter] += values[coefficient] * tipreds[predictor]
    end
    return subject_values
end

"""
    TIMissingRecipe(base, predictor_index, parameter_index)

How to rebuild one subject's TI-predictor row when some of its cells are
sampled rather than observed.

`base` is the subject's row of `tipred_data` exactly as R sent it -- a fixed
`Vector{Float64}`, real data at the observed positions and an unused filler
(`UNSET_PARAMETER`) at the missing ones, since those are overwritten on every
call and never read as data. `predictor_index` names which columns are
missing for this subject and `parameter_index` gives, for each of those
columns in the same order, the raw-parameter position that samples it.

A subject with no missing cells is never wrapped in this: its
`ContinuousEKFObjective` holds the plain `Vector{Float64}` it always has, so
nothing about that path changes. See `_ctsem_tipred_vector` for the dispatch
this exists to drive.
"""
struct TIMissingRecipe{V<:AbstractVector{<:Real}}
    base::V
    predictor_index::Vector{Int}
    parameter_index::Vector{Int}
end

"""
    _ctsem_tipred_vector(tipreds, p)

The TI-predictor vector a row evaluation should read, for the raw parameter
vector `p` of the current call.

Two methods, dispatched on the *stored* type, not a runtime flag:

  * `tipreds::AbstractVector` (today's only case, and every case for a model
    with no missing TI predictors) -- returned unchanged. This method is
    exactly the identity function and the compiler resolves it at the call
    site, so a model with nothing missing pays nothing for this feature
    existing: no branch, no allocation, no new code on its path.
  * `tipreds::TIMissingRecipe` -- a fresh `Vector{T}` (`T` = `eltype(p)`, so a
    `ForwardDiff.Dual` during a gradient) built from the fixed base row with
    each missing cell overwritten by `p[parameter_index]`. This is the
    "honest cost" the spec accepts for sampling a missing predictor: one
    small per-subject allocation per evaluation, paid only by subjects that
    actually have a missing cell.
"""
@inline _ctsem_tipred_vector(tipreds::AbstractVector, p::AbstractVector) = tipreds

@inline function _ctsem_tipred_vector(spec::TIMissingRecipe, p::AbstractVector{T}) where {T}
    buffer = Vector{T}(undef, length(spec.base))
    copyto!(buffer, spec.base)
    @inbounds for k in eachindex(spec.predictor_index)
        buffer[spec.predictor_index[k]] = p[spec.parameter_index[k]]
    end
    return buffer
end

"""
    _materialize_all_params!(all_params, values, sp)

Materialize the full transformed parameter vector for an EKF evaluation.

Mutable parameters are produced by `sp.regular_transforms(values)`, and fixed
values from `sp` are copied into their fixed positions.
"""
function _materialize_all_params!(all_params::AbstractVector, values::AbstractVector, sp::EKFParameters)
    @boundscheck begin
        nall = length(sp.mutables)
        length(all_params) == nall || throw(DimensionMismatch("all_params length must match parameter mask length"))

        nmut = count(sp.mutables)
        ntf = length(sp.regular_transforms)
        nmut == ntf || throw(DimensionMismatch("Number of mutable positions must match number of regular transforms"))
    end

    tf_idx = 1
    @inbounds for idx in eachindex(sp.mutables)
        if sp.mutables[idx]
            all_params[idx] = sp.regular_transforms[tf_idx](values)
            tf_idx += 1
        end
    end
    map_fixed_values!(all_params, sp.fixed_indices, sp.fixed_values)
    return all_params
end

"""
Materialising-transform derivative below which a raw parameter is treated as
saturated: its transform has stopped responding, so no gradient reaches it for
reasons that have nothing to do with the data.

Not a bound on the raw coordinate -- see `_ctsem_saturated_parameters` -- but
on the one quantity that actually says whether a transform has gone flat.
Every ctsem transform is either always non-saturating (identity: derivative
exactly its multiplier, at any raw magnitude) or asymptotically flat on
(at least) one side, and the two regimes are far apart wherever that side is
reached. Measured with `ForwardDiff` across the transform strings
`ctJuliaBackend.R` actually writes:

  transform                                  raw    derivative
  `2/(1+exp(-x))-1` (correlation)             10     9.1e-5
  `2/(1+exp(-x))-1` (correlation)             17     8.3e-8
  `-log1p_exp(x)` (drift, negative side)    -18.5    9.2e-9   (see
      `test_state_sampling.jl`, "a count model fits over the joint density":
      this is the drift coordinate that test's own comment already calls
      unidentified)
  `-(1e-6+2log1p_exp(-2x))` (drift diagonal)  10     8.2e-9
  `1e-10+5log1p_exp(2x)` (variance)            0     5.0       (never
      saturates for raw > 0; only for raw very negative, i.e. variance -> 0)

`1e-6` sits two orders of magnitude above the largest of the "gone flat"
figures and four below the smallest "still responding" one, so where exactly
it falls inside that gap does not change which of the measurements above it
classifies. It is deliberately not tied to `_CTSEM_SATURATION` (the retired
raw-magnitude threshold, `binary_measurement.jl`): that constant still guards
an unrelated decision (whether `ctsem_optimize`'s Hager-Zhang fallback is
worth attempting) and lowering it would not fix what was wrong with using a
raw magnitude for saturation in the first place -- see git history for the
false positives that motivated this.
"""
const _CTSEM_TRANSFORM_FLOOR = Ref(1e-6)

"""
    _ctsem_saturated_parameters(sp, values, range; threshold=_CTSEM_TRANSFORM_FLOOR[])

Which raw parameter indices in `range` have a materialising transform that has
stopped responding at `values`.

Each `sp.regular_transforms[tf_idx]` reads exactly one entry of `values` --
`sp.parnumber[tf_idx]` -- by construction on the R side, the same fact
`adjoint_parameters.jl` relies on (and verifies once, at adjoint-workspace
construction) to pull a cotangent back through this layer. So the sensitivity
of each materialised cell to its raw coordinate is one scalar derivative,
taken with a seeded `ForwardDiff.Dual` the same way `_ctsem_regular_pullback!`
does, rather than an `nmut x nvalues` Jacobian over the whole vector -- and
because it is evaluated at the fit's own `values`, not characterised for the
transform in the abstract, a raw parameter whose materialised scale itself
depends on other parameters or predictors is handled for free.

A raw index in `range` that no regular transform names -- a TI-predictor
coefficient, which enters only as a linear multiplier in
`_materialize_subject_values!`, or any other position this layer does not
materialise -- is left out of the result: there is no transform here for it
to have gone flat in, so it cannot saturate by this mechanism.

If a raw parameter feeds more than one materialised cell, the largest of
their derivatives decides: one cell still responding means the raw coordinate
still does something to the likelihood.
"""
function _ctsem_saturated_parameters(sp::EKFParameters, values::AbstractVector{T},
        range; threshold::Real=_CTSEM_TRANSFORM_FLOOR[]) where {T<:Real}
    isempty(range) && return Int[]
    D = ForwardDiff.Dual{Nothing,T,1}
    scratch = Vector{D}(undef, length(values))
    @inbounds for i in eachindex(values)
        scratch[i] = D(values[i], ForwardDiff.Partials((zero(T),)))
    end
    best = Dict{Int,T}()
    tf_idx = 0
    @inbounds for idx in eachindex(sp.mutables)
        sp.mutables[idx] || continue
        tf_idx += 1
        pn = sp.parnumber[tf_idx]
        pn in range || continue
        base = values[pn]
        scratch[pn] = _seed_dual(scratch[pn], base, true)
        derivative = abs(_partial1(sp.regular_transforms[tf_idx](scratch)))
        scratch[pn] = _seed_dual(scratch[pn], base, false)
        best[pn] = haskey(best, pn) ? max(best[pn], derivative) : derivative
    end
    flagged = [pn for (pn, derivative) in best if derivative < threshold]
    sort!(flagged)
    return flagged
end
