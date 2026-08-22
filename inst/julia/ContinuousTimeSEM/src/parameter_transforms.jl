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
