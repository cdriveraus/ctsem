using ComponentArrays

export EKFParameters

"""
    EKFParameters

Container for parameter metadata used by the continuous-time EKF routines.

The fields describe which flattened parameter positions are mutable or fixed,
which transforms should be applied before filtering, and the `ComponentArray`
axis used to view the flattened parameter vector as named model matrices.
"""
struct EKFParameters{RT,PT,UT,TT,AX,FV}
    mutables::BitVector
    transform_indices::Vector{Bool}
    predict_transforms_indices::Vector{Bool}
    update_transforms_indices::Vector{Bool}
    td_transforms_indices::Vector{Bool}
    regular_transforms::RT
    predict_transforms::PT
    update_transforms::UT
    td_transforms::TT
    parnumber::Vector{Int}
    parameter_axis::AX
    fixed_indices::Vector{Bool}
    fixed_values::FV
    ti_parameter_indices::Vector{Int}
    ti_predictor_indices::Vector{Int}
    ti_coefficient_indices::Vector{Int}
    diffusion_state_indices::Vector{Int}

    # The constructor ensures that the provided vectors are of the correct types and converts them if necessary.
    function EKFParameters(
        mutables,
        transform_indices,
        predict_transforms_indices,
        update_transforms_indices,
        td_transforms_indices,
        regular_transforms,
        predict_transforms,
        update_transforms,
        td_transforms,
        parnumber,
        parameter_axis,
        fixed_indices,
        fixed_values,
        ti_parameter_indices=Int[],
        ti_predictor_indices=Int[],
        ti_coefficient_indices=Int[],
        diffusion_state_indices=Int[],
    )
        regular_transforms_tuple = Tuple(regular_transforms)
        predict_transforms_tuple = Tuple(predict_transforms)
        update_transforms_tuple = Tuple(update_transforms)
        td_transforms_tuple = Tuple(td_transforms)
        fixed_values_vector = _concrete_float_vector(fixed_values)

        # Creation of the EKFParameters struct with the appropriate types for each field
        return new{
            typeof(regular_transforms_tuple),
            typeof(predict_transforms_tuple),
            typeof(update_transforms_tuple),
            typeof(td_transforms_tuple),
            typeof(parameter_axis),
            typeof(fixed_values_vector),
        }(
            BitVector(mutables),
            Vector{Bool}(transform_indices),
            Vector{Bool}(predict_transforms_indices),
            Vector{Bool}(update_transforms_indices),
            Vector{Bool}(td_transforms_indices),
            regular_transforms_tuple,
            predict_transforms_tuple,
            update_transforms_tuple,
            td_transforms_tuple,
            Vector{Int}(parnumber),
            parameter_axis,
            Vector{Bool}(fixed_indices),
            fixed_values_vector,
            Vector{Int}(ti_parameter_indices),
            Vector{Int}(ti_predictor_indices),
            Vector{Int}(ti_coefficient_indices),
            Vector{Int}(diffusion_state_indices),
        )
    end
end

# Compatibility constructor for existing callers without TD transform metadata.
function EKFParameters(mutables, transform_indices, predict_transforms_indices,
    update_transforms_indices, regular_transforms, predict_transforms,
    update_transforms, parnumber, parameter_axis, fixed_indices, fixed_values)
    td_transforms_indices = falses(length(mutables))
    return EKFParameters(mutables, transform_indices, predict_transforms_indices,
        update_transforms_indices, td_transforms_indices, regular_transforms,
        predict_transforms, update_transforms, Function[], parnumber,
        parameter_axis, fixed_indices, fixed_values)
end

function _concrete_float_vector(values)
    isempty(values) && return Float64[]
    T = mapreduce(typeof, promote_type, values)
    T <: AbstractFloat || throw(ArgumentError("fixed_values must contain floating-point values"))
    return Vector{T}(values)
end
