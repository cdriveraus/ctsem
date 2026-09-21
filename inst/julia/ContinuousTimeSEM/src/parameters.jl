using ComponentArrays

export EKFParameters
"""
    UNSET_PARAMETER

Fill value for a parameter slot that nothing has written yet.

The full parameter vector is materialised at *mutable* and at *fixed*
positions only, so a model-matrix cell owned solely by a state-dependent
transform is covered by neither: it holds whatever the buffer was filled
with until its transform group runs. Reading it before then is always a
bug, and the fill value decides whether that bug is visible.

Zero is the worst possible choice, because it is what a well-behaved
model matrix is mostly made of: it annihilates the products it enters and
reads back as a plausible fixed value. The row-1 defect fixed in
`kalman_filters.jl` survived precisely that way -- and one attempt to
reproduce it saw nothing, because the test model happened to fix the
initial latent mean at zero, so multiplying it by the unwritten cell hid
the fault a second time.

This follows the convention ctsem already uses for a missing predictor: a
number no real parameter can take, finite so that ForwardDiff partials
stay clean rather than turning into NaNs that travel silently, and large
enough that any likelihood it reaches is absurd on sight.
"""
const UNSET_PARAMETER = 99999.0

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
    # Discrete time is not a different filter, only a different discretization:
    # DRIFT, CINT and DIFFUSION are already the one-step quantities, so the
    # matrix exponential, the Lyapunov solve and the intercept solve that turn
    # continuous parameters into per-interval ones all collapse to identities.
    continuous_time::Bool
    # Which manifest variables are Gaussian (0), binary (1) or ordinal (2).
    # Last because the constructor takes it last: a field declared in one order
    # and supplied in another put a Bool into this slot and failed with a
    # convert error pointing at neither.
    #
    # Empty for the Gaussian-only models that were all this engine handled
    # before, so the filter can skip the branch rather than test a vector of
    # zeros on every row.
    manifesttype::Vector{Int}
    # Category count per manifest variable, used only by the ordinal ones: an
    # ordinal variable with `K` categories reads the first `K-1` columns of
    # THRESHOLDS and ignores the rest, so variables with different numbers of
    # categories can share one rectangular matrix.
    ncategories::Vector{Int}
    # Asymptote count per manifest variable, read only by the binary ones: 0
    # for a plain two parameter logistic, 1 when the lower asymptote is
    # estimated (three parameter) and 2 when both are (four parameter). It
    # says whether that row's THRESHOLDS cells are its asymptotes, which is a
    # property of the model rather than of the values, so it cannot be
    # recovered by looking at them -- a guessing parameter of zero is an
    # ordinary place for an optimizer to be.
    nasymptotes::Vector{Int}
    # Censoring limits per manifest variable, read only by the censored ones.
    # Known constants rather than parameters: a scale's floor and ceiling are
    # properties of the instrument, not things the data can inform. Empty means
    # no variable is censored, which is every model that does not ask for it.
    censormin::Vector{Float64}
    censormax::Vector{Float64}
    # Which covariance construction the model asked for, the same code the stan
    # path reads from `standata$choleskymats`: 0 for the unconstrained
    # correlation square root, 2 for covmattransform='z'. Declared last and
    # supplied last, for the reason the `manifesttype` comment above gives.
    covmatcode::Int
    # Which augmented states the population covariance occupies, in the order
    # of the RAWPOPVAR matrix rows: `population_indices[j]` is the state that
    # population row `j` describes. Those indices are not disjoint from the
    # model own latents -- an individually varying T0MEANS gets no carrier
    # state, so its own latent state is what its population row refers to.
    #
    # Empty means the model has no separate population covariance and T0VAR is
    # the whole of the initial covariance, which is every model that does not
    # use intoverpop.
    population_indices::Vector{Int}
    # Where that matrix sits in the flat parameter vector, as first and last
    # positions. A name lookup would not do: `pars` is a ComponentVector whose
    # axis is a type parameter, so `pars.RAWPOPVAR` has to compile for every
    # model, and for a model without the block it cannot. `0:-1` when absent.
    population_range::UnitRange{Int}
    # Which construction builds the population block, which need not be the one
    # the rest of the model's covariances use. A reduced-rank population
    # covariance is a factor -- code 1 -- while DIFFUSION and T0VAR keep the
    # model's own setting, and forcing all four to agree would change three
    # matrices to fix one. Defaults to `covmatcode`, so a model that has not
    # asked for anything different is exactly as it was.
    population_covmatcode::Int
    # The state-unit conversion for each population row, applied where the
    # block is placed. An appended carrier is 1, because its T0MEANS uses the
    # identity transform and the consuming cell does the scaling; an
    # individually varying T0MEANS carries its cell multiplier*meanscale,
    # because the carrier is the model latent itself and its covariance has to
    # be in that latent's units. Deliberately NOT folded into the block's own
    # transform: the same conversion is applied at the same point on the stan
    # path, and two implementations of one rule about this matrix is what the
    # previous arrangement cost.
    population_scale::Vector{Float64}
    # How many leading states are genuine dynamics rather than the static
    # coordinates a random effect augments the state with.
    #
    # The continuous form's local affine offset needs `JAx` inverted, and `JAx`
    # is exactly singular on a static coordinate -- zero row and zero column --
    # so that solve runs over this leading block. Restricting it is exact rather
    # than an approximation: order the states [d, s] and `JAx` is block
    # triangular with zero static rows while the offset is genuinely zero there,
    # so the sub-solve equals the full one, and the static-to-dynamic coupling
    # still arrives through the full matrix exponential.
    #
    # This is Stan's `1:nlatent` (`ctModelWriter.R`:
    # `mdivide_left(JAx[1:nlatent,1:nlatent], ...)`), NOT its `derrind`.
    # `derrind` additionally drops a latent that merely has no diffusion of its
    # own and no coupling to anything that has some -- whose offset is *not*
    # zero -- and using it took that state's CINT out of the likelihood
    # altogether: gradient exactly zero, so the optimizer left the parameter at
    # its starting value and the summary reported that as an estimate.
    #
    # 0 means the whole state vector. Declared last and supplied last, for the
    # reason the `manifesttype` comment above gives.
    affine_dim::Int

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
        continuous_time::Bool=true,
        manifesttype=Int[],
        ncategories=Int[],
        nasymptotes=Int[],
        censormin=Float64[],
        censormax=Float64[],
        covmatcode::Int=0,
        population_indices=Int[],
        population_range::UnitRange{Int}=0:-1,
        population_covmatcode::Union{Nothing,Integer}=nothing,
        population_scale=Float64[],
        affine_dim::Int=0,
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
            continuous_time,
            Vector{Int}(manifesttype),
            Vector{Int}(ncategories),
            Vector{Int}(nasymptotes),
            Vector{Float64}(censormin),
            Vector{Float64}(censormax),
            covmatcode,
            Vector{Int}(population_indices),
            population_range,
            population_covmatcode === nothing ? Int(covmatcode) :
                Int(population_covmatcode),
            isempty(population_scale) ?
                ones(Float64, length(population_indices)) :
                Vector{Float64}(population_scale),
            affine_dim,
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
