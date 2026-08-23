using ComponentArrays
using LinearAlgebra

################################################################################
# Model-implied parameter matrices
################################################################################
#
# The primitive the R side's summary and plot functions are built on. It runs
# inside the engine rather than re-deriving the transforms in R, because the
# engine already materializes every model matrix from the raw parameter vector
# on its way to a log likelihood; asking it for that same materialization is the
# only way to guarantee that what a summary reports is what the likelihood used.
#
# Two choices here are deliberate and surfaced rather than hidden:
#
#   * State-dependent cells have no single value -- they are functions of the
#     latent state. They are evaluated at a caller-supplied state, defaulting to
#     T0MEANS, with TD predictors fixed at zero. For a linear model this is
#     exact (no cell depends on the state); for a nonlinear one it is a value
#     conditional on that state, and `ctsem_parameter_layout` reports which
#     cells those are so a caller can say so.
#
#   * The derived matrices follow the generated Stan code's definitions,
#     including its restriction of the diffusion-related ones to the states
#     carrying their own diffusion.
#
# The interface is deliberately flat -- names and dimensions from one call, a
# plain `Matrix{Float64}` of stacked columns from another -- because the R
# bridge marshals nested containers element by element. A whole posterior of
# parameter matrices is one array transfer, not one call per sample.

export ctsem_parameter_layout, ctsem_parameter_matrices

const _CTSEM_DERIVED_MATRICES = (:DIFFUSIONcov, :MANIFESTcov, :T0cov,
    :asymDIFFUSIONcov, :asymCINT)

"""
    _ctsem_base_layout(sp)

Names, row counts, column counts and flat offsets of the engine's own model
matrices, in the order the parameter table defined them.
"""
function _ctsem_base_layout(sp::EKFParameters)
    template = ComponentVector(zeros(length(sp.mutables)), sp.parameter_axis)
    names = collect(propertynames(template))
    nrows = Int[]
    ncols = Int[]
    offsets = Int[]
    position = 0
    for name in names
        block = getproperty(template, name)
        r, c = size(block, 1), size(block, 2)
        push!(nrows, r)
        push!(ncols, c)
        push!(offsets, position)
        position += r * c
    end
    return (names, nrows, ncols, offsets, position)
end

"""
    ctsem_parameter_layout(objective)

Describe what `ctsem_parameter_matrices` returns.

Returns a named tuple of flat vectors: `matrix`, `nrow`, `ncol` and `offset`
(zero based) for each matrix, the total flat `size`, the latent and manifest
dimensions, and `statedep_matrix`/`statedep_row`/`statedep_col` naming the cells
whose values are conditional on the state they were evaluated at.
"""
function ctsem_parameter_layout(objective::CTSEMObjective)
    sp = objective.params
    names, nrows, ncols, offsets, nall = _ctsem_base_layout(sp)
    nlatent = size(getproperty(ComponentVector(zeros(nall), sp.parameter_axis), :DRIFT), 1)
    nmanifest = size(getproperty(ComponentVector(zeros(nall), sp.parameter_axis), :LAMBDA), 1)

    matrix = String.(names)
    nrow = copy(nrows)
    ncol = copy(ncols)
    offset = copy(offsets)
    position = nall
    for (name, r, c) in zip(_CTSEM_DERIVED_MATRICES,
        (nlatent, nmanifest, nlatent, nlatent, nlatent),
        (nlatent, nmanifest, nlatent, nlatent, 1))
        push!(matrix, String(name))
        push!(nrow, r)
        push!(ncol, c)
        push!(offset, position)
        position += r * c
    end

    statedep = sort!(unique(vcat(findall(sp.predict_transforms_indices),
        findall(sp.update_transforms_indices), findall(sp.td_transforms_indices))))
    sdmat = String[]
    sdrow = Int[]
    sdcol = Int[]
    for flat in statedep
        for k in eachindex(names)
            first = offsets[k] + 1
            last = offsets[k] + nrows[k] * ncols[k]
            (first <= flat <= last) || continue
            local_index = flat - offsets[k] - 1
            push!(sdmat, String(names[k]))
            push!(sdrow, local_index % nrows[k] + 1)
            push!(sdcol, local_index ÷ nrows[k] + 1)
            break
        end
    end

    return (matrix=matrix, nrow=nrow, ncol=ncol, offset=offset, size=position,
        nlatent=nlatent, nmanifest=nmanifest,
        statedep_matrix=sdmat, statedep_row=sdrow, statedep_col=sdcol)
end

"""
    _ctsem_asymptotics(DRIFT, DIFFUSIONcov, CINT, dyn)

Asymptotic diffusion covariance and asymptotic intercept over the states in
`dyn`, zero elsewhere.

A non-stationary or singular drift has no asymptotic form. That is a property of
the parameter value, not an error in the caller, so it comes back as `NaN`
rather than thrown: one non-stationary posterior sample should not abort the
summary of the other 199.
"""
function _ctsem_asymptotics(DRIFT, DIFFUSIONcov, CINT, dyn, nlatent)
    asym_diffusion = zeros(Float64, nlatent, nlatent)
    asym_cint = zeros(Float64, nlatent, 1)
    isempty(dyn) && return (asym_diffusion, asym_cint)

    A = Matrix{Float64}(DRIFT[dyn, dyn])
    Q = Matrix{Float64}(DIFFUSIONcov[dyn, dyn])
    k = length(dyn)
    solved = try
        X = zeros(Float64, k, k)
        ntri = (k * (k + 1)) ÷ 2
        ksolve!(X, A, Q, zeros(Float64, ntri, ntri), zeros(Float64, ntri),
            Vector{Int}(undef, ntri))
        all(isfinite, X) ? X : nothing
    catch
        nothing
    end
    if solved === nothing
        fill!(asym_diffusion, NaN)
    else
        asym_diffusion[dyn, dyn] .= solved
    end

    intercept = try
        value = (-A) \ Vector{Float64}(CINT[dyn, 1])
        all(isfinite, value) ? value : nothing
    catch
        nothing
    end
    if intercept === nothing
        fill!(asym_cint, NaN)
    else
        asym_cint[dyn, 1] .= intercept
    end
    return (asym_diffusion, asym_cint)
end

"""
    ctsem_parameter_matrices(objective, values; tipreds, state, time, dt)

Materialize every model matrix for one or many raw parameter vectors.

`values` is `npar` by `nsamples` -- a whole posterior costs one call rather than
one per sample -- and the returned `Matrix{Float64}` is `size` by `nsamples`,
with `size` and the block offsets given by `ctsem_parameter_layout`.
"""
function ctsem_parameter_matrices(objective::CTSEMObjective, values::AbstractMatrix;
    tipreds=Float64[], state=Float64[], time::Real=0.0, dt::Real=0.0)

    sp = objective.params
    layout = ctsem_parameter_layout(objective)
    names, nrows, ncols, offsets, nall = _ctsem_base_layout(sp)

    ntipred = length(sp.ti_predictor_indices) == 0 ? length(tipreds) :
              max(length(tipreds), maximum(sp.ti_predictor_indices))
    ti = zeros(Float64, ntipred)
    ti[eachindex(tipreds)] .= Float64.(tipreds)

    predict_indices = findall(sp.predict_transforms_indices)
    update_indices = findall(sp.update_transforms_indices)
    td_indices = findall(sp.td_transforms_indices)
    ntdpred = size(first(objective.subject_objectives).tdpreds, 1)

    out = zeros(Float64, layout.size, size(values, 2))
    subject_values = Vector{Float64}(undef, size(values, 1))
    all_params = Vector{Float64}(undef, nall)

    for column in axes(values, 2)
        raw = Vector{Float64}(view(values, :, column))
        _materialize_subject_values!(subject_values, raw, sp, ti)
        _materialize_all_params!(all_params, subject_values, sp)
        pars = ComponentVector(all_params, sp.parameter_axis)

        current = isempty(state) ? Vector{Float64}(vec(pars.T0MEANS)) : Vector{Float64}(state)
        tdzero = zeros(Float64, ntdpred)
        context = CTSEMRowContext(current, pars, tdzero, ti, Float64(time), Float64(dt), 1, 1)
        # All three groups, in filter order, so every Jacobian block is current.
        apply_complex_transforms_at_indices!(getdata(pars), predict_indices,
            sp.predict_transforms, context)
        apply_complex_transforms_at_indices!(getdata(pars), td_indices,
            sp.td_transforms, context)
        apply_complex_transforms_at_indices!(getdata(pars), update_indices,
            sp.update_transforms, context)

        out[1:nall, column] .= getdata(pars)

        nlatent = layout.nlatent
        dyn = isempty(sp.diffusion_state_indices) ? collect(1:nlatent) :
              sp.diffusion_state_indices
        diffusioncov = zeros(Float64, nlatent, nlatent)
        diffusioncov[dyn, dyn] .= Matrix(sdcovsqrt2cov(pars.DIFFUSION, 0))[dyn, dyn]
        manifestcov = Matrix(sdcovsqrt2cov(pars.MANIFESTVAR, 0))
        t0cov = Matrix(sdcovsqrt2cov(pars.T0VAR, 0))
        asym_diffusion, asym_cint = _ctsem_asymptotics(pars.DRIFT, diffusioncov,
            pars.CINT, dyn, nlatent)

        derived = (diffusioncov, manifestcov, t0cov, asym_diffusion, asym_cint)
        for (k, block) in enumerate(derived)
            slot = length(names) + k
            first = layout.offset[slot] + 1
            out[first:(first+length(block)-1), column] .= vec(block)
        end
    end
    return out
end

ctsem_parameter_matrices(objective::CTSEMObjective, values::AbstractVector; kwargs...) =
    ctsem_parameter_matrices(objective, reshape(collect(Float64, values), :, 1); kwargs...)
