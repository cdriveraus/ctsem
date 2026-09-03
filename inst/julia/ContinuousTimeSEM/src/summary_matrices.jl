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
#     conditional on that state, and `ctsem_state_dependent_cells` names those
#     cells so a caller can say so.
#
#   * The derived matrices follow the generated Stan code's definitions,
#     including its restriction of the diffusion-related ones to the states
#     carrying their own diffusion.
#
# The interface is deliberately flat -- names and dimensions from one call, a
# plain `Matrix{Float64}` of stacked columns from another -- because the R
# bridge marshals nested containers element by element. A whole posterior of
# parameter matrices is one array transfer, not one call per sample.

export ctsem_parameter_layout, ctsem_parameter_matrices, ctsem_state_dependent_cells

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
dimensions, and `n_statedep`, the number of cells whose values are conditional
on the state they were evaluated at.

The state-dependent cells themselves come from `ctsem_state_dependent_cells`
rather than from here, and only that count is reported. That is a bridge
constraint, not a design preference: JuliaConnectoR hangs marshalling a
zero-length vector, and a linear model has no state-dependent cells at all, so
returning those vectors unconditionally would deadlock every caller who has the
simplest kind of model.
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

    return (matrix=matrix, nrow=nrow, ncol=ncol, offset=offset, size=position,
        nlatent=nlatent, nmanifest=nmanifest,
        n_statedep=length(_ctsem_state_dependent_positions(sp)))
end

_ctsem_state_dependent_positions(sp::EKFParameters) =
    sort!(unique(vcat(findall(sp.predict_transforms_indices),
        findall(sp.update_transforms_indices), findall(sp.td_transforms_indices))))

"""
    ctsem_state_dependent_cells(objective)

Name the cells whose values are conditional on the state they were evaluated at,
as `(matrix, row, col)` rather than as flat engine offsets.

Only call this when `ctsem_parameter_layout(objective).n_statedep` is positive;
with no such cells the returned vectors are empty, and the R bridge hangs
marshalling those.
"""
function ctsem_state_dependent_cells(objective::CTSEMObjective)
    sp = objective.params
    names, nrows, ncols, offsets, _ = _ctsem_base_layout(sp)
    sdmat = String[]
    sdrow = Int[]
    sdcol = Int[]
    for flat in _ctsem_state_dependent_positions(sp)
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
    return (matrix=sdmat, row=sdrow, col=sdcol)
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
function _ctsem_asymptotics(DRIFT, DIFFUSIONcov, CINT, dyn, nlatent,
    continuous_time::Bool=true)
    asym_diffusion = zeros(Float64, nlatent, nlatent)
    asym_cint = zeros(Float64, nlatent, 1)
    isempty(dyn) && return (asym_diffusion, asym_cint)

    A = Matrix{Float64}(DRIFT[dyn, dyn])
    Q = Matrix{Float64}(DIFFUSIONcov[dyn, dyn])
    k = length(dyn)
    solved = try
        if continuous_time
            X = zeros(Float64, k, k)
            ntri = (k * (k + 1)) ÷ 2
            ksolve!(X, A, Q, zeros(Float64, ntri, ntri), zeros(Float64, ntri),
                Vector{Int}(undef, ntri))
            all(isfinite, X) ? X : nothing
        else
            # Discrete time solves X = A X A' + Q, i.e.
            # (I - A kron A) vec(X) = vec(Q).
            X = reshape((I - kron(A, A)) \ vec(Q), k, k)
            all(isfinite, X) ? X : nothing
        end
    catch
        nothing
    end
    if solved === nothing
        fill!(asym_diffusion, NaN)
    else
        asym_diffusion[dyn, dyn] .= solved
    end

    intercept = try
        # Continuous: -A x = c. Discrete: (I - A) x = c.
        system = continuous_time ? -A : (I - A)
        value = system \ Vector{Float64}(CINT[dyn, 1])
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

`rows` selects flat positions to return, and exists because of what happens
*after* this function: JuliaConnectoR moves about 1 MB/s, and a caller wanting
two cells of a 5-node quadrature was paying for the whole 82-by-1000 array five
times over -- 3.7s of a 16.7s fit, to keep two percent of what crossed. The
materialization itself costs under 0.02s either way, so this is purely about
what goes back over the bridge. Empty means all of them.
"""
function ctsem_parameter_matrices(objective::CTSEMObjective, values::AbstractMatrix;
    tipreds=Float64[], state=Float64[], time::Real=0.0, dt::Real=0.0, rows=Int[])

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
    # Sentinel-filled for the same reason as the EKF workspace: this path
    # runs all three transform groups, so a cell in none of them and
    # neither mutable nor fixed would otherwise be read as garbage.
    all_params = fill(Float64(UNSET_PARAMETER), nall)

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

        _ctsem_pack_matrices!(view(out, :, column), pars, sp, layout)
    end
    isempty(rows) && return out
    selected = Int.(rows)
    all(r -> 1 <= r <= layout.size, selected) || throw(BoundsError(out, selected))
    return out[selected, :]
end

"""
    _ctsem_pack_matrices!(column, pars, sp, layout)

Write one materialized parameter vector, plus the covariance and asymptotic
matrices derived from it, into one column of the flat layout.

Split out because there are two ways to arrive at a materialized parameter
vector -- transform it from a raw vector, or take the one a subject's filter
pass ended with (see `kalman_trace.jl`) -- and only the first half differs.
"""
function _ctsem_pack_matrices!(column, pars, sp::EKFParameters, layout)
    nall = length(getdata(pars))
    column[1:nall] .= getdata(pars)

    nlatent = layout.nlatent
    dyn = isempty(sp.diffusion_state_indices) ? collect(1:nlatent) :
          sp.diffusion_state_indices
    diffusioncov = zeros(Float64, nlatent, nlatent)
    diffusioncov[dyn, dyn] .= Matrix(sdcovsqrt2cov(pars.DIFFUSION, 0))[dyn, dyn]
    manifestcov = Matrix(sdcovsqrt2cov(pars.MANIFESTVAR, 0))
    t0cov = Matrix(sdcovsqrt2cov(pars.T0VAR, 0))
    asym_diffusion, asym_cint = _ctsem_asymptotics(pars.DRIFT, diffusioncov,
        pars.CINT, dyn, nlatent, sp.continuous_time)

    derived = (diffusioncov, manifestcov, t0cov, asym_diffusion, asym_cint)
    nbase = length(layout.matrix) - length(_CTSEM_DERIVED_MATRICES)
    for (k, block) in enumerate(derived)
        slot = nbase + k
        first = layout.offset[slot] + 1
        column[first:(first+length(block)-1)] .= vec(block)
    end
    return column
end

ctsem_parameter_matrices(objective::CTSEMObjective, values::AbstractVector; kwargs...) =
    ctsem_parameter_matrices(objective, reshape(collect(Float64, values), :, 1); kwargs...)
