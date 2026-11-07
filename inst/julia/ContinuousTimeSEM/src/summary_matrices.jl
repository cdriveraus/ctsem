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
    template = CTSEMPars(zeros(length(sp.mutables)), sp.parameter_axis)
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
    nlatent = size(getproperty(CTSEMPars(zeros(nall), sp.parameter_axis), :DRIFT), 1)
    nmanifest = size(getproperty(CTSEMPars(zeros(nall), sp.parameter_axis), :LAMBDA), 1)

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
    # The population covariance of the individually varying parameters on the
    # raw scale, before the state-unit conversion that places it in T0cov:
    # what Stan reports as `rawpopcov`, and what every population summary
    # reads (R's `.ctBackendPopulations()`). The Laplace route appends its
    # levels' own (`ctsem_parameter_layout(::CTSEMLaplaceObjective)`).
    k = length(sp.population_indices)
    if k > 0
        push!(matrix, "rawpopcov")
        push!(nrow, k)
        push!(ncol, k)
        push!(offset, position)
        position += k * k
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
    _ctsem_settling_states(A, continuous_time)

Which states of the drift `A` settle to an equilibrium whatever they start
from: those whose closure under driving -- the state and every state that
drives it, directly or through others (`A[i, j] != 0`) -- has a stable drift.
The closure is an autonomous system of its own, block triangular once ordered,
so its eigenvalues are those of the blocks it depends on, and a state settles
exactly when each of them does. A state with an auto effect of zero (one in
discrete time) and nothing pulling it back -- a trend, a trait, an integrator
-- does not settle, and
neither does any state it drives; a stable process beside it, or driving it,
still does. The same closure the stationary T0MEANS and T0VAR cells are solved
over (`_ctsem_closure_drift`), here from the values rather than the model.
"""
function _ctsem_settling_states(A::AbstractMatrix, continuous_time::Bool)
    k = size(A, 1)
    settles = fill(false, k)
    reach = fill(false, k)
    stack = Int[]
    for i in 1:k
        fill!(reach, false)
        reach[i] = true
        push!(stack, i)
        while !isempty(stack)
            r = pop!(stack)
            for j in 1:k
                if !reach[j] && A[r, j] != 0
                    reach[j] = true
                    push!(stack, j)
                end
            end
        end
        closure = findall(reach)
        block = A[closure, closure]
        all(isfinite, block) || continue
        ev = try
            eigvals(block)
        catch
            continue
        end
        settles[i] = continuous_time ? all(e -> real(e) < 0, ev) :
            all(e -> abs(e) < 1, ev)
    end
    return settles
end

"""
    _ctsem_noise_reach(A, Q, candidates)

The states among `candidates` the diffusion reaches: those with noise of their
own (`Q[i, i] != 0`) and every state they drive, directly or through others
(`A[i, j] != 0` for `j` driving `i`). Sorted. Closed under driving within the
set, so the noise-driven part of the process restricted to it is a system of
its own: every other state carries none of the noise.
"""
function _ctsem_noise_reach(A::AbstractMatrix, Q::AbstractMatrix, candidates::Vector{Int})
    k = size(A, 1)
    allowed = fill(false, k)
    allowed[candidates] .= true
    reach = fill(false, k)
    stack = Int[]
    for i in candidates
        if Q[i, i] != 0
            reach[i] = true
            push!(stack, i)
        end
    end
    while !isempty(stack)
        r = pop!(stack)
        for i in 1:k
            if allowed[i] && !reach[i] && A[i, r] != 0
                reach[i] = true
                push!(stack, i)
            end
        end
    end
    return findall(reach)
end

"""
    _ctsem_asymptotics(DRIFT, DIFFUSIONcov, CINT, dyn, nlatent, continuous_time, naff)

Asymptotic diffusion covariance and asymptotic intercept: the mean and the
diffusion's covariance each latent settles to whatever its start.

Both are over the latents of genuine dynamics (`1:naff`; the carrier states of
augmented random effects past them stay zero). The intercept is solved over
the latents that settle (`_ctsem_settling_states`) and is `NaN` for the rest:
a latent with a fixed auto effect of zero, or an unstable drift, or one driven
by such a latent, has no equilibrium mean, and solving with it anyway gives a
number -- a singular solve's error, or `CINT / 1e-6` once that zero is moved
off -- where there is none.

The covariance is zero wherever the diffusion does not reach
(`_ctsem_noise_reach`, within the diffusion states `dyn`), whatever the drift
there. Among the states it reaches it is solved where the reached states
driving each one are stable, and is `NaN` for the others. So a stable process
driven by a deterministic trend has a covariance though its mean has no
equilibrium, and the trend's covariance is zero; a random walk's is `NaN`.

Each solve is over a set closed under driving, so it is the subsystem's own.
A missing equilibrium is a property of the parameter value, not an error in the
caller, so it comes back as `NaN` rather than thrown: one non-stationary
posterior sample should not abort the summary of the other 199.
"""
function _ctsem_asymptotics(DRIFT, DIFFUSIONcov, CINT, dyn, nlatent,
    continuous_time::Bool=true, naff::Int=nlatent)
    asym_diffusion = zeros(Float64, nlatent, nlatent)
    asym_cint = zeros(Float64, nlatent, 1)
    naff >= 1 || return (asym_diffusion, asym_cint)

    A = Matrix{Float64}(DRIFT[1:naff, 1:naff])
    settles = _ctsem_settling_states(A, continuous_time)
    settled = findall(settles)

    asym_cint[1:naff, 1] .= NaN
    if !isempty(settled)
        intercept = try
            # Continuous: -A x = c. Discrete: (I - A) x = c.
            As = A[settled, settled]
            system = continuous_time ? -As : (I - As)
            value = system \ Vector{Float64}(CINT[settled, 1])
            all(isfinite, value) ? value : nothing
        catch
            nothing
        end
        intercept === nothing || (asym_cint[settled, 1] .= intercept)
    end

    isempty(dyn) && return (asym_diffusion, asym_cint)
    # The noise reaches the latents it enters and those they drive, and no
    # others: there the diffusion's covariance is zero whatever the drift.
    noisy = _ctsem_noise_reach(A, DIFFUSIONcov, [i for i in dyn if i <= naff])
    isempty(noisy) && return (asym_diffusion, asym_cint)
    asym_diffusion[noisy, noisy] .= NaN
    # Within that set, closed under driving, a latent's covariance settles
    # where the noisy states driving it are stable: a trend that no noise
    # reaches moves the mean of what it drives and adds nothing to its spread.
    K = noisy[_ctsem_settling_states(A[noisy, noisy], continuous_time)]
    isempty(K) && return (asym_diffusion, asym_cint)
    AK = Matrix{Float64}(DRIFT[K, K])
    Q = Matrix{Float64}(DIFFUSIONcov[K, K])
    k = length(K)
    solved = try
        if continuous_time
            X = zeros(Float64, k, k)
            # The filter's own solve, built per call.
            my_lyap!(X, AK, Q, LyapBuffer(Float64, k))
            all(isfinite, X) ? X : nothing
        else
            # Discrete time solves X = A X A' + Q, i.e.
            # (I - A kron A) vec(X) = vec(Q).
            X = reshape((I - kron(AK, AK)) \ vec(Q), k, k)
            all(isfinite, X) ? X : nothing
        end
    catch
        nothing
    end
    solved === nothing || (asym_diffusion[K, K] .= solved)
    return (asym_diffusion, asym_cint)
end

"""
    ctsem_parameter_matrices(objective, values; tipreds, state, time, dt)

Materialize every model matrix for one or many raw parameter vectors.

`values` is `npar` by `nsamples` -- a whole posterior costs one call rather than
one per sample -- and the returned `Matrix{Float64}` is `size` by `nsamples`,
with `size` and the block offsets given by `ctsem_parameter_layout`.

`state` is the evaluation point for state-dependent cells, given as the dynamic
processes -- the leading entries of the filter's state vector. Anything shorter
than the full state vector is filled out from the column's own `T0MEANS`, so
the carrier states that represent individually varying parameters always follow
the parameter vector being materialised. A full-length vector still overwrites
them, which is what reporting one person's matrices at a state that person's
own random effects produced needs; the R side spells that `filterstate=` and
the evaluation point `state=`, so the two cannot be confused for each other.

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

    # As wide as a subject's predictor row: a cell may read any predictor
    # (`ctx.tipreds[j]`), not only one with an effect. Unsupplied ones are zero,
    # where the population is reported.
    nsubjectti = isempty(objective.subject_objectives) ? 0 :
        _ctsem_tipred_length(first(objective.subject_objectives).tipreds)
    ntipred = max(length(tipreds), nsubjectti,
        maximum(sp.ti_predictor_indices; init=0))
    ti = zeros(Float64, ntipred)
    ti[eachindex(tipreds)] .= Float64.(tipreds)

    predict_indices = findall(sp.predict_transforms_indices)
    update_indices = findall(sp.update_transforms_indices)
    td_indices = findall(sp.td_transforms_indices)
    ntdpred = size(first(objective.subject_objectives).tdpreds, 1)

    out = zeros(Float64, layout.size, size(values, 2))
    # Single stationary elements are placed by the filter's own functions, on
    # a workspace of the summary's.
    elements = _stationary_elements_on(sp.stationary_elements) ?
        _init_continuous_ekf_workspace(Float64, sp) : nothing
    subject_values = Vector{Float64}(undef, size(values, 1))
    # Sentinel-filled for the same reason as the EKF workspace: this path
    # runs all three transform groups, so a cell in none of them and
    # neither mutable nor fixed would otherwise be read as garbage.
    all_params = fill(Float64(UNSET_PARAMETER), nall)

    for column in axes(values, 2)
        raw = Vector{Float64}(view(values, :, column))
        _materialize_subject_values!(subject_values, raw, sp, ti)
        _materialize_all_params!(all_params, subject_values, sp)
        pars = CTSEMPars(all_params, sp.parameter_axis)
        # A stationary model's T0MEANS is the stationary mean, so that is what
        # this column reports and where its state-dependent cells are evaluated
        # by default. The predict group runs first for the reason it does in
        # the filter: a DRIFT or CINT cell written by a transform has no value
        # until it has. Those cells read no state, so the placeholder state is
        # harmless.
        if sp.stationary || elements !== nothing
            first_state = Vector{Float64}(vec(pars.T0MEANS))
            apply_complex_transforms_at_indices!(getdata(pars), predict_indices,
                sp.predict_transforms, CTSEMRowContext(first_state, pars,
                    zeros(Float64, ntdpred), ti, Float64(time), Float64(dt), 1, 1))
            if sp.stationary
                _ctsem_stationary_mean!(pars, sp)
            elseif !_ctsem_stationary_elements_mean!(pars.T0MEANS, elements, pars,
                    sp.stationary_elements)
                pars.T0MEANS[sp.stationary_elements.means] .= NaN
            end
        end

        # `state` names where the real processes are, and nothing else. The
        # entries past them are carriers: an individually varying parameter is
        # represented as a latent state with no dynamics, so a cell reading one
        # *is* that parameter, and its value belongs to the column being
        # materialised rather than to the caller. Taking a whole state vector
        # from outside pinned those cells to whichever parameter vector built
        # that vector, so a posterior swept through here reported one draw's
        # dynamics for every draw -- zero width uncertainty, no error. Building
        # from this column's own T0MEANS and overwriting only what was supplied
        # makes that unreachable, and lets a caller pass the dynamic processes
        # alone, which is the only part it has any business choosing.
        current = Vector{Float64}(vec(pars.T0MEANS))
        if !isempty(state)
            length(state) <= length(current) || throw(DimensionMismatch(
                "state has $(length(state)) entries; the model has $(length(current))"))
            copyto!(current, 1, state, 1, length(state))
        end
        tdzero = zeros(Float64, ntdpred)
        context = CTSEMRowContext(current, pars, tdzero, ti, Float64(time), Float64(dt), 1, 1)
        # All three groups, in filter order, so every Jacobian block is current.
        apply_complex_transforms_at_indices!(getdata(pars), predict_indices,
            sp.predict_transforms, context)
        apply_complex_transforms_at_indices!(getdata(pars), td_indices,
            sp.td_transforms, context)
        apply_complex_transforms_at_indices!(getdata(pars), update_indices,
            sp.update_transforms, context)

        _ctsem_pack_matrices!(view(out, :, column), pars, sp, layout, elements)
    end
    isempty(rows) && return out
    selected = Int.(rows)
    all(r -> 1 <= r <= layout.size, selected) || throw(BoundsError(out, selected))
    return out[selected, :]
end

"""
    _ctsem_stationary_mean!(pars, sp)

Write the stationary mean, `-DRIFT \\ CINT` over the leading block of genuine
dynamics, over that block of `pars.T0MEANS`: the summary's copy of the mean
`_ctsem_stationary!` gives the filter. NaN when DRIFT is singular there.
"""
function _ctsem_stationary_mean!(pars, sp::EKFParameters)
    n = size(pars.DRIFT, 1)
    dyn = isempty(sp.diffusion_state_indices) ? collect(1:n) :
          sp.diffusion_state_indices
    naff = _ctsem_affine_dim(sp, dyn, n)
    A = Matrix{Float64}(pars.DRIFT[1:naff, 1:naff])
    c = Vector{Float64}(vec(pars.CINT)[1:naff])
    mean = try
        value = -(A \ c)
        all(isfinite, value) ? value : fill(NaN, naff)
    catch
        fill(NaN, naff)
    end
    @inbounds for i in 1:naff
        pars.T0MEANS[i, 1] = mean[i]
    end
    return pars
end

"""
    _ctsem_pack_matrices!(column, pars, sp, layout)

Write one materialized parameter vector, plus the covariance and asymptotic
matrices derived from it, into one column of the flat layout.

Split out because there are two ways to arrive at a materialized parameter
vector -- transform it from a raw vector, or take the one a subject's filter
pass ended with (see `kalman_trace.jl`) -- and only the first half differs.
"""
function _ctsem_pack_matrices!(column, pars, sp::EKFParameters, layout,
        elements=nothing)
    nall = length(getdata(pars))
    column[1:nall] .= getdata(pars)

    nlatent = layout.nlatent
    dyn = isempty(sp.diffusion_state_indices) ? collect(1:nlatent) :
          sp.diffusion_state_indices
    # `sp.covmatcode`, not a hardcoded 0: these are the covariances the
    # summaries report, and the filter builds its own with the requested
    # construction (`kalman_filters.jl` passes `ws.sp.covmatcode`). Hardcoding it
    # here meant a fit under 'z' or 'cholesky' was estimated with one
    # construction and reported through another -- T0cov came back 0.083 out
    # with a sign flip, and nothing failed.
    diffusioncov = zeros(Float64, nlatent, nlatent)
    diffusioncov[dyn, dyn] .=
        Matrix(sdcovsqrt2cov(pars.DIFFUSION, sp.covmatcode))[dyn, dyn]
    manifestcov = Matrix(sdcovsqrt2cov(pars.MANIFESTVAR, sp.covmatcode))
    t0cov = Matrix(sdcovsqrt2cov(pars.T0VAR, sp.covmatcode))
    # And the population block, through the same placement the filter uses.
    # Reporting T0VAR alone here was wrong the moment the population
    # covariance became its own matrix: the block came back zero while the
    # filter had it, which is the disagreement this whole separation exists to
    # remove.
    # The block's raw construction is kept as `rawpopcov`.
    rawpopcov = nothing
    if !isempty(sp.population_indices)
        k = length(sp.population_indices)
        buffer = _make_square_buffer(Float64, k)
        _place_population_block!(t0cov, getdata(pars), sp.population_indices,
            sp.population_range, sp.population_covmatcode, sp.population_scale,
            buffer)
        rawpopcov = Matrix{Float64}(buffer.out[1:k, 1:k])
    end
    asym_diffusion, asym_cint = _ctsem_asymptotics(pars.DRIFT, diffusioncov,
        pars.CINT, dyn, nlatent, sp.continuous_time,
        _ctsem_affine_dim(sp, dyn, nlatent))
    # The block `_ctsem_stationary!` places in the filter: the asymptotic
    # covariance over the diffusion states and zero across the rest of the
    # latent processes' rows and columns. NaN where there is none, as it is in
    # `asym_diffusion`.
    if sp.stationary
        naff = _ctsem_affine_dim(sp, dyn, nlatent)
        t0cov[1:naff, :] .= 0.0
        t0cov[:, 1:naff] .= 0.0
        t0cov[dyn, dyn] .= asym_diffusion[dyn, dyn]
    elseif _stationary_elements_on(sp.stationary_elements)
        # The marked blocks, by the filter's own congruence. Its workspace is
        # the caller's when there is one, and otherwise one of this column's.
        ws = elements === nothing ? _init_continuous_ekf_workspace(Float64, sp) :
            elements
        _ctsem_stationary_elements_cov!(t0cov, ws, pars, sp.stationary_elements) ||
            fill!(t0cov, NaN)
    end

    derived = rawpopcov === nothing ?
        (diffusioncov, manifestcov, t0cov, asym_diffusion, asym_cint) :
        (diffusioncov, manifestcov, t0cov, asym_diffusion, asym_cint, rawpopcov)
    nbase = length(_ctsem_base_layout(sp)[1])
    for (k, block) in enumerate(derived)
        slot = nbase + k
        first = layout.offset[slot] + 1
        column[first:(first+length(block)-1)] .= vec(block)
    end
    return column
end

ctsem_parameter_matrices(objective::CTSEMObjective, values::AbstractVector; kwargs...) =
    ctsem_parameter_matrices(objective, reshape(collect(Float64, values), :, 1); kwargs...)
