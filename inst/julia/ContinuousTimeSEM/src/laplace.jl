"""
Laplace-approximate marginal likelihood for subject-level parameter random
effects.

# What this replaces

ctsem's existing route for individual differences (`intoverpop`) *augments the
latent state*: each varying parameter becomes an extra, static latent state
carrying its own population variance, and the ordinary Kalman filter integrates
it out along with the dynamic states. That is exact for a linear-Gaussian model
and needs no new machinery, but every subject then pays for a state space of
dimension `nlatent + nindvarying` -- and the filter's cost is cubic in that
dimension, in the Lyapunov solve, the matrix exponential and its Frechet
derivative alike. Twelve random effects on a four-latent model is a sixteen-
dimensional system per row, for every subject, at every evaluation.

This file takes the other route. Each subject keeps the *small* `nlatent`
system, and the random effects are integrated out per subject by a Laplace
approximation over a `k = nindvarying` dimensional integral. The per-subject
inner problems are completely independent, so the curvature that has to be
factorized is `nsubjects` separate `k x k` blocks rather than one large system.

# The model

Exactly the one the generated Stan model already states for `intoverpop == 0`
(`ctModelWriter.R`, `rawindparams[indvaryingindex] += rawpopcovchol *
baseindparams[si]`):

    raw_i = theta + scatter(L * z_i) + TI-predictor effects,   z_i ~ N(0, I)

with `L` the Cholesky factor of the raw-scale population covariance, built from
`theta` by `_laplace_popchol` in the same parameterisation Stan uses. The
difference is only in how `z_i` is handled: Stan samples it, this integrates it
out. TI-predictor effects are untouched and stay where they are, inside the
engine's own per-subject parameter layer -- both contributions are additive on
the raw scale, so they compose without either knowing about the other.

# The approximation

Write the inner objective for subject `i`, with the standard-normal density's
own exponent folded in:

    g_i(z) = ll_i(theta, z) - z'z / 2

Then, because the `(2*pi)^(-k/2)` in `N(z|0,I)` cancels the `(2*pi)^(k/2)` the
Laplace approximation produces,

    log integral_i  =  g_i(zhat_i) - logdet(-H_i) / 2

with `zhat_i` the inner mode and `H_i` the inner curvature there. No constants
survive; that cancellation is asserted in `test_laplace.jl` against a
deliberately Gaussian integrand, where Laplace is exact.

# The outer gradient

`H_i` depends on `theta` both directly and through `zhat_i(theta)`, so an exact
outer gradient needs third derivatives of the process likelihood. That is the
requirement that stopped this feature's earlier Stan-based attempt: Stan Math's
higher-order primitives would not instantiate for the generated ctsem model.
Here it is simply `ForwardDiff` over the engine's existing forward-over-reverse
Hessian, which nests because every layer below is generic in its scalar type.

`zhat_i(theta)`'s own derivative comes from one Newton step taken in dual
arithmetic from the converged primal mode (`_laplace_dual_mode`). Since the
inner gradient's primal part is zero at the mode, that step contributes nothing
to the primal and exactly `-H^-1 dg/dtheta` to the dual -- which is the implicit
function theorem, obtained without ever forming `dH/dtheta` by hand.

"""

using LinearAlgebra
using ForwardDiff

export CTSEMLaplaceSpec, CTSEMLaplaceObjective, ctsem_laplace_objective,
    ctsem_laplace_evaluate, ctsem_laplace_modes, ctsem_laplace_optimize,
    ctsem_laplace_popcov

"""
    CTSEMLaplaceLevel(re_index, sd_index, cor_index, sd_scale, group, ngroups)

One level of the hierarchy: which raw parameters vary at it, which raw
parameters say how much, and which group each subject belongs to.

The subject level is the case `group == 1:nsubjects`; a study level assigns
several subjects the same group. Levels are ordered innermost first, so
`levels[1]` is the subject level whenever there is one.

All the index fields point into the *same* flat raw parameter vector the rest
of the engine works in; the population covariance parameters live in a tail of
that vector which no model matrix cell reads. That keeps one contiguous
parameter vector for the outer optimizer, the prior term and the Hessian, with
no packing anywhere.

  * `re_index[j]`: raw position of the `j`-th parameter varying at this level.
  * `sd_index[j]`: raw position of its population scale, before transformation.
  * `cor_index`: raw positions of the unconstrained correlation parameters, in
    column-major lower-triangular order -- the order Stan's own counter walks.
  * `sd_scale[j]`: the model's `sdscale` multiplier for that parameter.
  * `group[i]`: which group at this level subject `i` belongs to.
"""
struct CTSEMLaplaceLevel
    re_index::Vector{Int}
    sd_index::Vector{Int}
    cor_index::Vector{Int}
    sd_scale::Vector{Float64}
    group::Vector{Int}
    ngroups::Int

    function CTSEMLaplaceLevel(re_index, sd_index, cor_index, sd_scale, group,
        ngroups::Integer)
        re = Vector{Int}(collect(re_index))
        sd = Vector{Int}(collect(sd_index))
        cor = Vector{Int}(collect(cor_index))
        scale = Vector{Float64}(collect(sd_scale))
        grp = Vector{Int}(collect(group))
        k = length(re)
        length(sd) == k ||
            throw(DimensionMismatch("one population scale parameter per random effect is required"))
        length(scale) == k ||
            throw(DimensionMismatch("one sdscale per random effect is required"))
        expected = div(k * (k - 1), 2)
        length(cor) == expected || throw(DimensionMismatch(
            "expected $(expected) correlation parameters for $(k) random effects, got $(length(cor))"))
        allunique(re) || throw(ArgumentError("random-effect parameter indices must be distinct within a level"))
        isempty(grp) || (minimum(grp) >= 1 && maximum(grp) <= ngroups) ||
            throw(ArgumentError("group ids must lie in 1:ngroups"))
        return new(re, sd, cor, scale, grp, Int(ngroups))
    end
end

"""Number of random effects carried by one level."""
nrandomeffects(level::CTSEMLaplaceLevel) = length(level.re_index)

"""
    CTSEMLaplaceSpec(levels)

The whole hierarchy: one `CTSEMLaplaceLevel` per level, innermost first.

A single-level spec is the ordinary subject-random-effects case and is what the
four-vector constructor below builds.
"""
struct CTSEMLaplaceSpec
    levels::Vector{CTSEMLaplaceLevel}
end

"""Single-level convenience form: subject-level effects, one group per subject."""
function CTSEMLaplaceSpec(re_index, sd_index, cor_index, sd_scale,
    nsubjects::Integer=0)
    group = collect(1:Int(nsubjects))
    return CTSEMLaplaceSpec([CTSEMLaplaceLevel(re_index, sd_index, cor_index,
        sd_scale, group, Int(nsubjects))])
end

"""Total number of random effects across every level, per subject."""
nrandomeffects(spec::CTSEMLaplaceSpec) = sum(nrandomeffects(l) for l in spec.levels; init=0)

nlevels(spec::CTSEMLaplaceSpec) = length(spec.levels)

"""
    CTSEMLaplaceUnits

How the integral factorises.

Subjects sharing an outer random effect cannot be integrated separately: the
effect couples them. The *unit* is therefore the group at the outermost level
-- a study, when subjects are nested in studies -- and the latent vector `u` a
unit integrates over stacks one block per (level, group) it contains: the
study's own effects, and each of its subjects' effects.

With one level this degenerates to one subject per unit and `u == z`, which is
exactly the two-level-free case and costs nothing extra.

  * `members[U]`: subject indices in unit `U`.
  * `offsets[U][m][l]`: where level `l`'s block for `members[U][m]` starts in
    `u`, zero-based. Members of the same study share the study-level offset,
    which is precisely how the coupling is expressed.
  * `dims[U]`: length of `u` for unit `U`.
"""
struct CTSEMLaplaceUnits
    members::Vector{Vector{Int}}
    offsets::Vector{Vector{Vector{Int}}}
    dims::Vector{Int}
    # The blocks of `u`, and which members each one belongs to. A subject block
    # belongs to one member; a study block to every member of the study. This
    # is what makes the curvature sparse and, more importantly, what makes
    # *forming* it cheap -- see `_laplace_unit_hessian`.
    #
    # Each entry is `(offset, size, member positions)`, positions indexing into
    # `members[U]`.
    blocks::Vector{Vector{Tuple{Int,Int,Vector{Int}}}}
end

"""
    _laplace_build_units(spec, nsubjects)

Group subjects into integration units and lay out each unit's latent vector.

Strict nesting is assumed and checked by the R side; here the outermost level's
grouping simply defines the units, and every inner level's groups are placed
inside whichever unit its subjects fall in.
"""
function _laplace_build_units(spec::CTSEMLaplaceSpec, nsubjects::Integer)
    levels = spec.levels
    nlev = length(levels)
    if nlev == 0 || nsubjects == 0
        return CTSEMLaplaceUnits([[i] for i in 1:nsubjects],
            [[Int[] for _ in 1:1] for _ in 1:nsubjects], zeros(Int, nsubjects),
            [Tuple{Int,Int,Vector{Int}}[] for _ in 1:nsubjects])
    end
    outer = levels[end]
    nunits = outer.ngroups
    members = [Int[] for _ in 1:nunits]
    for i in 1:nsubjects
        push!(members[outer.group[i]], i)
    end
    offsets = Vector{Vector{Vector{Int}}}(undef, nunits)
    dims = zeros(Int, nunits)
    for U in 1:nunits
        cursor = 0
        seen = Dict{Tuple{Int,Int},Int}()
        offsets[U] = [zeros(Int, nlev) for _ in eachindex(members[U])]
        for (m, i) in enumerate(members[U])
            for l in 1:nlev
                k = nrandomeffects(levels[l])
                key = (l, levels[l].group[i])
                slot = get(seen, key, -1)
                if slot < 0
                    slot = cursor
                    seen[key] = slot
                    cursor += k
                end
                offsets[U][m][l] = slot
            end
        end
        dims[U] = cursor
    end

    # Invert the offset map: which member positions share each block.
    blocks = Vector{Vector{Tuple{Int,Int,Vector{Int}}}}(undef, nunits)
    for U in 1:nunits
        owners = Dict{Int,Vector{Int}}()
        sizes = Dict{Int,Int}()
        for (m, i) in enumerate(members[U])
            for l in 1:nlev
                k = nrandomeffects(levels[l])
                k == 0 && continue
                off = offsets[U][m][l]
                push!(get!(owners, off, Int[]), m)
                sizes[off] = k
            end
        end
        blocks[U] = [(off, sizes[off], sort(unique(owners[off])))
                     for off in sort(collect(keys(owners)))]
    end
    return CTSEMLaplaceUnits(members, offsets, dims, blocks)
end

"""
    CTSEMLaplaceObjective(objective, spec; inner_maxiter, inner_tol)

A `CTSEMObjective` plus the hierarchy to integrate out of it.

The inner modes are *state*, not output: they are retained between calls and
warm-start the next evaluation's Newton solve. Across an outer optimizer's
trajectory consecutive parameter vectors are close, so the inner solve usually
converges in one or two steps after the first evaluation.
"""
mutable struct CTSEMLaplaceObjective{O}
    objective::O
    spec::CTSEMLaplaceSpec
    units::CTSEMLaplaceUnits
    # One mode vector per unit; ragged, because units differ in size whenever
    # studies do.
    modes::Vector{Vector{Float64}}
    inner_maxiter::Int
    inner_tol::Float64
    # Adjoint workspaces: one dictionary per chunk of the unit loop, each
    # mapping a scalar type to its workspace. The engine's own cache lives on
    # the CTSEMObjective and holds one type at a time, which would thrash badly
    # here -- an exact outer gradient uses Float64 and two nested dual types
    # within a single evaluation.
    #
    # A dictionary *per chunk* rather than one shared dictionary keyed by chunk:
    # the workspaces were already private, but the container was not, and
    # concurrent inserts into one `Dict` corrupt it. Julia catches that
    # ("Multiple concurrent writes to Dict detected!") rather than silently
    # returning wrong numbers, which is how this was found.
    workspaces::Vector{Dict{Any,Any}}
    # Per unit, filled by the last evaluation; see `ctsem_laplace_diagnostics`.
    inner_iterations::Vector{Int}
    inner_gradient::Vector{Float64}
    inner_converged::Vector{Bool}
    hessian_repaired::Vector{Bool}
end

function CTSEMLaplaceObjective(objective::CTSEMObjective, spec::CTSEMLaplaceSpec;
    inner_maxiter::Integer=50, inner_tol::Real=1e-10)
    nsubjects = length(objective.subject_objectives)
    units = _laplace_build_units(spec, nsubjects)
    nunits = length(units.members)
    return CTSEMLaplaceObjective{typeof(objective)}(objective, spec, units,
        [zeros(Float64, units.dims[U]) for U in 1:nunits],
        Int(inner_maxiter), Float64(inner_tol), [Dict{Any,Any}()],
        zeros(Int, nunits), zeros(Float64, nunits), falses(nunits), falses(nunits))
end

"""
    ctsem_laplace_objective(objective; re_index, sd_index, cor_index, sd_scale, ...)

Build a `CTSEMLaplaceObjective` from flat index vectors, which is the form the
R side sends across the bridge.

The index vectors are *keyword* arguments with empty defaults because of that
bridge: JuliaConnectoR deadlocks marshalling a zero-length vector, so R must be
able to omit one rather than send it. `cor_index` is empty for exactly one real
model shape -- a single random effect at a level, which has no correlations --
and that shape is common enough that it cannot be an error.

For more than one level, the vectors are concatenated innermost level first and
split by `level_nre` (random effects per level). `group` is likewise
concatenated, `nsubjects` entries per level, and `level_ngroups` says how many
groups each level has. With one level all of that collapses to the
single-level form and none of it needs sending.
"""
function ctsem_laplace_objective(objective::CTSEMObjective; re_index=Int[],
    sd_index=Int[], cor_index=Int[], sd_scale=Float64[], level_nre=Int[],
    group=Int[], level_ngroups=Int[], inner_maxiter::Integer=50,
    inner_tol::Real=1e-10)
    nsubjects = length(objective.subject_objectives)
    counts = isempty(level_nre) ? [length(re_index)] : Vector{Int}(Int.(level_nre))
    ngroups = isempty(level_ngroups) ? [nsubjects] : Vector{Int}(Int.(level_ngroups))
    length(counts) == length(ngroups) || throw(DimensionMismatch(
        "level_nre and level_ngroups must describe the same number of levels"))
    groups = isempty(group) ? collect(1:nsubjects) : Vector{Int}(Int.(group))
    length(groups) == nsubjects * length(counts) || throw(DimensionMismatch(
        "group must hold one entry per subject per level"))

    levels = CTSEMLaplaceLevel[]
    re_at = 0; cor_at = 0
    for l in eachindex(counts)
        k = counts[l]
        ncor = div(k * (k - 1), 2)
        push!(levels, CTSEMLaplaceLevel(
            Int.(re_index[(re_at + 1):(re_at + k)]),
            Int.(sd_index[(re_at + 1):(re_at + k)]),
            Int.(cor_index[(cor_at + 1):(cor_at + ncor)]),
            Float64.(sd_scale[(re_at + 1):(re_at + k)]),
            groups[((l - 1) * nsubjects + 1):(l * nsubjects)],
            ngroups[l]))
        re_at += k; cor_at += ncor
    end
    return CTSEMLaplaceObjective(objective, CTSEMLaplaceSpec(levels);
        inner_maxiter=inner_maxiter, inner_tol=inner_tol)
end

"""Positional convenience form, for single-level calls written in Julia."""
ctsem_laplace_objective(objective::CTSEMObjective, re_index, sd_index, cor_index,
    sd_scale; kwargs...) =
    ctsem_laplace_objective(objective; re_index=re_index, sd_index=sd_index,
        cor_index=cor_index, sd_scale=sd_scale, kwargs...)

################################################################################
# The population covariance
################################################################################

"""
Largest unconstrained correlation coordinate allowed, and the correlation it
maps to: `2/(1+exp(-p)) - 1`, so `p = 5.2933` is a correlation of `0.99`.

A random-effect correlation can be driven to the boundary by a design that
cannot identify it -- a level estimates `k + k(k-1)/2` parameters from as many
draws as it has groups, and there are usually far fewer studies than subjects.
Left alone the coordinate runs off to infinity, the covariance approaches
singularity, and the curvature factorizations start failing for reasons that
have nothing to do with the model.

Capping keeps the fit alive and bounded. It is a real restriction of the
parameter space, not a numerical detail, so `ctsem_laplace_boundary` reports
which coordinates are sitting on it and the R side warns naming them. A prior
would also solve this, but it would change the estimator silently; a cap that
announces itself does not.
"""
const _LAPLACE_COR_CAP = Ref(5.2933)

export ctsem_set_correlation_cap!, ctsem_laplace_boundary
"""
    ctsem_set_correlation_cap!(p)

Set the largest unconstrained correlation coordinate. `0` removes the cap.
"""
function ctsem_set_correlation_cap!(p::Real)
    p >= 0 || throw(ArgumentError("correlation cap must be non-negative"))
    _LAPLACE_COR_CAP[] = Float64(p)
    return Float64(p)
end

"""Clamp one correlation coordinate, leaving other scalar types alone."""
@inline function _laplace_cap_correlation(x)
    cap = _LAPLACE_COR_CAP[]
    cap <= 0 && return x
    return clamp(x, -cap, cap)
end

"""
    ctsem_laplace_boundary(laplace, values)

Which unconstrained correlation coordinates are sitting on the cap, as
`(level, position, value)` triples. Empty when none are, which is the ordinary
case.
"""
function ctsem_laplace_boundary(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    cap = _LAPLACE_COR_CAP[]
    levels = Int[]; positions = Int[]; found = Float64[]
    cap <= 0 && return (level=levels, position=positions, value=found)
    theta = collect(Float64, values)
    for (l, level) in enumerate(laplace.spec.levels)
        for (t, idx) in enumerate(level.cor_index)
            if abs(theta[idx]) >= cap - 1e-8
                push!(levels, l); push!(positions, t); push!(found, theta[idx])
            end
        end
    end
    return (level=levels, position=positions, value=found)
end

"""
    _laplace_check_indices(laplace, npar)

Refuse a parameter vector too short to hold every index the spec references.

The R side checks the same invariant when it builds the specification, so this
should be unreachable; it exists because the failure it replaces is a
`BoundsError` raised deep inside a filter, which the optimizer's invalid-point
guard turns into a silently wrong fit rather than an error.
"""
function _laplace_check_indices(laplace::CTSEMLaplaceObjective, npar::Integer)
    for (l, level) in enumerate(laplace.spec.levels)
        for (what, index) in (("varying parameter", level.re_index),
                              ("population scale", level.sd_index),
                              ("correlation", level.cor_index))
            isempty(index) && continue
            maximum(index) <= npar || throw(ArgumentError(string(
                "level ", l, " references ", what, " ", maximum(index),
                " but the parameter vector holds ", npar)))
        end
    end
    return nothing
end

"""
    _laplace_popchol(values, level)

The Cholesky factor of one level's raw-scale population covariance.

Term for term the generated Stan model's construction (`ctModelWriter.R`
around `rawpopsd = ...` through `rawpopcovchol = cholesky_decompose(...)`),
including both of its numerical offsets and its exact correlation
parameterisation, so that a Laplace fit and a Stan fit of the same model are
describing the same population distribution rather than two similar ones.

`constraincorsqrt1` reads only the off-diagonal entries of its argument, so the
scales sitting on the diagonal of `base` are inert there; they are written in
anyway to keep the correspondence with Stan's `rawpopcovbase` literal.
"""
function _laplace_popchol(values::AbstractVector{T}, level::CTSEMLaplaceLevel) where {T}
    k = nrandomeffects(level)
    k == 0 && return zeros(T, 0, 0)
    scales = Vector{T}(undef, k)
    @inbounds for j in 1:k
        raw = values[level.sd_index[j]]
        scales[j] = log1p_exp(2 * raw - 1) * level.sd_scale[j] + 1e-10
    end
    base = zeros(T, k, k)
    counter = 0
    @inbounds for j in 1:k
        base[j, j] = scales[j]
        for i in 1:k
            if i > j
                counter += 1
                base[i, j] = 2 / (1 + exp(-_laplace_cap_correlation(
                    values[level.cor_index[counter]]))) - 1
            end
        end
    end
    corsqrt = constraincorsqrt1(base)
    correlation = corsqrt * transpose(corsqrt)
    scaled = scales .+ 1e-8
    covariance = (scaled .* correlation) .* transpose(scaled)
    symmetric = (covariance .+ transpose(covariance)) ./ 2
    return Matrix(cholesky(Symmetric(symmetric)).L)
end

"""Every level's Cholesky factor, innermost first."""
_laplace_popchols(values::AbstractVector{T}, spec::CTSEMLaplaceSpec) where {T} =
    [_laplace_popchol(values, level) for level in spec.levels]

"""Single-level shorthand, kept for the seeded gradient path."""
_laplace_popchol(values::AbstractVector, spec::CTSEMLaplaceSpec) =
    _laplace_popchol(values, spec.levels[1])

"""
    ctsem_laplace_popcov(laplace, values, level=1)

The raw-scale population covariance matrix implied by `values` at one level.
Reporting helper; the fit itself only ever needs the Cholesky factor.
"""
function ctsem_laplace_popcov(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    level::Integer=1)
    L = _laplace_popchol(collect(Float64, values), laplace.spec.levels[level])
    return L * transpose(L)
end

"""
    _laplace_member_values(values, spec, Ls, u, offsets)

The raw parameter vector one subject is filtered with: the population vector,
shifted at every level by that level's effects for the group this subject
belongs to.

The engine's own per-subject parameter layer then adds TI-predictor effects on
top, exactly as it does without random effects. Every shift is additive on the
raw scale, so none of them has to know about the others.
"""
function _laplace_member_values(values::AbstractVector{T}, spec::CTSEMLaplaceSpec,
    Ls::Vector{<:AbstractMatrix}, u::AbstractVector, offsets::Vector{Int}) where {T}
    S = promote_type(T, eltype(u), eltype(eltype(Ls)))
    shifted = Vector{S}(undef, length(values))
    copyto!(shifted, values)
    @inbounds for l in eachindex(spec.levels)
        level = spec.levels[l]
        k = nrandomeffects(level)
        k == 0 && continue
        base = offsets[l]
        L = Ls[l]
        for p in 1:k
            acc = zero(S)
            for q in 1:k
                acc += L[p, q] * u[base + q]
            end
            shifted[level.re_index[p]] += acc
        end
    end
    return shifted
end

"""
    _laplace_subject_values(values, spec, L, z)

Single-level shorthand: the population vector shifted by one subject's own
effects. Retained because the seeded gradient path is specialised to one level
and reads more clearly in those terms.
"""
function _laplace_subject_values(values::AbstractVector{T}, spec::CTSEMLaplaceSpec,
    L::AbstractMatrix, z::AbstractVector) where {T}
    S = promote_type(T, eltype(L), eltype(z))
    shifted = Vector{S}(undef, length(values))
    copyto!(shifted, values)
    level = spec.levels[1]
    offset = L * z
    @inbounds for j in eachindex(level.re_index)
        shifted[level.re_index[j]] += offset[j]
    end
    return shifted
end

"""
    _laplace_workspace!(laplace, ::Type{T}, nvalues, slot)

An adjoint workspace for scalar type `T`, cached on the Laplace object.

The engine's own cache holds a single type and rebuilds when it changes, which
is right for its callers (a Hessian call, then ordinary gradients again) and
wrong for this one: an exact outer gradient uses `Float64` for the primal inner
solve and nested duals for the outer sweep, interleaved per subject. Caching
per type turns that thrash into one build per type per model.
"""
function _laplace_workspace!(laplace::CTSEMLaplaceObjective, ::Type{T},
    nvalues::Integer, slot::Integer=1) where {T}
    # `slot` is the chunk index, not the thread id. A task can migrate between
    # threads at any yield point, so thread-indexed mutable scratch is a race
    # rather than an optimisation -- the same reason `_get_or_init_adjoint_
    # workspaces!` hands one workspace to each chunk.
    store = laplace.workspaces[slot]
    key = (T, Int(nvalues))
    cached = get(store, key, nothing)
    cached === nothing || return cached
    objective = laplace.objective
    ntdpred = isempty(objective.subject_objectives) ? 0 :
        size(first(objective.subject_objectives).tdpreds, 1)
    built = CTSEMAdjointWorkspace(T, objective.params, Int(nvalues), ntdpred)
    store[key] = built
    return built
end

"""
    _laplace_subject_value_gradient!(gradient, subject_objective, aws, values)

Subject `i`'s log likelihood at `values`, and its gradient with respect to
`values`, written into `gradient`.

This is one iteration of `_ctsem_subject_gradient_chunk!`'s loop, lifted out
because the Laplace layer needs each subject evaluated at a *different*
parameter vector -- which is exactly what the summed adjoint is built never to
have to do. The two shortcuts that make the summed gradient fast, the shared
parameter layer and the cross-subject Frechet batch, are per-subject-incorrect
for the same reason they are in `ctsem_subject_gradients`, so the Frechet batch
is flushed within the subject here too.
"""
function _laplace_subject_value_gradient!(gradient::AbstractVector{T},
    subject_objective, aws, values::AbstractVector{T}) where {T}
    ws = _get_or_init_objective_workspace!(subject_objective, T)
    tape = _tape_reset!(aws.tape)
    resize!(aws.tipreds, length(subject_objective.tipreds))
    copyto!(aws.tipreds, subject_objective.tipreds)
    aws.frechet_pending = false
    deferred = aws.defer_frechet
    aws.defer_frechet = false
    try
        loglik = _extended_kalman_filter_continuous!(ws, values,
            subject_objective.data, subject_objective.timesteps,
            subject_objective.params, subject_objective.tdpreds,
            subject_objective.tipreds, subject_objective.subject,
            subject_objective.max_timestep, tape)
        isfinite(loglik) || return loglik
        fill!(aws.theta_bar, zero(T))
        _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
        fill!(gradient, zero(T))
        _ctsem_parameter_layer!(gradient, aws.theta_bar, tape.subject_values,
            subject_objective.params, aws, subject_objective.tipreds)
        return loglik
    finally
        aws.defer_frechet = deferred
    end
end

"""
    _laplace_negate_definite(H)

`(-H, was_negative_definite)` with `-H` made positive definite if it was not.

`H` is the inner curvature, which the approximation requires to be negative
definite: `-H` is the precision matrix of the Gaussian being fitted to the
integrand, and its determinant is the approximation's normalizing constant. A
non-negative-definite `H` means the point is not a maximum, which happens on
the way to one and can happen at one for a badly identified model. Shifting the
diagonal is the standard repair; reporting that it happened is what keeps it
from being a silent change of objective.
"""
function _laplace_negate_definite(H::AbstractMatrix{T}) where {T}
    negated = -(H .+ transpose(H)) ./ 2
    size(negated, 1) == 0 && return (true, negated)
    factorization = cholesky(Symmetric(negated); check=false)
    issuccess(factorization) && return (true, negated)
    shift = sqrt(eps(real(float(one(T)))))
    scale = maximum(abs, diag(negated))
    scale = isfinite(scale) && scale > 0 ? scale : one(real(float(one(T))))
    for _ in 1:30
        candidate = negated + (shift * scale) * I
        if issuccess(cholesky(Symmetric(candidate); check=false))
            return (false, candidate)
        end
        shift *= 10
    end
    return (false, negated + (scale + one(scale)) * I)
end

################################################################################
# One unit, one latent vector
################################################################################

"""
    _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws, positions)

The summed *process* log likelihood of the given member positions and its
gradient with respect to `u` -- without the `-u'u/2` term, which belongs to the
unit as a whole rather than to any member.

Restricting the member set is what makes the blocked curvature below cheap. It
is exact, not an approximation: a member outside the set has no dependence on
the block being differentiated, so its second derivative with respect to that
block is structurally zero rather than merely small.
"""
function _laplace_unit_loglik_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws, positions) where {T}
    spec = laplace.spec
    units = laplace.units
    members = units.members[U]
    inner = zeros(T, length(u))
    total = zero(T)
    gradient = Vector{T}(undef, length(values))
    @inbounds for m in positions
        i = members[m]
        offsets = units.offsets[U][m]
        shifted = _laplace_member_values(values, spec, Ls, u, offsets)
        loglik = _laplace_subject_value_gradient!(gradient,
            laplace.objective.subject_objectives[i], aws, shifted)
        isfinite(loglik) || return (value=loglik, gradient=fill(T(NaN), length(u)))
        total += loglik
        for l in eachindex(spec.levels)
            level = spec.levels[l]
            k = nrandomeffects(level)
            k == 0 && continue
            base = offsets[l]
            L = Ls[l]
            for q in 1:k
                acc = zero(T)
                for pp in 1:k
                    acc += L[pp, q] * gradient[level.re_index[pp]]
                end
                inner[base + q] += acc
            end
        end
    end
    return (value=total, gradient=inner)
end

"""
    _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)

`(g_U(u), dg_U/du)` for unit `U`: the summed process log likelihood of its
subjects at their shifted parameter vectors, less `u'u/2`, and its gradient
with respect to `u`.

The chain rule through the shifts is one matrix-vector product per subject per
level: the likelihood's gradient with respect to that subject's shifted raw
vector, restricted to the level's varying positions, pushed back through that
level's Cholesky factor and accumulated into the level's block of `u`. Two
subjects in the same study accumulate into the *same* study block, which is
exactly the coupling that makes the study a single integration unit.
"""
function _laplace_unit_objective_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws) where {T}
    result = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws,
        eachindex(laplace.units.members[U]))
    isfinite(result.value) || return result
    inner = result.gradient
    @inbounds for a in eachindex(u)
        inner[a] -= u[a]
    end
    return (value=result.value - dot(u, u) / 2, gradient=inner)
end

"""Unit size at or above which the blocked curvature beats the dense one."""
const _LAPLACE_BLOCK_THRESHOLD = Ref(14)

export ctsem_set_block_threshold!
"""
    ctsem_set_block_threshold!(n)

Set the unit size at which curvature assembly switches from differentiating the
whole gradient to differentiating block by block. Exposed because the crossover
is a property of the machine and the model, not of the mathematics.
"""
function ctsem_set_block_threshold!(n::Integer)
    n >= 0 || throw(ArgumentError("threshold must be non-negative"))
    _LAPLACE_BLOCK_THRESHOLD[] = Int(n)
    return Int(n)
end

"""
    _laplace_unit_hessian(laplace, U, values, Ls, u, slot; dense=nothing)

`d2 g_U / du du` -- the unit's inner curvature.

Formed one *block* of `u` at a time rather than one column at a time, and that
distinction is the difference between quadratic and linear cost in study size.

Column block `b` of the curvature is `d(dg/du)/du_b`, and a member that does
not sit under `b` has no dependence on `u_b` at all -- so its contribution is
structurally zero and it need not be filtered at all. Seeding a subject block
therefore costs one subject sweep, not one per member of the study. Summed over
blocks the cost is `n * sum_of_k_over_levels`, against `(k_study + n*k_subject)
* n` for differentiating the whole gradient at once. On a 40-subject study with
two effects at each level that is 160 sweeps rather than 3280.

The subtracted identity is the `-u'u/2` term, applied once here rather than
inside each block.

Which of the two runs is chosen per unit, because the blocked route is not
universally better: it makes one `ForwardDiff.jacobian` call per block instead
of one for the whole matrix, and that per-call overhead costs more than the
saved sweeps until a unit is reasonably large. Measured on a one-latent model,
five waves, one random effect at each level, milliseconds per curvature:

    subjects/study   dim(u)   blocked   dense   ratio
                 4        5      0.34    0.20    0.6x
                 8        9      0.73    0.44    0.6x
                16       17      1.43    1.86    1.3x
                32       33      3.07    6.14    2.0x

Dense grows quadratically in study size and blocked linearly, exactly as the
structure says they should, but they cross over around `dim(u)` of 12 to 16.
`_LAPLACE_BLOCK_THRESHOLD` is that crossover, and it is an empirical constant
rather than a derived one -- set from the table above, on one machine. Passing
`dense` explicitly overrides the choice, which is what the test that compares
the two routes does.
"""
function _laplace_unit_hessian(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    slot::Integer=1; dense::Union{Nothing,Bool}=nothing) where {T}
    d = length(u)
    d == 0 && return zeros(T, 0, 0)
    blocks = laplace.units.blocks[U]
    # One block is the whole matrix, so the two routes are the same computation
    # and the blocked one only adds a layer.
    usedense = dense === nothing ?
        (length(blocks) <= 1 || d < _LAPLACE_BLOCK_THRESHOLD[]) : dense
    if usedense
        inner_of = function (uu)
            S = eltype(uu)
            ws = _laplace_workspace!(laplace, S, length(values), slot)
            vs = convert(Vector{S}, values)
            Lss = [convert(Matrix{S}, L) for L in Ls]
            return _laplace_unit_objective_gradient(laplace, U, vs, Lss, uu, ws).gradient
        end
        H = ForwardDiff.jacobian(inner_of, collect(u))
        return (H .+ transpose(H)) ./ 2
    end

    base = collect(u)
    H = zeros(T, d, d)
    for (offset, size, positions) in blocks
        columns = (offset + 1):(offset + size)
        block_of = function (ub)
            S = eltype(ub)
            ws = _laplace_workspace!(laplace, S, length(values), slot)
            vs = convert(Vector{S}, values)
            Lss = [convert(Matrix{S}, L) for L in Ls]
            uu = convert(Vector{S}, base)
            @inbounds for (t, c) in enumerate(columns)
                uu[c] = ub[t]
            end
            return _laplace_unit_loglik_gradient(laplace, U, vs, Lss, uu, ws,
                positions).gradient
        end
        H[:, columns] = ForwardDiff.jacobian(block_of, base[columns])
    end
    @inbounds for a in 1:d
        H[a, a] -= one(T)
    end
    return (H .+ transpose(H)) ./ 2
end

"""
    _laplace_solve_unit_mode!(laplace, U, values, Ls, slot)

Newton's method on `g_U`, warm-started from the retained mode for unit `U`.

`g_U` is a log likelihood minus a quadratic, so its curvature is negative
definite near the mode and Newton is the right method; away from the mode, and
for a nonlinear process model, it need not be. Two guards, both visible in the
diagnostics rather than silent: a curvature that is not negative definite is
shifted until it is and the unit flagged as repaired, and a step that does not
improve `g_U` is halved before the iteration gives up.
"""
function _laplace_solve_unit_mode!(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, slot::Integer=1)
    d = laplace.units.dims[U]
    u = copy(laplace.modes[U])
    aws = _laplace_workspace!(laplace, Float64, length(values), slot)
    repaired = false
    converged = false
    iterations = 0
    current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)
    if !isfinite(current.value)
        # A warm start can be stranded outside the support after a large outer
        # step. The origin is always inside it: u = 0 is the population mean.
        fill!(u, 0.0)
        current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)
    end
    for iteration in 1:laplace.inner_maxiter
        iterations = iteration
        if d == 0 || maximum(abs, current.gradient) < laplace.inner_tol
            converged = true
            break
        end
        H = _laplace_unit_hessian(laplace, U, values, Ls, u, slot)
        negative_definite, Hneg = _laplace_negate_definite(H)
        repaired |= !negative_definite
        step = Hneg \ current.gradient
        accepted = false
        scale = 1.0
        for _ in 1:20
            candidate = u .+ scale .* step
            trial = _laplace_unit_objective_gradient(laplace, U, values, Ls, candidate, aws)
            if isfinite(trial.value) && trial.value >= current.value - 1e-12
                u = candidate
                current = trial
                accepted = true
                break
            end
            scale /= 2
        end
        accepted || break
    end
    if d == 0 || maximum(abs, current.gradient) < laplace.inner_tol
        converged = true
    end
    laplace.modes[U] = u
    laplace.inner_iterations[U] = iterations
    laplace.inner_gradient[U] = d == 0 ? 0.0 : maximum(abs, current.gradient)
    laplace.inner_converged[U] = converged
    laplace.hessian_repaired[U] = repaired
    return (u=u, value=current.value, converged=converged)
end

"""
    _laplace_dual_unit_mode(laplace, U, values, Ls, uhat, Hneg, aws)

The unit's inner mode as a function of the outer parameters, to first order.

One Newton step from the converged primal mode, taken with dual parameters.
The inner gradient's *primal* part is zero there, so the step's primal part is
zero and its dual part is exactly `-H^-1 dg/dtheta`: the implicit function
theorem, without forming that cross-derivative explicitly. Using the primal
`Hneg` for the solve rather than a dual one is not an approximation for the
same reason -- any dual part of the inverse would multiply a zero primal
gradient.
"""
function _laplace_dual_unit_mode(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix},
    uhat::Vector{Float64}, Hneg::Matrix{Float64}, aws) where {T}
    isempty(uhat) && return Vector{T}(undef, 0)
    u0 = convert(Vector{T}, uhat)
    inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u0, aws)
    return u0 .+ (Hneg \ inner.gradient)
end

"""
    _laplace_unit_term(laplace, U, values, Ls, u, aws, slot)

`g_U(u) - logdet(-d2 g_U/du du) / 2`: unit `U`'s contribution to the
approximated log marginal likelihood.

The `(2*pi)^(d/2)` the Laplace approximation produces cancels the
`(2*pi)^(-d/2)` in the standard-normal density of `u`, so no dimension-dependent
constant appears here whatever the unit's size.
"""
function _laplace_unit_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws, slot::Integer=1) where {T}
    inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)
    isfinite(inner.value) || return inner.value
    isempty(u) && return inner.value
    H = _laplace_unit_hessian(laplace, U, values, Ls, u, slot)
    negated = -(H .+ transpose(H)) ./ 2
    factorization = cholesky(Symmetric(negated); check=false)
    issuccess(factorization) || return T(NaN)
    return inner.value - logdet(factorization) / 2
end

################################################################################
# The exact outer gradient, in k + 1 reverse sweeps per subject
################################################################################
#
# The straightforward way to differentiate the Laplace term is to run
# `ForwardDiff` over the whole of it, which is what `_laplace_nested_gradient`
# below still does. That costs `O(npar * k)` reverse sweeps per subject --
# `npar` outer forward directions, each carrying `k` more for the inner
# curvature -- and it is the dominant cost of a Laplace fit.
#
# It is avoidable, and the identity that avoids it is worth stating plainly.
# Write `M = -H` and `C = inv(M)`, and note that for *fixed* vectors `q_r` with
# `sum_r q_r q_r' = C`,
#
#     tr(C dH/dx) = sum_r  q_r' (dH/dx) q_r  =  sum_r d/dx ( q_r' H q_r )
#
# Holding `q_r` fixed is legitimate: the identity needs `C`'s value at the
# current point, not its derivative. And `q' H q` is a *directional* second
# derivative -- so seeding one direction with a second-order dual and running
# the engine's ordinary reverse pass returns its gradient with respect to every
# parameter at once, because reverse differentiation in the parameters and
# forward differentiation in the seed commute.
#
# Concretely, with `H = L' A L - I` for `A` the log likelihood's Hessian block
# on the varying positions, `psi := tr(C H) = tr(W A) - tr(C)` where
# `W = L C L'`. Factoring `W = Q Q'` puts the seeded directions in *parameter*
# space rather than in `z` space, which matters: a direction in `z` space moves
# with `L(theta)` while a parameter-space direction does not, and that is what
# keeps the bookkeeping below finite.
#
# So: `k` sweeps seeded along the columns of `Q`, one more along `L s` for the
# term carrying the mode's own dependence, and the whole gradient falls out.

struct _LaplaceSeedInner end
struct _LaplaceSeedOuter end

"""
    _laplace_directional_pass(laplace, i, base, direction, order)

One reverse sweep of subject `i` at `base`, with the parameter vector seeded
along `direction` by a dual of the given order.

Returns the seed-order coefficients of the gradient the sweep produces:

  * `d0` -- the ordinary gradient at `base`;
  * `d1` -- its first derivative along `direction`, i.e. the log likelihood's
    Hessian contracted with `direction`;
  * `d2` -- its second derivative along `direction` (order 2 only), i.e. the
    gradient of the directional second derivative.

`d2` is the point of the whole exercise. Obtaining it from a single sweep, for
every parameter simultaneously, is what replaces `npar` forward directions.

Nesting two one-dimensional duals rather than using a single second-order
partial is deliberate: with `x = a + d1 + d2`, the cross term `d1*d2` of `f(x)`
is exactly `f''(a)`, with no factorial to remember and no chunking to
configure.
"""
function _laplace_directional_pass(laplace::CTSEMLaplaceObjective, i::Integer,
    base::Vector{Float64}, direction::Vector{Float64}, order::Integer,
    slot::Integer=1)
    npar = length(base)
    subject = laplace.objective.subject_objectives[i]
    failure = (ok=false, d0=Float64[], d1=Float64[], d2=Float64[])
    if order == 2
        seed = ForwardDiff.Dual{_LaplaceSeedOuter}(
            ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0),
            ForwardDiff.Dual{_LaplaceSeedInner}(1.0, 0.0))
        S = typeof(seed)
        x = Vector{S}(undef, npar)
        @inbounds for m in 1:npar
            x[m] = base[m] + seed * direction[m]
        end
        aws = _laplace_workspace!(laplace, S, npar, slot)
        g = Vector{S}(undef, npar)
        loglik = _laplace_subject_value_gradient!(g, subject, aws, x)
        isfinite(ForwardDiff.value(ForwardDiff.value(loglik))) || return failure
        d0 = Vector{Float64}(undef, npar)
        d1 = Vector{Float64}(undef, npar)
        d2 = Vector{Float64}(undef, npar)
        @inbounds for m in 1:npar
            inner = ForwardDiff.value(g[m])
            d0[m] = ForwardDiff.value(inner)
            d1[m] = ForwardDiff.partials(inner)[1]
            d2[m] = ForwardDiff.partials(ForwardDiff.partials(g[m])[1])[1]
        end
        return (ok=true, d0=d0, d1=d1, d2=d2)
    end
    seed = ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0)
    S = typeof(seed)
    x = Vector{S}(undef, npar)
    @inbounds for m in 1:npar
        x[m] = base[m] + seed * direction[m]
    end
    aws = _laplace_workspace!(laplace, S, npar, slot)
    g = Vector{S}(undef, npar)
    loglik = _laplace_subject_value_gradient!(g, subject, aws, x)
    isfinite(ForwardDiff.value(loglik)) || return failure
    d0 = Vector{Float64}(undef, npar)
    d1 = Vector{Float64}(undef, npar)
    @inbounds for m in 1:npar
        d0[m] = ForwardDiff.value(g[m])
        d1[m] = ForwardDiff.partials(g[m])[1]
    end
    return (ok=true, d0=d0, d1=d1, d2=Float64[])
end

"""
    _laplace_popchol_derivatives(values, spec)

`(positions, dL)`: where the population covariance parameters sit in the raw
vector, and `dL[t]` the derivative of `L` with respect to `values[positions[t]]`.

Cheap whatever the model, because `L` is a `k x k` Cholesky of something built
only from those parameters -- no filter, no data, no process model. Every
subject uses the same derivatives, so this is computed once per evaluation
rather than once per subject.
"""
function _laplace_popchol_derivatives(values::AbstractVector{Float64},
    spec::CTSEMLaplaceSpec)
    positions = vcat(spec.levels[1].sd_index, spec.levels[1].cor_index)
    k = nrandomeffects(spec.levels[1])
    isempty(positions) && return (positions, Matrix{Float64}[])
    chol_of = function (p)
        v = convert(Vector{eltype(p)}, values)
        @inbounds for (slot, position) in enumerate(positions)
            v[position] = p[slot]
        end
        return vec(_laplace_popchol(v, spec))
    end
    jacobian = ForwardDiff.jacobian(chol_of, values[positions])
    dL = [Matrix{Float64}(reshape(view(jacobian, :, t), k, k))
          for t in eachindex(positions)]
    return (positions, dL)
end

"""
    _laplace_seeded_subject_gradient!(out, laplace, i, values, L, positions, dL, z, Mneg)

Accumulate subject `i`'s exact contribution to `dT/dtheta` into `out`.

With `v(theta, z) = theta + S(L(theta) z)`, `g = ll(v) - z'z/2`, and the mode
`zhat(theta)` defined by `dg/dz = 0`,

    dT/dtheta = dg/dtheta + [ dpsi/dtheta + s' B ] / 2

where `psi = tr(C H)`, `s = C dpsi/dz`, and `B = d2g/dz dtheta` carries the
mode's own dependence through `dzhat/dtheta = C B`. The envelope theorem is
what removes `zhat` from the first term and leaves it only inside `psi`.

Every likelihood-dependent piece comes from the `k + 1` sweeps: the envelope
gradient from the seed-order-zero coefficient, `dpsi/dtheta` and `dpsi/dz` from
the second-order ones, and the likelihood half of `s' B` from the extra sweep
along `L s`. The only terms not from a sweep are `L`'s own derivatives, which
are the `k x k` matrices in `dL`.

Returns `false` if a factorization or a sweep failed, leaving `out` untouched
in the caller's judgement -- the caller falls back rather than proceeding with
a partial answer.
"""
function _laplace_seeded_subject_gradient!(out::Vector{Float64},
    laplace::CTSEMLaplaceObjective, i::Integer, values::Vector{Float64},
    L::Matrix{Float64}, positions::Vector{Int}, dL::Vector{Matrix{Float64}},
    z::Vector{Float64}, Mneg::Matrix{Float64}, slot::Integer=1)

    spec = laplace.spec
    rho = spec.levels[1].re_index
    k = length(z)
    npar = length(values)
    base = _laplace_subject_values(values, spec, L, z)

    if k == 0
        pass = _laplace_directional_pass(laplace, i, base, zeros(Float64, npar), 1, slot)
        pass.ok || return false
        out .+= pass.d0
        return true
    end

    Mfact = cholesky(Symmetric(Mneg); check=false)
    issuccess(Mfact) || return false
    C = inv(Mfact); C = (C .+ transpose(C)) ./ 2
    W = L * C * transpose(L)
    Wfact = cholesky(Symmetric((W .+ transpose(W)) ./ 2); check=false)
    issuccess(Wfact) || return false
    Q = Matrix(Wfact.L)

    # k second-order sweeps, one per column of Q. Their d2 parts sum to
    # grad_v tr(W A), which is grad_v psi up to the constant tr(C).
    llv = Vector{Float64}(undef, npar)
    P = zeros(Float64, npar)
    direction = zeros(Float64, npar)
    for r in 1:k
        fill!(direction, 0.0)
        @inbounds for p in 1:k
            direction[rho[p]] = Q[p, r]
        end
        pass = _laplace_directional_pass(laplace, i, base, direction, 2, slot)
        pass.ok || return false
        r == 1 && copyto!(llv, pass.d0)
        P .+= pass.d2
    end

    Prho = Vector{Float64}(undef, k)
    llrho = Vector{Float64}(undef, k)
    @inbounds for p in 1:k
        Prho[p] = P[rho[p]]
        llrho[p] = llv[rho[p]]
    end
    s = C * (transpose(L) * Prho)

    # One more sweep, along L s: its d1 part is the likelihood factor of s' B.
    u = L * s
    fill!(direction, 0.0)
    @inbounds for p in 1:k
        direction[rho[p]] = u[p]
    end
    extra = _laplace_directional_pass(laplace, i, base, direction, 1, slot)
    extra.ok || return false
    Ds = extra.d1

    # dv/dtheta is the identity away from the population parameters, so for
    # every other parameter the contribution is a plain read-off.
    @inbounds for j in 1:npar
        out[j] += llv[j] + (P[j] + Ds[j]) / 2
    end

    # The population parameters move v through L as well, and move psi through
    # L explicitly. `H + I = L' A L`, so `tr(C (H+I) inv(L) dL)` gives the
    # explicit term without ever forming A.
    if !isempty(positions)
        K = C * (Matrix{Float64}(LinearAlgebra.I, k, k) .- Mneg)
        @inbounds for t in eachindex(positions)
            j = positions[t]
            shift = dL[t] * z
            chain = 0.0
            for p in 1:k
                chain += (llrho[p] + (Prho[p] + Ds[rho[p]]) / 2) * shift[p]
            end
            sshift = dL[t] * s
            bpart = 0.0
            for p in 1:k
                bpart += llrho[p] * sshift[p]
            end
            solved = L \ dL[t]
            trace = 0.0
            for a in 1:k, b in 1:k
                trace += K[a, b] * solved[b, a]
            end
            out[j] += chain + (2 * trace + bpart) / 2
        end
    end
    return true
end

"""
    ctsem_laplace_evaluate(laplace, values; gradient=true)

The approximated log marginal likelihood, and optionally its gradient.

The inner modes are always found in ordinary `Float64` arithmetic first: they
are the solution of an optimization problem, and an optimizer's iterates carry
no useful derivative information. Only the converged mode does, and it gets it
from `_laplace_dual_mode`.

The gradient differentiates the complete per-subject term, log determinant and
implicit mode dependence included.
"""
function ctsem_laplace_evaluate(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, nested_gradient::Bool=false)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    nsubjects = length(laplace.objective.subject_objectives)
    nunits = length(laplace.units.members)
    single_level = nlevels(laplace.spec) == 1

    # 1. Inner modes and the value at them, in primal arithmetic, warm-started
    #    from the last call. Each subject's term is its own approximated log
    #    marginal likelihood, which is the per-subject quantity that means the
    #    same thing here as `subject_loglik` does without random effects.
    #
    #    The curvature is computed once and used three times -- for the
    #    definiteness check, for the log determinant, and as the primal solve in
    #    `_laplace_dual_mode` below. It is the most expensive primal quantity
    #    here, so recomputing it for each of those would be a third of the
    #    primal pass thrown away.
    Ls = _laplace_popchols(theta, laplace.spec)
    L = Ls[1]
    primal_hessians = Vector{Matrix{Float64}}(undef, nunits)
    unit_loglik = zeros(Float64, nunits)
    subject_loglik = zeros(Float64, nsubjects)
    value = 0.0

    # The subject loop is the parallelism here, and it is the natural one: each
    # subject's mode solve, curvature and value are completely independent, and
    # nothing is shared but the read-only parameter vector. Chunk count comes
    # from `ctsem_set_max_chunks!`, which is what the R side sets from `cores`,
    # so it is the same control the non-Laplace path uses.
    nchunks = _ctsem_nchunks(nunits)
    # Grow the per-chunk workspace stores serially, before anything is spawned.
    while length(laplace.workspaces) < nchunks
        push!(laplace.workspaces, Dict{Any,Any}())
    end
    ranges = _ctsem_chunk_ranges(nunits, nchunks)
    chunk_ok = fill(true, nchunks)
    chunk_bad = fill(NaN, nchunks)
    run_primal = function (c)
        aws = _laplace_workspace!(laplace, Float64, length(theta), c)
        @inbounds for U in ranges[c]
            _laplace_solve_unit_mode!(laplace, U, theta, Ls, c)
            u = laplace.modes[U]
            H = isempty(u) ? zeros(Float64, 0, 0) :
                _laplace_unit_hessian(laplace, U, theta, Ls, u, c)
            _, negated = _laplace_negate_definite(H)
            primal_hessians[U] = negated
            inner = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
            term = if !isfinite(inner.value) || isempty(u)
                inner.value
            else
                factorization = cholesky(Symmetric(negated); check=false)
                issuccess(factorization) ? inner.value - logdet(factorization) / 2 : NaN
            end
            if !isfinite(term)
                chunk_ok[c] = false
                chunk_bad[c] = term
                return nothing
            end
            unit_loglik[U] = term
            # A unit's term is the marginal likelihood of all its subjects
            # jointly, and with a study level it does not decompose over them:
            # the study effect is shared. Attributing it to the first member
            # would be a fiction, so it is spread evenly and the honest
            # per-unit figures are what the diagnostics expose.
            share = term / length(laplace.units.members[U])
            for i in laplace.units.members[U]
                subject_loglik[i] = share
            end
        end
        return nothing
    end
    if nchunks <= 1
        run_primal(1)
    else
        Threads.@sync for c in 1:nchunks
            Threads.@spawn run_primal(c)
        end
    end
    @inbounds for c in 1:nchunks
        chunk_ok[c] || return (value=chunk_bad[c],
            gradient=gradient ? fill(NaN, length(theta)) : nothing,
            subject_loglik=subject_loglik, converged=all(laplace.inner_converged))
    end
    value = sum(unit_loglik) + _ctsem_log_prior(laplace.objective, theta)

    gradient || return (value=value, gradient=nothing,
        subject_loglik=subject_loglik, converged=all(laplace.inner_converged))

    # 3. The gradient, by `k + 1` seeded reverse sweeps per subject.
    #
    # `_laplace_nested_gradient` computes the same thing by running ForwardDiff
    # over the whole per-subject term. It is kept because `test_laplace.jl`
    # checks the two against each other: they share the primal and nothing
    # else, so agreement to machine precision is a real check on the seeded
    # assembly, which has a lot of chain rule in it. Set
    # `nested_gradient = true` to use it.
    grad = zeros(Float64, length(theta))
    # The seeded scheme is specialised to one level: it factors the inner
    # curvature per subject and seeds directions in that subject's parameter
    # space. With a study level the units couple subjects and that
    # factorisation is not the right one, so the general route is used until
    # the seeded one is generalised.
    if !nested_gradient && single_level
        positions, dL = _laplace_popchol_derivatives(theta, laplace.spec)
        # One accumulator per chunk rather than one shared vector: the subject
        # contributions are a sum, and summing per chunk and then across chunks
        # is the same sum in a different order.
        partials = [zeros(Float64, length(theta)) for _ in 1:nchunks]
        fill!(chunk_ok, true)
        run_gradient = function (c)
            @inbounds for i in ranges[c]
                z = laplace.modes[i]
                if !_laplace_seeded_subject_gradient!(partials[c], laplace, i, theta,
                        L, positions, dL, z, primal_hessians[i], c)
                    chunk_ok[c] = false
                    return nothing
                end
            end
            return nothing
        end
        if nchunks <= 1
            run_gradient(1)
        else
            Threads.@sync for c in 1:nchunks
                Threads.@spawn run_gradient(c)
            end
        end
        ok = all(chunk_ok)
        if ok
            for c in 1:nchunks
                grad .+= partials[c]
            end
            _ctsem_log_prior_gradient!(grad, laplace.objective, theta)
        else
            # A failed factorization or a non-finite sweep is not a licence to
            # return a partial sum: fall back to the route that does not need
            # those factorizations.
            fill!(grad, 0.0)
            grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_hessians)
        end
    else
        grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_hessians)
    end
    return (value=value, gradient=grad, subject_loglik=subject_loglik,
        converged=all(laplace.inner_converged))
end

"""
    _laplace_nested_gradient(laplace, theta, L, hessians)

The exact outer gradient by ForwardDiff over the whole per-subject term.

`O(npar * k)` reverse sweeps per subject, which is what
`_laplace_seeded_subject_gradient!` exists to avoid. Retained as the oracle the
seeded path is tested against, and as its fallback when a factorization the
seeded path needs is not available.
"""
function _laplace_nested_gradient(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}},
    hessians::Vector{Matrix{Float64}})
    nunits = length(laplace.units.members)
    total_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        accumulated = zero(S)
        for U in 1:nunits
            uhat = laplace.modes[U]
            ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, hessians[U], wsd)
            accumulated += _laplace_unit_term(laplace, U, x, Lsd, ud, wsd)
        end
        return accumulated + _ctsem_log_prior(laplace.objective, x)
    end
    return ForwardDiff.gradient(total_of, theta)
end

"""
    ctsem_laplace_subject_values(laplace, values)

Each subject's own raw parameter vector at the current inner modes: the
population vector shifted by that subject's random effects, and then by its
TI-predictor effects.

Both shifts come from the code that already applies them during fitting --
`_laplace_subject_values` and the engine's own `_materialize_subject_values!`
-- rather than being reconstructed by the caller. That matters because the
caller is the R side's subject-parameter reporting, and a second, slightly
different copy of "what parameters does this subject have" is exactly the kind
of divergence that shows up as a summary disagreeing with the fit.

Returns `nsubjects x length(values)`; push a row through the model's transforms
to get that subject's parameter matrices.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; from_level::Integer=1)
    theta = collect(Float64, values)
    Ls = _laplace_popchols(theta, laplace.spec)
    nsubjects = length(laplace.objective.subject_objectives)
    out = zeros(Float64, nsubjects, length(theta))
    buffer = Float64[]
    for U in eachindex(laplace.units.members)
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = _laplace_restrict_levels(laplace, U, laplace.modes[U], from_level)
        for (m, i) in enumerate(laplace.units.members[U])
            shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                laplace.units.offsets[U][m])
            subject = laplace.objective.subject_objectives[i]
            _materialize_subject_values!(buffer, shifted, subject.params, subject.tipreds)
            out[i, :] = buffer
        end
    end
    return out
end

export ctsem_laplace_subject_values

"""
    _laplace_restrict_levels(laplace, U, u, from_level)

A copy of the unit's latent vector with every level inside `from_level` zeroed.

Naming a level means "this level and everything outside it": `from_level = 1`
keeps all of them, which is the ordinary subject-level answer; a value one past
the last level keeps none, which is the population answer. Zeroing rather than
dropping keeps the vector's layout intact, so nothing downstream has to know a
restriction happened.
"""
function _laplace_restrict_levels(laplace::CTSEMLaplaceObjective, U::Integer,
    u::Vector{Float64}, from_level::Integer)
    from_level <= 1 && return copy(u)
    out = copy(u)
    units = laplace.units
    for (m, _) in enumerate(units.members[U])
        for l in 1:min(from_level - 1, nlevels(laplace.spec))
            k = nrandomeffects(laplace.spec.levels[l])
            k == 0 && continue
            base = units.offsets[U][m][l]
            @inbounds for q in 1:k
                out[base + q] = 0.0
            end
        end
    end
    return out
end

"""
    ctsem_kalman(laplace, values; from_level=1, subject_matrices=true)

Filter a Laplace fit with each subject at its own realized parameters.

`from_level` selects which levels of random effect are included: 1 (the
default) gives each subject its own estimate, 2 gives study-level effects only
so every subject in a study shares its mean trajectory, and one past the last
level gives the population trajectory with no random effects at all.

The random effects used are the *modes* -- estimated from all of a subject's
data at once. That makes these the smoothed-equivalent trajectories, not
filtered ones: an augmented fit's carrier states are updated observation by
observation, so its filtered output shows a random effect being learned, and
this cannot. The R side says so when it is used.
"""
function ctsem_kalman(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    from_level::Integer=1, subject_matrices::Bool=true)
    persubject = ctsem_laplace_subject_values(laplace, values; from_level=from_level)
    return ctsem_kalman(laplace.objective, persubject;
        subject_matrices=subject_matrices)
end

"""
    ctsem_laplace_mode_jacobian(laplace, values)

`dzhat/dtheta` for every unit: how each unit's inner mode moves when the
population parameters move.

The mode is defined implicitly by `dg/dz = 0`, so differentiating that gives
`dzhat/dtheta = C * d2g/dz dtheta`, and `_laplace_dual_mode` already produces
exactly that -- one Newton step from the converged mode taken in dual
arithmetic. Seeding `theta` with a full set of forward directions turns the
single directional answer into the whole Jacobian.

Costs one forward sweep per parameter chunk over each unit, once, which is why
it is worth computing here and reusing across every posterior draw rather than
re-solving a mode per draw.

Returns a vector of `dim(u) x npar` matrices, one per unit.
"""
function ctsem_laplace_mode_jacobian(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    nunits = length(laplace.units.members)
    Ls = _laplace_popchols(theta, laplace.spec)
    hessians = Vector{Matrix{Float64}}(undef, nunits)
    for U in 1:nunits
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = laplace.modes[U]
        H = isempty(u) ? zeros(Float64, 0, 0) :
            _laplace_unit_hessian(laplace, U, theta, Ls, u)
        _, negated = _laplace_negate_definite(H)
        hessians[U] = negated
    end
    out = Vector{Matrix{Float64}}(undef, nunits)
    for U in 1:nunits
        d = laplace.units.dims[U]
        if d == 0
            out[U] = zeros(Float64, 0, length(theta))
            continue
        end
        uhat = laplace.modes[U]
        mode_of = function (x)
            S = eltype(x)
            wsd = _laplace_workspace!(laplace, S, length(x))
            Lsd = _laplace_popchols(x, laplace.spec)
            return _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, hessians[U], wsd)
        end
        out[U] = ForwardDiff.jacobian(mode_of, theta)
    end
    return out
end

export ctsem_laplace_mode_jacobian

"""
    ctsem_laplace_subject_values(laplace, draws, thetahat)

Each subject's own raw parameter vector at every posterior draw.

The population Cholesky factors are rebuilt exactly at each draw, because they
are `k x k` and cost nothing. Only the *mode* is approximated, by

    uhat(theta) ~= uhat + (duhat/dtheta) (theta - thetahat)

rather than re-solved. That is deliberate and it is an approximation: solving a
mode per draw is a Newton iteration per unit per draw, where this is a
matrix-vector product, and the linearisation is accurate exactly where a
normal-approximation posterior puts its draws. It degrades for draws far from
the estimate, which is also where the normal approximation the draws come from
is itself least trustworthy.

Returns `ndraws x nsubjects x npar`.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    draws::AbstractMatrix, thetahat::AbstractVector)
    # The Jacobian is computed here rather than handed in: it is a vector of
    # matrices, and round-tripping one through the R bridge only to send it
    # straight back costs two marshalling steps for no gain.
    jacobians = ctsem_laplace_mode_jacobian(laplace, thetahat)
    ndraws = size(draws, 1)
    npar = size(draws, 2)
    nsubjects = length(laplace.objective.subject_objectives)
    centre = collect(Float64, thetahat)
    out = zeros(Float64, ndraws, nsubjects, npar)
    buffer = Float64[]
    for s in 1:ndraws
        theta = collect(Float64, view(draws, s, :))
        Ls = _laplace_popchols(theta, laplace.spec)
        step = theta .- centre
        for U in eachindex(laplace.units.members)
            u = laplace.units.dims[U] == 0 ? Float64[] :
                laplace.modes[U] .+ jacobians[U] * step
            for (m, i) in enumerate(laplace.units.members[U])
                shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                    laplace.units.offsets[U][m])
                subject = laplace.objective.subject_objectives[i]
                _materialize_subject_values!(buffer, shifted, subject.params,
                    subject.tipreds)
                out[s, i, :] = buffer
            end
        end
    end
    return out
end

"""
    ctsem_laplace_population(laplace, values)

Raw-scale population standard deviations and correlations, for a whole matrix
of raw parameter vectors at once (`values[s, :]` is one draw).

Batched because the caller is the summary, which has a posterior sample rather
than a point: one call for a thousand draws instead of a thousand calls. The
correlations come back in the same column-major lower-triangular order as the
parameters that produced them, which is the order `sdcovsqrt2cov`'s and Stan's
own correlation coordinates use.
"""
function ctsem_laplace_population(laplace::CTSEMLaplaceObjective, values::AbstractMatrix,
    level::Integer=1)
    nsamples = size(values, 1)
    lv = laplace.spec.levels[level]
    k = nrandomeffects(lv)
    noffdiagonals = div(k * (k - 1), 2)
    sd = zeros(Float64, nsamples, k)
    correlation = zeros(Float64, nsamples, noffdiagonals)
    for s in 1:nsamples
        L = _laplace_popchol(collect(Float64, view(values, s, :)), lv)
        covariance = L * transpose(L)
        scales = sqrt.(max.(diag(covariance), 0.0))
        sd[s, :] = scales
        counter = 0
        for j in 1:k, i in (j + 1):k
            counter += 1
            denominator = scales[i] * scales[j]
            correlation[s, counter] = denominator > 0 ? covariance[i, j] / denominator : 0.0
        end
    end
    return (sd=sd, correlation=correlation)
end

export ctsem_laplace_population

"""
    ctsem_laplace_modes(laplace, values)

The per-subject random-effect modes at `values`, on both the standardized
`z` scale and the raw parameter scale, with their conditional standard errors.

The conditional covariance of `z_i` is the inverse of the inner precision at
the mode, which is the same matrix the log determinant is taken of, so this
costs nothing that the objective did not already compute. On the raw scale it
is `L * cov * L'`: the same change of variables the model itself applies.
"""
function ctsem_laplace_modes(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    level::Integer=1)
    theta = collect(Float64, values)
    spec = laplace.spec
    lv = spec.levels[level]
    k = nrandomeffects(lv)
    Ls = _laplace_popchols(theta, spec)
    L = Ls[level]
    ngroups = lv.ngroups
    z = zeros(Float64, ngroups, k)
    raw = zeros(Float64, ngroups, k)
    zsd = zeros(Float64, ngroups, k)
    rawsd = zeros(Float64, ngroups, k)
    filled = falses(ngroups)
    for U in eachindex(laplace.units.members)
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = laplace.modes[U]
        covariance = if isempty(u)
            zeros(Float64, 0, 0)
        else
            H = _laplace_unit_hessian(laplace, U, theta, Ls, u)
            _, negated = _laplace_negate_definite(H)
            inv(Symmetric(negated))
        end
        for (m, i) in enumerate(laplace.units.members[U])
            g = lv.group[i]
            # Several subjects share a group at an outer level; its mode is one
            # vector, not one per subject, so it is written once.
            filled[g] && continue
            filled[g] = true
            k == 0 && continue
            base = laplace.units.offsets[U][m][level]
            slice = (base + 1):(base + k)
            zi = u[slice]
            z[g, :] = zi
            raw[g, :] = L * zi
            block = covariance[slice, slice]
            zsd[g, :] = sqrt.(max.(diag(block), 0.0))
            rawcov = L * block * transpose(L)
            rawsd[g, :] = sqrt.(max.(diag(rawcov), 0.0))
        end
    end
    return (z=z, raw=raw, z_sd=zsd, raw_sd=rawsd, parameter=lv.re_index,
        group=lv.group, ngroups=ngroups,
        converged=copy(laplace.inner_converged),
        iterations=copy(laplace.inner_iterations))
end

"""
    ctsem_laplace_diagnostics(laplace)

Inner-solve status from the last evaluation, per subject.
"""
ctsem_laplace_diagnostics(laplace::CTSEMLaplaceObjective) = (
    iterations=copy(laplace.inner_iterations),
    max_gradient=copy(laplace.inner_gradient),
    converged=copy(laplace.inner_converged),
    hessian_repaired=copy(laplace.hessian_repaired),
)

export ctsem_laplace_diagnostics

"""
    ctsem_laplace_optimize(laplace, start; ...)

Maximize with L-BFGS, mirroring `ctsem_optimize`'s contract so the R side can
treat the two the same way.

The Laplace value and the Laplace gradient are a consistent pair, so the line
search behaves and convergence means something.

There used to be a second mode here that maximized the value without its
log-determinant, using the cheaper envelope gradient. It is gone. Dropping the
log-determinant leaves the PQL-shaped objective, which is degenerate in the
variance components -- nothing in it penalises the population scales growing,
so they run away. A nine-replication recovery study put its population sd at
19.6 against a truth of 1.0, and one replication failed outright. It was not a
cheaper route to the same answer, so there is no version of it worth keeping.
"""
function ctsem_laplace_optimize(laplace::CTSEMLaplaceObjective, start::AbstractVector;
    maxiter::Integer=1000, g_tol::Real=1e-8, f_tol::Real=0.0, x_tol::Real=0.0,
    verbose::Bool=false, nested_gradient::Bool=false)
    start_values = collect(Float64, start)
    invalid_objective = floatmax(Float64) / 1e8
    gradient_limit = sqrt(floatmax(Float64))
    fg! = function (F, G, x)
        result = try
            ctsem_laplace_evaluate(laplace, x; gradient=G !== nothing,
                nested_gradient=nested_gradient)
        catch
            nothing
        end
        objective = result === nothing ? NaN : result.value
        valid = result !== nothing && isfinite(objective)
        if valid && G !== nothing
            valid = all(isfinite, result.gradient) &&
                all(abs(value) < gradient_limit for value in result.gradient)
        end
        if !valid
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        G !== nothing && (G .= -result.gradient)
        return F === nothing ? nothing : -objective
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol, f_reltol=f_tol,
        x_abstol=x_tol, show_trace=verbose, store_trace=false)
    if verbose
        chunks = ctsem_max_chunks()
        println("Laplace: ", length(laplace.objective.subject_objectives),
            " subjects in ", length(laplace.units.members), " unit(s), ",
            nlevels(laplace.spec), " level(s), ", nrandomeffects(laplace.spec),
            " random effects, ", min(max(chunks.max_chunks == 0 ? chunks.nthreads :
                chunks.max_chunks, 1), length(laplace.units.members)),
            " chunk(s) over ", chunks.nthreads, " thread(s)")
    end
    # Optim's default Hager-Zhang line search asserts its own bracketing
    # invariant (`B > A`) and *throws* when an evaluation it is handed is
    # invalid -- which happens here whenever a trial point makes an inner mode
    # solve or a curvature factorization fail, since those return a sentinel
    # objective with a zero gradient and a zero directional derivative breaks
    # the bracket. Backtracking makes no such assumption: it simply shrinks the
    # step. Falling back to it turns a crashed fit into a slower one, which is
    # the right trade, and `linesearch` on the result says which was used
    # rather than leaving it to be guessed.
    linesearch = "hagerzhang"
    result = try
        Optim.optimize(Optim.only_fg!(fg!), start_values, Optim.LBFGS(), options)
    catch err
        err isa InterruptException && rethrow()
        linesearch = "backtracking"
        verbose && println("Laplace: Hager-Zhang line search failed (",
            sprint(showerror, err), "); retrying with backtracking")
        Optim.optimize(Optim.only_fg!(fg!), start_values,
            Optim.LBFGS(linesearch=Optim.LineSearches.BackTracking()), options)
    end
    if verbose
        println("Laplace: inner modes ",
            count(laplace.inner_converged), "/", length(laplace.inner_converged),
            " converged, max |dg/dz| ",
            isempty(laplace.inner_gradient) ? 0.0 : maximum(laplace.inner_gradient),
            ", curvature repaired for ", count(laplace.hessian_repaired), " subject(s)")
    end
    minimizer = collect(Optim.minimizer(result))
    final = ctsem_laplace_evaluate(laplace, minimizer; gradient=true)
    return (
        minimizer=minimizer,
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        iterations=Optim.iterations(result),
        linesearch=linesearch,
        converged=Optim.converged(result),
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
        inner_converged=all(laplace.inner_converged),
        inner_iterations=copy(laplace.inner_iterations),
        hessian_repaired=copy(laplace.hessian_repaired),
    )
end

"""
    ctsem_laplace_hessian(laplace, values; step)

The outer Hessian of the approximated log marginal likelihood, for population
parameter standard errors.

A central difference of the *exact* outer gradient, not of the value: the
gradient is already third-order accurate, so differencing it once costs `2n`
gradient evaluations and inherits their accuracy, where differencing the value
twice would cost `O(n^2)` and lose half the digits. Nesting `ForwardDiff` once
more would be a fourth derivative of the process model; that is left until
there is evidence the difference matters, and this route is the one whose
error is at least bounded and reportable.
"""
function ctsem_laplace_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    step::Real=1e-4)
    x = collect(Float64, values)
    n = length(x)
    H = zeros(Float64, n, n)
    for j in 1:n
        h = step * max(1.0, abs(x[j]))
        plus = copy(x); plus[j] += h
        minus = copy(x); minus[j] -= h
        gp = ctsem_laplace_evaluate(laplace, plus; gradient=true).gradient
        gm = ctsem_laplace_evaluate(laplace, minus; gradient=true).gradient
        H[:, j] = (gp .- gm) ./ (2h)
    end
    return (H .+ transpose(H)) ./ 2
end

export ctsem_laplace_hessian

################################################################################
# Behaving like the objective it wraps
################################################################################
#
# The R side reaches the engine through a handful of entry points that all take
# "the objective", and a Laplace fit has to answer them too. They split cleanly
# in three.
#
#  1. Questions about the *model* -- its parameter layout, its matrices, which
#     cells are state-dependent. Integrating the random effects out does not
#     change any of them, so these forward to the wrapped objective unchanged.
#
#  2. Questions about the *fit* -- value, gradient, curvature, per-subject
#     scores. These have Laplace counterparts that mean the same thing about
#     the same model, so the generic entry points route to them. That is what
#     makes `ctOptimUncertainty` and the summary work without a Laplace branch
#     on the R side: they ask for the log posterior's curvature, and they get
#     the log *marginal* posterior's curvature, which is the right answer to
#     their question.
#
#  3. Questions that are genuinely not answered yet -- filtering and data
#     generation, both of which need a per-subject parameter vector rather than
#     one shared vector. These throw. Forwarding them to the wrapped objective
#     would run and return the population-level answer, silently, where the
#     caller asked for a subject-level one; a refusal is the honest result
#     until the subject-conditional versions exist.

for f in (:ctsem_parameter_layout, :ctsem_state_dependent_cells, :ctsem_parameter_matrices)
    @eval $f(laplace::CTSEMLaplaceObjective, args...; kwargs...) =
        $f(laplace.objective, args...; kwargs...)
end

"""
Evaluate a Laplace objective through the generic entry point.

Returns the approximated log *marginal* likelihood and its gradient. The
`gradient_method` names the process gradient elsewhere in the engine and has no
meaning here -- the Laplace gradient is a forward sweep over the reverse pass
regardless -- so it is accepted and ignored rather than rejected.
"""
function ctsem_evaluate(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    result = ctsem_laplace_evaluate(laplace, values; gradient=gradient)
    contributions || return (value=result.value, gradient=result.gradient)
    # No `row_loglik`: the integral is over a whole subject's trajectory, so a
    # single row has no marginal contribution to report. Returning the subject
    # terms and omitting the row ones is more honest than inventing a
    # decomposition that the approximation does not have.
    return (value=result.value, gradient=result.gradient,
        subject_loglik=result.subject_loglik)
end

"""The outer Hessian, for the generic curvature entry point."""
ctsem_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector; chunk::Integer=0) =
    ctsem_laplace_hessian(laplace, values)

"""The gradient of the approximated log marginal likelihood."""
function ctsem_adjoint_gradient(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    return (value=result.value, gradient=result.gradient)
end

"""
    ctsem_subject_gradients(laplace, values)

Per-subject scores: row `i` is the gradient of subject `i`'s own approximated
log marginal likelihood.

This is the quantity the OPG, sandwich and score-bootstrap uncertainty methods
consume, and for a marginal likelihood it is the marginal score rather than the
joint one -- the two differ by exactly the random-effect terms that have been
integrated out, so using the joint score would understate the uncertainty it is
there to measure.

It costs one sweep, not one per subject: the per-subject terms are assembled
into a vector and differentiated together, so the same forward directions serve
every row.
"""
function ctsem_subject_gradients(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    theta = collect(Float64, values)
    nsubjects = length(laplace.objective.subject_objectives)
    nunits = length(laplace.units.members)

    Ls = _laplace_popchols(theta, laplace.spec)
    primal_hessians = Vector{Matrix{Float64}}(undef, nunits)
    for U in 1:nunits
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = laplace.modes[U]
        H = isempty(u) ? zeros(Float64, 0, 0) :
            _laplace_unit_hessian(laplace, U, theta, Ls, u)
        _, negated = _laplace_negate_definite(H)
        primal_hessians[U] = negated
    end

    terms_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        out = Vector{S}(undef, nunits)
        for U in 1:nunits
            uhat = laplace.modes[U]
            ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, primal_hessians[U], wsd)
            out[U] = _laplace_unit_term(laplace, U, x, Lsd, ud, wsd)
        end
        return out
    end
    # Rows are *units*, and with a study level a unit is a study rather than a
    # subject: the study effect is shared, so the marginal likelihood does not
    # decompose over its members and there is no per-subject score to report.
    scores = ForwardDiff.jacobian(terms_of, theta)

    # Each subject carries its share of the prior, so the rows still sum to the
    # full posterior gradient -- the same convention `ctsem_subject_gradients`
    # uses without random effects.
    if !isempty(laplace.objective.prior_index) && nunits > 0
        share = 1 / nunits
        for U in 1:nunits
            _ctsem_log_prior_gradient!(view(scores, U, :), laplace.objective, theta, share)
        end
    end
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    value = sum(_laplace_unit_term(laplace, U, theta, Ls, laplace.modes[U], aws)
                for U in 1:nunits; init=0.0)
    return (value=value + _ctsem_log_prior(laplace.objective, theta), scores=scores)
end

for (f, what) in ((:ctsem_generate, "Data generation"),)
    @eval function $f(laplace::CTSEMLaplaceObjective, args...; kwargs...)
        throw(ArgumentError(string($what, " is not implemented for the Laplace ",
            "random-effect route yet. It needs a per-subject parameter vector ",
            "rather than one shared vector, and returning the population-level ",
            "answer instead would be a different quantity than the one asked ",
            "for. Use intoverpop=TRUE for this.")))
    end
end
