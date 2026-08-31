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
using Printf
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
    CTSEMLaplaceBlock(offset, size, members, level, ancestors)

One block of a unit's latent vector: where it sits, how big it is, which
members of the unit own it, which level it belongs to, and which blocks sit
above it in the hierarchy.

`ancestors` is what makes the curvature's sparsity usable. Two blocks couple
only if some member depends on both, and with strict nesting that happens only
when one contains the other -- so a block's nonzero couplings are exactly its
ancestors, and the elimination that exploits this produces no fill-in at all.
"""
struct CTSEMLaplaceBlock
    offset::Int
    size::Int
    members::Vector{Int}
    level::Int
    ancestors::Vector{Int}
end

"""
    CTSEMLaplaceUnits

How the integral factorises.

# Nesting is required, and is not an oversight

The hierarchy must be strictly nested: every group at an inner level belongs to
exactly one group at each outer level. Crossed designs -- subjects seen by
several raters, pupils in both a class and a neighbourhood that cut across each
other -- cannot be expressed here.

That is a consequence of what makes this fast rather than a missing feature. A
unit's curvature couples two blocks only when some member depends on both, and
with strict nesting that happens exactly when one block contains the other. The
sparsity pattern is then a *tree*, eliminating innermost-first produces no
fill-in, and the cost is `sum_b k_b^3` instead of `dim(u)^3` -- linear rather
than cubic in the members of a group. Crossing two levels puts a cycle in that
graph: the elimination fills in, the tree recursion in `quadrature.jl` no longer
enumerates the integral correctly, and the per-unit factorisation stops being
the right decomposition at all.

Supporting crossed effects means a different factorisation, not a relaxed check
here. Until then, the R side refuses such a model rather than silently treating
one grouping as nested inside the other.

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
    blocks::Vector{Vector{CTSEMLaplaceBlock}}
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
            [CTSEMLaplaceBlock[] for _ in 1:nsubjects])
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
    blocks = Vector{Vector{CTSEMLaplaceBlock}}(undef, nunits)
    for U in 1:nunits
        owners = Dict{Int,Vector{Int}}()
        sizes = Dict{Int,Int}()
        blocklevel = Dict{Int,Int}()
        for (m, i) in enumerate(members[U])
            for l in 1:nlev
                k = nrandomeffects(levels[l])
                k == 0 && continue
                off = offsets[U][m][l]
                push!(get!(owners, off, Int[]), m)
                sizes[off] = k
                blocklevel[off] = l
            end
        end
        # Innermost level first. The elimination below depends on that order:
        # a block is eliminated only once every block beneath it has been, which
        # is what keeps the factorization free of fill-in.
        ordered = sort(collect(keys(owners)); by = off -> (blocklevel[off], off))
        position = Dict(off => t for (t, off) in enumerate(ordered))
        blocks[U] = CTSEMLaplaceBlock[]
        for off in ordered
            members_here = sort(unique(owners[off]))
            l = blocklevel[off]
            # Ancestors: the blocks this one's members share at every outer
            # level. Strict nesting makes them the same for every member, so
            # reading them off the first is enough.
            m = members_here[1]
            ancestors = Int[]
            for outer in (l + 1):nlev
                nrandomeffects(levels[outer]) == 0 && continue
                push!(ancestors, position[offsets[U][m][outer]])
            end
            push!(blocks[U], CTSEMLaplaceBlock(off, sizes[off], members_here, l,
                ancestors))
        end
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
    # Repair *at the reported mode*, as distinct from `hessian_repaired`, which
    # is true if any Newton iterate needed shifting on the way there. Those mean
    # different things and only one of them is a reason to doubt the answer:
    # Newton away from a mode routinely passes through a point where the
    # curvature is not negative definite, and on a 40-subject model at the
    # generating parameters that happened for twelve subjects while every one of
    # their final curvatures was fine.
    mode_repaired::Vector{Bool}
end

function CTSEMLaplaceObjective(objective::CTSEMObjective, spec::CTSEMLaplaceSpec;
    inner_maxiter::Integer=50, inner_tol::Real=1e-10)
    nsubjects = length(objective.subject_objectives)
    units = _laplace_build_units(spec, nsubjects)
    nunits = length(units.members)
    return CTSEMLaplaceObjective{typeof(objective)}(objective, spec, units,
        [zeros(Float64, units.dims[U]) for U in 1:nunits],
        Int(inner_maxiter), Float64(inner_tol), [Dict{Any,Any}()],
        zeros(Int, nunits), zeros(Float64, nunits), falses(nunits), falses(nunits),
        falses(nunits))
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
    return _laplace_member_values!(shifted, values, spec, Ls, u, offsets)
end

"""
    _laplace_member_values!(shifted, values, spec, Ls, u, offsets)

The same shift, into a buffer the caller owns.

The sampler evaluates this once per subject per gradient and does hundreds of
thousands of gradients, so the allocating form's vector-per-subject is worth
removing there. Everywhere else the allocating form reads better and the
allocation is lost in the adjoint sweep that follows it.
"""
function _laplace_member_values!(shifted::AbstractVector, values::AbstractVector,
    spec::CTSEMLaplaceSpec, Ls::Vector{<:AbstractMatrix}, u::AbstractVector,
    offsets::Vector{Int})
    S = eltype(shifted)
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
    _laplace_ekf_workspace!(laplace, T, slot)

The filter's own workspace, per chunk rather than per subject.

`ContinuousEKFObjective` caches one workspace on itself, which is right when the
only parallelism partitions subjects -- the Laplace unit loop does, so no two
tasks ever filter the same subject at once. It is *wrong* for a sampler running
several chains at once, where every chain filters every subject: they would
share one set of buffers and quietly corrupt each other's filters.

The workspace depends only on the shared `EKFParameters`, never on the subject,
so one per chunk covers all of them.
"""
function _laplace_ekf_workspace!(laplace::CTSEMLaplaceObjective, ::Type{T},
    slot::Integer=1) where {T}
    store = laplace.workspaces[slot]
    key = (:ekf_workspace, T)
    cached = get(store, key, nothing)
    cached === nothing || return cached
    built = _init_continuous_ekf_workspace(T, laplace.objective.params)
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
    subject_objective, aws, values::AbstractVector{T};
    ekf_workspace=nothing) where {T}
    # `ekf_workspace` overrides the one cached on the subject. See
    # `_laplace_ekf_workspace!`: the cached one is shared between tasks that
    # filter the same subject, which only a sampler does.
    ws = ekf_workspace === nothing ?
        _get_or_init_objective_workspace!(subject_objective, T) : ekf_workspace
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

################################################################################
# Factorizing a unit's curvature without forming it
################################################################################
#
# A unit's curvature couples two blocks only when some member depends on both,
# and with strict nesting that happens only when one block contains the other.
# So the nonzero couplings of a block are exactly its ancestors, the sparsity
# pattern is a tree, and eliminating blocks innermost-first produces *no
# fill-in*: every Schur update lands on ancestor-ancestor couplings that were
# already nonzero.
#
# That gives `sum_b k_b^3` work instead of `dim(u)^3`, which for a study of `n`
# subjects is linear rather than cubic in `n`. The larger gain is memory: a
# study of 2000 subjects with two effects each has `dim(u) ~ 4000`, and a dense
# curvature for it is 128 MB *per unit*, held for every unit across the whole
# fit. The block form stores `sum_b k_b^2` plus the couplings, which is
# kilobytes.
#
# `CTSEMBlockMatrix` is symmetric and only the lower side is kept: `diag[b]` is
# the block's own `k_b x k_b` piece and `coupling[b][t]` is the coupling to
# `blocks[b].ancestors[t]`, stored as `k_b x k_a`.

"""
Symmetric block curvature for one unit, stored by block rather than densely.

Generic in its element type because the outer gradient differentiates through
the factorization: the nested route evaluates the log determinant with dual
numbers, so every step below has to work for those as readily as for `Float64`.
"""
struct CTSEMBlockMatrix{T}
    diag::Vector{Matrix{T}}
    coupling::Vector{Vector{Matrix{T}}}
end

"""An all-zero block matrix shaped for one unit."""
function CTSEMBlockMatrix(::Type{T}, blocks::Vector{CTSEMLaplaceBlock}) where {T}
    diag = [zeros(T, b.size, b.size) for b in blocks]
    coupling = [[zeros(T, b.size, blocks[a].size) for a in b.ancestors]
                for b in blocks]
    return CTSEMBlockMatrix{T}(diag, coupling)
end

CTSEMBlockMatrix(blocks::Vector{CTSEMLaplaceBlock}) =
    CTSEMBlockMatrix(Float64, blocks)

"""Reassemble a block matrix densely. Testing and small-unit use only."""
function _laplace_block_dense(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}, dim::Integer) where {T}
    out = zeros(T, dim, dim)
    for (b, block) in enumerate(blocks)
        rows = (block.offset + 1):(block.offset + block.size)
        out[rows, rows] .= M.diag[b]
        for (t, a) in enumerate(block.ancestors)
            cols = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
            out[rows, cols] .= M.coupling[b][t]
            out[cols, rows] .= transpose(M.coupling[b][t])
        end
    end
    return out
end

"""Take the block form of a dense symmetric matrix. Testing use only."""
function _laplace_block_of(dense::AbstractMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    M = CTSEMBlockMatrix(T, blocks)
    for (b, block) in enumerate(blocks)
        rows = (block.offset + 1):(block.offset + block.size)
        M.diag[b] .= dense[rows, rows]
        for (t, a) in enumerate(block.ancestors)
            cols = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
            M.coupling[b][t] .= dense[rows, cols]
        end
    end
    return M
end

"""The Cholesky type `_laplace_block_factor` produces for element type `T`."""
const _LaplaceCholesky{T} = LinearAlgebra.Cholesky{T,Matrix{T}}

"""One unit's factorized curvature: the per-block Choleskys and the eliminated
couplings, exactly what `_laplace_block_solve` and `_laplace_selected_inverse`
consume. Named because it is held one per unit for a whole evaluation and a
`Vector{Any}` of them makes every solve a dynamic dispatch."""
const _LaplaceFactorization{T} =
    Tuple{Vector{_LaplaceCholesky{T}},Vector{Vector{Matrix{T}}}}

"""
    _laplace_block_factor(M, blocks)

Eliminate the blocks innermost-first, returning `(ok, logdet, factors)`.

`factors` holds, per block, the Cholesky of its diagonal *after* every
descendant has been eliminated into it, and the couplings that were used --
enough to solve with. `ok` is false if any diagonal failed to factorize, which
is how a curvature that is not positive definite reports itself here.

The Schur update `A -= B' D^-1 B` is applied to the ancestor couplings, and
because a block's ancestors form a chain, every entry it touches is one that
was already nonzero. That is the no-fill-in property, and it is why this stays
linear in the number of members.
"""
function _laplace_block_factor(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    nb = length(blocks)
    diag = [copy(d) for d in M.diag]
    coupling = [[copy(c) for c in row] for row in M.coupling]
    # Concretely typed, not `Vector{Any}`: `factors[b] \ x` and `logdet(factors[b])`
    # are called once per block per solve, and through an `Any` element every one
    # of them is a dynamic dispatch. A unit with one block per subject makes that
    # two dispatches per subject per Newton step.
    factors = Vector{_LaplaceCholesky{T}}(undef, nb)
    total = zero(T)
    for b in 1:nb
        f = cholesky(Symmetric(_laplace_symmetrise(diag[b])); check=false)
        issuccess(f) || return (false, T(NaN), factors, coupling)
        # Numerically singular counts as failure, so the caller shifts it.
        #
        # `issuccess` asks whether the matrix is positive definite, not whether
        # it can be inverted usefully, and a block that is barely definite
        # factorizes happily and then gives an enormous inverse. The outer
        # gradient goes through exactly that inverse -- the implicit function
        # theorem step in `_laplace_dual_unit_mode` -- so the *value* at such a
        # point is fine while the *gradient* overflows to Inf. Trial points
        # were then rejected for a non-finite gradient, 67 of them in one
        # 25-subject fit, and the line search ran out of room and stopped 0.6
        # log units short with a gradient of 11.
        #
        # The test is on the Cholesky diagonal, which costs one pass over it,
        # and the threshold corresponds to a condition number around 1e15 --
        # so this fires where the inverse has no significant digits left, and
        # not on a merely awkward model.
        let dg = LinearAlgebra.diag(f.factors)
            lo, hi = extrema(abs, dg)
            hi > 0 && lo <= sqrt(eps(real(float(one(T))))) * hi &&
                return (false, T(NaN), factors, coupling)
        end
        factors[b] = f
        total += logdet(f)
        ancestors = blocks[b].ancestors
        isempty(ancestors) && continue
        # W = D^-1 B for this block's couplings, then push B' W into the
        # ancestor-ancestor entries.
        W = [f \ coupling[b][t] for t in eachindex(ancestors)]
        for t in eachindex(ancestors)
            a = ancestors[t]
            update = transpose(coupling[b][t]) * W[t]
            diag[a] .-= _laplace_symmetrise(update)
            for s in eachindex(ancestors)
                s == t && continue
                c = ancestors[s]
                # Where does the (a, c) coupling live? One of them is an
                # ancestor of the other; the update belongs on the lower one's
                # row.
                slot = findfirst(==(c), blocks[a].ancestors)
                if slot !== nothing
                    coupling[a][slot] .-= transpose(coupling[b][t]) * W[s]
                end
            end
        end
    end
    return (true, total, factors, coupling)
end

@inline _laplace_symmetrise(A) = (A .+ transpose(A)) ./ 2

"""
    _laplace_selected_inverse(factors, coupling, blocks)

The entries of `C = inv(M)` that lie in `M`'s own sparsity pattern: every
block's diagonal, and every block's coupling to each of its ancestors.

This is the Takahashi recursion. Writing the factorization as `M = L D L'` with
unit lower `L`, the below-diagonal entries for block `b` are
`L[a,b] = B[b,a]' inv(D_b)` over `a` in `anc(b)`, and then, taking blocks from
the outermost inwards,

    C[anc(b), b] = -C[anc(b), anc(b)] L[anc(b), b]
    C[b, b]      = inv(D_b) - L[anc(b), b]' C[anc(b), b]

The recursion closes on the stored pattern rather than spilling outside it: the
ancestors of a block form a chain, so for two of them one is an ancestor of the
other and their coupling is already a stored entry. Nothing dense is ever
formed, which is the whole point -- `C` itself is dense for an arrow matrix,
and only these selected entries are wanted.

They are wanted because `dH/dtheta` is block sparse, so

    tr(C dH/dtheta) = sum_b tr(C[b,b] dH[b,b]/dtheta)
                    + 2 sum_b sum_{a in anc(b)} tr(C[a,b] dH[b,a]/dtheta)

and every term needs only an entry this returns.

Returns `(diag, coupling)` shaped exactly like the `CTSEMBlockMatrix` it
inverts: `diag[b]` is `C[b,b]` and `coupling[b][t]` is `C[b, anc(b)[t]]`.
"""
function _laplace_selected_inverse(factors, elim, blocks::Vector{CTSEMLaplaceBlock})
    nb = length(blocks)
    Cdiag = [zeros(Float64, b.size, b.size) for b in blocks]
    Ccoup = [[zeros(Float64, b.size, blocks[a].size) for a in b.ancestors]
             for b in blocks]

    # Where block `a`'s coupling to block `c` is stored, and whether it needs
    # transposing to read as `C[a, c]`.
    function entry(a::Int, c::Int)
        a == c && return (:diag, a, 0, false)
        slot = findfirst(==(c), blocks[a].ancestors)
        slot !== nothing && return (:coup, a, slot, false)
        slot = findfirst(==(a), blocks[c].ancestors)
        slot !== nothing && return (:coup, c, slot, true)
        return (:none, 0, 0, false)
    end
    function readC(a::Int, c::Int)
        kind, i, t, flip = entry(a, c)
        kind === :diag && return Cdiag[i]
        kind === :coup && return flip ? transpose(Ccoup[i][t]) : Ccoup[i][t]
        # Outside the pattern. Two blocks with no ancestor relation never
        # couple in `M`, and no term of the trace above asks for them.
        return zeros(Float64, blocks[a].size, blocks[c].size)
    end

    for b in nb:-1:1
        ancestors = blocks[b].ancestors
        Dinv = inv(factors[b])
        if isempty(ancestors)
            Cdiag[b] .= _laplace_symmetrise(Dinv)
            continue
        end
        # L[a, b] = B[b,a]' inv(D_b), as a k_a x k_b block.
        Lb = [transpose(elim[b][t]) * Dinv for t in eachindex(ancestors)]
        # C[anc, b] = -C[anc, anc] L[anc, b], accumulated over the ancestor set.
        for t in eachindex(ancestors)
            a = ancestors[t]
            acc = zeros(Float64, blocks[a].size, blocks[b].size)
            for sidx in eachindex(ancestors)
                c = ancestors[sidx]
                acc .+= readC(a, c) * Lb[sidx]
            end
            Ccoup[b][t] .= .-transpose(acc)
        end
        # C[b,b] = inv(D_b) - L[anc,b]' C[anc,b]
        acc = copy(Dinv)
        for t in eachindex(ancestors)
            acc .-= transpose(Lb[t]) * transpose(Ccoup[b][t])
        end
        Cdiag[b] .= _laplace_symmetrise(acc)
    end
    return (Cdiag, Ccoup)
end

"""
    _laplace_block_solve(factors, coupling, blocks, rhs)

Solve `M x = rhs` using the factorization above, by forward substitution over
the elimination order and back substitution against it.
"""
function _laplace_block_solve(factors, coupling, blocks::Vector{CTSEMLaplaceBlock},
    rhs::Vector{T}) where {T}
    nb = length(blocks)
    y = [rhs[(blocks[b].offset + 1):(blocks[b].offset + blocks[b].size)] for b in 1:nb]
    # Forward: subtract each eliminated block's contribution from its ancestors.
    for b in 1:nb
        w = factors[b] \ y[b]
        for (t, a) in enumerate(blocks[b].ancestors)
            y[a] .-= transpose(coupling[b][t]) * w
        end
    end
    # Back: solve outermost-first, substituting into the blocks beneath.
    x = [zeros(T, blocks[b].size) for b in 1:nb]
    for b in nb:-1:1
        acc = copy(y[b])
        for (t, a) in enumerate(blocks[b].ancestors)
            acc .-= coupling[b][t] * x[a]
        end
        x[b] = factors[b] \ acc
    end
    out = zeros(T, sum(b.size for b in blocks; init=0))
    for b in 1:nb
        out[(blocks[b].offset + 1):(blocks[b].offset + blocks[b].size)] .= x[b]
    end
    return out
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

Nothing in the fitting path calls this any more: the production route is
`_laplace_unit_curvature`, which produces the same curvature negated and in
block form without materializing a dense matrix. It is kept because it computes
that quantity by a route knowing nothing about the block structure, which is
what makes it a usable reference for the tests that check the block assembly
and the block factorization.
"""
function _laplace_unit_hessian(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    slot::Integer=1; dense::Union{Nothing,Bool}=nothing) where {T}
    d = length(u)
    d == 0 && return zeros(T, 0, 0)
    blocks = laplace.units.blocks[U]
    usedense = dense === nothing ? true : dense
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
    for block in blocks
        offset = block.offset; size = block.size; positions = block.members
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
    _laplace_unit_curvature(laplace, U, values, Ls, u, slot)

`-d2g/du du` for a unit, in block form, assembled without ever building the
dense matrix.

Column block `b` of the log-likelihood's curvature is `d(dll/du)/du_b`, and a
member that does not sit under `b` has no dependence on `u_b`, so one Jacobian
per block gives that block's diagonal *and* its couplings to its ancestors --
which, by the tree sparsity, is every nonzero entry in its row. The `-u'u/2`
term contributes the identity, added once here.

Returned negated (`M = -H`) because that is the positive definite matrix
everything downstream wants: the precision of the Gaussian being fitted, whose
determinant is the approximation's normalizing constant.

There is one representation downstream, always the block one, but two ways of
filling it. Per-block assembly is `O(n)` in unit size where a single Jacobian
over the whole gradient is `O(n^2)`, but it makes one `ForwardDiff.jacobian`
call per block instead of one, and below roughly a dozen members that overhead
costs more than the sweeps it saves -- measured at two to three times more.
Small units are therefore assembled with one dense Jacobian and split into
blocks afterwards, which is the same matrix by a cheaper route at that size.
`_LAPLACE_BLOCK_THRESHOLD` is the crossover, and it is empirical rather than
derived; `ctsem_set_block_threshold!` moves it.
"""
function _laplace_unit_curvature(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    slot::Integer=1) where {T}
    blocks = laplace.units.blocks[U]
    base = collect(u)
    if length(blocks) <= 1 || length(u) < _LAPLACE_BLOCK_THRESHOLD[]
        gradient_of = function (uu)
            S = eltype(uu)
            ws = _laplace_workspace!(laplace, S, length(values), slot)
            vs = convert(Vector{S}, values)
            Lss = [convert(Matrix{S}, L) for L in Ls]
            return _laplace_unit_loglik_gradient(laplace, U, vs, Lss, uu, ws,
                eachindex(laplace.units.members[U])).gradient
        end
        A = ForwardDiff.jacobian(gradient_of, base)
        dense = Matrix{T}(LinearAlgebra.I, length(u), length(u)) .-
            _laplace_symmetrise(A)
        return _laplace_block_of(dense, blocks)
    end
    M = CTSEMBlockMatrix(T, blocks)
    for (b, block) in enumerate(blocks)
        columns = (block.offset + 1):(block.offset + block.size)
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
                block.members).gradient
        end
        J = ForwardDiff.jacobian(block_of, base[columns])
        # M = I - d2(sum ll)/du du, block by block.
        M.diag[b] .= Matrix{T}(LinearAlgebra.I, block.size, block.size) .-
            J[columns, :]
        for (t, a) in enumerate(block.ancestors)
            rows = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
            M.coupling[b][t] .= .-transpose(J[rows, :])
        end
    end
    return M
end

"""
    _laplace_repair_blocks!(M, blocks)

Shift the block diagonals until the whole matrix factorizes, reporting whether
anything had to be done.

The block equivalent of `_laplace_negate_definite`: the curvature has to be
negative definite for the approximation to mean anything, and it need not be on
the way to a mode or at one for a badly identified model. Shifting is the
standard repair; reporting it is what stops it being a silent change of
objective.
"""
function _laplace_repair_blocks!(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    return _laplace_factor_repaired!(M, blocks).repaired
end

"""
    _laplace_factor_repaired!(M, blocks)

Repair `M` if it needs it and return the factorization that proved it did not,
as `(repaired, ok, logdet, factors, coupling)`.

The same work as `_laplace_repair_blocks!` followed by `_laplace_block_factor`,
minus one whole factorization. The repair has to factorize to find out whether
anything is wrong, and every caller then factorized again to get the factors --
so the ordinary case, where nothing needs repairing, was paying for two
eliminations per Newton step and two more per evaluation.
"""
function _laplace_factor_repaired!(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}) where {T}
    ok, logdetM, factors, coupling = _laplace_block_factor(M, blocks)
    ok && return (repaired=false, ok=true, logdet=logdetM, factors=factors,
        coupling=coupling)
    scale = maximum((maximum(abs, d) for d in M.diag); init=one(real(T)))
    scale = isfinite(scale) && scale > 0 ? scale : one(real(T))
    shift = sqrt(eps(real(float(one(T)))))
    for _ in 1:30
        for d in M.diag
            for i in axes(d, 1); d[i, i] += shift * scale; end
        end
        ok, logdetM, factors, coupling = _laplace_block_factor(M, blocks)
        ok && return (repaired=true, ok=true, logdet=logdetM, factors=factors,
            coupling=coupling)
        shift *= 10
    end
    return (repaired=true, ok=false, logdet=T(NaN), factors=factors,
        coupling=coupling)
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
    aws = _laplace_workspace!(laplace, Float64, length(values), slot)
    warm = _laplace_newton_unit_mode(laplace, U, values, Ls, aws,
        copy(laplace.modes[U]), slot)
    best = warm
    # The retained mode is a *warm start*, not part of the definition of the
    # objective. If Newton reached the tolerance from it, the mode it found is
    # the mode, and where the previous evaluation happened to leave the warm
    # start cannot matter. If it did not, the value about to be returned does
    # depend on that history -- and then the objective the outer optimizer sees
    # is not a function of theta at all.
    #
    # That is not hypothetical. On a 40-subject model, an outer line search that
    # visited a distant point left the modes stranded, and the *same* parameter
    # vector then evaluated to -170856 where a fresh object gave -967. The
    # optimizer stopped at its starting values and reported convergence.
    #
    # So a warm start that fails is retried from the origin, which is the same
    # for every caller and always inside the support: u = 0 is the population
    # mean. The better of the two is kept, so a cold retry can only help.
    if !warm.converged && any(!iszero, laplace.modes[U])
        cold = _laplace_newton_unit_mode(laplace, U, values, Ls, aws,
            zeros(Float64, d), slot)
        if cold.converged || (isfinite(cold.value) &&
                (!isfinite(best.value) || cold.value > best.value))
            best = cold
        end
    end
    laplace.modes[U] = best.u
    laplace.inner_iterations[U] = best.iterations
    laplace.inner_gradient[U] = d == 0 ? 0.0 : maximum(abs, best.gradient)
    laplace.inner_converged[U] = best.converged
    laplace.hessian_repaired[U] = best.repaired
    return (u=best.u, value=best.value, converged=best.converged)
end

"""
    _laplace_inner_tolerance(laplace, value)

How small the inner gradient has to be for the mode to count as found.

`inner_tol` alone is an absolute bound, and an absolute bound on a gradient is
the same mistake the outer `g_tol` makes: it asks for a number of digits that
depends on how large the unit's objective happens to be. On a 25-subject
ordinal model one unit stalled at `1.009e-10` against a tolerance of `1e-10` --
missing by one percent -- and because a unit that misses makes the *whole*
trial point invalid, the outer line search lost 50 of its 89 evaluations to it
and stopped after three iterations, 129 log units short.

What actually has to be true is that the mode is precise enough for the value
and the gradient to be reproducible. At a maximum the value error is second
order in the mode error, so stopping at `|g| = ε` costs about `ε²/2λ` in the
value -- 2e-17 at `ε = 1e-8`, which is below the last bit of a log likelihood
of order 1e3. The outer gradient's error is first order, about `ε` times a
cross-derivative of order ten, so 1e-7 at the same `ε`; the outer tolerance
this feeds is `1e-6 * |value|`, around 1e-3 here, so that is four orders of
margin. A relative floor of `1e-10` with the absolute one kept underneath it
sits comfortably inside all of that.
"""
@inline _laplace_inner_tolerance(laplace::CTSEMLaplaceObjective, value::Real) =
    max(laplace.inner_tol, 1e-10 * (one(value) + abs(value)))

"""
    _laplace_newton_unit_mode(laplace, U, values, Ls, aws, start, slot)

Newton on `g_U` from one given starting point, reporting what happened rather
than writing anything back.

Split out of `_laplace_solve_unit_mode!` so that the same iteration can be run
twice from different starts -- see the cold retry there.
"""
function _laplace_newton_unit_mode(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    start::Vector{Float64}, slot::Integer)
    d = laplace.units.dims[U]
    u = start
    repaired = false
    converged = false
    iterations = 0
    current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)
    if !isfinite(current.value) && any(!iszero, u)
        fill!(u, 0.0)
        current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws)
    end
    for iteration in 1:laplace.inner_maxiter
        iterations = iteration
        if d == 0 || maximum(abs, current.gradient) <
                _laplace_inner_tolerance(laplace, current.value)
            converged = true
            break
        end
        M = _laplace_unit_curvature(laplace, U, values, Ls, u, slot)
        fac = _laplace_factor_repaired!(M, laplace.units.blocks[U])
        repaired |= fac.repaired
        fac.ok || break
        step = _laplace_block_solve(fac.factors, fac.coupling,
            laplace.units.blocks[U], current.gradient)
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
    if d == 0 || maximum(abs, current.gradient) <
            _laplace_inner_tolerance(laplace, current.value)
        converged = true
    end
    return (u=u, value=current.value, gradient=current.gradient,
        converged=converged, iterations=iterations, repaired=repaired)
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
    uhat::Vector{Float64}, curvature, aws) where {T}
    isempty(uhat) && return Vector{T}(undef, 0)
    u0 = convert(Vector{T}, uhat)
    inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u0, aws)
    factors, coupling = curvature
    return u0 .+ _laplace_block_solve(factors, coupling, laplace.units.blocks[U],
        inner.gradient)
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
    M = _laplace_unit_curvature(laplace, U, values, Ls, u, slot)
    ok, logdetM, _, _ = _laplace_block_factor(M, laplace.units.blocks[U])
    ok || return T(NaN)
    return inner.value - logdetM / 2
end

################################################################################
# The exact outer gradient, in O(members) reverse sweeps per unit
################################################################################
#
# The straightforward way to differentiate the Laplace term is to run
# `ForwardDiff` over the whole of it, which is what `_laplace_nested_gradient`
# below still does. That costs `O(npar * k)` reverse sweeps per unit -- `npar`
# outer forward directions, each carrying `k` more for the inner curvature --
# and it is the dominant cost of a Laplace fit.
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
# With one level that is the whole story. `C` is a `k x k` block per subject,
# `psi := tr(C H) = tr(W A) - tr(C)` for `W = L C L'` and `A` the log
# likelihood's Hessian on the varying positions, and factoring `W = Q Q'` gives
# `k` sweeps along the columns of `Q` plus one more along `L s` for the term
# carrying the mode's own dependence.
#
# Across a *unit* the piece that made that work does not survive unaltered. `C`
# is now the inverse of an arrow matrix and is dense: a direction built from it
# has mass in every member's block, so seeding one would cost a sweep per
# member and the saving would evaporate.
#
# What rescues it is that `dH/dtheta` is block *sparse* even though `C` is not.
# With `A = I - M` the log likelihood's curvature in `u` space,
#
#     psi := tr(C A) = sum_b tr(C[b,b] A[b,b]) + 2 sum_b sum_a tr(C[a,b] A[b,a])
#
# so only the selected entries of `C` appear -- exactly what
# `_laplace_selected_inverse` returns -- and every term involving block `b`
# needs only `b`'s own members. A member outside `b` has no dependence on `u_b`
# at all, so its contribution to both `A[b,b]` and `A[b,a]` is structurally
# zero. Seeding a subject block therefore costs one sweep however large the
# study is, and one level is just the case where every unit is a single block
# with no ancestors and the cross terms vanish.
#
# The cross terms need a *mixed* second derivative rather than a directional
# one. Polarisation identities would give it and are where sign errors live;
# seeding two *independent* duals gives it directly. With
# `x = base + d1*eps1 + d2*eps2` the `eps1*eps2` coefficient is the mixed
# derivative, and it degenerates to the pure one when the directions coincide,
# so one primitive covers both.
#
# Directions are taken in *parameter* space throughout, and that is not a
# detail: a direction in `u` space moves with `L(theta)` while a
# parameter-space one does not, which is what keeps the bookkeeping finite. The
# price is that `L`'s own derivative reappears as an explicit set of trace
# terms, handled at the end of the assembly.

# Distinct tags, so the two seeds are independent nilpotents rather than one
# nested perturbation: `eps1 * eps2` is then the mixed derivative and each
# `eps^2` still vanishes.
struct _LaplaceSeedInner end
struct _LaplaceSeedOuter end

"""
    _laplace_unit_seeded_gradient(laplace, U, values, Ls, u, members, d1, d2, order)

One reverse sweep over the given members, each evaluated at *its own* shifted
parameter vector with two independent seed directions added.

The per-member shift is the whole reason this cannot seed a single shared
vector: members of a unit differ precisely by their random effects. The seed
directions are in raw-parameter space and are the same for every member, which
is what lets one sweep serve them all.

Results are returned *per member* rather than summed, as `npar x length(members)`
matrices whose columns follow `members`. Summing is what the caller usually
wants, but not always: every quantity attached to a block needs the members
under that block and no others, and once summed they cannot be separated again.

`d0` is the ordinary gradient, `d1c` its derivative along the first direction,
and `d12` the mixed derivative along both -- the pure second directional
derivative when the two directions coincide. `order = 1` seeds only the first
direction and leaves `d12` empty.
"""
function _laplace_unit_seeded_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    members, d1::Vector{Float64}, d2::Vector{Float64}, order::Integer,
    slot::Integer=1)

    spec = laplace.spec
    units = laplace.units
    unitmembers = units.members[U]
    npar = length(values)
    nm = length(members)
    empty = zeros(Float64, 0, 0)
    failure = (ok=false, d0=empty, d1c=empty, d12=empty)

    d0 = zeros(Float64, npar, nm)
    d1c = zeros(Float64, npar, nm)

    if order == 1
        seed = ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0)
        S = typeof(seed)
        aws = _laplace_workspace!(laplace, S, npar, slot)
        gradient = Vector{S}(undef, npar)
        x = Vector{S}(undef, npar)
        for (c, m) in enumerate(members)
            shifted = _laplace_member_values(values, spec, Ls, u, units.offsets[U][m])
            @inbounds for t in 1:npar
                x[t] = shifted[t] + seed * d1[t]
            end
            loglik = _laplace_subject_value_gradient!(gradient,
                laplace.objective.subject_objectives[unitmembers[m]], aws, x)
            _laplace_finite(loglik) || return failure
            @inbounds for t in 1:npar
                d0[t, c] = ForwardDiff.value(gradient[t])
                d1c[t, c] = ForwardDiff.partials(gradient[t])[1]
            end
        end
        return (ok=true, d0=d0, d1c=d1c, d12=empty)
    end

    # Two independent nilpotents, one per tag, so the eps1*eps2 coefficient is
    # the mixed derivative directly rather than through a polarisation identity.
    e1 = ForwardDiff.Dual{_LaplaceSeedOuter}(
        ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0),
        ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 0.0))
    e2 = ForwardDiff.Dual{_LaplaceSeedOuter}(
        ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 0.0),
        ForwardDiff.Dual{_LaplaceSeedInner}(1.0, 0.0))
    S = typeof(e1)
    aws = _laplace_workspace!(laplace, S, npar, slot)
    gradient = Vector{S}(undef, npar)
    x = Vector{S}(undef, npar)
    d12 = zeros(Float64, npar, nm)
    for (c, m) in enumerate(members)
        shifted = _laplace_member_values(values, spec, Ls, u, units.offsets[U][m])
        @inbounds for t in 1:npar
            x[t] = shifted[t] + e1 * d1[t] + e2 * d2[t]
        end
        loglik = _laplace_subject_value_gradient!(gradient,
            laplace.objective.subject_objectives[unitmembers[m]], aws, x)
        _laplace_finite(loglik) || return failure
        @inbounds for t in 1:npar
            inner = ForwardDiff.value(gradient[t])
            d0[t, c] = ForwardDiff.value(inner)
            d1c[t, c] = ForwardDiff.partials(inner)[1]
            d12[t, c] = ForwardDiff.partials(ForwardDiff.partials(gradient[t])[1])[1]
        end
    end
    return (ok=true, d0=d0, d1c=d1c, d12=d12)
end

# Deep: `isfinite` on a dual tests only its value, and a sweep's whole point
# is the partials. Testing the value alone is what let a NaN derivative travel
# from a unit whose predicted variance had collapsed all the way into the
# assembled gradient, where it was indistinguishable from a bad trial point.
@inline _laplace_finite(x::Real) = isfinite(x)
@inline _laplace_finite(x::ForwardDiff.Dual) = _finite_deep(x)

"""
    _laplace_seeded_unit_gradient!(out, laplace, U, values, Ls, dL, M, factors, elim)

Accumulate unit `U`'s exact contribution to `dT/dtheta` into `out`, in a number
of sweeps proportional to the unit's members rather than to the parameter count.

    dT/dtheta = dg/dtheta + [ dpsi/dtheta + s' B ] / 2

with `psi = tr(C A)`, `s = C dpsi/du` and `B = d2g/du dtheta` carrying the
mode's own dependence through `duhat/dtheta = C B`. The envelope theorem is what
removes `uhat` from the first term and leaves it only inside `psi`.

Two facts make the assembly finite. Directions are taken in *parameter* space,
because a direction in `u` space moves with `L(theta)` and a parameter-space one
does not; the price is that `L`'s own derivative reappears as the explicit terms
at the end. And every quantity that belongs to a block is built from that
block's members alone -- a member outside block `b` has no dependence on `u_b`,
so it contributes nothing to `A[b,b]` or `A[b,a]`, and folding it in is simply
wrong rather than merely wasteful. That is why the sweeps return per-member
results: summed, they could not be taken apart again.

Returns `false` if any sweep or factorization failed, so the caller can fall
back rather than proceed on a partial answer.
"""
function _laplace_seeded_unit_gradient!(out::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, values::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dL::Vector{Vector{Matrix{Float64}}},
    M::CTSEMBlockMatrix{Float64}, factors, elim, slot::Integer=1)

    spec = laplace.spec
    units = laplace.units
    blocks = units.blocks[U]
    npar = length(values)
    d = units.dims[U]
    nmem = length(units.members[U])
    uhat = laplace.modes[U]
    zerodir = zeros(Float64, npar)

    # No random effects anywhere in this unit: there is no integral, the term
    # is the plain log likelihood, and one ordinary reverse sweep is the whole
    # gradient. Returning early without it would silently contribute nothing.
    if d == 0
        pass = _laplace_unit_seeded_gradient(laplace, U, values, Ls, uhat,
            1:nmem, zerodir, zerodir, 1, slot)
        pass.ok || return false
        @inbounds for m in 1:nmem, t in 1:npar
            out[t] += pass.d0[t, m]
        end
        return true
    end

    Cdiag, Ccoup = _laplace_selected_inverse(factors, elim, blocks)
    # The selected inverse of the curvature is where a badly conditioned unit
    # first shows, and every sweep below is scaled by it.
    if !(all(x -> all(isfinite, x), Cdiag) &&
            all(r -> all(x -> all(isfinite, x), r), Ccoup))
        return false
    end
    sweep = (mm, a1, a2, order) -> _laplace_unit_seeded_gradient(laplace, U,
        values, Ls, uhat, mm, a1, a2, order, slot)

    # Scatter a level-space vector onto the raw parameter vector. The local name
    # must not collide with anything in the enclosing scope: a Julia closure
    # rebinds an enclosing local rather than shadowing it, so assigning `out`
    # here would silently replace the caller's accumulator.
    function scatter(level::Int, w::AbstractVector)
        dir = zeros(Float64, npar)
        rho = spec.levels[level].re_index
        @inbounds for p in eachindex(rho); dir[rho[p]] = w[p]; end
        return dir
    end

    Pm = zeros(Float64, npar, nmem)     # d psi / d v_m
    llvm = zeros(Float64, npar, nmem)   # d loglik_m / d v_m
    Bsm = zeros(Float64, npar, nmem)    # member share of s' B
    seen = falses(nmem)

    for (b, block) in enumerate(blocks)
        l = block.level
        L = Ls[l]
        k = block.size
        k == 0 && continue
        # Diagonal term. tr(C[b,b] A[b,b]) = tr(W H) with W = L C[b,b] L', so a
        # Cholesky of W turns the trace into k pure second directional
        # derivatives whose directions live in parameter space.
        W = _laplace_symmetrise(L * Cdiag[b] * transpose(L))
        F = cholesky(Symmetric(W); check=false)
        issuccess(F) || return false
        Q = Matrix(F.L)
        for r in 1:k
            dir = scatter(l, Q[:, r])
            pass = sweep(block.members, dir, dir, 2)
            pass.ok || return false
            # The second-order sweep is where this fails when it fails: it is a
            # third derivative of the process model once the reverse pass is
            # counted, and a unit whose predicted variance has collapsed can
            # produce a NaN there with a perfectly finite log likelihood and
            # first derivative. Caught here so the caller can take the nested
            # route rather than carry the NaN into the sum.
            (all(isfinite, pass.d12) && all(isfinite, pass.d0)) || return false
            @inbounds for (c, m) in enumerate(block.members)
                for t in 1:npar
                    Pm[t, m] += pass.d12[t, c]
                    seen[m] || (llvm[t, m] = pass.d0[t, c])
                end
                seen[m] = true
            end
        end
        # Cross terms with each ancestor, twice over as the trace requires.
        # tr(C[a,b] A[b,a]) = tr(V H) with V = L_a C[a,b] L_b', which has
        # absorbed *both* Cholesky factors -- so the first direction is a bare
        # basis vector, and applying L to it again would count it twice.
        for (t, a) in enumerate(block.ancestors)
            la = blocks[a].level
            V = Ls[la] * transpose(Ccoup[b][t]) * transpose(L)
            for q in 1:k
                e = zeros(Float64, k); e[q] = 1.0
                pass = sweep(block.members, scatter(l, e), scatter(la, V[:, q]), 2)
                pass.ok || return false
                all(isfinite, pass.d12) || return false
                @inbounds for (c, m) in enumerate(block.members)
                    for tt in 1:npar
                        Pm[tt, m] += 2 * pass.d12[tt, c]
                    end
                end
            end
        end
    end

    # Any member no block covered, which happens only for empty blocks.
    if !all(seen)
        pass = sweep(1:nmem, zerodir, zerodir, 1)
        pass.ok || return false
        @inbounds for m in 1:nmem
            seen[m] && continue
            for t in 1:npar; llvm[t, m] = pass.d0[t, m]; end
        end
    end

    # dpsi/du_b runs through this block's members only, then s = C dpsi/du.
    gradu = zeros(Float64, d)
    for (b, block) in enumerate(blocks)
        l = block.level
        rho = spec.levels[l].re_index
        k = block.size
        k == 0 && continue
        acc = zeros(Float64, k)
        for m in block.members, p in 1:k
            acc[p] += Pm[rho[p], m]
        end
        contribution = transpose(Ls[l]) * acc
        @inbounds for q in 1:k
            gradu[block.offset + q] += contribution[q]
        end
    end
    s = _laplace_block_solve(factors, elim, blocks, gradu)
    all(isfinite, s) || return false

    # s' B: one sweep per block, along that block's share of L s. Summed over
    # blocks this is one directional derivative per member along the total shift
    # its own block and its ancestors impose, which is what the product needs.
    for (b, block) in enumerate(blocks)
        k = block.size
        k == 0 && continue
        sb = [s[block.offset + q] for q in 1:k]
        dir = scatter(block.level, Ls[block.level] * sb)
        pass = sweep(block.members, dir, dir, 1)
        pass.ok || return false
        all(isfinite, pass.d1c) || return false
        @inbounds for (c, m) in enumerate(block.members)
            for t in 1:npar; Bsm[t, m] += pass.d1c[t, c]; end
        end
    end

    # dv/dtheta is the identity away from the population parameters, so for
    # every other parameter the contribution is a plain read-off.
    (all(isfinite, llvm) && all(isfinite, Pm) && all(isfinite, Bsm)) ||
        return false
    @inbounds for j in 1:npar
        acc = 0.0
        for m in 1:nmem
            acc += llvm[j, m] + (Pm[j, m] + Bsm[j, m]) / 2
        end
        out[j] += acc
    end

    # The population parameters move every member's v through L, and move psi
    # through L explicitly. A[b,b] = I - M.diag[b] and A[b,a] = -M.coupling[b][t]
    # give the curvature blocks with no further sweeps, and writing each trace
    # through A rather than through the v-space Hessian needs only one
    # triangular solve X = L \ dL.
    Gb = Vector{Vector{Float64}}(undef, length(blocks))
    Fb = Vector{Vector{Float64}}(undef, length(blocks))
    for (b, block) in enumerate(blocks)
        k = block.size
        rho = spec.levels[block.level].re_index
        Gb[b] = zeros(Float64, k)
        Fb[b] = zeros(Float64, k)
        for m in block.members, p in 1:k
            Gb[b][p] += llvm[rho[p], m] + (Pm[rho[p], m] + Bsm[rho[p], m]) / 2
            Fb[b][p] += llvm[rho[p], m]
        end
    end

    for l in eachindex(spec.levels)
        isempty(dL[l]) && continue
        levelpositions = _laplace_level_positions(spec, l)
        for (t, j) in enumerate(levelpositions)
            X = Ls[l] \ dL[l][t]
            total = 0.0
            for (b, block) in enumerate(blocks)
                k = block.size
                if block.level == l && k > 0
                    ub = [uhat[block.offset + q] for q in 1:k]
                    sb = [s[block.offset + q] for q in 1:k]
                    shift = dL[l][t] * ub
                    sshift = dL[l][t] * sb
                    for p in 1:k
                        total += Gb[b][p] * shift[p] + Fb[b][p] * sshift[p] / 2
                    end
                    KB = Cdiag[b] * (Matrix{Float64}(LinearAlgebra.I, k, k) .- M.diag[b])
                    for x in 1:k, y in 1:k
                        total += KB[x, y] * X[y, x]
                    end
                    # tr(X' A[b,a] C[a,b]): this block's own factor moving.
                    for (tt, a) in enumerate(block.ancestors)
                        Z = (.-M.coupling[b][tt]) * transpose(Ccoup[b][tt])
                        for x in 1:k, y in 1:k
                            total += X[y, x] * Z[y, x]
                        end
                    end
                end
                # tr(C[a,b] A[b,a] X): an ancestor's factor moving.
                for (tt, a) in enumerate(block.ancestors)
                    blocks[a].level == l || continue
                    ka = blocks[a].size
                    Z = transpose(Ccoup[b][tt]) * (.-M.coupling[b][tt])
                    for x in 1:ka, y in 1:ka
                        total += Z[x, y] * X[y, x]
                    end
                end
            end
            out[j] += total
        end
    end
    return true
end

"""Raw positions of one level's population parameters, scales then correlations."""
_laplace_level_positions(spec::CTSEMLaplaceSpec, l::Integer) =
    vcat(spec.levels[l].sd_index, spec.levels[l].cor_index)

"""
    _laplace_level_chol_derivatives(values, spec)

`dL_l / d values[p]` for every level `l` and every population parameter `p` of
that level, in the order `_laplace_level_positions` gives them.

Cheap whatever the model: each `L_l` is a small Cholesky built only from that
level's own scales and correlations, with no data in it. Shared by every unit,
so it is computed once per evaluation rather than once per unit.
"""
function _laplace_level_chol_derivatives(values::AbstractVector{Float64},
    spec::CTSEMLaplaceSpec)
    out = Vector{Vector{Matrix{Float64}}}(undef, length(spec.levels))
    for l in eachindex(spec.levels)
        level = spec.levels[l]
        positions = _laplace_level_positions(spec, l)
        k = length(level.re_index)
        if isempty(positions) || k == 0
            out[l] = Matrix{Float64}[]
            continue
        end
        chol_of = function (p)
            v = convert(Vector{eltype(p)}, values)
            @inbounds for (slot, position) in enumerate(positions)
                v[position] = p[slot]
            end
            return vec(_laplace_popchol(v, level))
        end
        J = ForwardDiff.jacobian(chol_of, values[positions])
        out[l] = [Matrix{Float64}(reshape(collect(view(J, :, t)), k, k))
                  for t in eachindex(positions)]
    end
    return out
end


"""
    _laplace_unit_weights(laplace)

Roughly what each unit costs, for load balancing.

A unit's evaluation walks its members once per block of `u`, and the number of
blocks grows with the members, so the cost is superlinear in unit size. Rows
times members captures both without measuring anything: it is a relative
weight, and only its ordering and rough scale matter to the assignment.
"""
function _laplace_unit_weights(laplace::CTSEMLaplaceObjective)
    subjects = laplace.objective.subject_objectives
    return [begin
        members = laplace.units.members[U]
        rows = sum(size(subjects[i].data, 2) for i in members; init=0)
        Float64(rows) * max(1, length(members))
    end for U in eachindex(laplace.units.members)]
end

"""How often the seeded gradient assembly was abandoned for the nested route."""
const _CTSEM_LAPLACE_FALLBACKS = Ref(0)

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
    # Factorizations rather than dense curvatures. On a study of a few thousand
    # subjects a dense one is over a hundred megabytes, and one is held per unit
    # for the whole evaluation; the factors are kilobytes.
    primal_curvature = Vector{_LaplaceFactorization{Float64}}(undef, nunits)
    # The seeded assembly needs the curvature itself as well as its factors:
    # `A[b,b] = I - M.diag[b]` and `A[b,a] = -M.coupling[b][t]` are where the
    # explicit population-parameter terms come from, and recovering them from
    # the factors would cost more than keeping them. Block form, so this is
    # kilobytes per unit rather than the megabytes a dense one would be.
    primal_matrices = Vector{CTSEMBlockMatrix{Float64}}(undef, nunits)
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
    # Cost-weighted rather than contiguous. A unit's cost is roughly its
    # observations times its members -- the members enter twice, once through
    # the sweeps and once through the block count -- and with studies of
    # different sizes an equal-count split leaves one chunk holding most of the
    # work while the rest wait at the barrier.
    ranges = _ctsem_chunk_assignment(_laplace_unit_weights(laplace), nchunks)
    chunk_ok = fill(true, nchunks)
    chunk_bad = fill(NaN, nchunks)
    run_primal = function (c)
        aws = _laplace_workspace!(laplace, Float64, length(theta), c)
        @inbounds for U in ranges[c]
            _laplace_solve_unit_mode!(laplace, U, theta, Ls, c)
            u = laplace.modes[U]
            blocks = laplace.units.blocks[U]
            M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
                _laplace_unit_curvature(laplace, U, theta, Ls, u, c)
            fac = _laplace_factor_repaired!(M, blocks)
            laplace.mode_repaired[U] = fac.repaired
            ok, logdetM, factors, coupling = fac.ok, fac.logdet, fac.factors, fac.coupling
            primal_curvature[U] = (factors, coupling)
            primal_matrices[U] = M
            inner = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
            term = if !isfinite(inner.value) || isempty(u)
                inner.value
            elseif ok
                inner.value - logdetM / 2
            else
                NaN
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

    # 3. The gradient, by a number of seeded reverse sweeps proportional to a
    #    unit's members rather than to the parameter count.
    #
    # `_laplace_nested_gradient` computes the same thing by running ForwardDiff
    # over the whole per-unit term. It is kept because `test_laplace.jl` checks
    # the two against each other at one, two and three levels: they share the
    # primal and nothing else, so agreement to machine precision is a real
    # check on the seeded assembly, which has a lot of chain rule in it. Set
    # `nested_gradient = true` to use it.
    grad = zeros(Float64, length(theta))
    # One route at every depth. A one-level unit is a single block with no
    # ancestors, so the cross terms and the selected inverse degenerate to
    # nothing and the assembly reduces to the `k + 1` sweeps per subject that
    # the specialised single-level version used to do by hand.
    if !nested_gradient
        dLlevels = _laplace_level_chol_derivatives(theta, laplace.spec)
        # One accumulator per chunk rather than one shared vector: the unit
        # contributions are a sum, and summing per chunk and then across chunks
        # is the same sum in a different order.
        partials = [zeros(Float64, length(theta)) for _ in 1:nchunks]
        fill!(chunk_ok, true)
        run_gradient = function (c)
            @inbounds for U in ranges[c]
                factors, elim = primal_curvature[U]
                if !_laplace_seeded_unit_gradient!(partials[c], laplace, U, theta,
                        Ls, dLlevels, primal_matrices[U], factors, elim, c)
                    chunk_ok[c] = false
                    return nothing
                end
                # A sweep can return success and still have accumulated a
                # non-finite contribution: it reports whether its
                # factorizations worked, not whether the numbers that came out
                # of them are usable. One NaN here is the whole gradient, and
                # the trial point is then rejected with a perfectly good
                # objective value attached to it -- 54 of one fit's 267
                # evaluations went that way, and the fit stopped 0.14 log units
                # short with a gradient of 3.2.
                #
                # The fallback below exists for exactly this and was reachable
                # only through a failed factorization. It computes the same
                # quantity by ForwardDiff over the whole per-unit term, sharing
                # only the primal, so it is a genuinely different route rather
                # than a retry.
                if !all(isfinite, partials[c])
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
            grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_curvature)
            _CTSEM_LAPLACE_FALLBACKS[] += 1
        end
    else
        grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_curvature)
    end
    return (value=value, gradient=grad, subject_loglik=subject_loglik,
        converged=all(laplace.inner_converged))
end

"""
    _laplace_nested_gradient(laplace, theta, Ls, curvature)

The exact outer gradient by ForwardDiff over the whole per-unit term.

`O(npar * k)` reverse sweeps per unit, which is what
`_laplace_seeded_unit_gradient!` exists to avoid. Retained as the oracle the
seeded path is tested against at one, two and three levels, and as its fallback
when a factorization the seeded path needs is not available.
"""
function _laplace_nested_gradient(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}},
    curvature::AbstractVector)
    nunits = length(laplace.units.members)
    total_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        accumulated = zero(S)
        for U in 1:nunits
            uhat = laplace.modes[U]
            ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, curvature[U], wsd)
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
    ctsem_kalman(laplace, values; from_level=1, subject_matrices=true, subject_values=nothing)

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

`subject_values` supplies those per-subject vectors instead of solving for them
here, as an `nsubjects x length(values)` matrix in this objective's subject
order. It exists because a prediction call may filter over different rows than
the fit did -- a subset of subjects, an interpolated time grid, observations
withheld -- and a mode re-solved against those rows is not the fitted one. In
the extreme it is not even close: withhold every observation and there is
nothing left to condition on, so the mode collapses to zero and every subject
silently reverts to the population parameters. The caller holds the fitted
objective and can compute the modes against the data they were fitted to, so
that is where they come from.
"""
function ctsem_kalman(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    from_level::Integer=1, subject_matrices::Bool=true, fields=String[],
    subject_values::Union{Nothing,AbstractMatrix}=nothing)
    persubject = if subject_values === nothing
        ctsem_laplace_subject_values(laplace, values; from_level=from_level)
    else
        nsubjects = length(laplace.objective.subject_objectives)
        size(subject_values, 1) == nsubjects || throw(ArgumentError(
            "subject_values has $(size(subject_values, 1)) rows for " *
            "$nsubjects subjects."))
        size(subject_values, 2) == length(values) || throw(ArgumentError(
            "subject_values has $(size(subject_values, 2)) columns for " *
            "$(length(values)) parameters."))
        subject_values
    end
    return ctsem_kalman(laplace.objective, persubject;
        subject_matrices=subject_matrices, fields=fields)
end

"""
    _laplace_primal_curvature(laplace, theta, Ls)

Solve every unit's mode and factorize its curvature, returning the factors.

Shared by the places that need a primal solve later -- the nested gradient, the
per-unit scores, the mode Jacobian -- so that "solve the modes, then factorize"
is written once.
"""
function _laplace_primal_curvature(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}})
    nunits = length(laplace.units.members)
    out = Vector{_LaplaceFactorization{Float64}}(undef, nunits)
    for U in 1:nunits
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = laplace.modes[U]
        blocks = laplace.units.blocks[U]
        M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
            _laplace_unit_curvature(laplace, U, theta, Ls, u)
        fac = _laplace_factor_repaired!(M, blocks)
        laplace.mode_repaired[U] = fac.repaired
        out[U] = (fac.factors, fac.coupling)
    end
    return out
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
    hessians = _laplace_primal_curvature(laplace, theta, Ls)
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
        blocks = laplace.units.blocks[U]
        factors = nothing; coupling = nothing
        if !isempty(u)
            M = _laplace_unit_curvature(laplace, U, theta, Ls, u)
            fac = _laplace_factor_repaired!(M, blocks)
            factors = fac.factors; coupling = fac.coupling
        end
        # Only the diagonal block of the inverse is wanted, so it is solved for
        # a column at a time rather than by inverting the whole curvature --
        # which for a large study would be the dense matrix this file exists to
        # avoid forming.
        conditional = function (slice)
            out = zeros(Float64, length(slice), length(slice))
            for (t, column) in enumerate(slice)
                e = zeros(Float64, length(u)); e[column] = 1.0
                out[:, t] = _laplace_block_solve(factors, coupling, blocks, e)[slice]
            end
            return out
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
            block = conditional(slice)
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
    mode_repaired=copy(laplace.mode_repaired),
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
    verbose::Bool=false, nested_gradient::Bool=false, tune_chunks::Bool=true,
    lbfgs_memory::Integer=_CTSEM_LBFGS_MEMORY, progress_overwrite::Bool=true,
    progress_callback=nothing, progress::Bool=verbose)
    start_values = collect(Float64, start)
    invalid_objective = floatmax(Float64) / 1e8
    gradient_limit = sqrt(floatmax(Float64))
    # Why trial points were rejected, counted rather than guessed at. A fit
    # that stops short of a stationary point almost always did so because its
    # line search ran out of points it was allowed to accept, and these say
    # which of the three reasons was doing it.
    _CTSEM_LAPLACE_FALLBACKS[] = 0
    rejected_nonfinite = 0
    rejected_inner = 0
    rejected_gradient = 0
    accepted_calls = 0
    # `evaluated`, not `result`: Julia binds an assignment inside a closure to
    # the enclosing scope's local of the same name, and the outer Optim result
    # below is called `result`. Every objective call was overwriting it. That
    # was invisible while the outer one was only read after `Optim.optimize`
    # returned, and stopped being invisible the moment anything read it
    # between runs -- the restart loop below did, and `Optim.minimum` was
    # handed an evaluation named tuple.
    fg! = function (F, G, x)
        evaluated = try
            ctsem_laplace_evaluate(laplace, x; gradient=G !== nothing,
                nested_gradient=nested_gradient)
        catch
            nothing
        end
        objective = evaluated === nothing ? NaN : evaluated.value
        # A unit whose inner Newton did not reach `inner_tol` has not produced
        # the mode the term is defined at, so the number is not the objective
        # -- it is whatever the iteration happened to stop on. Treating the
        # trial point as invalid makes the line search shrink towards a point
        # where the inner problem is solvable, which is the right response;
        # accepting it lets the outer optimizer follow a function that is not
        # a function of theta.
        finite_value = evaluated !== nothing && isfinite(objective)
        inner_ok = evaluated !== nothing && evaluated.converged
        valid = finite_value && inner_ok
        if valid && G !== nothing
            valid = all(isfinite, evaluated.gradient) &&
                all(abs(value) < gradient_limit for value in evaluated.gradient)
            if !valid
                rejected_gradient += 1
                if verbose && rejected_gradient == 1
                    finite_part = filter(isfinite, evaluated.gradient)
                    println("Laplace probe: first gradient rejection, objective ",
                        objective, ", ", count(!isfinite, evaluated.gradient),
                        " of ", length(evaluated.gradient),
                        " entries non-finite, largest finite ",
                        isempty(finite_part) ? 0.0 : maximum(abs, finite_part))
                    println("Laplace probe: at x = ", collect(x))
                    println("Laplace probe: gradient = ", collect(evaluated.gradient))
                end
            end
        end
        finite_value || (rejected_nonfinite += 1)
        (finite_value && !inner_ok) && (rejected_inner += 1)
        valid && (accepted_calls += 1)
        if !valid
            # Zeros, and deliberately not the last valid gradient.
            #
            # Handing back the previous gradient at a new point was tried, on
            # the reasoning that it keeps the search direction pointing back
            # toward the feasible region. It does, and it also feeds L-BFGS a
            # secant pair whose gradient never belonged to that point, which
            # corrupts the curvature history and sends later directions
            # somewhere arbitrary: one fit in ten came back after a single
            # iteration with a NaN gradient, and the next spent half an hour
            # not finishing. The sentinel objective is what makes the line
            # search shrink, and it does that on the value alone.
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        G !== nothing && (G .= -evaluated.gradient)
        return F === nothing ? nothing : -objective
    end
    # See `ctsem_optimize` for why this is a callback rather than `show_trace`.
    # The inner mode count is worth reporting here and not there: a Laplace
    # iteration that is re-solving every unit's mode from scratch costs an order
    # of magnitude more than one that is warm-starting, and the difference shows
    # up as a stall that the objective alone does not explain.
    # See `ctsem_optimize`: progress is not verbosity.
    reporter = CTSEMProgress(progress; label="optimise",
        overwrite=progress_overwrite)
    # Recorded every iteration whatever `verbose` says; see `ctsem_optimize`.
    # `inner` is traced too, because a Laplace fit that stalls usually stalls
    # in the inner solve and the outer objective alone does not show it.
    trace = CTSEMTrace(:objective, :gradient_norm, :inner_converged)
    watcher = CTSEMCallback(progress_callback)
    watch = function (state)
        latest = state isa AbstractVector ? last(state) : state
        inner = count(laplace.inner_converged)
        _record!(trace, latest.iteration, -latest.value, latest.g_norm, inner)
        if _due(reporter)
            _progress_line(reporter, latest.iteration, Int(maxiter),
                @sprintf("logpost %11.2f", -latest.value),
                @sprintf("|g| %9.2e", latest.g_norm),
                @sprintf("inner %d/%d", inner,
                    length(laplace.inner_converged)))
        end
        # Its own cadence; see `ctsem_optimize`.
        _invoke_callback(watcher, latest.iteration, Int(maxiter),
            -latest.value, latest.g_norm)
        return false
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol, f_reltol=f_tol,
        x_abstol=x_tol, show_trace=false, store_trace=false, callback=watch,
        extended_trace=false)
    # Measure the chunk count rather than trusting `cores`. See
    # `ctsem_tune_chunks!`: on small models the wide split is slower than the
    # serial one, by up to 3.7x, and no rule from the model shape alone
    # predicts where the crossover is.
    tuning = tune_chunks ? ctsem_tune_chunks!(
        () -> ctsem_laplace_evaluate(laplace, start_values; gradient=true,
            nested_gradient=nested_gradient); verbose=verbose) : nothing
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
    # `alphaguess` for the same reason `ctsem_optimize` has it, which this
    # route was simply left out of. L-BFGS has no curvature history on its
    # first iteration, so it goes downhill with whatever the initial step guess
    # gives, and Optim's default `InitialStatic()` is an unscaled alpha of one
    # -- a first step as long as the gradient. On this objective the gradient
    # at the starting values is routinely in the hundreds, so the first step
    # was hundreds of units into a region where every ctsem transform is flat
    # to machine precision and the inner mode solve has nothing to work with.
    #
    # It shows up as extreme sensitivity to the starting draw, which is only
    # `rnorm(npar, 0, 0.01)` and cannot itself explain anything. Ten seeds on
    # one 25-subject ordinal model before this line: two converged, six stopped
    # short with gradients between 5.6 and 82, one threw, and one reported
    # convergence 146 log units below the answer.
    #
    # `scaled=true` divides alpha by the gradient norm, so the first step has
    # length one in parameter space however steep the objective is.
    lbfgs = Optim.LBFGS(m=Int(lbfgs_memory),
        alphaguess=Optim.LineSearches.InitialStatic(scaled=true))
    linesearch = "hagerzhang"
    run_from = function (from, method)
        try
            (Optim.optimize(Optim.only_fg!(fg!), from, method, options), true)
        catch err
            err isa InterruptException && rethrow()
            verbose && println("Laplace: Hager-Zhang line search failed (",
                sprint(showerror, err), "); retrying with backtracking")
            (Optim.optimize(Optim.only_fg!(fg!), from,
                Optim.LBFGS(m=Int(lbfgs_memory),
                    alphaguess=Optim.LineSearches.InitialStatic(scaled=true),
                    linesearch=Optim.LineSearches.BackTracking()), options),
                false)
        end
    end
    result, hz = run_from(start_values, lbfgs)
    hz || (linesearch = "backtracking")
    # Hager-Zhang has a second way to fail, and it is quieter than throwing: it
    # runs out of line search, returns whatever iterate it had reached, and
    # Optim reports that as a finished optimisation. Watched on a 60-subject
    # ordinal model, it stopped after two iterations having gone -2503.9 ->
    # -2500.2 -> -2501.3, uphill on the last one, with 68 objective evaluations
    # spent for those two steps and a final gradient of 475. Nothing was
    # rejected and every inner mode solve converged; the objective is simply
    # curved enough here that the Wolfe bracketing gives up. The verdict below
    # calls that not converged, correctly, but a correct verdict on a failed
    # fit is still a failed fit.
    #
    # The `catch` above already treats a Hager-Zhang failure as a reason to
    # switch line searches. This is the same failure without the exception, so
    # it gets the same answer. Backtracking asks only for sufficient decrease,
    # so it cannot fail to bracket -- on the eight fits that reproduced this
    # stall it converged every time, from the same starting values.
    #
    # It resumes from Hager-Zhang's own minimizer rather than the starting
    # values: that point is already downhill, and a fresh L-BFGS there has no
    # stale curvature history from the steps that went wrong.
    #
    # Backtracking is not simply made the default because it is worse when
    # nothing has gone wrong. Armijo alone imposes no curvature condition, so
    # it stops polishing sooner: over those same eight fits it took 47
    # iterations against Hager-Zhang's 36 and finished at gradients around
    # 1e-5 where Hager-Zhang reaches 1e-9. Fast and precise where that works,
    # robust where it does not.
    if hz
        reached = collect(Optim.minimizer(result))
        probe = ctsem_laplace_evaluate(laplace, reached; gradient=true)
        gnorm = isempty(probe.gradient) ? 0.0 : maximum(abs, probe.gradient)
        if !isfinite(gnorm) || gnorm > max(g_tol,
                1e-6 * max(one(gnorm), abs(probe.value)))
            verbose && println("Laplace: Hager-Zhang stopped after ",
                Optim.iterations(result), " iteration(s) with |g| ", gnorm,
                "; continuing with backtracking")
            retry = Optim.optimize(Optim.only_fg!(fg!), reached,
                Optim.LBFGS(m=Int(lbfgs_memory),
                    alphaguess=Optim.LineSearches.InitialStatic(scaled=true),
                    linesearch=Optim.LineSearches.BackTracking()), options)
            after = ctsem_laplace_evaluate(laplace,
                collect(Optim.minimizer(retry)); gradient=true)
            # Only if it actually helped. Keeping the better of the two points
            # means the fallback can never make a fit worse than not having it.
            if isfinite(after.value) && after.value >= probe.value
                result = retry
                linesearch = "hagerzhang+backtracking"
            end
        end
    end

    if verbose
        println("Laplace: ", accepted_calls, " objective evaluations accepted, ",
            rejected_nonfinite, " rejected as non-finite, ", rejected_inner,
            " for an inner mode solve that did not converge, ",
            rejected_gradient, " for the gradient; ",
            _CTSEM_LAPLACE_FALLBACKS[],
            " gradient(s) fell back to the nested route")
        println("Laplace: inner modes ",
            count(laplace.inner_converged), "/", length(laplace.inner_converged),
            " converged, max |dg/dz| ",
            isempty(laplace.inner_gradient) ? 0.0 : maximum(laplace.inner_gradient),
            ", curvature repaired at the mode for ", count(laplace.mode_repaired),
            " unit(s) (", count(laplace.hessian_repaired), " somewhere on the way)")
    end
    minimizer = collect(Optim.minimizer(result))
    final = ctsem_laplace_evaluate(laplace, minimizer; gradient=true)
    # `ctsem_optimize` has always closed its progress line and this route never
    # did, so an in-place update was left open and whatever R printed next
    # landed on the same line -- reported as "inner 100/100Computing exact
    # Hessian".
    progress && _progress_done(reporter,
        @sprintf("%d iterations", Optim.iterations(result)),
        @sprintf("logpost %.4f", final.value))
    # A fit that ends where it started, with a gradient nowhere near zero, has
    # not converged whatever Optim says. Optim's own verdict is the disjunction
    # of three criteria, and a line search that fails on its first try
    # satisfies the `f` one trivially: the objective did not change because
    # nothing was accepted. Two fits in thirteen returned their starting values
    # this way, reporting success, which in a simulation study is silently
    # wrong rather than loudly broken.
    moved = isempty(minimizer) ? 0.0 : maximum(abs, minimizer .- start_values)
    gradient_norm = isempty(final.gradient) ? 0.0 : maximum(abs, final.gradient)
    stalled = moved == 0 && (!isfinite(final.value) || gradient_norm > max(g_tol, 1e-6))
    # Optim's `g_tol` is an *absolute* bound on the gradient, and a log
    # likelihood of order 1e3 puts 1e-8 out of reach however good the fit is --
    # L-BFGS runs out of line search first and reports nothing converged. So
    # convergence is also allowed on a criterion that scales with the problem.
    #
    # This is deliberately not a loosening of the strict test, because that is
    # how the failure this replaced stayed hidden: on one 12-subject model the
    # old code stopped after a single iteration with a gradient of 1.2e9 and a
    # log likelihood of -7.7e6, and reported success. That point fails the
    # scaled criterion by eight orders of magnitude.
    scaled_tolerance = max(g_tol, 1e-6 * max(one(gradient_norm), abs(final.value)))
    # `isfinite(gradient_norm)` explicitly: `Optim.g_converged` can be true at
    # a point whose gradient is NaN, since it was set on an earlier iterate,
    # and `NaN <= tolerance` is false so the scaled test alone would not have
    # caught it. One draw in ten reported convergence with a NaN gradient.
    finite_gradient = isfinite(gradient_norm)
    # The same guard `ctsem_optimize` has, which this route was left out of.
    #
    # Every ctsem transform is flat to machine precision by |raw| ~ 20: the
    # exponential underflows and the derivative is *exactly* zero. A parameter
    # that walks out there therefore reports a zero gradient, satisfies any
    # tolerance, and is indistinguishable from an optimum -- while being pinned
    # by the transform's floating-point limit rather than by the data. A
    # variance going to zero is the usual way in.
    saturated = isempty(minimizer) ? false :
        maximum(abs, minimizer) >= _CTSEM_SATURATION[]
    converged_enough = isfinite(final.value) && finite_gradient &&
        gradient_norm <= scaled_tolerance
    # Convergence needs a small gradient, and nothing else counts as one.
    #
    # `Optim.converged` is the disjunction of its x, f and g criteria, and the
    # first two are satisfied by a line search that stops making progress: the
    # step went to zero, so x did not move and f did not change. That is the
    # signature of giving up, not of arriving. Measured on a 25-subject ordinal
    # model, one starting draw in ten stopped after three iterations with a
    # gradient of 681 and a log likelihood 146 units below the optimum, and
    # `Optim.converged` said true. `stalled` did not catch it because the
    # optimizer had moved.
    #
    # `g_converged` is kept as an alternative to the scaled test because it is
    # a genuine gradient criterion; it is just an absolute one, and `g_tol` is
    # out of reach on a log likelihood of order 1e3 however good the fit.
    verbose && stalled && println("Laplace: the optimizer made no progress from ",
        "its starting values; reporting this as not converged")
    verbose && saturated && println("Laplace: a parameter reached ",
        maximum(abs, minimizer), " on the unconstrained scale, where its ",
        "transform is flat to machine precision; reporting this as not converged")
    verbose && !stalled && !(finite_gradient &&
        (Optim.g_converged(result) || converged_enough)) &&
        println("Laplace: the optimizer stopped with a largest gradient of ",
            gradient_norm, " against a tolerance of ", scaled_tolerance,
            "; reporting this as not converged")
    return (
        minimizer=minimizer,
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        iterations=Optim.iterations(result),
        f_calls=Optim.f_calls(result),
        g_calls=Optim.g_calls(result),
        linesearch=linesearch,
        stalled=stalled,
        chunks=ctsem_max_chunks().max_chunks,
        chunk_timings=tuning === nothing ? Tuple{Int,Float64}[] : tuning.timings,
        gradient_norm=gradient_norm,
        scaled_tolerance=scaled_tolerance,
        saturated=saturated,
        converged=!stalled && !saturated && finite_gradient &&
            (Optim.g_converged(result) || converged_enough),
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
        inner_converged=all(laplace.inner_converged),
        inner_failures=count(!, laplace.inner_converged),
        rejected_nonfinite=rejected_nonfinite,
        rejected_inner=rejected_inner,
        rejected_gradient=rejected_gradient,
        accepted_calls=accepted_calls,
        inner_worst_gradient=isempty(laplace.inner_gradient) ? 0.0 :
            maximum(laplace.inner_gradient),
        inner_iterations=copy(laplace.inner_iterations),
        hessian_repaired=copy(laplace.hessian_repaired),
        mode_repaired=copy(laplace.mode_repaired),
        trace=_trace_result(trace),
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
    primal_curvature = _laplace_primal_curvature(laplace, theta, Ls)

    terms_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        out = Vector{S}(undef, nunits)
        for U in 1:nunits
            uhat = laplace.modes[U]
            ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, primal_curvature[U], wsd)
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

"""
    ctsem_laplace_effect_layout(laplace)

Where each entry of the flat random-effect vector comes from.

The sampler returns the effects as one vector per draw, laid out unit by unit
and, within a unit, by block offset. With a single level that degenerates to one
unit per subject with `k` effects each, which the R side can and does work out
for itself. With more than one level it cannot: a unit's blocks interleave a
group's own effects with its members', and the arrangement is decided here.

Rather than have R reimplement `_laplace_build_units` and keep the two in step,
the engine says what it did. One entry per position, giving the unit, the level
the block belongs to, the block's first member as a global subject index, how
many members the block covers, and which effect within the block it is -- from
which R can attach a subject or a group label and the parameter's own name.

`level` is 1 for the subject level and increases outwards. A block covering more
than one member is a grouping block, and belongs to the group those members
share at that level.
"""
function ctsem_laplace_effect_layout(laplace::CTSEMLaplaceObjective)
    units = laplace.units
    nunits = length(units.members)
    position = Int[]
    unit = Int[]
    level = Int[]
    first_member = Int[]
    nmembers = Int[]
    within = Int[]
    base = 0
    for U in 1:nunits
        for block in units.blocks[U]
            # `offset` is zero-based inside the unit and units concatenate in
            # order, which is the arithmetic every other consumer of these
            # blocks already does.
            for k in 1:block.size
                push!(position, base + block.offset + k)
                push!(unit, U)
                push!(level, block.level)
                push!(first_member, isempty(block.members) ? 0 :
                    units.members[U][block.members[1]])
                push!(nmembers, length(block.members))
                push!(within, k)
            end
        end
        base += units.dims[U]
    end
    order = sortperm(position)
    return (position=position[order], unit=unit[order], level=level[order],
        first_member=first_member[order], nmembers=nmembers[order],
        within=within[order])
end
