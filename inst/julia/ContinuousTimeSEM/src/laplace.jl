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
    CTSEMLaplaceLevel(re_index, sd_index, cor_index, sd_scale, group, ngroups;
                      covmatcode = 0)

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
  * `covmatcode`: which covariance construction this level's population matrix
    uses, in `sdcovsqrt2cov`'s own encoding -- 0 the row-normalised correlation
    square root, 1 a factor, 2 Fisher z inside a matrix exponential. Per level
    because each level's population covariance is its own matrix.
  * `rank`: how many dimensions this level's population covariance spans.
    Equal to `k` for the ordinary full-rank level, in which case `sd_index` and
    `cor_index` describe it as above and `load_index` is empty. Below `k` it is
    a reduced-rank level: `load_index` holds the `k*rank - rank*(rank-1)/2`
    entries of a `k x rank` loading matrix in column-major lower-triangular
    order, `sd_index` and `cor_index` are empty, and the covariance is
    `L * L'`.

    A reduced level keeps a `k`-dimensional latent block rather than a
    `rank`-dimensional one, and that is deliberate. The deviation a unit
    contributes is `L * u` for `u ~ N(0, I)` (see `_laplace_member_values!`),
    so the `k - rank` columns of `L` that are zero make those coordinates of
    `u` invisible to the likelihood. Their inner mode is the origin and their
    inner curvature is exactly the prior identity, so they integrate out
    against it contributing nothing -- no singular covariance is ever formed
    and nothing is inverted that cannot be. The cost is `k - rank` wasted
    coordinates per unit; the alternative, resizing every block, would touch
    the offsets, the ancestor sets and the curvature factorisation, where the
    same number means both "parameters this level moves" and "dimensions its
    block spans" and separating the two silently is how a wrong answer gets
    made here.
"""
struct CTSEMLaplaceLevel
    re_index::Vector{Int}
    sd_index::Vector{Int}
    cor_index::Vector{Int}
    load_index::Vector{Int}
    sd_scale::Vector{Float64}
    group::Vector{Int}
    ngroups::Int
    covmatcode::Int
    rank::Int

    function CTSEMLaplaceLevel(re_index, sd_index, cor_index, sd_scale, group,
        ngroups::Integer; covmatcode::Integer=0, rank::Integer=-1,
        load_index=Int[])
        re = Vector{Int}(collect(re_index))
        sd = Vector{Int}(collect(sd_index))
        cor = Vector{Int}(collect(cor_index))
        load = Vector{Int}(collect(load_index))
        scale = Vector{Float64}(collect(sd_scale))
        grp = Vector{Int}(collect(group))
        k = length(re)
        # A negative rank means "not stated", which is the full-rank level every
        # caller before this argument existed builds.
        r = rank < 0 ? k : Int(rank)
        (r >= 0 && r <= k) || throw(ArgumentError(
            "rank must lie in 0:$(k) for a level with $(k) random effects, got $(r)"))
        length(scale) == k ||
            throw(DimensionMismatch("one sdscale per random effect is required"))
        if r < k
            expected = k * r - div(r * (r - 1), 2)
            length(load) == expected || throw(DimensionMismatch(
                "expected $(expected) loading parameters for rank $(r) over $(k) random effects, got $(length(load))"))
            # Refused rather than ignored: a level carrying both descriptions
            # would leave which one the covariance came from decided by the
            # order of two branches, which is the kind of thing that stays
            # wrong for months because both produce a perfectly ordinary matrix.
            (isempty(sd) && isempty(cor)) || throw(ArgumentError(
                "a reduced-rank level is described by load_index alone; sd_index and cor_index must be empty"))
        else
            length(sd) == k ||
                throw(DimensionMismatch("one population scale parameter per random effect is required"))
            expected = div(k * (k - 1), 2)
            length(cor) == expected || throw(DimensionMismatch(
                "expected $(expected) correlation parameters for $(k) random effects, got $(length(cor))"))
            isempty(load) || throw(ArgumentError(
                "load_index applies to a reduced-rank level only; this level is full rank"))
        end
        allunique(re) || throw(ArgumentError("random-effect parameter indices must be distinct within a level"))
        isempty(grp) || (minimum(grp) >= 1 && maximum(grp) <= ngroups) ||
            throw(ArgumentError("group ids must lie in 1:ngroups"))
        return new(re, sd, cor, load, scale, grp, Int(ngroups), Int(covmatcode), r)
    end
end

"""True when this level's population covariance spans fewer dimensions than it
has random effects."""
isreducedrank(level::CTSEMLaplaceLevel) = level.rank < length(level.re_index)

"""Number of model parameters this level moves."""
nrandomeffects(level::CTSEMLaplaceLevel) = length(level.re_index)

"""
Dimensions this level's latent block spans.

The same as `nrandomeffects` for an ordinary level and equal to the rank for a
reduced one, and keeping the two apart is the whole of what makes a reduced
level cheap. They were one number until a rank could differ from an effect
count, and every site that used it meant one or the other: `re_index[p]` is
indexed by the parameter count, `u[offset + q]` by the block dimension. A site
that takes the wrong one still runs and still returns a number.
"""
nlatent(level::CTSEMLaplaceLevel) = level.rank

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
    nsubjects::Integer=0; covmatcode::Integer=0)
    group = collect(1:Int(nsubjects))
    return CTSEMLaplaceSpec([CTSEMLaplaceLevel(re_index, sd_index, cor_index,
        sd_scale, group, Int(nsubjects); covmatcode=covmatcode)])
end

"""Total number of random effects across every level, per subject."""
nrandomeffects(spec::CTSEMLaplaceSpec) = sum(nrandomeffects(l) for l in spec.levels; init=0)

nlevels(spec::CTSEMLaplaceSpec) = length(spec.levels)

"""Total latent dimensions across every level, per subject."""
nlatent(spec::CTSEMLaplaceSpec) = sum(nlatent(l) for l in spec.levels; init=0)

"""True when any level's covariance spans fewer dimensions than it has effects."""
hasreducedrank(spec::CTSEMLaplaceSpec) = any(isreducedrank, spec.levels)

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
                # Latent dimensions, not parameters: this is how wide the
                # block is in `u`.
                k = nlatent(levels[l])
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
                k = nlatent(levels[l])
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
                nlatent(levels[outer]) == 0 && continue
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

The inner modes are retained between calls, but they do not start the next
evaluation's Newton solve: that starts at the origin, so the value is a
function of theta alone. See `_laplace_solve_unit_mode!`, and
`ctsem_set_warm_start!` to retain them as starts and measure what that costs.
"""
mutable struct CTSEMLaplaceObjective{O} <: CTSEMOptimisable
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
    # Whether the prior floor bound for this unit at the last evaluation; see
    # `_laplace_prior_floor_logdet`. Not a fault -- a legitimately wide
    # posterior is floored too -- but the reported term is a bound rather than
    # the Laplace value there, and nothing else says so.
    logdet_floored::Vector{Bool}
    # The parameter vector of the last evaluation, so that report-time
    # diagnostics (`ctsem_laplace_conditioning`) describe the same point the
    # modes and flags above do. A copy the evaluation already made; nothing is
    # computed from it per evaluation.
    last_values::Vector{Float64}
    # The prior floor this objective evaluates under, `:total` or `:gated`,
    # and the gated rule's hand-off band; see `ctsem_set_laplace_floor!`. On
    # the objective rather than the session, so that no fit changes how
    # another evaluates.
    floor::Symbol
    gate_lo::Float64
    gate_hi::Float64
    # Units the gated rule flagged at the last evaluation.
    gated_units::Int
end

function CTSEMLaplaceObjective(objective::CTSEMObjective, spec::CTSEMLaplaceSpec;
    inner_maxiter::Integer=_LAPLACE_INNER_MAXITER[], inner_tol::Real=1e-10,
    floor=:total, gate_lo::Real=0.2, gate_hi::Real=0.7)
    Symbol(floor) in (:total, :gated) ||
        throw(ArgumentError("the Laplace floor must be :total or :gated"))
    0 < gate_lo < gate_hi || throw(ArgumentError("the gate needs 0 < lo < hi"))
    nsubjects = length(objective.subject_objectives)
    units = _laplace_build_units(spec, nsubjects)
    nunits = length(units.members)
    return CTSEMLaplaceObjective{typeof(objective)}(objective, spec, units,
        [zeros(Float64, units.dims[U]) for U in 1:nunits],
        Int(inner_maxiter), Float64(inner_tol), [Dict{Any,Any}()],
        zeros(Int, nunits), zeros(Float64, nunits), falses(nunits), falses(nunits),
        falses(nunits), falses(nunits), Float64[], Symbol(floor), Float64(gate_lo),
        Float64(gate_hi), 0)
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
    group=Int[], level_ngroups=Int[], level_covmatcode=Int[],
    level_rank=Int[], load_index=Int[],
    inner_maxiter::Integer=_LAPLACE_INNER_MAXITER[], inner_tol::Real=1e-10,
    floor="total", gate_lo::Real=0.2, gate_hi::Real=0.7)
    nsubjects = length(objective.subject_objectives)
    counts = isempty(level_nre) ? [length(re_index)] : Vector{Int}(Int.(level_nre))
    ngroups = isempty(level_ngroups) ? [nsubjects] : Vector{Int}(Int.(level_ngroups))
    length(counts) == length(ngroups) || throw(DimensionMismatch(
        "level_nre and level_ngroups must describe the same number of levels"))
    groups = isempty(group) ? collect(1:nsubjects) : Vector{Int}(Int.(group))
    length(groups) == nsubjects * length(counts) || throw(DimensionMismatch(
        "group must hold one entry per subject per level"))
    # The model's own construction unless a level names its own. Read from the
    # objective rather than defaulted to 0, so a covmattransform='z' model does
    # not get a Laplace population distribution built the other way -- which is
    # a wrong answer with no symptom, since both constructions return a
    # perfectly ordinary covariance matrix.
    codes = if isempty(level_covmatcode)
        fill(objective.params.covmatcode, length(counts))
    else
        Vector{Int}(Int.(level_covmatcode))
    end
    length(codes) == length(counts) || throw(DimensionMismatch(
        "level_covmatcode must hold one construction code per level"))

    # One rank per level, defaulting to that level's own `k`, which is the
    # full-rank description every caller before this existed sends.
    ranks = if isempty(level_rank)
        copy(counts)
    else
        Vector{Int}(Int.(level_rank))
    end
    length(ranks) == length(counts) || throw(DimensionMismatch(
        "level_rank must hold one rank per level"))

    levels = CTSEMLaplaceLevel[]
    # Four cursors rather than one reused for two things: a reduced level
    # consumes loadings and no scales or correlations, so the scale cursor and
    # the effect cursor stop advancing together the moment one level is reduced.
    re_at = 0; sd_at = 0; cor_at = 0; load_at = 0
    for l in eachindex(counts)
        k = counts[l]
        r = ranks[l]
        reduced = r < k
        ncor = reduced ? 0 : div(k * (k - 1), 2)
        nsd = reduced ? 0 : k
        nload = reduced ? k * r - div(r * (r - 1), 2) : 0
        push!(levels, CTSEMLaplaceLevel(
            Int.(re_index[(re_at + 1):(re_at + k)]),
            Int.(sd_index[(sd_at + 1):(sd_at + nsd)]),
            Int.(cor_index[(cor_at + 1):(cor_at + ncor)]),
            Float64.(sd_scale[(re_at + 1):(re_at + k)]),
            groups[((l - 1) * nsubjects + 1):(l * nsubjects)],
            ngroups[l]; covmatcode=codes[l], rank=r,
            load_index=Int.(load_index[(load_at + 1):(load_at + nload)])))
        re_at += k; sd_at += nsd; cor_at += ncor; load_at += nload
    end
    return CTSEMLaplaceObjective(objective, CTSEMLaplaceSpec(levels);
        inner_maxiter=inner_maxiter, inner_tol=inner_tol, floor=Symbol(floor),
        gate_lo=gate_lo, gate_hi=gate_hi)
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
                              ("correlation", level.cor_index),
                              ("loading", level.load_index))
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

`base` is exactly Stan's `rawpopcovbase`: the transformed scales on the
diagonal, the capped correlation coordinates in the lower triangle. Handing it
to `sdcovsqrt2cov` is what makes this the same construction the filter and the
stan program use rather than a second implementation of it -- this function
used to spell the code=0 route out by hand, and so silently ignored the model's
`covmattransform`.
"""
function _laplace_popchol(values::AbstractVector{T}, level::CTSEMLaplaceLevel) where {T}
    k = nrandomeffects(level)
    k == 0 && return zeros(T, 0, 0)
    isreducedrank(level) && return _laplace_poploading(values, level)
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
                # No squash here any more: `constraincorsqrt1` applies it.
                # The cap still bounds the correlation at 0.99, because it
                # bounds the coordinate the squash then maps.
                base[i, j] = _laplace_cap_correlation(
                    values[level.cor_index[counter]])
            end
        end
    end
    # `scales` already carries the 1e-10 floor its own transform applies, and
    # nothing is added on top of it here: one floor, in the diagonal element's
    # transform, matching the stan path.
    # `_sdcovsqrt2cov_uncached!` rather than `sdcovsqrt2cov`: the same
    # construction, without the content-keyed cache in front of it. This
    # function runs under nested ForwardDiff -- the outer gradient needs third
    # derivatives -- and a cache whose key comparison is `!=` on Duals compares
    # values and not partials. One construction per outer evaluation per level
    # is not a path worth caching anyway.
    buffer = _make_square_buffer(T, k)
    _sdcovsqrt2cov_uncached!(buffer, base, level.covmatcode, Val(k))
    return Matrix(cholesky(Symmetric(buffer.out, :L)).L)
end

"""
    _laplace_poploading(values, level)

A reduced-rank level's factor: `k x k`, with the loading in its first `rank`
columns and zeros in the rest, so that `L * L'` is the population covariance
and `L * u` is the deviation for `u ~ N(0, I)`.

The zero columns are what let the block stay `k`-dimensional; see the note on
`rank` in `CTSEMLaplaceLevel`. No Cholesky is taken, because there is nothing
to decompose -- the parameters *are* the factor, which is also why they need no
positivity transform and no correlation cap. `sd_scale` multiplies a row, so a
level still scales its own spread exactly as the full-rank form does, and a
loading is signed: negating a whole column leaves `L * L'` alone, so read the
covariance rather than the sign of one entry.
"""
function _laplace_poploading(values::AbstractVector{T}, level::CTSEMLaplaceLevel) where {T}
    k = nrandomeffects(level)
    r = level.rank
    L = zeros(T, k, r)
    counter = 0
    @inbounds for q in 1:r
        for p in q:k
            counter += 1
            L[p, q] = values[level.load_index[counter]] * level.sd_scale[p]
        end
    end
    return L
end

"""Every level's Cholesky factor, innermost first."""
_laplace_popchols(values::AbstractVector{T}, spec::CTSEMLaplaceSpec) where {T} =
    [_laplace_popchol(values, level) for level in spec.levels]

"""Single-level shorthand. The fitting path takes every level through
`_laplace_popchols`; only `test_laplace.jl` reaches for this one."""
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
    _laplace_saturated_parameters(laplace, values; threshold=_CTSEM_TRANSFORM_FLOOR[])

Which raw parameter indices have a materialising transform that has stopped
responding at `values`, over the whole raw vector `ctsem_laplace_optimize`
works in.

Two layers, because two different pieces of code materialise raw values here:
ordinary parameter-table coordinates (including each level's `re_index`) go
through `laplace.objective.params.regular_transforms` and are covered by
`_ctsem_saturated_parameters` (`parameter_transforms.jl`); each level's
`sd_index` and `cor_index` coordinates go through the population-scale and
-correlation formulas `_laplace_popchol` evaluates directly
(`log1p_exp(2x-1) * sdscale + 1e-10` and `2/(1+exp(-clamp(x))) - 1`), which sit
outside `regular_transforms` entirely -- "a tail of the vector no model matrix
cell reads", per `CTSEMLaplaceLevel`'s docstring -- so they need their own
derivative here. A correlation coordinate sitting past
`ctsem_set_correlation_cap!`'s cap is flagged by the same mechanism, because
the clamp itself has zero derivative past it: the same practical fact
`ctsem_laplace_boundary` already diagnoses separately, by name.

A raw index reached by neither layer -- a TI-predictor coefficient -- is never
flagged, as in `_ctsem_saturated_parameters`.
"""
function _laplace_saturated_parameters(laplace::CTSEMLaplaceObjective,
        values::AbstractVector{T}; threshold::Real=_CTSEM_TRANSFORM_FLOOR[]) where {T<:Real}
    found = Set(_ctsem_saturated_parameters(laplace.objective.params, values,
        eachindex(values); threshold=threshold))
    for level in laplace.spec.levels
        for j in eachindex(level.sd_index)
            idx = level.sd_index[j]
            scale = level.sd_scale[j]
            d = abs(ForwardDiff.derivative(
                raw -> log1p_exp(2 * raw - 1) * scale + 1e-10, values[idx]))
            d < threshold && push!(found, idx)
        end
        for idx in level.cor_index
            d = abs(ForwardDiff.derivative(
                raw -> 2 / (1 + exp(-_laplace_cap_correlation(raw))) - 1, values[idx]))
            d < threshold && push!(found, idx)
        end
    end
    return sort!(collect(found))
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
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        base = offsets[l]
        L = Ls[l]
        # `p` walks the parameters this level moves, `q` the dimensions of its
        # block. They are the same number only for a full-rank level.
        for p in 1:k
            acc = zero(S)
            for q in 1:r
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
effects. The fitting path shifts through `_laplace_member_values`, which takes
any number of levels; this one is reached only by `test_laplace.jl`, as the
readable statement of what that shift is at one level.
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

# The worker pool.
#
# Everything in this file that runs in parallel runs on one pool of at most
# `cores` workers, and a worker's identity is *ambient*: it lives in task-local
# storage as a base slot and a budget, and `_laplace_parallel` splits that budget
# among the workers it claims. Nothing is passed down and so nothing can be
# dropped on the way.
#
# That is the point. The previous arrangement threaded a slot and a count through
# eight functions as two separate arguments with independent defaults. A caller
# could supply half of a band, and the half whose absence is *safe* -- the count,
# whose default means "run serially" -- is the one that got dropped. Twice. It
# cost 4.6x down to 1.2x on the primal, with the answer correct throughout, so no
# test could report it. An ambient budget cannot be half-supplied, and a caller
# that forgets to touch it inherits its parent's, which is always within bounds.
#
# Nesting is bounded by construction: a region hands each worker `budget / nw`,
# so the workers alive at any moment never outnumber the pool however deep the
# nesting goes. That is also the `cores` ceiling the caller asked for, held by
# arithmetic rather than by a tuner multiplying two numbers and hoping.
"""Workers the pool may use: the `cores` ceiling, capped by the threads there are."""
@inline function _laplace_pool_size()
    requested = _CTSEM_MAX_CHUNKS[]
    return max(1, requested == 0 ? Threads.nthreads() :
        min(requested, Threads.nthreads()))
end


# Free workspace slots, as a stack behind a lock.
#
# A lock and a `Vector{Int}` rather than anything cleverer: acquisition happens
# a handful of times per parallel region, against subject filters that take
# milliseconds, so contention is irrelevant and being obviously correct is not.
const _LAPLACE_FREE_LOCK = ReentrantLock()
const _LAPLACE_FREE = Int[]

"""The pool size the free stack was built for; 0 until it is built."""
const _LAPLACE_POOL_N = Ref(0)

"""
    _laplace_slot_try_acquire()

A free slot, or `nothing` if none is free *right now*.

Never blocks, and that is the whole safety argument. A region always makes
progress on the slot its own task already holds, so a caller that gets nothing
here runs the work itself rather than waiting for a resource that a task
downstream of it may be holding. There is no cycle to deadlock on.
"""
function _laplace_slot_try_acquire()
    # An unlocked peek first. This reads a length while another task may be
    # pushing or popping, so the answer can be stale -- and both ways of being
    # wrong are harmless. A false "empty" costs a helper this region did not
    # have to take; a false "non-empty" still takes the lock and rechecks.
    # What it buys is the common case: once the pool is exhausted, every region
    # below asks and is refused, and those refusals should not serialise on a
    # lock.
    isempty(_LAPLACE_FREE) && return nothing
    lock(_LAPLACE_FREE_LOCK) do
        isempty(_LAPLACE_FREE) ? nothing : pop!(_LAPLACE_FREE)
    end
end

"""Hand a slot back. Called as soon as a worker runs out of items, not at the
end of the region -- see `_laplace_parallel` for why that timing is the point."""
function _laplace_slot_release(slot::Int)
    lock(_LAPLACE_FREE_LOCK) do
        push!(_LAPLACE_FREE, slot)
    end
    return nothing
end

"""How wide the pool is. An upper bound for sizing per-slot storage, never a
promise that this many workers are available."""
@inline _laplace_pool_width() = max(1, _laplace_pool_size())

"""The slot this task's scratch lives in. Slot 1 unless a pool worker set it."""
@inline _laplace_slot() = get(task_local_storage(), :ctsem_slot, 1)::Int

"""
Whether this task was handed its band by a parallel region, rather than taking
the whole pool as the outermost caller.

Set only by `_laplace_parallel` and `_laplace_partition`, on the tasks they
spawn, and never cleared -- a task is one or the other for its whole life.
"""
@inline _laplace_is_worker() =
    get(task_local_storage(), :ctsem_pool_worker, false)::Bool

"""
Refuse a slot the store does not have, by name.

Restoring a guard that was removed as "unnecessary under the pool". It was not:
the thing it catches is an arithmetic error in the claim, and the pool moved
that arithmetic rather than removing it. Without this the failure is a
`BoundsError` inside a spawned task, wrapped in a `TaskFailedException`, with
nothing in it naming the pool.
"""
@inline function _laplace_check_slot(laplace::CTSEMLaplaceObjective, slot::Integer)
    1 <= slot <= length(laplace.workspaces) && return nothing
    error("laplace worker pool: slot $slot asked for, " *
          "$(length(laplace.workspaces)) exist. The caller reached a parallel " *
          "region without `_laplace_ensure_pool!`, or claimed more workers " *
          "than its band holds.")
end

"""
    _laplace_ensure_pool!(laplace)

Fix the pool size and make sure there is a scratch store per worker.

Called by `_laplace_parallel` and `_laplace_partition`, which are the only two
things that spawn and so the only two that need the store. It used to be the
responsibility of every entry point instead, and three of them did not do it;
each was a slot past the end of a store nothing had sized, and each cost a
suite run to find. A few entry points still call it directly, which is
harmless -- it is idempotent -- and says at the top of the function what state
the rest of it assumes.

Idempotent, and a no-op for a pool worker: a worker returns the band it was
lent and never touches the store, so only the outermost task ever grows the
vector and there is no concurrent `push!` on an array its siblings are
indexing.
"""
function _laplace_ensure_pool!(laplace::CTSEMLaplaceObjective)
    # A worker returns the band it was lent and touches nothing. Only the
    # outermost task grows the vector, so there is no concurrent `push!` on it
    # -- which would be a data race on the backing array while its siblings are
    # indexing it.
    _laplace_is_worker() && return _laplace_pool_width()
    n = _laplace_pool_size()
    while length(laplace.workspaces) < n
        push!(laplace.workspaces, Dict{Any,Any}())
    end
    # The band too, not just the store. Reached only by a non-worker, which is
    # the outermost task by definition, so it takes the whole pool.
    #
    # "Has a band already" was the first test for that and is wrong: a band set
    # on the main task stays there for the rest of the session, so a later call
    # with a different `cores` sized the store to the new pool and kept the old
    # band -- in one direction a band wider than the store, which
    # `_laplace_check_slot` then refused. Being a worker is the thing actually
    # being asked about, so it is what gets recorded.
    task_local_storage(:ctsem_slot, 1)
    # Slot 1 belongs to this task, the root of the evaluation; the rest are the
    # pool.
    #
    # **Rebuilt only when the pool size changes**, and that condition is the
    # whole of the correctness argument. "Only a non-worker reaches here, so no
    # helper can hold a slot" was the first version and it is false: the root
    # task is not a worker, and it runs items itself, so an item whose body
    # opens a nested region brings the root back through here *while its own
    # helpers are still holding slots*. Refilling then hands the same slot to
    # two tasks, which surfaces as two subject evaluations sharing one adjoint
    # tape -- a `BoundsError` indexing a record vector, from inside a spawned
    # task, with nothing in it naming the pool.
    #
    # Slots come back through a `finally`, so the stack is full again between
    # evaluations and there is nothing to reclaim.
    if _LAPLACE_POOL_N[] != n
        lock(_LAPLACE_FREE_LOCK) do
            empty!(_LAPLACE_FREE)
            for slot in n:-1:2
                push!(_LAPLACE_FREE, slot)
            end
            _LAPLACE_POOL_N[] = n
        end
    end
    return n
end


"""
    _laplace_partition(f, laplace, n)

Call `f(mine, w)` once per worker, with a disjoint stride of `1:n` and the
worker's own index, so a reduction can keep one accumulator per worker and sum
them after the join.

For work whose per-item cost is small but whose *setup* is not. The seeded
assembly walks every block of the unit for the members it owns; giving it one
member at a time would be that same walk once per member.

Unlike `_laplace_parallel` this spends the whole budget on one axis, because
the members of a unit are the finest axis there is -- there is nothing below
for a sub-budget to buy.
"""
function _laplace_partition(f, laplace::CTSEMLaplaceObjective, n::Int)
    n <= 0 && return true
    # Sizing the store is this function's job, not its callers'.
    #
    # It used to be "called at the top of every entry point", and three entry
    # points did not: `_laplace_unit_hessian`, `ctsem_sample_start` and
    # `ctsem_sample_metric`, each of which is reached on an objective this
    # session may never have evaluated -- `ctSample()` on a reloaded fit builds
    # a fresh one. The band lives on the task and the store lives on the
    # object, so a task carrying a band from an earlier object spawned workers
    # that indexed past this one's empty store. Each was found by a separate
    # forty-minute suite run.
    #
    # Spawning is the only thing that needs the store, and these two functions
    # are the only things that spawn, so putting it here makes the omission
    # unrepresentable rather than merely detected.
    _laplace_ensure_pool!(laplace)
    n == 1 && return f(1:1:1, 1) !== false
    # Take what is idle, then divide among what was actually taken. The count
    # is discovered rather than computed, so nothing here has to know what an
    # item costs -- see `_laplace_parallel` for the argument.
    #
    # Slots are held for the whole region here, unlike there, because each
    # worker gets a fixed stride and works to completion. It is the *choice* of
    # how many workers that has become dynamic, not the division among them.
    slots = Int[]
    while length(slots) < n - 1
        slot = _laplace_slot_try_acquire()
        slot === nothing && break
        push!(slots, slot)
    end
    nw = length(slots) + 1
    nw == 1 && return f(1:1:n, 1) !== false
    oks = fill(true, nw)
    try
        Threads.@sync begin
            for (k, slot) in enumerate(slots)
                Threads.@spawn begin
                    task_local_storage(:ctsem_pool_worker, true)
                    task_local_storage(:ctsem_slot, slot)
                    oks[k + 1] = f((k + 1):nw:n, k + 1) !== false
                end
            end
            # The caller takes a stride too, on the slot it already holds.
            oks[1] = f(1:nw:n, 1) !== false
        end
    finally
        for slot in slots
            _laplace_slot_release(slot)
        end
    end
    return all(oks)
end

"""
    _laplace_parallel(f, laplace, items)

Run `f(item)` for every item, on the pool, and report whether all succeeded.

`f` returning `false` stops that worker and fails the whole region, which is
how the failure paths here already report. Items are taken from a shared
counter rather than dealt out in advance, so one slow item cannot leave the
rest of the pool idle.

A budget of one runs serially, and so does a single item -- which is what makes
nesting free: the innermost regions of a saturated pool cost nothing but a
branch.
"""
function _laplace_parallel(f, laplace::CTSEMLaplaceObjective, items)
    _laplace_ensure_pool!(laplace)      # see `_laplace_partition`
    return _laplace_parallel(f, items)
end

"""
    _laplace_parallel(f, items)

The same, for a region whose caller has already sized the store.

One caller: the chain runner in `sample_run.jl`, which is generic over a
density closure and has no objective to hand -- deliberately, since that is
what lets it run the same chains for the marginal and the state-explicit
targets. `ctsem_sample_marginal` sizes the store before any chain starts, so
the region below is never the outermost one.

This is the hole the objective argument closes everywhere else, kept open here
on purpose and watched by `_laplace_check_slot`, which turns a slot past the
end into a named error rather than a `BoundsError` inside a task.
"""
function _laplace_parallel(f, items)
    n = length(items)
    n == 0 && return true
    n == 1 && return f(@inbounds items[1]) !== false

    # Ask for one helper before building anything a lone worker would not need.
    #
    # The shared counter, the failure flag and the task vector are three heap
    # allocations and two atomic operations *per item*, and a region with no
    # helper needs none of them -- it is a plain loop. That distinction was
    # missing in the first version of this, which routed the alone case through
    # the same machinery: measured on a thirteen-study model it cost 6% at one
    # worker and 32% at five, which is the whole of this design's overhead and
    # none of its purpose.
    first_slot = _laplace_slot_try_acquire()
    if first_slot === nothing
        @inbounds for i in 1:n
            f(items[i]) === false && return false
        end
        return true
    end

    next = Threads.Atomic{Int}(0)
    # `Atomic`, not a plain `Bool`: several workers may set it, and a shared
    # `Bool` written by several tasks is a race even when they all write the
    # same value.
    failed = Threads.Atomic{Bool}(false)
    helpers = Task[]
    # Whatever is idle right now, and never wait for it.
    #
    # This is the whole change from a computed band. Nothing here predicts how
    # much work an item is, because nothing needs to: a worker that draws a
    # cheap item runs out of items, releases its slot, and the workers still
    # inside an expensive one pick it up at their next region. A study that
    # turns out to sit in an awkward part of the parameter space ends up with
    # more of the pool for exactly as long as it deserves it, and the cost
    # model that used to decide this in advance -- and was wrong by about a
    # factor of two on the thirteen-study affect model -- is gone.
    #
    # `n - 1` because this task is a worker too. Leaving the rest of the pool
    # alone matters: regions nested inside these items need slots, and a loop
    # that grabbed everything would starve them.
    # The slot goes in as an *argument*, so each task captures its own value.
    # A `while` loop reassigning one `slot` variable and closing over it would
    # give every task the same binding -- the last one written -- which is the
    # closure-rebinding trap this file has been caught by twice before, and
    # here it would hand one workspace to several workers.
    start_worker = function (mine::Int)
        Threads.@spawn begin
            task_local_storage(:ctsem_pool_worker, true)
            task_local_storage(:ctsem_slot, mine)
            try
                _laplace_pull!(f, items, n, next, failed)
            finally
                # As soon as the queue is empty, not at the end of the region.
                # That timing is what lets a finished worker's slot reach the
                # unit that is still going.
                _laplace_slot_release(mine)
            end
        end
    end
    push!(helpers, start_worker(first_slot))
    while length(helpers) < n - 1
        slot = _laplace_slot_try_acquire()
        slot === nothing && break
        push!(helpers, start_worker(slot))
    end
    # The caller works too, on the slot it already holds, so the region always
    # completes even when the pool is empty.
    _laplace_pull!(f, items, n, next, failed)
    foreach(wait, helpers)
    return !failed[]
end

"""Pull items off the shared counter until they run out or one fails."""
@inline function _laplace_pull!(f, items, n::Int, next, failed)
    while !failed[]
        i = Threads.atomic_add!(next, 1) + 1
        i > n && break
        if f(@inbounds items[i]) === false
            failed[] = true
            break
        end
    end
    return nothing
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
    nvalues::Integer, slot::Integer=_laplace_slot()) where {T}
    # `slot` is the chunk index, not the thread id. A task can migrate between
    # threads at any yield point, so thread-indexed mutable scratch is a race
    # rather than an optimisation -- the same reason `_get_or_init_adjoint_
    # workspaces!` hands one workspace to each chunk.
    _laplace_check_slot(laplace, slot)
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
    _laplace_scratch_vector!(laplace, V, n, slot, tag)

A reusable vector of length `n`, cached per slot under `tag`, like the adjoint
and filter workspaces beside it.

The member loops used to allocate their working vectors *per member*: the
subject gradient, the seeded input and the shifted parameter vector, three
fresh arrays for every subject of every sweep. That was invisible while the
LAPACK lock dominated threading and is not now. Measured on dev1 at 20 threads,
collection is 28% of a serial gradient and 51% of a sixteen-way one, and it
scales 2.5x where the work scales 6.7x -- so it is collection that sets the
ceiling, and allocation is what collection costs.

Per slot, not per thread, and never shared: a slot belongs to one task for the
duration of a call, which is the contract the worker pool enforces -- a task
holds a band of slots and lends disjoint sub-bands to anything it spawns, so no
two live tasks can name the same one.
"""
function _laplace_scratch_vector!(laplace::CTSEMLaplaceObjective, ::Type{V},
    n::Integer, tag::Symbol, slot::Integer=_laplace_slot()) where {V}
    _laplace_check_slot(laplace, slot)
    store = laplace.workspaces[slot]
    key = (tag, V, Int(n))
    cached = get(store, key, nothing)
    cached === nothing || return cached::Vector{V}
    built = Vector{V}(undef, Int(n))
    store[key] = built
    return built
end

"""
    _laplace_scratch_matrix!(laplace, V, nrow, ncol, tag)

A reusable `nrow` by `ncol` matrix, cached per slot under `tag`.

For results whose shape is fixed by the model rather than by the call --
principally the curvature's per-block jacobian, which is as tall as the whole
unit's random-effect dimension and as wide as one block.
"""
function _laplace_scratch_matrix!(laplace::CTSEMLaplaceObjective, ::Type{V},
    nrow::Integer, ncol::Integer, tag::Symbol,
    slot::Integer=_laplace_slot()) where {V}
    _laplace_check_slot(laplace, slot)
    store = laplace.workspaces[slot]
    key = (tag, V, Int(nrow), Int(ncol))
    cached = get(store, key, nothing)
    cached === nothing || return cached::Matrix{V}
    built = Matrix{V}(undef, Int(nrow), Int(ncol))
    store[key] = built
    return built
end

"""
    _laplace_scratch_levels!(laplace, S, Ls, tag)

`Ls` converted to scalar type `S`, in buffers cached per slot.

The curvature's per-block jacobian used to rebuild its whole dual working set on
every evaluation -- the parameter vector, the level factors, and `uu`, which is
as long as all the unit's random effects together. ForwardDiff calls that body
once per chunk of every block, so on a two-hundred-subject unit it was hundreds
of megabytes per Newton step, and the mode solve takes many. Measured at 918 MB
of the 1.44 GB a single gradient allocated, which is 63% of it, and collection
was 58% of the wall clock. Copying into a buffer costs a pass over the vector;
allocating one costs that plus the collector.
"""
function _laplace_scratch_levels!(laplace::CTSEMLaplaceObjective, ::Type{S},
    Ls::Vector{<:AbstractMatrix}, tag::Symbol) where {S}
    _laplace_check_slot(laplace, _laplace_slot())
    store = laplace.workspaces[_laplace_slot()]
    key = (tag, S, size.(Ls))
    cached = get(store, key, nothing)
    if cached === nothing
        cached = [Matrix{S}(undef, size(L)) for L in Ls]
        store[key] = cached
    end
    out = cached::Vector{Matrix{S}}
    @inbounds for i in eachindex(Ls)
        copyto!(out[i], Ls[i])
    end
    return out
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
    slot::Integer=_laplace_slot()) where {T}
    _laplace_check_slot(laplace, slot)
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
    # `aws.ekf_ws` rather than the subject's own cached workspace, and
    # `ekf_workspace` overrides both.
    #
    # A workspace depends only on the shared `EKFParameters` and the element
    # type, never on the subject, so the adjoint workspace's is the right one
    # for every subject this task filters -- and being a field of a typed struct
    # it arrives here *concrete*, where the accessor hands back an abstract
    # `ContinuousEKFWorkspace{T}` and makes the forward call below a dynamic
    # dispatch. It is also per task, where the subject's cached one is shared
    # between tasks that filter the same subject, which only a sampler does.
    #
    # An override must have the same type for the same reason, and it does: the
    # sampler builds it from the same `params` and the same `T`.
    ws = ekf_workspace === nothing ? aws.ekf_ws : ekf_workspace::typeof(aws.ekf_ws)
    tape = _tape_reset!(aws.tape)
    resize!(aws.tipreds, length(subject_objective.tipreds))
    copyto!(aws.tipreds, subject_objective.tipreds)
    aws.frechet_count = 0
    deferred = aws.defer_frechet
    aws.defer_frechet = false
    try
        local bf, br
        bf = _laplace_mark()
        loglik = _extended_kalman_filter_continuous!(ws, values,
            subject_objective.data, subject_objective.timesteps,
            subject_objective.params, subject_objective.tdpreds,
            subject_objective.tipreds, subject_objective.subject,
            subject_objective.max_timestep, tape)
        _laplace_charge!(_LAPLACE_BYTES_FWD, bf)
        isfinite(loglik) || return loglik
        br = _laplace_mark()
        fill!(aws.theta_bar, zero(T))
        _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
        fill!(gradient, zero(T))
        _ctsem_parameter_layer!(gradient, aws.theta_bar, tape.subject_values,
            subject_objective.params, aws, subject_objective.tipreds, values)
        _laplace_charge!(_LAPLACE_BYTES_REV, br)
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

"""The stand-in a failed or order-1 sweep returns where a matrix would go."""
const _LAPLACE_NO_MATRIX = zeros(Float64, 0, 0)

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

"""
    _laplace_scratch_blockmatrix!(laplace, T, U, blocks, tag)

A block matrix shaped for unit `U`, cached per slot under `tag`.

The curvature builds one of these per Newton step, and the mode solve takes
many steps per evaluation -- on two hundred subjects that is two hundred small
matrices plus their couplings, thrown away and rebuilt each time. Every entry
is written before the matrix is read, block by block over every level, so
reusing the storage cannot carry a value forward; if a block fails, the caller
abandons the whole curvature rather than reading a partial one.

Keyed by unit, because the shape is the unit's. A worker holds one per unit it
is given, and the pool strides units across workers deterministically, so the
total held is the unit count rather than units times workers.
"""
function _laplace_scratch_blockmatrix!(laplace::CTSEMLaplaceObjective, ::Type{T},
    U::Integer, blocks::Vector{CTSEMLaplaceBlock}, tag::Symbol) where {T}
    _laplace_check_slot(laplace, _laplace_slot())
    store = laplace.workspaces[_laplace_slot()]
    key = (tag, T, Int(U))
    cached = get(store, key, nothing)
    cached === nothing || return cached::CTSEMBlockMatrix{T}
    built = CTSEMBlockMatrix(T, blocks)
    store[key] = built
    return built
end

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

"""
Take the block form of a dense symmetric matrix.

The small-unit route in `_laplace_unit_curvature` returns through here, so this
is not test-only however much it reads like it.
"""
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
# The engine's own factorization, not `LinearAlgebra.Cholesky`: this runs once
# per block per Newton step per subject, inside the threaded subject loop, and
# LAPACK's per-call lock is what stopped that loop from threading (see
# small_linalg.jl). The blocks are the size of a subject's random effects.
const _LaplaceCholesky{T} = CTSEMCholesky{T,Matrix{T}}

"""One unit's factorized curvature: the per-block Choleskys and the eliminated
couplings, exactly what `_laplace_block_solve` and `_laplace_selected_inverse`
consume. Named because it is held one per unit for a whole evaluation and a
`Vector{Any}` of them makes every solve a dynamic dispatch."""
const _LaplaceFactorization{T} =
    Tuple{Vector{_LaplaceCholesky{T}},Vector{Vector{Matrix{T}}}}

"""
    _laplace_block_factor(M, blocks)

Eliminate the blocks innermost-first, returning
`(ok, logdet, factors, coupling)`.

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
        f = _ctsem_cholesky(Matrix{T}(_laplace_symmetrise(diag[b])), size(diag[b], 1))
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
    _laplace_selected_inverse(factors, elim, blocks)

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
    aws, positions, into::Union{Nothing,AbstractVector}=nothing) where {T}
    spec = laplace.spec
    units = laplace.units
    members = units.members[U]

    # One member's contribution, into the accumulators the caller hands it.
    # Unlike the seeded assembly these are *reductions*, not per-member columns:
    # every member of a block adds into that block's slice of `u`, which is
    # exactly the coupling that makes the unit one integral. So a task gets its
    # own `inner` and its own total, and they are summed after the join.
    one! = function (m, ws, grad, shift, acc_inner)
        local i, offsets, shifted, loglik, l, level, k, r, base, L, q, acc, pp
        i = members[m]
        offsets = units.offsets[U][m]
        shifted = _laplace_member_values!(shift, values, spec, Ls, u, offsets)
        loglik = _laplace_subject_value_gradient!(grad,
            laplace.objective.subject_objectives[i], ws, shifted)
        isfinite(loglik) || return loglik
        @inbounds for l in eachindex(spec.levels)
            level = spec.levels[l]
            k = nrandomeffects(level)
            r = nlatent(level)
            (k == 0 || r == 0) && continue
            base = offsets[l]
            L = Ls[l]
            for q in 1:r
                acc = zero(T)
                for pp in 1:k
                    acc += L[pp, q] * grad[level.re_index[pp]]
                end
                acc_inner[base + q] += acc
            end
        end
        return loglik
    end

    pos = collect(positions)
    npos = length(pos)
    # "Is anyone free", not "is the pool wide". Those differ whenever the pool
    # is fully committed above this call -- which is the normal state of a
    # nested region -- and asking the wrong one takes the parallel path with no
    # helpers: four per-slot arrays allocated, a partition entered, a single
    # worker doing all of it, and the reduction run anyway. Measured at 43% on
    # a thirteen-study model at five workers.
    if _laplace_pool_width() <= 1 || npos < 2 || isempty(_LAPLACE_FREE)
        # The caller may lend its accumulator. This vector is the whole unit's
        # random-effect dimension in dual arithmetic -- about 144 KB on two
        # hundred subjects -- and the curvature allocates one per block per
        # Newton step while writing nine of its entries. Measured at 318 MB of
        # a 1.21 GB gradient. The Newton path still allocates, because it holds
        # one gradient across a line search while computing the next.
        inner = into === nothing ? zeros(T, length(u)) : fill!(into, zero(T))
        total = zero(T)
        gradient = _laplace_scratch_vector!(laplace, T, length(values), :loglik_grad)
        shift = _laplace_scratch_vector!(laplace, T, length(values), :loglik_shift)
        for m in pos
            loglik = one!(m, aws, gradient, shift, inner)
            isfinite(loglik) ||
                return (value=loglik, gradient=fill(T(NaN), length(u)))
            total += loglik
        end
        return (value=total, gradient=inner)
    end

    # A reduction, not per-member columns: every member of a block adds into
    # that block's slice of `u`, which is exactly the coupling that makes the
    # unit one integral. So each worker keeps its own and they are summed after.
    # Keyed by the slot a worker holds, not by a dense worker index, and
    # taken from per-slot scratch rather than allocated here.
    #
    # Two things follow. The worker count no longer has to be known before the
    # region starts, which is what lets it be whatever the pool can spare. And
    # the accumulators stop being built per call: this used to allocate one
    # vector of the unit's whole random-effect dimension per worker, every
    # time, which is the shape of allocation the rest of this file spent a
    # night removing.
    #
    # `fill(false, ...)` and not `falses`: a `BitVector` packs eight flags to a
    # byte, so two workers setting different slots would read-modify-write the
    # same word and one would be lost.
    width = _laplace_pool_width()
    used = fill(false, width)
    totals = zeros(T, width)
    bad = zeros(T, width)
    okk = fill(true, width)
    _laplace_partition(laplace, npos) do mine, _w
        local ws, grad, shift, acc, tot, i, loglik, slot
        slot = _laplace_slot()
        # `::typeof(aws)` and not a bare fetch. `_laplace_workspace!` cannot
        # promise a concrete type -- the filter workspace's dimensions come
        # from the model at run time, so the best it can say is
        # `CTSEMAdjointWorkspace{T}` -- and calling `one!` with an abstractly
        # typed workspace is a dynamic dispatch per member, which boxes the
        # log likelihood it returns. `aws` is this function's own argument and
        # therefore concrete, and every slot's workspace for the same `T` and
        # the same objective has the same type, so naming it costs nothing and
        # throws loudly if that ever stops being true.
        ws = _laplace_workspace!(laplace, T, length(values))::typeof(aws)
        grad = _laplace_scratch_vector!(laplace, T, length(values), :loglik_grad)
        shift = _laplace_scratch_vector!(laplace, T, length(values), :loglik_shift)
        acc = _laplace_scratch_vector!(laplace, T, length(u), :primal_inner)
        fill!(acc, zero(T))
        used[slot] = true
        tot = zero(T)
        for i in mine
            loglik = one!(pos[i], ws, grad, shift, acc)
            if !isfinite(loglik)
                okk[slot] = false
                bad[slot] = loglik
                return false
            end
            tot += loglik
        end
        totals[slot] = tot
        return true
    end
    @inbounds for w in 1:width
        okk[w] || return (value=bad[w], gradient=fill(T(NaN), length(u)))
    end
    inner = into === nothing ? zeros(T, length(u)) : fill!(into, zero(T))
    @inbounds for w in 1:width
        used[w] || continue
        acc = _laplace_scratch_vector!(laplace, T, length(u), :primal_inner, w)
        for j in eachindex(inner)
            inner[j] += acc[j]
        end
    end
    return (value=sum(totals), gradient=inner)
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
    aws, into::Union{Nothing,AbstractVector}=nothing) where {T}
    result = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u, aws,
        eachindex(laplace.units.members[U]), into)
    isfinite(result.value) || return result
    inner = result.gradient
    @inbounds for a in eachindex(u)
        inner[a] -= u[a]
    end
    return (value=result.value - dot(u, u) / 2, gradient=inner)
end

"""
Whether a unit's inner Newton starts from the mode the last evaluation left, or
from the origin every time. Off, and see `_laplace_solve_unit_mode!` for why.
"""
const _LAPLACE_WARM_START = Ref(false)

"""
Whether to charge allocated bytes to a phase of the evaluation.

Off, and a `const false` rather than a `Ref`, so every `if _LAPLACE_COUNT_BYTES`
below folds away at compile time and the counters cost nothing at all. Flip it
to `true` and rebuild to get the table back.

It stays in the tree because the question it answers -- which phase allocates --
has no other instrument here. The sampled allocation profiler attributes 96.6%
of the bytes to "unattributed", because the engine is precompiled into a shared
library and its stack frames carry `.so` paths that the profiler cannot resolve
to source lines. That is not a setting to change; it is what a pkgimage is.

Read the counters serially. `Base.gc_bytes()` is process-wide, so a delta taken
around one worker's phase counts every other worker's allocation as well, and
the parts then sum to several times the whole -- 5.7 GB against a 1.2 GB total,
on the first attempt. `ctsem_set_max_chunks!(1)` first, and check the sum.
"""
const _LAPLACE_COUNT_BYTES = false

const _LAPLACE_BYTES_MODE = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_CURV = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_OBJ = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_GRAD = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_FWD = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_SWEEP = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_SELINV = Threads.Atomic{Int}(0)
const _LAPLACE_BYTES_REV = Threads.Atomic{Int}(0)
function ctsem_phase_bytes()
    _LAPLACE_COUNT_BYTES || error("phase byte counting is compiled out: set " *
        "`_LAPLACE_COUNT_BYTES = true` in src/laplace.jl and rebuild. " *
        "Reporting zeros here would read as an evaluation that allocated " *
        "nothing.")
    return _ctsem_phase_bytes()
end

_ctsem_phase_bytes() = (mode=_LAPLACE_BYTES_MODE[], curv=_LAPLACE_BYTES_CURV[],
    obj=_LAPLACE_BYTES_OBJ[], grad=_LAPLACE_BYTES_GRAD[],
    fwd=_LAPLACE_BYTES_FWD[], rev=_LAPLACE_BYTES_REV[],
    sweep=_LAPLACE_BYTES_SWEEP[], selinv=_LAPLACE_BYTES_SELINV[])
function ctsem_phase_reset!()
    _LAPLACE_BYTES_MODE[] = 0
    _LAPLACE_BYTES_CURV[] = 0
    _LAPLACE_BYTES_OBJ[] = 0
    _LAPLACE_BYTES_GRAD[] = 0
    _LAPLACE_BYTES_FWD[] = 0
    _LAPLACE_BYTES_REV[] = 0
    _LAPLACE_BYTES_SWEEP[] = 0
    _LAPLACE_BYTES_SELINV[] = 0
    return nothing
end
export ctsem_phase_bytes, ctsem_phase_reset!
@inline _laplace_mark() = _LAPLACE_COUNT_BYTES ? Int(Base.gc_bytes()) : 0

@inline function _laplace_charge!(counter, before)
    _LAPLACE_COUNT_BYTES || return nothing
    Threads.atomic_add!(counter, Int(Base.gc_bytes()) - before)
    return nothing
end

# Temporary: is the member loop actually dividing, and how wide?

"""Temporary: which check in the seeded assembly refused, readable from R."""
const _LAPLACE_DIAG = Ref(0)
ctsem_laplace_diag() = _LAPLACE_DIAG[]
ctsem_laplace_diag_reset!() = (_LAPLACE_DIAG[] = 0)
export ctsem_laplace_diag, ctsem_laplace_diag_reset!


"""
How many Newton iterations a unit's inner solve is allowed, by default.

Fifty was the figure under the warm start, where it was a *continuation*
budget: a unit that ran out carried its iterate into the next evaluation and
went on from there, so the cap was per call and the solve had no cap at all.
From the origin it is the whole solve, and one unit of the count model in
`test-julia-count.R` needs 53 -- it stalled at `|dg/du|` 2.4e-3 on the fiftieth
and reached 1.1e-12 three iterations later.
"""
const _LAPLACE_INNER_MAXITER = Ref(200)

"""
    ctsem_set_inner_maxiter!(n)

Set the default inner Newton iteration budget. Returns the previous value.
Objectives already built keep the budget they were built with.
"""
function ctsem_set_inner_maxiter!(n::Integer)
    previous = _LAPLACE_INNER_MAXITER[]
    _LAPLACE_INNER_MAXITER[] = Int(n)
    return previous
end

export ctsem_set_inner_maxiter!

"""
    ctsem_set_warm_start!(warm)

Start each unit's inner Newton from the retained mode (`true`) or from the
origin (`false`, the default). Returns the previous setting.

Here to measure what the warm start costs, not as a setting to fit with: it
makes the objective depend on the order the points were visited in, which is
what `_laplace_solve_unit_mode!` describes.
"""
function ctsem_set_warm_start!(warm::Bool)
    previous = _LAPLACE_WARM_START[]
    _LAPLACE_WARM_START[] = warm
    return previous
end

export ctsem_set_warm_start!

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
    _laplace_unit_hessian(laplace, U, values, Ls, u)

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
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T}
    ) where {T}
    # An entry point in its own right: tests and the dense curvature route both
    # reach it directly, and it divides inside. `_laplace_ensure_pool!` is a
    # no-op for a task that already holds a band, so calling it here is right
    # whether this is the outermost call or one worker of an outer region.
    #
    # The `slot::Integer=1` that used to sit here was doing nothing at all --
    # the body takes the ambient slot -- and an argument accepted and ignored is
    # worse than no argument, because a caller reads it as a promise.
    _laplace_ensure_pool!(laplace)
    d = length(u)
    d == 0 && return zeros(T, 0, 0)
    inner_of = function (uu)
        S = eltype(uu)
        ws = _laplace_workspace!(laplace, S, length(values))
        vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
        copyto!(vs, values)
        Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
        return _laplace_unit_objective_gradient(laplace, U, vs, Lss, uu, ws,
            _laplace_scratch_vector!(laplace, S, length(uu), :curv_inner)
            ).gradient
    end
    H = ForwardDiff.jacobian(inner_of, collect(u))
    return (H .+ transpose(H)) ./ 2
end

"""
    _laplace_unit_curvature(laplace, U, values, Ls, u)

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
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T}
    ) where {T}
    blocks = laplace.units.blocks[U]
    base = _laplace_scratch_vector!(laplace, T, length(u), :curv_base)
    copyto!(base, u)
    if length(blocks) <= 1 || length(u) < _LAPLACE_BLOCK_THRESHOLD[]
        gradient_of = function (uu)
            S = eltype(uu)
            ws = _laplace_workspace!(laplace, S, length(values))
            vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
            copyto!(vs, values)
            Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
            return _laplace_unit_loglik_gradient(laplace, U, vs, Lss, uu, ws,
                eachindex(laplace.units.members[U]),
                _laplace_scratch_vector!(laplace, S, length(uu), :curv_inner)
                ).gradient
        end
        A = ForwardDiff.jacobian(gradient_of, base)
        dense = Matrix{T}(LinearAlgebra.I, length(u), length(u)) .-
            _laplace_symmetrise(A)
        return _laplace_block_of(dense, blocks)
    end
    M = _laplace_scratch_blockmatrix!(laplace, T, U, blocks, :curv_M)
    # One block's fill. `M.diag[b]` and `M.coupling[b][*]` belong to that block
    # alone, so two blocks never write the same memory. What they can share is a
    # *subject*: an outer block's members are all of its descendants'. Blocks of
    # one level cover disjoint members, so a level at a time is safe and levels
    # are not.
    fill_block! = function (b::Int)
        # `local` for the reason spelled out in `_laplace_seeded_unit_gradient!`:
        # a name assigned here that is also a local of the enclosing function
        # would be shared by every task, silently.
        local block, columns, block_of, J, t, a, rows, xb
        block = blocks[b]
        columns = (block.offset + 1):(block.offset + block.size)
        block_of = function (ub)
            local S, ws, vs, Lss, uu, tt, c
            S = eltype(ub)
            ws = _laplace_workspace!(laplace, S, length(values))
            vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
            copyto!(vs, values)
            Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
            uu = _laplace_scratch_vector!(laplace, S, length(base), :curv_uu)
            copyto!(uu, base)
            @inbounds for (tt, c) in enumerate(columns)
                uu[c] = ub[tt]
            end
            return _laplace_unit_loglik_gradient(laplace, U, vs, Lss, uu, ws,
                block.members,
                _laplace_scratch_vector!(laplace, S, length(base), :curv_inner)
                ).gradient
        end
        # `jacobian!` into a cached result, and views out of it. The jacobian
        # is as tall as the unit's whole random-effect dimension because that
        # is what the gradient body returns, while the block reads its own
        # `block.size` rows and its ancestors' -- fifteen of eighteen hundred
        # on a two-hundred-subject unit. Allocating it fresh cost 130 KB per
        # block per Newton step, and the two row slices below copied again.
        # Computing the whole thing is unavoidable, since forward mode fills
        # every row of the dual result whether or not it is read; allocating
        # it repeatedly is not.
        xb = _laplace_scratch_vector!(laplace, T, block.size, :curv_ub)
        @inbounds for (t, c) in enumerate(columns)
            xb[t] = base[c]
        end
        J = _laplace_scratch_matrix!(laplace, T, length(base), block.size,
            :curv_jac)
        ForwardDiff.jacobian!(J, block_of, xb)
        # M = I - d2(sum ll)/du du, block by block.
        M.diag[b] .= .-(@view J[columns, :])
        @inbounds for t in 1:block.size
            M.diag[b][t, t] += one(T)
        end
        for (t, a) in enumerate(block.ancestors)
            rows = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
            M.coupling[b][t] .= .-transpose(@view J[rows, :])
        end
        return nothing
    end

    # A level with many blocks divides over blocks; a level with one divides
    # over that block's members instead. A unit's innermost level is one block
    # per subject, which the member axis cannot touch at all, and its outermost
    # is a single block spanning every member, which the block axis cannot touch
    # -- so neither axis alone is enough and the choice is per level.
    # A level at a time, because blocks of one level cover disjoint members and
    # blocks of different levels do not. Within a level the pool takes them;
    # a level holding a single block spends the whole budget on that block's
    # members instead, which is automatic -- the region below inherits whatever
    # this one did not claim. A unit's innermost level is one block per subject,
    # which only this axis can spread, and its outermost is a single block over
    # every member, which only the member axis can.
    at = 1
    while at <= length(blocks)
        stop = at
        while stop < length(blocks) && blocks[stop + 1].level == blocks[at].level
            stop += 1
        end
        _laplace_parallel(laplace, at:stop) do b
            fill_block!(b)
        end
        at = stop + 1
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
    _laplace_solve_unit_mode!(laplace, U, values, Ls)

Newton's method on `g_U`, from the origin -- or from the retained mode for
unit `U` when `ctsem_set_warm_start!` has asked for that.

`g_U` is a log likelihood minus a quadratic, so its curvature is negative
definite near the mode and Newton is the right method; away from the mode, and
for a nonlinear process model, it need not be. Two guards, both visible in the
diagnostics rather than silent: a curvature that is not negative definite is
shifted until it is and the unit flagged as repaired, and a step that does not
improve `g_U` is halved before the iteration gives up.
"""
function _laplace_solve_unit_mode!(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}})
    d = laplace.units.dims[U]
    aws = _laplace_workspace!(laplace, Float64, length(values))
    retained = _LAPLACE_WARM_START[] ? copy(laplace.modes[U]) : zeros(Float64, d)
    started_warm = any(!iszero, retained)
    warm = _laplace_newton_unit_mode(laplace, U, values, Ls, aws, retained)
    best = warm
    # The start is the origin, and that is not a detail of the iteration: it is
    # what makes the value a function of theta.
    #
    # Newton warm-started from the mode the last evaluation left is the obvious
    # economy, and it is sound exactly when `g_U` has one mode -- which it does
    # for a linear-Gaussian model, where the integrand is Gaussian in `u` and
    # the inner problem is concave. It is not, for a nonlinear one. Measured on
    # a 50-subject model with a random `-log1p_exp(-param)` drift, at one fixed
    # parameter vector, alternating with distant evaluations:
    #
    #     from the origin   -3234.201758   every time, |dg/du| 4.1e-9
    #     warm              -3234.201758
    #     warm              -3387.582975   same theta, |dg/du| 4.4e-6
    #     warm              -3354.568832
    #
    # Every warm value is a different stationary point and every one is worse.
    # An outer optimiser handed that is not maximising a function: its own trace
    # climbed to -3234.20 and the re-evaluation of that same minimizer returned
    # -3349.14, and over eight seeded starts three reached the optimum, one
    # stopped at a log likelihood of -1.5e7 with a gradient of 2.4e15. From the
    # origin, eight of eight, all to -2959.05.
    #
    # The origin is the population mean, it is the same point for every caller
    # and always inside the support, and it costs the inner iterations the warm
    # start would have saved: 16% on a 200-subject linear-Gaussian fit where the
    # warm start is safe -- identical estimate, identical outer iterations -- and
    # less than nothing on the model above, where a poisoned mode costs the outer
    # optimiser far more than the inner solve ever saved.
    #
    # `_LAPLACE_WARM_START` turns it back on, to measure that. When it is on, a
    # warm start that fails is retried from the origin and the better of the two
    # kept, which is the guard that used to stand here on its own.
    if !warm.converged && started_warm
        cold = _laplace_newton_unit_mode(laplace, U, values, Ls, aws,
            zeros(Float64, d))
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

So comfortably that it is worth knowing what loosening it would buy, since the
floor is the only lever on the mode solve's cost: it decides whether a unit
takes one Newton step or two, and the line search accepts the full step every
time, so a solve is one member sweep at the origin plus one per step. Where the
random effects enter through a transform, the first step leaves `|dg/du|`
around 5e-6 and the second drives it to 6e-15 -- ten million times past what is
asked -- for a third sweep over every member of every unit.

Measured on dev1, fitting a 200-subject one-latent model simulated from known
parameters to convergence, each setting twice in both orders, at 200 units of
one subject and at 20 units of ten. Speedup, then sweeps, then the largest
parameter change against the fit at the default:

    1e-9    1.03  0.96  3.7e-7      1e-8    1.11  0.88  2.5e-6
    1e-7    1.17  0.82  6.0e-5

About a ninth of a fit for a sixth decimal place. Not taken: the estimates move,
and a looser floor puts a discontinuity into the objective wherever a unit flips
between one step and two, which the outer optimiser feels as a different path --
at 1e-9 that path effect is larger than the saving and one shape came out
slower. Recorded here rather than made a setting, so that the next person to
look at the mode solve has the number.
"""
@inline _laplace_inner_tolerance(laplace::CTSEMLaplaceObjective, value::Real) =
    max(laplace.inner_tol, 1e-10 * (one(value) + abs(value)))

"""
    _laplace_stationary_gain(gradient, step)

How much objective one Newton step was predicted to gain: `g' M^-1 g / 2`.

`step` is `M^-1 g` and is already computed by the solve, so this costs a dot
product. It is the quantity `_laplace_newton_unit_mode` needs when its line
search comes back empty -- in the same units as the value, which is the only
comparison that survives a reparameterisation of `u`. `Inf` when it is not
finite, so an unusable number can never certify a mode.
"""
@inline function _laplace_stationary_gain(gradient::AbstractVector,
        step::AbstractVector)
    isempty(gradient) && return 0.0
    gain = dot(gradient, step) / 2
    isfinite(gain) ? abs(gain) : Inf
end

"""
    _laplace_newton_unit_mode(laplace, U, values, Ls, aws, start)

Newton on `g_U` from one given starting point, reporting what happened rather
than writing anything back.

Split out of `_laplace_solve_unit_mode!` so that the same iteration can be run
twice from different starts -- see the cold retry there.
"""
function _laplace_newton_unit_mode(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws,
    start::Vector{Float64})
    d = laplace.units.dims[U]
    u = start
    repaired = false
    converged = false
    iterations = 0
    current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
        )
    if !isfinite(current.value) && any(!iszero, u)
        fill!(u, 0.0)
        current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
            )
    end
    for iteration in 1:laplace.inner_maxiter
        iterations = iteration
        if d == 0 || maximum(abs, current.gradient) <
                _laplace_inner_tolerance(laplace, current.value)
            converged = true
            break
        end
        # The band, not just the slot. Each Newton step builds the curvature
        # over every block, and that is the bulk of the primal -- passing only
        # the slot alone left the count at its default, so the whole mode solve
        # ran serially while everything else divided. It cost 4.6x down to
        # 1.2x on the value-only path and no test could have reported it.
        M = _laplace_unit_curvature(laplace, U, values, Ls, u)
        fac = _laplace_factor_repaired!(M, laplace.units.blocks[U])
        repaired |= fac.repaired
        fac.ok || break
        step = _laplace_block_solve(fac.factors, fac.coupling,
            laplace.units.blocks[U], current.gradient)
        accepted = false
        # Whether any trial could be evaluated at all, which is what separates
        # "there is nothing left to gain here" from "the step went somewhere
        # the model has no value". See the `stationary` branch below.
        evaluable = false
        scale = 1.0
        for _ in 1:20
            candidate = u .+ scale .* step
            trial = _laplace_unit_objective_gradient(laplace, U, values, Ls,
                candidate, aws)
            evaluable |= isfinite(trial.value)
            if isfinite(trial.value) && trial.value >= current.value - 1e-12
                u = candidate
                current = trial
                accepted = true
                break
            end
            scale /= 2
        end
        if !accepted
            # Aim for the tolerance; do not fail for missing it when the
            # arithmetic is what stopped you.
            #
            # The curvature here is repaired to negative definite before the
            # solve, so the Newton direction is an ascent direction. If no
            # scale of it -- down to 2^-20 -- produces a value that is even
            # equal, then the improvement still available has fallen below the
            # objective's own roundoff, and that is what being at a mode looks
            # like in floating point. Measured on a count model: |g| 2.43e-4
            # against a curvature of 7158, so `g^2/2M` = 4.1e-12 of value left
            # to gain, while the first trial came back 4.1e-11 lower -- noise,
            # on a value of size 242.
            #
            # Reporting that as a failure is what killed the fit rather than
            # the unit: one unit of forty invalidates the whole trial point, so
            # the outer optimiser was offered nothing usable and stopped at its
            # starting values, six times in twenty random starts.
            #
            # What separates the two is not whether a trial could be
            # evaluated -- it is how much objective the step that failed was
            # predicted to gain. `step` is `M^-1 g`, so `g'step/2` is exactly
            # that, for a dot product, and it is in the units the value is in.
            # Below the objective's own tolerance the line search is measuring
            # roundoff, which is the case above; at `g'step/2` of 1e4 it is not,
            # and calling that a mode hands the outer optimiser a Laplace term
            # evaluated away from one. That is not a function of theta: the same
            # parameter vector then evaluates differently depending on which
            # trial points preceded it, the line search follows a surface that
            # is not the objective, and the fit walks somewhere worse than it
            # started. Measured on a 50-subject nonlinear-drift model: units
            # flagged converged at `|dg/du|` of 4e7, an outer fit that ran from
            # -4282 to -14192, and three correction restarts that could not
            # recover it.
            #
            # `inner_gradient` keeps what was actually reached, so a unit that
            # stopped four orders above the tolerance is visible in the
            # diagnostics rather than silently equated with one that did not.
            converged = evaluable &&
                _laplace_stationary_gain(current.gradient, step) <=
                    _laplace_inner_tolerance(laplace, current.value)
            break
        end
    end
    if d == 0 || maximum(abs, current.gradient) <
            _laplace_inner_tolerance(laplace, current.value)
        converged = true
    end
    return (u=u, value=current.value, gradient=current.gradient,
        converged=converged, iterations=iterations, repaired=repaired)
end

"""
    _laplace_dual_unit_mode(laplace, U, values, Ls, uhat, curvature, aws)

The unit's inner mode as a function of the outer parameters, to first order.

One Newton step from the converged primal mode, taken with dual parameters.
The inner gradient's *primal* part is zero there, so the step's primal part is
zero and its dual part is exactly `-H^-1 dg/dtheta`: the implicit function
theorem, without forming that cross-derivative explicitly. Using the primal
`curvature` for the solve rather than a dual one is not an approximation for
the same reason -- any dual part of the inverse would multiply a zero primal
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
Whether the approximating Gaussian may claim more volume than the prior.

On by default. `ctsem_set_prior_floor!(false)` restores the unfloored Laplace,
which is what the measurements comparing the two are made against.
"""
const _LAPLACE_PRIOR_FLOOR = Ref(true)

"""
    ctsem_set_prior_floor!(on)

Turn the prior floor on `_laplace_prior_floor_logdet` describes on or off.
Returns the previous setting.
"""
function ctsem_set_prior_floor!(on::Bool)
    previous = _LAPLACE_PRIOR_FLOOR[]
    _LAPLACE_PRIOR_FLOOR[] = on
    return previous
end
export ctsem_set_prior_floor!

"""
    _laplace_prior_floor_logdet(logdetM)

`logdet(M)` floored at zero: the approximating Gaussian may not claim more
volume than the prior it is approximating.

# Why the unfloored term is unbounded

In the standardised metric `u ~ N(0, I)` and `g_U(u) = log L_U(u) - u'u/2`, so

    integral exp(g_U) du = (2pi)^(d/2) E_{u~N(0,I)}[L_U(u)] <= (2pi)^(d/2) sup_u L_U(u)

For fixed data `L_U` is bounded in `u`, so the *exact* term is finite always.
The Laplace term is not: it is `g_U(uhat) - logdet(M)/2` with
`M = -d2 log L/du du + I`, and `-logdet(M)/2 -> +inf` as an eigenvalue of `M`
goes to zero. So the surrogate diverges where the truth does not, and an
optimiser maximising it is not solving a well posed problem -- measured on a
25-subject model with a random `-log1p_exp(-param)` drift: one unit's smallest
eigenvalue at 2.6e-6 against 664 for its largest, a term 3.1 nats above every
neighbour within 1e-7 of theta, an outer gradient of 1.7e10 that is the true
derivative of a diverging quantity, and an estimate 26.7 nats *above* the
honest optimum. See `CT-SEM/review/LAPLACE-singular-unit-curvature-2026-09-22.md`.

# Why zero is the floor, and not a tuned constant

The mass the approximation claims is `L(uhat) exp(-uhat'uhat/2) (2pi)^(d/2) /
sqrt(det M)`, and the bound above caps the truth at `sup L (2pi)^(d/2)`. Since
`L(uhat) <= sup L` and `exp(-uhat'uhat/2) <= 1`, requiring `det M >= 1` is
sufficient for the approximation to respect that cap. `det M >= 1` is
`logdet M >= 0`, which is this.

The floor therefore has a meaning rather than a calibration. `M >= I` is
exactly "the likelihood is concave in `u`", because the `I` is the prior's own
curvature, and it implies `det M >= 1` -- so where the likelihood is concave,
every well behaved unit, this floor is inactive and nothing changes.

The converse does not hold. `logdet M >= 0` is a statement about the product of
the eigenvalues, so with `d > 1` one eigenvalue can go to zero while the others
keep the total positive: the floor then does not engage, and that direction
still contributes `-log(lambda_small)/2`. What this bounds is the divergence --
the term can no longer exceed `g_U(uhat)` -- not the credit a unit collects by
letting one direction go convex. And an optimiser can climb to the edge of the
floor: on the same model with six observations per subject, converged fits
leave units at eigenvalues like (8.5e-4, 14.9, 79.0), `logdet = -2.5e-10`, each
2.2 to 3.1 nats above its exact integral. `ctsem_laplace_conditioning` reports
such units, and the `:gated` floor (`ctsem_set_laplace_floor!`) scores them by
quadrature along the soft direction instead. See
`CT-SEM/review/LAPLACE-eigenwise-floor-2026-09-23.md`.

Measured on the model above: at the honest optimum the smallest eigenvalue over
all 25 units is 1.122 and no unit is floored; at the spurious one seven units
are below 0.5 and one is at 2.6e-6.

# What it does not claim

Not that the floored term is closer to the truth, and not that it is a bound on
it -- only that it obeys the same upper bound the exact integral does, so it
cannot diverge, and that where it binds the reported term is a bound rather
than the Laplace value. `ctLaplaceCheck()` measures the remaining gap against
quadrature, which is the question "how good is this approximation" and is not
this function's to answer.

`max` rather than a smoothed one: the term is `C0` at the crossing, with a kink
where a line search is no worse off than at any other kink, against a jump of
several nats where it is not floored. A smooth blend would need a width, which
is the calibration this avoids.
"""
@inline function _laplace_prior_floor_logdet(logdetM::T) where {T}
    _LAPLACE_PRIOR_FLOOR[] || return logdetM
    return max(logdetM, zero(T))
end

"""
Largest unit dimension decomposed densely: by the gated floor to find its soft
direction, and by `ctsem_laplace_conditioning` to report its smallest
eigenvalue. A wider unit keeps the total floor under `:gated`, and is reported
as `NaN` by the conditioning, because a dense eigendecomposition of a
block-sparse unit curvature is `O(d^3)` and fills in what the block
elimination keeps sparse.
"""
const _LAPLACE_EIGEN_MAXDIM = Ref(64)

"""
    ctsem_set_laplace_floor!(laplace, floor; lo, hi)

The prior floor this objective evaluates under. `:total`, the default, floors
each unit's `logdet(M)` at zero (`_laplace_prior_floor_logdet`). `:gated`
scores a unit whose inner curvature has an eigenvalue below `hi` by the gated
soft-direction rule, handing off to `:total` between `lo` and `hi` (defaults
0.2 and 0.7; see `_laplace_gated_term`).

The floor belongs to the objective, not to the session, so fitting one model
under `:gated` changes nothing about how any other objective evaluates. From R
it is set when the objective is built, from the `laplace_floor` entry of
`optimcontrol`, and
travels on the fit's specification, so every post-fit route that rebuilds the
objective evaluates under the floor the fit used. Returns the previous setting
as `(floor, lo, hi)`.
"""
function ctsem_set_laplace_floor!(laplace::CTSEMLaplaceObjective, floor;
    lo::Real=laplace.gate_lo, hi::Real=laplace.gate_hi)
    f = Symbol(floor)
    f in (:total, :gated) ||
        throw(ArgumentError("the Laplace floor must be :total or :gated"))
    0 < lo < hi || throw(ArgumentError("the gate needs 0 < lo < hi"))
    previous = (floor=laplace.floor, lo=laplace.gate_lo, hi=laplace.gate_hi)
    laplace.floor = f
    laplace.gate_lo = Float64(lo)
    laplace.gate_hi = Float64(hi)
    return previous
end
export ctsem_set_laplace_floor!

@inline _laplace_deepvalue(x::Real) = x
@inline _laplace_deepvalue(x::ForwardDiff.Dual) =
    _laplace_deepvalue(ForwardDiff.value(x))

"""
    _laplace_exceeds_identity(M, blocks, shift = 1)

Whether `M - shift I` is positive definite, i.e. every eigenvalue of `M`
exceeds `shift`; at one, the log likelihood is strictly concave in `u` at the
mode. One more block elimination with `M`'s own sparsity, so it costs what the
factorization already did and fills in nothing. A `false` is conservative --
the factorization also refuses a nearly singular `M - shift I` -- and the
callers (the gated floor's gate, `ctsem_laplace_conditioning`) then decompose
the unit, which gives the exact answer there.
"""
function _laplace_exceeds_identity(M::CTSEMBlockMatrix{T},
    blocks::Vector{CTSEMLaplaceBlock}, shift::Real=1.0) where {T}
    shifted = CTSEMBlockMatrix{T}([copy(d) for d in M.diag],
        [[copy(c) for c in row] for row in M.coupling])
    for d in shifted.diag
        for i in axes(d, 1); d[i, i] -= shift; end
    end
    ok, _, _, _ = _laplace_block_factor(shifted, blocks)
    return ok
end

################################################################################
# The gated soft-direction rule, as a floor mode with an exact gradient
################################################################################
#
# `floor = :gated` on the objective (`ctsem_set_laplace_floor!`, or from R
# `optimcontrol$laplace_floor = 'gated'`). A unit whose inner curvature `M` has
# every eigenvalue above `hi` keeps the plain Laplace term, unfloored, and pays
# one extra block elimination (`M - hi I`) to be told so -- its eigenvalues
# bound the term, and flooring it put a corner at logdet 0 that fits climbed to
# and stalled on (see `_laplace_gated_term`). A unit below it gets
#
#     T = T_lap + w(lambda_1) (T_soft - T_lap)
#
# with `w` a C1 smoothstep from one below `lo` to zero at `hi`, and `T_soft` a
# 3-point Gauss-Hermite rule along the smallest eigenvector `v_1` at the
# prior's scale, the other coordinates moved at each node by one Newton step
# with that node's own conditional curvature (clipped at one) and integrated by
# Laplace there. Measured against an exact reference on 16 weak-data configs:
# mean |error| 0.03 to 0.06 per eigenvalue class, where the total floor is off
# by +2.7 on near-singular units; see
# `CT-SEM/review/LAPLACE-eigenwise-floor-2026-09-23.md`, third to fifth addenda.
#
# The gradient of a flagged unit is ForwardDiff over that unit's term alone,
# with the inner mode's first derivative from `_laplace_dual_unit_mode` (the
# implicit function theorem, exactly as the nested oracle has it) and the
# eigen-quantities differentiated to first order from the primal decomposition:
# eigenvalues as `v' dM v`, the soft eigenvector by the standard perturbation
# sum, and the clipped conditional curvatures by the Daleckii-Krein divided
# differences, which stay finite across a near-degenerate pair on the same side
# of the clip. The cost is confined to flagged units: every other unit goes
# through the seeded assembly unchanged, and above the gate the gradient is the
# unfloored Laplace one, which is `:total`'s wherever `logdet(M) >= 0`.
#
# Scope. Dense, so units wider than `_LAPLACE_EIGEN_MAXDIM` keep `T_total`.
# Multilevel units no wider than
# that are handled whole. A unit whose two smallest eigenvalues are within a
# relative 1e-3 of each other also keeps `T_total`: the soft direction is not
# defined there, and a rule that chose one would be discontinuous in theta. On
# the weak-data study the smallest gap was 0.24.


"""
    _laplace_eig_first_order(A)

Eigenvalues and eigenvectors of a symmetric matrix, exact to first order in
the dual parts of `A`: the decomposition of the primal part, then
`dlambda_i = v_i' dA v_i` and `dv_i = sum_{k != i} v_k (v_k' dA v_i) /
(lambda_i - lambda_k)`. For a `Float64` matrix this is `_ctsem_symeig`.
"""
function _laplace_eig_first_order(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    A0 = [Float64(_laplace_deepvalue(A[i, j])) for i in 1:n, j in 1:n]
    E = _ctsem_symeig(A0)
    T === Float64 && return (values=E.values, vectors=E.vectors, primal=E)
    D = A .- A0
    V = E.vectors
    Dt = transpose(V) * D * V
    values = [E.values[i] + Dt[i, i] for i in 1:n]
    vectors = Matrix{T}(undef, n, n)
    for i in 1:n
        col = convert(Vector{T}, V[:, i])
        for k in 1:n
            k == i && continue
            col = col .+ V[:, k] .* (Dt[k, i] / (E.values[i] - E.values[k]))
        end
        vectors[:, i] = col
    end
    return (values=values, vectors=vectors, primal=E)
end

"""
    _laplace_clip_first_order(B)

`V max(Lambda, 1) V'` for symmetric `B`, exact to first order in its dual
parts by Daleckii-Krein: in the primal eigenbasis the derivative's `(i, j)`
entry is `(f(l_i) - f(l_j)) / (l_i - l_j)` times `(V' dB V)_ij`, with `f'` on
the diagonal and wherever the pair coincides. No division by a small gap.
"""
function _laplace_clip_first_order(B::AbstractMatrix{T}) where {T}
    n = size(B, 1)
    B0 = [Float64(_laplace_deepvalue(B[i, j])) for i in 1:n, j in 1:n]
    E = _ctsem_symeig(B0)
    V = E.vectors
    l = E.values
    f(x) = max(x, 1.0)
    fp(x) = x > 1.0 ? 1.0 : 0.0
    Dt = transpose(V) * (B .- B0) * V
    C = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        coef = (i == j || abs(l[i] - l[j]) <= 1e-12 * max(1.0, abs(l[i]))) ?
            fp(l[i]) : (f(l[i]) - f(l[j])) / (l[i] - l[j])
        C[i, j] = coef * Dt[i, j] + (i == j ? f(l[i]) : 0.0)
    end
    return V * C * transpose(V)
end

@inline function _laplace_smoothstep_weight(lambda, lo::Real, hi::Real)
    lp = _laplace_deepvalue(lambda)
    lp <= lo && return one(lambda)
    lp >= hi && return zero(lambda)
    s = (lambda - lo) / (hi - lo)
    return 1 - (3 * s^2 - 2 * s^3)
end

"""
    _laplace_gated_term(laplace, U, values, Ls, u, aws; M, logdetM, inner)

Unit `U`'s term under the gated rule, generic in the element type: returns
`(value, flagged, floor)`. `M`, `logdetM` and `inner` may be passed when the
caller already has them at `u`.

`floor` says whether the unit still needs the total floor. Only the explicit
fallbacks do -- a unit wider than `_LAPLACE_EIGEN_MAXDIM` or with a degenerate
soft pair, where no eigenvalue is known. Everywhere else the eigenvalue is, and
it makes the floor unnecessary: above the band (`lambda_min >= hi`) the plain
Laplace term is bounded by `-d log(hi) / 2` and accurate there (+0.015 nats on
units with eigenvalues in 0.7 to 1, against the exact reference), and inside
the band it is blended with the soft rule, which is bounded by construction.

Flooring those units anyway was a kink. A unit with every eigenvalue in (hi,
1) has `logdet(M) < 0`, and `max(logdet, 0)` turns its term into a corner the
optimiser climbs to and then cannot leave: measured on AnomAuth (800 subjects)
from its spurious Laplace maximum, 69 to 77 units sat on it at once and L-BFGS
stopped with |g| = 144 and no ascent step, the continuation then spending an
hour at zero gain.
"""
function _laplace_gated_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws; M=nothing, logdetM=nothing, inner=nothing) where {T}
    blocks = laplace.units.blocks[U]
    inner === nothing &&
        (inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws))
    g = inner.value
    (isfinite(g) && !isempty(u)) || return (value=g, flagged=false, floor=true)
    M === nothing && (M = _laplace_unit_curvature(laplace, U, values, Ls, u))
    # A curvature that does not factor at the mode -- the inner solve stopped
    # at a point that is not a strict maximum, the primal repaired it, and the
    # repaired Laplace term is unbounded there -- is scored by the soft rule
    # alone (its smallest eigenvalue is below `lo`, so the weight is one). The
    # early returns below then report `flagged = false` with a NaN, and the
    # caller keeps its own term, as it did before this branch existed.
    indefinite = false
    if logdetM === nothing
        ok, logdetM, _, _ = _laplace_block_factor(M, blocks)
        ok || (indefinite = true)
    end
    # `total` only for the fallbacks, which have no eigenvalue to go on;
    # `lap`, unfloored, wherever there is one. See the docstring.
    total = indefinite ? T(NaN) : g - max(logdetM, zero(logdetM)) / 2
    lap = indefinite ? T(NaN) : g - logdetM / 2
    lo, hi = laplace.gate_lo, laplace.gate_hi
    _laplace_exceeds_identity(M, blocks, hi) &&
        return (value=lap, flagged=false, floor=false)
    d = length(u)
    d > _LAPLACE_EIGEN_MAXDIM[] && return (value=total, flagged=false, floor=true)
    dense = _laplace_block_dense(M, blocks, d)
    E = _laplace_eig_first_order(dense)
    lam = E.primal.values
    if d > 1 && (lam[2] - lam[1]) <= 1e-3 * max(abs(lam[2]), 1.0)
        return (value=total, flagged=false, floor=true)
    end
    lam1 = E.values[1]
    w = _laplace_smoothstep_weight(lam1, lo, hi)
    iszero(_laplace_deepvalue(w)) && return (value=lap, flagged=true, floor=false)
    v1 = E.vectors[:, 1]
    # The other directions as an orthonormal basis of v1's complement, by
    # Gram-Schmidt from the primal eigenvectors. Everything below depends on
    # that subspace only, not on the basis, so a degenerate pair among the
    # stiff directions costs nothing here.
    Vh = Matrix{T}(undef, d, d - 1)
    for k in 2:d
        c = convert(Vector{T}, E.primal.vectors[:, k])
        c = c .- v1 .* sum(v1 .* c)
        for j in 1:(k - 2)
            c = c .- Vh[:, j] .* sum(Vh[:, j] .* c)
        end
        Vh[:, k - 1] = c ./ sqrt(sum(c .* c))
    end
    lam1c = max(lam1, one(lam1))
    xs, ws = _gauss_hermite(3)
    terms = Vector{T}(undef, length(xs))
    for j in eachindex(xs)
        x = xs[j]
        local h
        if abs(x) < 1e-12 || d == 1
            if d == 1
                ub = u .+ (sqrt(2) * x / sqrt(lam1c)) .* v1
                h = abs(x) < 1e-12 ? g :
                    _laplace_unit_objective_gradient(laplace, U, values, Ls, ub, aws).value
            else
                Bc = _laplace_clip_first_order(transpose(Vh) * dense * Vh)
                F = _ctsem_cholesky(Matrix(_laplace_symmetrise(Bc)), d - 1)
                h = g - logdet(F) / 2
            end
        else
            ub = u .+ (sqrt(2) * x / sqrt(lam1c)) .* v1
            r = _laplace_unit_objective_gradient(laplace, U, values, Ls, ub, aws)
            Mb = _laplace_block_dense(_laplace_unit_curvature(laplace, U, values,
                Ls, ub), blocks, d)
            Bc = _laplace_clip_first_order(transpose(Vh) * Mb * Vh)
            F = _ctsem_cholesky(Matrix(_laplace_symmetrise(Bc)), d - 1)
            issuccess(F) || return (value=T(NaN), flagged=true, floor=false)
            z = F \ (transpose(Vh) * r.gradient)
            h = _laplace_unit_objective_gradient(laplace, U, values, Ls,
                ub .+ Vh * z, aws).value - logdet(F) / 2
        end
        terms[j] = h + log(ws[j]) + x^2
    end
    peak = maximum(t -> Float64(_laplace_deepvalue(t)), terms)
    isfinite(peak) || return (value=T(NaN), flagged=true, floor=false)
    soft = peak + log(sum(exp.(terms .- peak))) - log(lam1c) / 2 - log(pi) / 2
    indefinite && return (value=soft, flagged=true, floor=false)
    return (value=lap + w * (soft - lap), flagged=true, floor=false)
end

"""
    _laplace_gated_unit_gradient(laplace, U, theta, curvature)

The exact gradient of a flagged unit's gated term: ForwardDiff over the term,
the inner mode's first derivative from `_laplace_dual_unit_mode`.
"""
function _laplace_gated_unit_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, curvature)
    uhat = laplace.modes[U]
    # A repaired unit's factors are of a shifted curvature, which is not the
    # implicit function theorem's matrix; take the mode's derivative from the
    # unshifted one's eigendecomposition instead (it is nonsingular wherever
    # the mode moves smoothly, whatever its signature).
    dense_inverse = nothing
    if laplace.mode_repaired[U]
        Ls0 = _laplace_popchols(theta, laplace.spec)
        E = _ctsem_symeig(_laplace_block_dense(_laplace_unit_curvature(laplace, U,
            theta, Ls0, uhat), laplace.units.blocks[U], length(uhat)))
        dense_inverse = E.vectors * Diagonal(1 ./ E.values) * transpose(E.vectors)
    end
    term_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        ud = if dense_inverse === nothing
            _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, curvature, wsd)
        else
            u0 = convert(Vector{S}, uhat)
            u0 .+ dense_inverse * _laplace_unit_objective_gradient(laplace, U, x,
                Lsd, u0, wsd).gradient
        end
        return _laplace_gated_term(laplace, U, x, Lsd, ud, wsd).value
    end
    return ForwardDiff.gradient(term_of, theta)
end

"""
    _laplace_unit_term(laplace, U, values, Ls, u, aws)

`g_U(u) - logdet(-d2 g_U/du du) / 2`: unit `U`'s contribution to the
approximated log marginal likelihood.

The `(2*pi)^(d/2)` the Laplace approximation produces cancels the
`(2*pi)^(-d/2)` in the standard-normal density of `u`, so no dimension-dependent
constant appears here whatever the unit's size.
"""
function _laplace_unit_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws) where {T}
    inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
        )
    isfinite(inner.value) || return inner.value
    isempty(u) && return inner.value
    M = _laplace_unit_curvature(laplace, U, values, Ls, u)
    ok, logdetM, _, _ = _laplace_block_factor(M, laplace.units.blocks[U])
    ok || return T(NaN)
    if _LAPLACE_PRIOR_FLOOR[] && laplace.floor === :gated
        return _laplace_gated_term(laplace, U, values, Ls, u, aws; M=M,
            logdetM=logdetM, inner=inner).value
    end
    return inner.value - _laplace_prior_floor_logdet(logdetM) / 2
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
    _laplace_run_members(body, laplace, members)

Run `body(c, m)` for every member, over whatever the calling task's pool band
holds. Serial when that band is one slot wide or there is little to divide, in
which case every call uses the caller's own slot and nothing is spawned.
"""
function _laplace_run_members(body, laplace::CTSEMLaplaceObjective, members)
    return _laplace_parallel(laplace, 1:length(members)) do c
        body(c, @inbounds members[c])
    end
end

"""
    _laplace_unit_seeded_gradient(laplace, U, values, Ls, u, members, d1, d2,
                                  order)

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
    members, d1::Vector{Float64}, d2::Vector{Float64}, order::Integer)
    # The switch folds to a constant, so with counting off this is a plain call
    # and the try/finally below is not compiled at all. It is a sweep: the
    # assembly runs thousands per gradient.
    _LAPLACE_COUNT_BYTES || return _laplace_unit_seeded_gradient_(laplace, U,
        values, Ls, u, members, d1, d2, order)
    bs = _laplace_mark()
    try
        return _laplace_unit_seeded_gradient_(laplace, U, values, Ls, u,
            members, d1, d2, order)
    finally
        _laplace_charge!(_LAPLACE_BYTES_SWEEP, bs)
    end
end

function _laplace_unit_seeded_gradient_(laplace::CTSEMLaplaceObjective, U::Integer,
    values::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    members, d1::Vector{Float64}, d2::Vector{Float64}, order::Integer)

    spec = laplace.spec
    units = laplace.units
    unitmembers = units.members[U]
    npar = length(values)
    nm = length(members)
    failure = (ok=false, d0=_LAPLACE_NO_MATRIX, d1c=_LAPLACE_NO_MATRIX,
        d12=_LAPLACE_NO_MATRIX)

    # Not zeroed, and not freshly allocated. Every column is written by the
    # member that owns it before anything reads one, and a member that fails
    # abandons the whole sweep -- so there is no path that reads an entry this
    # call did not write. The one caller that *keeps* a result past the next
    # sweep copies it; see `levelsweeps` below, which is the reason this is
    # spelled out rather than left to be noticed.
    d0 = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d0)
    d1c = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d1c)

    if order == 1
        seed = ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0)
        S = typeof(seed)
        one_member! = function (c::Int, m)
            local aws, gradient, x, shifted, loglik, t
            aws = _laplace_workspace!(laplace, S, npar)
            gradient = _laplace_scratch_vector!(laplace, S, npar, :sweep_grad)
            x = _laplace_scratch_vector!(laplace, S, npar, :sweep_x)
            shifted = _laplace_member_values!(
                _laplace_scratch_vector!(laplace, Float64, npar, :sweep_shift),
                values, spec, Ls, u, units.offsets[U][m])
            @inbounds for t in 1:npar
                x[t] = shifted[t] + seed * d1[t]
            end
            loglik = _laplace_subject_value_gradient!(gradient,
                laplace.objective.subject_objectives[unitmembers[m]], aws, x)
            _laplace_finite(loglik) || return false
            @inbounds for t in 1:npar
                d0[t, c] = ForwardDiff.value(gradient[t])
                d1c[t, c] = ForwardDiff.partials(gradient[t])[1]
            end
            return true
        end
        _laplace_run_members(one_member!, laplace, members) || return failure
        return (ok=true, d0=d0, d1c=d1c, d12=_LAPLACE_NO_MATRIX)
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
    d12 = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d12)
    # One member. Writes only column `c` of the three outputs, so members never
    # collide however they are divided; the workspace comes from the ambient
    # slot, and the dual buffers are cached against it.
    one_member! = function (c::Int, m)
        local aws, gradient, x, shifted, loglik, t, inner
        aws = _laplace_workspace!(laplace, S, npar)
        gradient = _laplace_scratch_vector!(laplace, S, npar, :sweep_grad)
        x = _laplace_scratch_vector!(laplace, S, npar, :sweep_x)
        shifted = _laplace_member_values!(
            _laplace_scratch_vector!(laplace, Float64, npar, :sweep_shift),
            values, spec, Ls, u, units.offsets[U][m])
        @inbounds for t in 1:npar
            x[t] = shifted[t] + e1 * d1[t] + e2 * d2[t]
        end
        loglik = _laplace_subject_value_gradient!(gradient,
            laplace.objective.subject_objectives[unitmembers[m]], aws, x)
        _laplace_finite(loglik) || return false
        @inbounds for t in 1:npar
            inner = ForwardDiff.value(gradient[t])
            d0[t, c] = ForwardDiff.value(inner)
            d1c[t, c] = ForwardDiff.partials(inner)[1]
            d12[t, c] = ForwardDiff.partials(ForwardDiff.partials(gradient[t])[1])[1]
        end
        return true
    end
    _laplace_run_members(one_member!, laplace, members) || return failure
    return (ok=true, d0=d0, d1c=d1c, d12=d12)
end

# Deep: `isfinite` on a dual tests only its value, and a sweep's whole point
# is the partials. Testing the value alone is what let a NaN derivative travel
# from a unit whose predicted variance had collapsed all the way into the
# assembled gradient, where it was indistinguishable from a bad trial point.
@inline _laplace_finite(x::Real) = isfinite(x)
@inline _laplace_finite(x::ForwardDiff.Dual) = _finite_deep(x)

"""
    _laplace_seeded_unit_gradient!(out, laplace, U, values, Ls, dL, M, factors,
                                   elim)

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
    M::CTSEMBlockMatrix{Float64}, factors, elim)

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
            1:nmem, zerodir, zerodir, 1)
        pass.ok || return false
        @inbounds for m in 1:nmem, t in 1:npar
            out[t] += pass.d0[t, m]
        end
        return true
    end

    bsi = _laplace_mark()
    Cdiag, Ccoup = _laplace_selected_inverse(factors, elim, blocks)
    _laplace_charge!(_LAPLACE_BYTES_SELINV, bsi)
    # The selected inverse of the curvature is where a badly conditioned unit
    # first shows, and every sweep below is scaled by it.
    if !(all(x -> all(isfinite, x), Cdiag) &&
            all(r -> all(x -> all(isfinite, x), r), Ccoup))
        return false
    end
    # A sweep this function runs itself, on its own slot and serially: the
    # parallelism here is one level out, over members.
    sweep = (mm, a1, a2, order) -> _laplace_unit_seeded_gradient(laplace, U,
        values, Ls, uhat, mm, a1, a2, order)
    sweep_at = (mm, a1, a2, order) -> _laplace_unit_seeded_gradient(
        laplace, U, values, Ls, uhat, mm, a1, a2, order)

    # Members, partitioned across the pool: `mine` says which of the unit's
    # members a worker owns.
    #
    # Partitioning *members* rather than blocks is what makes this safe without
    # a single barrier. Every accumulator here -- `Pm`, `llvm`, `Bsm`, `seen` --
    # is indexed by member, so a worker writes only the columns it owns,
    # whatever block or level the contribution came from. An earlier attempt
    # partitioned blocks instead and had to argue that blocks of one level
    # cover disjoint members, which needs a barrier between levels.
    #
    # The cost of the choice is that every worker walks the whole block list
    # and skips blocks holding none of its members, and that a block spanning
    # many members has its small Cholesky factored once per worker. Both are
    # nothing beside one subject sweep.
    run_by_member! = function (body)
        return _laplace_partition(laplace, nmem) do mine, _w
            body(collect(mine))
        end
    end

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
    # `Vector{Bool}`, not `falses`. A `BitVector` packs 64 flags into a word, so
    # `seen[m] = true` is a read-modify-write of that word: two tasks setting
    # *different* members race on it and one write is simply lost.
    seen = fill(false, nmem)

    assemble! = function (mine::Vector{Int})
        # `local`, and not for tidiness. Assignment inside a closure binds the
        # *enclosing* function's variable whenever that name is already a local
        # there, and this function assigns `pass`, `dir`, `s`, `b`, `block` and
        # `l` in its own serial parts. Without these declarations all the tasks
        # share one `pass`: the size check above passes, and by the time the
        # loop reads it another task has replaced it with its own block's
        # result. It surfaces as a BoundsError indexing column 7 of a
        # one-column matrix, and would otherwise be a wrong gradient.
        local ownedA, b, block, l, L, k, subsA, FC, Q, j, dir, pass, c, m, t
        local kparl, la, nb, na, accumulate, Cab, Vcross, e, q, i, tt
        ownedA = fill(false, nmem)
        @inbounds for m in mine; ownedA[m] = true; end
        for (b, block) in enumerate(blocks)
            l = block.level
            L = Ls[l]
            k = block.size
            k == 0 && continue
            subsA = Int[m for m in block.members if ownedA[m]]
            isempty(subsA) && continue
        # Diagonal term. tr(C[b,b] A[b,b]) = tr(W H) with W = L C[b,b] L', so a
        # factor of W turns the trace into pure second directional derivatives
        # whose directions live in parameter space.
        #
        # The factor is taken of `C[b,b]`, not of `W`. `W` is parameters by
        # parameters with the rank of the block, so for a reduced level it is
        # singular and a Cholesky of it simply fails -- which returned `false`
        # here and sent every reduced fit to the nested gradient by the silent
        # fallback. `C[b,b]` is block by block and positive definite whatever
        # the rank, and `W = (L Cc)(L Cc)'` exactly, so `Q = L Cc` has one
        # column per block dimension, which is the number of directions the
        # trace needs. At full rank this is the same factorisation of the same
        # matrix, reached without forming it.
        FC = cholesky(Symmetric(_laplace_symmetrise(Cdiag[b])); check=false)
        issuccess(FC) || (_LAPLACE_DIAG[] = 11; return false)
        Q = L * Matrix(FC.L)
        # One sweep per direction, and it stays that way. This is the only site
        # whose sweep count grows with the number of random effects -- 100k
        # member-sweeps against 100 for every other site on a single-level
        # model -- so it is where the eye goes, and the obvious move is to seed
        # all `k` directions into one pass as `k` partials on each tag, reading
        # the diagonal of the block that comes back.
        #
        # Built and measured. It is correct -- the whole suite passes, oracles
        # included -- and it is slower. A width-`W` nested dual costs `(1 + W)^2`
        # components per scalar operation against `4W` for `W` separate passes,
        # and `(1 + W)^2 - 4W = (W - 1)^2`, so the wide pass always does more
        # arithmetic. Whether that matters depends on how much of a pass is
        # arithmetic rather than tape and control flow, and that depends on the
        # subject. Timing this engine's own filter at both widths, on dev1:
        #
        #                  1 latent, 5 rows      3 latents, 24 rows
        #     W = 2              1.63x                   1.18x
        #     W = 3              1.80                    1.01
        #     W = 4              1.87                    0.97
        #     W = 6              1.80                    0.79
        #
        # On a toy subject the fixed per-pass cost dominates and the wide pass
        # wins by most of a factor of two; on a subject the size anyone fits,
        # the arithmetic dominates and it loses from three directions up. End to
        # end on a 3-latent, 100-subject model the widened version ran 3.6%,
        # 2.5%, 13.8% and 17.7% slower at one, two, four and six random effects.
        #
        # The first probe of this used the one-latent fixture and reported 1.8x,
        # which is why the measurement above names its model. See
        # `review/LAPLACE-where-the-time-goes-2026-09-23.md`.
        for j in 1:k
            dir = scatter(l, Q[:, j])
            pass = sweep_at(subsA, dir, dir, 2)
            pass.ok || (_LAPLACE_DIAG[] = 12; return false)
            # The second-order sweep is where this fails when it fails: it is a
            # third derivative of the process model once the reverse pass is
            # counted, and a unit whose predicted variance has collapsed can
            # produce a NaN there with a perfectly finite log likelihood and
            # first derivative. Caught here so the caller can take the nested
            # route rather than carry the NaN into the sum.
            (all(isfinite, pass.d12) && all(isfinite, pass.d0)) ||
                (_LAPLACE_DIAG[] = 13; return false)
            size(pass.d12, 2) == length(subsA) || error(
                "diag sweep: block $b level $l j $j returned " *
                "$(size(pass.d12, 2)) columns for $(length(subsA)) members " *
                "(block has $(length(block.members)))")
            for (c, m) in enumerate(subsA)
                for t in 1:npar
                    Pm[t, m] += pass.d12[t, c]
                    seen[m] || (llvm[t, m] = pass.d0[t, c])
                end
                seen[m] = true
            end
        end
        # Cross terms with each ancestor, twice over as the trace requires.
        # `tr(C[a,b] A[b,a]) = tr(V H)` with `V = L_a C[a,b] L_b'`, and there
        # are two ways to spell that trace as directional derivatives.
        #
        # Over *parameters*: `V` absorbs both factors, so the first direction is
        # a bare basis vector and the second is `V`'s matching column. One sweep
        # per parameter of this level.
        #
        # Over *dimensions*: `tr(V H) = tr(C[a,b] L_b' H L_a)`, so sweeping
        # along the two loadings' own columns and weighting by `C[a,b]` gives
        # the same number in `rank_b * rank_a` sweeps. The two are equal by
        # linearity -- expand `V[:,q] = sum_j L_a[:,j] sum_i C[a,b][j,i] L_b[q,i]`
        # and collect over `q`.
        #
        # Neither dominates. A subject block of 18 parameters at rank 6 under a
        # rank-1 study is 6 sweeps against 18; a burst block of 6 parameters at
        # rank 3 under a rank-6 subject is 18 against 6. So the cheaper one is
        # chosen per pair, which at full rank is always the parameter form
        # because `rank == kpar` there and `rank^2 >= rank`.
        kparl = length(spec.levels[l].re_index)
        for (t, a) in enumerate(block.ancestors)
            la = blocks[a].level
            nb = block.size
            na = blocks[a].size
            accumulate = function (pass, subs, weight)
                pass.ok || (_LAPLACE_DIAG[] = 14; return false)
                all(isfinite, pass.d12) || (_LAPLACE_DIAG[] = 15; return false)
                size(pass.d12, 2) == length(subs) || error(
                    "seeded assembly: sweep returned $(size(pass.d12, 2)) " *
                    "columns for $(length(subs)) members")
                @inbounds for (c, m) in enumerate(subs)
                    for tt in 1:npar
                        Pm[tt, m] += weight * pass.d12[tt, c]
                    end
                end
                return true
            end
            if nb * na < kparl
                Cab = transpose(Ccoup[b][t])          # C[a,b], rank_a by rank_b
                for i in 1:nb, j in 1:na
                    accumulate(sweep_at(subsA, scatter(l, L[:, i]),
                        scatter(la, Ls[la][:, j]), 2), subsA,
                        2 * Cab[j, i]) || return false
                end
            else
                # `Vcross`, not `V`: the explicit terms below use a `V` of
                # their own at this function's scope, and this assignment would
                # land on it. Harmless only because that one is reallocated
                # before it is read -- a reordering, or lifting this loop into
                # a closure, makes it a silent wrong answer.
                Vcross = Ls[la] * transpose(Ccoup[b][t]) * transpose(L)
                for q in 1:kparl
                    e = zeros(Float64, kparl); e[q] = 1.0
                    accumulate(sweep_at(subsA, scatter(l, e),
                        scatter(la, Vcross[:, q]), 2), subsA, 2.0) || return false
                end
            end
        end
    end
    return true
    end
    run_by_member!(assemble!) || return false

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
        # `block.size` is the block's width in `u`; `length(rho)` is how many
        # parameters the level moves. Equal for a full-rank level and not
        # otherwise, so the accumulator is built in parameter space and `L'`
        # brings it back to block space.
        k = block.size
        kpar = length(rho)
        k == 0 && continue
        acc = zeros(Float64, kpar)
        for m in block.members, p in 1:kpar
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
    sB! = function (mine::Vector{Int})
        # See the note in `assemble!`: without this the tasks share `pass`.
        local ownedB, b, block, k, subsB, sb, dir, pass, c, m, t
        ownedB = fill(false, nmem)
        @inbounds for m in mine; ownedB[m] = true; end
        for (b, block) in enumerate(blocks)
            k = block.size
            k == 0 && continue
            subsB = Int[m for m in block.members if ownedB[m]]
            isempty(subsB) && continue
            sb = [s[block.offset + q] for q in 1:k]
            dir = scatter(block.level, Ls[block.level] * sb)
            pass = sweep_at(subsB, dir, dir, 1)
            pass.ok || (_LAPLACE_DIAG[] = 16; return false)
            all(isfinite, pass.d1c) || (_LAPLACE_DIAG[] = 17; return false)
            size(pass.d1c, 2) == length(subsB) || error(
                "seeded assembly: s'B sweep returned $(size(pass.d1c, 2)) " *
                "columns for $(length(subsB)) members")
            @inbounds for (c, m) in enumerate(subsB)
                for t in 1:npar; Bsm[t, m] += pass.d1c[t, c]; end
            end
        end
        return true
    end
    run_by_member!(sB!) || return false

    # dv/dtheta is the identity away from the population parameters, so for
    # every other parameter the contribution is a plain read-off.
    if !(all(isfinite, llvm) && all(isfinite, Pm) && all(isfinite, Bsm))
        _LAPLACE_DIAG[] = all(isfinite, llvm) ? (all(isfinite, Pm) ? 3 : 2) : 1
        return false
    end
    @inbounds for j in 1:npar
        acc = 0.0
        for m in 1:nmem
            acc += llvm[j, m] + (Pm[j, m] + Bsm[j, m]) / 2
        end
        out[j] += acc
    end

    # The explicit trace terms, written against `dL` rather than against
    # `X = L \ dL`.
    #
    # The `X` form needs no sweeps at all and is kept wherever it is valid, but
    # it exists only for a square `L`: it rewrites `tr(C L' H dL)` as
    # `tr(C L' H L X)`, which is the same thing only because `L X = dL` has a
    # solution. A reduced level's `L` is parameters by rank, so it has none.
    #
    # The other way round, each of the three traces is
    # `<H_{b,x} L_x C[x,b], dL_b>` for `x` over the block, its ancestors and
    # its descendants: one Frobenius product per population parameter against a
    # matrix that does not depend on which parameter it is. `H_{b,x} L_x` is a
    # directional second derivative along a *parameter*-space direction, which
    # is what one order-1 sweep along a column of `L_x` returns -- so this costs
    # one sweep per level dimension rather than one per level parameter.
    #
    # `H_{b,x}` runs over the members depending on both blocks, which under
    # strict nesting is the members of whichever is inner.
    nlev = length(spec.levels)
    needV = hasreducedrank(spec)
    levelsweeps = Vector{Vector{Matrix{Float64}}}(undef, nlev)
    if needV
        for l in 1:nlev
            rl = nlatent(spec.levels[l])
            levelsweeps[l] = Vector{Matrix{Float64}}(undef, rl)
            for q in 1:rl
                pass = sweep(1:nmem, scatter(l, Ls[l][:, q]), zerodir, 1)
                pass.ok || return false
                all(isfinite, pass.d1c) || return false
                # A copy, because the sweep's output buffer is reused by
                # the next `q` and these are read after the loop.
                levelsweeps[l][q] = copy(pass.d1c)
            end
        end
    end

    # `H_{target,partner} L_partner`, in the target level's parameter rows,
    # summed over the given members.
    HLmatrix = function (targetlevel::Int, partnerlevel::Int, memberset)
        rho = spec.levels[targetlevel].re_index
        rp = nlatent(spec.levels[partnerlevel])
        acc = zeros(Float64, length(rho), rp)
        for q in 1:rp
            D = levelsweeps[partnerlevel][q]
            @inbounds for m in memberset, pp in eachindex(rho)
                acc[pp, q] += D[rho[pp], m]
            end
        end
        return acc
    end

    V = Vector{Matrix{Float64}}(undef, length(blocks))
    if needV
        for (b, block) in enumerate(blocks)
            V[b] = zeros(Float64,
                length(spec.levels[block.level].re_index), block.size)
        end
        for (b, block) in enumerate(blocks)
            block.size == 0 && continue
            V[b] .+= HLmatrix(block.level, block.level, block.members) * Cdiag[b]
            # A coupling moves both ends: this block's factor against the
            # ancestor, and the ancestor's factor against this block. Both run
            # over this block's members, it being the inner of the two.
            for (tt, a) in enumerate(block.ancestors)
                blocks[a].size == 0 && continue
                V[b] .+= HLmatrix(block.level, blocks[a].level, block.members) *
                    transpose(Ccoup[b][tt])
                V[a] .+= HLmatrix(blocks[a].level, block.level, block.members) *
                    Ccoup[b][tt]
            end
        end
    end

    # The population parameters move every member's v through L, and move psi
    # through L explicitly. A[b,b] = I - M.diag[b] and A[b,a] = -M.coupling[b][t]
    # give the curvature blocks with no further sweeps, and writing each trace
    # through A rather than through the v-space Hessian needs only one
    # triangular solve X = L \ dL.
    Gb = Vector{Vector{Float64}}(undef, length(blocks))
    Fb = Vector{Vector{Float64}}(undef, length(blocks))
    for (b, block) in enumerate(blocks)
        rho = spec.levels[block.level].re_index
        kpar = length(rho)
        Gb[b] = zeros(Float64, kpar)
        Fb[b] = zeros(Float64, kpar)
        for m in block.members, p in 1:kpar
            Gb[b][p] += llvm[rho[p], m] + (Pm[rho[p], m] + Bsm[rho[p], m]) / 2
            Fb[b][p] += llvm[rho[p], m]
        end
    end

    for l in eachindex(spec.levels)
        isempty(dL[l]) && continue
        levelpositions = _laplace_level_positions(spec, l)
        for (t, j) in enumerate(levelpositions)
            X = needV ? zeros(Float64, 0, 0) : Ls[l] \ dL[l][t]
            total = 0.0
            for (b, block) in enumerate(blocks)
                k = block.size
                if block.level == l && k > 0
                    kpar = length(spec.levels[l].re_index)
                    ub = [uhat[block.offset + q] for q in 1:k]
                    sb = [s[block.offset + q] for q in 1:k]
                    # `dL` is parameter-by-block, so both shifts land in
                    # parameter space where `Gb` and `Fb` live.
                    shift = dL[l][t] * ub
                    sshift = dL[l][t] * sb
                    for p in 1:kpar
                        total += Gb[b][p] * shift[p] + Fb[b][p] * sshift[p] / 2
                    end
                    if needV
                        for x in axes(V[b], 1), y in axes(V[b], 2)
                            total += V[b][x, y] * dL[l][t][x, y]
                        end
                    else
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
                end
                # tr(C[a,b] A[b,a] X): an ancestor's factor moving.
                if !needV
                for (tt, a) in enumerate(block.ancestors)
                    blocks[a].level == l || continue
                    ka = blocks[a].size
                    Z = transpose(Ccoup[b][tt]) * (.-M.coupling[b][tt])
                    for x in 1:ka, y in 1:ka
                        total += Z[x, y] * X[y, x]
                    end
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
    isreducedrank(spec.levels[l]) ? copy(spec.levels[l].load_index) :
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
        r = nlatent(level)
        if isempty(positions) || k == 0 || r == 0
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
        out[l] = [Matrix{Float64}(reshape(collect(view(J, :, t)), k, r))
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

`unit_loglik` is each unit's own term -- the approximated log marginal
likelihood of all its members jointly. `subject_loglik` spreads that evenly
over the members, which is exact only when a unit is one subject; anything that
needs a marginal per independent block (leave-one-unit-out, for one) reads
`unit_loglik`.
"""
function ctsem_laplace_evaluate(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, nested_gradient::Bool=false)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    laplace.last_values = theta
    nsubjects = length(laplace.objective.subject_objectives)
    nunits = length(laplace.units.members)

    # 1. Inner modes and the value at them, in primal arithmetic, each solved
    #    from the origin. Each subject's term is its own approximated log
    #    marginal likelihood, which is the per-subject quantity that means the
    #    same thing here as `subject_loglik` does without random effects.
    #
    #    The curvature is computed once and used three times -- for the
    #    definiteness check, for the log determinant, and as the primal solve in
    #    `_laplace_dual_mode` below. It is the most expensive primal quantity
    #    here, so recomputing it for each of those would be a third of the
    #    primal pass thrown away.
    Ls = _laplace_popchols(theta, laplace.spec)
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
    # The gated floor's per-unit record: which units its rule scored.
    gated_floor = _LAPLACE_PRIOR_FLOOR[] && laplace.floor === :gated
    primal_gated = fill(false, nunits)
    unit_loglik = zeros(Float64, nunits)
    subject_loglik = zeros(Float64, nsubjects)
    value = 0.0

    # The subject loop is the parallelism here, and it is the natural one: each
    # subject's mode solve, curvature and value are completely independent, and
    # nothing is shared but the read-only parameter vector. Chunk count comes
    # from `ctsem_set_max_chunks!`, which is what the R side sets from `cores`,
    # so it is the same control the non-Laplace path uses.
    # One piece per unit, and deliberately *not* one per worker.
    #
    # `_ctsem_nchunks` returns `min(max_chunks, nunits)`, which ties the number
    # of pieces to the size of the pool. That is the wrong coupling and it is
    # worst exactly where balance matters most: with thirteen studies and five
    # workers it makes five pieces, each holding two or three studies, so every
    # worker takes one piece and there is nothing left to steal. A worker that
    # draws three small studies finishes and idles while the one holding the
    # big study grinds, and no amount of dynamic pulling can help because the
    # queue is empty.
    #
    # Pieces are for balance and workers are for hardware, so they are sized
    # separately. One unit each is the natural grain here -- units are the
    # outermost thing the integral factorises over, so they are independent by
    # construction -- and `_laplace_parallel` hands them out as workers come
    # free. Cost-weighted chunking is no longer needed for the same reason: the
    # queue does the balancing, and it does it against what the work actually
    # cost rather than against an estimate of it.
    nchunks = nunits
    ranges = [U:U for U in 1:nunits]
    _laplace_ensure_pool!(laplace)
    # Per *slot*, not per piece. Pieces are now one per unit so that the queue
    # can balance, and a model whose units are subjects has thousands of them;
    # anything sized by the piece count would then allocate thousands of
    # vectors per evaluation. Workers are bounded by the pool, so keying on the
    # slot bounds this whatever the granularity becomes.
    nslot = _laplace_pool_width()
    chunk_ok = fill(true, nslot)
    chunk_bad = fill(NaN, nslot)
    run_primal = function (c)
        aws = _laplace_workspace!(laplace, Float64, length(theta))
        @inbounds for U in ranges[c]
            local b0
            b0 = _laplace_mark()
            _laplace_solve_unit_mode!(laplace, U, theta, Ls)
            _laplace_charge!(_LAPLACE_BYTES_MODE, b0)
            u = laplace.modes[U]
            blocks = laplace.units.blocks[U]
            b0 = _laplace_mark()
            M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
                _laplace_unit_curvature(laplace, U, theta, Ls, u)
            _laplace_charge!(_LAPLACE_BYTES_CURV, b0)
            fac = _laplace_factor_repaired!(M, blocks)
            laplace.mode_repaired[U] = fac.repaired
            ok, logdetM, factors, coupling = fac.ok, fac.logdet, fac.factors, fac.coupling
            primal_curvature[U] = (factors, coupling)
            primal_matrices[U] = M
            b0 = _laplace_mark()
            inner = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
            _laplace_charge!(_LAPLACE_BYTES_OBJ, b0)
            laplace.logdet_floored[U] =
                ok && _LAPLACE_PRIOR_FLOOR[] && logdetM < 0
            term = if !isfinite(inner.value) || isempty(u)
                inner.value
            elseif ok
                inner.value - _laplace_prior_floor_logdet(logdetM) / 2
            else
                NaN
            end
            if gated_floor && ok && !isempty(u) && isfinite(inner.value)
                local gated
                # The factorization shifted a repaired `M` in place, so the
                # gated term gets the curvature afresh and factors it itself.
                gated = fac.repaired ?
                    _laplace_gated_term(laplace, U, theta, Ls, u, aws;
                        M=_laplace_unit_curvature(laplace, U, theta, Ls, u),
                        inner=inner) :
                    _laplace_gated_term(laplace, U, theta, Ls, u, aws; M=M,
                        logdetM=logdetM, inner=inner)
                if gated.flagged
                    term = gated.value
                    primal_gated[U] = true
                    # Scored by the gated rule, not floored: the diagnostics
                    # count `logdet_floored` as "the term is a bound".
                    laplace.logdet_floored[U] = false
                elseif !gated.floor && !fac.repaired
                    # Cleared by the gate: plain Laplace, no floor, and so the
                    # ordinary seeded gradient rather than the floored one.
                    term = inner.value - logdetM / 2
                    laplace.logdet_floored[U] = false
                end
            end
            if !isfinite(term)
                chunk_ok[_laplace_slot()] = false
                chunk_bad[_laplace_slot()] = term
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
    _laplace_parallel(laplace, 1:nchunks) do c
        run_primal(c)
        return true
    end
    @inbounds for c in 1:nslot
        chunk_ok[c] || return (value=chunk_bad[c],
            gradient=gradient ? fill(NaN, length(theta)) : nothing,
            subject_loglik=subject_loglik, unit_loglik=unit_loglik,
            converged=all(laplace.inner_converged))
    end
    value = sum(unit_loglik) + _ctsem_log_prior(laplace.objective, theta)
    laplace.gated_units = count(primal_gated)

    gradient || return (value=value, gradient=nothing,
        subject_loglik=subject_loglik, unit_loglik=unit_loglik,
        converged=all(laplace.inner_converged))

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
        # One accumulator per slot rather than one shared vector: the unit
        # contributions are a sum, and summing per slot and then across slots
        # is the same sum in a different order.
        #
        # Indexed by `_laplace_slot()` in *both* branches below, never by the
        # piece `c`. There are `nunits` pieces and `nslot` accumulators, so a
        # piece index past the pool width reads beyond the end of `partials`
        # -- under `@inbounds` that is not a BoundsError but a garbage array
        # handed to the floored gradient, and the process dies with
        # EXCEPTION_ACCESS_VIOLATION and no message on the R side. It also
        # races: piece `c` and the worker holding slot `c` would share one
        # accumulator.
        partials = [zeros(Float64, length(theta)) for _ in 1:nslot]
        fill!(chunk_ok, true)
        run_gradient = function (c)
            @inbounds for U in ranges[c]
                # A unit the gated floor scored is differentiated through its
                # own term; see `_laplace_gated_unit_gradient`.
                if primal_gated[U]
                    local gg
                    gg = _laplace_gated_unit_gradient(laplace, U, theta,
                        primal_curvature[U])
                    if !all(isfinite, gg)
                        chunk_ok[_laplace_slot()] = false
                        return nothing
                    end
                    partials[_laplace_slot()] .+= gg
                    continue
                end
                # A floored unit is differentiated whole; see
                # `_laplace_floored_unit_gradient!`.
                if laplace.logdet_floored[U]
                    if !_laplace_floored_unit_gradient!(
                            partials[_laplace_slot()], laplace, U,
                            theta, Ls, dLlevels)
                        chunk_ok[_laplace_slot()] = false
                        return nothing
                    end
                    continue
                end
                factors, elim = primal_curvature[U]
                local bg
                bg = _laplace_mark()
                if !_laplace_seeded_unit_gradient!(partials[_laplace_slot()],
                        laplace, U, theta,
                        Ls, dLlevels, primal_matrices[U], factors, elim)
                    chunk_ok[_laplace_slot()] = false
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
                _laplace_charge!(_LAPLACE_BYTES_GRAD, bg)
                if !all(isfinite, partials[_laplace_slot()])
                    chunk_ok[_laplace_slot()] = false
                    return nothing
                end
            end
            return nothing
        end
        _laplace_parallel(laplace, 1:nchunks) do c
            run_gradient(c)
            return true
        end
        # Not `ok`: `run_primal` assigns that name, and a closure binds an
        # enclosing local rather than shadowing it, so a function-level `ok`
        # here made every concurrent unit share one boxed `ok` between its
        # factorization and its term.
        gradient_ok = all(chunk_ok)
        if gradient_ok
            for w in 1:nslot
                grad .+= partials[w]
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
        unit_loglik=unit_loglik, converged=all(laplace.inner_converged))
end

"""
    _laplace_floored_unit_gradient!(out, laplace, U, values, Ls, dL)

The gradient of one floored unit's term.

Where the prior floor binds, the term is `g_U(uhat)` and carries no `logdet`,
so two things follow and both of them simplify this. The derivative of
`logdet` is not part of the gradient -- and the seeded assembly's trace terms
*are* that derivative, taken from the selected inverse of a curvature the floor
has just declared untrustworthy. And `dg_U/du = 0` at the mode, so by the
envelope theorem the mode's own movement contributes nothing either: what is
left is the partial derivative at fixed `uhat`.

That is one reverse sweep per member and the `dL/dtheta` chain, against the
`k + 1` sweeps and the selected inverse the seeded route needs -- cheaper than
the path it replaces rather than dearer. Differentiating the term with
ForwardDiff instead is correct and was measured: it costs `O(npar * k)` sweeps
per unit, and since the floor binds at *every* iteration while the optimiser is
in that region, it turned a fit of seconds into 350-480s.

# The chain

`_laplace_member_values!` builds member `m`'s parameter vector as

    shifted = values;  shifted[re_index[p]] += sum_q L_l[p, q] u[base_l + q]

so with `grad = d loglik_m / d shifted` from the reverse sweep,

    d loglik_m / d theta_t  =  grad[t]
        + sum_l sum_p grad[re_index[p]] sum_q (dL_l[p, q] / d theta_t) u[base_l + q]

The first term is every parameter's direct appearance, the TI coefficients
included, because those are entries of `shifted` and the sweep has already
accounted for them. The second is the only other route theta takes into the
term at fixed `u`: through the population Cholesky that scales the random
effects. `-u'u/2` is constant at fixed `u` and contributes nothing.
"""
function _laplace_floored_unit_gradient!(out::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Integer, values::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, dL::Vector{Vector{Matrix{Float64}}})
    spec = laplace.spec
    units = laplace.units
    members = units.members[U]
    u = laplace.modes[U]
    npar = length(values)
    aws = _laplace_workspace!(laplace, Float64, npar)
    grad = Vector{Float64}(undef, npar)
    shift = Vector{Float64}(undef, npar)
    # Once per level, not once per member per level: `_laplace_level_positions`
    # allocates, and it depends on the spec alone.
    levelpositions = [_laplace_level_positions(spec, l)
                      for l in eachindex(spec.levels)]
    @inbounds for m in eachindex(members)
        i = members[m]
        offsets = units.offsets[U][m]
        shifted = _laplace_member_values!(shift, values, spec, Ls, u, offsets)
        loglik = _laplace_subject_value_gradient!(grad,
            laplace.objective.subject_objectives[i], aws, shifted)
        isfinite(loglik) || return false
        all(isfinite, grad) || return false
        for t in 1:npar
            out[t] += grad[t]
        end
        for l in eachindex(spec.levels)
            level = spec.levels[l]
            k = nrandomeffects(level)
            r = nlatent(level)
            (k == 0 || r == 0) && continue
            positions = levelpositions[l]
            isempty(positions) && continue
            base = offsets[l]
            dLl = dL[l]
            for t in eachindex(positions)
                dLt = dLl[t]
                acc = 0.0
                for pp in 1:k
                    gp = grad[level.re_index[pp]]
                    iszero(gp) && continue
                    inner = 0.0
                    for q in 1:r
                        inner += dLt[pp, q] * u[base + q]
                    end
                    acc += gp * inner
                end
                out[positions[t]] += acc
            end
        end
    end
    return true
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
`_laplace_member_values` and the engine's own `_materialize_subject_values!`
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

"""
    ctsem_laplace_subject_values(laplace, values, effects)

Each subject's own raw parameter vector at the random effects it is *given*,
rather than at the conditional modes solved from the data.

The method above estimates each unit's effects as a mode. That is the right
answer for a fit that integrated them out, and the wrong one for a fit that
sampled them: a sampler already has a draw of every effect, and pairing draw
`s` of the parameters with draw `s` of the effects is what makes a posterior
predictive integrate over the random effects instead of conditioning on a
point estimate of them.

`effects` is the flat vector one draw of the sampler produces: unit by unit in
unit order, and within a unit by block offset, exactly the layout
`ctsem_laplace_effect_layout` reports and the sampler returns. Its length must
be the total latent dimension over all units.

Everything after the slicing is shared with the mode method -- the same
`_laplace_popchols`, `_laplace_member_values` and
`_materialize_subject_values!` -- so the two cannot drift in how a shifted
parameter vector is built.

Returns `nsubjects x length(values)`.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, effects::AbstractVector)
    theta = collect(Float64, values)
    Ls = _laplace_popchols(theta, laplace.spec)
    units = laplace.units
    total = sum(units.dims; init=0)
    length(effects) == total || throw(DimensionMismatch(string(
        "effects must have ", total, " entries for this design, got ",
        length(effects))))
    nsubjects = length(laplace.objective.subject_objectives)
    out = zeros(Float64, nsubjects, length(theta))
    buffer = Float64[]
    base = 0
    for U in eachindex(units.members)
        u = Vector{Float64}(view(effects, (base + 1):(base + units.dims[U])))
        for (m, i) in enumerate(units.members[U])
            shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                units.offsets[U][m])
            subject = laplace.objective.subject_objectives[i]
            _materialize_subject_values!(buffer, shifted, subject.params,
                subject.tipreds)
            out[i, :] = buffer
        end
        base += units.dims[U]
    end
    return out
end

export ctsem_laplace_subject_values

"""
    ctsem_state_dimension(laplace)

How many innovations the design needs. That does not depend on how the random
effects are handled -- it is a property of the subjects, their rows and their
substeps -- so a caller holding a Laplace objective can size an innovation
vector without unwrapping it.
"""
ctsem_state_dimension(laplace::CTSEMLaplaceObjective) =
    ctsem_state_dimension(laplace.objective)

"""
    ctsem_state_layout(laplace)

The inner objective's layout. A Laplace fit's random effects are not carrier
states, so the innovation vector describes the dynamic states alone.
"""
ctsem_state_layout(laplace::CTSEMLaplaceObjective) =
    ctsem_state_layout(laplace.objective)

"""
    ctsem_generate_states(laplace, values, z, base; effects, subject_values)

One dataset with each subject's trajectory drawn from *its own* model.

`values` is the population parameter vector. `effects` is one draw of the
random effects in the sampler's layout, from which each subject's own
parameter vector is built; `subject_values` supplies those vectors directly
instead. Given neither, each subject is put at its conditional mode, matching
what `ctsem_generate` does.

The states are drawn afresh here -- from the process, at that subject's own
parameters -- and never carried over from anything a fit sampled. What comes
from a fit is the individual differences; the trajectory is regenerated
conditional on the subject-specific model.
"""
function ctsem_generate_states(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, z::AbstractVector, base::AbstractMatrix;
    effects::Union{Nothing,AbstractVector}=nothing,
    subject_values::Union{Nothing,AbstractMatrix}=nothing, kwargs...)
    persubject = subject_values !== nothing ? subject_values :
        effects !== nothing ?
            ctsem_laplace_subject_values(laplace, values, effects) :
            ctsem_laplace_subject_values(laplace, values)
    return ctsem_generate_states(laplace.objective, persubject, z, base; kwargs...)
end

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
            k = nlatent(laplace.spec.levels[l])
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

    # Chunked, like the identical work inside `run_primal`.
    #
    # This loop is every unit's mode solve followed by its curvature and
    # factorization -- the same six lines `ctsem_laplace_evaluate` runs across
    # chunks. Here it ran serially, and worse, it omitted the `slot` argument
    # that both inner calls take, so every unit shared adjoint workspace 1.
    # That is what made it unthreadable rather than merely unthreaded:
    # `_laplace_solve_unit_mode!` and `_laplace_unit_curvature` both default
    # `slot` to 1, and two units on one slot corrupt each other's workspace.
    #
    # Its two callers, `ctsem_laplace_mode_jacobian` and
    # `ctsem_subject_gradients`, are top-level entry points rather than
    # anything reached from inside a chunk, so there is no nesting to worry
    # about. Each unit is written by exactly one chunk -- `out[U]`,
    # `laplace.modes[U]`, `laplace.mode_repaired[U]` -- so the only sharing was
    # the workspace, and a slot per chunk removes it.
    nchunks = _ctsem_nchunks(nunits)
    # Slot bands. This function sizes `laplace.workspaces`, so this function is
    # what decides how many slots each chunk may divide: chunk `c` owns
    # `(c-1)*width+1` upward for `width` slots, and hands that band to anything
    # below it that can spend it. Nothing downstream invents a slot, which is why
    # the routes that size this vector themselves -- quadrature, and the
    # sampler's per-chain blocks -- are unaffected by any of it.
    _laplace_ensure_pool!(laplace)
    # Written out rather than held in a closure. `width` is captured by the
    # per-chunk worker either way, and a captured variable that a closure also
    # reads is the shape Julia boxes -- which turns every call that takes it
    # into a dynamic one. Two integers are not worth a closure.
    ranges = nchunks > 1 ?
        _ctsem_chunk_assignment(_laplace_unit_weights(laplace), nchunks) :
        [1:nunits]

    run = function (c)
        @inbounds for U in ranges[c]
            _laplace_solve_unit_mode!(laplace, U, theta, Ls)
            u = laplace.modes[U]
            blocks = laplace.units.blocks[U]
            M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
                _laplace_unit_curvature(laplace, U, theta, Ls, u)
            fac = _laplace_factor_repaired!(M, blocks)
            laplace.mode_repaired[U] = fac.repaired
            out[U] = (fac.factors, fac.coupling)
        end
        return nothing
    end
    _laplace_parallel(laplace, 1:nchunks) do c
        run(c)
        return true
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
    ctsem_laplace_population(laplace, values, level=1)

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
    ctsem_laplace_modes(laplace, values, level=1)

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
    # `z` lives in the block's own space and `raw` in the parameter space the
    # loading maps it to; for a reduced level those have different widths.
    r = nlatent(lv)
    Ls = _laplace_popchols(theta, spec)
    L = Ls[level]
    ngroups = lv.ngroups
    z = zeros(Float64, ngroups, r)
    raw = zeros(Float64, ngroups, k)
    zsd = zeros(Float64, ngroups, r)
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
            (k == 0 || r == 0) && continue
            base = laplace.units.offsets[U][m][level]
            slice = (base + 1):(base + r)
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

"""Smallest eigenvalue of `M` below which a unit is reported near-singular."""
const _LAPLACE_NEAR_SINGULAR = 0.05

"""
    ctsem_laplace_conditioning(laplace)

How each unit's inner curvature `M = -d2 log L/du du + I` stands against the
prior's own, at the last evaluation's parameters and modes.

`M >= I` is the likelihood being concave in `u`, where the Laplace term is
well founded. An eigenvalue below one is a direction where it has gone convex:
mildly, which the Gaussian at that curvature still integrates well, or nearly
to zero, where the Gaussian extrapolates across a range the prior forbids and
the term is several nats too high. That second case is invisible otherwise --
the total prior floor does not engage while the other eigenvalues keep
`logdet M` positive, and an optimiser can park a unit exactly at `logdet = 0`
-- so it is counted here. See
`CT-SEM/review/LAPLACE-eigenwise-floor-2026-09-23.md`.

Report-time only, and it changes nothing: it rebuilds each unit's curvature
once, tests `M - I` with the block factorization (no fill-in), and
decomposes densely only the units that fail that test and are no wider than
`_LAPLACE_EIGEN_MAXDIM`. Returns

  * `min_eigenvalue`, per unit: the smallest eigenvalue where the test
    failed, `Inf` where it passed (every eigenvalue above one, not computed),
    and `NaN` for a failed unit too wide to decompose;
  * `below_one`, units with an eigenvalue below one (undecomposed failures
    included);
  * `near_singular`, units whose smallest eigenvalue is below 0.05.
"""
function ctsem_laplace_conditioning(laplace::CTSEMLaplaceObjective)
    nunits = length(laplace.units.members)
    mins = fill(Inf, nunits)
    theta = laplace.last_values
    if length(theta) != 0
        _laplace_ensure_pool!(laplace)
        Ls = _laplace_popchols(theta, laplace.spec)
        for U in 1:nunits
            u = laplace.modes[U]
            isempty(u) && continue
            blocks = laplace.units.blocks[U]
            M = _laplace_unit_curvature(laplace, U, theta, Ls, u)
            all(d -> all(isfinite, d), M.diag) || (mins[U] = NaN; continue)
            _laplace_exceeds_identity(M, blocks) && continue
            d = length(u)
            if d > _LAPLACE_EIGEN_MAXDIM[]
                mins[U] = NaN
                continue
            end
            mins[U] = _ctsem_symeig(_laplace_block_dense(M, blocks, d)).values[1]
        end
    end
    return (min_eigenvalue=mins,
        below_one=count(x -> isnan(x) || x < 1, mins),
        near_singular=count(x -> x < _LAPLACE_NEAR_SINGULAR, mins))
end

export ctsem_laplace_conditioning

"""
    ctsem_laplace_diagnostics(laplace)

Inner-solve status from the last evaluation, per subject, with the curvature
conditioning `ctsem_laplace_conditioning` reports.
"""
function ctsem_laplace_diagnostics(laplace::CTSEMLaplaceObjective)
    conditioning = ctsem_laplace_conditioning(laplace)
    return (
        iterations=copy(laplace.inner_iterations),
        max_gradient=copy(laplace.inner_gradient),
        converged=copy(laplace.inner_converged),
        hessian_repaired=copy(laplace.hessian_repaired),
        mode_repaired=copy(laplace.mode_repaired),
        logdet_floored=copy(laplace.logdet_floored),
        min_eigenvalue=conditioning.min_eigenvalue,
        below_one=conditioning.below_one,
        near_singular=conditioning.near_singular,
        floor=laplace.floor,
        gated_units=laplace.gated_units,
    )
end

export ctsem_laplace_diagnostics

"""Per-run state: the fallback counter this route keeps for its own report."""
function _ctsem_optimise_setup!(o::CTSEMLaplaceObjective)
    _CTSEM_LAPLACE_FALLBACKS[] = 0
    # The modes are a warm start *within* one optimisation, where consecutive
    # parameter vectors are close and the Newton solve lands in one or two
    # steps. Between optimisations they are nothing of the kind: R caches this
    # object by a hash of the model spec, so without this the next `ctFit()` of
    # the same model begins at wherever the last one left its units -- and a
    # previous fit that failed leaves them somewhere its successor cannot
    # recover from. Measured on a 40-subject count model: from zeros in a clean
    # session, -1471.584 in 28 iterations; from the same start after a failed
    # laplace fit, -6.4e86 in 2.
    #
    # Zero is what the constructor uses, so this is the state a fit would have
    # had if it were the first one in the session -- which is the property
    # being restored.
    for mode in o.modes
        fill!(mode, 0.0)
    end
    # The diagnostics describe the run that is about to happen, and
    # `_ctsem_optimise_verbose_shape` reads `inner_converged` before the first
    # evaluation of it -- so leaving the previous fit's values here reports one
    # fit's inner solve as another's.
    fill!(o.inner_iterations, 0)
    fill!(o.inner_gradient, 0.0)
    fill!(o.inner_converged, false)
    fill!(o.hessian_repaired, false)
    fill!(o.mode_repaired, false)
    fill!(o.logdet_floored, false)
    return nothing
end

_ctsem_optimise_label(::CTSEMLaplaceObjective) = "Laplace"

"""
Why trial points were rejected, counted rather than guessed at.

A fit that stops short of a stationary point almost always did so because its
line search ran out of points it was allowed to accept, and these say which of
the three reasons was doing it.
"""
mutable struct CTSEMLaplaceCallLog
    rejected_nonfinite::Int
    rejected_inner::Int
    rejected_gradient::Int
    accepted::Int
    verbose::Bool
    reported_gradient::Bool
end

_ctsem_optimise_log(::CTSEMLaplaceObjective, verbose::Bool) =
    CTSEMLaplaceCallLog(0, 0, 0, 0, verbose, false)

# The inner mode count is worth tracing and reporting here and not on the
# marginal route: an iteration re-solving every unit's mode from scratch costs
# an order of magnitude more than one warm-starting, and the difference shows up
# as a stall the objective alone does not explain.
_ctsem_optimise_trace_keys(::CTSEMLaplaceObjective) =
    (:objective, :gradient_norm, :inner_converged)
_ctsem_optimise_trace_values(o::CTSEMLaplaceObjective) =
    (count(o.inner_converged),)
_ctsem_optimise_progress_extra(o::CTSEMLaplaceObjective) =
    (@sprintf("inner %d/%d", count(o.inner_converged),
        length(o.inner_converged)),)

# A level of correlations can saturate together here, which the per-cell check
# cannot see. See `_laplace_saturated_parameters`.
_ctsem_saturated_for(o::CTSEMLaplaceObjective, minimizer) =
    _laplace_saturated_parameters(o, minimizer)

"""
A trial point on the laplace route, where finite is not the same as usable.

`converged` is the inner solve's verdict, and a point it failed at is invalid
however finite its value: the Laplace term is defined at the mode, so away from
one the number is not the objective. `ctsem_laplace_evaluate` is called rather
than the generic `ctsem_evaluate`, which deliberately drops that flag.

`gradient_method === :nested` is how `nested_gradient` reaches here through the
shared driver's one gradient-method argument.
"""
function _ctsem_optimise_trial(o::CTSEMLaplaceObjective, x, want_gradient::Bool,
        gradient_method, limit::Real, log)
    evaluated = try
        ctsem_laplace_evaluate(o, x; gradient=want_gradient,
            nested_gradient=(Symbol(gradient_method) === :nested))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    value = evaluated === nothing ? NaN : evaluated.value
    finite_value = evaluated !== nothing && isfinite(value)
    inner_ok = evaluated !== nothing && evaluated.converged
    valid = finite_value && inner_ok
    if valid && want_gradient
        valid = all(isfinite, evaluated.gradient) &&
            all(abs(entry) < limit for entry in evaluated.gradient)
        if !valid
            log.rejected_gradient += 1
            # Once, and only if asked: the probe prints the whole raw vector
            # and the whole gradient, which is what makes it worth having and
            # what makes it unwelcome on a quiet fit.
            if log.verbose && !log.reported_gradient
                log.reported_gradient = true
                finite_part = filter(isfinite, evaluated.gradient)
                println(_console(), "Laplace probe: first gradient rejection, objective ",
                    value, ", ", count(!isfinite, evaluated.gradient), " of ",
                    length(evaluated.gradient), " entries non-finite, largest finite ",
                    isempty(finite_part) ? 0.0 : maximum(abs, finite_part))
                println(_console(), "Laplace probe: at x = ", collect(x))
                println(_console(), "Laplace probe: gradient = ",
                    collect(evaluated.gradient))
            end
        end
    end
    finite_value || (log.rejected_nonfinite += 1)
    (finite_value && !inner_ok) && (log.rejected_inner += 1)
    valid && (log.accepted += 1)
    return (evaluated=evaluated, valid=valid)
end

"""
A probe point on the laplace route, where finite is not the same as usable.

The same predicate `_ctsem_optimise_trial` applies to a trial point, for the
same reason: the Laplace term is defined at the mode, so a value computed where
a unit's inner Newton did not reach one is not the objective. `_ctsem_overshot`
compares probe values against the estimate's and reports an improvement as
proof that the estimate is not a maximum -- so a number from a failed inner
solve there does not merely add noise, it manufactures that proof.

It is not hypothetical on this route: the probe deliberately evaluates far from
the estimate, pulling coordinates to zero, which is exactly where a warm inner
solve used to fail.
"""
function _ctsem_probe_value(o::CTSEMLaplaceObjective, x)
    evaluated = try
        ctsem_laplace_evaluate(o, x; gradient=false)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    evaluated === nothing && return -Inf
    (evaluated.converged && isfinite(evaluated.value)) ? evaluated.value : -Inf
end

"""The laplace route's own result fields, and no `row_loglik`.

The integral is over a whole subject's trajectory, so a single row has no
marginal contribution to report: returning the subject terms and omitting the
row ones is more honest than inventing a decomposition the approximation does
not have.
"""
_ctsem_optimise_result_extra(o::CTSEMLaplaceObjective, final, log) = (
    inner_converged=all(o.inner_converged),
    inner_failures=count(!, o.inner_converged),
    rejected_nonfinite=log.rejected_nonfinite,
    rejected_inner=log.rejected_inner,
    rejected_gradient=log.rejected_gradient,
    accepted_calls=log.accepted,
    inner_worst_gradient=isempty(o.inner_gradient) ? 0.0 :
        maximum(o.inner_gradient),
    inner_iterations=copy(o.inner_iterations),
    hessian_repaired=copy(o.hessian_repaired),
    mode_repaired=copy(o.mode_repaired),
    # At the final evaluation, which is the one just made at the minimizer;
    # flat rather than nested so it crosses the bridge as plain fields.
    unit_min_eigenvalue=ctsem_laplace_conditioning(o).min_eigenvalue,
    laplace_floor=String(o.floor),
    gated_units=o.gated_units,
    # How many gradients were computed twice. The seeded assembly returns
    # `false` on a failed factorization and the caller silently recomputes the
    # whole thing by the nested route, which is correct and much slower -- and
    # invisible, because the answer is right either way. A reduced-rank level
    # once failed on *every* gradient for a structural reason and the only
    # symptom was a fit that crawled. Counted already for the verbose report;
    # carried out here so a caller can see it without asking for one.
    gradient_fallbacks=_CTSEM_LAPLACE_FALLBACKS[])

"""
What is about to be fitted: the sizes that decide what the fit will cost.

The unit and level counts are what a misdeclared grouping shows up as -- one
unit per subject when a study level was meant -- and the chunk and thread counts
are the difference between a run that uses the machine and one that does not.
Printed after chunk tuning, which is what decides the chunk count.
"""
function _ctsem_optimise_verbose_shape(o::CTSEMLaplaceObjective)
    chunks = ctsem_max_chunks()
    println(_console(), "Laplace: ", length(o.objective.subject_objectives),
        " subjects in ", length(o.units.members), " unit(s), ",
        nlevels(o.spec), " level(s), ", nrandomeffects(o.spec),
        " random effects, ", min(max(chunks.max_chunks == 0 ? chunks.nthreads :
            chunks.max_chunks, 1), length(o.units.members)),
        " chunk(s) over ", chunks.nthreads, " thread(s)")
    return nothing
end

"""
What the run did, for a user who asked.

Two things that are otherwise invisible and both of which explain a fit that
stopped short. The first is why trial points were refused, which is the
difference between a line search that ran out of room and one that was never
offered a usable point. The second is the inner solve: the Laplace term is
defined at each unit's mode, so a run where modes went unconverged or needed
their curvature repaired is one whose objective was approximate in a way the
final numbers do not show.

Asserted in `test-julia-laplace.R`, and that is deliberate -- the engine suite
proves the inner solve works, and this is the only check anywhere that its
status reaches the user who asked for it.
"""
function _ctsem_optimise_verbose_report(o::CTSEMLaplaceObjective,
        log::CTSEMLaplaceCallLog)
    println(_console(), "Laplace: ", log.accepted, " objective evaluations accepted, ",
        log.rejected_nonfinite, " rejected as non-finite, ", log.rejected_inner,
        " for an inner mode solve that did not converge, ",
        log.rejected_gradient, " for the gradient; ",
        _CTSEM_LAPLACE_FALLBACKS[],
        " gradient(s) fell back to the nested route")
    nfloored = count(o.logdet_floored)
    nfloored > 0 && println(_console(), "Laplace: the prior floor bound for ",
        nfloored, "/", length(o.logdet_floored), " unit(s) at the last ",
        "evaluation -- their posterior is wider than the prior in some ",
        "direction, so the term reported for them is a bound rather than the ",
        "Laplace value. ctLaplaceCheck() measures the gap.")
    println(_console(), "Laplace: inner modes ",
        count(o.inner_converged), "/", length(o.inner_converged),
        " converged, max |dg/dz| ",
        isempty(o.inner_gradient) ? 0.0 : maximum(o.inner_gradient),
        ", curvature repaired at the mode for ", count(o.mode_repaired),
        " unit(s) (", count(o.hessian_repaired), " somewhere on the way)")
    conditioning = ctsem_laplace_conditioning(o)
    if conditioning.below_one > 0
        finite = filter(isfinite, conditioning.min_eigenvalue)
        println(_console(), "Laplace: ", conditioning.below_one, "/",
            length(conditioning.min_eigenvalue),
            " unit(s) with curvature below the prior's, ",
            conditioning.near_singular, " near-singular",
            isempty(finite) ? "" : string(" (smallest eigenvalue ",
                round(minimum(finite); sigdigits=3), ")"))
    end
    o.floor === :gated && println(_console(), "Laplace: gated floor, ",
        o.gated_units, " unit(s) scored by its soft-direction rule at the last ",
        "evaluation")
    return nothing
end

"""
    ctsem_laplace_optimize(laplace, start; ...)

Maximize a Laplace-approximated multilevel likelihood with L-BFGS.

The Laplace value and the Laplace gradient are a consistent pair, so the line
search behaves and convergence means something.

There used to be a second mode here that maximized the value without its
log-determinant, using the cheaper envelope gradient. It is gone. Dropping the
log-determinant leaves the PQL-shaped objective, which is degenerate in the
variance components -- nothing in it penalises the population scales growing,
so they run away. A nine-replication recovery study put its population sd at
19.6 against a truth of 1.0, and one replication failed outright. It was not a
cheaper route to the same answer, so there is no version of it worth keeping.

The driver is `ctsem_optimize`, which this route used to carry a copy of. The
copy cost three bugs -- the convergence verdict was fixed on one route and left
on the other, the saturation guard had to be copied across afterwards, and
`gap_tol` went into one and killed every fit here with a MethodError -- and what
it differed by is the eleven methods above.

`nested_gradient` travels as `gradient_method`, which is the shared driver's one
name for the same choice.
"""
function ctsem_laplace_optimize(laplace::CTSEMLaplaceObjective,
        start::AbstractVector; nested_gradient::Bool=false, kwargs...)
    ctsem_optimize(laplace, start;
        gradient_method=(nested_gradient ? :nested : :adjoint), kwargs...)
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

A column whose two evaluations did not both solve their inner modes is `NaN`,
for the reason `ctsem_laplace_optimize`'s `fg!` gives for rejecting such a
point: a unit that has not reached its mode gives whatever its iteration
stopped on, not the objective, so the difference is not a derivative of the
objective either. The caller sees a non-finite Hessian and falls back.

Before writing that `NaN`, the inner budget is raised and the point evaluated
again -- see `_laplace_hessian_point`. A whole Hessian of `NaN` columns is not
a usable fallback for anything, and the differencing points are `1e-4` from one
where the modes were found, so a unit that merely ran out of iterations there
will reach its mode given more. `retries` caps that at `4^retries` times the
objective's own budget.

`warm` (the default) starts every differencing point's inner modes from the
modes at `values`, solved there from the origin first, rather than from the
origin. It is the one place a warm start is sound: every point is `1e-4` from
the same base, so the start is fixed and the objective does not depend on the
order points are visited in -- which is what `_laplace_solve_unit_mode!`
refuses a warm start for elsewhere. Measured on dev1: warm and cold Hessians
agree to 1e-8 relative or better, standard errors to 2e-8, including on the
nonlinear fixture whose inner problem is multimodal away from the mode; inner
iterations fall from 3.0-3.7 to 2.0 per unit, the Hessian 1.3-1.6x faster on
ordinal and nonlinear models and unchanged on a linear-Gaussian one, which
already needed two. `warm = false` is the behaviour before, which the
inner-budget escalation test needs.
"""
function ctsem_laplace_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    step::Real=1e-4, retries::Integer=2, warm::Bool=true)
    x = collect(Float64, values)
    n = length(x)
    H = zeros(Float64, n, n)
    budget = laplace.inner_maxiter
    base = nothing
    if warm
        # The modes at `x`, from the origin, which every differencing point
        # then starts from.
        previous = ctsem_set_warm_start!(false)
        try
            ctsem_laplace_evaluate(laplace, x; gradient=false)
        finally
            ctsem_set_warm_start!(previous)
        end
        base = deepcopy(laplace.modes)
    end
    function restore!()
        base === nothing && return nothing
        for U in eachindex(base)
            laplace.modes[U] = copy(base[U])
        end
        return nothing
    end
    previous = ctsem_set_warm_start!(warm)
    unconverged = Int[]
    escalated = Int[]
    worst = 0.0
    units = 0
    try
        for j in 1:n
            h = step * max(1.0, abs(x[j]))
            plus = copy(x); plus[j] += h
            minus = copy(x); minus[j] -= h
            restore!()
            pp = _laplace_hessian_point(laplace, plus, budget, retries)
            restore!()
            pm = _laplace_hessian_point(laplace, minus, budget, retries)
            (pp.tries + pm.tries) > 0 && push!(escalated, j)
            if !(pp.out.converged && pm.out.converged)
                # The same predicate `fg!` applies to a trial point, applied to
                # the two points this column is differenced from. A gradient at
                # a not-quite-mode is a gradient of a different function, and
                # one that converged to 1e-9 instead of 1e-10 is finite, so
                # nothing downstream would notice. NaN makes the column visible
                # to the non-finite fallback callers already have.
                push!(unconverged, j)
                worst = max(worst, pp.worst, pm.worst)
                units = max(units, pp.units, pm.units)
                H[:, j] .= NaN
                continue
            end
            H[:, j] = (pp.out.gradient .- pm.out.gradient) ./ (2h)
        end
    finally
        # The budget is a field of a mutable objective the caller keeps, so an
        # escalation that threw would otherwise be inherited by every later
        # evaluation of it; the warm-start switch is global, and the same.
        laplace.inner_maxiter = budget
        ctsem_set_warm_start!(previous)
    end
    # A count, not a list of indices: on a model with a thousand parameters the
    # list is the whole message and says no more than the count does.
    isempty(escalated) || @info string("Laplace inner budget raised above ",
        budget, " for ", length(escalated), " of ", n, " Hessian columns.")
    isempty(unconverged) || @warn string("Laplace inner solve did not converge ",
        "at the Hessian step for parameter ", join(unconverged, ", "),
        "; those columns are NaN. ", units, " unit(s) short, worst inner ",
        "gradient ", worst, " against a budget of ", budget * 4^retries,
        " iterations.")
    return (H .+ transpose(H)) ./ 2
end

"""
    _laplace_hessian_point(laplace, x, budget, retries)

One of the two points a Hessian column is differenced from, evaluated with the
inner budget raised rather than the column abandoned.

Retried only while some unit actually *exhausted* its iterations. Every other
way the inner solve reports failure is deterministic from a fixed start -- the
start is the origin, so a repeat evaluation at the same `x` visits the same
points and stops the same way -- and the retry would buy nothing but the cost
of another pass. A unit that ran out of iterations is the one case where a
larger budget changes the answer.

Reports what it reached as well as whether it got there: `worst` is the largest
inner gradient over units and `units` how many are still short, which is what
tells a caller whether raising `inner_maxiter` would help. (R spells that
optimcontrol name `laplace_inner_maxiter`; a bare dollar sign interpolates in a
Julia docstring, which is why it is not written out here.)
"""
function _laplace_hessian_point(laplace::CTSEMLaplaceObjective,
    x::AbstractVector, budget::Integer, retries::Integer)
    capped() = any(>=(laplace.inner_maxiter), laplace.inner_iterations)
    reached() = isempty(laplace.inner_gradient) ? 0.0 :
        maximum(laplace.inner_gradient)
    short() = count(!, laplace.inner_converged)
    laplace.inner_maxiter = budget
    out = ctsem_laplace_evaluate(laplace, x; gradient=true)
    tries = 0
    while !out.converged && capped() && tries < retries
        tries += 1
        laplace.inner_maxiter = budget * 4^tries
        out = ctsem_laplace_evaluate(laplace, x; gradient=true)
    end
    result = (out=out, tries=tries, worst=reached(), units=short())
    laplace.inner_maxiter = budget
    return result
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

Deliberately without `converged`, which `ctsem_laplace_evaluate` does carry: a
caller through this entry point wants a value and a gradient, and the one thing
that needs the inner solve's verdict is the optimiser's own `fg!`, which calls
`ctsem_laplace_evaluate` directly. Merging the two drivers means giving the
shared one a validity hook rather than smuggling the flag through here.
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

# `ctsem_generate_states` and `ctsem_state_dimension` used to be on this list.
# Both now delegate to the wrapped objective: the first because generating a
# trajectory needs a per-subject parameter vector and there is one to give it,
# the second because the innovation count is a property of the design and not
# of how the random effects are handled. What remains is the state-explicit
# *target*, which really is not implemented for this route -- it would have to
# integrate the random effects and sample the states in one objective.
for (f, what) in ((:ctsem_joint_loglikelihood, "The state-explicit path"),
    (:ctsem_joint_evaluate, "The state-explicit path"))
    @eval function $f(laplace::CTSEMLaplaceObjective, args...; kwargs...)
        throw(ArgumentError(string($what, " is not implemented for the Laplace ",
            "random-effect route yet. It needs a per-subject parameter vector ",
            "rather than one shared vector, and returning the population-level ",
            "answer instead would be a different quantity than the one asked ",
            "for. Use intoverpop=TRUE for this.")))
    end
end

"""
    ctsem_generate(laplace, values, base; subject_values=nothing, seed=1)

One posterior-predictive dataset from a Laplace fit.

Each subject is generated at its own realized parameters -- the population
vector shifted by that subject's estimated random effects and TI-predictor
effects -- rather than at the shared population vector, using exactly the
per-subject values `ctsem_kalman(laplace, ...)` already computes for
prediction (`ctsem_laplace_subject_values`). A Laplace random effect is a
conditional mode estimated from the subject's whole record, not a fresh draw
from the population distribution, so the individual difference in the
generated data is the fitted one: this is the same smoothed-equivalent,
conditional-on-the-subject's-own-data quantity `ctsem_kalman(laplace, ...)`
already reports for residuals and predictions on this route, not a new
statistical convention introduced for generation.

`subject_values` supplies those per-subject vectors instead of solving for
them here, for the same reason `ctsem_kalman` takes it: a caller filtering
different rows than the fit did needs the fitted modes, not modes re-solved
against rows that may not condition on anything.
"""
function ctsem_generate(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    base::AbstractMatrix; subject_values::Union{Nothing,AbstractMatrix}=nothing,
    seed::Integer=1)
    persubject = subject_values === nothing ?
        ctsem_laplace_subject_values(laplace, values) : subject_values
    return ctsem_generate(laplace.objective, persubject, base; seed=seed)
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

"""The parameters a saturation check reads, for the wrapped objective."""
_ctsem_params(o::CTSEMLaplaceObjective) = _ctsem_params(o.objective)
