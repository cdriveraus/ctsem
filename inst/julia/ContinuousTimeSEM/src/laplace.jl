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
  * `sd_fixed`, `cor_fixed`: raw values for scales and correlation coordinates
    the model states rather than estimates. A stated entry has index 0 in
    `sd_index` or `cor_index` and its raw value here; a free one has its raw
    position there and `NaN` here. Empty means nothing is stated.
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
    `L * L'`. Its latent block is `nlatent(level)`-dimensional, not
    `k`-dimensional: the deviation a unit contributes is `F * u` for
    `u ~ N(0, I)` (see `_laplace_member_values!`) with `F` the level's root,
    `k x nlatent`, and every site that indexes `u` asks `nlatent`.
  * `unique_index[j]`: raw position of effect `j`'s own variance, or 0 when it
    has none; empty when no effect has one. Only on a reduced level (the rank
    may then be 0): the covariance becomes `L * L' + diag(psi.^2)`, the
    factor-analytic form, and the root `F = [L | Psi]` gains one column per
    own variance, a single entry `psi_j` in row `j`, transformed as a
    full-rank level's scale. An own variance on a full-rank level would not
    be identified, so it is refused.
  * `load_row`, `load_col`: the effect and dimension each loading fills, in
    `load_index` order, when the level states a pattern; empty for the
    triangle above. A pattern places dimensions in blocks, each loaded by
    the effects the model names for it and triangular in its own first rows,
    so an effect named in no block shares no dimension with the others
    (`ctFit(popblocks =)`). A patterned level is a loading level whatever its
    rank.
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
    unique_index::Vector{Int}
    # The order this level's effects take in the Cholesky root that standardises
    # them, `d = S u` (`_laplace_popchol`); empty for the natural order. Set by
    # the samplers for the run (`_laplace_order_roots!`) and cleared after it,
    # so nothing outside a sampling run sees permuted coordinates.
    order::Vector{Int}
    sd_fixed::Vector{Float64}
    cor_fixed::Vector{Float64}
    # Effects whose own variance the subject filter integrates rather than this
    # level (`CTSEMHybrid`): no column of the level's root, no coordinate of
    # its block. Empty when there are none.
    hybrid::Vector{Bool}
    # The effect and dimension of each loading, when the level states a
    # pattern (`_laplace_loading_entries`); empty for the triangle.
    load_row::Vector{Int}
    load_col::Vector{Int}

    function CTSEMLaplaceLevel(re_index, sd_index, cor_index, sd_scale, group,
        ngroups::Integer; covmatcode::Integer=0, rank::Integer=-1,
        load_index=Int[], sd_fixed=Float64[], cor_fixed=Float64[],
        unique_index=Int[], hybrid=Bool[], load_row=Int[], load_col=Int[])
        re = Vector{Int}(collect(re_index))
        sd = Vector{Int}(collect(sd_index))
        cor = Vector{Int}(collect(cor_index))
        load = Vector{Int}(collect(load_index))
        scale = Vector{Float64}(collect(sd_scale))
        grp = Vector{Int}(collect(group))
        lrow = Vector{Int}(collect(load_row))
        lcol = Vector{Int}(collect(load_col))
        patterned = !isempty(lrow)
        k = length(re)
        # A negative rank means "not stated", which is the full-rank level every
        # caller before this argument existed builds. A pattern may span more
        # dimensions than effects: an effect can load on several blocks.
        r = rank < 0 ? k : Int(rank)
        (r >= 0 && (patterned || r <= k)) || throw(ArgumentError(
            "rank must lie in 0:$(k) for a level with $(k) random effects, got $(r)"))
        length(scale) == k ||
            throw(DimensionMismatch("one sdscale per random effect is required"))
        if patterned
            (length(lrow) == length(lcol) == length(load)) || throw(DimensionMismatch(
                "load_row and load_col give one effect and one dimension per loading"))
            (all(p -> 1 <= p <= k, lrow) && all(q -> 1 <= q <= r, lcol)) ||
                throw(ArgumentError("a loading's effect lies in 1:$(k) and its dimension in 1:$(r)"))
            allunique(zip(lrow, lcol)) || throw(ArgumentError(
                "each effect loads on a dimension at most once"))
            all(q -> q in lcol, 1:r) || throw(ArgumentError(
                "every dimension of a patterned level needs a loading"))
            (isempty(sd) && isempty(cor)) || throw(ArgumentError(
                "a reduced-rank level is described by load_index alone; sd_index and cor_index must be empty"))
        elseif r < k
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
        uniq = Vector{Int}(collect(unique_index))
        if !isempty(uniq)
            length(uniq) == k || throw(DimensionMismatch(
                "unique_index holds one entry per random effect, 0 for none"))
            all(>=(0), uniq) || throw(ArgumentError("unique_index entries are raw positions or 0"))
            any(>(0), uniq) && r == k && !patterned && throw(ArgumentError(
                "own variances need a reduced-rank level: on a full-rank level they are not identified"))
            all(iszero, uniq) && (uniq = Int[])
        end
        r == 0 && isempty(uniq) && k > 0 && throw(ArgumentError(
            "a rank-0 level needs an own variance on some effect; otherwise it has no variation"))
        hyb = Vector{Bool}(collect(hybrid))
        if !isempty(hyb)
            length(hyb) == k || throw(DimensionMismatch(
                "hybrid holds one entry per random effect"))
            all(p -> !hyb[p] || (!isempty(uniq) && uniq[p] > 0), 1:k) ||
                throw(ArgumentError("an effect the filter integrates needs an own variance"))
            any(hyb) || (hyb = Bool[])
        end
        isempty(grp) || (minimum(grp) >= 1 && maximum(grp) <= ngroups) ||
            throw(ArgumentError("group ids must lie in 1:ngroups"))
        sdfix = isempty(sd_fixed) ? fill(NaN, length(sd)) : Vector{Float64}(collect(sd_fixed))
        corfix = isempty(cor_fixed) ? fill(NaN, length(cor)) : Vector{Float64}(collect(cor_fixed))
        for (what, index, fixed) in (("scale", sd, sdfix), ("correlation", cor, corfix))
            length(fixed) == length(index) || throw(DimensionMismatch(
                "one fixed $(what) value, or NaN, per $(what) entry is required"))
            all((index .== 0) .== isfinite.(fixed)) || throw(ArgumentError(
                "a stated $(what) has index 0 and a finite value; a free one a raw position and NaN"))
        end
        return new(re, sd, cor, load, scale, grp, Int(ngroups), Int(covmatcode), r, uniq,
            Int[], sdfix, corfix, hyb, lrow, lcol)
    end
end

"""True when this level's population covariance is a loading matrix: it spans
fewer dimensions than it has random effects, or it states a pattern."""
isreducedrank(level::CTSEMLaplaceLevel) =
    level.rank < length(level.re_index) || !isempty(level.load_row)

"""
    _laplace_loading_entries(level)

The effect and dimension each of a reduced-rank level's loadings fills, in
`load_index` order: the level's pattern (`load_row`, `load_col`), or the
triangle a level without one has, dimension `q` loaded by effects `q:k`,
dimension by dimension.
"""
function _laplace_loading_entries(level::CTSEMLaplaceLevel)
    isempty(level.load_row) || return (level.load_row, level.load_col)
    k = nrandomeffects(level)
    rows = Int[]
    cols = Int[]
    for q in 1:level.rank, p in q:k
        push!(rows, p)
        push!(cols, q)
    end
    return (rows, cols)
end

"""
    _laplace_loading_allowed(level)

Which entries `A[m, q]` a move `L -> L A` of a reduced-rank level's loadings
may use and keep its pattern: `m >= q`, and every effect loading on dimension
`m` loads on `q` too. The triangle allows its whole lower triangle; a pattern
of blocks, each block's. The allowed pairs are closed under products, so
`exp(E)` of a matrix on them stays on them.
"""
function _laplace_loading_allowed(level::CTSEMLaplaceLevel)
    r = level.rank
    rows, cols = _laplace_loading_entries(level)
    has = falses(nrandomeffects(level), r)
    for t in eachindex(rows)
        has[rows[t], cols[t]] = true
    end
    allowed = falses(r, r)
    for q in 1:r, m in q:r
        allowed[m, q] = all(p -> !has[p, m] || has[p, q], axes(has, 1))
    end
    return allowed
end

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
nlatent(level::CTSEMLaplaceLevel) = level.rank + nunique(level)

"""
Number of this level's effects with an own variance (`unique_index`) that the
level itself carries: one column of its root each. An own variance the filter
integrates (`hybrid`) is not counted -- it has no column here.
"""
nunique(level::CTSEMLaplaceLevel) = count(p -> _laplace_own(level, p),
    eachindex(level.unique_index))

"""Whether effect `p`'s own variance is a column of the level's root."""
@inline _laplace_own(level::CTSEMLaplaceLevel, p::Integer) =
    level.unique_index[p] > 0 && (isempty(level.hybrid) || !level.hybrid[p])

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
    # The coordinates of `u` that can reach unit `U`'s data, and the rest
    # (`inactive`), both sorted. A coordinate is inactive when every effect its
    # column of the root loads on is inert for the block's group -- it moves
    # only cells of manifests the group never observes (`inert` in
    # `_laplace_build_units`). Nothing in the likelihood depends on it, so its
    # curvature row is exactly the prior's, its mode exactly 0, and every sweep
    # along it exactly zero: `_laplace_unit_curvature` and
    # `_laplace_seeded_unit_gradient!` skip it, and nothing else changes --
    # the layout of `u`, the modes and what is reported stay as they were.
    active::Vector{Vector{Int}}
    inactive::Vector{Vector{Int}}
end

"""Units with every coordinate active: the layout before inert effects were
known, and what a spec without them builds."""
CTSEMLaplaceUnits(members, offsets, dims, blocks) = CTSEMLaplaceUnits(members,
    offsets, dims, blocks, [collect(1:d) for d in dims], [Int[] for _ in dims])

"""
    _laplace_coordinate_inert(level, q, inert)

Whether coordinate `q` of `level`'s block loads only on effects in `inert`
(sorted effect indices). Read from the root's structural pattern in its natural
order: a loading column, or a full-rank Cholesky column, has rows `q:k`; an own
variance's column has its one effect's row.
"""
function _laplace_coordinate_inert(level::CTSEMLaplaceLevel, q::Integer,
    inert::AbstractVector{Int})
    isempty(inert) && return false
    k = nrandomeffects(level)
    if isreducedrank(level) && q > level.rank
        t = q - level.rank
        seen = 0
        for p in eachindex(level.unique_index)
            _laplace_own(level, p) || continue
            seen += 1
            seen == t && return insorted(p, inert)
        end
        return false
    end
    # A pattern's dimension is loaded by the effects it names, not by `q:k`.
    if !isempty(level.load_row)
        for t in eachindex(level.load_col)
            level.load_col[t] == q || continue
            insorted(level.load_row[t], inert) || return false
        end
        return true
    end
    for p in q:k
        insorted(p, inert) || return false
    end
    return true
end

"""
    _laplace_build_units(spec, nsubjects)

Group subjects into integration units and lay out each unit's latent vector.

Strict nesting is assumed and checked by the R side; here the outermost level's
grouping simply defines the units, and every inner level's groups are placed
inside whichever unit its subjects fall in.

`inert`, when given, holds per level and per group the effects (sorted indices
into the level's `re_index`) that cannot reach that group's data, which the R
side reads off the model and the data: an effect whose every cell is a
measurement cell of a manifest none of the group's subjects observes. It
decides which coordinates are inactive (`CTSEMLaplaceUnits`).
"""
function _laplace_build_units(spec::CTSEMLaplaceSpec, nsubjects::Integer;
    inert::Union{Nothing,Vector{Vector{Vector{Int}}}}=nothing)
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
    inactive = [Int[] for _ in 1:nunits]
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
                    if inert !== nothing
                        groupinert = inert[l][levels[l].group[i]]
                        for q in 1:k
                            _laplace_coordinate_inert(levels[l], q, groupinert) &&
                                push!(inactive[U], slot + q)
                        end
                    end
                end
                offsets[U][m][l] = slot
            end
        end
        dims[U] = cursor
        sort!(inactive[U])
    end
    active = [setdiff(1:dims[U], inactive[U]) for U in 1:nunits]

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
    return CTSEMLaplaceUnits(members, offsets, dims, blocks, active, inactive)
end

"""
Whether units skip their inactive coordinates (`CTSEMLaplaceUnits`). On; off
evaluates every coordinate, which gives the same objective and gradient at more
cost and is there to show that it does. `ctsem_set_laplace_prune!` sets it.
"""
const _LAPLACE_PRUNE = Ref(true)

export ctsem_set_laplace_prune!
"""Set whether Laplace units skip coordinates that cannot reach their data."""
ctsem_set_laplace_prune!(on::Bool) = (_LAPLACE_PRUNE[] = on)

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
    # Whether every unit's inner problem is a concave quadratic: Gaussian
    # indicators, a linear model, and every effect entering a mean vector
    # through the affine transform (R's `.ctJuliaLaplaceQuadratic()`). Then it
    # has one maximum, Newton from the origin finds it, and the mode probes,
    # which look for another, are skipped (`_laplace_solve_unit_mode!`).
    quadratic::Bool
end

# `gate_lo = 0.2` and `gate_hi = 0.7` are measured, not derived: see the
# provenance paragraph at the gated rule's section further down.
function CTSEMLaplaceObjective(objective::CTSEMObjective, spec::CTSEMLaplaceSpec;
    inner_maxiter::Integer=_LAPLACE_INNER_MAXITER[], inner_tol::Real=1e-10,
    floor=:total, gate_lo::Real=0.2, gate_hi::Real=0.7,
    inert::Union{Nothing,Vector{Vector{Vector{Int}}}}=nothing, quadratic::Bool=false)
    Symbol(floor) in (:total, :gated) ||
        throw(ArgumentError("the Laplace floor must be :total or :gated"))
    0 < gate_lo < gate_hi || throw(ArgumentError("the gate needs 0 < lo < hi"))
    nsubjects = length(objective.subject_objectives)
    units = _laplace_build_units(spec, nsubjects; inert=inert)
    nunits = length(units.members)
    return CTSEMLaplaceObjective{typeof(objective)}(objective, spec, units,
        [zeros(Float64, units.dims[U]) for U in 1:nunits],
        Int(inner_maxiter), Float64(inner_tol), [Dict{Any,Any}()],
        zeros(Int, nunits), zeros(Float64, nunits), falses(nunits), falses(nunits),
        falses(nunits), falses(nunits), Float64[], Symbol(floor), Float64(gate_lo),
        Float64(gate_hi), 0, quadratic)
end

"""
    _laplace_active_columns(laplace, U)

Unit `U`'s active coordinates, or `nothing` when all of them are to be
evaluated: none is inactive, pruning is off, or a level's root is permuted for
a sampling run (`_laplace_order_roots!`), whose structural pattern the inactive
set was not computed for.
"""
@inline function _laplace_active_columns(laplace::CTSEMLaplaceObjective, U::Integer)
    _LAPLACE_PRUNE[] || return nothing
    isempty(laplace.units.inactive[U]) && return nothing
    any(level -> !isempty(level.order), laplace.spec.levels) && return nothing
    return laplace.units.active[U]
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

`sd_fixed` and `cor_fixed`, when sent, run parallel to `sd_index` and
`cor_index`: the raw value of a scale or correlation coordinate the model
states, at an index of 0, and `NaN` beside a free one.

`inert_level`, `inert_group` and `inert_effect`, when sent, are parallel: each
triple says that effect `inert_effect` (an index into that level's effects)
cannot reach group `inert_group`'s data at level `inert_level`
(`_laplace_build_units`).

`hybrid_effect`, `hybrid_kind`, `hybrid_target` and `hybrid_slope`, when sent,
are parallel too: the subject level's effects (indices into its effects) whose
own variance the subject filter integrates instead of this objective, the
kind of cell each one moves (1 a manifest intercept, 2 an initial state mean,
3 a continuous intercept), the manifest or latent it moves, and that cell's
transform's slope in the raw value (`CTSEMHybrid`). They leave the level's
root and its block, and fill the wrapped objective's `params.hybrid`, which
every subject's filter shares.

`quadratic`, when sent, says every unit's inner problem is a concave quadratic,
so the mode probes are skipped (`CTSEMLaplaceObjective`'s field of that name).
"""
function ctsem_laplace_objective(objective::CTSEMObjective; re_index=Int[],
    sd_index=Int[], cor_index=Int[], sd_scale=Float64[], level_nre=Int[],
    group=Int[], level_ngroups=Int[], level_covmatcode=Int[],
    level_rank=Int[], load_index=Int[], sd_fixed=Float64[], cor_fixed=Float64[],
    unique_index=Int[], inert_level=Int[], inert_group=Int[], inert_effect=Int[],
    hybrid_effect=Int[], hybrid_kind=Int[], hybrid_target=Int[], hybrid_slope=Float64[],
    level_nload=Int[], load_row=Int[], load_col=Int[],
    inner_maxiter::Integer=_LAPLACE_INNER_MAXITER[], inner_tol::Real=1e-10,
    floor="total", gate_lo::Real=0.2, gate_hi::Real=0.7, quadratic::Bool=false)
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
    # Loading patterns (`load_row`, `load_col`), when sent, cover every reduced
    # level, with `level_nload` saying how many loadings each has.
    patterns = !isempty(load_row)
    patterns && length(level_nload) != length(counts) && throw(DimensionMismatch(
        "load_row needs level_nload, one loading count per level"))
    for l in eachindex(counts)
        k = counts[l]
        r = ranks[l]
        patterned = patterns && level_nload[l] > 0
        reduced = r < k || patterned
        ncor = reduced ? 0 : div(k * (k - 1), 2)
        nsd = reduced ? 0 : k
        nload = patterns ? Int(level_nload[l]) :
            reduced ? k * r - div(r * (r - 1), 2) : 0
        push!(levels, CTSEMLaplaceLevel(
            Int.(re_index[(re_at + 1):(re_at + k)]),
            Int.(sd_index[(sd_at + 1):(sd_at + nsd)]),
            Int.(cor_index[(cor_at + 1):(cor_at + ncor)]),
            Float64.(sd_scale[(re_at + 1):(re_at + k)]),
            groups[((l - 1) * nsubjects + 1):(l * nsubjects)],
            ngroups[l]; covmatcode=codes[l], rank=r,
            load_index=Int.(load_index[(load_at + 1):(load_at + nload)]),
            sd_fixed=isempty(sd_fixed) ? Float64[] :
                Float64.(sd_fixed[(sd_at + 1):(sd_at + nsd)]),
            cor_fixed=isempty(cor_fixed) ? Float64[] :
                Float64.(cor_fixed[(cor_at + 1):(cor_at + ncor)]),
            unique_index=isempty(unique_index) ? Int[] :
                Int.(unique_index[(re_at + 1):(re_at + k)]),
            hybrid=(l == 1 && !isempty(hybrid_effect)) ?
                [p in hybrid_effect for p in 1:k] : Bool[],
            load_row=patterned ? Int.(load_row[(load_at + 1):(load_at + nload)]) : Int[],
            load_col=patterned ? Int.(load_col[(load_at + 1):(load_at + nload)]) : Int[]))
        re_at += k; sd_at += nsd; cor_at += ncor; load_at += nload
    end
    # The deviations the filter integrates, into the shared `EKFParameters`.
    hyb = objective.params.hybrid
    hyb.block = false
    if isempty(hybrid_effect)
        hyb.kind = Int[]; hyb.target = Int[]; hyb.psi_index = Int[]; hyb.effect = Int[]
        hyb.sd_scale = Float64[]; hyb.slope = Float64[]
    else
        (length(hybrid_target) == length(hybrid_effect) == length(hybrid_slope) ==
            length(hybrid_kind)) || throw(DimensionMismatch(
                "hybrid_effect, hybrid_kind, hybrid_target and hybrid_slope must be the same length"))
        all(k -> k in (1, 2, 3), hybrid_kind) || throw(ArgumentError(
            "hybrid_kind must be 1 (a manifest intercept), 2 (an initial state mean) or 3 (a continuous intercept)"))
        lv = levels[1]
        hyb.kind = Int.(hybrid_kind)
        hyb.target = Int.(hybrid_target)
        hyb.psi_index = [lv.unique_index[p] for p in hybrid_effect]
        hyb.effect = Int.(hybrid_effect)
        hyb.sd_scale = [lv.sd_scale[p] for p in hybrid_effect]
        hyb.slope = Float64.(hybrid_slope)
    end
    inert = nothing
    if !isempty(inert_level)
        (length(inert_group) == length(inert_level) == length(inert_effect)) ||
            throw(DimensionMismatch(
                "inert_level, inert_group and inert_effect must be the same length"))
        inert = [[Int[] for _ in 1:ngroups[l]] for l in eachindex(counts)]
        for (l, g, p) in zip(inert_level, inert_group, inert_effect)
            (1 <= l <= length(counts) && 1 <= g <= ngroups[l] && 1 <= p <= counts[l]) ||
                throw(ArgumentError("inert entry (level $l, group $g, effect $p) " *
                    "is outside the hierarchy"))
            push!(inert[l][g], Int(p))
        end
        foreach(lv -> foreach(x -> unique!(sort!(x)), lv), inert)
    end
    return CTSEMLaplaceObjective(objective, CTSEMLaplaceSpec(levels);
        inner_maxiter=inner_maxiter, inner_tol=inner_tol, floor=Symbol(floor),
        gate_lo=gate_lo, gate_hi=gate_hi, inert=inert, quadratic=quadratic)
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
                              ("loading", level.load_index),
                              ("own variance", filter(>(0), level.unique_index)))
            isempty(index) && continue
            maximum(index) <= npar || throw(ArgumentError(string(
                "level ", l, " references ", what, " ", maximum(index),
                " but the parameter vector holds ", npar)))
        end
    end
    return nothing
end

"""
    _laplace_extend(laplace, theta)

`theta` as everything below an entry point reads it: the free parameters, then
any stated population values (`CTSEMAnchors`), which the random effects shift
like any other raw parameter (a level's `re_index` names them). Every entry
point that takes a parameter vector from R, the optimiser or a sampler passes
it through here first, and what it returns per parameter -- a gradient, a
Hessian, a Jacobian's columns, scores -- goes back at the caller's length. A
vector that already carries the values, or a model with none, comes back as it
is, so an entry point reached from another passes through unchanged.
"""
function _laplace_extend(laplace::CTSEMLaplaceObjective, theta::AbstractVector{T}) where {T}
    sp = laplace.objective.params
    n = length(theta)
    na = _ctsem_anchor_append(sp, n)
    na == 0 && return theta
    out = Vector{T}(undef, n + na)
    copyto!(out, 1, theta, 1, n)
    return _ctsem_fill_anchors!(out, sp, n, na)
end

"""Every row of `draws` extended as `_laplace_extend` extends one vector."""
function _laplace_extend_rows(laplace::CTSEMLaplaceObjective,
    draws::AbstractMatrix{T}) where {T}
    sp = laplace.objective.params
    n = size(draws, 2)
    na = _ctsem_anchor_append(sp, n)
    na == 0 && return draws
    out = Matrix{T}(undef, size(draws, 1), n + na)
    out[:, 1:n] .= draws
    for j in 1:na
        out[:, n + j] .= sp.anchors.values[j]
    end
    return out
end

"""How many stated population values the model has (`CTSEMAnchors`); they
are the last entries of an extended vector."""
_laplace_nanchors(laplace::CTSEMLaplaceObjective) =
    length(laplace.objective.params.anchors.values)

"""An evaluation's result with its gradient cut back to the caller's `n`
parameters (`_laplace_extend`)."""
function _laplace_free_result(r, n::Integer)
    (r === nothing || r.gradient === nothing || length(r.gradient) == n) && return r
    return merge(r, (gradient=r.gradient[1:n],))
end

"""
    _laplace_popchol(values, level)

The Cholesky factor of one level's raw-scale population covariance.

Term for term the generated Stan model's construction (`ctModelWriter.R`
around `rawpopsd = ...` through `rawpopcovchol = cholesky_decompose(...)`),
including both of its numerical offsets and its exact correlation
parameterisation, so that a Laplace fit and a Stan fit of the same model are
describing the same population distribution rather than two similar ones.

`_laplace_popcov_full` builds the covariance. Its `base` is exactly Stan's
`rawpopcovbase`: the transformed scales on the diagonal, the correlation
coordinates in the lower triangle. Handing it
to `sdcovsqrt2cov` is what makes this the same construction the filter and the
stan program use rather than a second implementation of it -- this function
used to spell the code=0 route out by hand, and so silently ignored the model's
`covmattransform`.
"""
function _laplace_popchol(values::AbstractVector{T}, level::CTSEMLaplaceLevel) where {T}
    k = nrandomeffects(level)
    k == 0 && return zeros(T, 0, 0)
    isreducedrank(level) && return _laplace_poploading(values, level)
    Sigma = _laplace_popcov_full(values, level)
    order = level.order
    if length(order) == k
        Lp = cholesky(Symmetric(Sigma[order, order])).L
        S = similar(Sigma)
        S[order, :] = Lp
        return S
    end
    return Matrix(cholesky(Symmetric(Sigma, :L)).L)
end

"""
    _laplace_popcov_full(values, level)

A full-rank level's raw-scale population covariance, by the model's own
construction and without a factor: what `_laplace_popchol` factors, and what a
report reads. A report must not factor it -- a draw from an estimate's
uncertainty can fall far enough past the correlation bound that the matrix,
valid but near singular, no longer factors in floating point.
"""
_laplace_popcov_full(values::AbstractVector, level::CTSEMLaplaceLevel) =
    _level_popcov(values, level.sd_index, level.cor_index, level.sd_fixed,
        level.cor_fixed, level.sd_scale, level.covmatcode)

"""
    _level_popcov(values, sd_index, cor_index, sd_fixed, cor_fixed, sd_scale, covmatcode)

The construction itself, from a level's arrays: `_laplace_popcov_full` for a
level, and the subject filter's for the effects it integrates on one
(`CTSEMHybrid`'s `block`), so the two read one definition.
"""
function _level_popcov(values::AbstractVector{T}, sd_index::AbstractVector{<:Integer},
    cor_index::AbstractVector{<:Integer}, sd_fixed::AbstractVector,
    cor_fixed::AbstractVector, sd_scale::AbstractVector, covmatcode::Integer) where {T}
    k = length(sd_index)
    k == 0 && return zeros(T, 0, 0)
    scales = Vector{T}(undef, k)
    @inbounds for j in 1:k
        idx = sd_index[j]
        raw = idx == 0 ? sd_fixed[j] : values[idx]
        scales[j] = log1p_exp(2 * raw - 1) * sd_scale[j] + 1e-10
    end
    base = zeros(T, k, k)
    counter = 0
    @inbounds for j in 1:k
        base[j, j] = scales[j]
        for i in 1:k
            if i > j
                counter += 1
                # No squash here: `constraincorsqrt1` applies it, and `'z'`
                # has none. How nearly singular the matrix may become is the
                # objective's bound (`_ctsem_corbound`), the same on every
                # route; a clamp on the coordinates here once made this route
                # fit a different parameter space from the augmented one.
                idx = cor_index[counter]
                base[i, j] = idx == 0 ? cor_fixed[counter] : values[idx]
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
    _sdcovsqrt2cov_uncached!(buffer, base, Int(covmatcode), k)
    return Matrix(Symmetric(buffer.out, :L))
end

# Which square root of a full-rank level's covariance standardises its effects,
# `d = S u`, is a choice: every root gives the same population distribution, the
# same Laplace objective (det(I + S'HS) = det(I + H Sigma)) and the same joint
# posterior over the deviations, with no Jacobian between them. What it decides
# is which deviations move when a sampler moves theta with `u` held. A Cholesky
# writes each correlation into the later effect's row, so the order matters:
# with a strongly determined effect after a weakly determined one, moving their
# correlation shifts the strong effect's deviations, which the data pin -- on
# bigre (CINT effects 0.9 determined, DRIFT 0.1-0.25) every chain kept its own
# correlations, R-hat about 2. Strongest first puts each mixed correlation in
# the weaker effect's row; there, R-hat 1.11 at the same draws and three times
# the effective draws a second. The equivariant root of covmattransform 'z',
# `D rownormalise(exp(A/2))`, was tried as well and was no better than the
# Cholesky on the 20-model check (review/HANDOFF-sampler-checks-2026-10-03.md).

"""
    _laplace_determined_shares(laplace, theta)

Per full-rank level, the share of each effect's variation the data determine,
averaged over the level's units: `1 - Var(d_p | data) / Var(d_p)` on the
natural deviations, from each unit's curvature at its mode -- the quantity the
identifiability report's `determined` column gives. On the deviations, so it
does not depend on the root the curvature was formed in. Reduced-rank levels
are left at zero.
"""
function _laplace_determined_shares(laplace::CTSEMLaplaceObjective,
    theta::AbstractVector)
    theta = _laplace_extend(laplace, collect(Float64, theta))
    spec = laplace.spec
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, spec)
    for U in eachindex(laplace.units.members)
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    end
    totals = [zeros(nrandomeffects(lv)) for lv in spec.levels]
    counts = [zeros(Int, nrandomeffects(lv)) for lv in spec.levels]
    for U in eachindex(laplace.units.members)
        blocks = laplace.units.blocks[U]
        isempty(blocks) && continue
        Cdiag = try
            M = _laplace_unit_curvature(laplace, U, theta, Ls, laplace.modes[U])
            fac = _laplace_factor_repaired!(M, blocks)
            fac.ok ? first(_laplace_selected_inverse(fac.factors, fac.coupling,
                blocks)) : nothing
        catch err
            _ctsem_must_propagate(err) && rethrow()
            nothing
        end
        Cdiag === nothing && continue
        for (b, blk) in enumerate(blocks)
            l = blk.level
            isreducedrank(spec.levels[l]) && continue
            L = Ls[l]
            posterior = L * Matrix(Cdiag[b]) * transpose(L)
            for p in axes(L, 1)
                prior = sum(abs2, view(L, p, :))
                (prior > 0 && isfinite(posterior[p, p])) || continue
                totals[l][p] += 1 - posterior[p, p] / prior
                counts[l][p] += 1
            end
        end
    end
    return [totals[l] ./ max.(counts[l], 1) for l in eachindex(spec.levels)]
end

"""
    _laplace_order_roots!(laplace, theta)

Order each full-rank level's effects strongest first by
`_laplace_determined_shares` at `theta`, the order `_laplace_popchol` then takes
its Cholesky root in. The stored modes are reset, since their coordinates
change. Returns the orders.
"""
function _laplace_order_roots!(laplace::CTSEMLaplaceObjective, theta::AbstractVector)
    _laplace_clear_root_orders!(laplace)
    shares = _laplace_determined_shares(laplace, theta)
    for (l, level) in enumerate(laplace.spec.levels)
        k = nrandomeffects(level)
        (k < 2 || isreducedrank(level)) && continue
        order = sortperm(shares[l]; rev=true)
        resize!(level.order, k)
        copyto!(level.order, order)
    end
    for U in eachindex(laplace.modes)
        fill!(laplace.modes[U], 0.0)
    end
    return [copy(level.order) for level in laplace.spec.levels]
end

"""Every level back to the natural order, and the stored modes reset."""
function _laplace_clear_root_orders!(laplace::CTSEMLaplaceObjective)
    for level in laplace.spec.levels
        empty!(level.order)
    end
    for U in eachindex(laplace.modes)
        fill!(laplace.modes[U], 0.0)
    end
    return laplace
end

"""
    _with_root_orders(f, laplace, theta, natural)

`f()` with each level's root ordered by `_laplace_order_roots!` for the call,
the natural order restored afterwards whatever happens. `natural = true` -- a
caller that supplied its own starts or SAEM state, expressed in the natural
order -- runs `f()` unordered.
"""
function _with_root_orders(f, laplace::CTSEMLaplaceObjective, theta::AbstractVector,
    natural::Bool)
    natural && return f()
    _laplace_order_roots!(laplace, theta)
    try
        return f()
    finally
        _laplace_clear_root_orders!(laplace)
    end
end

"""
    _laplace_poploading(values, level)

A reduced-rank level's factor: `k x nlatent`, the loading itself followed by
one column per own variance (`unique_index`), so that `F * F'` is the
population covariance and `F * u` is the deviation for `u ~ N(0, I)` of
dimension `nlatent` (the rank plus the own variances).

The block is `nlatent`-dimensional, not `k`, which is what makes a reduced
level cheap; see `nlatent`. No Cholesky is taken, because there is nothing
to decompose -- the parameters *are* the factor, which is also why the loadings
need no positivity transform and no correlation cap (an own variance takes the
scales' transform, which also fixes the sign of its single-entry column). `sd_scale` multiplies a row, so a
level still scales its own spread exactly as the full-rank form does, and a
loading is signed: negating a whole column leaves `L * L'` alone, so read the
covariance rather than the sign of one entry.
"""
function _laplace_poploading(values::AbstractVector{T}, level::CTSEMLaplaceLevel) where {T}
    k = nrandomeffects(level)
    r = level.rank
    L = zeros(T, k, nlatent(level))
    if isempty(level.load_row)
        counter = 0
        @inbounds for q in 1:r
            for p in q:k
                counter += 1
                L[p, q] = values[level.load_index[counter]] * level.sd_scale[p]
            end
        end
    else
        # A pattern (`load_row`, `load_col`): each loading where it says.
        @inbounds for t in eachindex(level.load_row)
            p = level.load_row[t]
            L[p, level.load_col[t]] = values[level.load_index[t]] * level.sd_scale[p]
        end
    end
    # Own variances (`unique_index`): one column each, a single entry in the
    # effect's own row, transformed as a full-rank level's scale
    # (`_laplace_popchol`) so a raw value means the same spread either way.
    column = r
    @inbounds for p in eachindex(level.unique_index)
        _laplace_own(level, p) || continue
        idx = level.unique_index[p]
        column += 1
        L[p, column] = log1p_exp(2 * values[idx] - 1) * level.sd_scale[p] + 1e-10
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
    return _laplace_level_covariance(laplace.spec.levels[level],
        collect(Float64, values))
end

"""
    _laplace_level_covariance(level, x)

A level's raw-scale population covariance at `x`: its root's `L L'`, plus the
squared own sd of every effect whose deviation the subject filter integrates
(`hybrid`), which has no column of the root.
"""
function _laplace_level_covariance(lv::CTSEMLaplaceLevel, x::Vector{Float64})
    Sigma = if isreducedrank(lv)
        L = _laplace_poploading(x, lv)
        L * transpose(L)
    else
        _laplace_popcov_full(x, lv)
    end
    for p in eachindex(lv.hybrid)
        lv.hybrid[p] || continue
        psi = log1p_exp(2 * x[lv.unique_index[p]] - 1) * lv.sd_scale[p] + 1e-10
        Sigma[p, p] += psi^2
    end
    return Sigma
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
(`log1p_exp(2x-1) * sdscale + 1e-10`, and for the correlation square root
`2/(1+exp(-x)) - 1`; `'z'` has no squash to saturate), which sit outside
`regular_transforms` entirely -- "a tail of the vector no model matrix cell
reads", per `CTSEMLaplaceLevel`'s docstring -- so they need their own
derivative here.

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
            idx == 0 && continue
            scale = level.sd_scale[j]
            d = abs(ForwardDiff.derivative(
                raw -> log1p_exp(2 * raw - 1) * scale + 1e-10, values[idx]))
            d < threshold && push!(found, idx)
        end
        level.covmatcode == 2 && continue
        for idx in level.cor_index
            idx == 0 && continue
            d = abs(ForwardDiff.derivative(
                raw -> 2 / (1 + exp(-raw)) - 1, values[idx]))
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
    # Outside the loop's `@inbounds`: a vector short of an index a level
    # shifts -- a stated population value (`CTSEMAnchors`) an entry point did
    # not append (`_laplace_extend`) -- would otherwise be written past its end.
    for level in spec.levels
        isempty(level.re_index) || checkbounds(values, maximum(level.re_index))
    end
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
    # session may never have evaluated -- `ctFitUncertainty(fit, uncertainty =
    # 'sample')` on a reloaded fit builds a fresh one. The band lives on the
    # task and the store lives on the object, so a task carrying a band from
    # an earlier object spawned workers
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
                Threads.@spawn _ctsem_with_scratch() do
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
targets. It reaches here only when the caller passed `parallel = true`, which
`ctsem_sample` does after sizing the store, so the region below is never the
outermost one. `ctsem_sample_marginal` passes `false` and never reaches here:
its chains share one objective, and run one after another.

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
        Threads.@spawn _ctsem_with_scratch() do
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
    #
    # If the caller's share throws -- R asking the call to stop is the usual
    # reason, from a chain's progress checkpoint (interrupt.jl) -- the helpers
    # are told to stop and waited for before the exception goes on. Otherwise
    # they outlive the call that spawned them: its reply reaches R, R removes
    # the file they would have stopped at, and they sample on in the background
    # alongside whatever R asks for next.
    try
        _laplace_pull!(f, items, n, next, failed)
    catch
        failed[] = true
        for helper in helpers
            try
                wait(helper)
            catch
            end
        end
        rethrow()
    end
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

"""The calling slot's `CTSEMSymEigWork` for dimension `d`."""
function _laplace_scratch_symeig!(laplace::CTSEMLaplaceObjective, d::Integer)
    _laplace_check_slot(laplace, _laplace_slot())
    store = laplace.workspaces[_laplace_slot()]
    key = (:symeig, Int(d))
    cached = get(store, key, nothing)
    cached === nothing || return cached::CTSEMSymEigWork
    built = CTSEMSymEigWork(d)
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
    # Keyed by type alone and checked against the shapes, rather than keyed by
    # `size.(Ls)`, which built a vector on every curvature sweep. An objective's
    # levels keep their shapes, so the check passes after the first build.
    key = (tag, S)
    cached = get(store, key, nothing)
    if cached === nothing || !_laplace_same_shapes(cached::Vector{Matrix{S}}, Ls)
        cached = [Matrix{S}(undef, size(L)) for L in Ls]
        store[key] = cached
    end
    out = cached::Vector{Matrix{S}}
    @inbounds for i in eachindex(Ls)
        copyto!(out[i], Ls[i])
    end
    return out
end

function _laplace_same_shapes(a::Vector, b::Vector)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        size(a[i]) == size(b[i]) || return false
    end
    return true
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
    subject_objective, aws, values::AbstractVector{T},
    ekf_workspace=nothing, mode::Int=_HYBRID_PROFILE) where {T}
    # Positional, not a keyword: the sampler calls this with `aws` untyped, and
    # a dynamic call with a keyword builds its named tuple on the heap, the
    # workspace copied inside it -- 56 KB an evaluation on gC8.
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
    # The filter's own deviations (`CTSEMHybrid`): their profile unless `mode`
    # asks for something else, and back to the profile afterwards, so a direct
    # filter call on this workspace never inherits a caller's setting.
    ws.hyb.mode = mode
    tape = _tape_reset!(aws.tape)
    # The row this evaluation reads: a sampled TI predictor cell's value is
    # the raw parameter's at `values`, and the forward pass and the group
    # replay (`ctx.tipreds`) must both see it.
    tipreds = _ctsem_tipred_fill!(aws.tipreds, subject_objective.tipreds, values)
    aws.frechet_count = 0
    deferred = aws.defer_frechet
    aws.defer_frechet = false
    try
        local bf, br
        bf = _laplace_mark()
        loglik = _extended_kalman_filter_continuous!(ws, values,
            subject_objective.data, subject_objective.timesteps,
            subject_objective.params, subject_objective.tdpreds,
            tipreds, subject_objective.subject,
            subject_objective.max_timestep, tape)
        _laplace_charge!(_LAPLACE_BYTES_FWD, bf)
        isfinite(loglik) || return loglik
        if _growth_unstable(ws)
            _ctsem_forward_subject_gradient!(gradient, subject_objective, values,
                mode, ws.hyb.dfix)
            return loglik
        end
        br = _laplace_mark()
        fill!(aws.theta_bar, zero(T))
        _ctsem_tipred_bar_reset!(aws, subject_objective.tipreds, length(tipreds))
        _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
        fill!(gradient, zero(T))
        _ctsem_parameter_layer!(gradient, aws.theta_bar, tape.subject_values,
            subject_objective.params, aws, subject_objective.tipreds, values)
        _hybrid_scatter!(gradient, tape.hyb, subject_objective.params.hybrid)
        _ctsem_tipred_scatter!(gradient, aws, subject_objective.tipreds)
        _laplace_charge!(_LAPLACE_BYTES_REV, br)
        return loglik
    finally
        aws.defer_frechet = deferred
        ws.hyb.mode = _HYBRID_PROFILE
    end
end

"""
    _laplace_member_loglik(aws, so, shifted)

A member's log likelihood at `shifted`, value only, on this task's workspace,
with the deviations its filter integrates (`CTSEMHybrid`) integrated exactly --
their profile less `logdet(M) / 2` (`_HYBRID_FULL`) -- so that it is the
integrand over the unit's effects: what quadrature, the continuation and SAEM
integrate, where the Laplace layer takes the profile. Without such deviations
it is the subject's own value. A barrier on `aws`, which arrives untyped from
the slot's cache: read from it directly, its inline `ekf_ws` was copied to the
heap for every member, 700 bytes a call on gC8.
"""
@noinline function _laplace_member_loglik(aws, so, shifted)
    hb = aws.ekf_ws.hyb
    hb.mode = _HYBRID_FULL
    try
        return _ekf_run(aws.ekf_ws, so, shifted,
            _ctsem_tipred_fill!(aws.tipreds, so.tipreds, shifted), nothing)
    finally
        hb.mode = _HYBRID_PROFILE
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

"""A block matrix's copy, to keep past the slot scratch it was formed in."""
_laplace_copy_blockmatrix(M::CTSEMBlockMatrix{T}) where {T} =
    CTSEMBlockMatrix{T}([copy(b) for b in M.diag],
        [[copy(c) for c in row] for row in M.coupling])

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

Keyed by the unit's shape (`_laplace_shape_key`): a unit of one block shares
it with every unit of that size, and any other unit holds its own. A worker
holds one per unit it is given, and the pool strides units across workers
deterministically, so the total held is the unit count rather than units times
workers.
"""
function _laplace_scratch_blockmatrix!(laplace::CTSEMLaplaceObjective, ::Type{T},
    U::Integer, blocks::Vector{CTSEMLaplaceBlock}, tag::Symbol) where {T}
    _laplace_check_slot(laplace, _laplace_slot())
    store = laplace.workspaces[_laplace_slot()]
    key = _laplace_shape_key(tag, T, U, blocks)
    cached = get(store, key, nothing)
    if cached !== nothing
        c = cached::CTSEMBlockMatrix{T}
        # A layout other than the unit's own (the gated term's joint
        # coordinates) can share a key with it; what does not fit is replaced.
        _laplace_blocks_fit(c.diag, c.coupling, blocks) && return c
    end
    built = CTSEMBlockMatrix(T, blocks)
    store[key] = built
    return built
end

"""Reassemble a block matrix densely. Testing and small-unit use only."""
_laplace_block_dense(M::CTSEMBlockMatrix{T}, blocks::Vector{CTSEMLaplaceBlock},
    dim::Integer) where {T} = _laplace_block_dense!(zeros(T, dim, dim), M, blocks)

"""`_laplace_block_dense` into `out`, which is overwritten."""
function _laplace_block_dense!(out::AbstractMatrix, M::CTSEMBlockMatrix,
    blocks::Vector{CTSEMLaplaceBlock})
    fill!(out, zero(eltype(out)))
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
Storage for a factorization nobody keeps -- a Newton step's, used for one solve
and dropped: the per-block copies `_laplace_block_factor` factors in place, the
eliminated couplings, and the factors. One per slot and unit shape
(`_laplace_factor_store!`); the primal's factorization, which the gradient
phase reads, never uses one.
"""
struct _LaplaceFactorStore{T}
    diag::Vector{Matrix{T}}
    coupling::Vector{Vector{Matrix{T}}}
    factors::Vector{_LaplaceCholesky{T}}
end

function _LaplaceFactorStore{T}(blocks::Vector{CTSEMLaplaceBlock}) where {T}
    M = CTSEMBlockMatrix(T, blocks)
    return _LaplaceFactorStore{T}(M.diag, M.coupling,
        Vector{_LaplaceCholesky{T}}(undef, length(blocks)))
end

"""
A per-slot scratch key for something shaped like unit `U`: by size alone for a
unit of one block with nothing above it -- every unit of a single-level model --
so all units of that size share one, and by unit otherwise.
"""
_laplace_shape_key(tag::Symbol, ::Type{T}, U::Integer,
    blocks::Vector{CTSEMLaplaceBlock}) where {T} =
    length(blocks) == 1 && isempty(blocks[1].ancestors) ?
        (tag, T, :single, blocks[1].size) : (tag, T, :unit, Int(U))

# The working arrays `_laplace_block_factor` factors into: copies of `M`'s
# blocks, fresh or the store's. `factors` is concretely typed, not `Vector{Any}`:
# `factors[b] \ x` and `logdet(factors[b])` are called once per block per solve,
# and through an `Any` element every one of them is a dynamic dispatch.
function _laplace_factor_buffers(M::CTSEMBlockMatrix{T}, nb::Int, ::Nothing) where {T}
    return ([copy(d) for d in M.diag], [[copy(c) for c in row] for row in M.coupling],
        Vector{_LaplaceCholesky{T}}(undef, nb))
end

function _laplace_factor_buffers(M::CTSEMBlockMatrix{T}, nb::Int,
    store::_LaplaceFactorStore{T}) where {T}
    for b in 1:nb
        copyto!(store.diag[b], M.diag[b])
        for t in eachindex(M.coupling[b])
            copyto!(store.coupling[b][t], M.coupling[b][t])
        end
    end
    return (store.diag, store.coupling, store.factors)
end

"""
The calling slot's `_LaplaceFactorStore` for `blocks` -- unit `U`'s own unless
the caller factors another layout of it, as the gated term does in joint
coordinates (`_hybrid_joint_gated`), whose one block also spans the deviations
the filter integrates. A cached store that does not fit `blocks` is replaced.
"""
function _laplace_factor_store!(laplace::CTSEMLaplaceObjective, ::Type{T},
    U::Integer, blocks::Vector{CTSEMLaplaceBlock}=laplace.units.blocks[U]) where {T}
    _laplace_check_slot(laplace, _laplace_slot())
    store = laplace.workspaces[_laplace_slot()]
    key = _laplace_shape_key(:factor_store, T, U, blocks)
    cached = get(store, key, nothing)
    if cached !== nothing
        c = cached::_LaplaceFactorStore{T}
        _laplace_blocks_fit(c.diag, c.coupling, blocks) && return c
    end
    built = _LaplaceFactorStore{T}(blocks)
    store[key] = built
    return built
end

"""Whether block storage `diag`, `coupling` is shaped for `blocks`."""
function _laplace_blocks_fit(diag::Vector, coupling::Vector,
    blocks::Vector{CTSEMLaplaceBlock})
    length(diag) == length(blocks) == length(coupling) || return false
    @inbounds for (b, block) in enumerate(blocks)
        size(diag[b]) == (block.size, block.size) || return false
        length(coupling[b]) == length(block.ancestors) || return false
        for (t, a) in enumerate(block.ancestors)
            size(coupling[b][t]) == (block.size, blocks[a].size) || return false
        end
    end
    return true
end

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
    blocks::Vector{CTSEMLaplaceBlock};
    store::Union{Nothing,_LaplaceFactorStore{T}}=nothing) where {T}
    nb = length(blocks)
    # Copies of the blocks to factor in place: fresh where the caller keeps the
    # factorization (the primal's, read by the gradient phase), a store's where
    # it does not -- a Newton step uses its factors once. Assigned once, from a
    # function: the comprehension below captures `coupling`, and a captured
    # name assigned in two branches is boxed.
    diag, coupling, factors = _laplace_factor_buffers(M, nb, store)
    total = zero(T)
    for b in 1:nb
        # In place on this function's own copy of the block, which nothing
        # reads once it is factored -- the factor keeps it. Symmetrising into
        # a new matrix and then copying that was two more allocations per block
        # per evaluation (gC8: 2,400 a Laplace gradient).
        f = _ctsem_cholesky(_laplace_symmetrise!(diag[b]), size(diag[b], 1))
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
        # The diagonal's extremes by a loop: `diag(f.factors)` was a vector per
        # block per factorization.
        let F = f.factors, k = size(F, 1)
            if k > 0
                lo = hi = abs(F[1, 1])
                @inbounds for i in 2:k
                    lo = min(lo, abs(F[i, i]))
                    hi = max(hi, abs(F[i, i]))
                end
                hi > 0 && lo <= sqrt(eps(real(float(one(T))))) * hi &&
                    return (false, T(NaN), factors, coupling)
            end
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
            # `_ctsem_mulTN!`, not `transpose(B) * W`: that is OpenBLAS's
            # `dgemm_tn`, which takes a process-wide lock on every call.
            update = _ctsem_mulTN!(similar(W[t], size(coupling[b][t], 2), size(W[t], 2)),
                coupling[b][t], W[t])
            diag[a] .-= _laplace_symmetrise!(update)
            for s in eachindex(ancestors)
                s == t && continue
                c = ancestors[s]
                # Where does the (a, c) coupling live? One of them is an
                # ancestor of the other; the update belongs on the lower one's
                # row.
                slot = findfirst(==(c), blocks[a].ancestors)
                if slot !== nothing
                    _ctsem_mulTN!(coupling[a][slot], coupling[b][t], W[s], -1, true)
                end
            end
        end
    end
    return (true, total, factors, coupling)
end

@inline _laplace_symmetrise(A) = (A .+ transpose(A)) ./ 2

"""`A` overwritten by `(A + A') / 2`, and returned."""
function _laplace_symmetrise!(A::AbstractMatrix)
    n = size(A, 1)
    @inbounds for j in 1:n, i in (j + 1):n
        v = (A[i, j] + A[j, i]) / 2
        A[i, j] = v
        A[j, i] = v
    end
    return A
end

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
`_laplace_selected_inverse` for the seeded gradient, which reads it and drops
it: for a unit of one block, `C = sym(inv(D))` into the slot's storage by the
same column solves `inv` makes; any other unit takes the general route.
"""
function _laplace_selected_inverse_scratch!(laplace::CTSEMLaplaceObjective,
    U::Integer, factors, elim, blocks::Vector{CTSEMLaplaceBlock})
    (length(blocks) == 1 && isempty(blocks[1].ancestors)) ||
        return _laplace_selected_inverse(factors, elim, blocks)
    k = blocks[1].size
    C = _laplace_scratch_blockmatrix!(laplace, Float64, U, blocks, :sg_C)
    Dinv = _laplace_scratch_matrix!(laplace, Float64, k, k, :sg_Dinv)
    e = fill!(_laplace_scratch_vector!(laplace, Float64, k, :sg_unit), 0.0)
    @inbounds for j in 1:k
        e[j] = 1.0
        ldiv!(view(Dinv, :, j), factors[1], e)
        e[j] = 0.0
    end
    Cd = C.diag[1]
    @inbounds for j in 1:k, i in 1:k
        Cd[i, j] = (Dinv[i, j] + Dinv[j, i]) / 2
    end
    return (C.diag, C.coupling)
end

"""
    _laplace_block_solve(factors, coupling, blocks, rhs)

Solve `M x = rhs` using the factorization above, by forward substitution over
the elimination order and back substitution against it.
"""
function _laplace_block_solve(factors, coupling, blocks::Vector{CTSEMLaplaceBlock},
    rhs::Vector{T}) where {T}
    nb = length(blocks)
    # One block with nothing above it -- every single-level unit -- is one
    # solve; the elimination below is the same solve with five temporaries
    # round it, once per Newton step.
    nb == 1 && blocks[1].size == length(rhs) && return factors[1] \ rhs
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

"""`_laplace_block_solve` into `out`; for a single block, with nothing built."""
function _laplace_block_solve!(out::AbstractVector, factors, coupling,
    blocks::Vector{CTSEMLaplaceBlock}, rhs::Vector)
    length(blocks) == 1 && blocks[1].size == length(rhs) &&
        return ldiv!(out, factors[1], rhs)
    return copyto!(out, _laplace_block_solve(factors, coupling, blocks, rhs))
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
gradient with respect to `u`, and `hlogdet`, the summed `logdet(M) / 2` of the
deviations the members' filters integrate (`CTSEMHybrid`; zero without), whose
profile the log likelihood holds -- without the `-u'u/2` term, which belongs to the
unit as a whole rather than to any member.

Restricting the member set is what makes the blocked curvature below cheap. It
is exact, not an approximation: a member outside the set has no dependence on
the block being differentiated, so its second derivative with respect to that
block is structurally zero rather than merely small.
"""
function _laplace_unit_loglik_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws, positions, into::Union{Nothing,AbstractVector}=nothing,
    mode::Int=_HYBRID_PROFILE, hook=nothing) where {T}
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
        # `hook`: what a caller reads from, or hands to, each member's filter
        # around its pass (`_hybrid_before!`, `_hybrid_after!`).
        _hybrid_before!(hook, m, ws)
        loglik = _laplace_subject_value_gradient!(grad,
            laplace.objective.subject_objectives[i], ws, shifted, nothing, mode)
        isfinite(loglik) || return loglik
        _hybrid_after!(hook, m, ws, u)
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

    # A range or vector is indexed as it is; anything else is collected.
    pos = positions isa AbstractVector{<:Integer} ? positions : collect(positions)
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
        hsum = zero(T)
        gradient = _laplace_scratch_vector!(laplace, T, length(values), :loglik_grad)
        shift = _laplace_scratch_vector!(laplace, T, length(values), :loglik_shift)
        for m in pos
            loglik = one!(m, aws, gradient, shift, inner)
            isfinite(loglik) ||
                return (value=loglik, gradient=fill(_ctsem_nan(T), length(u)), hlogdet=zero(T))
            total += loglik
            hsum += aws.ekf_ws.hyb.half_logdet
        end
        return (value=total, gradient=inner, hlogdet=hsum)
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
    hls = zeros(T, width)
    bad = zeros(T, width)
    okk = fill(true, width)
    _laplace_partition(laplace, npos) do mine, _w
        local ws, grad, shift, acc, tot, htot, i, loglik, slot
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
        htot = zero(T)
        for i in mine
            loglik = one!(pos[i], ws, grad, shift, acc)
            if !isfinite(loglik)
                okk[slot] = false
                bad[slot] = loglik
                return false
            end
            tot += loglik
            htot += ws.ekf_ws.hyb.half_logdet
        end
        totals[slot] = tot
        hls[slot] = htot
        return true
    end
    @inbounds for w in 1:width
        okk[w] || return (value=bad[w], gradient=fill(_ctsem_nan(T), length(u)), hlogdet=zero(T))
    end
    inner = into === nothing ? zeros(T, length(u)) : fill!(into, zero(T))
    @inbounds for w in 1:width
        used[w] || continue
        acc = _laplace_scratch_vector!(laplace, T, length(u), :primal_inner, w)
        for j in eachindex(inner)
            inner[j] += acc[j]
        end
    end
    return (value=sum(totals), gradient=inner, hlogdet=sum(hls))
end

"""
`_laplace_unit_loglik_gradient(...).gradient`, for the curvature's sweeps. Their
workspace is a slot's, which `_laplace_workspace!` cannot type, so the call is a
dynamic dispatch: returning the named tuple boxed it on every sweep, where the
vector alone is a reference.
"""
_laplace_unit_loglik_gradient_vector(laplace, U, values, Ls, u, aws, positions,
    into, mode, hook) = _laplace_unit_loglik_gradient(laplace, U, values, Ls, u,
    aws, positions === nothing ? eachindex(laplace.units.members[U]) : positions,
    into, mode, hook).gradient

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
    return (value=result.value - dot(u, u) / 2, gradient=inner, hlogdet=result.hlogdet)
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

"""
Whether the inner solve looks for a better mode than the one it reaches from
the origin (`_laplace_probe_modes`). On by default; `ctsem_set_mode_probes!`.
"""
const _LAPLACE_MODE_PROBES = Ref(true)

# Whether the inner solve updates its curvature by SR1 between exact ones
# (`_laplace_newton_unit_mode`); off, every Newton step forms it exactly. For
# the tests, which compare the two.
const _LAPLACE_INNER_SR1 = Ref(true)

"""
Precision below which a direction of a mode's curvature is probed: the
continuation's `soft_tau`, the precision below which a direction can be
non-Gaussian enough to matter.
"""
const _LAPLACE_PROBE_SOFT = 3.5

"""
How far a probe starts from the mode, in posterior sds along its direction.
Two missed, on the pilot's person-persistence fit, a person's second mode two
units away past a valley at one, which a probe started at 2.5 and at 3 reached;
the point a probe starts at has to be past the valley, not merely towards it.
"""
const _LAPLACE_PROBE_REACH = 3.0

"""
    ctsem_set_mode_probes!(on)

Turn the inner solve's probes for a better mode (`_laplace_probe_modes`) on or
off. Returns the previous setting. Off is the objective as it was, Newton from
the origin alone, which is what a measurement of the probes is made against.
"""
function ctsem_set_mode_probes!(on::Bool)
    previous = _LAPLACE_MODE_PROBES[]
    _LAPLACE_MODE_PROBES[] = on
    return previous
end

export ctsem_set_mode_probes!

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
    H = _ctsem_jacobian(inner_of, collect(u), _CTSEM_CURVATURE_WIDTHS[])
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
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T};
    hook=nothing, mode::Int=_HYBRID_PROFILE,
    into::Union{Nothing,CTSEMBlockMatrix{T}}=nothing) where {T}
    # `into`: where a caller that uses the curvature once (a Newton step) wants
    # it, rather than a block matrix of its own. Only the small-unit route below
    # reads it; the blocked route fills the slot's per-unit matrix anyway.
    # `hook` is handed every member's pass (`_laplace_unit_loglik_gradient`):
    # a `CTSEMHybridCross` reads the cross curvature with the deviations the
    # filter integrates off the same sweeps. `mode` is what those passes return
    # for them (`CTSEMHybridBuffers`): their profile, or at fixed deviations,
    # which with a hook is the joint curvature's every block
    # (`_hybrid_joint_curvature`).
    blocks = laplace.units.blocks[U]
    base = _laplace_scratch_vector!(laplace, T, length(u), :curv_base)
    copyto!(base, u)
    act = _laplace_active_columns(laplace, U)
    act === nothing ||
        return _laplace_unit_curvature_pruned(laplace, U, values, Ls, base, act;
            hook=hook, mode=mode)
    if length(blocks) <= 1 || length(u) < _LAPLACE_BLOCK_THRESHOLD[]
        gradient_of = function (uu)
            S = eltype(uu)
            ws = _laplace_workspace!(laplace, S, length(values))
            vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
            copyto!(vs, values)
            Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
            return _laplace_unit_loglik_gradient_vector(laplace, U, vs, Lss, uu, ws,
                nothing,
                _laplace_scratch_vector!(laplace, S, length(uu), :curv_inner),
                mode, hook)
        end
        # Bucketed width: see `_CTSEM_CURVATURE_WIDTHS`. Into the slot's
        # scratch, padded input and all, as the blocked route below does: the
        # padding is `n + 1:nx` of `uu`, which no member reads, and the rows of
        # the gradient it adds are zero and not read either. Built fresh, the
        # padded input, the Jacobian, its leading block, the identity, the
        # symmetrised copy and the dense curvature were six arrays a Newton
        # step, most of what a small unit's mode solve allocated.
        n = length(base)
        width = _ctsem_dual_width(n, _CTSEM_CURVATURE_WIDTHS[])
        nx = max(n, width)
        xb = _laplace_scratch_vector!(laplace, T, nx, :curv_ub)
        @inbounds for t in 1:nx
            xb[t] = t <= n ? base[t] : zero(T)
        end
        Jw = _laplace_scratch_matrix!(laplace, T, nx, nx, :curv_jac)
        _ctsem_width_jacobian!(Jw, gradient_of, xb, n, width; whole=true,
            store=laplace.workspaces[_laplace_slot()])
        # M = I - (J + J') / 2, block by block. A matrix of its own unless the
        # caller handed one: `_laplace_unit_primal`'s is kept until the gradient
        # reads it.
        # Not `M`: the blocked route below assigns that name and its closures
        # capture it, and a captured name assigned twice is boxed.
        Mdense = into === nothing ? CTSEMBlockMatrix(T, blocks) : into
        for (b, block) in enumerate(blocks)
            ob = block.offset
            D = Mdense.diag[b]
            @inbounds for j in 1:block.size, i in 1:block.size
                D[i, j] = (i == j ? one(T) : zero(T)) -
                    (Jw[ob + i, ob + j] + Jw[ob + j, ob + i]) / 2
            end
            for (t, a) in enumerate(block.ancestors)
                oa = blocks[a].offset
                C = Mdense.coupling[b][t]
                @inbounds for j in 1:blocks[a].size, i in 1:block.size
                    C[i, j] = -(Jw[ob + i, oa + j] + Jw[oa + j, ob + i]) / 2
                end
            end
        end
        return Mdense
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
        local block, columns, block_of, J, t, a, rows, xb, width, nx, Jw
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
            return _laplace_unit_loglik_gradient_vector(laplace, U, vs, Lss, uu, ws,
                block.members,
                _laplace_scratch_vector!(laplace, S, length(base), :curv_inner),
                mode, hook)
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
        # At a bucketed width (see `_CTSEM_CURVATURE_WIDTHS`): a block smaller
        # than its width is padded with zeros `block_of` never reads, and the
        # padding's columns of the result are ignored.
        width = _ctsem_dual_width(block.size, _CTSEM_CURVATURE_WIDTHS[])
        nx = max(block.size, width)
        xb = _laplace_scratch_vector!(laplace, T, nx, :curv_ub)
        @inbounds for (t, c) in enumerate(columns)
            xb[t] = base[c]
        end
        @inbounds for t in (block.size + 1):nx
            xb[t] = zero(T)
        end
        Jw = _laplace_scratch_matrix!(laplace, T, length(base), nx, :curv_jac)
        _ctsem_width_jacobian!(Jw, block_of, xb, block.size, width; whole=true,
            store=laplace.workspaces[_laplace_slot()])
        J = view(Jw, :, 1:block.size)
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
    _laplace_unit_curvature_pruned(laplace, U, values, Ls, base, act)

`_laplace_unit_curvature` for a unit with inactive coordinates
(`CTSEMLaplaceUnits`): the Jacobian is taken over the active coordinates `act`
alone, and every inactive row and column is the prior's identity with no
coupling -- which is what it is exactly, since no member's likelihood depends
on an inactive coordinate. The same matrix the full Jacobian gives, at one
forward lane per active coordinate rather than per coordinate. `base` holds
`u`. A separate function rather than a branch inside the other, so a unit
without inactive coordinates runs, and compiles, exactly what it did.
"""
function _laplace_unit_curvature_pruned(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, base::AbstractVector{T},
    act::Vector{Int}; hook=nothing, mode::Int=_HYBRID_PROFILE) where {T}
    blocks = laplace.units.blocks[U]
    d = length(base)
    if length(blocks) <= 1 || d < _LAPLACE_BLOCK_THRESHOLD[]
        dense = Matrix{T}(LinearAlgebra.I, d, d)
        if !isempty(act)
            gradient_of = function (ua)
                local S, ws, vs, Lss, full, g, out
                S = eltype(ua)
                ws = _laplace_workspace!(laplace, S, length(values))
                vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
                copyto!(vs, values)
                Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
                full = _laplace_scratch_vector!(laplace, S, d, :curv_full)
                @inbounds for c in 1:d
                    full[c] = base[c]
                end
                @inbounds for (t, c) in enumerate(act)
                    full[c] = ua[t]
                end
                g = _laplace_unit_loglik_gradient(laplace, U, vs, Lss, full, ws,
                    eachindex(laplace.units.members[U]),
                    _laplace_scratch_vector!(laplace, S, d, :curv_inner),
                    mode, hook).gradient
                out = _laplace_scratch_vector!(laplace, S, length(act), :curv_act)
                @inbounds for (t, c) in enumerate(act)
                    out[t] = g[c]
                end
                return out
            end
            # Bucketed width: see `_CTSEM_CURVATURE_WIDTHS`.
            A = _laplace_symmetrise(_ctsem_width_jacobian(gradient_of, base[act],
                _CTSEM_CURVATURE_WIDTHS[]))
            @inbounds for (j, cj) in enumerate(act), (i, ci) in enumerate(act)
                dense[ci, cj] -= A[i, j]
            end
        end
        return _laplace_block_of(dense, blocks)
    end
    isactive = fill(false, d)
    @inbounds for c in act
        isactive[c] = true
    end
    M = _laplace_scratch_blockmatrix!(laplace, T, U, blocks, :curv_M)
    # As in `_laplace_unit_curvature`'s block route, over each block's active
    # columns only; `local` for the reason given there.
    fill_block! = function (b::Int)
        local block, off, columns, n, D, Cp, block_of, width, nx, xb, Jw, rows
        block = blocks[b]
        off = block.offset
        columns = Int[c for c in (off + 1):(off + block.size) if isactive[c]]
        n = length(columns)
        D = M.diag[b]
        fill!(D, zero(T))
        @inbounds for t in 1:block.size
            D[t, t] = one(T)
        end
        for t in eachindex(block.ancestors)
            fill!(M.coupling[b][t], zero(T))
        end
        n == 0 && return nothing
        block_of = function (ub)
            local S, ws, vs, Lss, uu
            S = eltype(ub)
            ws = _laplace_workspace!(laplace, S, length(values))
            vs = _laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
            copyto!(vs, values)
            Lss = _laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
            uu = _laplace_scratch_vector!(laplace, S, d, :curv_uu)
            copyto!(uu, base)
            @inbounds for (t, c) in enumerate(columns)
                uu[c] = ub[t]
            end
            return _laplace_unit_loglik_gradient_vector(laplace, U, vs, Lss, uu, ws,
                block.members, _laplace_scratch_vector!(laplace, S, d, :curv_inner),
                mode, hook)
        end
        width = _ctsem_dual_width(n, _CTSEM_CURVATURE_WIDTHS[])
        nx = max(n, width)
        xb = _laplace_scratch_vector!(laplace, T, nx, :curv_ub)
        @inbounds for (t, c) in enumerate(columns)
            xb[t] = base[c]
        end
        @inbounds for t in (n + 1):nx
            xb[t] = zero(T)
        end
        Jw = _laplace_scratch_matrix!(laplace, T, d, nx, :curv_jac)
        _ctsem_width_jacobian!(Jw, block_of, xb, n, width; whole=true,
            store=laplace.workspaces[_laplace_slot()])
        # M = I - d2(sum ll)/du du on the active columns, as the full route.
        @inbounds for (j, cj) in enumerate(columns), ci in columns
            D[ci - off, cj - off] -= Jw[ci, j]
        end
        for (t, a) in enumerate(block.ancestors)
            rows = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
            Cp = M.coupling[b][t]
            @inbounds for (j, cj) in enumerate(columns), (i, r) in enumerate(rows)
                Cp[cj - off, i] = -Jw[r, j]
            end
        end
        return nothing
    end
    # A level at a time, as `_laplace_unit_curvature`.
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
    _laplace_block_precision(laplace, U, theta, Ls, u, b, z; mode=_HYBRID_PROFILE)

Block `b`'s precision `I - d2 ll / du_b du_b`, with its own coordinates at `z`
and the rest of the unit's effects where `u` has them: what the adaptive
quadrature scales a leaf's nodes by (`_quadrature_leaf_rule!`) and the
continuation a soft rule's stiff complement by (`_continuation_block_functions`),
at points other than the mode. `mode` is what each member's filter returns for
the deviations it integrates (`CTSEMHybridBuffers`); the quadrature asks for
their exact integral, `_HYBRID_FULL`, the integrand it integrates.

A unit of one block takes it from `_laplace_unit_curvature` there, which for
such a unit is the same Jacobian of the same gradient at the same width, so
the two rules compile no dual type the Laplace fit's own curvature has not.
Each took its own `ForwardDiff.jacobian` of a closure of its own before, and a
closure's tag is a new dual type through the whole filter and its reverse
pass: about 75 s of compilation (local) in the first correction of a session
on ord4, for numbers that come out the same to the bit. A block of a nested
unit takes the Jacobian over its own coordinates, as both did, through this one
closure rather than two.
"""
function _laplace_block_precision(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64},
    b::Integer, z::AbstractVector; mode::Int=_HYBRID_PROFILE)
    blocks = laplace.units.blocks[U]
    block = blocks[b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    work = copy(u)
    @inbounds for (t, c) in enumerate(columns)
        work[c] = z[t]
    end
    length(blocks) == 1 &&
        return Matrix{Float64}(_laplace_unit_curvature(laplace, U, theta, Ls, work;
            mode=mode).diag[1])
    members = block.members
    gradient_of = function (zz)
        local S, ws, uu, r
        S = eltype(zz)
        ws = _laplace_workspace!(laplace, S, length(theta))
        uu = convert(Vector{S}, work)
        @inbounds for (t, c) in enumerate(columns)
            uu[c] = zz[t]
        end
        r = _laplace_unit_loglik_gradient(laplace, U, convert(Vector{S}, theta),
            [convert(Matrix{S}, L) for L in Ls], uu, ws, members, nothing, mode)
        return [r.gradient[c] for c in columns]
    end
    A = _ctsem_jacobian(gradient_of, collect(Float64, z), _CTSEM_CURVATURE_WIDTHS[])
    return Matrix{Float64}(LinearAlgebra.I, k, k) .- _laplace_symmetrise(A)
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
    blocks::Vector{CTSEMLaplaceBlock};
    store::Union{Nothing,_LaplaceFactorStore{T}}=nothing) where {T}
    ok, logdetM, factors, coupling = _laplace_block_factor(M, blocks; store=store)
    ok && return (repaired=false, ok=true, logdet=logdetM, factors=factors,
        coupling=coupling)
    scale = maximum((maximum(abs, d) for d in M.diag); init=one(real(T)))
    scale = isfinite(scale) && scale > 0 ? scale : one(real(T))
    shift = sqrt(eps(real(float(one(T)))))
    for _ in 1:30
        for d in M.diag
            for i in axes(d, 1); d[i, i] += shift * scale; end
        end
        ok, logdetM, factors, coupling = _laplace_block_factor(M, blocks; store=store)
        ok && return (repaired=true, ok=true, logdet=logdetM, factors=factors,
            coupling=coupling)
        shift *= 10
    end
    return (repaired=true, ok=false, logdet=T(NaN), factors=factors,
        coupling=coupling)
end

"""
    _laplace_solve_unit_mode!(laplace, U, values, Ls; retain)

Newton's method on `g_U`, from the origin -- or from the retained mode for
unit `U` when `retain` is true, which by default it is when
`ctsem_set_warm_start!` has asked for that -- and then from points along that
mode's soft directions and its reflection, where the line to one shows a second
maximum, keeping the mode with the largest Laplace term
(`_laplace_probe_modes`). At every width: a cap of five effects, taken from the
quadrature's (whose cost is exponential in the width, where probing's is not),
left a loading model of six effects per person, at the default start's
loading means of zero, at each person's stationary origin between two maxima
-- its objective 4800 below the true, its gradient off tenfold (dev1).

A solve that `ctsem_set_warm_start!` started from a retained mode is not
probed. Those warm starts are the differencing schemes' (the Hessian's
columns, the continuation's Hessian): their base point was solved from the
origin and probed, and every differenced point belongs in the basin the base
chose, since the derivative there is that mode's. A probe can only cost time
or, switching basin between a column's two points, difference across the
jump. On the Fisher example's final model (11 effects a person) the probes
made each gradient 1.6 times slower.

`g_U` is a log likelihood minus a quadratic, so its curvature is negative
definite near the mode and Newton is the right method; away from the mode, and
for a nonlinear process model, it need not be. Two guards, both visible in the
diagnostics rather than silent: a curvature that is not negative definite is
shifted until it is and the unit flagged as repaired, and a step that does not
improve `g_U` is halved before the iteration gives up.

`retain` per call is for a caller that chains related solves and owns what
`laplace.modes[U]` holds between them: the sequential pass starts each
truncated record from the mode of the record one row shorter, and the first
from the origin, so its result is still a function of theta alone.
"""
_laplace_solve_unit_mode!(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}};
    retain::Bool=_LAPLACE_WARM_START[]) =
    _laplace_solve_unit_mode!(laplace, U, values, Ls,
        _laplace_workspace!(laplace, Float64, length(values)); retain=retain)

# With the slot's workspace in hand: called with it concretely typed, Newton and
# the probes below are static calls, and nothing they return is boxed.
function _laplace_solve_unit_mode!(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws;
    retain::Bool=_LAPLACE_WARM_START[])
    d = laplace.units.dims[U]
    # Newton copies its start, so neither needs a vector of its own.
    origin = _laplace_scratch_vector!(laplace, Float64, d, :mode_origin)
    retained = retain ? laplace.modes[U] : fill!(origin, 0.0)
    started_warm = any(!iszero, retained)
    warm = _laplace_newton_unit_mode(laplace, U, values, Ls, aws, retained)
    best = warm
    fromwarm = started_warm
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
            fill!(origin, 0.0))
        if cold.converged || (isfinite(cold.value) &&
                (!isfinite(best.value) || cold.value > best.value))
            best = cold
            fromwarm = false
        end
    end
    # Not where every unit's inner problem is a concave quadratic: there is no
    # second maximum to find, and screening for one cost 25-50% of an
    # evaluation on random-intercept models, values unchanged
    # (scripts/gsm-uncertainty-2026-10/warm_reference.R).
    if _LAPLACE_MODE_PROBES[] && !laplace.quadratic && best.converged && 0 < d &&
            !(fromwarm && _LAPLACE_WARM_START[])
        best = _laplace_probe_modes(laplace, U, values, Ls, aws, best)
    end
    laplace.modes[U] = best.u
    laplace.inner_iterations[U] = best.iterations
    laplace.inner_gradient[U] = best.gradmax
    laplace.inner_converged[U] = best.converged
    laplace.hessian_repaired[U] = best.repaired
    return (u=best.u, value=best.value, converged=best.converged,
        final=best.final)
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
    start::AbstractVector{Float64})
    d = laplace.units.dims[U]
    blocks = laplace.units.blocks[U]
    # The current point and the trial, with their gradients, swapped when a step
    # is accepted; the curvature, its factorization and the step. All the
    # slot's: a step uses each once, and building them fresh was most of what a
    # small unit's mode solve allocated. `start` is read, not written. Only the
    # mode found is the caller's, copied out at the end, because it is kept
    # (`laplace.modes[U]`, and every mode the probes reach).
    u = copyto!(_laplace_scratch_vector!(laplace, Float64, d, :newton_u1), start)
    utry = _laplace_scratch_vector!(laplace, Float64, d, :newton_u2)
    gcur = _laplace_scratch_vector!(laplace, Float64, d, :newton_g1)
    gtry = _laplace_scratch_vector!(laplace, Float64, d, :newton_g2)
    step = _laplace_scratch_vector!(laplace, Float64, d, :newton_step)
    Mwork = _laplace_scratch_blockmatrix!(laplace, Float64, U, blocks, :newton_M)
    fstore = _laplace_factor_store!(laplace, Float64, U)
    repaired = false
    converged = false
    iterations = 0
    # The last exact curvature formed, before any repair, for
    # `_laplace_probe_modes`: one step from the mode, so its directions are the
    # mode's, never an SR1 update. The slot's scratch, so read it before the
    # next solve on this slot.
    curvature = nothing
    keep = _LAPLACE_MODE_PROBES[] && d > 0
    # The curvature is most of what a mode solve costs (60-70% of a Laplace
    # evaluation, dev1), and Newton formed it afresh at every iteration. Here it
    # is formed exactly at the first iterate -- so the solve stays a function of
    # the parameters alone, from the origin as before; started instead from the
    # prior's identity, which forms nothing there, solves were 5-15% slower on
    # wide2 and bigre (dev2), the first poor step forming it anyway after the
    # steps it wasted -- and then updated by
    # SR1 from each accepted step's change in the gradient, as the finish's walk
    # updates its Hessian, until a step shows the matrix no longer describes the
    # unit: one the line search had to shorten, or one that did not halve the
    # gradient. Then it is formed exactly again. The Laplace term at the mode
    # forms its own exact curvature (`_laplace_unit_primal`), so nothing the
    # objective reports rests on an updated matrix. The probes read the last
    # exact curvature, which every solve ends on (below). Not on a unit of
    # nested blocks, whose curvature is stored block-sparse and a rank-one
    # update would fill.
    #
    # Where the solve ends must not depend on which, and the gradient bound
    # alone did not say it: it is relative to the unit's value
    # (`_laplace_inner_tolerance`), so it says the gradient is small, not that
    # the mode is near, and a solve could stop anywhere up to the bound over
    # the curvature short of it -- 2e-8 for exact Newton on a 15-effect model,
    # more for an updated curvature, which converges only superlinearly and
    # stopped just under the bound. On the affect pilot (35 effects a unit)
    # SR1 moved the objective 1.7e-6, above the outer gap tolerance, and the
    # gradient 1e-4 in the loadings, against Newton; neither was the mode. So
    # the solve ends only when the Newton step on the exact curvature at the
    # point is under `inner_tol`, and that curvature is the term's too
    # (`final`). An SR1 solve reaches it by its own steps: when the updated
    # curvature's step is under `inner_tol`, the exact one is formed at the
    # same point to verify it. Neither the step bound nor a failed line search
    # is read from an updated matrix: like a repaired one, it does not measure
    # the distance to the mode.
    sr1 = _LAPLACE_INNER_SR1[] && d > 0 && length(blocks) == 1
    # Taken whether or not they are used: a local that is either `nothing` or
    # an array makes every call on it dynamic.
    B = _laplace_scratch_matrix!(laplace, Float64, d, d, :newton_sr1)
    sstep = _laplace_scratch_vector!(laplace, Float64, d, :newton_sr1_s)
    rvec = _laplace_scratch_vector!(laplace, Float64, d, :newton_sr1_r)
    exact = true
    # Whether the step that reached `u` was taken on an exact curvature: the
    # origin counts as one, as it always did. And whether `fstore` holds an
    # exact, unrepaired factorization, whose solve measures the distance left.
    stepped_exact = true
    factored_exact = false
    # The exact curvature at the point the solve ends on, formed there for the
    # last decision, which the term then reads instead of forming it again
    # (`_laplace_unit_primal`'s `given`): one curvature a unit fewer.
    final = nothing
    current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
        gcur)
    if !isfinite(current.value) && any(!iszero, u)
        fill!(u, 0.0)
        current = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
            gcur)
    end
    for iteration in 1:laplace.inner_maxiter
        iterations = iteration
        small = d == 0 || maximum(abs, current.gradient) <
            _laplace_inner_tolerance(laplace, current.value)
        # Nothing exact to measure the distance with -- the origin, or a
        # curvature the factorization repaired: the bound is all there is.
        if small && stepped_exact && !factored_exact
            converged = true
            break
        end
        # Under the bound the distance to the mode is measured on the exact
        # curvature at this point, which is then the term's curvature too, so
        # the measurement costs nothing more. An SR1 solve first walks on, a
        # gradient a step, until its own step is under `inner_tol` (below).
        small && !sr1 && (exact = true)
        # The band, not just the slot. Each Newton step builds the curvature
        # over every block, and that is the bulk of the primal -- passing only
        # the slot alone left the count at its default, so the whole mode solve
        # ran serially while everything else divided. It cost 4.6x down to
        # 1.2x on the value-only path and no test could have reported it.
        # Asserted: inference does not see through the curvature's routes, and
        # an untyped `M` made the factorization below a dynamic call that boxed
        # its result, every Newton step.
        M = if exact
            Mx = _laplace_unit_curvature(laplace, U, values, Ls, u;
                into=Mwork)::CTSEMBlockMatrix{Float64}
            sr1 && copyto!(B, Mx.diag[1])
            Mx
        else
            copyto!(Mwork.diag[1], B)
            Mwork
        end
        used_exact = exact
        keep && used_exact && (curvature = _laplace_block_dense!(_laplace_scratch_matrix!(laplace,
            Float64, d, d, :probe_curvature), M, laplace.units.blocks[U]))
        fac = _laplace_factor_repaired!(M, blocks; store=fstore)
        factored_exact = used_exact && fac.ok && !fac.repaired
        repaired |= fac.repaired
        fac.ok || break
        _laplace_block_solve!(step, fac.factors, fac.coupling, blocks,
            current.gradient)
        # The step, too, says the mode is found: `M^-1 g` is how far it is, in
        # the effects' own standardised units, whatever the curvature. The
        # gradient bound alone cannot say so where the curvature is large,
        # because roundoff in `g` grows with it. A random T0MEANS whose state's
        # remaining variance is 1e-6 under no measurement error has a
        # curvature of 1e7, and its unit sat at |g| 3.3e-10 against a
        # tolerance of 1.4e-10 for 200 iterations with a step of 3e-17, then
        # failed the trial point. Not from a repaired curvature, whose step
        # does not measure the distance to a maximum.
        if !fac.repaired && maximum(abs, step) < laplace.inner_tol
            if exact
                converged = true
                final = _laplace_copy_blockmatrix(M)
                break
            end
            exact = true
            continue
        end
        accepted = false
        # Whether any trial could be evaluated at all, which is what separates
        # "there is nothing left to gain here" from "the step went somewhere
        # the model has no value". See the `stationary` branch below.
        evaluable = false
        scale = 1.0
        for _ in 1:20
            utry .= u .+ scale .* step
            trial = _laplace_unit_objective_gradient(laplace, U, values, Ls,
                utry, aws, gtry)
            evaluable |= isfinite(trial.value)
            if isfinite(trial.value) && trial.value >= current.value - 1e-12
                u, utry = utry, u
                gcur, gtry = gtry, gcur
                current = trial
                accepted = true
                break
            end
            scale /= 2
        end
        accepted && (stepped_exact = used_exact)
        if accepted && sr1
            # `u` is the accepted point and `utry` the one it left; `gcur` and
            # `gtry` their gradients. The curvature is minus the gradient's
            # Jacobian, so its secant is the old gradient less the new.
            @inbounds for i in 1:d
                sstep[i] = u[i] - utry[i]
            end
            @inbounds for i in 1:d
                acc = 0.0
                for j in 1:d
                    acc += B[i, j] * sstep[j]
                end
                rvec[i] = (gtry[i] - gcur[i]) - acc
            end
            den = dot(rvec, sstep)
            if abs(den) > 1e-8 * norm(rvec) * norm(sstep)
                @inbounds for j in 1:d, i in 1:d
                    B[i, j] += rvec[i] * rvec[j] / den
                end
            end
            exact = scale < 1 ||
                maximum(abs, gcur) > maximum(abs, gtry) / 2
        end
        if !accepted && !exact
            # An updated curvature's step found nothing: form it exactly and
            # try again before any verdict.
            exact = true
            continue
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
            #
            # But a value that cannot rank the trials does not mean the mode is
            # found: on an exact, unrepaired curvature the Newton step is
            # reliable whatever the value can see, and the gradient still ranks
            # it. On the affect pilot, well-conditioned units (curvature 1.3 to
            # 134) stopped here one step short -- |g| 25-100 times the bound,
            # the step 1e3-1e4 times `inner_tol`, 1e-13 left to gain against a
            # value of -1000 -- and which units did depended on the CPU's last
            # bits: the objective differed by 1.1e-6 between Intel and AMD
            # nodes, through the log-determinant at those modes. So that step
            # is taken if it makes the gradient smaller, and only otherwise is
            # this the floor.
            stationary = evaluable &&
                _laplace_stationary_gain(current.gradient, step) <=
                    _laplace_inner_tolerance(laplace, current.value)
            if stationary && exact && !fac.repaired
                utry .= u .+ step
                trial = _laplace_unit_objective_gradient(laplace, U, values, Ls,
                    utry, aws, gtry)
                if isfinite(trial.value) &&
                        maximum(abs, gtry) < maximum(abs, current.gradient)
                    u, utry = utry, u
                    gcur, gtry = gtry, gcur
                    current = trial
                    stepped_exact = true
                    continue
                end
            end
            converged = stationary
            stationary && exact && !fac.repaired &&
                (final = _laplace_copy_blockmatrix(M))
            break
        end
    end
    # A solve that ran out of iterations has not measured its distance to the
    # mode, whatever its gradient: the bound is not the stopping rule above, so
    # counting it here let a starved budget stop short and still call the mode
    # found. The Hessian's escalation, which raises the budget only for a unit
    # that did not converge, then never ran, and its columns differed from the
    # full budget's by 3e-8 relative (test-julia-laplace.R, 2026-10-18).
    d == 0 && (converged = true)
    # The gradient's largest entry, the one thing its readers ask of it: the
    # vector is the slot's, and the next solve on this slot overwrites it.
    return (u=copy(u), value=current.value,
        gradmax=d == 0 ? 0.0 : maximum(abs, current.gradient),
        converged=converged, iterations=iterations, repaired=repaired,
        curvature=curvature, final=final)
end

"""
    _laplace_probe_modes(laplace, U, values, Ls, aws, found)

The unit's mode with the largest Laplace term among `found`, the one Newton
reached from the origin, and those Newton reaches from `_LAPLACE_PROBE_REACH`
posterior sds either way (at most three prior sds) along each of `found`'s soft
directions -- precision below `_LAPLACE_PROBE_SOFT`, a negative one included --
and from `-found.u`, its reflection through the origin.

Why. Where `g_U` has two maxima, Newton from the origin reaches whichever one
the origin's basin holds, and as the parameters move the basin's edge crosses
the origin: the unit's term jumps to the other mode's, and the objective with
it. On the affect pilot's person-persistence fit the optimiser climbed into
exactly such an edge and stalled there: one person with modes 5 log units apart
in `g_U`, a valley 0.15 below the lower, the term 7.7 lower on the far side of
the edge and central differences of the objective growing as one over the
step. The lower mode carried `exp(-7.7)` of the higher's mass, so the higher
mode's term was the Laplace approximation, and the jump only which basin the
origin fell in. A lower mode is often shallow along the way to a higher one, so
its soft directions are where to look.

Not always, and the reflection is for when not. On the pilot's
person-persistence fit (rank-3 loadings, three effects a person) three
people's better modes lay 4 to 6 posterior sds from the found ones, past a
valley, with the found mode pinned
in every direction in two of them (least precision 3.7 and 4.3): a person's
data explained on one loading dimension or on another, at -3.2 on the first or
+3.0 on the second. No soft-direction probe ran or reached them, and the
objective was 12.9 short there. Newton from `-found.u` found all three; probes
at 6 posterior sds along every direction found two, at ten times the cost; at
the reflection the found mode's quadratic is wrong by construction, so the
screen passes it, and near the origin, where the reflection is the found mode,
the screen stops it.

A probe is screened before Newton runs from it, at the cost of one gradient:
where `g_U` is the quadratic the found mode's curvature `M` describes, the
gradient at `p` is `-M (p - u)` and a Newton step from `p` lands on `u`. Newton
runs only where the actual gradient `a` would land it more than one posterior
sd away, `(a + M (p - u))' M^-1 (a + M (p - u)) > 1`, or where `M` does not
factor. Most probes return to the mode they started from, and this is what
they cost then.

Between two maxima the term (`_laplace_unit_primal`), not `g_U`, decides:
where their terms cross, the larger is the better single-mode approximation of
the integral and the objective stays continuous there, which choosing by height
does not. But only between maxima (`_laplace_modes_separated`): where `g_U`
climbs from one point to the other without a dip, the lower is a point on the
higher's slope, where Newton stopped on a nearly flat stretch, and the higher is
the mode. Comparing terms there kept the stopping point whenever its curvature
was soft enough for the gated rule to score it above plain Laplace at the mode
-- on the pilot, a person 0.18 below their mode in `g_U` and 0.51 above it in
term -- so the unit's term depended on where Newton stopped.

It is still a search: a better mode no probe reaches is not found, and where a
probe starts to reach one the objective changes by that mode's gain. A unit
with two modes of comparable mass is one no single-mode term approximates, and
the quadrature correction is the remedy there (`ctLaplaceCheck`).
"""
function _laplace_probe_modes(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws, found)
    d = length(found.u)
    M = found.curvature === nothing ?
        _laplace_block_dense(_laplace_unit_curvature(laplace, U, values, Ls,
            copy(found.u)), laplace.units.blocks[U], d) : found.curvature
    # A copy: `found.curvature` is the slot's scratch, and the probes'
    # own solves overwrite it. This copy, the decomposition and the factor are
    # the slot's too, under tags no solve uses: they run for every unit,
    # probed or not.
    Msym = _laplace_scratch_matrix!(laplace, Float64, d, d, :probe_M)
    @inbounds for j in 1:d, i in 1:d
        Msym[i, j] = (M[i, j] + M[j, i]) / 2
    end
    M = Msym
    E = _ctsem_symeig!(_laplace_scratch_symeig!(laplace, d), M)
    all(isfinite, E.values) || return found
    # The values ascend, so the soft directions are the first `nsoft`.
    nsoft = 0
    while nsoft < d && E.values[nsoft + 1] < _LAPLACE_PROBE_SOFT
        nsoft += 1
    end
    reflected = any(!iszero, found.u)
    (nsoft > 0 || reflected) || return found
    F = _ctsem_cholesky(copyto!(_laplace_scratch_matrix!(laplace, Float64, d, d,
        :probe_F), M), d)
    screened = issuccess(F)
    best = found
    bestterm = NaN
    seen = [found.u]
    # The screen's working vectors, the slot's: none outlives one probe.
    atgrad = _laplace_scratch_vector!(laplace, Float64, d, :probe_grad)
    miss = _laplace_scratch_vector!(laplace, Float64, d, :probe_miss)
    solved = _laplace_scratch_vector!(laplace, Float64, d, :probe_solve)
    # Each soft direction either way, then the reflection.
    for probe in 1:(2 * nsoft + (reflected ? 1 : 0))
        # The slot's: Newton reads its start and does not keep it.
        start = _laplace_scratch_vector!(laplace, Float64, d, :probe_start)
        if probe <= 2 * nsoft
            j = (probe + 1) >> 1
            reach = min(_LAPLACE_PROBE_REACH / sqrt(max(E.values[j], 1e-8)), 3.0)
            sign = isodd(probe) ? -1.0 : 1.0
            start .= found.u .+ (sign * reach) .* view(E.vectors, :, j)
        else
            start .= .-found.u
        end
        if screened
            at = _laplace_unit_objective_gradient(laplace, U, values, Ls, start,
                aws, atgrad)
            isfinite(at.value) || continue
            # miss = a + M (p - u)
            @inbounds for i in 1:d
                acc = at.gradient[i]
                for k in 1:d
                    acc += M[i, k] * (start[k] - found.u[k])
                end
                miss[i] = acc
            end
            dot(miss, ldiv!(solved, F, miss)) > 1 || continue
            # And only where a second maximum can lie on the line: past the
            # found one the line falls, and a maximum beyond a valley makes it
            # rise again, so the slope along it turns outward somewhere. The
            # screen above passes wherever the unit is not Gaussian along the
            # line, a maximum or not, and a Newton solve from every such start
            # made probing units of 15 effects cost four times their solve and
            # find nothing (dev1). The slope at the start, then at three
            # quarters, a half and a quarter of the way; Newton starts at the
            # first point that rises outward, which lies in the other basin.
            # A maximum narrower than a quarter of the line is missed.
            outward = dot(at.gradient, start) - dot(at.gradient, found.u) > 0
            if !outward
                ray = _laplace_scratch_vector!(laplace, Float64, d, :probe_ray)
                raygrad = _laplace_scratch_vector!(laplace, Float64, d,
                    :probe_raygrad)
                for t in (0.75, 0.5, 0.25)
                    @inbounds for i in 1:d
                        ray[i] = found.u[i] + t * (start[i] - found.u[i])
                    end
                    rg = _laplace_unit_objective_gradient(laplace, U, values, Ls,
                        ray, aws, raygrad)
                    isfinite(rg.value) || continue
                    if dot(rg.gradient, start) - dot(rg.gradient, found.u) > 0
                        copyto!(start, ray)
                        outward = true
                        break
                    end
                end
            end
            outward || continue
        end
        r = _laplace_newton_unit_mode(laplace, U, values, Ls, aws, start)
        (r.converged && isfinite(r.value)) || continue
        any(m -> _laplace_maxabsdiff(m, r.u) < 1e-4, seen) && continue
        push!(seen, r.u)
        if _laplace_modes_separated(laplace, U, values, Ls, aws, best, r)
            isnan(bestterm) && (bestterm = _laplace_unit_primal(laplace, U,
                values, Ls, copy(best.u), aws).term)
            term = _laplace_unit_primal(laplace, U, values, Ls, copy(r.u), aws).term
            (isfinite(term) && (!isfinite(bestterm) || term > bestterm)) || continue
        else
            r.value > best.value || continue
            term = NaN
        end
        best = (u=r.u, value=r.value, gradmax=r.gradmax, converged=true,
            iterations=found.iterations + r.iterations,
            repaired=found.repaired || r.repaired, curvature=nothing,
            final=r.final)
        bestterm = term
    end
    return best
end

"""`maximum(abs, a .- b)` without the difference vector."""
function _laplace_maxabsdiff(a::AbstractVector, b::AbstractVector)
    m = 0.0
    @inbounds for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end
    return m
end

"""
    _laplace_modes_separated(laplace, U, values, Ls, aws, a, b)

Whether two of unit `U`'s Newton end points are separate maxima of `g_U`: it
dips below the lower of the two somewhere on the segment between them (seven
interior points), or cannot be evaluated there. Otherwise the lower is a point
on the higher's slope. See `_laplace_probe_modes`.

Any dip, however shallow: a shoulder -- a maximum a valley of hundredths of a
nat parts from a higher one -- has a soft curvature, and its larger term is
the better one. On a bench model (gD3) two units of two effects had shoulders
0.009 and 0.02 nats deep, and against each unit's integral on a grid the
shoulder's term was 0.007 and 0.46 off, the higher maximum's 0.54 and 0.61.
A depth threshold would choose the higher maximum there. What does go wrong
with shoulders is the search: from the higher maximum no probe reached
either, so which term a unit got depended on which maximum Newton reached
first, and the SR1 and exact solves reached different ones.
"""
function _laplace_modes_separated(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{Float64}, Ls::Vector{Matrix{Float64}}, aws, a, b)
    lower = min(a.value, b.value)
    margin = 1e-8 * (1 + abs(lower))
    point = _laplace_scratch_vector!(laplace, Float64, length(a.u), :sep_point)
    grad = _laplace_scratch_vector!(laplace, Float64, length(a.u), :sep_grad)
    for k in 1:7
        t = k / 8
        point .= a.u .+ t .* (b.u .- a.u)
        g = _laplace_unit_objective_gradient(laplace, U, values, Ls, point, aws,
            grad).value
        (isfinite(g) && g >= lower - margin) || return true
    end
    return false
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
direction, and by `_laplace_min_eigenvalue` for the reports. A wider unit keeps
the total floor under `:gated` -- a dense eigendecomposition of a block-sparse
unit curvature is `O(d^3)` and fills in what the block elimination keeps sparse
-- and its smallest eigenvalue is found by bisection on the inertia instead.
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
    blocks::Vector{CTSEMLaplaceBlock}, shift::Real=1.0;
    shifted::Union{Nothing,CTSEMBlockMatrix{T}}=nothing,
    store::Union{Nothing,_LaplaceFactorStore{T}}=nothing) where {T}
    # `shifted` and `store`: a caller's storage for the shifted copy and its
    # factorization, neither of which outlives this call -- the gated floor's
    # gate passes the slot's, once per unit per evaluation.
    if shifted === nothing
        shifted = CTSEMBlockMatrix{T}([copy(d) for d in M.diag],
            [[copy(c) for c in row] for row in M.coupling])
    else
        for b in eachindex(M.diag)
            copyto!(shifted.diag[b], M.diag[b])
            for t in eachindex(M.coupling[b])
                copyto!(shifted.coupling[b][t], M.coupling[b][t])
            end
        end
    end
    for d in shifted.diag
        for i in axes(d, 1); d[i, i] -= shift; end
    end
    ok, _, _, _ = _laplace_block_factor(shifted, blocks; store=store)
    return ok
end

"""
    _laplace_count_below(M, blocks, shift)

How many eigenvalues of `M` lie below `shift`: the inertia of `M - shift I`,
by Sylvester's law of inertia, read off the same innermost-first elimination
`_laplace_block_factor` performs. The elimination is a congruence with a block
unit-triangular matrix, so `M - shift I` and the block diagonal of eliminated
pivots have the same inertia; each pivot is decomposed by `_ctsem_symeig`
rather than Cholesky-factored, so an indefinite one is counted instead of
refused. `nothing` when a pivot is singular to working precision -- `shift` is
an eigenvalue there, and the caller treats it as one.
"""
function _laplace_count_below(M::CTSEMBlockMatrix{Float64},
    blocks::Vector{CTSEMLaplaceBlock}, shift::Real)
    nb = length(blocks)
    diag = [Matrix{Float64}(d) for d in M.diag]
    for d in diag
        for i in axes(d, 1); d[i, i] -= shift; end
    end
    coupling = [[Matrix{Float64}(c) for c in row] for row in M.coupling]
    count = 0
    for b in 1:nb
        E = _ctsem_symeig(diag[b])
        scale = max(maximum(abs, E.values; init=0.0), 1.0)
        any(v -> abs(v) <= 1e-13 * scale, E.values) && return nothing
        count += Base.count(<(0.0), E.values)
        ancestors = blocks[b].ancestors
        isempty(ancestors) && continue
        Dinv = E.vectors * Diagonal(1 ./ E.values) * transpose(E.vectors)
        W = [Dinv * coupling[b][t] for t in eachindex(ancestors)]
        for t in eachindex(ancestors)
            a = ancestors[t]
            diag[a] .-= _laplace_symmetrise(transpose(coupling[b][t]) * W[t])
            for s in eachindex(ancestors)
                s == t && continue
                slot = findfirst(==(ancestors[s]), blocks[a].ancestors)
                slot === nothing && continue
                coupling[a][slot] .-= transpose(coupling[b][t]) * W[s]
            end
        end
    end
    return count
end

"""
    _laplace_min_eigenvalue(M, blocks, d; rtol)

The smallest eigenvalue of a unit's curvature, for the reports. A unit up to
`_LAPLACE_EIGEN_MAXDIM` wide is decomposed densely by `_ctsem_symeig`; a wider
one by bisection on `_laplace_count_below`, which needs only block eliminations
-- no fill-in, `O(sum_b k_b^3)` each, every width -- to relative precision
`rtol`. The bisection halves from 1 until no eigenvalue lies below (doubling
downwards first if one is negative), then narrows geometrically.

Wide units used to be reported as `NaN`, which the fit's print counts as below
the prior's but never as near-singular: on the SNSF pilot a study unit of 93
coordinates reached an eigenvalue of 2.5e-4 while the print said nothing
(`CT-SEM/review/LAPLACE-nested-spike-2026-10-02.md`).
"""
function _laplace_min_eigenvalue(M::CTSEMBlockMatrix{Float64},
    blocks::Vector{CTSEMLaplaceBlock}, d::Integer; rtol::Real=1e-4)
    d == 0 && return Inf
    d <= _LAPLACE_EIGEN_MAXDIM[] &&
        return _ctsem_symeig(_laplace_block_dense(M, blocks, d)).values[1]
    below(s) = (c = _laplace_count_below(M, blocks, s); c === nothing || c > 0)
    hi = 1.0
    below(hi) || return Inf
    lo = 0.0
    if below(0.0)
        # An eigenvalue at or below zero: walk down until none is below.
        step = 1.0
        lo = -step
        while below(lo)
            hi = lo
            step *= 2
            lo = -step
            step > 1e12 && return -Inf
        end
    else
        # Halve until nothing is below, so a tiny positive eigenvalue is
        # bracketed in relative terms before the geometric narrowing.
        s = hi / 2
        while below(s)
            hi = s
            s /= 2
            s < 1e-14 && return s
        end
        lo = s
    end
    for _ in 1:200
        (hi - lo) <= rtol * max(abs(hi), abs(lo)) && break
        mid = (lo > 0 && hi > 0) ? sqrt(lo * hi) : (lo + hi) / 2
        below(mid) ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
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
#
# Provenance of the constants, as Appendix B of
# `CT-SEM/review/OPTIM-consolidation-plan-2026-09-25.md` sorts them. The band
# (`lo = 0.2`, `hi = 0.7`) and its C1 hand-off, and the rule itself (3 nodes,
# one Newton step for the complement at each), are measured rather than
# derived, and they move where a fit ends. They were set on the weak-data
# families A, B and C -- 40 subjects, a random drift beside random T0MEANS and
# CINT (A), a random diffusion (B), a random measurement error (C) -- against
# the exact per-unit reference (`LAPLACE-eigenwise-floor-2026-09-23.md`,
# `LAPLACE-gated-gaps-2026-09-24.md` section 1). Categorical indicators
# (family D) and nested units (family N) were only checked afterwards, and
# the constants were not refitted there. The `1e-3` degenerate-pair gap is a
# numerical guard, and `_LAPLACE_EIGEN_MAXDIM` a cost bound.


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

`objective`, `curvature` and `blocks` evaluate the unit somewhere else than its
own coordinates -- `(value, gradient, hlogdet)` and the curvature at a point of
that space, and its block layout -- and `dcheck` is the dimension the
eigenvalue cap compares. With none given, an objective whose filter integrates
deviations (`CTSEMHybrid`) is scored by `_hybrid_gated_term`, which decides on
the joint curvature and comes back here in joint coordinates; `cross` is the
`CTSEMHybridCross` a caller read with `M`.

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
    aws; M=nothing, logdetM=nothing, inner=nothing, cross=nothing, objective=nothing,
    curvature=nothing, blocks=nothing, dcheck=nothing, flagged::Bool=false) where {T}
    objective === nothing && _hybrid_on(laplace.objective.params.hybrid) &&
        return _hybrid_gated_term(laplace, U, values, Ls, u, aws; M=M,
            logdetM=logdetM, inner=inner, cross=cross, flagged=flagged)
    objf = objective === nothing ?
        (z -> _laplace_unit_objective_gradient(laplace, U, values, Ls, z, aws)) :
        objective
    curvf = curvature === nothing ?
        (z -> _laplace_unit_curvature(laplace, U, values, Ls, z)) : curvature
    blocks === nothing && (blocks = laplace.units.blocks[U])
    inner === nothing && (inner = objf(u))
    g = inner.value
    (isfinite(g) && !isempty(u)) || return (value=g, flagged=false, floor=true)
    M === nothing && (M = curvf(u))
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
        # Joint, as a caller passing `logdetM` passes it: see `_laplace_unit_term`.
        logdetM += 2 * inner.hlogdet
    end
    # `total` only for the fallbacks, which have no eigenvalue to go on;
    # `lap`, unfloored, wherever there is one. See the docstring.
    total = indefinite ? T(NaN) : g - max(logdetM, zero(logdetM)) / 2
    lap = indefinite ? T(NaN) : g - logdetM / 2
    lo, hi = laplace.gate_lo, laplace.gate_hi
    _laplace_exceeds_identity(M, blocks, hi;
        shifted=_laplace_scratch_blockmatrix!(laplace, T, U, blocks, :gate_shifted),
        store=_laplace_factor_store!(laplace, T, U, blocks)) &&
        return (value=lap, flagged=false, floor=false)
    d = length(u)
    (dcheck === nothing ? d : dcheck) > _LAPLACE_EIGEN_MAXDIM[] &&
        return (value=total, flagged=false, floor=true)
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
                h = abs(x) < 1e-12 ? g : objf(ub).value
            else
                Bc = _laplace_clip_first_order(transpose(Vh) * dense * Vh)
                F = _ctsem_cholesky(Matrix(_laplace_symmetrise(Bc)), d - 1)
                h = g - logdet(F) / 2
            end
        else
            ub = u .+ (sqrt(2) * x / sqrt(lam1c)) .* v1
            r = objf(ub)
            Mb = _laplace_block_dense(curvf(ub), blocks, d)
            Bc = _laplace_clip_first_order(transpose(Vh) * Mb * Vh)
            F = _ctsem_cholesky(Matrix(_laplace_symmetrise(Bc)), d - 1)
            issuccess(F) || return (value=T(NaN), flagged=true, floor=false)
            z = F \ (transpose(Vh) * r.gradient)
            h = objf(ub .+ Vh * z).value - logdet(F) / 2
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
    _laplace_unit_reach(laplace, U, theta, Ls)

The raw positions unit `U`'s term can depend on, sorted, or `nothing` if a
member cannot be evaluated: every position some member's log likelihood reads,
from its gradient at the unit's mode -- a position the member's filter never
reads, another study's predictor effect at a predictor value of zero or an
item the member never answers, has a derivative of exactly zero there -- and
every population parameter of the model's levels, which shape the unit's
curvature whatever the gradient at one point says.
"""
function _laplace_unit_reach(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{<:AbstractMatrix})
    spec = laplace.spec
    units = laplace.units
    members = units.members[U]
    npar = length(theta)
    reach = fill(false, npar)
    aws = _laplace_workspace!(laplace, Float64, npar)
    grad = _laplace_scratch_vector!(laplace, Float64, npar, :reach_grad)
    shift = _laplace_scratch_vector!(laplace, Float64, npar, :reach_shift)
    u = laplace.modes[U]
    for m in eachindex(members)
        shifted = _laplace_member_values!(shift, theta, spec, Ls, u, units.offsets[U][m])
        loglik = _laplace_subject_value_gradient!(grad,
            laplace.objective.subject_objectives[members[m]], aws, shifted, nothing,
            _HYBRID_FULL)
        isfinite(loglik) || return nothing
        @inbounds for t in 1:npar
            iszero(grad[t]) || (reach[t] = true)
        end
    end
    for l in eachindex(spec.levels), t in _laplace_level_positions(spec, l)
        reach[t] = true
    end
    return findall(reach)
end

"""
    _laplace_gated_unit_gradient(laplace, U, theta, curvature)

The exact gradient of a flagged unit's gated term: ForwardDiff over the term,
the inner mode's first derivative from `_laplace_dual_unit_mode`, along the
positions the unit can depend on (`_laplace_unit_reach`) alone. The others'
derivatives are exactly zero, and on a model with per-study effects and items
a subject never answers they are most of the vector: each forward pass here
repeats the whole gated rule, so the count of them is the cost.
"""
function _laplace_gated_unit_gradient(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, curvature)
    # The reverse pass (laplace_gated_adjoint.jl), at a cost that does not grow
    # with the parameter count; what follows is the route it falls back to
    # where it does not apply, and the oracle it is tested against.
    if _LAPLACE_GATED_ADJOINT[]
        adjoint = _laplace_gated_unit_adjoint(laplace, U, theta)
        adjoint === nothing || return adjoint
    end
    uhat = laplace.modes[U]
    # A repaired unit's factors are of a shifted curvature, which is not the
    # implicit function theorem's matrix; take the mode's derivative from the
    # unshifted one's eigendecomposition instead (it is nonsingular wherever
    # the mode moves smoothly, whatever its signature).
    # One assignment: `term_of` captures it (a captured local assigned more
    # than once is a `Core.Box`).
    dense_inverse = if laplace.mode_repaired[U]
        Ls0 = _laplace_popchols(theta, laplace.spec)
        E = _ctsem_symeig(_laplace_block_dense(_laplace_unit_curvature(laplace, U,
            theta, Ls0, uhat), laplace.units.blocks[U], length(uhat)))
        E.vectors * Diagonal(1 ./ E.values) * transpose(E.vectors)
    else
        nothing
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
        # `flagged`: the primal already scored this unit by the gated rule,
        # so a filter-integrating objective need not decide the gate again.
        return _laplace_gated_term(laplace, U, x, Lsd, ud, wsd; flagged=true).value
    end
    reach = _laplace_unit_reach(laplace, U, theta, _laplace_popchols(theta, laplace.spec))
    reach === nothing && return _ctsem_gradient(term_of, theta)
    term_along = function (y)
        S = eltype(y)
        x = Vector{S}(undef, length(theta))
        @inbounds for t in eachindex(theta)
            x[t] = theta[t]
        end
        @inbounds for (k, t) in enumerate(reach)
            x[t] = y[k]
        end
        return term_of(x)
    end
    gradient = zeros(Float64, length(theta))
    isempty(reach) && return gradient
    partial = _ctsem_gradient(term_along, theta[reach])
    @inbounds for (k, t) in enumerate(reach)
        gradient[t] = partial[k]
    end
    return gradient
end

"""
    _laplace_unit_term(laplace, U, values, Ls, u, aws)

`g_U(u) - logdet(-d2 g_U/du du) / 2`: unit `U`'s contribution to the
approximated log marginal likelihood.

The `(2*pi)^(d/2)` the Laplace approximation produces cancels the
`(2*pi)^(-d/2)` in the standard-normal density of `u`, so no dimension-dependent
constant appears here whatever the unit's size.
"""
_laplace_unit_term(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws) where {T} = _laplace_unit_primal(laplace, U, values, Ls, u, aws).term

"""
One unit's mode and primal into the evaluation's records: its factorization
and curvature, which the gradient phase reads, its flags, and its term -- or,
when the term is not finite, the slot's failure. `false` then. A function of
its own because `aws` is a slot's and untyped, so the call is a dynamic
dispatch: made here once, everything inside is static, where before the mode
solve's and the primal's results each came back boxed through one.
"""
function _laplace_store_unit_primal!(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls, aws, primal_curvature, primal_matrices,
    primal_gated, unit_loglik, chunk_ok, chunk_bad)
    b0 = _laplace_mark()
    solved = _laplace_solve_unit_mode!(laplace, U, theta, Ls, aws)
    _laplace_charge!(_LAPLACE_BYTES_MODE, b0)
    primal = _laplace_unit_primal(laplace, U, theta, Ls, laplace.modes[U], aws;
        given=solved.final)
    laplace.mode_repaired[U] = primal.fac.repaired
    primal_curvature[U] = (primal.fac.factors, primal.fac.coupling)
    primal_matrices[U] = primal.M
    laplace.logdet_floored[U] = primal.floored
    primal_gated[U] = primal.gated
    term = primal.term
    if !isfinite(term)
        chunk_ok[_laplace_slot()] = false
        chunk_bad[_laplace_slot()] = term
        return false
    end
    unit_loglik[U] = term
    return true
end

"""
    _laplace_unit_primal(laplace, U, values, Ls, u, aws)

Unit `U`'s term at `u` as the objective scores it, and what the gradient reads
from the way there: `(term, M, fac, gated, floored)`, the curvature `M` (shifted
in place where `fac` repaired it), its factorization, whether the gated rule
scored the unit, and whether the term is the prior floor's bound. The one
definition: the evaluation's primal, `_laplace_unit_term`, and the inner
solve's choice between modes (`_laplace_probe_modes`) all take it from here, so
a unit cannot be scored one way when its mode is chosen and another when it is
summed. They were three before, and the third returned `NaN` for a curvature
the primal repairs and scores by the gated rule's soft directions.

Under the gated floor, deviations the filter integrates (`CTSEMHybrid`) need the
joint curvature's gate, whose cross terms the same curvature pass gives
(`CTSEMHybridCross`).
"""
function _laplace_unit_primal(laplace::CTSEMLaplaceObjective, U::Integer,
    values::AbstractVector{T}, Ls::Vector{<:AbstractMatrix}, u::AbstractVector{T},
    aws; given=nothing) where {T}
    blocks = laplace.units.blocks[U]
    gated_floor = _LAPLACE_PRIOR_FLOOR[] && laplace.floor === :gated
    b0 = _laplace_mark()
    cross = (gated_floor && !isempty(u) &&
             _hybrid_on(laplace.objective.params.hybrid)) ?
        _hybrid_cross!(laplace, T, length(u), length(laplace.units.members[U])) :
        nothing
    # Asserted, as in `_laplace_newton_unit_mode`.
    # `given`: the exact curvature the inner solve formed at `u` to end on
    # (`_laplace_newton_unit_mode`), not when the gated hybrid needs its cross
    # terms from the same pass.
    M = isempty(u) ? CTSEMBlockMatrix(T, blocks) :
        (T === Float64 && given !== nothing && cross === nothing) ?
            given::CTSEMBlockMatrix{T} :
        _laplace_unit_curvature(laplace, U, values, Ls, u;
            hook=cross)::CTSEMBlockMatrix{T}
    _laplace_charge!(_LAPLACE_BYTES_CURV, b0)
    fac = _laplace_factor_repaired!(M, blocks)
    ok, logdetM = fac.ok, fac.logdet
    b0 = _laplace_mark()
    # The gradient into the slot's: only this call and the gated term below
    # read it, and neither keeps it.
    inner = _laplace_unit_objective_gradient(laplace, U, values, Ls, u, aws,
        _laplace_scratch_vector!(laplace, T, length(u), :primal_grad))
    _laplace_charge!(_LAPLACE_BYTES_OBJ, b0)
    # The joint curvature's logdet: the profile's plus the filter-integrated
    # deviations' (`CTSEMHybrid`), before the floor or the gate reads it.
    ok && (logdetM += 2 * inner.hlogdet)
    floored = ok && _LAPLACE_PRIOR_FLOOR[] && logdetM < 0
    gated = false
    term = if !isfinite(inner.value)
        inner.value
    elseif isempty(u)
        inner.value - inner.hlogdet
    elseif ok
        inner.value - _laplace_prior_floor_logdet(logdetM) / 2
    else
        T(NaN)
    end
    if gated_floor && ok && !isempty(u) && isfinite(inner.value)
        # The factorization shifted a repaired `M` in place, so the gated term
        # gets the curvature afresh and factors it itself.
        g = fac.repaired ?
            _laplace_gated_term(laplace, U, values, Ls, u, aws;
                M=_laplace_unit_curvature(laplace, U, values, Ls, u), inner=inner) :
            _laplace_gated_term(laplace, U, values, Ls, u, aws; M=M,
                logdetM=logdetM, inner=inner, cross=cross)
        if g.flagged
            term = g.value
            gated = true
            # Scored by the gated rule, not floored: the diagnostics count
            # `logdet_floored` as "the term is a bound".
            floored = false
        elseif !g.floor && !fac.repaired
            # Cleared by the gate: plain Laplace, no floor, and so the ordinary
            # seeded gradient rather than the floored one.
            term = inner.value - logdetM / 2
            floored = false
        end
    end
    return (term=term, M=M, fac=fac, gated=gated, floored=floored)
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

# Their order among other tags, which ForwardDiff asks whenever duals of two
# tags meet (`a ≺ b`: `Dual{b}` goes outside `Dual{a}`). A bare tag type has
# none, so any meeting threw `DualMismatchError`: a categorical row's moment
# Jacobian (`CTSEMNestedTag`) reading a seeded value it did not seed itself,
# in the sweep that gives a Laplace fit's per-unit scores, and the seeds
# inside the forward Jacobian that is that sweep's fallback. The seeds go
# outside an ordinary tag (they perturb a point that tag's differentiation is
# at), inside the nested one (which differentiates within whatever reaches it),
# and the outer seed outside the inner, as `_laplace_seed_pair` builds them.
const _LaplaceSeedTag = Union{Type{_LaplaceSeedInner},Type{_LaplaceSeedOuter}}
ForwardDiff.:≺(::Type{_LaplaceSeedInner}, ::Type{_LaplaceSeedOuter}) = true
ForwardDiff.:≺(::Type{_LaplaceSeedOuter}, ::Type{_LaplaceSeedInner}) = false
ForwardDiff.:≺(::Type{<:ForwardDiff.Tag}, ::_LaplaceSeedTag) = true
ForwardDiff.:≺(::_LaplaceSeedTag, ::Type{<:ForwardDiff.Tag}) = false
ForwardDiff.:≺(::Type{ForwardDiff.Tag{CTSEMNestedTag,V}}, ::_LaplaceSeedTag) where {V} =
    false
ForwardDiff.:≺(::_LaplaceSeedTag, ::Type{ForwardDiff.Tag{CTSEMNestedTag,V}}) where {V} =
    true

"""The two seeds `(eps1, eps2)`, one per tag: `x + d1 * eps1 + d2 * eps2`
carries the mixed second derivative along `d1` and `d2` as the coefficient of
`eps1 * eps2`."""
@inline _laplace_seed_pair() = (
    ForwardDiff.Dual{_LaplaceSeedOuter}(ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 1.0),
        ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 0.0)),
    ForwardDiff.Dual{_LaplaceSeedOuter}(ForwardDiff.Dual{_LaplaceSeedInner}(0.0, 0.0),
        ForwardDiff.Dual{_LaplaceSeedInner}(1.0, 0.0)))

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
        # `S1`, not `S`: the order-2 path below assigns `S` too, and a local
        # assigned twice and captured by a closure is boxed, so both closures
        # saw an untyped `S`, untyped scratch vectors, and a dynamic dispatch
        # and a boxed dual on every element of every member's sweep -- 3,700
        # heap objects a gradient evaluation on gC8.
        S1 = typeof(seed)
        one_member! = function (c::Int, m)
            local aws, gradient, x, shifted, t
            aws = _laplace_workspace!(laplace, S1, npar)
            gradient = _laplace_scratch_vector!(laplace, S1, npar, :sweep_grad)
            x = _laplace_scratch_vector!(laplace, S1, npar, :sweep_x)
            shifted = _laplace_member_values!(
                _laplace_scratch_vector!(laplace, Float64, npar, :sweep_shift),
                values, spec, Ls, u, units.offsets[U][m])
            @inbounds for t in 1:npar
                x[t] = shifted[t] + seed * d1[t]
            end
            _laplace_member_pass_finite!(gradient,
                laplace.objective.subject_objectives[unitmembers[m]], aws, x) ||
                return false
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
    e1, e2 = _laplace_seed_pair()
    S = typeof(e1)
    d12 = _laplace_scratch_matrix!(laplace, Float64, npar, nm, :sweep_d12)
    # One member. Writes only column `c` of the three outputs, so members never
    # collide however they are divided; the workspace comes from the ambient
    # slot, and the dual buffers are cached against it.
    one_member! = function (c::Int, m)
        local aws, gradient, x, shifted, t, inner
        aws = _laplace_workspace!(laplace, S, npar)
        gradient = _laplace_scratch_vector!(laplace, S, npar, :sweep_grad)
        x = _laplace_scratch_vector!(laplace, S, npar, :sweep_x)
        shifted = _laplace_member_values!(
            _laplace_scratch_vector!(laplace, Float64, npar, :sweep_shift),
            values, spec, Ls, u, units.offsets[U][m])
        @inbounds for t in 1:npar
            x[t] = shifted[t] + e1 * d1[t] + e2 * d2[t]
        end
        _laplace_member_pass_finite!(gradient,
            laplace.objective.subject_objectives[unitmembers[m]], aws, x) ||
            return false
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

# Whether one member's pass came out finite, deep. The sweeps' workspace is a
# slot's, which `_laplace_workspace!` cannot type, so the pass is a dynamic call:
# returning its dual log likelihood boxed it, once per member per sweep, where a
# `Bool` comes back as itself.
_laplace_member_pass_finite!(gradient, subject_objective, aws, x) =
    _laplace_finite(_laplace_subject_value_gradient!(gradient, subject_objective,
        aws, x))

# Deep: `isfinite` on a dual tests only its value, and a sweep's whole point
# is the partials. Testing the value alone is what let a NaN derivative travel
# from a unit whose predicted variance had collapsed all the way into the
# assembled gradient, where it was indistinguishable from a bad trial point.
@inline _laplace_finite(x::Real) = isfinite(x)
@inline _laplace_finite(x::ForwardDiff.Dual) = _finite_deep(x)

"""
    _laplace_level_chol_solves(spec, Ls, dL)

`Ls[l] \\ dL[l][t]` for every level `l` and parameter `t`, as
`_laplace_seeded_unit_gradient!` takes it; `nothing` under reduced rank, whose
`L` is not square and whose units take the other route (`needV`) there. No unit
changes it, so it is solved once per evaluation: solved by each unit, the
general `\\` was an LAPACK factorization -- and its process-wide lock -- per
unit and parameter, nearly every lock sample left in a Laplace gradient or
Hessian at eight threads (dev1).
"""
_laplace_level_chol_solves(spec, Ls, dL) = hasreducedrank(spec) ? nothing :
    [[Ls[l] \ dL[l][t] for t in eachindex(dL[l])] for l in eachindex(dL)]

# `levelsweeps` for a unit that needs none: never written, so one is enough.
const _LAPLACE_NO_LEVELSWEEPS = Vector{Matrix{Float64}}[]

"""
    _laplace_seeded_unit_gradient!(out, laplace, U, values, Ls, dL, dLsolved, M,
                                   factors, elim)

Accumulate unit `U`'s exact contribution to `dT/dtheta` into `out`, in a number
of sweeps proportional to the unit's members rather than to the parameter count.
`dLsolved` is `_laplace_level_chol_solves(spec, Ls, dL)`.

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
    dLsolved, M::CTSEMBlockMatrix{Float64}, factors, elim)

    spec = laplace.spec
    units = laplace.units
    blocks = units.blocks[U]
    npar = length(values)
    d = units.dims[U]
    nmem = length(units.members[U])
    uhat = laplace.modes[U]
    # Every working array below is the calling slot's, cached by shape
    # (`_laplace_scratch_vector!`), or a worker's own where the member loops
    # divide: built afresh, they were most of what a Laplace gradient
    # allocated, once per unit, block and population parameter.
    zerodir = fill!(_laplace_scratch_vector!(laplace, Float64, npar, :sg_zero), 0.0)
    # Inactive coordinates (`CTSEMLaplaceUnits`): nothing depends on them, so
    # their curvature rows are the identity with no coupling, their columns of
    # the factor of C below are basis vectors, and a sweep along one moves only
    # effects no member here can see -- exactly zero, and skipped.
    act = _laplace_active_columns(laplace, U)
    isactive = fill!(_laplace_scratch_vector!(laplace, Bool, d, :sg_active), true)
    if act !== nothing
        fill!(isactive, false)
        @inbounds for c in act
            isactive[c] = true
        end
    end

    # No random effects anywhere in this unit: there is no integral, the term
    # is the plain log likelihood, and one ordinary reverse sweep is the whole
    # gradient. Returning early without it would silently contribute nothing.
    if d == 0
        # With deviations the filter integrates (`CTSEMHybrid`), the term is
        # their integral, so the gradient is the integral's.
        if _hybrid_on(laplace.objective.params.hybrid)
            full = _laplace_hybrid_full_gradients(laplace, U, values, Ls, uhat)
            full === nothing && return false
            @inbounds for m in 1:nmem, t in 1:npar
                out[t] += full[t, m]
            end
            return true
        end
        pass = _laplace_unit_seeded_gradient(laplace, U, values, Ls, uhat,
            1:nmem, zerodir, zerodir, 1)
        pass.ok || return false
        @inbounds for m in 1:nmem, t in 1:npar
            out[t] += pass.d0[t, m]
        end
        return true
    end

    bsi = _laplace_mark()
    Cdiag, Ccoup = _laplace_selected_inverse_scratch!(laplace, U, factors, elim,
        blocks)
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
            body(mine)
        end
    end

    # Scatter a level-space vector onto the raw parameter vector, into `dir`.
    # The local names must not collide with anything in the enclosing scope: a
    # Julia closure rebinds an enclosing local rather than shadowing it, so
    # assigning `out` here would silently replace the caller's accumulator.
    function scatter!(dir::Vector{Float64}, level::Int, w::AbstractVector)
        local rho, p
        fill!(dir, 0.0)
        rho = spec.levels[level].re_index
        @inbounds for p in eachindex(rho); dir[rho[p]] = w[p]; end
        return dir
    end
    # The same for `L w`, a level's loading applied to block coordinates.
    function scatter_loaded!(dir::Vector{Float64}, level::Int, L::AbstractMatrix,
        w::AbstractVector)
        local rho, p, q, acc
        fill!(dir, 0.0)
        rho = spec.levels[level].re_index
        @inbounds for p in eachindex(rho)
            acc = 0.0
            for q in eachindex(w)
                acc += L[p, q] * w[q]
            end
            dir[rho[p]] = acc
        end
        return dir
    end

    # d psi / d v_m, d loglik_m / d v_m, and the member share of s' B.
    Pm = fill!(_laplace_scratch_matrix!(laplace, Float64, npar, nmem, :sg_Pm), 0.0)
    llvm = fill!(_laplace_scratch_matrix!(laplace, Float64, npar, nmem, :sg_llvm), 0.0)
    Bsm = fill!(_laplace_scratch_matrix!(laplace, Float64, npar, nmem, :sg_Bsm), 0.0)
    # `Vector{Bool}`, not `falses`. A `BitVector` packs 64 flags into a word, so
    # `seen[m] = true` is a read-modify-write of that word: two tasks setting
    # *different* members race on it and one write is simply lost.
    seen = fill!(_laplace_scratch_vector!(laplace, Bool, nmem, :sg_seen), false)

    assemble! = function (mine)
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
        local Cs, dir1, dir2, T1, p, x, y, z, acc
        # This worker's own: `_laplace_scratch_*` is keyed by the calling slot,
        # and every worker of the partition holds a different one.
        ownedA = fill!(_laplace_scratch_vector!(laplace, Bool, nmem, :sg_ownedA), false)
        subsA = _laplace_scratch_vector!(laplace, Int, 0, :sg_subsA)
        dir1 = _laplace_scratch_vector!(laplace, Float64, npar, :sg_dir1)
        dir2 = _laplace_scratch_vector!(laplace, Float64, npar, :sg_dir2)
        @inbounds for m in mine; ownedA[m] = true; end
        for (b, block) in enumerate(blocks)
            l = block.level
            L = Ls[l]
            k = block.size
            k == 0 && continue
            empty!(subsA)
            for m in block.members
                ownedA[m] && push!(subsA, m)
            end
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
        #
        # The engine's own factor, `Cs = U'U` in place, rather than LAPACK's,
        # which took its process-wide lock once per block per unit on every
        # worker (`small_linalg.jl`); and `Q = L U'` by hand rather than a BLAS
        # product, which locks the same way.
        Cs = _laplace_scratch_matrix!(laplace, Float64, k, k, :sg_Cs)
        @inbounds for y in 1:k, x in 1:k
            Cs[x, y] = (Cdiag[b][x, y] + Cdiag[b][y, x]) / 2
        end
        _ctsem_cholesky!(Cs, k) || (_LAPLACE_DIAG[] = 11; return false)
        Q = _laplace_scratch_matrix!(laplace, Float64, size(L, 1), k, :sg_Q)
        @inbounds for j in 1:k, p in 1:size(L, 1)
            acc = 0.0
            for z in j:k
                acc += L[p, z] * Cs[j, z]
            end
            Q[p, j] = acc
        end
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
            isactive[block.offset + j] || continue
            dir = scatter!(dir1, l, view(Q, :, j))
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
                    # C[a,b] is zero in an inactive row or column.
                    (isactive[block.offset + i] && isactive[blocks[a].offset + j]) ||
                        continue
                    accumulate(sweep_at(subsA, scatter!(dir1, l, view(L, :, i)),
                        scatter!(dir2, la, view(Ls[la], :, j)), 2), subsA,
                        2 * Cab[j, i]) || return false
                end
            else
                # `Vcross`, not `V`: the explicit terms below use a `V` of
                # their own at this function's scope, and this assignment would
                # land on it. Harmless only because that one is reallocated
                # before it is read -- a reordering, or lifting this loop into
                # a closure, makes it a silent wrong answer.
                # Vcross = Ls[la] C[a,b] L', through `T1 = C[a,b] L'`, by the
                # engine's kernels rather than two locking BLAS products.
                T1 = _laplace_scratch_matrix!(laplace, Float64, na, kparl, :sg_T1)
                @inbounds for q in 1:kparl, j in 1:na
                    acc = 0.0
                    for i in 1:nb
                        acc += Ccoup[b][t][i, j] * L[q, i]
                    end
                    T1[j, q] = acc
                end
                Vcross = _laplace_scratch_matrix!(laplace, Float64,
                    size(Ls[la], 1), kparl, :sg_Vcross)
                _ctsem_mul!(Vcross, Ls[la], T1)
                for q in 1:kparl
                    # The basis direction `e_q` of this level, scattered.
                    fill!(dir1, 0.0)
                    dir1[spec.levels[l].re_index[q]] = 1.0
                    accumulate(sweep_at(subsA, dir1,
                        scatter!(dir2, la, view(Vcross, :, q)), 2), subsA,
                        2.0) || return false
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

    # Deviations the filter integrates (`CTSEMHybrid`): the joint logdet holds
    # each member's `logdet(M_m)` beside the profile curvature's, and enters
    # the term as `-1/2` of it, exactly as `-1/2 logdet` of the curvature does
    # through `psi / 2`. So its gradient joins `psi`'s, twice over: half of it
    # is added below as `dpsi/dtheta` is, and its `u` part goes into `s`,
    # carrying the mode's own dependence. The gradient of `-logdet(M_m) / 2` is
    # the full pass's less the profile's (`llvm`).
    if _hybrid_on(laplace.objective.params.hybrid)
        full = _laplace_hybrid_full_gradients(laplace, U, values, Ls, uhat)
        full === nothing && return false
        @inbounds for m in 1:nmem, t in 1:npar
            Pm[t, m] += 2 * (full[t, m] - llvm[t, m])
        end
    end

    # dpsi/du_b runs through this block's members only, then s = C dpsi/du.
    gradu = fill!(_laplace_scratch_vector!(laplace, Float64, d, :sg_gradu), 0.0)
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
        acc = fill!(_laplace_scratch_vector!(laplace, Float64, kpar, :sg_acc), 0.0)
        for m in block.members, p in 1:kpar
            acc[p] += Pm[rho[p], m]
        end
        @inbounds for q in 1:k
            contribution = 0.0
            for p in 1:kpar
                contribution += Ls[l][p, q] * acc[p]
            end
            gradu[block.offset + q] += contribution
        end
    end
    s = _laplace_block_solve!(_laplace_scratch_vector!(laplace, Float64, d, :sg_s),
        factors, elim, blocks, gradu)
    all(isfinite, s) || return false

    # s' B: one sweep per block, along that block's share of L s. Summed over
    # blocks this is one directional derivative per member along the total shift
    # its own block and its ancestors impose, which is what the product needs.
    sB! = function (mine)
        # See the note in `assemble!`: without this the tasks share `pass`.
        local ownedB, b, block, k, subsB, sb, dir, pass, c, m, t
        ownedB = fill!(_laplace_scratch_vector!(laplace, Bool, nmem, :sg_ownedB), false)
        subsB = _laplace_scratch_vector!(laplace, Int, 0, :sg_subsB)
        dir = _laplace_scratch_vector!(laplace, Float64, npar, :sg_dirB)
        @inbounds for m in mine; ownedB[m] = true; end
        for (b, block) in enumerate(blocks)
            k = block.size
            k == 0 && continue
            empty!(subsB)
            for m in block.members
                ownedB[m] && push!(subsB, m)
            end
            isempty(subsB) && continue
            sb = view(s, (block.offset + 1):(block.offset + k))
            scatter_loaded!(dir, block.level, Ls[block.level], sb)
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
    levelsweeps = needV ? Vector{Vector{Matrix{Float64}}}(undef, nlev) :
        _LAPLACE_NO_LEVELSWEEPS
    if needV
        for l in 1:nlev
            rl = nlatent(spec.levels[l])
            levelsweeps[l] = Vector{Matrix{Float64}}(undef, rl)
            for q in 1:rl
                # Only the members for whom column `q` is active: for the
                # others the direction moves nothing they can see, and their
                # columns are zero.
                if act === nothing
                    pass = sweep(1:nmem, scatter!(
                        _laplace_scratch_vector!(laplace, Float64, npar, :sg_dirV),
                        l, view(Ls[l], :, q)), zerodir, 1)
                    pass.ok || return false
                    all(isfinite, pass.d1c) || return false
                    # A copy, because the sweep's output buffer is reused by
                    # the next `q` and these are read after the loop.
                    levelsweeps[l][q] = copy(pass.d1c)
                    continue
                end
                seeing = Int[m for m in 1:nmem if isactive[units.offsets[U][m][l] + q]]
                if isempty(seeing)
                    levelsweeps[l][q] = zeros(Float64, npar, nmem)
                    continue
                end
                pass = sweep(seeing, scatter!(
                    _laplace_scratch_vector!(laplace, Float64, npar, :sg_dirV),
                    l, view(Ls[l], :, q)), zerodir, 1)
                pass.ok || return false
                all(isfinite, pass.d1c) || return false
                if length(seeing) == nmem
                    levelsweeps[l][q] = copy(pass.d1c)
                else
                    swept = zeros(Float64, npar, nmem)
                    @inbounds for (c, m) in enumerate(seeing), t in 1:npar
                        swept[t, m] = pass.d1c[t, c]
                    end
                    levelsweeps[l][q] = swept
                end
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

    if needV
        V = Vector{Matrix{Float64}}(undef, length(blocks))
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
    #
    # Block by block, and everything that does not depend on which population
    # parameter moves -- `G` and `F`, and the one matrix each trace contracts
    # `X` against -- is formed once per block, in the slot's scratch. Formed
    # inside the parameter loop, as it was, each was a few arrays and a BLAS
    # product per parameter per unit: a quarter of what a Laplace gradient
    # allocated, on a one-effect model.
    positions = _laplace_scratch_vector!(laplace, Int, 0, :sg_positions)
    for (b, block) in enumerate(blocks)
        l = block.level
        k = block.size
        if k > 0 && !isempty(dL[l])
            rho = spec.levels[l].re_index
            kpar = length(rho)
            G = fill!(_laplace_scratch_vector!(laplace, Float64, kpar, :sg_G), 0.0)
            F = fill!(_laplace_scratch_vector!(laplace, Float64, kpar, :sg_F), 0.0)
            for m in block.members, p in 1:kpar
                G[p] += llvm[rho[p], m] + (Pm[rho[p], m] + Bsm[rho[p], m]) / 2
                F[p] += llvm[rho[p], m]
            end
            # tr(C[b,b] A[b,b] X) and, for each ancestor, tr(X' A[b,a] C[a,b])
            # -- this block's own factor moving -- as one contraction
            # `sum(X .* W)`, `W = (C[b,b] A[b,b])' + sum_a A[b,a] C[a,b]'`.
            if !needV
                W = _laplace_scratch_matrix!(laplace, Float64, k, k, :sg_W)
                Cb = Cdiag[b]
                Md = M.diag[b]
                @inbounds for y in 1:k, x in 1:k
                    acc = Cb[x, y]
                    for z in 1:k
                        acc -= Cb[x, z] * Md[z, y]
                    end
                    W[y, x] = acc
                end
                for (tt, a) in enumerate(block.ancestors)
                    Mc = M.coupling[b][tt]
                    Cc = Ccoup[b][tt]
                    @inbounds for x in 1:k, y in 1:k
                        acc = 0.0
                        for z in axes(Mc, 2)
                            acc += Mc[y, z] * Cc[x, z]
                        end
                        W[y, x] -= acc
                    end
                end
            end
            _laplace_level_positions!(positions, spec, l)
            for (t, j) in enumerate(positions)
                # `dL` is parameter-by-block, so both shifts land in parameter
                # space where `G` and `F` live.
                D = dL[l][t]
                total = 0.0
                @inbounds for q in 1:k
                    uq = uhat[block.offset + q]
                    sq = s[block.offset + q]
                    for p in 1:kpar
                        total += D[p, q] * (G[p] * uq + F[p] * sq / 2)
                    end
                end
                if needV
                    Vb = V[b]
                    @inbounds for x in axes(Vb, 1), y in axes(Vb, 2)
                        total += Vb[x, y] * D[x, y]
                    end
                else
                    X = dLsolved[l][t]
                    @inbounds for x in 1:k, y in 1:k
                        total += X[y, x] * W[y, x]
                    end
                end
                out[j] += total
            end
        end
        # tr(C[a,b] A[b,a] X): an ancestor's factor moving, in the ancestor's
        # level's parameters, against `Z = -C[a,b] A[b,a]`.
        needV && continue
        for (tt, a) in enumerate(block.ancestors)
            la = blocks[a].level
            isempty(dL[la]) && continue
            ka = blocks[a].size
            Mc = M.coupling[b][tt]
            Cc = Ccoup[b][tt]
            Z = _laplace_scratch_matrix!(laplace, Float64, ka, ka, :sg_Z)
            @inbounds for y in 1:ka, x in 1:ka
                acc = 0.0
                for z in 1:k
                    acc += Cc[z, x] * Mc[z, y]
                end
                Z[x, y] = -acc
            end
            _laplace_level_positions!(positions, spec, la)
            for (t, j) in enumerate(positions)
                X = dLsolved[la][t]
                total = 0.0
                @inbounds for x in 1:ka, y in 1:ka
                    total += Z[x, y] * X[y, x]
                end
                out[j] += total
            end
        end
    end
    return true
end

"""
    _laplace_hybrid_full_gradients(laplace, U, values, Ls, u)

Each of unit `U`'s members' raw gradient at `u` with the deviations its filter
integrates (`CTSEMHybrid`) integrated -- `full`, logdet and all -- as an
`npar x members` matrix, or `nothing` if a member cannot be evaluated. In
parallel over members.
"""
function _laplace_hybrid_full_gradients(laplace::CTSEMLaplaceObjective, U::Integer,
    values::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64})
    units = laplace.units
    members = units.members[U]
    npar = length(values)
    # The calling slot's: read by the caller at once, before its next unit.
    out = _laplace_scratch_matrix!(laplace, Float64, npar, length(members), :hybrid_full)
    ok = _laplace_partition(laplace, length(members)) do mine, _w
        local aws, grad, shift, m, shifted, loglik
        aws = _laplace_workspace!(laplace, Float64, npar)
        grad = _laplace_scratch_vector!(laplace, Float64, npar, :hybrid_grad)
        shift = _laplace_scratch_vector!(laplace, Float64, npar, :hybrid_shift)
        for m in mine
            shifted = _laplace_member_values!(shift, values, laplace.spec, Ls, u,
                units.offsets[U][m])
            loglik = _laplace_subject_value_gradient!(grad,
                laplace.objective.subject_objectives[members[m]], aws, shifted,
                nothing, _HYBRID_FULL)
            (isfinite(loglik) && all(isfinite, grad)) || return false
            @inbounds for t in 1:npar
                out[t, m] = grad[t]
            end
        end
        return true
    end
    return ok ? out : nothing
end

"""Raw positions of one level's free population parameters, scales then
correlations. A stated entry (index 0) is not a parameter and has no position."""
_laplace_level_positions(spec::CTSEMLaplaceSpec, l::Integer) =
    _laplace_level_positions!(Int[], spec, l)

"""The same, into `dest`, for a caller that runs once per unit."""
function _laplace_level_positions!(dest::Vector{Int}, spec::CTSEMLaplaceSpec, l::Integer)
    level = spec.levels[l]
    empty!(dest)
    if isreducedrank(level)
        append!(dest, level.load_index)
        for i in level.unique_index
            i > 0 && push!(dest, i)
        end
    else
        for i in level.sd_index
            i != 0 && push!(dest, i)
        end
        for i in level.cor_index
            i != 0 && push!(dest, i)
        end
    end
    return dest
end

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
    return [_laplace_level_chol_derivatives(values, spec, l) for l in eachindex(spec.levels)]
end

"""
    _laplace_level_chol_derivatives(values, spec, l)

Level `l`'s entry of the above alone, for a caller that steps one level's
population parameters (SAEM's centred step).
"""
function _laplace_level_chol_derivatives(values::AbstractVector{Float64},
    spec::CTSEMLaplaceSpec, l::Integer)
    level = spec.levels[l]
    positions = _laplace_level_positions(spec, l)
    k = length(level.re_index)
    r = nlatent(level)
    (isempty(positions) || k == 0 || r == 0) && return Matrix{Float64}[]
    chol_of = function (p)
        v = convert(Vector{eltype(p)}, values)
        @inbounds for (slot, position) in enumerate(positions)
            v[position] = p[slot]
        end
        return vec(_laplace_popchol(v, level))
    end
    J = _ctsem_jacobian(chol_of, values[positions])
    return [Matrix{Float64}(reshape(collect(view(J, :, t)), k, r))
            for t in eachindex(positions)]
end

"""Row `p` of `D * u[base .+ (1:r)]`: a level's `dL` (or second derivative)
against its effects in `u`, which is how one of its effects' shifts moves with
that population parameter at fixed `u`."""
@inline function _laplace_dLu(D::AbstractMatrix, p::Integer, u::AbstractVector,
    base::Integer, r::Integer)
    inner = zero(promote_type(eltype(D), eltype(u)))
    @inbounds for q in 1:r
        inner += D[p, q] * u[base + q]
    end
    return inner
end

"""
    _laplace_chol_chain!(out, grad, spec, dL, positions, u, offsets)

Add into `out` the part of a member's gradient in `theta` that reaches the
population parameters through `L(theta) u` at fixed `u`: for each level and
each of its parameters (`positions`, `dL` in that order),
`sum_p grad[re_index[p]] (dL u)_p`, `grad` being the member's gradient in its
shifted parameters and `offsets` its places in `u`. The one statement of that
chain rule: at the mode (`_laplace_floored_unit_gradient!`), at a quadrature
node (`_continuation_member!`) and in the continuation's exact Hessian.
"""
function _laplace_chol_chain!(out::AbstractVector, grad::AbstractVector,
    spec::CTSEMLaplaceSpec, dL, positions, u::AbstractVector, offsets)
    @inbounds for l in eachindex(spec.levels)
        level = spec.levels[l]
        k = nrandomeffects(level)
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        pos = positions[l]
        isempty(pos) && continue
        base = offsets[l]
        for t in eachindex(pos)
            dLt = dL[l][t]
            acc = zero(eltype(out))
            for pp in 1:k
                gp = grad[level.re_index[pp]]
                iszero(gp) && continue
                acc += gp * _laplace_dLu(dLt, pp, u, base, r)
            end
            out[pos[t]] += acc
        end
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

"""
    CTSEMCoreset(units, weights)

A weighted subset of a Laplace objective's units: an evaluation given one
solves only these units' inner problems and returns the log prior plus
`sum(weights .* terms)` over them, in place of every unit's term. Importance
sampling uses it where a few units' terms, refitted as the draws move, follow
the full sum closely (R's `.ctCoresetLpg()`); `units` index
`laplace.units.members`, and a weight replaces the unit's sampling weight.
"""
struct CTSEMCoreset
    units::Vector{Int}
    weights::Vector{Float64}
end

"""A coreset over `laplace`'s units, checked: distinct units in range, finite
non-negative weights, one per unit."""
function _laplace_coreset(laplace::CTSEMLaplaceObjective, units::AbstractVector,
    weights::AbstractVector)
    nunits = length(laplace.units.members)
    u = Int[Int(x) for x in units]
    w = Float64[Float64(x) for x in weights]
    length(u) == length(w) || throw(ArgumentError(
        "a coreset needs one weight per unit: $(length(u)) units, $(length(w)) weights"))
    isempty(u) && throw(ArgumentError("a coreset needs at least one unit"))
    all(U -> 1 <= U <= nunits, u) || throw(ArgumentError(
        "coreset units must lie in 1:$nunits"))
    allunique(u) || throw(ArgumentError("coreset units must be distinct"))
    all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError(
        "coreset weights must be finite and non-negative"))
    return CTSEMCoreset(u, w)
end

"""Each unit's weight in the evaluation: the coreset's, zero off it; empty
without a coreset, when `_laplace_unit_weight` decides."""
function _laplace_coreset_weights(coreset::Union{Nothing,CTSEMCoreset}, nunits::Integer)
    coreset === nothing && return Float64[]
    w = zeros(Float64, nunits)
    @inbounds for k in eachindex(coreset.units)
        w[coreset.units[k]] = coreset.weights[k]
    end
    return w
end

"""NaN for every unit, and its subjects, that the coreset did not evaluate."""
function _laplace_coreset_mask!(unit_loglik::Vector{Float64},
    subject_loglik::Vector{Float64}, laplace::CTSEMLaplaceObjective,
    coreset::CTSEMCoreset)
    evaluated = falses(length(unit_loglik))
    evaluated[coreset.units] .= true
    @inbounds for U in eachindex(unit_loglik)
        evaluated[U] && continue
        unit_loglik[U] = NaN
        for i in laplace.units.members[U]
            subject_loglik[i] = NaN
        end
    end
    return nothing
end

"""Whether every evaluated unit's inner mode was found."""
_laplace_converged(laplace::CTSEMLaplaceObjective, coreset::Union{Nothing,CTSEMCoreset}) =
    coreset === nothing ? all(laplace.inner_converged) :
    all(U -> laplace.inner_converged[U], coreset.units)

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
    gradient::Bool=true, contributions::Bool=false, nested_gradient::Bool=false,
    scores::Union{Nothing,Matrix{Float64}}=nothing,
    coreset::Union{Nothing,CTSEMCoreset}=nothing)
    n = length(values)
    theta = _laplace_extend(laplace, collect(Float64, values))
    length(theta) == n && return _laplace_evaluate(laplace, theta; gradient=gradient,
        contributions=contributions, nested_gradient=nested_gradient, scores=scores,
        coreset=coreset)
    # With stated population values (`CTSEMAnchors`) the gradient and the
    # scores are formed over the extended vector and returned over the
    # caller's.
    full = scores === nothing ? nothing : zeros(Float64, size(scores, 1), length(theta))
    r = _laplace_evaluate(laplace, theta; gradient=gradient,
        contributions=contributions, nested_gradient=nested_gradient, scores=full,
        coreset=coreset)
    full === nothing || copyto!(scores, view(full, :, 1:n))
    return _laplace_free_result(r, n)
end

"""`ctsem_laplace_evaluate` at a vector already extended (`_laplace_extend`)."""
function _laplace_evaluate(laplace::CTSEMLaplaceObjective, theta::Vector{Float64};
    gradient::Bool=true, contributions::Bool=false, nested_gradient::Bool=false,
    scores::Union{Nothing,Matrix{Float64}}=nothing,
    coreset::Union{Nothing,CTSEMCoreset}=nothing)
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
    #
    # A coreset (`CTSEMCoreset`) is the same region over its own units only:
    # the rest are neither solved nor summed, so their entries come back NaN.
    ranges = coreset === nothing ? [U:U for U in 1:nunits] : [U:U for U in coreset.units]
    nchunks = length(ranges)
    coreset_w = _laplace_coreset_weights(coreset, nunits)
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
            _laplace_store_unit_primal!(laplace, U, theta, Ls, aws, primal_curvature,
                primal_matrices, primal_gated, unit_loglik, chunk_ok, chunk_bad) ||
                return nothing
            term = unit_loglik[U]
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
            converged=_laplace_converged(laplace, coreset))
    end
    value = (coreset === nothing ? _laplace_weighted_total(laplace, unit_loglik) :
        sum(coreset.weights[k] * unit_loglik[coreset.units[k]]
            for k in eachindex(coreset.units))) +
        _ctsem_log_prior(laplace.objective, theta)
    coreset === nothing || _laplace_coreset_mask!(unit_loglik, subject_loglik,
        laplace, coreset)
    laplace.gated_units = count(primal_gated)

    gradient || return (value=value, gradient=nothing,
        subject_loglik=subject_loglik, unit_loglik=unit_loglik,
        converged=_laplace_converged(laplace, coreset))

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
        dLsolved = _laplace_level_chol_solves(laplace.spec, Ls, dLlevels)
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
        # With `scores`, each unit's gradient is formed in a row buffer of the
        # slot's and then added both to the slot's sum and to the unit's row:
        # the units' own gradients are the score rows
        # (`ctsem_subject_gradients`). Without, the buffer is the sum itself.
        # With sampling weights too: each unit's gradient is then its weight
        # times the row.
        weighted = !isempty(laplace.objective.subject_weights)
        separate = scores !== nothing || weighted || coreset !== nothing
        rows = separate ? [zeros(Float64, length(theta)) for _ in 1:nslot] : partials
        fill!(chunk_ok, true)
        run_gradient = function (c)
            @inbounds for U in ranges[c]
                separate && fill!(rows[_laplace_slot()], 0.0)
                # One pass, so a unit handled early can `break` to the row below.
                for _ in 1:1
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
                        rows[_laplace_slot()] .+= gg
                        break
                    end
                    # A floored unit is differentiated whole; see
                    # `_laplace_floored_unit_gradient!`.
                    if laplace.logdet_floored[U]
                        if !_laplace_floored_unit_gradient!(
                                rows[_laplace_slot()], laplace, U,
                                theta, Ls, dLlevels)
                            chunk_ok[_laplace_slot()] = false
                            return nothing
                        end
                        break
                    end
                    factors, elim = primal_curvature[U]
                    local bg
                    bg = _laplace_mark()
                    if !_laplace_seeded_unit_gradient!(rows[_laplace_slot()],
                            laplace, U, theta,
                            Ls, dLlevels, dLsolved, primal_matrices[U], factors, elim)
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
                    if !all(isfinite, rows[_laplace_slot()])
                        chunk_ok[_laplace_slot()] = false
                        return nothing
                    end
                end
                if separate
                    w = coreset === nothing ? _laplace_unit_weight(laplace, U) : coreset_w[U]
                    partials[_laplace_slot()] .+= w .* rows[_laplace_slot()]
                    scores === nothing || (scores[U, :] .= w .* rows[_laplace_slot()])
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
            grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_curvature;
                unitweights=coreset_w)
            _CTSEM_LAPLACE_FALLBACKS[] += 1
            scores === nothing || fill!(scores, NaN)
        end
    else
        grad .= _laplace_nested_gradient(laplace, theta, Ls, primal_curvature;
                unitweights=coreset_w)
        scores === nothing || fill!(scores, NaN)
    end
    return (value=value, gradient=grad, subject_loglik=subject_loglik,
        unit_loglik=unit_loglik, converged=_laplace_converged(laplace, coreset))
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
        _laplace_chol_chain!(out, grad, spec, dL, levelpositions, u, offsets)
    end
    return true
end

"""Unit `U`'s sampling weight: its subjects', which the R side requires to be
one value within a unit; one when the objective has none."""
@inline function _laplace_unit_weight(laplace::CTSEMLaplaceObjective, U::Integer)
    weights = laplace.objective.subject_weights
    isempty(weights) && return 1.0
    members = laplace.units.members[U]
    return isempty(members) ? 1.0 : @inbounds(weights[first(members)])
end

"""The units' terms summed with their sampling weights."""
function _laplace_weighted_total(laplace::CTSEMLaplaceObjective, terms::AbstractVector)
    isempty(laplace.objective.subject_weights) && return sum(terms)
    total = zero(eltype(terms))
    for U in eachindex(terms)
        w = _laplace_unit_weight(laplace, U)
        w == 0 || (total += w * terms[U])
    end
    return total
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
    curvature::AbstractVector; unitweights::Vector{Float64}=Float64[])
    nunits = length(laplace.units.members)
    total_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Lsd = _laplace_popchols(x, laplace.spec)
        accumulated = zero(S)
        for U in 1:nunits
            w = isempty(unitweights) ? _laplace_unit_weight(laplace, U) :
                unitweights[U]
            w == 0 && continue
            uhat = laplace.modes[U]
            ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, curvature[U], wsd)
            accumulated += w * _laplace_unit_term(laplace, U, x, Lsd, ud, wsd)
        end
        return accumulated + _ctsem_log_prior(laplace.objective, x)
    end
    return _ctsem_gradient(total_of, theta)
end

"""
    _laplace_shifted_values(laplace, theta, Ls, effect_of)

Every subject's raw vector as a filter is handed it: the population vector
shifted by the random effects of each unit the subject belongs to, and *not*
by its TI-predictor effects. `effect_of(U)` returns unit `U`'s latent vector.

This is the vector the Laplace term itself filters -- `_laplace_member_values`,
then the filter, which adds TI-predictor effects on its first line -- so it is
the one every other filtering consumer on this route must be handed as well:
`ctsem_kalman`, `ctsem_generate`, `ctsem_generate_states`, the effect draws
behind leave-one-row-out, the substep mesh. A filter handed the vector with TI
effects already in it applies them a second time, and nothing fails: the
likelihood, the estimates and the gradients never pass through here, only the
trajectories, residuals and generated data do.

So a subject's vector is built here and nowhere else. The one step taken after
it is `_laplace_add_ti_effects`, for reporting what a subject's parameters are.

A subject no unit reached is `NaN`, so it would filter as not finite rather
than at whatever the allocation held.
"""
function _laplace_shifted_values(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, Ls::Vector{<:AbstractMatrix}, effect_of)
    units = laplace.units
    out = fill(NaN, length(laplace.objective.subject_objectives), length(theta))
    for U in eachindex(units.members)
        u = effect_of(U)
        for (m, i) in enumerate(units.members[U])
            out[i, :] = _laplace_member_values(theta, laplace.spec, Ls, u,
                units.offsets[U][m])
        end
    end
    return out
end

"""
    _laplace_add_ti_effects(laplace, shifted)

`shifted`, one row per subject, with each subject's TI-predictor effects added
by the filter's own `_materialize_subject_values!`: what that subject's raw
parameters are, for pushing through the transforms and reporting. Never for
handing to a filter, which adds them itself; see `_laplace_shifted_values`.
"""
function _laplace_add_ti_effects(laplace::CTSEMLaplaceObjective,
    shifted::AbstractMatrix)
    out = Matrix{Float64}(undef, size(shifted))
    buffer = Float64[]
    for (i, subject) in enumerate(laplace.objective.subject_objectives)
        row = view(shifted, i, :)
        _materialize_subject_values!(buffer, row, subject.params,
            _ctsem_tipred_vector(subject.tipreds, row))
        # Without the TI predictor tail (`EKFParameters.ti_tail`).
        out[i, :] = view(buffer, 1:size(out, 2))
    end
    return out
end

"""
    ctsem_laplace_subject_values(laplace, values; from_level=1, ti_effects=true)

Each subject's own raw parameter vector at the current inner modes: the
population vector shifted by that subject's random effects, and then by its
TI-predictor effects.

Both shifts come from the code that already applies them during fitting --
`_laplace_member_values` and the engine's own `_materialize_subject_values!`
-- rather than being reconstructed by the caller. That matters because the
caller is the R side's subject-parameter reporting, and a second, slightly
different copy of "what parameters does this subject have" is exactly the kind
of divergence that shows up as a summary disagreeing with the fit.

`ti_effects = false` stops before the second shift and returns what a filter is
handed instead, which is what the `subject_values` keyword of `ctsem_kalman`,
`ctsem_generate` and `ctsem_generate_states` takes: the filter adds each
subject's TI-predictor effects itself, so a vector that already carries them
has them twice.

Returns `nsubjects x length(values)`; push a row through the model's transforms
to get that subject's parameter matrices.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; from_level::Integer=1, ti_effects::Bool=true)
    theta = _laplace_extend(laplace, collect(Float64, values))
    Ls = _laplace_popchols(theta, laplace.spec)
    # Every unit's mode on the pool, as an evaluation solves them: each writes
    # only its own `laplace.modes[U]`, and the slot its task holds. Solved one
    # unit after another, the plug-in Kalman passes an adequacy check makes,
    # two a parameter, ran on one core -- 28 s each on the Fisher example's
    # six-latent, 11-effect model, 97 minutes in all.
    _laplace_parallel(laplace, eachindex(laplace.units.members)) do U
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        return true
    end
    shifted = _laplace_shifted_values(laplace, theta, Ls, U ->
        _laplace_restrict_levels(laplace, U, laplace.modes[U], from_level))
    # Deviations the subject filter integrates (`CTSEMHybrid`), at their
    # conditional modes given each subject's whole record and its effects'
    # modes: the subject level's own, so only where that level is kept.
    from_level <= 1 && _hybrid_add_modes!(shifted, laplace, theta, Ls)
    return ti_effects ? _laplace_add_ti_effects(laplace, shifted) : shifted
end

"""
    ctsem_laplace_subject_values(laplace, values, effects; ti_effects=true, from_level=1)

Each subject's own raw parameter vector at the random effects it is *given*,
rather than at the conditional modes solved from the data.

The method above estimates each unit's effects as a mode. That is the right
answer for a fit that integrated them out, and the wrong one for a fit that
sampled them: a sampler already has a draw of every effect, and pairing draw
`s` of the parameters with draw `s` of the effects is what makes a posterior
predictive integrate over the random effects instead of conditioning on a
point estimate of them.

`effects` is one draw's natural deviations as a sampled fit reports them: each
block's `S u` on its level's parameters, in `ctsem_laplace_deviation_layout`'s
layout. Deviations rather than the sampler's coordinates `u`, so a draw means
the same thing whichever square root of the population covariance it was made
under.

Everything after the slicing is shared with the mode method, including what
`ti_effects` and `from_level` mean, so the two cannot drift in how a shifted
parameter vector is built.

Returns `nsubjects x length(values)`.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, effects::AbstractVector; ti_effects::Bool=true,
    from_level::Integer=1)
    theta = _laplace_extend(laplace, collect(Float64, values))
    shifted = _laplace_shifted_values_from_deviations(laplace, theta, effects;
        from_level=from_level)
    return ti_effects ? _laplace_add_ti_effects(laplace, shifted) : shifted
end

"""
    ctsem_laplace_subject_values(laplace, draws, effects; ti_effects=true, from_level=1)

The method above for a whole sample at once: `draws[s, :]` is one draw of the
parameters and `effects[s, :]` the natural deviations drawn with it, so each
subject's vector at draw `s` is built from the pair. Batched because the caller
is the summary, which has a sample rather than a point, and one call a draw
would cost a bridge round trip each.

Returns `ndraws x nsubjects x npar`, as the linearised method does.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    draws::AbstractMatrix, effects::AbstractMatrix; ti_effects::Bool=true,
    from_level::Integer=1)
    size(effects, 1) == size(draws, 1) || throw(DimensionMismatch(string(
        "effects has ", size(effects, 1), " draws for ", size(draws, 1),
        " parameter draws")))
    draws = _laplace_extend_rows(laplace, draws)
    ndraws = size(draws, 1)
    nsubjects = length(laplace.objective.subject_objectives)
    out = zeros(Float64, ndraws, nsubjects, size(draws, 2))
    for s in 1:ndraws
        out[s, :, :] = ctsem_laplace_subject_values(laplace,
            collect(Float64, view(draws, s, :)), view(effects, s, :);
            ti_effects=ti_effects, from_level=from_level)
    end
    return out
end

"""
    _laplace_filter_values(laplace, values; from_level, effects, subject_values)

The per-subject matrix a Laplace method hands the wrapped objective's filter:
`subject_values` as the caller supplied it, else the vectors at `effects`, else
at the conditional modes with the levels inside `from_level` left out. Always
before TI-predictor effects (`_laplace_shifted_values`), including a supplied
`subject_values`, which is why the R side asks for those with
`ti_effects = false`.

One function for all three filtering entry points, so that what a filter is
handed on this route is decided once.
"""
function _laplace_filter_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; from_level::Integer=1,
    effects::Union{Nothing,AbstractVector}=nothing,
    subject_values::Union{Nothing,AbstractMatrix}=nothing)
    if subject_values !== nothing
        nsubjects = length(laplace.objective.subject_objectives)
        size(subject_values, 1) == nsubjects || throw(ArgumentError(
            "subject_values has $(size(subject_values, 1)) rows for " *
            "$nsubjects subjects."))
        # A subject's vector carries any stated population values
        # (`CTSEMAnchors`) after the free parameters, as every vector this
        # route builds does.
        width = length(values) + _ctsem_anchor_append(laplace.objective.params,
            length(values))
        size(subject_values, 2) == width || throw(ArgumentError(
            "subject_values has $(size(subject_values, 2)) columns for " *
            "$width parameters."))
        return subject_values
    end
    effects === nothing ||
        return ctsem_laplace_subject_values(laplace, values, effects; ti_effects=false)
    return ctsem_laplace_subject_values(laplace, values; from_level=from_level,
        ti_effects=false)
end

export ctsem_laplace_subject_values

"""
    ctsem_auto_substeps(laplace, values; tol, max_substeps, passes, floor_rule, fallback)

The substep mesh for a Laplace objective, with each subject measured at its
conditional random-effect modes at `values`.

A mesh is measured along the state path the filter takes, and on this route
that path depends on the random effects: a subject is filtered at the
population vector shifted by `L u`. Every evaluation of the Laplace term takes
its value and its curvature at the inner mode, so that is where the
linearisation has to be accurate, and that is where the mesh is measured. The
effects at zero would measure the population's system instead, which for a
subject whose effect is not small is a different path -- with a varying
manifest mean, the filter moves the latent state to absorb the subject's
offset -- and not one the fit ever evaluates.

The modes are solved under the mesh the objective already carries, since no
other exists before one is chosen; the R side chooses again at the optimum,
where that mesh is the fit's own. A unit whose inner problem is not finite is
measured at the origin, the population vector, and if the filter is not finite
there either, `finite` says so as it does for any objective. The solve leaves
its modes in `laplace.modes`, as `ctsem_laplace_subject_values` does.

The keywords, and what is returned, are the `CTSEMObjective` method's. Each
subject's shifted vector is one row of the matrix that method takes, before
TI-predictor effects, which the filter adds exactly as it does in the fit.
"""
function ctsem_auto_substeps(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    kwargs...)
    theta = _laplace_extend(laplace, collect(Float64, values))
    _laplace_check_indices(laplace, length(theta))
    Ls = _laplace_popchols(theta, laplace.spec)
    shifted = _laplace_shifted_values(laplace, theta, Ls, function (U)
        solved = _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        return isfinite(solved.value) ? solved.u :
            zeros(Float64, laplace.units.dims[U])
    end)
    return ctsem_auto_substeps(laplace.objective, shifted; kwargs...)
end

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
random effects as natural deviations, in `ctsem_laplace_deviation_layout`'s
layout, from which each subject's own parameter vector is built; `subject_values` supplies those vectors directly
instead, before TI-predictor effects as `ctsem_kalman` takes them. Given
neither, each subject is put at its conditional mode, matching what
`ctsem_generate` does.

The states are drawn afresh here -- from the process, at that subject's own
parameters -- and never carried over from anything a fit sampled. What comes
from a fit is the individual differences; the trajectory is regenerated
conditional on the subject-specific model.
"""
function ctsem_generate_states(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, z::AbstractVector, base::AbstractMatrix;
    effects::Union{Nothing,AbstractVector}=nothing,
    subject_values::Union{Nothing,AbstractMatrix}=nothing, kwargs...)
    persubject = _laplace_filter_values(laplace, values; effects=effects,
        subject_values=subject_values)
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

Those rows are *before* TI-predictor effects --
`ctsem_laplace_subject_values(fitted, values; ti_effects=false)` -- because the
filter adds each subject's effects itself, from this objective's predictors,
exactly as it does in the fit.
"""
function ctsem_kalman(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    from_level::Integer=1, subject_matrices::Bool=true, fields=String[], compact::Bool=false,
    subject_values::Union{Nothing,AbstractMatrix}=nothing)
    persubject = _laplace_filter_values(laplace, values; from_level=from_level,
        subject_values=subject_values)
    return ctsem_kalman(laplace.objective, persubject;
        subject_matrices=subject_matrices, fields=fields, compact=compact)
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

With deviations each member's filter integrates (`CTSEMHybrid`), a unit's rows
are its effects and then, member by member, one per hybrid item in the order
`_laplace_effect_posterior` gives -- each the item's conditional mode, which
moves with the parameters directly and through the effects' mode -- and an
item the member never observes is a row of zeros, its mode the prior's.
"""
function ctsem_laplace_mode_jacobian(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    n = length(values)
    J = _laplace_mode_jacobian(laplace,
        _laplace_extend(laplace, collect(Float64, values))).jacobians
    return [size(j, 2) == n ? j : j[:, 1:n] for j in J]
end

"""
    _laplace_mode_jacobian(laplace, theta)

`ctsem_laplace_mode_jacobian`'s Jacobians, and each unit's mode in the same
coordinates (`modes`) at `theta`.
"""
function _laplace_mode_jacobian(laplace::CTSEMLaplaceObjective, theta::Vector{Float64})
    _laplace_check_indices(laplace, length(theta))
    units = laplace.units
    nunits = length(units.members)
    Ls = _laplace_popchols(theta, laplace.spec)
    hessians = _laplace_primal_curvature(laplace, theta, Ls)
    hyb = laplace.objective.params.hybrid
    nh = _hybrid_on(hyb) ? _hybrid_items(hyb) : 0
    slot = nh > 0 ? invperm(sortperm(hyb.effect)) : Int[]
    out = Vector{Matrix{Float64}}(undef, nunits)
    modes = Vector{Vector{Float64}}(undef, nunits)
    for U in 1:nunits
        d = units.dims[U]
        nm = length(units.members[U])
        D = d + nh * nm
        if D == 0
            out[U] = zeros(Float64, 0, length(theta))
            modes[U] = Float64[]
            continue
        end
        uhat = laplace.modes[U]
        mode_of = function (x)
            local S, wsd, Lsd, u, w, hook, r
            S = eltype(x)
            wsd = _laplace_workspace!(laplace, S, length(x))
            Lsd = _laplace_popchols(x, laplace.spec)
            u = d == 0 ? S[] :
                _laplace_dual_unit_mode(laplace, U, x, Lsd, uhat, hessians[U], wsd)
            nh == 0 && return u
            w = zeros(S, D)
            w[1:d] .= u
            hook = _HybridModeHook{S}(nm)
            r = _laplace_unit_loglik_gradient(laplace, U, x, Lsd, u, wsd,
                eachindex(units.members[U]), nothing, _HYBRID_PROFILE, hook)
            isfinite(r.value) || return fill(S(NaN), D)
            for m in 1:nm, (c, k) in enumerate(hook.items[m])
                w[d + (m - 1) * nh + slot[k]] = hook.dhat[m][c]
            end
            return w
        end
        out[U] = _ctsem_jacobian(mode_of, theta)
        if nh == 0
            modes[U] = copy(uhat)
        else
            # The effects' mode as solved, and the deviations' conditional modes
            # given it.
            modes[U] = zeros(Float64, D)
            modes[U][1:d] .= uhat
            hook = _hybrid_member_modes(laplace, U, theta, Ls, copy(uhat),
                _laplace_workspace!(laplace, Float64, length(theta)))
            if hook !== nothing
                for m in 1:nm, (c, k) in enumerate(hook.items[m])
                    modes[U][d + (m - 1) * nh + slot[k]] = hook.dhat[m][c]
                end
            end
        end
    end
    return (jacobians=out, modes=modes)
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

`ti_effects` and `from_level` mean what they do for the other methods:
`from_level = l` leaves out the effects of the levels inside `l`, so a draw at
a grouping level is that group's own values. A unit's effects are drawn
jointly over its levels and the inner ones then dropped, which leaves a draw
from the group's marginal posterior; the same `seed` at two `from_level`s
makes the same joint draw, so their difference is one level's share of it.

Deviations each member's filter integrates (`CTSEMHybrid`) are linearised with
the effects (`ctsem_laplace_mode_jacobian`) and enter at each draw's own
variances, as the layer holding them as coordinates would have them.

With `conditional = true` each draw also carries a draw of every unit's
effects from their conditional posterior given its data, `N(0, M^-1)` about
the linearised mode in the unit's coordinates (`_laplace_effect_posterior`,
deviations included): plausible values, whose spread is the population
parameters' uncertainty and each person's own. Without it the spread is the
population parameters' alone, which for a person with little data is far too
narrow. The posterior's covariance is taken at `thetahat`, as the mode's
Jacobian is. `seed` makes those draws reproducible from the R side's seed. A
unit whose curvature cannot be factorized keeps its linearised mode, and the
count of such units is warned about.

Returns `ndraws x nsubjects x npar`.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    draws::AbstractMatrix, thetahat::AbstractVector; ti_effects::Bool=true,
    conditional::Bool=false, seed::Integer=1, from_level::Integer=1)
    # The Jacobian is computed here rather than handed in: it is a vector of
    # matrices, and round-tripping one through the R bridge only to send it
    # straight back costs two marshalling steps for no gain.
    centre = _laplace_extend(laplace, collect(Float64, thetahat))
    draws = _laplace_extend_rows(laplace, draws)
    linear = _laplace_mode_jacobian(laplace, centre)
    ndraws = size(draws, 1)
    npar = size(draws, 2)
    nsubjects = length(laplace.objective.subject_objectives)
    units = laplace.units
    hyb = laplace.objective.params.hybrid
    nh = _hybrid_on(hyb) ? _hybrid_items(hyb) : 0
    items = nh > 0 ? sortperm(hyb.effect) : Int[]
    lv = laplace.spec.levels[1]
    factors = conditional ? _laplace_conditional_factors(laplace, centre, linear, nh) :
        nothing
    rng = MersenneTwister(seed)
    out = zeros(Float64, ndraws, nsubjects, npar)
    for s in 1:ndraws
        theta = collect(Float64, view(draws, s, :))
        Ls = _laplace_popchols(theta, laplace.spec)
        step = theta .- centre
        w = [linear.modes[U] .+ linear.jacobians[U] * step for U in eachindex(units.members)]
        if factors !== nothing
            for U in eachindex(units.members)
                F = factors[U]
                F === nothing && continue
                w[U] .+= F * randn(rng, size(F, 1))
            end
        end
        shifted = _laplace_shifted_values(laplace, theta, Ls, U ->
            _laplace_restrict_levels(laplace, U, w[U][1:units.dims[U]], from_level))
        # The integrated deviations are the subject level's own.
        if nh > 0 && from_level <= 1
            psi = [_hybrid_item_sd(hyb, lv, theta, k) for k in items]
            for U in eachindex(units.members), (m, i) in enumerate(units.members[U])
                at = units.dims[U] + (m - 1) * nh
                for (c, k) in enumerate(items)
                    shifted[i, lv.re_index[hyb.effect[k]]] += psi[c] * w[U][at + c]
                end
            end
        end
        out[s, :, :] = ti_effects ? _laplace_add_ti_effects(laplace, shifted) : shifted
    end
    return out
end

"""
    _laplace_conditional_factors(laplace, theta, linear, nh)

Per unit, the lower Cholesky factor of its effects' conditional posterior
covariance at `theta` (`_laplace_effect_posterior`), in the coordinates of
`_laplace_mode_jacobian`'s modes -- effects, then each member's integrated
deviations -- or `nothing` for a unit with no coordinates, one whose inner
solve did not converge, or one whose curvature could not be factorized; the
last two are warned about once, by count.
"""
function _laplace_conditional_factors(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, linear, nh::Integer)
    units = laplace.units
    Ls = _laplace_popchols(theta, laplace.spec)
    out = Vector{Union{Nothing,Matrix{Float64}}}(nothing, length(units.members))
    failed = 0
    for U in eachindex(units.members)
        d = units.dims[U]
        D = d + nh * length(units.members[U])
        D == 0 && continue
        uhat = collect(Float64, view(linear.modes[U], 1:d))
        post = _laplace_effect_posterior(laplace, U, theta, Ls, uhat, nh)
        if !post.ok || !laplace.inner_converged[U]
            failed += 1
            continue
        end
        F = cholesky(Symmetric((post.cov .+ transpose(post.cov)) ./ 2); check=false)
        if !issuccess(F)
            failed += 1
            continue
        end
        out[U] = Matrix{Float64}(F.L)
    end
    failed > 0 && @warn string(failed, " unit(s) have no conditional draws: their ",
        "curvature could not be factorized or their mode was not found, so their ",
        "draws carry the population parameters' uncertainty only.")
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
        covariance = _laplace_level_covariance(lv, collect(Float64, view(values, s, :)))
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
    theta = _laplace_extend(laplace, collect(Float64, values))
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
    # Deviations the subject filter integrates (`CTSEMHybrid`): reported as the
    # layer holding them as coordinates reports them, from the cross curvature
    # read with the unit's own (`_hybrid_add_mode_rows!`).
    hyb = laplace.objective.params.hybrid
    hybrid = level == 1 && _hybrid_on(hyb)
    aws = hybrid ? _laplace_workspace!(laplace, Float64, length(theta)) : nothing
    for U in eachindex(laplace.units.members)
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        u = laplace.modes[U]
        blocks = laplace.units.blocks[U]
        cross = nothing
        fac = nothing
        if !isempty(u)
            hybrid && (cross = CTSEMHybridCross{Float64}(length(u),
                length(laplace.units.members[U])))
            M = _laplace_unit_curvature(laplace, U, theta, Ls, u; hook=cross)
            fac = _laplace_factor_repaired!(M, blocks)
        end
        # Assigned once: `conditional` captures both (a captured local assigned
        # more than once is a `Core.Box`).
        factors = fac === nothing ? nothing : fac.factors
        coupling = fac === nothing ? nothing : fac.coupling
        hmodes = hybrid ? _hybrid_member_modes(laplace, U, theta, Ls, u, aws) : nothing
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
            # A level with no coordinate of its own still reports deviations
            # the filter integrates.
            (k == 0 || (r == 0 && !hybrid)) && continue
            base = laplace.units.offsets[U][m][level]
            slice = (base + 1):(base + r)
            zi = u[slice]
            z[g, :] = zi
            raw[g, :] = L * zi
            block = conditional(slice)
            zsd[g, :] = sqrt.(max.(diag(block), 0.0))
            rawcov = L * block * transpose(L)
            if hybrid
                _hybrid_add_mode_rows!(raw, rawcov, g, m, slice, laplace, hyb, lv,
                    theta, hmodes, cross, factors, coupling, blocks, length(u))
            end
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
once, tests `M - I` with the block factorization (no fill-in), and finds the
smallest eigenvalue only of the units that fail that test, by
`_laplace_min_eigenvalue` -- densely up to `_LAPLACE_EIGEN_MAXDIM`, by
bisection on the inertia above it. Returns

  * `min_eigenvalue`, per unit: the smallest eigenvalue where the test
    failed, `Inf` where it passed (every eigenvalue above one, not computed),
    and `NaN` where the curvature is not finite;
  * `below_one`, units with an eigenvalue below one (non-finite ones
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
            mins[U] = _laplace_min_eigenvalue(M, blocks, length(u))
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
function _ctsem_optimise_refusal(o::CTSEMLaplaceObjective)
    failed = findall(!, o.inner_converged)
    isempty(failed) && return ""
    return string(": the inner mode did not converge for unit",
        length(failed) == 1 ? " " : "s ", join(failed, ", "), " of ",
        length(o.inner_converged))
end

_ctsem_optimise_progress_extra(o::CTSEMLaplaceObjective) =
    (@sprintf("inner %d/%d", count(o.inner_converged),
        length(o.inner_converged)),)

# A level of correlations can saturate together here, which the per-cell check
# cannot see. See `_laplace_saturated_parameters`.
_ctsem_saturated_for(o::CTSEMLaplaceObjective, minimizer) =
    _laplace_saturated_parameters(o, minimizer)

_ctsem_scores_supported(o::CTSEMLaplaceObjective) = _ctsem_scores_supported(o.objective)

# The Laplace route's scores come from the same evaluation as its gradient.
function _ctsem_optimise_trial_scores(o::CTSEMLaplaceObjective, x, gradient_method,
        limit::Real, log, store::CTSEMScoreStore)
    store.valid = false
    fill!(store.S, 0.0)
    trial = _ctsem_optimise_trial(o, x, true, gradient_method, limit, log;
        scores=store.S)
    if trial.valid && all(isfinite, store.S)
        copyto!(store.x, x); store.valid = true
    end
    return trial
end

"""
A trial point on the laplace route, where finite is not the same as usable.

`converged` is the inner solve's verdict, and a point it failed at is invalid
however finite its value: the Laplace term is defined at the mode, so away from
one the number is not the objective. `ctsem_laplace_evaluate` is called rather
than the generic `ctsem_evaluate`, which deliberately drops that flag.

`gradient_method === :nested` is how `nested_gradient` reaches here through the
shared driver's one gradient-method argument. `coreset`, a `CTSEMCoreset`, is a
coreset stage's trial (`CTSEMCoresetTarget`): its units alone, weighted.
"""
function _ctsem_optimise_trial(o::CTSEMLaplaceObjective, x, want_gradient::Bool,
        gradient_method, limit::Real, log; scores=nothing, coreset=nothing)
    evaluated = try
        ctsem_laplace_evaluate(o, x; gradient=want_gradient,
            nested_gradient=(Symbol(gradient_method) === :nested), scores=scores,
            coreset=coreset)
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

"""
The gradient counterpart of the laplace `_ctsem_probe_value`, refusing a point
whose inner solve did not converge for the same reason.
"""
function _ctsem_probe_value_gradient(o::CTSEMLaplaceObjective, x)
    evaluated = try
        ctsem_laplace_evaluate(o, x; gradient=true)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    evaluated === nothing && return (value=-Inf, gradient=nothing)
    ok = evaluated.converged && isfinite(evaluated.value) &&
        all(isfinite, evaluated.gradient)
    return ok ? (value=Float64(evaluated.value), gradient=evaluated.gradient) :
        (value=-Inf, gradient=nothing)
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

`localised` (the default) differences each column over the units whose terms
read its parameter, and over no others. A unit that does not read a parameter
runs the same computation at both differencing points, so its share of that
column is exactly zero and evaluating it is waste. Which units read which
parameter comes from their scores at `values`: a unit's score in a parameter it
never reads is exactly zero -- an item it never answered, a predictor that is
zero for it, a level it is not in. On a 697-parameter model of 4031 units
(affint's joint fit) the columns then need 12 to 20 per cent of the unit
passes. The result differs from `localised = false` only in rounding.

`localised` also takes a TI-predictor coefficient's column from its target's
(`_laplace_ti_derivable`): where every unit is one subject, a coefficient
`c` on predictor `k` moves a unit only through its target parameter `j`, by
the unit's `x_k`, so `d s_U / d c = x_k d s_U / d theta_j` for every unit's
score `s_U`, and the column is the units' score differences along `j`
weighted by their `x_k`, plus the prior's curvature at `c`. A smooth's
columns then cost nothing past their target's: on a model where a covariate's
spline reaches every varying parameter, most of the columns. The result
differs from differencing along `c` by the difference rule's truncation.

`progress`, when given, is called `(done, total)` after every column with the
gradient evaluations completed so far out of `2n` -- this loop is the one
place a finish or a certification can go quiet for minutes on a wide Laplace
model, so the caller reports through it rather than this function printing
anything itself; see `_progress_fraction` in progress.jl.
"""
function ctsem_laplace_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    step::Real=1e-4, retries::Integer=2, warm::Bool=true, progress=nothing,
    localised::Bool=true)
    H = _laplace_hessian_columns(laplace, values, 1:length(values); step=step,
        retries=retries, warm=warm, progress=progress, localised=localised)
    return (H .+ transpose(H)) ./ 2
end

"""
    _laplace_hessian_columns(laplace, values, columns; step, retries, warm,
        progress, localised)

`ctsem_laplace_hessian`'s columns `columns` alone, unsymmetrised, as an `n` by
`length(columns)` matrix: each differenced from the gradient at two points as
there, with the same inner-budget escalation and `NaN` for a column whose
points did not both solve their modes. With `localised`, a column whose
parameter some units do not read is differenced on the objective over the
units that do (`_ctsem_subset_objective`), one such objective for every set of
units that some column needs, its modes started where the full objective's
were.
"""
function _laplace_hessian_columns(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, columns; step::Real=1e-4, retries::Integer=2,
    warm::Bool=true, progress=nothing, localised::Bool=true)
    # Over the extended vector (`_laplace_extend`), and returned with the
    # caller's rows.
    nin = length(values)
    x = _laplace_extend(laplace, collect(Float64, values))
    n = length(x)
    cols = collect(Int, columns)
    k = length(cols)
    H = zeros(Float64, n, k)
    k == 0 && return H[1:nin, :]
    budget = laplace.inner_maxiter
    NU = length(laplace.units.members)
    # The modes at `x`, from the origin, which every differencing point then
    # starts from; and, for `localised`, each unit's score there, which says
    # which parameters it reads.
    scores = localised ? zeros(Float64, NU, n) : nothing
    if warm || localised
        previous = ctsem_set_warm_start!(false)
        try
            # Escalated as the points are (`_laplace_hessian_point`): the
            # points start from these modes, and modes a short budget stopped
            # short of sent them along other paths than the full budget's, to
            # other ends within the inner tolerance -- columns 3e-8 apart
            # relative where the escalation should have changed nothing.
            out = ctsem_laplace_evaluate(laplace, x; gradient=localised,
                scores=scores)
            tries = 0
            while !out.converged && tries < retries &&
                    any(>=(laplace.inner_maxiter), laplace.inner_iterations)
                tries += 1
                laplace.inner_maxiter = budget * 4^tries
                out = ctsem_laplace_evaluate(laplace, x; gradient=localised,
                    scores=scores)
            end
        finally
            laplace.inner_maxiter = budget
            ctsem_set_warm_start!(previous)
        end
    end
    base = warm ? deepcopy(laplace.modes) : nothing
    # Coefficient columns taken from their targets' score differences
    # (`_laplace_ti_derivable`): `derived[c]` is the TI effect entry of
    # column `c`, and the targets are differenced with their units' scores
    # kept, whether or not their own columns were asked for.
    derived = localised ? _laplace_ti_derivable(laplace, cols) : Dict{Int,Int}()
    sp = laplace.objective.params
    targets = sort!(unique([sp.ti_parameter_indices[e] for e in Base.values(derived)]))
    extra = setdiff(targets, cols)
    if length(cols) - length(derived) + length(extra) >= length(cols)
        empty!(derived)
        targets = Int[]
        extra = Int[]
    end
    work = vcat(cols, extra)
    targetscores = Dict{Int,Matrix{Float64}}()
    # The units each column is differenced over, grouped so that one subset
    # objective serves every column that reads the same units. All units, or
    # scores that are not finite, or a parameter no unit reads (only the prior
    # does), mean the full objective.
    groups = Dict{Vector{Int},Vector{Int}}()
    full = Int[]
    for c in eachindex(work)
        c <= k && haskey(derived, c) && continue
        reads = scores === nothing || !all(isfinite, scores) ? Int[] :
            findall(!iszero, view(scores, :, work[c]))
        if isempty(reads) || length(reads) == NU
            push!(full, c)
        else
            push!(get!(groups, reads, Int[]), c)
        end
    end
    plan = Tuple{Vector{Int},Vector{Int}}[(Int[], full)]
    for (reads, cs) in groups
        push!(plan, (reads, cs))
    end
    previous = ctsem_set_warm_start!(warm)
    unconverged = Int[]
    escalated = Int[]
    worst = 0.0
    units = 0
    done = 0
    try
        for (reads, cs) in plan
            isempty(cs) && continue
            obj = isempty(reads) ? laplace : _ctsem_subset_objective(laplace, reads)
            from = base === nothing ? nothing :
                isempty(reads) ? base : [base[U] for U in reads]
            for c in cs
                j = work[c]
                h = step * max(1.0, abs(x[j]))
                plus = copy(x); plus[j] += h
                minus = copy(x); minus[j] -= h
                target = j in targets
                nreads = isempty(reads) ? NU : length(reads)
                splus = target ? zeros(Float64, nreads, n) : nothing
                sminus = target ? zeros(Float64, nreads, n) : nothing
                pp = _laplace_hessian_point(obj, plus, budget, retries; from=from,
                    scores=splus)
                pm = _laplace_hessian_point(obj, minus, budget, retries; from=from,
                    scores=sminus)
                (pp.tries + pm.tries) > 0 && push!(escalated, j)
                # Two gradients just finished, converged or not: a column this
                # function gives up on (NaN, below) still cost the same two
                # inner solves, so it counts here too or the reported total
                # falls behind the real one.
                done += 1
                progress !== nothing && progress(2done, 2 * (length(work) - length(derived)))
                if target && pp.out.converged && pm.out.converged
                    # Every unit's score difference along `j`, by full unit.
                    D = zeros(Float64, NU, n)
                    rows = isempty(reads) ? (1:NU) : reads
                    D[rows, :] .= (splus .- sminus) ./ (2h)
                    targetscores[j] = D
                end
                c > k && continue
                if !(pp.out.converged && pm.out.converged)
                    # The same predicate `fg!` applies to a trial point, applied
                    # to the two points this column is differenced from. A
                    # gradient at a not-quite-mode is a gradient of a different
                    # function, and one that converged to 1e-9 instead of 1e-10
                    # is finite, so nothing downstream would notice. NaN makes
                    # the column visible to the non-finite fallback callers
                    # already have.
                    push!(unconverged, j)
                    worst = max(worst, pp.worst, pm.worst)
                    units = max(units, pp.units, pm.units)
                    H[:, c] .= NaN
                    continue
                end
                H[:, c] = (pp.out.gradient .- pm.out.gradient) ./ (2h)
            end
        end
        isempty(derived) || _laplace_ti_derived_columns!(H, laplace, x, cols,
            derived, targetscores)
    finally
        # The budget is a field of a mutable objective the caller keeps, so an
        # escalation that threw would otherwise be inherited by every later
        # evaluation of it; the warm-start switch is global, and the same.
        laplace.inner_maxiter = budget
        ctsem_set_warm_start!(previous)
        # Where the caller's objective's modes were left: at `x`, as the full
        # objective's last column would have left them near it.
        _laplace_restore_modes!(laplace, base)
    end
    # A count, not a list of indices: on a model with a thousand parameters the
    # list is the whole message and says no more than the count does.
    isempty(escalated) || @info string("Laplace inner budget raised above ",
        budget, " for ", length(escalated), " of ", k, " Hessian columns.")
    isempty(unconverged) || @warn string("Laplace inner solve did not converge ",
        "at the Hessian step for parameter ", join(unconverged, ", "),
        "; those columns are NaN. ", units, " unit(s) short, worst inner ",
        "gradient ", worst, " against a budget of ", budget * 4^retries,
        " iterations.")
    return n == nin ? H : H[1:nin, :]
end

"""
    _laplace_ti_derivable(laplace, cols)

The columns among `cols` that are TI-predictor coefficients whose Hessian
column can be taken from their target's (`_laplace_hessian_columns`), as
column position => TI effect entry: every unit one subject, so a coefficient
moves a unit through its own predictor value; no sampled predictor value,
which would make that value a parameter too; and a coefficient with one
target, not one named effect shared between parameters.
"""
function _laplace_ti_derivable(laplace::CTSEMLaplaceObjective, cols::AbstractVector{<:Integer})
    out = Dict{Int,Int}()
    sp = laplace.objective.params
    isempty(sp.ti_coefficient_indices) && return out
    all(m -> length(m) == 1, laplace.units.members) || return out
    isempty(laplace.objective.ti_missing_parameter) || return out
    all(s -> s.tipreds isa AbstractVector{<:Real},
        laplace.objective.subject_objectives) || return out
    coefficients = sp.ti_coefficient_indices
    for (c, j) in enumerate(cols)
        entries = findall(==(j), coefficients)
        length(entries) == 1 || continue
        e = entries[1]
        # The target's own column must be differenced, not derived.
        sp.ti_parameter_indices[e] in coefficients && continue
        out[c] = e
    end
    return out
end

"""
    _laplace_ti_derived_columns!(H, laplace, x, cols, derived, targetscores)

Fill each derived column of `H` from its target's units' score differences:
`sum_U x_Uk d s_U / d theta_j`, plus the prior's and penalty's curvature at the
coefficient itself. The units' scores carry no prior, which is diagonal in a
coefficient. A target whose differencing failed leaves its columns `NaN`.
"""
function _laplace_ti_derived_columns!(H::AbstractMatrix, laplace::CTSEMLaplaceObjective,
    x::AbstractVector, cols::AbstractVector{<:Integer}, derived::Dict{Int,Int},
    targetscores::Dict{Int,Matrix{Float64}})
    sp = laplace.objective.params
    subjects = laplace.objective.subject_objectives
    members = laplace.units.members
    precision = _ctsem_prior_precision(laplace.objective, length(x))
    for (c, e) in derived
        j = sp.ti_parameter_indices[e]
        k = sp.ti_predictor_indices[e]
        D = get(targetscores, j, nothing)
        if D === nothing
            H[:, c] .= NaN
            continue
        end
        column = zeros(Float64, size(H, 1))
        for U in eachindex(members)
            value = subjects[members[U][1]].tipreds[k]
            iszero(value) && continue
            @inbounds for m in eachindex(column)
                column[m] += value * D[U, m]
            end
        end
        column[cols[c]] -= precision[cols[c]]
        H[:, c] .= column
    end
    return H
end

"""
    _laplace_restore_modes!(laplace, modes)

Put `laplace`'s inner modes back to `modes` (one vector per unit, copied), the
start every differencing point is solved from; nothing when `modes` is
`nothing`.
"""
function _laplace_restore_modes!(laplace::CTSEMLaplaceObjective, modes)
    modes === nothing && return nothing
    for U in eachindex(modes)
        laplace.modes[U] = copy(modes[U])
    end
    return nothing
end

"""
    _laplace_hessian_point(laplace, x, budget, retries; from=nothing)

One of the two points a Hessian column is differenced from, evaluated with the
inner budget raised rather than the column abandoned. Each attempt starts from
`from`, the warm start (`nothing`: the origin): a retry from where the short
attempt stopped walks another path than the larger budget would have from the
start, so ends elsewhere within the inner tolerance, and the column is not the
one the larger budget gives.

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
    x::AbstractVector, budget::Integer, retries::Integer; from=nothing,
    scores::Union{Nothing,Matrix{Float64}}=nothing)
    capped() = any(>=(laplace.inner_maxiter), laplace.inner_iterations)
    reached() = isempty(laplace.inner_gradient) ? 0.0 :
        maximum(laplace.inner_gradient)
    short() = count(!, laplace.inner_converged)
    laplace.inner_maxiter = budget
    _laplace_restore_modes!(laplace, from)
    out = ctsem_laplace_evaluate(laplace, x; gradient=true, scores=scores)
    tries = 0
    while !out.converged && capped() && tries < retries
        tries += 1
        laplace.inner_maxiter = budget * 4^tries
        _laplace_restore_modes!(laplace, from)
        out = ctsem_laplace_evaluate(laplace, x; gradient=true, scores=scores)
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

"""The outer Hessian, for the generic curvature entry point. `progress`, when
given, reports the finite-difference loop's gradients as they finish -- see
`ctsem_laplace_hessian`."""
ctsem_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    chunk::Integer=0, progress=nothing) =
    ctsem_laplace_hessian(laplace, values; progress=progress)

"""The outer Hessian's `columns`, for `ctsem_hessian_columns`."""
ctsem_hessian_columns(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    columns; progress=nothing) =
    _laplace_hessian_columns(laplace, values, columns; progress=progress)

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

The rows are the units' own gradients from the evaluation's seeded reverse
sweep (`ctsem_laplace_evaluate`'s `scores`), so they cost one gradient and run
in parallel over units like it. They were a forward-mode Jacobian of the
per-unit terms (`_laplace_subject_gradients_forward`, still the fallback): one
dual pass per chunk of parameters over every unit, on one thread. The batched
optimiser takes every gradient of a batch stage from here, and on the SNSF
pilot (310 parameters, 493 single-member units) its first iterations then took
5-20 minutes each at one core, against seconds once the batch was full.

`prior = false` leaves each unit's share of the prior out, as for the plain
objective.
"""
function ctsem_subject_gradients(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; prior::Bool=true)
    # Over the extended vector (`_laplace_extend`), and returned with the
    # caller's columns.
    n = length(values)
    theta = _laplace_extend(laplace, collect(Float64, values))
    nunits = length(laplace.units.members)
    S = zeros(Float64, nunits, length(theta))
    r = ctsem_laplace_evaluate(laplace, theta; gradient=true, scores=S)
    isfinite(r.value) || return (value=r.value, scores=fill(NaN, nunits, n))
    if !all(isfinite, S)
        f = _laplace_subject_gradients_forward(laplace, theta; prior=prior)
        return (value=f.value, scores=size(f.scores, 2) == n ? f.scores : f.scores[:, 1:n])
    end
    # Each unit carries its share of the prior, as the forward route's rows do.
    prior && _ctsem_prior_scores!(S, laplace.objective, theta,
        _laplace_unit_of_subject(laplace))
    return (value=r.value, scores=size(S, 2) == n ? S : S[:, 1:n])
end

"""
    _laplace_subject_gradients_forward(laplace, values)

`ctsem_subject_gradients` by a forward-mode Jacobian of the per-unit terms: the
reference its seeded rows are tested against, and its fallback when the seeded
sweep fails. One sweep, not one per unit: the per-unit terms are assembled into
a vector and differentiated together, so the same forward directions serve
every row.
"""
function _laplace_subject_gradients_forward(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; prior::Bool=true)
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
    scores = _ctsem_jacobian(terms_of, theta)

    # Each unit carries its share of the prior, so the rows still sum to the
    # full posterior gradient -- the same convention `ctsem_subject_gradients`
    # uses without random effects.
    prior && _ctsem_prior_scores!(scores, laplace.objective, theta,
        _laplace_unit_of_subject(laplace))
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
    ctsem_generate(laplace, values, base; effects=nothing, subject_values=nothing, seed=1)

One posterior-predictive dataset from a Laplace fit.

Each subject is generated at its own realized parameters -- the population
vector shifted by that subject's estimated random effects and TI-predictor
effects -- rather than at the shared population vector, using exactly the
per-subject values `ctsem_kalman(laplace, ...)` already computes for
prediction (`_laplace_filter_values`). A Laplace random effect is a
conditional mode estimated from the subject's whole record, not a fresh draw
from the population distribution, so the individual difference in the
generated data is the fitted one: this is the same smoothed-equivalent,
conditional-on-the-subject's-own-data quantity `ctsem_kalman(laplace, ...)`
already reports for residuals and predictions on this route, not a new
statistical convention introduced for generation.

`effects` is one draw of the random effects as natural deviations, for a fit
that sampled them, as `ctsem_generate_states` takes it. `subject_values`
supplies the per-subject vectors instead of solving for them here, before
TI-predictor effects, for the same reason `ctsem_kalman` takes it: a caller
filtering different rows than the fit did needs the fitted modes, not modes
re-solved against rows that may not condition on anything.
"""
function ctsem_generate(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    base::AbstractMatrix; effects::Union{Nothing,AbstractVector}=nothing,
    subject_values::Union{Nothing,AbstractMatrix}=nothing, seed::Integer=1)
    persubject = _laplace_filter_values(laplace, values; effects=effects,
        subject_values=subject_values)
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

# The effects a sampler draws, `u`, are coordinates: each block's deviation of
# its level's parameters is `d = S u` for a square root `S` of the level's
# population covariance, and which root is a choice (the Cholesky here). What a
# fit reports and what generation takes are the deviations, which are the same
# whatever the root and are what a parameter's random effect is.

"""
    _laplace_deviation_offsets(laplace)

Each block's place in the flat vector of natural deviations: a zero-based start
per unit per block (in `units.blocks[U]`'s order), blocks laid out in offset
order and each taking its level's `k` entries -- one per parameter the level
moves. On a full-rank level that is the effects' own layout; a reduced-rank
block is wider, its `rank` effects moving `k` parameters. Returns the starts and
the total length.
"""
function _laplace_deviation_offsets(laplace::CTSEMLaplaceObjective)
    units = laplace.units
    starts = [zeros(Int, length(units.blocks[U])) for U in eachindex(units.members)]
    total = 0
    for U in eachindex(units.members)
        blocks = units.blocks[U]
        for b in sortperm([blk.offset for blk in blocks])
            starts[U][b] = total
            total += nrandomeffects(laplace.spec.levels[blocks[b].level])
        end
    end
    return starts, total
end

"""
    _laplace_deviation_blocks(laplace)

The natural deviations' blocks in their layout's order: each block's
zero-based start and width (its level's `k`). A block is one group's effects at
one level, so sums of products within blocks (`_effect_cross_add!`) give every
group's posterior covariance from a sampler's draws.
"""
function _laplace_deviation_blocks(laplace::CTSEMLaplaceObjective)
    starts, _ = _laplace_deviation_offsets(laplace)
    out = Tuple{Int,Int}[]
    for U in eachindex(laplace.units.members), (b, blk) in enumerate(laplace.units.blocks[U])
        push!(out, (starts[U][b], nrandomeffects(laplace.spec.levels[blk.level])))
    end
    sort!(out; by=first)
    return out
end

"""Length of the packed within-block products for `blocks`."""
_effect_cross_length(blocks) = sum(k * (k + 1) ÷ 2 for (_, k) in blocks; init=0)

"""
    _effect_cross_add!(cross, dev, blocks)

Add one draw's within-block products of the deviations `dev` to `cross`,
packed block by block, each column by column from the diagonal down.
"""
function _effect_cross_add!(cross::AbstractVector{Float64}, dev::AbstractVector,
    blocks::Vector{Tuple{Int,Int}})
    at = 0
    @inbounds for (s, k) in blocks
        for jb in 1:k, ja in jb:k
            at += 1
            cross[at] += dev[s + ja] * dev[s + jb]
        end
    end
    return cross
end

"""
    _effect_cross_index(blocks)

The two deviation positions (one-based) of each packed entry, in packing
order: what pooling covariances from several processes needs for the part
between them.
"""
function _effect_cross_index(blocks::Vector{Tuple{Int,Int}})
    rows = Int[]
    cols = Int[]
    for (s, k) in blocks, jb in 1:k, ja in jb:k
        push!(rows, s + ja)
        push!(cols, s + jb)
    end
    return rows, cols
end

"""
    _effect_cross_cov(cross, mean, n, blocks)

Each block's sample covariance, packed as `cross` is, from the sums of products
over `n` draws and the deviations' means.
"""
function _effect_cross_cov(cross::AbstractVector{Float64}, mean::AbstractVector{Float64},
    n::Integer, blocks::Vector{Tuple{Int,Int}})
    out = zeros(Float64, length(cross))
    at = 0
    for (s, k) in blocks, jb in 1:k, ja in jb:k
        at += 1
        out[at] = (cross[at] / n - mean[s + ja] * mean[s + jb]) * (n / max(n - 1, 1))
    end
    return out
end

"""
    _laplace_deviations!(dest, laplace, Ls, effects, starts)

One draw's natural deviations into `dest` (`_laplace_deviation_offsets`'s
layout) from its effects `effects` (the sampler's layout, unit by unit) and the
population roots `Ls` at that draw's parameters.
"""
function _laplace_deviations!(dest::AbstractVector, laplace::CTSEMLaplaceObjective,
    Ls::Vector{<:AbstractMatrix}, effects::AbstractVector, starts)
    units = laplace.units
    ubase = 0
    @inbounds for U in eachindex(units.members)
        for (b, blk) in enumerate(units.blocks[U])
            L = Ls[blk.level]
            s = starts[U][b]
            for p in axes(L, 1)
                acc = 0.0
                for q in axes(L, 2)
                    acc += L[p, q] * effects[ubase + blk.offset + q]
                end
                dest[s + p] = acc
            end
        end
        ubase += units.dims[U]
    end
    return dest
end

"""
    ctsem_laplace_deviation_layout(laplace)

`ctsem_laplace_effect_layout`'s account for the natural deviations a sampled
fit reports: one entry per position, giving the unit, the level, the block's
first member, how many members it covers, and which of the level's parameters
(`within`, 1 to `k`) the position is.
"""
function ctsem_laplace_deviation_layout(laplace::CTSEMLaplaceObjective)
    units = laplace.units
    starts, total = _laplace_deviation_offsets(laplace)
    unit = zeros(Int, total); level = zeros(Int, total)
    first_member = zeros(Int, total); nmembers = zeros(Int, total)
    within = zeros(Int, total)
    for U in eachindex(units.members), (b, blk) in enumerate(units.blocks[U])
        for p in 1:nrandomeffects(laplace.spec.levels[blk.level])
            i = starts[U][b] + p
            unit[i] = U
            level[i] = blk.level
            first_member[i] = isempty(blk.members) ? 0 : units.members[U][blk.members[1]]
            nmembers[i] = length(blk.members)
            within[i] = p
        end
    end
    return (position=collect(1:total), unit=unit, level=level,
        first_member=first_member, nmembers=nmembers, within=within)
end

"""
    ctsem_laplace_coordinate_dimension(laplace)

How many standard-normal coordinates `u` one draw of every random effect takes:
the samplers' layout, unit by unit, which on a reduced-rank level is narrower
than the deviations it builds.
"""
ctsem_laplace_coordinate_dimension(laplace::CTSEMLaplaceObjective) =
    sum(laplace.units.dims; init=0)

"""
    ctsem_laplace_population_deviations(laplace, values, u)

Natural deviations `S u` -- each block's deviation of its level's raw
parameters, in `ctsem_laplace_deviation_layout`'s layout -- for coordinates `u`
at the population covariance `values` implies. With `u` independent standard
normals this is a draw of new groups and subjects from the population, which is
what generation from a model, rather than from a fit's subjects, needs. The
deviations go to `ctsem_generate` or `ctsem_generate_states` as their `effects`.

The same construction the samplers report their draws through
(`_laplace_popchols`, `_laplace_deviations!`), so a deviation means one thing
however it was made.
"""
function ctsem_laplace_population_deviations(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, u::AbstractVector)
    n = ctsem_laplace_coordinate_dimension(laplace)
    length(u) == n || throw(DimensionMismatch(string("u has ", length(u),
        " coordinates where the random effects take ", n)))
    theta = _laplace_extend(laplace, collect(Float64, values))
    _laplace_check_indices(laplace, length(theta))
    Ls = _laplace_popchols(theta, laplace.spec)
    starts, ndev = _laplace_deviation_offsets(laplace)
    return _laplace_deviations!(zeros(Float64, ndev), laplace, Ls,
        collect(Float64, u), starts)
end

"""
    _laplace_shifted_values_from_deviations(laplace, theta, deviations; from_level=1)

Every subject's raw parameters before TI-predictor effects, from the population
vector and one draw's natural deviations: each of its blocks' deviations added
to that level's parameters, for the levels from `from_level` outward. No
population root is involved, so this is the same whichever root the draws were
made under.
"""
function _laplace_shifted_values_from_deviations(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, deviations::AbstractVector; from_level::Integer=1)
    units = laplace.units
    spec = laplace.spec
    starts, total = _laplace_deviation_offsets(laplace)
    length(deviations) == total || throw(DimensionMismatch(string(
        "deviations must have ", total, " entries for this design, got ",
        length(deviations))))
    out = fill(NaN, length(laplace.objective.subject_objectives), length(theta))
    for U in eachindex(units.members)
        blocks = units.blocks[U]
        for (m, i) in enumerate(units.members[U])
            shifted = copy(theta)
            for l in eachindex(spec.levels)
                l < from_level && continue
                level = spec.levels[l]
                k = nrandomeffects(level)
                (k == 0 || nlatent(level) == 0) && continue
                offset = units.offsets[U][m][l]
                b = findfirst(blk -> blk.level == l && blk.offset == offset, blocks)
                b === nothing && continue
                s = starts[U][b]
                for p in 1:k
                    shifted[level.re_index[p]] += deviations[s + p]
                end
            end
            out[i, :] = shifted
        end
    end
    return out
end

"""The parameters a saturation check reads, for the wrapped objective."""
_ctsem_params(o::CTSEMLaplaceObjective) = _ctsem_params(o.objective)
