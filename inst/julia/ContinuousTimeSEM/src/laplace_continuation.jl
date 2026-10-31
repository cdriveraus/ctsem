"""
The quadrature continuation: a short climb, from the Laplace optimum, of the
objective with the flagged units' integrals taken by quadrature.

# Why this exists

`ctsem_laplace_quadrature` scores a Laplace fit's units by adaptive
Gauss-Hermite quadrature, and at the Laplace optimum the two often disagree by
more than a fit can afford: on AnomAuth's spurious maximum by 28.8 nats, 444
of 800 units each over-credited by about 0.035. The default correction steps
towards the quadrature optimum by finite differences of that objective, `2
npar` quadratures per step, and stops a tenth of a standard error short on its
own fixture. What was missing was a gradient.

With the nodes held fixed, the quadrature objective has one that costs no
third derivatives and no mode adjoint. For a unit with fixed nodes `u_j` and
fixed log weights `w_j`,

    Q_U(theta) = log sum_j exp(w_j + g_U(u_j; theta))

and its derivative is `sum_j pi_j dg_U(u_j; theta)/dtheta` with `pi` the
normalised node weights -- the Fisher identity, a node-weighted average of the
per-node score. `dg_U(u; theta)/dtheta` at a fixed `u` is one reverse sweep
per member at the shifted parameters plus the chain rule through `L(theta)`,
which is what `_laplace_floored_unit_gradient!` already takes at the mode.
This is bigIRT's `laplaceRefine` design (`../bigIRT/package/src/agq_backend.cpp`,
`R/laplaceBackend.R`), adapted to units that are whole filters, that nest, and
that can be wide.

# The objective, as built

A *hybrid*: the units the screen flags carry fixed-node rules and every other
unit keeps its Laplace term and its seeded gradient (through a unit subset of
the Laplace objective, which also carries the prior, once). Flagged is decided
per unit, largest `|Q_U - T_U|` first, until what is left carries at most the
screen's own tolerance: the units left to Laplace are, together, what the
screen would have called exact. A unit once flagged stays flagged when the
nodes move, so no placement loses a unit an earlier one scored.

Nodes are held fixed in the *standardised* effects `u` -- the non-centred
coordinates the engine already works in, where a member's parameters are
`theta + L(theta) u` and the prior is `N(0, I)` whatever `theta` is. That
choice keeps the prior's density out of the gradient entirely, handles a
reduced-rank level (whose `L` has zero columns and no density in the raw effect
space) without a special case, and reuses the engine's fixed-`u` score as it
stands. Its cost is that a node fixed in `u` moves in the raw effect space as a
population scale changes, so the refined fit climbs an objective whose product
rules are placed at every point it evaluates (`laplace_placed.jl`), its frozen
draws placed again where their effective size has collapsed.

# The rule

Per block of the unit's tree, the rule `_quadrature_block` uses: a leaf
re-solves its conditional mode at each node of the blocks above it (held fixed
afterwards: fixed conditional maps), an outer block keeps the joint mode with
the eliminated precision's scale. At either, the block's point set
(`point_sets.jl`) is placed at that centre and the prior-clipped scale: product
Gauss-Hermite, `nodes^k` points, or lattice-t, `128` points and up whatever the
width. Which a unit takes, and how many points, is chosen by counted cost for
the accuracy at the screen and kept (`o.unit_rules`, `_continuation_refine!`):
there is no width at which a unit stops being scored. The placed objective
(`laplace_placed.jl`) places the same sets at every evaluation, so at the
centre it and the continuation agree.

# Where the constants come from

`nodes = 5` and `tolerance = 0.01` are the screen's (`.ctLaplaceCorrectDefaults`
on the R side). `draws_se = 0.05` is the error each flagged unit's rule is
refined to (`_continuation_refine!`). The lattice's constants are in
`point_sets.jl`. None of these was tuned on a fit.

# Units no rule converges on: frozen importance draws

A unit whose product and lattice rules both fail to converge
(`_continuation_refine!`), or every unit under `draws_all` (the random-effects
check, `ctsem_laplace_check`), gets draws in place of nodes
(`_continuation_importance_rule`):
held fixed in the standardised effects, as nodes are, from a defensive
mixture -- a Gaussian at the block's mode, at any other point the caller knows
for it (`points`: where the inner solve stopped at nearby parameters), at
draws that carried weight once moved a Newton step, and the prior -- and
weighted by that mixture. A nested unit gets them per block of its tree, as a
product rule is placed: each block's draws conditional on a draw of the blocks
above, an outer block's rule holding its children's rules at each of its
draws, so each integral is one block's effects. Its cost is the number of
draws whatever the width -- per member, the product of the draws along its
path -- and the same Fisher identity gives its gradient and the same Louis
identity its Hessian. Why: a unit whose posterior has a broad
shoulder beside its mode has a Laplace term that moves with where the inner
solve stops, by up to the shoulder's mass, and an optimiser climbing onto the
point where it stops on the shoulder stalls there with a gradient that is that
spike's (`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`: 0.33 nats
on an 8-effect person, the fit stalled for hours). A frozen rule does not read
the stopping point, so its term is smooth. Such a unit's posterior is not one
basin a placed rule can follow, which is why it is the one left on draws.

Its value carries Monte Carlo error, read from the spread of its replicate
lattices (`_CONTINUATION_REPLICATES`), which no number of draws a fit can
afford takes below the screen's tolerance. So the screen flags what it can
resolve (`_continuation_flags`): a unit whose own gap exceeds twice its error
by more than the tolerance, and, since small gaps of one sign sum to a
resolved total -- 60-80 per cent of the gap on acnonlin and gA1 sat in units
under twice their error (review/LAPLACE-frozen-draws-validation-2026-10-22.md)
-- units largest first until the unflagged units' signed sum is within the
tolerance plus twice its combined error. A unit the caller names (`force`)
is flagged whatever its gap: the optimiser
names a unit whose term jumped (`_ctsem_jump_units`), and at a stall the gap
at the stalled point can be small, since that is where the solve stopped on
the shoulder.

The constants. `draws_floor = 0.25` floors each component's precision
eigenvalues, so a soft or indefinite direction is sampled at up to twice the
prior's sd; `draws_prior = 0.2` is the prior's share of every batch. Both are
the 20000-draw reference's, which reached an ESS of 8000 on that person. The
first batch of `draws = 256`, doubling, to a standard error of `draws_se =
0.05`: on that person a rule placed from the stiff mode alone reached it at
1024 draws with seven components, its error against the reference at most
0.019 nats over four seeds and constant along the fold to 1e-3, and its
directional derivative within 1e-4 of its own differences. `draws_max = 16384`
is cost: a gradient is one reverse sweep per member per draw. Every component
but the prior of a unit's top block is a Student t on `draws_df = 4` degrees
of freedom (`CTSEMImportanceComponent` says why); `Inf` makes them Gaussian,
as the blocks beneath another are.
"""

using LinearAlgebra

export ctsem_laplace_continuation, ctsem_laplace_continuation_info,
    ctsem_laplace_continuation_evaluate, ctsem_laplace_continuation_se

"""
One block's fixed rule, for one configuration of the blocks above it.

`points[:, j]` is the block's own coordinates at node `j`, in the unit's
standardised `u`; `logweights[j]` carries everything about the node that does
not depend on the parameters (the rule's weight, the Gaussian factor of any
stiff complement, and the block's own `-z'z/2`), and `constant` what the block
adds to its log integral whatever the node. A leaf's value at a node is its
members' log likelihood there; an outer block's is the sum of its children's
rules, and `children[j]` holds those rules as they were placed with this block
at node `j`. `replicate` is an importance or lattice rule's replicate lattice
for each point (`_CONTINUATION_REPLICATES`, `_LATTICE_SHIFTS`) and
`replicate_logweights` each point's log weight within its replicate's own
estimate, both fixed with the points, which its standard error is read from
away from where it was placed (`ctsem_laplace_continuation_se`); empty for a
product rule.
"""
struct CTSEMFixedRule
    block::Int
    points::Matrix{Float64}
    logweights::Vector{Float64}
    constant::Float64
    children::Vector{Vector{CTSEMFixedRule}}
    replicate::Vector{Int}
    replicate_logweights::Vector{Float64}
end

"""A unit's fixed rule, one `CTSEMFixedRule` per root of its block tree, and
how many member evaluations one value of it costs."""
struct CTSEMFixedUnitRule
    unit::Int
    roots::Vector{CTSEMFixedRule}
    evaluations::Int
end

_continuation_evaluations(rule::CTSEMFixedRule, blocks) =
    isempty(rule.children) ?
        size(rule.points, 2) * length(blocks[rule.block].members) :
        sum(sum(_continuation_evaluations(c, blocks) for c in kids; init=0)
            for kids in rule.children; init=0)

"""The points a rule holds, at every block of its tree: an importance rule's
draws in all."""
_continuation_draws(rule::CTSEMFixedRule) = size(rule.points, 2) +
    sum(sum(_continuation_draws(c) for c in kids; init=0)
        for kids in rule.children; init=0)

"""
    CTSEMLaplaceContinuation

The hybrid objective with the flagged units' rules frozen where they were
placed, with the state it was placed from. See the module docstring. No run
optimises it: the refined fit climbs the placed objective built on it
(`laplace_placed.jl`), and the shift check evaluates it
(`ctsem_laplace_check`).
"""
mutable struct CTSEMLaplaceContinuation{L}
    laplace::L
    rest::L
    rest_units::Vector{Int}
    flagged::Vector{Int}
    rules::Vector{CTSEMFixedUnitRule}
    centre::Vector{Float64}
    nodes::Int
    tolerance::Float64
    # Every unit on importance draws, whatever its rules would do: the
    # random-effects check's continuation (`ctsem_laplace_check`).
    draws_all::Bool
    # The importance rule's options (see the module docstring), the points the
    # caller knows for each unit, and per unit at the centre: the rule's
    # standard error (0 for a product rule) and whether it is an importance rule.
    draws::Int
    draws_se::Float64
    draws_max::Int
    draws_floor::Float64
    draws_prior::Float64
    draws_df::Float64
    draws_seed::Int
    points::Dict{Int,Vector{Vector{Float64}}}
    centre_se::Vector{Float64}
    # Each flagged unit's rule error at the centre, in nats
    # (`_continuation_node_errors`); `NaN` for an unflagged unit.
    centre_error::Vector{Float64}
    importance::Vector{Bool}
    # Whether a placement scores every unit, or only the flagged ones (a fit
    # rescued from a jump, which needs the units the optimiser named and no
    # screen).
    screen::Bool
    # The units `force` named.
    forced::Vector{Int}
    # What the rule scored every unit at the centre, and what Laplace did.
    centre_quadrature::Vector{Float64}
    centre_laplace::Vector{Float64}
    rule_failed::Vector{Bool}
    centre_value::Float64
    # Accounting over the object's life.
    value_calls::Int
    gradient_calls::Int
    member_values::Int
    member_sweeps::Int
    recentres::Int
    # Wall seconds in each kind of call, so a correction's cost can be broken
    # down without a profiler: values, gradients, placements (the screen is
    # the first).
    seconds_values::Float64
    seconds_gradients::Float64
    seconds_placements::Float64
    # Whether the Laplace part's inner solves converged at the centre, so a
    # value read off the placement (`centre_value`) can be refused as the
    # evaluation refuses one.
    centre_converged::Bool
    # The units on importance draws because no rule reached the draws' target
    # at some placement (`_continuation_refine!`), which keep their draws at
    # every later one. Let go, a unit whose product rule's error sits
    # at the target switched between the two from one placement to the next,
    # and the value with it, by more than either's error: gN3's unit 3, its
    # rule's error 0.050 against a target of 0.05, cost 0.11-0.135 nats at
    # four of fourteen rounds, each then rejected, and the rounds stalled.
    drawn::Vector{Int}
    # Each unit's rule (`CTSEMRuleChoice`): the one it was screened on
    # (`_rule_default`), or the one the refine chose for it
    # (`_continuation_refine!`), kept at every later placement. Let go, a rule
    # whose cheap check sat near the target took the finer count at one
    # placement and the default at the next, and the value moved by the
    # difference between the two rules: on gN3, with `drawn` kept, the two
    # rounds still rejected were that switch.
    unit_rules::Vector{CTSEMRuleChoice}
end

"""The rule unit `U` of `o` is on (`o.unit_rules`)."""
_continuation_choice(o::CTSEMLaplaceContinuation, U::Integer) = o.unit_rules[U]


"""
    _continuation_block_functions(laplace, U, theta, Ls, b, u, aws)

Block `b`'s own log integrand and its precision, as functions of the block's
coordinates with the rest of `u` held where it is now: `value_gradient(z)` is
its members' log likelihood less `z'z/2` and that function's gradient, and
`precision(z)` is `I - d2 ll/dz dz`, dense (`_laplace_block_precision`) -- what
`_quadrature_leaf_rule!` uses for the same block, with the deviations each
member's filter integrates integrated exactly, as there.
"""
function _continuation_block_functions(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, aws)
    block = laplace.units.blocks[U][b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    members = block.members
    base = copy(u)
    value_gradient = function (z)
        work = copy(base)
        @inbounds for (t, c) in enumerate(columns); work[c] = z[t]; end
        r = _laplace_unit_loglik_gradient(laplace, U, theta, Ls, work, aws, members,
            nothing, _HYBRID_FULL)
        isfinite(r.value) || return (value=-Inf, gradient=fill(NaN, k))
        return (value=r.value - dot(z, z) / 2,
            gradient=[r.gradient[c] for c in columns] .- z)
    end
    precision = z -> _laplace_block_precision(laplace, U, theta, Ls, base, b, z;
        mode=_HYBRID_FULL)
    return (value_gradient=value_gradient, precision=precision)
end

"""
    _continuation_lattice_normals(d, n, rng)

`n` standard normal points in `d` dimensions, as columns, from a randomly
shifted rank-1 lattice (`frac(j sqrt(p_a) + s_a)`, `p_a` the `a`-th prime, one
shift `s` from `rng`) mapped through the normal quantile: the points
`_laplace_mixture_unit!` draws each component's particles from. Quasi-random,
because the integrands are smooth and their error falls faster than with
pseudo-random draws at the same number.
"""
function _continuation_lattice_normals(d::Integer, n::Integer, rng)
    alpha = [mod(sqrt(Float64(p)), 1.0) for p in _mixture_primes(Int(d))]
    shift = rand(rng, Int(d))
    Z = Matrix{Float64}(undef, Int(d), Int(n))
    @inbounds for j in 1:Int(n), a in 1:Int(d)
        Z[a, j] = sqrt(2.0) * erfinv(clamp(2 * mod(j * alpha[a] + shift[a], 1.0) - 1,
            -1 + 1e-15, 1 - 1e-15))
    end
    return Z
end

"""
A component of an importance rule's proposal in a block's own standardised
coordinates: centre `mu`, scale matrix the inverse of `vecs * Diagonal(lam) *
vecs'` -- the block's precision at the centre, its eigenvalues floored -- and
`df` degrees of freedom: a Student t, or a Gaussian at `Inf`, which the prior
component is. A t, because a Gaussian proposal at a posterior's mode is
lighter-tailed than the posterior, and every replicate lattice misses the
same tail: on the bench's gA1 the rule stopped on its replicates' spread read
0.6 to 1.0 nats low over 40 units with the spread saying 0.1, 0.2 to 0.35 on
the t at the same draws (review/LAPLACE-frozen-draws-validation-2026-10-22.md).
A t's tails cover what a Gaussian's cannot, at a cost in effective draws
where the posterior is Gaussian after all. Only a unit's top block takes it;
blocks beneath another keep Gaussian components (`_continuation_importance_rule`
says why).
"""
struct CTSEMImportanceComponent
    mu::Vector{Float64}
    vecs::Matrix{Float64}
    lam::Vector{Float64}
    df::Float64
end

function _continuation_component_logdensity(c::CTSEMImportanceComponent,
    z::AbstractVector{Float64}, work::Vector{Float64})
    k = length(c.mu)
    @inbounds for a in 1:k
        work[a] = z[a] - c.mu[a]
    end
    q = 0.0
    @inbounds for t in 1:k
        s = 0.0
        for a in 1:k
            s += c.vecs[a, t] * work[a]
        end
        q += c.lam[t] * s * s
    end
    isfinite(c.df) || return -(q + k * log(2 * pi)) / 2 + sum(log, c.lam) / 2
    nu = c.df
    return loggamma((nu + k) / 2) - loggamma(nu / 2) - k * log(nu * pi) / 2 +
        sum(log, c.lam) / 2 - (nu + k) / 2 * log1p(q / nu)
end

"""
    _continuation_t_scale(df, z)

The factor a t component's draw scales its Gaussian draw by, `sqrt(df / W)`
with `W` a chi-square on `df` degrees of freedom at the quantile the standard
normal `z` sits at: the draw's extra lattice coordinate, so a t draw is a
Gaussian draw stretched by a quasi-random amount. 1 at `df = Inf`.
"""
function _continuation_t_scale(df::Float64, z::Float64)
    isfinite(df) || return 1.0
    p = clamp(erfc(-z / sqrt(2.0)) / 2, 1e-15, 1 - 1e-15)
    w = 2 * gamma_inc_inv(df / 2, p, 1 - p)
    return sqrt(df / max(w, 1e-300))
end

"""
    _continuation_importance_component(fns, x, floor, df)

The component at block coordinates `x`: the block's precision there
(`fns.precision`) with its eigenvalues floored at `floor`, so a soft or an
indefinite direction is sampled at a bounded width rather than an infinite
one, on `df` degrees of freedom. `nothing` where the precision is not finite.
"""
function _continuation_importance_component(fns, x::Vector{Float64}, floor::Float64,
    df::Float64)
    P = fns.precision(x)
    all(isfinite, P) || return nothing
    E = _ctsem_symeig(_laplace_symmetrise(P))
    all(isfinite, E.values) || return nothing
    return CTSEMImportanceComponent(copy(x), Matrix(E.vectors), max.(E.values, floor),
        df)
end

"""Add the component at `x` unless one sits within 1e-3 of it; whether it was
added."""
function _continuation_add_component!(comps::Vector{CTSEMImportanceComponent}, fns,
    x::Vector{Float64}, floor::Float64, df::Float64)
    any(c -> maximum(abs, c.mu .- x) < 1e-3, comps) && return false
    c = _continuation_importance_component(fns, x, floor, df)
    c === nothing && return false
    push!(comps, c)
    return true
end

"""
    _continuation_newton_point(fns, x, floor; maxstep=3)

One damped Newton step on the block's log integrand from `x`, on its curvature
there with the eigenvalues floored, its length capped at `maxstep` in that
curvature's metric: where a component goes from a draw that carries weight, as
`.ctImisNewtonComponent()` moves the outer sampler's (after optimized
population Monte Carlo, Elvira and Chouzenoux 2022). `nothing` where the value
or the curvature is not finite.
"""
function _continuation_newton_point(fns, x::Vector{Float64}, floor::Float64;
    maxstep::Real=3.0)
    vg = fns.value_gradient(x)
    (isfinite(vg.value) && all(isfinite, vg.gradient)) || return nothing
    c = _continuation_importance_component(fns, x, floor, Inf)
    c === nothing && return nothing
    w = transpose(c.vecs) * vg.gradient
    len = sqrt(sum(abs2, w ./ sqrt.(c.lam)))
    scale = len > maxstep ? maxstep / len : 1.0
    return x .+ c.vecs * (scale .* w ./ c.lam)
end

"""The draws an importance rule has taken: their block coordinates, the
component each came from (0 the prior), the replicate lattice each belongs to
(`_CONTINUATION_REPLICATES` per batch, for the standard error), the log
integrand beneath each -- its members' log likelihood at a leaf, the sum of its
children's log integrals at an outer block -- the children's rules placed at
each draw (none at a leaf), and how many batches have been drawn."""
mutable struct _ContinuationDraws
    points::Matrix{Float64}
    origin::Vector{Int}
    replicate::Vector{Int}
    h::Vector{Float64}
    children::Vector{Vector{CTSEMFixedRule}}
    batches::Int
end

"""
How many independently shifted lattices each batch of an importance rule
draws, per component. The rule's standard error is the spread of their log
integrals over `sqrt` of their number: a shifted lattice's error falls faster
than independent draws', and the independent-draws formula `sqrt(1/ESS - 1/n)`
overstated it 1.7-fold at eight effects and 3- to 10-fold at one or two
(review/LAPLACE-frozen-draws-validation-2026-10-22.md), so batches doubled
needlessly and the screen called resolvable gaps noise. Eight: the spread of
eight has about a quarter's relative error itself, and a 256-draw batch still
gives each replicate 32 points over the components.
"""
const _CONTINUATION_REPLICATES = 8

"""`m` split as evenly as whole numbers allow over `J` parts, the larger
first."""
function _continuation_split(m::Integer, J::Integer)
    J <= 0 && return Int[]
    q, r = divrem(Int(m), Int(J))
    return [q + (t <= r ? 1 : 0) for t in 1:Int(J)]
end

"""`m` draws of block `b` from component `component` (`c`), scored at `theta`
with the blocks above where `u` has them: at a leaf each draw's members' log
likelihood, the deviations their filters integrate integrated exactly, as at a
product rule's node; at an outer block the sum of the log integrals of the
blocks beneath it, each by its own rule placed with this block at the draw
(`_continuation_importance_rule`), kept in `D.children`. In parallel over the
draws, each on its own copy of `u`.

The `m` are `_CONTINUATION_REPLICATES` lattices, each with its own shift from
`key`, the batch, the component and the replicate alone, so the same replicate
of the same component of the same batch draws the same standard normals at
every placement: common random numbers. A t component draws one lattice
coordinate more, which stretches each draw (`_continuation_t_scale`). A block
beneath takes the draw's index into its own `path`, and so into its key: every
draw's children have draws of their own, whose errors are independent across
the draws above, so the spread of this block's replicates counts them."""
function _continuation_draw!(D::_ContinuationDraws, laplace::CTSEMLaplaceObjective,
    U::Integer, theta::Vector{Float64}, Ls, b::Integer, u::Vector{Float64}, context,
    opts, c::CTSEMImportanceComponent, component::Int, m::Int, key,
    path::Vector{Int})
    m > 0 || return D
    block = laplace.units.blocks[U][b]
    columns = (block.offset + 1):(block.offset + block.size)
    kids = context.children[b]
    k = length(c.mu)
    R = _CONTINUATION_REPLICATES
    sizes = _continuation_split(m, R)
    kz = isfinite(c.df) ? k + 1 : k
    Z = Matrix{Float64}(undef, kz, m)
    reps = Vector{Int}(undef, m)
    at = 0
    for r in 1:R
        sizes[r] > 0 || continue
        Z[:, (at + 1):(at + sizes[r])] = _continuation_lattice_normals(kz, sizes[r],
            Random.Xoshiro(hash((key, D.batches, component, r))))
        reps[(at + 1):(at + sizes[r])] .= r
        at += sizes[r]
    end
    # Through the component's eigenvectors signed (`_ctsem_signed_columns!`),
    # so a placement whose eigensolve flips one does not reflect its draws.
    V = _ctsem_signed_columns!(copy(c.vecs))
    X = Matrix{Float64}(undef, k, m)
    for j in 1:m
        stretch = kz > k ? _continuation_t_scale(c.df, Z[kz, j]) : 1.0
        @inbounds for t in 1:k
            s = 0.0
            for a in 1:k
                s += V[t, a] * Z[a, j] / sqrt(c.lam[a])
            end
            X[t, j] = c.mu[t] + stretch * s
        end
    end
    before = size(D.points, 2)
    npar = length(theta)
    hs = fill(-Inf, m)
    placed = [CTSEMFixedRule[] for _ in 1:(isempty(kids) ? 0 : m)]
    score = function (j)
        # `local` throughout, as in `_continuation_flagged`.
        local uj, aws, inner, r
        uj = copy(u)
        @inbounds for (t, col) in enumerate(columns); uj[col] = X[t, j]; end
        aws = _laplace_workspace!(laplace, Float64, npar)
        inner = 0.0
        if isempty(kids)
            for mm in block.members
                inner += _laplace_member_loglik(aws, laplace.objective.subject_objectives[
                    laplace.units.members[U][mm]], _laplace_member_values(theta,
                    laplace.spec, Ls, uj, laplace.units.offsets[U][mm]))
            end
        else
            for child in kids
                r = _continuation_importance_rule(laplace, U, theta, Ls, child, uj,
                    context, aws, opts, vcat(path, before + j))
                # A block beneath that cannot be placed leaves this draw no
                # weight; it is dropped from the rule (`_continuation_importance_rule`).
                if r === nothing
                    inner = -Inf
                    break
                end
                push!(placed[j], r.rule)
                inner += r.value
            end
        end
        hs[j] = isfinite(inner) ? inner : -Inf
        return true
    end
    # Below 64 member evaluations on this thread, as `_continuation_block_rule`
    # places its nodes.
    m * length(block.members) >= 64 ? _laplace_parallel(score, laplace, 1:m) :
        foreach(score, 1:m)
    D.points = hcat(D.points, X)
    append!(D.origin, fill(component, m))
    append!(D.replicate, reps)
    append!(D.h, hs)
    append!(D.children, placed)
    return D
end

"""A batch of `m` draws: the prior's share, the rest split over the
components."""
function _continuation_batch!(D::_ContinuationDraws, laplace, U, theta, Ls, b, u,
    context, opts, comps, prior::CTSEMImportanceComponent, share::Float64, m::Int,
    key, path::Vector{Int})
    D.batches += 1
    mp = clamp(round(Int, share * m), 1, m)
    _continuation_draw!(D, laplace, U, theta, Ls, b, u, context, opts, prior, 0, mp,
        key, path)
    for (t, n) in enumerate(_continuation_split(m - mp, length(comps)))
        _continuation_draw!(D, laplace, U, theta, Ls, b, u, context, opts, comps[t],
            t, n, key, path)
    end
    return D
end

"""
The log weights of every draw under the balance heuristic -- each component
weighted by its share of the draws -- and the log integral, its standard error
and effective sample size, for a block of `k` effects. The standard error is
the spread of the replicate lattices' own log integrals (each weighted over
its own draws, by its own counts) over `sqrt` of their number
(`_CONTINUATION_REPLICATES`); where fewer than two replicates score, the
independent-draws formula `sqrt(1/ESS - 1/n)` stands in. `replw` is each
draw's log weight within its replicate's estimate.
"""
function _continuation_importance_weigh(D::_ContinuationDraws,
    comps::Vector{CTSEMImportanceComponent}, prior::CTSEMImportanceComponent, k::Int)
    n = size(D.points, 2)
    J = length(comps)
    dens = Matrix{Float64}(undef, J + 1, n)
    work = zeros(k)
    @inbounds for j in 1:n
        z = view(D.points, :, j)
        for t in 0:J
            dens[t + 1, j] = _continuation_component_logdensity(t == 0 ? prior : comps[t],
                z, work)
        end
    end
    # The log weights of a set of draws under the balance heuristic over the
    # set's own counts, and the set's log integral.
    function weigh(sel)
        m = length(sel)
        counts = zeros(Int, J + 1)
        for j in sel
            counts[D.origin[j] + 1] += 1
        end
        lw = Vector{Float64}(undef, m)
        logq = Vector{Float64}(undef, J + 1)
        @inbounds for (i, j) in enumerate(sel)
            for t in 1:(J + 1)
                logq[t] = counts[t] == 0 ? -Inf : log(counts[t] / m) + dens[t, j]
            end
            z = view(D.points, :, j)
            lw[i] = -dot(z, z) / 2 - _logsumexp(logq) - log(m)
        end
        terms = lw .+ D.h[sel]
        peak = maximum(terms)
        isfinite(peak) || return (lw=lw, terms=terms, value=NaN, ess=0.0)
        w = exp.(terms .- peak)
        sw = sum(w)
        return (lw=lw, terms=terms, value=peak + log(sw) - k * log(2 * pi) / 2,
            ess=sw^2 / sum(abs2, w))
    end
    all_ = weigh(1:n)
    replw = fill(-Inf, n)
    isfinite(all_.value) || return (lw=all_.lw, terms=all_.terms, value=NaN, se=Inf,
        ess=0.0, replw=replw)
    values = Float64[]
    for r in 1:_CONTINUATION_REPLICATES
        sel = findall(==(r), D.replicate)
        isempty(sel) && continue
        part = weigh(sel)
        replw[sel] = part.lw
        isfinite(part.value) && push!(values, part.value)
    end
    se = length(values) >= 2 ?
        sqrt(sum(abs2, values .- sum(values) / length(values)) /
             (length(values) - 1) / length(values)) :
        sqrt(max(1 / all_.ess - 1 / n, 0.0))
    return (lw=all_.lw, terms=all_.terms, value=all_.value, se=se, ess=all_.ess,
        replw=replw)
end

"""
    _continuation_outer_components(M, blocks, children, mode, floor, df)

The first component of every outer block's proposal, `nothing` for a leaf: a
t on `df` degrees of freedom for a block with nothing above it, a Gaussian
beneath another (`_continuation_importance_rule` says why), at the unit's
mode on the block's precision once every block beneath
it is eliminated, its eigenvalues floored at `floor`. The elimination is the
block factorization's (`_laplace_block_factor`, innermost first), each block
eliminated with its precision's eigenvalues clipped at the prior's, as the
quadrature clips a rule's scale (`_quadrature_clipped_scale`): a block's
posterior covariance is no wider than its prior's where its likelihood is
log-concave, while a soft curvature's inverse can be any width. Unclipped, the
eliminated precision -- what a product rule scales an outer block by -- is
where a near-singular child reappears (the SNSF pilot's subject 1174:
`CT-SEM/review/LAPLACE-nested-spike-2026-10-02.md`, section 2), and the
study's proposal is as wide as the floor allows along that direction. Where
`M - I` is positive definite nothing is clipped and this is the eliminated
precision exactly. `nothing` throughout where a block is not finite.
"""
function _continuation_outer_components(M::CTSEMBlockMatrix{Float64},
    blocks::Vector{CTSEMLaplaceBlock}, children, mode::Vector{Float64},
    floor::Float64, df::Float64)
    nb = length(blocks)
    out = Vector{Union{Nothing,CTSEMImportanceComponent}}(nothing, nb)
    diag = [_laplace_symmetrise(Matrix(d)) for d in M.diag]
    coupling = [[Matrix(c) for c in row] for row in M.coupling]
    for b in 1:nb
        all(isfinite, diag[b]) || return fill!(out, nothing)
        E = _ctsem_symeig(_laplace_symmetrise(diag[b]))
        all(isfinite, E.values) || return fill!(out, nothing)
        if !isempty(children[b])
            columns = (blocks[b].offset + 1):(blocks[b].offset + blocks[b].size)
            out[b] = CTSEMImportanceComponent(mode[columns], Matrix(E.vectors),
                max.(E.values, floor), isempty(blocks[b].ancestors) ? df : Inf)
        end
        ancestors = blocks[b].ancestors
        isempty(ancestors) && continue
        inverse = E.vectors * Diagonal(1 ./ max.(E.values, 1.0)) * transpose(E.vectors)
        W = [inverse * coupling[b][t] for t in eachindex(ancestors)]
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
    return out
end

"""
    _continuation_moment_component!(comps, D, current, floor, df)

Add the component matched to the weighted draws: their weighted mean, and
their weighted covariance shrunk towards the first component's scale by
`k + 2` draws against the draws' effective number (the fewest that keep it
proper whatever the weights), as its scale matrix, its precision's
eigenvalues floored at `floor`, on `df` degrees of freedom. How an outer
block's proposal adapts (adaptive multiple importance sampling, Cornuet et al.
2012): its marginal density, with the blocks beneath integrated out, has no
gradient or curvature to take a Newton step on, only the draws' weights.
Whether it was added.
"""
function _continuation_moment_component!(comps::Vector{CTSEMImportanceComponent},
    D::_ContinuationDraws, current, floor::Float64, df::Float64)
    k = size(D.points, 1)
    peak = maximum(current.terms)
    isfinite(peak) || return false
    w = exp.(current.terms .- peak)
    sw = sum(w)
    mu = (D.points * w) ./ sw
    X = D.points .- mu
    S = (X .* transpose(w)) * transpose(X) ./ sw
    first = comps[1]
    S0 = first.vecs * Diagonal(1 ./ first.lam) * transpose(first.vecs)
    S = (current.ess .* S .+ (k + 2) .* S0) ./ (current.ess + k + 2)
    E = _ctsem_symeig(_laplace_symmetrise(S))
    (all(isfinite, E.values) && all(>(0), E.values)) || return false
    push!(comps, CTSEMImportanceComponent(mu, Matrix(E.vectors),
        max.(1 ./ E.values, floor), df))
    return true
end

"""
    _continuation_draw_target(block, context, opts)

The standard error block `block`'s draws aim at and its first batch: for a
block with nothing above it `opts.draws_se` and `opts.draws`; beneath another,
its share of `draws_se` over its siblings -- their squared errors together at
most `draws_se`, each aiming at `sqrt(draws_se / siblings)`, never below
`draws_se` itself -- starting at the fewest draws that could reach it,
`(share / (1 - share)) / se^2`, the error of a proposal matching the block
exactly but for its prior share, at least 8 and at most `opts.draws`. See
`_continuation_importance_rule`. Every block keeps drawing until both its
replicates' spread (`_continuation_importance_weigh`, what it reports and the
screen reads) and its weights' spread, `sqrt(1/ESS - 1/n)`, are under that
target.

Why no looser, and why the weights too. A block's estimate is unbiased, but
its draws come from Gaussians lighter-tailed than these posteriors, and at a
few hundred draws it falls short: of the mass in the tails, which every
replicate misses alike, so their spread does not show the shortfall, while
the weights' spread does. At the top of a unit the same: on the bench's gA1
(review/LAPLACE-frozen-draws-validation-2026-10-22.md), stopped on the
replicates alone at about 1,100 draws a unit, the draws' total error against
exact was -1.2 to -2.1 nats with a reported error 0.72 of the truth, where
driven to 16,000 draws it was within 0.13. On the SNSF pilot's study unit at
it290 (`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`, section
10), with the study at 256 draws and every person's draws to 0.091, 0.04 and
0.02 by their weights, the unit came to -16745.691, -16745.657 and
-16745.654, against -16745.649 by Gauss-Hermite at both levels: a deficit of
about a sixth of the children's summed squared error. Tied to `draws_se`, it
stays under a sixth of the error the rule aims at. The first budget, a fixed
0.25 (what the prior's fifth of every batch costs any rule's weights), cost
that unit 0.04. Children stopped on their replicates' spread at 0.041 took
about 350 draws each, against 966, and left the unit 0.10 short, with a
standard error of 0.019.
"""
function _continuation_draw_target(block::CTSEMLaplaceBlock, context, opts)
    isempty(block.ancestors) && return (se=Float64(opts.draws_se), first=Int(opts.draws))
    siblings = length(context.children[block.ancestors[1]])
    se = max(Float64(opts.draws_se), sqrt(Float64(opts.draws_se) / siblings))
    share = Float64(opts.draws_prior)
    least = ceil(Int, share / (1 - share) / se^2)
    return (se=se, first=clamp(least, min(8, Int(opts.draws)), Int(opts.draws)))
end

"""
    _continuation_importance_rule(laplace, U, theta, Ls, b, u, context, aws,
        opts, path=Int[])

Block `b`'s frozen importance rule, and every block's beneath it (see the
module docstring): draws held fixed in the standardised effects, weighted by
the mixture they came from, so the block's term is

    log (1/n) sum_j L(z_j) exp(-z_j'z_j/2) / q(z_j) - k/2 log 2 pi

-- a `CTSEMFixedRule` whose log weights carry `-z'z/2 - log q - log n`, valued
and differentiated as a product rule's nodes are. At a leaf `L(z_j)` is its
members' likelihood. At an outer block it is the product of its children's
integrals, each by a rule of its own placed by this function with the block at
`z_j` (`children[j]`), as `_continuation_block_rule` places product rules: each
integral is one block's effects, never the whole unit's, since importance
sampling jointly over a nested unit's 93 effects would compound the proposal's
mismatch over every direction
(`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`, 3.2). `q` is the
balance heuristic over its components, each weighted by its share of the
draws:

  - at a leaf, a Gaussian at the block's mode on the block's precision there,
    its eigenvalues floored at `opts.draws_floor`: the unit's mode
    (`context.mode`) for a block with nothing above it, and beneath a draw of
    the blocks above the conditional mode there, re-solved from the unit's as a
    product rule's leaf is (`_quadrature_leaf_rule!`). At an outer block, a
    Gaussian at the unit's mode on its precision with the blocks beneath
    eliminated, their curvatures clipped at the prior's, floored
    (`_continuation_outer_components`, placed by `_continuation_unit_rule` as
    `context.outer`);
  - one at every point `opts.points` holds for the unit where this block is
    not at the unit's mode, on the block's precision there at a leaf, on the
    outer shape at an outer block;
  - the prior `N(0, I)`, `opts.draws_prior` of every batch, which bounds every
    weight by the likelihood's supremum over that share, so the estimate has a
    finite variance in any dimension (Hesterberg 1995);
  - where the standard error is not yet reached: at a leaf, where the weights
    show a gap -- a draw carrying more than `n / ESS*` times its share, `ESS* =
    1 / draws_se^2` -- a component at each of up to three such draws, moved a
    damped Newton step and shaped by the precision there
    (`_continuation_newton_point`); at an outer block, the Gaussian matched to
    the weighted draws (`_continuation_moment_component!`).

Draws come in batches, each later one as large as all before it,
quasi-random (`_continuation_lattice_normals`), until both the standard error
-- the spread of the replicate lattices' log integrals
(`_CONTINUATION_REPLICATES`, `_continuation_importance_weigh`) -- and the
weights' spread, which sees a tail the proposal misses where every replicate
misses it alike, are under a target (`_continuation_draw_target`), or
`opts.draws_max` are spent. A block with nothing above it starts at `opts.draws` and aims at
`opts.draws_se`. An outer block's standard error counts its children's errors
as well: they are independent across its draws (each draw's children are keyed
by it), so they spread its replicates' integrals as its own proposal's
mismatch does. And the product of its children's estimates, each unbiased for
its integral, is unbiased for their product, so the block's estimate of its
integral is unbiased as a one-block unit's is. So a block beneath another need
not reach `draws_se` itself: its siblings' squared errors together at most
`draws_se`, each aiming at `sqrt(draws_se / siblings)` (never below
`draws_se`) and starting at the fewest draws that could reach it
(`_continuation_draw_target`, which says why the share is no larger). Children
driven to `draws_se` each made an outer draw of the SNSF pilot's study unit
cost 30 people at 256 draws or more.

Each replicate's lattice in each batch is seeded by the unit, the block, its
`path` (the draws above it), the batch, the component and the replicate, the
same at every placement: common random numbers, so two placements differ by
where their proposals sit rather than by their draws' luck. With fresh draws at
each placement the correction's rounds (removed 2026-10-30), which kept a
step only if the objective with the rule placed again was no lower, were
rejected by the draws' error on the round-9 fit: 0.4 to 0.7 nats at radii down to 0.03 standard
errors.

`u` holds the blocks above where the caller placed them and is not changed.
Returns `(rule, value, se, ess, components)`, or `nothing` when no draw scores
or a mode cannot be found. A draw of an outer block whose children could not
all be placed has no weight, and is left out of the rule.
"""
function _continuation_importance_rule(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, context, aws, opts, path::Vector{Int}=Int[])
    block = laplace.units.blocks[U][b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    leaf = isempty(context.children[b])
    floor = Float64(opts.draws_floor)
    # A t only at the top of a unit. A block beneath another takes small
    # batches and stops on its weights' spread, and with a t's variable
    # weights -- large where the draws sit in the core, small in the tails it
    # samples -- that stopping picks the batches whose draws sit where the
    # weights are largest: on the linear nested fixtures, where Laplace is
    # exact, t leaves read 0.03 to 0.12 nats high at 3 to 16 standard errors,
    # unchanged by more draws, and Gaussian ones within their error. The top
    # block's first batch is 256 and its replicates' spread rules, where the
    # same selection is slight; and the SNSF study unit's t top block doubled
    # its draws for a stiff, near-Gaussian block (section 11 of the shoulders
    # note), which is the t's price there.
    df = isempty(path) ? Float64(opts.draws_df) : Inf
    joint = Float64[context.mode[c] for c in columns]
    comps = CTSEMImportanceComponent[]
    fns = nothing
    if leaf
        centre = joint
        if !isempty(block.ancestors)
            lr = _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws; start=joint)
            lr.ok || return nothing
            centre = lr.centre
        end
        fns = _continuation_block_functions(laplace, U, theta, Ls, b, u, aws)
        _continuation_add_component!(comps, fns, centre, floor, df)
    else
        first = context.outer[b]
        first === nothing && return nothing
        push!(comps, first)
    end
    isempty(comps) && return nothing
    for p in get(opts.points, Int(U), Vector{Vector{Float64}}())
        length(p) == length(u) || continue
        x = Float64[p[c] for c in columns]
        (all(isfinite, x) && maximum(abs, x .- joint) >= 1e-3) || continue
        if leaf
            _continuation_add_component!(comps, fns, x, floor, df)
        elseif !any(c -> maximum(abs, c.mu .- x) < 1e-3, comps)
            push!(comps, CTSEMImportanceComponent(x, comps[1].vecs, comps[1].lam, df))
        end
    end
    # The prior is the integrand's own Gaussian, kept as it is.
    prior = CTSEMImportanceComponent(zeros(k), Matrix{Float64}(I, k, k), ones(k), Inf)
    # A block with nothing above it keeps the key it had before blocks could
    # nest, so a one-block unit draws what it drew then.
    key = isempty(path) ? (Int(opts.draws_seed), Int(U), Int(b)) :
        (Int(opts.draws_seed), Int(U), Int(b), path)
    D = _ContinuationDraws(zeros(k, 0), Int[], Int[], Float64[],
        Vector{CTSEMFixedRule}[], 0)
    share = Float64(opts.draws_prior)
    target = _continuation_draw_target(block, context, opts)
    _continuation_batch!(D, laplace, U, theta, Ls, b, u, context, opts, comps, prior,
        share, target.first, key, path)
    current = _continuation_importance_weigh(D, comps, prior, k)
    target_ess = 1 / target.se^2
    # Every block stops on the larger of its replicates' spread and its
    # weights' (`_continuation_draw_target` says why); the replicates' is what
    # it reports and the screen reads.
    while (!isfinite(current.value) ||
            max(current.se, sqrt(max(1 / current.ess - 1 / size(D.points, 2), 0.0))) >
            target.se) && size(D.points, 2) < opts.draws_max
        n = size(D.points, 2)
        if isfinite(current.value)
            if leaf
                peak = maximum(current.terms)
                w = exp.(current.terms .- peak)
                ratio = w ./ (sum(w) / n)
                added = 0
                for j in sortperm(ratio; rev=true)
                    (added >= 3 || ratio[j] <= n / target_ess) && break
                    x = _continuation_newton_point(fns, D.points[:, j], floor)
                    x === nothing && continue
                    _continuation_add_component!(comps, fns, x, floor, df) &&
                        (added += 1)
                end
            else
                _continuation_moment_component!(comps, D, current, floor, df)
            end
        end
        _continuation_batch!(D, laplace, U, theta, Ls, b, u, context, opts, comps,
            prior, share, min(n, Int(opts.draws_max) - n), key, path)
        current = _continuation_importance_weigh(D, comps, prior, k)
    end
    isfinite(current.value) || return nothing
    points, logweights, children = D.points, current.lw, D.children
    replicate, replw = D.replicate, current.replw
    if !leaf
        keep = findall(isfinite, D.h)
        if length(keep) < length(D.h)
            points = points[:, keep]
            logweights = logweights[keep]
            children = children[keep]
            replicate = replicate[keep]
            replw = replw[keep]
        end
    end
    rule = CTSEMFixedRule(b, points, logweights, -k * log(2 * pi) / 2, children,
        replicate, replw)
    return (rule=rule, value=current.value, se=current.se, ess=current.ess,
        components=length(comps))
end

"""
    _continuation_block_rule(laplace, U, theta, Ls, b, u, context, aws, opts)

Block `b`'s fixed rule at `theta`, its ancestors held where `u` has them, and
its log integral there: the block's point set (`context.sets[b]`,
`point_sets.jl`) placed at a leaf's conditional mode and clipped scale
(`_quadrature_leaf_rule!`), or at an outer block at the unit's joint mode and
the eliminated precision's clipped scale (`_quadrature_clipped_scale`), as the
placed objective places them (`laplace_placed.jl`), so that at the centre the
two agree. A node's value is its members' log likelihood at a leaf, and the
sum of its children's rules placed with this block at the node at an outer
block. Returns `(rule, value, se)`, `se` the rule's error: a lattice rule's
shifts' spread and, at an outer block, its children's errors at its nodes in
quadrature, weighted as the nodes are (`_placed_nested_unit!` says why);
`nothing` where the rule cannot be placed.
"""
function _continuation_block_rule(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, context, aws, opts)
    blocks = laplace.units.blocks[U]
    block = blocks[b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    kids = context.children[b]
    leaf = isempty(kids)
    rule = if leaf
        _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws;
            start=Float64[context.mode[c] for c in columns])
    else
        clipped = _quadrature_clipped_scale(context.factors[b], k)
        (ok=true, centre=Float64[context.mode[c] for c in columns],
         scale=clipped.scale, logdetscale=clipped.logdetscale)
    end
    rule.ok || return nothing
    ps = context.sets[b]
    S = Matrix{Float64}(rule.scale)
    pl = _point_set_place(ps, rule.centre, S, Matrix(inv(UpperTriangular(S))),
        Float64(rule.logdetscale))
    points = pl.Z
    n = size(points, 2)
    logweights = Vector{Float64}(undef, n)
    for j in 1:n
        logweights[j] = pl.logw[j] - sum(abs2, view(points, :, j)) / 2
    end
    # A lattice rule's value is its ratio against its Gaussian
    # (`_point_set_gaussian`), whose sum is a constant of the points once
    # placed, as the placed objective reads it.
    gs = ps.family === :lattice ? _point_set_gaussian(pl,
        Matrix(inv(UpperTriangular(S))), ps.nrep) : nothing
    constant = gs === nothing ? -k * log(2 * pi) / 2 :
        Float64(rule.logdetscale) - gs.peak - log(gs.sum)
    # In parallel over the nodes, each on its own copy of `u`, as an importance
    # rule's draws are (`_continuation_draw!`): a nested unit whose outer block
    # holds 125 nodes, each placing every member's rule beneath it, was placed
    # on one thread whenever it was the only unit being placed. Not below 64
    # member evaluations (`_continuation_ranges`' bar): a region opened for
    # every small block made test_continuation.jl three times slower.
    children = [CTSEMFixedRule[] for _ in 1:(leaf ? 0 : n)]
    terms = fill(-Inf, n)
    childse = zeros(n)
    npar = length(theta)
    place_node = function (j)
        # `local` throughout, as in `_continuation_flagged`.
        local uj, awsj, inner, placed, v
        uj = copy(u)
        @inbounds for (t, c) in enumerate(columns); uj[c] = points[t, j]; end
        awsj = _laplace_workspace!(laplace, Float64, npar)
        inner = 0.0
        v = 0.0
        if leaf
            for m in block.members
                inner += _laplace_member_loglik(awsj, laplace.objective.subject_objectives[
                    laplace.units.members[U][m]], _laplace_member_values(theta,
                    laplace.spec, Ls, uj, laplace.units.offsets[U][m]))
            end
        else
            for c in kids
                placed = _continuation_block_rule(laplace, U, theta, Ls, c, uj,
                    context, awsj, opts)
                placed === nothing && return false
                push!(children[j], placed.rule)
                inner += placed.value
                v += placed.se^2
            end
        end
        terms[j] = isfinite(inner) ? inner + logweights[j] : -Inf
        childse[j] = sqrt(v)
        return true
    end
    placed_all = n * length(block.members) >= 64 ?
        _laplace_parallel(place_node, laplace, 1:n) : all(place_node, 1:n)
    placed_all || return nothing
    peak = maximum(terms)
    isfinite(peak) || return nothing
    acc = 0.0
    weighted = 0.0
    repsums = zeros(ps.nrep)
    for j in 1:n
        w = exp(terms[j] - peak)
        acc += w
        weighted += w * childse[j]
        ps.nrep > 1 && (repsums[pl.rep[j]] += w)
    end
    value = peak + log(acc) + constant
    se = sqrt(_point_set_error(gs === nothing ? repsums : repsums ./ gs.repsums,
        ps.nrep)^2 + (weighted / acc)^2)
    replicate = ps.nrep > 1 ? copy(pl.rep) : Int[]
    replw = ps.nrep > 1 ? logweights .+ log(ps.nrep) : Float64[]
    return (rule=CTSEMFixedRule(b, points, logweights, constant, children, replicate,
        replw), value=value, se=se)
end

"""
    _continuation_unit_rule(laplace, U, theta, Ls, aws, opts)

Unit `U`'s fixed rule at `theta` and its quadrature value there, placed as
`_quadrature_chunk!` places the adaptive rule: the unit's mode from the
origin, its curvature, the block factorization, the recursion over the tree,
each block on the point set of the unit's rule (`opts.choice` where given,
else `opts.unit_rules[U]`; `_rule_sets`), or on importance draws where the unit
is named for them. Returns `(rule, value, se, importance)`; `nothing` for a
unit with no random effects, or when the rule cannot be placed.
"""
function _continuation_unit_rule(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, aws, opts)
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return nothing
    _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    mode = copy(laplace.modes[U])
    tree = _quadrature_children(blocks)
    if opts.draws_all || U in get(opts, :draw_units, ())
        # Every unit on draws (the random-effects check), or named for draws
        # because no rule reached the draws' target there
        # (`_continuation_refine!`). A nested unit's outer blocks take their
        # first components from the curvature, here, before anything runs in
        # parallel: it is this slot's scratch.
        # Assigned once: the check below is a closure over it, and a closure
        # over a local assigned twice boxes it.
        outer = length(blocks) > 1 ?
            _continuation_outer_components(
                _laplace_unit_curvature(laplace, U, theta, Ls, mode), blocks,
                tree.children, mode, Float64(opts.draws_floor),
                Float64(opts.draws_df)) :
            Vector{Union{Nothing,CTSEMImportanceComponent}}(nothing, length(blocks))
        length(blocks) > 1 &&
            any(b -> !isempty(tree.children[b]) && outer[b] === nothing,
                eachindex(blocks)) && return nothing
        context = (children=tree.children, mode=mode, outer=outer)
        roots = CTSEMFixedRule[]
        total = 0.0
        variance = 0.0
        for root in tree.roots
            placed = _continuation_importance_rule(laplace, U, theta, Ls, root,
                copy(mode), context, aws, opts)
            placed === nothing && return nothing
            push!(roots, placed.rule)
            total += placed.value
            variance += placed.se^2
        end
        evaluations = sum(_continuation_evaluations(r, blocks) for r in roots; init=0)
        return (rule=CTSEMFixedUnitRule(U, roots, evaluations), value=total,
            se=sqrt(variance), importance=true)
    end
    M = _laplace_unit_curvature(laplace, U, theta, Ls, mode)
    factorization = _laplace_factor_repaired!(M, blocks)
    factorization.ok || return nothing
    choice = haskey(opts, :choice) ? opts.choice : opts.unit_rules[U]
    context = (children=tree.children, mode=mode, factors=factorization.factors,
        sets=_rule_sets(laplace, U, choice))
    u = copy(mode)
    roots = CTSEMFixedRule[]
    total = 0.0
    variance = 0.0
    for root in tree.roots
        placed = _continuation_block_rule(laplace, U, theta, Ls, root, u, context,
            aws, opts)
        placed === nothing && return nothing
        push!(roots, placed.rule)
        total += placed.value
        variance += placed.se^2
    end
    isfinite(total) || return nothing
    evaluations = sum(_continuation_evaluations(r, blocks) for r in roots; init=0)
    return (rule=CTSEMFixedUnitRule(U, roots, evaluations), value=total,
        se=sqrt(variance), importance=false)
end

"""
    _continuation_place(laplace, theta, opts)

Every unit's rule and quadrature value at `theta`, and every unit's Laplace
term there, in parallel over units. A unit whose rule could not be placed has
`nothing` for a rule and a `NaN` value; a unit with no random effects has no
rule and its plain log likelihood, which is what both terms are there.
"""
function _continuation_place(laplace::CTSEMLaplaceObjective, theta::Vector{Float64},
    opts)
    nunits = length(laplace.units.members)
    lap = ctsem_laplace_evaluate(laplace, theta; gradient=false)
    Ls = _laplace_popchols(theta, laplace.spec)
    _laplace_ensure_pool!(laplace)
    rules = Vector{Union{Nothing,CTSEMFixedUnitRule}}(nothing, nunits)
    values = fill(NaN, nunits)
    ses = zeros(nunits)
    importance = fill(false, nunits)
    # Each unit's rule, which the refine may change (`_continuation_refine!`).
    choices = copy(opts.unit_rules)
    only = get(opts, :only, nothing)
    _laplace_parallel(laplace, 1:nunits) do U
        (only === nothing || U in only) || return true
        placed = try
            _continuation_unit_rule(laplace, U, theta, Ls,
                _laplace_workspace!(laplace, Float64, length(theta)), opts)
        catch err
            # A point the model cannot evaluate places no rule; code that is
            # wrong is an error (`_ctsem_must_propagate`).
            _ctsem_must_propagate(err) && rethrow()
            nothing
        end
        if placed === nothing
            isempty(laplace.units.blocks[U]) && (values[U] = lap.unit_loglik[U])
        else
            rules[U] = placed.rule
            values[U] = placed.value
            ses[U] = placed.se
            importance[U] = placed.importance
        end
        return true
    end
    return (rules=rules, quadrature=values, laplace=copy(lap.unit_loglik),
        laplace_value=lap.value, converged=lap.converged, se=ses,
        importance=importance, choices=choices)
end

"""
The most nodes a product rule is refined to (`_continuation_refine!`): two
steps of two above the default 5, so 81 member evaluations per gradient on a
two-effect block against 25. Past that a unit whose rules have not converged
is one whose posterior over its effects is not one basin the nodes can cover.
"""
const _CONTINUATION_NODES_MAX = 9

"""
The most member evaluations a product rule over a unit's effects may cost a
member before the refine leaves it out (`_continuation_refine!`): four times
the lattice's most points (`_LATTICE_MAX` shifts' worth), so the product
family is tried wherever it could still be the cheaper one at the same error,
and never where its next step alone would dwarf the lattice's last.
"""
const _CONTINUATION_PRODUCT_COST_MAX = 4 * _LATTICE_SHIFTS * _LATTICE_MAX

"""
    _continuation_steps(laplace, U, opts)

The rules the refine may take for unit `U`, family by family, each a list of
`CTSEMRuleChoice` in the order the refine raises them
(`_continuation_refine!`): for a unit of one block, the product rule at
`nodes - 2`, `nodes`, ... to `_CONTINUATION_NODES_MAX` (the first only as the
check of the second) and the lattice at `_LATTICE_START` to `_LATTICE_MAX`
points a shift, doubling; for a tree of blocks, the two raised together
(`:mixed`, each block on the cheaper, `_rule_block_product`). A product
step past `_CONTINUATION_PRODUCT_COST_MAX` a member is left out, and the
product family with it if fewer than two steps are left.
"""
function _continuation_steps(laplace::CTSEMLaplaceObjective, U::Integer, opts)
    blocks = laplace.units.blocks[U]
    nm = sum((length(b.members) for b in blocks); init=0)
    nm = max(nm, 1)
    nodes = Int(opts.nodes)
    families = Vector{Vector{CTSEMRuleChoice}}()
    if length(blocks) == 1
        product = CTSEMRuleChoice[]
        for n in (nodes - 2):2:_CONTINUATION_NODES_MAX
            c = CTSEMRuleChoice(:product, n, 0)
            _rule_cost(laplace, U, c) / nm <= _CONTINUATION_PRODUCT_COST_MAX || break
            push!(product, c)
        end
        length(product) >= 2 && push!(families, product)
        lattice = CTSEMRuleChoice[]
        m = _LATTICE_START
        while m <= _LATTICE_MAX
            push!(lattice, CTSEMRuleChoice(:lattice, nodes, m))
            m *= 2
        end
        push!(families, lattice)
    else
        mixed = CTSEMRuleChoice[]
        n = nodes - 2
        m = _LATTICE_START ÷ 2
        while n <= _CONTINUATION_NODES_MAX && m <= _LATTICE_MAX
            c = CTSEMRuleChoice(:mixed, n, m)
            _rule_cost(laplace, U, c) / nm <= _CONTINUATION_PRODUCT_COST_MAX || break
            push!(mixed, c)
            n += 2
            m *= 2
        end
        isempty(mixed) || push!(families, mixed)
    end
    return families
end

"""
    _continuation_product_error(values, s)

Product step `s`'s error and the step before's, from the geometric tail its
last two node steps imply: with `d` the last step's change and `dp` the one
before, `r = d / dp`, the finer rule's error is `d max(1, r / (1 - r))` --
never less than `d`, one ratio being a noisy estimate of the next -- and the
coarser's `d` more. `nothing` before three steps or where the steps do not
shrink, unless both are within 1e-9 nats (an integrand the rules hold
exactly). One step read alone understated a slowly converging rule about
threefold: on gA14's one-effect units, whose integrands have heavy tails, 5
nodes stood 0.045 nats from the integral where their step to 7 said 0.017,
and the forty units' errors, all of one sign, summed to 2.2 nats
(`CT-SEM/scripts/unit-sampling-2026-10/unit_smc.R`).
"""
function _continuation_product_error(values::Vector{Float64}, s::Int)
    s >= 3 || return nothing
    (isfinite(values[s]) && isfinite(values[s - 1]) && isfinite(values[s - 2])) ||
        return nothing
    d = abs(values[s] - values[s - 1])
    dp = abs(values[s - 1] - values[s - 2])
    max(d, dp) <= 1e-9 && return (fine=d, coarse=d + dp)
    d < dp || return nothing
    r = d / dp
    fine = d * max(1.0, r / (1 - r))
    return (fine=fine, coarse=d + fine)
end

"""
    _continuation_count_step(values, s)

How far lattice step `s`'s value lies from the step before's, at half the
count: `Inf` where that step was not evaluated or failed. A lattice's shifts'
spread cannot see a tail every shift misses, and the count step does: on the
round-9 stall's unit 30 (eight effects, a shoulder beside a stiff mode) the
256-point lattice read 0.21 nats low with a spread of 0.013, and 512 points
read 0.2 higher, where a point reached the shoulder
(`CT-SEM/scripts/placed-mixture-2026-10/r9_placed.R`).
"""
function _continuation_count_step(values::Vector{Float64}, s::Int)
    (s >= 2 && isfinite(values[s - 1]) && isfinite(values[s])) || return Inf
    return abs(values[s] - values[s - 1])
end

"""
    _continuation_refine!(laplace, theta, opts, placed, flagged)

Each flagged unit's rule, chosen by counted cost for the accuracy, and every
flagged unit's rule error at the centre, in nats. What
`ctsem_laplace_continuation_info` sums into the hybrid's own error, the one
number a fit reports for how far its objective is from the exact marginal.

Each family of rule (`_continuation_steps`) is a sequence of counts. A
lattice step's error is the larger of its shifts' spread and its count step,
the difference between its value and the half count's
(`_continuation_count_step`); a product step's is the geometric tail its last
two node steps imply (`_continuation_product_error`), held to `opts.draws_se`
over the root of the flagged units, since product rules' errors share a sign
and add; a mixed step's is the difference between its value and the step
before's. Where the step after is within its target of a step, that step
stands with the finer rule's measure of it as its error. The refine evaluates,
at each turn, the family whose next step costs the fewest member evaluations
(`_rule_cost`), and stops at the first step whose error is within
`opts.draws_se`: the cheapest rule that reaches the target, at no more than a
few times its cost in rules tried. So which family a unit takes is decided by
what each costs it at that accuracy, never by its width: product
Gauss-Hermite where `nodes^k` is cheaper at the same error, which on the
fixtures was up to three effects, lattice-t past that. The screen's own rule
(`opts.unit_rules`, `_rule_default`) counts as already evaluated.

Where no step reaches the target, a product or mixed family whose last
difference `d` is smaller than the one before, `d_prev`, is closing on one
basin: it keeps its finest rule and reports as its error the geometric tail
that ratio implies, `d r / (1 - r)` with `r = d / d_prev`, and never less than
`d`, since one ratio is a noisy estimate of the next -- or the lattice's
finest stands, with its error, where that is smaller. A unit whose product
rules do not converge and whose lattice does not reach the target at its most
points is not one basin a placed rule can follow: it takes draws, which stop
at `draws_se` on their own standard error (`_continuation_importance_rule`),
or keeps whichever family's finest rule has the smaller error where draws are
off. Units on draws
already keep them. `placed`'s arrays are updated in place, its `choices` with
each unit's rule, which the continuation keeps for every later placement
(`o.unit_rules`).

Why the error and not draws: a unit whose rules converge has a known error,
which the fit reports, and draws cost far more than the rule -- a nested
unit's are conditional per block, so a member pays the product of the draws
along its path. On gN3 at its refined estimate the 9-node rule of unit 3
differed from the 7-node by 0.059 nats against a target of 0.05 (0.246 before
it, so converging), and the draws it took instead were 353,781: 364,440
member evaluations a value against 12,630 for the other seven units, 2,301 s
for the refined fit and its interval check against 68 s without them, at the
same estimate to 0.000 standard errors
(scripts/fit-review-2026-10/race_roundstops.R, refined_rules.R).

Why the refine at all: a hybrid is only as honest as its rules. On the bench's
gC8 the 5-node rules carried 6.9 nats of node error at the climb's end, and on
anomS2's spurious start 54 nats over 101 units: the climb on that hybrid
walked where the rules over-credited the Laplace basin and ended 11 and 15
nats below the baseline, certified on an objective whose own error was that
large.
"""
function _continuation_refine!(laplace::CTSEMLaplaceObjective, theta::Vector{Float64},
    opts, placed, flagged::Vector{Int})
    nunits = length(laplace.units.members)
    errors = fill(NaN, nunits)
    isempty(flagged) && return errors
    target = Float64(opts.draws_se)
    Ls = _laplace_popchols(theta, laplace.spec)
    todo = Int[]
    for U in flagged
        if placed.importance[U]
            errors[U] = placed.se[U]
        elseif placed.rules[U] !== nothing
            push!(todo, U)
        end
    end
    # A product rule's error has one sign at every unit, so the units' errors
    # add rather than average: each is held to the target over the root of
    # the units, for their sum to stand where independent errors would.
    systematic = target / sqrt(length(flagged))
    _laplace_parallel(laplace, 1:length(todo)) do i
        _continuation_refine_unit!(errors, laplace, todo[i], theta, Ls, opts, placed,
            target, systematic)
        return true
    end
    return errors
end

"""One unit of `_continuation_refine!`, as a function of its own so that the
task's locals are its own (a closure assigning a name its enclosing function
uses rebinds that function's)."""
function _continuation_refine_unit!(errors::Vector{Float64},
    laplace::CTSEMLaplaceObjective, U::Int, theta::Vector{Float64}, Ls, opts, placed,
    target::Float64, systematic::Float64=target)
    ws = _laplace_workspace!(laplace, Float64, length(theta))
    families = _continuation_steps(laplace, U, opts)
    nf = length(families)
    # Per family: the values and errors of its steps (NaN until evaluated),
    # the rules placed, and the next step to evaluate.
    values = [fill(NaN, length(f)) for f in families]
    ses = [fill(NaN, length(f)) for f in families]
    rules = [Vector{Any}(nothing, length(f)) for f in families]
    next = ones(Int, nf)
    done = falses(nf)
    # The screen's rule, where it is one of the steps, is evaluated already.
    current = placed.choices[U]
    for (f, steps) in enumerate(families), (s, c) in enumerate(steps)
        if c == current
            values[f][s] = placed.quadrature[U]
            ses[f][s] = placed.se[U]
            rules[f][s] = placed.rules[U]
        end
    end
    evaluate! = function (f, s)
        local r
        isnan(values[f][s]) || return isfinite(values[f][s])
        r = try
            _continuation_unit_rule(laplace, U, theta, Ls, ws,
                merge(opts, (choice=families[f][s],)))
        catch e
            _ctsem_must_propagate(e) && rethrow()
            nothing
        end
        if r === nothing || !isfinite(r.value)
            values[f][s] = -Inf
            return false
        end
        values[f][s] = r.value
        ses[f][s] = r.se
        rules[f][s] = r.rule
        return true
    end
    # Whether family `f` has a step under its target as it stands: the step
    # and its error, `(0, NaN)` where none.
    accepted = function (f)
        local steps, s, d, e
        steps = families[f]
        s = next[f] - 1
        s >= 1 || return (0, NaN)
        if steps[s].family === :lattice
            e = max(ses[f][s], _continuation_count_step(values[f], s))
            return e <= target ? (s, e) : (0, NaN)
        end
        if steps[s].family === :product
            # The tail two node steps imply, against the units' share of the
            # target: a product rule's error has one sign at every unit, so
            # the units' errors add (`_continuation_product_error`).
            e = _continuation_product_error(values[f], s)
            e === nothing && return (0, NaN)
            s - 1 >= 2 && e.coarse <= systematic && return (s - 1, e.coarse)
            return e.fine <= systematic ? (s, e.fine) : (0, NaN)
        end
        # A mixed rule's difference needs two steps; the first is only a check.
        s >= 2 || return (0, NaN)
        d = abs(values[f][s] - values[f][s - 1])
        d <= target || return (0, NaN)
        # The coarser rule stands, with the finer one's measure of it, unless
        # the coarser is only the check.
        return (s - 1 >= 2 ? s - 1 : s, d)
    end
    cost(f, s) = _rule_cost(laplace, U, families[f][s])
    # Start every family where the screen left it: the steps already
    # evaluated count as taken.
    for f in 1:nf
        while next[f] <= length(families[f]) && !isnan(values[f][next[f]]) &&
                (next[f] == 1 || families[f][next[f]].family === :lattice ||
                 !isnan(values[f][next[f] - 1]))
            next[f] += 1
        end
    end
    chosen = (f=0, s=0, err=NaN)
    while true
        for f in 1:nf
            done[f] && continue
            a = accepted(f)
            if a[1] != 0
                chosen = (f=f, s=a[1], err=a[2])
                break
            end
        end
        chosen.f == 0 || break
        # The family whose next step costs least: a product family's first
        # move takes its check and its first rule together.
        best = 0
        bestcost = Inf
        for f in 1:nf
            (done[f] || next[f] > length(families[f])) && continue
            c = cost(f, next[f])
            if families[f][next[f]].family !== :lattice && next[f] == 1 &&
                    length(families[f]) >= 2 && isnan(values[f][2])
                c += cost(f, 2)
            end
            if c < bestcost
                bestcost = c
                best = f
            end
        end
        best == 0 && break
        f = best
        if !evaluate!(f, next[f])
            done[f] = true
            continue
        end
        next[f] += 1
        if families[f][next[f] - 1].family !== :lattice && next[f] == 2 &&
                length(families[f]) >= 2
            if !evaluate!(f, 2)
                done[f] = true
                continue
            end
            next[f] = 3
        end
    end
    # Each family's finest rule and its error: a lattice's spread; a product
    # or mixed family's last difference, or where it is smaller than the one
    # before, the tail that ratio implies -- and whether it is converging.
    finest = [(s=0, err=NaN, converging=false) for _ in 1:nf]
    for f in 1:nf
        last = next[f] - 1
        (last >= 1 && isfinite(values[f][last])) || continue
        if families[f][last].family === :lattice
            cs = _continuation_count_step(values[f], last)
            finest[f] = (s=last, err=isfinite(cs) ? max(ses[f][last], cs) :
                ses[f][last], converging=false)
        elseif last >= 2 && isfinite(values[f][last - 1])
            dl = abs(values[f][last] - values[f][last - 1])
            dp = last >= 3 ? abs(values[f][last - 1] - values[f][last - 2]) : NaN
            if isfinite(dp) && dp > 0 && dl < dp
                ratio = dl / dp
                finest[f] = (s=last, err=dl * max(1.0, ratio / (1 - ratio)),
                    converging=true)
            else
                finest[f] = (s=last, err=dl, converging=false)
            end
        end
    end
    pick = function (ok)
        local best
        best = (f=0, s=0, err=NaN)
        for f in 1:nf
            (finest[f].s > 0 && isfinite(finest[f].err) && ok(f)) || continue
            (best.f == 0 || finest[f].err < best.err) &&
                (best = (f=f, s=finest[f].s, err=finest[f].err))
        end
        return best
    end
    if chosen.f == 0 && any(f -> finest[f].converging, 1:nf)
        # None reached the target, and the product rules converge: one basin
        # the nodes are closing on. The finest product rule with its tail, or
        # the lattice's finest where its spread is smaller.
        chosen = pick(f -> finest[f].converging ||
            families[f][finest[f].s].family === :lattice)
    end
    if chosen.f == 0 && opts.draws > 0
        r = try
            _continuation_unit_rule(laplace, U, theta, Ls, ws,
                merge(opts, (draw_units=(U,),)))
        catch e
            _ctsem_must_propagate(e) && rethrow()
            nothing
        end
        if r !== nothing
            placed.rules[U] = r.rule
            placed.quadrature[U] = r.value
            placed.se[U] = r.se
            placed.importance[U] = true
            errors[U] = r.se
            return nothing
        end
    end
    # Draws off or failed: each family's finest rule, the smaller error
    # standing.
    chosen.f == 0 && (chosen = pick(f -> true))
    if chosen.f == 0
        errors[U] = NaN
        return nothing
    end
    f, s = chosen.f, chosen.s
    placed.rules[U] = rules[f][s]
    placed.quadrature[U] = values[f][s]
    placed.se[U] = ses[f][s]
    placed.choices[U] = families[f][s]
    errors[U] = chosen.err
    return nothing
end

"""
    _continuation_resolved(gaps, se)

Each unit's gap beyond twice its rule's standard error, towards zero: what the
screen can tell of that unit alone from the draws' own error. The same gap
where the rule is a product rule, whose `se` is zero. What a set of units
leaves together is their signed sum against its combined error
(`_continuation_screen`), since gaps each below their own error sum to a
resolved total.
"""
_continuation_resolved(gaps::Vector{Float64}, se::Vector{Float64}) =
    [isfinite(g) ? sign(g) * max(abs(g) - 2 * se[U], 0.0) : g
     for (U, g) in enumerate(gaps)]

"""
    _continuation_screen(gaps, se)

What the screen can resolve of a set of units' gaps: the larger of the sum of
each unit's own resolved gap and the set's signed sum beyond twice its
combined error. Both are the plain absolute sum where every `se` is zero.
"""
function _continuation_screen(gaps::AbstractVector{Float64}, se::AbstractVector{Float64})
    own = 0.0; signed = 0.0; var = 0.0
    for (U, g) in enumerate(gaps)
        isfinite(g) || continue
        own += max(abs(g) - 2 * se[U], 0.0)
        signed += g
        var += se[U]^2
    end
    return max(own, abs(signed) - 2 * sqrt(var), 0.0)
end

"""
    _continuation_flags(gaps, tolerance; se=zeros(length(gaps)))

The units the continuation puts nodes on: largest `|gap|` first, the smallest
set whose removal leaves the rest carrying at most `tolerance` of what the
screen can resolve (`_continuation_screen`). That is the tolerance the screen
passes a whole fit on, so what the hybrid leaves to Laplace is what the screen
would have called exact. A unit whose rule could not be placed (`NaN`) is
never flagged. A negative `tolerance` flags every other unit, a zero gap
included: at `0` a unit whose gap is exactly zero, as on a linear Gaussian
model where Laplace is exact and the gaps are rounding, is left to Laplace, so
whether it is flagged would be a matter of the last digit.
"""
function _continuation_flags(gaps::Vector{Float64}, tolerance::Real;
    se::Vector{Float64}=zeros(length(gaps)))
    g = [isfinite(x) ? abs(x) : 0.0 for x in gaps]
    order = sortperm(g; rev=true)
    # What the unflagged units carry, as a sum of them rather than the total
    # less what has been flagged: subtracting the largest gaps from their total
    # can leave zero or less while a small one is still unflagged, which at
    # `tolerance = 0` dropped a unit on rounding alone.
    flagged = Int[]
    for (k, U) in enumerate(order)
        rest = order[k:end]
        _continuation_screen(gaps[rest], se[rest]) <= tolerance && break
        (g[U] > 0 || tolerance < 0) || break
        push!(flagged, U)
    end
    return sort!(flagged)
end

"""
    ctsem_laplace_continuation(laplace, values; nodes=5, tolerance=0.01,
        draws=0, draws_se=0.05, draws_max=16384, draws_floor=0.25,
        draws_prior=0.2, draws_df=4, draws_seed=20261021, force=nothing,
        screen=true, draws_all=false)

The hybrid objective with its nodes placed at `values`. Every unit is scored by
the rule there against its Laplace term, the units that differ beyond the
rule's own error are flagged (`_continuation_resolved`, `_continuation_flags`),
and those carry fixed nodes from then on while the rest keep their Laplace
term; the units in `force` are flagged whatever their gap. Every unit is
scored, on the rule its cost for the accuracy chooses (`_rule_default`,
`_continuation_refine!`); a flagged unit no rule converges on gets a frozen
importance rule, per block of its tree, when `draws > 0`, and with
`draws_all` every unit does. Each unit's known modes (`laplace.known`: the
modes a fit's jumps and searches saw, which R keeps on the fit's
specification) are components of its importance proposal
(`_continuation_known_points!`). `screen = false` scores only the units `force` names, the
rest keeping their Laplace terms unscored: what a fit climbing the hybrid
needs. See the module docstring for the rules and the constants.
"""
function ctsem_laplace_continuation(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, tolerance::Real=0.01,
    draws::Integer=0, draws_se::Real=0.05,
    draws_max::Integer=16384, draws_floor::Real=0.25, draws_prior::Real=0.2,
    draws_df::Real=4, draws_seed::Integer=20261021, force=nothing,
    screen::Bool=true, draws_all::Bool=false)
    nodes >= 3 || throw(ArgumentError("need at least three quadrature nodes"))
    draws_df > 0 || throw(ArgumentError("draws_df must be positive (Inf for Gaussian)"))
    draws_all && draws < 1 && throw(ArgumentError("draws_all needs draws"))
    draws >= 0 || throw(ArgumentError("draws must be at least zero"))
    draws_se > 0 || throw(ArgumentError("draws_se must be positive"))
    0 < draws_prior < 1 || throw(ArgumentError("draws_prior must be between 0 and 1"))
    draws_floor > 0 || throw(ArgumentError("draws_floor must be positive"))
    theta = _laplace_extend(laplace, collect(Float64, values))
    nunits = length(laplace.units.members)
    o = CTSEMLaplaceContinuation{typeof(laplace)}(laplace, laplace, collect(1:nunits),
        Int[], CTSEMFixedUnitRule[], theta, Int(nodes), Float64(tolerance),
        Bool(draws_all), Int(draws), Float64(draws_se),
        max(Int(draws_max), Int(draws)), Float64(draws_floor), Float64(draws_prior),
        Float64(draws_df), Int(draws_seed), Dict{Int,Vector{Vector{Float64}}}(),
        zeros(nunits), fill(NaN, nunits),
        fill(false, nunits), Bool(screen), Int[],
        fill(NaN, nunits), fill(NaN, nunits), fill(false, nunits), NaN,
        0, 0, 0, 0, 0, 0.0, 0.0, 0.0,
        false, Int[], CTSEMRuleChoice[_rule_default(laplace, U, nodes) for U in 1:nunits])
    _continuation_known_points!(o)
    _continuation_force!(o, force)
    _continuation_place!(o, theta)
    return o
end

"""Units `force` names, flagged from now on whatever their gap; out-of-range
entries are an error."""
function _continuation_force!(o::CTSEMLaplaceContinuation, force)
    force === nothing && return o
    nunits = length(o.laplace.units.members)
    units = Int[Int(U) for U in force]
    all(U -> 1 <= U <= nunits, units) || throw(ArgumentError(
        "force names a unit outside 1:$(nunits)"))
    o.flagged = sort!(union(o.flagged, units))
    o.forced = sort!(union(o.forced, units))
    return o
end

"""Each unit's known modes (`laplace.known`, those of the unit's length on this
objective) added to `o.points` as components for its importance proposal: the
modes a fit's jumps saw either side of where a unit's term jumped, which R
keeps on the fit's specification (`.ctJuliaKnownModes()`), as the placed
units' searches also start from them (`_placed_search!`). R once sent the
jump's points again beside them."""
function _continuation_known_points!(o::CTSEMLaplaceContinuation)
    laplace = o.laplace
    for U in eachindex(laplace.known)
        for x in laplace.known[U]
            (length(x) == laplace.units.dims[U] && all(isfinite, x)) || continue
            push!(get!(o.points, U, Vector{Vector{Float64}}()), copy(x))
        end
    end
    return o
end

"""The options every placement of `o`'s rules reads: the rules' constants,
which units it places, and the draws and node counts the units keep."""
_continuation_opts(o::CTSEMLaplaceContinuation) = (nodes=o.nodes,
    draws_all=o.draws_all, draws=o.draws, draws_se=o.draws_se,
    draws_max=o.draws_max, draws_floor=o.draws_floor, draws_prior=o.draws_prior,
    draws_df=o.draws_df, draws_seed=o.draws_seed, points=o.points,
    only=o.screen ? nothing : Set(o.flagged),
    draw_units=Set(o.drawn), unit_rules=o.unit_rules)

function _continuation_place!(o::CTSEMLaplaceContinuation, theta::Vector{Float64})
    started = time_ns()
    laplace = o.laplace
    opts = _continuation_opts(o)
    placed = _continuation_place(laplace, theta, opts)
    nunits = length(laplace.units.members)
    failed = [placed.rules[U] === nothing && !isempty(laplace.units.blocks[U])
              for U in 1:nunits]
    gaps = placed.quadrature .- placed.laplace
    flagged = sort!(union(o.flagged, _continuation_flags(gaps, o.tolerance;
        se=placed.se)))
    # A unit flagged before whose rule could not be placed here keeps its
    # Laplace term: an honest Laplace term beats a rule placed somewhere else.
    filter!(U -> placed.rules[U] !== nothing, flagged)
    errors = _continuation_refine!(laplace, theta, opts, placed, flagged)
    # A unit on draws took them in the refine's fallback, and keeps them
    # (`o.drawn`), unless every unit is on draws.
    o.draws_all || (o.drawn = sort!(union(o.drawn,
        [U for U in flagged if placed.importance[U]])))
    # And a unit keeps the rule the refine chose for it (`o.unit_rules`): a new
    # vector, since a saved placement holds the old one.
    o.unit_rules = placed.choices
    o.flagged = flagged
    o.rules = CTSEMFixedUnitRule[placed.rules[U] for U in flagged]
    # Not `setdiff(1:nunits, flagged)`: Base infers its working set only as
    # an `AbstractSet`, and Tables.jl, which the bridge's server loads, adds a
    # method that invalidates that instance and with it this function, the
    # coreset densities and DAIS above it -- 5 s of compile in a session's
    # first `ctFitUncertainty(fit, 'sample')` after a Laplace fit.
    isflagged = falses(nunits)
    isflagged[flagged] .= true
    o.rest_units = Int[U for U in 1:nunits if !isflagged[U]]
    o.rest = length(o.rest_units) == nunits ? laplace :
        _ctsem_subset_objective(laplace, o.rest_units)
    o.centre = copy(theta)
    o.centre_quadrature = placed.quadrature
    o.centre_laplace = placed.laplace
    o.centre_se = placed.se
    o.importance = placed.importance
    o.centre_error = errors
    o.rule_failed = failed
    prior = _ctsem_log_prior(laplace.objective, theta)
    o.centre_value = sum(_laplace_unit_weight(laplace, U) * placed.laplace[U]
        for U in o.rest_units; init=0.0) +
        sum(_laplace_unit_weight(laplace, U) * placed.quadrature[U]
        for U in flagged; init=0.0) + prior
    o.centre_converged = placed.converged
    o.recentres += 1
    o.seconds_placements += (time_ns() - started) / 1e9
    return o
end

"""The per-slot buffers one member's gradient at fixed `u` needs."""
_continuation_scratch(laplace::CTSEMLaplaceObjective, npar::Integer) = (
    shift=_laplace_scratch_vector!(laplace, Float64, npar, :cont_shift),
    grad=_laplace_scratch_vector!(laplace, Float64, npar, :cont_grad))

"""
    _continuation_member!(gh, laplace, U, m, theta, Ls, dL, positions, u, aws,
        scratch, want_gradient)

Member `m` of unit `U` at the fixed standardised effects `u`: its log
likelihood at `theta + L(theta) u`, and, when asked, its gradient in `theta`
added into `gh` -- one reverse sweep at the shifted parameters and the chain
through `L(theta)`, the partial `_laplace_floored_unit_gradient!` takes at the
mode, taken here at a node. Deviations the member's filter integrates are
integrated exactly (`_laplace_member_loglik`), their own variances read by the
filter directly.
"""
function _continuation_member!(gh, laplace::CTSEMLaplaceObjective, U::Integer,
    m::Integer, theta::Vector{Float64}, Ls, dL, positions, u::Vector{Float64},
    aws, scratch, want_gradient::Bool)
    spec = laplace.spec
    units = laplace.units
    i = units.members[U][m]
    offsets = units.offsets[U][m]
    shifted = _laplace_member_values!(scratch.shift, theta, spec, Ls, u, offsets)
    subject = laplace.objective.subject_objectives[i]
    want_gradient || return _laplace_member_loglik(aws, subject, shifted)
    grad = scratch.grad
    loglik = _laplace_subject_value_gradient!(grad, subject, aws, shifted, nothing,
        _HYBRID_FULL)
    isfinite(loglik) || return loglik
    all(isfinite, grad) || return NaN
    @inbounds for t in eachindex(grad)
        gh[t] += grad[t]
    end
    _laplace_chol_chain!(gh, grad, spec, dL, positions, u, offsets)
    return loglik
end

"""
    _continuation_rule_value(o, U, theta, Ls, dL, positions, rule, u, aws,
        scratch, want_gradient)

One block's log integral under its fixed rule at `theta`, with its gradient:
`log sum_j exp(w_j + h_j(theta)) + c`, whose derivative with the nodes held is
`sum_j pi_j dh_j/dtheta`, `pi` the normalised node weights -- the Fisher
identity. A rule's nodes are summed by `_continuation_node_sum`.
"""
function _continuation_rule_value(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, want_gradient::Bool)
    a = _continuation_node_sum(o, U, theta, Ls, dL, positions, rule, u, aws,
        scratch, want_gradient, 1:size(rule.points, 2))
    (a.failed || !isfinite(a.peak)) && return (value=NaN, gradient=a.gacc, ess=NaN)
    want_gradient && (a.gacc ./= a.acc)
    return (value=a.peak + log(a.acc) + rule.constant, gradient=a.gacc,
        ess=a.acc^2 / max(a.acc2, eps()))
end

"""
    _continuation_node_sum(o, U, theta, Ls, dL, positions, rule, u, aws,
        scratch, want_gradient, range)

A rule's nodes `range`, summed: their largest log term `peak`, the sum `acc`
of `exp(term - peak)` over them, and with the gradient the sum of each node's
score times its weight (`gacc`), unnormalised -- so a rule's nodes can be
summed in parts on different threads and the parts combined
(`_continuation_flagged`). One pass, the running maximum rescaling what has
been summed so far. A node's term is its
log weight plus, at a leaf, its members' log likelihood, and at an outer block
the log integrals of its children's rules as they were placed at that node
(`_continuation_rule_value`), the gradient theirs. A leaf node whose likelihood
is not finite carries no weight, as it does in `_quadrature_block`; a child
that cannot be scored fails the sum (`failed`), rather than silently removing
its node's mass. `peak` is `-Inf` where no node in the range scores.
"""
function _continuation_node_sum(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, want_gradient::Bool, range::UnitRange{Int})
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    columns = (block.offset + 1):(block.offset + block.size)
    leaf = isempty(rule.children)
    npar = length(theta)
    peak = -Inf
    acc = 0.0
    # The sum of the squared weights, for the nodes' effective number
    # `acc^2 / acc2` (`_continuation_flagged`), which says how far the rule has
    # drifted from the point it was placed at.
    acc2 = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    gh = want_gradient ? zeros(npar) : Float64[]
    for j in range
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        want_gradient && fill!(gh, 0.0)
        h = 0.0
        ok = true
        if leaf
            for m in block.members
                ll = _continuation_member!(gh, laplace, U, m, theta, Ls, dL,
                    positions, u, aws, scratch, want_gradient)
                if !isfinite(ll)
                    ok = false
                    break
                end
                h += ll
            end
        else
            for child in rule.children[j]
                r = _continuation_rule_value(o, U, theta, Ls, dL, positions, child,
                    u, aws, scratch, want_gradient)
                isfinite(r.value) || return (peak=NaN, acc=0.0, acc2=0.0,
                    gacc=gacc, failed=true)
                h += r.value
                want_gradient && (gh .+= r.gradient)
            end
        end
        ok || continue
        t = rule.logweights[j] + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                acc *= s
                acc2 *= s * s
                want_gradient && (gacc .*= s)
            end
            peak = t
        end
        w = exp(t - peak)
        acc += w
        acc2 += w * w
        want_gradient && (gacc .+= w .* gh)
    end
    return (peak=peak, acc=acc, acc2=acc2, gacc=gacc, failed=false)
end

"""
    _continuation_tasks(o, nslot)

The pieces `_continuation_flagged` hands the pool: each flagged unit whole, as
`(f, nothing)`, except that where there are fewer flagged units than slots a
unit whose rule has one root of at least two nodes and at least 64 member
evaluations -- a wide unit's draws, a nested unit's draws or nodes above its
leaves, a product rule over several effects -- is split into ranges of its
root's nodes, `(f, range)`, at least 32 evaluations a range and up to the
slots each unit's share leaves. Without that a fit whose one flagged unit has
4096 draws evaluated them on one thread while the rest of the pool waited.
"""
function _continuation_tasks(o::CTSEMLaplaceContinuation, nslot::Integer)
    nf = length(o.rules)
    share = nf < nslot ? cld(Int(nslot), nf) : 1
    tasks = Tuple{Int,Union{Nothing,UnitRange{Int}}}[]
    for (f, rule) in enumerate(o.rules)
        ranges = length(rule.roots) == 1 ? _continuation_ranges(
            size(rule.roots[1].points, 2), rule.evaluations, share) :
            UnitRange{Int}[]
        if length(ranges) > 1
            append!(tasks, ((f, range) for range in ranges))
        else
            push!(tasks, (f, nothing))
        end
    end
    return tasks
end

"""
    _continuation_task_runs(tasks)

The runs of `tasks` (`_continuation_tasks`) that belong to one flagged unit,
as `(f, run)` with `run` the tasks' indices: a unit's node ranges sit together
and in order, so both passes over the flagged units merge a run's pieces into
the unit's term (`_continuation_flagged`).
"""
function _continuation_task_runs(tasks::Vector{Tuple{Int,Union{Nothing,UnitRange{Int}}}})
    runs = Tuple{Int,UnitRange{Int}}[]
    at = 1
    while at <= length(tasks)
        f = tasks[at][1]
        stop = at
        while stop < length(tasks) && tasks[stop + 1][1] == f
            stop += 1
        end
        push!(runs, (f, at:stop))
        at = stop + 1
    end
    return runs
end

"""
    _continuation_ranges(n, evaluations, share)

A root's `n` nodes, costing `evaluations` member evaluations in all, as the
ranges `_continuation_tasks` splits them into for `share` slots: one range
unless `share > 1`, `n >= 2` and `evaluations >= 64`, otherwise at least 32
evaluations a range and at most `share` ranges.
"""
function _continuation_ranges(n::Integer, evaluations::Integer, share::Integer)
    (share > 1 && n >= 2 && evaluations >= 64) || return [1:Int(n)]
    ranges = UnitRange{Int}[]
    stop = 0
    for m in _continuation_split(n, min(Int(share), Int(n), cld(Int(evaluations), 32)))
        push!(ranges, (stop + 1):(stop + m))
        stop += m
    end
    return ranges
end

"""The flagged units' fixed-rule terms at `theta`, and their summed gradient,
in parallel over units and, where there are few of them, over a wide unit's
nodes (`_continuation_tasks`). Each piece's result is kept by the piece and
summed in order afterwards. With `scores` (`nunits x length(theta)`), each
flagged unit's row is its term's gradient times its sampling weight."""
function _continuation_flagged(o::CTSEMLaplaceContinuation, theta::Vector{Float64},
    want_gradient::Bool; scores::Union{Nothing,Matrix{Float64}}=nothing)
    laplace = o.laplace
    spec = laplace.spec
    nf = length(o.rules)
    npar = length(theta)
    values = fill(NaN, nf)
    nf == 0 && return (values=values, gradient=zeros(npar), ok=true)
    Ls = _laplace_popchols(theta, spec)
    dL = want_gradient ? _laplace_level_chol_derivatives(theta, spec) : nothing
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    _laplace_ensure_pool!(laplace)
    tasks = _continuation_tasks(o, _laplace_pool_width())
    ntask = length(tasks)
    peaks = fill(-Inf, ntask)
    accs = zeros(ntask)
    grads = [want_gradient ? zeros(npar) : Float64[] for _ in 1:ntask]
    ok = _laplace_parallel(laplace, 1:ntask) do k
        # `local` throughout: a name assigned here that is also a local of the
        # enclosing function would be one binding shared by every task.
        local f, range, rule, U, aws, scratch, u, total, r, failed, a
        f, range = tasks[k]
        rule = o.rules[f]
        U = rule.unit
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        failed = false
        try
            if range === nothing
                # The whole unit: its term as a peak with a sum of one, its
                # gradient already normalised.
                total = 0.0
                for root in rule.roots
                    r = _continuation_rule_value(o, U, theta, Ls, dL, positions,
                        root, u, aws, scratch, want_gradient)
                    if !isfinite(r.value)
                        failed = true
                        break
                    end
                    total += r.value
                    want_gradient && (grads[k] .+= r.gradient)
                end
                peaks[k] = total
                accs[k] = 1.0
            else
                a = _continuation_node_sum(o, U, theta, Ls, dL, positions,
                    rule.roots[1], u, aws, scratch, want_gradient, range)
                failed = a.failed
                peaks[k] = a.peak
                accs[k] = a.acc
                want_gradient && (grads[k] .= a.gacc)
            end
        catch err
            _ctsem_must_propagate(err) && rethrow()
            failed = true
        end
        return !failed
    end
    gradient = zeros(npar)
    ok || return (values=values, gradient=gradient, ok=false)
    g = want_gradient ? zeros(npar) : Float64[]
    for (f, run) in _continuation_task_runs(tasks)
        rule = o.rules[f]
        if tasks[first(run)][2] === nothing
            values[f] = peaks[first(run)]
            want_gradient && (g .= grads[first(run)])
        else
            P = maximum(@view peaks[run])
            isfinite(P) || return (values=values, gradient=gradient, ok=false)
            A = 0.0
            want_gradient && fill!(g, 0.0)
            for q in run
                e = exp(peaks[q] - P)
                A += accs[q] * e
                want_gradient && (g .+= e .* grads[q])
            end
            values[f] = P + log(A) + rule.roots[1].constant
            want_gradient && (g ./= A)
        end
        if want_gradient
            all(isfinite, g) || return (values=values, gradient=gradient, ok=false)
            gradient .+= _laplace_unit_weight(laplace, rule.unit) .* g
            scores === nothing ||
                (scores[rule.unit, :] .= _laplace_unit_weight(laplace, rule.unit) .* g)
        end
    end
    return (values=values, gradient=gradient, ok=true)
end

"""
    ctsem_laplace_continuation_evaluate(o, values; gradient=true)

The hybrid objective at `values`: the Laplace terms of the units left to it and
the prior, plus the fixed-rule terms of the flagged units. `converged` is false
when the Laplace part's inner solve did not converge or a flagged unit could
not be scored -- a finite value is not the objective there, as on the Laplace
route. `unit_loglik` is every unit's term in the fit's unit order, and
`subject_loglik` spreads a flagged unit's term evenly over its members, as
`ctsem_laplace_evaluate` does.
"""
function ctsem_laplace_continuation_evaluate(o::CTSEMLaplaceContinuation,
    values::AbstractVector; gradient::Bool=true)
    # Over the extended vector (`_laplace_extend`), the gradient returned at
    # the caller's length.
    n = length(values)
    theta = _laplace_extend(o.laplace, collect(Float64, values))
    return _laplace_free_result(_continuation_evaluate(o, theta, gradient), n)
end

# Copies of the vectors, so that neither a caller nor a cache can change what
# the other holds (`ctsem_laplace_placed_evaluate`'s).
_continuation_copy(r) = (value=r.value,
    gradient=r.gradient === nothing ? nothing : copy(r.gradient),
    unit_loglik=copy(r.unit_loglik), subject_loglik=copy(r.subject_loglik),
    converged=r.converged)

# With `scores` (`nunits x length(theta)`), every unit's score, as the Laplace
# route gives them (`ctsem_laplace_evaluate`): the rest's from it, the flagged
# units' from their rules.
function _continuation_evaluate(o::CTSEMLaplaceContinuation, theta::Vector{Float64},
    gradient::Bool; scores::Union{Nothing,Matrix{Float64}}=nothing)
    started = time_ns()
    rs = scores === nothing ? nothing : zeros(length(o.rest_units), length(theta))
    rest = ctsem_laplace_evaluate(o.rest, theta; gradient=gradient, scores=rs)
    flagged = _continuation_flagged(o, theta, gradient; scores=scores)
    if scores !== nothing
        for (position, U) in enumerate(o.rest_units)
            scores[U, :] .= view(rs, position, :)
        end
    end
    gradient ? (o.gradient_calls += 1) : (o.value_calls += 1)
    evaluations = sum(r.evaluations for r in o.rules; init=0)
    gradient ? (o.member_sweeps += evaluations) : (o.member_values += evaluations)
    value = flagged.ok ? rest.value + sum(_laplace_unit_weight(o.laplace, U) *
        flagged.values[f] for (f, U) in enumerate(o.flagged); init=0.0) : NaN
    grad = gradient ? (flagged.ok ? rest.gradient .+ flagged.gradient :
        fill(NaN, length(theta))) : nothing
    nunits = length(o.laplace.units.members)
    unit_loglik = zeros(nunits)
    for (position, U) in enumerate(o.rest_units)
        unit_loglik[U] = rest.unit_loglik[position]
    end
    for (f, U) in enumerate(o.flagged)
        unit_loglik[U] = flagged.values[f]
    end
    subject_loglik = copy(rest.subject_loglik)
    for (f, U) in enumerate(o.flagged)
        members = o.laplace.units.members[U]
        for i in members
            subject_loglik[i] = flagged.values[f] / length(members)
        end
    end
    elapsed = (time_ns() - started) / 1e9
    gradient ? (o.seconds_gradients += elapsed) : (o.seconds_values += elapsed)
    return (value=value, gradient=grad, unit_loglik=unit_loglik,
        subject_loglik=subject_loglik, converged=rest.converged && flagged.ok)
end

# A Laplace objective keeps the pairs its last run ended with where the run
# asked to (`lbfgs_pairs = :objective`, `ctsem_optimize`), in a weak table, so
# an objective that is let go takes them with it. The refined fit's placed
# objective starts from them where it holds none of its own
# (`_ctsem_lbfgs_held(::CTSEMLaplacePlaced)`).
const _CTSEM_LBFGS_HELD = WeakKeyDict{Any,CTSEMLBFGSMemory}()
_ctsem_lbfgs_held(o::CTSEMLaplaceObjective) = get(_CTSEM_LBFGS_HELD, o, nothing)
_ctsem_lbfgs_hold!(o::CTSEMLaplaceObjective, memory::CTSEMLBFGSMemory) =
    (_CTSEM_LBFGS_HELD[o] = memory; nothing)

"""
    _continuation_error(o)

The hybrid's own error against the exact marginal at its placement, in nats:
the frozen draws' and the lattice rules' standard errors in quadrature
(`draws`), each independent of the others', the product and mixed rules' node
errors summed (`nodes`), what the screen leaves on Laplace (`residual`), and
their sum (`total`). What `ctsem_laplace_continuation_info` reports.
"""
function _continuation_error(o::CTSEMLaplaceContinuation)
    gaps = o.centre_quadrature .- o.centre_laplace
    residual = _continuation_screen(gaps[o.rest_units], o.centre_se[o.rest_units])
    # Draws' and lattices' spreads are random and add in quadrature; product
    # rules' node errors, and a lattice's count step where it exceeds its
    # spread (`_continuation_count_step`), share a sign and add.
    random(U) = o.importance[U] || (o.unit_rules[U].family === :lattice &&
        !(o.centre_error[U] > o.centre_se[U]))
    draws = sqrt(sum(abs2(o.importance[U] ? o.centre_se[U] : o.centre_error[U])
        for U in o.flagged if random(U) && isfinite(o.centre_error[U]); init=0.0))
    nodes = sum(o.centre_error[U] for U in o.flagged
        if !random(U) && isfinite(o.centre_error[U]); init=0.0)
    return (draws=draws, nodes=nodes, residual=residual,
        total=draws + nodes + residual)
end

"""
    ctsem_laplace_continuation_info(o)

What the continuation's rule looks like at its current centre, for the R side.
Never an empty vector, since one deadlocks the bridge: `flagged` is `[0]` when
nothing is flagged, and `nflagged` says so.
"""
function ctsem_laplace_continuation_info(o::CTSEMLaplaceContinuation)
    nunits = length(o.laplace.units.members)
    gaps = o.centre_quadrature .- o.centre_laplace
    # The screen reads the gaps it can resolve (`_continuation_screen`); the
    # signed `gap` reported is the raw difference, and each unit's own
    # resolved gap sets the bar.
    finite = [isfinite(g) ? abs(g) : 0.0 for g in _continuation_resolved(gaps, o.centre_se)]
    flaggedgaps = [finite[U] for U in o.flagged]
    prior = _ctsem_log_prior(o.laplace.objective, o.centre)
    quadrature_units = [isfinite(o.centre_quadrature[U]) ? o.centre_quadrature[U] :
        o.centre_laplace[U] for U in 1:nunits]
    members = o.laplace.units.members
    nsubjects = length(o.laplace.objective.subject_objectives)
    subjects = zeros(nsubjects)
    for U in 1:nunits, i in members[U]
        subjects[i] = quadrature_units[U] / length(members[U])
    end
    # The hybrid's own error against the exact marginal, in nats: the frozen
    # rules' standard errors in quadrature, the product rules' node errors,
    # and what the screen leaves on Laplace (`residual`).
    err = _continuation_error(o)
    residual = err.residual
    error_draws = err.draws
    error_nodes = err.nodes
    return (flagged=isempty(o.flagged) ? [0] : copy(o.flagged),
        nflagged=length(o.flagged), nunits=nunits,
        screen=_continuation_screen(gaps, o.centre_se),
        gap=sum(gaps[isfinite.(gaps)]; init=0.0),
        bar=isempty(flaggedgaps) ? 0.0 : minimum(flaggedgaps),
        residual=residual,
        centre_error=copy(o.centre_error), error_draws=error_draws,
        error_nodes=error_nodes, error_residual=residual,
        error=error_draws + error_nodes + residual,
        centre=copy(o.centre), value=o.centre_value,
        quadrature=sum(quadrature_units) + prior, prior=prior,
        quadrature_units=quadrature_units, quadrature_subjects=subjects,
        laplace_units=copy(o.centre_laplace),
        laplace=sum(o.centre_laplace) + prior,
        rule_failures=count(o.rule_failed),
        nlattice=count(U -> !o.importance[U] &&
            o.unit_rules[U].family !== :product, o.flagged),
        importance=[o.importance[U] for U in 1:nunits], nimportance=count(o.importance),
        centre_se=copy(o.centre_se),
        importance_draws=sum(sum(_continuation_draws(root) for root in r.roots;
            init=0) for r in o.rules if o.importance[r.unit]; init=0),
        draws=o.draws, draws_se=o.draws_se,
        evaluations_per_value=sum(r.evaluations for r in o.rules; init=0),
        value_calls=o.value_calls, gradient_calls=o.gradient_calls,
        member_values=o.member_values, member_sweeps=o.member_sweeps,
        recentres=o.recentres,
        seconds_values=o.seconds_values, seconds_gradients=o.seconds_gradients,
        seconds_placements=o.seconds_placements,
        nodes=o.nodes, tolerance=o.tolerance)
end

"""
    _ctsem_jump_units(objective, jump)

The units whose terms jumped along the step the line search saw jump
(`CTSEMJump`), with their changes and where each one's inner solve stopped
either side, as `(units, changes, point_units, point_values)`; `nothing` on a
route with no units to name. On the Laplace route each unit's term is taken at
`x`, `x + alpha_lo s` and `x + alpha_hi s`, and its score at `x`: a unit whose
excess over its own linear change does not shrink with the square of the step
between the two (`_ctsem_jump_excess`), and carries a tenth of the jump, moved
with where its inner solve stopped rather than with the parameters
(`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`). Compared
whole, as they were while the two steps stood a hundredfold apart, a unit's
changes at steps four times apart keep a quarter of each other wherever its
term is linear in the step. `changes` are each named unit's excess at
`alpha_lo`.
Its modes at the three points are what a frozen importance rule for it starts
from (`ctsem_laplace_continuation`'s `point_units`): at one point the solve
returns one of them, and the rule needs both sides.
"""
_ctsem_jump_units(::Any, jump) = nothing
_ctsem_jump_capable(::CTSEMLaplaceObjective) = true

function _ctsem_jump_units(laplace::CTSEMLaplaceObjective, jump::CTSEMJump)
    return _laplace_jump_units(laplace, collect(1:length(laplace.units.members)),
        jump)
end

function _laplace_jump_units(laplace::CTSEMLaplaceObjective, numbering::Vector{Int},
    jump::CTSEMJump)
    nunits = length(numbering)
    # Each unit's score at `x`, weighted as its term is below, for its linear
    # change along `s`.
    S = zeros(nunits, length(jump.x))
    r0 = ctsem_laplace_evaluate(laplace, jump.x; gradient=true, scores=S)
    slope = all(isfinite, S) ? S * jump.s : fill(NaN, nunits)
    terms = Vector{Vector{Float64}}([copy(r0.unit_loglik)])
    modes = Vector{Vector{Vector{Float64}}}([[copy(m) for m in laplace.modes]])
    for th in (jump.x .+ jump.alpha_lo .* jump.s, jump.x .+ jump.alpha_hi .* jump.s)
        r = ctsem_laplace_evaluate(laplace, th; gradient=false)
        push!(terms, copy(r.unit_loglik))
        push!(modes, [copy(m) for m in laplace.modes])
    end
    units = Int[]
    changes = Float64[]
    point_units = Int[]
    columns = Vector{Vector{Float64}}()
    for p in 1:nunits
        w = _laplace_unit_weight(laplace, p)
        lo = w * (terms[2][p] - terms[1][p]) - jump.alpha_lo * slope[p]
        hi = w * (terms[3][p] - terms[1][p]) - jump.alpha_hi * slope[p]
        (isfinite(lo) && isfinite(hi)) || continue
        abs(lo) >= jump.excess / 10 || continue
        _ctsem_jump_excess(jump.alpha_lo, abs(lo), jump.alpha_hi, abs(hi)) ||
            continue
        push!(units, numbering[p])
        push!(changes, lo)
        for k in 1:3
            if !any(c -> c == modes[k][p], (columns[j] for j in eachindex(columns)
                    if point_units[j] == numbering[p]))
                push!(point_units, numbering[p])
                push!(columns, modes[k][p])
            end
        end
    end
    width = maximum(length, columns; init=1)
    values = fill(NaN, width, max(length(columns), 1))
    for (j, c) in enumerate(columns)
        values[1:length(c), j] = c
    end
    return (units=units, changes=changes,
        point_units=isempty(point_units) ? [0] : point_units, point_values=values)
end

"""
    ctsem_laplace_continuation_se(o, values)

Each flagged unit's importance rule's standard error at `values`, its draws
held where they were placed, read as at placement: the spread of each root's
replicate lattices' log integrals over `sqrt` of their number
(`_CONTINUATION_REPLICATES`), the roots' combined as independent errors; `0`
for a unit on a product rule. Each replicate's draws keep the log weights its
placement gave them, which do not move with the parameters, so at the point
of placement this is the standard error placed. A root above other blocks is weighted by its children's
integrals as they were placed at each of its draws, so its error counts
theirs (`_continuation_importance_rule`). Where fewer than two replicates
score, the independent-draws formula `sqrt(1/ESS - 1/n)` stands in, as at
placement. The draws are the proposal's of the point they were placed at, and
as the parameters move away from it their weights spread; the refined fit
places them again only once their effective size has collapsed
(`ctsem_laplace_placed_optimize`). In the order of
`ctsem_laplace_continuation_info`'s `flagged`; `[0.0]` when none is flagged.
In parallel over units and, as a value is, over a root's draws
(`_continuation_ranges`).
"""
function ctsem_laplace_continuation_se(o::CTSEMLaplaceContinuation,
    values::AbstractVector)
    nf = length(o.rules)
    nf == 0 && return [0.0]
    theta = _laplace_extend(o.laplace, collect(Float64, values))
    laplace = o.laplace
    spec = laplace.spec
    Ls = _laplace_popchols(theta, spec)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    npar = length(theta)
    _laplace_ensure_pool!(laplace)
    share = nf < _laplace_pool_width() ? cld(_laplace_pool_width(), nf) : 1
    # Every root's draws' log integrands beneath (its members' log likelihood,
    # or its children's log integrals), filled by pieces of ranges: (flagged
    # position, root, range).
    roots = Tuple{Int,Int}[]
    terms = Vector{Vector{Float64}}()
    tasks = Tuple{Int,Int,UnitRange{Int}}[]
    for (f, rule) in enumerate(o.rules)
        o.importance[rule.unit] || continue
        for (r, root) in enumerate(rule.roots)
            push!(roots, (f, r))
            push!(terms, fill(-Inf, size(root.points, 2)))
            for range in _continuation_ranges(size(root.points, 2),
                    _continuation_evaluations(root, laplace.units.blocks[rule.unit]),
                    length(rule.roots) == 1 ? share : 1)
                push!(tasks, (length(roots), r, range))
            end
        end
    end
    failed = fill(false, length(roots))
    _laplace_parallel(laplace, 1:length(tasks)) do q
        local i, r, range, rule, U, aws, scratch, u, ok
        i, r, range = tasks[q]
        rule = o.rules[roots[i][1]]
        U = rule.unit
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        ok = try
            _continuation_node_terms!(terms[i], o, U, theta, Ls, positions,
                rule.roots[r], u, aws, scratch, range)
        catch err
            _ctsem_must_propagate(err) && rethrow()
            false
        end
        ok || (failed[i] = true)
        return true
    end
    variance = zeros(nf)
    for (i, (f, r)) in enumerate(roots)
        root = o.rules[f].roots[r]
        h = terms[i]
        n = length(h)
        t = root.logweights .+ h
        peak = maximum(t; init=-Inf)
        if failed[i] || !isfinite(peak)
            variance[f] = Inf
            continue
        end
        reps = Float64[]
        if length(root.replicate) == n
            for rep in 1:_CONTINUATION_REPLICATES
                sel = findall(==(rep), root.replicate)
                isempty(sel) && continue
                tr = root.replicate_logweights[sel] .+ h[sel]
                pk = maximum(tr)
                isfinite(pk) && push!(reps, pk + log(sum(exp.(tr .- pk))))
            end
        end
        if length(reps) >= 2
            mu = sum(reps) / length(reps)
            variance[f] += sum(abs2, reps .- mu) / (length(reps) - 1) / length(reps)
        else
            w = exp.(t .- peak)
            variance[f] += max(sum(abs2, w) / sum(w)^2 - 1 / n, 0.0)
        end
    end
    return sqrt.(variance)
end

"""
    _continuation_node_terms!(terms, o, U, theta, Ls, positions, rule, u, aws,
        scratch, range)

Each of a rule's nodes `range`'s log integrand at `theta` into `terms`, its
weight left out: its members' log likelihood at a leaf (`-Inf` where that is
not finite, as `_continuation_node_sum` skips it), or its children's log
integrals at an outer block. `false` where a child cannot be scored.
"""
function _continuation_node_terms!(terms::Vector{Float64}, o::CTSEMLaplaceContinuation,
    U::Integer, theta::Vector{Float64}, Ls, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, range::UnitRange{Int})
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    columns = (block.offset + 1):(block.offset + block.size)
    leaf = isempty(rule.children)
    for j in range
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        h = 0.0
        ok = true
        if leaf
            for m in block.members
                ll = _continuation_member!(Float64[], laplace, U, m, theta, Ls, nothing,
                    positions, u, aws, scratch, false)
                if !isfinite(ll)
                    ok = false
                    break
                end
                h += ll
            end
        else
            for child in rule.children[j]
                r = _continuation_rule_value(o, U, theta, Ls, nothing, positions, child,
                    u, aws, scratch, false)
                isfinite(r.value) || return false
                h += r.value
            end
        end
        terms[j] = ok ? h : -Inf
    end
    return true
end

"""
    ctsem_laplace_check(laplace, values; draws=128, seeds=1:4, draws_df=4,
        unit_scores=true, expansion=true, expansion_gate=1.0)

The Laplace term checked at `values` on every unit, whatever its width
(review/PROPOSAL-laplace-check-every-unit-2026-10-27.md): at each seed every
unit on `draws` frozen importance draws (the continuation with `draws_all`
and every unit flagged), and the corrected objective's value and gradient
against the Laplace objective's. The seeds are independent replicates, so
their spread is the check's own noise: at ten effects a unit with a weakly
identified covariance one seed's predicted shift was mostly noise.

`laplace_unit` is each unit's Laplace term; `check_unit` and `check_se` (units
x seeds) its corrected term and that term's replicate standard error;
`gradient` (free parameters x seeds) the corrected minus the Laplace gradient;
with `unit_scores`, `unit_gradient` (units x free parameters) each unit's
corrected minus Laplace score, averaged over the seeds that scored; and
A seed whose rules could not score every unit leaves its columns NaN. With `expansion`, each unit's
Laplace term with its O(1/n) correction, a bound on that one's error and the
curvature's smallest eigenvalue (`ctsem_laplace_expansion`; `expansion_term`,
`expansion_correction`, `expansion_bound`, `expansion_mineig`; the
correction only where the eigenvalue reaches `expansion_gate`), from which the
caller takes the gap of every unit that passes its gate.
"""
function ctsem_laplace_check(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    draws::Integer=128, seeds=1:4, draws_df::Real=4, unit_scores::Bool=true,
    expansion::Bool=true, expansion_gate::Real=1.0)
    draws >= 1 || throw(ArgumentError("draws must be at least one"))
    started = time_ns()
    n = length(values)
    x = collect(Float64, values)
    theta = _laplace_extend(laplace, x)
    nunits = length(laplace.units.members)
    seeds = collect(Int, seeds)
    ns = length(seeds)
    lscores = unit_scores ? zeros(nunits, n) : nothing
    lap = ctsem_laplace_evaluate(laplace, x; gradient=true, scores=lscores)
    gradient = fill(NaN, n, ns)
    check_unit = fill(NaN, nunits, ns)
    check_se = fill(NaN, nunits, ns)
    unit_gradient = unit_scores ? zeros(nunits, n) : nothing
    cscores = unit_scores ? zeros(nunits, length(theta)) : nothing
    scored = 0
    for (j, seed) in enumerate(seeds)
        o = ctsem_laplace_continuation(laplace, x; draws=draws, draws_max=draws,
            draws_se=1.0, draws_all=true, tolerance=-1.0, draws_seed=seed,
            draws_df=draws_df)
        unit_scores && fill!(cscores, 0.0)
        r = _continuation_evaluate(o, theta, true; scores=cscores)
        (r.converged && all(isfinite, r.gradient)) || continue
        gradient[:, j] .= view(r.gradient, 1:n) .- lap.gradient
        check_unit[:, j] .= r.unit_loglik
        check_se[:, j] .= o.centre_se
        unit_scores && (unit_gradient .+= view(cscores, :, 1:n) .- lscores)
        scored += 1
    end
    unit_scores && scored > 0 && (unit_gradient ./= scored)
    ex = expansion ? ctsem_laplace_expansion(laplace, x; gate=expansion_gate) :
        (term=Float64[], correction=Float64[], bound=Float64[], mineig=Float64[])
    return (laplace_unit=collect(Float64, lap.unit_loglik), laplace_value=Float64(lap.value),
        check_unit=check_unit, check_se=check_se, gradient=gradient,
        unit_gradient=unit_gradient, scored=scored,
        expansion_term=ex.term, expansion_correction=ex.correction,
        expansion_bound=ex.bound, expansion_mineig=ex.mineig,
        seconds=(time_ns() - started) / 1e9)
end
