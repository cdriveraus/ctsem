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
nodes move, so no round loses a unit it was scored on.

Nodes are held fixed in the *standardised* effects `u` -- the non-centred
coordinates the engine already works in, where a member's parameters are
`theta + L(theta) u` and the prior is `N(0, I)` whatever `theta` is. That
choice keeps the prior's density out of the gradient entirely, handles a
reduced-rank level (whose `L` has zero columns and no density in the raw effect
space) without a special case, and reuses the engine's fixed-`u` score as it
stands. Its cost is that a node fixed in `u` moves in the raw effect space as a
population scale changes; the rounds below re-place the nodes before that
matters, and a round is only kept when the quadrature objective with its nodes
re-placed at the new point is higher.

# The rule

Per block of the unit's tree, the rule `_quadrature_block` uses: a leaf
re-solves its conditional mode at each node of the blocks above it (held fixed
afterwards: fixed conditional maps), an outer block keeps the joint mode with
the eliminated precision's scale. A block no wider than `product_maxdim` gets
the product Gauss-Hermite rule, so at the centre the fixed rule and
`ctsem_laplace_quadrature` agree to rounding. A wider one gets `nodes` points
along each of its softest directions -- precision eigenvalues below `soft_tau`,
at most `soft_maxdirs` of them and at least one -- at the prior-clipped scale,
and in the stiff complement the same product rule (`_continuation_stiff_rule`,
which says why nothing cheaper will do), centred at a leaf on the complement's
conditional mode at that soft node and scaled by its conditional precision
there, at an outer block on the joint mode, scaled by the eliminated precision.
So a soft rule costs `nodes^k` like a product rule; what it adds is the
prior-clipped scale on its soft directions and a complement that follows its
conditional mode along them.

# Where the constants come from

`nodes = 5` and `tolerance = 0.01` are the screen's (`.ctLaplaceCorrectDefaults`
on the R side). `product_maxdim = 2` is cost: a product rule costs `nodes^k`
reverse sweeps per member per gradient, 25 at `k = 2` and 125 at `k = 3`, where
a Laplace gradient is `k + 1` sweeps and a curvature. `soft_tau = 3.5` is the
softcut the gated-gaps job's exact reference settled on
(`CT-SEM/review/LAPLACE-gated-gaps-2026-09-24.md` 1d): a direction whose
precision is below it can be non-Gaussian enough to matter. `soft_maxdirs = 2`
bounds how many directions the complement's conditional mode is followed along.
None of these was tuned on a fit.

`maxdim = 5` is cost. A unit's rule costs `nodes^d` member evaluations per
value and reverse sweeps per gradient, `d` the widest root-to-leaf path of its
block tree, so each effect multiplies it by 5. Measured on dev1 at 8 threads,
150 subjects of two latents (dev/lapcontinue/highdim.R): at `d = 4` the
correction took 92 s beside a 343 s Laplace fit (the step correction 57 s), at
`d = 5` 727 s beside 184 s (the step 579 s), two thirds of it the Hessian. At
`d = 6` that is an hour, so the default stops at 5. Since the Hessian became
exact (`ctsem_laplace_continuation_hessian`) the `d = 5` correction takes
132 s there beside a 170-second fit, 30 s of it the Hessian; the cap has not
been revisited. A unit wider
than `maxdim` gets a frozen importance rule instead when `draws > 0` (below);
with `draws = 0` it is not scored and keeps its Laplace term (`wide` in the
info): the step correction, the alternative, pays the same `nodes^d` per
quadrature value, so falling back to it would save nothing.

# Units too wide for a product rule: frozen importance draws

A unit wider than `maxdim` gets draws in place of nodes
(`_continuation_importance_rule`) -- every such unit, or with `draws_scope =
:forced` only those the caller names, the rest keeping their Laplace terms:
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
the stopping point, so its term is smooth.

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
    ctsem_laplace_continuation_evaluate, ctsem_laplace_continuation_recentre!,
    ctsem_laplace_continuation_se,
    ctsem_laplace_continuation_revert!, ctsem_laplace_continuation_optimize,
    ctsem_laplace_continuation_hessian

"""
One block's fixed rule, for one configuration of the blocks above it.

`points[:, j]` is the block's own coordinates at node `j`, in the unit's
standardised `u`; `logweights[j]` carries everything about the node that does
not depend on the parameters (the rule's weight, the Gaussian factor of any
stiff complement, and the block's own `-z'z/2`), and `constant` what the block
adds to its log integral whatever the node. A leaf's value at a node is its
members' log likelihood there; an outer block's is the sum of its children's
rules, and `children[j]` holds those rules as they were placed with this block
at node `j`. `soft` is the number of soft directions a soft rule used, zero for
a product rule. `replicate` is an importance rule's replicate lattice for each
point (`_CONTINUATION_REPLICATES`) and `replicate_logweights` each point's log
weight within its replicate's own estimate, both fixed with the draws, which
its standard error is read from away from where it was placed
(`ctsem_laplace_continuation_se`); empty for a product rule.
"""
struct CTSEMFixedRule
    block::Int
    points::Matrix{Float64}
    logweights::Vector{Float64}
    constant::Float64
    children::Vector{Vector{CTSEMFixedRule}}
    soft::Int
    replicate::Vector{Int}
    replicate_logweights::Vector{Float64}
end

CTSEMFixedRule(block, points, logweights, constant, children, soft) =
    CTSEMFixedRule(block, points, logweights, constant, children, soft, Int[],
        Float64[])

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

_continuation_soft_blocks(rule::CTSEMFixedRule) = (rule.soft > 0 ? 1 : 0) +
    sum(sum(_continuation_soft_blocks(c) for c in kids; init=0)
        for kids in rule.children; init=0)

"""The points a rule holds, at every block of its tree: an importance rule's
draws in all."""
_continuation_draws(rule::CTSEMFixedRule) = size(rule.points, 2) +
    sum(sum(_continuation_draws(c) for c in kids; init=0)
        for kids in rule.children; init=0)

"""
    CTSEMLaplaceContinuation

The hybrid objective the continuation climbs, with the state it was placed
from. See the module docstring. `previous` holds the placement before the last
`ctsem_laplace_continuation_recentre!`, so a round the caller rejects can be
undone without placing the rule again.
"""
mutable struct CTSEMLaplaceContinuation{L} <: CTSEMOptimisable
    laplace::L
    rest::L
    rest_units::Vector{Int}
    flagged::Vector{Int}
    rules::Vector{CTSEMFixedUnitRule}
    centre::Vector{Float64}
    nodes::Int
    tolerance::Float64
    product_maxdim::Int
    soft_tau::Float64
    soft_maxdirs::Int
    maxdim::Int
    # Units wider than `maxdim` that no importance rule can take (`draws = 0`,
    # or not named under `draws_scope = :forced`): never scored, never flagged.
    wide::Vector{Bool}
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
    # Which units too wide for a product rule draws may take: every one of one
    # block (`:all`), or only those `force` named (`:forced`), and those units.
    draws_scope::Symbol
    forced::Vector{Int}
    # What the rule scored every unit at the centre, and what Laplace did.
    centre_quadrature::Vector{Float64}
    centre_laplace::Vector{Float64}
    rule_failed::Vector{Bool}
    centre_value::Float64
    previous::Any
    # Accounting over the object's life.
    value_calls::Int
    gradient_calls::Int
    member_values::Int
    member_sweeps::Int
    recentres::Int
    refused::Int
    # Wall seconds in each kind of call, so a correction's cost can be broken
    # down without a profiler: values, gradients, placements (the screen is
    # the first), Hessians.
    seconds_values::Float64
    seconds_gradients::Float64
    seconds_placements::Float64
    seconds_hessian::Float64
    # The last evaluation with a gradient and the last without, each at its
    # point, for the placement it was made under: a round's first evaluation
    # is at the point the previous round's stationarity check just evaluated,
    # and the uncertainty stage asks for the gradient at the estimate again
    # after the certification's probe has taken values elsewhere. One slot for
    # both let the probe's values evict that gradient, which was then paid
    # twice: a whole gradient per converged correction, one of the seven
    # ord4's took. Both cleared by every placement.
    cache::Any
    value_cache::Any
    # The curvature the last `ctsem_optimize` run on this objective ended with,
    # for the next to start from (`lbfgs_pairs = :objective`): a rescued fit's
    # stages, one after each placement. Kept across placements, which change
    # the objective but not its coordinates.
    lbfgs::Union{Nothing,CTSEMLBFGSMemory}
    # Per flagged rule, its effective number of nodes at the last evaluation
    # and at the first after its placement (`_continuation_record_ess!`),
    # which `_continuation_drifted` compares.
    ess::Vector{Float64}
    ess_placed::Vector{Float64}
end


"""Gauss-Hermite grids a rule of these options can ask for, filled before any
threaded region so that `_gh_grid` is a read there (see
`_quadrature_warm_caches`): a product rule's, for a block no wider than
`product_maxdim`, and the soft directions', at most `soft_maxdirs` of them. A
wider block's stiff complement is built where it is placed
(`_continuation_stiff_rule`), so its size is not warmed: warming every block
size built a 5^13-point grid for units of thirteen effects, which `maxdim`
then left out."""
function _continuation_warm_caches(laplace::CTSEMLaplaceObjective, nodes::Integer,
    soft_maxdirs::Integer, product_maxdim::Integer)
    _quadrature_warm_caches(laplace, nodes,
        1:max(1, Int(soft_maxdirs), Int(product_maxdim)))
    return nothing
end

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
    _continuation_complement(fns, base, Vh, laplace)

The stiff complement's conditional mode at one soft node, and its conditional
precision there with the eigenvalues clipped at the prior's: their log
determinant, and the eigenvectors and clipped eigenvalues that scale the
complement's nodes.

Newton on `h` for `G(h) = g(base + Vh h)` from zero, with step halving and the
precision's eigenvalues floored for the step, since the complement can go
convex away from the mode. Converged rather than one step: nothing
differentiates through this point -- the node is fixed once placed -- so a
converged mode costs only its accuracy, and a fixed number of steps is what
made the gated-gaps job's first trapezoid rule chaotic in theta.
"""
function _continuation_complement(fns, base::Vector{Float64}, Vh::Matrix{Float64},
    laplace::CTSEMLaplaceObjective)
    nh = size(Vh, 2)
    nh == 0 && return (h=Float64[], logdet=0.0, vectors=zeros(0, 0),
        values=Float64[])
    h = zeros(nh)
    z = copy(base)
    current = fns.value_gradient(z)
    isfinite(current.value) ||
        return (h=h, logdet=NaN, vectors=zeros(nh, nh), values=fill(NaN, nh))
    for _ in 1:50
        gh = transpose(Vh) * current.gradient
        maximum(abs, gh) < _laplace_inner_tolerance(laplace, current.value) && break
        E = _ctsem_symeig(transpose(Vh) * fns.precision(z) * Vh)
        lowest = 1e-8 * max(maximum(abs, E.values), 1.0)
        step = E.vectors * ((transpose(E.vectors) * gh) ./ max.(E.values, lowest))
        # A step at roundoff is a mode however large the gradient's roundoff
        # is; see `_laplace_newton_unit_mode`.
        all(>=(lowest), E.values) && maximum(abs, step) < laplace.inner_tol && break
        accepted = false
        scale = 1.0
        for _ in 1:20
            hn = h .+ scale .* step
            zn = base .+ Vh * hn
            trial = fns.value_gradient(zn)
            if isfinite(trial.value) && trial.value >= current.value - 1e-12
                h = hn; z = zn; current = trial; accepted = true
                break
            end
            scale /= 2
        end
        accepted || break
    end
    E = _ctsem_symeig(transpose(Vh) * fns.precision(z) * Vh)
    lam = max.(E.values, 1.0)
    return (h=h, logdet=sum(log, lam; init=0.0), vectors=Matrix(E.vectors),
        values=lam)
end

"""
    _continuation_stiff_rule(nh, nodes)

The rule a soft rule's stiff complement gets, in its whitened coordinates: the
`nodes`-point Gauss-Hermite product in `nh` dimensions for the standard normal,
as `(points, logweights)` with the weights summing to one. Built from
`_gauss_hermite`, whose cache is locked, so a thread may call it.

Why a full product rule. With the nodes held in the standardised effects, a
stiff direction's derivative in `theta` is a weighted average of per-node
scores that are large and cancel: the effect's posterior is narrow in `u`, so
a population scale moves a node by many posterior widths. How well they cancel
is how well the rule integrates, and two cheaper complements were measured on
the gated-gaps A14 config (three effects a subject, one soft direction each):

  - one node at the conditional mode. Its value is Laplace's, and right, but
    held fixed it leaves the complement's log determinant out of the
    derivative and its covariance out of the curvature: the continuation's
    standard errors came out at 0.07 to 0.45 of Laplace's.
  - the unscented transform's `2 nh + 1` points, exact to degree three and so
    on any Gaussian integrand, but not for the cross moments of a
    two-dimensional complement. At the Laplace optimum its derivative
    disagreed with that of the rule re-placed at every point by up to 280
    against 2.4, and the rounds stalled at once.
"""
function _continuation_stiff_rule(nh::Integer, nodes::Integer)
    nh == 0 && return (points=zeros(0, 1), logweights=[0.0])
    x, w = _gauss_hermite(nodes)
    m = Int(nodes)^Int(nh)
    points = Matrix{Float64}(undef, nh, m)
    logweights = Vector{Float64}(undef, m)
    for (col, index) in enumerate(Iterators.product(ntuple(_ -> 1:Int(nodes), Int(nh))...))
        for (d, i) in enumerate(index)
            points[d, col] = sqrt(2) * x[i]
        end
        logweights[col] = sum(log(w[i]) - log(pi) / 2 for i in index; init=0.0)
    end
    return (points=points, logweights=logweights)
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
each placement the continuation's rounds, which keep a step only if the
objective with the rule placed again is no lower, were rejected by the draws'
error on the round-9 fit: 0.4 to 0.7 nats at radii down to 0.03 standard
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
    rule = CTSEMFixedRule(b, points, logweights, -k * log(2 * pi) / 2, children, 0,
        replicate, replw)
    return (rule=rule, value=current.value, se=current.se, ess=current.ess,
        components=length(comps))
end

"""
    _continuation_block_rule(laplace, U, theta, Ls, b, u, context, aws, opts)

Place block `b`'s nodes, and those of every block beneath it, with its
ancestors at whatever `u` says, and score the block there. Returns
`(rule, value)` -- the fixed rule and its log integral at `theta`, which is the
quadrature value there -- or `nothing` when a mode or a precision could not be
found. See the module docstring for which rule a block gets.
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
    nsoft = 0
    if k <= opts.product_maxdim
        rule = if leaf
            _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws;
                start=Float64[context.mode[c] for c in columns])
        else
            clipped = _quadrature_clipped_scale(context.factors[b], k)
            (ok=true, centre=Float64[context.mode[c] for c in columns],
             scale=clipped.scale, logdetscale=clipped.logdetscale)
        end
        rule.ok || return nothing
        grid, lw = _gh_grid(k, opts.nodes)
        n = length(grid)
        points = Matrix{Float64}(undef, k, n)
        logweights = Vector{Float64}(undef, n)
        for j in 1:n
            z = rule.centre .+ sqrt(2) .* (rule.scale * grid[j])
            points[:, j] = z
            logweights[j] = lw[j] - dot(z, z) / 2
        end
        constant = rule.logdetscale + k * (log(2) - log(2 * pi)) / 2
    else
        fns = nothing
        if leaf
            lr = _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws;
                start=Float64[context.mode[c] for c in columns])
            lr.ok || return nothing
            centre = lr.centre
            fns = _continuation_block_functions(laplace, U, theta, Ls, b, u, aws)
            D = fns.precision(centre)
        else
            centre = Float64[context.mode[c] for c in columns]
            Lf = Matrix(context.factors[b].L)
            D = Lf * transpose(Lf)
        end
        E = _ctsem_symeig(D)
        all(isfinite, E.values) || return nothing
        lam = E.values
        nsoft = clamp(count(<(opts.soft_tau), lam), 1, min(opts.soft_maxdirs, k))
        Vs = E.vectors[:, 1:nsoft]
        Vh = E.vectors[:, (nsoft + 1):k]
        lsoft = max.(lam[1:nsoft], 1.0)
        lstiff = max.(lam[(nsoft + 1):k], 1.0)
        grid, lw = _gh_grid(nsoft, opts.nodes)
        stiff = _continuation_stiff_rule(k - nsoft, opts.nodes)
        m = size(stiff.points, 2)
        points = Matrix{Float64}(undef, k, length(grid) * m)
        logweights = Vector{Float64}(undef, length(grid) * m)
        col = 0
        for j in eachindex(grid)
            base = centre .+ Vs * (sqrt(2) .* grid[j] ./ sqrt.(lsoft))
            # Where the complement's nodes sit and the map from its whitened
            # coordinates: at a leaf its conditional mode and precision at this
            # soft node; at an outer block the joint mode and the eliminated
            # precision, whose eigenvectors `Vh` already are.
            mid, S, ld = if leaf
                c = _continuation_complement(fns, base, Vh, laplace)
                isfinite(c.logdet) || return nothing
                (base .+ Vh * c.h, Vh * c.vectors * Diagonal(1 ./ sqrt.(c.values)),
                 c.logdet)
            else
                (base, Vh * Diagonal(1 ./ sqrt.(lstiff)), sum(log, lstiff; init=0.0))
            end
            for i in 1:m
                col += 1
                xi = stiff.points[:, i]
                z = mid .+ S * xi
                points[:, col] = z
                logweights[col] = lw[j] + stiff.logweights[i] - ld / 2 +
                    dot(xi, xi) / 2 - dot(z, z) / 2
            end
        end
        # Per soft direction `sqrt(2 / lambda~)` from the change of variable,
        # and `(2 pi)^(-1/2)` from the prior's density. A complement node is
        # divided by the Gaussian it was placed by, `(2 pi)^(-nh/2) |P|^(1/2)
        # exp(-xi'xi/2)`: the `|P|` and `xi'xi` sit in its log weight, and the
        # `(2 pi)^(nh/2)` the division leaves cancels the rest of the prior's
        # density.
        constant = nsoft * (log(2) - log(2 * pi)) / 2 - sum(log, lsoft; init=0.0) / 2
    end
    n = size(points, 2)
    # Bound once for the closure below: `points` and `logweights` are assigned
    # in both branches above, and a closure over a local assigned twice boxes
    # it.
    nodes, nodeweights = points, logweights
    # In parallel over the nodes, each on its own copy of `u`, as an importance
    # rule's draws are (`_continuation_draw!`): a nested unit whose outer block
    # holds 125 nodes, each placing every member's rule beneath it, was placed
    # on one thread whenever it was the only unit being placed. Not below 64
    # member evaluations (`_continuation_ranges`' bar): a region opened for
    # every small block made test_continuation.jl three times slower.
    children = [CTSEMFixedRule[] for _ in 1:(leaf ? 0 : n)]
    terms = fill(-Inf, n)
    npar = length(theta)
    place_node = function (j)
        # `local` throughout, as in `_continuation_flagged`.
        local uj, awsj, inner, placed
        uj = copy(u)
        @inbounds for (t, c) in enumerate(columns); uj[c] = nodes[t, j]; end
        awsj = _laplace_workspace!(laplace, Float64, npar)
        inner = 0.0
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
            end
        end
        terms[j] = isfinite(inner) ? inner + nodeweights[j] : -Inf
        return true
    end
    placed_all = n * length(block.members) >= 64 ?
        _laplace_parallel(place_node, laplace, 1:n) : all(place_node, 1:n)
    placed_all || return nothing
    peak = maximum(terms)
    isfinite(peak) || return nothing
    value = peak + log(sum(t -> exp(t - peak), terms)) + constant
    return (rule=CTSEMFixedRule(b, points, logweights, constant, children, nsoft),
        value=value)
end

"""
    _continuation_unit_width(laplace, U)

The number of effects on the widest root-to-leaf path of unit `U`'s block tree:
the `d` whose `nodes^d` a rule over the unit costs per member.
"""
function _continuation_unit_width(laplace::CTSEMLaplaceObjective, U::Integer)
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return 0
    tree = _quadrature_children(blocks)
    return maximum((_continuation_block_width(blocks, tree, r) for r in tree.roots);
        init=0)
end
# Not a local recursive function: a closure that calls itself is a `Core.Box`.
_continuation_block_width(blocks, tree, b) = blocks[b].size +
    maximum((_continuation_block_width(blocks, tree, c) for c in tree.children[b]);
        init=0)

"""
    _continuation_unit_rule(laplace, U, theta, Ls, aws, opts)

Unit `U`'s fixed rule at `theta` and its quadrature value there, placed as
`_quadrature_chunk!` places the adaptive rule: the unit's mode from the
origin, its curvature, the block factorization, the recursion over the tree.
`nothing` for a unit with no random effects, or when the rule cannot be
placed.
"""
function _continuation_unit_rule(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, aws, opts)
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return nothing
    _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    mode = copy(laplace.modes[U])
    tree = _quadrature_children(blocks)
    if _continuation_unit_width(laplace, U) > opts.maxdim ||
            U in get(opts, :draw_units, ())
        # Too wide for a product rule, or named for draws because no product
        # rule reached the draws' target there (`_continuation_refine!`):
        # only `_continuation_place` sends such a unit here, and only when
        # draws are on. A nested unit's outer blocks take their first
        # components from the curvature, here, before anything runs in
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
    context = (children=tree.children, nodes=opts.nodes, mode=mode,
        factors=factorization.factors)
    u = copy(mode)
    roots = CTSEMFixedRule[]
    total = 0.0
    for root in tree.roots
        placed = _continuation_block_rule(laplace, U, theta, Ls, root, u, context,
            aws, opts)
        placed === nothing && return nothing
        push!(roots, placed.rule)
        total += placed.value
    end
    isfinite(total) || return nothing
    evaluations = sum(_continuation_evaluations(r, blocks) for r in roots; init=0)
    return (rule=CTSEMFixedUnitRule(U, roots, evaluations), value=total, se=0.0,
        importance=false)
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
    _continuation_warm_caches(laplace, opts.nodes, opts.soft_maxdirs, opts.product_maxdim)
    rules = Vector{Union{Nothing,CTSEMFixedUnitRule}}(nothing, nunits)
    values = fill(NaN, nunits)
    ses = zeros(nunits)
    importance = fill(false, nunits)
    # Wide and out of reach: a unit too wide for a product rule when draws are
    # off, or when they are kept for the units the caller named and it is not
    # one.
    drawable = get(opts, :drawable, nothing)
    wide = [_continuation_unit_width(laplace, U) > opts.maxdim &&
        !(opts.draws > 0 && (drawable === nothing || U in drawable)) for U in 1:nunits]
    only = get(opts, :only, nothing)
    _laplace_parallel(laplace, 1:nunits) do U
        wide[U] && return true
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
        laplace_value=lap.value, converged=lap.converged, wide=wide, se=ses,
        importance=importance)
end

"""
    _continuation_node_errors(laplace, theta, opts, placed, flagged)

Each flagged unit's rule error at the centre, in nats: a frozen rule's
standard error; for a product or soft rule the difference between its value
at `opts.nodes` and at two nodes fewer, placed again here (a mode solve and
a 3-node rule beside the 5-node rule already placed), `NaN` where that rule
cannot be placed or `nodes` is under 3. `NaN` for an unflagged unit, whose
term is Laplace's. What `ctsem_laplace_continuation_info` sums into the
hybrid's own error (`error`), the one number a fit reports for how far its
objective is from the exact marginal.
"""
function _continuation_node_errors(laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64}, opts, placed, flagged::Vector{Int})
    nunits = length(laplace.units.members)
    errors = fill(NaN, nunits)
    isempty(flagged) && return errors
    fewer = Int(opts.nodes) - 2
    opts3 = merge(opts, (nodes=max(fewer, 1),))
    Ls = _laplace_popchols(theta, laplace.spec)
    _continuation_warm_caches(laplace, opts3.nodes, opts.soft_maxdirs,
        opts.product_maxdim)
    _laplace_parallel(laplace, 1:length(flagged)) do i
        local U, r
        U = flagged[i]
        if placed.importance[U]
            errors[U] = placed.se[U]
        elseif fewer >= 1 && placed.rules[U] !== nothing
            r = try
                _continuation_unit_rule(laplace, U, theta, Ls,
                    _laplace_workspace!(laplace, Float64, length(theta)), opts3)
            catch err
                _ctsem_must_propagate(err) && rethrow()
                nothing
            end
            errors[U] = r === nothing ? NaN : abs(placed.quadrature[U] - r.value)
        end
        return true
    end
    return errors
end

"""
The most nodes a product or soft rule is refined to before a unit takes
draws (`_continuation_refine!`): two steps of two above the default 5, so
81 member evaluations per gradient on a two-effect block against 25. Past
that a unit whose rule has not converged is one whose posterior over its
effects is not one basin the nodes can cover, which draws from a defensive
mixture can.
"""
const _CONTINUATION_NODES_MAX = 9

"""
    _continuation_refine!(laplace, theta, opts, placed, flagged)

Each flagged unit's rule made honest where its nodes have not converged, and
every flagged unit's rule error at the centre (`_continuation_node_errors`).
A product or soft rule whose cheap check -- its value against the rule two
nodes fewer -- is within `opts.draws_se` keeps its rule and that error. One
beyond it is checked against the rule two nodes more: the finer rule's
difference from it is the rule's error, and where that too is beyond
`draws_se` the finer rule is taken and checked the same way, up to
`_CONTINUATION_NODES_MAX` nodes; a unit no rule up to there resolves takes
draws, which stop at `draws_se` on their own standard error
(`_continuation_importance_rule`), or keeps the finest rule with its error
where draws are off. `placed`'s arrays are updated in place.

Why: a hybrid is only as honest as its rules. On the bench's gC8 the 5-node
rules carried 6.9 nats of node error at the climb's end, and on anomS2's
spurious start 54 nats over 101 units: the climb on that hybrid walked where
the rules over-credited the Laplace basin and ended 11 and 15 nats below the
baseline, certified on an objective whose own error was that large. The
cheap check first, because on a unit whose posterior is one Gaussian basin
3 nodes already agree with 5, and the finer rule costs `(nodes + 2)^k`
evaluations a unit a placement.
"""
function _continuation_refine!(laplace::CTSEMLaplaceObjective, theta::Vector{Float64},
    opts, placed, flagged::Vector{Int})
    errors = _continuation_node_errors(laplace, theta, opts, placed, flagged)
    target = Float64(opts.draws_se)
    coarse = [U for U in flagged if !placed.importance[U] &&
        placed.rules[U] !== nothing && isfinite(errors[U]) && errors[U] > target]
    isempty(coarse) && return errors
    Ls = _laplace_popchols(theta, laplace.spec)
    # Every node count the loop below may ask for, warmed here: the grid
    # caches are pure reads inside the threads (`_quadrature_warm_caches`).
    for m in (Int(opts.nodes) + 2):2:_CONTINUATION_NODES_MAX
        _continuation_warm_caches(laplace, m, opts.soft_maxdirs, opts.product_maxdim)
    end
    _laplace_parallel(laplace, 1:length(coarse)) do i
        local U, n, finer, err, r, ws
        U = coarse[i]
        ws = _laplace_workspace!(laplace, Float64, length(theta))
        n = Int(opts.nodes)
        finer = false
        err = errors[U]
        while n + 2 <= _CONTINUATION_NODES_MAX
            r = try
                _continuation_unit_rule(laplace, U, theta, Ls, ws,
                    merge(opts, (nodes=n + 2,)))
            catch e
                _ctsem_must_propagate(e) && rethrow()
                nothing
            end
            r === nothing && break
            err = abs(r.value - placed.quadrature[U])
            if err <= target
                # The coarser rule stands, with the finer one's measure of it.
                errors[U] = err
                finer = false
                break
            end
            placed.rules[U] = r.rule
            placed.quadrature[U] = r.value
            errors[U] = err
            finer = true
            n += 2
        end
        if finer && errors[U] > target && opts.draws > 0
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
            end
        end
        return true
    end
    return errors
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
        product_maxdim=2, soft_tau=3.5, soft_maxdirs=2, maxdim=5, draws=0,
        draws_se=0.05, draws_max=16384, draws_floor=0.25, draws_prior=0.2,
        draws_df=4, draws_seed=20261021, force=nothing, point_units=nothing,
        point_values=nothing, screen=true, draws_scope=:all)

The hybrid objective with its nodes placed at `values`. Every unit is scored by
the rule there against its Laplace term, the units that differ beyond the
rule's own error are flagged (`_continuation_resolved`, `_continuation_flags`),
and those carry fixed nodes from then on while the rest keep their Laplace
term; the units in `force` are flagged whatever their gap. A unit wider than
`maxdim` effects along any path of its block tree gets a frozen importance
rule, per block of its tree, when `draws > 0`, and is not scored at all
otherwise, keeping its Laplace term. `point_units` and `point_values` are
points the caller knows for some units, column `j` of `point_values` for unit
`point_units[j]` (its first `dims` entries), each a component of that unit's
importance proposal. `screen = false` scores only the units `force` names, the
rest keeping their Laplace terms unscored: what a fit climbing the hybrid
needs. `draws_scope = :forced` gives draws only to the units `force` names,
every other unit too wide for a product rule keeping its Laplace term
unscored, as with `draws = 0`. See the module docstring for the rules and the
constants.
"""
function ctsem_laplace_continuation(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, tolerance::Real=0.01,
    product_maxdim::Integer=2, soft_tau::Real=3.5, soft_maxdirs::Integer=2,
    maxdim::Integer=5, draws::Integer=0, draws_se::Real=0.05,
    draws_max::Integer=16384, draws_floor::Real=0.25, draws_prior::Real=0.2,
    draws_df::Real=4, draws_seed::Integer=20261021, force=nothing,
    point_units=nothing, point_values=nothing, screen::Bool=true, draws_scope=:all)
    nodes >= 1 || throw(ArgumentError("need at least one quadrature node"))
    draws_df > 0 || throw(ArgumentError("draws_df must be positive (Inf for Gaussian)"))
    Symbol(draws_scope) in (:all, :forced) || throw(ArgumentError(
        "draws_scope must be :all or :forced"))
    draws >= 0 || throw(ArgumentError("draws must be at least zero"))
    draws_se > 0 || throw(ArgumentError("draws_se must be positive"))
    0 < draws_prior < 1 || throw(ArgumentError("draws_prior must be between 0 and 1"))
    draws_floor > 0 || throw(ArgumentError("draws_floor must be positive"))
    theta = _laplace_extend(laplace, collect(Float64, values))
    nunits = length(laplace.units.members)
    o = CTSEMLaplaceContinuation{typeof(laplace)}(laplace, laplace, collect(1:nunits),
        Int[], CTSEMFixedUnitRule[], theta, Int(nodes), Float64(tolerance),
        Int(product_maxdim), Float64(soft_tau), Int(soft_maxdirs), Int(maxdim),
        fill(false, nunits), Int(draws), Float64(draws_se),
        max(Int(draws_max), Int(draws)), Float64(draws_floor), Float64(draws_prior),
        Float64(draws_df), Int(draws_seed), Dict{Int,Vector{Vector{Float64}}}(),
        zeros(nunits), fill(NaN, nunits),
        fill(false, nunits), Bool(screen), Symbol(draws_scope), Int[],
        fill(NaN, nunits), fill(NaN, nunits), fill(false, nunits), NaN, nothing,
        0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, nothing, nothing, nothing,
        Float64[], Float64[])
    _continuation_add_points!(o, point_units, point_values)
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

"""Points the caller knows for some units (`point_units[j]`, `point_values[:,
j]`), added to `o.points` as components for their importance proposals."""
function _continuation_add_points!(o::CTSEMLaplaceContinuation, point_units,
    point_values)
    (point_units === nothing || point_values === nothing) && return o
    V = Matrix{Float64}(point_values)
    units = Int[Int(U) for U in point_units]
    size(V, 2) == length(units) || throw(DimensionMismatch(
        "point_values needs one column per entry of point_units"))
    for (j, U) in enumerate(units)
        1 <= U <= length(o.laplace.units.members) || continue
        d = o.laplace.units.dims[U]
        d <= size(V, 1) || continue
        x = V[1:d, j]
        all(isfinite, x) || continue
        push!(get!(o.points, U, Vector{Vector{Float64}}()), x)
    end
    return o
end

const _CONTINUATION_STATE = (:rest, :rest_units, :flagged, :rules, :centre,
    :centre_quadrature, :centre_laplace, :rule_failed, :centre_value, :centre_se,
    :centre_error, :importance, :ess, :ess_placed)

"""
    ctsem_laplace_continuation_recentre!(o, values; force=nothing,
        point_units=nothing, point_values=nothing, screen=nothing)

Re-place `o`'s nodes at `values`: every unit is scored again, the flagged set
grows by any unit whose gap now matters and by any `force` names -- it never
shrinks, so no later round loses a unit an earlier one was scored on -- and
every flagged unit's rule moves to its mode and curvature there, an importance
rule drawing afresh with the points it knows and any given here. The placement
it replaces is kept for `ctsem_laplace_continuation_revert!`. Returns
`ctsem_laplace_continuation_info`.
"""
function ctsem_laplace_continuation_recentre!(o::CTSEMLaplaceContinuation,
    values::AbstractVector; force=nothing, point_units=nothing, point_values=nothing,
    screen=nothing)
    o.previous = NamedTuple{_CONTINUATION_STATE}(
        Tuple(getfield(o, f) for f in _CONTINUATION_STATE))
    # `screen` given turns the screen on or off from here: a fit that climbed
    # with the screen off places every unit once at its end.
    screen === nothing || (o.screen = Bool(screen))
    _continuation_add_points!(o, point_units, point_values)
    _continuation_force!(o, force)
    _continuation_place!(o, _laplace_extend(o.laplace, collect(Float64, values)))
    return ctsem_laplace_continuation_info(o)
end

"""
    ctsem_laplace_continuation_revert!(o)

Put back the placement the last `ctsem_laplace_continuation_recentre!`
replaced: the round that moved there was rejected, so the rule goes back to
the point it was placed at. Nothing is evaluated.
"""
function ctsem_laplace_continuation_revert!(o::CTSEMLaplaceContinuation)
    o.previous === nothing && return ctsem_laplace_continuation_info(o)
    for f in _CONTINUATION_STATE
        setfield!(o, f, getfield(o.previous, f))
    end
    o.previous = nothing
    o.cache = nothing
    o.value_cache = nothing
    return ctsem_laplace_continuation_info(o)
end

function _continuation_place!(o::CTSEMLaplaceContinuation, theta::Vector{Float64})
    started = time_ns()
    laplace = o.laplace
    opts = (nodes=o.nodes, product_maxdim=o.product_maxdim, soft_tau=o.soft_tau,
        maxdim=o.maxdim,
        soft_maxdirs=o.soft_maxdirs, draws=o.draws, draws_se=o.draws_se,
        draws_max=o.draws_max, draws_floor=o.draws_floor, draws_prior=o.draws_prior,
        draws_df=o.draws_df, draws_seed=o.draws_seed, points=o.points,
        only=o.screen ? nothing : Set(o.flagged),
        drawable=o.draws_scope === :all ? nothing : Set(o.forced))
    placed = _continuation_place(laplace, theta, opts)
    nunits = length(laplace.units.members)
    failed = [placed.rules[U] === nothing && !isempty(laplace.units.blocks[U]) &&
              !placed.wide[U] for U in 1:nunits]
    o.wide = placed.wide
    gaps = placed.quadrature .- placed.laplace
    flagged = sort!(union(o.flagged, _continuation_flags(gaps, o.tolerance;
        se=placed.se)))
    # A unit flagged before whose rule could not be placed here keeps its
    # Laplace term: an honest Laplace term beats a rule placed somewhere else.
    filter!(U -> placed.rules[U] !== nothing, flagged)
    errors = _continuation_refine!(laplace, theta, opts, placed, flagged)
    o.flagged = flagged
    o.rules = CTSEMFixedUnitRule[placed.rules[U] for U in flagged]
    o.rest_units = setdiff(1:nunits, flagged)
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
    o.recentres += 1
    o.ess = fill(NaN, length(o.rules))
    o.ess_placed = fill(NaN, length(o.rules))
    o.cache = nothing
    o.value_cache = nothing
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
    _continuation_record_ess!(o, f, ess)

Keep flagged rule `f`'s effective number of nodes at the last evaluation
(`o.ess`), and the first one after its placement (`o.ess_placed`), which
`_continuation_drifted` compares: a rule whose nodes
were placed at one point carries less and less of the integrand's mass as
the parameters move, and its weights concentrate on the nodes nearest the
new mode. Written in the serial part of `_continuation_flagged`.
"""
function _continuation_record_ess!(o::CTSEMLaplaceContinuation, f::Integer, ess::Real)
    f <= length(o.ess) || return nothing
    o.ess[f] = Float64(ess)
    isnan(o.ess_placed[f]) && isfinite(ess) && (o.ess_placed[f] = Float64(ess))
    return nothing
end

"""
    _continuation_drifted(o, ratio)

Whether any flagged rule's effective number of nodes has fallen under `ratio`
of what it was at its placement (`_continuation_record_ess!`): nodes fixed
where they were placed carry less of the integrand as the parameters move,
and read far from home they can say anything. A round of the correction ends
on it (`ctsem_laplace_continuation_optimize`), which is as far as a run on
frozen rules can be trusted: unbounded, a first round on the bench's gA14
ran 40 standard errors out and was kept 0.3 exact nats below where bounded
rounds ended.
"""
_continuation_drifted(o::CTSEMLaplaceContinuation, ratio::Real) =
    any(f -> isfinite(o.ess[f]) && isfinite(o.ess_placed[f]) &&
        o.ess[f] < ratio * o.ess_placed[f], eachindex(o.rules))

# A round's stop: sets `flag` so the round can say it ended on drift.
_continuation_drift_stop!(o::CTSEMLaplaceContinuation, ratio::Real, flag::Base.RefValue{Bool}) =
    _continuation_drifted(o, ratio) ? (flag[] = true; true) : false

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
summed in order afterwards."""
function _continuation_flagged(o::CTSEMLaplaceContinuation, theta::Vector{Float64},
    want_gradient::Bool)
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
    _continuation_warm_caches(laplace, o.nodes, o.soft_maxdirs, o.product_maxdim)
    tasks = _continuation_tasks(o, _laplace_pool_width())
    ntask = length(tasks)
    peaks = fill(-Inf, ntask)
    accs = zeros(ntask)
    # Per task the sum of squared weights (a range) or the effective number
    # itself (a whole unit, the smallest over its roots).
    acc2s = zeros(ntask)
    esss = fill(NaN, ntask)
    grads = [want_gradient ? zeros(npar) : Float64[] for _ in 1:ntask]
    ok = _laplace_parallel(laplace, 1:ntask) do k
        # `local` throughout: a name assigned here that is also a local of the
        # enclosing function would be one binding shared by every task.
        local f, range, rule, U, aws, scratch, u, total, r, failed, a, essmin
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
                essmin = Inf
                for root in rule.roots
                    r = _continuation_rule_value(o, U, theta, Ls, dL, positions,
                        root, u, aws, scratch, want_gradient)
                    if !isfinite(r.value)
                        failed = true
                        break
                    end
                    total += r.value
                    essmin = min(essmin, r.ess)
                    want_gradient && (grads[k] .+= r.gradient)
                end
                peaks[k] = total
                accs[k] = 1.0
                esss[k] = essmin
            else
                a = _continuation_node_sum(o, U, theta, Ls, dL, positions,
                    rule.roots[1], u, aws, scratch, want_gradient, range)
                failed = a.failed
                peaks[k] = a.peak
                accs[k] = a.acc
                acc2s[k] = a.acc2
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
    at = 1
    while at <= ntask
        f = tasks[at][1]
        last = at
        while last < ntask && tasks[last + 1][1] == f
            last += 1
        end
        rule = o.rules[f]
        if tasks[at][2] === nothing
            values[f] = peaks[at]
            _continuation_record_ess!(o, f, esss[at])
            want_gradient && (g .= grads[at])
        else
            P = maximum(@view peaks[at:last])
            isfinite(P) || return (values=values, gradient=gradient, ok=false)
            A = 0.0
            B = 0.0
            want_gradient && fill!(g, 0.0)
            for q in at:last
                e = exp(peaks[q] - P)
                A += accs[q] * e
                B += acc2s[q] * e * e
                want_gradient && (g .+= e .* grads[q])
            end
            values[f] = P + log(A) + rule.roots[1].constant
            _continuation_record_ess!(o, f, A^2 / max(B, eps()))
            want_gradient && (g ./= A)
        end
        if want_gradient
            all(isfinite, g) || return (values=values, gradient=gradient, ok=false)
            gradient .+= _laplace_unit_weight(laplace, rule.unit) .* g
        end
        at = last + 1
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
    for hit in (o.cache, o.value_cache)
        if hit !== nothing && hit.theta == theta &&
                (!gradient || hit.result.gradient !== nothing)
            return _laplace_free_result(_continuation_copy(hit.result), n)
        end
    end
    result = _continuation_evaluate(o, theta, gradient)
    entry = (theta=theta, result=_continuation_copy(result))
    gradient ? (o.cache = entry) : (o.value_cache = entry)
    return _laplace_free_result(result, n)
end

# Copies of the vectors, so that neither a caller nor the cache can change
# what the other holds.
_continuation_copy(r) = (value=r.value,
    gradient=r.gradient === nothing ? nothing : copy(r.gradient),
    unit_loglik=copy(r.unit_loglik), subject_loglik=copy(r.subject_loglik),
    converged=r.converged)

function _continuation_evaluate(o::CTSEMLaplaceContinuation, theta::Vector{Float64},
    gradient::Bool)
    started = time_ns()
    rest = ctsem_laplace_evaluate(o.rest, theta; gradient=gradient)
    flagged = _continuation_flagged(o, theta, gradient)
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

# The optimisable protocol (see ctsem_backend.jl), so `ctsem_optimize` and the
# flat probe take the hybrid as they take any objective.

function ctsem_evaluate(o::CTSEMLaplaceContinuation, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    r = ctsem_laplace_continuation_evaluate(o, values; gradient=gradient)
    contributions || return (value=r.value, gradient=r.gradient)
    return (value=r.value, gradient=r.gradient, subject_loglik=r.subject_loglik)
end

_ctsem_optimise_label(::CTSEMLaplaceContinuation) = "Laplace continuation"
# The finish walks the hybrid on its own Hessian, as it walks the Laplace
# objective, where the optimiser takes the hybrid whole: chord, at the cost
# of the Laplace part's differences.
_ctsem_finish_curvature(::CTSEMLaplaceContinuation) = :chord
_ctsem_hessian_cost(::CTSEMLaplaceContinuation, npar::Integer) = 2 * Int(npar)
# A continuation built where a Laplace run ended starts from that run's
# pairs: the hybrid with nothing flagged is the Laplace objective, and with
# units flagged it differs by their terms, so the pairs still give descent
# directions. The Laplace objective keeps a run's pairs in a weak table, so
# an objective that is let go takes them with it.
const _CTSEM_LBFGS_HELD = WeakKeyDict{Any,CTSEMLBFGSMemory}()
_ctsem_lbfgs_held(o::CTSEMLaplaceObjective) = get(_CTSEM_LBFGS_HELD, o, nothing)
_ctsem_lbfgs_hold!(o::CTSEMLaplaceObjective, memory::CTSEMLBFGSMemory) =
    (_CTSEM_LBFGS_HELD[o] = memory; nothing)
_ctsem_lbfgs_held(o::CTSEMLaplaceContinuation) =
    o.lbfgs === nothing ? _ctsem_lbfgs_held(o.laplace) : o.lbfgs
_ctsem_lbfgs_hold!(o::CTSEMLaplaceContinuation, memory::CTSEMLBFGSMemory) =
    (o.lbfgs = memory; nothing)
ctsem_hessian(o::CTSEMLaplaceContinuation, values::AbstractVector;
    chunk::Integer=0, progress=nothing) = ctsem_laplace_continuation_hessian(o, values)
_ctsem_params(o::CTSEMLaplaceContinuation) = _ctsem_params(o.laplace)
_ctsem_saturated_for(o::CTSEMLaplaceContinuation, minimizer) =
    _laplace_saturated_parameters(o.laplace, minimizer)
# The flagged count, for a run's record.
_ctsem_optimise_result_extra(o::CTSEMLaplaceContinuation, final, log) =
    (continuation_flagged=length(o.flagged),)

function _ctsem_optimise_trial(o::CTSEMLaplaceContinuation, x, want_gradient::Bool,
    gradient_method, limit::Real, log)
    evaluated = _ctsem_try_evaluate(() ->
        ctsem_laplace_continuation_evaluate(o, x; gradient=want_gradient))
    valid = evaluated !== nothing && isfinite(evaluated.value) && evaluated.converged &&
        (!want_gradient || _ctsem_gradient_usable(evaluated.gradient, limit))
    return (evaluated=evaluated, valid=valid)
end

# The estimate read the same way, so a probe's gain is never the placement
# difference (`_ctsem_probe_baseline`).
_ctsem_probe_baseline(o::CTSEMLaplaceContinuation, minimizer, value) =
    _ctsem_probe_value(o, minimizer)

# A probe asks whether the objective is better somewhere else -- the
# overshoot probe's pullbacks, the stall watch's, the finish's flat and
# saddle ladders -- and on the hybrid the objective at another point is its
# value with the rules placed there, not the frozen rules' value far from
# their placement: a rule's nodes fixed at one point carry less and less of
# the integrand as the parameters move, and read far away they can say
# anything. On the bench's anomS2 the climb reached the baseline's optimum,
# the probe on its frozen rules found a gain of 0.016 nats toward the basin
# the Laplace terms over-credit, and the escape taken on it ended the fit 2.5
# exact nats and 6.5 reported nats off. So the rules are placed at the point
# and the placement put back, around the correction rounds' own saved
# placement (`o.previous`), which a probe inside a round must not disturb.
function _ctsem_probe_value(o::CTSEMLaplaceContinuation, x)
    theta = try
        _laplace_extend(o.laplace, collect(Float64, x))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        return -Inf
    end
    saved = NamedTuple{_CONTINUATION_STATE}(
        Tuple(getfield(o, f) for f in _CONTINUATION_STATE))
    previous = o.previous
    # Not a placement of the correction's own: the count of placements the
    # record reports stays as it was.
    recentres = o.recentres
    value = -Inf
    try
        _continuation_place!(o, theta)
        # What the hybrid can vouch for there: its value less its own error at
        # that placement (`_continuation_error`), as the estimate is read
        # (`_ctsem_probe_baseline`). Read whole, a probe four raw units from
        # round 9's estimate found 13000 nats on a unit whose draws had
        # degenerated there, and the fit was certified not converged on it.
        value = o.centre_value - _continuation_error(o).total
    catch err
        _ctsem_must_propagate(err) && rethrow()
        value = -Inf
    finally
        for f in _CONTINUATION_STATE
            setfield!(o, f, getfield(saved, f))
        end
        o.previous = previous
        o.recentres = recentres
        o.cache = nothing
        o.value_cache = nothing
    end
    return isfinite(value) ? value : -Inf
end

"""
    CTSEMContinuationStep(o, centre, basis, radius)

One round's problem: the hybrid `o` over `x = centre + basis * y`, with `y`
confined to the ball `|y| <= radius`. The caller builds `basis` from the
Laplace fit's curvature -- its identified directions, each scaled to one
standard error -- so the round's L-BFGS starts with that curvature as its
metric, which is what makes the hand-over from the Laplace fit a short climb,
and the ball is a trust region measured in standard errors. A trial point
outside the ball is refused as invalid, and the line search shrinks. The
directions `basis` leaves out are the ones the Laplace fit says the data do
not identify, held where the fit left them.
"""
struct CTSEMContinuationStep{C} <: CTSEMOptimisable
    cont::C
    centre::Vector{Float64}
    basis::Matrix{Float64}
    radius::Float64
end

_continuation_point(s::CTSEMContinuationStep, y) = s.centre .+ s.basis * collect(Float64, y)
_continuation_inside(s::CTSEMContinuationStep, y) =
    s.radius <= 0 || sqrt(sum(abs2, y)) <= s.radius * (1 + 1e-12)

function ctsem_evaluate(s::CTSEMContinuationStep, y::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    r = ctsem_evaluate(s.cont, _continuation_point(s, y); gradient=gradient,
        contributions=contributions)
    g = gradient ? transpose(s.basis) * r.gradient : nothing
    contributions || return (value=r.value, gradient=g)
    return (value=r.value, gradient=g, subject_loglik=r.subject_loglik)
end

_ctsem_optimise_label(::CTSEMContinuationStep) = "Laplace continuation"
_ctsem_params(s::CTSEMContinuationStep) = _ctsem_params(s.cont)
# Nothing in `y` is a transform's coordinate; saturation is the caller's to
# judge, in the parameters, at the point the round ends.
_ctsem_saturation_range(::CTSEMContinuationStep, minimizer) = Int[]
_ctsem_saturated_for(::CTSEMContinuationStep, minimizer) = Int[]
_ctsem_optimise_result_extra(::CTSEMContinuationStep, final, log) = NamedTuple()

function _ctsem_optimise_trial(s::CTSEMContinuationStep, y, want_gradient::Bool,
    gradient_method, limit::Real, log)
    if !_continuation_inside(s, y)
        s.cont.refused += 1
        return (evaluated=nothing, valid=false)
    end
    trial = _ctsem_optimise_trial(s.cont, _continuation_point(s, y), want_gradient,
        gradient_method, limit, log)
    evaluated = trial.evaluated
    if evaluated !== nothing && want_gradient && evaluated.gradient !== nothing
        evaluated = merge(evaluated, (gradient=transpose(s.basis) * evaluated.gradient,))
    end
    return (evaluated=evaluated, valid=trial.valid)
end

_ctsem_probe_value(s::CTSEMContinuationStep, y) = _continuation_inside(s, y) ?
    _ctsem_probe_value(s.cont, _continuation_point(s, y)) : -Inf

"""
    ctsem_laplace_continuation_optimize(o, values, basis, radius; maxiter=100,
        tol=1e-6, stationary_only=false)

One round: maximise the hybrid, its nodes where they are, over
`values + basis * y` with `|y| <= radius`, by the engine's own `ctsem_optimize`.
The first step is the whitened Newton step `basis' g` (under the curvature
`basis` was built from), cut to the radius. Returns the point reached in the
parameters, the hybrid's value there and at the start, the optimiser's counts,
`moved` (`|y|`) and whether the round ended on the boundary of its region.

`start_gain` is `|basis' g|^2 / 2`, the gain that Newton step predicts, and a
round whose `start_gain` is below `tol` is not run: `stationary` says the point
is already the answer on these nodes. `stationary_only` asks only that.
"""
# One round of the correction: L-BFGS on the frozen rules from `values`,
# in the directions `basis` spans, within `radius` of the start (0 for no
# bound), and only while the rules have not drifted (`drift`, the share of a
# rule's placed effective node count under which the round ends,
# `_continuation_drifted`; 0 for never). `drifted` in the result says it
# ended so.
function ctsem_laplace_continuation_optimize(o::CTSEMLaplaceContinuation,
    values::AbstractVector, basis::AbstractMatrix, radius::Real;
    maxiter::Integer=100, tol::Real=1e-6, stationary_only::Bool=false,
    drift::Real=0.5)
    x0 = collect(Float64, values)
    B = Matrix{Float64}(basis)
    size(B, 1) == length(x0) || throw(DimensionMismatch(
        "the basis must have one row per parameter"))
    start = ctsem_laplace_continuation_evaluate(o, x0; gradient=true)
    usable = size(B, 2) > 0 && isfinite(start.value) && start.converged &&
        all(isfinite, start.gradient)
    gy = usable ? transpose(B) * start.gradient : zeros(0)
    start_gain = usable ? sum(abs2, gy) / 2 : NaN
    stationary = usable && start_gain < tol
    none = (minimizer=x0, value=start.value, start_value=start.value, iterations=0,
        f_calls=0, g_calls=1, converged=stationary, stopped_by_gap=false, moved=0.0,
        boundary=false, drifted=false, gradient_norm=NaN, refused=0,
        start_gain=start_gain, stationary=stationary)
    (usable && !stationary && !stationary_only) || return none
    newton = sqrt(sum(abs2, gy))
    alpha = radius > 0 ? min(newton, Float64(radius)) : newton
    step = CTSEMContinuationStep(o, x0, B, Float64(radius))
    refused = o.refused
    drifted = Ref(false)
    ratio = Float64(drift)
    stop = ratio > 0 ? (st -> _continuation_drift_stop!(o, ratio, drifted)) : nothing
    res = ctsem_optimize(step, zeros(size(B, 2)); maxiter=Int(maxiter), g_tol=1e-8,
        f_tol=0.0, x_tol=0.0, gap_tol=Float64(tol), converge_tol=Float64(tol),
        overshoot_probe=:off, stall_window=0, batch=false, newton=false,
        precondition=nothing, lbfgs_diagonal=false,
        # Its basis is whitened by the Hessian, so its first step is Newton's
        # and its prediction is the gain left.
        lbfgs_scaled=true,
        initial_alpha=max(alpha, 1e-12), progress=false, stop=stop)
    y = collect(Float64, res.minimizer)
    moved = sqrt(sum(abs2, y))
    return (minimizer=_continuation_point(step, y), value=res.maximum_loglik,
        start_value=start.value, iterations=res.iterations, f_calls=res.f_calls,
        g_calls=res.g_calls + 1, converged=res.converged,
        stopped_by_gap=res.stopped_by_gap, moved=moved,
        boundary=radius > 0 && moved >= 0.9 * radius, drifted=drifted[],
        gradient_norm=res.gradient_norm,
        refused=o.refused - refused, start_gain=start_gain, stationary=false)
end

"""
    ctsem_laplace_continuation_hessian(o, values; step=1e-4, scheme=:exact,
        width=0)

The Hessian of the hybrid at `values`, its nodes where they are.

`scheme = :exact`, the default, takes the flagged units' part exactly, by the
Louis identity over their fixed nodes, and the part the other units and the
prior carry by central differences of that part's own gradient
(`_continuation_exact_hessian`). Where the flagged part cannot be formed that
way it falls back to `:forward`. `width` is the dual width of its sweeps, `0`
to choose (`_continuation_hessian_width`); the result does not depend on it.

`:forward` and `:central` difference the hybrid's exact gradient, the Laplace
part's inner modes warm-started from their modes at `values` for the reason
`ctsem_laplace_hessian` gives: `:forward` one gradient per column and one at
`values`, `:central` two per column -- `npar + 1` and `2 npar` reverse sweeps
per member per node, where the exact scheme takes `ceil(active / width)` dual
ones. Measured on dev1, each after a first call: on the bench cell ord4 (27
parameters, 49 of 50 units flagged, 625 nodes each; one thread) exact 175 s,
forward 825 s, central 1593 s, the exact standard errors within 2e-7 of the
central ones; on dev/lapcontinue/highdim.R (150 subjects; 8 threads) at five
effects a unit exact 8.7 s against forward 125 s, and at four, where the width
stays one, 10.6 s against 20.3 s. A session's first exact Hessian adds its
compilation: 30 s and 15 s in all there. Forward was the default before this
scheme (e258ed65), central before that; on the optimiser bench's Laplace fits
whose continuation reported this Hessian (32 fits of 9 models;
review/bench/2026-09-26-baseline.md) forward moved no standard error by more
than 0.24% against central. The Laplace Hessian is no substitute: at the
Laplace optimum or at the continuation's estimate it moved some standard
error by more than 20% on three of the nine, by up to 69%.

A column whose points did not all evaluate is `NaN`, under either difference
scheme and in the exact scheme's differenced part, which the caller sees as an
unusable Hessian.
"""
function ctsem_laplace_continuation_hessian(o::CTSEMLaplaceContinuation,
    values::AbstractVector; step::Real=1e-4, scheme=:exact, width::Integer=0)
    scheme = scheme isa Symbol ? scheme : Symbol(scheme)
    scheme in (:exact, :forward, :central) || throw(ArgumentError(
        "scheme must be :exact, :forward or :central, got " * repr(scheme)))
    # Over the extended vector (`_laplace_extend`), returned at the caller's.
    nin = length(values)
    x = _laplace_extend(o.laplace, collect(Float64, values))
    started = time_ns()
    H = scheme === :exact ? _continuation_exact_hessian(o, x; step=step, width=width) :
        nothing
    if H === nothing
        # A difference scheme's gradients are in `seconds_gradients` as well.
        H = _continuation_difference_hessian(o, x, step,
            scheme === :exact ? :forward : scheme)
    end
    o.seconds_hessian += (time_ns() - started) / 1e9
    return size(H, 1) == nin ? H : H[1:nin, 1:nin]
end

"""The `:forward` and `:central` schemes of `ctsem_laplace_continuation_hessian`."""
function _continuation_difference_hessian(o::CTSEMLaplaceContinuation,
    x::Vector{Float64}, step::Real, scheme::Symbol)
    n = length(x)
    H = fill(NaN, n, n)
    previous = ctsem_set_warm_start!(false)
    try
        ctsem_laplace_evaluate(o.rest, x; gradient=false)
    finally
        ctsem_set_warm_start!(previous)
    end
    base = deepcopy(o.rest.modes)
    restore!() = (for U in eachindex(base); o.rest.modes[U] = copy(base[U]); end)
    previous = ctsem_set_warm_start!(true)
    try
        # The gradient at `values` itself, which the forward scheme differences
        # every column against.
        g0 = if scheme === :forward
            # Fresh, not cached: every difference starts from the same modes.
            restore!()
            _continuation_evaluate(o, x, true)
        else
            nothing
        end
        g0 === nothing || g0.converged || return H
        for j in 1:n
            h = step * max(1.0, abs(x[j]))
            plus = copy(x); plus[j] += h
            restore!()
            gp = _continuation_evaluate(o, plus, true)
            if g0 !== nothing
                gp.converged || continue
                H[:, j] = (gp.gradient .- g0.gradient) ./ h
            else
                minus = copy(x); minus[j] -= h
                restore!()
                gm = _continuation_evaluate(o, minus, true)
                (gp.converged && gm.converged) || continue
                H[:, j] = (gp.gradient .- gm.gradient) ./ (2h)
            end
        end
    finally
        ctsem_set_warm_start!(previous)
        restore!()
    end
    return (H .+ transpose(H)) ./ 2
end

# The exact Hessian
#
# A flagged unit's term is `Q(theta) = log sum_j exp(w_j + h_j(theta)) + c`
# with the nodes held, so with `pi` the normalised node weights
#
#     d2Q/dtheta2 = sum_j pi_j d2h_j + sum_j pi_j (dh_j - dQ)(dh_j - dQ)'
#
# -- the Louis identity: the node-weighted mean of the per-node Hessians, plus
# the node-weighted covariance of the per-node gradients, which are the Fisher
# identity's scores that the gradient already takes. An outer block's `h_j` is
# the sum of its children's rules as they were placed at node `j`, so the same
# identity applies at every level of the unit's block tree, children first.
#
# A leaf's `h_j` is its members' log likelihoods at `theta + L(theta) u_j`. For
# one member, with `s(theta)` that shifted vector, `g` and `Hs` the gradient and
# Hessian of its log likelihood in `s`,
#
#     d2h/dtheta2 = J' Hs J + sum_p g_rho(p) d2(L u)_p/dtheta2,    J = ds/dtheta
#
# `J` is the identity plus, in the column of each population parameter, the
# vector `dL u` placed at the level's `re_index`; the last term has entries
# only between population parameters of one level. The filter never reads a
# population parameter, so `Hs` is needed only over the others (`active`), and
# it comes from one forward-over-reverse sweep: the member's value-and-gradient
# with dual inputs, each seeded along one active parameter, `width` of them to
# a sweep (`_continuation_member_hessian!`). Its value parts are the plain
# sweep's, so the gradient comes with it, and a node costs
# `ceil(active / width)` dual sweeps where differencing the gradient costs
# `npar + 1` or `2 npar` plain ones.
#
# The units left to Laplace, and the prior, which that part carries, have no
# exact Hessian in the engine: their term's gradient is itself built from
# third derivatives at an inner mode. They are differenced as the Laplace fit
# differences them, centrally, over those units alone (`ctsem_laplace_hessian`
# on the subset objective), which is where every unit is when nothing is
# flagged and nowhere much when most are.

"""
Which raw parameters a member's filter reads (`active`, the ones a sweep
seeds) and which are population parameters (`pop`, reaching it only through
`L(theta)`), with each level's `re_index` in the active numbering (`reloc`) and
one `(level, position, parameter)` entry per population parameter that moves a
member (`entries`), and the `positions` it was built from. `nothing` when a
random effect sits on a population parameter, which the assembly does not
handle.

`filtered` are the level parameters a member's filter reads itself
(`_hybrid_filtered`): own variances of deviations the filter integrates, no
part of `L`, whose entries' `dL` is zero; or a block's scales and
correlations, which reach the member through `L` as well. They are seeded
with the active ones.
"""
function _continuation_hessian_layout(spec::CTSEMLaplaceSpec, npar::Integer,
    positions::Vector{Vector{Int}}, dL, filtered::Vector{Int}=Int[])
    pop = setdiff!(sort!(unique!(reduce(vcat, positions; init=Int[]))), filtered)
    active = Int[]
    local_of = zeros(Int, npar)
    for t in 1:npar
        insorted(t, pop) && continue
        push!(active, t)
        local_of[t] = length(active)
    end
    reloc = Vector{Vector{Int}}(undef, length(spec.levels))
    entries = Tuple{Int,Int,Int}[]
    kmax = 0
    for l in eachindex(spec.levels)
        level = spec.levels[l]
        reloc[l] = Int[local_of[t] for t in level.re_index]
        any(iszero, reloc[l]) && return nothing
        # A level whose `dL` is empty moves nothing (no effects, or no
        # population parameters): `_continuation_member!` skips it too.
        isempty(dL[l]) && continue
        kmax = max(kmax, nrandomeffects(level))
        for t in eachindex(positions[l])
            push!(entries, (l, t, positions[l][t]))
        end
    end
    return (active=active, pop=pop, reloc=reloc, entries=entries, kmax=kmax,
        positions=positions)
end

"""
    _continuation_level_chol_second_derivatives(values, spec, positions)

`d2L_l / d values[p] d values[q]` for every level `l` and every pair of its
population parameters, indexed `[l][t][t2]` in the order
`_laplace_level_positions` gives them.

One evaluation of `_laplace_popchol` per pair, seeded along the pair's two
parameters with the seeded sweeps' own pair (`_laplace_seed_pair`): the
coefficient of the product of the two nilpotents is the mixed second
derivative, and the pure one when the two coincide. Not a nested jacobian of
`_laplace_level_chol_derivatives`'s construction: that compiles the factor
again for each count of population parameters, and on gD1 (local) took 2.95 s
of the first Hessian where this takes 0.1 s, since this type's width never
changes.
"""
function _continuation_level_chol_second_derivatives(values::Vector{Float64},
    spec::CTSEMLaplaceSpec, positions::Vector{Vector{Int}})
    e1, e2 = _laplace_seed_pair()
    E = typeof(e1)
    out = Vector{Vector{Vector{Matrix{Float64}}}}(undef, length(spec.levels))
    v = Vector{E}(undef, length(values))
    for l in eachindex(spec.levels)
        level = spec.levels[l]
        pos = positions[l]
        k = length(level.re_index)
        r = nlatent(level)
        np = length(pos)
        if np == 0 || k == 0 || r == 0
            out[l] = Vector{Matrix{Float64}}[]
            continue
        end
        out[l] = [[zeros(k, r) for _ in 1:np] for _ in 1:np]
        for t1 in 1:np, t2 in t1:np
            @inbounds for i in eachindex(values)
                v[i] = convert(E, values[i])
            end
            v[pos[t1]] += e1
            v[pos[t2]] += e2
            L = _laplace_popchol(v, level)
            D = out[l][t1][t2]
            @inbounds for q in 1:r, p in 1:k
                D[p, q] = ForwardDiff.partials(ForwardDiff.partials(L[p, q])[1])[1]
            end
            t1 == t2 || copyto!(out[l][t2][t1], D)
        end
    end
    return out
end

"""
    _continuation_hessian_scratch(laplace, S, npar, layout, depth)

The exact Hessian's buffers for this pool slot at dual type `S`, from the
pool's cached scratch as `_continuation_scratch` is, whose two it shares (the
shifted parameters, and the member's plain gradient `g`): the member's dual
input and gradient, its active Hessian `Hs`, the `dL u` columns `C` and
`V = Hs C`, and per level of the block tree the running node-weighted mean
gradient, mean Hessian and gradient scatter and the node's own gradient and
Hessian; and a unit's sum over its roots. Nothing in the node loop allocates.
"""
function _continuation_hessian_scratch(laplace::CTSEMLaplaceObjective, ::Type{S},
    npar::Integer, layout, depth::Integer) where {S}
    n = Int(npar)
    na = length(layout.active)
    ne = length(layout.entries)
    # No default arguments: they make a local function call itself, which is a
    # `Core.Box`.
    vbuf(tag, T) = _laplace_scratch_vector!(laplace, T, n, tag)
    mbuf(tag, nrow, ncol) = _laplace_scratch_matrix!(laplace, Float64, nrow, ncol, tag)
    plain = _continuation_scratch(laplace, n)
    return (x=vbuf(:cont_hess_x, S), gradient=vbuf(:cont_hess_gradient, S),
        shift=plain.shift, g=plain.grad, Hs=mbuf(:cont_hess_Hs, na, na),
        C=mbuf(:cont_hess_C, layout.kmax, ne), V=mbuf(:cont_hess_V, na, ne),
        delta=vbuf(:cont_hess_delta, Float64),
        mean=[vbuf(Symbol(:cont_hess_mean, d), Float64) for d in 1:depth],
        hess=[mbuf(Symbol(:cont_hess_hess, d), n, n) for d in 1:depth],
        scatter=[mbuf(Symbol(:cont_hess_scatter, d), n, n) for d in 1:depth],
        gnode=[vbuf(Symbol(:cont_hess_gnode, d), Float64) for d in 1:depth],
        hnode=[mbuf(Symbol(:cont_hess_hnode, d), n, n) for d in 1:depth],
        unit_gradient=vbuf(:cont_hess_unit_gradient, Float64),
        unit_hessian=mbuf(:cont_hess_unit_hessian, n, n))
end

"""
    _continuation_member_hessian!(gnode, hnode, laplace, U, m, theta, Ls, dL, d2L,
        layout, u, aws, sc)

Member `m` of unit `U` at the fixed standardised effects `u`: adds its gradient
in `theta` into `gnode` and its Hessian into `hnode`, and returns
`(loglik, status)`. Status 1 is a node the gradient path skips (the value or
the gradient not finite, as `_continuation_member!` reports it); 2 is a
Hessian that cannot be formed here: a partial not finite where the gradient
is, or a population parameter the filter reads after all.
"""
function _continuation_member_hessian!(gnode::Vector{Float64},
    hnode::Matrix{Float64}, laplace::CTSEMLaplaceObjective, U::Integer,
    m::Integer, theta::Vector{Float64}, Ls, dL, d2L, layout, u::Vector{Float64},
    aws, sc)
    S = eltype(sc.x)
    N = ForwardDiff.npartials(S)
    spec = laplace.spec
    units = laplace.units
    offsets = units.offsets[U][m]
    subject = laplace.objective.subject_objectives[units.members[U][m]]
    npar = length(theta)
    active = layout.active
    na = length(active)
    shifted = _laplace_member_values!(sc.shift, theta, spec, Ls, u, offsets)
    x = sc.x
    gd = sc.gradient
    g = sc.g
    Hs = sc.Hs
    loglik = NaN
    zero_seed = zero(ForwardDiff.Partials{N,Float64})
    for lo in 0:N:(na - 1)
        ncol = min(N, na - lo)
        @inbounds for t in 1:npar
            x[t] = S(shifted[t], zero_seed)
        end
        @inbounds for c in 1:ncol
            t = active[lo + c]
            x[t] = S(shifted[t], ForwardDiff.Partials{N,Float64}(
                ntuple(q -> q == c ? 1.0 : 0.0, Val(N))))
        end
        ll = ForwardDiff.value(_laplace_subject_value_gradient!(gd, subject, aws, x,
            nothing, _HYBRID_FULL))
        isfinite(ll) || return (ll, 1)
        if lo == 0
            loglik = ll
            @inbounds for t in 1:npar
                g[t] = ForwardDiff.value(gd[t])
            end
            all(isfinite, g) || return (NaN, 1)
            @inbounds for a in layout.pop
                iszero(g[a]) || return (NaN, 2)
            end
        end
        @inbounds for c in 1:ncol, r in 1:na
            Hs[r, lo + c] = ForwardDiff.partials(gd[active[r]], c)
        end
    end
    all(isfinite, Hs) || return (NaN, 2)
    # The gradient, J' g: what `_continuation_member!` adds.
    @inbounds for t in 1:npar
        gnode[t] += g[t]
    end
    _laplace_chol_chain!(gnode, g, spec, dL, layout.positions, u, offsets)
    # dL u for every population parameter: `J`'s extra column, on the
    # level's effects.
    C = sc.C
    entries = layout.entries
    @inbounds for (e, (l, t, a)) in enumerate(entries)
        level = spec.levels[l]
        for p in 1:nrandomeffects(level)
            C[p, e] = _laplace_dLu(dL[l][t], p, u, offsets[l], nlatent(level))
        end
    end
    # J' Hs J, block by block: active with active, then each population
    # parameter's column `V = Hs dL u` against the active parameters and, through
    # its own `dL u`, against every other population parameter.
    @inbounds for q in 1:na
        tq = active[q]
        for r in 1:na
            hnode[active[r], tq] += Hs[r, q]
        end
    end
    isempty(entries) && return (loglik, 0)
    V = sc.V
    @inbounds for (e, (l, t, a)) in enumerate(entries)
        rl = layout.reloc[l]
        for r in 1:na
            acc = 0.0
            for p in eachindex(rl)
                acc += Hs[r, rl[p]] * C[p, e]
            end
            V[r, e] = acc
            hnode[active[r], a] += acc
            hnode[a, active[r]] += acc
        end
    end
    @inbounds for (e2, (l2, t2, a2)) in enumerate(entries)
        for (e1, (l1, t1, a1)) in enumerate(entries)
            rl = layout.reloc[l1]
            acc = 0.0
            for p in eachindex(rl)
                acc += C[p, e1] * V[rl[p], e2]
            end
            # The last term: the second derivative of `L u`, within a level.
            if l1 == l2
                level = spec.levels[l1]
                rho = level.re_index
                for p in eachindex(rho)
                    acc += g[rho[p]] * _laplace_dLu(d2L[l1][t1][t2], p, u,
                        offsets[l1], nlatent(level))
                end
            end
            hnode[a1, a2] += acc
        end
    end
    return (loglik, 0)
end

"""
    _continuation_rule_hessian!(o, U, theta, Ls, dL, d2L, layout, rule, u, aws,
        sc, depth; range=1:size(rule.points, 2), finalise=true)

One block's log integral under its fixed rule at `theta`, with its gradient
left in `sc.mean[depth]` and its Hessian in `sc.hess[depth]`: the Louis
identity over the rule's nodes, the node weights normalised in one pass as
`_continuation_rule_value` does, and the mean and scatter of the node
gradients accumulated by West's weighted update, so that scores which are
large and nearly equal -- a stiff direction's, see `_continuation_stiff_rule`
-- do not cancel in a difference of two sums. Skips the nodes the gradient
path skips, and fails a block whose child cannot be scored, as it does.
Returns `(value, failed, peak, W)`, `failed` meaning the Hessian cannot be
formed.

Over the nodes `range` alone and with `finalise = false`, the parts a root's
nodes summed on different threads are merged from
(`_continuation_merge_moments!`): `peak`, the weights' sum `W` relative to it,
and in `sc` the weighted mean gradient, the weighted mean of the node Hessians
and the gradients' scatter, the last not yet divided by `W` and added.

The running-maximum rescale is written out here as it is in
`_continuation_rule_value` and `_quadrature_block`, which carry one moment
and none: an accumulator shared across the three would need a closure or a
buffer per moment in a node loop that allocates nothing.
"""
function _continuation_rule_hessian!(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, d2L, layout, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, sc, depth::Int;
    range::UnitRange{Int}=1:size(rule.points, 2), finalise::Bool=true)
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    npar = length(theta)
    n = size(rule.points, 2)
    leaf = isempty(rule.children)
    mean = sc.mean[depth]
    hess = sc.hess[depth]
    scatter = sc.scatter[depth]
    gnode = sc.gnode[depth]
    hnode = sc.hnode[depth]
    delta = sc.delta
    fill!(mean, 0.0)
    fill!(hess, 0.0)
    fill!(scatter, 0.0)
    peak = -Inf
    W = 0.0
    for j in range
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        fill!(gnode, 0.0)
        fill!(hnode, 0.0)
        h = 0.0
        ok = true
        if leaf
            for m in block.members
                ll, status = _continuation_member_hessian!(gnode, hnode, laplace, U,
                    m, theta, Ls, dL, d2L, layout, u, aws, sc)
                status == 2 && return (value=NaN, failed=true, peak=NaN, W=0.0)
                if status == 1 || !isfinite(ll)
                    ok = false
                    break
                end
                h += ll
            end
        else
            for child in rule.children[j]
                r = _continuation_rule_hessian!(o, U, theta, Ls, dL, d2L, layout,
                    child, u, aws, sc, depth + 1)
                r.failed && return r
                isfinite(r.value) || return (value=NaN, failed=false, peak=NaN, W=0.0)
                h += r.value
                gnode .+= sc.mean[depth + 1]
                hnode .+= sc.hess[depth + 1]
            end
        end
        ok || continue
        t = rule.logweights[j] + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                W *= s
                scatter .*= s
            end
            peak = t
        end
        w = exp(t - peak)
        Wn = W + w
        rho = w / Wn
        f = w * (1 - rho)
        @inbounds for a in 1:npar
            delta[a] = gnode[a] - mean[a]
        end
        @inbounds for b in 1:npar
            db = delta[b]
            for a in 1:npar
                scatter[a, b] += f * delta[a] * db
                hess[a, b] += rho * (hnode[a, b] - hess[a, b])
            end
            mean[b] += rho * db
        end
        W = Wn
    end
    isfinite(peak) || return (value=NaN, failed=false, peak=peak, W=W)
    if finalise
        invW = 1 / W
        @inbounds for b in 1:npar, a in 1:npar
            hess[a, b] += scatter[a, b] * invW
        end
    end
    return (value=peak + log(W) + rule.constant, failed=false, peak=peak, W=W)
end

"""
    _continuation_merge_moments!(mean, hess, parts)

The Louis identity's moments over a root's nodes from those of its parts, each
`(peak, W, mean, hess, scatter)` as `_continuation_rule_hessian!` leaves them
unfinalised over a range of the nodes: the parts' weights brought to the
largest peak, their means and mean Hessians averaged by them, and the scatter
about the whole mean -- each part's own plus its mean's distance from the
whole's (Chan, Golub and LeVeque 1979) -- divided by the total weight and
added to the mean Hessian. Writes `mean` and `hess`, and returns `(peak, W)`.
"""
function _continuation_merge_moments!(mean::Vector{Float64}, hess::Matrix{Float64},
    parts)
    P = maximum(p.peak for p in parts)
    W = sum(p.W * exp(p.peak - P) for p in parts)
    fill!(mean, 0.0)
    fill!(hess, 0.0)
    for p in parts
        r = p.W * exp(p.peak - P) / W
        mean .+= r .* p.mean
        hess .+= r .* p.hess
    end
    for p in parts
        e = exp(p.peak - P)
        r = p.W * e / W
        d = p.mean .- mean
        hess .+= (e / W) .* p.scatter .+ r .* (d .* transpose(d))
    end
    return (peak=P, W=W)
end

# Width `wide` has been run on a model of this shape in this session, so a
# sweep at it compiles nothing (`_continuation_hessian_width`). One entry per
# model shape and width a session runs; like the transform cache, never
# emptied.
const _CONTINUATION_WIDE_SEEN = Set{Any}()
const _CONTINUATION_WIDE_LOCK = ReentrantLock()
# Seconds of width-one sweeps above which a wider dual pays for its compile.
const _CONTINUATION_WIDE_AFTER = Ref(30.0)

"""
    _continuation_hessian_width(o, theta, layout)

How many directions one dual sweep of the exact Hessian carries.

A sweep's cost grows much more slowly than its width: on the ordinal fixture
of the bench cell ord4 (two latents, 30 rows, 17 active parameters; local
machine, one thread) a plain sweep took 0.79 ms and a dual one 1.27 ms at
width 1, 1.39 ms at 8, 1.87 ms at 12 and 2.81 ms at 17. But a width the model
has not been swept at compiles the whole filter and its reverse pass again for
that type, which there took 18 s at width 2, 26 s at 8 and 45 s at 17, while
width one is the type the Laplace fit's own seeded sweeps compiled. So the
width is ForwardDiff's chunk for the active count where one has already run
on a model of this shape in this session, or where width one is forecast to
take more than `_CONTINUATION_WIDE_AFTER[]` seconds (one timed sweep, times
the sweeps the flagged rules need, over the workers), and one otherwise.
"""
function _continuation_hessian_width(o::CTSEMLaplaceContinuation,
    theta::Vector{Float64}, layout)
    na = length(layout.active)
    wide = ForwardDiff.pickchunksize(max(na, 1))
    (wide <= 1 || isempty(o.rules)) && return 1
    key = (typeof(o.laplace.objective.params), wide)
    lock(() -> key in _CONTINUATION_WIDE_SEEN, _CONTINUATION_WIDE_LOCK) && return wide
    laplace = o.laplace
    npar = length(theta)
    S = ForwardDiff.Dual{_LaplaceSeedInner,Float64,1}
    subject = laplace.objective.subject_objectives[laplace.units.members[o.rules[1].unit][1]]
    aws = _laplace_workspace!(laplace, S, npar)
    x = S[S(theta[t], ForwardDiff.Partials((t == layout.active[1] ? 1.0 : 0.0,)))
          for t in 1:npar]
    gd = similar(x)
    best = Inf
    for _ in 1:3
        started = time_ns()
        _laplace_subject_value_gradient!(gd, subject, aws, x, nothing, _HYBRID_FULL)
        best = min(best, (time_ns() - started) / 1e9)
    end
    sweeps = sum(r.evaluations for r in o.rules; init=0)
    workers = min(_laplace_pool_width(), length(o.rules))
    return best * na * sweeps / workers > _CONTINUATION_WIDE_AFTER[] ? wide : 1
end

"""
    _continuation_flagged_hessian(o, theta; width=0)

The flagged units' terms at `theta`, their summed gradient and their summed
Hessian, exactly, in parallel over units; `nothing` when the Hessian cannot be
formed this way (`_continuation_member_hessian!`'s status 2, or a unit that
cannot be scored), so that the caller differences instead.
"""
function _continuation_flagged_hessian(o::CTSEMLaplaceContinuation,
    theta::Vector{Float64}; width::Integer=0)
    laplace = o.laplace
    spec = laplace.spec
    npar = length(theta)
    isempty(o.rules) && return (values=Float64[], gradient=zeros(npar),
        hessian=zeros(npar, npar), width=0)
    Ls = _laplace_popchols(theta, spec)
    dL = _laplace_level_chol_derivatives(theta, spec)
    positions = Vector{Int}[_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    hyb = laplace.objective.params.hybrid
    layout = _continuation_hessian_layout(spec, npar, positions, dL,
        _hybrid_on(hyb) ? _hybrid_filtered(hyb) : Int[])
    (layout === nothing || isempty(layout.active)) && return nothing
    d2L = _continuation_level_chol_second_derivatives(theta, spec, positions)
    _laplace_ensure_pool!(laplace)
    _continuation_warm_caches(laplace, o.nodes, o.soft_maxdirs, o.product_maxdim)
    N = width > 0 ? Int(width) : _continuation_hessian_width(o, theta, layout)
    out = _continuation_flagged_hessian(o, theta, Ls, dL, d2L, layout, Val(N))
    if out !== nothing && N > 1
        key = (typeof(laplace.objective.params), N)
        lock(() -> push!(_CONTINUATION_WIDE_SEEN, key), _CONTINUATION_WIDE_LOCK)
    end
    return out
end

function _continuation_flagged_hessian(o::CTSEMLaplaceContinuation,
    theta::Vector{Float64}, Ls, dL, d2L, layout, ::Val{N}) where {N}
    laplace = o.laplace
    S = ForwardDiff.Dual{_LaplaceSeedInner,Float64,N}
    nf = length(o.rules)
    npar = length(theta)
    depth = length(laplace.spec.levels)
    nslot = _laplace_pool_width()
    gparts = [zeros(npar) for _ in 1:nslot]
    hparts = [zeros(npar, npar) for _ in 1:nslot]
    values = fill(NaN, nf)
    good = fill(true, nslot)
    # A unit whole on one slot, or, where there are fewer units than slots,
    # its root's nodes in ranges (`_continuation_tasks`), each range's moments
    # kept by the task and merged below.
    tasks = _continuation_tasks(o, nslot)
    ntask = length(tasks)
    parts = Vector{Any}(nothing, ntask)
    _laplace_parallel(laplace, 1:ntask) do q
        # `local` throughout, as in `_continuation_flagged`.
        local f, range, rule, U, aws, sc, u, total, failed, r, slot
        slot = _laplace_slot()
        f, range = tasks[q]
        rule = o.rules[f]
        U = rule.unit
        aws = _laplace_workspace!(laplace, S, npar)
        sc = _continuation_hessian_scratch(laplace, S, npar, layout, depth)
        u = zeros(laplace.units.dims[U])
        failed = false
        if range !== nothing
            try
                r = _continuation_rule_hessian!(o, U, theta, Ls, dL, d2L, layout,
                    rule.roots[1], u, aws, sc, 1; range=range, finalise=false)
                failed = r.failed || isnan(r.peak)
                failed || (parts[q] = (peak=r.peak, W=r.W, mean=copy(sc.mean[1]),
                    hess=copy(sc.hess[1]), scatter=copy(sc.scatter[1])))
            catch err
                _ctsem_must_propagate(err) && rethrow()
                failed = true
            end
            failed && (good[slot] = false)
            return !failed
        end
        fill!(sc.unit_gradient, 0.0)
        fill!(sc.unit_hessian, 0.0)
        total = 0.0
        try
            for root in rule.roots
                r = _continuation_rule_hessian!(o, U, theta, Ls, dL, d2L, layout,
                    root, u, aws, sc, 1)
                if r.failed || !isfinite(r.value)
                    failed = true
                    break
                end
                total += r.value
                sc.unit_gradient .+= sc.mean[1]
                sc.unit_hessian .+= sc.hess[1]
            end
        catch err
            _ctsem_must_propagate(err) && rethrow()
            failed = true
        end
        if failed || !all(isfinite, sc.unit_gradient) || !all(isfinite, sc.unit_hessian)
            good[slot] = false
            return false
        end
        values[f] = total
        w = _laplace_unit_weight(laplace, U)
        gparts[slot] .+= w .* sc.unit_gradient
        hparts[slot] .+= w .* sc.unit_hessian
        return true
    end
    all(good) || return nothing
    gradient = zeros(npar)
    hessian = zeros(npar, npar)
    for s in 1:nslot
        gradient .+= gparts[s]
        hessian .+= hparts[s]
    end
    mean = zeros(npar)
    hess = zeros(npar, npar)
    at = 1
    while at <= ntask
        f, range = tasks[at]
        last = at
        while last < ntask && tasks[last + 1][1] == f
            last += 1
        end
        if range !== nothing
            pieces = [parts[q] for q in at:last if isfinite(parts[q].peak)]
            isempty(pieces) && return nothing
            m = _continuation_merge_moments!(mean, hess, pieces)
            (all(isfinite, mean) && all(isfinite, hess)) || return nothing
            rule = o.rules[f]
            values[f] = m.peak + log(m.W) + rule.roots[1].constant
            w = _laplace_unit_weight(laplace, rule.unit)
            gradient .+= w .* mean
            hessian .+= w .* hess
        end
        at = last + 1
    end
    # Dual sweeps, counted as the reverse sweeps they are.
    o.member_sweeps += sum(r.evaluations for r in o.rules; init=0) *
        cld(length(layout.active), N)
    return (values=values, gradient=gradient, hessian=hessian, width=N)
end

"""
    _continuation_exact_hessian(o, x; step=1e-4, width=0)

The exact scheme of `ctsem_laplace_continuation_hessian`: the flagged units'
part by `_continuation_flagged_hessian`, and the part the other units and the
prior carry by `ctsem_laplace_hessian` on the subset objective that holds
them, whose inner modes are put back as they were. `nothing` when the flagged
part cannot be formed exactly.
"""
function _continuation_exact_hessian(o::CTSEMLaplaceContinuation, x::Vector{Float64};
    step::Real=1e-4, width::Integer=0)
    flagged = _continuation_flagged_hessian(o, x; width=width)
    flagged === nothing && return nothing
    saved = [copy(m) for m in o.rest.modes]
    rest = try
        ctsem_laplace_hessian(o.rest, x; step=step)
    finally
        for U in eachindex(saved)
            o.rest.modes[U] = saved[U]
        end
    end
    H = rest .+ flagged.hessian
    return (H .+ transpose(H)) ./ 2
end

"""
    _continuation_error(o)

The hybrid's own error against the exact marginal at its placement, in nats:
the frozen rules' standard errors in quadrature (`draws`), the product rules'
node errors (`nodes`), what the screen leaves on Laplace (`residual`), and
their sum (`total`). What `ctsem_laplace_continuation_info` reports, and what
a probe's value is read net of (`_ctsem_probe_value`).
"""
function _continuation_error(o::CTSEMLaplaceContinuation)
    gaps = o.centre_quadrature .- o.centre_laplace
    residual = _continuation_screen(gaps[o.rest_units], o.centre_se[o.rest_units])
    draws = sqrt(sum(abs2(o.centre_se[U]) for U in o.flagged if o.importance[U];
        init=0.0))
    nodes = sum(o.centre_error[U] for U in o.flagged
        if !o.importance[U] && isfinite(o.centre_error[U]); init=0.0)
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
    softblocks = sum(sum(_continuation_soft_blocks(r) for r in rule.roots; init=0)
        for rule in o.rules; init=0)
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
        rule_failures=count(o.rule_failed), soft_blocks=softblocks,
        nwide=count(o.wide), maxdim=o.maxdim,
        importance=[o.importance[U] for U in 1:nunits], nimportance=count(o.importance),
        centre_se=copy(o.centre_se),
        importance_draws=sum(sum(_continuation_draws(root) for root in r.roots;
            init=0) for r in o.rules if o.importance[r.unit]; init=0),
        draws=o.draws, draws_se=o.draws_se,
        evaluations_per_value=sum(r.evaluations for r in o.rules; init=0),
        value_calls=o.value_calls, gradient_calls=o.gradient_calls,
        member_values=o.member_values, member_sweeps=o.member_sweeps,
        recentres=o.recentres, refused=o.refused,
        seconds_values=o.seconds_values, seconds_gradients=o.seconds_gradients,
        seconds_placements=o.seconds_placements, seconds_hessian=o.seconds_hessian,
        nodes=o.nodes,
        product_maxdim=o.product_maxdim, soft_tau=o.soft_tau,
        soft_maxdirs=o.soft_maxdirs, tolerance=o.tolerance)
end

"""
    _ctsem_jump_units(objective, jump)

The units whose terms jumped along the step the line search saw jump
(`CTSEMJump`), with their changes and where each one's inner solve stopped
either side, as `(units, changes, point_units, point_values)`; `nothing` on a
route with no units to name. On the Laplace route each unit's term is taken at
`x`, `x + alpha_lo s` and `x + alpha_hi s`: a unit whose change over the short
step keeps a quarter of its change over the long one, a hundred times longer,
and carries a tenth of the jump, moved with where its inner solve stopped
rather than with the parameters (`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`).
Its modes at the three points are what a frozen importance rule for it starts
from (`ctsem_laplace_continuation`'s `point_units`): at one point the solve
returns one of them, and the rule needs both sides.
"""
_ctsem_jump_units(::Any, jump) = nothing
_ctsem_jump_capable(::Union{CTSEMLaplaceObjective,CTSEMLaplaceContinuation}) = true

function _ctsem_jump_units(laplace::CTSEMLaplaceObjective, jump::CTSEMJump)
    return _laplace_jump_units(laplace, collect(1:length(laplace.units.members)),
        jump)
end

# On the hybrid only its Laplace units can jump; the flagged ones are frozen
# rules, smooth by construction. Named in the full objective's numbering.
_ctsem_jump_units(o::CTSEMLaplaceContinuation, jump::CTSEMJump) =
    _laplace_jump_units(o.rest, o.rest_units, jump)

function _laplace_jump_units(laplace::CTSEMLaplaceObjective, numbering::Vector{Int},
    jump::CTSEMJump)
    thetas = [jump.x, jump.x .+ jump.alpha_lo .* jump.s,
        jump.x .+ jump.alpha_hi .* jump.s]
    terms = Vector{Vector{Float64}}()
    modes = Vector{Vector{Vector{Float64}}}()
    for th in thetas
        r = ctsem_laplace_evaluate(laplace, th; gradient=false)
        push!(terms, copy(r.unit_loglik))
        push!(modes, [copy(m) for m in laplace.modes])
    end
    nunits = length(numbering)
    units = Int[]
    changes = Float64[]
    point_units = Int[]
    columns = Vector{Vector{Float64}}()
    for p in 1:nunits
        w = _laplace_unit_weight(laplace, p)
        lo = w * (terms[2][p] - terms[1][p])
        hi = w * (terms[3][p] - terms[1][p])
        (isfinite(lo) && isfinite(hi)) || continue
        abs(lo) >= jump.excess / 10 || continue
        abs(lo) >= abs(hi) / 4 || continue
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
    values = zeros(width, max(length(columns), 1))
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
as the parameters move away from it their weights spread; a fit rescued from a
jump places them again only once this has drifted past what they were placed
to (`.ctLaplaceAnytimeRecord()`, the correction's rounds). In the order of
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
