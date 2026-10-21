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
than `maxdim` and of one block gets a frozen importance rule instead when
`draws > 0` (below); a nested one is not scored and keeps its Laplace term
(`wide` in the info): the step correction, the alternative, pays the same
`nodes^d` per quadrature value, so falling back to it would save nothing.

# Units too wide for a product rule: frozen importance draws

A unit of one block wider than `maxdim` gets draws in place of nodes
(`_continuation_importance_rule`): held fixed in the standardised effects, as
nodes are, from a defensive mixture -- a Gaussian at the unit's mode, at any
other point the caller knows for it (`points`: where the inner solve stopped
at nearby parameters), at draws that carried weight once moved a Newton step,
and the prior -- and weighted by that mixture. Its cost is the number of draws
whatever the width, and the same Fisher identity gives its gradient and the
same Louis identity its Hessian. Why: a unit whose posterior has a broad
shoulder beside its mode has a Laplace term that moves with where the inner
solve stops, by up to the shoulder's mass, and an optimiser climbing onto the
point where it stops on the shoulder stalls there with a gradient that is that
spike's (`CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md`: 0.33 nats
on an 8-effect person, the fit stalled for hours). A frozen rule does not read
the stopping point, so its term is smooth.

Its value carries Monte Carlo error, `sqrt(1/ESS - 1/n)` estimated from its
weights, which no number of draws a fit can afford takes below the screen's
tolerance. So a unit's gap counts towards flagging only beyond twice that
error (`_continuation_place!`): the screen flags what it can resolve, and a
unit Laplace gets right to within the draws' error keeps its Laplace term. A
unit the caller names (`force`) is flagged whatever its gap: the optimiser
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
is cost: a gradient is one reverse sweep per member per draw.
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
a product rule.
"""
struct CTSEMFixedRule
    block::Int
    points::Matrix{Float64}
    logweights::Vector{Float64}
    constant::Float64
    children::Vector{Vector{CTSEMFixedRule}}
    soft::Int
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

_continuation_soft_blocks(rule::CTSEMFixedRule) = (rule.soft > 0 ? 1 : 0) +
    sum(sum(_continuation_soft_blocks(c) for c in kids; init=0)
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
    # Units wider than `maxdim` that no importance rule can take (nested, or
    # `draws = 0`): never scored, never flagged.
    wide::Vector{Bool}
    # The importance rule's options (see the module docstring), the points the
    # caller knows for each unit, and per unit at the centre: the rule's
    # standard error (0 for a product rule) and whether it is an importance rule.
    draws::Int
    draws_se::Float64
    draws_max::Int
    draws_floor::Float64
    draws_prior::Float64
    draws_seed::Int
    points::Dict{Int,Vector{Vector{Float64}}}
    centre_se::Vector{Float64}
    importance::Vector{Bool}
    # Whether a placement scores every unit, or only the flagged ones (a fit
    # rescued from a jump, which needs the units the optimiser named and no
    # screen).
    screen::Bool
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
A Gaussian component of an importance rule's proposal in a block's own
standardised coordinates: centre `mu`, precision `vecs * Diagonal(lam) * vecs'`
-- the block's precision at the centre, its eigenvalues floored.
"""
struct CTSEMImportanceComponent
    mu::Vector{Float64}
    vecs::Matrix{Float64}
    lam::Vector{Float64}
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
    return -(q + k * log(2 * pi)) / 2 + sum(log, c.lam) / 2
end

"""
    _continuation_importance_component(fns, x, floor)

The component at block coordinates `x`: the block's precision there
(`fns.precision`) with its eigenvalues floored at `floor`, so a soft or an
indefinite direction is sampled at a bounded width rather than an infinite
one. `nothing` where the precision is not finite.
"""
function _continuation_importance_component(fns, x::Vector{Float64}, floor::Float64)
    P = fns.precision(x)
    all(isfinite, P) || return nothing
    E = _ctsem_symeig(_laplace_symmetrise(P))
    all(isfinite, E.values) || return nothing
    return CTSEMImportanceComponent(copy(x), Matrix(E.vectors), max.(E.values, floor))
end

"""Add the component at `x` unless one sits within 1e-3 of it; whether it was
added."""
function _continuation_add_component!(comps::Vector{CTSEMImportanceComponent}, fns,
    x::Vector{Float64}, floor::Float64)
    any(c -> maximum(abs, c.mu .- x) < 1e-3, comps) && return false
    c = _continuation_importance_component(fns, x, floor)
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
    c = _continuation_importance_component(fns, x, floor)
    c === nothing && return nothing
    w = transpose(c.vecs) * vg.gradient
    len = sqrt(sum(abs2, w ./ sqrt.(c.lam)))
    scale = len > maxstep ? maxstep / len : 1.0
    return x .+ c.vecs * (scale .* w ./ c.lam)
end

"""The draws an importance rule has taken: their block coordinates, the
component each came from (0 the prior), its members' log likelihood, and how
many batches have been drawn."""
mutable struct _ContinuationDraws
    points::Matrix{Float64}
    origin::Vector{Int}
    h::Vector{Float64}
    batches::Int
end

"""`m` split as evenly as whole numbers allow over `J` parts, the larger
first."""
function _continuation_split(m::Integer, J::Integer)
    J <= 0 && return Int[]
    q, r = divrem(Int(m), Int(J))
    return [q + (t <= r ? 1 : 0) for t in 1:Int(J)]
end

"""`m` draws from component `component` (`c`), scored at `theta`: each draw's
members' log likelihood, the deviations their filters integrate integrated
exactly, as at a product rule's node. The lattice's shift comes from `key`, the
batch and the component alone, so the same component of the same batch draws
the same standard normals at every placement: common random numbers."""
function _continuation_draw!(D::_ContinuationDraws, laplace::CTSEMLaplaceObjective,
    U::Integer, theta::Vector{Float64}, Ls, block, columns, u::Vector{Float64}, aws,
    c::CTSEMImportanceComponent, component::Int, m::Int, key)
    m > 0 || return D
    k = length(c.mu)
    Z = _continuation_lattice_normals(k, m,
        Random.Xoshiro(hash((key, D.batches, component))))
    X = Matrix{Float64}(undef, k, m)
    hs = Vector{Float64}(undef, m)
    for j in 1:m
        @inbounds for t in 1:k
            s = 0.0
            for a in 1:k
                s += c.vecs[t, a] * Z[a, j] / sqrt(c.lam[a])
            end
            X[t, j] = c.mu[t] + s
        end
        @inbounds for (t, col) in enumerate(columns); u[col] = X[t, j]; end
        inner = 0.0
        for mm in block.members
            shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                laplace.units.offsets[U][mm])
            inner += _laplace_member_loglik(aws, laplace.objective.subject_objectives[
                laplace.units.members[U][mm]], shifted)
        end
        hs[j] = isfinite(inner) ? inner : -Inf
    end
    D.points = hcat(D.points, X)
    append!(D.origin, fill(component, m))
    append!(D.h, hs)
    return D
end

"""A batch of `m` draws: the prior's share, the rest split over the
components."""
function _continuation_batch!(D::_ContinuationDraws, laplace, U, theta, Ls, block,
    columns, u, aws, comps, prior::CTSEMImportanceComponent, share::Float64,
    m::Int, key)
    D.batches += 1
    mp = clamp(round(Int, share * m), 1, m)
    _continuation_draw!(D, laplace, U, theta, Ls, block, columns, u, aws, prior, 0,
        mp, key)
    for (t, n) in enumerate(_continuation_split(m - mp, length(comps)))
        _continuation_draw!(D, laplace, U, theta, Ls, block, columns, u, aws,
            comps[t], t, n, key)
    end
    return D
end

"""
The log weights of every draw under the balance heuristic -- each component
weighted by its share of the draws -- and the log integral, its standard error
`sqrt(1/ESS - 1/n)` and effective sample size, for a block of `k` effects.
"""
function _continuation_importance_weigh(D::_ContinuationDraws,
    comps::Vector{CTSEMImportanceComponent}, prior::CTSEMImportanceComponent, k::Int)
    n = size(D.points, 2)
    J = length(comps)
    counts = zeros(Int, J + 1)
    for t in D.origin
        counts[t + 1] += 1
    end
    lw = Vector{Float64}(undef, n)
    dens = Vector{Float64}(undef, J + 1)
    work = zeros(k)
    @inbounds for j in 1:n
        z = view(D.points, :, j)
        for t in 0:J
            dens[t + 1] = counts[t + 1] == 0 ? -Inf : log(counts[t + 1] / n) +
                _continuation_component_logdensity(t == 0 ? prior : comps[t], z, work)
        end
        lw[j] = -dot(z, z) / 2 - _logsumexp(dens) - log(n)
    end
    terms = lw .+ D.h
    peak = maximum(terms)
    isfinite(peak) || return (lw=lw, terms=terms, value=NaN, se=Inf, ess=0.0)
    w = exp.(terms .- peak)
    sw = sum(w)
    ess = sw^2 / sum(abs2, w)
    return (lw=lw, terms=terms, value=peak + log(sw) - k * log(2 * pi) / 2,
        se=sqrt(max(1 / ess - 1 / n, 0.0)), ess=ess)
end

"""
    _continuation_importance_rule(laplace, U, theta, Ls, b, u, context, aws, opts)

Block `b`'s frozen importance rule (see the module docstring): draws held fixed
in the standardised effects, weighted by the mixture they came from, so the
block's term is

    log (1/n) sum_j L(z_j) exp(-z_j'z_j/2) / q(z_j) - k/2 log 2 pi

-- a `CTSEMFixedRule` whose log weights carry `-z'z/2 - log q - log n`, valued
and differentiated as a product rule's nodes are. `q` is the balance heuristic
over its components, each weighted by its share of the draws:

  - a Gaussian at the block's mode (`context.mode`) and at every point
    `opts.points` holds for the unit, each on the block's precision there with
    its eigenvalues floored at `opts.draws_floor`;
  - the prior `N(0, I)`, `opts.draws_prior` of every batch, which bounds every
    weight by the likelihood's supremum over that share, so the estimate has a
    finite variance in any dimension (Hesterberg 1995);
  - where the weights show a gap -- a draw carrying more than `n / ESS*` times
    its share, `ESS* = 1 / draws_se^2` -- a component at each of up to three
    such draws, moved a damped Newton step and shaped by the precision there
    (`_continuation_newton_point`).

Draws come in batches, the first of `opts.draws` and each later as large as
all before it, quasi-random (`_continuation_lattice_normals`), until the
standard error is under `opts.draws_se` or `opts.draws_max` are spent. Each
component's lattice in each batch is seeded by the unit, block, batch and
component, the same at every placement: common random numbers, so two
placements differ by where their proposals sit rather than by their draws'
luck. With fresh draws at each placement the continuation's rounds, which keep
a step only if the objective with the rule placed again is no lower, were
rejected by the draws' error on the round-9 fit: 0.4 to 0.7 nats at radii
down to 0.03 standard errors. Returns `(rule, value, se, ess,
components)`, or `nothing` when no draw scores.
"""
function _continuation_importance_rule(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, b::Integer,
    u::Vector{Float64}, context, aws, opts)
    block = laplace.units.blocks[U][b]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    fns = _continuation_block_functions(laplace, U, theta, Ls, b, u, aws)
    floor = Float64(opts.draws_floor)
    comps = CTSEMImportanceComponent[]
    _continuation_add_component!(comps, fns, Float64[context.mode[c] for c in columns],
        floor)
    for p in get(opts.points, Int(U), Vector{Vector{Float64}}())
        length(p) == length(u) || continue
        x = Float64[p[c] for c in columns]
        all(isfinite, x) && _continuation_add_component!(comps, fns, x, floor)
    end
    isempty(comps) && return nothing
    prior = CTSEMImportanceComponent(zeros(k), Matrix{Float64}(I, k, k), ones(k))
    key = (Int(opts.draws_seed), Int(U), Int(b))
    D = _ContinuationDraws(zeros(k, 0), Int[], Float64[], 0)
    share = Float64(opts.draws_prior)
    _continuation_batch!(D, laplace, U, theta, Ls, block, columns, u, aws, comps,
        prior, share, Int(opts.draws), key)
    current = _continuation_importance_weigh(D, comps, prior, k)
    target_ess = 1 / Float64(opts.draws_se)^2
    while (!isfinite(current.value) || current.se > opts.draws_se) &&
            size(D.points, 2) < opts.draws_max
        n = size(D.points, 2)
        if isfinite(current.value)
            peak = maximum(current.terms)
            w = exp.(current.terms .- peak)
            ratio = w ./ (sum(w) / n)
            added = 0
            for j in sortperm(ratio; rev=true)
                (added >= 3 || ratio[j] <= n / target_ess) && break
                x = _continuation_newton_point(fns, D.points[:, j], floor)
                x === nothing && continue
                _continuation_add_component!(comps, fns, x, floor) && (added += 1)
            end
        end
        _continuation_batch!(D, laplace, U, theta, Ls, block, columns, u, aws, comps,
            prior, share, min(n, Int(opts.draws_max) - n), key)
        current = _continuation_importance_weigh(D, comps, prior, k)
    end
    isfinite(current.value) || return nothing
    # The block's own coordinates back where the caller had them.
    @inbounds for c in columns; u[c] = context.mode[c]; end
    rule = CTSEMFixedRule(b, D.points, current.lw, -k * log(2 * pi) / 2,
        Vector{Vector{CTSEMFixedRule}}(), 0)
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
    children = Vector{Vector{CTSEMFixedRule}}()
    terms = Vector{Float64}(undef, n)
    for j in 1:n
        @inbounds for (t, c) in enumerate(columns); u[c] = points[t, j]; end
        inner = 0.0
        if leaf
            for m in block.members
                shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                    laplace.units.offsets[U][m])
                inner += _laplace_member_loglik(aws, laplace.objective.subject_objectives[
                    laplace.units.members[U][m]], shifted)
            end
        else
            kidrules = CTSEMFixedRule[]
            for c in kids
                placed = _continuation_block_rule(laplace, U, theta, Ls, c, u,
                    context, aws, opts)
                placed === nothing && return nothing
                push!(kidrules, placed.rule)
                inner += placed.value
            end
            push!(children, kidrules)
        end
        terms[j] = isfinite(inner) ? inner + logweights[j] : -Inf
    end
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
    if _continuation_unit_width(laplace, U) > opts.maxdim
        # Too wide for a product rule: only `_continuation_place` sends such a
        # unit here, and only one of one block when draws are on.
        tree = _quadrature_children(blocks)
        placed = _continuation_importance_rule(laplace, U, theta, Ls, tree.roots[1],
            copy(mode), (children=tree.children, mode=mode), aws, opts)
        placed === nothing && return nothing
        evaluations = size(placed.rule.points, 2) * length(blocks[1].members)
        return (rule=CTSEMFixedUnitRule(U, [placed.rule], evaluations),
            value=placed.value, se=placed.se, importance=true)
    end
    M = _laplace_unit_curvature(laplace, U, theta, Ls, mode)
    factorization = _laplace_factor_repaired!(M, blocks)
    factorization.ok || return nothing
    tree = _quadrature_children(blocks)
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
    # Wide and out of reach: a unit too wide for a product rule that is nested,
    # or any such unit when draws are off.
    wide = [_continuation_unit_width(laplace, U) > opts.maxdim &&
        !(opts.draws > 0 && length(laplace.units.blocks[U]) == 1) for U in 1:nunits]
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
    _continuation_resolved(gaps, se)

Each unit's gap beyond twice its rule's standard error, towards zero: what the
screen can tell from the draws' own error (see the module docstring). The same
gap where the rule is a product rule, whose `se` is zero.
"""
_continuation_resolved(gaps::Vector{Float64}, se::Vector{Float64}) =
    [isfinite(g) ? sign(g) * max(abs(g) - 2 * se[U], 0.0) : g
     for (U, g) in enumerate(gaps)]

"""
    _continuation_flags(gaps, tolerance)

The units the continuation puts nodes on: largest `|gap|` first, the smallest
set whose removal leaves the rest carrying at most `tolerance` in total. That
is the tolerance the screen passes a whole fit on, so what the hybrid leaves to
Laplace is what the screen would have called exact. A unit whose rule could
not be placed (`NaN`) is never flagged. A negative `tolerance` flags every
other unit, a zero gap included: at `0` a unit whose gap is exactly zero, as on
a linear Gaussian model where Laplace is exact and the gaps are rounding, is
left to Laplace, so whether it is flagged would be a matter of the last digit.
"""
function _continuation_flags(gaps::Vector{Float64}, tolerance::Real)
    g = [isfinite(x) ? abs(x) : 0.0 for x in gaps]
    order = sortperm(g; rev=true)
    # What the unflagged units carry, as a sum of them rather than the total
    # less what has been flagged: subtracting the largest gaps from their total
    # can leave zero or less while a small one is still unflagged, which at
    # `tolerance = 0` dropped a unit on rounding alone.
    remaining = zeros(length(order) + 1)
    for k in length(order):-1:1
        remaining[k] = remaining[k + 1] + g[order[k]]
    end
    flagged = Int[]
    for (k, U) in enumerate(order)
        remaining[k] <= tolerance && break
        (g[U] > 0 || tolerance < 0) || break
        push!(flagged, U)
    end
    return sort!(flagged)
end

"""
    ctsem_laplace_continuation(laplace, values; nodes=5, tolerance=0.01,
        product_maxdim=2, soft_tau=3.5, soft_maxdirs=2, maxdim=5, draws=0,
        draws_se=0.05, draws_max=16384, draws_floor=0.25, draws_prior=0.2,
        draws_seed=20261021, force=nothing, point_units=nothing,
        point_values=nothing, screen=true)

The hybrid objective with its nodes placed at `values`. Every unit is scored by
the rule there against its Laplace term, the units that differ beyond the
rule's own error are flagged (`_continuation_resolved`, `_continuation_flags`),
and those carry fixed nodes from then on while the rest keep their Laplace
term; the units in `force` are flagged whatever their gap. A unit wider than
`maxdim` effects along any path of its block tree gets a frozen importance
rule when it is one block and `draws > 0`, and is not scored at all
otherwise, keeping its Laplace term. `point_units` and `point_values` are
points the caller knows for some units, column `j` of `point_values` for unit
`point_units[j]` (its first `dims` entries), each a component of that unit's
importance proposal. `screen = false` scores only the units `force` names, the
rest keeping their Laplace terms unscored: what a fit rescued from a jump
needs. See the module docstring for the rules and the constants.
"""
function ctsem_laplace_continuation(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, tolerance::Real=0.01,
    product_maxdim::Integer=2, soft_tau::Real=3.5, soft_maxdirs::Integer=2,
    maxdim::Integer=5, draws::Integer=0, draws_se::Real=0.05,
    draws_max::Integer=16384, draws_floor::Real=0.25, draws_prior::Real=0.2,
    draws_seed::Integer=20261021, force=nothing, point_units=nothing,
    point_values=nothing, screen::Bool=true)
    nodes >= 1 || throw(ArgumentError("need at least one quadrature node"))
    draws >= 0 || throw(ArgumentError("draws must be at least zero"))
    draws_se > 0 || throw(ArgumentError("draws_se must be positive"))
    0 < draws_prior < 1 || throw(ArgumentError("draws_prior must be between 0 and 1"))
    draws_floor > 0 || throw(ArgumentError("draws_floor must be positive"))
    theta = _laplace_extend(laplace, collect(Float64, values))
    _laplace_check_indices(laplace, length(theta))
    nunits = length(laplace.units.members)
    o = CTSEMLaplaceContinuation{typeof(laplace)}(laplace, laplace, collect(1:nunits),
        Int[], CTSEMFixedUnitRule[], theta, Int(nodes), Float64(tolerance),
        Int(product_maxdim), Float64(soft_tau), Int(soft_maxdirs), Int(maxdim),
        fill(false, nunits), Int(draws), Float64(draws_se),
        max(Int(draws_max), Int(draws)), Float64(draws_floor), Float64(draws_prior),
        Int(draws_seed), Dict{Int,Vector{Vector{Float64}}}(), zeros(nunits),
        fill(false, nunits), Bool(screen),
        fill(NaN, nunits), fill(NaN, nunits), fill(false, nunits), NaN, nothing,
        0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, nothing, nothing)
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
    :importance)

"""
    ctsem_laplace_continuation_recentre!(o, values; force=nothing,
        point_units=nothing, point_values=nothing)

Re-place `o`'s nodes at `values`: every unit is scored again, the flagged set
grows by any unit whose gap now matters and by any `force` names -- it never
shrinks, so no later round loses a unit an earlier one was scored on -- and
every flagged unit's rule moves to its mode and curvature there, an importance
rule drawing afresh with the points it knows and any given here. The placement
it replaces is kept for `ctsem_laplace_continuation_revert!`. Returns
`ctsem_laplace_continuation_info`.
"""
function ctsem_laplace_continuation_recentre!(o::CTSEMLaplaceContinuation,
    values::AbstractVector; force=nothing, point_units=nothing, point_values=nothing)
    o.previous = NamedTuple{_CONTINUATION_STATE}(
        Tuple(getfield(o, f) for f in _CONTINUATION_STATE))
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
        draws_seed=o.draws_seed, points=o.points,
        only=o.screen ? nothing : Set(o.flagged))
    placed = _continuation_place(laplace, theta, opts)
    nunits = length(laplace.units.members)
    failed = [placed.rules[U] === nothing && !isempty(laplace.units.blocks[U]) &&
              !placed.wide[U] for U in 1:nunits]
    o.wide = placed.wide
    gaps = _continuation_resolved(placed.quadrature .- placed.laplace, placed.se)
    flagged = sort!(union(o.flagged, _continuation_flags(gaps, o.tolerance)))
    # A unit flagged before whose rule could not be placed here keeps its
    # Laplace term: an honest Laplace term beats a rule placed somewhere else.
    filter!(U -> placed.rules[U] !== nothing, flagged)
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
    o.rule_failed = failed
    prior = _ctsem_log_prior(laplace.objective, theta)
    o.centre_value = sum(_laplace_unit_weight(laplace, U) * placed.laplace[U]
        for U in o.rest_units; init=0.0) +
        sum(_laplace_unit_weight(laplace, U) * placed.quadrature[U]
        for U in flagged; init=0.0) + prior
    o.recentres += 1
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
    _continuation_rule_value(o, U, theta, Ls, dL, positions, rule, u, aws,
        scratch, want_gradient)

One block's log integral under its fixed rule at `theta`, with its gradient:
`log sum_j exp(w_j + h_j(theta)) + c`, whose derivative with the nodes held is
`sum_j pi_j dh_j/dtheta`, `pi` the normalised node weights -- the Fisher
identity. One pass, the running maximum rescaling what has been summed so far.
A node whose likelihood is not finite carries no weight, as it does in
`_quadrature_block`. A leaf's nodes are summed by `_continuation_leaf_sum`.
"""
function _continuation_rule_value(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, want_gradient::Bool)
    if isempty(rule.children)
        a = _continuation_leaf_sum(o, U, theta, Ls, dL, positions, rule, u, aws,
            scratch, want_gradient, 1:size(rule.points, 2))
        isfinite(a.peak) || return (value=NaN, gradient=a.gacc)
        want_gradient && (a.gacc ./= a.acc)
        return (value=a.peak + log(a.acc) + rule.constant, gradient=a.gacc)
    end
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    npar = length(theta)
    n = size(rule.points, 2)
    peak = -Inf
    acc = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    gh = want_gradient ? zeros(npar) : Float64[]
    for j in 1:n
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        want_gradient && fill!(gh, 0.0)
        h = 0.0
        for child in rule.children[j]
            r = _continuation_rule_value(o, U, theta, Ls, dL, positions, child,
                u, aws, scratch, want_gradient)
            # A child that cannot be scored fails the block, rather than
            # silently removing this node's mass.
            isfinite(r.value) || return (value=NaN, gradient=gacc)
            h += r.value
            want_gradient && (gh .+= r.gradient)
        end
        t = rule.logweights[j] + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                acc *= s
                want_gradient && (gacc .*= s)
            end
            peak = t
        end
        w = exp(t - peak)
        acc += w
        want_gradient && (gacc .+= w .* gh)
    end
    isfinite(peak) || return (value=NaN, gradient=gacc)
    want_gradient && (gacc ./= acc)
    return (value=peak + log(acc) + rule.constant, gradient=gacc)
end

"""
    _continuation_leaf_sum(o, U, theta, Ls, dL, positions, rule, u, aws,
        scratch, want_gradient, range)

A leaf rule's nodes `range`, summed: their largest log term `peak`, the sum
`acc` of `exp(term - peak)` over them, and with the gradient the same sum of
each node's score times its weight (`gacc`), unnormalised -- so a rule's nodes
can be summed in parts on different threads and the parts combined
(`_continuation_flagged`). `peak` is `-Inf` where no node in the range scores.
"""
function _continuation_leaf_sum(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, want_gradient::Bool, range::UnitRange{Int})
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    columns = (block.offset + 1):(block.offset + block.size)
    npar = length(theta)
    peak = -Inf
    acc = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    gh = want_gradient ? zeros(npar) : Float64[]
    for j in range
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        want_gradient && fill!(gh, 0.0)
        h = 0.0
        ok = true
        for m in block.members
            ll = _continuation_member!(gh, laplace, U, m, theta, Ls, dL,
                positions, u, aws, scratch, want_gradient)
            if !isfinite(ll)
                ok = false
                break
            end
            h += ll
        end
        ok || continue
        t = rule.logweights[j] + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                acc *= s
                want_gradient && (gacc .*= s)
            end
            peak = t
        end
        w = exp(t - peak)
        acc += w
        want_gradient && (gacc .+= w .* gh)
    end
    return (peak=peak, acc=acc, gacc=gacc)
end

"""
    _continuation_tasks(o, nslot)

The pieces `_continuation_flagged` hands the pool: each flagged unit whole, as
`(f, nothing)`, except that where there are fewer flagged units than slots a
unit whose rule is one leaf of at least 64 nodes -- a wide unit's draws, a
product rule over several effects -- is split into ranges of its nodes,
`(f, range)`, at least 32 a range and up to the slots each unit's share
leaves. Without that a fit whose one flagged unit has 4096 draws evaluated
them on one thread while the rest of the pool waited.
"""
function _continuation_tasks(o::CTSEMLaplaceContinuation, nslot::Integer)
    nf = length(o.rules)
    share = nf < nslot ? cld(Int(nslot), nf) : 1
    tasks = Tuple{Int,Union{Nothing,UnitRange{Int}}}[]
    for (f, rule) in enumerate(o.rules)
        root = rule.roots[1]
        n = size(root.points, 2)
        if share > 1 && length(rule.roots) == 1 && isempty(root.children) && n >= 64
            stop = 0
            for m in _continuation_split(n, min(share, cld(n, 32)))
                push!(tasks, (f, (stop + 1):(stop + m)))
                stop += m
            end
        else
            push!(tasks, (f, nothing))
        end
    end
    return tasks
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
                a = _continuation_leaf_sum(o, U, theta, Ls, dL, positions,
                    rule.roots[1], u, aws, scratch, want_gradient, range)
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
            want_gradient && (g .= grads[at])
        else
            P = maximum(@view peaks[at:last])
            isfinite(P) || return (values=values, gradient=gradient, ok=false)
            A = 0.0
            want_gradient && fill!(g, 0.0)
            for q in at:last
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
# objective, where the optimiser takes the hybrid whole (a fit rescued from a
# jump, `.ctLaplaceFrozenRounds()`): chord, at the cost of the Laplace part's
# differences.
_ctsem_finish_curvature(::CTSEMLaplaceContinuation) = :chord
_ctsem_hessian_cost(::CTSEMLaplaceContinuation, npar::Integer) = 2 * Int(npar)
ctsem_hessian(o::CTSEMLaplaceContinuation, values::AbstractVector;
    chunk::Integer=0, progress=nothing) = ctsem_laplace_continuation_hessian(o, values)
_ctsem_params(o::CTSEMLaplaceContinuation) = _ctsem_params(o.laplace)
_ctsem_saturated_for(o::CTSEMLaplaceContinuation, minimizer) =
    _laplace_saturated_parameters(o.laplace, minimizer)
_ctsem_optimise_result_extra(o::CTSEMLaplaceContinuation, final, log) =
    (continuation_flagged=length(o.flagged),)

function _ctsem_optimise_trial(o::CTSEMLaplaceContinuation, x, want_gradient::Bool,
    gradient_method, limit::Real, log)
    evaluated = try
        ctsem_laplace_continuation_evaluate(o, x; gradient=want_gradient)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    valid = evaluated !== nothing && isfinite(evaluated.value) && evaluated.converged
    if valid && want_gradient
        valid = all(isfinite, evaluated.gradient) &&
            all(abs(entry) < limit for entry in evaluated.gradient)
    end
    return (evaluated=evaluated, valid=valid)
end

function _ctsem_probe_value(o::CTSEMLaplaceContinuation, x)
    evaluated = try
        ctsem_laplace_continuation_evaluate(o, x; gradient=false)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    evaluated === nothing && return -Inf
    (evaluated.converged && isfinite(evaluated.value)) ? evaluated.value : -Inf
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
function ctsem_laplace_continuation_optimize(o::CTSEMLaplaceContinuation,
    values::AbstractVector, basis::AbstractMatrix, radius::Real;
    maxiter::Integer=100, tol::Real=1e-6, stationary_only::Bool=false)
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
        boundary=false, gradient_norm=NaN, refused=0, start_gain=start_gain,
        stationary=stationary)
    (usable && !stationary && !stationary_only) || return none
    newton = sqrt(sum(abs2, gy))
    alpha = radius > 0 ? min(newton, Float64(radius)) : newton
    step = CTSEMContinuationStep(o, x0, B, Float64(radius))
    refused = o.refused
    res = ctsem_optimize(step, zeros(size(B, 2)); maxiter=Int(maxiter), g_tol=1e-8,
        f_tol=0.0, x_tol=0.0, gap_tol=Float64(tol), converge_tol=Float64(tol),
        overshoot_probe=:off, stall_window=0, batch=false, newton=false,
        precondition=nothing, lbfgs_diagonal=false,
        # Its basis is whitened by the Hessian, so its first step is Newton's
        # and its prediction is the gain left.
        lbfgs_scaled=true,
        initial_alpha=max(alpha, 1e-12), progress=false)
    y = collect(Float64, res.minimizer)
    moved = sqrt(sum(abs2, y))
    return (minimizer=_continuation_point(step, y), value=res.maximum_loglik,
        start_value=start.value, iterations=res.iterations, f_calls=res.f_calls,
        g_calls=res.g_calls + 1, converged=res.converged,
        stopped_by_gap=res.stopped_by_gap, moved=moved,
        boundary=radius > 0 && moved >= 0.9 * radius, gradient_norm=res.gradient_norm,
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

`filtered` are the own variances of deviations the filter integrates
(`CTSEMHybrid`'s `psi_index`): a level's parameters, but read by the filter
and no part of `L`, so they are seeded with the active ones, and their
entries' `dL` is zero.
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
        sc, depth)

One block's log integral under its fixed rule at `theta`, with its gradient
left in `sc.mean[depth]` and its Hessian in `sc.hess[depth]`: the Louis
identity over the rule's nodes, the node weights normalised in one pass as
`_continuation_rule_value` does, and the mean and scatter of the node
gradients accumulated by West's weighted update, so that scores which are
large and nearly equal -- a stiff direction's, see `_continuation_stiff_rule`
-- do not cancel in a difference of two sums. Skips the nodes the gradient
path skips, and fails a block whose child cannot be scored, as it does.
Returns `(value, failed)`, `failed` meaning the Hessian cannot be formed.

The running-maximum rescale is written out here as it is in
`_continuation_rule_value` and `_quadrature_block`, which carry one moment
and none: an accumulator shared across the three would need a closure or a
buffer per moment in a node loop that allocates nothing.
"""
function _continuation_rule_hessian!(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, d2L, layout, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, sc, depth::Int)
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
    for j in 1:n
        @inbounds for (t, c) in enumerate(columns); u[c] = rule.points[t, j]; end
        fill!(gnode, 0.0)
        fill!(hnode, 0.0)
        h = 0.0
        ok = true
        if leaf
            for m in block.members
                ll, status = _continuation_member_hessian!(gnode, hnode, laplace, U,
                    m, theta, Ls, dL, d2L, layout, u, aws, sc)
                status == 2 && return (value=NaN, failed=true)
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
                isfinite(r.value) || return (value=NaN, failed=false)
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
    isfinite(peak) || return (value=NaN, failed=false)
    invW = 1 / W
    @inbounds for b in 1:npar, a in 1:npar
        hess[a, b] += scatter[a, b] * invW
    end
    return (value=peak + log(W) + rule.constant, failed=false)
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
        _hybrid_on(hyb) ? hyb.psi_index : Int[])
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
    _laplace_parallel(laplace, 1:nf) do f
        # `local` throughout, as in `_continuation_flagged`.
        local rule, U, aws, sc, u, total, failed, r, slot
        slot = _laplace_slot()
        rule = o.rules[f]
        U = rule.unit
        aws = _laplace_workspace!(laplace, S, npar)
        sc = _continuation_hessian_scratch(laplace, S, npar, layout, depth)
        u = zeros(laplace.units.dims[U])
        fill!(sc.unit_gradient, 0.0)
        fill!(sc.unit_hessian, 0.0)
        total = 0.0
        failed = false
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
    ctsem_laplace_continuation_info(o)

What the continuation's rule looks like at its current centre, for the R side.
Never an empty vector, since one deadlocks the bridge: `flagged` is `[0]` when
nothing is flagged, and `nflagged` says so.
"""
function ctsem_laplace_continuation_info(o::CTSEMLaplaceContinuation)
    nunits = length(o.laplace.units.members)
    gaps = o.centre_quadrature .- o.centre_laplace
    # The screen reads the gaps it can resolve (`_continuation_resolved`); the
    # signed `gap` reported is the raw difference.
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
    return (flagged=isempty(o.flagged) ? [0] : copy(o.flagged),
        nflagged=length(o.flagged), nunits=nunits,
        screen=sum(finite; init=0.0), gap=sum(gaps[isfinite.(gaps)]; init=0.0),
        bar=isempty(flaggedgaps) ? 0.0 : minimum(flaggedgaps),
        residual=sum(finite[U] for U in o.rest_units; init=0.0),
        centre=copy(o.centre), value=o.centre_value,
        quadrature=sum(quadrature_units) + prior, prior=prior,
        quadrature_units=quadrature_units, quadrature_subjects=subjects,
        laplace_units=copy(o.centre_laplace),
        laplace=sum(o.centre_laplace) + prior,
        rule_failures=count(o.rule_failed), soft_blocks=softblocks,
        nwide=count(o.wide), maxdim=o.maxdim,
        importance=[o.importance[U] for U in 1:nunits], nimportance=count(o.importance),
        centre_se=copy(o.centre_se),
        importance_draws=sum(size(r.roots[1].points, 2) for r in o.rules
            if o.importance[r.unit]; init=0),
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
held where they were placed: `sqrt(1/ESS - 1/n)` from the weights there, `0`
for a unit on a product rule. The draws are the proposal's of the point they
were placed at, and as the parameters move away from it their weights spread;
a fit rescued from a jump places them again only once this has drifted past
what they were placed to (`.ctLaplaceFrozenRescue()`). In the order of
`ctsem_laplace_continuation_info`'s `flagged`; `[0.0]` when none is flagged.
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
    out = zeros(nf)
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, 1:nf) do f
        local rule, U, root, block, columns, aws, scratch, u, n, terms, peak, w
        rule = o.rules[f]
        U = rule.unit
        o.importance[U] || return true
        root = rule.roots[1]
        block = laplace.units.blocks[U][root.block]
        columns = (block.offset + 1):(block.offset + block.size)
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        n = size(root.points, 2)
        terms = fill(-Inf, n)
        for j in 1:n
            @inbounds for (t, c) in enumerate(columns); u[c] = root.points[t, j]; end
            h = 0.0
            for m in block.members
                h += _continuation_member!(Float64[], laplace, U, m, theta, Ls,
                    nothing, positions, u, aws, scratch, false)
            end
            isfinite(h) && (terms[j] = root.logweights[j] + h)
        end
        peak = maximum(terms)
        if !isfinite(peak)
            out[f] = Inf
            return true
        end
        w = exp.(terms .- peak)
        out[f] = sqrt(max(sum(abs2, w) / sum(w)^2 - 1 / n, 0.0))
        return true
    end
    return out
end
