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
`d = 6` that is an hour, so the default stops at 5. A unit wider
than `maxdim` is not scored and keeps its Laplace term (`wide` in the info):
the step correction, the alternative, pays the same `nodes^d` per quadrature
value, so falling back to it would save nothing.
"""

using LinearAlgebra

export ctsem_laplace_continuation, ctsem_laplace_continuation_info,
    ctsem_laplace_continuation_evaluate, ctsem_laplace_continuation_recentre!,
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
    # Units wider than `maxdim`: never scored, never flagged.
    wide::Vector{Bool}
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
end

"""Gauss-Hermite grids a rule of these options can ask for, filled before any
threaded region so that `_gh_grid` is a read there (see
`_quadrature_warm_caches`)."""
function _continuation_warm_caches(laplace::CTSEMLaplaceObjective, nodes::Integer,
    soft_maxdirs::Integer)
    _quadrature_warm_caches(laplace, nodes)
    for k in 1:max(1, Int(soft_maxdirs))
        _gh_grid(k, nodes)
    end
    return nothing
end

"""
    _continuation_block_functions(laplace, U, theta, Ls, b, u, aws)

Block `b`'s own log integrand and its precision, as functions of the block's
coordinates with the rest of `u` held where it is now: `value_gradient(z)` is
its members' log likelihood less `z'z/2` and that function's gradient, and
`precision(z)` is `I - d2 ll/dz dz`, dense -- the closures
`_quadrature_leaf_rule!` builds for the same block.
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
        r = _laplace_unit_loglik_gradient(laplace, U, theta, Ls, work, aws, members)
        isfinite(r.value) || return (value=-Inf, gradient=fill(NaN, k))
        return (value=r.value - dot(z, z) / 2,
            gradient=[r.gradient[c] for c in columns] .- z)
    end
    loglik_gradient = function (z)
        S = eltype(z)
        ws = _laplace_workspace!(laplace, S, length(theta))
        work = convert(Vector{S}, base)
        @inbounds for (t, c) in enumerate(columns); work[c] = z[t]; end
        r = _laplace_unit_loglik_gradient(laplace, U, convert(Vector{S}, theta),
            [convert(Matrix{S}, L) for L in Ls], work, ws, members)
        return [r.gradient[c] for c in columns]
    end
    precision = function (z)
        A = ForwardDiff.jacobian(loglik_gradient, collect(Float64, z))
        return Matrix{Float64}(LinearAlgebra.I, k, k) .- _laplace_symmetrise(A)
    end
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
                inner += laplace.objective.subject_objectives[
                    laplace.units.members[U][m]](shifted)
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
    width(b) = blocks[b].size + maximum((width(c) for c in tree.children[b]); init=0)
    return maximum((width(r) for r in tree.roots); init=0)
end

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
    return (rule=CTSEMFixedUnitRule(U, roots, evaluations), value=total)
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
    _continuation_warm_caches(laplace, opts.nodes, opts.soft_maxdirs)
    rules = Vector{Union{Nothing,CTSEMFixedUnitRule}}(nothing, nunits)
    values = fill(NaN, nunits)
    wide = [_continuation_unit_width(laplace, U) > opts.maxdim for U in 1:nunits]
    _laplace_parallel(laplace, 1:nunits) do U
        wide[U] && return true
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
        end
        return true
    end
    return (rules=rules, quadrature=values, laplace=copy(lap.unit_loglik),
        laplace_value=lap.value, converged=lap.converged, wide=wide)
end

"""
    _continuation_flags(gaps, tolerance)

The units the continuation puts nodes on: largest `|gap|` first, the smallest
set whose removal leaves the rest carrying at most `tolerance` in total. That
is the tolerance the screen passes a whole fit on, so what the hybrid leaves to
Laplace is what the screen would have called exact. A unit whose rule could
not be placed (`NaN`) is never flagged.
"""
function _continuation_flags(gaps::Vector{Float64}, tolerance::Real)
    g = [isfinite(x) ? abs(x) : 0.0 for x in gaps]
    order = sortperm(g; rev=true)
    remaining = sum(g; init=0.0)
    flagged = Int[]
    for U in order
        remaining <= tolerance && break
        g[U] > 0 || break
        push!(flagged, U)
        remaining -= g[U]
    end
    return sort!(flagged)
end

"""
    ctsem_laplace_continuation(laplace, values; nodes=5, tolerance=0.01,
        product_maxdim=2, soft_tau=3.5, soft_maxdirs=2, maxdim=5)

The hybrid objective with its nodes placed at `values`. Every unit is scored by
the rule there against its Laplace term, the units that differ are flagged
(`_continuation_flags`), and those carry fixed nodes from then on while the
rest keep their Laplace term. A unit wider than `maxdim` effects along any
path of its block tree is not scored at all and keeps its Laplace term. See the
module docstring for the rule and the constants.
"""
function ctsem_laplace_continuation(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5, tolerance::Real=0.01,
    product_maxdim::Integer=2, soft_tau::Real=3.5, soft_maxdirs::Integer=2,
    maxdim::Integer=5)
    nodes >= 1 || throw(ArgumentError("need at least one quadrature node"))
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    nunits = length(laplace.units.members)
    o = CTSEMLaplaceContinuation{typeof(laplace)}(laplace, laplace, collect(1:nunits),
        Int[], CTSEMFixedUnitRule[], theta, Int(nodes), Float64(tolerance),
        Int(product_maxdim), Float64(soft_tau), Int(soft_maxdirs), Int(maxdim),
        fill(false, nunits),
        fill(NaN, nunits), fill(NaN, nunits), fill(false, nunits), NaN, nothing,
        0, 0, 0, 0, 0, 0)
    _continuation_place!(o, theta)
    return o
end

const _CONTINUATION_STATE = (:rest, :rest_units, :flagged, :rules, :centre,
    :centre_quadrature, :centre_laplace, :rule_failed, :centre_value)

"""
    ctsem_laplace_continuation_recentre!(o, values)

Re-place `o`'s nodes at `values`: every unit is scored again, the flagged set
grows by any unit whose gap now matters -- it never shrinks, so no later round
loses a unit an earlier one was scored on -- and every flagged unit's rule
moves to its mode and curvature there. The placement it replaces is kept for
`ctsem_laplace_continuation_revert!`. Returns `ctsem_laplace_continuation_info`.
"""
function ctsem_laplace_continuation_recentre!(o::CTSEMLaplaceContinuation,
    values::AbstractVector)
    o.previous = NamedTuple{_CONTINUATION_STATE}(
        Tuple(getfield(o, f) for f in _CONTINUATION_STATE))
    _continuation_place!(o, collect(Float64, values))
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
    return ctsem_laplace_continuation_info(o)
end

function _continuation_place!(o::CTSEMLaplaceContinuation, theta::Vector{Float64})
    laplace = o.laplace
    opts = (nodes=o.nodes, product_maxdim=o.product_maxdim, soft_tau=o.soft_tau,
        maxdim=o.maxdim,
        soft_maxdirs=o.soft_maxdirs)
    placed = _continuation_place(laplace, theta, opts)
    nunits = length(laplace.units.members)
    failed = [placed.rules[U] === nothing && !isempty(laplace.units.blocks[U]) &&
              !placed.wide[U] for U in 1:nunits]
    o.wide = placed.wide
    gaps = placed.quadrature .- placed.laplace
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
    o.rule_failed = failed
    prior = _ctsem_log_prior(laplace.objective, theta)
    o.centre_value = sum(placed.laplace[U] for U in o.rest_units; init=0.0) +
        sum(placed.quadrature[U] for U in flagged; init=0.0) + prior
    o.recentres += 1
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
mode, taken here at a node.
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
    want_gradient || return subject(shifted)
    grad = scratch.grad
    loglik = _laplace_subject_value_gradient!(grad, subject, aws, shifted)
    isfinite(loglik) || return loglik
    all(isfinite, grad) || return NaN
    @inbounds for t in eachindex(grad)
        gh[t] += grad[t]
    end
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
            gh[pos[t]] += acc
        end
    end
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
`_quadrature_block`.
"""
function _continuation_rule_value(o::CTSEMLaplaceContinuation, U::Integer,
    theta::Vector{Float64}, Ls, dL, positions, rule::CTSEMFixedRule,
    u::Vector{Float64}, aws, scratch, want_gradient::Bool)
    laplace = o.laplace
    block = laplace.units.blocks[U][rule.block]
    k = block.size
    columns = (block.offset + 1):(block.offset + k)
    npar = length(theta)
    n = size(rule.points, 2)
    leaf = isempty(rule.children)
    peak = -Inf
    acc = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    gh = want_gradient ? zeros(npar) : Float64[]
    for j in 1:n
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
                # A child that cannot be scored fails the block, rather than
                # silently removing this node's mass.
                isfinite(r.value) || return (value=NaN, gradient=gacc)
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

"""The flagged units' fixed-rule terms at `theta`, and their summed gradient,
in parallel over units."""
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
    _continuation_warm_caches(laplace, o.nodes, o.soft_maxdirs)
    nslot = _laplace_pool_width()
    partials = [zeros(npar) for _ in 1:(want_gradient ? nslot : 0)]
    good = fill(true, nslot)
    _laplace_parallel(laplace, 1:nf) do f
        # `local` throughout: a name assigned here that is also a local of the
        # enclosing function would be one binding shared by every task.
        local rule, U, aws, scratch, u, total, grad, r, failed
        rule = o.rules[f]
        U = rule.unit
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        total = 0.0
        grad = want_gradient ? zeros(npar) : Float64[]
        failed = false
        try
            for root in rule.roots
                r = _continuation_rule_value(o, U, theta, Ls, dL, positions, root,
                    u, aws, scratch, want_gradient)
                if !isfinite(r.value)
                    failed = true
                    break
                end
                total += r.value
                want_gradient && (grad .+= r.gradient)
            end
        catch err
            _ctsem_must_propagate(err) && rethrow()
            failed = true
        end
        if failed || (want_gradient && !all(isfinite, grad))
            good[_laplace_slot()] = false
            return false
        end
        values[f] = total
        want_gradient && (partials[_laplace_slot()] .+= grad)
        return true
    end
    ok = all(good)
    gradient = zeros(npar)
    if want_gradient && ok
        for p in partials
            gradient .+= p
        end
    end
    return (values=values, gradient=gradient, ok=ok)
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
    theta = collect(Float64, values)
    rest = ctsem_laplace_evaluate(o.rest, theta; gradient=gradient)
    flagged = _continuation_flagged(o, theta, gradient)
    gradient ? (o.gradient_calls += 1) : (o.value_calls += 1)
    evaluations = sum(r.evaluations for r in o.rules; init=0)
    gradient ? (o.member_sweeps += evaluations) : (o.member_values += evaluations)
    value = flagged.ok ? rest.value + sum(flagged.values; init=0.0) : NaN
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
        tune_chunks=false, precondition=nothing,
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
    ctsem_laplace_continuation_hessian(o, values; step=1e-4)

The Hessian of the hybrid at `values`, its nodes where they are: a central
difference of its exact gradient, as `ctsem_laplace_hessian` takes the Laplace
one, the Laplace part's inner modes warm-started from their modes at `values`
for the reason given there. A column whose two points did not both evaluate is
`NaN`, which the caller sees as an unusable Hessian.
"""
function ctsem_laplace_continuation_hessian(o::CTSEMLaplaceContinuation,
    values::AbstractVector; step::Real=1e-4)
    x = collect(Float64, values)
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
        for j in 1:n
            h = step * max(1.0, abs(x[j]))
            plus = copy(x); plus[j] += h
            minus = copy(x); minus[j] -= h
            restore!()
            gp = ctsem_laplace_continuation_evaluate(o, plus; gradient=true)
            restore!()
            gm = ctsem_laplace_continuation_evaluate(o, minus; gradient=true)
            (gp.converged && gm.converged) || continue
            H[:, j] = (gp.gradient .- gm.gradient) ./ (2h)
        end
    finally
        ctsem_set_warm_start!(previous)
        restore!()
    end
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
    finite = [isfinite(g) ? abs(g) : 0.0 for g in gaps]
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
        evaluations_per_value=sum(r.evaluations for r in o.rules; init=0),
        value_calls=o.value_calls, gradient_calls=o.gradient_calls,
        member_values=o.member_values, member_sweeps=o.member_sweeps,
        recentres=o.recentres, refused=o.refused, nodes=o.nodes,
        product_maxdim=o.product_maxdim, soft_tau=o.soft_tau,
        soft_maxdirs=o.soft_maxdirs, tolerance=o.tolerance)
end
