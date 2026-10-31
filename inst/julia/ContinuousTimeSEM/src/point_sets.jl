"""
The fixed point sets a unit's rule is built from: product Gauss-Hermite and
lattice-t, for the quadrature continuation's screen (`laplace_continuation.jl`)
and the refined fit's rules placed at every evaluation (`laplace_placed.jl`).

# The form

A block's integral over its standardised effects is taken at its mode `c` and
clipped scale `S` (`_quadrature_leaf_rule!`, `_quadrature_clipped_scale`) as

    log sum_j exp(lw_j + f(c + S y_j)) + log|S| - k log(2 pi) / 2

`f` the members' log likelihood less `z'z/2`. A point set is the `y_j`, in the
coordinates `S` whitens, and the `lw_j`, the log weights of a rule for
`integral g(y) dy`. Nothing in the placement or its gradient reads more than
that, so the placed rule's derivative (`_placed_leaf!`) is the same for every
set. A set built for a block of one unit is fixed for that unit: the same
points at every placement, so the placed value is a smooth function of the
parameters wherever the unit's posterior is one basin.

# Product Gauss-Hermite

`y_j = sqrt(2) x_j` on the Gauss-Hermite grid, `lw_j` the grid's log weights
with `x_j'x_j` and `k log(2) / 2` folded in: exact for a Gaussian integrand,
and on a near-Gaussian one far more accurate than its `nodes^k` points suggest.
Its cost grows by `nodes` with every effect.

# Lattice-t

A rank-1 lattice, `frac(i z / n + Delta_r)`, `i = 0, ..., n - 1`, for each of
`_LATTICE_SHIFTS` shifts `Delta_r` seeded by the unit, the block and the shift
(`_point_set_lattice`), never by where the rule sits. Each coordinate goes
through the quantile of a Student t on `_LATTICE_DF` degrees of freedom, and
`lw_j` is `-log(n R)` less the log of that product of t densities. The
generating vector `z` is built component by component (`_lattice_cbc`). Its
cost is `n R` whatever the width.

Why this mapping. A multivariate t, one more lattice coordinate for its
chi-square, puts every coordinate's tail at once on the normal quantile, whose
derivative grows as `exp(y^2/2)` towards the cube's faces, faster than the
integrand falls: the integrand in the cube then has unbounded derivatives and
the lattice converged as `n^-0.6`, little better than independent draws. A t
quantile's derivative grows as a power of `y`, so the integrand vanishes
smoothly at every face and the rule converges as about `1/n`: on w4's units
(four effects) 6.9e-5 nats at 16,384 points against 6.3e-4
(`scripts/anywidth-2026-10/proto.jl`). With the integrand vanishing at the
faces it is periodic already, so no tent transform: one halved the
resolution for nothing. Gaussian marginals in place of t are lighter-tailed
than a unit's posterior: the rms of their true error over their stated one
was 1.5-2.5 on w4, against 0.7-1.5 on t. At ten effects (k5's units, against
the same rule at 16,384 points) a lattice of 512 points was within 0.015
nats, where the 3-node product rule needs 59,049 for the same.

The error is the spread over the shifts: `R` independent estimates of the same
integral, so `sd(I_r) / (sqrt(R) mean(I_r))` on the log scale.

# A defensive prior component

A t at the mode and the clipped curvature misses whatever mass sits where it
puts little: a shoulder beside a stiff mode, and every shift misses it alike,
so their spread does not show it. So a lattice rule also takes
`_LATTICE_PRIOR_SHARE` of its points from the prior, `N(0, I)` in the unit's
standardised effects, held fixed there (`prior`), and every point is weighted
by the balance heuristic, `1 / (N (b1 rho_placed(z) + b0 phi(z)))`
(`_point_set_place`). That weight moves with the placement, so the placed
gradient reads it (`_placed_leaf!`).
"""

using LinearAlgebra
using Random

"""
    CTSEMPointSet

A fixed rule over a block's `k` standardised coordinates: `points[:, j] = y_j`,
placed as `z = c + S y_j`, with `integral g(y) dy ~ sum_j exp(logweights[j])
g(y_j)`. `family` is `:product` or `:lattice`, `count` the nodes per
coordinate or the lattice points per shift; `replicate[j]` is the shift point
`j` belongs to (empty for a product rule) and `nrep` their number. A lattice
rule's defensive component (see the module docstring): `prior`, its points in
the standardised effects themselves, held there, `prior_replicate` the shift
each belongs to, and `logrho` each placed point's log t density; empty for a
product rule. `_point_set_place` places either kind.
"""
struct CTSEMPointSet
    family::Symbol
    count::Int
    points::Matrix{Float64}
    logweights::Vector{Float64}
    replicate::Vector{Int}
    nrep::Int
    prior::Matrix{Float64}
    prior_replicate::Vector{Int}
    logrho::Vector{Float64}
end

"""How many shifts a lattice rule takes. Its error is their spread, which over
eight has about a quarter's relative error itself: on w4's twelve units the
rms of the true error over the stated one was 0.7-1.5 at every size."""
const _LATTICE_SHIFTS = 8

"""Degrees of freedom of the lattice rule's t marginals: tails heavier than a
unit's posterior's. As the importance draws' `draws_df`."""
const _LATTICE_DF = 4.0

"""The share of a lattice rule's points from its defensive prior component.
On anomS2's six shoulder units at its spurious start, the lattice-t alone at
1,024 points read them 0.09 low in sum, its worst unit 0.055 low at a stated
0.031; with a quarter of the points at the prior, 0.017 high and 0.013 at
worst, but the summed absolute error over all 800 units rose from 0.48 to
0.95, the quarter spent on easy units
(`review/MIXTURE-multibasin-measure-2026-10-30.md`). A tenth, as phase B
suggests for wide units."""
const _LATTICE_PRIOR_SHARE = 0.1

"""Lattice points a shift a rule starts at, 128 points in all: on w4 that met
the 0.05-nat target on every unit, at a fifth of the 5-node product rule's
625; at ten effects (k5) the target took 256 to 512."""
const _LATTICE_START = 16

"""The most lattice points a shift a rule is refined to, 16,384 points in all,
the importance draws' `draws_max`."""
const _LATTICE_MAX = 2048

"""Points in the set, its prior's included."""
_point_set_size(ps::CTSEMPointSet) = size(ps.points, 2) + size(ps.prior, 2)

"""
    _point_set_product(k, nodes)

The product Gauss-Hermite set over `k` coordinates at `nodes` a coordinate
(see the module docstring), from `_gauss_hermite`, built each call: the
callers keep the set for the units that use it.
"""
function _point_set_product(k::Integer, nodes::Integer)
    big(nodes)^k > _CTSEM_GH_GRID_MAXPOINTS && throw(ArgumentError(string(
        "a product quadrature rule over ", k, " random effects at ", nodes,
        " nodes has ", big(nodes)^k, " points, too many to evaluate")))
    x, w = _gauss_hermite(nodes)
    m = Int(nodes)
    n = m^Int(k)
    Y = Matrix{Float64}(undef, Int(k), n)
    lw = Vector{Float64}(undef, n)
    half = Int(k) * log(2.0) / 2
    for (col, index) in enumerate(Iterators.product(ntuple(_ -> 1:m, Int(k))...))
        s = half
        @inbounds for (a, i) in enumerate(index)
            Y[a, col] = sqrt(2.0) * x[i]
            s += log(w[i]) + x[i]^2
        end
        lw[col] = s
    end
    return CTSEMPointSet(:product, m, Y, lw, Int[], 1, zeros(Int(k), 0), Int[], Float64[])
end

"""
    _lattice_cbc(d, n)

The generating vector of an `n`-point rank-1 lattice in `d` dimensions, `n` a
power of two, built component by component (Sloan, Kuo and Joe 2002): each
component the odd number that minimises the lattice's worst-case error in the
Korobov space of smoothness two, with product weights of one half, given the
components before it. `O(d n^2)` arithmetic, once per `(d, n)` a session
(`_lattice_vector`).
"""
function _lattice_cbc(d::Integer, n::Integer)
    (n >= 2 && ispow2(n)) || throw(ArgumentError("a lattice size must be a power of two"))
    gamma = 0.5
    omega = Vector{Float64}(undef, n)
    @inbounds for i in 0:(n - 1)
        x = i / n
        omega[i + 1] = 2 * pi^2 * (x * x - x + 1 / 6)
    end
    p = ones(n)
    z = Vector{Int}(undef, Int(d))
    for s in 1:Int(d)
        # Every odd number scores the same as the first component, a
        # permutation of the points; one, rather than a tie broken by rounding.
        if s == 1
            z[1] = 1
            @inbounds for i in 0:(n - 1)
                p[i + 1] *= 1 + gamma * omega[i + 1]
            end
            continue
        end
        best = Inf
        bz = 1
        for c in 1:2:(n - 1)
            e = 0.0
            @inbounds for i in 0:(n - 1)
                e += p[i + 1] * (1 + gamma * omega[mod(i * c, n) + 1])
            end
            if e < best
                best = e
                bz = c
            end
        end
        z[s] = bz
        @inbounds for i in 0:(n - 1)
            p[i + 1] *= 1 + gamma * omega[mod(i * bz, n) + 1]
        end
    end
    return z
end

const _LATTICE_VECTORS = Dict{Tuple{Int,Int},Vector{Int}}()
const _LATTICE_VECTORS_LOCK = ReentrantLock()

"""`_lattice_cbc(d, n)`, kept for the session. Locked: a placement reaches it
from the unit loop's threads."""
function _lattice_vector(d::Integer, n::Integer)
    key = (Int(d), Int(n))
    lock(_LATTICE_VECTORS_LOCK) do
        cached = get(_LATTICE_VECTORS, key, nothing)
        cached === nothing || return cached
        z = _lattice_cbc(d, n)
        _LATTICE_VECTORS[key] = z
        return z
    end
end

"""
    _t4_quantile(p)

The quantile of a Student t on four degrees of freedom, in closed form
(`alpha = 4 p (1 - p)`, `q = cos(acos(sqrt(alpha)) / 3) / sqrt(alpha)`, `t =
sign(p - 1/2) 2 sqrt(q - 1)`).
"""
function _t4_quantile(p::Float64)
    alpha = 4 * p * (1 - p)
    q = cos(acos(sqrt(alpha)) / 3) / sqrt(alpha)
    t = 2 * sqrt(max(q - 1, 0.0))
    return p < 0.5 ? -t : t
end

"""The log density of a Student t on four degrees of freedom: `(3/8) (1 +
y^2/4)^(-5/2)`."""
_t4_logdensity(y::Float64) = log(3 / 8) - 2.5 * log1p(y * y / 4)

"""
    _point_set_lattice(k, n, key)

The lattice-t set over `k` coordinates, `n` points a shift and
`_LATTICE_SHIFTS` shifts (see the module docstring); `key` seeds the shifts,
each from `hash((key, r))`. `_LATTICE_DF` is four, the closed form's. Each
shift also takes `ceil(n s / (1 - s))` points of the prior, `s` the
`_LATTICE_PRIOR_SHARE` (`prior`), from a shifted Kronecker lattice through
the normal quantile (`_continuation_lattice_normals`) seeded by `key` and the
shift; `logweights` are then the t's alone, as `share = 0` would have them,
and `_point_set_place` weighs every point by the balance heuristic.
"""
function _point_set_lattice(k::Integer, n::Integer, key; share::Real=_LATTICE_PRIOR_SHARE)
    z = _lattice_vector(k, n)
    R = _LATTICE_SHIFTS
    N = Int(n) * R
    n0 = share > 0 ? ceil(Int, n * share / (1 - share)) : 0
    Y = Matrix{Float64}(undef, Int(k), N)
    lw = Vector{Float64}(undef, N)
    lrho = Vector{Float64}(undef, N)
    rep = Vector{Int}(undef, N)
    base = -log(Float64(N))
    shift = Vector{Float64}(undef, Int(k))
    P0 = Matrix{Float64}(undef, Int(k), n0 * R)
    rep0 = Vector{Int}(undef, n0 * R)
    col = 0
    for r in 1:R
        rand!(Random.Xoshiro(hash((key, r))), shift)
        for i in 0:(Int(n) - 1)
            col += 1
            s = 0.0
            @inbounds for a in 1:Int(k)
                x = clamp(mod(i * z[a] / n + shift[a], 1.0), 1e-15, 1 - 1e-15)
                y = _t4_quantile(x)
                Y[a, col] = y
                s += _t4_logdensity(y)
            end
            lrho[col] = s
            lw[col] = base - s
            rep[col] = r
        end
        if n0 > 0
            cols = ((r - 1) * n0 + 1):(r * n0)
            P0[:, cols] = _continuation_lattice_normals(Int(k), n0,
                Random.Xoshiro(hash((key, :prior, r))))
            rep0[cols] .= r
        end
    end
    return CTSEMPointSet(:lattice, Int(n), Y, lw, rep, R, P0, rep0, lrho)
end

"""
    _point_set_place(ps, c, S, U, logdetS)

The set `ps` placed at centre `c` and scale `S = U^-1` (`log|S| = logdetS`):
each point's position `Z[:, j]` in the standardised effects, its coordinates
`Y[:, j] = U (z - c)` in the scale's, its log weight `logw[j]` as a rule for
`integral g(z) dz`, its responsibility `r[j]` under the placed component, and
whether it is one of the prior's (`prior[j]`) and which shift it belongs to
(`rep`). A product rule's, or a lattice's with no prior points, weighs a point
`lw_j + log|S|`, `r = 1`; with prior points every point is weighed
`-log N - log(b1 rho(y) / |S| + b0 phi(z))` (the balance heuristic, `b1`,
`b0` the shares of points), and `r` is the placed component's share of that
mixture density there.
"""
function _point_set_place(ps::CTSEMPointSet, c::AbstractVector{Float64},
    S::AbstractMatrix{Float64}, U::AbstractMatrix{Float64}, logdetS::Float64)
    k = length(c)
    n1 = size(ps.points, 2)
    n0 = size(ps.prior, 2)
    N = n1 + n0
    Z = Matrix{Float64}(undef, k, N)
    Y = Matrix{Float64}(undef, k, N)
    logw = Vector{Float64}(undef, N)
    r = ones(N)
    prior = fill(false, N)
    rep = ones(Int, N)
    y = zeros(k)
    sy = zeros(k)
    zc = zeros(k)
    mix = n0 > 0
    lb1 = mix ? log(n1 / N) : 0.0
    lb0 = mix ? log(n0 / N) : -Inf
    lN = log(Float64(N))
    lphi0 = -k * log(2 * pi) / 2
    for j in 1:n1
        @inbounds for t in 1:k
            y[t] = ps.points[t, j]
        end
        _ctsem_mulvec!(sy, S, y)
        q = 0.0
        @inbounds for t in 1:k
            Z[t, j] = c[t] + sy[t]
            Y[t, j] = y[t]
            q += Z[t, j]^2
        end
        if mix
            lt = lb1 + ps.logrho[j] - logdetS
            lp = lb0 + lphi0 - q / 2
            lq = max(lt, lp) + log1p(exp(-abs(lt - lp)))
            logw[j] = -lN - lq
            r[j] = exp(lt - lq)
        else
            logw[j] = ps.logweights[j] + logdetS
        end
        isempty(ps.replicate) || (rep[j] = ps.replicate[j])
    end
    for i in 1:n0
        col = n1 + i
        q = 0.0
        @inbounds for t in 1:k
            Z[t, col] = ps.prior[t, i]
            zc[t] = ps.prior[t, i] - c[t]
            q += ps.prior[t, i]^2
        end
        _ctsem_mulvec!(y, U, zc)
        lrho = 0.0
        @inbounds for t in 1:k
            Y[t, col] = y[t]
            lrho += _t4_logdensity(y[t])
        end
        lt = lb1 + lrho - logdetS
        lp = lb0 + lphi0 - q / 2
        lq = max(lt, lp) + log1p(exp(-abs(lt - lp)))
        logw[col] = -lN - lq
        r[col] = exp(lt - lq)
        prior[col] = true
        rep[col] = ps.prior_replicate[i]
    end
    return (Z=Z, Y=Y, logw=logw, r=r, prior=prior, rep=rep)
end

"""The seed of a lattice rule's shifts: fixed, so that every placement of every
unit's rule, at the screen and in the climb, reads the same points."""
const _LATTICE_SEED = 20261030

"""The key a lattice rule for block `b` of unit `U` seeds its shifts by;
`alt` names an independent set of shifts for the same rule, which the
refined fit's slope check reads its gradient's error from
(`_placed_slope_check`)."""
_lattice_key(U::Integer, b::Integer, alt::Integer=0) =
    (_LATTICE_SEED + Int(alt), Int(U), Int(b))

"""
    _point_set_error(repsums, nrep)

A lattice rule's standard error on the log scale from its shifts' sums
`repsums` -- each `sum_{j in r} exp(t_j - peak)`, on one `peak` -- as
`sd(I_r) / (sqrt(R) mean(I_r))`; zero for a product rule (`nrep` 1).
"""
function _point_set_error(repsums::AbstractVector{Float64}, nrep::Integer)
    nrep < 2 && return 0.0
    m = sum(repsums) / nrep
    m > 0 || return Inf
    v = 0.0
    @inbounds for r in 1:nrep
        v += (repsums[r] / m - 1)^2
    end
    return sqrt(v / (nrep - 1) / nrep)
end

"""
    CTSEMRuleChoice

The rule a unit's blocks take: `family` `:product` (every block on product
Gauss-Hermite at `nodes`), `:lattice` (every block on lattice-t at `points` a
shift), or `:mixed`, each block the cheaper of the two at those counts
(`_rule_block_family`). See `_continuation_refine!` for how one is chosen.
"""
struct CTSEMRuleChoice
    family::Symbol
    nodes::Int
    points::Int
end

"""Whether a block of `k` effects takes the product rule under `choice`: under
`:mixed`, where its `nodes^k` points are no more than the lattice's."""
function _rule_block_product(choice::CTSEMRuleChoice, k::Integer)
    choice.family === :product && return true
    choice.family === :lattice && return false
    return big(choice.nodes)^k <= _LATTICE_SHIFTS * choice.points
end

"""
    _rule_default(laplace, U, nodes)

The rule a unit is screened on: for one block, the product rule at `nodes`
where its `nodes^k` points are no more than the lattice's first `128`, the
lattice there otherwise; for a tree of blocks, each block so (`:mixed`).
"""
function _rule_default(laplace::CTSEMLaplaceObjective, U::Integer, nodes::Integer)
    blocks = laplace.units.blocks[U]
    length(blocks) == 1 || return CTSEMRuleChoice(:mixed, Int(nodes), _LATTICE_START)
    k = blocks[1].size
    return big(nodes)^k <= _LATTICE_SHIFTS * _LATTICE_START ?
        CTSEMRuleChoice(:product, Int(nodes), 0) :
        CTSEMRuleChoice(:lattice, Int(nodes), _LATTICE_START)
end

"""
    _rule_sets(laplace, U, choice)

Unit `U`'s point set per block under `choice`, a lattice block's seeded by
`_lattice_key(U, b, alt)`.
"""
function _rule_sets(laplace::CTSEMLaplaceObjective, U::Integer, choice::CTSEMRuleChoice;
    alt::Integer=0)
    blocks = laplace.units.blocks[U]
    sets = Vector{CTSEMPointSet}(undef, length(blocks))
    for (b, block) in enumerate(blocks)
        k = block.size
        sets[b] = _rule_block_product(choice, k) ?
            _point_set_product(k, choice.nodes) :
            _point_set_lattice(k, choice.points, _lattice_key(U, b, alt))
    end
    return sets
end

"""
    _rule_cost(laplace, U, choice)

Member evaluations one value of unit `U` costs under `choice`: each leaf's
members once at every combination of its path's points. What the refine compares rules
by (`_continuation_refine!`), and the count a gradient's reverse sweeps scale
with.
"""
function _rule_cost(laplace::CTSEMLaplaceObjective, U::Integer, choice::CTSEMRuleChoice)
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return 0.0
    npts = Float64[_rule_block_product(choice, b.size) ?
        Float64(choice.nodes)^b.size : Float64(_LATTICE_SHIFTS * (choice.points +
            ceil(Int, choice.points * _LATTICE_PRIOR_SHARE / (1 - _LATTICE_PRIOR_SHARE))))
        for b in blocks]
    tree = _quadrature_children(blocks)
    total = 0.0
    for (b, block) in enumerate(blocks)
        isempty(tree.children[b]) || continue
        path = npts[b]
        for a in block.ancestors
            path *= npts[a]
        end
        total += path * length(block.members)
    end
    return total
end
