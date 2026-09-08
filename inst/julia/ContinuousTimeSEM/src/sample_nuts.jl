"""
No-U-Turn sampling, with a block-structured metric.

# Why this is written here rather than taken from a package

The sampler itself is standard and could have come from AdvancedHMC. The metric
could not. This target has dimension `npar + sum_U dim(u_U)` -- a hundred
subjects with two effects each is a little over two hundred -- and its geometry
is neither diagonal nor dense:

  * a *diagonal* metric throws away the correlation between a subject's own
    effects and between the population parameters, which is where the
    conditioning actually is;
  * a *dense* metric is a 200x200 covariance to estimate from warmup draws and
    factorize every iteration, almost all of it structural zeros -- two
    subjects' effects are conditionally independent given the population
    parameters, and no amount of warmup will discover that as precisely as the
    model already states it.

The right structure is block diagonal, one block per population vector and one
per random-effect block, which is exactly the sparsity `laplace.jl` already
factorizes. Both packages express diagonal and dense; neither expresses this.

# The metric comes from the Laplace fit

`ctsem_sample_metric` reads the initial metric straight off a Laplace fit: the
outer Hessian gives the population block, and each unit's curvature at its mode
-- through the selected inverse the Laplace path already computes -- gives that
unit's blocks. So the chain starts with a metric that is correct to the accuracy
of the Gaussian approximation, and warmup refines it rather than discovering it.
That is the whole performance argument for sampling *after* optimising rather
than instead of it, and it is why `ctsem_sample` takes a fitted object.

# What is implemented

Stan's variant: multinomial sampling from the trajectory, the generalised
U-turn criterion with the extra subtree checks, dual-averaging step size with
Nesterov's smoothing, and windowed metric adaptation. Divergences, tree depth
saturation and E-BFMI are reported rather than silently absorbed -- on this model
they mean something specific (a population scale near zero pinching the
non-centred funnel) and a user needs to see them.
"""

using LinearAlgebra
using Random

export CTSEMMetric, ctsem_sample_metric

################################################################################
# Metric
################################################################################

"""
    CTSEMMetric(ranges, factors)

A block-diagonal inverse mass matrix -- an estimate of the posterior covariance,
which is what a metric should be.

`ranges` partitions `1:ndim`; `factors[b]` is the lower Cholesky factor `C_b` of
that block's covariance `Sigma_b = C_b C_b'`. A diagonal metric is the special
case of 1x1 blocks and needs no separate code path.

Three operations, all block-local:

  * velocity `Sigma p`, as `C (C' p)`;
  * kinetic energy `p'Sigma p / 2`, as `||C' p||^2 / 2`, sharing the same
    intermediate;
  * a momentum draw with covariance `Sigma^-1`, as the solution of `C' p = xi`
    for standard normal `xi`.
"""
struct CTSEMMetric
    ranges::Vector{UnitRange{Int}}
    factors::Vector{Matrix{Float64}}
    ndim::Int
end

function CTSEMMetric(ranges::Vector{UnitRange{Int}}, factors::Vector{Matrix{Float64}})
    length(ranges) == length(factors) ||
        throw(DimensionMismatch("one factor per block is required"))
    ndim = isempty(ranges) ? 0 : maximum(last, ranges)
    for (r, C) in zip(ranges, factors)
        size(C, 1) == size(C, 2) == length(r) ||
            throw(DimensionMismatch("block factor does not match its range"))
    end
    return CTSEMMetric(ranges, factors, ndim)
end

"""An identity metric over `ndim` coordinates, as one 1x1 block each."""
ctsem_identity_metric(ndim::Integer) = CTSEMMetric(
    [i:i for i in 1:Int(ndim)], [ones(1, 1) for _ in 1:Int(ndim)])

"""`v .= Sigma * p`, and return the kinetic energy `p'Sigma p / 2`."""
function _metric_velocity!(v::Vector{Float64}, metric::CTSEMMetric,
    p::Vector{Float64}, scratch::Vector{Float64})
    kinetic = 0.0
    @inbounds for b in eachindex(metric.ranges)
        r = metric.ranges[b]
        C = metric.factors[b]
        k = length(r)
        # w = C' p
        for j in 1:k
            acc = 0.0
            for i in j:k
                acc += C[i, j] * p[r[i]]
            end
            scratch[j] = acc
            kinetic += acc * acc
        end
        # v = C w
        for i in 1:k
            acc = 0.0
            for j in 1:i
                acc += C[i, j] * scratch[j]
            end
            v[r[i]] = acc
        end
    end
    return kinetic / 2
end

"""Draw `p` with covariance `Sigma^-1`, by solving `C' p = xi`."""
function _metric_momentum!(p::Vector{Float64}, metric::CTSEMMetric, rng::AbstractRNG)
    @inbounds for b in eachindex(metric.ranges)
        r = metric.ranges[b]
        C = metric.factors[b]
        k = length(r)
        # Back substitution on the upper triangular C', innermost index last.
        for i in k:-1:1
            acc = randn(rng)
            for j in (i + 1):k
                acc -= C[j, i] * p[r[j]]
            end
            p[r[i]] = acc / C[i, i]
        end
    end
    return p
end

"""
    _metric_from_covariances(ranges, covariances; jitter)

Factorize each block, falling back to its diagonal when it is not positive
definite.

A block estimated from too few draws, or a curvature that needed repair, can
fail to factorize. That is not a reason to stop: a diagonal metric for that
block is worse than the intended one and far better than none, and the sampler
stays correct either way -- the metric affects efficiency, never the stationary
distribution.
"""
function _metric_from_covariances(ranges::Vector{UnitRange{Int}},
    covariances::Vector{Matrix{Float64}}; jitter::Real=1e-10)
    factors = Vector{Matrix{Float64}}(undef, length(ranges))
    for b in eachindex(ranges)
        S = covariances[b]
        k = size(S, 1)
        S = (S .+ transpose(S)) ./ 2
        for i in 1:k
            S[i, i] += jitter
        end
        F = cholesky(Symmetric(S); check=false)
        factors[b] = if issuccess(F)
            Matrix(F.L)
        else
            d = [sqrt(max(S[i, i], jitter)) for i in 1:k]
            Matrix(Diagonal(d))
        end
    end
    return CTSEMMetric(ranges, factors)
end

"""
    _bounded_inverse(information; rtol)

Invert a symmetric information matrix along the directions it actually
identifies, flooring the rest.

A plain `inv` here is not merely inaccurate, it is unusable. On a weakly
identified fit the outer Hessian has eigenvalues at 1e-83 and, from numerical
error, sometimes small *positive* ones where it should be negative definite --
measured on a six-subject fixture whose population scales had collapsed. The
inverse then has variances of 1e83, the initial draw lands there, and the chain
cannot start.

Flooring the eigenvalues at `rtol` times the largest gives a metric that is
merely uninformative in those directions rather than catastrophic. That is the
right failure mode: the metric only ever affects efficiency, so a bad one costs
time, while an infinite one costs the run.
"""
function _bounded_inverse(information::Symmetric{Float64}; rtol::Real=1e-8)
    n = size(information, 1)
    n == 0 && return zeros(0, 0)
    decomposition = try
        eigen(information)
    catch err
        err isa InterruptException && rethrow()
        return Matrix(1.0I, n, n)
    end
    lambda = decomposition.values
    scale = maximum(abs, lambda)
    (!isfinite(scale) || scale <= 0) && return Matrix(1.0I, n, n)
    floor = rtol * scale
    inverted = [1.0 / max(l, floor) for l in lambda]
    out = decomposition.vectors * Diagonal(inverted) * transpose(decomposition.vectors)
    all(isfinite, out) || return Matrix(1.0I, n, n)
    return (out .+ transpose(out)) ./ 2
end

"""
    _conditional_population_covariance(sampler, theta, Ls, fallback)

The population block of the metric: `(-H_theta,theta)^-1` for the *joint* log
density, at the Laplace mode.

Central differences of the joint gradient's population part, which costs
`2 * npar` joint gradient evaluations -- a few dozen milliseconds, once, against
a sampler that will spend minutes. Differencing the gradient rather than the
value keeps the accuracy the adjoint already has.

`fallback` is the marginal Hessian, used only if the difference fails. It is the
wrong matrix for this purpose, as the caller explains, but a wrong metric costs
efficiency where no metric costs the run.
"""
function _conditional_population_covariance(sampler::CTSEMSampler,
    theta::Vector{Float64}, Ls::Vector{<:AbstractMatrix},
    fallback::Union{Nothing,AbstractMatrix}; step::Real=1e-4)
    n = sampler.npar
    laplace = sampler.laplace
    # Evaluated at the joint mode, not at zero effects: the conditional
    # curvature of `theta` is a local quantity and the modes are where the
    # posterior actually is.
    x = ctsem_sample_start(sampler, theta)
    for U in 1:sampler.nunits
        mode = laplace.modes[U]
        length(mode) == sampler.udims[U] || continue
        @inbounds for q in eachindex(mode)
            x[sampler.uoffsets[U] + q] = mode[q]
        end
    end

    H = zeros(n, n)
    gplus = zeros(sampler.ndim)
    gminus = zeros(sampler.ndim)
    ok = true
    for j in 1:n
        original = x[j]
        x[j] = original + step
        vplus = ctsem_sample_density!(gplus, sampler, x)
        x[j] = original - step
        vminus = ctsem_sample_density!(gminus, sampler, x)
        x[j] = original
        if !isfinite(vplus) || !isfinite(vminus)
            ok = false
            break
        end
        @inbounds for i in 1:n
            H[i, j] = (gplus[i] - gminus[i]) / (2 * step)
        end
    end
    if !ok
        fallback === nothing && return Matrix(1.0I, n, n)
        F = Matrix(fallback)
        return _bounded_inverse(Symmetric((-(F .+ transpose(F))) ./ 2))
    end
    return _bounded_inverse(Symmetric((-(H .+ transpose(H))) ./ 2))
end

"""
    ctsem_sample_metric(sampler, values; regularize)

The initial metric, read off the Laplace approximation at `values`.

The population block is `inv(-H)` for the outer Hessian `H`; each random-effect
block is the corresponding diagonal block of the inverse unit curvature, which
is that block's conditional covariance under the Gaussian approximation and is
already computed by the Laplace factorization.

This is the difference between a chain that starts well conditioned and one that
spends its warmup finding out what the model could have told it.

The per-unit modes are re-solved at `values` first, so that what comes back
describes `values` and not wherever the objective was last evaluated. The body
records what reading them instead was measured to cost.
"""
function ctsem_sample_metric(sampler::CTSEMSampler, values::AbstractVector;
    hessian::Union{Nothing,AbstractMatrix}=nothing, regularize::Real=1e-8)
    laplace = sampler.laplace
    theta = collect(Float64, values)[1:sampler.npar]
    ranges = UnitRange{Int}[]
    covariances = Matrix{Float64}[]

    # Population block, from the *conditional* curvature and not the marginal.
    #
    # This distinction is the whole difference between a metric that works here
    # and one that does not. `ctsem_laplace_hessian` gives the curvature of the
    # marginal posterior, with the effects integrated out; a block-diagonal
    # metric has no term coupling `theta` to `u`, so what each block needs is
    # the curvature *conditional* on the others. The two differ by exactly the
    # coupling the block form drops -- the marginal is the conditional minus
    # `H_tu H_uu^-1 H_ut`, so it is always the flatter of the two and using it
    # over-scales every population direction. Measured on a 30-subject model,
    # the marginal metric gave R-hat 1.11 and 21 effective draws from 1200.
    #
    # It is the same mistake, in the same place, that the nested quadrature made
    # by scaling its outer blocks with the marginal covariance rather than the
    # eliminated diagonal.
    Ls = _laplace_popchols(theta, laplace.spec)

    # The per-unit modes are solved here, at `theta`, rather than read as they
    # stand -- and this is the whole difference between a chain that starts at
    # the fit and one that does not come back.
    #
    # `laplace.modes` is retained *state*: a warm start for the next inner
    # Newton solve, holding whatever the objective's last call left there (see
    # `CTSEMLaplaceObjective`). Three functions downstream of here read it as
    # though it described `theta` -- this metric's effect blocks, the
    # conditional population block below, which evaluates the joint curvature
    # at those effects, and `ctsem_sample_start`, which places the chain's
    # starting effects at them. A freshly constructed objective has them all at
    # zero, and one last evaluated elsewhere has that elsewhere's modes.
    #
    # Both states are reached by ordinary use, because the objective a chain
    # samples need not be the one that was optimised. `ctSample()` on a fit
    # reloaded into a new session builds a fresh one. Worse, every chain of
    # `ctFit(optimize = FALSE)` run as a process does: the worker pool is warmed
    # with one gradient at the *pre-optimisation start values*, deliberately, so
    # that the engine compiles while the optimisation it overlaps is still
    # running -- and the chain then samples with `theta` at the estimate and
    # every subject's effects at whatever was modal for the start.
    #
    # Measured on a 100-subject, 50-occasion, 311-dimension fit whose joint
    # density at the mode is -5297, reading the modes rather than solving them:
    #
    #     objective state                      centre     chain 1 start
    #     last evaluated at the estimate      -5296.99         -5471.27
    #     freshly built (all modes zero)    -411899.24        -4.70e+09
    #     last evaluated at the start values -114649.76        -4.00e+14
    #
    # Nothing announced any of it. The metric was plausible, the density finite,
    # the gradient finite, and the chain went nowhere: reported from a real run
    # as two chains that between them had made twelve transitions in forty
    # minutes, one sitting at logp -2.3e6 and the other at -3e4, while the
    # optimisation they started from had reported -5834.
    #
    # The solve costs one inner Newton pass per unit, once per sample run,
    # against the thousands of joint gradients that follow. It warm-starts from
    # whatever is there and retries from the origin when that fails, so a bad
    # warm start costs a few iterations rather than the answer.
    for U in 1:sampler.nunits
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    end

    popcov = _conditional_population_covariance(sampler, theta, Ls, hessian)
    push!(ranges, 1:sampler.npar)
    push!(covariances, popcov)

    # One block per random-effect block, from the unit curvature at the mode --
    # already conditional, being the inverse of that block's own curvature with
    # everything else held.
    for U in 1:sampler.nunits
        blocks = laplace.units.blocks[U]
        isempty(blocks) && continue
        u = laplace.modes[U]
        Cdiag = try
            M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
                _laplace_unit_curvature(laplace, U, theta, Ls, u, 1)
            fac = _laplace_factor_repaired!(M, blocks)
            fac.ok ? first(_laplace_selected_inverse(fac.factors, fac.coupling, blocks)) :
                nothing
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        for (b, block) in enumerate(blocks)
            block.size == 0 && continue
            start = sampler.uoffsets[U] + block.offset + 1
            push!(ranges, start:(start + block.size - 1))
            candidate = Cdiag === nothing ? Matrix(1.0I, block.size, block.size) :
                Matrix(Cdiag[b])
            # Same guard as the population block: a repaired curvature can give
            # a conditional covariance that is not usable as a metric.
            if !all(isfinite, candidate) || any(candidate[i, i] <= 0 for i in 1:block.size)
                candidate = Matrix(1.0I, block.size, block.size)
            end
            push!(covariances, candidate)
        end
    end

    # Any coordinate the blocks above did not cover -- there should be none, but
    # a metric that silently omits a dimension would give it zero momentum and
    # freeze it, which is a far worse failure than an identity block.
    covered = falses(sampler.ndim)
    for r in ranges, i in r
        covered[i] = true
    end
    for i in 1:sampler.ndim
        covered[i] && continue
        push!(ranges, i:i)
        push!(covariances, ones(1, 1))
    end

    order = sortperm(first.(ranges))
    return _metric_from_covariances(ranges[order], covariances[order];
        jitter=regularize)
end

################################################################################
# Leapfrog and trees
################################################################################

"""One point of a trajectory: position, momentum, gradient, log density."""
mutable struct _NUTSPoint
    x::Vector{Float64}
    p::Vector{Float64}
    g::Vector{Float64}
    logp::Float64
end

_nuts_point(ndim::Int) = _NUTSPoint(zeros(ndim), zeros(ndim), zeros(ndim), -Inf)

function _nuts_copy!(dest::_NUTSPoint, src::_NUTSPoint)
    copyto!(dest.x, src.x); copyto!(dest.p, src.p); copyto!(dest.g, src.g)
    dest.logp = src.logp
    return dest
end

"""
    _leapfrog!(point, logdensity!, eps, metric, velocity, scratch)

One leapfrog step of size `eps`, in place. Returns the kinetic energy at the new
point, or `NaN` if the density there is not finite.
"""
function _leapfrog!(point::_NUTSPoint, logdensity!, eps::Float64,
    metric::CTSEMMetric, velocity::Vector{Float64}, scratch::Vector{Float64})
    n = length(point.x)
    @inbounds for i in 1:n
        point.p[i] += (eps / 2) * point.g[i]
    end
    _metric_velocity!(velocity, metric, point.p, scratch)
    @inbounds for i in 1:n
        point.x[i] += eps * velocity[i]
    end
    point.logp = logdensity!(point.g, point.x)
    isfinite(point.logp) || return NaN
    @inbounds for i in 1:n
        point.p[i] += (eps / 2) * point.g[i]
    end
    return _metric_velocity!(velocity, metric, point.p, scratch)
end

"""
    _uturn(pminus, pplus, xminus, xplus, metric, velocity, scratch)

The generalised no-U-turn criterion: has the trajectory started to double back?

Compares the position difference against the *velocities* at both ends, which is
the metric-aware form -- the plain momentum version is only correct for an
identity metric and quietly under-runs trajectories otherwise.
"""
function _uturn(pminus::Vector{Float64}, pplus::Vector{Float64},
    xminus::Vector{Float64}, xplus::Vector{Float64}, metric::CTSEMMetric,
    velocity::Vector{Float64}, scratch::Vector{Float64})
    _metric_velocity!(velocity, metric, pminus, scratch)
    back = 0.0
    @inbounds for i in eachindex(xplus)
        back += velocity[i] * (xplus[i] - xminus[i])
    end
    back < 0 && return true
    _metric_velocity!(velocity, metric, pplus, scratch)
    forward = 0.0
    @inbounds for i in eachindex(xplus)
        forward += velocity[i] * (xplus[i] - xminus[i])
    end
    return forward < 0
end

"""Workspace for one chain, so a transition allocates nothing."""
struct _NUTSWorkspace
    velocity::Vector{Float64}
    scratch::Vector{Float64}
    proposal::_NUTSPoint
    subproposal::_NUTSPoint
    minus::_NUTSPoint
    plus::_NUTSPoint
    working::_NUTSPoint
    # One candidate buffer *per depth*. A single shared one is a subtle
    # aliasing bug: the recursion at depth `d` holds a candidate while its own
    # children run, and a child writing to the same buffer overwrites it. The
    # symptom is not a crash but an over-dispersed posterior -- measured at
    # sd 1.03-1.10 and E[x^4] 3.7-4.9 on a standard normal, where 1 and 3 are
    # the answers.
    candidates::Vector{_NUTSPoint}
    # The *near* endpoint of each depth's subtree -- the state its first
    # leapfrog produced, not the state it leapfrogged from. The distinction is
    # the whole correctness of the stopping rule: the criterion has to be a
    # function of the subtree alone, and reaching one state back outside it
    # breaks the reversibility argument NUTS rests on. Measured cost of getting
    # this wrong on an isotropic 4-D normal: posterior sd 0.974 instead of 1,
    # and it does not shrink with more draws.
    xnear::Vector{Vector{Float64}}
    pnear::Vector{Vector{Float64}}
end

function _NUTSWorkspace(ndim::Int, maxdepth::Int)
    return _NUTSWorkspace(zeros(ndim), zeros(ndim),
        _nuts_point(ndim), _nuts_point(ndim), _nuts_point(ndim),
        _nuts_point(ndim), _nuts_point(ndim),
        [_nuts_point(ndim) for _ in 0:maxdepth],
        [zeros(ndim) for _ in 0:(maxdepth + 1)],
        [zeros(ndim) for _ in 0:(maxdepth + 1)])
end

"""
    _build_tree!(...)

One doubling of the trajectory, recursively.

Returns `(valid, logweight, nleapfrog, sumaccept, divergent)`. `logweight` is the
log of the summed multinomial weight over the subtree, which is what makes this
multinomial rather than slice NUTS: a proposal is drawn from the whole
trajectory in proportion to `exp(-H)`, so every state contributes rather than
only those above a slice.
"""
function _build_tree!(ws::_NUTSWorkspace, logdensity!, metric::CTSEMMetric,
    rng::AbstractRNG, depth::Int, eps::Float64, direction::Int, h0::Float64,
    maxdelta::Float64, tip::_NUTSPoint, out::_NUTSPoint)

    if depth == 0
        kinetic = _leapfrog!(tip, logdensity!, direction * eps, metric,
            ws.velocity, ws.scratch)
        if !isfinite(kinetic) || !isfinite(tip.logp)
            return (false, -Inf, 1, 0.0, true)
        end
        h = -tip.logp + kinetic
        divergent = !isfinite(h) || (h - h0) > maxdelta
        _nuts_copy!(out, tip)
        # A single-state subtree: it is its own near end.
        copyto!(ws.xnear[1], tip.x)
        copyto!(ws.pnear[1], tip.p)
        # `min(1, exp(h0 - h))` accumulated over the trajectory is Stan's
        # accept_stat, which is what dual averaging targets -- not the
        # acceptance of the final proposal, which for multinomial NUTS is 1.
        accept = min(1.0, exp(h0 - h))
        return (!divergent, h0 - h, 1, accept, divergent)
    end

    left = _build_tree!(ws, logdensity!, metric, rng, depth - 1, eps, direction,
        h0, maxdelta, tip, out)
    left[1] || return left
    # This subtree's near end is its left half's near end. Copied up now,
    # because the right half is about to reuse the level below.
    copyto!(ws.xnear[depth + 1], ws.xnear[depth])
    copyto!(ws.pnear[depth + 1], ws.pnear[depth])

    candidate = ws.candidates[depth]
    right = _build_tree!(ws, logdensity!, metric, rng, depth - 1, eps, direction,
        h0, maxdelta, tip, candidate)
    nleapfrog = left[3] + right[3]
    sumaccept = left[4] + right[4]
    divergent = left[5] || right[5]
    right[1] || return (false, left[2], nleapfrog, sumaccept, divergent)

    # Multinomial progressive sampling: take the right half's candidate with
    # probability equal to its share of the combined weight.
    logweight = _logaddexp(left[2], right[2])
    if log(rand(rng)) < right[2] - logweight
        _nuts_copy!(out, candidate)
    end

    # `tip` is the far end; `ws.xnear[depth + 1]` the near one. Which is
    # "minus" and which "plus" depends on the direction the trajectory was
    # extended in, because the criterion is stated in trajectory order.
    near_x = ws.xnear[depth + 1]
    near_p = ws.pnear[depth + 1]
    valid = if direction > 0
        !_uturn(near_p, tip.p, near_x, tip.x, metric, ws.velocity, ws.scratch)
    else
        !_uturn(tip.p, near_p, tip.x, near_x, metric, ws.velocity, ws.scratch)
    end
    return (valid, logweight, nleapfrog, sumaccept, divergent)
end

@inline _logaddexp(a::Float64, b::Float64) =
    a == -Inf ? b : (b == -Inf ? a : (a > b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b))))

"""
    _nuts_transition!(ws, logdensity!, metric, rng, x, g, logp, eps, maxdepth)

One NUTS transition from `x`, updating `x`, `g` and `logp` in place.
"""
function _nuts_transition!(ws::_NUTSWorkspace, logdensity!, metric::CTSEMMetric,
    rng::AbstractRNG, x::Vector{Float64}, g::Vector{Float64}, logp::Float64,
    eps::Float64, maxdepth::Int, maxdelta::Float64)

    _metric_momentum!(ws.minus.p, metric, rng)
    copyto!(ws.minus.x, x); copyto!(ws.minus.g, g); ws.minus.logp = logp
    _nuts_copy!(ws.plus, ws.minus)
    _nuts_copy!(ws.proposal, ws.minus)

    kinetic = _metric_velocity!(ws.velocity, metric, ws.minus.p, ws.scratch)
    h0 = -logp + kinetic
    logweight = 0.0            # log of exp(h0 - h0)
    nleapfrog = 0
    sumaccept = 0.0
    divergent = false
    depth = 0

    while depth < maxdepth
        direction = rand(rng, Bool) ? 1 : -1
        tip = direction > 0 ? ws.plus : ws.minus
        _nuts_copy!(ws.working, tip)
        result = _build_tree!(ws, logdensity!, metric, rng, depth, eps, direction,
            h0, maxdelta, ws.working, ws.subproposal)
        _nuts_copy!(tip, ws.working)
        nleapfrog += result[3]
        sumaccept += result[4]
        divergent |= result[5]
        result[1] || break

        # Between-tree multinomial step. Stan biases towards the *new* subtree
        # here, which is what makes the trajectory expand rather than dwell:
        # accept with probability min(1, w_new / w_old) instead of the share.
        if log(rand(rng)) < min(0.0, result[2] - logweight)
            _nuts_copy!(ws.proposal, ws.subproposal)
        end
        logweight = _logaddexp(logweight, result[2])
        depth += 1

        _uturn(ws.minus.p, ws.plus.p, ws.minus.x, ws.plus.x, metric,
            ws.velocity, ws.scratch) && break
    end

    copyto!(x, ws.proposal.x)
    copyto!(g, ws.proposal.g)
    energy = -ws.proposal.logp +
        _metric_velocity!(ws.velocity, metric, ws.proposal.p, ws.scratch)
    return (logp=ws.proposal.logp, accept=nleapfrog == 0 ? 0.0 : sumaccept / nleapfrog,
        divergent=divergent, depth=depth, nleapfrog=nleapfrog, energy=energy)
end
