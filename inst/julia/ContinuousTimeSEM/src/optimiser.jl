# The optimiser `ctsem_optimize` drives.
#
# Written here rather than taken from Optim because two of Optim's defaults,
# combined with what the engine passes it, were costing most of every fit, and
# because what comes next -- a batch that grows under the optimiser, and a
# Newton finish -- needs control no library interface offers.
#
# The two defaults, both measured on dev1 against the models in
# `dev/stochopt/models.R`:
#
# * `InitialStatic(alpha = 0.1, scaled = true)` caps EVERY iteration's trial step
#   at 0.1 raw units (LineSearches' initialguess.jl: alpha = min(0.1, |s|)/|s|),
#   not only the first. L-BFGS never tried its own unit step until it was
#   already within 0.1 of the optimum.
# * Supplying a preconditioner switches off the secant rescaling of the initial
#   inverse Hessian (Optim's l_bfgs.jl: `scaleinvH0 = P === nothing`), so H0 was
#   the bare transform metric, which does not grow with the data while the
#   curvature does. The line search then backtracked on most iterations: 863
#   objective calls for 280 iterations on a 1000-subject panel.
#
# Here the metric sets the SHAPE of H0 and the secant ratio sets its SCALE,
# `H0 = gamma D^-1` with `gamma = s'y / y'D^-1 y`, and only the first step (no
# curvature history yet) is shortened. The same panel: 25 iterations to a
# predicted gain of 0.1 against the engine's 280 to its stop.

using Random

"""
The state an iteration callback sees, shaped like the Optim state it replaced,
plus the predicted gain of the step that produced it where the step knows it
exactly (a Newton step does; an L-BFGS step leaves it to `CTSEMDirectional`).
"""
struct CTSEMIterate
    iteration::Int
    value::Float64
    g_norm::Float64
    gain::Float64
end
CTSEMIterate(iteration, value, g_norm) = CTSEMIterate(iteration, value, g_norm, NaN)

"""What `_ctsem_lbfgs` returns: the point and the reasons it stopped."""
struct CTSEMLBFGSResult
    minimizer::Vector{Float64}
    minimum::Float64
    gradient::Vector{Float64}
    iterations::Int
    f_calls::Int
    g_calls::Int
    g_converged::Bool
    f_converged::Bool
    x_converged::Bool
    linesearch_failed::Bool
    stopped_by_callback::Bool
    batch_sizes::Vector{Int}
    batch_iterations::Vector{Int}
end

"""Curvature pairs for the two-loop recursion, on the MINIMISED objective."""
mutable struct CTSEMLBFGSMemory
    S::Vector{Vector{Float64}}
    Y::Vector{Vector{Float64}}
    rho::Vector{Float64}
    m::Int
end
CTSEMLBFGSMemory(m::Integer) = CTSEMLBFGSMemory(Vector{Float64}[],
    Vector{Float64}[], Float64[], Int(m))

function _ctsem_lbfgs_push!(M::CTSEMLBFGSMemory, s::Vector{Float64},
        y::Vector{Float64})
    sy = dot(s, y)
    # A pair that is not a curvature (sy <= 0, possible after an Armijo step
    # with no curvature condition) would break positive definiteness; skip it.
    sy > 1e-12 * norm(s) * norm(y) || return false
    push!(M.S, s); push!(M.Y, y); push!(M.rho, 1 / sy)
    if length(M.S) > M.m
        popfirst!(M.S); popfirst!(M.Y); popfirst!(M.rho)
    end
    true
end

_ctsem_lbfgs_reset!(M::CTSEMLBFGSMemory) =
    (empty!(M.S); empty!(M.Y); empty!(M.rho); M)

# H * q for the inverse-Hessian approximation, with `dinv` the metric's inverse
# diagonal (all ones for no metric). With no pairs, returns `h0 * D^-1 q`.
function _ctsem_lbfgs_hmul(M::CTSEMLBFGSMemory, q::AbstractVector, h0::Float64,
        dinv::Vector{Float64})
    k = length(M.S)
    k == 0 && return h0 .* dinv .* q
    r = collect(Float64, q)
    a = Vector{Float64}(undef, k)
    @inbounds for i in k:-1:1
        a[i] = M.rho[i] * dot(M.S[i], r)
        r .-= a[i] .* M.Y[i]
    end
    y = M.Y[k]
    gamma = dot(M.S[k], y) / sum(i -> y[i]^2 * dinv[i], eachindex(y))
    r .= gamma .* dinv .* r
    @inbounds for i in 1:k
        b = M.rho[i] * dot(M.Y[i], r)
        r .+= (a[i] - b) .* M.S[i]
    end
    r
end

"""
    _ctsem_lbfgs(fg!, x0; ...)

L-BFGS with Armijo backtracking, minimising through `fg!(F, G, x)` -- the same
closure contract Optim's `only_fg!` used: `G === nothing` asks for the value
alone, `F === nothing` for the gradient alone, and an invalid trial point comes
back as a huge value with a zero gradient.

`callback(state::CTSEMIterate)` is called at iteration 0 and after every
accepted step, and stops the run by returning `true`. `directional` receives the
directional derivative `g's` of each step taken, which is what
`_ctsem_predicted_gain` reads.

`batch`, when given, is a `CTSEMBatch` whose `fg!` and `scores` replace the
full-data ones while it is smaller than the data; see `_ctsem_batch_step!`.
"""
function _ctsem_lbfgs(fg!, x0::AbstractVector; memory::Integer=20,
        metric=nothing, initial_alpha::Real=0.1, maxiter::Integer=1000,
        g_tol::Real=1e-8, f_tol::Real=0.0, x_tol::Real=0.0,
        callback=nothing, directional=nothing, batch=nothing,
        c1::Real=1e-4, maxbacktrack::Integer=40, iteration0::Integer=0)
    n = length(x0)
    x = collect(Float64, x0)
    dinv = metric === nothing ? ones(n) : begin
        d = collect(Float64, diag(metric))
        [isfinite(v) && v > 0 ? 1 / v : 1.0 for v in d]
    end
    M = CTSEMLBFGSMemory(memory)
    G = zeros(n)
    fcalls = 0; gcalls = 0
    sizes = Int[]; its = Int[]
    # The objective in force: the full one, or the current batch.
    evaluate!(F, Gout, y) = batch === nothing ? fg!(F, Gout, y) :
        _ctsem_batch_fg!(batch, F, Gout, y)
    f = evaluate!(0.0, G, x); fcalls += 1; gcalls += 1
    # Iteration 0 only on a run that is a run of its own; a continuation's
    # starting point is the last row its predecessor already recorded.
    stopped = callback !== nothing && iteration0 == 0 &&
        callback(CTSEMIterate(0, f, maximum(abs, G; init=0.0))) === true
    # Length `initial_alpha` in the metric's norm: a short first step when
    # there is no curvature history, measured so that a unit means the same
    # amount of model in every coordinate.
    metric_norm(v) = sqrt(sum(i -> v[i]^2 * dinv[i], eachindex(v)))
    h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
    iteration = 0
    gconv = maximum(abs, G; init=0.0) <= g_tol
    fconv = false; xconv = false; lsfail = false
    retried = false
    while !stopped && !gconv && iteration < maxiter
        s = -_ctsem_lbfgs_hmul(M, G, h0, dinv)
        # With no curvature pairs the step is the metric's alone, and at a start
        # where the transforms are flat the metric says a raw unit is worth
        # almost nothing -- so 0.1 in model units is an enormous raw step. One
        # raw unit at most, then: a long first step is how a fit once went from
        # raw 0 to raw 20.9 and ended where every transform is flat.
        if isempty(M.S)
            len = norm(s)
            len > 1 && (s .*= 1 / len)
        end
        dphi = dot(G, s)
        if !(dphi < 0) || !isfinite(dphi)
            # Not a descent direction: the memory has gone bad. Start it again.
            isempty(M.S) && (lsfail = true; break)
            _ctsem_lbfgs_reset!(M); h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
            continue
        end
        # The first trial carries the gradient too: most steps are accepted
        # whole, and then the iteration has cost one gradient and nothing else.
        alpha = 1.0
        xn = x .+ s
        Gn = similar(G)
        fn = evaluate!(0.0, Gn, xn); fcalls += 1; gcalls += 1
        have_gradient = true
        accepted = isfinite(fn) && fn <= f + c1 * alpha * dphi
        k = 0
        while !accepted && k < maxbacktrack
            k += 1
            # Quadratic interpolation of phi(alpha), safeguarded to [0.1, 0.5].
            denom = 2 * (fn - f - dphi * alpha)
            trial = isfinite(fn) && denom > 0 ? -dphi * alpha^2 / denom : 0.5alpha
            alpha = clamp(trial, 0.1alpha, 0.5alpha)
            xn = x .+ alpha .* s
            fn = evaluate!(0.0, nothing, xn); fcalls += 1
            have_gradient = false
            accepted = isfinite(fn) && fn <= f + c1 * alpha * dphi
        end
        if !accepted
            # A stale memory is the usual cause; drop it once, then give up.
            if !retried && !isempty(M.S)
                retried = true
                _ctsem_lbfgs_reset!(M)
                h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
                continue
            end
            lsfail = true
            break
        end
        retried = false
        if !have_gradient
            evaluate!(nothing, Gn, xn); gcalls += 1
        end
        directional === nothing || (directional.dphi0 = alpha * dphi)
        step = xn .- x
        iteration += 1
        fconv = abs(fn - f) <= f_tol * abs(fn)
        xconv = maximum(abs, step; init=0.0) <= x_tol
        _ctsem_lbfgs_push!(M, step, Gn .- G)
        x = xn; f = fn; G = Gn
        gconv = maximum(abs, G; init=0.0) <= g_tol
        stopped = callback !== nothing && callback(CTSEMIterate(
            iteration0 + iteration, f, maximum(abs, G; init=0.0))) === true
        # A growing batch changes the objective under the optimiser. The
        # curvature pairs are kept -- the batch objective is scaled to the full
        # data, so they estimate the same curvature -- but no pair spans the
        # change, and the value and gradient are re-read on the new batch. Once
        # it is the whole data the batch steps aside and the caller's own
        # objective is used from then on.
        if batch !== nothing && !stopped &&
                _ctsem_batch_step!(batch, x, G, q -> _ctsem_lbfgs_hmul(M, q, h0, dinv), iteration)
            if _ctsem_batch_full(batch)
                sizes = batch.sizes; its = batch.iterations
                batch = nothing
            end
            f = evaluate!(0.0, G, x); fcalls += 1; gcalls += 1
            gconv = maximum(abs, G; init=0.0) <= g_tol
        end
        (fconv && f_tol > 0) && break
        (xconv && x_tol > 0) && break
    end
    if batch !== nothing
        sizes = batch.sizes; its = batch.iterations
    end
    CTSEMLBFGSResult(x, f, G, iteration, fcalls, gcalls, gconv, fconv, xconv,
        lsfail, stopped, sizes, its)
end

# ------------------------------------------------------------------ batches
#
# Progressive batching over independent units. Units are permuted once and the
# batch is always a prefix of that permutation, so growing it ADDS units: within
# a stage the objective is deterministic and the line search is ordinary. The
# batch objective is scaled to the full data, (N/n) * loglik_batch, so curvature
# pairs from a small batch estimate the full curvature and the memory survives
# growth.
#
# Growth is decided in the optimiser's own metric and in objective units. With
# H the inverse-Hessian approximation and G the scaled batch gradient, a step
# predicts a gain of G'HG/2, and E[G'HG] = G'HG + tr(H Cov G). When the noise
# term reaches `theta` times the predicted gain, the batch cannot tell which way
# is up, and it grows by the factor that would bring the ratio back to `theta`
# (Byrd, Chin, Nocedal & Wu 2012's norm test, in this metric). A model with too
# few units to batch never starts one, so it is plain L-BFGS by construction.
#
# The prior is a global term and must not be scaled with the likelihood, so it
# is separated out: scaled = (N/n) (batch - prior) + prior, and likewise for the
# gradient. Both routes add exactly `_ctsem_log_prior` of the marginal
# objective to their value, and their score rows each carry 1/n of its
# gradient, which is what makes the separation exact. (ctFit's default
# `priors = 'randomCorr'` puts a prior on random-effect correlations even in a
# maximum likelihood fit, so refusing a prior would refuse every model with
# random effects.) A sampled TI predictor's imputation density is also global
# but is not attributable per subject, so that one still disables batching.
#
# In a batch stage the gradient comes from the score sweep itself -- the rows
# sum to it -- so the growth test costs nothing beyond the gradient. Taking
# the gradient from the ordinary trial and the rows from a second sweep, as a
# first version here did, doubled every batch iteration (tripled on laplace,
# where the score sweep costs two gradients) and made batching a loss.

mutable struct CTSEMBatch{O}
    full::O
    stage::Any
    perm::Vector{Int}
    m::Int
    nunits::Int
    nsubjects::Int
    theta::Float64
    fg!::Any              # (objective, F, G, x) -> F, the caller's trial path
    sizes::Vector{Int}
    iterations::Vector{Int}
    scores::Union{Nothing,Matrix{Float64}}   # rows at `scores_x`, from the last gradient
    scores_x::Vector{Float64}
end

_ctsem_nunits(o::CTSEMObjective) = length(o.subject_objectives)
_ctsem_nunits(o) = 0
_ctsem_batch_nsubjects(o::CTSEMObjective) = length(o.subject_objectives)

"""A CTSEMObjective over some of its subjects; the subject objectives are self-contained."""
_ctsem_subset_objective(o::CTSEMObjective, idx::AbstractVector{<:Integer}) =
    CTSEMObjective(o.params, o.subject_objectives[idx], nothing,
        o.prior_index, o.prior_scale, o.prior_weight,
        o.ti_missing_parameter, o.ti_missing_mu, o.ti_missing_sigma)
_ctsem_subset_objective(o, idx) = nothing

_ctsem_batchable(o::CTSEMObjective) = isempty(o.ti_missing_parameter)
_ctsem_batchable(o::CTSEMLaplaceObjective) = _ctsem_batchable(o.objective)
_ctsem_batchable(o) = false

# Laplace: whole UNITS (outer-level groups). Every evaluation path goes
# units.members[U] -> subject_objectives[i], so restricting `units` restricts
# the objective and the full subject list underneath is never touched for an
# excluded unit. Pairing the full spec with a subset *objective* instead would
# silently misassign groups: `_laplace_build_units` reads `group[i]` for the
# first n subjects without checking the length.
#
# Built field by field by name rather than positionally: the struct gains
# fields (the floor mode and its gate did, after this was first written, and a
# positional call then fails on every batched laplace fit). Per-unit and per-run
# state starts fresh, as the constructor starts it; every other field is the
# full objective's, shared. A field this does not know that looks per-unit is
# refused rather than copied at the wrong length.
function _ctsem_subset_objective(L::CTSEMLaplaceObjective, keep::AbstractVector{<:Integer})
    u = L.units
    units = CTSEMLaplaceUnits(u.members[keep], u.offsets[keep], u.dims[keep],
        u.blocks[keep])
    n = length(keep)
    NU = length(u.members)
    fresh = Dict{Symbol,Any}(
        :units => units,
        :modes => [zeros(units.dims[U]) for U in 1:n],
        :workspaces => [Dict{Any,Any}() for _ in eachindex(L.workspaces)],
        :inner_iterations => zeros(Int, n), :inner_gradient => zeros(n),
        :inner_converged => falses(n), :hessian_repaired => falses(n),
        :mode_repaired => falses(n), :logdet_floored => falses(n),
        # per-run: where the objective was last evaluated, and a counter
        :last_values => Float64[], :gated_units => 0)
    args = map(fieldnames(typeof(L))) do name
        haskey(fresh, name) && return fresh[name]
        v = getfield(L, name)
        if v isa AbstractVector && NU > 1 && length(v) == NU
            error("_ctsem_subset_objective: CTSEMLaplaceObjective field `$(name)` ",
                "has one entry per unit and is not known here; add it to `fresh`.")
        end
        v
    end
    typeof(L)(args...)
end
_ctsem_nunits(o::CTSEMLaplaceObjective) = length(o.units.members)
_ctsem_prior_objective(o::CTSEMObjective) = o
_ctsem_prior_objective(o::CTSEMLaplaceObjective) = o.objective
_ctsem_batch_nsubjects(o::CTSEMLaplaceObjective) = sum(length, o.units.members; init=0)

"""
    _ctsem_batch_plan(objective; theta, n0, rng)

A batch for `objective`, or `nothing` when it cannot or should not batch: a
route with no subset, a prior or sampled TI predictor, or too few units -- a
first batch of at least 20 units and a quarter of the data at most.
"""
function _ctsem_batch_plan(objective, fg!; theta::Real=0.25, n0::Integer=0,
        rng=Random.MersenneTwister(20260923))
    _ctsem_batchable(objective) || return nothing
    NU = _ctsem_nunits(objective)
    m = n0 > 0 ? Int(n0) : max(20, cld(NU, 32))
    m > NU ÷ 4 && return nothing
    perm = randperm(rng, NU)
    stage = _ctsem_subset_objective(objective, sort(perm[1:m]))
    stage === nothing && return nothing
    CTSEMBatch(objective, stage, perm, m, NU, _ctsem_batch_nsubjects(objective),
        Float64(theta), fg!, [m], [0], nothing, Float64[])
end

_ctsem_batch_full(b::CTSEMBatch) = b.m >= b.nunits
_ctsem_batch_scale(b::CTSEMBatch) = b.nsubjects / _ctsem_batch_nsubjects(b.stage)

# The stage objective on the minimised scale, with its likelihood scaled to the
# full data and its prior left whole. A value alone goes through the caller's
# trial, which decides validity; a gradient comes from the score sweep, whose
# rows the growth test then reads at the same point.
const _CTSEM_BATCH_INVALID = floatmax(Float64) / 1e8

function _ctsem_batch_fg!(b::CTSEMBatch, F, G, x)
    c = _ctsem_batch_scale(b)
    po = _ctsem_prior_objective(b.full)
    xv = collect(Float64, x)
    prior = _ctsem_log_prior(po, xv)
    if G === nothing
        v = b.fg!(b.stage, F, nothing, x)
        v === nothing && return nothing
        v >= _CTSEM_BATCH_INVALID && return v
        return -(c * (-v - prior) + prior)
    end
    r = try
        ctsem_subject_gradients(b.stage, xv)
    catch err
        # A point the model cannot evaluate is no data; code that is wrong is
        # an error, and an interrupt stops the fit (`_ctsem_must_propagate`).
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    if r === nothing || !isfinite(r.value) || !all(isfinite, r.scores)
        fill!(G, 0.0); b.scores = nothing
        return F === nothing ? nothing : _CTSEM_BATCH_INVALID
    end
    S = Matrix{Float64}(r.scores)
    gprior = _ctsem_log_prior_gradient!(zeros(length(xv)), po, xv)
    g = vec(sum(S; dims=1))
    G .= -(c .* (g .- gprior) .+ gprior)
    b.scores = S; b.scores_x = xv
    F === nothing ? nothing : -(c * (r.value - prior) + prior)
end

"""
    _ctsem_batch_step!(batch, x, G, hmul, iteration)

Apply the growth test at `x`, where `G` is the scaled batch gradient (of the
minimised objective) and `hmul` applies the current inverse-Hessian
approximation. Grows the batch and returns `true` when the noise of the
predicted gain exceeds `theta` times the gain; `false` otherwise.
"""
function _ctsem_batch_step!(b::CTSEMBatch, x, G, hmul, iteration)
    _ctsem_batch_full(b) && return false
    # The rows from the gradient just taken at `x`; a fresh sweep only if the
    # last gradient was somewhere else (it is not, on the path through
    # `_ctsem_lbfgs`, which always ends an iteration with a gradient at `x`).
    S = b.scores !== nothing && b.scores_x == x ? b.scores : try
        Matrix{Float64}(ctsem_subject_gradients(b.stage, collect(Float64, x)).scores)
    catch err
        # A point the model cannot evaluate is no data; code that is wrong is
        # an error, and an interrupt stops the fit (`_ctsem_must_propagate`).
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    grow = false
    newm = b.m
    if S === nothing || !all(isfinite, S)
        # No spread to read: the batch is not telling us anything reliable.
        grow = true; newm = min(b.nunits, 2b.m)
    else
        n = size(S, 1)
        d = hmul(G)
        gain = 0.5 * dot(G, d)
        Sbar = vec(sum(S; dims=1)) ./ n
        acc = 0.0
        for i in 1:n
            di = S[i, :] .- Sbar
            acc += dot(di, hmul(di))
        end
        V = acc / max(n - 1, 1)
        c = _ctsem_batch_scale(b)
        noise = 0.5 * c^2 * n * (1 - b.m / b.nunits) * V
        if !(gain > 0) || noise > b.theta * gain
            grow = true
            factor = gain > 0 ? noise / (b.theta * gain) : 2.0
            newm = min(b.nunits, max(2b.m, ceil(Int, b.m * min(factor, 1e6))))
        end
    end
    grow || return false
    b.m = newm
    b.stage = newm >= b.nunits ? b.full :
        _ctsem_subset_objective(b.full, sort(b.perm[1:newm]))
    push!(b.sizes, newm); push!(b.iterations, iteration)
    true
end

# ------------------------------------------------------------------ Newton finish
#
# Once L-BFGS has come close -- its predicted gain below `switch` -- a damped
# Newton iteration on the exact Hessian finishes in a handful of steps what
# L-BFGS would take tens or hundreds to polish, and the exact Hessian at the
# final point is the one the certification needs anyway, so it is returned for
# reuse rather than computed again. Only where a Hessian is cheap: the marginal
# route's forward-over-adjoint Hessian is a few gradients (6 on a 1000-subject
# panel, dev1); the laplace route's is finite differences, 2p gradients, and
# there L-BFGS to the end is faster (measured: 0.6x with the finish, 1.3-1.9x
# without).

_ctsem_cheap_hessian(::CTSEMObjective) = true
_ctsem_cheap_hessian(::Any) = false

"""
    _ctsem_newton_finish(objective, x, fg!; tol, maxit, curvature, contraction)

Damped Newton from `x` on the minimised objective behind `fg!`. What the
steps are taken against is `curvature`:

- `:exact`  the exact Hessian, refreshed whenever the predicted gain is not
  contracting by `contraction` per step;
- `:chord`  the exact Hessian at `x`, kept for every step (the chord, or
  simplified Newton, method: linear convergence at the rate the Hessian
  changes between `x` and the optimum, which from a hand-over within ~0.1
  nats is fast);
- `:subset` the likelihood Hessian of a random `subset` share of the units
  (at least `subset_min`), scaled up to the data, with the prior's curvature
  kept whole -- a chord Hessian at a fraction of the cost.

Whatever the steps used, a failed line search is answered with the exact
Hessian, and the run always ends with the exact Hessian at the final point,
which is the certification's matrix and is returned for it. Returns the point,
its value and gradient, that Hessian of the MAXIMISED objective (or `nothing`),
the steps taken, the predicted gain, and how many full and subset Hessians
were formed.
"""
function _ctsem_newton_finish(objective, x0, f0, G0, fg!; tol::Real=1e-8,
        maxit::Integer=30, contraction::Real=0.25, callback=nothing,
        iteration0::Integer=0, curvature::Symbol=:exact,
        subset::Real=0.125, subset_min::Integer=200,
        rng=Random.MersenneTwister(20260924))
    curvature in (:exact, :chord, :subset) || throw(ArgumentError(
        "newton curvature must be exact, chord or subset, got $(curvature)"))
    x = collect(Float64, x0); f = f0; G = collect(Float64, G0)
    full_hessians = 0; subset_hessians = 0
    hessof(o, y) = try
        H = Matrix{Float64}(ctsem_hessian(o, y))
        all(isfinite, H) ? H : nothing
    catch err
        # A point the model cannot evaluate is no data; code that is wrong is
        # an error, and an interrupt stops the fit (`_ctsem_must_propagate`).
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    # Of the minimised objective, as every step below expects.
    # `local`: an assignment to `H` in a closure would rebind the H below.
    hess(y) = (full_hessians += 1; local Hx = hessof(objective, y);
        Hx === nothing ? nothing : -Hx)
    NU = _ctsem_nunits(objective)
    m = min(NU, max(Int(subset_min), ceil(Int, subset * NU)))
    sub = (curvature === :subset && 2m <= NU) ?
        _ctsem_subset_objective(objective, sort(randperm(rng, NU)[1:m])) : nothing
    function subhess(y)
        sub === nothing && return hess(y)
        subset_hessians += 1
        Hs = hessof(sub, y)
        Hs === nothing && return hess(y)
        # The prior's curvature is diagonal and global: take it out of the
        # subset's, scale the likelihood part, and put it back once.
        po = _ctsem_prior_objective(objective)
        Hp = zeros(length(y))
        for k in eachindex(po.prior_index)
            Hp[po.prior_index[k]] -= po.prior_weight / po.prior_scale[k]^2
        end
        local c = _ctsem_batch_nsubjects(objective) / _ctsem_batch_nsubjects(sub)
        -(c .* (Hs .- Diagonal(Hp)) .+ Diagonal(Hp))
    end
    step_hessian(y) = curvature === :subset ? subhess(y) : hess(y)
    H = step_hessian(x)
    H === nothing && return (x=x, f=f, G=G, hessian=nothing, steps=0,
        gain=Inf, full_hessians=full_hessians, subset_hessians=subset_hessians,
        fcalls=0, gcalls=0)
    exact_at_x = curvature !== :subset
    fcalls = 0; gcalls = 0
    steps = 0; mu = 0.0; prevgain = Inf; gain = Inf
    # The gain the undamped step predicts, flat directions included at their
    # floored curvature. Judging convergence over the trusted directions alone
    # stopped the finish at the start of a nearly flat ray: on a fixture with a
    # diffusion correlation flat below raw -6, 3.5e-6 nats short of the point a
    # profile then found, so the estimate was not the maximum. A truly flat
    # direction has no gradient and adds nothing; a nearly flat one is walked
    # while the step still promises more than `tol`, as L-BFGS used to.
    function newton(H, G, mu)
        E = eigen(Symmetric(H))
        lmax = maximum(abs, E.values; init=0.0)
        lmax > 0 || return nothing
        floor = 1e-8 * lmax
        c = E.vectors' * G
        floored = max.(E.values, floor)
        gain = 0.5 * sum(abs2.(c) ./ floored; init=0.0)
        lam = floored .+ mu * lmax
        (step=-(E.vectors * (c ./ lam)), gain=gain)
    end
    at_x = exact_at_x    # whether H is the exact Hessian at the current x
    while steps < maxit
        nt = newton(H, G, mu)
        nt === nothing && break
        gain = nt.gain
        gain < tol && break
        dphi = dot(G, nt.step)
        alpha = 1.0
        xn = x .+ nt.step
        fn = fg!(0.0, nothing, xn); fcalls += 1
        k = 0
        while !(isfinite(fn) && fn <= f + 1e-4 * alpha * dphi) && k < 30
            k += 1; alpha /= 2
            xn = x .+ alpha .* nt.step
            fn = fg!(0.0, nothing, xn); fcalls += 1
        end
        if !(isfinite(fn) && fn <= f + 1e-4 * alpha * dphi)
            if !at_x
                H = hess(x); at_x = true
                H === nothing && break
            else
                mu = mu == 0 ? 1e-4 : 10mu
                mu > 1e2 && break
            end
            continue
        end
        Gn = similar(G)
        fg!(nothing, Gn, xn); gcalls += 1
        x = xn; f = fn; G = Gn; steps += 1; at_x = false
        mu = alpha == 1 ? mu / 10 : mu
        mu < 1e-8 && (mu = 0.0)
        callback === nothing || callback(CTSEMIterate(iteration0 + steps, f,
            maximum(abs, G; init=0.0), gain))
        # Only the exact variant refreshes on slow contraction; the chord and
        # the subset keep their matrix, which is the point of them. A step
        # that fails outright still gets the exact Hessian (above).
        if curvature === :exact && gain / prevgain > contraction && steps > 1
            H = hess(x); at_x = true
            H === nothing && break
            prevgain = Inf
        else
            prevgain = gain
        end
    end
    # The certification's Hessian: at the final point, exactly. A kept or a
    # subset Hessian judges the gain against itself, so the loop above can
    # stop where the exact curvature still predicts more than `tol`; then the
    # finish carries on with the exact Hessian until it agrees. Otherwise the
    # certification would find the point unfinished and resume L-BFGS with its
    # stopping rule switched off -- measured: 1045 iterations to the cap on a
    # 30-subject model the exact finish closes in 10 steps.
    for _ in 1:5
        if !at_x
            H = hess(x); at_x = true
        end
        H === nothing && break
        nt = newton(H, G, 0.0)
        gain = nt === nothing ? Inf : nt.gain
        (nt === nothing || gain < tol || steps >= maxit + 5) && break
        dphi = dot(G, nt.step)
        alpha = 1.0
        xn = x .+ nt.step
        fn = fg!(0.0, nothing, xn); fcalls += 1
        k = 0
        while !(isfinite(fn) && fn <= f + 1e-4 * alpha * dphi) && k < 30
            k += 1; alpha /= 2
            xn = x .+ alpha .* nt.step
            fn = fg!(0.0, nothing, xn); fcalls += 1
        end
        (isfinite(fn) && fn <= f + 1e-4 * alpha * dphi) || break
        Gn = similar(G)
        fg!(nothing, Gn, xn); gcalls += 1
        x = xn; f = fn; G = Gn; steps += 1; at_x = false
        callback === nothing || callback(CTSEMIterate(iteration0 + steps, f,
            maximum(abs, G; init=0.0), gain))
    end
    (x=x, f=f, G=G, hessian=H === nothing ? nothing : -H, steps=steps,
     gain=gain, full_hessians=full_hessians, subset_hessians=subset_hessians,
     fcalls=fcalls, gcalls=gcalls)
end
