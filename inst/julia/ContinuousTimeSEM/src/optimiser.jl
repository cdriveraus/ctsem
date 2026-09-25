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
route with no subset, a sampled TI predictor (the one global term that is not
attributable to a subject -- a prior is separated out, see above, and does not
stop batching), or too few units -- a first batch of at least 20 units and a
quarter of the data at most.

Provenance: `theta = 0.25` and the first batch of `max(20, NU/32)` units were
set when this optimiser replaced Optim's (92499af5, 2026-09-24), on the seven
regimes of `dev/stochopt/models.R` -- panel, panel5k, long, ordinal Laplace,
nonlin, small, bigp -- and have not been swept since (the consolidation plan's
Appendix B). Batching changes the path to the optimum, not the optimum.
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

# ------------------------------------------------------------------ the endgame
#
# Once L-BFGS has come close, a damped Newton iteration finishes in a handful of
# steps what L-BFGS would take tens or hundreds to polish, and the curvature it
# ends on is the one the certification needs. So everything the certification
# reads is produced here and handed back -- the final Hessian and where it was
# evaluated, the escape from a saddle, and what the flat directions are worth --
# and R reads its verdict off those numbers (`.ctBackendCorrectResult()`) without
# another engine call. Until 2026-09-25 R took its own damped Newton step and its
# own negative-curvature step, on a Hessian it computed again every round: one
# computation in two languages, and on a Laplace fit `2 npar` gradients a round.
#
# What the steps are taken against depends on what a Hessian costs:
#
# * the marginal route's is forward-over-adjoint, a few gradients (6 on a
#   1000-subject panel, dev1), so L-BFGS hands over early and every step uses
#   the exact Hessian, refreshed while the steps contract slowly;
# * the Laplace route's is central differences of the gradient, `2 npar`
#   gradients, and an exact finish there made fits slower (0.6x with it, 1.3 to
#   1.9x without, dev1). So L-BFGS runs to its own stopping rule there, and the
#   finish takes one Hessian at the hand-over, reuses it for its steps (the
#   chord) and keeps it as the final one unless the steps moved the estimate by
#   more than `_CTSEM_HESSIAN_REUSE_SE` standard errors: one Hessian per fit in
#   the common case.
#
# Provenance: the hand-over at a predicted gain of 0.1, the cap of 30 steps and
# the contraction of 0.25 that triggers a refresh were set when this optimiser
# replaced Optim's (92499af5, 2026-09-24), on the seven regimes of
# `dev/stochopt/models.R` (panel, panel5k, long, ordinal Laplace, nonlin, small,
# bigp), and have not been swept since. The eigenvalue floor at 1e-8 of the
# largest and the Levenberg start at 1e-4 are the usual safeguards, not tuned.
# The consolidation plan's Appendix B has the list.

_ctsem_cheap_hessian(::CTSEMObjective) = true
_ctsem_cheap_hessian(::Any) = false

"""
The route's own curvature for the finish: `:exact` where a Hessian costs a few
gradients, `:chord` where it costs `2 npar` of them, and `nothing` for a route
the finish does not run on -- a pinned profile point, whose wrapper has no
Hessian of its own, and the state-explicit target, which is never certified.
"""
_ctsem_finish_curvature(::CTSEMObjective) = :exact
_ctsem_finish_curvature(::CTSEMLaplaceObjective) = :chord
_ctsem_finish_curvature(::Any) = nothing

"""
How far the estimate may be from where a Hessian was evaluated, in the standard
errors that Hessian implies, for that Hessian still to count as the curvature at
the estimate. A displacement of a hundredth of a standard error moves the
curvature by a hundredth of what it changes over a whole standard error, which
in a regular model is itself small. Chosen, not measured (decision 3 of
`review/OPTIM-consolidation-plan-2026-09-25.md`); the optimiser bench is what
says whether it matters. R's `.ctBackendHessianReuse()` is the same number and
is what the fit passes in.
"""
const _CTSEM_HESSIAN_REUSE_SE = 0.01

"""
The relative curvature at or below which a direction is flat: excluded from the
certification's gap and probed instead. R's `.ctFlatDirectionRtol()`, which the
fit passes in, so the certification, the intervals and the probe share one rule.
"""
const _CTSEM_FLAT_RTOL = 1e-12

"""
How negative a direction's curvature has to be, relative to the largest, for
the point to count as a saddle: the `negative` of R's
`.ctBackendInformationSplit()`.
"""
const _CTSEM_NEGATIVE_RTOL = 1e-8

"""
The lengths, in raw units, the finish tries along the most negative curvature at
a saddle, both signs, the gradient's first. They stop at 3 because every ctsem
transform is nearly flat a few units from its centre, and a longer rung measures
the plateau rather than the likelihood. Set on the rank-deficient mvmix fixture
(`test-julia-multivariate-mixed.R`), where L-BFGS resumed from the saddle crept
0.02 nats in 1000 iterations along this direction while a maximum 2 nats higher
sat in the other basin (8ad83c19); moved here from the R-side
negative-curvature step on 2026-09-25.
"""
const _CTSEM_SADDLE_LADDER = (0.03, 0.1, 0.3, 1.0, 3.0)

"""
The lengths, in raw units, the flat-direction probe steps along the gradient's
component outside the trusted directions. Raw units are standardised by
construction -- the priors are normal(0, 1) and each transform carries its own
scale -- so a quarter, one and four span a small move to a large one without
knowing the model. Chosen by that argument when the certification came in
(18a1872b), not measured; moved here from the R-side probe on 2026-09-25.
"""
const _CTSEM_FLAT_PROBE_LENGTHS = (0.25, 1.0, 4.0)

# What a probe that had nothing to probe reports. A function rather than a
# constant, so no caller shares a mutable vector.
_ctsem_no_probe() = (gain=0.0, length=0.0, longest=0.0, direction=Float64[],
    evaluations=0)

"""
    _ctsem_information_split(H; rtol, negative)

The eigen-decomposition of `H`, the Hessian of the MINIMISED objective -- the
observed information -- with its trusted and its negative directions marked by
the rule R's `.ctBackendInformationSplit()` applies, so the saddle test and the
flat probe here see the directions the certification will. `nothing` when it
cannot be decomposed.
"""
function _ctsem_information_split(H::AbstractMatrix; rtol::Real=_CTSEM_FLAT_RTOL,
        negative::Real=_CTSEM_NEGATIVE_RTOL)
    (size(H, 1) == size(H, 2) && size(H, 1) > 0 && all(isfinite, H)) ||
        return nothing
    E = try
        eigen(Symmetric((H .+ transpose(H)) ./ 2))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        return nothing
    end
    values = E.values
    scale = maximum(values)
    if !(isfinite(scale) && scale > 0)
        # Not a maximum in any direction, so nothing is trusted, and saying so
        # is the finding.
        scale = maximum(abs, values)
        (isfinite(scale) && scale > 0) || return nothing
        return (values=values, vectors=E.vectors, scale=scale,
            trusted=falses(length(values)), negative=values .< -negative * scale)
    end
    (values=values, vectors=E.vectors, scale=scale,
        trusted=values .> rtol * scale, negative=values .< -negative * scale)
end

"""
    _ctsem_hessian_distance(split, from, to)

How far `to` is from `from`, in the standard errors the curvature `split`
describes: the largest `abs(to[i] - from[i]) / se[i]`, with `se` the square
root of the diagonal of the inverse information over its trusted directions --
the covariance a fit reports. A coordinate the trusted curvature says nothing
about has a standard error of zero there, so moving it at all is infinitely far
and the Hessian is not reused, which is the safe side. R's
`.ctBackendHessianDistance()` is the same arithmetic.
"""
function _ctsem_hessian_distance(split, from::AbstractVector, to::AbstractVector)
    d = collect(Float64, to) .- collect(Float64, from)
    all(iszero, d) && return 0.0
    split === nothing && return Inf
    keep = findall(split.trusted)
    worst = 0.0
    for i in eachindex(d)
        iszero(d[i]) && continue
        variance = 0.0
        for k in keep
            variance += split.vectors[i, k]^2 / split.values[k]
        end
        worst = max(worst, variance > 0 ? abs(d[i]) / sqrt(variance) : Inf)
    end
    worst
end

"""
    _ctsem_flat_residual(split, ascent)

The part of the ascent gradient outside the trusted directions -- what R's
`.ctBackendOptimGap()` calls the residual. The gap cannot see it, so it is
probed.
"""
function _ctsem_flat_residual(split, ascent::AbstractVector)
    V = split.vectors[:, split.trusted]
    collect(Float64, ascent .- V * (transpose(V) * ascent))
end

"""
    _ctsem_flat_probe(value_at, x, value, direction; lengths)

What the directions the curvature does not trust are worth, measured rather
than predicted: step along `direction` by each of `lengths` raw units and keep
the best improvement on `value`. A norm of the gradient there cannot say -- it
has units -- and in a flat direction the quadratic model that would turn it into
a likelihood is exactly what does not hold.

`value_at(y)` is the MAXIMISED objective at `y`, or `-Inf` where the route
refuses the point: on the Laplace route, a point where a unit's inner solve did
not converge, which R's probe could not see. Returns the gain, the length that
gave it, the longest length that could be evaluated at all, the unit direction
(empty when there was nothing to probe) and the evaluations spent. The logic of
the R-side probe it replaced.
"""
function _ctsem_flat_probe(value_at, x::AbstractVector, value::Real,
        direction::AbstractVector; lengths=_CTSEM_FLAT_PROBE_LENGTHS)
    magnitude = norm(direction)
    (isfinite(magnitude) && magnitude > 0) || return _ctsem_no_probe()
    unit = collect(Float64, direction) ./ magnitude
    best = 0.0; at = 0.0; longest = 0.0; evaluations = 0
    for len in lengths
        got = value_at(x .+ len .* unit)
        evaluations += 1
        isfinite(got) || continue
        longest = max(longest, Float64(len))
        if got - value > best
            best = got - value
            at = Float64(len)
        end
    end
    (gain=best, length=at, longest=longest, direction=unit, evaluations=evaluations)
end

"""
    _ctsem_saddle_ladder(value_at, x, value, split, ascent; ladder)

At a saddle, the best point along the most negative curvature: the unit
eigenvector, both signs -- the gradient's first, since along a saddle direction
the gradient says which way is up -- at each rung of `ladder` raw units. The
Newton step cannot find this: it lives in the trusted directions, where the
point is already a maximum. `value_at` is as for `_ctsem_flat_probe`. Returns
the best improvement on `value` as `(point, value, length)`, or `nothing`, and
the evaluations spent. The logic of the R-side negative-curvature step it
replaced.
"""
function _ctsem_saddle_ladder(value_at, x::AbstractVector, value::Real, split,
        ascent::AbstractVector; ladder=_CTSEM_SADDLE_LADDER)
    (split === nothing || !any(split.negative)) &&
        return (best=nothing, evaluations=0)
    v = split.vectors[:, argmin(split.values)]
    first = dot(ascent, v) >= 0 ? 1.0 : -1.0
    best = nothing
    evaluations = 0
    for side in (first, -first), len in ladder
        point = x .+ (side * len) .* v
        got = value_at(point)
        evaluations += 1
        if isfinite(got) && got > value && (best === nothing || got > best.value)
            best = (point=point, value=got, length=side * len)
        end
    end
    (best=best, evaluations=evaluations)
end

"""
    _ctsem_newton_finish(objective, x, f, G, fg!; ...)

The endgame, from `x` on the minimised objective behind `fg!` (`f` and `G` its
value and gradient there): damped Newton steps, the final Hessian, an escape
from a saddle, and the flat-direction probe.

What the steps are taken against is `curvature`:

- `:exact`  the exact Hessian, refreshed whenever the predicted gain is not
  contracting by `contraction` per step;
- `:chord`  the exact Hessian at `x`, kept for every step (the chord, or
  simplified Newton, method: linear convergence at the rate the Hessian
  changes between `x` and the optimum, which from a hand-over this close is
  fast);
- `:subset` the likelihood Hessian of a random `subset` share of the units
  (at least `subset_min`), scaled up to the data, with the prior's curvature
  kept whole -- a chord Hessian at a fraction of the cost.

The final Hessian is the exact one at the final point, with one exception: a
chord Hessian whose steps converged within `reuse_se` standard errors of where
it was evaluated is kept (`hessian_at` and `distance` say where, and how far
that is -- `_ctsem_hessian_distance`). Whatever the steps used, a failed line
search is answered with the exact Hessian, and a subset Hessian, or a chord one
the steps moved too far from, is replaced by the exact Hessian at the final
point and the steps continued until it agrees.

At a saddle -- negative curvature in the final Hessian -- the ladder along the
most negative curvature is tried (`_ctsem_saddle_ladder`), and the finish
continues from the better point with the exact Hessian there, at most
`max_escapes` times: it never returns a point with a direction of negative
curvature it did not try, short of that cap. Then, when `probe`, the
flat-direction probe (`_ctsem_flat_probe`) at the final point.

`take_steps = false` gives the certification's numbers alone: the Hessian at `x`
and the probe, with no step, no refresh and no ladder, for a point the optimiser
left without a finish (`ctsem_endgame`).

`value_at(y)`, the maximised objective or `-Inf` for the ladder and the probe,
defaults to the route's own predicate (`_ctsem_probe_value`); `fg!` is the
optimiser's trial path, for the steps.

Returns the point and its value and gradient (minimised); the Hessian of the
MAXIMISED objective, or `nothing`, with `hessian_at` and `distance`; the steps
taken, including escapes, and the predicted gain at the end; the full and
subset Hessians formed and the calls made; the saddle record (`escapes`,
`saddle`, `ladder_tried`, `ladder_gain`); the probe record; and the history of
the steps, one entry each (`kind` is "newton", "exact" or "saddle").
"""
function _ctsem_newton_finish(objective, x0, f0, G0, fg!; tol::Real=1e-8,
        maxit::Integer=30, contraction::Real=0.25, callback=nothing,
        iteration0::Integer=0, curvature::Symbol=:exact,
        subset::Real=0.125, subset_min::Integer=200,
        rng=Random.MersenneTwister(20260924), take_steps::Bool=true,
        probe::Bool=true, reuse_se::Real=_CTSEM_HESSIAN_REUSE_SE,
        flat_rtol::Real=_CTSEM_FLAT_RTOL, negative::Real=_CTSEM_NEGATIVE_RTOL,
        ladder=_CTSEM_SADDLE_LADDER, probe_lengths=_CTSEM_FLAT_PROBE_LENGTHS,
        max_escapes::Integer=3, value_at=nothing)
    curvature in (:exact, :chord, :subset) || throw(ArgumentError(
        "newton curvature must be exact, chord or subset, got $(curvature)"))
    valueof = value_at === nothing ? (y -> _ctsem_probe_value(objective, y)) :
        value_at
    x = collect(Float64, x0); f = Float64(f0); G = collect(Float64, G0)
    full_hessians = 0; subset_hessians = 0
    hessof(o, y) = try
        local Hy = Matrix{Float64}(ctsem_hessian(o, y))
        all(isfinite, Hy) ? Hy : nothing
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
        local Hs = hessof(sub, y)
        Hs === nothing && return hess(y)
        # The prior's curvature is diagonal and global: take it out of the
        # subset's, scale the likelihood part, and put it back once.
        local po = _ctsem_prior_objective(objective)
        local Hp = zeros(length(y))
        for k in eachindex(po.prior_index)
            Hp[po.prior_index[k]] -= po.prior_weight / po.prior_scale[k]^2
        end
        local c = _ctsem_batch_nsubjects(objective) / _ctsem_batch_nsubjects(sub)
        -(c .* (Hs .- Diagonal(Hp)) .+ Diagonal(Hp))
    end
    step_hessian(y) = curvature === :subset ? subhess(y) : hess(y)
    history = (kind=String[], gain=Float64[], alpha=Float64[], value=Float64[])
    remember!(kind, g, a) = (push!(history.kind, kind); push!(history.gain, g);
        push!(history.alpha, a); push!(history.value, -f); nothing)
    fcalls = 0; gcalls = 0
    # One Armijo backtracking search along `step` from the current point. The
    # names inside are `local` so they cannot rebind the finish's own.
    #
    # Accepted only on a decrease the objective can represent, and the halving
    # stops once the first-order gain falls below that same resolution: Armijo
    # alone scales with the step, so at a small enough `alpha` an increase of
    # 1e-14 satisfies it, and accepting that spends a step on a point no
    # different from the one it left. The rule R's damped step applied before
    # the steps moved here, now applied to every step.
    function search(step, dphi)
        local floor = max(abs(f), 1.0) * eps()
        local alpha = 1.0
        local xn = x .+ step
        local fn = fg!(0.0, nothing, xn)
        fcalls += 1
        local k = 0
        while !(isfinite(fn) && fn <= f + 1e-4 * alpha * dphi && f - fn > floor) &&
                k < 30 && alpha * abs(dphi) > floor
            k += 1; alpha /= 2
            xn = x .+ alpha .* step
            fn = fg!(0.0, nothing, xn)
            fcalls += 1
        end
        (ok=isfinite(fn) && fn <= f + 1e-4 * alpha * dphi && f - fn > floor,
            x=xn, f=fn, alpha=alpha)
    end
    # The gain the undamped step predicts, flat directions included at their
    # floored curvature. Judging convergence over the trusted directions alone
    # stopped the finish at the start of a nearly flat ray: on a fixture with a
    # diffusion correlation flat below raw -6, 3.5e-6 nats short of the point a
    # profile then found, so the estimate was not the maximum. A truly flat
    # direction has no gradient and adds nothing; a nearly flat one is walked
    # while the step still promises more than `tol`, as L-BFGS used to.
    function newton(Hm, Gv, mu)
        local E = eigen(Symmetric(Hm))
        local lmax = maximum(abs, E.values; init=0.0)
        lmax > 0 || return nothing
        local floored = max.(E.values, 1e-8 * lmax)
        local c = E.vectors' * Gv
        local g = 0.5 * sum(abs2.(c) ./ floored; init=0.0)
        local lam = floored .+ mu * lmax
        (step=-(E.vectors * (c ./ lam)), gain=g)
    end
    report(g) = callback === nothing || callback(CTSEMIterate(iteration0 + steps,
        f, maximum(abs, G; init=0.0), g))
    H = step_hessian(x)
    hat = copy(x)                       # where `H` was evaluated
    exact = curvature !== :subset       # `H` is the exact Hessian at `hat`
    steps = 0; escapes = 0; gain = Inf
    saddle = false; ladder_tried = false; ladder_gain = 0.0
    if H === nothing
        return (x=x, f=f, G=G, hessian=nothing, hessian_at=hat, distance=NaN,
            steps=0, gain=Inf, full_hessians=full_hessians,
            subset_hessians=subset_hessians, fcalls=0, gcalls=0, escapes=0,
            saddle=false, ladder_tried=false, ladder_gain=0.0,
            probe=_ctsem_no_probe(), history=history)
    end
    at_x = exact                        # `H` is the exact Hessian at `x`
    while true                          # once, and again after each escape
        converged = false
        mu = 0.0; prevgain = Inf
        while take_steps && steps < maxit
            nt = newton(H, G, mu)
            nt === nothing && break
            gain = nt.gain
            if gain < tol
                converged = true
                break
            end
            trial = search(nt.step, dot(G, nt.step))
            if !trial.ok
                if !at_x
                    H = hess(x); hat = copy(x); exact = true; at_x = true
                    H === nothing && break
                else
                    mu = mu == 0 ? 1e-4 : 10mu
                    mu > 1e2 && break
                end
                continue
            end
            Gn = similar(G)
            fg!(nothing, Gn, trial.x); gcalls += 1
            x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
            remember!("newton", gain, trial.alpha)
            mu = trial.alpha == 1 ? mu / 10 : mu
            mu < 1e-8 && (mu = 0.0)
            report(gain)
            # Only the exact variant refreshes on slow contraction; the chord
            # and the subset keep their matrix, which is the point of them. A
            # step that fails outright still gets the exact Hessian (above).
            if curvature === :exact && gain / prevgain > contraction && steps > 1
                H = hess(x); hat = copy(x); exact = true; at_x = true
                H === nothing && break
                prevgain = Inf
            else
                prevgain = gain
            end
        end
        H === nothing && break
        # The final Hessian. A chord Hessian the steps converged on, taken
        # within `reuse_se` standard errors of here, is kept: the exact one at
        # this point would differ from it by less than the certification can
        # see, and on the Laplace route it would cost another `2 npar`
        # gradients. Anything else is replaced by the exact Hessian here, and a
        # kept or subset Hessian that judged the gain against itself is
        # followed by steps on the exact one until they agree -- otherwise the
        # certification finds the point unfinished and resumes the optimiser
        # (measured before the finish carried on: 1045 iterations to the cap
        # on a 30-subject model the exact finish closes in 10 steps).
        if !at_x
            keep = curvature === :chord && exact && converged &&
                _ctsem_hessian_distance(_ctsem_information_split(H;
                    rtol=flat_rtol, negative=negative), hat, x) <= reuse_se
            if !keep
                for _ in 1:5
                    if !at_x
                        H = hess(x); hat = copy(x); exact = true; at_x = true
                    end
                    H === nothing && break
                    take_steps || break
                    nt = newton(H, G, 0.0)
                    gain = nt === nothing ? Inf : nt.gain
                    (nt === nothing || gain < tol || steps >= maxit + 5) && break
                    trial = search(nt.step, dot(G, nt.step))
                    trial.ok || break
                    Gn = similar(G)
                    fg!(nothing, Gn, trial.x); gcalls += 1
                    x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
                    remember!("exact", gain, trial.alpha)
                    report(gain)
                end
                # A step on the last round leaves the Hessian one point behind.
                if !at_x && H !== nothing
                    H = hess(x); hat = copy(x); exact = true; at_x = true
                end
            end
        end
        H === nothing && break
        # A saddle: the ascent is along the negative curvature, which the Newton
        # step, living in the trusted directions, cannot see. Measured on the
        # rank-deficient mvmix laplace fixture before this existed: the Newton
        # step gained its predicted 2.2e-7 and L-BFGS resumed from there crept
        # 0.02 nats in 1000 iterations along the same direction, while a
        # maximum 2 nats higher sat in the other basin.
        take_steps || break
        split = _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
        (split === nothing || !any(split.negative)) && break
        saddle = true
        escapes >= max_escapes && break
        ladder_tried = true
        tried = _ctsem_saddle_ladder(valueof, x, -f, split, -G; ladder=ladder)
        fcalls += tried.evaluations
        tried.best === nothing && break
        Gn = similar(G)
        fn = fg!(0.0, Gn, tried.best.point); fcalls += 1; gcalls += 1
        # The point the ladder measured, through the optimiser's own trial path:
        # a point the route refuses a gradient at is not one to continue from.
        (isfinite(fn) && fn < f) || break
        ladder_gain += f - fn
        x = collect(Float64, tried.best.point); f = fn; G = Gn
        steps += 1; escapes += 1
        remember!("saddle", NaN, tried.best.length)
        report(NaN)
        H = hess(x); hat = copy(x); exact = true; at_x = true
        saddle = false; ladder_tried = false
        H === nothing && break
    end
    split = H === nothing ? nothing :
        _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
    if H !== nothing
        nt = newton(H, G, 0.0)
        gain = nt === nothing ? Inf : nt.gain
    end
    distance = H === nothing ? NaN : _ctsem_hessian_distance(split, hat, x)
    # What the directions the curvature does not trust are worth, at the point
    # being returned, against its own value: the probe the certification
    # reads (`.ctBackendCertify()`).
    probed = _ctsem_no_probe()
    if probe && split !== nothing && all(isfinite, G)
        residual = _ctsem_flat_residual(split, -G)
        if norm(residual) > 0
            probed = _ctsem_flat_probe(valueof, x, -f, residual;
                lengths=probe_lengths)
            fcalls += probed.evaluations
        end
    end
    (x=x, f=f, G=G, hessian=H === nothing ? nothing : -H, hessian_at=hat,
     distance=distance, steps=steps, gain=gain, full_hessians=full_hessians,
     subset_hessians=subset_hessians, fcalls=fcalls, gcalls=gcalls,
     escapes=escapes, saddle=saddle, ladder_tried=ladder_tried,
     ladder_gain=ladder_gain, probe=probed, history=history)
end

"""
    _ctsem_trial_closure(objective, gradient_method)

The optimiser's trial path (`fg!`) over `objective`, for a caller that is not
`ctsem_optimize`: the route decides what a usable point is
(`_ctsem_optimise_trial`), and an unusable one comes back as the sentinel value
with a zero gradient.
"""
function _ctsem_trial_closure(objective, gradient_method)
    invalid = floatmax(Float64) / 1e8
    limit = sqrt(floatmax(Float64))
    log = _ctsem_optimise_log(objective, false)
    return function (F, G, x)
        trial = _ctsem_optimise_trial(objective, x, G !== nothing,
            gradient_method, limit, log)
        if !trial.valid
            G !== nothing && fill!(G, 0.0)
            return F === nothing ? nothing : invalid
        end
        G !== nothing && (G .= -trial.evaluated.gradient)
        return F === nothing ? nothing : -trial.evaluated.value
    end
end

"""
    ctsem_endgame(objective, values; gradient_method, flat_rtol, probe_lengths)

The certification's numbers at a point the optimiser left without its finish --
an iteration cap, a stall, a finish that could not form a Hessian -- so R can
certify it the way it certifies a finished one: the exact Hessian at `values`
and what the flat directions are worth there. No step is taken. Whether to go
on from here is R's decision, and it goes on by resuming the optimiser, whose
own finish takes the steps.

The Hessian is of the MAXIMISED objective, and a 1x1 zero when it could not be
formed: an empty matrix would hang the bridge, as would an empty probe
direction, which is `[0.0]` when `probe_ran` is false. The counts are named as
`ctsem_optimize` names its own, so whatever tallies engine runs counts this
one too.
"""
function ctsem_endgame(objective::CTSEMOptimisable, values::AbstractVector;
        gradient_method=:adjoint, flat_rtol::Real=_CTSEM_FLAT_RTOL,
        probe_lengths=_CTSEM_FLAT_PROBE_LENGTHS)
    x = collect(Float64, values)
    fg! = _ctsem_trial_closure(objective, gradient_method)
    G = zeros(length(x))
    f = fg!(0.0, G, x)
    out = _ctsem_newton_finish(objective, x, f, G, fg!; take_steps=false,
        probe=true, curvature=:exact, flat_rtol=flat_rtol,
        probe_lengths=collect(Float64, probe_lengths))
    probed = out.probe
    return (minimizer=x, maximum_loglik=-f, gradient=-G,
        hessian=out.hessian === nothing ? zeros(1, 1) : out.hessian,
        iterations=0, f_calls=1 + out.fcalls, g_calls=1 + out.gcalls,
        newton_steps=0, newton_hessians=out.full_hessians,
        probe_ran=!isempty(probed.direction), probe_gain=probed.gain,
        probe_length=probed.length, probe_longest=probed.longest,
        probe_direction=isempty(probed.direction) ? [0.0] : probed.direction)
end

export ctsem_endgame

"""
    ctsem_flat_probe(objective, values, direction; lengths)

The flat-direction probe (`_ctsem_flat_probe`) at `values` along `direction`,
with the value at `values` taken by the same predicate as the probe points --
for a certification R assembles from a Hessian it already holds
(`.ctBackendCertification()`). On the Laplace route that predicate refuses a
point whose inner solve did not converge, which the R-side probe this replaced
could not see. The fields are named as `ctsem_optimize` names its probe's, so R
reads both one way; `probe_ran` false and a `[0.0]` direction when there was
nothing to probe.
"""
function ctsem_flat_probe(objective, values::AbstractVector,
        direction::AbstractVector; lengths=_CTSEM_FLAT_PROBE_LENGTHS)
    x = collect(Float64, values)
    value_at = y -> _ctsem_probe_value(objective, y)
    base = value_at(x)
    isfinite(base) || return (probe_ran=false, probe_gain=0.0,
        probe_length=0.0, probe_longest=0.0, probe_direction=[0.0])
    out = _ctsem_flat_probe(value_at, x, base, collect(Float64, direction);
        lengths=collect(Float64, lengths))
    return (probe_ran=!isempty(out.direction), probe_gain=out.gain,
        probe_length=out.length, probe_longest=out.longest,
        probe_direction=isempty(out.direction) ? [0.0] : out.direction)
end

export ctsem_flat_probe
