# Prototype optimiser endgames, kept outside the engine so that editing them
# costs a re-include rather than a precompile. Loaded into Main of a session
# whose ContinuousTimeSEM is already loaded (ctsem does that on first use):
#
#   include("dev/stochopt/stochopt.jl")      # re-include freely
#
# Everything maximises the engine's objective (a log likelihood or log
# posterior); `H` below is always the Hessian of that objective, so negative
# definite at a maximum.
module StochOpt

using LinearAlgebra, Printf, Random
import ContinuousTimeSEM as CT

# ---------------------------------------------------------------- cost ledger
# Every engine call goes through here, so a run's cost is counted rather than
# inferred from a wall clock on a shared machine. Units are full-data
# gradients' worth of subjects: a gradient over m of N subjects adds m/N to
# `grad`, and the conversion to wall time uses per-primitive timings measured
# once per model (see `primitive_times`).
mutable struct Ledger
    nsubjects::Int
    grad::Float64
    value::Float64
    hess::Float64
    scores::Float64
    seconds::Float64
end
Ledger(n::Integer) = Ledger(Int(n), 0.0, 0.0, 0.0, 0.0, 0.0)

nsubjects(o::CT.CTSEMObjective) = length(o.subject_objectives)
# Subjects actually evaluated: a unit-subset keeps the full list underneath.
nsubjects(o::CT.CTSEMLaplaceObjective) = sum(length, o.units.members; init=0)

function evalvalue(led::Ledger, o, x)
    t = @elapsed r = try
        CT.ctsem_evaluate(o, x; gradient=false)
    catch
        nothing
    end
    led.seconds += t
    led.value += nsubjects(o) / led.nsubjects
    r === nothing && return -Inf
    v = Float64(r.value)
    isfinite(v) ? v : -Inf
end

function evalgrad(led::Ledger, o, x)
    t = @elapsed r = try
        CT.ctsem_evaluate(o, x; gradient=true)
    catch
        nothing
    end
    led.seconds += t
    led.grad += nsubjects(o) / led.nsubjects
    r === nothing && return (-Inf, fill(NaN, length(x)))
    v = Float64(r.value)
    g = collect(Float64, r.gradient)
    (isfinite(v) && all(isfinite, g)) ? (v, g) : (-Inf, g)
end

function evalhess(led::Ledger, o, x)
    t = @elapsed H = try
        Matrix{Float64}(CT.ctsem_hessian(o, collect(Float64, x)))
    catch
        nothing
    end
    led.seconds += t
    led.hess += nsubjects(o) / led.nsubjects
    H
end

function evalscores(led::Ledger, o, x)
    t = @elapsed r = try
        CT.ctsem_subject_gradients(o, collect(Float64, x))
    catch
        nothing
    end
    led.seconds += t
    led.scores += nsubjects(o) / led.nsubjects
    r === nothing ? nothing : Matrix{Float64}(r.scores)
end

# ---------------------------------------------------------------- subsets
# A CTSEMObjective over some of its subjects. The subject objectives are
# self-contained; the prior is global and so is kept whole -- the caller
# scales the *likelihood* part, not the prior, when it rescales curvature.
subset_objective(o::CT.CTSEMObjective, idx::AbstractVector{<:Integer}) =
    CT.CTSEMObjective(o.params, o.subject_objectives[idx], nothing,
        o.prior_index, o.prior_scale, o.prior_weight,
        o.ti_missing_parameter, o.ti_missing_mu, o.ti_missing_sigma)
# Laplace: a subset of whole UNITS (outer-level groups), sharing everything
# else. Every evaluation path goes units.members[U] -> subject_objectives[i],
# so restricting `units` restricts the objective; the full subject list stays
# in the wrapped objective and is never touched for excluded units. Pairing the
# full spec with a *subset objective* instead would silently misassign groups.
function subset_objective(L::CT.CTSEMLaplaceObjective, keep::AbstractVector{<:Integer})
    u = L.units
    units = CT.CTSEMLaplaceUnits(u.members[keep], u.offsets[keep], u.dims[keep],
        u.blocks[keep])
    n = length(keep)
    CT.CTSEMLaplaceObjective{typeof(L.objective)}(L.objective, L.spec, units,
        [zeros(units.dims[U]) for U in 1:n], L.inner_maxiter, L.inner_tol,
        [Dict{Any,Any}() for _ in eachindex(L.workspaces)], zeros(Int, n),
        zeros(n), falses(n), falses(n), falses(n), falses(n))
end
subset_objective(o, idx) = nothing   # joint objective: not subsettable here

# Units, not subjects, are what a Laplace subset draws; for a single-level
# model they are the same thing.
nunits(o::CT.CTSEMObjective) = nsubjects(o)
nunits(o::CT.CTSEMLaplaceObjective) = length(o.units.members)

# ---------------------------------------------------------------- Newton step
# The step over the trusted subspace of -H: directions whose curvature is below
# `rtol` of the largest are left out, as the certification does, so a flat
# ridge does not produce an enormous step. Returns the step, the predicted gain
# 1/2 g' step, and how many directions were excluded or negative.
function newton_step(H::AbstractMatrix, g::AbstractVector; rtol=1e-8)
    E = eigen(Symmetric(-Matrix(H)))
    lmax = maximum(abs, E.values)
    keep = E.values .> rtol * lmax
    V = E.vectors[:, keep]
    c = V' * g
    step = V * (c ./ E.values[keep])
    (step=step, gain=0.5 * dot(g, step), nneg=count(<(0), E.values),
     nflat=count(!, keep))
end

# Armijo backtracking on value only (a value costs no reverse pass).
function linesearch(led, o, x, f, step, gain; c1=1e-4, maxhalve=30)
    alpha = 1.0
    for _ in 1:maxhalve
        xn = x .+ alpha .* step
        fn = evalvalue(led, o, xn)
        if isfinite(fn) && fn >= f + c1 * alpha * 2 * gain
            return (x=xn, f=fn, alpha=alpha)
        end
        alpha /= 2
        alpha * 2 * gain < eps(abs(f)) && break
    end
    nothing
end

# ---------------------------------------------------------------- endgames
"""
    endgame(o, x0; curvature, tol, maxit, submax, rng)

From `x0` (typically where a loosely-stopped L-BFGS left off), iterate Newton
steps on the FULL objective with curvature from `curvature`:

- `:exact`    exact Hessian, refreshed whenever a step is not taken whole
- `:chord`    one exact Hessian at x0, reused; refreshed only on failure
- `:subset`   Hessian of a random `submax`-subject subset, likelihood scaled
               by N/m; refreshed (new subset) only on failure
- `:bhhh`     -(scores' scores), refreshed every iteration (one sweep each)

Always finishes with one exact full Hessian at the final point -- the one a fit
computes anyway to certify and for standard errors -- and one more Newton step
from it if that Hessian still predicts more than `tol`.
"""
function endgame(o, x0::AbstractVector; curvature::Symbol=:chord, tol=1e-6,
    maxit::Integer=100, submax::Integer=0, rng=Random.default_rng())
    N = nsubjects(o)
    led = Ledger(N)
    x = collect(Float64, x0)
    f, g = evalgrad(led, o, x)
    isfinite(f) || return (status="nonfinite start", led=led, x=x, f=f)

    function curvature_at(y, gy)
        if curvature in (:exact, :chord)
            return evalhess(led, o, y)
        elseif curvature === :subset
            NU = nunits(o)
            m = min(NU, submax)
            if m >= NU
                return evalhess(led, o, y)
            end
            idx = sort(randperm(rng, NU)[1:m])
            so = subset_objective(o, idx)
            so === nothing && return evalhess(led, o, y)
            Hs = evalhess(led, so, y)
            Hs === nothing && return nothing
            # Prior curvature is not a sum over subjects; with priors off (as in
            # every ML fit) the whole Hessian is likelihood and scales by the
            # share of subjects evaluated.
            return Hs .* (N / nsubjects(so))
        elseif curvature === :bhhh
            S = evalscores(led, o, y)
            S === nothing && return nothing
            return -(S' * S)
        end
        error("unknown curvature $(curvature)")
    end

    H = curvature_at(x, g)
    H === nothing && return (status="no curvature", led=led, x=x, f=f)
    iters = 0
    refreshes = 1
    status = "maxit"
    lastgain = Inf
    for k in 1:maxit
        iters = k
        s = newton_step(H, g)
        lastgain = s.gain
        if s.gain < tol
            status = "converged"
            break
        end
        ls = linesearch(led, o, x, f, s.step, s.gain)
        if ls === nothing
            # The quadratic model is wrong here. A fresh exact Hessian is the
            # remedy for a stale or approximate one; for :exact it is not.
            if curvature === :exact
                status = "linesearch failed"
                break
            end
            H = curvature === :bhhh ? evalhess(led, o, x) : curvature_at(x, g)
            refreshes += 1
            H === nothing && (status = "no curvature"; break)
            continue
        end
        x = ls.x
        f, g = evalgrad(led, o, x)
        isfinite(f) || (status = "nonfinite"; break)
        if curvature === :bhhh || (curvature === :exact && ls.alpha < 1)
            H = curvature_at(x, g)
            refreshes += 1
            H === nothing && (status = "no curvature"; break)
        end
    end

    # The certification Hessian, which every arrangement pays.
    Hc = evalhess(led, o, x)
    final_gain = NaN
    extra = 0
    if Hc !== nothing
        s = newton_step(Hc, g)
        final_gain = s.gain
        if s.gain > tol
            ls = linesearch(led, o, x, f, s.step, s.gain)
            if ls !== nothing
                x = ls.x
                f, g = evalgrad(led, o, x)
                extra = 1
                final_gain = newton_step(Hc, g).gain
            end
        end
    end
    (status=status, x=x, f=f, iterations=iters, refreshes=refreshes,
     final_gain=final_gain, extra_step=extra, last_gain=lastgain,
     grad=led.grad, value=led.value, hess=led.hess, scores=led.scores,
     seconds=led.seconds)
end

# ---------------------------------------------------------------- timings
"""Minimum-of-`reps` seconds for one full value, gradient, Hessian, scores."""
function primitive_times(o, x; reps::Integer=3)
    x = collect(Float64, x)
    tmin(f) = minimum(begin
        t = @elapsed try f() catch end
        t
    end for _ in 1:reps)
    tv = tmin(() -> CT.ctsem_evaluate(o, x; gradient=false))
    tg = tmin(() -> CT.ctsem_evaluate(o, x; gradient=true))
    th = tmin(() -> CT.ctsem_hessian(o, x))
    ts = tmin(() -> CT.ctsem_subject_gradients(o, x))
    (value=tv, grad=tg, hess=th, scores=ts)
end

# Subset-Hessian timing, for the cost model of :subset at size m.
function subset_hessian_time(o, x, m::Integer; reps::Integer=2)
    so = subset_objective(o, collect(1:min(m, nsubjects(o))))
    so === nothing && return NaN
    x = collect(Float64, x)
    minimum(begin
        t = @elapsed try CT.ctsem_hessian(so, x) catch end
        t
    end for _ in 1:reps)
end

end # module
