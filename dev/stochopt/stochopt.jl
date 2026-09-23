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

A reused curvature is refreshed exactly when the gain stops contracting by
`contraction` per step, and a failed step raises a Levenberg damping rather
than giving up. Always finishes with an exact full Hessian at the final point
-- the one a fit computes anyway to certify and for standard errors -- and
keeps iterating with exact curvature while that Hessian predicts more than
`tol`.
"""
function endgame(o, x0::AbstractVector; curvature::Symbol=:chord, tol=1e-6,
    maxit::Integer=100, submax::Integer=0, rng=Random.default_rng(),
    contraction::Real=0.25, finalit::Integer=20)
    N = nsubjects(o)
    led = Ledger(N)
    x = collect(Float64, x0)
    f, g = evalgrad(led, o, x)
    isfinite(f) || return (status="nonfinite start", x=x, f=f, iterations=0,
        refreshes=0, final_gain=NaN, extra_step=0, last_gain=NaN,
        grad=led.grad, value=led.value, hess=led.hess, scores=led.scores,
        seconds=led.seconds)

    function curvature_at(kind, y)
        if kind in (:exact, :chord)
            return evalhess(led, o, y)
        elseif kind === :subset
            NU = nunits(o)
            m = min(NU, submax)
            m >= NU && return evalhess(led, o, y)
            idx = sort(randperm(rng, NU)[1:m])
            so = subset_objective(o, idx)
            so === nothing && return evalhess(led, o, y)
            Hs = evalhess(led, so, y)
            Hs === nothing && return nothing
            # Priors off (every ML fit): the whole Hessian is likelihood and
            # scales by the share of subjects evaluated.
            return Hs .* (N / nsubjects(so))
        elseif kind === :bhhh
            S = evalscores(led, o, y)
            S === nothing && return nothing
            return -(S' * S)
        end
        error("unknown curvature $(kind)")
    end

    # Damped Newton on the maximisation. `mu` is a Levenberg damping in units
    # of the largest curvature: zero while full steps succeed, raised when a
    # step fails, so a saddle or a boundary gets a shorter, gradient-leaning
    # step instead of a failed one. Negative and flat directions are floored
    # at `rtol` of the largest curvature rather than dropped, so the gradient
    # there is still followed, gently.
    function damped_step(H, g, mu; rtol=1e-8)
        E = eigen(Symmetric(-Matrix(H)))
        lmax = maximum(abs, E.values)
        lam = max.(E.values, rtol * lmax) .+ mu * lmax
        c = E.vectors' * g
        step = E.vectors * (c ./ lam)
        # the undamped predicted gain, which is what convergence is judged on
        keep = E.values .> rtol * lmax
        gain_undamped = 0.5 * sum(abs2.(c[keep]) ./ E.values[keep])
        (step=step, gain=0.5 * dot(g, step), gain_undamped=gain_undamped)
    end

    # One Newton loop with a given curvature source. Returns when the undamped
    # gain under the current curvature is below `tol`.
    function run!(kind, H, maxsteps)
        mu = 0.0
        prevgain = Inf
        steps = 0; refreshes = 0; st = "maxit"
        while steps < maxsteps
            s = damped_step(H, g, mu)
            if s.gain_undamped < tol
                st = "converged"; break
            end
            ls = linesearch(led, o, x, f, s.step, s.gain)
            if ls === nothing
                if mu == 0 && kind !== :exact
                    # a stale or approximate curvature: refresh it exactly
                    H = curvature_at(:exact, x); refreshes += 1
                    kind = :exact
                    H === nothing && (st = "no curvature"; break)
                else
                    mu = mu == 0 ? 1e-4 : 10mu
                    mu > 1e2 && (st = "linesearch failed"; break)
                end
                continue
            end
            steps += 1
            x .= ls.x
            f, g = evalgrad(led, o, x)
            isfinite(f) || (st = "nonfinite"; break)
            mu = ls.alpha == 1 ? mu / 10 : mu
            mu < 1e-8 && (mu = 0.0)
            # A reused curvature converges linearly at the rate it is wrong
            # by; when the gain is not shrinking fast enough, refresh.
            ratio = s.gain_undamped / prevgain
            prevgain = s.gain_undamped
            if kind === :bhhh
                H = curvature_at(:bhhh, x); refreshes += 1
            elseif ratio > contraction && steps > 1
                H = curvature_at(:exact, x); refreshes += 1
                kind = :exact
                prevgain = Inf
                H === nothing && (st = "no curvature"; break)
            end
        end
        (H=H, steps=steps, refreshes=refreshes, status=st)
    end

    H = curvature_at(curvature, x)
    H === nothing && return (status="no curvature", x=x, f=f, iterations=0,
        refreshes=0, final_gain=NaN, extra_step=0, last_gain=NaN,
        grad=led.grad, value=led.value, hess=led.hess, scores=led.scores,
        seconds=led.seconds)
    r1 = run!(curvature, H, maxit)

    # The certification Hessian every fit computes. If it disagrees -- the
    # fit is not finished under the exact curvature -- keep going with it,
    # exactly, and certify again at the new point.
    Hc = evalhess(led, o, x)
    extra = 0
    final_gain = NaN
    status = r1.status
    for round in 1:3
        Hc === nothing && break
        final_gain = damped_step(Hc, g, 0.0).gain_undamped
        final_gain < tol && break
        r2 = run!(:exact, Hc, finalit)
        extra += r2.steps
        status = r1.status * "+" * r2.status
        r2.steps == 0 && break
        Hc = evalhess(led, o, x)
    end
    (status=status, x=x, f=f, iterations=r1.steps, refreshes=r1.refreshes,
     final_gain=final_gain, extra_step=extra, last_gain=final_gain,
     grad=led.grad, value=led.value, hess=led.hess, scores=led.scores,
     seconds=led.seconds)
end

# ---------------------------------------------------------------- phase 2
# Progressive-batch L-BFGS. Units are permuted once and the batch is always a
# prefix of that permutation, so growing it adds units rather than resampling
# -- within a stage the objective is deterministic and the line search is
# ordinary. The batch objective is scaled to full-data units, (N/n) * loglik_B,
# so curvature pairs gathered on a small batch estimate the full curvature and
# the memory is carried across growth rather than thrown away.
#
# Growth is decided in the optimiser's own metric, in objective units. With H
# the current inverse-Hessian approximation and G the scaled batch gradient,
# the step predicts a gain of G'HG/2; the sampling noise of that prediction is
# tr(H Cov(G))/2 with Cov(G) = N^2 (1 - n/N) S/n from the per-unit scores. When
# noise reaches `theta` times the predicted gain the batch is too small to say
# which way is up, and it grows by the factor that would bring the ratio back
# to `theta` (Byrd, Chin, Nocedal & Wu 2012's norm test, in this metric).

struct LBFGSMemory
    S::Vector{Vector{Float64}}
    Y::Vector{Vector{Float64}}
    m::Int
end
LBFGSMemory(m::Integer) = LBFGSMemory(Vector{Float64}[], Vector{Float64}[], Int(m))

function pushpair!(M::LBFGSMemory, s, y)
    sy = dot(s, y)
    sy > 1e-12 * norm(s) * norm(y) || return false   # keep H positive definite
    push!(M.S, copy(s)); push!(M.Y, copy(y))
    length(M.S) > M.m && (popfirst!(M.S); popfirst!(M.Y))
    true
end

# H*q for the inverse Hessian of the NEGATIVE objective (so a positive
# definite metric); `h0` is the initial scaling when memory is empty.
function hmul(M::LBFGSMemory, q::AbstractVector, h0::Float64)
    k = length(M.S)
    k == 0 && return h0 .* q
    q = copy(q)
    a = zeros(k); rho = [1 / dot(M.Y[i], M.S[i]) for i in 1:k]
    for i in k:-1:1
        a[i] = rho[i] * dot(M.S[i], q); q .-= a[i] .* M.Y[i]
    end
    gamma = dot(M.S[k], M.Y[k]) / dot(M.Y[k], M.Y[k])
    r = gamma .* q
    for i in 1:k
        b = rho[i] * dot(M.Y[i], r); r .+= (a[i] - b) .* M.S[i]
    end
    r
end

# Scaled batch gradient and per-unit score rows, from one sweep. E[G'HG] is
# G'HG + tr(H Cov G), which is why the growth test compares the two.
function scaled_scores(led::Ledger, so, x, N)
    S = evalscores(led, so, x)
    S === nothing && return nothing
    n = size(S, 1)
    nsub = nsubjects(so)
    # rows carry the batch's prior share; with priors off that is zero
    G = vec(sum(S; dims=1)) .* (N / nsub)
    (G=G, S=S, n=n, nsub=nsub)
end

"""
    pbatch(o, x0; n0, theta, memory, tol_switch, maxit, grow)

Progressive-batch L-BFGS on the maximisation of `o`. Returns when the batch is
the full data and the predicted gain G'HG/2 is below `tol_switch`, so the
endgame can take over. `grow=false` is the control: the same optimiser on the
full data throughout.
"""
function pbatch(o, x0::AbstractVector; n0::Integer=0, theta::Real=0.5,
    memory::Integer=10, tol_switch::Real=1e-1, maxit::Integer=2000,
    grow::Bool=true, rng=Random.default_rng(), c1=1e-4)
    N = nsubjects(o)
    NU = nunits(o)
    led = Ledger(N)
    perm = randperm(rng, NU)
    m = grow ? clamp(n0 > 0 ? Int(n0) : max(20, cld(NU, 32)), 1, NU) : NU
    # Too few units to batch: this is plain L-BFGS, by construction.
    m >= NU ÷ 4 && (m = NU)
    batch() = m == NU ? o : subset_objective(o, sort(perm[1:m]))
    so = batch()
    x = collect(Float64, x0)
    ev = scaled_scores(led, so, x, N)
    ev === nothing && return (status="no scores", x=x)
    f = evalvalue_scaled(led, so, x, N)
    G = ev.G
    M = LBFGSMemory(memory)
    h0 = 0.1 / max(norm(G), 1e-12)          # first step of length 0.1
    history = Tuple{Int,Int}[]               # (iteration, batch size) at growth
    status = "maxit"
    it = 0
    for k in 1:maxit
        it = k
        d = hmul(M, G, h0)                   # ascent direction
        gain = 0.5 * dot(G, d)
        if m < NU
            # noise of the predicted gain from the per-unit spread
            Sbar = vec(sum(ev.S; dims=1)) ./ ev.n
            acc = 0.0
            for i in 1:ev.n
                di = ev.S[i, :] .- Sbar
                acc += dot(di, hmul(M, di, h0))
            end
            V = acc / max(ev.n - 1, 1)                   # per-unit, unscaled
            scale = (N / ev.nsub)^2 * ev.n * (1 - m / NU)  # Cov(G) = scale * V-ish
            noise = 0.5 * scale * V
            if noise > theta * gain
                newm = min(NU, max(2m, ceil(Int, m * noise / (theta * gain))))
                push!(history, (k, newm))
                m = newm
                so = batch()
                ev = scaled_scores(led, so, x, N)
                ev === nothing && (status = "no scores"; break)
                f = evalvalue_scaled(led, so, x, N)
                G = ev.G
                continue                     # no pair across a change of objective
            end
        elseif gain < tol_switch
            status = "switch"
            break
        end
        # backtracking on the (scaled) batch value
        alpha = 1.0
        xn = x; fn = -Inf; ok = false
        for _ in 1:40
            xn = x .+ alpha .* d
            fn = evalvalue_scaled(led, so, xn, N)
            if isfinite(fn) && fn >= f + c1 * alpha * 2 * gain
                ok = true; break
            end
            alpha /= 2
        end
        if !ok
            # a stale memory is the usual cause; drop it once before giving up
            if !isempty(M.S)
                empty!(M.S); empty!(M.Y); h0 = 0.1 / max(norm(G), 1e-12)
                continue
            end
            status = m < NU ? "linesearch (grow)" : "linesearch"
            if m < NU
                m = NU; so = o; push!(history, (k, m))
                ev = scaled_scores(led, so, x, N); f = evalvalue_scaled(led, so, x, N)
                G = ev.G; continue
            end
            break
        end
        evn = scaled_scores(led, so, xn, N)
        evn === nothing && (status = "no scores"; break)
        # maximising: the pair for the negative objective is (s, -(Gn - G))
        pushpair!(M, xn .- x, -(evn.G .- G))
        x = xn; f = fn; ev = evn; G = evn.G
    end
    (status=status, x=x, f=f, iterations=it, final_batch=m,
     # never empty: a zero-length vector deadlocks the JuliaConnectoR bridge
     growth=vcat(-1, [h[1] for h in history]), sizes=vcat(-1, [h[2] for h in history]),
     grad=led.grad, value=led.value, hess=led.hess, scores=led.scores,
     seconds=led.seconds)
end

function evalvalue_scaled(led, so, x, N)
    v = evalvalue(led, so, x)
    v * (N / nsubjects(so))
end

# ---------------------------------------------------------------- timings
"""Minimum-of-`reps` seconds for one full value, gradient, Hessian, scores."""
function primitive_times(o, x, reps::Integer=3)
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
function subset_hessian_time(o, x, m::Integer, reps::Integer=2)
    so = subset_objective(o, collect(1:min(m, nsubjects(o))))
    so === nothing && return NaN
    x = collect(Float64, x)
    minimum(begin
        t = @elapsed try CT.ctsem_hessian(so, x) catch end
        t
    end for _ in 1:reps)
end

end # module
