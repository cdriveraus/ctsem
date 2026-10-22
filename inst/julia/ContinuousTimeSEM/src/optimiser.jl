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
exactly (a Newton step does; an L-BFGS step leaves it to `CTSEMDirectional`),
and the point itself where the stepper holds it apart from the objective's
evaluations: the finish's evaluations include Hessian columns and rejected
trials, so the last point evaluated is not its iterate. `nothing` for L-BFGS,
whose caller reads the last accepted trial instead.
"""
struct CTSEMIterate
    iteration::Int
    value::Float64
    g_norm::Float64
    gain::Float64
    x::Union{Nothing,Vector{Float64}}
end
CTSEMIterate(iteration, value, g_norm, gain=NaN) =
    CTSEMIterate(iteration, value, g_norm, gain, nothing)

"""
A jump the line search saw (`_ctsem_lbfgs`'s `jump`): from `x` along `s`, a
trial at `alpha_lo` kept at least a quarter of the excess over the line
search's linear model, `excess`, that a trial at `alpha_hi`, a hundred times
longer or more, had. A smooth objective's excess shrinks with the square of
the step, so that is a jump of about `excess` between the two, not curvature:
on the Laplace route, a unit whose term moved with where its inner solve
stopped (`_ctsem_jump_units`). `seen` is false until one is.
"""
mutable struct CTSEMJump
    seen::Bool
    x::Vector{Float64}
    s::Vector{Float64}
    alpha_hi::Float64
    alpha_lo::Float64
    excess::Float64
end
CTSEMJump() = CTSEMJump(false, Float64[], Float64[], NaN, NaN, NaN)

function _ctsem_jump_record!(jump::CTSEMJump, x, s, alpha_hi, alpha_lo, excess)
    jump.seen = true
    jump.x = copy(x)
    jump.s = copy(s)
    jump.alpha_hi = alpha_hi
    jump.alpha_lo = alpha_lo
    jump.excess = excess
    return jump
end

"""Whether the route can name the units a jump came from
(`_ctsem_jump_units`), and so whether its line searches watch for one: the
Laplace objective and the hybrid built on it. Elsewhere a line search behaves
exactly as it did before the watch."""
_ctsem_jump_capable(::Any) = false

"""Whether the objective screens and flags its own units as the climb goes
(`_ctsem_anytime!`, the Laplace continuation with `anytime` on): the loop
then gives it each accepted step's predicted gain and any jump the watch
names, and re-reads the objective where it changed."""
_ctsem_anytime_capable(::Any) = false

# At or above this a trial's value is the sentinel a route returns for a point
# it refuses (`ctsem_optimize`'s `invalid_objective`, `floatmax / 1e8`), not a
# measurement of the objective.
const _CTSEM_REFUSED = floatmax(Float64) / 1e9

"""Curvature pairs for the two-loop recursion, on the MINIMISED objective, and
the learned diagonal of the initial matrix (`_ctsem_lbfgs_diagonal!`'s `B`,
`nothing` until a run learns one)."""
mutable struct CTSEMLBFGSMemory
    S::Vector{Vector{Float64}}
    Y::Vector{Vector{Float64}}
    rho::Vector{Float64}
    m::Int
    diagonal::Union{Nothing,Vector{Float64}}
end
CTSEMLBFGSMemory(m::Integer) = CTSEMLBFGSMemory(Vector{Float64}[],
    Vector{Float64}[], Float64[], Int(m), nothing)

"""What `_ctsem_lbfgs` returns: the point, the reasons it stopped, and the
curvature it ended with, for a later run to start from (its `pairs`)."""
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
    memory::CTSEMLBFGSMemory
end

"""
    _ctsem_lbfgs_carry(pairs, m, n)

The memory a run of `n` parameters keeping `m` pairs starts from: a fresh one
for `pairs === nothing`, else a copy of `pairs`' newest `m`, so the run never
writes into the memory it was handed. A memory of another length is refused.
The pairs' vectors are shared, not copied: nothing writes into one once it is
stored.
"""
function _ctsem_lbfgs_carry(pairs::Union{Nothing,CTSEMLBFGSMemory}, m::Integer,
        n::Integer)
    pairs === nothing && return CTSEMLBFGSMemory(m)
    _ctsem_lbfgs_check(pairs, n)
    keep = max(1, length(pairs.S) - Int(m) + 1):length(pairs.S)
    CTSEMLBFGSMemory(pairs.S[keep], pairs.Y[keep], pairs.rho[keep], Int(m),
        pairs.diagonal === nothing ? nothing : copy(pairs.diagonal))
end

"""Refuses a memory whose vectors are not of length `n`."""
function _ctsem_lbfgs_check(pairs::CTSEMLBFGSMemory, n::Integer)
    sized = all(v -> length(v) == n, pairs.S) && all(v -> length(v) == n, pairs.Y) &&
        (pairs.diagonal === nothing || length(pairs.diagonal) == n)
    sized || throw(DimensionMismatch(string("an L-BFGS memory made for another ",
        "number of parameters cannot start a run of ", n)))
    return nothing
end

# The curvature an objective keeps between `ctsem_optimize` runs on it
# (`lbfgs_pairs = :objective`), or `nothing` where it keeps none and each run
# starts cold; `_ctsem_lbfgs_hold!` leaves a run's there.
_ctsem_lbfgs_held(::Any) = nothing
_ctsem_lbfgs_hold!(::Any, memory) = nothing

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

"""Drop the `k` newest curvature pairs."""
function _ctsem_lbfgs_drop!(M::CTSEMLBFGSMemory, k::Integer)
    for _ in 1:min(k, length(M.S))
        pop!(M.S); pop!(M.Y); pop!(M.rho)
    end
    return M
end

_ctsem_lbfgs_reset!(M::CTSEMLBFGSMemory) =
    (empty!(M.S); empty!(M.Y); empty!(M.rho); M.diagonal = nothing; M)

# H * q for the inverse-Hessian approximation, with `dinv` the metric's inverse
# diagonal (all ones for no metric). With no pairs, returns `h0 * D^-1 q`.
# `diagonal`, when given, is the initial inverse Hessian itself, one scale per
# coordinate (`_ctsem_lbfgs_diagonal!`), in place of the metric's shape times
# the latest pair's one secant ratio; or a whole matrix the caller supplied
# (`_ctsem_lbfgs`'s `inverse`), taken as it is, with or without pairs.
function _ctsem_lbfgs_hmul(M::CTSEMLBFGSMemory, q::AbstractVector, h0::Float64,
        dinv::Vector{Float64}, diagonal=nothing)
    k = length(M.S)
    k == 0 && return diagonal isa AbstractMatrix ? diagonal * q : h0 .* dinv .* q
    r = collect(Float64, q)
    a = Vector{Float64}(undef, k)
    @inbounds for i in k:-1:1
        a[i] = M.rho[i] * dot(M.S[i], r)
        r .-= a[i] .* M.Y[i]
    end
    if diagonal === nothing
        y = M.Y[k]
        gamma = dot(M.S[k], y) / sum(i -> y[i]^2 * dinv[i], eachindex(y))
        r .= gamma .* dinv .* r
    elseif diagonal isa AbstractMatrix
        r = diagonal * r
    else
        r .= diagonal .* r
    end
    @inbounds for i in 1:k
        b = M.rho[i] * dot(M.Y[i], r)
        r .+= (a[i] - b) .* M.S[i]
    end
    r
end

# `f(lo) + ... + f(hi)` in the order Base's `sum` adds terms it cannot
# vectorise: one at a time within blocks of 1024, the halves of a longer range
# added pairwise. `_ctsem_lbfgs_diagonal!`'s scale summed in this order while
# its `B` was boxed (untyped); typed, `sum` would vectorise from 16 terms on and
# round differently, and on a weakly identified model that last bit took the
# fit down another path to the same maximum. The fixed order keeps fits as they
# were.
function _ctsem_sum_inorder(f, lo::Int, hi::Int)
    hi == lo && return f(lo)
    if hi - lo < 1024
        v = f(lo) + f(lo + 1)
        for i in (lo + 2):hi
            v += f(i)
        end
        return v
    end
    mid = lo + (hi - lo) >> 1
    return _ctsem_sum_inorder(f, lo, mid) + _ctsem_sum_inorder(f, mid + 1, hi)
end

"""
    _ctsem_lbfgs_diagonal!(B, s, y, dinv)

The diagonal of the BFGS update of a diagonal Hessian approximation `B` by the
curvature pair `(s, y)` -- Gilbert & Lemaréchal's (1989) diagonal update, the
one M1QN3 uses for L-BFGS's initial matrix -- in place, returning the inverse
diagonal the two-loop recursion takes. `B === nothing` starts it at the
metric's shape `D`. Before each update `B` is rescaled by the Oren-Spedicato
factor `y'B^-1 y / s'y`, as theirs is, so its overall scale follows the latest
pair the way the scalar rule's does and only the shape across coordinates is
learned; from `D` that is the matrix the scalar rule would have used.

    B <- (y'B^-1 y / s'y) B
    B_i <- B_i - (B_i s_i)^2 / s'Bs + y_i^2 / s'y

Each term keeps `B_i` positive for a pair with `s'y > 0` (the first two are
`B_i (1 - B_i s_i^2 / s'Bs) >= 0`), so no coordinate can turn negative; a floor
at `1e-12` of the largest stops a coordinate the pairs never touch from
reaching zero. Why a diagonal: one secant ratio sets the scale of every
coordinate from the latest step, which is the stiffest curvature the step saw,
so on a model whose curvatures span orders of magnitude the weakly determined
coordinates take steps a thousandth of their own scale; the diagonal learns each
coordinate's scale from the pairs, as a per-parameter step size would, while the
two-loop recursion keeps the curvature between coordinates.
"""
function _ctsem_lbfgs_diagonal!(B0, s::Vector{Float64}, y::Vector{Float64},
        dinv::Vector{Float64})
    sy = dot(s, y)
    B = B0 === nothing ? [1 / d for d in dinv] : B0
    scale = _ctsem_sum_inorder(i -> y[i]^2 / B[i], 1, length(y)) / sy
    isfinite(scale) && scale > 0 && (B .*= scale)
    Bs = B .* s
    sBs = dot(s, Bs)
    if sBs > 0 && isfinite(sBs)
        @inbounds for i in eachindex(B)
            B[i] = B[i] - Bs[i]^2 / sBs + y[i]^2 / sy
        end
    end
    lo = 1e-12 * maximum(B; init=0.0)
    @inbounds for i in eachindex(B)
        (isfinite(B[i]) && B[i] > lo) || (B[i] = lo > 0 ? lo : 1.0)
    end
    B
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

`coreset`, when given, is a `CTSEMCoresetPlan`: once no batch is in force,
stretches of the run go on a coreset of the units, each checked by a full
gradient (`_ctsem_coreset_step!`), and the run returns the full objective's
value and gradient at its end.

`inverse`, when given, is the initial inverse Hessian, a whole matrix: the
curvature a caller already holds near `x0`, as a fit's covariance is for a
refit to slightly different data. The pairs correct it as they would the
learned diagonal, which it replaces (`diagonal` is then ignored). From an
inverse that is close, the first step is nearly the Newton step and what each
step predicts it gains is nearly what is left, which is what the gap rule reads
(`_ctsem_predicted_gain`); learned from the pairs alone, that prediction was a
hundredth of what was left on a 433-parameter model.

`jump`, a `CTSEMJump`, asks the line search to watch for a jump. A trial
whose excess over the linear model `alpha g's` does not shrink as the step does
-- a quarter of it left at a hundredth of an earlier trial's length -- may be a
discontinuity, where backtracking only reaches a step too short to leave the
point: a Laplace fit whose person's term jumped by 0.33 nats as its inner solve
changed where it stopped accepted steps below the resolution of the parameters
for hours, each after 40 backtracks
(CT-SEM/review/LAPLACE-shoulders-frozen-draws-2026-10-21.md). It may also be a
plateau the search gets past -- a run of refused trials, a profile point's
boundary -- so the search backtracks as before and the jump counts only if it
then finds no step, the retry with a fresh memory included, or finds one only
below the jump's scale three searches running. Then the run stops, the line
search failed and the jump recorded in `jump`.
An accepted step that leaves the point where it was fails the line search
too: it is no step.

`pairs`, a `CTSEMLBFGSMemory`, is the curvature to start from: what an earlier
run on the same coordinates ended with (`CTSEMLBFGSResult`'s `memory`), as a
rescued fit's stages hand on across a placement of the draws. It is copied,
its newest `memory` pairs kept; one for another number of parameters is
refused. The result's `memory` is what this run ended with.
"""
function _ctsem_lbfgs(fg!, x0::AbstractVector; memory::Integer=20,
        pairs::Union{Nothing,CTSEMLBFGSMemory}=nothing,
        metric=nothing, initial_alpha::Real=0.1, maxiter::Integer=1000,
        g_tol::Real=1e-8, f_tol::Real=0.0, x_tol::Real=0.0,
        callback=nothing, directional=nothing, batch=nothing,
        c1::Real=1e-4, maxbacktrack::Integer=40, iteration0::Integer=0,
        diagonal::Bool=false, nonmonotone::Real=0.0, gll::Integer=0,
        inverse::Union{Nothing,Matrix{Float64}}=nothing, scaled::Bool=false,
        structured=nothing, coreset=nothing,
        jump::Union{Nothing,CTSEMJump}=nothing, checkpoint=nothing)
    n = length(x0)
    x = collect(Float64, x0)
    dinv = metric === nothing ? ones(n) : begin
        d = collect(Float64, diag(metric))
        [isfinite(v) && v > 0 ? 1 / v : 1.0 for v in d]
    end
    # Pairs carried from another run, perhaps on an objective that has since
    # changed (a placement of the draws): the two-loop recursion stays positive
    # definite for any stored pairs, since `_ctsem_lbfgs_push!` keeps only
    # pairs with s'y > 0 and the learned diagonal is positive by construction,
    # so they still give a descent direction on the new objective. Where they
    # lead nowhere the line search fails and the retry below starts afresh.
    M = _ctsem_lbfgs_carry(pairs, memory, n)
    G = zeros(n)
    fcalls = 0; gcalls = 0
    sizes = Int[]; its = Int[]
    # The objective in force: the full one, or the current batch. The closures
    # here take what changes as arguments rather than capturing it: a captured
    # local assigned more than once is a `Core.Box`, untyped wherever it is
    # read, and an untyped trial point made inference walk the whole objective
    # again for `Any` -- most of what the package image held and what loading
    # it cost (review/PRECOMPILE-workload-2026-10-12.md).
    # A coreset stage, when one is in force, before either.
    evaluate!(b, F, Gout, y) = _ctsem_coreset_live(coreset) ?
        coreset.fg!(coreset.target, F, Gout, y) :
        b === nothing ? fg!(F, Gout, y) : _ctsem_batch_fg!(b, F, Gout, y)
    coreset === nothing || (coreset.g_tol = Float64(g_tol))
    f = evaluate!(batch, 0.0, G, x); fcalls += 1; gcalls += 1
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
    # Searches in a row that found a step only below a jump the watch saw.
    crept = 0
    # `diagonal`: the initial inverse Hessian learned per coordinate from the
    # pairs (`_ctsem_lbfgs_diagonal!`), `B` its inverse, rather than the
    # metric's shape times one secant ratio. `nonmonotone`: Zhang & Hager's
    # (2004) line search, which accepts a step against `C`, a running average
    # of the values the steps reached weighted by `eta`, rather than against the
    # last one -- so a step may rise above where it started while it stays
    # below the recent level, which lets a step cross a curved valley that a
    # monotone search would cut short. `eta = 0` makes `C` the last value: the
    # search is then exactly the monotone Armijo one.
    learn = diagonal && inverse === nothing
    B = learn && M.diagonal !== nothing ? M.diagonal : nothing
    Dinv = B === nothing ? inverse : 1 ./ B
    eta = clamp(Float64(nonmonotone), 0.0, 1.0)
    C = f; Q = 1.0; recent = [f]
    # `gll`: Grippo, Lampariello & Lucidi's (1986) rule instead -- a step is
    # judged against the worst of the last `gll` values the steps reached, the
    # acceptance sgd() uses (`_ctsem_sgd`). `recent` is set beside C.
    reference(recent, C) = gll > 0 ? maximum(@view recent[max(1, end - gll + 1):end]) : C
    while !stopped && !gconv && iteration < maxiter
        # `structured` (R's `lbfgs_scores`): a second model of the curvature,
        # whose step is taken while it predicts better than this one
        # (`CTSEMStructured`).
        sdir = structured === nothing ? nothing :
            _ctsem_structured_direction(structured, x, G)
        took_structured = sdir !== nothing
        s = took_structured ? sdir : -_ctsem_lbfgs_hmul(M, G, h0, dinv, Dinv)
        # Whether `s` comes from a model of the curvature -- the structured
        # model, pairs, an inverse the caller holds, or coordinates the caller
        # has scaled by it (`scaled`) -- rather than the gradient alone. Only
        # such a step predicts anything of what is left to gain.
        informed = took_structured || scaled || !isempty(M.S) || Dinv !== nothing
        # With no curvature pairs the step is the metric's alone, and at a start
        # where the transforms are flat the metric says a raw unit is worth
        # almost nothing -- so 0.1 in model units is an enormous raw step. One
        # raw unit at most, then: a long first step is how a fit once went from
        # raw 0 to raw 20.9 and ended where every transform is flat.
        if isempty(M.S) && !took_structured
            len = norm(s)
            len > 1 && (s .*= 1 / len)
        end
        dphi = dot(G, s)
        if (!(dphi < 0) || !isfinite(dphi)) && took_structured
            structured.active = false
            continue
        end
        if !(dphi < 0) || !isfinite(dphi)
            # Not a descent direction: the memory has gone bad. Start it again.
            isempty(M.S) && (lsfail = true; break)
            _ctsem_lbfgs_reset!(M); h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
            B = nothing; Dinv = inverse; C = f; Q = 1.0; recent = [f]
            continue
        end
        # The first trial carries the gradient too: most steps are accepted
        # whole, and then the iteration has cost one gradient and nothing else.
        alpha = 1.0
        xn = x .+ s
        Gn = similar(G)
        fn = evaluate!(batch, 0.0, Gn, xn); fcalls += 1; gcalls += 1
        have_gradient = true
        accepted = isfinite(fn) && fn <= reference(recent, C) + c1 * alpha * dphi
        k = 0
        # The trials' lengths and their excess over the linear model, for the
        # jump watch (`jump`); empty when nothing watches. A trial the route
        # refused comes back as a huge sentinel and is no measurement of the
        # objective, so it is not one of them. Not in a coreset stage, whose
        # objective weights a few units many times over: a spike there is one
        # unit's, scaled, and a search that fails in a stage leaves it for the
        # full objective, where the watch sees the jump if there is one.
        watching = jump !== nothing && !_ctsem_coreset_live(coreset)
        trial_alpha = Float64[]
        trial_excess = Float64[]
        found = false
        found_hi = NaN; found_lo = NaN; found_excess = NaN
        while !accepted && k < maxbacktrack
            if watching && !found && isfinite(fn) && fn < _CTSEM_REFUSED
                excess = fn - f - alpha * dphi
                if excess > 1e-8 * (1 + abs(f))
                    for (a, e) in zip(trial_alpha, trial_excess)
                        if a >= 100 * alpha && excess >= e / 4
                            found = true
                            found_hi = a; found_lo = alpha; found_excess = excess
                            break
                        end
                    end
                    if !found
                        push!(trial_alpha, alpha)
                        push!(trial_excess, excess)
                    end
                end
            end
            k += 1
            # Quadratic interpolation of phi(alpha), safeguarded to [0.1, 0.5].
            denom = 2 * (fn - f - dphi * alpha)
            trial = isfinite(fn) && denom > 0 ? -dphi * alpha^2 / denom : 0.5alpha
            alpha = clamp(trial, 0.1alpha, 0.5alpha)
            xn = x .+ alpha .* s
            fn = evaluate!(batch, 0.0, nothing, xn); fcalls += 1
            have_gradient = false
            accepted = isfinite(fn) && fn <= reference(recent, C) + c1 * alpha * dphi
        end
        # A step accepted on its value alone still needs a usable gradient.
        # When the gradient evaluation is refused (see `trial_fg!`), the step
        # is treated as too long: halved, and tried again with both.
        if accepted && !have_gradient
            evaluate!(batch, nothing, Gn, xn); gcalls += 1
            have_gradient = true
            while !all(isfinite, Gn)
                accepted = false
                k < maxbacktrack || break
                k += 1
                alpha *= 0.5
                xn = x .+ alpha .* s
                fn = evaluate!(batch, 0.0, Gn, xn); fcalls += 1; gcalls += 1
                accepted = isfinite(fn) && fn <= reference(recent, C) + c1 * alpha * dphi &&
                    all(isfinite, Gn)
                accepted && break
            end
        end
        # A step that leaves the point where it was is none, however the
        # Armijo test reads it: below the parameters' resolution `fn == f`, and
        # `f + c1 alpha dphi` rounds to `f`.
        if accepted && xn == x
            accepted = false
        end
        if !accepted && _ctsem_coreset_live(coreset)
            # A coreset whose line search finds nothing has stopped leading
            # anywhere the full objective goes: check it and leave the stage.
            hm = let h0 = h0, Dinv = Dinv
                q -> _ctsem_lbfgs_hmul(M, q, h0, dinv, Dinv)
            end
            r = _ctsem_coreset_check!(coreset, x, f, G, hm,
                _ctsem_lbfgs_scale(Dinv, dinv); leave=true)
            r.act === :revert && _ctsem_lbfgs_drop!(M, coreset.pairs)
            x, f, G, nev = _ctsem_coreset_switch!(coreset, r, evaluate!, batch)
            fcalls += nev; gcalls += nev
            gconv = maximum(abs, G; init=0.0) <= g_tol
            C = f; Q = 1.0; recent = [f]
            continue
        end
        # A jump the watch saw counts only once backtracking has confirmed it:
        # the search, retried as it always is, found no step at all (below), or
        # found one only below the jump's scale for the third search running
        # -- creeping along a spike's top. A plateau backtracking gets past is
        # no jump, and the search goes on as it always did.
        if accepted && found && alpha < found_lo
            crept += 1
            if crept >= 3
                _ctsem_jump_record!(jump, x, s, found_hi, found_lo, found_excess)
                # A route that takes the jump in (`_ctsem_anytime!`) changes
                # its objective under the point, which stays where it is: the
                # step is given up and the objective read again there.
                if checkpoint !== nothing && checkpoint(x, 0.0, :jump)
                    jump.seen = false; crept = 0; retried = false
                    f = evaluate!(batch, 0.0, G, x); fcalls += 1; gcalls += 1
                    gconv = maximum(abs, G; init=0.0) <= g_tol
                    C = f; Q = 1.0; recent = [f]
                    continue
                end
                lsfail = true
                break
            end
        elseif accepted
            crept = 0
        end
        if !accepted && took_structured
            # The structured model's step found nothing: plain L-BFGS from
            # here, until the structured model predicts better again.
            structured.active = false
            continue
        end
        if !accepted
            # A stale memory is the usual cause; drop it once, then give up.
            if !retried && !isempty(M.S)
                retried = true
                _ctsem_lbfgs_reset!(M)
                h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
                B = nothing; Dinv = inverse; C = f; Q = 1.0; recent = [f]
                continue
            end
            found && _ctsem_jump_record!(jump, x, s, found_hi, found_lo, found_excess)
            if found && checkpoint !== nothing && checkpoint(x, 0.0, :jump)
                jump.seen = false; crept = 0; retried = false
                f = evaluate!(batch, 0.0, G, x); fcalls += 1; gcalls += 1
                gconv = maximum(abs, G; init=0.0) <= g_tol
                C = f; Q = 1.0; recent = [f]
                continue
            end
            lsfail = true
            break
        end
        retried = false
        directional === nothing ||
            (directional.dphi0 = informed ? alpha * dphi : NaN)
        # What the step was predicted to gain, as the gap rule reads it
        # (`_ctsem_predicted_gain`), for a coreset stage to stop where the full
        # sum would.
        stepgain = informed ? abs(alpha * dphi) / 2 : Inf
        step = xn .- x
        iteration += 1
        fconv = abs(fn - f) <= f_tol * abs(fn)
        xconv = maximum(abs, step; init=0.0) <= x_tol
        dG = Gn .- G
        structured === nothing || _ctsem_structured_update!(structured, x, xn,
            G, Gn, f, fn, step, alpha, took_structured)
        pushed = _ctsem_lbfgs_push!(M, step, dG)
        pushed && _ctsem_coreset_live(coreset) && (coreset.pairs += 1)
        if pushed && learn
            B = _ctsem_lbfgs_diagonal!(B, step, dG, dinv)
            Dinv = 1 ./ B
        end
        x = xn; f = fn; G = Gn
        push!(recent, f)
        Qn = eta * Q + 1
        C = (eta * Q * C + f) / Qn; Q = Qn
        gconv = maximum(abs, G; init=0.0) <= g_tol
        # A coreset stage's value is shown on the full objective's scale, from
        # where the stage began (`offset`), so the trace and the progress line
        # read on across the switches.
        shown = _ctsem_coreset_live(coreset) ? f + coreset.offset : f
        stopped = callback !== nothing && callback(CTSEMIterate(
            iteration0 + iteration, shown, maximum(abs, G; init=0.0))) === true
        # The route's own checkpoint, on the full objective only: given the
        # step's predicted gain it may screen its units and change itself,
        # and the objective is then read again at the point, the pairs kept
        # (they stay positive definite, as carried pairs do).
        if checkpoint !== nothing && !stopped && batch === nothing &&
                !_ctsem_coreset_live(coreset) && checkpoint(x, stepgain, :step)
            f = evaluate!(batch, 0.0, G, x); fcalls += 1; gcalls += 1
            gconv = maximum(abs, G; init=0.0) <= g_tol
            C = f; Q = 1.0; recent = [f]
        end
        # A growing batch changes the objective under the optimiser. The
        # curvature pairs are kept -- the batch objective is scaled to the full
        # data, so they estimate the same curvature -- but no pair spans the
        # change, and the value and gradient are re-read on the new batch. Once
        # it is the whole data the batch steps aside and the caller's own
        # objective is used from then on.
        if batch !== nothing && !stopped &&
                _ctsem_batch_step!(batch, x, G, let h0 = h0, Dinv = Dinv
                    q -> _ctsem_lbfgs_hmul(M, q, h0, dinv, Dinv)
                end, iteration)
            if _ctsem_batch_full(batch)
                sizes = batch.sizes; its = batch.iterations
                batch = nothing
            end
            f = evaluate!(batch, 0.0, G, x); fcalls += 1; gcalls += 1
            gconv = maximum(abs, G; init=0.0) <= g_tol
            # The reference level belongs to the objective it averaged.
            C = f; Q = 1.0; recent = [f]
        end
        # Coreset stages, once no batch is in force: the step may begin one,
        # check one, or end one (`_ctsem_coreset_step!`). The pairs are kept
        # whichever way it switches; the value, gradient and reference level
        # are those of the objective now in force.
        if coreset !== nothing && !stopped && batch === nothing
            hm = let h0 = h0, Dinv = Dinv
                q -> _ctsem_lbfgs_hmul(M, q, h0, dinv, Dinv)
            end
            r = _ctsem_coreset_step!(coreset, x, f, G, hm,
                _ctsem_lbfgs_scale(Dinv, dinv); stepgain=stepgain)
            if r.act !== :none
                r.act === :revert && _ctsem_lbfgs_drop!(M, coreset.pairs)
                x, f, G, nev = _ctsem_coreset_switch!(coreset, r, evaluate!, batch)
                fcalls += nev; gcalls += nev
                gconv = maximum(abs, G; init=0.0) <= g_tol
                C = f; Q = 1.0; recent = [f]
            end
        end
        (fconv && f_tol > 0) && break
        (xconv && x_tol > 0) && break
    end
    if batch !== nothing
        sizes = batch.sizes; its = batch.iterations
    end
    # A run that ends inside a coreset stage ends on the full objective.
    if _ctsem_coreset_live(coreset)
        hm = let h0 = h0, Dinv = Dinv
            q -> _ctsem_lbfgs_hmul(M, q, h0, dinv, Dinv)
        end
        r = _ctsem_coreset_check!(coreset, x, f, G, hm, _ctsem_lbfgs_scale(Dinv, dinv);
            leave=true)
        r.act === :revert && _ctsem_lbfgs_drop!(M, coreset.pairs)
        x, f, G, nev = _ctsem_coreset_switch!(coreset, r, evaluate!, batch)
        fcalls += nev; gcalls += nev
        gconv = maximum(abs, G; init=0.0) <= g_tol
    end
    M.diagonal = learn ? B : nothing
    CTSEMLBFGSResult(x, f, G, iteration, fcalls, gcalls, gconv, fconv, xconv,
        lsfail, stopped, sizes, its, M)
end

"""
    _ctsem_sgd(fg!, x0; ...)

The early phase for large, ill-conditioned models: ctsem 3.11's `sgd()` (the
stan path's `stochastic = TRUE`), reduced to the mechanisms that measured as
carrying it. One gradient per iteration, no line search.

  * A step per coordinate, adapted by how often that coordinate's gradient
    changes sign: two exponential averages of the flip rate -- of the mean of
    consecutive gradients (aimed at 0.5) and of the momentum (aimed at 0.1,
    pushed back hard above it) -- each scaling the step towards its target, by
    at most a fifth an iteration. A trust radius per coordinate: a flip says the
    last step overshot that coordinate's optimum, a run of one sign says it fell
    short. This is the information a diagonal curvature estimate carries, read
    from signs, which stay meaningful where the curvature measured at one point
    says little about the next.
  * A direction that is a moving average of the gradient, the momentum (its
    weight rises from 0.6 to 0.8 over the first 20 iterations), which carries the
    steps along a curved valley that a single gradient crosses.
  * Its magnitude compressed to `sign(m) sqrt|m|`: the gradients of
    exponential and logistic transforms span orders of magnitude and change by
    orders as a coordinate moves, faster than a multiplicative step can follow,
    and the square root halves that range. Each coordinate's step is capped at
    `cap` raw units.
  * Non-monotone acceptance, Grippo, Lampariello & Lucidi's (1986) rule with a
    slack: a step is kept if the objective stays below the worst of the last 20
    values plus one; the k-th rejection in a row divides every step by e^k.

Measured against `sgd()` on the 715-parameter model it was written for (dev2,
500 gradients from one start), removing the per-coordinate steps cost 2500 log
posterior, the compression 1900, the non-monotone acceptance 300 and the
momentum 200. `sgd()`'s other mechanisms -- a look-ahead term, a global step
rule from how often the objective falls, a pull back to the best point, random
walks on its momentum and target, an inflection rule, a warm-up check -- were
each worth nothing measurable there and on a 116-parameter analogue, and
together they were worth nothing (the analogue, two data sets: the same end
points to 0.4); they are not here. A single sign rule in place of the two
averages (Rprop's 1.2 and 0.5) was tried and is worse with momentum: the
momentum's sign rarely flips, the steps compound, and a third of the trials are
rejected.

The first direction is the gradient's sign at the size `sqrt|f0|`, as `sgd()`'s
was, which makes the first steps a fixed small move in every coordinate rather
than one proportional to gradients that can be ten thousand at a cold start.

Stops at `maxiter`, when the callback says so, after `maxreject` rejections in
a row, when the best value has improved by less than `itertol` an iteration
over the last 30 (`sgd()`'s own rule), or, with `progress > 0`, when those 30
iterations gained less than `progress` of everything gained since the start --
where steps of the curvature's own size, L-BFGS's, do better. Minimises through
the same `fg!` contract as `_ctsem_lbfgs` and returns the same result type, at
the best point reached; `iteration0` offsets the iteration numbers the callback
sees.
"""
function _ctsem_sgd(fg!, x0::AbstractVector; maxiter::Integer=1000,
        step0::Real=1e-3, cap::Real=0.5, window::Integer=20,
        maxreject::Integer=20, itertol::Real=1e-3, progress::Real=0.0,
        callback=nothing, iteration0::Integer=0)
    n = length(x0)
    x = collect(Float64, x0)
    G = zeros(n)
    f = fg!(0.0, G, x)
    fcalls = 1; gcalls = 1
    stopped = callback !== nothing && iteration0 == 0 &&
        callback(CTSEMIterate(0, f, maximum(abs, G; init=0.0))) === true
    warmup = 20
    # Ascent quantities, for the maximised objective -f. `g` holds the previous
    # gradient when a point is accepted; before the first it is the scaled sign
    # direction, as sgd()'s was.
    first = sqrt(abs(f))
    gsmooth = [-sign(v) * first for v in G]
    gmid = copy(gsmooth)
    g = copy(gsmooth)
    step = fill(Float64(step0), n)
    groughness = fill(0.5, n)
    gsmoothroughness = fill(0.1, n)
    # Refs: `accept!` updates both and the loop below reads `gmemory2` (a
    # captured local assigned more than once is a `Core.Box`, untyped).
    gmemory = Ref(0.6)
    gmemory2 = Ref(0.0)
    history = Float64[]
    bestx = copy(x); bestf = f; bestG = copy(G)
    xn = similar(x); Gn = similar(G)
    oldgmid = similar(x); oldgsmooth = similar(x)
    mod(r, t) = t / (r + t) - 0.5
    iteration = 0; rejected = 0
    # The start is the first accepted point; `accept!` is everything sgd() does
    # with one.
    accept! = function (fnew, gnew, i)
        push!(history, fnew)
        gmemory2[] = gmemory[] * min(i / warmup, 1)^(1 / 8)
        rm2 = 0.9 * min(i / warmup, 1)^(1 / 8)
        oldgmid .= gmid
        oldgsmooth .= gsmooth
        @inbounds for k in 1:n
            gmid[k] = (g[k] + gnew[k]) / 2
            gsmooth[k] = gsmooth[k] * gmemory2[] + (1 - gmemory2[]) * gnew[k]
            groughness[k] = groughness[k] * rm2 + (1 - rm2) *
                (sign(gmid[k]) != sign(oldgmid[k]))
            gsmoothroughness[k] = gsmoothroughness[k] * rm2 + (1 - rm2) *
                (sign(gsmooth[k]) != sign(oldgsmooth[k]))
            if i > warmup
                a = mod(gsmoothroughness[k], 0.1)
                step[k] *= 1 + 0.4 * sign(a) * a^4 / 0.5^4
                step[k] *= 1 + 0.12 * mod(groughness[k], 0.5)
            end
            g[k] = gnew[k]
        end
        # The momentum's weight settles at 0.8 once warmed up.
        (i > 25 && i % 20 == 0) && (gmemory[] = clamp(gmemory[], 0.8, 0.95))
        return nothing
    end
    accept!(-f, -G, 1)
    iteration = 1
    while !stopped && iteration < maxiter
        i = iteration + 1
        accepted = false
        tries = 0
        local fn
        while !accepted
            tries += 1
            tries > maxreject && break
            _ctsem_interrupt_check()
            @inbounds for k in 1:n
                d = step[k] * sign(gsmooth[k]) * sqrt(abs(gsmooth[k]))
                xn[k] = x[k] + clamp(d, -cap, cap)
            end
            fn = fg!(0.0, Gn, xn); fcalls += 1; gcalls += 1
            lpn = -fn
            # sgd()'s slack: one, and one more for every try so far.
            bar = minimum(@view history[max(1, end - window + 1):end]) - tries - 1
            accepted = isfinite(lpn) && all(isfinite, Gn) && lpn > bar
            if !accepted
                if i > warmup
                    @inbounds for k in 1:n
                        gsmooth[k] = gsmooth[k] * gmemory2[]^2 + (1 - gmemory2[]^2) * g[k]
                    end
                end
                step ./= exp(tries)
            end
        end
        accepted || (rejected = tries; break)
        iteration = i
        x .= xn; f = fn; G .= Gn
        accept!(-fn, -Gn, i)
        if f < bestf
            bestf = f; bestx .= x; bestG .= G
        end
        stopped = callback !== nothing && callback(CTSEMIterate(
            iteration0 + iteration, f, maximum(abs, G; init=0.0))) === true
        # sgd()'s own rule, and the relative one: the best value over the last
        # 30 iterations against the best before them.
        if iteration > 30 && maximum(@view history[end-29:end]) == maximum(history)
            before = maximum(@view history[1:end-30])
            gained = maximum(history) - before
            (gained / 30 < itertol && gained > 0) && break
            total = maximum(history) - history[1]
            (progress > 0 && total > 0 && gained < progress * total) && break
        end
    end
    CTSEMLBFGSResult(bestx, bestf, bestG, iteration, fcalls, gcalls,
        maximum(abs, bestG; init=0.0) <= 0, false, false, rejected > 0,
        stopped, Int[], Int[], CTSEMLBFGSMemory(0))
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
# gradient, which is what makes the separation exact. (`priors = 'randomCorr'`
# puts a prior on random-effect correlations even in a maximum likelihood fit,
# and refusing a prior would refuse those fits.) A sampled TI predictor's
# imputation density is also global
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

"""A CTSEMObjective over some of its subjects; the subject objectives are
self-contained, and the sampled TI predictor cells are those subjects' own,
renumbered with them."""
function _ctsem_subset_objective(o::CTSEMObjective, idx::AbstractVector{<:Integer})
    renumber = Dict(Int(s) => j for (j, s) in enumerate(idx))
    keep = [k for k in eachindex(o.ti_missing_subject) if haskey(renumber, o.ti_missing_subject[k])]
    CTSEMObjective(o.params, o.subject_objectives[idx], nothing,
        o.prior_index, o.prior_scale, o.prior_weight,
        o.ti_missing_parameter[keep], o.ti_missing_mu[keep], o.ti_missing_sigma[keep],
        Int[renumber[o.ti_missing_subject[k]] for k in keep],
        isempty(o.subject_weights) ? Float64[] : o.subject_weights[idx], o.corbound,
        o.penalty_index, o.penalty_precision)
end
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
        u.blocks[keep], u.active[keep], u.inactive[keep])
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

# ------------------------------------------------------------------ coreset stages
#
# On the Laplace route a gradient is a sum over units, each an inner mode solve,
# and along an optimiser's path the units' scores vary in about as many
# dimensions as there are parameters, however many units there are. So a few of
# them with non-negative weights (`CTSEMCoreset`) can follow how the full
# gradient changes over a stretch of the path, and L-BFGS can run on their
# weighted sum at their share of the cost (review/OPTIM-coreset-2026-10-21.md).
#
# Every full gradient on the way already holds every unit's score: the trial
# takes them from the same reverse sweep (`CTSEMScoreStore`). The last few of
# those points are the training set. A stage is ordinary L-BFGS on the coreset
# objective, with the curvature pairs kept across the switch both ways, and one
# full gradient checks each stretch of it; the full sum finishes, and the Newton
# finish and its Hessian, which certify the estimate and give its standard
# errors, never see a coreset.
#
# Four rules, each for a reason that does not depend on the model:
#
# - A stage's gradient is exact where it begins. The coreset's objective is
#   anchored there (`CTSEMCoresetTarget`): its gradient at the anchor is
#   replaced by the full one and the difference carried as a linear term, so
#   that `g(x) = g_full(x_a) + g_c(x) - g_c(x_a)`, a control variate as in SVRG
#   (Johnson & Zhang, 2013). A coreset then has only to follow how the gradient
#   changes, which is the curvature, and its error grows from zero with the
#   distance from the anchor instead of standing at a floor. To first order a
#   stage's maximum misses the full one by the coreset's relative curvature
#   error times the gradient it began with, so stages chained from check to
#   check close on the full maximum, as inexact Newton steps do. Unanchored, a
#   coreset's error stood while the gradient shrank: a stage led to the
#   coreset's own maximum, and none passed near the full one.
# - A coreset is tested where it was not fitted: its change in gradient over
#   every step between training points, each held out once, must be within
#   `tol` of the full change, relative, in L-BFGS's inverse-curvature norm, the
#   norm in which a gradient error becomes an error in the step and in what it
#   gains (`_CTSEM_CORESET_TOL`) -- the mean of those errors plus its standard
#   error. Every step, because a change in gradient along one direction says
#   nothing of the others: tested on the newest step alone, coresets of one to
#   five units of 250 passed and two stages in three made the full objective
#   worse; tested on the newest quarter of the steps, coresets of one to five
#   units of 400 passed early on, while there were two or three steps, and a
#   random effect's sd went with them to its floor. The standard error,
#   because a few small errors are no evidence: with one step there is no
#   test, and with few the coreset must do clearly better. And the steps
#   need not have moved every parameter, so the coreset must also carry, in
#   every parameter, enough of the units' spread of scores where the stage
#   begins -- the diagonal of their covariance, which estimates the
#   information (`_ctsem_coreset_shortfall`) -- that a stage's step there is
#   off by at most `tol`, as the steps' test asks of the steps. Too little
#   makes a step too long; too much only shortens it, so the test is one-sided.
#   At the contraction boundary instead (half the information), coresets
#   chosen for their spread sat just above it, steps along two sds ran up to
#   twice too long stage after stage, and both collapsed, 20 below the
#   optimum. A random effect's sd has its curvature
#   from the spread of the units and from little else: on bigre from refit
#   starts, coresets of 3-8 of 400 units passed the steps' test, and their
#   stages took an sd to its floor, 8 and 12 below the full sum's optimum.
#   Asked to be within the tolerance both ways, the coresets took most of the
#   units and saved nothing. Entering and
#   carrying on alike, and always the smallest coreset that passes. An
#   accuracy and not a price: the stages share one curvature memory, and a
#   coreset whose gradients differ more feeds it pairs from another function.
#   Tested in sample, or priced by a model of what a coarser coreset would
#   still save, coresets were let in whose stages made the full objective
#   worse, and runs cost more. And it is trusted only where it was fitted: a
#   stage is checked once it is farther from where it began than the farthest
#   training point, in L-BFGS's scaled metric, or once it has moved any
#   parameter farther than the training points reach in that parameter, and
#   each check is a training point there, so the region grows as fast as it is
#   confirmed. A few units stand for a population's spread only near where
#   they were chosen; unbounded, stages on coresets of 4-6 of 500 and 37-338 of
#   400 units walked a random effect's sd to the far end of its transform, and
#   the fits ended on maxima 157 and 8 log units below the full sum's. Bounded
#   by distance alone, a stage could spend the whole reach of steps that moved
#   every parameter a little along one sd, and on bigre from refit starts sds
#   collapsed in about five runs of eight, 8-20 below the optimum the full sum
#   reached from every one.
# - A stage spends only what it can save. `n` coreset iterations at a share `s`
#   of a full evaluation, plus the full check and the full iteration after the
#   stage, beat `n` full iterations while `n s + 2 < n`, so a stage planned for
#   `n` iterations may use units costing up to `1 - 2/n` of the work
#   (`_ctsem_coreset_budget`), and it carries on only if it paid: its iterations
#   times the full objective's gain over the coreset's own against their share
#   of the work and the check. The first stage is planned for 10 iterations and
#   each that pays for twice the last, the cheapest search for how long a
#   coreset stays good: its checks cost at most twice those of the longest stage
#   that would have paid.
# - A stage that led somewhere worse leaves nothing behind. The point goes back
#   to where it began, and the curvature pairs it gave L-BFGS go with it; and the
#   stall watch judges only windows of full iterations after the last stage
#   (`ctsem_optimize`), since a window holding a stage measures the coreset's
#   path, and one holding a stage sent back reads as no progress.
#
# Entering: after a full iteration, the steps between training points are
# split into up to four contiguous folds; for each, coresets of growing size,
# from one unit, are fitted to reproduce the change in gradient over the other
# steps, by forward selection with non-negative weights, and tried on the
# fold's, and the first that passes within the budget is its size. The coreset
# is the largest of those sizes, fitted again on every step, if the errors met
# pass together. At a check, the check is a training point whatever it says,
# and a stage carried on from it is chosen the same way, so that its tests
# hold the change over the stage just run. Within a stage, a check comes when
# the stage has run as planned, or
# when what L-BFGS's model says is left to gain on the coreset is no more than
# its error -- its held-out relative error, squared, times what was left to gain
# where the stage began -- since what it would go on to gain is then its error.
# No stage is carried on once the full objective gained less than the fit's
# gap tolerance an iteration over the stage just run -- the bar a full step's
# predicted gain meets in the gap rule, read from what the stage measurably
# did: the end is the full sum's, and its stopping rules'. And a stage whose
# step is predicted to gain less than the tolerance, where the gap rule would
# stop a fit on the full sum, is checked and left with nothing carried on. Carried on while
# L-BFGS's model of what was left exceeded the tolerance, stages climbed on
# along AnomAuth's flat directions where the full sum's own gap rule had
# stopped it, 8-14 log units up an over-credited Laplace term, and the
# certification paid twice the fit's cost to come back. Nor is one entered
# after a full iteration that gained less than the tolerance: stages are for
# climbing. Entered after every full iteration on bigre's collapsed-sd plateau
# from a refit start, stage after stage drifted along the sd, flat for the full
# objective, wherever its coreset's error pointed, from raw -2 to -9.5, and the
# fit ended there, 7.8 below the optimum the full sum reached.
#
# Looking for a coreset backs off: after a search that found none, the next
# waits twice as many full iterations as the last wait, and one that finds one
# starts the count again, so a run where nothing passes searches about `log2`
# of its full iterations' times. Counted in iterations, not timed: decisions
# that read the clock made a fit's path depend on the machine's load, and two
# identical fits at one core ended at different points.

"""
    CTSEMCoresetTarget(laplace, coreset, anchor, correction)

What a coreset stage optimises: `laplace`'s log prior plus the coreset's
weighted unit terms, through the Laplace route's own trial
(`_ctsem_optimise_trial`, with its validity rule), solving only the coreset's
units on the pool a full evaluation uses; plus `correction' (x - anchor)`,
where `correction` is the full gradient less the coreset's at `anchor`
(maximised), so that the stage's gradient there is the full one.
`_ctsem_coreset_switch!` sets `correction` once it has read the coreset's
gradient at the anchor.
"""
struct CTSEMCoresetTarget
    laplace::CTSEMLaplaceObjective
    coreset::CTSEMCoreset
    anchor::Vector{Float64}
    correction::Vector{Float64}
end

function _ctsem_optimise_trial(t::CTSEMCoresetTarget, x, want_gradient::Bool,
        gradient_method, limit::Real, log)
    trial = _ctsem_optimise_trial(t.laplace, x, want_gradient, gradient_method,
        limit, log; coreset=t.coreset)
    trial.valid || return trial
    e = trial.evaluated
    shift = 0.0
    @inbounds for i in eachindex(t.anchor)
        shift += t.correction[i] * (x[i] - t.anchor[i])
    end
    value = e.value + shift
    gradient = e.gradient === nothing ? nothing : e.gradient .+ t.correction
    return (evaluated=merge(e, (value=value, gradient=gradient)),
        valid=isfinite(value))
end

"""The first stage's planned length; a stage that pays plans the next twice as
long."""
const _CTSEM_CORESET_FIRST_STAGE = 10

"""The relative error a coreset may have in the change of the gradient, in
L-BFGS's inverse-curvature norm; see the head of this section for why it is an
accuracy and not a price. A stage run to its maximum leaves about its square
of what was left to gain where it began, so stages contract towards the full
maximum while it is below one; 0.3 leaves about a tenth, with room for the
error of measuring it along one step."""
const _CTSEM_CORESET_TOL = 0.3

"""The folds the steps between training points are split into to choose a
coreset (`_ctsem_coreset_choose`): every step is held out once, and a search
fits at most this many paths and one more."""
const _CTSEM_CORESET_FOLDS = 4

"""
    CTSEMCoresetPlan

The state of the coreset stages in one L-BFGS run (`_ctsem_lbfgs`'s
`coreset`): the full objective and the caller's trial (`fg!(target, F, G, x)`),
the store the full trials leave their units' scores in, each unit's sampling
weight (which its scores carry) and share of a full evaluation's work, the
training points (each point, its units' scores without the prior, and the full
gradient of the maximised objective there), the stage in force -- its target,
where it began on the full objective (`x0`, `f0`, `G0`, minimised), the
coreset's value there less the full one (`offset`), its length so far and as
planned, and the coreset's error, in objective units, `e' H e / 2`, where the
stage began -- the full iterations since the last search and how many to wait
before the next, the record, how many curvature pairs the stage in force has
given L-BFGS, which a stage that is sent back takes with it, the region its
coreset was fitted on: the farthest training point's distance from where the
stage began, in the metric `metric` scales by; and the held-out relative error
of the last coreset chosen, squared (`rel`), the fit's gap tolerance: no
stage is carried on from one that gained less than it an iteration (`gap_tol`,
set by `ctsem_optimize`); per parameter, the farthest the training points lie
from where the stage began (`reach`), which the stage may not pass; and the
accepted iterations so far, full or not (`seen`), with the one each stage
began after (`starts`), which places the stages on the run's trace; and the
full objective's value at the last full point (`f_last`, minimised), from
which a full iteration's gain is measured.
"""
mutable struct CTSEMCoresetPlan
    full::CTSEMLaplaceObjective
    fg!::Any
    store::CTSEMScoreStore
    sampling::Vector{Float64}
    cost::Vector{Float64}
    keep::Int
    g_tol::Float64
    X::Vector{Vector{Float64}}
    S::Vector{Matrix{Float64}}
    gfull::Vector{Vector{Float64}}
    live::Bool
    target::Union{Nothing,CTSEMCoresetTarget}
    x0::Vector{Float64}
    f0::Float64
    G0::Vector{Float64}
    offset::Float64
    iterations::Int
    stretch::Int
    error::Float64
    since::Int
    wait::Int
    sizes::Vector{Int}
    shares::Vector{Float64}
    entered::Vector{Int}
    stage_iterations::Vector{Int}
    bought::Vector{Float64}
    total_iterations::Int
    checks::Int
    reverted::Int
    attempts::Int
    select_seconds::Float64
    pairs::Int
    radius::Float64
    metric::Vector{Float64}
    rel::Float64
    gap_tol::Float64
    reach::Vector{Float64}
    seen::Int
    starts::Vector{Int}
    f_last::Float64
end

"""
    _ctsem_coreset_plan(objective, fg!, store)

Coreset stages for `objective`'s L-BFGS run, or `nothing` where they cannot
apply: a route other than Laplace, units whose scores the trial does not keep,
or a single unit. The training set holds the last 30 full gradients, fewer
where their scores would pass 2.5e7 numbers.
"""
_ctsem_coreset_plan(objective, fg!, store) = nothing
function _ctsem_coreset_plan(L::CTSEMLaplaceObjective, fg!, store)
    store === nothing && return nothing
    N, p = size(store.S)
    N >= 2 || return nothing
    w = _laplace_unit_weights(L)
    cost = w ./ sum(w)
    keep = clamp(floor(Int, 2.5e7 / (N * p)), 3, 30)
    sampling = Float64[_laplace_unit_weight(L, U) for U in 1:N]
    return CTSEMCoresetPlan(L, fg!, store, sampling, cost, keep, 0.0,
        Vector{Float64}[], Matrix{Float64}[], Vector{Float64}[], false, nothing,
        Float64[], NaN, Float64[], 0.0, 0, _CTSEM_CORESET_FIRST_STAGE, 0.0, 0,
        0, Int[], Float64[], Int[], Int[], Float64[], 0, 0, 0, 0, 0.0, 0,
        0.0, Float64[], 0.0, 0.0, Float64[], 0, Int[], NaN)
end

_ctsem_coreset_live(c) = c !== nothing && c.live

"""The largest share of a full evaluation's work a stage planned for `n`
iterations can spend on each: with its check and the full iteration after it,
`n` iterations at share `s` cost less than `n` full ones only while
`s < 1 - 2/n`."""
_ctsem_coreset_budget(n::Integer) = 1 - 2 / n

"""What a coreset costs, as a share of a full evaluation's work."""
_ctsem_coreset_share(c::CTSEMCoresetPlan, units) =
    sum(c.cost[U] for U in units; init=0.0)

"""Keep the full gradient at `x` (`G` minimised) and its units' scores, which
the trial just left in the store, as a training point."""
function _ctsem_coreset_record!(c::CTSEMCoresetPlan, x::AbstractVector, G::AbstractVector)
    st = c.store
    (st.valid && st.x == x) || return false
    push!(c.X, collect(Float64, x))
    push!(c.S, copy(st.S))
    push!(c.gfull, -collect(Float64, G))
    while length(c.S) > c.keep
        popfirst!(c.X); popfirst!(c.S); popfirst!(c.gfull)
    end
    return true
end

"""
    _ctsem_nnls(AtA, Atb)

Lawson and Hanson's non-negative least squares, from the normal equations:
the `x >= 0` minimising `|A x - b|^2`, given `A'A` and `A'b`. Each passive set
is solved by a Cholesky factorization of its block, with a ridge of `1e-12` of
the diagonal's mean so that collinear units do not stop it.
"""
function _ctsem_nnls(AtA::AbstractMatrix, Atb::AbstractVector)
    n = length(Atb)
    x = zeros(n)
    n == 0 && return x
    P = falses(n)
    w = collect(Float64, Atb)
    tol = 1e-10 * max(1.0, maximum(abs, w))
    ridge = 1e-12 * max(sum(AtA[i, i] for i in 1:n) / n, eps())
    solve(idx) = begin
        B = AtA[idx, idx]
        for i in eachindex(idx)
            B[i, i] += ridge
        end
        F = cholesky!(Symmetric(B); check=false)
        issuccess(F) ? F \ Atb[idx] : pinv(AtA[idx, idx]) * Atb[idx]
    end
    for _ in 1:(3n)
        j = 0; best = tol
        for i in 1:n
            (!P[i] && w[i] > best) && (best = w[i]; j = i)
        end
        j == 0 && break
        P[j] = true
        z = zeros(n)
        for _ in 1:(3n)
            idx = findall(P)
            z = zeros(n)
            z[idx] = solve(idx)
            all(z[i] > 0 for i in idx) && break
            alpha = Inf
            for i in idx
                z[i] <= 0 && (alpha = min(alpha, x[i] / (x[i] - z[i])))
            end
            x .+= alpha .* (z .- x)
            for i in idx
                x[i] <= 1e-12 && (P[i] = false; x[i] = 0.0)
            end
        end
        x = z
        w = Atb .- AtA * x
    end
    return x
end

"""
    _ctsem_coreset_information(S, scale)

Each parameter's spread of the units' scores (`S`, units by parameters) about
their mean `m`, `sum_U (s_U - m)^2` -- the diagonal of the score covariance,
which estimates the information -- with the means and the scores, and which
parameters count: those whose spread, in L-BFGS's metric (`scale`), is under
`1e-8` of the largest are left out, since no unit moves them.
"""
function _ctsem_coreset_information(S::AbstractMatrix, scale::AbstractVector)
    N, p = size(S)
    full = zeros(p)
    mean = zeros(p)
    @inbounds for i in 1:p
        m = 0.0
        for U in 1:N
            m += S[U, i]
        end
        m /= N
        d = 0.0
        for U in 1:N
            d += (S[U, i] - m)^2
        end
        full[i] = d
        mean[i] = m
    end
    scaled = [full[i] * max(scale[i], 0.0) for i in 1:p]
    top = maximum(scaled; init=0.0)
    return (spread=full, mean=mean, S=S,
        counts=[top > 0 && v > 1e-8 * top for v in scaled])
end

"""
    _ctsem_coreset_shortfall(S, units, w, info)

How far a coreset's information falls short, at worst over the parameters:
`(d - d_c) / d_c`, where `d` is the units' spread of scores
(`_ctsem_coreset_information`, `info`) and `d_c` the coreset's weighted spread
about its weighted mean. To first order and in each parameter alone, a stage's
step is off by that factor: the stages close on the full maximum only while it
is below one, and `_ctsem_coreset_choose` asks that it be within
`_CTSEM_CORESET_TOL`. Carrying more only shortens a step. `Inf` where the
coreset carries none.
"""
function _ctsem_coreset_shortfall(S::AbstractMatrix, units, w::AbstractVector, info)
    W = sum(w)
    (W > 0 && isfinite(W)) || return Inf
    worst = -Inf
    @inbounds for i in eachindex(info.spread)
        info.counts[i] || continue
        mc = 0.0
        for (j, U) in enumerate(units)
            mc += w[j] * S[U, i]
        end
        mc /= W
        dc = 0.0
        for (j, U) in enumerate(units)
            dc += w[j] * (S[U, i] - mc)^2
        end
        dc > 0 || return Inf
        worst = max(worst, (info.spread[i] - dc) / dc)
    end
    return worst
end

"""The sizes tried, 1.25 times apart, from one unit to `N`."""
function _ctsem_coreset_sizes(N::Integer)
    sizes = Int[]
    m = 1
    while m <= N
        push!(sizes, m)
        m = max(m + 1, round(Int, 1.25m))
    end
    return sizes
end

"""
    _ctsem_coreset_path(c, steps, budget, value, scale; size)

Forward selection of units for the training `steps`; the coreset that
`value(units, weights)` scores highest, if it scores above zero, as
`(units, weights, error)`, or `nothing`. Step `k` runs from training point `k`
to `k + 1`, and its rows are the change in the units' scores over it, scaled
per coordinate by `sqrt(scale)` -- L-BFGS's current inverse-curvature diagonal,
so that a coordinate counts as much as a step along it -- and per step by the
length of the full gradient's change over it in that metric, so that each
step's relative error counts alike; the target is their sum over units, which
every unit at weight one fits exactly. Each step adds the unit whose scores,
orthogonalised against those in, best explain what is left, among those that
help with a positive weight, while the units cost at most `budget` of a full
evaluation; the weights at each size tried (`_ctsem_coreset_sizes`) are
non-negative least squares (`_ctsem_nnls`). `value` returns the score, a bound
on any larger coreset's, and the error to keep; the search stops once the bound
is no better than the best. With `size`, the path goes to that many units and
is taken there, unscored.

With `information` (`_ctsem_coreset_information`), each parameter that counts
adds a row: the units' squared deviations of score from the mean over their
spread, which every unit at weight one fits exactly, so that the units chosen
carry the information as well as the changes in gradient. The rows together
weigh as one step's: the coreset's spread about the units' mean is linear in
the weights, and `_ctsem_coreset_shortfall`, which a coreset must pass, reads
it about the coreset's own.
"""
function _ctsem_coreset_path(c::CTSEMCoresetPlan, steps, budget::Real, value,
        scale::AbstractVector; size::Integer=0, information=nothing)
    root = sqrt.(max.(collect(Float64, scale), 0.0))
    N, p = Base.size(c.S[1])
    counted = information === nothing ? Int[] : findall(information.counts)
    rows = p * length(steps) + length(counted)
    A = Matrix{Float64}(undef, rows, N)
    b = Vector{Float64}(undef, rows)
    for (j, k) in enumerate(steps)
        nk = norm(root .* (c.gfull[k + 1] .- c.gfull[k]))
        nk = nk > 0 && isfinite(nk) ? nk : 1.0
        r0 = (j - 1) * p
        S0 = c.S[k]; S1 = c.S[k + 1]
        @inbounds for U in 1:N, i in 1:p
            A[r0 + i, U] = root[i] * (S1[U, i] - S0[U, i]) / nk
        end
    end
    if !isempty(counted)
        r0 = p * length(steps)
        scale_rows = 1 / sqrt(length(counted))
        SI = information.S
        @inbounds for (j, i) in enumerate(counted), U in 1:N
            A[r0 + j, U] = scale_rows * (SI[U, i] - information.mean[i])^2 /
                information.spread[i]
        end
    end
    sum!(reshape(b, rows, 1), A)
    Ap = copy(A)
    nrm = vec(sum(abs2, Ap; dims=1))
    floor0 = 1e-14 * maximum(nrm; init=0.0)
    r = copy(b)
    chosen = Int[]
    spent = 0.0
    along = zeros(N)
    t = zeros(N)
    best = nothing
    bestv = 0.0
    for target in (size > 0 ? [min(Int(size), N)] : _ctsem_coreset_sizes(N))
        while length(chosen) < target
            mul!(along, transpose(Ap), r)
            j = 0; top = 0.0
            @inbounds for U in 1:N
                a = along[U]
                if a > 0 && nrm[U] > floor0
                    s = a^2 / nrm[U]
                    s > top && (top = s; j = U)
                end
            end
            j == 0 && break
            spent + c.cost[j] > budget && break
            push!(chosen, j)
            spent += c.cost[j]
            q = Ap[:, j]
            q ./= norm(q)
            mul!(t, transpose(Ap), q)
            LinearAlgebra.BLAS.ger!(-1.0, q, t, Ap)
            # The lengths by downdating, and afresh every eighth unit, so that
            # rounding cannot leave a spent column looking long.
            if length(chosen) % 8 == 0
                sum!(abs2, reshape(nrm, 1, N), Ap)
            else
                @inbounds for U in 1:N
                    nrm[U] = max(nrm[U] - t[U]^2, 0.0)
                end
            end
            for U in chosen
                nrm[U] = 0.0
            end
            r .-= q .* dot(q, r)
        end
        isempty(chosen) && return nothing
        As = A[:, chosen]
        w = _ctsem_nnls(transpose(As) * As, transpose(As) * b)
        size > 0 && return (units=copy(chosen), weights=w, error=NaN)
        v, bound, err = value(chosen, w)
        if v > bestv
            best = (units=copy(chosen), weights=w, error=err)
            bestv = v
        end
        # No larger coreset can do better, or there is none.
        (bound <= bestv || length(chosen) < target) && break
    end
    return best
end

"""
    _ctsem_coreset_choose(c, scale, hmul; refit)

A coreset a stage of `c.stretch` iterations can pay for
(`_ctsem_coreset_budget`) whose changes in gradient over the steps between
training points, each held out once, are within `_CTSEM_CORESET_TOL` of the
full ones -- the mean squared relative error plus its standard error at most
the tolerance squared -- or `nothing` if none is; `c.rel` becomes that bound.

The steps are split into up to `_CTSEM_CORESET_FOLDS` contiguous folds, the
newest last -- for a stage carried on, ending with the stage just run. For
each, the path is fitted on the other steps and the first size within the
tolerance on the fold is that fold's; the coreset is fitted again on every
step and on the units' spread of scores at the newest training point, at the
first size from the largest of them up whose shortfall of that spread is within
the tolerance in every parameter (`_ctsem_coreset_shortfall`). Fitted in sample, a coreset carried on near the
optimum passed on an error a fraction of the one it met, its stage ran past
its floor, and the full objective came out worse.
"""
function _ctsem_coreset_choose(c::CTSEMCoresetPlan, scale::AbstractVector, hmul;
        refit::Bool)
    nsteps = length(c.S) - 1
    nsteps >= 2 || return nothing
    t0 = time()
    budget = _ctsem_coreset_budget(c.stretch)
    gg = [let dn = c.gfull[k + 1] .- c.gfull[k]; dot(dn, hmul(dn)) end for k in 1:nsteps]
    all(v -> isfinite(v) && v > 0, gg) || return nothing
    tol2 = _CTSEM_CORESET_TOL^2
    nfolds = min(nsteps, _CTSEM_CORESET_FOLDS)
    edges = round.(Int, range(0, nsteps; length=nfolds + 1))
    largest = 0
    rels = Float64[]
    for f in 1:nfolds
        held = (edges[f] + 1):edges[f + 1]
        D = [c.S[k + 1] .- c.S[k] for k in held]
        lik = [vec(sum(d; dims=1)) for d in D]
        gh = gg[held]
        # The first size within the tolerance on the fold: scored one there,
        # after which no larger can score more. The errors kept are the
        # fold's steps' own.
        value = (units, w) -> begin
            r = Vector{Float64}(undef, length(D))
            for j in eachindex(D)
                e = transpose(view(D[j], units, :)) * w .- lik[j]
                r[j] = dot(e, hmul(e)) / gh[j]
            end
            local mean_r = sum(r) / length(r)
            passes = isfinite(mean_r) && mean_r <= tol2
            (passes ? 1.0 : 0.0, passes ? 1.0 : 2.0, r)
        end
        found = _ctsem_coreset_path(c, setdiff(1:nsteps, held), budget, value, scale)
        if found === nothing
            empty!(rels)
            break
        end
        largest = max(largest, length(found.units))
        append!(rels, found.error)
    end
    result = nothing
    if length(rels) == nsteps
        m = sum(rels) / nsteps
        bound = m + sqrt(sum(abs2, rels .- m) / (nsteps - 1) / nsteps)
        if isfinite(bound) && bound <= tol2
            SK = c.S[end]
            info = _ctsem_coreset_information(SK, scale)
            tol = _CTSEM_CORESET_TOL
            need = largest
            spread = (units, w) -> begin
                ok = length(units) >= need &&
                    _ctsem_coreset_shortfall(SK, units, w, info) <= tol
                (ok ? 1.0 : 0.0, ok ? 1.0 : 2.0, nothing)
            end
            again = _ctsem_coreset_path(c, 1:nsteps, budget, spread, scale;
                information=info)
            again === nothing ||
                (result = (units=again.units, weights=again.weights, error=bound))
        end
    end
    c.attempts += 1
    c.select_seconds += time() - t0
    result === nothing && return nothing
    keep = findall(>(0), result.weights)
    isempty(keep) && return nothing
    c.rel = result.error
    units = result.units[keep]
    return CTSEMCoreset(units, result.weights[keep] .* c.sampling[units])
end

"""
    _ctsem_coreset_step!(c, x, f, G, hmul, scale)

After an accepted L-BFGS step to `x`, with value `f` and gradient `G` (of the
minimised objective in force). On the full objective: keep the point, and,
where the step gained more than `gap_tol`, look for a coreset
(`_ctsem_coreset_choose`) once the full iterations since the last search
number `wait`, which doubles after a search that finds none and is one after
one that finds one. In a
stage: check it with a full gradient (`_ctsem_coreset_check!`) once it has run
its planned iterations, once what L-BFGS's model says is left to gain on the
coreset (`G' H G / 2`) is no more than the coreset's error, once its gradient
is below `g_tol`, or once it has gone farther from where it began than the
coreset's training points reach, as a whole or in any one parameter
(`_ctsem_coreset_begin!`); and check it and leave it, carrying nothing on,
once a step is predicted to gain less than `gap_tol` (`stepgain`, as the gap
rule reads it), where the full sum's gap rule would stop the fit.
Neither the coreset's error nor `G' H G` is the full objective's inside a
stage; both are the stage's own measures of where it stands.

Returns `(act, x, f, G)`: `:none`; `:enter` or `:refit`, a stage that begins
at `x` on `c.target`, with the full value and gradient there; `:leave`, back to
the full objective at `x`; `:revert`, back to the full objective at the start
of the stage that just failed its check.
"""
function _ctsem_coreset_step!(c::CTSEMCoresetPlan, x, f, G, hmul, scale;
        stepgain::Real=Inf)
    none = (act=:none, x=x, f=f, G=G)
    c.seen += 1
    if !c.live
        gained = isfinite(c.f_last) ? c.f_last - f : Inf
        c.f_last = f
        _ctsem_coreset_record!(c, x, G) || return none
        gained > c.gap_tol || return none
        c.since += 1
        c.since >= c.wait || return none
        c.since = 0
        c.stretch = _CTSEM_CORESET_FIRST_STAGE
        cs = _ctsem_coreset_choose(c, scale, hmul; refit=false)
        if cs === nothing
            c.wait = max(1, 2c.wait)
            return none
        end
        c.wait = 1
        push!(c.entered, c.total_iterations)
        _ctsem_coreset_begin!(c, cs, x, f, G, scale, hmul)
        return (act=:enter, x=x, f=f, G=G)
    end
    c.iterations += 1
    c.total_iterations += 1
    left = dot(G, hmul(G)) / 2
    # A step predicted to gain less than the gap tolerance ends the stage, and
    # no stage is carried on from it: where the full sum's gap rule would stop
    # the fit, the full sum decides.
    quiet = stepgain < c.gap_tol
    due = quiet || c.iterations >= c.stretch || left <= c.error ||
        maximum(abs, G; init=0.0) <= c.g_tol ||
        norm((x .- c.x0) .* c.metric) > c.radius ||
        any(i -> abs(x[i] - c.x0[i]) > c.reach[i], eachindex(c.reach))
    due || return none
    return _ctsem_coreset_check!(c, x, f, G, hmul, scale; leave=quiet)
end

"""A stage on `cs` anchored at `x` (minimised value `f`, gradient `G`), trusted
within the training points' reach of `x` in the metric that divides each
coordinate by `sqrt(scale)` -- L-BFGS's inverse-curvature diagonal, so that a
unit is about one standard deviation of its quadratic model in every
coordinate -- and within their reach in each parameter, and whose error is the coreset's held-out relative error times what
L-BFGS's model (`hmul`) says is left to gain at `x`."""
function _ctsem_coreset_begin!(c::CTSEMCoresetPlan, cs::CTSEMCoreset, x, f, G,
        scale, hmul)
    c.live = true
    c.target = CTSEMCoresetTarget(c.full, cs, collect(Float64, x), zeros(length(x)))
    c.error = c.rel * dot(G, hmul(G)) / 2
    c.x0 = collect(Float64, x); c.f0 = f; c.G0 = collect(Float64, G)
    c.iterations = 0
    c.pairs = 0
    c.metric = [v > 0 && isfinite(v) ? 1 / sqrt(v) : 1.0 for v in scale]
    c.radius = maximum((norm((xk .- c.x0) .* c.metric) for xk in c.X); init=0.0)
    c.reach = [maximum((abs(xk[i] - c.x0[i]) for xk in c.X); init=0.0)
        for i in eachindex(c.x0)]
    push!(c.sizes, length(cs.units))
    push!(c.starts, c.seen)
    push!(c.shares, _ctsem_coreset_share(c, cs.units))
    push!(c.stage_iterations, 0)
    push!(c.bought, NaN)
    return nothing
end

"""
    _ctsem_coreset_check!(c, x, fc, Gc, hmul, scale; leave)

The full gradient at `x`, where the stage's coreset has value `fc`: back to the
stage's start if the full value is worse there or unusable; otherwise a
training point, and, unless `leave`, a stage carried on from `x`, planned twice
as long, if this one paid -- its iterations times the full objective's gain
over the coreset's predicted one, against its iterations' share of the work and
its check -- the full objective gained more than `gap_tol` an iteration over
it, and some coreset the longer stage can pay for is accurate enough; the full
objective otherwise.
"""
function _ctsem_coreset_check!(c::CTSEMCoresetPlan, x, fc, Gc, hmul, scale;
        leave::Bool=false)
    c.checks += 1
    n = c.iterations
    c.stage_iterations[end] = n
    Gf = similar(collect(Float64, x))
    ff = c.fg!(c.full, 0.0, Gf, x)
    ok = ff !== nothing && isfinite(ff) && ff < _CTSEM_BATCH_INVALID && all(isfinite, Gf)
    c.live = false
    c.since = 0
    # A usable check is a training point whatever it says.
    ok && _ctsem_coreset_record!(c, x, Gf)
    predicted = (c.f0 - c.offset) - fc
    actual = ok ? c.f0 - ff : -Inf
    fidelity = predicted > 0 ? actual / predicted : 0.0
    c.bought[end] = fidelity
    if !ok || ff > c.f0
        c.reverted += 1
        c.f_last = c.f0
        return (act=:revert, x=c.x0, f=c.f0, G=c.G0)
    end
    c.f_last = ff
    paid = n * fidelity > n * c.shares[end] + 1
    if !leave && paid && actual > c.gap_tol * max(n, 1)
        stretch = 2c.stretch
        c.stretch = stretch
        cs = _ctsem_coreset_choose(c, scale, hmul; refit=true)
        if cs !== nothing
            _ctsem_coreset_begin!(c, cs, x, ff, Gf, scale, hmul)
            return (act=:refit, x=x, f=ff, G=Gf)
        end
    end
    return (act=:leave, x=x, f=ff, G=Gf)
end

"""The objective in force after a coreset action `r` (`_ctsem_coreset_step!`):
for a stage that begins, the coreset's value at its point, read through
`evaluate!` (which reads the plan), and the full gradient, which the stage's
correction makes its own there; the full ones if the coreset cannot be
evaluated there, and otherwise. Returns the point, value, gradient and the
number of evaluations it took."""
function _ctsem_coreset_switch!(c::CTSEMCoresetPlan, r, evaluate!, batch)
    r.act === :none && return (r.x, r.f, r.G, 0)
    if r.act === :enter || r.act === :refit
        t = c.target
        fill!(t.correction, 0.0)
        Gc = similar(r.G)
        fc = evaluate!(batch, 0.0, Gc, r.x)
        if fc !== nothing && isfinite(fc) && fc < _CTSEM_BATCH_INVALID && all(isfinite, Gc)
            # Maximised: the full gradient less the coreset's, from minimised
            # arrays.
            t.correction .= Gc .- r.G
            c.offset = r.f - fc
            return (r.x, Float64(fc), copy(r.G), 1)
        end
        c.live = false
        return (r.x, r.f, r.G, 1)
    end
    return (r.x, r.f, r.G, 1)
end

"""A field of the plan's record for the bridge: `[0]` when there is no plan or
the field is empty."""
function _ctsem_coreset_report(c, field::Symbol)
    c === nothing && return field in (:shares, :bought) ? [0.0] : [0]
    v = getfield(c, field)
    return isempty(v) ? (eltype(v) <: AbstractFloat ? [0.0] : [0]) : copy(v)
end

"""The per-coordinate scale of L-BFGS's initial inverse Hessian: the learned
diagonal, the diagonal of a caller's inverse, or the metric's."""
_ctsem_lbfgs_scale(Dinv::AbstractVector, dinv) = Dinv
_ctsem_lbfgs_scale(Dinv::AbstractMatrix, dinv) = diag(Dinv)
_ctsem_lbfgs_scale(::Nothing, dinv) = dinv

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
# What a Hessian costs (`_ctsem_hessian_cost`) decides when one is formed, never
# a size:
#
# * the marginal route's is forward-over-adjoint, one lane per parameter at
#   about a gradient each: a few gradients on a small model (6 on a
#   1000-subject panel, dev1), about 750, a quarter of an hour, on a
#   715-parameter one. L-BFGS hands over early, on its own predicted gain, only
#   once the stage has spent at least that many gradients, since forming and
#   discarding the Hessian then at most doubles the fit; before that it runs
#   on to its own rule. Until 2026-10-17 the comparison was with the finish's
#   step budget, `maxit`, which made 30 parameters a boundary where the
#   hand-over changed for no change in cost;
# * the Laplace route's is central differences of the gradient, `2 npar`
#   gradients, and early hand-overs there made fits slower (0.6x with them,
#   1.3 to 1.9x without, dev1). So L-BFGS runs to its own stopping rule there.
#
# On either route the finish takes one Hessian where L-BFGS stopped and walks
# it (the chord), updated by SR1 after each step, and buys a fresh one only
# when the walk stops paying for it (the marginal-value rule in
# `_ctsem_newton_finish`, which weighs the steps' gains against
# `_ctsem_hessian_cost`); the final one is kept unless the steps moved the
# estimate by more than `_CTSEM_HESSIAN_REUSE_SE` standard errors. Where a
# Hessian costs a few gradients that rule refreshes it as readily as Newton
# would; where it costs hours it does not. Until 2026-10-17 the marginal route
# refreshed the exact Hessian at every step that contracted slowly, for up to
# `maxit` steps, which nothing bounded on a large model: on the bench's small
# marginal cells the chord matched it in every score and certification, at
# 629 gradients' cost against 647 (dev1, 2026-10-17). `newton = 'exact'` still
# asks for it.
#
# The steps are the model's minimiser within a trust region (`newton` inside
# `_ctsem_newton_finish`), so a Hessian formed at a point where its quadratic
# model holds only nearby -- negative curvature, or directions so weakly
# determined that their Newton step runs for raw units -- is still walked on,
# its flattest directions damped first.
#
# The early hand-over is on L-BFGS's own predicted gain, a proxy that after a
# few iterations knows little of the curvature, so where a Hessian is cheap the
# finish checks it: a first Newton step the objective does not take whole, or
# one that does not contract the predicted gain by `contraction`, is a point
# outside Newton's region, and the finish hands back to L-BFGS, which runs on
# to its own stopping rule before the finish runs again. Measured on the noise
# fixture of test-backend-summary.R from the prior warm-up's start: L-BFGS
# handed over at its third iteration on a proxy of 0.024 with an exact gain of
# 0.94 there, and the finish spent its cap of 30 steps and 13 Hessians walking
# 4.7 nats to a point it still called unfinished. Where a Hessian is dear the
# finish walks on it instead: until 2026-09-30 it handed back there too, and on
# a 715-parameter model threw away a quarter of an hour of Hessian twice, after
# which L-BFGS gained 0.02 log likelihood an iteration. The check also asked
# contraction of damped steps then, which a step damped to `alpha` cannot give
# below `(1 - alpha)^2` however right the model: on the synthetic analogue of
# that model (116 parameters, from its hand-over point), 3 negative directions
# and weak ones carried 77% of the predicted gain of 28.8, the line search took
# 1/16 of the step, and the check read a contraction of 1.63.
#
# Provenance: the hand-over at a predicted gain of 0.1, the cap of 30 steps and
# the contraction of 0.25 that triggers a refresh were set when this optimiser
# replaced Optim's (92499af5, 2026-09-24), on the seven regimes of
# `dev/stochopt/models.R` (panel, panel5k, long, ordinal Laplace, nonlin, small,
# bigp), and have not been swept since. The eigenvalue floor at 1e-8 of the
# largest, for the directions the certification does not trust, is the usual
# safeguard, not tuned; the trust region's constants are the textbook ones
# (Nocedal & Wright's Algorithm 4.1) and its acceptance the Armijo constant of
# every line search here. The consolidation plan's Appendix B has the list.

_ctsem_cheap_hessian(::CTSEMObjective) = true
_ctsem_cheap_hessian(::Any) = false

"""
    _ctsem_flat_escape(objective)

Whether the finish may leave a point along a direction its curvature calls flat
(`_ctsem_flat_ladder`); `_ctsem_flat_allowed` judges the point it would go to.
On the marginal objective, the likelihood itself, a gain along such a direction
is a gain. On the Laplace objective it is one only if the integral gains too:
a flat direction is where the approximation is least reliable, and on AnomAuth
S1 (bench, default start) the ladder raised the random drift's population sd
from raw -3.5 to 0.7, which the Laplace objective scored 25 nats up and the
quadrature correction 0.3 down (2026-10-18; decision 1 of
review/OPTIM-consolidation-plan-2026-09-25 recorded it as a scale collapsing).
That case kept the escape off the Laplace objective altogether until then, and
a count model's fit from a supplied start (test-julia-count.R) stopped where
the certification called it not stationary, 14 nats short along a population
sd that the correction agrees about. Not on the state-explicit objective,
which is sampled, never optimised.
"""
_ctsem_flat_escape(::Union{CTSEMObjective,CTSEMLaplaceObjective}) = true
_ctsem_flat_escape(::Any) = false

"""
    _ctsem_flat_allowed(objective, from, to, tol)

Whether the finish may take the flat ladder's point `to` from `from`: on the
Laplace objective, when the corrected objective gains more than `tol` there --
the continuation's value with its rules placed at each point
(`ctsem_laplace_continuation`, `centre_value`: quadrature for the units whose
Laplace term it finds off, the Laplace term for the rest and for units wider
than its `maxdim`). A point either value cannot be had at is refused.
Everywhere else, always.
"""
_ctsem_flat_allowed(::Any, from, to, tol) = true
function _ctsem_flat_allowed(laplace::CTSEMLaplaceObjective, from::AbstractVector,
        to::AbstractVector, tol::Real)
    corrected(p) = try
        local o = ctsem_laplace_continuation(laplace, collect(Float64, p))
        any(o.rule_failed) ? NaN : o.centre_value
    catch err
        _ctsem_must_propagate(err) && rethrow()
        NaN
    end
    a = corrected(from)
    isfinite(a) || return false
    b = corrected(to)
    return isfinite(b) && b - a > tol
end

"""
    _ctsem_hessian_cost(objective, npar)

What one Hessian costs, in gradients: `2 npar` on the Laplace route (central
differences of the gradient), and about `npar` on the marginal route, whose
forward-over-reverse sweeps carry one lane per parameter at about a gradient
each (measured: 750 gradients' time for 715 parameters). The unit the finish
weighs a Hessian against its steps in, each of which costs a gradient.
"""
_ctsem_hessian_cost(::CTSEMLaplaceObjective, npar::Integer) = 2 * Int(npar)
_ctsem_hessian_cost(::Any, npar::Integer) = Int(npar)

"""
    _ctsem_finish_curvature(objective)

The route's own curvature for the finish: `:chord` on the marginal and Laplace
routes, one Hessian walked with SR1 updates and refreshed when the walk stops
paying for it, whatever it costs; and `nothing` for a route the finish does not
run on: a pinned profile point, whose wrapper has no Hessian of its own, and
the state-explicit target, which is never certified.
"""
_ctsem_finish_curvature(::Union{CTSEMObjective,CTSEMLaplaceObjective}) = :chord
_ctsem_finish_curvature(::Any) = nothing

"""
How far the estimate may be from where a Hessian was evaluated, in the standard
errors that Hessian implies, for that Hessian still to give the standard errors
(R's `hessianreuse` control). A Wald standard error treats the curvature as
constant over about a standard error around the estimate, and carries an error
of however much it is not -- on the bench, 20 to 50 per cent over one standard
error on well identified models and 150 to 300 on weakly identified Laplace
ones (review/bench/2026-10-17-reuse-probe.csv). A Hessian `d` standard errors
away is off by about `d` of that, so a tenth adds a tenth to the error the
standard error already carries, whatever the model. R's
`.ctBackendHessianReuse()` is the same default and is what the fit passes in.
"""
const _CTSEM_HESSIAN_REUSE_SE = 0.1

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
    _ctsem_saddle_ladder(value_at, x, value, split, ascent; ladder, enough)

At a saddle, the best point along the directions of negative curvature: each
unit eigenvector, most negative first, both signs -- the gradient's first, since
along a saddle direction the gradient says which way is up -- at each rung of
`ladder` raw units. The Newton step cannot find this: it lives in the trusted
directions, where the point is already a maximum. `value_at` is as for
`_ctsem_flat_probe`. The walk stops after a direction that gained more than
`enough`, the most an escape is not taken for; `Inf` walks them all.

Returns the best improvement on `value` as `(point, value, length)`, or
`nothing`; the evaluations spent; and `walked`, the directions measured, which
is what lets a certification read a ladder that found nothing as a finding: a
direction counts once its shortest rung was evaluated on both sides, where a
curvature `c` would show as a gain of `c * ladder[1]^2 / 2` on at least one. A
route that refuses those points measured nothing along it. The sign of a
curvature that is zero by symmetry is rounding -- a population loading the data
fix the length of but not the direction among dimensions nothing else informs,
SNSF pilot: -0.016 against a largest eigenvalue of 1.2e6, under the bar by a
third, while its three twins came out positive -- and the ladder is what tells
that from a saddle.
"""
function _ctsem_saddle_ladder(value_at, x::AbstractVector, value::Real, split,
        ascent::AbstractVector; ladder=_CTSEM_SADDLE_LADDER, enough::Real=0.0)
    (split === nothing || !any(split.negative)) &&
        return (best=nothing, evaluations=0, walked=0)
    order = sort!(findall(split.negative); by=k -> split.values[k])
    best = nothing
    evaluations = 0
    walked = 0
    for k in order
        v = split.vectors[:, k]
        first = dot(ascent, v) >= 0 ? 1.0 : -1.0
        shortest = 0
        for side in (first, -first), len in ladder
            point = x .+ (side * len) .* v
            got = value_at(point)
            evaluations += 1
            isfinite(got) && len == ladder[1] && (shortest += 1)
            if isfinite(got) && got > value && (best === nothing || got > best.value)
                best = (point=point, value=got, length=side * len)
            end
        end
        shortest == 2 && (walked += 1)
        best !== nothing && best.value - value > enough && break
    end
    (best=best, evaluations=evaluations, walked=walked)
end

"""
    _ctsem_flat_ladder(fg!, x, f, G, split; lengths, corrections)

At a point with no negative curvature, the best point along the directions the
curvature calls flat, followed as a ridge: a step of each of `lengths` raw units
along the gradient's part outside the trusted directions, both signs, the
gradient's first -- along a flat direction its sign can be noise -- then up to
`corrections` Newton steps in the trusted directions, so the point rides the
ridge rather than leaving it. Each correction takes its direction from the
curvature at `x` and its length from the secant of the gradients at its two
ends: the Hessian at `x` describes the ridge less well the further along it the
rung is, and on a ridge curving as sharply as a parabola's four units out its
curvature was 2.5 times too strong, so uncorrected lengths each recovered only
two fifths of the way back. Returns the best point below `f` (minimised) as
`(point, f, G, length)`, or `nothing`, and the calls spent.

Why followed, not straight: a ridge flat in one direction is curved in raw
coordinates, so a straight step leaves it. On test-julia-convergence.R's random
walks fitted with an OU model, the fit can stop on the white-noise plateau,
where drift and diffusion trade along a ridge flat to 1e-5 over five raw units
and the Hessian is negative definite. Four raw units along it the straight
step has lost 1.04 nats and the followed one gained 0.02, and the maximum, 14.7
higher, lies further up the same ridge. The probe the certification reads
(`_ctsem_flat_probe`) steps straight, so it could not see this, and certified
the plateau.

With no gradient at all outside the trusted directions, the flattest direction.
`fg!` is the optimiser's trial path: an unusable point comes back as the
sentinel value, which no correction improves and no comparison accepts.
"""
function _ctsem_flat_ladder(fg!, x::AbstractVector, f::Real, G::AbstractVector,
        split; lengths=_CTSEM_FLAT_PROBE_LENGTHS, corrections::Integer=3)
    flat = findall(.!split.trusted)
    isempty(flat) && return (best=nothing, fcalls=0, gcalls=0)
    r = _ctsem_flat_residual(split, -G)
    if !(all(isfinite, r) && norm(r) > 0)
        r = split.vectors[:, flat[argmin(abs.(split.values[flat]))]]
    end
    u = r ./ norm(r)
    V = split.vectors[:, split.trusted]
    curvature = split.values[split.trusted]
    best = nothing; fcalls = 0; gcalls = 0
    for side in (1.0, -1.0), len in lengths
        y = x .+ (side * len) .* u
        g = similar(G)
        fy = fg!(0.0, g, y); fcalls += 1; gcalls += 1
        isfinite(fy) || continue
        for _ in 1:corrections
            d = -(V * ((transpose(V) * g) ./ curvature))
            slope = dot(d, g)
            slope < 0 || break
            z = y .+ d
            gz = similar(G)
            fz = fg!(0.0, gz, z); fcalls += 1; gcalls += 1
            # The secant of the slope along `d` between its ends: the step's
            # length where the slope vanishes, exact on a quadratic along it.
            ends = isfinite(fz) ? dot(d, gz) : NaN
            if isfinite(ends) && ends > slope
                t = slope / (slope - ends)
                if 0 < t <= 10 && abs(t - 1) > 0.05
                    zt = y .+ t .* d
                    gt = similar(G)
                    ft = fg!(0.0, gt, zt); fcalls += 1; gcalls += 1
                    if isfinite(ft) && !(isfinite(fz) && fz <= ft)
                        z = zt; fz = ft; gz = gt
                    end
                end
            end
            (isfinite(fz) && fz < fy) || break
            y = z; fy = fz; g = gz
        end
        if fy < f && (best === nothing || fy < best.f)
            best = (point=y, f=fy, G=g, length=side * len)
        end
    end
    (best=best, fcalls=fcalls, gcalls=gcalls)
end

"""
    _ctsem_sr1(H, s, y)

`H` updated by the symmetric rank-one formula so that it maps the step `s` to
the change in gradient `y` over it: `H + r r' / (r's)` with `r = y - H s`, which
gives `H s = y` and changes nothing orthogonal to `r`. What lets the chord's
copy of its Hessian follow a curvature that changes as it walks, from the
gradients the steps compute anyway. Unlike BFGS it does not keep `H` positive
definite, which the finish's model does not need: a direction the steps find
convex the other way is taken as the saddle-free step takes negative
curvature. Skipped -- `H` returned -- where `r's` is below `1e-8 |r| |s|`, the
usual safeguard (Nocedal & Wright, eq. 6.26), or not finite.
"""
function _ctsem_sr1(H::AbstractMatrix, s::AbstractVector, y::AbstractVector)
    r = y .- H * s
    denom = dot(r, s)
    (isfinite(denom) && abs(denom) > 1e-8 * norm(r) * norm(s)) || return H
    H .+ (r * transpose(r)) ./ denom
end

"""
    _ctsem_trust_multiplier(curvature, c, radius)

The Levenberg multiplier `mu >= 0` that puts the step `-c ./ (curvature .+ mu)`
-- a step in the eigenbasis of a model Hessian whose `curvature`s are all
positive, `c` the gradient there -- on the boundary of a trust region of
`radius`: zero when the Newton step (`mu = 0`) is inside it already. The step's
length falls monotonically in `mu`, and its reciprocal is nearly linear in it,
so Newton's method on the reciprocal converges from zero in a few iterations,
from below, never overshooting (Moré & Sorensen 1983). In the eigenbasis each
iteration is a sum, so this is the exact trust-region step at no cost beyond
the decomposition the finish already takes.
"""
function _ctsem_trust_multiplier(curvature::AbstractVector, c::AbstractVector,
        radius::Real)
    steplength(mu) = sqrt(sum(i -> (c[i] / (curvature[i] + mu))^2, eachindex(c);
        init=0.0))
    len = steplength(0.0)
    (isfinite(len) && len > radius) || return 0.0
    mu = 0.0
    for _ in 1:100
        len = steplength(mu)
        len <= radius * (1 + 1e-6) && break
        q = let mu = mu
            sum(i -> c[i]^2 / (curvature[i] + mu)^3, eachindex(c); init=0.0)
        end
        mu += (len^2 / q) * (len - radius) / radius
    end
    mu
end

"""
    _ctsem_newton_finish(objective, x, f, G, fg!; ...)

The endgame, from `x` on the minimised objective behind `fg!` (`f` and `G` its
value and gradient there): Newton steps within a trust region, the final
Hessian, an escape from a saddle, and the flat-direction probe.

Each step minimises the quadratic model within `radius` raw units (`newton`
inside, `_ctsem_trust_multiplier`): the Newton step while it fits, otherwise a
step with its flattest directions damped most. A Hessian's first step is its
Newton step; one the objective does not bear out is re-solved at a quarter of
its length, and the radius then follows how well the model predicted each step
(Nocedal & Wright's Algorithm 4.1).

What the steps are taken against is `curvature`:

- `:exact`  the exact Hessian, refreshed whenever a Newton step does not
  contract the predicted gain by `contraction`, or a step the trust region held
  back gained less than a quarter of what the model predicted;
- `:chord`  the exact Hessian at `x`, kept for every step (the chord, or
  simplified Newton, method: linear convergence at the rate the Hessian
  changes between `x` and the optimum, which from a hand-over this close is
  fast), with the steps' copy of it updated after every step by the step and
  its change in gradient (`_ctsem_sr1`), until the steps are not worth a
  Hessian: at the rate the last half of the walk gained,
  closing the gain its model predicts over the trusted directions would take
  more steps than a Hessian costs in gradients, and the point has moved more
  than `reuse_se` standard errors from where the Hessian was taken;
- `:subset` the likelihood Hessian of a random `subset` share of the units
  (at least `subset_min`), scaled up to the data, with the prior's curvature
  kept whole -- a chord Hessian at a fraction of the cost.

The final Hessian is the exact one at the final point, with one exception: a
chord Hessian whose steps converged within `reuse_se` standard errors of where
it was evaluated is kept (`hessian_at` and `distance` say where, and how far
that is -- `_ctsem_hessian_distance`). Whatever the steps used, a step the
trust region cannot make on a stale matrix is answered with the exact Hessian.
A subset Hessian, or a chord one whose steps converged farther away, is
replaced by the exact Hessian at the final point and the chord walked again on
that, as is a chord that ended without converging, while the finish has spent
less than six Hessians and walks of as many steps as one costs (counted in
gradients, `_ctsem_hessian_cost`); on the exact variant the steps continue on a
fresh Hessian each until it agrees.

At a saddle -- negative curvature in the final Hessian -- the ladder along the
negative curvature is tried (`_ctsem_saddle_ladder`), and the finish
continues from the better point with the exact Hessian there, at most
`max_escapes` times: it never returns a point with a direction of negative
curvature it did not try, short of that cap. What the last ladder found at the
returned point is returned with it (`saddle_walked`, `saddle_gain`), so a
direction whose negative curvature is rounding, walked and found flat, is a
measurement the certification reads rather than a saddle. With no negative curvature but a
direction the curvature calls flat, or a saddle whose ladder found nothing
worth having, the flat ladder (`_ctsem_flat_ladder`) is tried, when
`flat_escape` (the fit sets it where `_ctsem_flat_escape` allows), along what
`flat_allowed(from, to)` accepts its point (`_ctsem_flat_allowed`), within the same cap, the steps continuing on the exact Hessian there with a fresh step
budget. Either escape is taken only when it gains more than `certify_tol`.
Then, when `probe`, the flat-direction probe (`_ctsem_flat_probe`) at the final
point.

The two tolerances divide by what a decision reads. `certify_tol`, the
certification's own (the fit sets it; 0 turns off what follows), is for what
is exact: a gain read off the exact Hessian at the point ends the steps once
it is below it, and an escape has to gain more. `tol`, the fit's inner
tolerance (a hundredth of the other by default), is for what is not: a gain
predicted by a matrix taken elsewhere, the chord's, or a refreshed one's
first step, has to fall below it, as L-BFGS's own prediction does, since such
a prediction can be optimistic.

`hessian0`, the Hessian of the MAXIMISED objective at `hessian0_at`, is one the
caller holds -- the certification's, for a stage it resumed -- and is walked on
in place of forming one at `x`, on the chord and the subset; the exact variant
forms its own, as it is cheap.

`take_steps = false` gives the certification's numbers alone: the Hessian at `x`,
the ladder's measurement along any negative curvature there and the probe, with
no step and no refresh, for a point the optimiser left without a finish
(`ctsem_endgame`).

`handback = true` is for a hand-over L-BFGS made on its own predicted gain
rather than at its stopping rule, where the Hessian is cheap enough to discard
(`ctsem_optimize` decides): when the objective does not take the first Newton
step whole, or that step does not contract the predicted gain by
`contraction`, the point is outside Newton's region and the finish returns
there, with `handback` set and no final Hessian, for L-BFGS to go on from.

`value_at(y)`, the maximised objective or `-Inf` for the ladder and the probe,
defaults to the route's own predicate (`_ctsem_probe_value`); `fg!` is the
optimiser's trial path, for the steps.

Returns the point and its value and gradient (minimised); the Hessian of the
MAXIMISED objective, or `nothing`, with `hessian_at` and `distance`; the steps
taken, including escapes, and the predicted gain at the end; the full and
subset Hessians formed and the calls made; the saddle record (`escapes`,
`saddle`, `ladder_tried`, `ladder_gain`, and at the returned point
`saddle_walked` and `saddle_gain`, 0 and `NaN` when no ladder ran there); the
probe record; whether it handed
back (`handback`); and the history of the steps, one entry each (`kind` is "newton", "exact", "saddle" or "flat", and `alpha`
the share of its Newton step's length the step took).
"""
function _ctsem_newton_finish(objective, x0, f0, G0, fg!; tol::Real=1e-8,
        maxit::Integer=30, contraction::Real=0.25, callback=nothing,
        iteration0::Integer=0, curvature::Symbol=:exact,
        subset::Real=0.125, subset_min::Integer=200,
        rng=Random.MersenneTwister(20260924), take_steps::Bool=true,
        probe::Bool=true, reuse_se::Real=_CTSEM_HESSIAN_REUSE_SE,
        flat_rtol::Real=_CTSEM_FLAT_RTOL, negative::Real=_CTSEM_NEGATIVE_RTOL,
        ladder=_CTSEM_SADDLE_LADDER, probe_lengths=_CTSEM_FLAT_PROBE_LENGTHS,
        max_escapes::Integer=3, value_at=nothing, handback::Bool=false,
        certify_tol::Real=0.0, flat_escape::Bool=false,
        flat_allowed=(from, to) -> true,
        reporter=nothing, hessian0=nothing, hessian0_at=nothing)
    curvature in (:exact, :chord, :subset) || throw(ArgumentError(
        "newton curvature must be exact, chord or subset, got $(curvature)"))
    valueof = value_at === nothing ? (y -> _ctsem_probe_value(objective, y)) :
        value_at
    x = collect(Float64, x0); f = Float64(f0); G = collect(Float64, G0)
    # The closures below count into Refs and are handed the current point, as
    # `_ctsem_lbfgs`'s are: a captured local assigned more than once is a
    # `Core.Box`, and an untyped point sent inference through the objective
    # for `Any`.
    full_hessians = Ref(0); subset_hessians = Ref(0)
    # Reports a Hessian this finish forms as it goes, rather than leaving the
    # caller's progress line frozen for as long as forming one takes -- up to
    # 80s on the Laplace fixture that found this, and 15 minutes on a
    # 715-parameter marginal model, all of it inside one `ctsem_hessian` call
    # the line above had no way to see into. The Laplace method reports per
    # gradient of its finite differences, the marginal one per forward sweep.
    # `nothing` when there is nothing to report through (see
    # `_ctsem_hessian_progress`), which costs nothing: each method calls back
    # only when it is not `nothing`.
    hessian_progress = _ctsem_hessian_progress(reporter)
    hessof(o, y) = try
        local Hy = Matrix{Float64}(ctsem_hessian(o, y; progress=hessian_progress))
        all(isfinite, Hy) ? Hy : nothing
    catch err
        # A point the model cannot evaluate is no data; code that is wrong is
        # an error, and an interrupt stops the fit (`_ctsem_must_propagate`).
        _ctsem_must_propagate(err) && rethrow()
        nothing
    end
    # Of the minimised objective, as every step below expects.
    # `local`: an assignment to `H` in a closure would rebind the H below.
    hess(y) = (full_hessians[] += 1; local Hx = hessof(objective, y);
        Hx === nothing ? nothing : -Hx)
    NU = _ctsem_nunits(objective)
    m = min(NU, max(Int(subset_min), ceil(Int, subset * NU)))
    sub = (curvature === :subset && 2m <= NU) ?
        _ctsem_subset_objective(objective, sort(randperm(rng, NU)[1:m])) : nothing
    function subhess(y)
        sub === nothing && return hess(y)
        subset_hessians[] += 1
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
    remember!(kind, g, a, f) = (push!(history.kind, kind); push!(history.gain, g);
        push!(history.alpha, a); push!(history.value, -f); nothing)
    fcalls = Ref(0); gcalls = 0
    # The model the steps are taken on: the eigen-decomposition of the matrix
    # they use, with the curvatures described below. Kept until that matrix
    # changes, since a step the trust region refuses is re-solved on it.
    modelled = Ref{Any}(nothing)
    function model(Hm)
        local M = modelled[]
        (M !== nothing && M.matrix === Hm) && return M
        local E = eigen(Symmetric(Hm))
        local lmax = maximum(abs, E.values; init=0.0)
        local top = maximum(E.values)
        M = (matrix=Hm, vectors=E.vectors, usable=lmax > 0,
            curvature=[top > 0 && v > flat_rtol * top ? v :
                max(abs(v), 1e-8 * lmax) for v in E.values],
            trusted=[top > 0 && v > flat_rtol * top for v in E.values])
        modelled[] = M
        M
    end
    # The gain the undamped step predicts, flat directions included at their
    # floored curvature. Judging convergence over the trusted directions alone
    # stopped the finish at the start of a nearly flat ray: on a fixture with a
    # diffusion correlation flat below raw -6, 3.5e-6 nats short of the point a
    # profile then found, so the estimate was not the maximum. A truly flat
    # direction has no gradient and adds nothing; a nearly flat one is walked
    # while the step still promises more than `tol`, as L-BFGS used to.
    #
    # A direction of negative curvature is taken at the size of its curvature,
    # not at the floor (the saddle-free Newton step): downhill along it by a
    # length the curvature can justify. Floored, its component was up to 1e8
    # times too long, the line search shrank the whole step to tame it, and the
    # trusted directions' part of the step -- the part that closes the gap --
    # shrank with it; the chord never refreshes, so every later step did the
    # same. On gated-gaps config B8 (bench, seed 1) seven chord steps from a
    # point 10.8 nats short, with one direction of negative curvature, moved
    # the log posterior by less than 1e-4: inferred to be this, from those
    # numbers. At a saddle proper the gradient along that direction is zero
    # either way, and the ladder below is what leaves it.
    #
    # And a direction the certification trusts -- curvature above `flat_rtol`
    # of the largest, `_ctsem_information_split`'s rule -- is taken at its own
    # curvature however small, so the finish closes the gap the certification
    # will measure; only the rest are floored. Floored as well, a trusted
    # direction below 1e-8 of the largest was walked at a fraction of its
    # Newton step, gaining a few percent of the missing likelihood a step: on
    # AnomAuth S1 (bench, default start) a random-effect sd on a ray toward
    # zero, at relative curvature 9e-11 and a Newton step of 0.25 raw units,
    # took 35 chord steps of 1e-7 nats each against a gap of 1.1e-5, and the
    # certification resumed the fit twice for it: 21 Hessians.
    #
    # The step itself is the minimiser of that model within a trust region of
    # `radius` raw units (`_ctsem_trust_multiplier`): the Newton step while it
    # fits, and otherwise the step whose Levenberg multiplier `mu` puts it on
    # the boundary, which damps each direction by `curvature / (curvature + mu)`
    # -- the flattest most, the stiff ones hardly at all. `gain` is the undamped
    # step's, the convergence test's; `predicted` is what the model promises
    # for the step taken, which the trust region judges it by; `known` is
    # `gain` over the trusted directions alone, what the model can say is left
    # -- a floored direction can make `gain` large on a model that cannot say
    # where along it the maximum is.
    function newton(Hm, Gv, radius)
        local M = model(Hm)
        M.usable || return nothing
        local lam = M.curvature
        local c = transpose(M.vectors) * Gv
        local mu = _ctsem_trust_multiplier(lam, c, radius)
        local d = c ./ (lam .+ mu)
        local known = 0.0
        for i in eachindex(c)
            M.trusted[i] && (known += abs2(c[i]) / lam[i])
        end
        (step=-(M.vectors * d), gain=0.5 * sum(abs2.(c) ./ lam; init=0.0),
            known=0.5 * known,
            predicted=sum(abs2.(c) .* (lam .+ 2mu) ./ (2 .* (lam .+ mu) .^ 2); init=0.0),
            length=norm(d), whole=mu == 0,
            fraction=mu == 0 ? 1.0 : norm(d) / norm(c ./ lam))
    end
    # A step on `Hm` from the current point, starting from `nt` (its step at
    # `radius`): accepted when the objective falls by at least 1e-4 of what the
    # model predicted -- the Armijo constant of every line search here -- and
    # by more than it can represent, and otherwise re-solved at a quarter of
    # its length (Nocedal & Wright's Algorithm 4.1), until the model promises
    # less than that resolution. Armijo alone scales with the step, so a small
    # enough step would take an increase of 1e-14, spending a step on a point
    # no different from the one it left. The names inside are `local` so they
    # cannot rebind the finish's own.
    function trust(Hm, nt, radius, x, f, G)
        local floor = max(abs(f), 1.0) * eps()
        local at = nt
        local r = radius
        local k = 0
        while true
            local xn = x .+ at.step
            local fn = fg!(0.0, nothing, xn)
            fcalls[] += 1
            local rho = isfinite(fn) ? (f - fn) / at.predicted : -Inf
            if isfinite(fn) && rho > 1e-4 && f - fn > floor
                return (ok=true, x=xn, f=fn, rho=rho, nt=at, radius=r)
            end
            (at.predicted > floor && k < 100) ||
                return (ok=false, x=x, f=f, rho=rho, nt=at, radius=r)
            k += 1
            r = at.length / 4
            at = newton(Hm, G, r)
        end
    end
    report(g, steps, f, G, x) = callback === nothing || callback(CTSEMIterate(iteration0 + steps,
        f, maximum(abs, G; init=0.0), g, x))
    # A Hessian the caller holds, where Hessians are dear: `hessian0` (of the
    # MAXIMISED objective, as the finish returns it) at `hessian0_at`, walked
    # on rather than formed again. A stage resumed after a certification
    # starts where the certification's Hessian was taken and hands over within
    # a few L-BFGS steps: on a 529-parameter marginal fit (ECLS-K, dev2) its
    # finish formed a new one there, 47 minutes, every resume. The walk judges
    # it as it judges its own (below), and the final Hessian is still the
    # exact one at the final point unless the steps ended within `reuse_se`
    # standard errors of `hessian0_at`.
    seeded = curvature !== :exact && hessian0 !== nothing &&
        hessian0_at !== nothing && size(hessian0) == (length(x), length(x)) &&
        length(hessian0_at) == length(x) && all(isfinite, hessian0) &&
        all(isfinite, hessian0_at)
    H = seeded ? -Matrix{Float64}(hessian0) : step_hessian(x)
    # What the steps are taken against: `H` itself, except that the chord
    # corrects its copy along the direction it walks (see the refresh below),
    # while `H`, exact where it was evaluated, is what the reuse test and the
    # certification read.
    Hs = H
    hat = seeded ? collect(Float64, hessian0_at) : copy(x)   # where `H` was evaluated
    exact = seeded || curvature !== :subset   # `H` is the exact Hessian at `hat`
    steps = 0; escapes = 0; gain = Inf
    saddle = false; ladder_tried = false; ladder_gain = 0.0
    saddle_walked = 0; saddle_gain = NaN
    firstgain = Inf; firstwhole = true; handed = false
    if H === nothing
        return (x=x, f=f, G=G, hessian=nothing, hessian_at=hat, distance=NaN,
            steps=0, gain=Inf, full_hessians=full_hessians[],
            subset_hessians=subset_hessians[], fcalls=0, gcalls=0, escapes=0,
            saddle=false, ladder_tried=false, ladder_gain=0.0, saddle_walked=0,
            saddle_gain=NaN, probe=_ctsem_no_probe(), history=history,
            handback=false)
    end
    at_x = exact && hat == x            # `H` is the exact Hessian at `x`
    # The trust region, in raw units. Each Hessian starts without one -- its
    # first step is the Newton step -- and the region is set by what the steps
    # on it then show.
    radius = Inf
    # The steps a Hessian is walked for: `maxit` on the exact variant. On the
    # chord and the subset a walk goes on while its steps are worth what a
    # Hessian costs (the rule below), and the finish as a whole may spend, in
    # gradients, what six Hessians and walks of as many steps as one costs
    # did under the count this replaced (`finish_budget`): a step costs about
    # a gradient and a Hessian `hessian_cost` of them (`_ctsem_hessian_cost`).
    # Counted as a number of rewalks, walks cut short for a fresh Hessian used
    # the allowance up sooner, and a finish that could not close its gap ended
    # having covered less ground. A walk that ends without converging, or
    # converges somewhere its Hessian no longer describes, is walked again on
    # the exact Hessian there (below) while the budget lasts. An escape's new
    # Hessian is walked for `per_hessian` steps whatever is left.
    hessian_cost = _ctsem_hessian_cost(objective, length(x))
    subset_cost = sub === nothing ? Float64(hessian_cost) : hessian_cost * m / NU
    per_hessian = curvature === :exact ? Int(maxit) : max(Int(maxit), hessian_cost)
    finish_budget = 6 * (hessian_cost + per_hessian)
    budget = curvature === :exact ? per_hessian :
        max(0, floor(Int, finish_budget - hessian_cost * full_hessians[] -
            subset_cost * subset_hessians[]))
    while true                          # once, and again after an escape or a rewalk
        converged = false
        prevgain = Inf; prevwhole = true; prevrho = 1.0
        # This walk's own record: where it started, the objective after each
        # of its steps, and the curvature its Hessian implies, for `worth`.
        walk_start = steps; walk_values = [f]
        walk_split = curvature === :exact ? nothing :
            _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
        while take_steps && steps < budget
            nt = newton(Hs, G, radius)
            nt === nothing && break
            gain = nt.gain
            # A gain read off the exact Hessian at the point (`at_x`) is the
            # number the certification computes, so within its tolerance the
            # point is certified as it stands; `tol`'s margin is for
            # predictions that can be optimistic, which this is not. Stepping
            # on toward `tol` moved the point off its Hessian, and past
            # `reuse_se` standard errors another was needed to certify where
            # it went: on ECLS-K at full N (529 parameters, 3 h a Hessian, a
            # tolerance of 2) the finish re-formed its Hessian after 5 steps,
            # then 1, then 1, each time to close part of what was already
            # within the tolerance.
            if at_x && gain < certify_tol
                converged = true
                break
            end
            if gain < tol
                # On the exact variant a gap is closed when the Hessian here
                # says so. One from an earlier point can call it closed on a
                # plateau that is still rising, and the steps then stopped:
                # test-julia-particle.R's linear fit, from a saddle, reached a
                # plateau where drift runs to minus infinity, converged there
                # on a two-step-old Hessian, and was certified at -39.879 on
                # 6.4e-7 of predicted gain after a capped check of five steps
                # whose gains grew. The exact profile over drift rises from
                # there to -25.189 at the boundary, and steps on a fresh
                # Hessian climb it. So the steps go on under the budget.
                if curvature === :exact && !at_x
                    H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                    H === nothing && break
                    prevgain = Inf; prevwhole = true; prevrho = 1.0; radius = Inf
                    continue
                end
                converged = true
                break
            end
            # Where Hessians are dear, a walk ends for a fresh Hessian (the
            # rewalk below) once its steps are not worth one, judged at the rate
            # the last half of the walk gained, per step and so per gradient:
            #
            # - that rate is below the walk's average, its Hessian's cost
            #   included -- the point of diminishing returns, past which the
            #   steps earn less than a Hessian and its walk did (the marginal
            #   value rule). A fresh walk whose first steps are slow has not
            #   yet earned its Hessian, so it is not judged on them: without
            #   this a refresh followed every refresh two steps later on a
            #   curved valley, and the rewalks ran out short of the optimum;
            # - and closing what its model says is left (`known`, over the
            #   trusted directions) would take more steps than a Hessian costs.
            #   A walk still converging never meets this -- what is left
            #   shrinks with its steps' gains -- and one whose Hessian no
            #   longer describes the valley it walks meets it soon;
            # - and the point is more than `reuse_se` standard errors from where
            #   the Hessian was taken, or the fresh one would be the same.
            #
            # Each walk of a 529-parameter marginal fit (ECLS-K, dev2) gained
            # nine tenths of its total in its first 25 steps and then walked on
            # for 500 more, as many as its Hessian cost, a few hundredths of a
            # nat a step: seven walks and a resume in 19 hours.
            k = steps - walk_start
            if curvature !== :exact && k >= 2
                half = k ÷ 2
                recent = (walk_values[half + 1] - f) / (k - half)
                if recent * (hessian_cost + k) < walk_values[1] - f &&
                        nt.known > hessian_cost * recent &&
                        _ctsem_hessian_distance(walk_split, hat, x) > reuse_se
                    break
                end
            end
            # The check `handback` asks for, once: in Newton's region the
            # objective takes the whole Newton step, and the step contracts the
            # predicted gain -- the hand-over Hessian's prediction at the point
            # it reached -- by `contraction`. A first step it did not take whole
            # says the point is outside that region by itself, and contraction
            # is asked only of a whole one: a step damped to `alpha` contracts
            # the gain to no less than `(1 - alpha)^2` of itself however right
            # the model, so asking it of a damped step, as this did until
            # 2026-09-30, found every damped first step outside the region.
            if handback && steps == 1 && escapes == 0 &&
                    (!firstwhole || gain > contraction * firstgain)
                handed = true
                break
            end
            steps == 0 && (firstgain = gain)
            trial = trust(Hs, nt, radius, x, f, G)
            if !trial.ok
                # The model promises nothing the objective can represent at
                # any radius: a stale matrix is replaced, and on the exact
                # one here the steps are over.
                at_x && break
                H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                H === nothing && break
                prevgain = Inf; prevwhole = true; prevrho = 1.0; radius = Inf
                continue
            end
            Gn = similar(G)
            fg!(nothing, Gn, trial.x); gcalls += 1
            # A refused gradient at the step's end (see `trial_fg!`): the
            # finish stops where it was rather than step to a point it cannot
            # differentiate.
            all(isfinite, Gn) || break
            dx = trial.x .- x; dg = Gn .- G
            x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
            push!(walk_values, f)
            steps == 1 && (firstwhole = trial.nt.whole)
            remember!("newton", gain, trial.nt.fraction, f)
            report(gain, steps, f, G, x)
            # The radius follows how well the model predicted the step's gain:
            # a quarter of the step when it predicted less than a quarter of
            # it, twice the radius when it predicted three quarters and the
            # region held the step back (Nocedal & Wright's Algorithm 4.1).
            radius = trial.rho < 0.25 ? trial.nt.length / 4 :
                (trial.rho > 0.75 && !trial.nt.whole) ? 2 * trial.radius :
                trial.radius
            # Only the exact variant refreshes on slow contraction; the chord
            # and the subset keep their matrix, which is the point of them. A
            # step that fails outright still gets the exact Hessian (above).
            #
            # Also after a step that was on a fresh one, where slow contraction
            # means the quadratic model does not hold rather than a stale
            # matrix. Stopping the refreshes there was measured and it is worse:
            # on the noise fixture of test-backend-summary.R from twelve random
            # starts (cores = 1), 68 Hessians in all against 222, and four fits
            # left on the ridge at -210.13, three of them notstationary, that
            # the refreshed steps had climbed off -- to -207.02 twice, -207.29
            # and -205.83.
            #
            # The chord's refresh is a free one: after every step its copy of
            # the matrix is updated by the step and the change in gradient over
            # it, which the steps compute anyway (`_ctsem_sr1`), so the curvature
            # it walks on follows the walk rather than staying where the
            # Hessian was taken. Without any update the chord walked a ray whose
            # curvature decays along it -- a random-effect sd heading toward
            # zero -- at the curvature of where it started, every step shorter
            # than the last: on AnomAuth S1 from its default start (bench,
            # dev1), 30 chord steps for gains of 4e-5 down to 3e-7, then five
            # exact Hessians to finish, seven in the stage. A correction along
            # the step alone, on slow steps, mended that ray and left the rest
            # of the matrix fixed: a curved valley bends it in other directions
            # too, and on a 529-parameter marginal fit (ECLS-K, dev2) each walk
            # gained nine tenths of its total in 25 steps and crawled on for
            # 500 more.
            #
            # A step the trust region held back cannot contract the gain as a
            # Newton step does -- the directions it damped keep their share of
            # it -- so contraction is asked of whole steps only, and a held
            # step counts as slow when the model predicted its gain poorly, the
            # trust region's own measure of a model that no longer holds.
            # Asking contraction of a held step too refreshed an accurate
            # Hessian after every step the region shortened.
            slow = steps > 1 &&
                (prevwhole ? gain / prevgain > contraction : prevrho < 0.25)
            prevwhole = trial.nt.whole; prevrho = trial.rho
            if curvature === :exact && slow
                H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                H === nothing && break
                prevgain = Inf; prevwhole = true; prevrho = 1.0; radius = Inf
            else
                curvature !== :exact && (Hs = _ctsem_sr1(Hs, dx, dg))
                prevgain = gain
            end
        end
        # Handed back: the point the step reached, and no final Hessian, since
        # L-BFGS goes on from here and the finish runs again where it stops.
        handed && return (x=x, f=f, G=G, hessian=nothing, hessian_at=hat,
            distance=NaN, steps=steps, gain=gain, full_hessians=full_hessians[],
            subset_hessians=subset_hessians[], fcalls=fcalls[], gcalls=gcalls,
            escapes=0, saddle=false, ladder_tried=false, ladder_gain=0.0,
            saddle_walked=0, saddle_gain=NaN, probe=_ctsem_no_probe(),
            history=history, handback=true)
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
        #
        # Where Hessians are dear (the chord and the subset), the new one is
        # walked on as the first was -- a rewalk -- rather than replaced after
        # every step: each step there used to cost a Hessian of its own, five
        # of them in the worst case, and on the synthetic analogue of a
        # 715-parameter model those were most of a stage's time. So is a chord
        # that used up its steps without converging, which is walking on a
        # Hessian that no longer describes where it is: on that analogue a
        # stale chord predicted a sixtieth of the gap the certification then
        # found, and its steps each gained about that.
        if !at_x
            keep = curvature !== :exact && exact && converged &&
                _ctsem_hessian_distance(_ctsem_information_split(H;
                    rtol=flat_rtol, negative=negative), hat, x) <= reuse_se
            if !keep && curvature !== :exact
                H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                H === nothing && break
                radius = Inf
                left = finish_budget - hessian_cost * full_hessians[] -
                    subset_cost * subset_hessians[] - steps
                if take_steps && left >= 1
                    budget = steps + floor(Int, left)
                    continue
                end
            elseif !keep
                for _ in 1:5
                    if !at_x
                        H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                        radius = Inf
                    end
                    H === nothing && break
                    take_steps || break
                    nt = newton(H, G, radius)
                    gain = nt === nothing ? Inf : nt.gain
                    (nt === nothing || gain < tol || steps >= maxit + 5) && break
                    trial = trust(H, nt, radius, x, f, G)
                    trial.ok || break
                    Gn = similar(G)
                    fg!(nothing, Gn, trial.x); gcalls += 1
                    all(isfinite, Gn) || break
                    x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
                    remember!("exact", gain, trial.nt.fraction, f)
                    report(gain, steps, f, G, x)
                end
                # A step on the last round leaves the Hessian one point behind.
                if !at_x && H !== nothing
                    H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
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
        split = _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
        split === nothing && break
        # The certification's numbers alone: what the negative curvature is
        # worth, measured along every such direction, and no step.
        if !take_steps
            if any(split.negative)
                saddle = true; ladder_tried = true
                tried = _ctsem_saddle_ladder(valueof, x, -f, split, -G;
                    ladder=ladder, enough=Inf)
                fcalls[] += tried.evaluations
                saddle_walked = tried.walked
                saddle_gain = tried.best === nothing ? 0.0 : tried.best.value + f
            end
            break
        end
        # An escape has to gain more than `certify_tol`, the certification's
        # tolerance: less is a move the verdict cannot tell from staying, and it
        # spends one of the `max_escapes`. On test-julia-convergence.R's random
        # walks the Hessian on the white-noise plateau showed a noise-level
        # negative curvature, and three saddle escapes of 0.03 raw units, gaining
        # nothing at eight digits, used the cap before anything else was tried.
        if any(split.negative)
            saddle = true
            escapes >= max_escapes && break
            ladder_tried = true
            tried = _ctsem_saddle_ladder(valueof, x, -f, split, -G; ladder=ladder,
                enough=certify_tol)
            fcalls[] += tried.evaluations
            saddle_walked = tried.walked
            saddle_gain = tried.best === nothing ? 0.0 : tried.best.value + f
            if tried.best !== nothing
                Gn = similar(G)
                fn = fg!(0.0, Gn, tried.best.point); fcalls[] += 1; gcalls += 1
                # The point the ladder measured, through the optimiser's own
                # trial path: a point the route refuses a gradient at is not one
                # to continue from.
                if isfinite(fn) && f - fn > certify_tol
                    ladder_gain += f - fn
                    x = collect(Float64, tried.best.point); f = fn; G = Gn
                    steps += 1; escapes += 1
                    remember!("saddle", NaN, tried.best.length, f)
                    report(NaN, steps, f, G, x)
                    H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                    saddle = false; ladder_tried = false
                    saddle_walked = 0; saddle_gain = NaN
                    H === nothing && break
                    # A new point and a new Hessian: walked as the first one
                    # was, with its own budget and no trust region. Without the
                    # budget an escape after the steps had spent theirs was
                    # returned unwalked: a 529-parameter marginal fit (ECLS-K,
                    # dev2, 2026-10-16) ended its finish on the ladder's point
                    # with |g| 1.3e3 against 7 before it, its certification
                    # found 2 nats left, and the resume began again with
                    # L-BFGS and the Hessian it had just formed.
                    radius = Inf
                    budget = max(budget, steps + per_hessian)
                    continue
                end
            end
        end
        # No saddle to leave, or none worth leaving: a direction the curvature
        # does not trust may still have somewhere better along it, which neither
        # the Newton step (no gradient there) nor the certification's straight
        # probe can see -- `_ctsem_flat_ladder`, where the route allows it.
        (flat_escape && escapes < max_escapes) || break
        walked = _ctsem_flat_ladder(fg!, x, f, G, split; lengths=probe_lengths)
        fcalls[] += walked.fcalls; gcalls += walked.gcalls
        (walked.best !== nothing && f - walked.best.f > certify_tol &&
            flat_allowed(x, walked.best.point)) || break
        x = collect(Float64, walked.best.point); f = walked.best.f
        G = walked.best.G
        steps += 1; escapes += 1
        remember!("flat", NaN, walked.best.length, f)
        report(NaN, steps, f, G, x)
        H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
        saddle = false; ladder_tried = false
        saddle_walked = 0; saddle_gain = NaN
        H === nothing && break
        # A new point and a new Hessian: walked as the first one was.
        radius = Inf
        budget = max(budget, steps + per_hessian)
    end
    split = H === nothing ? nothing :
        _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
    if H !== nothing
        nt = newton(H, G, Inf)
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
            fcalls[] += probed.evaluations
        end
    end
    (x=x, f=f, G=G, hessian=H === nothing ? nothing : -H, hessian_at=hat,
     distance=distance, steps=steps, gain=gain, full_hessians=full_hessians[],
     subset_hessians=subset_hessians[], fcalls=fcalls[], gcalls=gcalls,
     escapes=escapes, saddle=saddle, ladder_tried=ladder_tried,
     ladder_gain=ladder_gain, saddle_walked=saddle_walked,
     saddle_gain=saddle_gain, probe=probed, history=history, handback=false)
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
    ctsem_endgame(objective, values; gradient_method, flat_rtol, probe_lengths,
        progress, progress_overwrite, progress_sink, progress_label)

The certification's numbers at a point the optimiser left without its finish --
an iteration cap, a stall, a finish that could not form a Hessian -- so R can
certify it the way it certifies a finished one: the exact Hessian at `values`,
what any negative curvature there is worth along the ladder, and what the flat
directions are worth. No step is taken. Whether to go
on from here is R's decision, and it goes on by resuming the optimiser, whose
own finish takes the steps.

The Hessian is of the MAXIMISED objective, and a 1x1 zero when it could not be
formed: an empty matrix would hang the bridge, as would an empty probe
direction, which is `[0.0]` when `probe_ran` is false. The counts are named as
`ctsem_optimize` names its own, so whatever tallies engine runs counts this
one too.

This is exactly the path a certification takes when the optimiser stopped
without a finish having run -- so on the Laplace route it is the whole `2 npar`
gradient loop with nothing to show for it until it returns, the same freeze a
finish already reports through `hessian_progress`. `progress`/`progress_sink`/
`progress_overwrite`/`progress_label`/`progress_every` mirror `ctsem_optimize`'s
own vocabulary, for the same `CTSEMProgress` and the same sink on the R side
(`.ctBackendEndgameAt()`); `progress = false` (the default) reports nothing, as
every caller before this did. `progress_every` is exposed mainly for a
fit-free test to force a cadence the certification's own small models would
not otherwise reach in the time they take.
"""
function ctsem_endgame(objective::CTSEMOptimisable, values::AbstractVector;
        gradient_method=:adjoint, flat_rtol::Real=_CTSEM_FLAT_RTOL,
        probe_lengths=_CTSEM_FLAT_PROBE_LENGTHS,
        progress::Bool=false, progress_overwrite::Bool=true, progress_sink=nothing,
        progress_label::AbstractString="certify", progress_every::Real=0.0)
    x = collect(Float64, values)
    fg! = _ctsem_trial_closure(objective, gradient_method)
    G = zeros(length(x))
    f = fg!(0.0, G, x)
    reporter = _ctsem_progress_reporter(progress, progress_label,
        progress_overwrite, progress_sink, progress_every)
    out = _ctsem_newton_finish(objective, x, f, G, fg!; take_steps=false,
        probe=true, curvature=:exact, flat_rtol=flat_rtol,
        probe_lengths=collect(Float64, probe_lengths), reporter=reporter)
    # Only when something was actually shown: a certification is usually fast
    # (the cheap, non-Laplace Hessian is one chunked ForwardDiff call), and a
    # "done" line for every one of those would be the paragraph this whole
    # effort exists to avoid rather than the phrase. `reporter` is always a
    # struct now (`_ctsem_progress_reporter`), so `.lines > 0` alone already
    # implies it was enabled -- `_due` gates every increment on that -- but
    # the explicit check says so rather than leaning on it.
    reporter.enabled && reporter.lines > 0 &&
        _progress_done(reporter, "hessian formed")
    probed = out.probe
    return (minimizer=x, maximum_loglik=-f, gradient=-G,
        hessian=out.hessian === nothing ? zeros(1, 1) : out.hessian,
        iterations=0, f_calls=1 + out.fcalls, g_calls=1 + out.gcalls,
        newton_steps=0, newton_hessians=out.full_hessians,
        newton_saddle_walked=out.saddle_walked,
        newton_saddle_gain=out.saddle_gain,
        probe_ran=!isempty(probed.direction), probe_gain=probed.gain,
        probe_length=probed.length, probe_longest=probed.longest,
        probe_direction=isempty(probed.direction) ? [0.0] : probed.direction)
end

export ctsem_endgame

"""
    ctsem_hessian_progress(objective, values; forward, chunk, progress,
        progress_sink, progress_overwrite, progress_label, progress_every)

`ctsem_hessian` (or, with `forward = true`, `ctsem_hessian_forward`), reporting
through the same `CTSEMProgress`/`_progress_fraction` line
`_ctsem_newton_finish`'s own Hessian already uses, for a caller with no
enclosing optimiser or certification to lend it a reporter. Today that is only
`.ctBackendHessian()` (R/ctBackendUncertainty.R), the standalone
`ctFitUncertainty()`/certification Hessian reached outside any running fit: it
used to build its own R-side rate limiter and a hand-written copy of this
line's text, one more way of saying the same thing `ctsem_endgame` already
says for the same reason (no finish to borrow a reporter from).

The `progress*` keywords are `ctsem_optimize`/`ctsem_endgame`'s own vocabulary
(`_ctsem_progress_reporter`); `progress = false`, the default, reports
nothing, exactly as a bare `ctsem_hessian(objective, values)` call always has.
`forward` selects `ctsem_hessian_forward` for the one route that needs it (a
model with a sampled TI predictor value, forced to forward mode) -- a single
ForwardDiff.hessian call, so `progress` is accepted and has nothing to report,
same as calling that function directly.
"""
function ctsem_hessian_progress(objective, values::AbstractVector;
        forward::Bool=false, chunk::Integer=0,
        progress::Bool=false, progress_sink=nothing, progress_overwrite::Bool=true,
        progress_label::AbstractString="hessian", progress_every::Real=0.0)
    reporter = _ctsem_progress_reporter(progress, progress_label,
        progress_overwrite, progress_sink, progress_every)
    cb = _ctsem_hessian_progress(reporter)
    forward ? ctsem_hessian_forward(objective, values; chunk=chunk, progress=cb) :
        ctsem_hessian(objective, values; chunk=chunk, progress=cb)
end

export ctsem_hessian_progress

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

"""
    ctsem_saddle_probe(objective, values, hessian; ascent, ladder)

The saddle ladder (`_ctsem_saddle_ladder`) along every direction of negative
curvature in `hessian` -- of the MAXIMISED objective, as R holds it -- at
`values`, measured and no step taken: for a certification R assembles from a
Hessian it already holds (`.ctBackendCertification()`, the Laplace
continuation). Its split is `_ctsem_information_split`'s, the rule R's
`.ctBackendInformationSplit()` applies. `ascent`, the gradient there, orders
the two sides of each direction and nothing else. Named as `ctsem_optimize`
names the finish's own, so R reads both one way: `newton_saddle_walked` the
directions measured, `newton_saddle_gain` the best gain found (0 for none).
"""
function ctsem_saddle_probe(objective, values::AbstractVector,
        hessian::AbstractMatrix; ascent=nothing, ladder=_CTSEM_SADDLE_LADDER)
    x = collect(Float64, values)
    split = _ctsem_information_split(-Matrix{Float64}(hessian))
    value_at = y -> _ctsem_probe_value(objective, y)
    base = (split === nothing || !any(split.negative)) ? -Inf : value_at(x)
    isfinite(base) || return (newton_saddle_walked=0, newton_saddle_gain=NaN)
    a = ascent === nothing ? zeros(length(x)) : collect(Float64, ascent)
    tried = _ctsem_saddle_ladder(value_at, x, base, split, a;
        ladder=collect(Float64, ladder), enough=Inf)
    return (newton_saddle_walked=tried.walked,
        newton_saddle_gain=tried.best === nothing ? 0.0 : tried.best.value - base)
end

export ctsem_saddle_probe
