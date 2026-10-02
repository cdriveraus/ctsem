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
# `diagonal`, when given, is the initial inverse Hessian itself, one scale per
# coordinate (`_ctsem_lbfgs_diagonal!`), in place of the metric's shape times
# the latest pair's one secant ratio.
function _ctsem_lbfgs_hmul(M::CTSEMLBFGSMemory, q::AbstractVector, h0::Float64,
        dinv::Vector{Float64}, diagonal=nothing)
    k = length(M.S)
    k == 0 && return h0 .* dinv .* q
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
    else
        r .= diagonal .* r
    end
    @inbounds for i in 1:k
        b = M.rho[i] * dot(M.Y[i], r)
        r .+= (a[i] - b) .* M.S[i]
    end
    r
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
function _ctsem_lbfgs_diagonal!(B, s::Vector{Float64}, y::Vector{Float64},
        dinv::Vector{Float64})
    sy = dot(s, y)
    B === nothing && (B = [1 / d for d in dinv])
    scale = sum(i -> y[i]^2 / B[i], eachindex(y)) / sy
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
"""
function _ctsem_lbfgs(fg!, x0::AbstractVector; memory::Integer=20,
        metric=nothing, initial_alpha::Real=0.1, maxiter::Integer=1000,
        g_tol::Real=1e-8, f_tol::Real=0.0, x_tol::Real=0.0,
        callback=nothing, directional=nothing, batch=nothing,
        c1::Real=1e-4, maxbacktrack::Integer=40, iteration0::Integer=0,
        diagonal::Bool=false, nonmonotone::Real=0.0, gll::Integer=0)
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
    # `diagonal`: the initial inverse Hessian learned per coordinate from the
    # pairs (`_ctsem_lbfgs_diagonal!`), `B` its inverse, rather than the
    # metric's shape times one secant ratio. `nonmonotone`: Zhang & Hager's
    # (2004) line search, which accepts a step against `C`, a running average
    # of the values the steps reached weighted by `eta`, rather than against the
    # last one -- so a step may rise above where it started while it stays
    # below the recent level, which lets a step cross a curved valley that a
    # monotone search would cut short. `eta = 0` makes `C` the last value: the
    # search is then exactly the monotone Armijo one.
    B = nothing; Dinv = nothing
    eta = clamp(Float64(nonmonotone), 0.0, 1.0)
    C = f; Q = 1.0; recent = [f]
    # `gll`: Grippo, Lampariello & Lucidi's (1986) rule instead -- a step is
    # judged against the worst of the last `gll` values the steps reached, the
    # acceptance sgd() uses (`_ctsem_sgd`). `recent` is set beside C.
    reference() = gll > 0 ? maximum(@view recent[max(1, end - gll + 1):end]) : C
    while !stopped && !gconv && iteration < maxiter
        s = -_ctsem_lbfgs_hmul(M, G, h0, dinv, Dinv)
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
            B = nothing; Dinv = nothing; C = f; Q = 1.0; recent = [f]
            continue
        end
        # The first trial carries the gradient too: most steps are accepted
        # whole, and then the iteration has cost one gradient and nothing else.
        alpha = 1.0
        xn = x .+ s
        Gn = similar(G)
        fn = evaluate!(0.0, Gn, xn); fcalls += 1; gcalls += 1
        have_gradient = true
        accepted = isfinite(fn) && fn <= reference() + c1 * alpha * dphi
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
            accepted = isfinite(fn) && fn <= reference() + c1 * alpha * dphi
        end
        if !accepted
            # A stale memory is the usual cause; drop it once, then give up.
            if !retried && !isempty(M.S)
                retried = true
                _ctsem_lbfgs_reset!(M)
                h0 = Float64(initial_alpha) / max(metric_norm(G), eps())
                B = nothing; Dinv = nothing; C = f; Q = 1.0; recent = [f]
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
        dG = Gn .- G
        if _ctsem_lbfgs_push!(M, step, dG) && diagonal
            B = _ctsem_lbfgs_diagonal!(B, step, dG, dinv)
            Dinv = 1 ./ B
        end
        x = xn; f = fn; G = Gn
        push!(recent, f)
        Qn = eta * Q + 1
        C = (eta * Q * C + f) / Qn; Q = Qn
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
                _ctsem_batch_step!(batch, x, G, q -> _ctsem_lbfgs_hmul(M, q, h0, dinv, Dinv), iteration)
            if _ctsem_batch_full(batch)
                sizes = batch.sizes; its = batch.iterations
                batch = nothing
            end
            f = evaluate!(0.0, G, x); fcalls += 1; gcalls += 1
            gconv = maximum(abs, G; init=0.0) <= g_tol
            # The reference level belongs to the objective it averaged.
            C = f; Q = 1.0; recent = [f]
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
    gmemory = 0.6
    gmemory2 = 0.0
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
        gmemory2 = gmemory * min(i / warmup, 1)^(1 / 8)
        rm2 = 0.9 * min(i / warmup, 1)^(1 / 8)
        oldgmid .= gmid
        oldgsmooth .= gsmooth
        @inbounds for k in 1:n
            gmid[k] = (g[k] + gnew[k]) / 2
            gsmooth[k] = gsmooth[k] * gmemory2 + (1 - gmemory2) * gnew[k]
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
        (i > 25 && i % 20 == 0) && (gmemory = clamp(gmemory, 0.8, 0.95))
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
                        gsmooth[k] = gsmooth[k] * gmemory2^2 + (1 - gmemory2^2) * g[k]
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
        stopped, Int[], Int[])
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
# * the marginal route's is forward-over-adjoint, one lane per parameter at
#   about a gradient each: a few gradients on a small model (6 on a
#   1000-subject panel, dev1), so L-BFGS hands over early, and there every step
#   uses the exact Hessian, refreshed while the steps contract slowly. On a
#   large one it is not few: about 750 gradients, a quarter of an hour, on a
#   715-parameter model. A refresh pays only while a Hessian costs less than
#   the steps it saves, which are at most the finish's budget of `maxit`, so
#   above `maxit` parameters the marginal route takes the chord below
#   (`_ctsem_finish_curvature`);
# * the Laplace route's is central differences of the gradient, `2 npar`
#   gradients, and an exact finish there made fits slower (0.6x with it, 1.3 to
#   1.9x without, dev1). So L-BFGS runs to its own stopping rule there, and the
#   finish takes one Hessian at the hand-over, reuses it for its steps (the
#   chord) and keeps it as the final one unless the steps moved the estimate by
#   more than `_CTSEM_HESSIAN_REUSE_SE` standard errors: one Hessian per fit in
#   the common case.
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
(`_ctsem_flat_ladder`). On the marginal objective, the likelihood itself, a
gain along such a direction is a gain. Not on the Laplace objective: there the
approximation is least reliable exactly where a direction goes flat -- a
random-effect scale collapsing, a unit's curvature going singular -- and moving
along one walked AnomAuth S1 from its optimum into the spurious basin, 25
Laplace nats up and 1.8 exact nats down (decision 1 of
review/OPTIM-consolidation-plan-2026-09-25). Nor on the state-explicit
objective, which is sampled, never optimised.
"""
_ctsem_flat_escape(::CTSEMObjective) = true
_ctsem_flat_escape(::Any) = false

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
    _ctsem_finish_curvature(objective[, npar, maxit])

The route's own curvature for the finish: `:exact` where a Hessian costs no
more than the `maxit` steps a refresh could save (`_ctsem_hessian_cost`),
`:chord` where it costs more -- always on the Laplace route, and on the
marginal route above `maxit` parameters -- and `nothing` for a route the
finish does not run on: a pinned profile point, whose wrapper has no Hessian
of its own, and the state-explicit target, which is never certified.
"""
_ctsem_finish_curvature(::CTSEMObjective) = :exact
_ctsem_finish_curvature(::CTSEMLaplaceObjective) = :chord
_ctsem_finish_curvature(::Any) = nothing
_ctsem_finish_curvature(o, npar::Integer, maxit::Integer) = _ctsem_finish_curvature(o)
_ctsem_finish_curvature(o::CTSEMObjective, npar::Integer, maxit::Integer) =
    _ctsem_hessian_cost(o, npar) <= maxit ? :exact : :chord

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
    _ctsem_secant_along(H, s, y)

`H` with its curvature along the step `s` replaced by the secant's, `s'y / s's`
(`y` the change in gradient over the step), and every other direction's kept:
a symmetric rank-one correction, which is what lets the chord's copy of its
Hessian follow a curvature that changes along the direction it walks. Returned
unchanged when the secant is not a positive, finite curvature -- noise, or a
step into negative curvature, which the chord then keeps from its Hessian.
"""
function _ctsem_secant_along(H::AbstractMatrix, s::AbstractVector, y::AbstractVector)
    ss = dot(s, s)
    (isfinite(ss) && ss > 0) || return H
    kappa = dot(s, y) / ss
    (isfinite(kappa) && kappa > 0) || return H
    current = dot(s, H * s) / ss
    isfinite(current) || return H
    H .+ ((kappa - current) / ss) .* (s * transpose(s))
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
        q = sum(i -> c[i]^2 / (curvature[i] + mu)^3, eachindex(c); init=0.0)
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
  fast), with the steps' copy of it corrected along any step that is slow in
  that sense by the secant curvature there (`_ctsem_secant_along`);
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
that, as is a chord that used up its steps -- `maxit`, or as many as its
Hessian cost in gradients (`_ctsem_hessian_cost`) if more -- at most five
times; on the exact variant the steps continue on a fresh Hessian each until it
agrees.

At a saddle -- negative curvature in the final Hessian -- the ladder along the
most negative curvature is tried (`_ctsem_saddle_ladder`), and the finish
continues from the better point with the exact Hessian there, at most
`max_escapes` times: it never returns a point with a direction of negative
curvature it did not try, short of that cap. With no negative curvature but a
direction the curvature calls flat, or a saddle whose ladder found nothing
worth having, the flat ladder (`_ctsem_flat_ladder`) is tried, when
`flat_escape` (the fit sets it where `_ctsem_flat_escape` allows), within the
same cap, the steps continuing on the exact Hessian there with a fresh step
budget. Either escape is taken only when it gains more than `escape_gain`,
which the fit sets to its certification tolerance. Then, when `probe`, the
flat-direction probe (`_ctsem_flat_probe`) at the final point.

`take_steps = false` gives the certification's numbers alone: the Hessian at `x`
and the probe, with no step, no refresh and no ladder, for a point the optimiser
left without a finish (`ctsem_endgame`).

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
`saddle`, `ladder_tried`, `ladder_gain`); the probe record; whether it handed
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
        escape_gain::Real=0.0, flat_escape::Bool=false,
        reporter=nothing)
    curvature in (:exact, :chord, :subset) || throw(ArgumentError(
        "newton curvature must be exact, chord or subset, got $(curvature)"))
    valueof = value_at === nothing ? (y -> _ctsem_probe_value(objective, y)) :
        value_at
    x = collect(Float64, x0); f = Float64(f0); G = collect(Float64, G0)
    full_hessians = 0; subset_hessians = 0
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
                max(abs(v), 1e-8 * lmax) for v in E.values])
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
    # for the step taken, which the trust region judges it by.
    function newton(Hm, Gv, radius)
        local M = model(Hm)
        M.usable || return nothing
        local lam = M.curvature
        local c = transpose(M.vectors) * Gv
        local mu = _ctsem_trust_multiplier(lam, c, radius)
        local d = c ./ (lam .+ mu)
        (step=-(M.vectors * d), gain=0.5 * sum(abs2.(c) ./ lam; init=0.0),
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
    function trust(Hm, nt, radius)
        local floor = max(abs(f), 1.0) * eps()
        local at = nt
        local r = radius
        local k = 0
        while true
            local xn = x .+ at.step
            local fn = fg!(0.0, nothing, xn)
            fcalls += 1
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
    report(g) = callback === nothing || callback(CTSEMIterate(iteration0 + steps,
        f, maximum(abs, G; init=0.0), g))
    H = step_hessian(x)
    # What the steps are taken against: `H` itself, except that the chord
    # corrects its copy along the direction it walks (see the refresh below),
    # while `H`, exact where it was evaluated, is what the reuse test and the
    # certification read.
    Hs = H
    hat = copy(x)                       # where `H` was evaluated
    exact = curvature !== :subset       # `H` is the exact Hessian at `hat`
    steps = 0; escapes = 0; gain = Inf
    saddle = false; ladder_tried = false; ladder_gain = 0.0
    firstgain = Inf; firstwhole = true; handed = false
    if H === nothing
        return (x=x, f=f, G=G, hessian=nothing, hessian_at=hat, distance=NaN,
            steps=0, gain=Inf, full_hessians=full_hessians,
            subset_hessians=subset_hessians, fcalls=0, gcalls=0, escapes=0,
            saddle=false, ladder_tried=false, ladder_gain=0.0,
            probe=_ctsem_no_probe(), history=history, handback=false)
    end
    at_x = exact                        # `H` is the exact Hessian at `x`
    # The trust region, in raw units. Each Hessian starts without one -- its
    # first step is the Newton step -- and the region is set by what the steps
    # on it then show.
    radius = Inf
    # The steps a Hessian is walked for: `maxit`, and on the chord and the
    # subset as many as the Hessian cost (`_ctsem_hessian_cost`), since a
    # step there costs a gradient and a new Hessian is not worth buying
    # before the steps have spent as much as the last one did. A chord that
    # uses them up, or converges somewhere its Hessian no longer describes,
    # is walked again on the exact Hessian there (below), at most
    # `max_rewalks` times.
    per_hessian = curvature === :exact ? Int(maxit) :
        max(Int(maxit), _ctsem_hessian_cost(objective, length(x)))
    budget = per_hessian; rewalks = 0; max_rewalks = 5
    while true                          # once, and again after an escape or a rewalk
        converged = false
        prevgain = Inf; prevwhole = true; prevrho = 1.0
        while take_steps && steps < budget
            nt = newton(Hs, G, radius)
            nt === nothing && break
            gain = nt.gain
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
            trial = trust(Hs, nt, radius)
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
            dx = trial.x .- x; dg = Gn .- G
            x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
            steps == 1 && (firstwhole = trial.nt.whole)
            remember!("newton", gain, trial.nt.fraction)
            report(gain)
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
            # The chord's refresh is a free one: its copy of the matrix takes
            # the secant's curvature along the step it just took, which the
            # gradients at both ends measure, and keeps the exact Hessian's
            # everywhere else (`_ctsem_secant_along`). Without it the chord
            # walked a ray whose curvature decays along it -- a random-effect sd
            # heading toward zero -- at the curvature of where it started, every
            # step shorter than the last: on AnomAuth S1 from its default start
            # (bench, dev1), 30 chord steps for gains of 4e-5 down to 3e-7, then
            # five exact Hessians to finish, seven in the stage.
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
                curvature === :chord && slow && (Hs = _ctsem_secant_along(Hs, dx, dg))
                prevgain = gain
            end
        end
        # Handed back: the point the step reached, and no final Hessian, since
        # L-BFGS goes on from here and the finish runs again where it stops.
        handed && return (x=x, f=f, G=G, hessian=nothing, hessian_at=hat,
            distance=NaN, steps=steps, gain=gain, full_hessians=full_hessians,
            subset_hessians=subset_hessians, fcalls=fcalls, gcalls=gcalls,
            escapes=0, saddle=false, ladder_tried=false, ladder_gain=0.0,
            probe=_ctsem_no_probe(), history=history, handback=true)
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
                if take_steps && rewalks < max_rewalks
                    rewalks += 1
                    budget = steps + per_hessian
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
                    trial = trust(H, nt, radius)
                    trial.ok || break
                    Gn = similar(G)
                    fg!(nothing, Gn, trial.x); gcalls += 1
                    x = trial.x; f = trial.f; G = Gn; steps += 1; at_x = false
                    remember!("exact", gain, trial.nt.fraction)
                    report(gain)
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
        take_steps || break
        split = _ctsem_information_split(H; rtol=flat_rtol, negative=negative)
        split === nothing && break
        # An escape has to gain more than `escape_gain`, the certification's
        # tolerance: less is a move the verdict cannot tell from staying, and it
        # spends one of the `max_escapes`. On test-julia-convergence.R's random
        # walks the Hessian on the white-noise plateau showed a noise-level
        # negative curvature, and three saddle escapes of 0.03 raw units, gaining
        # nothing at eight digits, used the cap before anything else was tried.
        if any(split.negative)
            saddle = true
            escapes >= max_escapes && break
            ladder_tried = true
            tried = _ctsem_saddle_ladder(valueof, x, -f, split, -G; ladder=ladder)
            fcalls += tried.evaluations
            if tried.best !== nothing
                Gn = similar(G)
                fn = fg!(0.0, Gn, tried.best.point); fcalls += 1; gcalls += 1
                # The point the ladder measured, through the optimiser's own
                # trial path: a point the route refuses a gradient at is not one
                # to continue from.
                if isfinite(fn) && f - fn > escape_gain
                    ladder_gain += f - fn
                    x = collect(Float64, tried.best.point); f = fn; G = Gn
                    steps += 1; escapes += 1
                    remember!("saddle", NaN, tried.best.length)
                    report(NaN)
                    H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
                    saddle = false; ladder_tried = false
                    H === nothing && break
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
        fcalls += walked.fcalls; gcalls += walked.gcalls
        (walked.best !== nothing && f - walked.best.f > escape_gain) || break
        x = collect(Float64, walked.best.point); f = walked.best.f
        G = walked.best.G
        steps += 1; escapes += 1
        remember!("flat", NaN, walked.best.length)
        report(NaN)
        H = hess(x); Hs = H; hat = copy(x); exact = true; at_x = true
        saddle = false; ladder_tried = false
        H === nothing && break
        # A new point and a new Hessian: walked as the first one was.
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
            fcalls += probed.evaluations
        end
    end
    (x=x, f=f, G=G, hessian=H === nothing ? nothing : -H, hessian_at=hat,
     distance=distance, steps=steps, gain=gain, full_hessians=full_hessians,
     subset_hessians=subset_hessians, fcalls=fcalls, gcalls=gcalls,
     escapes=escapes, saddle=saddle, ladder_tried=ladder_tried,
     ladder_gain=ladder_gain, probe=probed, history=history, handback=false)
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
certify it the way it certifies a finished one: the exact Hessian at `values`
and what the flat directions are worth there. No step is taken. Whether to go
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
