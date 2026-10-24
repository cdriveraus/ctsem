"""
Doubly adaptive importance sampling: the normal moved to the posterior's
moments.

`uncertainty = 'checked'` tests the Hessian's normal with gradients at draws
from it (R's `.ctOptimNormalCheck()`), takes one step of the posterior's mean
from the pairs' intercept and one of its widths from their odd and even
parts, and hands the directions it names to importance sampling with the
rest held at that normal. Where the posterior's mass is far from the mode --
a 75-parameter time series whose bulk sits 95 normal sds from it, 41 of its
directions three times wider than the Hessian says -- one step is not
enough, importance sampling in 67 named directions kept 2 draws in 2,000,
and what was reported was the Hessian's normal with a warning.

This iterates that step with damping until the normal sits on the
posterior's moments: van den Boom, Cremaschi and Thiery, "Doubly adaptive
importance sampling", JCGS 35 (2026), arXiv:2404.18556. In the normal's
whitened coordinates `z` (`theta = centre + A z`, `A A'` the covariance),
with `q = N(mu, Gamma)` the current normal there, each iteration draws `S`
points from `q` in antithetic pairs, takes their log densities and
gradients, and forms `Phi = log p - log q`. The damping `gamma` is the
largest value in (0, 1] whose weights `exp(gamma Phi)` keep an effective
sample of half the draws (bisection on the same draws). Stein's lemma for
the Gaussian `q` then gives the damped target `q exp(gamma Phi)` its moments
from the weighted gradients,

    mu'    = mu    + gamma E_w[Gamma grad Phi]
    Gamma' = Gamma + gamma Cov_w(Gamma grad Phi, z),

with `Gamma grad Phi(z) = Gamma grad log p(z) + (z - mu)`, and `q` moves to
them. At `gamma = 1` that is the posterior's mean and covariance by
importance sampling; at `gamma -> 0` it is natural-gradient Gaussian
variational inference (Khan and Rue 2023's Bayesian learning rule; the fixed
point `E_q[grad log p] = 0`, `Gamma^-1 = E_q[-Hessian]` of Opper and
Archambeau 2009). Below one, which is where it stays in many dimensions --
a Gaussian cannot weight a skewed 75-dimensional posterior to half its draws
-- the fixed point is the bridge `q^(1 - gamma) p^gamma`'s moments: the mean
the posterior's, the sds a little under it (review/SAMPLING-smc-2026-10-23.md:
on the time series every mean within 0.4 reference sds where the check's
normal was 2.5 off, the sds 0.75 of the posterior's).

Three things the paper's 10^5 draws an iteration never met at 64-256:

  * The iterate is noisy where it is not wrong: its per-coordinate sds come
    out right while its eigenstructure is noise. So the answer is the
    average of the iterates over a window at the fixed point (Polyak).
  * When to stop. The answer is held to stated precisions
    (`_DAIS_PRECISION`): every parameter's mean within `mean_sd` (0.2) of
    its sd and its variance within `var` (40%) of itself, all parameters at
    once at `level` (5%), and the averaged covariance's Monte Carlo noise
    distorting no direction's width by more than `width` (a third), which
    the normal needs to serve as a sampler's metric. Each iteration takes the
    shortest trailing window of iterates whose draws could show a departure
    from those (`_dais_window`): effective draws at least `k / width^2` for
    `k` free dimensions, and for every parameter a remaining distance of its
    mean or variance of that size visible above the random walk its damped
    steps' noise makes. It stops where, over that window, no parameter's mean
    or variance has drifted beyond the bound for testing them all at once,
    and the last natural gradient `E_w[Gamma grad Phi]` is within twice its
    own. The answer is the window's average. So a run that starts near its
    fixed point stops as soon as its draws can say so, and one still
    travelling goes on. A test that cannot be read (a batch with no spread)
    is unmet.

    The precision is per parameter because an aggregate hides a direction:
    held to a tenth of an sd a coordinate as a norm over the time series'
    75, the rule allowed most of 0.9 sd in one, stopped at 87 iterations
    where the iterate was still moving, and its largest mean error was 0.86
    sd where the run left to go on had 0.43. The ELBO summarises every
    coordinate too, and it no longer decides. The level is a fifth of an sd
    because that is about what DAIS delivers (`_DAIS_PRECISION`): a tenth
    asked for more draws than either hard posterior gave in 150 iterations,
    for no better answer.

    Until 2026-10-24 the windows were fixed at 20 iterations, with three
    passes in a row and 20 more iterations averaged: 61 before any stop,
    however early the iterate had settled. The 266-parameter scale fit
    stopped at exactly 61. Four rules before that failed on the time series:
    stopping when `gamma` reached 1, which it never does in many dimensions;
    when the damped step was small, which it is whenever `gamma` is; when the
    undamped step was within a crude noise floor, which it always was; and on
    the ELBO's drift over 5 iterations, which hid 0.4 nats an iteration. Each
    stopped where too few draws could see the departure. The windows here
    are as long as seeing it takes, and no longer.
  * A covariance step that is not positive definite halves `gamma` on the
    same draws until it is.

A direction the data do not identify has no fixed point: the gradient along
it is zero, so each step multiplies its variance by `1 + gamma`. So the
directions the flat screen confirmed (`hold`) are held at the normal's width
-- the moves and the covariance are projected off them -- and a covariance
that still grows past `maxgrowth` times the start's in any direction ends the
run as diverged, `used = false`, and the check's normal stands. A draw the
density cannot use carries no weight and no gradient (zero weight times an
infinite gradient was NaN, and stopped a fit); a batch with fewer than two
usable draws ends the run the same way.

Where `'checked'` uses a coreset (`coreset = true`, on a Laplace objective),
DAIS's draws go through one of its own (`CTSEMISCoreset`,
coreset_sampler.jl): its gradient batch runs on the coreset where that is
expected to cost less, checked at a few draws in full every batch, and an
iteration evaluated in full teaches it, since a value batch at the same draws
reads the unit terms the gradients computed, for nothing, and is where it
learns them and chooses. Its own, not the importance sampler's: DAIS's draws
vary in every direction and the sampler's in the few the check named, so
DAIS's sit far below the sampler's in density, and a coreset weighs its
training draws by their distance below the best. Fed DAIS's draws, the
sampler's coreset priced nearly every choice at a failure rate of one half,
took the cheapest size and failed: 27 of 28 choices at 1,000 units on the
strong shell. Elsewhere the draws are evaluated at once, on the pool or on
copies of the objective (`_each_point!`).

Moved here from R (2026-10-24). In R each iteration crossed the bridge for
its batch and did its arithmetic and the coreset's refits in R: on the
266-parameter scale fit the R process used nine hours of CPU against the
engine's four, and `'checked'` had not finished in nine and a half. Here
nothing crosses until the answer does.
"""

using LinearAlgebra
using Random

export ctsem_dais

"""Draws an iteration: the dimension, at least 64 and at most 256, even."""
_dais_draws(k::Int) = 2 * ceil(Int, min(256, max(64, k)) / 2)

"""
What DAIS's answer is held to (see the file comment): each parameter's mean
within `mean_sd` of its sd and its variance within `var` of itself, every
parameter at once at `level`, and no direction's width distorted by more than
`width` by the averaged covariance's Monte Carlo noise. The last is the
spread of a sample covariance's eigenvalues, `(1 +- sqrt(k / n))^2` for `n`
effective draws in `k` dimensions, so it asks for `k / width^2` of them; it is
what the normal needs to serve as a sampler's metric.

A fifth of an sd is about what DAIS delivers: its damped fixed point sits up
to 0.4 sd from the posterior's mean on the time series, and 'checked''s
recheck and importance sampler correct the directions that matter after it.
From one seed and one start (local, 2026-10-24), against the old rule's 142
iterations on the time series and 85 on the 250-unit shell, both with the
largest mean error 0.43 and 0.77 sd: at 0.1, neither settled in 150, for
0.40 and 0.89; at 0.2, 123 and 112 iterations, 0.42 and 0.73; at 0.3, 105
and 54, 0.52 and 1.03, the shell's smallest sd ratio 0.64 against 0.83.
"""
const _DAIS_PRECISION = (mean_sd=0.2, var=0.4, level=0.05, width=1 / 3)

"""
    _dais_window(ess, gamma, se2p, se2v, need, z, tol_mean, tol_var)

The shortest trailing window of iterates, at least three, whose draws could
show a departure from what the answer is held to: effective draws (`ess`) at
least `need`, and, for every parameter, a remaining distance of its mean of
`tol_mean` of its sd, and of its variance of `tol_var` of itself, visible
above the random walk its damped steps' noise makes. A remaining distance
`d` moves the mean `d` times the window's summed damping (`gamma`), and the
steps' noise (`se2p`, `se2v`: a vector per iteration, each step's own, in
the parameter's sd units) walks `sqrt(sum(gamma^2 se2))`; visible means
beyond `z + 1` of that walk, `z` the bound a parameter's drift must stay
within when every parameter is tested at once. Returns the window's length
(0 where the run is not yet long enough) and each parameter's walk over it.
"""
function _dais_window(ess::Vector{Float64}, gamma::Vector{Float64},
    se2p::Vector{Vector{Float64}}, se2v::Vector{Vector{Float64}}, need::Real, z::Real,
    tol_mean::Real, tol_var::Real)
    n = length(ess)
    np = n == 0 ? 0 : length(se2p[1])
    walkp = zeros(np)
    walkv = zeros(np)
    draws = 0.0
    damp = 0.0
    for m in 1:n
        t = n - m + 1
        draws += ess[t]
        damp += gamma[t]
        g2 = gamma[t]^2
        walkp .+= g2 .* se2p[t]
        walkv .+= g2 .* se2v[t]
        m >= 3 || continue
        draws >= need || continue
        (z + 1) * sqrt(maximum(walkp; init=0.0)) <= tol_mean * damp || continue
        (z + 1) * sqrt(maximum(walkv; init=0.0)) <= tol_var * damp || continue
        return (m=m, walkp=sqrt.(walkp), walkv=sqrt.(walkv))
    end
    return (m=0, walkp=walkp, walkv=walkv)
end

"""
    _dais_evaluate!(vals, grads, objective, coreset, X, mode)

Each column of `X`'s log density and gradient into `vals` and the columns of
`grads`, `-Inf` and zeros where a draw cannot be used, and the full
evaluations spent, a gradient as two. Through DAIS's coreset where there is
one: its gradient batch, checked at a few draws in full every time (`check`),
since DAIS has no value batches to check it with; and where that ran in full,
a value batch at the same draws, which reads the batch's unit terms for
nothing and is where the coreset learns them and is chosen. Otherwise on the
objective, the draws at once (`_each_point!`, `mode`).
"""
function _dais_evaluate!(vals::Vector{Float64}, grads::Matrix{Float64}, objective,
    ::Nothing, X::Matrix{Float64}, mode::Symbol)
    fill!(vals, -Inf)
    fill!(grads, 0.0)
    _each_point!(objective, size(X, 2), mode) do i, obj, _
        r = _ctsem_probe_value_gradient(obj, X[:, i])
        if r.gradient !== nothing && isfinite(r.value) && all(isfinite, r.gradient)
            vals[i] = r.value
            grads[:, i] .= r.gradient
        end
        return nothing
    end
    return 2.0 * size(X, 2)
end

function _dais_evaluate!(vals::Vector{Float64}, grads::Matrix{Float64}, objective,
    cs::CTSEMISCoreset, X::Matrix{Float64}, mode::Symbol)
    fill!(vals, -Inf)
    fill!(grads, 0.0)
    before = cs.cost
    r = ctsem_is_coreset_gradbatch!(cs, X; check=true)
    _is_cached(cs, X) && ctsem_is_coreset_batch!(cs, X)
    for i in axes(X, 2)
        v = r.value[i]
        g = view(r.gradient, :, i)
        if isfinite(v) && v > -1e100 && all(isfinite, g)
            vals[i] = v
            grads[:, i] .= g
        end
    end
    return cs.cost - before
end

"""
    ctsem_dais(objective, centre, cov; coreset, hold, draws, maxiter, seed,
        rtol, maxgrowth, parallel, verbose, progress_sink, progress_overwrite)

Move the normal `N(centre, cov)` to the posterior's moments by doubly
adaptive importance sampling (see the file comment). `coreset = true` sends
the draws through a coreset of DAIS's own where `objective` is a sum of units
(`_smc_coreset_capable`). `hold`, a matrix whose columns are directions in
the parameters, holds them at the normal's width. `draws` (0: `_dais_draws`
of the dimension) and `maxiter` (150) bound the run. Returns `used`, the
settled normal (`mean`, `cov`), and the record: iterations, how many were
averaged, whether it settled, the last damping and effective sample, the
unusable draws, the draws an iteration, the gradients asked for, the full
evaluations they cost (`cost`, a gradient as two), whether a coreset served
them, with the coreset's last size, how many it chose and how many of its
checks failed, how the draws were evaluated (`parallel`: `_points_mode`'s
serial, pool or copies, or `coreset`), and the trace, one row an iteration
(`trace_names`); `used = false` with `diverged` (and `growth`) or with no
usable draws, and `cov = nothing`.
"""
function ctsem_dais(objective, centre::AbstractVector, cov::AbstractMatrix;
    coreset::Bool=false, kwargs...)
    # DAIS normalises its weights within each iteration's batch, so a draw's
    # relevance is measured from the batch's best (`_is_top`).
    cs = coreset && _smc_coreset_capable(objective) ?
        ctsem_is_coreset(objective; top=:batch) : nothing
    return _dais(objective, cs, centre, cov; kwargs...)
end

function _dais(objective, coreset, centre::AbstractVector, cov::AbstractMatrix;
    hold::Union{Nothing,AbstractMatrix}=nothing, draws::Integer=0, maxiter::Integer=150,
    seed::Integer=20261024, rtol::Real=1e-10, maxgrowth::Real=1e4, parallel::Bool=true,
    verbose::Bool=false, progress_sink=nothing, progress_overwrite::Bool=true)
    # One body whatever keywords were passed: with `hold` absent and present it
    # was compiled twice, and a session's first default fit compiled for 2 s
    # the copy the package image did not hold.
    c0 = collect(Float64, centre)
    held = hold === nothing ? zeros(Float64, length(c0), 0) : Matrix{Float64}(hold)
    return _ctsem_dais(objective, coreset, c0, cov, held, Int(draws), Int(maxiter), Int(seed),
        Float64(rtol), Float64(maxgrowth), parallel, verbose, progress_sink,
        progress_overwrite)
end

function _ctsem_dais(objective, coreset, c0::Vector{Float64}, cov::AbstractMatrix,
    hold::Matrix{Float64}, draws::Int, maxiter::Int, seed::Int, rtol::Float64,
    maxgrowth::Float64, parallel::Bool, verbose::Bool, progress_sink,
    progress_overwrite::Bool)
    p = length(c0)
    C = Matrix{Float64}((cov .+ transpose(cov)) ./ 2)
    size(C) == (p, p) || throw(DimensionMismatch("cov must be $p x $p"))
    e = eigen(Symmetric(C))
    top = maximum(e.values)
    keep = findall(v -> v > rtol * top, e.values)
    k = length(keep)
    k > 0 || return _dais_stopped(:empty, NaN, 0, Vector{Float64}[], 0, 0, 0.0,
        coreset !== nothing)
    V = e.vectors[:, keep]
    sq = sqrt.(e.values[keep])
    A = V .* transpose(sq)
    # The held directions in the whitened coordinates, orthonormal: `A^+ h`
    # for each, `A^+ = D^-1/2 V'` on the kept eigenvectors.
    U = _dais_held_basis(V, sq, hold)
    offhold(M) = U === nothing ? M : M .- (M * U) * transpose(U)
    S = draws > 0 ? max(4, 2 * ((Int(draws) + 1) ÷ 2)) : _dais_draws(k)
    essfrac = 0.5
    kfree = k - (U === nothing ? 0 : size(U, 2))
    P = _DAIS_PRECISION
    need = kfree / P.width^2
    # The bound a parameter's drift over the window must stay within, every
    # parameter tested at once at `level` (Bonferroni).
    zp = sqrt(2) * erfinv(1 - P.level / p)
    mode = coreset === nothing ? _points_mode(objective, parallel) : :coreset
    rng = Random.Xoshiro(UInt64(seed))
    reporter = CTSEMProgress(verbose; label="DAIS", overwrite=progress_overwrite,
        sink=progress_sink)

    mu = zeros(k)
    Gamma = Matrix{Float64}(I, k, k)
    # One entry an iteration: the mean after it (with the start first), and
    # what the stopping rule reads. The covariances after the last `cap`
    # iterations, which bounds what a window's average can hold to 256 MB.
    mus = Vector{Vector{Float64}}([copy(mu)])
    cap = clamp(floor(Int, 2^28 / (8 * k^2)), 2, max(2, Int(maxiter)))
    Gs = Vector{Matrix{Float64}}(undef, cap)
    gammas = Float64[]
    esss = Float64[]
    # Each step's noise in each parameter's mean and variance, in its sd
    # units, and each parameter's variance after each iteration.
    se2ps = Vector{Vector{Float64}}()
    se2vs = Vector{Vector{Float64}}()
    varsv = Vector{Vector{Float64}}([vec(sum(abs2, A; dims=2))])
    window = 0
    gradients = 0
    cost = 0.0
    unusable = 0
    gamma = NaN
    essw = NaN
    trace = Vector{Vector{Float64}}()
    vals = Vector{Float64}(undef, S)
    grads = Matrix{Float64}(undef, p, S)
    half = S ÷ 2
    for iter in 1:Int(maxiter)
        _ctsem_interrupt_check()
        Lc = cholesky(Symmetric(Gamma)).L
        E = randn(rng, half, k)
        E = vcat(E, -E)
        Z = E * transpose(Lc) .+ transpose(mu)
        lq = vec(-0.5 .* sum(abs2, E; dims=2)) .- sum(log, diag(Lc)) .- 0.5 * k * log(2pi)
        X = c0 .+ A * transpose(Z)
        cost += _dais_evaluate!(vals, grads, objective, coreset, X, mode)
        gradients += S
        ok = isfinite.(vals)
        nok = count(ok)
        unusable += S - nok
        if nok < 2
            verbose && _due(reporter, true) && _emit(reporter,
                "  DAIS: $nok of $S draws usable at iteration $iter; the check's normal stands")
            _progress_done(reporter)
            return _dais_stopped(:unusable, NaN, unusable, trace, S, gradients, cost,
                coreset !== nothing)
        end
        G = transpose(grads) * A
        G[.!ok, :] .= 0.0
        Phi = [ok[i] ? vals[i] - lq[i] : -Inf for i in 1:S]
        # Off the held directions, where the normal stands as it is.
        GP = offhold(G * Gamma .+ (Z .- transpose(mu)))
        gamma = 1.0
        if _importance_ess(Phi) < essfrac * S
            lo, hi = 0.0, 1.0
            for _ in 1:60
                mid = (lo + hi) / 2
                _importance_ess(mid .* Phi) >= essfrac * S ? (lo = mid) : (hi = mid)
            end
            gamma = max(lo, 1e-12)
        end
        w = zeros(S)
        munew = mu
        Gnew = Gamma
        gbar = zeros(k)
        zbar = zeros(k)
        for attempt in 1:40
            lw = gamma .* Phi
            m = maximum(lw)
            w .= exp.(lw .- m)
            w[.!ok] .= 0.0
            w ./= sum(w)
            gbar = vec(sum(w .* GP; dims=1))
            zbar = vec(sum(w .* Z; dims=1))
            munew = mu .+ gamma .* gbar
            Cv = transpose(w .* (GP .- transpose(gbar))) * offhold(Z .- transpose(zbar))
            Gtry = Gamma .+ gamma .* (Cv .+ transpose(Cv)) ./ 2
            Gtry = (Gtry .+ transpose(Gtry)) ./ 2
            if isposdef(Symmetric(Gtry))
                Gnew = Gtry
                break
            end
            gamma /= 2
            if attempt == 40
                Gnew = Gamma
                munew = mu
            end
        end
        essw = _importance_ess(gamma .* Phi)
        shift = norm(Lc \ (munew .- mu))
        natural = norm(Lc \ gbar)
        W = Lc \ transpose(GP .- transpose(gbar))
        se2 = sum(w .^ 2 .* vec(sum(abs2, W; dims=1)))
        noise = 2 * sqrt(se2)
        lpf = vals[ok]
        elbo = sum(lpf) / length(lpf) + sum(log, diag(Lc))
        # Each parameter's step noise, in the sds of the normal the draws came
        # from: the mean step's weighted spread, and the variance step's.
        sig2 = vec(sum(abs2, A * Lc; dims=2))
        sig2 .= max.(sig2, eps() * maximum(sig2))
        Gp = (GP .- transpose(gbar)) * transpose(A)
        Zp = offhold(Z .- transpose(zbar)) * transpose(A)
        w2 = w .^ 2
        push!(se2ps, vec(sum(w2 .* Gp .^ 2; dims=1)) ./ sig2)
        push!(se2vs, vec(sum(w2 .* (Gp .* Zp) .^ 2; dims=1)) ./ sig2 .^ 2)
        push!(gammas, gamma)
        push!(esss, essw)
        mu = munew
        Gamma = Gnew
        push!(mus, copy(mu))
        push!(varsv, vec(sum((A * Gamma) .* A; dims=2)))
        Gs[mod1(iter, cap)] = copy(Gamma)
        win = _dais_window(esss, gammas, se2ps, se2vs, need, zp, P.mean_sd, P.var)
        m = win.m
        push!(trace, [gamma, essw, natural, noise, elbo, shift, 0.0, Float64(nok),
            Float64(m)])
        # A variance still growing past `maxgrowth` times the start's is a
        # direction with no fixed point: the run ends and the caller keeps the
        # check's normal.
        growth = maximum(eigvals(Symmetric(Gamma)))
        if growth > maxgrowth
            verbose && _due(reporter, true) && _emit(reporter, @sprintf(
                "  DAIS: the covariance grew %.0f times the normal's after %d iterations; the check's normal stands",
                growth, iter))
            _progress_done(reporter)
            return _dais_stopped(:diverged, growth, unusable, trace, S, gradients, cost,
                coreset !== nothing)
        end
        if _due(reporter)
            _emit(reporter, @sprintf("  DAIS %d: damping %.3g, ESS %.0f of %d, step %.2f (noise %.2f), ELBO %.2f%s",
                iter, gamma, essw, S, natural, noise, elbo,
                m > 0 ? ", window $m" : ""))
        end
        # The window's tests: no parameter's mean or variance has drifted over
        # it beyond `zp` times the random walk of its steps' noise, and the
        # last natural gradient is within twice its own.
        if m > 0 && isfinite(natural) && isfinite(noise) && natural <= noise
            sd = sqrt.(max.(varsv[end], eps() * maximum(varsv[end])))
            drift = abs.(A * (mus[end] .- mus[end - m])) ./ sd
            vdrift = abs.(varsv[end] .- varsv[end - m]) ./ sd .^ 2
            if all(drift .<= zp .* win.walkp .+ 1e-12) &&
                all(vdrift .<= zp .* win.walkv .+ 1e-12)
                window = m
                break
            end
        end
    end
    # The answer: the window's average where it stopped, the last half's where
    # the budget ended first.
    n = length(trace)
    settled = window > 0
    nav = min(settled ? window : max(1, n ÷ 2), cap, n)
    mu = sum(mus[(n - nav + 2):(n + 1)]) ./ nav
    Gamma = sum(Gs[mod1(t, cap)] for t in (n - nav + 1):n) ./ nav
    for t in (n - nav + 1):n
        trace[t][7] = 1.0
    end
    Gamma = (Gamma .+ transpose(Gamma)) ./ 2
    covraw = A * Gamma * transpose(A)
    _progress_done(reporter, @sprintf("%d iterations, %d averaged, %d gradients, %.0f full evaluations",
        n, nav, gradients, cost))
    return (used=true, diverged=false, growth=NaN, mean=c0 .+ A * mu,
        cov=(covraw .+ transpose(covraw)) ./ 2, iterations=n, averaged=nav,
        settled=settled, gamma=gamma, ess=essw, unusable=unusable, draws=S,
        gradients=gradients, cost=cost, _dais_coreset_record(coreset)...,
        parallel=String(mode), trace=_dais_trace(trace), trace_names=_DAIS_TRACE)
end

"""What DAIS's coreset did: whether there was one, its last size (0 where
none held), the coresets it chose and the checks that failed, and its whole
record as the importance sampler's coreset reports it (`ctsem_is_coreset_record`:
each choice's size and price among it)."""
_dais_coreset_record(::Nothing) = (coreset=false, coreset_size=0, coreset_chosen=0,
    coreset_failed=0)
_dais_coreset_record(cs::CTSEMISCoreset) = (coreset=true,
    coreset_size=cs.phase == :coreset ? _is_active_n(cs) : 0, coreset_chosen=cs.nchoose,
    coreset_failed=cs.fails_total, coreset_record=ctsem_is_coreset_record(cs))

"""The held directions `hold` (columns, in the parameters) in the whitened
coordinates of the kept eigenvectors `V` and root eigenvalues `sq`,
orthonormal: `A^+ h` for each, `A^+ = D^-1/2 V'`; `nothing` without any."""
function _dais_held_basis(V::AbstractMatrix, sq::AbstractVector, hold)
    (hold === nothing || size(hold, 2) == 0) && return nothing
    f = svd((transpose(V) * Matrix{Float64}(hold)) ./ sq)
    r = count(>(1e-7 * maximum(f.S; init=0.0)), f.S)
    return r > 0 ? f.U[:, 1:r] : nothing
end

"""A run that ended without a normal: diverged (with the growth that ended
it), with no usable draws, or with no direction to move. The check's normal
stands."""
_dais_stopped(why::Symbol, growth::Real, unusable::Int, trace, S::Int, gradients::Int,
    cost::Real, coreset::Bool) =
    (used=false, diverged=why == :diverged, growth=Float64(growth), unusable=unusable,
        iterations=length(trace), cov=nothing, draws=S, gradients=gradients,
        cost=Float64(cost), coreset=coreset, coreset_size=0, coreset_chosen=0,
        coreset_failed=0, trace=_dais_trace(trace), trace_names=_DAIS_TRACE)

"""The trace's columns, one row an iteration: `averaging` marks the iterates
the answer averages, `window` the stopping rule's window at that iteration (0
before the run could show what the answer is held to)."""
const _DAIS_TRACE = ["gamma", "ess", "natural", "noise", "elbo", "shift", "averaging",
    "usable", "window"]

_dais_trace(rows::Vector{Vector{Float64}}) =
    isempty(rows) ? zeros(0, length(_DAIS_TRACE)) : Matrix(transpose(reduce(hcat, rows)))
