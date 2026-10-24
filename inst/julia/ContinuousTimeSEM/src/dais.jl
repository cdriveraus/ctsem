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
    out right while its eigenstructure is noise. So once settled, `average`
    more iterations are run and their iterates averaged (Polyak), and the
    average is the answer; the last `average` are averaged where the budget
    ends first.
  * When to stop. Four rules failed on the time series; what holds is three
    tests together over a window of `window` iterations, `settle` times
    running: the natural gradient `E_w[Gamma grad Phi]` within twice its own
    weighted standard error; the ELBO's last window within two standard
    errors of the one before; and the mean's displacement over the window
    within twice the random walk its damped steps' noise would make. A test
    that cannot be read (a batch with no spread) is unmet.
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

Moved here from R (2026-10-24). In R each iteration crossed the bridge for
its batch and did its arithmetic and the coreset's refits in R: on the
266-parameter scale fit the R process used nine hours of CPU against the
engine's four, and `'checked'` had not finished in nine and a half. Here the
draws are evaluated at once on the pool or on copies (`_each_point!`), and
nothing crosses until the answer does.
"""

using LinearAlgebra
using Random

export ctsem_dais

"""Draws an iteration: the dimension, at least 64 and at most 256, even."""
_dais_draws(k::Int) = 2 * ceil(Int, min(256, max(64, k)) / 2)

"""The effective sample size of log weights (0 where none is finite)."""
function _dais_ess(lw::AbstractVector{Float64})
    m = maximum(lw)
    isfinite(m) || return 0.0
    w = exp.(lw .- m)
    return sum(w)^2 / sum(abs2, w)
end

"""
    ctsem_dais(objective, centre, cov; hold, draws, maxiter, seed, rtol,
        maxgrowth, parallel, verbose, progress_sink, progress_overwrite)

Move the normal `N(centre, cov)` to the posterior's moments by doubly
adaptive importance sampling (see the file comment). `hold`, a matrix whose
columns are directions in the parameters, holds them at the normal's width.
`draws` (0: `_dais_draws` of the dimension) and `maxiter` (150) bound the
run. Returns `used`, the settled normal (`mean`, `cov`), and the record:
iterations, how many were averaged, whether it settled, the last damping and
effective sample, the unusable draws, the draws an iteration, the gradients
spent, how the draws were evaluated (`parallel`: `_points_mode`'s serial, pool
or copies), and the trace, one row an iteration (`trace_names`); `used = false`
with `diverged` (and `growth`) or with no usable draws, and `cov = nothing`.
"""
function ctsem_dais(objective, centre::AbstractVector, cov::AbstractMatrix;
    hold::Union{Nothing,AbstractMatrix}=nothing, draws::Integer=0, maxiter::Integer=150,
    seed::Integer=20261024, rtol::Real=1e-10, maxgrowth::Real=1e4, parallel::Bool=true,
    verbose::Bool=false, progress_sink=nothing, progress_overwrite::Bool=true)

    c0 = collect(Float64, centre)
    p = length(c0)
    C = Matrix{Float64}((cov .+ transpose(cov)) ./ 2)
    size(C) == (p, p) || throw(DimensionMismatch("cov must be $p x $p"))
    e = eigen(Symmetric(C))
    top = maximum(e.values)
    keep = findall(v -> v > rtol * top, e.values)
    k = length(keep)
    k > 0 || return (used=false, diverged=false, unusable=0, iterations=0, cov=nothing,
        draws=0, gradients=0, trace=zeros(0, 8), trace_names=_DAIS_TRACE)
    V = e.vectors[:, keep]
    sq = sqrt.(e.values[keep])
    A = V .* transpose(sq)
    # The held directions in the whitened coordinates, orthonormal: `A^+ h`
    # for each, `A^+ = D^-1/2 V'` on the kept eigenvectors.
    U = _dais_held_basis(V, sq, hold)
    offhold(M) = U === nothing ? M : M .- (M * U) * transpose(U)
    S = draws > 0 ? max(4, 2 * ((Int(draws) + 1) ÷ 2)) : _dais_draws(k)
    essfrac = 0.5
    settle = 3
    window = 20
    average = 20
    mode = _points_mode(objective, parallel)
    rng = Random.Xoshiro(UInt64(seed))
    reporter = CTSEMProgress(verbose; label="DAIS", overwrite=progress_overwrite,
        sink=progress_sink)

    mu = zeros(k)
    Gamma = Matrix{Float64}(I, k, k)
    elbos = Float64[]
    mus = Vector{Vector{Float64}}()
    stepnoise = Float64[]
    recent = Vector{Tuple{Vector{Float64},Matrix{Float64}}}()
    averaging = 0
    settled = 0
    musum = zeros(k)
    Gsum = zeros(k, k)
    gradients = 0
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
        fill!(vals, -Inf)
        fill!(grads, 0.0)
        Xs = X
        valsv = vals
        gradsm = grads
        _each_point!(objective, S, mode) do i, obj, _
            r = _ctsem_probe_value_gradient(obj, Xs[:, i])
            if r.gradient !== nothing && isfinite(r.value) && all(isfinite, r.gradient)
                valsv[i] = r.value
                gradsm[:, i] .= r.gradient
            end
            return nothing
        end
        gradients += S
        ok = isfinite.(vals)
        nok = count(ok)
        unusable += S - nok
        if nok < 2
            verbose && _due(reporter, true) && _emit(reporter,
                "  DAIS: $nok of $S draws usable at iteration $iter; the check's normal stands")
            _progress_done(reporter)
            return _dais_stopped(:unusable, NaN, unusable, trace, S, gradients)
        end
        G = transpose(grads) * A
        G[.!ok, :] .= 0.0
        Phi = [ok[i] ? vals[i] - lq[i] : -Inf for i in 1:S]
        # Off the held directions, where the normal stands as it is.
        GP = offhold(G * Gamma .+ (Z .- transpose(mu)))
        gamma = 1.0
        if _dais_ess(Phi) < essfrac * S
            lo, hi = 0.0, 1.0
            for _ in 1:60
                mid = (lo + hi) / 2
                _dais_ess(mid .* Phi) >= essfrac * S ? (lo = mid) : (hi = mid)
            end
            gamma = max(lo, 1e-12)
        end
        w = zeros(S)
        munew = mu
        Gnew = Gamma
        gbar = zeros(k)
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
        essw = _dais_ess(gamma .* Phi)
        shift = norm(Lc \ (munew .- mu))
        natural = norm(Lc \ gbar)
        W = Lc \ transpose(GP .- transpose(gbar))
        se2 = sum(w .^ 2 .* vec(sum(abs2, W; dims=1)))
        noise = 2 * sqrt(se2)
        lpf = vals[ok]
        elbo = sum(lpf) / length(lpf) + sum(log, diag(Lc))
        push!(elbos, elbo)
        elbo_se = length(lpf) > 1 ? sqrt(sum(abs2, lpf .- sum(lpf) / length(lpf)) /
            (length(lpf) - 1)) / sqrt(length(lpf)) : NaN
        n_e = length(elbos)
        elbo_settled = n_e >= 2window &&
            (sum(elbos[(n_e - window + 1):n_e]) - sum(elbos[(n_e - 2window + 1):(n_e - window)])) /
            window <= 2 * elbo_se * sqrt(2 / window)
        push!(stepnoise, gamma * sqrt(se2))
        travel_settled = length(mus) >= window &&
            norm(Lc \ (munew .- mus[end - window + 1])) <=
            2 * sqrt(sum(abs2, stepnoise[max(1, end - window + 1):end]))
        mu = munew
        Gamma = Gnew
        push!(mus, copy(mu))
        push!(recent, (copy(mu), copy(Gamma)))
        length(recent) > average && popfirst!(recent)
        push!(trace, [gamma, essw, natural, noise, elbo, shift, averaging > 0 ? 1.0 : 0.0,
            Float64(nok)])
        # A variance still growing past `maxgrowth` times the start's is a
        # direction with no fixed point: the run ends and the caller keeps the
        # check's normal.
        growth = maximum(eigvals(Symmetric(Gamma)))
        if growth > maxgrowth
            verbose && _due(reporter, true) && _emit(reporter, @sprintf(
                "  DAIS: the covariance grew %.0f times the normal's after %d iterations; the check's normal stands",
                growth, iter))
            _progress_done(reporter)
            return _dais_stopped(:diverged, growth, unusable, trace, S, gradients)
        end
        if _due(reporter)
            _emit(reporter, @sprintf("  DAIS %d: damping %.3g, ESS %.0f of %d, step %.2f (noise %.2f), ELBO %.2f%s",
                iter, gamma, essw, S, natural, noise, elbo, averaging > 0 ? ", averaging" : ""))
        end
        # `averaging` counts the iterates summed and is the divisor.
        if averaging == 0
            settled = (isfinite(natural) && isfinite(noise) && natural <= noise &&
                elbo_settled === true && travel_settled === true) ? settled + 1 : 0
        end
        if averaging > 0 || settled >= settle
            averaging += 1
            musum .+= mu
            Gsum .+= Gamma
            averaging >= average && break
        end
    end
    if averaging > 0
        mu = musum ./ averaging
        Gamma = Gsum ./ averaging
    elseif !isempty(recent)
        mu = sum(first.(recent)) ./ length(recent)
        Gamma = sum(last.(recent)) ./ length(recent)
    end
    Gamma = (Gamma .+ transpose(Gamma)) ./ 2
    covraw = A * Gamma * transpose(A)
    _progress_done(reporter, @sprintf("%d iterations, %d averaged, %d gradients",
        length(trace), averaging, gradients))
    return (used=true, diverged=false, growth=NaN, mean=c0 .+ A * mu,
        cov=(covraw .+ transpose(covraw)) ./ 2, iterations=length(trace),
        averaged=averaging, settled=averaging > 0, gamma=gamma, ess=essw,
        unusable=unusable, draws=S, gradients=gradients, parallel=String(mode),
        trace=_dais_trace(trace), trace_names=_DAIS_TRACE)
end

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
it) or with no usable draws. The check's normal stands."""
_dais_stopped(why::Symbol, growth::Real, unusable::Int, trace, S::Int, gradients::Int) =
    (used=false, diverged=why == :diverged, growth=Float64(growth), unusable=unusable,
        iterations=length(trace), cov=nothing, draws=S, gradients=gradients,
        trace=_dais_trace(trace), trace_names=_DAIS_TRACE)

"""The trace's columns, one row an iteration."""
const _DAIS_TRACE = ["gamma", "ess", "natural", "noise", "elbo", "shift", "averaging",
    "usable"]

_dais_trace(rows::Vector{Vector{Float64}}) =
    isempty(rows) ? zeros(0, length(_DAIS_TRACE)) : Matrix(transpose(reduce(hcat, rows)))
