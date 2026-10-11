# The sequential predictions of a Laplace fit as the model's own predictive
# distribution: each row predicted from the subject's earlier rows with the
# random effects integrated against their posterior given those rows, by
# importance sampling, and the deviations the filter integrates (`CTSEMHybrid`)
# integrated exactly within each draw.
#
# `ctsem_laplace_sequential` re-solves the effects' mode on every truncated
# record and adds their uncertainty by the delta method: exact where the
# effects enter the predictions linearly, a little narrow and a little off
# centre where they reach the dynamics or the variances, and O(T^2) per
# subject. Its errors are then not quite martingale differences under the
# model, which a misfit test with the estimation removed sees: on data
# generated from a pilot-shaped model (random drift and diffusion) the row
# distances ran 1.5 per cent high and their mean rejected the correct model at
# z = 6. Here, for a fixed set of draws u_k from a proposal covering both the
# prior and the whole-record posterior, one filter pass per draw over the whole
# record gives every row's prediction given u_k, and weighting the draws by the
# likelihood of the rows so far gives the effects' posterior at every row at
# once: O(K T) per subject. Two passes per draw -- the first for the weights,
# the second to accumulate the weighted moments -- so nothing per draw is kept
# but its row log densities.
#
# The draws are fixed per unit (seeded) and placed by the whole-record mode and
# curvature, so the result is a smooth function of the parameters, as a
# derivative of statistics built on it needs.

"""
Per-row output of one filter pass at fixed random effects, with the
deviations the filter integrates (`CTSEMHybrid`) marginalised at every row from
the pass's own augmented quantities: given the rows so far the standardised
deviations are `N(M^-1 G s, M^-1)`, `M = I + G S G`, and they move the state by
`X G` and the observations by `(Jy X + E) G` (`hybrid_filter.jl`).
"""
mutable struct CTSEMMixtureTrace
    n::Int
    m::Int
    nrows::Int
    data::Matrix{Float64}
    values::Vector{Float64}
    hybrid::CTSEMHybrid
    g::Vector{Float64}
    etaP::Matrix{Float64}
    etaPcov::Array{Float64,3}
    etaU::Matrix{Float64}
    etaUcov::Array{Float64,3}
    yP::Matrix{Float64}
    yPcov::Array{Float64,3}
    yU::Matrix{Float64}
    yUcov::Array{Float64,3}
    llgauss::Vector{Float64}
    llfilter::Vector{Float64}
    # Each non-Gaussian observation's predictive distribution function bounds
    # and log predictive probability under this draw (item_pit.jl).
    pitlower::Matrix{Float64}
    pitupper::Matrix{Float64}
    pitlogz::Matrix{Float64}
    pitwanted::Bool
    # scratch
    Mq::Matrix{Float64}
    Cq::Matrix{Float64}
    bq::Vector{Float64}
    dq::Vector{Float64}
    eq::Vector{Float64}
    XG::Matrix{Float64}
    VG::Matrix{Float64}
    Soo::Matrix{Float64}
    r::Vector{Float64}
    idx::Vector{Int}
    AC::Matrix{Float64}
    JP::Matrix{Float64}
end

function CTSEMMixtureTrace(n::Int, m::Int, nrows::Int, hybrid::CTSEMHybrid)
    CTSEMMixtureTrace(n, m, nrows, zeros(m, nrows), Float64[], hybrid, Float64[],
        zeros(nrows, n), zeros(nrows, n, n), zeros(nrows, n), zeros(nrows, n, n),
        zeros(nrows, m), zeros(nrows, m, m), zeros(nrows, m), zeros(nrows, m, m),
        zeros(nrows), zeros(nrows),
        fill(NaN, nrows, m), fill(NaN, nrows, m), fill(NaN, nrows, m), false,
        zeros(0, 0), zeros(0, 0), Float64[], Float64[], Float64[],
        zeros(n, 0), zeros(m, 0), zeros(m, m), zeros(m), zeros(Int, m),
        zeros(max(n, m), 0), zeros(m, n))
end

# Every hook this trace does not want, declined as the Kalman trace declines them.
@inline _record_subject_values!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_init!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_theta!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_group!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_td!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_update!(::CTSEMMixtureTrace, args...) = nothing
@inline _begin_predict!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_predict!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_stationary!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_transition!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_td_transition!(::CTSEMMixtureTrace, args...) = nothing

# A row's non-Gaussian observations: each one's distribution function bounds
# and log predictive probability, as the update applies it.
function _record_binary!(tr::CTSEMMixtureTrace, ws, pars, data, obs_col, rows, n)
    (tr.pitwanted && !isempty(rows)) || return nothing
    return CTSEMPITRecord(tr.pitlower, tr.pitupper, tr.pitlogz, Int(obs_col),
        collect(Int, rows), 0)
end

"""Each hybrid column's scale (an item's own sd times its intercept's slope),
as `_hybrid_marginal!` computes it, and the scratch sized to `q`."""
function _mixture_hybrid_scale!(tr::CTSEMMixtureTrace, hb)
    q = hb.q
    h = tr.hybrid
    length(tr.g) == q || (tr.g = zeros(q))
    @inbounds for c in 1:q
        k = hb.item[c]
        raw = tr.values[h.psi_index[k]]
        tr.g[c] = h.slope[k] * (log1p_exp(2 * raw - 1) * h.sd_scale[k] + 1e-10)
    end
    if size(tr.Mq, 1) != q
        tr.Mq = zeros(q, q); tr.Cq = zeros(q, q)
        tr.bq = zeros(q); tr.dq = zeros(q); tr.eq = zeros(q)
        tr.XG = zeros(tr.n, q); tr.VG = zeros(tr.m, q)
        tr.AC = zeros(max(tr.n, tr.m), q)
    end
    return nothing
end

"""The deviations' posterior given the rows the buffers have seen: mean into
`tr.dq`, covariance into `tr.Cq`; `false` where `M` does not factor."""
function _mixture_hybrid_posterior!(tr::CTSEMMixtureTrace, hb)
    q = hb.q
    g = tr.g
    @inbounds for b in 1:q, a in 1:q
        tr.Mq[a, b] = g[a] * hb.S[a, b] * g[b] + (a == b ? 1.0 : 0.0)
    end
    F = _ctsem_cholesky(tr.Mq, q)
    issuccess(F) || return false
    @inbounds for c in 1:q
        tr.bq[c] = g[c] * hb.s[c]
    end
    ldiv!(tr.dq, F, tr.bq)
    @inbounds for c in 1:q
        fill!(tr.eq, 0.0)
        tr.eq[c] = 1.0
        ldiv!(view(tr.Cq, :, c), F, tr.eq)
    end
    return true
end

"""`mean[row, :] += A d` and `cov[row, :, :] += A C A'`, `A` the first `k` rows
of a `k x q` loading, through `AC = A C` (`k q^2 + k^2 q`, not `k^2 q^2`)."""
function _mixture_add!(mean, cov, A, d, C, k::Int, q::Int, row::Int, AC)
    @inbounds for i in 1:k
        acc = 0.0
        for c in 1:q
            acc += A[i, c] * d[c]
        end
        mean[row, i] += acc
    end
    @inbounds for b in 1:q, i in 1:k
        acc = 0.0
        for a in 1:q
            acc += A[i, a] * C[a, b]
        end
        AC[i, b] = acc
    end
    @inbounds for j in 1:k, i in 1:k
        acc = 0.0
        for b in 1:q
            acc += AC[i, b] * A[j, b]
        end
        cov[row, i, j] += acc
    end
    return nothing
end

"""The state's and the observations' loadings on the deviations, `X G` and
`(Jy X + E) G`, from the buffers as they stand."""
function _mixture_loadings!(tr::CTSEMMixtureTrace, hb, pars)
    n, m, q = tr.n, tr.m, hb.q
    @inbounds for c in 1:q
        for r in 1:n
            tr.XG[r, c] = hb.X[r, c] * tr.g[c]
        end
        for j in 1:m
            acc = hb.col[j] == c ? 1.0 : 0.0
            for k in 1:n
                acc += pars.Jy[j, k] * hb.X[k, c]
            end
            tr.VG[j, c] = acc * tr.g[c]
        end
    end
    return nothing
end

"""The measurement model at `state` and `statecov` into `y[row, :]` and
`ycov[row, :, :]`, as `_kalman_store!` computes it, through `JP = Jy P`."""
function _mixture_measure!(y, ycov, ws, pars, state, statecov, n::Int, m::Int, row::Int, JP)
    Θ = ws.bufferΘ.out
    @inbounds for i in 1:m
        value = pars.MANIFESTMEANS[i, 1]
        for j in 1:n
            value += pars.LAMBDA[i, j] * state[j]
        end
        y[row, i] = value
    end
    @inbounds for b in 1:n, i in 1:m
        acc = 0.0
        for a in 1:n
            acc += pars.Jy[i, a] * statecov[a, b]
        end
        JP[i, b] = acc
    end
    @inbounds for j in 1:m, i in 1:m
        value = Θ[i, j]
        for b in 1:n
            value += JP[i, b] * pars.Jy[j, b]
        end
        ycov[row, i, j] = value
    end
    return nothing
end

function _record_row_prior!(tr::CTSEMMixtureTrace, ws, pars, row::Int)
    n, m = tr.n, tr.m
    state = ws.state
    P = ws.P_predict.data
    @inbounds for j in 1:n, i in 1:n
        tr.etaPcov[row, i, j] = P[i, j]
    end
    @inbounds for i in 1:n
        tr.etaP[row, i] = state[i]
    end
    _mixture_measure!(tr.yP, tr.yPcov, ws, pars, state, P, n, m, row, tr.JP)
    hb = ws.hyb
    if hb.active && hb.q > 0
        row == 1 && _mixture_hybrid_scale!(tr, hb)
        if !_mixture_hybrid_posterior!(tr, hb)
            tr.llgauss[row] = NaN
            return nothing
        end
        _mixture_loadings!(tr, hb, pars)
        _mixture_add!(tr.etaP, tr.etaPcov, tr.XG, tr.dq, tr.Cq, n, hb.q, row, tr.AC)
        _mixture_add!(tr.yP, tr.yPcov, tr.VG, tr.dq, tr.Cq, m, hb.q, row, tr.AC)
    end
    tr.llgauss[row] = _mixture_row_density(tr, row)
    return nothing
end

function _record_row_update!(tr::CTSEMMixtureTrace, ws, pars, sp, context, row::Int, row_ll)
    n, m = tr.n, tr.m
    state = ws.state
    P = ws.P_update.data
    @inbounds for j in 1:n, i in 1:n
        tr.etaUcov[row, i, j] = P[i, j]
    end
    @inbounds for i in 1:n
        tr.etaU[row, i] = state[i]
    end
    _mixture_measure!(tr.yU, tr.yUcov, ws, pars, state, P, n, m, row, tr.JP)
    hb = ws.hyb
    if hb.active && hb.q > 0 && _mixture_hybrid_posterior!(tr, hb)
        _mixture_loadings!(tr, hb, pars)
        _mixture_add!(tr.etaU, tr.etaUcov, tr.XG, tr.dq, tr.Cq, n, hb.q, row, tr.AC)
        _mixture_add!(tr.yU, tr.yUcov, tr.VG, tr.dq, tr.Cq, m, hb.q, row, tr.AC)
    end
    tr.llfilter[row] = row_ll === nothing ? NaN : Float64(row_ll)
    return nothing
end

"""The observed items' Gaussian log density at `row`, from the predicted mean
and covariance the trace holds."""
function _mixture_row_density(tr::CTSEMMixtureTrace, row::Int)
    k = 0
    @inbounds for j in 1:tr.m
        _ctsem_observed(tr.data[j, row]) || continue
        k += 1
        tr.idx[k] = j
        tr.r[k] = tr.data[j, row] - tr.yP[row, j]
    end
    k == 0 && return 0.0
    @inbounds for b in 1:k, a in 1:k
        tr.Soo[a, b] = tr.yPcov[row, tr.idx[a], tr.idx[b]]
    end
    F = _ctsem_cholesky(tr.Soo, k)
    issuccess(F) || return NaN
    U = F.U
    logdet = 0.0
    quad = 0.0
    # w = U^-T r by forward substitution on U', in place; quad = w'w.
    @inbounds for i in 1:k
        acc = tr.r[i]
        for j in 1:(i - 1)
            acc -= U[j, i] * tr.r[j]
        end
        tr.r[i] = acc / U[i, i]
        quad += tr.r[i]^2
        logdet += 2 * log(U[i, i])
    end
    return -0.5 * (k * log(2pi) + logdet + quad)
end

"""The row log densities a pass gives the weights: the marginal over the
integrated deviations where there are any, the filter's own otherwise."""
_mixture_rowll(tr::CTSEMMixtureTrace, hb) = hb.active && hb.q > 0 ? tr.llgauss : tr.llfilter

"""The first `d` primes."""
function _mixture_primes(d::Int)
    out = Int[]
    c = 2
    while length(out) < d
        all(c % p != 0 for p in out if p * p <= c) && push!(out, c)
        c += 1
    end
    return out
end

function _logsumexp(x::AbstractVector{Float64})
    mx = maximum(x)
    isfinite(mx) || return mx
    s = 0.0
    @inbounds for v in x
        s += exp(v - mx)
    end
    return mx + log(s)
end

"""
    ctsem_laplace_sequential_mixture(laplace, values; particles=128, scale=1.15,
        seed=20261011)

The output of `ctsem_laplace_sequential` -- predicted (leading index 1) and
updated (2) states and observations with their covariances, and each row's log
predictive density -- with the random effects integrated against their
posterior given the subject's earlier rows by importance sampling rather than
by the mode of each truncated record and the delta method. `particles` draws a
unit, quasi-random, shared equally between the prior and the Laplace posteriors
given the record cut at rows T/8, T/4, T/2 and T (covariances scaled by
`scale^2`), so some component sits near every row's posterior, and weighted by
their mixture (deterministic multiple importance sampling, the balance
heuristic). `ess` is each unit's smallest effective sample size over its rows.
Single-level fits. With non-Gaussian indicators, each one's predictive
distribution function just below and at what was observed, in the filter's
order (`pitlower`, `pitupper`), for the probability integral transforms the
misfit tests use in place of their errors (item_pit.jl).
"""
function ctsem_laplace_sequential_mixture(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; particles::Integer=128, scale::Real=1.15,
    seed::Integer=20261011, compact::Bool=false)
    theta = collect(Float64, values)
    units = laplace.units
    all(length(m) == 1 for m in units.members) || throw(ArgumentError(
        "sequential prediction needs every subject to be its own unit: a single " *
        "level of random effects"))
    particles >= 2 || throw(ArgumentError("need at least two particles"))
    sp = laplace.objective.params
    ws0 = _init_continuous_ekf_workspace(Float64, sp)
    n = _val(ws0.state_dim)
    m = _val(ws0.manifest_dim)
    types = ws0.manifesttype
    Ls = _laplace_popchols(theta, laplace.spec)
    subjects = laplace.objective.subject_objectives
    rows = [size(s.data, 2) for s in subjects]
    first_row = cumsum(vcat(0, rows[1:(end - 1)]))
    nrows = sum(rows)
    out = (eta=fill(NaN, 2, nrows, n), etacov=fill(NaN, 2, nrows, n, n),
        y=fill(NaN, 2, nrows, m), ycov=fill(NaN, 2, nrows, m, m),
        llrow=fill(NaN, nrows), subject=zeros(Int, nrows),
        ess=fill(NaN, length(units.members)),
        pitlower=fill(NaN, nrows, m), pitupper=fill(NaN, nrows, m))
    for (i, r) in enumerate(rows)
        out.subject[(first_row[i] + 1):(first_row[i] + r)] .= i
    end
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, eachindex(units.members)) do U
        _laplace_mixture_unit!(out, laplace, U, theta, Ls,
            first_row[units.members[U][1]], Int(particles), Float64(scale), Int(seed))
        return true
    end
    return compact ? _ctsem_kalman_compact(out) : out
end

function _laplace_mixture_unit!(out, laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, offset::Integer,
    particles::Int, scale::Float64, seed::Int)
    spec = laplace.spec
    units = laplace.units
    i = units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    sp = laplace.objective.params
    T = size(so.data, 2)
    d = units.dims[U]
    # The proposal's components: the prior, and the Laplace posterior given the
    # record cut at rows T/8, T/4, T/2 and T, each solved as the objective
    # solves it, so that one sits near every row's posterior. The unit's record
    # and retained state are put back afterwards.
    saved = (mode=copy(laplace.modes[U]), iterations=laplace.inner_iterations[U],
        gradient=laplace.inner_gradient[U], converged=laplace.inner_converged[U],
        hrepaired=laplace.hessian_repaired[U], mrepaired=laplace.mode_repaired[U],
        floored=laplace.logdet_floored[U])
    full = (data=so.data, timesteps=so.timesteps, tdpreds=so.tdpreds,
        max_timestep=so.max_timestep)
    comps = [(mu=zeros(d), vecs=Matrix{Float64}(I, d, d), lam=ones(d), c=1.0)]
    if d > 0
        try
            for t in unique(clamp.(ceil.(Int, T .* (0.125, 0.25, 0.5, 1.0)), 1, T))
                _laplace_truncate!(so, full, t)
                _laplace_solve_unit_mode!(laplace, U, theta, Ls)
                mu = copy(laplace.modes[U])
                Sig = _laplace_sequential_inverse(_laplace_block_dense(
                    _laplace_unit_curvature(laplace, U, theta, Ls, mu), units.blocks[U], d))
                E = _ctsem_symeig(Sig)
                push!(comps, (mu=mu, vecs=E.vectors, lam=max.(E.values, 1e-12), c=scale))
            end
        finally
            so.data = full.data
            so.timesteps = full.timesteps
            so.tdpreds = full.tdpreds
            so.max_timestep = full.max_timestep
            laplace.modes[U] = saved.mode
            laplace.inner_iterations[U] = saved.iterations
            laplace.inner_gradient[U] = saved.gradient
            laplace.inner_converged[U] = saved.converged
            laplace.hessian_repaired[U] = saved.hrepaired
            laplace.mode_repaired[U] = saved.mrepaired
            laplace.logdet_floored[U] = saved.floored
        end
    end

    # The draws: a randomly shifted lattice (Kronecker, frac(k sqrt(p_j) + s),
    # p_j the j-th prime) mapped to normals, the same number in each component,
    # fixed for the unit; each weighted by the prior over the components'
    # mixture. Quasi-random points, because the integrands are smooth: at the
    # same number of points they put far less error into each row's moments
    # than pseudo-random draws, whose relative error in the predictive variance
    # was 3 per cent at 512 draws a unit.
    J = length(comps)
    perj = d == 0 ? 0 : max(2, fld(particles, J))
    K = d == 0 ? 1 : perj * J
    rng = Random.Xoshiro(hash((seed, Int(U))))
    draws = Vector{Vector{Float64}}(undef, K)
    logw0 = zeros(K)
    if d == 0
        draws[1] = Float64[]
    else
        alpha = [mod(sqrt(Float64(p)), 1.0) for p in _mixture_primes(d)]
        k = 0
        for cmp in comps
            shift = rand(rng, d)
            for j in 1:perj
                z = [sqrt(2.0) * erfinv(clamp(2 * mod(j * alpha[a] + shift[a], 1.0) - 1,
                    -1 + 1e-15, 1 - 1e-15)) for a in 1:d]
                k += 1
                draws[k] = cmp.mu .+ cmp.c .* (cmp.vecs * (sqrt.(cmp.lam) .* z))
            end
        end
        logdens = function (cmp, u)
            w = cmp.vecs' * (u .- cmp.mu)
            return -0.5 * (sum(w .^ 2 ./ (cmp.c^2 .* cmp.lam)) + d * log(2pi)) -
                sum(log.(cmp.lam)) / 2 - d * log(cmp.c)
        end
        for k in 1:K
            u = draws[k]
            logprior = -0.5 * (dot(u, u) + d * log(2pi))
            logq = _logsumexp([logdens(cmp, u) for cmp in comps]) - log(J)
            logw0[k] = logprior - logq
        end
    end

    n = size(out.eta, 3)
    m = size(out.y, 3)
    ws = _init_continuous_ekf_workspace(Float64, sp)
    ws.hyb.mode = _HYBRID_FULL
    tr = CTSEMMixtureTrace(n, m, T, sp.hybrid)
    tr.data = map(v -> ismissing(v) ? NaN : Float64(v), so.data)
    timesteps = collect(so.timesteps)
    pass! = function (k)
        x = _laplace_member_values(theta, spec, Ls, draws[k], units.offsets[U][1])
        tr.values = collect(Float64, x)
        ll = _extended_kalman_filter_continuous!(ws, x, so.data, timesteps, sp,
            so.tdpreds, _ctsem_tipred_vector(so.tipreds, x), i, so.max_timestep, tr)
        return isfinite(ll)
    end

    # First pass: each draw's row log densities, for the weights, and with
    # non-Gaussian indicators each observation's distribution function bounds
    # and log probability under the draw.
    types = ws.manifesttype
    cats = [j for j in 1:m if _ekf_is_categorical(types, j)]
    nc = length(cats)
    PL = nc > 0 ? fill(NaN, K, T, nc) : zeros(0, 0, 0)
    PH = similar(PL); PZ = similar(PL)
    LL = fill(-Inf, K, T)
    ok = falses(K)
    tr.pitwanted = nc > 0
    for k in 1:K
        if nc > 0
            fill!(tr.pitlower, NaN); fill!(tr.pitupper, NaN); fill!(tr.pitlogz, NaN)
        end
        ok[k] = pass!(k)
        ok[k] || continue
        rl = _mixture_rowll(tr, ws.hyb)
        @inbounds for t in 1:T
            LL[k, t] = isfinite(rl[t]) ? rl[t] : -Inf
            for (c, j) in enumerate(cats)
                PL[k, t, c] = tr.pitlower[t, j]
                PH[k, t, c] = tr.pitupper[t, j]
                PZ[k, t, c] = tr.pitlogz[t, j]
            end
        end
    end
    tr.pitwanted = false
    any(ok) || return nothing
    # Weights before each row (its prediction) and after it (its update).
    WP = zeros(K, T)
    WU = zeros(K, T)
    cum = copy(logw0)
    cum[.!ok] .= -Inf
    minESS = Inf
    @inbounds for t in 1:T
        a = _logsumexp(cum)
        for k in 1:K
            WP[k, t] = exp(cum[k] - a)
        end
        # The row's log predictive density under the mixture.
        out.llrow[offset + t] = _logsumexp(cum .+ LL[:, t]) - a
        for k in 1:K
            cum[k] += LL[k, t]
        end
        b = _logsumexp(cum)
        for k in 1:K
            WU[k, t] = exp(cum[k] - b)
        end
        minESS = min(minESS, 1 / sum(abs2, view(WP, :, t)))
    end
    out.ess[U] = minESS

    # Each non-Gaussian observation's predictive distribution function under
    # the mixture, in the order the filter applies them: the draws weighted by
    # the row's prediction weights times each draw's probability of the row's
    # earlier non-Gaussian observations, so that it conditions on what the
    # filter had applied before it (item_pit.jl).
    if nc > 0
        logw = zeros(K)
        @inbounds for t in 1:T
            for k in 1:K
                logw[k] = ok[k] && WP[k, t] > 0 ? log(WP[k, t]) : -Inf
            end
            for c in 1:nc
                any(k -> ok[k] && isfinite(PL[k, t, c]), 1:K) || continue
                a = _logsumexp(logw)
                lo = 0.0; hi = 0.0
                for k in 1:K
                    isfinite(logw[k]) || continue
                    w = exp(logw[k] - a)
                    lo += w * PL[k, t, c]
                    hi += w * PH[k, t, c]
                end
                out.pitlower[offset + t, cats[c]] = lo
                out.pitupper[offset + t, cats[c]] = hi
                for k in 1:K
                    logw[k] += isfinite(PZ[k, t, c]) ? PZ[k, t, c] : -Inf
                end
            end
        end
    end

    # Second pass: the weighted moments, with the spread between draws.
    sP = zeros(T, n); sPP = zeros(T, n, n); sU = zeros(T, n); sUU = zeros(T, n, n)
    yP = zeros(T, m); yPP = zeros(T, m, m); yU = zeros(T, m); yUU = zeros(T, m, m)
    for k in 1:K
        ok[k] || continue
        pass!(k) || continue
        @inbounds for t in 1:T
            wp = WP[k, t]
            wu = WU[k, t]
            for a in 1:n
                sP[t, a] += wp * tr.etaP[t, a]
                sU[t, a] += wu * tr.etaU[t, a]
                for b in 1:n
                    sPP[t, a, b] += wp * (tr.etaPcov[t, a, b] + tr.etaP[t, a] * tr.etaP[t, b])
                    sUU[t, a, b] += wu * (tr.etaUcov[t, a, b] + tr.etaU[t, a] * tr.etaU[t, b])
                end
            end
            for a in 1:m
                yP[t, a] += wp * tr.yP[t, a]
                yU[t, a] += wu * tr.yU[t, a]
                for b in 1:m
                    yPP[t, a, b] += wp * (tr.yPcov[t, a, b] + tr.yP[t, a] * tr.yP[t, b])
                    yUU[t, a, b] += wu * (tr.yUcov[t, a, b] + tr.yU[t, a] * tr.yU[t, b])
                end
            end
        end
    end
    @inbounds for t in 1:T
        r = offset + t
        for a in 1:n
            out.eta[1, r, a] = sP[t, a]
            out.eta[2, r, a] = sU[t, a]
            for b in 1:n
                out.etacov[1, r, a, b] = sPP[t, a, b] - sP[t, a] * sP[t, b]
                out.etacov[2, r, a, b] = sUU[t, a, b] - sU[t, a] * sU[t, b]
            end
        end
        for a in 1:m
            out.y[1, r, a] = yP[t, a]
            out.y[2, r, a] = yU[t, a]
            for b in 1:m
                out.ycov[1, r, a, b] = yPP[t, a, b] - yP[t, a] * yP[t, b]
                out.ycov[2, r, a, b] = yUU[t, a, b] - yU[t, a] * yU[t, b]
            end
        end
    end
    return nothing
end

export ctsem_laplace_sequential_mixture
