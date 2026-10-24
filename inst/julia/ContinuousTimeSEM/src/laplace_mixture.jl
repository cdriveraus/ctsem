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
deviations are `N(M^-1 G' s, M^-1)`, `M = I + G' S G`, and they move the state
by `X G` and the observations by `(Jy X + E) G` (`hybrid_filter.jl`), `G` the
filter's own factor (`_hybrid_factor!`).
"""
mutable struct CTSEMMixtureTrace
    n::Int
    m::Int
    nrows::Int
    data::Matrix{Float64}
    values::Vector{Float64}
    hybrid::CTSEMHybrid
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
    CTSEMMixtureTrace(n, m, nrows, zeros(m, nrows), Float64[], hybrid,
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
@inline _record_stationary_elements!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_transition!(::CTSEMMixtureTrace, args...) = nothing
@inline _record_td_transition!(::CTSEMMixtureTrace, args...) = nothing

# A row's non-Gaussian observations: each one's distribution function bounds
# and log predictive probability, as the update applies it.
function _record_binary!(tr::CTSEMMixtureTrace, ws, pars, data, obs_col, rows, n)
    (tr.pitwanted && !isempty(rows)) || return nothing
    return CTSEMPITRecord(tr.pitlower, tr.pitupper, tr.pitlogz, Int(obs_col),
        collect(Int, rows), 0)
end

"""The columns' factor `G` into the buffers, as `_hybrid_marginal!` computes
it (`_hybrid_factor!`), and the scratch sized to its `rd` columns; false where
it does not factor."""
function _mixture_hybrid_scale!(tr::CTSEMMixtureTrace, hb)
    rd = _hybrid_factor!(hb, tr.hybrid, tr.values)
    rd > 0 || return false
    if size(tr.Mq, 1) != rd
        tr.Mq = zeros(rd, rd); tr.Cq = zeros(rd, rd)
        tr.bq = zeros(rd); tr.dq = zeros(rd); tr.eq = zeros(rd)
        tr.XG = zeros(tr.n, rd); tr.VG = zeros(tr.m, rd)
        tr.AC = zeros(max(tr.n, tr.m), rd)
    end
    return true
end

"""The deviations' posterior given the rows the buffers have seen: mean into
`tr.dq`, covariance into `tr.Cq`; `false` where `M` does not factor."""
function _mixture_hybrid_posterior!(tr::CTSEMMixtureTrace, hb)
    q, rd = hb.q, hb.rd
    # The filter's own closing (`_hybrid_closing!`), `Cq` its scratch until it
    # takes the covariance (`q <= rd`, so `rd x q` fits).
    _hybrid_closing!(tr.Mq, tr.bq, tr.Cq, hb.G, hb.S, hb.s, q, rd)
    F = _ctsem_cholesky(tr.Mq, rd)
    issuccess(F) || return false
    ldiv!(tr.dq, F, tr.bq)
    @inbounds for a in 1:rd
        fill!(tr.eq, 0.0)
        tr.eq[a] = 1.0
        ldiv!(view(tr.Cq, :, a), F, tr.eq)
    end
    return true
end

"""`mean[row, :] += A d` and `cov[row, :, :] += A C A'`, `A` the first `k` rows
of a `k x q` loading (`q` the deviations' count), through `AC = A C` (`k q^2 + k^2 q`, not `k^2 q^2`)."""
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
    n, m, q, rd, G = tr.n, tr.m, hb.q, hb.rd, hb.G
    fill!(tr.XG, 0.0)
    fill!(tr.VG, 0.0)
    @inbounds for c in 1:q
        for j in 1:m
            v = hb.col[j] == c ? 1.0 : 0.0
            for k in 1:n
                v += pars.Jy[j, k] * hb.X[k, c]
            end
            for a in 1:rd
                tr.VG[j, a] += v * G[c, a]
            end
        end
        for r in 1:n
            x = hb.X[r, c]
            for a in 1:rd
                tr.XG[r, a] += x * G[c, a]
            end
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
        if (row == 1 && !_mixture_hybrid_scale!(tr, hb)) || !_mixture_hybrid_posterior!(tr, hb)
            tr.llgauss[row] = NaN
            return nothing
        end
        _mixture_loadings!(tr, hb, pars)
        _mixture_add!(tr.etaP, tr.etaPcov, tr.XG, tr.dq, tr.Cq, n, hb.rd, row, tr.AC)
        _mixture_add!(tr.yP, tr.yPcov, tr.VG, tr.dq, tr.Cq, m, hb.rd, row, tr.AC)
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
        _mixture_add!(tr.etaU, tr.etaUcov, tr.XG, tr.dq, tr.Cq, n, hb.rd, row, tr.AC)
        _mixture_add!(tr.yU, tr.yUcov, tr.VG, tr.dq, tr.Cq, m, hb.rd, row, tr.AC)
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
heuristic). `ess` is each unit's smallest effective sample size over its rows;
`Inf` for a unit with no coordinate, whose one draw is exact (the filter
integrating every effect it has). Single-level fits. With non-Gaussian indicators, each one's predictive
distribution function just below and at what was observed, in the filter's
order (`pitlower`, `pitupper`), for the probability integral transforms the
misfit tests use in place of their errors (item_pit.jl).

For the Gaussian items, the mixture's own distribution function at what was
observed (`_mixture_gauss_pit!`): `pitgauss`, each item given the row's other
Gaussian items, and `pitseq`, each given those before it in item order. The
mixture is not normal, so an error standardised by its mean and covariance is
not either: its fourth moments and the row distance's spread are not the
normal's, and checks reading them rejected a correct model more often the more
people there were. Through the normal quantile these are standard normal under
the model, and with one component they are the standardised errors.
"""
function ctsem_laplace_sequential_mixture(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; particles::Integer=128, scale::Real=1.15,
    seed::Integer=20261011, compact::Bool=false)
    theta = _laplace_extend(laplace, collect(Float64, values))
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
        pitlower=fill(NaN, nrows, m), pitupper=fill(NaN, nrows, m),
        pitgauss=fill(NaN, nrows, m), pitseq=fill(NaN, nrows, m))
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
    saved = _laplace_unit_saved(laplace, U, so)
    full = saved.full
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
            _laplace_unit_restore!(laplace, U, so, saved)
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
    # Nothing sampled is nothing thin: read as 1, a unit the filter integrates
    # whole made ctFitAdequacy report its errors' sampling as thin.
    out.ess[U] = d == 0 ? Inf : minESS

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

    # Second pass: the weighted moments, with the spread between draws, and the
    # Gaussian items' distribution functions under the mixture.
    sP = zeros(T, n); sPP = zeros(T, n, n); sU = zeros(T, n); sUU = zeros(T, n, n)
    yP = zeros(T, m); yPP = zeros(T, m, m); yU = zeros(T, m); yUU = zeros(T, m, m)
    gauss = [j for j in 1:m if !_ekf_is_categorical(types, j)]
    pit = _CTSEMMixturePIT(T, m)
    for k in 1:K
        ok[k] || continue
        pass!(k) || continue
        isempty(gauss) || _mixture_gauss_pit!(pit, tr, WP, k, gauss)
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
    @inbounds for t in 1:T, j in gauss
        pit.gden[t, j] > 0 && (out.pitgauss[offset + t, j] = pit.gnum[t, j] / pit.gden[t, j])
        pit.sden[t, j] > 0 && (out.pitseq[offset + t, j] = pit.snum[t, j] / pit.sden[t, j])
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

# Running sums, over the draws, of each row's Gaussian items' conditional
# distribution functions at what was observed and of their weights, each pair
# on the scale of the largest log weight so far (log-sum-exp as it goes): `g`
# each item given the row's other Gaussian items, `s` given those before it.
struct _CTSEMMixturePIT
    gmx::Matrix{Float64}
    gnum::Matrix{Float64}
    gden::Matrix{Float64}
    smx::Matrix{Float64}
    snum::Matrix{Float64}
    sden::Matrix{Float64}
    U::Matrix{Float64}
    Ui::Matrix{Float64}
    r::Vector{Float64}
    obs::Vector{Int}
end

_CTSEMMixturePIT(T::Int, m::Int) = _CTSEMMixturePIT(fill(-Inf, T, m), zeros(T, m),
    zeros(T, m), fill(-Inf, T, m), zeros(T, m), zeros(T, m), zeros(m, m), zeros(m, m),
    zeros(m), zeros(Int, m))

@inline function _mixture_pit_add!(mx, num, den, t, j, lw, p)
    isfinite(lw) || return nothing
    @inbounds if lw > mx[t, j]
        s = exp(mx[t, j] - lw)
        num[t, j] = num[t, j] * s + p
        den[t, j] = den[t, j] * s + 1
        mx[t, j] = lw
    else
        e = exp(lw - mx[t, j])
        num[t, j] += e * p
        den[t, j] += e
    end
    return nothing
end

# Draw k's share of each row's Gaussian items' distribution functions under the
# mixture. Given the person's earlier rows, the draw predicts the row's
# observed Gaussian items as N(mu, S) (the trace's prior prediction, as the
# standardised errors use); with U'U = S, r = U'^-1 (y - mu) are the items
# standardised given those before them, and (S^-1 (y - mu))_j / sqrt(S^-1_jj)
# each item given the others. A draw's weight for an item is its prediction
# weight times its density of what the item is conditioned on: the items
# before it, or the others -- the row's density less the item's conditional
# one. Hand-written triangular arithmetic: these run per row and draw inside a
# threaded loop, where LAPACK's lock would serialise them.
function _mixture_gauss_pit!(pit::_CTSEMMixturePIT, tr::CTSEMMixtureTrace,
    WP::Matrix{Float64}, k::Int, gauss::Vector{Int})
    T = size(WP, 2)
    U = pit.U; Ui = pit.Ui; r = pit.r; obs = pit.obs
    half_log_2pi = 0.5 * log(2pi)
    @inbounds for t in 1:T
        wp = WP[k, t]
        wp > 0 || continue
        g = 0
        for j in gauss
            if isfinite(tr.data[j, t])
                g += 1
                obs[g] = j
            end
        end
        g == 0 && continue
        for a in 1:g, b in 1:g
            U[a, b] = tr.yPcov[t, obs[a], obs[b]]
        end
        _ctsem_cholesky!(U, g) || continue
        logdet = 0.0
        rr = 0.0
        for a in 1:g
            s = tr.data[obs[a], t] - tr.yP[t, obs[a]]
            for c in 1:(a - 1)
                s -= U[c, a] * r[c]
            end
            r[a] = s / U[a, a]
            logdet += log(U[a, a])
            rr += r[a]^2
        end
        isfinite(rr) || continue
        lwp = log(wp)
        llrow = -0.5 * rr - logdet - g * half_log_2pi
        lpre = 0.0
        for a in 1:g
            _mixture_pit_add!(pit.smx, pit.snum, pit.sden, t, obs[a], lwp + lpre,
                _ctsem_normal_cdf(r[a]))
            lpre += -0.5 * r[a]^2 - log(U[a, a]) - half_log_2pi
        end
        # U^-1, upper triangular, column by column.
        for c in 1:g
            for a in c:-1:1
                s = a == c ? 1.0 : 0.0
                for i in (a + 1):c
                    s -= U[a, i] * Ui[i, c]
                end
                Ui[a, c] = s / U[a, a]
            end
        end
        for a in 1:g
            Pjj = 0.0
            Pv = 0.0
            for i in a:g
                Pjj += Ui[a, i]^2
                Pv += Ui[a, i] * r[i]
            end
            z = Pv / sqrt(Pjj)
            logcond = -0.5 * z^2 + 0.5 * log(Pjj) - half_log_2pi
            _mixture_pit_add!(pit.gmx, pit.gnum, pit.gden, t, obs[a], lwp + llrow - logcond,
                _ctsem_normal_cdf(z))
        end
    end
    return nothing
end

export ctsem_laplace_sequential_mixture
