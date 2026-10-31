# Reference integrals for the bench, included into ContinuousTimeSEM at run
# time (harness.R: Core.eval(ContinuousTimeSEM, :(include(path)))). Not part of
# the engine and not loaded by it.
#
# Copied from the gaps job's probes (session 00fd0b41, g2/probe.jl and
# g2/probe_is.jl; review/LAPLACE-gated-gaps-2026-09-24.md sections 1d and 3a),
# renamed with a bench_ prefix so nothing here can collide with an engine name.
# They reach into engine internals (_laplace_unit_curvature and friends): if a
# build renames those, the harness records the error and the cell carries no
# reference rather than a wrong one.
#
# What each computes, for one unit U at theta:
#   log int exp(g_U(u)) du / (2pi)^(d/2)  =  log E_{N(0,I)} L_U(u)
# which is the unit's exact log marginal likelihood. The penalised exact log
# likelihood of a fit is the sum over units plus the objective's own prior
# term (harness.R adds it).
#
# bench_probe_reference: directions whose eigenvalue of M at the mode is below
# `softcut` are integrated by a trapezoid over a wide grid; for each grid point
# the stiff coordinates' conditional mode is re-solved by Newton and integrated
# by adaptive Gauss-Hermite. Nothing reuses the Laplace term. softcut 3.5 is
# the setting the gaps note settled on (at 1.5 it missed mass on non-Gaussian
# directions, by up to 3.9 nats per estimate); it agreed with importance
# sampling to a median 0.009 and at most 0.18 over the units checked.
# Cost grows as (2 halfwidth / h)^nsoft, so it is used for d <= 3 only.
#
# bench_probe_is: importance sampling with a multivariate t (df 4) proposal at
# the mode, scale inflate / sqrt(max(eigenvalue, 1)) per eigendirection, drawn
# through the symmetric root so the sample is a continuous function of M.
using LinearAlgebra, Random

function bench_probe_units(laplace, theta::Vector{Float64})
    _laplace_ensure_pool!(laplace)
    ctsem_laplace_evaluate(laplace, theta; gradient=false)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    nunits = length(laplace.units.members)
    dmax = maximum(laplace.units.dims)
    out = fill(NaN, nunits, 6 + dmax)
    for U in 1:nunits
        u = laplace.modes[U]
        d = length(u)
        blocks = laplace.units.blocks[U]
        M = _laplace_unit_curvature(laplace, U, theta, Ls, u)
        D = _laplace_block_dense(M, blocks, d)
        ev = eigvals(Symmetric(D))
        g = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws).value
        ld = all(>(0), ev) ? sum(log.(ev)) : NaN
        out[U, 1] = U
        out[U, 2] = g
        out[U, 3] = ld
        out[U, 4] = isnan(ld) ? NaN : g - max(ld, 0.0) / 2
        out[U, 5] = g - sum(log.(max.(ev, 1.0))) / 2
        out[U, 6] = g - ld / 2
        out[U, 7:(6 + d)] .= ev
    end
    return out
end

function _bench_probe_g(laplace, U, theta, Ls, u, aws)
    return _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws)
end

function bench_probe_reference(laplace, theta::Vector{Float64}, U::Integer;
        softcut=3.5, nstiff=9, h=0.1, halfwidth=6.0, resolve=true)
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    uhat = copy(laplace.modes[U])
    d = length(uhat)
    blocks = laplace.units.blocks[U]
    M0 = _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta, Ls, uhat), blocks, d)
    E = eigen(Symmetric(M0))
    soft = findall(<(softcut), E.values)
    stiff = setdiff(1:d, soft)
    Vs = E.vectors[:, soft]
    Vh = E.vectors[:, stiff]
    ns = length(soft)
    nh = length(stiff)
    x, w = _gauss_hermite(nstiff)
    function stiff_integral(ubase::Vector{Float64}, zstart::Vector{Float64})
        local z, uu, r, gz, Mf, Mz, F, step, t, g0, g1, ez, lam, S, terms, xi, zz, gv, lw, pk, val, accepted
        nh == 0 && return (_bench_probe_g(laplace, U, theta, Ls, ubase, aws).value, zstart)
        z = copy(zstart)
        for it in 1:100
            uu = ubase .+ Vh * z
            r = _bench_probe_g(laplace, U, theta, Ls, uu, aws)
            gz = transpose(Vh) * r.gradient
            maximum(abs, gz) < 1e-9 && break
            Mf = _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta, Ls, uu), blocks, d)
            Mz = Symmetric(transpose(Vh) * Mf * Vh)
            F = cholesky(Mz; check=false)
            if !issuccess(F)
                Mz = Symmetric(Matrix(Mz) + (abs(eigmin(Mz)) + 1.0) * I)
                F = cholesky(Mz)
            end
            step = F \ gz
            t = 1.0
            g0 = r.value
            accepted = false
            for _ in 1:30
                g1 = _bench_probe_g(laplace, U, theta, Ls, ubase .+ Vh * (z .+ t .* step), aws).value
                if isfinite(g1) && g1 >= g0 - 1e-12
                    accepted = true
                    break
                end
                t /= 2
            end
            accepted || break
            z .+= t .* step
        end
        uu = ubase .+ Vh * z
        Mf = _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta, Ls, uu), blocks, d)
        Mz = Symmetric(transpose(Vh) * Mf * Vh)
        ez = eigen(Mz)
        lam = max.(ez.values, 1e-8)
        S = ez.vectors * Diagonal(1 ./ sqrt.(lam))
        terms = Float64[]
        for idx in Iterators.product(ntuple(_ -> 1:nstiff, nh)...)
            xi = [x[i] for i in idx]
            zz = z .+ sqrt(2) .* (S * xi)
            gv = _bench_probe_g(laplace, U, theta, Ls, ubase .+ Vh * zz, aws).value
            lw = sum(log(w[i]) + x[i]^2 for i in idx)
            push!(terms, isfinite(gv) ? gv + lw : -Inf)
        end
        pk = maximum(terms)
        val = pk + log(sum(exp.(terms .- pk))) - sum(log.(lam)) / 2 +
            nh * log(2) / 2
        return (val, z)
    end
    if ns == 0
        v, _ = stiff_integral(uhat, zeros(nh))
        return (value=v - d * log(2pi) / 2, nsoft=0, eig=E.values)
    end
    p = transpose(Vs) * uhat
    ranges = [collect((min(0.0, -p[j]) - halfwidth):h:(max(0.0, -p[j]) + halfwidth)) for j in 1:ns]
    terms = Float64[]
    zstart = zeros(nh)
    for idx in Iterators.product(ranges...)
        t = collect(idx)
        ubase = uhat .+ Vs * t
        # Two starts, the continuation and the unit mode's own stiff
        # coordinates, keeping the larger: a continuation that has walked into
        # a poor basin otherwise carries every later node with it.
        v, zz = stiff_integral(ubase, resolve ? zstart : zeros(nh))
        if resolve
            v0, z0 = stiff_integral(ubase, zeros(nh))
            if !isfinite(v) || (isfinite(v0) && v0 > v)
                v, zz = v0, z0
            end
        end
        resolve && all(isfinite, zz) && (zstart = zz)
        push!(terms, isfinite(v) ? v : -Inf)
    end
    pk = maximum(terms)
    val = pk + log(sum(exp.(terms .- pk))) + ns * log(h)
    return (value=val - d * log(2pi) / 2, nsoft=ns, eig=E.values)
end

function bench_probe_is(laplace, theta::Vector{Float64}, U::Integer; n::Integer=100_000,
        df::Real=4.0, inflate::Real=1.5, seed::Integer=1)
    rng = MersenneTwister(seed)
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    uhat = copy(laplace.modes[U]); d = length(uhat)
    M = _laplace_block_dense(_laplace_unit_curvature(laplace, U, theta, Ls, uhat),
        laplace.units.blocks[U], d)
    E = eigen(Symmetric(M))
    s = inflate ./ sqrt.(max.(E.values, 1.0))
    # The symmetric root, a continuous function of M. Laid along the
    # eigenvectors as the solver returns them, a sign or a rotation it chose
    # from M's last bits re-drew the whole sample: the nested cells' references
    # moved 0.01-0.03 nats between builds whose curvatures differed only in
    # rounding (rescue grid, juliaFit 3a780671 against 70f7d34a).
    S = E.vectors * Diagonal(s) * transpose(E.vectors)
    logdetS = sum(log, s)
    c = _bench_lgamma((df + d) / 2) - _bench_lgamma(df / 2) - d / 2 * log(df * pi)
    lw = Vector{Float64}(undef, n)
    for i in 1:n
        z = randn(rng, d); wv = sum(abs2, randn(rng, Int(df))) / df
        xv = z ./ sqrt(wv)
        u = uhat .+ S * xv
        g = _laplace_unit_objective_gradient(laplace, U, theta, Ls, u, aws).value
        lq = c - (df + d) / 2 * log1p(sum(abs2, xv) / df) - logdetS
        lw[i] = isfinite(g) ? g - lq : -Inf
    end
    pk = maximum(lw)
    wts = exp.(lw .- pk)
    m = sum(wts) / n
    se_log = sqrt(sum(abs2, wts .- m) / (n - 1)) / sqrt(n) / m
    ess = sum(wts)^2 / sum(abs2, wts)
    return (value=pk + log(m) - d * log(2pi) / 2, se=se_log, ess=ess, eig=E.values)
end

# The contamination control: seconds per value-and-gradient evaluation at `x`,
# timed inside Julia so the bridge's round trip is not in it. One untimed call,
# one call to size the batch at about `target` seconds, then `reps` batches.
function bench_time_gradient(objective, x::Vector{Float64}; target::Real=0.25,
        reps::Integer=3, maxk::Integer=5000)
    ctsem_evaluate(objective, x; gradient=true)
    t1 = @elapsed ctsem_evaluate(objective, x; gradient=true)
    k = clamp(ceil(Int, target / max(t1, 1e-6)), 1, maxk)
    out = zeros(reps)
    for r in 1:reps
        out[r] = (@elapsed for _ in 1:k
            ctsem_evaluate(objective, x; gradient=true)
        end) / k
    end
    return (seconds=out, batch=k, first=t1)
end

# Seconds a call of `f` takes, over `k` calls.
function _bench_per(f, k::Integer)
    t = time_ns()
    for _ in 1:k
        f()
    end
    return (time_ns() - t) / 1e9 / k
end

# What the pricers' unit of cost is worth on this objective at `x`: seconds a
# value and a gradient take, in batches of about `target` seconds interleaved
# so the machine's load moves both alike -- searched, as a sampler evaluates
# outside an optimisation, and anchored, as an optimisation evaluates, at
# points a short step either side of the last, the anchors moved after each as
# an accepted step moves them; the operation counts of one of each; one
# Hessian's seconds where it is expected to take under `hessian_cap`; and on a
# sum over units, single-unit coreset evaluations of `nunits` units spread
# over the unit weights the pricers read (`_laplace_unit_weights`), for how
# well rows times members predicts a unit's share.
function bench_time_costs(objective, x::Vector{Float64}; target::Real=0.25,
        reps::Integer=3, nunits::Integer=40, hessian_cap::Real=120.0)
    x = copy(x)
    ev(g, at=x) = ctsem_evaluate(objective, at; gradient=g)
    ev(false); ev(true)
    t1 = @elapsed ev(true)
    k = clamp(ceil(Int, target / max(t1, 1e-6)), 1, 2000)
    tv = zeros(reps); tg = zeros(reps)
    for r in 1:reps
        tv[r] = _bench_per(() -> ev(false), k)
        tg[r] = _bench_per(() -> ev(true), k)
    end
    ctsem_reset_opcounts!(); ev(false); ov = collect(values(ctsem_opcounts()))
    ctsem_reset_opcounts!(); ev(true); og = collect(values(ctsem_opcounts()))
    opnames = string.(collect(keys(ctsem_opcounts())))
    # Anchored: a fixed direction of length 0.05 in raw units, alternated.
    av = fill(NaN, reps); ag = fill(NaN, reps)
    _ctsem_anchor_start!(objective, x)
    anchored = objective isa CTSEMLaplaceObjective && !isempty(objective.anchors)
    if anchored
        try
            d = [sin(1.618 * i) for i in eachindex(x)]
            d .*= 0.05 / max(sqrt(sum(abs2, d)), eps())
            pts = (x .+ d, x .- d)
            j = Ref(0)
            step!(g) = (j[] += 1; ev(g, pts[mod1(j[], 2)]); _ctsem_anchor_here!(objective))
            step!(false); step!(true)
            for r in 1:reps
                av[r] = _bench_per(() -> step!(false), k)
                ag[r] = _bench_per(() -> step!(true), k)
            end
        finally
            ctsem_laplace_clear_anchors!(objective)
        end
    end
    # One Hessian, the finish's unit of cost in gradients.
    th = NaN
    if 2 * length(x) * median(tg) < hessian_cap
        th = try
            ctsem_hessian(objective, x)
            @elapsed ctsem_hessian(objective, x)
        catch
            NaN
        end
    end
    # Units: a coreset of one, value and gradient.
    L = _ctsem_unit_sum(objective)
    uw = Float64[]; usel = Int[]; uv = Float64[]; ug = Float64[]; nU = 0
    if L !== nothing
        w = _laplace_unit_weights(L)
        nU = length(w)
        if nU >= 2
            order = sortperm(w)
            usel = unique(order[round.(Int, range(1, nU; length=min(nunits, nU)))])
            X = reshape(copy(x), :, 1)
            for U in usel
                ctsem_coreset_batch(L, X, [U], [1.0])
                ctsem_coreset_batch_gradient(L, X, [U], [1.0])
                t = @elapsed ctsem_coreset_batch(L, X, [U], [1.0])
                ku = clamp(ceil(Int, 0.02 / max(t, 1e-7)), 1, 500)
                push!(uv, _bench_per(() -> ctsem_coreset_batch(L, X, [U], [1.0]), ku))
                push!(ug, _bench_per(() -> ctsem_coreset_batch_gradient(L, X, [U], [1.0]), ku))
                push!(uw, w[U])
            end
        end
    end
    return (value=tv, gradient=tg, batch=k, anchored=anchored,
        value_anchored=av, gradient_anchored=ag, hessian=th, npar=length(x),
        opnames=opnames, ops_value=Float64.(ov), ops_gradient=Float64.(og),
        nunits=nU, unit_index=isempty(usel) ? [0] : usel,
        unit_weight=isempty(uw) ? [NaN] : uw, unit_value=isempty(uv) ? [NaN] : uv,
        unit_gradient=isempty(ug) ? [NaN] : ug)
end

# log Gamma at positive integers and half-integers, by recursion.
function _bench_lgamma(x::Real)
    r = 0.0
    while x > 1.25
        x -= 1; r += log(x)
    end
    return r + (isapprox(x, 0.5) ? log(sqrt(pi)) : 0.0)
end

