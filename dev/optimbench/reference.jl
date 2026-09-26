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
# the mode, scale inflate / sqrt(max(eigenvalue, 1)) per eigendirection.
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
    S = E.vectors * Diagonal(inflate ./ sqrt.(max.(E.values, 1.0)))
    logdetS = sum(log, inflate ./ sqrt.(max.(E.values, 1.0)))
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

# log Gamma at positive integers and half-integers, by recursion.
function _bench_lgamma(x::Real)
    r = 0.0
    while x > 1.25
        x -= 1; r += log(x)
    end
    return r + (isapprox(x, 0.5) ? log(sqrt(pi)) : 0.0)
end

# The quadrature continuation's Hessian by FORWARD differences of the hybrid's
# exact gradient: `ctsem_laplace_continuation_hessian` (laplace_continuation.jl)
# with the same steps, warm starts and restores, but one gradient per column
# plus one at `values` instead of two per column. For the speed job's
# comparison of the standard errors each gives (harness.R section 5b).
# Untyped on purpose: an annotation naming CTSEMLaplaceContinuation would stop
# this file from loading on a build without that type, and with it every
# reference integral above.
function bench_continuation_hessian_forward(o, values::AbstractVector; step::Real=1e-4)
    x = collect(Float64, values)
    n = length(x)
    H = fill(NaN, n, n)
    previous = ctsem_set_warm_start!(false)
    try
        ctsem_laplace_evaluate(o.rest, x; gradient=false)
    finally
        ctsem_set_warm_start!(previous)
    end
    base = deepcopy(o.rest.modes)
    restore!() = (for U in eachindex(base); o.rest.modes[U] = copy(base[U]); end)
    previous = ctsem_set_warm_start!(true)
    try
        restore!()
        g0 = ctsem_laplace_continuation_evaluate(o, x; gradient=true)
        g0.converged || return H
        for j in 1:n
            h = step * max(1.0, abs(x[j]))
            plus = copy(x); plus[j] += h
            restore!()
            gp = ctsem_laplace_continuation_evaluate(o, plus; gradient=true)
            gp.converged || continue
            H[:, j] = (gp.gradient .- g0.gradient) ./ h
        end
    finally
        ctsem_set_warm_start!(previous)
        restore!()
    end
    return (H .+ transpose(H)) ./ 2
end
