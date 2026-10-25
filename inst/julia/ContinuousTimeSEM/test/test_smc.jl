# Tempered SMC (sample_smc.jl) on closed-form targets, and its coreset on a
# Laplace fixture.
#
# A wrong SMC does not error: it returns a plausible cloud. Two things a
# closed form can check that a reference sampler cannot. The final moments
# recover whatever happened in the stages, because the last stages' moves are
# invariant for the posterior; but the stages' log mean weights sum to the
# target's log constant only if the particles sat at each stage's density, so
# a lag shows there first (the first prototype's came out 17 nats low with
# fewer particles than dimensions). And a start that is narrow and off
# centre is the case the sampler exists for.

using LinearAlgebra, Logging, Random, Statistics

isdefined(@__MODULE__, :_fresh_linear) || include("laplace_fixtures.jl")

# A Gaussian target as an objective: `ctsem_evaluate` is all the sampler asks of
# one. Log density -(x - m)' P (x - m) / 2, whose constant is
# (d / 2) log(2 pi) - log det(P) / 2.
struct _SMCGauss
    m::Vector{Float64}
    P::Matrix{Float64}
end
function ContinuousTimeSEM.ctsem_evaluate(o::_SMCGauss, x::AbstractVector;
    gradient::Bool=true, kwargs...)
    r = x .- o.m
    Pr = o.P * r
    return (value=-dot(r, Pr) / 2, gradient=gradient ? -Pr : nothing)
end

# The first k coordinates of a rotated frame log-gamma(a) (mean digamma(a),
# variance trigamma(a)), the rest standard normal: skewed, with exact moments.
struct _SMCSkew
    Q::Matrix{Float64}
    a::Vector{Float64}
end
function ContinuousTimeSEM.ctsem_evaluate(o::_SMCSkew, x::AbstractVector;
    gradient::Bool=true, kwargs...)
    y = transpose(o.Q) * x
    k = length(o.a)
    v = 0.0
    gy = similar(y)
    @inbounds for j in eachindex(y)
        if j <= k
            e = exp(y[j])
            v += o.a[j] * y[j] - e
            gy[j] = o.a[j] - e
        else
            v -= y[j]^2 / 2
            gy[j] = -y[j]
        end
    end
    return (value=v, gradient=gradient ? o.Q * gy : nothing)
end

function _smc_gauss(d::Int, seed::Int)
    rng = Random.Xoshiro(seed)
    A = randn(rng, d, d)
    R = A * transpose(A) ./ d + 0.2I
    sd = sqrt.(diag(R))
    R = R ./ (sd * transpose(sd))
    sc = exp.(range(-1, 1; length=d))
    S = R .* (sc * transpose(sc))
    return (m=randn(rng, d), S=Symmetric(S), P=inv(Symmetric(S)))
end

@testset "tempered SMC recovers a correlated normal from a narrow, off-centre start" begin
    for (d, np) in ((20, 200), (100, 150))
        g = _smc_gauss(d, 100 + d)
        sd = sqrt.(diag(g.S))
        target = _SMCGauss(g.m, Matrix(g.P))
        # The start: 0.3 of the sds, no correlation, 3 sds off along three
        # coordinates; handed in as the Hessian the metric would be read from.
        start = g.m .+ vcat([3.0, -3.0, 3.0], zeros(d - 3)) .* sd
        hessian = -Diagonal(1 ./ (0.09 .* sd .^ 2))
        r = ContinuousTimeSEM.ctsem_sample_smc(target, start; hessian=Matrix(hessian),
            widen=1.0, nparticles=np, nmoves=3, maxdepth=6, seed=7, ndraws=10,
            min_ess=200.0, max_draws=40)
        X = r.draws
        @test size(X, 1) == d
        # The final phase starts as many chains as its target needs at the
        # draws an effective draw took in the last stage's moves, and more
        # where they fall short, at most one a particle.
        @test r.final_tau >= 1
        @test r.nchains >= ContinuousTimeSEM._smc_final_chains(r.final_tau, 200.0, 10, np)
        @test r.nchains <= np
        @test size(X, 2) == r.nchains * r.ndraws
        @test r.sampler == "smc"
        @test r.nstages >= 2
        @test all(isfinite, X)
        err = abs.(vec(mean(X; dims=2)) .- g.m) ./ sd
        ratio = vec(std(X; dims=2)) ./ sd
        @test maximum(err) < 0.25
        @test all(0.8 .< ratio .< 1.2)
        # The stages' log mean weights sum to the log constant where the
        # particles kept up with the stages. Fewer particles than dimensions
        # (the second case) is where the unshrunk metric lost 17 nats.
        logz_exact = d / 2 * log(2pi) - logdet(g.P) / 2
        @test abs(sum(r.stage_logz) - logz_exact) < 2.5
        # The final phase met the effective-size target and reports per chain.
        @test r.min_ess >= 200 * 0.8
        @test length(r.stepsize) == r.nchains
        @test all(r.stepsize .== r.stepsize[1])
        @test length(r.accept) == r.nchains * r.ndraws
        @test r.stage_beta[end] == 1.0
        @test issorted(r.stage_beta)
        # No coreset off the Laplace route, and no correction.
        @test !r.coreset
        @test all(r.stage_coreset_size .== 0)
        @test r.final_correction_ess == 1.0
    end
end

@testset "the final phase's chains are as many as the moves' autocorrelation needs" begin
    # Positions a lag apart with correlation 0.125 across the cloud: an AR(1)
    # at 0.5 a draw over three draws, whose autocorrelation time is 3.
    rng = Random.Xoshiro(41)
    n = 20000
    Xfrom = randn(rng, 3, n)
    Xto = 0.125 .* Xfrom .+ sqrt(1 - 0.125^2) .* randn(rng, 3, n)
    Xto[3, :] .= -Xfrom[3, :]
    @test ContinuousTimeSEM._smc_move_tau(Xfrom, Xto, 3) ≈ 3 rtol = 0.1
    # A negative correlation is not credited, and a coordinate that does not
    # vary is passed over.
    Xto[1:2, :] .= -Xfrom[1:2, :]
    Xfrom[2, :] .= 0.0
    @test ContinuousTimeSEM._smc_move_tau(Xfrom, Xto, 3) == 1.0
    # 1.25 times the target at tau draws each, between 16 and the particles.
    chains = ContinuousTimeSEM._smc_final_chains
    @test chains(3.0, 200.0, 20, 200) == 38
    @test chains(1.0, 200.0, 20, 200) == 16
    @test chains(1.0, 200.0, 20, 10) == 10
    @test chains(100.0, 200.0, 20, 200) == 200
    @test chains(3.0, 0.0, 20, 200) == 200
    # Every particle once, in an order whose prefixes spread over the cloud.
    o = ContinuousTimeSEM._smc_spread_order(200)
    @test sort(o) == collect(1:200)
    @test o[1:4] == [1, 101, 51, 151]
    @test maximum(diff(sort(o[1:38]))) <= 2 * cld(200, 38)
    @test sort(ContinuousTimeSEM._smc_spread_order(7)) == collect(1:7)
end

@testset "a stage's moves set the next stage's step" begin
    next = ContinuousTimeSEM._smc_next_step
    # At the target acceptance the step stands; below it the step shrinks and
    # above it grows, by at most a factor of two.
    @test next(0.3, 0.8, 0.8) ≈ 0.3
    @test 0.15 <= next(0.3, 0.6, 0.8) < 0.3
    @test 0.3 < next(0.3, 0.95, 0.8) <= 0.6
    @test next(0.3, 1.0, 0.8) == 0.6
    @test next(0.3, 0.0, 0.8) == 0.15
    @test next(0.3, NaN, 0.8) == 0.3
    # The energy error's sd grows as the step's square: at an acceptance of
    # 2 Phi(-s / 2), the step for 0.8 scales by sqrt(sd(0.8) / s).
    s = 1.0
    a = 2 * (0.5 * (1 + ContinuousTimeSEM.SpecialFunctions.erf(-s / 2 / sqrt(2))))
    sd8 = -2 * sqrt(2) * ContinuousTimeSEM.SpecialFunctions.erfinv(0.8 - 1)
    @test next(0.3, a, 0.8) ≈ 0.3 * sqrt(sd8 / s)
end

@testset "tempered SMC recovers a skewed target from its Laplace normal" begin
    d, a = 30, [0.5, 0.5, 1.0, 2.0, 4.0]
    k = length(a)
    rng = Random.Xoshiro(200)
    Q = Matrix(qr(randn(rng, d, d)).Q)
    target = _SMCSkew(Q, a)
    # digamma and trigamma at these points, exact to 1e-12.
    dig = [-1.9635100260214235, -1.9635100260214235, -0.5772156649015329,
        0.42278433509846713, 1.2561176684318005]
    trig = [4.934802200544679, 4.934802200544679, 1.6449340668482264,
        0.6449340668482264, 0.2838229557371153]
    m_true = Q * vcat(dig, zeros(d - k))
    sd_true = sqrt.(diag(Q * Diagonal(vcat(trig, ones(d - k))) * transpose(Q)))
    mode = Q * vcat(log.(a), zeros(d - k))
    lapcov = Q * Diagonal(vcat(1 ./ a, ones(d - k))) * transpose(Q)
    # The Laplace normal's mean is off by more than half an sd somewhere.
    @test maximum(abs.(mode .- m_true) ./ sd_true) > 0.5
    r = ContinuousTimeSEM.ctsem_sample_smc(target, mode; hessian=-inv(Symmetric(lapcov)),
        widen=2.0, nparticles=200, nmoves=3, maxdepth=6, seed=11, ndraws=10,
        min_ess=200.0, max_draws=40)
    X = r.draws
    err = abs.(vec(mean(X; dims=2)) .- m_true) ./ sd_true
    ratio = vec(std(X; dims=2)) ./ sd_true
    @test maximum(err) < 0.25
    @test all(0.8 .< ratio .< 1.2)
    @test r.ndivergent == 0
end

@testset "tempered SMC refuses what it cannot run" begin
    g = _smc_gauss(4, 3)
    target = _SMCGauss(g.m, Matrix(g.P))
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_sample_smc(target, g.m;
        hessian=-Matrix(g.P), nparticles=1)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_sample_smc(target, g.m;
        hessian=-Matrix(g.P), nmoves=0)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_sample_smc(target, g.m;
        hessian=-Matrix(g.P), ess_frac=1.0)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_sample_smc(target, g.m;
        hessian=-Matrix(g.P), coreset_tol=0.0)
end

@testset "tempered SMC reports its stages to the caller's sink" begin
    # A front end or a knitr run listens on a sink (progress.jl's `_deliver!`);
    # the stage lines once went straight to the console and it heard nothing.
    g = _smc_gauss(4, 3)
    target = _SMCGauss(g.m, Matrix(g.P))
    heard = Tuple{String,String}[]
    sink = (text, kind) -> (push!(heard, (String(text), String(kind))); nothing)
    ContinuousTimeSEM.ctsem_sample_smc(target, g.m; hessian=-Matrix(g.P), widen=2.0,
        nparticles=20, nmoves=1, maxdepth=4, seed=1, ndraws=5, verbose=true,
        progress_sink=sink)
    @test any(h -> occursin("tempered SMC: 20 particles", h[1]), heard)
    @test any(h -> occursin("tempering stage 1", h[1]) && h[2] == "update", heard)
    @test any(h -> h[2] == "done" && occursin("stages", h[1]), heard)
    # And nothing without verbose.
    empty!(heard)
    ContinuousTimeSEM.ctsem_sample_smc(target, g.m; hessian=-Matrix(g.P), widen=2.0,
        nparticles=20, nmoves=1, maxdepth=4, seed=1, ndraws=5, verbose=false,
        progress_sink=sink)
    @test isempty(heard)
end

# ---- the coreset ------------------------------------------------------------

@testset "a coreset's fit follows a sum of unit terms, and holds out" begin
    rng = Random.Xoshiro(5)
    N, n = 60, 240
    # Units' terms that vary in four dimensions over the cloud, each at its own
    # level, with a little of their own: a few units carry the sum.
    Z = randn(rng, 4, n)
    B = randn(rng, N, 4)
    terms = B * Z .+ 0.02 .* randn(rng, N, n) .+ randn(rng, N)
    total = vec(sum(terms; dims=1))
    cost = fill(1 / N, N)
    fit = ContinuousTimeSEM._smc_coreset_fit(terms, total, cost, 0.3, 0.5, Random.Xoshiro(1))
    @test fit !== nothing
    @test 4 <= length(fit.units) <= 12
    @test allunique(fit.units)
    @test all(fit.weights .> 0)
    @test fit.error <= 0.3
    @test fit.share <= 0.5
    predicted = vec(transpose(fit.weights) * terms[fit.units, :])
    e = predicted .- total
    e .-= mean(e)
    @test sqrt(mean(abs2, e)) <= 0.3
    # A particle the route could not use is left out of the fit, not fitted
    # to.
    terms2 = copy(terms); total2 = copy(total)
    terms2[3, 7] = NaN; total2[11] = -Inf
    fit2 = ContinuousTimeSEM._smc_coreset_fit(terms2, total2, cost, 0.3, 0.5, Random.Xoshiro(1))
    @test fit2 !== nothing && fit2.error <= 0.3
    # Terms that share nothing with the total are not carried within a small
    # budget.
    @test ContinuousTimeSEM._smc_coreset_fit(randn(rng, N, n), randn(rng, n), cost,
        0.1, 0.05, Random.Xoshiro(1)) === nothing
end

@testset "particles over copies of a one-subject filter objective draw what the serial run draws" begin
    # A filter objective with fewer subjects than the pool runs its particles
    # on deep copies, one task each (`_smc_each!`'s `:copies`). Each particle
    # draws from its own seed, and one subject's arithmetic does not depend on
    # the layout, so the draws are the serial run's. A shared workspace
    # between the tasks once left them moving on each other's buffers: early
    # ends, a posterior 0.6 of its width, and no error.
    if Threads.nthreads() > 1
        starts, times, data = _laplace_test_data(1, 1; seed=4)
        one = ctsem_objective(_LAPLACE_LINEAR_OBJECTIVE.params, starts, times, data;
            prior_index=collect(1:5), prior_scale=ones(5))
        x0 = _LAPLACE_LINEAR_VALUES[1:5]
        H = ContinuousTimeSEM.ctsem_hessian(one, x0)
        held = ContinuousTimeSEM._CTSEM_MAX_CHUNKS[]
        ContinuousTimeSEM._CTSEM_MAX_CHUNKS[] = min(4, Threads.nthreads())
        try
            common = (hessian=H, widen=2.0, nparticles=24, nmoves=2, maxdepth=5, seed=9,
                ndraws=5, min_ess=0.0, max_draws=5)
            ser = ContinuousTimeSEM.ctsem_sample_smc(one, x0; parallel=false, common...)
            par = ContinuousTimeSEM.ctsem_sample_smc(one, x0; parallel=true, common...)
            @test ser.parallel == "serial"
            @test par.parallel == "copies"
            @test par.nstages == ser.nstages
            @test par.stage_depth ≈ ser.stage_depth
            @test maximum(abs, par.draws .- ser.draws) < 1e-8
            @test all(ser.stage_moved .> 0)
        finally
            ContinuousTimeSEM._CTSEM_MAX_CHUNKS[] = held
        end
    end
end

@testset "a move on the coreset owes its error at both ends" begin
    mc = ContinuousTimeSEM._smc_move_correction
    # beta ((lp(y) - lc(y)) - (lp(x) - lc(x))); one for a particle that stayed.
    @test mc(0.5, -10.0, -12.0, -11.0, -11.5, true) ≈ 0.5 * (2.0 - 0.5)
    @test mc(0.5, -10.0, -12.0, -11.0, -11.5, false) == 0.0
    @test mc(1.0, -Inf, -12.0, -11.0, -11.5, true) == -Inf
    @test mc(1.0, -Inf, -12.0, -11.0, -11.5, false) == 0.0
end

# A sum of Gaussian unit terms the sampler can move on a coreset of, through
# the hooks the Laplace route answers (`_smc_coreset_capable` and the rest):
# exact moments and log constant, and a few units stand for the sum poorly.
struct _SMCUnits
    a::Matrix{Float64}            # d x N: each unit's centre
    B::Vector{Matrix{Float64}}    # each unit's precision
    P0::Matrix{Float64}           # the prior's precision
end
_smcu_term(o::_SMCUnits, U::Int, x) = (r = x .- view(o.a, :, U); -dot(r, o.B[U] * r) / 2)
_smcu_prior(o::_SMCUnits, x) = -dot(x, o.P0 * x) / 2
function _smcu_density(o::_SMCUnits, units, weights, x)
    v = _smcu_prior(o, x)
    g = -(o.P0 * x)
    for (k, U) in enumerate(units)
        r = x .- view(o.a, :, U)
        Br = o.B[U] * r
        v -= weights[k] * dot(r, Br) / 2
        g .-= weights[k] .* Br
    end
    return (value=v, gradient=g)
end
function ContinuousTimeSEM.ctsem_evaluate(o::_SMCUnits, x::AbstractVector;
    gradient::Bool=true, kwargs...)
    r = _smcu_density(o, 1:size(o.a, 2), ones(size(o.a, 2)), collect(Float64, x))
    return (value=r.value, gradient=gradient ? r.gradient : nothing)
end
ContinuousTimeSEM._smc_coreset_capable(::_SMCUnits) = true
ContinuousTimeSEM._smc_unit_cost(o::_SMCUnits) = fill(1 / size(o.a, 2), size(o.a, 2))
function ContinuousTimeSEM._smc_unit_terms(o::_SMCUnits, X::AbstractMatrix{Float64})
    terms = [_smcu_term(o, U, view(X, :, j)) for U in axes(o.a, 2), j in axes(X, 2)]
    total = vec(sum(terms; dims=1))
    return (value=[_smcu_prior(o, view(X, :, j)) + total[j] for j in axes(X, 2)],
        total=total, terms=terms)
end
ContinuousTimeSEM._smc_coreset_values(o::_SMCUnits, X::AbstractMatrix{Float64}, cs) =
    [_smcu_density(o, cs.units, cs.weights, collect(Float64, view(X, :, j))).value
        for j in axes(X, 2)]
ContinuousTimeSEM._smc_coreset_eval(o::_SMCUnits, cs, x::AbstractVector) =
    _smcu_density(o, cs.units, cs.weights, collect(Float64, x))

@testset "a poor coreset slows SMC but leaves its posterior and constant exact" begin
    # Forty units whose precisions are rank one along random directions: their
    # sum is a well-conditioned normal in four dimensions, and a dozen of them
    # stand for it with errors of a log unit and more, so the correction does
    # real work (its effective size falls below a tenth at some stage). One
    # move a stage, where a weight that assumed the moves had mixed to the
    # coreset's density is furthest from right: before the backward kernel's
    # weights the stages' log constant was 270-650 nats off on this target;
    # with them, within a nat (2026-10-23).
    d, N = 4, 40
    rng = Random.Xoshiro(3)
    a = 1.5 .* randn(rng, d, N)
    B = [begin
        v = randn(rng, d)
        (2.0 * (0.5 + rand(rng))) .* (v * transpose(v)) ./ dot(v, v)
    end for _ in 1:N]
    o = _SMCUnits(a, B, Matrix(0.05I, d, d))
    P = o.P0 .+ sum(B)
    b = sum(B[U] * a[:, U] for U in 1:N)
    S = inv(Symmetric(P))
    m = S * b
    sd = sqrt.(diag(S))
    logz = -sum(dot(a[:, U], B[U] * a[:, U]) for U in 1:N) / 2 + dot(b, S * b) / 2 +
        d / 2 * log(2pi) - logdet(Symmetric(P)) / 2
    start = m .+ vcat([2.0, -2.0], zeros(d - 2)) .* sd
    H = -Matrix(Diagonal(1 ./ (0.25 .* sd .^ 2)))
    # Its final chains keep about half their effective size, which warns.
    r = @test_logs min_level=Logging.Error ContinuousTimeSEM.ctsem_sample_smc(o, start;
        hessian=H, widen=1.0, nparticles=400, nmoves=1, maxdepth=5, seed=2, ndraws=10,
        min_ess=0.0, max_draws=10, coreset=true, coreset_tol=1.0, parallel=false)
    @test r.coreset
    @test all(r.stage_coreset_size .> 0)
    @test minimum(r.stage_correction_ess) < 0.5
    X = r.draws
    @test maximum(abs.(vec(mean(X; dims=2)) .- m) ./ sd) < 0.2
    @test all(0.85 .< vec(std(X; dims=2)) ./ sd .< 1.15)
    @test abs(sum(r.stage_logz) - logz) < 1.5
end

@testset "tempered SMC on a coreset samples the same posterior as the full Laplace objective" begin
    # Forty subjects of the shared linear model, a standard normal prior on
    # every raw parameter (without one the scales run to zero), two random
    # effects: a sum of forty units' terms, which a coreset of a few follows.
    starts, times, data = _laplace_test_data(40, 1; seed=3)
    prior = ctsem_objective(_LAPLACE_LINEAR_OBJECTIVE.params, starts, times, data;
        prior_index=collect(1:8), prior_scale=ones(8))
    laplace = ctsem_laplace_objective(prior, [1, 2], [6, 7], [8], [1.0, 1.0])
    x0 = copy(_LAPLACE_LINEAR_VALUES)
    H = ContinuousTimeSEM.ctsem_hessian(laplace, x0)
    common = (hessian=H, widen=1.5, nparticles=100, nmoves=2, maxdepth=5, seed=3,
        ndraws=10, min_ess=80.0, max_draws=40, parallel=false)
    full = ContinuousTimeSEM.ctsem_sample_smc(laplace, x0; coreset=false, common...)
    core = ContinuousTimeSEM.ctsem_sample_smc(laplace, x0; coreset=true, common...)
    @test !full.coreset
    @test all(full.stage_coreset_size .== 0)
    @test core.coreset
    @test core.nunits == 40
    # Some stage moved on a coreset, and none spent more than the budget's
    # share of the units.
    @test any(core.stage_coreset_size .> 0)
    # A stage run in full records a share of one.
    @test all(i -> core.stage_coreset_size[i] == 0 ||
        core.stage_coreset_share[i] <= 0.5 + 1e-12, eachindex(core.stage_coreset_size))
    @test all(s -> s == 0 || s <= 20, core.stage_coreset_size)
    # The correction's effective size is reported a stage, in (0, 1].
    @test all(0 .< core.stage_correction_ess .<= 1)
    @test 0 < core.final_correction_ess <= 1
    @test core.parallel == "serial"
    # The same posterior, within the Monte Carlo error of two runs at an
    # effective size of about 80.
    m1 = vec(mean(full.draws; dims=2)); s1 = vec(std(full.draws; dims=2))
    m2 = vec(mean(core.draws; dims=2)); s2 = vec(std(core.draws; dims=2))
    @test maximum(abs.(m1 .- m2) ./ s1) < 0.45
    @test all(0.65 .< s2 ./ s1 .< 1.5)
    # With the backward kernel's weights the stages' log mean weights still
    # estimate the log constant on a coreset: the weights before the backward
    # kernel were hundreds of nats off. How close is the Monte Carlo error of
    # two runs of 100 particles: over seeds 3-8 the two differed by -4.3 to
    # +1.9 (sd 2.2, mean -1.2 +- 0.9), so the bar is three of those sds. The
    # constant itself is tested against its closed form above.
    @test abs(sum(core.stage_logz) - sum(full.stage_logz)) < 7.0
    if Threads.nthreads() > 1
        # Particles spread over the pool on copies: the same seeds per
        # particle, so the same posterior, laid out differently.
        par = ContinuousTimeSEM.ctsem_sample_smc(laplace, x0; coreset=true,
            merge(common, (parallel=true,))...)
        @test par.parallel == "pool"
        m3 = vec(mean(par.draws; dims=2)); s3 = vec(std(par.draws; dims=2))
        @test maximum(abs.(m3 .- m1) ./ s1) < 0.45
        @test all(0.65 .< s3 ./ s1 .< 1.5)
    end
end

# ---- NUTS on a coreset (sample_surrogate.jl) ---------------------------------

# Two dimensions, forty rank-one units: a few stand for the sum over the
# posterior, imperfectly, so the endpoint test has work to do.
function _surrogate_units(seed)
    d, N = 2, 40
    rng = Random.Xoshiro(seed)
    a = 1.5 .* randn(rng, d, N)
    B = [begin
        v = randn(rng, d)
        (2.0 * (0.5 + rand(rng))) .* (v * transpose(v)) ./ dot(v, v)
    end for _ in 1:N]
    o = _SMCUnits(a, B, Matrix(0.05I, d, d))
    P = o.P0 .+ sum(B)
    S = inv(Symmetric(P))
    return (o=o, P=P, m=S * sum(B[U] * a[:, U] for U in 1:N), sd=sqrt.(diag(S)))
end

@testset "a draw on a coreset is tested against the full density at both ends" begin
    t = _surrogate_units(5)
    d = ContinuousTimeSEM._SurrogateDensity(t.o, :adjoint, 2, 1)
    d.coreset = ContinuousTimeSEM.CTSEMCoreset([1, 2, 3], [10.0, 10.0, 10.0])
    full(x) = ContinuousTimeSEM.ctsem_evaluate(t.o, x).value
    core(x) = _smcu_density(t.o, d.coreset.units, d.coreset.weights, x).value
    x0 = t.m .+ [0.5, -0.5] .* t.sd
    y = t.m .+ [-0.5, 1.0] .* t.sd
    d.lp = full(x0); d.lc = core(x0)
    g = zeros(2)
    ContinuousTimeSEM._surrogate_begin!(d, x0, g)
    # pi(y) c(x0) / (pi(x0) c(y)): accepted where a uniform falls under it.
    logratio = (full(y) - core(y)) - (full(x0) - core(x0))
    rng = Random.Xoshiro(11)
    accept = log(rand(copy(rng))) < logratio
    x = copy(y); g .= 1.0
    lp = ContinuousTimeSEM._surrogate_correct!(d, x, g, core(y), rng)
    @test d.proposed == 1
    @test d.accepted == (accept ? 1 : 0)
    @test x == (accept ? y : x0)
    @test lp == (accept ? core(y) : core(x0))
    accept && @test d.lp == full(y)
    # A transition that did not move owes nothing and asks nothing.
    ContinuousTimeSEM._surrogate_begin!(d, x, g)
    @test ContinuousTimeSEM._surrogate_correct!(d, x, g, lp, rng) == lp
    @test d.proposed == 1
    # A plain density has nothing to test.
    @test ContinuousTimeSEM._surrogate_correct!((g, x) -> 0.0, x, g, -3.0, rng) == -3.0
end

@testset "NUTS on a coreset samples the full posterior" begin
    t = _surrogate_units(5)
    common = (hessian=-t.P, nchains=2, nwarmup=300, ndraws=3000, maxdepth=6,
        seed=4, min_ess=0.0)
    core = ContinuousTimeSEM.ctsem_sample_marginal(t.o, t.m; coreset=true, common...)
    full = ContinuousTimeSEM.ctsem_sample_marginal(t.o, t.m; coreset=false, common...)
    @test core.coreset
    @test !full.coreset
    @test core.nunits == 40
    # Each chain fitted its own coreset from its warm-up, within the budget's
    # share of the units, and the full density accepted most of its moves.
    @test all(0 .< core.coreset_size .<= 20)
    @test all(core.coreset_share .<= 0.5 + 1e-12)
    rate = core.coreset_accepted ./ core.coreset_proposed
    @test all(0.3 .< rate .<= 1)
    @test !any(core.coreset_abandoned)
    # The full posterior, exactly: within the Monte Carlo error of 6000 draws.
    X = core.draws
    @test maximum(abs.(vec(mean(X; dims=2)) .- t.m) ./ t.sd) < 0.1
    @test all(0.9 .< vec(std(X; dims=2)) ./ t.sd .< 1.1)
end
