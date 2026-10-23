# Tempered SMC (sample_smc.jl) on closed-form targets.
#
# A wrong SMC does not error: it returns a plausible cloud. Two things a
# closed form can check that a reference sampler cannot. The final moments
# recover whatever happened in the stages, because the last stages' moves are
# invariant for the posterior; but the stages' log mean weights sum to the
# target's log constant only if the particles sat at each stage's density, so
# a lag shows there first (the first prototype's came out 17 nats low with
# fewer particles than dimensions). And a start that is narrow and off
# centre is the case the sampler exists for.

using LinearAlgebra, Random, Statistics

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
        @test r.nchains == np
        @test size(X, 2) == np * r.ndraws
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
        @test length(r.stepsize) == np
        @test all(r.stepsize .== r.stepsize[1])
        @test length(r.accept) == np * r.ndraws
        @test r.stage_beta[end] == 1.0
        @test issorted(r.stage_beta)
    end
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
end
