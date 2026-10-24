# Doubly adaptive importance sampling (dais.jl) on targets with exact
# moments: a rotated log-gamma skew whose Laplace normal is off, a normal it
# should leave alone, a direction the density does not depend on, and a
# density that refuses a tail.

using LinearAlgebra, Random, Statistics

# A target as an objective: `ctsem_evaluate` is all DAIS asks of one.
struct _DaisTarget{F}
    f::F    # x -> (value, gradient)
end
function ContinuousTimeSEM.ctsem_evaluate(o::_DaisTarget, x::AbstractVector;
    gradient::Bool=true, kwargs...)
    v, g = o.f(collect(Float64, x))
    return (value=v, gradient=gradient ? g : nothing)
end

# The first k coordinates of a rotated frame are log-gamma(a) (density
# exp(a y - exp(y)); mean digamma(a), variance trigamma(a)), the rest
# standard normal.
function _dais_skew(d::Int, a::Vector{Float64}; seed::Int=11)
    k = length(a)
    Q = Matrix(qr(randn(Random.Xoshiro(seed), d, d)).Q)
    f = x -> begin
        y = transpose(Q) * x
        e = exp.(y[1:k])
        (sum(a .* y[1:k] .- e) - sum(abs2, y[(k + 1):end]) / 2,
            Q * vcat(a .- e, -y[(k + 1):end]))
    end
    dig = Dict(0.5 => -1.9635100260214235, 1.0 => -0.5772156649015329,
        2.0 => 0.42278433509846713)
    trig = Dict(0.5 => 4.934802200544679, 1.0 => 1.6449340668482264,
        2.0 => 0.6449340668482264)
    return (target=_DaisTarget(f), mode=Q * vcat(log.(a), zeros(d - k)),
        lapcov=Q * Diagonal(vcat(1 ./ a, ones(d - k))) * transpose(Q),
        mean=Q * vcat([dig[x] for x in a], zeros(d - k)),
        sd=sqrt.(diag(Q * Diagonal(vcat([trig[x] for x in a], ones(d - k))) * transpose(Q))))
end

@testset "DAIS moves a normal from the mode to a skewed posterior's moments" begin
    t = _dais_skew(20, [0.5, 0.5, 1.0, 2.0])
    # The Laplace normal's mean is off by over half a posterior sd somewhere.
    @test maximum(abs.(t.mode .- t.mean) ./ t.sd) > 0.5
    r = ctsem_dais(t.target, t.mode, t.lapcov; draws=64, seed=12)
    @test r.used
    @test r.settled
    @test r.averaged > 0
    @test maximum(abs.(r.mean .- t.mean) ./ t.sd) < 0.15
    # The damped fixed point's sds sit a little under the posterior's: over
    # ten seeds the smallest ratio was 0.82-0.86, the median 0.94-0.96.
    ratio = sqrt.(diag(r.cov)) ./ t.sd
    @test all(0.78 .< ratio .< 1.1)
    # Even draws an iteration, counted; the trace carries one row an iteration.
    @test r.draws == 64
    @test r.gradients == 64 * r.iterations
    @test size(r.trace) == (r.iterations, length(r.trace_names))
    @test all(r.trace[:, findfirst(==("usable"), r.trace_names)] .== 64)
    # A normal target from its own normal: nothing to move.
    P = inv(Symmetric(t.lapcov))
    g = _DaisTarget(x -> (-dot(x, P * x) / 2, -(P * x)))
    r0 = ctsem_dais(g, zeros(20), t.lapcov; draws=64, seed=13)
    @test r0.settled
    @test maximum(abs.(r0.mean) ./ sqrt.(diag(t.lapcov))) < 0.1
    @test isapprox(r0.cov, t.lapcov; rtol=0.15)
end

@testset "DAIS holds a direction the data do not identify, or ends as diverged" begin
    # Skewed in two coordinates and flat in a third: a moment-matching step
    # multiplies that variance by 1 + gamma every iteration.
    a = [0.5, 1.0]
    target = _DaisTarget(x -> begin
        e = exp.(x[1:2])
        (sum(a .* x[1:2] .- e), vcat(a .- e, 0.0))
    end)
    mode = vcat(log.(a), 0.0)
    cov0 = Matrix(Diagonal(vcat(1 ./ a, 1.0)))
    held = ctsem_dais(target, mode, cov0; draws=64, seed=18, hold=reshape([0.0, 0.0, 1.0], 3, 1))
    @test held.used
    dig = [-1.9635100260214235, -0.5772156649015329]
    trig = [4.934802200544679, 1.6449340668482264]
    @test maximum(abs.(held.mean[1:2] .- dig) ./ sqrt.(trig)) < 0.15
    @test held.mean[3] ≈ 0 atol = 1e-12
    @test held.cov[3, 3] ≈ 1 atol = 1e-10
    @test all(abs.(held.cov[1:2, 3]) .< 1e-10)
    free = ctsem_dais(target, mode, cov0; draws=64, seed=18)
    @test !free.used
    @test free.diverged
    @test free.growth > 1e4
    @test free.cov === nothing
end

@testset "DAIS gives the draws a density cannot use no weight and no gradient" begin
    # The skewed target refusing every draw past 1.5 on its first coordinate:
    # a value of -Inf and an infinite gradient, as a filter overflowing in a
    # tail returns. Zero weight times an infinite gradient was NaN once.
    t = _dais_skew(6, [0.5, 1.0])
    refusing = _DaisTarget(x -> begin
        v, g = t.target.f(x)
        x[1] > 1.5 ? (-Inf, fill(Inf, length(x))) : (v, g)
    end)
    r = ctsem_dais(refusing, t.mode, t.lapcov; draws=64, seed=19)
    @test r.used
    @test r.unusable > 0
    @test all(isfinite, r.mean)
    @test all(isfinite, r.cov)
    # And a normal none of whose draws can be used leaves the check's.
    none = ctsem_dais(_DaisTarget(x -> (-Inf, fill(Inf, length(x)))), t.mode, t.lapcov;
        draws=64, seed=20)
    @test !none.used
    @test none.unusable == 64
    @test none.cov === nothing
end

@testset "DAIS reports to the caller's sink" begin
    t = _dais_skew(6, [0.5, 1.0])
    heard = String[]
    ctsem_dais(t.target, t.mode, t.lapcov; draws=64, seed=21, maxiter=5, verbose=true,
        progress_sink=(text, kind) -> (push!(heard, String(text)); nothing))
    @test any(occursin("DAIS", h) for h in heard)
end
