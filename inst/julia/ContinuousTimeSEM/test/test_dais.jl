# Doubly adaptive importance sampling (dais.jl) on targets with exact
# moments: a rotated log-gamma skew whose Laplace normal is off, a normal it
# should leave alone, a direction the density does not depend on, and a
# density that refuses a tail.

using LinearAlgebra, Random, Statistics
isdefined(@__MODULE__, :_fresh_linear) || include("laplace_fixtures.jl")

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
    # It stops once its draws can show every parameter is held to its
    # precision, not after a fixed count: over twelve seeds 7-12 iterations,
    # where the windows fixed at 20 stopped at 61 at the earliest.
    @test r.iterations < 20
    # Every parameter's remaining distance within 0.2 sd, its Monte Carlo
    # error a fraction of that, and the damped fixed point's own bias: the
    # largest of twenty was 0.05-0.22 over twelve seeds.
    @test maximum(abs.(r.mean .- t.mean) ./ t.sd) < 0.25
    # The damped fixed point's sds sit a little under the posterior's: over
    # twelve seeds the smallest ratio was 0.78-0.91, the median 0.92-0.96.
    ratio = sqrt.(diag(r.cov)) ./ t.sd
    @test all(0.75 .< ratio .< 1.1)
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
    @test r0.iterations < 30
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
    # Within the 0.2 sd DAIS stops at, and its noise: over 24 seeds the
    # larger of the two was 0.01-0.23.
    @test maximum(abs.(held.mean[1:2] .- dig) ./ sqrt.(trig)) < 0.25
    @test held.mean[3] ≈ 0 atol = 1e-12
    @test held.cov[3, 3] ≈ 1 atol = 1e-10
    @test all(abs.(held.cov[1:2, 3]) .< 1e-10)
    free = ctsem_dais(target, mode, cov0; draws=64, seed=18)
    @test !free.used
    @test free.diverged
    @test free.growth > 1e4
    @test free.cov === nothing
end

@testset "DAIS draws the same for one seed whatever eigenbasis LAPACK returns" begin
    # Five equal eigenvalues, whose eigenvectors LAPACK picks by roundoff, as
    # it does every column's sign. Drawn through the eigenbasis, one seed's
    # draws turned with a 1e-13 change of the covariance.
    t = _dais_skew(6, [0.5, 1.0])
    noise = 1e-13 .* randn(Random.Xoshiro(23), 6, 6)
    cov2 = t.lapcov .+ (noise .+ transpose(noise)) ./ 2
    r1 = ctsem_dais(t.target, t.mode, t.lapcov; draws=64, seed=22)
    r2 = ctsem_dais(t.target, t.mode, cov2; draws=64, seed=22)
    @test r1.used && r2.used
    @test r1.iterations == r2.iterations
    @test maximum(abs.(r1.mean .- r2.mean) ./ t.sd) < 1e-6
    @test isapprox(r1.cov, r2.cov; rtol=1e-6)
    # Where some directions are left out, the kept span's basis is the
    # span's, whichever orthonormal basis of it, signs included, it comes in.
    V = Matrix(qr(randn(Random.Xoshiro(24), 6, 3)).Q)[:, 1:3]
    turn = Matrix(qr(randn(Random.Xoshiro(25), 3, 3)).Q)[:, 1:3] .* [1.0 -1.0 1.0]
    Q1 = ContinuousTimeSEM._dais_span_basis(V)
    Q2 = ContinuousTimeSEM._dais_span_basis(V * turn)
    @test Q1 ≈ Q2 atol = 1e-12
    @test Q1 * transpose(Q1) ≈ V * transpose(V) atol = 1e-12
    @test transpose(Q1) * Q1 ≈ Matrix(1.0I, 3, 3) atol = 1e-12
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

@testset "DAIS draws on the pool what it draws one point after another" begin
    # On the Laplace route each worker task evaluates its draws on a copy of
    # the objective of its own. Written as the region's closure, the task's
    # locals were the enclosing function's, so the tasks shared one variable:
    # a task could run on another's copy, racing on the units' modes, or find
    # `nothing` there and stop the fit.
    starts, times, data = _laplace_test_data(40, 1; seed=3)
    prior = ctsem_objective(_LAPLACE_LINEAR_OBJECTIVE.params, starts, times, data;
        prior_index=collect(1:8), prior_scale=ones(8))
    laplace = ctsem_laplace_objective(prior, [1, 2], [6, 7], [8], [1.0, 1.0])
    x0 = copy(_LAPLACE_LINEAR_VALUES)
    C = inv(Symmetric(-ContinuousTimeSEM.ctsem_hessian(laplace, x0)))
    one = ctsem_dais(laplace, x0, C; draws=32, maxiter=8, seed=22, parallel=false)
    @test one.used
    @test one.parallel == "serial"
    if Threads.nthreads() > 1 && ContinuousTimeSEM._laplace_pool_width() > 1
        pool = ctsem_dais(laplace, x0, C; draws=32, maxiter=8, seed=22, parallel=true)
        @test pool.parallel == "pool"
        @test pool.iterations == one.iterations
        # The same draws and the same terms; only the order of sums over a
        # unit's members can differ.
        @test isapprox(pool.mean, one.mean; rtol=1e-8)
        @test isapprox(pool.cov, one.cov; rtol=1e-8)
    end
end

@testset "DAIS on a coreset of its own lands where the full objective does, for less" begin
    # Its gradient batches run on a coreset once one is expected to cost
    # less, checked every batch, and an iteration in full teaches it for
    # nothing beyond the gradients. Unchecked once trusted, a coreset of 3-6
    # units of 200 chosen from the first draws put the answer 0.9 sd off.
    starts, times, data = _laplace_test_data(40, 1; seed=3)
    prior = ctsem_objective(_LAPLACE_LINEAR_OBJECTIVE.params, starts, times, data;
        prior_index=collect(1:8), prior_scale=ones(8))
    laplace = ctsem_laplace_objective(prior, [1, 2], [6, 7], [8], [1.0, 1.0])
    x0 = copy(_LAPLACE_LINEAR_VALUES)
    C = inv(Symmetric(-ContinuousTimeSEM.ctsem_hessian(laplace, x0)))
    full = ctsem_dais(laplace, x0, C; draws=64, seed=23, parallel=false)
    viacs = ctsem_dais(laplace, x0, C; draws=64, seed=23, coreset=true)
    @test full.used && viacs.used
    @test !full.coreset
    @test viacs.coreset
    @test viacs.parallel == "coreset"
    @test viacs.coreset_chosen >= 1
    # In full, a gradient costs two; the coreset costs less: on three seeds
    # 1,242-1,562 full evaluations against 3,968-4,480.
    @test full.cost == 2 * full.gradients
    @test viacs.cost < 0.7 * full.cost
    # The same normal, within two runs' Monte Carlo error (each mean to 0.05
    # sd) and the coreset's (0.1 log units RMS at the draws): 0.04-0.05 sd on
    # three seeds.
    sdf = sqrt.(diag(full.cov))
    @test maximum(abs.(viacs.mean .- full.mean) ./ sdf) < 0.3
    @test all(0.75 .< sqrt.(diag(viacs.cov)) ./ sdf .< 1.33)
end

@testset "DAIS reports to the caller's sink" begin
    t = _dais_skew(6, [0.5, 1.0])
    heard = String[]
    ctsem_dais(t.target, t.mode, t.lapcov; draws=64, seed=21, maxiter=5, verbose=true,
        progress_sink=(text, kind) -> (push!(heard, String(text)); nothing))
    @test any(occursin("DAIS", h) for h in heard)
end
