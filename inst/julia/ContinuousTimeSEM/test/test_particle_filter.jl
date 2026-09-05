using LinearAlgebra
using DataFrames
using Random

# The particle filter is the reference the filter is checked against, so it is
# itself checked where an exact answer exists: a linear Gaussian model, where
# the exponential particle step is exact and the EKF likelihood is the truth;
# and a one-latent linear process with any measurement density, where a
# Chapman-Kolmogorov recursion on a fine grid is the truth and the EKF is not.
# Everything else is consistency: two transitions on a fine mesh agree, the
# per-row increments sum to the total, a seed reproduces.

function _pf_dataframe(; drift, jax, cint, diffusion, lambda, jy, manifestmeans,
    manifestvar, t0var, t0means, pars=zeros(1, 1), free=Dict(), predict=Dict(),
    update=Dict())
    matrices = Symbol[]; rows = Int[]; cols = Int[]
    parnumber = Union{Missing,Int}[]; value = Union{Missing,Float64}[]
    transform = Union{Missing,String}[]; predicttransform = Union{Missing,String}[]
    updatetransform = Union{Missing,String}[]
    function addmat!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            key = (name, i, j)
            push!(matrices, name); push!(rows, i); push!(cols, j)
            if haskey(free, key)
                pn, tf = free[key]
                push!(parnumber, pn); push!(value, missing); push!(transform, tf)
            else
                push!(parnumber, missing); push!(value, Float64(mat[i, j])); push!(transform, missing)
            end
            push!(predicttransform, get(predict, key, missing))
            push!(updatetransform, get(update, key, missing))
        end
    end
    addmat!(:DRIFT, drift); addmat!(:JAx, jax); addmat!(:CINT, cint)
    addmat!(:DIFFUSION, diffusion); addmat!(:LAMBDA, lambda); addmat!(:Jy, jy)
    addmat!(:MANIFESTMEANS, manifestmeans); addmat!(:MANIFESTVAR, manifestvar)
    addmat!(:T0VAR, t0var); addmat!(:T0MEANS, t0means); addmat!(:PARS, pars)
    DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber, value=value,
        transform=transform, predicttransform=predicttransform, updatetransform=updatetransform)
end

# The shared `ekf_from_data_frame` has no indicator-type arguments.
function _pf_spec(df; manifesttype=Int[], ncategories=Int[])
    ContinuousTimeSEM.ekf_from_columns(df.matrix, df.row, df.col,
        _tbl_int(df.parnumber), _tbl_float(df.value), _tbl_string(df.transform),
        _tbl_string(df.predicttransform), _tbl_string(df.updatetransform),
        fill("", nrow(df)); manifesttype=manifesttype, ncategories=ncategories)
end

# The variance a standard deviation becomes, from the transform rather than by
# squaring: `sdcovsqrt2cov` carries a 1e-5 that does not cancel.
function _pf_variance(sd)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 1)
    ContinuousTimeSEM.sdcovsqrt2cov!(buffer, [sd;;], 0, Val(1))
    return buffer.out[1, 1]
end

# Exact, to grid resolution, marginal log likelihood of one subject of a
# one-latent linear process under any measurement density: the
# Chapman-Kolmogorov recursion on a grid. `rowlogpdf(x, t)` is log p(y_t | x)
# over the grid vector `x` for the subject's row `t`. The transition is the
# exact Gaussian kernel of the linear SDE, so this shares no approximation with
# the EKF's measurement update, which is what the tests below need it for.
function _pf_grid_loglik(rowlogpdf, times; drift, cint, diffusion_sd, t0mean, t0sd,
    grid=range(-8.0, 8.0; length=1601))
    q = _pf_variance(diffusion_sd)
    x = collect(grid)
    p = exp.(-(x .- t0mean) .^ 2 ./ (2 * _pf_variance(t0sd)))
    p ./= sum(p)
    kernels = Dict{Float64,Matrix{Float64}}()
    total = 0.0
    for t in eachindex(times)
        if t > 1
            dt = times[t] - times[t - 1]
            K = get!(kernels, dt) do
                a = exp(drift * dt)
                b = cint * (a - 1) / drift
                qd = q * (a^2 - 1) / (2 * drift)
                K = [exp(-(xi - a * xj - b)^2 / (2 * qd)) for xi in x, xj in x]
                K ./ sum(K; dims=1)
            end
            p = K * p
        end
        l = exp.(rowlogpdf(x, t))
        s = sum(p .* l)
        total += log(s)
        p .*= l ./ s
    end
    return total
end

_PF_LINEAR = ekf_from_data_frame(_pf_dataframe(
    drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3], cint=[0.1; 0.0;;],
    diffusion=[0.4 0.0; 0.1 0.3], lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
    manifestmeans=[0.0; 0.2;;], manifestvar=[0.3 0.0; 0.0 0.4],
    t0var=[1.0 0.0; 0.2 1.0], t0means=[0.5; 0.0;;],
    free=Dict((:DRIFT, 1, 2) => (1, "param[1]"), (:JAx, 1, 2) => (1, "param[1]"))))

# Cubic damping, stable everywhere: f(x) = a x (1 + c x^2).
_PF_CUBIC = ekf_from_data_frame(_pf_dataframe(
    drift=[0.0;;], jax=[0.0;;], cint=[0.0;;], diffusion=[0.5;;], lambda=[1.0;;], jy=[1.0;;],
    manifestmeans=[0.0;;], manifestvar=[0.3;;], t0var=[1.0;;], t0means=[2.0;;],
    pars=[0.5;;],
    free=Dict((:PARS, 1, 1) => (1, "log1p_exp(param[1])")),
    predict=Dict((:DRIFT, 1, 1) => "-0.6 * (1 + PARS[1,1] * state[1]^2)",
        (:JAx, 1, 1) => "-0.6 * (1 + 3 * PARS[1,1] * state[1]^2)")))

# Two latents, the cross effect of the first on the second weakening with the
# first's square: f_2 = 0.2 (1 - c x_1^2) x_1 - 0.4 x_2. The strength `c` is
# the free parameter, so c = 0 is the linear model and the EKF is exact there.
_PF_CROSS = ekf_from_data_frame(_pf_dataframe(
    drift=[-0.5 0.3; 0.2 -0.4], jax=[-0.5 0.3; 0.2 -0.4], cint=[0.1; 0.0;;],
    diffusion=[0.4 0.0; 0.1 0.3], lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
    manifestmeans=[0.0; 0.2;;], manifestvar=[0.3 0.0; 0.0 0.4],
    t0var=[1.0 0.0; 0.2 1.0], t0means=[0.5; 0.0;;],
    pars=[0.0;;],
    free=Dict((:PARS, 1, 1) => (1, "param[1]")),
    predict=Dict((:DRIFT, 2, 1) => "0.2 * (1 - PARS[1,1] * state[1]^2)",
        (:JAx, 2, 1) => "0.2 * (1 - 3 * PARS[1,1] * state[1]^2)")))

# One latent, a binary and a count indicator, nothing Gaussian to measure it.
_PF_CATEGORICAL = _pf_spec(_pf_dataframe(
    drift=[-0.5;;], jax=[-0.5;;], cint=[0.1;;], diffusion=[0.5;;],
    lambda=[1.0; 0.8;;], jy=[1.0; 0.8;;], manifestmeans=[0.2; 0.5;;],
    manifestvar=zeros(2, 2), t0var=[1.0;;], t0means=[0.3;;],
    free=Dict((:DRIFT, 1, 1) => (1, "param[1]"), (:JAx, 1, 1) => (1, "param[1]")));
    manifesttype=[ContinuousTimeSEM.CTSEM_OBS_BINARY, ContinuousTimeSEM.CTSEM_OBS_COUNT],
    ncategories=[0, 0])

# One latent measured through a loading that depends on the state it loads on:
# y = 0.1 + (1 + c x) x + e. `c` free, so c = 0 is again the exact linear case.
_PF_LOADING = ekf_from_data_frame(_pf_dataframe(
    drift=[-0.6;;], jax=[-0.6;;], cint=[0.0;;], diffusion=[0.5;;], lambda=[1.0;;], jy=[1.0;;],
    manifestmeans=[0.1;;], manifestvar=[0.3;;], t0var=[1.0;;], t0means=[0.5;;],
    pars=[0.0;;],
    free=Dict((:PARS, 1, 1) => (1, "param[1]")),
    update=Dict((:LAMBDA, 1, 1) => "1 + PARS[1,1] * state[1]",
        (:Jy, 1, 1) => "1 + 2 * PARS[1,1] * state[1]")))

_pf_times(nsub, nrow, dt) = ([1 + (s - 1) * nrow for s in 1:nsub],
    repeat(collect(0.0:dt:(dt * (nrow - 1))), nsub))

@testset "particle filter reproduces the exact likelihood of a linear model" begin
    starts, times = _pf_times(4, 6, 1.0)
    data = 0.5 .* randn(MersenneTwister(1), 2, length(times))
    obj = ContinuousTimeSEM.ctsem_objective(_PF_LINEAR, starts, times, data)
    θ = [0.2]
    exact = obj(θ)
    # One exponential step per interval is the exact transition for a linear
    # model, so the only error is Monte Carlo.
    pf = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=1, seed=3)
    @test isfinite(pf.loglik) && pf.se < 0.5
    @test abs(pf.loglik - exact) < 4 * pf.se + 0.05
    @test length(pf.row_loglik) == length(times)
    @test sum(pf.row_loglik) ≈ pf.loglik atol=1e-9
    @test pf.ess_min > 0
    # More substeps change nothing but the draws for a linear model.
    pf8 = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=8, seed=3)
    @test abs(pf8.loglik - exact) < 4 * pf8.se + 0.05
    # A seed reproduces; another seed is another draw of the same quantity.
    again = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=1, seed=3)
    @test again.loglik == pf.loglik
    other = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=1, seed=4)
    @test other.loglik != pf.loglik
    @test abs(other.loglik - exact) < 4 * other.se + 0.05
end

@testset "the two particle transitions agree on a nonlinear model at a fine mesh" begin
    starts, times = _pf_times(5, 6, 1.0)
    rng = MersenneTwister(2)
    data = reshape(2.0 .* exp.(-0.5 .* repeat(0:5, 5)) .+ 0.3 .* randn(rng, 30), 1, :)
    obj = ContinuousTimeSEM.ctsem_objective(_PF_CUBIC, starts, times, data)
    θ = [-0.5]
    expo = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=20, seed=5)
    euler = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=80,
        transition=:euler, seed=5)
    @test isfinite(expo.loglik) && isfinite(euler.loglik)
    @test abs(expo.loglik - euler.loglik) < 4 * (expo.se + euler.se) + 0.2
    @test expo.transition == :exponential && euler.transition == :euler
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; transition=:bogus)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=1)
end

@testset "multivariate state-dependent drift: exact at zero strength, consistent otherwise" begin
    starts, times = _pf_times(4, 6, 1.0)
    data = 0.5 .* randn(MersenneTwister(8), 2, length(times))
    obj = ContinuousTimeSEM.ctsem_objective(_PF_CROSS, starts, times, data)
    # The state-dependent cells are present and evaluated per particle, but at
    # zero strength they equal the linear model's, so the EKF is exact.
    exact = obj([0.0])
    pf = ContinuousTimeSEM.ctsem_particle_loglik(obj, [0.0]; particles=4000, substeps=1, seed=9)
    @test isfinite(exact) && isfinite(pf.loglik)
    @test abs(pf.loglik - exact) < 4 * pf.se + 0.05
    # At nonzero strength the two transitions, which share nothing over a step,
    # agree on a fine mesh; the exponential one needs far fewer substeps.
    θ = [0.3]
    expo = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=2000, substeps=8, seed=10)
    euler = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=2000, substeps=40,
        transition=:euler, seed=10)
    @test isfinite(expo.loglik) && isfinite(euler.loglik)
    @test abs(expo.loglik - euler.loglik) < 4 * (expo.se + euler.se) + 0.2
    @test sum(expo.row_loglik) ≈ expo.loglik atol=1e-9
    @test length(expo.subject_loglik) == 4
end

@testset "binary and count indicators: the particle likelihood matches a grid filter" begin
    nsub, nrow = 3, 5
    starts, times = _pf_times(nsub, nrow, 1.0)
    rng = MersenneTwister(11)
    data = vcat(Float64.(rand(rng, 0:1, 1, nsub * nrow)), Float64.(rand(rng, 0:4, 1, nsub * nrow)))
    obj = ContinuousTimeSEM.ctsem_objective(_PF_CATEGORICAL, starts, times, data)
    θ = [-0.5]
    # The same densities the engine uses: logistic for the binary indicator,
    # Poisson with a log link for the count.
    grid = sum(1:nsub) do s
        cols = starts[s]:(starts[s] + nrow - 1)
        _pf_grid_loglik(times[cols]; drift=-0.5, cint=0.1, diffusion_sd=0.5,
            t0mean=0.3, t0sd=1.0) do x, t
            y1, y2 = data[1, cols[t]], data[2, cols[t]]
            η1 = 0.2 .+ x
            η2 = 0.5 .+ 0.8 .* x
            binary = y1 > 0.5 ? -log1p.(exp.(-η1)) : -log1p.(exp.(η1))
            count = y2 .* η2 .- exp.(η2) .- sum(log, 1:Int(y2); init=0.0)
            binary .+ count
        end
    end
    @test isfinite(grid)
    pf = ContinuousTimeSEM.ctsem_particle_loglik(obj, θ; particles=4000, substeps=1, seed=12)
    @test isfinite(pf.loglik) && pf.ess_min > 0
    @test abs(pf.loglik - grid) < 4 * pf.se + 0.05
    @test sum(pf.row_loglik) ≈ pf.loglik atol=1e-9
    # The EKF's categorical update is an approximation here; the point of the
    # tool is that the particle estimate is not.
    @test isfinite(obj(θ))
end

@testset "state-dependent loading: exact at zero strength, matches the grid filter otherwise" begin
    nsub, nrow = 3, 5
    starts, times = _pf_times(nsub, nrow, 1.0)
    data = 0.6 .* randn(MersenneTwister(13), 1, nsub * nrow) .+ 0.3
    obj = ContinuousTimeSEM.ctsem_objective(_PF_LOADING, starts, times, data)
    exact = obj([0.0])
    pf0 = ContinuousTimeSEM.ctsem_particle_loglik(obj, [0.0]; particles=4000, substeps=1, seed=14)
    @test abs(pf0.loglik - exact) < 4 * pf0.se + 0.05
    c = 0.3
    v = _pf_variance(0.3)
    grid = sum(1:nsub) do s
        cols = starts[s]:(starts[s] + nrow - 1)
        _pf_grid_loglik(times[cols]; drift=-0.6, cint=0.0, diffusion_sd=0.5,
            t0mean=0.5, t0sd=1.0) do x, t
            μ = 0.1 .+ (1 .+ c .* x) .* x
            -(data[1, cols[t]] .- μ) .^ 2 ./ (2v) .- 0.5 * log(2π * v)
        end
    end
    pf = ContinuousTimeSEM.ctsem_particle_loglik(obj, [c]; particles=4000, substeps=1, seed=15)
    @test isfinite(pf.loglik) && isfinite(grid)
    @test abs(pf.loglik - grid) < 4 * pf.se + 0.05
    # And the loading really is state dependent: the filter's answer differs.
    @test obj([c]) != exact
end

@testset "state-explicit path takes an Euler transition" begin
    starts, times = _pf_times(3, 4, 1.0)
    rng = MersenneTwister(6)
    data = reshape(2.0 .* exp.(-0.5 .* repeat(0:3, 3)) .+ 0.3 .* randn(rng, 12), 1, :)
    obj = ContinuousTimeSEM.ctsem_objective(_PF_CUBIC, starts, times, data,
        zeros(0, length(times)), zeros(3, 0), 0.25)
    npar = 1
    joint_exp = ContinuousTimeSEM.ctsem_joint_objective(obj, npar)
    joint_eul = ContinuousTimeSEM.ctsem_joint_objective(obj, npar; transition="euler")
    @test joint_exp.transition == :exponential && joint_eul.transition == :euler
    @test ContinuousTimeSEM.ctsem_joint_dimension(joint_exp) == ContinuousTimeSEM.ctsem_joint_dimension(joint_eul)
    x = vcat([-0.5], 0.3 .* randn(MersenneTwister(7), ContinuousTimeSEM.ctsem_joint_dimension(joint_exp) - npar))
    a = joint_exp(x); b = joint_eul(x)
    @test isfinite(a) && isfinite(b) && a != b
    # The same joint density through the exported function that takes z directly
    # (both include the innovation prior; this model has no parameter prior).
    z = x[(npar + 1):end]
    @test ContinuousTimeSEM.ctsem_joint_loglikelihood(obj, [-0.5], z; transition=:euler) ≈ b atol=1e-9
    @test ContinuousTimeSEM.ctsem_joint_loglikelihood(obj, [-0.5], z) ≈ a atol=1e-9
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_joint_objective(obj, npar; transition=:bogus)
end

@testset "batch evaluation matches single runs and separates the prior from the likelihood" begin
    starts, times = _pf_times(3, 4, 1.0)
    data = 0.5 .* randn(MersenneTwister(21), 2, length(times))
    obj = ContinuousTimeSEM.ctsem_objective(_PF_LINEAR, starts, times, data)
    thetas = [0.1 0.3]                      # one parameter, two columns
    b = ContinuousTimeSEM.ctsem_particle_batch(obj, thetas; particles=500, substeps=1, seed=4)
    @test length(b.particle) == 2
    for j in 1:2
        single = ContinuousTimeSEM.ctsem_particle_loglik(obj, thetas[:, j]; particles=500,
            substeps=1, seed=4)
        @test b.particle[j] == single.loglik
        @test b.se[j] == single.se
        @test b.ess_min[j] == single.ess_min
        @test b.filter[j] ≈ obj(thetas[:, j]) atol=1e-10      # no prior on this model
        @test b.posterior[j] == b.filter[j]
    end
    # With a prior the posterior column carries it and the filter column does not.
    objp = ContinuousTimeSEM.ctsem_objective(_PF_LINEAR, starts, times, data;
        prior_index=[1], prior_scale=[2.0])
    bp = ContinuousTimeSEM.ctsem_particle_batch(objp, thetas; particles=200, substeps=1, seed=4)
    for j in 1:2
        @test bp.filter[j] ≈ b.filter[j] atol=1e-10
        @test bp.posterior[j] - bp.filter[j] ≈
            ContinuousTimeSEM._ctsem_log_prior(objp, thetas[:, j]) atol=1e-10
        @test bp.posterior[j] != bp.filter[j]
    end
    # One seed per column: the second column then equals a single run at its seed.
    bs = ContinuousTimeSEM.ctsem_particle_batch(obj, thetas; particles=500, substeps=1, seed=[4, 5])
    @test bs.particle[1] == b.particle[1]
    @test bs.particle[2] == ContinuousTimeSEM.ctsem_particle_loglik(obj, thetas[:, 2];
        particles=500, substeps=1, seed=5).loglik
    @test bs.particle[2] != b.particle[2]
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_particle_batch(obj, thetas; seed=[4, 5, 6])
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_particle_batch(obj, zeros(1, 0))
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_particle_loglik(obj, [0.1]; seed=-1)
    # Splitting the subjects across chunks does not change the answer: every
    # subject has its own stream. (Vacuous on one thread, real on several.)
    previous = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    ContinuousTimeSEM.ctsem_set_max_chunks!(2)
    two = ContinuousTimeSEM.ctsem_particle_loglik(obj, [0.1]; particles=500, substeps=1, seed=4)
    ContinuousTimeSEM.ctsem_set_max_chunks!(previous)
    one = ContinuousTimeSEM.ctsem_particle_loglik(obj, [0.1]; particles=500, substeps=1, seed=4)
    @test two.loglik == one.loglik
    @test two.row_loglik == one.row_loglik
end
