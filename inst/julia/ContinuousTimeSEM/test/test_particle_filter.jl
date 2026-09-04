using LinearAlgebra
using DataFrames
using Random

# The particle filter is the reference the filter is checked against, so it is
# itself checked where an exact answer exists: a linear model, where the
# exponential particle step is exact and the EKF likelihood is the truth.
# Everything else is consistency: two transitions on a fine mesh agree, the
# per-row increments sum to the total, a seed reproduces.

function _pf_dataframe(; drift, jax, cint, diffusion, lambda, jy, manifestmeans,
    manifestvar, t0var, t0means, pars=zeros(1, 1), free=Dict(), predict=Dict())
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
            push!(updatetransform, missing)
        end
    end
    addmat!(:DRIFT, drift); addmat!(:JAx, jax); addmat!(:CINT, cint)
    addmat!(:DIFFUSION, diffusion); addmat!(:LAMBDA, lambda); addmat!(:Jy, jy)
    addmat!(:MANIFESTMEANS, manifestmeans); addmat!(:MANIFESTVAR, manifestvar)
    addmat!(:T0VAR, t0var); addmat!(:T0MEANS, t0means); addmat!(:PARS, pars)
    DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber, value=value,
        transform=transform, predicttransform=predicttransform, updatetransform=updatetransform)
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
