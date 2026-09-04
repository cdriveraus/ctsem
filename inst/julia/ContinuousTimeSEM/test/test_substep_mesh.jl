using LinearAlgebra
using DataFrames
using Random

# A substep mesh -- one count per row -- travels in the same slot as the
# `maxtimestep` rule. These tests pin that the two agree where they must, that
# the automatic mesh refines exactly the models that need it, and that every
# consumer of the objective (adjoint, state-explicit dimension) honours a mesh.

function _mesh_test_dataframe(; drift, jax, cint, diffusion, lambda, jy, manifestmeans,
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

# Linear two-latent model, one cross effect free.
_MESH_LINEAR = ekf_from_data_frame(_mesh_test_dataframe(
    drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3], cint=[0.0; 0.0;;],
    diffusion=[0.2 0.0; 0.0 0.15], lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
    manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
    t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
    free=Dict((:DRIFT, 1, 2) => (1, "param[1]"), (:JAx, 1, 2) => (1, "param[1]"))))

# One latent whose drift strengthens with the state: f(x) = p (1 + c x) x, so the
# Jacobian p (1 + 2 c x) moves as the state does and one step per interval is
# not exact. `c` is PARS[1,1].
_MESH_NONLINEAR = ekf_from_data_frame(_mesh_test_dataframe(
    drift=[0.0;;], jax=[0.0;;], cint=[0.0;;], diffusion=[0.5;;], lambda=[1.0;;], jy=[1.0;;],
    manifestmeans=[0.0;;], manifestvar=[0.05;;], t0var=[1.0;;], t0means=[2.0;;],
    pars=[0.6;;],
    free=Dict((:PARS, 1, 1) => (1, "param[1]")),
    predict=Dict((:DRIFT, 1, 1) => "-0.8 * (1 + PARS[1,1] * state[1])",
        (:JAx, 1, 1) => "-0.8 * (1 + 2 * PARS[1,1] * state[1])")))

_mesh_times(nsub, nrow, dt) = ([1 + (s - 1) * nrow for s in 1:nsub],
    repeat(collect(0.0:dt:(dt * (nrow - 1))), nsub))

@testset "substep rule and mesh agree where they should" begin
    @test ContinuousTimeSEM._ctsem_substeps(1.0, Inf, 3) == 1
    @test ContinuousTimeSEM._ctsem_substeps(1.0, 0.3, 3) == 4
    @test ContinuousTimeSEM._ctsem_substeps(1.0, [1, 7, 2], 2) == 7
    starts, times = _mesh_times(3, 5, 1.0)
    data = 0.3 .* randn(MersenneTwister(1), 2, length(times))
    by_rule = ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts, times, data,
        zeros(0, length(times)), zeros(3, 0), 0.5)
    mesh = fill(2, length(times))
    by_mesh = ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts, times, data,
        zeros(0, length(times)), zeros(3, 0), mesh)
    @test by_rule([0.2]) ≈ by_mesh([0.2]) atol=1e-12
    # A linear model is exact at any step, so the mesh can move the answer only
    # through the Stan-matching 1e-10 ridge the filter adds to the covariance
    # before *each* propagation: two substeps ridge twice. Measured at 1.4e-7 on
    # this model; the same holds for `maxtimestep`, and is not new here.
    one_step = ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts, times, data)
    @test one_step([0.2]) ≈ by_mesh([0.2]) atol=1e-5
    # Length and positivity are checked.
    @test_throws DimensionMismatch ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts,
        times, data, zeros(0, length(times)), zeros(3, 0), fill(2, length(times) - 1))
    bad = fill(2, length(times)); bad[3] = 0
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts,
        times, data, zeros(0, length(times)), zeros(3, 0), bad)
end

@testset "a linear model gets the floor without a filter pass" begin
    starts, times = _mesh_times(4, 6, 0.7)
    data = 0.3 .* randn(MersenneTwister(2), 2, length(times))
    obj = ContinuousTimeSEM.ctsem_objective(_MESH_LINEAR, starts, times, data)
    out = ContinuousTimeSEM.ctsem_auto_substeps(obj, [0.2]; tol=1e-6)
    @test out.finite && out.passes == 0
    @test out.refined == 0 && out.total == out.intervals == 4 * 5
    @test all(==(1), out.mesh)
    # With a maxtimestep floor the mesh is the rule's counts.
    out2 = ContinuousTimeSEM.ctsem_auto_substeps(obj, [0.2]; tol=1e-6, floor_rule=0.3)
    @test all(==(3), out2.mesh[2:6]) && out2.refined == 0
end

@testset "a nonlinear model is refined where the state is large, and converges" begin
    starts, times = _mesh_times(20, 8, 1.0)
    rng = MersenneTwister(3)
    # States start near 2 (T0MEANS) and decay; data that keep some subjects high
    # and let others fall so the indicator varies across rows.
    data = reshape(vcat([2.0 .* exp.(-0.3 .* (0:7)) .+ 0.2 .* randn(rng, 8) .+ (s % 3 == 0 ? 1.5 : 0.0)
                          for s in 1:20]...), 1, :)
    obj = ContinuousTimeSEM.ctsem_objective(_MESH_NONLINEAR, starts, times, data)
    θ = [0.6]
    loose = ContinuousTimeSEM.ctsem_auto_substeps(obj, θ; tol=0.5)
    tight = ContinuousTimeSEM.ctsem_auto_substeps(obj, θ; tol=0.01)
    @test loose.finite && tight.finite
    @test tight.refined > 0
    @test tight.total > loose.total
    @test tight.max_substeps <= 64
    # Each tightening of the tolerance brings the likelihood closer to a fine
    # mesh: the mechanism, not a calibration -- the tolerance's meaning in
    # likelihood units is what review/INTEGRATOR-substeps measures on simulated
    # data, and it depends on the model.
    with_mesh(mesh) = ContinuousTimeSEM.ctsem_objective(_MESH_NONLINEAR, starts, times, data,
        zeros(0, length(times)), zeros(20, 0), mesh)(θ)
    ll_fine = with_mesh(fill(64, length(times)))
    err_one = abs(obj(θ) - ll_fine)
    err_loose = abs(with_mesh(loose.mesh) - ll_fine)
    err_tight = abs(with_mesh(tight.mesh) - ll_fine)
    @test err_tight < err_loose < err_one
    @test err_tight < 0.25 * err_one
    # A mesh at the tolerance re-measures as within tolerance (fixed point).
    again = ContinuousTimeSEM.ctsem_auto_substeps(
        ContinuousTimeSEM.ctsem_objective(_MESH_NONLINEAR, starts, times, data,
            zeros(0, length(times)), zeros(20, 0), tight.mesh), θ; tol=0.01)
    @test again.passes >= 1
end

@testset "adjoint and state-explicit dimension honour a mesh" begin
    starts, times = _mesh_times(3, 5, 1.0)
    rng = MersenneTwister(4)
    data = reshape(2.0 .* exp.(-0.3 .* repeat(0:4, 3)) .+ 0.2 .* randn(rng, 15), 1, :)
    mesh = ones(Int, length(times))
    mesh[[2, 4, 7, 9, 13]] .= 3
    obj = ContinuousTimeSEM.ctsem_objective(_MESH_NONLINEAR, starts, times, data,
        zeros(0, length(times)), zeros(3, 0), mesh)
    θ = [0.6]
    result = ContinuousTimeSEM.ctsem_validate_forward_gradient(obj, θ)
    @test result.relative_error < 1e-4
    @test result.adjoint_relative_error < 1e-9
    # One innovation per substep: 1 (initial) + sum of counts over intervals per subject.
    expected = 3 * 1 + sum(mesh[t] for t in eachindex(mesh) if (t - 1) % 5 != 0)
    @test ContinuousTimeSEM.ctsem_state_dimension(obj) == expected
end
