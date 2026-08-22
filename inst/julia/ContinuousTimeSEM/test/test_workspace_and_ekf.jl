using ComponentArrays
using DataFrames
using ForwardDiff
using LinearAlgebra

# Build the smallest fixed-parameter model that still exercises the structured
# EKF parameter axis and all matrices expected by the filtering code.
function fixed_one_dimensional_ekf_parameters()
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    values = Float64[]

    function add_matrix!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name)
            push!(rows, i)
            push!(cols, j)
            push!(values, mat[i, j])
        end
    end

    add_matrix!(:DRIFT, [-0.5;;])
    add_matrix!(:JAx, [-0.5;;])
    add_matrix!(:CINT, reshape([0.0], :, 1))
    add_matrix!(:DIFFUSION, [0.2;;])
    add_matrix!(:LAMBDA, [1.0;;])
    add_matrix!(:Jy, [1.0;;])
    add_matrix!(:MANIFESTMEANS, reshape([0.0], :, 1))
    add_matrix!(:MANIFESTVAR, [0.4;;])
    add_matrix!(:T0VAR, [0.3;;])
    add_matrix!(:T0MEANS, reshape([0.0], :, 1))
    add_matrix!(:PARS, reshape([0.0], :, 1))

    axis = retrieve_axes(DataFrame(matrix = matrices, row = rows, col = cols))
    n = length(values)

    return ContinuousTimeSEM.EKFParameters(
        falses(n),
        fill(false, n),
        fill(false, n),
        fill(false, n),
        Function[],
        Function[],
        Function[],
        Int[],
        axis,
        fill(true, n),
        AbstractFloat[values...],
    )
end

function mutable_one_dimensional_ekf_parameters()
    sp = fixed_one_dimensional_ekf_parameters()
    all_values = Vector{Float64}(sp.fixed_values)
    mutables = copy(sp.mutables)
    mutables[1] = true   # DRIFT[1, 1]
    mutables[2] = true   # JAx[1, 1]
    fixed_indices = collect(.!mutables)

    return ContinuousTimeSEM.EKFParameters(
        mutables,
        collect(mutables),
        fill(false, length(mutables)),
        fill(false, length(mutables)),
        Function[p -> -p[1], p -> -p[1]],
        Function[],
        Function[],
        Int[],
        sp.parameter_axis,
        fixed_indices,
        AbstractFloat[all_values[fixed_indices]...],
    )
end

function fixed_rectangular_measurement_ekf_parameters()
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    values = Float64[]

    function add_matrix!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name)
            push!(rows, i)
            push!(cols, j)
            push!(values, mat[i, j])
        end
    end

    drift = [-0.5 0.1; 0.0 -0.2]
    add_matrix!(:DRIFT, drift)
    add_matrix!(:JAx, drift)
    add_matrix!(:CINT, reshape([0.0, 0.0], :, 1))
    add_matrix!(:DIFFUSION, [0.2 0.0; 0.0 0.1])
    add_matrix!(:LAMBDA, [1.0 0.5])
    add_matrix!(:Jy, [1.0 0.5])
    add_matrix!(:MANIFESTMEANS, reshape([0.0], :, 1))
    add_matrix!(:MANIFESTVAR, [0.4;;])
    add_matrix!(:T0VAR, [0.3 0.0; 0.0 0.2])
    add_matrix!(:T0MEANS, reshape([0.0, 0.0], :, 1))
    add_matrix!(:PARS, reshape([0.0], :, 1))

    axis = retrieve_axes(DataFrame(matrix = matrices, row = rows, col = cols))
    n = length(values)

    return ContinuousTimeSEM.EKFParameters(
        falses(n),
        fill(false, n),
        fill(false, n),
        fill(false, n),
        Function[],
        Function[],
        Function[],
        Int[],
        axis,
        fill(true, n),
        AbstractFloat[values...],
    )
end

# Workspace caching is important for ForwardDiff-heavy likelihood evaluation;
# this checks dimensions and that compatible calls reuse the same object.
@testset "Continuous EKF workspace" begin
    sp = fixed_one_dimensional_ekf_parameters()
    ws = ContinuousTimeSEM._init_continuous_ekf_workspace(Float64, sp)

    @test length(ws.all_params) == length(sp.mutables)
    @test size(ws.K) == (1, 1)
    @test size(ws.P_predict.data) == (1, 1)

    ws_ref = Ref{Any}(nothing)
    cached = ContinuousTimeSEM._get_or_init_continuous_ekf_workspace!(ws_ref, Float64[], sp)
    @test cached === ws_ref[]
    @test ContinuousTimeSEM._get_or_init_continuous_ekf_workspace!(ws_ref, Float64[], sp) === cached
end

# A tiny fixed-parameter likelihood path catches integration breakage across
# covariance transforms, discretization, update steps, and log-likelihood code.
@testset "Continuous EKF fixed-parameter likelihood" begin
    sp = fixed_one_dimensional_ekf_parameters()
    data = reshape([0.1, -0.2, 0.05], 1, :)
    timesteps = [0.0, 0.5, 1.5]

    ll = ContinuousTimeSEM.extended_kalman_log_likelihood_continuous(Float64[], sp, timesteps, data)
    @test isfinite(ll)

    trajectory_sp = mutable_one_dimensional_ekf_parameters()
    trajectory_values = [0.5]
    workspace = ContinuousTimeSEM.cts_benchmark_subject_workspace(
        trajectory_values,
        trajectory_sp,
        [timesteps],
        [data],
    )
    trajectory = ContinuousTimeSEM.cts_benchmark_predict_trajectory(
        trajectory_values,
        workspace,
        timesteps,
    )
    @test size(trajectory.states) == (1, length(timesteps))
    @test size(trajectory.fitted_manifest) == size(data)
    @test all(isfinite, trajectory.states)
    @test all(isfinite, trajectory.fitted_manifest)

    substep_workspace = ContinuousTimeSEM.cts_benchmark_substep_workspace(
        trajectory_values,
        workspace,
        0.05,
    )
    substep_result = ContinuousTimeSEM.cts_benchmark_res_and_grad_subjects_cached(
        substep_workspace,
        trajectory_values,
    )
    @test isfinite(substep_result.value)
    @test all(isfinite, substep_result.gradient)
end

@testset "Continuous EKF input validation" begin
    sp = fixed_one_dimensional_ekf_parameters()
    data = reshape([0.1, -0.2, 0.05], 1, :)

    decreasing_time_error = try
        ContinuousTimeSEM.ContinuousEKFGradientWorkspace(
            Float64[],
            sp,
            [0.0, 1.0, 0.0],
            data,
        )
        nothing
    catch err
        err
    end
    @test decreasing_time_error isa ArgumentError
    @test occursin("separate subjects", sprint(showerror, decreasing_time_error))

    @test_throws DimensionMismatch ContinuousTimeSEM.extended_kalman_log_likelihood_continuous(
        Float64[],
        sp,
        [0.0, 1.0],
        data,
    )
end

# Rectangular measurement matrices occur when a model has more latent states
# than manifest variables, e.g. intervention accumulator states.
@testset "Continuous EKF rectangular measurement likelihood" begin
    sp = fixed_rectangular_measurement_ekf_parameters()
    data = reshape([0.1, -0.2, 0.05], 1, :)
    timesteps = [0.0, 0.5, 1.5]

    ll = ContinuousTimeSEM.extended_kalman_log_likelihood_continuous(Float64[], sp, timesteps, data)
    @test isfinite(ll)
end

@testset "Continuous EKF cached gradient workspace" begin
    sp = mutable_one_dimensional_ekf_parameters()
    data = reshape([0.1, -0.2, 0.05], 1, :)
    timesteps = [0.0, 0.5, 1.5]
    values = [0.5]

    workspace = ContinuousTimeSEM.ContinuousEKFGradientWorkspace(values, sp, timesteps, data)
    result = zeros(length(values))

    @test @inferred(workspace.objective(values)) isa Float64
    dual_values = ForwardDiff.Dual.(values, one.(values))
    @test @inferred(workspace.objective(dual_values)) isa ForwardDiff.Dual

    cached_gradient = ContinuousTimeSEM.grad_log_likelihood_ekf_continuous!(result, values, workspace)
    public_gradient = ContinuousTimeSEM.grad_log_likelihood_ekf_continuous(values, sp, timesteps, data)

    @test cached_gradient === result
    @test result ≈ public_gradient

    workspace_gradient = ContinuousTimeSEM.grad_log_likelihood_ekf_continuous!(values, workspace)
    @test workspace_gradient === workspace.result
    @test workspace_gradient ≈ public_gradient

    cached_res = ContinuousTimeSEM.res_and_grad_likelihood_ekf_continuous!(workspace, values)
    public_res = ContinuousTimeSEM.res_and_grad_likelihood_ekf_continuous(values, sp, timesteps, data)

    @test cached_res.value ≈ public_res.value
    @test cached_res.gradient ≈ public_res.gradient
end
