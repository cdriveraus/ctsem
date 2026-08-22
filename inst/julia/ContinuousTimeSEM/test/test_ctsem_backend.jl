using DataFrames

# Keep this test independent of filename order in runtests.jl.
function _ctsem_backend_parameters(; manifest_variance=0.4, t0_variance=0.3)
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    values = Float64[]
    function add_matrix!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name); push!(rows, i); push!(cols, j)
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
    add_matrix!(:MANIFESTVAR, [manifest_variance;;])
    add_matrix!(:T0VAR, [t0_variance;;])
    add_matrix!(:T0MEANS, reshape([0.0], :, 1))
    add_matrix!(:TDPREDEFFECT, [1.0;;])
    add_matrix!(:Jtd, [1.0;;])
    add_matrix!(:PARS, reshape([0.0], :, 1))
    axis = retrieve_axes(DataFrame(matrix=matrices, row=rows, col=cols))
    n = length(values)
    ContinuousTimeSEM.EKFParameters(falses(n), fill(false, n), fill(false, n),
        fill(false, n), fill(false, n), Function[], Function[], Function[], Function[],
        Int[], axis, fill(true, n), AbstractFloat[values...], Int[], Int[], Int[], [1])
end

function _ctsem_static_state_parameters()
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    values = Float64[]
    function add_matrix!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            push!(matrices, name); push!(rows, i); push!(cols, j)
            push!(values, mat[i, j])
        end
    end
    add_matrix!(:DRIFT, [-0.5 0.0; 0.0 0.0])
    add_matrix!(:JAx, [-0.5 10.0; 0.0 0.0])
    add_matrix!(:CINT, reshape([0.0, 0.0], :, 1))
    add_matrix!(:DIFFUSION, [0.2 0.0; 0.0 0.0])
    add_matrix!(:LAMBDA, [1.0 0.0])
    add_matrix!(:Jy, [1.0 0.0])
    add_matrix!(:MANIFESTMEANS, reshape([0.0], :, 1))
    add_matrix!(:MANIFESTVAR, [0.4;;])
    add_matrix!(:T0VAR, [0.3 0.0; 0.0 3.0])
    add_matrix!(:T0MEANS, reshape([0.0, 0.0], :, 1))
    add_matrix!(:PARS, reshape([0.0], :, 1))
    axis = retrieve_axes(DataFrame(matrix=matrices, row=rows, col=cols))
    n = length(values)
    ContinuousTimeSEM.EKFParameters(falses(n), fill(false, n), fill(false, n),
        fill(false, n), fill(false, n), Function[], Function[], Function[], Function[],
        Int[], axis, fill(true, n), AbstractFloat[values...], Int[], Int[], Int[], [1])
end

@testset "TD impulses apply before each row measurement" begin
    sp = _ctsem_backend_parameters()
    subject_starts = [1]
    times = [0.0, 1.0]
    data = reshape([1.0, 1.0], 1, :)
    no_impulse = ContinuousTimeSEM.ctsem_objective(sp, subject_starts, times, data,
        zeros(0, 2), zeros(1, 0))
    impulse = ContinuousTimeSEM.ctsem_objective(sp, subject_starts, times, data,
        reshape([1.0, 0.0], 1, :), zeros(1, 0))
    @test impulse(Float64[]) > no_impulse(Float64[])
end

@testset "ctsem backend multi-subject objective" begin
    sp = _ctsem_backend_parameters()
    subject_starts = [1, 3]
    times = [0.0, 0.5, 0.0, 0.5, 1.0]
    data = reshape([0.1, -0.2, 0.0, 0.2, -0.1], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, subject_starts, times, data)
    result = ContinuousTimeSEM.ctsem_evaluate(objective, Float64[]; contributions=true)
    @test isfinite(result.value)
    @test sum(result.subject_loglik) ≈ result.value
    @test sum(result.row_loglik) ≈ result.value
    @test length(result.row_loglik) == length(times)
    @test result.gradient == Float64[]
end

@testset "ctsem backend validates R-prepared subject starts" begin
    sp = _ctsem_backend_parameters()
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_objective(
        sp, [2], [0.0, 0.0, 1.0], reshape([0.0, 0.1, 0.2], 1, :),
    )
end

@testset "static augmented states have finite likelihoods" begin
    sp = _ctsem_static_state_parameters()
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 1.0],
        reshape([0.0, 0.0], 1, :))
    @test isfinite(objective(Float64[]))
end

@testset "adjoint primitives agree with directional finite differences" begin
    A = [-0.4 0.1; 0.0 -0.2]
    E = [0.2 -0.1; 0.05 0.3]
    h = 1e-6
    expected = (exp(A + h * E) - exp(A - h * E)) / (2h)
    @test ContinuousTimeSEM._ctsem_exp_frechet_block(A, E) ≈ expected rtol=1e-5
end
