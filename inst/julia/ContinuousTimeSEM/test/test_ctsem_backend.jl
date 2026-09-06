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

# `_ctsem_overshot` decides which half of "saturated" a fit is in, and
# `converged` is keyed on its answer. A mock objective is used rather than a
# fitted model because the two cases differ only in the shape of the objective
# along the flagged coordinate, and that is exactly what a mock can state.
struct _OvershotMock
    peak::Float64
end
ContinuousTimeSEM.ctsem_evaluate(m::_OvershotMock, x::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint) =
    (value = -sum(abs2, x .- m.peak), gradient = nothing)

@testset "a saturated coordinate that is not a maximum is an overshoot" begin
    # The optimizer overstepped: the objective peaks at 1 and the reported
    # estimate is at 20, where the transform is flat. Pulling the coordinate
    # back improves the objective, so this is not a maximum. This is the shape
    # measured on a binary model -- one L-BFGS iteration to raw 20.9, log
    # likelihood 16 units below the profile peak.
    mock = _OvershotMock(1.0)
    value = -sum(abs2, [20.0] .- 1.0)
    out = ContinuousTimeSEM._ctsem_overshot(mock, [20.0], [1], value, 1e-3)
    @test out.overshot
    @test out.gain > 100

    # And the other half: the data do not identify the coordinate, so the
    # optimizer ran it to the edge from a maximum that really is there. Nothing
    # a pullback can do improves it, so the fit converged.
    atpeak = _OvershotMock(0.0)
    out2 = ContinuousTimeSEM._ctsem_overshot(atpeak, [0.0], [1],
        -sum(abs2, [0.0]), 1e-3)
    @test !out2.overshot
    @test out2.gain == 0.0

    # Nothing flagged means nothing evaluated and nothing claimed.
    out3 = ContinuousTimeSEM._ctsem_overshot(mock, [20.0], Int[], value, 1e-3)
    @test !out3.overshot
    @test out3.gain == 0.0
end
