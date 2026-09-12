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


# The gradient bar and the objective floor were one quantity, `1e-6 * max(1,
# |value|)`, used both as a gradient threshold and as an objective one. These
# are about what survives a change of scale, because that is where the old one
# failed and no measurement on a single model could have said so.
@testset "how big the objective is, not what it equals" begin
    scale = ContinuousTimeSEM._ctsem_objective_scale

    # The ordinary case: every subject contributes the same sign, and this is
    # the total's magnitude. That is what keeps it the same bar the tolerances
    # were measured against.
    @test scale(fill(-37.0, 40)) ≈ 1480.0
    @test scale(fill(37.0, 40)) ≈ 1480.0

    # The cases `abs(sum)` gets wrong. A log likelihood is a sum of log
    # *densities*, positive wherever a density exceeds one, so contributions can
    # cancel: forty subjects whose total is zero still have forty subjects'
    # worth of data, and a bar built on the total would collapse to its floor.
    mixed = vcat(fill(37.0, 20), fill(-37.0, 20))
    @test abs(sum(mixed)) == 0.0
    @test scale(mixed) ≈ 1480.0

    # Bounded below by one, and defined for a route that reports no per-subject
    # split at all.
    @test scale(Float64[]) == 1.0
    @test scale([0.0, 0.0]) == 1.0
end

@testset "the gradient bar asks for a gradient per unit of objective" begin
    bar = ContinuousTimeSEM._ctsem_gradient_tolerance

    # Ten times the data is ten times the objective and ten times the gradient,
    # so the bar grows with it and the same fit is judged the same way. That is
    # what makes one tolerance usable on a 6-subject model and a 600-subject
    # one.
    @test bar(0.0, 10 * 1483.0) ≈ 10 * bar(0.0, 1483.0)

    # Never relative to the gradient's own history, which was tried: a fit that
    # starts somewhere terrible would then be judged by a bar as large as the
    # gradient it started with. Stated as a property of a *bad* fit, because
    # that is the case that matters -- a largest gradient of 2.7e4 at an
    # objective of size 3870 must fail, and under a worst-gradient bar it
    # passed.
    @test 2.7e4 > bar(0.0, 3870.0)

    # Bounded below, so a tiny objective does not get a bar of zero, and the
    # absolute tolerance is a floor that is never overridden downward.
    @test bar(0.0, 0.0) == bar(0.0, 1.0)
    @test bar(1e-3, 1.0) == 1e-3
    @test bar(1e-8, 1483.0) > 1e-8
end

@testset "the objective bar is in objective units, at materiality" begin
    bar = ContinuousTimeSEM._ctsem_objective_tolerance

    # Proportional to the value, so it means the same thing on a log likelihood
    # of -30 and one of -3e6, and blind to the sign.
    @test bar(2e6) ≈ 100 * bar(2e4)
    # Never smaller than the bar for a scale of one, so a tiny objective does
    # not get a tolerance of zero.
    @test bar(0.0) == bar(1.0)
    @test bar(1e-30) == bar(0.5)

    # Between the arithmetic and the finding, with room on both sides: the
    # overstep this guards against was worth 16 log likelihood units and the
    # flat transform it must not flag is worth exactly zero. Deliberately not
    # `sqrt(eps)`, which answers a different question -- a false positive here
    # tells a user a converged fit is not a maximum.
    @test bar(2000.0) > 2000.0 * sqrt(eps(Float64))
    @test bar(2000.0) < 1.0
end
