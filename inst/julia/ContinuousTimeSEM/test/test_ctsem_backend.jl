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

    # The saturation flag is not what selects the probe. The default probe
    # orders coordinates by magnitude and never reads the flag, so an empty
    # flag changes nothing about what it finds.
    out3 = ContinuousTimeSEM._ctsem_overshot(mock, [20.0], Int[], value, 1e-3)
    @test out3.overshot
    @test out3.coordinates == [1]

    # Under `:saturation` it is the flag, and an empty one means nothing
    # evaluated and nothing claimed -- the behaviour the switch preserves.
    previous = ContinuousTimeSEM.ctsem_set_overshoot_probe!(:saturation)
    try
        @test ContinuousTimeSEM._ctsem_overshot(mock, [20.0], [1], value, 1e-3).overshot
        flagless = ContinuousTimeSEM._ctsem_overshot(mock, [20.0], Int[], value, 1e-3)
        @test !flagless.overshot
        @test flagless.gain == 0.0
    finally
        ContinuousTimeSEM.ctsem_set_overshoot_probe!(previous)
    end

    # And `:off` claims nothing whatever is handed to it.
    previous = ContinuousTimeSEM.ctsem_set_overshoot_probe!(:off)
    try
        @test !ContinuousTimeSEM._ctsem_overshot(mock, [20.0], [1], value, 1e-3).overshot
    finally
        ContinuousTimeSEM.ctsem_set_overshoot_probe!(previous)
    end
    @test ContinuousTimeSEM._CTSEM_OVERSHOOT_PROBE[] === :magnitude
end

# The reason the probe is not one coordinate at a time. A degenerate corner is
# left by moving a whole block together, and every single-coordinate move out of
# one is worse than staying -- which is what a collapsed population scale and
# its correlations look like, and what made the old probe report those fits as
# maxima. See `_ctsem_overshot`.
struct _JointMock end
ContinuousTimeSEM.ctsem_evaluate(::_JointMock, x::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint) =
    (value = -(x[1] - x[2])^2 - 0.01 * sum(abs2, x), gradient = nothing)

@testset "the move out of a degenerate corner is not along a coordinate" begin
    mock = _JointMock()
    estimate = [10.0, 10.0]
    value = ContinuousTimeSEM.ctsem_evaluate(mock, estimate; gradient=false).value
    @test value ≈ -2.0
    # Either coordinate alone is far worse; both together gain 2.
    @test ContinuousTimeSEM.ctsem_evaluate(mock, [0.0, 10.0]; gradient=false).value < -100
    @test ContinuousTimeSEM.ctsem_evaluate(mock, [10.0, 0.0]; gradient=false).value < -100
    @test ContinuousTimeSEM.ctsem_evaluate(mock, [0.0, 0.0]; gradient=false).value ≈ 0.0

    out = ContinuousTimeSEM._ctsem_overshot(mock, estimate, [1], value, 1e-3)
    @test out.overshot
    @test out.coordinates == [1, 2]
    # 1.5, not the 2.0 available at the origin: the probe stops at the first
    # improvement, which here is the half-way pullback. The gain is a lower
    # bound on what is left, and it is reported as one -- the question it
    # answers is whether this is a maximum, and any improvement settles that.
    @test out.gain ≈ 1.5 atol = 1e-8

    # The old probe, on the same point, says maximum. Both coordinates are
    # flagged, so this is not a matter of flagging the right one.
    previous = ContinuousTimeSEM.ctsem_set_overshoot_probe!(:saturation)
    try
        singly = ContinuousTimeSEM._ctsem_overshot(mock, estimate, [1, 2], value, 1e-3)
        @test !singly.overshot
    finally
        ContinuousTimeSEM.ctsem_set_overshoot_probe!(previous)
    end
end

@testset "the pullback sets are magnitude-ordered prefixes" begin
    sets = ContinuousTimeSEM._ctsem_pullback_sets([0.5, -9.0, 2.0, 0.1, -4.0])
    # Ordered by |raw|: 2 (9), 5 (4), 3 (2), 1 (0.5), 4 (0.1).
    @test length.(sets) == [1, 2, 3, 4, 5]
    @test sets[1] == [2]
    @test sets[2] == [2, 5]
    @test sets[3] == [2, 5, 3]
    @test sort(sets[5]) == [1, 2, 3, 4, 5]
    # Every set is a prefix of the next. Every size is present, which is what
    # a geometric ladder gave up and what dataset 12's escaping set of three
    # needed -- see `_ctsem_pullback_sets`.
    for i in 1:(length(sets) - 1)
        @test sets[i] == sets[i + 1][1:length(sets[i])]
    end
    # One coordinate, and none.
    @test ContinuousTimeSEM._ctsem_pullback_sets([3.0]) == [[1]]
    @test isempty(ContinuousTimeSEM._ctsem_pullback_sets(Float64[]))
end


# The convergence criterion. `converge_tol` is in nats and what is tested
# against it is the objective still available, so these are about units and
# about what the verdict refuses to look at, not about one model's numbers.
#
# `_ctsem_optimise_verdict` needs no usable objective here. The probe does now
# evaluate whatever it is handed -- it no longer returns early on an empty
# saturation flag -- but `nothing` has no `ctsem_evaluate` method, so every
# probe point is refused and no gain is claimed. That is the same guard a real
# objective's non-finite point gets, exercised on the cheapest possible
# objective, and it leaves the convergence rule testable on its own.
@testset "convergence is judged on the objective still available" begin
    verdict(gain, last = Inf; g = 1.0e3, tol = 1.0e-6) =
        ContinuousTimeSEM._ctsem_optimise_verdict(nothing, [1.0], [0.0],
            -1483.0, g, gain, last, Int[], 1e-8, tol)

    # Below the tolerance is converged, above it is not -- and the gradient is
    # a thousand in both, which is the point: it is not consulted.
    @test verdict(1.0e-9).converged_enough
    @test !verdict(1.0e-3).converged_enough
    @test verdict(1.0e-9; g = 1.0e9).converged_enough

    # Either condition suffices, and the second is the one that carries a fit
    # whose metric has been corrupted by a flat direction: measured on a
    # saturated fit, `1/2 g'Bg` read 1.002 where the exact gap was 2.1e-15,
    # while its last iteration gained nothing at all.
    @test verdict(1.002, 0.0).converged_enough
    @test !verdict(1.002, 0.093).converged_enough

    # The tolerance is absolute, in nats. A log likelihood difference is a
    # likelihood ratio, so it means the same thing at any N and needs no scale
    # -- which is why nothing here divides by the objective.
    @test verdict(2.0e-6; tol = 1.0e-6).converged_enough == false
    @test verdict(2.0e-6; tol = 1.0e-5).converged_enough

    # A run with neither quantity yet has `Inf` for both, so a fit that took no
    # step cannot pass on an uninitialised number.
    @test !verdict(Inf, Inf).converged_enough

    # A NaN gradient is refused even with nothing left to gain: `Optim` can set
    # its own flag on an earlier iterate, and `NaN <= tolerance` is false, so
    # this is tested rather than inferred. One draw in ten reported convergence
    # this way.
    @test !verdict(1.0e-9; g = NaN).converged_enough
end


