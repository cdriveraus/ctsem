using DataFrames, ForwardDiff, LinearAlgebra, Random

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

# The end-of-run verdict skips the probe where every coordinate is small and
# nothing is flagged (`_CTSEM_OVERSHOOT_MIN_RAW`). The mock peaks at zero and
# the estimate sits at 1.5, so a probe WOULD find a gain there: what is tested
# is that it is not asked, and that a flag or one large coordinate asks it.
@testset "the end-of-run probe is skipped only where every coordinate is small" begin
    mock = _OvershotMock(0.0)
    verdict(x, flagged) = ContinuousTimeSEM._ctsem_optimise_verdict(mock, x,
        zeros(length(x)), -sum(abs2, x), 1e-9, 1e-9, 0.0, flagged, 1e-8, 1e-6)
    @test !verdict([1.5, -0.5], Int[]).overshot
    @test verdict([1.5, -0.5], [1]).overshot
    @test verdict([2.5, -0.5], Int[]).overshot
    @test ContinuousTimeSEM._ctsem_overshoot_skippable([1.9, -1.9], Int[])
    @test !ContinuousTimeSEM._ctsem_overshoot_skippable([NaN, 0.0], Int[])
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

# One free drift and one free diffusion, with the transforms that matter here:
# `-log1p_exp` saturates in one direction and not the other, and `log1p_exp`
# does the reverse. Built locally rather than borrowed from
# `test_adjoint_gradient_validation.jl`, which has the same shape -- a fixture
# reached across files makes a test depend on the suite's include order, and
# that is the defect `test-julia-convergence.R` was rewritten to remove.
function _backend_saturating_parameters()
    matrices = Symbol[]; rows = Int[]; cols = Int[]
    parnumber = Union{Missing,Int}[]; value = Union{Missing,Float64}[]
    transform = Union{Missing,String}[]
    add!(name, i, j, v, pn, tf) = begin
        push!(matrices, name); push!(rows, i); push!(cols, j)
        push!(parnumber, pn); push!(value, v); push!(transform, tf)
    end
    add!(:DRIFT, 1, 1, missing, 1, "-log1p_exp(param[1])")
    add!(:JAx, 1, 1, missing, 1, "-log1p_exp(param[1])")
    add!(:DIFFUSION, 1, 1, missing, 2, "log1p_exp(param[2])")
    add!(:CINT, 1, 1, 0.0, missing, missing)
    add!(:LAMBDA, 1, 1, 1.0, missing, missing)
    add!(:Jy, 1, 1, 1.0, missing, missing)
    add!(:MANIFESTMEANS, 1, 1, 0.0, missing, missing)
    add!(:MANIFESTVAR, 1, 1, 0.3, missing, missing)
    add!(:T0VAR, 1, 1, 0.5, missing, missing)
    add!(:T0MEANS, 1, 1, 0.0, missing, missing)
    add!(:PARS, 1, 1, 0.0, missing, missing)
    df = DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber,
        value=value, transform=transform)
    ekf_from_data_frame(df)
end

# The relative flatness measure, on a real parameter object rather than a mock,
# because what it has to get right is a *transform* -- and the fixture's drift
# is `-log1p_exp(param[1])`, which is exactly the shape that stalled the fit
# this was written for.
@testset "flatness is measured against a transform's own live value" begin
    sp = _backend_saturating_parameters()
    ratios(v) = ContinuousTimeSEM._ctsem_flat_ratios(sp, v, eachindex(v))

    # At the origin every ratio is 1 by construction: that is the reference.
    at_zero = ratios([0.0, 0.0])
    @test all(isapprox(r, 1.0; atol=1e-12) for r in values(at_zero))

    # Run the drift out to where its transform dies and the ratio collapses,
    # while the diffusion's, which has not moved, stays at one. An absolute
    # floor cannot make that distinction without knowing which transform it is
    # looking at.
    #
    # Negative, not positive. `-log1p_exp(param)` has derivative
    # `-sigmoid(param)`, which goes to zero as `param` goes to *minus*
    # infinity and to one as it grows -- so raw -12 is the dead end and raw +12
    # is the live one, where the ratio is 2 rather than small. Which direction
    # a transform dies in is exactly what a magnitude rule cannot know and this
    # measure does not need to.
    far = ratios([-12.0, 0.0])
    @test far[1] < 1e-4
    @test isapprox(far[2], 1.0; atol=1e-12)
    @test ratios([12.0, 0.0])[1] > 1

    @test ContinuousTimeSEM._ctsem_flat_coordinates(sp, [-12.0, 0.0],
        eachindex([12.0, 0.0]); ratio=1e-3) == [1]
    @test isempty(ContinuousTimeSEM._ctsem_flat_coordinates(sp, [12.0, 0.0],
        eachindex([12.0, 0.0]); ratio=1e-3))
end

# A mock whose maximum is exactly where the flat coordinate sits, so nothing a
# pullback can reach is better. This is `test_state_sampling.jl`'s count model
# in miniature: a fit converging *into* a flat region, which stalled-and-flat
# alone stopped before it arrived.
struct _AtPeakMock end
ContinuousTimeSEM.ctsem_evaluate(::_AtPeakMock, x::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint) =
    (value = -sum(abs2, x .- [-12.0, 0.0]), gradient = nothing)

# The conjunction and its hysteresis. Three conditions, and the third is the
# one an earlier version lacked.
@testset "a fit is stopped only when there is somewhere better to go" begin
    sp = _backend_saturating_parameters()
    trace = ContinuousTimeSEM.CTSEMTrace(:objective, :gradient_norm)
    # The flat tail has to be longer than the window, or the window reaches
    # back into the climb and the progress test correctly declines to fire.
    climb = collect(range(-1000.0, -1.0; length = 121))
    for (i, v) in enumerate(vcat(climb, fill(-1.0, 100)))
        ContinuousTimeSEM._record!(trace, i, v, 1.0)
    end
    @test ContinuousTimeSEM._ctsem_stalled(trace, 80, 1e-2)
    watch() = ContinuousTimeSEM.CTSEMStallWatch(window=80, fraction=1e-2,
        cooldown=30, tighten=0.1, tightenings=2)

    # Stalled, nothing flat: slow rather than stuck. Not stopped, and asked
    # less readily next time.
    w = watch()
    before = w.fraction
    @test !ContinuousTimeSEM._ctsem_stall_verdict!(w, trace, 221, _AtPeakMock(),
        sp, [0.0, 0.0], 1:2, 1e-3)
    @test w.triggers == 1
    @test w.fraction == before * 0.1
    @test w.quiet_until == 251

    # Inside the cooldown it does not even look.
    @test !ContinuousTimeSEM._ctsem_stall_verdict!(w, trace, 240, _AtPeakMock(),
        sp, [0.0, 0.0], 1:2, 1e-3)
    @test w.triggers == 1

    # Tightening is capped, so the bar cannot walk away to nothing.
    for iteration in (251, 300, 400, 500)
        ContinuousTimeSEM._ctsem_stall_verdict!(w, trace, iteration,
            _AtPeakMock(), sp, [0.0, 0.0], 1:2, 1e-3)
    end
    @test w.fraction >= before * 0.01 - 1e-18

    # Stalled AND flat, but the objective is at its supremum there: still not
    # stopped. Without this the count model over the joint density -- whose
    # drift arrives at raw -18.5 on a flat transform and is a genuine maximum
    # -- was stopped before it arrived and reported as a failed fit.
    peak = watch()
    @test !ContinuousTimeSEM._ctsem_stall_verdict!(peak, trace, 221,
        _AtPeakMock(), sp, [-12.0, 0.0], 1:2, 1e-3)
    @test isempty(peak.flat)
    @test peak.quiet_until == 251          # treated as the slow case

    # Stalled, flat, and somewhere better: stopped, with the coordinate named
    # and the point to resume from in hand.
    stuck = watch()
    @test ContinuousTimeSEM._ctsem_stall_verdict!(stuck, trace, 221,
        _OvershotMock(1.0), sp, [-12.0, 0.0], 1:2, 1e-3)
    @test stuck.flat == [1]
    @test length(stuck.point) == 2
    @test stuck.gain > 0
    # The point is a pullback of the flat coordinate, not of everything.
    @test stuck.point[2] == 0.0
    @test abs(stuck.point[1]) < 12.0

    # No parameter object is not half a conjunction.
    idle = watch()
    @test !ContinuousTimeSEM._ctsem_stall_verdict!(idle, trace, 221,
        _OvershotMock(1.0), nothing, [-12.0, 0.0], 1:2, 1e-3)
end

# `ctsem_pullback` is what a caller acts on. The mock's maximum is at the
# origin and the estimate is out at [10, 10], where no single coordinate
# improves and both together do -- the degenerate-corner shape.
@testset "the pullback hands back a point, or says there is none" begin
    mock = _JointMock()
    out = ContinuousTimeSEM.ctsem_pullback(mock, [10.0, 10.0]; tolerance=1e-3)
    @test out.found
    @test out.point == [5.0, 5.0]
    @test out.coordinates == [1, 2]
    @test out.gain > 1

    # At a maximum there is nothing to hand back, and the point it returns is
    # the one it was given rather than a half-formed candidate.
    at_peak = ContinuousTimeSEM.ctsem_pullback(mock, [0.0, 0.0]; tolerance=1e-3)
    @test !at_peak.found
    @test at_peak.point == [0.0, 0.0]
end

# The ladder reaches past zero. A contraction cannot change a sign, and one
# measured escape needed a population correlation to go from -2.358 to +0.327.
# A pinned objective: the primitive a profile point and an escape attempt share.
#
# It has to hold the coordinate *and* keep being the objective it wraps. The
# second half is the one that can break silently: the laplace route overrides
# most of the optimisable protocol -- its own trial predicate, trace columns and
# verbose report -- and a wrapper that let any of that fall back to the generic
# default would still run, still converge, and quietly stop refusing points
# whose inner Newton had not converged.
# A `CTSEMOptimisable`, because that is what `ctsem_pin` requires and the
# requirement is the point: the wrapper forwards the whole optimisable
# protocol, so anything it can wrap has to have one.
struct _PinMock <: ContinuousTimeSEM.CTSEMOptimisable
    seen::Vector{Vector{Float64}}
end
_PinMock() = _PinMock(Vector{Float64}[])
function ContinuousTimeSEM.ctsem_evaluate(m::_PinMock, x::AbstractVector;
        gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    push!(m.seen, collect(Float64, x))
    value = -sum(abs2, x) / 2
    # `row_loglik` and `subject_loglik` because the result assembly reads both
    # -- the filter routes decompose the likelihood by row and by subject, and
    # this mock has neither, so they are empty rather than absent. The full set
    # `ctsem_optimize` reads off a final evaluation is value, gradient,
    # row_loglik and subject_loglik.
    return (value = value, gradient = gradient ? -collect(Float64, x) : nothing,
        row_loglik = Float64[], subject_loglik = Float64[])
end
ContinuousTimeSEM._ctsem_params(::_PinMock) = nothing
# No transforms behind this mock, so nothing can saturate. Defined because the
# generic `_ctsem_saturated_for` would go looking for `EKFParameters`, and
# because the pinned wrapper forwards this rather than rebuilding it -- which
# is what the forwarding is for.
ContinuousTimeSEM._ctsem_saturated_for(::_PinMock, minimizer) = Int[]

@testset "a pinned objective holds its coordinates and stays itself" begin
    inner = _PinMock()
    pinned = ContinuousTimeSEM.ctsem_pin(inner, [2], [5.0])

    # The inner objective is evaluated at the pinned value whatever the
    # optimiser hands in, so the pin cannot be walked past.
    out = ContinuousTimeSEM.ctsem_evaluate(pinned, [1.0, -3.0])
    @test inner.seen[end] == [1.0, 5.0]
    @test out.value ≈ -(1.0^2 + 5.0^2) / 2

    # And the gradient comes back zero there, which is what actually holds it:
    # L-BFGS builds its direction from gradients and secant pairs, and a
    # coordinate contributing zero to both keeps what it started with.
    @test out.gradient[1] ≈ -1.0
    @test out.gradient[2] == 0.0

    # A gradient-free evaluation is not given one.
    @test ContinuousTimeSEM.ctsem_evaluate(pinned, [1.0, -3.0];
        gradient=false).gradient === nothing

    # The probe path expands too, or the pullback would measure the objective
    # at a point the fit can never reach.
    @test ContinuousTimeSEM._ctsem_probe_value(pinned, [1.0, -3.0]) ≈
        ContinuousTimeSEM._ctsem_probe_value(inner, [1.0, 5.0])

    # A pinned coordinate is out of the saturation range: it cannot saturate or
    # overshoot, because it cannot move, and leaving it in would let the
    # pullback report a gain from moving something this objective holds still.
    @test ContinuousTimeSEM._ctsem_saturation_range(pinned, [0.0, 0.0]) == [1]
    @test collect(ContinuousTimeSEM._ctsem_saturation_range(inner, [0.0, 0.0])) == [1, 2]

    # The route's own methods are forwarded rather than reimplemented. The
    # label is the one thing deliberately changed, so a progress line says
    # which stage this is.
    @test occursin("pinned", ContinuousTimeSEM._ctsem_optimise_label(pinned))
    @test ContinuousTimeSEM._ctsem_params(pinned) === ContinuousTimeSEM._ctsem_params(inner)

    # Validated at construction, not at the first evaluation: a misspelled
    # index should cost nothing rather than a whole optimisation that pinned
    # the wrong coordinate.
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_pin(inner, [1, 1], [0.0, 0.0])
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_pin(inner, [1], [0.0, 1.0])
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_pin(inner, [1], [NaN])
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_pin(inner, [0], [0.0])
end

@testset "a pinned optimisation is a profile point" begin
    # The definition: maximise over everything else with one coordinate fixed.
    # The constrained maximum is below the free one and the gap is what the
    # profile reports -- here exactly `p^2 / 2`, because the coordinates are
    # independent in this objective.
    inner = _PinMock()
    for at in (0.5, 1.0, 2.0)
        pinned = ContinuousTimeSEM.ctsem_pin(inner, [2], [at])
        result = ContinuousTimeSEM.ctsem_optimize(pinned, [0.0, 0.0];
            maxiter=200, verbose=false, progress=false)
        best = collect(Float64, result.minimizer)
        # The free coordinate finds its own optimum...
        @test isapprox(best[1], 0.0; atol=1e-6)
        # ...and the drop from the free maximum is the pinned coordinate's own
        # contribution, which is what a profile measures.
        @test isapprox(result.maximum_loglik, -at^2 / 2; atol=1e-6)
    end
end

@testset "the pullback fractions include reflections" begin
    fractions = ContinuousTimeSEM._CTSEM_PULLBACK_FRACTIONS
    @test any(f -> f < 0, fractions)
    @test -1.0 in fractions          # a pure sign flip
    @test 0.0 in fractions
    @test issorted(fractions; rev=true)
end

# The stall rule, on the shape it was written for. A trace rather than a fit,
# because what the rule reads is a sequence of objective values and the whole
# question is which sequences it calls stalled -- a fit would take a quarter of
# an hour to produce one of them.
@testset "a run that has stopped converging is stalled" begin
    trace_of(values) = begin
        t = ContinuousTimeSEM.CTSEMTrace(:objective, :gradient_norm)
        for (i, v) in enumerate(values)
            ContinuousTimeSEM._record!(t, i, v, 1.0)
        end
        t
    end
    stalled(values, window = 80, fraction = 1e-5) =
        ContinuousTimeSEM._ctsem_stalled(trace_of(values), window, fraction)

    # The measured fit: 4749 nats climbed, then eighty iterations that moved
    # 0.0155. Its own numbers, so a change to the rule that stopped catching
    # this fails here rather than in a quarter of an hour of wall clock.
    climb = collect(range(-8470.0, -3721.05; length = 241))
    plateau = collect(range(-3721.04874, -3721.03329; length = 80))
    @test stalled(vcat(climb, plateau))
    # And not while it was still climbing: at iteration 240 the window covers
    # the climb, whose share is 7.9e-3.
    @test !stalled(climb)

    # A run still gaining a real share of its own progress is not stalled,
    # however small the absolute numbers are.
    @test !stalled(collect(range(-1000.0, 0.0; length = 300)))

    # Fewer rows than the window: nothing to say yet.
    @test !stalled(collect(1.0:50.0))

    # No progress at all leaves the share undefined rather than zero, and this
    # is not the rule that catches it -- `stalled` in the verdict is, which is
    # why the guard is here and not a division.
    @test !stalled(fill(-5.0, 200))

    # Off, by the control the R side exposes as `optimcontrol$stallwindow = 0`.
    @test !stalled(vcat(climb, plateau), 0)

    # The window is what it says: the same plateau is not stalled when the rule
    # is asked to look further back than the plateau is long.
    @test !stalled(vcat(climb, plateau), 200)
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

@testset "the diagonal L-BFGS scaling learns each coordinate's curvature" begin
    C = ContinuousTimeSEM
    rng = Random.MersenneTwister(3)
    a = exp10.(range(-3, 3; length=12))      # curvatures six orders apart
    # On a separable quadratic the true diagonal is a fixed point of the
    # rescaled update, whatever the pair.
    for _ in 1:5
        s = randn(rng, 12)
        B = C._ctsem_lbfgs_diagonal!(copy(a), s, a .* s, ones(12))
        @test B ≈ a rtol=1e-10
    end
    # From the metric's shape it stays positive. (It does not approach the
    # truth quickly from random pairs -- a factor of ~1000 off after 200 --
    # which is not what it is for: the pairs L-BFGS feeds it follow the path.)
    B = nothing
    for _ in 1:200
        s = randn(rng, 12)
        B = C._ctsem_lbfgs_diagonal!(B, s, a .* s, ones(12))
        @test all(>(0), B)
    end
    # L-BFGS with it reaches the minimum of the ill-conditioned quadratic in
    # fewer iterations than with one secant ratio: measured 55 against 133 at
    # the default memory, and 67 against no convergence in 5000 at memory 3.
    fg! = function (F, G, x)
        G === nothing || (G .= a .* x)
        F === nothing ? nothing : 0.5 * sum(a .* x .^ 2)
    end
    x0 = fill(1.0, 12)
    plain = C._ctsem_lbfgs(fg!, x0; memory=20, maxiter=2000, g_tol=1e-8)
    diag = C._ctsem_lbfgs(fg!, x0; memory=20, maxiter=2000, g_tol=1e-8, diagonal=true)
    @test plain.g_converged && diag.g_converged
    @test diag.iterations < plain.iterations
    @test maximum(abs, diag.minimizer) < 1e-4
end

@testset "the sgd phase climbs an ill-conditioned problem and hands over on its progress rule" begin
    C = ContinuousTimeSEM
    # Minimised: a quadratic with curvatures five orders apart, offset so the
    # start is far from the optimum in every coordinate.
    a = exp10.(range(-2, 3; length=10))
    c = fill(3.0, 10)
    fg! = function (F, G, x)
        G === nothing || (G .= a .* (x .- c))
        F === nothing ? nothing : 0.5 * sum(a .* (x .- c) .^ 2)
    end
    f0 = 0.5 * sum(a .* c .^ 2)
    seen = Int[]
    r = C._ctsem_sgd(fg!, zeros(10); maxiter=3000,
        callback=st -> (push!(seen, st.iteration); false))
    # It climbed most of the way, returned its best point with that point's
    # value, and the callback saw the start and then each iteration in order.
    @test r.minimum < 1e-3 * f0
    @test r.minimum ≈ 0.5 * sum(a .* (r.minimizer .- c) .^ 2)
    @test seen == [0; 2:r.iterations]
    @test r.f_calls == r.g_calls
    # No coordinate moved more than the cap in one step.
    steps = Float64[]
    prev = zeros(10)
    tracked = function (F, G, x)
        push!(steps, maximum(abs, x .- prev)); prev .= x
        fg!(F, G, x)
    end
    C._ctsem_sgd(tracked, zeros(10); maxiter=200)
    @test maximum(steps) <= 0.5 + 1e-12
    # The progress rule ends the phase far sooner than sgd()'s own absolute one,
    # well short of the optimum, which is left to L-BFGS.
    early = C._ctsem_sgd(fg!, zeros(10); maxiter=3000, progress=1e-2)
    @test early.iterations < r.iterations
    @test early.minimum < 0.1 * f0
    # An objective that refuses every step ends the phase rather than looping.
    refuse = function (F, G, x)
        G === nothing || (G .= a .* (x .- c))
        F === nothing ? nothing : (x == zeros(10) ? f0 : 1e300)
    end
    stuck = C._ctsem_sgd(refuse, zeros(10); maxiter=100)
    @test stuck.linesearch_failed
    @test stuck.minimizer == zeros(10)
end

@testset "beyond the dense prefixes the pullback sets double" begin
    C = ContinuousTimeSEM
    n = 100
    x = collect(Float64, n:-1:1)   # already in |raw| order: coordinate k is rank k
    sets = C._ctsem_pullback_sets(x; dense=4)
    @test length.(sets) == [1, 2, 3, 4, 8, 16, 32, 64, 100]
    for i in 1:(length(sets) - 1)
        @test sets[i] == sets[i + 1][1:length(sets[i])]
    end
    # A flagged coordinate inside the dense prefixes adds nothing; one ranked
    # beyond them is tried as a set of its own, with the others flagged.
    @test C._ctsem_pullback_sets(x, [2]; dense=4) == sets
    withflag = C._ctsem_pullback_sets(x, [70, 3]; dense=4)
    @test withflag[5] == [3, 70]
    @test withflag[[1:4; 6:end]] == sets
    # At the default a 16-parameter model keeps every prefix, and a 715-parameter
    # one is probed in 22 sets rather than 715.
    @test length(C._ctsem_pullback_sets(randn(16))) == 16
    @test length(C._ctsem_pullback_sets(randn(715))) == 16 + 5 + 1
end

@testset "the overshoot probe reports each set and stops on an interrupt" begin
    C = ContinuousTimeSEM
    seen = Tuple{Int,Int}[]
    # `nothing` refuses every point, so every set is tried and nothing is found.
    out = C._ctsem_overshot(nothing, collect(Float64, 40:-1:1), Int[], 0.0, 1e-3;
        progress=(d, t) -> push!(seen, (d, t)))
    @test !out.overshot
    total = length(C._ctsem_pullback_sets(zeros(40)))
    @test seen == [(i, total) for i in 0:(total - 1)]
    # R's request to stop is a file; set it by hand rather than through
    # `ctsem_set_interrupt!`, which also detaches the process from its console.
    flag = tempname(); touch(flag)
    saved = (C._CTSEM_INTERRUPT_FILE[], C._CTSEM_PARENT_PID[])
    C._CTSEM_INTERRUPT_FILE[] = flag; C._CTSEM_PARENT_PID[] = 0
    C._CTSEM_INTERRUPT_NEXT[] = 0.0
    try
        @test_throws InterruptException C._ctsem_overshot(nothing,
            collect(Float64, 40:-1:1), Int[], 0.0, 1e-3)
    finally
        C._CTSEM_INTERRUPT_FILE[], C._CTSEM_PARENT_PID[] = saved
        C._CTSEM_INTERRUPT_NEXT[] = 0.0; C._CTSEM_INTERRUPT_SEEN[] = false
        rm(flag; force=true)
    end
end


# The convergence criterion. `converge_tol` is in nats and what is tested
# against it is the objective still available, so these are about units and
# about what the verdict refuses to look at, not about one model's numbers.
#
# `_ctsem_optimise_verdict` needs no usable objective here. The probe does now
# evaluate whatever it is handed -- it no longer returns early on an empty
# saturation flag, and the estimate is put at raw 3 so that the magnitude gate
# (`_CTSEM_OVERSHOOT_MIN_RAW`) does not skip it -- but `nothing` has no
# `ctsem_evaluate` method, so every probe point is refused and no gain is
# claimed. That is the same guard a real objective's non-finite point gets,
# exercised on the cheapest possible objective, and it leaves the convergence
# rule testable on its own.
@testset "convergence is judged on the objective still available" begin
    verdict(gain, last = Inf; g = 1.0e3, tol = 1.0e-6) =
        ContinuousTimeSEM._ctsem_optimise_verdict(nothing, [3.0], [0.0],
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



# A trial point that throws for a numerical reason is rejected. One that throws
# because the code is wrong must not be: a caught `MethodError` once turned a
# refactor into a NaN quadrature gap and an identity sampler metric, and the
# same catch sat around every optimiser trial point.
struct _ThrowMock <: ContinuousTimeSEM.CTSEMOptimisable
    err::Any
end
ContinuousTimeSEM.ctsem_evaluate(m::_ThrowMock, x::AbstractVector; kwargs...) =
    throw(m.err)

@testset "a trial point rejects numerical failures and rethrows bugs" begin
    trial(err) = ContinuousTimeSEM._ctsem_optimise_trial(_ThrowMock(err),
        [0.0], true, :adjoint, 1e10, nothing)

    numerical = trial(DomainError(-1.0, "log of a negative"))
    @test numerical.evaluated === nothing
    @test numerical.valid == false
    @test trial(ArgumentError("matrix is not positive definite")).valid == false

    @test_throws MethodError trial(MethodError(+, ("a",)))
    @test_throws UndefVarError trial(UndefVarError(:nowhere))
    @test_throws BoundsError trial(BoundsError([1], 2))
    @test_throws InterruptException trial(InterruptException())

    # From a spawned region the error arrives wrapped, and the wrapper is not
    # what decides.
    task = @task throw(MethodError(+, ("a",)))
    schedule(task)
    wrapped = try
        wait(task)
        nothing
    catch err
        err
    end
    @test wrapped isa TaskFailedException
    @test ContinuousTimeSEM._ctsem_must_propagate(wrapped)
    @test !ContinuousTimeSEM._ctsem_must_propagate(DomainError(-1.0))
end


# ------------------------------------------------------------------ the endgame
#
# The finish (`_ctsem_newton_finish`) on objectives whose shape is known, so each
# property is asserted where it can only hold for the reason given. Several of
# these moved here from R's test-backend-optimgap.R on 2026-09-25 with the steps
# they tested -- the flat-direction probe, the ladder along negative curvature,
# the damped step's line search -- which moved into the engine then.

# The maximised objective, its gradient and its Hessian, each given, so a test
# can hand the finish curvature that is wrong on purpose. Counts what the value
# path and the Hessian cost.
struct _EndgameMock <: ContinuousTimeSEM.CTSEMOptimisable
    value::Function
    gradient::Function
    hessian::Function
    calls::Base.RefValue{Int}
    hessians::Base.RefValue{Int}
end
_endgame_mock(f; g = x -> ForwardDiff.gradient(f, x),
    h = x -> ForwardDiff.hessian(f, x)) = _EndgameMock(f, g, h, Ref(0), Ref(0))
function ContinuousTimeSEM.ctsem_evaluate(m::_EndgameMock, x::AbstractVector;
        gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    m.calls[] += 1
    y = collect(Float64, x)
    return (value = m.value(y), gradient = gradient ? m.gradient(y) : nothing,
        row_loglik = Float64[], subject_loglik = Float64[])
end
function ContinuousTimeSEM.ctsem_hessian(m::_EndgameMock, x::AbstractVector;
        chunk::Integer=0, progress=nothing)
    m.hessians[] += 1
    return Matrix{Float64}(m.hessian(collect(Float64, x)))
end
ContinuousTimeSEM._ctsem_params(::_EndgameMock) = nothing
ContinuousTimeSEM._ctsem_saturated_for(::_EndgameMock, minimizer) = Int[]
# A route whose Hessian is cheap and exact, so `ctsem_optimize` finishes on it
# as it does on the marginal route.
ContinuousTimeSEM._ctsem_finish_curvature(::_EndgameMock) = :exact
ContinuousTimeSEM._ctsem_cheap_hessian(::_EndgameMock) = true

# The finish from `x0`, on the optimiser's own trial path.
function _endgame_run(m, x0; kwargs...)
    fg! = ContinuousTimeSEM._ctsem_trial_closure(m, :adjoint)
    x = collect(Float64, x0)
    G = zeros(length(x))
    f = fg!(0.0, G, x)
    return ContinuousTimeSEM._ctsem_newton_finish(m, x, f, G, fg!; kwargs...)
end

@testset "the flat probe reports the best actual gain, and zero when there is none" begin
    probe = ContinuousTimeSEM._ctsem_flat_probe
    # A likelihood rising along the direction: the probe has to find it and say
    # how far it went, because that number is what a norm of the gradient
    # cannot give.
    rising = x -> 3 * x[2]
    found = probe(rising, [0.0, 0.0], 0.0, [0.0, 1.0])
    @test found.gain == 12          # the longest step, 4, times 3
    @test found.length == 4
    @test probe(x -> -3 * x[2], [0.0, 0.0], 0.0, [0.0, 1.0]).gain == 0
    # A direction of length zero is nothing to probe, not an error.
    none = probe(rising, [0.0, 0.0], 0.0, [0.0, 0.0])
    @test none.gain == 0
    @test isempty(none.direction)
    @test none.evaluations == 0
    # A point the route refuses is a step not taken, not a failed
    # certification.
    @test probe(x -> -Inf, [0.0, 0.0], 0.0, [0.0, 1.0]).gain == 0
end

@testset "a flat direction that still gains says what the probe measured" begin
    # AnomAuth S1 from default starts, in miniature: a gain of 5.5e-05 at a
    # quarter of a unit and less further out. The length that gave it, the
    # longest the objective could be evaluated at, and the unit direction are
    # what the verdict's message reads (`.ctBackendFlatGainReason()` in R).
    bump = x -> x[2] <= 0.25 ? 2.2e-4 * x[2] : 5.5e-5 - 1e-5 * (x[2] - 0.25)
    out = ContinuousTimeSEM._ctsem_flat_probe(bump, [0.0, 0.0], 0.0, [0.0, 3.0])
    @test out.gain ≈ 5.5e-5
    @test out.length == 0.25
    @test out.longest == 4
    @test out.direction == [0.0, 1.0]
    # The longest length is the longest the objective could be evaluated at,
    # not the longest asked for: refused beyond 2, the probe looked no further
    # than 1.
    short = ContinuousTimeSEM._ctsem_flat_probe(x -> x[2] > 2 ? -Inf : bump(x),
        [0.0, 0.0], 0.0, [0.0, 1.0])
    @test short.longest == 1
end

@testset "at a saddle the ascent is looked for along the negative curvature" begin
    # f(x, y) = x^2/2 - y^2/2, maximised, has a saddle at the origin: a maximum
    # in y, a minimum in x. The Newton step over the trusted direction (y) has
    # nothing to offer there, which is how a rank-deficient laplace fit sat in
    # R's correction loop for hours; the ascent is along x, and the side the
    # gradient leans to is the one to try first.
    saddle = p -> 0.5 * p[1]^2 - 0.5 * p[2]^2
    split = ContinuousTimeSEM._ctsem_information_split([-1.0 0.0; 0.0 1.0])
    @test count(split.negative) == 1
    at = [0.01, 0.0]
    out = ContinuousTimeSEM._ctsem_saddle_ladder(saddle, at, saddle(at), split,
        [0.01, 0.0])
    @test out.best !== nothing
    @test out.best.value > saddle(at)
    @test out.best.point[1] > at[1]        # uphill, on the gradient's side
    @test out.best.point[2] == 0.0
    # At a maximum there is no negative curvature, so nothing is proposed and
    # the objective is never called.
    called = Ref(0)
    bowl = p -> (called[] += 1; -sum(abs2, p))
    none = ContinuousTimeSEM._ctsem_saddle_ladder(bowl, [0.0, 0.0], 0.0,
        ContinuousTimeSEM._ctsem_information_split([2.0 0.0; 0.0 2.0]),
        [0.0, 0.0])
    @test none.best === nothing
    @test called[] == 0
end

@testset "the finish at a saddle escapes along the negative curvature" begin
    # Maximised: -x^2/2 + y^2/2 - y^4/4. The origin is a saddle (a minimum in y)
    # with a zero gradient, so the Newton step has nothing to take and a
    # gradient-based optimiser nowhere to go; the maxima are at y = +-1, a
    # quarter above it.
    m = _endgame_mock(p -> -0.5 * p[1]^2 + 0.5 * p[2]^2 - 0.25 * p[2]^4)
    out = _endgame_run(m, [0.0, 0.0])
    @test out.escapes == 1
    @test -out.f ≈ 0.25 atol = 1e-10
    @test abs(out.x[2]) ≈ 1.0 atol = 1e-6
    # It ended on a maximum, so no saddle is left, and on the exact Hessian
    # there: one at the start, one after the escape.
    @test !out.saddle
    @test out.full_hessians == 2
    @test out.hessian_at == out.x
    @test all(eigvals(Symmetric(-out.hessian)) .> 0)
    @test "saddle" in out.history.kind

    # A saddle the route will not let it leave -- every point off it refused,
    # as a Laplace point is when a unit's inner solve fails -- ends there, and
    # says the ladder was tried: a verdict the certification can read, rather
    # than a point returned with a direction nobody looked along.
    refused = y -> y[2] == 0 ? m.value(y) : -Inf
    stuck = _endgame_run(m, [0.0, 0.0]; value_at = refused)
    @test stuck.escapes == 0
    @test stuck.saddle
    @test stuck.ladder_tried
    @test stuck.x == [0.0, 0.0]
end

@testset "negative curvature does not starve the step that closes the gap" begin
    # Maximised: -(x - 3)^2/2 + y^2/2 - y^4/4, from beside the saddle in y with
    # three units still to go in x -- a gap of 4.5 in the trusted direction.
    # Flooring the negative curvature made the y part of the step 1e8 times too
    # long; the line search shrank the whole step to tame it, so the x part
    # shrank too, and the chord, which keeps its Hessian, did the same every
    # step. Taken at the size of its curvature, the first step is whole.
    m = _endgame_mock(p -> -0.5 * (p[1] - 3)^2 + 0.5 * p[2]^2 - 0.25 * p[2]^4)
    out = _endgame_run(m, [0.0, 0.05]; curvature = :chord)
    @test out.history.alpha[1] == 1.0
    @test abs(out.x[1] - 3) < 1e-6
    @test abs(abs(out.x[2]) - 1) < 1e-3
    @test -out.f ≈ 0.25 atol = 1e-7
    @test out.escapes == 0
    @test out.steps < 20
end

@testset "a chord finish keeps its Hessian when the steps moved the estimate little" begin
    # Maximised: a quartic bowl, so the curvature changes with the point and a
    # Hessian from elsewhere is not the one here. Standard errors near one.
    c = [0.3, -0.2, 0.5]
    f = p -> -0.5 * sum(abs2, p .- c) - 0.1 * sum(q -> q^4, p .- c)
    # From a thousandth of a standard error away: the chord's steps converge
    # and move the estimate far less than a hundredth of one, so the Hessian
    # taken at the hand-over is the final one, and says where it was taken.
    m = _endgame_mock(f)
    x0 = c .+ 1e-3
    near = _endgame_run(m, x0; curvature = :chord)
    @test near.full_hessians == 1
    @test m.hessians[] == 1
    @test near.hessian_at == x0
    @test 0 < near.distance <= ContinuousTimeSEM._CTSEM_HESSIAN_REUSE_SE
    @test near.gain < 1e-8
    @test maximum(abs, near.x .- c) < 1e-8
    # From a fifth of a standard error, the steps move it too far for the
    # hand-over Hessian to stand, and the final one is exact at the estimate.
    far = _endgame_run(_endgame_mock(f), c .+ 0.2; curvature = :chord)
    @test far.full_hessians >= 2
    @test far.hessian_at == far.x
    @test far.distance == 0
    # Walked on the new one as on the first, not one Hessian a step.
    @test all(==("newton"), far.history.kind)
    # And a chord that uses up its steps is walked again on the Hessian where
    # it stopped, rather than given no more steps or one Hessian a step.
    short = _endgame_run(_endgame_mock(f), c .+ 1.0; curvature = :chord, maxit = 2)
    @test short.gain < 1e-8
    @test maximum(abs, short.x .- c) < 1e-3
    @test all(==("newton"), short.history.kind)
    @test short.full_hessians >= 2
    # The exact finish never keeps one from elsewhere.
    exact = _endgame_run(_endgame_mock(f), x0; curvature = :exact)
    @test exact.hessian_at == exact.x
end

@testset "a direction the certification trusts is stepped at its own curvature" begin
    # Maximised: one stiff direction and one at relative curvature 1e-10 --
    # trusted by the certification's rule (above 1e-12) and below the finish's
    # floor (1e-8) -- with a quarter of a unit to go along it, as a random-effect
    # sd on a ray toward zero had on AnomAuth S1. Floored, each chord step went
    # a hundredth of the way and gained a fiftieth of the gap, and the finish ran
    # out of steps; at its own curvature the first step closes it, and the
    # Hessian is kept, since a quarter of a unit is a four-hundredth of a
    # standard error there.
    f = p -> -0.5e6 * p[1]^2 - 0.5e-4 * (p[2] - 0.25)^2
    trusted = _endgame_run(_endgame_mock(f), [0.0, 0.0]; curvature = :chord)
    @test trusted.steps <= 2
    @test trusted.x[2] ≈ 0.25 atol = 1e-8
    @test trusted.full_hessians == 1
    @test trusted.gain < 1e-8
    # Below the certification's rule the direction is flat, the probe's to
    # measure, and the step along it is still floored: a relative curvature of
    # 1e-14 moves the estimate by nothing a step could see.
    g = p -> -0.5e6 * p[1]^2 - 0.5e-8 * (p[2] - 0.25)^2
    flat = _endgame_run(_endgame_mock(g), [0.0, 0.0]; curvature = :chord)
    @test abs(flat.x[2]) < 1e-6
end

@testset "the chord follows a curvature that decays along its walk" begin
    # Maximised: one stiff direction and an exponential tail, the shape of a
    # random-effect sd on a ray toward zero -- the likelihood rises without
    # bound in raw units and the curvature decays with it, so every Hessian is
    # stale one step later. At the curvature of the hand-over each chord step
    # was shorter than the last and the finish ran out of steps, then formed
    # exact Hessians to finish (seven in the stage on AnomAuth S1); the
    # secant along each slow step keeps the steps the length the ray needs.
    b = 1e-4
    f = p -> -0.5e4 * p[1]^2 - b * exp(2 * p[2])
    m = _endgame_mock(f)
    out = _endgame_run(m, [0.0, 0.0]; curvature = :chord, probe = false)
    @test out.gain < 1e-8
    @test out.steps < 30
    @test out.full_hessians <= 2
    @test m.hessians[] == out.full_hessians
    # Walked as far as the tolerance needs: what is left along the ray is the
    # objective's own remaining gain, b exp(2 p2).
    @test b * exp(2 * out.x[2]) < 1e-7
    # The secant itself: exact along the step for a quadratic, and nothing
    # changed when the step measures no positive curvature.
    H = [4.0 1.0; 1.0 3.0]
    s = [1.0, 0.0]
    @test ContinuousTimeSEM._ctsem_secant_along(H, s, [2.0, 7.0])[1, 1] ≈ 2.0
    @test ContinuousTimeSEM._ctsem_secant_along(H, s, [2.0, 7.0])[2, 2] == 3.0
    @test ContinuousTimeSEM._ctsem_secant_along(H, s, [-1.0, 0.0]) === H
end

@testset "a first step that does not contract hands the point back" begin
    # A quadratic whose Hessian is reported a hundred times too large, as the
    # curvature at a point outside Newton's region misleads: the first step goes
    # a hundredth of the way, the gain predicted there has barely contracted,
    # and a finish asked to check (`handback`) gives the point back after it.
    c = [1.0, -2.0]
    f = p -> -0.5 * sum(abs2, p .- c)
    wrong = p -> -100.0 * Matrix{Float64}(I, 2, 2)
    back = _endgame_run(_endgame_mock(f; h = wrong), [4.0, 3.0]; handback = true)
    @test back.handback
    @test back.steps == 1
    @test back.hessian === nothing
    @test -back.f > f([4.0, 3.0])
    # A first step the objective would not take whole is outside the region by
    # itself, even one that then contracts the gain: on a tenth of the true
    # curvature the Newton step overshoots tenfold, the trust region's third
    # try is taken, and the gain there is a seventh of the first.
    short = _endgame_mock(p -> -100 * (p[1] - 0.1)^2; h = p -> fill(-20.0, 1, 1))
    held = _endgame_run(short, [0.0]; handback = true, probe = false)
    @test held.handback
    @test held.steps == 1
    @test held.history.alpha[1] < 1
    # Not asked, the same finish walks on; with the right curvature the first
    # step closes the gap and there is nothing to hand back.
    @test !_endgame_run(_endgame_mock(f; h = wrong), [4.0, 3.0]).handback
    right = _endgame_run(_endgame_mock(f), [4.0, 3.0]; handback = true)
    @test !right.handback
    @test maximum(abs, right.x .- c) < 1e-8
end

# A quadratic whose Hessian at the first point it is asked for is reported a
# hundred times too large, as the curvature at a point outside Newton's region
# misleads, and is right everywhere after.
function _endgame_misled(c)
    asked = Ref(0)
    f = p -> -0.5 * sum(abs2, p .- c)
    n = length(c)
    h = p -> (asked[] += 1) == 1 ? -100.0 * Matrix{Float64}(I, n, n) :
        -Matrix{Float64}(I, n, n)
    return f, h
end

@testset "a finish that may not hand back walks on a misleading first Hessian" begin
    # The first step on it goes a hundredth of the way and contracts the gain
    # by almost nothing. Where the Hessian is too dear to throw away, that used
    # to be thrown away too; the exact finish refreshes on the slow step
    # instead, and the chord corrects its copy along the step by the secant,
    # and both close the gap.
    c = [1.0, -2.0]
    f, h = _endgame_misled(c)
    exact = _endgame_run(_endgame_mock(f; h = h), [4.0, 3.0])
    @test maximum(abs, exact.x .- c) < 1e-8
    @test exact.hessian !== nothing
    @test exact.hessian_at == exact.x
    @test exact.full_hessians <= 3
    f, h = _endgame_misled(c)
    chord = _endgame_run(_endgame_mock(f; h = h), [4.0, 3.0]; curvature = :chord)
    @test maximum(abs, chord.x .- c) < 1e-8
    @test chord.full_hessians <= 2
    @test !chord.handback
    @test !exact.handback
end

@testset "the trust region holds back the direction the model cannot carry" begin
    # One stiff direction a thousandth from its optimum, and one so weakly
    # curved at the start that its Newton step runs 130 raw units, past an
    # optimum 5 away (a pseudo-Huber valley: curvature falls as the distance
    # grows). Damping the whole step, as a line search along it does, left the
    # stiff direction short by as much as the weak one was cut; the trust
    # region damps by curvature, so its first step closes the stiff direction
    # while the weak one is held to what the model bears.
    f = p -> -0.5e4 * p[1]^2 - 1e-2 * (sqrt(1 + (p[2] - 5)^2) - 1)
    norms = Float64[]
    out = _endgame_run(_endgame_mock(f), [1e-3, 0.0]; curvature = :chord,
        probe = false, callback = st -> (push!(norms, st.g_norm); false))
    @test out.history.alpha[1] < 0.1           # the weak direction held back
    # The stiff one taken whole: its gradient, 10 at the start, is gone after
    # the first step, and what is left is the weak direction's, about 0.01.
    @test norms[1] < 0.02
    # Walked on, it reaches the valley floor.
    @test abs(out.x[2] - 5) < 1e-4
    @test abs(out.x[1]) < 1e-8
    # The multiplier puts a step outside the region on its boundary, and
    # leaves one inside it alone.
    mult = ContinuousTimeSEM._ctsem_trust_multiplier
    lam = [1e4, 7.5e-5]; cc = [10.0, -0.0098]
    mu = mult(lam, cc, 2.0)
    @test mu > 0
    @test sqrt(sum(abs2, cc ./ (lam .+ mu))) ≈ 2.0 rtol = 1e-5
    @test mult(lam, cc, 200.0) == 0
end

@testset "the finish probes the directions its Hessian does not trust" begin
    # A direction with no curvature at the estimate and a live gradient, which
    # the gap cannot see: stepping along it gains, so the certification must
    # hear it. And no probe unless a certification asked for one.
    f = p -> -0.5 * p[1]^2 + 1e-14 * p[2] + 8e-4 * p[2]^3
    out = _endgame_run(_endgame_mock(f), [0.0, 0.0])
    @test out.probe.length == 4
    @test out.probe.gain ≈ f([0.0, 4.0]) - f([0.0, 0.0])
    @test out.probe.direction ≈ [0.0, 1.0]
    @test isempty(_endgame_run(_endgame_mock(f), [0.0, 0.0];
        probe = false).probe.direction)
end

@testset "the finish's line search backtracks as far as the arithmetic allows" begin
    # Curvature wrong by a factor of sixteen, as a trusted curvature spanning
    # nine orders once made a Newton step fifty units long: the undamped step
    # overshoots and the search has to halve four times. Nothing is tuned to
    # that -- halving stops at the objective's own resolution, so a differently
    # conditioned model simply takes a different number.
    peak = 1 / 16
    f = p -> -100 * (p[1] - peak)^2
    m = _endgame_mock(f; h = p -> fill(-200 / 16, 1, 1))
    out = _endgame_run(m, [0.0]; curvature = :chord)
    @test out.history.kind[1] == "newton"
    @test out.history.alpha[1] <= 0.125
    @test -out.f > f([0.0])
    @test out.x[1] ≈ peak atol = 1e-6
end

@testset "an increase the objective cannot represent is not a step" begin
    # Armijo alone does not rule this out: sufficient increase scales with the
    # step, so an increase of 1e-20 satisfies it once alpha is small enough, and
    # taking it spends a step on a point no different from the one it left. The
    # gradient and curvature here promise half a nat; the value gains nothing it
    # can represent.
    crumbs = _endgame_mock(p -> p[1] > 0 ? -2705.0 + 1e-20 : -2705.0;
        g = p -> [1.0], h = p -> fill(-1.0, 1, 1))
    out = _endgame_run(crumbs, [0.0]; curvature = :exact)
    @test out.steps == 0
    @test out.x == [0.0]
    # And it ends on the objective's resolution, not on a rung count fitted to
    # some model: a bounded number of evaluations.
    @test crumbs.calls[] < 400
    # A representable, Armijo-sufficient increase is taken.
    real = _endgame_mock(p -> p[1] > 0 ? -2705.0 + 1e-3 : -2705.0;
        g = p -> [1.0], h = p -> fill(-1.0, 1, 1))
    taken = _endgame_run(real, [0.0]; curvature = :exact)
    @test taken.steps >= 1
    @test -taken.f == -2705.0 + 1e-3
end

# A resumed stage starts where the last one stopped, so its own trace shows no
# progress and the share the watch asks for is undefined. Carrying the fit's
# progress in fixes that, and on a resumed stage the watch may then stop on
# progress alone, once the hysteresis every fit runs under is spent: stalled at
# a share of 1e-2, of 1e-3 thirty iterations later, and of 1e-4 thirty after
# that. It is the one rule for a resume that has stopped gaining: it replaced
# resuming with the predicted-gain stop switched off and the cap quadrupled.
@testset "a resumed stage stops on progress alone once its hysteresis is spent" begin
    C = ContinuousTimeSEM
    function first_stop(values; carried, alone = true)
        watch = C.CTSEMStallWatch(window = 80, fraction = 1e-2, cooldown = 30,
            tighten = 0.1, tightenings = 2, carried = carried, alone = alone)
        t = C.CTSEMTrace(:objective, :gradient_norm)
        for (i, v) in enumerate(values)
            C._record!(t, i, v, 1.0)
            C._ctsem_stall_verdict!(watch, t, i, nothing, nothing, Float64[],
                nothing, 1e-3) && return (i, watch.exhausted)
        end
        return (0, watch.exhausted)
    end
    # AnomAuth's resume: flat at -2900.7749 for its whole budget, after a fit
    # that had gained something before it.
    ridge = fill(-2900.7749, 300)
    @test first_stop(ridge; carried = 5.0) == (141, true)
    # Without the carried progress the share is undefined and it never stops,
    # which is the 57 minutes.
    @test first_stop(ridge; carried = 0.0) == (0, false)
    # A first stage keeps the conjunction: no transform layer, no stop.
    @test first_stop(ridge; carried = 5.0, alone = false) == (0, false)
    # A resume that is still gaining is left alone: 108 iterations for 2.9
    # nats, the case switching the predicted-gain rule off was measured on.
    @test first_stop(collect(range(-100.0, -97.1; length = 108));
        carried = 50.0) == (0, false)
end

@testset "ctsem_optimize hands the endgame's numbers back" begin
    # From beside the saddle above: L-BFGS converges into it (the gradient in y
    # is exactly zero along y = 0) and hands over, and the finish escapes.
    m = _endgame_mock(p -> -0.5 * p[1]^2 + 0.5 * p[2]^2 - 0.25 * p[2]^4)
    r = ContinuousTimeSEM.ctsem_optimize(m, [0.5, 0.0]; maxiter=200,
        tune_chunks=false, gap_tol=1e-8, newton=true, certify=true,
        progress=false, verbose=false, overshoot_probe=:off)
    @test r.newton_escapes == 1
    @test abs(r.minimizer[2]) ≈ 1 atol = 1e-6
    @test r.maximum_loglik ≈ 0.25 atol = 1e-10
    @test !r.newton_saddle
    @test size(r.hessian) == (2, 2)
    @test r.hessian_evaluated_at == r.minimizer
    @test r.hessian_distance == 0
    @test r.stop_reason in ("gap", "gradient", "linesearch")
    @test length(r.newton_history_kind) == r.newton_steps
    @test "saddle" in r.newton_history_kind
    # With nothing to certify, no probe; the steps are still taken where the
    # Hessian is cheap, as they always were.
    bare = ContinuousTimeSEM.ctsem_optimize(m, [0.5, 0.0]; maxiter=200,
        tune_chunks=false, gap_tol=1e-8, newton=true, certify=false,
        progress=false, verbose=false, overshoot_probe=:off)
    @test bare.newton_steps >= 1
    @test !bare.probe_ran
end

@testset "ctsem_optimize goes back to L-BFGS when the hand-over was early" begin
    # The quadratic with the hundredfold Hessian, and the switch set so high
    # that L-BFGS hands over after its first iteration: the finish hands back,
    # L-BFGS runs to its own stopping rule, and the finish runs again there.
    # Everything either finish did is counted.
    c = [1.0, -2.0]
    f = p -> -0.5 * sum(abs2, p .- c)
    m = _endgame_mock(f; h = p -> -100.0 * Matrix{Float64}(I, 2, 2))
    r = ContinuousTimeSEM.ctsem_optimize(m, [4.0, 3.0]; maxiter=200,
        tune_chunks=false, gap_tol=1e-8, newton=true, newton_switch=1e6,
        certify=false, progress=false, verbose=false, overshoot_probe=:off)
    @test r.newton_handback
    @test maximum(abs, r.minimizer .- c) < 1e-6
    @test r.newton_steps >= 1
    @test length(r.newton_history_kind) == r.newton_steps
    @test r.newton_hessians == m.hessians[] == 2
    @test r.iterations > r.newton_steps
    @test r.stop_reason in ("gap", "gradient", "linesearch")
    # With the Hessian right, the early hand-over is borne out and kept.
    kept = ContinuousTimeSEM.ctsem_optimize(_endgame_mock(f), [4.0, 3.0];
        maxiter=200, tune_chunks=false, gap_tol=1e-8, newton=true,
        newton_switch=1e6, certify=false, progress=false, verbose=false,
        overshoot_probe=:off)
    @test !kept.newton_handback
    @test maximum(abs, kept.minimizer .- c) < 1e-8
end


@testset "the cap the closing line reports is L-BFGS's, not the finish's steps" begin
    # The finish's steps count in the line but not against `maxiter`: a dear
    # Hessian is walked for as many steps as it cost, so on a large model they
    # alone can pass the cap in a fit that converged. The exponential tail
    # needs about seven steps after a hand-over at L-BFGS's first iteration,
    # on a chord too dear to hand back (a budget of one step).
    m = _endgame_mock(p -> -0.5e4 * p[1]^2 - 1e-4 * exp(2 * p[2]))
    lines = String[]
    r = ContinuousTimeSEM.ctsem_optimize(m, [0.01, 0.0]; maxiter=3,
        tune_chunks=false, gap_tol=1e-8, newton=true, newton_switch=1e6,
        newton_maxit=1, newton_curvature=:chord,
        certify=false, progress=true, progress_every=1e-9,
        progress_sink=(text, kind) -> (push!(lines, text); nothing),
        verbose=false, overshoot_probe=:off)
    @test r.iterations > 3
    @test r.stop_reason != "cap"
    @test any(occursin("iterations", l) for l in lines)
    @test !any(occursin("ITERATION CAP", l) for l in lines)
end

@testset "ctsem_optimize finishes on a Hessian too dear to hand back" begin
    # The quadratic whose first Hessian is a hundredfold wrong, handed over
    # after L-BFGS's first iteration, with a step budget below what a Hessian
    # costs (4 parameters, `newton_maxit = 3`): the first step fails the check
    # a hand-back would make, but the finish may not throw the Hessian away,
    # so it keeps what it formed, refreshes on the slow step and closes the
    # gap itself, and everything it did is counted.
    c = [1.0, -2.0, 0.5, 3.0]
    f, h = _endgame_misled(c)
    m = _endgame_mock(f; h = h)
    r = ContinuousTimeSEM.ctsem_optimize(m, [4.0, 3.0, -1.0, 0.0]; maxiter=200,
        tune_chunks=false, gap_tol=1e-8, newton=true, newton_switch=1e6,
        newton_maxit=3, certify=false, progress=false, verbose=false,
        overshoot_probe=:off)
    @test maximum(abs, r.minimizer .- c) < 1e-8
    @test r.newton_steps >= 1
    @test length(r.newton_history_kind) == r.newton_steps
    @test r.newton_hessians == m.hessians[] <= 3
    @test r.iterations - r.newton_steps <= 2    # L-BFGS did not have to go on
    @test r.stop_reason in ("gap", "gradient", "linesearch")
    @test !r.newton_handback
end

# Shared with `test_laplace.jl`; see `laplace_fixtures.jl`.
isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

# The nonlinear Laplace fixture's model -- a random initial level, a drift that
# grows with the state, free diffusion and manifest mean -- over data simulated
# from a process of that shape, fifty subjects of eight waves, so its optimum
# is interior and determined in every direction. The shared linear fixture is
# not, twice over: it frees T0MEANS, CINT and MANIFESTMEANS together, three
# means of which the data determine two, and its deterministic sinusoids put
# the optimum out along flat transforms (raw drift near 35, diffusion near -47).
# There the floored Newton step walks a nearly flat ray far in raw units however
# little it moves in standard errors -- measured on that fixture and on OU data
# fitted by the same model: a thousandth of a standard error came back as a
# second Hessian. Built on the shared nonlinear fixture's parameter object, so
# no new model type is compiled.
function _endgame_laplace_objective(; nsubjects=50, nobs=8, seed=20260925)
    rng = Random.MersenneTwister(seed)
    starts = Int[]; times = Float64[]; ys = Float64[]
    a = -0.5; dt = 0.5; substeps = 20; h = dt / substeps
    position = 1
    for s in 1:nsubjects
        push!(starts, position)
        state = 0.5 + 0.8 * randn(rng) + 0.4 * randn(rng)
        for t in 1:nobs
            if t > 1
                for _ in 1:substeps
                    state += a * (1 + 0.15 * state) * state * h +
                        0.3 * sqrt(h) * randn(rng)
                end
            end
            push!(times, dt * (t - 1))
            push!(ys, state + 0.2 + 0.3 * randn(rng))
            position += 1
        end
    end
    objective = ctsem_objective(_LAPLACE_NONLINEAR_OBJECTIVE.params, starts,
        times, reshape(ys, 1, :))
    return ctsem_laplace_objective(objective, [1], [5], Int[], [1.0])
end

@testset "a Laplace fit's finish ends on one Hessian when its steps are small" begin
    # The Laplace Hessian is 2 npar warm-started gradients, so the finish takes
    # one at the hand-over, reuses it for its steps (the chord) and keeps it as
    # the final one when they moved the estimate less than a hundredth of a
    # standard error. Deterministic below: from the optimum, and from the
    # optimum moved along its best-determined direction by a thousandth of a
    # standard error and by half of one.
    laplace = _endgame_laplace_objective()
    start = [0.5, -0.4, -1.0, 0.2, 0.0]
    fit = ctsem_laplace_optimize(laplace, start; maxiter=500,
        tune_chunks=false, gap_tol=1e-8, newton=true, certify=true,
        progress=false)
    @test fit.converged
    @test fit.newton_hessians >= 1
    @test size(fit.hessian) == (length(start), length(start))
    @test length(fit.hessian_evaluated_at) == length(start)
    @test fit.hessian_distance <= ContinuousTimeSEM._CTSEM_HESSIAN_REUSE_SE
    best = collect(Float64, fit.minimizer)
    split = ContinuousTimeSEM._ctsem_information_split(-fit.hessian)
    # Interior and determined in every direction, or this tests something else.
    @test all(split.trusted)
    k = argmax(split.values)
    along(amount) = best .+ (amount / sqrt(split.values[k])) .* split.vectors[:, k]
    fg! = ContinuousTimeSEM._ctsem_trial_closure(laplace, :adjoint)
    function finish_from(x0)
        G = zeros(length(x0))
        f = fg!(0.0, G, x0)
        ContinuousTimeSEM._ctsem_newton_finish(laplace, x0, f, G, fg!;
            curvature = :chord)
    end
    here = finish_from(best)
    @test here.full_hessians == 1
    @test here.hessian_at == best
    @test here.distance <= ContinuousTimeSEM._CTSEM_HESSIAN_REUSE_SE
    near = finish_from(along(1e-3))
    @test near.steps >= 1
    @test near.full_hessians == 1
    @test near.hessian_at == along(1e-3)
    @test 0 < near.distance <= ContinuousTimeSEM._CTSEM_HESSIAN_REUSE_SE
    @test near.gain < 1e-8
    far = finish_from(along(0.5))
    @test far.full_hessians >= 2
    @test far.hessian_at == far.x
    @test far.distance == 0
    # And without a certification to read it, the Laplace route runs no finish
    # at all: its Hessian would be one nobody asked for.
    bare = ctsem_laplace_optimize(_endgame_laplace_objective(), start;
        maxiter=500, tune_chunks=false, gap_tol=1e-8, newton=true,
        certify=false, progress=false)
    @test bare.newton_hessians == 0
end

@testset "the marginal route takes the chord where a Hessian outprices its steps" begin
    # A forward-over-reverse Hessian costs about a gradient per parameter, and
    # a refresh can save at most the finish's step budget: above that many
    # parameters the marginal route walks one Hessian (the chord), as the
    # Laplace route always has.
    laplace = _endgame_laplace_objective(; nsubjects=5)
    marginal = laplace.objective
    @test ContinuousTimeSEM._ctsem_finish_curvature(marginal, 30, 30) === :exact
    @test ContinuousTimeSEM._ctsem_finish_curvature(marginal, 31, 30) === :chord
    @test ContinuousTimeSEM._ctsem_finish_curvature(laplace, 5, 30) === :chord
    @test ContinuousTimeSEM._ctsem_finish_curvature(
        _endgame_mock(p -> -sum(abs2, p)), 100, 30) === :exact
end
