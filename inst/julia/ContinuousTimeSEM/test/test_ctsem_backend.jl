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


