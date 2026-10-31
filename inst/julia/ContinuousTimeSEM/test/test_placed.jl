using ForwardDiff, LinearAlgebra, Statistics

# The placed objective (`laplace_placed.jl`): the hybrid with each flagged
# unit's rule -- product Gauss-Hermite or lattice-t (`point_sets.jl`) -- placed
# at its mode and curvature at every evaluation, and its exact gradient
# through the placement.
#
# What has to hold, each checked on its own:
#
#  1. At the continuation's centre the placed value is the continuation's
#     placed value, unit by unit: the same rule at the same point.
#  2. The gradient is the derivative of the value: central differences of the
#     placed value, at the centre and away from it, on a nonlinear one-level
#     fixture, a nonlinear two-level one (a leaf's conditional mode chained
#     through the root's nodes), a linear two-level one, and on lattice
#     rules: four effects a subject, and a four-effect study over subjects.
#  3. On a Gaussian integrand the rule is exact wherever it sits, so the
#     placement terms vanish and the gradient is the frozen rule's.
#  4. The Hessian, localised to each column's readers, is the plain central
#     difference of the exact gradient.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

# Subjects nested in two studies on the nonlinear objective: a T0MEANS effect
# a subject, whose state-dependent drift makes the integrand non-Gaussian, and
# a MANIFESTMEANS effect a study. Raw layout: 1-4 the model parameters, 5 the
# subject scale, 6 the study scale.
_fresh_placed_nested() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 4], sd_index=[5, 6], cor_index=Int[], sd_scale=[1.0, 1.0],
    level_nre=[1, 1], group=vcat(1:5, [1, 1, 2, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.1])

# Four effects a subject on the nonlinear objective: T0MEANS, the
# state-dependent drift, DIFFUSION and MANIFESTMEANS. Raw layout: 1-4 the model
# parameters, 5-8 the scales, 9-14 the correlations.
_fresh_placed_wide() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 2, 3, 4], sd_index=[5, 6, 7, 8], cor_index=collect(9:14),
    sd_scale=[1.0, 1.0, 1.0, 1.0]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.6, -0.6, -0.4, 0.1, 0.0, -0.1, 0.05, 0.0, 0.1])

# Every subject's T0MEANS under a study's four effects: each study a root of
# four effects, on the lattice, over leaves of one, on product rules. Raw
# layout: 1-4 the model parameters, 5 the subject scale, 6-9 the study scales,
# 10-15 their correlations.
_fresh_placed_wide_nested() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 1, 2, 3, 4], sd_index=[5, 6, 7, 8, 9], cor_index=collect(10:15),
    sd_scale=[1.0, 1.0, 1.0, 1.0, 1.0], level_nre=[1, 4],
    group=vcat(1:5, [1, 1, 2, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.5, -0.6, -0.6, -0.4, 0.1, 0.0, -0.1, 0.05,
     0.0, 0.1])

# Every unit flagged (a negative tolerance), each on its screened rule
# (`draws_se` out of reach, so the refine keeps it).
_placed_from(laplace, theta) = ctsem_laplace_placed(ctsem_laplace_continuation(
    laplace, theta; tolerance=-1.0, draws_se=1e9))

@testset "the placed value is the continuation's placed value at its centre" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("nested", _fresh_placed_nested),
                           ("two levels, linear", _fresh_twolevel),
                           ("lattice", _fresh_placed_wide),
                           ("lattice root", _fresh_placed_wide_nested))
        laplace, values = fresh()
        theta = collect(Float64, values)
        o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0, draws_se=1e9)
        P = ctsem_laplace_placed(o)
        @test (label, length(P.placed)) == (label, length(o.flagged))
        r = ctsem_laplace_placed_evaluate(P, theta; gradient=false)
        @test (label, isapprox(r.value, o.centre_value; rtol=1e-10)) == (label, true)
        @test (label, maximum(abs, r.unit_loglik[o.flagged] .-
            o.centre_quadrature[o.flagged]) < 1e-9) == (label, true)
    end
end

@testset "the placed gradient is the derivative of the placed value" begin
    for (label, fresh) in (("one level", _fresh_nonlinear),
                           ("nested", _fresh_placed_nested),
                           ("two levels, linear", _fresh_twolevel),
                           ("lattice", _fresh_placed_wide),
                           ("lattice root", _fresh_placed_wide_nested))
        laplace, values = fresh()
        theta = collect(Float64, values)
        P = _placed_from(laplace, theta)
        lattice = label in ("lattice", "lattice root")
        @test (label, lattice == any(U -> any(ps -> ps.family === :lattice,
            P.sets[U]), P.placed)) == (label, true)
        n = length(theta)
        for shift in (zeros(n), fill(0.05, n), [0.1 * (-1)^j for j in 1:n])
            r = ctsem_laplace_placed_check(P, theta .+ shift; step=1e-6)
            @test (label, all(isfinite, r.gradient)) == (label, true)
            @test (label, all(isapprox.(r.gradient, r.differences; rtol=1e-5,
                atol=1e-6))) == (label, true)
        end
    end
end

@testset "on a Gaussian integrand the placement terms vanish" begin
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    o = ctsem_laplace_continuation(laplace, theta; tolerance=-1.0)
    P = ctsem_laplace_placed(o)
    placed = ctsem_laplace_placed_evaluate(P, theta; gradient=true).gradient
    frozen = ctsem_laplace_continuation_evaluate(o, theta; gradient=true).gradient
    @test all(isapprox.(placed, frozen; rtol=1e-7, atol=1e-8))
end

@testset "a lattice rule past its target is raised, and the error is read where it sits" begin
    C = ContinuousTimeSEM
    laplace, values = _fresh_placed_wide()
    theta = collect(Float64, values)
    P = _placed_from(laplace, theta)
    ctsem_laplace_placed_evaluate(P, theta; gradient=false)
    U = P.placed[1]
    @test P.sets[U][1].family === :lattice && P.sets[U][1].count == C._LATTICE_START
    before = P.se[U]
    @test isfinite(before) && before > 0
    # Under a target the start's spread exceeds, the unit is rough and its
    # points double until the spread is within it, or the lattice's most.
    P.cont.draws_se = before / 3
    @test U in C._placed_rough(P, 1.0)
    @test C._placed_raise!(P, theta, [U]) == 1
    @test P.sets[U][1].count > C._LATTICE_START
    @test P.cont.unit_rules[U].points == P.sets[U][1].count
    @test P.se[U] < before
    @test P.se[U] <= P.cont.draws_se || P.sets[U][1].count == C._LATTICE_MAX
    # The error at the last evaluation: the lattice units' spread there in
    # quadrature.
    ctsem_laplace_placed_evaluate(P, theta .+ 0.05; gradient=false)
    e = ctsem_laplace_placed_error(P)
    lat = [V for V in P.placed if P.sets[V][1].family === :lattice]
    @test isapprox(e.error_draws, sqrt(sum(abs2, P.se[lat])); rtol=1e-12)
    @test e.centre_error[lat] == P.se[lat]
    @test isapprox(e.error, e.error_draws + e.error_nodes + e.error_residual; rtol=1e-12)
    # The slope check: the lattice units' gradients on two independent sets
    # of shifts, through the curvature, in standard errors. Every unit a
    # bar of zero leaves carrying something, none under an infinite one, and
    # the other sets are lent and taken back.
    H = ctsem_laplace_placed_hessian(P, theta)
    sets = copy(P.sets)
    loose = C._placed_slope_check(P, theta, H; bar=Inf)
    @test isfinite(loose.sd) && loose.sd > 0
    @test isempty(loose.units)
    tight = C._placed_slope_check(P, theta, H; bar=0.0)
    @test tight.sd == loose.sd
    @test !isempty(tight.units) && issubset(tight.units, lat)
    @test loose.doublings == 0 && 1 <= tight.doublings <= 5
    @test all(V -> P.sets[V] === sets[V], P.placed)
    # Without a curvature there is nothing to carry the error through.
    @test isempty(C._placed_slope_check(P, theta, zeros(1, 1)).units)
    # The choice by the estimate's error: under no bar each unit keeps the
    # cheapest rule it is measured on, the first lattice; under any bar every
    # unit's rule is one of its families' and its sets are that rule's.
    Q = _placed_from(laplace, theta)
    @test C._placed_choose!(Q, theta, H; bar=Inf) == 0
    @test all(V -> Q.cont.unit_rules[V] == C.CTSEMRuleChoice(:lattice, 5, C._LATTICE_START),
        Q.placed)
    C._placed_choose!(Q, theta, H; bar=1e-3)
    for V in Q.placed
        c = Q.cont.unit_rules[V]
        @test c.family in (:product, :lattice)
        @test Q.sets[V][1].family === c.family
        @test Q.sets[V][1].count == (c.family === :product ? c.nodes : c.points)
    end
    @test C._placed_choose!(Q, theta, zeros(1, 1)) == 0
end

@testset "the localised Hessian is the difference of the exact gradient" begin
    laplace, values = _fresh_placed_nested()
    theta = collect(Float64, values)
    P = _placed_from(laplace, theta)
    H = ctsem_laplace_placed_hessian(P, theta; step=1e-4)
    @test all(isfinite, H)
    @test H == transpose(H)
    home = deepcopy(P.anchors)
    n = length(theta)
    D = zeros(n, n)
    for t in 1:n
        h = 1e-4 * max(1.0, abs(theta[t]))
        xp = copy(theta); xp[t] += h
        xm = copy(theta); xm[t] -= h
        P.anchors = deepcopy(home); P.cache = nothing; P.value_cache = nothing
        gp = ctsem_laplace_placed_evaluate(P, xp; gradient=true).gradient
        P.anchors = deepcopy(home); P.cache = nothing; P.value_cache = nothing
        gm = ctsem_laplace_placed_evaluate(P, xm; gradient=true).gradient
        D[:, t] .= (gp .- gm) ./ (2h)
    end
    D = (D .+ transpose(D)) ./ 2
    @test all(isapprox.(H, D; rtol=1e-6, atol=1e-8))
end

@testset "within an optimisation the anchors move only at accepted points" begin
    C = ContinuousTimeSEM
    laplace, values = _fresh_nonlinear()
    theta = collect(Float64, values)
    P = _placed_from(laplace, theta)
    C._ctsem_anchor_start!(P, theta)
    @test P.running
    v0 = ctsem_laplace_placed_evaluate(P, theta; gradient=false).value
    home = deepcopy(P.anchors)
    far = theta .+ 0.5
    # A trial moves no anchor, so the value at the accepted point does not
    # depend on what was tried after it.
    ctsem_laplace_placed_evaluate(P, far; gradient=false)
    @test P.anchors == home
    P.cache = nothing; P.value_cache = nothing
    @test ctsem_laplace_placed_evaluate(P, theta; gradient=false).value ≈ v0 atol=1e-12
    # An accepted step does.
    ctsem_laplace_placed_evaluate(P, far; gradient=false)
    C._ctsem_anchor_here!(P)
    @test P.anchors != home
    C.ctsem_laplace_clear_anchors!(P)
    @test !P.running
    # Outside a run every evaluation moves them.
    ctsem_laplace_placed_evaluate(P, theta; gradient=false)
    @test maximum(maximum(abs, P.anchors[U] .- home[U]) for U in P.placed) < 1e-8
end

# A loading whose sign the data cannot fix (test-julia-laplace-modes.R's
# fixture): y = lambda eta with eta centred, so each subject's likelihood
# depends on lambda^2 and its integrand in the loading's effect has a maximum
# either side of minus the population mean over its sd. Raw layout: 1 the
# loading's mean, 2 its population scale.
function _build_placed_sign_objective(; nsubjects=6)
    df = _laplace_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.8;;],
        lambda=[0.0;;], jy=[0.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[1.0;;], t0means=[0.0;;],
        free=Dict((:LAMBDA, 1, 1) => (1, "param[1]"), (:Jy, 1, 1) => (1, "param[1]")))
    sp = ekf_from_data_frame(df)
    subject_starts, times, data = _laplace_test_data(nsubjects, 1; seed=7)
    return ctsem_objective(sp, subject_starts, times, 4 .* data)
end
const _PLACED_SIGN_OBJECTIVE = _build_placed_sign_objective()
_fresh_placed_sign() = ctsem_laplace_objective(_PLACED_SIGN_OBJECTIVE, [1], [2], Int[],
    [1.0])

@testset "a unit with two maxima averages its rules placed at each" begin
    C = ContinuousTimeSEM
    # Mirror images at a population mean of zero: the two rules are one rule
    # reflected, and the average is either.
    theta = [0.0, 0.3]
    P = _placed_from(_fresh_placed_sign(), theta)
    @test all(U -> length(P.alt_anchors[U]) == 1, P.placed)
    both = ctsem_laplace_placed_evaluate(P, theta; gradient=false)
    @test all(U -> P.live[U] == 2, P.placed)
    previous = C.ctsem_set_mode_sum!(false)
    one = try
        ctsem_laplace_placed_evaluate(_placed_from(_fresh_placed_sign(), theta), theta;
            gradient=false)
    finally
        C.ctsem_set_mode_sum!(previous)
    end
    @test both.unit_loglik[P.placed] ≈ one.unit_loglik[P.placed] atol=1e-8
    # Away from it the two rules differ and the weights move with the
    # parameters; the gradient is still the value's. The far maximum is
    # followed from there, as an optimisation follows it: a search started
    # out here does not reach it.
    P = _placed_from(_fresh_placed_sign(), [0.0, 0.3])
    for (i, m) in enumerate(range(0.1, 1.2; length=12))
        ctsem_laplace_placed_evaluate(P, [m, 0.3]; gradient=false)
        i in (4, 12) || continue
        @test all(U -> P.live[U] == 2, P.placed)
        r = ctsem_laplace_placed_check(P, [m, 0.3]; step=1e-6)
        @test all(isfinite, r.gradient)
        @test r.gradient ≈ r.differences rtol=1e-5 atol=1e-6
    end
    # Along the mean the maximum on the far side falls away, its weight gone
    # before it is dropped (here past `_LAPLACE_MODE_SPAN` below the near one):
    # the term has no step where it goes, and the gradient is the value's on
    # both sides. Followed from one point to the next, as an optimisation's
    # accepted points follow them.
    P = _placed_from(_fresh_placed_sign(), [0.0, 0.3])
    grid = collect(range(0.0, 3.0; length=301))
    h = step(range(0.0, 3.0; length=301))
    value = Float64[]
    live = Int[]
    for m in grid
        push!(value, ctsem_laplace_placed_evaluate(P, [m, 0.3]; gradient=false).value)
        push!(live, sum(P.live[P.placed]))
    end
    @test live[1] == 2 * length(P.placed)
    @test live[end] < live[1]
    d2 = abs.(diff(diff(value)))
    @test maximum(d2) < 50 * h^2 * (1 + maximum(abs, diff(value)) / h)
    for m in (grid[findfirst(<(live[1]), live)] - h, grid[findfirst(<(live[1]), live)])
        Q = _placed_from(_fresh_placed_sign(), [m, 0.3])
        r = ctsem_laplace_placed_check(Q, [m, 0.3]; step=1e-7)
        @test r.gradient ≈ r.differences rtol=1e-4 atol=1e-5
    end
end
