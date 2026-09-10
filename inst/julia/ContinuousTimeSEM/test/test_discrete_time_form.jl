using ComponentArrays
using LinearAlgebra

# The one-dimensional case has closed-form discrete drift, intercept, and
# process covariance, making it a compact regression test for discretization.
@testset "One-dimensional continuous-to-discrete form" begin
    a = 0.7
    q = 0.2
    c = 0.3
    Δt = 0.4

    pars = ComponentVector(
        DRIFT = [-a;;],
        JAx = [-a;;],
        CINT = [c],
    )
    discrete_ca = ContinuousTimeSEM._make_discrete_ca_buffer(Float64, 1)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 1)
    exp_buffer = ContinuousTimeSEM.ExpBuffer(pars.DRIFT)
    lyap_buffer = ContinuousTimeSEM.LyapBuffer(pars.DRIFT)

    ContinuousTimeSEM._compute_discrete_time_form!(
        discrete_ca,
        buffer,
        [q;;],
        pars,
        Δt,
        exp_buffer,
        lyap_buffer,
    )

    expected_drift = exp(-a * Δt)
    expected_intercept = c * (1 - expected_drift) / a
    expected_diffusion = q / (2a) * (1 - exp(-2a * Δt))

    @test discrete_ca.dDRIFT[1, 1] ≈ expected_drift atol = 1e-12
    @test discrete_ca.eJAx[1, 1] ≈ expected_drift atol = 1e-12
    @test discrete_ca.dINT[1] ≈ expected_intercept atol = 1e-12
    @test discrete_ca.dDIFFUSION[1, 1] ≈ expected_diffusion atol = 1e-12
end

@testset "Local affine correction follows Stan propagation" begin
    Δt = 0.4
    state = [2.0]
    pars = ComponentVector(DRIFT=[-0.5;;], JAx=[-0.7;;], CINT=[1.0])
    discrete_ca = ContinuousTimeSEM._make_discrete_ca_buffer(Float64, 1)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 1)
    exp_buffer = ContinuousTimeSEM.ExpBuffer(pars.DRIFT)
    lyap_buffer = ContinuousTimeSEM.LyapBuffer(pars.DRIFT)
    diffusion_buffer = ContinuousTimeSEM._make_square_buffer(Float64, 1)
    discretization_buffer = ContinuousTimeSEM._make_discretization_buffer(Float64, 1)

    ContinuousTimeSEM._compute_discrete_time_form!(discrete_ca, buffer, [0.2;;],
        pars, Δt, exp_buffer, lyap_buffer, state, [1], diffusion_buffer,
        discretization_buffer, Val(1))

    jacobian_transition = exp(-0.7 * Δt)
    f0 = -0.5 * state[1] + 1.0
    affine_intercept = f0 - (-0.7 * state[1])
    expected = jacobian_transition * state[1] +
        (jacobian_transition - 1) / -0.7 * affine_intercept
    @test discrete_ca.dDRIFT[1, 1] ≈ jacobian_transition atol = 1e-12
    @test discrete_ca.dDRIFT[1, 1] * state[1] + discrete_ca.dINT[1] ≈ expected atol = 1e-12
end

# A random effect augments the state with a static coordinate, and DRIFT must
# be the true one-step map on that row: diagonal 1, because a static state
# carries forward whole. `.ctJuliaAugmentRandomEffects()` guarantees it for a
# discrete-time model, matching what `ctJacobian()` does to the copy of DRIFT
# it builds JAx from. This asserts the invariant from the engine's side --
# given that DRIFT, the offset cancels on the static row and the state
# survives the step.
#
# When the padding left that diagonal at 0 the offset came out as
# `(0 - 1) x[i] = -x[i]`, cancelling the state the transition had just
# carried: `x_next[i] = 1*x[i] - x[i] = 0`. Every random-effect coordinate was
# zeroed at every step, so a subject's deviation reached the first transition
# and nothing after it, with a finite likelihood throughout. On the tutorial
# handbook's two-variable ESM model that read as an autoregression of .993
# instead of .709 -- a near-unit root absorbing the between-person differences
# in level the intercept could no longer carry.
#
# State 2 is the static coordinate, and state 1 reads it the way an
# `intoverpop` CINT does: JAx[1,2] is the multiplier, CINT[1] the value.
@testset "One-step form carries static augmented coordinates forward" begin
    drift = -0.6
    multiplier = 10.0
    state = [2.0, 0.3]
    pars = ComponentVector(
        DRIFT = [drift 0.0; 0.0 1.0],
        JAx = [drift multiplier; 0.0 1.0],
        CINT = [multiplier * state[2], 0.0],
    )
    discrete_ca = ContinuousTimeSEM._make_discrete_ca_buffer(Float64, 2)
    Qc = [0.2 0.0; 0.0 0.0]

    ContinuousTimeSEM._compute_one_step_form!(discrete_ca, Qc, pars, state,
        [1], Val(2))

    x_next = discrete_ca.dDRIFT * state .+ discrete_ca.dINT
    # Stan's discrete-time map: the dynamic block gets DRIFT*x + CINT, the
    # augmented coordinate is left alone (`state[1:nlatent] *= DRIFT'`).
    @test x_next[1] ≈ drift * state[1] + multiplier * state[2] atol = 1e-12
    @test x_next[2] ≈ state[2] atol = 1e-12
    @test discrete_ca.dINT[2] ≈ 0.0 atol = 1e-12
    # Only the diffusing state gets process noise.
    @test discrete_ca.dDIFFUSION[1, 1] ≈ 0.2 atol = 1e-12
    @test discrete_ca.dDIFFUSION[2, 2] ≈ 0.0 atol = 1e-12

    # And the transition the covariance is propagated with must be the
    # Jacobian of exactly that map, or the mean and the covariance describe
    # different systems.
    @test discrete_ca.dDRIFT ≈ pars.JAx atol = 1e-12
end

# The offset must run over EVERY row, not just the diffusing ones. An isolated
# deterministic latent -- no diffusion of its own and no coupling to anything
# that has some -- is excluded from `.ctJuliaDerrind()`'s set, so restricting
# the offset to that set drops its CINT entirely and the parameter stops
# affecting the likelihood at all (gradient exactly zero). Stan applies it:
# `state[CINTnonzero] += CINT[CINTnonzero,1]` is over every row.
#
# State 2 here is that latent: zero process noise, no coupling, its own drift
# and its own intercept, and it is deliberately NOT in the dynamic set passed
# to the engine.
@testset "One-step form keeps CINT on a non-diffusing state" begin
    state = [2.0, 0.3]
    pars = ComponentVector(
        DRIFT = [-0.6 0.0; 0.0 0.8],
        JAx = [-0.6 0.0; 0.0 0.8],
        CINT = [0.4, 0.7],
    )
    discrete_ca = ContinuousTimeSEM._make_discrete_ca_buffer(Float64, 2)
    Qc = [0.2 0.0; 0.0 0.0]

    ContinuousTimeSEM._compute_one_step_form!(discrete_ca, Qc, pars, state,
        [1], Val(2))

    # A linear model has JAx == DRIFT, so the offset is exactly CINT.
    @test discrete_ca.dINT[1] ≈ 0.4 atol = 1e-12
    @test discrete_ca.dINT[2] ≈ 0.7 atol = 1e-12
    x_next = discrete_ca.dDRIFT * state .+ discrete_ca.dINT
    @test x_next[2] ≈ 0.8 * state[2] + 0.7 atol = 1e-12
end
