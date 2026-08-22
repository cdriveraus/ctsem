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
