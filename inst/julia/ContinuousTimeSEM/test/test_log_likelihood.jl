using LinearAlgebra

# The Kalman likelihood helper should match the direct Gaussian formula while
# reusing its solve buffer for `S \ y`.
@testset "Kalman innovation log-likelihood" begin
    Σ = [2.0 0.3; 0.3 1.5]
    ỹ = [0.5, -1.0]
    chol = cholesky(Symmetric(copy(Σ)))
    ll_buffer = zeros(2)

    ll = ContinuousTimeSEM._kalman_loglikelihood_cholesky!(ll_buffer, chol, ỹ, log(2π))
    expected = -0.5 * (length(ỹ) * log(2π) + logdet(Σ) + dot(ỹ, Σ \ ỹ))

    @test ll ≈ expected atol = 1e-12
    @test ll_buffer ≈ Σ \ ỹ atol = 1e-12
end
