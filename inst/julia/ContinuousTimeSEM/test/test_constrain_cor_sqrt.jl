using LinearAlgebra

# The loop, vectorized, and buffered correlation transforms should agree. The
# resulting factor should also produce a positive-definite correlation matrix.
@testset "Correlation square-root constraints" begin
    mat = [
        1.2 0.0 0.0
        0.2 0.9 0.0
       -0.1 0.3 1.1
    ]

    loop_factor = ContinuousTimeSEM.constraincorsqrt1(mat)
    vector_factor = ContinuousTimeSEM.constraincorsqrt1_vec(Symmetric(mat, :L))

    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 3)
    ContinuousTimeSEM.constraincorsqrt1_vec!(buffer, mat)
    buffered_factor = copy(buffer.out)

    @test loop_factor ≈ vector_factor atol = 1e-12
    @test buffered_factor ≈ vector_factor atol = 1e-12

    corr = buffered_factor * buffered_factor'
    # The transform intentionally leaves an epsilon margin on the diagonal.
    @test diag(corr) ≈ fill(1.0 + 1e-5, 3) atol = 1e-12
    @test ContinuousTimeSEM.is_positive_definite(Matrix(corr))
end

# The allocating and buffered covariance paths are both used by higher-level
# code, so they should produce the same covariance from the same parameters.
@testset "Covariance conversion" begin
    mat = [
        1.2 0.0 0.0
        0.2 0.9 0.0
       -0.1 0.3 1.1
    ]

    cov = ContinuousTimeSEM.sdcovsqrt2cov(mat, 0)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 3)
    ContinuousTimeSEM.sdcovsqrt2cov!(buffer, mat, 0)

    @test Matrix(cov) ≈ buffer.out atol = 1e-12
    @test ContinuousTimeSEM.is_positive_definite(Matrix(cov))
end
