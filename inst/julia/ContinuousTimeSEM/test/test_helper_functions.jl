using ForwardDiff
using LinearAlgebra

# `add_diag!` is used in numerical kernels, so test both scalar and vector
# updates without relying on broader matrix routines.
@testset "Diagonal updates" begin
    mat = [1.0 2.0; 3.0 4.0]
    ContinuousTimeSEM.add_diag!(mat, 10.0)
    @test mat == [11.0 2.0; 3.0 14.0]

    ContinuousTimeSEM.add_diag!(mat, [1.0, -2.0])
    @test mat == [12.0 2.0; 3.0 12.0]
end

# These helpers replace `mul!` in paths involving views or generic scalar
# types. Compare them directly with the mathematical products they implement.
@testset "Matrix-vector and transpose products" begin
    A = [1.0 2.0 3.0; 4.0 5.0 6.0]
    x = [0.5, -1.0, 2.0]
    y = zeros(2)

    ContinuousTimeSEM._matvec_mul!(y, A, x)
    @test y ≈ A * x

    ContinuousTimeSEM._matvec_mul!(y, A, reshape(x, :, 1))
    @test y ≈ A * x

    B = [2.0 -1.0 0.5; 0.0 1.5 -2.0]
    C = zeros(2, 2)
    ContinuousTimeSEM._mul_right_transpose!(C, A, B)
    @test C ≈ A * B'

    Abig = BigFloat.(A)
    Bbig = BigFloat.(B)
    Cbig = zeros(BigFloat, 2, 2)
    ContinuousTimeSEM._mul_right_transpose!(Cbig, Abig, Bbig)
    @test Cbig == Abig * Bbig'
end

# Approximation helpers decide whether matrix exponentials can be reused. Check
# primal-value handling for dual numbers and the intended tolerance behavior.
@testset "Approximation helpers" begin
    A = [1.0 0.0; 0.0 1.0]
    B = [1.0 + 1e-12 0.0; 0.0 1.0 - 1e-12]
    C = [1.0 + 1e-3 0.0; 0.0 1.0 - 1e-3]

    @test ContinuousTimeSEM._custom_abs(-2.0) == 2.0
    @test ContinuousTimeSEM._custom_abs(ForwardDiff.Dual(-2.0, 1.0)) == 2.0
    @test ContinuousTimeSEM._isapprox_default_rtol(Float64) == sqrt(eps(Float64))
    @test ContinuousTimeSEM._isapprox_matrix_noalloc(A, B)
    @test !ContinuousTimeSEM._isapprox_matrix_noalloc(A, C)
    @test ContinuousTimeSEM._can_reuse_same_exponential(A, B)
    @test !ContinuousTimeSEM._can_reuse_same_exponential(A, C)
end
