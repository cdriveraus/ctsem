using LinearAlgebra

# The square-system solver is shared by the exponential and intercept kernels,
# so verify the overwritten right-hand side solves the original system.
@testset "Square linear solve" begin
    A = [2.0 0.5; -0.1 1.7]
    b = [1.0, 2.0]
    A_work = copy(A)
    b_work = copy(b)

    ContinuousTimeSEM._solve_square_system!(A_work, b_work, Vector{Int}(undef, 2))
    @test A * b_work ≈ b
end
