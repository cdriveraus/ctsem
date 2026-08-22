using LinearAlgebra

# `ksolve!` works on packed upper-triangular symmetric matrices internally; this
# protects the packing order that the generated linear system expects.
@testset "Packed symmetric helpers" begin
    Q = [1.0 0.2; 0.2 2.0]
    triQ = zeros(3)
    ContinuousTimeSEM._ksolve_pack_upper!(triQ, Q)
    @test triQ == [1.0, 0.2, 2.0]

    unpacked = zeros(2, 2)
    ContinuousTimeSEM._ksolve_unpack_upper!(unpacked, triQ)
    @test unpacked == Q
end

# The square-system solver is shared by the exponential and Lyapunov kernels,
# so verify the overwritten right-hand side solves the original system.
@testset "Square linear solve" begin
    A = [2.0 0.5; -0.1 1.7]
    b = [1.0, 2.0]
    A_work = copy(A)
    b_work = copy(b)

    ContinuousTimeSEM._solve_square_system!(A_work, b_work, Vector{Int}(undef, 2))
    @test A * b_work ≈ b
end

# Validate the public Lyapunov solver by its residual rather than by the exact
# packed-system coefficients.
@testset "Lyapunov ksolve residual" begin
    A = [-0.8 0.2; -0.1 -1.1]
    Q = [0.3 0.05; 0.05 0.4]
    X = zeros(2, 2)
    O = zeros(3, 3)
    triQ = zeros(3)
    piv = Vector{Int}(undef, 3)

    ContinuousTimeSEM.ksolve!(X, A, Q, O, triQ, piv)
    @test X ≈ X'
    @test A * X + X * A' + Q ≈ zeros(2, 2) atol = 1e-10
end
