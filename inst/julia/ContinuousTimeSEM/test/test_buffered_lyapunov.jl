using LinearAlgebra
using ForwardDiff

# Buffer dispatch chooses different algorithms by size and scalar type. This
# guards the intended split between the direct packed solver and Schur solver.
@testset "Lyapunov buffer selection" begin
    @test ContinuousTimeSEM.gauss_summation(5) == 15
    ksolve_buffer = ContinuousTimeSEM.LyapBuffer(Float64, 2)
    @test ksolve_buffer isa ContinuousTimeSEM.LyapKsolveBuffer{Float64}
    @test ksolve_buffer.ksolve_dim isa Val{2}
    @test ksolve_buffer.ksolve_system_dim isa Val{3}
    @test ContinuousTimeSEM.LyapBuffer(Float64, 5) isa ContinuousTimeSEM.LyapSchurBuffer{Float64}
end

# Rather than compare against one implementation detail, validate the defining
# Lyapunov equation residual for both solver paths.
@testset "Buffered Lyapunov residuals" begin
    for n in (2, 5)
        A = -Matrix(Diagonal(collect(1.0:n)))
        Q = [1.0 / (i + j) for i in 1:n, j in 1:n]
        Q = Matrix(Symmetric(Q))
        X = zeros(n, n)
        buffer = ContinuousTimeSEM.LyapBuffer(Float64, n)

        ContinuousTimeSEM.my_lyap!(X, A, Q, buffer)
        @test X ≈ X' atol = 1e-12
        @test A * X + X * A' + Q ≈ zeros(n, n) atol = 1e-10
    end
end

@testset "Buffered Lyapunov ForwardDiff path" begin
    function lyap_sum(x)
        T = eltype(x)
        A = T[-x[1] 0.2; -0.1 -1.1]
        Q = T[0.3 0.05; 0.05 0.4]
        X = zeros(T, 2, 2)
        buffer = ContinuousTimeSEM.LyapBuffer(T, 2)

        ContinuousTimeSEM.my_lyap!(X, A, Q, buffer)
        return sum(X)
    end

    x0 = [0.8]
    grad = ForwardDiff.gradient(lyap_sum, x0)
    ε = sqrt(eps())
    finite_diff = (lyap_sum([x0[1] + ε]) - lyap_sum([x0[1] - ε])) / (2ε)

    @test grad[1] ≈ finite_diff rtol = 1e-6 atol = 1e-8
end
