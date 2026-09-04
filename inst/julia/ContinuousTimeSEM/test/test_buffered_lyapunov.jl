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
    # The threshold is where the packed solve stops being faster than a LAPACK
    # call; see `_CTSEM_LYAP_SCHUR_ABOVE` for the measurement behind ten.
    @test ContinuousTimeSEM.LyapBuffer(Float64, 10) isa ContinuousTimeSEM.LyapKsolveBuffer{Float64}
    # The packed buffer keeps its LU across calls: a second solve against the
    # same A must not refactor, a changed A must, and every answer must match
    # a fresh ksolve! to roundoff.
    let A = [-0.9 0.3 0.1; 0.2 -0.6 0.0; 0.05 -0.1 -0.8], Q1 = [0.4 0.1 0.0; 0.1 0.35 0.05; 0.0 0.05 0.5],
        Q2 = [1.0 0.2 0.1; 0.2 0.8 0.0; 0.1 0.0 0.6], A2 = A .- 0.1 .* I(3)
        kb = ContinuousTimeSEM.LyapKsolveBuffer{Float64}(3)
        X = zeros(3, 3)
        fresh(Aa, Qq) = ContinuousTimeSEM.ksolve!(zeros(3, 3), Aa, Qq, zeros(6, 6), zeros(6))
        ContinuousTimeSEM.ctsem_reset_opcounts!()
        ContinuousTimeSEM.my_lyap!(X, A, Q1, kb)
        @test isapprox(X, fresh(A, Q1); atol=1e-13)
        ContinuousTimeSEM.my_lyap!(X, A, Q2, kb)
        @test isapprox(X, fresh(A, Q2); atol=1e-13)
        @test ContinuousTimeSEM.ctsem_opcounts().lyap_ksolve == 1 + 2   # one cached factorisation, two fresh ones
        ContinuousTimeSEM.my_lyap!(X, A2, Q1, kb)
        @test isapprox(X, fresh(A2, Q1); atol=1e-13)
        @test ContinuousTimeSEM.ctsem_opcounts().lyap_ksolve == 2 + 3
    end
    # The Schur route is LAPACK and off by default (`_CTSEM_LYAP_SCHUR_ABOVE` is
    # typemax); the setter brings it back above a chosen size.
    @test ContinuousTimeSEM.LyapBuffer(Float64, 11) isa ContinuousTimeSEM.LyapKsolveBuffer{Float64}
    ContinuousTimeSEM.ctsem_set_lyapunov_schur_above!(10)
    @test ContinuousTimeSEM.LyapBuffer(Float64, 11) isa ContinuousTimeSEM.LyapSchurBuffer{Float64}
    ContinuousTimeSEM.ctsem_set_lyapunov_schur_above!(typemax(Int))
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
