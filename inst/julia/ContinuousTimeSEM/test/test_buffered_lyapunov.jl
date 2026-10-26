using LinearAlgebra
using ForwardDiff
using Random

# Buffer construction chooses the route by size and scalar type. This guards
# the intended split between the packed solve and the Schur route, and that
# both answer the same equation.
@testset "Lyapunov buffer selection" begin
    CT = ContinuousTimeSEM
    @test CT.gauss_summation(5) == 15
    ksolve_buffer = CT.LyapBuffer(Float64, 2)
    @test !ksolve_buffer.schur
    # Run-time integers, not `Val`s: a size in a type compiles again for every
    # model (CLAUDE.md, "A model's types do not depend on the model").
    @test ksolve_buffer.ksolve_dim === 2
    @test ksolve_buffer.ksolve_system_dim === 3
    # One concrete type on both sides of the threshold, for the same reason:
    # the workspace's type must not depend on the diffusion block's size.
    threshold = CT._CTSEM_LYAP_SCHUR_ABOVE[]
    small = CT.LyapBuffer(Float64, threshold)
    large = CT.LyapBuffer(Float64, threshold + 1)
    @test !small.schur && large.schur
    @test typeof(small) === typeof(large)
    # The packed system is not allocated on the Schur route, nor the factors
    # on the packed one.
    @test isempty(large.ksolve_O) && isempty(small.S)
    # Every Dual stays on the packed route.
    D = ForwardDiff.Dual{Nothing,Float64,2}
    @test !CT.LyapBuffer(D, threshold + 5).schur
    # The setter moves the line for buffers built after it.
    CT.ctsem_set_lyapunov_schur_above!(3)
    @test CT.LyapBuffer(Float64, 4).schur
    CT.ctsem_set_lyapunov_schur_above!(threshold)
    @test !CT.LyapBuffer(Float64, 4).schur
    # Both routes keep their factors across calls: a second solve against the
    # same A must not refactor, a changed A must, and every answer must match
    # a fresh ksolve! to roundoff.
    let A = [-0.9 0.3 0.1; 0.2 -0.6 0.0; 0.05 -0.1 -0.8], Q1 = [0.4 0.1 0.0; 0.1 0.35 0.05; 0.0 0.05 0.5],
        Q2 = [1.0 0.2 0.1; 0.2 0.8 0.0; 0.1 0.0 0.6], A2 = A .- 0.1 .* I(3)
        fresh(Aa, Qq) = CT.ksolve!(zeros(3, 3), Aa, Qq, zeros(6, 6), zeros(6))
        for above in (typemax(Int), 2)
            CT.ctsem_set_lyapunov_schur_above!(above)
            buffer = CT.LyapBuffer(Float64, 3)
            CT.ctsem_set_lyapunov_schur_above!(threshold)
            @test buffer.schur == (above == 2)
            X = zeros(3, 3)
            CT.ctsem_reset_opcounts!()
            CT.my_lyap!(X, A, Q1, buffer)
            @test isapprox(X, fresh(A, Q1); atol=1e-13)
            CT.my_lyap!(X, A, Q2, buffer)
            @test isapprox(X, fresh(A, Q2); atol=1e-13)
            CT.my_lyap!(X, A2, Q1, buffer)
            @test isapprox(X, fresh(A2, Q1); atol=1e-13)
            counts = CT.ctsem_opcounts()
            # The fresh ksolve! calls count too: three of them, plus the
            # buffer's own two factorisations.
            @test (buffer.schur ? counts.lyap_schur : counts.lyap_ksolve - 3) == 2
        end
    end
end

@testset "Schur and packed routes agree across the threshold" begin
    CT = ContinuousTimeSEM
    rng = MersenneTwister(31)
    threshold = CT._CTSEM_LYAP_SCHUR_ABOVE[]
    try
        for n in (3, 8, threshold, threshold + 1, 16, 21)
            for rep in 1:3
                A = 0.4 .* randn(rng, n, n) .- 0.8 .* I(n)
                B = randn(rng, n, n)
                Q = B * B' + 0.1I
                CT.ctsem_set_lyapunov_schur_above!(typemax(Int))
                packed = CT.LyapBuffer(Float64, n)
                CT.ctsem_set_lyapunov_schur_above!(1)
                schur = CT.LyapBuffer(Float64, n)
                CT.ctsem_set_lyapunov_schur_above!(threshold)
                X1 = zeros(n, n); X2 = zeros(n, n)
                CT.my_lyap!(X1, A, Q, packed)
                CT.my_lyap!(X2, A, Q, schur)
                @test isapprox(X2, X1; rtol=1e-11)
                @test X2 == X2'
                @test isapprox(A * X2 + X2 * A' + Q, zeros(n, n); atol=1e-11 * opnorm(Q, 1))
                # The transposed solve on A's own factors, against the packed
                # route's solve with the transpose.
                W1 = zeros(n, n); W2 = zeros(n, n)
                CT.my_lyap_transposed!(W1, A, Q, packed, zeros(n, n))
                CT.ctsem_reset_opcounts!()
                CT.my_lyap_transposed!(W2, A, Q, schur, zeros(n, n))
                @test CT.ctsem_opcounts().lyap_schur == 0    # factors reused
                @test isapprox(W2, W1; rtol=1e-11)
                @test isapprox(A' * W2 + W2 * A + Q, zeros(n, n); atol=1e-11 * opnorm(Q, 1))
                # The series route's conditioning reads the factors: the
                # smallest real part of an eigenvalue sum.
                ev = eigvals(A)
                @test isapprox(CT._series_lyap_conditioning(schur),
                    minimum(abs(real(a) + real(b)) for a in ev, b in ev); rtol=1e-8)
            end
        end
    finally
        CT.ctsem_set_lyapunov_schur_above!(threshold)
    end
end

@testset "Schur route: blocks larger than the solve, failure, no allocation" begin
    CT = ContinuousTimeSEM
    threshold = CT._CTSEM_LYAP_SCHUR_ABOVE[]
    n = threshold + 2
    rng = MersenneTwister(32)
    A = 0.3 .* randn(rng, n, n) .- 0.7 .* I(n)
    B = randn(rng, n, n); Q = B * B' + I
    buffer = CT.LyapBuffer(Float64, n)
    @test buffer.schur
    # Buffers larger than the block: only the leading block is read and
    # written, as the packed route's unpacking does.
    Xbig = fill(7.0, n + 2, n + 3)
    Qbig = fill(NaN, n + 1, n + 1); Qbig[1:n, 1:n] = Q
    CT.my_lyap!(Xbig, A, Qbig, buffer)
    @test isapprox(A * Xbig[1:n, 1:n] + Xbig[1:n, 1:n] * A' + Q, zeros(n, n); atol=1e-10)
    @test all(==(7.0), Xbig[(n + 1):end, :]) && all(==(7.0), Xbig[:, (n + 1):end])
    # A view of a larger matrix as A, as the discretisation passes it.
    big = randn(rng, n + 3, n + 3); big[1:n, 1:n] = A
    X = zeros(n, n)
    CT.my_lyap!(X, view(big, 1:n, 1:n), Q, buffer)
    @test isapprox(A * X + X * A' + Q, zeros(n, n); atol=1e-10)
    # Nothing allocated per call, refactoring or not.
    A2 = A .- 0.01 .* I(n)
    W = zeros(n, n); At = zeros(n, n)
    CT.my_lyap!(X, A2, Q, buffer); CT.my_lyap!(X, A, Q, buffer)
    CT.my_lyap_transposed!(W, A2, Q, buffer, At)
    @test (@allocated CT.my_lyap!(X, A, Q, buffer)) == 0
    @test (@allocated CT.my_lyap!(X, A, Q, buffer)) == 0
    @test (@allocated CT.my_lyap_transposed!(W, A2, Q, buffer, At)) == 0
    # A non-finite A gives NaN, and conditioning NaN, which sends the interval
    # to the series route, rather than an error or a loop.
    Abad = copy(A); Abad[2, 3] = Inf
    CT.my_lyap!(X, Abad, Q, buffer)
    @test all(isnan, X)
    @test isnan(CT._series_lyap_conditioning(buffer))
    @test CT._series_needed(CT._series_lyap_conditioning(buffer), 1.0)
    # And the next finite A is solved again.
    CT.my_lyap!(X, A, Q, buffer)
    @test isapprox(A * X + X * A' + Q, zeros(n, n); atol=1e-10)
end

# Rather than compare against one implementation detail, validate the defining
# Lyapunov equation residual for both solver paths.
@testset "Buffered Lyapunov residuals" begin
    for n in (2, 5, 14)
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
