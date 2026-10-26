using LinearAlgebra
using ForwardDiff
using Random

# The Lyapunov solve is the Schur route at every size and element type. Its
# reference here is the Kronecker system `(I ⊗ A + A ⊗ I) vec(X) = -vec(Q)`,
# solved by LinearAlgebra's generic LU, which takes Float64, BigFloat and
# Duals alike: an independent derivation, so agreement says the factors, the
# quasi triangular solve and the tangents are right rather than that two of
# the engine's routes agree.
_kron_lyap(A, Q) = (n = size(A, 1); reshape(-((kron(Matrix{eltype(A)}(I, n, n), A) +
    kron(A, Matrix{eltype(A)}(I, n, n))) \ vec(Q)), n, n))

# Drifts the subject loop meets: a random stable one, a repeated eigenvalue
# (diagonal, already in Schur form), a defective, non-normal one (a Jordan
# block), and an oscillator (a 2x2 block, complex pair).
function _lyap_drifts(rng, n)
    out = Any[0.4 .* randn(rng, n, n) .- 0.8 .* I(n)]
    push!(out, -Matrix(Diagonal(fill(0.7, n))))
    J = -0.6 .* Matrix(1.0I, n, n)
    for i in 1:(n - 1)
        J[i, i + 1] = 1.3
    end
    push!(out, J)
    if n >= 2
        O = 0.2 .* randn(rng, n, n) .- 0.5 .* I(n)
        O[1:2, 1:2] = [-0.05 2.0; -2.0 -0.05]
        push!(out, O)
    end
    return out
end

@testset "Lyapunov buffer: one type at every size, factors kept" begin
    CT = ContinuousTimeSEM
    # Run-time sizes, one concrete type: a size in a type compiles the filter
    # again for every model (CLAUDE.md, "A model's types do not depend on the
    # model").
    @test typeof(CT.LyapBuffer(Float64, 2)) === typeof(CT.LyapBuffer(Float64, 30))
    D = ForwardDiff.Dual{Nothing,Float64,2}
    @test typeof(CT.LyapBuffer(D, 2)) === typeof(CT.LyapBuffer(D, 30))
    @test CT.LyapBuffer(D, 3).tangents.primal isa CT.LyapBuffer{Float64}
    @test CT.LyapBuffer(Float64, 0).dim == 0     # no diffusing state
    # A second solve against the same A does not refactor, a changed A does.
    A = [-0.9 0.3 0.1; 0.2 -0.6 0.0; 0.05 -0.1 -0.8]
    Q1 = [0.4 0.1 0.0; 0.1 0.35 0.05; 0.0 0.05 0.5]
    Q2 = [1.0 0.2 0.1; 0.2 0.8 0.0; 0.1 0.0 0.6]
    A2 = A .- 0.1 .* I(3)
    buffer = CT.LyapBuffer(Float64, 3)
    X = zeros(3, 3)
    CT.ctsem_reset_opcounts!()
    CT.my_lyap!(X, A, Q1, buffer)
    @test isapprox(X, _kron_lyap(A, Q1); atol=1e-14)
    CT.my_lyap!(X, A, Q2, buffer)
    @test isapprox(X, _kron_lyap(A, Q2); atol=1e-14)
    CT.my_lyap!(X, A2, Q1, buffer)
    @test isapprox(X, _kron_lyap(A2, Q1); atol=1e-14)
    @test CT.ctsem_opcounts().lyap_schur == 2
end

@testset "Schur solve against the Kronecker system, every size" begin
    CT = ContinuousTimeSEM
    rng = MersenneTwister(31)
    for n in (1, 2, 3, 5, 8, 9, 16, 21), A in _lyap_drifts(rng, n)
        B = randn(rng, n, n)
        Q = B * B' + 0.1I
        buffer = CT.LyapBuffer(Float64, n)
        X = zeros(n, n)
        CT.my_lyap!(X, A, Q, buffer)
        @test isapprox(X, _kron_lyap(A, Q); rtol=1e-10)
        @test X == X'
        # A backward-stable solve's residual is relative to the terms it
        # cancels, `|A||X|`, not to `Q` alone.
        scale(Y) = n * (2 * opnorm(A, 1) * opnorm(Y, 1) + opnorm(Q, 1))
        @test opnorm(A * X + X * A' + Q, 1) <= 1e-13 * scale(X)
        # The transposed solve on A's own factors.
        W = zeros(n, n)
        CT.ctsem_reset_opcounts!()
        CT.my_lyap_transposed!(W, A, Q, buffer)
        @test CT.ctsem_opcounts().lyap_schur == 0    # factors reused
        @test isapprox(W, _kron_lyap(Matrix(A'), Q); rtol=1e-10)
        @test opnorm(A' * W + W * A + Q, 1) <= 1e-13 * scale(W)
        # The series route's conditioning reads the factors: the smallest
        # real part of an eigenvalue sum.
        ev = eigvals(A)
        @test isapprox(CT._series_lyap_conditioning(buffer),
            minimum(abs(real(a) + real(b)) for a in ev, b in ev); rtol=1e-6)
    end
end

@testset "Schur route: blocks larger than the solve, failure, no allocation" begin
    CT = ContinuousTimeSEM
    for n in (2, 10)
        rng = MersenneTwister(32 + n)
        A = 0.3 .* randn(rng, n, n) .- 0.7 .* I(n)
        B = randn(rng, n, n); Q = B * B' + I
        buffer = CT.LyapBuffer(Float64, n)
        # Buffers larger than the block: only the leading block is read and
        # written.
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
        W = zeros(n, n)
        CT.my_lyap!(X, A2, Q, buffer); CT.my_lyap!(X, A, Q, buffer)
        CT.my_lyap_transposed!(W, A2, Q, buffer)
        @test (@allocated CT.my_lyap!(X, A, Q, buffer)) == 0
        @test (@allocated CT.my_lyap!(X, A, Q, buffer)) == 0
        @test (@allocated CT.my_lyap_transposed!(W, A2, Q, buffer)) == 0
        # A non-finite A gives NaN, and conditioning NaN, which sends the
        # interval to the series route, rather than an error or a loop.
        Abad = copy(A); Abad[2, 1] = Inf
        CT.my_lyap!(X, Abad, Q, buffer)
        @test all(isnan, X)
        @test isnan(CT._series_lyap_conditioning(buffer))
        @test CT._series_needed(CT._series_lyap_conditioning(buffer), 1.0)
        # And the next finite A is solved again.
        CT.my_lyap!(X, A, Q, buffer)
        @test isapprox(A * X + X * A' + Q, zeros(n, n); atol=1e-10)
    end
end

@testset "Duals: tangents on the value's factors" begin
    CT = ContinuousTimeSEM
    for n in (1, 2, 3, 11)
        rng = MersenneTwister(33 + n)
        A0 = 0.3 .* randn(rng, n, n) .- 0.7 .* I(n)
        B = randn(rng, n, n); Q0 = B * B' + I
        dA = [randn(rng, n, n) for _ in 1:3]
        dQ = [let C = randn(rng, n, n); C + C' end for _ in 1:3]
        D = ForwardDiff.Dual{Nothing,Float64,3}
        A = [D(A0[i, j], ForwardDiff.Partials(ntuple(p -> dA[p][i, j], 3))) for i in 1:n, j in 1:n]
        Q = [D(Q0[i, j], ForwardDiff.Partials(ntuple(p -> dQ[p][i, j], 3))) for i in 1:n, j in 1:n]
        buffer = CT.LyapBuffer(D, n)
        for transposed in (false, true)
            X = zeros(D, n, n)
            CT.ctsem_reset_opcounts!()
            transposed ? CT.my_lyap_transposed!(X, A, Q, buffer) : CT.my_lyap!(X, A, Q, buffer)
            # One factorisation of the values serves all three partials, and
            # the transposed solve reuses the forward one's.
            @test CT.ctsem_opcounts().lyap_schur == (transposed ? 0 : 1)
            # ForwardDiff through the generic LU of the Kronecker system.
            ref = _kron_lyap(transposed ? Matrix(transpose(A)) : A, Q)
            @test isapprox(ForwardDiff.value.(X), ForwardDiff.value.(ref); rtol=1e-10)
            for p in 1:3
                @test isapprox(ForwardDiff.partials.(X, p), ForwardDiff.partials.(ref, p);
                    rtol=1e-9, atol=1e-11)
            end
        end
    end
    # A nested Dual, as a Hessian over the adjoint runs: second derivatives
    # of a scalar function of the solution against differences of first
    # derivatives.
    n = 4
    rng = MersenneTwister(40)
    A0 = 0.3 .* randn(rng, n, n) .- 0.7 .* I(n)
    B = randn(rng, n, n); Q0 = B * B' + I
    dA1 = randn(rng, n, n); dQ1 = let C = randn(rng, n, n); C + C' end
    f(t) = begin
        T = typeof(t)
        At = T.(A0) .+ t .* T.(dA1)
        Qt = T.(Q0) .+ t .* T.(dQ1)
        X = zeros(T, n, n)
        CT.my_lyap!(X, At, Qt, CT.LyapBuffer(T, n))
        sum(abs2, X)
    end
    d2 = ForwardDiff.derivative(t -> ForwardDiff.derivative(f, t), 0.0)
    h = 1e-5
    fd = (ForwardDiff.derivative(f, h) - ForwardDiff.derivative(f, -h)) / (2h)
    @test isapprox(d2, fd; rtol=1e-6)
end

# The adjoint's BigFloat checks (CLAUDE.md, "An adjoint can be exactly right
# and numerically wrong") run the filter in BigFloat, so the route has to
# solve to that precision, not to Float64's.
@testset "BigFloat: the Schur route at its own precision" begin
    CT = ContinuousTimeSEM
    setprecision(BigFloat, 256) do
        rng = MersenneTwister(41)
        for n in (1, 2, 5), A64 in _lyap_drifts(rng, n)
            A = BigFloat.(A64)
            B = BigFloat.(randn(rng, n, n)); Q = B * B' + I
            X = zeros(BigFloat, n, n)
            CT.my_lyap!(X, A, Q, CT.LyapBuffer(BigFloat, n))
            @test maximum(abs, X - _kron_lyap(A, Q)) <= 1e-60 * maximum(abs, X)
            @test maximum(abs, A * X + X * A' + Q) <=
                1e-60 * (maximum(abs, A) * maximum(abs, X) + maximum(abs, Q))
        end
        # A Dual over BigFloat: tangents on BigFloat factors.
        n = 3
        A0 = BigFloat.(0.3 .* randn(rng, n, n) .- 0.7 .* I(n))
        B = BigFloat.(randn(rng, n, n)); Q0 = B * B' + I
        dA = BigFloat.(randn(rng, n, n))
        D = ForwardDiff.Dual{Nothing,BigFloat,1}
        A = [D(A0[i, j], ForwardDiff.Partials((dA[i, j],))) for i in 1:n, j in 1:n]
        Q = D.(Q0)
        X = zeros(D, n, n)
        CT.my_lyap!(X, A, Q, CT.LyapBuffer(D, n))
        ref = _kron_lyap(A, Q)
        @test maximum(abs, ForwardDiff.partials.(X, 1) - ForwardDiff.partials.(ref, 1)) <=
            1e-60 * maximum(abs, ForwardDiff.partials.(ref, 1))
    end
end

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
