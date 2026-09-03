using LinearAlgebra
using ForwardDiff
using Random

# `my_exp_frechet!` is the Al-Mohy & Higham recurrence; the block identity
# `exp([A E; 0 A])` in `_ctsem_exp_frechet_block_reference` is an independent
# route to the same Padé approximant, so the two must agree to roundoff -- and
# at every size class the linear solve dispatches on (hand-written LU at or
# below `_CTSEM_SMALL_CHOLESKY`, LAPACK above it), and with enough norm to force
# several squarings.

function _frechet_pair(A, E)
    n = size(A, 1)
    T = promote_type(eltype(A), eltype(E))
    Y = Matrix{T}(undef, n, n)
    L = Matrix{T}(undef, n, n)
    ContinuousTimeSEM.my_exp_frechet!(Y, L, A, E, ContinuousTimeSEM.ExpFrechetBuffer{T}(n))
    return Y, L
end

@testset "Fréchet recurrence matches the block identity" begin
    rng = MersenneTwister(2026)
    for n in (1, 2, 3, 6, 9, 20), scale in (0.05, 1.0, 30.0)
        A = scale .* (randn(rng, n, n) .- 0.5 .* I(n))
        E = randn(rng, n, n)
        Y, L = _frechet_pair(A, E)
        Yref = zeros(n, n)
        ContinuousTimeSEM.my_exp!(Yref, copy(A), zeros(n, n), ContinuousTimeSEM.ExpBuffer{Float64}(n), Val(n))
        Lref = ContinuousTimeSEM._ctsem_exp_frechet_block_reference(A, E)
        tol = 1e-11 * max(1.0, opnorm(Yref, 1))
        @test isapprox(Y, Yref; atol=tol, rtol=1e-11)
        @test isapprox(L, Lref; atol=1e-11 * max(1.0, opnorm(Lref, 1)), rtol=1e-11)
    end
end

@testset "Fréchet recurrence: linearity and transpose identity" begin
    rng = MersenneTwister(7)
    A = randn(rng, 4, 4) .- 0.7 .* I(4)
    E1 = randn(rng, 4, 4); E2 = randn(rng, 4, 4)
    _, L1 = _frechet_pair(A, E1)
    _, L2 = _frechet_pair(A, E2)
    _, L12 = _frechet_pair(A, 2.0 .* E1 .- 3.0 .* E2)
    @test isapprox(L12, 2.0 .* L1 .- 3.0 .* L2; atol=1e-12, rtol=1e-12)
    # <Ȳ, L(A, E)> = <L(A', Ȳ), E>, which is the identity the adjoint rests on.
    Ȳ = randn(rng, 4, 4)
    _, Lt = _frechet_pair(transpose(A), Ȳ)
    @test isapprox(dot(Ȳ, L1), dot(Lt, E1); rtol=1e-12)
    # The engine's wrapper is the same computation.
    @test isapprox(ContinuousTimeSEM._ctsem_exp_frechet_block(A, E1), L1; atol=1e-13)
    @test isapprox(ContinuousTimeSEM._ctsem_exp_frechet_adjoint(A, Ȳ), Lt; atol=1e-13)
end

@testset "Fréchet recurrence differentiates under ForwardDiff" begin
    # `ctsem_hessian` runs the adjoint at Dual, so the kernel must carry
    # partials through the LU and the squarings. Compare the derivative of
    # L(A + εU, E) in ε against a central difference of the Float64 kernel.
    rng = MersenneTwister(11)
    A = randn(rng, 3, 3) .- 0.5 .* I(3)
    U = randn(rng, 3, 3)
    E = randn(rng, 3, 3)
    f(ε) = ContinuousTimeSEM._ctsem_exp_frechet_block(A .+ ε .* U, E)
    dual = ForwardDiff.derivative(f, 0.0)
    h = 1e-6
    fd = (f(h) .- f(-h)) ./ (2h)
    @test isapprox(dual, fd; atol=1e-7, rtol=1e-7)
    # Dual entries in the direction too: the partials of E must reach L.
    g(ε) = ContinuousTimeSEM._ctsem_exp_frechet_block(A, E .+ ε .* U)
    @test isapprox(ForwardDiff.derivative(g, 0.0), ContinuousTimeSEM._ctsem_exp_frechet_block(A, U);
        atol=1e-12, rtol=1e-12)
end

@testset "Fréchet recurrence: non-finite input gives NaN, not an error" begin
    A = [NaN 0.1; 0.2 -0.5]
    Y, L = _frechet_pair(A, ones(2, 2))
    @test all(isnan, Y) && all(isnan, L)
    Y, L = _frechet_pair([-0.5 0.1; 0.2 -0.5], [Inf 0.0; 0.0 0.0])
    @test all(isnan, L)
end
