using LinearAlgebra
using Random

# `my_exp!` and `my_exp_frechet!` both need the thirteenth-order Padé numerator
# coefficients, and used to carry their own literal tuple of them. The second
# said "as in `my_exp!`" -- a claim nothing checked. `_pade13_b` is now the one
# definition; these tests check the claim the comment used to make, which is
# that both routines really do compute the same approximant.
#
# Two independent references, because a shared constant is only as good as the
# values in it: `exp` from LinearAlgebra for the values themselves, and the
# published b0..b13 for the tuple.

@testset "Padé 13 coefficients are Higham (2008)'s" begin
    # Table 10.4 of Higham, Functions of Matrices (2008). Integers, so an exact
    # comparison is the right one -- these are not approximations of anything.
    expected = (
        64764752532480000, 32382376266240000, 7771770303897600,
        1187353796428800, 129060195264000, 10559470521600,
        670442572800, 33522128640, 1323241920, 40840800,
        960960, 16380, 182, 1)

    b = ContinuousTimeSEM._pade13_b(Float64)
    @test length(b) == 14
    @test all(b[i] == Float64(expected[i]) for i in 1:14)

    # The tuple has to arrive in the caller's own working type, because the
    # Dual paths use it too. A Float32 or a Dual must not silently get Float64s
    # back and widen the arithmetic around them.
    @test eltype(ContinuousTimeSEM._pade13_b(Float32)) === Float32
    @test eltype(ContinuousTimeSEM._pade13_b(Float64)) === Float64
    @test all(ContinuousTimeSEM._pade13_b(Float32)[i] == Float32(expected[i])
        for i in 1:14)
end

@testset "both exponential routines compute the same approximant" begin
    # This is the assertion the comment in frechet_exponential.jl used to make
    # in prose. Across sizes and across enough norm to force several squarings,
    # since the coefficients enter before the squaring loop and a mismatch
    # would show up scaled.
    rng = MersenneTwister(20260912)
    for n in (1, 2, 3, 5, 9), scale in (0.01, 1.0, 8.0)
        A = scale .* randn(rng, n, n)

        Y1 = Matrix{Float64}(undef, n, n)
        W1 = Matrix{Float64}(undef, n, n)   # scratch, per my_exp!'s signature
        ContinuousTimeSEM.my_exp!(Y1, A, W1, ContinuousTimeSEM.ExpBuffer{Float64}(n))

        # The Fréchet routine's primal output is the same exponential.
        Y2 = Matrix{Float64}(undef, n, n)
        L = Matrix{Float64}(undef, n, n)
        E = randn(rng, n, n)
        ContinuousTimeSEM.my_exp_frechet!(Y2, L, A, E,
            ContinuousTimeSEM.ExpFrechetBuffer{Float64}(n))

        @test isapprox(Y1, Y2; rtol = 1e-12, atol = 1e-14)
        # And both against an implementation that shares no code with either.
        @test isapprox(Y1, exp(A); rtol = 1e-9, atol = 1e-12)
    end
end
