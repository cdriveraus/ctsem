using LinearAlgebra
using Random

# The subject loop is LAPACK-free by default: the hand-written Cholesky and LU
# at every size, the packed Lyapunov solve at every size, `my_exp!` for every
# exponential. These pin the defaults and the three pieces the Laplace and
# quadrature routes now take from the engine's own factorization.

@testset "LAPACK-free defaults" begin
    @test ContinuousTimeSEM._CTSEM_SMALL_CHOLESKY[] == typemax(Int)
    @test ContinuousTimeSEM._CTSEM_LYAP_SCHUR_ABOVE[] == typemax(Int)
    @test ContinuousTimeSEM.LyapBuffer(Float64, 40) isa ContinuousTimeSEM.LyapKsolveBuffer{Float64}
end

@testset "CTSEMCholesky solves, inverse and U inverse match LinearAlgebra" begin
    rng = MersenneTwister(9)
    for d in (1, 3, 7)
        B = randn(rng, d, d); S = B * B' .+ d .* I(d)
        F = ContinuousTimeSEM._ctsem_cholesky(copy(S), d)
        @test issuccess(F)
        ref = cholesky(Symmetric(S))
        b = randn(rng, d); M = randn(rng, d, 4)
        @test isapprox(F \ b, ref \ b; rtol=1e-12)
        @test isapprox(F \ M, ref \ M; rtol=1e-12)
        @test isapprox(inv(F), inv(ref); rtol=1e-11)
        @test isapprox(ContinuousTimeSEM._ctsem_cholesky_uinv(F), inv(ref.U); rtol=1e-12)
        @test isapprox(logdet(F), logdet(ref); rtol=1e-12)
    end
    # A matrix that is not positive definite reports failure rather than throwing.
    notpd = [1.0 2.0; 2.0 1.0]
    @test !issuccess(ContinuousTimeSEM._ctsem_cholesky(copy(notpd), 2))
end

@testset "_ctsem_expm takes the buffered kernel at every size" begin
    rng = MersenneTwister(10)
    for n in (3, 20, 30)
        A = 0.2 .* randn(rng, n, n) .- 0.5 .* I(n)
        @test isapprox(ContinuousTimeSEM._ctsem_expm(A), exp(A); rtol=1e-11, atol=1e-13)
    end
end
