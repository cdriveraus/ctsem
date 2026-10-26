using LinearAlgebra
using ForwardDiff
using Random

# `CTSEMTraceJet` carries a value, `W` directional derivatives and the sum of
# the second directional derivatives along them. The references here are
# ForwardDiff's: the gradient's projections for the first parts and
# `sum_j q_j' H q_j` from the Hessian for the summed second part.

const _JT = ContinuousTimeSEM.CTSEMTraceJet
const _JSF = ContinuousTimeSEM.SpecialFunctions

# The jet of `f` at `x0` along the columns of `Q` (at most the width).
function _jet_eval(f, x0::Vector{Float64}, Q::Matrix{Float64}, W::Int)
    x = [_JT{W}(x0[i], ntuple(j -> j <= size(Q, 2) ? Q[i, j] : 0.0, Val(W)), 0.0)
         for i in eachindex(x0)]
    return f(x)
end

@testset "Trace jet: scalar rules against ForwardDiff" begin
    rng = MersenneTwister(71)
    CT = ContinuousTimeSEM
    fs = [
        x -> x[1] * x[2] + x[3] / x[1] - sqrt(x[2]),
        x -> exp(x[1] - x[2]) * log(x[3]) + log1p(x[1]^2) - expm1(-x[2]),
        x -> CT.log1p_exp(x[1] + 2x[2]) - CT.inv_logit(x[3] - x[1]),
        x -> x[1]^3 - x[2]^2.5 + x[3]^-2 + inv(x[1] + x[3]) + abs(x[1] - x[2]),
        x -> _JSF.erfc(x[1] - x[3]) + _JSF.loggamma(x[2] + 1) + tanh(x[1]) * cos(x[2]),
        x -> x[2]^x[1] + 2.0^x[3] - (x[1] * x[2] * x[3]) / (1 + x[2]^2),
    ]
    x0 = [0.7, 1.3, 0.9]
    for f in fs, nd in (1, 3, 8)
        Q = randn(rng, 3, nd)
        r = _jet_eval(f, x0, Q, 8)
        g = ForwardDiff.gradient(f, x0)
        H = ForwardDiff.hessian(f, x0)
        @test r.v ≈ f(x0) rtol = 1e-13
        for j in 1:8
            @test r.d[j] ≈ (j <= nd ? dot(g, Q[:, j]) : 0.0) atol = 1e-12 rtol = 1e-12
        end
        @test r.s ≈ sum(dot(Q[:, j], H * Q[:, j]) for j in 1:nd) rtol = 1e-10 atol = 1e-12
    end
    # Equality reads every component, as the engine's caches need; a constant
    # jet equals its value.
    a = _JT{2}(1.0, (0.5, 0.0), 0.0)
    @test a != _JT{2}(1.0, (0.0, 0.0), 0.0)
    @test _JT{2}(1.0) == 1.0
    @test !(a == 1.0)
    @test !iszero(_JT{2}(0.0, (1.0, 0.0), 0.0))
    @test CT._finite_deep(_JT{2}(1.0, (NaN, 0.0), 0.0)) == false
    @test CT._primal(a) == 1.0
    @test all(isnan, (CT._ctsem_nan(_JT{2}).v, CT._ctsem_nan(_JT{2}).s))
end

@testset "Trace jet: the Lyapunov solve" begin
    CT = ContinuousTimeSEM
    rng = MersenneTwister(72)
    kron_lyap(A, Q) = (n = size(A, 1); reshape(-((kron(Matrix{eltype(A)}(I, n, n), A) +
        kron(A, Matrix{eltype(A)}(I, n, n))) \ vec(Q)), n, n))
    for n in (1, 3, 6), transposed in (false, true)
        A0 = 0.3 .* randn(rng, n, n) .- 0.8 .* I(n)
        B = randn(rng, n, n); Q0 = B * B' + I
        W = 8
        nd = 3
        dA = [randn(rng, n, n) for _ in 1:nd]
        dQ = [let C = randn(rng, n, n); C + C' end for _ in 1:nd]
        sA = randn(rng, n, n)
        sQ = let C = randn(rng, n, n); C + C' end
        J = _JT{W}
        A = [J(A0[i, j], ntuple(p -> p <= nd ? dA[p][i, j] : 0.0, Val(W)), sA[i, j])
             for i in 1:n, j in 1:n]
        Q = [J(Q0[i, j], ntuple(p -> p <= nd ? dQ[p][i, j] : 0.0, Val(W)), sQ[i, j])
             for i in 1:n, j in 1:n]
        X = zeros(J, n, n)
        buffer = CT.LyapBuffer(J, n)
        transposed ? CT.my_lyap_transposed!(X, A, Q, buffer) : CT.my_lyap!(X, A, Q, buffer)
        # Reference: along each direction p the path A(t) = A0 + t dA_p + t^2 sA/(2 nd),
        # likewise Q, so that the summed second derivative of the inputs is sA
        # and sQ; X's summed second derivative is the sum over p of each path's.
        Aref(p, t) = (M = A0 .+ t .* dA[p] .+ (t^2 / (2nd)) .* sA; transposed ? Matrix(M') : M)
        Qref(p, t) = Q0 .+ t .* dQ[p] .+ (t^2 / (2nd)) .* sQ
        X0 = kron_lyap(Aref(1, 0.0), Q0)
        @test isapprox(getfield.(X, :v), X0; rtol=1e-10)
        S = zeros(n, n)
        for p in 1:nd
            Xt(t) = kron_lyap(Aref(p, t), Qref(p, t))
            d1 = ForwardDiff.derivative(Xt, 0.0)
            d2 = ForwardDiff.derivative(t -> ForwardDiff.derivative(Xt, t), 0.0)
            @test isapprox([x.d[p] for x in X], d1; rtol=1e-9, atol=1e-11)
            S .+= d2
        end
        @test isapprox(getfield.(X, :s), S; rtol=1e-9, atol=1e-11)
    end
end
