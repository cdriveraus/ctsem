using LinearAlgebra
using ForwardDiff

# The buffered exponential is a custom allocation-conscious kernel; compare it
# against Julia's dense `exp` on a nontrivial matrix and against the identity
# case, where regressions are easy to spot.
@testset "Buffered matrix exponential" begin
    A = [0.0 1.0; -2.0 -3.0]
    Y = zeros(2, 2)
    W1 = zeros(2, 2)
    buffer = ContinuousTimeSEM.ExpBuffer(A)

    ContinuousTimeSEM.my_exp!(Y, A, W1, buffer)
    @test Y ≈ exp(A) atol = 1e-12 rtol = 1e-12

    ContinuousTimeSEM.my_exp!(Y, zeros(2, 2), W1, buffer)
    @test Y ≈ Matrix{Float64}(I, 2, 2)
end

@testset "Buffered matrix exponential ForwardDiff path" begin
    function exp_sum(x)
        T = eltype(x)
        A = T[0.0 x[1]; -2.0 -3.0]
        Y = zeros(T, 2, 2)
        W1 = zeros(T, 2, 2)
        buffer = ContinuousTimeSEM.ExpBuffer(A)

        ContinuousTimeSEM.my_exp!(Y, A, W1, buffer)
        return sum(Y)
    end

    x0 = [1.0]
    grad = ForwardDiff.gradient(exp_sum, x0)
    ε = sqrt(eps())
    finite_diff = (exp_sum([x0[1] + ε]) - exp_sum([x0[1] - ε])) / (2ε)

    @test grad[1] ≈ finite_diff rtol = 1e-6 atol = 1e-8
end
