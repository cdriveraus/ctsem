# The matrix-exponential covariance construction: correctness of the forward,
# of its pullback, and of the default route being left alone.
#
# The pullback is checked against ForwardDiff rather than by inspection. A wrong
# adjoint returns a plausible wrong gradient, which is the failure mode that has
# been most expensive in this engine, and the first version of this one was
# wrong by a factor of two in the standard-deviation term while looking right.

using ForwardDiff

const _COVEXPM_BUF = (k, T) -> (s = Vector{T}(undef, k), ss = Vector{T}(undef, k),
    r = Vector{T}(undef, k), row_sq = Vector{T}(undef, k),
    out = Matrix{T}(undef, k, k), intermediate = Matrix{T}(undef, k, k))

_covexpm_lower_n(k) = k * (k - 1) ÷ 2

function _covexpm_unpack(v, k, ::Type{T}) where {T}
    M = zeros(T, k, k)
    for i in 1:k
        M[i, i] = v[i]
    end
    idx = k
    for i in 1:k, j in 1:(i - 1)
        idx += 1
        M[i, j] = v[idx]
    end
    M
end

function _covexpm_pack(M, k)
    v = zeros(eltype(M), k + _covexpm_lower_n(k))
    for i in 1:k
        v[i] = M[i, i]
    end
    idx = k
    for i in 1:k, j in 1:(i - 1)
        idx += 1
        v[idx] = M[i, j]
    end
    v
end

# Sigma through whichever route is armed, as a function of the packed coordinates
function _covexpm_sigma(v, k)
    T = eltype(v)
    b = _COVEXPM_BUF(k, T)
    ContinuousTimeSEM.sdcovsqrt2cov!(b, _covexpm_unpack(v, k, T), 0, Val(k))
    copy(b.out)
end

function _covexpm_coords(k)
    v = zeros(k + _covexpm_lower_n(k))
    for i in 1:k
        v[i] = 0.7 + 0.3 * i / k
    end
    for m in (k + 1):length(v)
        v[m] = 0.55 * sinpi(0.6 * m / length(v))
    end
    v
end

# An asymmetric cotangent, so the symmetrisation convention is exercised rather
# than hidden by a symmetric one.
_covexpm_cotangent(k) = [0.03 * (1 + i) + 0.011 * j for i in 1:k, j in 1:k]

previous_flag = ContinuousTimeSEM.ctsem_cov_expm!(false)
try
    @testset "expm covariance forward" begin
        ContinuousTimeSEM.ctsem_cov_expm!(true)
        for k in (3, 6)
            v = _covexpm_coords(k)
            S = _covexpm_sigma(v, k)
            @test diag(S) ≈ v[1:k] .^ 2
            @test S ≈ S'
            @test minimum(eigvals(Symmetric(S))) > 0
        end
        # A lone coordinate is exactly Fisher's z of the correlation: the 2x2
        # block [[0,a],[a,0]] exponentiates to [[cosh a, sinh a],[sinh a, cosh a]],
        # which normalises to tanh(a).
        for k in (3, 6, 12), a in (0.2, 0.8, 2.0)
            v = zeros(k + _covexpm_lower_n(k))
            fill!(view(v, 1:k), 1.0)
            v[k + 1] = a
            S = _covexpm_sigma(v, k)
            @test S[2, 1] / sqrt(S[1, 1] * S[2, 2]) ≈ tanh(a) atol = 1e-12
        end
    end

    @testset "expm covariance pullback against ForwardDiff" begin
        ContinuousTimeSEM.ctsem_cov_expm!(true)
        for k in (2, 3, 6, 12)
            v = _covexpm_coords(k)
            cb = _covexpm_cotangent(k)
            analytic = ForwardDiff.gradient(x -> sum(cb .* _covexpm_sigma(x, k)), v)
            mb = zeros(k, k)
            ContinuousTimeSEM._sdcovexpm2cov_pullback!(mb,
                _covexpm_unpack(v, k, Float64), cb, k)
            @test _covexpm_pack(mb, k) ≈ analytic rtol = 1e-9
        end
    end

    @testset "the route is generic in the element type" begin
        ContinuousTimeSEM.ctsem_cov_expm!(true)
        # Float64 and BigFloat must agree. The construction calls no LAPACK, so
        # there is no fast path that only some element types take -- this pins
        # that, and would fail if one were reintroduced for BLAS floats alone.
        for k in (3, 6)
            v = _covexpm_coords(k)
            cb = _covexpm_cotangent(k)
            fast = zeros(k, k)
            ContinuousTimeSEM._sdcovexpm2cov_pullback!(fast,
                _covexpm_unpack(v, k, Float64), cb, k)
            slow = zeros(BigFloat, k, k)
            ContinuousTimeSEM._sdcovexpm2cov_pullback!(slow,
                _covexpm_unpack(BigFloat.(v), k, BigFloat), BigFloat.(cb), k)
            @test Float64.(slow) ≈ fast rtol = 1e-8
        end
    end

    @testset "the cache serves several matrices at one size" begin
        ContinuousTimeSEM.ctsem_cov_expm!(true)
        # T0VAR, DIFFUSION and MANIFESTVAR of the same size share one scratch
        # entry. With a single cache slot they evicted each other; this cycles
        # three distinct matrices repeatedly and checks every answer.
        k = 6
        mats = [_covexpm_unpack(_covexpm_coords(k) .* f, k, Float64)
                for f in (1.0, 0.6, 1.4)]
        expected = map(m -> begin
            b = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b, m, 0, Val(k))
            copy(b.out)
        end, mats)
        for _ in 1:5, (m, want) in zip(mats, expected)
            b = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b, m, 0, Val(k))
            @test b.out ≈ want
        end
        # and a matrix that has fallen out of the table still comes back right
        many = [_covexpm_unpack(_covexpm_coords(k) .* (0.5 + 0.1 * i), k, Float64)
                for i in 1:10]
        for m in many
            b = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b, m, 0, Val(k))
            b2 = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovexpm2cov!(b2, m, Val(k))
            @test b.out ≈ b2.out
        end
    end

    @testset "the default route is untouched" begin
        ContinuousTimeSEM.ctsem_cov_expm!(false)
        for k in (3, 6)
            v = _covexpm_coords(k)
            M = _covexpm_unpack(v, k, Float64)
            b = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b, M, 0, Val(k))
            b2 = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.constraincorsqrt1_vec!(b2, M, 1e-5, Val(k))
            D = Diagonal([M[i, i] for i in 1:k])
            @test b.out ≈ (D * b2.out) * (D * b2.out)'
            # and the existing pullback still agrees with ForwardDiff, so the
            # branch added to it did not disturb the arithmetic
            cb = _covexpm_cotangent(k)
            analytic = ForwardDiff.gradient(x -> sum(cb .* _covexpm_sigma(x, k)), v)
            mb = zeros(k, k)
            ContinuousTimeSEM._sdcovsqrt2cov_pullback!(mb, M, cb, k)
            @test _covexpm_pack(mb, k) ≈ analytic rtol = 1e-9
        end
    end

    @testset "choleskymats selects the construction" begin
        # `choleskymats == 2` is covmattransform='z', the same code the stan
        # path reads from its data block. This pins that the argument is
        # honoured rather than accepted and ignored, which is the failure mode
        # this engine has had before with covmattransform.
        ContinuousTimeSEM.ctsem_cov_expm!(false)
        for k in (3, 6)
            v = _covexpm_coords(k)
            M = _covexpm_unpack(v, k, Float64)
            b0 = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b0, M, 0, Val(k))
            b2 = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(b2, M, 2, Val(k))
            @test !isapprox(b0.out, b2.out)
            # 0 still reproduces constraincorsqrt1 exactly
            bc = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.constraincorsqrt1_vec!(bc, M, 1e-5, Val(k))
            D = Diagonal([M[i, i] for i in 1:k])
            @test b0.out ≈ (D * bc.out) * (D * bc.out)'
            # and 2 matches what the flag forces
            ContinuousTimeSEM.ctsem_cov_expm!(true)
            bf = _COVEXPM_BUF(k, Float64)
            ContinuousTimeSEM.sdcovsqrt2cov!(bf, M, 0, Val(k))
            ContinuousTimeSEM.ctsem_cov_expm!(false)
            @test b2.out ≈ bf.out
        end
    end

    @testset "flag toggling reports the previous value" begin
        @test ContinuousTimeSEM.ctsem_cov_expm() == false
        @test ContinuousTimeSEM.ctsem_cov_expm!(true) == false
        @test ContinuousTimeSEM.ctsem_cov_expm() == true
        @test ContinuousTimeSEM.ctsem_cov_expm!(false) == true
    end
finally
    # Leaking this flag would silently reroute every later test file's
    # covariance construction.
    ContinuousTimeSEM.ctsem_cov_expm!(previous_flag)
end
