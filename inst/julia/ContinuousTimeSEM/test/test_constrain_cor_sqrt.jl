using LinearAlgebra

# The loop, vectorized, and buffered correlation transforms should agree. The
# resulting factor should also produce a positive-definite correlation matrix.
@testset "Correlation square-root constraints" begin
    mat = [
        1.2 0.0 0.0
        0.2 0.9 0.0
       -0.1 0.3 1.1
    ]

    loop_factor = ContinuousTimeSEM.constraincorsqrt1(mat)
    vector_factor = ContinuousTimeSEM.constraincorsqrt1_vec(Symmetric(mat, :L))

    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 3)
    ContinuousTimeSEM.constraincorsqrt1_vec!(buffer, mat)
    buffered_factor = copy(buffer.out)

    @test loop_factor ≈ vector_factor atol = 1e-12
    @test buffered_factor ≈ vector_factor atol = 1e-12

    corr = buffered_factor * buffered_factor'
    # The transform intentionally leaves an epsilon margin on the diagonal.
    @test diag(corr) ≈ fill(1.0 + 1e-5, 3) atol = 1e-12
    @test ContinuousTimeSEM.is_positive_definite(Matrix(corr))
end

# J9/F2: `_ctsem_corrsqrt_row` (adjoint_primitives.jl) is a hand-copied,
# statement-for-statement mirror of this file's loop body, kept separate so
# the reverse-mode pullback differentiates a plain function of one row rather
# than the buffered in-place form. Its docstring claims this agreement is
# tested; before this test it was not (test_adjoint_primitives.jl only pinned
# the pullback to its own mirror, never the mirror to the primal it claims to
# track -- see review/J9-duplicate-edge-blocks.md F2). `constrain_cor_sqrt.jl`
# carries several "TODO: rewrite for performance" markers on the code being
# mirrored, so this is a plausible future divergence site.
@testset "corrsqrt row mirror matches the buffered primal it claims to track" begin
    for (label, mat) in (
        ("moderate", [1.2 0.0 0.0; 0.2 0.9 0.0; -0.1 0.3 1.1]),
        ("near-zero off-diagonal", [1.0 0.0 0.0; 1e-6 1.0 0.0; -1e-6 1e-6 1.0]),
        ("large magnitude", [3.5 0.0 0.0; 2.4 4.1 0.0; -3.9 2.7 5.0]),
    )
        d = size(mat, 1)
        buffer = ContinuousTimeSEM._make_square_buffer(Float64, d)
        ContinuousTimeSEM.constraincorsqrt1_vec!(buffer, mat)
        for i in 1:d
            v = zeros(d)
            ContinuousTimeSEM._ctsem_symmetric_row!(v, mat, i, d)
            mirrored = ContinuousTimeSEM._ctsem_corrsqrt_row(v, i, 1e-5)
            @test isapprox(mirrored, buffer.out[i, :]; atol=1e-12, rtol=1e-12) ||
                error("corrsqrt row mirror mismatch for $label row $i")
        end
    end
end

# The allocating and buffered covariance paths are both used by higher-level
# code, so they should produce the same covariance from the same parameters.
@testset "Covariance conversion" begin
    mat = [
        1.2 0.0 0.0
        0.2 0.9 0.0
       -0.1 0.3 1.1
    ]

    cov = ContinuousTimeSEM.sdcovsqrt2cov(mat, 0)
    buffer = ContinuousTimeSEM._make_square_buffer(Float64, 3)
    ContinuousTimeSEM.sdcovsqrt2cov!(buffer, mat, 0)

    @test Matrix(cov) ≈ buffer.out atol = 1e-12
    @test ContinuousTimeSEM.is_positive_definite(Matrix(cov))
end
