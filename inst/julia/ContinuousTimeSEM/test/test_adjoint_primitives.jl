using ForwardDiff
using ChainRulesCore
using LinearAlgebra

# The primitive layer (`src/adjoint_primitives.jl`) carries the only
# hand-derived reverse-mode mathematics in the package, so it gets the most
# direct possible check: for each primitive, contract its pullback with a
# cotangent and compare against the same contraction computed by ForwardDiff
# through the primal. A `<v, J u>` vs `<J' v, u>` agreement over random `u`
# and `v` catches transposition and sign errors that a symmetric test case
# would hide, which is why the test matrices below are deliberately
# non-symmetric and the cotangents random.

const _PRIM_TOL = 1e-9

"""
Compare a primitive's pullback against ForwardDiff.

`primal(args...)` must return the primitive's output. The pullback is
contracted with `cotangent`, and each returned input cotangent is checked
against the corresponding ForwardDiff directional derivative: for input `k`
and random direction `u`, `<cotangent, d primal / d args[k] * u>` must equal
`<pullback_k, u>`.
"""
function _check_pullback(primal, args, cotangent; tol=_PRIM_TOL, symmetric_args=())
    _, pullback = ChainRulesCore.rrule(primal, args...)
    pulled = pullback(cotangent)[2:end]
    for k in eachindex(args)
        u = randn(size(args[k]))
        # Arguments the primal only reads symmetrically (the Lyapunov `Q`)
        # must be perturbed symmetrically: the derivative in a non-symmetric
        # direction is a property of the packing convention inside `ksolve!`,
        # not of the mathematical map, and every caller in the filter supplies
        # a covariance.
        k in symmetric_args && (u = (u .+ transpose(u)) ./ 2)
        directional = ForwardDiff.derivative(
            ε -> dot(cotangent, primal(ntuple(j -> j == k ? args[j] .+ ε .* u : args[j], length(args))...)),
            0.0,
        )
        @test isapprox(dot(pulled[k], u), directional; atol=tol, rtol=tol)
    end
end

@testset "matrix exponential adjoint" begin
    A = [-0.7 0.4 0.1; 0.2 -0.5 -0.3; -0.1 0.15 -0.9]
    Ȳ = randn(3, 3)
    _check_pullback(ContinuousTimeSEM._ctsem_expm, (A,), Ȳ)

    # exp(A) itself must agree with the package's buffered forward kernel,
    # otherwise the pullback would be differentiating a different function
    # from the one the filter evaluates.
    Y = zeros(3, 3)
    ContinuousTimeSEM.my_exp!(Y, copy(A), zeros(3, 3), ContinuousTimeSEM.ExpBuffer{Float64}(3), Val(3))
    @test isapprox(Y, ContinuousTimeSEM._ctsem_expm(A); atol=1e-12, rtol=1e-12)
end

@testset "Lyapunov solve adjoint" begin
    # Stable A (negative real parts) so the Lyapunov solution exists and is
    # well conditioned, and a genuinely non-diagonal Q.
    A = [-0.9 0.3; 0.2 -0.6]
    Q = [0.4 0.12; 0.12 0.35]
    X = ContinuousTimeSEM._ctsem_lyap(A, Q)
    @test isapprox(A * X + X * transpose(A) + Q, zeros(2, 2); atol=1e-12)
    @test isapprox(X, transpose(X); atol=1e-14)

    # `X` is symmetric by construction, so only the symmetric part of a
    # cotangent is meaningful; the pullback symmetrises internally, and this
    # asserts it does (an asymmetric cotangent must give the same answer as
    # its symmetric part, not a different one).
    cotangent = [0.7 -0.3; 0.15 0.9]
    symmetric = (cotangent .+ transpose(cotangent)) ./ 2
    _check_pullback(ContinuousTimeSEM._ctsem_lyap, (A, Q), symmetric; symmetric_args=(2,))
    asym_pull = ChainRulesCore.rrule(ContinuousTimeSEM._ctsem_lyap, A, Q)[2](cotangent)
    sym_pull = ChainRulesCore.rrule(ContinuousTimeSEM._ctsem_lyap, A, Q)[2](symmetric)
    @test isapprox(asym_pull[2], sym_pull[2]; atol=1e-12)
    @test isapprox(asym_pull[3], sym_pull[3]; atol=1e-12)
end

@testset "linear solve adjoint" begin
    A = [1.4 0.3 -0.2; 0.1 1.1 0.4; -0.3 0.2 1.7]
    B = [0.5 -0.2; 1.1 0.3; -0.4 0.9]
    @test isapprox(A * ContinuousTimeSEM._ctsem_linsolve(A, B), B; atol=1e-12)
    _check_pullback(ContinuousTimeSEM._ctsem_linsolve, (A, B), randn(3, 2))

    b = [0.5, 1.1, -0.4]
    _check_pullback(ContinuousTimeSEM._ctsem_linsolve, (A, b), randn(3))
end

@testset "Frechet block identity" begin
    # L(A, E) is linear in E, and L(A, I) for diagonal A reduces to
    # exp(A) elementwise on the diagonal -- a case that can be checked by hand.
    A = [-0.5 0.0; 0.0 -0.25]
    L = ContinuousTimeSEM._ctsem_exp_frechet_block(A, Matrix{Float64}(I, 2, 2))
    @test isapprox(diag(L), exp.(diag(A)); atol=1e-12)

    B = [-0.5 0.3; 0.1 -0.8]
    E1, E2 = randn(2, 2), randn(2, 2)
    @test isapprox(
        ContinuousTimeSEM._ctsem_exp_frechet_block(B, 2.0 .* E1 .- 3.0 .* E2),
        2.0 .* ContinuousTimeSEM._ctsem_exp_frechet_block(B, E1) .-
        3.0 .* ContinuousTimeSEM._ctsem_exp_frechet_block(B, E2);
        atol=1e-10,
    )
end

@testset "correlation-sqrt row pullback matches AD" begin
    # `_ctsem_corrsqrt_row_pullback!` is hand-derived for speed; this pins it
    # to `_ctsem_corrsqrt_row`, the readable reference implementation it
    # replaced, via ForwardDiff. Random rows (including a near-degenerate one
    # with a large row sum, which is where the softening term matters) are
    # checked for every row index.
    for (label, mat) in (("moderate", [0.0 0.0 0.0; 0.35 0.0 0.0; -0.2 0.45 0.0]),
                         ("large row sum", [0.0 0.0 0.0; 2.4 0.0 0.0; 1.9 2.7 0.0]),
                         ("near zero", [0.0 0.0 0.0; 1e-6 0.0 0.0; -1e-6 1e-6 0.0]))
        d = size(mat, 1)
        for i in 1:d
            v = zeros(d)
            ContinuousTimeSEM._ctsem_symmetric_row!(v, mat, i, d)
            obar = randn(d)
            vbar = zeros(d)
            ContinuousTimeSEM._ctsem_corrsqrt_row_pullback!(vbar, v, obar, i, 1e-5)
            reference = ForwardDiff.gradient(
                u -> dot(obar, ContinuousTimeSEM._ctsem_corrsqrt_row(u, i, 1e-5)), v)
            @test isapprox(vbar, reference; atol=1e-9, rtol=1e-9) ||
                error("row pullback mismatch for $label row $i")
        end
    end
end

@testset "cache guard depends on ForwardDiff comparing partials" begin
    # `_blocks_identical` (discrete_time_form.jl) decides "these inputs are
    # unchanged, reuse the cached factorization" using plain `==`. Under
    # ForwardDiff those entries are `Dual`s, so the caches are only correct
    # because ForwardDiff defines `==` on same-tag duals as
    # `value == value && partials == partials`.
    #
    # That is a dependency on another package's semantics, and if it ever
    # became value-only the caches would silently reuse a matrix exponential
    # whose derivative belongs to a different point -- a wrong gradient on the
    # default ForwardDiff path, with nothing to signal it. This asserts the
    # property directly, so such a change fails loudly here.
    D = ForwardDiff.Dual{Nothing,Float64,2}
    a = D(1.5, ForwardDiff.Partials((1.0, 0.0)))
    b = D(1.5, ForwardDiff.Partials((0.0, 1.0)))
    same = D(1.5, ForwardDiff.Partials((1.0, 0.0)))

    @test ForwardDiff.value(a) == ForwardDiff.value(b)   # same value...
    @test a != b                                          # ...but not equal
    @test a == same

    # ...and through the block comparison the caches actually call.
    A = [a a; a a]
    B = [a a; a b]
    @test ContinuousTimeSEM._blocks_identical(A, copy(A), 2, 2)
    @test !ContinuousTimeSEM._blocks_identical(A, B, 2, 2)
    # Differing only outside the compared block must still count as identical.
    @test ContinuousTimeSEM._blocks_identical(A, B, 1, 2)
end

