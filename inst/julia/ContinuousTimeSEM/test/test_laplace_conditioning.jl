using DataFrames, LinearAlgebra

# The smallest eigenvalue of a unit's curvature, as the reports take it.
#
# Units wider than `_LAPLACE_EIGEN_MAXDIM` used to be reported as `NaN`, which
# hid a near-singular study unit of 93 coordinates on the SNSF pilot. They are
# now found by bisection on the inertia of `M - s I`, read off the block
# elimination. Checked here against a dense decomposition on the nested
# fixtures, with the dense cap lowered so the bisection runs on units small
# enough to decompose, and on curvatures shifted to be indefinite, which is the
# case the elimination's pivots must count rather than refuse.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

const _CT = ContinuousTimeSEM

"""The unit's curvature at its mode, its blocks, and the dense matrix."""
function _conditioning_unit(laplace, values, U)
    ctsem_laplace_evaluate(laplace, values; gradient=false)
    theta = collect(Float64, values)
    Ls = _CT._laplace_popchols(theta, laplace.spec)
    M = nothing
    _CT._laplace_parallel(laplace, [U]) do U
        M = _CT._laplace_unit_curvature(laplace, U, theta, Ls, laplace.modes[U])
        true
    end
    blocks = laplace.units.blocks[U]
    d = laplace.units.dims[U]
    return M, blocks, d, _CT._laplace_block_dense(M, blocks, d)
end

"""`M - c I` in block form."""
function _shifted_blocks(M, c)
    out = _CT.CTSEMBlockMatrix{Float64}([copy(d) for d in M.diag],
        [[copy(x) for x in row] for row in M.coupling])
    for d in out.diag
        for i in axes(d, 1); d[i, i] -= c; end
    end
    return out
end

@testset "inertia count and smallest eigenvalue, nested units" begin
    old = _CT._LAPLACE_EIGEN_MAXDIM[]
    try
        for (name, fresh) in (("two levels", _fresh_twolevel),
                              ("three levels", _fresh_threelevel))
            laplace, values = fresh()
            for U in eachindex(laplace.units.members)
                M, blocks, d, dense = _conditioning_unit(laplace, values, U)
                lam = eigvals(Symmetric(dense))
                # Counts at shifts between, below and above the spectrum, on
                # the curvature and on indefinite versions of it.
                for c in (0.0, lam[1] + 0.5, (lam[1] + lam[end]) / 2)
                    Mc = _shifted_blocks(M, c)
                    lc = lam .- c
                    for s in (-1.0, minimum(lc) - 0.1,
                              (lc[1] + lc[min(2, d)]) / 2, maximum(lc) + 0.1)
                        any(v -> abs(v - s) < 1e-8, lc) && continue
                        @test _CT._laplace_count_below(Mc, blocks, s) == count(<(s), lc)
                    end
                end
                # The smallest eigenvalue by bisection (cap forced to zero)
                # against the dense one, and on the indefinite shift.
                _CT._LAPLACE_EIGEN_MAXDIM[] = 0
                for c in (0.0, lam[1] - 1e-3, (lam[1] + lam[min(2, d)]) / 2)
                    Mc = _shifted_blocks(M, c)
                    target = lam[1] - c
                    got = _CT._laplace_min_eigenvalue(Mc, blocks, d)
                    # At exactly one (a direction no data inform, M = I there)
                    # the pivot is singular at the shift and the bisection
                    # brackets one rather than reporting Inf; both are right.
                    if target > 1 + 1e-6
                        @test got == Inf
                    elseif target >= 1 - 1e-6
                        @test got == Inf || isapprox(got, 1.0; rtol=1e-3)
                    else
                        @test isapprox(got, target; rtol=1e-3, atol=1e-9)
                    end
                end
                # The dense route still answers for a unit within the cap.
                _CT._LAPLACE_EIGEN_MAXDIM[] = old
                @test _CT._laplace_min_eigenvalue(M, blocks, d) ≈ lam[1]
            end
        end
    finally
        _CT._LAPLACE_EIGEN_MAXDIM[] = old
    end
end

@testset "conditioning reports a wide near-singular unit" begin
    old = _CT._LAPLACE_EIGEN_MAXDIM[]
    try
        laplace, values = _fresh_twolevel()
        ctsem_laplace_evaluate(laplace, values; gradient=false)
        dense_report = ctsem_laplace_conditioning(laplace)
        # Every unit "wide": the report must come from the bisection and agree.
        _CT._LAPLACE_EIGEN_MAXDIM[] = 0
        wide_report = ctsem_laplace_conditioning(laplace)
        @test !any(isnan, wide_report.min_eigenvalue)
        for U in eachindex(dense_report.min_eigenvalue)
            a = dense_report.min_eigenvalue[U]
            b = wide_report.min_eigenvalue[U]
            isinf(a) ? (@test isinf(b)) : (@test isapprox(a, b; rtol=1e-3))
        end
        @test wide_report.below_one == dense_report.below_one
        @test wide_report.near_singular == dense_report.near_singular
    finally
        _CT._LAPLACE_EIGEN_MAXDIM[] = old
    end
end
