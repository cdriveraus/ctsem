using LinearAlgebra

# A nested unit's softest eigenpairs, and draws with covariance its inverse,
# from its block storage alone (`_laplace_block_lowest`, `_laplace_block_draw`,
# SAEM's `_saem_block_proposal`), checked against dense decompositions on block
# matrices wide enough for the structured route: persons under a study, and
# persons under groups under a study, as the mode search and SAEM's start meet
# them on data with many members a unit.

const _CB = ContinuousTimeSEM

"""A random symmetric positive definite matrix of size `k`, well conditioned."""
function _blocks_spd(k, rng_state)
    X = [sin(rng_state * 1.3 + 0.7 * i * j + 0.11 * i * i) for i in 1:k, j in 1:k]
    return Matrix(Symmetric(X' * X + 0.5I))
end

"""
A unit of `nleaf` blocks of `leaf` under a root of `root` (two levels), or, with
`nmid > 0`, the leaves split evenly under `nmid` blocks of `mid` under the root
(three levels): its blocks innermost first and its block matrix.
"""
function _blocks_nested(nleaf; leaf=4, root=4, nmid=0, mid=2, coupling=0.15)
    blocks = _CB.CTSEMLaplaceBlock[]
    nb = nleaf + nmid + 1
    rootat = nb
    off = 0
    for i in 1:nleaf
        anc = nmid > 0 ? [nleaf + 1 + (i - 1) % nmid, rootat] : [rootat]
        push!(blocks, _CB.CTSEMLaplaceBlock(off, leaf, [i], 1, anc))
        off += leaf
    end
    for m in 1:nmid
        push!(blocks, _CB.CTSEMLaplaceBlock(off, mid, Int[], 2, [rootat]))
        off += mid
    end
    push!(blocks, _CB.CTSEMLaplaceBlock(off, root, Int[], nmid > 0 ? 3 : 2, Int[]))
    M = _CB.CTSEMBlockMatrix(Float64, blocks)
    for (b, block) in enumerate(blocks)
        M.diag[b] .= _blocks_spd(block.size, b) +
            (b == nb ? nleaf * coupling^2 * 4 : 0.0) * I
        for (t, a) in enumerate(block.ancestors)
            M.coupling[b][t] .= [coupling * cos(0.9 * b + 1.7 * i + 2.3 * j + t)
                                 for i in 1:block.size, j in 1:blocks[a].size]
        end
    end
    return M, blocks, off + root
end

"""`M - c I` in block form."""
function _blocks_shift(M, c)
    out = _CB.CTSEMBlockMatrix{Float64}([copy(d) for d in M.diag],
        [[copy(x) for x in row] for row in M.coupling])
    for d in out.diag
        for i in axes(d, 1); d[i, i] -= c; end
    end
    return out
end

@testset "a nested unit's softest eigenpairs from its blocks" begin
    for (name, M0, blocks, d) in (("two levels", _blocks_nested(30)...),
                                  ("three levels", _blocks_nested(24; leaf=3, nmid=4)...))
        dense0 = _CB._laplace_block_dense(M0, blocks, d)
        lam0 = eigvals(Symmetric(dense0))
        @test lam0[1] > 0
        # Shifted so a few eigenvalues sit low, and one so the matrix is
        # indefinite, which the count and the iteration must take as they come.
        for c in (lam0[1] - 0.05, lam0[2] + 1e-3, (lam0[3] + lam0[4]) / 2)
            M = _blocks_shift(M0, c)
            E = eigen(Symmetric(dense0 - c * I))
            for s in (E.values[1] - 0.1, (E.values[2] + E.values[3]) / 2,
                      (E.values[5] + E.values[6]) / 2, 0.0, 1.0)
                any(v -> abs(v - s) < 1e-6, E.values) && continue
                k = count(<(s), E.values)
                low = _CB._laplace_block_lowest(M, blocks, s)
                @test low !== nothing
                low === nothing && continue
                @test length(low.values) == k
                @test isapprox(low.values, E.values[1:k]; rtol=1e-8, atol=1e-10)
                V = low.vectors
                @test size(V) == (d, k)
                @test isapprox(V' * V, Matrix(1.0I, k, k); atol=1e-10)
                # The subspace, to what the residual allows over the gap at
                # `threshold` (two clustered leaves' pairs sit 4e-4 apart here).
                Vd = E.vectors[:, 1:k]
                gap = k == 0 ? 1.0 : E.values[k + 1] - E.values[k]
                @test maximum(abs.(V * V' .- Vd * Vd'); init=0.0) <=
                    10 * sqrt(max(k, 1)) * 1e-9 / min(gap, 1.0)
                # Each a column of `M`'s eigenvectors, to the residual it stops at.
                for j in 1:k
                    r = _CB._soft_block_mulvec(M, blocks, V[:, j]) .- low.values[j] .* V[:, j]
                    @test norm(r) <= 1e-8 * max(1.0, abs(low.values[j]))
                end
            end
        end
        # The oracle declines, so the callers decompose densely.
        previous = ctsem_set_block_eigen!(false)
        try
            @test _CB._laplace_block_lowest(M0, blocks, 1.0) === nothing
        finally
            ctsem_set_block_eigen!(previous)
        end
    end
end

@testset "draws from a nested unit's blocks have the inverse covariance" begin
    for (name, M0, blocks, d) in (("two levels", _blocks_nested(30)...),
                                  ("three levels", _blocks_nested(24; leaf=3, nmid=4)...))
        dense0 = _CB._laplace_block_dense(M0, blocks, d)
        lam0 = eigvals(Symmetric(dense0))
        # One with no eigenvalue under the prior's, and one with several.
        for c in (lam0[1] - 1.5, lam0[1] - 0.2)
            M = _blocks_shift(M0, c)
            dense = dense0 - c * I
            F = _CB._soft_shifted_factor(M, blocks, 0.0)
            @test F !== nothing
            Id = Matrix{Float64}(I, d, d)
            # Through the identity's columns, so `X X'` is the covariance.
            X = reduce(hcat, [_CB._laplace_block_draw(F.factors, F.coupling, blocks,
                Id[:, j]) for j in 1:d])
            @test isapprox(X * X', inv(Symmetric(dense)); rtol=1e-9)
            P = _CB._saem_block_proposal(M, blocks)
            @test P !== nothing
            Y = reduce(hcat, [_CB._saem_block_draw(P, blocks, Id[:, j]) for j in 1:d])
            E = eigen(Symmetric(dense))
            clipped = E.vectors * Diagonal(1 ./ max.(E.values, 1.0)) * E.vectors'
            @test isapprox(Y * Y', clipped; rtol=1e-8)
            # What SAEM's dense start draws from, to the same covariance.
            Fd = _CB._saem_proposal_factor(dense)
            @test isapprox(Fd * Fd', clipped; rtol=1e-8)
        end
    end
end

@testset "a nested unit's indefinite elimination counts, solves and inverts" begin
    for (name, M0, blocks, d) in (("two levels", _blocks_nested(30)...),
                                  ("three levels", _blocks_nested(24; leaf=3, nmid=4)...))
        dense0 = _CB._laplace_block_dense(M0, blocks, d)
        lam0 = eigvals(Symmetric(dense0))
        for c in (lam0[1] - 0.3, lam0[2] + 1e-3, (lam0[5] + lam0[6]) / 2)
            M = _blocks_shift(M0, c)
            A = Symmetric(dense0 - c * I)
            for s in (0.0, 0.2)
                F = _CB._laplace_block_ldl(M, blocks, s)
                @test F !== nothing
                @test F.negative == count(<(s), eigvals(A))
                @test isapprox(F.logabsdet, logabsdet(Matrix(A) - s * I)[1]; rtol=1e-10)
                x = [sin(0.3 * i) for i in 1:d]
                @test isapprox(_CB._laplace_block_ldl_solve(F, blocks, x),
                    (Matrix(A) - s * I) \ x; rtol=1e-9)
            end
            # A node's factor: Cholesky where definite, the elimination's
            # pivots otherwise, which the block solve and the selected inverse
            # take alike.
            F = _CB._soft_node_factor(M, blocks)
            Ainv = inv(A)
            x = [cos(0.7 * i) for i in 1:d]
            @test isapprox(_CB._soft_solve(F, blocks, x), Ainv * x; rtol=1e-9)
            @test isapprox(F.logdet, logabsdet(Matrix(A))[1]; rtol=1e-10)
            Cd, Cc = _CB._laplace_selected_inverse(F.factors, F.coupling, blocks)
            for (b, blk) in enumerate(blocks)
                rb = (blk.offset + 1):(blk.offset + blk.size)
                @test isapprox(Cd[b], Ainv[rb, rb]; rtol=1e-8, atol=1e-11)
                for (t, a) in enumerate(blk.ancestors)
                    ra = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
                    @test isapprox(Cc[b][t], Ainv[rb, ra]; rtol=1e-8, atol=1e-11)
                end
            end
        end
    end
end

@testset "a soft node's floored complement and its high-space solves" begin
    top = _CB._laplace_soft_node_curvature_floor_top()
    for (name, M0, blocks, d) in (("two levels", _blocks_nested(30)...),
                                  ("three levels", _blocks_nested(24; leaf=3, nmid=4)...))
        dense0 = _CB._laplace_block_dense(M0, blocks, d)
        lam0 = eigvals(Symmetric(dense0))
        v = normalize([1 + 0.5 * sin(1.3 * i) for i in 1:d])
        Vh = nullspace(transpose(v))
        # Definite with a few pairs under the top, and indefinite.
        for c in (lam0[1] - 0.05, lam0[2] + 1e-3)
            M = _blocks_shift(M0, c)
            A = Symmetric(dense0 - c * I)
            E = eigen(Symmetric(transpose(Vh) * A * Vh))
            k = count(<(top), E.values)
            low = _CB._soft_complement_lowest(M, blocks, v, top)
            @test low !== nothing
            low === nothing && continue
            @test length(low.values) == k
            @test isapprox(low.values, E.values[1:k]; rtol=1e-8, atol=1e-10)
            W = low.vectors
            @test maximum(abs, transpose(W) * v; init=0.0) < 1e-10
            for j in 1:k
                r = _CB._soft_off(Vector(A * W[:, j]), v) .- low.values[j] .* W[:, j]
                @test norm(r) <= 1e-8 * max(1.0, abs(low.values[j]))
            end
            # (C_H - mu)^-1 on the complement of v and the low pairs.
            Z = hcat(v, W)
            Ph = I - Z * transpose(Z)
            e = [cos(0.4 * i) for i in 1:d]
            for j in 1:k
                r = _CB._soft_high_solve(M, blocks, Z, low.values[j], e)
                @test r !== nothing
                @test maximum(abs, transpose(Z) * r; init=0.0) < 1e-9
                @test isapprox(Ph * (Matrix(A) * r .- low.values[j] .* r), Ph * e;
                    rtol=1e-8, atol=1e-10)
            end
        end
    end
end

@testset "the map's root from a nested unit's block elimination" begin
    for (name, M, blocks, d) in (("two levels", _blocks_nested(30)...),
                                 ("three levels", _blocks_nested(24; leaf=3, nmid=4)...))
        dense = _CB._laplace_block_dense(M, blocks, d)
        ok, _, factors, coupling = _CB._laplace_block_factor(M, blocks)
        @test ok
        Id = Matrix{Float64}(I, d, d)
        # Q = K K': K^-T by columns, K^-1 its transpose, K' its inverse.
        X = reduce(hcat, [_CB._laplace_block_draw(factors, coupling, blocks, Id[:, j])
            for j in 1:d])
        Y = reduce(hcat, [_CB._laplace_block_lsolve(factors, coupling, blocks, Id[:, j])
            for j in 1:d])
        R = reduce(hcat, [_CB._laplace_block_rmul(factors, coupling, blocks, Id[:, j])
            for j in 1:d])
        @test isapprox(X * transpose(X), inv(Symmetric(dense)); rtol=1e-9)
        @test isapprox(Y, transpose(X); rtol=1e-9, atol=1e-12)
        @test isapprox(R * X, Id; atol=1e-9)
        @test isapprox(transpose(R) * R, dense; rtol=1e-10)
        @test isapprox(_CB._laplace_block_halflogdet(factors, blocks),
            logdet(Symmetric(dense)) / 2; rtol=1e-12)
    end
end

@testset "the soft rule's centre takes a binding floor from the unit's blocks" begin
    # A nested unit whose curvature at the mode has distinct eigenvalues under
    # the floor's top -- one leaf, then two, given low ones -- on a cubic
    # integrand: the factored rule's value against the dense rule's, which
    # floors the centre's complement eigenvalue by eigenvalue. The objective's
    # laplace object only lends its gate.
    isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
        include(joinpath(@__DIR__, "laplace_fixtures.jl"))
    laplace, _ = _fresh_twolevel()
    ctsem_set_soft_directions!(laplace, :quadrature; lo=0.2, hi=0.7)
    top = _CB._laplace_soft_node_curvature_floor_top()
    Q4 = Matrix(qr([sin(1.1 * i + 0.7 * j) for i in 1:4, j in 1:4]).Q)
    for lows in ([0.08, 0.2], [0.08, 0.2, 0.26])
        Ab, blocks, d = _blocks_nested(30)
        Ab.diag[1] .= Q4 * Diagonal([lows[1], lows[2], 1.5, 2.5]) * transpose(Q4)
        if length(lows) > 2
            Ab.diag[2] .= Q4 * Diagonal([lows[3], 1.2, 1.8, 2.2]) * transpose(Q4)
        end
        A = Symmetric(_CB._laplace_block_dense(Ab, blocks, d))
        E = eigen(A)
        lam = E.values
        nlow = count(<(top), lam)
        @test lam[1] > 0
        @test nlow == length(lows)
        @test lam[2] - lam[1] > 0.05
        u = zeros(d)
        k3 = [0.02 * sin(0.9 * i) for i in 1:d]
        objf = z -> (value=-dot(z, A * z) / 2 + sum(k3 .* z .^ 3),
            gradient=-(A * z) .+ 3 .* k3 .* z .^ 2, hlogdet=0.0)
        curvf = z -> begin
            M = _blocks_shift(Ab, 0.0)
            for (b, blk) in enumerate(blocks), i in 1:blk.size
                M.diag[b][i, i] -= 6 * k3[blk.offset + i] * z[blk.offset + i]
            end
            M
        end
        ok, logdetM, _, _ = _CB._laplace_block_factor(Ab, blocks)
        @test ok
        centre = _CB._soft_centre_floor(Ab, blocks, E.vectors[:, 1], lam[1], top)
        @test centre !== nothing
        centre === nothing && continue
        @test isapprox(centre.values, lam[2:nlow]; rtol=1e-8)
        fv = _CB._laplace_soft_quadrature_factored(laplace, Ab, blocks, logdetM, 0.0, u,
            objf, curvf)
        @test fv !== nothing
        @test fv !== nothing && !fv.logdet_floor
        previous = ctsem_set_soft_quadrature_factored!(false)
        dv = try
            _CB._laplace_soft_quadrature_term(laplace, 1, Float64[], Matrix{Float64}[], u,
                nothing; M=Ab, logdetM=logdetM, inner=objf(u), objective=objf,
                curvature=curvf, blocks=blocks)
        finally
            ctsem_set_soft_quadrature_factored!(previous)
        end
        fv === nothing || @test isapprox(fv.value, dv.value; rtol=1e-11)
    end
end
