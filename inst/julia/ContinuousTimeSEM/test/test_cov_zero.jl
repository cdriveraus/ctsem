# Zero covariances under covmattransform 'z' (cov_zero.jl): a cell fixed at
# zero between connected variables is solved so that the covariance is zero.
#
# The forward is checked against what it promises (exact zeros, plain 'z'
# where there is nothing to solve, the same matrix whatever the order of the
# variables), its derivatives against differences of the Float64 forward, and
# the reverse pass against forward mode, at one and two orders.

using LinearAlgebra
using ForwardDiff

const CT = ContinuousTimeSEM

_covzero_buf(k, T) = (s = Vector{T}(undef, k), ss = Vector{T}(undef, k),
    r = Vector{T}(undef, k), row_sq = Vector{T}(undef, k),
    out = Matrix{T}(undef, k, k), intermediate = Matrix{T}(undef, k, k))

# The open cells' positions for a pattern: every cell below the diagonal that
# is not a zero, down each column.
function _covzero_open(k, cells)
    closed = Set(zip(cells.rows, cells.cols))
    [(i, j) for j in 1:k for i in (j + 1):k if !((i, j) in closed)]
end

# The matrix of cells from a packed vector: k sds, then the open cells.
function _covzero_unpack(v, k, cells)
    T = eltype(v)
    M = zeros(T, k, k)
    for i in 1:k
        M[i, i] = v[i]
    end
    for (n, (i, j)) in enumerate(_covzero_open(k, cells))
        M[i, j] = v[k + n]
    end
    M
end

function _covzero_sigma(v, k, cells)
    T = eltype(v)
    b = _covzero_buf(k, T)
    CT.sdcovsqrt2cov!(b, _covzero_unpack(v, k, cells), CT.CTSEMCovCode(2, cells), k)
    copy(b.out)
end

function _covzero_pull(v, cb, k, cells)
    T = promote_type(eltype(v), eltype(cb))
    out = zeros(T, k, k)
    CT._sdcovsqrt2cov_pullback!(out, T.(_covzero_unpack(v, k, cells)), T.(cb), k;
        covmatcode=CT.CTSEMCovCode(2, cells))
    g = zeros(T, length(v))
    for i in 1:k
        g[i] = out[i, i]
    end
    for (n, (i, j)) in enumerate(_covzero_open(k, cells))
        g[k + n] = out[i, j]
    end
    g
end

# Every zero of a pattern whose variables the open cells connect.
function _covzero_cells(k, zeros_)
    closed = Set(zeros_)
    open = [(i, j) for j in 1:k for i in (j + 1):k if !((i, j) in closed)]
    CT._ctsem_linked_zeros(k, first.(zeros_), last.(zeros_), first.(open), last.(open))
end

function _covzero_point(k, cells; scale=0.6)
    n = length(_covzero_open(k, cells))
    v = zeros(k + n)
    for i in 1:k
        v[i] = 0.7 + 0.5 * i / k
    end
    for m in 1:n
        v[k + m] = scale * sinpi(0.37 * m + 0.2)
    end
    v
end

const _COVZERO_PATTERNS = [
    (3, [(2, 1)]),
    (3, [(3, 2)]),
    (4, [(3, 1), (4, 2)]),
    (5, [(2, 1), (4, 1), (5, 3), (5, 4)]),
    (6, [(3, 1), (4, 1), (5, 1), (6, 1), (4, 2), (5, 2), (6, 2), (5, 3), (6, 3), (6, 4)]),
]

@testset "which zeros are solved" begin
    # y1 ~~ y2 and y2 ~~ y3: the zero between 1 and 3 is linked.
    c = _covzero_cells(3, [(3, 1)])
    @test (c.rows, c.cols) == ([3], [1])
    @test c.id != 0
    # Two separate pairs: nothing to solve.
    c = _covzero_cells(4, [(3, 1), (4, 1), (3, 2), (4, 2)])
    @test isempty(c.rows) && c.id == 0
    # A pair and a chain of three: only the chain's zero.
    c = _covzero_cells(5, [(3, 1), (4, 1), (5, 1), (3, 2), (4, 2), (5, 2), (5, 3)])
    @test (c.rows, c.cols) == ([5], [3])
    @test_throws ArgumentError CT.CTSEMZeroCells([1], [2])
end

@testset "a solved cell is a zero covariance" begin
    for (k, zs) in _COVZERO_PATTERNS
        cells = _covzero_cells(k, zs)
        v = _covzero_point(k, cells)
        S = _covzero_sigma(v, k, cells)
        @test all(S[i, j] == 0 && S[j, i] == 0 for (i, j) in zip(cells.rows, cells.cols))
        @test diag(S) ≈ v[1:k] .^ 2
        @test S ≈ S'
        @test minimum(eigvals(Symmetric(S))) > 0
        # Held at 0 the same cells are not zero covariances.
        b = _covzero_buf(k, Float64)
        CT.sdcovsqrt2cov!(b, _covzero_unpack(v, k, cells), 2, k)
        @test maximum(abs(b.out[i, j]) for (i, j) in zip(cells.rows, cells.cols)) > 1e-3
    end
    # Three variables at unit sd, the open cells at 0.5.
    cells = _covzero_cells(3, [(2, 1)])
    S = _covzero_sigma([1.0, 1.0, 1.0, 0.5, 0.5], 3, cells)
    @test S[3, 1] ≈ 0.4306862 atol = 1e-6
    @test S[3, 2] ≈ S[3, 1] atol = 1e-12
end

@testset "nothing to solve is plain 'z'" begin
    k = 5
    none = CT.CTSEMZeroCells()
    v = _covzero_point(k, none)
    b = _covzero_buf(k, Float64)
    CT.sdcovsqrt2cov!(b, _covzero_unpack(v, k, none), 2, k)
    @test _covzero_sigma(v, k, none) == b.out
    # Zeros between separate blocks: the engine finds none linked.
    zs = [(3, 1), (4, 1), (5, 1), (3, 2), (4, 2), (5, 2)]
    @test isempty(_covzero_cells(k, zs).rows)
end

@testset "the order of the variables does not matter" begin
    k = 5
    zs = [(2, 1), (4, 1), (5, 3), (5, 4)]
    cells = _covzero_cells(k, zs)
    M = _covzero_unpack(_covzero_point(k, cells), k, cells)
    full = M + M' - Diagonal(diag(M))
    S = _covzero_sigma(_covzero_point(k, cells), k, cells)
    perm = [4, 2, 5, 1, 3]
    Mp = LowerTriangular(full[perm, perm]) |> Matrix
    zp = [(max(a, b), min(a, b)) for (a, b) in
          ((findfirst(==(i), perm), findfirst(==(j), perm)) for (i, j) in zs)]
    cellsp = _covzero_cells(k, zp)
    b = _covzero_buf(k, Float64)
    CT.sdcovsqrt2cov!(b, Mp, CT.CTSEMCovCode(2, cellsp), k)
    @test b.out ≈ S[perm, perm] atol = 1e-11
end

@testset "derivatives of the forward, against differences" begin
    for (k, zs) in _COVZERO_PATTERNS
        cells = _covzero_cells(k, zs)
        v = _covzero_point(k, cells)
        f = x -> vec(_covzero_sigma(x, k, cells))
        J = ForwardDiff.jacobian(f, v)
        h = 1e-6
        fd = reduce(hcat, [(f(v .+ h .* (1:length(v) .== l)) .-
            f(v .- h .* (1:length(v) .== l))) ./ (2h) for l in eachindex(v)])
        @test J ≈ fd atol = 1e-8
        # The solved covariances stay zero as the open cells move.
        for (i, j) in zip(cells.rows, cells.cols)
            @test all(iszero, J[(j - 1) * k + i, :])
        end
        # Second order: a nested dual against differences of the first.
        w = [0.3 + 0.1 * i - 0.05 * j for i in 1:k, j in 1:k]
        g = x -> ForwardDiff.gradient(y -> sum(w .* _covzero_sigma(y, k, cells)), x)
        H = ForwardDiff.jacobian(g, v)
        fdH = reduce(hcat, [(g(v .+ h .* (1:length(v) .== l)) .-
            g(v .- h .* (1:length(v) .== l))) ./ (2h) for l in eachindex(v)])
        @test H ≈ fdH atol = 1e-7
    end
end

@testset "the reverse pass, against forward mode" begin
    for (k, zs) in _COVZERO_PATTERNS
        cells = _covzero_cells(k, zs)
        v = _covzero_point(k, cells)
        # Asymmetric, and nonzero at the solved cells, where it must not count.
        cb = [0.03 * (1 + i) + 0.011 * j for i in 1:k, j in 1:k]
        forward = ForwardDiff.gradient(x -> sum(cb .* _covzero_sigma(x, k, cells)), v)
        @test _covzero_pull(v, cb, k, cells) ≈ forward rtol = 1e-9
        # At a dual: the reverse pass differentiated, against the forward's
        # second derivative.
        Hrev = ForwardDiff.jacobian(x -> _covzero_pull(x, cb, k, cells), v)
        Hfwd = ForwardDiff.hessian(x -> sum(cb .* _covzero_sigma(x, k, cells)), v)
        @test Hrev ≈ Hfwd rtol = 1e-7 atol = 1e-9
    end
end

@testset "at a trace jet" begin
    # A jet's first parts and its summed second part, against forward mode
    # along the same directions.
    k = 4
    cells = _covzero_cells(k, [(3, 1), (4, 2)])
    v = _covzero_point(k, cells)
    dirs = [[0.4 * sinpi(0.3 * i + 0.21 * l) for i in eachindex(v)] for l in 1:2]
    J = CT.CTSEMTraceJet{2}
    vj = [J(v[i], (dirs[1][i], dirs[2][i]), 0.0) for i in eachindex(v)]
    S = _covzero_sigma(vj, k, cells)
    val = _covzero_sigma(v, k, cells)
    along(d) = t -> vec(_covzero_sigma(v .+ t .* d, k, cells))
    first_ = [ForwardDiff.derivative(along(d), 0.0) for d in dirs]
    second = sum(ForwardDiff.derivative(t -> ForwardDiff.derivative(along(d), t), 0.0)
        for d in dirs)
    @test [s.v for s in vec(S)] ≈ vec(val) atol = 1e-12
    @test [s.d[1] for s in vec(S)] ≈ first_[1] atol = 1e-10
    @test [s.d[2] for s in vec(S)] ≈ first_[2] atol = 1e-10
    @test [s.s for s in vec(S)] ≈ second atol = 1e-9
end

@testset "the reverse pass at a trace jet" begin
    # The Laplace gradient's trace terms carry jets through the reverse pass:
    # its first parts and summed second part, against forward mode over the
    # reverse pass along the same directions.
    k = 4
    cells = _covzero_cells(k, [(3, 1), (4, 2)])
    v = _covzero_point(k, cells)
    cb = [0.03 * (1 + i) + 0.011 * j for i in 1:k, j in 1:k]
    dirs = [[0.4 * sinpi(0.3 * i + 0.21 * l) for i in eachindex(v)] for l in 1:2]
    J = CT.CTSEMTraceJet{2}
    vj = [J(v[i], (dirs[1][i], dirs[2][i]), 0.0) for i in eachindex(v)]
    g = _covzero_pull(vj, J.(cb), k, cells)
    along(d) = t -> _covzero_pull(v .+ t .* d, cb, k, cells)
    first_ = [ForwardDiff.derivative(along(d), 0.0) for d in dirs]
    second = sum(ForwardDiff.derivative(t -> ForwardDiff.derivative(along(d), t), 0.0)
        for d in dirs)
    @test [x.v for x in g] ≈ _covzero_pull(v, cb, k, cells) atol = 1e-12
    @test [x.d[1] for x in g] ≈ first_[1] atol = 1e-9
    @test [x.d[2] for x in g] ≈ first_[2] atol = 1e-9
    @test [x.s for x in g] ≈ second atol = 1e-8
end

@testset "BigFloat carries the solution to its own precision" begin
    k = 3
    cells = _covzero_cells(k, [(2, 1)])
    v = BigFloat.(_covzero_point(k, cells))
    b = _covzero_buf(k, BigFloat)
    M = _covzero_unpack(v, k, cells)
    CT.sdcovsqrt2cov!(b, M, CT.CTSEMCovCode(2, cells), k)
    @test Float64.(b.out) ≈ _covzero_sigma(Float64.(v), k, cells) rtol = 1e-13
    # The solved matrix's exponential is zero there far below double precision.
    sc = CT._expm_cov_scratch(BigFloat, k)
    A = zeros(BigFloat, k, k)
    CT._fill_hollow!(A, M, k, BigFloat)
    slot = CT._expm_cov_slot!(sc, A, k, cells)
    Y = sc.Ytab[slot]
    @test abs(Y[2, 1]) / sqrt(Y[1, 1] * Y[2, 2]) < 1e-40
end

@testset "no solution is NaN, not a wrong matrix" begin
    # Coordinates far past any correlation bound: the correlation matrix is
    # singular to machine precision and the zero cells cannot be solved.
    k = 5
    cells = _covzero_cells(k, [(2, 1), (4, 1), (5, 3), (5, 4)])
    v = _covzero_point(k, cells; scale=400.0)
    S = _covzero_sigma(v, k, cells)
    @test all(isnan, S) || all(S[i, j] == 0 for (i, j) in zip(cells.rows, cells.cols))
    cb = ones(k, k)
    g = _covzero_pull(v, cb, k, cells)
    @test all(isnan, S) == any(isnan, g)
end

@testset "the cache keeps patterns apart" begin
    # Equal cells under different zero patterns are different matrices.
    k = 3
    a = _covzero_cells(k, [(2, 1)])
    M = [1.0 0 0; 0 1.0 0; 0.5 0.5 1.0]
    for _ in 1:3
        b1 = _covzero_buf(k, Float64)
        CT.sdcovsqrt2cov!(b1, M, CT.CTSEMCovCode(2, a), k)
        b2 = _covzero_buf(k, Float64)
        CT.sdcovsqrt2cov!(b2, M, 2, k)
        @test b1.out[2, 1] == 0
        @test b2.out[2, 1] ≈ 0.1152759 atol = 1e-6
    end
end

@testset "the correlation bound reads the solved matrix" begin
    # A block with a solved cell: the wall's terms are those of the matrix the
    # construction builds, and its gradient is the derivative of its value.
    k = 3
    index = [0, 1, 2]                     # cells (2,1), (3,1), (3,2)
    block = CT.CTSEMCorBound(k, 2, index, [0.0, 0.0, 0.0], CT.CTSEMZeroCells([2], [1]))
    values = [2.9, 2.7]
    terms = CT._ctsem_corbound_zterms(block, values)
    cells = CT.CTSEMZeroCells([2], [1])
    S = _covzero_sigma([1.0, 1.0, 1.0, values...], k, cells)
    R = S ./ sqrt.(diag(S) * diag(S)')
    @test terms[1] ≈ [log(inv(R)[j, j]) for j in 1:k] rtol = 1e-8
    limit = CT._ctsem_corbound_limit()
    kappa = CT._CTSEM_CORBOUND_STEEPNESS
    wall = x -> begin
        logs = CT._ctsem_corbound_zterms(block, x)[1]
        -kappa * sum(max(l - limit, zero(l))^3 for l in logs)
    end
    @test wall(values) < 0
    g = zeros(2)
    CT._ctsem_corbound_zgradient!(g, block, values, limit, kappa)
    @test g ≈ ForwardDiff.gradient(wall, values) rtol = 1e-7
end

@testset "a model's table says which cells are solved" begin
    # DIFFUSION of three processes: (2,1) fixed at zero, (3,1) and (3,2) free.
    matrix = String[]; row = Int[]; col = Int[]; par = Int[]; value = Float64[]
    np = 0
    for j in 1:3, i in 1:3
        push!(matrix, "DIFFUSION"); push!(row, i); push!(col, j)
        if i == j || (i, j) in ((3, 1), (3, 2))
            np += 1
            push!(par, np); push!(value, NaN)
        else
            push!(par, 0); push!(value, 0.0)
        end
    end
    none = fill("", length(matrix))
    cz = CT.CTSEMCovZeros()
    CT._ctsem_set_cov_zeros!(cz, matrix, row, col, par, value, none, none, none, none)
    @test (cz.diffusion.rows, cz.diffusion.cols) == ([2], [1])
    @test isempty(cz.t0var.rows) && isempty(cz.manifestvar.rows)
    # A state dependent cell at (2,1) is open, so nothing is fixed at zero.
    pt = copy(none)
    pt[findfirst(i -> row[i] == 2 && col[i] == 1, eachindex(row))] = "state[1]"
    CT._ctsem_set_cov_zeros!(cz, matrix, row, col, par, value, none, pt, none, none)
    @test isempty(cz.diffusion.rows)
end

@testset "thirty variables" begin
    k = 30
    # Twenty zeros scattered through an otherwise free matrix.
    zs = [(i, j) for j in 1:k for i in (j + 1):k if (7 * i + 13 * j) % 21 == 0]
    cells = _covzero_cells(k, zs)
    @test length(cells.rows) == length(zs) >= 15
    v = _covzero_point(k, cells; scale=0.25)
    S = _covzero_sigma(v, k, cells)
    @test all(S[i, j] == 0 for (i, j) in zip(cells.rows, cells.cols))
    @test minimum(eigvals(Symmetric(S))) > 0
    cb = [0.01 * (1 + i) + 0.003 * j for i in 1:k, j in 1:k]
    forward = ForwardDiff.gradient(x -> sum(cb .* _covzero_sigma(x, k, cells)), v)
    @test _covzero_pull(v, cb, k, cells) ≈ forward rtol = 1e-8
    # Only neighbours covary: 406 zeros, all linked.
    zs = [(i, j) for j in 1:k for i in (j + 2):k]
    cells = _covzero_cells(k, zs)
    @test length(cells.rows) == 406
    v = _covzero_point(k, cells; scale=0.5)
    S = _covzero_sigma(v, k, cells)
    @test all(S[i, j] == 0 for (i, j) in zip(cells.rows, cells.cols))
    @test minimum(eigvals(Symmetric(S))) > 0
end
