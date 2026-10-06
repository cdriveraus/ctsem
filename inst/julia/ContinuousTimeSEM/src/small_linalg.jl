"""
Cholesky factorization and solves for the small dense matrices this engine
actually works with, written out rather than delegated to LAPACK.

# Why

LAPACK acquires an internal buffer from a process-global pool on every call,
under a lock, and `BLAS.set_num_threads(1)` does not change that. For a
`potrf` on a 4x4 that lock *is* the cost -- there are about thirty flops to
amortize it over. Measured on 23 cores, 20,000 factorizations split across
threads:

    n                     1      2      4      8
    LAPACK, 23 threads  0.30x  0.11x  0.14x  0.20x
    this file, 23       5.86x  8.47x  6.74x  10.87x
    this file, serial   11x    3.2x   2.1x   1.7x   faster than LAPACK

So it is not a trade: the unblocked loop is faster on one core *and* it scales,
because it touches nothing outside its own arguments. Above roughly sixteen
states LAPACK's blocking starts to earn its overhead back and this stops being
the right call -- `_CTSEM_SMALL_CHOLESKY` is where that line is drawn.

This is the reason the engine's threading was a net loss. The filter does at
least one Cholesky per row and one matrix exponential per prediction substep,
so a subject sweep is hundreds of LAPACK calls, and 23 threads spent their time
queueing for that lock rather than filtering.

# What

`_ctsem_cholesky!` factorizes the leading `d x d` block of `A` in place into an
upper triangular `U` with `U' U = A`, reading and writing only the upper
triangle -- the same convention `LinearAlgebra.cholesky!(::Matrix)` uses, so
the factor is laid out where the rest of the code already looks for it.
`CTSEMCholesky` wraps it and answers `ldiv!`, `rdiv!`, `issuccess`, `logdet`,
`.factors` and `.L`, which is everything the filter asks a factorization for.
"""

using LinearAlgebra

"""
State dimension up to which the hand-written factorization is used.

Above it LAPACK's blocking has enough arithmetic to amortize its call overhead
and its lock, and a plain triple loop starts losing badly -- ctsem models run to
a hundred latent states in extreme cases, and an unblocked `n^3` loop at that
size is not what anyone wants. The default is nevertheless *never LAPACK*,
because the lock serialises every thread and the price of the loops is small:
measured on dev1 (single thread, minimum of 3000 repeats, microseconds)

    n            4     8    16    24    32    48    64
    Cholesky hand 0.07 0.13  0.46  1.18  2.56  8.39 20.6
    Cholesky LAPACK 0.12 0.22 0.63 1.22 2.09 6.75 14.8
    LU solve hand 0.10 0.34  1.54  3.97  7.92 21.6  49.1
    LU solve LAPACK 0.20 0.44 1.26 2.51 4.53 9.78 18.8

so the unblocked loops win to about 24 states and cost at most 1.4x
(Cholesky) and 2.6x (LU) at 64 -- tens of microseconds per row.
`ctsem_set_small_linalg!(dimension = 24)` brings LAPACK back above that size
for a session that runs one fit on one core and wants it.
"""
const _CTSEM_SMALL_CHOLESKY = Ref(typemax(Int))

"""
Arithmetic budget above which a matrix product goes to BLAS `gemm` rather than
the engine's own kernels, in `rows * cols * inner`. Never, by default.

OpenBLAS's `gemm` takes a process-wide lock for its work buffer on every call
on most CPUs, so threads filtering different subjects queue for it. Measured
with eight threads calling at once, one BLAS thread (scripts/
sampler-checks-2026-10/kernels.jl): a 4x4 to 32x32 `gemm` costs about 2.2 us a
call locally (i9, AVX2) and 8.5 us on dev2 (EPYC 7702, AVX2) whatever its size
-- 30 and 50 times its serial cost at 4x4. Only CPUs where OpenBLAS has its
small-matrix kernels skip the lock, and then for `A*B` and `A*B'` only: dev1
(EPYC 9654, AVX-512) does, which is where an earlier threshold of 64 was
measured and looked right. 3x3 and smaller never reach OpenBLAS.

The engine's kernels (`_ctsem_mul!` and its transposed forms) take no lock and
come within about 1.5x of `gemm`'s serial speed from 8x8 up, so the threshold
is a setter for a single-threaded run on a large model, not a default.
"""
const _CTSEM_SMALL_PRODUCT = Ref(typemax(Int))

"""
Arithmetic budget above which a matrix-vector product goes to BLAS `gemv`, in
`length(y) * inner`. BLAS's matrix-vector products take no buffer and no lock:
with eight threads calling at once they cost what they cost serially (12x12:
33 against 35 ns locally, 67 against 101 on dev2), so this stays where the
hand loop stops being quicker.
"""
const _CTSEM_SMALL_MATVEC = Ref(64)

export ctsem_set_small_linalg!
"""
    ctsem_set_small_linalg!(; dimension, product, matvec)

Move the thresholds separating the engine's own small-matrix kernels from
LAPACK and BLAS. Exposed because the crossover is a property of the machine, not
of the mathematics -- the same reason `ctsem_set_block_threshold!` exists.
"""
function ctsem_set_small_linalg!(; dimension::Integer=_CTSEM_SMALL_CHOLESKY[],
    product::Integer=_CTSEM_SMALL_PRODUCT[], matvec::Integer=_CTSEM_SMALL_MATVEC[])
    dimension >= 0 && product >= 0 && matvec >= 0 ||
        throw(ArgumentError("thresholds must be non-negative"))
    _CTSEM_SMALL_CHOLESKY[] = Int(dimension)
    _CTSEM_SMALL_PRODUCT[] = Int(product)
    _CTSEM_SMALL_MATVEC[] = Int(matvec)
    return (dimension=Int(dimension), product=Int(product), matvec=Int(matvec))
end

"""
    _ctsem_cholesky!(A, d)

Factorize the leading `d x d` block of `A` in place, upper triangular, and
report whether it was positive definite. `A` is read as symmetric through its
upper triangle and is left partly overwritten when the answer is `false`.
"""
function _ctsem_cholesky!(A::AbstractMatrix{T}, d::Int) where {T}
    @inbounds for j in 1:d
        s = A[j, j]
        for k in 1:(j - 1)
            s -= A[k, j] * A[k, j]
        end
        # `>` and not `>=`: a zero pivot is a singular matrix, and dividing by
        # it below would put Inf in the factor rather than reporting failure.
        s > zero(real(T)) || return false
        u = sqrt(s)
        A[j, j] = u
        for i in (j + 1):d
            t = A[j, i]
            for k in 1:(j - 1)
                t -= A[k, j] * A[k, i]
            end
            A[j, i] = t / u
        end
    end
    return true
end

"""
    CTSEMCholesky(U, d, ok)

An upper triangular factor `U` with `U' U = S`, over the leading `d x d` block.

Immutable and holding only a view, so constructing one per row costs nothing:
it never escapes the function that makes it.
"""
struct CTSEMCholesky{T,M<:AbstractMatrix{T}}
    U::M
    d::Int
    ok::Bool
end

"""Factorize `A`'s leading `d x d` block and wrap the result."""
@inline function _ctsem_cholesky(A::AbstractMatrix, d::Integer)
    ok = _ctsem_cholesky!(A, Int(d))
    return CTSEMCholesky(A, Int(d), ok)
end

LinearAlgebra.issuccess(F::CTSEMCholesky) = F.ok

function Base.getproperty(F::CTSEMCholesky, name::Symbol)
    name === :factors && return getfield(F, :U)
    # The filter only asks for `.L` when generating data, which is not the path
    # this file exists for; a triangular wrapper over the transpose is the right
    # answer and the cost of getting there does not matter.
    name === :L && return LowerTriangular(transpose(getfield(F, :U)))
    return getfield(F, name)
end

function LinearAlgebra.logdet(F::CTSEMCholesky{T}) where {T}
    total = zero(real(T))
    U = F.U
    @inbounds for i in 1:F.d
        total += log(U[i, i])
    end
    return 2 * total
end

"""`y .= S \\ b`, by forward then back substitution through `U`."""
function LinearAlgebra.ldiv!(y::AbstractVector, F::CTSEMCholesky, b::AbstractVector)
    U = F.U
    d = F.d
    @inbounds for i in 1:d
        t = b[i]
        for k in 1:(i - 1)
            t -= U[k, i] * y[k]
        end
        y[i] = t / U[i, i]
    end
    @inbounds for i in d:-1:1
        t = y[i]
        for k in (i + 1):d
            t -= U[i, k] * y[k]
        end
        y[i] = t / U[i, i]
    end
    return y
end

LinearAlgebra.ldiv!(F::CTSEMCholesky, b::AbstractVector) = ldiv!(b, F, b)

"""`F \\ b`, allocating: the solve the Laplace block elimination asks for."""
function Base.:\(F::CTSEMCholesky, b::AbstractVector)
    y = similar(b, promote_type(eltype(F.U), eltype(b)))
    return ldiv!(y, F, b)
end

function Base.:\(F::CTSEMCholesky, B::AbstractMatrix)
    Y = similar(B, promote_type(eltype(F.U), eltype(B)))
    @inbounds for j in axes(B, 2)
        ldiv!(view(Y, :, j), F, view(B, :, j))
    end
    return Y
end

"""The inverse of the factorized matrix, dense: `d` solves against the identity."""
Base.inv(F::CTSEMCholesky{T}) where {T} = F \ Matrix{T}(LinearAlgebra.I, F.d, F.d)

"""
    _ctsem_cholesky_uinv(F)

`U^-1` for the upper factor, by back substitution column by column -- what
`inv(F.U)` would give through LAPACK's `trtri`, without the call.
"""
_ctsem_cholesky_uinv(F::CTSEMCholesky{T}) where {T} =
    _ctsem_cholesky_uinv!(zeros(T, F.d, F.d), F)

"""`X` set to `U^-1` (upper triangular, zero below) for `F = U'U`; `X` at least `d x d`."""
function _ctsem_cholesky_uinv!(X::AbstractMatrix{T}, F::CTSEMCholesky{T}) where {T}
    d = F.d
    U = F.U
    @inbounds for j in 1:d, i in (j + 1):d
        X[i, j] = zero(T)
    end
    @inbounds for j in 1:d
        X[j, j] = one(T) / U[j, j]
        for i in (j - 1):-1:1
            acc = zero(T)
            for k in (i + 1):j
                acc -= U[i, k] * X[k, j]
            end
            X[i, j] = acc / U[i, i]
        end
    end
    return X
end

"""`X .= X / S`, one row at a time: `X U^-1` then `X U^-T`."""
function LinearAlgebra.rdiv!(X::AbstractMatrix, F::CTSEMCholesky)
    U = F.U
    d = F.d
    @inbounds for r in axes(X, 1)
        for i in 1:d
            t = X[r, i]
            for k in 1:(i - 1)
                t -= U[k, i] * X[r, k]
            end
            X[r, i] = t / U[i, i]
        end
        for i in d:-1:1
            t = X[r, i]
            for k in (i + 1):d
                t -= U[i, k] * X[r, k]
            end
            X[r, i] = t / U[i, i]
        end
    end
    return X
end

################################################################################
# Small matrix products
################################################################################
#
# Matrix products in per-row and per-substep code do not call BLAS: `gemm`
# locks on most CPUs (`_CTSEM_SMALL_PRODUCT` has the measurements), and a
# transposed view is worse still, because Julia cannot hand it to `gemm` and
# falls back to a copy or the generic kernel. These touch nothing outside their
# own arguments, so they thread.
#
# Two kernels. Below an inner dimension of `_CTSEM_COLUMN_KERNEL`, a dot product
# per element, which is quickest when everything is a few registers. From there
# up, each column of `C` built as a running sum of columns of `A`, four of them
# a pass (`_ctsem_colmul!`): contiguous in the row index, so it vectorises, and
# a quarter of the passes over `C`. At 12x12 that is 241 ns against 497 for the
# dot form and 167 for an unlocked `gemm` (local), 430 against 883 and 465 on
# dev2. The dot form read `A` along a row, which is strided.
#
# All six follow `mul!(C, A, B, alpha, beta)`: `C = alpha * op(A) * op(B) +
# beta * C`, with `beta = 0` overwriting rather than reading `C` -- so an
# uninitialised buffer is safe.

"""
Is this product to be done by the engine's kernels rather than by BLAS?

`inner` is the contracted dimension, which differs between the transposed forms
-- hence passing it rather than reading it off one operand.
"""
@inline _ctsem_small_product(C, inner::Integer) =
    length(C) * inner <= _CTSEM_SMALL_PRODUCT[]

"""The matrix-vector form of `_ctsem_small_product`, against `_CTSEM_SMALL_MATVEC`."""
@inline _ctsem_small_matvec(y, inner::Integer) =
    length(y) * inner <= _CTSEM_SMALL_MATVEC[]

"""Inner dimension from which the column kernel (`_ctsem_colmul!`) is used."""
const _CTSEM_COLUMN_KERNEL = 8

@inline _ctsem_bkj(B, k, j, ::Val{false}) = @inbounds B[k, j]
@inline _ctsem_bkj(B, k, j, ::Val{true}) = @inbounds B[j, k]
@inline _ctsem_bcols(B, ::Val{false}) = size(B, 2)
@inline _ctsem_bcols(B, ::Val{true}) = size(B, 1)

"""
    _ctsem_colmul!(C, A, B, alpha, beta, Val(transB))

`C = alpha * A * op(B) + beta * C`, with `op(B)` = `B` or `B'`, column by
column: `C[:, j]` is a sum of columns of `A` weighted by `op(B)[:, j]`, four
columns a pass. Contiguous in the row index, so it vectorises.

Writes the leading `size(A, 1)` by `size(op(B), 2)` block of `C` and nothing
else, as the dot-product loops always did: callers hand over a buffer sized for
the largest case and use its leading block -- the measurement update, for the
variables a row actually observed. Iterating over `C`'s own columns instead
wrote past that block and, under `@inbounds`, read past the end of `B`: a
wrong likelihood with no error, on one study of a 13-study fit whose other
units agreed to 1e-6 (juliaFit d568f112).

`@noinline`: inlined, its four unrolled loops were compiled into every caller
instance -- each of
some 130 product sites, at Float64 and every dual type the routes use, against
the 428 filter instances the precompile workload leaves -- which took the
engine's precompile from 811 s to about 1400 s (local) and the CI runners past
their limits. The wrappers stay inline: their dot loops are small code, and a
call per tiny product cost 13% of a small model's Laplace gradient.
"""
@noinline function _ctsem_colmul!(C, A, B, alpha, beta, tb::Val)
    m = size(A, 1)
    K = size(A, 2)
    T = eltype(C)
    @inbounds for j in 1:_ctsem_bcols(B, tb)
        if iszero(beta)
            @simd for i in 1:m
                C[i, j] = zero(T)
            end
        elseif !isone(beta)
            @simd for i in 1:m
                C[i, j] *= beta
            end
        end
        k = 1
        while k + 3 <= K
            b1 = alpha * _ctsem_bkj(B, k, j, tb)
            b2 = alpha * _ctsem_bkj(B, k + 1, j, tb)
            b3 = alpha * _ctsem_bkj(B, k + 2, j, tb)
            b4 = alpha * _ctsem_bkj(B, k + 3, j, tb)
            @simd for i in 1:m
                C[i, j] = muladd(A[i, k], b1, muladd(A[i, k + 1], b2,
                    muladd(A[i, k + 2], b3, muladd(A[i, k + 3], b4, C[i, j]))))
            end
            k += 4
        end
        r = K - k + 1
        if r == 3
            b1 = alpha * _ctsem_bkj(B, k, j, tb)
            b2 = alpha * _ctsem_bkj(B, k + 1, j, tb)
            b3 = alpha * _ctsem_bkj(B, k + 2, j, tb)
            @simd for i in 1:m
                C[i, j] = muladd(A[i, k], b1, muladd(A[i, k + 1], b2,
                    muladd(A[i, k + 2], b3, C[i, j])))
            end
        elseif r == 2
            b1 = alpha * _ctsem_bkj(B, k, j, tb)
            b2 = alpha * _ctsem_bkj(B, k + 1, j, tb)
            @simd for i in 1:m
                C[i, j] = muladd(A[i, k], b1, muladd(A[i, k + 1], b2, C[i, j]))
            end
        elseif r == 1
            b1 = alpha * _ctsem_bkj(B, k, j, tb)
            @simd for i in 1:m
                C[i, j] = muladd(A[i, k], b1, C[i, j])
            end
        end
    end
    return C
end

"""`C = alpha * A * B + beta * C`."""
@inline function _ctsem_mul!(C, A, B, alpha=true, beta=false)
    _ctsem_small_product(C, size(A, 2)) ||
        return (_ctsem_barrier(mul!, C, A, B, alpha, beta); C)
    size(A, 2) >= _CTSEM_COLUMN_KERNEL &&
        return _ctsem_colmul!(C, A, B, alpha, beta, Val(false))
    @inbounds for j in axes(B, 2), i in axes(A, 1)
        acc = zero(eltype(C))
        for k in axes(A, 2)
            acc += A[i, k] * B[k, j]
        end
        C[i, j] = iszero(beta) ? alpha * acc : alpha * acc + beta * C[i, j]
    end
    return C
end

"""Task-local-storage key for `_ctsem_transpose_buffer`'s tables."""
struct _CTSEMTransposeKey{T} end

"""
    _ctsem_transpose_buffer(A)

A `size(A, 2) x size(A, 1)` matrix of `A`'s element type, kept in the calling
task's local storage and reused by every later call of that shape on the same
task. Task-local rather than per thread: a task can move between threads, and
the pool's workers are tasks.

The task's store is keyed by a singleton per element type and holds a typed
table keyed by shape. Keying the store by `(name, T, rows, cols)` directly
built that tuple on the heap every call -- it holds a type and two runtime
integers -- 48 bytes a call and three times as slow (83 ns against 30, dev1),
on a path that exists to keep the threads from waiting on the collector.
"""
function _ctsem_transpose_buffer(A::AbstractMatrix{T}) where {T}
    table = get!(Dict{Tuple{Int,Int},Matrix{T}}, _ctsem_scratch(),
        _CTSEMTransposeKey{T}())::Dict{Tuple{Int,Int},Matrix{T}}
    return get!(() -> Matrix{T}(undef, size(A, 2), size(A, 1)), table,
        (size(A, 2), size(A, 1)))
end

"""`C = alpha * A' * B + beta * C`.

From the column kernel's size up, `A'` is formed explicitly and the product is
the column kernel's: the copy is O(nk) against the O(nkm) product, into a
buffer the task keeps per element type and shape (`_ctsem_transpose_buffer`)
-- a fresh matrix every call, from `permutedims`, made eight threads stop for
the collector instead. Above `_CTSEM_SMALL_PRODUCT` the same copy goes to an
ordinary `gemm`, never to `dgemm_tn`, which locks even where OpenBLAS's other
products do not: a quarter of a 12-state gradient's samples at 8 threads on
dev1.
"""
@inline function _ctsem_mulTN!(C, A, B, alpha=true, beta=false)
    _ctsem_small_product(C, size(A, 1)) || return (_ctsem_barrier(mul!, C,
        transpose!(_ctsem_transpose_buffer(A), A), B, alpha, beta); C)
    size(A, 1) >= _CTSEM_COLUMN_KERNEL && return _ctsem_colmul!(C,
        transpose!(_ctsem_transpose_buffer(A), A), B, alpha, beta, Val(false))
    @inbounds for j in axes(B, 2), i in axes(A, 2)
        acc = zero(eltype(C))
        for k in axes(A, 1)
            acc += A[k, i] * B[k, j]
        end
        C[i, j] = iszero(beta) ? alpha * acc : alpha * acc + beta * C[i, j]
    end
    return C
end

"""`C = alpha * A * B' + beta * C`."""
@inline function _ctsem_mulNT!(C, A, B, alpha=true, beta=false)
    _ctsem_small_product(C, size(A, 2)) ||
        return (_ctsem_barrier(mul!, C, A, transpose(B), alpha, beta); C)
    size(A, 2) >= _CTSEM_COLUMN_KERNEL &&
        return _ctsem_colmul!(C, A, B, alpha, beta, Val(true))
    @inbounds for j in axes(B, 1), i in axes(A, 1)
        acc = zero(eltype(C))
        for k in axes(A, 2)
            acc += A[i, k] * B[j, k]
        end
        C[i, j] = iszero(beta) ? alpha * acc : alpha * acc + beta * C[i, j]
    end
    return C
end

"""`y = alpha * A * x + beta * y`."""
@inline function _ctsem_mulvec!(y, A, x, alpha=true, beta=false)
    _ctsem_small_matvec(y, size(A, 2)) ||
        return (_ctsem_barrier(mul!, y, A, x, alpha, beta); y)
    @inbounds for i in axes(A, 1)
        acc = zero(eltype(y))
        for k in axes(A, 2)
            acc += A[i, k] * x[k]
        end
        y[i] = iszero(beta) ? alpha * acc : alpha * acc + beta * y[i]
    end
    return y
end

"""`y = alpha * A' * x + beta * y`."""
@inline function _ctsem_mulTvec!(y, A, x, alpha=true, beta=false)
    _ctsem_small_matvec(y, size(A, 1)) ||
        return (_ctsem_barrier(mul!, y, transpose(A), x, alpha, beta); y)
    @inbounds for i in axes(A, 2)
        acc = zero(eltype(y))
        for k in axes(A, 1)
            acc += A[k, i] * x[k]
        end
        y[i] = iszero(beta) ? alpha * acc : alpha * acc + beta * y[i]
    end
    return y
end

"""`C = alpha * a * b' + beta * C`, the outer product of two vectors."""
@inline function _ctsem_outer!(C, a, b, alpha=true, beta=false)
    @inbounds for j in eachindex(b), i in eachindex(a)
        acc = a[i] * b[j]
        C[i, j] = iszero(beta) ? alpha * acc : alpha * acc + beta * C[i, j]
    end
    return C
end

"""
    _ctsem_symeig(A)

Eigenvalues (ascending) and orthonormal eigenvectors of a small symmetric
`Float64` matrix, by cyclic Jacobi rotations -- for the same reason as the
Cholesky above: LAPACK's `syevr` takes the process-global lock, and this is
called from inside the threaded unit loop.

Reads `A` as symmetric through both triangles averaged, and does not modify it.
Cyclic Jacobi converges quadratically once the off-diagonal is small, and is
accurate to a few ulps of the largest eigenvalue in *absolute* terms, which is
all its one caller (the eigenwise prior floor, where only eigenvalues near 1
matter) needs. Cost is about `6 n^3` flops per sweep and typically 5 to 8
sweeps, so it is for the dimension of one unit's random effects, not a model's.
"""
function _ctsem_symeig(A::AbstractMatrix{Float64}; maxsweeps::Integer=60)
    n = size(A, 1)
    a = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        a[i, j] = (A[i, j] + A[j, i]) / 2
    end
    V = Matrix{Float64}(LinearAlgebra.I, n, n)
    scale = 0.0
    @inbounds for j in 1:n, i in 1:n
        scale += a[i, j]^2
    end
    tol = (eps(Float64) * sqrt(scale))^2
    @inbounds for _ in 1:maxsweeps
        off = 0.0
        for q in 2:n, p in 1:(q - 1)
            off += a[p, q]^2
        end
        off <= tol && break
        for q in 2:n, p in 1:(q - 1)
            apq = a[p, q]
            iszero(apq) && continue
            theta = (a[q, q] - a[p, p]) / (2 * apq)
            t = (theta >= 0 ? 1.0 : -1.0) / (abs(theta) + sqrt(theta^2 + 1))
            c = 1 / sqrt(t^2 + 1)
            s = t * c
            for k in 1:n
                akp = a[k, p]; akq = a[k, q]
                a[k, p] = c * akp - s * akq
                a[k, q] = s * akp + c * akq
            end
            for k in 1:n
                apk = a[p, k]; aqk = a[q, k]
                a[p, k] = c * apk - s * aqk
                a[q, k] = s * apk + c * aqk
            end
            for k in 1:n
                vkp = V[k, p]; vkq = V[k, q]
                V[k, p] = c * vkp - s * vkq
                V[k, q] = s * vkp + c * vkq
            end
        end
    end
    values = [a[i, i] for i in 1:n]
    order = sortperm(values)
    return (values=values[order], vectors=V[:, order])
end
