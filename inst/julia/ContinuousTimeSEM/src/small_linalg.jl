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
size is not what anyone wants. Sixteen is where the two measured even on this
machine; `ctsem_set_small_linalg!` moves it.
"""
const _CTSEM_SMALL_CHOLESKY = Ref(16)

"""
Arithmetic budget below which a product is done by hand rather than by `gemm`.

The same trade as the factorization threshold, in the units a matrix product is
naturally measured in: `rows * cols * inner`. `16^3 = 4096` is a cube of side
sixteen, so a product stays hand-written exactly while a Cholesky of the same
dimension would.
"""
const _CTSEM_SMALL_PRODUCT = Ref(4096)

export ctsem_set_small_linalg!
"""
    ctsem_set_small_linalg!(; dimension, product)

Move the thresholds separating the engine's own small-matrix kernels from
LAPACK and BLAS. Exposed because the crossover is a property of the machine, not
of the mathematics -- the same reason `ctsem_set_block_threshold!` exists.
"""
function ctsem_set_small_linalg!(; dimension::Integer=_CTSEM_SMALL_CHOLESKY[],
    product::Integer=_CTSEM_SMALL_PRODUCT[])
    dimension >= 0 && product >= 0 ||
        throw(ArgumentError("thresholds must be non-negative"))
    _CTSEM_SMALL_CHOLESKY[] = Int(dimension)
    _CTSEM_SMALL_PRODUCT[] = Int(product)
    return (dimension=Int(dimension), product=Int(product))
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
# `mul!` on plain matrices scales across threads and is fine. `mul!` with a
# *transposed view* operand is not: measured at 23 threads on one BLAS thread,
#
#     operands                n=1     n=2     n=4
#     plain Matrix           4.49x   7.15x   4.71x
#     views of a buffer      4.39x   3.51x   6.29x
#     transpose(view)        0.09x   4.37x   0.11x
#
# and it is slower serially too, because Julia cannot hand a transposed
# non-contiguous view to `gemm` and falls back to a copy or to the generic
# kernel. The buffered reverse pass is made of exactly that shape -- every
# temporary is a view, and half the products transpose one side.
#
# At these sizes the answer is not to shuffle operands into a form BLAS likes.
# It is to not call BLAS: a hand-written multiply is about five times faster
# than `gemm` at n = 1 and comparable at n = 4, and it threads because it
# touches nothing outside its own arguments.
#
# All six follow `mul!(C, A, B, alpha, beta)`: `C = alpha * op(A) * op(B) +
# beta * C`, with `beta = 0` overwriting rather than reading `C` -- so an
# uninitialised buffer is safe.

"""
Is this product small enough to be worth doing by hand?

`inner` is the contracted dimension, which differs between the transposed forms
-- hence passing it rather than reading it off one operand.
"""
@inline _ctsem_small_product(C, inner::Integer) =
    length(C) * inner <= _CTSEM_SMALL_PRODUCT[]

"""`C = alpha * A * B + beta * C`."""
@inline function _ctsem_mul!(C, A, B, alpha=true, beta=false)
    _ctsem_small_product(C, size(A, 2)) || return mul!(C, A, B, alpha, beta)
    @inbounds for j in axes(B, 2), i in axes(A, 1)
        acc = zero(eltype(C))
        for k in axes(A, 2)
            acc += A[i, k] * B[k, j]
        end
        C[i, j] = iszero(beta) ? alpha * acc : alpha * acc + beta * C[i, j]
    end
    return C
end

"""`C = alpha * A' * B + beta * C`."""
@inline function _ctsem_mulTN!(C, A, B, alpha=true, beta=false)
    _ctsem_small_product(C, size(A, 1)) || return mul!(C, transpose(A), B, alpha, beta)
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
    _ctsem_small_product(C, size(A, 2)) || return mul!(C, A, transpose(B), alpha, beta)
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
    _ctsem_small_product(y, size(A, 2)) || return mul!(y, A, x, alpha, beta)
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
    _ctsem_small_product(y, size(A, 1)) || return mul!(y, transpose(A), x, alpha, beta)
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
