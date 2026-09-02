using LinearAlgebra
using ForwardDiff

################################################################################
# Helper functions 
################################################################################
"""
    _custom_abs(x)

Return `abs(x)` for ordinary numeric values.

This method is paired with a `ForwardDiff.Dual` specialization so pivoting and
approximate comparisons can work on primal values during automatic
differentiation.
"""
@inline _custom_abs(x) = abs(x)

"""
    _finite_deep(x)

Whether `x` is finite *including its derivatives*.

`isfinite` on a `ForwardDiff.Dual` tests only the value, so a NaN in the
partials passes every validity guard in this package and surfaces much later as
a rejected trial point with a finite objective and an unusable gradient. This
walks the whole dual tree instead.
"""
_finite_deep(x::Real) = isfinite(x)
_finite_deep(x::ForwardDiff.Dual) = _finite_deep(ForwardDiff.value(x)) &&
    all(_finite_deep, ForwardDiff.partials(x))
_finite_deep(x::AbstractArray) = all(_finite_deep, x)

"""
    _custom_abs(x::ForwardDiff.Dual)

Return the absolute value of the primal part of a dual number.
"""
@inline _custom_abs(x::ForwardDiff.Dual) = abs(ForwardDiff.value(x))

@inline _val(::Val{N}) where {N} = N

"""
    _isapprox_default_rtol(T)

Return the default relative tolerance used by `_isapprox_matrix_noalloc`.

For non-floating scalar types this falls back to `sqrt(eps(Float64))`.
"""
@inline _isapprox_default_rtol(::Type) = sqrt(eps(Float64))

"""
    _isapprox_default_rtol(::Type{T}) where {T<:AbstractFloat}

Return `sqrt(eps(T))` for floating-point matrix element types.
"""
@inline _isapprox_default_rtol(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(T))

"""
    _isapprox_default_rtol(::Type{<:ForwardDiff.Dual})

Return a relative tolerance based on the dual number's floating primal type.
"""
@inline _isapprox_default_rtol(::Type{<:ForwardDiff.Dual{Tag,V,N}}) where {Tag,V<:AbstractFloat,N} = sqrt(eps(V))

"""
    add_diag!(mat, val::Number)

Add the scalar `val` to each diagonal entry of `mat` in place.
"""
function add_diag!(mat, val::Number)
    return add_diag!(mat, val, Val(min(size(mat, 1), size(mat, 2))))
end

function add_diag!(mat, val::Number, ::Val{d}) where {d}
    @boundscheck begin
        size(mat, 1) >= d || throw(DimensionMismatch("matrix row count must be at least d"))
        size(mat, 2) >= d || throw(DimensionMismatch("matrix column count must be at least d"))
    end
    @inbounds for idx in 1:d
        mat[idx, idx] += val
    end
    return mat
end

"""
    add_diag!(mat, vals::AbstractVector)

Add `vals[i]` to the `i`th diagonal entry of `mat` in place.
"""
function add_diag!(mat, vals::AbstractVector)
    return add_diag!(mat, vals, Val(length(vals)))
end

function add_diag!(mat, vals::AbstractVector, ::Val{d}) where {d}
    @boundscheck begin
        length(vals) == d || throw(DimensionMismatch("vals length must match d"))
        size(mat, 1) >= d || throw(DimensionMismatch("matrix row count must be at least d"))
        size(mat, 2) >= d || throw(DimensionMismatch("matrix column count must be at least d"))
    end
    @inbounds for idx in 1:d
        mat[idx, idx] += vals[idx]
    end
    return mat
end

"""
    _matvec_mul!(y, A, x::AbstractVector)

Compute `y .= A * x` without allocating temporaries.

This manual implementation is used for `SubArray`-backed EKF buffers where
`mul!` is not suitable.
"""
function _matvec_mul!(y::AbstractVector, A::AbstractMatrix, x::AbstractVector)
    return _matvec_mul!(y, A, x, Val(size(A, 1)), Val(size(A, 2)))
end

function _matvec_mul!(
    y::AbstractVector,
    A::AbstractMatrix,
    x::AbstractVector,
    ::Val{rows},
    ::Val{cols},
) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && length(y) == rows || throw(DimensionMismatch("A*y output size mismatch"))
        size(A, 2) == cols && length(x) == cols || throw(DimensionMismatch("A*x input size mismatch"))
    end
    @inbounds for row in 1:rows
        acc = zero(eltype(y))
        for col in 1:cols
            acc += A[row, col] * x[col]
        end
        y[row] = acc
    end
    return y
end

"""
    _matvec_mul!(y, A, x::AbstractMatrix)

Compute `y .= A * x[:, 1]` without allocating temporaries.

The matrix input must have exactly one column.
"""
function _matvec_mul!(y::AbstractVector, A::AbstractMatrix, x::AbstractMatrix)
    return _matvec_mul!(y, A, x, Val(size(A, 1)), Val(size(A, 2)))
end

function _matvec_mul!(
    y::AbstractVector,
    A::AbstractMatrix,
    x::AbstractMatrix,
    ::Val{rows},
    ::Val{cols},
) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && length(y) == rows || throw(DimensionMismatch("A*y output size mismatch"))
        size(x, 2) == 1 || throw(DimensionMismatch("x must be a column matrix"))
        size(A, 2) == cols && size(x, 1) == cols || throw(DimensionMismatch("A*x input size mismatch"))
    end
    @inbounds for row in 1:rows
        acc = zero(eltype(y))
        for col in 1:cols
            acc += A[row, col] * x[col, 1]
        end
        y[row] = acc
    end
    return y
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` for BLAS-compatible strided matrices.
"""
function _mul_right_transpose!(C::StridedMatrix{T}, A::StridedMatrix{T}, B::StridedMatrix{T}) where {T<:LinearAlgebra.BlasFloat}
    return _mul_right_transpose!(C, A, B, Val(size(C, 1)), Val(size(C, 2)), Val(size(A, 2)))
end

function _mul_right_transpose!(
    C::StridedMatrix{T},
    A::StridedMatrix{T},
    B::StridedMatrix{T},
    ::Val{rows},
    ::Val{cols},
    ::Val{inner},
) where {T<:LinearAlgebra.BlasFloat, rows, cols, inner}
    @boundscheck begin
        size(C, 1) == rows && size(A, 1) == rows || throw(DimensionMismatch("A row count must match C row count"))
        size(C, 2) == cols && size(B, 1) == cols || throw(DimensionMismatch("B row count must match C column count"))
        size(A, 2) == inner && size(B, 2) == inner || throw(DimensionMismatch("A and B must have the same column count"))
    end
    LinearAlgebra.BLAS.gemm!('N', 'T', one(T), A, B, zero(T), C)
    return C
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` with explicit loops for generic numeric matrices.
"""
function _mul_right_transpose!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T<:Number}
    return _mul_right_transpose!(C, A, B, Val(size(C, 1)), Val(size(C, 2)), Val(size(A, 2)))
end

function _mul_right_transpose!(
    C::AbstractMatrix{T},
    A::AbstractMatrix{T},
    B::AbstractMatrix{T},
    ::Val{rows},
    ::Val{cols},
    ::Val{inner},
) where {T<:Number, rows, cols, inner}
    @boundscheck begin
        size(C, 1) == rows && size(A, 1) == rows || throw(DimensionMismatch("A row count must match C row count"))
        size(C, 2) == cols && size(B, 1) == cols || throw(DimensionMismatch("B row count must match C column count"))
        size(A, 2) == inner && size(B, 2) == inner || throw(DimensionMismatch("A and B must have the same column count"))
    end
    @inbounds for j in 1:cols
        for i in 1:rows
            acc = zero(T)
            for k in 1:inner
                acc += A[i, k] * B[j, k]
            end
            C[i, j] = acc
        end
    end
    return C
end

"""
    _isapprox_matrix_noalloc(A, B; atol=0.0, rtol=...)

Return whether two matrices are approximately equal without allocating.

Dual-number entries are compared using their primal values through
`_custom_abs`.
"""
function _isapprox_matrix_noalloc(
    A::AbstractMatrix,
    B::AbstractMatrix;
    atol::Real=0.0,
    rtol::Real=_isapprox_default_rtol(promote_type(eltype(A), eltype(B))),
)
    size(A) == size(B) || return false
    return _isapprox_matrix_noalloc(A, B, Val(size(A, 1)), Val(size(A, 2)); atol=atol, rtol=rtol)
end

function _isapprox_matrix_noalloc(
    A::AbstractMatrix,
    B::AbstractMatrix,
    ::Val{rows},
    ::Val{cols};
    atol::Real=0.0,
    rtol::Real=_isapprox_default_rtol(promote_type(eltype(A), eltype(B))),
) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && size(B, 1) == rows || return false
        size(A, 2) == cols && size(B, 2) == cols || return false
    end
    @inbounds for j in 1:cols, i in 1:rows
        aij = A[i, j]
        bij = B[i, j]
        diff = _custom_abs(aij - bij)
        scale = max(_custom_abs(aij), _custom_abs(bij))
        if diff > atol + rtol * scale
            return false
        end
    end
    return true
end

function _matrix_equal_noalloc(A::AbstractMatrix, B::AbstractMatrix, ::Val{rows}, ::Val{cols}) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && size(B, 1) == rows || return false
        size(A, 2) == cols && size(B, 2) == cols || return false
    end
    @inbounds for j in 1:cols, i in 1:rows
        A[i, j] == B[i, j] || return false
    end
    return true
end

"""
    _can_reuse_same_exponential(A, B)

Return whether the exponential computed for `A` can be reused for `B`.

For dual-number matrices, exact equality is required so derivative information
is preserved.
"""
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:ForwardDiff.Dual}, B::AbstractMatrix{<:ForwardDiff.Dual}) = A == B
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:ForwardDiff.Dual}, B::AbstractMatrix{<:ForwardDiff.Dual}, dim::Val{d}) where {d} =
    _matrix_equal_noalloc(A, B, dim, dim)

"""
    _can_reuse_same_exponential(A, B)

Return whether two ordinary matrices are close enough to share an exponential.
"""
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix) = _isapprox_matrix_noalloc(A, B)
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix, dim::Val{d}) where {d} =
    _isapprox_matrix_noalloc(A, B, dim, dim)
