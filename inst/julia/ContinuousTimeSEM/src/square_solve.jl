using LinearAlgebra

################################################################################
# Square linear solves without LAPACK
################################################################################

"""
    _solve_square_system!(A, B, piv)

Solve `A * X = B` in place for BLAS-compatible scalar types.

`A` is overwritten by its LU factorization, `B` is overwritten by the solution,
and `piv` supplies reusable pivot storage.
"""
function _solve_square_system!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T}, piv::AbstractVector{Int}) where {T<:LinearAlgebra.BlasFloat}
    # Uses LAPACK's LU factorization and solve routines, which are optimized for BLAS-compatible types.
    LinearAlgebra.LAPACK.getrf!(A, piv; check=false)
    LinearAlgebra.LAPACK.getrs!('N', A, piv, B)
    return B
end

"""
Size-dispatched form. Below `_CTSEM_SMALL_CHOLESKY` the generic LU above is
used even for BLAS element types, and that is not a compromise: `getrf!` on a
matrix this size is almost entirely the process-global lock OpenBLAS takes to
acquire its scratch buffer, so the plain loop is both faster on one core and
the difference between a subject loop that threads and one that does not. See
`small_linalg.jl`, which makes the same trade for the Cholesky.
"""
function _solve_square_system!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T},
    piv::AbstractVector{Int}, size_hint::Int) where {T<:LinearAlgebra.BlasFloat}
    n = size_hint
    n <= _CTSEM_SMALL_CHOLESKY[] &&
        return _solve_square_system_generic!(A, B, piv, size_hint)
    return _solve_square_system!(A, B, piv)
end

"""Every other scalar type -- `ForwardDiff.Dual` above all -- has no LAPACK
route at any size, so it goes straight to the generic loops."""
_solve_square_system!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T},
    piv::AbstractVector{Int}, size_hint::Int) where {T} =
    _solve_square_system_generic!(A, B, piv, size_hint)

"""
    _solve_square_system_generic!(A, B, piv, n)

Solve `A X = B` in place, `B` a vector or a matrix, by the engine's own LU
with partial pivoting (`_lu_factor_generic!`, then `_lu_solve_generic!`): the
route for `ForwardDiff.Dual` and the other types LAPACK cannot take, and for
BLAS types below `_CTSEM_SMALL_CHOLESKY`. `A` is left holding the factors as
`_lu_factor_generic!` lays them out. A caller solving against one matrix
several times factors it once and calls `_lu_solve_generic!` per right-hand
side instead.
"""
function _solve_square_system_generic!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T},
    piv::AbstractVector{Int}, n::Int) where {T}
    @boundscheck begin
        size(A, 2) == n || throw(DimensionMismatch("A must be square"))
        size(B, 1) == n || throw(DimensionMismatch("B row count must match A size"))
        length(piv) >= n || throw(DimensionMismatch("Pivot buffer too small"))
    end
    _lu_factor_generic!(A, piv, n)
    _lu_solve_generic!(A, piv, B, n)
    return B
end

"""
    _solve_square_system!(A, B, piv)

Solve `A * X = B` in place for generic scalar types using Julia's LU solver.

The `piv` argument is accepted for API compatibility with the specialized
methods.
"""
function _solve_square_system!(A::AbstractMatrix, B::AbstractVecOrMat, piv::AbstractVector{Int})
    ldiv!(lu!(A), B)
    return B
end

function _solve_square_system!(A::AbstractMatrix, B::AbstractVecOrMat, piv::AbstractVector{Int}, ::Int)
    return _solve_square_system!(A, B, piv)
end

"""
    _lu_factor_generic!(A, piv, n)

LU-factorize the leading `n × n` block of `A` in place with partial pivoting,
recording the row interchanges in `piv`, without a right-hand side.

Step `k` swaps only columns `k:n`, so a multiplier stays where the row it
belongs to stood when it was formed, and `_lu_solve_generic!` replays each
interchange just before eliminating its column. Every right-hand-side entry
then sees the same operations in the same order as an elimination carrying
the right-hand side along, so factoring apart costs nothing in the result
(bitwise, for Float64 and Duals) and the solve takes every right-hand side in
one pass. Unlike LAPACK's `getrf!`, the multipliers are not in permuted
order.
"""
function _lu_factor_generic!(A::AbstractMatrix{T}, piv::AbstractVector{Int}, n::Int) where {T}
    @inbounds for k in 1:n
        _, rel = findmax(_custom_abs, @view A[k:n, k])
        pivot = k + rel - 1
        piv[k] = pivot
        if pivot != k
            for col in k:n
                A[k, col], A[pivot, col] = A[pivot, col], A[k, col]
            end
        end
        akk = A[k, k]
        for i in (k + 1):n
            lik = A[i, k] / akk
            A[i, k] = lik
            for col in (k + 1):n
                A[i, col] -= lik * A[k, col]
            end
        end
    end
    return A
end

"""
    _lu_solve_generic!(LU, piv, B, n)

Solve `A X = B` in place, `B` a vector or a matrix, given the factors
`_lu_factor_generic!` left in `LU` and `piv`.
"""
function _lu_solve_generic!(LU::AbstractMatrix{T}, piv::AbstractVector{Int},
    B::AbstractVecOrMat{T}, n::Int) where {T}
    nrhs = size(B, 2)
    @inbounds for k in 1:n
        p = piv[k]
        if p != k
            for r in 1:nrhs
                B[k, r], B[p, r] = B[p, r], B[k, r]
            end
        end
        for i in (k + 1):n
            lik = LU[i, k]
            for r in 1:nrhs
                B[i, r] -= lik * B[k, r]
            end
        end
    end
    @inbounds for i in n:-1:1
        aii = LU[i, i]
        for r in 1:nrhs
            acc = B[i, r]
            for col in (i + 1):n
                acc -= LU[i, col] * B[col, r]
            end
            B[i, r] = acc / aii
        end
    end
    return B
end
