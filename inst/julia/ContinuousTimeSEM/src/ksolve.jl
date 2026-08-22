using LinearAlgebra

################################################################################
# ksolve
################################################################################

"""
    _ksolve_system_matrix!(O, A)

Fill `O` with the packed linear system for the continuous Lyapunov equation.

The rows and columns of `O` correspond to upper-triangular entries of a
symmetric unknown `X` in the equation `A * X + X * A'`.
"""
function _ksolve_system_matrix!(O::AbstractMatrix{T}, A::AbstractMatrix{T}, ::Val{d}) where {T<:Number, d}
    # d = size(A, 1)
    fill!(O, zero(T))

    if d == 1
        O[1, 1] = 2 * A[1, 1]
        return O
    end

    # For symmetric X packed by upper triangle, row (i,j) of O corresponds to:
    # (A*X + X*A')_ij = sum_k A[i,k] * X[k,j] + sum_k X[i,k] * A[j,k]
    @inbounds for j in 1:d
        base_j = ((j - 1) * j) ÷ 2

        for k in 1:d
            ajk = A[j, k]
            tri_k = ((k - 1) * k) ÷ 2
            y1 = k <= j ? base_j + k : tri_k + j
            m = min(j, k)

            # i <= k  => packed(i,k) = tri_k + i
            for i in 1:m
                z = base_j + i
                O[z, y1] += A[i, k]
                O[z, tri_k + i] += ajk
            end

            # i > k   => packed(k,i) = base_i + k
            base_i = (m * (m + 1)) ÷ 2
            for i in (m + 1):j
                z = base_j + i
                O[z, y1] += A[i, k]
                O[z, base_i + k] += ajk
                base_i += i
            end
        end
    end
    return O
end

_ksolve_system_matrix!(O::AbstractMatrix{T}, A::AbstractMatrix{T}) where {T<:Number} =
    _ksolve_system_matrix!(O, A, Val(size(A, 1)))

"""
    _ksolve_pack_upper!(triQ, Q)

Pack the upper triangle of the symmetric matrix `Q` into `triQ`.

The packing order matches the system matrix produced by
`_ksolve_system_matrix!`.
"""
function _ksolve_pack_upper!(triQ::AbstractVector{T}, Q::AbstractMatrix{T}, ::Val{d}) where {T<:Number, d}
    # d = size(Q, 1)
    z = 0
    @inbounds for j in 1:d
        for i in 1:j
            z += 1
            triQ[z] = Q[i, j]
        end
    end
    return triQ
end

_ksolve_pack_upper!(triQ::AbstractVector{T}, Q::AbstractMatrix{T}) where {T<:Number} =
    _ksolve_pack_upper!(triQ, Q, Val(size(Q, 1)))

"""
    _ksolve_unpack_upper!(AQ, triQ)

Unpack an upper-triangular packed vector into the symmetric matrix `AQ`.
"""
function _ksolve_unpack_upper!(AQ::AbstractMatrix{T}, triQ::AbstractVector{T}, ::Val{d}) where {T<:Number, d}
    # d = size(AQ, 1)
    z = 0
    @inbounds for j in 1:d
        for i in 1:j
            z += 1
            AQ[i, j] = triQ[z]
            AQ[j, i] = triQ[z]
        end
    end
    return AQ
end

_ksolve_unpack_upper!(AQ::AbstractMatrix{T}, triQ::AbstractVector{T}) where {T<:Number} =
    _ksolve_unpack_upper!(AQ, triQ, Val(size(AQ, 1)))

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

function _solve_square_system!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T}, piv::AbstractVector{Int}, ::Val) where {T<:LinearAlgebra.BlasFloat}
    return _solve_square_system!(A, B, piv)
end

"""
    _solve_square_system!(A, B::AbstractMatrix, piv)

Solve `A * X = B` in place for `ForwardDiff.Dual` matrix right-hand sides.

This custom LU path avoids BLAS/LAPACK calls that do not support dual numbers.
"""
function _solve_square_system!(A::AbstractMatrix{T}, B::AbstractMatrix{T}, piv::AbstractVector{Int}, ::Val{n}) where {T<:ForwardDiff.Dual, n}
    # n = size(A, 1)
    nrhs = size(B, 2)
    @boundscheck begin
        size(A, 2) == n || throw(DimensionMismatch("A must be square"))
        size(B, 1) == n || throw(DimensionMismatch("B row count must match A size"))
        length(piv) >= n || throw(DimensionMismatch("Pivot buffer too small"))
    end

    @inbounds for k in 1:n
        _, rel = findmax(_custom_abs, @view A[k:n, k])  # includes k, never empty
        pivot = k + rel - 1
        piv[k] = pivot


        # Row swap in A and B (manual loops avoid temporary index vectors).
        if pivot != k
            for col in k:n
                A[k, col], A[pivot, col] = A[pivot, col], A[k, col]
            end
            for col in 1:nrhs
                B[k, col], B[pivot, col] = B[pivot, col], B[k, col]
            end
        end

        # Elimination step below pivot.
        akk = A[k, k]
        for i in (k + 1):n
            lik = A[i, k] / akk
            A[i, k] = lik
            for col in (k + 1):n
                A[i, col] -= lik * A[k, col]
            end
            for col in 1:nrhs
                B[i, col] -= lik * B[k, col]
            end
        end
    end

    # Back substitution.
    @inbounds for i in n:-1:1
        aii = A[i, i]
        for rhs in 1:nrhs
            acc = B[i, rhs]
            for col in (i + 1):n
                acc -= A[i, col] * B[col, rhs]
            end
            B[i, rhs] = acc / aii
        end
    end
    return B
end

"""
    _solve_square_system!(A, b::AbstractVector, piv)

Solve `A * x = b` in place for `ForwardDiff.Dual` vector right-hand sides.
"""
function _solve_square_system!(A::AbstractMatrix{T}, B::AbstractVector{T}, piv::AbstractVector{Int}, ::Val{n}) where {T<:ForwardDiff.Dual, n}
    # n = size(A, 1)
    @boundscheck begin
        size(A, 2) == n || throw(DimensionMismatch("A must be square"))
        length(B) == n || throw(DimensionMismatch("B length must match A size"))
        length(piv) >= n || throw(DimensionMismatch("Pivot buffer too small"))
    end

    # Custom LU factorization with partial pivoting, adapted from the standard algorithm but using _custom_abs for comparisons to handle Dual numbers.
    @inbounds for k in 1:n
        _, rel = findmax(_custom_abs, @view A[k:n, k])  # includes k, never empty
        pivot = k + rel - 1
        piv[k] = pivot

        if pivot != k
            for col in k:n
                A[k, col], A[pivot, col] = A[pivot, col], A[k, col]
            end
            B[k], B[pivot] = B[pivot], B[k]
        end

        akk = A[k, k]
        for i in (k + 1):n
            lik = A[i, k] / akk
            A[i, k] = lik
            for col in (k + 1):n
                A[i, col] -= lik * A[k, col]
            end
            B[i] -= lik * B[k]
        end
    end

    @inbounds for i in n:-1:1
        acc = B[i]
        for col in (i + 1):n
            acc -= A[i, col] * B[col]
        end
        B[i] = acc / A[i, i]
    end
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

function _solve_square_system!(A::AbstractMatrix, B::AbstractVecOrMat, piv::AbstractVector{Int}, ::Val)
    return _solve_square_system!(A, B, piv)
end

"""
    ksolve!(AQ, A, Q, O, triQ, piv)

Solve the continuous Lyapunov equation `A * X + X * A' + Q = 0`.

The result is written to `AQ`. The arguments `O`, `triQ`, and `piv` are reusable
workspace buffers.
"""
function ksolve!(
    AQ::AbstractMatrix{T},
    A::AbstractMatrix{T},
    Q::AbstractMatrix{T},
    O::AbstractMatrix{T},
    triQ::AbstractVector{T},
    piv::AbstractVector{Int},
    dim::Val{d},
    system_dim::Val{ntri},
) where {T<:Number, d, ntri}
    @boundscheck begin
        # d = size(A, 1)
        size(A, 1) == d && size(A, 2) == d || throw(DimensionMismatch("A must be square"))
        size(Q, 1) == d && size(Q, 2) == d || throw(DimensionMismatch("Q must match size(A)"))
        size(AQ, 1) == d && size(AQ, 2) == d || throw(DimensionMismatch("AQ must match size(A)"))

        ntri == (d * (d + 1)) ÷ 2 || throw(DimensionMismatch("system_dim must match the packed size for A"))
        size(O, 1) == ntri && size(O, 2) == ntri || throw(DimensionMismatch("O must be ntri×ntri"))
        length(triQ) == ntri || throw(DimensionMismatch("triQ length must be ntri"))
        length(piv) >= ntri || throw(DimensionMismatch("Pivot buffer too small"))
    end

    _ksolve_system_matrix!(O, A, dim)
    _ksolve_pack_upper!(triQ, Q, dim)
    _solve_square_system!(O, triQ, piv, system_dim)
    rmul!(triQ, -one(T))
    _ksolve_unpack_upper!(AQ, triQ, dim)
    return AQ
end

function ksolve!(
    AQ::AbstractMatrix{T},
    A::AbstractMatrix{T},
    Q::AbstractMatrix{T},
    O::AbstractMatrix{T},
    triQ::AbstractVector{T},
    piv::AbstractVector{Int},
    dim::Val{d},
) where {T<:Number, d}
    return ksolve!(AQ, A, Q, O, triQ, piv, dim, Val((d * (d + 1)) ÷ 2))
end

function ksolve!(
    AQ::AbstractMatrix{T},
    A::AbstractMatrix{T},
    Q::AbstractMatrix{T},
    O::AbstractMatrix{T},
    triQ::AbstractVector{T},
    piv::AbstractVector{Int},
) where {T<:Number}
    return ksolve!(AQ, A, Q, O, triQ, piv, Val(size(A, 1)))
end

"""
    ksolve!(AQ, A, Q, O, triQ)

Solve the continuous Lyapunov equation using internally allocated pivot storage.
"""
function ksolve!(AQ::AbstractMatrix{T}, A::AbstractMatrix{T}, Q::AbstractMatrix{T}, O::AbstractMatrix{T}, triQ::AbstractVector{T}) where {T<:Number}
    piv = Vector{Int}(undef, length(triQ))
    return ksolve!(AQ, A, Q, O, triQ, piv, Val(size(A, 1)))
end
