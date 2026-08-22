using LinearAlgebra

################################################################################
# Lyapunov related functions
################################################################################
"""
    gauss_summation(n)

Return the triangular number `n * (n + 1) ÷ 2`.
"""
@inline gauss_summation(n::Int) = (n * (n + 1)) ÷ 2

abstract type AbstractLyapBuffer{TYPE<:Number} end

"""
    LyapSchurBuffer{T}(n)

Workspace for Schur-based Lyapunov solves with BLAS-compatible scalar types.
"""
mutable struct LyapSchurBuffer{TYPE<:LinearAlgebra.BlasFloat} <: AbstractLyapBuffer{TYPE}
    S::Matrix{TYPE}
    U::Matrix{TYPE}
    tmp::Matrix{TYPE}
    # Last matrix factorized, so a repeated solve against the same `A` can
    # reuse the factors. See `_lyap_factorize_schur!`.
    last_A::Matrix{TYPE}
    valid::Bool
    LyapSchurBuffer{TYPE}(n::Int) where {TYPE<:LinearAlgebra.BlasFloat} = begin
        new(zeros(TYPE, n, n), zeros(TYPE, n, n), zeros(TYPE, n, n),
            zeros(TYPE, n, n), false)
    end
end

"""
    LyapKsolveBuffer{T}(n)

Workspace for direct packed-system Lyapunov solves.

This buffer is used for small systems and scalar types such as
`ForwardDiff.Dual` that are not supported by LAPACK.
"""
mutable struct LyapKsolveBuffer{TYPE<:Number,D,NTRI} <: AbstractLyapBuffer{TYPE}
    ksolve_O::Matrix{TYPE}
    ksolve_triQ::Vector{TYPE}
    ksolve_piv::Vector{Int}
    ksolve_dim::Val{D}
    ksolve_system_dim::Val{NTRI}
    LyapKsolveBuffer{TYPE}(n::Int) where {TYPE <: Number} = begin
        tri_n = gauss_summation(n)
        new{TYPE,n,tri_n}(
            zeros(TYPE, tri_n, tri_n),
            zeros(TYPE, tri_n),
            zeros(Int, tri_n),
            Val(n),
            Val(tri_n),
        )
    end
end

"""
    LyapBuffer(T, n)

Create a Lyapunov workspace for `n × n` matrices with scalar type `T`.

Large BLAS-compatible systems use a Schur buffer; other cases use the direct
packed `ksolve!` buffer.
"""
function LyapBuffer(::Type{TYPE}, n::Int) where {TYPE<:Number}
    if TYPE <: LinearAlgebra.BlasFloat && n > 4
        return LyapSchurBuffer{TYPE}(n)
    else
        return LyapKsolveBuffer{TYPE}(n)
    end
end

"""
    LyapBuffer(mat::AbstractMatrix)

Create a Lyapunov workspace sized and typed from `mat`.
"""
LyapBuffer(mat::AbstractMatrix) = LyapBuffer(eltype(mat), size(mat, 1))

"""
    _lyap_factorize_schur!(buffer, A)

Compute the Schur factorization of `A` into reusable Lyapunov workspace.

The input matrix `A` is preserved.
"""
function _lyap_factorize_schur!(buffer::LyapSchurBuffer{TYPE}, A::AbstractMatrix{TYPE}) where {TYPE<:LinearAlgebra.BlasFloat}
    # Reuse the factors when `A` has not changed. `schur!` allocates fresh
    # work arrays and a `Schur` object on every call -- it was the single
    # largest allocation source in the filter -- and the reverse pass in
    # particular solves against the *same* `A` at every row while only the
    # right-hand side varies. `_lyap_solve_factorized!` treats `S` and `U` as
    # read-only (LAPACK's `trsyl!` does not modify its `A`/`B` arguments), so
    # the cached factors stay valid across solves.
    n = size(A, 1)
    if buffer.valid && size(buffer.last_A, 1) == n && _blocks_identical(buffer.last_A, A, n, n)
        return nothing
    end

    # Preserve A by factorizing a workspace copy. `schur!` factorizes in
    # place and returns its `T` factor as the very array it was handed, so
    # `buffer.S` aliases itself across calls and needs no fresh storage.
    copyto!(buffer.S, A)
    F = schur!(buffer.S)

    # Keep direct references to Schur outputs (no elementwise copy).
    buffer.S = F.T
    buffer.U = F.Z
    _copy_block!(buffer.last_A, A, n, n)
    buffer.valid = true
    return nothing
end

"""
    my_lyap!(X, A, Q, buffer::LyapSchurBuffer)

Solve `A * X + X * A' + Q = 0` using a Schur factorization.

The solution is written to `X`.
"""
function my_lyap!(X::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, Q::AbstractMatrix{TYPE}, buffer::LyapSchurBuffer{TYPE}) where {TYPE<:LinearAlgebra.BlasFloat}
    # For larger BLAS-compatible matrices, use the Schur factorization approach, which is more efficient and numerically stable than the direct ksolve approach.
    _lyap_factorize_schur!(buffer, A)
    _lyap_solve_factorized!(X, buffer.S, buffer.U, Q, buffer.tmp)
    return X
end

"""
    my_lyap!(X, A, Q, buffer::LyapKsolveBuffer)

Solve `A * X + X * A' + Q = 0` using the direct packed `ksolve!` path.

The solution is written to `X`.
"""
function my_lyap!(X::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, Q::AbstractMatrix{TYPE}, buffer::LyapKsolveBuffer{TYPE}) where {TYPE<:Number}
    # Small BLAS matrices and non-BLAS element types use ksolve!.
    # This includes Dual-number paths (e.g., ForwardDiff).
    ksolve!(
        X,
        A,
        Q,
        buffer.ksolve_O,
        buffer.ksolve_triQ,
        buffer.ksolve_piv,
        buffer.ksolve_dim,
        buffer.ksolve_system_dim,
    )
    return X
end

"""
    _lyap_solve_factorized!(X, S, U, Q, tmp)

Solve a Lyapunov equation from a precomputed real Schur factorization.

`S` and `U` are the Schur factors of `A`; `tmp` is scratch storage. The solution
is written to `X`.
"""
function _lyap_solve_factorized!(X::AbstractMatrix, S::AbstractMatrix, U::AbstractMatrix, Q::AbstractMatrix, tmp::AbstractMatrix)
    # Transform RHS into Schur coordinates: X = -U' * Q * U
    mul!(tmp, adjoint(U), Q)
    mul!(X, tmp, U)
    rmul!(X, -one(eltype(X)))

    # Solve Schur-space Sylvester equation in-place: S*Y + Y*S' = X
    _, scale = LinearAlgebra.LAPACK.trsyl!('N', 'T', S, S, X)

    # This is a scaling factor to prevent overflow/underflow in the solution of the Sylvester equation. If the solution is scaled, we need to back-transform it by multiplying with the inverse of the scale factor.
    if scale != one(scale)
        rmul!(X, inv(scale))
    end

    # Back-transform
    mul!(tmp, U, X)
    mul!(X, tmp, adjoint(U))
    return X
end
