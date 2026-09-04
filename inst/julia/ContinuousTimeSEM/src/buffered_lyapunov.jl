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
    # The `A` whose packed system is currently factored in `ksolve_O`, so a
    # repeated solve against the same `A` reuses the factors -- the same
    # economy `LyapSchurBuffer` makes. The reverse pass solves against one
    # `A` at every substep of a linear model; without this, moving the Schur
    # threshold up handed those models a fresh LU per substep and made the
    # 6-latent reverse pass 8% slower (dev1).
    last_A::Matrix{TYPE}
    valid::Bool
    LyapKsolveBuffer{TYPE}(n::Int) where {TYPE <: Number} = begin
        tri_n = gauss_summation(n)
        new{TYPE,n,tri_n}(
            zeros(TYPE, tri_n, tri_n),
            zeros(TYPE, tri_n),
            zeros(Int, tri_n),
            Val(n),
            Val(tri_n),
            zeros(TYPE, n, n),
            false,
        )
    end
end

"""
Diffusion-block size above which a `Float64` Lyapunov solve goes to LAPACK's
Schur route rather than the packed `ksolve!`.

`ksolve!` is O(k^6) and Schur is O(k^3), so the crossover is real, but it
sits far higher than the operation counts suggest: the packed solve is a
plain LU with no library call, while `schur!` takes LAPACK's process-global
lock and allocates fresh work arrays every time. Measured on dev1 (23-core
EPYC, single thread, minimum of 2000 repeats, a fresh factorisation each
call as a state-dependent model pays at every substep):

    k          2     4     6     8    10    11    12
    ksolve µs  0.09  0.65  2.7   6.8  14.6  22.2  30.8
    schur  µs  1.2   3.9   6.6  11.6  17.5  21.6  26.1
    schur allocates 1-5 KB per call; ksolve none

The two agree to about 1e-14 relative throughout. Beyond that table the packed
solve loses fast when it has to refactor -- k = 16: 357 vs 56; k = 20: 1344 vs
80; k = 24: 5628 vs 124 -- but a solve against a *cached* factor costs 13, 33
and 69 microseconds at those sizes, and a linear model's reverse pass solves
against one `A` throughout.

The default is therefore *never Schur*: no LAPACK call, no lock, at any size.
The one case that pays is state-dependent drift with more than about twelve
diffusing states, which refactors at every substep; a session fitting such a
model on one core can call `ctsem_set_lyapunov_schur_above!(12)`.
"""
const _CTSEM_LYAP_SCHUR_ABOVE = Ref(typemax(Int))

"""Set the diffusion-block size above which the Schur Lyapunov route is used."""
function ctsem_set_lyapunov_schur_above!(k::Integer)
    k >= 1 || throw(ArgumentError("threshold must be at least 1"))
    _CTSEM_LYAP_SCHUR_ABOVE[] = Int(k)
    return _CTSEM_LYAP_SCHUR_ABOVE[]
end

"""
    LyapBuffer(T, n)

Create a Lyapunov workspace for `n × n` matrices with scalar type `T`.

`Float64` systems above `_CTSEM_LYAP_SCHUR_ABOVE` use a Schur buffer; other
cases, including every `ForwardDiff.Dual` system, use the direct packed
`ksolve!` buffer.
"""
function LyapBuffer(::Type{TYPE}, n::Int) where {TYPE<:Number}
    if TYPE <: LinearAlgebra.BlasFloat && n > _CTSEM_LYAP_SCHUR_ABOVE[]
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
    _CTSEM_OPCOUNT.lyap_schur[] += 1
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
function my_lyap!(X::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, Q::AbstractMatrix{TYPE},
    buffer::LyapKsolveBuffer{TYPE,d,ntri}) where {TYPE<:Number,d,ntri}
    # Small BLAS matrices and non-BLAS element types (every ForwardDiff.Dual
    # path) take the packed system. Its LU is kept across calls and rebuilt
    # only when `A` changes; `_blocks_identical` compares Dual partials too,
    # so a factor never serves a point it was not built at.
    O = buffer.ksolve_O
    if !(buffer.valid && _blocks_identical(buffer.last_A, A, d, d))
        _CTSEM_OPCOUNT.lyap_ksolve[] += 1
        _ksolve_system_matrix!(O, A, buffer.ksolve_dim)
        _lu_factor_generic!(O, buffer.ksolve_piv, ntri)
        _copy_block!(buffer.last_A, A, d, d)
        buffer.valid = true
    end
    triQ = buffer.ksolve_triQ
    _ksolve_pack_upper!(triQ, Q, buffer.ksolve_dim)
    _lu_solve_generic!(O, buffer.ksolve_piv, triQ, ntri)
    rmul!(triQ, -one(TYPE))
    _ksolve_unpack_upper!(X, triQ, buffer.ksolve_dim)
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
