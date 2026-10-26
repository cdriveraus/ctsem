using LinearAlgebra

################################################################################
# Lyapunov related functions
################################################################################
"""
    gauss_summation(n)

Return the triangular number `n * (n + 1) ÷ 2`.
"""
@inline gauss_summation(n::Int) = (n * (n + 1)) ÷ 2

"""
Diffusion-block size above which a `Float64` Lyapunov solve takes the Schur
route (`real_schur.jl`) rather than the packed `ksolve!`.

`ksolve!` solves a `k(k+1)/2`-square dense system, O(k^6) to factor and O(k^4)
to solve against its factors; the Schur route is O(k^3) for both. Both are
written out, with no LAPACK call and no lock, so the choice is arithmetic
alone. Microseconds a solve, single-threaded, minimum of ten batches, with a
fresh factorisation each call as a state-dependent drift pays at every row,
and against cached factors as a linear model's rows solve
(`scripts/stangrad-2026-10`'s bench_lyap.jl; local i9-12900KS / dev2 EPYC 7702):

    k                   6     8    10    11    12    16     20      30
    packed, fresh     1.9   8.1  25.7  42.9    70   358   1280   23211   local
    Schur, fresh      2.9   7.0  10.2  13.3  14.6  26.2   43.4     114
    packed, fresh     3.3  12.3  36.5  61.0    97   525   1931   38700   dev2
    Schur, fresh      5.4  12.2  18.5  24.3  27.4  42.7   76.7     210
    packed, cached    0.3   0.7   1.5   2.1   2.9   8.4   18.9    89.8   local
    Schur, cached     0.7   1.6   2.6   3.5   3.5   5.9    9.8    25.4
    packed, cached    0.6   1.6   3.4   4.8   6.1  18.6   41.9     190   dev2
    Schur, cached     1.5   2.8   4.5   5.8   5.4  10.3   17.3    48.8

Fresh, the Schur route wins from k = 8 on both machines; against cached
factors from 12-14, losing about a microsecond a solve below that. Ten sits
between them. Eight tasks each refactoring on their own buffer ran 4.9-5.1
times one task's throughput on the Schur route and 3.2-4.6 times on the packed
one (local), so neither waits on the other. A buffer's route is fixed when it
is made (`LyapBuffer(T, n)`), so `ctsem_set_lyapunov_schur_above!` reaches the
workspaces built after it.
"""
const _CTSEM_LYAP_SCHUR_ABOVE = Ref(10)

"""Set the diffusion-block size above which the Schur Lyapunov route is used."""
function ctsem_set_lyapunov_schur_above!(k::Integer)
    k >= 1 || throw(ArgumentError("threshold must be at least 1"))
    _CTSEM_LYAP_SCHUR_ABOVE[] = Int(k)
    return _CTSEM_LYAP_SCHUR_ABOVE[]
end

"""
    LyapBuffer{T}

Workspace for the Lyapunov solve `A X + X A' + Q = 0` on a `dim x dim` block,
by the packed `ksolve!` system or by the real Schur factors of `A`
(`real_schur.jl`), whichever `LyapBuffer(T, n)` chose: Schur above
`_CTSEM_LYAP_SCHUR_ABOVE` for `Float64` and `Float32`, the packed system
otherwise and for every `ForwardDiff.Dual`.

One type for both routes, with the route a run-time field, so that a
workspace's type does not depend on the model's size: two types here would
compile the filter again for every model above the threshold (CLAUDE.md, "A
model's types do not depend on the model"). The fields of the route not taken
are empty.

Either route keeps its factors across calls and rebuilds them only when `A`
changes (`last_A`; `_blocks_identical` compares Dual partials too, so a factor
never serves a point it was not built at). The reverse pass solves against one
`A` at every substep of a linear model; without the cache the packed route
made the 6-latent reverse pass 8% slower (dev1).
"""
mutable struct LyapBuffer{TYPE<:Number}
    dim::Int
    schur::Bool
    # The packed route: the system's LU, the right-hand side, the pivots.
    ksolve_O::Matrix{TYPE}
    ksolve_triQ::Vector{TYPE}
    ksolve_piv::Vector{Int}
    ksolve_dim::Int
    ksolve_system_dim::Int
    # The Schur route: `A = U S U'`, the right-hand side and then the solution
    # in Schur coordinates (`C`), product scratch, the reduction's reflector
    # (and a block solve's right-hand side) in `work`, a block's Kronecker
    # system in `small`, and whether the QR iteration converged.
    S::Matrix{TYPE}
    U::Matrix{TYPE}
    C::Matrix{TYPE}
    tmp::Matrix{TYPE}
    work::Vector{TYPE}
    small::Matrix{TYPE}
    converged::Bool
    last_A::Matrix{TYPE}
    valid::Bool
end

"""
    LyapBuffer(T, n)

A Lyapunov workspace for `n × n` matrices with scalar type `T`, on the Schur
route for `Float64` and `Float32` above `_CTSEM_LYAP_SCHUR_ABOVE` and on the
packed route otherwise.
"""
function LyapBuffer(::Type{TYPE}, n::Int) where {TYPE<:Number}
    schur = TYPE <: LinearAlgebra.BlasReal && n > _CTSEM_LYAP_SCHUR_ABOVE[]
    tri_n = schur ? 0 : gauss_summation(n)
    s = schur ? n : 0
    return LyapBuffer{TYPE}(n, schur,
        zeros(TYPE, tri_n, tri_n), zeros(TYPE, tri_n), zeros(Int, tri_n), n, tri_n,
        zeros(TYPE, s, s), zeros(TYPE, s, s), zeros(TYPE, s, s), zeros(TYPE, s, s),
        zeros(TYPE, schur ? max(2n, 4) : 0), zeros(TYPE, schur ? 4 : 0, schur ? 4 : 0),
        false, zeros(TYPE, n, n), false)
end

"""
    LyapBuffer(mat::AbstractMatrix)

Create a Lyapunov workspace sized and typed from `mat`.
"""
LyapBuffer(mat::AbstractMatrix) = LyapBuffer(eltype(mat), size(mat, 1))

"""
    my_lyap!(X, A, Q, buffer::LyapBuffer)

Solve `A * X + X * A' + Q = 0` on the leading `buffer.dim` block, by the
buffer's route. The solution is written to that block of `X`.
"""
function my_lyap!(X::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, Q::AbstractMatrix{TYPE},
    buffer::LyapBuffer{TYPE}) where {TYPE<:Number}
    buffer.schur && return _lyap_schur!(X, A, Q, buffer, false)
    d = buffer.ksolve_dim
    ntri = buffer.ksolve_system_dim
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
    my_lyap_transposed!(W, A, Q, buffer, At)

Solve `A' W + W A + Q = 0`, the equation a Lyapunov pullback solves. The Schur
route solves it on `A`'s own factors (`A' = U S' U'`), so the pullback
factorises only where `A` changes, as the forward pass does; the packed route
solves against the transpose, which it copies into the scratch `At`.
"""
function my_lyap_transposed!(W::AbstractMatrix{TYPE}, A::AbstractMatrix, Q::AbstractMatrix{TYPE},
    buffer::LyapBuffer{TYPE}, At::AbstractMatrix{TYPE}) where {TYPE<:Number}
    buffer.schur && return _lyap_schur!(W, A, Q, buffer, true)
    n = buffer.dim
    # A transposed copy rather than `adjoint(A)`: the packed route compares
    # and factors the matrix it is handed.
    @inbounds for j in 1:n, i in 1:n
        At[j, i] = A[i, j]
    end
    return my_lyap!(W, At, Q, buffer)
end

"""
    _lyap_schur!(X, A, Q, buffer, transposed)

The Schur route: factor `A` where it has changed, then `C = -U' Q U`, the
quasi triangular solve (`_ctsem_quasitri_lyap!`, `S Y + Y S' = C`, or
`S' Y + Y S = C` when `transposed`) and `X = U Y U'`, symmetrised. Every
product goes through the engine's kernels (`small_linalg.jl`). Where the QR
iteration did not converge, or `A` was not finite, `X` is NaN, as a
non-finite `A` leaves the packed route's.

Defined for `Float64` and `Float32` only, which are the only buffers
`LyapBuffer(T, n)` puts on this route. The route is chosen by a run-time
field, so without the fallback below every Dual type a model runs would
compile the whole route as well, although it never takes it.
"""
_lyap_schur!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix, buffer::LyapBuffer,
    transposed::Bool) =
    throw(ArgumentError("the Schur Lyapunov route is for Float64 and Float32 buffers"))

function _lyap_schur!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix,
    buffer::LyapBuffer{T}, transposed::Bool) where {T<:LinearAlgebra.BlasReal}
    n = buffer.dim
    S = buffer.S
    U = buffer.U
    C = buffer.C
    tmp = buffer.tmp
    if !(buffer.valid && _blocks_identical(buffer.last_A, A, n, n))
        _CTSEM_OPCOUNT.lyap_schur[] += 1
        _copy_block!(S, A, n, n)
        buffer.converged = _ctsem_real_schur!(S, U, buffer.work, n)
        _copy_block!(buffer.last_A, A, n, n)
        buffer.valid = true
    end
    if !buffer.converged
        @inbounds for j in 1:n, i in 1:n
            X[i, j] = T(NaN)
        end
        return X
    end
    _copy_block!(C, Q, n, n)
    _ctsem_mulTN!(tmp, U, C)
    _ctsem_mul!(C, tmp, U, -one(T), false)
    _ctsem_quasitri_lyap!(C, S, n, buffer.work, buffer.small, transposed)
    _ctsem_mul!(tmp, U, C)
    _ctsem_mulNT!(C, tmp, U)
    @inbounds for j in 1:n, i in 1:j
        v = (C[i, j] + C[j, i]) / 2
        X[i, j] = v
        X[j, i] = v
    end
    return X
end
