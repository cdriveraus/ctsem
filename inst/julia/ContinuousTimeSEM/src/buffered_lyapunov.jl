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
factors from 12-14, losing a microsecond or two a solve below that. The costs
are not symmetric -- a fresh packed factorisation grows as k^6, a cached
solve's loss stays a few microseconds -- so the line is drawn where fresh
factorisations win: Schur above eight. Eight tasks each refactoring on their
own buffer ran 4.9-5.1 times one task's throughput on the Schur route and
3.2-4.6 times on the packed one (local), so neither waits on the other. A
buffer's route is fixed when it is made (`LyapBuffer(T, n)`), so
`ctsem_set_lyapunov_schur_above!` reaches the workspaces built after it.
"""
const _CTSEM_LYAP_SCHUR_ABOVE = Ref(8)

"""
Diffusion-block size above which a `ForwardDiff.Dual` Lyapunov solve takes
the tangent route (`LyapTangents`): its values solved on whichever route their
own type takes at that size, then one solve per partial on the same factors,
where the packed route pushes the Duals through the LU. Microseconds a solve
at dual width 4 (a Laplace unit's curvature) and 8, local i9-12900KS,
`bench_lyap_dual2.jl`:

    k                         3      5      8     10     12
    w4 packed Dual, fresh   0.7    7.2     31    258    707
    w4 tangents, fresh      1.0    3.5   13.5     37     54
    w4 packed Dual, cached  0.1    1.5    2.3     20     39
    w4 tangents, cached     0.8    2.3    7.9     17     32
    w8 packed Dual, fresh   0.8    7.9     90    300    808
    w8 tangents, fresh      1.7    5.5     27     58     84
    w8 packed Dual, cached  0.3    1.6    9.5     22     44
    w8 tangents, cached     1.5    4.4   15.4     34     58

Tangents win fresh from five states, cached from ten to twelve; by the same
asymmetry as the values' line, above four.
"""
const _CTSEM_LYAP_TANGENTS_ABOVE = Ref(4)

"""
    ctsem_set_lyapunov_schur_above!(k; tangents = k)

Set the diffusion-block sizes above which a Lyapunov solve takes the Schur
route (`Float64`, `Float32`) and above which a Dual's takes tangents.
"""
function ctsem_set_lyapunov_schur_above!(k::Integer; tangents::Integer=k)
    (k >= 1 && tangents >= 1) || throw(ArgumentError("thresholds must be at least 1"))
    _CTSEM_LYAP_SCHUR_ABOVE[] = Int(k)
    _CTSEM_LYAP_TANGENTS_ABOVE[] = Int(tangents)
    return _CTSEM_LYAP_SCHUR_ABOVE[]
end

"""
    LyapTangents{V, B}

The Schur route at a `ForwardDiff.Dual` element type: the value parts of `A`
and `Q` (`A`, `Q`), solved on a buffer of the value type (`primal`, itself on
the Schur route; a nested Dual recurses), then each partial as the same
equation on the same factors -- differentiating `A X + X A' + Q = 0`,

    A dX + dX A' + (dA X + X dA' + dQ) = 0,

so `dX` is a Lyapunov solve with `Q` replaced by `dA X + X dA' + dQ`. Exact,
and a solve against cached factors per partial, where pushing the Duals
through a factorisation costs it once per partial and, on the packed route,
O(k^6) each time. `X`, `R` (that right-hand side), `dA` and `tmp` are scratch;
`parts` holds the partials' solutions.
"""
struct LyapTangents{V,B}
    primal::B
    A::Matrix{V}
    Q::Matrix{V}
    X::Matrix{V}
    R::Matrix{V}
    dA::Matrix{V}
    tmp::Matrix{V}
    parts::Vector{Matrix{V}}
end

"""
    LyapBuffer{T}

Workspace for the Lyapunov solve `A X + X A' + Q = 0` on a `dim x dim` block,
by the packed `ksolve!` system or by the real Schur factors of `A`
(`real_schur.jl`), whichever `LyapBuffer(T, n)` chose: Schur above
`_CTSEM_LYAP_SCHUR_ABOVE` for `Float64`, `Float32` and every
`ForwardDiff.Dual` over them (`LyapTangents`), the packed system otherwise.

One type per element type for both routes, with the route a run-time field,
so that a workspace's type does not depend on the model's size: two types
here would compile the filter again for every model above the threshold
(CLAUDE.md, "A model's types do not depend on the model"). The fields of the
route not taken are empty; `tangents` is `nothing` unless the element type is
a Dual.

Either route keeps its factors across calls and rebuilds them only when `A`
changes (`last_A`; `_blocks_identical` compares Dual partials too, so a factor
never serves a point it was not built at). The reverse pass solves against one
`A` at every substep of a linear model; without the cache the packed route
made the 6-latent reverse pass 8% slower (dev1).
"""
mutable struct LyapBuffer{TYPE<:Number,TG}
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
    tangents::TG
end

"""
    LyapBuffer(T, n)

A Lyapunov workspace for `n × n` matrices with scalar type `T`, on the Schur
route above `_CTSEM_LYAP_SCHUR_ABOVE` for `Float64`, `Float32` and Duals over
them, and on the packed route otherwise.
"""
function LyapBuffer(::Type{TYPE}, n::Int) where {TYPE<:Number}
    schur = _lyap_schur_route(TYPE, n)
    tri_n = schur ? 0 : gauss_summation(n)
    s = schur ? n : 0
    return LyapBuffer{TYPE,Nothing}(n, schur,
        zeros(TYPE, tri_n, tri_n), zeros(TYPE, tri_n), zeros(Int, tri_n), n, tri_n,
        zeros(TYPE, s, s), zeros(TYPE, s, s), zeros(TYPE, s, s), zeros(TYPE, s, s),
        zeros(TYPE, schur ? max(2n, 4) : 0), zeros(TYPE, schur ? 4 : 0, schur ? 4 : 0),
        false, zeros(TYPE, n, n), false, nothing)
end

function LyapBuffer(::Type{D}, n::Int) where {Tg,V,N,D<:ForwardDiff.Dual{Tg,V,N}}
    schur = _lyap_schur_route(D, n)
    s = schur ? n : 0
    tangents = LyapTangents(LyapBuffer(V, s),
        zeros(V, s, s), zeros(V, s, s), zeros(V, s, s), zeros(V, s, s), zeros(V, s, s),
        zeros(V, s, s), [zeros(V, s, s) for _ in 1:N])
    tri_n = schur ? 0 : gauss_summation(n)
    return LyapBuffer{D,typeof(tangents)}(n, schur,
        zeros(D, tri_n, tri_n), zeros(D, tri_n), zeros(Int, tri_n), n, tri_n,
        zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0),
        zeros(D, 0, 0), false, zeros(D, n, n), false, tangents)
end

"""
Whether `LyapBuffer(T, n)` leaves the packed route: for `Float64` and
`Float32` the Schur route, above `_CTSEM_LYAP_SCHUR_ABOVE`; for a Dual the
tangent route, above `_CTSEM_LYAP_TANGENTS_ABOVE`, its values then solved as
their own type's rule says (a Dual over BigFloat, as the adjoint's BigFloat
checks run, solves its values packed).
"""
_lyap_schur_route(::Type{T}, n::Int) where {T<:Number} =
    T <: LinearAlgebra.BlasReal && n > _CTSEM_LYAP_SCHUR_ABOVE[]
_lyap_schur_route(::Type{<:ForwardDiff.Dual}, n::Int) = n > _CTSEM_LYAP_TANGENTS_ABOVE[]

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

Defined for `Float64` and `Float32`, and for Duals over them by tangents on
the value type's factors (`LyapTangents`), which are the only buffers
`LyapBuffer(T, n)` puts on this route. The route is chosen by a run-time
field, so without the fallback below any other element type (BigFloat) would
compile the whole route as well, although it never takes it.
"""
_lyap_schur!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix, buffer::LyapBuffer,
    transposed::Bool) =
    throw(ArgumentError("the Schur Lyapunov route is for Float64, Float32 and Duals over them"))

function _lyap_schur!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix,
    buffer::LyapBuffer{D}, transposed::Bool) where {Tg,V,N,D<:ForwardDiff.Dual{Tg,V,N}}
    n = buffer.dim
    tg = buffer.tangents
    Av = tg.A
    Qv = tg.Q
    Xv = tg.X
    R = tg.R
    dA = tg.dA
    tmp = tg.tmp
    @inbounds for j in 1:n, i in 1:n
        Av[i, j] = ForwardDiff.value(A[i, j])
        Qv[i, j] = ForwardDiff.value(Q[i, j])
    end
    # `dA` is free scratch for the packed route's transpose here and below.
    _lyap_route!(Xv, Av, Qv, tg.primal, dA, transposed)
    parts = tg.parts
    @inbounds for p in 1:N
        for j in 1:n, i in 1:n
            dA[i, j] = ForwardDiff.partials(A[i, j], p)
        end
        # `dA X + X dA'` is `tmp + tmp'` with `tmp = dA X`, `X` symmetric;
        # transposed, `dA' X + X dA`.
        transposed ? _ctsem_mulTN!(tmp, dA, Xv) : _ctsem_mul!(tmp, dA, Xv)
        for j in 1:n, i in 1:n
            R[i, j] = tmp[i, j] + tmp[j, i] + ForwardDiff.partials(Q[i, j], p)
        end
        _lyap_route!(parts[p], Av, R, tg.primal, dA, transposed)
    end
    @inbounds for j in 1:n, i in 1:n
        X[i, j] = D(Xv[i, j], ForwardDiff.Partials(ntuple(p -> parts[p][i, j], Val(N))))
    end
    return X
end

"""`A X + X A' + Q = 0`, or `A' X + X A + Q = 0` when `transposed`, on `buffer`."""
@inline _lyap_route!(X, A, Q, buffer, At, transposed::Bool) =
    transposed ? my_lyap_transposed!(X, A, Q, buffer, At) : my_lyap!(X, A, Q, buffer)

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
