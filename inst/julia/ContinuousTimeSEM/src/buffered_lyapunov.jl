using LinearAlgebra

################################################################################
# Lyapunov related functions
################################################################################

"""
    LyapTangents{V, B}

The Lyapunov solve at a `ForwardDiff.Dual` element type: the value parts of
`A` and `Q` (`A`, `Q`), solved on a buffer of the value type (`primal`; a
nested Dual recurses), then each partial as the same equation on the same
factors -- differentiating `A X + X A' + Q = 0`,

    A dX + dX A' + (dA X + X dA' + dQ) = 0,

so `dX` is a Lyapunov solve with `Q` replaced by `dA X + X dA' + dQ`. Exact,
and a solve against cached factors per partial, where pushing the Duals
through the factorisation would repeat it once per partial. `X`, `R` (that
right-hand side), `dA` and `tmp` are scratch; `parts` holds the partials'
solutions.
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
by the real Schur factors of `A` (`real_schur.jl`, Bartels-Stewart), at every
size and every floating-point element type; a `ForwardDiff.Dual` solves its
values so and each partial on the same factors (`tangents`, `LyapTangents`),
and is `nothing` otherwise.

There is one route. A packed solve over the `k(k+1)/2` symmetric unknowns
served small blocks until 2026-10: O(k^6) to factor against Schur's O(k^3),
it won only below about eight states, by a microsecond or two a solve
(local i9-12900KS, `scripts/stangrad-2026-10/bench_lyap.jl`), on models whose
gradients take milliseconds, and it measured the operator's conditioning for
the series route by its LU pivots where Schur reads the eigenvalue sums the
solve divides by.

The factors are kept across calls and rebuilt only when `A` changes
(`last_A`; `_blocks_identical` compares Dual partials too, so a factor never
serves a point it was not built at). The reverse pass solves against one `A`
at every substep of a linear model, and the pullback's transposed solve uses
the forward factors.
"""
mutable struct LyapBuffer{TYPE<:Number,TG}
    dim::Int
    # `A = U S U'`, the right-hand side and then the solution in Schur
    # coordinates (`C`), product scratch, the reduction's reflector (and a
    # block solve's right-hand side) in `work`, a block's Kronecker system in
    # `small`, and whether the QR iteration converged. Empty at a Dual type,
    # whose factors are its value buffer's.
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

A Lyapunov workspace for `n × n` matrices with scalar type `T`: a
floating-point type, or a `ForwardDiff.Dual` over one.
"""
function LyapBuffer(::Type{T}, n::Int) where {T<:AbstractFloat}
    return LyapBuffer{T,Nothing}(n,
        zeros(T, n, n), zeros(T, n, n), zeros(T, n, n), zeros(T, n, n),
        zeros(T, max(2n, 4)), zeros(T, 4, 4), false, zeros(T, n, n), false, nothing)
end

function LyapBuffer(::Type{D}, n::Int) where {Tg,V,N,D<:ForwardDiff.Dual{Tg,V,N}}
    tangents = LyapTangents(LyapBuffer(V, n),
        zeros(V, n, n), zeros(V, n, n), zeros(V, n, n), zeros(V, n, n), zeros(V, n, n),
        zeros(V, n, n), [zeros(V, n, n) for _ in 1:N])
    return LyapBuffer{D,typeof(tangents)}(n,
        zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0, 0), zeros(D, 0),
        zeros(D, 0, 0), false, zeros(D, 0, 0), false, tangents)
end

"""
    LyapBuffer(mat::AbstractMatrix)

Create a Lyapunov workspace sized and typed from `mat`.
"""
LyapBuffer(mat::AbstractMatrix) = LyapBuffer(eltype(mat), size(mat, 1))

"""
    my_lyap!(X, A, Q, buffer::LyapBuffer)

Solve `A * X + X * A' + Q = 0` on the leading `buffer.dim` block. The
solution is written to that block of `X`, symmetric; `Q` must be symmetric,
as every covariance the filter hands it is, since the solve works on one
triangle and mirrors it.
"""
my_lyap!(X::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, Q::AbstractMatrix{TYPE},
    buffer::LyapBuffer{TYPE}) where {TYPE<:Number} = _lyap_solve!(X, A, Q, buffer, false)

"""
    my_lyap_transposed!(W, A, Q, buffer)

Solve `A' W + W A + Q = 0`, the equation a Lyapunov pullback solves, on `A`'s
own factors (`A' = U S' U'`), so the pullback factorises only where `A`
changes, as the forward pass does.
"""
my_lyap_transposed!(W::AbstractMatrix{TYPE}, A::AbstractMatrix, Q::AbstractMatrix{TYPE},
    buffer::LyapBuffer{TYPE}) where {TYPE<:Number} = _lyap_solve!(W, A, Q, buffer, true)

"""
    _lyap_solve!(X, A, Q, buffer, transposed)

Factor `A` where it has changed, then `C = -U' Q U`, the quasi triangular
solve (`_ctsem_quasitri_lyap!`, `S Y + Y S' = C`, or `S' Y + Y S = C` when
`transposed`) and `X = U Y U'`, symmetrised. Every product goes through the
engine's kernels (`small_linalg.jl`), so nothing takes a library lock. Where
the QR iteration did not converge, or `A` was not finite, `X` is NaN.
"""
function _lyap_solve!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix,
    buffer::LyapBuffer{T}, transposed::Bool) where {T<:AbstractFloat}
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

function _lyap_solve!(X::AbstractMatrix, A::AbstractMatrix, Q::AbstractMatrix,
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
    _lyap_solve!(Xv, Av, Qv, tg.primal, transposed)
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
        _lyap_solve!(parts[p], Av, R, tg.primal, transposed)
    end
    @inbounds for j in 1:n, i in 1:n
        X[i, j] = D(Xv[i, j], ForwardDiff.Partials(ntuple(p -> parts[p][i, j], Val(N))))
    end
    return X
end
