using LinearAlgebra
using ChainRulesCore

################################################################################
# Adjoint primitive layer
################################################################################
#
# This file holds most of the hand-derived reverse-mode mathematics in the
# package -- the exception is the discrete-intercept solve pulled back inline
# at `adjoint_ekf.jl:591-599` (same math as the `linsolve` derivation below,
# but run at every prediction substep, so it stays hand-inlined against the
# scratch buffers rather than routed through a pure primitive here; see
# `test_adjoint_primitives.jl` for its direct test). Everything else in the
# adjoint path is expressed as a composition of these primitives, which is
# what `docs/src/adjoint-roadmap.md` calls "primitive-level rrules composed by
# replay": if the forward filter changes shape (extra substeps, a new
# predictor kind, a different missingness rule) but keeps expressing its math
# as calls to these kernels, the reverse pass stays correct without a matching
# hand-edit.
#
# Two kernels need hand-written pullbacks here because their forward
# implementations bottom out in LAPACK calls (`dgetrf`, `dtrsyl`) that no
# reverse-mode AD in the Julia ecosystem differentiates today:
#
#   * the matrix exponential   `exp(A)`          -- Padé + scaling/squaring,
#     whose inner Padé solve is an LU factorization;
#   * the continuous Lyapunov solve `A X + X A' + Q = 0` -- either the packed
#     `ksolve!` LU or the Schur/`trsyl!` path, depending on size.
#
# (The plain square solve `A \ B` is the same LAPACK situation, but see the
# note above -- its one caller inlines the pullback rather than using a
# primitive from this file.)
#
# Each is given here as a *pure* function plus a `ChainRulesCore.rrule`. The
# pure form matters: the forward filter's kernels (`my_exp!`, `my_lyap!`,
# `_solve_square_system!`) are in-place and destroy their inputs as documented
# scratch, which is fine for a primal but has no well-defined cotangent. The
# pure wrappers give the reverse pass a value-semantics boundary to attach the
# derivative to.
#
# Derivations (all standard, but written out so they can be checked rather
# than trusted):
#
#   exp:      d/dA exp(A) applied to a direction E is the Fréchet derivative
#             L(A, E), computable as the top-right block of
#             exp([A E; 0 A]). Its adjoint satisfies <Ȳ, L(A,E)> = <L(A', Ȳ), E>,
#             so Ā = L(A', Ȳ) -- one block exponential of twice the dimension.
#
#   lyap:     F(A,Q,X) = A X + X A' + Q = 0. Differentiating,
#             A dX + dX A' = -(dA X + X dA' + dQ). With L(Z) = A Z + Z A',
#             <X̄, dX> = <-L^{-*}(X̄), dA X + X dA' + dQ>, so setting
#             W = -L^{-*}(X̄) -- i.e. W solves A' W + W A + X̄ = 0, which is
#             another Lyapunov solve with A replaced by A' -- gives
#             Ā = W X' + W' X and Q̄ = W. Note the pullback of a Lyapunov solve
#             is a Lyapunov solve, so it costs the same as one forward solve.
#
#   linsolve: X = A^{-1} B, so B̄ = A^{-T} X̄ and Ā = -B̄ X'. Reusing the
#             factorization from the forward solve makes the pullback cheap.
#             (Hand-inlined at `adjoint_ekf.jl:591-599`, not a primitive here
#             -- see the note at the top of this file.)
#
################################################################################

"""
    _ctsem_exp_frechet_block_reference(A, E)

The Fréchet derivative `L(A, E)` of the matrix exponential at `A` in direction
`E`, as the top-right block of `exp([A E; 0 A])`. A reference route, kept so
tests can check `my_exp_frechet!` against something independent of it; the
engine calls `_ctsem_exp_frechet_block` below.

This is the textbook block-exponential identity (Higham, *Functions of
Matrices*, Thm. 10.13). It costs one exponential of a `2n × 2n` matrix, about
eight times the arithmetic of one of size `n`, where the recurrence costs
roughly three.
"""
function _ctsem_exp_frechet_block_reference(A::AbstractMatrix, E::AbstractMatrix)
    n = size(A, 1)
    T = promote_type(eltype(A), eltype(E))

    # Normalise the direction before forming the block, and undo it afterwards.
    #
    # `L(A, E)` is *linear* in `E`, so this is exact -- but it is not cosmetic.
    # `Base.exp` chooses its Padé degree and squaring count from the norm of the
    # matrix it is given, and the block `[A E; 0 A]` inherits the norm of `E`.
    # Cotangents arriving here routinely have norms orders of magnitude above
    # the drift's, so an unnormalised block gets squared several extra times for
    # no reason. Scaling `E` to unit max-norm makes the block's norm that of `A`
    # alone, which is what actually governs how hard this exponential is.
    scale = maximum(abs, E)
    iszero(scale) && return zeros(T, n, n)

    block = zeros(T, 2n, 2n)
    @views begin
        block[1:n, 1:n] .= A
        block[1:n, (n + 1):2n] .= E ./ scale
        block[(n + 1):2n, (n + 1):2n] .= A
    end
    # `_ctsem_expm`, not `Base.exp`: for Float64 the two are the same call, but
    # `Base.exp` has no method for a matrix of duals, which is what this sees
    # when the gradient itself is being differentiated (see `ctsem_hessian`).
    return Matrix(@view _ctsem_expm(block)[1:n, (n + 1):2n]) .* scale
end

"""
    _ctsem_exp_frechet_block(A, E)

The Fréchet derivative `L(A, E)` of the matrix exponential at `A` in direction
`E`, by the Al-Mohy & Higham recurrence in `my_exp_frechet!`. Allocates its
workspace; the reverse pass holds one and calls the kernel directly.
"""
function _ctsem_exp_frechet_block(A::AbstractMatrix, E::AbstractMatrix)
    _CTSEM_OPCOUNT.frechet[] += 1
    T = promote_type(eltype(A), eltype(E))
    n = size(A, 1)
    Y = Matrix{T}(undef, n, n)
    L = Matrix{T}(undef, n, n)
    my_exp_frechet!(Y, L, A, E, ExpFrechetBuffer{T}(n))
    return L
end

"""
    _ctsem_exp_frechet_adjoint(A, Ȳ)

Return `Ā` such that `<Ȳ, d exp(A)> = <Ā, dA>`, i.e. the Fréchet derivative
evaluated at `A'` in direction `Ȳ`.
"""
@inline _ctsem_exp_frechet_adjoint(A::AbstractMatrix, Ȳ::AbstractMatrix) =
    _ctsem_exp_frechet_block(adjoint(A), Ȳ)

################################################################################
# Pure primitives
################################################################################

"""
    _ctsem_expm(A)

Matrix exponential of `A`, as a pure function with a hand-written pullback.

Every scalar type goes through the package's buffered `my_exp!`, `Float64`
included: `Base.exp` bottoms out in a LAPACK solve, and this can run inside
the subject loop. The validation oracle's `ForwardDiff.Dual` path computes
the same Padé approximant by the same code.
"""
function _ctsem_expm(A::Matrix{Float64})
    # `Base.exp` picks its Pade degree adaptively from the matrix norm and is
    # the faster of the two on one core -- but it bottoms out in LAPACK, and
    # LAPACK takes a process-global lock per call (see `small_linalg.jl`). This
    # runs once per Frechet flush in the reverse pass, which is why the reverse
    # tape was the last part of the subject loop that would not thread.
    #
    # `my_exp!` reaches the same Pade approximant through `mul!` and the
    # engine's own LU, neither of which contends, so above one chunk it wins by
    # more than the fixed degree costs.
    n = size(A, 1)
    Y = Matrix{Float64}(undef, n, n)
    scratch = Matrix{Float64}(undef, n, n)
    my_exp!(Y, copy(A), scratch, ExpBuffer{Float64}(n), Val(n))
    return Y
end

function _ctsem_expm(A::AbstractMatrix)
    T = eltype(A)
    n = size(A, 1)
    Y = Matrix{T}(undef, n, n)
    scratch = Matrix{T}(undef, n, n)
    my_exp!(Y, Matrix{T}(A), scratch, ExpBuffer{T}(n), Val(n))
    return Y
end

"""
    _ctsem_lyap(A, Q)

Solve the continuous Lyapunov equation `A X + X A' + Q = 0` for `X`.

`Q` must be symmetric, and `X` is symmetric by construction. Both hold for
every use in the filter, where `Q` is a covariance.

Uses `LyapBuffer`, i.e. exactly the same size-dependent choice the primal
filter makes: a Schur solve for `Float64` systems above
`_CTSEM_LYAP_SCHUR_ABOVE` (ten states), the packed
`ksolve!` otherwise (including every `ForwardDiff.Dual` case, which LAPACK
cannot take).

**This matters far more than it looks.** `ksolve!` solves an
`n(n+1)/2`-square dense system, so its cost is `O(n^6)` against the Schur
path's `O(n^3)`: at 20 latents that is a 210x210 LU instead of a 20x20 Schur
decomposition, per prediction substep. An earlier revision pinned this to
`ksolve!` -- a workaround for a Mooncake evaluation that was then abandoned --
and profiling showed the resulting LU accounting for ~58% of the entire
reverse pass on a 20-latent model.
"""
function _ctsem_lyap(A::AbstractMatrix, Q::AbstractMatrix,
    buffer=LyapBuffer(promote_type(eltype(A), eltype(Q)), size(A, 1)))
    T = promote_type(eltype(A), eltype(Q))
    X = Matrix{T}(undef, size(Q, 1), size(Q, 2))
    my_lyap!(X, Matrix{T}(A), Matrix{T}(Q), buffer)
    return X
end

################################################################################
# Pullbacks
################################################################################

function ChainRulesCore.rrule(::typeof(_ctsem_expm), A::AbstractMatrix)
    Y = _ctsem_expm(A)
    function _ctsem_expm_pullback(Ȳ)
        return (NoTangent(), _ctsem_exp_frechet_adjoint(A, unthunk(Ȳ)))
    end
    return Y, _ctsem_expm_pullback
end

"""
    _ctsem_lyap_pullback(A, X, X̄, buffer)

Return `(Ā, Q̄)` for `X = _ctsem_lyap(A, Q)`.

`Q̄` is the solution `W` of the transposed Lyapunov equation
`A' W + W A + X̄ = 0`, and `Ā = W X' + W' X`.

**The incoming cotangent is symmetrised first, and it must be.** `ksolve!`
solves the Lyapunov equation over *symmetric* unknowns only (it packs the
system into `n(n+1)/2` triangular unknowns), so handing it a non-symmetric
right-hand side does not fail -- it quietly returns the symmetric matrix
nearest to solving it, which is not the `W` this pullback needs. Symmetrising
is also the mathematically correct thing to do rather than a workaround:
`_ctsem_lyap` always returns a symmetric `X`, so `dX` is symmetric and
`<X̄, dX> = <sym(X̄), dX>` for any `X̄`.
"""
function _ctsem_lyap_pullback(A::AbstractMatrix, X::AbstractMatrix, X̄::AbstractMatrix,
    buffer=LyapBuffer(promote_type(eltype(A), eltype(X), eltype(X̄)), size(A, 1)))
    T = promote_type(eltype(A), eltype(X), eltype(X̄))
    n = size(A, 1)
    return _ctsem_lyap_pullback!(Matrix{T}(undef, n, n), Matrix{T}(undef, n, n),
        Matrix{T}(undef, n, n), Matrix{T}(undef, n, n), A, X, X̄, buffer)
end

"""
    _ctsem_lyap_pullback!(Ā, W, rhs, At, A, X, X̄, buffer)

The pullback above into caller-owned storage: `Ā` and `W` receive the two
cotangents, `rhs` and `At` are scratch. Four fresh matrices per prediction
substep is what the allocating form cost (`Profile.Allocs`, dev1).
"""
function _ctsem_lyap_pullback!(Ā::AbstractMatrix{T}, W::AbstractMatrix{T},
    rhs::AbstractMatrix{T}, At::AbstractMatrix{T},
    A::AbstractMatrix, X::AbstractMatrix, X̄::AbstractMatrix, buffer) where {T}
    n = size(A, 1)
    @inbounds for j in 1:n, i in 1:n
        rhs[i, j] = (X̄[i, j] + X̄[j, i]) / 2
    end
    # A transposed copy rather than `adjoint(A)`: `my_lyap!` wants a matrix it
    # can compare and factor, and a product with an `adjoint` operand does not
    # reach `gemm` and does not thread (0.18x on 23 threads at this size),
    # which left this the last serialised piece of the reverse tape. See
    # `small_linalg.jl`.
    @inbounds for j in 1:n, i in 1:n
        At[j, i] = A[i, j]
    end
    my_lyap!(W, At, rhs, buffer)
    _ctsem_mulNT!(Ā, W, X)
    _ctsem_mulTN!(Ā, W, X, one(T), one(T))
    return Ā, W
end

function ChainRulesCore.rrule(::typeof(_ctsem_lyap), A::AbstractMatrix, Q::AbstractMatrix)
    X = _ctsem_lyap(A, Q)
    function _ctsem_lyap_rrule_pullback(X̄)
        Ā, Q̄ = _ctsem_lyap_pullback(A, X, unthunk(X̄))
        return (NoTangent(), Ā, Q̄)
    end
    return X, _ctsem_lyap_rrule_pullback
end

################################################################################
# SD/correlation-square-root -> covariance
################################################################################
#
# `sdcovsqrt2cov!` is the fourth primitive needing a hand-written pullback, but
# for a different reason from the three above: it is differentiable by AD, it
# is just expensive to differentiate naively. The map is d² -> d², so a
# ForwardDiff Jacobian costs O(d^5) -- fine for a 3-latent model, ruinous for
# the 20-60 state augmented models that motivate having an adjoint at all.
#
# The saving structure is that `constraincorsqrt1_vec!` is **row-separable**:
# row `i` of its output depends only on row/column `i` of the (symmetric)
# input, because the row sums `s[i]`, `ss[i]`, the scale `r[i]` and the
# diagonal repair `sqrt(1 - rowsum + eps)` all read one row. So its Jacobian is
# block diagonal with d blocks of size d×d, and the pullback costs O(d³) --
# the same order as the covariance product it feeds.
#
# The composition is
#     O   = constraincorsqrt1_vec(mat)        (row-separable, nonlinear)
#     B   = Diagonal(diag(mat)) * O           (row scaling by the SDs)
#     cov = B * B'
# so, given a covariance cotangent `C̄`,
#     B̄        = (C̄ + C̄') * B
#     Ō[i,j]   = mat[i,i] * B̄[i,j]
#     mat̄[i,i] += Σ_j B̄[i,j] * O[i,j]        (the SD path)
# and each row of `Ō` is pushed back through one row of the row-separable map.

"""
    _ctsem_corrsqrt_row(v, i, epsilon)

Row `i` of `constraincorsqrt1_vec!`'s output, as a pure function of `v`, the
`i`th symmetric row of the parameter matrix (`v[j] == mat[i,j]` for `j < i`
and `mat[j,i]` for `j > i`; `v[i]` is read but cancels, since the diagonal is
excluded from every sum below and the diagonal output is recomputed).

This mirrors `constraincorsqrt1_vec!` statement for statement; it exists so the
pullback differentiates the function the filter actually evaluates rather than
an independently rederived formula. `test_constrain_cor_sqrt.jl`'s "corrsqrt
row mirror matches the buffered primal it claims to track" asserts the two
agree, row by row, against `constraincorsqrt1_vec!`'s own buffered output.
"""
_ctsem_corrsqrt_row(v::AbstractVector{T}, i::Int, epsilon) where {T} =
    _ctsem_corrsqrt_row!(similar(v), v, i, epsilon)

"""In-place form of `_ctsem_corrsqrt_row`: the row is written into `out`."""
function _ctsem_corrsqrt_row!(out::AbstractVector{T}, v::AbstractVector{T}, i::Int, epsilon) where {T}
    d = length(v)
    e = convert(T, epsilon)
    si = e
    ssi = e
    @inbounds for j in 1:d
        j == i && continue
        si += v[j]
        ssi += v[j] * v[j]
    end
    abs_si = abs(si)
    tmp = sqrt(log1p(exp(2 * (abs_si - si - one(T)) - 4)))
    r = sqrt(ssi + (tmp * (abs_si / sqrt(ssi) - one(T)) + one(T)) * tmp + one(T))
    sq = zero(T)
    @inbounds for j in 1:d
        if j == i
            out[j] = zero(T)
        else
            o = v[j] / r
            out[j] = o
            sq += o * o
        end
    end
    @inbounds out[i] = sqrt(one(T) - sq + e)
    return out
end

"""
    _ctsem_symmetric_row!(v, mat, i, d)

Fill `v` with the `i`th row of `mat` read symmetrically from its lower
triangle, matching `_sym_lower_get`.
"""
@inline function _ctsem_symmetric_row!(v::AbstractVector, mat::AbstractMatrix, i::Int, d::Int)
    @inbounds for j in 1:d
        v[j] = j >= i ? mat[j, i] : mat[i, j]
    end
    return v
end

"""
    CTSEMCovSqrtScratch(T, d)

Working storage for `_sdcovsqrt2cov_pullback!`, sized for the largest `d` it
will see. The pullback runs at every prediction substep (DIFFUSION) and
every observed row (MANIFESTVAR), and allocated seven arrays each time --
about half the reverse pass's bytes on a small linear model
(`Profile.Allocs`, dev1) -- so callers on the hot path hand it one of these.
"""
struct CTSEMCovSqrtScratch{T}
    O::Matrix{T}
    B::Matrix{T}
    Csym::Matrix{T}
    Bbar::Matrix{T}
    v::Vector{T}
    Obar_row::Vector{T}
    vbar::Vector{T}
end

CTSEMCovSqrtScratch(::Type{T}, d::Int) where {T} = CTSEMCovSqrtScratch{T}(
    zeros(T, d, d), zeros(T, d, d), zeros(T, d, d), zeros(T, d, d),
    zeros(T, d), zeros(T, d), zeros(T, d))

# With no scratch, allocate as the pullback always did (tests, the rrule).
_covsqrt_scratch(::Nothing, ::Type{T}, d::Int) where {T} = (
    Matrix{T}(undef, d, d), Matrix{T}(undef, d, d), Matrix{T}(undef, d, d),
    Matrix{T}(undef, d, d), Vector{T}(undef, d), Vector{T}(undef, d), Vector{T}(undef, d))

function _covsqrt_scratch(sc::CTSEMCovSqrtScratch{T}, ::Type{T}, d::Int) where {T}
    return (view(sc.O, 1:d, 1:d), view(sc.B, 1:d, 1:d), view(sc.Csym, 1:d, 1:d),
        view(sc.Bbar, 1:d, 1:d), view(sc.v, 1:d), view(sc.Obar_row, 1:d), view(sc.vbar, 1:d))
end

"""
    _sdcovsqrt2cov_pullback!(mat_bar, mat, cov_bar, d; epsilon=1e-5, scratch=nothing)

Accumulate into `mat_bar` the cotangent of `mat` for
`cov = sdcovsqrt2cov(mat)`, given the covariance cotangent `cov_bar`.

`mat_bar` is *added to*, not overwritten, so several uses of the same
parameter matrix across rows accumulate correctly. Only the lower triangle and
diagonal of `mat_bar` are written, matching where the free parameters live.
"""
function _sdcovsqrt2cov_pullback!(mat_bar::AbstractMatrix, mat::AbstractMatrix,
    cov_bar::AbstractMatrix, d::Int; epsilon::Real=1e-5, scratch=nothing)
    d == 0 && return mat_bar
    # A zero cotangent happens routinely -- e.g. the manifest-covariance
    # cotangent on a fully missing row, where no measurement update ran -- and
    # this map is O(d^3) to pull back, so it is worth not doing.
    all(iszero, cov_bar) && return mat_bar
    T = promote_type(eltype(mat), eltype(cov_bar))

    # Recompute O and B. The forward pass overwrites its correlation-factor
    # buffer with the finished covariance, so O is not recoverable from the
    # workspace; recomputing it is one O(d²) pass and keeps the trace small.
    O, B, Csym, Bbar, v, Obar_row, vbar = _covsqrt_scratch(scratch, T, d)
    @inbounds for i in 1:d
        _ctsem_symmetric_row!(v, mat, i, d)
        _ctsem_corrsqrt_row!(Obar_row, v, i, epsilon)
        for j in 1:d
            O[i, j] = Obar_row[j]
        end
    end
    @inbounds for j in 1:d, i in 1:d
        B[i, j] = mat[i, i] * O[i, j]
    end

    # B̄ = (C̄ + C̄') B. cov is symmetric, so this is the correct pullback
    # whether or not the incoming cotangent has been symmetrized.
    @inbounds for j in 1:d, i in 1:d
        Csym[i, j] = cov_bar[i, j] + cov_bar[j, i]
    end
    mul!(Bbar, Csym, B)

    # SD path: mat[i,i] scales the whole of row i of O.
    @inbounds for i in 1:d
        acc = zero(T)
        for j in 1:d
            acc += Bbar[i, j] * O[i, j]
        end
        mat_bar[i, i] += acc
    end

    # Correlation path, one row at a time (the block-diagonal structure).
    @inbounds for i in 1:d
        for j in 1:d
            Obar_row[j] = mat[i, i] * Bbar[i, j]
        end
        _ctsem_symmetric_row!(v, mat, i, d)
        _ctsem_corrsqrt_row_pullback!(vbar, v, Obar_row, i, epsilon)
        # Scatter the symmetric row cotangent back into the lower triangle.
        for j in 1:d
            j == i && continue
            if j < i
                mat_bar[i, j] += vbar[j]
            else
                mat_bar[j, i] += vbar[j]
            end
        end
    end
    return mat_bar
end

"""
    _ctsem_corrsqrt_row_scale(s, ss)

The row scale `r` of `constraincorsqrt1_vec!`, isolated as a function of just
the row sum `s` and the row sum-of-squares `ss`.

Everything awkward in that kernel (the `sqrt(log1p(exp(...)))` softening term)
depends on the row only through these two scalars, which is what makes the
row pullback below cheap: two partial derivatives, not `d` of them.
"""
@inline function _ctsem_corrsqrt_row_scale(s, ss)
    one_ = one(s)
    abs_s = abs(s)
    tmp = sqrt(log1p(exp(2 * (abs_s - s - one_) - 4)))
    return sqrt(ss + (tmp * (abs_s / sqrt(ss) - one_) + one_) * tmp + one_)
end

"""
    _ctsem_corrsqrt_row_pullback!(vbar, v, obar, i, epsilon)

Write into `vbar` the cotangent of `v` for one row of
`constraincorsqrt1_vec!`, given the row's output cotangent `obar`.

Hand-derived, and `O(d)` rather than the `O(d^2)`-with-AD-overhead of taking a
`ForwardDiff` gradient of the whole row map -- profiling put that at ~21% of
the reverse pass. Only the scale `r` needs AD, through two scalar partials.

Forward, writing `S = {j : j != i}`:

    s   = eps + sum_{j in S} v[j]
    ss  = eps + sum_{j in S} v[j]^2
    r   = R(s, ss)
    o[j] = v[j] / r         (j in S)
    sq  = sum_{j in S} o[j]^2
    o[i] = sqrt(1 - sq + eps)

Reverse, in that order backwards. `test_adjoint_primitives.jl` checks this
against `ForwardDiff` applied to `_ctsem_corrsqrt_row`, which is kept as the
readable reference implementation.
"""
function _ctsem_corrsqrt_row_pullback!(vbar::AbstractVector{T}, v::AbstractVector,
    obar::AbstractVector, i::Int, epsilon) where {T}
    d = length(v)
    e = convert(T, epsilon)
    s = e
    ss = e
    @inbounds for j in 1:d
        j == i && continue
        s += v[j]
        ss += v[j] * v[j]
    end

    # r and its two partials, from a single dual evaluation.
    Dual2 = ForwardDiff.Dual{Nothing,T,2}
    r_dual = _ctsem_corrsqrt_row_scale(
        Dual2(s, ForwardDiff.Partials((one(T), zero(T)))),
        Dual2(ss, ForwardDiff.Partials((zero(T), one(T)))))
    r = ForwardDiff.value(r_dual)
    dr_ds = ForwardDiff.partials(r_dual, 1)
    dr_dss = ForwardDiff.partials(r_dual, 2)
    inv_r = inv(r)

    # o[i] = sqrt(1 - sq + eps): pull its cotangent back onto sq, then onto
    # the off-diagonal outputs it was built from.
    sq = zero(T)
    @inbounds for j in 1:d
        j == i && continue
        o = v[j] * inv_r
        sq += o * o
    end
    diagonal = sqrt(one(T) - sq + e)
    sq_bar = -obar[i] * T(0.5) / diagonal

    # o[j] = v[j] / r, for j != i.
    r_bar = zero(T)
    @inbounds for j in 1:d
        if j == i
            vbar[j] = zero(T)
        else
            o = v[j] * inv_r
            o_bar = obar[j] + 2 * sq_bar * o
            vbar[j] = o_bar * inv_r
            r_bar -= o_bar * o * inv_r
        end
    end

    # r = R(s, ss); s and ss are plain sums over j != i.
    s_bar = r_bar * dr_ds
    ss_bar = r_bar * dr_dss
    @inbounds for j in 1:d
        j == i && continue
        vbar[j] += s_bar + 2 * ss_bar * v[j]
    end
    return vbar
end
