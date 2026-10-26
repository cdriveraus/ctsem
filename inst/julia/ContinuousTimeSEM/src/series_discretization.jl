using LinearAlgebra

################################################################################
# The discretisation where the closed forms divide by the drift
################################################################################
#
# `_compute_discrete_time_form!` builds the discrete intercept and process noise
# of one interval from closed forms:
#
#     dINT  = J^-1 (e^{J dt} - I) c
#     dDIFF = X - e^{J dt} X e^{J' dt},   J X + X J' + Q = 0
#
# Both are exact, and both are removable singularities: the first divides by
# the eigenvalues of J, the second by their pairwise sums. The quantities
# themselves are smooth everywhere -- they are the integrals
#
#     Phi(dt) = int_0^dt e^{J s} ds,          dINT  = Phi c
#     V(dt)   = int_0^dt e^{J s} Q e^{J' s} ds, dDIFF = V
#
# -- but computed through the closed forms they lose accuracy as a pivot
# approaches zero, relative to the interval. With x the smallest pivot times
# dt, the value is good to about eps / x, the gradient to eps / x^2, a
# curvature to eps / x^3. At x = 1e-5 a curvature has no correct digits, at
# x = 1e-8 neither has the gradient, and at x = 0 the closed forms divide by
# zero.
#
# Nothing about the model is wrong there. A drift with a unit root is a random
# walk, a drift with a positive eigenvalue an explosive process, and the
# likelihood is smooth through the boundary between them. It was found through
# a Laplace fit whose random cross effect was symmetric: one unit's DRIFT went
# singular at its inner mode, a1 * a2 - c^2 crossing zero, and the inner Newton
# climbed the rounding spike the closed forms produce there until the
# Hessian's differencing points failed. Short intervals reach the same place
# from the other side, since x scales with dt: a substep of 1e-3 on a drift of
# -0.5 has x = 5e-4 whatever the conditioning.
#
# So when the pivot is small relative to the interval, the discretisation takes
# the integrals directly instead, by a truncated Taylor series on a scaled-down
# interval and then doubling:
#
#     Phi(2h) = Phi(h) + e^{J h} Phi(h)
#     V(2h)   = V(h) + e^{J h} V(h) e^{J' h}
#
# Neither recurrence divides by anything, and both are generic in the element
# type, so the forward-mode paths differentiate them as they stand. The
# adjoint differentiates them by hand, below, recomputing the forward sweep
# from the recorded inputs rather than taping it.
#
# Every product here goes through `_ctsem_mul!` and its transposed forms: no
# LAPACK, and no BLAS below the small-product threshold, for the reasons in
# `small_linalg.jl`.

"""
Conditioning, in units of the interval, below which the discretisation takes
the series route.

Compared against the smallest LU pivot of the intercept block times `dt`, and
separately against the smallest pivot of the packed Lyapunov operator times
`dt`, so each piece switches on its own. At 0.01 the closed forms are good to
about 2e-14 in value, 2e-12 in gradient and 2e-10 in curvature where they hand
over, which is also the size of the step the switch makes: small enough for a
central-difference Hessian with steps of 1e-4 to see about 0.02 in the worst
case, where a curvature that feeds the objective -- a Laplace log determinant
-- straddles the switch. `ctsem_set_series_discretization_below!` moves it; 0
restores the closed forms everywhere and `Inf` takes the series everywhere.
"""
const _CTSEM_SERIES_BELOW = Ref(0.01)

"""Set the conditioning below which the discretisation takes the series route."""
function ctsem_set_series_discretization_below!(x::Real)
    x >= 0 || throw(ArgumentError("threshold must be non-negative"))
    _CTSEM_SERIES_BELOW[] = Float64(x)
    return _CTSEM_SERIES_BELOW[]
end

# The scaled interval's norm bound, and the most terms the series may take.
# With || J h || at most 0.25 the noise series, whose operator norm is twice
# that, needs 18 terms for its third derivatives to be accurate to 1e-17; see
# `_series_terms`.
const _CTSEM_SERIES_THETA = 0.25
const _CTSEM_SERIES_MAXTERMS = 20
# Doublings are capped so that an absurd trial point cannot grow the tapes
# without bound. Past || J dt || = 0.25 * 2^64, about 5e18, the scaled step no
# longer meets the Taylor bound and the result is inaccurate -- a point no fit
# has any business evaluating.
const _CTSEM_SERIES_MAXDOUBLINGS = 64

"""
    SeriesKernelBuffer{T}(n)

Scratch and tapes for the series intercept and noise kernels.

One buffer serves both kernels, one at a time, on the leading `p x p` block of
`n x n` storage. The tapes -- Horner iterates, and the exponentials and
integrals at each doubling -- grow on first use to the number of terms and
doublings an interval needs, and are reused after that.

`intercept` and `noise` record which kernels the last discretisation used, so
that the forward filter's tape can say which pullback its substep needs. `Q`
holds the noise kernel's diffusion block for the same reason: the reverse pass
needs the input, where the closed form needed only the Lyapunov solution.
"""
mutable struct SeriesKernelBuffer{T}
    n::Int
    X::Matrix{T}
    Q::Matrix{T}
    out::Matrix{T}
    horner::Vector{Matrix{T}}
    taylor::Vector{Matrix{T}}
    Es::Vector{Matrix{T}}
    Ps::Vector{Matrix{T}}
    b1::Matrix{T}
    b2::Matrix{T}
    b3::Matrix{T}
    Ebar::Matrix{T}
    Xbar::Matrix{T}
    s::Int
    m::Int
    h0::Float64
    intercept::Bool
    noise::Bool
end

function SeriesKernelBuffer(::Type{T}, n::Int) where {T}
    z() = zeros(T, n, n)
    return SeriesKernelBuffer{T}(n, z(), z(), z(), Matrix{T}[], Matrix{T}[],
        Matrix{T}[], Matrix{T}[], z(), z(), z(), z(), z(), 0, 0, 0.0, false, false)
end

@inline function _series_grow!(tape::Vector{Matrix{T}}, count::Int, n::Int) where {T}
    while length(tape) < count
        push!(tape, zeros(T, n, n))
    end
    return tape
end

@inline _series_block(A::AbstractMatrix, p::Int) = view(A, 1:p, 1:p)

"""
    CTSEMSeriesChain{T}()

One substep's series kernel as its pullback reads it: the plan (`s`, `m`,
`h0`), the doubling chain's exponentials and integrals (`Es`, `Ps`) and the
result (`out`), each `p x p`. The adjoint tape keeps one per kernel a substep
used (`CTSEMPredictRecord`), so the reverse pass loads the chain the forward
pass computed (`_series_chain_load!`) where it ran the kernel again: that was
a third of the chain's products in a gradient, and the chains were a third of
a Laplace value and gradient on the affect pilot.
"""
mutable struct CTSEMSeriesChain{T}
    s::Int
    m::Int
    h0::Float64
    Es::Vector{Matrix{T}}
    Ps::Vector{Matrix{T}}
    out::Matrix{T}
end

CTSEMSeriesChain(::Type{T}) where {T} = CTSEMSeriesChain{T}(0, 0, 0.0, Matrix{T}[],
    Matrix{T}[], zeros(T, 0, 0))

"""Copy the leading `p x p` blocks of `buf`'s last kernel into `chain`."""
function _series_chain_store!(chain::CTSEMSeriesChain{T}, buf::SeriesKernelBuffer{T},
    p::Int) where {T}
    chain.s, chain.m, chain.h0 = buf.s, buf.m, buf.h0
    _series_grow!(chain.Es, buf.s, p)
    _series_grow!(chain.Ps, buf.s, p)
    size(chain.out, 1) == p || (chain.out = zeros(T, p, p))
    @inbounds for k in 1:buf.s
        _copy_block!(chain.Es[k], buf.Es[k], p, p)
        _copy_block!(chain.Ps[k], buf.Ps[k], p, p)
    end
    _copy_block!(chain.out, buf.out, p, p)
    return chain
end

"""Load `chain` into `buf` as its kernel left it, for the pullbacks."""
function _series_chain_load!(buf::SeriesKernelBuffer{T}, chain::CTSEMSeriesChain{T},
    p::Int) where {T}
    buf.s, buf.m, buf.h0 = chain.s, chain.m, chain.h0
    _series_grow!(buf.Es, chain.s, buf.n)
    _series_grow!(buf.Ps, chain.s, buf.n)
    @inbounds for k in 1:chain.s
        _copy_block!(buf.Es[k], chain.Es[k], p, p)
        _copy_block!(buf.Ps[k], chain.Ps[k], p, p)
    end
    _copy_block!(buf.out, chain.out, p, p)
    return buf
end

"""The 1-norm of the leading `p x p` block of `J`, on primal values."""
@inline function _series_norm1(J::AbstractMatrix, p::Int)
    nrm = 0.0
    @inbounds for j in 1:p
        col = 0.0
        for i in 1:p
            col += abs(Float64(_primal(J[i, j])))
        end
        nrm = max(nrm, col)
    end
    return nrm
end

"""
Terms of the Taylor series for a scaled norm `z`.

The smallest `m` for which the tail's third derivative, about
`z^(m-2) / (m-2)!`, is below 1e-17: third because that is the deepest the
engine differentiates the likelihood (the Laplace route's outer gradient,
through its inner curvature).
"""
function _series_terms(z::Float64)
    k = 2
    term = z^2 / 2
    while term > 1e-17 && k < _CTSEM_SERIES_MAXTERMS - 2
        k += 1
        term *= z / k
    end
    return k + 2
end

"""Doublings and scaled step for `h * J`, and the Taylor terms for `factor * z`."""
function _series_plan!(buf::SeriesKernelBuffer, J::AbstractMatrix, h, p::Int,
    factor::Float64)
    hval = Float64(_primal(h))
    nrm = _series_norm1(J, p) * hval
    s = 0
    if isfinite(nrm) && nrm > _CTSEM_SERIES_THETA
        s = min(ceil(Int, log2(nrm / _CTSEM_SERIES_THETA)), _CTSEM_SERIES_MAXDOUBLINGS)
    end
    scale = 2.0^s
    h0 = hval / scale
    z = isfinite(nrm) ? factor * nrm / scale : 0.0
    buf.s = s
    buf.m = _series_terms(z)
    buf.h0 = h0
    n = buf.n
    _series_grow!(buf.horner, buf.m, n)
    _series_grow!(buf.Es, s, n)
    _series_grow!(buf.Ps, s, n)
    X = _series_block(buf.X, p)
    @inbounds for j in 1:p, i in 1:p
        X[i, j] = J[i, j] * h0
    end
    return buf
end

@inline function _series_add_identity!(A::AbstractMatrix, p::Int)
    @inbounds for i in 1:p
        A[i, i] += one(eltype(A))
    end
    return A
end

"""
    _series_intercept!(buf, J, h, p)

`Phi(h) = int_0^h e^{J s} ds` for the leading `p x p` block of `J`, into the
same block of `buf.out`, taping what `_series_intercept_pullback!` needs.

At the scaled interval `h0`, `Phi = h0 G` and `e^{J h0} = I + X G`, with
`X = J h0` and `G` the series of `phi_1(X) = sum_j X^j / (j+1)!` by Horner's
rule: `G_m = I + X / (m+1)`, `G_j = I + X G_{j+1} / (j+1)`, `G = G_1`.
"""
function _series_intercept!(buf::SeriesKernelBuffer{T}, J::AbstractMatrix, h,
    p::Int) where {T}
    _CTSEM_OPCOUNT.series_intercept[] += 1
    _series_plan!(buf, J, h, p, 1.0)
    s, m, h0 = buf.s, buf.m, buf.h0
    X = _series_block(buf.X, p)
    Gm = _series_block(buf.horner[m], p)
    @inbounds for j in 1:p, i in 1:p
        Gm[i, j] = X[i, j] / (m + 1)
    end
    _series_add_identity!(Gm, p)
    for j in (m - 1):-1:1
        Gj = _series_block(buf.horner[j], p)
        _ctsem_mul!(Gj, X, _series_block(buf.horner[j + 1], p), inv(T(j + 1)))
        _series_add_identity!(Gj, p)
    end
    G = _series_block(buf.horner[1], p)
    out = _series_block(buf.out, p)
    first = s == 0 ? out : _series_block(buf.Ps[1], p)
    @inbounds for j in 1:p, i in 1:p
        first[i, j] = G[i, j] * h0
    end
    s == 0 && return out
    E0 = _series_block(buf.Es[1], p)
    _ctsem_mul!(E0, X, G)
    _series_add_identity!(E0, p)
    for k in 1:s
        E = _series_block(buf.Es[k], p)
        P = _series_block(buf.Ps[k], p)
        next = k == s ? out : _series_block(buf.Ps[k + 1], p)
        _ctsem_mul!(next, E, P)
        next .+= P
        k < s && _ctsem_mul!(_series_block(buf.Es[k + 1], p), E, E)
    end
    return out
end

"""
    _series_intercept_pullback!(Jbar, buf, p, Phibar)

Add to `Jbar` the pullback of `_series_intercept!` for cotangent `Phibar`.
`buf` must hold that kernel's forward sweep, unchanged since.
"""
function _series_intercept_pullback!(Jbar::AbstractMatrix, buf::SeriesKernelBuffer{T},
    p::Int, Phibar::AbstractMatrix) where {T}
    s, m, h0 = buf.s, buf.m, buf.h0
    X = _series_block(buf.X, p)
    t1 = _series_block(buf.b2, p)
    t2 = _series_block(buf.b3, p)
    Pb, Eb = _series_intercept_doubling_pullback!(buf, p, Phibar)
    # Phi_0 = h0 G and E_0 = I + X G.
    G = _series_block(buf.horner[1], p)
    Gb = t2
    Xb = _series_block(buf.Xbar, p)
    @inbounds for j in 1:p, i in 1:p
        Gb[i, j] = Pb[i, j] * h0
    end
    if s > 0
        _ctsem_mulTN!(Gb, X, Eb, one(T), one(T))
        _ctsem_mulNT!(Xb, Eb, G)
    else
        fill!(Xb, zero(T))
    end
    # G_j = I + X G_{j+1} / (j+1), with G_{m+1} = I.
    for j in 1:m
        c = inv(T(j + 1))
        if j < m
            _ctsem_mulNT!(Xb, Gb, _series_block(buf.horner[j + 1], p), c, one(T))
            _ctsem_mulTN!(t1, X, Gb, c)
            copyto!(Gb, t1)
        else
            @inbounds for jj in 1:p, ii in 1:p
                Xb[ii, jj] += c * Gb[ii, jj]
            end
        end
    end
    @inbounds for j in 1:p, i in 1:p
        Jbar[i, j] += Xb[i, j] * h0
    end
    return Jbar
end

"""
    _series_intercept_doubling_pullback!(buf, p, Phibar)

The doubling chain of an intercept series reversed, for cotangent `Phibar` on
its result: the cotangents of the first scaled integral `Phi(h0)` (in `buf.b1`)
and of the first exponential `e^{J h0}` (in `buf.Ebar`, zero without doublings).
`buf` must hold the forward sweep, by either `_series_intercept!` or
`_series_intercept_table!`; `buf.b2` and `buf.b3` are its scratch.
"""
function _series_intercept_doubling_pullback!(buf::SeriesKernelBuffer{T}, p::Int,
    Phibar::AbstractMatrix) where {T}
    s = buf.s
    Pb = _series_block(buf.b1, p)
    Eb = _series_block(buf.Ebar, p)
    t1 = _series_block(buf.b2, p)
    t2 = _series_block(buf.b3, p)
    copyto!(Pb, Phibar)
    fill!(Eb, zero(T))
    # Phi_k = Phi_{k-1} + E_{k-1} Phi_{k-1}, and E_k = E_{k-1}^2 while k < s.
    for k in s:-1:1
        E = _series_block(buf.Es[k], p)
        P = _series_block(buf.Ps[k], p)
        if k < s
            _ctsem_mulNT!(t1, Eb, E)
            _ctsem_mulTN!(t1, E, Eb, one(T), one(T))
        else
            fill!(t1, zero(T))
        end
        _ctsem_mulNT!(t1, Pb, P, one(T), one(T))
        _ctsem_mulTN!(t2, E, Pb)
        Pb .+= t2
        copyto!(Eb, t1)
    end
    return Pb, Eb
end

"""
    _series_noise!(buf, J, Q, h, p)

`V(h) = int_0^h e^{J s} Q e^{J' s} ds` for the leading `p x p` blocks of `J`
and the symmetric `Q`, into the same block of `buf.out`, taping what
`_series_noise_pullback!` needs. `Q` is copied into `buf.Q` first.

At the scaled interval `h0`, `V = h0 W` with `W` the series
`sum_j h0^j / (j+1)! L^j(Q)`, `L(W) = J W + W J'`, by Horner's rule on the
operator: `W_{m+1} = Q`, `W_j = Q + (X W_{j+1} + W_{j+1} X') / (j+1)`. The
exponential the doubling needs, `e^X`, is its own Taylor series, taped the same
way in `buf.taylor`.
"""
function _series_noise!(buf::SeriesKernelBuffer{T}, J::AbstractMatrix,
    Q::AbstractMatrix, h, p::Int) where {T}
    _CTSEM_OPCOUNT.series_noise[] += 1
    _series_plan!(buf, J, h, p, 2.0)
    s, m, h0 = buf.s, buf.m, buf.h0
    X = _series_block(buf.X, p)
    Qb = _series_block(buf.Q, p)
    @inbounds for j in 1:p, i in 1:p
        Qb[i, j] = Q[i, j]
    end
    t1 = _series_block(buf.b1, p)
    for j in m:-1:1
        prev = j == m ? Qb : _series_block(buf.horner[j + 1], p)
        W = _series_block(buf.horner[j], p)
        _ctsem_mul!(t1, X, prev)
        c = inv(T(j + 1))
        @inbounds for jj in 1:p, ii in 1:p
            W[ii, jj] = Qb[ii, jj] + c * (t1[ii, jj] + t1[jj, ii])
        end
    end
    W1 = _series_block(buf.horner[1], p)
    out = _series_block(buf.out, p)
    first = s == 0 ? out : _series_block(buf.Ps[1], p)
    @inbounds for j in 1:p, i in 1:p
        first[i, j] = W1[i, j] * h0
    end
    if s > 0
        _series_grow!(buf.taylor, m, buf.n)
        Fm = _series_block(buf.taylor[m], p)
        @inbounds for j in 1:p, i in 1:p
            Fm[i, j] = X[i, j] / m
        end
        _series_add_identity!(Fm, p)
        for j in (m - 1):-1:1
            Fj = _series_block(buf.taylor[j], p)
            _ctsem_mul!(Fj, X, _series_block(buf.taylor[j + 1], p), inv(T(j)))
            _series_add_identity!(Fj, p)
        end
        copyto!(_series_block(buf.Es[1], p), _series_block(buf.taylor[1], p))
        for k in 1:s
            E = _series_block(buf.Es[k], p)
            P = _series_block(buf.Ps[k], p)
            next = k == s ? out : _series_block(buf.Ps[k + 1], p)
            _ctsem_mul!(t1, E, P)
            _ctsem_mulNT!(next, t1, E)
            next .+= P
            k < s && _ctsem_mul!(_series_block(buf.Es[k + 1], p), E, E)
        end
    end
    # Symmetric by construction up to rounding; made exactly so, since the
    # filter adds it to a covariance.
    @inbounds for j in 1:p, i in 1:(j - 1)
        v = (out[i, j] + out[j, i]) / 2
        out[i, j] = v
        out[j, i] = v
    end
    return out
end

"""
    _series_noise_pullback!(Jbar, Qbar, buf, p, Vbar)

Add to `Jbar` and `Qbar` the pullback of `_series_noise!` for cotangent `Vbar`,
which is symmetrised first. `buf` must hold that kernel's forward sweep,
unchanged since.
"""
function _series_noise_pullback!(Jbar::AbstractMatrix, Qbar::AbstractMatrix,
    buf::SeriesKernelBuffer{T}, p::Int, Vbar::AbstractMatrix) where {T}
    s, m, h0 = buf.s, buf.m, buf.h0
    X = _series_block(buf.X, p)
    Qb = _series_block(buf.Q, p)
    t1 = _series_block(buf.b2, p)
    Xb = _series_block(buf.Xbar, p)
    Vb, Eb = _series_noise_doubling_pullback!(buf, p, Vbar)
    fill!(Xb, zero(T))
    # E_0 = F_1, with F_{m+1} = I and F_j = I + X F_{j+1} / j.
    if s > 0
        Fb = Eb
        for j in 1:m
            c = inv(T(j))
            if j < m
                _ctsem_mulNT!(Xb, Fb, _series_block(buf.taylor[j + 1], p), c, one(T))
                _ctsem_mulTN!(t1, X, Fb, c)
                copyto!(Fb, t1)
            else
                @inbounds for jj in 1:p, ii in 1:p
                    Xb[ii, jj] += c * Fb[ii, jj]
                end
            end
        end
    end
    # V_0 = h0 W_1, and W_j = Q + (X W_{j+1} + W_{j+1} X') / (j+1), W_{m+1} = Q.
    Wb = Vb
    Wb .*= h0
    for j in 1:m
        c = inv(T(j + 1))
        next = j < m ? _series_block(buf.horner[j + 1], p) : Qb
        @inbounds for jj in 1:p, ii in 1:p
            Qbar[ii, jj] += Wb[ii, jj]
        end
        _ctsem_mul!(Xb, Wb, next, 2c, one(T))
        _ctsem_mulTN!(t1, X, Wb, c)
        _ctsem_mul!(t1, Wb, X, c, one(T))
        copyto!(Wb, t1)
    end
    @inbounds for jj in 1:p, ii in 1:p
        Qbar[ii, jj] += Wb[ii, jj]
    end
    @inbounds for j in 1:p, i in 1:p
        Jbar[i, j] += Xb[i, j] * h0
    end
    return Jbar, Qbar
end

"""
    _series_noise_doubling_pullback!(buf, p, Vbar)

The doubling chain of a noise series reversed, for cotangent `Vbar` on its
result (symmetrised first): the cotangents of the first scaled integral `V(h0)`
(in `buf.b1`) and of the first exponential `e^{J h0}` (in `buf.Ebar`, zero
without doublings). `buf` must hold the forward sweep, by either
`_series_noise!` or `_series_noise_table!`; `buf.b2` and `buf.b3` are its
scratch.
"""
function _series_noise_doubling_pullback!(buf::SeriesKernelBuffer{T}, p::Int,
    Vbar::AbstractMatrix) where {T}
    s = buf.s
    Vb = _series_block(buf.b1, p)
    Eb = _series_block(buf.Ebar, p)
    t1 = _series_block(buf.b2, p)
    t2 = _series_block(buf.b3, p)
    @inbounds for j in 1:p, i in 1:p
        Vb[i, j] = (Vbar[i, j] + Vbar[j, i]) / 2
    end
    fill!(Eb, zero(T))
    # V_k = V_{k-1} + E V_{k-1} E', and E_k = E^2 while k < s, E = E_{k-1}.
    for k in s:-1:1
        E = _series_block(buf.Es[k], p)
        P = _series_block(buf.Ps[k], p)
        if k < s
            _ctsem_mulNT!(t1, Eb, E)
            _ctsem_mulTN!(t1, E, Eb, one(T), one(T))
        else
            fill!(t1, zero(T))
        end
        _ctsem_mul!(t2, Vb, E)                     # Vbar E
        _ctsem_mul!(t1, t2, P, T(2), one(T))       # + 2 Vbar E V
        copyto!(Eb, t1)
        _ctsem_mulTN!(t1, E, t2)                   # E' Vbar E
        Vb .+= t1
    end
    return Vb, Eb
end

################################################################################
# The series from powers built once per drift
################################################################################
#
# The kernels above sum their Taylor series by Horner's rule at every
# interval: about 2m products of the block for the noise integral and its
# exponential, and m for the intercept, with m 15-18 at the scaled step, and
# the reverse pass sums them again before pulling back through every term.
# Where the drift and diffusion do not change from row to row -- every linear
# model, within a subject's filter pass -- the series' matrices do not change
# either; only an interval's scalar weights on them do. So a `SeriesTable`
# holds the powers `J^j` and, for the noise, `L^j(Q)` with
# `L(W) = J W + W J'`, built once for its `(J, Q)` and extended when a longer
# series is asked for, and an interval's series is their weighted sum: O(m p^2)
# where Horner's rule was O(m p^3). Only the doublings stay products. The
# truncation is `_series_plan!`'s, the same as the kernels above, so the two
# agree to rounding.
#
# The reverse pass adds each interval's cotangents onto the powers
# (`SeriesTableBar`) and pulls them through the power chains once for the run
# of intervals that share a table (`_series_table_pullback!`, released by
# `_flush_series!` in adjoint_ekf.jl as the Fréchet and diffusion cotangents
# are). A drift that changes every row, as a state-dependent one does, builds
# and releases a table each row, which costs about what Horner's rule did.
#
# On the affect pilot -- a slow tier whose small drift puts almost every
# interval on the series route -- the series was 56% of a Laplace value and
# gradient (8 threads, local, 2026-10-25). Per interval, a value and gradient
# through the tables cost 2.3 times less than Horner's on a 10-state block
# with 30 intervals to a drift (scripts/pilot-curvature-2026-10/
# series_bench.jl, local, one thread); the doublings are most of what is left.
# A longer scaled step trades doublings for terms, which the tables make
# cheap (3.0 times at a step of norm 1), but it costs digits where the series
# matters most: with the series forced, the explosive drift over a 43-day gap
# in test_adjoint_gradient_validation.jl read 1.2e-5 off BigFloat, over the
# test's 1e-5. So the tables keep Horner's step.

"""
    SeriesTable{T}(n, noise)

Powers of the leading `p x p` block of a drift `J`, for the series kernels:
`Jp[j] = J^j` for `j = 1:count`, and for a noise table also
`Lp[j + 1] = L^j(Q)` for `j = 0:count`. `J` and `Q` are copies of the block
the powers belong to; `_series_table_key!` compares against them.
"""
mutable struct SeriesTable{T}
    p::Int
    noise::Bool
    valid::Bool
    count::Int
    J::Matrix{T}
    Q::Matrix{T}
    Jp::Vector{Matrix{T}}
    Lp::Vector{Matrix{T}}
    tmp::Matrix{T}
end

SeriesTable(::Type{T}, n::Int, noise::Bool) where {T} =
    SeriesTable{T}(0, noise, false, 0, zeros(T, n, n), zeros(T, n, n),
        Matrix{T}[], Matrix{T}[], zeros(T, n, n))

"""
    SeriesTableBar{T}(n)

Cotangents on a `SeriesTable`'s powers, laid out as its `Jp` and `Lp`,
accumulated over the intervals that read the table; `count` is the highest
power written and `pending` whether anything has been.
"""
mutable struct SeriesTableBar{T}
    pending::Bool
    count::Int
    Jp::Vector{Matrix{T}}
    Lp::Vector{Matrix{T}}
end

SeriesTableBar(::Type{T}) where {T} = SeriesTableBar{T}(false, 0, Matrix{T}[],
    Matrix{T}[])

"""
    _series_table_key!(tab, J, Q, p)

Whether `tab` holds the powers of the leading `p x p` block of `J` (and of
`Q`, for a noise table). If not, it is reset to that block with no powers
built, and the answer is `false`. The comparison is exact, and under
`ForwardDiff.Dual` it compares partials too (`_blocks_identical`).
"""
function _series_table_key!(tab::SeriesTable, J::AbstractMatrix, Q, p::Int)
    _series_table_matches(tab, J, Q, p) && return true
    _series_table_reset!(tab, J, Q, p)
    return false
end

"""Whether `tab` holds the powers of the leading `p x p` block of `J` (and `Q`)."""
@inline _series_table_matches(tab::SeriesTable, J::AbstractMatrix, Q, p::Int) =
    tab.valid && tab.p == p && _blocks_identical(tab.J, J, p, p) &&
        (!tab.noise || _blocks_identical(tab.Q, Q, p, p))

"""Reset `tab` to the leading `p x p` block of `J` (and `Q`), no powers built."""
function _series_table_reset!(tab::SeriesTable, J::AbstractMatrix, Q, p::Int)
    _CTSEM_OPCOUNT.series_tables[] += 1
    _copy_block!(tab.J, J, p, p)
    tab.noise && _copy_block!(tab.Q, Q, p, p)
    tab.p = p
    tab.count = 0
    tab.valid = true
    return tab
end

"""
    CTSEMSeriesReverse{T}(n, k)

The reverse pass's series tables (`SeriesTable`), one for the noise over the
`k` diffusing states and one for the intercept over the affine block, with the
cotangents waiting on their powers, the raw DIFFUSION matrix the noise table's
intervals defer their diffusion cotangent against, and scratch for releasing
them (`_flush_series!` in adjoint_ekf.jl).
"""
struct CTSEMSeriesReverse{T}
    noise::SeriesTable{T}
    noise_bar::SeriesTableBar{T}
    intercept::SeriesTable{T}
    intercept_bar::SeriesTableBar{T}
    diffusion::Matrix{T}
    Jbar::Matrix{T}
    Qbar::Matrix{T}
    Qfull::Matrix{T}
end

CTSEMSeriesReverse(::Type{T}, n::Int, k::Int) where {T} = CTSEMSeriesReverse{T}(
    SeriesTable(T, k, true), SeriesTableBar(T), SeriesTable(T, n, false),
    SeriesTableBar(T), zeros(T, n, n), zeros(T, n, n), zeros(T, n, n), zeros(T, n, n))

"""Extend `tab`'s powers to `count`, building only those it lacks."""
function _series_table_extend!(tab::SeriesTable{T}, count::Int) where {T}
    count <= tab.count && return tab
    p = tab.p
    n = size(tab.J, 1)
    _series_grow!(tab.Jp, count, n)
    tab.noise && _series_grow!(tab.Lp, count + 1, n)
    J = _series_block(tab.J, p)
    t = _series_block(tab.tmp, p)
    if tab.count == 0 && tab.noise
        _copy_block!(tab.Lp[1], tab.Q, p, p)
    end
    for j in (tab.count + 1):count
        Jj = _series_block(tab.Jp[j], p)
        if j == 1
            copyto!(Jj, J)
        else
            _ctsem_mul!(Jj, J, _series_block(tab.Jp[j - 1], p))
        end
        if tab.noise
            # L^j(Q) = J L^{j-1}(Q) + its transpose, symmetric as Q is.
            _ctsem_mul!(t, J, _series_block(tab.Lp[j], p))
            Lj = _series_block(tab.Lp[j + 1], p)
            @inbounds for jj in 1:p, ii in 1:p
                Lj[ii, jj] = t[ii, jj] + t[jj, ii]
            end
        end
    end
    tab.count = count
    return tab
end

"""`A .+= c .* B` over the leading `p x p` block."""
@inline function _series_axpy!(A::AbstractMatrix, c, B::AbstractMatrix, p::Int)
    @inbounds for j in 1:p, i in 1:p
        A[i, j] += c * B[i, j]
    end
    return A
end

"""
    _series_doublings!(buf, out, p, noise)

The doubling chain both table kernels end with, from `buf.Es[1]` and
`buf.Ps[1]` to `out`: `P_{k+1} = P_k + E_k P_k E_k'` for the noise,
`P_{k+1} = P_k + E_k P_k` for the intercept, `E_{k+1} = E_k^2`; the same
recurrences, in the same order, as `_series_noise!` and `_series_intercept!`.
"""
function _series_doublings!(buf::SeriesKernelBuffer{T}, out::AbstractMatrix, p::Int,
    noise::Bool) where {T}
    s = buf.s
    t1 = _series_block(buf.b1, p)
    for k in 1:s
        E = _series_block(buf.Es[k], p)
        P = _series_block(buf.Ps[k], p)
        next = k == s ? out : _series_block(buf.Ps[k + 1], p)
        if noise
            _ctsem_mul!(t1, E, P)
            _ctsem_mulNT!(next, t1, E)
        else
            _ctsem_mul!(next, E, P)
        end
        next .+= P
        k < s && _ctsem_mul!(_series_block(buf.Es[k + 1], p), E, E)
    end
    return out
end

"""
    _series_noise_table!(buf, tab, h, p)

`_series_noise!` for the `J` and `Q` that `tab` holds, from its powers:
`V(h0) = sum_{j=0}^m h0^{j+1} / (j+1)! L^j(Q)` and
`e^{J h0} = sum_{j=0}^m h0^j / j! J^j`, then the same doublings. Tapes what
`_series_noise_table_pullback!` needs in `buf`.
"""
function _series_noise_table!(buf::SeriesKernelBuffer{T}, tab::SeriesTable{T}, h,
    p::Int) where {T}
    _CTSEM_OPCOUNT.series_noise[] += 1
    _series_plan!(buf, tab.J, h, p, 2.0)
    s, m, h0 = buf.s, buf.m, buf.h0
    _series_table_extend!(tab, m)
    out = _series_block(buf.out, p)
    first = s == 0 ? out : _series_block(buf.Ps[1], p)
    c = h0
    _copy_block!(first, tab.Lp[1], p, p)
    first .*= c
    for j in 1:m
        c *= h0 / (j + 1)
        _series_axpy!(first, c, tab.Lp[j + 1], p)
    end
    if s > 0
        E0 = _series_block(buf.Es[1], p)
        fill!(E0, zero(T))
        _series_add_identity!(E0, p)
        e = 1.0
        for j in 1:m
            e *= h0 / j
            _series_axpy!(E0, e, tab.Jp[j], p)
        end
        _series_doublings!(buf, out, p, true)
    end
    @inbounds for j in 1:p, i in 1:(j - 1)
        v = (out[i, j] + out[j, i]) / 2
        out[i, j] = v
        out[j, i] = v
    end
    return out
end

"""
    _series_intercept_table!(buf, tab, h, p)

`_series_intercept!` for the `J` that `tab` holds, from its powers:
`Phi(h0) = sum_{j=0}^m h0^{j+1} / (j+1)! J^j` and
`e^{J h0} = sum_{i=0}^{m+1} h0^i / i! J^i`, then the same doublings. Tapes what
`_series_intercept_table_pullback!` needs in `buf`.
"""
function _series_intercept_table!(buf::SeriesKernelBuffer{T}, tab::SeriesTable{T}, h,
    p::Int) where {T}
    _CTSEM_OPCOUNT.series_intercept[] += 1
    _series_plan!(buf, tab.J, h, p, 1.0)
    s, m, h0 = buf.s, buf.m, buf.h0
    _series_table_extend!(tab, s > 0 ? m + 1 : m)
    out = _series_block(buf.out, p)
    first = s == 0 ? out : _series_block(buf.Ps[1], p)
    c = h0
    fill!(first, zero(T))
    @inbounds for i in 1:p
        first[i, i] = c
    end
    for j in 1:m
        c *= h0 / (j + 1)
        _series_axpy!(first, c, tab.Jp[j], p)
    end
    s == 0 && return out
    E0 = _series_block(buf.Es[1], p)
    fill!(E0, zero(T))
    _series_add_identity!(E0, p)
    e = 1.0
    for i in 1:(m + 1)
        e *= h0 / i
        _series_axpy!(E0, e, tab.Jp[i], p)
    end
    return _series_doublings!(buf, out, p, false)
end

"""Grow `bar` to hold cotangents for powers up to `count` of `tab`'s shape."""
function _series_bar_grow!(bar::SeriesTableBar{T}, tab::SeriesTable{T},
    count::Int) where {T}
    n = size(tab.J, 1)
    p = tab.p
    while length(bar.Jp) < count
        push!(bar.Jp, zeros(T, n, n))
    end
    if tab.noise
        while length(bar.Lp) < count + 1
            push!(bar.Lp, zeros(T, n, n))
        end
    end
    # Nothing pending is a fresh start, whatever an earlier walk left; past
    # that, the powers not yet written are zeroed as they are first reached.
    if !bar.pending
        bar.count = 0
        tab.noise && fill!(_series_block(bar.Lp[1], p), zero(T))
    end
    for j in (bar.count + 1):count
        fill!(_series_block(bar.Jp[j], p), zero(T))
        tab.noise && fill!(_series_block(bar.Lp[j + 1], p), zero(T))
    end
    bar.count = max(bar.count, count)
    bar.pending = true
    return bar
end

"""
    _series_noise_table_pullback!(bar, buf, tab, p, Vbar)

Add to `bar` the cotangents, on `tab`'s powers, of the interval
`_series_noise_table!` last computed in `buf`, for cotangent `Vbar` on it.
"""
function _series_noise_table_pullback!(bar::SeriesTableBar{T}, buf::SeriesKernelBuffer{T},
    tab::SeriesTable{T}, p::Int, Vbar::AbstractMatrix) where {T}
    s, m, h0 = buf.s, buf.m, buf.h0
    Vb, Eb = _series_noise_doubling_pullback!(buf, p, Vbar)
    _series_bar_grow!(bar, tab, m)
    c = h0
    _series_axpy!(bar.Lp[1], c, Vb, p)
    for j in 1:m
        c *= h0 / (j + 1)
        _series_axpy!(bar.Lp[j + 1], c, Vb, p)
    end
    if s > 0
        e = 1.0
        for j in 1:m
            e *= h0 / j
            _series_axpy!(bar.Jp[j], e, Eb, p)
        end
    end
    return bar
end

"""
    _series_intercept_table_pullback!(bar, buf, tab, p, Phibar)

Add to `bar` the cotangents, on `tab`'s powers, of the interval
`_series_intercept_table!` last computed in `buf`, for cotangent `Phibar`.
"""
function _series_intercept_table_pullback!(bar::SeriesTableBar{T},
    buf::SeriesKernelBuffer{T}, tab::SeriesTable{T}, p::Int,
    Phibar::AbstractMatrix) where {T}
    s, m, h0 = buf.s, buf.m, buf.h0
    Pb, Eb = _series_intercept_doubling_pullback!(buf, p, Phibar)
    _series_bar_grow!(bar, tab, s > 0 ? m + 1 : m)
    c = h0
    for j in 1:m
        c *= h0 / (j + 1)
        _series_axpy!(bar.Jp[j], c, Pb, p)
    end
    if s > 0
        e = 1.0
        for i in 1:(m + 1)
            e *= h0 / i
            _series_axpy!(bar.Jp[i], e, Eb, p)
        end
    end
    return bar
end

"""
    _series_table_pullback!(Jbar, Qbar, tab, bar)

Pull `bar`'s cotangents on `tab`'s powers back through the power chains,
adding to the leading `p x p` blocks of `Jbar` (and `Qbar`, for a noise
table), and release `bar`. `J^j = J J^{j-1}`;
`L^j = J L^{j-1} + L^{j-1} J'` with `L^0 = Q`.
"""
function _series_table_pullback!(Jbar::AbstractMatrix, Qbar, tab::SeriesTable{T},
    bar::SeriesTableBar{T}) where {T}
    bar.pending || return Jbar
    p = tab.p
    c = bar.count
    J = _series_block(tab.J, p)
    Jb = _series_block(Jbar, p)
    if tab.noise
        for j in c:-1:1
            Lb = _series_block(bar.Lp[j + 1], p)
            A = _series_block(tab.Lp[j], p)
            prev = _series_block(bar.Lp[j], p)
            _ctsem_mulNT!(Jb, Lb, A, one(T), one(T))
            _ctsem_mulTN!(Jb, Lb, A, one(T), one(T))
            _ctsem_mulTN!(prev, J, Lb, one(T), one(T))
            _ctsem_mul!(prev, Lb, J, one(T), one(T))
        end
        _series_axpy!(Qbar, one(T), bar.Lp[1], p)
    end
    for j in c:-1:2
        Pjb = _series_block(bar.Jp[j], p)
        _ctsem_mulNT!(Jb, Pjb, _series_block(tab.Jp[j - 1], p), one(T), one(T))
        _ctsem_mulTN!(_series_block(bar.Jp[j - 1], p), J, Pjb, one(T), one(T))
    end
    c >= 1 && _series_axpy!(Jb, one(T), bar.Jp[1], p)
    bar.pending = false
    bar.count = 0
    return Jbar
end

"""The smallest diagonal magnitude of the leading `p x p` block, on primal values."""
@inline function _series_min_pivot(LU::AbstractMatrix, p::Int)
    smallest = Inf
    @inbounds for i in 1:p
        smallest = min(smallest, abs(Float64(_primal(LU[i, i]))))
    end
    return smallest
end

"""
The conditioning of the Lyapunov operator a buffer last solved against: the
smallest pivot of the packed system's LU, or, on the Schur route, the smallest
`|s_ii + s_jj|` over the Schur factor's diagonal -- the eigenvalue sums the
operator divides by, a 2x2 block contributing its real part twice.
"""
_series_lyap_conditioning(buffer::LyapKsolveBuffer{T}) where {T} =
    _series_min_pivot(buffer.ksolve_O, buffer.ksolve_system_dim)

function _series_lyap_conditioning(buffer::LyapSchurBuffer)
    S = buffer.S
    k = size(S, 1)
    smallest = Inf
    @inbounds for j in 1:k, i in 1:j
        smallest = min(smallest, abs(Float64(S[i, i] + S[j, j])))
    end
    return smallest
end

"""Whether a pivot `conditioning` is small against the interval `h`."""
@inline _series_needed(conditioning::Float64, h) =
    !(conditioning * Float64(_primal(h)) >= _CTSEM_SERIES_BELOW[])
