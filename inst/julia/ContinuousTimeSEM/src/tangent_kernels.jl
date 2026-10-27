# The hand-coded second order's kernels, on `tangent_primitives.jl`: each the
# engine's own computation, the same operations in the same order, each paired
# with its tangent, in Float64 with the W directions innermost -- what
# `ctsem_hessian`'s dual route computes for it, against which
# test_tangent_hessian.jl checks them.
#
# The Gaussian measurement update's reverse (`_reverse_update!`, a row whose
# update is not relinearised and integrates no deviations).

struct TanUpdate{W}
    n::Int
    m::Int
    # The record and the cotangents coming in.
    P_in::TanMat{W}; H::TanMat{W}; Lambda::TanMat{W}; R::TanMat{W}
    x::TanVec{W}; mu::TanVec{W}; y::TanVec{W}
    xbar::TanVec{W}; Pbar::TanMat{W}
    # Out: the parameter-side cotangents.
    Rbar::TanMat{W}; Hbar::TanMat{W}; Lbar::TanMat{W}; ybarbar::TanVec{W}
    # Intermediates.
    Pr::TanMat{W}; PHt::TanMat{W}; mm1::TanMat{W}; S::TanMat{W}; Ui::TanMat{W}; Xc::TanMat{W}; Tc::TanMat{W}
    yt::TanVec{W}; e::TanVec{W}; a::TanVec{W}
    Kb::TanMat{W}; K::TanMat{W}; Mm::TanMat{W}; Wh::TanMat{W}; nn1::TanMat{W}; Pp::TanMat{W}; nm1::TanMat{W}
    Ps::TanMat{W}
    Kx::TanVec{W}; We::TanVec{W}; HKx::TanVec{W}; Ky::TanVec{W}; Ppx::TanVec{W}; gv::TanVec{W}; cv::TanVec{W}
    Pbar_new::TanMat{W}; Kc::TanMat{W}; Mc::TanMat{W}; Z::TanMat{W}; Y::TanMat{W}; KcS::TanMat{W}
    mm2::TanMat{W}; S0::TanMat{W}
    xbar_new::TanVec{W}
    F::Matrix{Float64}
end

function TanUpdate(n::Int, m::Int, W::Int)
    T = TanMat{W}; V = TanVec{W}
    nn() = T(n, n); nm() = T(n, m); mn() = T(m, n); mm() = T(m, m)
    return TanUpdate{W}(n, m,
        nn(), mn(), mn(), mm(), V(n), V(m), V(m), V(n), nn(),
        mm(), mn(), mn(), V(m),
        nn(), nm(), mm(), mm(), mm(), mm(), mm(),
        V(m), V(m), V(m),
        nm(), nm(), nn(), mn(), nn(), nn(), nm(), nn(),
        V(m), V(n), V(n), V(n), V(n), V(n), V(n),
        nn(), nm(), nn(), mn(), mm(), nm(), mm(), mm(),
        V(n),
        zeros(m, m))
end

"""
The reverse update with tangents, reading the record and the incoming `xbar`,
`Pbar` from `u`, and leaving the outgoing `xbar`, `Pbar` there and the
parameter-side cotangents in `u.Rbar`, `u.Hbar`, `u.Lbar` and `u.ybarbar`
(MANIFESTMEANS' is its negative). `false` where `S` does not factor.
"""
function _tan_reverse_update!(u::TanUpdate{W}) where {W}
    n, m = u.n, u.m
    # Pr = P_in + ridge, PHt = Pr H', S = sym(H PHt + R) + ridge.
    _tan_copy!(u.Pr, u.P_in)
    @inbounds for i in 1:n; u.Pr.v[i, i] += _CTSEM_RIDGE; end
    _tan_mulNT!(u.PHt, u.Pr, u.H)
    _tan_copy!(u.mm1, u.R)
    _tan_mul!(u.mm1, u.H, u.PHt, true, true, Val(:full), Val(true))
    _tan_sym!(u.S, u.mm1)
    @inbounds for i in 1:m; u.S.v[i, i] += _CTSEM_RIDGE; end
    _tan_cholinv!(u.Ui, u.S, u.Xc, u.Tc, u.F) || return false
    Ui = u.Ui
    # yt = y - (Lambda x + mu).
    _tan_copy!(u.yt, u.mu)
    _tan_mulvec!(u.yt, u.Lambda, u.x, true, true)
    _tan_lincomb!(u.yt, true, u.y, -1.0, u.yt)
    e, a, Kb, K, M, Wh = u.e, u.a, u.Kb, u.K, u.Mm, u.Wh
    _tan_mulTvec!(e, Ui, u.yt, true, false, Val(:AupperT))  # e = Ui' yt
    _tan_mulvec!(a, Ui, e, true, false, Val(:Aupper))  # a = Ui e
    _tan_mul!(Kb, u.PHt, Ui, true, false, Val(:Bupper))  # Kb = PHt Ui
    _tan_mulNT!(K, Kb, Ui, true, false, Val(:BupperT))  # K = Kb Ui'
    _tan_mul!(M, K, u.H, -1.0, false)  # M = I - K H
    @inbounds for i in 1:n; M.v[i, i] += 1.0; end
    _tan_mulTN!(Wh, Ui, u.H, true, false, Val(:AupperT))  # Wh = Ui' H
    _tan_mul!(u.nn1, M, u.Pr)  # Pp = M Pr M' + K R K' + eps K K'
    _tan_mulNT!(u.Pp, u.nn1, M, true, false, Val(:full), Val(true))
    _tan_mul!(u.nm1, K, u.R)
    _tan_mulNT!(u.Pp, u.nm1, K, true, true, Val(:full), Val(true))
    _tan_mulNT!(u.Pp, K, K, _CTSEM_RIDGE, true, Val(:full), Val(true))
    _tan_sym!(u.Ps, u.Pbar)
    _tan_mulTvec!(u.Kx, K, u.xbar)
    _tan_mulTvec!(u.We, Wh, e)
    _tan_mulTvec!(u.HKx, u.H, u.Kx)
    _tan_mulvec!(u.Ky, Kb, e)
    _tan_mulvec!(u.Ppx, u.Pp, u.xbar)

    # Pbar of the prediction.
    Pn = u.Pbar_new
    _tan_mulTN!(Pn, Wh, Wh, -0.5, false, Val(:full), Val(true))
    # + 1/2 (We We' + sym(xbar We') - sym(HKx We')) = 1/2 (We c' + g We'),
    # g = xbar - HKx, c = We + g.
    We, xb, HKx, g, c = u.We, u.xbar, u.HKx, u.gv, u.cv
    _tan_lincomb!(g, true, xb, -1.0, HKx)
    _tan_lincomb!(c, true, g, true, We)
    _tan_ger!(Pn, 0.5, We, c)
    _tan_ger!(Pn, 0.5, g, We)
    _tan_mul!(u.nn1, u.Ps, M)
    _tan_mulTN!(Pn, M, u.nn1, true, true, Val(:full), Val(true))  # + M' Ps M

    # The ridge correction.
    eps = _CTSEM_RIDGE
    Kc, Mc = u.Kc, u.Mc
    _tan_mul!(Kc, u.Ps, K, -2eps, false)
    _tan_mul!(Mc, u.Ps, M, -2eps, false)
    _tan_mulNT!(Kc, Mc, u.H, -1.0, true)
    _tan_mul!(u.nm1, Kc, Ui, true, false, Val(:Bupper))
    _tan_mul!(u.nn1, u.nm1, Wh)
    _tan_add!(Pn, 0.5, u.nn1)
    _tan_add!(Pn, 0.5, u.nn1, true, Val(:T))

    # Hbar.
    Z, Ky, Hbar = u.Z, u.Ky, u.Hbar
    _tan_add!(Z, true, Kb, false, Val(:T))  # Z = Kb' - e Ky'
    _tan_ger!(Z, -1.0, e, Ky)
    _tan_mul!(Hbar, Ui, Z, -1.0, false, Val(:Aupper))
    Ppx, Kx = u.Ppx, u.Kx
    _tan_ger!(Hbar, true, a, Ppx)
    _tan_ger!(Hbar, -1.0, Kx, Ky)
    _tan_mul!(u.nn1, u.Ps, u.Pp)
    _tan_mulTN!(Hbar, K, u.nn1, -2.0, true)

    # Rbar.
    Y, Rbar = u.Y, u.Rbar
    _tan_add!(Y, true, Ui, false, Val(:T))  # Y = Ui' - e a'
    _tan_ger!(Y, -1.0, e, a)
    _tan_mul!(Rbar, Ui, Y, -0.5, false, Val(:Aupper), Val(true))
    _tan_ger!(Rbar, -0.5, Kx, a)
    _tan_ger!(Rbar, -0.5, a, Kx)
    _tan_mul!(u.nm1, u.Ps, K)
    _tan_mulTN!(Rbar, K, u.nm1, true, true, Val(:full), Val(true))

    # The ridge correction's Hbar, Rbar and Pbar through S.
    _tan_mulTN!(Hbar, K, Mc, -1.0, true)
    _tan_mul!(u.nm1, Kc, Ui, true, false, Val(:Bupper))
    _tan_mulNT!(u.KcS, u.nm1, Ui, true, false, Val(:BupperT))  # Kc S^-1
    _tan_mul!(u.nm1, u.Pr, u.KcS)
    _tan_add!(Hbar, true, u.nm1, true, Val(:T))
    _tan_mulTN!(u.mm2, K, u.KcS, -1.0, false)
    _tan_sym!(u.S0, u.mm2)
    Lbar = u.Lbar
    _tan_mul!(Lbar, u.S0, u.H)
    _tan_mulTN!(Pn, u.H, Lbar, true, true, Val(:full), Val(true))
    _tan_mul!(Hbar, Lbar, u.Pr, 2.0, true)
    _tan_add!(Rbar, true, u.S0)

    # yt through x+ = x + K yt and the log likelihood.
    yb = u.ybarbar
    _tan_lincomb!(yb, true, Kx, -1.0, a)
    _tan_copy!(u.xbar_new, xb)
    _tan_ger!(Lbar, -1.0, yb, u.x, false)
    _tan_mulTvec!(u.xbar_new, u.Lambda, yb, -1.0, true)
    _tan_copy!(u.xbar, u.xbar_new)
    _tan_copy!(u.Pbar, Pn)
    return true
end

# The exponential and its Fréchet derivative (`my_exp_frechet!`, Al-Mohy &
# Higham 2009), for the flush of the queued Fréchet directions -- the reverse
# pass's largest cost -- except that the Padé denominator is factored once, its
# values only, where the dual factors it twice with the directions riding along.

struct TanFrechet{W}
    n::Int
    As::TanMat{W}; Es::TanMat{W}; A2::TanMat{W}; A4::TanMat{W}; A6::TanMat{W}
    M2::TanMat{W}; M4::TanMat{W}; M6::TanMat{W}
    W1::TanMat{W}; Wm::TanMat{W}; Z1::TanMat{W}; V::TanMat{W}; U::TanMat{W}
    Lw::TanMat{W}; Lu::TanMat{W}; Lv::TanMat{W}; T1::TanMat{W}; T2::TanMat{W}
    Y::TanMat{W}; L::TanMat{W}
    F::Matrix{Float64}
    piv::Vector{Int}
end
TanFrechet(n::Int, W::Int) =
    TanFrechet{W}(n, (TanMat{W}(n, n) for _ in 1:20)..., zeros(n, n), zeros(Int, n))

"""
`f.Y = exp(A)` and `f.L = L(A, E)`, the Fréchet derivative of the exponential
at `A` in direction `E`, with their tangents from `A`'s and `E`'s.
"""
function _tan_exp_frechet!(f::TanFrechet{W}, A::TanMat{W}, E::TanMat{W}) where {W}
    n = f.n
    As, Es, Y, L = f.As, f.Es, f.Y, f.L
    _tan_copy!(As, A)
    _tan_copy!(Es, E)
    a1 = _opnorm1_noalloc(As.v, n)
    e1 = _opnorm1_noalloc(Es.v, n)
    if !(isfinite(a1) && isfinite(e1))
        for X in (Y, L)
            fill!(X.v, NaN)
            fill!(X.d, NaN)
            X.nz[] = true
        end
        return f
    end
    theta = _FRECHET_THETA13
    s = a1 <= theta ? 0 : ceil(Int, log2(a1 / theta))
    factor = inv(2.0^s)
    _tan_scale!(As, factor)
    _tan_scale!(Es, factor)

    A2, A4, A6, M2, M4, M6 = f.A2, f.A4, f.A6, f.M2, f.M4, f.M6
    _tan_mul!(A2, As, As)
    _tan_mul!(M2, As, Es)
    _tan_mul!(M2, Es, As, true, true)
    _tan_mul!(A4, A2, A2)
    _tan_mul!(M4, A2, M2)
    _tan_mul!(M4, M2, A2, true, true)
    _tan_mul!(A6, A4, A2)
    _tan_mul!(M6, A4, M2)
    _tan_mul!(M6, M4, A2, true, true)

    b = _pade13_b(Float64)
    W1, Wm, Z1, V, U = f.W1, f.Wm, f.Z1, f.V, f.U
    Lw, Lu, Lv, T1, T2 = f.Lw, f.Lu, f.Lv, f.T1, f.T2

    # Odd part U = A (A6 W1 + W2), even part V = A6 Z1 + Z2.
    _tan_comb!(W1, ((b[14], A6), (b[12], A4), (b[10], A2)))
    _tan_comb!(Z1, ((b[13], A6), (b[11], A4), (b[9], A2)))
    _tan_mul!(Wm, A6, W1)
    _tan_comb!(Wm, ((b[8], A6), (b[6], A4), (b[4], A2)), true)
    @inbounds for i in 1:n; Wm.v[i, i] += b[2]; end
    _tan_mul!(V, A6, Z1)
    _tan_comb!(V, ((b[7], A6), (b[5], A4), (b[3], A2)), true)
    @inbounds for i in 1:n; V.v[i, i] += b[1]; end
    _tan_mul!(U, As, Wm)

    # Their Fréchet derivatives, term by term.
    _tan_comb!(T1, ((b[14], M6), (b[12], M4), (b[10], M2)))  # L(W1)
    _tan_mul!(Lw, A6, T1)
    _tan_mul!(Lw, M6, W1, true, true)
    _tan_comb!(Lw, ((b[8], M6), (b[6], M4), (b[4], M2)), true)  # L(W)
    _tan_mul!(Lu, As, Lw)
    _tan_mul!(Lu, Es, Wm, true, true)  # L(U) = A L(W) + E W
    _tan_comb!(T1, ((b[13], M6), (b[11], M4), (b[9], M2)))  # L(Z1)
    _tan_mul!(Lv, A6, T1)
    _tan_mul!(Lv, M6, Z1, true, true)
    _tan_comb!(Lv, ((b[7], M6), (b[5], M4), (b[3], M2)), true)  # L(V)

    # R = (V - U)^-1 (V + U), and L = (V - U)^-1 (LU + LV + (LU - LV) R): one
    # factorization of the values of V - U serves both solves.
    _tan_comb!(T1, ((1.0, V), (-1.0, U)))
    _tan_comb!(Y, ((1.0, V), (1.0, U)))
    F, piv = f.F, f.piv
    copyto!(F, T1.v)
    _lu_factor_generic!(F, piv, n)
    _tan_solve!(Y, F, piv, T1)
    _tan_comb!(T2, ((1.0, Lu), (-1.0, Lv)))
    _tan_mul!(L, T2, Y)
    _tan_comb!(L, ((1.0, Lu), (1.0, Lv)), true)
    _tan_solve!(L, F, piv, T1)

    # Undo the scaling: L(A^2, E) = A L + L A at each doubling.
    for _ in 1:s
        _tan_mul!(T1, Y, L)
        _tan_mul!(T1, L, Y, true, true)
        _tan_copy!(L, T1)
        _tan_mul!(T1, Y, Y)
        _tan_copy!(Y, T1)
    end
    return f
end


# The Lyapunov solve (`my_lyap!`, `my_lyap_transposed!`) and its pullback
# (`_ctsem_lyap_pullback!`).
# The values are the engine's Float64 solve, on its cached Schur factors
# `A = U S U'`; the tangents solve, differentiating `A X + X A' + Q = 0`,
#
#     A dX + dX A' + (dA X + X dA' + dQ) = 0
#
# (transposed: `dA' X + X dA`) on the same factors, as the engine's dual route
# does -- but all W directions at once, as tuples through one pass of the
# transforms and the quasi-triangular back substitution, each small Kronecker
# block factored once where the dual route factors it once per direction.

struct TanLyap{W}
    n::Int
    buf::LyapBuffer{Float64,Nothing}
    Ut::TanMat{W}  # U, values only
    Xv::TanMat{W}  # the solution's values, values only
    R::TanMat{W}; tmp::TanMat{W}; C::TanMat{W}; Y::TanMat{W}; rhs::TanMat{W}
    M::Matrix{Float64}  # a block's Kronecker system, 4 x 4
end
function TanLyap(n::Int, W::Int)
    Ut, Xv = TanMat{W}(n, n), TanMat{W}(n, n)
    _tan_setnz!(Ut, false)
    _tan_setnz!(Xv, false)
    return TanLyap{W}(n, LyapBuffer(Float64, n), Ut, Xv,
        (TanMat{W}(n, n) for _ in 1:5)..., zeros(4, 4))
end

"""
`A X + X A' + Q = 0` (`transposed`: `A' X + X A + Q = 0`) for `X`, with
tangents from `A`'s and `Q`'s. `Q` symmetric, as the engine's solve requires.
"""
function _tan_lyap!(X::TanMat{W}, A::TanMat{W}, Q::TanMat{W}, h::TanLyap{W}, transposed::Bool) where {W}
    n, buf = h.n, h.buf
    _lyap_solve!(X.v, A.v, Q.v, buf, transposed)
    (_tan_nz(A) || _tan_nz(Q)) || return _tan_setnz!(X, false)
    if !buf.converged
        fill!(X.d, NaN)
        X.nz[] = true
        return X
    end
    Ut, Xv, R, tmp, C, Y = h.Ut, h.Xv, h.R, h.tmp, h.C, h.Y
    copyto!(Ut.v, buf.U)
    copyto!(Xv.v, X.v)
    # R = dA X + X dA' + dQ (transposed: dA' X + X dA + dQ), X symmetric.
    transposed ? _tan_mulTN!(tmp, A, Xv) : _tan_mul!(tmp, A, Xv)
    _tan_add!(R, true, tmp, false)
    _tan_add!(R, true, tmp, true, Val(:T))
    _tan_add!(R, true, Q, true)
    # C = -U' R U; the quasi-triangular solve; Y = U C U', symmetrised.
    _tan_mulTN!(tmp, Ut, R)
    _tan_mul!(C, tmp, Ut, -1.0, false)
    _tan_quasitri_lyap!(C.d, buf.S, n, h.M, transposed, Val(W))
    _tan_mul!(tmp, Ut, C)
    _tan_mulNT!(Y, tmp, Ut)
    @inbounds for j in 1:n, i in 1:j
        t = _tan_axpby(0.5, _tan_col(Y, i, j), 0.5, _tan_col(Y, j, i))
        _tan_put!(X, i, j, t)
        _tan_put!(X, j, i, t)
    end
    X.nz[] = true
    return X
end

"""
`(Abar, Wt)` for `X = lyap(A, Q)` with tangents: `Wt` solves
`A' W + W A + sym(Xbar) = 0` and `Abar = W X' + W' X`, as
`_ctsem_lyap_pullback!` computes them.
"""
function _tan_lyap_pullback!(Abar::TanMat{W}, Wt::TanMat{W}, h::TanLyap{W}, A::TanMat{W}, X::TanMat{W},
    Xbar::TanMat{W}) where {W}
    _tan_sym!(h.rhs, Xbar)
    _tan_lyap!(Wt, A, h.rhs, h, true)
    _tan_mulNT!(Abar, Wt, X)
    _tan_mulTN!(Abar, Wt, X, true, true)
    return Abar, Wt
end

# The block of the tangent array `d` at (ii, jj) for the k-th unknown of the
# block (is.., js..) with p rows, column-major as the engine's `R` vector.
@inline _tan_blk(d, V, is, js, p, k) = _tan_off(d, V, is + (k - 1) % p, js + (k - 1) ÷ p)

"""
`_ctsem_quasitri_lyap!` on every direction's right-hand side at once: `d` is
the `W x n x n` tangent array of `C`, overwritten with the solutions. Each
block's right-hand side is written in place over its own entries, which no
other entry of the block reads.
"""
function _tan_quasitri_lyap!(d::Array{Float64,3}, S::Matrix{Float64}, n::Int,
    M::Matrix{Float64}, transposed::Bool, V::Val{W}) where {W}
    if !transposed
        je = n
        @inbounds while je >= 1
            js = _ctsem_block_top(S, je)
            for jj in js:je, ii in (je + 1):n
                _tan_st!(d, _tan_off(d, V, ii, jj), _tan_ld(d, _tan_off(d, V, jj, ii), V))
            end
            ie = je
            while ie >= 1
                is = _ctsem_block_top(S, ie)
                for jj in js:je, ii in is:ie
                    o = _tan_off(d, V, ii, jj)
                    t = _tan_ld(d, o, V)
                    for l in (ie + 1):n
                        t = _tan_axpy(t, -S[ii, l], d, l, jj)
                    end
                    for l in (je + 1):n
                        t = _tan_axpy(t, -S[jj, l], d, ii, l)
                    end
                    _tan_st!(d, o, t)
                end
                _tan_small_sylvester!(d, S, is, ie, js, je, M, false, V)
                ie = is - 1
            end
            je = js - 1
        end
    else
        js = 1
        @inbounds while js <= n
            je = _ctsem_block_end(S, js, n)
            for jj in js:je, ii in 1:(js - 1)
                _tan_st!(d, _tan_off(d, V, ii, jj), _tan_ld(d, _tan_off(d, V, jj, ii), V))
            end
            is = js
            while is <= n
                ie = _ctsem_block_end(S, is, n)
                for jj in js:je, ii in is:ie
                    o = _tan_off(d, V, ii, jj)
                    t = _tan_ld(d, o, V)
                    for l in 1:(is - 1)
                        t = _tan_axpy(t, -S[l, ii], d, l, jj)
                    end
                    for l in 1:(js - 1)
                        t = _tan_axpy(t, -S[l, jj], d, ii, l)
                    end
                    _tan_st!(d, o, t)
                end
                _tan_small_sylvester!(d, S, is, ie, js, je, M, true, V)
                is = ie + 1
            end
            js = je + 1
        end
    end
    return d
end

"""
`_ctsem_small_sylvester!` with every direction's right-hand side, held in
place in `d`: the Kronecker system is built and eliminated once, and its row
operations applied to the W-tuples.
"""
function _tan_small_sylvester!(d::Array{Float64,3}, S::Matrix{Float64}, is::Int, ie::Int,
    js::Int, je::Int, M::Matrix{Float64}, transposed::Bool, V::Val{W}) where {W}
    p = ie - is + 1
    q = je - js + 1
    m = p * q
    @inbounds begin
        if m == 1
            o = _tan_off(d, V, is, js)
            _tan_st!(d, o, _tan_times(inv(S[is, is] + S[js, js]), _tan_ld(d, o, V)))
            return d
        end
        for c in 1:m, r in 1:m
            M[r, c] = 0.0
        end
        for b in 1:q, a in 1:p, a2 in 1:p
            M[a + (b - 1) * p, a2 + (b - 1) * p] +=
                transposed ? S[is + a2 - 1, is + a - 1] : S[is + a - 1, is + a2 - 1]
        end
        for b in 1:q, b2 in 1:q, a in 1:p
            M[a + (b - 1) * p, a + (b2 - 1) * p] +=
                transposed ? S[js + b2 - 1, js + b - 1] : S[js + b - 1, js + b2 - 1]
        end
        for c in 1:(m - 1)
            piv = c
            big = abs(M[c, c])
            for r in (c + 1):m
                if abs(M[r, c]) > big
                    big = abs(M[r, c])
                    piv = r
                end
            end
            if piv != c
                for cc in c:m
                    M[c, cc], M[piv, cc] = M[piv, cc], M[c, cc]
                end
                oc, op = _tan_blk(d, V, is, js, p, c), _tan_blk(d, V, is, js, p, piv)
                tc = _tan_ld(d, oc, V)
                _tan_st!(d, oc, _tan_ld(d, op, V))
                _tan_st!(d, op, tc)
            end
            oc = _tan_blk(d, V, is, js, p, c)
            for r in (c + 1):m
                f = M[r, c] / M[c, c]
                for cc in (c + 1):m
                    M[r, cc] -= f * M[c, cc]
                end
                orr = _tan_blk(d, V, is, js, p, r)
                _tan_st!(d, orr, _tan_axpyo(_tan_ld(d, orr, V), -f, d, oc))
            end
        end
        for r in m:-1:1
            orr = _tan_blk(d, V, is, js, p, r)
            t = _tan_ld(d, orr, V)
            for cc in (r + 1):m
                t = _tan_axpyo(t, -M[r, cc], d, _tan_blk(d, V, is, js, p, cc))
            end
            _tan_st!(d, orr, _tan_times(inv(M[r, r]), t))
        end
    end
    return d
end

