# The deviations the subject filter integrates (`CTSEMHybrid`,
# `hybrid_filter.jl`) with their tangents by hand: the columns `X`, the
# information `S` and score `s` they gather through each update, the closing
# term `1/2 b' M^-1 b - w/2 logdet M` (`M = I + G' S G`, `b = G' s`), and the
# reverse of each -- the operations of `_hybrid_update!`, `_hybrid_predict!`,
# `_reverse_hybrid_marginal!`, `_hybrid_reverse_update!` and
# `_hybrid_reverse_predict!`, each paired with its tangent.
#
# Every hybrid item has a column here, where the filter gives a subject columns
# only for the manifest intercepts it observes: an unobserved item's column
# starts at zero and stays there, its rows of `S` and `s` are zero, so it adds
# exact zeros to `M` and `b`, and `log 1` to the logdet -- the same numbers,
# at sizes fixed by the model. The factor `G` (`nh x nh`, an own variance's
# diagonal or a block's rows) is the parameter layer's, global, and so is the
# cotangent it gathers (`TanAcc.Gbar`), which `tangent_hessian.jl` pushes back.

"""`v = A[:, c]`, with tangents."""
function _tan_getcol!(v::TanVec{W}, A::TanMat{W}, c::Int) where {W}
    n = length(v.v)
    on = _tan_nz(A)
    @inbounds for r in 1:n
        v.v[r] = A.v[r, c]
        on && _tan_put!(v, r, _tan_col(A, r, c))
    end
    on ? (v.nz[] = true) : _tan_setnz!(v, false)
    return v
end

"""`A[:, c] += alpha v`, with tangents."""
function _tan_addcol!(A::TanMat{W}, c::Int, alpha::Float64, v::TanVec{W}) where {W}
    n = length(v.v)
    on = _tan_nz(v)
    @inbounds for r in 1:n
        A.v[r, c] = muladd(alpha, v.v[r], A.v[r, c])
        on && _tan_put!(A, r, c, _tan_axpby(1.0, _tan_col(A, r, c), alpha, _tan_col(v, r)))
    end
    on && (A.nz[] = true)
    return A
end

"""A row's share of the hybrid's scratch, for a row observing `mo` manifests."""
struct TanHybRow{W}
    V::TanMat{W}; Wm::TanMat{W}; Z::TanMat{W}; W2::TanMat{W}; Vbar::TanMat{W}
    Kbar::TanMat{W}; Rbar::TanMat{W}; Lbar::TanMat{W}; ybar::TanVec{W}
end
TanHybRow(n::Int, mo::Int, nh::Int, W::Int) = TanHybRow{W}(
    (TanMat{W}(mo, nh) for _ in 1:5)..., TanMat{W}(n, mo), TanMat{W}(mo, mo),
    TanMat{W}(mo, n), TanVec{W}(mo))

"""
A subject's hybrid state and its closing's scratch: the columns `X` and the
cotangent `Xbar` they carry back; `S`, `s`; and `Sbar`, `sbar`, `Gbar`, the
closing term's cotangents, which the reverse walk reads at every row.
"""
struct TanHybState{W}
    X::TanMat{W}; Xt::TanMat{W}; Xbar::TanMat{W}
    S::TanMat{W}; s::TanVec{W}
    GS::TanMat{W}; M::TanMat{W}; UiM::TanMat{W}; Minv::TanMat{W}; Mbar::TanMat{W}
    GM::TanMat{W}; Xc::TanMat{W}; Tc::TanMat{W}; F::Matrix{Float64}
    b::TanVec{W}; beta::TanVec{W}
    Sbar::TanMat{W}; sbar::TanVec{W}; Gbar::TanMat{W}
    col::TanVec{W}; inc::TanVec{W}; sv::TanVec{W}
end
function TanHybState(n::Int, nh::Int, W::Int)
    T, V = TanMat{W}, TanVec{W}
    return TanHybState{W}(T(n, nh), T(n, nh), T(n, nh), T(nh, nh), V(nh),
        (T(nh, nh) for _ in 1:8)..., zeros(nh, nh), V(nh), V(nh), T(nh, nh), V(nh),
        T(nh, nh), V(n), V(n), V(n))
end

"""
A subject's start: every column zero, an initial mean's (`kind` 2) the unit
vector at its latent; `S` and `s` zero.
"""
function _tan_hybrid_start!(hs::TanHybState{W}, h::CTSEMHybrid) where {W}
    _tan_zero!(hs.X)
    _tan_zero!(hs.S)
    _tan_zero!(hs.s)
    @inbounds for k in eachindex(h.kind)
        h.kind[k] == 2 && (hs.X.v[h.target[k], k] = 1.0)
    end
    return hs
end

"""
An update's share (`_hybrid_update!`): `V = L X + E`, `W = Ui' V`, `S += W'W`,
`s += W' e`, `X -= K V`; `L` the observed rows of LAMBDA, `Ui`, `e = Ui' y`
and `K` the update's own.
"""
function _tan_hybrid_update!(hs::TanHybState{W}, hr::TanHybRow{W}, h::CTSEMHybrid,
    o::AbstractVector{Int}, L::TanMat{W}, Ui::TanMat{W}, e::TanVec{W},
    K::TanMat{W}) where {W}
    V = hr.V
    _tan_mul!(V, L, hs.X)
    _tan_hybrid_selector!(V, h, o)
    _tan_mulTN!(hr.Wm, Ui, V, true, false, Val(:AupperT))
    _tan_mulTN!(hs.S, hr.Wm, hr.Wm, true, true, Val(:full), Val(true))
    _tan_mulTvec!(hs.s, hr.Wm, e, true, true)
    _tan_mul!(hs.X, K, V, -1.0, true)
    return hs
end

"""`V += E`: a manifest intercept's column (`kind` 1) at its item's row."""
function _tan_hybrid_selector!(V::TanMat, h::CTSEMHybrid, o::AbstractVector{Int})
    @inbounds for k in eachindex(h.kind)
        h.kind[k] == 1 || continue
        for i in eachindex(o)
            o[i] == h.target[k] && (V.v[i, k] += 1.0)
        end
    end
    return V
end

"""
A prediction's share (`_hybrid_predict!`): `X = A X`, and a continuous
intercept's column (`kind` 3) gains `J^-1 (A - I) e_l` (`Xinc`, recorded), by
the drift's factors (`FJ`) as the intercept is, or `Phi e_l` where the
interval's intercept was taken as a series (`Phi` not `nothing`).
"""
function _tan_hybrid_predict!(hs::TanHybState{W}, Xh::TanMat{W}, Xinc::TanMat{W},
    h::CTSEMHybrid, A::TanMat{W}, JAx::TanMat{W}, FJ::Matrix{Float64},
    pivJ::Vector{Int}, Phi::Union{Nothing,TanMat{W}}) where {W}
    _tan_copy!(Xh, hs.X)
    _tan_mul!(hs.X, A, Xh)
    _tan_zero!(Xinc)
    @inbounds for k in eachindex(h.kind)
        h.kind[k] == 3 || continue
        l = h.target[k]
        inc = hs.inc
        if Phi === nothing
            _tan_getcol!(inc, A, l)
            inc.v[l] -= 1.0
            _tan_solve!(inc, FJ, pivJ, JAx)
        else
            _tan_getcol!(inc, Phi, l)
        end
        _tan_addcol!(Xinc, k, 1.0, inc)
        _tan_addcol!(hs.X, k, 1.0, inc)
    end
    return hs
end

"""
The closing term and its reverse (`_hybrid_marginal!`,
`_reverse_hybrid_marginal!`): adds `1/2 b' beta - w/2 logdet M` and its
tangents to `ll`, or with deviations held at `dfix` (`_HYBRID_FIXED`, its
tangents zero) `d' b - 1/2 d' (M - I) d`; and leaves `Sbar`, `sbar` for the
walk and `Gbar` for the layer. `false` where `M` does not factor.
"""
function _tan_hybrid_closing!(hs::TanHybState{W}, G::TanMat{W}, w::Float64,
    ll::Vector{Float64}, dfix::Union{Nothing,TanVec{W}}=nothing) where {W}
    nh = size(G.v, 1)
    M, b, beta = hs.M, hs.b, hs.beta
    _tan_mulTN!(hs.GS, G, hs.S)
    _tan_mul!(M, hs.GS, G)
    _tan_diag!(M, 1.0)
    _tan_sym!(hs.Xc, M)
    _tan_copy!(M, hs.Xc)
    _tan_mulTvec!(b, G, hs.s)
    _tan_cholinv!(hs.UiM, M, hs.Xc, hs.Tc, hs.F) || return false
    UiM = hs.UiM
    _tan_mulNT!(hs.Minv, UiM, UiM, true, false, Val(:full), Val(true))
    _tan_mulvec!(beta, hs.Minv, b)
    dfix === nothing || return _tan_hybrid_fixed!(hs, G, ll, dfix)
    # ll += 1/2 b' beta - w/2 logdet M, logdet M = -2 sum log UiM_ii.
    quad = 0.0
    ld = 0.0
    @inbounds for a in 1:nh
        quad += b.v[a] * beta.v[a]
        ld -= 2 * log(UiM.v[a, a])
    end
    ll[1] += quad / 2 - w * ld / 2
    @inbounds for t in 1:W
        dq = 0.0
        dld = 0.0
        for a in 1:nh
            dq += b.d[t, a] * beta.v[a] + b.v[a] * beta.d[t, a]
            dld -= 2 * UiM.d[t, a, a] / UiM.v[a, a]
        end
        ll[1 + t] += dq / 2 - w * dld / 2
    end
    # Mbar = -(w M^-1 + beta beta')/2; Sbar = G Mbar G', sbar = G beta,
    # Gbar = 2 S G Mbar + s beta'.
    _tan_comb!(hs.Mbar, ((-w / 2, hs.Minv),))
    _tan_ger!(hs.Mbar, -0.5, beta, beta)
    _tan_mul!(hs.GM, G, hs.Mbar)
    _tan_mulNT!(hs.Sbar, hs.GM, G)
    _tan_mulvec!(hs.sbar, G, beta)
    _tan_mul!(hs.Gbar, hs.S, hs.GM, 2.0, false)
    _tan_ger!(hs.Gbar, true, hs.s, beta)
    _tan_zero!(hs.Xbar)
    return true
end

"""
The closing term at fixed deviations `d` (`_tan_hybrid_closing!`): adds
`d' b - 1/2 d' (M - I) d` and its tangents to `ll`; `Mbar = -d d' / 2` and
`bbar = d`, so `Sbar = G Mbar G'`, `sbar = G d`, `Gbar = 2 S G Mbar + s d'`.
"""
function _tan_hybrid_fixed!(hs::TanHybState{W}, G::TanMat{W}, ll::Vector{Float64},
    d::TanVec{W}) where {W}
    nh = size(G.v, 1)
    M, b = hs.M, hs.b
    val = 0.0
    @inbounds for a in 1:nh
        val += d.v[a] * b.v[a]
        for c in 1:nh
            val -= d.v[a] * (M.v[a, c] - (a == c ? 1.0 : 0.0)) * d.v[c] / 2
        end
    end
    ll[1] += val
    @inbounds for t in 1:W
        dv = 0.0
        for a in 1:nh
            dv += d.v[a] * b.d[t, a]
            for c in 1:nh
                dv -= d.v[a] * M.d[t, a, c] * d.v[c] / 2
            end
        end
        ll[1 + t] += dv
    end
    _tan_ger!(hs.Mbar, -0.5, d, d, false)
    _tan_mul!(hs.GM, G, hs.Mbar)
    _tan_mulNT!(hs.Sbar, hs.GM, G)
    _tan_mulvec!(hs.sbar, G, d)
    _tan_mul!(hs.Gbar, hs.S, hs.GM, 2.0, false)
    _tan_ger!(hs.Gbar, true, hs.s, d)
    _tan_zero!(hs.Xbar)
    return true
end

"""
An update's reverse share (`_hybrid_reverse_update!`), before the base reverse
uses its intermediates: from the carried `Xbar` and the closing's `Sbar`,
`sbar`, the cotangents of the gain (`Kbar`), the innovation covariance
(`Rbar`), the innovation (`ybar`) and the observed LAMBDA rows (`Lbar`), which
the base reverse adds where it forms its own; and `Xbar` carried back to
before the update. `Ui`, `K`, `a = S^-1 y` as the base reverse recomputed them.
"""
function _tan_hybrid_reverse_update!(hs::TanHybState{W}, hr::TanHybRow{W}, h::CTSEMHybrid,
    o::AbstractVector{Int}, Xh::TanMat{W}, L::TanMat{W}, Ui::TanMat{W}, K::TanMat{W},
    a::TanVec{W}) where {W}
    V, Wm, Z, W2, Vbar = hr.V, hr.Wm, hr.Z, hr.W2, hr.Vbar
    _tan_mul!(V, L, Xh)
    _tan_hybrid_selector!(V, h, o)
    _tan_mulTN!(Wm, Ui, V, true, false, Val(:AupperT))  # W = Ui' V
    _tan_mul!(Z, Ui, Wm, true, false, Val(:Aupper))  # Z = Ui W = S^-1 V
    Xbar = hs.Xbar
    # X+ = X - K V: Kbar = -Xbar V', Vbar = -K' Xbar.
    _tan_mulNT!(hr.Kbar, Xbar, V, -1.0, false)
    _tan_mul!(W2, Z, hs.Sbar)
    # S += V' S^-1 V: Vbar += 2 Z Sbar. s += V' S^-1 y: Vbar += a sbar'.
    _tan_mulTN!(Vbar, K, Xbar, -1.0, false)
    _tan_comb!(Vbar, ((2.0, W2),), true)
    _tan_ger!(Vbar, true, a, hs.sbar)
    # ybar = Z sbar; Rbar = -Z Sbar Z' - sym(Z sbar a').
    _tan_mulvec!(hr.ybar, Z, hs.sbar)
    _tan_mulNT!(hr.Rbar, W2, Z, -1.0, false)
    _tan_ger!(hr.Rbar, -0.5, hr.ybar, a)
    _tan_ger!(hr.Rbar, -0.5, a, hr.ybar)
    # V = L X + E: Lbar = Vbar X', Xbar(before) = Xbar(after) + L' Vbar.
    _tan_mulNT!(hr.Lbar, Vbar, Xh)
    _tan_mulTN!(Xbar, L, Vbar, true, true)
    return hr
end

"""
A prediction's reverse share (`_hybrid_reverse_predict!`): a continuous
intercept's column first, `sbar = J^-T Xbar[:, c]`, `JAxbar -= sbar Xinc[:, c]'`,
`Abar[:, l] += sbar` -- or, where the interval's intercept was a series,
`Phibar[:, l] += Xbar[:, c]`, which the caller pulls back with its own; then
`Abar += Xbar X'` and `Xbar = A' Xbar`. `JT`, `FT` are JAx' and its factors.
"""
function _tan_hybrid_reverse_predict!(hs::TanHybState{W}, h::CTSEMHybrid, Xh::TanMat{W},
    Xinc::TanMat{W}, A::TanMat{W}, Abar::TanMat{W}, JAxbar::TanMat{W}, JT::TanMat{W},
    FT::Matrix{Float64}, pivT::Vector{Int}, Phibar::Union{Nothing,TanMat{W}}) where {W}
    Xbar = hs.Xbar
    @inbounds for k in eachindex(h.kind)
        h.kind[k] == 3 || continue
        l = h.target[k]
        sv, inc = hs.sv, hs.inc
        _tan_getcol!(sv, Xbar, k)
        if Phibar !== nothing
            _tan_addcol!(Phibar, l, 1.0, sv)
            continue
        end
        _tan_solve!(sv, FT, pivT, JT)
        _tan_getcol!(inc, Xinc, k)
        _tan_ger!(JAxbar, -1.0, sv, inc)
        _tan_addcol!(Abar, l, 1.0, sv)
    end
    _tan_mulNT!(Abar, Xbar, Xh, true, true)
    _tan_copy!(hs.Xt, Xbar)
    _tan_mulTN!(Xbar, A, hs.Xt)
    return hs
end
