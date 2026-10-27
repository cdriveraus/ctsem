# The series discretisation (`series_discretization.jl`) with its tangents by
# hand: where a pivot is small against the interval the filter takes the
# intercept's integral `Phi(dt) = int_0^dt e^{J s} ds` and the noise's
# `V(dt) = int_0^dt e^{J s} Q e^{J' s} ds` as truncated Taylor series on a
# scaled interval and then doublings. Each kernel and pullback here is the one
# there -- `_series_intercept!`, `_series_noise!` and their pullbacks, Horner's
# rule rather than the power tables, which agree to rounding -- with every
# product carrying its tangents. The plan (doublings, terms, scaled step) is
# the primal's, as there. The reverse runs the forward sweep again from the
# record's inputs rather than taping its chain.

"""A series kernel's sweep, with tangents (`SeriesKernelBuffer`)."""
mutable struct TanSeries{W}
    n::Int
    X::TanMat{W}; Q::TanMat{W}; out::TanMat{W}
    horner::Vector{TanMat{W}}; taylor::Vector{TanMat{W}}
    Es::Vector{TanMat{W}}; Ps::Vector{TanMat{W}}
    b1::TanMat{W}; b2::TanMat{W}; b3::TanMat{W}; Ebar::TanMat{W}; Xbar::TanMat{W}
    s::Int
    m::Int
    h0::Float64
end
TanSeries(n::Int, W::Int) = TanSeries{W}(n, (TanMat{W}(n, n) for _ in 1:3)...,
    TanMat{W}[], TanMat{W}[], TanMat{W}[], TanMat{W}[], (TanMat{W}(n, n) for _ in 1:5)...,
    0, 0, 0.0)

@inline function _tan_series_grow!(tape::Vector{TanMat{W}}, count::Int, n::Int) where {W}
    while length(tape) < count
        push!(tape, TanMat{W}(n, n))
    end
    return tape
end

"""Doublings and scaled step for `h J`, Taylor terms for `factor` times its
norm, and `X = J h0` (`_series_plan!`)."""
function _tan_series_plan!(sb::TanSeries{W}, J::TanMat{W}, h::Float64,
    factor::Float64) where {W}
    n = sb.n
    nrm = _series_norm1(J.v, n) * h
    s = 0
    if isfinite(nrm) && nrm > _CTSEM_SERIES_THETA
        s = min(ceil(Int, log2(nrm / _CTSEM_SERIES_THETA)), _CTSEM_SERIES_MAXDOUBLINGS)
    end
    scale = 2.0^s
    h0 = h / scale
    z = isfinite(nrm) ? factor * nrm / scale : 0.0
    sb.s = s
    sb.m = _series_terms(z)
    sb.h0 = h0
    _tan_series_grow!(sb.horner, sb.m, n)
    _tan_series_grow!(sb.Es, s, n)
    _tan_series_grow!(sb.Ps, s, n)
    _tan_comb!(sb.X, ((h0, J),))
    return sb
end

"""`Phi(h)` into `sb.out` (`_series_intercept!`): `Phi = h0 G`, `G` the series
of `phi_1(X)` by Horner's rule, then the doublings `Phi(2h) = Phi + e^{Jh} Phi`."""
function _tan_series_intercept!(sb::TanSeries{W}, J::TanMat{W}, h::Float64) where {W}
    _tan_series_plan!(sb, J, h, 1.0)
    s, m, h0 = sb.s, sb.m, sb.h0
    X = sb.X
    Gm = sb.horner[m]
    _tan_comb!(Gm, ((1 / (m + 1), X),))
    _tan_diag!(Gm, 1.0)
    for j in (m - 1):-1:1
        Gj = sb.horner[j]
        _tan_mul!(Gj, X, sb.horner[j + 1], 1 / (j + 1), false)
        _tan_diag!(Gj, 1.0)
    end
    G = sb.horner[1]
    out = sb.out
    _tan_comb!(s == 0 ? out : sb.Ps[1], ((h0, G),))
    s == 0 && return out
    E0 = sb.Es[1]
    _tan_mul!(E0, X, G)
    _tan_diag!(E0, 1.0)
    for k in 1:s
        E, P = sb.Es[k], sb.Ps[k]
        next = k == s ? out : sb.Ps[k + 1]
        _tan_mul!(next, E, P)
        _tan_comb!(next, ((1.0, P),), true)
        k < s && _tan_mul!(sb.Es[k + 1], E, E)
    end
    return out
end

"""The intercept's doubling chain reversed (`_series_intercept_doubling_pullback!`):
the cotangents of `Phi(h0)` (`b1`) and of `e^{J h0}` (`Ebar`)."""
function _tan_series_intercept_doubling_pullback!(sb::TanSeries{W}, Phibar::TanMat{W}) where {W}
    Pb, Eb, t1, t2 = sb.b1, sb.Ebar, sb.b2, sb.b3
    _tan_copy!(Pb, Phibar)
    _tan_zero!(Eb)
    for k in sb.s:-1:1
        E, P = sb.Es[k], sb.Ps[k]
        if k < sb.s
            _tan_mulNT!(t1, Eb, E)
            _tan_mulTN!(t1, E, Eb, true, true)
        else
            _tan_zero!(t1)
        end
        _tan_mulNT!(t1, Pb, P, true, true)
        _tan_mulTN!(t2, E, Pb)
        _tan_comb!(Pb, ((1.0, t2),), true)
        _tan_copy!(Eb, t1)
    end
    return Pb, Eb
end

"""`Jbar += ` the intercept series' pullback for `Phibar`
(`_series_intercept_pullback!`); `sb` holds its forward sweep."""
function _tan_series_intercept_pullback!(Jbar::TanMat{W}, sb::TanSeries{W},
    Phibar::TanMat{W}) where {W}
    s, m, h0 = sb.s, sb.m, sb.h0
    X, t1 = sb.X, sb.b2
    Pb, Eb = _tan_series_intercept_doubling_pullback!(sb, Phibar)
    G = sb.horner[1]
    Gb = sb.b3
    _tan_comb!(Gb, ((h0, Pb),))
    Xb = sb.Xbar
    if s > 0
        _tan_mulTN!(Gb, X, Eb, true, true)
        _tan_mulNT!(Xb, Eb, G)
    else
        _tan_zero!(Xb)
    end
    for j in 1:m
        c = 1 / (j + 1)
        if j < m
            _tan_mulNT!(Xb, Gb, sb.horner[j + 1], c, true)
            _tan_mulTN!(t1, X, Gb, c)
            _tan_copy!(Gb, t1)
        else
            _tan_comb!(Xb, ((c, Gb),), true)
        end
    end
    _tan_comb!(Jbar, ((h0, Xb),), true)
    return Jbar
end

"""`V(h)` into `sb.out` (`_series_noise!`): `V = h0 W`, `W` by Horner's rule on
`L(W) = X W + W X'`; `e^X` its own series for the doublings
`V(2h) = V + e^{Jh} V e^{J'h}`; made exactly symmetric."""
function _tan_series_noise!(sb::TanSeries{W}, J::TanMat{W}, Q::TanMat{W},
    h::Float64) where {W}
    _tan_series_plan!(sb, J, h, 2.0)
    s, m, h0 = sb.s, sb.m, sb.h0
    X, Qb, t1 = sb.X, sb.Q, sb.b1
    _tan_copy!(Qb, Q)
    for j in m:-1:1
        prev = j == m ? Qb : sb.horner[j + 1]
        Wj = sb.horner[j]
        _tan_mul!(t1, X, prev)
        c = 1 / (j + 1)
        _tan_copy!(Wj, Qb)
        _tan_add!(Wj, c, t1, true)
        _tan_add!(Wj, c, t1, true, Val(:T))
    end
    out = sb.out
    _tan_comb!(s == 0 ? out : sb.Ps[1], ((h0, sb.horner[1]),))
    if s > 0
        _tan_series_grow!(sb.taylor, m, sb.n)
        Fm = sb.taylor[m]
        _tan_comb!(Fm, ((1 / m, X),))
        _tan_diag!(Fm, 1.0)
        for j in (m - 1):-1:1
            Fj = sb.taylor[j]
            _tan_mul!(Fj, X, sb.taylor[j + 1], 1 / j, false)
            _tan_diag!(Fj, 1.0)
        end
        _tan_copy!(sb.Es[1], sb.taylor[1])
        for k in 1:s
            E, P = sb.Es[k], sb.Ps[k]
            next = k == s ? out : sb.Ps[k + 1]
            _tan_mul!(t1, E, P)
            _tan_mulNT!(next, t1, E)
            _tan_comb!(next, ((1.0, P),), true)
            k < s && _tan_mul!(sb.Es[k + 1], E, E)
        end
    end
    _tan_sym!(t1, out)
    _tan_copy!(out, t1)
    return out
end

"""The noise's doubling chain reversed (`_series_noise_doubling_pullback!`),
`Vbar` symmetrised first: the cotangents of `V(h0)` (`b1`) and `e^{J h0}`
(`Ebar`)."""
function _tan_series_noise_doubling_pullback!(sb::TanSeries{W}, Vbar::TanMat{W}) where {W}
    Vb, Eb, t1, t2 = sb.b1, sb.Ebar, sb.b2, sb.b3
    _tan_sym!(Vb, Vbar)
    _tan_zero!(Eb)
    for k in sb.s:-1:1
        E, P = sb.Es[k], sb.Ps[k]
        if k < sb.s
            _tan_mulNT!(t1, Eb, E)
            _tan_mulTN!(t1, E, Eb, true, true)
        else
            _tan_zero!(t1)
        end
        _tan_mul!(t2, Vb, E)
        _tan_mul!(t1, t2, P, 2.0, true)
        _tan_copy!(Eb, t1)
        _tan_mulTN!(t1, E, t2)
        _tan_comb!(Vb, ((1.0, t1),), true)
    end
    return Vb, Eb
end

"""`Jbar += `, `Qbar += ` the noise series' pullback for `Vbar`
(`_series_noise_pullback!`); `sb` holds its forward sweep."""
function _tan_series_noise_pullback!(Jbar::TanMat{W}, Qbar::TanMat{W}, sb::TanSeries{W},
    Vbar::TanMat{W}) where {W}
    s, m, h0 = sb.s, sb.m, sb.h0
    X, Qb, t1, Xb = sb.X, sb.Q, sb.b2, sb.Xbar
    Vb, Eb = _tan_series_noise_doubling_pullback!(sb, Vbar)
    _tan_zero!(Xb)
    if s > 0
        Fb = Eb
        for j in 1:m
            c = 1 / j
            if j < m
                _tan_mulNT!(Xb, Fb, sb.taylor[j + 1], c, true)
                _tan_mulTN!(t1, X, Fb, c)
                _tan_copy!(Fb, t1)
            else
                _tan_comb!(Xb, ((c, Fb),), true)
            end
        end
    end
    Wb = Vb
    _tan_scale!(Wb, h0)
    for j in 1:m
        c = 1 / (j + 1)
        next = j < m ? sb.horner[j + 1] : Qb
        _tan_comb!(Qbar, ((1.0, Wb),), true)
        _tan_mul!(Xb, Wb, next, 2c, true)
        _tan_mulTN!(t1, X, Wb, c)
        _tan_mul!(t1, Wb, X, c, true)
        _tan_copy!(Wb, t1)
    end
    _tan_comb!(Qbar, ((1.0, Wb),), true)
    _tan_comb!(Jbar, ((h0, Xb),), true)
    return Jbar, Qbar
end
