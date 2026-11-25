# The series discretisation (`series_discretization.jl`) with its tangents by
# hand: where a pivot is small against the interval the filter takes the
# intercept's integral `Phi(dt) = int_0^dt e^{J s} ds` and the noise's
# `V(dt) = int_0^dt e^{J s} Q e^{J' s} ds` as truncated Taylor series on a
# scaled interval and then doublings. Each kernel and pullback here is the one
# there, with every product carrying its tangents: `_series_intercept!`,
# `_series_noise!` and their pullbacks by Horner's rule, where the drift can
# change by row, and the power tables (`SeriesTable`, `TanSeriesTable` below)
# where it holds over a subject's rows. The plan (doublings, terms, scaled step)
# is the primal's, as there. The reverse runs an interval's forward again from
# the record's inputs rather than taping its chain.

"""A series kernel's sweep, with tangents (`SeriesKernelBuffer`)."""
mutable struct TanSeries{W}
    n::Int
    X::TanMat{W}; Q::TanMat{W}; out::TanMat{W}
    horner::Vector{TanMat{W}}; taylor::Vector{TanMat{W}}
    Es::Vector{TanMat{W}}; Ps::Vector{TanMat{W}}
    b1::TanMat{W}; b2::TanMat{W}; b3::TanMat{W}; Ebar::TanMat{W}; Xbar::TanMat{W}
    A::TanMat{W}            # e^{J h} off the chain (`_tan_series_exp!`)
    s::Int
    m::Int
    h0::Float64
end
TanSeries(n::Int, W::Int) = TanSeries{W}(n, (TanMat{W}(n, n) for _ in 1:3)...,
    TanMat{W}[], TanMat{W}[], TanMat{W}[], TanMat{W}[], (TanMat{W}(n, n) for _ in 1:6)...,
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

# The doublings by their nonzero terms. A drift is usually sparse -- processes
# that do not act on each other, a level a faster process reverts to, an
# average read off one process -- and what the doublings multiply is then
# sparse too: `e^{J h}` and the intercept's integral are nonzero only where
# one process reaches another (the pattern's transitive closure), the noise's
# integral only between processes that noise reaches together. Each product of
# a doubling is listed once for a drift's pattern (`TanProdList`) and summed
# over its terms; the reverse's cotangents are taken as dense, so its products
# keep the terms of their primal factor. On the affect pilot's tied model (13
# processes, 37 of the exponential's 169 cells) the doublings were most of a
# Laplace value and gradient, 5.2 of them a row (2026-11-18, dev2).
# `ctsem_set_series_lists!(false)` is the oracle.
"""
The doubling products' terms for a drift's and a diffusion's patterns (`mJ`,
`mQ`): `EE` is `E E` and `E Phi`, `EV` is `E V`, `TEt` is `(E V) E'`; of the
reverse, with `D` a dense cotangent, `DEt` is `D E'` and `D Phi'`, `EtD` is
`E' D`, `DE` is `D E`, `ED` is `E D`, `DV` is `D V`, `DVt` is `D V'` and `DtV`
is `D' V`. `on` where the pattern is sparse enough to pay.
"""
mutable struct TanSeriesLists
    on::Bool
    built::Bool             # `mJ` and `mQ` are a drift's and a diffusion's
    ready::Bool             # and the lists are theirs (sparse enough to list)
    mJ::Matrix{Bool}
    mQ::Matrix{Bool}
    sJ::Matrix{Bool}        # scratch for the patterns in hand
    sQ::Matrix{Bool}
    EE::TanProdList
    EV::TanProdList
    TEt::TanProdList
    DEt::TanProdList
    EtD::TanProdList
    DE::TanProdList
    ED::TanProdList
    DV::TanProdList
    DVt::TanProdList
    DtV::TanProdList
    eE::Vector{Int32}       # the exponential's entries, and the noise integral's
    eV::Vector{Int32}
end
TanSeriesLists(n::Int) = TanSeriesLists(false, false, false, fill(false, n, n), fill(false, n, n),
    fill(false, n, n), fill(false, n, n),
    _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST,
    _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST, _TAN_NO_LIST, Int32[], Int32[])

"""Listed doublings on or off; off is the dense products, the oracle."""
const _TAN_SERIES_LISTS = Ref(true)
ctsem_set_series_lists!(on::Bool) = (_TAN_SERIES_LISTS[] = on)
"""The largest share of the exponential's cells that may be nonzero for the
doublings to go by lists, and the fewest processes."""
const _TAN_SERIES_LIST_SHARE = 0.5
const _TAN_SERIES_LIST_MIN = 4

# `m` with its diagonal, closed under composition: where one process reaches
# another.
function _tan_reach(m::Matrix{Bool})
    n = size(m, 1)
    R = copy(m)
    for i in 1:n
        R[i, i] = true
    end
    changed = true
    while changed
        changed = false
        for j in 1:n, i in 1:n
            R[i, j] && continue
            for q in 1:n
                if R[i, q] && R[q, j]
                    R[i, j] = true
                    changed = true
                    break
                end
            end
        end
    end
    return R
end
_tan_bmul(A::Matrix{Bool}, B::Matrix{Bool}) =
    Bool[any(A[i, q] && B[q, j] for q in axes(A, 2)) for i in axes(A, 1), j in axes(B, 2)]

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
    return _tan_series_doublings!(sb, out, false)
end

"""The doubling chain both series end with (`_series_doublings!`), from
`sb.Es[1]` and `sb.Ps[1]` into `out`: `P_{k+1} = P_k + E_k P_k E_k'` for the
noise, `P_k + E_k P_k` for the intercept, `E_{k+1} = E_k^2`."""
function _tan_series_doublings!(sb::TanSeries{W}, out::TanMat{W}, noise::Bool,
    lists::Union{Nothing,TanSeriesLists}=nothing) where {W}
    s, t1 = sb.s, sb.b1
    listed = lists !== nothing && lists.on
    for k in 1:s
        E, P = sb.Es[k], sb.Ps[k]
        next = k == s ? out : sb.Ps[k + 1]
        if listed
            if noise
                _tan_mul_list!(t1, E, P, lists.EV)
                _tan_mul_list!(next, t1, E, lists.TEt)
            else
                _tan_mul_list!(next, E, P, lists.EE)
            end
        elseif noise
            _tan_mul!(t1, E, P)
            _tan_mulNT!(next, t1, E)
        else
            _tan_mul!(next, E, P)
        end
        _tan_comb!(next, ((1.0, P),), true)
        if k < s
            listed ? _tan_mul_list!(sb.Es[k + 1], E, E, lists.EE) :
                _tan_mul!(sb.Es[k + 1], E, E)
        end
    end
    return out
end

"""The intercept's doubling chain reversed (`_series_intercept_doubling_pullback!`):
the cotangents of `Phi(h0)` (`b1`) and of `e^{J h0}` (`Ebar`)."""
function _tan_series_intercept_doubling_pullback!(sb::TanSeries{W}, Phibar::TanMat{W},
    lists::Union{Nothing,TanSeriesLists}=nothing,
    Abar::Union{Nothing,TanMat{W}}=nothing) where {W}
    Pb, Eb, t1, t2 = sb.b1, sb.Ebar, sb.b2, sb.b3
    _tan_copy!(Pb, Phibar)
    # `Abar`: the cotangent of `e^{J h}`, the square of the chain's last
    # exponential (`_tan_series_exp!`), which is what `Eb` is of at every
    # doubling below the last.
    Abar === nothing ? _tan_zero!(Eb) : _tan_copy!(Eb, Abar)
    listed = lists !== nothing && lists.on
    for k in sb.s:-1:1
        E, P = sb.Es[k], sb.Ps[k]
        if k == sb.s && Abar === nothing
            _tan_zero!(t1)
        elseif listed
            _tan_mul_list!(t1, Eb, E, lists.DEt)
            _tan_mul_list!(t1, E, Eb, lists.EtD, 1.0, 1.0)
        else
            _tan_mulNT!(t1, Eb, E)
            _tan_mulTN!(t1, E, Eb, true, true)
        end
        if listed
            _tan_mul_list!(t1, Pb, P, lists.DEt, 1.0, 1.0)
            _tan_mul_list!(t2, E, Pb, lists.EtD)
        else
            _tan_mulNT!(t1, Pb, P, true, true)
            _tan_mulTN!(t2, E, Pb)
        end
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
        _tan_series_doublings!(sb, out, true)
    end
    _tan_sym!(t1, out)
    _tan_copy!(out, t1)
    return out
end

"""The noise's doubling chain reversed (`_series_noise_doubling_pullback!`),
`Vbar` symmetrised first: the cotangents of `V(h0)` (`b1`) and `e^{J h0}`
(`Ebar`)."""
function _tan_series_noise_doubling_pullback!(sb::TanSeries{W}, Vbar::TanMat{W},
    lists::Union{Nothing,TanSeriesLists}=nothing) where {W}
    Vb, Eb, t1, t2 = sb.b1, sb.Ebar, sb.b2, sb.b3
    _tan_sym!(Vb, Vbar)
    _tan_zero!(Eb)
    listed = lists !== nothing && lists.on
    for k in sb.s:-1:1
        E, P = sb.Es[k], sb.Ps[k]
        if k == sb.s
            _tan_zero!(t1)
        elseif listed
            _tan_mul_list!(t1, Eb, E, lists.DEt)
            _tan_mul_list!(t1, E, Eb, lists.EtD, 1.0, 1.0)
        else
            _tan_mulNT!(t1, Eb, E)
            _tan_mulTN!(t1, E, Eb, true, true)
        end
        if listed
            _tan_mul_list!(t2, Vb, E, lists.DE)
            _tan_mul_list!(t1, t2, P, lists.DV, 2.0, 1.0)
            _tan_copy!(Eb, t1)
            _tan_mul_list!(t1, E, t2, lists.EtD)
        else
            _tan_mul!(t2, Vb, E)
            _tan_mul!(t1, t2, P, 2.0, true)
            _tan_copy!(Eb, t1)
            _tan_mulTN!(t1, E, t2)
        end
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

# The tables. One interval's series differs from another's only in its scalar
# weights on the drift's powers, so where the drift holds over a subject's rows
# the powers are built once and an interval's series is their weighted sum,
# O(m n^2) where Horner's rule is O(m n^3); only the doublings stay products.
# The reverse adds each interval's cotangents onto the powers, already
# weighted, and pulls them through the power chains once
# (`_tan_series_release!`), where the exponential's queued directions are
# released (`_tan_flush_frechet!`) -- which every reverse sweep does before its
# cotangents are read. On the affect pilot, whose slow tier puts nearly every
# interval on the series, Horner's rule was half of a Laplace value and
# gradient by hand.

"""
The powers of a drift `J` with tangents (`SeriesTable`): `Jp[j] = J^j`, and
for a noise table `Lp[j + 1] = L^j(Q)`, `L(V) = J V + V J'`; `J` and `Q` the
copies they are of (`_tan_table_for!`). `Jpbar` and `Lpbar` hold the
reverse's cotangents on the powers while `pending`, up to power `bcount`;
`Jbar` and `Qbar` what a release pulls them back to.
"""
mutable struct TanSeriesTable{W}
    n::Int
    noise::Bool
    valid::Bool
    count::Int
    J::TanMat{W}
    Q::TanMat{W}
    Jp::Vector{TanMat{W}}
    Lp::Vector{TanMat{W}}
    tmp::TanMat{W}
    pending::Bool
    bcount::Int
    Jpbar::Vector{TanMat{W}}
    Lpbar::Vector{TanMat{W}}
    Jbar::TanMat{W}
    Qbar::TanMat{W}
    lists::TanSeriesLists
end
TanSeriesTable(n::Int, W::Int, noise::Bool) = TanSeriesTable{W}(n, noise, false, 0,
    TanMat{W}(n, n), TanMat{W}(n, n), TanMat{W}[], TanMat{W}[], TanMat{W}(n, n), false, 0,
    TanMat{W}[], TanMat{W}[], TanMat{W}(n, n), TanMat{W}(n, n), TanSeriesLists(n))

"""`tab`'s lists for the drift and diffusion it holds: kept where their
patterns are the ones the lists were built for, as a model's subjects' are."""
function _tan_series_lists!(tab::TanSeriesTable{W}) where {W}
    L = tab.lists
    n = tab.n
    if !_TAN_SERIES_LISTS[] || n < _TAN_SERIES_LIST_MIN
        L.on = false
        return L
    end
    # The patterns in hand, into scratch; swapped with the lists' own where
    # they differ.
    mJ = _tan_mask!(L.sJ, tab.J)
    mQ = tab.noise ? _tan_mask!(L.sQ, tab.Q) : fill!(L.sQ, false)
    if L.built && mJ == L.mJ && mQ == L.mQ
        L.on = L.ready
        return L
    end
    L.sJ, L.sQ = L.mJ, L.mQ
    L.mJ, L.mQ, L.built, L.ready = mJ, mQ, true, false
    mE = _tan_reach(mJ)
    if count(mE) > _TAN_SERIES_LIST_SHARE * n * n
        L.on = false
        return L
    end
    dense = fill(true, n, n)
    L.EE = _tan_prod_list(mE, false, mE, false)
    L.DEt = _tan_prod_list(dense, false, mE, true)
    L.EtD = _tan_prod_list(mE, true, dense, false)
    L.DE = _tan_prod_list(dense, false, mE, false)
    L.ED = _tan_prod_list(mE, false, dense, false)
    L.eE = _tan_entries(mE)
    if tab.noise
        mV = _tan_bmul(_tan_bmul(mE, mQ), permutedims(mE))
        L.EV = _tan_prod_list(mE, false, mV, false)
        L.TEt = _tan_prod_list(_tan_bmul(mE, mV), false, mE, true)
        L.DV = _tan_prod_list(dense, false, mV, false)
        L.DVt = _tan_prod_list(dense, false, mV, true)
        L.eV = _tan_entries(mV)
        L.DtV = _tan_prod_list(dense, true, mV, false)
    end
    L.ready = true
    L.on = true
    return L
end

"""Whether `A` and `B` are the same, tangents and second parts included."""
function _tan_identical(A::TanMat, B::TanMat)
    A.v == B.v || return false
    _tan_nz(A) == _tan_nz(B) || return false
    _tan_nz(A) && !(A.d == B.d) && return false
    _tan_snz(A) == _tan_snz(B) || return false
    _tan_snz(A) && !(A.s == B.s) && return false
    return true
end

@noinline _tan_table_stale() = throw(ErrorException(string("series cotangents ",
    "pending on a drift the sweep has left: released at the flush, or dropped at ",
    "the sweep's start")))

"""`tab`, holding the powers of `J` (and `Q`, a noise table's): reset to them
where it holds another drift's, built as the series asks (`_tan_table_extend!`).
Cotangents pending on the drift it leaves are an error: a reverse sweep
releases them before its drift changes, and a sweep starts with none
(`_tan_series_drop!`)."""
function _tan_table_for!(tab::TanSeriesTable{W}, J::TanMat{W}, Q) where {W}
    tab.valid && _tan_identical(tab.J, J) && (!tab.noise || _tan_identical(tab.Q, Q)) &&
        return tab
    tab.pending && _tan_table_stale()
    _tan_copy!(tab.J, J)
    tab.noise && _tan_copy!(tab.Q, Q)
    tab.count = 0
    tab.valid = true
    _tan_series_lists!(tab)
    return tab
end

"""Extend `tab`'s powers to `count`, building only those it lacks."""
function _tan_table_extend!(tab::TanSeriesTable{W}, count::Int) where {W}
    count <= tab.count && return tab
    n = tab.n
    _tan_series_grow!(tab.Jp, count, n)
    tab.noise && _tan_series_grow!(tab.Lp, count + 1, n)
    tab.count == 0 && tab.noise && _tan_copy!(tab.Lp[1], tab.Q)
    for j in (tab.count + 1):count
        listed = tab.lists.on
        if j == 1
            _tan_copy!(tab.Jp[1], tab.J)
        elseif listed
            _tan_mul_list!(tab.Jp[j], tab.J, tab.Jp[j - 1], tab.lists.EE)
        else
            _tan_mul!(tab.Jp[j], tab.J, tab.Jp[j - 1])
        end
        if tab.noise
            # L^j(Q) = J L^{j-1}(Q) + its transpose.
            listed ? _tan_mul_list!(tab.tmp, tab.J, tab.Lp[j], tab.lists.EV) :
                _tan_mul!(tab.tmp, tab.J, tab.Lp[j])
            _tan_comb!(tab.Lp[j + 1], ((1.0, tab.tmp),))
            _tan_add!(tab.Lp[j + 1], 1.0, tab.tmp, true, Val(:T))
        end
    end
    tab.count = count
    return tab
end

"""`Phi(h)` into `sb.out` from `tab`'s powers (`_series_intercept_table!`):
`Phi(h0) = sum_{j=0}^m h0^{j+1} / (j+1)! J^j`, `e^{J h0} =
sum_{i=0}^{m+1} h0^i / i! J^i`, then the doublings, on Horner's plan."""
function _tan_series_intercept_table!(sb::TanSeries{W}, tab::TanSeriesTable{W},
    h::Float64) where {W}
    _tan_series_plan!(sb, tab.J, h, 1.0)
    s, m, h0 = sb.s, sb.m, sb.h0
    _tan_table_extend!(tab, s > 0 ? m + 1 : m)
    out = sb.out
    first = s == 0 ? out : sb.Ps[1]
    listed = tab.lists.on
    _tan_zero!(first)
    _tan_diag!(first, h0)
    c = h0
    for j in 1:m
        c *= h0 / (j + 1)
        listed ? _tan_axpy_list!(first, c, tab.Jp[j], tab.lists.eE) :
            _tan_comb!(first, ((c, tab.Jp[j]),), true)
    end
    s == 0 && return out
    E0 = sb.Es[1]
    _tan_zero!(E0)
    _tan_diag!(E0, 1.0)
    e = 1.0
    for i in 1:(m + 1)
        e *= h0 / i
        listed ? _tan_axpy_list!(E0, e, tab.Jp[i], tab.lists.eE) :
            _tan_comb!(E0, ((e, tab.Jp[i]),), true)
    end
    return _tan_series_doublings!(sb, out, false, tab.lists)
end

"""`V(h)` into `sb.out` from `tab`'s powers (`_series_noise_table!`):
`V(h0) = sum_{j=0}^m h0^{j+1} / (j+1)! L^j(Q)`, `e^{J h0} =
sum_{j=0}^m h0^j / j! J^j`, then the doublings; made exactly symmetric."""
function _tan_series_noise_table!(sb::TanSeries{W}, tab::TanSeriesTable{W},
    h::Float64) where {W}
    _tan_series_plan!(sb, tab.J, h, 2.0)
    s, m, h0 = sb.s, sb.m, sb.h0
    _tan_table_extend!(tab, m)
    out = sb.out
    first = s == 0 ? out : sb.Ps[1]
    listed = tab.lists.on
    c = h0
    _tan_comb!(first, ((c, tab.Lp[1]),))
    for j in 1:m
        c *= h0 / (j + 1)
        listed ? _tan_axpy_list!(first, c, tab.Lp[j + 1], tab.lists.eV) :
            _tan_comb!(first, ((c, tab.Lp[j + 1]),), true)
    end
    if s > 0
        E0 = sb.Es[1]
        _tan_zero!(E0)
        _tan_diag!(E0, 1.0)
        e = 1.0
        for j in 1:m
            e *= h0 / j
            listed ? _tan_axpy_list!(E0, e, tab.Jp[j], tab.lists.eE) :
                _tan_comb!(E0, ((e, tab.Jp[j]),), true)
        end
        _tan_series_doublings!(sb, out, true, tab.lists)
    end
    _tan_sym!(sb.b1, out)
    _tan_copy!(out, sb.b1)
    return out
end

"""Grow `tab`'s cotangents to powers up to `count`, those first reached zeroed;
a release left none pending, which is a fresh start."""
function _tan_table_bar_grow!(tab::TanSeriesTable{W}, count::Int) where {W}
    n = tab.n
    _tan_series_grow!(tab.Jpbar, count, n)
    tab.noise && _tan_series_grow!(tab.Lpbar, count + 1, n)
    if !tab.pending
        tab.bcount = 0
        tab.noise && _tan_zero!(tab.Lpbar[1])
    end
    for j in (tab.bcount + 1):count
        _tan_zero!(tab.Jpbar[j])
        tab.noise && _tan_zero!(tab.Lpbar[j + 1])
    end
    tab.bcount = max(tab.bcount, count)
    tab.pending = true
    return tab
end

"""Add to `tab`'s cotangents, times `w`, those of the interval
`_tan_series_intercept_table!` last computed in `sb`, for `Phibar`
(`_series_intercept_table_pullback!`)."""
function _tan_series_intercept_table_pullback!(tab::TanSeriesTable{W}, sb::TanSeries{W},
    Phibar::TanMat{W}, w::Float64, Abar::Union{Nothing,TanMat{W}}=nothing) where {W}
    s, m, h0 = sb.s, sb.m, sb.h0
    Pb, Eb = _tan_series_intercept_doubling_pullback!(sb, Phibar, tab.lists, Abar)
    # The exponential's own series is read where the chain doubles, or where
    # `e^{J h}` was handed out with no doubling (`_tan_series_exp!`).
    expo = s > 0 || Abar !== nothing
    _tan_table_bar_grow!(tab, expo ? m + 1 : m)
    c = h0
    for j in 1:m
        c *= h0 / (j + 1)
        _tan_comb!(tab.Jpbar[j], ((w * c, Pb),), true)
    end
    if expo
        e = 1.0
        for i in 1:(m + 1)
            e *= h0 / i
            _tan_comb!(tab.Jpbar[i], ((w * e, Eb),), true)
        end
    end
    return tab
end

"""
`e^{J h}` for the interval `_tan_series_intercept_table!` last computed in
`sb`, off its chain: the square of the last doubling's exponential, or with no
doubling the series itself. A prediction on the series takes it from here
(`_tan_predict!`) rather than as a Pade approximant and its Frechet derivative
beside the chain: on the affect pilot's tied model those were a quarter of a
Laplace gradient, the same squarings a third time (2026-11-18, dev2). Its
cotangent goes back with the intercept's (`_tan_series_intercept_table_pullback!`).
"""
function _tan_series_exp!(sb::TanSeries{W}, tab::TanSeriesTable{W}) where {W}
    A = sb.A
    s, m, h0 = sb.s, sb.m, sb.h0
    if s > 0
        E = sb.Es[s]
        tab.lists.on ? _tan_mul_list!(A, E, E, tab.lists.EE) : _tan_mul!(A, E, E)
        return A
    end
    _tan_table_extend!(tab, m + 1)
    _tan_zero!(A)
    _tan_diag!(A, 1.0)
    e = 1.0
    for i in 1:(m + 1)
        e *= h0 / i
        tab.lists.on ? _tan_axpy_list!(A, e, tab.Jp[i], tab.lists.eE) :
            _tan_comb!(A, ((e, tab.Jp[i]),), true)
    end
    return A
end

"""Add to `tab`'s cotangents, times `w`, those of the interval
`_tan_series_noise_table!` last computed in `sb`, for `Vbar`
(`_series_noise_table_pullback!`)."""
function _tan_series_noise_table_pullback!(tab::TanSeriesTable{W}, sb::TanSeries{W},
    Vbar::TanMat{W}, w::Float64) where {W}
    s, m, h0 = sb.s, sb.m, sb.h0
    Vb, Eb = _tan_series_noise_doubling_pullback!(sb, Vbar, tab.lists)
    _tan_table_bar_grow!(tab, m)
    c = h0
    _tan_comb!(tab.Lpbar[1], ((w * c, Vb),), true)
    for j in 1:m
        c *= h0 / (j + 1)
        _tan_comb!(tab.Lpbar[j + 1], ((w * c, Vb),), true)
    end
    if s > 0
        e = 1.0
        for j in 1:m
            e *= h0 / j
            _tan_comb!(tab.Jpbar[j], ((w * e, Eb),), true)
        end
    end
    return tab
end

"""Pull `tab`'s pending cotangents back through its power chains into
`tab.Jbar` (and `tab.Qbar`) and clear them (`_series_table_pullback!`):
`J^j = J J^{j-1}`, `L^j = J L^{j-1} + L^{j-1} J'`, `L^0 = Q`. `false` where none
were pending."""
function _tan_table_release!(tab::TanSeriesTable{W}) where {W}
    tab.pending || return false
    c, J, Jb = tab.bcount, tab.J, tab.Jbar
    _tan_zero!(Jb)
    L = tab.lists
    listed = L.on
    if tab.noise
        for j in c:-1:1
            Lb, A, prev = tab.Lpbar[j + 1], tab.Lp[j], tab.Lpbar[j]
            if listed
                _tan_mul_list!(Jb, Lb, A, L.DVt, 1.0, 1.0)
                _tan_mul_list!(Jb, Lb, A, L.DtV, 1.0, 1.0)
                _tan_mul_list!(prev, J, Lb, L.EtD, 1.0, 1.0)
                _tan_mul_list!(prev, Lb, J, L.DE, 1.0, 1.0)
            else
                _tan_mulNT!(Jb, Lb, A, true, true)
                _tan_mulTN!(Jb, Lb, A, true, true)
                _tan_mulTN!(prev, J, Lb, true, true)
                _tan_mul!(prev, Lb, J, true, true)
            end
        end
        _tan_copy!(tab.Qbar, tab.Lpbar[1])
    end
    for j in c:-1:2
        if listed
            _tan_mul_list!(Jb, tab.Jpbar[j], tab.Jp[j - 1], L.DEt, 1.0, 1.0)
            _tan_mul_list!(tab.Jpbar[j - 1], J, tab.Jpbar[j], L.EtD, 1.0, 1.0)
        else
            _tan_mulNT!(Jb, tab.Jpbar[j], tab.Jp[j - 1], true, true)
            _tan_mulTN!(tab.Jpbar[j - 1], J, tab.Jpbar[j], true, true)
        end
    end
    c >= 1 && _tan_comb!(Jb, ((1.0, tab.Jpbar[1]),), true)
    tab.pending = false
    tab.bcount = 0
    return true
end

"""The intercept's series over `dt` for `J` into `sb.out`: from `tab` where the
drift holds over the rows, by Horner's rule where it may change by row
(`drows`, `TanModel.drows`)."""
@inline _tan_intercept_series!(sb::TanSeries{W}, tab::TanSeriesTable{W}, J::TanMat{W},
    dt::Float64, drows::Bool) where {W} = drows ? _tan_series_intercept!(sb, J, dt) :
    _tan_series_intercept_table!(sb, _tan_table_for!(tab, J, nothing), dt)

"""The noise's series over `dt` for `J` and `Q` into `sb.out`, as
`_tan_intercept_series!`."""
@inline _tan_noise_series!(sb::TanSeries{W}, tab::TanSeriesTable{W}, J::TanMat{W},
    Q::TanMat{W}, dt::Float64, drows::Bool) where {W} = drows ?
    _tan_series_noise!(sb, J, Q, dt) : _tan_series_noise_table!(sb, _tan_table_for!(tab, J, Q), dt)
