# A row's categorical observations by hand, as the filter applies them
# (`_ekf_categorical_rows!`) and its adjoint undoes them (`_reverse_binary!`):
# each an exact scalar conditioning on its linear predictor,
#
#     c = P λ,  b = λ'c (+ a count's σ²),  a = λ'x + μ,
#     (logZ, m, v) = moments(a, b, extras),  x += c m/b,  P -= c c' (1 - v/b)/b,
#
# `λ` the measurement Jacobian's row and `μ` the intercept linearised at the
# row's prior (`_ekf_categorical_offset`). The matrix algebra is the sweep's,
# with tangents; the moment function is the engine's own
# (`_binary_moment_partials!`), its values and its Jacobian `J` in `(a, b,` the
# extras the observation reads`)` taken at a dual of five partials, one an
# input, whose partials are `J`'s derivatives, contracted with the lanes here --
# where the filter at the Hessian's dual carried every lane through the rule --
# or in a jet sweep at the jet of the sweep's width, whose second parts are
# `J`'s too. The extras -- an ordinal row's thresholds, a binary row's
# asymptotes, a count's dispersion, a censored row's sd and response scale --
# are cells of the layer (`TanModel.XC`, `_tan_cat_each`), assembled for each
# row as `_ordinal_thresholds!` assembles them (`TanCatPlan`).

"""A scalar with `W` tangents and a jet's second part, as a `TanVec`'s entry."""
struct TanS{W}
    v::Float64
    d::NTuple{W,Float64}
    s::Float64
end
@inline TanS{W}(v::Real) where {W} = TanS{W}(Float64(v), _tan_zeros(Val(W)), 0.0)
@inline _ts_dd(a::NTuple{W,Float64}, b::NTuple{W,Float64}) where {W} =
    sum(ntuple(w -> a[w] * b[w], Val(W)))
@inline _ts_add(a::TanS{W}, b::TanS{W}) where {W} =
    TanS{W}(a.v + b.v, ntuple(w -> a.d[w] + b.d[w], Val(W)), a.s + b.s)
@inline _ts_sub(a::TanS{W}, b::TanS{W}) where {W} =
    TanS{W}(a.v - b.v, ntuple(w -> a.d[w] - b.d[w], Val(W)), a.s - b.s)
@inline _ts_scale(a::TanS{W}, k::Float64) where {W} =
    TanS{W}(k * a.v, _tan_times(k, a.d), k * a.s)
@inline _ts_mul(a::TanS{W}, b::TanS{W}) where {W} =
    TanS{W}(a.v * b.v, ntuple(w -> a.d[w] * b.v + a.v * b.d[w], Val(W)),
        a.s * b.v + a.v * b.s + 2 * _ts_dd(a.d, b.d))
"""`q = a / b`: `dq = (da - q db) / b`, `s(q) = (s(a) - q s(b) - 2 sum_w dq_w db_w) / b`."""
@inline function _ts_div(a::TanS{W}, b::TanS{W}) where {W}
    q = a.v / b.v
    d = ntuple(w -> (a.d[w] - q * b.d[w]) / b.v, Val(W))
    return TanS{W}(q, d, (a.s - q * b.s - 2 * _ts_dd(d, b.d)) / b.v)
end

# Entries of the sweep's quantities as scalars, and back. A quantity's
# tangents and second parts are stored zero where its flags say none, so they
# are read as they are; a write marks them live only where it brings some.
@inline _ts_get(x::TanVec{W}, i::Int) where {W} =
    TanS{W}(@inbounds(x.v[i]), _tan_col(x, i), x.jet ? @inbounds(x.s[i]) : 0.0)
@inline _ts_get(A::TanMat{W}, i::Int, j::Int) where {W} =
    TanS{W}(@inbounds(A.v[i, j]), _tan_col(A, i, j), A.jet ? @inbounds(A.s[i, j]) : 0.0)
@inline function _ts_sput!(X, k::Int, s::Float64)
    if X.snz[]
        @inbounds X.s[k] = s
    elseif !iszero(s)
        X.snz[] = true
        @inbounds X.s[k] = s
    end
    return nothing
end
@inline function _ts_set!(x::TanVec{W}, i::Int, a::TanS{W}) where {W}
    @inbounds x.v[i] = a.v
    if x.nz[] || any(!iszero, a.d)
        x.nz[] = true
        _tan_put!(x, i, a.d)
    end
    x.jet && _ts_sput!(x, i, a.s)
    return x
end
@inline function _ts_set!(A::TanMat{W}, i::Int, j::Int, a::TanS{W}) where {W}
    @inbounds A.v[i, j] = a.v
    if A.nz[] || any(!iszero, a.d)
        A.nz[] = true
        _tan_put!(A, i, j, a.d)
    end
    A.jet && _ts_sput!(A, i + size(A.v, 1) * (j - 1), a.s)
    return A
end

"""`acc + x'y`, summed in order from `acc`."""
@inline function _ts_dot(x::TanVec{W}, y::TanVec{W}, acc::TanS{W}) where {W}
    @inbounds for i in eachindex(x.v)
        acc = _ts_add(acc, _ts_mul(_ts_get(x, i), _ts_get(y, i)))
    end
    return acc
end
"""`x[i] += a`."""
@inline _ts_addat!(x::TanVec{W}, i::Int, a::TanS{W}) where {W} =
    _ts_set!(x, i, _ts_add(_ts_get(x, i), a))
"""`y += a x`."""
function _ts_vaxpy!(y::TanVec{W}, a::TanS{W}, x::TanVec{W}) where {W}
    @inbounds for i in eachindex(x.v)
        _ts_set!(y, i, _ts_add(_ts_get(y, i), _ts_mul(a, _ts_get(x, i))))
    end
    return y
end
"""`y = x a` (`y` not `x`)."""
function _ts_vscale!(y::TanVec{W}, x::TanVec{W}, a::TanS{W}) where {W}
    @inbounds for i in eachindex(x.v)
        _ts_set!(y, i, _ts_mul(_ts_get(x, i), a))
    end
    return y
end
"""`M[i, :] += k a x'`."""
function _ts_rowaxpy!(M::TanMat{W}, i::Int, k::Float64, a::TanS{W}, x::TanVec{W}) where {W}
    ka = _ts_scale(a, k)
    @inbounds for j in eachindex(x.v)
        _ts_set!(M, i, j, _ts_add(_ts_get(M, i, j), _ts_mul(ka, _ts_get(x, j))))
    end
    return M
end
"""`M[i, :] += k x'`."""
function _ts_rowadd!(M::TanMat{W}, i::Int, k::Float64, x::TanVec{W}) where {W}
    @inbounds for j in eachindex(x.v)
        _ts_set!(M, i, j, _ts_add(_ts_get(M, i, j), _ts_scale(_ts_get(x, j), k)))
    end
    return M
end

"""`v = A[i, :]`, with tangents and second parts."""
function _tan_rowvec!(v::TanVec{W}, A::TanMat{W}, i::Int) where {W}
    on, son = _tan_nz(A), _tan_snz(A)
    _tan_jet_check(v, son)
    @inbounds for j in eachindex(v.v)
        v.v[j] = A.v[i, j]
        on && _tan_put!(v, j, _tan_col(A, i, j))
        son && (v.s[j] = A.s[i, j])
    end
    on ? (v.nz[] = true) : _tan_setnz!(v, false)
    son ? (v.snz[] = true) : _tan_setsnz!(v, false)
    return v
end

"""
    _tan_cat_each(f, sp)

`f(k, cell, i)` for each cell the categorical extras read, `k` its place on
the layer's `XC`, `cell` its place in the cells (`all_params`), `i` the
manifest that reads it, in the order `_ordinal_thresholds!` reads them: an
ordinal row's THRESHOLDS cells, a binary row's ASYMPTOTES row, a count's and
a censored row's MANIFESTVAR diagonal (a standard deviation before the
covariance is assembled), and a censored row's RESPONSE cells. The number of
cells, none for a model without categorical indicators.
"""
function _tan_cat_each(f::F, sp) where {F}
    types = sp.manifesttype
    isempty(types) && return 0
    l = sp.parameter_axis
    k = 0
    @inbounds for i in eachindex(types)
        t = types[i]
        if t == CTSEM_OBS_ORDINAL
            s = _layout_index(l, :THRESHOLDS)
            s == 0 && continue
            ncat = i <= length(sp.ncategories) ? sp.ncategories[i] : 0
            nr = l.nrows[s]
            for j in 1:min(max(ncat - 1, 0), l.ncols[s])
                k += 1
                f(k, l.first[s] + (i - 1) + nr * (j - 1), i)
            end
        elseif t == CTSEM_OBS_BINARY
            r = i <= length(sp.asymptoterow) ? sp.asymptoterow[i] : 0
            s = _layout_index(l, :ASYMPTOTES)
            (r > 0 && s != 0) || continue
            nr = l.nrows[s]
            for j in 1:2
                k += 1
                f(k, l.first[s] + (r - 1) + nr * (j - 1), i)
            end
        elseif t == CTSEM_OBS_COUNT || t == CTSEM_OBS_CENSORED
            s = _layout_index(l, :MANIFESTVAR)
            k += 1
            f(k, l.first[s] + (i - 1) + l.nrows[s] * (i - 1), i)
            t == CTSEM_OBS_CENSORED || continue
            s = _layout_index(l, :RESPONSE)
            s == 0 && continue
            nr = l.nrows[s]
            for j in 1:2
                k += 1
                f(k, l.first[s] + (i - 1) + nr * (j - 1), i)
            end
        end
    end
    return k
end

"""The number of the categorical extras' cells (`_tan_cat_each`)."""
_tan_cat_ncells(sp) = _tan_cat_each((k, cell, i) -> nothing, sp)

"""
Each manifest's categorical kind (0 Gaussian) and where its extras are: its
first cell on `XC` (0-based), how many, how many extras they assemble into
(`_ordinal_thresholds!`: an ordinal row's cumulated thresholds, a censored
row's limits `lo`, `hi` before its cells), the cells in all and the most
extras a row has.
"""
struct TanCatPlan
    kind::Vector{Int}
    first::Vector{Int}
    ncell::Vector{Int}
    ne::Vector{Int}
    lo::Vector{Float64}
    hi::Vector{Float64}
    nxc::Int
    emax::Int
end
TanCatPlan() = TanCatPlan(Int[], Int[], Int[], Int[], Float64[], Float64[], 0, 0)

function _tan_cat_plan(sp, m::Int)
    types = sp.manifesttype
    any(t -> t != 0, types) || return TanCatPlan()
    kind, first, ncell, ne = zeros(Int, m), zeros(Int, m), zeros(Int, m), zeros(Int, m)
    lo, hi = fill(-Inf, m), fill(Inf, m)
    for i in 1:min(m, length(types))
        kind[i] = types[i]
    end
    nxc = _tan_cat_each(sp) do k, cell, i
        ncell[i] == 0 && (first[i] = k - 1)
        ncell[i] += 1
    end
    for i in 1:m
        if kind[i] == CTSEM_OBS_CENSORED
            ne[i] = 2 + ncell[i]
            i <= length(sp.censormin) && (lo[i] = sp.censormin[i])
            i <= length(sp.censormax) && (hi[i] = sp.censormax[i])
        else
            ne[i] = ncell[i]
        end
    end
    return TanCatPlan(kind, first, ncell, ne, lo, hi, nxc, maximum(ne; init=0))
end

"""
Manifest `i`'s extras into `tau` from the layer's `XC`, as
`_ordinal_thresholds!` assembles them; how many.
"""
function _tan_cat_tau!(tau::TanVec{W}, XC::TanVec{W}, plan::TanCatPlan, i::Int) where {W}
    f, e, kind = plan.first[i], plan.ne[i], plan.kind[i]
    if kind == CTSEM_OBS_CENSORED
        _ts_set!(tau, 1, TanS{W}(plan.lo[i]))
        _ts_set!(tau, 2, TanS{W}(plan.hi[i]))
        for j in 3:e
            _ts_set!(tau, j, _ts_get(XC, f + j - 2))
        end
    elseif kind == CTSEM_OBS_ORDINAL
        for j in 1:e
            _ts_set!(tau, j, j == 1 ? _ts_get(XC, f + 1) :
                _ts_add(_ts_get(tau, j - 1), _ts_get(XC, f + j)))
        end
    else
        for j in 1:e
            _ts_set!(tau, j, _ts_get(XC, f + j))
        end
    end
    return e
end

"""
The cotangents `tb` of manifest `i`'s `e` assembled extras onto their cells in
`XCbar`, times `w` (`_extras_cotangent!`): a cumulated threshold's to every
cell at or below it, an asymptote's and a censored row's to their own.
"""
function _tan_cat_extras_bar!(XCbar::TanVec{W}, plan::TanCatPlan, i::Int,
    tb::Vector{TanS{W}}, e::Int, w::Float64) where {W}
    f, kind = plan.first[i], plan.kind[i]
    if kind == CTSEM_OBS_ORDINAL
        running = TanS{W}(0.0)
        for j in e:-1:1
            running = _ts_add(running, tb[j])
            _ts_addat!(XCbar, f + j, _ts_scale(running, w))
        end
    elseif kind == CTSEM_OBS_CENSORED
        for j in 3:e
            _ts_addat!(XCbar, f + j - 2, _ts_scale(tb[j], w))
        end
    else
        for j in 1:e
            _ts_addat!(XCbar, f + j, _ts_scale(tb[j], w))
        end
    end
    return XCbar
end

struct _TanCatTag end
"""The dual the moment function runs at: one partial for each input an
observation reads, `(a, b)` and at most three extras."""
const _TanCatDual = ForwardDiff.Dual{ForwardDiff.Tag{_TanCatTag,Float64},Float64,5}

"""
One categorical observation as the forward applied it, for the reverse: its
manifest, kind and value, the extras it reads (`lo` on, `nr` of them), the
state and covariance it was applied to, `c = P λ`, `a`, `b`, the moments, and
their Jacobian in `(a, b,` the extras read`)` with its tangents and second
parts.
"""
mutable struct TanCatObs{W}
    i::Int
    kind::Int
    y::Float64
    lo::Int
    nr::Int
    x::TanVec{W}
    P::TanMat{W}
    c::TanVec{W}
    a::TanS{W}
    b::TanS{W}
    logZ::TanS{W}
    m::TanS{W}
    v::TanS{W}
    Jv::Matrix{Float64}     # 3 x 5
    Jd::Array{Float64,3}    # W x 3 x 5
    Js::Matrix{Float64}
end
TanCatObs(n::Int, ::Val{W}) where {W} = TanCatObs{W}(0, 0, 0.0, 1, 0, TanVec{W}(n),
    TanMat{W}(n, n), TanVec{W}(n), TanS{W}(0.0), TanS{W}(0.0), TanS{W}(0.0), TanS{W}(0.0),
    TanS{W}(0.0), zeros(3, 5), zeros(W, 3, 5), zeros(3, 5))

"""A row's categorical observations: how many, each, and the row's prior,
which their intercepts read, and its cotangent."""
mutable struct TanCatRec{W}
    nobs::Int
    obs::Vector{TanCatObs{W}}
    prior::TanVec{W}
    priorbar::TanVec{W}
end
TanCatRec(n::Int, ::Val{W}) where {W} = TanCatRec{W}(0, TanCatObs{W}[], TanVec{W}(n),
    TanVec{W}(n))

"""The dual of `_tan_cat_moments_lanes!`: its partials two lanes, the most a
unit's coordinate or two carry."""
const _TanCatDual2 = ForwardDiff.Dual{ForwardDiff.Tag{_TanCatTag,Float64},Float64,2}

"""A sweep's scratch for its categorical rows, the extras' grown to the model
at first use (`_tan_cat_scratch!`)."""
mutable struct TanCatScratch{W}
    lam::TanVec{W}; dl::TanVec{W}; cbar::TanVec{W}; v1::TanVec{W}; v2::TanVec{W}
    tau::TanVec{W}
    taud::Vector{_TanCatDual}
    Jd::Matrix{_TanCatDual}
    taud2::Vector{_TanCatDual2}
    Jd2::Matrix{_TanCatDual2}
    taub::Vector{TanS{W}}
    # A jet sweep's extras and Jacobian as jets, at the sweep's width and at
    # two (`_tan_cat_moments_jet!`).
    jet::Any
    jet2::Any
end
TanCatScratch(n::Int, ::Val{W}) where {W} = TanCatScratch{W}(TanVec{W}(n), TanVec{W}(n),
    TanVec{W}(n), TanVec{W}(n), TanVec{W}(n), TanVec{W}(0), _TanCatDual[],
    zeros(_TanCatDual, 3, 2), _TanCatDual2[], zeros(_TanCatDual2, 3, 2), TanS{W}[],
    nothing, nothing)

function _tan_cat_scratch!(cs::TanCatScratch{W}, emax::Int) where {W}
    if length(cs.tau.v) < emax
        cs.tau = TanVec{W}(emax)
        resize!(cs.taud, emax)
        cs.Jd = zeros(_TanCatDual, 3, 2 + emax)
        resize!(cs.taud2, emax)
        cs.Jd2 = zeros(_TanCatDual2, 3, 2 + emax)
        resize!(cs.taub, emax)
        cs.jet = nothing
        cs.jet2 = nothing
    end
    return cs
end

"""
The last lane along which observation `o`'s inputs move -- `a`, `b`, the
extras it reads (`tau`) -- 0 where none: a unit with one coordinate has one,
and its moments are cheaper along that lane alone than in every input.
"""
function _tan_cat_live(o::TanCatObs{W}, tau::TanVec{W}) where {W}
    a, b = o.a, o.b
    live = 0
    @inbounds for w in 1:W
        moves = !iszero(a.d[w]) || !iszero(b.d[w])
        for t in o.lo:(o.lo + o.nr - 1)
            moves |= !iszero(tau.d[w, t])
        end
        moves && (live = w)
    end
    return live
end

"""
`_tan_cat_moments!` where at most two lanes move the inputs
(`_tan_cat_live`): the moment function at `_TanCatDual2` seeded with those
lanes, whose partials are then the moments' and the Jacobian's tangents
themselves. Needs `W >= 2`.
"""
function _tan_cat_moments_lanes!(o::TanCatObs{W}, cs::TanCatScratch{W}, tau::TanVec{W},
    e::Int) where {W}
    D = _TanCatDual2
    p = 2 + o.nr
    a, b = o.a, o.b
    taud = cs.taud2
    @inbounds for t in 1:e
        taud[t] = o.lo <= t < o.lo + o.nr ?
            D(tau.v[t], ForwardDiff.Partials((tau.d[1, t], tau.d[2, t]))) : D(tau.v[t])
    end
    Jd = cs.Jd2
    fill!(Jd, zero(D))
    nodes, weights = _binary_rule()
    logZ, m, v, handled = _binary_moment_partials!(Jd,
        D(a.v, ForwardDiff.Partials((a.d[1], a.d[2]))),
        D(b.v, ForwardDiff.Partials((b.d[1], b.d[2]))), o.y, nodes, weights,
        view(taud, 1:e), o.kind)
    handled || return _TAN_REFUSED
    isfinite(ForwardDiff.value(logZ)) || return _TAN_INVALID
    lanes(x) = ntuple(w -> w <= 2 ? ForwardDiff.partials(x, w) : 0.0, Val(W))
    o.logZ = TanS{W}(ForwardDiff.value(logZ), lanes(logZ), 0.0)
    o.m = TanS{W}(ForwardDiff.value(m), lanes(m), 0.0)
    o.v = TanS{W}(ForwardDiff.value(v), lanes(v), 0.0)
    @inbounds for cc in 1:p, r in 1:3
        x = Jd[r, cc <= 2 ? cc : o.lo + cc - 1]
        o.Jv[r, cc] = ForwardDiff.value(x)
        for w in 1:W
            o.Jd[w, r, cc] = w <= 2 ? ForwardDiff.partials(x, w) : 0.0
        end
        o.Js[r, cc] = 0.0
    end
    return _TAN_OK
end

"""
Observation `o`'s moments and their Jacobian at its `a`, `b` and the row's
`e` assembled extras `tau`: the moment function at `_TanCatDual`, an input it
reads a partial, its moments' tangents along the lanes and its Jacobian's
contracted from the partials. A status code: `_TAN_REFUSED` where the
partials are not the rule's own (`_binary_moment_partials!`: a variance at
the floor, a mode solve that used its budget), `_TAN_INVALID` where the
observation has no support.
"""
function _tan_cat_moments!(o::TanCatObs{W}, cs::TanCatScratch{W}, tau::TanVec{W},
    e::Int) where {W}
    D = _TanCatDual
    p = 2 + o.nr
    taud = cs.taud
    @inbounds for t in 1:e
        q = t - o.lo + 3
        taud[t] = 3 <= q <= p ?
            D(tau.v[t], ForwardDiff.Partials(ntuple(j -> j == q ? 1.0 : 0.0, Val(5)))) :
            D(tau.v[t])
    end
    Jd = cs.Jd
    fill!(Jd, zero(D))
    nodes, weights = _binary_rule()
    ad = D(o.a.v, ForwardDiff.Partials(ntuple(j -> j == 1 ? 1.0 : 0.0, Val(5))))
    bd = D(o.b.v, ForwardDiff.Partials(ntuple(j -> j == 2 ? 1.0 : 0.0, Val(5))))
    logZ, m, v, handled = _binary_moment_partials!(Jd, ad, bd, o.y, nodes, weights,
        view(taud, 1:e), o.kind)
    handled || return _TAN_REFUSED
    isfinite(ForwardDiff.value(logZ)) || return _TAN_INVALID
    # The inputs' tangents: `a`'s, `b`'s, the extras read.
    a, b = o.a, o.b
    du(j::Int, w::Int) = j == 1 ? a.d[w] : j == 2 ? b.d[w] : @inbounds(tau.d[w, o.lo + j - 3])
    o.logZ = _tan_cat_contract(logZ, du, p, Val(W))
    o.m = _tan_cat_contract(m, du, p, Val(W))
    o.v = _tan_cat_contract(v, du, p, Val(W))
    @inbounds for cc in 1:p, r in 1:3
        x = Jd[r, cc <= 2 ? cc : o.lo + cc - 1]
        o.Jv[r, cc] = ForwardDiff.value(x)
        h = ForwardDiff.partials(x)
        for w in 1:W
            acc = 0.0
            for j in 1:p
                acc = muladd(h[j], du(j, w), acc)
            end
            o.Jd[w, r, cc] = acc
        end
        o.Js[r, cc] = 0.0
    end
    return _TAN_OK
end

"""A moment at the dual as a scalar: its value, and its partials contracted
with the inputs' tangents `du`."""
@inline function _tan_cat_contract(x::_TanCatDual, du::F, p::Int, ::Val{W}) where {F,W}
    h = ForwardDiff.partials(x)
    d = ntuple(Val(W)) do w
        acc = 0.0
        for j in 1:p
            acc = muladd(h[j], du(j, w), acc)
        end
        acc
    end
    return TanS{W}(ForwardDiff.value(x), d, 0.0)
end

"""
`_tan_cat_moments!` in a jet sweep: the moment function at a jet of `JW`
directions, the sweep's width or two where no lane past the second moves the
inputs, seeded with the inputs' tangents and second parts, whose moments and
Jacobian then carry theirs. Behind a barrier, so that only a sweep that runs
it compiles it.
"""
function _tan_cat_moments_jet!(o::TanCatObs{W}, cs::TanCatScratch{W}, tau::TanVec{W},
    e::Int, ::Val{JW}) where {W,JW}
    J = CTSEMTraceJet{JW}
    sc = JW == W ? cs.jet : cs.jet2
    if !(sc isa Tuple{Vector{J},Matrix{J}}) || length(sc[1]) < e
        sc = (zeros(J, max(e, 1)), zeros(J, 3, 2 + max(e, 1)))
        JW == W ? (cs.jet = sc) : (cs.jet2 = sc)
    end
    tauj, Jm = sc::Tuple{Vector{J},Matrix{J}}
    p = 2 + o.nr
    # The lanes past `JW` move nothing here (`_tan_cat_live`).
    jet(v, d::NTuple{W,Float64}, s) = J(v, ntuple(j -> d[j], Val(JW)), s)
    @inbounds for t in 1:e
        tauj[t] = o.lo <= t < o.lo + o.nr ?
            jet(tau.v[t], _tan_col(tau, t), tau.jet ? tau.s[t] : 0.0) : J(tau.v[t])
    end
    fill!(Jm, zero(J))
    nodes, weights = _binary_rule()
    a, b = o.a, o.b
    logZ, m, v, handled = _binary_moment_partials!(Jm, jet(a.v, a.d, a.s),
        jet(b.v, b.d, b.s), o.y, nodes, weights, view(tauj, 1:e), o.kind)
    handled || return _TAN_REFUSED
    isfinite(logZ.v) || return _TAN_INVALID
    lanes(x) = ntuple(w -> w <= JW ? x.d[w] : 0.0, Val(W))
    o.logZ = TanS{W}(logZ.v, lanes(logZ), logZ.s)
    o.m = TanS{W}(m.v, lanes(m), m.s)
    o.v = TanS{W}(v.v, lanes(v), v.s)
    @inbounds for cc in 1:p, r in 1:3
        x = Jm[r, cc <= 2 ? cc : o.lo + cc - 1]
        o.Jv[r, cc] = x.v
        for w in 1:W
            o.Jd[w, r, cc] = w <= JW ? x.d[w] : 0.0
        end
        o.Js[r, cc] = x.s
    end
    return _TAN_OK
end

"""
Row `t`'s categorical observations, in the order the filter applies them,
from the state `x` and covariance `P` (`sw.x`, `sw.P`), recorded in
`sw.cats[t]`; their log likelihood's share into `sw.ll`. A status code
(`_tan_cat_moments!`'s).
"""
function _tan_cat_row!(x::TanVec{W}, P::TanMat{W}, sw, mod, y::AbstractVector{Float64},
    t::Int) where {W}
    plan = mod.cat::TanCatPlan
    t > length(sw.cats) && push!(sw.cats, TanCatRec(length(x.v), Val(W)))
    rec = sw.cats[t]::TanCatRec{W}
    rec.nobs = 0
    kind = plan.kind
    any(k -> isfinite(y[k]) && kind[k] > 0, eachindex(y)) || return _TAN_OK
    cs = _tan_cat_scratch!(sw.cat::TanCatScratch{W}, plan.emax)
    _tan_copy!(rec.prior, x)
    for k in eachindex(y)
        (isfinite(y[k]) && kind[k] > 0) || continue
        st = _tan_cat_obs!(x, P, rec, cs, mod, k, y[k], sw.ll)
        st == _TAN_OK || return st
    end
    return _TAN_OK
end

"""λ = Jy[i, :] into `cs.lam`, and `LAMBDA[i, :] - λ`, which the intercept's
prior term reads, into `cs.dl`."""
function _tan_cat_loading!(cs::TanCatScratch, mod, i::Int)
    _tan_rowvec!(cs.lam, mod.Jy, i)
    _tan_rowvec!(cs.dl, mod.LAMBDA, i)
    _tan_lincomb!(cs.dl, true, cs.dl, -1.0, cs.lam)
    return nothing
end

"""Manifest `i`'s observation `yi` (`_ekf_binary_update!`), recorded as `rec`'s
next. A status code."""
function _tan_cat_obs!(x::TanVec{W}, P::TanMat{W}, rec::TanCatRec{W},
    cs::TanCatScratch{W}, mod, i::Int, yi::Float64, ll::Vector{Float64}) where {W}
    plan = mod.cat::TanCatPlan
    n = length(x.v)
    rec.nobs += 1
    rec.nobs > length(rec.obs) && push!(rec.obs, TanCatObs(n, Val(W)))
    o = rec.obs[rec.nobs]
    kind = plan.kind[i]
    o.i, o.kind, o.y = i, kind, yi
    _tan_cat_loading!(cs, mod, i)
    lam = cs.lam
    # The intercept, MANIFESTMEANS[i] + (LAMBDA[i, :] - λ)' prior.
    mu = _ts_dot(cs.dl, rec.prior, _ts_get(mod.MM, i))
    _tan_copy!(o.x, x)
    _tan_copy!(o.P, P)
    c = o.c
    _tan_mulvec!(c, P, lam)
    b = _ts_dot(lam, c, TanS{W}(0.0))
    a = _ts_dot(lam, x, mu)
    b.v < 0 && (b = TanS{W}(0.0))
    e = _tan_cat_tau!(cs.tau, mod.XC, plan, i)
    if kind == CTSEM_OBS_COUNT
        # A count's dispersion joins the predictor's variance
        # (`_count_dispersion`).
        sig = _ts_get(cs.tau, 1)
        b = _ts_add(b, _ts_mul(sig, sig))
    end
    o.a, o.b = a, b
    # At the floor the filter moves nothing and adds the likelihood at `a`.
    b.v > _CTSEM_MIN_VARIANCE[] || return _TAN_REFUSED
    read = _extras_read(yi, view(cs.tau.v, 1:e), kind)
    o.lo, o.nr = first(read), length(read)
    # Along at most two lanes, the moments at two partials; otherwise in every
    # input they read, or along the jet sweep's lanes.
    narrow = W >= 2 && _tan_cat_live(o, cs.tau) <= 2
    st = if x.jet
        narrow ? _ctsem_barrier(_tan_cat_moments_jet!, o, cs, cs.tau, e, Val(2))::Int :
            _ctsem_barrier(_tan_cat_moments_jet!, o, cs, cs.tau, e, Val(W))::Int
    else
        narrow ? _tan_cat_moments_lanes!(o, cs, cs.tau, e) : _tan_cat_moments!(o, cs, cs.tau, e)
    end
    st == _TAN_OK || return st
    ll[1] += o.logZ.v
    @inbounds for w in 1:W
        ll[1 + w] += o.logZ.d[w]
    end
    x.jet && (ll[2 + W] += o.logZ.s)
    # x += c m/b; P -= c c' (1 - v/b)/b, on the lower triangle and mirrored.
    shift = _ts_div(o.m, b)
    shrink = _ts_div(_ts_sub(TanS{W}(1.0), _ts_div(o.v, b)), b)
    @inbounds for k in 1:n
        _ts_set!(x, k, _ts_add(_ts_get(x, k), _ts_mul(_ts_get(c, k), shift)))
    end
    @inbounds for j in 1:n, k in j:n
        _ts_set!(P, k, j, _ts_sub(_ts_get(P, k, j),
            _ts_mul(_ts_mul(shrink, _ts_get(c, k)), _ts_get(c, j))))
    end
    _tan_lower_to_upper!(P)
    return _TAN_OK
end

"""
Undo row `t`'s categorical observations (`_reverse_binary!`), last first:
`xbar`, `Pbar` carried back, the measurement's and the extras' cotangents
into `acc` times `w`; the intercepts' share of the prior's cotangent joins
`xbar` once the row is walked back.
"""
function _tan_reverse_cat_row!(xbar::TanVec{W}, Pbar::TanMat{W}, sw, mod, acc, w::Float64,
    t::Int) where {W}
    rec = sw.cats[t]::TanCatRec{W}
    rec.nobs == 0 && return nothing
    cs = sw.cat::TanCatScratch{W}
    _tan_zero!(rec.priorbar)
    for j in rec.nobs:-1:1
        _tan_reverse_cat_obs!(xbar, Pbar, rec.obs[j], rec, cs, mod, acc, w)
    end
    _tan_lincomb!(xbar, true, xbar, true, rec.priorbar)
    return nothing
end

function _tan_reverse_cat_obs!(xbar::TanVec{W}, Pbar::TanMat{W}, o::TanCatObs{W},
    rec::TanCatRec{W}, cs::TanCatScratch{W}, mod, acc, w::Float64) where {W}
    plan = mod.cat::TanCatPlan
    i = o.i
    _tan_cat_loading!(cs, mod, i)
    lam, c, b, m, v = cs.lam, o.c, o.b, o.m, o.v
    one_ = TanS{W}(1.0)
    shift = _ts_div(m, b)
    shrink = _ts_div(_ts_sub(one_, _ts_div(v, b)), b)
    # The cotangents of the two things the step wrote: the shift, `x̄'c`, and
    # `c`'s, `x̄ shift - shrink (P̄ + P̄') c` (P̄ need not be symmetric:
    # `c = P λ` below adds `c̄ λ'`); the shrink's, `-c'P̄c`.
    shiftbar = _ts_dot(xbar, c, TanS{W}(0.0))
    cbar = cs.cbar
    @inbounds for k in eachindex(c.v)
        _ts_set!(cbar, k, _ts_mul(_ts_get(xbar, k), shift))
    end
    Pc, PTc = cs.v1, cs.v2
    _tan_mulvec!(Pc, Pbar, c)
    _tan_mulTvec!(PTc, Pbar, c)
    shrinkbar = _ts_scale(_ts_dot(c, Pc, TanS{W}(0.0)), -1.0)
    _tan_lincomb!(Pc, true, Pc, true, PTc)
    @inbounds for k in eachindex(c.v)
        _ts_set!(cbar, k, _ts_sub(_ts_get(cbar, k), _ts_mul(shrink, _ts_get(Pc, k))))
    end
    # Through shift = m/b and shrink = 1/b - v/b².
    bb = _ts_mul(b, b)
    mbar = _ts_div(shiftbar, b)
    bbar = _ts_scale(_ts_div(_ts_mul(shiftbar, m), bb), -1.0)
    vbar = _ts_scale(_ts_div(shrinkbar, bb), -1.0)
    bbar = _ts_add(bbar, _ts_mul(shrinkbar, _ts_add(_ts_scale(_ts_div(one_, bb), -1.0),
        _ts_div(_ts_scale(v, 2.0), _ts_mul(bb, b)))))
    # Through the moment function, the likelihood's cotangent one.
    Jg(r, cc) = @inbounds TanS{W}(o.Jv[r, cc], ntuple(q -> o.Jd[q, r, cc], Val(W)),
        o.Js[r, cc])
    abar = _ts_add(_ts_add(Jg(1, 1), _ts_mul(mbar, Jg(2, 1))), _ts_mul(vbar, Jg(3, 1)))
    bbar = _ts_add(bbar,
        _ts_add(_ts_add(Jg(1, 2), _ts_mul(mbar, Jg(2, 2))), _ts_mul(vbar, Jg(3, 2))))
    f = plan.first[i]
    if o.kind == CTSEM_OBS_COUNT
        # The dispersion enters through `b` alone: `db/dσ = 2σ`.
        sig = _ts_get(mod.XC, f + 1)
        _ts_addat!(acc.XC, f + 1, _ts_scale(_ts_mul(bbar, sig), 2.0 * w))
    elseif o.nr > 0
        e = plan.ne[i]
        tb = cs.taub
        @inbounds for t in 1:e
            tb[t] = TanS{W}(0.0)
        end
        @inbounds for q in 1:o.nr
            tb[o.lo + q - 1] = _ts_add(_ts_add(Jg(1, 2 + q), _ts_mul(mbar, Jg(2, 2 + q))),
                _ts_mul(vbar, Jg(3, 2 + q)))
        end
        _tan_cat_extras_bar!(acc.XC, plan, i, tb, e, w)
    end
    # a = λ'x + μ: the state, the direction, and the intercept's cells and
    # prior (`_intercept_cotangent!`).
    _ts_vaxpy!(xbar, abar, lam)
    _ts_rowaxpy!(acc.Jy, i, w, abar, o.x)
    _ts_addat!(acc.MM, i, _ts_scale(abar, w))
    _ts_rowaxpy!(acc.LAMBDA, i, w, abar, rec.prior)
    _ts_rowaxpy!(acc.Jy, i, -w, abar, rec.prior)
    _ts_vaxpy!(rec.priorbar, abar, cs.dl)
    # b = λ'P λ, directly; then c = P λ.
    _ts_rowaxpy!(acc.Jy, i, 2.0 * w, bbar, c)
    _ts_vscale!(PTc, lam, bbar)
    _tan_ger!(Pbar, 1.0, PTc, lam, true)
    _tan_ger!(Pbar, 1.0, cbar, lam, true)
    _tan_mulvec!(Pc, o.P, cbar)
    _ts_rowadd!(acc.Jy, i, w, Pc)
    return nothing
end
