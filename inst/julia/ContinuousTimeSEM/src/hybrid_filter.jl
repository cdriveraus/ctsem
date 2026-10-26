# Deviations the subject filter integrates itself (`CTSEMHybrid`): the
# forward recursion, the closing term, the tape records and their reverse.
#
# Forward, per subject pass, with `q` columns (the hybrid items the subject
# carries: a manifest intercept's where it observes the item at all, a state
# mean's always) and the filter's own gains and innovation covariances:
#
#   start:    X = 0, except an initial mean's column, e_l
#   update:   V = Lambda[o,:] X + E          (E selects each observed item's column)
#             S += V' S_t^-1 V,  s += V' S_t^-1 y_t      (through U^-T, never S^-1)
#             X -= K V
#   predict:  X = A X                         (A the discrete drift)
#             + J^-1 (A - I) e_l  for a continuous intercept's column, or
#             + e_l               for a discrete-time model's intercept
#   impulse:  X[l, c] += u_p                  a predictor's effect, where the
#                                             state takes TDPREDEFFECT u
#   close:    ll += 1/2 b' M^-1 b  [ - 1/2 logdet(M) ],   M = I + G' S G,  b = G' s
#
# with the deviations `G u`, `u ~ N(0, I)` of dimension `rd` (`_hybrid_factor!`):
# `G = diag(g)`, `g` each item's own sd times its cell's slope, for own
# variances; the slopes times the trailing block of a full-rank level's
# ordered covariance factor for a block (`CTSEMHybrid`). The columns are per
# unit of the cell's value, as a manifest intercept's `E` is, so the slope
# enters through `G` alone. For fixed
# deviations a linear Gaussian filter's covariances and gains do not depend on
# them, so the log likelihood is exactly quadratic in them, and this is its
# integral against their N(0, I) prior (de Jong's augmented filter). Nothing
# here propagates a covariance for the deviations: a column carries no dynamics
# but `A X`, and a continuous intercept's the discrete intercept's Jacobian in
# it, which the filter forms for the intercept itself (`_compute_discrete_time_form!`:
# the affine block's solve, or the series where its pivots are small).
#
# By default the pass returns the profile, `1/2 b' M^-1 b` without the logdet:
# the joint log density at the deviations' conditional mode. That is what
# reproduces Laplace over the random effects *and* the deviations exactly: the
# joint mode in `u` is the profile's, the profile's curvature is the Schur
# complement of the joint curvature, and the joint logdet is that one's plus
# `logdet(M)` -- which the Laplace layer adds at the mode
# (`_laplace_unit_loglik_gradient`'s `hlogdet`). The other modes are in
# `CTSEMHybridBuffers`.
#
# Every buffer here and on the tape is sized by capacity (items, manifests,
# states) and only grows, so after a model's first pass nothing allocates; the
# loops read the leading block a subject or row needs. Products and solves are
# loops: nothing takes a BLAS or LAPACK lock.

"""The number of deviations `u` the last pass integrated: zero where it
carried none."""
_hybrid_ndev(hb::CTSEMHybridBuffers) = hb.active && hb.q > 0 ? hb.rd : 0

"""
`out[base + a] (+)= b - (M - I) u`, the fixed deviations' own gradient, from
the pass just run in `_HYBRID_FIXED` (`u = hb.dfix`); `add` accumulates.
"""
function _hybrid_fixed_gradient!(out::AbstractVector, base::Integer,
    hb::CTSEMHybridBuffers; add::Bool = false)
    rd = _hybrid_ndev(hb)
    @inbounds for a in 1:rd
        acc = hb.bv[a] + hb.dfix[a]
        for b in 1:rd
            acc -= hb.Mh[a, b] * hb.dfix[b]
        end
        out[base + a] = add ? out[base + a] + acc : acc
    end
    return out
end

"""Start a pass: which hybrid items this subject carries, and their state."""
function _hybrid_start!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid, data,
    n::Int, m::Int) where {T}
    hb.active = _hybrid_on(h) && hb.mode != _HYBRID_OFF
    hb.q = 0
    hb.half_logdet = zero(T)
    hb.active || return nothing
    nh = _hybrid_items(h)
    length(hb.col) == m || (hb.col = zeros(Int, m))
    fill!(hb.col, 0)
    if length(hb.item) < nh
        hb.item = zeros(Int, nh)
        hb.uitem = zeros(Int, nh)
        hb.ckind = zeros(Int, nh)
        hb.clat = zeros(Int, nh)
        hb.cslope = zeros(Float64, nh)
        hb.cpred = zeros(Int, nh)
    end
    q = 0
    @inbounds for k in 1:nh
        kind = h.kind[k]
        if kind == 1
            # A manifest intercept's column, where the subject observes the item.
            j = h.target[k]
            seen = false
            for c in axes(data, 2)
                if _ctsem_observed(data[j, c])
                    seen = true
                    break
                end
            end
            seen || continue
            q += 1
            hb.col[j] = q
        else
            q += 1
        end
        hb.item[q] = k
        hb.ckind[q] = kind
        hb.clat[q] = kind == 1 ? 0 : h.target[k]
        hb.cslope[q] = h.slope[k]
        hb.cpred[q] = kind == 4 ? h.pred[k] : 0
    end
    hb.q = q
    q == 0 && return nothing
    if size(hb.X, 1) < n || size(hb.X, 2) < nh
        hb.X = zeros(T, n, nh)
        hb.Xt = zeros(T, n, nh)
        hb.Xinc = zeros(T, n, nh)
    end
    if size(hb.V, 1) < m || size(hb.V, 2) < nh
        hb.V = zeros(T, m, nh)
        hb.W = zeros(T, m, nh)
        hb.e = zeros(T, m)
    end
    if size(hb.S, 1) < nh
        hb.S = zeros(T, nh, nh)
        hb.s = zeros(T, nh)
        hb.Mh = zeros(T, nh, nh)
        hb.Mw = zeros(T, nh, nh)
        hb.gv = zeros(T, nh)
        hb.dgv = zeros(T, nh)
        hb.bv = zeros(T, nh)
        hb.dhat = zeros(T, nh)
        hb.r = zeros(T, nh)
        hb.G = zeros(T, nh, nh)
    end
    if h.block && size(hb.Lb, 1) < length(h.perm)
        hb.Lb = zeros(T, length(h.perm), length(h.perm))
        hb.lversion = -1
    end
    @inbounds for c in 1:q
        for r in 1:n
            hb.X[r, c] = zero(T)
            hb.Xinc[r, c] = zero(T)
        end
        # An initial mean's column: a unit of the cell's value shifts the first
        # state mean by one (the slope is in `G`, as for `E`).
        hb.ckind[c] == 2 && (hb.X[hb.clat[c], c] = one(T))
        for a in 1:q
            hb.S[a, c] = zero(T)
        end
        hb.s[c] = zero(T)
    end
    return nothing
end

"""
One measurement update's share: `V`, the accumulation into `S` and `s`, and
`X -= K V`. `L` is LAMBDA's observed rows (the mean path, as the innovation's),
`factor` the innovation covariance's Cholesky, `K` the gain, `y` the
innovation.
"""
function _hybrid_update!(hb::CTSEMHybridBuffers{T}, L, factor, K, y,
    observed, n::Int) where {T}
    q = hb.q
    m = length(observed)
    X, V, W, S, s, e = hb.X, hb.V, hb.W, hb.S, hb.s, hb.e
    U = factor.U
    @inbounds for c in 1:q, i in 1:m
        acc = hb.col[observed[i]] == c ? one(T) : zero(T)
        for k in 1:n
            acc += L[i, k] * X[k, c]
        end
        V[i, c] = acc
    end
    # W = U^-T V and e = U^-T y, by forward substitution on U'.
    @inbounds for c in 1:q, i in 1:m
        acc = V[i, c]
        for k in 1:(i - 1)
            acc -= U[k, i] * W[k, c]
        end
        W[i, c] = acc / U[i, i]
    end
    @inbounds for i in 1:m
        acc = y[i]
        for k in 1:(i - 1)
            acc -= U[k, i] * e[k]
        end
        e[i] = acc / U[i, i]
    end
    _ctsem_mulTN!(_rs(S, q, q), _rs(W, m, q), _rs(W, m, q), true, true)
    @inbounds for b in 1:q
        acc = zero(T)
        for i in 1:m
            acc += W[i, b] * e[i]
        end
        s[b] += acc
    end
    _ctsem_mul!(_rs(X, n, q), K, _rs(V, m, q), -1, true)
    return nothing
end

"""
One prediction substep's share: `X = A X`, into the spare and swapped in, and
for a continuous intercept's column (`ckind` 3) the discrete intercept's
Jacobian in a unit of the cell's value, `J^-1 (A - I) e_l` over the affine
block `naff` by the engine's own solve, or `Phi e_l` where the substep took the
series (`series.intercept`, set by `_compute_discrete_time_form!` just
before, which left `Phi` in `series.out`: `series` is the intercept kernel's
buffer); in a discrete-time model, whose intercept is added as it stands
(`_compute_one_step_form!`), `e_l`. Kept in `Xinc` for the tape.
"""
function _hybrid_predict!(hb::CTSEMHybridBuffers{T}, A, JAx, series, naff::Int,
    dt, n::Int, continuous_time::Bool) where {T}
    q = hb.q
    X, Xt = hb.X, hb.Xt
    _ctsem_mul!(_rs(Xt, n, q), A, _rs(X, n, q))
    hb.X, hb.Xt = Xt, X
    X, Xinc = hb.X, hb.Xinc
    phi = nothing
    # Every continuous intercept's column solves against the same `JAx`: it is
    # factored once, at the first.
    factored = false
    @inbounds for c in 1:q
        hb.ckind[c] == 3 || continue
        l = hb.clat[c]
        if !continuous_time
            for r in 1:n
                Xinc[r, c] = r == l ? one(T) : zero(T)
                X[r, c] += Xinc[r, c]
            end
            continue
        end
        l <= naff || throw(ArgumentError(
            "a continuous intercept's deviation on latent $(l) lies outside the affine block of $(naff)"))
        if series.intercept
            phi === nothing && (phi = _series_block(series.out, naff))
            for r in 1:naff
                Xinc[r, c] = phi[r, l]
            end
        else
            if !factored
                if size(hb.Jc, 1) < naff
                    hb.Jc = zeros(T, naff, naff)
                    hb.jv = zeros(T, naff)
                    hb.piv = zeros(Int, naff)
                end
                for j in 1:naff, i in 1:naff
                    hb.Jc[i, j] = JAx[i, j]
                end
                _lu_factor_generic!(hb.Jc, hb.piv, naff)
                factored = true
            end
            jv = hb.jv
            for i in 1:naff
                jv[i] = A[i, l] - (i == l ? one(T) : zero(T))
            end
            _lu_solve_generic!(hb.Jc, hb.piv, jv, naff)
            for r in 1:naff
                Xinc[r, c] = jv[r]
            end
        end
        for r in (naff + 1):n
            Xinc[r, c] = zero(T)
        end
        for r in 1:n
            X[r, c] += Xinc[r, c]
        end
    end
    return nothing
end

"""
A time dependent predictor's impulse share, where the state takes
`TDPREDEFFECT u` (`_apply_td_impulse!`): a predictor effect's column gains the
predictor's value at its latent, in a unit of the cell's value as every column
is. The impulse's own Jacobian is the identity in a linear model, the only kind
the filter integrates for, so nothing else moves; and its reverse is the
identity on `X`'s cotangent, which needs no record.
"""
function _hybrid_td_impulse!(hb::CTSEMHybridBuffers, tdpreds)
    @inbounds for c in 1:hb.q
        hb.ckind[c] == 4 || continue
        hb.X[hb.clat[c], c] += tdpreds[hb.cpred[c]]
    end
    return nothing
end

"""
The columns' factor `G` (`q x rd`) at `values`, into `hb.G`, and `rd`, the
dimension of `u`: for own variances the diagonal of the items' sds times their
slopes (`gv`, and `dgv` the derivative in the raw own sd); for a block the
slopes times the rows of the subject's columns in the trailing block of the
level's covariance factor in `h.perm` (`_hybrid_block_factor!`), every column
kept. `rd` is 0 where a block's covariance does not factor.
"""
function _hybrid_factor!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid,
    values::AbstractVector) where {T}
    q = hb.q
    G = hb.G
    if !h.block
        @inbounds for c in 1:q
            k = hb.item[c]
            raw = values[h.psi_index[k]]
            g = h.slope[k] * (log1p_exp(2 * raw - 1) * h.sd_scale[k] + 1e-10)
            hb.gv[c] = g
            hb.uitem[c] = k
            # d log1p_exp(z)/dz is the logistic function.
            hb.dgv[c] = h.slope[k] * h.sd_scale[k] * 2 / (1 + exp(-(2 * raw - 1)))
            for a in 1:q
                G[a, c] = a == c ? g : zero(T)
            end
        end
        hb.rd = q
        return q
    end
    nh = _hybrid_items(h)
    _hybrid_block_cached!(hb, h, values) || (hb.rd = 0; return 0)
    nL = length(h.perm) - nh
    Lb = hb.Lb
    @inbounds for f in 1:nh
        hb.uitem[f] = f
    end
    # `_ctsem_cholesky!` leaves the upper factor: the lower one's [r, f] is
    # Lb[f, r] for r >= f, and zero above.
    @inbounds for c in 1:q
        j = hb.item[c]
        sl = T(h.slope[j])
        for f in 1:nh
            G[c, f] = f <= j ? sl * Lb[nL + f, nL + j] : zero(T)
        end
    end
    hb.rd = nh
    return nh
end

"""
Whether two numbers are the same value with the same partials, at every depth
of dual: a cached quantity built from one is the other's. `==` on duals reads
the values alone.
"""
_hybrid_same(a::T, b::T) where {T<:AbstractFloat} = a == b
function _hybrid_same(a::ForwardDiff.Dual{G,V,N}, b::ForwardDiff.Dual{G,V,N}) where {G,V,N}
    _hybrid_same(ForwardDiff.value(a), ForwardDiff.value(b)) || return false
    @inbounds for i in 1:N
        _hybrid_same(ForwardDiff.partials(a, i), ForwardDiff.partials(b, i)) || return false
    end
    return true
end
_hybrid_same(a, b) = false

"""
The block's factor into `hb.Lb` (`_hybrid_block_factor!`), kept while the
level's parameters are the ones it was built at, value and partials alike, and
the spec is the one it was built for: a level's parameters are the same for
every subject a task filters in one pass, so the factor is built once a pass
rather than once a subject -- per subject it was most of a block's cost, and
63 MB an evaluation on AnomAuth's 800. Returns whether it factored.
"""
function _hybrid_block_cached!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid,
    values::AbstractVector) where {T}
    np = length(h.lparams)
    if hb.lversion == h.version && length(hb.lkey) == np
        same = true
        @inbounds for j in 1:np
            if !_hybrid_same(hb.lkey[j], values[h.lparams[j]])
                same = false
                break
            end
        end
        same && return hb.lok
    end
    length(hb.lkey) == np || (hb.lkey = Vector{T}(undef, np))
    @inbounds for j in 1:np
        hb.lkey[j] = values[h.lparams[j]]
    end
    hb.lversion = h.version
    hb.lok = _hybrid_block_factor!(hb.Lb, h, values)
    return hb.lok
end

"""
A block's level covariance at `values` in the order `h.perm`, factored in
place into `Lb` (`_ctsem_cholesky!`, the upper factor); false where it does
not factor. The construction is the layer's own (`_level_popcov`).
"""
function _hybrid_block_factor!(Lb::AbstractMatrix, h::CTSEMHybrid, values::AbstractVector)
    Sigma = _level_popcov(values, h.lv_sd, h.lv_cor, h.lv_sdfix, h.lv_corfix,
        h.lv_scale, h.lv_code)
    k = length(h.perm)
    perm = h.perm
    @inbounds for j in 1:k, i in 1:k
        Lb[i, j] = Sigma[perm[i], perm[j]]
    end
    return _ctsem_cholesky!(Lb, k)
end

"""
`M = I + G' S G` into `M` (`rd x rd`, both triangles) and `b = G' s` into `b`
over the subject's `q` columns, through `GS = G' S` (`rd x q`, scratch): the
closing term's matrices, which the pass, its reverse and the sequential
mixture (laplace_mixture.jl) all form. `GS` once, then `M` from it, is
`O(rd q^2 + rd^2 q)` where `G' S` per column of `M` was `O(rd^2 q^2)`.
"""
function _hybrid_closing!(M::AbstractMatrix, b::AbstractVector, GS::AbstractMatrix,
    G::AbstractMatrix, S::AbstractMatrix, s::AbstractVector, q::Int, rd::Int)
    T = eltype(M)
    @inbounds for e in 1:q, a in 1:rd
        acc = zero(T)
        for c in 1:q
            acc += G[c, a] * S[c, e]
        end
        GS[a, e] = acc
    end
    @inbounds for j in 1:rd, a in 1:j
        acc = a == j ? one(T) : zero(T)
        for e in 1:q
            acc += GS[a, e] * G[e, j]
        end
        M[a, j] = acc
        M[j, a] = acc
    end
    @inbounds for a in 1:rd
        acc = zero(T)
        for c in 1:q
            acc += G[c, a] * s[c]
        end
        b[a] = acc
    end
    return nothing
end

"""
The closing term, by `hb.mode`: the profile `1/2 b' M^-1 b`, less
`logdet(M) / 2` for the full integral, or at the caller's `u` (`hb.dfix`)
`u' b - 1/2 u' (M - I) u`; with `M = I + G' S G` and `b = G' s` over the `rd`
dimensions of `u` (`_hybrid_factor!`). What `CTSEMHybridBuffers` lists is left
in `hb` whatever the mode, and the term is recorded for the reverse when
`trace` is a tape. `values` is the raw parameter vector the filter was given,
which holds the own variances or the level's scales and correlations.
"""
function _hybrid_marginal!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid,
    values::AbstractVector, trace) where {T}
    q = hb.q
    if q == 0
        hb.rd = 0
        _record_hybrid!(trace, hb, h, values)
        return zero(T)
    end
    rd = _hybrid_factor!(hb, h, values)
    rd > 0 || return T(NaN)
    G, S, s = hb.G, hb.S, hb.s
    Mh, Mw, bv, z = hb.Mh, hb.Mw, hb.bv, hb.dhat
    # `M` twice: `Mw` is factored in place (`_ctsem_cholesky`), and the cross
    # score and the fixed-deviation value read `M` itself, from `Mh`. `Mw` is
    # the closing's scratch until it takes its copy.
    _hybrid_closing!(Mh, bv, Mw, G, S, s, q, rd)
    @inbounds for j in 1:rd, a in 1:rd
        Mw[a, j] = Mh[a, j]
    end
    F = _ctsem_cholesky(Mw, rd)
    issuccess(F) || return T(NaN)
    ldiv!(z, F, bv)
    half_logdet = zero(T)
    quad = zero(T)
    @inbounds for a in 1:rd
        half_logdet += log(Mw[a, a])
        quad += bv[a] * z[a]
    end
    hb.half_logdet = half_logdet
    fixed = hb.mode == _HYBRID_FIXED
    fixed && length(hb.dfix) < rd && throw(ArgumentError(string(
        "the fixed deviations hold ", length(hb.dfix), " values for ", rd, " dimensions")))
    # `r = b - M u` at the `u` the pass is at -- the conditional mode, or the
    # fixed one -- held fixed in the outermost seed's derivative.
    dref = fixed ? hb.dfix : z
    @inbounds for a in 1:rd
        acc = bv[a]
        for b in 1:rd
            acc -= Mh[a, b] * _hybrid_constant(T, _hybrid_strip(dref[b]))
        end
        hb.r[a] = acc
    end
    _record_hybrid!(trace, hb, h, values)
    if fixed
        d = hb.dfix
        value = zero(T)
        @inbounds for a in 1:rd
            value += d[a] * bv[a]
            for b in 1:rd
                value -= d[a] * (Mh[a, b] - (a == b ? one(T) : zero(T))) * d[b] / 2
            end
        end
        return value
    end
    return (hb.mode == _HYBRID_FULL ? -half_logdet : zero(T)) + quad / 2
end

# One dual layer off: the conditional mode held fixed in the derivative the
# outermost seed takes, which is what makes `r`'s partials the cross curvature.
@inline _hybrid_strip(x::ForwardDiff.Dual) = ForwardDiff.value(x)
@inline _hybrid_strip(x) = x
# ... and lifted back to `T` as a constant of the outermost seed, so the
# product is between two `T`s: the seeded sweeps' two tags have no order, and
# multiplying a dual by its own inner layer is a `DualMismatchError` there.
@inline _hybrid_constant(::Type{T}, x) where {T} = convert(T, x)
@inline _hybrid_constant(::Type{ForwardDiff.Dual{G,V,N}}, x) where {G,V,N} =
    ForwardDiff.Dual{G,V,N}(convert(V, x), zero(ForwardDiff.Partials{N,V}))

# Tape records. Every hook has a fallback that does nothing, so the plain
# primal and ForwardDiff paths compile to what they did.
_hybrid_record_X!(record, hb, n::Int) = nothing
function _hybrid_record_X!(record::Union{CTSEMUpdateRecord{T},CTSEMPredictRecord{T}},
    hb::CTSEMHybridBuffers{T}, n::Int) where {T}
    q = hb.active ? hb.q : 0
    q == 0 && return nothing
    Xh = record.Xh
    if size(Xh, 1) != n || size(Xh, 2) < q
        Xh = zeros(T, n, size(hb.X, 2))
        record.Xh = Xh
    end
    X = hb.X
    @inbounds for c in 1:q, r in 1:n
        Xh[r, c] = X[r, c]
    end
    return nothing
end

"""A prediction record's `Xinc`: what the substep added to each column."""
_hybrid_record_Xinc!(record, hb, n::Int) = nothing
function _hybrid_record_Xinc!(record::CTSEMPredictRecord{T},
    hb::CTSEMHybridBuffers{T}, n::Int) where {T}
    q = hb.active ? hb.q : 0
    q == 0 && return nothing
    Xi = record.Xinc
    if size(Xi, 1) != n || size(Xi, 2) < q
        Xi = zeros(T, n, size(hb.Xinc, 2))
        record.Xinc = Xi
    end
    Xinc = hb.Xinc
    @inbounds for c in 1:q, r in 1:n
        Xi[r, c] = Xinc[r, c]
    end
    return nothing
end

_record_hybrid!(trace, hb, h, values) = nothing
function _record_hybrid!(tape::CTSEMAdjointTape{T}, hb::CTSEMHybridBuffers{T},
    h::CTSEMHybrid, values::AbstractVector) where {T}
    rec = tape.hyb
    q = hb.q
    rd = hb.rd
    rec.active = true
    rec.mode = hb.mode
    rec.q = q
    rec.rd = rd
    rec.block = h.block
    q == 0 && return _tape_push!(tape, :hybrid, 0)
    nh = length(hb.gv)
    if length(rec.g) < nh
        rec.g = zeros(T, nh)
        rec.dg = zeros(T, nh)
        rec.s = zeros(T, nh)
        rec.dfix = zeros(T, nh)
        rec.S = zeros(T, nh, nh)
        rec.G = zeros(T, nh, nh)
    end
    length(rec.item) >= q || resize!(rec.item, nh)
    length(rec.col) == length(hb.col) || resize!(rec.col, length(hb.col))
    copyto!(rec.col, hb.col)
    @inbounds for a in 1:q
        rec.item[a] = hb.item[a]
        rec.g[a] = hb.gv[a]
        rec.dg[a] = hb.dgv[a]
        rec.s[a] = hb.s[a]
        for b in 1:q
            rec.S[a, b] = hb.S[a, b]
        end
        for b in 1:rd
            rec.G[a, b] = hb.G[a, b]
        end
    end
    if rec.mode == _HYBRID_FIXED
        @inbounds for a in 1:rd
            rec.dfix[a] = hb.dfix[a]
        end
    end
    if h.block
        np = length(h.lparams)
        length(rec.lvals) == np || (rec.lvals = zeros(T, np))
        @inbounds for j in 1:np
            rec.lvals[j] = values[h.lparams[j]]
        end
    end
    _tape_push!(tape, :hybrid, 0)
    return nothing
end

"""
Reverse of the closing term, run first in the walk: the cotangents of `S` and
`s`, of `G` (`Gbar`), and from `Gbar` each item's own-sd cotangent into
`psibar` (own variances, indexed by hybrid item; only `G`'s diagonal moves
with them) or the level's free scales' and correlations' into `lpbar` (a
block, `_hybrid_block_pullback!`); and the carried cotangent of `X` started at
zero. With `ll = 1/2 b' M^-1 b - w/2 logdet M`, `M = I + G' S G`, `b = G' s`:
`Mbar = -(w M^-1 + beta beta')/2`, `beta = M^-1 b`, so `Sbar = G Mbar G'`,
`sbar = G beta`, `Gbar = 2 S G Mbar + s beta'`; at fixed `u`,
`ll = u' G' s - 1/2 u' G' S G u`, `Sbar = -(G u)(G u)'/2`, `sbar = G u`,
`Gbar = (s - S G u) u'`. `nh` is the model's hybrid item count, the capacity
every per-item array is sized by.
"""
function _reverse_hybrid_marginal!(rec::CTSEMHybridTape{T}, n::Int, nh::Int,
    h::CTSEMHybrid) where {T}
    q = rec.q
    rd = rec.rd
    length(rec.psibar) == nh || (rec.psibar = zeros(T, nh))
    fill!(rec.psibar, zero(T))
    np = length(h.lparams)
    length(rec.lpbar) == np || (rec.lpbar = zeros(T, np))
    fill!(rec.lpbar, zero(T))
    if size(rec.Sbar, 1) < nh || size(rec.Xbar, 1) != n
        rec.Sbar = zeros(T, nh, nh)
        rec.sbar = zeros(T, nh)
        rec.Xbar = zeros(T, n, nh)
        rec.Mw = zeros(T, nh, nh)
        rec.C = zeros(T, nh, nh)
        rec.wv = zeros(T, nh)
        rec.uv = zeros(T, nh)
        rec.beta = zeros(T, nh)
        rec.Gbar = zeros(T, nh, nh)
    end
    q == 0 && return nothing
    G, S, s = rec.G, rec.S, rec.s
    Sbar, sbar, Xbar, Gbar = rec.Sbar, rec.sbar, rec.Xbar, rec.Gbar
    @inbounds for c in 1:q
        for r in 1:n
            Xbar[r, c] = zero(T)
        end
    end
    if rec.mode == _HYBRID_FIXED
        d = rec.dfix
        gd, sgd = rec.wv, rec.beta
        @inbounds for c in 1:q
            acc = zero(T)
            for a in 1:rd
                acc += G[c, a] * d[a]
            end
            gd[c] = acc
        end
        @inbounds for c in 1:q
            acc = zero(T)
            for e in 1:q
                acc += S[c, e] * gd[e]
            end
            sgd[c] = acc
        end
        @inbounds for e in 1:q, c in 1:q
            Sbar[c, e] = -gd[c] * gd[e] / 2
        end
        @inbounds for c in 1:q
            sbar[c] = gd[c]
            for a in 1:rd
                Gbar[c, a] = (s[c] - sgd[c]) * d[a]
            end
        end
    else
        Mw, C, col, unit, beta = rec.Mw, rec.C, rec.wv, rec.uv, rec.beta
        # `C` is the closing's scratch until it takes `M^-1`, and `col` holds
        # `b` until beta = M^-1 b is solved.
        _hybrid_closing!(Mw, col, C, G, S, s, q, rd)
        F = _ctsem_cholesky(Mw, rd)
        ldiv!(beta, F, col)
        @inbounds for b in 1:rd
            for a in 1:rd
                unit[a] = a == b ? one(T) : zero(T)
            end
            ldiv!(col, F, unit)
            for a in 1:rd
                C[a, b] = col[a]
            end
        end
        # C becomes Mbar, then Mw (factored, no longer needed) holds G Mbar.
        w = rec.mode == _HYBRID_FULL ? one(T) : zero(T)
        @inbounds for b in 1:rd, a in 1:rd
            C[a, b] = -(w * C[a, b] + beta[a] * beta[b]) / 2
        end
        GM = Mw
        @inbounds for a in 1:rd, c in 1:q
            acc = zero(T)
            for b in 1:rd
                acc += G[c, b] * C[b, a]
            end
            GM[c, a] = acc
        end
        @inbounds for e in 1:q, c in 1:q
            acc = zero(T)
            for a in 1:rd
                acc += GM[c, a] * G[e, a]
            end
            Sbar[c, e] = acc
        end
        @inbounds for c in 1:q
            acc = zero(T)
            for a in 1:rd
                acc += G[c, a] * beta[a]
            end
            sbar[c] = acc
            for a in 1:rd
                t = zero(T)
                for e in 1:q
                    t += S[c, e] * GM[e, a]
                end
                Gbar[c, a] = 2 * t + s[c] * beta[a]
            end
        end
    end
    if rec.block
        _hybrid_block_pullback!(rec, h, Gbar, q, rd)
    else
        @inbounds for c in 1:q
            rec.psibar[rec.item[c]] += Gbar[c, c] * rec.dg[c]
        end
    end
    return nothing
end

# A tag of its own for the block factor's derivative in each level parameter,
# one at a time: one dual type whatever the model (CLAUDE.md, "A model's types
# do not depend on the model").
struct _HybridFactorTag end

"""
`lpbar[j] = sum(Gbar .* dG/dtheta_j)` over the level's free scales and
correlations (`h.lparams`, raw values `rec.lvals`), `G` a block's factor
(`_hybrid_factor!`) for the subject's columns `rec.item[1:q]`, from the
factor's derivatives (`_hybrid_block_derivatives!`).
"""
function _hybrid_block_pullback!(rec::CTSEMHybridTape{T}, h::CTSEMHybrid,
    Gbar::AbstractMatrix{T}, q::Int, rd::Int) where {T}
    np = length(h.lparams)
    np == 0 && return rec.lpbar
    _hybrid_block_derivatives!(rec, h)
    nL = length(h.perm) - rd
    dLb = rec.dLb
    @inbounds for j in 1:np
        rec.dok[j] || continue
        acc = zero(T)
        for c in 1:q
            it = rec.item[c]
            sl = T(h.slope[it])
            for f in 1:min(it, rd)
                acc += Gbar[c, f] * sl * dLb[nL + f, nL + it, j]
            end
        end
        rec.lpbar[j] = acc
    end
    return rec.lpbar
end

"""
The block factor's derivative in each level parameter into `rec.dLb[:, :, j]`
(`rec.dok[j]` false where it does not factor), the factor rebuilt at a dual of
width one per parameter -- one tag, so no dual type depends on the model --
and kept while the level parameters are the values they were taken at
(`_hybrid_same`), as the forward factor is.
"""
function _hybrid_block_derivatives!(rec::CTSEMHybridTape{T}, h::CTSEMHybrid) where {T}
    np = length(h.lparams)
    lvals = rec.lvals
    if rec.dversion == h.version && length(rec.dkey) == np
        same = true
        @inbounds for j in 1:np
            if !_hybrid_same(rec.dkey[j], lvals[j])
                same = false
                break
            end
        end
        same && return nothing
    end
    k = length(h.perm)
    D = ForwardDiff.Dual{_HybridFactorTag,T,1}
    size(rec.dLb) == (k, k, np) || (rec.dLb = zeros(T, k, k, np))
    length(rec.dok) == np || (rec.dok = fill(false, np))
    vd = zeros(D, maximum(h.lparams))
    Lb = zeros(D, k, k)
    for j in 1:np
        for t in 1:np
            vd[h.lparams[t]] = D(lvals[t], ForwardDiff.Partials{1,T}((t == j ? one(T) : zero(T),)))
        end
        rec.dok[j] = _hybrid_block_factor!(Lb, h, vd)
        rec.dok[j] || continue
        @inbounds for c in 1:k, r in 1:k
            rec.dLb[r, c, j] = ForwardDiff.partials(Lb[r, c])[1]
        end
    end
    length(rec.dkey) == np || (rec.dkey = Vector{T}(undef, np))
    copyto!(rec.dkey, lvals)
    rec.dversion = h.version
    return nothing
end

"""
One measurement update's reverse, before the base `_reverse_update!` uses its
intermediates: from the carried `Xbar` (the cotangent of `X` after the update)
and the closing term's `Sbar`, `sbar`, the cotangents of the gain (`Kbar`), of
the innovation covariance (`Rbar`), of the innovation (`ybar`) and of the
observed LAMBDA rows (`Lbar`), which the base reverse folds into its own
(`_hybrid_accumulate!`, over the leading block); and `Xbar` carried back to
before the update. `Ui` is `U^-1`, `K` the gain, `a` `S^-1 y`, all as the base
reverse recomputed them.
"""
function _hybrid_reverse_update!(rec::CTSEMHybridTape{T}, record::CTSEMUpdateRecord{T},
    Ui, K, a, n::Int) where {T}
    q = rec.q
    o = record.observed
    m = length(o)
    Xh = record.Xh
    L = record.Lambda
    nh = size(rec.Xbar, 2)
    if size(rec.V, 1) < m || size(rec.V, 2) < nh
        rows = max(m, size(rec.V, 1))
        rec.V = zeros(T, rows, nh)
        rec.W = zeros(T, rows, nh)
        rec.Z = zeros(T, rows, nh)
        rec.Vbar = zeros(T, rows, nh)
    end
    if size(rec.Kbar, 1) != n || size(rec.Kbar, 2) < m
        cols = max(m, size(rec.Kbar, 2))
        rec.Kbar = zeros(T, n, cols)
        rec.Rbar = zeros(T, cols, cols)
        rec.Lbar = zeros(T, cols, n)
        rec.ybar = zeros(T, cols)
    end
    size(rec.Xt) == size(rec.Xbar) || (rec.Xt = zeros(T, n, nh))
    V, W, Z, Vbar = rec.V, rec.W, rec.Z, rec.Vbar
    Kbar, Rbar, Lbar, ybar = rec.Kbar, rec.Rbar, rec.Lbar, rec.ybar
    Xbar, Sbar, sbar, Xt = rec.Xbar, rec.Sbar, rec.sbar, rec.Xt
    @inbounds for c in 1:q, i in 1:m
        acc = rec.col[o[i]] == c ? one(T) : zero(T)
        for k in 1:n
            acc += L[i, k] * Xh[k, c]
        end
        V[i, c] = acc
    end
    # W = Ui' V, Z = Ui W = S^-1 V; Ui is upper triangular.
    @inbounds for c in 1:q, i in 1:m
        acc = zero(T)
        for k in 1:i
            acc += Ui[k, i] * V[k, c]
        end
        W[i, c] = acc
    end
    @inbounds for c in 1:q, i in 1:m
        acc = zero(T)
        for k in i:m
            acc += Ui[i, k] * W[k, c]
        end
        Z[i, c] = acc
    end
    # X+ = X - K V:  Kbar = -Xbar V',  Vbar = -K' Xbar.
    _ctsem_mulNT!(_rs(Kbar, n, m), _rs(Xbar, n, q), _rs(V, m, q), -1)
    # `W` is spent: it holds `Z Sbar` from here.
    _ctsem_mul!(_rs(W, m, q), _rs(Z, m, q), _rs(Sbar, q, q))
    @inbounds for c in 1:q, i in 1:m
        acc = zero(T)
        for r in 1:n
            acc -= K[r, i] * Xbar[r, c]
        end
        # S += V' S^-1 V:  Vbar += 2 Z Sbar.  s += V' S^-1 y:  Vbar += a sbar'.
        Vbar[i, c] = acc + 2 * W[i, c] + a[i] * sbar[c]
    end
    # Rbar (the innovation covariance) = -Z Sbar Z' - sym(Z sbar a');
    # ybar = Z sbar.
    @inbounds for i in 1:m
        zs = zero(T)
        for c in 1:q
            zs += Z[i, c] * sbar[c]
        end
        ybar[i] = zs
    end
    @inbounds for j in 1:m, i in 1:m
        acc = zero(T)
        for c in 1:q
            acc -= W[i, c] * Z[j, c]
        end
        Rbar[i, j] = acc - (ybar[i] * a[j] + a[i] * ybar[j]) / 2
    end
    # V = L X + E:  Lbar = Vbar X',  Xbar(before) = Xbar(after) + L' Vbar.
    _ctsem_mulNT!(_rs(Lbar, m, n), _rs(Vbar, m, q), _rs(Xh, n, q))
    @inbounds for c in 1:q, k in 1:n
        acc = Xbar[k, c]
        for i in 1:m
            acc += L[i, k] * Vbar[i, c]
        end
        Xt[k, c] = acc
    end
    rec.Xbar, rec.Xt = Xt, Xbar
    return nothing
end

"""
One prediction substep's reverse: `X+ = A X + Xinc`, so `Abar += Xbar X'` into
the transition's cotangent and `Xbar = A' Xbar`; and first, while `Xbar` is
still the cotangent of `X+`, each continuous intercept's column, whose
increment `Xinc[:, c] = J^-1 (A - I) e_l` is reversed as the discrete
intercept is in `_reverse_predict!`: `sbar = J^-T Xbar[1:naff, c]`,
`JAx_bar -= sbar Xinc[1:naff, c]'`, `Abar[:, l] += sbar`; or on the
series route the series pullback of `Xbar[1:naff, c] e_l'`, which reads
no exponential. In a discrete-time model (`continuous_time` false) the
increment is `e_l` and reads no parameter. Scratch from `sc` (`av1`, `aa1`,
`piv`, `series`), used before the intercept's own reverse reuses them.
"""
function _hybrid_reverse_predict!(Abar, JAx_bar, A, JAx, rec::CTSEMHybridTape{T},
    h::CTSEMHybrid, record::CTSEMPredictRecord{T}, n::Int, naff::Int, sc,
    dt; continuous_time::Bool=true) where {T}
    q = rec.q
    q == 0 && return nothing
    Xh = record.Xh
    Xbar = rec.Xbar
    Xinc = record.Xinc
    series_done = false
    sbar = _rs(sc.av1, naff)
    aa1 = _rs(sc.aa1, naff, naff)
    # `JAx'` is factored once, at the first column that solves against it; the
    # series route, which uses `aa1` as scratch instead, holds for the whole
    # record.
    factored = false
    @inbounds for c in 1:q
        continuous_time || break
        k = rec.item[c]
        h.kind[k] == 3 || continue
        l = h.target[k]
        if record.series_intercept
            if !series_done
                _series_intercept!(sc.series, JAx, dt, naff)
                series_done = true
            end
            for j in 1:naff, i in 1:naff
                aa1[i, j] = j == l ? Xbar[i, c] : zero(T)
            end
            _series_intercept_pullback!(JAx_bar, sc.series, naff, aa1)
        else
            for i in 1:naff
                sbar[i] = Xbar[i, c]
            end
            if !factored
                for j in 1:naff, i in 1:naff
                    aa1[i, j] = JAx[j, i]
                end
                _lu_factor_generic!(aa1, sc.piv, naff)
                factored = true
            end
            _lu_solve_generic!(aa1, sc.piv, sbar, naff)
            for j in 1:naff, i in 1:naff
                JAx_bar[i, j] -= sbar[i] * Xinc[j, c]
            end
            for i in 1:naff
                Abar[i, l] += sbar[i]
            end
        end
    end
    size(rec.Xt) == size(Xbar) || (rec.Xt = zeros(T, size(Xbar)...))
    Xt = rec.Xt
    _ctsem_mulNT!(view(Abar, 1:n, 1:n), _rs(Xbar, n, q), _rs(Xh, n, q), true, true)
    _ctsem_mulTN!(_rs(Xt, n, q), view(A, 1:n, 1:n), _rs(Xbar, n, q))
    rec.Xbar, rec.Xt = Xt, Xbar
    return nothing
end

"""`dest += src` over `dest`'s extent: `src` is a capacity-sized buffer."""
function _hybrid_accumulate!(dest::AbstractMatrix, src::AbstractMatrix)
    @inbounds for j in axes(dest, 2), i in axes(dest, 1)
        dest[i, j] += src[i, j]
    end
    return dest
end
function _hybrid_accumulate!(dest::AbstractVector, src::AbstractVector)
    @inbounds for i in eachindex(dest)
        dest[i] += src[i]
    end
    return dest
end

"""Add the own-sd cotangents, or a block's level ones, to the raw gradient,
after the parameter layer."""
function _hybrid_scatter!(gradient::AbstractVector, rec::CTSEMHybridTape,
    h::CTSEMHybrid)
    rec.active || return gradient
    if rec.block
        @inbounds for j in eachindex(h.lparams)
            j <= length(rec.lpbar) || break
            gradient[h.lparams[j]] += rec.lpbar[j]
        end
        return gradient
    end
    @inbounds for k in eachindex(h.psi_index)
        k <= length(rec.psibar) || break
        gradient[h.psi_index[k]] += rec.psibar[k]
    end
    return gradient
end
