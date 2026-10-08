# Deviations the subject filter integrates itself (`CTSEMHybrid`): the
# forward recursion, the closing term, the tape records and their reverse.
#
# Forward, per subject pass, with `q` columns (the hybrid items the subject
# observes at all) and the filter's own gains and innovation covariances:
#
#   update:   V = Lambda[o,:] X + E          (E selects each observed item's column)
#             S += V' S_t^-1 V,  s += V' S_t^-1 y_t      (through U^-T, never S^-1)
#             X -= K V
#   predict:  X = A X                         (A the discrete drift)
#   close:    ll += 1/2 b' M^-1 b  [ - 1/2 logdet(M) ],   M = I + G S G,  b = G s
#
# with `G = diag(g)`, `g` each item's own sd times its intercept's slope. For
# fixed deviations a linear Gaussian filter's covariances and gains do not
# depend on them, so the log likelihood is exactly quadratic in them, and this
# is its integral against their N(0, I) prior (de Jong's augmented filter).
# Nothing here propagates a covariance for the deviations: they carry no
# dynamics but `A X`, and no discretisation at all.
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

"""Start a pass: which hybrid items this subject observes, and zeroed state."""
function _hybrid_start!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid, data,
    n::Int, m::Int, continuous_time::Bool) where {T}
    hb.active = _hybrid_on(h)
    hb.q = 0
    hb.half_logdet = zero(T)
    hb.active || return nothing
    continuous_time || throw(ArgumentError(
        "deviations integrated in the filter need a continuous-time model"))
    nh = length(h.manifests)
    length(hb.col) == m || (hb.col = zeros(Int, m))
    fill!(hb.col, 0)
    length(hb.item) >= nh || (hb.item = zeros(Int, nh))
    q = 0
    @inbounds for (k, j) in enumerate(h.manifests)
        seen = false
        for c in axes(data, 2)
            if _ctsem_observed(data[j, c])
                seen = true
                break
            end
        end
        if seen
            q += 1
            hb.col[j] = q
            hb.item[q] = k
        end
    end
    hb.q = q
    q == 0 && return nothing
    if size(hb.X, 1) < n || size(hb.X, 2) < nh
        hb.X = zeros(T, n, nh)
        hb.Xt = zeros(T, n, nh)
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
    end
    @inbounds for c in 1:q
        for r in 1:n
            hb.X[r, c] = zero(T)
        end
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
    @inbounds for b in 1:q
        for a in 1:b
            acc = zero(T)
            for i in 1:m
                acc += W[i, a] * W[i, b]
            end
            S[a, b] += acc
            a == b || (S[b, a] += acc)
        end
        acc = zero(T)
        for i in 1:m
            acc += W[i, b] * e[i]
        end
        s[b] += acc
    end
    @inbounds for c in 1:q, r in 1:n
        acc = zero(T)
        for i in 1:m
            acc += K[r, i] * V[i, c]
        end
        X[r, c] -= acc
    end
    return nothing
end

"""One prediction substep's share: `X = A X`, into the spare and swapped in."""
function _hybrid_predict!(hb::CTSEMHybridBuffers{T}, A, n::Int) where {T}
    q = hb.q
    X, Xt = hb.X, hb.Xt
    @inbounds for c in 1:q, r in 1:n
        acc = zero(T)
        for k in 1:n
            acc += A[r, k] * X[k, c]
        end
        Xt[r, c] = acc
    end
    hb.X, hb.Xt = Xt, X
    return nothing
end

"""
The closing term, by `hb.mode`: the profile `1/2 b' M^-1 b`, less
`logdet(M) / 2` for the full integral, or at the caller's deviations `d`
(`hb.dfix`) `d' b - 1/2 d' (M - I) d`; with `M = I + G S G` and `b = G s`.
What `CTSEMHybridBuffers` lists is left in `hb` whatever the mode, and the term
is recorded for the reverse when `trace` is a tape. `values` is the raw
parameter vector the filter was given, which holds the own variances.
"""
function _hybrid_marginal!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid,
    values::AbstractVector, trace) where {T}
    q = hb.q
    if q == 0
        _record_hybrid!(trace, hb)
        return zero(T)
    end
    g, dg, Mh, Mw, bv, z = hb.gv, hb.dgv, hb.Mh, hb.Mw, hb.bv, hb.dhat
    @inbounds for c in 1:q
        k = hb.item[c]
        raw = values[h.psi_index[k]]
        g[c] = h.slope[k] * (log1p_exp(2 * raw - 1) * h.sd_scale[k] + 1e-10)
        # d log1p_exp(z)/dz is the logistic function.
        dg[c] = h.slope[k] * h.sd_scale[k] * 2 / (1 + exp(-(2 * raw - 1)))
    end
    # `M` twice: `Mw` is factored in place (`_ctsem_cholesky`), and the cross
    # score and the fixed-deviation value read `M` itself, from `Mh`.
    @inbounds for b in 1:q, a in 1:q
        v = g[a] * hb.S[a, b] * g[b] + (a == b ? one(T) : zero(T))
        Mh[a, b] = v
        Mw[a, b] = v
    end
    F = _ctsem_cholesky(Mw, q)
    issuccess(F) || return T(NaN)
    @inbounds for c in 1:q
        bv[c] = g[c] * hb.s[c]
    end
    ldiv!(z, F, bv)
    half_logdet = zero(T)
    quad = zero(T)
    @inbounds for c in 1:q
        half_logdet += log(Mw[c, c])
        quad += bv[c] * z[c]
    end
    hb.half_logdet = half_logdet
    fixed = hb.mode == _HYBRID_FIXED
    fixed && length(hb.dfix) < q && throw(ArgumentError(string(
        "the fixed deviations hold ", length(hb.dfix), " values for ", q, " columns")))
    # `r = b - M d` at the deviations the pass is at -- the conditional mode,
    # or the fixed ones -- held fixed in the outermost seed's derivative.
    dref = fixed ? hb.dfix : z
    @inbounds for a in 1:q
        acc = bv[a]
        for b in 1:q
            acc -= Mh[a, b] * _hybrid_constant(T, _hybrid_strip(dref[b]))
        end
        hb.r[a] = acc
    end
    _record_hybrid!(trace, hb)
    if fixed
        d = hb.dfix
        value = zero(T)
        @inbounds for a in 1:q
            value += d[a] * bv[a]
            for b in 1:q
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

_record_hybrid!(trace, hb) = nothing
function _record_hybrid!(tape::CTSEMAdjointTape{T}, hb::CTSEMHybridBuffers{T}) where {T}
    rec = tape.hyb
    q = hb.q
    rec.active = true
    rec.mode = hb.mode
    rec.q = q
    q == 0 && return _tape_push!(tape, :hybrid, 0)
    nh = length(hb.gv)
    if length(rec.g) < nh
        rec.g = zeros(T, nh)
        rec.dg = zeros(T, nh)
        rec.s = zeros(T, nh)
        rec.dfix = zeros(T, nh)
        rec.S = zeros(T, nh, nh)
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
    end
    if rec.mode == _HYBRID_FIXED
        @inbounds for a in 1:q
            rec.dfix[a] = hb.dfix[a]
        end
    end
    _tape_push!(tape, :hybrid, 0)
    return nothing
end

"""
Reverse of the closing term, run first in the walk: the cotangents of `S` and
`s`, each item's own-sd cotangent into `psibar` (indexed by hybrid item), and
the carried cotangent of `X` started at zero. `nh` is the model's hybrid item
count, the capacity every per-item array is sized by.
"""
function _reverse_hybrid_marginal!(rec::CTSEMHybridTape{T}, n::Int, nh::Int) where {T}
    q = rec.q
    length(rec.psibar) == nh || (rec.psibar = zeros(T, nh))
    fill!(rec.psibar, zero(T))
    if size(rec.Sbar, 1) < nh || size(rec.Xbar, 1) != n
        rec.Sbar = zeros(T, nh, nh)
        rec.sbar = zeros(T, nh)
        rec.Xbar = zeros(T, n, nh)
        rec.Mw = zeros(T, nh, nh)
        rec.C = zeros(T, nh, nh)
        rec.wv = zeros(T, nh)
        rec.uv = zeros(T, nh)
        rec.beta = zeros(T, nh)
    end
    q == 0 && return nothing
    g, S, s = rec.g, rec.S, rec.s
    Sbar, sbar, Xbar = rec.Sbar, rec.sbar, rec.Xbar
    @inbounds for c in 1:q
        for r in 1:n
            Xbar[r, c] = zero(T)
        end
    end
    if rec.mode == _HYBRID_FIXED
        # d' G s - 1/2 d' G S G d at the caller's `d`.
        d = rec.dfix
        @inbounds for b in 1:q, a in 1:q
            Sbar[a, b] = -g[a] * d[a] * d[b] * g[b] / 2
        end
        @inbounds for a in 1:q
            sbar[a] = g[a] * d[a]
            gbar = d[a] * s[a]
            for b in 1:q
                gbar -= d[a] * S[a, b] * g[b] * d[b]
            end
            rec.psibar[rec.item[a]] += gbar * rec.dg[a]
        end
        return nothing
    end
    Mw, C, col, unit, beta = rec.Mw, rec.C, rec.wv, rec.uv, rec.beta
    @inbounds for b in 1:q, a in 1:q
        Mw[a, b] = g[a] * S[a, b] * g[b] + (a == b ? one(T) : zero(T))
    end
    F = _ctsem_cholesky(Mw, q)
    @inbounds for b in 1:q
        for a in 1:q
            unit[a] = a == b ? one(T) : zero(T)
        end
        ldiv!(col, F, unit)
        for a in 1:q
            C[a, b] = col[a]
        end
    end
    @inbounds for a in 1:q
        acc = zero(T)
        for b in 1:q
            acc += C[a, b] * g[b] * s[b]
        end
        beta[a] = acc
    end
    # dll/dM = -1/2 (w C + beta beta'), dll/db = beta, with `w` 1 when the
    # pass returned the logdet and 0 for the profile.
    w = rec.mode == _HYBRID_FULL ? one(T) : zero(T)
    @inbounds for b in 1:q, a in 1:q
        Mbar = -(w * C[a, b] + beta[a] * beta[b]) / 2
        Sbar[a, b] = g[a] * Mbar * g[b]
    end
    @inbounds for a in 1:q
        sbar[a] = g[a] * beta[a]
        gbar = beta[a] * s[a]
        for b in 1:q
            gbar -= (w * C[a, b] + beta[a] * beta[b]) * S[a, b] * g[b]
        end
        rec.psibar[rec.item[a]] += gbar * rec.dg[a]
    end
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
    @inbounds for i in 1:m, r in 1:n
        acc = zero(T)
        for c in 1:q
            acc -= Xbar[r, c] * V[i, c]
        end
        Kbar[r, i] = acc
    end
    # `W` is spent: it holds `Z Sbar` from here.
    @inbounds for c in 1:q, i in 1:m
        acc = zero(T)
        for b in 1:q
            acc += Z[i, b] * Sbar[b, c]
        end
        W[i, c] = acc
    end
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
    @inbounds for k in 1:n, i in 1:m
        acc = zero(T)
        for c in 1:q
            acc += Vbar[i, c] * Xh[k, c]
        end
        Lbar[i, k] = acc
    end
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
One prediction substep's reverse: `X+ = A X`, so `Abar += Xbar X'` into the
transition's cotangent and `Xbar = A' Xbar`.
"""
function _hybrid_reverse_predict!(Abar, A, rec::CTSEMHybridTape{T},
    record::CTSEMPredictRecord{T}, n::Int) where {T}
    q = rec.q
    q == 0 && return nothing
    Xh = record.Xh
    Xbar = rec.Xbar
    size(rec.Xt) == size(Xbar) || (rec.Xt = zeros(T, size(Xbar)...))
    Xt = rec.Xt
    @inbounds for j in 1:n, i in 1:n
        acc = zero(T)
        for c in 1:q
            acc += Xbar[i, c] * Xh[j, c]
        end
        Abar[i, j] += acc
    end
    @inbounds for c in 1:q, j in 1:n
        acc = zero(T)
        for i in 1:n
            acc += A[i, j] * Xbar[i, c]
        end
        Xt[j, c] = acc
    end
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

"""Add the own-sd cotangents to the raw gradient, after the parameter layer."""
function _hybrid_scatter!(gradient::AbstractVector, rec::CTSEMHybridTape,
    h::CTSEMHybrid)
    rec.active || return gradient
    @inbounds for k in eachindex(h.psi_index)
        k <= length(rec.psibar) || break
        gradient[h.psi_index[k]] += rec.psibar[k]
    end
    return gradient
end
