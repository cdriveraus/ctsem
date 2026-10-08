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
# (`_laplace_unit_loglik_gradient`'s `hlogdet`). With `full` set the pass
# returns the integral itself, profile less `logdet(M) / 2`.

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
    if size(hb.X, 1) < n || size(hb.X, 2) < q
        hb.X = zeros(T, n, nh)
        hb.Xt = zeros(T, n, nh)
    end
    if size(hb.V, 1) < m || size(hb.V, 2) < q
        hb.V = zeros(T, m, nh)
        hb.W = zeros(T, m, nh)
        hb.e = zeros(T, m)
    end
    if size(hb.S, 1) < q
        hb.S = zeros(T, nh, nh)
        hb.s = zeros(T, nh)
    end
    fill!(view(hb.X, 1:n, 1:q), zero(T))
    fill!(view(hb.S, 1:q, 1:q), zero(T))
    fill!(view(hb.s, 1:q), zero(T))
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

"""One prediction substep's share: `X = A X`."""
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
    @inbounds for c in 1:q, r in 1:n
        X[r, c] = Xt[r, c]
    end
    return nothing
end

"""
The closing term, `1/2 b' M^-1 b`, less `logdet(M) / 2` when `hb.full`, with
`M = I + G S G` and `b = G s`; `logdet(M) / 2` left in `hb.half_logdet` either
way, and the term recorded for the reverse when `trace` is a tape. `values` is
the raw parameter vector the filter was given, which holds the own variances.
"""
function _hybrid_marginal!(hb::CTSEMHybridBuffers{T}, h::CTSEMHybrid,
    values::AbstractVector, trace) where {T}
    q = hb.q
    if q == 0
        _record_hybrid!(trace, hb, T[], T[])
        return zero(T)
    end
    g = Vector{T}(undef, q)
    dg = Vector{T}(undef, q)
    @inbounds for c in 1:q
        k = hb.item[c]
        raw = values[h.psi_index[k]]
        g[c] = h.slope[k] * (log1p_exp(2 * raw - 1) * h.sd_scale[k] + 1e-10)
        # d log1p_exp(z)/dz is the logistic function.
        dg[c] = h.slope[k] * h.sd_scale[k] * 2 / (1 + exp(-(2 * raw - 1)))
    end
    M = Matrix{T}(undef, q, q)
    @inbounds for b in 1:q, a in 1:q
        M[a, b] = g[a] * hb.S[a, b] * g[b] + (a == b ? one(T) : zero(T))
    end
    F = _ctsem_cholesky(M, q)
    issuccess(F) || return T(NaN)
    bvec = Vector{T}(undef, q)
    @inbounds for c in 1:q
        bvec[c] = g[c] * hb.s[c]
    end
    z = similar(bvec)
    ldiv!(z, F, bvec)
    half_logdet = zero(T)
    @inbounds for c in 1:q
        half_logdet += log(F.U[c, c])
    end
    hb.half_logdet = half_logdet
    _record_hybrid!(trace, hb, g, dg)
    return (hb.full ? -half_logdet : zero(T)) + dot(bvec, z) / 2
end

# Tape records. Every hook has a fallback that does nothing, so the plain
# primal and ForwardDiff paths compile to what they did.
_hybrid_record_X!(record, hb, n::Int) = nothing
function _hybrid_record_X!(record::Union{CTSEMUpdateRecord{T},CTSEMPredictRecord{T}},
    hb::CTSEMHybridBuffers{T}, n::Int) where {T}
    q = hb.active ? hb.q : 0
    if size(record.Xh) != (n, q)
        record.Xh = zeros(T, n, q)
    end
    q == 0 && return nothing
    @inbounds for c in 1:q, r in 1:n
        record.Xh[r, c] = hb.X[r, c]
    end
    return nothing
end

_record_hybrid!(trace, hb, g, dg) = nothing
function _record_hybrid!(tape::CTSEMAdjointTape{T}, hb::CTSEMHybridBuffers{T},
    g::Vector{T}, dg::Vector{T}) where {T}
    rec = tape.hyb
    q = hb.q
    rec.active = true
    rec.full = hb.full
    rec.q = q
    resize!(rec.item, q)
    copyto!(rec.item, 1, hb.item, 1, q)
    rec.col = copy(hb.col)
    rec.S = Matrix(view(hb.S, 1:q, 1:q))
    rec.s = Vector(view(hb.s, 1:q))
    rec.g = copy(g)
    rec.dg = copy(dg)
    _tape_push!(tape, :hybrid, 0)
    return nothing
end

"""
Reverse of the closing term, run first in the walk: the cotangents of `S` and
`s`, each item's own-sd cotangent into `psibar` (indexed by hybrid item), and
the carried cotangent of `X` started at zero.
"""
function _reverse_hybrid_marginal!(rec::CTSEMHybridTape{T}, n::Int, nh::Int) where {T}
    q = rec.q
    length(rec.psibar) == nh || (rec.psibar = zeros(T, nh))
    fill!(rec.psibar, zero(T))
    rec.Sbar = zeros(T, q, q)
    rec.sbar = zeros(T, q)
    rec.Xbar = zeros(T, n, q)
    q == 0 && return nothing
    g, S, s = rec.g, rec.S, rec.s
    M = Matrix{T}(undef, q, q)
    @inbounds for b in 1:q, a in 1:q
        M[a, b] = g[a] * S[a, b] * g[b] + (a == b ? one(T) : zero(T))
    end
    F = _ctsem_cholesky(M, q)
    C = Matrix{T}(undef, q, q)
    col = zeros(T, q)
    unit = zeros(T, q)
    @inbounds for b in 1:q
        fill!(unit, zero(T))
        unit[b] = one(T)
        ldiv!(col, F, unit)
        for a in 1:q
            C[a, b] = col[a]
        end
    end
    beta = zeros(T, q)
    @inbounds for b in 1:q, a in 1:q
        beta[a] += C[a, b] * g[b] * s[b]
    end
    # dll/dM = -1/2 (w C + beta beta'), dll/db = beta, with `w` 1 when the
    # pass returned the logdet and 0 for the profile.
    w = rec.full ? one(T) : zero(T)
    @inbounds for b in 1:q, a in 1:q
        Mbar = -(w * C[a, b] + beta[a] * beta[b]) / 2
        rec.Sbar[a, b] = g[a] * Mbar * g[b]
    end
    @inbounds for a in 1:q
        rec.sbar[a] = g[a] * beta[a]
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
observed LAMBDA rows (`Lbar`), which the base reverse folds into its own; and
`Xbar` carried back to before the update. `Ui` is `U^-1`, `K` the gain, `a`
`S^-1 y`, all as the base reverse recomputed them.
"""
function _hybrid_reverse_update!(rec::CTSEMHybridTape{T}, record::CTSEMUpdateRecord{T},
    Ui, K, a, n::Int) where {T}
    q = rec.q
    o = record.observed
    m = length(o)
    Xh = record.Xh
    L = record.Lambda
    for (name, rows, cols) in ((:V, m, q), (:W, m, q), (:Z, m, q), (:Vbar, m, q),
            (:Kbar, n, m), (:Rbar, m, m), (:Lbar, m, n), (:Xt, n, q))
        A = getfield(rec, name)
        size(A) == (rows, cols) || setfield!(rec, name, zeros(T, rows, cols))
    end
    length(rec.ybar) == m || (rec.ybar = zeros(T, m))
    V, W, Z, Vbar = rec.V, rec.W, rec.Z, rec.Vbar
    Xbar, Sbar, sbar = rec.Xbar, rec.Sbar, rec.sbar
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
        rec.Kbar[r, i] = acc
    end
    @inbounds for c in 1:q, i in 1:m
        acc = zero(T)
        for r in 1:n
            acc -= K[r, i] * Xbar[r, c]
        end
        # S += V' S^-1 V:  Vbar += 2 Z Sbar.  s += V' S^-1 y:  Vbar += a sbar'.
        for b in 1:q
            acc += 2 * Z[i, b] * Sbar[b, c]
        end
        Vbar[i, c] = acc + a[i] * sbar[c]
    end
    # Rbar (the innovation covariance) = -Z Sbar Z' - sym(Z sbar a');
    # ybar = Z sbar.
    @inbounds for i in 1:m
        zs = zero(T)
        for c in 1:q
            zs += Z[i, c] * sbar[c]
        end
        rec.ybar[i] = zs
    end
    @inbounds for j in 1:m, i in 1:m
        acc = zero(T)
        for c in 1:q
            zc = zero(T)
            for b in 1:q
                zc += Sbar[c, b] * Z[j, b]
            end
            acc -= Z[i, c] * zc
        end
        rec.Rbar[i, j] = acc - (rec.ybar[i] * a[j] + a[i] * rec.ybar[j]) / 2
    end
    # V = L X + E:  Lbar = Vbar X',  Xbar(before) = Xbar(after) + L' Vbar.
    @inbounds for k in 1:n, i in 1:m
        acc = zero(T)
        for c in 1:q
            acc += Vbar[i, c] * Xh[k, c]
        end
        rec.Lbar[i, k] = acc
    end
    @inbounds for c in 1:q, k in 1:n
        acc = Xbar[k, c]
        for i in 1:m
            acc += L[i, k] * Vbar[i, c]
        end
        rec.Xt[k, c] = acc
    end
    @inbounds for c in 1:q, k in 1:n
        Xbar[k, c] = rec.Xt[k, c]
    end
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
    size(rec.Xt) == (n, q) || (rec.Xt = zeros(T, n, q))
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
        rec.Xt[j, c] = acc
    end
    @inbounds for c in 1:q, j in 1:n
        Xbar[j, c] = rec.Xt[j, c]
    end
    return nothing
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
