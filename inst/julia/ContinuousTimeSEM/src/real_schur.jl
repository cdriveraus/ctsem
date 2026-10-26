################################################################################
# Real Schur decomposition and the Lyapunov solve on its factor, by hand
################################################################################
#
# The Lyapunov equation `A X + X A' + Q = 0` is O(k^3) through the real Schur
# form of `A` (Bartels-Stewart): `A = U S U'` with `S` quasi upper triangular,
# then `S Y + Y S' = -U' Q U` by back substitution over its blocks, and
# `X = U Y U'`. LAPACK's `schur!` and `trsyl!` take OpenBLAS's process-wide lock
# on every call and allocate, so the subject loop once used a packed solve over
# the k(k+1)/2 symmetric unknowns instead, which is O(k^6): at k = 20 a 210x210
# LU per factorisation, refactored per subject wherever the drift's Jacobian
# reads the state. These are the Schur route written out as the rest of
# `small_linalg.jl` is: no library call, no allocation, nothing touched outside
# the arguments, so it threads. Generic in the element type, so BigFloat checks
# take it too (`eps(T)`, `floatmin(T)`).
#
# The decomposition follows LAPACK's unblocked path, which is what LAPACK
# itself runs at these sizes: Householder reduction to Hessenberg form (as
# `dgehd2`, accumulating the reflectors), the Francis double-shift QR iteration
# of `dlahqr` (its Ahues-Kressner deflation test, exceptional shifts every ten
# iterations without a deflation, and its iteration limit), and `dlanv2`'s
# standard form for each 2x2 block: equal diagonal elements and off-diagonal
# elements of opposite sign, so a 2x2 block always holds a complex pair and
# its diagonal is the pair's real part. The scaling LAPACK applies against
# overflow and underflow is left out: the matrices here are drift Jacobians.

"""
    _ctsem_hessenberg!(H, Z, w, n)

Reduce the leading `n x n` block of `H` to upper Hessenberg form in place by
Householder reflections, and set `Z` to their product, so that the input is
`Z H Z'`. `w` is scratch of length at least `2n`.
"""
function _ctsem_hessenberg!(H::AbstractMatrix{T}, Z::AbstractMatrix{T},
    w::AbstractVector{T}, n::Int) where {T}
    @inbounds for j in 1:n, i in 1:n
        Z[i, j] = i == j ? one(T) : zero(T)
    end
    @inbounds for k in 1:(n - 2)
        alpha = H[k + 1, k]
        xnorm2 = zero(T)
        for i in (k + 2):n
            xnorm2 += H[i, k]^2
        end
        # Nothing below the subdiagonal: the reflector is the identity.
        iszero(xnorm2) && continue
        beta = -copysign(sqrt(alpha^2 + xnorm2), alpha)
        tau = (beta - alpha) / beta
        scale = inv(alpha - beta)
        # `P = I - tau v v'` with `v = [1; x / (alpha - beta)]` maps the column
        # below the diagonal to `beta e1` (`dlarfg`'s convention).
        w[k + 1] = one(T)
        for i in (k + 2):n
            w[i] = H[i, k] * scale
        end
        # From the left, columns k+1:n (column k is set below).
        for j in (k + 1):n
            s = zero(T)
            for i in (k + 1):n
                s += w[i] * H[i, j]
            end
            s *= tau
            for i in (k + 1):n
                H[i, j] -= s * w[i]
            end
        end
        # From the right, on H and on Z: row sums first, so both passes run
        # down columns.
        for M in (H, Z)
            for i in 1:n
                w[n + i] = zero(T)
            end
            for j in (k + 1):n
                wj = w[j]
                for i in 1:n
                    w[n + i] += M[i, j] * wj
                end
            end
            for j in (k + 1):n
                wj = tau * w[j]
                for i in 1:n
                    M[i, j] -= w[n + i] * wj
                end
            end
        end
        H[k + 1, k] = beta
        for i in (k + 2):n
            H[i, k] = zero(T)
        end
    end
    return H
end

"""
    _ctsem_lanv2(a, b, c, d)

The Schur factorisation of the real 2x2 matrix `[a b; c d]` in standard form,
LAPACK's `dlanv2`: returns `(aa, bb, cc, dd, cs, sn)` with
`[a b; c d] = R [aa bb; cc dd] R'`, `R = [cs -sn; sn cs]`. Real eigenvalues
give `cc = 0`; a complex pair gives `aa = dd` and `bb cc < 0`.
"""
function _ctsem_lanv2(a::T, b::T, c::T, d::T) where {T}
    epsT = eps(T)
    if iszero(c)
        cs = one(T); sn = zero(T)
    elseif iszero(b)
        # Swap rows and columns.
        cs = zero(T); sn = one(T)
        a, d = d, a
        b = -c; c = zero(T)
    elseif iszero(a - d) && sign(b) != sign(c)
        cs = one(T); sn = zero(T)
    else
        temp = a - d
        p = temp / 2
        bcmax = max(abs(b), abs(c))
        bcmis = min(abs(b), abs(c)) * sign(b) * sign(c)
        scale = max(abs(p), bcmax)
        z = (p / scale) * p + (bcmax / scale) * bcmis
        if z >= 4 * epsT
            # Real eigenvalues: compute a and d, then b and the rotation.
            z = p + copysign(sqrt(scale) * sqrt(z), p)
            a = d + z
            d = d - (bcmax / z) * bcmis
            tau = hypot(c, z)
            cs = z / tau
            sn = c / tau
            b = b - c
            c = zero(T)
        else
            # A complex pair, or real eigenvalues nearly equal: make the
            # diagonal elements equal.
            sigma = b + c
            tau = hypot(sigma, temp)
            cs = sqrt((one(T) + abs(sigma) / tau) / 2)
            sn = -(p / (tau * cs)) * (sigma >= 0 ? one(T) : -one(T))
            aa = a * cs + b * sn
            bb = -a * sn + b * cs
            cc = c * cs + d * sn
            dd = -c * sn + d * cs
            a = aa * cs + cc * sn
            b = bb * cs + dd * sn
            c = -aa * sn + cc * cs
            d = -bb * sn + dd * cs
            temp = (a + d) / 2
            a = temp
            d = temp
            if !iszero(c)
                if !iszero(b)
                    if sign(b) == sign(c)
                        # Real eigenvalues after all: reduce to triangular.
                        sab = sqrt(abs(b))
                        sac = sqrt(abs(c))
                        p = copysign(sab * sac, c)
                        tau = inv(sqrt(abs(b + c)))
                        a = temp + p
                        d = temp - p
                        b = b - c
                        c = zero(T)
                        cs1 = sab * tau
                        sn1 = sac * tau
                        temp = cs * cs1 - sn * sn1
                        sn = cs * sn1 + sn * cs1
                        cs = temp
                    end
                else
                    b = -c
                    c = zero(T)
                    temp = cs
                    cs = -sn
                    sn = temp
                end
            end
        end
    end
    return a, b, c, d, cs, sn
end

"""
    _ctsem_hqr!(H, Z, n) -> Bool

The Francis double-shift QR iteration of LAPACK's `dlahqr` on the upper
Hessenberg leading `n x n` block of `H`, with the full Schur form and the
transformations accumulated into `Z`: on success `H` is quasi upper triangular
with standardised 2x2 blocks (`_ctsem_lanv2`) and `Z` has been multiplied on
the right by the orthogonal factor. Returns false if some eigenvalue has not
converged within `30 max(10, n)` iterations, as `dlahqr` reports `info > 0`.
"""
function _ctsem_hqr!(H::AbstractMatrix{T}, Z::AbstractMatrix{T}, n::Int) where {T}
    n <= 1 && return true
    @inbounds for j in 1:(n - 3)
        H[j + 2, j] = zero(T)
        H[j + 3, j] = zero(T)
    end
    n >= 3 && (@inbounds H[n, n - 2] = zero(T))
    ulp = eps(T)
    smlnum = floatmin(T) * (T(n) / ulp)
    itmax = 30 * max(10, n)
    i = n
    @inbounds while i >= 1
        l = 1
        converged = false
        for its in 0:itmax
            # A single negligible subdiagonal element splits the problem.
            k = i
            while k > l
                hk = abs(H[k, k - 1])
                hk <= smlnum && break
                tst = abs(H[k - 1, k - 1]) + abs(H[k, k])
                if iszero(tst)
                    k - 2 >= 1 && (tst += abs(H[k - 1, k - 2]))
                    k + 1 <= n && (tst += abs(H[k + 1, k]))
                end
                if hk <= ulp * tst
                    # Ahues and Kressner's conservative deflation criterion.
                    hu = abs(H[k - 1, k])
                    ab = max(hk, hu)
                    ba = min(hk, hu)
                    dd = abs(H[k - 1, k - 1] - H[k, k])
                    aa = max(abs(H[k, k]), dd)
                    bb = min(abs(H[k, k]), dd)
                    s = aa + ab
                    ba * (ab / s) <= max(smlnum, ulp * (bb * (aa / s))) && break
                end
                k -= 1
            end
            l = k
            l > 1 && (H[l, l - 1] = zero(T))
            # A 1x1 or 2x2 block has split off.
            if l >= i - 1
                converged = true
                break
            end
            kdefl = its + 1
            if kdefl % 20 == 0
                s = abs(H[i, i - 1]) + abs(H[i - 1, i - 2])
                h11 = T(0.75) * s + H[i, i]
                h12 = T(-0.4375) * s
                h21 = s
                h22 = h11
            elseif kdefl % 10 == 0
                s = abs(H[l + 1, l]) + abs(H[l + 2, l + 1])
                h11 = T(0.75) * s + H[l, l]
                h12 = T(-0.4375) * s
                h21 = s
                h22 = h11
            else
                h11 = H[i - 1, i - 1]
                h21 = H[i, i - 1]
                h12 = H[i - 1, i]
                h22 = H[i, i]
            end
            s = abs(h11) + abs(h12) + abs(h21) + abs(h22)
            if iszero(s)
                rt1r = zero(T); rt1i = zero(T); rt2r = zero(T); rt2i = zero(T)
            else
                h11 /= s; h21 /= s; h12 /= s; h22 /= s
                tr = (h11 + h22) / 2
                det = (h11 - tr) * (h22 - tr) - h12 * h21
                rtdisc = sqrt(abs(det))
                if det >= 0
                    # A complex conjugate pair of shifts.
                    rt1r = tr * s
                    rt2r = rt1r
                    rt1i = rtdisc * s
                    rt2i = -rt1i
                else
                    # Real shifts: both the one nearer h22.
                    rt1r = tr + rtdisc
                    rt2r = tr - rtdisc
                    if abs(rt1r - h22) <= abs(rt2r - h22)
                        rt1r *= s
                        rt2r = rt1r
                    else
                        rt2r *= s
                        rt1r = rt2r
                    end
                    rt1i = zero(T); rt2i = zero(T)
                end
            end
            # Two consecutive small subdiagonal elements: start the bulge there.
            m = i - 2
            v1 = zero(T); v2 = zero(T); v3 = zero(T)
            while true
                h21s = H[m + 1, m]
                s = abs(H[m, m] - rt2r) + abs(rt2i) + abs(h21s)
                h21s = H[m + 1, m] / s
                v1 = h21s * H[m, m + 1] + (H[m, m] - rt1r) * ((H[m, m] - rt2r) / s) -
                    rt1i * (rt2i / s)
                v2 = h21s * (H[m, m] + H[m + 1, m + 1] - rt1r - rt2r)
                v3 = h21s * H[m + 2, m + 1]
                s = abs(v1) + abs(v2) + abs(v3)
                v1 /= s; v2 /= s; v3 /= s
                m == l && break
                h00 = abs(H[m, m - 1]) * (abs(v2) + abs(v3))
                h01 = abs(v1) * (abs(H[m - 1, m - 1]) + abs(H[m, m]) + abs(H[m + 1, m + 1]))
                h00 <= ulp * h01 && break
                m -= 1
            end
            # The double-shift QR sweep, chasing the bulge from m to i.
            for k in m:(i - 1)
                nr = min(3, i - k + 1)
                if k > m
                    v1 = H[k, k - 1]
                    v2 = H[k + 1, k - 1]
                    nr == 3 && (v3 = H[k + 2, k - 1])
                end
                # `dlarfg` on (v1; v2[; v3]): reflector `I - t1 [1 v2 v3]'[1 v2 v3]`.
                xnorm = nr == 3 ? hypot(v2, v3) : abs(v2)
                if iszero(xnorm)
                    t1 = zero(T)
                    beta = v1
                else
                    beta = -copysign(hypot(v1, xnorm), v1)
                    t1 = (beta - v1) / beta
                    sc = inv(v1 - beta)
                    v2 *= sc
                    v3 *= sc
                end
                if k > m
                    H[k, k - 1] = beta
                    H[k + 1, k - 1] = zero(T)
                    k < i - 1 && (H[k + 2, k - 1] = zero(T))
                elseif m > l
                    # As `dlahqr`: not `-H[k, k-1]`, which goes wrong when v2
                    # and v3 underflow.
                    H[k, k - 1] *= (one(T) - t1)
                end
                t2 = t1 * v2
                if nr == 3
                    t3 = t1 * v3
                    for j in k:n
                        sm = H[k, j] + v2 * H[k + 1, j] + v3 * H[k + 2, j]
                        H[k, j] -= sm * t1
                        H[k + 1, j] -= sm * t2
                        H[k + 2, j] -= sm * t3
                    end
                    for j in 1:min(k + 3, i)
                        sm = H[j, k] + v2 * H[j, k + 1] + v3 * H[j, k + 2]
                        H[j, k] -= sm * t1
                        H[j, k + 1] -= sm * t2
                        H[j, k + 2] -= sm * t3
                    end
                    for j in 1:n
                        sm = Z[j, k] + v2 * Z[j, k + 1] + v3 * Z[j, k + 2]
                        Z[j, k] -= sm * t1
                        Z[j, k + 1] -= sm * t2
                        Z[j, k + 2] -= sm * t3
                    end
                else
                    for j in k:n
                        sm = H[k, j] + v2 * H[k + 1, j]
                        H[k, j] -= sm * t1
                        H[k + 1, j] -= sm * t2
                    end
                    for j in 1:i
                        sm = H[j, k] + v2 * H[j, k + 1]
                        H[j, k] -= sm * t1
                        H[j, k + 1] -= sm * t2
                    end
                    for j in 1:n
                        sm = Z[j, k] + v2 * Z[j, k + 1]
                        Z[j, k] -= sm * t1
                        Z[j, k + 1] -= sm * t2
                    end
                end
            end
        end
        converged || return false
        if l == i - 1
            # A 2x2 block: put it in standard form and carry the rotation
            # through the rest of H and into Z.
            a, b, c, d, cs, sn = _ctsem_lanv2(H[i - 1, i - 1], H[i - 1, i], H[i, i - 1], H[i, i])
            H[i - 1, i - 1] = a
            H[i - 1, i] = b
            H[i, i - 1] = c
            H[i, i] = d
            for j in (i + 1):n
                x = H[i - 1, j]
                y = H[i, j]
                H[i - 1, j] = cs * x + sn * y
                H[i, j] = cs * y - sn * x
            end
            for j in 1:(i - 2)
                x = H[j, i - 1]
                y = H[j, i]
                H[j, i - 1] = cs * x + sn * y
                H[j, i] = cs * y - sn * x
            end
            for j in 1:n
                x = Z[j, i - 1]
                y = Z[j, i]
                Z[j, i - 1] = cs * x + sn * y
                Z[j, i] = cs * y - sn * x
            end
        end
        i = l - 1
    end
    return true
end

"""
    _ctsem_real_schur!(S, U, w, n) -> Bool

The real Schur decomposition of the leading `n x n` block of `S`, in place:
on success `S` holds the quasi upper triangular factor (standardised 2x2
blocks) and `U` the orthogonal factor, so that the input was `U S U'`. `w` is
scratch of length at least `2n`. False where the QR iteration did not
converge, or the input was not finite.
"""
function _ctsem_real_schur!(S::AbstractMatrix{T}, U::AbstractMatrix{T},
    w::AbstractVector{T}, n::Int) where {T}
    @inbounds for j in 1:n, i in 1:n
        isfinite(S[i, j]) || return false
    end
    _ctsem_hessenberg!(S, U, w, n)
    return _ctsem_hqr!(S, U, n)
end

"""The first index of the block of quasi triangular `S` that ends at `e`."""
@inline _ctsem_block_top(S, e::Int) = (e > 1 && !iszero(@inbounds S[e, e - 1])) ? e - 1 : e

"""The last index of the block of quasi triangular `S` that starts at `s`."""
@inline _ctsem_block_end(S, s::Int, n::Int) =
    (s < n && !iszero(@inbounds S[s + 1, s])) ? s + 1 : s

"""
    _ctsem_small_sylvester!(C, S, is, ie, js, je, R, M, transposed)

Solve one block equation of the quasi triangular Lyapunov solve and write it
into `C[is:ie, js:je]`: `S1 Y + Y S2' = R` (or `S1' Y + Y S2 = R` when
`transposed`), with `S1 = S[is:ie, is:ie]` and `S2 = S[js:je, js:je]`, each
1x1 or 2x2. `R` holds the right-hand side column-major on entry; `M` is 4x4
scratch. Its Kronecker form, at most 4x4, is solved by Gaussian elimination
with partial pivoting. A singular block -- an eigenvalue sum of zero -- gives
non-finite values rather than an error, as the packed solve's zero pivot does.
"""
function _ctsem_small_sylvester!(C::AbstractMatrix{T}, S::AbstractMatrix{T},
    is::Int, ie::Int, js::Int, je::Int, R::AbstractVector{T}, M::AbstractMatrix{T},
    transposed::Bool) where {T}
    p = ie - is + 1
    q = je - js + 1
    m = p * q
    @inbounds begin
        if m == 1
            C[is, js] = R[1] / (S[is, is] + S[js, js])
            return C
        end
        # vec(S1 Y) = (I ⊗ S1) vec(Y) and vec(Y S2') = (S2 ⊗ I) vec(Y);
        # transposed, S1' and S2' in their places.
        for c in 1:m, r in 1:m
            M[r, c] = zero(T)
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
                R[c], R[piv] = R[piv], R[c]
            end
            for r in (c + 1):m
                f = M[r, c] / M[c, c]
                for cc in (c + 1):m
                    M[r, cc] -= f * M[c, cc]
                end
                R[r] -= f * R[c]
            end
        end
        for r in m:-1:1
            acc = R[r]
            for cc in (r + 1):m
                acc -= M[r, cc] * R[cc]
            end
            R[r] = acc / M[r, r]
        end
        for b in 1:q, a in 1:p
            C[is + a - 1, js + b - 1] = R[a + (b - 1) * p]
        end
    end
    return C
end

"""
    _ctsem_quasitri_lyap!(C, S, n, R, M, transposed)

Bartels-Stewart on the Schur factor: overwrite the symmetric right-hand side
`C` (leading `n x n` block) with the symmetric `Y` solving `S Y + Y S' = C`,
or `S' Y + Y S = C` when `transposed`, for quasi upper triangular `S` with
standardised 2x2 blocks. Only one triangle's blocks are solved; the other is
copied across, so `Y` comes back symmetric. `R` (length 4) and `M` (4x4) are
scratch. The transposed form serves the pullback from the same factors: if
`A = U S U'` then `A' = U S' U'`.
"""
function _ctsem_quasitri_lyap!(C::AbstractMatrix{T}, S::AbstractMatrix{T}, n::Int,
    R::AbstractVector{T}, M::AbstractMatrix{T}, transposed::Bool) where {T}
    if !transposed
        # Block (i, j) reads the blocks below it in its column and to its right
        # in its row: columns from the last, rows upwards from the diagonal.
        je = n
        @inbounds while je >= 1
            js = _ctsem_block_top(S, je)
            for jj in js:je, ii in (je + 1):n
                C[ii, jj] = C[jj, ii]
            end
            ie = je
            while ie >= 1
                is = _ctsem_block_top(S, ie)
                for jj in js:je, ii in is:ie
                    r = C[ii, jj]
                    for l in (ie + 1):n
                        r -= S[ii, l] * C[l, jj]
                    end
                    for l in (je + 1):n
                        r -= C[ii, l] * S[jj, l]
                    end
                    R[(ii - is + 1) + (jj - js) * (ie - is + 1)] = r
                end
                _ctsem_small_sylvester!(C, S, is, ie, js, je, R, M, false)
                ie = is - 1
            end
            je = js - 1
        end
    else
        # `S'` is lower quasi triangular: columns from the first, rows
        # downwards from the diagonal.
        js = 1
        @inbounds while js <= n
            je = _ctsem_block_end(S, js, n)
            for jj in js:je, ii in 1:(js - 1)
                C[ii, jj] = C[jj, ii]
            end
            is = js
            while is <= n
                ie = _ctsem_block_end(S, is, n)
                for jj in js:je, ii in is:ie
                    r = C[ii, jj]
                    for l in 1:(is - 1)
                        r -= S[l, ii] * C[l, jj]
                    end
                    for l in 1:(js - 1)
                        r -= C[ii, l] * S[l, jj]
                    end
                    R[(ii - is + 1) + (jj - js) * (ie - is + 1)] = r
                end
                _ctsem_small_sylvester!(C, S, is, ie, js, je, R, M, true)
                is = ie + 1
            end
            js = je + 1
        end
    end
    return C
end
