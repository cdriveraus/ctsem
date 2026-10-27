# Tangent-aware dense primitives in Float64 for the hand-coded second order
# (`tangent_kernels.jl`, `tangent_sweep.jl`, `tangent_hessian.jl`): each
# quantity carries its value and W tangents, a matrix's as a `W x r x c` array,
# the directions innermost and W a type parameter, so that every tangent
# operation is a W-wide fused multiply-add over contiguous memory -- what a
# dual does per scalar, without the dual type. A product's value and tangents
# are computed in one pass over (j, i, q), the tangents summed in a register
# tuple and stored once; transposed operands are only another index order.
#
# Each quantity also says whether its tangents may be nonzero (`nz`); where it
# says not, they are zero, and a product skips the terms they would feed. A
# Hessian sweep's directions are one chunk of parameters, and most quantities
# do not depend on most chunks -- the measurement's on none but its own -- which
# a dual carries as zeros through every multiply. Nothing here allocates once
# built.
#
# Three ways of writing this lost to the dual, each found by timing one
# primitive against the same loop on duals: directions outermost (four times
# slower), a tuple stored by a loop indexing it by `w` (spilled to the stack,
# its vector reloads stalling: four to five times a primitive), and a column
# addressed `d[w, i, j]` (the vectoriser left its first and last directions
# scalar). Transposition, a triangular factor and a symmetric result are
# run-time arguments of one product kernel: as type parameters they compiled a
# 12-wide kernel per product signature, most of this code's image.

struct TanMat{W}
    v::Matrix{Float64}
    d::Array{Float64,3}     # W x r x c
    nz::Base.RefValue{Bool}
end
TanMat{W}(r::Int, c::Int) where {W} = TanMat{W}(zeros(r, c), zeros(W, r, c), Ref(true))

struct TanVec{W}
    v::Vector{Float64}
    d::Matrix{Float64}      # W x len
    nz::Base.RefValue{Bool}
end
TanVec{W}(r::Int) where {W} = TanVec{W}(zeros(r), zeros(W, r), Ref(true))

@inline _tan_nz(A) = A.nz[]
"""Mark whether `A`'s tangents may be nonzero, zeroing them where not."""
@inline function _tan_setnz!(A, on::Bool)
    on || (A.nz[] && fill!(A.d, 0.0))
    A.nz[] = on
    return A
end
"""Mark `A` from its tangents themselves (an input)."""
_tan_setnz_from!(A) = (A.nz[] = any(!iszero, A.d); A)

# A plain loop, not `copyto!`, whose call costs more than these copies.
function _tan_copyto!(dst::Array{Float64}, src::Array{Float64})
    @inbounds @simd for k in eachindex(dst, src)
        dst[k] = src[k]
    end
    return nothing
end
function _tan_copy!(C::Union{TanMat,TanVec}, A::Union{TanMat,TanVec})
    _tan_copyto!(C.v, A.v)
    if _tan_nz(A)
        _tan_copyto!(C.d, A.d)
        C.nz[] = true
    else
        _tan_setnz!(C, false)
    end
    return C
end

# A tangent column as a register tuple, and its store. Elementwise work reads
# its operands' columns into tuples and stores the result once: a `@simd` loop
# over W storing into one array while reading others is alias-checked at every
# entry, which at W = 12 cost three times the arithmetic.
#
# Column (i, j) is `d[o + 1], ..., d[o + W]` from its offset `o`, and every
# access goes through the offset: indexed `d[w, i, j]`, the first direction's
# address was computed apart from the rest, and the vectoriser left the first
# and last directions scalar.
@inline _tan_off(d::Array{Float64,3}, ::Val{W}, i, j) where {W} =
    W * ((i - 1) + size(d, 2) * (j - 1))
@inline _tan_off(d::Matrix{Float64}, ::Val{W}, i) where {W} = W * (i - 1)
@inline _tan_ld(d, o, ::Val{W}) where {W} = ntuple(w -> @inbounds(d[o + w]), Val(W))
# Stores are unrolled by `ntuple`, so each reads the tuple at a constant
# index: a loop indexing it by `w` puts the tuple on the stack, and its vector
# reloads of the scalar stores stall (4-5 times the dual's time a primitive).
@inline function _tan_st!(d, o, t::NTuple{W,Float64}) where {W}
    ntuple(w -> (@inbounds(d[o + w] = t[w]); nothing), Val(W))
    return nothing
end
@inline _tan_col(A::TanMat{W}, i, j) where {W} = _tan_ld(A.d, _tan_off(A.d, Val(W), i, j), Val(W))
@inline _tan_col(x::TanVec{W}, i) where {W} = _tan_ld(x.d, _tan_off(x.d, Val(W), i), Val(W))
@inline _tan_put!(A::TanMat{W}, i, j, t::NTuple{W,Float64}) where {W} =
    _tan_st!(A.d, _tan_off(A.d, Val(W), i, j), t)
@inline _tan_put!(x::TanVec{W}, i, t::NTuple{W,Float64}) where {W} =
    _tan_st!(x.d, _tan_off(x.d, Val(W), i), t)
@inline _tan_zeros(::Val{W}) where {W} = ntuple(_ -> 0.0, Val(W))
@inline _tan_times(a, x::NTuple{W,Float64}) where {W} = ntuple(w -> a * x[w], Val(W))
@inline _tan_axpby(a, x::NTuple{W,Float64}, b, y::NTuple{W,Float64}) where {W} =
    ntuple(w -> muladd(a, x[w], b * y[w]), Val(W))
# acc + s d[:, i, j] (or the column at offset `o`), as a function so that a
# loop's reassigned accumulator is never captured by a closure (which would
# box it).
@inline _tan_axpyo(acc::NTuple{W,Float64}, s, d, o) where {W} =
    ntuple(w -> @inbounds(muladd(s, d[o + w], acc[w])), Val(W))
@inline _tan_axpy(acc::NTuple{W,Float64}, s, d, i, j) where {W} =
    _tan_axpyo(acc, s, d, _tan_off(d, Val(W), i, j))
@inline _tan_axpyv(acc::NTuple{W,Float64}, s, d, i) where {W} =
    _tan_axpyo(acc, s, d, _tan_off(d, Val(W), i))

# Store a register accumulator into the tangent column `d[:, i, j]` (or
# `d[:, i]`) as `alpha acc + beta d`.
@inline function _tan_store!(d, o, acc::NTuple{W,Float64}, alpha, beta) where {W}
    t = iszero(beta) ? _tan_times(alpha, acc) : _tan_axpby(alpha, acc, beta, _tan_ld(d, o, Val(W)))
    return _tan_st!(d, o, t)
end

# `acc + dA[:, ia, ja] b + a dB[:, ib, jb]` as a register tuple, either term
# left out where its tangents are zero, from the columns' offsets `oa`, `ob`. A
# function of its own: a closure over the loop's accumulator, which the loop
# reassigns, would box it.
@inline _tan_fma(acc::NTuple{W,Float64}, ::Val{true}, ::Val{true}, dA, oa, b, a, dB,
    ob) where {W} =
    ntuple(w -> @inbounds(muladd(dA[oa + w], b, muladd(a, dB[ob + w], acc[w]))), Val(W))
@inline _tan_fma(acc::NTuple{W,Float64}, ::Val{true}, ::Val{false}, dA, oa, b, a, dB,
    ob) where {W} =
    ntuple(w -> @inbounds(muladd(dA[oa + w], b, acc[w])), Val(W))
@inline _tan_fma(acc::NTuple{W,Float64}, ::Val{false}, ::Val{true}, dA, oa, b, a, dB,
    ob) where {W} =
    ntuple(w -> @inbounds(muladd(a, dB[ob + w], acc[w])), Val(W))

# A product's structure, as run-time arguments: whether a triangular factor
# (an upper triangle, its tangent too: `Ui` and `dUi`) bounds the inner index,
# and where (`tri`), and whether the result is symmetric (`sym`: its upper
# triangle computed, then mirrored). As type parameters, with the
# transpositions and the scalars' types, they compiled one 12-wide kernel per
# product signature, and those were most of the hand-coded kernels' image.
@inline _tan_tricode(::Val{:full}) = 0
@inline _tan_tricode(::Val{:Bupper}) = 1     # A U: U[q, j], q <= j
@inline _tan_tricode(::Val{:BupperT}) = 2    # A U': U[j, q], q >= j
@inline _tan_tricode(::Val{:Aupper}) = 3     # U B: U[i, q], q >= i
@inline _tan_tricode(::Val{:AupperT}) = 4    # U' B: U[q, i], q <= i
@inline _tan_symflag(::Val{S}) where {S} = S::Bool
@inline function _tan_qbounds(tri::Int, i::Int, j::Int, k::Int)
    tri == 0 && return (1, k)
    tri == 1 && return (1, min(j, k))
    tri == 2 && return (j, k)
    tri == 3 && return (i, k)
    return (1, min(i, k))
end
function _tan_mirror!(C::TanMat{W}) where {W}
    n = size(C.v, 1)
    @inbounds for j in 1:n, i in (j + 1):n
        C.v[i, j] = C.v[j, i]
        _tan_put!(C, i, j, _tan_col(C, j, i))
    end
    return nothing
end

"""
The product kernel: `C = alpha op(A) op(B) + beta C` with tangents. `ta`, `tb`
say which factor is transposed, read through strides; with `tri` and `sym`
they are run-time arguments, so one compiled kernel serves every product. Only
whether each factor's tangents enter (`FA`, `FB`) is a type parameter, since it
takes terms out of the inner loop.
"""
function _tan_prod!(C::TanMat{W}, A::TanMat{W}, B::TanMat{W}, alpha::Float64, beta::Float64,
    tri::Int, sym::Bool, ta::Bool, tb::Bool, fa::Val{FA}, fb::Val{FB}) where {W,FA,FB}
    lda, ldb, ldc = size(A.v, 1), size(B.v, 1), size(C.v, 1)
    r = ta ? size(A.v, 2) : lda
    k = ta ? lda : size(A.v, 2)
    c = tb ? ldb : size(B.v, 2)
    # op(A)[i, q] is A.v[1 + (i - 1) sai + (q - 1) saq], op(B)[q, j] likewise;
    # a tangent column sits at W times its entry's linear offset.
    sai, saq = ta ? (lda, 1) : (1, lda)
    sbq, sbj = tb ? (ldb, 1) : (1, ldb)
    Av, Bv, Cv, dA, dB, dC = A.v, B.v, C.v, A.d, B.d, C.d
    bz = iszero(beta)
    @inbounds for j in 1:c, i in 1:(sym ? j : r)
        q0, q1 = _tan_qbounds(tri, i, j, k)
        acc = 0.0
        dacc = _tan_zeros(Val(W))
        for q in q0:q1
            la = (i - 1) * sai + (q - 1) * saq
            lb = (q - 1) * sbq + (j - 1) * sbj
            a = Av[la + 1]
            b = Bv[lb + 1]
            acc = muladd(a, b, acc)
            if FA || FB
                dacc = _tan_fma(dacc, fa, fb, dA, W * la, b, a, dB, W * lb)
            end
        end
        lc = (i - 1) + (j - 1) * ldc
        Cv[lc + 1] = bz ? alpha * acc : muladd(alpha, acc, beta * Cv[lc + 1])
        (FA || FB) && _tan_store!(dC, W * lc, dacc, alpha, beta)
    end
    sym && _tan_mirror!(C)
    return C
end

function _tan_prod_flags!(C::TanMat{W}, A::TanMat{W}, B::TanMat{W}, alpha::Float64, beta::Float64,
    tri::Int, sym::Bool, ta::Bool, tb::Bool) where {W}
    fa, fb = _tan_nz(A), _tan_nz(B)
    if fa && fb
        _tan_prod!(C, A, B, alpha, beta, tri, sym, ta, tb, Val(true), Val(true))
    elseif fa
        _tan_prod!(C, A, B, alpha, beta, tri, sym, ta, tb, Val(true), Val(false))
    elseif fb
        _tan_prod!(C, A, B, alpha, beta, tri, sym, ta, tb, Val(false), Val(true))
    else
        _tan_prod!(C, A, B, alpha, beta, tri, sym, ta, tb, Val(false), Val(false))
        # No tangent entered: what is left is beta times C's own.
        if !iszero(beta) && _tan_nz(C)
            isone(beta) || (C.d .*= beta)
        else
            _tan_setnz!(C, false)
        end
        return C
    end
    # The stores above read `beta d` where `beta` is nonzero: zeros where C's
    # tangents were zero, so they are right either way.
    (iszero(beta) || !_tan_nz(C)) && (C.nz[] = true)
    return C
end

"""C = alpha A B + beta C, with tangents."""
@inline _tan_mul!(C::TanMat, A::TanMat, B::TanMat, alpha=true, beta=false, tri=Val(:full), sym=Val(false)) =
    _tan_prod_flags!(C, A, B, Float64(alpha), Float64(beta), _tan_tricode(tri), _tan_symflag(sym),
        false, false)
"""C = alpha A B' + beta C, with tangents."""
@inline _tan_mulNT!(C::TanMat, A::TanMat, B::TanMat, alpha=true, beta=false, tri=Val(:full), sym=Val(false)) =
    _tan_prod_flags!(C, A, B, Float64(alpha), Float64(beta), _tan_tricode(tri), _tan_symflag(sym),
        false, true)
"""C = alpha A' B + beta C, with tangents."""
@inline _tan_mulTN!(C::TanMat, A::TanMat, B::TanMat, alpha=true, beta=false, tri=Val(:full), sym=Val(false)) =
    _tan_prod_flags!(C, A, B, Float64(alpha), Float64(beta), _tan_tricode(tri), _tan_symflag(sym),
        true, false)

"""The matrix-vector kernel: `y = alpha op(A) x + beta y`, with tangents."""
function _tan_prodvec!(y::TanVec{W}, A::TanMat{W}, x::TanVec{W}, alpha::Float64, beta::Float64,
    tri::Int, ta::Bool, fa::Val{FA}, fx::Val{FX}) where {W,FA,FX}
    lda = size(A.v, 1)
    r = ta ? size(A.v, 2) : lda
    k = ta ? lda : size(A.v, 2)
    sai, saq = ta ? (lda, 1) : (1, lda)
    Av, xv, yv, dA, dx, dy = A.v, x.v, y.v, A.d, x.d, y.d
    bz = iszero(beta)
    @inbounds for i in 1:r
        q0, q1 = _tan_qbounds(tri, i, 1, k)
        acc = 0.0
        dacc = _tan_zeros(Val(W))
        for q in q0:q1
            la = (i - 1) * sai + (q - 1) * saq
            a = Av[la + 1]
            b = xv[q]
            acc = muladd(a, b, acc)
            if FA || FX
                dacc = _tan_fma(dacc, fa, fx, dA, W * la, b, a, dx, W * (q - 1))
            end
        end
        yv[i] = bz ? alpha * acc : muladd(alpha, acc, beta * yv[i])
        (FA || FX) && _tan_store!(dy, W * (i - 1), dacc, alpha, beta)
    end
    return y
end

function _tan_prodvec_flags!(y::TanVec{W}, A::TanMat{W}, x::TanVec{W}, alpha::Float64, beta::Float64,
    tri::Int, ta::Bool) where {W}
    fa, fx = _tan_nz(A), _tan_nz(x)
    if fa && fx
        _tan_prodvec!(y, A, x, alpha, beta, tri, ta, Val(true), Val(true))
    elseif fa
        _tan_prodvec!(y, A, x, alpha, beta, tri, ta, Val(true), Val(false))
    elseif fx
        _tan_prodvec!(y, A, x, alpha, beta, tri, ta, Val(false), Val(true))
    else
        _tan_prodvec!(y, A, x, alpha, beta, tri, ta, Val(false), Val(false))
        if !iszero(beta) && _tan_nz(y)
            isone(beta) || (y.d .*= beta)
        else
            _tan_setnz!(y, false)
        end
        return y
    end
    (iszero(beta) || !_tan_nz(y)) && (y.nz[] = true)
    return y
end

"""y = alpha A x + beta y, with tangents."""
@inline _tan_mulvec!(y::TanVec, A::TanMat, x::TanVec, alpha=true, beta=false, tri=Val(:full)) =
    _tan_prodvec_flags!(y, A, x, Float64(alpha), Float64(beta), _tan_tricode(tri), false)
"""y = alpha A' x + beta y, with tangents."""
@inline _tan_mulTvec!(y::TanVec, A::TanMat, x::TanVec, alpha=true, beta=false, tri=Val(:full)) =
    _tan_prodvec_flags!(y, A, x, Float64(alpha), Float64(beta), _tan_tricode(tri), true)

"""S = (A + A') / 2, with tangents."""
function _tan_sym!(S::TanMat{W}, A::TanMat{W}) where {W}
    n = size(A.v, 1)
    Av, Sv = A.v, S.v
    @inbounds for j in 1:n, i in 1:j
        s = (Av[i, j] + Av[j, i]) / 2
        Sv[i, j] = s
        Sv[j, i] = s
    end
    _tan_nz(A) || return _tan_setnz!(S, false)
    @inbounds for j in 1:n, i in 1:j
        t = _tan_axpby(0.5, _tan_col(A, i, j), 0.5, _tan_col(A, j, i))
        _tan_put!(S, i, j, t)
        _tan_put!(S, j, i, t)
    end
    S.nz[] = true
    return S
end

"""C = alpha op(A) + beta C, with tangents; `tr` `Val(:N)` or `Val(:T)`."""
@inline _tan_add!(C::TanMat, alpha, A::TanMat, beta=true, tr=Val(:N)) =
    _tan_add_impl!(C, Float64(alpha), A, Float64(beta), tr === Val(:T))
function _tan_add_impl!(C::TanMat{W}, alpha::Float64, A::TanMat{W}, beta::Float64, tr::Bool) where {W}
    r, c = size(C.v)
    lda = size(A.v, 1)
    # op(A)[i, j] is A.v[1 + (i - 1) sai + (j - 1) saj].
    sai, saj = tr ? (lda, 1) : (1, lda)
    Cv, Av, dA, dC = C.v, A.v, A.d, C.d
    bz = iszero(beta)
    @inbounds for j in 1:c, i in 1:r
        la = (i - 1) * sai + (j - 1) * saj
        Cv[i, j] = bz ? alpha * Av[la + 1] : muladd(alpha, Av[la + 1], beta * Cv[i, j])
    end
    if _tan_nz(A)
        keep = !bz && _tan_nz(C)
        @inbounds for j in 1:c, i in 1:r
            a = _tan_ld(dA, W * ((i - 1) * sai + (j - 1) * saj), Val(W))
            o = W * ((i - 1) + (j - 1) * r)
            _tan_st!(dC, o, keep ? _tan_axpby(alpha, a, beta, _tan_ld(dC, o, Val(W))) : _tan_times(alpha, a))
        end
        C.nz[] = true
    elseif bz || !_tan_nz(C)
        _tan_setnz!(C, false)
    elseif !isone(beta)
        C.d .*= beta
    end
    return C
end

"""C = alpha x y' + beta C, with tangents."""
@inline _tan_ger!(C::TanMat, alpha, x::TanVec, y::TanVec, beta=true) =
    _tan_ger_flags!(C, Float64(alpha), x, y, Float64(beta))
function _tan_ger_flags!(C::TanMat{W}, alpha::Float64, x::TanVec{W}, y::TanVec{W}, beta::Float64) where {W}
    xv, yv, Cv = x.v, y.v, C.v
    bz = iszero(beta)
    @inbounds for j in eachindex(yv), i in eachindex(xv)
        Cv[i, j] = bz ? alpha * xv[i] * yv[j] : muladd(alpha * xv[i], yv[j], beta * Cv[i, j])
    end
    fx, fy = _tan_nz(x), _tan_nz(y)
    keep = !bz && _tan_nz(C)
    if fx && fy
        _tan_ger_impl!(C, alpha, x, y, beta, keep, Val(true), Val(true))
    elseif fx
        _tan_ger_impl!(C, alpha, x, y, beta, keep, Val(true), Val(false))
    elseif fy
        _tan_ger_impl!(C, alpha, x, y, beta, keep, Val(false), Val(true))
    else
        keep || return _tan_setnz!(C, false)
        isone(beta) || (C.d .*= beta)
        return C
    end
    C.nz[] = true
    return C
end
function _tan_ger_impl!(C::TanMat{W}, alpha::Float64, x::TanVec{W}, y::TanVec{W}, beta::Float64, keep::Bool, ::Val{FX},
    ::Val{FY}) where {W,FX,FY}
    xv, yv = x.v, y.v
    @inbounds for j in eachindex(yv), i in eachindex(xv)
        t = keep ? _tan_times(beta, _tan_col(C, i, j)) : _tan_zeros(Val(W))
        FX && (t = _tan_axpyv(t, alpha * yv[j], x.d, i))
        FY && (t = _tan_axpyv(t, alpha * xv[i], y.d, j))
        _tan_put!(C, i, j, t)
    end
    return nothing
end

"""A *= s, with tangents."""
function _tan_scale!(A::Union{TanMat,TanVec}, s)
    A.v .*= s
    _tan_nz(A) && (A.d .*= s)
    return A
end

"""
C = sum of `a A` over `terms`, a tuple of `(a, A)` pairs, plus `beta C`, with
tangents. The operands share C's shape, so one offset serves every column.
"""
function _tan_comb!(C::TanMat{W}, terms::Tuple, beta=false) where {W}
    r, c = size(C.v)
    Cv = C.v
    bz = iszero(beta)
    @inbounds for j in 1:c, i in 1:r
        Cv[i, j] = _tan_comb_val(bz ? 0.0 : beta * Cv[i, j], terms, i, j)
    end
    keep = !bz && _tan_nz(C)
    if _tan_comb_nz(terms)
        d = C.d
        @inbounds for j in 1:c, i in 1:r
            o = _tan_off(d, Val(W), i, j)
            t = keep ? _tan_times(beta, _tan_ld(d, o, Val(W))) : _tan_zeros(Val(W))
            _tan_st!(d, o, _tan_comb_tan(t, terms, o))
        end
        C.nz[] = true
    elseif !keep
        _tan_setnz!(C, false)
    elseif !isone(beta)
        C.d .*= beta
    end
    return C
end
@inline _tan_comb_val(v, ::Tuple{}, i, j) = v
@inline _tan_comb_val(v, terms::Tuple, i, j) =
    _tan_comb_val(muladd(terms[1][1], @inbounds(terms[1][2].v[i, j]), v), Base.tail(terms), i, j)
@inline _tan_comb_nz(::Tuple{}) = false
@inline _tan_comb_nz(terms::Tuple) = _tan_nz(terms[1][2]) || _tan_comb_nz(Base.tail(terms))
@inline _tan_comb_tan(t, ::Tuple{}, o) = t
@inline function _tan_comb_tan(t, terms::Tuple, o)
    a, A = terms[1]
    return _tan_comb_tan(_tan_nz(A) ? _tan_axpyo(t, a, A.d, o) : t, Base.tail(terms), o)
end

"""
X = Q \\ X in place, with tangents, given Q's value factors `F`, `piv` from
`_lu_factor_generic!`: `dX = Q^-1 (dX - dQ X)`, the right-hand side's
tangents solved as W-tuples by the same elimination as the values, on one
factorization of the values -- where a dual factors its own matrix, the
directions riding every operation of the elimination.
"""
function _tan_solve!(X::Union{TanMat{W},TanVec{W}}, F::Matrix{Float64}, piv::Vector{Int},
    Q::TanMat{W}) where {W}
    n = size(X.v, 1)
    r = size(X.v, 2)
    _lu_solve_generic!(F, piv, X.v, n)
    fq = _tan_nz(Q)
    (fq || _tan_nz(X)) || return X
    # Entry (i, c) of X, a matrix or a vector (c = 1), by its linear offset.
    d, Xv, dQ, V = X.d, X.v, Q.d, Val(W)
    off(i, c) = W * ((i - 1) + n * (c - 1))
    if fq
        @inbounds for c in 1:r, i in 1:n
            o = off(i, c)
            t = _tan_ld(d, o, V)
            for q in 1:n
                t = _tan_axpyo(t, -Xv[q + n * (c - 1)], dQ, W * ((i - 1) + n * (q - 1)))
            end
            _tan_st!(d, o, t)
        end
    end
    @inbounds for k in 1:n
        p = piv[k]
        if p != k
            for c in 1:r
                ok, op = off(k, c), off(p, c)
                tk = _tan_ld(d, ok, V)
                _tan_st!(d, ok, _tan_ld(d, op, V))
                _tan_st!(d, op, tk)
            end
        end
        for i in (k + 1):n
            lik = F[i, k]
            for c in 1:r
                oi = off(i, c)
                _tan_st!(d, oi, _tan_axpyo(_tan_ld(d, oi, V), -lik, d, off(k, c)))
            end
        end
    end
    @inbounds for i in n:-1:1
        ainv = inv(F[i, i])
        for c in 1:r
            oi = off(i, c)
            t = _tan_ld(d, oi, V)
            for col in (i + 1):n
                t = _tan_axpyo(t, -F[i, col], d, off(col, c))
            end
            _tan_st!(d, oi, _tan_times(ainv, t))
        end
    end
    X.nz[] = true
    return X
end

"""Copy C's lower triangle onto its upper, tangents too, as the filter does."""
function _tan_lower_to_upper!(C::TanMat{W}) where {W}
    n = size(C.v, 1)
    on = _tan_nz(C)
    @inbounds for j in 1:n, i in (j + 1):n
        C.v[j, i] = C.v[i, j]
        on && _tan_put!(C, j, i, _tan_col(C, i, j))
    end
    return C
end

"""C = A[rows, cols], with tangents."""
function _tan_gather!(C::TanMat{W}, A::TanMat{W}, rows::AbstractVector{Int},
    cols::AbstractVector{Int}) where {W}
    on = _tan_nz(A)
    @inbounds for (jj, j) in enumerate(cols), (ii, i) in enumerate(rows)
        C.v[ii, jj] = A.v[i, j]
        on && _tan_put!(C, ii, jj, _tan_col(A, i, j))
    end
    on ? (C.nz[] = true) : _tan_setnz!(C, false)
    return C
end
"""c = a[rows], with tangents."""
function _tan_gather!(c::TanVec{W}, a::TanVec{W}, rows::AbstractVector{Int}) where {W}
    on = _tan_nz(a)
    @inbounds for (ii, i) in enumerate(rows)
        c.v[ii] = a.v[i]
        on && _tan_put!(c, ii, _tan_col(a, i))
    end
    on ? (c.nz[] = true) : _tan_setnz!(c, false)
    return c
end
"""A[rows, cols] += alpha C, with tangents."""
function _tan_scatter!(A::TanMat{W}, alpha::Float64, C::TanMat{W}, rows::AbstractVector{Int},
    cols::AbstractVector{Int}) where {W}
    on = _tan_nz(C)
    @inbounds for (jj, j) in enumerate(cols), (ii, i) in enumerate(rows)
        A.v[i, j] = muladd(alpha, C.v[ii, jj], A.v[i, j])
        on && _tan_put!(A, i, j, _tan_axpby(1.0, _tan_col(A, i, j), alpha, _tan_col(C, ii, jj)))
    end
    on && (A.nz[] = true)
    return A
end
"""a[rows] += alpha c, with tangents."""
function _tan_scatter!(a::TanVec{W}, alpha::Float64, c::TanVec{W}, rows::AbstractVector{Int}) where {W}
    on = _tan_nz(c)
    @inbounds for (ii, i) in enumerate(rows)
        a.v[i] = muladd(alpha, c.v[ii], a.v[i])
        on && _tan_put!(a, i, _tan_axpby(1.0, _tan_col(a, i), alpha, _tan_col(c, ii)))
    end
    on && (a.nz[] = true)
    return a
end

"""C += a I (values: a constant has no tangent)."""
function _tan_diag!(C::TanMat, a)
    @inbounds for i in 1:min(size(C.v)...)
        C.v[i, i] += a
    end
    return C
end

"""z = a x + b y, with tangents (`z` may be `x` or `y`)."""
@inline _tan_lincomb!(z::TanVec, a, x::TanVec, b, y::TanVec) = _tan_lincomb_impl!(z, Float64(a), x, Float64(b), y)
function _tan_lincomb_impl!(z::TanVec{W}, a::Float64, x::TanVec{W}, b::Float64, y::TanVec{W}) where {W}
    zv, xv, yv = z.v, x.v, y.v
    on = _tan_nz(x) || _tan_nz(y)
    @inbounds for i in eachindex(zv)
        zv[i] = muladd(a, xv[i], b * yv[i])
        on && _tan_put!(z, i, _tan_axpby(a, _tan_col(x, i), b, _tan_col(y, i)))
    end
    on ? (z.nz[] = true) : _tan_setnz!(z, false)
    return z
end

"""
`Ui = U^-1` for `S = U'U` (upper, zeros below) and its tangents:
`dUi = -Ui triu_half(Ui' dS Ui)`, `triu_half` the upper triangle with the
diagonal halved. `X` and `T` are scratch of `S`'s shape. `false` where `S`
does not factor.
"""
function _tan_cholinv!(Ui::TanMat{W}, S::TanMat{W}, X::TanMat{W}, T::TanMat{W}, F::Matrix{Float64}) where {W}
    m = size(S.v, 1)
    copyto!(F, S.v)
    _ctsem_cholesky!(F, m) || return false
    @inbounds for j in 1:m, i in (j + 1):m
        F[i, j] = 0.0
    end
    U = Ui.v
    _ctsem_cholesky_uinv!(U, CTSEMCholesky(F, m, true))
    if !_tan_nz(S)
        _tan_setnz!(Ui, false)
        return true
    end
    dS, dT, dX = S.d, T.d, X.d
    V = Val(W)
    # T = dS Ui, then X = Ui' T, per direction; Ui is upper triangular.
    @inbounds for j in 1:m, i in 1:m
        acc = _tan_zeros(V)
        for q in 1:j
            acc = _tan_axpy(acc, U[q, j], dS, i, q)
        end
        _tan_put!(T, i, j, acc)
    end
    @inbounds for j in 1:m, i in 1:j
        acc = _tan_zeros(V)
        for q in 1:i
            acc = _tan_axpy(acc, U[q, i], dT, q, j)
        end
        _tan_put!(X, i, j, i == j ? _tan_times(0.5, acc) : acc)
    end
    # dUi = -Ui triu_half(X): upper triangular, rows of Ui from the diagonal.
    @inbounds for j in 1:m, i in 1:m
        acc = _tan_zeros(V)
        for q in i:j
            acc = _tan_axpy(acc, -U[i, q], dX, q, j)
        end
        _tan_put!(Ui, i, j, acc)
    end
    Ui.nz[] = true
    return true
end
