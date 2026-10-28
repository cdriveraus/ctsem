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
#
# A quantity built under `_TAN_JETS` (a jet: `jet`) also carries `s`, the sum
# over the W directions of its second derivatives along them -- what a
# `CTSEMTraceJet` carries beside its first parts (`trace_jet.jl`), so a sweep
# of jets is that number's forward Laplacian without the filter at its type:
#
#     s(A B) = s(A) B + A s(B) + 2 sum_w dA_w dB_w
#
# and a linear operation's is its own on `s`. `snz` says whether `s` may be
# nonzero, as `nz` does for the tangents. Each primitive computes its output's
# `s` where the output is a jet; one that is not, handed an operand whose `s`
# is live, throws rather than drop it. The second parts are each a function of
# their own, never inlined: inlined, their locals enlarged every primitive's
# frame, which is cleared at each call, and the curvature of a two-latent model,
# which carries no second parts, ran 60 per cent slower, a seventh of it in
# that clearing.

"""Whether quantities built now carry second parts (`TanMat.jet`): set around
building and running a jet sweep (`with(_TAN_JETS => true)`)."""
const _TAN_JETS = Base.ScopedValues.ScopedValue(false)

struct TanMat{W}
    v::Matrix{Float64}
    d::Array{Float64,3}     # W x r x c
    nz::Base.RefValue{Bool}
    s::Matrix{Float64}      # r x c where a jet, 0 x 0 otherwise
    snz::Base.RefValue{Bool}
    jet::Bool
end
function TanMat{W}(r::Int, c::Int) where {W}
    j = _TAN_JETS[]
    return TanMat{W}(zeros(r, c), zeros(W, r, c), Ref(true), j ? zeros(r, c) : zeros(0, 0),
        Ref(false), j)
end

struct TanVec{W}
    v::Vector{Float64}
    d::Matrix{Float64}      # W x len
    nz::Base.RefValue{Bool}
    s::Vector{Float64}
    snz::Base.RefValue{Bool}
    jet::Bool
end
function TanVec{W}(r::Int) where {W}
    j = _TAN_JETS[]
    return TanVec{W}(zeros(r), zeros(W, r), Ref(true), j ? zeros(r) : zeros(0), Ref(false), j)
end

@inline _tan_nz(A) = A.nz[]
"""Mark whether `A`'s tangents may be nonzero, zeroing them where not."""
@inline function _tan_setnz!(A, on::Bool)
    on || (A.nz[] && fill!(A.d, 0.0))
    A.nz[] = on
    return A
end
"""Mark `A` from its tangents themselves (an input)."""
_tan_setnz_from!(A) = (A.nz[] = any(!iszero, A.d); A)

@inline _tan_snz(A) = A.snz[]
"""Mark whether `A`'s second parts may be nonzero, zeroing them where not."""
@inline function _tan_setsnz!(A, on::Bool)
    on || (A.snz[] && fill!(A.s, 0.0))
    A.snz[] = on
    return A
end
"""Mark `A` from its second parts themselves (an input)."""
_tan_setsnz_from!(A) = (A.snz[] = A.jet && any(!iszero, A.s); A)

@noinline _tan_jet_mismatch() = throw(ArgumentError(string("a quantity with ",
    "second parts written into one without: a jet sweep's scratch built outside ",
    "`_TAN_JETS`")))
"""Throws where `C` is no jet but `live` says an operand's second parts are."""
@inline _tan_jet_check(C, live::Bool) = (!C.jet && live && _tan_jet_mismatch(); nothing)

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
    _tan_jet_check(C, _tan_snz(A))
    C.jet && _tan_copy_s!(C, A)
    return C
end
@noinline function _tan_copy_s!(C::Union{TanMat,TanVec}, A::Union{TanMat,TanVec})
    if _tan_snz(A)
        _tan_copyto!(C.s, A.s)
        C.snz[] = true
    else
        _tan_setsnz!(C, false)
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
# and last directions scalar. A second-parts matrix is a width-one array of
# the same shape, its entry (i, j) at offset `(i - 1) + r (j - 1)`.
@inline _tan_off(d::Array{Float64,3}, ::Val{W}, i, j) where {W} =
    W * ((i - 1) + size(d, 2) * (j - 1))
@inline _tan_off(d::Matrix{Float64}, ::Val{W}, i) where {W} = W * (i - 1)
@inline _tan_off(d::Matrix{Float64}, ::Val{1}, i, j) = (i - 1) + size(d, 1) * (j - 1)
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
"""`sum_w dA[oa + w] dB[ob + w]`: two tangent columns' product summed over the
directions, a second part's cross term."""
@inline _tan_dot(dA, oa, dB, ob, ::Val{W}) where {W} =
    sum(ntuple(w -> @inbounds(dA[oa + w] * dB[ob + w]), Val(W)))

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
"""Copy `C`'s second parts' upper triangle onto its lower."""
function _tan_mirror_s!(C::TanMat)
    n = size(C.s, 1)
    @inbounds for j in 1:n, i in (j + 1):n
        C.s[i, j] = C.s[j, i]
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

"""
The product's second parts: `s(C) = alpha (op(sA) op(B) + op(A) op(sB) +
2 sum_w op(dA_w) op(dB_w)) + beta s(C)`, each term where its operands carry it,
from the operands as the product reads them -- so before the product where
`C` is one of them. A triangular factor's second parts are triangular as its
tangents are.
"""
@noinline function _tan_prod_s!(C::TanMat{W}, A::TanMat{W}, B::TanMat{W}, alpha::Float64,
    beta::Float64, tri::Int, sym::Bool, ta::Bool, tb::Bool) where {W}
    sa, sb, x = _tan_snz(A), _tan_snz(B), _tan_nz(A) && _tan_nz(B)
    keep = !iszero(beta) && _tan_snz(C)
    if !(sa || sb || x)
        keep ? (isone(beta) || (C.s .*= beta)) : _tan_setsnz!(C, false)
        return C
    end
    lda, ldb, ldc = size(A.v, 1), size(B.v, 1), size(C.v, 1)
    r = ta ? size(A.v, 2) : lda
    k = ta ? lda : size(A.v, 2)
    c = tb ? ldb : size(B.v, 2)
    sai, saq = ta ? (lda, 1) : (1, lda)
    sbq, sbj = tb ? (ldb, 1) : (1, ldb)
    Av, Bv, As, Bs, Cs, dA, dB = A.v, B.v, A.s, B.s, C.s, A.d, B.d
    @inbounds for j in 1:c, i in 1:(sym ? j : r)
        q0, q1 = _tan_qbounds(tri, i, j, k)
        acc = 0.0
        for q in q0:q1
            la = (i - 1) * sai + (q - 1) * saq
            lb = (q - 1) * sbq + (j - 1) * sbj
            sa && (acc = muladd(As[la + 1], Bv[lb + 1], acc))
            sb && (acc = muladd(Av[la + 1], Bs[lb + 1], acc))
            x && (acc = muladd(2.0, _tan_dot(dA, W * la, dB, W * lb, Val(W)), acc))
        end
        lc = (i - 1) + (j - 1) * ldc
        Cs[lc + 1] = keep ? muladd(alpha, acc, beta * Cs[lc + 1]) : alpha * acc
    end
    sym && _tan_mirror_s!(C)
    C.snz[] = true
    return C
end

function _tan_prod_flags!(C::TanMat{W}, A::TanMat{W}, B::TanMat{W}, alpha::Float64, beta::Float64,
    tri::Int, sym::Bool, ta::Bool, tb::Bool) where {W}
    fa, fb = _tan_nz(A), _tan_nz(B)
    _tan_jet_check(C, _tan_snz(A) || _tan_snz(B))
    C.jet && _tan_prod_s!(C, A, B, alpha, beta, tri, sym, ta, tb)
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

"""The matrix-vector product's second parts, as `_tan_prod_s!`'s."""
@noinline function _tan_prodvec_s!(y::TanVec{W}, A::TanMat{W}, x::TanVec{W}, alpha::Float64,
    beta::Float64, tri::Int, ta::Bool) where {W}
    sa, sx, c = _tan_snz(A), _tan_snz(x), _tan_nz(A) && _tan_nz(x)
    keep = !iszero(beta) && _tan_snz(y)
    if !(sa || sx || c)
        keep ? (isone(beta) || (y.s .*= beta)) : _tan_setsnz!(y, false)
        return y
    end
    lda = size(A.v, 1)
    r = ta ? size(A.v, 2) : lda
    k = ta ? lda : size(A.v, 2)
    sai, saq = ta ? (lda, 1) : (1, lda)
    Av, xv, As, xs, ys, dA, dx = A.v, x.v, A.s, x.s, y.s, A.d, x.d
    @inbounds for i in 1:r
        q0, q1 = _tan_qbounds(tri, i, 1, k)
        acc = 0.0
        for q in q0:q1
            la = (i - 1) * sai + (q - 1) * saq
            sa && (acc = muladd(As[la + 1], xv[q], acc))
            sx && (acc = muladd(Av[la + 1], xs[q], acc))
            c && (acc = muladd(2.0, _tan_dot(dA, W * la, dx, W * (q - 1), Val(W)), acc))
        end
        ys[i] = keep ? muladd(alpha, acc, beta * ys[i]) : alpha * acc
    end
    y.snz[] = true
    return y
end

function _tan_prodvec_flags!(y::TanVec{W}, A::TanMat{W}, x::TanVec{W}, alpha::Float64, beta::Float64,
    tri::Int, ta::Bool) where {W}
    fa, fx = _tan_nz(A), _tan_nz(x)
    _tan_jet_check(y, _tan_snz(A) || _tan_snz(x))
    y.jet && _tan_prodvec_s!(y, A, x, alpha, beta, tri, ta)
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
    _tan_jet_check(S, _tan_snz(A))
    S.jet && _tan_sym_s!(S, A)
    _tan_nz(A) || return _tan_setnz!(S, false)
    @inbounds for j in 1:n, i in 1:j
        t = _tan_axpby(0.5, _tan_col(A, i, j), 0.5, _tan_col(A, j, i))
        _tan_put!(S, i, j, t)
        _tan_put!(S, j, i, t)
    end
    S.nz[] = true
    return S
end

@noinline function _tan_sym_s!(S::TanMat, A::TanMat)
    _tan_snz(A) || return _tan_setsnz!(S, false)
    n = size(A.v, 1)
    As, Ss = A.s, S.s
    @inbounds for j in 1:n, i in 1:j
        s = (As[i, j] + As[j, i]) / 2
        Ss[i, j] = s
        Ss[j, i] = s
    end
    S.snz[] = true
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
    _tan_jet_check(C, _tan_snz(A))
    C.jet && _tan_add_s!(C, alpha, A, beta, tr)
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

@noinline function _tan_add_s!(C::TanMat, alpha::Float64, A::TanMat, beta::Float64,
    tr::Bool)
    r, c = size(C.v)
    lda = size(A.v, 1)
    sai, saj = tr ? (lda, 1) : (1, lda)
    keep = !iszero(beta) && _tan_snz(C)
    if _tan_snz(A)
        As, Cs = A.s, C.s
        @inbounds for j in 1:c, i in 1:r
            la = (i - 1) * sai + (j - 1) * saj
            Cs[i, j] = keep ? muladd(alpha, As[la + 1], beta * Cs[i, j]) : alpha * As[la + 1]
        end
        C.snz[] = true
    elseif !keep
        _tan_setsnz!(C, false)
    elseif !isone(beta)
        C.s .*= beta
    end
    return C
end

"""C = alpha x y' + beta C, with tangents."""
@inline _tan_ger!(C::TanMat, alpha, x::TanVec, y::TanVec, beta=true) =
    _tan_ger_flags!(C, Float64(alpha), x, y, Float64(beta))
function _tan_ger_flags!(C::TanMat{W}, alpha::Float64, x::TanVec{W}, y::TanVec{W}, beta::Float64) where {W}
    xv, yv, Cv = x.v, y.v, C.v
    bz = iszero(beta)
    _tan_jet_check(C, _tan_snz(x) || _tan_snz(y))
    C.jet && _tan_ger_s!(C, alpha, x, y, beta)
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
"""The outer product's second parts: `alpha (s(x) y' + x s(y)' + 2 sum_w dx_w
dy_w') + beta s(C)`."""
@noinline function _tan_ger_s!(C::TanMat{W}, alpha::Float64, x::TanVec{W}, y::TanVec{W},
    beta::Float64) where {W}
    sx, sy, c = _tan_snz(x), _tan_snz(y), _tan_nz(x) && _tan_nz(y)
    keep = !iszero(beta) && _tan_snz(C)
    if !(sx || sy || c)
        keep ? (isone(beta) || (C.s .*= beta)) : _tan_setsnz!(C, false)
        return C
    end
    xv, yv, xs, ys, Cs, dx, dy = x.v, y.v, x.s, y.s, C.s, x.d, y.d
    @inbounds for j in eachindex(yv), i in eachindex(xv)
        acc = 0.0
        sx && (acc = muladd(xs[i], yv[j], acc))
        sy && (acc = muladd(xv[i], ys[j], acc))
        c && (acc = muladd(2.0, _tan_dot(dx, W * (i - 1), dy, W * (j - 1), Val(W)), acc))
        Cs[i, j] = keep ? muladd(alpha, acc, beta * Cs[i, j]) : alpha * acc
    end
    C.snz[] = true
    return C
end

"""A *= s, with tangents."""
function _tan_scale!(A::Union{TanMat,TanVec}, s)
    A.v .*= s
    _tan_nz(A) && (A.d .*= s)
    _tan_snz(A) && (A.s .*= s)
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
    _tan_jet_check(C, _tan_comb_snz(terms))
    C.jet && _tan_comb_s!(C, terms, Float64(beta))
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
@noinline function _tan_comb_s!(C::TanMat, terms::Tuple, beta::Float64)
    r, c = size(C.v)
    keep = !iszero(beta) && _tan_snz(C)
    if _tan_comb_snz(terms)
        Cs = C.s
        @inbounds for j in 1:c, i in 1:r
            Cs[i, j] = _tan_comb_sval(keep ? beta * Cs[i, j] : 0.0, terms, i, j)
        end
        C.snz[] = true
    elseif !keep
        _tan_setsnz!(C, false)
    elseif !isone(beta)
        C.s .*= beta
    end
    return C
end
@inline _tan_comb_val(v, ::Tuple{}, i, j) = v
@inline _tan_comb_val(v, terms::Tuple, i, j) =
    _tan_comb_val(muladd(terms[1][1], @inbounds(terms[1][2].v[i, j]), v), Base.tail(terms), i, j)
@inline _tan_comb_nz(::Tuple{}) = false
@inline _tan_comb_nz(terms::Tuple) = _tan_nz(terms[1][2]) || _tan_comb_nz(Base.tail(terms))
@inline _tan_comb_snz(::Tuple{}) = false
@inline _tan_comb_snz(terms::Tuple) = _tan_snz(terms[1][2]) || _tan_comb_snz(Base.tail(terms))
@inline _tan_comb_sval(v, ::Tuple{}, i, j) = v
@inline function _tan_comb_sval(v, terms::Tuple, i, j)
    a, A = terms[1]
    return _tan_comb_sval(_tan_snz(A) ? muladd(a, @inbounds(A.s[i, j]), v) : v,
        Base.tail(terms), i, j)
end
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
directions riding every operation of the elimination. A jet's second parts,
`s(X) = Q^-1 (s(B) - s(Q) X - 2 sum_w dQ_w dX_w)`, on the same factors.
"""
function _tan_solve!(X::Union{TanMat{W},TanVec{W}}, F::Matrix{Float64}, piv::Vector{Int},
    Q::TanMat{W}) where {W}
    n = size(X.v, 1)
    r = size(X.v, 2)
    _lu_solve_generic!(F, piv, X.v, n)
    fq = _tan_nz(Q)
    _tan_jet_check(X, _tan_snz(Q))
    if fq || _tan_nz(X)
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
    end
    X.jet && _tan_solve_s!(X, F, piv, Q)
    return X
end

"""`_tan_solve!`'s second parts, after its tangents: `X.s` holds the
right-hand side's, `X.d` the solution's tangents."""
@noinline function _tan_solve_s!(X::Union{TanMat{W},TanVec{W}}, F::Matrix{Float64},
    piv::Vector{Int}, Q::TanMat{W}) where {W}
    n = size(X.v, 1)
    r = size(X.v, 2)
    sb, sq, c = _tan_snz(X), _tan_snz(Q), _tan_nz(Q) && _tan_nz(X)
    (sb || sq || c) || return _tan_setsnz!(X, false)
    Xs, Xv, Qs, dQ, dX = X.s, X.v, Q.s, Q.d, X.d
    sb || fill!(Xs, 0.0)
    @inbounds for cc in 1:r, i in 1:n
        acc = Xs[i + n * (cc - 1)]
        for q in 1:n
            sq && (acc -= Qs[i, q] * Xv[q + n * (cc - 1)])
            c && (acc -= 2 * _tan_dot(dQ, W * ((i - 1) + n * (q - 1)), dX,
                W * ((q - 1) + n * (cc - 1)), Val(W)))
        end
        Xs[i + n * (cc - 1)] = acc
    end
    _lu_solve_generic!(F, piv, Xs, n)
    X.snz[] = true
    return X
end

"""Copy C's lower triangle onto its upper, tangents too, as the filter does."""
function _tan_lower_to_upper!(C::TanMat{W}) where {W}
    n = size(C.v, 1)
    on = _tan_nz(C)
    son = _tan_snz(C)
    @inbounds for j in 1:n, i in (j + 1):n
        C.v[j, i] = C.v[i, j]
        on && _tan_put!(C, j, i, _tan_col(C, i, j))
        son && (C.s[j, i] = C.s[i, j])
    end
    return C
end

"""C = A[rows, cols], with tangents."""
function _tan_gather!(C::TanMat{W}, A::TanMat{W}, rows::AbstractVector{Int},
    cols::AbstractVector{Int}) where {W}
    on = _tan_nz(A)
    son = _tan_snz(A)
    _tan_jet_check(C, son)
    @inbounds for (jj, j) in enumerate(cols), (ii, i) in enumerate(rows)
        C.v[ii, jj] = A.v[i, j]
        on && _tan_put!(C, ii, jj, _tan_col(A, i, j))
        son && (C.s[ii, jj] = A.s[i, j])
    end
    on ? (C.nz[] = true) : _tan_setnz!(C, false)
    son ? (C.snz[] = true) : _tan_setsnz!(C, false)
    return C
end
"""c = a[rows], with tangents."""
function _tan_gather!(c::TanVec{W}, a::TanVec{W}, rows::AbstractVector{Int}) where {W}
    on = _tan_nz(a)
    son = _tan_snz(a)
    _tan_jet_check(c, son)
    @inbounds for (ii, i) in enumerate(rows)
        c.v[ii] = a.v[i]
        on && _tan_put!(c, ii, _tan_col(a, i))
        son && (c.s[ii] = a.s[i])
    end
    on ? (c.nz[] = true) : _tan_setnz!(c, false)
    son ? (c.snz[] = true) : _tan_setsnz!(c, false)
    return c
end
"""A[rows, cols] += alpha C, with tangents."""
function _tan_scatter!(A::TanMat{W}, alpha::Float64, C::TanMat{W}, rows::AbstractVector{Int},
    cols::AbstractVector{Int}) where {W}
    on = _tan_nz(C)
    son = _tan_snz(C)
    _tan_jet_check(A, son)
    @inbounds for (jj, j) in enumerate(cols), (ii, i) in enumerate(rows)
        A.v[i, j] = muladd(alpha, C.v[ii, jj], A.v[i, j])
        on && _tan_put!(A, i, j, _tan_axpby(1.0, _tan_col(A, i, j), alpha, _tan_col(C, ii, jj)))
        son && (A.s[i, j] = muladd(alpha, C.s[ii, jj], A.s[i, j]))
    end
    on && (A.nz[] = true)
    son && (A.snz[] = true)
    return A
end
"""a[rows] += alpha c, with tangents."""
function _tan_scatter!(a::TanVec{W}, alpha::Float64, c::TanVec{W}, rows::AbstractVector{Int}) where {W}
    on = _tan_nz(c)
    son = _tan_snz(c)
    _tan_jet_check(a, son)
    @inbounds for (ii, i) in enumerate(rows)
        a.v[i] = muladd(alpha, c.v[ii], a.v[i])
        on && _tan_put!(a, i, _tan_axpby(1.0, _tan_col(a, i), alpha, _tan_col(c, ii)))
        son && (a.s[i] = muladd(alpha, c.s[ii], a.s[i]))
    end
    on && (a.nz[] = true)
    son && (a.snz[] = true)
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
    _tan_jet_check(z, _tan_snz(x) || _tan_snz(y))
    # The second parts first: `z` may be `x` or `y`.
    z.jet && _tan_lincomb_s!(z, a, x, b, y)
    @inbounds for i in eachindex(zv)
        zv[i] = muladd(a, xv[i], b * yv[i])
        on && _tan_put!(z, i, _tan_axpby(a, _tan_col(x, i), b, _tan_col(y, i)))
    end
    on ? (z.nz[] = true) : _tan_setnz!(z, false)
    return z
end
@noinline function _tan_lincomb_s!(z::TanVec, a::Float64, x::TanVec, b::Float64, y::TanVec)
    sx, sy = _tan_snz(x), _tan_snz(y)
    (sx || sy) || return _tan_setsnz!(z, false)
    @inbounds for i in eachindex(z.s)
        z.s[i] = muladd(a, sx ? x.s[i] : 0.0, b * (sy ? y.s[i] : 0.0))
    end
    z.snz[] = true
    return z
end

"""
`Ui = U^-1` for `S = U'U` (upper, zeros below) and its tangents:
`dUi = -Ui triu_half(Ui' dS Ui)`, `triu_half` the upper triangle with the
diagonal halved. `X` and `T` are scratch of `S`'s shape. `false` where `S`
does not factor. A jet's second parts follow from `Ui' S Ui = I`
differentiated twice: with `Z_w = Ui^-1 dUi_w = -triu_half(Ui' dS_w Ui)`,
`s(Ui) = -Ui triu_half(Y)`, `Y = Ui' s(S) Ui - 2 sum_w (Z_w' Z_w + Z_w Z_w +
(Z_w Z_w)')`.
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
    _tan_jet_check(Ui, _tan_snz(S))
    if !_tan_nz(S)
        _tan_setnz!(Ui, false)
        Ui.jet && _tan_cholinv_s!(Ui, S, X, T, false)
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
    Ui.jet && _tan_cholinv_s!(Ui, S, X, T, true)
    return true
end

"""
`_tan_cholinv!`'s second parts, after its tangents: `X.d` holds `triu_half(Ui'
dS_w Ui)` on and above its diagonal (`-Z_w` there; `Z_w` is zero below), and
`X.s`, `T.s` are scratch. `cross` where the tangents are live. Each term of
`Y`'s sum is a product of two `Z_w`, so `X.d` serves for `Z_w` as it stands.
"""
@noinline function _tan_cholinv_s!(Ui::TanMat{W}, S::TanMat{W}, X::TanMat{W}, T::TanMat{W},
    cross::Bool) where {W}
    m = size(S.v, 1)
    ss = _tan_snz(S)
    (ss || cross) || return _tan_setsnz!(Ui, false)
    U, Y, Ts, dX = Ui.v, X.s, T.s, X.d
    V = Val(W)
    # Y = Ui' s(S) Ui, its upper triangle.
    @inbounds for j in 1:m, i in 1:j
        Y[i, j] = 0.0
    end
    if ss
        Ss = S.s
        @inbounds for j in 1:m, i in 1:m
            acc = 0.0
            for q in 1:j
                acc = muladd(Ss[i, q], U[q, j], acc)
            end
            Ts[i, j] = acc
        end
        @inbounds for j in 1:m, i in 1:j
            acc = 0.0
            for q in 1:i
                acc = muladd(U[q, i], Ts[q, j], acc)
            end
            Y[i, j] = acc
        end
    end
    # - 2 sum_w (Z' Z + Z Z + (Z Z)'), on the upper triangle: Z upper, so
    # (Z' Z)[i, j] runs over q <= i, (Z Z)[i, j] over i <= q <= j, and
    # (Z Z)'[i, j] above the diagonal is zero.
    if cross
        @inbounds for j in 1:m, i in 1:j
            acc = 0.0
            for q in 1:i
                acc += _tan_dot(dX, _tan_off(dX, V, q, i), dX, _tan_off(dX, V, q, j), V)
            end
            for q in i:j
                acc += _tan_dot(dX, _tan_off(dX, V, i, q), dX, _tan_off(dX, V, q, j), V)
            end
            i == j && (acc += _tan_dot(dX, _tan_off(dX, V, i, i), dX,
                _tan_off(dX, V, i, i), V))
            Y[i, j] -= 2 * acc
        end
    end
    # s(Ui) = Ui Z2, Z2 = -triu_half(Y).
    Uis = Ui.s
    @inbounds for j in 1:m, i in 1:m
        acc = 0.0
        for q in i:j
            acc -= U[i, q] * (q == j ? Y[q, j] / 2 : Y[q, j])
        end
        Uis[i, j] = acc
    end
    Ui.snz[] = true
    return Ui
end
