################################################################################
# A number carrying the trace of second directional derivatives
################################################################################
#
# The Laplace term's gradient needs, for every parameter, the sum over a
# block's directions `q_j` of the second directional derivative of the
# gradient: `sum_j d2(grad f)[q_j, q_j]` (`_laplace_seeded_unit_gradient!`).
# A hyper-dual sweep gives one direction's, with a value, two first-order parts
# that are equal and the second-order one, so the sum took one sweep of the
# whole subject a direction.
#
# The sum needs less. For `x = x0 + sum_j t_j q_j`, write `x'_j` for the
# derivative along `t_j` and `S_x` for `sum_j` of the second ones. Then
#
#     S_{ab}   = S_a b + a S_b + 2 sum_j a'_j b'_j
#     S_{g(a)} = g''(a) sum_j a'_j^2 + g'(a) S_a
#
# so a number that holds the value, the `W` first parts and the one summed
# second part closes under every operation (the "forward Laplacian"). Run
# through the engine's reverse pass, it returns every parameter's gradient
# with that sum in one sweep of `W` directions: `W + 2` components, about
# `4W + 3` multiplies a product against `9W` for `W` hyper-dual sweeps, and one
# sweep's fixed cost instead of `W`.
#
# Unary functions without a rule of their own take `g'` and `g''` from one
# evaluation at a nested scalar dual (`_jet_derivatives`), so anything the
# engine evaluates at ForwardDiff's duals works here.

"""
    CTSEMTraceJet{W}

The value `v`, the derivatives `d` along `W` directions and `s`, the sum of
the second derivatives along them (see the note above). Equality, `iszero` and
finiteness read every component, as the engine's caches and validity checks
read a dual's.
"""
struct CTSEMTraceJet{W} <: Real
    v::Float64
    d::NTuple{W,Float64}
    s::Float64
end

"""The width every jet sweep runs at: one width, so one compiled type."""
const _CTSEM_JET_WIDTH = 8

const _JetReal = Union{AbstractFloat,Integer,Rational,AbstractIrrational}
# What the mixed operations take directly. An irrational reaches them by
# promotion, and `Bool` has its own `*`: Base's `Bool * x` methods and its
# irrationals' would otherwise be ambiguous with these.
const _JetMixed = Union{AbstractFloat,Integer,Rational}

@inline _jet_zeros(::Val{W}) where {W} = ntuple(_ -> 0.0, Val(W))
@inline CTSEMTraceJet{W}(x::_JetReal) where {W} =
    CTSEMTraceJet{W}(Float64(x), _jet_zeros(Val(W)), 0.0)
@inline CTSEMTraceJet{W}(x::CTSEMTraceJet{W}) where {W} = x

Base.convert(::Type{CTSEMTraceJet{W}}, x::CTSEMTraceJet{W}) where {W} = x
Base.convert(::Type{CTSEMTraceJet{W}}, x::_JetReal) where {W} = CTSEMTraceJet{W}(x)
Base.promote_rule(::Type{CTSEMTraceJet{W}}, ::Type{<:_JetReal}) where {W} =
    CTSEMTraceJet{W}

Base.zero(::Type{CTSEMTraceJet{W}}) where {W} = CTSEMTraceJet{W}(0.0)
Base.one(::Type{CTSEMTraceJet{W}}) where {W} = CTSEMTraceJet{W}(1.0)
Base.zero(x::CTSEMTraceJet) = zero(typeof(x))
Base.one(x::CTSEMTraceJet) = one(typeof(x))
Base.float(x::CTSEMTraceJet) = x
Base.eps(::Type{<:CTSEMTraceJet}) = eps(Float64)
Base.eps(x::CTSEMTraceJet) = eps(x.v)
Base.floatmin(::Type{<:CTSEMTraceJet}) = floatmin(Float64)
Base.floatmax(::Type{<:CTSEMTraceJet}) = floatmax(Float64)
Base.typemin(::Type{CTSEMTraceJet{W}}) where {W} = CTSEMTraceJet{W}(-Inf)
Base.typemax(::Type{CTSEMTraceJet{W}}) where {W} = CTSEMTraceJet{W}(Inf)
Base.rtoldefault(::Type{<:CTSEMTraceJet}) = Base.rtoldefault(Float64)
ForwardDiff.value(x::CTSEMTraceJet) = x.v

Base.isfinite(x::CTSEMTraceJet) = isfinite(x.v)
Base.isnan(x::CTSEMTraceJet) = isnan(x.v)
Base.isinf(x::CTSEMTraceJet) = isinf(x.v)
Base.signbit(x::CTSEMTraceJet) = signbit(x.v)
Base.isinteger(x::CTSEMTraceJet) = isinteger(x.v) && _jet_constant(x)
Base.iszero(x::CTSEMTraceJet) = iszero(x.v) && _jet_constant(x)
Base.isone(x::CTSEMTraceJet) = isone(x.v) && _jet_constant(x)
@inline _jet_constant(x::CTSEMTraceJet) = iszero(x.s) && all(iszero, x.d)

for f in (:floor, :ceil, :round, :trunc)
    @eval Base.$f(::Type{I}, x::CTSEMTraceJet) where {I<:Integer} = $f(I, x.v)
end
Base.Int(x::CTSEMTraceJet) = Int(x.v)

Base.:(==)(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} =
    a.v == b.v && a.d == b.d && a.s == b.s
Base.:(==)(a::CTSEMTraceJet, b::_JetMixed) = a.v == b && _jet_constant(a)
Base.:(==)(a::_JetMixed, b::CTSEMTraceJet) = b == a
Base.isequal(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} =
    isequal(a.v, b.v) && isequal(a.d, b.d) && isequal(a.s, b.s)
Base.hash(x::CTSEMTraceJet, h::UInt) = hash(x.s, hash(x.d, hash(x.v, h)))
for op in (:<, :<=, :isless)
    @eval begin
        Base.$op(a::CTSEMTraceJet, b::CTSEMTraceJet) = $op(a.v, b.v)
        Base.$op(a::CTSEMTraceJet, b::_JetMixed) = $op(a.v, b)
        Base.$op(a::_JetMixed, b::CTSEMTraceJet) = $op(a, b.v)
    end
end

# A pairwise sum: as a running sum the dot product was one chain of `W`
# dependent multiply-adds in every jet product, and the jet's products ran at a
# seventh of a Dual's speed a component (`bench_jetmul.jl`).
@inline _jet_tree(t::NTuple{1,Float64}) = @inbounds t[1]
@inline _jet_tree(t::NTuple{2,Float64}) = @inbounds t[1] + t[2]
@inline function _jet_tree(t::NTuple{N,Float64}) where {N}
    h = N >> 1
    return _jet_tree(ntuple(i -> @inbounds(t[i]), Val(h))) +
        _jet_tree(ntuple(i -> @inbounds(t[h + i]), Val(N - h)))
end
@inline _jet_dot(a::NTuple{W,Float64}, b::NTuple{W,Float64}) where {W} =
    _jet_tree(ntuple(i -> @inbounds(a[i] * b[i]), Val(W)))

# g(x) given g, g' and g'' at the value.
@inline _jet_chain(x::CTSEMTraceJet{W}, f0, f1, f2) where {W} =
    CTSEMTraceJet{W}(f0, ntuple(i -> f1 * x.d[i], Val(W)),
        f2 * _jet_dot(x.d, x.d) + f1 * x.s)

Base.:+(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} =
    CTSEMTraceJet{W}(a.v + b.v, ntuple(i -> a.d[i] + b.d[i], Val(W)), a.s + b.s)
Base.:-(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} =
    CTSEMTraceJet{W}(a.v - b.v, ntuple(i -> a.d[i] - b.d[i], Val(W)), a.s - b.s)
Base.:-(a::CTSEMTraceJet{W}) where {W} =
    CTSEMTraceJet{W}(-a.v, ntuple(i -> -a.d[i], Val(W)), -a.s)
Base.:+(a::CTSEMTraceJet) = a
Base.:*(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} =
    CTSEMTraceJet{W}(a.v * b.v, ntuple(i -> muladd(a.v, b.d[i], b.v * a.d[i]), Val(W)),
        muladd(a.s, b.v, muladd(a.v, b.s, 2 * _jet_dot(a.d, b.d))))
function Base.:/(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W}
    q = a.v / b.v
    d = ntuple(i -> (a.d[i] - q * b.d[i]) / b.v, Val(W))
    return CTSEMTraceJet{W}(q, d, (a.s - q * b.s - 2 * _jet_dot(d, b.d)) / b.v)
end

Base.:+(a::CTSEMTraceJet{W}, r::_JetMixed) where {W} = CTSEMTraceJet{W}(a.v + r, a.d, a.s)
Base.:+(r::_JetMixed, a::CTSEMTraceJet) = a + r
Base.:-(a::CTSEMTraceJet{W}, r::_JetMixed) where {W} = CTSEMTraceJet{W}(a.v - r, a.d, a.s)
Base.:-(r::_JetMixed, a::CTSEMTraceJet) = -a + r
Base.:*(a::CTSEMTraceJet{W}, r::_JetMixed) where {W} =
    CTSEMTraceJet{W}(a.v * r, ntuple(i -> a.d[i] * r, Val(W)), a.s * r)
Base.:*(r::_JetMixed, a::CTSEMTraceJet) = a * r
Base.:/(a::CTSEMTraceJet{W}, r::_JetMixed) where {W} =
    CTSEMTraceJet{W}(a.v / r, ntuple(i -> a.d[i] / r, Val(W)), a.s / r)
Base.:/(r::_JetMixed, a::CTSEMTraceJet) = r * inv(a)

Base.:*(b::Bool, a::CTSEMTraceJet) = b ? a : zero(a)
Base.:*(a::CTSEMTraceJet, b::Bool) = b * a

# Fused, as ForwardDiff's is: the engine's product kernels are nests of these.
@inline function Base.muladd(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W},
    c::CTSEMTraceJet{W}) where {W}
    av = a.v
    bv = b.v
    return CTSEMTraceJet{W}(muladd(av, bv, c.v),
        ntuple(i -> @inbounds(muladd(av, b.d[i], muladd(bv, a.d[i], c.d[i]))), Val(W)),
        muladd(a.s, bv, muladd(av, b.s, muladd(2.0, _jet_dot(a.d, b.d), c.s))))
end
Base.fma(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}, c::CTSEMTraceJet{W}) where {W} =
    muladd(a, b, c)

function Base.inv(x::CTSEMTraceJet)
    i = inv(x.v)
    return _jet_chain(x, i, -i * i, 2 * i * i * i)
end
function Base.exp(x::CTSEMTraceJet)
    e = exp(x.v)
    return _jet_chain(x, e, e, e)
end
function Base.expm1(x::CTSEMTraceJet)
    e = exp(x.v)
    return _jet_chain(x, expm1(x.v), e, e)
end
Base.log(x::CTSEMTraceJet) = _jet_chain(x, log(x.v), inv(x.v), -inv(x.v)^2)
Base.log1p(x::CTSEMTraceJet) = _jet_chain(x, log1p(x.v), inv(1 + x.v), -inv(1 + x.v)^2)
function Base.sqrt(x::CTSEMTraceJet)
    r = sqrt(x.v)
    return _jet_chain(x, r, 0.5 / r, -0.25 / (r * x.v))
end
Base.abs(x::CTSEMTraceJet) = signbit(x.v) ? -x : x
Base.abs2(x::CTSEMTraceJet) = x * x
function Base.:^(x::CTSEMTraceJet, n::Integer)
    n == 0 && return one(x)
    n == 1 && return x
    n == 2 && return x * x
    v = x.v
    return _jet_chain(x, v^n, n * v^(n - 1), n * (n - 1) * v^(n - 2))
end
function Base.:^(x::CTSEMTraceJet, p::AbstractFloat)
    isinteger(p) && abs(p) < 64 && return x^Int(p)
    v = x.v
    return _jet_chain(x, v^p, p * v^(p - 1), p * (p - 1) * v^(p - 2))
end
Base.:^(x::CTSEMTraceJet, p::Rational) = x^Float64(p)
Base.:^(x::CTSEMTraceJet, y::CTSEMTraceJet) = exp(y * log(x))
Base.:^(r::_JetMixed, y::CTSEMTraceJet) = exp(y * log(r))
# Literal powers (`x^2` in an expression) reach `^` through `literal_pow`.
Base.literal_pow(::typeof(^), x::CTSEMTraceJet, ::Val{p}) where {p} = x^p

# The rest from one evaluation at a nested scalar dual: `g(v + e1 + e2)` holds
# `g`, `g'` (either first part) and `g''` (the `e1 e2` part).
struct _JetTagInner end
struct _JetTagOuter end
ForwardDiff.:≺(::Type{_JetTagInner}, ::Type{_JetTagOuter}) = true
ForwardDiff.:≺(::Type{_JetTagOuter}, ::Type{_JetTagInner}) = false

@inline function _jet_derivatives(g::F, v::Float64) where {F}
    z = ForwardDiff.Dual{_JetTagOuter}(ForwardDiff.Dual{_JetTagInner}(v, 1.0),
        ForwardDiff.Dual{_JetTagInner}(1.0, 0.0))
    y = g(z)
    inner = ForwardDiff.value(y)
    return ForwardDiff.value(inner), ForwardDiff.partials(inner, 1),
        ForwardDiff.partials(ForwardDiff.partials(y, 1), 1)
end
@inline _jet_unary(g::F, x::CTSEMTraceJet) where {F} =
    _jet_chain(x, _jet_derivatives(g, x.v)...)

for g in (:sin, :cos, :tan, :sinh, :cosh, :tanh, :atan, :asin, :acos, :cbrt,
    :log2, :log10, :exp2, :exp10)
    @eval Base.$g(x::CTSEMTraceJet) = _jet_unary($g, x)
end
for g in (:erf, :erfc, :erfcx, :logerfc, :loggamma, :gamma, :digamma, :trigamma)
    @eval SpecialFunctions.$g(x::CTSEMTraceJet) = _jet_unary(SpecialFunctions.$g, x)
end

function Base.show(io::IO, x::CTSEMTraceJet{W}) where {W}
    print(io, "CTSEMTraceJet{", W, "}(", x.v, ", ", x.d, ", ", x.s, ")")
end

"""
    _jet_seed(v, dirs, t, W)

The jet for coordinate `t` of `v + sum_j t_j dirs[j]`: value `v[t]` and first
parts `dirs[j][t]`, the second part zero. `dirs` holds the sweep's directions,
at most `W`; the rest of the width carries zeros.
"""
@inline _jet_seed(v::Float64, dirs, t::Int, ::Val{W}) where {W} =
    CTSEMTraceJet{W}(v, ntuple(j -> j <= length(dirs) ? @inbounds(dirs[j][t]) : 0.0,
        Val(W)), 0.0)

# The engine's dual-aware helpers, for jets: every component read where a
# dual's partials are, and the value where a dual's primal is.
_finite_deep(x::CTSEMTraceJet) = isfinite(x.v) && isfinite(x.s) && all(isfinite, x.d)
_zero_deep(x::CTSEMTraceJet) = iszero(x)
@inline _custom_abs(x::CTSEMTraceJet) = abs(x.v)
_ctsem_nan(::Type{CTSEMTraceJet{W}}) where {W} =
    CTSEMTraceJet{W}(NaN, ntuple(_ -> NaN, Val(W)), NaN)
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:CTSEMTraceJet},
    B::AbstractMatrix{<:CTSEMTraceJet}) = A == B
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:CTSEMTraceJet},
    B::AbstractMatrix{<:CTSEMTraceJet}, dim) = _matrix_equal_noalloc(A, B, dim, dim)
@inline _primal(x::CTSEMTraceJet) = x.v
_carries_derivatives(::Type{<:CTSEMTraceJet}) = true
_primal_type(::Type{<:CTSEMTraceJet}) = Float64
@inline _laplace_deepvalue(x::CTSEMTraceJet) = x.v
@inline _laplace_finite(x::CTSEMTraceJet) = _finite_deep(x)
_hybrid_same(a::CTSEMTraceJet{W}, b::CTSEMTraceJet{W}) where {W} = a == b
@inline _hybrid_strip(x::CTSEMTraceJet) = x.v
