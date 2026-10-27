# ---- State-dependent transforms as data ---------------------------------------
#
# A state-dependent cell -- a random effect the augmented route carries as a
# state (`10*state[3]`), a nonlinear DRIFT, a TD predictor in a LAMBDA -- used
# to be `eval`ed into its own closure, and the closures' types were type
# parameters of `EKFParameters`: every model with a different set of such
# cells compiled the filter, its reverse pass and every dual layer again. A
# `CTSEMExpr` is the same expression held as data and evaluated by one compiled
# interpreter, generic in the number type, so every model's transforms are one
# type: `Vector{CTSEMExpr}`.
#
# The vocabulary is what ctsem writes and what users write in a cell:
# constants, `state[k]`, a model matrix's cell, `tdpreds[j]`, `tipreds[j]`,
# the arithmetic operators and the scalar functions below. An expression using
# anything else returns `nothing` from `_ctsem_expr` and keeps its closure.

const _XC = 0x01      # constant
const _XS = 0x02      # state[a]
const _XP = 0x03      # parameter vector entry a (a model matrix cell)
const _XTD = 0x04     # tdpreds[a]
const _XTI = 0x05     # tipreds[a]
const _XNEG = 0x10
const _XEXP = 0x11
const _XLOG = 0x12
const _XSQRT = 0x13
const _XABS = 0x14
const _XLOG1PEXP = 0x15
const _XINVLOGIT = 0x16
const _XLOGIT = 0x17
const _XTANH = 0x18
const _XLOG1P = 0x19
const _XEXPM1 = 0x1a
const _XSIN = 0x1b
const _XCOS = 0x1c
const _XINV = 0x1d
const _XSIGN = 0x1e   # the derivative ctJacobian writes for abs, min and max
const _XPOWI = 0x2f   # a ^ b, b an integer literal: b stored in `b`
const _XADD = 0x30
const _XSUB = 0x31
const _XMUL = 0x32
const _XDIV = 0x33
const _XPOW = 0x34
const _XMIN = 0x35
const _XMAX = 0x36

# What a node's value depends on: nothing (a constant, folded at build), the
# row's data (TD or TI predictors), or the model (a state or a parameter).
const _XK_CONST = 0x00
const _XK_DATA = 0x01
const _XK_MODEL = 0x02

const _CTSEM_EXPR_UNARY = Dict{Symbol,UInt8}(:exp => _XEXP, :log => _XLOG,
    :sqrt => _XSQRT, :abs => _XABS, :log1p_exp => _XLOG1PEXP,
    :inv_logit => _XINVLOGIT, :logit => _XLOGIT, :tanh => _XTANH,
    :log1p => _XLOG1P, :expm1 => _XEXPM1, :sin => _XSIN, :cos => _XCOS,
    :inv => _XINV, :sign => _XSIGN)
const _CTSEM_EXPR_BINARY = Dict{Symbol,UInt8}(:+ => _XADD, :- => _XSUB,
    :* => _XMUL, :/ => _XDIV, :^ => _XPOW, :min => _XMIN, :max => _XMAX)

"""
    CTSEMExpr

A state-dependent transform held as a tree of operations: node `i` has
operation `op[i]`, operands `a[i]` and `b[i]` (child nodes, or a leaf's index),
`val[i]` for a constant, and `kind[i]`, what its value depends on. Called as
`e(ctx)` on a `CTSEMRowContext`, like the closure it replaces.

A subtree reading the model is evaluated in the promoted element type of the
context's state and parameters, so a seeded dual context differentiates it; a
constant stays a `Float64` and a subtree reading only predictors stays in their
type, and either meets the model's values through the same mixed arithmetic
the closure would have done. That keeps the result identical to the closure's,
derivatives included: `p / s` and `p * inv(s)` differ in the last bit, and
promoting `s` to a dual first is the second.
"""
struct CTSEMExpr
    op::Vector{UInt8}
    a::Vector{Int}
    b::Vector{Int}
    val::Vector{Float64}
    kind::Vector{UInt8}
    root::Int
    # The states it reads, in order: the only ones the reverse pass needs its
    # derivative against (`_group_reverse_walk!`).
    states::Vector{Int}
    # A shape evaluated directly rather than walked: `_XFORM_SCALED_STATE`,
    # `fc * state[fa]` -- a random effect the augmented route carries as a
    # state, which is nearly every state-dependent cell a model has -- or
    # `_XFORM_TREE`.
    form::UInt8
    fa::Int
    fc::Float64
    source::String
end

const _XFORM_TREE = 0x00
const _XFORM_SCALED_STATE = 0x01

_ctsem_expr_states(op, a) = sort!(unique(Int[a[i] for i in eachindex(op) if op[i] === _XS]))

# The direct form of a built tree, as `(form, fa, fc)`.
function _ctsem_expr_form(op, a, b, val, root)
    if op[root] === _XMUL
        l = a[root]
        r = b[root]
        op[l] === _XC && op[r] === _XS && return (_XFORM_SCALED_STATE, a[r], val[l])
        op[r] === _XC && op[l] === _XS && return (_XFORM_SCALED_STATE, a[l], val[r])
    end
    return (_XFORM_TREE, 0, 0.0)
end

function _ctsem_expr_make(bld, root, states, source)
    form, fa, fc = _ctsem_expr_form(bld.op, bld.a, bld.b, bld.val, root)
    return CTSEMExpr(bld.op, bld.a, bld.b, bld.val, bld.kind, root, states,
        form, fa, fc, String(source))
end

Base.show(io::IO, e::CTSEMExpr) = print(io, "CTSEMExpr(", repr(e.source), ")")

@inline _ctsem_expr_type(ctx) = promote_type(eltype(ctx.state), eltype(ctx.pars))
@inline _ctsem_expr_datatype(ctx) =
    promote_type(Float64, eltype(ctx.tdpreds), eltype(ctx.tipreds))

(e::CTSEMExpr)(ctx) = _ctsem_expr_eval(e, ctx, nothing)

"""
    _ctsem_expr_eval(e, ctx, scratch)

Evaluate `e` in `ctx`. `scratch` is the pair of task buffers a long expression
evaluates in (`_ctsem_expr_scratches(ctx)`), or `nothing` to look them up here.
A loop over a group's expressions fetches them once (`_ctsem_expr_any_long`):
the lookup is two `IdDict` reads per element type, and per expression it was
about 5% of a nonlinear Laplace gradient.
"""
function _ctsem_expr_eval(e::CTSEMExpr, ctx, scratch)
    T = _ctsem_expr_type(ctx)
    # `fc * x` with `x` in `T`, as `10 * state[3]` multiplies in the closure.
    e.form === _XFORM_SCALED_STATE && return convert(T, e.fc * convert(T, @inbounds ctx.state[e.fa]))
    r = e.root
    @inbounds k = e.kind[r]
    k === _XK_CONST && return convert(T, @inbounds e.val[r])
    S = _ctsem_expr_datatype(ctx)
    n = length(e.op)
    # Up to four nodes in tuples, which live on the stack; past that in the
    # task's buffers (`_ctsem_expr_scratch`), since a write into a tuple
    # selects across it. Five to eight nodes go to tuples only where nobody
    # has fetched the buffers, whose lookup then costs more than the selects:
    # with them in hand the buffers were faster, up to four times at Dual
    # types (`scripts/stangrad-2026-10/bench_expr.jl`, local).
    n <= 4 && return _ctsem_expr_flat(e, ctx, T, S,
        ntuple(_ -> zero(T), Val(4)), ntuple(_ -> zero(S), Val(4)))
    scratch === nothing && n <= 8 && return _ctsem_expr_flat(e, ctx, T, S,
        ntuple(_ -> zero(T), Val(8)), ntuple(_ -> zero(S), Val(8)))
    scT, scS = scratch === nothing ? _ctsem_expr_scratches(ctx) : scratch
    return _ctsem_expr_flat(e, ctx, T, S, _ctsem_expr_buffer(scT.a, n),
        _ctsem_expr_buffer(scS.b, n))
end

"""The task buffers `ctx`'s long expressions evaluate in, one per element type."""
@inline _ctsem_expr_scratches(ctx) =
    (_ctsem_expr_scratch(_ctsem_expr_type(ctx)), _ctsem_expr_scratch(_ctsem_expr_datatype(ctx)))

"""Whether any of a group's expressions is long enough to use the task buffers."""
@inline function _ctsem_expr_any_long(transforms::Vector{CTSEMExpr})
    @inbounds for e in transforms
        length(e.op) > 4 && return true
    end
    return false
end

# Nodes are stored children first (`_ctsem_expr_node!`), so one pass in index
# order evaluates the tree, each node's value held in a slot of `vt` (model
# values, in `T`) or `vs` (predictor values, in `S`): tuples or vectors. No call
# per node: a recursive walk cost a call and its frame for every node, about
# 9 ns, and made a long expression -- the derivative R writes for a long
# nonlinear cell can run to a hundred nodes -- several times its closure.

@inline _ctsem_expr_put(t::Tuple, v, i::Int) = Base.setindex(t, v, i)
@inline _ctsem_expr_put(a::Vector, v, i::Int) = (@inbounds a[i] = v; a)

"""
A task's buffers for evaluating and differentiating a long expression, one
pair per element type, in the task-local store the engine's other per-task
caches use (`_ctsem_scratch`). An expression's evaluation calls no other, and
its reverse sweep evaluates none, so one pair per type is enough; with `T`
and `S` both Float64 the forward pass takes `a` for one and `b` for the other.
"""
struct _CTSEMExprScratch{T}
    a::Vector{T}
    b::Vector{T}
end
struct _CTSEMExprScratchKey{T} end

@inline _ctsem_expr_scratch(::Type{T}) where {T} =
    get!(() -> _CTSEMExprScratch{T}(T[], T[]), _ctsem_scratch(),
        _CTSEMExprScratchKey{T}())::_CTSEMExprScratch{T}

@inline _ctsem_expr_buffer(v::Vector, n::Int) = (length(v) < n && resize!(v, n); v)

function _ctsem_expr_flat(e::CTSEMExpr, ctx, ::Type{T}, ::Type{S}, vt, vs) where {T,S}
    n = length(e.op)
    @inbounds for i in 1:n
        op = e.op[i]
        op === _XC && continue
        if e.kind[i] === _XK_MODEL
            v = if op === _XP
                convert(T, ctx.pars[e.a[i]])
            elseif op === _XS
                convert(T, ctx.state[e.a[i]])
            elseif op >= _XADD
                # Each operand in the type it was computed in, so a constant or
                # a predictor meets the model's values as in the closure.
                a = e.a[i]
                b = e.b[i]
                ka = e.kind[a]
                kb = e.kind[b]
                if ka === _XK_MODEL
                    x = vt[a]
                    kb === _XK_MODEL ? convert(T, _ctsem_expr_binary(op, x, vt[b])) :
                        kb === _XK_DATA ? convert(T, _ctsem_expr_binary(op, x, vs[b])) :
                        convert(T, _ctsem_expr_binary(op, x, e.val[b]))
                elseif ka === _XK_DATA
                    convert(T, _ctsem_expr_binary(op, vs[a], vt[b]))
                else
                    convert(T, _ctsem_expr_binary(op, e.val[a], vt[b]))
                end
            elseif op === _XPOWI
                convert(T, _ctsem_expr_powi(vt[e.a[i]], e.b[i]))
            else
                convert(T, _ctsem_expr_unary(op, vt[e.a[i]]))
            end
            vt = _ctsem_expr_put(vt, v, i)
        else
            w = if op === _XTD
                convert(S, ctx.tdpreds[e.a[i]])
            elseif op === _XTI
                convert(S, ctx.tipreds[e.a[i]])
            elseif op >= _XADD
                a = e.a[i]
                b = e.b[i]
                if e.kind[a] === _XK_DATA
                    e.kind[b] === _XK_DATA ? convert(S, _ctsem_expr_binary(op, vs[a], vs[b])) :
                        convert(S, _ctsem_expr_binary(op, vs[a], e.val[b]))
                else
                    convert(S, _ctsem_expr_binary(op, e.val[a], vs[b]))
                end
            elseif op === _XPOWI
                convert(S, _ctsem_expr_powi(vs[e.a[i]], e.b[i]))
            else
                convert(S, _ctsem_expr_unary(op, vs[e.a[i]]))
            end
            vs = _ctsem_expr_put(vs, w, i)
        end
    end
    r = e.root
    return @inbounds e.kind[r] === _XK_MODEL ? vt[r] : convert(T, vs[r])
end

# `x ^ n` for an integer literal `n`, as Julia evaluates it in source:
# `literal_pow`, so the value matches the closure this replaces.
@inline function _ctsem_expr_powi(x, n::Int)
    n == 2 && return Base.literal_pow(^, x, Val(2))
    n == 3 && return Base.literal_pow(^, x, Val(3))
    n == 1 && return Base.literal_pow(^, x, Val(1))
    n == 0 && return Base.literal_pow(^, x, Val(0))
    n == -1 && return Base.literal_pow(^, x, Val(-1))
    n == -2 && return Base.literal_pow(^, x, Val(-2))
    n == 4 && return Base.literal_pow(^, x, Val(4))
    return x^n
end

@inline function _ctsem_expr_binary(op::UInt8, x, y)
    op === _XMUL && return x * y
    op === _XADD && return x + y
    op === _XSUB && return x - y
    op === _XDIV && return x / y
    op === _XPOW && return x^y
    op === _XMIN && return min(x, y)
    return max(x, y)
end

@inline function _ctsem_expr_unary(op::UInt8, x)
    op === _XNEG && return -x
    op === _XEXP && return exp(x)
    op === _XLOG && return log(x)
    op === _XSQRT && return sqrt(x)
    op === _XABS && return abs(x)
    op === _XLOG1PEXP && return log1p_exp(x)
    op === _XINVLOGIT && return inv_logit(x)
    op === _XLOGIT && return logit(x)
    op === _XTANH && return tanh(x)
    op === _XLOG1P && return log1p(x)
    op === _XEXPM1 && return expm1(x)
    op === _XSIN && return sin(x)
    op === _XCOS && return cos(x)
    # In `x`'s own type with zero partials, as Base's sign gives a Dual away
    # from zero; `zero(x) - one(x)`, since `-one(x)`'s partials are -0.0.
    op === _XSIGN && return x > 0 ? one(x) : x < 0 ? zero(x) - one(x) : zero(x)
    return inv(x)
end

# ---- Building one from source -------------------------------------------------

mutable struct _CTSEMExprBuilder
    op::Vector{UInt8}
    a::Vector{Int}
    b::Vector{Int}
    val::Vector{Float64}
    kind::Vector{UInt8}
    layout::Union{Nothing,CTSEMLayout}
end

function _ctsem_expr_node!(bld::_CTSEMExprBuilder, op::UInt8, a::Int=0, b::Int=0,
    val::Float64=0.0)
    kind = op === _XC ? _XK_CONST :
        (op === _XS || op === _XP) ? _XK_MODEL :
        (op === _XTD || op === _XTI) ? _XK_DATA :
        op >= _XADD ? max(bld.kind[a], bld.kind[b]) : bld.kind[a]
    push!(bld.op, op); push!(bld.a, a); push!(bld.b, b); push!(bld.val, val)
    push!(bld.kind, kind)
    return _ctsem_expr_fold!(bld, length(bld.op))
end

# An integer literal, as Julia's parser gives it.
_ctsem_expr_intlit(x) = x isa Integer && !(x isa Bool) ? Int(x) : nothing

"""
The flat parameter index of `ctx.pars.<name>[i, j]` (or `[i]`), or `nothing`
when the layout has no such cell. Without a layout -- `ctsem_transforms_cached`
only asks whether an expression is in the vocabulary -- any cell is index 1.
"""
function _ctsem_expr_cell(bld::_CTSEMExprBuilder, name::Symbol, idx)
    all(i -> _ctsem_expr_intlit(i) !== nothing, idx) || return nothing
    (length(idx) == 1 || length(idx) == 2) || return nothing
    l = bld.layout
    l === nothing && return 1
    m = _layout_index(l, name)
    m == 0 && return nothing
    nr = l.nrows[m]
    nc = l.ncols[m]
    if length(idx) == 1
        k = Int(idx[1])
        1 <= k <= nr * nc || return nothing
        return l.first[m] + k - 1
    end
    r = Int(idx[1]); c = Int(idx[2])
    (1 <= r <= nr && 1 <= c <= nc) || return nothing
    return l.first[m] + (c - 1) * nr + (r - 1)
end

# `ctx.<field>` as the parser gives it.
_ctsem_expr_ctxfield(x, f::Symbol) = x isa Expr && x.head === :. &&
    x.args[1] === :ctx && x.args[2] == QuoteNode(f)

function _ctsem_expr_build!(bld::_CTSEMExprBuilder, x)
    if x isa Real && !(x isa Bool)
        return _ctsem_expr_node!(bld, _XC, 0, 0, Float64(x))
    end
    x isa Expr || return nothing
    if x.head === :ref
        target = x.args[1]
        idx = x.args[2:end]
        if _ctsem_expr_ctxfield(target, :param)
            # A regular transform's raw parameter (`_ctsem_regular_expr`).
            length(idx) == 1 || return nothing
            k = _ctsem_expr_intlit(idx[1])
            (k === nothing || k < 1) && return nothing
            return _ctsem_expr_node!(bld, _XP, k)
        elseif _ctsem_expr_ctxfield(target, :state)
            k = _ctsem_expr_intlit(idx[1])
            (k === nothing || k < 1) && return nothing
            length(idx) == 1 || (length(idx) == 2 && idx[2] == 1) || return nothing
            return _ctsem_expr_node!(bld, _XS, k)
        elseif _ctsem_expr_ctxfield(target, :tdpreds) || _ctsem_expr_ctxfield(target, :tipreds)
            length(idx) == 1 || return nothing
            k = _ctsem_expr_intlit(idx[1])
            (k === nothing || k < 1) && return nothing
            return _ctsem_expr_node!(bld,
                _ctsem_expr_ctxfield(target, :tdpreds) ? _XTD : _XTI, k)
        elseif target isa Expr && target.head === :. && _ctsem_expr_ctxfield(target.args[1], :pars) &&
               target.args[2] isa QuoteNode
            flat = _ctsem_expr_cell(bld, target.args[2].value, idx)
            flat === nothing && return nothing
            return _ctsem_expr_node!(bld, _XP, flat)
        end
        return nothing
    end
    x.head === :call || return nothing
    f = x.args[1]
    f isa Symbol || return nothing
    args = x.args[2:end]
    if f === :^ && length(args) == 2 && _ctsem_expr_intlit(args[2]) !== nothing
        base = _ctsem_expr_build!(bld, args[1])
        base === nothing && return nothing
        return _ctsem_expr_node!(bld, _XPOWI, base, Int(args[2]))
    end
    if length(args) == 1
        child = _ctsem_expr_build!(bld, args[1])
        child === nothing && return nothing
        f === :+ && return child
        op = f === :- ? _XNEG : get(_CTSEM_EXPR_UNARY, f, 0x00)
        op === 0x00 && return nothing
        return _ctsem_expr_node!(bld, op, child)
    end
    op = get(_CTSEM_EXPR_BINARY, f, 0x00)
    op === 0x00 && return nothing
    # `+` and `*` take any number of operands and fold them from the left, as
    # Base does for numbers; the others are binary.
    (length(args) == 2 || ((f === :+ || f === :*) && length(args) > 2)) || return nothing
    acc = _ctsem_expr_build!(bld, args[1])
    acc === nothing && return nothing
    for k in 2:length(args)
        rhs = _ctsem_expr_build!(bld, args[k])
        rhs === nothing && return nothing
        acc = _ctsem_expr_node!(bld, op, acc, rhs)
    end
    return acc
end

# A node whose operands are all constants is replaced by its value, computed
# in Float64 as the closure would have.
function _ctsem_expr_fold!(bld::_CTSEMExprBuilder, i::Int)
    op = bld.op[i]
    (op >= _XNEG && bld.kind[i] === _XK_CONST) || return i
    v = if op >= _XADD
        Float64(_ctsem_expr_binary(op, bld.val[bld.a[i]], bld.val[bld.b[i]]))
    elseif op === _XPOWI
        Float64(_ctsem_expr_powi(bld.val[bld.a[i]], bld.b[i]))
    else
        Float64(_ctsem_expr_unary(op, bld.val[bld.a[i]]))
    end
    bld.op[i] = _XC; bld.a[i] = 0; bld.b[i] = 0; bld.val[i] = v
    return i
end

"""
    _ctsem_expr(source, layout)

The `CTSEMExpr` for a state-dependent transform written as ctsem writes them
(`state[3]`, `PARS[1,2]`, `ctx.tdpreds[1]`), or `nothing` when it uses
something outside the vocabulary, which then keeps its compiled closure.
`layout === nothing` asks only whether the expression is in the vocabulary.
"""
function _ctsem_expr(source::AbstractString, layout::Union{Nothing,CTSEMLayout})
    text = generate_complex_transform_string(String(source))
    parsed = try
        Meta.parse(text)
    catch
        return nothing
    end
    (parsed isa Expr && parsed.head === :-> && parsed.args[1] === :ctx) || return nothing
    body = parsed.args[2]
    if body isa Expr && body.head === :block
        stmts = filter(s -> !(s isa LineNumberNode), body.args)
        length(stmts) == 1 || return nothing
        body = stmts[1]
    end
    bld = _CTSEMExprBuilder(UInt8[], Int[], Int[], Float64[], UInt8[], layout)
    root = try
        _ctsem_expr_build!(bld, body)
    catch
        nothing
    end
    root === nothing && return nothing
    return _ctsem_expr_make(bld, root, _ctsem_expr_states(bld.op, bld.a), source)
end

# ---- Regular transforms --------------------------------------------------------
#
# A regular transform maps the subject's raw parameter vector to one cell
# (`10 * log1p_exp(2 * param[3]) + 1e-10`). The templates ctsem writes by
# default stay compiled closures, one group each (`_group_regular_transforms`);
# any other -- a transform a user wrote, a scale that is not the default, a
# T0MEANS cell composed from several parameters -- is held as a `CTSEMExpr`
# reading `param[k]`, in one group every model carries, rather than as a
# closure type of its own that compiled the model again.

const _CTSEM_EXPR_NOSTATE = Float64[]

# Called on the raw parameter vector, as the closure is.
(e::CTSEMExpr)(param::AbstractVector) = e((state = _CTSEM_EXPR_NOSTATE, pars = param,
    tdpreds = _CTSEM_EXPR_NOSTATE, tipreds = _CTSEM_EXPR_NOSTATE))

"""
    _ctsem_regular_expr(text)

The `CTSEMExpr` for a regular transform written as ctsem writes them, with
`param[k]` reading raw parameter `k`; `nothing` outside the vocabulary.
"""
function _ctsem_regular_expr(text::AbstractString)
    lifted = replace(String(text), r"\bparam\s*\[" => "ctx.param[")
    occursin(r"\b(state|tdpreds|tipreds)\b", lifted) && return nothing
    parsed = try
        Meta.parse("ctx -> " * lifted)
    catch
        return nothing
    end
    (parsed isa Expr && parsed.head === :->) || return nothing
    body = parsed.args[2]
    if body isa Expr && body.head === :block
        stmts = filter(s -> !(s isa LineNumberNode), body.args)
        length(stmts) == 1 || return nothing
        body = stmts[1]
    end
    bld = _CTSEMExprBuilder(UInt8[], Int[], Int[], Float64[], UInt8[], nothing)
    root = try
        _ctsem_expr_build!(bld, body)
    catch
        nothing
    end
    root === nothing && return nothing
    # A transform reading no parameter is a rendering fault the support probe
    # reports; leave it to the closure, which reaches that check the same way.
    any(==(_XP), bld.op) || return nothing
    return _ctsem_expr_make(bld, root, Int[], text)
end

const _CTSEM_DEFAULT_TEMPLATE_KEYS = Ref{Union{Nothing,Set{String}}}(nothing)

"""Whether a regular transform's template is one ctsem writes by default."""
function _ctsem_default_template(text::AbstractString)
    keys = _CTSEM_DEFAULT_TEMPLATE_KEYS[]
    if keys === nothing
        keys = Set(ctsem_default_template_keys())
        _CTSEM_DEFAULT_TEMPLATE_KEYS[] = keys
    end
    return _regular_transform_key(text) in keys
end

"""
The container `EKFParameters` holds a group's state-dependent transforms in:
`Vector{CTSEMExpr}` when every one is an expression -- one type for every
model, the empty group included -- and otherwise a tuple of closures, which
compiles per model as it always has.
"""
function _ctsem_transform_container(transforms)
    all(t -> t isa CTSEMExpr, transforms) && return CTSEMExpr[t for t in transforms]
    return Tuple(transforms)
end
