using LinearAlgebra
using ForwardDiff

################################################################################
# Helper functions 
################################################################################
"""
    _custom_abs(x)

Return `abs(x)` for ordinary numeric values.

This method is paired with a `ForwardDiff.Dual` specialization so pivoting and
approximate comparisons can work on primal values during automatic
differentiation.
"""
@inline _custom_abs(x) = abs(x)

"""
    _finite_deep(x)

Whether `x` is finite *including its derivatives*.

`isfinite` on a `ForwardDiff.Dual` tests only the value, so a NaN in the
partials passes every validity guard in this package and surfaces much later as
a rejected trial point with a finite objective and an unusable gradient. This
walks the whole dual tree instead.
"""
_finite_deep(x::Real) = isfinite(x)
_finite_deep(x::ForwardDiff.Dual) = _finite_deep(ForwardDiff.value(x)) &&
    all(_finite_deep, ForwardDiff.partials(x))
_finite_deep(x::AbstractArray) = all(_finite_deep, x)

"""
    _ctsem_nan(T)

NaN of type `T` in the value and in every partial, at every level of a dual
number: the sentinel for an invalid point. `T(NaN)` sets the value alone, and a
Jacobian reads only the partials, so a gradient filled with it differentiated to
an exact zero -- `ctsem_hessian` returned a zero matrix at a point whose dual
pass had failed, which every check downstream took for a Hessian.
"""
_ctsem_nan(::Type{T}) where {T} = T(NaN)
_ctsem_nan(::Type{ForwardDiff.Dual{G,V,N}}) where {G,V,N} = ForwardDiff.Dual{G}(
    _ctsem_nan(V), ForwardDiff.Partials(ntuple(_ -> _ctsem_nan(V), Val(N))))

"""
    _zero_deep(x)

Whether `x` is exactly zero *including its derivatives*: the value and every
partial at every nesting level. A dual whose value is zero still moves when its
partials do not vanish, so a structural shortcut for a zero input asks this,
not `iszero` on the value.
"""
_zero_deep(x::Real) = iszero(x)
_zero_deep(x::ForwardDiff.Dual) = _zero_deep(ForwardDiff.value(x)) &&
    all(_zero_deep, ForwardDiff.partials(x))

"""
    _custom_abs(x::ForwardDiff.Dual)

Return the absolute value of the primal part of a dual number.
"""
@inline _custom_abs(x::ForwardDiff.Dual) = abs(ForwardDiff.value(x))


"""
    _isapprox_default_rtol(T)

Return the default relative tolerance used by `_isapprox_matrix_noalloc`.

For non-floating scalar types this falls back to `sqrt(eps(Float64))`.
"""
@inline _isapprox_default_rtol(::Type) = sqrt(eps(Float64))

"""
    _isapprox_default_rtol(::Type{T}) where {T<:AbstractFloat}

Return `sqrt(eps(T))` for floating-point matrix element types.
"""
@inline _isapprox_default_rtol(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(T))

"""
    _isapprox_default_rtol(::Type{<:ForwardDiff.Dual})

Return a relative tolerance based on the dual number's floating primal type.
"""
@inline _isapprox_default_rtol(::Type{<:ForwardDiff.Dual{Tag,V,N}}) where {Tag,V<:AbstractFloat,N} = sqrt(eps(V))

"""
    add_diag!(mat, val::Number)

Add the scalar `val` to each diagonal entry of `mat` in place.
"""
function add_diag!(mat, val::Number)
    return add_diag!(mat, val, min(size(mat, 1), size(mat, 2)))
end

function add_diag!(mat, val::Number, __d::Int)
    d = __d
    @boundscheck begin
        size(mat, 1) >= d || throw(DimensionMismatch("matrix row count must be at least d"))
        size(mat, 2) >= d || throw(DimensionMismatch("matrix column count must be at least d"))
    end
    @inbounds for idx in 1:d
        mat[idx, idx] += val
    end
    return mat
end

"""
    add_diag!(mat, vals::AbstractVector)

Add `vals[i]` to the `i`th diagonal entry of `mat` in place.
"""
function add_diag!(mat, vals::AbstractVector)
    return add_diag!(mat, vals, length(vals))
end

function add_diag!(mat, vals::AbstractVector, __d::Int)
    d = __d
    @boundscheck begin
        length(vals) == d || throw(DimensionMismatch("vals length must match d"))
        size(mat, 1) >= d || throw(DimensionMismatch("matrix row count must be at least d"))
        size(mat, 2) >= d || throw(DimensionMismatch("matrix column count must be at least d"))
    end
    @inbounds for idx in 1:d
        mat[idx, idx] += vals[idx]
    end
    return mat
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` for BLAS-compatible strided matrices.
"""
function _mul_right_transpose!(C::StridedMatrix{T}, A::StridedMatrix{T}, B::StridedMatrix{T}) where {T<:LinearAlgebra.BlasFloat}
    return _mul_right_transpose!(C, A, B, size(C, 1), size(C, 2), size(A, 2))
end

function _mul_right_transpose!(
    C::StridedMatrix{T},
    A::StridedMatrix{T},
    B::StridedMatrix{T},
    __rows::Int,
    __cols::Int,
    __inner::Int,
) where {T<:LinearAlgebra.BlasFloat}
    rows = __rows
    cols = __cols
    inner = __inner
    @boundscheck begin
        size(C, 1) == rows && size(A, 1) == rows || throw(DimensionMismatch("A row count must match C row count"))
        size(C, 2) == cols && size(B, 1) == cols || throw(DimensionMismatch("B row count must match C column count"))
        size(A, 2) == inner && size(B, 2) == inner || throw(DimensionMismatch("A and B must have the same column count"))
    end
    # Through `_ctsem_mulNT!`: the engine's kernels, BLAS only above
    # `_CTSEM_SMALL_PRODUCT`, which is never by default.
    return _ctsem_mulNT!(C, A, B)
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` with explicit loops for generic numeric matrices.
"""
function _mul_right_transpose!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T<:Number}
    return _mul_right_transpose!(C, A, B, size(C, 1), size(C, 2), size(A, 2))
end

function _mul_right_transpose!(
    C::AbstractMatrix{T},
    A::AbstractMatrix{T},
    B::AbstractMatrix{T},
    __rows::Int,
    __cols::Int,
    __inner::Int,
) where {T<:Number}
    rows = __rows
    cols = __cols
    inner = __inner
    @boundscheck begin
        size(C, 1) == rows && size(A, 1) == rows || throw(DimensionMismatch("A row count must match C row count"))
        size(C, 2) == cols && size(B, 1) == cols || throw(DimensionMismatch("B row count must match C column count"))
        size(A, 2) == inner && size(B, 2) == inner || throw(DimensionMismatch("A and B must have the same column count"))
    end
    @inbounds for j in 1:cols
        for i in 1:rows
            acc = zero(T)
            for k in 1:inner
                acc += A[i, k] * B[j, k]
            end
            C[i, j] = acc
        end
    end
    return C
end

"""
    _isapprox_matrix_noalloc(A, B; atol=0.0, rtol=...)

Return whether two matrices are approximately equal without allocating.

Dual-number entries are compared using their primal values through
`_custom_abs`.
"""
function _isapprox_matrix_noalloc(
    A::AbstractMatrix,
    B::AbstractMatrix;
    atol::Real=0.0,
    rtol::Real=_isapprox_default_rtol(promote_type(eltype(A), eltype(B))),
)
    size(A) == size(B) || return false
    return _isapprox_matrix_noalloc(A, B, size(A, 1), size(A, 2); atol=atol, rtol=rtol)
end

function _isapprox_matrix_noalloc(
    A::AbstractMatrix,
    B::AbstractMatrix,
    __rows::Int,
    __cols::Int;
    atol::Real=0.0,
    rtol::Real=_isapprox_default_rtol(promote_type(eltype(A), eltype(B))),
)
    rows = __rows
    cols = __cols
    @boundscheck begin
        size(A, 1) == rows && size(B, 1) == rows || return false
        size(A, 2) == cols && size(B, 2) == cols || return false
    end
    @inbounds for j in 1:cols, i in 1:rows
        aij = A[i, j]
        bij = B[i, j]
        diff = _custom_abs(aij - bij)
        scale = max(_custom_abs(aij), _custom_abs(bij))
        if diff > atol + rtol * scale
            return false
        end
    end
    return true
end

function _matrix_equal_noalloc(A::AbstractMatrix, B::AbstractMatrix, __rows::Int, __cols::Int)
    rows = __rows
    cols = __cols
    @boundscheck begin
        size(A, 1) == rows && size(B, 1) == rows || return false
        size(A, 2) == cols && size(B, 2) == cols || return false
    end
    @inbounds for j in 1:cols, i in 1:rows
        A[i, j] == B[i, j] || return false
    end
    return true
end

"""
    _can_reuse_same_exponential(A, B)

Return whether the exponential computed for `A` can be reused for `B`.

For dual-number matrices, exact equality is required so derivative information
is preserved.
"""
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:ForwardDiff.Dual}, B::AbstractMatrix{<:ForwardDiff.Dual}) = A == B
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:ForwardDiff.Dual}, B::AbstractMatrix{<:ForwardDiff.Dual}, dim) =
    _matrix_equal_noalloc(A, B, dim, dim)

"""
    _can_reuse_same_exponential(A, B)

Return whether two ordinary matrices are close enough to share an exponential.
"""
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix) = _isapprox_matrix_noalloc(A, B)
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix, dim) =
    _isapprox_matrix_noalloc(A, B, dim, dim)


"""
    _ctsem_must_propagate(err)

Whether an exception caught while evaluating a trial point must be rethrown
rather than converted into an invalid-point sentinel.

The catches around trial evaluations exist for facts about a parameter vector:
a curvature that will not factorize, a matrix exponential whose scaling step
takes `ceil(Int, NaN)`. The engine's contract is that such a point reports NaN
or is rejected. These are not that. `MethodError`, `UndefVarError`,
`UndefKeywordError`, `BoundsError` and `TypeError` mean the code is wrong --
when the worker pool dropped the `slot` argument, two callers kept passing it,
and their catches turned the `MethodError` into a NaN quadrature gap and an
identity sampler metric, with nothing reported. A `ForwardDiff.DualMismatchError`
is the same kind of fact about the code: two dual tags met in an order
ForwardDiff cannot resolve (see `CTSEMNestedTag` below), which no parameter
vector causes or avoids, and scoring it as an invalid point ends a fit early
with nothing said. An `InterruptException` is the user stopping the run. A
`TaskFailedException` or `CompositeException` from a spawned region is
unwrapped, since the error that matters is inside it.
"""
function _ctsem_must_propagate(err)
    err isa TaskFailedException && return _ctsem_must_propagate(err.task.result)
    err isa CompositeException && return any(_ctsem_must_propagate, err.exceptions)
    return err isa MethodError || err isa UndefVarError ||
        err isa UndefKeywordError || err isa BoundsError || err isa TypeError ||
        err isa ForwardDiff.DualMismatchError || err isa InterruptException
end

"""
    _ctsem_try_evaluate(f)

`f()`, or `nothing` where it throws what a point the model cannot evaluate
throws; code that is wrong, and an interrupt, still propagate
(`_ctsem_must_propagate`). Every trial and probe point is read through it.
"""
function _ctsem_try_evaluate(f)
    try
        return f()
    catch err
        _ctsem_must_propagate(err) && rethrow()
        return nothing
    end
end

"""Whether a gradient may steer the optimiser: every entry finite, and under
`limit` in size."""
_ctsem_gradient_usable(g, limit::Real) =
    all(isfinite, g) && all(abs(v) < limit for v in g)


"""
Differentiation nested inside code that may itself be differentiated.

ForwardDiff decides which of two tags is the outer one by `tagcount`: a
`@generated` function that numbers a tag type from a counter when that type is
first compiled. The numbers therefore follow compilation order -- across the
package image and the session -- and not nesting. Nothing guarantees that an
inner tag outranks the outer one, and when it does not, the two meet in a
conversion and ForwardDiff throws `DualMismatchError`.

Measured on a censored model: `ctsem_hessian` differentiates the adjoint
gradient, whose measurement update takes its own jacobian of the quadrature
(`_binary_moment_jacobian`), and it threw exactly that when the censored
row's standard deviation -- an outer dual, since MANIFESTVAR is differentiated
-- was converted into the inner dual type. The Hessian was unavailable for
every censored fit, so certification and standard errors had none. The
forward-over-forward Hessian, which nests nothing, was finite and agreed with
finite differences. Why the numbers came out inverted for that model and not
for binary or ordinal ones was not established; the fix does not depend on it.

An inner differentiation uses this tag. Its number is above any ordinary tag's,
and rises by one for each level of dual nesting in its input, so it is inner to
everything it can be nested in, whatever was compiled when.
"""
struct CTSEMNestedTag end

_ctsem_dual_depth(::Type) = 0
_ctsem_dual_depth(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N} =
    1 + _ctsem_dual_depth(V)

ForwardDiff.tagcount(::Type{ForwardDiff.Tag{CTSEMNestedTag,V}}) where {V} =
    (typemax(UInt) >> 1) + UInt(_ctsem_dual_depth(V))

_ctsem_nested_tag(x::AbstractArray) = ForwardDiff.Tag{CTSEMNestedTag,eltype(x)}()

"""
    _ctsem_nested_seed(value, i, Val(N))

`value` as a dual under `CTSEMNestedTag`, seeded along direction `i` of `N`:
one input of an inner differentiation that takes all `N` directions in one
chunk, seeded by hand. `_binary_moment_jacobian`, the one such
differentiation, went through `ForwardDiff.jacobian` before, whose config, the
closure's result array and the result matrix were over half of what each of
its calls allocated -- one call per categorical observation per pass. The
seeds are the ones `ForwardDiff.jacobian` builds, so the numbers are the same.
"""
@inline _ctsem_nested_seed(value::T, i::Int, ::Val{N}) where {T,N} =
    ForwardDiff.Dual{ForwardDiff.Tag{CTSEMNestedTag,T},T,N}(value,
        ForwardDiff.Partials{N,T}(ntuple(j -> j == i ? one(T) : zero(T), Val(N))))

"""
    _ctsem_gradient(f, x[, widths]), _ctsem_jacobian(f, x[, widths]),
    _ctsem_dual_hessian(f, x[, chunk; widths])

`ForwardDiff.gradient`, `jacobian` and `hessian` of `f` at `x`, at a dual
width from `widths` (`_ctsem_dual_width`) rather than `length(x)`, behind
`_ctsem_barrier`. `widths` is `_CTSEM_HESSIAN_WIDTHS[]` for an input as long
as the parameters; a caller over a unit's random effects passes
`_CTSEM_CURVATURE_WIDTHS[]`. An input shorter than its width is padded with
zeros `f` never reads (`_CTSEMLeading`, which also gives one `f` one tag type
at every width), and the padding's rows and columns are dropped from the
result: the same numbers. An explicit `chunk` to `_ctsem_dual_hessian` is used
as it is.

A width is a type, so ForwardDiff's own choice, the input length up to 12,
compiled everything these differentiate -- the filter, its reverse pass, the
Laplace unit term with its inner solve -- again for every parameter count. In
the engine's tests, whose ForwardDiff oracles run over fixtures of several
sizes, compiling was nearly all of the time: test_laplace.jl ran in 3.5 s
once compiled, after 1711 s of compiling.

The barrier, because the width is a run-time value: a direct call hands its
worker a config whose type inference cannot know, and inference then walks the
whole dual pipeline below it for an unknown width. None of that code ever runs,
but it was two thirds of a Laplace model's first-fit compile, and the package
image kept it (review/PRECOMPILE-workload-2026-10-12.md). It costs one dynamic
dispatch a call, which the call made anyway.
"""
function _ctsem_gradient(f::F, x::AbstractVector,
        widths::Vector{Int}=_CTSEM_HESSIAN_WIDTHS[]) where {F}
    n = length(x)
    n == 0 && return zeros(eltype(x), 0)
    w = _ctsem_dual_width(n, widths)
    g = _CTSEMLeading(f, n)
    w <= n && return _ctsem_barrier(_ctsem_gradient_at, g, x,
        ForwardDiff.Chunk{w}())::Vector{eltype(x)}
    xp = vcat(x, zeros(eltype(x), w - n))
    full = _ctsem_barrier(_ctsem_gradient_at, g, xp,
        ForwardDiff.Chunk{w}())::Vector{eltype(x)}
    return full[1:n]
end
_ctsem_jacobian(f::F, x::AbstractVector,
    widths::Vector{Int}=_CTSEM_HESSIAN_WIDTHS[]) where {F} =
    _ctsem_width_jacobian(f, x, widths)
function _ctsem_dual_hessian(f::F, x::AbstractVector, chunk=nothing;
        widths::Vector{Int}=_CTSEM_HESSIAN_WIDTHS[]) where {F}
    chunk === nothing || return _ctsem_barrier(_ctsem_hessian_at, f, x,
        chunk)::Matrix{eltype(x)}
    n = length(x)
    n == 0 && return zeros(eltype(x), 0, 0)
    w = _ctsem_dual_width(n, widths)
    g = _CTSEMLeading(f, n)
    w <= n && return _ctsem_barrier(_ctsem_hessian_at, g, x,
        ForwardDiff.Chunk{w}())::Matrix{eltype(x)}
    xp = vcat(x, zeros(eltype(x), w - n))
    full = _ctsem_barrier(_ctsem_hessian_at, g, xp,
        ForwardDiff.Chunk{w}())::Matrix{eltype(x)}
    return full[1:n, 1:n]
end

_ctsem_gradient_at(f::F, x, chunk::ForwardDiff.Chunk) where {F} =
    ForwardDiff.gradient(f, x, ForwardDiff.GradientConfig(f, x, chunk))
_ctsem_hessian_at(f::F, x, chunk::ForwardDiff.Chunk) where {F} =
    ForwardDiff.hessian(f, x, ForwardDiff.HessianConfig(f, x, chunk))

"""
Dual widths for the forward-mode Jacobians whose input length is a property of
the model: the Hessian's (the parameter count) and the Laplace unit
curvature's (the unit's random-effect dimension).

`ForwardDiff.pickchunksize(n)` is `n` itself up to 12, and a dual width is a
type, so every parameter count compiled the whole filter and reverse pass
again: about 30 s on dev2 for one more free parameter on a model the session,
or the package image, had already compiled. Rounding the width up to the next
of a few buckets makes every count in a bucket one type, and past the largest
the width is the largest, in more sweeps. The extra lanes of a padded input are
zero seeds on inputs the function ignores, so the result is the same
Jacobian. An empty list is `pickchunksize`. The two lists are the `Ref`s
`_CTSEM_HESSIAN_WIDTHS` and `_CTSEM_CURVATURE_WIDTHS`.

The Hessian takes one width, 12, for every parameter count. A padded lane
costs 5-12% of a small model's Hessian, and a binary model's 5 parameters at
width 8 cost 31% more (dev1, job H2: 7.4 -> 9.8 ms) -- milliseconds, a few
times a fit -- where each width the session had not run compiled the filter
and its reverse pass again, about 37 s on dev2 for a model whose every other
type the package image already holds, since the model's other types no longer
depend on its size. The curvature takes one width too, 4, whatever the unit's
size: a unit of more effects takes more sweeps. Its width is the unit's dual
type for everything the unit term differentiates -- the filter, its reverse
pass, a categorical item's moments -- so widths 8 and 12 for units of 5 to 24
effects, which the package image did not hold, compiled 1800 s of it on the
affect pilot's first evaluation (units of about 60 effects, interval-censored
items; review/PRECOMPILE-curvature-width-2026-10-26.md). Width 4 costs a unit
of 60 effects about 15 per cent more on its curvature (15 sweeps of 5 lanes'
work against 5 of 13), and one of one effect 2-8 per cent of an evaluation
against no padding.
"""
const _CTSEM_HESSIAN_WIDTHS = Ref(Int[12])
const _CTSEM_CURVATURE_WIDTHS = Ref(Int[4])

"""
    _ctsem_hessian_dual(objective)

The dual type `ctsem_hessian`'s sweeps run `objective`'s filter at: the tag of
`_CTSEMAdjointGradientOf` through `_CTSEMLeading`, at the Hessian's width. A
forward pass seeded in parameter directions that is never nested in another
derivative -- a row's scores -- runs at it rather than at a type of its own,
since every type is a compile of the filter and its reverse pass, and the
package image holds this one for any model a fit forms a Hessian for.
"""
_ctsem_hessian_dual(objective) = ForwardDiff.Dual{ForwardDiff.Tag{
    typeof(_CTSEMLeading(_CTSEMAdjointGradientOf(objective, 0), 0)),Float64},Float64,
    last(_CTSEM_HESSIAN_WIDTHS[])}

"""
The width a forward-mode Jacobian over `n` inputs takes from `widths`: the
smallest at least `pickchunksize(n)`, or the largest where none is.
"""
function _ctsem_dual_width(n::Integer, widths::Vector{Int})
    w0 = ForwardDiff.pickchunksize(n)
    isempty(widths) && return w0
    for w in widths
        w >= w0 && return w
    end
    return last(widths)
end

"""
`f` applied to the first `n` entries of its input: the padded Jacobian's function.

`whole` when `f` reads only those entries and ignores the padding after them,
so it can be handed the padded input itself rather than a copy of its head --
a new vector of duals on every sweep otherwise. A field and not a type
parameter: this type names the Jacobian's tag, so a second type would be a
second dual type for the filter to compile.
"""
struct _CTSEMLeading{F}
    f::F
    n::Int
    whole::Bool
end
_CTSEMLeading(f, n::Integer) = _CTSEMLeading(f, Int(n), false)
(g::_CTSEMLeading)(y::AbstractVector) =
    g.f(g.whole || length(y) == g.n ? y : y[1:g.n])

"""
    _ctsem_width_jacobian!(J, f, x, n, width)

The Jacobian of `f` over the first `n` entries of `x`, by forward mode at
`width` lanes a sweep, written into `J`, which has `length(x)` columns.
`length(x)` is `max(n, width)`: when `width` exceeds `n`, the entries past `n`
are padding, zero, and `J`'s columns past `n` come back zero. Always through
`_CTSEMLeading`, padded or not, so one `f` has one tag type at every width.
`whole`: `f` takes the padded input as it is (`_CTSEMLeading`).

`store`, a caller's per-task dictionary, keeps the `JacobianConfig` -- a vector
of duals, built per call otherwise -- by its type and the input's length. Its
tag is the one built here either way, so a cached config changes no dual type.
"""
function _ctsem_width_jacobian!(J::AbstractMatrix, f::F, x::AbstractVector,
    n::Integer, width::Integer; whole::Bool=false,
    store::Union{Nothing,AbstractDict}=nothing) where {F}
    # The width is a run-time value, so this call dispatches at run time
    # anyway; the barrier stops inference from also walking the whole sweep
    # for an unknown width, which the package image used to keep. `J` is
    # filled in place, so it is what is returned.
    _ctsem_barrier(_ctsem_width_jacobian!, J, _CTSEMLeading(f, Int(n), whole), x,
        Val(Int(width)), store)
    return J
end

function _ctsem_width_jacobian!(J::AbstractMatrix, g::G, x::AbstractVector{V},
    ::Val{W}, store) where {G,V,W}
    (store === nothing || !(x isa Vector)) && return ForwardDiff.jacobian!(J, g, x,
        ForwardDiff.JacobianConfig(g, x, ForwardDiff.Chunk{W}()))
    T = ForwardDiff.Tag{G,V}
    C = ForwardDiff.JacobianConfig{T,V,W,Vector{ForwardDiff.Dual{T,V,W}}}
    key = (:width_jacobian, C, length(x))
    cfg = get(store, key, nothing)
    if cfg === nothing
        cfg = ForwardDiff.JacobianConfig(g, x, ForwardDiff.Chunk{W}())
        store[key] = cfg
    end
    return ForwardDiff.jacobian!(J, g, x, cfg::C)
end

"""
`ForwardDiff.jacobian(f, x)` at the width `widths` gives `length(x)`.
`progress`, when given, is called as `progress(done, n)` after each sweep,
`done` the columns finished so far (`_ctsem_reported_jacobian`).
"""
function _ctsem_width_jacobian(f::F, x::AbstractVector, widths::Vector{Int};
        progress=nothing) where {F}
    n = length(x)
    width = _ctsem_dual_width(n, widths)
    # Behind the barrier for the reason `_ctsem_width_jacobian!` gives.
    width <= n && return _ctsem_barrier(_ctsem_reported_jacobian, _CTSEMLeading(f, n),
        x, ForwardDiff.Chunk{width}(), progress, n)::Matrix{eltype(x)}
    xp = vcat(x, zeros(eltype(x), width - n))
    J = _ctsem_barrier(_ctsem_reported_jacobian, _CTSEMLeading(f, n), xp,
        ForwardDiff.Chunk{width}(), progress, n)::Matrix{eltype(x)}
    return J[:, 1:n]
end

"""`g`, counting its calls and reporting after each; see `_ctsem_reported_jacobian`."""
struct _CTSEMSweepReporter{G,P}
    g::G
    progress::P
    width::Int
    n::Int
    calls::Base.RefValue{Int}
end
function (r::_CTSEMSweepReporter)(y::AbstractVector)
    out = r.g(y)
    r.calls[] += 1
    r.progress(min(r.calls[] * r.width, r.n), r.n)
    return out
end

"""
    _ctsem_reported_jacobian(g, x, chunk, progress, n)

`ForwardDiff.jacobian(g, x)` at `chunk`'s width, calling `progress(done, n)`
after each sweep when `progress` is not `nothing` -- `done` the columns of the
first `n` finished so far. ForwardDiff's chunk mode calls `g` once a sweep, so
a wrapper that counts its calls is the whole mechanism. It runs under `g`'s own
tag, with tag checking off (the one check the wrapper would otherwise fail), so
every dual the reverse pass sees is the type it has already compiled for: a
wrapper with a tag of its own would compile the filter and its reverse pass
again for a new dual type, tens of seconds on a large model.
"""
function _ctsem_reported_jacobian(g::G, x::AbstractVector, ::ForwardDiff.Chunk{W},
        progress, n::Integer) where {G,W}
    progress === nothing && return ForwardDiff.jacobian(g, x,
        ForwardDiff.JacobianConfig(g, x, ForwardDiff.Chunk{W}()))
    h = _CTSEMSweepReporter(g, progress, W, Int(n), Ref(0))
    cfg = ForwardDiff.JacobianConfig(h, x, ForwardDiff.Chunk{W}(),
        ForwardDiff.Tag(g, eltype(x)))
    return ForwardDiff.jacobian(h, x, cfg, Val{false}())
end

"""
    _ctsem_barrier(f, args...)

Call `f(args...)` without letting inference look into `f`: the call dispatches
at run time, and `f` is compiled only when a call actually reaches it.

For code a model can reach and seldom runs: forward mode for the categorical
moments' Jacobian (`_binary_moment_dual!`), which `_binary_moment_jacobian`
falls back to only for an asymptote item, a variance too small to integrate or
a mode solve that used its budget. The price is one dynamic dispatch per call and a result that
inference cannot see, so a caller that uses the result asserts its type.

Code that a whole class of models never reaches is better gated on a type, as
the categorical update is (`_ekf_categorical_call`): that folds at compile
time and costs nothing per call, where this barrier in front of it cost one
dispatch per categorical row on each pass.
"""
@inline _ctsem_barrier(f::F, args::Vararg{Any,N}) where {F,N} =
    Base.inferencebarrier(f)(args...)
