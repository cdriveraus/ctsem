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
    _custom_abs(x::ForwardDiff.Dual)

Return the absolute value of the primal part of a dual number.
"""
@inline _custom_abs(x::ForwardDiff.Dual) = abs(ForwardDiff.value(x))

@inline _val(::Val{N}) where {N} = N

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
    return add_diag!(mat, val, Val(min(size(mat, 1), size(mat, 2))))
end

function add_diag!(mat, val::Number, ::Val{d}) where {d}
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
    return add_diag!(mat, vals, Val(length(vals)))
end

function add_diag!(mat, vals::AbstractVector, ::Val{d}) where {d}
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
    _matvec_mul!(y, A, x::AbstractVector)

Compute `y .= A * x` without allocating temporaries.

This manual implementation is used for `SubArray`-backed EKF buffers where
`mul!` is not suitable.
"""
function _matvec_mul!(y::AbstractVector, A::AbstractMatrix, x::AbstractVector)
    return _matvec_mul!(y, A, x, Val(size(A, 1)), Val(size(A, 2)))
end

function _matvec_mul!(
    y::AbstractVector,
    A::AbstractMatrix,
    x::AbstractVector,
    ::Val{rows},
    ::Val{cols},
) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && length(y) == rows || throw(DimensionMismatch("A*y output size mismatch"))
        size(A, 2) == cols && length(x) == cols || throw(DimensionMismatch("A*x input size mismatch"))
    end
    @inbounds for row in 1:rows
        acc = zero(eltype(y))
        for col in 1:cols
            acc += A[row, col] * x[col]
        end
        y[row] = acc
    end
    return y
end

"""
    _matvec_mul!(y, A, x::AbstractMatrix)

Compute `y .= A * x[:, 1]` without allocating temporaries.

The matrix input must have exactly one column.
"""
function _matvec_mul!(y::AbstractVector, A::AbstractMatrix, x::AbstractMatrix)
    return _matvec_mul!(y, A, x, Val(size(A, 1)), Val(size(A, 2)))
end

function _matvec_mul!(
    y::AbstractVector,
    A::AbstractMatrix,
    x::AbstractMatrix,
    ::Val{rows},
    ::Val{cols},
) where {rows, cols}
    @boundscheck begin
        size(A, 1) == rows && length(y) == rows || throw(DimensionMismatch("A*y output size mismatch"))
        size(x, 2) == 1 || throw(DimensionMismatch("x must be a column matrix"))
        size(A, 2) == cols && size(x, 1) == cols || throw(DimensionMismatch("A*x input size mismatch"))
    end
    @inbounds for row in 1:rows
        acc = zero(eltype(y))
        for col in 1:cols
            acc += A[row, col] * x[col, 1]
        end
        y[row] = acc
    end
    return y
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` for BLAS-compatible strided matrices.
"""
function _mul_right_transpose!(C::StridedMatrix{T}, A::StridedMatrix{T}, B::StridedMatrix{T}) where {T<:LinearAlgebra.BlasFloat}
    return _mul_right_transpose!(C, A, B, Val(size(C, 1)), Val(size(C, 2)), Val(size(A, 2)))
end

function _mul_right_transpose!(
    C::StridedMatrix{T},
    A::StridedMatrix{T},
    B::StridedMatrix{T},
    ::Val{rows},
    ::Val{cols},
    ::Val{inner},
) where {T<:LinearAlgebra.BlasFloat, rows, cols, inner}
    @boundscheck begin
        size(C, 1) == rows && size(A, 1) == rows || throw(DimensionMismatch("A row count must match C row count"))
        size(C, 2) == cols && size(B, 1) == cols || throw(DimensionMismatch("B row count must match C column count"))
        size(A, 2) == inner && size(B, 2) == inner || throw(DimensionMismatch("A and B must have the same column count"))
    end
    # Small products by hand, BLAS above the threshold: see `_ctsem_mul!`
    # and the note at the predict step in kalman_filters.jl.
    return _ctsem_mulNT!(C, A, B)
end

"""
    _mul_right_transpose!(C, A, B)

Compute `C .= A * B'` with explicit loops for generic numeric matrices.
"""
function _mul_right_transpose!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T<:Number}
    return _mul_right_transpose!(C, A, B, Val(size(C, 1)), Val(size(C, 2)), Val(size(A, 2)))
end

function _mul_right_transpose!(
    C::AbstractMatrix{T},
    A::AbstractMatrix{T},
    B::AbstractMatrix{T},
    ::Val{rows},
    ::Val{cols},
    ::Val{inner},
) where {T<:Number, rows, cols, inner}
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
    return _isapprox_matrix_noalloc(A, B, Val(size(A, 1)), Val(size(A, 2)); atol=atol, rtol=rtol)
end

function _isapprox_matrix_noalloc(
    A::AbstractMatrix,
    B::AbstractMatrix,
    ::Val{rows},
    ::Val{cols};
    atol::Real=0.0,
    rtol::Real=_isapprox_default_rtol(promote_type(eltype(A), eltype(B))),
) where {rows, cols}
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

function _matrix_equal_noalloc(A::AbstractMatrix, B::AbstractMatrix, ::Val{rows}, ::Val{cols}) where {rows, cols}
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
@inline _can_reuse_same_exponential(A::AbstractMatrix{<:ForwardDiff.Dual}, B::AbstractMatrix{<:ForwardDiff.Dual}, dim::Val{d}) where {d} =
    _matrix_equal_noalloc(A, B, dim, dim)

"""
    _can_reuse_same_exponential(A, B)

Return whether two ordinary matrices are close enough to share an exponential.
"""
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix) = _isapprox_matrix_noalloc(A, B)
@inline _can_reuse_same_exponential(A::AbstractMatrix, B::AbstractMatrix, dim::Val{d}) where {d} =
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

"""`ForwardDiff.gradient(f, x)` under `CTSEMNestedTag`; see there."""
_ctsem_nested_gradient(f, x::AbstractArray) = ForwardDiff.gradient(f, x,
    ForwardDiff.GradientConfig(f, x, ForwardDiff.Chunk(x), _ctsem_nested_tag(x)),
    Val{false}())

"""
Dual widths for the forward-mode Jacobians whose input length is a property of
the model: the Hessian's (the parameter count) and the Laplace unit
curvature's (the unit's random-effect dimension).

`ForwardDiff.pickchunksize(n)` is `n` itself up to 12, and a dual width is a
type, so every parameter count compiled the whole filter and reverse pass
again: about 30 s on dev2 for one more free parameter on a model the session,
or the package image, had already compiled. Rounding the width up to the next
of a few buckets makes every count in a bucket one type. Above 12,
`pickchunksize`'s width is rounded up the same way, which never adds a sweep.
The extra lanes of a padded input are zero seeds on inputs the function
ignores, so the result is the same Jacobian. An empty list is `pickchunksize`.
`ctsem_set_dual_widths!` sets both lists.

The Hessian's list is empty by default: measured on dev1 (job H2), a padded
lane costs 5-12% of a small model's Hessian, and a binary model's 5 parameters
at width 8 cost 31% more (7.4 -> 9.8 ms; even buckets 2, 4, ..., 12: +8%, and
+12% for 7 Gaussian parameters at 8), so no bucketing both kept a Hessian
within 10% and made every parameter count free. The curvature's buckets cost
a Laplace evaluation 2% (one or two random effects, Gaussian) to 8% (one,
ordinal) and make one to four random effects per unit a single type.
"""
const _CTSEM_HESSIAN_WIDTHS = Ref(Int[])
const _CTSEM_CURVATURE_WIDTHS = Ref(Int[4, 8, 12])

export ctsem_set_dual_widths!
function ctsem_set_dual_widths!(; hessian=nothing, curvature=nothing)
    hessian === nothing || (_CTSEM_HESSIAN_WIDTHS[] = sort!(Int.(collect(hessian))))
    curvature === nothing || (_CTSEM_CURVATURE_WIDTHS[] = sort!(Int.(collect(curvature))))
    return (hessian = copy(_CTSEM_HESSIAN_WIDTHS[]),
        curvature = copy(_CTSEM_CURVATURE_WIDTHS[]))
end

"""The width a forward-mode Jacobian over `n` inputs takes from `widths`."""
function _ctsem_dual_width(n::Integer, widths::Vector{Int})
    w0 = ForwardDiff.pickchunksize(n)
    for w in widths
        w >= w0 && return w
    end
    return w0
end

"""`f` applied to the first `n` entries of its input: the padded Jacobian's function."""
struct _CTSEMLeading{F}
    f::F
    n::Int
end
(g::_CTSEMLeading)(y::AbstractVector) = g.f(length(y) == g.n ? y : y[1:g.n])

"""
    _ctsem_width_jacobian!(J, f, x, n, width)

The Jacobian of `f` over the first `n` entries of `x`, by forward mode at
`width` lanes a sweep, written into `J`, which has `length(x)` columns.
`length(x)` is `max(n, width)`: when `width` exceeds `n`, the entries past `n`
are padding, zero, and `J`'s columns past `n` come back zero. Always through
`_CTSEMLeading`, padded or not, so one `f` has one tag type at every width.
"""
function _ctsem_width_jacobian!(J::AbstractMatrix, f::F, x::AbstractVector,
    n::Integer, width::Integer) where {F}
    g = _CTSEMLeading(f, Int(n))
    cfg = ForwardDiff.JacobianConfig(g, x, ForwardDiff.Chunk{Int(width)}())
    return ForwardDiff.jacobian!(J, g, x, cfg)
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
    width <= n && return _ctsem_reported_jacobian(_CTSEMLeading(f, n), x,
        ForwardDiff.Chunk{width}(), progress, n)
    xp = vcat(x, zeros(eltype(x), width - n))
    J = _ctsem_reported_jacobian(_CTSEMLeading(f, n), xp,
        ForwardDiff.Chunk{width}(), progress, n)
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
falls back to only for a censored row, an asymptote item or a mode solve that
used its budget. The price is one dynamic dispatch per call and a result that
inference cannot see, so a caller that uses the result asserts its type.

Code that a whole class of models never reaches is better gated on a type, as
the categorical update is (`_ekf_categorical_call`): that folds at compile
time and costs nothing per call, where this barrier in front of it cost one
dispatch per categorical row on each pass.
"""
@inline _ctsem_barrier(f::F, args::Vararg{Any,N}) where {F,N} =
    Base.inferencebarrier(f)(args...)
