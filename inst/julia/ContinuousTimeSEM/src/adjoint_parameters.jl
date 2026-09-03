using ForwardDiff
using ComponentArrays

################################################################################
# Parameter-layer adjoint
################################################################################
#
# The reverse EKF pass (adjoint_ekf.jl) produces a cotangent for the
# *materialized* model matrices -- DRIFT, JAx, CINT, LAMBDA, Jy, T0MEANS, and
# the covariance square-root matrices -- laid out exactly like `ws.all_params`.
# This file turns that into a gradient with respect to the free parameter
# vector `values` that R hands in.
#
# There are three layers to undo, in the order the forward pass applies them:
#
#   1. `_materialize_subject_values!`: subject_values = values, plus TI
#      predictor effects `subject_values[par] += values[coef] * tipred[pred]`.
#      Linear, so its pullback is exact and trivial.
#
#   2. `_materialize_all_params!`: for each mutable position `idx`,
#      `all_params[idx] = regular_transform_idx(subject_values)`. By
#      construction on the R side (`ctJuliaBackend.R` rewrites every bare
#      `param` in a transform string to `param[parnumber]`) each of these
#      reads exactly one entry of `subject_values`, so the Jacobian is
#      diagonal-in-parameter and one scalar derivative per position suffices.
#      That "by construction" claim is *verified*, not assumed -- see
#      `_ctsem_regular_transform_supports` below, which throws if it fails.
#
#   3. `apply_complex_transforms_at_indices!`: state-dependent cells,
#      re-evaluated at every row (and every prediction substep). These read the
#      current state and other already-materialized parameter cells, and write
#      back into the same `all_params` buffer they read from. That makes them
#      genuine in-place overwrites, and the reverse pass has to treat them as
#      such: walking a group backwards, each transform *consumes* the cotangent
#      sitting on the cell it wrote (zeroing it, because the value that was
#      there before the write never reached the likelihood) and scatters it to
#      the cells and state entries it read. Anything still sitting on a cell
#      when the reverse pass reaches the start of the filter belongs to layer 2.
#
# The read sets in layers 2 and 3 are discovered exactly, by evaluating each
# transform once against a parameter vector that records which indices are
# read, rather than by parsing expression strings or by probing derivatives
# numerically (which would silently miss a dependency whose derivative happens
# to vanish at the probe point).

"""
    RecordingVector(parent)

An `AbstractVector` that forwards to `parent` and records which linear indices
have been read.

Used only during adjoint-workspace construction, to discover exactly which
parameter cells each transform expression touches. Wrapping it in a
`ComponentVector` works because every `ctx.pars.DRIFT[i,j]`-style access in a
transform expression ultimately bottoms out in a linear `getindex` on the
backing vector.
"""
mutable struct RecordingVector{T} <: AbstractVector{T}
    parent::Vector{T}
    touched::Vector{Bool}
end

RecordingVector(parent::Vector{T}) where {T} =
    RecordingVector{T}(parent, fill(false, length(parent)))

Base.size(v::RecordingVector) = size(v.parent)
Base.IndexStyle(::Type{<:RecordingVector}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::RecordingVector, i::Int)
    v.touched[i] = true
    return v.parent[i]
end
Base.@propagate_inbounds Base.setindex!(v::RecordingVector, x, i::Int) =
    (v.parent[i] = x)

_recorded_indices(v::RecordingVector) = findall(v.touched)
_reset_recording!(v::RecordingVector) = (fill!(v.touched, false); v)

"""
    _ctsem_regular_transform_supports(sp, nvalues)

Return, for each mutable position in `sp`, the single `subject_values` index
its regular transform reads.

Throws if any transform reads something other than exactly its own
`parnumber`. That would not be a Julia bug but an R-side change to how
transform strings are rendered, and it must fail loudly here rather than
produce a gradient that silently drops a term.
"""
function _ctsem_regular_transform_supports(sp::EKFParameters, nvalues::Integer)
    probe = RecordingVector(collect(range(0.11, step=0.017, length=Int(nvalues))))
    supports = Vector{Int}(undef, length(sp.regular_transforms))
    tf_idx = 0
    for idx in eachindex(sp.mutables)
        sp.mutables[idx] || continue
        tf_idx += 1
        _reset_recording!(probe)
        sp.regular_transforms[tf_idx](probe)
        read = _recorded_indices(probe)
        expected = sp.parnumber[tf_idx]
        if read != [expected]
            throw(ArgumentError(string(
                "adjoint: regular transform for flattened parameter position ", idx,
                " reads subject_values indices ", read, ", expected exactly [",
                expected, "]. The adjoint's parameter-layer pullback assumes one ",
                "free parameter per transform, which is how ctJuliaBackend.R ",
                "renders them; a multi-parameter transform needs the pullback in ",
                "adjoint_parameters.jl generalised before it can be used.")))
        end
        supports[tf_idx] = expected
    end
    return supports
end

"""
    _ctsem_complex_transform_supports(transforms, indices, sp, state_dim, ntdpred)

Return a `Vector{Vector{Int}}` giving, for each state-dependent transform, the
`all_params` indices it reads.

The state dependence is handled separately (every complex transform is
differentiated with respect to the whole state vector, which is short), so
only the parameter-cell reads need discovering here.
"""
function _ctsem_complex_transform_supports(transforms, indices::AbstractVector{Int},
    sp::EKFParameters, state_dim::Int, ntdpred::Int)
    isempty(indices) && return Vector{Int}[]
    nall = length(sp.mutables)
    # Probe values are distinct and non-degenerate so that an expression like
    # `PARS[1,1] * pars.X` cannot be read-free by accident. Reads are recorded
    # structurally, so the actual numbers only need to avoid domain errors.
    probe = RecordingVector(collect(range(0.31, step=0.013, length=nall)))
    pars = ComponentVector(probe, sp.parameter_axis)
    state = reshape(collect(range(0.23, step=0.011, length=state_dim)), state_dim, 1)
    tdpreds = collect(range(0.19, step=0.007, length=ntdpred))
    tipreds = Float64[]
    ctx = CTSEMRowContext(state, pars, tdpreds, tipreds, 0.5, 0.25, 1, 2)
    supports = Vector{Vector{Int}}(undef, length(indices))
    for k in eachindex(indices)
        _reset_recording!(probe)
        transforms[k](ctx)
        supports[k] = _recorded_indices(probe)
    end
    return supports
end

################################################################################
# Layer 1 and 2 pullbacks
################################################################################

"""
    _ctsem_regular_pullback!(subject_values_bar, all_params_bar, subject_values, sp, supports, scratch)

Push the cotangent on `all_params` back through the regular transforms into
`subject_values_bar`.

`scratch` is a `Vector{ForwardDiff.Dual{...,1}}` mirror of `subject_values`
whose entries all carry a zero partial; each transform is evaluated with a
single seeded entry, which is why one scalar derivative per mutable position
is enough.
"""
function _ctsem_regular_pullback!(subject_values_bar::AbstractVector,
    all_params_bar::AbstractVector, subject_values::AbstractVector,
    sp::EKFParameters, supports::AbstractVector{Int}, scratch::AbstractVector)
    tf_idx = 0
    @inbounds for idx in eachindex(sp.mutables)
        sp.mutables[idx] || continue
        tf_idx += 1
        cotangent = all_params_bar[idx]
        iszero(cotangent) && continue
        pn = supports[tf_idx]
        base = subject_values[pn]
        scratch[pn] = _seed_dual(scratch[pn], base, true)
        derivative = _partial1(sp.regular_transforms[tf_idx](scratch))
        scratch[pn] = _seed_dual(scratch[pn], base, false)
        subject_values_bar[pn] += cotangent * derivative
    end
    return subject_values_bar
end

@inline _seed_dual(::D, value, seed::Bool) where {D<:ForwardDiff.Dual} =
    D(value, ForwardDiff.Partials((seed ? one(value) : zero(value),)))

"""
    _ctsem_ti_pullback!(values_bar, subject_values_bar, sp, tipreds, values)

Undo `_materialize_subject_values!`: the identity copy plus each TI predictor
effect `subject_values[par] += values[coef] * tipred[pred]`.

Two methods, dispatched on the *stored* `tipreds` type -- exactly mirroring
`_ctsem_tipred_vector` in parameter_transforms.jl:

  * `tipreds::AbstractVector` -- every predictor value here is data (a
    constant as far as the gradient is concerned), so the product rule only
    has one variable factor: the coefficient.
  * `tipreds::TIMissingRecipe` -- one or more predictor values are
    *themselves* raw parameters (`_ctsem_tipred_vector`'s other method
    substitutes `values[parameter_index[k]]` at each missing cell before the
    row is ever read). The product rule then has two variable factors for
    every TI effect that reads one of those cells:
    `d(coef * tipred)/d(values) = tipred * d(coef) + coef * d(tipred)`, so
    this method accumulates into `values_bar[coefficient]` exactly as the
    plain method does *and* into `values_bar[parameter_index[k]]`, using
    `values[coefficient]` -- the raw trial point, not `subject_values`, which
    only equals it when no TI effect's coefficient is itself the target of
    another TI effect. Two or more TI effects reading the same missing cell
    (`predictor` shared across several `i`) each contribute their own term to
    the same `values_bar[parameter_index[k]]` slot, and the `+=` below sums
    them rather than overwriting, as it must.
"""
function _ctsem_ti_pullback!(values_bar::AbstractVector,
    subject_values_bar::AbstractVector, sp::EKFParameters, tipreds::AbstractVector,
    values::AbstractVector)
    @inbounds for i in eachindex(values_bar)
        values_bar[i] += subject_values_bar[i]
    end
    @inbounds for i in eachindex(sp.ti_parameter_indices)
        parameter = sp.ti_parameter_indices[i]
        predictor = sp.ti_predictor_indices[i]
        coefficient = sp.ti_coefficient_indices[i]
        values_bar[coefficient] += subject_values_bar[parameter] * tipreds[predictor]
    end
    return values_bar
end

function _ctsem_ti_pullback!(values_bar::AbstractVector,
    subject_values_bar::AbstractVector, sp::EKFParameters, spec::TIMissingRecipe,
    values::AbstractVector)
    tipred_vec = _ctsem_tipred_vector(spec, values)
    @inbounds for i in eachindex(values_bar)
        values_bar[i] += subject_values_bar[i]
    end
    @inbounds for i in eachindex(sp.ti_parameter_indices)
        parameter = sp.ti_parameter_indices[i]
        predictor = sp.ti_predictor_indices[i]
        coefficient = sp.ti_coefficient_indices[i]
        cotangent = subject_values_bar[parameter]
        values_bar[coefficient] += cotangent * tipred_vec[predictor]
        # Second product-rule term, present only when this TI effect's own
        # predictor column is one of this subject's sampled cells.
        k = findfirst(==(predictor), spec.predictor_index)
        if k !== nothing
            values_bar[spec.parameter_index[k]] += cotangent * values[coefficient]
        end
    end
    return values_bar
end

################################################################################
# Layer 3 pullback
################################################################################

"""
    _ctsem_complex_group_pullback!(all_params_bar, state_bar, transforms, indices,
                                   supports, ctx, relevant, dual_ctx)

Reverse one group of state-dependent transforms (the `predict`, `td`, or
`update` group at one row).

Walks the group backwards, because within a group a later transform may read a
cell an earlier one has already written. For each transform this:

  1. takes the cotangent currently on the cell it wrote,
  2. zeroes that cotangent -- the cell was *overwritten*, so whatever value was
     there beforehand did not reach the likelihood through this path,
  3. adds `cotangent * ∂f/∂input` to every parameter cell and state entry the
     transform reads (including the written cell itself, which is legitimate:
     a self-referencing expression reads the pre-write value, and step 2 has
     already cleared the post-write cotangent).

Derivatives come from `dual_ctx`, a `ForwardDiff` mirror of `ctx` with a
single seeded input at a time. One-at-a-time seeding keeps this simple and
exact; each transform costs `state_dim + length(support)` scalar evaluations
of a short expression.

**Each transform's derivative is evaluated at the values that transform
actually saw**, which are not the group's final values. Transform `k` runs
after `1..k-1` have already written their cells but before `k+1..K` have
written theirs, so it reads the current row's values for the former and the
*previous* row's for the latter. The group is therefore replayed forward once
from `params_before` to recover what each transform wrote, and the reverse walk
then unwinds those writes one at a time as it descends -- which costs `O(1)`
per step rather than a fresh full-length sync.

Note the cotangent *routing* needs no such care and is correct either way: a
contribution that transform `k` sends to a cell written by a later transform
lands on that cell after it has already been consumed and zeroed, so it flows
onward to earlier tape entries -- which is exactly right, because the value
`k` read came from earlier in time.
"""
function _ctsem_complex_group_pullback!(all_params_bar::AbstractVector,
    state_bar::AbstractVector, transforms, indices::AbstractVector{Int},
    supports::AbstractVector{Vector{Int}}, ctx::CTSEMRowContext,
    relevant::AbstractVector{Int}, dual_ctx)
    isempty(indices) && return nothing
    nk = length(indices)
    working = getdata(ctx.pars)

    # Keep the pre-group values of the cells this group writes, because the
    # replay below overwrites them and the reverse walk needs to put them back
    # one at a time.
    saved = _group_saved_buffer(dual_ctx, nk)
    @inbounds for k in 1:nk
        saved[k] = working[indices[k]]
    end

    # Replay the group forward over `ctx`'s working buffer (seeded from the
    # compact record by the caller) to recover what each transform wrote.
    written = _group_write_buffer(dual_ctx, nk)
    @inbounds for k in 1:nk
        value = transforms[k](ctx)
        written[k] = value
        working[indices[k]] = value
    end

    # Start the reverse walk at transform K, whose inputs are the pre-group
    # values with every earlier transform's write applied. Only the relevant
    # entries are synced: a transform cannot read anything else by
    # construction, and syncing the whole vector per group per row is what
    # this compaction exists to avoid.
    # `working` holds *post*-replay values now, so every written cell has to be
    # set explicitly rather than relying on the sync: cells 1..K-1 to what
    # their transforms wrote, and cell K back to its pre-group value, since
    # transform K ran before its own write landed.
    _sync_dual_context!(dual_ctx, ctx, working, relevant)
    @inbounds for k in 1:(nk - 1)
        _set_dual_value!(dual_ctx, indices[k], written[k])
    end
    @inbounds _set_dual_value!(dual_ctx, indices[nk], saved[nk])

    @inbounds for k in nk:-1:1
        idx = indices[k]
        cotangent = all_params_bar[idx]
        all_params_bar[idx] = zero(cotangent)
        if !iszero(cotangent)
            transform = transforms[k]
            for j in supports[k]
                all_params_bar[j] += cotangent * _dual_param_derivative(dual_ctx, transform, j)
            end
            for s in eachindex(state_bar)
                state_bar[s] += cotangent * _dual_state_derivative(dual_ctx, transform, s)
            end
        end
        # Step back one transform: undo transform k-1's write so that
        # transform k-1 sees the inputs it actually had.
        k > 1 && _set_dual_value!(dual_ctx, indices[k - 1], saved[k - 1])
    end
    return nothing
end

"""
    CTSEMDualContext

A `ForwardDiff.Dual`-typed mirror of a `CTSEMRowContext`, reused across rows.

Holds a full-length dual parameter buffer and a dual state so that seeding one
input costs `O(1)` rather than rebuilding the whole context. Entries outside a
transform's recorded read set are never consulted, so only the read set and the
state need to be kept in sync with the primal.
"""
mutable struct CTSEMDualContext{D,PD,ST,T}
    data::Vector{D}
    pars::PD
    state::ST
    row_context::Any
    # Scratch for the forward replay in `_ctsem_complex_group_pullback!`,
    # grown on demand to the largest group seen. Typed `T` -- the scalar the
    # surrounding reverse pass works in -- and not `Float64`. It used to be
    # `Float64` on the reasoning that the adjoint is only ever entered with
    # doubles from R, which stopped being true the moment the gradient itself
    # became something to differentiate: `ctsem_hessian` runs this whole pass
    # at `T = ForwardDiff.Dual`, and parameter values round-trip through here.
    written::Vector{T}
    # Pre-group values of the cells a group writes, for the reverse walk.
    saved::Vector{T}
end

function CTSEMDualContext(::Type{T}, sp::EKFParameters, state_like) where {T}
    D = ForwardDiff.Dual{Nothing,T,1}
    zero_dual = D(zero(T), ForwardDiff.Partials((zero(T),)))
    data = fill(zero_dual, length(sp.mutables))
    pars = ComponentVector(data, sp.parameter_axis)
    state = similar(state_like, D)
    fill!(state, zero_dual)
    return CTSEMDualContext{D,typeof(pars),typeof(state),T}(data, pars, state, nothing,
        T[], T[])
end

"""Scratch vector of at least `n` slots for one group's replayed writes."""
function _group_write_buffer(dual_ctx::CTSEMDualContext, n::Int)
    length(dual_ctx.written) < n && resize!(dual_ctx.written, n)
    return dual_ctx.written
end

"""Scratch vector of at least `n` slots for one group's pre-group values."""
function _group_saved_buffer(dual_ctx::CTSEMDualContext, n::Int)
    length(dual_ctx.saved) < n && resize!(dual_ctx.saved, n)
    return dual_ctx.saved
end

"""Set the value (keeping a zero partial) of one entry of the dual mirror."""
@inline function _set_dual_value!(dual_ctx::CTSEMDualContext, j::Int, value)
    dual_ctx.data[j] = _seed_dual(dual_ctx.data[j], value, false)
    return nothing
end

"""
Refresh the dual mirror's values (with zero partials) from the primal context,
and rebuild the dual row context so the seeded derivative calls below see this
row's time, dt, and predictors.
"""
function _sync_dual_context!(dual_ctx::CTSEMDualContext, ctx::CTSEMRowContext,
    primal_data::AbstractVector=getdata(ctx.pars),
    relevant::AbstractVector{Int}=eachindex(dual_ctx.data))
    @inbounds for i in relevant
        dual_ctx.data[i] = _seed_dual(dual_ctx.data[i], primal_data[i], false)
    end
    @inbounds for i in eachindex(dual_ctx.state)
        dual_ctx.state[i] = _seed_dual(dual_ctx.state[i], ctx.state[i], false)
    end
    dual_ctx.row_context = CTSEMRowContext(dual_ctx.state, dual_ctx.pars, ctx.tdpreds,
        ctx.tipreds, ctx.time, ctx.dt, ctx.subject, ctx.row)
    return dual_ctx
end

# A transform that happens not to depend on any dual input returns a plain
# real; treat that as a zero partial rather than letting `partials` throw.
@inline _partial1(x::ForwardDiff.Dual) = ForwardDiff.partials(x, 1)
@inline _partial1(x::Real) = zero(x)

function _dual_param_derivative(dual_ctx::CTSEMDualContext, transform, j::Int)
    saved = dual_ctx.data[j]
    dual_ctx.data[j] = _seed_dual(saved, ForwardDiff.value(saved), true)
    derivative = _partial1(transform(dual_ctx.row_context))
    dual_ctx.data[j] = saved
    return derivative
end

function _dual_state_derivative(dual_ctx::CTSEMDualContext, transform, s::Int)
    saved = dual_ctx.state[s]
    dual_ctx.state[s] = _seed_dual(saved, ForwardDiff.value(saved), true)
    derivative = _partial1(transform(dual_ctx.row_context))
    dual_ctx.state[s] = saved
    return derivative
end
