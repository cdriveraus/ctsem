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
#      `all_params[idx] = regular_transform_idx(subject_values)`. Almost every
#      one of these reads a single entry of `subject_values` -- R renders the
#      cell's transform with its own `param[parnumber]` substituted in -- so
#      the Jacobian is diagonal-in-parameter and one scalar derivative per
#      position suffices. The exception is a T0MEANS or T0VAR cell written as
#      an expression over several PARS parameters, which
#      `.ctJuliaResolveStaticRefs` composes into one transform reading each of
#      them. The read set is *discovered*, not assumed -- see
#      `_ctsem_regular_transform_supports` below -- and the pullback takes one
#      scalar derivative per index in it.
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

Return, for each mutable position in `sp`, the `subject_values` indices its
regular transform reads.

Discovered with a `RecordingVector` rather than assumed, so a cell composed
from several raw parameters (`.ctJuliaResolveStaticRefs` builds those for
T0MEANS and T0VAR) reports all of them and `_ctsem_regular_pullback!` pushes a
term back through each.

Two things still throw, because both would be R-side rendering mistakes that
the gradient would otherwise absorb silently:

  * a transform that reads *nothing*, which means the cell was marked mutable
    but its text has no `param[...]` in it at all, so its raw coordinate has
    no path to the likelihood; and
  * a transform whose own `parnumber` is not among the indices it reads, which
    is the part of the old one-parameter assertion that catches a substitution
    landing on the wrong slot.
"""
function _ctsem_regular_transform_supports(sp::EKFParameters, nvalues::Integer)
    probe = RecordingVector(collect(range(0.11, step=0.017, length=Int(nvalues))))
    supports = Vector{Vector{Int}}(undef, length(sp.regular_transforms))
    tf_idx = 0
    for idx in eachindex(sp.mutables)
        sp.mutables[idx] || continue
        tf_idx += 1
        _reset_recording!(probe)
        sp.regular_transforms[tf_idx](probe)
        read = _recorded_indices(probe)
        expected = sp.parnumber[tf_idx]
        if isempty(read) || !(expected in read)
            throw(ArgumentError(string(
                "adjoint: regular transform for flattened parameter position ", idx,
                " reads subject_values indices ", read, ", which does not include ",
                "its own parameter number ", expected, ". A transform is rendered ",
                "by ctJuliaBackend.R from the cell's own parameter, plus any PARS ",
                "parameters a composed T0MEANS/T0VAR expression references; a ",
                "transform that reads neither is a rendering fault, and the ",
                "parameter-layer pullback would silently drop the term.")))
        end
        supports[tf_idx] = read
    end
    return supports
end

"""
A cell may read the subject's time independent predictors (`ctx.tipreds[j]`),
so the support probe hands it this: a plain number at whatever index the cell
names, without needing to know how many predictors the model has, and a record
of which it read. A predictor is data to the cell unless the subject's value
of it is sampled (`TIMissingRecipe`), and then the reverse pass needs the
cell's derivative with respect to it (`_group_reverse_walk!`). It had an empty
vector, and a cell reading a predictor failed there.
"""
struct _CTSEMProbePredictors <: AbstractVector{Float64}
    read::Vector{Int}
end
_CTSEMProbePredictors() = _CTSEMProbePredictors(Int[])
Base.size(::_CTSEMProbePredictors) = (0,)
function Base.getindex(p::_CTSEMProbePredictors, i::Int)
    i in p.read || push!(p.read, i)
    return 0.17 + 0.003 * i
end

"""
    _ctsem_complex_transform_supports(transforms, indices, sp, state_dim, ntdpred)

For each state-dependent transform, the `all_params` indices it reads and the
TI predictor indices it reads, as two `Vector{Vector{Int}}`.

The state dependence is handled separately (every complex transform is
differentiated with respect to the whole state vector, which is short), so
only the parameter-cell and predictor reads need discovering here.
"""
function _ctsem_complex_transform_supports(transforms, indices::AbstractVector{Int},
    sp::EKFParameters, state_dim::Int, ntdpred::Int)
    isempty(indices) && return Vector{Int}[], Vector{Int}[]
    nall = length(sp.mutables)
    # Probe values are distinct and non-degenerate so that an expression like
    # `PARS[1,1] * pars.X` cannot be read-free by accident. Reads are recorded
    # structurally, so the actual numbers only need to avoid domain errors.
    probe = RecordingVector(collect(range(0.31, step=0.013, length=nall)))
    pars = CTSEMPars(probe, sp.parameter_axis)
    state = reshape(collect(range(0.23, step=0.011, length=state_dim)), state_dim, 1)
    tdpreds = collect(range(0.19, step=0.007, length=ntdpred))
    predictors = _CTSEMProbePredictors()
    ctx = CTSEMRowContext(state, pars, tdpreds, predictors, 0.5, 0.25, 1, 2)
    supports = Vector{Vector{Int}}(undef, length(indices))
    tipred_supports = Vector{Vector{Int}}(undef, length(indices))
    for k in eachindex(indices)
        _reset_recording!(probe)
        empty!(predictors.read)
        transforms[k](ctx)
        supports[k] = _recorded_indices(probe)
        tipred_supports[k] = sort(predictors.read)
    end
    return supports, tipred_supports
end

################################################################################
# Layer 1 and 2 pullbacks
################################################################################

"""
    _ctsem_regular_pullback!(subject_values_bar, all_params_bar, subject_values, sp, supports, scratch)

Push the cotangent on `all_params` back through the regular transforms into
`subject_values_bar`.

`scratch` is a `Vector{ForwardDiff.Dual{...,1}}` mirror of `subject_values`,
seeded one entry at a time: the transform is evaluated once per index in its
support, and each evaluation reads a partial of 1 at that index and 0
everywhere else. A cell materialised from a single raw parameter -- all of them
except a T0MEANS or T0VAR cell composed from several PARS parameters -- has a
one-element support, so the common case is the one seeded evaluation this
always was, and the loop is neither a per-cell allocation nor a Jacobian.

**`scratch` carries no values between calls.** Seeding an entry writes its
*value* as well as its partial, so a single-parameter transform never needed
any other entry of `scratch` to be current -- and none of them were. A
composed transform reads every index in its support on each of its
evaluations, so all of them are refreshed up front. Without that, the
unseeded ones were read at whatever was last left in the buffer, which is a
gradient wrong by a term rather than an error: on the composed-cell scenario
in `test_adjoint_gradient_validation.jl` it disagreed with ForwardDiff by 4.5%
while every other gate still passed.
"""
function _ctsem_regular_pullback!(subject_values_bar::AbstractVector,
    all_params_bar::AbstractVector, subject_values::AbstractVector,
    sp::EKFParameters, supports::AbstractVector{Vector{Int}}, scratch::AbstractVector)
    _regular_pullback_groups!(subject_values_bar, all_params_bar, subject_values,
        sp.regular_groups, sp.regular_group_cells, sp.regular_group_tf,
        supports, scratch, 1)
    return subject_values_bar
end

# Grouped for the same reason the forward pass is -- see
# `_apply_regular_groups!` -- and it matters more here, since this runs once per
# index in each cell's support rather than once per cell. Cells are visited in
# template order rather than in cell order, so the terms landing on one raw
# coordinate are summed in a different order than before; nothing else about the
# result changes.
@inline _regular_pullback_groups!(subject_values_bar, all_params_bar,
    subject_values, ::Tuple{}, cells, tfs, supports, scratch, g::Int) = nothing

@inline function _regular_pullback_groups!(subject_values_bar, all_params_bar,
    subject_values, groups::Tuple, cells, tfs, supports, scratch, g::Int)
    _regular_pullback_group!(subject_values_bar, all_params_bar, subject_values,
        groups[1], cells[g], tfs[g], supports, scratch)
    return _regular_pullback_groups!(subject_values_bar, all_params_bar,
        subject_values, Base.tail(groups), cells, tfs, supports, scratch, g + 1)
end

function _regular_pullback_group!(subject_values_bar::AbstractVector,
    all_params_bar::AbstractVector, subject_values::AbstractVector,
    group::Vector{F}, cells::Vector{Int}, tfs::Vector{Int},
    supports::AbstractVector{Vector{Int}}, scratch::AbstractVector) where {F}
    @inbounds for k in eachindex(group)
        cotangent = all_params_bar[cells[k]]
        iszero(cotangent) && continue
        transform = group[k]
        support = supports[tfs[k]]
        if length(support) > 1
            for pn in support
                scratch[pn] = _seed_dual(scratch[pn], subject_values[pn], false)
            end
        end
        for pn in support
            base = subject_values[pn]
            scratch[pn] = _seed_dual(scratch[pn], base, true)
            derivative = _partial1(transform(scratch))
            scratch[pn] = _seed_dual(scratch[pn], base, false)
            subject_values_bar[pn] += cotangent * derivative
        end
    end
    return nothing
end

@inline _seed_dual(::D, value, seed::Bool) where {D<:ForwardDiff.Dual} =
    D(value, ForwardDiff.Partials((seed ? one(value) : zero(value),)))

"""
    _ctsem_ti_pullback!(values_bar, subject_values_bar, sp, tipreds, values)

Undo `_materialize_subject_values!`: the identity copy plus each TI predictor
effect `subject_values[par] += values[coef] * tipred[pred]`. A fixed effect
(coefficient index 0) has no coefficient to receive a term.

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
    them rather than overwriting, as it must. The second term uses the effect,
    `_ctsem_ti_effect`, so it holds for a fixed one too; each predictor value
    is read from the recipe rather than from a row built for the purpose.
"""
function _ctsem_ti_pullback!(values_bar::AbstractVector,
    subject_values_bar::AbstractVector, sp::EKFParameters, tipreds::AbstractVector,
    values::AbstractVector)
    @inbounds for i in eachindex(values_bar)
        values_bar[i] += subject_values_bar[i]
    end
    @inbounds for i in eachindex(sp.ti_parameter_indices)
        coefficient = sp.ti_coefficient_indices[i]
        coefficient == 0 && continue
        values_bar[coefficient] += subject_values_bar[sp.ti_parameter_indices[i]] *
            tipreds[sp.ti_predictor_indices[i]]
    end
    return values_bar
end

function _ctsem_ti_pullback!(values_bar::AbstractVector,
    subject_values_bar::AbstractVector, sp::EKFParameters, spec::TIMissingRecipe,
    values::AbstractVector)
    @inbounds for i in eachindex(values_bar)
        values_bar[i] += subject_values_bar[i]
    end
    @inbounds for i in eachindex(sp.ti_parameter_indices)
        predictor = sp.ti_predictor_indices[i]
        coefficient = sp.ti_coefficient_indices[i]
        cotangent = subject_values_bar[sp.ti_parameter_indices[i]]
        # Present only when this effect's predictor is one of this subject's
        # sampled cells, whose value is then that raw parameter.
        if coefficient != 0
            values_bar[coefficient] += cotangent * _ctsem_tipred_value(spec, predictor, values)
        end
        k = findfirst(==(predictor), spec.predictor_index)
        if k !== nothing
            values_bar[spec.parameter_index[k]] += cotangent * _ctsem_ti_effect(sp, values, i)
        end
    end
    # A T0MEANS or T0VAR cell that reads a sampled predictor reads it from the
    # tail (`EKFParameters.ti_tail`), and its cotangent lands there.
    if !isempty(sp.ti_tail)
        @inbounds for k in eachindex(spec.predictor_index)
            values_bar[spec.parameter_index[k]] +=
                subject_values_bar[first(sp.ti_tail) - 1 + spec.predictor_index[k]]
        end
    end
    return values_bar
end

"""
    _ctsem_tipred_bar_reset!(aws, tipreds, n)

Before a subject's reverse pass: its TI predictor cotangent zeroed, and which
of its `n` predictors are sampled -- a `TIMissingRecipe`'s cells, the only ones
whose cotangent `_group_reverse_walk!` computes. A plain row has none.
"""
function _ctsem_tipred_bar_reset!(aws, tipreds, n::Integer)
    resize!(aws.tipreds_bar, n)
    fill!(aws.tipreds_bar, zero(eltype(aws.tipreds_bar)))
    resize!(aws.tipred_sampled, n)
    fill!(aws.tipred_sampled, false)
    if tipreds isa TIMissingRecipe
        @inbounds for j in tipreds.predictor_index
            aws.tipred_sampled[j] = true
        end
    end
    return nothing
end

"""
    _ctsem_tipred_scatter!(gradient, aws, tipreds)

After a subject's reverse pass: what cells reading a sampled TI predictor sent
back to it (`aws.tipreds_bar`), onto that cell's raw parameter. Per subject and
straight onto raw positions, as `_hybrid_scatter!` does, so a parameter layer
unwound once for all subjects cannot drop it. An effect's route to the same
parameter is `_ctsem_ti_pullback!`'s.
"""
_ctsem_tipred_scatter!(gradient, aws, tipreds) = gradient
function _ctsem_tipred_scatter!(gradient, aws, tipreds::TIMissingRecipe)
    @inbounds for k in eachindex(tipreds.predictor_index)
        gradient[tipreds.parameter_index[k]] += aws.tipreds_bar[tipreds.predictor_index[k]]
    end
    return gradient
end

################################################################################
# Layer 3 pullback
################################################################################

"""
    _ctsem_complex_group_pullback!(all_params_bar, state_bar, transforms, indices,
                                   supports, ctx, relevant, dual_ctx, ti)

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
     already cleared the post-write cotangent), and to every sampled TI
     predictor it reads (`ti`, see `_group_reverse_walk!`).

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

Both walks recurse over the transform tuple rather than index it, so each
transform is called at its own concrete type with a concretely typed context.
Indexed with a runtime `k`, the element is a union of the group's closure
types, which the compiler splits only up to four; past that, or with an
untyped context, every call boxed its context and its result -- on a
one-effect augmented model 90% of what a gradient allocated, per row.
"""
function _ctsem_complex_group_pullback!(all_params_bar::AbstractVector,
    state_bar::AbstractVector, transforms::Union{Tuple,AbstractVector},
    indices::AbstractVector{Int},
    supports::AbstractVector{Vector{Int}}, ctx::CTSEMRowContext,
    relevant::AbstractVector{Int}, dual_ctx, ti)
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
    _group_replay!(written, working, indices, transforms, ctx, 1)

    # Start the reverse walk at transform K, whose inputs are the pre-group
    # values with every earlier transform's write applied. Only the relevant
    # entries are synced: a transform cannot read anything else by
    # construction, and syncing the whole vector per group per row is what
    # this compaction exists to avoid.
    # `working` holds *post*-replay values now, so every written cell has to be
    # set explicitly rather than relying on the sync: cells 1..K-1 to what
    # their transforms wrote, and cell K back to its pre-group value, since
    # transform K ran before its own write landed.
    dctx = _sync_dual_context!(dual_ctx, ctx, working, relevant)
    @inbounds for k in 1:(nk - 1)
        _set_dual_value!(dual_ctx, indices[k], written[k])
    end
    @inbounds _set_dual_value!(dual_ctx, indices[nk], saved[nk])

    _group_reverse_walk!(all_params_bar, state_bar, transforms, indices, supports,
        saved, dual_ctx, dctx, ti, nk)
    return nothing
end

@inline _group_replay!(written, working, indices, ::Tuple{}, ctx, k::Int) = nothing

function _group_replay!(written, working, indices, transforms::Tuple, ctx, k::Int)
    value = first(transforms)(ctx)
    @inbounds written[k] = value
    @inbounds working[indices[k]] = value
    return _group_replay!(written, working, indices, Base.tail(transforms), ctx, k + 1)
end

# A group of `CTSEMExpr`s is one type, so it is walked by index.
function _group_replay!(written, working, indices, transforms::Vector{CTSEMExpr}, ctx,
    k::Int)
    @inbounds for j in k:length(transforms)
        value = transforms[j](ctx)
        written[j] = value
        working[indices[j]] = value
    end
    return nothing
end

@inline _group_reverse_walk!(all_params_bar, state_bar, ::Tuple{}, indices, supports,
    saved, dual_ctx, dctx, ti, k::Int) = nothing

# Transform `k` is the tuple's last element; the walk then steps back to the
# tuple without it. `ti` is the group's predictor reads (`ti.supports`), the
# subject's predictor cotangent (`ti.bar`), and which of its predictors are
# sampled (`ti.sampled`): a read of a sampled one sends the cell's cotangent
# back to it, where `_ctsem_tipred_scatter!` collects it for the parameter.
function _group_reverse_walk!(all_params_bar, state_bar, transforms::Tuple, indices,
    supports, saved, dual_ctx, dctx, ti, k::Int)
    transform = last(transforms)
    @inbounds begin
        idx = indices[k]
        cotangent = all_params_bar[idx]
        all_params_bar[idx] = zero(cotangent)
        if !iszero(cotangent)
            for j in supports[k]
                all_params_bar[j] += cotangent *
                    _dual_param_derivative(dual_ctx, transform, dctx, j)
            end
            for s in eachindex(state_bar)
                state_bar[s] += cotangent *
                    _dual_state_derivative(dual_ctx, transform, dctx, s)
            end
            for j in ti.supports[k]
                ti.sampled[j] || continue
                ti.bar[j] += cotangent *
                    _dual_tipred_derivative(dual_ctx, transform, dctx, j)
            end
        end
        # Step back one transform: undo transform k-1's write so that
        # transform k-1 sees the inputs it actually had.
        k > 1 && _set_dual_value!(dual_ctx, indices[k - 1], saved[k - 1])
    end
    return _group_reverse_walk!(all_params_bar, state_bar, Base.front(transforms),
        indices, supports, saved, dual_ctx, dctx, ti, k - 1)
end

# An expression knows the states it reads, so it is differentiated against
# those alone: a random effect's `10*state[k]` once rather than once per state.
function _group_reverse_walk!(all_params_bar, state_bar, transforms::Vector{CTSEMExpr},
    indices, supports, saved, dual_ctx, dctx, ti, k::Int)
    @inbounds for m in k:-1:1
        transform = transforms[m]
        idx = indices[m]
        cotangent = all_params_bar[idx]
        all_params_bar[idx] = zero(cotangent)
        if !iszero(cotangent)
            if transform.form === _XFORM_SCALED_STATE
                # d(c x)/dx = c, exactly what the seeded dual returns.
                state_bar[transform.fa] += cotangent * transform.fc
            else
                V = eltype(state_bar)
                n = length(transform.op)
                if n <= 8
                    _ctsem_expr_pullback!(all_params_bar, state_bar, transform, dctx, ti,
                        cotangent, V, ntuple(_ -> zero(V), Val(8)),
                        Base.setindex(ntuple(_ -> zero(V), Val(8)), one(V), transform.root))
                else
                    sc = _ctsem_expr_scratch(V)
                    g = _ctsem_expr_buffer(sc.b, n)
                    fill!(view(g, 1:n), zero(V))
                    g[transform.root] = one(V)
                    _ctsem_expr_pullback!(all_params_bar, state_bar, transform, dctx, ti,
                        cotangent, V, _ctsem_expr_buffer(sc.a, n), g)
                end
            end
        end
        m > 1 && _set_dual_value!(dual_ctx, indices[m - 1], saved[m - 1])
    end
    return nothing
end

"""
    _ctsem_expr_pullback!(all_params_bar, state_bar, e, dctx, ti, cotangent, V, x, g)

Every partial of expression `e` at the values `dctx` holds, by one pass forward
over its nodes and one back, in the record's scalar type `V`, each added times
`cotangent` to the parameter, state and sampled TI predictor it belongs to.
What the seeded duals (`_dual_param_derivative` and its siblings) compute, at
one evaluation per input each; the two agree to rounding, not bit for bit,
since a chain of products is multiplied in the other order. `x` holds the
node values and `g` their adjoints, tuples or vectors at least as long as the
expression; `g` arrives zero but for one at the root.
"""
function _ctsem_expr_pullback!(all_params_bar, state_bar, e::CTSEMExpr, dctx, ti,
    cotangent, ::Type{V}, x, g) where {V}
    n = length(e.op)
    @inbounds for i in 1:n
        op = e.op[i]
        v = if op === _XC
            convert(V, e.val[i])
        elseif op === _XP
            convert(V, ForwardDiff.value(dctx.pars[e.a[i]]))
        elseif op === _XS
            convert(V, ForwardDiff.value(dctx.state[e.a[i]]))
        elseif op === _XTD
            convert(V, ForwardDiff.value(dctx.tdpreds[e.a[i]]))
        elseif op === _XTI
            convert(V, ForwardDiff.value(dctx.tipreds[e.a[i]]))
        elseif op >= _XADD
            convert(V, _ctsem_expr_binary(op, x[e.a[i]], x[e.b[i]]))
        elseif op === _XPOWI
            convert(V, _ctsem_expr_powi(x[e.a[i]], e.b[i]))
        else
            convert(V, _ctsem_expr_unary(op, x[e.a[i]]))
        end
        x = _ctsem_expr_put(x, v, i)
    end
    @inbounds for i in n:-1:1
        gi = g[i]
        iszero(gi) && continue
        op = e.op[i]
        op === _XC && continue
        if op === _XP
            all_params_bar[e.a[i]] += cotangent * gi
        elseif op === _XS
            state_bar[e.a[i]] += cotangent * gi
        elseif op === _XTI
            j = e.a[i]
            ti.sampled[j] && (ti.bar[j] += cotangent * gi)
        elseif op === _XTD
            continue
        elseif op >= _XADD
            a = e.a[i]
            b = e.b[i]
            xa = x[a]
            xb = x[b]
            da, db = _ctsem_expr_binary_partials(op, xa, xb, x[i])
            g = _ctsem_expr_put(g, g[a] + gi * da, a)
            g = _ctsem_expr_put(g, g[b] + gi * db, b)
        else
            a = e.a[i]
            g = _ctsem_expr_put(g, g[a] + gi * _ctsem_expr_unary_partial(op, x[a], x[i], e.b[i]), a)
        end
    end
    return nothing
end

# Partials of a binary node with operands `a`, `b` and value `y`.
@inline function _ctsem_expr_binary_partials(op::UInt8, a, b, y)
    op === _XMUL && return (b, a)
    op === _XADD && return (one(a), one(b))
    op === _XSUB && return (one(a), -one(b))
    op === _XDIV && return (inv(b), -y / b)
    op === _XPOW && return (b * a^(b - one(b)), iszero(a) ? zero(a) : y * log(a))
    # min and max route to the operand they returned, as Base picks it.
    op === _XMIN && return isless(b, a) ? (zero(a), one(b)) : (one(a), zero(b))
    return isless(b, a) ? (one(a), zero(b)) : (zero(a), one(b))
end

# Partial of a unary node with operand `a` and value `y`; `n` is `_XPOWI`'s power.
@inline function _ctsem_expr_unary_partial(op::UInt8, a, y, n::Int)
    op === _XPOWI && return n * a^(n - 1)
    op === _XNEG && return -one(a)
    op === _XEXP && return y
    op === _XLOG && return inv(a)
    op === _XSQRT && return inv(2 * y)
    op === _XABS && return signbit(a) ? -one(a) : one(a)
    op === _XLOG1PEXP && return inv_logit(a)
    op === _XINVLOGIT && return y * (1 - y)
    op === _XLOGIT && return inv(a * (1 - a))
    op === _XTANH && return 1 - y^2
    op === _XLOG1P && return inv(1 + a)
    op === _XEXPM1 && return y + 1
    op === _XSIN && return cos(a)
    op === _XCOS && return -sin(a)
    return -inv(a^2)
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
    # The row's TD predictors lifted into `D`, when they arrive as a dual `T`;
    # see `_dual_context_data`.
    tdpreds::Vector{D}
    # The subject's TI predictors, lifted the same way: a cell can read them.
    tipreds::Vector{D}
end

function CTSEMDualContext(::Type{T}, sp::EKFParameters, state_like) where {T}
    D = ForwardDiff.Dual{Nothing,T,1}
    zero_dual = D(zero(T), ForwardDiff.Partials((zero(T),)))
    data = fill(zero_dual, length(sp.mutables))
    pars = CTSEMPars(data, sp.parameter_axis)
    state = similar(state_like, D)
    fill!(state, zero_dual)
    return CTSEMDualContext{D,typeof(pars),typeof(state),T}(data, pars, state,
        T[], T[], D[], D[])
end

# The row's data -- TD predictors, time, dt -- reach the reverse pass in its
# scalar `T`, which inside `ctsem_hessian` is the outer `ForwardDiff.Dual`.
# The mirror's parameters are `Dual{Nothing,T,1}`, and ForwardDiff cannot
# order an untagged dual against a tagged one, so an expression combining a
# predictor with a parameter threw `DualMismatchError` -- and one reading only
# a predictor returned an outer dual that `_partial1` took for the inner one.
# Lifted into `D` with a zero inner partial, they keep their outer derivatives
# and every operand shares one type. Plain floats mix with any dual and pass
# through, so the ordinary gradient copies nothing.
@inline _dual_context_data(::Type{D}, x::AbstractFloat) where {D} = x
@inline _dual_context_data(::Type{D}, x::Real) where {D} =
    D(x, ForwardDiff.Partials((zero(x),)))
_dual_context_data(::Type{D}, x::AbstractVector{<:AbstractFloat}, buffer) where {D} = x
_dual_context_data(::Type{D}, x, buffer) where {D} = x
function _dual_context_data(::Type{D}, x::AbstractVector{<:Real}, buffer) where {D}
    resize!(buffer, length(x))
    @inbounds for (i, v) in enumerate(x)
        buffer[i] = _dual_context_data(D, v)
    end
    return buffer
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
and return the dual row context the seeded derivative calls below evaluate the
transforms at: this row's time, dt and predictors over the mirror's state and
parameters. Returned rather than kept on `dual_ctx`, where its type -- which
depends on the row's predictors -- could only be `Any`, and every transform
called with it was a dynamic dispatch boxing it.
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
    D = eltype(dual_ctx.data)
    return CTSEMRowContext(dual_ctx.state, dual_ctx.pars,
        _dual_context_data(D, ctx.tdpreds, dual_ctx.tdpreds),
        _dual_context_data(D, ctx.tipreds, dual_ctx.tipreds),
        _dual_context_data(D, ctx.time), _dual_context_data(D, ctx.dt),
        ctx.subject, ctx.row)
end

# A transform that happens not to depend on any dual input returns a plain
# real; treat that as a zero partial rather than letting `partials` throw.
# Only the untagged dual is this layer's seed: a tagged one is an outer
# derivative (`ctsem_hessian`'s) and its first partial is not the one asked for.
@inline _partial1(x::ForwardDiff.Dual{Nothing}) = ForwardDiff.partials(x, 1)
@inline _partial1(x::Real) = zero(x)

# The first partial of a transform's result, as the seed's value type `V`: a
# transform reading only data returns a plain real, whose partial is zero.
@inline function _partial1_as(r, ::Type{ForwardDiff.Dual{Nothing,V,N}}) where {V,N}
    r isa ForwardDiff.Dual{Nothing,V,N} && return ForwardDiff.partials(r, 1)::V
    return convert(V, _partial1(r))::V
end

# Specialised on the transform and on `dctx`, both concrete where the group
# walk calls these, so the transform call is static and nothing is boxed.
@inline function _dual_param_derivative(dual_ctx::CTSEMDualContext{D}, transform::F,
    dctx::CTSEMRowContext, j::Int) where {D,F}
    saved = dual_ctx.data[j]
    dual_ctx.data[j] = _seed_dual(saved, ForwardDiff.value(saved), true)
    derivative = _partial1_as(transform(dctx), D)
    dual_ctx.data[j] = saved
    return derivative
end

@inline function _dual_state_derivative(dual_ctx::CTSEMDualContext{D}, transform::F,
    dctx::CTSEMRowContext, s::Int) where {D,F}
    saved = dual_ctx.state[s]
    dual_ctx.state[s] = _seed_dual(saved, ForwardDiff.value(saved), true)
    derivative = _partial1_as(transform(dctx), D)
    dual_ctx.state[s] = saved
    return derivative
end

# The derivative with respect to TI predictor `j`, which the row context
# carries as data -- plain floats on an ordinary gradient, passed through
# unlifted by `_sync_dual_context!` -- so the row is lifted into the mirror's
# predictor buffer here, on the one path that needs it: a cell reading a
# predictor whose value this subject samples. Under `ctsem_hessian` the row
# arrives lifted already, as that same buffer, and is left as it is.
@inline function _dual_tipred_derivative(dual_ctx::CTSEMDualContext{D}, transform::F,
    dctx::CTSEMRowContext, j::Int) where {D,F}
    lifted = _dual_lift!(dual_ctx.tipreds, dctx.tipreds)
    saved = lifted[j]
    lifted[j] = _seed_dual(saved, ForwardDiff.value(saved), true)
    seeded = CTSEMRowContext(dctx.state, dctx.pars, dctx.tdpreds, lifted, dctx.time,
        dctx.dt, dctx.subject, dctx.row)
    derivative = _partial1_as(transform(seeded), D)
    lifted[j] = saved
    return derivative
end

function _dual_lift!(buffer::Vector{D}, x::AbstractVector) where {D<:ForwardDiff.Dual}
    x === buffer && return buffer
    resize!(buffer, length(x))
    @inbounds for i in eachindex(x)
        buffer[i] = _dual_lift(D, x[i])
    end
    return buffer
end
@inline _dual_lift(::Type{D}, v::D) where {D<:ForwardDiff.Dual} = v
@inline _dual_lift(::Type{D}, v::Real) where {D<:ForwardDiff.Dual} =
    D(v, ForwardDiff.Partials((zero(v),)))
