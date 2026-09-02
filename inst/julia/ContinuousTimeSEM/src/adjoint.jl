"""
Adjoint-gradient integration boundary.

The hand-derived reverse-mode mathematics lives one level down, in
`adjoint_primitives.jl` (generic matrix kernels), `adjoint_parameters.jl` (the
transform layer) and `adjoint_ekf.jl` (the reverse filter). This file holds the
public entry points: the workspace that caches everything discovered once per
model, `ctsem_adjoint_gradient`, the gating status, and the validation oracle
that cross-checks the adjoint against ForwardDiff and finite differences.
"""

using LinearAlgebra
using ForwardDiff

export ctsem_adjoint_status, ctsem_validate_forward_gradient, ctsem_adjoint_gradient

"""
    CTSEMAdjointWorkspace

Per-model state for the reverse pass, built once and reused across evaluations.

Holds the things that are discovered by inspecting the model rather than by
running it: which `values` index each regular transform reads, which parameter
cells each state-dependent transform reads, and the dual-number scratch those
lookups are evaluated through. `tipreds` is the one genuinely per-subject
field, refreshed by `ctsem_adjoint_gradient` before each subject's reverse
pass.
"""
mutable struct CTSEMAdjointWorkspace{T,SP,LB}
    sp::SP
    n::Int
    m::Int
    diffusion_state_indices::Vector{Int}
    predict_indices::Vector{Int}
    update_indices::Vector{Int}
    td_indices::Vector{Int}
    predict_supports::Vector{Vector{Int}}
    update_supports::Vector{Vector{Int}}
    td_supports::Vector{Vector{Int}}
    regular_supports::Vector{Int}
    regular_dual_scratch::Vector{ForwardDiff.Dual{Nothing,T,1}}
    dual_context::Any
    tipreds::Vector{T}
    # Persistent working copy of `all_params` used to replay a transform group
    # forward during the reverse pass; see `_ctsem_complex_group_pullback!`.
    # Initialised to `NaN`: only cells in a group's relevant set are ever
    # written, so anything a transform reads that was not recorded stays `NaN`
    # and poisons the gradient visibly instead of quietly using a stale value.
    group_scratch::Vector{T}
    # Reused across every reverse prediction substep. Sized to the *dynamic*
    # state block, which is what the Lyapunov solve actually operates on.
    lyap_buffer::LB
    # Reverse-pass accumulators, reused across subjects rather than allocated
    # per subject. The reverse pass wraps `theta_bar` in a ComponentVector
    # view locally rather than caching one here -- see the note in
    # `_ctsem_reverse_tape!` for why that matters.
    theta_bar::Vector{T}
    x_bar::Vector{T}
    P_bar::Matrix{T}
    manifest_cov_bar::Matrix{T}
    subject_values_bar::Vector{T}
    # True when the parameter layer is identical for every subject, so it can
    # be unwound once at the end instead of once per subject.
    parameter_layer_shareable::Bool
    # --- deferred matrix-exponential Frechet derivative (see `_flush_frechet!`)
    # `A = exp(JAx * dt)` is the single most expensive step in the reverse pass,
    # and `L(A, E)` is linear in `E`, so directions belonging to the same
    # `JAx * dt` are accumulated and pushed through one block exponential
    # instead of one each.
    frechet_pending::Bool
    frechet_dt::T
    frechet_A::Matrix{T}
    frechet_accum::Matrix{T}
    # Frechet contributions held back so the batch can span subjects too;
    # `ctsem_adjoint_gradient` pushes these through the parameter layer once.
    jax_bar_deferred::Matrix{T}
    # Whether that cross-subject deferral is legal for this model.
    defer_frechet::Bool
    # Whether any state-dependent transform group writes a JAx cell. Such a
    # group's reverse *zeroes* that cell's cotangent, so it is the one thing
    # that can consume an outstanding Frechet contribution mid-tape.
    groups_write_jax::Bool
    # Flat `all_params` positions belonging to the JAx component.
    jax_positions::Vector{Int}
    # Working storage for the reverse pass; see `CTSEMReverseScratch` for why
    # the temporaries come from here rather than from the heap.
    reverse_scratch::CTSEMReverseScratch{T}
    tape::CTSEMAdjointTape{T}
end

function CTSEMAdjointWorkspace(::Type{T}, sp::EKFParameters, nvalues::Integer,
    ntdpred::Integer) where {T}
    ws = _init_continuous_ekf_workspace(T, sp)
    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)

    predict_indices = findall(sp.predict_transforms_indices)
    update_indices = findall(sp.update_transforms_indices)
    td_indices = findall(sp.td_transforms_indices)

    zero_dual = ForwardDiff.Dual{Nothing,T,1}(zero(T), ForwardDiff.Partials((zero(T),)))
    theta_bar = zeros(T, length(sp.mutables))
    lyap_buffer = LyapBuffer(T, length(ws.diffusion_state_indices))

    predict_supports = _ctsem_complex_transform_supports(sp.predict_transforms, predict_indices, sp, n, Int(ntdpred))
    update_supports = _ctsem_complex_transform_supports(sp.update_transforms, update_indices, sp, n, Int(ntdpred))
    td_supports = _ctsem_complex_transform_supports(sp.td_transforms, td_indices, sp, n, Int(ntdpred))
    group_relevant = [
        _ctsem_group_relevant(predict_supports, predict_indices),
        _ctsem_group_relevant(td_supports, td_indices),
        _ctsem_group_relevant(update_supports, update_indices),
    ]
    group_scratch = fill(T(NaN), length(sp.mutables))

    jax_positions = _ctsem_jax_positions(sp)
    jax_set = Set(jax_positions)
    groups_write_jax = any(i -> i in jax_set,
        Iterators.flatten((predict_indices, td_indices, update_indices)))
    defer_frechet = isempty(sp.ti_parameter_indices) && !groups_write_jax

    return CTSEMAdjointWorkspace{T,typeof(sp),typeof(lyap_buffer)}(
        sp, n, m,
        collect(ws.diffusion_state_indices),
        predict_indices, update_indices, td_indices,
        predict_supports,
        update_supports,
        td_supports,
        _ctsem_regular_transform_supports(sp, nvalues),
        fill(zero_dual, Int(nvalues)),
        CTSEMDualContext(T, sp, zeros(T, n)),
        T[],
        group_scratch,
        lyap_buffer,
        theta_bar,
        zeros(T, n),
        zeros(T, n, n),
        zeros(T, m, m),
        Vector{T}(undef, Int(nvalues)),
        _ctsem_parameter_layer_shareable(sp, predict_indices, update_indices, td_indices),
        false, zero(T), zeros(T, n, n), zeros(T, n, n), zeros(T, n, n),
        defer_frechet, groups_write_jax, jax_positions,
        CTSEMReverseScratch(T, n, m, length(ws.diffusion_state_indices)),
        CTSEMAdjointTape(T, group_relevant),
    )
end

"""
    _ctsem_jax_positions(sp)

The flat `all_params` positions belonging to the `JAx` component.

Discovered by marking the component through a `ComponentVector` view rather
than by arithmetic on the axis, so it stays correct whatever order the R-side
parameter table happens to list matrices in.
"""
function _ctsem_jax_positions(sp::EKFParameters)
    marker = ComponentVector(zeros(Int, length(sp.mutables)), sp.parameter_axis)
    marker.JAx .= 1
    return findall(!iszero, getdata(marker))
end

"""
    _ctsem_group_relevant(supports, indices)

The `all_params` indices one transform group can read or write: the union of
its transforms' recorded read sets with the cells they write.

This is exactly what a group record has to snapshot. Everything else in
`all_params` is unreachable from that group's expressions, so recording it
would be dead weight -- ~23 KB per group per row on a 20-latent model.
"""
function _ctsem_group_relevant(supports::AbstractVector{Vector{Int}},
    indices::AbstractVector{Int})
    isempty(indices) && return Int[]
    relevant = Set{Int}(indices)
    for support in supports
        union!(relevant, support)
    end
    return sort!(collect(relevant))
end

"""
    _ctsem_parameter_layer_shareable(sp, predict_indices, update_indices, td_indices)

Whether every subject sees the same parameter layer, so the cotangent on
`all_params` can be accumulated across subjects and unwound once.

Two things break it:

  * **TI predictor effects.** `_materialize_subject_values!` shifts a
    parameter by `values[coef] * tipred[pred]`, so `subject_values` -- and
    therefore each regular transform's derivative -- differs per subject.
  * **State-dependent transforms.** Their reverse *overwrites* (zeroes) the
    cotangent on the cells they write, because the forward pass overwrote the
    values. Running subject 2's reverse over a buffer still holding subject
    1's accumulated cotangent would destroy it. This is a genuine correctness
    constraint, not a conservatism.

When neither applies -- an ordinary linear multi-subject model, which is the
common ctsem shape -- the parameter layer is `O(number of parameter cells)`
work that would otherwise be repeated once per subject for no reason. On a
20-latent, 20-subject model that is the difference between doing it 20 times
and once.
"""
function _ctsem_parameter_layer_shareable(sp::EKFParameters, predict_indices,
    update_indices, td_indices)
    isempty(sp.ti_parameter_indices) || return false
    return isempty(predict_indices) && isempty(update_indices) && isempty(td_indices)
end

"""
    ctsem_adjoint_gradient(objective, values)

Return `(value = ..., gradient = ...)` for a prepared `CTSEMObjective`, using
the reverse-mode adjoint rather than ForwardDiff.

Each subject is filtered forward once with tracing on, then its tape is
replayed backwards and its contribution accumulated into the shared gradient.
The cost is therefore independent of the number of free parameters, which is
the whole point: `ForwardDiff` needs one dual pass per chunk of parameters, so
its cost grows with the parameter count that large `indvarying` models
inflate.

An invalid trial point -- one where the primal returns a non-finite
log-likelihood, e.g. because a Cholesky factorization failed -- propagates as a
non-finite value *and* a non-finite gradient, matching what the ForwardDiff
path does and what `ctsem_optimize`'s `fg!` guards expect. There is deliberately
no silent fallback to ForwardDiff.
"""
function ctsem_adjoint_gradient(objective::CTSEMObjective, values::AbstractVector{T}) where {T}
    subjects = objective.subject_objectives
    nsubjects = length(subjects)
    nchunks = _ctsem_nchunks(nsubjects)
    workspaces = _get_or_init_adjoint_workspaces!(objective, T, length(values), nchunks)
    ranges = _ctsem_chunk_ranges(nsubjects, nchunks)

    # Per-chunk accumulators, summed at the end. Each chunk is independently a
    # valid sub-objective: the parameter layer and the deferred Frechet
    # contribution are both unwound within the chunk that produced them, so
    # summing chunk gradients is exactly summing subject gradients, only in a
    # different order.
    gradients = [zeros(T, length(values)) for _ in 1:nchunks]
    totals = zeros(T, nchunks)
    valid = fill(true, nchunks)
    badvalue = fill(T(NaN), nchunks)

    if nchunks <= 1
        _ctsem_adjoint_chunk!(gradients[1], totals, valid, badvalue, 1,
            ranges[1], subjects, objective.params, workspaces[1], values)
    else
        Threads.@sync for c in 1:nchunks
            Threads.@spawn _ctsem_adjoint_chunk!(gradients[c], totals, valid,
                badvalue, c, ranges[c], subjects, objective.params,
                workspaces[c], values)
        end
    end

    @inbounds for c in 1:nchunks
        if !valid[c]
            # One invalid subject invalidates the whole evaluation, exactly as
            # in the serial path: a non-finite value *and* a non-finite
            # gradient, so an optimizer's finiteness guards reject the trial
            # point rather than accepting a partial gradient.
            return (value=badvalue[c], gradient=fill(T(NaN), length(values)))
        end
    end

    total = zero(T)
    gradient = gradients[1]
    @inbounds for c in 1:nchunks
        total += totals[c]
        if c > 1
            gradient .+= gradients[c]
        end
    end
    # The prior is a closed-form function of the raw parameters alone, so it is
    # added once here rather than inside a chunk.
    total += _ctsem_log_prior(objective, values)
    _ctsem_log_prior_gradient!(gradient, objective, values)
    return (value=total, gradient=gradient)
end

export ctsem_subject_gradients
"""
    ctsem_subject_gradients(objective, values)

Return `(value = ..., scores = ...)` where `scores[i, :]` is subject `i`'s own
gradient contribution.

This is the score matrix the R side's `scorecalc()` produces for the Stan
backend, and it is what the OPG, sandwich and score-bootstrap uncertainty
methods consume. The adjoint already computes exactly this and then adds it up,
so producing it costs one ordinary gradient evaluation rather than the
per-subject model re-initialisation the Stan path needs.

Two shortcuts of the summed gradient have to be switched off for the rows to be
individually correct, and both are per-subject flags rather than changes to the
reverse pass:

  * the shared parameter layer, which unwinds the transform layer once for a
    whole chunk of subjects, is unwound per subject instead;
  * the deferred matrix-exponential Frechet contribution, which is normally
    batched across subjects, is flushed within each subject.

`sum(scores, dims=1)` therefore equals `ctsem_adjoint_gradient`'s gradient, and
`test_subject_gradients.jl` asserts that -- which is the natural check, since
the two routes share every primitive but differ in where they accumulate.
"""
function ctsem_subject_gradients(objective::CTSEMObjective, values::AbstractVector{T}) where {T}
    subjects = objective.subject_objectives
    nsubjects = length(subjects)
    npars = length(values)
    scores = zeros(T, nsubjects, npars)
    nchunks = _ctsem_nchunks(nsubjects)
    workspaces = _get_or_init_adjoint_workspaces!(objective, T, npars, nchunks)
    ranges = _ctsem_chunk_ranges(nsubjects, nchunks)
    totals = zeros(T, nchunks)
    valid = fill(true, nchunks)
    badvalue = fill(T(NaN), nchunks)

    if nchunks <= 1
        _ctsem_subject_gradient_chunk!(scores, totals, valid, badvalue, 1,
            ranges[1], subjects, workspaces[1], values)
    else
        Threads.@sync for c in 1:nchunks
            Threads.@spawn _ctsem_subject_gradient_chunk!(scores, totals, valid,
                badvalue, c, ranges[c], subjects, workspaces[c], values)
        end
    end

    @inbounds for c in 1:nchunks
        valid[c] || return (value=badvalue[c], scores=fill(T(NaN), nsubjects, npars))
    end

    # Each subject's row carries 1/nsubjects of the prior, which is what the R
    # side's `scorecalc()` does for the Stan backend (it sets
    # `standata$priormod = 1/nsubjects` before taking per-subject gradients).
    # The rows then still sum to the full posterior gradient.
    if !isempty(objective.prior_index) && nsubjects > 0
        share = 1 / nsubjects
        @inbounds for i in 1:nsubjects
            _ctsem_log_prior_gradient!(view(scores, i, :), objective, values, share)
        end
    end
    return (value=sum(totals) + _ctsem_log_prior(objective, values), scores=scores)
end

"""Per-subject gradients for one contiguous chunk, written into `scores`."""
function _ctsem_subject_gradient_chunk!(scores::Matrix{T}, totals::Vector{T},
    valid::Vector{Bool}, badvalue::Vector{T}, c::Int, range::UnitRange{Int},
    subjects, aws, values::AbstractVector{T}) where {T}

    deferred = aws.defer_frechet
    aws.defer_frechet = false
    total = zero(T)
    try
        @inbounds for i in range
            subject_objective = subjects[i]
            ws = _get_or_init_objective_workspace!(subject_objective, T)
            tape = _tape_reset!(aws.tape)
            resize!(aws.tipreds, length(subject_objective.tipreds))
            copyto!(aws.tipreds, subject_objective.tipreds)
            aws.frechet_pending = false

            loglik = _extended_kalman_filter_continuous!(ws, values,
                subject_objective.data, subject_objective.timesteps,
                subject_objective.params, subject_objective.tdpreds,
                subject_objective.tipreds, subject_objective.subject,
                subject_objective.max_timestep, tape)
            # _finite_deep, not isfinite: under ctsem_hessian this runs at Dual
            # and isfinite tests the value alone, so a NaN partial would pass.
            if !_finite_deep(loglik)
                valid[c] = false
                badvalue[c] = loglik
                return nothing
            end
            total += loglik

            fill!(aws.theta_bar, zero(T))
            _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
            _ctsem_parameter_layer!(view(scores, i, :), aws.theta_bar,
                tape.subject_values, subject_objective.params, aws,
                subject_objective.tipreds)
        end
    finally
        aws.defer_frechet = deferred
    end
    totals[c] = total
    return nothing
end

"""
    _ctsem_adjoint_chunk!(gradient, totals, valid, badvalue, c, range, subjects,
                          params, aws, values)

Accumulate one contiguous chunk of subjects into its own gradient and workspace.

This is the serial adjoint loop, scoped to `range`. Everything it touches --
`aws` and its tape, `gradient`, and slot `c` of the shared result vectors -- is
private to this chunk, so no synchronisation is needed beyond the enclosing
`@sync`.
"""
function _ctsem_adjoint_chunk!(gradient::Vector{T}, totals::Vector{T},
    valid::Vector{Bool}, badvalue::Vector{T}, c::Int, range::UnitRange{Int},
    subjects, params, aws, values::AbstractVector{T}) where {T}

    shared = aws.parameter_layer_shareable
    shared && fill!(aws.theta_bar, zero(T))
    # A previous call that bailed out on an invalid trial point can leave a
    # queued Frechet direction behind; it belongs to that call's cotangent.
    aws.frechet_pending = false
    fill!(aws.jax_bar_deferred, zero(T))
    last_subject_values = nothing
    last_tipreds = nothing
    total = zero(T)

    @inbounds for i in range
        subject_objective = subjects[i]
        ws = _get_or_init_objective_workspace!(subject_objective, T)
        tape = _tape_reset!(aws.tape)
        resize!(aws.tipreds, length(subject_objective.tipreds))
        copyto!(aws.tipreds, subject_objective.tipreds)

        loglik = _extended_kalman_filter_continuous!(ws, values,
            subject_objective.data, subject_objective.timesteps,
            subject_objective.params, subject_objective.tdpreds,
            subject_objective.tipreds, subject_objective.subject,
            subject_objective.max_timestep, tape)
        # _finite_deep, not isfinite: under ctsem_hessian this runs at Dual and
        # isfinite tests the value alone, so a NaN partial would pass.
        if !_finite_deep(loglik)
            valid[c] = false
            badvalue[c] = loglik
            return nothing
        end
        total += loglik

        shared || fill!(aws.theta_bar, zero(T))
        _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
        # Kept for the shared parameter layer below, and for the deferred
        # Frechet pass; with either of those enabled every subject sees the
        # same parameter layer, so any subject's values do.
        last_subject_values = tape.subject_values
        last_tipreds = subject_objective.tipreds
        if !shared
            _ctsem_parameter_layer!(gradient, aws.theta_bar, tape.subject_values,
                subject_objective.params, aws, subject_objective.tipreds)
        end
    end

    if shared && last_subject_values !== nothing
        # Every subject in this chunk shares one parameter layer; unwind once.
        _ctsem_parameter_layer!(gradient, aws.theta_bar, last_subject_values,
            params, aws, last_tipreds)
    end

    # The deferred matrix-exponential Frechet contribution (see
    # `_flush_frechet!`), pushed through the parameter layer in one extra pass.
    # It only ever touches JAx cells, and `defer_frechet` guarantees every
    # subject shares the same parameter layer for those, so one pass covers all
    # of this chunk's subjects. `jax_positions` lists the JAx component's flat
    # positions in the same column-major order the matrix is stored in.
    if aws.defer_frechet && last_subject_values !== nothing
        _flush_frechet!(aws)
        fill!(aws.theta_bar, zero(T))
        @inbounds for k in eachindex(aws.jax_positions)
            aws.theta_bar[aws.jax_positions[k]] = aws.jax_bar_deferred[k]
        end
        _ctsem_parameter_layer!(gradient, aws.theta_bar, last_subject_values,
            params, aws, last_tipreds)
    end

    totals[c] = total
    return nothing
end

"""Return a cached adjoint workspace for `objective`, building it on first use."""
function _get_or_init_adjoint_workspaces!(objective::CTSEMObjective, ::Type{T},
    nvalues::Integer, nchunks::Integer) where {T}
    cached = objective.adjoint_ws
    if cached isa Vector{Any} && length(cached) >= nchunks &&
        all(w -> w isa CTSEMAdjointWorkspace{T} &&
            length(w.regular_dual_scratch) == nvalues, view(cached, 1:nchunks))
        return cached
    end
    ntdpred = isempty(objective.subject_objectives) ? 0 :
        size(first(objective.subject_objectives).tdpreds, 1)
    # One workspace per chunk, owned by that chunk's task. Not indexed by
    # `threadid()`: a task can migrate between threads at any yield point, so
    # thread-indexed mutable scratch is a race rather than an optimisation.
    built = Any[CTSEMAdjointWorkspace(T, objective.params, nvalues, ntdpred)
                for _ in 1:nchunks]
    objective.adjoint_ws = built
    return built
end

"""
Report the state of the adjoint implementation.

`available == true` since the acceptance gates in
`test/test_adjoint_gradient_validation.jl` (adjoint vs. ForwardDiff vs. central
finite differences, across linear, cross-effect, state/PARS-expression, free
covariance/mean, TD/TI-predictor, bounded-substep, partially missing, fully
missing and unequal-length multi-subject models) became part of the always-run
suite and went green. `fallback == false` is a promise, not a description:
selecting the adjoint never silently reverts to ForwardDiff.

It is not yet the *default* gradient -- see the performance table in
`docs/src/adjoint-roadmap.md`; below roughly 40 free parameters ForwardDiff is
still faster.
"""
ctsem_adjoint_status() = (
    available=true,
    phase="tape-replay reverse EKF",
    reference="primitive pullbacks in adjoint_primitives.jl; reverse filter in adjoint_ekf.jl",
    fallback=false,
)

"""
    ctsem_validate_forward_gradient(objective, values; step=1e-6, adjoint=true)

Three-way gradient comparison: ForwardDiff, central finite differences, and
(when `adjoint` is true) the reverse-mode adjoint.

This is the acceptance-gate oracle the adjoint roadmap requires. Finite
differences are the independent referee: ForwardDiff and the adjoint are both
"clever" and could in principle share a misunderstanding of the primal, whereas
a central difference only knows the primal's outputs.
"""
function ctsem_validate_forward_gradient(objective, values::AbstractVector;
    step::Real=1e-6, adjoint::Bool=true)
    forward = ForwardDiff.gradient(objective, values)
    finite = similar(forward)
    for i in eachindex(values)
        plus = copy(values); minus = copy(values)
        plus[i] += step; minus[i] -= step
        finite[i] = (objective(plus) - objective(minus)) / (2step)
    end
    denom = max(norm(forward), eps(eltype(values)))
    adjoint_gradient = (adjoint && objective isa CTSEMObjective) ?
        ctsem_adjoint_gradient(objective, collect(values)).gradient : nothing
    return (
        forward=forward,
        finite=finite,
        adjoint=adjoint_gradient,
        relative_error=norm(forward - finite) / denom,
        adjoint_relative_error=adjoint_gradient === nothing ? nothing :
            norm(adjoint_gradient - forward) / denom,
    )
end

export ctsem_hessian

"""
    ctsem_hessian(objective, values; chunk=0)

The Hessian of the log posterior, by forward-mode differentiation *of the
reverse-mode gradient*.

The alternative the R side used before this existed is a central finite
difference of the same gradient, which costs `2 * npar` reverse sweeps and is
accurate to roughly the square root of machine precision -- and only if the
step happens to suit the parameter's scale, which one global step cannot do for
a vector mixing log standard deviations with unconstrained correlations. Nesting
forward over reverse costs `ceil(npar / chunksize)` sweeps instead of `2 * npar`,
and is exact to machine precision: no step to choose, and nothing to tune.

Forward-over-reverse rather than forward-over-forward because the reverse pass
is where this engine's work already is. `ForwardDiff.hessian` would need
`O(npar^2 / chunksize)` primal passes; differentiating the adjoint needs
`O(npar / chunksize)`, each one a single traced forward sweep plus its reverse.

The nesting works because every layer below is generic in its element type: the
workspace is `CTSEMAdjointWorkspace{T}`, built for whatever `T` the values
arrive as, and the transform layer's own dual scratch is
`Dual{Nothing,T,1}` -- so with a dual `T` it simply becomes a nested dual. The
workspace cache is keyed on `T`, so a Hessian call rebuilds it once and the
next ordinary gradient rebuilds it back; that is one model-inspection pass, not
a per-evaluation cost.

The result is symmetrised. The exact Hessian is symmetric, and each entry is
computed once, so the two triangles differ only by floating-point association
order; averaging them is free and keeps the matrix usable by a Cholesky.
"""
function ctsem_hessian(objective::CTSEMObjective, values::AbstractVector;
    chunk::Integer=0)
    x = collect(Float64, values)
    n = length(x)
    n == 0 && return zeros(Float64, 0, 0)
    gradient_of = y -> ctsem_adjoint_gradient(objective, y).gradient
    chunksize = chunk > 0 ? min(Int(chunk), n) : ForwardDiff.pickchunksize(n)
    config = ForwardDiff.JacobianConfig(gradient_of, x, ForwardDiff.Chunk{chunksize}())
    hessian = ForwardDiff.jacobian(gradient_of, x, config)
    return (hessian .+ transpose(hessian)) ./ 2
end
