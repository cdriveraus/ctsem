"""
Laplace-approximate marginal likelihood for subject-level parameter random
effects.

# What this replaces

ctsem's existing route for individual differences (`intoverpop`) *augments the
latent state*: each varying parameter becomes an extra, static latent state
carrying its own population variance, and the ordinary Kalman filter integrates
it out along with the dynamic states. That is exact for a linear-Gaussian model
and needs no new machinery, but every subject then pays for a state space of
dimension `nlatent + nindvarying` -- and the filter's cost is cubic in that
dimension, in the Lyapunov solve, the matrix exponential and its Frechet
derivative alike. Twelve random effects on a four-latent model is a sixteen-
dimensional system per row, for every subject, at every evaluation.

This file takes the other route. Each subject keeps the *small* `nlatent`
system, and the random effects are integrated out per subject by a Laplace
approximation over a `k = nindvarying` dimensional integral. The per-subject
inner problems are completely independent, so the curvature that has to be
factorized is `nsubjects` separate `k x k` blocks rather than one large system.

# The model

Exactly the one the generated Stan model already states for `intoverpop == 0`
(`ctModelWriter.R`, `rawindparams[indvaryingindex] += rawpopcovchol *
baseindparams[si]`):

    raw_i = theta + scatter(L * z_i) + TI-predictor effects,   z_i ~ N(0, I)

with `L` the Cholesky factor of the raw-scale population covariance, built from
`theta` by `_laplace_popchol` in the same parameterisation Stan uses. The
difference is only in how `z_i` is handled: Stan samples it, this integrates it
out. TI-predictor effects are untouched and stay where they are, inside the
engine's own per-subject parameter layer -- both contributions are additive on
the raw scale, so they compose without either knowing about the other.

# The approximation

Write the inner objective for subject `i`, with the standard-normal density's
own exponent folded in:

    g_i(z) = ll_i(theta, z) - z'z / 2

Then, because the `(2*pi)^(-k/2)` in `N(z|0,I)` cancels the `(2*pi)^(k/2)` the
Laplace approximation produces,

    log integral_i  =  g_i(zhat_i) - logdet(-H_i) / 2

with `zhat_i` the inner mode and `H_i` the inner curvature there. No constants
survive; that cancellation is asserted in `test_laplace.jl` against a
deliberately Gaussian integrand, where Laplace is exact.

# The outer gradient

`H_i` depends on `theta` both directly and through `zhat_i(theta)`, so an exact
outer gradient needs third derivatives of the process likelihood. That is the
requirement that stopped this feature's earlier Stan-based attempt: Stan Math's
higher-order primitives would not instantiate for the generated ctsem model.
Here it is simply `ForwardDiff` over the engine's existing forward-over-reverse
Hessian, which nests because every layer below is generic in its scalar type.

`zhat_i(theta)`'s own derivative comes from one Newton step taken in dual
arithmetic from the converged primal mode (`_laplace_dual_mode`). Since the
inner gradient's primal part is zero at the mode, that step contributes nothing
to the primal and exactly `-H^-1 dg/dtheta` to the dual -- which is the implicit
function theorem, obtained without ever forming `dH/dtheta` by hand.

`gradient_method = :approximate` drops the log-determinant's dependence on
`theta` and keeps only the envelope term, which costs one reverse sweep per
subject instead of one per parameter chunk per subject.

It is a warm-start and exploration tool, not a cheaper route to the same
answer. The dropped term is `-tr(H^-1 dH/dtheta)/2`, and `H` depends on the
population scales directly, so those are precisely the parameters whose
gradient is most wrong without it. On the linear test model in
`tests/testthat/test-julia-laplace.R` the approximate route settles some 35 log
likelihood units short of the exact one. The *value* reported is the true
Laplace value in both cases -- only the gradient differs -- so the two are
directly comparable, and which was used is carried through to the fit rather
than left to be inferred.
"""

using LinearAlgebra
using ForwardDiff

export CTSEMLaplaceSpec, CTSEMLaplaceObjective, ctsem_laplace_objective,
    ctsem_laplace_evaluate, ctsem_laplace_modes, ctsem_laplace_optimize,
    ctsem_laplace_popcov

"""
    CTSEMLaplaceSpec(re_index, sd_index, cor_index, sd_scale)

Which raw parameters vary between subjects, and which raw parameters describe
how much they vary.

All four fields index the *same* flat raw parameter vector the rest of the
engine already works in; the population covariance parameters simply live in a
tail of that vector which no model matrix cell reads. Keeping them there rather
than in a separate vector means the outer optimizer, the prior term, the
gradient and the Hessian all keep operating on one contiguous parameter vector,
with no packing and unpacking anywhere.

  * `re_index[j]`: raw position of the `j`-th varying parameter. This is Stan's
    `indvaryingindex`.
  * `sd_index[j]`: raw position of that parameter's population scale parameter,
    before transformation (Stan's `rawpopsdbase`).
  * `cor_index`: raw positions of the unconstrained correlation parameters
    (Stan's `sqrtpcov`), in *column-major lower-triangular* order -- the order
    Stan's own `counter` walks, so the two parameter vectors correspond
    element for element.
  * `sd_scale[j]`: the model's `sdscale` multiplier for that parameter.
"""
struct CTSEMLaplaceSpec
    re_index::Vector{Int}
    sd_index::Vector{Int}
    cor_index::Vector{Int}
    sd_scale::Vector{Float64}

    function CTSEMLaplaceSpec(re_index, sd_index, cor_index, sd_scale)
        re = Vector{Int}(collect(re_index))
        sd = Vector{Int}(collect(sd_index))
        cor = Vector{Int}(collect(cor_index))
        scale = Vector{Float64}(collect(sd_scale))
        k = length(re)
        length(sd) == k ||
            throw(DimensionMismatch("one population scale parameter per random effect is required"))
        length(scale) == k ||
            throw(DimensionMismatch("one sdscale per random effect is required"))
        expected = div(k * (k - 1), 2)
        length(cor) == expected || throw(DimensionMismatch(
            "expected $(expected) correlation parameters for $(k) random effects, got $(length(cor))"))
        allunique(re) || throw(ArgumentError("random-effect parameter indices must be distinct"))
        return new(re, sd, cor, scale)
    end
end

"""Number of random effects per subject."""
nrandomeffects(spec::CTSEMLaplaceSpec) = length(spec.re_index)

"""
    CTSEMLaplaceObjective(objective, spec; inner_maxiter, inner_tol)

A `CTSEMObjective` plus the random-effect structure to integrate out of it.

The inner modes are *state*, not output: they are retained between calls and
warm-start the next evaluation's Newton solve. Across an outer optimizer's
trajectory consecutive parameter vectors are close, so the inner solve usually
converges in one or two steps after the first evaluation.
"""
mutable struct CTSEMLaplaceObjective{O}
    objective::O
    spec::CTSEMLaplaceSpec
    # k x nsubjects. Column `i` is subject `i`'s current inner mode.
    modes::Matrix{Float64}
    inner_maxiter::Int
    inner_tol::Float64
    # Adjoint workspaces, one per scalar type in flight. The engine's own cache
    # lives on the CTSEMObjective and holds one type at a time, which would
    # thrash badly here: an exact outer gradient alternates between Float64 and
    # two levels of nested Dual within a single evaluation.
    workspaces::Dict{Any,Any}
    # Filled by the last evaluation; see `ctsem_laplace_diagnostics`.
    inner_iterations::Vector{Int}
    inner_gradient::Vector{Float64}
    inner_converged::Vector{Bool}
    hessian_repaired::Vector{Bool}
end

function CTSEMLaplaceObjective(objective::CTSEMObjective, spec::CTSEMLaplaceSpec;
    inner_maxiter::Integer=50, inner_tol::Real=1e-10)
    nsubjects = length(objective.subject_objectives)
    k = nrandomeffects(spec)
    return CTSEMLaplaceObjective{typeof(objective)}(objective, spec,
        zeros(Float64, k, nsubjects), Int(inner_maxiter), Float64(inner_tol),
        Dict{Any,Any}(), zeros(Int, nsubjects), zeros(Float64, nsubjects),
        falses(nsubjects), falses(nsubjects))
end

"""
    ctsem_laplace_objective(objective; re_index, sd_index, cor_index, sd_scale, ...)

Build a `CTSEMLaplaceObjective` from flat index vectors, which is the form the
R side sends across the bridge.

The index vectors are *keyword* arguments with empty defaults because of that
bridge: JuliaConnectoR deadlocks marshalling a zero-length vector, so R must be
able to omit one rather than send it. `cor_index` is empty for exactly one real
model shape -- a single random effect, which has no correlations -- and that
shape is common enough that it cannot be an error.
"""
function ctsem_laplace_objective(objective::CTSEMObjective; re_index=Int[],
    sd_index=Int[], cor_index=Int[], sd_scale=Float64[],
    inner_maxiter::Integer=50, inner_tol::Real=1e-10)
    spec = CTSEMLaplaceSpec(Int.(re_index), Int.(sd_index), Int.(cor_index),
        Float64.(sd_scale))
    return CTSEMLaplaceObjective(objective, spec; inner_maxiter=inner_maxiter,
        inner_tol=inner_tol)
end

"""Positional convenience form, for calls written in Julia."""
ctsem_laplace_objective(objective::CTSEMObjective, re_index, sd_index, cor_index,
    sd_scale; kwargs...) =
    ctsem_laplace_objective(objective; re_index=re_index, sd_index=sd_index,
        cor_index=cor_index, sd_scale=sd_scale, kwargs...)

################################################################################
# The population covariance
################################################################################

"""
    _laplace_popchol(values, spec)

The Cholesky factor of the raw-scale population covariance.

Term for term the generated Stan model's construction (`ctModelWriter.R`
around `rawpopsd = ...` through `rawpopcovchol = cholesky_decompose(...)`),
including both of its numerical offsets and its exact correlation
parameterisation, so that a Laplace fit and a Stan fit of the same model are
describing the same population distribution rather than two similar ones.

`constraincorsqrt1` reads only the off-diagonal entries of its argument, so the
scales sitting on the diagonal of `base` are inert there; they are written in
anyway to keep the correspondence with Stan's `rawpopcovbase` literal.
"""
function _laplace_popchol(values::AbstractVector{T}, spec::CTSEMLaplaceSpec) where {T}
    k = nrandomeffects(spec)
    scales = Vector{T}(undef, k)
    @inbounds for j in 1:k
        raw = values[spec.sd_index[j]]
        scales[j] = log1p_exp(2 * raw - 1) * spec.sd_scale[j] + 1e-10
    end
    base = zeros(T, k, k)
    counter = 0
    @inbounds for j in 1:k
        base[j, j] = scales[j]
        for i in 1:k
            if i > j
                counter += 1
                base[i, j] = 2 / (1 + exp(-values[spec.cor_index[counter]])) - 1
            end
        end
    end
    corsqrt = constraincorsqrt1(base)
    correlation = corsqrt * transpose(corsqrt)
    scaled = scales .+ 1e-8
    covariance = (scaled .* correlation) .* transpose(scaled)
    symmetric = (covariance .+ transpose(covariance)) ./ 2
    return Matrix(cholesky(Symmetric(symmetric)).L)
end

"""
    ctsem_laplace_popcov(laplace, values)

The raw-scale population covariance matrix implied by `values`. Reporting
helper; the fit itself only ever needs its Cholesky factor.
"""
function ctsem_laplace_popcov(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    L = _laplace_popchol(collect(Float64, values), laplace.spec)
    return L * transpose(L)
end

"""
    _laplace_subject_values(values, spec, L, z)

The raw parameter vector subject `i` is filtered with: the population vector,
shifted on the varying positions by `L * z`.

The engine's own per-subject parameter layer then adds TI-predictor effects to
this on top, exactly as it does without random effects. Both shifts are
additive on the raw scale, so neither has to know about the other.
"""
function _laplace_subject_values(values::AbstractVector{T}, spec::CTSEMLaplaceSpec,
    L::AbstractMatrix, z::AbstractVector) where {T}
    S = promote_type(T, eltype(L), eltype(z))
    shifted = Vector{S}(undef, length(values))
    copyto!(shifted, values)
    offset = L * z
    @inbounds for j in eachindex(spec.re_index)
        shifted[spec.re_index[j]] += offset[j]
    end
    return shifted
end

################################################################################
# One subject, one parameter vector
################################################################################

"""
    _laplace_workspace!(laplace, ::Type{T}, nvalues)

An adjoint workspace for scalar type `T`, cached on the Laplace object.

The engine's own cache holds a single type and rebuilds when it changes, which
is right for its callers (a Hessian call, then ordinary gradients again) and
wrong for this one: an exact outer gradient uses `Float64` for the primal inner
solve and nested duals for the outer sweep, interleaved per subject. Caching
per type turns that thrash into one build per type per model.
"""
function _laplace_workspace!(laplace::CTSEMLaplaceObjective, ::Type{T},
    nvalues::Integer) where {T}
    key = (T, Int(nvalues))
    cached = get(laplace.workspaces, key, nothing)
    cached === nothing || return cached
    objective = laplace.objective
    ntdpred = isempty(objective.subject_objectives) ? 0 :
        size(first(objective.subject_objectives).tdpreds, 1)
    built = CTSEMAdjointWorkspace(T, objective.params, Int(nvalues), ntdpred)
    laplace.workspaces[key] = built
    return built
end

"""
    _laplace_subject_value_gradient!(gradient, subject_objective, aws, values)

Subject `i`'s log likelihood at `values`, and its gradient with respect to
`values`, written into `gradient`.

This is one iteration of `_ctsem_subject_gradient_chunk!`'s loop, lifted out
because the Laplace layer needs each subject evaluated at a *different*
parameter vector -- which is exactly what the summed adjoint is built never to
have to do. The two shortcuts that make the summed gradient fast, the shared
parameter layer and the cross-subject Frechet batch, are per-subject-incorrect
for the same reason they are in `ctsem_subject_gradients`, so the Frechet batch
is flushed within the subject here too.
"""
function _laplace_subject_value_gradient!(gradient::AbstractVector{T},
    subject_objective, aws, values::AbstractVector{T}) where {T}
    ws = _get_or_init_objective_workspace!(subject_objective, T)
    tape = _tape_reset!(aws.tape)
    resize!(aws.tipreds, length(subject_objective.tipreds))
    copyto!(aws.tipreds, subject_objective.tipreds)
    aws.frechet_pending = false
    deferred = aws.defer_frechet
    aws.defer_frechet = false
    try
        loglik = _extended_kalman_filter_continuous!(ws, values,
            subject_objective.data, subject_objective.timesteps,
            subject_objective.params, subject_objective.tdpreds,
            subject_objective.tipreds, subject_objective.subject,
            subject_objective.max_timestep, tape)
        isfinite(loglik) || return loglik
        fill!(aws.theta_bar, zero(T))
        _ctsem_reverse_tape!(tape, subject_objective.params, aws, aws.n, aws.m)
        fill!(gradient, zero(T))
        _ctsem_parameter_layer!(gradient, aws.theta_bar, tape.subject_values,
            subject_objective.params, aws, subject_objective.tipreds)
        return loglik
    finally
        aws.defer_frechet = deferred
    end
end

"""
    _laplace_inner_objective_gradient(laplace, i, values, L, z, aws)

`(g_i(z), dg_i/dz)` for subject `i`: the process log likelihood at the shifted
parameter vector less `z'z/2`, and its gradient with respect to `z`.

The chain rule through the shift is one matrix-vector product: the likelihood's
gradient with respect to the shifted raw vector, restricted to the varying
positions, pushed back through `L`.
"""
function _laplace_inner_objective_gradient(laplace::CTSEMLaplaceObjective, i::Integer,
    values::AbstractVector{T}, L::AbstractMatrix{T}, z::AbstractVector{T},
    aws) where {T}
    spec = laplace.spec
    shifted = _laplace_subject_values(values, spec, L, z)
    gradient = Vector{T}(undef, length(shifted))
    loglik = _laplace_subject_value_gradient!(gradient, laplace.objective.subject_objectives[i],
        aws, shifted)
    isfinite(loglik) || return (value=loglik, gradient=fill(T(NaN), length(z)))
    restricted = Vector{T}(undef, length(spec.re_index))
    @inbounds for j in eachindex(spec.re_index)
        restricted[j] = gradient[spec.re_index[j]]
    end
    inner = transpose(L) * restricted .- z
    value = loglik - dot(z, z) / 2
    return (value=value, gradient=inner)
end

"""
    _laplace_inner_hessian(laplace, i, values, L, z, ...)

`d2 g_i / dz dz` -- the `k x k` inner curvature, by forward-mode differentiation
of the inner gradient above.

Only `k` forward directions are seeded, not one per model parameter: the
integral is over `z` alone, so this is the small block the whole design exists
to keep small. Its cost is `k` chunked dual sweeps of one subject, whatever the
model's parameter count.
"""
function _laplace_inner_hessian(laplace::CTSEMLaplaceObjective, i::Integer,
    values::AbstractVector{T}, L::AbstractMatrix{T}, z::AbstractVector{T}) where {T}
    k = length(z)
    k == 0 && return zeros(T, 0, 0)
    inner_of = function (zz)
        S = eltype(zz)
        aws = _laplace_workspace!(laplace, S, length(values))
        vs = convert(Vector{S}, values)
        Ls = convert(Matrix{S}, L)
        return _laplace_inner_objective_gradient(laplace, i, vs, Ls, zz, aws).gradient
    end
    H = ForwardDiff.jacobian(inner_of, collect(z))
    return (H .+ transpose(H)) ./ 2
end

################################################################################
# The inner mode
################################################################################

"""
    _laplace_solve_mode!(laplace, i, values, L)

Newton's method on `g_i`, warm-started from the retained mode for subject `i`.

`g_i` is a log likelihood minus a quadratic, so its curvature is negative
definite near the mode and Newton is the right method; away from the mode, and
for a nonlinear process model, it need not be. Two guards, both visible in the
diagnostics rather than silent:

  * a curvature that is not negative definite is shifted until it is, and the
    subject is flagged as repaired;
  * a step that does not improve `g_i` is halved, up to a fixed number of
    times, before the iteration gives up.

The mode is written back into `laplace.modes` so the next outer evaluation
starts from it.
"""
function _laplace_solve_mode!(laplace::CTSEMLaplaceObjective, i::Integer,
    values::AbstractVector{Float64}, L::AbstractMatrix{Float64})
    spec = laplace.spec
    k = nrandomeffects(spec)
    z = Vector{Float64}(laplace.modes[:, i])
    aws = _laplace_workspace!(laplace, Float64, length(values))
    repaired = false
    converged = false
    iterations = 0
    current = _laplace_inner_objective_gradient(laplace, i, values, L, z, aws)
    if !isfinite(current.value)
        # A warm start can be stranded outside the support after a large outer
        # step. The origin is always inside it: z = 0 is the population mean.
        fill!(z, 0.0)
        current = _laplace_inner_objective_gradient(laplace, i, values, L, z, aws)
    end
    H = zeros(Float64, k, k)
    for iteration in 1:laplace.inner_maxiter
        iterations = iteration
        if maximum(abs, current.gradient) < laplace.inner_tol
            converged = true
            break
        end
        H = _laplace_inner_hessian(laplace, i, values, L, z)
        negative_definite, H = _laplace_negate_definite(H)
        repaired |= !negative_definite
        step = H \ current.gradient        # H here is already -curvature
        accepted = false
        scale = 1.0
        for _ in 1:20
            candidate = z .+ scale .* step
            trial = _laplace_inner_objective_gradient(laplace, i, values, L, candidate, aws)
            if isfinite(trial.value) && trial.value >= current.value - 1e-12
                z = candidate
                current = trial
                accepted = true
                break
            end
            scale /= 2
        end
        accepted || break
    end
    if maximum(abs, current.gradient) < laplace.inner_tol
        converged = true
    end
    laplace.modes[:, i] = z
    laplace.inner_iterations[i] = iterations
    laplace.inner_gradient[i] = k == 0 ? 0.0 : maximum(abs, current.gradient)
    laplace.inner_converged[i] = converged
    laplace.hessian_repaired[i] = repaired
    return (z=z, value=current.value, converged=converged)
end

"""
    _laplace_negate_definite(H)

`(-H, was_negative_definite)` with `-H` made positive definite if it was not.

`H` is the inner curvature, which the approximation requires to be negative
definite: `-H` is the precision matrix of the Gaussian being fitted to the
integrand, and its determinant is the approximation's normalizing constant. A
non-negative-definite `H` means the point is not a maximum, which happens on
the way to one and can happen at one for a badly identified model. Shifting the
diagonal is the standard repair; reporting that it happened is what keeps it
from being a silent change of objective.
"""
function _laplace_negate_definite(H::AbstractMatrix{T}) where {T}
    negated = -(H .+ transpose(H)) ./ 2
    size(negated, 1) == 0 && return (true, negated)
    factorization = cholesky(Symmetric(negated); check=false)
    issuccess(factorization) && return (true, negated)
    shift = sqrt(eps(real(float(one(T)))))
    scale = maximum(abs, diag(negated))
    scale = isfinite(scale) && scale > 0 ? scale : one(real(float(one(T))))
    for _ in 1:30
        candidate = negated + (shift * scale) * I
        if issuccess(cholesky(Symmetric(candidate); check=false))
            return (false, candidate)
        end
        shift *= 10
    end
    return (false, negated + (scale + one(scale)) * I)
end

"""
    _laplace_dual_mode(laplace, i, values, zhat, Hneg, L, aws)

The inner mode as a function of the outer parameters, to first order.

One Newton step from the converged primal mode, taken with dual parameters.
The inner gradient's *primal* part is zero there, so the step's primal part is
zero and its dual part is exactly `-H^-1 dg/dtheta`: the implicit function
theorem, without forming that cross-derivative explicitly. Using the primal
`Hneg` for the solve rather than a dual one is not an approximation for the
same reason -- any dual part of the inverse would multiply a zero primal
gradient.
"""
function _laplace_dual_mode(laplace::CTSEMLaplaceObjective, i::Integer,
    values::AbstractVector{T}, zhat::AbstractVector{Float64},
    Hneg::AbstractMatrix{Float64}, L::AbstractMatrix{T}, aws) where {T}
    isempty(zhat) && return Vector{T}(undef, 0)
    z0 = convert(Vector{T}, zhat)
    inner = _laplace_inner_objective_gradient(laplace, i, values, L, z0, aws)
    return z0 .+ (Hneg \ inner.gradient)
end

################################################################################
# The Laplace objective
################################################################################

"""
    _laplace_subject_term(laplace, i, values, L, z)

`g_i(z) - logdet(-d2 g_i/dz dz) / 2`: subject `i`'s contribution to the
approximated log marginal likelihood.

The `(2*pi)^(k/2)` the Laplace approximation produces cancels the
`(2*pi)^(-k/2)` in the standard-normal density of `z`, so no dimension-dependent
constant appears here. `test_laplace.jl` pins that against a Gaussian integrand
whose integral is known in closed form.
"""
function _laplace_subject_term(laplace::CTSEMLaplaceObjective, i::Integer,
    values::AbstractVector{T}, L::AbstractMatrix{T}, z::AbstractVector{T},
    aws) where {T}
    inner = _laplace_inner_objective_gradient(laplace, i, values, L, z, aws)
    isfinite(inner.value) || return inner.value
    isempty(z) && return inner.value
    H = _laplace_inner_hessian(laplace, i, values, L, z)
    negated = -(H .+ transpose(H)) ./ 2
    factorization = cholesky(Symmetric(negated); check=false)
    issuccess(factorization) || return T(NaN)
    return inner.value - logdet(factorization) / 2
end

"""
    _laplace_popchol_jacobian(values, spec)

`(positions, jacobian)`: where the population covariance parameters sit in the
raw vector, and the derivative of `vec(L)` with respect to them.

Cheap regardless of the model, because `L` is a `k x k` Cholesky of something
built only from those parameters -- no filter, no data, no process model. This
is what lets the envelope gradient stay at one reverse sweep per subject: the
random effects' dependence on the population parameters is `(dL/dp) * z`, and
`dL/dp` is the same for every subject.
"""
function _laplace_popchol_jacobian(values::AbstractVector{Float64},
    spec::CTSEMLaplaceSpec)
    positions = vcat(spec.sd_index, spec.cor_index)
    isempty(positions) && return (positions, zeros(Float64, 0, 0))
    chol_of = function (p)
        v = convert(Vector{eltype(p)}, values)
        @inbounds for (slot, position) in enumerate(positions)
            v[position] = p[slot]
        end
        return vec(_laplace_popchol(v, spec))
    end
    return (positions, ForwardDiff.jacobian(chol_of, values[positions]))
end

"""
    _laplace_envelope_gradient(laplace, values, L, hessians)

The approximate outer gradient: the log determinant's parameter dependence
dropped, the mode held fixed, everything else exact.

    dL/dtheta ~= sum_i [ dll_i/dv  +  (dL/dp * zhat_i)' restricted(dll_i/dv) ]

The first term is the engine's own reverse pass at subject `i`'s shifted
parameter vector; the second is how that subject's shift moves when the
population scales and correlations move, which is `dL/dp` -- computed once for
all subjects -- contracted with the same gradient. One reverse sweep per
subject, and no forward directions over the model parameters at all.
"""
function _laplace_envelope_gradient(laplace::CTSEMLaplaceObjective,
    values::AbstractVector{Float64}, L::AbstractMatrix{Float64},
    hessians::Vector{Matrix{Float64}})
    spec = laplace.spec
    k = nrandomeffects(spec)
    npar = length(values)
    nsubjects = length(laplace.objective.subject_objectives)
    aws = _laplace_workspace!(laplace, Float64, npar)
    positions, chol_jacobian = _laplace_popchol_jacobian(values, spec)

    total = zeros(Float64, npar)
    subject_gradient = Vector{Float64}(undef, npar)
    for i in 1:nsubjects
        z = Vector{Float64}(laplace.modes[:, i])
        shifted = _laplace_subject_values(values, spec, L, z)
        loglik = _laplace_subject_value_gradient!(subject_gradient,
            laplace.objective.subject_objectives[i], aws, shifted)
        isfinite(loglik) || return fill(NaN, npar)
        total .+= subject_gradient
        k == 0 && continue
        # (dL/dp_j * z) . restricted gradient, for each population parameter.
        @inbounds for slot in eachindex(positions)
            derivative = reshape(view(chol_jacobian, :, slot), k, k)
            accumulated = 0.0
            for a in 1:k
                shift = 0.0
                for b in 1:k
                    shift += derivative[a, b] * z[b]
                end
                accumulated += subject_gradient[spec.re_index[a]] * shift
            end
            total[positions[slot]] += accumulated
        end
    end
    _ctsem_log_prior_gradient!(total, laplace.objective, values)
    return total
end

"""
    ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:exact)

The approximated log marginal likelihood, and optionally its gradient.

The inner modes are always found in ordinary `Float64` arithmetic first: they
are the solution of an optimization problem, and an optimizer's iterates carry
no useful derivative information. Only the converged mode does, and it gets it
from `_laplace_dual_mode`.

`gradient_method`:

  * `:exact` differentiates the complete per-subject term, log determinant and
    implicit mode dependence included. This is the third-order path.
  * `:approximate` keeps the log determinant out of the differentiation, so the
    gradient omits `-tr(H^-1 dH/dtheta)/2` and the mode is held fixed (which
    costs nothing extra: the envelope theorem makes that term vanish anyway).
    The *value* returned is the same true Laplace value in both cases; only the
    gradient differs, and `approximate = true` is returned alongside it so no
    caller has to infer which one it got. It converges somewhere materially
    different -- see the note at the top of this file -- so it is for warm
    starts and exploration, not for final estimates.
"""
function ctsem_laplace_evaluate(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    gradient::Bool=true, gradient_method=:exact, contributions::Bool=false)
    method = Symbol(gradient_method)
    method in (:exact, :approximate) ||
        throw(ArgumentError("gradient_method must be :exact or :approximate, got :$(method)"))
    theta = collect(Float64, values)
    nsubjects = length(laplace.objective.subject_objectives)
    k = nrandomeffects(laplace.spec)

    # 1. Inner modes and the value at them, in primal arithmetic, warm-started
    #    from the last call. Each subject's term is its own approximated log
    #    marginal likelihood, which is the per-subject quantity that means the
    #    same thing here as `subject_loglik` does without random effects.
    #
    #    The curvature is computed once and used three times -- for the
    #    definiteness check, for the log determinant, and as the primal solve in
    #    `_laplace_dual_mode` below. It is the most expensive primal quantity
    #    here, so recomputing it for each of those would be a third of the
    #    primal pass thrown away.
    L = _laplace_popchol(theta, laplace.spec)
    primal_hessians = Vector{Matrix{Float64}}(undef, nsubjects)
    subject_loglik = zeros(Float64, nsubjects)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    value = 0.0
    for i in 1:nsubjects
        _laplace_solve_mode!(laplace, i, theta, L)
        z = Vector{Float64}(laplace.modes[:, i])
        H = k == 0 ? zeros(Float64, 0, 0) : _laplace_inner_hessian(laplace, i, theta, L, z)
        _, negated = _laplace_negate_definite(H)
        primal_hessians[i] = negated
        inner = _laplace_inner_objective_gradient(laplace, i, theta, L, z, aws)
        term = if !isfinite(inner.value)
            inner.value
        elseif k == 0
            inner.value
        else
            factorization = cholesky(Symmetric(negated); check=false)
            issuccess(factorization) ? inner.value - logdet(factorization) / 2 : NaN
        end
        isfinite(term) || return (value=term, gradient=gradient ? fill(NaN, length(theta)) : nothing,
            subject_loglik=subject_loglik, approximate=method === :approximate,
            converged=all(laplace.inner_converged))
        subject_loglik[i] = term
        value += term
    end
    value += _ctsem_log_prior(laplace.objective, theta)

    gradient || return (value=value, gradient=nothing, subject_loglik=subject_loglik,
        approximate=method === :approximate, converged=all(laplace.inner_converged))

    # 3. The gradient.
    #
    # The approximate route does not need the forward sweep at all. With the
    # mode held fixed -- which the envelope theorem says costs nothing, since
    # `dg/dz` is zero there -- what is left is `dll_i/dtheta` at the shifted
    # parameter vector, which is exactly what the engine's reverse pass already
    # returns. That is *one reverse sweep per subject* rather than one per
    # parameter chunk per subject, and on a 4-latent model with 44 parameters
    # the difference is the whole reason the cheap mode exists.
    if method === :approximate
        return (value=value,
            gradient=_laplace_envelope_gradient(laplace, theta, L, primal_hessians),
            subject_loglik=subject_loglik, approximate=true,
            converged=all(laplace.inner_converged))
    end

    total_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Ld = _laplace_popchol(x, laplace.spec)
        accumulated = zero(S)
        for i in 1:nsubjects
            zhat = Vector{Float64}(laplace.modes[:, i])
            zd = _laplace_dual_mode(laplace, i, x, zhat, primal_hessians[i], Ld, wsd)
            accumulated += _laplace_subject_term(laplace, i, x, Ld, zd, wsd)
        end
        return accumulated + _ctsem_log_prior(laplace.objective, x)
    end
    grad = ForwardDiff.gradient(total_of, theta)
    return (value=value, gradient=grad, subject_loglik=subject_loglik,
        approximate=method === :approximate, converged=all(laplace.inner_converged))
end

"""
    ctsem_laplace_subject_values(laplace, values)

Each subject's own raw parameter vector at the current inner modes: the
population vector shifted by that subject's random effects, and then by its
TI-predictor effects.

Both shifts come from the code that already applies them during fitting --
`_laplace_subject_values` and the engine's own `_materialize_subject_values!`
-- rather than being reconstructed by the caller. That matters because the
caller is the R side's subject-parameter reporting, and a second, slightly
different copy of "what parameters does this subject have" is exactly the kind
of divergence that shows up as a summary disagreeing with the fit.

Returns `nsubjects x length(values)`; push a row through the model's transforms
to get that subject's parameter matrices.
"""
function ctsem_laplace_subject_values(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    theta = collect(Float64, values)
    L = _laplace_popchol(theta, laplace.spec)
    nsubjects = length(laplace.objective.subject_objectives)
    out = zeros(Float64, nsubjects, length(theta))
    buffer = Float64[]
    for i in 1:nsubjects
        _laplace_solve_mode!(laplace, i, theta, L)
        z = Vector{Float64}(laplace.modes[:, i])
        shifted = _laplace_subject_values(theta, laplace.spec, L, z)
        subject = laplace.objective.subject_objectives[i]
        _materialize_subject_values!(buffer, shifted, subject.params, subject.tipreds)
        out[i, :] = buffer
    end
    return out
end

export ctsem_laplace_subject_values

"""
    ctsem_laplace_population(laplace, values)

Raw-scale population standard deviations and correlations, for a whole matrix
of raw parameter vectors at once (`values[s, :]` is one draw).

Batched because the caller is the summary, which has a posterior sample rather
than a point: one call for a thousand draws instead of a thousand calls. The
correlations come back in the same column-major lower-triangular order as the
parameters that produced them, which is the order `sdcovsqrt2cov`'s and Stan's
own correlation coordinates use.
"""
function ctsem_laplace_population(laplace::CTSEMLaplaceObjective, values::AbstractMatrix)
    nsamples = size(values, 1)
    k = nrandomeffects(laplace.spec)
    noffdiagonals = div(k * (k - 1), 2)
    sd = zeros(Float64, nsamples, k)
    correlation = zeros(Float64, nsamples, noffdiagonals)
    for s in 1:nsamples
        L = _laplace_popchol(collect(Float64, view(values, s, :)), laplace.spec)
        covariance = L * transpose(L)
        scales = sqrt.(max.(diag(covariance), 0.0))
        sd[s, :] = scales
        counter = 0
        for j in 1:k, i in (j + 1):k
            counter += 1
            denominator = scales[i] * scales[j]
            correlation[s, counter] = denominator > 0 ? covariance[i, j] / denominator : 0.0
        end
    end
    return (sd=sd, correlation=correlation)
end

export ctsem_laplace_population

"""
    ctsem_laplace_modes(laplace, values)

The per-subject random-effect modes at `values`, on both the standardized
`z` scale and the raw parameter scale, with their conditional standard errors.

The conditional covariance of `z_i` is the inverse of the inner precision at
the mode, which is the same matrix the log determinant is taken of, so this
costs nothing that the objective did not already compute. On the raw scale it
is `L * cov * L'`: the same change of variables the model itself applies.
"""
function ctsem_laplace_modes(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    theta = collect(Float64, values)
    nsubjects = length(laplace.objective.subject_objectives)
    k = nrandomeffects(laplace.spec)
    L = _laplace_popchol(theta, laplace.spec)
    z = zeros(Float64, nsubjects, k)
    raw = zeros(Float64, nsubjects, k)
    zsd = zeros(Float64, nsubjects, k)
    rawsd = zeros(Float64, nsubjects, k)
    for i in 1:nsubjects
        _laplace_solve_mode!(laplace, i, theta, L)
        zi = Vector{Float64}(laplace.modes[:, i])
        z[i, :] = zi
        raw[i, :] = L * zi
        if k > 0
            H = _laplace_inner_hessian(laplace, i, theta, L, zi)
            _, negated = _laplace_negate_definite(H)
            covariance = inv(Symmetric(negated))
            zsd[i, :] = sqrt.(max.(diag(covariance), 0.0))
            rawcov = L * covariance * transpose(L)
            rawsd[i, :] = sqrt.(max.(diag(rawcov), 0.0))
        end
    end
    return (z=z, raw=raw, z_sd=zsd, raw_sd=rawsd,
        parameter=laplace.spec.re_index,
        converged=copy(laplace.inner_converged),
        iterations=copy(laplace.inner_iterations))
end

"""
    ctsem_laplace_diagnostics(laplace)

Inner-solve status from the last evaluation, per subject.
"""
ctsem_laplace_diagnostics(laplace::CTSEMLaplaceObjective) = (
    iterations=copy(laplace.inner_iterations),
    max_gradient=copy(laplace.inner_gradient),
    converged=copy(laplace.inner_converged),
    hessian_repaired=copy(laplace.hessian_repaired),
)

export ctsem_laplace_diagnostics

"""
    ctsem_laplace_optimize(laplace, start; ...)

Maximize the approximated log marginal likelihood with L-BFGS, mirroring
`ctsem_optimize`'s contract so the R side can treat the two the same way.
"""
function ctsem_laplace_optimize(laplace::CTSEMLaplaceObjective, start::AbstractVector;
    maxiter::Integer=1000, g_tol::Real=1e-8, f_tol::Real=0.0, x_tol::Real=0.0,
    verbose::Bool=false, gradient_method=:exact)
    start_values = collect(Float64, start)
    invalid_objective = floatmax(Float64) / 1e8
    gradient_limit = sqrt(floatmax(Float64))
    fg! = function (F, G, x)
        result = try
            ctsem_laplace_evaluate(laplace, x; gradient=G !== nothing,
                gradient_method=gradient_method)
        catch
            nothing
        end
        valid = result !== nothing && isfinite(result.value)
        if valid && G !== nothing
            valid = all(isfinite, result.gradient) &&
                all(abs(value) < gradient_limit for value in result.gradient)
        end
        if !valid
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        G !== nothing && (G .= -result.gradient)
        return F === nothing ? nothing : -result.value
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol, f_reltol=f_tol,
        x_abstol=x_tol, show_trace=verbose, store_trace=false)
    result = Optim.optimize(Optim.only_fg!(fg!), start_values, Optim.LBFGS(), options)
    minimizer = collect(Optim.minimizer(result))
    final = ctsem_laplace_evaluate(laplace, minimizer; gradient=true,
        gradient_method=gradient_method)
    return (
        minimizer=minimizer,
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        approximate=final.approximate,
        iterations=Optim.iterations(result),
        converged=Optim.converged(result),
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
        inner_converged=all(laplace.inner_converged),
        inner_iterations=copy(laplace.inner_iterations),
        hessian_repaired=copy(laplace.hessian_repaired),
    )
end

"""
    ctsem_laplace_hessian(laplace, values; step)

The outer Hessian of the approximated log marginal likelihood, for population
parameter standard errors.

A central difference of the *exact* outer gradient, not of the value: the
gradient is already third-order accurate, so differencing it once costs `2n`
gradient evaluations and inherits their accuracy, where differencing the value
twice would cost `O(n^2)` and lose half the digits. Nesting `ForwardDiff` once
more would be a fourth derivative of the process model; that is left until
there is evidence the difference matters, and this route is the one whose
error is at least bounded and reportable.
"""
function ctsem_laplace_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    step::Real=1e-4, gradient_method=:exact)
    x = collect(Float64, values)
    n = length(x)
    H = zeros(Float64, n, n)
    for j in 1:n
        h = step * max(1.0, abs(x[j]))
        plus = copy(x); plus[j] += h
        minus = copy(x); minus[j] -= h
        gp = ctsem_laplace_evaluate(laplace, plus; gradient=true, gradient_method=gradient_method).gradient
        gm = ctsem_laplace_evaluate(laplace, minus; gradient=true, gradient_method=gradient_method).gradient
        H[:, j] = (gp .- gm) ./ (2h)
    end
    return (H .+ transpose(H)) ./ 2
end

export ctsem_laplace_hessian

################################################################################
# Behaving like the objective it wraps
################################################################################
#
# The R side reaches the engine through a handful of entry points that all take
# "the objective", and a Laplace fit has to answer them too. They split cleanly
# in three.
#
#  1. Questions about the *model* -- its parameter layout, its matrices, which
#     cells are state-dependent. Integrating the random effects out does not
#     change any of them, so these forward to the wrapped objective unchanged.
#
#  2. Questions about the *fit* -- value, gradient, curvature, per-subject
#     scores. These have Laplace counterparts that mean the same thing about
#     the same model, so the generic entry points route to them. That is what
#     makes `ctOptimUncertainty` and the summary work without a Laplace branch
#     on the R side: they ask for the log posterior's curvature, and they get
#     the log *marginal* posterior's curvature, which is the right answer to
#     their question.
#
#  3. Questions that are genuinely not answered yet -- filtering and data
#     generation, both of which need a per-subject parameter vector rather than
#     one shared vector. These throw. Forwarding them to the wrapped objective
#     would run and return the population-level answer, silently, where the
#     caller asked for a subject-level one; a refusal is the honest result
#     until the subject-conditional versions exist.

for f in (:ctsem_parameter_layout, :ctsem_state_dependent_cells, :ctsem_parameter_matrices)
    @eval $f(laplace::CTSEMLaplaceObjective, args...; kwargs...) =
        $f(laplace.objective, args...; kwargs...)
end

"""
Evaluate a Laplace objective through the generic entry point.

Returns the approximated log *marginal* likelihood and its gradient. The
`gradient_method` names the process gradient elsewhere in the engine and has no
meaning here -- the Laplace gradient is a forward sweep over the reverse pass
either way -- so it is accepted and ignored rather than rejected, and the
choice that does matter is `laplace_gradient`.
"""
function ctsem_evaluate(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint,
    laplace_gradient=:exact)
    result = ctsem_laplace_evaluate(laplace, values; gradient=gradient,
        gradient_method=laplace_gradient)
    contributions || return (value=result.value, gradient=result.gradient)
    # No `row_loglik`: the integral is over a whole subject's trajectory, so a
    # single row has no marginal contribution to report. Returning the subject
    # terms and omitting the row ones is more honest than inventing a
    # decomposition that the approximation does not have.
    return (value=result.value, gradient=result.gradient,
        subject_loglik=result.subject_loglik)
end

"""The outer Hessian, for the generic curvature entry point."""
ctsem_hessian(laplace::CTSEMLaplaceObjective, values::AbstractVector; chunk::Integer=0) =
    ctsem_laplace_hessian(laplace, values)

"""The gradient of the approximated log marginal likelihood."""
function ctsem_adjoint_gradient(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    return (value=result.value, gradient=result.gradient)
end

"""
    ctsem_subject_gradients(laplace, values)

Per-subject scores: row `i` is the gradient of subject `i`'s own approximated
log marginal likelihood.

This is the quantity the OPG, sandwich and score-bootstrap uncertainty methods
consume, and for a marginal likelihood it is the marginal score rather than the
joint one -- the two differ by exactly the random-effect terms that have been
integrated out, so using the joint score would understate the uncertainty it is
there to measure.

It costs one sweep, not one per subject: the per-subject terms are assembled
into a vector and differentiated together, so the same forward directions serve
every row.
"""
function ctsem_subject_gradients(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    theta = collect(Float64, values)
    nsubjects = length(laplace.objective.subject_objectives)
    k = nrandomeffects(laplace.spec)

    L = _laplace_popchol(theta, laplace.spec)
    primal_hessians = Vector{Matrix{Float64}}(undef, nsubjects)
    for i in 1:nsubjects
        _laplace_solve_mode!(laplace, i, theta, L)
        z = Vector{Float64}(laplace.modes[:, i])
        H = k == 0 ? zeros(Float64, 0, 0) : _laplace_inner_hessian(laplace, i, theta, L, z)
        _, negated = _laplace_negate_definite(H)
        primal_hessians[i] = negated
    end

    terms_of = function (x)
        S = eltype(x)
        wsd = _laplace_workspace!(laplace, S, length(x))
        Ld = _laplace_popchol(x, laplace.spec)
        out = Vector{S}(undef, nsubjects)
        for i in 1:nsubjects
            zhat = Vector{Float64}(laplace.modes[:, i])
            zd = _laplace_dual_mode(laplace, i, x, zhat, primal_hessians[i], Ld, wsd)
            out[i] = _laplace_subject_term(laplace, i, x, Ld, zd, wsd)
        end
        return out
    end
    scores = ForwardDiff.jacobian(terms_of, theta)

    # Each subject carries its share of the prior, so the rows still sum to the
    # full posterior gradient -- the same convention `ctsem_subject_gradients`
    # uses without random effects.
    if !isempty(laplace.objective.prior_index) && nsubjects > 0
        share = 1 / nsubjects
        for i in 1:nsubjects
            _ctsem_log_prior_gradient!(view(scores, i, :), laplace.objective, theta, share)
        end
    end
    value = sum(_laplace_subject_term(laplace, i, theta, L,
        Vector{Float64}(laplace.modes[:, i]),
        _laplace_workspace!(laplace, Float64, length(theta))) for i in 1:nsubjects)
    return (value=value + _ctsem_log_prior(laplace.objective, theta), scores=scores)
end

for (f, what) in ((:ctsem_kalman, "Filtering"), (:ctsem_generate, "Data generation"))
    @eval function $f(laplace::CTSEMLaplaceObjective, args...; kwargs...)
        throw(ArgumentError(string($what, " is not implemented for the Laplace ",
            "random-effect route yet. It needs a per-subject parameter vector ",
            "rather than one shared vector, and returning the population-level ",
            "answer instead would be a different quantity than the one asked ",
            "for. Use intoverpop=TRUE for this.")))
    end
end
