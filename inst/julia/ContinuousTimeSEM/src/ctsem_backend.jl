"""
    CTSEMObjective(params, subject_starts, timesteps, data)

Prepared continuous-time likelihood for the R `ctsem` backend. `data` is
stored as manifest variables by observation. Subject starts are prepared by R,
which owns data ordering and predictor preprocessing. This type retains only
the numerical per-subject EKF workspaces needed for repeated evaluation.
"""
mutable struct CTSEMObjective{P,O}
    params::P
    subject_objectives::O
    # Lazily built reverse-mode workspace (see adjoint.jl). Cached here so the
    # one-off model inspection it does -- discovering each transform's read set
    # -- happens once per objective rather than once per gradient evaluation.
    adjoint_ws::Any
end

CTSEMObjective(params, subject_objectives) =
    CTSEMObjective(params, subject_objectives, nothing)

export CTSEMObjective, ctsem_objective, ctsem_evaluate, ctsem_optimize

function _ctsem_row_loglikelihood(objective::CTSEMObjective, values::AbstractVector)
    contributions = Vector{eltype(values)}()
    for subject_objective in objective.subject_objectives
        previous = zero(eltype(values))
        nrows = size(subject_objective.data, 2)
        for row in 1:nrows
            prefix = ContinuousEKFObjective(objective.params,
                Matrix(subject_objective.data[:, 1:row]),
                collect(subject_objective.timesteps[1:row]);
                tdpreds=Matrix(subject_objective.tdpreds[:, 1:row]),
                tipreds=subject_objective.tipreds, subject=subject_objective.subject,
                max_timestep=subject_objective.max_timestep)
            current = prefix(values)
            push!(contributions, current - previous)
            previous = current
        end
    end
    contributions
end

function _ctsem_subject_ranges(subject_starts::AbstractVector, timesteps::AbstractVector,
    data::AbstractMatrix)
    n = length(timesteps)
    n > 0 || throw(ArgumentError("data must contain at least one observation"))
    n == size(data, 2) || throw(DimensionMismatch("timesteps length must equal data columns"))
    starts = Int.(subject_starts)
    !isempty(starts) || throw(ArgumentError("subject_starts must contain the first observation"))
    starts[1] == 1 || throw(ArgumentError("subject_starts must begin at one"))
    all((1 .<= starts) .& (starts .<= n)) || throw(ArgumentError("subject_starts are outside the observation range"))
    all(diff(starts) .> 0) || throw(ArgumentError("subject_starts must be strictly increasing"))
    stops = vcat(starts[2:end] .- 1, n)
    out = UnitRange{Int}[start:stop for (start, stop) in zip(starts, stops)]
    for r in out
        _validate_continuous_ekf_inputs(collect(view(timesteps, r)), Matrix(view(data, :, r)))
    end
    out
end

function CTSEMObjective(params::EKFParameters, subject_starts::AbstractVector,
    timesteps::AbstractVector, data::AbstractMatrix,
    tdpred_data::AbstractMatrix=zeros(eltype(data), 0, size(data, 2)),
    tipred_data::AbstractMatrix=zeros(eltype(data), length(subject_starts), 0),
    max_timestep::Real=Inf)
    ranges = _ctsem_subject_ranges(subject_starts, timesteps, data)
    size(tdpred_data, 2) == size(data, 2) || throw(DimensionMismatch("TD predictor columns must match observations"))
    size(tipred_data, 1) == length(ranges) || throw(DimensionMismatch("TI predictor rows must match subjects"))
    # Copy each subject once. This avoids R proxy/view lifetime issues and makes
    # the objective safe to retain in a Julia session.
    objects = Any[
        ContinuousEKFObjective(params, Matrix(view(data, :, r)), collect(view(timesteps, r));
            tdpreds=Matrix(view(tdpred_data, :, r)), tipreds=vec(tipred_data[i, :]), subject=i,
            max_timestep=max_timestep)
        for (i, r) in enumerate(ranges)
    ]
    return CTSEMObjective(params, objects)
end

ctsem_objective(params::EKFParameters, subject_starts, timesteps, data,
    tdpred_data=zeros(eltype(data), 0, size(data, 2)),
    tipred_data=zeros(eltype(data), length(subject_starts), 0), max_timestep::Real=Inf) =
    CTSEMObjective(params, subject_starts, timesteps, data, tdpred_data, tipred_data, max_timestep)

function (objective::CTSEMObjective)(values::AbstractVector)
    total = zero(eltype(values))
    # Deliberately serial for now: each subject objective owns mutable primal and
    # dual workspaces. Parallel reduction is added only with isolated worker
    # workspaces, rather than sharing these caches across threads.
    @inbounds for subject_objective in objective.subject_objectives
        total += subject_objective(values)
    end
    total
end

"""
Evaluate a prepared objective and optionally return its gradient.

`gradient_method` selects how the gradient is computed:

  * `:adjoint` -- the reverse-mode pass in `adjoint.jl`, and the default. One
    traced forward sweep plus one reverse sweep per subject, independent of
    the parameter count.
  * `:forward` -- ForwardDiff. Cost scales with the number of free parameters,
    since ForwardDiff needs one dual pass per chunk of them. Kept for
    cross-checking, and marginally quicker on very small nonlinear models.

The two are checked against each other (and against finite differences) by
`test/test_adjoint_gradient_validation.jl`. `:adjoint` never silently falls
back to `:forward`: an unsupported model or an invalid trial point produces a
non-finite gradient or a thrown error, not a quietly different answer.
"""
function ctsem_evaluate(objective::CTSEMObjective, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    # Accepts a String as well as a Symbol: R passes this across the
    # JuliaConnectoR boundary, which marshals character vectors to `String`.
    method = Symbol(gradient_method)
    method in (:forward, :adjoint) ||
        throw(ArgumentError("gradient_method must be :forward or :adjoint, got :$(method)"))
    if gradient && method === :adjoint
        result = ctsem_adjoint_gradient(objective, collect(values))
        value = result.value
        grad = result.gradient
    else
        value = objective(values)
        grad = gradient ? ForwardDiff.gradient(objective, values) : nothing
    end
    if contributions
        subject = [obj(values) for obj in objective.subject_objectives]
        return (value=value, gradient=grad, subject_loglik=subject,
            row_loglik=_ctsem_row_loglikelihood(objective, values))
    end
    return (value=value, gradient=grad)
end

"""Optimize a prepared likelihood entirely within Julia using L-BFGS."""
function ctsem_optimize(objective::CTSEMObjective, start::AbstractVector;
    maxiter::Integer=1000, g_tol::Real=1e-8, f_tol::Real=0.0,
    x_tol::Real=0.0, verbose::Bool=false, gradient_method=:adjoint)
    start_values = collect(start)
    invalid_objective = floatmax(eltype(start_values)) / 1e8
    gradient_limit = sqrt(floatmax(eltype(start_values)))
    fg! = function (F, G, x)
        result = try
            ctsem_evaluate(objective, x; gradient=G !== nothing,
                gradient_method=gradient_method)
        catch
            nothing
        end
        valid = result !== nothing && isfinite(result.value)
        if valid && G !== nothing
            valid = all(isfinite, result.gradient) && all(abs(value) < gradient_limit for value in result.gradient)
        end
        if !valid
            G !== nothing && fill!(G, zero(eltype(G)))
            return F === nothing ? nothing : invalid_objective
        end
        if G !== nothing
            G .= -result.gradient
        end
        return F === nothing ? nothing : -result.value
    end
    options = Optim.Options(iterations=Int(maxiter), g_tol=g_tol,
        f_reltol=f_tol, x_abstol=x_tol, show_trace=verbose, store_trace=false)
    result = Optim.optimize(Optim.only_fg!(fg!), start_values, Optim.LBFGS(), options)
    final = ctsem_evaluate(objective, Optim.minimizer(result); gradient=true,
        contributions=true, gradient_method=gradient_method)
    return (
        minimizer=collect(Optim.minimizer(result)),
        maximum_loglik=final.value,
        gradient=collect(final.gradient),
        subject_loglik=collect(final.subject_loglik),
        row_loglik=final.row_loglik,
        iterations=Optim.iterations(result),
        converged=Optim.converged(result),
        g_converged=Optim.g_converged(result),
        f_converged=Optim.f_converged(result),
        x_converged=Optim.x_converged(result),
    )
end
