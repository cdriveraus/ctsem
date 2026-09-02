using Optim

# Execute a Julia function repeatedly behind a single JuliaConnectoR call. The
# last result is returned so the compiler cannot discard the calls. Keeping the
# loop on the Julia side amortizes the fixed R-to-Julia bridge latency.
function cts_benchmark_repeat(f, iterations::Integer, args...)
    iterations > 0 || throw(ArgumentError("iterations must be positive"))

    result = f(args...)
    for _ in 2:iterations
        result = f(args...)
    end
    return result
end

# Build reusable per-subject workspaces. This removes setup cost from cached
# gradient benchmarks and mirrors the lower-level Julia API.
function cts_benchmark_subject_workspace(values, params, timesteps_by_subject, data_by_subject)
    workspaces = [
        ContinuousTimeSEM.ContinuousEKFGradientWorkspace(
            values,
            params,
            timesteps_by_subject[idx],
            data_by_subject[idx],
        )
        for idx in eachindex(timesteps_by_subject)
    ]
    return (workspaces = workspaces, gradient = similar(values))
end

# Objective/workspace variant that refreshes state-dependent transforms several
# times within each observation interval. This is important for strongly
# nonlinear dynamics such as Lotka--Volterra oscillations.
mutable struct CTSBenchmarkSubstepObjective{P,D,TS,S}
    value_ws::Any
    dual_ws::Any
    params::P
    data::D
    timesteps::TS
    max_step::S
end

function _cts_benchmark_substep_likelihood!(ws, values, objective)
    params = objective.params
    data = objective.data
    timesteps = objective.timesteps
    _materialize_subject_values!(ws.subject_values, values, params, Float64[])
    _materialize_all_params!(ws.all_params, ws.subject_values, params)
    pars = ws.pars
    all_params = getdata(pars)
    all_params .= ws.all_params

    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.T0VAR, 0, ws.state_dim)
    copyto!(ws.P_predict.data, ws.bufferQ.out)
    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0, ws.manifest_dim)
    copyto!(ws.state, pars.T0MEANS)
    first_context = CTSEMRowContext(ws.state, pars, Float64[], Float64[],
        timesteps[1], zero(eltype(values)), 1, 1)
    apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
        params.update_transforms, first_context)
    _ekf_update_step!(ws, pars, data, 1)
    log2π = log(2π)
    likelihood = _kalman_loglikelihood_cholesky!(
        ws.ll_buffer, ws.S, ws.ỹ, log2π, ws.manifest_dim,
    )

    previous_time = timesteps[1]
    @inbounds for observation in 2:length(timesteps)
        interval = timesteps[observation] - previous_time
        n_substeps = max(1, ceil(Int, interval / objective.max_step))
        substep = interval / n_substeps
        for _ in 1:n_substeps
            context = CTSEMRowContext(ws.state, pars, Float64[], Float64[],
                timesteps[observation], substep, 1, observation)
            apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
                params.predict_transforms, context)
            _ekf_predict_step!(ws, pars, substep)
            copyto!(ws.P_update.data, ws.P_predict.data)
        end

        context = CTSEMRowContext(ws.state, pars, Float64[], Float64[],
            timesteps[observation], interval, 1, observation)
        apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
            params.update_transforms, context)
        ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0, ws.manifest_dim)
        _ekf_update_step!(ws, pars, data, observation)
        likelihood += _kalman_loglikelihood_cholesky!(
            ws.ll_buffer, ws.S, ws.ỹ, log2π, ws.manifest_dim,
        )
        previous_time = timesteps[observation]
    end
    return likelihood
end

function (objective::CTSBenchmarkSubstepObjective)(values::AbstractVector{T}) where {T}
    workspace_ref = T <: ForwardDiff.Dual ? :dual_ws : :value_ws
    ws = getfield(objective, workspace_ref)
    if ws === nothing || eltype(ws.all_params) != T
        ws = _init_continuous_ekf_workspace(T, objective.params)
        setfield!(objective, workspace_ref, ws)
    end
    return _cts_benchmark_substep_likelihood!(ws, values, objective)::T
end

struct CTSBenchmarkSubstepGradientWorkspace{O,C,R,DR}
    objective::O
    cfg::C
    result::R
    diff_result::DR
end

function CTSBenchmarkSubstepGradientWorkspace(values, params, timesteps, data, max_step)
    objective = CTSBenchmarkSubstepObjective(
        nothing, nothing, params, data, timesteps, max_step,
    )
    return CTSBenchmarkSubstepGradientWorkspace(
        objective,
        ForwardDiff.GradientConfig(objective, values),
        similar(values),
        DiffResults.GradientResult(values),
    )
end

function grad_log_likelihood_ekf_continuous!(result, values, workspace::CTSBenchmarkSubstepGradientWorkspace)
    ForwardDiff.gradient!(result, workspace.objective, values, workspace.cfg)
    return result
end

function res_and_grad_likelihood_ekf_continuous!(workspace::CTSBenchmarkSubstepGradientWorkspace, values)
    ForwardDiff.gradient!(
        workspace.diff_result, workspace.objective, values, workspace.cfg,
    )
    workspace.result .= DiffResults.gradient(workspace.diff_result)
    return (
        value = DiffResults.value(workspace.diff_result),
        gradient = workspace.result,
    )
end

function cts_benchmark_substep_workspace(values, benchmark_workspace, max_step = 0.02)
    length(benchmark_workspace.workspaces) == 1 || throw(ArgumentError(
        "benchmark_workspace must contain exactly one subject",
    ))
    original = only(benchmark_workspace.workspaces).objective
    workspace = CTSBenchmarkSubstepGradientWorkspace(
        values, original.params, original.timesteps, original.data, max_step,
    )
    return (workspaces = [workspace], gradient = similar(values))
end

# Sum subject log-likelihoods using the reusable workspaces.
function cts_benchmark_loglikelihood_subjects_cached(benchmark_workspace, values)
    total = zero(eltype(values))
    @inbounds for workspace in benchmark_workspace.workspaces
        total += workspace.objective(values)
    end
    return total
end

# Sum subject gradients using the reusable workspaces.
function cts_benchmark_grad_subjects_cached(benchmark_workspace, values)
    result = benchmark_workspace.gradient
    fill!(result, zero(eltype(result)))
    @inbounds for workspace in benchmark_workspace.workspaces
        result .+= ContinuousTimeSEM.grad_log_likelihood_ekf_continuous!(
            workspace.result,
            values,
            workspace,
        )
    end
    return result
end

# Sum subject likelihoods and gradients using the reusable workspaces.
function cts_benchmark_res_and_grad_subjects_cached(benchmark_workspace, values)
    result = benchmark_workspace.gradient
    fill!(result, zero(eltype(result)))
    total = zero(eltype(values))
    @inbounds for workspace in benchmark_workspace.workspaces
        subject_result = ContinuousTimeSEM.res_and_grad_likelihood_ekf_continuous!(
            workspace,
            values,
        )
        total += subject_result.value
        result .+= subject_result.gradient
    end
    return (value = total, gradient = result)
end

# Optim minimizes, while the SEM target is a log-likelihood to maximize. This
# fg! callback therefore returns the negative summed likelihood and writes the
# negative summed gradient into the Optim gradient buffer.
function cts_benchmark_negative_loglikelihood_fg_cached!(F, G, values, benchmark_workspace)
    if G === nothing
        if F === nothing
            return nothing
        end
        return -cts_benchmark_res_and_grad_subjects_cached(benchmark_workspace, values).value
    end

    if F === nothing
        gradient = cts_benchmark_grad_subjects_cached(benchmark_workspace, values)
        @. G = -gradient
        return nothing
    end

    result = cts_benchmark_res_and_grad_subjects_cached(benchmark_workspace, values)
    @. G = -result.gradient
    return -result.value
end

# Construct the reusable Optim objective for the cached multi-subject EKF
# likelihood. The initial vector defines the parameter dimension and element
# type used by Optim.
function cts_benchmark_cached_objective(initial_values, benchmark_workspace)
    fg! = let benchmark_workspace = benchmark_workspace
        (F, G, values) -> cts_benchmark_negative_loglikelihood_fg_cached!(
            F,
            G,
            values,
            benchmark_workspace,
        )
    end
    return Optim.OnceDifferentiable(Optim.only_fg!(fg!), initial_values)
end

# Optimize the summed Julia EKF log-likelihood with LBFGS and return plain
# values that JuliaConnectoR can expose cleanly to R.
function cts_benchmark_optimize_cached(initial_values, benchmark_workspace, max_iterations::Integer = 5000)
    starting_values = copy(initial_values)
    objective = cts_benchmark_cached_objective(starting_values, benchmark_workspace)
    result = Optim.optimize(
        objective,
        starting_values,
        Optim.LBFGS(),
        Optim.Options(iterations = max_iterations),
    )

    return (
        minimizer = Optim.minimizer(result),
        minimum = Optim.minimum(result),
        maximum = -Optim.minimum(result),
        iterations = Optim.iterations(result),
        converged = Optim.converged(result),
    )
end

# Time the full Julia optimization procedure. This is the Julia-side analogue of
# the system.time(ctStanFit(..., optimize = TRUE)) call below.
function cts_benchmark_time_optimize_cached(initial_values, benchmark_workspace, max_iterations::Integer = 5000)
    GC.gc()
    start_time = time_ns()
    result = cts_benchmark_optimize_cached(initial_values, benchmark_workspace, max_iterations)
    elapsed = (time_ns() - start_time) / 1e9

    return (
        elapsed = elapsed,
        minimizer = result.minimizer,
        minimum = result.minimum,
        maximum = result.maximum,
        iterations = result.iterations,
        converged = result.converged,
    )
end

# Integrate the fitted deterministic state equation without measurement
# updates. A fourth-order Runge--Kutta scheme makes this useful for smooth
# diagnostic curves of nonlinear models instead of filtered interpolation.
function cts_benchmark_predict_trajectory(values, benchmark_workspace, timesteps, max_step = 0.02)
    length(benchmark_workspace.workspaces) == 1 || throw(ArgumentError(
        "benchmark_workspace must contain exactly one subject",
    ))
    objective = only(benchmark_workspace.workspaces).objective
    params = objective.params
    ws = _init_continuous_ekf_workspace(eltype(values), params)
    _materialize_subject_values!(ws.subject_values, values, params, Float64[])
    _materialize_all_params!(ws.all_params, ws.subject_values, params)

    pars = ws.pars
    all_params = getdata(pars)
    all_params .= ws.all_params
    state = collect(vec(pars.T0MEANS))
    state_dim = length(state)
    manifest_dim = size(pars.LAMBDA, 1)
    states = Matrix{eltype(values)}(undef, state_dim, length(timesteps))
    fitted_manifest = Matrix{eltype(values)}(undef, manifest_dim, length(timesteps))
    k1, k2, k3, k4, trial = ntuple(_ -> similar(state), 5)

    function derivative!(output, at_state)
        context = CTSEMRowContext(at_state, pars, Float64[], Float64[],
            zero(eltype(values)), zero(eltype(values)), 1, 1)
        apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
            params.predict_transforms, context)
        mul!(output, pars.DRIFT, at_state)
        output .+= vec(pars.CINT)
        return output
    end

    function record!(index)
        states[:, index] .= state
        context = CTSEMRowContext(state, pars, Float64[], Float64[],
            zero(eltype(values)), zero(eltype(values)), 1, index)
        apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
            params.update_transforms, context)
        mul!(@view(fitted_manifest[:, index]), pars.LAMBDA, state)
        @views fitted_manifest[:, index] .+= vec(pars.MANIFESTMEANS)
        return nothing
    end

    record!(1)
    previous_time = timesteps[1]
    @inbounds for index in 2:length(timesteps)
        interval = timesteps[index] - previous_time
        interval >= 0 || throw(ArgumentError("timesteps must be nondecreasing"))
        n_steps = max(1, ceil(Int, interval / max_step))
        step = interval / n_steps
        for _ in 1:n_steps
            derivative!(k1, state)
            @. trial = state + (step / 2) * k1
            derivative!(k2, trial)
            @. trial = state + (step / 2) * k2
            derivative!(k3, trial)
            @. trial = state + step * k3
            derivative!(k4, trial)
            @. state += (step / 6) * (k1 + 2k2 + 2k3 + k4)
        end
        record!(index)
        previous_time = timesteps[index]
    end

    return (states = states, fitted_manifest = fitted_manifest)
end
