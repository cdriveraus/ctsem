"""
particle_filter.jl

A ground-truth likelihood for the nonlinear model: a bootstrap particle filter.

# Why it exists

The filter this engine fits with is an assumed-density approximation. It
carries a Gaussian for the latent state, linearises the drift over each
prediction step, and projects every non-Gaussian measurement update back onto
a Gaussian. For a linear model with Gaussian indicators it is exact; for
anything else its likelihood is an approximation whose error nobody can see
from inside it. This file computes the same marginal likelihood without any of
those approximations, so the two can be compared at the fitted values, row by
row.

# What it does

Per subject, `N` particles are drawn from the initial state distribution and
carried through the observation intervals. Between rows each particle steps
through the process with the drift, diffusion and transforms evaluated at *its
own* state, so state dependence is exact; the step is either the filter's
exponential discretisation applied at the particle (`:exponential`, exact for
the linear parts of the model, second order otherwise) or plain Euler-Maruyama
(`:euler`, nothing linearised over the step at all). At each row every particle
is weighted by the conditional density of the observations given its state --
the same `_ctsem_state_row!` the state-explicit path uses, so Gaussian,
binary, ordinal, count and censored indicators, time-dependent predictors and
state-dependent measurement parameters all come for free -- and the particles
are resampled when the weights degenerate. The log of the mean weight, row by
row, is an unbiased estimate of the marginal likelihood (Del Moral 2004), with
Monte Carlo error that falls as `1 / sqrt(N)`.

Random effects carried as augmented states ride along as static particle
coordinates and are integrated exactly by the same mechanism. The Laplace
route's random effects are not handled here.

# What it is for

A diagnostic at fixed parameters, not an estimator: the estimate is not smooth
in the parameters, and it costs `N` transform evaluations per particle step.
`ctParticleLik()` on the R side runs it at a fit's estimates and reports the
difference to the filter's likelihood, in total and per row, with the Monte
Carlo standard error beside it. The substep count is the reference mesh, chosen
by the caller and independent of the fit's; doubling it and running again is
the convergence check.
"""

using Random

export ctsem_particle_loglik

"""
    ctsem_particle_loglik(objective, values; particles=2000, substeps=20,
                          transition=:exponential, seed=1, resample_threshold=0.5)

The bootstrap particle filter estimate of the marginal log likelihood of
`objective`'s data at the parameter values `values`.

Returns a NamedTuple: `loglik`; `se`, an approximate Monte Carlo standard error
from the effective sample sizes (`sqrt(sum_t (N / ESS_t - 1) / N)`, a rough
guide -- independent seeds are the honest one); `row_loglik`, one increment per
row of the data in data order; `subject_loglik`; `ess_min`, the smallest
effective sample size any row saw; and the settings. `substeps` is the number
of transition steps per observation interval, floored by the objective's own
`maxtimestep` rule or mesh. `seed` fixes every draw.
"""
function ctsem_particle_loglik(objective::CTSEMObjective, values::AbstractVector;
    particles::Integer=2000, substeps::Integer=20, transition=:exponential,
    seed::Integer=1, resample_threshold::Real=0.5)
    N = Int(particles)
    N >= 2 || throw(ArgumentError("particles must be at least 2"))
    substeps >= 1 || throw(ArgumentError("substeps must be at least 1"))
    0 < resample_threshold <= 1 || throw(ArgumentError("resample_threshold must be in (0, 1]"))
    trans = _ctsem_transition(transition)
    rng = MersenneTwister(Int(seed))
    x = Vector{Float64}(values)
    subjects = objective.subject_objectives
    row_loglik = Float64[]
    subject_loglik = zeros(Float64, length(subjects))
    varsum = 0.0
    ess_min = Inf
    for (i, sub) in enumerate(subjects)
        ll, rows, v, e = _ctsem_particle_subject!(sub, objective.params, x, N,
            Int(substeps), trans, Float64(resample_threshold), rng)
        subject_loglik[i] = ll
        append!(row_loglik, rows)
        varsum += v
        ess_min = min(ess_min, e)
    end
    return (loglik=sum(subject_loglik), se=sqrt(varsum), row_loglik=row_loglik,
        subject_loglik=subject_loglik, ess_min=ess_min, particles=N,
        substeps=Int(substeps), transition=trans)
end

"""
Fold this row's log weights `logw` into the running normalised log weights
`lw`. Returns the row's log-likelihood increment and the effective sample size
after the update; a row no particle can explain returns `-Inf`.
"""
function _ctsem_pf_weigh!(lw::Vector{Float64}, logw::Vector{Float64}, N::Int)
    mx = -Inf
    @inbounds for p in 1:N
        v = lw[p] + (isnan(logw[p]) ? -Inf : logw[p])
        lw[p] = v
        mx = max(mx, v)
    end
    isfinite(mx) || return (-Inf, 0.0)
    s = 0.0
    @inbounds for p in 1:N
        s += exp(lw[p] - mx)
    end
    increment = mx + log(s)
    s2 = 0.0
    @inbounds for p in 1:N
        lw[p] -= increment
        s2 += exp(2 * lw[p])
    end
    return (increment, 1 / s2)
end

"""Systematic resampling of the particle columns of `X` by the weights `exp.(lw)`."""
function _ctsem_pf_resample!(X::Matrix{Float64}, Xnew::Matrix{Float64},
    lw::Vector{Float64}, rng, N::Int)
    n = size(X, 1)
    u = rand(rng) / N
    c = exp(lw[1])
    i = 1
    @inbounds for p in 1:N
        target = u + (p - 1) / N
        while c < target && i < N
            i += 1
            c += exp(lw[i])
        end
        for r in 1:n
            Xnew[r, p] = X[r, i]
        end
    end
    copyto!(X, Xnew)
    fill!(lw, -log(N))
    return nothing
end

function _ctsem_particle_subject!(sub, sp, x::Vector{Float64}, N::Int, nsubsteps::Int,
    transition::Symbol, threshold::Float64, rng)
    T = Float64
    ws = _get_or_init_objective_workspace!(sub, T)
    tipreds = _ctsem_tipred_vector(sub.tipreds, x)
    data = sub.data
    ts = sub.timesteps
    tdpreds = sub.tdpreds
    subject = sub.subject
    _materialize_subject_values!(ws.subject_values, x, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = ws.pars
    all_params = getdata(pars)
    all_params .= ws.all_params

    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)
    indices = ws.diffusion_state_indices
    k = length(indices)
    log2pi = log(2 * T(pi))

    factor = Matrix{T}(undef, n, n)
    qfactor = Matrix{T}(undef, k, k)
    theta = Matrix{T}(undef, m, m)
    chol = Matrix{T}(undef, m, m)
    pred = Vector{T}(undef, m)
    resid = Vector{T}(undef, m)
    scratch = Vector{T}(undef, m)
    gaussian = Vector{Int}(undef, m)
    X = Matrix{T}(undef, n, N)
    Xnew = Matrix{T}(undef, n, N)
    logw = Vector{T}(undef, N)
    lw = fill(-log(N), N)
    z = Vector{T}(undef, max(n, k, 1))
    nrows = min(length(ts), size(data, 2))
    rows = fill(-Inf, nrows)
    varsum = 0.0
    ess_min = Inf

    # Initial particles: T0MEANS + the T0VAR factor applied to standard normals.
    _ctsem_sdcor_factor!(factor, pars.T0VAR, ws.bufferQ, ws.state_dim)
    @inbounds for p in 1:N
        randn!(rng, view(z, 1:n))
        for i in 1:n
            acc = T(pars.T0MEANS[i])
            for j in 1:n
                acc += factor[i, j] * z[j]
            end
            X[i, p] = acc
        end
    end

    # Row 1: the three transform groups and the TD impulse at each particle's
    # state, in the order the filter and the state-explicit pass use, then the
    # conditional density of the observations.
    @inbounds for p in 1:N
        for i in 1:n; ws.state[i] = X[i, p]; end
        ctx = CTSEMRowContext(ws.state, pars, view(tdpreds, :, 1), tipreds,
            ts[1], zero(T), subject, 1)
        apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
            sp.predict_transforms, ctx)
        apply_complex_transforms_at_indices!(all_params, ws.td_param_indices,
            sp.td_transforms, ctx)
        _ctsem_td_impulse!(ws, pars, ctx.tdpreds)
        apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
            sp.update_transforms, ctx)
        logw[p] = _ctsem_state_row!(ws, pars, data, 1, pred, theta, chol, resid,
            scratch, gaussian, log2pi, nothing)
        for i in 1:n; X[i, p] = ws.state[i]; end
    end
    increment, ess = _ctsem_pf_weigh!(lw, logw, N)
    rows[1] = increment
    isfinite(increment) || return (-Inf, rows, varsum, 0.0)
    varsum += (N / ess - 1) / N
    ess_min = min(ess_min, ess)
    ess < threshold * N && _ctsem_pf_resample!(X, Xnew, lw, rng, N)

    prev = ts[1]
    @inbounds for t in 2:nrows
        dt = ts[t] - prev
        nsub = max(nsubsteps, _ctsem_substeps(dt, sub.max_timestep, t))
        h = dt / nsub
        for p in 1:N
            for i in 1:n; ws.state[i] = X[i, p]; end
            for substep in 1:nsub
                ctx = CTSEMRowContext(ws.state, pars, view(tdpreds, :, t), tipreds,
                    prev + substep * h, h, subject, t)
                apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
                    sp.predict_transforms, ctx)
                ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.DIFFUSION, 0, ws.state_dim)
                randn!(rng, view(z, 1:k))
                if transition === :euler && ws.continuous_time
                    _matvec_mul!(ws.bufferQ.r, pars.DRIFT, ws.state, ws.state_dim, ws.state_dim)
                    for i in 1:n
                        ws.state[i] += (ws.bufferQ.r[i] + pars.CINT[i]) * h
                    end
                    for j in 1:k, i in 1:k
                        qfactor[i, j] = ws.bufferQ.out[indices[i], indices[j]]
                    end
                    _ctsem_lower_chol!(qfactor, k)
                    sqrt_h = sqrt(h)
                    for i in 1:k
                        acc = zero(T)
                        for j in 1:i
                            acc += qfactor[i, j] * z[j]
                        end
                        ws.state[indices[i]] += sqrt_h * acc
                    end
                else
                    if ws.continuous_time
                        _compute_discrete_time_form!(ws.discrete_ca, ws.bufferQ, ws.bufferQ.out,
                            pars, h, ws.exp_buffer, ws.lyap_buffer, ws.state, indices,
                            ws.diffusion_buffer, ws.discretization_buffer, ws.state_dim,
                            ws.discretization_cache)
                    else
                        _compute_one_step_form!(ws.discrete_ca, ws.bufferQ.out, pars, ws.state,
                            indices, ws.state_dim)
                    end
                    _matvec_mul!(ws.bufferQ.r, ws.discrete_ca.dDRIFT, ws.state, ws.state_dim,
                        ws.state_dim)
                    for i in 1:n
                        ws.state[i] = ws.bufferQ.r[i] + ws.discrete_ca.dINT[i]
                    end
                    for j in 1:k, i in 1:k
                        qfactor[i, j] = ws.discrete_ca.dDIFFUSION[indices[i], indices[j]]
                    end
                    _ctsem_lower_chol!(qfactor, k)
                    for i in 1:k
                        acc = zero(T)
                        for j in 1:i
                            acc += qfactor[i, j] * z[j]
                        end
                        ws.state[indices[i]] += acc
                    end
                end
            end
            ctx = CTSEMRowContext(ws.state, pars, view(tdpreds, :, t), tipreds, ts[t], dt,
                subject, t)
            apply_complex_transforms_at_indices!(all_params, ws.td_param_indices,
                sp.td_transforms, ctx)
            _ctsem_td_impulse!(ws, pars, ctx.tdpreds)
            apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
                sp.update_transforms, ctx)
            logw[p] = _ctsem_state_row!(ws, pars, data, t, pred, theta, chol, resid,
                scratch, gaussian, log2pi, nothing)
            for i in 1:n; X[i, p] = ws.state[i]; end
        end
        increment, ess = _ctsem_pf_weigh!(lw, logw, N)
        rows[t] = increment
        isfinite(increment) || return (-Inf, rows, varsum, 0.0)
        varsum += (N / ess - 1) / N
        ess_min = min(ess_min, ess)
        ess < threshold * N && _ctsem_pf_resample!(X, Xnew, lw, rng, N)
        prev = ts[t]
    end
    return (sum(rows), rows, varsum, ess_min)
end
