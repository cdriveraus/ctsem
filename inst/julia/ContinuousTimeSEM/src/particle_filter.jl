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

export ctsem_particle_loglik, ctsem_particle_batch

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

Each subject has its own random stream, seeded from `(seed, subject)`, and the
subjects are split across threads the way the filter's own subject loop is
(`ctsem_set_max_chunks!` caps it). So the result is the same whatever the
thread count, and a difference between two runs is a difference in `seed`.
"""
function ctsem_particle_loglik(objective::CTSEMObjective, values::AbstractVector;
    particles::Integer=2000, substeps::Integer=20, transition=:exponential,
    seed::Integer=1, resample_threshold::Real=0.5)
    N = Int(particles)
    N >= 2 || throw(ArgumentError("particles must be at least 2"))
    substeps >= 1 || throw(ArgumentError("substeps must be at least 1"))
    seed >= 0 || throw(ArgumentError("seed must be non-negative"))
    0 < resample_threshold <= 1 || throw(ArgumentError("resample_threshold must be in (0, 1]"))
    trans = _ctsem_transition(transition)
    x = Vector{Float64}(values)
    subjects = objective.subject_objectives
    nsubjects = length(subjects)
    results = Vector{Tuple{Float64,Vector{Float64},Float64,Float64}}(undef, nsubjects)
    nchunks = _ctsem_nchunks(nsubjects)
    if nchunks <= 1
        for i in 1:nsubjects
            _ctsem_particle_run!(results, i, objective, x, N, Int(substeps), trans,
                Float64(resample_threshold), seed)
        end
    else
        ranges = _ctsem_chunk_ranges(nsubjects, nchunks)
        Threads.@sync for c in 1:nchunks
            Threads.@spawn for i in ranges[c]
                _ctsem_particle_run!(results, i, objective, x, N, Int(substeps), trans,
                    Float64(resample_threshold), seed)
            end
        end
    end
    row_loglik = Float64[]
    subject_loglik = zeros(Float64, nsubjects)
    varsum = 0.0
    ess_min = Inf
    @inbounds for i in 1:nsubjects
        ll, rows, v, e = results[i]
        subject_loglik[i] = ll
        append!(row_loglik, rows)
        varsum += v
        ess_min = min(ess_min, e)
    end
    return (loglik=sum(subject_loglik), se=sqrt(varsum), row_loglik=row_loglik,
        subject_loglik=subject_loglik, ess_min=ess_min, particles=N,
        substeps=Int(substeps), transition=trans)
end

# One subject into its slot. A function rather than a closure in the spawn, so
# nothing in the loop body can rebind a caller's local.
function _ctsem_particle_run!(results, i::Int, objective, x, N, nsubsteps, trans,
    threshold, seed)
    results[i] = _ctsem_particle_subject!(objective.subject_objectives[i],
        objective.params, x, N, nsubsteps, trans, threshold,
        _ctsem_particle_rng(seed, i))
    return nothing
end

# The stream for one subject: seeded from the two halves of `seed` and the
# subject index, so no two (seed, subject) pairs share a stream and the
# assignment of subjects to threads cannot change a result.
function _ctsem_particle_rng(seed::Integer, i::Integer)
    s = UInt64(seed)
    return MersenneTwister(UInt32[UInt32(s & 0xffffffff), UInt32(s >> 32), UInt32(i)])
end

"""
    ctsem_particle_batch(objective, values::AbstractMatrix; particles=1000,
                         substeps=10, transition=:exponential, seed=1,
                         resample_threshold=0.5)

The particle log likelihood at each column of `values`, beside the filter's.

Returns a NamedTuple of vectors, one entry per column: `particle`, the particle
estimate; `se` and `ess_min` as in `ctsem_particle_loglik`; `filter`, the
filter's log likelihood summed over subjects; and `posterior`, the objective's
full value at the column, which is `filter` plus the prior and any term for
missing time-independent predictors. So `posterior - filter + particle` is the
log posterior with the filter's likelihood replaced by the particle one, which
is what an importance weight against the particle posterior needs.

`seed` is one integer for every column, or a vector with one per column. The
importance weights downstream want the latter: the estimates must be
independent across draws for the self-normalised weights to be consistent. A
shared seed was tried first and is wrong for that use -- it makes the filter's
Monte Carlo error a fixed random function of the parameters, which tilts the
weighted posterior by that seed's error surface instead of averaging out, and
more draws do not cure it. Measured on a two-parameter linear fit: a shared
seed put the corrected mean 0.7 standard errors from the exact posterior with
an effective sample size of 158; independent seeds removed it.
"""
function ctsem_particle_batch(objective::CTSEMObjective, values::AbstractMatrix;
    particles::Integer=1000, substeps::Integer=10, transition=:exponential,
    seed=1, resample_threshold::Real=0.5)
    ncol = size(values, 2)
    ncol >= 1 || throw(ArgumentError("values must have at least one column"))
    seeds = seed isa AbstractVector ? seed : fill(seed, ncol)
    length(seeds) == ncol ||
        throw(ArgumentError("seed must be one integer or one per column of values"))
    particle = Vector{Float64}(undef, ncol)
    se = Vector{Float64}(undef, ncol)
    ess_min = Vector{Float64}(undef, ncol)
    filter = Vector{Float64}(undef, ncol)
    posterior = Vector{Float64}(undef, ncol)
    for j in 1:ncol
        x = Vector{Float64}(view(values, :, j))
        pf = ctsem_particle_loglik(objective, x; particles=particles, substeps=substeps,
            transition=transition, seed=Int(seeds[j]), resample_threshold=resample_threshold)
        particle[j] = pf.loglik
        se[j] = pf.se
        ess_min[j] = pf.ess_min
        full = objective(x)
        posterior[j] = full
        filter[j] = full - _ctsem_log_prior(objective, x) - _ctsem_ti_missing_loglik(objective, x)
    end
    return (particle=particle, se=se, ess_min=ess_min, filter=filter, posterior=posterior)
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
                ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.DIFFUSION, ws.covmatcode, ws.state_dim)
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
                            ws.discretization_cache, ws.affine_buffer)
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
