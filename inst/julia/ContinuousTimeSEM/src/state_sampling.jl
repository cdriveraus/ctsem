"""
state_sampling.jl

The state-explicit path: latent states are *sampled* rather than integrated
out, which is what ctsem's `intoverstates=FALSE` asks for.

# What this is an alternative to

Everything else in this engine marginalises the latent states as it goes. The
EKF carries a Gaussian `(mean, covariance)` and, at each row, replaces it with
the moments of the posterior given that row's observation. For a Gaussian
indicator that update is exact. For a categorical one it is an *assumed
density* step: the true conditional is not Gaussian, so the filter integrates
the observation against the predicted state (`binary_measurement.jl`) and then
projects the answer back onto a Gaussian. The result is a legitimate
approximate marginal likelihood, and generation through it draws from that
approximation -- the density the fit maximises, which is why the round trip in
`kalman_trace.jl` holds exactly.

It is not, however, the model. A count observation far into its own tail moves
the projected state a long way, the next row's rate is larger still, and the
walk can run away; bounded indicators (binary, ordinal) hide it because their
support does not change with the state. Nothing here can do that, because
nothing here conditions the state on an observation at all.

# The parameterisation

States are built from standard normal *innovations*, not sampled directly:

    eta_1 = T0MEANS + M0 z_1,             M0 M0' = T0VAR
    eta_t = A_dt eta_{t-1} + b_dt + L z_t, L L'  = Q_dt

so the joint density of everything the model does not observe is

    log p(y, z | theta) = sum_t log p(y_t | eta_t(z, theta))
                          - z'z / 2 - dim(z) / 2 log 2pi

with `p(y_t | eta_t)` the *conditional* observation model: Gaussian for a
Gaussian indicator, Bernoulli for a binary one, cumulative logit for an
ordinal, Poisson for a count, Tobit for a censored one. Every one of those
already exists as `_category_loglikelihood`, evaluated at a known linear
predictor -- which is exactly what a sampled state provides.

This is the non-centred form, and it is the right one for three separate
reasons. It needs no inverse and no determinant of the process noise, so a
singular one -- every `intoverpop='augmented'` model has one, and so does any
model with a DIFFUSION entry fixed at zero -- is not a special case but simply
a direction the innovation cannot move. It is the parameterisation a
Hamiltonian sampler needs, for the reason `sample_density.jl` gives about the
random effects. And it is what generation wants anyway: draw the innovations,
walk the trajectory forward, draw the observations given the states.

It is also, deliberately, the object the generated Stan model builds for
`intoverstates=0`, where `etaupdbasestates` is standard normal and the states
are reconstructed from it.

# One innovation per substep, not per row

`max_timestep` splits a long interval into bounded substeps and re-materialises
the local affine model at each one. For a nonlinear model those substeps *are*
the process, so each gets its own innovation: putting a single innovation at
the end of the interval would run a deterministic trajectory through a
state-dependent drift and perturb it afterwards, which is a different process.
The count is deterministic given the data, so `ctsem_state_dimension` reports
it and the caller draws that many standard normals -- rather than both sides
recomputing the substep arithmetic and eventually disagreeing about it.

With no `maxtimestep` set, which is the ordinary case, this is one innovation
per row and the arithmetic disappears.

# What is drawn where

Only the caller has a random number generator. Every standard normal arrives
from outside, as it does for `ctsem_generate`, so `set.seed()` in R governs
generation completely.
"""

using LinearAlgebra, DiffResults

export ctsem_state_dimension, ctsem_joint_loglikelihood, ctsem_joint_evaluate,
    ctsem_generate_states

"""
Poisson rate above which a count is drawn from its normal approximation rather
than by inverting its distribution function.

Inversion walks upward one value at a time, so its cost is the value drawn, and
its first term `exp(-lambda)` underflows to zero at `lambda = 745` -- which
makes the walk not merely slow but wrong. The normal approximation has skewness
`1/sqrt(lambda)`, so at five hundred it is 0.045, far finer than the gap
between two adjacent counts out there.

Both forms are driven by the *same* standard normal, one through `Phi` and one
directly, so the value drawn moves continuously as a rate crosses the boundary
instead of jumping.
"""
const _CTSEM_POISSON_NORMAL_RATE = Ref(500.0)

"""
    _ctsem_normal_cdf(z)

`Phi(z)`, from `erfc` rather than from the rational approximation
`_standard_normal_cdf` uses.

That one is accurate to 8e-8 absolute, which is ample against a Bernoulli
threshold and is not ample for inverting a count distribution: once the rate is
large enough that individual masses are smaller than 8e-8, that much misplaced
probability is a whole category. `erfc` is already a dependency here
(`_norm_logcdf` is written with `logerfc`), so the accurate form is free.
"""
@inline _ctsem_normal_cdf(z::Real) = erfc(-z / sqrt(oftype(float(z), 2))) / 2

"""
    _ctsem_substeps(dt, rule, t)

Number of prediction substeps for the interval `dt` ending at row `t`.

`rule` is either `maxtimestep` (a positive real: the interval is split into
steps no longer than it) or a mesh, an integer vector with one entry per row
giving the count directly. The mesh is what `ctsem_auto_substeps` produces;
it travels in the same slot as the rule so every consumer of the objective --
the filter, the tape, the state-explicit dimension count, generation -- sees
one policy without a second code path.
"""
@inline _ctsem_substeps(dt, rule::Real, t::Int) = max(1, ceil(Int, dt / rule))
@inline _ctsem_substeps(dt, mesh::AbstractVector, t::Int) = Int(mesh[t])

"""NaN, for the reason `_invalid_ekf_loglikelihood` is NaN."""
@inline _ctsem_invalid(::Type{T}) where {T} = -one(T) * NaN

"""
    _ctsem_lower_chol!(L, d)

Overwrite the leading `d` by `d` block of `L` with a lower-triangular factor
satisfying `L L' = A`, and report whether every pivot was positive.

Unlike `_ctsem_cholesky!` this does not give up on a non-positive pivot: it
zeroes that column and carries on, which produces a valid factor of a
*positive semi-definite* matrix. That is what the state path needs. A process
noise covariance is legitimately singular whenever a state has no diffusion of
its own -- every augmented random-effect coordinate, and any model with a
DIFFUSION entry fixed at zero -- and the right behaviour there is an innovation
direction with no effect, not a failure.

The flag still matters, because a *density* over a singular Gaussian block does
not exist: the observation model uses it to refuse, while generation ignores it
and draws from the degenerate distribution, which is perfectly well defined.

Only the lower triangle is written and only the lower triangle is read.
"""
function _ctsem_lower_chol!(L::AbstractMatrix{T}, d::Int) where {T}
    definite = true
    @inbounds for j in 1:d
        s = L[j, j]
        for p in 1:(j - 1)
            s -= L[j, p] * L[j, p]
        end
        # Value-only comparison, as everywhere else in this engine: a Dual
        # orders on its value, so the branch a gradient takes is the branch the
        # primal took.
        if s > zero(real(T))
            u = sqrt(s)
            L[j, j] = u
            for i in (j + 1):d
                v = L[i, j]
                for p in 1:(j - 1)
                    v -= L[i, p] * L[j, p]
                end
                L[i, j] = v / u
            end
        else
            definite = false
            L[j, j] = zero(T)
            for i in (j + 1):d
                L[i, j] = zero(T)
            end
        end
    end
    return definite
end

"""
    _ctsem_mvn_logpdf(L, resid, scratch, d, log2pi)

`log N(resid; 0, L L')` by forward substitution, for lower-triangular `L`.
"""
@inline function _ctsem_mvn_logpdf(L::AbstractMatrix{T}, resid, scratch, d::Int,
    log2pi) where {T}
    quad = zero(T)
    logdet_half = zero(T)
    @inbounds for i in 1:d
        acc = resid[i]
        for j in 1:(i - 1)
            acc -= L[i, j] * scratch[j]
        end
        scratch[i] = acc / L[i, i]
        quad += scratch[i] * scratch[i]
        logdet_half += log(L[i, i])
    end
    return -quad / 2 - logdet_half - d * log2pi / 2
end

"""
    _ctsem_sdcor_factor!(dest, mat, buffer, dim)

The factor `M` with `M M'` equal to `sdcovsqrt2cov(mat)`, written to `dest`.

The same two lines `sdcovsqrt2cov!` uses to *form* the covariance, stopped one
step earlier: it builds `M` from the standard deviations on the diagonal and
the constrained correlation square root, then returns `M M'`. So this is not an
approximation of that covariance nor a re-factorisation of it -- it is the
factor the covariance was built from, which means a state drawn through it has
exactly the covariance the filter would have carried, and a zero standard
deviation gives a zero row rather than a failed Cholesky.

Reading it out of `sdcovsqrt2cov!`'s own scratch buffer would work today and
would break silently the first time that function reuses the buffer.
"""
@inline function _ctsem_sdcor_factor!(dest, mat, buffer, dim::Val{d}) where {d}
    constraincorsqrt1_vec!(buffer, mat, 1e-5, dim)
    @inbounds for j in 1:d, i in 1:d
        dest[i, j] = mat[i, i] * buffer.out[i, j]
    end
    return dest
end

"""State-only form of `_apply_td_impulse!`: there is no covariance to carry."""
@inline function _ctsem_td_impulse!(ws, pars, tdpreds)
    isempty(tdpreds) && return nothing
    _matvec_mul!(ws.bufferQ.r, pars.TDPREDEFFECT, tdpreds, ws.state_dim,
        Val(length(tdpreds)))
    ws.state .+= ws.bufferQ.r
    return nothing
end


################################################################################
# Layout
################################################################################

"""
    CTSEMStateLayout

Where each subject's innovations and data columns sit, and how many there are
in total.

The innovation count is a property of the *design* -- the observation times and
`maxtimestep` -- and not of the parameters, so it is computed once and holds
for every evaluation of a given objective.
"""
struct CTSEMStateLayout
    nlatent::Int
    ndiffusion::Int
    nmanifest::Int
    nrows::Int
    ndim::Int
    zoffsets::Vector{Int}
    rowoffsets::Vector{Int}
end

function _ctsem_subject_innovations(sub, nlatent::Int, ndiffusion::Int)
    timesteps = sub.timesteps
    nsteps = min(length(timesteps), size(sub.data, 2))
    total = nlatent
    prev = timesteps[1]
    @inbounds for t in 2:nsteps
        dt = timesteps[t] - prev
        total += ndiffusion * _ctsem_substeps(dt, sub.max_timestep, t)
        prev = timesteps[t]
    end
    return total
end

function _ctsem_state_layout(objective::CTSEMObjective)
    sp = objective.params
    pars = ComponentVector(zeros(Float64, length(sp.mutables)), sp.parameter_axis)
    nlatent = size(pars.DRIFT, 1)
    nmanifest = size(pars.LAMBDA, 1)
    ndiffusion = isempty(sp.diffusion_state_indices) ? nlatent :
        length(sp.diffusion_state_indices)
    subjects = objective.subject_objectives
    zoffsets = Vector{Int}(undef, length(subjects))
    rowoffsets = Vector{Int}(undef, length(subjects))
    at = 0
    rows = 0
    for (i, sub) in enumerate(subjects)
        zoffsets[i] = at
        rowoffsets[i] = rows
        at += _ctsem_subject_innovations(sub, nlatent, ndiffusion)
        rows += size(sub.data, 2)
    end
    return CTSEMStateLayout(nlatent, ndiffusion, nmanifest, rows, at, zoffsets,
        rowoffsets)
end

"""
    ctsem_state_dimension(objective)

How many standard normal innovations the state path needs for this design.

`nlatent` for each subject's first row, and one per diffusing state per
bounded substep thereafter. The caller draws exactly this many and hands them
back; nothing on this side has an RNG.
"""
ctsem_state_dimension(objective::CTSEMObjective) =
    _ctsem_state_layout(objective).ndim


################################################################################
# Generation record
################################################################################

"""
    CTSEMStateGenerate(base, out, states, llrow)

The generation half of the state pass.

`base` is `nmanifest` by `nrows` standard normals, used exactly as
`CTSEMGenerateSpec` uses them -- one entry per manifest cell, turned into a
uniform where the draw needs one. `out` receives the generated observations and
keeps `NaN` wherever the skeleton had no observation, `states` receives the
sampled latent state at each row, and `llrow` each row's *conditional* log
likelihood given that state.

`llrow` is therefore not comparable with the `llrow` `ctsem_generate` returns:
that one is a marginal one-step-ahead density, this one conditions on the state
that produced the row. It is `NaN` for a row whose Gaussian block is singular,
which is a real answer rather than an error -- with no measurement error a
Gaussian indicator is a deterministic function of the state, so the draw exists
and its density with respect to Lebesgue measure does not.
"""
mutable struct CTSEMStateGenerate
    base::Matrix{Float64}
    out::Matrix{Float64}
    states::Matrix{Float64}
    llrow::Vector{Float64}
    offset::Int
end

CTSEMStateGenerate(base, out, states, llrow) =
    CTSEMStateGenerate(base, out, states, llrow, 0)

"""
    _ctsem_draw_count(rate, u, z)

One Poisson draw, by inverting the distribution function below
`_CTSEM_POISSON_NORMAL_RATE` and from the normal approximation above it.

Inversion is exact and needs no rejection step, which is what keeps a draw
reproducible from the caller's standard normal alone.
"""
function _ctsem_draw_count(rate::T, u::T, z::T) where {T}
    if !(rate < T(_CTSEM_POISSON_NORMAL_RATE[]))
        return max(zero(T), round(rate + sqrt(rate) * z))
    end
    p = exp(-rate)
    cumulative = p
    k = 0
    limit = _CTSEM_COUNT_GENERATE_MAX[]
    while cumulative < u && k < limit
        k += 1
        p *= rate / k
        cumulative += p
    end
    return T(k)
end

"""
    _ctsem_draw_categorical(gen, eta, row, col, thresholds, kind)

One non-Gaussian observation drawn from its conditional distribution given the
state, through the standard normal the caller supplied for this cell.

Inverting the *conditional* distribution, not the marginal one: the state is
known here, so there is nothing to integrate and no quadrature is involved.
That is the whole difference from `_generate_binary!`, which has to integrate
the state's own uncertainty out before it can invert anything.
"""
function _ctsem_draw_categorical(gen::CTSEMStateGenerate, eta::T, row::Int,
    col::Int, thresholds, kind::Int) where {T}
    z = T(gen.base[row, col])
    u = _ctsem_normal_cdf(z)
    if kind == CTSEM_OBS_COUNT
        return _ctsem_draw_count(exp(min(eta, T(_CTSEM_COUNT_MAX_LOG_RATE[]))),
            u, z)
    elseif kind == CTSEM_OBS_CENSORED
        lower, upper, sd = _censor_limits(thresholds, T)
        return min(max(eta + sd * z, lower), upper)
    elseif kind == CTSEM_OBS_BINARY || isempty(thresholds)
        return u < inv(one(T) + exp(-eta)) ? one(T) : zero(T)
    end
    # Cumulative logit, the same identity `_category_loglikelihood` evaluates:
    # P(y <= k) is the logistic CDF at `thresholds[k] - eta`. The thresholds
    # arrive already cumulated and so are increasing, which makes this walk a
    # single pass.
    K = length(thresholds) + 1
    @inbounds for k in 1:(K - 1)
        u < inv(one(T) + exp(eta - thresholds[k])) && return T(k)
    end
    return T(K)
end


################################################################################
# The pass
################################################################################

"""
    _ctsem_state_row!(ws, pars, data, col, buffers..., gen)

One row: the conditional log likelihood of its observations given the sampled
state, drawing them first when generating.

Gaussian indicators are handled as one block, because their measurement errors
are correlated through MANIFESTVAR; every other kind is scalar given the state
and is handled on its own. Only the observed entries take part, and generation
writes only where an observation already exists -- the same contract
`_generate_row!` follows, so a generated dataset keeps the missingness pattern
it was given.
"""
function _ctsem_state_row!(ws, pars, data::AbstractMatrix, col::Int,
    pred::AbstractVector{T}, theta::AbstractMatrix{T}, chol::AbstractMatrix{T},
    resid::AbstractVector{T}, scratch::AbstractVector{T},
    gaussian::AbstractVector{Int}, log2pi, gen) where {T}

    m = _val(ws.manifest_dim)
    n = _val(ws.state_dim)
    types = ws.manifesttype
    row = gen === nothing ? col : gen.offset + col

    # The linear predictor. `pars.LAMBDA` rather than `pars.Jy`, matching the
    # filter's own distinction: LAMBDA evaluated at the current state gives the
    # manifest mean, and the Jacobian exists only to propagate a covariance --
    # of which there is none here.
    @inbounds for i in 1:m
        acc = pars.MANIFESTMEANS[i]
        for j in 1:n
            acc += pars.LAMBDA[i, j] * ws.state[j]
        end
        pred[i] = acc
    end

    ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferΘ, pars.MANIFESTVAR, 0,
        ws.manifest_dim)
    @inbounds for j in 1:m, i in 1:m
        theta[i, j] = ws.bufferΘ.out[i, j]
    end

    total = zero(T)
    ngauss = 0
    @inbounds for i in 1:m
        _ctsem_observed(data[i, col]) || continue
        kind = i <= length(types) ? types[i] : 0
        if kind == 0
            ngauss += 1
            gaussian[ngauss] = i
            continue
        end
        thresholds = _ordinal_thresholds!(ws, pars, i)
        y = gen === nothing ? data[i, col] :
            _ctsem_draw_categorical(gen, pred[i], i, row, thresholds, kind)
        gen === nothing || (gen.out[i, row] = y)
        contribution = _category_loglikelihood(pred[i], y, thresholds, kind)
        # Generating, an unusable contribution is reported and stepped over:
        # the draw itself is valid whatever its density came out as, and the
        # remaining indicators in this row still have to be drawn. Evaluating,
        # it ends the pass, exactly as a bad row ends the filter's.
        if gen === nothing && !isfinite(contribution)
            return _ctsem_invalid(T)
        end
        total += contribution
    end

    ngauss == 0 && return total

    @inbounds for j in 1:ngauss, i in 1:ngauss
        chol[i, j] = theta[gaussian[i], gaussian[j]]
    end
    definite = _ctsem_lower_chol!(chol, ngauss)
    if gen === nothing
        # No measurement error means no density: with the state known, the
        # indicator is a deterministic function of it, and a point mass has no
        # log likelihood with respect to Lebesgue measure to contribute. So a
        # Gaussian indicator with MANIFESTVAR fixed at zero -- perfectly
        # ordinary under the marginal filter, where the state's own uncertainty
        # supplies the variance -- has no joint density here at all, and this
        # is an invalid evaluation rather than a small number.
        definite || return _ctsem_invalid(T)
        @inbounds for i in 1:ngauss
            resid[i] = data[gaussian[i], col] - pred[gaussian[i]]
        end
        return total + _ctsem_mvn_logpdf(chol, resid, scratch, ngauss, log2pi)
    end

    # The draw is `mean + L z` with `z` indexed by the Gaussian rows only, so
    # each manifest cell's standard normal is consumed by exactly one draw --
    # the same discipline `_generate_row!` keeps, and what stops a Gaussian
    # indicator and a categorical one in the same row from sharing a deviate
    # and coming out correlated.
    @inbounds for i in 1:ngauss
        acc = pred[gaussian[i]]
        for j in 1:i
            acc += chol[i, j] * T(gen.base[gaussian[j], row])
        end
        gen.out[gaussian[i], row] = acc
        resid[i] = acc - pred[gaussian[i]]
    end
    definite || return T(NaN)
    return total + _ctsem_mvn_logpdf(chol, resid, scratch, ngauss, log2pi)
end

"""
    _ctsem_state_pass!(ws, params, data, timesteps, sp, tdpreds, tipreds,
                       subject, max_timestep, z, zoffset, gen)

One subject's forward pass with the states built from `z`, returning the total
conditional log likelihood of its observations.

The row contract is the filter's, step for step -- predict, time-dependent
impulse, measurement, with the bounded substeps in between and the same three
groups of state-dependent transforms materialised at the same three points. The
only difference is what the state is: here it is a draw, so the transforms see
a sampled state rather than a filtered mean, and the measurement conditions on
it rather than updating it.

Keeping that ordering identical is what makes the two paths comparable. A
nonlinear model's parameters are functions of the state, so a pass that
materialised them at a different point would be fitting a different model and
would still return a perfectly plausible number.
"""
function _ctsem_state_pass!(ws, params::AbstractVector{T}, data::AbstractMatrix,
    timesteps, sp, tdpreds::AbstractMatrix, tipreds::AbstractVector,
    subject::Int, max_timestep, z::AbstractVector, zoffset::Int, gen) where {T}

    _materialize_subject_values!(ws.subject_values, params, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = ws.pars
    all_params = getdata(pars)
    all_params .= ws.all_params

    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)
    indices = ws.diffusion_state_indices
    k = length(indices)
    log2pi = log(2 * T(pi))

    # Allocated once per subject, not once per row: a handful of small matrices
    # against a pass over every observation. The filter's own row loop
    # allocates its observed-index vector on the same terms.
    factor = Matrix{T}(undef, n, n)
    qfactor = Matrix{T}(undef, k, k)
    theta = Matrix{T}(undef, m, m)
    chol = Matrix{T}(undef, m, m)
    pred = Vector{T}(undef, m)
    resid = Vector{T}(undef, m)
    scratch = Vector{T}(undef, m)
    gaussian = Vector{Int}(undef, m)

    # Initial state: T0MEANS plus the T0VAR factor applied to the first block
    # of innovations.
    _ctsem_sdcor_factor!(factor, pars.T0VAR, ws.bufferQ, ws.state_dim)
    at = zoffset
    @inbounds for i in 1:n
        acc = T(pars.T0MEANS[i])
        for j in 1:n
            acc += factor[i, j] * T(z[at + j])
        end
        ws.state[i] = acc
    end
    at += n

    # All three groups, in the loop's own order, even though row 1 has no
    # interval to predict over: a group supplies *values* as well as a
    # prediction. PARS is in the predict group, and an update-group or td-group
    # cell -- LAMBDA, MANIFESTMEANS, MANIFESTVAR, Jy, TDPREDEFFECT, Jtd -- may
    # be written by a transform that reads one. Running td and update only left
    # that read pointing at a parameter slot nothing had written, so it took
    # the zero the buffer is filled with: 6.0 log units on one row of a
    # LAMBDA-reads-PARS model, and this is the path `intoverstates` takes for a
    # categorical model. The filter's row 1 carries the same three groups in
    # the same order, which is what the docstring above means by materialising
    # them at the same three points.
    #
    # One context serves all three, as it already served two. Its interval is
    # zero, which is the only thing about it a predict-group transform could
    # object to, and none can: generate_complex_transform_string substitutes
    # only state, PARS and the model matrices, so no transform expression can
    # reference ctx.dt or ctx.time at all.
    first_context = CTSEMRowContext(ws.state, pars, view(tdpreds, :, 1), tipreds,
        timesteps[1], zero(T), subject, 1)
    apply_complex_transforms_at_indices!(all_params, ws.predict_param_indices,
        sp.predict_transforms, first_context)
    apply_complex_transforms_at_indices!(all_params, ws.td_param_indices,
        sp.td_transforms, first_context)
    _ctsem_td_impulse!(ws, pars, first_context.tdpreds)
    update_context = CTSEMRowContext(ws.state, pars, first_context.tdpreds,
        tipreds, timesteps[1], zero(T), subject, 1)
    apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
        sp.update_transforms, update_context)

    gen === nothing || _ctsem_record_state!(gen, ws, 1, n)
    total = _ctsem_state_row!(ws, pars, data, 1, pred, theta, chol, resid,
        scratch, gaussian, log2pi, gen)
    if gen === nothing
        isfinite(total) || return _ctsem_invalid(T)
    else
        gen.llrow[gen.offset + 1] = total
    end

    prev = timesteps[1]
    nsteps = min(length(timesteps), size(data, 2))
    @inbounds for t in 2:nsteps
        dt = timesteps[t] - prev
        nsub = _ctsem_substeps(dt, max_timestep, t)
        substep_dt = dt / nsub
        for substep in 1:nsub
            substep_time = prev + substep * substep_dt
            predict_context = CTSEMRowContext(ws.state, pars,
                view(tdpreds, :, t), tipreds, substep_time, substep_dt,
                subject, t)
            apply_complex_transforms_at_indices!(all_params,
                ws.predict_param_indices, sp.predict_transforms, predict_context)

            ContinuousTimeSEM.sdcovsqrt2cov!(ws.bufferQ, pars.DIFFUSION, 0,
                ws.state_dim)
            if ws.continuous_time
                _compute_discrete_time_form!(ws.discrete_ca, ws.bufferQ,
                    ws.bufferQ.out, pars, substep_dt, ws.exp_buffer,
                    ws.lyap_buffer, ws.state, indices, ws.diffusion_buffer,
                    ws.discretization_buffer, ws.state_dim,
                    ws.discretization_cache)
            else
                _compute_one_step_form!(ws.discrete_ca, ws.bufferQ.out, pars,
                    ws.state, indices, ws.state_dim)
            end

            # Deterministic part first, into scratch, because the noise reads
            # the state it is added to only after the whole product is formed.
            _matvec_mul!(ws.bufferQ.r, ws.discrete_ca.dDRIFT, ws.state,
                ws.state_dim, ws.state_dim)
            for i in 1:n
                ws.state[i] = ws.bufferQ.r[i] + ws.discrete_ca.dINT[i]
            end

            # Process noise, over the diffusing states only. Everything else is
            # a static coordinate whose covariance the filter propagates through
            # the transition alone, and which correspondingly gets no
            # innovation of its own here.
            for j in 1:k, i in 1:k
                qfactor[i, j] = ws.discrete_ca.dDIFFUSION[indices[i], indices[j]]
            end
            _ctsem_lower_chol!(qfactor, k)
            for i in 1:k
                acc = zero(T)
                for j in 1:i
                    acc += qfactor[i, j] * T(z[at + j])
                end
                ws.state[indices[i]] += acc
            end
            at += k
        end

        td_context = CTSEMRowContext(ws.state, pars, view(tdpreds, :, t),
            tipreds, timesteps[t], dt, subject, t)
        apply_complex_transforms_at_indices!(all_params, ws.td_param_indices,
            sp.td_transforms, td_context)
        _ctsem_td_impulse!(ws, pars, td_context.tdpreds)

        measurement_context = CTSEMRowContext(ws.state, pars, td_context.tdpreds,
            tipreds, timesteps[t], dt, subject, t)
        apply_complex_transforms_at_indices!(all_params, ws.update_param_indices,
            sp.update_transforms, measurement_context)

        gen === nothing || _ctsem_record_state!(gen, ws, t, n)
        row_ll = _ctsem_state_row!(ws, pars, data, t, pred, theta, chol, resid,
            scratch, gaussian, log2pi, gen)
        if gen === nothing
            isfinite(row_ll) || return _ctsem_invalid(T)
        else
            gen.llrow[gen.offset + t] = row_ll
        end
        total += row_ll
        prev = timesteps[t]
    end
    return total
end

@inline function _ctsem_record_state!(gen::CTSEMStateGenerate, ws, col::Int,
    n::Int)
    @inbounds for i in 1:n
        gen.states[i, gen.offset + col] = Float64(ws.state[i])
    end
    return nothing
end


################################################################################
# The joint density
################################################################################

"""
    ctsem_joint_loglikelihood(objective, values, z)

`log p(y, z | theta)`: the observations given the states the innovations `z`
build, plus the standard normal density of `z` itself, plus whatever prior the
objective carries on `theta`.

This is a proper joint density, not a profile or a penalised likelihood: it
integrates to the marginal likelihood over `z`, which is what makes it usable
as a sampling target and what makes a mode of it interpretable. Nothing is
approximated anywhere in it -- the Gaussian assumption the EKF makes about the
state posterior has no counterpart here.

`values` and `z` are promoted together, so a gradient may be taken with respect
to either or both.

Serial over subjects, unlike the filter objective. Each subject's innovations
are a disjoint block and each already owns its workspace, so the same
cost-weighted chunking would apply unchanged -- it is left out because nothing
calls this in a loop yet. A fit over this target is what would want it.
"""
function ctsem_joint_loglikelihood(objective::CTSEMObjective,
    values::AbstractVector, z::AbstractVector)
    T = promote_type(eltype(values), eltype(z), Float64)
    parameters = collect(T, values)
    innovations = collect(T, z)
    layout = _ctsem_state_layout(objective)
    length(innovations) == layout.ndim || throw(DimensionMismatch(
        "z must have $(layout.ndim) entries for this design, got $(length(innovations))"))

    subjects = objective.subject_objectives
    sp = objective.params
    total = zero(T)
    @inbounds for (i, sub) in enumerate(subjects)
        ws = _get_or_init_objective_workspace!(sub, T)
        contribution = _ctsem_state_pass!(ws, parameters, sub.data,
            sub.timesteps, sp, sub.tdpreds, sub.tipreds, sub.subject,
            sub.max_timestep, innovations, layout.zoffsets[i], nothing)
        isfinite(contribution) || return _ctsem_invalid(T)
        total += contribution
    end

    # The innovations' own density. Written out rather than folded into the
    # subject loop so that a subject whose pass fails cannot leave a partial
    # prior behind.
    prior = zero(T)
    @inbounds for value in innovations
        prior -= value * value / 2
    end
    prior -= layout.ndim * log(2 * T(pi)) / 2
    return total + prior + _ctsem_log_prior(objective, parameters)
end

"""
    ctsem_joint_evaluate(objective, values, z; gradient=true)

The joint density and, optionally, its gradient with respect to `[values; z]`.

Forward mode, because there is no adjoint for this pass yet and because the
gradient a sampler needs is with respect to a vector whose length is dominated
by `z` -- for which reverse mode is the only sensible answer eventually. The
value alone is cheap and is what generation and testing use.

The pieces are returned alongside the total: `observation` is the conditional
log likelihood of the data, `state_prior` the innovations' own density, and
`parameter_prior` whatever prior the objective carries. They sum to `value`.
"""
function ctsem_joint_evaluate(objective::CTSEMObjective, values::AbstractVector,
    z::AbstractVector; gradient::Bool=true)
    npar = length(values)
    parameters = collect(Float64, values)
    innovations = collect(Float64, z)
    packed = vcat(parameters, innovations)
    target = x -> ctsem_joint_loglikelihood(objective, view(x, 1:npar),
        view(x, (npar + 1):length(x)))
    value = target(packed)
    grad = gradient ? ForwardDiff.gradient(target, packed) : nothing

    layout = _ctsem_state_layout(objective)
    state_prior = -sum(abs2, innovations) / 2 - layout.ndim * log(2pi) / 2
    parameter_prior = _ctsem_log_prior(objective, parameters)
    return (value=value, gradient=grad,
        observation=value - state_prior - parameter_prior,
        state_prior=state_prior, parameter_prior=parameter_prior,
        ndim=layout.ndim, npar=npar)
end


################################################################################
# Generation
################################################################################

"""
    ctsem_generate_states(objective, values, z, base)

One dataset drawn from the model itself: the trajectory the innovations `z`
build, and then each observation from its conditional distribution given the
state at its row.

`z` is `ctsem_state_dimension(objective)` standard normals and `base` is
`nmanifest` by `nrows` more, one per manifest cell. Returns the generated
observations `Y` (`NaN` wherever the input had no observation), the sampled
`states`, each row's conditional log likelihood `llrow`, and the per-subject
totals.

# How this differs from `ctsem_generate`

That one draws each row from the filter's one-step-ahead predictive and then
lets the filter condition on what it drew, which makes the dataset an exact
draw from the density the filter maximises. Here the state is drawn from the
process and the observation from the state, which makes the dataset a draw from
the *model* -- and for a Gaussian model those are the same thing, since the
filter's predictive is then exact.

For a non-Gaussian one they are not, and the difference is not a subtlety. The
filter's categorical update is an assumed-density projection: it can move the
state a long way on one improbable observation, and with an unbounded indicator
the row after that is drawn from a rate that has already moved. Nothing of the
kind can happen here, because no observation ever moves a state.
"""
function ctsem_generate_states(objective::CTSEMObjective,
    values::AbstractVector, z::AbstractVector, base::AbstractMatrix)

    sp = objective.params
    parameters = collect(Float64, values)
    innovations = collect(Float64, z)
    layout = _ctsem_state_layout(objective)
    length(innovations) == layout.ndim || throw(DimensionMismatch(
        "z must have $(layout.ndim) entries for this design, got $(length(innovations))"))
    size(base) == (layout.nmanifest, layout.nrows) || throw(DimensionMismatch(
        "base must be $(layout.nmanifest) by $(layout.nrows)"))

    gen = CTSEMStateGenerate(Matrix{Float64}(base),
        fill(NaN, layout.nmanifest, layout.nrows),
        zeros(Float64, layout.nlatent, layout.nrows),
        fill(NaN, layout.nrows))
    subjects = objective.subject_objectives
    loglik = zeros(Float64, length(subjects))

    for (i, sub) in enumerate(subjects)
        gen.offset = layout.rowoffsets[i]
        ws = _get_or_init_objective_workspace!(sub, Float64)
        loglik[i] = _ctsem_state_pass!(ws, parameters, sub.data, sub.timesteps,
            sp, sub.tdpreds, sub.tipreds, sub.subject, sub.max_timestep,
            innovations, layout.zoffsets[i], gen)
    end
    return (Y=gen.out, states=gen.states, llrow=gen.llrow,
        subject_loglik=loglik)
end


################################################################################
# Fitting over the joint density
################################################################################
#
# The optimiser and the sampler both want the same thing from an objective: a
# value and a gradient at a vector. `CTSEMJointObjective` supplies them over
# `x = [theta; z]`, so `ctsem_optimize` and `ctsem_sample_marginal` drive the
# state-explicit target with no changes at all -- the first is typed on
# `CTSEMOptimisable`, which this is, and the second was already duck-typed.
#
# What that buys, for `optimize=FALSE`, is the estimator this path is actually
# for: NUTS over parameters *and* states, which is the exact posterior with no
# Gaussian assumption about the state anywhere. `optimize=TRUE` gives the joint
# mode instead -- the same thing Stan's optimiser does at `intoverstates=0`,
# and carrying the same caveat, which ctFit already states: maximising over
# `z` rather than integrating it out biases the variance parameters downward,
# because a variance whose realisations are also being chosen can always be
# made to look smaller.

export CTSEMJointObjective, ctsem_joint_objective, ctsem_joint_blocks,
    ctsem_joint_hessian, ctsem_joint_states, ctsem_joint_dimension

"""
    CTSEMJointObjective(objective, npar)

The joint density as something to optimise or sample: `x = [theta; z]`, the
`npar` population parameters first and then every subject's innovations in
subject order.

`npar` is supplied rather than inferred. The engine cannot infer it -- the raw
vector's length is decided by ctsem's parameter layout on the R side, which
counts population scales and TI-predictor coefficients the parameter table
alone does not reach -- and inferring it from the largest `parnumber` would
silently shorten the vector for exactly the models where that is wrong.
"""
struct CTSEMJointObjective{O} <: CTSEMOptimisable
    objective::O
    layout::CTSEMStateLayout
    npar::Int
end

function ctsem_joint_objective(objective::CTSEMObjective, npar::Integer)
    npar = Int(npar)
    npar >= 0 || throw(ArgumentError("npar must be non-negative"))
    return CTSEMJointObjective(objective, _ctsem_state_layout(objective), npar)
end

"""Total dimension optimised or sampled: parameters plus innovations."""
ctsem_joint_dimension(o::CTSEMJointObjective) = o.npar + o.layout.ndim

"""
    ctsem_joint_blocks(objective)

Where each block of `x` lives: the population parameters, then one block per
subject.

The joint Hessian is arrow-shaped -- dense in the population block, block
diagonal in the innovations, coupled between the two -- because a subject's
innovations enter no other subject's likelihood. Everything below exploits
that, and a caller building a sampler metric wants the same partition.
"""
function ctsem_joint_blocks(o::CTSEMJointObjective)
    layout = o.layout
    nsubjects = length(layout.zoffsets)
    starts = Vector{Int}(undef, nsubjects)
    stops = Vector{Int}(undef, nsubjects)
    for i in 1:nsubjects
        starts[i] = o.npar + layout.zoffsets[i] + 1
        stops[i] = i < nsubjects ? o.npar + layout.zoffsets[i + 1] :
            o.npar + layout.ndim
    end
    return (npar=o.npar, ndim=layout.ndim, start=starts, stop=stops)
end

"""One subject's slice of the innovation vector, as a range into `x`."""
@inline function _ctsem_joint_range(o::CTSEMJointObjective, i::Int)
    layout = o.layout
    first = o.npar + layout.zoffsets[i] + 1
    last = i < length(layout.zoffsets) ? o.npar + layout.zoffsets[i + 1] :
        o.npar + layout.ndim
    return first:last
end

"""
    _ctsem_joint_subject(o, i, theta, z, offset)

One subject's contribution: the conditional log likelihood of its observations
plus the density of its own innovations.

The innovation prior is split across subjects rather than added once at the
end, so a subject's term is self-contained and its gradient can be taken
against `[theta; z_i]` alone.
"""
function _ctsem_joint_subject(o::CTSEMJointObjective, i::Int,
    theta::AbstractVector{T}, z::AbstractVector{T}, offset::Int) where {T}
    sub = o.objective.subject_objectives[i]
    ws = _get_or_init_objective_workspace!(sub, T)
    value = _ctsem_state_pass!(ws, theta, sub.data, sub.timesteps,
        o.objective.params, sub.tdpreds, sub.tipreds, sub.subject,
        sub.max_timestep, z, offset, nothing)
    isfinite(value) || return value
    prior = zero(T)
    @inbounds for k in eachindex(z)
        prior -= z[k] * z[k] / 2
    end
    return value + prior - length(z) * log(2 * T(pi)) / 2
end

function (o::CTSEMJointObjective)(x::AbstractVector{T}) where {T}
    length(x) == ctsem_joint_dimension(o) || throw(DimensionMismatch(
        "x must have $(ctsem_joint_dimension(o)) entries, got $(length(x))"))
    theta = view(x, 1:o.npar)
    total = zero(T)
    @inbounds for i in eachindex(o.layout.zoffsets)
        range = _ctsem_joint_range(o, i)
        contribution = _ctsem_joint_subject(o, i, theta, view(x, range), 0)
        isfinite(contribution) || return _ctsem_invalid(T)
        total += contribution
    end
    return total + _ctsem_log_prior(o.objective, theta)
end

"""
    ctsem_evaluate(objective::CTSEMJointObjective, x; gradient=true, ...)

The joint density and its gradient with respect to `[theta; z]`.

Forward mode, but *per subject*, over `[theta; z_i]` rather than over the whole
vector at once. Both are exact and they differ by an enormous constant factor.
A subject's innovations enter no other subject's likelihood, so a dual pass
seeded with the whole vector carries `dim(z)` partials through every subject to
compute the `dim(z_i)` of them that are not structurally zero. On thirty
subjects with twenty innovations each and ten parameters that is 610 seeded
partials through thirty subjects, against thirty passes of thirty seeded
partials through one subject -- twenty times less arithmetic, and the ratio
grows with the subject count without bound.

`gradient_method` is accepted and ignored: there is no reverse pass for this
target. The filter's adjoint is built around the covariance recursion, which
this path does not have.
"""
function ctsem_evaluate(o::CTSEMJointObjective, x::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:forward)
    ndim = ctsem_joint_dimension(o)
    length(x) == ndim || throw(DimensionMismatch(
        "x must have $(ndim) entries, got $(length(x))"))
    values = collect(Float64, x)
    theta = collect(Float64, view(values, 1:o.npar))
    nsubjects = length(o.layout.zoffsets)
    subject = zeros(Float64, nsubjects)
    grad = gradient ? zeros(Float64, ndim) : nothing
    total = 0.0
    valid = true

    for i in 1:nsubjects
        range = _ctsem_joint_range(o, i)
        if gradient
            # `[theta; z_i]`, so the population partials come back from the
            # same pass and are accumulated across subjects.
            packed = vcat(theta, values[range])
            target = function (v)
                return _ctsem_joint_subject(o, i, view(v, 1:o.npar),
                    view(v, (o.npar + 1):length(v)), 0)
            end
            result = DiffResults.GradientResult(packed)
            result = ForwardDiff.gradient!(result, target, packed)
            value = DiffResults.value(result)
            partials = DiffResults.gradient(result)
            subject[i] = value
            if !isfinite(value) || !all(isfinite, partials)
                valid = false
                break
            end
            total += value
            @inbounds for k in 1:o.npar
                grad[k] += partials[k]
            end
            @inbounds for (position, k) in enumerate(range)
                grad[k] = partials[o.npar + position]
            end
        else
            value = _ctsem_joint_subject(o, i, theta, view(values, range), 0)
            subject[i] = value
            if !isfinite(value)
                valid = false
                break
            end
            total += value
        end
    end

    if !valid
        # The same contract `ctsem_optimize`'s `fg!` expects of the marginal
        # objective: a non-finite value is a rejected trial point, not an error.
        nan = _ctsem_invalid(Float64)
        gradient && fill!(grad, nan)
        return contributions ?
            (value=nan, gradient=grad, subject_loglik=subject,
                row_loglik=Float64[]) :
            (value=nan, gradient=grad)
    end

    total += _ctsem_log_prior(o.objective, theta)
    if gradient && !isempty(o.objective.prior_index)
        _ctsem_log_prior_gradient!(view(grad, 1:o.npar), o.objective, theta, 1.0)
    end
    contributions || return (value=total, gradient=grad)
    # Row likelihoods, for the summaries that read them. Conditional on the
    # sampled state, like `ctsem_generate_states`' own `llrow` and unlike the
    # filter's -- so they sum to the observation term rather than to a marginal
    # log likelihood, which is the honest thing for them to sum to here.
    return (value=total, gradient=grad, subject_loglik=subject,
        row_loglik=_ctsem_joint_record(o, values).llrow)
end

"""
    _ctsem_joint_record(o, x)

One recording pass at `x`: the states it implies and each row's conditional log
likelihood.

The generation record is used as a *recorder* here -- `_ctsem_state_row!`
writes `llrow` and `states` only on that path -- and the data it would draw is
written into a buffer that is thrown away. The alternative is a second copy of
the row loop that records instead of drawing, which is the kind of duplicate
that drifts.
"""
function _ctsem_joint_record(o::CTSEMJointObjective, x::AbstractVector)
    layout = o.layout
    theta = collect(Float64, view(x, 1:o.npar))
    innovations = collect(Float64, view(x, (o.npar + 1):length(x)))
    gen = CTSEMStateGenerate(zeros(Float64, layout.nmanifest, layout.nrows),
        fill(NaN, layout.nmanifest, layout.nrows),
        zeros(Float64, layout.nlatent, layout.nrows), fill(NaN, layout.nrows))
    for (i, sub) in enumerate(o.objective.subject_objectives)
        gen.offset = layout.rowoffsets[i]
        ws = _get_or_init_objective_workspace!(sub, Float64)
        _ctsem_state_pass!(ws, theta, sub.data, sub.timesteps,
            o.objective.params, sub.tdpreds, sub.tipreds, sub.subject,
            sub.max_timestep, innovations, layout.zoffsets[i], gen)
    end
    return gen
end

"""
    ctsem_joint_states(objective, x)

The latent states `x` implies, `nlatent` by `nrows`.

What a fit over this target has that a marginal one does not: the trajectory
itself, at the estimate, with no smoother needed to recover it.
"""
ctsem_joint_states(o::CTSEMJointObjective, x::AbstractVector) =
    _ctsem_joint_record(o, x).states

"""
    ctsem_joint_hessian(objective, x; profile=true)

The second derivatives of the joint density at `x`.

Built subject by subject, over `[theta; z_i]`, and scattered: the innovation
blocks of two different subjects are *structurally* zero, not merely small, so
a dense Hessian over the whole vector would spend almost all of its arithmetic
confirming that.

With `profile=true` (the default) the answer is the `npar` by `npar` Schur
complement

    H_pp - H_pz inv(H_zz) H_zp

rather than the whole matrix. That is the curvature of the *profile* function
`max_z log p(y, z | theta)`, and it is strictly better than reading `H_pp`
alone, which reports the curvature at fixed states and so treats a trajectory
that was estimated as though it had been observed.

It is **not** the observed information for `theta`, and must not be presented
as one. The Laplace approximation to the marginal adds `-log det(-H_zz)/2`,
whose own curvature in `theta` this omits -- and that omitted term is what
identifies the parameters. With roughly as many innovations as observations the
profile is nearly flat without it: measured on fifteen subjects of six Gaussian
rows, four free parameters, the largest eigenvalue of the profiled curvature was
0.05, against 22.8 for a single well-determined parameter in a count model. The
R side does not report standard errors from an optimised state-explicit fit for
that reason; sampling the same density is what gives usable intervals.

The inner block is inverted per subject, which is what makes this affordable:
`H_zz` is block diagonal by subject, so the complement is a sum of per-subject
terms and no matrix of the full innovation dimension is ever formed.

Returns `nothing` if any subject's inner block is not negative definite, which
a caller must treat as "no Hessian available" rather than as a number.
"""
function ctsem_joint_hessian(o::CTSEMJointObjective, x::AbstractVector;
    profile::Bool=true)
    ndim = ctsem_joint_dimension(o)
    length(x) == ndim || throw(DimensionMismatch(
        "x must have $(ndim) entries, got $(length(x))"))
    values = collect(Float64, x)
    theta = collect(Float64, view(values, 1:o.npar))
    nsubjects = length(o.layout.zoffsets)
    npar = o.npar

    full = profile ? zeros(Float64, npar, npar) : zeros(Float64, ndim, ndim)
    for i in 1:nsubjects
        range = _ctsem_joint_range(o, i)
        k = length(range)
        packed = vcat(theta, values[range])
        target = function (v)
            return _ctsem_joint_subject(o, i, view(v, 1:npar),
                view(v, (npar + 1):length(v)), 0)
        end
        block = ForwardDiff.hessian(target, packed)
        block = (block .+ transpose(block)) ./ 2
        all(isfinite, block) || return nothing
        Hpp = view(block, 1:npar, 1:npar)
        Hpz = view(block, 1:npar, (npar + 1):(npar + k))
        Hzz = view(block, (npar + 1):(npar + k), (npar + 1):(npar + k))
        if profile
            # `-Hzz` is positive definite at a maximum, so the solve is a
            # Cholesky on it. A subject whose inner curvature is singular has a
            # direction of its trajectory the data does not determine, and no
            # arithmetic here produces a standard error from that.
            factor = cholesky(Symmetric(-Matrix(Hzz)), check=false)
            issuccess(factor) || return nothing
            full .+= Hpp .+ Hpz * (factor \ transpose(Hpz))
        else
            full[1:npar, 1:npar] .+= Hpp
            full[1:npar, range] .= Hpz
            full[range, 1:npar] .= transpose(Hpz)
            full[range, range] .= Hzz
        end
    end
    # The prior is over the parameters alone and is not a per-subject term, so
    # it is added once. Its second derivative is `-weight / scale^2` on the
    # diagonal, matching `_ctsem_log_prior`.
    @inbounds for (position, index) in enumerate(o.objective.prior_index)
        1 <= index <= npar || continue
        scale = o.objective.prior_scale[position]
        full[index, index] -= o.objective.prior_weight / (scale * scale)
    end
    return full
end

"""
`ctsem_hessian` for the joint target: the whole arrow-shaped matrix, which is
what a sampler metric over `[theta; z]` needs.

`ctsem_joint_hessian(o, x)` is the profiled `npar` by `npar` form, which is
what a standard error needs. They are different objects for different callers,
and neither is a default for the other.
"""
function ctsem_hessian(o::CTSEMJointObjective, values::AbstractVector;
    chunk::Integer=0)
    result = ctsem_joint_hessian(o, values; profile=false)
    result === nothing && throw(ArgumentError(
        "the joint Hessian was not finite at this point"))
    return result
end

"""
Only the population block can saturate: the innovations carry no transform.

See `_ctsem_saturation_range` in `ctsem_backend.jl` for what the check is for.
Without this an ordinary trajectory excursion -- five standard deviations is
one row in a few hundred -- would be reported as a parameter pinned at the
floating-point limit of its transform, and the fit declared not converged.
"""
_ctsem_saturation_range(o::CTSEMJointObjective, minimizer) = 1:o.npar

"""
Unwrap to the `CTSEMObjective` the joint objective is built on -- see
`_ctsem_params` in `ctsem_backend.jl`. The population block's raw coordinates
are exactly that objective's, so no separate lookup is needed.
"""
_ctsem_params(o::CTSEMJointObjective) = _ctsem_params(o.objective)
