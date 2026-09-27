################################################################################
# How much of each random effect a subject's own data determine
################################################################################
#
# A random effect is only as good as the data behind each subject's value of
# it. When each subject contributes a handful of observations of a quantity
# that a handful of observations cannot pin down -- a subject's own rate of
# change from three waves, say -- the population sd of that effect is barely
# determined, the likelihood of many subjects goes convex in it, and the
# Laplace term over-credits them. Every symptom of that arrives at the end of
# a fit, and obliquely: a flat direction in the information matrix, a unit
# count in `ctsem_laplace_conditioning`. What says it directly, and before the
# fit, is the fraction of each subject's effect its own data determine,
#
#     1 - posterior variance / population variance,
#
# the complement of what population pharmacokinetics calls shrinkage. Zero is
# a subject whose data say nothing about its effect; one is a subject whose
# data determine it outright. It can go negative, where the likelihood is
# convex in the effect and the posterior is wider than the population.
#
# Both routes compute the same quantity, each from the object its own
# objective already builds, at a point the caller names. One pass over the
# subjects, serial, and no LAPACK: the Laplace side reuses the unit curvature
# and its block factorization, and the augmented side's backward pass solves
# through the engine's own Cholesky.

"""
    ctsem_effect_information(objective, values)

How much of each random effect each group's own data determine, at the raw
parameter vector `values`: `1 - posterior variance / population variance` of
the effect, per group and effect.

On the Laplace route (`CTSEMLaplaceObjective`) the posterior covariance of a
unit's standardised effects is the inverse of its inner curvature `M` at the
mode -- the matrix whose log determinant the Laplace term takes -- so the
posterior variance of effect `j` on the raw scale is `[L M^-1 L']_jj`, against
the population variance `[L L']_jj`. Returns, per level, `determined`
(`ngroups` by effects) and `popvar` (one per effect); per unit, the smallest
eigenvalue of `M` as `ctsem_laplace_conditioning` reports it and whether the
inner solve converged. The modes are solved at `values` from the origin, as an
evaluation solves them.

On the augmented route (`CTSEMObjective`) each effect is a state of the
filter (`population_indices`), its population variance the variance the
filter starts it at, and its posterior variance that state's smoothed
variance at the subject's first row: a carrier state has no dynamics, so that
is its variance given every observation, and for an individually varying
initial state it is the initial state's. Returns `determined` and `priorvar`,
subjects by effects.

A level or route with no random effects returns a one-by-one `NaN`
placeholder rather than an empty array, which the R bridge cannot carry.
"""
function ctsem_effect_information(laplace::CTSEMLaplaceObjective,
    values::AbstractVector)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    spec = laplace.spec
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, spec)
    nlev = length(spec.levels)
    determined = Vector{Matrix{Float64}}(undef, nlev)
    popvar = Vector{Vector{Float64}}(undef, nlev)
    filled = Vector{Vector{Bool}}(undef, nlev)
    for l in 1:nlev
        lv = spec.levels[l]
        k = nrandomeffects(lv)
        L = Ls[l]
        determined[l] = fill(NaN, max(lv.ngroups, 1), max(k, 1))
        popvar[l] = k == 0 ? [NaN] : [sum(abs2, view(L, j, :)) for j in 1:k]
        filled[l] = fill(false, lv.ngroups)
    end
    nunits = length(laplace.units.members)
    mins = fill(Inf, nunits)
    converged = fill(true, nunits)
    # The modes below are solved at `theta`, so the diagnostics that read
    # `last_values` describe the point they were left at.
    laplace.last_values = theta
    for U in 1:nunits
        d = laplace.units.dims[U]
        d == 0 && continue
        converged[U] = _laplace_solve_unit_mode!(laplace, U, theta, Ls).converged
        u = laplace.modes[U]
        blocks = laplace.units.blocks[U]
        M = _laplace_unit_curvature(laplace, U, theta, Ls, u)
        if !all(b -> all(isfinite, b), M.diag)
            mins[U] = NaN
            continue
        end
        # The smallest eigenvalue exactly as `ctsem_laplace_conditioning`
        # takes it, and before the factorization below, which shifts the
        # diagonal of a curvature it has to repair.
        if !_laplace_exceeds_identity(M, blocks)
            mins[U] = d > _LAPLACE_EIGEN_MAXDIM[] ? NaN :
                _ctsem_symeig(_laplace_block_dense(M, blocks, d)).values[1]
        end
        factored = _laplace_factor_repaired!(M, blocks)
        factored.ok || continue
        e = zeros(Float64, d)
        for (m, i) in enumerate(laplace.units.members[U])
            for l in 1:nlev
                lv = spec.levels[l]
                k = nrandomeffects(lv)
                L = Ls[l]
                r = size(L, 2)
                (k == 0 || r == 0) && continue
                # An outer group's effect is one vector shared by its members,
                # so it is computed once.
                g = lv.group[i]
                filled[l][g] && continue
                filled[l][g] = true
                base = laplace.units.offsets[U][m][l]
                # The group's block of M^-1, a column at a time, as
                # `ctsem_laplace_modes` takes it: never the whole inverse.
                block = zeros(Float64, r, r)
                for t in 1:r
                    fill!(e, 0.0)
                    e[base + t] = 1.0
                    x = _laplace_block_solve(factored.factors, factored.coupling,
                        blocks, e)
                    for q in 1:r
                        block[q, t] = x[base + q]
                    end
                end
                for j in 1:k
                    posterior = 0.0
                    for a in 1:r, b in 1:r
                        posterior += L[j, a] * block[a, b] * L[j, b]
                    end
                    population = popvar[l][j]
                    determined[l][g, j] = population > 0 ?
                        1 - posterior / population : NaN
                end
            end
        end
    end
    return (determined=determined, popvar=popvar, min_eigenvalue=mins,
        converged=converged,
        nrandom=[nrandomeffects(lv) for lv in spec.levels])
end

function ctsem_effect_information(objective::CTSEMObjective, values::AbstractVector)
    sp = objective.params
    states = sp.population_indices
    k = length(states)
    subjects = objective.subject_objectives
    nsub = length(subjects)
    determined = fill(NaN, max(nsub, 1), max(k, 1))
    priorvar = fill(NaN, max(nsub, 1), max(k, 1))
    (k == 0 || nsub == 0) &&
        return (determined=determined, priorvar=priorvar, nrandom=k)
    raw = collect(Float64, values)
    ws = _init_continuous_ekf_workspace(Float64, sp)
    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)
    # One subject's rows at a time, so the trace is sized to the longest
    # subject rather than to the data.
    maxrows = maximum(size(sub.data, 2) for sub in subjects)
    trace = CTSEMKalmanTrace(Float64, n, m, maxrows, 1, length(sp.mutables))
    for (i, sub) in enumerate(subjects)
        nobs = size(sub.data, 2)
        trace.offset = 0
        trace.current_subject = 1
        value = _extended_kalman_filter_continuous!(ws, raw, sub.data,
            collect(sub.timesteps), sp, sub.tdpreds, sub.tipreds, i,
            sub.max_timestep, trace)
        isfinite(value) || continue
        smoothed = _effect_smoothed_initial_cov(trace, nobs, n)
        for (j, s) in enumerate(states)
            prior = trace.etacov[_CTSEM_KALMAN_PRIOR, 1, s, s]
            priorvar[i, j] = prior
            determined[i, j] = prior > 0 ? 1 - smoothed[s, s] / prior : NaN
        end
    end
    return (determined=determined, priorvar=priorvar, nrandom=k)
end

export ctsem_effect_information

"""
    _effect_smoothed_initial_cov(trace, nobs, n)

The smoothed state covariance at a subject's first row: the covariance half of
the backward pass `_kalman_smooth!` makes, with the same ridge on the matrix it
solves against, but through the engine's own Cholesky rather than LAPACK's, and
without the means and the measurement quantities nothing here reads.
"""
function _effect_smoothed_initial_cov(trace::CTSEMKalmanTrace, nobs::Int, n::Int)
    smoothed = trace.etacov[_CTSEM_KALMAN_UPD, nobs, :, :]
    P = Matrix{Float64}(undef, n, n)
    gain = Matrix{Float64}(undef, n, n)
    scratch = Matrix{Float64}(undef, n, n)
    for r in (nobs - 1):-1:1
        prior = trace.etacov[_CTSEM_KALMAN_PRIOR, r + 1, :, :]
        updated = trace.etacov[_CTSEM_KALMAN_UPD, r, :, :]
        copyto!(P, prior)
        _symmetrize_and_ridge!(P, n)
        factor = _ctsem_cholesky(P, n)
        issuccess(factor) || return fill(NaN, n, n)
        # gain = P_upd A' P_prior^-1
        _ctsem_mulNT!(gain, updated, trace.transition[r + 1])
        rdiv!(gain, factor)
        # P_sm[r] = P_upd + gain (P_sm[r+1] - P_prior[r+1]) gain'
        smoothed .-= prior
        _ctsem_mul!(scratch, gain, smoothed)
        _ctsem_mulNT!(updated, scratch, gain, true, true)
        smoothed = updated
    end
    return smoothed
end
