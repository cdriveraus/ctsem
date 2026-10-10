# Per-row scores and conditional effect likelihoods, for specification tests.
#
# Two quantities the R side's ctScoreTest() and ctFitAdequacy() test with, and
# that nothing on the fitting path computes:
#
#   ctsem_row_scores             each data row's score: the derivative of its
#                                log predictive density (`llrow`) along a set
#                                of directions in the raw parameters. Under the
#                                model the row scores are a martingale
#                                difference sequence within each subject, which
#                                is what the dynamic information matrix test
#                                (White 1987) and the parameter-stability
#                                tests (Nyblom 1989) read.
#   ctsem_laplace_effect_loglik  each unit's log likelihood with its random
#                                effects fixed at given values, the deviations
#                                its filter integrates integrated exactly: the
#                                integrand the gradient function of Verbeke and
#                                Molenberghs (2013) compares with the marginal.
#
# Neither touches a fitting path. The derivatives are forward mode, seeded in
# at most `_ROW_SCORE_WIDTH` directions a pass, so one dual type is compiled
# for a model shape whatever the number of directions; a call asking for more
# makes more passes rather than a wider type.

struct _RowScoreTag end

"""The partials a row-score pass carries: directions are taken this many at a
time, and a narrower set is padded with zero directions, so a model compiles
one dual type for any number of them."""
const _ROW_SCORE_WIDTH = 8

const _RowScoreDual = ForwardDiff.Dual{_RowScoreTag,Float64,_ROW_SCORE_WIDTH}

"""The raw parameters as duals seeded in `directions[:, cols]`, padded with
zero directions to the pass's width."""
function _row_score_seed(theta::Vector{Float64}, directions::AbstractMatrix,
    cols::UnitRange{Int})
    k = length(cols)
    return [_RowScoreDual(theta[j], ForwardDiff.Partials{_ROW_SCORE_WIDTH,Float64}(
        ntuple(a -> a <= k ? Float64(directions[j, cols[a]]) : 0.0,
            Val(_ROW_SCORE_WIDTH)))) for j in eachindex(theta)]
end

"""
    ctsem_row_scores(objective, values, directions)

Each data row's score along each column of `directions` (`npar` by `K`): the
directional derivative of the row's log one-step-ahead predictive density at
`values`, rows in the order `ctsem_kalman` reports them. A subject's rows sum
to the directional derivative of its log likelihood. Rows the filter cannot
evaluate, and every row of a subject whose filter fails, are `NaN`.
"""
function ctsem_row_scores(objective::CTSEMObjective, values::AbstractVector,
    directions::AbstractMatrix)
    theta = collect(Float64, values)
    size(directions, 1) == length(theta) || throw(DimensionMismatch(
        "directions needs one row per raw parameter"))
    K = size(directions, 2)
    sp = objective.params
    subjects = objective.subject_objectives
    rows = [size(s.data, 2) for s in subjects]
    out = fill(NaN, sum(rows; init=0), K)
    K == 0 && return out
    ws = _init_continuous_ekf_workspace(_RowScoreDual, sp)
    ws.hyb.mode = _HYBRID_OFF
    n = _val(ws.state_dim)
    m = _val(ws.manifest_dim)
    nall = length(sp.mutables)
    for c0 in 1:_ROW_SCORE_WIDTH:K
        cols = c0:min(c0 + _ROW_SCORE_WIDTH - 1, K)
        x = _row_score_seed(theta, directions, cols)
        offset = 0
        for (i, sub) in enumerate(subjects)
            nobs = rows[i]
            trace = CTSEMKalmanTrace(_RowScoreDual, n, m, nobs, 1, nall)
            trace.offset = 0
            trace.current_subject = 1
            ll = _extended_kalman_filter_continuous!(ws, x, sub.data,
                collect(sub.timesteps), sp, sub.tdpreds,
                _ctsem_tipred_vector(sub.tipreds, x), i, sub.max_timestep, trace)
            if isfinite(ForwardDiff.value(ll))
                @inbounds for r in 1:nobs, (a, col) in enumerate(cols)
                    out[offset + r, col] = ForwardDiff.partials(trace.llrow[r], a)
                end
            end
            offset += nobs
        end
    end
    return out
end

"""
    ctsem_row_scores(laplace, values, directions)

The same for a Laplace fit, from the predictions `ctsem_laplace_sequential`
makes: a row's log predictive density is the difference of the subject's
Laplace terms on its record up to the row and up to the row before, each with
the effects' mode solved on that truncated record, so its score is the
difference of those terms' derivatives. Each term's derivative is the dual
term at the dual mode (`_laplace_dual_unit_mode`, the implicit function
theorem by one Newton step), seeded in the directions. Single-level fits only,
as for the sequential predictions.
"""
function ctsem_row_scores(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    directions::AbstractMatrix)
    theta = collect(Float64, values)
    size(directions, 1) == length(theta) || throw(DimensionMismatch(
        "directions needs one row per raw parameter"))
    units = laplace.units
    all(length(mm) == 1 for mm in units.members) || throw(ArgumentError(
        "row scores need every subject to be its own unit: a single level of " *
        "random effects"))
    K = size(directions, 2)
    subjects = laplace.objective.subject_objectives
    rows = [size(s.data, 2) for s in subjects]
    first_row = cumsum(vcat(0, rows[1:(end - 1)]))
    out = fill(NaN, sum(rows; init=0), K)
    K == 0 && return out
    Ls = _laplace_popchols(theta, laplace.spec)
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, eachindex(units.members)) do U
        _laplace_row_scores_unit!(out, laplace, U, theta, Ls, directions,
            first_row[units.members[U][1]])
        return true
    end
    return out
end

function _laplace_row_scores_unit!(out, laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, directions::AbstractMatrix,
    offset::Integer)
    spec = laplace.spec
    units = laplace.units
    i = units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    full = (data=so.data, timesteps=so.timesteps, tdpreds=so.tdpreds,
        max_timestep=so.max_timestep)
    T = size(full.data, 2)
    K = size(directions, 2)
    saved = (mode=copy(laplace.modes[U]), iterations=laplace.inner_iterations[U],
        gradient=laplace.inner_gradient[U], converged=laplace.inner_converged[U],
        hrepaired=laplace.hessian_repaired[U], mrepaired=laplace.mode_repaired[U],
        floored=laplace.logdet_floored[U])
    npar = length(theta)
    blocks = units.blocks[U]
    awsd = _laplace_workspace!(laplace, _RowScoreDual, npar)
    previous = zeros(Float64, K)
    current = zeros(Float64, K)
    try
        for t in 1:T
            _laplace_truncate!(so, full, t)
            _laplace_solve_unit_mode!(laplace, U, theta, Ls)
            u = copy(laplace.modes[U])
            M = isempty(u) ? CTSEMBlockMatrix(Float64, blocks) :
                _laplace_unit_curvature(laplace, U, theta, Ls, u)
            fac = _laplace_factor_repaired!(M, blocks)
            curvature = (fac.factors, fac.coupling)
            ok = true
            for c0 in 1:_ROW_SCORE_WIDTH:K
                cols = c0:min(c0 + _ROW_SCORE_WIDTH - 1, K)
                x = _row_score_seed(theta, directions, cols)
                Lsd = _laplace_popchols(x, spec)
                ud = _laplace_dual_unit_mode(laplace, U, x, Lsd, u, curvature, awsd)
                term = _laplace_unit_term(laplace, U, x, Lsd, ud, awsd)
                if !isfinite(ForwardDiff.value(term))
                    ok = false
                    break
                end
                for (a, col) in enumerate(cols)
                    current[col] = ForwardDiff.partials(term, a)
                end
            end
            if ok
                @inbounds for col in 1:K
                    out[offset + t, col] = current[col] - previous[col]
                end
                copyto!(previous, current)
            else
                # The next row's difference has no defined start either.
                fill!(previous, NaN)
            end
        end
    finally
        so.data = full.data
        so.timesteps = full.timesteps
        so.tdpreds = full.tdpreds
        so.max_timestep = full.max_timestep
        laplace.modes[U] = saved.mode
        laplace.inner_iterations[U] = saved.iterations
        laplace.inner_gradient[U] = saved.gradient
        laplace.inner_converged[U] = saved.converged
        laplace.hessian_repaired[U] = saved.hrepaired
        laplace.mode_repaired[U] = saved.mrepaired
        laplace.logdet_floored[U] = saved.floored
    end
    return nothing
end

"""
    ctsem_laplace_effect_loglik(laplace, values, points)

Each unit's log likelihood with its random effects fixed at each column of
`points` (`d` by `G`, in the unit's standardised coordinates, whose prior is
standard normal), with the deviations its members' filters integrate
(`CTSEMHybrid`) integrated exactly -- the integrand quadrature integrates
(`_laplace_member_loglik`). `nunits` by `G`. Single-level fits only: there a
unit is one subject and its coordinates are the level's.
"""
function ctsem_laplace_effect_loglik(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, points::AbstractMatrix)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    units = laplace.units
    all(length(mm) == 1 for mm in units.members) || throw(ArgumentError(
        "effect likelihoods need every subject to be its own unit: a single " *
        "level of random effects"))
    nunits = length(units.members)
    G = size(points, 2)
    out = fill(NaN, nunits, G)
    Ls = _laplace_popchols(theta, laplace.spec)
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, 1:nunits) do U
        d = units.dims[U]
        size(points, 1) == d || throw(DimensionMismatch(
            "points needs one row per random effect of a unit"))
        aws = _laplace_workspace!(laplace, Float64, length(theta))
        so = laplace.objective.subject_objectives[units.members[U][1]]
        u = zeros(Float64, d)
        for g in 1:G
            @inbounds for k in 1:d
                u[k] = points[k, g]
            end
            shifted = _laplace_member_values(theta, laplace.spec, Ls, u,
                units.offsets[U][1])
            out[U, g] = _laplace_member_loglik(aws, so, shifted)
        end
        return true
    end
    return out
end

"""
    ctsem_laplace_unit_dims(laplace)

Each unit's number of random-effect coordinates, the rows
`ctsem_laplace_effect_loglik` takes.
"""
ctsem_laplace_unit_dims(laplace::CTSEMLaplaceObjective) = copy(laplace.units.dims)

export ctsem_row_scores, ctsem_laplace_effect_loglik, ctsem_laplace_unit_dims
