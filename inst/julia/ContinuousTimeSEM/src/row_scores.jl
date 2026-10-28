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
#   ctsem_laplace_posterior_nodes  each unit's posterior over its standardised
#                                random effects, as the adaptive rule the
#                                quadrature integrates with: for a test of the
#                                effects' distribution, whose posterior
#                                expectations average to their prior ones under
#                                the model.
#
# Neither touches a fitting path. On the filter routes the derivatives are the
# hand-coded forward sweep's (`_tan_row_scores!`); on the Laplace route each
# truncated record's unit gradient, as the evaluation forms it.

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
    rows = [size(s.data, 2) for s in objective.subject_objectives]
    out = fill(NaN, sum(rows; init=0), K)
    K == 0 && return out
    _tan_row_scores!(out, objective, theta, directions, rows)
    return out
end

"""
    ctsem_row_scores(laplace, values, directions)

The same for a Laplace fit, from the predictions `ctsem_laplace_sequential`
makes: a row's log predictive density is the difference of the subject's
Laplace terms on its record up to the row and up to the row before, each with
the effects' mode solved on that truncated record, so its score is the
difference of those terms' derivatives. Each term's derivative is its unit's
gradient as the evaluation forms it (`_laplace_unit_gradient_route!`), along
the directions. Single-level fits only, as for the sequential predictions.
"""
function ctsem_row_scores(laplace::CTSEMLaplaceObjective, values::AbstractVector,
    directions::AbstractMatrix)
    size(directions, 1) == length(values) || throw(DimensionMismatch(
        "directions needs one row per raw parameter"))
    # A stated population value (`CTSEMAnchors`, `_laplace_extend`) is not a
    # direction anyone asks along: a zero row for each, after the caller's.
    theta = _laplace_extend(laplace, collect(Float64, values))
    na = length(theta) - length(values)
    along = na == 0 ? Matrix{Float64}(directions) :
        vcat(Matrix{Float64}(directions), zeros(Float64, na, size(directions, 2)))
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
    dLlevels = _laplace_level_chol_derivatives(theta, laplace.spec)
    dLsolved = _laplace_level_chol_solves(laplace.spec, Ls, dLlevels)
    _laplace_ensure_pool!(laplace)
    # The evaluation's per-unit records (`_laplace_store_unit_primal!`).
    nunits = length(units.members)
    nslot = _laplace_pool_width()
    records = (Vector{_LaplaceFactorization{Float64}}(undef, nunits),
        Vector{CTSEMBlockMatrix{Float64}}(undef, nunits), fill(false, nunits),
        zeros(Float64, nunits), fill(true, nslot), fill(NaN, nslot))
    _laplace_parallel(laplace, eachindex(units.members)) do U
        _laplace_row_scores_unit!(out, laplace, U, theta, Ls, dLlevels, dLsolved,
            records, along, first_row[units.members[U][1]])
        return true
    end
    return out
end

function _laplace_row_scores_unit!(out, laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, dLlevels, dLsolved, records,
    directions::Matrix{Float64}, offset::Integer)
    units = laplace.units
    i = units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    saved = _laplace_unit_saved(laplace, U, so)
    full = saved.full
    T = size(full.data, 2)
    K = size(directions, 2)
    npar = length(theta)
    curvatures, matrices, gated, terms, chunk_ok, chunk_bad = records
    aws = _laplace_workspace!(laplace, Float64, npar)
    row = zeros(Float64, npar)
    previous = zeros(Float64, K)
    current = zeros(Float64, K)
    try
        for t in 1:T
            _laplace_truncate!(so, full, t)
            # From the previous row's mode, as the sequential pass solves.
            ok = _laplace_store_unit_primal!(laplace, U, theta, Ls, aws, curvatures,
                matrices, gated, terms, chunk_ok, chunk_bad; retain=t > 1) &&
                first(_laplace_unit_gradient_route!(row, laplace, U, theta, Ls,
                    dLlevels, dLsolved, gated[U], matrices[U], curvatures[U]))
            if ok
                @inbounds for col in 1:K
                    s = 0.0
                    for j in 1:npar
                        s += directions[j, col] * row[j]
                    end
                    current[col] = s
                    out[offset + t, col] = s - previous[col]
                end
                copyto!(previous, current)
            else
                # The next row's difference has no defined start either.
                fill!(previous, NaN)
            end
        end
    finally
        _laplace_unit_restore!(laplace, U, so, saved)
    end
    return nothing
end

"""
    ctsem_laplace_posterior_nodes(laplace, values; nodes=5)

Each unit's posterior over its random effects, in the unit's standardised
coordinates (prior standard normal), as the adaptive Gauss-Hermite rule
`ctsem_laplace_quadrature` integrates with: `points` (`nunits` by `nodes^d` by
`d`), normalised `weights` (`nunits` by `nodes^d`), and `logmarginal`, each
unit's quadrature term. `sum_j weights[i, j] * g(points[i, j, :])` is then
E[g(u) | y_i] for a smooth `g`. The nodes are the quadrature's own, through
`_quadrature_block`'s sink. Single-level fits whose units all have the same
number of effects; a unit the rule cannot be built for is `NaN`.

Where every unit's term is a concave quadratic in its effects
(`laplace.quadratic`), the posterior is the normal at the mode with the
curvature there, and the rule's normalised weights are the Gauss-Hermite
weights themselves: the nodes are placed without evaluating any, and the
marginal is the Laplace term, which is exact there.
"""
function ctsem_laplace_posterior_nodes(laplace::CTSEMLaplaceObjective,
    values::AbstractVector; nodes::Integer=5)
    nodes >= 1 || throw(ArgumentError("need at least one quadrature node"))
    theta = _laplace_extend(laplace, collect(Float64, values))
    _laplace_check_indices(laplace, length(theta))
    units = laplace.units
    all(length(mm) == 1 for mm in units.members) || throw(ArgumentError(
        "posterior nodes need every subject to be its own unit: a single level " *
        "of random effects"))
    nunits = length(units.members)
    d = nunits == 0 ? 0 : units.dims[1]
    all(==(d), units.dims) || throw(ArgumentError(
        "posterior nodes need every unit to have the same number of effects"))
    Q = Int(nodes)^d
    points = fill(NaN, nunits, Q, d)
    weights = fill(NaN, nunits, Q)
    logmarginal = fill(NaN, nunits)
    (d == 0 || nunits == 0) && return (points=points, weights=weights,
        logmarginal=logmarginal)
    Ls = _laplace_popchols(theta, laplace.spec)
    _laplace_ensure_pool!(laplace)
    _quadrature_warm_caches(laplace, nodes)
    _laplace_parallel(laplace, 1:nunits) do U
        try
            _laplace_posterior_nodes_unit!(points, weights, logmarginal, laplace, U,
                theta, Ls, Int(nodes))
        catch err
            # As in the quadrature: a point the model cannot evaluate leaves the
            # unit NaN; a bug propagates.
            _ctsem_must_propagate(err) && rethrow()
        end
        return true
    end
    return (points=points, weights=weights, logmarginal=logmarginal)
end

function _laplace_posterior_nodes_unit!(points, weights, logmarginal,
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, nodes::Int)
    blocks = laplace.units.blocks[U]
    length(blocks) == 1 || return nothing
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    _laplace_solve_unit_mode!(laplace, U, theta, Ls)
    mode = copy(laplace.modes[U])
    M = _laplace_unit_curvature(laplace, U, theta, Ls, mode)
    factorization = _laplace_factor_repaired!(M, blocks)
    factorization.ok || return nothing
    laplace.quadratic && return _laplace_posterior_nodes_normal!(points, weights,
        logmarginal, laplace, U, theta, Ls, mode, factorization, nodes, aws)
    tree = _quadrature_children(blocks)
    context = (children=tree.children, nodes=nodes, mode=mode,
        factors=factorization.factors)
    terms = fill(-Inf, size(weights, 2))
    sink = function (j, z, term)
        terms[j] = term
        @inbounds for k in eachindex(z)
            points[U, j, k] = z[k]
        end
        return nothing
    end
    total = _quadrature_block(laplace, U, theta, Ls, tree.roots[1], copy(mode),
        context, aws; sink=sink)
    isfinite(total) || return nothing
    peak = maximum(terms)
    mass = 0.0
    @inbounds for j in eachindex(terms)
        mass += exp(terms[j] - peak)
    end
    @inbounds for j in eachindex(terms)
        weights[U, j] = exp(terms[j] - peak) / mass
    end
    logmarginal[U] = total
    return nothing
end

"""
The posterior nodes of a unit whose term is a concave quadratic: the rule
`_quadrature_block` builds at the mode (the curvature's clipping does nothing,
since such a curvature is the prior's identity plus a positive semidefinite
part), and each node's term there is the mode's less `|x|^2`, which the grid's
log weight adds back, so the normalised weights are the Gauss-Hermite ones.
"""
function _laplace_posterior_nodes_normal!(points, weights, logmarginal,
    laplace::CTSEMLaplaceObjective, U::Integer, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, mode::Vector{Float64}, factorization, nodes::Int, aws)
    k = length(mode)
    rule = _quadrature_clipped_scale(factorization.factors[1], k)
    grid, logweights = _gh_grid(k, nodes)
    mass = 0.0
    @inbounds for j in eachindex(grid)
        x = grid[j]
        z = mode .+ sqrt(2) .* (rule.scale * x)
        for c in 1:k
            points[U, j, c] = z[c]
        end
        w = exp(logweights[j] - dot(x, x))
        weights[U, j] = w
        mass += w
    end
    @inbounds for j in eachindex(grid)
        weights[U, j] /= mass
    end
    logmarginal[U] = _laplace_unit_term(laplace, U, theta, Ls, mode, aws)
    return nothing
end

"""
    ctsem_laplace_unit_dims(laplace)

Each unit's number of random-effect coordinates.
"""
ctsem_laplace_unit_dims(laplace::CTSEMLaplaceObjective) = copy(laplace.units.dims)

export ctsem_row_scores, ctsem_laplace_posterior_nodes, ctsem_laplace_unit_dims
