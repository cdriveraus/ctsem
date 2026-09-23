# Refit-free leave-one-out for Laplace fits ----------------------------------
#
# Two batched evaluations that the R side's `ctLOO(method = 'psis')` smooths
# with Pareto-smoothed importance sampling. Both exist to put a whole set of
# draws through one bridge call: a round trip costs ~41 ms, and a per-draw loop
# from R would spend most of its time there.
#
#   ctsem_laplace_unit_terms     the Laplace objective at each column of a
#                                matrix of population parameter draws, with
#                                each unit's own term -- leave one *unit* out.
#   ctsem_laplace_effect_draws   draws of each unit's random effects from the
#                                Gaussian the Laplace approximation fits at the
#                                mode, and every row's one-step-ahead log
#                                likelihood at each -- leave one *row* out,
#                                conditional on the population parameters.

"""
    ctsem_laplace_unit_terms(laplace, draws)

The Laplace objective without gradient at every column of `draws`
(`npar x ndraws`), returning

- `value`: the log posterior (every unit's term plus the log prior),
- `unit_loglik`: `nunits x ndraws`, each unit's approximated log marginal,
- `converged`: whether every inner mode was found at that draw,
- `unit`: the unit each subject belongs to.

A draw at which the model cannot be evaluated -- a factorization that fails, a
non-finite term -- comes back with a `NaN` value rather than an error, so the
caller can drop it and say how many it dropped. Only the linear algebra and
domain failures a proposal draw in the tails can provoke are treated that way;
anything else is a bug and is rethrown.
"""
function ctsem_laplace_unit_terms(laplace::CTSEMLaplaceObjective, draws::AbstractMatrix)
    ndraws = size(draws, 2)
    nunits = length(laplace.units.members)
    value = fill(NaN, ndraws)
    terms = fill(NaN, nunits, ndraws)
    converged = fill(false, ndraws)
    for s in 1:ndraws
        result = try
            ctsem_laplace_evaluate(laplace, collect(Float64, view(draws, :, s));
                gradient=false)
        catch err
            err isa Union{DomainError,LinearAlgebra.PosDefException,
                LinearAlgebra.SingularException,LinearAlgebra.LAPACKException} ||
                rethrow()
            nothing
        end
        result === nothing && continue
        value[s] = result.value
        terms[:, s] .= result.unit_loglik
        converged[s] = result.converged
    end
    return (value=value, unit_loglik=terms, converged=converged,
        unit=_laplace_unit_of_subject(laplace))
end

"""The unit each subject belongs to, in subject order."""
function _laplace_unit_of_subject(laplace::CTSEMLaplaceObjective)
    nsubjects = length(laplace.objective.subject_objectives)
    out = zeros(Int, nsubjects)
    for (U, members) in enumerate(laplace.units.members)
        for i in members
            out[i] = U
        end
    end
    return out
end

"""
    ctsem_laplace_effect_draws(laplace, values, ndraws; seed, scale)

Draws of every unit's random effects from `N(uhat_U, scale^2 M_U^{-1})` at the
population parameters `values`, and each data row's log likelihood at each.

`uhat_U` is the unit's inner mode and `M_U` its curvature there -- the same
Gaussian the Laplace approximation integrates, so for a linear Gaussian model
it is the exact conditional posterior of the effects. The covariance is formed
densely from the block factorization, one solve per column; a unit's effect
dimension is small next to its data, so this is not where the cost is.

Returns

- `llrow`: `nrows x ndraws`, each row's one-step-ahead log likelihood, from the
  filter, with every subject at the parameters its drawn effects imply,
- `logq`: `nunits x ndraws`, the log density of each unit's draw under the
  proposal,
- `logprior`: `nunits x ndraws`, the log density of each unit's draw under its
  standard normal prior -- the other half of `g_U` besides the rows,
- `converged`: per unit, whether its mode was found and its curvature
  factorized; the draws for a unit that fails are the mode repeated and should
  not be used,
- `unit`: the unit each subject belongs to.

`seed` makes the draws reproducible from the R side's own seed. `scale` widens
the proposal; see the R side (`.ctBackendLOOPsis`) for why it is above one.
"""
function ctsem_laplace_effect_draws(laplace::CTSEMLaplaceObjective,
    values::AbstractVector, ndraws::Integer; seed::Integer=1, scale::Real=1.0)
    theta = collect(Float64, values)
    _laplace_check_indices(laplace, length(theta))
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    units = laplace.units
    nunits = length(units.members)
    total = sum(units.dims; init=0)
    z = randn(MersenneTwister(seed), total, ndraws)
    effects = zeros(Float64, total, ndraws)
    logq = zeros(Float64, nunits, ndraws)
    logprior = zeros(Float64, nunits, ndraws)
    converged = fill(true, nunits)
    base = 0
    for U in 1:nunits
        d = units.dims[U]
        slice = (base + 1):(base + d)
        base += d
        d == 0 && continue
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        uhat = copy(laplace.modes[U])
        blocks = units.blocks[U]
        M = _laplace_unit_curvature(laplace, U, theta, Ls, uhat)
        fac = _laplace_factor_repaired!(M, blocks)
        C = zeros(Float64, d, d)
        if fac.ok
            for j in 1:d
                e = zeros(Float64, d); e[j] = 1.0
                C[:, j] = _laplace_block_solve(fac.factors, fac.coupling, blocks, e)
            end
        end
        F = fac.ok ? cholesky(Symmetric((C + transpose(C)) ./ 2); check=false) : nothing
        if !laplace.inner_converged[U] || F === nothing || !issuccess(F)
            converged[U] = false
            for s in 1:ndraws
                effects[slice, s] .= uhat
            end
            continue
        end
        halflogdet = sum(log, diag(F.L)) + d * log(scale)
        constant = -d / 2 * log(2pi)
        for s in 1:ndraws
            zs = view(z, slice, s)
            u = uhat .+ scale .* (F.L * zs)
            effects[slice, s] .= u
            logq[U, s] = constant - halflogdet - dot(zs, zs) / 2
            logprior[U, s] = constant - dot(u, u) / 2
        end
    end
    # Sized from the first trace rather than from the subjects' data, so the
    # rows are exactly the ones the filter reports.
    llrow = Matrix{Float64}(undef, 0, 0)
    for s in 1:ndraws
        persubject = ctsem_laplace_subject_values(laplace, theta, view(effects, :, s))
        trace = ctsem_kalman(laplace.objective, persubject; subject_matrices=false,
            fields=["llrow"])
        s == 1 && (llrow = fill(NaN, length(trace.llrow), ndraws))
        llrow[:, s] .= trace.llrow
    end
    return (llrow=llrow, logq=logq, logprior=logprior, converged=converged,
        unit=_laplace_unit_of_subject(laplace))
end

export ctsem_laplace_unit_terms, ctsem_laplace_effect_draws
