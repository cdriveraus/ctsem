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
    ctsem_laplace_unit_terms(laplace, draws; subjects=nothing)

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

`subjects`, when given, scores only the units holding those subjects, by the
same mode solve and term the evaluation uses; `value` is then `NaN` and every
other unit's term `NaN`, and `converged` covers the scored units. A
leave-one-subject-out fold reads its held-out subject's unit and nothing else.
"""
function ctsem_laplace_unit_terms(laplace::CTSEMLaplaceObjective, draws::AbstractMatrix;
    subjects=nothing)
    ndraws = size(draws, 2)
    nunits = length(laplace.units.members)
    value = fill(NaN, ndraws)
    terms = fill(NaN, nunits, ndraws)
    converged = fill(false, ndraws)
    unit = _laplace_unit_of_subject(laplace)
    # A bare integer when R sends one subject, so iterate rather than index.
    chosen = subjects === nothing ? Int[] : unique(Int[unit[Int(i)] for i in subjects])
    for s in 1:ndraws
        theta = collect(Float64, view(draws, :, s))
        ok = try
            if subjects === nothing
                result = ctsem_laplace_evaluate(laplace, theta; gradient=false)
                value[s] = result.value
                terms[:, s] .= result.unit_loglik
                converged[s] = result.converged
            else
                converged[s] = _laplace_unit_terms_of!(view(terms, :, s), laplace,
                    theta, chosen)
            end
            true
        catch err
            err isa Union{DomainError,LinearAlgebra.PosDefException,
                LinearAlgebra.SingularException,LinearAlgebra.LAPACKException} ||
                rethrow()
            false
        end
        if !ok
            terms[:, s] .= NaN
            converged[s] = false
        end
    end
    return (value=value, unit_loglik=terms, converged=converged, unit=unit)
end

"""
Units `chosen`'s terms at `theta` into `terms`, by the evaluation's own mode
solve and term (`_laplace_unit_term`), and whether every one's mode was found.
"""
function _laplace_unit_terms_of!(terms::AbstractVector{Float64},
    laplace::CTSEMLaplaceObjective, theta::Vector{Float64}, chosen::Vector{Int})
    _laplace_check_indices(laplace, length(theta))
    _laplace_ensure_pool!(laplace)
    Ls = _laplace_popchols(theta, laplace.spec)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    ok = true
    for U in chosen
        _laplace_solve_unit_mode!(laplace, U, theta, Ls, aws)
        terms[U] = _laplace_unit_term(laplace, U, theta, Ls, laplace.modes[U], aws)
        ok &= laplace.inner_converged[U]
    end
    return ok
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

Deviations each member's filter integrates (`CTSEMHybrid`) are drawn with the
effects, from the same approximation over both (`_laplace_effect_posterior`):
the coordinates the layer holding them would draw, in its order, so that the
draws are the layer's.

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
    hyb = laplace.objective.params.hybrid
    nh = _hybrid_on(hyb) ? length(hyb.manifests) : 0
    # A unit's coordinates: its effects, then each member's deviations, one per
    # hybrid item whether the member observes it or not (`_laplace_effect_posterior`).
    widths = [units.dims[U] + nh * length(units.members[U]) for U in 1:nunits]
    z = randn(MersenneTwister(seed), sum(widths; init=0), ndraws)
    effects = zeros(Float64, sum(units.dims; init=0), ndraws)
    # Each subject's drawn deviations on the raw scale, `psi * d`, per item.
    nsubjects = length(laplace.objective.subject_objectives)
    shifts = zeros(Float64, nsubjects, nh, ndraws)
    lv = laplace.spec.levels[1]
    psi = [_hybrid_item_sd(hyb, lv, theta, k) for k in 1:nh]
    logq = zeros(Float64, nunits, ndraws)
    logprior = zeros(Float64, nunits, ndraws)
    converged = fill(true, nunits)
    zbase = 0
    ubase = 0
    for U in 1:nunits
        d = units.dims[U]
        D = widths[U]
        zslice = (zbase + 1):(zbase + D)
        uslice = (ubase + 1):(ubase + d)
        zbase += D
        ubase += d
        D == 0 && continue
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        uhat = copy(laplace.modes[U])
        post = _laplace_effect_posterior(laplace, U, theta, Ls, uhat, nh)
        F = post.ok ? cholesky(Symmetric((post.cov + transpose(post.cov)) ./ 2);
            check=false) : nothing
        if !laplace.inner_converged[U] || F === nothing || !issuccess(F)
            converged[U] = false
            for s in 1:ndraws
                effects[uslice, s] .= uhat
            end
            continue
        end
        halflogdet = sum(log, diag(F.L)) + D * log(scale)
        constant = -D / 2 * log(2pi)
        for s in 1:ndraws
            zs = view(z, zslice, s)
            w = post.mean .+ scale .* (F.L * zs)
            effects[uslice, s] .= view(w, 1:d)
            for (m, i) in enumerate(units.members[U]), (c, k) in enumerate(post.items)
                shifts[i, k, s] = psi[k] * w[d + (m - 1) * nh + c]
            end
            logq[U, s] = constant - halflogdet - dot(zs, zs) / 2
            logprior[U, s] = constant - dot(w, w) / 2
        end
    end
    # Sized from the first trace rather than from the subjects' data, so the
    # rows are exactly the ones the filter reports.
    # The draws are coordinates u; the filter takes natural deviations S u.
    llrow = Matrix{Float64}(undef, 0, 0)
    starts, ndev = _laplace_deviation_offsets(laplace)
    deviations = zeros(Float64, ndev)
    for s in 1:ndraws
        _laplace_deviations!(deviations, laplace, Ls, view(effects, :, s), starts)
        persubject = _laplace_filter_values(laplace, theta; effects=deviations)
        for i in 1:nsubjects, k in 1:nh
            persubject[i, lv.re_index[hyb.effect[k]]] += shifts[i, k, s]
        end
        trace = ctsem_kalman(laplace.objective, persubject; subject_matrices=false,
            fields=["llrow"])
        s == 1 && (llrow = fill(NaN, length(trace.llrow), ndraws))
        llrow[:, s] .= trace.llrow
    end
    return (llrow=llrow, logq=logq, logprior=logprior, converged=converged,
        unit=_laplace_unit_of_subject(laplace))
end

"""
    _laplace_effect_posterior(laplace, U, theta, Ls, uhat, nh)

The Laplace approximation to unit `U`'s effects' posterior at its mode `uhat`:
`(ok, mean, cov, items)`, `cov` the inverse of the curvature as the fit
factorizes it (`_laplace_factor_repaired!`). With `nh` deviations each member's
filter integrates, the coordinates are the effects and then, member by member,
one per hybrid item, `items` their order -- the items by the intercept each
shifts, as the layer holding them as coordinates orders its own variances. An
item the member observes takes its conditional mode and the joint curvature
with the effects (`_hybrid_joint_curvature`); one it never observes keeps its
prior, `N(0, 1)`, independent of the rest.
"""
function _laplace_effect_posterior(laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, uhat::Vector{Float64},
    nh::Integer)
    d = length(uhat)
    blocks = laplace.units.blocks[U]
    members = laplace.units.members[U]
    failed = (ok=false, mean=Float64[], cov=zeros(0, 0), items=Int[])
    if nh == 0
        M = _laplace_unit_curvature(laplace, U, theta, Ls, uhat)
        fac = _laplace_factor_repaired!(M, blocks)
        fac.ok || return failed
        cov = zeros(Float64, d, d)
        for j in 1:d
            e = zeros(Float64, d); e[j] = 1.0
            cov[:, j] = _laplace_block_solve(fac.factors, fac.coupling, blocks, e)
        end
        return (ok=true, mean=copy(uhat), cov=cov, items=Int[])
    end
    hyb = laplace.objective.params.hybrid
    items = sortperm(hyb.effect)
    slot = invperm(items)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    modes = _hybrid_member_modes(laplace, U, theta, Ls, uhat, aws)
    modes === nothing && return failed
    layout = _hybrid_layout(d, modes.q)
    D = d + nh * length(members)
    # Where each coordinate of the joint layout sits among the unit's.
    at = collect(1:layout.n)
    zj = zeros(Float64, layout.n)
    zj[1:d] .= uhat
    mean = zeros(Float64, D)
    mean[1:d] .= uhat
    for m in eachindex(members), (c, k) in enumerate(modes.items[m])
        a = layout.offs[m] + c
        at[a] = d + (m - 1) * nh + slot[k]
        zj[a] = modes.dhat[m][c]
        mean[at[a]] = modes.dhat[m][c]
    end
    cov = Matrix{Float64}(LinearAlgebra.I, D, D)
    if layout.n > 0
        Mj = d > 0 ? _hybrid_joint_curvature(laplace, U, theta, Ls, zj, layout) :
            CTSEMBlockMatrix{Float64}([copy(modes.A[1])], [Matrix{Float64}[]])
        jblocks = [CTSEMLaplaceBlock(0, layout.n, collect(eachindex(members)), 1, Int[])]
        fac = _laplace_factor_repaired!(Mj, jblocks)
        fac.ok || return failed
        for j in 1:layout.n
            e = zeros(Float64, layout.n); e[j] = 1.0
            col = _laplace_block_solve(fac.factors, fac.coupling, jblocks, e)
            for i in 1:layout.n
                cov[at[i], at[j]] = col[i]
            end
        end
    end
    return (ok=true, mean=mean, cov=cov, items=items)
end

export ctsem_laplace_unit_terms, ctsem_laplace_effect_draws
