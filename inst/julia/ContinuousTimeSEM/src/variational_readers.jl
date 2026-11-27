# Readers of a variational fit's objective per unit -----------------------------
#
# What reads the fit's objective unit by unit -- leave-one-out's unit terms,
# the score tests' row scores and Hessian columns, ctLaplaceCheck's gap --
# reads each unit's bound, the term the fit
# maximised. None of it moves the approximations the fit holds: each
# evaluation starts from them and leaves them as it found them (`keep =
# false`).

"""
    ctsem_unit_terms(V::CTSEMVariational, draws; subjects = nothing)

`ctsem_unit_terms` on the variational objective: at every column of `draws`
(`npar x ndraws`) the bound summed with the prior (`value`) and each unit's
bound (`unit_loglik`), as the Laplace objective's method returns its terms.
`subjects` scores only the units holding them, `value` then `NaN`. A unit
with every row withheld scores 0, its approximation the prior. A draw where
the bound cannot be evaluated comes back `NaN`, for the caller to drop.
"""
function ctsem_unit_terms(V::CTSEMVariational, draws::AbstractMatrix; subjects=nothing)
    laplace = V.laplace
    ndraws = size(draws, 2)
    nunits = length(V.layouts)
    value = fill(NaN, ndraws)
    terms = fill(NaN, nunits, ndraws)
    converged = fill(false, ndraws)
    unit = _laplace_unit_of_subject(laplace)
    # A bare integer when R sends one subject, so iterate rather than index.
    held = subjects === nothing ? Int[] : unique(Int[unit[Int(i)] for i in subjects])
    chosen = subjects === nothing ? nothing : _laplace_coreset(laplace, held, ones(length(held)))
    for s in 1:ndraws
        r = try
            ctsem_variational_evaluate(V, collect(Float64, view(draws, :, s));
                gradient=false, coreset=chosen, keep=false)
        catch err
            err isa Union{DomainError,LinearAlgebra.PosDefException,
                LinearAlgebra.SingularException,LinearAlgebra.LAPACKException} ||
                rethrow()
            nothing
        end
        (r === nothing || !r.converged) && continue
        terms[:, s] .= r.unit_loglik
        converged[s] = true
        chosen === nothing && (value[s] = r.value)
    end
    return (value=value, unit_loglik=terms, converged=converged, unit=unit)
end

"""
    ctsem_row_scores(V::CTSEMVariational, values, directions)

Each data row's score on a variational fit, as the Laplace route's are formed
(`ctsem_row_scores(laplace, ...)`) with the unit's bound in place of its
Laplace term: a row's log predictive density is the difference of the
subject's bounds on its record up to the row and up to the row before, each
solved on that truncated record (from the previous row's approximation), so
its score is the difference of those bounds' gradients along the
directions. A subject's rows sum to its unit's score on the bound, which the
fit's estimate sets to zero over the units. Single-level fits only.
"""
function ctsem_row_scores(V::CTSEMVariational, values::AbstractVector,
    directions::AbstractMatrix)
    laplace = V.laplace
    npar = length(values)
    size(directions, 1) == npar || throw(DimensionMismatch(
        "directions needs one row per raw parameter"))
    units = laplace.units
    all(length(mm) == 1 for mm in units.members) || throw(ArgumentError(
        "row scores need every subject to be its own unit: a single level of " *
        "random effects"))
    K = size(directions, 2)
    subjects = laplace.objective.full_subjects
    rows = [size(s.data, 2) for s in subjects]
    first_row = cumsum(vcat(0, rows[1:(end - 1)]))
    out = fill(NaN, sum(rows; init=0), K)
    K == 0 && return out
    ctx = _vf_ctx(laplace, values)
    _laplace_ensure_pool!(laplace)
    D = Matrix{Float64}(directions)
    _laplace_parallel(laplace, eachindex(units.members)) do U
        _vf_row_scores_unit!(out, V, U, ctx, D, first_row[units.members[U][1]])
        return true
    end
    return out
end

function _vf_row_scores_unit!(out, V::CTSEMVariational, U::Integer, ctx::_VCtx,
    directions::Matrix{Float64}, offset::Integer)
    laplace = V.laplace
    i = laplace.units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    own = _laplace_record(so)
    full = _laplace_record(laplace.objective.full_subjects[i])
    T = size(full.data, 2)
    K = size(directions, 2)
    npar = size(directions, 1)
    previous = zeros(K)
    current = zeros(K)
    lam = zeros(V.layouts[U].p)
    try
        for t in 1:T
            _laplace_truncate!(so, full, t)
            s = try
                _vf_solve(V, Int(U), ctx, lam)
            catch err
                _ctsem_must_propagate(err) && rethrow()
                nothing
            end
            if s !== nothing && isfinite(s.e.B)
                lam = s.lam
                for col in 1:K
                    current[col] = dot(view(directions, :, col), view(s.e.GE, 1:npar))
                    out[offset + t, col] = current[col] - previous[col]
                end
                copyto!(previous, current)
            else
                # The next row's difference has no defined start either.
                fill!(previous, NaN)
            end
        end
    finally
        _laplace_set_record!(so, own)
    end
    return nothing
end

"""
    ctsem_hessian_columns(V::CTSEMVariational, values, columns; step = 1e-4)

The fit's Hessian's columns `columns` at `values`, as `ctsem_hessian(V, ...)`
forms the whole: per unit reading column `j`, `B_tt[:, j]` and `B_lt[:, j]`
by central differences in parameter `j` with the unit's approximation held
where it is solved at `values`, and the approximation's response `B_tl M^-1
B_lt[:, j]` as the difference of the parameters' gradient along `v = M^-1
B_lt[:, j]` in the approximation's parameters. Four of a unit's term
evaluations a column, whatever the number of parameters: the score tests'
columns along new terms of an extended model. Leaves the approximations the
objective holds, and their derivatives, where they were.
"""
function ctsem_hessian_columns(V::CTSEMVariational, values::AbstractVector, columns;
    step::Real=1e-4, kwargs...)
    x = collect(Float64, values)
    n = length(x)
    cols = collect(Int, columns)
    H = zeros(n, length(cols))
    isempty(cols) && return H
    h = Float64(step)
    laplace = V.laplace
    nunits = length(V.layouts)
    kept = (copy(V.lam_trial), copy(V.shares_trial), V.has_trial, copy(V.unit_value),
        copy(V.theta_trial))
    r0 = ctsem_variational_evaluate(V, x; gradient=false)
    lam = [copy(l) for l in V.lam_trial]
    V.lam_trial, V.shares_trial, V.has_trial, V.unit_value, V.theta_trial = kept
    isfinite(r0.value) || throw(ErrorException(
        "the variational bound cannot be evaluated where the Hessian's columns were asked for"))
    ctx = _vf_ctx(laplace, x)
    nx = length(ctx.theta)
    ctxy = _vf_ctx(laplace, _laplace_reader_point(x, h))
    wts = [_laplace_unit_weight(laplace, U) for U in 1:nunits]
    facs = Vector{Any}(nothing, nunits)
    G0 = zeros(nx, nunits)
    Gy = zeros(nx, nunits)
    ok = _laplace_parallel(laplace, 1:nunits) do U
        lay = V.layouts[U]
        e = _vf_unit_terms(V, U, ctx, lam[U], lay.p > 0)
        isfinite(e.B) || return false
        G0[:, U] .= e.GE
        if lay.p > 0
            fac, cpl, _ = _vf_block_factor(_vf_assemble(lay, [-M for M in e.Hm], e.Hneg),
                lay.pblocks)
            fac === nothing && return false
            facs[U] = (fac, cpl)
        end
        ey = _vf_unit_terms(V, U, ctxy, lam[U], false)
        isfinite(ey.B) ? (Gy[:, U] .= ey.GE) : fill!(view(Gy, :, U), 1.0)
        return true
    end
    ok || throw(ErrorException(
        "a unit's curvature in its approximation does not factor at these values"))
    col = zeros(nx, nunits)
    gpp = zeros(nx)
    gpm = zeros(nx)
    for (k, j) in enumerate(cols)
        xp = copy(x)
        xp[j] += h
        xm = copy(x)
        xm[j] -= h
        ctxp = _vf_ctx(laplace, xp)
        ctxm = _vf_ctx(laplace, xm)
        us = [U for U in 1:nunits if G0[j, U] != 0 || Gy[j, U] != 0]
        okj = _laplace_parallel(laplace, us) do U
            ep = _vf_unit_terms(V, U, ctxp, lam[U], false)
            em = _vf_unit_terms(V, U, ctxm, lam[U], false)
            (isfinite(ep.B) && isfinite(em.B)) || return false
            col[:, U] .= (ep.GE .- em.GE) ./ (2 * h)
            if facs[U] !== nothing
                fac, cpl = facs[U]
                lay = V.layouts[U]
                v = _laplace_block_solve(fac, cpl, lay.pblocks, (ep.g .- em.g) ./ (2 * h))
                size_v = norm(v)
                if size_v > 0
                    eps = h / size_v
                    e1 = _vf_unit_terms(V, U, ctx, lam[U] .+ eps .* v, false)
                    e2 = _vf_unit_terms(V, U, ctx, lam[U] .- eps .* v, false)
                    (isfinite(e1.B) && isfinite(e2.B)) || return false
                    col[:, U] .+= (e1.GE .- e2.GE) ./ (2 * eps)
                end
            end
            col[:, U] .*= wts[U]
            return true
        end
        okj || throw(ErrorException(
            "a unit's bound cannot be evaluated beside these values"))
        for U in us
            H[:, k] .+= view(col, 1:n, U)
        end
        fill!(gpp, 0.0)
        fill!(gpm, 0.0)
        _ctsem_log_prior_gradient!(gpp, laplace.objective, ctxp.theta)
        _ctsem_log_prior_gradient!(gpm, laplace.objective, ctxm.theta)
        H[:, k] .+= (view(gpp, 1:n) .- view(gpm, 1:n)) ./ (2 * h)
    end
    return H
end

# ---- ctLaplaceCheck ----------------------------------------------------------------
#
# No coreset of the units' bounds for the samplers' densities: with one the
# interval check took as long or longer, on gD3 81-98 s against 55-81 without
# and on gD1 15-17 against 16-17, since its draws cost what they did
# (review/VARIATIONAL-downstream-2026-11-27.md).

"""
    ctsem_laplace_correction(V::CTSEMVariational, values, hessian; nodes, step,
        laplace_gradient)

`ctLaplaceCheck`'s gap and correction on a variational fit: the units'
integral by product quadrature on the Laplace objective the bound is built on,
against the bound the fit maximised, and the step `-H^-1` times the gap's
gradient. With `laplace_gradient` (a fit not at the bound's optimum) the
bound's own gradient is added, so the step is the quadrature objective's.
"""
function ctsem_laplace_correction(V::CTSEMVariational, values::AbstractVector,
    hessian::AbstractMatrix; nodes::Integer=5, step::Real=1e-3,
    laplace_gradient::Bool=false)
    theta = collect(Float64, values)
    b = ctsem_variational_evaluate(V, theta; keep=false)
    isfinite(b.value) || throw(ErrorException(
        "the variational bound cannot be evaluated where the check was asked for"))
    gap = _quadrature_gap_gradient(V.laplace, theta, nodes, step) .+
        ctsem_laplace_evaluate(V.laplace, theta; gradient=true).gradient
    laplace_gradient || (gap .-= b.gradient)
    H = Symmetric((hessian .+ transpose(hessian)) ./ 2)
    delta, dropped = _correction_step(H, gap)
    quadrature = ctsem_laplace_quadrature(V.laplace, theta; nodes=nodes).value
    return (delta=delta, corrected=theta .+ delta, gap_gradient=gap,
        quadrature=quadrature, laplace=b.value, gap=quadrature - b.value,
        nodes=Int(nodes), dropped_directions=dropped)
end

ctsem_laplace_quadrature(V::CTSEMVariational, args...; kwargs...) =
    ctsem_laplace_quadrature(V.laplace, args...; kwargs...)
