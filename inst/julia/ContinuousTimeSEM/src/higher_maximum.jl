################################################################################
# A climb from a candidate: is there a point above the estimate?
################################################################################
#
# An estimate's own checks are local. Its certification, its flat-direction
# walk and its interval check read the objective's curvature at the estimate
# and its value a few standard errors out, and a higher maximum behind a dip
# shows in none of them: on the affect pilot's tied fit (819 parameters) the
# estimate was certified with every direction trusted, its Hessian put a
# maximum 41.6 nats higher 1,899 nats below it, and on the straight line
# between the two the objective fell 69 nats before it rose
# (review/OPEN.md, 2026-11-21). What showed it was the joint sampler, whose
# chains left the estimate's normal within their warm-up, and a climb of the
# fit's objective from a chain's mean, which passed the estimate's value in a
# third of an hour where the fit had taken seventeen.
#
# So the test of a candidate is a climb from it. Its value says nothing: the
# chains' means were 77 and 114 nats below the estimate, and every point of the
# line short of four fifths of the way was below it too. The climb is
# `ctsem_optimize`'s L-BFGS stage alone -- no finish, no Hessian, no
# certification -- stopped by the optimiser's own `stop` once it is a margin
# above the reference, which is all the question needs: what is done with a
# higher point (a fit from it) is the caller's.

# The objective at `x` as a run that starts there reads it: each unit's inner
# solve searched at `x` (`_ctsem_anchor_start!`), which is how the climbs'
# values are read too. NaN where the route refuses the point.
function _ctsem_searched_value(objective::CTSEMOptimisable, x::Vector{Float64},
    gradient_method::Symbol)
    _ctsem_optimise_setup!(objective)
    trial = try
        _ctsem_anchor_start!(objective, x)
        _ctsem_optimise_trial(objective, x, false, gradient_method,
            sqrt(floatmax(Float64)), _ctsem_optimise_log(objective, false))
    finally
        ctsem_laplace_clear_anchors!(objective)
    end
    return trial.valid ? Float64(trial.evaluated.value) : NaN
end

# One climb, in a function of its own: its `stop` closure captures `seen`
# and `bar`, each assigned once here.
function _ctsem_climb(objective::CTSEMOptimisable, start::Vector{Float64}, bar::Float64,
    maxiter::Int, gap_tol::Float64, verbose::Bool, progress_overwrite::Bool,
    progress_sink, gradient_method::Symbol)
    seen = Ref(-Inf)
    t = @elapsed res = ctsem_optimize(objective, start; maxiter=maxiter,
        gap_tol=gap_tol, verbose=verbose, progress_overwrite=progress_overwrite,
        progress_sink=progress_sink, progress_label="climb",
        gradient_method=gradient_method, overshoot_probe=:off,
        stop=st -> (seen[] = max(seen[], -Float64(st.value)); seen[] > bar))
    value = max(seen[], Float64(res.maximum_loglik))
    return (point=collect(Float64, res.minimizer), value=value,
        iterations=Int(res.iterations), secs=t)
end

"""
    ctsem_climb_above(objective, starts, reference; value, margin, maxiter, ...)

Whether `objective` has a point more than `margin` above its value at
`reference`, looked for by climbing from the columns of `starts`.

The reference's value is the larger of `value`, the caller's (a fit's own, read
with the modes its run ended with), and the objective searched at `reference`.
The starts are climbed in the order of their own values, highest first, each
by L-BFGS with no finish and no Hessian, stopped once it is past the bar, by
its own gap rule (`gap_tol`, a tenth of a standard error: the end is no
estimate) or at `maxiter`; the first to pass ends the search. A start the route
refuses is left out.

Returns `higher`, the reference's value and the bar, the highest end's
`value` and `point`, and per start, in the order of `starts`, its value, its
climb's end value, iterations and seconds (NaN, 0 where it was not climbed).
"""
function ctsem_climb_above(objective::CTSEMOptimisable, starts::AbstractMatrix,
    reference::AbstractVector; value::Real=-Inf, margin::Real=1.92,
    maxiter::Integer=1000, gap_tol::Real=_CTSEM_RESOLUTION_CAP, verbose::Bool=false,
    progress_overwrite::Bool=true, progress_sink=nothing, gradient_method=:adjoint)
    x0 = collect(Float64, reference)
    S = Matrix{Float64}(starts)
    size(S, 1) == length(x0) || throw(ArgumentError(string("starts has ", size(S, 1),
        " rows and the reference ", length(x0), " parameters")))
    method = Symbol(gradient_method)
    n = size(S, 2)
    searched = _ctsem_searched_value(objective, x0, method)
    ref = isfinite(searched) ? max(Float64(value), searched) : Float64(value)
    isfinite(ref) || throw(ArgumentError(
        "the objective could not be evaluated at the reference, and no value was given"))
    bar = ref + Float64(margin)
    startvalues = [_ctsem_searched_value(objective, S[:, j], method) for j in 1:n]
    endvalues = fill(NaN, n)
    iterations = zeros(Int, n)
    secs = zeros(Float64, n)
    best = -Inf
    point = copy(x0)
    higher = false
    order = sortperm([isfinite(v) ? v : -Inf for v in startvalues]; rev=true)
    for j in order
        isfinite(startvalues[j]) || continue
        verbose && (println(_console(), @sprintf(
            "climb: from start %d of %d, %.2f below the reference; the bar is %.2f above it",
            j, n, ref - startvalues[j], Float64(margin))); flush(_console()))
        climb = _ctsem_climb(objective, S[:, j], bar, Int(maxiter), Float64(gap_tol),
            verbose, progress_overwrite, progress_sink, method)
        endvalues[j] = climb.value
        iterations[j] = climb.iterations
        secs[j] = climb.secs
        if climb.value > best
            best = climb.value
            point = climb.point
        end
        verbose && (println(_console(), @sprintf(
            "climb: start %d ended %+.2f against the reference after %d iterations, %.0f s",
            j, climb.value - ref, climb.iterations, climb.secs)); flush(_console()))
        if climb.value > bar
            higher = true
            break
        end
    end
    return (higher=higher, reference=ref, searched=searched, bar=bar, value=best,
        point=point, start_values=startvalues, end_values=endvalues,
        iterations=iterations, secs=secs)
end

export ctsem_climb_above
