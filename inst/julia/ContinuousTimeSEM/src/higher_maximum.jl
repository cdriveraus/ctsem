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
# solve searched at `x` (`_ctsem_anchor_start!`), from the origin, the probes
# and the objective's known modes. NaN where the route refuses the point.
#
# With `modes`, also the maxima that search ended on (`_ctsem_anchor_modes`,
# `nothing` on a route without them): what a fit started at `x` takes as its
# opening search's further starts (`resume_units`), so that it reads `x` at
# this value. A run's found modes stay on the objective it ran on
# (`_ctsem_merge_modes!` puts the maxima it followed that a search missed in
# the known modes), and a fit from its end is another objective, built from
# the model: on the affect pilot's tied fit a fit from a climb's end read the
# point 6.7 nats below the climb's value, and a second climb's end read 14
# lower searched on a new objective than on its own (review/OPEN.md,
# 2026-11-21).
function _ctsem_searched(objective::CTSEMOptimisable, x::Vector{Float64},
    gradient_method::Symbol; modes::Bool=false)
    _ctsem_optimise_setup!(objective)
    held = nothing
    trial = try
        _ctsem_anchor_start!(objective, x)
        t = _ctsem_optimise_trial(objective, x, false, gradient_method,
            sqrt(floatmax(Float64)), _ctsem_optimise_log(objective, false))
        modes && (held = _ctsem_anchor_modes(objective; known=false))
        t
    finally
        ctsem_laplace_clear_anchors!(objective)
    end
    return (value=trial.valid ? Float64(trial.evaluated.value) : NaN, modes=held)
end
_ctsem_searched_value(objective::CTSEMOptimisable, x::Vector{Float64},
    gradient_method::Symbol) = _ctsem_searched(objective, x, gradient_method).value

# One climb, in a function of its own: what its `stop` closure captures is
# assigned once here.
#
# `inverse` is the reference's covariance, L-BFGS's first metric
# (`lbfgs_inverse`). Without one the climb learns each coordinate's scale from
# its own steps, and its gap rule reads that model: on the affect pilot's tied
# fit (819 parameters, 2026-11-21) such a climb took 427 iterations to pass the
# estimate's value from 124 nats below it, and its gap rule read under 0.005
# nats left at its 174th, 24 below, where only a mode search that moved the
# objective sent it on.
#
# With `finish` it does not stop at the bar: it ends by its own gap rule, and
# `passed`, `passed_secs` say when it went past.
function _ctsem_climb(objective::CTSEMOptimisable, start::Vector{Float64}, bar::Float64,
    maxiter::Int, gap_tol::Float64, verbose::Bool, progress_overwrite::Bool,
    progress_sink, gradient_method::Symbol, inverse::Union{Nothing,Matrix{Float64}},
    finish::Bool)
    seen = Ref(-Inf)
    passed = Ref(0)
    passed_secs = Ref(NaN)
    t0 = time()
    t = @elapsed res = ctsem_optimize(objective, start; maxiter=maxiter,
        gap_tol=gap_tol, verbose=verbose, progress_overwrite=progress_overwrite,
        progress_sink=progress_sink, progress_label="climb",
        gradient_method=gradient_method, overshoot_probe=:off,
        lbfgs_inverse=inverse,
        stop=st -> begin
            seen[] = max(seen[], -Float64(st.value))
            if seen[] > bar && passed[] == 0
                passed[] = max(1, Int(st.iteration))
                passed_secs[] = time() - t0
            end
            !finish && seen[] > bar
        end)
    value = max(seen[], Float64(res.maximum_loglik))
    return (point=collect(Float64, res.minimizer), value=value,
        iterations=Int(res.iterations), secs=t, passed=passed[],
        passed_secs=passed_secs[])
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
refuses is left out. `inverse` is the reference's covariance, each climb's
first metric (`lbfgs_inverse`); one that is not positive definite is not used.

With `finish` every start is climbed to its own gap rule's end, past the bar
or not, and the highest end is the point: what a caller that goes on to fit
from the point wants, since the climb's steps are then the fit's first, and
since the first start to pass need not lead to the highest maximum.

An end's value is the objective searched there, as the reference's and the
starts' are, on this objective, which holds what the climbs found; `higher` is
judged on it, and the maxima that search ended on come back with the point
(`mode_units`, `mode_values`) for a fit from it to start its own search from
(`_ctsem_searched`). A climb stopped at the bar ended without a search, so its
own value can rest on maxima a search does not find: where its end is not past
the bar searched, it goes on to its own end before it is judged.

Returns `higher`, the reference's value and the bar, the highest end's
`value`, `point` and modes (`mode_units` `[0]` where there are none), and per
start, in the order of `starts`, its value, its end's value (`end_values`,
searched) and its climb's own there (`climb_values`), end point (`ends`),
iterations and seconds (NaN, 0 where it was not climbed), and the iteration
and seconds at which the climb's own value first passed the bar (0, NaN where
it did not).
"""
function ctsem_climb_above(objective::CTSEMOptimisable, starts::AbstractMatrix,
    reference::AbstractVector; value::Real=-Inf, margin::Real=1.92,
    maxiter::Integer=1000, gap_tol::Real=_CTSEM_RESOLUTION_CAP, verbose::Bool=false,
    progress_overwrite::Bool=true, progress_sink=nothing, gradient_method=:adjoint,
    inverse=nothing, finish::Bool=false)
    x0 = collect(Float64, reference)
    S = Matrix{Float64}(starts)
    size(S, 1) == length(x0) || throw(ArgumentError(string("starts has ", size(S, 1),
        " rows and the reference ", length(x0), " parameters")))
    method = Symbol(gradient_method)
    n = size(S, 2)
    metric = _ctsem_lbfgs_inverse(inverse, length(x0))
    passed = zeros(Int, n)
    passed_secs = fill(NaN, n)
    ends = fill(NaN, length(x0), n)
    searched = _ctsem_searched_value(objective, x0, method)
    ref = isfinite(searched) ? max(Float64(value), searched) : Float64(value)
    isfinite(ref) || throw(ArgumentError(
        "the objective could not be evaluated at the reference, and no value was given"))
    bar = ref + Float64(margin)
    startvalues = [_ctsem_searched_value(objective, S[:, j], method) for j in 1:n]
    endvalues = fill(NaN, n)
    climbvalues = fill(NaN, n)
    modes = nothing
    iterations = zeros(Int, n)
    secs = zeros(Float64, n)
    best = -Inf
    point = copy(x0)
    higher = false
    order = sortperm([isfinite(v) ? v : -Inf for v in startvalues]; rev=true)
    for j in order
        isfinite(startvalues[j]) || continue
        # A start already past the bar needs no climb.
        if startvalues[j] > bar && !finish
            endvalues[j] = startvalues[j]
            climbvalues[j] = startvalues[j]
            ends[:, j] .= view(S, :, j)
            best = startvalues[j]
            point = S[:, j]
            modes = _ctsem_searched(objective, point, method; modes=true).modes
            higher = true
            break
        end
        verbose && (println(_console(), @sprintf(
            "climb: from start %d of %d, %.2f below the reference; the bar is %.2f above it",
            j, n, ref - startvalues[j], Float64(margin))); flush(_console()))
        climb = _ctsem_climb(objective, S[:, j], bar, Int(maxiter), Float64(gap_tol),
            verbose, progress_overwrite, progress_sink, method, metric, finish)
        endpoint = climb.point
        followed = climb.value
        iterations[j] = climb.iterations
        secs[j] = climb.secs
        passed[j] = climb.passed
        passed_secs[j] = climb.passed_secs
        read = _ctsem_searched(objective, endpoint, method; modes=true)
        if !finish && followed > bar && !(read.value > bar)
            more = _ctsem_climb(objective, endpoint, bar, Int(maxiter), Float64(gap_tol),
                verbose, progress_overwrite, progress_sink, method, metric, true)
            endpoint = more.point
            followed = more.value
            iterations[j] += more.iterations
            secs[j] += more.secs
            read = _ctsem_searched(objective, endpoint, method; modes=true)
        end
        endvalues[j] = read.value
        climbvalues[j] = followed
        ends[:, j] .= endpoint
        reached = isfinite(read.value) ? read.value : -Inf
        if reached > best
            best = reached
            point = endpoint
            modes = read.modes
        end
        verbose && (println(_console(), @sprintf(
            "climb: start %d ended %+.2f against the reference (%+.2f on the climb's own modes) after %d iterations, %.0f s%s",
            j, read.value - ref, followed - ref, iterations[j], secs[j],
            climb.passed > 0 ? @sprintf("; its own value past the bar at iteration %d, %.0f s",
                climb.passed, climb.passed_secs) : "")); flush(_console()))
        if reached > bar
            higher = true
            finish || break
        end
    end
    # None as `[0]` and a one-entry matrix, as `_ctsem_found_columns`.
    mode_units, mode_values = modes === nothing ? ([0], zeros(1, 1)) :
        (collect(Int, modes.units), Matrix{Float64}(modes.values))
    return (higher=higher, reference=ref, searched=searched, bar=bar, value=best,
        point=point, start_values=startvalues, end_values=endvalues,
        climb_values=climbvalues, iterations=iterations, secs=secs, passed=passed,
        passed_secs=passed_secs, ends=ends, metric=metric !== nothing,
        mode_units=mode_units, mode_values=mode_values)
end

export ctsem_climb_above
