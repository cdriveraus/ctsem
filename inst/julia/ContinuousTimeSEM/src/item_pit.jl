# Probability integral transforms of non-Gaussian observations, for misfit
# tests that need prediction errors which are standard normal under the model.
#
# A Gaussian item's one-step error, standardised, is N(0, 1) and independent
# of every earlier row under the model; a binary, ordinal, count or censored
# item's is not, so the error families had nothing valid to work on for them.
# What is: the item's predictive distribution function at what was observed,
# F(y), which is uniform under the model when it is evaluated given everything
# the filter has conditioned on before it -- the subject's earlier rows and,
# within the row, the items the filter has already applied (the categorical
# items one at a time, in item order, and then the Gaussian block). That is the
# Rosenblatt transform, so every row's values are independent of each other
# and of every earlier row. Where y has probability mass -- a category, a
# count, a censored value, a value recorded to an interval -- F jumps there,
# and a uniform drawn between the two sides (on the R side, from fixed draws,
# so the result stays smooth in the parameters) is uniform again.
#
# The bounds are taken in the update itself, from the same predictive the
# update integrates against (`ηbar`, `s2` in `_ekf_binary_update!`), through a
# record the Kalman trace hands the categorical update when they are wanted;
# nothing changes for any other pass.

"""
    CTSEMPITRecord

One row's categorical observations, collecting each one's distribution
function bounds as the update applies it, in `rows` order.
"""
mutable struct CTSEMPITRecord
    lower::Matrix{Float64}
    upper::Matrix{Float64}
    # Each observation's log predictive probability (or density), where a
    # caller needs it: empty otherwise.
    logz::Matrix{Float64}
    row::Int
    items::Vector{Int}
    k::Int
end

# The Kalman trace asks for the bounds when `pitwanted` is set: `ctsem_kalman`
# with the field `pitlower` or `pitupper`.
function _record_binary!(trace::CTSEMKalmanTrace, ws, pars, data, obs_col, rows, n)
    trace.pitwanted || return nothing
    isempty(rows) && return nothing
    return CTSEMPITRecord(trace.pitlower, trace.pitupper, zeros(0, 0),
        trace.offset + Int(obs_col), collect(Int, rows), 0)
end

function _record_binary_step!(record::CTSEMPITRecord, ws, c, ηbar, s2, y, nodes,
    weights, thresholds, kind, n)
    record.k += 1
    i = record.items[record.k]
    lo, hi = _ctsem_item_cdf(Float64(_primal(ηbar)), Float64(_primal(s2)), y,
        thresholds, kind, nodes, weights)
    @inbounds record.lower[record.row, i] = lo
    @inbounds record.upper[record.row, i] = hi
    moments = _binary_moments(ηbar, sqrt(s2), y, nodes, weights, thresholds, kind)
    isempty(record.logz) || (@inbounds record.logz[record.row, i] = Float64(_primal(first(moments))))
    return moments
end

"""
    _ctsem_item_cdf(ηbar, s2, y, thresholds, kind, nodes, weights)

`(F(y-), F(y))` for one non-Gaussian observation whose linear predictor is
`N(ηbar, s2)`: the predictive probability of a value below `y`, and of one at
most `y`; equal where `y` has a density. `s2` is the variance the update
integrates against (with a count's dispersion already in it).

Censored values use the closed form the update does (`_censored_moments`): the
latent response is `N(ηbar, s2 + sd^2)`, a value recorded exactly is a point of
it, one at a limit or recorded to an interval (`_response_interval`) is that
interval. The other kinds sum the probabilities the update's own rule gives
each outcome (`_binary_moments`), normalised over a finite support.
"""
function _ctsem_item_cdf(ηbar::T, s2::T, y::Real, thresholds, kind::Int, nodes,
    weights) where {T}
    s2 = max(s2, zero(T))
    if kind == CTSEM_OBS_CENSORED
        lower, upper, sd = _censor_limits(thresholds, T)
        half, power = _censor_response(thresholds, T)
        rootV = sqrt(sd * sd + s2)
        atlower = _censored_at(y, lower, false)
        atupper = _censored_at(y, upper, true)
        if !atlower && !atupper && !(_primal(half) > 0)
            z = (_response_transform(T(y), _response_origin(lower), power) - ηbar) / rootV
            u = exp(_norm_logcdf(z))
            return (u, u)
        end
        lo, hi = _response_interval(y, lower, upper, half, power, atlower, atupper)
        Flo = lo === nothing ? zero(T) : exp(_norm_logcdf((lo - ηbar) / rootV))
        # One minus the upper tail keeps the digits a value near the top needs.
        Fhi = hi === nothing ? one(T) : -expm1(_norm_logcdf(-(hi - ηbar) / rootV))
        return (Flo, Fhi)
    end
    s = sqrt(s2)
    prob(v) = exp(first(_binary_moments(ηbar, s, v, nodes, weights, thresholds, kind)))
    if kind == CTSEM_OBS_BINARY
        p0 = prob(0); p1 = prob(1)
        f0 = p0 / (p0 + p1)
        return y > 0.5 ? (f0, one(T)) : (zero(T), f0)
    elseif kind == CTSEM_OBS_ORDINAL
        k = Int(y)
        below = zero(T); at = zero(T); total = zero(T)
        for v in 1:(length(thresholds) + 1)
            p = prob(v)
            total += p
            v < k && (below += p)
            v == k && (at = p)
        end
        return (below / total, (below + at) / total)
    elseif kind == CTSEM_OBS_COUNT
        k = Int(round(y))
        below = zero(T)
        for v in 0:(k - 1)
            below += prob(v)
        end
        return (min(below, one(T)), min(below + prob(k), one(T)))
    end
    return (T(NaN), T(NaN))
end
