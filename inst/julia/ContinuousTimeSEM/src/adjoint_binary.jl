"""
The reverse pass for binary observations.

The forward update (see `binary_measurement.jl`) applies each binary
observation as an exact scalar conditioning:

    c = Pλ,  b = λ'Pλ,  a = λ'x + μ
    (logZ, m, v) = moments(a, b, y)   # m is the mean's offset from a
    x ← x + c m/b
    P ← P - c c' (1 - v/b)/b

so the reverse is the differential of that, with the moment function's own
partials supplied by `_binary_moment_jacobian`.

# Why the record holds what the forward saw, not only the prior

The moments are the expensive part: a 21-node rule per observation, and their
Jacobian is that rule again with its partials accumulated beside it. Recording
only the prior meant the reverse replayed the chain to find each observation's
`(x, P)` and evaluated the rule there, so a traced pass paid for the rule twice
per observation -- plainly in the forward, with partials in the replay. A
traced forward now evaluates the Jacobian itself (`_record_binary_step!`),
whose moments are the ones the plain rule would have returned, and keeps it
with the state the observation was applied to; the reverse reads them. One evaluation per
observation per pass instead of two: an ordinal gradient sweep about a third
cheaper on ord4 (local), and the rule's Jacobian is now taken at exactly the
point the forward used rather than at a replayed one that could differ from it
in the last bit. An untraced pass evaluates the plain rule, as before.

# Why the Gaussian block needed no changes

The forward applies binary observations *before* the Gaussian ones, so if the
update record is taken after the binary block, its `state_in`/`P_in` are exactly
the Gaussian block's own inputs and its existing reverse is untouched. The
binary chain then reverses after it, from the true prior. Ordering the other way
would have meant teaching `_reverse_update!` about a preceding step.
"""

"""One row's binary observations, and what the forward pass saw at each."""
mutable struct CTSEMBinaryRecord{T}
    rows::Vector{Int}         # manifest indices, in application order
    Lambda::Matrix{T}         # LAMBDA[rows, :]
    manifestmeans::Vector{T}  # MANIFESTMEANS[rows]
    y::Vector{T}              # the observations themselves
    # Cumulated thresholds per observation, empty for a binary one. Copied
    # rather than referenced: `_ordinal_thresholds!` hands back a view into a
    # workspace scratch that the next observation overwrites.
    thresholds::Vector{Vector{T}}
    # Which kind of observation each row is, in the same order. The forward
    # pass reads it from the workspace; the reverse has no workspace to read,
    # and inferring it from `thresholds` being empty would make a count look
    # like a Bernoulli -- the one confusion this whole argument exists to stop.
    kinds::Vector{Int}
    # Written by the forward pass as it applies each observation
    # (`_record_binary_step!`): how many it has applied, and for each the
    # state and covariance it was applied to, `c = P λ`, the predictor's mean
    # `a` and variance `b` (a count's dispersion included), and the moments'
    # offset, variance and Jacobian there. `states[1]` and `covs[1]` are the
    # row's prior. Sized on first use and reused by every later pass.
    applied::Int
    states::Vector{Vector{T}}
    covs::Vector{Matrix{T}}
    cs::Vector{Vector{T}}
    ab::Vector{Tuple{T,T}}
    moments::Vector{Tuple{T,T,Matrix{T}}}
    # The reverse pass's cotangent of `c`, one observation at a time.
    cbar::Vector{T}
end

"""
    _extras_cotangent!(θ̄ca, row, τ, kind, cot)

Push a row's extras cotangents back onto the matrix cells they came from.

`cot(i)` is the cotangent with respect to the i-th *assembled* extra -- what
`_ordinal_thresholds!` wrote -- and the cells are what the parameter vector
holds, which is not the same thing for any kind that accumulates. One function
because the two reverse passes below both need it and had a copy each, and the
copies were the place a new kind's rule could be added to one and not the
other.

The three rules:

  * censored, whose extras are two constant limits and a standard deviation
    that belongs to MANIFESTVAR;
  * binary with asymptotes, where cell 1 is `c` and cell 2 is a gap `g` with
    `d = c + (1-c)g`, so `dc/dcell1 = 1`, `dd/dcell1 = 1-g` and
    `dd/dcell2 = 1-c`;
  * ordinal, where threshold `k` is the sum of the first `k` cells, so a
    cell's cotangent is the running sum of every threshold at or above it.

Getting the second one wrong is not visible in a likelihood -- the forward pass
is untouched by it -- and shows up only as an optimizer walking off a cliff.
Measured before this existed, with the ordinal rule applied to a three
parameter logistic: the fit stopped at a gradient norm of 6536 where a
converged one is 1e-3, and every downstream diagnostic then described a point
that was not a mode.
"""
@inline function _extras_cotangent!(θ̄ca, row::Int, τ, kind::Int, cot::F
    ) where {F}
    if kind == CTSEM_OBS_CENSORED
        @inbounds θ̄ca.MANIFESTVAR[row, row] += cot(3)
        return nothing
    end
    if kind == CTSEM_OBS_BINARY
        @inbounds begin
            c = τ[1]
            room = one(c) - c
            g = room > zero(room) ? (τ[2] - c) / room : zero(c)
            d1 = cot(1)
            d2 = cot(2)
            θ̄ca.THRESHOLDS[row, 1] += d1 + d2 * (one(g) - g)
            θ̄ca.THRESHOLDS[row, 2] += d2 * room
        end
        return nothing
    end
    running = zero(cot(1))
    @inbounds for i in length(τ):-1:1
        running += cot(i)
        θ̄ca.THRESHOLDS[row, i] += running
    end
    return nothing
end

"""
    _record_binary!(tape, ws, pars, data, obs_col, rows, n)

Start the record of one row's categorical observations, `rows`, before the
forward pass applies them, and return it for `_record_binary_step!` to fill as
it does; `nothing` when nothing is being taped.
"""
function _record_binary!(tape, ws, pars, data, obs_col, rows, n)
    # Two different tracing mechanisms reach this: the adjoint tape, and the
    # Kalman trace that `ctKalman`/`ctPredict` use to collect filtered states.
    # Only the first wants a record, and only the first has anywhere to put one.
    # Tested by field rather than by type because this file is included before
    # the tape's own -- the tape holds a vector of these records, so the record
    # type has to exist first, and naming `CTSEMAdjointTape` here would make
    # that circular. Dropping the annotation entirely was the first attempt and
    # it turned a dispatch error into a `FieldError` on the Kalman trace, which
    # is worse: it surfaced in `ctKalman`, `ctPredict`, `ctLOO`,
    # `ctACFresiduals` and `ctPostPredPlots` all at once, far from the cause.
    tape === nothing && return nothing
    hasproperty(tape, :nbinaries) || return nothing
    isempty(rows) && return nothing
    # Only a model with categorical indicators gets here; see `_ctsem_barrier`.
    # The type asserted so the caller, which hands the record on to the
    # forward update, stays inferable.
    return _ctsem_barrier(_record_binary_rows!, tape, ws, pars, data, obs_col,
        rows, n)::CTSEMBinaryRecord{eltype(ws.state)}
end

function _record_binary_rows!(tape, ws, pars, data, obs_col, rows, n)
    T = eltype(ws.state)
    n = Int(n)
    obs_col = Int(obs_col)
    index = (tape.nbinaries += 1)
    if index > length(tape.binaries)
        push!(tape.binaries, CTSEMBinaryRecord{T}(Int[], zeros(T, 0, n), T[], T[],
            Vector{T}[], Int[], 0, Vector{T}[], Matrix{T}[], Vector{T}[],
            Tuple{T,T}[], Tuple{T,T,Matrix{T}}[], T[]))
    end
    # Refilled in place, as the tape's other records are (`_tape_fill!`): the
    # arrays of the last pass at this index already have the row's shape.
    record = tape.binaries[index]
    _tape_fill!(record.rows, rows)
    record.Lambda = _tape_gather!(record.Lambda, pars.LAMBDA, rows, n)
    _tape_gather!(record.manifestmeans, pars.MANIFESTMEANS, rows)
    resize!(record.y, length(rows))
    resize!(record.kinds, length(rows))
    length(record.thresholds) < length(rows) &&
        resize!(record.thresholds, length(rows))
    @inbounds for (k, i) in enumerate(rows)
        record.y[k] = data[i, obs_col]
        record.kinds[k] = Int(ws.manifesttype[i])
        τ = _ordinal_thresholds!(ws, pars, i)
        if isassigned(record.thresholds, k)
            _tape_fill!(record.thresholds[k], τ)
        else
            record.thresholds[k] = collect(T, τ)
        end
    end
    record.applied = 0
    push!(tape.program, (:binary, index))
    return record
end

"""
    _record_binary_step!(record, ws, c, a, b, y, nodes, weights, thresholds,
        kind, n)

The moments of the next observation of `record`'s row, `(logZ, offset,
variance)`, for a traced forward pass: the state and covariance it is applied
to (`ws.state`, `ws.P_predict`), `c`, `a` and `b` are kept, and where `b` is
above the variance floor so are the moments' offset, variance and Jacobian,
from `_binary_moment_jacobian`, whose moments are what `_binary_moments`
returns. At or below the floor the forward adds only the likelihood at `a` and
moves nothing, and the reverse differentiates that instead, so the plain rule
serves.
"""
function _record_binary_step!(record::CTSEMBinaryRecord{T}, ws, c, a::T, b::T,
    y::Real, nodes, weights, thresholds, kind::Int, n::Int) where {T}
    j = (record.applied += 1)
    if length(record.states) < j
        push!(record.states, Vector{T}(undef, n))
        push!(record.covs, Matrix{T}(undef, n, n))
        push!(record.cs, Vector{T}(undef, n))
        push!(record.ab, (a, b))
        push!(record.moments, (zero(T), zero(T), zeros(T, 0, 0)))
    end
    x = record.states[j]
    length(x) == n || resize!(x, n)
    P = record.covs[j]
    size(P) == (n, n) || (P = record.covs[j] = Matrix{T}(undef, n, n))
    cj = record.cs[j]
    length(cj) == n || resize!(cj, n)
    source = ws.P_predict.data
    @inbounds for q in 1:n
        x[q] = ws.state[q]
        cj[q] = c[q]
        for p in 1:n
            P[p, q] = source[p, q]
        end
    end
    record.ab[j] = (a, b)
    b > T(_CTSEM_MIN_VARIANCE[]) ||
        return _binary_moments(a, sqrt(b), y, nodes, weights, thresholds, kind)
    logZ, m, v, J = _binary_moment_jacobian(a, b, y, nodes, weights, thresholds,
        kind)
    record.moments[j] = (m, v, J)
    return (logZ, m, v)
end

"""
    _reverse_binary!(x̄, P̄, θ̄ca, record, n)

Reverse one row's chain of binary conditionings, last observation first, at
the state, covariance and moments the forward pass kept for each
(`_record_binary_step!`).
"""
function _reverse_binary!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMBinaryRecord{T}, n::Int) where {T}
    k = length(record.rows)
    k == 0 && return nothing
    # A reverse pass exists only for a forward pass that applied the whole row:
    # one that stopped on an observation it could not score was invalid.
    record.applied == k || error("binary record: the forward pass applied ",
        record.applied, " of ", k, " observations")
    cbar = record.cbar
    length(cbar) == n || resize!(cbar, n)

    for j in k:-1:1
        λ = view(record.Lambda, j, :)
        x0 = record.states[j]
        P0 = record.covs[j]
        c = record.cs[j]
        a, b = record.ab[j]
        τ = record.thresholds[j]
        row = record.rows[j]
        if !(b > T(_CTSEM_MIN_VARIANCE[]))
            # The forward left the state and covariance alone here and added
            # `log P(y | a)`. That term is not nothing: it depends on the
            # linear predictor, and so on the state, LAMBDA, MANIFESTMEANS and
            # the thresholds. Skipping the whole observation -- which is what
            # `continue` did -- drops all of it, and drops it exactly where a
            # variance is collapsing, which is where an optimiser most needs
            # the gradient to point somewhere sensible.
            score, _ = _category_score(a, record.y[j], τ, record.kinds[j])
            @inbounds for i in 1:n
                x̄[i] += score * λ[i]
                θ̄ca.LAMBDA[row, i] += score * x0[i]
            end
            θ̄ca.MANIFESTMEANS[row] += score
            # A count's extra is excluded: it is not a parameter of
            # `_category_loglikelihood` at all -- see `_count_dispersion` -- so
            # differentiating that against it gives zero, and the `else` below
            # would then write the zero into a THRESHOLDS matrix a count model
            # has no reason to have. Its cotangent rides `b`, which is
            # degenerate here.
            if !isempty(τ) && record.kinds[j] != CTSEM_OBS_COUNT
                dτ = _ctsem_nested_gradient(
                    t -> _category_loglikelihood(a, record.y[j], t,
                        record.kinds[j]),
                    collect(T, τ))
                _extras_cotangent!(θ̄ca, row, τ, record.kinds[j],
                    i -> dτ[i])
            end
            continue
        end
        m, v, J = record.moments[j]
        dlogZ_da, dlogZ_db = J[1, 1], J[1, 2]
        dm_da, dm_db = J[2, 1], J[2, 2]
        dv_da, dv_db = J[3, 1], J[3, 2]
        shift = m / b
        shrink = (one(T) - v / b) / b

        # Cotangents of the two things this step wrote.
        shiftbar = zero(T)
        @inbounds for i in 1:n
            shiftbar += x̄[i] * c[i]
        end
        @inbounds for i in 1:n
            cbar[i] = x̄[i] * shift
        end
        shrinkbar = zero(T)
        @inbounds for jj in 1:n, ii in 1:n
            shrinkbar -= P̄[ii, jj] * c[ii] * c[jj]
        end
        # P -= shrink c c' contributes -shrink (P̄ + P̄') c to c̄.
        #
        # Written out rather than as -2 shrink P̄ c, because P̄ is *not*
        # symmetric: `c = Pλ` two blocks below adds `c̄ λ'` to it, which is a
        # rank-one term with no reason to be. With one latent state that makes
        # no difference and every gradient test the binary path had used one;
        # with two -- a second latent, or a random effect under
        # `intoverpop='augmented'`, which expands the state by one per varying
        # parameter -- the adjoint came back about 1% wrong on DRIFT and
        # DIFFUSION and 12% on the random effect's own variance, against
        # forward mode. Small enough to look like quadrature error and quite
        # large enough to stop an optimiser short.
        @inbounds for i in 1:n
            acc = zero(T)
            for jj in 1:n
                acc += (P̄[i, jj] + P̄[jj, i]) * c[jj]
            end
            cbar[i] -= shrink * acc
        end

        # Through shift = m/b and shrink = 1/b - v/b². `m` is the posterior
        # mean's *offset* from `a` now, so the explicit `-a` that used to sit
        # here is gone: its derivative lives inside `dm_da` instead.
        mbar = shiftbar / b
        abar = zero(T)
        bbar = -shiftbar * m / (b * b)
        vbar = -shrinkbar / (b * b)
        bbar += shrinkbar * (-one(T) / (b * b) + T(2) * v / (b * b * b))

        # Through the moment function, including the likelihood it returned.
        logZbar = one(T)
        abar += logZbar * dlogZ_da + mbar * dm_da + vbar * dv_da
        bbar += logZbar * dlogZ_db + mbar * dm_db + vbar * dv_db

        # A count's dispersion enters only through `b`, so its cotangent is
        # `bbar` times `d b/d σ = 2σ` and nothing else moves: `b`'s own flow to
        # LAMBDA and to P below is unaltered, because the term `b` gained
        # depends on neither.
        if record.kinds[j] == CTSEM_OBS_COUNT && !isempty(τ)
            @inbounds θ̄ca.MANIFESTVAR[row, row] += bbar * T(2) * τ[1]
        end

        # Thresholds, when this observation has any. THRESHOLDS holds gaps and
        # the forward pass cumulates them, so the cotangent on gap `i` is the
        # sum of the cotangents on every threshold at or after it.
        if !isempty(τ) && record.kinds[j] != CTSEM_OBS_COUNT
            _extras_cotangent!(θ̄ca, row, τ, record.kinds[j],
                i -> logZbar * J[1, 2 + i] + mbar * J[2, 2 + i] +
                    vbar * J[3, 2 + i])
        end

        @inbounds for i in 1:n
            x̄[i] += abar * λ[i]
            θ̄ca.LAMBDA[row, i] += abar * x0[i]
        end
        θ̄ca.MANIFESTMEANS[row] += abar

        # b = λ'Pλ, taken directly rather than through c so nothing is counted
        # twice: c's own cotangent below carries only the state and covariance
        # updates.
        @inbounds for i in 1:n
            θ̄ca.LAMBDA[row, i] += T(2) * bbar * c[i]
            for jj in 1:n
                P̄[i, jj] += bbar * λ[i] * λ[jj]
            end
        end

        # c = Pλ.
        @inbounds for i in 1:n
            for jj in 1:n
                P̄[i, jj] += cbar[i] * λ[jj]
            end
            acc = zero(T)
            for jj in 1:n
                acc += P0[i, jj] * cbar[jj]
            end
            θ̄ca.LAMBDA[row, i] += acc
        end
    end
    return nothing
end
