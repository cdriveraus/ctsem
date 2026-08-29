"""
The reverse pass for binary observations.

The forward update (see `binary_measurement.jl`) applies each binary
observation as an exact scalar conditioning:

    c = Pλ,  b = λ'Pλ,  a = λ'x + μ
    (logZ, m, v) = moments(a, b, y)   # m is the mean's offset from a
    x ← x + c m/b
    P ← P - c c' (1 - v/b)/b

so the reverse is the differential of that, with the moment function's own
partials supplied in closed form by `_binary_moment_derivatives`.

# Why the record holds the prior rather than the intermediates

Every quantity above is reconstructible from `(x, P, λ, μ, y)`, and the forward
pass overwrites `x` and `P` in place. Storing the prior and replaying is both
smaller and less brittle than storing six intermediates per observation per
row -- and it is what `_reverse_update!` already does for the Gaussian block.

# Why the Gaussian block needed no changes

The forward applies binary observations *before* the Gaussian ones, so if the
update record is taken after the binary block, its `state_in`/`P_in` are exactly
the Gaussian block's own inputs and its existing reverse is untouched. The
binary chain then reverses after it, from the true prior. Ordering the other way
would have meant teaching `_reverse_update!` about a preceding step.
"""

"""One row's binary observations, and the state they were applied to."""
mutable struct CTSEMBinaryRecord{T}
    rows::Vector{Int}         # manifest indices, in application order
    state_in::Vector{T}       # prior mean, before any binary conditioning
    P_in::Matrix{T}           # prior covariance, likewise
    Lambda::Matrix{T}         # LAMBDA[rows, :]
    manifestmeans::Vector{T}  # MANIFESTMEANS[rows]
    y::Vector{T}              # the observations themselves
    # Cumulated thresholds per observation, empty for a binary one. Copied
    # rather than referenced: `_ordinal_thresholds!` hands back a view into a
    # workspace scratch that the next observation overwrites.
    thresholds::Vector{Vector{T}}
end

# Deliberately untyped in `tape`: this file is included before the tape's own,
# because the tape holds a vector of these records and so needs the type to
# exist first. Annotating `::CTSEMAdjointTape` here would make that circular.
# The `::Nothing` method above is what keeps an untraced pass free of dispatch.
function _record_binary!(tape, ws, pars, data, obs_col, rows, state_in, P_in, n)
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
    T = eltype(state_in)
    n = Int(n)
    obs_col = Int(obs_col)
    index = (tape.nbinaries += 1)
    r = collect(rows)
    if index > length(tape.binaries)
        push!(tape.binaries, CTSEMBinaryRecord{T}(
            r, collect(vec(state_in)), Matrix(P_in),
            Matrix(pars.LAMBDA[r, 1:n]), collect(pars.MANIFESTMEANS[r]),
            T[data[i, obs_col] for i in r],
            Vector{T}[collect(T, _ordinal_thresholds!(ws, pars, i)) for i in r]))
    else
        record = tape.binaries[index]
        record.rows = r
        record.state_in = collect(vec(state_in))
        record.P_in = Matrix(P_in)
        record.Lambda = Matrix(pars.LAMBDA[r, 1:n])
        record.manifestmeans = collect(pars.MANIFESTMEANS[r])
        record.y = T[data[i, obs_col] for i in r]
        record.thresholds =
            Vector{T}[collect(T, _ordinal_thresholds!(ws, pars, i)) for i in r]
    end
    push!(tape.program, (:binary, index))
    return nothing
end

"""
    _reverse_binary!(x̄, P̄, θ̄ca, record, n)

Reverse one row's chain of binary conditionings.

Replays the chain forward first, keeping each step's `(x, P)` so the reverse can
be taken at the values that step actually saw, then walks it backwards.
"""
function _reverse_binary!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMBinaryRecord{T}, n::Int) where {T}
    k = length(record.rows)
    k == 0 && return nothing

    # Forward replay, storing the state each observation was applied to.
    states = Vector{Vector{T}}(undef, k)
    covs = Vector{Matrix{T}}(undef, k)
    x = copy(record.state_in)
    P = copy(record.P_in)
    for j in 1:k
        states[j] = copy(x)
        covs[j] = copy(P)
        λ = view(record.Lambda, j, :)
        a = record.manifestmeans[j]
        c = P * λ
        b = zero(T)
        @inbounds for i in 1:n
            b += λ[i] * c[i]
            a += λ[i] * x[i]
        end
        b <= zero(T) && continue
        g = _binary_moment_derivatives(a, b, record.y[j], record.thresholds[j])
        m, v = g[2], g[3]
        shift = m / b
        shrink = (one(T) - v / b) / b
        @inbounds for i in 1:n
            x[i] += c[i] * shift
        end
        @inbounds for jj in 1:n, ii in 1:n
            P[ii, jj] -= shrink * c[ii] * c[jj]
        end
    end

    # Reverse, last observation first.
    for j in k:-1:1
        λ = view(record.Lambda, j, :)
        x0 = states[j]
        P0 = covs[j]
        a = record.manifestmeans[j]
        c = P0 * λ
        b = zero(T)
        @inbounds for i in 1:n
            b += λ[i] * c[i]
            a += λ[i] * x0[i]
        end
        b <= zero(T) && continue
        τ = record.thresholds[j]
        g = _binary_moment_derivatives(a, b, record.y[j], τ)
        m, v = g[2], g[3]
        dlogZ_da, dlogZ_db = g[4], g[5]
        dm_da, dm_db = g[6], g[7]
        dv_da, dv_db = g[8], g[9]
        shift = m / b
        shrink = (one(T) - v / b) / b

        # Cotangents of the two things this step wrote.
        shiftbar = zero(T)
        @inbounds for i in 1:n
            shiftbar += x̄[i] * c[i]
        end
        cbar = Vector{T}(undef, n)
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

        # Thresholds, when this observation has any. THRESHOLDS holds gaps and
        # the forward pass cumulates them, so the cotangent on gap `i` is the
        # sum of the cotangents on every threshold at or after it.
        row = record.rows[j]
        if !isempty(τ)
            Jτ = _binary_threshold_derivatives(a, b, record.y[j], τ)
            running = zero(T)
            @inbounds for i in length(τ):-1:1
                running += logZbar * Jτ[1, i] + mbar * Jτ[2, i] + vbar * Jτ[3, i]
                θ̄ca.THRESHOLDS[row, i] += running
            end
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
