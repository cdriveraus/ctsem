################################################################################
# SAEM: the exact marginal posterior mode, by stochastic approximation
################################################################################
#
# Why this exists. The Laplace route maximises an approximation to each unit's
# integral over its random effects, and that approximation fails where a
# unit's posterior is not near-Gaussian: at a flat-topped posterior the
# curvature at the mode goes to zero, `-logdet(M)/2` grows without bound, and
# an optimiser climbs the error rather than the likelihood. SAEM does not
# approximate the integral. It samples each unit's random effects from their
# conditional distribution given the data and the current parameters, and moves
# the parameters along the Fisher identity,
#
#     d/dtheta log p(y | theta) = E[ d/dtheta log p(y, u | theta) | y, theta ],
#
# estimated at the draws -- so its fixed point is the mode of the exact
# marginal posterior, whatever shape each unit's posterior has. The states
# stay integrated by the filter, as on the Laplace route
# (Rao-Blackwellisation); only the random effects are sampled. See
# CT-SEM/review/LAPLACE-nested-spike-2026-10-02.md.
#
# The draws. Metropolis-within-Gibbs over each unit's block tree, the tree
# `CTSEMLaplaceUnits` already lays out, in `chains` independent chains per
# unit:
#   * every leaf block (a subject) takes a random-walk move shaped by its
#     conditional curvature, clipped at the prior's so no proposal is wider
#     than the prior in any direction. Only that block's members see it, so a
#     proposal costs their filters and nothing else, and the leaves of a unit
#     are conditionally independent given their ancestors -- they move in
#     parallel;
#   * every block above a leaf (a study, a region) takes a collapsed move: it
#     shifts by Delta and each block beneath it by its Gaussian linear
#     response `A_c Delta = -clip(M_cc)^-1 M_ca Delta`. A fixed linear shear
#     has unit Jacobian and is reversed by -Delta, so the proposal is
#     symmetric; it moves a study together with its people. Its members'
#     filters run in parallel.
# Chains: a model with few units has few independent draws behind each
# population parameter, so it is replicated, as Monolix replicates a small
# data set, until units times chains reaches fifty. They run in parallel with
# the units, which also fills the workers a few-unit model leaves idle.
#
# The parameters. `theta += gamma .* (P \ g)`, with `g` the complete-data score
# at the draws averaged over chains, plus the prior's gradient -- the fixed-u
# sweep `_laplace_floored_unit_gradient!` takes -- and `P` the averaged outer
# product of the per-member complete-data scores plus the prior precision: the
# complete-data information, larger than the marginal one wherever the random
# effects carry missing information.
#
# No phases. Each parameter has its own step size, by Kesten's rule (1958;
# Delyon and Juditsky 1993 for the multivariate case): `gamma_i = (1 +
# flips_i)^(-2/3)`, where `flips_i` counts the sign changes of that
# parameter's step. A parameter still travelling keeps its sign and its full
# step; one oscillating about its optimum flips and has its step shrink, at
# the rate (exponent 2/3, inside Robbins-Monro's (1/2, 1]) the SAEM literature
# uses. So the move from "getting there" to "averaging" happens per parameter,
# on that parameter's own evidence, with nothing to tune and no switch to
# misfire -- the plateau test this replaces fired a thousand iterations early
# on the SNSF pilot. The estimate is the average of the last half of the
# iterates (suffix averaging), which forgets the transient without being told
# where it ended.
#
# Stopping, the one check that matters: when the estimate has stopped moving
# relative to its own uncertainty. The averages of the last two quarters of
# the run differ by `Delta`; the run stops once `Delta' P Delta / npar < 0.01`,
# a root-mean-square drift below a tenth of a standard error per parameter in
# the complete-data metric, which overstates the information and so errs
# strict. Monte Carlo noise in the averages keeps that quantity up, so a noisy
# problem runs longer on its own. The averaged point is then handed to the
# optimiser's finish and certification (`ctsem_optimize`).
#
# Reproducibility. Every move draws from a stream seeded by (seed, iteration,
# unit, chain, block), never from a shared generator, and every quantity a task
# writes belongs to its own unit and chain, so the draws do not depend on how
# the work is split. Member scores land in per-unit, per-chain matrices and are
# reduced in a fixed order after the join, so a fixed seed and `cores`
# reproduce the run exactly.
#
# Constants that remain, each a choice of scale rather than of phase: the step
# cap (0.25 raw, a trust radius in coordinates whose priors are standard
# normal), the information's averaging rate (0.1), the shape refresh interval
# (25), two sweeps and two collapsed moves per block per iteration, fifty
# units' worth of chains, and the tenth of a standard error.

mutable struct CTSEMSAEMState
    theta::Vector{Float64}
    # Per unit, per chain: the latent vector and each member's log likelihood
    # at it.
    u::Vector{Vector{Vector{Float64}}}
    ll::Vector{Vector{Vector{Float64}}}
    # Per unit and block, shared by its chains: the proposal's Cholesky factor
    # (lower), and for a block with ancestors its linear response to each of
    # them, aligned with `blocks[b].ancestors`.
    chol::Vector{Vector{Matrix{Float64}}}
    response::Vector{Vector{Vector{Matrix{Float64}}}}
    # Per unit, per chain, per block: proposal scale and its counts. Per chain
    # so that two chains of one unit, running at once, never write the same
    # number.
    logscale::Vector{Vector{Vector{Float64}}}
    accepted::Vector{Vector{Vector{Int}}}
    proposed::Vector{Vector{Vector{Int}}}
    # Per unit: the blocks with no descendants, the rest innermost first, and
    # each block's descendants.
    leaves::Vector{Vector{Int}}
    uppers::Vector{Vector{Int}}
    descendants::Vector{Vector{Vector{Int}}}
    # Per unit, per chain: the members' complete-data scores at the last M-step.
    scores::Vector{Vector{Matrix{Float64}}}
    info::Matrix{Float64}
    # Kesten's counts and the last step's signs, per parameter.
    flips::Vector{Int}
    lastsign::Vector{Int8}
    # Cumulative sums of the iterates, for averages over any window.
    csum::Vector{Vector{Float64}}
    chains::Int
    iteration::Int
    seed::UInt64
end

_saem_rng(st::CTSEMSAEMState, tags...) = Random.Xoshiro(hash((st.seed, st.iteration, tags...)))

"""
    _saem_member_ll(laplace, U, m, theta, Ls, u)

Member `m` of unit `U`'s log likelihood at the unit's latent vector `u`: the
same filter the inner objective runs, value only, on this task's workspace.
`-Inf` where the point cannot be evaluated; an interrupt or a bug propagates.
"""
function _saem_member_ll(laplace::CTSEMLaplaceObjective, U::Integer, m::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, u::Vector{Float64})
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    shift = _laplace_scratch_vector!(laplace, Float64, length(theta), :saem_shift)
    shifted = _laplace_member_values!(shift, theta, laplace.spec, Ls, u,
        laplace.units.offsets[U][m])
    so = laplace.objective.subject_objectives[laplace.units.members[U][m]]
    value = try
        _extended_kalman_filter_continuous!(aws.ekf_ws, shifted, so.data,
            so.timesteps, so.params, so.tdpreds, so.tipreds, so.subject,
            so.max_timestep)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        -Inf
    end
    return isfinite(value) ? Float64(value) : -Inf
end

"""Member positions' log likelihoods into `dest`, in parallel over members."""
function _saem_members_ll!(dest::Vector{Float64}, laplace::CTSEMLaplaceObjective,
    U::Integer, positions::AbstractVector{<:Integer}, theta::Vector{Float64},
    Ls::Vector{Matrix{Float64}}, u::Vector{Float64})
    n = length(positions)
    if n <= 2
        for (j, m) in enumerate(positions)
            dest[j] = _saem_member_ll(laplace, U, m, theta, Ls, u)
        end
        return dest
    end
    _laplace_partition(laplace, n) do mine, _w
        for j in mine
            dest[j] = _saem_member_ll(laplace, U, positions[j], theta, Ls, u)
        end
        nothing
    end
    return dest
end

function _saem_clip(P::AbstractMatrix)
    E = _ctsem_symeig(Matrix{Float64}(P))
    return E.vectors * Diagonal(max.(E.values, 1.0)) * transpose(E.vectors)
end

# Lower Cholesky factor of the inverse of a precision clipped at the prior's.
function _saem_proposal_factor(P::AbstractMatrix)
    E = _ctsem_symeig(Matrix{Float64}(P))
    C = E.vectors * Diagonal(1 ./ max.(E.values, 1.0)) * transpose(E.vectors)
    k = size(C, 1)
    F = _ctsem_cholesky(Matrix(_laplace_symmetrise(C)), k)
    return Matrix(LowerTriangular(transpose(UpperTriangular(F.U[1:k, 1:k]))))
end

"""
    _saem_refresh!(st, laplace, U, Ls)

Proposal shapes for unit `U`, shared by its chains, from its curvature at the
first chain's draw: a leaf's from its own diagonal block, a block above from
its diagonal less what its descendants explain (each through its own clipped
diagonal), and every block's linear response to each of its ancestors. A
curvature that is not finite keeps the shapes it had.
"""
function _saem_refresh!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, Ls::Vector{Matrix{Float64}})
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return st
    M = _laplace_unit_curvature(laplace, U, st.theta, Ls, st.u[U][1])
    all(d -> all(isfinite, d), M.diag) || return st
    all(row -> all(c -> all(isfinite, c), row), M.coupling) || return st
    clipped = [_saem_clip(M.diag[b]) for b in eachindex(blocks)]
    for (b, blk) in enumerate(blocks)
        st.response[U][b] = [-(clipped[b] \ Matrix(M.coupling[b][t]))
                             for t in eachindex(blk.ancestors)]
    end
    for b in eachindex(blocks)
        P = Matrix{Float64}(M.diag[b])
        for c in st.descendants[U][b]
            t = findfirst(==(b), blocks[c].ancestors)
            B = Matrix(M.coupling[c][t])
            P .-= transpose(B) * (clipped[c] \ B)
        end
        st.chol[U][b] = _saem_proposal_factor(_laplace_symmetrise(P))
    end
    return st
end

_saem_target_acceptance(k::Integer) = 0.234 + 0.2 / k

function _saem_record_move!(st::CTSEMSAEMState, U, c, b, accepted::Bool, adapt::Float64)
    st.proposed[U][c][b] += 1
    st.accepted[U][c][b] += accepted
    k = size(st.chol[U][b], 1)
    st.logscale[U][c][b] += adapt * (accepted - _saem_target_acceptance(k))
    return nothing
end

"""
    _saem_sweep!(st, laplace, U, c, Ls, sweep, adapt; nupper)

One sweep of chain `c` of unit `U`: every leaf once, in parallel, then
`nupper` collapsed moves of every block above, innermost first.
"""
function _saem_sweep!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, c::Integer, Ls::Vector{Matrix{Float64}}, sweep::Integer,
    adapt::Float64; nupper::Integer=2)
    blocks = laplace.units.blocks[U]
    u = st.u[U][c]; ll = st.ll[U][c]; theta = st.theta
    _laplace_parallel(laplace, st.leaves[U]) do b
        blk = blocks[b]
        cols = (blk.offset + 1):(blk.offset + blk.size)
        rng = _saem_rng(st, U, c, b, sweep, 1)
        old = u[cols]
        step = exp(st.logscale[U][c][b]) .* (st.chol[U][b] * randn(rng, blk.size))
        # Writing `u` from several tasks at once is safe here: each leaf owns
        # its columns, and a member reads only its own leaf and the ancestors,
        # which no leaf move touches.
        u[cols] .= old .+ step
        newll = [_saem_member_ll(laplace, U, m, theta, Ls, u) for m in blk.members]
        logr = sum(newll) - sum(@view ll[blk.members]) -
            (sum(abs2, @view u[cols]) - sum(abs2, old)) / 2
        accepted = isfinite(logr) && log(rand(rng)) < logr
        accepted ? (ll[blk.members] .= newll) : (u[cols] .= old)
        _saem_record_move!(st, U, c, b, accepted, adapt)
        true
    end
    for a in st.uppers[U]
        blk = blocks[a]
        acols = (blk.offset + 1):(blk.offset + blk.size)
        moved = vcat(collect(acols), [collect((blocks[d].offset + 1):(blocks[d].offset +
            blocks[d].size)) for d in st.descendants[U][a]]...)
        newll = Vector{Float64}(undef, length(blk.members))
        for r in 1:nupper
            rng = _saem_rng(st, U, c, a, sweep, 2, r)
            old = u[moved]
            delta = exp(st.logscale[U][c][a]) .* (st.chol[U][a] * randn(rng, blk.size))
            u[acols] .+= delta
            for d in st.descendants[U][a]
                t = findfirst(==(a), blocks[d].ancestors)
                dcols = (blocks[d].offset + 1):(blocks[d].offset + blocks[d].size)
                u[dcols] .+= st.response[U][d][t] * delta
            end
            _saem_members_ll!(newll, laplace, U, blk.members, theta, Ls, u)
            logr = sum(newll) - sum(@view ll[blk.members]) -
                (sum(abs2, @view u[moved]) - sum(abs2, old)) / 2
            accepted = isfinite(logr) && log(rand(rng)) < logr
            accepted ? (ll[blk.members] .= newll) : (u[moved] .= old)
            _saem_record_move!(st, U, c, a, accepted, adapt)
        end
    end
    return st
end

"""
    _saem_unit_scores!(st, laplace, U, c, Ls, dL, positions)

Each member's complete-data score at chain `c`'s draw of unit `U`, into the
columns of `st.scores[U][c]`: one reverse sweep per member and the chain rule
through the level loadings, the sweep `_laplace_floored_unit_gradient!` takes
at a mode. In parallel over members. False where a member cannot be evaluated.
"""
function _saem_unit_scores!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, c::Integer, Ls::Vector{Matrix{Float64}}, dL, positions)
    members = laplace.units.members[U]
    npar = length(st.theta)
    u = st.u[U][c]
    S = st.scores[U][c]
    ok = Threads.Atomic{Bool}(true)
    _laplace_partition(laplace, length(members)) do mine, _w
        aws = _laplace_workspace!(laplace, Float64, npar)
        grad = _laplace_scratch_vector!(laplace, Float64, npar, :saem_grad)
        shift = _laplace_scratch_vector!(laplace, Float64, npar, :saem_shift)
        s = _laplace_scratch_vector!(laplace, Float64, npar, :saem_score)
        for m in mine
            offsets = laplace.units.offsets[U][m]
            shifted = _laplace_member_values!(shift, st.theta, laplace.spec, Ls, u, offsets)
            value = _laplace_subject_value_gradient!(grad,
                laplace.objective.subject_objectives[members[m]], aws, shifted)
            if !(isfinite(value) && all(isfinite, grad))
                ok[] = false
                return nothing
            end
            copyto!(s, grad)
            _laplace_chol_chain!(s, grad, laplace.spec, dL, positions, u, offsets)
            @views S[:, m] .= s
        end
        nothing
    end
    return ok[]
end

"""Fifty units' worth of chains: `cld(50, nunits)`, between one and eight."""
_saem_default_chains(nunits::Integer) = clamp(cld(50, max(nunits, 1)), 1, 8)

"""
    ctsem_saem_init(laplace, theta; seed, chains)

SAEM's state at `theta`: every chain of each unit starts at the unit's Laplace
mode (from the origin, as the objective solves it), and the proposal shapes
come from the curvature there. `chains = 0` takes `_saem_default_chains`.
"""
function ctsem_saem_init(laplace::CTSEMLaplaceObjective, theta::AbstractVector;
    seed::Integer=1, chains::Integer=0)
    x = collect(Float64, theta)
    _laplace_check_indices(laplace, length(x))
    ctsem_laplace_evaluate(laplace, x; gradient=false)
    nunits = length(laplace.units.members)
    K = chains > 0 ? Int(chains) : _saem_default_chains(nunits)
    npar = length(x)
    Ls = _laplace_popchols(x, laplace.spec)
    blocks = laplace.units.blocks
    descendants = [[[d for d in eachindex(blocks[U]) if b in blocks[U][d].ancestors]
                    for b in eachindex(blocks[U])] for U in 1:nunits]
    leaves = [[b for b in eachindex(blocks[U]) if isempty(descendants[U][b])]
              for U in 1:nunits]
    # Innermost first: a block with more ancestors is deeper.
    uppers = [sort([b for b in eachindex(blocks[U]) if !isempty(descendants[U][b])];
                   by=b -> -length(blocks[U][b].ancestors)) for U in 1:nunits]
    nmem(U) = length(laplace.units.members[U])
    st = CTSEMSAEMState(x,
        [[copy(laplace.modes[U]) for _ in 1:K] for U in 1:nunits],
        [[zeros(nmem(U)) for _ in 1:K] for U in 1:nunits],
        [[Matrix{Float64}(I, b.size, b.size) for b in blocks[U]] for U in 1:nunits],
        [[Matrix{Float64}[] for _ in blocks[U]] for U in 1:nunits],
        [[[log(2.38 / sqrt(b.size)) for b in blocks[U]] for _ in 1:K] for U in 1:nunits],
        [[zeros(Int, length(blocks[U])) for _ in 1:K] for U in 1:nunits],
        [[zeros(Int, length(blocks[U])) for _ in 1:K] for U in 1:nunits],
        leaves, uppers, descendants,
        [[zeros(npar, nmem(U)) for _ in 1:K] for U in 1:nunits],
        zeros(npar, npar), zeros(Int, npar), zeros(Int8, npar),
        Vector{Float64}[], K, 0, UInt64(seed))
    ok = _laplace_parallel(laplace, 1:nunits) do U
        _saem_refresh!(st, laplace, U, Ls)
        for c in 1:K
            _saem_members_ll!(st.ll[U][c], laplace, U, eachindex(laplace.units.members[U]),
                x, Ls, st.u[U][c])
        end
        all(c -> all(isfinite, st.ll[U][c]), 1:K)
    end
    ok || throw(ArgumentError("SAEM cannot start: a unit's log likelihood is not " *
        "finite at its Laplace mode"))
    return st
end

export ctsem_saem_init

"""The complete-data log posterior at the state's draws, averaged over chains."""
function _saem_logpost(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective)
    total = 0.0
    for U in eachindex(st.u), c in 1:st.chains
        total += sum(st.ll[U][c]; init=0.0) - sum(abs2, st.u[U][c]; init=0.0) / 2
    end
    return total / st.chains + _ctsem_log_prior(laplace.objective, st.theta)
end

function _saem_acceptance(st::CTSEMSAEMState)
    a = 0; p = 0
    for U in eachindex(st.accepted), c in 1:st.chains
        a += sum(st.accepted[U][c]; init=0); p += sum(st.proposed[U][c]; init=0)
    end
    return p == 0 ? NaN : a / p
end

"""Prior precision on the raw parameters, as a vector."""
function _saem_prior_precision(laplace::CTSEMLaplaceObjective, npar::Integer)
    prec = zeros(npar)
    obj = laplace.objective
    for j in eachindex(obj.prior_index)
        prec[obj.prior_index[j]] += obj.prior_weight / obj.prior_scale[j]^2
    end
    return prec
end

"""
    _saem_curvature(st, prec)

`P`, the M-step's metric: the averaged complete-data information plus the
prior precision, with a ridge relative to its own diagonal.
"""
_saem_curvature(st::CTSEMSAEMState, prec::Vector{Float64}) =
    Symmetric(st.info) + Diagonal(prec .+ 1e-8 .* (1 .+ diag(st.info)))

"""
    ctsem_saem_step!(st, laplace; sweeps, nupper, maxstep, refresh, info_rate,
        adapt, prec)

One SAEM iteration: `sweeps` E-step sweeps of every chain of every unit (units
and chains in parallel, and within one its leaves and its members), then the
M-step with each parameter's Kesten step size. Returns the complete-data log
posterior before the step, the score norm, the largest coordinate of the step
taken and the median step size.
"""
function ctsem_saem_step!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective;
    sweeps::Integer=2, nupper::Integer=2, maxstep::Real=0.25,
    refresh::Integer=25, info_rate::Real=0.1, adapt::Real=-1.0,
    prec::Vector{Float64}=_saem_prior_precision(laplace, length(st.theta)))
    st.iteration += 1
    k = st.iteration
    nunits = length(laplace.units.members)
    npar = length(st.theta)
    K = st.chains
    rate = adapt < 0 ? 1.0 / (1 + k)^0.6 : Float64(adapt)
    Ls = _laplace_popchols(st.theta, laplace.spec)
    if k % refresh == 0
        _laplace_parallel(laplace, 1:nunits) do U
            isempty(st.u[U][1]) || _saem_refresh!(st, laplace, U, Ls)
            true
        end
    end
    items = [(U, c) for U in 1:nunits for c in 1:K]
    ok = _laplace_parallel(laplace, items) do item
        U, c = item
        # The members' log likelihoods move with theta -- in a unit with no
        # random effects too, whose members still count in the log posterior.
        _saem_members_ll!(st.ll[U][c], laplace, U, eachindex(laplace.units.members[U]),
            st.theta, Ls, st.u[U][c])
        all(isfinite, st.ll[U][c]) || return false
        isempty(st.u[U][c]) && return true
        for sweep in 1:sweeps
            _saem_sweep!(st, laplace, U, c, Ls, sweep, rate; nupper=nupper)
        end
        true
    end
    ok || throw(DomainError(st.theta,
        "SAEM: a unit's log likelihood is not finite at the current parameters"))
    # Between the halves as well as after them: on a large nested model one
    # iteration's E-step takes seconds, and Escape should not wait for the
    # M-step's sweeps too.
    _ctsem_interrupt_check()
    dL = _laplace_level_chol_derivatives(st.theta, laplace.spec)
    positions = [_laplace_level_positions(laplace.spec, l)
                 for l in eachindex(laplace.spec.levels)]
    ok = _laplace_parallel(laplace, items) do item
        U, c = item
        isempty(laplace.units.members[U]) && return true
        _saem_unit_scores!(st, laplace, U, c, Ls, dL, positions)
    end
    ok || throw(DomainError(st.theta,
        "SAEM: a member's score is not finite at the current parameters"))
    # Reduced in unit and chain order, after the join: the same sum whatever
    # the split.
    g = zeros(npar)
    B = zeros(npar, npar)
    for U in 1:nunits, c in 1:K
        S = st.scores[U][c]
        size(S, 2) == 0 && continue
        for m in axes(S, 2)
            @views g .+= S[:, m]
        end
        BLAS.syrk!('U', 'N', 1.0 / K, S, 1.0, B)
    end
    g ./= K
    logpost = _saem_logpost(st, laplace)
    _ctsem_log_prior_gradient!(g, laplace.objective, st.theta)
    B = Matrix(Symmetric(B, :U))
    w = k == 1 ? 1.0 : Float64(info_rate)
    st.info .= (1 - w) .* st.info .+ w .* B
    P = _saem_curvature(st, prec)
    F = cholesky(P; check=false)
    d = issuccess(F) ? F \ g : g ./ (diag(P) .+ 1.0)
    # Kesten, per parameter: a sign change of this parameter's step is the
    # evidence that it has reached its optimum and is oscillating about it.
    step = similar(d)
    for i in eachindex(d)
        s = Int8(sign(d[i]))
        (s != 0 && st.lastsign[i] != 0 && s != st.lastsign[i]) && (st.flips[i] += 1)
        s != 0 && (st.lastsign[i] = s)
        gamma = (1 + st.flips[i])^(-2 / 3)
        step[i] = clamp(gamma * d[i], -maxstep, maxstep)
    end
    st.theta .+= step
    push!(st.csum, isempty(st.csum) ? copy(st.theta) : st.csum[end] .+ st.theta)
    return (logpost=logpost, gradient_norm=norm(g),
        step=maximum(abs, step; init=0.0),
        gamma=sort((1 .+ st.flips) .^ (-2 / 3))[cld(npar, 2)])
end

export ctsem_saem_step!

"""Average of the iterates `a+1:b` from the cumulative sums."""
_saem_window(st::CTSEMSAEMState, a::Integer, b::Integer) =
    (a == 0 ? st.csum[b] : st.csum[b] .- st.csum[a]) ./ (b - a)

"""
    _saem_drift(st, P)

How far the estimate moved between the last two quarters of the run, as the
mean squared change per parameter in the metric `P`, in units of squared
standard errors; `NaN` before forty iterations, when a quarter is too short
to average anything.
"""
function _saem_drift(st::CTSEMSAEMState, P)
    k = length(st.csum)
    k < 40 && return NaN
    q = k ÷ 4
    delta = _saem_window(st, k - q, k) .- _saem_window(st, k - 2q, k - q)
    return dot(delta, P * delta) / length(delta)
end

"""
    ctsem_saem(laplace, start; maxiter, tol, seed, chains, ...)

SAEM from `start` until the estimate -- the average of the last half of the
iterates -- has stopped moving (`_saem_drift` below `tol`, a squared tenth of a
standard error per parameter by default), or `maxiter` iterations. Every
parameter's step follows Kesten's rule (see the file's header); there are no
phases. Returns the estimate (`minimizer`), the iterations, whether it stopped
on the drift rule (`settled`) and the drift there, the chains, the mean
acceptance rate and a trace. `progress`, `progress_*` and `callback` behave as
on `ctsem_optimize`; the callback receives the iteration, `maxiter`, the
complete-data log posterior, the score norm and the current estimate.
"""
function ctsem_saem(laplace::CTSEMLaplaceObjective, start::AbstractVector;
    maxiter::Integer=10000, tol::Real=0.01, seed::Integer=1, chains::Integer=0,
    sweeps::Integer=2, nupper::Integer=2, maxstep::Real=0.25, refresh::Integer=25,
    info_rate::Real=0.1, progress::Bool=false, progress_overwrite::Bool=true,
    progress_sink=nothing, progress_every::Real=0.0, callback=nothing)
    maxiter >= 1 || throw(ArgumentError("SAEM needs at least one iteration"))
    st = ctsem_saem_init(laplace, start; seed=seed, chains=chains)
    prec = _saem_prior_precision(laplace, length(st.theta))
    reporter = _ctsem_progress_reporter(progress, "saem", progress_overwrite,
        progress_sink, progress_every)
    watcher = CTSEMCallback(callback)
    trace = CTSEMTrace(:logpost_complete, :gradient_norm, :gamma, :step,
        :acceptance, :drift)
    settled = false
    drift = NaN
    for k in 1:Int(maxiter)
        out = ctsem_saem_step!(st, laplace; sweeps=sweeps, nupper=nupper,
            maxstep=maxstep, refresh=refresh, info_rate=info_rate, prec=prec)
        drift = _saem_drift(st, _saem_curvature(st, prec))
        acceptance = _saem_acceptance(st)
        _record!(trace, k, out.logpost, out.gradient_norm, out.gamma, out.step,
            acceptance, drift)
        estimate = _saem_window(st, k - max(1, k ÷ 2), k)
        if _due(reporter)
            _progress_optimise(reporter, k, Int(maxiter),
                @sprintf("logpost (complete) %11.2f", out.logpost),
                @sprintf("drift %8.2e", drift),
                @sprintf("median step %.2f", out.gamma),
                @sprintf("accept %.2f", acceptance))
        end
        _invoke_callback(watcher, k, Int(maxiter), out.logpost, out.gradient_norm,
            estimate)
        if isfinite(drift) && drift < tol
            settled = true
            break
        end
    end
    k = st.iteration
    minimizer = _saem_window(st, k - max(1, k ÷ 2), k)
    _progress_done(reporter, @sprintf("%d iterations", k),
        settled ? @sprintf("settled, drift %.2e", drift) :
            @sprintf("not settled, drift %.2e", drift),
        @sprintf("%d chain%s, accept %.2f", st.chains, st.chains == 1 ? "" : "s",
            _saem_acceptance(st)))
    return (minimizer=minimizer, iterations=k, settled=settled, drift=drift,
        chains=st.chains, acceptance=_saem_acceptance(st),
        trace=_trace_result(trace), state=st)
end

export ctsem_saem

# The phase `ctsem_optimize` runs before L-BFGS when asked: SAEM on a Laplace
# objective, nothing on any other route (the R side refuses those first).
_ctsem_saem_phase(objective, start; kwargs...) = nothing
_ctsem_saem_phase(laplace::CTSEMLaplaceObjective, start; kwargs...) =
    ctsem_saem(laplace, start; kwargs...)
