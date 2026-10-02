################################################################################
# SAEM: the exact marginal posterior mode, by stochastic approximation
################################################################################
#
# Why this exists. The Laplace route maximises an approximation to each unit's
# integral over its random effects, and that approximation fails where a
# unit's posterior is not near-Gaussian: at a flat-topped posterior the
# curvature at the mode goes to zero, `-logdet(M)/2` grows without bound, and
# an optimiser climbs the error rather than the likelihood. The gated floor
# repairs that for units it can decompose, and not for wide nested units. SAEM
# does not approximate the integral at all. It samples each unit's random
# effects from their conditional distribution given the data and the current
# parameters, and moves the parameters along the Fisher identity,
#
#     d/dtheta log p(y | theta) = E[ d/dtheta log p(y, u | theta) | y, theta ],
#
# estimated at the draws -- so its fixed point is the mode of the exact
# marginal posterior, whatever shape each unit's posterior has. The states
# stay integrated by the filter, as on the Laplace route
# (Rao-Blackwellisation); only the random effects are sampled.
#
# Measured on the SNSF pilot (380 people in 13 studies, 319 parameters), where
# the Laplace fit climbed a near-singular person's spike: an iteration costs
# 0.25 s against the Laplace route's 14 s, and in about 20 minutes it reached a
# point 416 nats above the Laplace fit's iteration 290 *on the Laplace
# objective itself*, with no unit near-singular there. See
# CT-SEM/review/LAPLACE-nested-spike-2026-10-02.md and
# CT-SEM/review/stochopt-saem-design.md.
#
# The draws. Metropolis-within-Gibbs over each unit's block tree, the tree
# `CTSEMLaplaceUnits` already lays out:
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
#     symmetric; it moves a study together with its people, which a plain
#     Gibbs step on the study block cannot do when the two are tightly coupled
#     (and they are: every person's data inform the study coordinates). Its
#     members' filters run in parallel.
# Shapes are refreshed from the curvature at the current draw every `refresh`
# iterations; scales adapt toward an acceptance rate that falls with the block
# size, at a diminishing rate.
#
# The parameters. `theta += gamma P^-1 g`, with `g` the complete-data score at
# the draws plus the prior's gradient -- exactly the fixed-u sweep
# `_laplace_floored_unit_gradient!` takes -- and `P` the averaged outer product
# of the per-member complete-data scores plus the prior precision. That is the
# complete-data information, larger than the marginal one wherever the random
# effects carry missing information, so the steps are conservative there; the
# marginal alternative (outer products of averaged scores) was measured to
# diverge from a far start, its estimate dominated by the scores' mean. Steps
# are capped in raw units.
#
# The schedule. Step size one until the complete-data log posterior stops
# rising -- the mean of the last quarter of the burn-in no longer exceeds the
# quarter before by its own standard error -- then `gamma = j^-alpha` with
# Polyak averaging, for `averaging` iterations. PX-SAEM's expansion step
# (`_saem_expand!`) is available for the burn-in and off by default: it was
# measured worse on the model this was built for. The averaged point is handed
# to L-BFGS on the Laplace objective (`ctsem_optimize`), whose finish and
# certification apply unchanged: where every unit is regular at SAEM's point,
# Laplace is accurate there and the polish is short.
#
# Reproducibility. Every move draws from a stream seeded by (seed, iteration,
# unit, block), never from a shared generator, so the draws do not depend on
# how the work is split. Member scores land in a per-unit matrix and are
# reduced in a fixed order after the join, so a fixed seed and `cores`
# reproduce the run exactly; across widths the BLAS reduction and the chunked
# filter differ in the last digits, as every route here does.
#
# Constants, measured or chosen, not derived: the window (50), the averaging
# length, the step cap (0.25 raw), the information's averaging rate (0.1), the
# refresh interval (25), two sweeps and two collapsed moves per block per
# iteration. Set on the SNSF pilot subset and the engine's linear fixtures;
# none was tuned on a fit's estimates.

mutable struct CTSEMSAEMState
    theta::Vector{Float64}
    u::Vector{Vector{Float64}}
    # Each member's log likelihood at the unit's current draw.
    ll::Vector{Vector{Float64}}
    # Per unit and block: the proposal's Cholesky factor (lower), and for a
    # block with ancestors its linear response to each of them, aligned with
    # `blocks[b].ancestors`.
    chol::Vector{Vector{Matrix{Float64}}}
    response::Vector{Vector{Vector{Matrix{Float64}}}}
    logscale::Vector{Vector{Float64}}
    accepted::Vector{Vector{Int}}
    proposed::Vector{Vector{Int}}
    # Per unit: the blocks with no descendants, the rest innermost first, and
    # each block's descendants.
    leaves::Vector{Vector{Int}}
    uppers::Vector{Vector{Int}}
    descendants::Vector{Vector{Vector{Int}}}
    # Per unit: the members' complete-data scores at the last M-step.
    scores::Vector{Matrix{Float64}}
    info::Matrix{Float64}
    thetabar::Vector{Float64}
    nbar::Int
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

Proposal shapes for unit `U` from its curvature at the current draw: a leaf's
from its own diagonal block, a block above from its diagonal less what its
descendants explain (each through its own clipped diagonal), and every
block's linear response to each of its ancestors. A curvature that is not
finite keeps the shapes it had.
"""
function _saem_refresh!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, Ls::Vector{Matrix{Float64}})
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return st
    M = _laplace_unit_curvature(laplace, U, st.theta, Ls, st.u[U])
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

function _saem_record_move!(st::CTSEMSAEMState, U, b, accepted::Bool, adapt::Float64)
    st.proposed[U][b] += 1
    st.accepted[U][b] += accepted
    k = size(st.chol[U][b], 1)
    st.logscale[U][b] += adapt * (accepted - _saem_target_acceptance(k))
    return nothing
end

"""
    _saem_sweep!(st, laplace, U, Ls, sweep, adapt; nupper)

One sweep of unit `U`: every leaf once, in parallel, then `nupper` collapsed
moves of every block above, innermost first.
"""
function _saem_sweep!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, Ls::Vector{Matrix{Float64}}, sweep::Integer, adapt::Float64;
    nupper::Integer=2)
    blocks = laplace.units.blocks[U]
    u = st.u[U]; ll = st.ll[U]; theta = st.theta
    _laplace_parallel(laplace, st.leaves[U]) do b
        blk = blocks[b]
        cols = (blk.offset + 1):(blk.offset + blk.size)
        rng = _saem_rng(st, U, b, sweep, 1)
        old = u[cols]
        step = exp(st.logscale[U][b]) .* (st.chol[U][b] * randn(rng, blk.size))
        # Writing `u` from several tasks at once is safe here: each leaf owns
        # its columns, and a member reads only its own leaf and the ancestors,
        # which no leaf move touches.
        u[cols] .= old .+ step
        newll = [_saem_member_ll(laplace, U, m, theta, Ls, u) for m in blk.members]
        logr = sum(newll) - sum(@view ll[blk.members]) -
            (sum(abs2, @view u[cols]) - sum(abs2, old)) / 2
        accepted = isfinite(logr) && log(rand(rng)) < logr
        accepted ? (ll[blk.members] .= newll) : (u[cols] .= old)
        _saem_record_move!(st, U, b, accepted, adapt)
        true
    end
    for a in st.uppers[U]
        blk = blocks[a]
        acols = (blk.offset + 1):(blk.offset + blk.size)
        moved = vcat(collect(acols), [collect((blocks[c].offset + 1):(blocks[c].offset +
            blocks[c].size)) for c in st.descendants[U][a]]...)
        newll = Vector{Float64}(undef, length(blk.members))
        for r in 1:nupper
            rng = _saem_rng(st, U, a, sweep, 2, r)
            old = u[moved]
            delta = exp(st.logscale[U][a]) .* (st.chol[U][a] * randn(rng, blk.size))
            u[acols] .+= delta
            for c in st.descendants[U][a]
                t = findfirst(==(a), blocks[c].ancestors)
                ccols = (blocks[c].offset + 1):(blocks[c].offset + blocks[c].size)
                u[ccols] .+= st.response[U][c][t] * delta
            end
            _saem_members_ll!(newll, laplace, U, blk.members, theta, Ls, u)
            logr = sum(newll) - sum(@view ll[blk.members]) -
                (sum(abs2, @view u[moved]) - sum(abs2, old)) / 2
            accepted = isfinite(logr) && log(rand(rng)) < logr
            accepted ? (ll[blk.members] .= newll) : (u[moved] .= old)
            _saem_record_move!(st, U, a, accepted, adapt)
        end
    end
    return st
end

"""
    _saem_unit_scores!(st, laplace, U, Ls, dL, positions)

Each member's complete-data score at the unit's draw, into `st.scores[U]`'s
columns: one reverse sweep per member and the chain rule through the level
loadings, the sweep `_laplace_floored_unit_gradient!` takes at a mode. In
parallel over members. False where a member cannot be evaluated.
"""
function _saem_unit_scores!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, Ls::Vector{Matrix{Float64}}, dL, positions)
    members = laplace.units.members[U]
    npar = length(st.theta)
    u = st.u[U]
    S = st.scores[U]
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

"""
    ctsem_saem_init(laplace, theta; seed)

SAEM's state at `theta`: each unit's draw starts at its Laplace mode (from the
origin, as the objective solves it) and the proposal shapes come from the
curvature there.
"""
function ctsem_saem_init(laplace::CTSEMLaplaceObjective, theta::AbstractVector;
    seed::Integer=1)
    x = collect(Float64, theta)
    _laplace_check_indices(laplace, length(x))
    ctsem_laplace_evaluate(laplace, x; gradient=false)
    nunits = length(laplace.units.members)
    npar = length(x)
    Ls = _laplace_popchols(x, laplace.spec)
    blocks = laplace.units.blocks
    descendants = [[[c for c in eachindex(blocks[U]) if b in blocks[U][c].ancestors]
                    for b in eachindex(blocks[U])] for U in 1:nunits]
    leaves = [[b for b in eachindex(blocks[U]) if isempty(descendants[U][b])]
              for U in 1:nunits]
    # Innermost first: a block with more ancestors is deeper.
    uppers = [sort([b for b in eachindex(blocks[U]) if !isempty(descendants[U][b])];
                   by=b -> -length(blocks[U][b].ancestors)) for U in 1:nunits]
    st = CTSEMSAEMState(x, [copy(laplace.modes[U]) for U in 1:nunits],
        [zeros(length(laplace.units.members[U])) for U in 1:nunits],
        [[Matrix{Float64}(I, b.size, b.size) for b in blocks[U]] for U in 1:nunits],
        [[Matrix{Float64}[] for _ in blocks[U]] for U in 1:nunits],
        [[log(2.38 / sqrt(b.size)) for b in blocks[U]] for U in 1:nunits],
        [zeros(Int, length(blocks[U])) for U in 1:nunits],
        [zeros(Int, length(blocks[U])) for U in 1:nunits],
        leaves, uppers, descendants,
        [zeros(npar, length(laplace.units.members[U])) for U in 1:nunits],
        zeros(npar, npar), zeros(npar), 0, 0, UInt64(seed))
    ok = _laplace_parallel(laplace, 1:nunits) do U
        _saem_refresh!(st, laplace, U, Ls)
        _saem_members_ll!(st.ll[U], laplace, U, eachindex(laplace.units.members[U]),
            x, Ls, st.u[U])
        all(isfinite, st.ll[U])
    end
    ok || throw(ArgumentError("SAEM cannot start: a unit's log likelihood is not " *
        "finite at its Laplace mode"))
    return st
end

export ctsem_saem_init

"""
    _saem_expand!(st, laplace)

The parameter-expansion step of PX-SAEM (Lavielle and Meza 2007), between
the E-step and the M-step during burn-in. Each level's draws are recentred
and, on a reduced-rank level, rescaled, with the population parameters
moved to compensate, so that every member's parameters `theta + L u` -- and
so the likelihood -- are unchanged while the draws' empirical mean and
second moment are put back to the prior's zero and identity (each shrunk toward no
change by `r + 1` pseudo-groups, so a level with few groups is not jolted by
its own sampling noise):

  * location, every level: with `m` the mean draw over the level's groups,
    `theta[re] += L m` and `u -= m`;
  * scale, a reduced-rank level, whose loadings are raw values
    (`_laplace_poploading`): with `S = A A'` the centred draws' second
    moment, the loadings `R -> R A` and `u -> A^-1 u`. `A` is lower
    triangular, so `R A` keeps the loadings' structural zeros. A full-rank
    level's scales and correlations are not linear in its Cholesky factor,
    so it gets the location step only.

Why: with an effect per subject and plenty of data each, a subject's draw is
pinned by its own data relative to the population parameters, so the
ordinary M-step and the next E-step chase each other -- EM's slow rate,
fastest to see in a level's scale and in loadings that few groups inform.
The expansion moves the population parameters to where the draws say they
are in one step. Burn-in only: with a prior on the raw parameters the step
does not vanish exactly at the posterior mode, so the averaged phase, whose
fixed point must be that mode, runs without it. The proposal factors and
linear responses are carried into the new coordinates.

Off by default (`expand` on `ctsem_saem`), because it did not help where it
was meant to. On the SNSF pilot subset (380 people in 13 studies, rank 3 at
both levels; same start and seed, 3000 iterations, local), the Laplace
objective at the averaged point was -297809.9 without it, -297921.3 with it on
every level and -298038.6 with it on the subject level alone (`min_groups =
50`). On the engine's six-subject linear fixture it brought the average from
0.23 to 0.12 standard errors of the exact mode. Kept for models where the
scale directions are what is slow and the groups are many.
"""
function _saem_expand!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective;
    min_groups::Integer=0)
    spec = laplace.spec
    for (l, level) in enumerate(spec.levels)
        k = nrandomeffects(level)
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        sites = Tuple{Int,Int}[]
        for U in eachindex(st.u), (b, blk) in enumerate(laplace.units.blocks[U])
            blk.level == l && push!(sites, (U, b))
        end
        G = length(sites)
        (G == 0 || G < min_groups) && continue
        cols(U, b) = (laplace.units.blocks[U][b].offset + 1):(laplace.units.blocks[U][b].offset +
            laplace.units.blocks[U][b].size)
        m = zeros(r)
        for (U, b) in sites
            m .+= @view st.u[U][cols(U, b)]
        end
        # Shrunk toward no change by r + 1 pseudo-groups at the prior's mean
        # and identity: negligible at hundreds of groups, and it halves the
        # step's own noise at the dozen a study level may have.
        kappa = r + 1
        m ./= (G + kappa)
        L = _laplace_popchol(st.theta, level)
        shift = L * m
        for p in 1:k
            st.theta[level.re_index[p]] += shift[p]
        end
        for (U, b) in sites
            st.u[U][cols(U, b)] .-= m
        end
        (isreducedrank(level) && G > r) || continue
        S = zeros(r, r)
        for (U, b) in sites
            v = st.u[U][cols(U, b)]
            S .+= v * transpose(v)
        end
        S = (S .+ kappa .* Matrix{Float64}(I, r, r)) ./ (G + kappa)
        F = _ctsem_cholesky(Matrix(_laplace_symmetrise(S)), r)
        issuccess(F) || continue
        A = Matrix(LowerTriangular(transpose(UpperTriangular(F.U[1:r, 1:r]))))
        Ainv = inv(LowerTriangular(A))
        R = zeros(k, r)
        counter = 0
        for q in 1:r, p in q:k
            counter += 1
            R[p, q] = st.theta[level.load_index[counter]]
        end
        R = R * A
        counter = 0
        for q in 1:r, p in q:k
            counter += 1
            st.theta[level.load_index[counter]] = R[p, q]
        end
        for (U, b) in sites
            c = cols(U, b)
            st.u[U][c] = Ainv * st.u[U][c]
            st.chol[U][b] = Matrix(Ainv * st.chol[U][b])
        end
        # A block's response to an ancestor maps the ancestor's coordinates to
        # its own: new = (own level's A^-1) old (ancestor level's A).
        for U in eachindex(st.u)
            blocks = laplace.units.blocks[U]
            for (c, blk) in enumerate(blocks), t in eachindex(blk.ancestors)
                a = blk.ancestors[t]
                isempty(st.response[U][c]) && continue
                blk.level == l && (st.response[U][c][t] = Ainv * st.response[U][c][t])
                blocks[a].level == l && (st.response[U][c][t] = st.response[U][c][t] * A)
            end
        end
    end
    return st
end

"""The complete-data log posterior at the state's draws."""
_saem_logpost(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective) =
    sum(sum, st.ll) - sum(u -> sum(abs2, u; init=0.0), st.u) / 2 +
    _ctsem_log_prior(laplace.objective, st.theta)

function _saem_acceptance(st::CTSEMSAEMState)
    a = 0; p = 0
    for U in eachindex(st.accepted)
        a += sum(st.accepted[U]; init=0); p += sum(st.proposed[U]; init=0)
    end
    return p == 0 ? NaN : a / p
end

"""
    ctsem_saem_step!(st, laplace; gamma, sweeps, nupper, maxstep, refresh,
        info_rate, adapt)

One SAEM iteration: `sweeps` E-step sweeps of every unit (units in parallel,
and within a unit its leaves and its members), the expansion step if
`expand` (`_saem_expand!`), then the M-step. Returns the
complete-data log posterior before the step, the score norm, and the largest
coordinate of the step taken.
"""
function ctsem_saem_step!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective;
    gamma::Real=1.0, sweeps::Integer=2, nupper::Integer=2, maxstep::Real=0.25,
    refresh::Integer=25, info_rate::Real=0.1, adapt::Real=-1.0, expand::Bool=false,
    expand_min_groups::Integer=0)
    st.iteration += 1
    k = st.iteration
    nunits = length(laplace.units.members)
    npar = length(st.theta)
    rate = adapt < 0 ? 1.0 / (1 + k)^0.6 : Float64(adapt)
    Ls = _laplace_popchols(st.theta, laplace.spec)
    ok = _laplace_parallel(laplace, 1:nunits) do U
        # The members' log likelihoods move with theta -- in a unit with no
        # random effects too, whose members still count in the log posterior.
        _saem_members_ll!(st.ll[U], laplace, U, eachindex(laplace.units.members[U]),
            st.theta, Ls, st.u[U])
        all(isfinite, st.ll[U]) || return false
        isempty(st.u[U]) && return true
        (k % refresh == 0) && _saem_refresh!(st, laplace, U, Ls)
        for sweep in 1:sweeps
            _saem_sweep!(st, laplace, U, Ls, sweep, rate; nupper=nupper)
        end
        true
    end
    ok || throw(DomainError(st.theta,
        "SAEM: a unit's log likelihood is not finite at the current parameters"))
    # Between the halves as well as after them: on a large nested model one
    # iteration's E-step takes seconds, and Escape should not wait for the
    # M-step's sweeps too.
    _ctsem_interrupt_check()
    expand && _saem_expand!(st, laplace; min_groups=expand_min_groups)
    dL = _laplace_level_chol_derivatives(st.theta, laplace.spec)
    positions = [_laplace_level_positions(laplace.spec, l)
                 for l in eachindex(laplace.spec.levels)]
    ok = _laplace_parallel(laplace, 1:nunits) do U
        isempty(laplace.units.members[U]) && return true
        _saem_unit_scores!(st, laplace, U, Ls, dL, positions)
    end
    ok || throw(DomainError(st.theta,
        "SAEM: a member's score is not finite at the current parameters"))
    # Reduced in unit order, after the join: the same sum whatever the split.
    g = zeros(npar)
    B = zeros(npar, npar)
    for U in 1:nunits
        S = st.scores[U]
        size(S, 2) == 0 && continue
        for m in axes(S, 2)
            @views g .+= S[:, m]
        end
        BLAS.syrk!('U', 'N', 1.0, S, 1.0, B)
    end
    logpost = _saem_logpost(st, laplace)
    _ctsem_log_prior_gradient!(g, laplace.objective, st.theta)
    B = Matrix(Symmetric(B, :U))
    w = k == 1 ? 1.0 : Float64(info_rate)
    st.info .= (1 - w) .* st.info .+ w .* B
    prec = zeros(npar)
    obj = laplace.objective
    for j in eachindex(obj.prior_index)
        prec[obj.prior_index[j]] += obj.prior_weight / obj.prior_scale[j]^2
    end
    P = Symmetric(st.info) + Diagonal(prec .+ 1e-8 .* (1 .+ diag(st.info)))
    F = cholesky(P; check=false)
    step = issuccess(F) ? F \ g : g ./ (diag(P) .+ 1.0)
    step .*= gamma
    smax = maximum(abs, step; init=0.0)
    smax > maxstep && (step .*= maxstep / smax)
    st.theta .+= step
    return (logpost=logpost, gradient_norm=norm(g), step=min(smax, Float64(maxstep)))
end

export ctsem_saem_step!

"""
    _saem_plateaued(history, window)

Whether the complete-data log posterior has stopped rising: the mean of the
last quarter of the burn-in so far exceeds the mean of the quarter before it
by less than their combined standard error. Never before four windows.

Quarters rather than fixed windows, because the test has to see a slow climb:
two windows of 50 against each other passed a climb of 0.36 an iteration on
the SNSF pilot at iteration 800 as a plateau (the noise in the series is tens
of nats), and the burn-in ended a thousand iterations early. Over quarters the
rise it can miss falls like `q^-1.5`. The series is autocorrelated, so the
standard error is an underestimate, which errs toward a longer burn-in.
"""
function _saem_plateaued(history::Vector{Float64}, window::Integer)
    n = length(history)
    n < 4 * window && return false
    q = n ÷ 4
    a = @view history[(n - 2 * q + 1):(n - q)]
    b = @view history[(n - q + 1):n]
    ma = sum(a) / q; mb = sum(b) / q
    va = sum(abs2, a .- ma) / (q - 1); vb = sum(abs2, b .- mb) / (q - 1)
    return mb - ma < sqrt((va + vb) / q)
end

"""
    ctsem_saem(laplace, start; maxiter, burnin_max, averaging, window, alpha,
        seed, ...)

SAEM from `start`: step size one, with the expansion step (`_saem_expand!`,
off by default) if `expand`, until the complete-data log posterior plateaus
(`_saem_plateaued`) or `burnin_max` iterations, then `averaging` iterations
at `gamma = j^-alpha` with Polyak averaging, all within `maxiter`.
Returns the averaged point (`minimizer`), the iterations, where the burn-in
ended and whether it ended on a plateau rather than its cap (`settled`), the
mean acceptance rate and a trace. `progress`, `progress_*` and
`callback` behave as on `ctsem_optimize`; the callback receives the iteration,
`maxiter`, the complete-data log posterior, the score norm and the current
estimate (the average once averaging has begun).
"""
function ctsem_saem(laplace::CTSEMLaplaceObjective, start::AbstractVector;
    maxiter::Integer=3000, burnin_max::Integer=-1, averaging::Integer=-1,
    window::Integer=50, alpha::Real=0.7, seed::Integer=1, sweeps::Integer=2,
    nupper::Integer=2, maxstep::Real=0.25, refresh::Integer=25,
    info_rate::Real=0.1, expand::Bool=false, expand_min_groups::Integer=0,
    progress::Bool=false,
    progress_overwrite::Bool=true, progress_sink=nothing, progress_every::Real=0.0,
    callback=nothing)
    maxiter >= 1 || throw(ArgumentError("SAEM needs at least one iteration"))
    burnin_cap = burnin_max < 0 ? max(2 * window, (3 * Int(maxiter)) ÷ 4) : Int(burnin_max)
    st = ctsem_saem_init(laplace, start; seed=seed)
    reporter = _ctsem_progress_reporter(progress, "saem", progress_overwrite,
        progress_sink, progress_every)
    watcher = CTSEMCallback(callback)
    trace = CTSEMTrace(:logpost_complete, :gradient_norm, :gamma, :step, :acceptance)
    history = Float64[]
    burnin = 0
    averaged = 0
    target = 0
    settled = false
    for k in 1:Int(maxiter)
        averaging_now = burnin > 0
        gamma = averaging_now ? (k - burnin)^(-Float64(alpha)) : 1.0
        out = ctsem_saem_step!(st, laplace; gamma=gamma, sweeps=sweeps,
            nupper=nupper, maxstep=maxstep, refresh=refresh, info_rate=info_rate,
            expand=expand && !averaging_now, expand_min_groups=expand_min_groups)
        push!(history, out.logpost)
        if averaging_now
            st.nbar += 1
            st.thetabar .+= (st.theta .- st.thetabar) ./ st.nbar
            averaged += 1
        end
        acceptance = _saem_acceptance(st)
        _record!(trace, k, out.logpost, out.gradient_norm, gamma, out.step, acceptance)
        estimate = averaging_now ? st.thetabar : st.theta
        if _due(reporter)
            _progress_optimise(reporter, k, Int(maxiter),
                averaging_now ? @sprintf("averaging %d/%d", averaged, target) : "burn-in",
                @sprintf("logpost (complete) %11.2f", out.logpost),
                @sprintf("step %8.2e", out.step),
                @sprintf("accept %.2f", acceptance))
        end
        _invoke_callback(watcher, k, Int(maxiter), out.logpost, out.gradient_norm,
            estimate)
        if !averaging_now && (k >= burnin_cap ||
                (k % window == 0 && _saem_plateaued(history, window)))
            settled = k < burnin_cap
            burnin = k
            target = averaging < 0 ? clamp(k ÷ 2, 2 * window, 1000) : Int(averaging)
            target = min(target, Int(maxiter) - k)
            target <= 0 && break
        end
        averaging_now && averaged >= target && break
    end
    minimizer = st.nbar > 0 ? copy(st.thetabar) : copy(st.theta)
    _progress_done(reporter, @sprintf("%d iterations", st.iteration),
        burnin > 0 ? @sprintf("burn-in %d, averaged %d", burnin, averaged) :
            "burn-in not finished",
        @sprintf("accept %.2f", _saem_acceptance(st)))
    return (minimizer=minimizer, iterations=st.iteration, burnin=burnin,
        settled=settled, averaged=averaged, acceptance=_saem_acceptance(st),
        trace=_trace_result(trace), state=st)
end

export ctsem_saem

# The phase `ctsem_optimize` runs before L-BFGS when asked: SAEM on a Laplace
# objective, nothing on any other route (the R side refuses those first).
_ctsem_saem_phase(objective, start; kwargs...) = nothing
_ctsem_saem_phase(laplace::CTSEMLaplaceObjective, start; kwargs...) =
    ctsem_saem(laplace, start; kwargs...)
