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
# With `proposal = :laplace` each leaf first takes an independence move from
# its conditional Laplace approximation -- its mode and curvature (clipped at
# the prior's) given its ancestors' current values, from the quadrature's leaf
# rule (`_quadrature_leaf_rule!`) -- before its random walk: the f-SAEM kernel
# (Karimi, Lavielle and Moulines 2020), which jumps straight to where the
# current parameters put the subject instead of walking there. Its acceptance
# is kept apart from the random walk's, which alone adapts the walk's scale.
#
# Chains: a model with few units has few independent draws behind each
# population parameter, so it is replicated, as Monolix replicates a small
# data set, until units times chains reaches fifty. They run in parallel with
# the units, which also fills the workers a few-unit model leaves idle.
#
# Every chain starts at a draw from its unit's Laplace approximation rather
# than at the mode, so the first iterations' draws are not under-dispersed
# (`_saem_disperse!`).
#
# The parameters: two steps an iteration, each a full one, under two
# augmentations of the same data (AECM; `_saem_centre!` explains why both).
#   * Holding the standardised effects `u` fixed: `theta += P \ g`, `g` the
#     complete-data score at the draws averaged over chains, plus the prior's
#     gradient -- the fixed-u sweep `_laplace_floored_unit_gradient!` takes --
#     and `P` the outer product of the per-member scores plus the prior
#     precision, each coordinate capped at `maxstep`. That moves every
#     parameter, and is the only step for those without random effects.
#   * Holding the effects themselves fixed: each level's population mean and
#     scale are a Gaussian's, given the effects, and step to that Gaussian's
#     maximum, the draws re-expressed so no member's likelihood moves.
# `P` averages the outer products of earlier iterations, at rate `info_rate`,
# and never this iteration's own. With this iteration's, the step divides the
# score by a matrix built from the same draws -- a self-normalised estimate
# whose fixed point is not the score's zero, 0.6 standard errors off on the
# rank-one test fixture. A slow average lags as theta moves, and with the
# centred step moving theta fast early on, the lag collapsed a variance on the
# bench's config A1 at rate 0.1; 0.3 did not. The remaining offset, about a
# tenth of a standard error on the fixtures, is the same at every rate from
# 0.3 to 0.01, so it is the full step's, not the average's.
#
# No phases and no step-size schedule. The steps never shrink; the estimate is
# the average of the last half of the iterates (suffix averaging), which
# forgets the transient without being told where it ended and averages the
# Monte Carlo noise away after it. The run stops once the iterates have
# stopped travelling (`_saem_trend`): every parameter's change between the
# last two quarters of the run is within its own Monte Carlo noise, measured
# from batch means, or else under a tenth of a standard error, so a parameter
# creeping along a flat direction toward a boundary does not hold the run
# open. A full step leaves the stationary average a little off the exact mode
# where the score is nonlinear in the parameters; the optimiser's finish and
# certification (`ctsem_optimize`), which follow, remove that and are the one
# serious check of convergence.
#
# What was built and measured before this, and removed. A burn-in ended by a
# plateau test fired a thousand iterations early on the SNSF pilot. Kesten's
# rule (a parameter's step shrinking with its sign changes) shrank every step
# from the first iterations, because Monte Carlo noise flips the signs as often
# as crossing the optimum does, and froze 12 nats short on config A1. Constant
# steps in the fixed-u augmentation alone wandered there for 10000 iterations:
# 98 per cent of the information about the CINT mean is missing in it, so it
# moves at EM's pace. Louis' identity for the marginal curvature failed for the
# same reason -- the marginal information is then the difference of two
# matrices a hundred times larger, and its Monte Carlo error was as large as
# itself.
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
# normal), the information's averaging rate (0.3), the shape refresh interval
# (25), two sweeps and two collapsed moves
# per block per iteration, fifty units' worth of chains, ten batches of at
# least five iterates per quarter for the trend, and the tenth of a standard
# error below which a change does not count.

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
    # Per unit, per chain: the independence move's counts (`proposal = :laplace`).
    indep_accepted::Vector{Vector{Int}}
    indep_proposed::Vector{Vector{Int}}
    # Per unit: the blocks with no descendants, the rest innermost first, and
    # each block's descendants.
    leaves::Vector{Vector{Int}}
    uppers::Vector{Vector{Int}}
    descendants::Vector{Vector{Vector{Int}}}
    # Per unit, per chain: the members' complete-data scores at the last M-step.
    scores::Vector{Vector{Matrix{Float64}}}
    # The complete-data information, averaged over iterations: the step's
    # metric (taken before this iteration's is added) and the trend's
    # yardstick for what counts as a negligible change.
    info::Matrix{Float64}
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
    _saem_independence_move!(st, laplace, U, c, b, Ls)

Leaf `b` of chain `c` of unit `U`: one Metropolis-Hastings move proposing from
the leaf's conditional Laplace approximation given its ancestors' current
values -- the mode and the prior-clipped curvature `_quadrature_leaf_rule!`
places, Newton from the current draw. Leaves the draw where it was when the
rule cannot be placed.
"""
function _saem_independence_move!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, c::Integer, b::Integer, Ls::Vector{Matrix{Float64}})
    blocks = laplace.units.blocks[U]
    blk = blocks[b]
    u = st.u[U][c]; ll = st.ll[U][c]; theta = st.theta
    cols = (blk.offset + 1):(blk.offset + blk.size)
    aws = _laplace_workspace!(laplace, Float64, length(theta))
    rule = _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws; start=u[cols])
    rule.ok || return false
    R = UpperTriangular(rule.scale)
    logq(x) = -sum(abs2, R \ (x .- rule.centre)) / 2
    rng = _saem_rng(st, U, c, b, 0, 3)
    old = u[cols]
    z = rule.centre .+ rule.scale * randn(rng, blk.size)
    u[cols] .= z
    newll = [_saem_member_ll(laplace, U, m, theta, Ls, u) for m in blk.members]
    logr = sum(newll) - sum(@view ll[blk.members]) -
        (sum(abs2, z) - sum(abs2, old)) / 2 - (logq(z) - logq(old))
    accepted = isfinite(logr) && log(rand(rng)) < logr
    accepted ? (ll[blk.members] .= newll) : (u[cols] .= old)
    # Counted per unit and chain, summed over its leaves: leaves of one chain
    # run in parallel, so the increment is atomic.
    _saem_count!(st.indep_proposed[U], c, 1)
    accepted && _saem_count!(st.indep_accepted[U], c, 1)
    return accepted
end

const _SAEM_COUNT_LOCK = ReentrantLock()
_saem_count!(v::Vector{Int}, c::Integer, n::Integer) =
    lock(() -> (v[c] += n), _SAEM_COUNT_LOCK)

"""
    _saem_sweep!(st, laplace, U, c, Ls, sweep, adapt; nupper)

One sweep of chain `c` of unit `U`: every leaf once, in parallel, then
`nupper` collapsed moves of every block above, innermost first.
"""
function _saem_sweep!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, c::Integer, Ls::Vector{Matrix{Float64}}, sweep::Integer,
    adapt::Float64; nupper::Integer=2, independence::Bool=false)
    blocks = laplace.units.blocks[U]
    u = st.u[U][c]; ll = st.ll[U][c]; theta = st.theta
    independence && _laplace_parallel(laplace, st.leaves[U]) do b
        _saem_independence_move!(st, laplace, U, c, b, Ls)
        true
    end
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

"""
    _saem_centre!(st, laplace, prec)

The centred half of an iteration (AECM: Meng and van Dyk 1997; the
interweaving of Yu and Meng 2011). The M-step before it moves theta holding
the standardised effects `u` fixed, and in that augmentation a population
parameter whose effects the data determine well carries almost none of its
information: shifting a mean at fixed `u` shifts every subject, which the data
refuse, so the step is tiny and EM crawls (on the bench's config A1, 98 per
cent of the information about the CINT mean is missing). Holding the effects
themselves fixed instead, `b = theta + L u`, the same parameter is a plain
Gaussian mean and scale of the `b`, and its step is the full one. Each
augmentation is fast exactly where the other is slow, so alternating them
moves every population parameter quickly whatever the data say about the
subjects.

Per level, over its groups and every chain (a chain's draws weigh 1/chains):

* the mean: `theta[re] += L delta` and every `u -= delta`, `delta` maximising
  `sum log N(u - delta; 0, I) + log prior(theta + L delta)` -- in closed form;
  along `L`'s columns, which is all a reduced-rank level allows;
* the scale, full rank: the level's scales and correlations by Fisher scoring
  on `sum log N(L u; 0, L(phi) L(phi)') + log prior(phi)` with the deviations
  `L u` held fixed, in raw coordinates through `_laplace_popchol` and its
  derivatives, so whatever the covariance transform, and the effects
  re-expressed as `L(phi)^-1 L u`;
* the loadings, reduced rank: `L <- L A`, which keeps the column space, and
  `u <- A^-1 u`, with `A = chol(S)` for `S` the mean of `u u'`, or with the
  loadings' prior when they carry one (`_saem_centre_loadings!`).

Each moves only part of the way: by `I - Lambda` of its full step, `Lambda` the
share of each direction the data leave to the prior (see the code). Every
member's shifted parameters, so its likelihood, are left where they were; the
proposal shapes are carried into the new coordinates.
"""
function _saem_centre!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    prec::Vector{Float64})
    spec = laplace.spec
    K = st.chains
    for (l, level) in enumerate(spec.levels)
        k = nrandomeffects(level)
        r = nlatent(level)
        (k == 0 || r == 0) && continue
        sites = Tuple{Int,Int}[]
        for U in eachindex(st.u), (b, blk) in enumerate(laplace.units.blocks[U])
            blk.level == l && push!(sites, (U, b))
        end
        G = length(sites)
        G == 0 && continue
        cols(U, b) = (laplace.units.blocks[U][b].offset + 1):(laplace.units.blocks[U][b].offset +
            laplace.units.blocks[U][b].size)
        # How far to centre. Lambda, the mean over the level's groups of the
        # posterior covariance of `u` (the proposal shapes' prior-clipped
        # inverse curvature), is the share of each direction the data leave to
        # the prior. The centred augmentation is the efficient one where the
        # data pin the effects down (Lambda near 0) and the fixed-u one where
        # they barely inform them (Lambda near I): there the centred step reads
        # a population scale off draws that are mostly prior, follows their
        # noise with almost nothing pulling it back, and walks -- on the
        # bench's bigre a drift effect's scale walked from -1.0 to -2.4 raw.
        # So the centred step moves each direction by `I - Lambda` of its full
        # step, the weight with which a partially non-centred
        # parameterisation (Papaspiliopoulos, Roberts and Skold 2003) leaves a
        # Gaussian mean no missing information; fixed points are unchanged.
        Lambda = zeros(r, r)
        for (U, b) in sites
            C = st.chol[U][b]
            Lambda .+= C * transpose(C)
        end
        Lambda ./= G
        E = eigen(Symmetric(Lambda))
        keep = clamp.(1 .- E.values, 0.0, 1.0)
        W = E.vectors * Diagonal(keep) * transpose(E.vectors)
        Whalf = E.vectors * Diagonal(sqrt.(keep)) * transpose(E.vectors)
        # The mean.
        L = _laplace_popchol(st.theta, level)
        usum = zeros(r)
        for (U, b) in sites, c in 1:K
            usum .+= @view st.u[U][c][cols(U, b)]
        end
        usum ./= K
        Dm = prec[level.re_index]
        mu = st.theta[level.re_index]
        A = Symmetric(G .* Matrix{Float64}(I, r, r) .+ transpose(L) * (Dm .* L))
        delta = W * (A \ (usum .- transpose(L) * (Dm .* mu)))
        st.theta[level.re_index] .+= L * delta
        for (U, b) in sites, c in 1:K
            st.u[U][c][cols(U, b)] .-= delta
        end
        # The scale, from the effects' second moment with its departure from
        # the identity -- what the full step would act on -- shrunk the same
        # way: `G I + Whalf (S - G I) Whalf`.
        S = zeros(r, r)
        for (U, b) in sites, c in 1:K
            v = st.u[U][c][cols(U, b)]
            S .+= v * transpose(v)
        end
        S ./= K
        S = Matrix(Symmetric(G .* Matrix{Float64}(I, r, r) .+
            Whalf * (S .- G .* Matrix{Float64}(I, r, r)) * Whalf))
        T = isreducedrank(level) ? _saem_centre_loadings!(st, laplace, l, S, G, prec) :
            _saem_centre_scale!(st, laplace, l, S, G, prec)
        T === nothing && continue
        Tinv = inv(T)
        for (U, b) in sites
            for c in 1:K
                st.u[U][c][cols(U, b)] = T * st.u[U][c][cols(U, b)]
            end
            st.chol[U][b] = Matrix(T * st.chol[U][b])
        end
        # A block's response maps its ancestor's coordinates to its own.
        for U in eachindex(st.u)
            blocks = laplace.units.blocks[U]
            for (c, blk) in enumerate(blocks), t in eachindex(blk.ancestors)
                isempty(st.response[U][c]) && continue
                a = blk.ancestors[t]
                blk.level == l && (st.response[U][c][t] = T * st.response[U][c][t])
                blocks[a].level == l && (st.response[U][c][t] = st.response[U][c][t] * Tinv)
            end
        end
    end
    return st
end

"""
    _saem_centre_scale!(st, laplace, l, S, G, prec)

Full-rank level `l`'s scales and correlations by Fisher scoring on the
Gaussian log likelihood of `G` groups' deviations `d = L u`, held fixed, whose
summed second moment in the current coordinates is `S` (`sum u u'`, averaged
over chains), plus the prior; with a halving line search on that same
objective. With `M_p = L^-1 dL_p` and `Sw = L^-1 L0 S L0' L^-T`, the gradient
is `tr(M_p (Sw - G I))` and the information
`G/2 tr((M_p + M_p')(M_q + M_q'))`. Returns the map `L_new^-1 L_old` that
re-expresses the effects, or `nothing` where the level has no scale
parameters or nothing moved.
"""
function _saem_centre_scale!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, S::Matrix{Float64}, G::Integer, prec::Vector{Float64})
    spec = laplace.spec
    level = spec.levels[l]
    positions = _laplace_level_positions(spec, l)
    isempty(positions) && return nothing
    L0 = _laplace_popchol(st.theta, level)
    D2 = Matrix(Symmetric(L0 * S * transpose(L0)))
    Dp = prec[positions]
    function second_moment(L)
        X = L \ D2
        return Matrix(Symmetric(transpose(L \ transpose(X))))
    end
    function objective(theta)
        L = LowerTriangular(_laplace_popchol(theta, level))
        phi = theta[positions]
        return -G * sum(log, diag(L)) - tr(second_moment(L)) / 2 - sum(Dp .* phi .^ 2) / 2
    end
    theta = copy(st.theta)
    q = objective(theta)
    np = length(positions)
    moved = false
    for _ in 1:25
        L = LowerTriangular(_laplace_popchol(theta, level))
        dL = _laplace_level_chol_derivatives(theta, spec, l)
        Sw = second_moment(L)
        M = [Matrix(L \ dL[t]) for t in 1:np]
        grad = [tr(M[t] * (Sw - G * I)) for t in 1:np] .- Dp .* theta[positions]
        F = zeros(np, np)
        for t1 in 1:np, t2 in t1:np
            F[t1, t2] = F[t2, t1] = G / 2 * tr((M[t1] + transpose(M[t1])) *
                (M[t2] + transpose(M[t2])))
        end
        F .+= Diagonal(Dp .+ 1e-8 .* (1 .+ diag(F)))
        step = Symmetric(F) \ grad
        all(isfinite, step) || break
        maximum(abs, step) < 1e-8 && break
        accepted = false
        for _ in 1:30
            trial = copy(theta)
            trial[positions] .+= step
            qt = objective(trial)
            if isfinite(qt) && qt >= q
                theta = trial; q = qt; accepted = true; moved = true
                break
            end
            step ./= 2
        end
        accepted || break
    end
    moved || return nothing
    st.theta[positions] .= theta[positions]
    L1 = _laplace_popchol(st.theta, level)
    return Matrix(LowerTriangular(L1) \ L0)
end

"""
    _saem_centre_loadings!(st, laplace, l, S, G, prec)

Reduced-rank level `l`: `L <- L A`, `A` lower triangular, which keeps the
loadings' column space, and `u <- A^-1 u`, with `A` maximising the Gaussian
log likelihood of `G` groups' deviations, held fixed, whose summed second
moment in the current coordinates is `S`, plus the loadings' prior. Without a
prior that is `A A' = S / G`; with one, Newton on `A`'s entries from there
(`_saem_loading_factor`). Returns `A^-1`.
"""
function _saem_centre_loadings!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, S::Matrix{Float64}, G::Integer, prec::Vector{Float64})
    level = laplace.spec.levels[l]
    k = nrandomeffects(level)
    r = nlatent(level)
    G > r || return nothing
    F = cholesky(Symmetric(S ./ G); check=false)
    issuccess(F) || return nothing
    A = Matrix(F.L)
    R = zeros(k, r)
    D = zeros(k, r)
    counter = 0
    for q in 1:r, p in q:k
        counter += 1
        R[p, q] = st.theta[level.load_index[counter]]
        D[p, q] = prec[level.load_index[counter]]
    end
    if any(>(0), D)
        A = _saem_loading_factor(A, S, R, D, G)
        A === nothing && return nothing
    end
    R = R * A
    counter = 0
    for q in 1:r, p in q:k
        counter += 1
        st.theta[level.load_index[counter]] = R[p, q]
    end
    return Matrix(inv(LowerTriangular(A)))
end

"""
    _saem_loading_factor(A0, S, R, D, G)

The lower-triangular `A` maximising
`-G log|det A| - tr(A^-1 S A^-T) / 2 - sum(D .* (R A).^2) / 2`: the
log likelihood of `G` groups' deviations `L u` (summed second moment `S` of
the `u`) under the loadings `L A`, plus a Gaussian prior of precision `D` on
the raw loadings `R A`. Newton on `A`'s `r (r + 1) / 2` entries with a halving
line search, from `A0`. `nothing` where no step could be taken from the start.
"""
function _saem_loading_factor(A0::Matrix{Float64}, S::Matrix{Float64},
    R::Matrix{Float64}, D::Matrix{Float64}, G::Integer)
    r = size(A0, 1)
    idx = [(p, q) for q in 1:r for p in q:r]
    function unpack(x)
        A = zeros(eltype(x), r, r)
        for (j, (p, q)) in enumerate(idx)
            A[p, q] = x[j]
        end
        return A
    end
    function f(x)
        A = unpack(x)
        any(iszero, diag(A)) && return oftype(x[1], -Inf)
        Ainv = inv(LowerTriangular(A))
        return -G * sum(a -> log(abs(a)), diag(A)) - tr(Ainv * S * transpose(Ainv)) / 2 -
            sum(D .* (R * A) .^ 2) / 2
    end
    x = [A0[p, q] for (p, q) in idx]
    fx = f(x)
    isfinite(fx) || return nothing
    for _ in 1:25
        g = ForwardDiff.gradient(f, x)
        H = ForwardDiff.hessian(f, x)
        F = cholesky(Symmetric(-H); check=false)
        step = issuccess(F) ? F \ g : g ./ (abs.(diag(H)) .+ 1.0)
        slope = dot(g, step)
        (all(isfinite, step) && slope > 0) || break
        t = 1.0
        moved = false
        for _ in 1:30
            trial = x .+ t .* step
            ft = f(trial)
            if isfinite(ft) && ft >= fx + 1e-4 * t * slope
                x = trial; fx = ft; moved = true
                break
            end
            t /= 2
        end
        (moved && maximum(abs, t .* step) > 1e-10) || break
    end
    return Matrix{Float64}(unpack(x))
end

"""Fifty units' worth of chains: `cld(50, nunits)`, between one and eight."""
_saem_default_chains(nunits::Integer) = clamp(cld(50, max(nunits, 1)), 1, 8)

"""
    ctsem_saem_init(laplace, theta; seed, chains, disperse)

SAEM's state at `theta`: every chain of each unit starts at a draw from the
unit's Laplace approximation (`_saem_disperse!`; at the mode itself with
`disperse = false`), and the proposal shapes come from the curvature at the
mode. `chains = 0` takes `_saem_default_chains`.
"""
function ctsem_saem_init(laplace::CTSEMLaplaceObjective, theta::AbstractVector;
    seed::Integer=1, chains::Integer=0, disperse::Bool=true)
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
        [zeros(Int, K) for U in 1:nunits], [zeros(Int, K) for U in 1:nunits],
        leaves, uppers, descendants,
        [[zeros(npar, nmem(U)) for _ in 1:K] for U in 1:nunits],
        zeros(npar, npar), Vector{Float64}[], K, 0, UInt64(seed))
    ok = _laplace_parallel(laplace, 1:nunits) do U
        _saem_refresh!(st, laplace, U, Ls)
        disperse && _saem_disperse!(st, laplace, U, Ls)
        for c in 1:K
            _saem_members_ll!(st.ll[U][c], laplace, U, eachindex(laplace.units.members[U]),
                x, Ls, st.u[U][c])
            # A draw the likelihood refuses goes back to the mode.
            if disperse && !all(isfinite, st.ll[U][c])
                st.u[U][c] .= laplace.modes[U]
                _saem_members_ll!(st.ll[U][c], laplace, U,
                    eachindex(laplace.units.members[U]), x, Ls, st.u[U][c])
            end
        end
        all(c -> all(isfinite, st.ll[U][c]), 1:K)
    end
    ok || throw(ArgumentError("SAEM cannot start: a unit's log likelihood is not " *
        "finite at its Laplace mode"))
    return st
end

"""
    _saem_disperse!(st, laplace, U, Ls)

Unit `U`'s chains start at draws from its Laplace approximation, `N(mode,
M^-1)` with `M` its curvature at the mode clipped at the prior's, not at the
mode itself. A mode is shrunk toward zero -- across units its second moment
falls short of the identity by the posterior variance -- so draws started
there are under-dispersed for the first sweeps, and the centred step, which
reads the population scale off them, would take that shrinkage for the
population's.
"""
function _saem_disperse!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, Ls::Vector{Matrix{Float64}})
    blocks = laplace.units.blocks[U]
    isempty(blocks) && return st
    mode = laplace.modes[U]
    M = _laplace_unit_curvature(laplace, U, st.theta, Ls, mode)
    dense = _laplace_block_dense(M, blocks, length(mode))
    all(isfinite, dense) || return st
    F = _saem_proposal_factor(_laplace_symmetrise(Matrix{Float64}(dense)))
    for c in 1:st.chains
        rng = _saem_rng(st, U, c, 0, 0, 4)
        st.u[U][c] .= mode .+ F * randn(rng, length(mode))
    end
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
    _saem_curvature(B, prec)

A complete-data information `B` plus the prior precision, with a ridge
relative to its own diagonal.
"""
_saem_curvature(B::AbstractMatrix, prec::Vector{Float64}) =
    Symmetric(B) + Diagonal(prec .+ 1e-8 .* (1 .+ diag(B)))

"""
    ctsem_saem_step!(st, laplace; sweeps, nupper, maxstep, refresh, info_rate,
        adapt, prec, proposal, centre, mstep)

One SAEM iteration: `sweeps` E-step sweeps of every chain of every unit (units
and chains in parallel, and within one its leaves and its members), then the
two M-steps of the file's header -- the full step at fixed `u`, each
coordinate capped at `maxstep`, and with `centre` the centred step at fixed
effects. `mstep = false` leaves theta alone. Returns the complete-data log
posterior before the step, the score norm, the largest coordinate of the
first step, and the gradient and information it took.
"""
function ctsem_saem_step!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective;
    sweeps::Integer=2, nupper::Integer=2, maxstep::Real=0.25,
    refresh::Integer=25, info_rate::Real=0.3, adapt::Real=-1.0,
    prec::Vector{Float64}=_saem_prior_precision(laplace, length(st.theta)),
    proposal::Symbol=:rw, centre::Bool=true, mstep::Bool=true)
    proposal in (:rw, :laplace) ||
        throw(ArgumentError("SAEM proposal must be :rw or :laplace"))
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
            _saem_sweep!(st, laplace, U, c, Ls, sweep, rate; nupper=nupper,
                independence=(proposal === :laplace && sweep == 1))
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
    # The step's information comes from earlier draws than its score: this
    # iteration's outer product moves with this iteration's score, and
    # dividing one by the other is a self-normalised estimate whose fixed point
    # is not the score's zero.
    P = _saem_curvature(k == 1 ? B : st.info, prec)
    w = k == 1 ? 1.0 : Float64(info_rate)
    st.info .= (1 - w) .* st.info .+ w .* B
    step = zeros(npar)
    if mstep
        F = cholesky(P; check=false)
        d = issuccess(F) ? F \ g : g ./ (diag(P) .+ 1.0)
        step = clamp.(d, -maxstep, maxstep)
        # Capping coordinates one by one can turn an ascent direction into
        # one that is not; scaled whole, it stays one.
        dot(g, step) > 0 || (step = d .* min(1.0, maxstep / maximum(abs, d; init=0.0)))
        step = _saem_ascend!(st, laplace, items, step, g)
        centre && _saem_centre!(st, laplace, prec)
    end
    push!(st.csum, isempty(st.csum) ? copy(st.theta) : st.csum[end] .+ st.theta)
    return (logpost=logpost, gradient_norm=norm(g),
        step=maximum(abs, step; init=0.0), gradient=g, information=B)
end

export ctsem_saem_step!

"""
    _saem_ascend!(st, laplace, items, step, g; halvings)

Take `step` from theta, or the largest of its halvings that raises the
complete-data log posterior at the current draws by at least a sliver of its
predicted gain (`g' step`, Armijo), or none: generalised EM's condition on an
M-step. The comparison is at fixed draws, so it carries no Monte Carlo noise;
it costs one likelihood pass over every member per trial, and the first trial
is usually taken. Without it the step's information, an average over earlier
iterations, could point the step anywhere while the centred step moved theta
fast -- on the bench's cf_gaussian it drove drift and diffusion the full cap
every iteration into a region where both transforms are flat, and they froze
there.
Returns the step taken.
"""
function _saem_ascend!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective, items,
    step::Vector{Float64}, g::Vector{Float64}; halvings::Integer=10)
    K = st.chains
    objective = laplace.objective
    q0 = sum(sum(st.ll[U][c]; init=0.0) for (U, c) in items; init=0.0) / K +
        _ctsem_log_prior(objective, st.theta)
    slope = dot(g, step)
    (isfinite(q0) && slope > 0) || return zero(step)
    trial = [[similar(st.ll[U][c]) for c in 1:K] for U in eachindex(st.ll)]
    theta0 = copy(st.theta)
    alpha = 1.0
    for _ in 0:halvings
        st.theta .= theta0 .+ alpha .* step
        Ls = _laplace_popchols(st.theta, laplace.spec)
        _laplace_parallel(laplace, items) do item
            U, c = item
            _saem_members_ll!(trial[U][c], laplace, U, eachindex(laplace.units.members[U]),
                st.theta, Ls, st.u[U][c])
            true
        end
        q = sum(sum(trial[U][c]; init=0.0) for (U, c) in items; init=0.0) / K +
            _ctsem_log_prior(objective, st.theta)
        isfinite(q) && q >= q0 + 1e-4 * alpha * slope && return alpha .* step
        alpha /= 2
    end
    st.theta .= theta0
    return zero(step)
end

"""Average of the iterates `a+1:b` from the cumulative sums."""
_saem_window(st::CTSEMSAEMState, a::Integer, b::Integer) =
    (a == 0 ? st.csum[b] : st.csum[b] .- st.csum[a]) ./ (b - a)

"""
    _saem_trend(st, P; batches, negligible)

Whether the iterates are still travelling. For each parameter, the change
`Delta` between the averages of the last two quarters of the run is divided by
its own Monte Carlo standard error, the spread of `batches` batch means within
each quarter -- so iterates correlated from one iteration to the next are not
taken for independent ones (with batches shorter than that correlation the
error comes out small and the rule errs toward running on), and no
information matrix says what "small" means. A parameter whose change is under
a tenth of a standard error in the complete-data information `P`
(`Delta^2 P_jj < negligible`), which overstates the marginal information and so
errs strict, counts as still: one creeping along a flat direction toward a
boundary moves the fit by nothing. Returns the mean of the squared ratios over
the parameters, near `_saem_trend_null()` once the two quarters are draws of
one stationary distribution and larger while the run still has somewhere to
go; `NaN` until each batch holds five iterates.
"""
function _saem_trend(st::CTSEMSAEMState, P::AbstractMatrix; batches::Integer=10,
    negligible::Real=0.01)
    k = length(st.csum)
    q = k ÷ 4
    bl = q ÷ batches
    bl < 5 && return NaN
    npar = length(st.theta)
    function quarter(stop)
        means = [_saem_window(st, stop - (batches - j + 1) * bl, stop - (batches - j) * bl)
                 for j in 1:batches]
        m = sum(means) ./ batches
        v = sum(x -> (x .- m) .^ 2, means) ./ ((batches - 1) * batches)
        return m, v
    end
    mA, vA = quarter(k - q)
    mB, vB = quarter(k)
    total = 0.0
    for j in 1:npar
        delta = mB[j] - mA[j]
        v = vA[j] + vB[j]
        (v > 0 && delta^2 * P[j, j] >= negligible) || continue
        total += delta^2 / v
    end
    return total / npar
end

"""The trend's expectation for a stationary run in which every parameter
counts: a squared difference over its batch-means variance, on
`2 (batches - 1)` degrees of freedom."""
_saem_trend_null(batches::Integer=10) = (2 * (batches - 1)) / (2 * (batches - 1) - 2)

"""
    ctsem_saem(laplace, start; maxiter, seed, chains, centre, proposal, ...)

SAEM from `start` until the iterates have stopped travelling (`_saem_trend` at
or below its stationary expectation, `_saem_trend_null`), or `maxiter`
iterations. Every step is a full one; there are no phases and no schedule (see
the file's header). Returns the estimate (`minimizer`, the average of the last
half of the iterates), the iterations, whether it stopped on the trend rule
(`settled`) and the trend there, the chains, the mean acceptance rate and a
trace. `centre = false` drops the centred step. `progress`, `progress_*` and
`callback` behave as on `ctsem_optimize`; the callback receives the
iteration, `maxiter`, the complete-data log posterior, the score norm and the
current estimate.
"""
function ctsem_saem(laplace::CTSEMLaplaceObjective, start::AbstractVector;
    maxiter::Integer=10000, seed::Integer=1, chains::Integer=0,
    sweeps::Integer=2, nupper::Integer=2, maxstep::Real=0.25, refresh::Integer=25,
    info_rate::Real=0.3, proposal=:rw, centre::Bool=true, progress::Bool=false,
    progress_overwrite::Bool=true, progress_sink=nothing, progress_every::Real=0.0,
    callback=nothing)
    proposal = Symbol(proposal)
    maxiter >= 1 || throw(ArgumentError("SAEM needs at least one iteration"))
    st = ctsem_saem_init(laplace, start; seed=seed, chains=chains)
    prec = _saem_prior_precision(laplace, length(st.theta))
    reporter = _ctsem_progress_reporter(progress, "saem", progress_overwrite,
        progress_sink, progress_every)
    watcher = CTSEMCallback(callback)
    trace = CTSEMTrace(:logpost_complete, :gradient_norm, :step, :acceptance,
        :trend)
    settled = false
    trend = NaN
    null = _saem_trend_null()
    for k in 1:Int(maxiter)
        out = ctsem_saem_step!(st, laplace; sweeps=sweeps, nupper=nupper,
            maxstep=maxstep, refresh=refresh, info_rate=info_rate, prec=prec,
            proposal=proposal, centre=centre)
        trend = _saem_trend(st, _saem_curvature(st.info, prec))
        acceptance = _saem_acceptance(st)
        _record!(trace, k, out.logpost, out.gradient_norm, out.step, acceptance,
            trend)
        estimate = _saem_window(st, k - max(1, k ÷ 2), k)
        if _due(reporter)
            _progress_optimise(reporter, k, Int(maxiter),
                @sprintf("logpost (complete) %11.2f", out.logpost),
                @sprintf("trend %6.2f", trend),
                @sprintf("accept %.2f", acceptance))
        end
        _invoke_callback(watcher, k, Int(maxiter), out.logpost, out.gradient_norm,
            estimate)
        if isfinite(trend) && trend <= null
            settled = true
            break
        end
    end
    k = st.iteration
    minimizer = _saem_window(st, k - max(1, k ÷ 2), k)
    _progress_done(reporter, @sprintf("%d iterations", k),
        settled ? @sprintf("settled, trend %.2f", trend) :
            @sprintf("not settled, trend %.2f", trend),
        @sprintf("%d chain%s, accept %.2f", st.chains, st.chains == 1 ? "" : "s",
            _saem_acceptance(st)))
    indep = sum(sum, st.indep_proposed; init=0)
    return (minimizer=minimizer, iterations=k, settled=settled, trend=trend,
        chains=st.chains, acceptance=_saem_acceptance(st),
        independence_acceptance=indep == 0 ? NaN :
            sum(sum, st.indep_accepted; init=0) / indep,
        trace=_trace_result(trace), state=st)
end

export ctsem_saem

# The phase `ctsem_optimize` runs before L-BFGS when asked: SAEM on a Laplace
# objective, nothing on any other route (the R side refuses those first).
_ctsem_saem_phase(objective, start; kwargs...) = nothing
_ctsem_saem_phase(laplace::CTSEMLaplaceObjective, start; kwargs...) =
    ctsem_saem(laplace, start; kwargs...)
