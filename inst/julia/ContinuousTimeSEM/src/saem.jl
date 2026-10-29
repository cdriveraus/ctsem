################################################################################
# SAEM: the exact marginal posterior mode, by stochastic approximation
################################################################################
#
# Two algorithms share this file and its state. `ctsem_saem` is SAEM, an
# optimiser: a phase `ctsem_optimize` runs first when asked
# (`optimcontrol$saem`), and what places the joint samplers' chains
# (`_saem_placement`). `ctsem_saem_sample`, in the second half, is not SAEM:
# it is a data-augmentation sampler that draws the effects with SAEM's E-step
# and draws theta where SAEM steps it (its own header, below).
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
# Where this departs from textbook SAEM, and what each departure costs.
#
# The M-step steps theta rather than solving for it. Textbook SAEM (Delyon,
# Lavielle and Moulines 1999) keeps a stochastic average of the complete-data
# sufficient statistics and maximises in closed form, which needs an
# exponential-family complete-data model. Drift and diffusion reach the
# likelihood through the filter's matrix exponentials and are not that, so
# theta itself takes a scoring step on the complete-data score (above), and
# the step cap, the lagged information and the ascent condition
# (`_saem_ascend!`) come with it. Only the level means and scales, a
# Gaussian's given the effects, keep a closed form, in the centred half.
#
# No step-size schedule. Textbook SAEM takes full steps through a burn-in and
# then decreasing ones (gamma_k -> 0, sum gamma_k infinite, sum gamma_k^2
# finite), and the decrease is what makes its own iterate converge. Here the
# steps never shrink, and the estimate is the average of the last half of the
# iterates (suffix averaging), which forgets the transient without being told
# where it ended and averages the Monte Carlo noise away after it. The
# decrease has to begin after the transient -- begun early, it leaves the
# iterate unable to travel the rest of the way -- and both automatic switches
# tried began it early. A burn-in ended by a plateau test fired a thousand
# iterations early on the SNSF pilot. Kesten's rule (a parameter's step
# shrinking with its sign changes) shrank every step from the first
# iterations, because Monte Carlo noise flips the signs as often as crossing
# the optimum does, and froze 12 nats short on config A1. What constant steps
# cost: the average loses the noise but not a full step's offset -- the
# stationary average sits a little off the exact mode where the score is
# nonlinear in the parameters, about a tenth of a standard error on the
# fixtures -- and the stop below bounds the iterates' travel, not the
# average's Monte Carlo error, which along a weakly identified direction can
# be large enough to decide between maxima
# (review/SAEM-fell-short-2026-10-05.md). The phased version this replaced,
# with a decaying step, also got further in fewer iterations on the SNSF
# pilot. So the estimate is a start for the finish, not a final point, and a
# second phase begun where the stop below fires did not make it one, in any of
# three forms measured on the SAEM-fell-short model
# (review/SAEM-smooth-2026-10-13.md; branch saem-smooth). A decreasing gain on
# this step froze the directions where information is missing, leaving the
# estimate 1.8 nats short where the first phase had. Stochastic Newton on the
# Laplace Hessian moved them, but multiplied the fixed-u score's noise by the
# inverse marginal information and ended 12-71 nats short. Accelerating by the
# rate read off the iterates' lag-one autocorrelation ended 48-137 short: the
# E-step chain's own persistence dominates that autocorrelation, so nearly
# every direction looked slow. Moving those directions quickly means moving
# them on the one augmentation that informs them, whose noise is the missing
# information itself.
#
# The stop. The run stops once the iterates have stopped travelling
# (`_saem_trend`): every parameter's change between the last two quarters of
# the run is within its own Monte Carlo noise, measured from batch means, or
# else under a tenth of a standard error, so a parameter creeping along a flat
# direction toward a boundary does not hold the run open.
#
# The finish is on the Laplace objective. Nothing evaluates the exact
# marginal posterior deterministically for every model -- estimating it is
# SAEM's job -- and the optimiser's finish (Newton on the exact Hessian), its
# certification and the Hessian's standard errors each need an objective that
# is deterministic and differentiable. The Laplace objective with its gated
# floor is the one every model has; from SAEM's point the finish removes the
# step's offset and noise, and it is the one serious check of convergence.
# What it costs: where the Laplace mode and the exact one differ, the finish
# returns to the Laplace mode. The quadrature correction after it
# (`laplace_refit`) continues the fit with quadrature in place of Laplace
# for each unit of at most five effects whose values differ; a wider unit -- a
# nested study with its people -- keeps the Laplace term, and there SAEM's own
# estimate, better on the exact likelihood than the Laplace optimum on the
# bench's gN1 and gN3, is given back (review/SAEM-aecm-2026-10-02.md). So in
# the optimiser SAEM chooses which maximum the finish climbs, and ctFit checks
# that choice against a finish from the fit's own start, on the
# quadrature-corrected objective. An estimate free of the Laplace error is the
# sampler's to give.
#
# Two more measured failures. Constant steps in the fixed-u augmentation alone
# wandered on config A1 for 10000 iterations: 98 per cent of the information
# about the CINT mean is missing in it, so it moves at EM's pace. Louis'
# identity for the marginal curvature failed for the same reason -- the
# marginal information is then the difference of two matrices a hundred times
# larger, and its Monte Carlo error was as large as itself.
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
        _laplace_member_loglik(aws, so, shifted)
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

`theta` evaluates somewhere other than the chain's own parameters, which the
sampler's trajectories do; `values`, when given, receives each member's log
likelihood, which the same sweep computes anyway.
"""
function _saem_unit_scores!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    U::Integer, c::Integer, Ls::Vector{Matrix{Float64}}, dL, positions;
    theta::Vector{Float64}=st.theta, values=nothing)
    members = laplace.units.members[U]
    npar = length(theta)
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
            shifted = _laplace_member_values!(shift, theta, laplace.spec, Ls, u, offsets)
            # The complete-data score integrates the filter's own deviations,
            # as `_laplace_member_loglik` does.
            value = _laplace_subject_value_gradient!(grad,
                laplace.objective.subject_objectives[members[m]], aws, shifted,
                nothing, _HYBRID_FULL)
            if !(isfinite(value) && all(isfinite, grad))
                ok[] = false
                return nothing
            end
            values === nothing || (values[m] = Float64(value))
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
  along `L`'s columns, which is all a reduced-rank level allows. With it, the
  TI coefficients on the level's parameters (`_saem_level_ti`), as the
  regression they are given the effects: `u_g -= N_g eta` for a change `eta`
  (`_saem_mean_system`). In the fixed-u augmentation a coefficient on a
  well-determined parameter moves every subject the data pin down, as the
  mean does;
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
        # The mean, and the TI coefficients on the level's parameters with it:
        # each coefficient's direction in `u` is scaled by the share the data
        # take there, averaged over the groups by the size of their shift.
        ti = _saem_level_ti(laplace, l)
        system = _saem_mean_system(st, laplace, l, sites, prec, K, ti)
        system === nothing && (ti = _saem_level_ti_none;
            system = _saem_mean_system(st, laplace, l, sites, prec, K, ti))
        A, rhs, L, N, Z = system
        full = isempty(rhs) ? rhs : Symmetric(A) \ rhs
        rz = Z === nothing ? r : size(Z, 2)
        # Confined (`_saem_mean_confine`), the shrinkage acts within the
        # subspace the move may take.
        delta = Z === nothing ? W * full[1:r] :
            Z * ((transpose(Z) * W * Z) * full[1:rz])
        eta = full[(rz + 1):end]
        for c in eachindex(eta)
            num = 0.0
            den = 0.0
            for Ng in N
                v = @view Ng[:, c]
                num += dot(v, Lambda * v)
                den += dot(v, v)
            end
            eta[c] *= den > 0 ? clamp(1 - num / den, 0.0, 1.0) : 0.0
        end
        _saem_mean_apply!(st, laplace, l, sites, K, ti, L, N, delta, eta)
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
    # Any square root of the covariance (`_laplace_popchol`), so the log
    # determinant is taken whole rather than off a triangular diagonal.
    function objective(theta)
        L = _laplace_popchol(theta, level)
        phi = theta[positions]
        return -G * first(logabsdet(L)) - tr(second_moment(L)) / 2 - sum(Dp .* phi .^ 2) / 2
    end
    theta = copy(st.theta)
    q = objective(theta)
    np = length(positions)
    moved = false
    for _ in 1:25
        L = _laplace_popchol(theta, level)
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
    return Matrix(L1 \ L0)
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
    # The move rotates the triangular loadings; own-variance columns are not
    # part of that triangle, so a level with any keeps its coordinates.
    nunique(level) > 0 && return nothing
    k = nrandomeffects(level)
    r = nlatent(level)
    G > r || return nothing
    # `A` keeps to the entries the level's pattern allows
    # (`_laplace_loading_allowed`): the whole lower triangle for the triangle,
    # and each block's for a pattern of blocks, where the optimum without a
    # prior is each block's own `chol(S_bb / G)`.
    allowed = _laplace_loading_allowed(level)
    A, whole = _saem_block_cholesky(S ./ G, allowed)
    A === nothing && return nothing
    rows, cols = _laplace_loading_entries(level)
    R = zeros(k, r)
    D = zeros(k, r)
    for t in eachindex(rows)
        R[rows[t], cols[t]] = st.theta[level.load_index[t]]
        D[rows[t], cols[t]] = prec[level.load_index[t]]
    end
    # Newton where a prior shapes the optimum, or where some block's triangle
    # is not wholly allowed and its Cholesky factor would leave the pattern;
    # it moves the allowed entries alone.
    if any(>(0), D) || !whole
        A = _saem_loading_factor(A, S, R, D, G, allowed)
        A === nothing && return nothing
    end
    R = R * A
    for t in eachindex(rows)
        st.theta[level.load_index[t]] = R[rows[t], cols[t]]
    end
    return Matrix(inv(LowerTriangular(A)))
end

"""
    _saem_block_cholesky(M, allowed)

The lower Cholesky factor of each block of `M` that `allowed` connects
(`_laplace_loading_allowed`), zero between blocks, and whether every pair
within each block is allowed, so that the factors stay on the pattern.
`(nothing, false)` where a block is not positive definite. One block, the
whole of `M`, for the triangle.
"""
function _saem_block_cholesky(M::AbstractMatrix, allowed::AbstractMatrix{Bool})
    r = size(M, 1)
    block = collect(1:r)
    root(i) = (while block[i] != i; i = block[i]; end; i)
    for q in 1:r, m in (q + 1):r
        allowed[m, q] || continue
        a, b = root(m), root(q)
        a == b || (block[max(a, b)] = min(a, b))
    end
    A = zeros(r, r)
    whole = true
    for b in unique(root.(1:r))
        members = [i for i in 1:r if root(i) == b]
        F = cholesky(Symmetric(M[members, members]); check=false)
        issuccess(F) || return (nothing, false)
        A[members, members] .= F.L
        for j in eachindex(members), m in members[(j + 1):end]
            allowed[m, members[j]] || (whole = false)
        end
    end
    return (A, whole)
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
    R::Matrix{Float64}, D::Matrix{Float64}, G::Integer,
    allowed::AbstractMatrix{Bool}=trues(size(A0)))
    r = size(A0, 1)
    # Only the entries the level's pattern allows (`_laplace_loading_allowed`).
    idx = [(p, q) for q in 1:r for p in q:r if allowed[p, q]]
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
        g = _ctsem_gradient(f, x)
        H = _ctsem_dual_hessian(f, x)
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
    # The state carries any stated population values (`CTSEMAnchors`) after
    # the free parameters, where the mean moves and the prior's precision
    # index them; the steps hold them still.
    x = _laplace_extend(laplace, collect(Float64, theta))
    _laplace_check_indices(laplace, length(x))
    # Its scale step reads a full-rank level's covariance off sampled
    # deviations, and a block's integrated ones are not sampled.
    any(_laplace_block, laplace.spec.levels) && throw(ArgumentError(
        "SAEM needs every effect of a full-rank level as a coordinate; this " *
        "objective integrates some in the subject filter (build it without the hybrid)"))
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
_saem_prior_precision(laplace::CTSEMLaplaceObjective, npar::Integer) =
    _ctsem_prior_precision(laplace.objective, npar)

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
    _saem_hold_anchors!(g, B, laplace)
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
        step[(npar - _laplace_nanchors(laplace) + 1):npar] .= 0.0
        step = _saem_ascend!(st, laplace, items, step, g)
        centre && _saem_centre!(st, laplace, prec)
    end
    push!(st.csum, isempty(st.csum) ? copy(st.theta) : st.csum[end] .+ st.theta)
    return (logpost=logpost, gradient_norm=norm(g),
        step=maximum(abs, step; init=0.0), gradient=g, information=B)
end

export ctsem_saem_step!

"""
    _saem_hold_anchors!(g, B, laplace)

Stated population values (`CTSEMAnchors`), the state's last entries, are not
parameters: no score, and an information row and column that leave the
step's solve decoupled from them, so it does not move them.
"""
function _saem_hold_anchors!(g::Vector{Float64}, B::Matrix{Float64},
    laplace::CTSEMLaplaceObjective)
    na = _laplace_nanchors(laplace)
    n = length(g)
    for j in (n - na + 1):n
        g[j] = 0.0
        B[j, :] .= 0.0
        B[:, j] .= 0.0
        B[j, j] = 1.0
    end
    return nothing
end

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
    negligible::Real=0.01, npar::Integer=length(st.theta))
    # `npar`: the free parameters, the first entries of the state; any stated
    # population values after them never move.
    k = length(st.csum)
    q = k ÷ 4
    bl = q ÷ batches
    bl < 5 && return NaN
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
        trend = _saem_trend(st, _saem_curvature(st.info, prec);
            npar=length(st.theta) - _laplace_nanchors(laplace))
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
            estimate[1:length(start)])
        if isfinite(trend) && trend <= null
            settled = true
            break
        end
    end
    k = st.iteration
    # The free parameters, as `start` gave them.
    minimizer = _saem_window(st, k - max(1, k ÷ 2), k)[1:length(start)]
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

################################################################################
# Data augmentation on SAEM's kernel: draws of the joint posterior, placed and
# stopped as NUTS is
################################################################################
#
# This is not SAEM, though it runs on SAEM's state and sweeps and keeps the
# name in the code (`ctsem_saem_sample`, `sampler = 'saem'`). SAEM's iterates
# are not posterior draws -- their spread is the algorithm's Monte Carlo
# noise -- but its draws of the effects are conditional posterior draws, and
# drawing theta given the effects instead of stepping it toward the maximum
# turns the iteration into data augmentation (Tanner and Wong 1987), the
# hierarchical Gibbs sampler, of the joint posterior of theta and u (states
# filtered). It is to SAEM what data augmentation is to EM, and its moves
# pair with the optimiser's the same way: the centred and fixed-u steps of
# AECM become the interweaving of Yu and Meng (2011), a step to a maximum
# becomes a draw. NONMEM's METHOD=BAYES has the same shape (effects by
# Metropolis-Hastings, population parameters drawn given them) and is commonly
# run from an earlier SAEM or importance-sampling estimate, as this is placed
# by SAEM. There is no step-size question here: no stochastic approximation is
# left, and the chains' spread is the posterior's.
# Every move below leaves that posterior invariant:
#
#   * the effects: SAEM's Metropolis-within-Gibbs sweeps (`_saem_sweep!`),
#     their scales adapted in warmup and fixed after it;
#   * theta given u: one NUTS transition (`_nuts_transition!`) on that
#     conditional, metered by its complete-data information at the draws
#     actually visited -- the members' score outer products plus the prior's
#     precision, from each chain's start and then averaged over Stan's warmup
#     windows (`_adapt_windows`). Not the joint sampler's metric: its
#     population block is measured once at the modes, where a population scale
#     direction can curve upward, and on mvmix and gD1 it let NUTS take steps
#     of 0.004 and 0.03 where the information allows 0.25 and 1 -- the chains
#     did not move, R-hat 3.5-4.1. One MALA step, which this took before, is
#     too short a move: on six bench models whose chains had not mixed in 4000
#     draws, five leapfrog steps brought every one to NUTS-level error with
#     seeds agreeing (review/POSTERIOR-race-2026-10-02.md). NUTS chooses the
#     length itself;
#   * each level's mean given the effects b = theta + L u: an exact Gaussian
#     draw along L's columns, the effects re-expressed (the centred half of
#     the interweaving, ASIS: Yu and Meng 2011);
#   * a full-rank level's scales and correlations, two ways. Centred: the
#     deviations L u held and the effects re-expressed, random-walk Metropolis
#     on their own conditional, at no likelihood cost. Non-centred: u held and
#     every member's likelihood evaluated. The first is fast where the data pin
#     the deviations down; the second where they do not -- in a variance's
#     funnel neck the centred move can shrink a variance only as fast as the
#     deviations shrink, which on a study level with eight groups took bursts
#     a few thousand iterations apart (gN3 in the same note);
#   * a reduced-rank level's loadings: L <- L exp(E), E lower triangular and
#     symmetric about zero, the effects mapped by exp(-E), accepted with the
#     map's Jacobian (loadings prod_q A_qq^(k - q + 1), effects |det A|^-G).
#
# It is placed exactly as `ctsem_sample` is: the same per-unit modes, the same
# metric, and each chain's start drawn from the same Laplace approximation of
# the joint posterior (`_sample_initial_point`) on the same stream, so the two
# samplers start chain `c` at the same point. An approximation's whole job here
# is to place the chains. And it stops by the same rule
# (`_sample_until_target`), so a target and a budget mean the same thing
# whichever sampler drew the run.

"""One chain of the SAEM sampler: its state, its tuning, and its draws."""
mutable struct _SAEMChain
    st::CTSEMSAEMState
    # Theta's NUTS transitions and the level moves. The effects' sweeps draw
    # from `st`'s own seeded streams, as SAEM's do.
    rng::Random.Xoshiro
    ws::_NUTSWorkspace
    g::Vector{Float64}
    # Per unit, every member's log likelihood at a trajectory's point; the
    # chain's own stay in `st.ll`.
    trial::Vector{Vector{Float64}}
    eps::Float64
    da::_DualAverage
    # Theta's metric, and the information summed since the current warmup
    # window began.
    metric::CTSEMMetric
    info_sum::Matrix{Float64}
    info_n::Int
    # Per level: the centred scale (or loading) move's and the non-centred
    # move's log step, and the centred move's metric factor where it has one.
    logscale::Vector{Float64}
    logncp::Vector{Float64}
    fisher::Vector{Any}
    draws::Vector{Vector{Float64}}
    effect_sum::Vector{Float64}
    effect_sq::Vector{Float64}
    effect_cross::Vector{Float64}
    accept::Vector{Float64}
    divergent::Vector{Bool}
    depth::Vector{Int}
    energy::Vector{Float64}
    warmup_divergent::Int
    scale_accept::Float64
    nscale::Int
    ncp_accept::Float64
    nncp::Int
end

"""Every member's log likelihood at `theta` and the chain's effects, into
`dest` (per unit); false where one is not finite, or where `theta` has no
population Cholesky factor."""
function _saem_sample_ll!(dest, st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    theta::Vector{Float64})
    Ls = try
        _laplace_popchols(theta, laplace.spec)
    catch err
        _ctsem_must_propagate(err) && rethrow()
        return false
    end
    return _laplace_parallel(laplace, collect(eachindex(st.u))) do U
        _saem_members_ll!(dest[U], laplace, U, eachindex(laplace.units.members[U]),
            theta, Ls, st.u[U][1])
        all(isfinite, dest[U])
    end
end

"""
    _saem_conditional!(g, st, laplace, theta, ll)

`log p(theta | u)` up to a constant -- every member's log likelihood at the
chain's effects, plus the prior -- with its gradient in `g`, from one reverse
sweep per member, and each member's log likelihood in `ll` (one vector per
unit). `-Inf` where a member or the population Cholesky factors cannot be
evaluated: a trajectory routinely steps somewhere the filter refuses, and that
is a rejection, not an error.
"""
function _saem_conditional!(g::Vector{Float64}, st::CTSEMSAEMState,
    laplace::CTSEMLaplaceObjective, free::Vector{Float64}, ll)
    # NUTS moves the free parameters (`g`'s length); the members read them with
    # any stated population values after them (`_laplace_extend`). One
    # assignment: the unit loop's closure captures it.
    theta = _laplace_extend(laplace, free)
    spec = laplace.spec
    Ls, dL = try
        (_laplace_popchols(theta, spec), _laplace_level_chol_derivatives(theta, spec))
    catch err
        _ctsem_must_propagate(err) && rethrow()
        return -Inf
    end
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    ok = _laplace_parallel(laplace, collect(eachindex(st.u))) do U
        isempty(laplace.units.members[U]) && return true
        _saem_unit_scores!(st, laplace, U, 1, Ls, dL, positions; theta=theta,
            values=ll[U])
    end
    ok || return -Inf
    fill!(g, 0.0)
    value = 0.0
    ng = length(g)
    for U in eachindex(st.u)
        S = st.scores[U][1]
        for m in axes(S, 2)
            @views g .+= S[1:ng, m]
        end
        value += sum(ll[U]; init=0.0)
    end
    _ctsem_log_prior_gradient!(g, laplace.objective, theta)
    value += _ctsem_log_prior(laplace.objective, theta)
    return isfinite(value) && all(isfinite, g) ? value : -Inf
end

"""
    _saem_sample_ncp!(st, laplace, l, scale, Rf, rng)

A full-rank level's scales and correlations by random-walk Metropolis with the
standardised effects `u` held -- the non-centred counterpart of
`_saem_sample_scale!`. Every member's likelihood moves, so a proposal costs one
pass over all of them; the proposal is `scale * Rf^-1 z`, `Rf` the centred
move's metric factor. Returns the acceptance probability, NaN for a level with
nothing to move.
"""
function _saem_sample_ncp!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, scale::Float64, Rf::UpperTriangular, rng)
    positions = _laplace_level_positions(laplace.spec, l)
    isempty(positions) && return NaN
    trial = copy(st.theta)
    trial[positions] .+= scale .* (Rf \ randn(rng, length(positions)))
    new = [similar(st.ll[U][1]) for U in eachindex(st.u)]
    _saem_sample_ll!(new, st, laplace, trial) || return 0.0
    logr = sum(sum(new[U]; init=0.0) for U in eachindex(st.u); init=0.0) -
        sum(sum(st.ll[U][1]; init=0.0) for U in eachindex(st.u); init=0.0) +
        _ctsem_log_prior(laplace.objective, trial) -
        _ctsem_log_prior(laplace.objective, st.theta)
    alpha = isfinite(logr) ? min(1.0, exp(logr)) : 0.0
    if rand(rng) < alpha
        st.theta .= trial
        for U in eachindex(st.u)
            st.ll[U][1] .= new[U]
        end
    end
    return alpha
end

"""The level's sites (unit, block) and each one's columns of `u`."""
function _saem_level_sites(laplace::CTSEMLaplaceObjective, l::Integer)
    sites = Tuple{Int,Int}[]
    for U in eachindex(laplace.units.blocks), (b, blk) in enumerate(laplace.units.blocks[U])
        blk.level == l && push!(sites, (U, b))
    end
    return sites
end

_saem_cols(laplace, U, b) = (laplace.units.blocks[U][b].offset + 1):(laplace.units.blocks[U][b].offset +
    laplace.units.blocks[U][b].size)

"""
    _saem_sample_mean!(st, laplace, l, prec, rng)

Level `l`'s mean, and the TI coefficients on its parameters, given its effects
`b`: their conditional (`_saem_mean_system`) is Gaussian; draw from it, move
theta by `L delta` and the coefficients by `eta`, and every group's effects by
`-delta - N_g eta`. Members' likelihoods do not move.
"""
function _saem_sample_mean!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, prec::Vector{Float64}, rng)
    level = laplace.spec.levels[l]
    r = nlatent(level)
    sites = _saem_level_sites(laplace, l)
    G = length(sites)
    (G == 0 || r == 0) && return nothing
    ti = _saem_level_ti(laplace, l)
    system = _saem_mean_system(st, laplace, l, sites, prec, 1, ti)
    system === nothing && (ti = _saem_level_ti_none;
        system = _saem_mean_system(st, laplace, l, sites, prec, 1, ti))
    A, rhs, L, N, Z = system
    isempty(rhs) && return nothing
    F = cholesky(Symmetric(A))
    z = (F \ rhs) .+ (F.U \ randn(rng, length(rhs)))
    rz = Z === nothing ? r : size(Z, 2)
    delta = Z === nothing ? z[1:r] : Z * z[1:rz]
    _saem_mean_apply!(st, laplace, l, sites, 1, ti, L, N, delta, z[(rz + 1):end])
    return nothing
end

const _saem_level_ti_none = Tuple{Int,Vector{Tuple{Int,Int}}}[]

"""
    _saem_level_ti(laplace, l)

The TI coefficients a centred move can carry with level `l`'s mean, each as
`(coefficient, [(predictor, row), ...])`: its raw index, and each effect it
has, `row` the parameter's place in `re_index`. A coefficient moves its
parameters by its change times each group's predictor, and the move takes that
back out of the group's effects -- which needs the shift of one parameter to
be one of `u` (a full-rank level with no own variances), one predictor value
per group (a level whose groups are single subjects), and every parameter the
coefficient reaches to be one of the level's, or that parameter's members
would move. Fixed effects have no coefficient. Empty where nothing qualifies.
"""
function _saem_level_ti(laplace::CTSEMLaplaceObjective, l::Integer)
    sp = laplace.objective.params
    out = Tuple{Int,Vector{Tuple{Int,Int}}}[]
    isempty(sp.ti_parameter_indices) && return out
    level = laplace.spec.levels[l]
    (!isreducedrank(level) && nlatent(level) == nrandomeffects(level) == level.rank) ||
        return out
    for U in eachindex(laplace.units.blocks), blk in laplace.units.blocks[U]
        blk.level == l && length(blk.members) != 1 && return out
    end
    for coef in unique(sp.ti_coefficient_indices)
        coef == 0 && continue
        effects = Tuple{Int,Int}[]
        for i in eachindex(sp.ti_coefficient_indices)
            sp.ti_coefficient_indices[i] == coef || continue
            row = findfirst(==(sp.ti_parameter_indices[i]), level.re_index)
            row === nothing && (empty!(effects); break)
            push!(effects, (sp.ti_predictor_indices[i], row))
        end
        isempty(effects) || push!(out, (coef, effects))
    end
    return out
end

"""
    _saem_mean_system(st, laplace, l, sites, prec, chains, ti)

The Gaussian that level `l`'s mean and the TI coefficients `ti`
(`_saem_level_ti`) have given the effects, in the centred move's coordinates:
`delta` along `L`'s columns for the mean and `eta` for the coefficients, which
shift group `g`'s effects by `-delta - N_g eta` (column `c` of `N_g` is
`L^-1` applied to coefficient `c`'s shift of the parameters for one unit of
it, at that group's predictor values), so no member's parameters move. The
effects are averaged over the first `chains` chains. Returns `(A, b, L, N)`,
`A` the precision and `A^-1 b` the mean, or `nothing` where `A` is singular --
a coefficient its predictor cannot tell from the mean, with no prior.
"""
function _saem_mean_system(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, sites, prec::Vector{Float64}, chains::Integer, ti)
    level = laplace.spec.levels[l]
    r = nlatent(level)
    C = length(ti)
    G = length(sites)
    L = _laplace_popchol(st.theta, level)
    Dm = prec[level.re_index]
    mu = st.theta[level.re_index]
    A = zeros(r + C, r + C)
    b = zeros(r + C)
    A[1:r, 1:r] .= transpose(L) * (Dm .* L)
    for q in 1:r
        A[q, q] += G
    end
    b[1:r] .= .-(transpose(L) * (Dm .* mu))
    N = Vector{Matrix{Float64}}(undef, C == 0 ? 0 : G)
    Linv = C == 0 ? L : inv(L)
    ubar = zeros(r)
    for (g, (U, blk)) in enumerate(sites)
        fill!(ubar, 0.0)
        for c in 1:chains
            ubar .+= @view st.u[U][c][_saem_cols(laplace, U, blk)]
        end
        ubar ./= chains
        @views b[1:r] .+= ubar
        C == 0 && continue
        member = laplace.units.members[U][laplace.units.blocks[U][blk].members[1]]
        tipreds = laplace.objective.subject_objectives[member].tipreds
        Ng = zeros(r, C)
        for (c, (_, effects)) in enumerate(ti), (predictor, row) in effects
            x = _ctsem_tipred_value(tipreds, predictor, st.theta)
            @views Ng[:, c] .+= x .* Linv[:, row]
        end
        N[g] = Ng
        @views A[1:r, (r + 1):end] .+= Ng
        @views A[(r + 1):end, (r + 1):end] .+= transpose(Ng) * Ng
        @views b[(r + 1):end] .+= transpose(Ng) * ubar
    end
    for (c, (coef, _)) in enumerate(ti)
        A[r + c, r + c] += prec[coef]
        b[r + c] -= prec[coef] * st.theta[coef]
    end
    C == 0 && return _saem_mean_confine(laplace, level, A, b, L, N)
    A[(r + 1):end, 1:r] .= transpose(A[1:r, (r + 1):end])
    issuccess(cholesky(Symmetric(A); check=false)) || return nothing
    return _saem_mean_confine(laplace, level, A, b, L, N)
end

"""
    _saem_mean_confine(laplace, level, A, b, L, N)

The mean system `(A, b, L, N)` with `Z` beside it: `nothing` where none of the
level's parameters has a stated population value (`CTSEMAnchors`), the system
unchanged; otherwise an orthonormal basis of the moves `delta` whose shift
`L delta` leaves every such parameter where it is, and the system in the
coordinates `(y, eta)` of `delta = Z y`. A Gaussian confined to that subspace
is the one in `y`, so the centred step and the draw are each taken there.
"""
function _saem_mean_confine(laplace::CTSEMLaplaceObjective, level, A, b, L, N)
    anchors = laplace.objective.params.anchors
    nf = isempty(anchors.values) ? typemax(Int) : anchors.nfree
    held = findall(>(nf), level.re_index)
    isempty(held) && return (A, b, L, N, nothing)
    r = size(L, 2)
    C = size(A, 1) - r
    Z = nullspace(L[held, :])
    T = zeros(r + C, size(Z, 2) + C)
    T[1:r, 1:size(Z, 2)] .= Z
    for c in 1:C
        T[r + c, size(Z, 2) + c] = 1.0
    end
    return (transpose(T) * A * T, transpose(T) * b, L, N, Z)
end

"""Move level `l`'s mean by `L delta`, the TI coefficients `ti` by `eta`, and
the first `chains` chains' effects at every site by `-delta - N_g eta`."""
function _saem_mean_apply!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, sites, chains::Integer, ti, L, N, delta, eta)
    level = laplace.spec.levels[l]
    st.theta[level.re_index] .+= L * delta
    for (c, (coef, _)) in enumerate(ti)
        st.theta[coef] += eta[c]
    end
    # A confined move (`_saem_mean_confine`) leaves stated population values
    # where they were to rounding; put them back exactly.
    na = _laplace_nanchors(laplace)
    if na > 0
        n = length(st.theta)
        st.theta[(n - na + 1):n] .= laplace.objective.params.anchors.values
    end
    for (g, (U, blk)) in enumerate(sites)
        shift = isempty(ti) ? delta : delta .+ N[g] * eta
        for c in 1:chains
            st.u[U][c][_saem_cols(laplace, U, blk)] .-= shift
        end
    end
    return st
end

"""
    _saem_sample_scale!(st, laplace, l, prec, scale, Rf, rng)

A full-rank level's scales and correlations `phi` given the deviations
`d = L u`: one random-walk Metropolis step on their conditional
`log p(phi) - G sum log diag L(phi) - sum |L(phi)^-1 d|^2 / 2`, proposal
`scale * Rf^-1 z` (`Rf` the upper factor of the conditional's Fisher
information, frozen), the effects re-expressed as `L(phi')^-1 d` when accepted.
Returns the acceptance probability.
"""
function _saem_sample_scale!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, prec::Vector{Float64}, scale::Float64, Rf::UpperTriangular, rng)
    spec = laplace.spec
    level = spec.levels[l]
    positions = _laplace_level_positions(spec, l)
    sites = _saem_level_sites(laplace, l)
    G = length(sites)
    (isempty(positions) || G == 0) && return NaN
    L0 = _laplace_popchol(st.theta, level)
    k = nrandomeffects(level)
    d = zeros(k, G)
    for (j, (U, b)) in enumerate(sites)
        d[:, j] = L0 * st.u[U][1][_saem_cols(laplace, U, b)]
    end
    Dp = prec[positions]
    function f(theta)
        L = _laplace_popchol(theta, level)
        all(isfinite, L) || return -Inf
        logdet, sign = logabsdet(L)
        (sign == 0 || !isfinite(logdet)) && return -Inf
        phi = theta[positions]
        return -G * logdet - sum(abs2, L \ d) / 2 - sum(Dp .* phi .^ 2) / 2
    end
    trial = copy(st.theta)
    trial[positions] .+= scale .* (Rf \ randn(rng, length(positions)))
    f0 = f(st.theta)
    f1 = f(trial)
    alpha = isfinite(f1) ? min(1.0, exp(f1 - f0)) : 0.0
    if rand(rng) < alpha
        st.theta .= trial
        L1 = _laplace_popchol(st.theta, level)
        for (j, (U, b)) in enumerate(sites)
            st.u[U][1][_saem_cols(laplace, U, b)] = L1 \ d[:, j]
        end
    end
    return alpha
end

"""
    _saem_sample_loadings!(st, laplace, l, prec, scale, rng)

A reduced-rank level's loadings: propose `A = exp(E)`, `E` lower triangular
with `N(0, scale^2)` entries, so the reverse move is `exp(-E)` with the same
density; map loadings `R -> R A` and every effect `u -> A^-1 u`, which leaves
each deviation `L u` and so every member's likelihood alone, and accept with
the prior and effects' density ratio times the map's Jacobian. Returns the
acceptance probability.
"""
function _saem_sample_loadings!(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, prec::Vector{Float64}, scale::Float64, rng)
    level = laplace.spec.levels[l]
    # As `_saem_centre_loadings!`: the move is for the triangular loadings only.
    nunique(level) > 0 && return NaN
    sites = _saem_level_sites(laplace, l)
    G = length(sites)
    k = nrandomeffects(level)
    r = nlatent(level)
    (G == 0 || r == 0) && return NaN
    # `E` on the entries the level's pattern allows (`_laplace_loading_allowed`),
    # so `exp(E)` keeps the loadings on their pattern.
    allowed = _laplace_loading_allowed(level)
    E = zeros(r, r)
    for q in 1:r, p in q:r
        allowed[p, q] || continue
        E[p, q] = scale * randn(rng)
    end
    A = exp(E)
    Ainv = exp(-E)
    rows, cols = _laplace_loading_entries(level)
    Rl = zeros(k, r)
    D = zeros(k, r)
    for t in eachindex(rows)
        Rl[rows[t], cols[t]] = st.theta[level.load_index[t]]
        D[rows[t], cols[t]] = prec[level.load_index[t]]
    end
    Rn = Rl * A
    uold = [st.u[U][1][_saem_cols(laplace, U, b)] for (U, b) in sites]
    unew = [Ainv * v for v in uold]
    # The map's Jacobian on the loadings: column `q`'s entries scale by
    # `A[q, q]`, one factor each (`k - q + 1` of them in the triangle).
    ncol = zeros(Int, r)
    for q in cols
        ncol[q] += 1
    end
    logratio = -sum(D .* Rn .^ 2) / 2 + sum(D .* Rl .^ 2) / 2 -
        sum(sum(abs2, v) for v in unew) / 2 + sum(sum(abs2, v) for v in uold) / 2 +
        sum(ncol[q] * E[q, q] for q in 1:r) - G * sum(E[q, q] for q in 1:r)
    alpha = isfinite(logratio) ? min(1.0, exp(logratio)) : 0.0
    if rand(rng) < alpha
        for t in eachindex(rows)
            st.theta[level.load_index[t]] = Rn[rows[t], cols[t]]
        end
        for (j, (U, b)) in enumerate(sites)
            st.u[U][1][_saem_cols(laplace, U, b)] = unew[j]
        end
    end
    return alpha
end

"""
    _saem_scale_fisher(st, laplace, l, prec)

The upper Cholesky factor of a full-rank level's scale conditional's Fisher
information at the chain's theta (`G/2 tr((M_p + M_p')(M_q + M_q'))` plus the
prior), the random-walk metric of `_saem_sample_scale!`.
"""
function _saem_scale_fisher(st::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    l::Integer, prec::Vector{Float64})
    spec = laplace.spec
    positions = _laplace_level_positions(spec, l)
    G = length(_saem_level_sites(laplace, l))
    np = length(positions)
    L = _laplace_popchol(st.theta, spec.levels[l])
    dL = _laplace_level_chol_derivatives(st.theta, spec, l)
    M = [Matrix(L \ dL[t]) for t in 1:np]
    F = zeros(np, np)
    for t1 in 1:np, t2 in t1:np
        F[t1, t2] = F[t2, t1] = G / 2 * tr((M[t1] + transpose(M[t1])) * (M[t2] + transpose(M[t2])))
    end
    F .+= Diagonal(prec[positions] .+ 1e-6 .* (1 .+ diag(F)))
    return cholesky(Symmetric(F)).U
end

"""
    _saem_new_chain(base, laplace, sampler, x, c, seed, rng, setup, maxdepth)

Chain `c`, started at the joint point `x` (`[theta; u_1; u_2; ...]`, the
sampler's layout): a one-chain copy of `base` with `x`'s effects in it, its
members' log likelihoods there, its own streams, and its centred moves'
metrics measured where it starts.
"""
function _saem_new_chain(base::CTSEMSAEMState, laplace::CTSEMLaplaceObjective,
    sampler::CTSEMSampler, x::Vector{Float64}, c::Integer, seed::Integer,
    rng::Random.Xoshiro, setup, maxdepth::Integer)
    npar = sampler.npar
    st = deepcopy(base)
    st.theta = _laplace_extend(laplace, x[1:npar])
    for U in 1:sampler.nunits
        st.u[U][1] .= view(x, _sample_urange(sampler, U))
    end
    st.seed = UInt64(hash((seed + c, :saem_sample)))
    st.iteration = 0
    for U in eachindex(st.u)
        fill!(st.accepted[U][1], 0)
        fill!(st.proposed[U][1], 0)
    end
    _saem_sample_ll!([st.ll[U][1] for U in eachindex(st.u)], st, laplace, st.theta) ||
        throw(DomainError(st.theta, "SAEM sampler: a member's log likelihood is " *
            "not finite where chain $c starts"))
    levels = setup.levels
    fisher = Any[(setup.active[l] && !setup.reduced[l] &&
                  !isempty(_laplace_level_positions(laplace.spec, l))) ?
                 _saem_scale_fisher(st, laplace, l, setup.prec) : nothing for l in levels]
    logscale = [setup.reduced[l] ?
        log(1 / sqrt(2 * max(1, length(_saem_level_sites(laplace, l))))) : log(0.5)
        for l in levels]
    neffects = last(_laplace_deviation_offsets(laplace))
    return _SAEMChain(st, rng, _NUTSWorkspace(npar, Int(maxdepth)), zeros(npar),
        [similar(st.ll[U][1]) for U in eachindex(st.u)], 0.0, _DualAverage(1.0, 0.8),
        ctsem_identity_metric(npar), zeros(length(st.theta), length(st.theta)), 0,
        logscale, fill(log(0.5), length(levels)), fisher, Vector{Float64}[],
        zeros(neffects), zeros(neffects),
        zeros(_effect_cross_length(_laplace_deviation_blocks(laplace))),
        Float64[], Bool[], Int[], Float64[], 0, 0.0, 0, 0.0, 0)
end

"""
    _saem_one_chain(base, c)

A one-chain copy of SAEM's state `base`, holding its chain `c`'s draws, tuned
scales and counts (cycling when there are fewer chains than asked for); `base`
itself when it has one chain.
"""
function _saem_one_chain(base::CTSEMSAEMState, c::Integer)
    base.chains == 1 && return base
    kc = mod1(c, base.chains)
    st = deepcopy(base)
    st.u = [[copy(base.u[U][kc])] for U in eachindex(base.u)]
    st.ll = [[copy(base.ll[U][kc])] for U in eachindex(base.ll)]
    st.logscale = [[copy(base.logscale[U][kc])] for U in eachindex(base.logscale)]
    st.accepted = [[copy(base.accepted[U][kc])] for U in eachindex(base.accepted)]
    st.proposed = [[copy(base.proposed[U][kc])] for U in eachindex(base.proposed)]
    st.indep_accepted = [zeros(Int, 1) for _ in eachindex(base.u)]
    st.indep_proposed = [zeros(Int, 1) for _ in eachindex(base.u)]
    st.scores = [[copy(base.scores[U][kc])] for U in eachindex(base.scores)]
    st.chains = 1
    return st
end

"""The joint point `[theta; u_1; u_2; ...]` from `theta` and a one-chain
state's effects."""
function _saem_joint_point(sampler::CTSEMSampler, theta::AbstractVector, st::CTSEMSAEMState)
    x = zeros(sampler.ndim)
    x[1:sampler.npar] .= view(theta, 1:sampler.npar)
    for U in 1:sampler.nunits
        x[_sample_urange(sampler, U)] .= st.u[U][1]
    end
    return x
end

"""
    _saem_placement(laplace, values; nchains, seed, nestep, init_scale)

Where the joint samplers start when asked to (`saem = true`): SAEM from
`values` to its own stop, on the exact marginal posterior rather than the
Laplace approximation to it, with at least as many chains as the sampler will
run so each starts from different effects. Its state's draws are then settled
by `nestep` E-steps at the averaged estimate -- the run ends on its last
iterate, not on the average it returns. Each chain then starts at the estimate
jittered by `init_scale` times one draw from the state's complete-data
information (plus the prior), with its own SAEM chain's effects -- the
samplers' `init_scale`, which would otherwise have nothing to scale once the
chains are placed. Returns the estimate, the state,
the starts (the joint vectors and their parameter parts) and SAEM's
summary.
"""
function _saem_placement(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    nchains::Integer=4, seed::Integer=1, nestep::Integer=50, init_scale::Real=1.0)
    t0 = time()
    npar = length(values)
    nunits = length(laplace.units.members)
    run = ctsem_saem(laplace, collect(Float64, values); seed=seed,
        chains=max(Int(nchains), _saem_default_chains(nunits)))
    st = run.state
    theta = collect(Float64, run.minimizer)
    st.theta .= _laplace_extend(laplace, theta)
    prec = _saem_prior_precision(laplace, length(st.theta))
    for _ in 1:nestep
        ctsem_saem_step!(st, laplace; prec=prec, mstep=false)
    end
    sampler = ctsem_sampler(laplace, npar)
    # The free parameters' block: any stated population values after them are
    # held, not drawn.
    R = cholesky(_saem_curvature(st.info[1:npar, 1:npar], prec[1:npar])).U
    rng = Random.Xoshiro(hash((seed, :saem_placement)))
    joint = zeros(sampler.ndim, nchains)
    pars = zeros(npar, nchains)
    for c in 1:nchains
        th = theta .+ Float64(init_scale) .* (R \ randn(rng, npar))
        pars[:, c] = th
        joint[:, c] = _saem_joint_point(sampler, th, _saem_one_chain(st, c))
    end
    return (theta=theta, state=st, joint=joint, pars=pars, iterations=run.iterations,
        settled=run.settled, trend=run.trend, secs=time() - t0)
end

"""The chain's own members' log likelihoods, one vector per unit (the arrays
in `st.ll`, not copies)."""
_saem_own_ll(st::CTSEMSAEMState) = [st.ll[U][1] for U in eachindex(st.u)]

"""The members' score outer products at the chain's draw, summed: the
complete-data information about theta given u, from the scores the last
`_saem_conditional!` left in `st.scores`."""
function _saem_score_information(st::CTSEMSAEMState)
    npar = length(st.theta)
    B = zeros(npar, npar)
    for U in eachindex(st.u)
        S = st.scores[U][1]
        size(S, 2) == 0 && continue
        mul!(B, S, transpose(S), 1.0, 1.0)
    end
    return B
end

"""Meter theta's transitions by the information `B` plus the prior's
precision (`_saem_curvature`)."""
function _saem_set_metric!(ch::_SAEMChain, B::AbstractMatrix, prec::Vector{Float64})
    # The free parameters NUTS moves (`ch.g`'s length), the first of the
    # state's; any stated population values after them are held.
    n = length(ch.g)
    ch.metric = _metric_from_covariances([1:n],
        [Matrix(inv(_saem_curvature(B[1:n, 1:n], prec[1:n])))])
    return ch
end

"""
    _saem_chain_iterate!(ch, laplace, setup, k, warm)

One iteration of chain `ch`: theta given u, then the effects, then each
level's mean and scales. `k` counts warmup iterations for the adaptation, which
stops once `warm` is false. Returns the theta transition's summary.
"""
function _saem_chain_iterate!(ch::_SAEMChain, laplace::CTSEMLaplaceObjective,
    setup, k::Integer, warm::Bool)
    st = ch.st
    st.iteration += 1
    own = _saem_own_ll(st)
    # The effects moved since the last transition, so the start's density and
    # gradient are evaluated afresh; the same sweep refreshes `st.ll`.
    logp = _saem_conditional!(ch.g, st, laplace, st.theta, own)
    isfinite(logp) || throw(DomainError(st.theta, "SAEM sampler: the parameters' " *
        "conditional density is not finite at the chain's own point"))
    # The same sweep left the members' scores here, which is the information
    # the metric is re-estimated from at each warmup window's end.
    if warm && k > _ADAPT_INIT_BUFFER
        ch.info_sum .+= _saem_score_information(st)
        ch.info_n += 1
    end
    # The free parameters (`ch.g`'s length); any stated population values
    # after them in the state are held.
    np = length(ch.g)
    x = st.theta[1:np]
    trial = ch.trial
    step = _nuts_transition!(ch.ws, (g, y) -> _saem_conditional!(g, st, laplace, y, trial),
        ch.metric, ch.rng, x, ch.g, logp, ch.eps, setup.maxdepth, setup.maxdelta)
    if x != view(st.theta, 1:np)
        # The transition chose one point of its trajectory and its members'
        # likelihoods went with the others; one forward pass recovers them.
        st.theta[1:np] .= x
        _saem_sample_ll!(own, st, laplace, st.theta) ||
            throw(DomainError(st.theta, "SAEM sampler: a member's log likelihood " *
                "is not finite at the point the transition accepted"))
    end
    if warm
        ch.eps = _dual_update!(ch.da, step.accept)
        ch.warmup_divergent += step.divergent
    end

    Ls = _laplace_popchols(st.theta, laplace.spec)
    rate = warm ? 1.0 / (1 + k)^0.6 : 0.0
    ok = _laplace_parallel(laplace, collect(eachindex(st.u))) do U
        isempty(st.u[U][1]) && return true
        for sweep in 1:setup.sweeps
            _saem_sweep!(st, laplace, U, 1, Ls, sweep, rate; nupper=setup.nupper)
        end
        true
    end
    ok || throw(DomainError(st.theta, "SAEM sampler: a likelihood is not finite"))

    for l in setup.levels
        setup.active[l] || continue
        _saem_sample_mean!(st, laplace, l, setup.prec, ch.rng)
        for _ in 1:setup.nscale
            s = exp(ch.logscale[l])
            a = setup.reduced[l] ? _saem_sample_loadings!(st, laplace, l, setup.prec, s, ch.rng) :
                ch.fisher[l] === nothing ? NaN :
                _saem_sample_scale!(st, laplace, l, setup.prec, s, ch.fisher[l], ch.rng)
            isfinite(a) || continue
            if warm
                ch.logscale[l] += (a - 0.3) / k^0.6
            else
                ch.scale_accept += a
                ch.nscale += 1
            end
        end
        (setup.reduced[l] || ch.fisher[l] === nothing) && continue
        for _ in 1:setup.nncp
            a = _saem_sample_ncp!(st, laplace, l, exp(ch.logncp[l]), ch.fisher[l], ch.rng)
            isfinite(a) || continue
            if warm
                ch.logncp[l] += (a - 0.3) / k^0.6
            else
                ch.ncp_accept += a
                ch.nncp += 1
            end
        end
    end
    return step
end

"""Keep the chain's current draw: theta (and the effects, with
`save_effects`), the effects' running sums, and the transition's summary."""
function _saem_chain_record!(ch::_SAEMChain, sampler::CTSEMSampler, step,
    save_effects::Bool)
    st = ch.st
    npar = sampler.npar
    laplace = sampler.laplace
    # Natural deviations, as `ctsem_sample` reports them (`_laplace_deviations!`).
    devstarts, ndev = _laplace_deviation_offsets(laplace)
    effects = reduce(vcat, (st.u[U][1] for U in 1:sampler.nunits); init=Float64[])
    dev = _laplace_deviations!(zeros(ndev), laplace,
        _laplace_popchols(st.theta, laplace.spec), effects, devstarts)
    x = save_effects ? zeros(npar + ndev) : st.theta[1:npar]
    save_effects && (x[1:npar] .= view(st.theta, 1:npar))
    @inbounds for j in 1:ndev
        ch.effect_sum[j] += dev[j]
        ch.effect_sq[j] += dev[j] * dev[j]
        save_effects && (x[npar + j] = dev[j])
    end
    _effect_cross_add!(ch.effect_cross, dev, _laplace_deviation_blocks(laplace))
    push!(ch.draws, x)
    push!(ch.accept, step.accept)
    push!(ch.divergent, step.divergent)
    push!(ch.depth, step.depth)
    push!(ch.energy, step.energy)
    return ch
end

"""
    _saem_chain_run!(ch, laplace, sampler, setup, niter, nwarm, progress,
        callback, save_effects)

`niter` iterations of one chain, the first `nwarm` of them adapting and
discarded. At each of Stan's warmup windows' ends (`_adapt_windows`) theta's
metric is re-estimated from the information averaged over the window and the
step size's averaging restarts around it; the averaged step size takes over
when warmup ends.
"""
function _saem_chain_run!(ch::_SAEMChain, laplace::CTSEMLaplaceObjective,
    sampler::CTSEMSampler, setup, niter::Int, nwarm::Int,
    progress::CTSEMProgress, callback::CTSEMCallback, save_effects::Bool)
    ndiv = 0
    windows = nwarm > 0 ? _adapt_windows(nwarm) : Int[]
    for k in 1:niter
        warm = k <= nwarm
        step = _saem_chain_iterate!(ch, laplace, setup, k, warm)
        if warm
            if k in windows && ch.info_n > 0
                _saem_set_metric!(ch, ch.info_sum ./ ch.info_n, setup.prec)
                fill!(ch.info_sum, 0.0)
                ch.info_n = 0
                _dual_restart!(ch.da, ch.eps)
            end
            k == nwarm && (ch.eps = _dual_final(ch.da))
        else
            _saem_chain_record!(ch, sampler, step, save_effects)
            ndiv += step.divergent
        end
        phase = warm ? "warmup" : "sampling"
        done, total = warm ? (k, nwarm) : (k - nwarm, niter - nwarm)
        # The rate is per phase, as `_sample_draws!` reports it.
        (nwarm > 0 && k == nwarm + 1) && (progress.started = time())
        if _due(progress)
            progress.label = phase
            _progress_line(progress, done, total, @sprintf("logp|u %11.2f", step.logp),
                @sprintf("depth %d", step.depth), @sprintf("div %d", ndiv))
        end
        _invoke_callback(callback, phase, done, total, step.logp, ndiv)
    end
    return ch
end

"""
    ctsem_saem_sample(laplace, values; kwargs...)

Sample the joint posterior over population parameters and random effects by
data augmentation on SAEM's kernel (see the header above): the effects by
SAEM's sweeps, theta given them by NUTS, the population means and scales by
moves that re-express the effects. The counterpart of `ctsem_sample`, which runs NUTS on the whole joint
vector at once; this one never forms that vector's trajectory, so its cost per
iteration grows with the subjects' filters rather than with a joint tree.

Arguments and result follow `ctsem_sample`: `values` (a Laplace estimate)
places the chains and the metric, `nchains`/`nwarmup`/`ndraws`/`seed`/
`init_scale`/`hessian`/`save_effects` mean what they mean there, `maxdepth`,
`target_accept`, `maxdelta` and `stepsize` tune theta's NUTS transitions, and
`min_ess`/`mean_ess`/`max_draws`/`rhat_target` are the same targets, stopped
on by the same rule. `sweeps` and `nupper` are the effects' sweeps per
iteration and collapsed moves per block above a leaf, `nscale` the centred
scale moves and `nncp` the non-centred ones per level per iteration.

`starts`, as for `ctsem_sample`, gives each chain its own starting point.
`state`, a SAEM state (`ctsem_saem`'s `state`, or `ctsem_saem_init` after
E-steps), continues SAEM's chains instead of starting fresh: chain `c` takes
the state's chain `c`'s effects (cycling), its tuned proposal scales and
shapes, and the state's averaged complete-data information as theta's first
metric, with theta from `starts` when given (parameters only) and `values`
otherwise. `saem = true` makes that state itself first (`_saem_placement`):
SAEM from `values`, then its chains continued, so nothing from the Laplace
approximation's optimum is needed to place the sampler.

The result has `ctsem_sample`'s fields, with `sampler = "saem"`; `depth`,
`divergent`, `energy` and `stepsize` describe theta's transitions, and `ebfmi`
is NaN, since an energy taken at a different `u` every iteration does not
measure what E-BFMI asks.
"""
function _ctsem_saem_sample(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    npar::Integer=length(values), nchains::Integer=4, nwarmup::Integer=500,
    ndraws::Integer=500, maxdepth::Integer=10, target_accept::Real=0.8,
    maxdelta::Real=1000.0, seed::Integer=20260828, init_scale::Real=1.0,
    save_effects::Bool=false, hessian::Union{Nothing,AbstractMatrix}=nothing,
    verbose::Bool=false, min_ess::Real=0.0, mean_ess::Real=0.0,
    max_draws::Integer=0, rhat_target::Real=1.01, progress_overwrite::Bool=true,
    progress_callback=nothing, progress_sink=nothing, stepsize::Real=0.0,
    sweeps::Integer=2, nupper::Integer=2, nscale::Integer=5, nncp::Integer=2,
    starts::Union{Nothing,AbstractMatrix}=nothing, state=nothing, saem::Bool=false,
    saem_nestep::Integer=50, coordinate=nothing)

    t0 = time()
    nchains = Int(nchains); nwarmup = Int(nwarmup); ndraws = Int(ndraws)
    placement = nothing
    if saem
        state === nothing || throw(ArgumentError(
            "saem = true makes its own state; pass state or saem, not both"))
        placement = _saem_placement(laplace, collect(Float64, values)[1:npar];
            nchains=nchains, seed=seed, nestep=saem_nestep, init_scale=init_scale)
        values = placement.theta
        state = placement.state
        starts === nothing && (starts = placement.pars)
    end
    _check_sample_arguments(nchains, ndraws, target_accept; nwarmup=nwarmup)

    # Placed as `ctsem_sample` places its chains, in the same order and for
    # the same reasons (see there): the metric solves the modes the start
    # reads.
    sampler = ctsem_sampler(laplace, npar)
    start = collect(Float64, values)
    metric = ctsem_sample_metric(sampler, start; hessian=hessian)
    centre = ctsem_sample_start(sampler, start)
    npar = sampler.npar
    base = state === nothing ?
        ctsem_saem_init(laplace, centre[1:npar]; seed=seed, chains=1, disperse=false) : state

    levels = collect(eachindex(laplace.spec.levels))
    setup = (
        maxdepth=Int(maxdepth), maxdelta=Float64(maxdelta),
        sweeps=Int(sweeps), nupper=Int(nupper), nscale=Int(nscale), nncp=Int(nncp),
        prec=_saem_prior_precision(laplace, npar + _laplace_nanchors(laplace)),
        levels=levels,
        active=[nrandomeffects(laplace.spec.levels[l]) > 0 &&
                nlatent(laplace.spec.levels[l]) > 0 for l in levels],
        reduced=[isreducedrank(laplace.spec.levels[l]) for l in levels])

    gjoint = zeros(sampler.ndim)
    density = (g, x) -> ctsem_sample_density!(g, sampler, x)
    chains = Vector{_SAEMChain}(undef, nchains)
    jointstarts = _sample_joint_starts(starts, centre, npar)
    for c in 1:nchains
        rng = Random.Xoshiro(UInt64(seed) + UInt64(c))
        chainbase = _saem_one_chain(base, c)
        x = if state !== nothing && (starts === nothing || size(starts, 1) == npar)
            _saem_joint_point(sampler, starts === nothing ? centre[1:npar] :
                view(starts, :, mod1(c, size(starts, 2))), chainbase)
        elseif jointstarts !== nothing
            jointstarts[mod1(c, length(jointstarts))]
        else
            first(_sample_initial_point(centre, metric, rng, Float64(init_scale),
                density, gjoint))
        end
        ch = _saem_new_chain(chainbase, laplace, sampler, x, c, seed, rng, setup, maxdepth)
        logp = _saem_conditional!(ch.g, ch.st, laplace, ch.st.theta, _saem_own_ll(ch.st))
        isfinite(logp) || throw(DomainError(ch.st.theta, "SAEM sampler: the " *
            "parameters' conditional density is not finite where chain $c starts"))
        _saem_set_metric!(ch, (state === nothing || all(iszero, state.info)) ?
            _saem_score_information(ch.st) : state.info, setup.prec)
        ch.eps = stepsize > 0 ? Float64(stepsize) :
            _init_stepsize((g, y) -> _saem_conditional!(g, ch.st, laplace, y, ch.trial),
                ch.metric, rng, ch.st.theta[1:npar], copy(ch.g), logp, ch.ws)
        ch.da = _DualAverage(ch.eps, Float64(target_accept))
        chains[c] = ch
    end

    width = _laplace_ensure_pool!(laplace)
    parallel = nchains > 1 && Threads.nthreads() > 1
    verbose && println(_console(), "Sampling (SAEM kernel): ", nchains, " chain(s), ",
        sampler.npar, " population parameters by NUTS given ", sampler.ndim - sampler.npar,
        " effects, ", width, " worker(s) shared")
    run_all = function (niter, nwarm)
        body = function (c)
            reporter = CTSEMProgress(verbose && c == 1;
                label=nwarm > 0 ? "warmup" : "sampling", overwrite=progress_overwrite,
                sink=progress_sink)
            watcher = CTSEMCallback(c == 1 ? progress_callback : nothing)
            _saem_chain_run!(chains[c], laplace, sampler, setup, niter, nwarm,
                reporter, watcher, save_effects)
            _progress_done(reporter, nwarm > 0 ?
                @sprintf("%d warmup + %d draws", nwarm, niter - nwarm) :
                @sprintf("%d draws", niter))
            return true
        end
        if parallel
            _laplace_parallel(body, laplace, collect(1:nchains))
        else
            for c in 1:nchains
                body(c)
            end
        end
        return nothing
    end
    run_all(nwarmup + ndraws, nwarmup)

    pool = function ()
        n = length(chains[1].draws)
        P = Matrix{Float64}(undef, npar, nchains * n)
        for c in 1:nchains, t in 1:n
            @inbounds P[:, (c - 1) * n + t] .= view(chains[c].draws[t], 1:npar)
        end
        return P
    end
    run = _sample_until_target(pool, (wanted, attempt) -> run_all(wanted, 0),
        nchains, ndraws, Float64(min_ess), Float64(mean_ess),
        max(Int(max_draws), ndraws), Float64(rhat_target), verbose; t0=t0,
        coordinate=coordinate)
    total = run.total

    neffects = last(_laplace_deviation_offsets(laplace))
    kept = save_effects ? npar + neffects : npar
    draws = Matrix{Float64}(undef, kept, nchains * total)
    for c in 1:nchains, t in 1:total
        @inbounds draws[:, (c - 1) * total + t] .= chains[c].draws[t]
    end
    nall = nchains * total
    effect_sum = sum(ch.effect_sum for ch in chains; init=zeros(neffects))
    effect_sq = sum(ch.effect_sq for ch in chains; init=zeros(neffects))
    blocks = _laplace_deviation_blocks(laplace)
    effect_cross = sum(ch.effect_cross for ch in chains;
        init=zeros(_effect_cross_length(blocks)))
    pooled = (draws=draws, accept=reduce(vcat, [ch.accept for ch in chains]),
        divergent=reduce(vcat, [ch.divergent for ch in chains]),
        depth=reduce(vcat, [ch.depth for ch in chains]),
        energy=reduce(vcat, [ch.energy for ch in chains]))
    # E-BFMI is not computed: an energy taken at different effects every
    # iteration does not measure it.
    return merge(_sample_result(pooled; nchains=nchains, ndraws=total, npar=npar,
            ndim=npar + neffects, maxdepth=Int(maxdepth),
            stepsize=[ch.eps for ch in chains],
            warmup_divergent=[ch.warmup_divergent for ch in chains], trace=run.trace,
            saved_effects=save_effects,
            effects=_effect_summary(effect_sum, effect_sq, effect_cross, nall, blocks),
            ebfmi=fill(NaN, nchains)),
        (sampler="saem",
            scale_accept=[ch.nscale == 0 ? NaN : ch.scale_accept / ch.nscale for ch in chains],
            ncp_accept=[ch.nncp == 0 ? NaN : ch.ncp_accept / ch.nncp for ch in chains],
            placement=_saem_placement_summary(placement)))
end

"""
    ctsem_saem_sample(laplace, values; kwargs...)

`_ctsem_saem_sample` (documented there) with each full-rank level's effects
standardised in a Cholesky ordered strongest-determined first for the run
(`_with_root_orders`). A caller's own `starts` or `state` are in the natural
order and keep it.
"""
ctsem_saem_sample(laplace::CTSEMLaplaceObjective, values::AbstractVector; kwargs...) =
    _with_sample_roots(() -> _ctsem_saem_sample(laplace, values; kwargs...), laplace,
        values, kwargs, (:starts, :state))

"""What a result reports of an SAEM placement: where it put the sampler, and
how SAEM got there; `nothing` when there was none."""
_saem_placement_summary(p) = p === nothing ? nothing :
    (theta=p.theta, saem_iterations=p.iterations, saem_settled=p.settled,
     saem_trend=p.trend, saem_secs=p.secs)

export ctsem_saem_sample
