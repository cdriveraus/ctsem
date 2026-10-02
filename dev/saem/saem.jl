# SAEM for `intoverpop = 'laplace'` objectives: a research prototype.
#
# The random effects of each unit are sampled; the states stay integrated by
# the filter, as on the Laplace route. The population parameters follow a
# Robbins-Monro recursion on the Fisher-identity gradient
#
#     d/dtheta log p(y | theta) = E[ d/dtheta log p(y, u | theta) | y, theta ],
#
# estimated at the current draws, so the fixed point is the mode of the exact
# marginal posterior -- no Laplace term, so nothing to over-credit where a
# unit's curvature goes singular. See CT-SEM/review/stochopt-saem-design.md
# and CT-SEM/review/LAPLACE-nested-spike-2026-10-02.md.
#
# Reaches engine internals by name through the objective's own module, so it
# is pinned to the engine it was written against (juliaFit f32bb179) and is
# not package code. Include into Main and drive it from R through the bridge.
#
# The E-step is Metropolis-within-Gibbs over each unit's block tree:
#   * every leaf block (a person) takes a random-walk move shaped by its
#     conditional curvature clipped at the prior's -- one filter evaluation per
#     proposal, since only that block's members see it;
#   * every non-leaf block (a study) takes a collapsed move: the block shifts
#     by Delta and each leaf beneath it by its Gaussian linear response
#     A_b Delta, A_b = -D_b^-1 M_ba. The shear has unit Jacobian and the
#     reverse move is -Delta, so the proposal is symmetric; it moves the study
#     and its people together, which a plain Gibbs step on the study block
#     cannot do when they are tightly coupled.
# Two-level trees only (a leaf's ancestors are at most one block).
#
# The M-step preconditions the gradient with the averaged outer product of the
# per-member complete-data scores plus the prior precision. That is the
# complete-data information, larger than the marginal one wherever the random
# effects carry missing information, so the steps are conservative there.

using LinearAlgebra, Random

_saem_ct(o) = parentmodule(typeof(o))

mutable struct SAEMState
    theta::Vector{Float64}
    u::Vector{Vector{Float64}}
    ll::Vector{Vector{Float64}}                    # per unit, per member position
    leafchol::Vector{Vector{Matrix{Float64}}}      # per unit, per block (leaves used)
    rootchol::Vector{Vector{Matrix{Float64}}}      # per unit, per block (non-leaves used)
    response::Vector{Vector{Matrix{Float64}}}      # per unit, per block: A for leaves
    logscale::Vector{Vector{Float64}}
    accrate::Vector{Vector{Float64}}
    bhhh::Matrix{Float64}
    thetabar::Vector{Float64}
    nbar::Int
    iteration::Int
    seed::UInt64
    sbar::Vector{Matrix{Float64}}                  # per unit: averaged member scores
end

_saem_isleaf(blocks, b) = !any(c -> b in c.ancestors, blocks)

# log p(y_i | theta, u) for member position m of unit U; -Inf if it fails.
function _saem_member_ll(o, U, m, theta, Ls, u)
    CT = _saem_ct(o)
    shifted = CT._laplace_member_values(theta, o.spec, Ls, u, o.units.offsets[U][m])
    v = try
        o.objective.subject_objectives[o.units.members[U][m]](shifted)
    catch err
        err isa InterruptException && rethrow()
        -Inf
    end
    return isfinite(v) ? Float64(v) : -Inf
end

# Clip a symmetric precision's eigenvalues at `floor` and return the lower
# Cholesky factor of its inverse: a proposal covariance no wider than the
# prior in any direction.
function _saem_proposal_chol(P::AbstractMatrix; floor::Float64=1.0)
    E = eigen(Symmetric(Matrix(P)))
    vals = max.(E.values, floor)
    C = E.vectors * Diagonal(1 ./ vals) * transpose(E.vectors)
    return Matrix(cholesky(Symmetric(C)).L)
end

# Proposal shapes from the unit's curvature at its current draw.
function _saem_refresh_shapes!(st::SAEMState, o, U, Ls)
    CT = _saem_ct(o)
    blocks = o.units.blocks[U]
    u = st.u[U]
    M = CT._laplace_unit_curvature(o, U, st.theta, Ls, u)
    for (b, blk) in enumerate(blocks)
        if _saem_isleaf(blocks, b)
            D = Symmetric(Matrix(M.diag[b]))
            st.leafchol[U][b] = _saem_proposal_chol(D)
            if !isempty(blk.ancestors)
                Dc = Matrix(Symmetric(Matrix(_saem_clip(D))))
                st.response[U][b] = -(Dc \ Matrix(M.coupling[b][1]))
            end
        end
    end
    for (a, blk) in enumerate(blocks)
        _saem_isleaf(blocks, a) && continue
        S = Matrix(M.diag[a])
        for (b, c) in enumerate(blocks)
            isempty(c.ancestors) && continue
            c.ancestors[1] == a || continue
            Bba = Matrix(M.coupling[b][1])
            Dc = Matrix(_saem_clip(Symmetric(Matrix(M.diag[b]))))
            S .-= transpose(Bba) * (Dc \ Bba)
        end
        st.rootchol[U][a] = _saem_proposal_chol(Symmetric((S .+ transpose(S)) ./ 2))
    end
    return st
end

function _saem_clip(P::Symmetric; floor::Float64=1.0)
    E = eigen(P)
    return Symmetric(E.vectors * Diagonal(max.(E.values, floor)) * transpose(E.vectors))
end

"""
    saem_init(o, theta; seed)

State at `theta`: each unit's draw starts at its Laplace mode, and the
proposal shapes come from the curvature there.
"""
function saem_init(o, theta::Vector{Float64}; seed::Integer=1)
    CT = _saem_ct(o)
    CT.ctsem_laplace_evaluate(o, theta; gradient=false)
    nunits = length(o.units.members)
    Ls = CT._laplace_popchols(theta, o.spec)
    npar = length(theta)
    st = SAEMState(copy(theta), [copy(o.modes[U]) for U in 1:nunits],
        [zeros(length(o.units.members[U])) for U in 1:nunits],
        [[zeros(b.size, b.size) for b in o.units.blocks[U]] for U in 1:nunits],
        [[zeros(b.size, b.size) for b in o.units.blocks[U]] for U in 1:nunits],
        [[zeros(b.size, 0) for b in o.units.blocks[U]] for U in 1:nunits],
        [[log(2.38 / sqrt(b.size)) for b in o.units.blocks[U]] for U in 1:nunits],
        [[0.3 for b in o.units.blocks[U]] for U in 1:nunits],
        zeros(npar, npar), zeros(npar), 0, 0, UInt64(seed),
        [zeros(npar, length(o.units.members[U])) for U in 1:nunits])
    for U in 1:nunits
        all(b -> length(b.ancestors) <= 1, o.units.blocks[U]) ||
            error("SAEM prototype handles two-level trees only")
    end
    CT._laplace_parallel(o, 1:nunits) do U
        _saem_refresh_shapes!(st, o, U, Ls)
        for m in eachindex(o.units.members[U])
            st.ll[U][m] = _saem_member_ll(o, U, m, theta, Ls, st.u[U])
        end
        true
    end
    return st
end

# One E-step sweep of unit U: every leaf once, then `nroot` collapsed moves of
# each non-leaf block. Returns acceptance counts.
function _saem_sweep!(st::SAEMState, o, U, Ls, rng, adapt::Float64; nroot::Int=2)
    blocks = o.units.blocks[U]
    u = st.u[U]; ll = st.ll[U]; theta = st.theta
    for (b, blk) in enumerate(blocks)
        _saem_isleaf(blocks, b) || continue
        cols = (blk.offset + 1):(blk.offset + blk.size)
        old = u[cols]
        step = exp(st.logscale[U][b]) .* (st.leafchol[U][b] * randn(rng, blk.size))
        u[cols] .= old .+ step
        newll = [_saem_member_ll(o, U, m, theta, Ls, u) for m in blk.members]
        logr = sum(newll) - sum(ll[blk.members]) - (sum(abs2, u[cols]) - sum(abs2, old)) / 2
        acc = isfinite(logr) && log(rand(rng)) < logr
        if acc
            ll[blk.members] .= newll
        else
            u[cols] .= old
        end
        _saem_adapt!(st, U, b, acc, adapt)
    end
    for (a, blk) in enumerate(blocks)
        _saem_isleaf(blocks, a) && continue
        kids = [b for (b, c) in enumerate(blocks) if !isempty(c.ancestors) && c.ancestors[1] == a]
        acols = (blk.offset + 1):(blk.offset + blk.size)
        for _ in 1:nroot
            old = copy(u)
            delta = exp(st.logscale[U][a]) .* (st.rootchol[U][a] * randn(rng, blk.size))
            u[acols] .+= delta
            for b in kids
                c = blocks[b]
                u[(c.offset + 1):(c.offset + c.size)] .+= st.response[U][b] * delta
            end
            newll = [_saem_member_ll(o, U, m, theta, Ls, u) for m in blk.members]
            logr = sum(newll) - sum(ll[blk.members]) - (sum(abs2, u) - sum(abs2, old)) / 2
            acc = isfinite(logr) && log(rand(rng)) < logr
            if acc
                ll[blk.members] .= newll
            else
                u .= old
            end
            _saem_adapt!(st, U, a, acc, adapt)
        end
    end
    return st
end

function _saem_adapt!(st, U, b, acc, adapt)
    st.accrate[U][b] = 0.95 * st.accrate[U][b] + 0.05 * acc
    adapt > 0 && (st.logscale[U][b] += adapt * (acc - 0.3))
    return nothing
end

# Complete-data score per member at the current draws, the unit's sum, and the
# member scores' outer product summed. The same sweep as
# `_laplace_floored_unit_gradient!`, kept per member.
function _saem_unit_scores!(gU::Vector{Float64}, BU::Matrix{Float64}, SU::Matrix{Float64},
        st, o, U, Ls, dL, positions)
    CT = _saem_ct(o)
    theta = st.theta; u = st.u[U]; npar = length(theta)
    aws = CT._laplace_workspace!(o, Float64, npar)
    grad = zeros(npar); shift = zeros(npar); s = zeros(npar)
    fill!(gU, 0.0); fill!(BU, 0.0)
    for m in eachindex(o.units.members[U])
        i = o.units.members[U][m]
        offsets = o.units.offsets[U][m]
        shifted = CT._laplace_member_values!(shift, theta, o.spec, Ls, u, offsets)
        v = CT._laplace_subject_value_gradient!(grad, o.objective.subject_objectives[i], aws, shifted)
        (isfinite(v) && all(isfinite, grad)) || return false
        s .= grad
        CT._laplace_chol_chain!(s, grad, o.spec, dL, positions, u, offsets)
        gU .+= s
        SU[:, m] .= s
        BLAS.syr!('U', 1.0, s, BU)
    end
    return true
end

"""
    saem_run!(st, o; iterations, burnin, alpha, sweeps, maxstep, refresh, every)

`iterations` SAEM iterations from the state's current one. Step size 1 for
iterations up to `burnin`, then `(k - burnin)^-alpha`; Polyak averaging over
the iterations after `burnin`. Returns a trace.
"""
function saem_run!(st::SAEMState, o; iterations::Int=25, burnin::Int=200,
        alpha::Float64=0.7, sweeps::Int=2, nroot::Int=2, maxstep::Float64=0.05,
        refresh::Int=25, ridge::Float64=1e-6, bhhh_rate::Float64=0.1,
        step_scale::Float64=1.0, verbose::Bool=true, precond::Symbol=:complete,
        gamma_burn::Float64=1.0, sbar_rate::Float64=0.05, average_from::Int=burnin)
    CT = _saem_ct(o)
    nunits = length(o.units.members)
    npar = length(st.theta)
    prec = zeros(npar)
    obj = o.objective
    for k in eachindex(obj.prior_index)
        prec[obj.prior_index[k]] += obj.prior_weight / obj.prior_scale[k]^2
    end
    trace = (iteration=Int[], gamma=Float64[], logpost=Float64[], gnorm=Float64[],
        stepmax=Float64[], acc_leaf=Float64[], acc_root=Float64[], seconds=Float64[],
        theta=Vector{Float64}[])
    gUs = [zeros(npar) for _ in 1:nunits]
    BUs = [zeros(npar, npar) for _ in 1:nunits]
    SUs = [zeros(npar, length(o.units.members[U])) for U in 1:nunits]
    for _ in 1:iterations
        t0 = time()
        st.iteration += 1
        k = st.iteration
        gamma = k <= burnin ? gamma_burn : gamma_burn * (k - burnin)^(-alpha)
        adapt = 1.0 / (1 + k)^0.6
        Ls = CT._laplace_popchols(st.theta, o.spec)
        # Members' log likelihood moves with theta: recompute before sampling.
        ok = CT._laplace_parallel(o, 1:nunits) do U
            rng = Random.Xoshiro(hash((st.seed, k, U)))
            for m in eachindex(o.units.members[U])
                st.ll[U][m] = _saem_member_ll(o, U, m, st.theta, Ls, st.u[U])
            end
            (k % refresh == 0) && _saem_refresh_shapes!(st, o, U, Ls)
            for _ in 1:sweeps
                _saem_sweep!(st, o, U, Ls, rng, adapt; nroot=nroot)
            end
            true
        end
        ok || error("E-step failed at iteration $k")
        dL = CT._laplace_level_chol_derivatives(st.theta, o.spec)
        positions = [CT._laplace_level_positions(o.spec, l) for l in eachindex(o.spec.levels)]
        ok = CT._laplace_parallel(o, 1:nunits) do U
            _saem_unit_scores!(gUs[U], BUs[U], SUs[U], st, o, U, Ls, dL, positions)
        end
        ok || error("M-step scores failed at iteration $k")
        g = sum(gUs)
        CT._ctsem_log_prior_gradient!(g, obj, st.theta)
        if precond === :marginal
            # Fisher identity per member: the averaged member score estimates
            # E[s_i | y], and its outer product the marginal information.
            rho = st.iteration == 1 ? 1.0 : sbar_rate
            fill!(st.bhhh, 0.0)
            for U in 1:nunits
                st.sbar[U] .= (1 - rho) .* st.sbar[U] .+ rho .* SUs[U]
                BLAS.syrk!('U', 'N', 1.0, st.sbar[U], 1.0, st.bhhh)
            end
            st.bhhh .= Matrix(Symmetric(st.bhhh, :U))
        else
            B = Symmetric(sum(BUs), :U)
            rate = k == 1 ? 1.0 : (k <= burnin ? bhhh_rate : max(gamma, bhhh_rate / 10))
            st.bhhh .= (1 - rate) .* st.bhhh .+ rate .* Matrix(B)
        end
        P = Symmetric(st.bhhh) + Diagonal(prec .+ ridge)
        step = gamma * step_scale .* (cholesky(P) \ g)
        smax = maximum(abs, step)
        smax > maxstep && (step .*= maxstep / smax)
        logpost = sum(sum, st.ll) - sum(u -> sum(abs2, u), st.u) / 2 +
            CT._ctsem_log_prior(obj, st.theta)
        st.theta .+= step
        if k > average_from
            st.nbar += 1
            st.thetabar .+= (st.theta .- st.thetabar) ./ st.nbar
        end
        accl = mean_leaf_root(st, o)
        push!(trace.iteration, k); push!(trace.gamma, gamma); push!(trace.logpost, logpost)
        push!(trace.gnorm, norm(g)); push!(trace.stepmax, smax)
        push!(trace.acc_leaf, accl[1]); push!(trace.acc_root, accl[2])
        push!(trace.seconds, time() - t0); push!(trace.theta, copy(st.theta))
    end
    return trace
end

function mean_leaf_root(st, o)
    leaf = Float64[]; root = Float64[]
    for U in eachindex(o.units.members)
        blocks = o.units.blocks[U]
        for b in eachindex(blocks)
            push!(_saem_isleaf(blocks, b) ? leaf : root, st.accrate[U][b])
        end
    end
    return (isempty(leaf) ? NaN : sum(leaf) / length(leaf),
        isempty(root) ? NaN : sum(root) / length(root))
end

# For R: the state's vectors, and a check that the member log likelihoods sum
# to what the Laplace inner objective says at the same u.
saem_theta(st::SAEMState) = copy(st.theta)
saem_thetabar(st::SAEMState) = st.nbar == 0 ? copy(st.theta) : copy(st.thetabar)
saem_units(st::SAEMState) = [copy(u) for u in st.u]
function saem_check(st::SAEMState, o, U::Int)
    CT = _saem_ct(o)
    Ls = CT._laplace_popchols(st.theta, o.spec)
    ref = Ref(0.0)
    CT._laplace_parallel(o, [U]) do U
        aws = CT._laplace_workspace!(o, Float64, length(st.theta))
        r = CT._laplace_unit_objective_gradient(o, U, st.theta, Ls, st.u[U], aws)
        ref[] = r.value + sum(abs2, st.u[U]) / 2
        true
    end
    mine = sum(_saem_member_ll(o, U, m, st.theta, Ls, st.u[U]) for m in eachindex(o.units.members[U]))
    return (engine=ref[], saem=mine)
end

# Globals for driving from R without handing the state across the bridge.
function saem_global_init!(o, theta::Vector{Float64}; seed::Integer=1)
    global __SAEM_O = o
    global __SAEM_ST = saem_init(o, theta; seed=seed)
    return length(__SAEM_ST.u)
end
function saem_global_run!(; precond="complete", kwargs...)
    tr = saem_run!(__SAEM_ST, __SAEM_O; precond=Symbol(precond), kwargs...)
    return (iteration=tr.iteration, gamma=tr.gamma, logpost=tr.logpost, gnorm=tr.gnorm,
        stepmax=tr.stepmax, acc_leaf=tr.acc_leaf, acc_root=tr.acc_root,
        seconds=tr.seconds, theta=reduce(hcat, tr.theta))
end
saem_global_theta() = saem_theta(__SAEM_ST)
saem_global_thetabar() = saem_thetabar(__SAEM_ST)
saem_global_check(U::Integer) = saem_check(__SAEM_ST, __SAEM_O, Int(U))
