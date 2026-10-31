"""
The refined fit's rules placed at every evaluation, and its gradient taken
through the placement.

# Why this exists

The continuation (`laplace_continuation.jl`) holds each flagged unit's rule
fixed in the standardised effects where it was placed, climbs that frozen
objective, and places the rules again between rounds. Its gradient is the
Fisher identity at fixed nodes, which leaves out how the rule's own error
moves with its placement. So the rounds' fixed point is not the maximum of the
objective with the rules placed at the point, the keep rule reads that error
moving as gain or loss, and the density the uncertainty draws on returns a
value placed at the point with a gradient that is not that value's
(`review/SURVEY-quadrature-2026-10-30.md`, sections 4 and 6B).

Here a flagged unit's product rule is placed at its mode and curvature at
every evaluation, as SAS NLMIXED and lme4's `nAGQ` place theirs, and the
gradient is the derivative of that value. One function, smooth wherever the
unit's posterior is one basin, that L-BFGS and the Newton finish climb as they
climb the Laplace objective, and that a sampler reads with a gradient that is
its value's.

# The gradient

A unit's term, its block's effects `u` (`k` of them), its log integrand
`f(u) = sum_m ll_m(theta + L u) - u'u/2` and its mode `c(theta)` with
precision `M = -d2f/du2` there, is

    Q = log sum_j exp(lw_j + f(c + S y_j)) + log|S| - k log(2 pi) / 2

with `y_j` and `lw_j` the unit's fixed point set (`point_sets.jl`: product
Gauss-Hermite or lattice-t), and `S = U^-1` for `M~ = U'U`, `M~` the precision
with its eigenvalues clipped at the prior's (`_quadrature_clipped_scale`). With
`pi_j` the normalised node weights,

    dQ/dtheta = sum_j pi_j df/dtheta(z_j)        Fisher, the frozen rule's
              + a' dc/dtheta                     the centre moving
              + <Mbar, dM/dtheta>                the scale and log|S| moving

where `a = sum_j pi_j grad_u f(z_j)` and `Mbar` is what `tr(dS B) - logdet(M~)/2`
asks of `dM~` (`B = sum_j pi_j y_j grad_u f(z_j)'`), carried back
through the Cholesky factor and the clip (`_placed_mbar`). `dM/dtheta` is
total: `M` moves at a fixed centre and with it, so its centre part joins
`a`, `b_i = <Mbar, dM/du_i>`, and by the implicit function theorem
`dc/dtheta = M^-1 d(grad_u f)/dtheta`. With `v = M^-1 (a + b)`:

  - `v' d(grad_u f)/dtheta`: one seeded sweep at the mode along `L v`, the
    derivative of the gradient along it, and the chain through `L(theta)`.
  - `<Mbar, dM/dtheta>` at a fixed centre, `M = I - L'HL`: the sum over
    `Mbar`'s eigendirections of the third derivative along each, which is the
    trace-jet sweep the Laplace log determinant's gradient runs
    (`_laplace_unit_jet_gradient`), a pass per sign, with `Mbar` in place of
    `M^-1`; and `2 <dL/dtheta, H L Mbar>`, `L`'s own derivative, from the
    seeded sweep's further lanes along `L Mbar`.

None of this reads what the point set is, so every family of rule shares it.
For a Gaussian integrand on a product rule the three placement terms cancel,
as they must: the rule is then exact wherever it sits. A lattice rule is not
exact there, but its error is then the same function of the fixed points at
every placement, so it does not move with the parameters.

# What is placed, and what is not

A flagged unit of one block, on the rule the continuation chose for it
(`o.unit_rules`: product or lattice, and its count, by cost for the accuracy,
`_continuation_refine!`), placed as `_continuation_unit_rule` places it: the
leaf's mode by Newton (`_quadrature_leaf_rule!`), from its maximum at the last
accepted point of an optimisation, or outside one at the last evaluation
(`anchors`), where the continuation searches from the Laplace mode. Within an
optimisation the units are anchored as the Laplace route's are: searched at a
run's start and where L-BFGS stops, followed between, so a trial's value does
not depend on the trials before it. A unit whose search found more than one
maximum follows each and averages its rule placed at each, weighted as the
Laplace route's mode sum weights its terms (`_placed_mode_sum!`), so its term
stays continuous through a fold. A unit of one root block over leaves is placed as the continuation
places it, the root at the joint mode and eliminated precision and each leaf
at its conditional mode at every root node (`_placed_nested_unit!`), which
chains a second implicit derivative; each block on its own set, the same set
for a leaf at every root node. A unit on importance draws keeps its frozen
draws, as does a deeper tree for now: their rules are the continuation's as
placed, with its Fisher gradient. The other units keep their Laplace terms
(`o.rest`), and the prior is theirs.

# Deviations the filter integrates

Where each subject's filter integrates deviations of its own (`CTSEMHybrid`),
the integrand at a node is the members' log likelihood with them integrated
exactly (`_HYBRID_FULL`), as the continuation's is, and the leaf's mode and
precision are that integrand's. So are the sweeps at the mode: seeded and jet
sweeps take that mode (`_tan_laplace_context`, `_tan_jet_context`), the closing
term's log determinant weighted one.
"""

using LinearAlgebra

export ctsem_laplace_placed, ctsem_laplace_placed_evaluate, ctsem_laplace_placed_check,
    ctsem_laplace_placed_hessian, ctsem_laplace_placed_info,
    ctsem_laplace_placed_optimize, ctsem_laplace_placed_quadrature,
    ctsem_laplace_placed_error

"""
    CTSEMLaplacePlaced

The hybrid with its flagged units' rules placed at every evaluation; see the
module docstring. `cont` is the continuation it was built from, which supplies
the flagged set, the Laplace part (`cont.rest`), the frozen rules and each
unit's rule (`cont.unit_rules`). `placed` are the units placed at every
evaluation, `frozen` positions in `cont.rules` of the units kept on their
frozen rules, `sets` each placed unit's point set per block (`point_sets.jl`;
empty for a unit not placed), `se` each placed unit's rule error at the last
evaluation (a lattice rule's shifts' spread, zero for a product rule's), and
`anchors` each placed unit's effects at its maximum where the next
evaluation's Newton starts, with the other maxima it follows (`alt_anchors`,
`_placed_mode_sum!`). Within an optimisation (`running`) they are the maxima
at the last accepted point, moved by the optimiser's per-step callbacks
(`_ctsem_anchor_here!`) from those the last evaluation reached (`modes`,
`alt_modes`); outside one, every evaluation moves them. `live` is how many
maxima carried weight in each unit's term at the last evaluation.
"""
mutable struct CTSEMLaplacePlaced{C} <: CTSEMOptimisable
    cont::C
    placed::Vector{Int}
    frozen::Vector{Int}
    sets::Vector{Vector{CTSEMPointSet}}
    se::Vector{Float64}
    anchors::Vector{Vector{Float64}}
    alt_anchors::Vector{Vector{Vector{Float64}}}
    modes::Vector{Vector{Float64}}
    alt_modes::Vector{Vector{Vector{Float64}}}
    live::Vector{Int}
    running::Bool
    cache::Any
    value_cache::Any
    value_calls::Int
    gradient_calls::Int
    seconds_values::Float64
    seconds_gradients::Float64
    seconds_hessian::Float64
    lbfgs::Union{Nothing,CTSEMLBFGSMemory}
    # Per frozen rule, its effective number of nodes or draws at the last
    # evaluation and at the first after its placement (`_placed_drifted`).
    ess::Vector{Float64}
    ess_placed::Vector{Float64}
    # Per placed unit, its rule's effective points at the last evaluation and
    # at the first after its sets were made (`_placed_rough`).
    unit_ess::Vector{Float64}
    unit_ess_placed::Vector{Float64}
end

"""
    ctsem_laplace_placed(o)

The placed objective from the continuation `o`, its flagged set and rules as
`o` holds them: every flagged unit of one block, or of one root block over
leaves (`_placed_two_level`), whose rule is not draws is placed at every
evaluation, on the point sets of the rule `o` chose for it (`o.unit_rules`,
`_rule_sets`); the rest of the flagged units keep `o`'s frozen rules. Each
placed unit's anchors start at the maxima a search at `o`'s centre finds
(`_laplace_solve_unit_mode!`), its mode first, so at the centre a unit with
one maximum takes `o`'s placed value.
"""
function ctsem_laplace_placed(o::CTSEMLaplaceContinuation)
    laplace = o.laplace
    nunits = length(laplace.units.members)
    placed = Int[]
    frozen = Int[]
    sets = [CTSEMPointSet[] for _ in 1:nunits]
    for (f, U) in enumerate(o.flagged)
        blocks = laplace.units.blocks[U]
        if !o.importance[U] && ((length(blocks) == 1 && blocks[1].size >= 1) ||
                _placed_two_level(laplace, U))
            push!(placed, U)
            sets[U] = _rule_sets(laplace, U, _continuation_choice(o, U))
        else
            push!(frozen, f)
        end
    end
    none() = [Vector{Float64}[] for _ in 1:nunits]
    P = CTSEMLaplacePlaced{typeof(o)}(o, placed, frozen, sets, fill(NaN, nunits),
        [Float64[] for _ in 1:nunits], none(), [Float64[] for _ in 1:nunits], none(),
        zeros(Int, nunits), false, nothing, nothing, 0, 0, 0.0, 0.0, 0.0, nothing,
        fill(NaN, length(frozen)), fill(NaN, length(frozen)), fill(NaN, nunits),
        fill(NaN, nunits))
    _placed_search!(P, copy(o.centre); join=false)
    return P
end

"""
    _placed_search!(P, theta; join=true)

A search at `theta` for each placed unit's maxima (`_laplace_solve_unit_mode!`:
Newton from the origin, the probes, the known modes, and past a maximum near
its fold the one it hands over to), each unit anchored on every maximum it
finds and, with `join`, on every one its anchors held (followed to `theta`
by Newton) that the search missed: what a run's start and its stops take
(`_ctsem_anchor_start!`, `_ctsem_take_search_modes!`), as the Laplace route's
merge does (`_ctsem_merge_modes!`). Each maximum one side held and the other
missed goes into the unit's known modes, and `found`, so that later searches
reach it too. Returns how many units' maxima changed; the evaluation caches
are dropped where any did.
"""
function _placed_search!(P::CTSEMLaplacePlaced, theta::Vector{Float64}; join::Bool=true,
    found::CTSEMFoundModes=CTSEMFoundModes())
    laplace = P.cont.laplace
    Ls = _laplace_popchols(theta, laplace.spec)
    np = length(P.placed)
    sets = Vector{Vector{Vector{Float64}}}(undef, np)
    missed = [Vector{Float64}[] for _ in 1:np]
    changed = falses(np)
    _laplace_parallel(laplace, 1:np) do i
        local U, aws, held, searched, d
        U = P.placed[i]
        d = laplace.units.dims[U]
        aws = _laplace_workspace!(laplace, Float64, length(theta))
        held = Vector{Float64}[]
        if join
            for a in vcat(Vector{Float64}[P.anchors[U]], P.alt_anchors[U])
                length(a) == d || continue
                r = _laplace_newton_unit_mode(laplace, U, theta, Ls, aws, a)
                (r.converged && isfinite(r.value)) || continue
                any(h -> _laplace_maxabsdiff(h, r.u) < 1e-4, held) ||
                    push!(held, collect(Float64, r.u))
            end
        end
        _laplace_solve_unit_mode!(laplace, U, theta, Ls, aws)
        searched = vcat(Vector{Float64}[copy(laplace.modes[U])],
            Vector{Float64}[copy(m) for m in laplace.alt_modes[U]])
        if isempty(held)
            sets[i] = searched
            changed[i] = true
            return true
        end
        new = Vector{Float64}[m for m in searched
            if !any(h -> _laplace_maxabsdiff(h, m) < 1e-4, held)]
        missed[i] = vcat(new, Vector{Float64}[h for h in held
            if !any(m -> _laplace_maxabsdiff(h, m) < 1e-4, searched)])
        sets[i] = vcat(held, new)
        changed[i] = !isempty(new)
        return true
    end
    n = 0
    for (i, U) in enumerate(P.placed)
        for m in missed[i]
            _ctsem_known_add!(laplace, U, m, found)
        end
        changed[i] || continue
        P.anchors[U] = sets[i][1]
        P.alt_anchors[U] = sets[i][2:end]
        n += 1
    end
    P.modes = deepcopy(P.anchors)
    P.alt_modes = deepcopy(P.alt_anchors)
    if n > 0
        P.cache = nothing
        P.value_cache = nothing
    end
    return n
end

"""
    _placed_scale(M, k)

The rule's scale at precision `M` exactly as the leaf rule takes it
(`_quadrature_leaf_rule!`, `_quadrature_clipped_scale`), with what its
derivative needs: `S = U^-1` for `M~ = U'U`, the log scale `-logdet(M~)/2`,
whether the eigenvalues were clipped at one, and if so the eigenvectors and
eigenvalues of the matrix that was clipped.
"""
function _placed_scale(M::Matrix{Float64}, k::Int)
    F = _ctsem_cholesky(copy(M), k)
    D = M
    if issuccess(F)
        Lf = Matrix(F.L)
        D = Lf * transpose(Lf)
        if issuccess(_ctsem_cholesky(D - Matrix{Float64}(LinearAlgebra.I, k, k), k))
            return (S=_ctsem_cholesky_uinv(F), U=Matrix(UpperTriangular(F.U[1:k, 1:k])),
                logdetscale=-logdet(F) / 2, clipped=false, vectors=zeros(0, 0),
                values=Float64[], ok=true)
        end
    end
    E = _ctsem_symeig(D)
    all(isfinite, E.values) || return (S=zeros(0, 0), U=zeros(0, 0), logdetscale=NaN,
        clipped=true, vectors=zeros(0, 0), values=Float64[], ok=false)
    clipped = E.vectors * Diagonal(max.(E.values, 1.0)) * transpose(E.vectors)
    G = _ctsem_cholesky(Matrix(_laplace_symmetrise(clipped)), k)
    issuccess(G) || return (S=zeros(0, 0), U=zeros(0, 0), logdetscale=NaN,
        clipped=true, vectors=zeros(0, 0), values=Float64[], ok=false)
    return (S=_ctsem_cholesky_uinv(G), U=Matrix(UpperTriangular(G.U[1:k, 1:k])),
        logdetscale=-logdet(G) / 2, clipped=true, vectors=Matrix(E.vectors),
        values=collect(E.values), ok=true)
end

"""
    _placed_mbar(sc, B, k)

`Mbar`, the symmetric `k x k` matrix with `<Mbar, dM>` the change in
`tr(dS B) + logdet logdet(S)` for a change `dM` of the unclipped precision
(`logdet` one but under a defensive component, `_placed_leaf!`), `S = U^-1`
and `M~ = U'U` (`_placed_scale`). Through the factor: `dS = -S dU S`, so
`tr(dS B) = <Ubar, dU>` with `Ubar = -triu((S B S)')`, and the upper Cholesky's
reverse, `Mbar = sym(S Phi(Ubar U') S')` with `Phi` the upper triangle and half
the diagonal; `logdet(S) = -logdet(M~)/2` adds `-S S'/2`. Through the clip,
`M~ = V max(Lambda, 1) V'`, by the divided differences of `max(., 1)` on the
eigenvalues (Daleckii-Krein).
"""
function _placed_mbar(sc, B::Matrix{Float64}, k::Int; logdet::Real=1.0)
    S = sc.S
    U = sc.U
    SBS = _placed_mm(_placed_mm(S, B), S)
    Ubar = zeros(k, k)
    @inbounds for j in 1:k, i in 1:j
        Ubar[i, j] = -SBS[j, i]
    end
    X = zeros(k, k)
    _ctsem_mulNT!(X, Ubar, U)
    Phi = zeros(k, k)
    @inbounds for j in 1:k, i in 1:j
        Phi[i, j] = i == j ? X[i, i] / 2 : X[i, j]
    end
    W = zeros(k, k)
    _ctsem_mulNT!(W, _placed_mm(S, Phi), S)
    SS = zeros(k, k)
    _ctsem_mulNT!(SS, S, S)
    Mbar = (W .+ transpose(W)) ./ 2 .- (logdet / 2) .* SS
    sc.clipped || return Mbar
    V = sc.vectors
    lam = sc.values
    G = zeros(k, k)
    _ctsem_mulTN!(G, V, _placed_mm(Mbar, V))
    @inbounds for j in 1:k, i in 1:k
        li, lj = lam[i], lam[j]
        K = if abs(li - lj) > 1e-12 * max(abs(li), abs(lj), 1.0)
            (max(li, 1.0) - max(lj, 1.0)) / (li - lj)
        else
            (li + lj) / 2 > 1.0 ? 1.0 : 0.0
        end
        G[i, j] *= K
    end
    out = zeros(k, k)
    _ctsem_mulNT!(out, _placed_mm(V, G), V)
    return (out .+ transpose(out)) ./ 2
end

"""`A * B` for the rule's small matrices, by the engine's kernel rather than a
BLAS product, which takes OpenBLAS's process-wide lock (`small_linalg.jl`)."""
function _placed_mm(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})
    C = zeros(size(A, 1), size(B, 2))
    _ctsem_mul!(C, A, B)
    return C
end

"""`A * x`, by the engine's kernel."""
function _placed_mv(A::AbstractMatrix{Float64}, x::AbstractVector{Float64})
    y = zeros(size(A, 1))
    _ctsem_mulvec!(y, A, x)
    return y
end

"""Scatter a level-space vector `w` (one entry per effect of level `l`) onto
the raw parameters, into `dir`; `add = true` adds it to what `dir` holds, as a
direction moving effects at two levels that share a raw parameter needs."""
function _placed_scatter!(dir::AbstractVector{Float64}, spec, l::Int, w::AbstractVector;
    add::Bool=false)
    add || fill!(dir, 0.0)
    rho = spec.levels[l].re_index
    @inbounds for p in eachindex(rho)
        dir[rho[p]] += w[p]
    end
    return dir
end

"""
    _placed_scale_factor(F, k)

`_placed_scale` from a factor already formed, exactly as
`_quadrature_clipped_scale` takes it: an outer block's eliminated precision
(`_laplace_factor_repaired!`'s `factors`).
"""
function _placed_scale_factor(F::CTSEMCholesky, k::Int)
    Lf = Matrix(F.L)
    D = Lf * transpose(Lf)
    if issuccess(_ctsem_cholesky(D - Matrix{Float64}(LinearAlgebra.I, k, k), k))
        return (S=_ctsem_cholesky_uinv(F), U=Matrix(UpperTriangular(F.U[1:k, 1:k])),
            logdetscale=-logdet(F) / 2, clipped=false, vectors=zeros(0, 0),
            values=Float64[], ok=true, D=D)
    end
    E = _ctsem_symeig(D)
    clipped = E.vectors * Diagonal(max.(E.values, 1.0)) * transpose(E.vectors)
    G = _ctsem_cholesky(Matrix(_laplace_symmetrise(clipped)), k)
    issuccess(G) || return (S=zeros(0, 0), U=zeros(0, 0), logdetscale=NaN,
        clipped=true, vectors=zeros(0, 0), values=Float64[], ok=false, D=D)
    return (S=_ctsem_cholesky_uinv(G), U=Matrix(UpperTriangular(G.U[1:k, 1:k])),
        logdetscale=-logdet(G) / 2, clipped=true, vectors=Matrix(E.vectors),
        values=collect(E.values), ok=true, D=D)
end

"""The largest gradient a leaf's centre may keep (`_quadrature_leaf_rule!`'s
`gradmax`) for the placement to be differentiated there."""
const _PLACED_STATIONARY = 1e-6

"""
    _placed_leaf!(row, wv, P, U, b, u, start, theta, Ls, dL, positions, want_gradient)

Leaf block `b` of unit `U`'s placed term, its ancestors held where `u` has
them: the conditional mode by Newton from `start` (`_quadrature_leaf_rule!`),
its precision, the block's point set (`sets[b]`, the unit's own by default)
placed there; see the
module docstring. With
`want_gradient`, its gradient in `theta` at those ancestors into `row`, and
into `wv` the raw-space vector whose loadings give its derivative in an
ancestor's effects (`L_a' wv[rho_a]`): the node-weighted mean of the members'
raw gradients, plus the gradient's derivative along the mode's direction `L v`,
less the jets' `grad <Mbar, L'HL>` -- the same three terms as the `theta`
gradient, without the chain through `L(theta)`, since an ancestor's effects
reach the members only through their raw values. `u`'s leaf columns are left at
the conditional mode. Returns `(ok, value, se)`, `se` the rule's error there
(`_point_set_error`).
"""
function _placed_leaf!(row::Vector{Float64}, wv::Vector{Float64}, P::CTSEMLaplacePlaced,
    U::Int, b::Int, u::Vector{Float64}, start::Vector{Float64},
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool;
    sets::Vector{CTSEMPointSet}=P.sets[U])
    laplace = P.cont.laplace
    spec = laplace.spec
    units = laplace.units
    block = units.blocks[U][b]
    k = block.size
    l = block.level
    columns = (block.offset + 1):(block.offset + k)
    members = block.members
    nm = length(members)
    npar = length(theta)
    d = units.dims[U]
    failure = (ok=false, value=NaN, se=NaN, ess=NaN)
    aws = _laplace_workspace!(laplace, Float64, npar)
    lr = _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws; start=start)
    # The implicit function theorem is the centre's derivative only at a
    # stationary point; a solve stopped short of one has none to give.
    (lr.ok && lr.gradmax <= _PLACED_STATIONARY) || return failure
    c = lr.centre
    M = Matrix{Float64}(lr.precision)
    sc = _placed_scale(M, k)
    sc.ok || return failure
    ps = sets[b]
    pl = _point_set_place(ps, c, sc.S, sc.U, sc.logdetscale)
    Z = pl.Z
    Yp = pl.Y
    lw = pl.logw
    scratch = _continuation_scratch(laplace, npar)
    peak = -Inf
    acc = 0.0
    acc2 = 0.0
    racc = 0.0
    repacc = zeros(ps.nrep)
    gacc = want_gradient ? zeros(npar) : Float64[]
    vacc = want_gradient ? zeros(npar) : Float64[]
    aacc = zeros(k)
    Bacc = zeros(k, k)
    gh = zeros(npar)
    gv = zeros(npar)
    gu = zeros(d)
    z = zeros(k)
    ey = zeros(k)
    gy = zeros(k)
    for j in 1:size(Z, 2)
        @inbounds for t in 1:k
            z[t] = Z[t, j]
            u[columns[t]] = z[t]
        end
        if want_gradient
            fill!(gh, 0.0)
            fill!(gv, 0.0)
            fill!(gu, 0.0)
        end
        h = 0.0
        ok = true
        for m in members
            ll = _continuation_member!(gh, laplace, U, m, theta, Ls, dL, positions, u,
                aws, scratch, want_gradient)
            if !isfinite(ll)
                ok = false
                break
            end
            h += ll
            if want_gradient
                gv .+= scratch.grad
                _laplace_loading_chain!(gu, scratch.grad, spec, Ls, units.offsets[U][m])
            end
        end
        ok || continue
        t = lw[j] - dot(z, z) / 2 + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                acc *= s
                acc2 *= s * s
                racc *= s
                repacc .*= s
                if want_gradient
                    gacc .*= s
                    vacc .*= s
                    aacc .*= s
                    Bacc .*= s
                end
            end
            peak = t
        end
        w = exp(t - peak)
        acc += w
        acc2 += w * w
        racc += w * pl.r[j]
        ps.nrep > 1 && (repacc[pl.rep[j]] += w)
        if want_gradient
            gacc .+= w .* gh
            vacc .+= w .* gv
            # What the point's term asks of the centre and the scale: a placed
            # point moves with them, `grad f + (1 - r) z`; a prior point is
            # held, and only its weight moves, `r U' grad log rho(y)`.
            if pl.prior[j]
                @inbounds for t in 1:k
                    yv = Yp[t, j]
                    gy[t] = -1.25 * yv / (1 + yv * yv / 4)
                end
                _ctsem_mulTvec!(ey, sc.U, gy)
                ey .*= pl.r[j]
            else
                @inbounds for i in 1:k
                    ey[i] = gu[columns[i]] - pl.r[j] * z[i]
                end
            end
            @inbounds for i in 1:k
                aacc[i] += w * ey[i]
                for r in 1:k
                    Bacc[r, i] += w * Yp[r, j] * ey[i]
                end
            end
        end
    end
    isfinite(peak) || return failure
    # A lattice rule is read as a ratio against its Gaussian
    # (`_point_set_gaussian`); a product rule integrates that exactly.
    ratio = ps.family === :lattice
    gs = ratio ? _point_set_gaussian(pl, sc.U, ps.nrep) : nothing
    value = ratio ? peak + log(acc) - gs.peak - log(gs.sum) + sc.logdetscale :
        peak + log(acc) - k * log(2 * pi) / 2
    se = ratio ? _point_set_error(repacc ./ gs.repsums, ps.nrep) :
        _point_set_error(repacc, ps.nrep)
    ess = acc^2 / acc2
    @inbounds for (t, col) in enumerate(columns)
        u[col] = c[t]
    end
    want_gradient || return (ok=true, value=value, se=se, ess=ess)

    fisher = gacc ./ acc
    a = aacc ./ acc
    B = Bacc ./ acc
    gamma = racc / acc
    if ratio
        # Less the Gaussian's sum's own moves, and Laplace's log|S| beside it.
        aD, BD, gD = _placed_gaussian_moments(gs, pl, k)
        a .-= aD
        B .-= BD
        gamma += 1 - gD
    end
    Mbar = _placed_mbar(sc, B, k; logdet=gamma)
    all(isfinite, Mbar) || return failure
    L = Ls[l]
    offsets = units.offsets[U][first(members)]
    dir = zeros(npar)
    # The integrand's own derivatives: with deviations the filter integrates,
    # their exact integral's, as at the nodes (`_continuation_member!`).
    mode = _hybrid_on(laplace.objective.params.hybrid) ? _HYBRID_FULL : _HYBRID_PROFILE

    # <Mbar, L'HL> by its eigendirections, a jet pass a sign.
    E = _ctsem_symeig(Mbar)
    g = zeros(npar)
    dSnet = zeros(npar)
    have_g = false
    for sgn in (1.0, -1.0)
        picked = [r for r in 1:k if sgn * E.values[r] > 0]
        isempty(picked) && continue
        for first_ in 1:_CTSEM_JET_WIDTH:length(picked)
            chunk = picked[first_:min(first_ + _CTSEM_JET_WIDTH - 1, length(picked))]
            D = zeros(npar, length(chunk))
            for (col, r) in enumerate(chunk)
                _placed_scatter!(dir, spec, l, _placed_mv(L,
                    sqrt(abs(E.values[r])) .* E.vectors[:, r]))
                D[:, col] .= dir
            end
            pass = _laplace_unit_jet_gradient(laplace, U, theta, Ls, u, members, D,
                length(chunk); mode=mode)
            pass.ok || return failure
            (all(isfinite, pass.dS) && all(isfinite, pass.d0)) || return failure
            @inbounds for m in 1:nm, t in 1:npar
                dSnet[t] += sgn * pass.dS[t, m]
                have_g || (g[t] += pass.d0[t, m])
            end
            have_g = true
        end
    end

    # b = <Mbar, dM/du> and the mode's direction v = M^-1 (a + b).
    bfull = zeros(d)
    _laplace_loading_chain!(bfull, dSnet, spec, Ls, offsets)
    rhs = Float64[a[i] - bfull[columns[i]] for i in 1:k]
    FM = _ctsem_cholesky(copy(M), k)
    v = if issuccess(FM)
        FM \ rhs
    else
        EM = _ctsem_symeig(M)
        all(x -> abs(x) > 1e-12, EM.values) || return failure
        _placed_mv(EM.vectors, _placed_mv(Matrix(transpose(EM.vectors)), rhs) ./
            EM.values)
    end
    all(isfinite, v) || return failure

    # H along L v and along L Mbar's columns, from seeded sweeps at the mode.
    ndir = 1 + k
    Dm = zeros(npar, ndir)
    _placed_scatter!(dir, spec, l, _placed_mv(L, v))
    Dm[:, 1] .= dir
    for q in 1:k
        _placed_scatter!(dir, spec, l, _placed_mv(L, Mbar[:, q]))
        Dm[:, 1 + q] .= dir
    end
    HD = _placed_sweep!(g, have_g, laplace, U, theta, Ls, u, members, Dm, mode)
    HD === nothing && return failure

    # The assembly: Fisher, the centre moving, the scale moving.
    vfull = zeros(d)
    @inbounds for (t, col) in enumerate(columns)
        vfull[col] = v[t]
    end
    @inbounds for t in 1:npar
        row[t] = fisher[t] + HD[t, 1] - dSnet[t]
        wv[t] = vacc[t] / acc + HD[t, 1] - dSnet[t]
    end
    Hd = HD[:, 1]
    _laplace_chol_chain!(row, Hd, spec, dL, positions, u, offsets)
    _laplace_chol_chain!(row, g, spec, dL, positions, vfull, offsets)
    _laplace_chol_chain!(row, .-dSnet, spec, dL, positions, u, offsets)
    eq = zeros(d)
    for q in 1:k
        fill!(eq, 0.0)
        eq[columns[q]] = 1.0
        _laplace_chol_chain!(row, -2 .* HD[:, 1 + q], spec, dL, positions, eq, offsets)
    end
    (all(isfinite, row) && all(isfinite, wv)) || return failure
    return (ok=true, value=value, se=se, ess=ess)
end

"""
    _placed_gaussian_moments(gs, pl, k)

What the Gaussian's sum (`_point_set_gaussian`) asks of the centre, the scale
and `log|S|`, normalised over its own weights: `sum pi_j e_j`, `sum pi_j y_j
e_j'` and `sum pi_j r_j`, as the integrand's sum's `a`, `B` and `gamma` are
taken in `_placed_leaf!`.
"""
function _placed_gaussian_moments(gs, pl, k::Int)
    N = length(gs.t)
    a = zeros(k)
    B = zeros(k, k)
    g = 0.0
    for j in 1:N
        w = exp(gs.t[j] - gs.peak) / gs.sum
        w == 0 && continue
        g += w * pl.r[j]
        @inbounds for i in 1:k
            ei = w * gs.e[i, j]
            a[i] += ei
            for r in 1:k
                B[r, i] += pl.Y[r, j] * ei
            end
        end
    end
    return a, B, g
end

"""
    _placed_sweep!(g, have_g, laplace, U, theta, Ls, u, members, D, mode)

The members' gradients' derivatives along `D`'s columns at `u`, summed over
the members (`npar x size(D, 2)`), by seeded sweeps as wide as the curvature's
lanes; the members' summed gradient into `g` where `have_g` is false.
`nothing` where a sweep fails.
"""
function _placed_sweep!(g::Vector{Float64}, have_g::Bool, laplace, U::Int,
    theta::Vector{Float64}, Ls, u::Vector{Float64}, members, D::Matrix{Float64},
    mode::Int)
    npar, ndir = size(D)
    nm = length(members)
    lanes = _laplace_seeded_lanes(laplace)
    HD = zeros(npar, ndir)
    for first_ in 1:lanes:ndir
        stop = min(first_ + lanes - 1, ndir)
        nd = stop - first_ + 1
        pass = _laplace_unit_seeded_gradient(laplace, U, theta, Ls, u, members,
            Matrix{Float64}(D[:, first_:stop]), nd; mode=mode)
        pass.ok || return nothing
        all(isfinite, pass.d1c) || return nothing
        @inbounds for w in 1:nd, m in 1:nm, t in 1:npar
            HD[t, first_ + w - 1] += pass.d1c[t, (w - 1) * nm + m]
        end
        if !have_g
            @inbounds for m in 1:nm, t in 1:npar
                g[t] += pass.d0[t, m]
            end
            have_g = true
        end
    end
    return HD
end

"""
    _placed_unit!(row, P, U, theta, Ls, dL, positions, want_gradient)

Unit `U`'s placed term at `theta`, and with `want_gradient` its gradient into
`row` (unweighted). A unit with one maximum is its rule placed there
(`_placed_one!`), from its anchor; one whose anchors hold more
(`P.alt_anchors`) follows each from its own and averages their rules
(`_placed_mode_sum!`), until all but one have gone past a fold. Returns `(ok,
value, se, ess, modes, live)`: `se` the rule's error, `ess` its effective
points, `modes` the maxima to follow from here (the one of largest weight
first) and `live` how many carried weight.
"""
function _placed_unit!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool;
    sets::Vector{CTSEMPointSet}=P.sets[U])
    start = P.anchors[U]
    if _LAPLACE_MODE_SUM[] && !isempty(P.alt_anchors[U])
        laplace = P.cont.laplace
        aws = _laplace_workspace!(laplace, Float64, length(theta))
        held = _laplace_newton_unit_mode(laplace, U, theta, Ls, aws, start)
        kept = _laplace_follow_modes(laplace, U, theta, Ls, aws, held, start,
            P.alt_anchors[U])
        length(kept) >= 2 && return _placed_mode_sum!(row, P, U, theta, Ls, dL,
            positions, want_gradient, sets, kept, aws)
        # One maximum left, or none: that one, or the anchor's, on its own.
        isempty(kept) || (start = collect(Float64, kept[1].u))
    end
    r = _placed_one!(row, P, U, start, theta, Ls, dL, positions, want_gradient, sets)
    return (ok=r.ok, value=r.value, se=r.se, ess=r.ess,
        modes=Vector{Float64}[r.anchor], live=1)
end

"""
    _placed_one!(row, P, U, start, theta, Ls, dL, positions, want_gradient, sets)

Unit `U`'s rule placed at the maximum Newton reaches from `start`: a unit of
one block is one leaf (`_placed_leaf!`), a unit of a root over leaves an outer
rule over them (`_placed_nested_unit!`). Returns `(ok, value, se, ess,
anchor)`, `anchor` the unit's effects at that maximum.
"""
function _placed_one!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    start::Vector{Float64}, theta::Vector{Float64}, Ls, dL, positions,
    want_gradient::Bool, sets::Vector{CTSEMPointSet})
    blocks = P.cont.laplace.units.blocks[U]
    length(blocks) == 1 || return _placed_nested_unit!(row, P, U, theta, Ls, dL,
        positions, want_gradient; sets=sets, start=start)
    block = blocks[1]
    u = copy(start)
    r = _placed_leaf!(row, zeros(length(theta)), P, U, 1, u,
        Float64[u[c] for c in (block.offset + 1):(block.offset + block.size)], theta,
        Ls, dL, positions, want_gradient; sets=sets)
    return (ok=r.ok, value=r.value, se=r.se, ess=r.ess, anchor=r.ok ? u : start)
end

"""
    _placed_mode_sum!(row, P, U, theta, Ls, dL, positions, want_gradient, sets,
        kept, aws)

Unit `U`'s term where it follows more than one maximum (`kept`, Newton results
from its anchors, `_laplace_follow_modes`): its rule placed at each, averaged
by the weights the Laplace route's mode sum gives their Laplace terms
(`_laplace_mode_sum`),

    Q = sum_k pi_k Q_k,   pi_k proportional to f_k exp(t_k),

`Q_k` the unit's rule placed at maximum `k` (`_placed_one!`: its own centre,
clipped curvature, set and guards), `t_k` its Laplace term and `f_k` its fades
in the curvature's smallest eigenvalue and in the valley's depth. Each `Q_k`
estimates the unit's whole integral, so their average estimates it too, and
no worse than the worst of them; no masses are summed, so a fit has no reason
to grow maxima, the reason sums of masses failed (gD1 s7). A maximum's weight
fades to zero before its fold, so the term is continuous wherever the mode
sum is. Only the maxima carrying weight are placed. The gradient is

    sum_k pi_k grad Q_k + sum_k (Q_k - Q) grad pi_k,

the second the mode sum's own gradient through its terms and fades
(`_laplace_mode_sum_gradient!`) with `pi_k (Q_k - Q)` where it takes
`pi_k (t_k - T)` (`_laplace_mode_seeds`): small, the rules' disagreement, but
what keeps the gradient the value's. The error is the larger of the weighted
rules' errors and the rules' weighted spread about `Q`.
"""
function _placed_mode_sum!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool,
    sets::Vector{CTSEMPointSet}, kept, aws)
    laplace = P.cont.laplace
    npar = length(theta)
    failure = (ok=false, value=NaN, se=NaN, ess=NaN,
        modes=Vector{Float64}[P.anchors[U]], live=0)
    laplace.modes[U] = collect(Float64, kept[1].u)
    laplace.alt_modes[U] = Vector{Float64}[collect(Float64, m.u) for m in kept[2:end]]
    primal = _laplace_unit_primal(laplace, U, theta, Ls, laplace.modes[U], aws)
    isfinite(primal.term) || return failure
    s = _laplace_mode_sum(laplace, U, theta, Ls, aws, primal)
    # The maxima to follow from here, as the sum recorded them: the one of
    # largest weight, then the rest within its span.
    modes = vcat(Vector{Float64}[copy(laplace.modes[U])],
        Vector{Float64}[copy(m) for m in laplace.alt_modes[U]])
    live = Int[k for k in eachindex(s.pi) if s.pi[k] > 0]
    K = length(live)
    rows = [want_gradient ? zeros(npar) : Float64[] for _ in 1:K]
    Q = zeros(K)
    se = zeros(K)
    ess = zeros(K)
    for (i, k) in enumerate(live)
        r = _placed_one!(rows[i], P, U, collect(Float64, s.modes[k]), theta, Ls, dL,
            positions, want_gradient, sets)
        r.ok || return failure
        Q[i] = r.value
        se[i] = isfinite(r.se) ? r.se : 0.0
        ess[i] = r.ess
    end
    p = s.pi[live]
    value = sum(p .* Q)
    err = max(sum(p .* se), sqrt(sum(p .* (Q .- value) .^ 2)))
    lead = argmax(p)
    if want_gradient
        fill!(row, 0.0)
        for i in 1:K
            row .+= p[i] .* rows[i]
        end
        if K > 1
            c = zeros(length(s.pi))
            for (i, k) in enumerate(live)
                c[k] = p[i] * (Q[i] - value)
            end
            seeds = _laplace_mode_seeds(s.fades, c)
            weighs = merge(s, (route=c, gbar=seeds.gbar, Mbar=seeds.Mbar, pts=seeds.pts))
            part = zeros(npar)
            _laplace_mode_sum_gradient!(part, laplace, U, theta, Ls, dL,
                _laplace_level_chol_solves(laplace.spec, Ls, dL), weighs) || return failure
            row .+= part
        end
        all(isfinite, row) || return failure
    end
    return (ok=true, value=value, se=err, ess=ess[lead], modes=modes, live=K)
end

"""Whether unit `U`'s block tree is one root over leaves, the shape
`_placed_nested_unit!` places."""
function _placed_two_level(laplace::CTSEMLaplaceObjective, U::Integer)
    blocks = laplace.units.blocks[U]
    length(blocks) > 1 || return false
    tree = _quadrature_children(blocks)
    length(tree.roots) == 1 || return false
    root = tree.roots[1]
    blocks[root].size >= 1 || return false
    return all(c -> isempty(tree.children[c]) && blocks[c].size >= 1 &&
        blocks[c].ancestors == [root], tree.children[root]) &&
        length(tree.children[root]) == length(blocks) - 1
end

"""
    _placed_nested_unit!(row, P, U, theta, Ls, dL, positions, want_gradient)

A unit of one root block over leaves, placed as the continuation places it
(`_continuation_block_rule`): the root's rule at the unit's joint mode, scaled
by the root's eliminated precision `D0 = M00 - sum_c M0c Mcc^-1 Mc0`, and at
each of its nodes each leaf's rule at its conditional mode
(`_placed_leaf!`). The root's term is the single-block one with the leaves'
summed placed terms for its members' log likelihood, so its gradient is the
same three parts: the node-weighted mean of the leaves' gradients at fixed
root effects (Fisher), the joint mode moving, and `D0` moving. The joint mode
moves by `M^-1 dgrad/dtheta` over the whole unit; `<Mbar0, dD0> = <W, dM>` with
`W = E Mbar0 E'`, `E = [I; -Mcc^-1 Mc0]`, so its jets run along `E`'s columns
-- for a leaf's members, the root's effects and that leaf's together -- and
`v = M^-1 (a0 + b_W)` is one block solve. A leaf's derivative in the root's
effects is `L0' wv[rho0]` (`_placed_leaf!`). Its error is the root rule's and,
in quadrature, the leaves' at the root's nodes, each weighted as its node is:
a leaf's set is the same at every root node, so its errors there are not
independent and do not average out. The joint mode's Newton starts from
`start`, the unit's anchor by default.
"""
function _placed_nested_unit!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool;
    sets::Vector{CTSEMPointSet}=P.sets[U], start::Vector{Float64}=P.anchors[U])
    laplace = P.cont.laplace
    spec = laplace.spec
    units = laplace.units
    blocks = units.blocks[U]
    tree = _quadrature_children(blocks)
    r0 = tree.roots[1]
    kids = tree.children[r0]
    root = blocks[r0]
    k0 = root.size
    l0 = root.level
    cols0 = (root.offset + 1):(root.offset + k0)
    npar = length(theta)
    d = units.dims[U]
    failure = (ok=false, value=NaN, se=NaN, ess=NaN, anchor=start)
    aws = _laplace_workspace!(laplace, Float64, npar)
    jm = _laplace_newton_unit_mode(laplace, U, theta, Ls, aws, copy(start))
    (jm.converged && isfinite(jm.value)) || return failure
    cj = collect(Float64, jm.u)
    M = _laplace_copy_blockmatrix(_laplace_unit_curvature(laplace, U, theta, Ls, cj))
    fz = _laplace_factor_repaired!(_laplace_copy_blockmatrix(M), blocks)
    (fz.ok && !fz.repaired) || return failure
    sc = _placed_scale_factor(fz.factors[r0], k0)
    sc.ok || return failure
    ps = sets[r0]
    pl = _point_set_place(ps, Float64[cj[c] for c in cols0], sc.S, sc.U, sc.logdetscale)
    Z = pl.Z
    Yp = pl.Y
    lw = pl.logw
    peak = -Inf
    acc = 0.0
    acc2 = 0.0
    racc = 0.0
    repacc = zeros(ps.nrep)
    leafse = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    aacc = zeros(k0)
    Bacc = zeros(k0, k0)
    rowc = zeros(npar)
    wvc = zeros(npar)
    G = zeros(k0)
    gd = zeros(d)
    z = zeros(k0)
    ey = zeros(k0)
    gy = zeros(k0)
    rowsum = zeros(npar)
    for j in 1:size(Z, 2)
        u = copy(cj)
        @inbounds for t in 1:k0
            z[t] = Z[t, j]
            u[cols0[t]] = z[t]
        end
        want_gradient && (fill!(rowsum, 0.0); fill!(G, 0.0))
        h = 0.0
        hse = 0.0
        ok = true
        for cb in kids
            ccols = (blocks[cb].offset + 1):(blocks[cb].offset + blocks[cb].size)
            leaf = _placed_leaf!(rowc, wvc, P, U, cb, u, Float64[cj[c] for c in ccols],
                theta, Ls, dL, positions, want_gradient; sets=sets)
            if !leaf.ok
                ok = false
                break
            end
            h += leaf.value
            hse += leaf.se^2
            if want_gradient
                rowsum .+= rowc
                fill!(gd, 0.0)
                _laplace_loading_chain!(gd, wvc, spec, Ls,
                    units.offsets[U][first(blocks[cb].members)])
                @inbounds for t in 1:k0
                    G[t] += gd[cols0[t]]
                end
            end
        end
        # A leaf that could not be scored fails the rule, as in
        # `_continuation_node_sum`, rather than dropping its node.
        ok || return failure
        t = lw[j] - dot(z, z) / 2 + h
        if t > peak
            if isfinite(peak)
                s = exp(peak - t)
                acc *= s
                acc2 *= s * s
                racc *= s
                repacc .*= s
                leafse *= s
                if want_gradient
                    gacc .*= s
                    aacc .*= s
                    Bacc .*= s
                end
            end
            peak = t
        end
        w = exp(t - peak)
        acc += w
        acc2 += w * w
        racc += w * pl.r[j]
        leafse += w * sqrt(hse)
        ps.nrep > 1 && (repacc[pl.rep[j]] += w)
        if want_gradient
            gacc .+= w .* rowsum
            # As at a leaf (`_placed_leaf!`), with the leaves' gradient in the
            # root's effects for the members'.
            if pl.prior[j]
                @inbounds for t in 1:k0
                    yv = Yp[t, j]
                    gy[t] = -1.25 * yv / (1 + yv * yv / 4)
                end
                _ctsem_mulTvec!(ey, sc.U, gy)
                ey .*= pl.r[j]
            else
                @inbounds for i in 1:k0
                    ey[i] = G[i] - pl.r[j] * z[i]
                end
            end
            @inbounds for i in 1:k0
                aacc[i] += w * ey[i]
                for r in 1:k0
                    Bacc[r, i] += w * Yp[r, j] * ey[i]
                end
            end
        end
    end
    isfinite(peak) || return failure
    ratio = ps.family === :lattice
    gs = ratio ? _point_set_gaussian(pl, sc.U, ps.nrep) : nothing
    value = ratio ? peak + log(acc) - gs.peak - log(gs.sum) + sc.logdetscale :
        peak + log(acc) - k0 * log(2 * pi) / 2
    se = sqrt((ratio ? _point_set_error(repacc ./ gs.repsums, ps.nrep) :
        _point_set_error(repacc, ps.nrep))^2 + (leafse / acc)^2)
    ess = acc^2 / acc2
    want_gradient || return (ok=true, value=value, se=se, ess=ess, anchor=cj)

    fisher = gacc ./ acc
    a0 = aacc ./ acc
    B0 = Bacc ./ acc
    gamma0 = racc / acc
    if ratio
        aD, BD, gD = _placed_gaussian_moments(gs, pl, k0)
        a0 .-= aD
        B0 .-= BD
        gamma0 += 1 - gD
    end
    Mbar0 = _placed_mbar(sc, B0, k0; logdet=gamma0)
    all(isfinite, Mbar0) || return failure
    mode = _hybrid_on(laplace.objective.params.hybrid) ? _HYBRID_FULL : _HYBRID_PROFILE
    E = _ctsem_symeig(Mbar0)
    # E's leaf blocks, -Mcc^-1 Mc0, one a leaf.
    Ec = Dict{Int,Matrix{Float64}}()
    for cb in kids
        Fc = _ctsem_cholesky(copy(M.diag[cb]), blocks[cb].size)
        issuccess(Fc) || return failure
        Ec[cb] = .-(Fc \ M.coupling[cb][1])
    end
    dir = zeros(npar)
    leafdir! = function (dirv, cb, y)
        _placed_scatter!(dirv, spec, l0, _placed_mv(Ls[l0], y))
        _placed_scatter!(dirv, spec, blocks[cb].level,
            _placed_mv(Ls[blocks[cb].level], _placed_mv(Ec[cb], y)); add=true)
        return dirv
    end
    # Per leaf: the jets along E q_r for its members, at the joint mode.
    dS = Dict{Int,Vector{Float64}}()
    gl = Dict{Int,Vector{Float64}}()
    bfull = zeros(d)
    for cb in kids
        members = blocks[cb].members
        nm = length(members)
        dSnet = zeros(npar)
        g = zeros(npar)
        have_g = false
        for sgn in (1.0, -1.0)
            picked = [r for r in 1:k0 if sgn * E.values[r] > 0]
            isempty(picked) && continue
            for first_ in 1:_CTSEM_JET_WIDTH:length(picked)
                chunk = picked[first_:min(first_ + _CTSEM_JET_WIDTH - 1, length(picked))]
                D = zeros(npar, length(chunk))
                for (col, r) in enumerate(chunk)
                    D[:, col] .= leafdir!(dir, cb, sqrt(abs(E.values[r])) .* E.vectors[:, r])
                end
                pass = _laplace_unit_jet_gradient(laplace, U, theta, Ls, cj, members, D,
                    length(chunk); mode=mode)
                pass.ok || return failure
                (all(isfinite, pass.dS) && all(isfinite, pass.d0)) || return failure
                @inbounds for m in 1:nm, t in 1:npar
                    dSnet[t] += sgn * pass.dS[t, m]
                    have_g || (g[t] += pass.d0[t, m])
                end
                have_g = true
            end
        end
        _laplace_loading_chain!(bfull, dSnet, spec, Ls, units.offsets[U][first(members)])
        dS[cb] = dSnet
        gl[cb] = have_g ? g : fill(NaN, npar)
    end
    rhs = .-bfull
    @inbounds for t in 1:k0
        rhs[cols0[t]] += a0[t]
    end
    v = _laplace_block_solve!(zeros(d), fz.factors, fz.coupling, blocks, rhs)
    all(isfinite, v) || return failure

    row .= fisher
    eq = zeros(d)
    for cb in kids
        members = blocks[cb].members
        ccols = (blocks[cb].offset + 1):(blocks[cb].offset + blocks[cb].size)
        kc = length(ccols)
        offsets = units.offsets[U][first(members)]
        # Along J v, and along J W e_q for the root's and this leaf's columns.
        qcols = vcat(collect(cols0), collect(ccols))
        Dm = zeros(npar, 1 + length(qcols))
        _placed_scatter!(dir, spec, l0, _placed_mv(Ls[l0], v[cols0]))
        _placed_scatter!(dir, spec, blocks[cb].level,
            _placed_mv(Ls[blocks[cb].level], v[ccols]); add=true)
        Dm[:, 1] .= dir
        for (i, q) in enumerate(qcols)
            # E' e_q: the root's own coordinate, or this leaf's row of Ec.
            e = i <= k0 ? Float64[t == i ? 1.0 : 0.0 for t in 1:k0] :
                Ec[cb][i - k0, :]
            Dm[:, 1 + i] .= leafdir!(dir, cb, _placed_mv(Mbar0, e))
        end
        g = gl[cb]
        have_g = all(isfinite, g)
        have_g || fill!(g, 0.0)
        HD = _placed_sweep!(g, have_g, laplace, U, theta, Ls, cj, members, Dm, mode)
        HD === nothing && return failure
        dSnet = dS[cb]
        @inbounds for t in 1:npar
            row[t] += HD[t, 1] - dSnet[t]
        end
        _laplace_chol_chain!(row, HD[:, 1], spec, dL, positions, cj, offsets)
        _laplace_chol_chain!(row, g, spec, dL, positions, v, offsets)
        _laplace_chol_chain!(row, .-dSnet, spec, dL, positions, cj, offsets)
        for (i, q) in enumerate(qcols)
            fill!(eq, 0.0)
            eq[q] = 1.0
            _laplace_chol_chain!(row, -2 .* HD[:, 1 + i], spec, dL, positions, eq, offsets)
        end
    end
    all(isfinite, row) || return failure
    return (ok=true, value=value, se=se, ess=ess, anchor=cj)
end

"""
    _placed_frozen(P, theta, Ls, dL, positions, want_gradient)

The frozen rules' terms (`P.frozen`) at `theta`, and their weighted gradient:
the continuation's rules as placed, by the Fisher identity
(`_continuation_rule_value`), a unit a task.
"""
function _placed_frozen(P::CTSEMLaplacePlaced, theta::Vector{Float64}, Ls, dL,
    positions, want_gradient::Bool)
    o = P.cont
    laplace = o.laplace
    npar = length(theta)
    nf = length(P.frozen)
    values = fill(NaN, nf)
    rows = [want_gradient ? zeros(npar) : Float64[] for _ in 1:nf]
    nf == 0 && return (values=values, rows=rows, ok=true)
    ok = _laplace_parallel(laplace, 1:nf) do i
        local rule, U, aws, scratch, u, total, r, ess
        rule = o.rules[P.frozen[i]]
        U = rule.unit
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        total = 0.0
        ess = Inf
        for root in rule.roots
            r = _continuation_rule_value(o, U, theta, Ls, dL, positions, root, u, aws,
                scratch, want_gradient)
            isfinite(r.value) || return false
            total += r.value
            ess = min(ess, r.ess)
            want_gradient && (rows[i] .+= r.gradient)
        end
        values[i] = total
        # Each task its own entries.
        P.ess[i] = ess
        isnan(P.ess_placed[i]) && isfinite(ess) && (P.ess_placed[i] = ess)
        return true
    end
    return (values=values, rows=rows, ok=ok)
end

"""
    ctsem_laplace_placed_evaluate(P, values; gradient=true)

The placed objective at `values`: the Laplace terms of the units left to it
and the prior (`P.cont.rest`), every placed unit's rule placed there, and the
frozen rules' terms. `converged` is false where the Laplace part's inner solve
did not converge or a flagged unit could not be scored, as on the continuation.
"""
function ctsem_laplace_placed_evaluate(P::CTSEMLaplacePlaced, values::AbstractVector;
    gradient::Bool=true)
    o = P.cont
    n = length(values)
    theta = _laplace_extend(o.laplace, collect(Float64, values))
    for hit in (P.cache, P.value_cache)
        if hit !== nothing && hit.theta == theta &&
                (!gradient || hit.result.gradient !== nothing)
            # The maxima that evaluation reached, as if it ran again: an
            # accepted step's anchors are the last evaluation's.
            _placed_restore_modes!(P, hit.modes)
            return _laplace_free_result(_continuation_copy(hit.result), n)
        end
    end
    result = _placed_evaluate(P, theta, gradient)
    entry = (theta=theta, result=_continuation_copy(result), modes=_placed_mode_snapshot(P))
    gradient ? (P.cache = entry) : (P.value_cache = entry)
    return _laplace_free_result(result, n)
end

function _placed_evaluate(P::CTSEMLaplacePlaced, theta::Vector{Float64}, gradient::Bool)
    started = time_ns()
    o = P.cont
    laplace = o.laplace
    spec = laplace.spec
    npar = length(theta)
    rest = ctsem_laplace_evaluate(o.rest, theta; gradient=gradient)
    Ls = _laplace_popchols(theta, spec)
    dL = gradient ? _laplace_level_chol_derivatives(theta, spec) : nothing
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    _laplace_ensure_pool!(laplace)
    np = length(P.placed)
    pvalues = fill(NaN, np)
    pse = fill(NaN, np)
    pess = fill(NaN, np)
    plive = zeros(Int, np)
    prows = [gradient ? zeros(npar) : Float64[] for _ in 1:np]
    pmodes = Vector{Vector{Vector{Float64}}}(undef, np)
    placed_ok = _laplace_parallel(laplace, 1:np) do i
        local U, r
        U = P.placed[i]
        r = try
            _placed_unit!(prows[i], P, U, theta, Ls, dL, positions, gradient)
        catch err
            _ctsem_must_propagate(err) && rethrow()
            (ok=false, value=NaN, se=NaN, ess=NaN,
                modes=Vector{Float64}[P.anchors[U]], live=0)
        end
        pvalues[i] = r.value
        pse[i] = r.se
        pess[i] = r.ess
        plive[i] = r.live
        pmodes[i] = r.modes
        return r.ok
    end
    if placed_ok
        for (i, U) in enumerate(P.placed)
            P.modes[U] = pmodes[i][1]
            P.alt_modes[U] = pmodes[i][2:end]
            P.live[U] = plive[i]
            P.se[U] = pse[i]
            P.unit_ess[U] = pess[i]
            isnan(P.unit_ess_placed[U]) && (P.unit_ess_placed[U] = pess[i])
        end
        # Outside an optimisation every evaluation moves the anchors; within
        # one only an accepted step does (`_ctsem_anchor_here!`), so that a
        # trial's value does not depend on the trials before it.
        P.running || _placed_anchor_here!(P)
    end
    frozen = _placed_frozen(P, theta, Ls, dL, positions, gradient)
    ok = placed_ok && frozen.ok
    nunits = length(laplace.units.members)
    unit_loglik = zeros(nunits)
    for (position, U) in enumerate(o.rest_units)
        unit_loglik[U] = rest.unit_loglik[position]
    end
    total = rest.value
    grad = gradient ? copy(rest.gradient) : nothing
    for (i, U) in enumerate(P.placed)
        w = _laplace_unit_weight(laplace, U)
        unit_loglik[U] = pvalues[i]
        total += w * pvalues[i]
        gradient && (grad .+= w .* prows[i])
    end
    for (i, f) in enumerate(P.frozen)
        U = o.rules[f].unit
        w = _laplace_unit_weight(laplace, U)
        unit_loglik[U] = frozen.values[i]
        total += w * frozen.values[i]
        gradient && (grad .+= w .* frozen.rows[i])
    end
    value = ok ? total : NaN
    gradient && !ok && (grad = fill(NaN, npar))
    subject_loglik = copy(rest.subject_loglik)
    for U in Iterators.flatten((P.placed, (o.rules[f].unit for f in P.frozen)))
        members = laplace.units.members[U]
        for i in members
            subject_loglik[i] = unit_loglik[U] / length(members)
        end
    end
    elapsed = (time_ns() - started) / 1e9
    if gradient
        P.gradient_calls += 1
        P.seconds_gradients += elapsed
    else
        P.value_calls += 1
        P.seconds_values += elapsed
    end
    return (value=value, gradient=grad, unit_loglik=unit_loglik,
        subject_loglik=subject_loglik, converged=rest.converged && ok)
end

"""What the last evaluation left besides its value: the maxima it reached, the
placed units' and the Laplace part's (`_laplace_mode_snapshot`), and the
rules' errors and effective sizes, which a run's stop and its end read
(`_placed_rough`, `_placed_drifted`); copied, and put back by
`_placed_restore_modes!` where a cached evaluation is returned. Without them
a run's end read a later trial's effective sizes at the point it returned to,
found nothing to raise there, and handed back a run its stop had ended."""
_placed_mode_snapshot(P::CTSEMLaplacePlaced) = (modes=deepcopy(P.modes),
    alt=deepcopy(P.alt_modes), live=copy(P.live), se=copy(P.se),
    unit_ess=copy(P.unit_ess), ess=copy(P.ess),
    rest=_laplace_mode_snapshot(P.cont.rest))
function _placed_restore_modes!(P::CTSEMLaplacePlaced, s)
    P.modes = deepcopy(s.modes)
    P.alt_modes = deepcopy(s.alt)
    P.live = copy(s.live)
    P.se = copy(s.se)
    P.unit_ess = copy(s.unit_ess)
    length(s.ess) == length(P.ess) && (P.ess = copy(s.ess))
    _laplace_restore_modes!(P.cont.rest, s.rest)
    P.running || _placed_anchor_here!(P)
    return nothing
end

"""The placed units' anchors moved to the maxima the last evaluation reached."""
function _placed_anchor_here!(P::CTSEMLaplacePlaced)
    P.anchors = deepcopy(P.modes)
    P.alt_anchors = deepcopy(P.alt_modes)
    return nothing
end

"""
    ctsem_laplace_placed_quadrature(P, values)

What a refined fit reports at `values`: every unit's quadrature, the flagged
units' by their rules placed there (as `ctsem_laplace_placed_evaluate` scores
them) and every other unit's by the rule the screen would place for it there
(`_continuation_unit_rule`, on the rule it was screened on); the prior; and
each unit's term spread over its members. The objective the
estimate maximises holds the unflagged units' Laplace terms instead, which the
screen found within its tolerance of these.
"""
function ctsem_laplace_placed_quadrature(P::CTSEMLaplacePlaced, values::AbstractVector)
    o = P.cont
    laplace = o.laplace
    theta = _laplace_extend(laplace, collect(Float64, values))
    r = _placed_evaluate(P, theta, false)
    units = copy(r.unit_loglik)
    Ls = _laplace_popchols(theta, laplace.spec)
    opts = _continuation_opts(o)
    _laplace_parallel(laplace, eachindex(o.rest_units)) do i
        local U, placed
        U = o.rest_units[i]
        isempty(laplace.units.blocks[U]) && return true
        placed = try
            _continuation_unit_rule(laplace, U, theta, Ls,
                _laplace_workspace!(laplace, Float64, length(theta)), opts)
        catch err
            _ctsem_must_propagate(err) && rethrow()
            nothing
        end
        placed === nothing || (units[U] = placed.value)
        return true
    end
    value = _laplace_weighted_total(laplace, units) +
        _ctsem_log_prior(laplace.objective, theta)
    members = laplace.units.members
    subjects = zeros(length(laplace.objective.subject_objectives))
    for U in eachindex(members), i in members[U]
        subjects[i] = units[U] / length(members[U])
    end
    return (value=value, unit_loglik=units, subject_loglik=subjects,
        converged=r.converged)
end

"""
    ctsem_laplace_placed_check(P, values; step=1e-5, index=nothing)

The placed objective's gradient at `values` against central differences of its
value, per parameter (`index`, every one by default): the oracle the gradient
is tested against. Each difference starts every unit's Newton from its mode at
`values`, so the two values differ by the placement and not by where a solve
stopped.
"""
function ctsem_laplace_placed_check(P::CTSEMLaplacePlaced, values::AbstractVector;
    step::Real=1e-5, index=nothing)
    x = collect(Float64, values)
    r = ctsem_laplace_placed_evaluate(P, x; gradient=true)
    home = (deepcopy(P.anchors), deepcopy(P.alt_anchors))
    homed!() = (P.anchors = deepcopy(home[1]); P.alt_anchors = deepcopy(home[2]))
    idx = index === nothing ? collect(eachindex(x)) : collect(Int, index)
    fd = fill(NaN, length(idx))
    for (j, t) in enumerate(idx)
        h = step * max(1.0, abs(x[t]))
        xp = copy(x); xp[t] += h
        xm = copy(x); xm[t] -= h
        homed!()
        P.cache = nothing; P.value_cache = nothing
        vp = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xp), false).value
        homed!()
        vm = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xm), false).value
        fd[j] = (vp - vm) / (2 * h)
    end
    homed!()
    P.cache = nothing; P.value_cache = nothing
    return (value=r.value, gradient=r.gradient[idx], differences=fd, index=idx)
end

"""
    _placed_term!(row, P, i, theta, Ls, dL, positions)

Flagged term `i`'s value and unweighted gradient at `theta` into `row`: the
`i`-th placed unit (`_placed_unit!`) from its anchor, or past them the
corresponding frozen rule (`_continuation_rule_value`). What a Hessian column
differences, a term at a time.
"""
function _placed_term!(row::Vector{Float64}, P::CTSEMLaplacePlaced, i::Int,
    theta::Vector{Float64}, Ls, dL, positions)
    np = length(P.placed)
    i <= np && return _placed_unit!(row, P, P.placed[i], theta, Ls, dL, positions,
        true).ok
    o = P.cont
    laplace = o.laplace
    rule = o.rules[P.frozen[i - np]]
    U = rule.unit
    npar = length(theta)
    fill!(row, 0.0)
    aws = _laplace_workspace!(laplace, Float64, npar)
    scratch = _continuation_scratch(laplace, npar)
    u = zeros(laplace.units.dims[U])
    for root in rule.roots
        r = _continuation_rule_value(o, U, theta, Ls, dL, positions, root, u, aws,
            scratch, true)
        isfinite(r.value) || return false
        row .+= r.gradient
    end
    return true
end

"""The flagged terms' units, placed then frozen, as `_placed_term!` numbers them."""
_placed_term_units(P::CTSEMLaplacePlaced) =
    vcat(P.placed, Int[P.cont.rules[f].unit for f in P.frozen])

"""
    ctsem_laplace_placed_hessian(P, values; step=1e-4)

The placed objective's Hessian at `values`: the Laplace part's and the prior's
(`ctsem_laplace_hessian` on `P.cont.rest`, localised to the units that read
each parameter), and the flagged terms' by central differences of their exact
gradients, each column over the terms that read it -- a nonzero gradient at
`values` or at a point moved in every coordinate, since a reader's gradient can
be exactly zero at a point of symmetry (`_laplace_hessian_readers`) -- and
columns whose readers do not overlap moved together
(`_ctsem_local_colours`). Every unit's Newton starts from its mode at `values`.
"""
function ctsem_laplace_placed_hessian(P::CTSEMLaplacePlaced, values::AbstractVector;
    step::Real=1e-4)
    started = time_ns()
    o = P.cont
    laplace = o.laplace
    spec = laplace.spec
    x = collect(Float64, values)
    n = length(x)
    ctsem_laplace_placed_evaluate(P, x; gradient=false)
    H = Matrix{Float64}(ctsem_laplace_hessian(o.rest, x; step=step))
    theta = _laplace_extend(laplace, x)
    npar = length(theta)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    units = _placed_term_units(P)
    nt = length(units)
    np = length(P.placed)
    weights = [_laplace_unit_weight(laplace, U) for U in units]
    rows_at = function (th)
        local Ls, dL, rows, ok
        Ls = _laplace_popchols(th, spec)
        dL = _laplace_level_chol_derivatives(th, spec)
        rows = [zeros(npar) for _ in 1:nt]
        ok = _laplace_parallel(laplace, 1:nt) do i
            _placed_term!(rows[i], P, i, th, Ls, dL, positions)
        end
        return ok ? rows : nothing
    end
    base = rows_at(theta)
    base === nothing && return fill(NaN, n, n)
    moved = copy(theta)
    for j in 1:n
        h = step * max(1.0, abs(x[j])) * (1 + mod(j * 0.6180339887498949, 1.0))
        moved[j] += isodd(j) ? h : -h
    end
    shifted = rows_at(moved)
    reads = [Int[i for i in 1:nt if base[i][t] != 0 ||
        (shifted !== nothing && shifted[i][t] != 0)] for t in 1:n]
    for colour in _ctsem_local_colours(reads, nt)
        hs = [step * max(1.0, abs(x[t])) for t in colour]
        plus = copy(theta)
        minus = copy(theta)
        for (c, t) in enumerate(colour)
            plus[t] += hs[c]
            minus[t] -= hs[c]
        end
        terms = sort!(unique(reduce(vcat, (reads[t] for t in colour); init=Int[])))
        isempty(terms) && continue
        rp = [zeros(npar) for _ in terms]
        rm = [zeros(npar) for _ in terms]
        Lp = _laplace_popchols(plus, spec)
        Lm = _laplace_popchols(minus, spec)
        dLp = _laplace_level_chol_derivatives(plus, spec)
        dLm = _laplace_level_chol_derivatives(minus, spec)
        # A unit following several maxima writes its Laplace records while its
        # rules are averaged (`_placed_mode_sum!`), so its two sides run in one
        # task (side 0); any other term's run apart.
        jobs = Tuple{Int,Int}[]
        for a in eachindex(terms)
            i = terms[a]
            if i <= np && !isempty(P.alt_anchors[P.placed[i]])
                push!(jobs, (a, 0))
            else
                push!(jobs, (a, 1))
                push!(jobs, (a, -1))
            end
        end
        ok = _laplace_parallel(laplace, eachindex(jobs)) do q
            local a, side
            a, side = jobs[q]
            side >= 0 && !_placed_term!(rp[a], P, terms[a], plus, Lp, dLp, positions) &&
                return false
            side <= 0 && !_placed_term!(rm[a], P, terms[a], minus, Lm, dLm, positions) &&
                return false
            return true
        end
        for (c, t) in enumerate(colour)
            if !ok
                H[:, t] .= NaN
                continue
            end
            for (a, i) in enumerate(terms)
                i in reads[t] || continue
                @inbounds for s in 1:n
                    H[s, t] += weights[i] * (rp[a][s] - rm[a][s]) / (2 * hs[c])
                end
            end
        end
    end
    P.seconds_hessian += (time_ns() - started) / 1e9
    return (H .+ transpose(H)) ./ 2
end

"""
    _placed_drifted(P, ratio)

Whether any frozen rule's effective number of nodes or draws has fallen under
`ratio` of what it was at its placement (`_placed_frozen` records both): a
property of the draws, not of a value compared across placements.
"""
_placed_drifted(P::CTSEMLaplacePlaced, ratio::Real) =
    any(i -> isfinite(P.ess[i]) && isfinite(P.ess_placed[i]) &&
        P.ess[i] < ratio * P.ess_placed[i], eachindex(P.frozen))

"""
    _placed_replace!(P, values, ratio)

Place again, at `values`, the frozen rules that have drifted
(`_placed_drifted`), each as its unit keeps it -- its rule, or draws
(`_continuation_unit_rule`). Returns how many were placed.
"""
function _placed_replace!(P::CTSEMLaplacePlaced, values::AbstractVector, ratio::Real)
    o = P.cont
    laplace = o.laplace
    theta = _laplace_extend(laplace, collect(Float64, values))
    Ls = _laplace_popchols(theta, laplace.spec)
    opts = _continuation_opts(o)
    which = [i for i in eachindex(P.frozen) if isfinite(P.ess[i]) &&
        isfinite(P.ess_placed[i]) && P.ess[i] < ratio * P.ess_placed[i]]
    placed = 0
    for i in which
        f = P.frozen[i]
        U = o.rules[f].unit
        uopts = o.importance[U] ? merge(opts, (draw_units=(U,),)) : opts
        r = _continuation_unit_rule(laplace, U, theta, Ls,
            _laplace_workspace!(laplace, Float64, length(theta)), uopts)
        r === nothing && continue
        o.rules[f] = r.rule
        P.ess_placed[i] = NaN
        placed += 1
    end
    P.cache = nothing
    P.value_cache = nothing
    return placed
end

"""
How far past its target a placed rule's error may grow before the run raises
its count (`_placed_rough`): four times `draws_se`, a rule the climb has
taken somewhere it plainly does not serve. A rule's count is chosen where the
continuation was placed (`_continuation_refine!`), and a climb can walk to
where it serves worse: on w6 the 128-point lattices chosen at the Laplace
estimate, 0.18 nats in all there, stood at 1.16 at the estimate the climb
reached, one unit's at 0.85. What decides whether a rule serves is the slope
check at the climb's end (`_placed_slope_check`); raised at twice the target,
w6 restarted its climb eight times and ended unconverged.
"""
const _PLACED_ROUGH = 4.0

"""
The standard deviation, in an estimate's standard errors, that the placed
rules' error may move any estimate by (`_placed_slope_check`): the shift
check's `material` bar.
"""
const _PLACED_SLOPE_BAR = 0.1

"""
    _placed_slope_check(P, values, H; bar=_PLACED_SLOPE_BAR)

How far the placed rules' error can move the estimate at `values`, in its
standard errors, and the units that carry it. What moves an estimate is a
rule's error's slope, not its level, and a lattice rule's value error is
smooth in the parameters but its slope need not be small: with each unit's
128-point lattice within 0.05 nats at the Laplace estimate, the climb on w4s
and w6 walked 1.4 and 1.9 se from a reference fit's estimate, to where the
rules together over-credited the integral by 1.4 and 3.4 nats. Product rules
are read too: on gA14 thirty-nine one-effect units at five nodes, each within
its target, together left the estimate 0.64 se from a reference fit's.

Each lattice unit's gradient on its own points against its gradient on twice
as many, the count step, and each product unit's against the rule two nodes
finer, its node step (`_placed_finer`): the difference sees the rule's
systematic error as well as its noise. Two
independent sets of shifts at the same count (`_rule_sets`'s `alt`), whose
difference over `sqrt(2)` has the variance of either's error, see the noise
alone: read so, w4s's estimate stood 0.63 se from a reference fit where the
check said 0.16, the ratio's small-sample bias shared by both sets. A unit
already at the most points is read on independent shifts.
Carried through the inverse curvature `-H` on its directions with curvature,
as the standard errors are, and divided by them, it is that unit's
contribution to the estimate's error in standard errors. A product rule's
error has one sign at every unit, so the product and mixed units'
contributions add as vectors, signed; the lattice units' are independent and
add in quadrature; per
coordinate the estimate's error from the rules is the first's magnitude plus
the second's root, as the value's error adds them
(`ctsem_laplace_placed_error`). Returns `(sd, units, doublings)`: the largest
coordinate's; the units taken largest contribution first until what is left is
under `bar`; and how many times their points are to double (or a product
rule's nodes to step by two) for the whole to come under `bar`, were their
errors to fall as `n^-0.7`, the slowest rate the lattice showed (at ten
effects; `point_sets.jl`) -- one to five, and checked again after -- and
whether it read every flagged unit (`complete`: none on frozen draws or at a
product rule's most nodes).
"""
function _placed_slope_check(P::CTSEMLaplacePlaced, values::AbstractVector,
    H::AbstractMatrix; bar::Real=_PLACED_SLOPE_BAR)
    o = P.cont
    laplace = o.laplace
    spec = laplace.spec
    x = collect(Float64, values)
    n = length(x)
    units = Int[U for U in P.placed if any(ps -> ps.family === :lattice, P.sets[U]) ||
        _placed_finer(o.unit_rules[U]) !== nothing]
    # Whether every flagged unit's error is read here: not a frozen rule's, nor
    # a product rule's at the most nodes, which has no finer rule to read
    # against.
    complete = isempty(P.frozen) && length(units) == length(P.placed)
    metric = isempty(units) ? nothing : _placed_metric(H, n)
    metric === nothing && return (sd=0.0, units=Int[], doublings=0, complete=false)
    Cinv, se = metric.Cinv, metric.se
    theta = _laplace_extend(laplace, x)
    npar = length(theta)
    Ls = _laplace_popchols(theta, spec)
    dL = _laplace_level_chol_derivatives(theta, spec)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    nu = length(units)
    own = [zeros(npar) for _ in 1:nu]
    other = [zeros(npar) for _ in 1:nu]
    finer = [_placed_finer(o.unit_rules[U]) for U in units]
    alt = [f === nothing ? _rule_sets(laplace, U, o.unit_rules[U]; alt=1) :
        _rule_sets(laplace, U, f) for (U, f) in zip(units, finer)]
    ok = _laplace_parallel(laplace, 1:nu) do i
        local U
        U = units[i]
        return _placed_unit!(own[i], P, U, theta, Ls, dL, positions, true).ok &&
            _placed_unit!(other[i], P, U, theta, Ls, dL, positions, true;
                sets=alt[i]).ok
    end
    ok || return (sd=NaN, units=Int[], doublings=0, complete=false)
    contribution = zeros(nu, n)
    for i in 1:nu
        w = _laplace_unit_weight(laplace, units[i])
        d = ((finer[i] === nothing ? 1 / sqrt(2) : 1.0) * w) .*
            (own[i][1:n] .- other[i][1:n])
        contribution[i, :] .= (Cinv * d) ./ se
    end
    signed = Bool[o.unit_rules[U].family !== :lattice for U in units]
    total = _placed_slope_measure(contribution, signed, trues(nu))
    sd = maximum(total)
    order = sortperm([maximum(abs, contribution[i, :]) for i in 1:nu]; rev=true)
    carrying = Int[]
    kept = trues(nu)
    left = copy(total)
    for i in order
        maximum(left) <= bar && break
        push!(carrying, units[i])
        kept[i] = false
        left = _placed_slope_measure(contribution, signed, kept)
    end
    need = 1.0
    for k in 1:n
        total[k] > bar || continue
        room = bar - left[k]
        share = total[k] - left[k]
        (room > 0 && share > 0) || continue
        need = max(need, log2(share / room) / 0.7)
    end
    return (sd=sd, units=carrying, doublings=isempty(carrying) ? 0 :
        clamp(ceil(Int, need), 1, 5), complete=complete)
end

"""Per coordinate, the estimate's error from the units `kept` of `contribution`'s
rows (`_placed_slope_check`): the `signed` units' rows summed, its magnitude,
plus the root of the others' squares."""
function _placed_slope_measure(contribution::Matrix{Float64}, signed::Vector{Bool},
    kept::AbstractVector{Bool})
    n = size(contribution, 2)
    s = zeros(n)
    r = zeros(n)
    for i in axes(contribution, 1)
        kept[i] || continue
        if signed[i]
            @inbounds for k in 1:n
                s[k] += contribution[i, k]
            end
        else
            @inbounds for k in 1:n
                r[k] += contribution[i, k]^2
            end
        end
    end
    return abs.(s) .+ sqrt.(r)
end

"""The rule a placed unit's error is read against (`_placed_slope_check`) and
the one it is raised to (`_placed_raise!`): a lattice rule with its points
doubled, its count step; a product rule two nodes finer, its node step; a
mixed rule a step along the refine's mixed sequence (`_continuation_steps`),
two nodes finer and its points doubled, each to its most. `nothing` where
nothing can grow. Doubling a mixed rule's points alone left a nested unit
whose blocks all take product rules on the same rule, and the slope check read
its error as zero (gN1, gN3)."""
function _placed_finer(c::CTSEMRuleChoice)
    if c.family === :product
        c.nodes + 2 <= _CONTINUATION_NODES_MAX || return nothing
        return CTSEMRuleChoice(:product, c.nodes + 2, 0)
    end
    if c.family === :mixed
        nodes = min(c.nodes + 2, _CONTINUATION_NODES_MAX)
        points = min(2 * c.points, _LATTICE_MAX)
        (nodes, points) == (c.nodes, c.points) && return nothing
        return CTSEMRuleChoice(:mixed, nodes, points)
    end
    c.points < _LATTICE_MAX || return nothing
    return CTSEMRuleChoice(c.family, c.nodes, 2 * c.points)
end

"""
    _placed_metric(H, n)

The inverse curvature `-H^-1` on the directions with curvature, as the
standard errors are taken, and those standard errors; `nothing` where `H` is
not an `n x n` matrix with curvature.
"""
function _placed_metric(H::AbstractMatrix, n::Integer)
    (size(H) == (n, n) && all(isfinite, H)) || return nothing
    info = -(Matrix{Float64}(H) .+ transpose(Matrix{Float64}(H))) ./ 2
    E = eigen(Symmetric(info))
    top = maximum(E.values)
    top > 0 || return nothing
    keep = E.values .> 1e-8 * top
    V = E.vectors[:, keep]
    Cinv = V * Diagonal(1 ./ E.values[keep]) * transpose(V)
    se = sqrt.(max.(diag(Cinv), 0.0))
    se[se .<= 0] .= Inf
    return (Cinv=Cinv, se=se)
end

"""
    _placed_choose!(P, values, H; bar=_PLACED_SLOPE_BAR)

Each placed unit of one block put on the cheapest rule whose own share of the
estimate's error is under `bar / sqrt(units)`, at `values` on the curvature
`H`: the slope check's measure (`_placed_slope_check`) taken unit by unit and
rule by rule. A product rule's share is its gradient against the rule two
nodes finer, the finer one's measure of it; a lattice rule's, its gradient
against the lattice of twice its points, likewise (at the most points,
against independent shifts, over `sqrt(2)`). A rule must also meet the
value's targets the refine holds it to (`_continuation_refine!`): a product
rule its two node steps' tail within the flagged units' share, a lattice its
spread and count step within `draws_se`; read on one node step alone, the
gradient's share put gA14's heavy-tailed one-effect units back on five nodes
the refine had rejected, two nats off. Rules are
tried in the order of what a pair of them costs (`_rule_cost`), so a unit
takes the product rule wherever its exactness on a near-Gaussian integrand
makes it the cheaper at the accuracy an estimate needs. The continuation
chose by the value's error (`_continuation_refine!`), which the lattice
meets far sooner than its gradient: on w4s every unit took a 128-point
lattice there, and the estimate needed about 2,000 points a unit, where the
625 of the 5-node product rule did as well. A unit no rule takes there keeps
the one with the smallest share. Returns how many units changed rule.
"""
function _placed_choose!(P::CTSEMLaplacePlaced, values::AbstractVector,
    H::AbstractMatrix; bar::Real=_PLACED_SLOPE_BAR)
    o = P.cont
    laplace = o.laplace
    spec = laplace.spec
    x = collect(Float64, values)
    n = length(x)
    units = Int[U for U in P.placed if length(laplace.units.blocks[U]) == 1]
    metric = isempty(units) ? nothing : _placed_metric(H, n)
    metric === nothing && return 0
    budget = Float64(bar) / sqrt(length(units))
    theta = _laplace_extend(laplace, x)
    Ls = _laplace_popchols(theta, spec)
    dL = _laplace_level_chol_derivatives(theta, spec)
    positions = [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)]
    chosen = Vector{CTSEMRuleChoice}(undef, length(units))
    _laplace_parallel(laplace, eachindex(units)) do i
        chosen[i] = _placed_choose_unit(P, units[i], theta, Ls, dL, positions,
            metric, n, budget)
        return true
    end
    rules = copy(o.unit_rules)
    changed = 0
    for (i, U) in enumerate(units)
        chosen[i] == rules[U] && continue
        rules[U] = chosen[i]
        P.sets[U] = _rule_sets(laplace, U, chosen[i])
        P.unit_ess_placed[U] = NaN
        changed += 1
    end
    o.unit_rules = rules
    P.cache = nothing
    P.value_cache = nothing
    return changed
end

"""One unit of `_placed_choose!`: the cheapest rule of its families
(`_continuation_steps`) whose share is within `budget`, by the cost of the
pair of rules that measures it."""
function _placed_choose_unit(P::CTSEMLaplacePlaced, U::Int, theta::Vector{Float64},
    Ls, dL, positions, metric, n::Int, budget::Float64)
    o = P.cont
    laplace = o.laplace
    npar = length(theta)
    w = _laplace_unit_weight(laplace, U)
    families = _continuation_steps(laplace, U, _continuation_opts(o))
    # Each candidate: the rule, the rule that measures it (itself on other
    # shifts for a lattice), a product rule's coarser step (`steps[j - 1]`,
    # for its value's tail), and the pair's cost; their sets built in turn.
    cands = Tuple{CTSEMRuleChoice,CTSEMRuleChoice,CTSEMRuleChoice,Float64}[]
    for steps in families
        for (j, c) in enumerate(steps)
            if c.family === :lattice
                f = _placed_finer(c)
                push!(cands, (c, f === nothing ? c : f, c,
                    _rule_cost(laplace, U, c) + _rule_cost(laplace, U, f === nothing ? c : f)))
            elseif j >= 2 && j < length(steps)
                push!(cands, (c, steps[j + 1], steps[j - 1],
                    _rule_cost(laplace, U, c) + _rule_cost(laplace, U, steps[j + 1])))
            end
        end
    end
    sort!(cands; by=last)
    # The value's targets the refine holds a rule to (`_continuation_refine!`):
    # a rule that misses them is no choice, whatever its gradient's share.
    systematic = o.draws_se / sqrt(max(length(o.flagged), 1))
    best = o.unit_rules[U]
    bestshare = Inf
    a = zeros(npar)
    b = zeros(npar)
    scratch = zeros(npar)
    for (c, m, prev, _) in cands
        lattice = c.family === :lattice
        sa = _rule_sets(laplace, U, c)
        same = lattice && m == c
        sb = same ? _rule_sets(laplace, U, m; alt=1) : _rule_sets(laplace, U, m)
        scale = same ? 1 / sqrt(2) : 1.0
        fill!(a, 0.0)
        fill!(b, 0.0)
        ra = _placed_unit!(a, P, U, theta, Ls, dL, positions, true; sets=sa)
        rb = _placed_unit!(b, P, U, theta, Ls, dL, positions, true; sets=sb)
        (ra.ok && rb.ok) || continue
        if c.family === :product
            rp = _placed_unit!(scratch, P, U, theta, Ls, dL, positions, false;
                sets=_rule_sets(laplace, U, prev))
            rp.ok || continue
            e = _continuation_product_error([rp.value, ra.value, rb.value], 3)
            (e !== nothing && e.coarse <= systematic) || continue
        elseif lattice
            max(ra.se, same ? 0.0 : abs(ra.value - rb.value)) <= o.draws_se || continue
        end
        d = (w * scale) .* (a[1:n] .- b[1:n])
        share = maximum(abs, (metric.Cinv * d) ./ metric.se)
        if share < bestshare
            best = c
            bestshare = share
        end
        share <= budget && return c
    end
    return best
end

"""
    _placed_rough(P, factor)

The placed lattice units whose points can still be doubled and whose rule, at
the last evaluation, has an error (`P.se`, the shifts' spread) past `factor`
times the target, or effective points under half what they were where its
sets were made, as a frozen rule's drift is read (`_placed_drifted`). The
spread is the late signal: on unit 30 of the round-9 stall a placed t
component at the mode was climbed where it was weakest, its effective points
falling from 87 to 1.5 of 1,280 while its spread moved from 0.100 to 0.112
by the time the effective size had halved, and ended at 0.82 against a true
error of 1.53 (`review/MIXTURE-multibasin-measure-2026-10-30.md`, section 3).
"""
_placed_rough(P::CTSEMLaplacePlaced, factor::Real) = Int[U for U in P.placed if
    any(ps -> ps.family === :lattice && ps.count < _LATTICE_MAX, P.sets[U]) &&
    ((isfinite(P.se[U]) && P.se[U] > factor * P.cont.draws_se) ||
     (isfinite(P.unit_ess[U]) && isfinite(P.unit_ess_placed[U]) &&
      P.unit_ess[U] < P.unit_ess_placed[U] / 2))]

"""
    _placed_raise!(P, values, units)

Each of `units`' lattice points doubled, at `values`, until its error there is
within the target or its points are the most a lattice takes
(`_LATTICE_MAX`), or with `times` doubled that many times (to the most); a
product rule's nodes are stepped by two `times` times, once at least, to
`_CONTINUATION_NODES_MAX`. The continuation keeps the new rule
(`o.unit_rules`), so counts only rise and a unit never switches back. Returns
how many were raised.
"""
function _placed_raise!(P::CTSEMLaplacePlaced, values::AbstractVector,
    units::Vector{Int}; times::Integer=0)
    o = P.cont
    laplace = o.laplace
    theta = _laplace_extend(laplace, collect(Float64, values))
    Ls = _laplace_popchols(theta, laplace.spec)
    positions = [_laplace_level_positions(laplace.spec, l)
        for l in eachindex(laplace.spec.levels)]
    rules = copy(o.unit_rules)
    row = zeros(length(theta))
    for U in units
        choice = rules[U]
        # Each step's change in the unit's value, the coarser rule's node or
        # count step, stands as the new rule's systematic error
        # (`o.centre_error`, `ctsem_laplace_placed_error`).
        previous = _placed_unit!(row, P, U, theta, Ls, nothing, positions, false)
        if choice.family === :product || times > 0
            for _ in 1:max(Int(times), 1)
                finer = _placed_finer(choice)
                finer === nothing && break
                choice = finer
            end
            P.sets[U] = _rule_sets(laplace, U, choice)
            r = _placed_unit!(row, P, U, theta, Ls, nothing, positions, false)
            (previous.ok && r.ok) && (o.centre_error[U] = abs(r.value - previous.value))
            r.ok && choice.family !== :product && (P.se[U] = r.se)
            P.unit_ess_placed[U] = NaN
        else
            while choice.points < _LATTICE_MAX
                choice = CTSEMRuleChoice(choice.family, choice.nodes, 2 * choice.points)
                P.sets[U] = _rule_sets(laplace, U, choice)
                P.unit_ess_placed[U] = NaN
                r = _placed_unit!(row, P, U, theta, Ls, nothing, positions, false)
                r.ok || continue
                previous.ok && (o.centre_error[U] = abs(r.value - previous.value))
                P.se[U] = r.se
                max(r.se, previous.ok ? o.centre_error[U] : Inf) <= o.draws_se && break
                previous = r
            end
        end
        rules[U] = choice
    end
    o.unit_rules = rules
    P.cache = nothing
    P.value_cache = nothing
    return length(units)
end

"""
    ctsem_laplace_placed_optimize(P, start; drift=0.5, replacements=8, kwargs...)

`ctsem_optimize` on the placed objective from `start`, its keywords passed on.
Its units are anchored as the Laplace route's are: a search at each run's
start (`_ctsem_anchor_start!`), each unit's Newton from its maxima at the last
accepted point, and a search where L-BFGS stops by its own rule, whose new
maxima join the ones followed (`_ctsem_take_search_modes!`), the Laplace
part's units with them. A unit with more than one maximum averages their
rules (`_placed_mode_sum!`).
The run stops where a frozen rule has drifted (`_placed_drifted`) or a placed
rule's error has grown past `_PLACED_ROUGH` times its target (`_placed_rough`),
and again where it ends with either, since its Newton finish does not consult
the stop; the drifted rules are placed again there (`_placed_replace!`), the
rough ones raised (`_placed_raise!`), and the run goes on from there, with
its L-BFGS memory dropped where frozen rules were placed again; a run its
stop ended goes on even where its end finds nothing to do. Where it ends
on a Hessian, the slope check (`_placed_slope_check`) reads how far the
lattice rules' error can move the estimate; the units carrying more than
`_PLACED_SLOPE_BAR` standard errors have their points raised as many
doublings as it says, and the run goes on with its memory, the objective
having moved by less than the rules' error. At most `replacements` times in
all, the last of them run to its end with no stop. With `hessian`, a curvature
at `start` (the Laplace fit's), each unit of one block is first put on the
cheapest rule whose share of the estimate's error is within the bar
(`_placed_choose!`) and the whole is checked there, so that the run starts
with the rules the estimate needs: checked only where each climb ended, w4s
and w6 restarted eight times. The result is the last run's, with
`replacements` the count, `raised` the units raised, `multimodal` the units
whose term averaged more than one maximum at the end, and `slope_sd` the slope
check's last sd (`NaN` where none ran or it did not read every flagged unit).
"""
function ctsem_laplace_placed_optimize(P::CTSEMLaplacePlaced, start::AbstractVector;
    drift::Real=0.5, replacements::Integer=8, hessian=nothing, kwargs...)
    x = collect(Float64, start)
    count = 0
    raised = 0
    pairs = get(kwargs, :lbfgs_pairs, nothing)
    rest = Dict{Symbol,Any}(k => v for (k, v) in kwargs if k !== :lbfgs_pairs)
    ratio = Float64(drift)
    slope = NaN
    if hessian !== nothing
        # Each unit on the cheapest rule its share of the estimate's error
        # allows, then the whole checked as it will be at the end.
        H0 = Matrix{Float64}(hessian)
        raised += _placed_choose!(P, x, H0)
        ctsem_laplace_placed_evaluate(P, x; gradient=false)
        check = _placed_slope_check(P, x, H0)
        slope = check.complete ? check.sd : NaN
        steep = Int[U for U in check.units if
            _placed_finer(P.cont.unit_rules[U]) !== nothing]
        isempty(steep) || (raised += _placed_raise!(P, x, steep; times=check.doublings))
    end
    try
        while true
            last = count >= replacements
            res = ctsem_optimize(P, x; rest..., lbfgs_pairs=pairs,
                stop=last ? nothing : _PlacedStop(P, ratio))
            x = collect(Float64, res.minimizer)
            # The last evaluation need not have been at the point returned.
            ctsem_laplace_placed_evaluate(P, x; gradient=false)
            drifted = !isempty(P.frozen) && ratio > 0 && _placed_drifted(P, ratio)
            rough = _placed_rough(P, _PLACED_ROUGH)
            # A run its stop ended has not converged, whatever its end finds
            # to do: its end can have moved off the point the stop read.
            stopped = string(get(res, :stop_reason, "")) == "stopped"
            steep = Int[]
            doublings = 0
            if !drifted && isempty(rough) && !stopped
                check = _placed_slope_check(P, x, res.hessian)
                slope = check.complete ? check.sd : NaN
                doublings = check.doublings
                steep = Int[U for U in check.units if
                    _placed_finer(P.cont.unit_rules[U]) !== nothing]
            end
            if last || !(drifted || !isempty(rough) || !isempty(steep) || stopped)
                return merge(res, (replacements=count, raised=raised,
                    multimodal=Base.count(U -> P.live[U] > 1, P.placed), slope_sd=slope))
            end
            drifted && _placed_replace!(P, x, ratio)
            isempty(rough) || (raised += _placed_raise!(P, x, rough))
            isempty(steep) || (raised += _placed_raise!(P, x, steep; times=doublings))
            count += 1
            # Frozen draws placed again change the objective by a rule's
            # error; a raised rule only by its own, which the pairs outlive.
            drifted && (pairs = nothing)
        end
    finally
        # A run an interrupt or an error ended leaves its anchors behind.
        ctsem_laplace_clear_anchors!(P)
    end
end

"""The stop `ctsem_laplace_placed_optimize` hands its run: a frozen rule
drifted, or a placed rule's error past `_PLACED_ROUGH` times its target. A
callable rather than a closure, so that nothing it reads is a captured local."""
struct _PlacedStop{P}
    placed::P
    ratio::Float64
end
(s::_PlacedStop)(st) = (!isempty(s.placed.frozen) && s.ratio > 0 &&
    _placed_drifted(s.placed, s.ratio)) || !isempty(_placed_rough(s.placed, _PLACED_ROUGH))

"""
    ctsem_laplace_placed_error(P)

The refined objective's error against the exact marginal at its last
evaluation, in `ctsem_laplace_continuation_info`'s fields (`error`,
`error_draws`, `error_nodes`, `error_residual`, `centre_error`): each placed
lattice or mixed unit's spread there (`P.se`) in quadrature with the frozen
draws' errors, and summed with them the product rules' node errors, a
lattice unit's count step where it exceeds its spread (`o.centre_error`, from
where its count was set) and what the screen left on Laplace as at the
placement. A lattice rule's spread is free at every evaluation and moves with
where the rule sits; read at the placement alone it said 0.18 nats on w6 where
the estimate the climb reached carried 1.16.
"""
function ctsem_laplace_placed_error(P::CTSEMLaplacePlaced)
    o = P.cont
    err = _continuation_error(o)
    units = copy(o.centre_error)
    random = 0.0
    nodes = 0.0
    for U in P.placed
        isfinite(P.se[U]) || continue
        if any(ps -> ps.family === :lattice, P.sets[U])
            # Its spread where it sits, or its count step where larger: the
            # systematic error the spread cannot see, from where its count was
            # set (`_continuation_count_step`, `_placed_raise!`).
            cs = isfinite(o.centre_error[U]) ? o.centre_error[U] : 0.0
            units[U] = max(P.se[U], cs)
            cs > P.se[U] ? (nodes += cs) : (random += P.se[U]^2)
        elseif P.live[U] > 1
            # A product rule's node error beside its maxima's disagreement
            # (`_placed_mode_sum!`), the larger.
            units[U] = max(P.se[U], isfinite(o.centre_error[U]) ? o.centre_error[U] : 0.0)
            nodes += units[U]
        end
    end
    for U in o.flagged
        o.importance[U] && (random += o.centre_se[U]^2; continue)
        U in P.placed && (any(ps -> ps.family === :lattice, P.sets[U]) ||
            P.live[U] > 1) && continue
        isfinite(o.centre_error[U]) && (nodes += o.centre_error[U])
    end
    draws = sqrt(random)
    return (error=draws + nodes + err.residual, error_draws=draws, error_nodes=nodes,
        error_residual=err.residual, centre_error=units)
end

"""
    ctsem_laplace_placed_info(P)

The placed objective's accounting: units placed and frozen, calls, seconds;
per placed unit its rule's family of the root block (`lattice` 1 where it is
lattice-t, 0 product), its count (`count`: nodes a coordinate, or lattice
points a shift), its points (`points`, member evaluations a value,
`_rule_cost`), its error at the last evaluation (`se`) and how many maxima
its term averaged there (`live`, `_placed_mode_sum!`). A list with no unit
is `[0]`, for the bridge (an empty vector hangs it).
"""
function ctsem_laplace_placed_info(P::CTSEMLaplacePlaced)
    nonempty(v) = isempty(v) ? [0] : v
    laplace = P.cont.laplace
    roots = [_quadrature_children(laplace.units.blocks[U]).roots[1] for U in P.placed]
    lattice = Int[P.sets[U][r].family === :lattice ? 1 : 0
        for (U, r) in zip(P.placed, roots)]
    counts = Int[P.sets[U][r].count for (U, r) in zip(P.placed, roots)]
    points = Float64[_rule_cost(laplace, U, _continuation_choice(P.cont, U))
        for U in P.placed]
    return (placed=nonempty(copy(P.placed)),
        frozen=nonempty([P.cont.rules[f].unit for f in P.frozen]),
        lattice=nonempty(lattice), count=nonempty(counts),
        points=isempty(points) ? [0.0] : points,
        se=isempty(P.placed) ? [0.0] : P.se[P.placed],
        live=isempty(P.placed) ? [0] : P.live[P.placed],
        value_calls=P.value_calls, gradient_calls=P.gradient_calls,
        seconds_values=P.seconds_values, seconds_gradients=P.seconds_gradients,
        seconds_hessian=P.seconds_hessian)
end

# The optimisable protocol (see ctsem_backend.jl), as the continuation takes it.

function ctsem_evaluate(P::CTSEMLaplacePlaced, values::AbstractVector;
    gradient::Bool=true, contributions::Bool=false, gradient_method=:adjoint)
    r = ctsem_laplace_placed_evaluate(P, values; gradient=gradient)
    value = r.converged ? r.value : NaN
    contributions || return (value=value, gradient=r.gradient)
    return (value=value, gradient=r.gradient, subject_loglik=r.subject_loglik)
end

_ctsem_optimise_label(::CTSEMLaplacePlaced) = "refined fit, rules placed"
_ctsem_finish_curvature(::CTSEMLaplacePlaced) = :chord
_ctsem_hessian_cost(::CTSEMLaplacePlaced, npar::Integer) = 2 * Int(npar)
_ctsem_lbfgs_held(P::CTSEMLaplacePlaced) =
    P.lbfgs === nothing ? _ctsem_lbfgs_held(P.cont.laplace) : P.lbfgs
_ctsem_lbfgs_hold!(P::CTSEMLaplacePlaced, memory::CTSEMLBFGSMemory) =
    (P.lbfgs = memory; nothing)
ctsem_hessian(P::CTSEMLaplacePlaced, values::AbstractVector; chunk::Integer=0,
    progress=nothing) = ctsem_laplace_placed_hessian(P, values)
_ctsem_params(P::CTSEMLaplacePlaced) = _ctsem_params(P.cont.laplace)
_ctsem_saturated_for(P::CTSEMLaplacePlaced, minimizer) =
    _laplace_saturated_parameters(P.cont.laplace, minimizer)
_ctsem_optimise_result_extra(P::CTSEMLaplacePlaced, final, log) =
    (placed_units=length(P.placed), frozen_units=length(P.frozen))

# Anchored as the Laplace route's optimisations are (`_ctsem_anchor_start!`):
# the placed units by their own anchors, the Laplace part's units by its.
function _ctsem_anchor_start!(P::CTSEMLaplacePlaced, x::AbstractVector;
    resume_units=nothing, resume_values=nothing)
    _LAPLACE_ANCHORS[] || return nothing
    xv = collect(Float64, x)
    _ctsem_anchor_start!(P.cont.rest, xv)
    _placed_search!(P, _laplace_extend(P.cont.laplace, xv))
    P.running = true
    return nothing
end

function _ctsem_anchor_here!(P::CTSEMLaplacePlaced)
    P.running || return nothing
    _placed_anchor_here!(P)
    _ctsem_anchor_here!(P.cont.rest)
    return nothing
end

_ctsem_anchor_finish!(P::CTSEMLaplacePlaced, x) = ctsem_laplace_clear_anchors!(P)

function ctsem_laplace_clear_anchors!(P::CTSEMLaplacePlaced)
    P.running = false
    ctsem_laplace_clear_anchors!(P.cont.rest)
    return nothing
end

"""
    _ctsem_take_search_modes!(P, x, found)

Where an anchored run on the placed objective stopped by its own rule at `x`:
the Laplace part's units joined with a search there as on the Laplace route
(`_ctsem_take_search_modes!`), and each placed unit's maxima with its own
(`_placed_search!`), so that a maximum born on the run's way, or one past a
fold, joins the ones the unit averages (`_placed_mode_sum!`). Returns the
change this makes to the objective at `x` and how many units' maxima changed,
`NaN` and 0 where none did or the joined objective cannot be evaluated there;
the run carries on from `x` where the change passes `_CTSEM_MODE_RESUME_TOL`.
"""
function _ctsem_take_search_modes!(P::CTSEMLaplacePlaced, x::AbstractVector,
    found::CTSEMFoundModes)
    P.running || return (change=NaN, units=0)
    o = P.cont
    xv = collect(Float64, x)
    before = ctsem_laplace_placed_evaluate(P, xv; gradient=false)
    isfinite(before.value) || return (change=NaN, units=0)
    rest_found = CTSEMFoundModes()
    r = _ctsem_take_search_modes!(o.rest, xv, rest_found)
    for (U, m) in zip(rest_found.units, rest_found.modes)
        push!(found.units, o.rest === o.laplace ? U : o.rest_units[U])
        push!(found.modes, m)
    end
    held = (deepcopy(P.anchors), deepcopy(P.alt_anchors))
    n = _placed_search!(P, _laplace_extend(o.laplace, xv); found=found)
    units = n + (isfinite(r.change) ? r.units : 0)
    units == 0 && return (change=NaN, units=0)
    P.cache = nothing
    P.value_cache = nothing
    after = ctsem_laplace_placed_evaluate(P, xv; gradient=false)
    if !isfinite(after.value)
        P.anchors, P.alt_anchors = held
        P.modes, P.alt_modes = deepcopy(held[1]), deepcopy(held[2])
        P.cache = nothing
        P.value_cache = nothing
        return (change=NaN, units=0)
    end
    _placed_anchor_here!(P)
    return (change=Float64(after.value) - Float64(before.value), units=units)
end

function _ctsem_optimise_trial(P::CTSEMLaplacePlaced, x, want_gradient::Bool,
    gradient_method, limit::Real, log)
    evaluated = _ctsem_try_evaluate(() ->
        ctsem_laplace_placed_evaluate(P, x; gradient=want_gradient))
    valid = evaluated !== nothing && isfinite(evaluated.value) && evaluated.converged &&
        (!want_gradient || _ctsem_gradient_usable(evaluated.gradient, limit))
    return (evaluated=evaluated, valid=valid)
end
