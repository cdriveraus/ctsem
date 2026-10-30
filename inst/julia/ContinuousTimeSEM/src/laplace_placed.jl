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

    Q = log sum_j exp(lw_j + f(c + sqrt(2) S x_j)) + log|S| + k (log 2 - log 2 pi) / 2

with `x_j` the Gauss-Hermite grid, `lw_j` its log weights with `x_j'x_j`
folded in, and `S = U^-1` for `M~ = U'U`, `M~` the precision with its
eigenvalues clipped at the prior's (`_quadrature_clipped_scale`). With `pi_j`
the normalised node weights,

    dQ/dtheta = sum_j pi_j df/dtheta(z_j)        Fisher, the frozen rule's
              + a' dc/dtheta                     the centre moving
              + <Mbar, dM/dtheta>                the scale and log|S| moving

where `a = sum_j pi_j grad_u f(z_j)` and `Mbar` is what `tr(dS B) - logdet(M~)/2`
asks of `dM~` (`B = sqrt(2) sum_j pi_j x_j grad_u f(z_j)'`), carried back
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

For a Gaussian integrand the three placement terms cancel, as they must: the
rule is then exact wherever it sits.

# What is placed, and what is not

A flagged unit of one block on a product rule, at the node count the
continuation keeps for it (`o.unit_nodes`), placed as `_continuation_unit_rule`
places it: the leaf's mode by Newton (`_quadrature_leaf_rule!`), from the mode
it had at the last evaluation (`anchors`), where the continuation searches
from the Laplace mode. A unit that the continuation gave a soft rule (more
than `product_maxdim` effects) gets the product rule here, whose placement
this differentiates. A unit of one root block over leaves is placed as the
continuation places it, the root at the joint mode and eliminated precision and
each leaf at its conditional mode at every root node (`_placed_nested_unit!`),
which chains a second implicit derivative. A unit on importance draws keeps its
frozen draws, as does a deeper tree for now: their rules are the
continuation's as placed, with its Fisher gradient. The other units keep their
Laplace terms (`o.rest`), and the prior is theirs.

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
    ctsem_laplace_placed_hessian, ctsem_laplace_placed_info

"""
    CTSEMLaplacePlaced

The hybrid with its flagged units' product rules placed at every evaluation;
see the module docstring. `cont` is the continuation it was built from, which
supplies the flagged set, the Laplace part (`cont.rest`), the frozen rules and
the node counts. `placed` are the units placed at every evaluation, `frozen`
positions in `cont.rules` of the units kept on their frozen rules, `nodes`
each unit's node count (zero for a unit not placed), and `anchors` each placed
unit's effects at its mode at the last evaluation, where the next one's
Newton starts.
"""
mutable struct CTSEMLaplacePlaced{C} <: CTSEMOptimisable
    cont::C
    placed::Vector{Int}
    frozen::Vector{Int}
    nodes::Vector{Int}
    anchors::Vector{Vector{Float64}}
    cache::Any
    value_cache::Any
    value_calls::Int
    gradient_calls::Int
    seconds_values::Float64
    seconds_gradients::Float64
    seconds_hessian::Float64
    lbfgs::Union{Nothing,CTSEMLBFGSMemory}
end

"""
    ctsem_laplace_placed(o)

The placed objective from the continuation `o`, its flagged set and node counts
as `o` holds them: every flagged unit of one block, or of one root block over
leaves (`_placed_two_level`), whose rule is not draws is placed at every
evaluation; the rest of the flagged units keep `o`'s frozen rules. Each placed unit's anchor starts at its mode at `o`'s centre, from
`o`'s own search (`_laplace_solve_unit_mode!`), so at the centre the value is
`o`'s placed value wherever `o`'s rules are product rules.
"""
function ctsem_laplace_placed(o::CTSEMLaplaceContinuation)
    laplace = o.laplace
    nunits = length(laplace.units.members)
    placed = Int[]
    frozen = Int[]
    nodes = zeros(Int, nunits)
    for (f, U) in enumerate(o.flagged)
        blocks = laplace.units.blocks[U]
        if !o.importance[U] && ((length(blocks) == 1 && blocks[1].size >= 1) ||
                _placed_two_level(laplace, U))
            push!(placed, U)
            nodes[U] = max(o.nodes, o.unit_nodes[U])
        else
            push!(frozen, f)
        end
    end
    theta = copy(o.centre)
    Ls = _laplace_popchols(theta, laplace.spec)
    anchors = [Float64[] for _ in 1:nunits]
    _laplace_parallel(laplace, eachindex(placed)) do i
        local U
        U = placed[i]
        _laplace_solve_unit_mode!(laplace, U, theta, Ls)
        anchors[U] = copy(laplace.modes[U])
        return true
    end
    return CTSEMLaplacePlaced{typeof(o)}(o, placed, frozen, nodes, anchors, nothing,
        nothing, 0, 0, 0.0, 0.0, 0.0, nothing)
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
`tr(dS B) + logdet(S)` for a change `dM` of the unclipped precision, `S = U^-1`
and `M~ = U'U` (`_placed_scale`). Through the factor: `dS = -S dU S`, so
`tr(dS B) = <Ubar, dU>` with `Ubar = -triu((S B S)')`, and the upper Cholesky's
reverse, `Mbar = sym(S Phi(Ubar U') S')` with `Phi` the upper triangle and half
the diagonal; `logdet(S) = -logdet(M~)/2` adds `-S S'/2`. Through the clip,
`M~ = V max(Lambda, 1) V'`, by the divided differences of `max(., 1)` on the
eigenvalues (Daleckii-Krein).
"""
function _placed_mbar(sc, B::Matrix{Float64}, k::Int)
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
    Mbar = (W .+ transpose(W)) ./ 2 .- SS ./ 2
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
its precision, the product rule there; see the module docstring. With
`want_gradient`, its gradient in `theta` at those ancestors into `row`, and
into `wv` the raw-space vector whose loadings give its derivative in an
ancestor's effects (`L_a' wv[rho_a]`): the node-weighted mean of the members'
raw gradients, plus the gradient's derivative along the mode's direction `L v`,
less the jets' `grad <Mbar, L'HL>` -- the same three terms as the `theta`
gradient, without the chain through `L(theta)`, since an ancestor's effects
reach the members only through their raw values. `u`'s leaf columns are left at
the conditional mode. Returns `(ok, value)`.
"""
function _placed_leaf!(row::Vector{Float64}, wv::Vector{Float64}, P::CTSEMLaplacePlaced,
    U::Int, b::Int, u::Vector{Float64}, start::Vector{Float64},
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool)
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
    failure = (ok=false, value=NaN)
    aws = _laplace_workspace!(laplace, Float64, npar)
    lr = _quadrature_leaf_rule!(laplace, U, theta, Ls, b, u, aws; start=start)
    # The implicit function theorem is the centre's derivative only at a
    # stationary point; a solve stopped short of one has none to give.
    (lr.ok && lr.gradmax <= _PLACED_STATIONARY) || return failure
    c = lr.centre
    M = Matrix{Float64}(lr.precision)
    sc = _placed_scale(M, k)
    sc.ok || return failure
    grid, lw = _gh_grid(k, P.nodes[U])
    scratch = _continuation_scratch(laplace, npar)
    peak = -Inf
    acc = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    vacc = want_gradient ? zeros(npar) : Float64[]
    aacc = zeros(k)
    Bacc = zeros(k, k)
    gh = zeros(npar)
    gv = zeros(npar)
    gu = zeros(d)
    z = zeros(k)
    sx = zeros(k)
    for j in eachindex(grid)
        x = grid[j]
        _ctsem_mulvec!(sx, sc.S, x)
        @inbounds for t in 1:k
            z[t] = c[t] + sqrt(2) * sx[t]
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
        if want_gradient
            gacc .+= w .* gh
            vacc .+= w .* gv
            @inbounds for i in 1:k
                fz = gu[columns[i]] - z[i]
                aacc[i] += w * fz
                for r in 1:k
                    Bacc[r, i] += w * x[r] * fz
                end
            end
        end
    end
    isfinite(peak) || return failure
    value = peak + log(acc) + sc.logdetscale + k * (log(2) - log(2 * pi)) / 2
    @inbounds for (t, col) in enumerate(columns)
        u[col] = c[t]
    end
    want_gradient || return (ok=true, value=value)

    fisher = gacc ./ acc
    a = aacc ./ acc
    B = (sqrt(2) / acc) .* Bacc
    Mbar = _placed_mbar(sc, B, k)
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
    return (ok=true, value=value)
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
`row` (unweighted): a unit of one block is one leaf (`_placed_leaf!`) from its
anchor, a unit of a root over leaves an outer rule over them
(`_placed_nested_unit!`). Returns `(ok, value, anchor)`, `anchor` the unit's
effects at the mode it found.
"""
function _placed_unit!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool)
    laplace = P.cont.laplace
    length(laplace.units.blocks[U]) == 1 ||
        return _placed_nested_unit!(row, P, U, theta, Ls, dL, positions, want_gradient)
    block = laplace.units.blocks[U][1]
    columns = (block.offset + 1):(block.offset + block.size)
    u = copy(P.anchors[U])
    r = _placed_leaf!(row, zeros(length(theta)), P, U, 1, u,
        Float64[u[c] for c in columns], theta, Ls, dL, positions, want_gradient)
    return (ok=r.ok, value=r.value, anchor=r.ok ? u : P.anchors[U])
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
effects is `L0' wv[rho0]` (`_placed_leaf!`).
"""
function _placed_nested_unit!(row::Vector{Float64}, P::CTSEMLaplacePlaced, U::Int,
    theta::Vector{Float64}, Ls, dL, positions, want_gradient::Bool)
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
    failure = (ok=false, value=NaN, anchor=P.anchors[U])
    aws = _laplace_workspace!(laplace, Float64, npar)
    jm = _laplace_newton_unit_mode(laplace, U, theta, Ls, aws, copy(P.anchors[U]))
    (jm.converged && isfinite(jm.value)) || return failure
    cj = collect(Float64, jm.u)
    M = _laplace_copy_blockmatrix(_laplace_unit_curvature(laplace, U, theta, Ls, cj))
    fz = _laplace_factor_repaired!(_laplace_copy_blockmatrix(M), blocks)
    (fz.ok && !fz.repaired) || return failure
    sc = _placed_scale_factor(fz.factors[r0], k0)
    sc.ok || return failure
    grid, lw = _gh_grid(k0, P.nodes[U])
    peak = -Inf
    acc = 0.0
    gacc = want_gradient ? zeros(npar) : Float64[]
    aacc = zeros(k0)
    Bacc = zeros(k0, k0)
    rowc = zeros(npar)
    wvc = zeros(npar)
    G = zeros(k0)
    gd = zeros(d)
    z = zeros(k0)
    sx = zeros(k0)
    rowsum = zeros(npar)
    for j in eachindex(grid)
        x = grid[j]
        _ctsem_mulvec!(sx, sc.S, x)
        u = copy(cj)
        @inbounds for t in 1:k0
            z[t] = cj[cols0[t]] + sqrt(2) * sx[t]
            u[cols0[t]] = z[t]
        end
        want_gradient && (fill!(rowsum, 0.0); fill!(G, 0.0))
        h = 0.0
        ok = true
        for cb in kids
            ccols = (blocks[cb].offset + 1):(blocks[cb].offset + blocks[cb].size)
            leaf = _placed_leaf!(rowc, wvc, P, U, cb, u, Float64[cj[c] for c in ccols],
                theta, Ls, dL, positions, want_gradient)
            if !leaf.ok
                ok = false
                break
            end
            h += leaf.value
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
        if want_gradient
            gacc .+= w .* rowsum
            @inbounds for i in 1:k0
                fz0 = G[i] - z[i]
                aacc[i] += w * fz0
                for r in 1:k0
                    Bacc[r, i] += w * x[r] * fz0
                end
            end
        end
    end
    isfinite(peak) || return failure
    value = peak + log(acc) + sc.logdetscale + k0 * (log(2) - log(2 * pi)) / 2
    want_gradient || return (ok=true, value=value, anchor=cj)

    fisher = gacc ./ acc
    a0 = aacc ./ acc
    B0 = (sqrt(2) / acc) .* Bacc
    Mbar0 = _placed_mbar(sc, B0, k0)
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
    return (ok=true, value=value, anchor=cj)
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
        local rule, U, aws, scratch, u, total, r
        rule = o.rules[P.frozen[i]]
        U = rule.unit
        aws = _laplace_workspace!(laplace, Float64, npar)
        scratch = _continuation_scratch(laplace, npar)
        u = zeros(laplace.units.dims[U])
        total = 0.0
        for root in rule.roots
            r = _continuation_rule_value(o, U, theta, Ls, dL, positions, root, u, aws,
                scratch, want_gradient)
            isfinite(r.value) || return false
            total += r.value
            want_gradient && (rows[i] .+= r.gradient)
        end
        values[i] = total
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
            return _laplace_free_result(_continuation_copy(hit.result), n)
        end
    end
    result = _placed_evaluate(P, theta, gradient)
    entry = (theta=theta, result=_continuation_copy(result))
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
    sizes = unique(b.size for U in P.placed for b in laplace.units.blocks[U])
    for m in unique(P.nodes[P.placed])
        _quadrature_warm_caches(laplace, m, sizes)
    end
    np = length(P.placed)
    pvalues = fill(NaN, np)
    prows = [gradient ? zeros(npar) : Float64[] for _ in 1:np]
    anchors = Vector{Vector{Float64}}(undef, np)
    placed_ok = _laplace_parallel(laplace, 1:np) do i
        local U, r
        U = P.placed[i]
        r = try
            _placed_unit!(prows[i], P, U, theta, Ls, dL, positions, gradient)
        catch err
            _ctsem_must_propagate(err) && rethrow()
            (ok=false, value=NaN, anchor=P.anchors[U])
        end
        pvalues[i] = r.value
        anchors[i] = r.anchor
        return r.ok
    end
    if placed_ok
        for (i, U) in enumerate(P.placed)
            P.anchors[U] = anchors[i]
        end
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
    home = deepcopy(P.anchors)
    idx = index === nothing ? collect(eachindex(x)) : collect(Int, index)
    fd = fill(NaN, length(idx))
    for (j, t) in enumerate(idx)
        h = step * max(1.0, abs(x[t]))
        xp = copy(x); xp[t] += h
        xm = copy(x); xm[t] -= h
        P.anchors = deepcopy(home)
        P.cache = nothing; P.value_cache = nothing
        vp = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xp), false).value
        P.anchors = deepcopy(home)
        vm = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xm), false).value
        fd[j] = (vp - vm) / (2 * h)
    end
    P.anchors = home
    P.cache = nothing; P.value_cache = nothing
    return (value=r.value, gradient=r.gradient[idx], differences=fd, index=idx)
end

"""
    ctsem_laplace_placed_hessian(P, values; step=1e-4)

The placed objective's Hessian at `values` by central differences of its exact
gradient, two gradients a column, every unit's Newton starting from its mode at
`values`.
"""
function ctsem_laplace_placed_hessian(P::CTSEMLaplacePlaced, values::AbstractVector;
    step::Real=1e-4)
    started = time_ns()
    x = collect(Float64, values)
    n = length(x)
    ctsem_laplace_placed_evaluate(P, x; gradient=false)
    home = deepcopy(P.anchors)
    H = fill(NaN, n, n)
    for t in 1:n
        h = step * max(1.0, abs(x[t]))
        xp = copy(x); xp[t] += h
        xm = copy(x); xm[t] -= h
        P.anchors = deepcopy(home)
        gp = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xp), true)
        P.anchors = deepcopy(home)
        gm = _placed_evaluate(P, _laplace_extend(P.cont.laplace, xm), true)
        (gp.converged && gm.converged) || continue
        H[:, t] .= (gp.gradient[1:n] .- gm.gradient[1:n]) ./ (2 * h)
    end
    P.anchors = home
    P.cache = nothing
    P.value_cache = nothing
    H = (H .+ transpose(H)) ./ 2
    P.seconds_hessian += (time_ns() - started) / 1e9
    return H
end

"""
    ctsem_laplace_placed_info(P)

The placed objective's accounting: units placed and frozen, calls, seconds;
`soft` the placed units wider than the continuation's `product_maxdim`, which
the rounds give a soft rule and this a product rule. A list with no unit is
`[0]`, for the bridge (an empty vector hangs it).
"""
function ctsem_laplace_placed_info(P::CTSEMLaplacePlaced)
    nonempty(v) = isempty(v) ? [0] : v
    soft = [U for U in P.placed if any(b -> b.size > P.cont.product_maxdim,
        P.cont.laplace.units.blocks[U])]
    return (placed=nonempty(copy(P.placed)),
        frozen=nonempty([P.cont.rules[f].unit for f in P.frozen]),
        nodes=nonempty(P.nodes[P.placed]), soft=nonempty(soft),
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
    P.lbfgs === nothing ? _ctsem_lbfgs_held(P.cont) : P.lbfgs
_ctsem_lbfgs_hold!(P::CTSEMLaplacePlaced, memory::CTSEMLBFGSMemory) =
    (P.lbfgs = memory; nothing)
ctsem_hessian(P::CTSEMLaplacePlaced, values::AbstractVector; chunk::Integer=0,
    progress=nothing) = ctsem_laplace_placed_hessian(P, values)
_ctsem_params(P::CTSEMLaplacePlaced) = _ctsem_params(P.cont.laplace)
_ctsem_saturated_for(P::CTSEMLaplacePlaced, minimizer) =
    _laplace_saturated_parameters(P.cont.laplace, minimizer)
_ctsem_optimise_result_extra(P::CTSEMLaplacePlaced, final, log) =
    (placed_units=length(P.placed), frozen_units=length(P.frozen))

function _ctsem_optimise_trial(P::CTSEMLaplacePlaced, x, want_gradient::Bool,
    gradient_method, limit::Real, log)
    evaluated = _ctsem_try_evaluate(() ->
        ctsem_laplace_placed_evaluate(P, x; gradient=want_gradient))
    valid = evaluated !== nothing && isfinite(evaluated.value) && evaluated.converged &&
        (!want_gradient || _ctsem_gradient_usable(evaluated.gradient, limit))
    return (evaluated=evaluated, valid=valid)
end
