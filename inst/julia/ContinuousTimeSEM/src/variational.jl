"""
The variational objective: each unit's random effects integrated by a lower
bound on its marginal likelihood, fitted at every evaluation, in place of the
Laplace term (review/DESIGN-unit-bound-2026-11-26.md).

A unit's approximation maps a standard normal `z` to its standardised
effects `u` on the hybrid's coordinates (the deviations the filter integrates
stay integrated exactly), `u = m + L x(z)`, and its bound is

    B_U = sum_members E[log L_m(u)] - E[|u|^2] / 2 + E[log |du/dz|] + dim / 2

less a weak prior on the family's shape (`kappa`). Each member's term is an
integral over its own chain of the unit's coordinates alone -- its own
block's and its ancestors' -- on nodes over that chain. The other two terms
involve no model: for the Gaussian family both are in closed form; for the
shaped family the entropy is on fine one-dimensional rules, so the family
cannot gain by stretching its map between the nodes, and `-E|u|^2 / 2` on the
unit's own nodes, whose points and Jacobians are there already. The approximation maximises `B_U`, so the bound's gradient in the
parameters is the partial with the nodes held: first derivatives of the
members' likelihoods at their nodes.

Two families:

- The Gaussian family on the unit's block tree (`_vf_tree_layout`), for every
  unit by default. A coordinate's row of `L` reaches its own block's earlier
  coordinates and every coordinate of its block's ancestors, nothing else, so
  the parameters grow linearly in the members. That is exact for a unit whose
  likelihood is Gaussian in its effects: the posterior precision has no
  coupling between sibling blocks, so in root-first order the covariance's
  factor has none either. Rules over a member's chain: product Gauss-Hermite,
  3 nodes a coordinate, up to three coordinates; the rotated `2k + 1` rule
  beyond (unbiased over rotations, not a guaranteed bound member by member).
- The shaped family (`CTSEMVFamily`) for a unit of one block and at most three
  coordinates: each coordinate skewed (`y = sinh((asinh z + e) / d)`), and at
  two or three bending with the earlier ones, `x_j = y_j exp(sum_i c_ji
  s(y_i)) + sum_i b_ji (s(y_i)^2 - 1)`, `s(y) = 5 tanh(y / 5)`. One coordinate
  is skewed from the start (7 nodes); a unit of two or three starts Gaussian
  and takes the bending family where it gains (`_ctsem_take_search_modes!`, 5
  nodes a coordinate).

Every rule has positive weights. A rule with negative weights is no bound: on
the cut rule (each axis by Gauss-Hermite, the rest by the degree-5 rule's
off-axis points) one of gA1's three-coordinate units moved its approximation's
mean 100 prior sds out, where the negative-weight nodes read a likelihood of
-2e7, for a "bound" of +1e6.

The inner solve is Newton on the bound's Hessian in the family's parameters,
which has the unit's tree: a member's term couples its own block's parameters
with its ancestors' and never with a sibling's. So the step is the Laplace
route's block elimination (`_laplace_block_factor`, `_laplace_block_solve`)
over the parameters' blocks, a Levenberg shift where the matrix is not
definite. The members' shares of the curvature need their likelihood's
curvature at every node, the dear part, so each unit holds them from where
they were last formed and steps on gradient passes alone, each member's share
updated by SR1 from its own gradient's change (as `_LaplacePSR1` holds the
Laplace route's); the curvature is formed afresh only where the held shares do
not serve (`_vf_solve_held`). The prior and entropy terms are exact at every
step.

Within an optimisation (`running`) each evaluation starts every unit from the
approximation the last accepted point left, and the optimiser's per-step
callback (`_ctsem_anchor_here!`) makes the last evaluation's approximations
the accepted ones, as the Laplace route moves its anchors. Outside one, an
evaluation starts from the accepted approximations and keeps nothing, so a
value is the parameters' given the held state.

The log likelihood a fit reports is not this file's: R takes the
random-effects check's integral at the estimate (`ctsem_laplace_check` on the
Laplace objective underneath), which integrates nested units block by block.
"""

export ctsem_variational, ctsem_variational_info, ctsem_variational_reset!,
    ctsem_variational_choose!

# ---- the shaped family (one block, at most three coordinates) ----------------------

"""The shaped family: `d` coordinates, which axes are skewed, whether it bends,
and where each parameter sits in its vector (`m`, `log diag L`, the strict
lower triangle of `L` by column, `(e_j, log d_j)` for each skewed axis, then
`(b_ji, c_ji)` for i < j where it bends). Its first `2d + d(d - 1)/2`
parameters are laid out as the Gaussian tree family's on one block, so a
Gaussian approximation carries over by appending zeros."""
struct CTSEMVFamily
    d::Int
    skew::Vector{Bool}
    bend::Bool
    im::UnitRange{Int}
    ild::UnitRange{Int}
    ilo::Matrix{Int}
    ie::Vector{Int}
    idl::Vector{Int}
    ib::Matrix{Int}
    ic::Matrix{Int}
    p::Int
    shape::Vector{Int}
end

function CTSEMVFamily(d::Integer, skew::Bool; bend::Bool=false)
    d = Int(d)
    at = 0
    im = (at + 1):(at + d)
    at += d
    ild = (at + 1):(at + d)
    at += d
    ilo = zeros(Int, d, d)
    for j in 1:d, i in (j + 1):d
        at += 1
        ilo[i, j] = at
    end
    ie = zeros(Int, d)
    idl = zeros(Int, d)
    if skew
        for j in 1:d
            at += 1
            ie[j] = at
            at += 1
            idl[j] = at
        end
    end
    bnd = bend && d > 1
    ib = zeros(Int, d, d)
    ic = zeros(Int, d, d)
    if bnd
        for j in 2:d, i in 1:(j - 1)
            at += 1
            ib[j, i] = at
            at += 1
            ic[j, i] = at
        end
    end
    shape = vcat(skew ? vcat(ie, idl) : Int[], filter(>(0), vec(ib)), filter(>(0), vec(ic)))
    return CTSEMVFamily(d, fill(skew, d), bnd, im, ild, ilo, ie, idl, ib, ic, at, shape)
end

# The shaped family's map at any element type, for its derivatives: u and
# log |du/dz|. Small and free of the filter, so differentiated by ForwardDiff
# on the callable structs below.
function _vf_map_generic(f::CTSEMVFamily, lam::AbstractVector, z::AbstractVector)
    T = promote_type(eltype(lam), eltype(z))
    d = f.d
    y = Vector{T}(undef, d)
    lj = zero(T)
    for j in 1:d
        if f.skew[j]
            dl = lam[f.idl[j]]
            a = (asinh(z[j]) + lam[f.ie[j]]) / exp(dl)
            y[j] = sinh(a)
            lj += _vf_lcosh(a) - dl - log1p(z[j]^2) / 2
        else
            y[j] = z[j]
        end
    end
    x = Vector{T}(undef, d)
    for j in 1:d
        if f.bend && j > 1
            loc = zero(T)
            sj = zero(T)
            for i in 1:(j - 1)
                si = 5 * tanh(y[i] / 5)
                loc += lam[f.ib[j, i]] * (si^2 - 1)
                sj += lam[f.ic[j, i]] * si
            end
            x[j] = y[j] * exp(sj) + loc
            lj += sj
        else
            x[j] = y[j]
        end
    end
    u = Vector{T}(undef, d)
    for i in 1:d
        li = lam[f.ild[i]]
        acc = lam[f.im[i]] + exp(li) * x[i]
        for j in 1:(i - 1)
            acc += lam[f.ilo[i, j]] * x[j]
        end
        u[i] = acc
        lj += li
    end
    return u, lj
end

struct _VFMapAt
    f::CTSEMVFamily
    z::Vector{Float64}
end
(m::_VFMapAt)(lam) = (q = _vf_map_generic(m.f, lam, m.z); vcat(q[1], q[2]))

# r'u alone, for the likelihood term's second derivatives through the map.
struct _VFUAt
    f::CTSEMVFamily
    z::Vector{Float64}
    r::Vector{Float64}
end
(m::_VFUAt)(lam) = dot(m.r, _vf_map_generic(m.f, lam, m.z)[1])

# log cosh, smooth at any element type.
_vf_lcosh(a) = abs(a) > 20 ? abs(a) - log(2.0) + log1p(exp(-2 * abs(a))) : log(cosh(a))

# One axis's transform at node value z: y and its derivatives in (e, log d).
function _vf_axis(f::CTSEMVFamily, lam::Vector{Float64}, z::Float64, j::Int)
    f.skew[j] || return (y=z, ye=0.0, yd=0.0, yee=0.0, yed=0.0, ydd=0.0)
    del = exp(lam[f.idl[j]])
    a = (asinh(z) + lam[f.ie[j]]) / del
    sh = sinh(a)
    ch = cosh(a)
    return (y=sh, ye=ch / del, yd=-a * ch, yee=sh / del^2, yed=-(a * sh + ch) / del,
        ydd=a * ch + a^2 * sh)
end

function _vf_L(f::CTSEMVFamily, lam::Vector{Float64})
    d = f.d
    L = zeros(d, d)
    @inbounds for i in 1:d
        L[i, i] = exp(lam[f.ild[i]])
        for j in 1:(i - 1)
            L[i, j] = lam[f.ilo[i, j]]
        end
    end
    return L
end

"""The nodes' points `U` (d x K), their log Jacobians, and each node's
Jacobian of `[u; log |du/dz|]` in the family's parameters."""
function _vf_nodes(f::CTSEMVFamily, lam::Vector{Float64}, Z::Matrix{Float64};
    jacobian::Bool=true)
    f.bend && return _vf_nodes_bend(f, lam, Z, jacobian)
    d, K = size(Z)
    L = _vf_L(f, lam)
    U = zeros(d, K)
    lj = zeros(K)
    J = Vector{Matrix{Float64}}(undef, K)
    ldL = 0.0
    for i in 1:d
        ldL += lam[f.ild[i]]
    end
    y = zeros(d)
    for k in 1:K
        Jk = zeros(d + 1, f.p)
        l = ldL
        for j in 1:d
            z = Z[j, k]
            ax = _vf_axis(f, lam, z, j)
            y[j] = ax.y
            if f.skew[j]
                dl = lam[f.idl[j]]
                a = (asinh(z) + lam[f.ie[j]]) / exp(dl)
                l += _vf_lcosh(a) - dl - log1p(z^2) / 2
                th = tanh(a)
                for i in j:d
                    Jk[i, f.ie[j]] = L[i, j] * ax.ye
                    Jk[i, f.idl[j]] = L[i, j] * ax.yd
                end
                Jk[d + 1, f.ie[j]] = th / exp(dl)
                Jk[d + 1, f.idl[j]] = -a * th - 1
            end
        end
        for i in 1:d
            acc = lam[f.im[i]]
            for j in 1:i
                acc += L[i, j] * y[j]
            end
            U[i, k] = acc
            Jk[i, f.im[i]] = 1.0
            Jk[i, f.ild[i]] = L[i, i] * y[i]
            Jk[d + 1, f.ild[i]] = 1.0
            for j in 1:(i - 1)
                Jk[i, f.ilo[i, j]] = y[j]
            end
        end
        lj[k] = l
        J[k] = Jk
    end
    return (U=U, lj=lj, J=J)
end

function _vf_nodes_bend(f::CTSEMVFamily, lam::Vector{Float64}, Z::Matrix{Float64},
    jacobian::Bool)
    d, K = size(Z)
    U = zeros(d, K)
    lj = zeros(K)
    J = Vector{Matrix{Float64}}(undef, K)
    for k in 1:K
        z = Z[:, k]
        u, l = _vf_map_generic(f, lam, z)
        U[:, k] .= u
        lj[k] = l
        J[k] = jacobian ? ForwardDiff.jacobian(_VFMapAt(f, z), lam) : zeros(0, 0)
    end
    return (U=U, lj=lj, J=J)
end

# The Hessian in the family's parameters of r'u at one node, the skewed
# family by hand (the bending family's is ForwardDiff's, `_VFUAt`).
function _vf_second!(H::Matrix{Float64}, wk::Float64, f::CTSEMVFamily,
    lam::Vector{Float64}, L::Matrix{Float64}, z::AbstractVector, r::AbstractVector)
    d = f.d
    for j in 1:d
        ax = _vf_axis(f, lam, z[j], j)
        dj = f.ild[j]
        H[dj, dj] += wk * r[j] * L[j, j] * ax.y
        f.skew[j] || continue
        e = f.ie[j]
        g = f.idl[j]
        t = wk * r[j] * L[j, j]
        H[dj, e] += t * ax.ye
        H[e, dj] += t * ax.ye
        H[dj, g] += t * ax.yd
        H[g, dj] += t * ax.yd
        for i in (j + 1):d
            o = f.ilo[i, j]
            H[o, e] += wk * r[i] * ax.ye
            H[e, o] += wk * r[i] * ax.ye
            H[o, g] += wk * r[i] * ax.yd
            H[g, o] += wk * r[i] * ax.yd
        end
        rl = 0.0
        for i in j:d
            rl += r[i] * L[i, j]
        end
        H[e, e] += wk * rl * ax.yee
        H[e, g] += wk * rl * ax.yed
        H[g, e] += wk * rl * ax.yed
        H[g, g] += wk * rl * ax.ydd
    end
    return H
end

"""E[log |du/dz|] under the standard normal, at any element type: log |L|,
each skewed axis's term, and each bending coefficient times the mean of
`s(y_i)`, every one an integral over one coordinate, on 96-node
Gauss-Hermite."""
function _vf_entropy(f::CTSEMVFamily, lam::AbstractVector)
    z1, w1 = _vf_gh(96)
    T = eltype(lam)
    out = zero(T)
    for i in 1:f.d
        out += lam[f.ild[i]]
    end
    (any(f.skew) || f.bend) || return out
    for i in 1:f.d
        sy = zero(T)
        for k in eachindex(z1)
            if f.skew[i]
                dl = lam[f.idl[i]]
                a = (asinh(z1[k]) + lam[f.ie[i]]) / exp(dl)
                out += w1[k] * (_vf_lcosh(a) - dl - log1p(z1[k]^2) / 2)
                f.bend && (sy += w1[k] * 5 * tanh(sinh(a) / 5))
            elseif f.bend
                sy += w1[k] * 5 * tanh(z1[k] / 5)
            end
        end
        if f.bend
            for j in (i + 1):f.d
                out += lam[f.ic[j, i]] * sy
            end
        end
    end
    return out
end

struct _VFEntropyAt
    f::CTSEMVFamily
end
(m::_VFEntropyAt)(lam) = _vf_entropy(m.f, lam)

# ---- rules --------------------------------------------------------------------------

"""Gauss-Hermite for the standard normal, `n` nodes: the engine's rule
(`_gauss_hermite`, for weight `exp(-x^2)`, cached) rescaled."""
function _vf_gh(n::Integer)
    x, w = _gauss_hermite(n)
    return sqrt(2.0) .* x, w ./ sqrt(pi)
end

"""Product Gauss-Hermite over `d` coordinates for the standard normal, `n` a
coordinate: the placed rules' product set (`_point_set_product`, which
refuses one too large to evaluate), its weights read against the normal."""
function _vf_rule_product(d::Integer, n::Integer)
    d == 0 && return zeros(0, 1), [1.0]
    ps = _point_set_product(d, n)
    Y = ps.points
    w = [exp(ps.logweights[k] - sum(abs2, view(Y, :, k)) / 2 - d * log(2pi) / 2)
         for k in axes(Y, 2)]
    return Y, w ./ sum(w)
end

"""The centre and +-sqrt(d + 2) along `R` random rotations, drawn once
(`seed`) and kept: exact to degree 3, unbiased to degree 5 over rotations,
`2Rd + 1` nodes, positive weights."""
function _vf_rule_rot(d::Integer, R::Integer, seed::Integer)
    d == 0 && return zeros(0, 1), [1.0]
    rng = Random.Xoshiro(seed)
    Z = zeros(d, 1)
    w = [2 / (d + 2)]
    for _ in 1:R
        Q = Matrix(qr(randn(rng, d, d)).Q)[:, 1:d]
        Z = hcat(Z, sqrt(d + 2) .* Q, -sqrt(d + 2) .* Q)
        append!(w, fill(1 / (2 * R * (d + 2)), 2 * d))
    end
    return Z, w
end

"""A Gaussian member's rule over its chain of `k` coordinates: product, 3 a
coordinate, to three; rotated beyond, seeded by the unit and member."""
_vf_gauss_rule(k::Int, rule::Symbol, seed::Int) =
    (rule === :product || (rule === :auto && k <= 3)) ? _vf_rule_product(k, 3) :
    _vf_rule_rot(k, 1, seed)

# ---- a unit's layout: its family's parameters on its block tree --------------------

"""One member of a unit: its chain of coordinates (outermost level first), the
family's parameters its term reads (ascending), and for each chain
coordinate the local index of its mean and log diagonal and its row's
`(chain position, local parameter)` pairs (Gaussian family), with its rule."""
struct _VFMember
    member::Int
    chain::Vector{Int}
    plist::Vector{Int}
    jm::Vector{Int}
    jld::Vector{Int}
    jrow::Vector{Vector{Tuple{Int,Int}}}
    Z::Matrix{Float64}
    w::Vector{Float64}
end

"""A unit's family: Gaussian on its block tree, or shaped on one block; the
family's parameters' blocks (the unit's tree over them, `pblocks`, with each
parameter's block and place in it), each coordinate's mean, log diagonal and
row parameters (Gaussian), and its members."""
struct CTSEMVLayout
    shaped::Union{Nothing,CTSEMVFamily}
    d::Int
    p::Int
    pblocks::Vector{CTSEMLaplaceBlock}
    pblock::Vector{Int}
    plocal::Vector{Int}
    im::Vector{Int}
    ild::Vector{Int}
    rowp::Vector{Vector{Int}}
    members::Vector{_VFMember}
    rule::Symbol
end

# Each member's chain of coordinates, outermost level first.
function _vf_chains(laplace::CTSEMLaplaceObjective, U::Int)
    spec = laplace.spec
    out = Vector{Vector{Int}}()
    for mm in eachindex(laplace.units.members[U])
        offs = laplace.units.offsets[U][mm]
        chain = Int[]
        for l in reverse(eachindex(spec.levels))
            r = nlatent(spec.levels[l])
            r == 0 && continue
            append!(chain, offs[l] .+ (1:r))
        end
        push!(out, chain)
    end
    return out
end

"""The Gaussian family on unit `U`'s block tree. Parameters by block (the
unit's blocks, innermost first): each coordinate's mean, then log diagonal,
the block's own strict lower triangle by column, then its rows' entries on
each ancestor's coordinates."""
function _vf_tree_layout(laplace::CTSEMLaplaceObjective, U::Int, rule::Symbol, seed::Int)
    blocks = laplace.units.blocks[U]
    d = laplace.units.dims[U]
    im = zeros(Int, d)
    ild = zeros(Int, d)
    rowp = [Int[] for _ in 1:d]
    rowc = [Int[] for _ in 1:d]
    pblocks = CTSEMLaplaceBlock[]
    pblock = Int[]
    plocal = Int[]
    at = 0
    for (b, blk) in enumerate(blocks)
        start = at
        coords = (blk.offset + 1):(blk.offset + blk.size)
        k = blk.size
        for i in coords
            at += 1
            im[i] = at
        end
        for i in coords
            at += 1
            ild[i] = at
        end
        for jl in 1:k, il in (jl + 1):k
            at += 1
            push!(rowp[coords[il]], at)
            push!(rowc[coords[il]], coords[jl])
        end
        for a in blk.ancestors
            ab = blocks[a]
            for i in coords, j in (ab.offset + 1):(ab.offset + ab.size)
                at += 1
                push!(rowp[i], at)
                push!(rowc[i], j)
            end
        end
        push!(pblocks, CTSEMLaplaceBlock(start, at - start, blk.members, blk.level,
            blk.ancestors))
        append!(pblock, fill(b, at - start))
        append!(plocal, 1:(at - start))
    end
    members = _VFMember[]
    for (mm, chain) in enumerate(_vf_chains(laplace, U))
        plist = Int[]
        for i in chain
            push!(plist, im[i], ild[i])
            append!(plist, rowp[i])
        end
        sort!(plist)
        local_of(q) = searchsortedfirst(plist, q)
        pos = Dict(c => a for (a, c) in enumerate(chain))
        jrow = [[(pos[rowc[i][t]], local_of(rowp[i][t])) for t in eachindex(rowp[i])]
                for i in chain]
        Z, w = _vf_gauss_rule(length(chain), rule, Int(hash((seed, U, mm)) >> 2))
        push!(members, _VFMember(mm, chain, plist, [local_of(im[i]) for i in chain],
            [local_of(ild[i]) for i in chain], jrow, Z, w))
    end
    return CTSEMVLayout(nothing, d, at, pblocks, pblock, plocal, im, ild, rowp, members,
        rule)
end

"""The shaped family on unit `U`, one block of `f.d` coordinates: every
member's chain is the whole unit, and its nodes are the unit's product rule."""
function _vf_shaped_layout(laplace::CTSEMLaplaceObjective, U::Int, f::CTSEMVFamily,
    nodes::Int)
    blocks = laplace.units.blocks[U]
    (length(blocks) == 1 && f.d == laplace.units.dims[U]) ||
        throw(ArgumentError("the shaped family is for a unit of one block"))
    Z, w = _vf_rule_product(f.d, nodes)
    all1 = collect(1:f.p)
    members = [_VFMember(mm, collect(1:f.d), all1, Int[], Int[], Vector{Tuple{Int,Int}}[],
        Z, w) for mm in eachindex(laplace.units.members[U])]
    pblocks = [CTSEMLaplaceBlock(0, f.p, blocks[1].members, blocks[1].level, Int[])]
    return CTSEMVLayout(f, f.d, f.p, pblocks, fill(1, f.p), collect(1:f.p), Int[], Int[],
        Vector{Int}[], members, :product)
end

# ---- a unit's terms --------------------------------------------------------------------

struct _VCtx
    theta::Vector{Float64}
    Ls::Vector{Matrix{Float64}}
    dL::Any
    positions::Vector{Vector{Int}}
end

function _vf_ctx(laplace::CTSEMLaplaceObjective, values::AbstractVector)
    theta = Vector{Float64}(_laplace_extend(laplace, collect(Float64, values)))
    spec = laplace.spec
    Ls = [Matrix{Float64}(M) for M in _laplace_popchols(theta, spec)]
    return _VCtx(theta, Ls, _laplace_level_chol_derivatives(theta, spec),
        [_laplace_level_positions(spec, l) for l in eachindex(spec.levels)])
end

"""The prior and entropy terms, `-E|u|^2 / 2 + E[log |du/dz|] + d / 2`: value,
gradient, and the negated Hessian. For the Gaussian tree family a vector, its
diagonal: `-E|u|^2 / 2 = -sum_i (m_i^2 + exp(2 l_i) + sum_j L_ij^2) / 2`
and the entropy `sum_i l_i`. For the shaped family a matrix: `-|u|^2 / 2` on
the unit's own nodes (`nd`, whose points and Jacobians the members' terms
use), the entropy on its fine one-dimensional rule (`_vf_entropy`) and the
shape prior; the Hessian only where `hessian`."""
function _vf_prior_entropy(lay::CTSEMVLayout, lam::Vector{Float64}, kappa::Float64,
    hessian::Bool, nd=nothing)
    if lay.shaped !== nothing
        f = lay.shaped
        ent = _VFEntropyAt(f)
        B = f.d / 2 + ent(lam)
        g = ForwardDiff.gradient(ent, lam)
        Hneg = hessian ? -ForwardDiff.hessian(ent, lam) : zeros(0, 0)
        for j in f.shape
            B -= kappa * lam[j]^2 / 2
            g[j] -= kappa * lam[j]
            hessian && (Hneg[j, j] += kappa)
        end
        Z = lay.members[1].Z
        w = lay.members[1].w
        L = hessian ? _vf_L(f, lam) : zeros(0, 0)
        for k in eachindex(w)
            u = view(nd.U, :, k)
            Ju = view(nd.J[k], 1:f.d, :)
            B -= w[k] * sum(abs2, u) / 2
            mul!(g, transpose(Ju), u, -w[k], 1.0)
            hessian || continue
            _ctsem_mulTN!(Hneg, Ju, Ju, w[k], true)
            if f.bend
                Hneg .+= w[k] .* ForwardDiff.hessian(_VFUAt(f, Z[:, k], Vector(u)), lam)
            else
                _vf_second!(Hneg, w[k], f, lam, L, view(Z, :, k), u)
            end
        end
        hessian && (Hneg .= (Hneg .+ transpose(Hneg)) ./ 2)
        return (B=B, g=g, Hneg=Hneg)
    end
    B = lay.d / 2
    g = zeros(lay.p)
    Hneg = zeros(lay.p)
    for i in 1:lay.d
        q = lay.im[i]
        B -= lam[q]^2 / 2
        g[q] = -lam[q]
        Hneg[q] = 1.0
        q = lay.ild[i]
        e2 = exp(2 * lam[q])
        B += lam[q] - e2 / 2
        g[q] = 1 - e2
        Hneg[q] = 2 * e2
        for q in lay.rowp[i]
            B -= lam[q]^2 / 2
            g[q] = -lam[q]
            Hneg[q] = 1.0
        end
    end
    return (B=B, g=g, Hneg=Hneg)
end

"""A Gaussian member's chain at node `z`: its points and their Jacobian in the
member's parameters (`u_a = m_a + exp(l_a) z_a + sum_j L_aj z_j`)."""
function _vf_tree_point!(uc::Vector{Float64}, Jn::Matrix{Float64}, mem::_VFMember,
    lam::Vector{Float64}, z::AbstractVector)
    fill!(Jn, 0.0)
    @inbounds for a in eachindex(mem.chain)
        qm = mem.plist[mem.jm[a]]
        ql = mem.plist[mem.jld[a]]
        Laa = exp(lam[ql])
        acc = lam[qm] + Laa * z[a]
        Jn[a, mem.jm[a]] = 1.0
        Jn[a, mem.jld[a]] = Laa * z[a]
        for (b, t) in mem.jrow[a]
            acc += lam[mem.plist[t]] * z[b]
            Jn[a, t] = z[b]
        end
        uc[a] = acc
    end
    return uc
end

"""Unit `U`'s bound at the family's parameters `lam`: value, gradient in them,
each member's gradient in its own parameters, and with `curvature` each
member's Hessian in them; the prior and entropy terms' negated Hessian; the
gradient in the population parameters (`GE`, the extended vector); counts
of gradient and curvature passes. `B = -Inf` where a member cannot be
evaluated at some node."""
function _vf_unit_terms(V, U::Int, ctx::_VCtx, lam::Vector{Float64}, curvature::Bool)
    lay = V.layouts[U]
    laplace = V.laplace
    nx = length(ctx.theta)
    aws = _laplace_workspace!(laplace, Float64, nx)
    scratch = _continuation_scratch(laplace, nx)
    f = lay.shaped
    nd = f === nothing ? nothing : _vf_nodes(f, lam, lay.members[1].Z)
    # The shaped family's prior-and-entropy Hessian is formed with the members'
    # curvature and held with their shares (`_vf_solve_held`); the Gaussian
    # family's is a diagonal in closed form.
    pe = _vf_prior_entropy(lay, lam, V.kappa, curvature || f === nothing, nd)
    B = pe.B
    g = copy(pe.g)
    GE = zeros(nx)
    gth = zeros(nx)
    u = zeros(lay.d)
    gu = zeros(lay.d)
    gm = [zeros(length(mem.plist)) for mem in lay.members]
    Hm = curvature ? [zeros(length(mem.plist), length(mem.plist)) for mem in lay.members] :
        Matrix{Float64}[]
    Lsh = f === nothing ? zeros(0, 0) : _vf_L(f, lam)
    mode = _laplace_hmode(laplace)
    gp = 0
    hp = 0
    invalid = (B=-Inf, g=g, gpe=pe.g, gm=gm, Hm=Hm, Hneg=pe.Hneg, GE=GE, gpasses=gp,
        hpasses=hp)
    for (mi, mem) in enumerate(lay.members)
        k = length(mem.chain)
        np = length(mem.plist)
        offs = laplace.units.offsets[U][mem.member]
        uc = zeros(k)
        Jt = zeros(k, np)
        r = zeros(k)
        T1 = curvature ? zeros(k, np) : zeros(0, 0)
        for kk in eachindex(mem.w)
            z = view(mem.Z, :, kk)
            Jn = if f === nothing
                _vf_tree_point!(uc, Jt, mem, lam, z)
                Jt
            else
                uc .= view(nd.U, :, kk)
                view(nd.J[kk], 1:k, :)
            end
            @inbounds for (a, c) in enumerate(mem.chain)
                u[c] = uc[a]
                gu[c] = 0.0
            end
            fill!(gth, 0.0)
            ll = try
                _continuation_member!(gth, laplace, U, mem.member, ctx.theta, ctx.Ls, ctx.dL,
                    ctx.positions, u, aws, scratch, true)
            catch err
                err isa InterruptException && rethrow()
                -Inf
            end
            gp += 1
            isfinite(ll) || return invalid
            _laplace_loading_chain!(gu, scratch.grad, laplace.spec, ctx.Ls, offs)
            @inbounds for (a, c) in enumerate(mem.chain)
                r[a] = gu[c]
            end
            wk = mem.w[kk]
            B += wk * ll
            GE .+= wk .* gth
            mul!(gm[mi], transpose(Jn), r, wk, 1.0)
            curvature || continue
            Jw = try
                _laplace_columns_jacobian(laplace, U, ctx.theta, ctx.Ls, u, mem.chain,
                    [mem.member], mode, nothing)
            catch err
                err isa InterruptException && rethrow()
                return invalid
            end
            hp += 1
            Hll = zeros(k, k)
            @inbounds for b in 1:k, a in 1:k
                Hll[a, b] = (Jw[mem.chain[a], b] + Jw[mem.chain[b], a]) / 2
            end
            all(isfinite, Hll) || return invalid
            # The engine's small products: BLAS's matrix product takes a
            # process-wide lock on every call (CLAUDE.md, "Julia threads").
            _ctsem_mul!(T1, Hll, Jn)
            _ctsem_mulTN!(Hm[mi], Jn, T1, wk, true)
            if f === nothing
                @inbounds for a in 1:k
                    t = mem.jld[a]
                    Hm[mi][t, t] += wk * r[a] * exp(lam[mem.plist[t]]) * z[a]
                end
            elseif f.bend
                Hm[mi] .+= wk .* ForwardDiff.hessian(_VFUAt(f, Vector(z), copy(r)), lam)
            else
                _vf_second!(Hm[mi], wk, f, lam, Lsh, z, r)
            end
        end
        @inbounds for (t, q) in enumerate(mem.plist)
            g[q] += gm[mi][t]
        end
        curvature && (Hm[mi] .= (Hm[mi] .+ transpose(Hm[mi])) ./ 2)
    end
    (isfinite(B) && all(isfinite, g)) || return invalid
    return (B=B, g=g, gpe=pe.g, gm=gm, Hm=Hm, Hneg=pe.Hneg, GE=GE, gpasses=gp,
        hpasses=hp)
end

"""The negated Hessian in the family's parameters on the unit's tree
(`CTSEMBlockMatrix` over `pblocks`): each member's negated share, scattered
into its blocks and their couplings, plus the prior and entropy terms'."""
function _vf_assemble(lay::CTSEMVLayout, shares::Vector{Matrix{Float64}}, Hneg)
    M = CTSEMBlockMatrix(Float64, lay.pblocks)
    for (mi, mem) in enumerate(lay.members)
        S = shares[mi]
        for (tb, qb) in enumerate(mem.plist), (ta, qa) in enumerate(mem.plist)
            ba = lay.pblock[qa]
            bb = lay.pblock[qb]
            if ba == bb
                M.diag[ba][lay.plocal[qa], lay.plocal[qb]] += S[ta, tb]
            else
                slot = findfirst(==(bb), lay.pblocks[ba].ancestors)
                slot === nothing && continue
                M.coupling[ba][slot][lay.plocal[qa], lay.plocal[qb]] += S[ta, tb]
            end
        end
    end
    if Hneg isa AbstractVector
        for q in eachindex(Hneg)
            M.diag[lay.pblock[q]][lay.plocal[q], lay.plocal[q]] += Hneg[q]
        end
    elseif !isempty(Hneg)
        M.diag[1] .+= Hneg
    end
    return M
end

"""A Newton step on `M` (the negated Hessian on the parameters' tree): the
block elimination, with a Levenberg shift where `M` is not definite.
`(step, shifted)`, `step === nothing` where no shift made it factor."""
function _vf_block_step(M::CTSEMBlockMatrix{Float64}, pblocks, g::Vector{Float64})
    F = _laplace_block_factor(M, pblocks)
    shifted = false
    if !F[1]
        top = 1.0
        for D in M.diag, i in axes(D, 1)
            top = max(top, abs(D[i, i]))
        end
        mu = 1e-8 * top
        for _ in 1:30
            Ms = _laplace_copy_blockmatrix(M)
            for D in Ms.diag, i in axes(D, 1)
                D[i, i] += mu
            end
            F = _laplace_block_factor(Ms, pblocks)
            F[1] && break
            mu *= 10
        end
        F[1] || return (nothing, true)
        shifted = true
    end
    return (_laplace_block_solve(F[3], F[4], pblocks, g), shifted)
end

# A held share's SR1 update from a step `s` and its term's gradient's change
# `y` (negated: the share is the negated Hessian), skipped where the update's
# denominator is small against the step and the residual.
function _vf_sr1!(S::Matrix{Float64}, s::Vector{Float64}, y::Vector{Float64})
    rr = y .- S * s
    den = dot(rr, s)
    abs(den) > 1e-8 * norm(s) * norm(rr) || return S
    S .+= (rr * transpose(rr)) ./ den
    return S
end

"""Where the shaped family's Newton matrix is not definite and its step has
nothing left to gain, the direction of its most negative curvature: a
saddle of the bound, which a Levenberg step does not leave. `nothing`
otherwise, and on the tree family."""
function _vf_saddle_direction(lay::CTSEMVLayout, M::CTSEMBlockMatrix{Float64})
    lay.shaped === nothing && return nothing
    E = _ctsem_symeig(Matrix(Symmetric(M.diag[1])))
    i = argmin(E.values)
    E.values[i] < -1e-6 || return nothing
    return E.vectors[:, i]
end

"""Whether unit `U`'s members can be evaluated at their prior mean: a unit that
cannot is beyond any approximation at these parameters (an explosive drift,
a variance below zero), said at the cost of one point a member."""
function _vf_finite_at_zero(V, U::Int, ctx::_VCtx)
    laplace = V.laplace
    nx = length(ctx.theta)
    aws = _laplace_workspace!(laplace, Float64, nx)
    scratch = _continuation_scratch(laplace, nx)
    u = zeros(V.layouts[U].d)
    for mem in V.layouts[U].members
        ll = try
            _continuation_member!(Float64[], laplace, U, mem.member, ctx.theta, ctx.Ls,
                nothing, ctx.positions, u, aws, scratch, false)
        catch err
            err isa InterruptException && rethrow()
            -Inf
        end
        isfinite(ll) || return false
    end
    return true
end

"""A start every node of which is finite: the approximation shrunk towards the
prior, on gradient passes, or `nothing`."""
function _vf_finite_start(V, U::Int, ctx::_VCtx)
    lay = V.layouts[U]
    _vf_finite_at_zero(V, U, ctx) || return nothing
    for s in 1:12
        lam = zeros(lay.p)
        ild = lay.shaped === nothing ? lay.ild : lay.shaped.ild
        lam[ild] .= -s * log(2.0)
        e = _vf_unit_terms(V, U, ctx, lam, false)
        isfinite(e.B) && return lam
    end
    return nothing
end

"""Newton on the bound's exact Hessian (the members' curvature at every node
of every point it evaluates), the step capped, Armijo backtracking;
converged where half the decrement is under `tol` on a definite matrix. A
shaped unit at a saddle leaves along its most negative curvature."""
function _vf_solve(V, U::Int, ctx::_VCtx, lam0::Vector{Float64})
    lay = V.layouts[U]
    lam = copy(lam0)
    gp = 0
    hp = 0
    e = _vf_unit_terms(V, U, ctx, lam, true)
    gp += e.gpasses
    hp += e.hpasses
    if !isfinite(e.B)
        start = _vf_finite_start(V, U, ctx)
        start === nothing && return (lam=lam, e=e, converged=false, gpasses=gp, hpasses=hp,
            shares=Matrix{Float64}[])
        lam = start
        e = _vf_unit_terms(V, U, ctx, lam, true)
        gp += e.gpasses
        hp += e.hpasses
        isfinite(e.B) || return (lam=lam, e=e, converged=false, gpasses=gp, hpasses=hp,
            shares=Matrix{Float64}[])
    end
    converged = lay.p == 0
    it = 0
    while !converged
        M = _vf_assemble(lay, [-H for H in e.Hm], e.Hneg)
        D, shifted = _vf_block_step(M, lay.pblocks, e.g)
        D === nothing && break
        dec = dot(e.g, D)
        escape = false
        if dec / 2 < V.tol
            if !shifted
                converged = true
                break
            end
            v = _vf_saddle_direction(lay, M)
            v === nothing && (converged = true; break)
            D = 0.5 .* (dot(v, e.g) < 0 ? -v : v)
            escape = true
        end
        it >= V.maxit && break
        it += 1
        big = maximum(abs, D)
        big > 5 && (D .*= 5 / big)
        slope = dot(e.g, D)
        alpha = 1.0
        accepted = false
        e1 = e
        lam1 = lam
        for _ in 1:30
            lam1 = lam .+ alpha .* D
            e1 = _vf_unit_terms(V, U, ctx, lam1, true)
            gp += e1.gpasses
            hp += e1.hpasses
            need = escape ? 1e-12 * max(1.0, abs(e.B)) : 1e-4 * alpha * slope
            if isfinite(e1.B) && e1.B >= e.B + need
                accepted = true
                break
            end
            alpha /= 2
        end
        if !accepted
            converged = escape || dec / 2 < 1e-6
            break
        end
        lam, e = lam1, e1
    end
    return (lam=lam, e=e, converged=converged, gpasses=gp, hpasses=hp,
        shares=_vf_shares(lay, e))
end

# The shares a solve hands to the held one: each member's negated curvature,
# and for the shaped family the prior-and-entropy term's last.
_vf_shares(lay::CTSEMVLayout, e) = lay.shaped === nothing ? [-H for H in e.Hm] :
    vcat([-H for H in e.Hm], isempty(e.Hneg) ? Matrix{Float64}[] : [copy(e.Hneg)])

"""The inner solve with the members' shares of the curvature held: from where
they were last formed, stepping on gradient passes alone, each share updated
by SR1 from its member's gradient's change (as `_LaplacePSR1` holds the
Laplace route's; damped BFGS on shares made definite formed the curvature
more often on gD3, 67% of unit-evaluations against 48%); the prior and
entropy terms exact for the Gaussian family, and for the shaped family a
share of their own, held and updated alike. The curvature is formed afresh (`_vf_solve`) where the start's
decrement is above `refresh`, a step had to be cut below a quarter, or
`quasi` steps did not converge."""
function _vf_solve_held(V, U::Int, ctx::_VCtx, lam0::Vector{Float64},
    shares0::Vector{Matrix{Float64}})
    lay = V.layouts[U]
    lam = copy(lam0)
    shares = [copy(S) for S in shares0]
    gp = 0
    e = _vf_unit_terms(V, U, ctx, lam, false)
    gp += e.gpasses
    if isfinite(e.B)
        it = 0
        while true
            shaped = lay.shaped !== nothing
            shaped && length(shares) != length(lay.members) + 1 && break
            M = _vf_assemble(lay, shares, shaped ? shares[end] : e.Hneg)
            D, shifted = _vf_block_step(M, lay.pblocks, e.g)
            D === nothing && break
            dec = dot(e.g, D)
            if dec / 2 < V.tol && !shifted
                return (lam=lam, e=e, converged=true, gpasses=gp, hpasses=0, shares=shares,
                    formed=false)
            end
            (it == 0 && dec / 2 > V.refresh) && break
            it >= V.quasi && break
            it += 1
            big = maximum(abs, D)
            big > 5 && (D .*= 5 / big)
            slope = dot(e.g, D)
            alpha = 1.0
            accepted = false
            lam1 = lam
            e1 = e
            for _ in 1:8
                lam1 = lam .+ alpha .* D
                e1 = _vf_unit_terms(V, U, ctx, lam1, false)
                gp += e1.gpasses
                if isfinite(e1.B) && e1.B >= e.B + 1e-4 * alpha * slope
                    accepted = true
                    break
                end
                alpha /= 2
            end
            (accepted && alpha >= 0.25) || break
            for (mi, mem) in enumerate(lay.members)
                _vf_sr1!(shares[mi], lam1[mem.plist] .- lam[mem.plist], e.gm[mi] .- e1.gm[mi])
            end
            shaped && _vf_sr1!(shares[end], lam1 .- lam, e.gpe .- e1.gpe)
            lam, e = lam1, e1
        end
    end
    r = _vf_solve(V, U, ctx, lam)
    return (lam=r.lam, e=r.e, converged=r.converged, gpasses=r.gpasses + gp,
        hpasses=r.hpasses, shares=r.shares, formed=true)
end

# ---- the objective ---------------------------------------------------------------------

"""
    CTSEMVariational

The variational objective over a Laplace objective's units (`laplace`, whose
spec, units, pool and hybrid hooks it uses): each unit's layout (its family
on its block tree), and the approximation and held curvature shares
accepted (`lam`, `shares`) and from the last evaluation (`*_trial`). See the
file's docstring.
"""
mutable struct CTSEMVariational{L} <: CTSEMOptimisable
    laplace::L
    layouts::Vector{CTSEMVLayout}
    lam::Vector{Vector{Float64}}
    shares::Vector{Vector{Matrix{Float64}}}
    lam_trial::Vector{Vector{Float64}}
    shares_trial::Vector{Vector{Matrix{Float64}}}
    accepted::Bool
    has_trial::Bool
    running::Bool
    tol::Float64
    kappa::Float64
    quasi::Int
    refresh::Float64
    maxit::Int
    unit_value::Vector{Float64}
    evaluations::Int
    gradient_passes::Int
    curvature_passes::Int
    curvature_formed::Int
    lbfgs::Union{Nothing,CTSEMLBFGSMemory}
    # Units whose family is chosen by need (`family = "auto"`, one block of two
    # or three coordinates): Gaussian, tried once with the bending family where
    # a run stops (`_ctsem_take_search_modes!`), kept bending where that gains
    # more than `need` nats.
    pending::Vector{Bool}
    need::Float64
    # The layouts and choices a new fit starts from (`ctsem_variational_reset!`).
    start::Tuple{Vector{CTSEMVLayout},Vector{Bool}}
end

"""
    ctsem_variational(laplace; family = "auto", rule = "auto", seed = 1,
        tol = 1e-8, kappa = 1e-3, quasi = 10, refresh = 5.0, maxit = 50,
        need = 0.005)

The variational objective over `laplace`'s units. `family = "auto"`: the
Gaussian family on each unit's block tree (`_vf_tree_layout`), except a unit
of one block and one coordinate, which is skewed (product Gauss-Hermite, 7
nodes), and a unit of one block and two or three coordinates, which is tried
with the bending family (5 nodes a coordinate) where a run stops and kept
bending where its bound gains more than `need` nats. "gauss" for the Gaussian
family everywhere; "skew" or "bend" for the shaped family on units of one
block and at most three coordinates (Gaussian elsewhere). `rule` = "product"
or "rot" fixes the Gaussian members' rule.
"""
function ctsem_variational(laplace::CTSEMLaplaceObjective; family::AbstractString="auto",
    rule::AbstractString="auto", seed::Integer=1, tol::Real=1e-8, kappa::Real=1e-3,
    quasi::Integer=10, refresh::Real=5.0, maxit::Integer=50, need::Real=0.005)
    family in ("auto", "gauss", "skew", "bend") ||
        throw(ArgumentError("family: auto, gauss, skew or bend"))
    rule in ("auto", "product", "rot") || throw(ArgumentError("rule: auto, product or rot"))
    _laplace_ensure_pool!(laplace)
    n = length(laplace.units.dims)
    layouts = CTSEMVLayout[]
    pending = fill(false, n)
    for U in 1:n
        d = laplace.units.dims[U]
        oneblock = length(laplace.units.blocks[U]) == 1 && 1 <= d <= 3
        shaped = oneblock && (family in ("skew", "bend") || (family == "auto" && d == 1))
        lay = if shaped
            f = CTSEMVFamily(d, true; bend=family == "bend")
            _vf_shaped_layout(laplace, U, f, d == 1 ? 7 : 5)
        else
            _vf_tree_layout(laplace, U, Symbol(rule), Int(seed))
        end
        pending[U] = family == "auto" && oneblock && d >= 2
        push!(layouts, lay)
    end
    return CTSEMVariational(laplace, layouts, [zeros(l.p) for l in layouts],
        [Matrix{Float64}[] for _ in layouts], [zeros(l.p) for l in layouts],
        [Matrix{Float64}[] for _ in layouts], false, false, false, Float64(tol),
        Float64(kappa), Int(quasi), Float64(refresh), Int(maxit), fill(NaN, n), 0, 0, 0, 0,
        nothing, pending, Float64(need), (copy(layouts), copy(pending)))
end

function _vf_unit!(V::CTSEMVariational, ctx::_VCtx, U::Int, unit::Vector{Float64},
    GE::Matrix{Float64}, gps::Vector{Int}, hps::Vector{Int}, formed::Vector{Bool})
    r = if isempty(V.shares[U])
        s = _vf_solve(V, U, ctx, V.lam[U])
        (lam=s.lam, e=s.e, gpasses=s.gpasses, hpasses=s.hpasses, shares=s.shares,
            formed=true)
    else
        _vf_solve_held(V, U, ctx, V.lam[U], V.shares[U])
    end
    gps[U] = r.gpasses
    hps[U] = r.hpasses
    formed[U] = r.formed
    isfinite(r.e.B) || return false
    V.lam_trial[U] = r.lam
    V.shares_trial[U] = r.shares
    unit[U] = r.e.B
    GE[:, U] .= r.e.GE
    return true
end

function _vf_commit!(V::CTSEMVariational)
    V.has_trial || return nothing
    for U in eachindex(V.lam)
        V.lam[U] = copy(V.lam_trial[U])
        V.shares[U] = V.shares_trial[U]
    end
    V.accepted = true
    return nothing
end

"""
    ctsem_variational_evaluate(V, values; gradient = true, units = false)

The bound summed over units with the prior, and its gradient; per unit the
bounds (`unit_loglik`) and, with `units`, each unit's gradient (its weight
applied, the prior not). Every unit is solved from the accepted
approximations; the result is the trial until the optimiser accepts it
(`_ctsem_anchor_here!`). Outside an optimisation nothing is kept, except that
a first evaluation's approximations become the accepted ones.
"""
function ctsem_variational_evaluate(V::CTSEMVariational, values::AbstractVector;
    gradient::Bool=true, units::Bool=false)
    laplace = V.laplace
    npar = length(values)
    nunits = length(V.layouts)
    V.evaluations += 1
    bad = (value=-Inf, gradient=zeros(npar), unit_loglik=fill(NaN, nunits),
        unit_gradient=zeros(npar, units ? nunits : 0), converged=false)
    ctx = try
        _vf_ctx(laplace, values)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    ctx === nothing && return bad
    nx = length(ctx.theta)
    unit = fill(NaN, nunits)
    GE = zeros(nx, nunits)
    gps = zeros(Int, nunits)
    hps = zeros(Int, nunits)
    formed = fill(false, nunits)
    ok = _laplace_parallel(laplace, 1:nunits) do U
        return _vf_unit!(V, ctx, U, unit, GE, gps, hps, formed)
    end
    V.gradient_passes += sum(gps)
    V.curvature_passes += sum(hps)
    V.curvature_formed += count(formed)
    (ok && all(isfinite, unit)) || return bad
    V.has_trial = true
    V.unit_value = unit
    (!V.running && !V.accepted) && _vf_commit!(V)
    wts = [_laplace_unit_weight(laplace, U) for U in 1:nunits]
    prior = _ctsem_log_prior(laplace.objective, ctx.theta)
    gp = zeros(nx)
    _ctsem_log_prior_gradient!(gp, laplace.objective, ctx.theta)
    G = GE .* transpose(wts)
    return (value=dot(wts, unit) + prior,
        gradient=(vec(sum(G; dims=2)) .+ gp)[1:npar], unit_loglik=unit,
        unit_gradient=units ? G[1:npar, :] : zeros(npar, 0), converged=true)
end

"""Counts (evaluations, gradient and curvature passes a member's node, how
many times a unit's curvature was formed) and each unit's shape: its
coordinates, its members' nodes in all, its family and rule."""
function ctsem_variational_info(V::CTSEMVariational)
    lay = V.layouts
    return (evaluations=V.evaluations, gradient_passes=V.gradient_passes,
        curvature_passes=V.curvature_passes, curvature_formed=V.curvature_formed,
        units=length(lay), dims=[l.d for l in lay], parameters=[l.p for l in lay],
        blocks=[length(l.pblocks) for l in lay],
        nodes=[sum(length(m.w) for m in l.members; init=0) for l in lay],
        skewed=[l.shaped !== nothing && any(l.shaped.skew) for l in lay],
        bending=[l.shaped !== nothing && l.shaped.bend for l in lay],
        pending=copy(V.pending), rules=[String(l.rule) for l in lay])
end

"""
    _ctsem_take_search_modes!(V::CTSEMVariational, x, found)

Where a run stops by its own rule at `x`: each unit whose family is chosen by
need and still Gaussian fits the bending family there from its Gaussian
approximation (whose parameters are the bending family's first ones, so the
bound starts where the Gaussian's is), and keeps it where its bound gains
more than `V.need` nats. The change to the objective at `x` sends the run on
from there with its curvature pairs, as the Laplace route's search for other
maxima does, up to the run's mode rounds. A unit that does not gain stays a
candidate at the run's later stops: tried once only, gD3's units were tried
where L-BFGS first stopped, far from the optimum, and 25 of 60 that bend at
the optimum stayed Gaussian (0.94 nats). On gD3 every unit gained at the
optimum (median 0.05 nats); on gD1 none did (at most 0.005), and its fit
stays at the Gaussian family's cost.
"""
function _ctsem_take_search_modes!(V::CTSEMVariational, x::AbstractVector, found)
    todo = findall(V.pending)
    isempty(todo) && return (change=NaN, units=0, joined=0)
    laplace = V.laplace
    ctx = _vf_ctx(laplace, x)
    gains = fill(-Inf, length(V.layouts))
    cand = Vector{Any}(undef, length(V.layouts))
    _laplace_parallel(laplace, todo) do U
        lay0 = V.layouts[U]
        base = isempty(V.shares[U]) ? _vf_solve(V, U, ctx, V.lam[U]) :
            _vf_solve_held(V, U, ctx, V.lam[U], V.shares[U])
        isfinite(base.e.B) || return true
        f = CTSEMVFamily(lay0.d, true; bend=true)
        lay = _vf_shaped_layout(laplace, U, f, 5)
        lam0 = vcat(base.lam, zeros(f.p - lay0.p))
        # Solved on the bending layout without touching the unit's own.
        trial = _VFOne(V, U, lay)
        s = _vf_solve(trial, U, ctx, lam0)
        isfinite(s.e.B) || return true
        gains[U] = s.e.B - base.e.B
        cand[U] = (lay=lay, lam=s.lam, shares=s.shares)
        return true
    end
    change = 0.0
    units = 0
    for U in todo
        gains[U] > V.need || continue
        V.pending[U] = false
        c = cand[U]
        V.layouts[U] = c.lay
        V.lam[U] = copy(c.lam)
        V.shares[U] = c.shares
        V.lam_trial[U] = copy(c.lam)
        V.shares_trial[U] = c.shares
        change += _laplace_unit_weight(laplace, U) * gains[U]
        units += 1
    end
    units == 0 && return (change=0.0, units=0, joined=0)
    return (change=change, units=units, joined=units)
end

_ctsem_search_note(::CTSEMVariational, taken) = @sprintf(
    "%d unit(s) took the bending family where the run stopped, moving the objective by %+.3g; the run goes on",
    taken.units, taken.change)

# A view of the objective with one unit's layout replaced, for solving a
# candidate family without touching the unit's own (`_vf_unit_terms` and the
# solves read `layouts`, `laplace` and the settings alone).
struct _VFOne{L}
    laplace::L
    layouts::Dict{Int,CTSEMVLayout}
    tol::Float64
    kappa::Float64
    quasi::Int
    refresh::Float64
    maxit::Int
end
_VFOne(V::CTSEMVariational, U::Int, lay::CTSEMVLayout) = _VFOne(V.laplace,
    Dict(U => lay), V.tol, V.kappa, V.quasi, V.refresh, V.maxit)

"""
    ctsem_variational_reset!(V)

Forgets every unit's approximation, held curvature and family choice, so the
next optimisation starts each unit from the prior on its starting family: a
new fit on an objective R had cached.
"""
function ctsem_variational_reset!(V::CTSEMVariational)
    layouts, pending = V.start
    V.layouts = copy(layouts)
    V.pending = copy(pending)
    for U in eachindex(V.layouts)
        V.lam[U] = zeros(V.layouts[U].p)
        V.shares[U] = Matrix{Float64}[]
        V.lam_trial[U] = zeros(V.layouts[U].p)
        V.shares_trial[U] = Matrix{Float64}[]
    end
    V.accepted = false
    V.has_trial = false
    V.running = false
    V.lbfgs = nothing
    V.evaluations = V.gradient_passes = V.curvature_passes = V.curvature_formed = 0
    return nothing
end

"""
    ctsem_variational_choose!(V, values)

The family choice a run makes where it stops (`_ctsem_take_search_modes!`),
made at `values`; for a caller outside an optimisation.
"""
ctsem_variational_choose!(V::CTSEMVariational, values::AbstractVector) =
    _ctsem_take_search_modes!(V, collect(Float64, values), nothing)

# ---- the optimisable protocol -----------------------------------------------------------

function ctsem_evaluate(V::CTSEMVariational, values::AbstractVector; gradient::Bool=true,
    contributions::Bool=false, gradient_method=:adjoint)
    r = ctsem_variational_evaluate(V, values; gradient=gradient)
    value = r.converged ? r.value : NaN
    contributions || return (value=value, gradient=r.gradient)
    return (value=value, gradient=r.gradient,
        subject_loglik=_vf_subject_values(V, r.unit_loglik))
end

# Each subject's share of its unit's term, spread evenly over the unit's
# members as the Laplace route spreads its own (`ctsem_laplace_evaluate`): a
# unit with a study level does not decompose over its subjects.
function _vf_subject_values(V::CTSEMVariational, unit::Vector{Float64})
    out = fill(NaN, length(V.laplace.objective.subject_objectives))
    for (U, members) in enumerate(V.laplace.units.members)
        share = unit[U] / length(members)
        for i in members
            out[i] = share
        end
    end
    return out
end

"""The fit's Hessian: central differences of the exact gradient, every unit
solved from the accepted approximations at each point."""
function ctsem_hessian(V::CTSEMVariational, values::AbstractVector; chunk::Integer=0,
    progress=nothing, step::Real=1e-4)
    x = collect(Float64, values)
    n = length(x)
    H = zeros(n, n)
    # The differences leave the trial where it was: the optimiser's next
    # accepted point is not one of them.
    kept = (copy(V.lam_trial), copy(V.shares_trial), V.has_trial, copy(V.unit_value))
    for j in 1:n
        xp = copy(x)
        xp[j] += step
        xm = copy(x)
        xm[j] -= step
        gpl = ctsem_variational_evaluate(V, xp).gradient
        gmi = ctsem_variational_evaluate(V, xm).gradient
        H[:, j] .= (gpl .- gmi) ./ (2 * step)
    end
    V.lam_trial, V.shares_trial, V.has_trial, V.unit_value = kept
    return (H .+ transpose(H)) ./ 2
end

"""Each unit's gradient (its weight applied) with its share of the prior, as
the Laplace route's rows carry theirs."""
function ctsem_subject_gradients(V::CTSEMVariational, values::AbstractVector;
    prior::Bool=true, sink::Union{Nothing,AbstractVector}=nothing)
    n = length(values)
    r = ctsem_variational_evaluate(V, values; units=true)
    sink === nothing || copyto!(sink, r.unit_loglik)
    isfinite(r.value) || return (value=r.value, scores=fill(NaN, length(V.layouts), n))
    S = Matrix(transpose(r.unit_gradient))
    if prior
        theta = _laplace_extend(V.laplace, collect(Float64, values))
        Sx = zeros(size(S, 1), length(theta))
        Sx[:, 1:n] .= S
        _ctsem_prior_scores!(Sx, V.laplace.objective, theta,
            _laplace_unit_of_subject(V.laplace))
        S = Sx[:, 1:n]
    end
    return (value=r.value, scores=S)
end

_ctsem_optimise_label(::CTSEMVariational) = "variational"
_ctsem_finish_curvature(::CTSEMVariational) = :chord
_ctsem_hessian_cost(::CTSEMVariational, npar::Integer) = 2 * Int(npar)
_ctsem_gradient_values(V::CTSEMVariational) = _ctsem_gradient_values(V.laplace)
_ctsem_lbfgs_held(V::CTSEMVariational) = V.lbfgs
_ctsem_lbfgs_hold!(V::CTSEMVariational, memory::CTSEMLBFGSMemory) =
    (V.lbfgs = memory; nothing)
_ctsem_params(V::CTSEMVariational) = _ctsem_params(V.laplace)
_ctsem_saturated_for(V::CTSEMVariational, minimizer) =
    _laplace_saturated_parameters(V.laplace, minimizer)
_ctsem_optimise_result_extra(V::CTSEMVariational, final, log) =
    (variational=ctsem_variational_info(V),)

# Within an optimisation each evaluation starts from the approximations the
# last accepted point left; the per-step callback accepts the last
# evaluation's.
function _ctsem_anchor_start!(V::CTSEMVariational, x::AbstractVector;
    resume_units=nothing, resume_values=nothing)
    V.accepted || ctsem_variational_evaluate(V, x)
    V.running = true
    return nothing
end

function _ctsem_anchor_here!(V::CTSEMVariational)
    V.running || return nothing
    _vf_commit!(V)
    return nothing
end

_ctsem_anchor_finish!(V::CTSEMVariational, x) = (V.running = false; nothing)

function ctsem_laplace_clear_anchors!(V::CTSEMVariational)
    V.running = false
    return nothing
end

# The Laplace route's entry points R calls on the fit's objective, for a
# variational fit's (`intoverpop = 'variational'`).
ctsem_laplace_optimize(V::CTSEMVariational, start::AbstractVector;
    nested_gradient::Bool=false, kwargs...) = ctsem_optimize(V, start; kwargs...)

ctsem_laplace_hessian(V::CTSEMVariational, values::AbstractVector; step::Real=1e-4,
    kwargs...) = ctsem_hessian(V, values; step=step)

ctsem_laplace_evaluate(V::CTSEMVariational, values::AbstractVector; gradient::Bool=true,
    kwargs...) = ctsem_variational_evaluate(V, values; gradient=gradient)

# The readers of the effects, the plans made before a fit (substeps,
# relinearisation), generation, and the random-effects check whose integral a
# variational fit reports, take the Laplace objective the bound is built on:
# a variational fit reports each unit's effects at their conditional modes
# there, at the fit's estimate. What reads the objective's value per unit or
# per row -- LOO's unit terms, the score tests' Hessian columns and row
# scores, the samplers -- has no method here yet and fails by name.
for f in (:ctsem_laplace_subject_values, :ctsem_laplace_sequential,
    :ctsem_laplace_sequential_mixture, :ctsem_laplace_population_deviations,
    :ctsem_laplace_deviation_layout, :ctsem_parameter_layout,
    :ctsem_state_dependent_cells, :ctsem_laplace_conditioning,
    :ctsem_auto_substeps, :ctsem_relin_plan, :ctsem_effect_information,
    :ctsem_generate, :ctsem_generate_draws, :ctsem_generate_states, :ctsem_kalman,
    :ctsem_laplace_coordinate_dimension, :ctsem_laplace_effect_draws,
    :ctsem_laplace_popcov, :ctsem_laplace_posterior_nodes, :ctsem_laplace_unit_dims,
    :ctsem_parameter_matrices, :ctsem_state_dimension, :ctsem_state_layout,
    :ctsem_state_trajectories, :ctsem_laplace_check)
    @eval $f(V::CTSEMVariational, args...; kwargs...) = $f(V.laplace, args...; kwargs...)
end
