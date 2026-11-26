"""
The variational objective: each unit's random effects integrated by a lower
bound on its marginal likelihood, fitted at every evaluation, in place of the
Laplace term (review/DESIGN-unit-bound-2026-11-26.md).

A unit's approximation maps a standard normal `z` to its standardised
effects `u` on the hybrid's coordinates (the deviations the filter integrates
stay integrated exactly):

    y_j = sinh((asinh(z_j) + e_j) / d_j)      skewed axes; y_j = z_j otherwise
    u   = m + L y                             L lower triangular

and on fixed nodes `z_k` with weights `w_k` the bound is

    B_U = sum_k w_k h(u_k) + E[log |du/dz|] + dim / 2 - kappa |shape|^2 / 2

with `h = log L_U - |u|^2 / 2`. The second term, the approximation's entropy
less a constant, does not involve the model, and for these families is a sum
of functions of one coordinate each; it is integrated on a fine
one-dimensional rule (`_vf_entropy`), so the family cannot gain by
stretching its map between the nodes, and the rule's error enters through
the likelihood term alone. The approximation maximises `B_U`, so the
bound's gradient in the parameters is the partial with the nodes held: first
derivatives of the members' likelihoods at the nodes. Any approximation gives
a lower bound up to the rule's error; the shape prior (`kappa`, on the skew
and tail parameters) keeps a unit whose posterior is the prior cut off at a
wall from creeping towards a limit at infinity, and does not depend on the
parameters.

The inner solve is Newton on the bound's Hessian in the family's parameters.
Its curvature term needs the unit's curvature at each node, the dear part, so
each unit holds the Hessian from where it was last formed and steps by it on
gradient passes alone, updating it by BFGS; the curvature is formed afresh
only where the held matrix does not serve (`_vf_solve_held`).

Within an optimisation (`running`) each evaluation starts every unit from the
approximation the last accepted point left, and the optimiser's per-step
callback (`_ctsem_anchor_here!`) makes the last evaluation's approximations
the accepted ones, as the Laplace route moves its anchors. Outside one, an
evaluation starts from the accepted approximations and keeps nothing, so a
value is the parameters' given the held state.

Rules: product Gauss-Hermite for the skewed and bending families (7 nodes at
one coordinate, 5 a coordinate at two or three); the rotated `2d + 1` rule
for the Gaussian family beyond (unbiased over rotations, not a guaranteed
bound unit by unit). The skewed family on any rule that misses the
transform's higher terms exploits it: on 3 nodes a dimension, the rotated
rule or the symmetric degree-5 rule its "bound" came out above the integral.
A rule with negative weights is no bound at all: on the cut rule (each axis
by Gauss-Hermite, the rest by the degree-5 rule's off-axis points), one of
gA1's three-coordinate units moved its approximation's mean 100 prior sds
out, where the negative-weight nodes' likelihood was -2e7, for a "bound" of
+1e6.
"""

export ctsem_variational, ctsem_variational_info, ctsem_variational_value,
    ctsem_variational_reset!

# ---- families -------------------------------------------------------------------

"""A unit's family: `d` coordinates, which axes are skewed, whether it bends,
and where each parameter sits in its vector (`m`, `log diag L`, the strict
lower triangle of `L` by column, `(e_j, log d_j)` for each skewed axis, then
`(b_ji, c_ji)` for i < j where it bends: `x_j = y_j exp(sum_i c_ji s(y_i)) +
sum_i b_ji (s(y_i)^2 - 1)`, `s(y) = 5 tanh(y / 5)`, `u = m + L x`)."""
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

# The bending family's map at any element type, for its derivatives: u and
# log |du/dz|. Small and free of the filter, so differentiated by ForwardDiff
# on the callable structs below, at units of two or three coordinates.
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
            lj += abs(a) + log1p(exp(-2 * abs(a))) - log(2.0) - dl - log1p(z[j]^2) / 2
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

struct _VFPhiAt
    f::CTSEMVFamily
    z::Vector{Float64}
    r::Vector{Float64}
end
(m::_VFPhiAt)(lam) = (q = _vf_map_generic(m.f, lam, m.z); dot(m.r, q[1]) + q[2])

# r'u alone, for the likelihood term's second derivatives through the map.
struct _VFUAt
    f::CTSEMVFamily
    z::Vector{Float64}
    r::Vector{Float64}
end
(m::_VFUAt)(lam) = dot(m.r, _vf_map_generic(m.f, lam, m.z)[1])

# log cosh, smooth at any element type.
_vf_lcosh(a) = abs(a) > 20 ? abs(a) - log(2.0) + log1p(exp(-2 * abs(a))) : log(cosh(a))

@inline _vf_logcosh(a::Float64) = abs(a) + log1p(exp(-2 * abs(a))) - log(2.0)

# One axis's transform at node value z: y and its derivatives in (e, log d),
# and the log Jacobian's.
function _vf_axis(f::CTSEMVFamily, lam::Vector{Float64}, z::Float64, j::Int)
    if !f.skew[j]
        return (y=z, lj=0.0, ye=0.0, yd=0.0, yee=0.0, yed=0.0, ydd=0.0, le=0.0, ld=0.0,
            lee=0.0, led=0.0, ldd=0.0)
    end
    dl = lam[f.idl[j]]
    del = exp(dl)
    a = (asinh(z) + lam[f.ie[j]]) / del
    sh = sinh(a)
    ch = cosh(a)
    th = tanh(a)
    se2 = 1 / ch^2
    return (y=sh, lj=_vf_logcosh(a) - dl - log1p(z^2) / 2, ye=ch / del, yd=-a * ch,
        yee=sh / del^2, yed=-(a * sh + ch) / del, ydd=a * ch + a^2 * sh, le=th / del,
        ld=-a * th - 1, lee=se2 / del^2, led=(-a * se2 - th) / del, ldd=a * th + a^2 * se2)
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
            ax = _vf_axis(f, lam, Z[j, k], j)
            y[j] = ax.y
            l += ax.lj
            if f.skew[j]
                for i in j:d
                    Jk[i, f.ie[j]] = L[i, j] * ax.ye
                    Jk[i, f.idl[j]] = L[i, j] * ax.yd
                end
                Jk[d + 1, f.ie[j]] = ax.le
                Jk[d + 1, f.idl[j]] = ax.ld
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

"""The bound, its gradient in the family's parameters and, where `hH` (h's
Hessian in u at each node) is given, its Hessian there. A non-finite part
makes the point invalid."""
function _vf_parts(f::CTSEMVFamily, lam::Vector{Float64}, Z::Matrix{Float64},
    w::Vector{Float64}, nd, h::Vector{Float64}, hg::Matrix{Float64}, hH,
    kappa::Float64)
    d, K = size(Z)
    p = f.p
    invalid = (B=-Inf, g=zeros(p), H=zeros(0, 0))
    all(isfinite, h) || return invalid
    B = d / 2
    g = zeros(p)
    want = hH !== nothing
    H = want ? zeros(p, p) : zeros(0, 0)
    for j in f.shape
        B -= kappa * lam[j]^2 / 2
        g[j] -= kappa * lam[j]
        want && (H[j, j] -= kappa)
    end
    L = want ? _vf_L(f, lam) : zeros(0, 0)
    ent = _VFEntropyAt(f)
    B += ent(lam)
    if any(f.skew) || f.bend
        g .+= ForwardDiff.gradient(ent, lam)
        want && (H .+= ForwardDiff.hessian(ent, lam))
    else
        g[f.ild] .+= 1.0
    end
    r = zeros(d)
    for k in 1:K
        B += w[k] * h[k]
        Jk = nd.J[k]
        for i in 1:d
            r[i] = hg[i, k]
        end
        for a in 1:p
            s = 0.0
            for i in 1:d
                s += Jk[i, a] * r[i]
            end
            g[a] += w[k] * s
        end
        if want
            Ju = view(Jk, 1:d, :)
            H .+= w[k] .* (transpose(Ju) * hH[k] * Ju)
            if f.bend
                H .+= w[k] .* ForwardDiff.hessian(_VFUAt(f, Z[:, k], copy(r)), lam)
            else
                _vf_second!(H, w[k], f, lam, L, view(Z, :, k), r)
            end
        end
    end
    want && (H .= (H .+ transpose(H)) ./ 2)
    (isfinite(B) && all(isfinite, g) && all(isfinite, H)) || return invalid
    return (B=B, g=g, H=H)
end

function _vf_pd(Mx::AbstractMatrix)
    E = eigen(Symmetric(Matrix(Mx)))
    top = max(maximum(abs, E.values), 1.0)
    vals = max.(abs.(E.values), 1e-8 * top)
    return Matrix(Symmetric(E.vectors * Diagonal(vals) * transpose(E.vectors)))
end

# ---- rules ------------------------------------------------------------------------

# Gauss-Hermite for the standard normal, n nodes (Golub-Welsch).
function _vf_gh(n::Integer)
    n == 1 && return [0.0], [1.0]
    J = zeros(n, n)
    for i in 1:(n - 1)
        J[i, i + 1] = J[i + 1, i] = sqrt(i)
    end
    E = eigen(Symmetric(J))
    return E.values, E.vectors[1, :] .^ 2
end

const _VF_FINE = _vf_gh(96)

"""E[log |du/dz|] under the standard normal, at any element type: log |L|,
each skewed axis's term, and each bending coefficient times the mean of
`s(y_i)`, every one an integral over one coordinate, on 96-node
Gauss-Hermite."""
function _vf_entropy(f::CTSEMVFamily, lam::AbstractVector)
    z1, w1 = _VF_FINE
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

function _vf_rule_product(d::Integer, n::Integer)
    d == 0 && return zeros(0, 1), [1.0]
    z1, w1 = _vf_gh(n)
    K = Int(n)^d
    Z = zeros(d, K)
    w = ones(K)
    k = 0
    for idx in CartesianIndices(ntuple(_ -> Int(n), Int(d)))
        k += 1
        for j in 1:d
            Z[j, k] = z1[idx[j]]
            w[k] *= w1[idx[j]]
        end
    end
    return Z, w ./ sum(w)
end

"""The cut rule: each axis by `n1`-node Gauss-Hermite (exact for any sum of
one-coordinate functions to degree 2 n1 - 1), what is left by the degree-5
symmetric rule's off-axis points; exact to degree 5, `2d^2 + (n1 - 1)d + 1`
nodes. The centre's and the axis points' at `+-s` weights are negative, on
points where the remainder of an additive function is zero."""
function _vf_rule_cut(d::Integer; n1::Integer=7)
    d == 0 && return zeros(0, 1), [1.0]
    z1, w1 = _vf_gh(n1)
    k0 = argmin(abs.(z1))
    pts = Vector{Vector{Float64}}([zeros(d)])
    w = [0.0]
    for j in 1:d, k in eachindex(z1)
        if k == k0
            w[1] += w1[k]
            continue
        end
        z = zeros(d)
        z[j] = z1[k]
        push!(pts, z)
        push!(w, w1[k])
    end
    w[1] -= d - 1
    if d > 1
        s = sqrt((d + 2) / 2)
        W2 = 1 / (d + 2)^2
        for i in 1:d, sg in (-1, 1)
            z = zeros(d)
            z[i] = sg * s
            push!(pts, z)
            push!(w, -2 * (d - 1) * W2)
        end
        w[1] += 2 * d * (d - 1) * W2
        for i in 1:d, j in (i + 1):d, a in (-1, 1), b in (-1, 1)
            z = zeros(d)
            z[i] = a * s
            z[j] = b * s
            push!(pts, z)
            push!(w, W2)
        end
    end
    return reduce(hcat, pts), w
end

"""The centre and +-sqrt(d + 2) along `R` random rotations, drawn once for
the unit (`seed`) and kept: exact to degree 3, unbiased to degree 5 over
rotations, `2Rd + 1` nodes, positive weights."""
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

# ---- a unit at its nodes ----------------------------------------------------------

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

"""h = log L_U - |u|^2 / 2 at the points `U`, its gradient in u, its gradient
in theta (each times nothing: the rule's weights apply later), and where
`curv` its Hessian in u. A point the engine cannot evaluate, or throws at, is
an invalid point."""
function _vf_eval(laplace::CTSEMLaplaceObjective, ctx::_VCtx, U::Int,
    Un::Matrix{Float64}, curv::Bool)
    d, K = size(Un)
    nx = length(ctx.theta)
    aws = _laplace_workspace!(laplace, Float64, nx)
    scratch = _continuation_scratch(laplace, nx)
    members = laplace.units.members[U]
    h = fill(-Inf, K)
    hg = zeros(d, K)
    GT = zeros(nx, K)
    gu = zeros(d)
    gth = zeros(nx)
    u = zeros(d)
    for k in 1:K
        for i in 1:d
            u[i] = Un[i, k]
        end
        fill!(gth, 0.0)
        fill!(gu, 0.0)
        hk = 0.0
        ok = true
        try
            for mm in eachindex(members)
                ll = _continuation_member!(gth, laplace, U, mm, ctx.theta, ctx.Ls, ctx.dL,
                    ctx.positions, u, aws, scratch, true)
                if !isfinite(ll)
                    ok = false
                    break
                end
                _laplace_loading_chain!(gu, scratch.grad, laplace.spec, ctx.Ls,
                    laplace.units.offsets[U][mm])
                hk += ll
            end
        catch err
            err isa InterruptException && rethrow()
            ok = false
        end
        ok || continue
        h[k] = hk - sum(abs2, u) / 2
        for i in 1:d
            hg[i, k] = gu[i] - u[i]
        end
        for t in 1:nx
            GT[t, k] = gth[t]
        end
    end
    hH = nothing
    if curv
        blocks = laplace.units.blocks[U]
        hH = Vector{Matrix{Float64}}(undef, K)
        for k in 1:K
            hH[k] = zeros(d, d)
            isfinite(h[k]) || continue
            try
                M = _laplace_unit_curvature(laplace, U, ctx.theta, ctx.Ls, Un[:, k])
                Md = length(M.diag) == 1 ? copy(M.diag[1]) : _laplace_block_dense(M, blocks, d)
                hH[k] = -(Md .+ transpose(Md)) ./ 2
            catch err
                err isa InterruptException && rethrow()
                h[k] = -Inf
            end
        end
    end
    return (h=h, hg=hg, hH=hH, GT=GT)
end

# ---- the inner solve ----------------------------------------------------------------

function _vf_step(Mx::Matrix{Float64}, g::Vector{Float64})
    E = eigen(Symmetric(Mx))
    top = max(maximum(abs, E.values), 1.0)
    vals = max.(abs.(E.values), 1e-8 * top)
    D = E.vectors * ((transpose(E.vectors) * g) ./ vals)
    return D, dot(g, D), minimum(E.values), E
end

"""Newton on the bound's exact Hessian (the curvature at every node it
evaluates), saddle-free, the step capped, Armijo backtracking; converged where
half the decrement is under `tol` on a definite Hessian, and at a saddle it
leaves along the positive-curvature direction."""
function _vf_solve(laplace, ctx::_VCtx, U::Int, f::CTSEMVFamily, lam0::Vector{Float64},
    Z::Matrix{Float64}, w::Vector{Float64}, tol::Float64, maxit::Int, kappa::Float64)
    K = size(Z, 2)
    lam = copy(lam0)
    hp = 0
    nd = _vf_nodes(f, lam, Z)
    ev = _vf_eval(laplace, ctx, U, nd.U, true)
    hp += K
    pt = _vf_parts(f, lam, Z, w, nd, ev.h, ev.hg, ev.hH, kappa)
    shrink = 0
    while !isfinite(pt.B) && shrink < 30 && f.d > 0
        shrink += 1
        lam = zeros(f.p)
        lam[f.ild] .= -shrink * log(2.0)
        nd = _vf_nodes(f, lam, Z)
        ev = _vf_eval(laplace, ctx, U, nd.U, true)
        hp += K
        pt = _vf_parts(f, lam, Z, w, nd, ev.h, ev.hg, ev.hH, kappa)
    end
    isfinite(pt.B) || return (lam=lam, B=-Inf, ev=ev, converged=false, iterations=0,
        gpasses=0, hpasses=hp, H=pt.H)
    it = 0
    converged = f.p == 0
    while !converged
        D, dec, lmin, E = _vf_step(-pt.H, pt.g)
        escape = false
        if dec / 2 < tol
            if lmin > -1e-6
                converged = true
                break
            end
            v = E.vectors[:, argmin(E.values)]
            D = 0.5 .* (dot(v, pt.g) < 0 ? -v : v)
            escape = true
        end
        it >= maxit && break
        it += 1
        big = maximum(abs, D)
        big > 5 && (D .*= 5 / big)
        slope = dot(pt.g, D)
        alpha = 1.0
        accepted = false
        lam1 = lam
        nd1 = nd
        ev1 = ev
        pt1 = pt
        for _ in 1:30
            lam1 = lam .+ alpha .* D
            nd1 = _vf_nodes(f, lam1, Z)
            ev1 = _vf_eval(laplace, ctx, U, nd1.U, true)
            hp += K
            pt1 = _vf_parts(f, lam1, Z, w, nd1, ev1.h, ev1.hg, ev1.hH, kappa)
            need = escape ? 1e-12 * max(1.0, abs(pt.B)) : 1e-4 * alpha * slope
            if isfinite(pt1.B) && pt1.B >= pt.B + need
                accepted = true
                break
            end
            alpha /= 2
        end
        if !accepted
            converged = escape || dec / 2 < 1e-6
            break
        end
        lam, nd, ev, pt = lam1, nd1, ev1, pt1
    end
    return (lam=lam, B=pt.B, ev=ev, converged=converged, iterations=it, gpasses=0,
        hpasses=hp, H=pt.H)
end

"""The inner solve with the curvature held: the negated Hessian `M0` from
where it was last formed steps the unit on gradient passes alone and is
updated by BFGS; the curvature is formed afresh (`_vf_solve`) where the start's
decrement on `M0` is above `refresh`, a step had to be cut below a quarter, or
`quasi` steps did not converge. Returns the solve and the matrix to hold."""
function _vf_solve_held(laplace, ctx::_VCtx, U::Int, f::CTSEMVFamily,
    lam0::Vector{Float64}, Z::Matrix{Float64}, w::Vector{Float64}, M0::Matrix{Float64},
    tol::Float64, maxit::Int, kappa::Float64, quasi::Int, refresh::Float64)
    K = size(Z, 2)
    lam = copy(lam0)
    gp = 0
    nd = _vf_nodes(f, lam, Z)
    ev = _vf_eval(laplace, ctx, U, nd.U, false)
    gp += K
    pt = _vf_parts(f, lam, Z, w, nd, ev.h, ev.hg, nothing, kappa)
    M = copy(M0)
    if isfinite(pt.B)
        it = 0
        while true
            F = cholesky(Symmetric(M); check=false)
            issuccess(F) || break
            D = F \ pt.g
            dec = dot(pt.g, D)
            if dec / 2 < tol
                return (lam=lam, B=pt.B, ev=ev, converged=true, iterations=it, gpasses=gp,
                    hpasses=0, M=M, formed=false)
            end
            (it == 0 && dec / 2 > refresh) && break
            it >= quasi && break
            it += 1
            big = maximum(abs, D)
            big > 5 && (D .*= 5 / big)
            slope = dot(pt.g, D)
            alpha = 1.0
            accepted = false
            lam1 = lam
            nd1 = nd
            ev1 = ev
            pt1 = pt
            for _ in 1:8
                lam1 = lam .+ alpha .* D
                nd1 = _vf_nodes(f, lam1, Z)
                ev1 = _vf_eval(laplace, ctx, U, nd1.U, false)
                gp += K
                pt1 = _vf_parts(f, lam1, Z, w, nd1, ev1.h, ev1.hg, nothing, kappa)
                if isfinite(pt1.B) && pt1.B >= pt.B + 1e-4 * alpha * slope
                    accepted = true
                    break
                end
                alpha /= 2
            end
            (accepted && alpha >= 0.25) || break
            s = lam1 .- lam
            y = pt.g .- pt1.g
            sy = dot(s, y)
            if sy > 1e-12 * dot(y, y)
                Ms = M * s
                M = M .- (Ms * transpose(Ms)) ./ dot(s, Ms) .+ (y * transpose(y)) ./ sy
                M = (M .+ transpose(M)) ./ 2
            end
            lam, nd, ev, pt = lam1, nd1, ev1, pt1
        end
    end
    r = _vf_solve(laplace, ctx, U, f, lam, Z, w, tol, maxit, kappa)
    Mn = isempty(r.H) ? M : _vf_pd(-r.H)
    return (lam=r.lam, B=r.B, ev=r.ev, converged=r.converged, iterations=r.iterations,
        gpasses=r.gpasses + gp, hpasses=r.hpasses, M=Mn, formed=true)
end

# ---- the objective -------------------------------------------------------------------

"""
    CTSEMVariational

The variational objective over a Laplace objective's units (`laplace`, whose
spec, units, pool and hybrid hooks it uses): each unit's family, rule, and
the approximation and held curvature accepted (`lam`, `M`) and from the last
evaluation (`lam_trial`, `M_trial`). See the file's docstring.
"""
mutable struct CTSEMVariational{L} <: CTSEMOptimisable
    laplace::L
    fams::Vector{CTSEMVFamily}
    Z::Vector{Matrix{Float64}}
    w::Vector{Vector{Float64}}
    rules::Vector{Symbol}
    lam::Vector{Vector{Float64}}
    M::Vector{Matrix{Float64}}
    lam_trial::Vector{Vector{Float64}}
    M_trial::Vector{Matrix{Float64}}
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
end

"""
    ctsem_variational(laplace; family = "auto", rule = "auto", seed = 1,
        tol = 1e-8, kappa = 1e-3, quasi = 10, refresh = 5.0)

The variational objective over `laplace`'s units. `family = "auto"`: the
skewed family at one coordinate on product Gauss-Hermite (7 nodes), the
bending family at two or three on product Gauss-Hermite (5 nodes a
coordinate), the Gaussian family on the rotated `2d + 1` rule beyond; "bend",
"skew" or "gauss" for every unit; `nodes` sets the product rule's nodes a
coordinate. Every default rule has positive weights:
on a rule with negative weights (the cut rule, by name) the bound is not a
bound, and an approximation can climb it without limit by moving the
negative-weight nodes to where the likelihood is very low. `rule` = "auto", "product", "cut" or "rot" overrides the rule.
"""
function ctsem_variational(laplace::CTSEMLaplaceObjective; family::AbstractString="auto",
    rule::AbstractString="auto", seed::Integer=1, tol::Real=1e-8, kappa::Real=1e-3,
    quasi::Integer=10, refresh::Real=5.0, maxit::Integer=50, nodes::Integer=0)
    _laplace_ensure_pool!(laplace)
    dims = laplace.units.dims
    n = length(dims)
    fams = CTSEMVFamily[]
    Zs = Matrix{Float64}[]
    ws = Vector{Float64}[]
    rules = Symbol[]
    for U in 1:n
        d = dims[U]
        skew = family in ("skew", "bend") ? true : family == "gauss" ? false : d <= 3
        bend = family == "bend" || (family == "auto" && 2 <= d <= 3)
        push!(fams, CTSEMVFamily(d, skew; bend=bend))
        r = rule != "auto" ? Symbol(rule) : (d <= 1 || skew) ? :product : :rot
        n1 = nodes > 0 ? Int(nodes) : !skew ? 3 : d == 1 ? 7 : 5
        Z, w = r === :product ? _vf_rule_product(d, n1) :
            r === :cut ? _vf_rule_cut(d; n1=skew ? 7 : 3) :
            r === :rot ? _vf_rule_rot(d, 1, Int(seed) * 100003 + U) :
            throw(ArgumentError("rule: auto, product, cut or rot"))
        push!(Zs, Z)
        push!(ws, w)
        push!(rules, r)
    end
    return CTSEMVariational(laplace, fams, Zs, ws, rules, [zeros(f.p) for f in fams],
        [zeros(0, 0) for _ in fams], [zeros(f.p) for f in fams], [zeros(0, 0) for _ in fams],
        false, false, false, Float64(tol), Float64(kappa), Int(quasi), Float64(refresh),
        Int(maxit), fill(NaN, n), 0, 0, 0, 0, nothing)
end

function _vf_unit!(V::CTSEMVariational, ctx::_VCtx, U::Int, unit::Vector{Float64},
    GE::Matrix{Float64}, gps::Vector{Int}, hps::Vector{Int}, formed::Vector{Bool})
    f = V.fams[U]
    Z = V.Z[U]
    w = V.w[U]
    laplace = V.laplace
    if f.p == 0
        ev = _vf_eval(laplace, ctx, U, Z, false)
        isfinite(ev.h[1]) || return false
        unit[U] = ev.h[1]
        GE[:, U] .= view(ev.GT, :, 1)
        gps[U] = 1
        return true
    end
    r = if isempty(V.M[U])
        s = _vf_solve(laplace, ctx, U, f, V.lam[U], Z, w, V.tol, V.maxit, V.kappa)
        (lam=s.lam, B=s.B, ev=s.ev, gpasses=s.gpasses, hpasses=s.hpasses,
            M=isempty(s.H) ? Matrix{Float64}(I, f.p, f.p) : _vf_pd(-s.H), formed=true)
    else
        _vf_solve_held(laplace, ctx, U, f, V.lam[U], Z, w, V.M[U], V.tol, V.maxit,
            V.kappa, V.quasi, V.refresh)
    end
    gps[U] = r.gpasses
    hps[U] = r.hpasses
    formed[U] = r.formed
    isfinite(r.B) || return false
    V.lam_trial[U] = r.lam
    V.M_trial[U] = r.M
    unit[U] = r.B
    GE[:, U] .= r.ev.GT * w
    return true
end

function _vf_commit!(V::CTSEMVariational)
    V.has_trial || return nothing
    for U in eachindex(V.lam)
        V.lam[U] = copy(V.lam_trial[U])
        V.M[U] = V.M_trial[U]
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
    nunits = length(V.fams)
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

"""
    ctsem_variational_value(V, values; draws = 2048, df = 4, seed = 11)

The log likelihood at `values` by importance sampling on each unit's accepted
approximation, its normal draws replaced by t (`df`) draws for heavier tails:
what a fit reports, beside the bound it climbed. Returns the total with the
prior, per unit values and standard errors, and each unit's effective draws.
"""
function ctsem_variational_value(V::CTSEMVariational, values::AbstractVector;
    draws::Integer=2048, df::Real=4, seed::Integer=11)
    laplace = V.laplace
    ctx = _vf_ctx(laplace, values)
    nunits = length(V.fams)
    val = fill(NaN, nunits)
    se = fill(NaN, nunits)
    ess = fill(NaN, nunits)
    _laplace_parallel(laplace, 1:nunits) do U
        return _vf_value_unit!(V, ctx, U, Int(draws), Float64(df), Int(seed), val, se, ess)
    end
    wts = [_laplace_unit_weight(laplace, U) for U in 1:nunits]
    prior = _ctsem_log_prior(laplace.objective, ctx.theta)
    return (value=dot(wts, val) + prior, loglik=sum(val), loglik_weighted=dot(wts, val),
        unit_loglik=val, subject_loglik=_vf_subject_values(V, val), unit_se=se, ess=ess,
        se=sqrt(sum(abs2, wts .* se)))
end

function _vf_lt(z::Vector{Float64}, df::Float64)
    d = length(z)
    c = _vf_lgamma((df + d) / 2) - _vf_lgamma(df / 2) - d / 2 * log(df * pi)
    return c - (df + d) / 2 * log1p(sum(abs2, z) / df)
end

function _vf_lgamma(x::Float64)
    r = 0.0
    while x > 1
        x -= 1
        r += log(x)
    end
    return r + (isapprox(x, 0.5) ? log(sqrt(pi)) : 0.0)
end

function _vf_value_unit!(V::CTSEMVariational, ctx::_VCtx, U::Int, draws::Int, df::Float64,
    seed::Int, val, se, ess)
    f = V.fams[U]
    d = f.d
    if d == 0
        ev = _vf_eval(V.laplace, ctx, U, zeros(0, 1), false)
        val[U] = ev.h[1]
        se[U] = 0.0
        ess[U] = 1.0
        return true
    end
    rng = Random.Xoshiro(seed + 7919 * U)
    half = draws ÷ 2
    Zt = zeros(d, 2 * half)
    for k in 1:half
        gsc = sqrt(sum(abs2, randn(rng, round(Int, df))) / df)
        z = randn(rng, d) ./ gsc
        Zt[:, k] .= z
        Zt[:, half + k] .= -z
    end
    nd = _vf_nodes(f, V.lam[U], Zt; jacobian=false)
    ev = _vf_eval(V.laplace, ctx, U, nd.U, false)
    lw = zeros(2 * half)
    for k in 1:(2 * half)
        z = Zt[:, k]
        # log N(u; 0, I) L_U(u) - log q(u), q the family pushed from a t:
        # h already holds -|u|^2 / 2; the normal's constant cancels the t's
        # standard part only partly, so both are written out.
        lw[k] = ev.h[k] - d / 2 * log(2pi) + nd.lj[k] - _vf_lt(z, df)
    end
    pk = maximum(lw)
    isfinite(pk) || return false
    p = exp.(lw .- pk)
    m = sum(p) / length(p)
    val[U] = pk + log(m)
    se[U] = sqrt(sum(abs2, p .- m) / (length(p) - 1)) / sqrt(length(p)) / m
    ess[U] = sum(p)^2 / sum(abs2, p)
    return true
end

"""Counts: evaluations, gradient and curvature passes a unit, and how many
times a unit's curvature was formed."""
ctsem_variational_info(V::CTSEMVariational) = (evaluations=V.evaluations,
    gradient_passes=V.gradient_passes, curvature_passes=V.curvature_passes,
    curvature_formed=V.curvature_formed, units=length(V.fams),
    nodes=[size(Z, 2) for Z in V.Z], dims=[f.d for f in V.fams],
    skewed=[any(f.skew) for f in V.fams], rules=String.(V.rules))

# ---- the optimisable protocol -------------------------------------------------------

function ctsem_evaluate(V::CTSEMVariational, values::AbstractVector; gradient::Bool=true,
    contributions::Bool=false, gradient_method=:adjoint)
    r = ctsem_variational_evaluate(V, values; gradient=gradient)
    value = r.converged ? r.value : NaN
    contributions || return (value=value, gradient=r.gradient)
    return (value=value, gradient=r.gradient,
        subject_loglik=_vf_subject_values(V, r.unit_loglik))
end

function _vf_subject_values(V::CTSEMVariational, unit::Vector{Float64})
    nsub = length(V.laplace.objective.subject_objectives)
    out = fill(NaN, nsub)
    for (U, members) in enumerate(V.laplace.units.members)
        length(members) == 1 && (out[members[1]] = unit[U])
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
    kept = (copy(V.lam_trial), copy(V.M_trial), V.has_trial, copy(V.unit_value))
    for j in 1:n
        xp = copy(x)
        xp[j] += step
        xm = copy(x)
        xm[j] -= step
        gpl = ctsem_variational_evaluate(V, xp).gradient
        gmi = ctsem_variational_evaluate(V, xm).gradient
        H[:, j] .= (gpl .- gmi) ./ (2 * step)
    end
    V.lam_trial, V.M_trial, V.has_trial, V.unit_value = kept
    return (H .+ transpose(H)) ./ 2
end

"""Each unit's gradient (its weight applied) with its share of the prior, as
the Laplace route's rows carry theirs."""
function ctsem_subject_gradients(V::CTSEMVariational, values::AbstractVector;
    prior::Bool=true, sink::Union{Nothing,AbstractVector}=nothing)
    n = length(values)
    r = ctsem_variational_evaluate(V, values; units=true)
    sink === nothing || copyto!(sink, r.unit_loglik)
    isfinite(r.value) || return (value=r.value, scores=fill(NaN, length(V.fams), n))
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

"""
    ctsem_variational_reset!(V)

Forgets every unit's approximation and held curvature, so the next
optimisation starts each unit from the prior: a new fit on an objective R
had cached.
"""
function ctsem_variational_reset!(V::CTSEMVariational)
    for U in eachindex(V.fams)
        V.lam[U] = zeros(V.fams[U].p)
        V.M[U] = zeros(0, 0)
        V.lam_trial[U] = zeros(V.fams[U].p)
        V.M_trial[U] = zeros(0, 0)
    end
    V.accepted = false
    V.has_trial = false
    V.running = false
    V.lbfgs = nothing
    V.evaluations = V.gradient_passes = V.curvature_passes = V.curvature_formed = 0
    return nothing
end

# The readers of the effects, the plans made before a fit (substeps,
# relinearisation) and generation take the Laplace objective the bound is
# built on: a variational fit reports each unit's effects at their conditional
# modes there, at the fit's estimate. `ctsem_laplace_evaluate` on the fit's
# objective is its value. What reads the objective's value per unit or per
# row -- LOO's unit terms, the score tests' Hessian columns and row scores,
# the samplers -- has no method here yet and fails by name.
for f in (:ctsem_laplace_subject_values, :ctsem_laplace_sequential,
    :ctsem_laplace_sequential_mixture, :ctsem_laplace_population_deviations,
    :ctsem_laplace_deviation_layout, :ctsem_parameter_layout,
    :ctsem_state_dependent_cells, :ctsem_laplace_conditioning,
    :ctsem_auto_substeps, :ctsem_relin_plan, :ctsem_effect_information,
    :ctsem_generate, :ctsem_generate_draws, :ctsem_generate_states, :ctsem_kalman,
    :ctsem_laplace_coordinate_dimension, :ctsem_laplace_effect_draws,
    :ctsem_laplace_popcov, :ctsem_laplace_posterior_nodes, :ctsem_laplace_unit_dims,
    :ctsem_parameter_matrices, :ctsem_state_dimension, :ctsem_state_layout,
    :ctsem_state_trajectories)
    @eval $f(V::CTSEMVariational, args...; kwargs...) = $f(V.laplace, args...; kwargs...)
end

ctsem_laplace_evaluate(V::CTSEMVariational, values::AbstractVector; gradient::Bool=true,
    kwargs...) = ctsem_variational_evaluate(V, values; gradient=gradient)
