# Forward mode through the filter: the routes the hand-coded second order
# (`tangent_*.jl`) replaced, which the engine no longer runs, kept here as the
# oracles it is tested against. Each is the engine's own generic code at a dual
# or a trace jet, at the engine's tags and widths (`_ctsem_hessian_dual`,
# `_laplace_curvature_dual`, `_LaplaceSeedInner`), so the tests' compile cache
# (test/cache) holds them. Exact in the arithmetic, so they agree with the hand
# sweeps to rounding: finite differences could not tell a wrong second-order
# term of 1e-7 from their own error.

const _DO = ContinuousTimeSEM
# The engine's, where the session has none of its own (an R session's Julia).
if !isdefined(@__MODULE__, :ForwardDiff)
    const ForwardDiff = _DO.ForwardDiff
end

"""The Hessian: forward mode over the adjoint gradient at the Hessian's dual."""
function _oracle_hessian(objective, values; hybrid_mode=_DO._HYBRID_PROFILE)
    x = collect(Float64, values)
    H = _DO._ctsem_width_jacobian(_DO._CTSEMAdjointGradientOf(objective, Int(hybrid_mode)),
        x, _DO._CTSEM_HESSIAN_WIDTHS[])
    return (H .+ transpose(H)) ./ 2
end

"""Each data row's score along `directions`: the filter at the Hessian's dual,
the deviations it integrates off, as `ctsem_row_scores` has them."""
function _oracle_row_scores(objective, values, directions)
    theta = collect(Float64, values)
    S = _DO._ctsem_hessian_dual(objective)
    W = ForwardDiff.npartials(S)
    K = size(directions, 2)
    sp = objective.params
    subjects = objective.subject_objectives
    rows = [size(s.data, 2) for s in subjects]
    out = fill(NaN, sum(rows; init=0), K)
    ws = _DO._init_continuous_ekf_workspace(S, sp)
    ws.hyb.mode = _DO._HYBRID_OFF
    n, m = ws.state_dim, ws.manifest_dim
    nall = length(sp.mutables)
    for c0 in 1:W:K
        cols = c0:min(c0 + W - 1, K)
        x = S[S(theta[j], ForwardDiff.Partials{W,Float64}(ntuple(a -> a <= length(cols) ?
            Float64(directions[j, cols[a]]) : 0.0, Val(W)))) for j in eachindex(theta)]
        offset = 0
        for (i, sub) in enumerate(subjects)
            nobs = rows[i]
            trace = _DO.CTSEMKalmanTrace(S, n, m, nobs, 1, nall)
            trace.offset = 0
            trace.current_subject = 1
            ll = _DO._extended_kalman_filter_continuous!(ws, x, sub.data,
                collect(sub.timesteps), sp, sub.tdpreds,
                _DO._ctsem_tipred_vector(sub.tipreds, x), i, sub.max_timestep, trace,
                nothing, nothing, sub.relin)
            if isfinite(ForwardDiff.value(ll))
                for r in 1:nobs, (a, col) in enumerate(cols)
                    out[offset + r, col] = ForwardDiff.partials(trace.llrow[r], a)
                end
            end
            offset += nobs
        end
    end
    return out
end

"""`d2 g_U / du du`, the unit objective's Hessian in its coordinates: forward
mode over its whole gradient at the curvature's widths, knowing nothing of the
blocks, so a reference for the blocked assembly and the hand-coded sweeps."""
function _oracle_unit_hessian(laplace, U, values, Ls, u)
    _DO._laplace_ensure_pool!(laplace)
    length(u) == 0 && return zeros(0, 0)
    inner_of = function (uu)
        S = eltype(uu)
        ws = _DO._laplace_workspace!(laplace, S, length(values))
        vs = _DO._laplace_scratch_vector!(laplace, S, length(values), :curv_vs)
        copyto!(vs, values)
        Lss = _DO._laplace_scratch_levels!(laplace, S, Ls, :curv_levels)
        return _DO._laplace_unit_objective_gradient(laplace, U, vs, Lss, uu, ws,
            _DO._laplace_scratch_vector!(laplace, S, length(uu), :curv_inner)).gradient
    end
    H = _DO._ctsem_jacobian(inner_of, collect(u), _DO._CTSEM_CURVATURE_WIDTHS[])
    return (H .+ transpose(H)) ./ 2
end

"""The unit curvature `I - (J + J') / 2`, dense, `J` the members' gradient's
Jacobian in the unit's coordinates at the curvature's dual
(`_LaplaceColumnsGradient`), whose passes read and fill `hook` (a
`CTSEMHybridCross`, the gated floor's cross curvature): the route the hand-coded
sweeps replaced."""
function _oracle_unit_curvature(laplace, U, values, Ls, u; hook=nothing,
    mode=_DO._HYBRID_PROFILE)
    n = length(u)
    block_of = _DO._LaplaceColumnsGradient(laplace, Int(U), values, Ls,
        collect(Float64, u), collect(1:n), collect(eachindex(laplace.units.members[U])),
        Int(mode), hook)
    width = _DO._ctsem_dual_width(n, _DO._CTSEM_CURVATURE_WIDTHS[])
    nx = max(n, width)
    xb = zeros(nx)
    xb[1:n] .= u
    Jw = zeros(n, nx)
    _DO._ctsem_width_jacobian!(Jw, block_of, xb, n, width; whole=true,
        store=laplace.workspaces[_DO._laplace_slot()])
    J = Jw[:, 1:n]
    return [Float64(i == j) for i in 1:n, j in 1:n] .- (J .+ transpose(J)) ./ 2
end

"""The joint curvature at fixed deviations, minus the Hessian of the unit's
joint objective (`_hybrid_joint_objective_gradient`) in its joint coordinates
`z`, dense, by forward mode over its gradient: a reference for
`_hybrid_joint_curvature` sharing none of its assembly."""
function _oracle_joint_curvature(laplace, U, values, Ls, z, layout)
    inner_of = function (zz)
        S = eltype(zz)
        ws = _DO._laplace_workspace!(laplace, S, length(values))
        return _DO._hybrid_joint_objective_gradient(laplace, U, S.(values),
            [S.(L) for L in Ls], zz, ws, layout).gradient
    end
    H = _DO._ctsem_jacobian(inner_of, collect(Float64, z), _DO._CTSEM_CURVATURE_WIDTHS[])
    return -(H .+ transpose(H)) ./ 2
end

"""`_laplace_unit_seeded_gradient`: each member's gradient (`d0`) and its
derivative along `d1` (`d1c`), the member's pass at a width-one dual."""
function _oracle_seeded_gradient(laplace, U, values, Ls, u, members, d1)
    spec = laplace.spec
    units = laplace.units
    unitmembers = units.members[U]
    npar = length(values)
    nm = length(members)
    d0 = zeros(npar, nm)
    d1c = zeros(npar, nm)
    seed = ForwardDiff.Dual{_DO._LaplaceSeedInner}(0.0, 1.0)
    S1 = typeof(seed)
    for (c, m) in enumerate(members)
        aws = _DO._laplace_workspace!(laplace, S1, npar)
        gradient = zeros(S1, npar)
        shifted = _DO._laplace_member_values!(zeros(npar), values, spec, Ls, u,
            units.offsets[U][m])
        x = S1[shifted[t] + seed * d1[t] for t in 1:npar]
        _DO._laplace_member_pass_finite!(gradient,
            laplace.objective.subject_objectives[unitmembers[m]], aws, x) ||
            return (ok=false, d0=d0, d1c=d1c)
        for t in 1:npar
            d0[t, c] = ForwardDiff.value(gradient[t])
            d1c[t, c] = ForwardDiff.partials(gradient[t])[1]
        end
    end
    return (ok=true, d0=d0, d1c=d1c)
end

"""`_laplace_unit_jet_gradient`: each member's gradient (`d0`) and the sum of
its second directional derivatives along `D`'s first `ndirs` columns (`dS`),
the member's pass at a trace jet of width `W`."""
function _oracle_jet_gradient(laplace, U, values, Ls, u, members, D, ndirs, W)
    spec = laplace.spec
    units = laplace.units
    unitmembers = units.members[U]
    npar = length(values)
    nm = length(members)
    J = _DO.CTSEMTraceJet{W}
    d0 = zeros(npar, nm)
    dS = zeros(npar, nm)
    for (c, m) in enumerate(members)
        aws = _DO._laplace_workspace!(laplace, J, npar)
        gradient = zeros(J, npar)
        shifted = _DO._laplace_member_values!(zeros(npar), values, spec, Ls, u,
            units.offsets[U][m])
        x = J[J(shifted[t], ntuple(j -> j <= ndirs ? D[t, j] : 0.0, Val(W)), 0.0)
            for t in 1:npar]
        _DO._laplace_member_pass_finite!(gradient,
            laplace.objective.subject_objectives[unitmembers[m]], aws, x) ||
            return (ok=false, d0=d0, dS=dS)
        for t in 1:npar
            d0[t, c] = gradient[t].v
            dS[t, c] = gradient[t].s
        end
    end
    return (ok=true, d0=d0, dS=dS)
end

"""`_laplace_mode_cross`: `d/dz_a` of the unit objective's raw gradient at
`z`, row `a`, by a width-one dual in each coordinate."""
function _oracle_mode_cross(laplace, U, theta, Ls, dL, positions, z, layout)
    D = length(z)
    npar = length(theta)
    H = zeros(D, npar)
    seed = ForwardDiff.Dual{_DO._LaplaceSeedInner}(0.0, 1.0)
    for a in 1:D
        zz = [z[b] + seed * (b == a ? 1.0 : 0.0) for b in 1:D]
        r = _DO._laplace_unit_theta_gradient(laplace, U, theta, Ls, dL, positions, zz,
            layout)
        r === nothing && return nothing
        for t in 1:npar
            H[a, t] = ForwardDiff.partials(r.theta[t])[1]
        end
    end
    return H
end

"""`_gated_directional`: the gated rule's gradients and their derivatives along
`beta`, at a width-one dual."""
function _oracle_gated_directional(G, z, beta)
    seed = ForwardDiff.Dual{_DO._LaplaceSeedInner}(0.0, 1.0)
    zz = [z[a] + seed * beta[a] for a in eachindex(z)]
    r = G(zz)
    r === nothing && return nothing
    return (theta=[ForwardDiff.partials(x)[1] for x in r.theta],
        z=[ForwardDiff.partials(x)[1] for x in r.z])
end

"""`_laplace_gated_jet_pass!`'s second parts along `chunk`'s eigendirections of
`E`, at a trace jet of width `W`: `(theta, z)`, before its sign."""
function _oracle_gated_jet(G, z, E, chunk, W)
    D = length(z)
    zz = [_DO.CTSEMTraceJet{W}(z[a], ntuple(j -> j <= length(chunk) ?
        sqrt(abs(E.values[chunk[j]])) * E.vectors[a, chunk[j]] : 0.0, Val(W)), 0.0)
        for a in 1:D]
    r = G(zz)
    r === nothing && return nothing
    return (theta=[x.s for x in r.theta], z=[x.s for x in r.z])
end

"""The plug-in Kalman trace of subject `i` at the duals `x`, the filter at
their type with the deviations it integrates off: what `_tan_laplace_trace!`
writes by hand."""
function _oracle_trace!(trace, laplace, i, x)
    sp = laplace.objective.params
    so = laplace.objective.subject_objectives[i]
    ws = _DO._init_continuous_ekf_workspace(eltype(x), sp)
    ws.hyb.mode = _DO._HYBRID_OFF
    ll = _DO._extended_kalman_filter_continuous!(ws, x, so.data, collect(so.timesteps),
        sp, so.tdpreds, _DO._ctsem_tipred_vector(so.tipreds, x), i, so.max_timestep,
        trace, nothing, nothing, so.relin)
    return _DO._laplace_finite(ll)
end
