# Predicted and updated states and observations of a Laplace fit, each
# conditional on the subject's own record up to that row only.
#
# `ctsem_kalman` on a Laplace fit filters every subject at its random effects'
# modes, solved from the whole record, so its predicted and updated passes see
# the subject's later data through those effects. Here every row's prediction
# is the fitted model's own: the subject's Laplace posterior given its rows so
# far -- the inner mode solved on the truncated record as the objective solves
# it, probes and all, but started from the mode one row earlier (the first row
# from the origin), and the curvature there -- with the
# deviations its filter integrates (`CTSEMHybrid`) at their conditional
# posterior given the same rows. The mean is the filter's at the mode, and
# the uncertainty of the effects is added to the filter's own by the delta
# method, the Jacobian taken by one dual pass seeded in the effects. Each row's
# log predictive density is the difference of the subject's Laplace terms with
# and without it, so a subject's rows sum to its term in the objective.
#
# A unit with no coordinate, whose effects the filter integrates whole, needs
# none of that: its model is linear and Gaussian in the deviations, so one pass
# gives every row (`_laplace_sequential_filtered!`), where the general path
# filters every prefix again and costs O(T^2) a subject.

struct _LaplaceSequentialTag end

@inline function _laplace_sequential_seed(::Type{ForwardDiff.Dual{G,V,N}}, v,
    a::Integer) where {G,V,N}
    return ForwardDiff.Dual{G,V,N}(V(v), ForwardDiff.Partials{N,V}(
        ntuple(j -> j == a ? one(V) : zero(V), Val(N))))
end

"""A truncated copy of a subject's record, rows `1:t`, in place on `so`."""
function _laplace_truncate!(so, full, t::Integer)
    so.data = full.data[:, 1:t]
    so.timesteps = full.timesteps[1:t]
    so.tdpreds = full.tdpreds[:, 1:t]
    so.max_timestep = full.max_timestep isa AbstractVector ?
        full.max_timestep[1:t] : full.max_timestep
    return so
end

"""
What a pass over a prefix of a unit's record overwrites: its retained mode and
inner-solve flags, and (in `full`) the member's record that
`_laplace_truncate!` cuts. `_laplace_unit_restore!` puts both back, in a
`finally`, so a unit is as the objective left it whatever the pass did.
"""
function _laplace_unit_saved(laplace, U::Integer, so)
    return (mode=copy(laplace.modes[U]), iterations=laplace.inner_iterations[U],
        gradient=laplace.inner_gradient[U], converged=laplace.inner_converged[U],
        hrepaired=laplace.hessian_repaired[U], mrepaired=laplace.mode_repaired[U],
        floored=laplace.logdet_floored[U],
        full=(data=so.data, timesteps=so.timesteps, tdpreds=so.tdpreds,
            max_timestep=so.max_timestep))
end

function _laplace_unit_restore!(laplace, U::Integer, so, saved)
    full = saved.full
    so.data = full.data
    so.timesteps = full.timesteps
    so.tdpreds = full.tdpreds
    so.max_timestep = full.max_timestep
    laplace.modes[U] = saved.mode
    laplace.inner_iterations[U] = saved.iterations
    laplace.inner_gradient[U] = saved.gradient
    laplace.inner_converged[U] = saved.converged
    laplace.hessian_repaired[U] = saved.hrepaired
    laplace.mode_repaired[U] = saved.mrepaired
    laplace.logdet_floored[U] = saved.floored
    return nothing
end

"""Inverse of a symmetric curvature by its eigenvalues, each floored at
`1e-8`: a covariance even where the mode is not a strict maximum."""
function _laplace_sequential_inverse(M::AbstractMatrix)
    d = size(M, 1)
    d == 0 && return zeros(Float64, 0, 0)
    E = _ctsem_symeig(Matrix{Float64}(M))
    out = zeros(Float64, d, d)
    @inbounds for k in 1:d
        w = 1 / max(E.values[k], 1e-8)
        for j in 1:d, i in 1:d
            out[i, j] += E.vectors[i, k] * w * E.vectors[j, k]
        end
    end
    return out
end

"""
    ctsem_laplace_sequential(laplace, values)

For every data row, the state and observation predicted from the subject's
rows before it (leading index 1) and updated with it (2), each with the
subject's random effects at their Laplace posterior given those rows alone, as
the top of this file describes; `llrow` each row's log predictive density.
Single-level fits only: in a nested unit a person's truncated record moves the
effects shared with the other members.
"""
function ctsem_laplace_sequential(laplace::CTSEMLaplaceObjective, values::AbstractVector;
    compact::Bool=false)
    theta = _laplace_extend(laplace, collect(Float64, values))
    units = laplace.units
    all(length(m) == 1 for m in units.members) || throw(ArgumentError(
        "sequential prediction needs every subject to be its own unit: a single " *
        "level of random effects"))
    Ls = _laplace_popchols(theta, laplace.spec)
    sp = laplace.objective.params
    subjects = laplace.objective.subject_objectives
    ws0 = _init_continuous_ekf_workspace(Float64, sp)
    n = ws0.state_dim
    m = ws0.manifest_dim
    rows = [size(s.data, 2) for s in subjects]
    first_row = cumsum(vcat(0, rows[1:(end - 1)]))
    nrows = sum(rows)
    # The effects' own posterior at each row, on the raw scale: the deviation
    # of each varying parameter and its covariance, given the rows before it
    # (1) and with it (2) -- what an augmented fit's carrier states are.
    k = nrandomeffects(laplace.spec.levels[1])
    out = (pitlower=fill(NaN, 0, 0), pitupper=fill(NaN, 0, 0),
        eta=fill(NaN, 2, nrows, n), etacov=fill(NaN, 2, nrows, n, n),
        y=fill(NaN, 2, nrows, m), ycov=fill(NaN, 2, nrows, m, m),
        llrow=fill(NaN, nrows), subject=zeros(Int, nrows),
        effect=fill(NaN, 2, nrows, max(k, 1)),
        effectcov=fill(NaN, 2, nrows, max(k, 1), max(k, 1)))
    for (i, r) in enumerate(rows)
        out.subject[(first_row[i] + 1):(first_row[i] + r)] .= i
    end
    hyb = sp.hybrid
    d = maximum(units.dims; init=0)
    W = max(d + (_hybrid_on(hyb) ? _hybrid_items(hyb) : 0), 1)
    S = ForwardDiff.Dual{_LaplaceSequentialTag,Float64,W}
    _laplace_ensure_pool!(laplace)
    _laplace_parallel(laplace, eachindex(units.members)) do U
        _laplace_sequential_unit!(out, laplace, U, theta, Ls,
            first_row[units.members[U][1]], S)
        return true
    end
    return compact ? _ctsem_kalman_compact(out) : out
end

function _laplace_sequential_unit!(out, laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, offset::Integer,
    ::Type{S}) where {S}
    laplace.units.dims[U] == 0 &&
        return _laplace_sequential_filtered!(out, laplace, U, theta, Ls, offset, S)
    spec = laplace.spec
    units = laplace.units
    i = units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    saved = _laplace_unit_saved(laplace, U, so)
    full = saved.full
    T = size(full.data, 2)
    sp = laplace.objective.params
    hyb = sp.hybrid
    hybrid = _hybrid_on(hyb)
    lv = spec.levels[1]
    B = hybrid ? _hybrid_dev_map(hyb, lv, theta) : zeros(0, 0)
    blocks = units.blocks[U]
    d = length(saved.mode)
    npar = length(theta)
    W = ForwardDiff.npartials(S)
    aws = _laplace_workspace!(laplace, Float64, npar)
    ws = _init_continuous_ekf_workspace(S, sp)
    ws.hyb.mode = _HYBRID_OFF
    n = ws.state_dim
    m = ws.manifest_dim
    nall = length(sp.mutables)
    prevterm = 0.0
    try
        for t in 0:T
            # The posterior of the effects given rows 1:t: mode, the observed
            # hybrid items' conditional modes, and the joint covariance.
            items = Int[]
            dhat = Float64[]
            if t == 0
                u = zeros(Float64, d)
                Sigma = Matrix{Float64}(LinearAlgebra.I, d, d)
            else
                _laplace_truncate!(so, full, t)
                # Started from the previous row's mode: one more row moves the
                # posterior a little, and Newton from there converges in a step
                # or two. The same values as from the origin where the
                # posterior has one mode (5e-9 on the pilot-shaped probe), and
                # the probes still choose among maxima where it has more.
                _laplace_solve_unit_mode!(laplace, U, theta, Ls; retain=t > 1)
                u = copy(laplace.modes[U])
                term = _laplace_unit_term(laplace, U, theta, Ls, u, aws)
                out.llrow[offset + t] = term - prevterm
                prevterm = term
                if hybrid
                    modes = _hybrid_member_modes(laplace, U, theta, Ls, u, aws)
                    modes === nothing && continue
                    items = modes.items[1]
                    dhat = modes.dhat[1]
                    layout = _hybrid_layout(d, modes.q)
                    Mj = d > 0 ?
                        _hybrid_joint_curvature(laplace, U, theta, Ls, vcat(u, dhat),
                            layout).diag[1] :
                        modes.A[1]
                else
                    Mj = d > 0 ? _laplace_block_dense(_laplace_unit_curvature(laplace, U,
                        theta, Ls, u), blocks, d) : zeros(Float64, 0, 0)
                end
                Sigma = _laplace_sequential_inverse(Mj)
            end
            # One plug-in pass over rows 1:tt, seeded in the effects, gives row
            # t's updated values and row t+1's predicted ones, both given rows
            # 1:t. Its coordinates: `u`, then every hybrid item, each at its
            # posterior mean given rows 1:t -- zero, with the prior's unit
            # variance, for one not observed by then, which a prediction of
            # that item still carries.
            tt = min(t + 1, T)
            _laplace_truncate!(so, full, tt)
            pass_items = hybrid ? collect(1:_hybrid_items(hyb)) : Int[]
            D = d + length(pass_items)
            z0 = zeros(Float64, D)
            SigmaP = zeros(Float64, D, D)
            z0[1:d] .= u
            known = zeros(Int, length(pass_items))
            for (c, k) in enumerate(pass_items)
                j = findfirst(==(k), items)
                if j === nothing
                    SigmaP[d + c, d + c] = 1.0
                else
                    known[c] = j
                    z0[d + c] = dhat[j]
                end
            end
            # The covariance's blocks over `u` and the known items, from `Sigma`
            # (ordered `u`, then `items`).
            for a in 1:D, b in 1:D
                ia = a <= d ? a : (known[a - d] == 0 ? 0 : d + known[a - d])
                ib = b <= d ? b : (known[b - d] == 0 ? 0 : d + known[b - d])
                (ia == 0 || ib == 0) && continue
                SigmaP[a, b] = Sigma[ia, ib]
            end
            _laplace_sequential_effects!(out, offset, t, T, Ls[1], B,
                pass_items, z0, SigmaP, d, D)
            zd = [_laplace_sequential_seed(S, z0[a], a) for a in 1:D]
            ud = S[zd[a] for a in 1:d]
            x = _laplace_member_values(theta, spec, Ls, ud, units.offsets[U][1])
            for (c, k) in enumerate(pass_items), p in axes(B, 1)
                b = B[p, k]
                iszero(b) || (x[lv.re_index[p]] += b * zd[d + c])
            end
            trace = CTSEMKalmanTrace(S, n, m, tt, 1, nall)
            trace.offset = 0
            trace.current_subject = 1
            ll = _extended_kalman_filter_continuous!(ws, x, so.data, collect(so.timesteps),
                sp, so.tdpreds, _ctsem_tipred_vector(so.tipreds, x), i, so.max_timestep,
                trace, nothing, nothing, so.relin)
            _laplace_finite(ll) || continue
            t >= 1 && _laplace_sequential_record!(out, 2, offset + t, trace, t, SigmaP, D)
            t + 1 <= T && _laplace_sequential_record!(out, 1, offset + t + 1, trace,
                t + 1, SigmaP, D)
        end
    finally
        _laplace_unit_restore!(laplace, U, so, saved)
    end
    return nothing
end

"""
    _laplace_sequential_filtered!(out, laplace, U, theta, Ls, offset, S)

The rows of a unit with no coordinate, from one filter pass. Every effect it
has is integrated by its member's filter (`CTSEMHybrid`), which needs a linear
Gaussian model, or it has none. At deviations `z` (one standard normal
coordinate per integrated item, mapped by `_hybrid_dev_map`) the plug-in
filter's prediction of row `r` is then `value_r + V_r z`, with a covariance
`S_r` that does not depend on `z`, and its log density is quadratic in `z`. So
one pass seeded at `z = 0` gives each row's `V_r`, `S_r`, value and score, and
the posterior of `z` given rows `1:t` is Gaussian: precision
`M_t = I + sum V_r' S_r^-1 V_r` over the observed manifests, mean
`M_t^-1 g_t`, `g_t` the summed scores. Rows are recorded at that mean, with
`J Sigma J'` added as the general path adds it; the record's marginal is
`sum value + g' M^-1 g / 2 - logdet M / 2`, and each row's log predictive
density is its difference with and without the row. The same numbers as the
general path, which filters every prefix again.
"""
function _laplace_sequential_filtered!(out, laplace::CTSEMLaplaceObjective, U::Integer,
    theta::Vector{Float64}, Ls::Vector{Matrix{Float64}}, offset::Integer,
    ::Type{S}) where {S}
    spec = laplace.spec
    i = laplace.units.members[U][1]
    so = laplace.objective.subject_objectives[i]
    T = size(so.data, 2)
    sp = laplace.objective.params
    hyb = sp.hybrid
    lv = spec.levels[1]
    D = _hybrid_on(hyb) ? _hybrid_items(hyb) : 0
    B = D > 0 ? _hybrid_dev_map(hyb, lv, theta) : zeros(Float64, nrandomeffects(lv), 0)
    items = collect(1:D)
    ws = _init_continuous_ekf_workspace(S, sp)
    ws.hyb.mode = _HYBRID_OFF
    n = ws.state_dim
    m = ws.manifest_dim
    x = _laplace_member_values(theta, spec, Ls, S[], laplace.units.offsets[U][1])
    for c in 1:D
        zc = _laplace_sequential_seed(S, 0.0, c)
        for p in axes(B, 1)
            b = B[p, c]
            iszero(b) || (x[lv.re_index[p]] += b * zc)
        end
    end
    trace = CTSEMKalmanTrace(S, n, m, T, 1, length(sp.mutables))
    trace.offset = 0
    trace.current_subject = 1
    ll = _extended_kalman_filter_continuous!(ws, x, so.data, collect(so.timesteps),
        sp, so.tdpreds, _ctsem_tipred_vector(so.tipreds, x), i, so.max_timestep, trace,
        nothing, nothing, so.relin)
    _laplace_finite(ll) || return nothing
    # The posterior given rows 1:t, in information form: H the rows' summed
    # precisions in z, g their summed scores at zero, c their summed values.
    H = zeros(Float64, D, D)
    g = zeros(Float64, D)
    c = 0.0
    zhat = zeros(Float64, D)
    Sigma = Matrix{Float64}(LinearAlgebra.I, D, D)
    prevterm = 0.0
    V = zeros(Float64, m, D)
    SV = zeros(Float64, m, D)
    observed = zeros(Int, m)
    for t in 0:T
        if t >= 1
            lr = trace.llrow[t]
            c += ForwardDiff.value(lr)
            p = 0
            for j in 1:m
                _ctsem_observed(so.data[j, t]) && (p += 1; observed[p] = j)
            end
            if D > 0
                @inbounds for a in 1:D
                    g[a] += ForwardDiff.partials(lr)[a]
                end
                if p > 0
                    Sr = zeros(Float64, p, p)
                    @inbounds for b in 1:p, a in 1:p
                        Sr[a, b] = ForwardDiff.value(trace.ycov[1, t, observed[a], observed[b]])
                    end
                    @inbounds for a in 1:p, e in 1:D
                        V[a, e] = ForwardDiff.partials(trace.y[1, t, observed[a]])[e]
                    end
                    F = _ctsem_cholesky(Sr, p)
                    issuccess(F) || return nothing
                    @inbounds for e in 1:D
                        ldiv!(view(SV, 1:p, e), F, view(V, 1:p, e))
                    end
                    @inbounds for f in 1:D, e in 1:D, a in 1:p
                        H[e, f] += V[a, e] * SV[a, f]
                    end
                end
                M = Matrix{Float64}(LinearAlgebra.I, D, D) .+ H
                FM = _ctsem_cholesky(M, D)
                issuccess(FM) || return nothing
                ldiv!(zhat, FM, g)
                Sigma = inv(FM)
                term = c + dot(g, zhat) / 2 - logdet(FM) / 2
            else
                term = c
            end
            out.llrow[offset + t] = term - prevterm
            prevterm = term
            _laplace_sequential_record!(out, 2, offset + t, trace, t, Sigma, D, zhat)
        end
        _laplace_sequential_effects!(out, offset, t, T, Ls[1], B, items, zhat, Sigma, 0, D)
        t + 1 <= T && _laplace_sequential_record!(out, 1, offset + t + 1, trace,
            t + 1, Sigma, D, zhat)
    end
    return nothing
end

"""
    _laplace_sequential_effects!(out, offset, t, T, L, B, items, z0, Sigma, d, D)

The effects' posterior given rows `1:t` on the raw scale, recorded as row `t`'s
updated value and row `t + 1`'s predicted one: each varying parameter's
deviation `A z0` and covariance `A Sigma A'`, `A` the level's loading for the
effects and `B`'s column (`_hybrid_dev_map`) for each integrated item.
"""
function _laplace_sequential_effects!(out, offset::Integer, t::Integer, T::Integer,
    L::AbstractMatrix, B::AbstractMatrix, items::Vector{Int},
    z0::Vector{Float64}, Sigma::Matrix{Float64}, d::Integer, D::Integer)
    k = size(L, 1)
    k == 0 && return out
    A = zeros(Float64, k, D)
    @inbounds for a in 1:min(d, size(L, 2)), j in 1:k
        A[j, a] = L[j, a]
    end
    for (c, item) in enumerate(items), p in axes(B, 1)
        A[p, d + c] += B[p, item]
    end
    # The engine's own products: this runs per row on every thread, where a BLAS
    # matrix product would queue on OpenBLAS's lock.
    centre = zeros(Float64, k)
    @inbounds for a in 1:D, j in 1:k
        centre[j] += A[j, a] * z0[a]
    end
    AS = zeros(Float64, k, D)
    _ctsem_mul!(AS, A, Sigma)
    spread = zeros(Float64, k, k)
    _ctsem_mulNT!(spread, AS, A)
    if t >= 1
        out.effect[2, offset + t, 1:k] .= centre
        out.effectcov[2, offset + t, 1:k, 1:k] .= spread
    end
    if t + 1 <= T
        out.effect[1, offset + t + 1, 1:k] .= centre
        out.effectcov[1, offset + t + 1, 1:k, 1:k] .= spread
    end
    return out
end

"""One row's means and covariances from a seeded pass: the filter's values at
the effects' posterior means, its covariance plus `J Sigma J'`. A pass seeded
at zero rather than at the means moves its values to them, `+ J centre`."""
function _laplace_sequential_record!(out, k::Integer, row::Integer, trace, r::Integer,
    Sigma::Matrix{Float64}, D::Integer, centre=nothing)
    for (mean, cov, means, covs) in ((trace.eta, trace.etacov, out.eta, out.etacov),
            (trace.y, trace.ycov, out.y, out.ycov))
        p = size(mean, 3)
        J = zeros(Float64, p, D)
        @inbounds for a in 1:p
            e = mean[k, r, a]
            means[k, row, a] = ForwardDiff.value(e)
            for b in 1:D
                J[a, b] = ForwardDiff.partials(e)[b]
            end
            if centre !== nothing
                for b in 1:D
                    means[k, row, a] += J[a, b] * centre[b]
                end
            end
        end
        JS = zeros(Float64, p, D)
        @inbounds for b in 1:D, c in 1:D
            s = Sigma[c, b]
            iszero(s) && continue
            for a in 1:p
                JS[a, b] += J[a, c] * s
            end
        end
        @inbounds for bb in 1:p, a in 1:p
            acc = ForwardDiff.value(cov[k, r, a, bb])
            for c in 1:D
                acc += JS[a, c] * J[bb, c]
            end
            covs[k, row, a, bb] = acc
        end
    end
    return nothing
end

export ctsem_laplace_sequential
