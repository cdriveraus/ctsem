# The marginal route's Hessian by the hand-coded second order: each subject's
# filter sweep carries W = 12 directions as Float64 tangents
# (`tangent_sweep.jl`), where `ctsem_hessian`'s other route runs the adjoint
# gradient at a 12-wide dual. The parameter layer -- subject values, the
# regular transforms, the covariance constructions, their pullbacks and the
# prior -- is the engine's own code, run once per column at a width-1 dual:
# forward along the column for the parameter matrices' tangents, and back with
# their cotangents carrying the column's tangent, so the result's partial is
# the column. Without TI predictors that layer is every subject's, run once a
# column; with them each subject's, run once a column a subject. What the
# sweep does not cover goes to the dual (`_tan_refusal`), as does a point where
# an interval needs the series. The dual route is the oracle
# (`ctsem_set_tangent_hessian!(false)`, test_tangent_hessian.jl). Deviations the
# filter integrates
# (`tangent_hybrid.jl`) are part of the sweep; their factor `G` is the layer's,
# own sds in closed form, a block's factor through `_hybrid_block_factor!` at
# the width-1 dual.

const _CTSEM_TANGENT_HESSIAN = Ref(true)

"""
    ctsem_set_tangent_hessian!(on)

Whether `ctsem_hessian` and `ctsem_hessian_columns` take the hand-coded second
order where it covers the model (`true`, the default) or forward mode over the
adjoint gradient at its dual (`false`): the oracle the hand-coded route is
checked against.
"""
function ctsem_set_tangent_hessian!(on::Bool)
    _CTSEM_TANGENT_HESSIAN[] = on
    return on
end
export ctsem_set_tangent_hessian!

const _TAN_WIDTH = 12

"""
Why the hand-coded Hessian does not serve `objective` in `hybrid_mode`, or
`nothing`: what its sweep covers is a linear Gaussian model, in continuous or
discrete time, with every state diffusing, TI predictors observed, TD
predictors entering as impulses, initial moments stated or stationary (in
continuous time, without integrated deviations), no population block,
deviations integrated whole or profiled (not at fixed values), and no
forward-mode subjects.
"""
function _tan_refusal(objective::CTSEMObjective, hybrid_mode::Integer=_HYBRID_PROFILE)
    subjects = objective.subject_objectives
    isempty(subjects) && return "no subjects"
    sp = objective.params
    any(sp.predict_transforms_indices) && return "per-row predict transforms"
    any(sp.update_transforms_indices) && return "per-row update transforms"
    any(sp.td_transforms_indices) && return "per-row TD transforms"
    any(t -> t != 0, sp.manifesttype) && return "categorical indicators"
    _stationary_elements_on(sp.stationary_elements) && return "stationary elements"
    sp.stationary && !sp.continuous_time && return "stationary moments in discrete time"
    sp.stationary && _hybrid_on(sp.hybrid) && hybrid_mode != _HYBRID_OFF &&
        return "stationary moments beside integrated deviations"
    isempty(sp.population_indices) || return "a population block"
    isempty(objective.ti_missing_parameter) || return "sampled TI predictor values"
    _hybrid_on(sp.hybrid) && hybrid_mode == _HYBRID_FIXED &&
        return "deviations at fixed values"
    isfinite(_CTSEM_ADJOINT_GROWTH[]) && return "forward-mode subjects"
    for so in subjects
        so.tipreds isa AbstractVector || return "sampled TI predictor values"
    end
    return nothing
end

"""The parameter matrices at `xv` (any element type), through the filter's layer."""
function _tan_layer_mats(aws, sp, xv::AbstractVector, tipreds)
    ws = aws.ekf_ws
    _materialize_subject_values!(ws.subject_values, xv, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = CTSEMPars(ws.all_params, sp.parameter_axis)
    n, m = ws.state_dim, ws.manifest_dim
    # As the filter builds them: MANIFESTVAR's and DIFFUSION's covariances as
    # the construction leaves them, T0VAR's with its lower triangle copied up.
    cov(mat, d, buf) = (sdcovsqrt2cov!(buf, mat, sp.covmatcode, d);
        Matrix(view(buf.out, 1:d, 1:d)))
    T0V = cov(pars.T0VAR, n, ws.bufferQ)
    _copy_lower_to_upper!(T0V, n)
    # The TD predictors' effects and the impulse's Jacobian, where the model
    # has them.
    S = eltype(ws.all_params)
    TDE = hasproperty(pars, :TDPREDEFFECT) ? Matrix(pars.TDPREDEFFECT) : zeros(S, n, 0)
    Jtd = hasproperty(pars, :Jtd) ? Matrix(view(pars.Jtd, 1:n, 1:n)) : Matrix{S}(I, n, n)
    return (JAx=Matrix(view(pars.JAx, 1:n, 1:n)), DRIFT=Matrix(view(pars.DRIFT, 1:n, 1:n)),
        CINT=collect(view(pars.CINT, 1:n)), Qc=cov(pars.DIFFUSION, n, ws.bufferQ),
        LAM=Matrix(view(pars.LAMBDA, 1:m, 1:n)), Jy=Matrix(view(pars.Jy, 1:m, 1:n)),
        MM=collect(view(pars.MANIFESTMEANS, 1:m)), R=cov(pars.MANIFESTVAR, m, ws.bufferΘ),
        T0M=collect(view(pars.T0MEANS, 1:n)), T0V=T0V, TDE=TDE, Jtd=Jtd)
end

"""
The raw gradient of `<C, matrices(xv)>`, the matrices' cotangents `C` pushed
back through the layer: the adjoint's own pullbacks (`_sdcovsqrt2cov_pullback!`,
`_ctsem_regular_pullback!`, `_ctsem_ti_pullback!`), at `xv`'s element type.
"""
function _tan_layer_reverse(aws, sp, xv::AbstractVector{S}, C, tipreds) where {S}
    ws = aws.ekf_ws
    n, m = ws.state_dim, ws.manifest_dim
    _materialize_subject_values!(ws.subject_values, xv, sp, tipreds)
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    pars = CTSEMPars(ws.all_params, sp.parameter_axis)
    tb = zeros(S, length(sp.mutables))
    tc = CTSEMPars(tb, sp.parameter_axis)
    view(tc.JAx, 1:n, 1:n) .+= C.JAx
    view(tc.CINT, 1:n) .+= C.CINT
    view(tc.LAMBDA, 1:m, 1:n) .+= C.LAM
    view(tc.MANIFESTMEANS, 1:m) .+= C.MM
    view(tc.T0MEANS, 1:n) .+= C.T0M
    hasproperty(tc, :TDPREDEFFECT) && (tc.TDPREDEFFECT .+= C.TDE)
    hasproperty(tc, :Jtd) && (view(tc.Jtd, 1:n, 1:n) .+= C.Jtd)
    pullback!(target, mat, cbar, d) = _sdcovsqrt2cov_pullback!(target, mat, cbar, d;
        scratch=aws.covsqrt_scratch, covmatcode=sp.covmatcode)
    pullback!(tc.DIFFUSION, pars.DIFFUSION, C.Qc, n)
    pullback!(tc.MANIFESTVAR, pars.MANIFESTVAR, C.R, m)
    t0 = zeros(S, n, n)
    pullback!(t0, pars.T0VAR, C.T0V, n)
    view(tc.T0VAR, 1:n, 1:n) .+= t0
    svb = zeros(S, length(ws.subject_values))
    _ctsem_regular_pullback!(svb, tb, ws.subject_values, sp, aws.regular_supports,
        aws.regular_dual_scratch)
    vb = zeros(S, length(xv))
    _ctsem_ti_pullback!(vb, svb, sp, tipreds, xv)
    return vb
end

const _TAN_LAYER_NAMES = (:JAx, :CINT, :Qc, :LAM, :MM, :R, :T0M, :T0V, :TDE, :Jtd)
_tan_tag(::Type{<:ForwardDiff.Dual{T}}) where {T} = T
@inline _tan_layer_field(X, k::Symbol) = getproperty(X, k === :LAM ? :LAMBDA : k)

"""A worker's sweep, kept in its scratch store between calls."""
function _tan_sweep_scratch(n::Int, m::Int, ntd::Int, nh::Int, ::Val{W}) where {W}
    store = _ctsem_scratch()
    key = (:tangent_sweep, W, n, m, ntd, nh)
    sw = get(store, key, nothing)
    if sw === nothing
        sw = TanSweep(n, m, ntd, nh, W)
        store[key] = sw
    end
    return sw::TanSweep{W}
end

"""
The Hessian's columns `cols` into `H` (`length(x)` by `length(cols)`) by the
hand-coded second order. `false` where it does not serve this point -- the
drift is not the Jacobian, or an interval needs the series -- and the caller
takes the dual route; `_tan_refusal(objective)` must be `nothing`.
"""
_tan_hessian_columns!(H::Matrix{Float64}, objective::CTSEMObjective,
    x::Vector{Float64}, cols::Vector{Int}, progress,
    hybrid_mode::Integer=_HYBRID_PROFILE)::Bool =
    _tan_hessian_columns!(H, objective, x, cols, progress, Int(hybrid_mode),
        Val(_TAN_WIDTH))

"""
The hybrid's own sds: each item's `g = slope (log1p_exp(2 r - 1) sd + 1e-10)`
at its raw `r` (`_hybrid_factor!`), and the first two derivatives in `r`.
"""
function _tan_hybrid_own(h::CTSEMHybrid, x::Vector{Float64}, k::Int)
    r = x[h.psi_index[k]]
    z = 2 * r - 1
    sig = 1 / (1 + exp(-z))
    sl, sd = h.slope[k], h.sd_scale[k]
    return (sl * (log1p_exp(z) * sd + 1e-10), sl * sd * 2 * sig, sl * sd * 4 * sig * (1 - sig))
end

"""
A block's factor, the rows of the trailing block of the level's covariance
factor the hybrid reads (`_hybrid_factor!`): `G[c, f] = slope_c Lb[nL+f, nL+c]`
for `f <= c`, from `Lb` at any element type.
"""
function _tan_hybrid_block_G(h::CTSEMHybrid, Lb::AbstractMatrix, c::Int, f::Int)
    nh = _hybrid_items(h)
    nL = length(h.perm) - nh
    return f <= c ? h.slope[c] * Lb[nL + f, nL + c] : zero(eltype(Lb))
end

"""
The matrices of the layer at `x` into `mod`'s values, and their tangents along
each seeded copy of `x` in `xs` into its tangents, for the subject whose TI
predictors are `tipreds`. `false` where the drift is not its own Jacobian.
"""
function _tan_layer_fill!(mod::TanModel{W}, aws0, aws1, sp, x::Vector{Float64},
    xs::Vector, tipreds) where {W}
    m0 = _tan_layer_mats(aws0, sp, x, tipreds)
    (m0.JAx == m0.DRIFT && m0.Jy == m0.LAM) || return false
    for k in _TAN_LAYER_NAMES
        X = _tan_layer_field(mod, k)
        X.v .= getproperty(m0, k)
        fill!(X.d, 0.0)
    end
    for (w, xw) in enumerate(xs)
        mw = _tan_layer_mats(aws1, sp, xw, tipreds)
        (mw.JAx == mw.DRIFT && mw.Jy == mw.LAM) || return false
        for k in _TAN_LAYER_NAMES
            X = _tan_layer_field(mod, k)
            P = ForwardDiff.partials.(getproperty(mw, k), 1)
            X.d isa Matrix ? (X.d[w, :] .= P) : (X.d[w, :, :] .= P)
        end
    end
    foreach(k -> _tan_setnz_from!(_tan_layer_field(mod, k)), _TAN_LAYER_NAMES)
    return true
end

"""
The layer's reverse along each seeded copy of `x` in `xs`, from the matrices'
cotangents in `acc` with their tangents: the `w`-th column into `Hc[:, w]`,
added.
"""
function _tan_layer_columns!(Hc::AbstractMatrix{Float64}, acc::TanAcc{W}, aws1, sp,
    xs::Vector{Vector{D1}}, tipreds) where {W,D1}
    for (w, xw) in enumerate(xs)
        C = NamedTuple{_TAN_LAYER_NAMES}(map(_TAN_LAYER_NAMES) do k
            X = _tan_layer_field(acc, k)
            dv = X.d isa Matrix ? X.d[w, :] : X.d[w, :, :]
            D1.(X.v, ForwardDiff.Partials.(tuple.(dv)))
        end)
        g = _tan_layer_reverse(aws1, sp, xw, C, tipreds)
        @inbounds for i in eachindex(g)
            Hc[i, w] += ForwardDiff.partials(g[i], 1)
        end
    end
    return Hc
end

"""A model with `mod`'s shape and its hybrid's factor, for a subject's own layer."""
_tan_model_like(mod::TanModel{W}) where {W} = TanModel{W}(mod.n, mod.m,
    mod.ct, mod.stat, TanMat{W}(mod.n, mod.n), TanVec{W}(mod.n), TanMat{W}(mod.n, mod.n),
    TanMat{W}(mod.m, mod.n), TanVec{W}(mod.m), TanMat{W}(mod.m, mod.m), TanVec{W}(mod.n),
    TanMat{W}(mod.n, mod.n), mod.ntd, TanMat{W}(mod.n, mod.ntd), TanMat{W}(mod.n, mod.n),
    mod.nh, mod.hyb, mod.G, mod.wlog)

"""
One block's subjects swept into `acc`. With a shared layer (`L === nothing`)
their cotangents are summed there for the caller's reverse; otherwise each
subject's layer is built into `L.mod`, its cotangents summed in `L.acc`, and
pulled back through its own layer into `Hb`, the block's columns; the
hybrid's factor's cotangent is summed in `acc` either way. A status code.
"""
function _tan_block!(acc::TanAcc{W}, Hb::Matrix{Float64}, sw::TanSweep{W},
    mod::TanModel{W}, L, objective::CTSEMObjective, range, x::Vector{Float64},
    xs::Vector) where {W}
    subjects, weights, sp = objective.subject_objectives, objective.subject_weights,
        objective.params
    _tan_acc_zero!(acc)
    if L === nothing
        _tan_sweep_begin!(sw, mod)
        for i in range
            so = subjects[i]
            wt = isempty(weights) ? 1.0 : weights[i]
            wt == 0 && continue
            st = _tan_sweep!(acc, sw, mod, so.timesteps, so.data, so.tdpreds,
                so.max_timestep, wt)
            st == _TAN_OK || return st
        end
        _tan_flush_frechet!(acc, sw.fw, mod, sw.sT)
        return _TAN_OK
    end
    fill!(Hb, 0.0)
    for i in range
        so = subjects[i]
        wt = isempty(weights) ? 1.0 : weights[i]
        wt == 0 && continue
        _tan_layer_fill!(L.mod, L.aws0, L.aws1, sp, x, xs, so.tipreds) || return _TAN_REFUSED
        _tan_sweep_begin!(sw, L.mod)
        sa = _tan_acc_zero!(L.acc)
        st = _tan_sweep!(sa, sw, L.mod, so.timesteps, so.data, so.tdpreds,
            so.max_timestep, wt)
        st == _TAN_OK || return st
        _tan_flush_frechet!(sa, sw.fw, L.mod, sw.sT)
        L.mod.nh > 0 && _tan_comb!(acc.Gbar, ((1.0, sa.Gbar),), true)
        _tan_layer_columns!(Hb, sa, L.aws1, sp, xs, so.tipreds)
    end
    return _TAN_OK
end

function _tan_hessian_columns!(H::Matrix{Float64}, objective::CTSEMObjective,
    x::Vector{Float64}, cols::Vector{Int}, progress, hybrid_mode::Int,
    ::Val{W})::Bool where {W}
    sp = objective.params
    npar = length(x)
    # Width 1 at the Hessian's own tag: the layer's own seeded duals nest inside
    # it as they nest inside the 12-wide one.
    D1 = ForwardDiff.Dual{_tan_tag(_ctsem_hessian_dual(objective)),Float64,1}
    subjects = objective.subject_objectives
    # Without TI predictors every subject's layer is the same, and any
    # subject's predictor row gives it; with them each subject has its own.
    shared = isempty(sp.ti_parameter_indices) && isempty(sp.ti_tail)
    nsub = length(subjects)
    nworkers = _ctsem_nchunks(nsub)
    nblocks = _ctsem_nblocks(nsub, nworkers)
    ranges = _ctsem_chunk_ranges(nsub, nblocks)
    aws0 = CTSEMAdjointWorkspace(Float64, sp, npar, 0)
    aws1 = CTSEMAdjointWorkspace(D1, sp, npar, 0)
    m0 = _tan_layer_mats(aws0, sp, x, subjects[1].tipreds)
    n, m, ntd = size(m0.JAx, 1), size(m0.LAM, 1), size(m0.TDE, 2)
    # Every state diffusing, none static.
    (aws0.diffusion_state_indices == 1:n && aws0.affine_dim == n) || return false
    seed(c) = [D1(x[i], ForwardDiff.Partials((i == c ? 1.0 : 0.0,))) for i in 1:npar]
    # The deviations the filter integrates, whole (`_HYBRID_FULL`) or profiled.
    h = sp.hybrid
    nh = _hybrid_on(h) && hybrid_mode != _HYBRID_OFF ? _hybrid_items(h) : 0
    mod = TanModel(n, m, W; ct=sp.continuous_time, stat=sp.stationary, ntd=ntd, hyb=h,
        nh=nh, wlog=hybrid_mode == _HYBRID_FULL ? 1.0 : 0.0)
    Lb = h.block && nh > 0 ? zeros(length(h.perm), length(h.perm)) : zeros(0, 0)
    if nh > 0 && h.block
        if !_hybrid_block_factor!(Lb, h, x)
            fill!(H, NaN)
            return true
        end
    end
    # A block's factor's derivatives in its level's parameters, once: the
    # columns read them by arithmetic.
    hd = nh > 0 && h.block ? _tan_hybrid_block_derivs(h, x, D1) : nothing
    accs = [TanAcc(n, m, ntd, nh, W) for _ in 1:nblocks]
    # Each worker's own layer and the cotangents of one subject, and each
    # block's columns, where the layer is the subject's.
    layers = shared ? nothing : [(aws0=CTSEMAdjointWorkspace(Float64, sp, npar, 0),
        aws1=CTSEMAdjointWorkspace(D1, sp, npar, 0), mod=_tan_model_like(mod),
        acc=TanAcc(n, m, ntd, nh, W)) for _ in 1:nworkers]
    Hbs = [zeros(shared ? 0 : npar, shared ? 0 : W) for _ in 1:nblocks]
    status = zeros(Int, nblocks)
    for start in 1:W:length(cols)
        chunk = start:min(start + W - 1, length(cols))
        xs = [seed(cols[p]) for p in chunk]
        # The matrices' tangents along the chunk's columns.
        shared && (_tan_layer_fill!(mod, aws0, aws1, sp, x, xs, subjects[1].tipreds) ||
            return false)
        # The hybrid's factor, and its tangents along the chunk's columns.
        if nh > 0
            G = mod.G
            fill!(G.v, 0.0)
            fill!(G.d, 0.0)
            if !h.block
                for k in 1:nh
                    g, g1, _ = _tan_hybrid_own(h, x, k)
                    G.v[k, k] = g
                    for (w, p) in enumerate(chunk)
                        cols[p] == h.psi_index[k] && (G.d[w, k, k] = g1)
                    end
                end
            else
                for c in 1:nh, f in 1:c
                    G.v[c, f] = _tan_hybrid_block_G(h, Lb, c, f)
                end
                for (w, p) in enumerate(chunk)
                    t = hd.pos[cols[p]]
                    t == 0 && continue
                    for c in 1:nh, f in 1:c
                        G.d[w, c, f] = hd.dG[c, f, t]
                    end
                end
            end
            _tan_setnz_from!(G)
        end
        # The subjects, in the gradient's blocks on its workers; each block's
        # exponential directions flushed once for all its subjects where the
        # layer is shared, once a subject where it is not.
        fill!(status, _TAN_OK)
        _ctsem_pull_blocks(nblocks, nworkers) do b, wk
            sw = _tan_sweep_scratch(n, m, ntd, nh, Val(W))
            status[b] = _tan_block!(accs[b], Hbs[b], sw, mod,
                layers === nothing ? nothing : layers[wk], objective, ranges[b], x, xs)
            return nothing
        end
        any(==(_TAN_REFUSED), status) && return false
        if any(==(_TAN_INVALID), status)
            # The gradient there is an invalid point's, NaN; so is its Hessian.
            fill!(H, NaN)
            return true
        end
        tot = accs[1]
        for b in 2:nblocks
            _tan_acc_add!(tot, accs[b])
        end
        # Each column: the layer's reverse at a width-1 dual seeded along it,
        # the cotangents carrying their tangent along it, or the subjects'
        # own, summed; the prior; and the hybrid's factor.
        Hc = view(H, :, chunk)
        fill!(Hc, 0.0)
        if shared
            _tan_layer_columns!(Hc, tot, aws1, sp, xs, subjects[1].tipreds)
        else
            for b in 1:nblocks
                Hc .+= view(Hbs[b], :, 1:length(chunk))
            end
        end
        for (w, p) in enumerate(chunk)
            gp = zeros(D1, npar)
            _ctsem_log_prior_gradient!(gp, objective, xs[w])
            @inbounds for i in 1:npar
                H[i, p] += ForwardDiff.partials(gp[i], 1)
            end
            # The hybrid's factor's cotangent, straight onto its raw positions
            # (`_hybrid_scatter!`): `Gbar` with its tangent along the column,
            # times G's first derivative, and `Gbar` times its second.
            nh > 0 && _tan_hybrid_layer!(H, p, w, cols[p], tot.Gbar, h, x, hd)
        end
        progress === nothing || progress(last(chunk), length(cols), "parameters")
    end
    return true
end

"""
A block's factor's derivatives in its level's parameters at `x`: `dG[:, :, j]`
along `h.lparams[j]` and `d2G[:, :, j, t]` along it and `h.lparams[t]`, from
`_hybrid_block_factor!` at a dual of its own tag nested in the width-1 one
(`D1`), one pass a pair; and `pos`, each raw parameter's place in `h.lparams`,
0 for none.
"""
function _tan_hybrid_block_derivs(h::CTSEMHybrid, x::Vector{Float64},
    ::Type{D1}) where {D1}
    nh = _hybrid_items(h)
    nl = length(h.lparams)
    k = length(h.perm)
    D2 = ForwardDiff.Dual{_HybridFactorTag,D1,1}
    pos = zeros(Int, length(x))
    for (j, i) in enumerate(h.lparams)
        pos[i] = j
    end
    dG = zeros(nh, nh, nl)
    d2G = zeros(nh, nh, nl, nl)
    vd = Vector{D2}(undef, length(x))
    Lb = zeros(D2, k, k)
    for t in 1:nl, j in t:nl
        @inbounds for i in eachindex(x)
            vd[i] = D2(D1(x[i], ForwardDiff.Partials((i == h.lparams[t] ? 1.0 : 0.0,))),
                ForwardDiff.Partials((i == h.lparams[j] ? one(D1) : zero(D1),)))
        end
        _hybrid_block_factor!(Lb, h, vd) || continue
        for c in 1:nh, f in 1:c
            dj = ForwardDiff.partials(_tan_hybrid_block_G(h, Lb, c, f), 1)
            t == 1 && (dG[c, f, j] = ForwardDiff.value(dj))
            d2 = ForwardDiff.partials(dj, 1)
            d2G[c, f, j, t] = d2
            d2G[c, f, t, j] = d2
        end
    end
    return (pos=pos, dG=dG, d2G=d2G)
end

"""
The hybrid factor's cotangent `Gbar` (summed over the subjects, with its
tangents) pushed onto its raw positions, column `p` of `H` along raw direction
`c` (its `w`-th tangent): `Gbar`'s tangent times G's first derivative, and
`Gbar` times its second along `c`; own sds in closed form, a block's from
`hd` (`_tan_hybrid_block_derivs`).
"""
function _tan_hybrid_layer!(H::Matrix{Float64}, p::Int, w::Int, c::Int,
    Gbar::TanMat, h::CTSEMHybrid, x::Vector{Float64}, hd)
    nh = _hybrid_items(h)
    if !h.block
        @inbounds for k in 1:nh
            idx = h.psi_index[k]
            _, g1, g2 = _tan_hybrid_own(h, x, k)
            H[idx, p] += Gbar.d[w, k, k] * g1 + (c == idx ? Gbar.v[k, k] * g2 : 0.0)
        end
        return H
    end
    t = hd.pos[c]
    @inbounds for j in eachindex(h.lparams)
        acc = 0.0
        for cc in 1:nh, f in 1:cc
            acc += Gbar.d[w, cc, f] * hd.dG[cc, f, j]
            t == 0 || (acc += Gbar.v[cc, f] * hd.d2G[cc, f, j, t])
        end
        H[h.lparams[j], p] += acc
    end
    return H
end
