# The marginal route's Hessian by the hand-coded second order: each subject's
# filter sweep carries W = 12 directions as Float64 tangents
# (`tangent_sweep.jl`), where `ctsem_hessian`'s other route runs the adjoint
# gradient at a 12-wide dual. The parameter layer -- subject values, the
# regular transforms, the covariance constructions, their pullbacks and the
# prior -- is the engine's own code, run once a chunk of columns at the
# Hessian's own dual, each column a lane: forward for the parameter matrices'
# tangents, and back with their cotangents carrying the tangents, so the
# result's partials are the columns. Without TI predictors that layer is every
# subject's, run once a chunk; with them each subject's. The filter, which is
# where the dual's compile and memory went, never runs at it. What the
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
`nothing`: what its sweep covers is a Gaussian model, in continuous or
discrete time, with every state diffusing, TI predictors observed, TD
predictors entering as impulses, initial moments stated or stationary (in
continuous time, without integrated deviations or cells evaluated by row), no
population block, deviations integrated whole or profiled (not at fixed
values), and no forward-mode subjects. Cells evaluated by row are covered,
with the iterated update a nonlinear measurement takes and the moment-matched
prediction a bending drift takes, and categorical indicators, but neither
beside integrated deviations.
"""
function _tan_refusal(objective::CTSEMObjective, hybrid_mode::Integer=_HYBRID_PROFILE)
    subjects = objective.subject_objectives
    isempty(subjects) && return "no subjects"
    sp = objective.params
    if _tan_rows_model(sp)
        ws = _init_continuous_ekf_workspace(Float64, sp)
        _hybrid_on(sp.hybrid) && hybrid_mode != _HYBRID_OFF &&
            return "integrated deviations beside cells evaluated by row"
    end
    # The engine itself refuses stationary moments in discrete time.
    sp.stationary && !sp.continuous_time && return "stationary moments in discrete time"
    ws0 = _init_continuous_ekf_workspace(Float64, sp)
    blocks = ws0.affine_buffer.dim < ws0.state_dim ||
        ws0.diffusion_state_indices != 1:ws0.state_dim
    _hybrid_on(sp.hybrid) && blocks && hybrid_mode != _HYBRID_OFF &&
        return "integrated deviations beside static states"
    sp.stationary && _hybrid_on(sp.hybrid) && hybrid_mode != _HYBRID_OFF &&
        return "stationary moments beside integrated deviations"
    _hybrid_on(sp.hybrid) && hybrid_mode == _HYBRID_FIXED &&
        return "deviations at fixed values"
    isfinite(_CTSEM_ADJOINT_GROWTH[]) && return "forward-mode subjects"
    return nothing
end

"""
Whether the subjects' sweeps run by row (`tangent_rows.jl`): a cell evaluated
by row, a categorical indicator, or stationary initial conditions the hand
sweep's own do not cover -- single elements, or beside static states -- each of
which runs as the engine's own at the Hessian's dual.
"""
function _tan_rows_model(sp)
    (any(sp.predict_transforms_indices) || any(sp.update_transforms_indices) ||
        any(sp.td_transforms_indices) || any(t -> t != 0, sp.manifesttype) ||
        _stationary_elements_on(sp.stationary_elements)) && return true
    sp.stationary || return false
    ws = _init_continuous_ekf_workspace(Float64, sp)
    return ws.affine_buffer.dim < ws.state_dim || ws.diffusion_state_indices != 1:ws.state_dim
end

"""
The subject values and cells at `xv` into `aws`'s filter workspace; the cells.
A subject's sampled TI predictor values (`TIMissingRecipe`) are `xv`'s.
"""
function _tan_layer_values!(aws, sp, xv::AbstractVector, tipreds)
    ws = aws.ekf_ws
    _materialize_subject_values!(ws.subject_values, xv, sp,
        _ctsem_tipred_vector(tipreds, xv))
    _materialize_all_params!(ws.all_params, ws.subject_values, sp)
    return ws.all_params
end

"""The parameter matrices at `xv` (any element type), through the filter's layer."""
_tan_layer_mats(aws, sp, xv::AbstractVector, tipreds) =
    _tan_mats_of(aws, sp, _tan_layer_values!(aws, sp, xv, tipreds))

"""Matrices of element type `S` shaped as the layer's (`_TAN_LAYER_NAMES`)."""
_tan_mats_alloc(::Type{S}, n::Int, m::Int, ntd::Int) where {S} =
    NamedTuple{_TAN_LAYER_NAMES}((zeros(S, n, n), zeros(S, n, n), zeros(S, n),
        zeros(S, n, n), zeros(S, m, n), zeros(S, m, n), zeros(S, m), zeros(S, m, m),
        zeros(S, n), zeros(S, n, n), zeros(S, n, ntd), zeros(S, n, n)))

"""The parameter matrices of the cells `all_params`, as the filter builds them."""
function _tan_mats_of(aws, sp, all_params::AbstractVector)
    pars = CTSEMPars(all_params, sp.parameter_axis)
    ws = aws.ekf_ws
    ntd = hasproperty(pars, :TDPREDEFFECT) ? size(pars.TDPREDEFFECT, 2) : 0
    return _tan_mats_of!(_tan_mats_alloc(eltype(all_params), ws.state_dim,
        ws.manifest_dim, ntd), aws, sp, all_params)
end

"""`_tan_mats_of` into the matrices `dst` (`_tan_mats_alloc`)."""
function _tan_mats_of!(dst, aws, sp, all_params::AbstractVector)
    ws = aws.ekf_ws
    pars = CTSEMPars(all_params, sp.parameter_axis)
    n, m = ws.state_dim, ws.manifest_dim
    # As the filter builds them: MANIFESTVAR's and DIFFUSION's covariances as
    # the construction leaves them, T0VAR's with its lower triangle copied up.
    cov!(out, mat, d, buf) = (sdcovsqrt2cov!(buf, mat, sp.covmatcode, d);
        copyto!(out, view(buf.out, 1:d, 1:d)))
    T0V = dst.T0V
    cov!(T0V, pars.T0VAR, n, ws.bufferQ)
    _copy_lower_to_upper!(T0V, n)
    # The population block, added (`_place_population_block!`).
    pidx = sp.population_indices
    if !isempty(pidx)
        k = length(pidx)
        raw = reshape(view(all_params, sp.population_range), k, k)
        sdcovsqrt2cov!(ws.population_buffer, raw, sp.population_covmatcode,
            _square_dim(ws.population_buffer, k))
        @inbounds for b in 1:k, a in 1:k
            T0V[pidx[a], pidx[b]] += sp.population_scale[a] * sp.population_scale[b] *
                ws.population_buffer.out[a, b]
        end
    end
    # The TD predictors' effects and the impulse's Jacobian, where the model
    # has them.
    S = eltype(all_params)
    hasproperty(pars, :TDPREDEFFECT) && copyto!(dst.TDE, pars.TDPREDEFFECT)
    if hasproperty(pars, :Jtd)
        copyto!(dst.Jtd, view(pars.Jtd, 1:n, 1:n))
    else
        @inbounds for j in 1:n, i in 1:n
            dst.Jtd[i, j] = i == j ? one(S) : zero(S)
        end
    end
    copyto!(dst.JAx, view(pars.JAx, 1:n, 1:n))
    copyto!(dst.DRIFT, view(pars.DRIFT, 1:n, 1:n))
    copyto!(dst.CINT, view(pars.CINT, 1:n))
    cov!(dst.Qc, pars.DIFFUSION, n, ws.bufferQ)
    copyto!(dst.LAM, view(pars.LAMBDA, 1:m, 1:n))
    copyto!(dst.Jy, view(pars.Jy, 1:m, 1:n))
    copyto!(dst.MM, view(pars.MANIFESTMEANS, 1:m))
    cov!(dst.R, pars.MANIFESTVAR, m, ws.bufferΘ)
    copyto!(dst.T0M, view(pars.T0MEANS, 1:n))
    return dst
end

"""
A member layer's buffers at element type `S` (`_tan_layer_fill!`,
`_tan_layer_reverse`): the matrices, their cotangents, and the cells', subject
values' and raw values' cotangents. Handed in, the reverse also reads the cells
and subject values the fill materialised on the same workspace, which nothing
between the two touches, rather than materialising them twice more.
"""
function _tan_layer_buffers(::Type{S}, aws, sp, nv::Int, ntd::Int) where {S}
    ws = aws.ekf_ws
    n, m = ws.state_dim, ws.manifest_dim
    return (mats=_tan_mats_alloc(S, n, m, ntd), C=_tan_mats_alloc(S, n, m, ntd),
        tb=zeros(S, length(sp.mutables)), svb=zeros(S, length(ws.subject_values)),
        vb=zeros(S, nv), t0=zeros(S, n, n))
end

"""
The matrices' cotangents `C` added into `tb`, the cells' (`all_params`'), the
covariances pulled back through their constructions at the cells `rawD`,
`rawR` and `rawT0` (DIFFUSION's, MANIFESTVAR's and T0VAR's).
"""
function _tan_bar_of!(tb::AbstractVector{S}, aws, sp, C, rawD, rawR, cells,
    t0=nothing) where {S}
    ws = aws.ekf_ws
    n, m = ws.state_dim, ws.manifest_dim
    tc = CTSEMPars(tb, sp.parameter_axis)
    view(tc.JAx, 1:n, 1:n) .+= C.JAx
    view(tc.DRIFT, 1:n, 1:n) .+= C.DRIFT
    view(tc.CINT, 1:n) .+= C.CINT
    view(tc.LAMBDA, 1:m, 1:n) .+= C.LAM
    view(tc.Jy, 1:m, 1:n) .+= C.Jy
    view(tc.MANIFESTMEANS, 1:m) .+= C.MM
    view(tc.T0MEANS, 1:n) .+= C.T0M
    hasproperty(tc, :TDPREDEFFECT) && (tc.TDPREDEFFECT .+= C.TDE)
    hasproperty(tc, :Jtd) && (view(tc.Jtd, 1:n, 1:n) .+= C.Jtd)
    pullback!(target, mat, cbar, d) = _sdcovsqrt2cov_pullback!(target, mat, cbar, d;
        scratch=aws.covsqrt_scratch, covmatcode=sp.covmatcode)
    pullback!(tc.DIFFUSION, rawD, C.Qc, n)
    pullback!(tc.MANIFESTVAR, rawR, C.R, m)
    T0Vbar = C.T0V
    pidx = sp.population_indices
    if !isempty(pidx) && any(!iszero, T0Vbar)
        # The population block's cotangent, scaled as it was placed, through its
        # own construction (`_pullback_population!`); then its rows and columns
        # masked out of T0VAR's, as the adjoint does.
        k = length(pidx)
        sc = sp.population_scale
        popsym = [(sc[a] * sc[b] * T0Vbar[pidx[a], pidx[b]] +
            sc[b] * sc[a] * T0Vbar[pidx[b], pidx[a]]) / 2 for a in 1:k, b in 1:k]
        popbar = zeros(S, k, k)
        _sdcovsqrt2cov_pullback!(popbar, reshape(view(cells, sp.population_range), k, k),
            popsym, k; scratch=aws.covsqrt_scratch, covmatcode=sp.population_covmatcode)
        reshape(view(tb, sp.population_range), k, k) .+= popbar
        T0Vbar = copy(T0Vbar)
        T0Vbar[pidx, :] .= zero(S)
        T0Vbar[:, pidx] .= zero(S)
    end
    t0 = t0 === nothing ? zeros(S, n, n) : fill!(t0, zero(S))
    pullback!(t0, view(CTSEMPars(cells, sp.parameter_axis).T0VAR, 1:n, 1:n), T0Vbar, n)
    view(tc.T0VAR, 1:n, 1:n) .+= t0
    return tb
end

"""
The raw gradient from the cells' cotangent `tb`: the regular transforms' and
TI effects' pullbacks (`_ctsem_regular_pullback!`, `_ctsem_ti_pullback!`) at
`xv`'s element type; with `buf` (`_tan_layer_buffers`), on the subject values
the fill left and into its buffers.
"""
function _tan_values_reverse(aws, sp, xv::AbstractVector{S}, tb, tipreds,
    buf=nothing) where {S}
    ws = aws.ekf_ws
    buf === nothing && _materialize_subject_values!(ws.subject_values, xv, sp,
        _ctsem_tipred_vector(tipreds, xv))
    # Sized here: a workspace's subject values have their length once first
    # materialised, which a buffer built beside it may precede.
    svb = buf === nothing ? zeros(S, length(ws.subject_values)) :
        fill!(resize!(buf.svb, length(ws.subject_values)), zero(S))
    _ctsem_regular_pullback!(svb, tb, ws.subject_values, sp, aws.regular_supports,
        aws.regular_dual_scratch)
    vb = buf === nothing ? zeros(S, length(xv)) : fill!(buf.vb, zero(S))
    _ctsem_ti_pullback!(vb, svb, sp, tipreds, xv)
    return vb
end

"""
The raw gradient of `<C, matrices(xv)>`, the matrices' cotangents `C` pushed
back through the layer: the adjoint's own pullbacks (`_sdcovsqrt2cov_pullback!`,
`_ctsem_regular_pullback!`, `_ctsem_ti_pullback!`), at `xv`'s element type;
with `buf`, at the cells the fill left (`_tan_layer_buffers`).
"""
function _tan_layer_reverse(aws, sp, xv::AbstractVector{S}, C, tipreds,
    buf=nothing) where {S}
    cells = buf === nothing ? _tan_layer_values!(aws, sp, xv, tipreds) : aws.ekf_ws.all_params
    pars = CTSEMPars(cells, sp.parameter_axis)
    tb = buf === nothing ? zeros(S, length(sp.mutables)) : fill!(buf.tb, zero(S))
    _tan_bar_of!(tb, aws, sp, C, pars.DIFFUSION, pars.MANIFESTVAR, cells,
        buf === nothing ? nothing : buf.t0)
    return _tan_values_reverse(aws, sp, xv, tb, tipreds, buf)
end

const _TAN_LAYER_NAMES = (:JAx, :DRIFT, :CINT, :Qc, :LAM, :Jy, :MM, :R, :T0M, :T0V, :TDE,
    :Jtd)
_tan_tag(::Type{<:ForwardDiff.Dual{T}}) where {T} = T
@inline _tan_layer_field(X, k::Symbol) = getproperty(X, k === :LAM ? :LAMBDA : k)

"""
`f(X's field, N's field)` for each of the layer's matrices (`_TAN_LAYER_NAMES`),
`X` a model or accumulator, `N` the matrices by name: unrolled, so each name is
a constant. A loop over the names reads each field at a run-time symbol, a
dynamic call each.
"""
@inline function _tan_each_layer(f::F, X, N) where {F}
    ntuple(Val(length(_TAN_LAYER_NAMES))) do i
        k = _TAN_LAYER_NAMES[i]
        f(_tan_layer_field(X, k), getproperty(N, k))
        nothing
    end
    return nothing
end

"""The accumulator's cotangents as duals of type `DW`, by name."""
_tan_layer_duals(acc::TanAcc, ::Type{DW}) where {DW} = NamedTuple{_TAN_LAYER_NAMES}(
    ntuple(i -> _tan_duals(_tan_layer_field(acc, _TAN_LAYER_NAMES[i]), DW),
        Val(length(_TAN_LAYER_NAMES))))

"""A worker's sweep, kept in its scratch store between calls."""
function _tan_sweep_scratch(n::Int, m::Int, ntd::Int, nh::Int, ::Val{W}, na::Int=n,
    k::Int=n) where {W}
    store = _ctsem_scratch()
    key = (:tangent_sweep, W, n, m, ntd, nh, na, k)
    sw = get(store, key, nothing)
    if sw === nothing
        sw = TanSweep(n, m, ntd, nh, W, na, k)
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
function _tan_hybrid_own(h::CTSEMHybrid, x::AbstractVector{Float64}, k::Int)
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

"""`X`'s values and tangents from the duals `M`, lane `w` its `w`-th tangent,
those past `nl` zero."""
function _tan_lanes!(X::Union{TanMat{W},TanVec{W}}, M::AbstractArray, nl::Int) where {W}
    @inbounds for idx in eachindex(M)
        X.v[idx] = ForwardDiff.value(M[idx])
        p = ForwardDiff.partials(M[idx])
        for w in 1:W
            X.d[w + W * (idx - 1)] = w <= nl ? p[w] : 0.0
        end
    end
    _tan_setnz_from!(X)
    return X
end

"""`X` as duals of type `DW`: its values, its tangents their partials."""
function _tan_duals(X::Union{TanMat{W},TanVec{W}}, ::Type{DW}) where {W,DW}
    out = similar(X.v, DW)
    @inbounds for idx in eachindex(X.v)
        o = W * (idx - 1)
        out[idx] = DW(X.v[idx], ForwardDiff.Partials(ntuple(w -> X.d[o + w], Val(W))))
    end
    return out
end

"""Whether the duals `a` and `b` are the same, partials included: `==` on a
dual reads its value alone."""
_tan_same(a::AbstractArray, b::AbstractArray) = size(a) == size(b) &&
    all(i -> ForwardDiff.value(a[i]) == ForwardDiff.value(b[i]) &&
        ForwardDiff.partials(a[i]) == ForwardDiff.partials(b[i]), eachindex(a))

"""
The matrices of the layer into `mod`, values and tangents, for the subject
whose TI predictors are `tipreds`: one pass at `xw`, the raw values seeded
with the chunk's `nl` columns as lanes. `false` where the drift is not its own
Jacobian. `buf` (`_tan_layer_buffers`) holds the matrices where given.
"""
function _tan_layer_fill!(mod::TanModel{W}, awsW, sp, xw::Vector, nl::Int,
    tipreds, buf=nothing) where {W}
    mw = buf === nothing ? _tan_layer_mats(awsW, sp, xw, tipreds) :
        _tan_mats_of!(buf.mats, awsW, sp, _tan_layer_values!(awsW, sp, xw, tipreds))
    # Where no cell changes by row, the drift and the measurement must be
    # their own Jacobians; otherwise the offset and the gain read each.
    mod.rows || (_tan_same(mw.JAx, mw.DRIFT) && _tan_same(mw.Jy, mw.LAM)) || return false
    _tan_each_layer((X, M) -> _tan_lanes!(X, M, nl), mod, mw)
    mod.rows && _tan_comb!(mod.DmJ, ((1.0, mod.DRIFT), (-1.0, mod.JAx)))
    _tan_aug_refresh!(mod)
    return true
end

"""
The layer's reverse at `xw`, from the matrices' cotangents in `acc` with their
tangents as partials: lane `w`'s result into `Hc[:, w]`, added, for the
chunk's `nl` columns.
"""
function _tan_layer_columns!(Hc::AbstractMatrix{Float64}, acc::TanAcc{W}, awsW, sp,
    xw::Vector{DW}, nl::Int, tipreds) where {W,DW}
    g = _tan_layer_reverse(awsW, sp, xw, _tan_layer_duals(acc, DW), tipreds)
    @inbounds for w in 1:nl, i in eachindex(g)
        Hc[i, w] += ForwardDiff.partials(g[i], w)
    end
    return Hc
end

"""A model with `mod`'s shape and its hybrid's factor, for a subject's own layer."""
_tan_model_like(mod::TanModel{W}) where {W} = TanModel{W}(mod.n, mod.m,
    mod.ct, mod.stat, mod.rows, TanMat{W}(mod.n, mod.n), TanMat{W}(mod.n, mod.n),
    TanMat{W}(mod.n, mod.n), TanVec{W}(mod.n), TanMat{W}(mod.n, mod.n),
    TanMat{W}(mod.m, mod.n), TanMat{W}(mod.m, mod.n), TanVec{W}(mod.m),
    TanMat{W}(mod.m, mod.m), TanVec{W}(mod.n), TanMat{W}(mod.n, mod.n), TanVec{W}(mod.n),
    Ref(false), mod.gmask, mod.na, mod.dyn, mod.aug, TanMat{W}(mod.na, mod.na),
    TanMat{W}(length(mod.dyn), length(mod.dyn)), TanMat{W}(length(mod.dyn), length(mod.dyn)),
    mod.ntd, TanMat{W}(mod.n, mod.ntd), TanMat{W}(mod.n, mod.n),
    mod.nh, mod.hyb, mod.G, mod.wlog, mod.fixed, mod.dfix)

"""
One block's subjects swept into `acc`. With a shared layer (`L === nothing`)
their cotangents are summed there for the caller's reverse; otherwise each
subject's layer is built into `L.mod`, its cotangents summed in `L.acc`, and
pulled back through its own layer into `Hb`, the block's columns; the
hybrid's factor's cotangent is summed in `acc` either way. A status code.
"""
function _tan_block!(acc::TanAcc{W}, Hb::Matrix{Float64}, sw::TanSweep{W},
    mod::TanModel{W}, L, objective::CTSEMObjective, range, xw::Vector,
    nl::Int) where {W}
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
        _tan_layer_fill!(L.mod, L.awsW, sp, xw, nl, so.tipreds) || return _TAN_REFUSED
        if L.rows !== nothing
            st = _tan_sweep_rows!(L.acc, Hb, sw, L.mod, L.rows, objective, i, wt, xw, nl)
            st == _TAN_OK || return st
            continue
        end
        _tan_sweep_begin!(sw, L.mod)
        sa = _tan_acc_zero!(L.acc)
        st = _tan_sweep!(sa, sw, L.mod, so.timesteps, so.data, so.tdpreds,
            so.max_timestep, wt)
        st == _TAN_OK || return st
        _tan_flush_frechet!(sa, sw.fw, L.mod, sw.sT)
        L.mod.nh > 0 && _tan_comb!(acc.Gbar, ((1.0, sa.Gbar),), true)
        _tan_layer_columns!(Hb, sa, L.awsW, sp, xw, nl, so.tipreds)
    end
    return _TAN_OK
end

function _tan_hessian_columns!(H::Matrix{Float64}, objective::CTSEMObjective,
    x::Vector{Float64}, cols::Vector{Int}, progress, hybrid_mode::Int,
    ::Val{W})::Bool where {W}
    sp = objective.params
    npar = length(x)
    # The layer at the Hessian's own dual, its lanes the chunk's columns; a
    # block's factor at width 1 of the same tag, nested.
    DW = _ctsem_hessian_dual(objective)
    ForwardDiff.npartials(DW) == W || return false
    D1 = ForwardDiff.Dual{_tan_tag(DW),Float64,1}
    subjects = objective.subject_objectives
    # Without TI predictors every subject's layer is the same, and any
    # subject's predictor row gives it; with them each subject has its own,
    # and where cells change by row each subject's cells are its own too.
    rows = _tan_rows_model(sp)
    shared = !rows && isempty(sp.ti_parameter_indices) && isempty(sp.ti_tail)
    nsub = length(subjects)
    nworkers = _ctsem_nchunks(nsub)
    nblocks = _ctsem_nblocks(nsub, nworkers)
    ranges = _ctsem_chunk_ranges(nsub, nblocks)
    aws0 = CTSEMAdjointWorkspace(Float64, sp, npar, 0)
    awsW = CTSEMAdjointWorkspace(DW, sp, npar, 0)
    m0 = _tan_layer_mats(aws0, sp, x, subjects[1].tipreds)
    n, m, ntd = size(m0.JAx, 1), size(m0.LAM, 1), size(m0.TDE, 2)
    # The leading block of genuine dynamics and the diffusing states.
    na, dyn = aws0.affine_dim, collect(aws0.diffusion_state_indices)
    lanes(chunk) = [DW(x[i], ForwardDiff.Partials(ntuple(w -> w <= length(chunk) &&
        cols[chunk[w]] == i ? 1.0 : 0.0, Val(W)))) for i in 1:npar]
    # The deviations the filter integrates, whole (`_HYBRID_FULL`) or profiled.
    h = sp.hybrid
    nh = _hybrid_on(h) && hybrid_mode != _HYBRID_OFF ? _hybrid_items(h) : 0
    gmask = any(t -> t != 0, sp.manifesttype) ?
        Bool[k > length(sp.manifesttype) || sp.manifesttype[k] == 0 for k in 1:m] : Bool[]
    mod = TanModel(n, m, W; ct=sp.continuous_time, stat=sp.stationary, rows=rows,
        gmask=gmask, na=na, dyn=dyn, ntd=ntd, hyb=h, nh=nh,
        wlog=hybrid_mode == _HYBRID_FULL ? 1.0 : 0.0)
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
    layers = shared ? nothing : map(1:nworkers) do _
        awsWk = CTSEMAdjointWorkspace(DW, sp, npar, ntd)
        (awsW=awsWk, mod=_tan_model_like(mod), acc=TanAcc(n, m, ntd, nh, W),
            rows=rows ? TanRows(awsWk, sp, n, m, ntd, Val(W), DW) : nothing)
    end
    Hbs = [zeros(shared ? 0 : npar, shared ? 0 : W) for _ in 1:nblocks]
    status = zeros(Int, nblocks)
    for start in 1:W:length(cols)
        chunk = start:min(start + W - 1, length(cols))
        nl = length(chunk)
        xw = lanes(chunk)
        # The matrices' tangents along the chunk's columns.
        shared && (_tan_layer_fill!(mod, awsW, sp, xw, nl, subjects[1].tipreds) ||
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
            sw = _tan_sweep_scratch(n, m, ntd, nh, Val(W), na, length(dyn))
            status[b] = _tan_block!(accs[b], Hbs[b], sw, mod,
                layers === nothing ? nothing : layers[wk], objective, ranges[b], xw, nl)
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
        # The columns: the layer's reverse, the cotangents carrying their
        # tangents, or the subjects' own, summed; the prior; and the hybrid's
        # factor.
        Hc = view(H, :, chunk)
        fill!(Hc, 0.0)
        if shared
            _tan_layer_columns!(Hc, tot, awsW, sp, xw, nl, subjects[1].tipreds)
        else
            for b in 1:nblocks
                Hc .+= view(Hbs[b], :, 1:nl)
            end
        end
        gp = zeros(DW, npar)
        _ctsem_log_prior_gradient!(gp, objective, xw)
        @inbounds for w in 1:nl, i in 1:npar
            Hc[i, w] += ForwardDiff.partials(gp[i], w)
        end
        for (w, p) in enumerate(chunk)
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
