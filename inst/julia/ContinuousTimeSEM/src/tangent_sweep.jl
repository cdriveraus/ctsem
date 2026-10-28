# A subject's filter sweep with its tangents by hand, on `tangent_primitives.jl`
# and `tangent_kernels.jl`: the forward pass, recorded, then the reverse walk,
# for a linear Gaussian continuous-time model (DRIFT = JAx, LAMBDA = Jy), every
# state diffusing (`dyn = 1:n`, `naff = n`), any manifests missing, and the
# discretisation's closed forms (an interval needing the series refuses). The
# forward is the operations of `my_exp!`, `_compute_discrete_time_form!`,
# `_ekf_predict_step!` and `_ekf_masked_update_step!` in their order, the
# reverse those of `_reverse_update!` and `_reverse_predict!`, each paired with
# its tangent: what the adjoint gradient's subject pass computes at the
# Hessian's dual, up to the parameter layer (`tangent_hessian.jl`), into the
# parameter matrices' cotangents and their tangents.

"""`f.Y = exp(A)` with tangents, as `my_exp!` computes it (Padé 13, squaring
in pairs); `f`'s other buffers are scratch."""
function _tan_exp!(f::TanFrechet{W}, A::TanMat{W}) where {W}
    n = f.n
    As, Y = f.As, f.Y
    a1 = _opnorm1_noalloc(A.v, n)
    if !isfinite(a1)
        fill!(Y.v, NaN)
        fill!(Y.d, NaN)
        Y.nz[] = true
        Y.jet && (fill!(Y.s, NaN); Y.snz[] = true)
        return Y
    end
    theta = 5.371920351148152
    s = a1 <= theta ? 0 : ceil(Int, log2(a1 / theta))
    _tan_copy!(As, A)
    _tan_scale!(As, inv(2.0^s))
    A2, A4, A6, W1, Wm, Z1, V, U, T1 = f.A2, f.A4, f.A6, f.W1, f.Wm, f.Z1, f.V, f.U, f.T1
    _tan_mul!(A2, As, As)
    _tan_mul!(A4, A2, A2)
    _tan_mul!(A6, A4, A2)
    b = _pade13_b(Float64)
    _tan_comb!(W1, ((b[14], A6), (b[12], A4), (b[10], A2)))
    _tan_mul!(Wm, A6, W1)
    _tan_comb!(Wm, ((b[8], A6), (b[6], A4), (b[4], A2)), true)
    _tan_diag!(Wm, b[2])
    _tan_mul!(U, As, Wm)
    _tan_comb!(Z1, ((b[13], A6), (b[11], A4), (b[9], A2)))
    _tan_mul!(V, A6, Z1)
    _tan_comb!(V, ((b[7], A6), (b[5], A4), (b[3], A2)), true)
    _tan_diag!(V, b[1])
    _tan_comb!(T1, ((1.0, V), (-1.0, U)))
    _tan_comb!(Y, ((1.0, V), (1.0, U)))
    copyto!(f.F, T1.v)
    _lu_factor_generic!(f.F, f.piv, n)
    _tan_solve!(Y, f.F, f.piv, T1)
    for _ in 1:div(s, 2)
        _tan_mul!(W1, Y, Y)
        _tan_mul!(Y, W1, W1)
    end
    if isodd(s)
        _tan_mul!(W1, Y, Y)
        _tan_copy!(Y, W1)
    end
    return Y
end

"""What the reverse of a prediction reads (`CTSEMPredictRecord`), with tangents."""
struct TanPredRec{W}
    A::TanMat{W}            # exp(JAx dt)
    x::TanVec{W}            # the state in
    P::TanMat{W}            # the covariance in
    affine::TanVec{W}       # CINT: DRIFT = JAx, so the offset reads no state
    dINT::TanVec{W}
    dt::Base.RefValue{Float64}
    Xh::TanMat{W}           # the hybrid's columns in, and the intercepts' increment
    Xinc::TanMat{W}
    sint::Base.RefValue{Bool}   # the intercept's integral taken as a series
    snoise::Base.RefValue{Bool} # and the noise's
end
TanPredRec(n::Int, nh::Int, W::Int) = TanPredRec{W}(TanMat{W}(n, n), TanVec{W}(n),
    TanMat{W}(n, n), TanVec{W}(n), TanVec{W}(n), Ref(0.0), TanMat{W}(n, nh), TanMat{W}(n, nh),
    Ref(false), Ref(false))

"""
What the reverse of an update reads (`CTSEMUpdateRecord`), with tangents: the
manifests observed (`o`, none where the row updates nothing) and their values.
An iterated update's (`relin` 1 at its first step, 2 after): the measurement
linearised at `xl`, its loading Jy, its intercept `c = MANIFESTMEANS +
(LAMBDA - Jy) xl` on the observed rows, and those rows of LAMBDA (`Lam`).
"""
mutable struct TanUpdRec{W}
    o::Vector{Int}
    P::TanMat{W}            # the predicted covariance, unridged
    x::TanVec{W}            # the predicted state
    y::TanVec{W}            # the observed values (no tangent)
    Xh::TanMat{W}           # the hybrid's columns in
    relin::Int
    c::TanVec{W}
    xl::TanVec{W}
    Lam::TanMat{W}
end
TanUpdRec(n::Int, m::Int, nh::Int, W::Int) = TanUpdRec{W}(Int[], TanMat{W}(n, n),
    TanVec{W}(n), TanVec{W}(m), TanMat{W}(n, nh), 0, TanVec{W}(0), TanVec{W}(n),
    TanMat{W}(0, n))

"""
An update's scratch for a row observing `mo` manifests: the gathered loadings,
means and measurement covariance, the forward's intermediates, and the
reverse's workspace (`TanUpdate`).
"""
struct TanRowScratch{W}
    H::TanMat{W}; L::TanMat{W}; MM::TanVec{W}; R::TanMat{W}
    Pr::TanMat{W}; PHt::TanMat{W}; S::TanMat{W}; Ui::TanMat{W}; Xc::TanMat{W}; Tc::TanMat{W}; F::Matrix{Float64}
    yt::TanVec{W}; e::TanVec{W}; a::TanVec{W}; Kb::TanMat{W}; K::TanMat{W}; Mm::TanMat{W}; nn1::TanMat{W}; nm1::TanMat{W}
    u::TanUpdate{W}
    hr::TanHybRow{W}
end
function TanRowScratch(n::Int, mo::Int, nh::Int, W::Int)
    T, V = TanMat{W}, TanVec{W}
    return TanRowScratch{W}(T(mo, n), T(mo, n), V(mo), T(mo, mo),
        T(n, n), T(n, mo), T(mo, mo), T(mo, mo), T(mo, mo), T(mo, mo), zeros(mo, mo),
        V(mo), V(mo), V(mo), T(n, mo), T(n, mo), T(n, n), T(n, n), T(n, mo),
        TanUpdate(n, mo, W), TanHybRow(n, mo, nh, W))
end

"""
The parameter matrices of a Gaussian model, with tangents, `ntd` TD
predictors' among them; whether it is in continuous time (`ct`), starts from
its stationary moments (`stat`), and has cells evaluated by row (`rows`,
`tangent_rows.jl`), where DRIFT need not be the drift's Jacobian JAx, nor the
measurement's Jacobian Jy its loadings: the local affine model's offset
reads `DmJ = DRIFT - JAx`, the update's gain `Jy`, and where the prediction
is moment-matched (`momon`) the offset gains the cubature's mean correction
`mc` (`_ekf_moment_cells!`). Where the filter
integrates deviations, their spec, their factor `G` and the closing term's
weight on the logdet (1 for the integral, 0 for the profile); or, `fixed`,
the deviations held at `dfix` (`_HYBRID_FIXED`), one per item.
"""
struct TanModel{W}
    n::Int
    m::Int
    ct::Bool
    stat::Bool
    rows::Bool
    JAx::TanMat{W}; DRIFT::TanMat{W}; DmJ::TanMat{W}; CINT::TanVec{W}; Qc::TanMat{W}
    LAMBDA::TanMat{W}; Jy::TanMat{W}; MM::TanVec{W}; R::TanMat{W}
    T0M::TanVec{W}; T0V::TanMat{W}
    mc::TanVec{W}; momon::Base.RefValue{Bool}
    gmask::Vector{Bool}         # the Gaussian manifests (empty: every one)
    # The leading block of genuine dynamics (`na`) and the diffusing states
    # (`dyn`), which the intercept's solve and the noise's Lyapunov solve run
    # over; `aug` where they are not every state, as with the random effects'
    # carriers of the augmented route. JAx's and Qc's blocks.
    na::Int
    dyn::Vector{Int}
    aug::Bool
    JAa::TanMat{W}; JAd::TanMat{W}; Qd::TanMat{W}
    ntd::Int
    TDE::TanMat{W}; Jtd::TanMat{W}
    nh::Int
    hyb::CTSEMHybrid
    G::TanMat{W}
    wlog::Float64
    fixed::Bool
    dfix::TanVec{W}
end
TanModel(n::Int, m::Int, W::Int; ct::Bool=true, stat::Bool=false, rows::Bool=false,
    gmask::Vector{Bool}=Bool[], na::Int=n, dyn::Vector{Int}=collect(1:n), ntd::Int=0,
    hyb::CTSEMHybrid=CTSEMHybrid(), nh::Int=0, wlog::Float64=1.0, fixed::Bool=false) =
    TanModel{W}(n, m, ct, stat, rows, TanMat{W}(n, n), TanMat{W}(n, n), TanMat{W}(n, n),
    TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(m, n), TanMat{W}(m, n), TanVec{W}(m),
    TanMat{W}(m, m), TanVec{W}(n), TanMat{W}(n, n), TanVec{W}(n), Ref(false), gmask,
    na, dyn, na < n || dyn != 1:n, TanMat{W}(na, na), TanMat{W}(length(dyn), length(dyn)),
    TanMat{W}(length(dyn), length(dyn)), ntd, TanMat{W}(n, ntd), TanMat{W}(n, n), nh, hyb,
    TanMat{W}(nh, nh), wlog, fixed, _tan_setnz!(TanVec{W}(nh), false))

"""
What the blocks of a model with states outside them need (`TanModel.aug`):
the Lyapunov solve and its solution over the diffusing states, JAx's factors
over the leading block of genuine dynamics and its transpose's for the
reverse, the two series at the blocks' sizes, and scratch.
"""
struct TanAugFw{W}
    lya::TanLyap{W}; Xd::TanMat{W}
    FJa::Matrix{Float64}; pivJa::Vector{Int}
    JTa::TanMat{W}; FTa::Matrix{Float64}; pivTa::Vector{Int}
    seria::TanSeries{W}; serna::TanSeries{W}
    Aa::TanMat{W}; Ad::TanMat{W}; kk1::TanMat{W}; kk2::TanMat{W}; Qb::TanMat{W}
    Xbar::TanMat{W}; Adbar::TanMat{W}; JAdbar::TanMat{W}; Qdbar::TanMat{W}
    JAabar::TanMat{W}; Phibar::TanMat{W}
    ra::TanVec{W}; dINTa::TanVec{W}; sbar::TanVec{W}; affa::TanVec{W}
    dD::TanMat{W}
end
TanAugFw(n::Int, na::Int, k::Int, W::Int) = TanAugFw{W}(TanLyap(k, W), TanMat{W}(k, k),
    zeros(na, na), zeros(Int, na), TanMat{W}(na, na), zeros(na, na), zeros(Int, na),
    TanSeries(na, W), TanSeries(k, W), TanMat{W}(na, na),
    (TanMat{W}(k, k) for _ in 1:8)..., TanMat{W}(na, na), TanMat{W}(na, na),
    (TanVec{W}(na) for _ in 1:4)..., TanMat{W}(n, n))

"""JAx's and Qc's blocks into the model, where it has blocks."""
function _tan_aug_refresh!(mod::TanModel)
    mod.aug || return mod
    _tan_gather!(mod.JAa, mod.JAx, 1:mod.na, 1:mod.na)
    _tan_gather!(mod.JAd, mod.JAx, mod.dyn, mod.dyn)
    _tan_gather!(mod.Qd, mod.Qc, mod.dyn, mod.dyn)
    return mod
end

"""
A block's forward state: the Lyapunov solution and the drift's factors, which
every prediction of a linear model shares, with the two numbers the series
test reads (`_series_needed`); the exponentials by interval, a ring of
`_CTSEM_EXP_TABLE` as the filter's `ExpTable`; and an update's scratch by the
number of manifests a row observes.
"""
struct TanForward{W}
    n::Int
    m::Int
    nh::Int
    aug::Union{Nothing,TanAugFw{W}}     # the blocks' own, where there are blocks
    fr::TanFrechet{W}
    ly::TanLyap{W}
    X::TanMat{W}                        # lyap(JAx, Qc)
    FJ::Matrix{Float64}                 # LU of JAx's values
    pivJ::Vector{Int}
    minpivot::Base.RefValue{Float64}    # JAx's smallest pivot
    conditioning::Base.RefValue{Float64}  # the Lyapunov operator's smallest sum
    dts::Vector{Float64}
    exps::Vector{TanMat{W}}
    count::Base.RefValue{Int}
    next::Base.RefValue{Int}
    As::TanMat{W}; nn1::TanMat{W}; nn2::TanMat{W}; v1::TanVec{W}
    rows::Vector{Union{Nothing,TanRowScratch{W}}}
    seri::TanSeries{W}      # the intercept's series, and the noise's
    sern::TanSeries{W}
end
TanForward(n::Int, m::Int, nh::Int, W::Int, na::Int=n, k::Int=n) = TanForward{W}(n, m, nh,
    na < n || k < n ? TanAugFw(n, na, k, W) : nothing, TanFrechet(n, W),
    TanLyap(n, W),
    TanMat{W}(n, n), zeros(n, n), zeros(Int, n), Ref(0.0), Ref(0.0),
    zeros(_CTSEM_EXP_TABLE), [TanMat{W}(n, n) for _ in 1:_CTSEM_EXP_TABLE], Ref(0), Ref(1),
    TanMat{W}(n, n), TanMat{W}(n, n), TanMat{W}(n, n), TanVec{W}(n),
    Union{Nothing,TanRowScratch{W}}[nothing for _ in 1:m], TanSeries(n, W), TanSeries(n, W))

"""
The scratch for a row observing `mo` manifests, built at first use. Its type
stated: built from `W` as a value, it is not inferred, and every row's update
then read its fields dynamically, a quarter of a sweep.
"""
function _tan_row(fw::TanForward{W}, mo::Int)::TanRowScratch{W} where {W}
    r = fw.rows[mo]
    r === nothing || return r
    r = TanRowScratch(fw.n, mo, fw.nh, W)
    fw.rows[mo] = r
    return r
end

"""
What every prediction of a block shares: X = lyap(JAx, Qc) and JAx's factors,
and the series test's two numbers from them; the exponential ring emptied.
"""
function _tan_forward_begin!(fw::TanForward{W}, mod::TanModel{W}) where {W}
    fw.count[] = 0
    fw.next[] = 1
    # A discrete-time step reads none of them.
    mod.ct || return fw
    if mod.aug
        fa = fw.aug::TanAugFw{W}
        _tan_lyap!(fa.Xd, mod.JAd, mod.Qd, fa.lya, false)
        fw.conditioning[] = _series_lyap_conditioning(fa.lya.buf)
        copyto!(fa.FJa, mod.JAa.v)
        _lu_factor_generic!(fa.FJa, fa.pivJa, mod.na)
        fw.minpivot[] = _series_min_pivot(fa.FJa, mod.na)
        return fw
    end
    _tan_lyap!(fw.X, mod.JAx, mod.Qc, fw.ly, false)
    fw.conditioning[] = _series_lyap_conditioning(fw.ly.buf)
    copyto!(fw.FJ, mod.JAx.v)
    _lu_factor_generic!(fw.FJ, fw.pivJ, fw.n)
    fw.minpivot[] = _series_min_pivot(fw.FJ, fw.n)
    return fw
end

"""exp(JAx dt) with tangents, from the ring or computed into it."""
function _tan_exp_cached!(fw::TanForward{W}, mod::TanModel{W}, dt::Float64) where {W}
    @inbounds for e in 1:fw.count[]
        fw.dts[e] == dt && return fw.exps[e]
    end
    _tan_copy!(fw.As, mod.JAx)
    _tan_scale!(fw.As, dt)
    E = _tan_exp!(fw.fr, fw.As)
    slot = fw.next[]
    fw.dts[slot] = dt
    _tan_copy!(fw.exps[slot], E)
    fw.count[] = max(fw.count[], slot)
    fw.next[] = slot == length(fw.dts) ? 1 : slot + 1
    return fw.exps[slot]
end

"""
One prediction over `dt`: `x <- A x + dINT`, `P <- A (P + eps I) A' + dDIFF`,
`A = exp(JAx dt)`, `dINT = JAx \\ (A CINT - CINT)`, `dDIFF = X - A X A'`;
recorded into `rec`. Where a pivot is small against the interval
(`_series_needed`) the filter takes `dINT = Phi CINT` and `dDIFF = V` as
series instead, and so does this (`tangent_series.jl`).
"""
function _tan_predict!(x::TanVec{W}, P::TanMat{W}, rec::TanPredRec{W}, fw::TanForward{W},
    mod::TanModel{W}, dt::Float64, hs::TanHybState{W}) where {W}
    mod.ct || return _tan_predict_dt!(x, P, rec, fw, mod, hs)
    mod.aug && return _tan_predict_aug!(x, P, rec, fw, mod, dt)
    sint = _series_needed(fw.minpivot[], dt)
    snoise = _series_needed(fw.conditioning[], dt)
    rec.sint[] = sint
    rec.snoise[] = snoise
    A = _tan_exp_cached!(fw, mod, dt)
    _tan_copy!(rec.A, A)
    _tan_copy!(rec.x, x)
    _tan_copy!(rec.P, P)
    rec.dt[] = dt
    # The intercept: s = A affine - affine, dINT = JAx \ s; or Phi affine.
    # affine = CINT + (DRIFT - JAx) x, the local affine model's offset.
    _tan_copy!(rec.affine, mod.CINT)
    mod.rows && _tan_mulvec!(rec.affine, mod.DmJ, x, true, true)
    mod.momon[] && _tan_lincomb!(rec.affine, true, rec.affine, true, mod.mc)
    dINT = rec.dINT
    Phi = nothing
    if sint
        Phi = _tan_series_intercept!(fw.seri, mod.JAx, dt)
        _tan_mulvec!(dINT, Phi, rec.affine)
    else
        _tan_mulvec!(dINT, A, rec.affine)
        _tan_lincomb!(dINT, true, dINT, -1.0, rec.affine)
        _tan_solve!(dINT, fw.FJ, fw.pivJ, mod.JAx)
    end
    # The mean.
    _tan_mulvec!(fw.v1, A, x)
    _tan_lincomb!(x, true, fw.v1, true, dINT)
    # The covariance, ridged before the propagation as the filter does.
    _tan_diag!(P, 1e-10)
    _tan_lower_to_upper!(P)
    _tan_mul!(fw.nn1, A, P)
    _tan_mulNT!(fw.nn2, fw.nn1, A)  # A P A'
    # dDIFF = X - (A X) A', or V, added in.
    if snoise
        _tan_comb!(P, ((1.0, fw.nn2), (1.0, _tan_series_noise!(fw.sern, mod.JAx, mod.Qc, dt))))
    else
        _tan_mul!(fw.nn1, A, fw.X)
        _tan_comb!(P, ((1.0, fw.nn2), (1.0, fw.X)))
        _tan_mulNT!(P, fw.nn1, A, -1.0, true)
    end
    _tan_lower_to_upper!(P)
    mod.nh > 0 && _tan_hybrid_predict!(hs, rec.Xh, rec.Xinc, mod.hyb, A, mod.JAx, fw.FJ,
        fw.pivJ, Phi)
    return true
end

"""
One discrete-time step (`_compute_one_step_form!`): `x <- JAx x + CINT`,
`P <- JAx (P + eps I) JAx' + Qc`, recorded into `rec`; the hybrid's columns
`X <- JAx X`, an intercept's (`kind` 3) gaining `e_l`.
"""
function _tan_predict_dt!(x::TanVec{W}, P::TanMat{W}, rec::TanPredRec{W},
    fw::TanForward{W}, mod::TanModel{W}, hs::TanHybState{W}) where {W}
    A = mod.JAx
    _tan_copy!(rec.x, x)
    _tan_copy!(rec.P, P)
    rec.dt[] = 0.0
    # The offset, CINT + (DRIFT - JAx) x.
    _tan_copy!(rec.affine, mod.CINT)
    mod.rows && _tan_mulvec!(rec.affine, mod.DmJ, x, true, true)
    mod.momon[] && _tan_lincomb!(rec.affine, true, rec.affine, true, mod.mc)
    _tan_mulvec!(fw.v1, A, x)
    _tan_lincomb!(x, true, fw.v1, true, rec.affine)
    _tan_diag!(P, 1e-10)
    _tan_lower_to_upper!(P)
    _tan_mul!(fw.nn1, A, P)
    _tan_mulNT!(P, fw.nn1, A)
    _tan_comb!(P, ((1.0, mod.Qc),), true)
    _tan_lower_to_upper!(P)
    if mod.nh > 0
        h = mod.hyb
        _tan_copy!(rec.Xh, hs.X)
        _tan_mul!(hs.X, A, rec.Xh)
        @inbounds for k in eachindex(h.kind)
            h.kind[k] == 3 && (hs.X.v[h.target[k], k] += 1.0)
        end
    end
    return true
end

"""
The stationary moments as the first row's prior (`_ctsem_stationary!`):
`x = -JAx^-1 CINT`, `P = lyap(JAx, Qc)`, the forward's own `X`; the mean
recorded. `false` where they do not exist, as the filter's invalid point.
"""
function _tan_stationary!(x::TanVec{W}, P::TanMat{W}, sw, mod::TanModel{W}) where {W}
    fw = sw.fw
    _tan_copy!(x, mod.CINT)
    _tan_scale!(x, -1.0)
    _tan_solve!(x, fw.FJ, fw.pivJ, mod.JAx)
    _tan_copy!(sw.x0, x)
    _tan_copy!(P, fw.X)
    all(isfinite, x.v) || return false
    return _ctsem_semidefinite(fw.X.v, sw.Lsd, sw.n)
end

"""
One measurement update on the manifests `y` observes (finite), recorded into
`rec`; the row's log likelihood and its tangents into `ll` (a value and W
tangents). A row observing nothing leaves the state as it is. `false` where S
does not factor.
"""
function _tan_update!(x::TanVec{W}, P::TanMat{W}, rec::TanUpdRec{W}, fw::TanForward{W},
    mod::TanModel{W}, y::AbstractVector{Float64}, ll::Vector{Float64},
    hs::TanHybState{W}) where {W}
    n = fw.n
    o = rec.o
    empty!(o)
    # The Gaussian manifests observed; a categorical row's took its own update.
    for i in eachindex(y)
        isfinite(y[i]) && (isempty(mod.gmask) || mod.gmask[i]) && push!(o, i)
    end
    mo = length(o)
    mo == 0 && return true
    rs = _tan_row(fw, mo)
    length(rec.y.v) == mo || (rec.y = TanVec{W}(mo))
    @inbounds for k in 1:mo
        rec.y.v[k] = y[o[k]]
    end
    _tan_setnz!(rec.y, false)
    _tan_copy!(rec.P, P)
    _tan_copy!(rec.x, x)
    mod.nh > 0 && _tan_copy!(rec.Xh, hs.X)
    allstates = 1:n
    H, MM, R = rs.H, rs.MM, rs.R
    # The gain's loading, Jy, which is LAMBDA in a linear model; the mean's,
    # LAMBDA, or Jy too where the update is iterated, its intercept then the
    # linearised one.
    _tan_gather!(H, mod.rows ? mod.Jy : mod.LAMBDA, o, allstates)
    L = H
    if mod.rows && rec.relin == 0
        L = rs.L
        _tan_gather!(L, mod.LAMBDA, o, allstates)
    end
    if rec.relin > 0
        _tan_copy!(MM, rec.c)
    else
        _tan_gather!(MM, mod.MM, o)
    end
    _tan_gather!(R, mod.R, o, o)
    # yt = y - (Lambda x + mu).
    yt = rs.yt
    _tan_copy!(yt, MM)
    _tan_mulvec!(yt, L, x, true, true)
    _tan_lincomb!(yt, true, rec.y, -1.0, yt)
    # S = sym(H Pr H' + R) + ridge, Pr = P + ridge.
    Pr = rs.Pr
    _tan_copy!(Pr, P)
    _tan_diag!(Pr, 1e-10)
    _tan_mulNT!(rs.PHt, Pr, H)
    _tan_mul!(rs.S, H, rs.PHt)
    _tan_comb!(rs.S, ((1.0, R),), true)
    _tan_sym!(rs.Xc, rs.S)
    _tan_copy!(rs.S, rs.Xc)
    _tan_diag!(rs.S, 1e-10)
    _tan_cholinv!(rs.Ui, rs.S, rs.Xc, rs.Tc, rs.F) || return false
    Ui, e, a = rs.Ui, rs.e, rs.a
    # The log likelihood: -1/2 (mo log 2pi + logdet S + e'e), e = Ui' yt,
    # logdet S = -2 sum log Ui_ii.
    _tan_mulTvec!(e, Ui, yt, true, false, Val(:AupperT))
    q = 0.0
    ld = 0.0
    @inbounds for i in 1:mo
        q += e.v[i]^2
        ld -= 2 * log(Ui.v[i, i])
    end
    ll[1] += -0.5 * (mo * log(2pi) + ld + q)
    if _tan_nz(e) || _tan_nz(Ui)
        @inbounds for w in 1:W
            dq = 0.0
            dld = 0.0
            for i in 1:mo
                dq += 2 * e.v[i] * e.d[w, i]
                dld -= 2 * Ui.d[w, i, i] / Ui.v[i, i]
            end
            ll[1 + w] += -0.5 * (dld + dq)
        end
    end
    e.jet && _tan_ll_s!(ll, e, Ui, mo)
    # x <- x + PHt S^-1 yt, S^-1 yt = Ui e.
    _tan_mulvec!(a, Ui, e, true, false, Val(:Aupper))
    _tan_mulvec!(x, rs.PHt, a, true, true)
    # P <- M P M' + K R K', K = PHt S^-1, M = I - K H, on the unridged P.
    _tan_mul!(rs.Kb, rs.PHt, Ui, true, false, Val(:Bupper))
    _tan_mulNT!(rs.K, rs.Kb, Ui, true, false, Val(:BupperT))
    Mm = rs.Mm
    _tan_mul!(Mm, rs.K, H, -1.0, false)
    _tan_diag!(Mm, 1.0)
    _tan_mul!(rs.nn1, Mm, P)
    _tan_mulNT!(P, rs.nn1, Mm)
    _tan_mul!(rs.nm1, rs.K, R)
    _tan_mulNT!(P, rs.nm1, rs.K, true, true)
    _tan_lower_to_upper!(P)
    mod.nh > 0 && _tan_hybrid_update!(hs, rs.hr, mod.hyb, o, H, Ui, e, rs.K)
    return true
end

"""
A row's log likelihood's second parts into `ll[2 + W]`, a jet's: of `-1/2
(e'e - 2 sum_i log Ui_ii)`, `s(e'e) = 2 sum_i (e_i s(e_i) + sum_w de_iw^2)` and
`s(log Ui_ii) = s(Ui_ii) / Ui_ii - sum_w (dUi_iiw / Ui_ii)^2`.
"""
@noinline function _tan_ll_s!(ll::Vector{Float64}, e::TanVec{W}, Ui::TanMat{W}, mo::Int) where {W}
    se, su, de, du = _tan_snz(e), _tan_snz(Ui), _tan_nz(e), _tan_nz(Ui)
    (se || su || de || du) || return nothing
    sq = 0.0
    sl = 0.0
    @inbounds for i in 1:mo
        se && (sq += 2 * e.v[i] * e.s[i])
        de && (sq += 2 * _tan_dot(e.d, W * (i - 1), e.d, W * (i - 1), Val(W)))
        u = Ui.v[i, i]
        su && (sl += Ui.s[i, i] / u)
        if du
            o = _tan_off(Ui.d, Val(W), i, i)
            sl -= _tan_dot(Ui.d, o, Ui.d, o, Val(W)) / (u * u)
        end
    end
    # ll = -1/2 (logdet S + q), logdet S = -2 sum log Ui_ii.
    ll[2 + W] += -0.5 * (-2 * sl + sq)
    return nothing
end

"""
The parameter matrices' cotangents with their tangents, and the exponential's
directions queued by interval (`_frechet_enqueue!`), which a block flushes once
for all its subjects (`_tan_flush_frechet!`).
"""
struct TanAcc{W}
    JAx::TanMat{W}; DRIFT::TanMat{W}; CINT::TanVec{W}; Qc::TanMat{W}
    LAMBDA::TanMat{W}; Jy::TanMat{W}; MM::TanVec{W}; R::TanMat{W}
    T0M::TanVec{W}; T0V::TanMat{W}
    TDE::TanMat{W}; Jtd::TanMat{W}
    Gbar::TanMat{W}             # the hybrid's factor's
    queue::Dict{Float64,Int}    # interval -> its slot in Abars
    Abars::Vector{TanMat{W}}
end
TanAcc(n::Int, m::Int, ntd::Int, nh::Int, W::Int) = TanAcc{W}(TanMat{W}(n, n),
    TanMat{W}(n, n), TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(m, n), TanMat{W}(m, n),
    TanVec{W}(m), TanMat{W}(m, m), TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(n, ntd),
    TanMat{W}(n, n), TanMat{W}(nh, nh), Dict{Float64,Int}(), TanMat{W}[])

function _tan_zero!(X)
    fill!(X.v, 0.0)
    _tan_setnz!(X, false)
    _tan_setsnz!(X, false)
    return X
end

"""Every cotangent zero and the queue empty. Each written out: a loop over
them, matrices and vectors, called `_tan_zero!` dynamically."""
function _tan_acc_zero!(acc::TanAcc)
    _tan_zero!(acc.JAx); _tan_zero!(acc.DRIFT); _tan_zero!(acc.CINT); _tan_zero!(acc.Qc)
    _tan_zero!(acc.LAMBDA); _tan_zero!(acc.Jy); _tan_zero!(acc.MM); _tan_zero!(acc.R)
    _tan_zero!(acc.T0M); _tan_zero!(acc.T0V); _tan_zero!(acc.TDE); _tan_zero!(acc.Jtd)
    _tan_zero!(acc.Gbar)
    empty!(acc.queue)
    return acc
end

"""`acc += b`, the queues merged by interval; both flushed or both not."""
function _tan_acc_add!(acc::TanAcc{W}, b::TanAcc{W}) where {W}
    _tan_comb!(acc.JAx, ((1.0, b.JAx),), true)
    _tan_comb!(acc.DRIFT, ((1.0, b.DRIFT),), true)
    _tan_comb!(acc.Jy, ((1.0, b.Jy),), true)
    _tan_lincomb!(acc.CINT, true, acc.CINT, true, b.CINT)
    _tan_comb!(acc.Qc, ((1.0, b.Qc),), true)
    _tan_comb!(acc.LAMBDA, ((1.0, b.LAMBDA),), true)
    _tan_lincomb!(acc.MM, true, acc.MM, true, b.MM)
    _tan_comb!(acc.R, ((1.0, b.R),), true)
    _tan_lincomb!(acc.T0M, true, acc.T0M, true, b.T0M)
    _tan_comb!(acc.T0V, ((1.0, b.T0V),), true)
    _tan_comb!(acc.TDE, ((1.0, b.TDE),), true)
    _tan_comb!(acc.Jtd, ((1.0, b.Jtd),), true)
    _tan_comb!(acc.Gbar, ((1.0, b.Gbar),), true)
    for (dt, k) in b.queue
        _tan_enqueue!(acc, dt, b.Abars[k], 1.0)
    end
    return acc
end

"""Queue `w Abar`, the exponential's direction at interval `dt`."""
function _tan_enqueue!(acc::TanAcc{W}, dt::Float64, Abar::TanMat{W}, w::Float64) where {W}
    k = get(acc.queue, dt, 0)
    if k == 0
        k = length(acc.queue) + 1
        acc.queue[dt] = k
        length(acc.Abars) < k && push!(acc.Abars, TanMat{W}(size(Abar.v)...))
        _tan_comb!(acc.Abars[k], ((w, Abar),))
    else
        _tan_comb!(acc.Abars[k], ((w, Abar),), true)
    end
    return acc
end

"""The queued directions into JAx's cotangent: `dt L(dt JAx', Abar)` each."""
function _tan_flush_frechet!(acc::TanAcc{W}, fw::TanForward{W}, mod::TanModel{W},
    sT::TanMat{W}) where {W}
    for (dt, k) in acc.queue
        _tan_add!(sT, dt, mod.JAx, false, Val(:T))
        _tan_exp_frechet!(fw.fr, sT, acc.Abars[k])
        _tan_comb!(acc.JAx, ((dt, fw.fr.L),), true)
    end
    empty!(acc.queue)
    return acc
end

"""The reverse of a prediction's scratch, and JAx' with its factors."""
struct TanPredRev{W}
    Abar::TanMat{W}; JAxbar::TanMat{W}; Pnew::TanMat{W}; Ps::TanMat{W}; Pt::TanMat{W}
    nn1::TanMat{W}; nn2::TanMat{W}; Xbar::TanMat{W}; Adbar::TanMat{W}; kk6::TanMat{W}
    JAxdbar::TanMat{W}; Qcdbar::TanMat{W}
    xnew::TanVec{W}; sbar::TanVec{W}; affbar::TanVec{W}
    JT::TanMat{W}
    FT::Matrix{Float64}
    pivT::Vector{Int}
    Phibar::TanMat{W}           # the intercept series' cotangent
end
TanPredRev(n::Int, W::Int) = TanPredRev{W}((TanMat{W}(n, n) for _ in 1:12)...,
    (TanVec{W}(n) for _ in 1:3)..., TanMat{W}(n, n), zeros(n, n), zeros(Int, n),
    TanMat{W}(n, n))

"""
Undo one prediction (`_reverse_predict!`, continuous time, closed forms, every
state diffusing, DRIFT = JAx): `xbar`, `Pbar` in and out, the matrices'
cotangents into `acc` with the subject's weight `w`, the exponential's
direction queued by interval.
"""
function _tan_reverse_predict!(xbar::TanVec{W}, Pbar::TanMat{W}, rec::TanPredRec{W},
    rv::TanPredRev{W}, fw::TanForward{W}, mod::TanModel{W}, acc::TanAcc{W},
    w::Float64, hs::TanHybState{W}) where {W}
    mod.ct || return _tan_reverse_predict_dt!(xbar, Pbar, rec, rv, mod, acc, w, hs)
    mod.aug && return _tan_reverse_predict_aug!(xbar, Pbar, rec, rv, fw, mod, acc, w)
    A, JAx, X = rec.A, mod.JAx, fw.X
    Abar, JAxbar, Ps = rv.Abar, rv.JAxbar, rv.Ps
    sint, snoise, dt = rec.sint[], rec.snoise[], rec.dt[]
    _tan_zero!(JAxbar)
    sint && _tan_zero!(rv.Phibar)
    # The mean: x+ = A x + dINT; and the hybrid's columns, X+ = A X + Xinc.
    _tan_ger!(Abar, true, xbar, rec.x, false)
    # Each branch a call of its own: an argument that is a matrix or `nothing`
    # made it a dynamic call every substep.
    if mod.nh > 0
        if sint
            _tan_hybrid_reverse_predict!(hs, mod.hyb, rec.Xh, rec.Xinc, A, Abar, JAxbar,
                rv.JT, rv.FT, rv.pivT, rv.Phibar)
        else
            _tan_hybrid_reverse_predict!(hs, mod.hyb, rec.Xh, rec.Xinc, A, Abar, JAxbar,
                rv.JT, rv.FT, rv.pivT, nothing)
        end
    end
    _tan_mulTvec!(rv.xnew, A, xbar)
    # The covariance: P+ = A (P + eps I) A' + dDIFF.
    _tan_sym!(Ps, Pbar)
    _tan_copy!(rv.Pt, rec.P)
    _tan_diag!(rv.Pt, 1e-10)
    _tan_mul!(rv.nn1, Ps, A)
    _tan_mul!(Abar, rv.nn1, rv.Pt, 2.0, true)
    _tan_mulTN!(rv.nn2, A, Ps)
    _tan_mul!(rv.Pnew, rv.nn2, A)
    if snoise
        # dDIFF = V, its series run again and pulled back; it reads no A.
        _tan_series_noise!(fw.sern, JAx, mod.Qc, dt)
        _tan_zero!(rv.JAxdbar)
        _tan_zero!(rv.Qcdbar)
        _tan_series_noise_pullback!(rv.JAxdbar, rv.Qcdbar, fw.sern, Ps)
    else
        # dDIFF = X - A X A', its cotangent Ps (symmetric already).
        _tan_mulTN!(rv.kk6, A, Ps)
        _tan_mul!(rv.Xbar, rv.kk6, A)
        _tan_comb!(rv.Xbar, ((1.0, Ps),), -1.0)  # Xbar = Qb - A' Qb A
        _tan_mul!(rv.kk6, Ps, A)
        _tan_mulNT!(rv.Adbar, rv.kk6, X)
        _tan_mulTN!(rv.kk6, Ps, A)
        _tan_mul!(rv.Adbar, rv.kk6, X, true, true)
        _tan_scale!(rv.Adbar, -1.0)  # -(Qb A X' + Qb' A X)
        _tan_comb!(Abar, ((1.0, rv.Adbar),), true)
        # X = lyap(JAx, Qc).
        _tan_lyap_pullback!(rv.JAxdbar, rv.Qcdbar, fw.ly, JAx, X, rv.Xbar)
    end
    if sint
        # dINT = Phi affine, the series run again: affine's cotangent Phi'
        # xbar, and Phi's, with the hybrid's columns', pulled back.
        Phi = _tan_series_intercept!(fw.seri, JAx, dt)
        _tan_mulTvec!(rv.affbar, Phi, xbar)
        _tan_ger!(rv.Phibar, true, xbar, rec.affine, true)
        _tan_series_intercept_pullback!(JAxbar, fw.seri, rv.Phibar)
    else
        # dINT = JAx \ s, s = A affine - affine: sbar = JAx' \ dINTbar.
        _tan_copy!(rv.sbar, xbar)
        _tan_solve!(rv.sbar, rv.FT, rv.pivT, rv.JT)
        _tan_ger!(JAxbar, -1.0, rv.sbar, rec.dINT, true)
        _tan_mulTvec!(rv.affbar, A, rv.sbar)
        _tan_lincomb!(rv.affbar, true, rv.affbar, -1.0, rv.sbar)
        _tan_ger!(Abar, true, rv.sbar, rec.affine, true)
    end
    # affine = CINT + (DRIFT - JAx) x: CINT's cotangent; DRIFT's and JAx's
    # parts cancel where DRIFT is JAx, and otherwise go to each, and to x.
    _tan_lincomb!(acc.CINT, true, acc.CINT, w, rv.affbar)
    if mod.rows
        _tan_ger!(acc.DRIFT, w, rv.affbar, rec.x, true)
        _tan_ger!(JAxbar, -1.0, rv.affbar, rec.x, true)
        _tan_mulTvec!(rv.xnew, mod.DmJ, rv.affbar, true, true)
    end
    _tan_comb!(JAxbar, ((1.0, rv.JAxdbar),), true)
    _tan_comb!(acc.JAx, ((w, JAxbar),), true)
    _tan_comb!(acc.Qc, ((w, rv.Qcdbar),), true)
    _tan_enqueue!(acc, rec.dt[], Abar, w)
    _tan_copy!(xbar, rv.xnew)
    _tan_copy!(Pbar, rv.Pnew)
    return xbar, Pbar
end

"""
Undo one discrete-time step: `JAx`'s, `CINT`'s and `Qc`'s cotangents into
`acc` with the subject's weight `w`; `xbar`, `Pbar` and the hybrid's `Xbar`
carried back.
"""
function _tan_reverse_predict_dt!(xbar::TanVec{W}, Pbar::TanMat{W}, rec::TanPredRec{W},
    rv::TanPredRev{W}, mod::TanModel{W}, acc::TanAcc{W}, w::Float64,
    hs::TanHybState{W}) where {W}
    A, Abar, Ps = mod.JAx, rv.Abar, rv.Ps
    # The offset's cotangent, which a moment-matched step's reverse reads.
    _tan_copy!(rv.affbar, xbar)
    # The mean: x+ = A x + CINT + (DRIFT - JAx) x; the hybrid's columns,
    # X+ = A X (+ e_l).
    _tan_ger!(Abar, true, xbar, rec.x, false)
    _tan_lincomb!(acc.CINT, true, acc.CINT, w, xbar)
    if mod.rows
        _tan_ger!(acc.DRIFT, w, xbar, rec.x, true)
        _tan_ger!(Abar, -1.0, xbar, rec.x, true)
    end
    if mod.nh > 0
        _tan_mulNT!(Abar, hs.Xbar, rec.Xh, true, true)
        _tan_copy!(hs.Xt, hs.Xbar)
        _tan_mulTN!(hs.Xbar, A, hs.Xt)
    end
    _tan_mulTvec!(rv.xnew, A, xbar)
    mod.rows && _tan_mulTvec!(rv.xnew, mod.DmJ, xbar, true, true)
    # The covariance: P+ = A (P + eps I) A' + Qc.
    _tan_sym!(Ps, Pbar)
    _tan_copy!(rv.Pt, rec.P)
    _tan_diag!(rv.Pt, 1e-10)
    _tan_mul!(rv.nn1, Ps, A)
    _tan_mul!(Abar, rv.nn1, rv.Pt, 2.0, true)
    _tan_mulTN!(rv.nn2, A, Ps)
    _tan_mul!(rv.Pnew, rv.nn2, A)
    _tan_comb!(acc.Qc, ((w, Ps),), true)
    _tan_comb!(acc.JAx, ((w, Abar),), true)
    _tan_copy!(xbar, rv.xnew)
    _tan_copy!(Pbar, rv.Pnew)
    return xbar, Pbar
end

"""
Undo the stationary moments: `x = -JAx^-1 CINT` gives `CINT` `-JAx^-T xbar`
and `JAx` that times `x'`; `P = lyap(JAx, Qc)` its pullback, by the forward's
factors. Into `acc`, times `w`.
"""
function _tan_reverse_stationary!(xbar::TanVec{W}, Pbar::TanMat{W}, sw,
    mod::TanModel{W}, acc::TanAcc{W}, w::Float64) where {W}
    rv, fw = sw.rv, sw.fw
    _tan_copy!(rv.sbar, xbar)
    _tan_solve!(rv.sbar, rv.FT, rv.pivT, rv.JT)
    _tan_lincomb!(acc.CINT, true, acc.CINT, -w, rv.sbar)
    _tan_ger!(acc.JAx, -w, rv.sbar, sw.x0, true)
    _tan_sym!(rv.Xbar, Pbar)
    _tan_lyap_pullback!(rv.JAxdbar, rv.Qcdbar, fw.ly, mod.JAx, fw.X, rv.Xbar)
    _tan_comb!(acc.JAx, ((w, rv.JAxdbar),), true)
    _tan_comb!(acc.Qc, ((w, rv.Qcdbar),), true)
    return nothing
end

const _TAN_OK = 0
const _TAN_INVALID = 1     # the filter's invalid point: S did not factor, or worse
const _TAN_REFUSED = 2     # an interval the filter takes as a series

"""A worker's sweep: its records, the walk's state, its scratch."""
struct TanSweep{W}
    n::Int
    m::Int
    fw::TanForward{W}
    rv::TanPredRev{W}
    preds::Vector{TanPredRec{W}}
    upds::Vector{TanUpdRec{W}}
    x::TanVec{W}; P::TanMat{W}; xbar::TanVec{W}; Pbar::TanMat{W}
    sT::TanMat{W}               # dt JAx', the flush's argument
    ll::Vector{Float64}         # the log likelihood, its W tangents, a jet's second part
    first::Vector{Int}          # each row's first prediction record
    nh::Int
    hs::TanHybState{W}
    tdP::Vector{TanMat{W}}      # each row's covariance before its TD impulse
    tdv::TanVec{W}              # a row's TD predictors (no tangent)
    tdn1::TanMat{W}; tdn2::TanMat{W}
    x0::TanVec{W}               # the stationary mean
    Lsd::Matrix{Float64}        # the stationary covariance's check
    # Each row's log likelihood and its tangents, (1 + W) a row, for the row
    # scores (`_tan_row_scores!`), and the running total they difference.
    llrow::Vector{Float64}
    llprev::Vector{Float64}
end
TanSweep(n::Int, m::Int, ntd::Int, nh::Int, W::Int, na::Int=n, k::Int=n) = TanSweep{W}(n, m,
    TanForward(n, m, nh, W, na, k), TanPredRev(n, W),
    TanPredRec{W}[], TanUpdRec{W}[],
    TanVec{W}(n), TanMat{W}(n, n), TanVec{W}(n), TanMat{W}(n, n), TanMat{W}(n, n), zeros(2 + W),
    Int[], nh, TanHybState(n, nh, W), TanMat{W}[], TanVec{W}(ntd), TanMat{W}(n, n),
    TanMat{W}(n, n), TanVec{W}(n), zeros(n, n), Float64[], zeros(2 + W))

"""Row `t`'s log likelihood and its tangents, what the total gained since the
last row, kept in `sw.llrow`."""
function _tan_llrow!(sw::TanSweep{W}, t::Int) where {W}
    o = (t - 1) * (1 + W)
    length(sw.llrow) < o + 1 + W && resize!(sw.llrow, o + 1 + W)
    @inbounds for a in 1:(1 + W)
        sw.llrow[o + a] = sw.ll[a] - sw.llprev[a]
        sw.llprev[a] = sw.ll[a]
    end
    return sw
end

"""
A row's TD impulse (`_apply_td_impulse!`): `x += TDE u`, `P = Jtd P Jtd'`, with
`u` the row's predictors `td`, the covariance before it recorded for row `t`;
and a predictor effect the filter integrates (`kind` 4) gains `u` in its
column at its latent (`_hybrid_td_impulse!`).
"""
function _tan_td!(x::TanVec{W}, P::TanMat{W}, sw::TanSweep{W}, mod::TanModel{W},
    td::AbstractVector{Float64}, t::Int) where {W}
    tdv = sw.tdv
    @inbounds for j in eachindex(td)
        tdv.v[j] = td[j]
    end
    _tan_setnz!(tdv, false)
    t > length(sw.tdP) && push!(sw.tdP, TanMat{W}(sw.n, sw.n))
    _tan_copy!(sw.tdP[t], P)
    _tan_mulvec!(x, mod.TDE, tdv, true, true)
    _tan_mul!(sw.tdn1, mod.Jtd, P)
    _tan_mulNT!(P, sw.tdn1, mod.Jtd)
    if mod.nh > 0
        h = mod.hyb
        @inbounds for k in eachindex(h.kind)
            h.kind[k] == 4 && (sw.hs.X.v[h.target[k], k] += td[h.pred[k]])
        end
    end
    return nothing
end

"""
Undo row `t`'s TD impulse (`_reverse_td!`): `TDE`'s and `Jtd`'s cotangents into
`acc`, times `w`; `xbar` passes through, `Pbar <- Jtd' Ps Jtd`; the hybrid's
columns' cotangent passes through too.
"""
function _tan_reverse_td!(xbar::TanVec{W}, Pbar::TanMat{W}, sw::TanSweep{W},
    mod::TanModel{W}, acc::TanAcc{W}, w::Float64, td::AbstractVector{Float64},
    t::Int) where {W}
    tdv = sw.tdv
    @inbounds for j in eachindex(td)
        tdv.v[j] = td[j]
    end
    _tan_setnz!(tdv, false)
    _tan_ger!(acc.TDE, w, xbar, tdv, true)
    Ps = sw.rv.Ps
    _tan_sym!(Ps, Pbar)
    # Jtd_bar = Ps Jtd P' + Ps' Jtd P, P the covariance before the impulse.
    _tan_mul!(sw.tdn1, Ps, mod.Jtd)
    _tan_mulNT!(sw.tdn2, sw.tdn1, sw.tdP[t])
    _tan_mul!(sw.tdn2, sw.tdn1, sw.tdP[t], true, true)
    _tan_comb!(acc.Jtd, ((w, sw.tdn2),), true)
    _tan_mulTN!(sw.tdn1, mod.Jtd, Ps)
    _tan_mul!(Pbar, sw.tdn1, mod.Jtd)
    return nothing
end

"""
Undo row `rec`'s update, if it observed anything: on the manifests it observed,
its gathered loadings (the gain's `Jy` and the mean's `LAMBDA`), means and
measurement covariance, and their cotangents scattered back into `acc` with
the weight `w`. `false` where S does not factor.
"""
function _tan_reverse_row_update!(xbar::TanVec{W}, Pbar::TanMat{W}, rec::TanUpdRec{W},
    fw::TanForward{W}, mod::TanModel{W}, acc::TanAcc{W}, w::Float64,
    hs::TanHybState{W}) where {W}
    o = rec.o
    isempty(o) && return true
    allstates = 1:fw.n
    rs = _tan_row(fw, length(o))
    u = rs.u
    rel = rec.relin > 0
    _tan_copy!(u.P_in, rec.P)
    _tan_gather!(u.H, mod.rows ? mod.Jy : mod.LAMBDA, o, allstates)
    if mod.rows && !rel
        _tan_gather!(u.Lambda, mod.LAMBDA, o, allstates)
    else
        _tan_copy!(u.Lambda, u.H)
    end
    _tan_gather!(u.R, mod.R, o, o)
    _tan_copy!(u.x, rec.x)
    if rel
        _tan_copy!(u.mu, rec.c)
    else
        _tan_gather!(u.mu, mod.MM, o)
    end
    _tan_copy!(u.y, rec.y)
    _tan_copy!(u.xbar, xbar)
    _tan_copy!(u.Pbar, Pbar)
    # The hook in a branch of its own: a hook that is a NamedTuple or `nothing`
    # made the reverse a dynamic call each row, and its bytes a third of a
    # hybrid member sweep's.
    ok = mod.nh > 0 ? _tan_reverse_update!(u,
        (hs=hs, hr=rs.hr, h=mod.hyb, o=o, Xh=rec.Xh)) : _tan_reverse_update!(u, nothing)
    ok || return false
    _tan_copy!(xbar, u.xbar)
    _tan_copy!(Pbar, u.Pbar)
    _tan_scatter!(acc.R, w, u.Rbar, o, o)
    _tan_scatter!(acc.MM, -w, u.ybarbar, o)
    if rel
        # Iterated: both loadings were Jy, and the intercept
        # c = MANIFESTMEANS + (LAMBDA - Jy) xl, whose cotangent is -ybarbar.
        _tan_scatter!(acc.Jy, w, u.Hbar, o, allstates)
        _tan_scatter!(acc.Jy, w, u.Lbar, o, allstates)
        G = rs.L
        _tan_ger!(G, true, u.ybarbar, rec.xl, false)
        _tan_scatter!(acc.LAMBDA, -w, G, o, allstates)
        _tan_scatter!(acc.Jy, w, G, o, allstates)
    else
        _tan_scatter!(mod.rows ? acc.Jy : acc.LAMBDA, w, u.Hbar, o, allstates)
        _tan_scatter!(acc.LAMBDA, w, u.Lbar, o, allstates)
    end
    return true
end

"""What a block's subjects share: `_tan_forward_begin!`, and JAx' and its factors."""
function _tan_sweep_begin!(sw::TanSweep{W}, mod::TanModel{W}) where {W}
    _tan_forward_begin!(sw.fw, mod)
    mod.ct || return sw
    if mod.aug
        fa = sw.fw.aug::TanAugFw{W}
        _tan_add!(fa.JTa, true, mod.JAa, false, Val(:T))
        copyto!(fa.FTa, fa.JTa.v)
        _lu_factor_generic!(fa.FTa, fa.pivTa, mod.na)
        return sw
    end
    rv = sw.rv
    _tan_add!(rv.JT, true, mod.JAx, false, Val(:T))
    copyto!(rv.FT, rv.JT.v)
    _lu_factor_generic!(rv.FT, rv.pivT, sw.n)
    return sw
end

"""
One subject, weight `w`: `times`, the observations `Y` (m x T, NaN where
missing), its TD predictors `TD` (ntd x T) and the filter's substep rule
(`_ctsem_substeps`). Adds its cotangents, times `w`, into `acc` and queues its
exponential directions there; leaves its log likelihood in `sw.ll`.
`_tan_sweep_begin!` first, once a block. `forward` stops before the reverse,
each row's log likelihood kept; `rec`, with it, records each row's predicted
and updated moments (`TanTraceRecorder`).
"""
function _tan_sweep!(acc::TanAcc{W}, sw::TanSweep{W}, mod::TanModel{W},
    times::AbstractVector{Float64}, Y::AbstractMatrix{Float64},
    TD::AbstractMatrix{Float64}, rule, w::Float64; forward::Bool=false,
    rec=nothing) where {W}
    n, fw, rv = sw.n, sw.fw, sw.rv
    T = length(times)
    # The predictions, row by row: a row's substeps in order, then its update.
    npred = 0
    first = resize!(sw.first, T + 1)
    x, P = sw.x, sw.P
    if mod.stat
        _tan_stationary!(x, P, sw, mod) || return _TAN_INVALID
    else
        _tan_copy!(x, mod.T0M)
        _tan_copy!(P, mod.T0V)
    end
    fill!(sw.ll, 0.0)
    fill!(sw.llprev, 0.0)
    hs = sw.hs
    mod.nh > 0 && _tan_hybrid_start!(hs, mod.hyb)
    for t in 1:T
        first[t] = npred + 1
        if t > 1
            dtrow = times[t] - times[t - 1]
            k = _ctsem_substeps(dtrow, rule, t)
            h = dtrow / k
            for _ in 1:k
                npred += 1
                npred > length(sw.preds) && push!(sw.preds, TanPredRec(n, sw.nh, W))
                _tan_predict!(x, P, sw.preds[npred], fw, mod, h, hs) || return _TAN_REFUSED
            end
        end
        mod.ntd > 0 && _tan_td!(x, P, sw, mod, view(TD, :, t), t)
        t > length(sw.upds) && push!(sw.upds, TanUpdRec(n, sw.m, sw.nh, W))
        sw.upds[t].relin = 0
        rec === nothing || _tan_record!(rec, 1, t, x, P, sw)
        _tan_update!(x, P, sw.upds[t], fw, mod, view(Y, :, t), sw.ll, hs) ||
            return _TAN_INVALID
        forward && _tan_llrow!(sw, t)
        rec === nothing || _tan_record!(rec, 2, t, x, P, sw)
    end
    first[T + 1] = npred + 1
    # The row scores need the forward alone.
    forward && return isfinite(sw.ll[1]) ? _TAN_OK : _TAN_INVALID
    # The closing term, and its reverse, which starts the walk.
    if mod.nh > 0
        (mod.fixed ? _tan_hybrid_closing!(hs, mod.G, mod.wlog, sw.ll, mod.dfix) :
            _tan_hybrid_closing!(hs, mod.G, mod.wlog, sw.ll)) || return _TAN_INVALID
        _tan_comb!(acc.Gbar, ((w, hs.Gbar),), true)
    end
    isfinite(sw.ll[1]) || return _TAN_INVALID

    # The reverse walk.
    xbar, Pbar = _tan_zero!(sw.xbar), _tan_zero!(sw.Pbar)
    for t in T:-1:1
        _tan_reverse_row_update!(xbar, Pbar, sw.upds[t], fw, mod, acc, w, hs) ||
            return _TAN_INVALID
        mod.ntd > 0 && _tan_reverse_td!(xbar, Pbar, sw, mod, acc, w, view(TD, :, t), t)
        for j in (first[t + 1] - 1):-1:first[t]
            _tan_reverse_predict!(xbar, Pbar, sw.preds[j], rv, fw, mod, acc, w, hs)
        end
    end
    if mod.stat
        _tan_reverse_stationary!(xbar, Pbar, sw, mod, acc, w)
    else
        _tan_lincomb!(acc.T0M, true, acc.T0M, w, xbar)
        _tan_sym!(rv.Ps, Pbar)
        _tan_comb!(acc.T0V, ((w, rv.Ps),), true)
    end
    return _TAN_OK
end

"""
`_tan_predict!` where the model has states outside its blocks
(`_compute_discrete_time_form!`'s `naff` and `dyn`): the exponential over every
state, the intercept's solve over the leading block of genuine dynamics,
`dINT = JAa \\ (Aa affine_a - affine_a)` or its series, zero beyond it, and the
noise over the diffusing states, `X - Ad X Ad'` or its series, zero elsewhere.
"""
function _tan_predict_aug!(x::TanVec{W}, P::TanMat{W}, rec::TanPredRec{W},
    fw::TanForward{W}, mod::TanModel{W}, dt::Float64) where {W}
    fa = fw.aug::TanAugFw{W}
    na, dyn = mod.na, mod.dyn
    sint = _series_needed(fw.minpivot[], dt)
    snoise = _series_needed(fw.conditioning[], dt)
    rec.sint[] = sint
    rec.snoise[] = snoise
    A = _tan_exp_cached!(fw, mod, dt)
    _tan_copy!(rec.A, A)
    _tan_copy!(rec.x, x)
    _tan_copy!(rec.P, P)
    rec.dt[] = dt
    _tan_copy!(rec.affine, mod.CINT)
    mod.rows && _tan_mulvec!(rec.affine, mod.DmJ, x, true, true)
    mod.momon[] && _tan_lincomb!(rec.affine, true, rec.affine, true, mod.mc)
    ra, da = fa.ra, fa.dINTa
    _tan_gather!(ra, rec.affine, 1:na)
    if sint
        _tan_mulvec!(da, _tan_series_intercept!(fa.seria, mod.JAa, dt), ra)
    else
        _tan_gather!(fa.Aa, A, 1:na, 1:na)
        _tan_mulvec!(da, fa.Aa, ra)
        _tan_lincomb!(da, true, da, -1.0, ra)
        _tan_solve!(da, fa.FJa, fa.pivJa, mod.JAa)
    end
    dINT = rec.dINT
    _tan_zero!(dINT)
    _tan_scatter!(dINT, 1.0, da, 1:na)
    _tan_mulvec!(fw.v1, A, x)
    _tan_lincomb!(x, true, fw.v1, true, dINT)
    _tan_diag!(P, 1e-10)
    _tan_lower_to_upper!(P)
    _tan_mul!(fw.nn1, A, P)
    _tan_mulNT!(fw.nn2, fw.nn1, A)
    D = fa.dD
    _tan_zero!(D)
    if snoise
        _tan_scatter!(D, 1.0, _tan_series_noise!(fa.serna, mod.JAd, mod.Qd, dt), dyn, dyn)
    else
        _tan_gather!(fa.Ad, A, dyn, dyn)
        _tan_mul!(fa.kk1, fa.Ad, fa.Xd)
        _tan_copy!(fa.kk2, fa.Xd)
        _tan_mulNT!(fa.kk2, fa.kk1, fa.Ad, -1.0, true)
        _tan_scatter!(D, 1.0, fa.kk2, dyn, dyn)
    end
    _tan_comb!(P, ((1.0, fw.nn2), (1.0, D)))
    _tan_lower_to_upper!(P)
    return true
end

"""Undo `_tan_predict_aug!`, as `_tan_reverse_predict!` undoes the full-state one."""
function _tan_reverse_predict_aug!(xbar::TanVec{W}, Pbar::TanMat{W}, rec::TanPredRec{W},
    rv::TanPredRev{W}, fw::TanForward{W}, mod::TanModel{W}, acc::TanAcc{W},
    w::Float64) where {W}
    fa = fw.aug::TanAugFw{W}
    na, dyn = mod.na, mod.dyn
    A = rec.A
    Abar, JAxbar, Ps = rv.Abar, rv.JAxbar, rv.Ps
    sint, snoise, dt = rec.sint[], rec.snoise[], rec.dt[]
    _tan_zero!(JAxbar)
    # The mean, x+ = A x + dINT.
    _tan_ger!(Abar, true, xbar, rec.x, false)
    _tan_mulTvec!(rv.xnew, A, xbar)
    # The covariance, P+ = A (P + eps I) A' + dDIFF.
    _tan_sym!(Ps, Pbar)
    _tan_copy!(rv.Pt, rec.P)
    _tan_diag!(rv.Pt, 1e-10)
    _tan_mul!(rv.nn1, Ps, A)
    _tan_mul!(Abar, rv.nn1, rv.Pt, 2.0, true)
    _tan_mulTN!(rv.nn2, A, Ps)
    _tan_mul!(rv.Pnew, rv.nn2, A)
    # The noise over the diffusing states.
    Qb = fa.Qb
    _tan_gather!(Qb, Ps, dyn, dyn)
    _tan_zero!(fa.JAdbar)
    _tan_zero!(fa.Qdbar)
    if snoise
        _tan_series_noise!(fa.serna, mod.JAd, mod.Qd, dt)
        _tan_series_noise_pullback!(fa.JAdbar, fa.Qdbar, fa.serna, Qb)
    else
        Ad, Xd = fa.Ad, fa.Xd
        _tan_gather!(Ad, A, dyn, dyn)
        _tan_mulTN!(fa.kk1, Ad, Qb)
        _tan_mul!(fa.Xbar, fa.kk1, Ad)
        _tan_comb!(fa.Xbar, ((1.0, Qb),), -1.0)  # Xbar = Qb - Ad' Qb Ad
        _tan_mul!(fa.kk1, Qb, Ad)
        _tan_mulNT!(fa.Adbar, fa.kk1, Xd)
        _tan_mulTN!(fa.kk1, Qb, Ad)
        _tan_mul!(fa.Adbar, fa.kk1, Xd, true, true)
        _tan_scale!(fa.Adbar, -1.0)
        _tan_scatter!(Abar, 1.0, fa.Adbar, dyn, dyn)
        _tan_lyap_pullback!(fa.JAdbar, fa.Qdbar, fa.lya, mod.JAd, Xd, fa.Xbar)
    end
    _tan_scatter!(JAxbar, 1.0, fa.JAdbar, dyn, dyn)
    _tan_scatter!(acc.Qc, w, fa.Qdbar, dyn, dyn)
    # The intercept over the leading block.
    ra, xa, affa = fa.ra, fa.sbar, fa.affa
    _tan_gather!(ra, rec.affine, 1:na)
    _tan_gather!(xa, xbar, 1:na)
    if sint
        Phi = _tan_series_intercept!(fa.seria, mod.JAa, dt)
        _tan_mulTvec!(affa, Phi, xa)
        _tan_ger!(fa.Phibar, true, xa, ra, false)
        _tan_zero!(fa.JAabar)
        _tan_series_intercept_pullback!(fa.JAabar, fa.seria, fa.Phibar)
        _tan_scatter!(JAxbar, 1.0, fa.JAabar, 1:na, 1:na)
    else
        # dINT_a = JAa \ s, s = Aa affine_a - affine_a: sbar = JAa' \ dINTbar_a.
        _tan_solve!(xa, fa.FTa, fa.pivTa, fa.JTa)
        _tan_gather!(fa.dINTa, rec.dINT, 1:na)
        _tan_ger!(fa.JAabar, -1.0, xa, fa.dINTa, false)
        _tan_scatter!(JAxbar, 1.0, fa.JAabar, 1:na, 1:na)
        _tan_gather!(fa.Aa, A, 1:na, 1:na)
        _tan_mulTvec!(affa, fa.Aa, xa)
        _tan_lincomb!(affa, true, affa, -1.0, xa)
        _tan_ger!(fa.Phibar, true, xa, ra, false)
        _tan_scatter!(Abar, 1.0, fa.Phibar, 1:na, 1:na)
    end
    affbar = rv.affbar
    _tan_zero!(affbar)
    _tan_scatter!(affbar, 1.0, affa, 1:na)
    _tan_lincomb!(acc.CINT, true, acc.CINT, w, affbar)
    if mod.rows
        _tan_ger!(acc.DRIFT, w, affbar, rec.x, true)
        _tan_ger!(JAxbar, -1.0, affbar, rec.x, true)
        _tan_mulTvec!(rv.xnew, mod.DmJ, affbar, true, true)
    end
    _tan_comb!(acc.JAx, ((w, JAxbar),), true)
    _tan_enqueue!(acc, dt, Abar, w)
    _tan_copy!(xbar, rv.xnew)
    _tan_copy!(Pbar, rv.Pnew)
    return xbar, Pbar
end
