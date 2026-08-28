using ComponentArrays
using LinearAlgebra

################################################################################
# Workspace struct and initialization functions for continuous-time EKF
################################################################################
struct SquareBuffer{T,N}
    s::Vector{T}
    ss::Vector{T}
    r::Vector{T}
    row_sq::Vector{T}
    intermediate::Matrix{T}
    out::Matrix{T}
    dim::Val{N}
    SquareBuffer{T}(n::Int) where {T} = new{T,n}(
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n, n),
        zeros(T, n, n),
        Val(n),
    )
end

"""
    _make_square_buffer(T, n)

Create reusable vector and square-matrix scratch storage with element type `T`.

The returned `SquareBuffer` is used by covariance, correlation, and EKF matrix
operations.
"""
@inline function _make_square_buffer(::Type{T}, n::Int) where {T}
    return SquareBuffer{T}(n)
end

"""
    _make_discrete_ca_buffer(T, n)

Create reusable storage for continuous-to-discrete EKF quantities.

The returned `ComponentVector` contains discrete drift, intercept, asymptotic
diffusion, and discrete diffusion matrices.
"""
@inline function _make_discrete_ca_buffer(::Type{T}, n::Int) where {T}
    return ComponentVector(
        eJAx = zeros(T, n, n),
        dDRIFT = zeros(T, n, n),
        dINT = zeros(T, n),
        asym_DIFFUSION = zeros(T, n, n),
        dDIFFUSION = zeros(T, n, n),
    )
end

"""Reusable workspace for exact block-exponential discretization."""
struct DiscretizationBuffer{T,D,EB}
    input::Matrix{T}
    output::Matrix{T}
    scratch::Matrix{T}
    exp_buffer::EB
    dim::Val{D}
end

function _make_discretization_buffer(::Type{T}, n::Int) where {T}
    dim = 2 * n
    return DiscretizationBuffer(zeros(T, dim, dim), zeros(T, dim, dim),
        zeros(T, dim, dim), ExpBuffer{T}(dim), Val(dim))
end

"""
    ContinuousEKFWorkspace

Reusable workspace for continuous-time extended Kalman filter evaluations.

The workspace holds materialized parameters, structured parameter views, matrix
factorizations, covariance buffers, and log-likelihood scratch storage for one
scalar type.
"""
struct ContinuousEKFWorkspace{T, N, M, PARS, BQ, BTHETA, DCA, EBUF, LBUF, DIFBUF, DBUF, DSI, ST, DCACHE}
    all_params::Vector{T}
    subject_values::Vector{T}
    pars::PARS
    predict_param_indices::Vector{Int}
    update_param_indices::Vector{Int}
    td_param_indices::Vector{Int}
    state_dim::Val{N}
    manifest_dim::Val{M}
    bufferQ::BQ
    bufferΘ::BTHETA
    discrete_ca::DCA
    exp_buffer::EBUF
    lyap_buffer::LBUF
    diffusion_buffer::DIFBUF
    discretization_buffer::DBUF
    discretization_cache::DCACHE
    diffusion_state_indices::DSI
    continuous_time::Bool
    ỹ::Vector{T}
    S::Cholesky{T, Matrix{T}}
    K::Matrix{T}
    KR::Matrix{T}
    state::ST
    P_update::Symmetric{T, Matrix{T}}
    P_predict::Symmetric{T, Matrix{T}}
    ll_buffer::Vector{T}
    # 1 for a binary manifest variable, 0 for Gaussian. Copied from the
    # `EKFParameters` at construction because the filter's `pars` is the
    # evaluated matrix ComponentVector -- LAMBDA, DRIFT and so on -- and has
    # nowhere for model-level metadata to live. Concrete, so it costs no extra
    # type parameter. Empty means every variable is Gaussian.
    manifesttype::Vector{Int}
end

"""
    _init_continuous_ekf_workspace(T, sp)

Allocate a `ContinuousEKFWorkspace` for scalar type `T` and parameter metadata
`sp`.
"""
function _init_continuous_ekf_workspace(::Type{T}, sp::EKFParameters) where {T}
    # Storage for the fully materialized parameter vector (mutable + fixed entries).
    all_params = Vector{T}(undef, length(sp.mutables))
    subject_values = Vector{T}(undef, maximum(vcat(sp.parnumber, sp.ti_coefficient_indices, [0])))

    # Structured parameter view used throughout the EKF loop.
    pars = ComponentVector(zeros(T, length(sp.mutables)), sp.parameter_axis)

    # Indices of transforms applied in predict/update phases.
    predict_param_indices = findall(sp.predict_transforms_indices)
    update_param_indices = findall(sp.update_transforms_indices)
    td_param_indices = findall(sp.td_transforms_indices)

    # Latent-state and manifest dimensions.
    n = size(pars.DRIFT, 1)
    m = size(pars.LAMBDA, 1)
    diffusion_state_indices = isempty(sp.diffusion_state_indices) ?
        collect(1:n) : sp.diffusion_state_indices
    all((1 .<= diffusion_state_indices) .& (diffusion_state_indices .<= n)) ||
        throw(ArgumentError("diffusion-state indices are outside the latent-state range"))
    length(unique(diffusion_state_indices)) == length(diffusion_state_indices) ||
        throw(ArgumentError("diffusion-state indices must be unique"))

    # Reusable matrix/vector work buffers.
    bufferQ = _make_square_buffer(T, n)
    bufferΘ = _make_square_buffer(T, m)
    discrete_ca = _make_discrete_ca_buffer(T, n)
    exp_buffer = ExpBuffer(pars.DIFFUSION)
    lyap_buffer = LyapBuffer(T, length(diffusion_state_indices))
    diffusion_buffer = _make_square_buffer(T, length(diffusion_state_indices))
    discretization_buffer = _make_discretization_buffer(T, n)
    discretization_cache = DiscretizationCache(T, n, length(diffusion_state_indices))

    # EKF state for innovation, gain, and covariance updates.
    ỹ = zeros(T, m)
    S = cholesky(Symmetric(Matrix{T}(LinearAlgebra.I, m, m)))
    K = zeros(T, n, m)
    KR = zeros(T, n, m)
    state = similar(pars.T0MEANS)
    fill!(state, zero(T))
    P_update = Symmetric(zeros(T, n, n), :L)
    P_predict = Symmetric(zeros(T, n, n), :L)

    # Scratch vector for log-likelihood solve.
    ll_buffer = zeros(T, m)
    manifesttype = isdefined(sp, :manifesttype) ? copy(sp.manifesttype) : Int[]

    return ContinuousEKFWorkspace(
        all_params,
        subject_values,
        pars,
        predict_param_indices,
        update_param_indices,
        td_param_indices,
        Val(n),
        Val(m),
        bufferQ,
        bufferΘ,
        discrete_ca,
        exp_buffer,
        lyap_buffer,
        diffusion_buffer,
        discretization_buffer,
        discretization_cache,
        diffusion_state_indices,
        sp.continuous_time,
        ỹ,
        S,
        K,
        KR,
        state,
        P_update,
        P_predict,
        ll_buffer,
        manifesttype,
    )
end

"""
    _get_or_init_continuous_ekf_workspace!(ws_ref, params, sp)

Return a cached continuous EKF workspace compatible with `params`.

The cache is refreshed when it is empty or when the parameter scalar type has
changed, for example between `Float64` and `ForwardDiff.Dual`.
"""
function _get_or_init_continuous_ekf_workspace!(ws_ref::Base.RefValue{Any}, params::AbstractVector, sp::EKFParameters)
    ws = ws_ref[]
    T = eltype(params)

    # Reinitialize only when cache is missing, of the wrong workspace type,
    # or built for a different scalar type (e.g., Float64 vs ForwardDiff.Dual).
    if ws === nothing || !(ws isa ContinuousEKFWorkspace) || eltype(ws.all_params) != T
        ws = _init_continuous_ekf_workspace(T, sp)
        ws_ref[] = ws
    end
    return ws::ContinuousEKFWorkspace{T}
end
