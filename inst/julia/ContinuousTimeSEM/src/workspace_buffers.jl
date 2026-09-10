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
    piv::Vector{Int}
    dim::Val{N}
    SquareBuffer{T}(n::Int) where {T} = new{T,n}(
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(Int, n),
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
struct ContinuousEKFWorkspace{T, N, M, PARS, BQ, BTHETA, DCA, EBUF, LBUF, DIFBUF, DBUF, DSI, ST, DCACHE, AFBUF}
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
    # 0 Gaussian, 1 binary, 2 ordinal, 3 count, 4 censored. Copied from the
    # `EKFParameters` at
    # construction because the filter's `pars` is the evaluated matrix
    # ComponentVector -- LAMBDA, DRIFT and so on -- and has nowhere for
    # model-level metadata to live. Concrete, so it costs no extra type
    # parameter. Empty means every variable is Gaussian.
    manifesttype::Vector{Int}
    ncategories::Vector{Int}
    # Scratch for one ordinal variable's cumulated thresholds. Length is the
    # widest THRESHOLDS row in the model, zero when there is no such matrix.
    thresholds::Vector{T}
    censormin::Vector{Float64}
    censormax::Vector{Float64}
    # Which covariance construction the model asked for, copied from the
    # `EKFParameters` for the same reason `manifesttype` is: `pars` is the
    # evaluated matrix ComponentVector and has nowhere for model-level metadata.
    # 0 is the unconstrained correlation square root, 2 is covmattransform='z'.
    covmatcode::Int
    # Scratch for the continuous form's intercept solve, sized to the leading
    # block of genuine dynamics (`sp.affine_dim`) rather than to the diffusion
    # subset. Its `dim` field carries that size as a `Val`, so the solve needs
    # nothing else passed alongside it. Separate from `diffusion_buffer`
    # because the two blocks are different sets and the Lyapunov solve wants
    # its own: see `affine_dim` in `EKFParameters`.
    affine_buffer::AFBUF
end

"""
    _init_continuous_ekf_workspace(T, sp)

Allocate a `ContinuousEKFWorkspace` for scalar type `T` and parameter metadata
`sp`.
"""
function _init_continuous_ekf_workspace(::Type{T}, sp::EKFParameters) where {T}
    # Storage for the fully materialized parameter vector (mutable + fixed entries).
    # UNSET_PARAMETER, not undef and not zero: this buffer is filled at mutable
    # and at fixed positions only, so a cell owned solely by a transform is
    # covered by neither. undef makes a read before that transform runs differ
    # between runs and between machines; zero makes it invisible. The sentinel
    # makes it deterministic *and* loud. See `UNSET_PARAMETER` for why 99999.
    all_params = fill(T(UNSET_PARAMETER), length(sp.mutables))
    subject_values = Vector{T}(undef, maximum(vcat(sp.parnumber, sp.ti_coefficient_indices, [0])))

    # Structured parameter view used throughout the EKF loop. Every evaluation
    # overwrites this wholesale from `all_params`, so the fill is only ever
    # read if some path forgets to -- which is the reason it is a sentinel.
    pars = ComponentVector(fill(T(UNSET_PARAMETER), length(sp.mutables)), sp.parameter_axis)

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
    # The leading block of genuine dynamics, over which the intercept solve
    # runs. Every fit from R states it. A caller that does not -- the engine's
    # own direct `EKFParameters` constructions, and anything predating the
    # argument -- falls back to the diffusion block, which is what the offset
    # used before this existed. That matters: such a caller building an
    # `intoverpop`-style model has static carrier states, and running the
    # solve over them inverts a structurally singular `JAx` and returns NaN.
    # The fallback only applies when the diffusion block is itself a leading
    # block, which is exactly the augmented shape; otherwise the whole state
    # vector is the safe reading.
    stated = isdefined(sp, :affine_dim) ? sp.affine_dim : 0
    affine_dim = if stated != 0
        stated
    elseif !isempty(diffusion_state_indices) &&
            diffusion_state_indices == 1:length(diffusion_state_indices)
        length(diffusion_state_indices)
    else
        n
    end
    1 <= affine_dim <= n ||
        throw(ArgumentError("affine_dim must be between 1 and the latent-state dimension"))

    # Reusable matrix/vector work buffers.
    bufferQ = _make_square_buffer(T, n)
    bufferΘ = _make_square_buffer(T, m)
    discrete_ca = _make_discrete_ca_buffer(T, n)
    exp_buffer = ExpBuffer(pars.DIFFUSION)
    lyap_buffer = LyapBuffer(T, length(diffusion_state_indices))
    diffusion_buffer = _make_square_buffer(T, length(diffusion_state_indices))
    discretization_buffer = _make_discretization_buffer(T, n)
    discretization_cache = DiscretizationCache(T, n, length(diffusion_state_indices))
    affine_buffer = _make_square_buffer(T, affine_dim)

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
    ncategories = isdefined(sp, :ncategories) ? copy(sp.ncategories) : Int[]
    censormin = isdefined(sp, :censormin) ? copy(sp.censormin) : Float64[]
    censormax = isdefined(sp, :censormax) ? copy(sp.censormax) : Float64[]
    # At least three, because a censored row borrows this same scratch to carry
    # its lower limit, upper limit and standard deviation -- one row at a time,
    # exactly as an ordinal row borrows it for its cumulated thresholds.
    thresholds = zeros(T, max(hasproperty(pars, :THRESHOLDS) ?
        size(pars.THRESHOLDS, 2) : 0, any(==(4), manifesttype) ? 3 : 0))

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
        ncategories,
        thresholds,
        censormin,
        censormax,
        sp.covmatcode,
        affine_buffer,
    )
end
