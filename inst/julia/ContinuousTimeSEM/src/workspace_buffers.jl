using ComponentArrays
using LinearAlgebra

################################################################################
# Workspace struct and initialization functions for continuous-time EKF
################################################################################
struct SquareBuffer{T}
    s::Vector{T}
    ss::Vector{T}
    r::Vector{T}
    row_sq::Vector{T}
    intermediate::Matrix{T}
    out::Matrix{T}
    piv::Vector{Int}
    dim::Int
    SquareBuffer{T}(n::Int) where {T} = new{T}(
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(Int, n),
        n,
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
struct DiscreteCA{T}
    eJAx::Matrix{T}
    dDRIFT::Matrix{T}
    dINT::Vector{T}
    asym_DIFFUSION::Matrix{T}
    dDIFFUSION::Matrix{T}
end

@inline function _make_discrete_ca_buffer(::Type{T}, n::Int) where {T}
    return DiscreteCA{T}(zeros(T, n, n), zeros(T, n, n), zeros(T, n),
        zeros(T, n, n), zeros(T, n, n))
end

"""Reusable workspace for exact block-exponential discretization."""
struct DiscretizationBuffer{T,EB}
    input::Matrix{T}
    output::Matrix{T}
    scratch::Matrix{T}
    exp_buffer::EB
    dim::Int
    # The series kernels for a pivot small against the interval, and the
    # record of which of them the last discretisation used; see
    # `series_discretization.jl`.
    series::SeriesKernelBuffer{T}
end

function _make_discretization_buffer(::Type{T}, n::Int) where {T}
    dim = 2 * n
    return DiscretizationBuffer(zeros(T, dim, dim), zeros(T, dim, dim),
        zeros(T, dim, dim), ExpBuffer{T}(dim), dim, SeriesKernelBuffer(T, n))
end

"""
    CTSEMGrowth

How far the filter's transitions have carried the state since the last row
with an observation: the 1-norm of the product of every transition since then
(substeps and unobserved rows included), and the largest such norm over the
pass. A row's own interval is not the measure, because nothing conditions the
covariance between observations -- three unobserved weekly rows grow it as one
three-week interval does. Primal values only. See `_CTSEM_EXPLOSIVE_GROWTH`
and `_CTSEM_ADJOINT_GROWTH`.
"""
mutable struct CTSEMGrowth
    max::Float64
    steps::Int
    phi::Matrix{Float64}
    tmp::Matrix{Float64}
end

CTSEMGrowth(n::Int) = CTSEMGrowth(0.0, 0, zeros(n, n), zeros(n, n))

# What a pass returns for the deviations its filter integrates; see
# `CTSEMHybridBuffers`.
const _HYBRID_PROFILE = 0
const _HYBRID_FULL = 1
const _HYBRID_FIXED = 2
# Not integrated at all: a plug-in filter at parameters that already carry the
# deviations (`ctsem_kalman`, `ctsem_generate` on a Laplace fit's subjects).
const _HYBRID_OFF = 3

"""
    CTSEMHybridBuffers{T}

One pass's state for the deviations the filter integrates (`CTSEMHybrid`):
`q` columns, one per hybrid item the subject carries -- a manifest intercept's
where the subject observes the item at all, a state mean's always (`col[j]`
the column of manifest `j`, 0 for none; `item[c]` the hybrid item of column
`c`; `ckind[c]`, `clat[c]`, `cslope[c]` and `cpred[c]` its kind, latent,
slope and, for a predictor's effect, predictor); `X`, the state mean's
dependence on them (`n x q`); `Xinc`, what the last prediction substep added to
a continuous intercept's column, which the tape records; the accumulated `S = sum V' S_t^-1 V` and `s = sum V' S_t^-1 y_t`;
row scratch; and `Jc`, `jv`, `piv` for the intercept Jacobian's solve. Mutable so
the arrays can grow when a model needs more than the empty default -- the
workspace may be built before the Laplace objective fills `CTSEMHybrid`.

`mode` says what the pass returns for them, and is the caller's to set
(`_laplace_subject_value_gradient!`): `_HYBRID_PROFILE`, the default, the
profile -- the log likelihood at the deviations' conditional mode -- which is
what Laplace's mode and curvature are taken of; `_HYBRID_FULL` the profile less
`half_logdet`, their exact integral; `_HYBRID_FIXED` the log likelihood at the
deviations `dfix` (one per column, the caller's), which is how a unit is
evaluated in joint coordinates. Every pass leaves `half_logdet`,
`logdet(M) / 2`; `G`, the `q x rd` factor of the columns' deviations,
`G u` with `u ~ N(0, I)` of dimension `rd` (`gv`, `dgv` its diagonal and
the diagonal's derivative in the raw own sd, for own variances); the
conditional mode `dhat` of `u`; `Mh`, `M` itself; and `r = b - M value(dhat)`,
zero in value, whose partials along a seeded direction are minus that
direction's cross curvature with `u` (`_hybrid_after!` on a
`CTSEMHybridCross`). `dfix`, `dhat`, `r`, `bv` and `M` are over `rd`, which is
`q` for own variances and the level's integrated count for a block; `uitem`
is the hybrid item each coordinate of `u` belongs to, and `_hybrid_ndev` reads
`rd` for a pass that may not have run.
`Lb` is a block's factor, kept with the level parameters it was built at
(`lkey`, `lok`, `lversion`; `_hybrid_block_cached!`).
"""
mutable struct CTSEMHybridBuffers{T}
    active::Bool
    mode::Int
    half_logdet::T
    dfix::Vector{T}
    dhat::Vector{T}
    gv::Vector{T}
    dgv::Vector{T}
    bv::Vector{T}
    r::Vector{T}
    Mh::Matrix{T}
    Mw::Matrix{T}
    G::Matrix{T}
    rd::Int
    uitem::Vector{Int}
    Lb::Matrix{T}
    lkey::Vector{T}
    lok::Bool
    lversion::Int
    q::Int
    col::Vector{Int}
    item::Vector{Int}
    ckind::Vector{Int}
    clat::Vector{Int}
    cslope::Vector{Float64}
    cpred::Vector{Int}
    X::Matrix{T}
    Xt::Matrix{T}
    Xinc::Matrix{T}
    V::Matrix{T}
    W::Matrix{T}
    S::Matrix{T}
    s::Vector{T}
    e::Vector{T}
    Jc::Matrix{T}
    jv::Vector{T}
    piv::Vector{Int}
end
CTSEMHybridBuffers{T}() where {T} = CTSEMHybridBuffers{T}(false, _HYBRID_PROFILE,
    zero(T), T[], T[], T[], T[], T[], T[], zeros(T, 0, 0), zeros(T, 0, 0),
    zeros(T, 0, 0), 0, Int[], zeros(T, 0, 0), T[], false, -1, 0,
    Int[], Int[], Int[], Int[], Float64[], Int[],
    zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0),
    zeros(T, 0, 0), T[], T[], zeros(T, 0, 0), T[], Int[])

"""
    ContinuousEKFWorkspace

Reusable workspace for continuous-time extended Kalman filter evaluations.

The workspace holds materialized parameters, structured parameter views, matrix
factorizations, covariance buffers, and log-likelihood scratch storage for one
scalar type.
"""
struct ContinuousEKFWorkspace{T, PARS, BQ, BTHETA, DCA, EBUF, LBUF, DIFBUF, DBUF, DSI, ST, DCACHE, BPOP, AFBUF, CAT, SP}
    all_params::Vector{T}
    subject_values::Vector{T}
    pars::PARS
    # The model's specification, which the filter reads its model-level
    # metadata from -- the indicators' kinds, categories, asymptote rows and
    # censoring limits, the covariance constructions, the population block --
    # as the adjoint workspace does (`aws.sp`): `pars` is the evaluated
    # matrices, and has nowhere for them.
    sp::SP
    predict_param_indices::Vector{Int}
    update_param_indices::Vector{Int}
    td_param_indices::Vector{Int}
    state_dim::Int
    manifest_dim::Int
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
    # Scratch for one categorical row's extras: an ordinal variable's cumulated
    # thresholds, a binary one's asymptotes, a censored one's limits and sd, a
    # count's dispersion. Long enough for the widest of them in the model.
    thresholds::Vector{T}
    # The row's prior state, where a categorical observation's measurement is
    # linearised (`_ekf_categorical_offset`): the state moves as the row's
    # categorical observations are applied one at a time, and each must be
    # linearised where the Gaussian update linearises the row. Empty for a
    # model with no categorical indicator.
    catprior::Vector{T}
    # Scratch for the population construction, sized to the population block
    # and not to the state dimension. Those are not the same number -- the
    # block is as wide as the number of random effects, and every population
    # row names a distinct state, so it is at most the state dimension and
    # usually less. `sdcovsqrt2cov!` requires its buffer to match the matrix
    # it is given, and the cache in front of it enforces that; sharing
    # `bufferQ` here was a size mismatch.
    population_buffer::BPOP
    # Scratch for the continuous form's intercept solve, sized to the leading
    # block of genuine dynamics (`sp.affine_dim`) rather than to the diffusion
    # subset. Its `dim` field carries that size as a `Val`, so the solve needs
    # nothing else passed alongside it. Separate from `diffusion_buffer`
    # because the two blocks are different sets and the Lyapunov solve wants
    # its own: see `affine_dim` in `EKFParameters`.
    affine_buffer::AFBUF

    # Row scratch for the observed-subset indices: which manifest rows this
    # row has, and the categorical and Gaussian parts of that. One row is
    # live at a time within a subject filter, and a workspace belongs to one
    # task, so one buffer each is enough. Fixed length, used through a view
    # of the part that is filled -- the count varies row to row but the
    # ceiling does not, which is what the old comment here missed when it
    # said a workspace field could not hold this. Half the rows of a real ESM
    # dataset are partially observed (32,468 of 64,784 on the thirteen-study
    # affect model), so this ran on every other row of every evaluation.
    observed_buf::Vector{Int}
    binary_buf::Vector{Int}
    gaussian_buf::Vector{Int}
    # Whether the model has a non-Gaussian indicator (`manifesttype` is not
    # empty), as a type, so that the categorical update, its tape record and
    # its reverse pass are called directly for a model that has one and never
    # compiled for a model that does not; see `_ekf_categorical_call`.
    categorical::Val{CAT}
    # Transition growth since the last observed row; see `CTSEMGrowth`.
    growth::CTSEMGrowth
    # True when no transform group writes a MANIFESTVAR cell, so the manifest
    # covariance is the same on every row of a pass: the filter builds and
    # records it on the first row only (`_extended_kalman_filter_continuous!`),
    # rather than rebuilding a manifest-by-manifest matrix per row and taping
    # it per row for the reverse to compare and accumulate. The reverse's
    # MANIFESTVAR cotangent accumulates over rows until a `:theta` record
    # (`_ctsem_reverse_tape!`), so one record pushes the whole pass back once.
    manifestvar_static::Bool
    # The same for DIFFUSION: when no transform group writes it, the subject
    # filter builds its covariance into `diffusion_cov` once per pass and hands
    # that to every prediction, rather than each prediction rebuilding it.
    diffusion_static::Bool
    diffusion_cov::Matrix{T}
    # The deviations this filter integrates itself, when the model has any
    # (`CTSEMHybrid`, `_hybrid_start!`).
    hyb::CTSEMHybridBuffers{T}
    # The predictor row of a subject with sampled TI predictor cells, built here
    # each evaluation (`_ctsem_tipred_vector(tipreds, p, buffer)`); unused, and
    # empty, for every other subject.
    tipred_row::Vector{T}
    # The last pass's log likelihood, which `_ekf_run` returns through here
    # rather than out of its function barrier; see there.
    loglik::Base.RefValue{T}
    # Whether a Jy cell reads the state, so that the measurement is nonlinear
    # in it and each row's Gaussian update is iterated to the posterior mode
    # (`_ekf_iterated_rows!`). False for a linear measurement, which runs the
    # update it always did; `ctsem_set_relinearise!` is the oracle.
    relinearise::Bool
    # That iteration's prior mean, next iterate and measurement offset at the
    # linearisation point; empty where it does not run.
    xprior::Vector{T}
    xnext::Vector{T}
    mm_eff::Vector{T}
    # The row's iterates, a column each, while a plan is learnt
    # (`_ekf_relin_need`); empty where the update does not iterate.
    relin_hist::Matrix{T}
    # Whether a predict-group transform that writes a JAx cell reads the
    # state, so that the dynamics are nonlinear in it and each prediction
    # substep is moment-matched (`_ekf_moment_cells!`); `ctsem_set_moment_predict!`
    # is the oracle. Whether this substep's prediction was, for its record.
    moment::Bool
    moment_now::Base.RefValue{Bool}
    # The predict group's JAx, DRIFT and CINT cells, averaged over the
    # cubature points, and its other cells, which keep their values at the
    # mean; then the mean, the ridged covariance's upper Cholesky factor, the
    # points' offsets from the mean (a column each), the drift's mean
    # correction, scratch for the averages and the other cells' values, and
    # the drift at each point (`n x n` blocks side by side), which the reverse
    # reads. Empty where the dynamics are linear.
    moment_avg::Vector{Int}
    moment_oth::Vector{Int}
    moment_m::Vector{T}
    moment_U::Matrix{T}
    moment_d::Matrix{T}
    moment_c::Vector{T}
    moment_acc::Vector{T}
    moment_saved::Vector{T}
    moment_drifts::Matrix{T}
end

"""
    _ekf_has_categorical(ws)

Whether the workspace's model has a non-Gaussian indicator: a constant
wherever the workspace's type is known.
"""
@inline _ekf_has_categorical(ws::ContinuousEKFWorkspace) = _val_bool(ws.categorical)
@inline _ekf_has_categorical(ws) = !isempty(ws.sp.manifesttype)
@inline _val_bool(::Val{C}) where {C} = C

"""
    _ekf_categorical_call(ws, f, args...)

`f(args...)`, the categorical code a workspace's model needs, called by what
the workspace's type says: directly when it has a non-Gaussian indicator,
not at all when it has none, and through `_ctsem_barrier` when the type is
not known, which is only ever inference looking at a caller that holds the
workspace abstractly.

Generated for that last case. A plain branch on `_ekf_has_categorical` folds
for a concrete workspace, but inference also walks callers whose workspace is
only `ContinuousEKFWorkspace{T}`, and through a branch it cannot fold it
specialised the categorical files for a Gaussian model: 188 method instances
on a Gaussian augmented fit (bench cf_gaussian, local, no package image),
against 22 behind the barrier alone. A generator is not run for a type that is
not concrete, so inference there stops at the call. And the barrier alone, in
front of every categorical row, cost one dynamic dispatch per row on each of
these paths: 7-10% of ord4's Laplace and continuation gradients (local
profile), where this costs nothing.
"""
@generated function _ekf_categorical_call(ws::ContinuousEKFWorkspace, f::F,
    args...) where {F}
    flag = fieldtype(ws, :categorical)
    flag === Val{true} && return :(f(args...))
    flag === Val{false} && return :(nothing)
    return :(_ctsem_barrier(f, args...))
end
_ekf_categorical_call(ws, f::F, args...) where {F} = _ctsem_barrier(f, args...)

"""
    _ctsem_affine_dim(sp, diffusion_state_indices, n)

The leading block of genuine dynamics, over which the intercept solve runs.

Every fit from R states it. A caller that does not -- the engine's own direct
`EKFParameters` constructions, and anything predating the argument -- falls
back to the diffusion block, which is what the offset used before this existed.
That matters: such a caller building an `intoverpop`-style model has static
carrier states, and running the solve over them inverts a structurally singular
`JAx` and returns NaN. The fallback only applies when the diffusion block is
itself a leading block, which is exactly the augmented shape; otherwise the
whole state vector is the safe reading.

A function rather than inline in the workspace constructor because the summary
needs the same block for the stationary mean it reports.
"""
function _ctsem_affine_dim(sp::EKFParameters, diffusion_state_indices, n::Int)
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
    return affine_dim
end

"""
    _ctsem_check_stationary(sp, diffusion_state_indices, affine_dim)

Refuse a stationary specification the engine cannot honour, at construction
rather than as a wrong likelihood later. R refuses all three first and says
what to do instead; this is the backstop for a caller that does not.
"""
function _ctsem_check_stationary(sp::EKFParameters, diffusion_state_indices,
        affine_dim::Int)
    se = sp.stationary_elements
    elements = _stationary_elements_on(se)
    (sp.stationary || elements) || return nothing
    sp.stationary && elements && throw(ArgumentError(
        "stationary = TRUE starts every latent from its stationary moments; " *
        "it takes no single stationary elements as well"))
    sp.continuous_time || throw(ArgumentError(
        "stationary initial conditions are implemented for continuous time only"))
    all(<=(affine_dim), diffusion_state_indices) || throw(ArgumentError(
        "stationary initial conditions need every diffusing state among the " *
        "leading $(affine_dim) latent processes"))
    if elements
        length(se.inclosure) <= affine_dim || throw(ArgumentError(
            "stationary elements are for the leading $(affine_dim) latent processes"))
        # A random initial mean is a population row on its latent; a stationary
        # mean leaves nothing to vary, and a stationary variance is not one the
        # population block adds to.
        marked = vcat(se.means, reduce(vcat, se.blocks; init=Int[]))
        any(in(marked), sp.population_indices) && throw(ArgumentError(
            "a stationary initial condition's latent has a random initial mean"))
        return nothing
    end
    # An individually varying T0MEANS is a population row on a latent itself.
    # Stationarity leaves no T0MEANS to vary, and the placement would write
    # over the block that row describes.
    any(<=(affine_dim), sp.population_indices) && throw(ArgumentError(
        "stationary initial conditions leave no T0MEANS to vary; a population " *
        "row refers to a latent process"))
    return nothing
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
    pars = CTSEMPars(fill(T(UNSET_PARAMETER), length(sp.mutables)), sp.parameter_axis)

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
    affine_dim = _ctsem_affine_dim(sp, diffusion_state_indices, n)
    _ctsem_check_stationary(sp, diffusion_state_indices, affine_dim)

    # Reusable matrix/vector work buffers.
    bufferQ = _make_square_buffer(T, n)
    population_buffer = _make_square_buffer(T, length(sp.population_indices))
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
    manifesttype = sp.manifesttype
    # At least three, because a censored row borrows this same scratch to carry
    # its lower limit, upper limit and standard deviation (five with a response
    # scale: its half-width and power too) -- one row at a time,
    # exactly as an ordinal row borrows it for its cumulated thresholds. At
    # least one for a count, which borrows it for its dispersion.
    thresholds = zeros(T, max(hasproperty(pars, :THRESHOLDS) ?
        size(pars.THRESHOLDS, 2) : 0,
        any(==(4), manifesttype) ? (hasproperty(pars, :RESPONSE) ? 5 : 3) : 0,
        any(==(3), manifesttype) ? 1 : 0, any(>(0), sp.asymptoterow) ? 2 : 0))
    relinearise = _ctsem_measurement_nonlinear(sp, update_param_indices, n)
    moment = _ctsem_dynamics_nonlinear(sp, predict_param_indices, n)
    moment_avg, moment_oth = moment ? _ctsem_moment_cells(sp, predict_param_indices) :
        (Int[], Int[])
    nm = moment ? n : 0

    return ContinuousEKFWorkspace(
        all_params,
        subject_values,
        pars,
        sp,
        predict_param_indices,
        update_param_indices,
        td_param_indices,
        n,
        m,
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
        thresholds,
        zeros(T, isempty(manifesttype) ? 0 : n),
        population_buffer,
        affine_buffer,
        Vector{Int}(undef, m),
        Vector{Int}(undef, m),
        Vector{Int}(undef, m),
        Val(!isempty(manifesttype)),
        CTSEMGrowth(n),
        _ctsem_static_component(sp, :MANIFESTVAR, predict_param_indices,
            update_param_indices, td_param_indices),
        _ctsem_static_component(sp, :DIFFUSION, predict_param_indices,
            update_param_indices, td_param_indices),
        zeros(T, n, n),
        CTSEMHybridBuffers{T}(),
        T[],
        Ref(zero(T)),
        relinearise,
        zeros(T, relinearise ? n : 0),
        zeros(T, relinearise ? n : 0),
        zeros(T, relinearise ? m : 0),
        zeros(T, relinearise ? n : 0, relinearise ? _CTSEM_RELIN_MAXIT + 1 : 0),
        moment,
        Ref(false),
        moment_avg,
        moment_oth,
        zeros(T, nm),
        zeros(T, nm, nm),
        zeros(T, nm, 2nm),
        zeros(T, nm),
        zeros(T, length(moment_avg)),
        zeros(T, length(moment_oth)),
        zeros(T, nm, 2nm * nm),
    )
end

"""
    _ctsem_static_component(sp, name, predict, update, td)

Whether no transform group -- the per-row predict, update and TD-predictor
groups -- writes a cell of matrix `name`, so that matrix is the same on every
row of a pass; see `ContinuousEKFWorkspace`'s `manifestvar_static` and
`diffusion_static`. The test `CTSEMAdjointWorkspace` makes for
`groups_write_manifestvar` and `groups_write_diffusion`, negated.
"""
function _ctsem_static_component(sp::EKFParameters, name::Symbol, predict, update, td)
    isempty(predict) && isempty(update) && isempty(td) && return true
    positions = Set(_ctsem_component_positions(sp, name))
    return !any(i -> i in positions, Iterators.flatten((predict, update, td)))
end
