using LinearAlgebra
using ComponentArrays

################################################################################
# Reverse-mode EKF
################################################################################
#
# The forward filter (`_extended_kalman_filter_continuous!`) optionally records
# a *tape* while it runs; this file replays that tape backwards.
#
# The tape is deliberately an ordered program of small typed records rather
# than a fixed-shape array-of-checkpoints. The forward pass is not a fixed
# sequence: the number of prediction substeps depends on `max_timestep` and the
# gap between observations, TD impulses only happen when there are TD
# predictors, state-dependent transform groups only exist for nonlinear models,
# and a fully missing row skips the measurement update entirely. Recording the
# order the forward pass actually took, rather than assuming a canonical one,
# is what lets one reverse pass cover every case the tests check --
# linear, state/PARS-expression, partially missing, fully missing, and
# unequal-length multi-subject models -- without special-casing any of them.
#
# Conventions used throughout the reverse pass:
#
#   * `x̄` and `P̄` are the running cotangents of "the state and covariance at
#     this point in the forward pass". Every reverse step consumes them as the
#     cotangent of its outputs and leaves them as the cotangent of its inputs.
#     Identity copies in the primal (`P_update .= P_predict` between substeps,
#     and on a fully missing row) therefore need no reverse code at all.
#   * `θ̄` is the cotangent of `all_params`, in exactly the layout `ws.pars`
#     views. Matrix cotangents are accumulated into it through a
#     `ComponentVector` view, so `θ̄ca.DRIFT` and friends name the right cells.
#   * `Θ̄` is the running cotangent of the *manifest covariance* Θ, which is
#     accumulated by measurement updates and consumed by the `:theta` tape
#     entry at the point where the forward pass built Θ from MANIFESTVAR. It is
#     kept separate from `θ̄` precisely because that build happens at a
#     different point in the row for the first row than for later rows.
#   * Covariance cotangents are symmetrised before use. The primal maintains
#     symmetry explicitly (`_copy_lower_to_upper!`, `_symmetrize_and_ridge!`),
#     and those operations are the identity on symmetric perturbations, so
#     keeping the cotangents symmetric is what makes it valid to treat them as
#     such.
#   * The `+1e-10` ridges Stan applies (and this port matches) are affine, so
#     they are invisible to the reverse pass except through the *value* at
#     which downstream derivatives are evaluated -- which is why the reverse
#     recomputes `Pr = P + 1e-10 I` rather than reusing the unridged `P`.

const _CTSEM_RIDGE = 1e-10

# Distinct `JAx * dt` keys the reverse pass holds Fréchet directions for at
# once. A shared wave schedule with nineteen distinct intervals batched to one
# Fréchet evaluation per substep under a last-value rule -- 3800 per gradient
# on 200 subjects of 20 rows -- and batches to nineteen with a table. Fully
# irregular times fill the table every 32 substeps and flush it, which costs
# what the last-value rule did plus a scan of 32 scalars per substep.
const _CTSEM_FRECHET_TABLE = 32

################################################################################
# Tape records
################################################################################

"""One prediction substep: everything the reverse pass needs to undo it."""
mutable struct CTSEMPredictRecord{T}
    state_in::Vector{T}      # state entering the substep
    P_in::Matrix{T}          # posterior covariance entering the substep (unridged)
    A::Matrix{T}             # eJAx = exp(JAx * dt), also the discrete drift
    JAx::Matrix{T}
    DRIFT::Matrix{T}
    DIFFUSION::Matrix{T}     # raw SD/correlation-sqrt parameters
    Xlyap::Matrix{T}         # Lyapunov solution on the dynamic block
    affine::Vector{T}        # `r` in the derivation: the local affine correction
    dINT_dynamic::Vector{T}  # solved discrete intercept on the dynamic block
    dt::T
    # Which pieces took the series route (`series_discretization.jl`), and the
    # diffusion block the noise series integrated, which its pullback needs
    # where the closed form needed only `Xlyap`, and each series kernel's
    # doubling chain, which its pullback reads (`CTSEMSeriesChain`). Empty
    # until a substep uses them, so a model that never takes the series
    # carries no extra storage.
    series_intercept::Bool
    series_noise::Bool
    Qd::Matrix{T}
    noise_chain::CTSEMSeriesChain{T}
    icpt_chain::CTSEMSeriesChain{T}
    # `X` entering the substep (`CTSEMHybrid`), `n x q`, and what the substep
    # added to a continuous intercept's column; `0 x 0` without.
    Xh::Matrix{T}
    Xinc::Matrix{T}
    # Whether the substep was moment-matched (`_ekf_moment_cells!`): its JAx,
    # DRIFT and CINT are the cubature's averages and `affine` carries the
    # mean's correction, whose cotangent the reverse leaves for
    # `_reverse_moment!` (`_moment_cbar!`).
    moment::Bool
end

"""
One moment-matched prediction substep's cubature (`_ekf_moment_cells!`): the
tape index of its first point's group record (the points' records follow it,
off the program, since this record's reverse walks them), how many points,
the ridged covariance's upper Cholesky factor (its upper triangle), the
points' offsets from the mean (a column each), and the drift at each point
(`n x n` blocks side by side).
"""
mutable struct CTSEMMomentRecord{T}
    first::Int
    npoints::Int
    U::Matrix{T}
    D::Matrix{T}
    DRIFTs::Matrix{T}
end

"""
Scratch for `_reverse_moment!`, sized to the state where the model's dynamics
are moment-matched and empty elsewhere, except `cbar`, the mean correction's
cotangent, which the prediction reverse writes on every moment-matched
substep. `drift_rc` lists the averaged cells that are DRIFT's, as (position
among them, row, column).
"""
struct CTSEMMomentScratch{T}
    cbar::Vector{T}
    cellbar::Vector{T}
    savedoth::Vector{T}
    xs::Vector{T}
    dbar::Vector{T}
    mbar::Vector{T}
    Ubar::Matrix{T}
    B::Matrix{T}
    Z::Matrix{T}
    drift_rc::Vector{NTuple{3,Int}}
end

function CTSEMMomentScratch(::Type{T}, sp, ws) where {T}
    n = ws.state_dim
    nm = ws.moment ? n : 0
    avg = ws.moment_avg
    drift = _ctsem_component_positions(sp, :DRIFT)
    nr = size(ws.pars.DRIFT, 1)
    rc = NTuple{3,Int}[]
    for (k, p) in enumerate(avg)
        q = findfirst(==(p), drift)
        q === nothing && continue
        push!(rc, (k, (q - 1) % nr + 1, (q - 1) ÷ nr + 1))
    end
    return CTSEMMomentScratch{T}(zeros(T, n), zeros(T, length(avg)),
        zeros(T, length(ws.moment_oth)), zeros(T, nm), zeros(T, nm), zeros(T, nm),
        zeros(T, nm, nm), zeros(T, nm, nm), zeros(T, nm, nm), rc)
end

"""One TD-predictor impulse."""
mutable struct CTSEMTDRecord{T}
    P_in::Matrix{T}
    Jtd::Matrix{T}
    tdpreds::Vector{T}
end

"""
One measurement update.

Stores only inputs, not intermediates: the innovation, gain and Cholesky
factor are cheap to recompute from these and doing so keeps the primal's
`_ekf_masked_update_step!` free of tracing code (it reuses `ws.K` as scratch
partway through, so its intermediates are not all live at the end anyway).
"""
mutable struct CTSEMUpdateRecord{T}
    observed::Vector{Int}
    state_in::Vector{T}
    P_in::Matrix{T}          # prior covariance, unridged
    Lambda::Matrix{T}        # LAMBDA[observed, :]
    manifestmeans::Vector{T} # MANIFESTMEANS[observed]
    H::Matrix{T}             # Jy[observed, :]
    R::Matrix{T}             # Θ[observed, observed]
    y::Vector{T}             # observed data for this row
    # `X` entering the update (`CTSEMHybrid`), `n x q`; `0 x 0` without.
    Xh::Matrix{T}
    # A relinearised update (`_ekf_iterated_rows!`): 0 for the EKF's, 1 when
    # linearised at the prior mean `state_in` itself, 2 at a later iterate.
    # `Lambda` and `manifestmeans` are then the loading and intercept it used,
    # Jy and `MANIFESTMEANS + (LAMBDA - Jy) x_lin`, and `Lambda_lin` is LAMBDA
    # at the linearisation point `x_lin`.
    relin::Int
    x_lin::Vector{T}
    Lambda_lin::Matrix{T}
end

"""
One step of an iterated update (`_ekf_relin_mean!`): the mean of the update
linearised at `x_lin`, from the row's prior `state_in` and `P_in` (unridged),
over the observed rows -- `H` the Jacobian there, `offset` the intercept
`MANIFESTMEANS + (LAMBDA - H) x_lin`, `Lambda` LAMBDA there -- and the step
`x_lin + alpha (mean - x_lin)` taken towards it. `first` when `x_lin` is the
prior mean itself.
"""
mutable struct CTSEMRelinRecord{T}
    observed::Vector{Int}
    state_in::Vector{T}
    P_in::Matrix{T}
    H::Matrix{T}
    offset::Vector{T}
    Lambda::Matrix{T}
    x_lin::Vector{T}
    R::Matrix{T}
    y::Vector{T}
    first::Bool
    alpha::Float64
end

"""
    CTSEMHybridTape{T}

The deviations the filter integrated (`CTSEMHybrid`), for the reverse pass:
what the pass ended with -- `q`, the columns' items and manifests (`item`,
`col`), the accumulated `S` and `s`, the factor `G` (`q x rd`), for own
variances its diagonal `g` and `dg/draw`, for a block the level's free raw
values `lvals` -- and the reverse state: the cotangents of `S` and `s` from
the closing term, the cotangent of `X` carried backwards through the rows,
`Gbar`, and from it each item's own-sd cotangent (`psibar`) or the level's
(`lpbar`, over `CTSEMHybrid`'s `lparams`), which the caller adds to the raw
gradient; and a block factor's derivatives in the level parameters, kept with
the values they were taken at (`dkey`, `dLb`, `dok`, `dversion`;
`_hybrid_block_derivatives!`). Scratch for
the per-row work beside it. Every array is sized by capacity -- the model's
hybrid item count, the manifest count, the state count -- and grows only, so a
pass after the first allocates nothing; loops read the leading block they need.
"""
mutable struct CTSEMHybridTape{T}
    active::Bool
    mode::Int
    dfix::Vector{T}
    G::Matrix{T}
    rd::Int
    block::Bool
    lvals::Vector{T}
    Gbar::Matrix{T}
    lpbar::Vector{T}
    dkey::Vector{T}
    dLb::Array{T,3}
    dok::Vector{Bool}
    dversion::Int
    q::Int
    item::Vector{Int}
    col::Vector{Int}
    S::Matrix{T}
    s::Vector{T}
    g::Vector{T}
    dg::Vector{T}
    Sbar::Matrix{T}
    sbar::Vector{T}
    Xbar::Matrix{T}
    psibar::Vector{T}
    V::Matrix{T}
    W::Matrix{T}
    Z::Matrix{T}
    Vbar::Matrix{T}
    Kbar::Matrix{T}
    Rbar::Matrix{T}
    ybar::Vector{T}
    Lbar::Matrix{T}
    Xt::Matrix{T}
    Mw::Matrix{T}
    C::Matrix{T}
    wv::Vector{T}
    uv::Vector{T}
    beta::Vector{T}
end
CTSEMHybridTape{T}() where {T} = CTSEMHybridTape{T}(false, 0, T[], zeros(T, 0, 0), 0,
    false, T[], zeros(T, 0, 0), T[], T[], zeros(T, 0, 0, 0), Bool[], -1, 0, Int[], Int[],
    zeros(T, 0, 0), T[], T[], T[], zeros(T, 0, 0), T[], zeros(T, 0, 0), T[],
    zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0),
    zeros(T, 0, 0), T[], zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0),
    zeros(T, 0, 0), T[], T[], T[])

"""
One application of a group of state-dependent transforms.

`params_before` holds the parameter values from *before* the group ran, and
both halves of that matter.

**Before, not after.** Within a group the transforms are applied in ascending
flattened-index order into the same buffer they read from, so a transform can
read a cell that a later transform in the same group is about to overwrite --
and when it does, it sees the value from the *previous* row, not this one.
(`PARS` sorts last in the parameter axis, so a state-dependent `DRIFT`
expression referencing a state-dependent `PARS` cell is exactly this case.)
Recording the pre-group state lets the reverse pass reconstruct, for each
transform, the values that transform actually saw; recording the post-group
state would silently evaluate some derivatives at the wrong point.

**Only the relevant cells, not the whole vector.** A group's transforms read
only the indices in their recorded read sets and write only their own
`indices`; every other entry of `all_params` is dead as far as this group is
concerned. Storing the lot cost ~23 KB per group per row on a 20-latent model
and made state-dependent models allocate ~50 MB per gradient. `params_before`
therefore holds just the values at `tape.group_relevant[group]`, in that
order.
"""
mutable struct CTSEMGroupRecord{T}
    group::Int               # 1 = predict, 2 = td, 3 = update
    params_before::Vector{T} # values at `group_relevant[group]`, pre-group
    state::Vector{T}
    tdpreds::Vector{T}
    time::T
    dt::T
    row::Int
end

"""Construction of Θ from MANIFESTVAR."""
mutable struct CTSEMThetaRecord{T}
    MANIFESTVAR::Matrix{T}
end

"""The `t = 1` prior: state from T0MEANS, covariance from T0VAR."""
mutable struct CTSEMInitRecord{T}
    T0VAR::Matrix{T}
    # The raw population covariance matrix, at the point the forward pass
    # constructed it. `0` by `0` for a model whose population covariance is not
    # a matrix of its own, which is every model that does not use intoverpop.
    POPRAW::Matrix{T}
end

"""
The stationary initial moments, `_ctsem_stationary!`: the drift over the
leading block of genuine dynamics, the mean it gave, the covariance over the
diffusion states, and the raw DIFFUSION, which keys the deferred diffusion
pullback exactly as a prediction's does. For single stationary elements
(`_ctsem_stationary_elements_cov!`) also the initial covariance before them,
and the packed factors and congruence they used; `0` by `0` otherwise.
"""
mutable struct CTSEMStationaryRecord{T}
    DRIFT::Matrix{T}
    mean::Vector{T}
    X::Matrix{T}
    DIFFUSION::Matrix{T}
    Pold::Matrix{T}
    LX::Matrix{T}
    LP::Matrix{T}
    M::Matrix{T}
end

"""
    CTSEMAdjointTape

Ordered record of one subject's forward pass.

`program` lists `(kind, index)` pairs in forward order; the reverse pass walks
it backwards and dispatches on `kind` into the matching record vector.
"""
mutable struct CTSEMAdjointTape{T}
    program::Vector{Tuple{Symbol,Int}}
    predicts::Vector{CTSEMPredictRecord{T}}
    tds::Vector{CTSEMTDRecord{T}}
    updates::Vector{CTSEMUpdateRecord{T}}
    groups::Vector{CTSEMGroupRecord{T}}
    thetas::Vector{CTSEMThetaRecord{T}}
    inits::Vector{CTSEMInitRecord{T}}
    binaries::Vector{CTSEMBinaryRecord{T}}
    stationaries::Vector{CTSEMStationaryRecord{T}}
    # How many of each vector the *current* pass has written. The vectors are
    # never emptied, so a subject after the first writes into records that
    # already exist and already have arrays of the right shape -- see
    # `_tape_reset!`. `program` is the exception: `empty!` keeps a Vector's
    # capacity, so pushing into it again allocates nothing.
    npredicts::Int
    ntds::Int
    nupdates::Int
    ngroups::Int
    nthetas::Int
    ninits::Int
    nbinaries::Int
    nstationaries::Int
    subject_values::Vector{T}
    # The integrated deviations' record and reverse state; see `CTSEMHybridTape`.
    hyb::CTSEMHybridTape{T}
    # For each transform group (1 = predict, 2 = td, 3 = update), the
    # `all_params` indices that group can read or write: the union of its
    # transforms' recorded read sets and their written indices. This is what
    # a group record snapshots instead of the whole parameter vector.
    group_relevant::Vector{Vector{Int}}
    # An iterated update's steps (`CTSEMRelinRecord`), and the row's prior
    # mean's cotangent, which its steps' reverses add to until the first of
    # them hands it on (`_relin_route!`).
    relins::Vector{CTSEMRelinRecord{T}}
    nrelins::Int
    xp_bar::Vector{T}
    # A moment-matched prediction's cubature (`CTSEMMomentRecord`).
    moments::Vector{CTSEMMomentRecord{T}}
    nmoments::Int
end

CTSEMAdjointTape(::Type{T}, group_relevant=[Int[], Int[], Int[]]) where {T} =
    CTSEMAdjointTape{T}(
        Tuple{Symbol,Int}[], CTSEMPredictRecord{T}[], CTSEMTDRecord{T}[],
        CTSEMUpdateRecord{T}[], CTSEMGroupRecord{T}[], CTSEMThetaRecord{T}[],
        CTSEMInitRecord{T}[], CTSEMBinaryRecord{T}[],
        CTSEMStationaryRecord{T}[], 0, 0, 0, 0, 0, 0, 0, 0,
        T[], CTSEMHybridTape{T}(), group_relevant, CTSEMRelinRecord{T}[], 0, T[],
        CTSEMMomentRecord{T}[], 0)

"""
Transition growth above which a subject's gradient is taken by forward mode
rather than by the reverse pass. The growth is that of the product of every
transition since the last observed row (`CTSEMGrowth`).

Not needed for accuracy since juliaFit 21b7a1a2. The reverse pass used to be
unstable in double precision once a prediction's covariance was very large in
some direction: the update's covariance cotangent came from an explicit `S^-1`
whose rounding (about 1e-6 on an innovation covariance of 1e10) was far larger
than its true value (about 1e-10), and the reverse through `e^{JAx dt}`
multiplied that direction by up to the square of the growth. Found on the SNSF
pilot: a person's drift with eigenvalues -1.05 and +0.28 over a 43-day gap,
growth 2e5, the adjoint gradient off by up to 1.7e4 (relative) against
BigFloat. The update's reverse is now stable (`_reverse_update!`), and on the
regression fixture the adjoint is closer to BigFloat than Float64 forward mode
is. This stays as an opt-in cross-check.

`Inf`, the reverse pass always, by default: forward mode needs a dual type the
fit has not otherwise used, so the first subject that takes it compiles the
filter again -- inside the Laplace route's seeded sweeps at a nested type. On a
one-latent model, 55 s on the plain route and 140 s on the Laplace one, once
per session (dev2); on the SNSF pilot's 32-effect model a single iteration
then ran 35 minutes on one core with the memory climbing past 25 GB. Opted
into per fit by optimcontrol's `explosive_forward`, at
`_CTSEM_EXPLOSIVE_GROWTH`. `ctsem_set_adjoint_growth!` moves it.
"""
const _CTSEM_ADJOINT_GROWTH = Ref(Inf)

"""
Transition growth past which a subject's filter pass counts as explosive: what
the progress line counts, what the fit's warning at the estimate and the
sampler's check of its draws report. Measured and reported whatever
`_CTSEM_ADJOINT_GROWTH` says, since the measurement costs nothing.
"""
const _CTSEM_EXPLOSIVE_GROWTH = Ref(100.0)

# Filter passes that crossed `_CTSEM_EXPLOSIVE_GROWTH`, since the session
# started; read as differences (`ctsem_explosive_passes`).
const _CTSEM_EXPLOSIVE_PASSES = Threads.Atomic{Int}(0)

"""Filter passes so far that crossed `_CTSEM_EXPLOSIVE_GROWTH`."""
ctsem_explosive_passes() = _CTSEM_EXPLOSIVE_PASSES[]

"""Set the transition growth above which subjects take the forward-mode gradient."""
function ctsem_set_adjoint_growth!(x::Real)
    x > 0 || throw(ArgumentError("threshold must be positive"))
    _CTSEM_ADJOINT_GROWTH[] = Float64(x)
    return _CTSEM_ADJOINT_GROWTH[]
end

"""
    _track_growth!(growth, E, n, subject)

Fold the transition `E[1:n, 1:n]` just applied into `growth` (`CTSEMGrowth`),
count the pass when it first exceeds `_CTSEM_EXPLOSIVE_GROWTH`, and note
`subject` then while `ctsem_explosive_subjects` is recording. A hand loop rather than
`_ctsem_mul!`: `E` may hold duals and only its primal values are wanted.
"""
function _track_growth!(growth::CTSEMGrowth, E::AbstractMatrix, n::Int, subject)
    phi, tmp = growth.phi, growth.tmp
    @inbounds if growth.steps == 0
        for j in 1:n, i in 1:n
            phi[i, j] = Float64(_primal(E[i, j]))
        end
    else
        for j in 1:n, i in 1:n
            acc = 0.0
            for k in 1:n
                acc += Float64(_primal(E[i, k])) * phi[k, j]
            end
            tmp[i, j] = acc
        end
        growth.phi, growth.tmp = tmp, phi
        phi = tmp
    end
    largest = 0.0
    @inbounds for j in 1:n
        column = 0.0
        for i in 1:n
            column += abs(phi[i, j])
        end
        largest = max(largest, column)
    end
    growth.steps += 1
    if largest > growth.max
        threshold = _CTSEM_EXPLOSIVE_GROWTH[]
        if growth.max <= threshold < largest
            Threads.atomic_add!(_CTSEM_EXPLOSIVE_PASSES, 1)
            _CTSEM_EXPLOSIVE_RECORDING[] && _note_explosive(subject)
        end
        growth.max = largest
    end
    return nothing
end

"""Whether the subject just filtered on `ws` needs the forward-mode gradient."""
@inline _growth_unstable(ws) = ws.growth.max > _CTSEM_ADJOINT_GROWTH[]

# Forward-mode subject gradients taken since the last reset, for the progress
# line (`_ctsem_forward_progress`) and the fit's result.
const _CTSEM_FORWARD_GRADIENTS = Threads.Atomic{Int}(0)

"""Forward-mode subject gradients taken since the last reset; `reset` zeroes the count."""
function ctsem_forward_gradients(; reset::Bool=false)
    return reset ? Threads.atomic_xchg!(_CTSEM_FORWARD_GRADIENTS, 0) :
        _CTSEM_FORWARD_GRADIENTS[]
end

const _CTSEM_EXPLOSIVE_RECORDING = Ref(false)
const _CTSEM_EXPLOSIVE_SUBJECTS = Set{Int}()
const _CTSEM_EXPLOSIVE_LOCK = ReentrantLock()

_note_explosive(subject) =
    lock(() -> push!(_CTSEM_EXPLOSIVE_SUBJECTS, Int(subject)), _CTSEM_EXPLOSIVE_LOCK)

"""
    ctsem_explosive_subjects(f)

Run `f()` and return its value with the sorted subject indices whose filter,
during it, exceeded `_CTSEM_EXPLOSIVE_GROWTH` -- the subjects whose predictions
grow that much between observations at the points `f` evaluates. Not
reentrant.
"""
function ctsem_explosive_subjects(f)
    lock(() -> empty!(_CTSEM_EXPLOSIVE_SUBJECTS), _CTSEM_EXPLOSIVE_LOCK)
    _CTSEM_EXPLOSIVE_RECORDING[] = true
    value = try
        f()
    finally
        _CTSEM_EXPLOSIVE_RECORDING[] = false
    end
    return value, lock(() -> sort!(collect(_CTSEM_EXPLOSIVE_SUBJECTS)),
        _CTSEM_EXPLOSIVE_LOCK)
end

"""
    _with_forward_count(f)

`f()`'s named tuple with `forward_gradients`, the forward-mode subject gradients
taken during it, and `explosive_passes`, the filter passes that crossed
`_CTSEM_EXPLOSIVE_GROWTH`, added: how a sampler reports them.
"""
function _with_forward_count(f)
    n0 = ctsem_forward_gradients()
    e0 = ctsem_explosive_passes()
    result = f()
    return merge(result, (forward_gradients=ctsem_forward_gradients() - n0,
        explosive_passes=ctsem_explosive_passes() - e0))
end

"""
    ctsem_explosive_draws(objective, draws)

Evaluate `objective` at each column of `draws` and report how many columns had
a subject whose filter exceeded `_CTSEM_EXPLOSIVE_GROWTH`, and which subjects
(`[0]` for none). A Laplace objective filters each subject at its effects' mode
given the column.
"""
function ctsem_explosive_draws(objective, draws::AbstractMatrix)
    hits = 0
    subjects = Set{Int}()
    for j in axes(draws, 2)
        _, s = ctsem_explosive_subjects(() ->
            ctsem_evaluate(objective, collect(Float64, view(draws, :, j)); gradient=false))
        isempty(s) && continue
        hits += 1
        union!(subjects, s)
    end
    return (checked=size(draws, 2), explosive=hits,
        subjects=isempty(subjects) ? [0] : sort!(collect(subjects)))
end

export ctsem_explosive_draws

"""
The progress line's counts since `forward0` and `explosive0`, each when any:
filter passes past `_CTSEM_EXPLOSIVE_GROWTH`, and forward-mode subject
gradients (optimcontrol's `explosive_forward`).
"""
function _ctsem_forward_progress(forward0::Int, explosive0::Int)
    e = ctsem_explosive_passes() - explosive0
    n = ctsem_forward_gradients() - forward0
    e > 0 || n > 0 || return ()
    n > 0 || return (@sprintf("explosive %d", e),)
    return (@sprintf("explosive %d", e), @sprintf("fwd-grad %d", n))
end

"""
    _ctsem_forward_subject_gradient!(g, subject_objective, values)

`g` set to the gradient of one subject's log likelihood at `values`, by forward
mode, on a workspace of its own (the subject's cached one can be shared by
sampler chains filtering the same subject). Behind a barrier, so the filter is
compiled for the new dual type only when a subject actually needs this. Works
at any element type, the nested duals of the Hessian and the Laplace sweeps
included: ForwardDiff's tag for this call is created after theirs.
"""
function _ctsem_forward_subject_gradient!(g::AbstractVector, so, values::AbstractVector,
    mode::Int=0, dfix::AbstractVector=Float64[])
    Threads.atomic_add!(_CTSEM_FORWARD_GRADIENTS, 1)
    _ctsem_barrier(_ctsem_forward_subject_gradient_impl!, g, so, values, mode, dfix)
    return g
end

# `mode` and `dfix`: what the pass returns for deviations the filter
# integrates, and the fixed ones it is handed; see `CTSEMHybridBuffers`.
function _ctsem_forward_subject_gradient_impl!(g, so, values, mode::Int=0,
    dfix::AbstractVector=Float64[])
    f = function (x)
        ws = _init_continuous_ekf_workspace(eltype(x), so.params)
        ws.hyb.mode = mode
        ws.hyb.dfix = convert(Vector{eltype(x)}, collect(dfix))
        return _ekf_run(ws, so, x, _ctsem_tipred_vector(so.tipreds, x), nothing)
    end
    # At a bucketed width (`_ctsem_gradient`), not ForwardDiff's own, which
    # compiled the filter again for each parameter count.
    copyto!(g, _ctsem_gradient(f, values))
    return nothing
end

"""
    _tape_reset!(tape)

Begin a new pass, keeping the records the last one built.

Emptying the record vectors is correct and was what this did; it is also why a
24-row subject of a one-latent model allocated 1,263 arrays in its forward
sweep alone, every one of them a 1x1 `Matrix{Float64}` that the previous
subject had already allocated at exactly that shape. The tape is per chunk and
a chunk's subjects share a model, so after the first subject every shape is
already there.

Resetting counts instead means the records are *reused*: `_record_predict!` and
friends copy into `tape.predicts[i]` when `i` is within reach and only build
when the tape has to grow. Nothing outside the recording hooks reads past the
count, and the reverse pass reads records by the index `program` recorded, so a
stale record beyond the count is unreachable rather than merely unread.
"""
function _tape_reset!(tape::CTSEMAdjointTape)
    empty!(tape.program)
    tape.npredicts = 0; tape.ntds = 0; tape.nupdates = 0
    tape.ngroups = 0; tape.nthetas = 0; tape.ninits = 0
    tape.nbinaries = 0; tape.nstationaries = 0; tape.nrelins = 0
    tape.nmoments = 0
    tape.hyb.active = false
    return tape
end

@inline _tape_push!(tape::CTSEMAdjointTape, kind::Symbol, index::Int) =
    push!(tape.program, (kind, index))

"""
    _tape_fill!(destination, source)

Overwrite a pooled record field, reallocating only when the shape changed.

A vector is `resize!`d, which keeps its capacity across resets and so is free
after the first pass. A matrix cannot be reshaped in place, so a genuine shape
change -- which for these records means a row with a different set of observed
manifests -- builds a new one and the caller stores it back. That is why the
records are mutable.
"""
@inline function _tape_fill!(destination::Vector{T}, source) where {T}
    resize!(destination, length(source))
    copyto!(destination, source)
    return destination
end

@inline function _tape_fill!(destination::Matrix{T}, source::AbstractMatrix) where {T}
    size(destination) == size(source) || return Matrix{T}(source)
    copyto!(destination, source)
    return destination
end

# Gathering copies, for the rows a record keeps.
#
# `view(pars.LAMBDA, o, 1:n)` with `o` a *vector* of indices builds a SubArray
# that is not strided, and unlike a range-indexed view it does not stay on the
# stack: line-level tracking put 96 KB per 300 subject evaluations on one of
# these, and there are four per row. The view existed only to avoid
# materialising `pars.LAMBDA[o, 1:n]`, which is the right instinct -- copying
# straight out of the parent with the index achieves it without the view.
@inline function _tape_gather!(destination::Vector{T}, source, rows) where {T}
    resize!(destination, length(rows))
    @inbounds for k in eachindex(rows)
        destination[k] = source[rows[k]]
    end
    return destination
end

@inline function _tape_gather!(destination::Matrix{T}, source, rows, ncol::Int) where {T}
    size(destination) == (length(rows), ncol) ||
        return T[source[rows[k], j] for k in eachindex(rows), j in 1:ncol]
    @inbounds for j in 1:ncol, k in eachindex(rows)
        destination[k, j] = source[rows[k], j]
    end
    return destination
end

@inline function _tape_gather_square!(destination::Matrix{T}, source, rows) where {T}
    m = length(rows)
    size(destination) == (m, m) ||
        return T[source[rows[i], rows[j]] for i in 1:m, j in 1:m]
    @inbounds for j in 1:m, i in 1:m
        destination[i, j] = source[rows[i], rows[j]]
    end
    return destination
end

################################################################################
# Recording hooks
################################################################################
#
# Every hook has a `::Nothing` method so the untraced primal pays nothing: with
# `trace === nothing` the calls have no arguments the compiler cannot see
# through and disappear entirely. This is why the forward filter has exactly
# one loop rather than a primal loop and a parallel traced copy -- the reverse
# pass cannot drift out of step with a forward pass it is literally recorded
# from.

@inline _record_subject_values!(::Nothing, args...) = nothing
@inline _record_init!(::Nothing, args...) = nothing
@inline _record_theta!(::Nothing, args...) = nothing
@inline _record_group!(::Nothing, args...) = nothing
@inline _record_td!(::Nothing, args...) = nothing
@inline _record_update!(::Nothing, args...) = nothing
# Every trace but the adjoint's ignores an iterated update's steps: what they
# record of a row is its final update.
@inline _record_relin!(trace, args...) = nothing
@inline _begin_predict!(::Nothing, args...) = nothing
@inline _record_predict!(::Nothing, args...) = nothing
@inline _record_stationary!(::Nothing, args...) = nothing
@inline _record_stationary_elements!(::Nothing, args...) = nothing

function _record_subject_values!(tape::CTSEMAdjointTape, subject_values)
    resize!(tape.subject_values, length(subject_values))
    copyto!(tape.subject_values, subject_values)
    return nothing
end

function _record_init!(tape::CTSEMAdjointTape{T}, pars, n::Int,
        popraw=nothing) where {T}
    index = (tape.ninits += 1)
    source = view(pars.T0VAR, 1:n, 1:n)
    if index <= length(tape.inits)
        record = tape.inits[index]
        record.T0VAR = _tape_fill!(record.T0VAR, source)
        record.POPRAW = _tape_fill_raw!(record.POPRAW, popraw)
    else
        push!(tape.inits, CTSEMInitRecord{T}(Matrix{T}(source),
            _tape_fill_raw!(zeros(T, 0, 0), popraw)))
    end
    _tape_push!(tape, :init, index)
    return nothing
end

# The raw population matrix, or `0 x 0` for a model without one. By method
# rather than by a local holding either, which would be a union of two views
# and boxed.
@inline _tape_fill_raw!(destination::Matrix{T}, ::Nothing) where {T} =
    size(destination) == (0, 0) ? destination : zeros(T, 0, 0)
@inline _tape_fill_raw!(destination::Matrix, popraw::AbstractMatrix) =
    _tape_fill!(destination, popraw)

"""
    _record_stationary!(tape, ws, pars)

Record `_ctsem_stationary!` just after it ran: the mean is in `ws.state` and the
covariance still in `ws.diffusion_buffer.out`, where it left them.
"""
function _record_stationary!(tape::CTSEMAdjointTape{T}, ws, pars) where {T}
    index = (tape.nstationaries += 1)
    n = ws.state_dim
    naff = ws.affine_buffer.dim
    k = length(ws.diffusion_state_indices)
    drift = view(pars.DRIFT, 1:naff, 1:naff)
    mean = view(ws.state, 1:naff)
    X = view(ws.diffusion_buffer.out, 1:k, 1:k)
    diffusion = view(pars.DIFFUSION, 1:n, 1:n)
    if index <= length(tape.stationaries)
        record = tape.stationaries[index]
        record.DRIFT = _tape_fill!(record.DRIFT, drift)
        record.mean = _tape_fill!(record.mean, mean)
        record.X = _tape_fill!(record.X, X)
        record.DIFFUSION = _tape_fill!(record.DIFFUSION, diffusion)
    else
        push!(tape.stationaries, CTSEMStationaryRecord{T}(Matrix{T}(drift),
            Vector{T}(mean), Matrix{T}(X), Matrix{T}(diffusion), zeros(T, 0, 0),
            zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0)))
    end
    _tape_push!(tape, :stationary, index)
    return nothing
end

"""
    _record_stationary_elements!(tape, ws, pars)

Record the single stationary elements just after they were placed: the mean
in `ws.affine_buffer.s`, the covariance in `ws.diffusion_buffer.out`, and the
initial covariance before them, the factors and the congruence in
`ws.discretization_buffer`, where `_ctsem_stationary_elements_cov!` left them.
What a model without means or blocks never computed is recorded and never read.
"""
function _record_stationary_elements!(tape::CTSEMAdjointTape{T}, ws, pars) where {T}
    index = (tape.nstationaries += 1)
    n = ws.state_dim
    naff = ws.affine_buffer.dim
    k = length(ws.diffusion_state_indices)
    buf = ws.discretization_buffer
    drift = view(pars.DRIFT, 1:naff, 1:naff)
    mean = view(ws.affine_buffer.s, 1:naff)
    X = view(ws.diffusion_buffer.out, 1:k, 1:k)
    diffusion = view(pars.DIFFUSION, 1:n, 1:n)
    pold = view(buf.input, 1:n, (n + 1):(2n))
    lx = view(buf.output, 1:n, 1:n)
    lp = view(buf.output, 1:n, (n + 1):(2n))
    m = view(buf.scratch, 1:n, 1:n)
    if index <= length(tape.stationaries)
        record = tape.stationaries[index]
        record.DRIFT = _tape_fill!(record.DRIFT, drift)
        record.mean = _tape_fill!(record.mean, mean)
        record.X = _tape_fill!(record.X, X)
        record.DIFFUSION = _tape_fill!(record.DIFFUSION, diffusion)
        record.Pold = _tape_fill!(record.Pold, pold)
        record.LX = _tape_fill!(record.LX, lx)
        record.LP = _tape_fill!(record.LP, lp)
        record.M = _tape_fill!(record.M, m)
    else
        push!(tape.stationaries, CTSEMStationaryRecord{T}(Matrix{T}(drift),
            Vector{T}(mean), Matrix{T}(X), Matrix{T}(diffusion), Matrix{T}(pold),
            Matrix{T}(lx), Matrix{T}(lp), Matrix{T}(m)))
    end
    _tape_push!(tape, :stationary_elements, index)
    return nothing
end

function _record_theta!(tape::CTSEMAdjointTape{T}, pars, m::Int) where {T}
    index = (tape.nthetas += 1)
    source = view(pars.MANIFESTVAR, 1:m, 1:m)
    if index <= length(tape.thetas)
        record = tape.thetas[index]
        record.MANIFESTVAR = _tape_fill!(record.MANIFESTVAR, source)
    else
        push!(tape.thetas, CTSEMThetaRecord{T}(Matrix{T}(source)))
    end
    _tape_push!(tape, :theta, index)
    return nothing
end

"""
Record a transform group. **Must be called before the group is applied** --
see `CTSEMGroupRecord` for why the pre-group parameter values are the ones
the reverse pass needs.
"""
function _record_group!(tape::CTSEMAdjointTape{T}, group::Int,
    indices::AbstractVector{Int}, all_params, ctx::CTSEMRowContext) where {T}
    index = _record_group_slot!(tape, group, indices, all_params, ctx)
    index == 0 || _tape_push!(tape, :group, index)
    return nothing
end

"""
Record a transform group into the tape's next group slot without a program
entry, and return the slot (0 for an empty group): `_record_group!`'s record,
and a cubature point's (`_record_moment_group!`), which its moment record's
reverse walks instead.
"""
function _record_group_slot!(tape::CTSEMAdjointTape{T}, group::Int,
    indices::AbstractVector{Int}, all_params, ctx::CTSEMRowContext) where {T}
    # Linear models have no state-dependent cells at all; skipping the
    # snapshot there keeps the tape out of the common case entirely.
    isempty(indices) && return 0
    relevant = tape.group_relevant[group]
    index = (tape.ngroups += 1)
    # Explicit `{T}` for the same reason as in `_record_update!`: the TD
    # predictors and the time/interval are data, so they arrive as Float64
    # whatever the tape's element type is, and the record requires every field
    # to share it.
    if index <= length(tape.groups)
        record = tape.groups[index]
        record.group = group
        # `params_before` is as long as the group's relevant set, and the three
        # groups have different sets -- so this is the one pooled field whose
        # length genuinely varies between neighbouring records. `resize!` keeps
        # the capacity, so it still costs nothing after the first pass.
        resize!(record.params_before, length(relevant))
        @inbounds for j in eachindex(relevant)
            record.params_before[j] = all_params[relevant[j]]
        end
        _tape_fill!(record.state, vec(ctx.state))
        _tape_fill!(record.tdpreds, ctx.tdpreds)
        record.time = T(ctx.time)
        record.dt = T(ctx.dt)
        record.row = ctx.row
    else
        compact = Vector{T}(undef, length(relevant))
        @inbounds for j in eachindex(relevant)
            compact[j] = all_params[relevant[j]]
        end
        push!(tape.groups, CTSEMGroupRecord{T}(group, compact, collect(vec(ctx.state)),
            T[x for x in ctx.tdpreds], T(ctx.time), T(ctx.dt), ctx.row))
    end
    return index
end

# A cubature point's group evaluation (`_ekf_moment_cells!`): on the tape, a
# slot its moment record's reverse walks; nothing on any other trace.
@inline _record_moment_group!(::Any, args...) = 0
@inline _record_moment_group!(tape::CTSEMAdjointTape, group::Int,
    indices::AbstractVector{Int}, all_params, ctx::CTSEMRowContext) =
    _record_group_slot!(tape, group, indices, all_params, ctx)

@inline _record_moment!(::Any, args...) = nothing

"""Record a moment-matched substep's cubature (`CTSEMMomentRecord`)."""
function _record_moment!(tape::CTSEMAdjointTape{T}, ws, first::Int, npoints::Int,
    n::Int) where {T}
    index = (tape.nmoments += 1)
    if index > length(tape.moments)
        push!(tape.moments, CTSEMMomentRecord{T}(0, 0, zeros(T, 0, 0), zeros(T, 0, 0),
            zeros(T, 0, 0)))
    end
    record = tape.moments[index]
    record.first = first
    record.npoints = npoints
    record.U = _tape_fill!(record.U, view(ws.moment_U, 1:n, 1:n))
    record.D = _tape_fill!(record.D, view(ws.moment_d, 1:n, 1:npoints))
    record.DRIFTs = _tape_fill!(record.DRIFTs, view(ws.moment_drifts, 1:n, 1:(n * npoints)))
    _tape_push!(tape, :moment, index)
    return nothing
end

function _record_td!(tape::CTSEMAdjointTape{T}, ws, pars, tdpreds, n::Int) where {T}
    # Mirrors `_apply_td_impulse!`'s own early return: with no TD predictors
    # the impulse is a no-op and contributes nothing to reverse.
    isempty(tdpreds) && return nothing
    index = (tape.ntds += 1)
    P_in = view(ws.P_predict.data, 1:n, 1:n)
    Jtd = view(pars.Jtd, 1:n, 1:n)
    if index <= length(tape.tds)
        record = tape.tds[index]
        record.P_in = _tape_fill!(record.P_in, P_in)
        record.Jtd = _tape_fill!(record.Jtd, Jtd)
        _tape_fill!(record.tdpreds, tdpreds)
    else
        push!(tape.tds, CTSEMTDRecord{T}(Matrix{T}(P_in), Matrix{T}(Jtd),
            T[x for x in tdpreds]))
    end
    _tape_push!(tape, :td, index)
    return nothing
end

function _record_update!(tape::CTSEMAdjointTape{T}, ws, pars, data, obs_col::Int,
    observed::AbstractVector{Int}, state_in, P_in, n::Int, relin::Int=0) where {T}
    index = (tape.nupdates += 1)
    # `T[...]` rather than `[float(...)]`: the observed data is a constant with
    # respect to the parameters, so it is Float64 whatever the tape's element
    # type is. Every other field is `T`, and `CTSEMUpdateRecord{T}` requires
    # them all to agree -- which they did as long as `T` was `Float64` too.
    # Differentiating this gradient (see `ctsem_hessian`) makes `T` a dual, and
    # the record then has to carry the data as a dual with zero partials.
    if index > length(tape.updates)
        o = collect(observed)
        push!(tape.updates, CTSEMUpdateRecord{T}(
            o, collect(vec(state_in)), Matrix(P_in),
            Matrix(pars.LAMBDA[o, 1:n]), collect(pars.MANIFESTMEANS[o]),
            Matrix(pars.Jy[o, 1:n]), Matrix(ws.bufferΘ.out[o, o]),
            T[data[i, obs_col] for i in o], zeros(T, 0, 0), 0, T[], zeros(T, 0, 0)))
        _hybrid_record_X!(tape.updates[index], ws.hyb, n)
        _record_update_relin!(tape.updates[index], ws, pars, relin, n)
        _tape_push!(tape, :update, index)
        return nothing
    end
    record = tape.updates[index]
    o = _tape_fill!(record.observed, observed)
    # Views rather than slices throughout: `pars.LAMBDA[o, 1:n]` materialises a
    # matrix only to copy it and throw it away, and this runs once per row.
    _tape_fill!(record.state_in, state_in)
    record.P_in = _tape_fill!(record.P_in, P_in)
    record.Lambda = _tape_gather!(record.Lambda, pars.LAMBDA, o, n)
    _tape_gather!(record.manifestmeans, pars.MANIFESTMEANS, o)
    record.H = _tape_gather!(record.H, pars.Jy, o, n)
    record.R = _tape_gather_square!(record.R, ws.bufferΘ.out, o)
    resize!(record.y, length(o))
    @inbounds for j in eachindex(o)
        record.y[j] = T(data[o[j], obs_col])
    end
    _hybrid_record_X!(record, ws.hyb, n)
    _record_update_relin!(record, ws, pars, relin, n)
    _tape_push!(tape, :update, index)
    return nothing
end

# A relinearised update (`_ekf_iterated_rows!`) loads on Jy with the intercept
# at its linearisation point, `ws.state` when it is recorded; the record keeps
# LAMBDA there and the point, for the reverse.
function _record_update_relin!(record::CTSEMUpdateRecord{T}, ws, pars, relin::Int,
    n::Int) where {T}
    record.relin = relin
    relin == 0 && return nothing
    o = record.observed
    record.Lambda_lin = _tape_gather!(record.Lambda_lin, pars.LAMBDA, o, n)
    record.Lambda = _tape_gather!(record.Lambda, pars.Jy, o, n)
    _tape_fill!(record.manifestmeans, view(ws.mm_eff, 1:length(o)))
    _tape_fill!(record.x_lin, vec(ws.state))
    return nothing
end

function _record_relin!(tape::CTSEMAdjointTape{T}, ws, pars, data, obs_col::Int,
    observed::AbstractVector{Int}, first::Bool, alpha::Float64) where {T}
    n = ws.state_dim
    m = length(observed)
    index = (tape.nrelins += 1)
    index > length(tape.relins) && push!(tape.relins, CTSEMRelinRecord{T}(Int[], T[],
        zeros(T, 0, 0), zeros(T, 0, 0), T[], zeros(T, 0, 0), T[], zeros(T, 0, 0), T[],
        first, alpha))
    record = tape.relins[index]
    o = _tape_fill!(record.observed, observed)
    _tape_fill!(record.state_in, ws.xprior)
    record.P_in = _tape_fill!(record.P_in, view(ws.P_predict.data, 1:n, 1:n))
    record.H = _tape_gather!(record.H, pars.Jy, o, n)
    _tape_fill!(record.offset, view(ws.mm_eff, 1:m))
    record.Lambda = _tape_gather!(record.Lambda, pars.LAMBDA, o, n)
    _tape_fill!(record.x_lin, vec(ws.state))
    record.R = _tape_gather_square!(record.R, ws.bufferΘ.out, o)
    resize!(record.y, m)
    @inbounds for j in 1:m
        record.y[j] = T(data[o[j], obs_col])
    end
    record.first = first
    record.alpha = alpha
    _tape_push!(tape, :relin, index)
    return nothing
end

"""An all-zero predict record shaped for `n` states and `k` dynamic states."""
# `k` sizes the Lyapunov block (the diffusing states), `naff` the affine one
# (the leading genuine dynamics). They are different sets; see `affine_dim`.
_empty_predict_record(::Type{T}, n::Int, k::Int, naff::Int=k) where {T} =
    CTSEMPredictRecord{T}(zeros(T, n), zeros(T, n, n), zeros(T, n, n),
        zeros(T, n, n), zeros(T, n, n), zeros(T, n, n), zeros(T, k, k),
        zeros(T, naff), zeros(T, naff), zero(T), false, false, zeros(T, 0, 0),
        CTSEMSeriesChain(T), CTSEMSeriesChain(T), zeros(T, 0, 0), zeros(T, 0, 0),
        false)

"""
Snapshot the substep inputs, which `_ekf_predict_step!` overwrites.

Returns the pooled record the matching `_record_predict!` will finish, rather
than a fresh tuple -- the two are called in lockstep around one
`_ekf_predict_step!`, so the slot is known here and the snapshot can go straight
into it. Growing the tape here rather than there keeps the return type a plain
`CTSEMPredictRecord`, so the filter's local stays concrete.
"""
function _begin_predict!(tape::CTSEMAdjointTape{T}, ws, n::Int) where {T}
    index = tape.npredicts + 1
    if index > length(tape.predicts)
        push!(tape.predicts,
            _empty_predict_record(T, n, length(ws.diffusion_state_indices),
                ws.affine_buffer.dim))
    end
    record = tape.predicts[index]
    _tape_fill!(record.state_in, ws.state)
    record.P_in = _tape_fill!(record.P_in, view(ws.P_update.data, 1:n, 1:n))
    _hybrid_record_X!(record, ws.hyb, n)
    return record
end

function _record_predict!(tape::CTSEMAdjointTape{T}, ws, pars,
    record::CTSEMPredictRecord{T}, Δt, n::Int) where {T}
    dyn = ws.diffusion_state_indices
    k = length(dyn)
    index = (tape.npredicts += 1)
    record.A = _tape_fill!(record.A, view(ws.discrete_ca.eJAx, 1:n, 1:n))
    record.JAx = _tape_fill!(record.JAx, view(pars.JAx, 1:n, 1:n))
    record.DRIFT = _tape_fill!(record.DRIFT, view(pars.DRIFT, 1:n, 1:n))
    record.DIFFUSION = _tape_fill!(record.DIFFUSION, view(pars.DIFFUSION, 1:n, 1:n))
    record.Xlyap = _tape_fill!(record.Xlyap, view(ws.diffusion_buffer.out, 1:k, 1:k))
    # The affine offset and the solved intercept live on the affine block, not
    # the diffusion one, so they come from `affine_buffer` and from the leading
    # entries of `dINT`. (Unused by the discrete-time reverse pass, which has
    # neither a solve nor a Lyapunov term.)
    naff = ws.affine_buffer.dim
    _tape_fill!(record.affine, view(ws.affine_buffer.r, 1:naff))
    _tape_fill!(record.dINT_dynamic, view(ws.discrete_ca.dINT, 1:naff))
    record.dt = T(Δt)
    series = ws.discretization_buffer.series
    icpt = ws.discretization_buffer.series_icpt
    record.series_intercept = icpt.intercept
    record.series_noise = series.noise
    if series.noise
        record.Qd = _tape_fill!(record.Qd, view(series.Q, 1:k, 1:k))
        _series_chain_store!(record.noise_chain, series, k)
    end
    icpt.intercept && _series_chain_store!(record.icpt_chain, icpt, naff)
    _hybrid_record_Xinc!(record, ws.hyb, n)
    record.moment = ws.moment_now[]
    _tape_push!(tape, :predict, index)
    return nothing
end

################################################################################
# Reverse pass
################################################################################

"""
    _flush_frechet!(aws)

Push any deferred matrix-exponential Fréchet derivative into the JAx cotangent.

`A = exp(JAx * dt)` is the single most expensive step in the reverse pass: its
adjoint is a matrix exponential of *twice* the state dimension. Profiling the
equivalent C++ port put it at ~66% of the whole reverse pass on a 20-latent
model -- ~540 µs per prediction substep against ~30 µs for everything else in
that substep combined.

But the Fréchet derivative is **linear in its direction argument**, and a
balanced panel design hands it the same `JAx * dt` at every substep of every
subject. So instead of `Σ_steps dt * L(A', Ā_step)`, accumulate the directions
and evaluate `dt * L(A', Σ_steps Ā_step)` once. That is exact, not an
approximation, and it turns an O(rows) count of block exponentials into
O(distinct `(JAx, dt)` pairs) -- one, for a linear model on an equally spaced
panel.

The key, in `_frechet_enqueue!`, is an exact comparison against the actual
`JAx * dt` that produced each pending direction, so nothing here assumes
linearity: a state-dependent model changes `JAx` per row, matches nothing and
pays only the comparison, while a schedule with few distinct intervals
batches to one evaluation per distinct interval however its rows are ordered.

Correctness depends on nothing consuming the JAx cotangent while a flush is
outstanding. Two things can: a state-dependent transform group that *writes* a
JAx cell (whose reverse zeroes that cell's cotangent), and the parameter layer.
Both flush first.
"""
function _flush_frechet!(aws)
    count = aws.frechet_count
    count == 0 && return nothing
    aws.frechet_count = 0
    n = aws.n
    L = aws.frechet_L
    @inbounds for e in 1:count
        # `L(A', Ā)` by the Al-Mohy-Higham recurrence, into the workspace's
        # own buffers; the `_ctsem_exp_frechet_adjoint` wrapper does the same
        # computation but allocates, which at one flush per substep is what
        # the GC share of the reverse pass was made of.
        _CTSEM_OPCOUNT.frechet[] += 1
        my_exp_frechet!(aws.frechet_Y, L, transpose(aws.frechet_As[e]),
            aws.frechet_accums[e], aws.frechet_buffer)
        dt = aws.frechet_dts[e]
        if aws.defer_frechet
            _frechet_add!(aws.jax_bar_deferred, L, dt, n)
        else
            _frechet_add!(CTSEMPars(aws.theta_bar, aws.sp.parameter_axis).JAx, L, dt, n)
        end
    end
    return nothing
end

"""`target .+= dt .* L` over the leading `n × n` block, without a temporary."""
@inline function _frechet_add!(target, L, dt, n::Int)
    @inbounds for j in 1:n, i in 1:n
        target[i, j] += dt * L[i, j]
    end
    return nothing
end

"""
Queue the direction `Ā` for the exponential of `scaled = JAx * dt`.

A direction whose key is already in the table is added to that entry; a new
key takes the next free slot, and a full table is flushed first. `dt` is
compared before the matrix, so a schedule of irregular intervals pays one
scalar comparison per live entry and no matrix comparison at all.
"""
function _frechet_enqueue!(aws, dt, scaled::AbstractMatrix, Ā::AbstractMatrix)
    @inbounds for e in 1:aws.frechet_count
        if aws.frechet_dts[e] == dt && aws.frechet_As[e] == scaled
            aws.frechet_accums[e] .+= Ā
            return nothing
        end
    end
    aws.frechet_count == length(aws.frechet_dts) && _flush_frechet!(aws)
    e = aws.frechet_count + 1
    aws.frechet_dts[e] = dt
    copyto!(aws.frechet_As[e], scaled)
    copyto!(aws.frechet_accums[e], Ā)
    aws.frechet_count = e
    return nothing
end

"""
    _reverse_predict!(x̄, P̄, θ̄ca, record, dyn, n, lyap_buffer, aws)

Undo one prediction substep.

Forward, on the dynamic index subset `D = dyn` of size `k`:

    A            = exp(JAx * dt)
    affine[i]    = CINT[i] + Σⱼ (DRIFT[i,j] - JAx[i,j]) x[j]      (i ≤ naff)
    s[i]         = -affine[i] + Σ_q A[i,q] affine[q]              (q ≤ naff)
    dINT[1:naff] = JAx[1:naff,1:naff] \\ s
    X            = lyap(JAx[D,D], Qc[D,D])
    dDIFF[D,D]   = X - A[D,D] X A[D,D]'
    x⁺           = A x + dINT
    P⁺           = A (P + εI) A' + dDIFF

Note `A` is the matrix exponential itself: this port sets `dDRIFT = eJAx`, so
the discrete drift and the exponential are the same object and share one
cotangent.

Where a pivot was small against the interval the forward pass took either of
two lines as the integral it is instead -- `dINT[1:naff] = Φ affine` with
`Φ = ∫₀^dt e^{JAx s} ds`, and `dDIFF[D,D] = ∫₀^dt e^{JAx s} Qc e^{JAx' s} ds`
-- and the record says which (`series_intercept`, `series_noise`). Neither
reads `A`, so on those lines nothing reaches `Ā`; see
`series_discretization.jl`.

Two different sub-blocks appear, and conflating them was a bug. `D = dyn` is
`diffusion_state_indices`, the states with their own diffusion; the Lyapunov
solve and the `dDIFF` term live there. `1:naff` is the leading block of genuine
dynamics -- everything that is not a static random-effect carrier -- and the
intercept solve and the affine offset live there. `dyn` is a subset of
`1:naff`, equal to it for any model with a full DIFFUSION, and strictly smaller
when a latent has no diffusion of its own and no coupling to one that has some.
Running the affine offset over `dyn` dropped such a state's CINT from the
likelihood entirely. `A` spans the whole augmented state throughout.
"""
function _reverse_predict!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMPredictRecord{T}, dyn::AbstractVector{Int}, n::Int,
    lyap_buffer, aws, hyb::Union{Nothing,CTSEMHybridTape{T}}=nothing) where {T}
    A = record.A
    JAx = record.JAx
    x = record.state_in
    k = length(dyn)

    # The discrete-time reverse pass is the continuous one with its three hard
    # pieces removed rather than a second implementation of it: A is JAx (no
    # Frechet derivative of an exponential), dDIFFUSION is the diffusion
    # covariance itself (no Lyapunov pullback), and dINT is the affine offset
    # (no linear solve). Written out separately rather than branched inside the
    # shared recursion, which would put a test in every step of it.
    if !aws.sp.continuous_time
        return _reverse_predict_discrete!(x̄, P̄, θ̄ca, record, dyn, n, aws, hyb)
    end

    # Every temporary is a view of a scratch buffer; see `CTSEMReverseScratch`.
    # The derivation below is unchanged -- what changed is only where the memory
    # comes from. Buffers are shared with `_reverse_update!`, which is safe
    # because the two are called at different points of the tape walk and
    # neither is ever live while the other runs.
    sc = aws.reverse_scratch
    Ā          = _rs(sc.Abar, n, n)
    JAx_bar    = _rs(sc.nn3, n, n)
    P̄_new      = _rs(sc.Pbar_new, n, n)
    Ps         = _rs(sc.Ps, n, n)
    Ptilde     = _rs(sc.Pr, n, n)
    nn1        = _rs(sc.nn1, n, n)
    nn2        = _rs(sc.nn2, n, n)
    scaled     = _rs(sc.nn4, n, n)
    Qc_bar     = _rs(sc.nn5, n, n)
    Ad         = _rs(sc.kk1, k, k)
    JAxd       = _rs(sc.kk2, k, k)
    Qb         = _rs(sc.kk3, k, k)
    X̄          = _rs(sc.kk4, k, k)
    Ād         = _rs(sc.kk5, k, k)
    kk6        = _rs(sc.kk6, k, k)
    x̄_new      = _rs(sc.xbar_new, n)
    dINT_bar   = _rs(sc.nv1, n)
    # The intercept solve and the affine offset run over the leading `naff`
    # states, a different and larger block than `dyn`; both are live at once,
    # so these cannot share the `kk`/`kv` slots.
    naff       = aws.affine_dim
    s̄          = _rs(sc.av1, naff)
    affine_bar = _rs(sc.av2, naff)
    aa1        = _rs(sc.aa1, naff, naff)

    @inbounds for j in 1:k, i in 1:k
        Ad[i, j] = A[dyn[i], dyn[j]]
        JAxd[i, j] = JAx[dyn[i], dyn[j]]
    end
    fill!(Ā, zero(T))
    fill!(JAx_bar, zero(T))

    # --- mean: x⁺ = A x + dINT
    _ctsem_outer!(Ā, x̄, x, one(T), one(T))
    # and the integrated deviations' X⁺ = A X + Xinc (`CTSEMHybrid`), a
    # continuous intercept's increment reversed through `JAx_bar` and `Ā` as
    # the intercept's own is below.
    hyb === nothing || _hybrid_reverse_predict!(Ā, JAx_bar, A, JAx, hyb,
        aws.sp.hybrid, record, n, aws.affine_dim, sc, record.dt)
    copyto!(dINT_bar, x̄)
    _ctsem_mulTvec!(x̄_new, A, x̄)

    # --- covariance: P⁺ = A (P + εI) A' + dDIFF
    _symmetrize_into!(Ps, P̄)
    copyto!(Ptilde, record.P_in)
    _ridge_diagonal!(Ptilde, n, _CTSEM_RIDGE)
    # d(A P̃ A')/dA contracted with P̄ is P̄ A P̃' + P̄' A P̃; both P̄ (symmetrised
    # just above) and P̃ are symmetric, so that collapses to twice one term.
    _ctsem_mul!(nn1, Ps, A)
    _ctsem_mul!(Ā, nn1, Ptilde, T(2), one(T))
    _ctsem_mulTN!(nn2, A, Ps)
    _ctsem_mul!(P̄_new, nn2, A)
    # `dDIFF_bar` is `Ps`, so `Ps` has to survive until `Qb` is taken from it.

    # --- dDIFF[D,D] = X - Ad X Ad'
    @inbounds for j in 1:k, i in 1:k
        Qb[i, j] = (Ps[dyn[i], dyn[j]] + Ps[dyn[j], dyn[i]]) / 2
    end
    JAxd_bar = sc.kk7
    Qcd_bar = sc.kk8
    if record.series_noise
        # --- dDIFF[D,D] = V(JAx[D,D], Qc[D,D], dt), the integral itself; see
        # `series_discretization.jl`. It reads no exponential, so Ād is zero.
        # Its cotangents wait on the table's powers until `_flush_series!`, so
        # nothing reaches JAx or DIFFUSION from this row directly.
        fill!(Ād, zero(T))
        fill!(JAxd_bar, zero(T))
        fill!(Qcd_bar, zero(T))
        sr = aws.series_reverse
        if !(_series_table_matches(sr.noise, JAxd, record.Qd, k) &&
                _ctsem_same_matrix(sr.diffusion, record.DIFFUSION, n))
            _flush_series_noise!(aws, θ̄ca, n)
            _series_table_reset!(sr.noise, JAxd, record.Qd, k)
            _copy_block!(sr.diffusion, record.DIFFUSION, n, n)
        end
        # The chain the forward pass computed, and the powers its pullback
        # adds onto, which the flush pulls through.
        _series_chain_load!(sc.series, record.noise_chain, k)
        _series_table_extend!(sr.noise, sc.series.m)
        _series_noise_table_pullback!(sr.noise_bar, sc.series, sr.noise, k, Qb)
    else
        X = record.Xlyap
        _ctsem_mulTN!(kk6, Ad, Qb)
        _ctsem_mul!(X̄, kk6, Ad)
        X̄ .= Qb .- X̄                                   # X̄ = Qb - Ad' Qb Ad
        _ctsem_mul!(kk6, Qb, Ad)
        _ctsem_mulNT!(Ād, kk6, X)
        _ctsem_mulTN!(kk6, Qb, Ad)
        _ctsem_mul!(Ād, kk6, X, one(T), one(T))
        Ād .= .-Ād                                     # Ād = -(Qb Ad X' + Qb' Ad X)

        # --- X = lyap(JAx[D,D], Qc[D,D])
        _ctsem_lyap_pullback!(JAxd_bar, Qcd_bar, sc.kk9, sc.kk10, JAxd, X, X̄, lyap_buffer)
    end

    if record.series_intercept
        # --- dINT[1:naff] = Phi affine, Phi = int_0^dt e^{JAx s} ds over the
        # leading block, taken as a series; see `series_discretization.jl`.
        # Like the noise series it reads no exponential, so nothing reaches Ā,
        # and its cotangent waits on the table's powers (`_flush_series!`).
        sr = aws.series_reverse
        if !_series_table_matches(sr.intercept, JAx, nothing, naff)
            _flush_series_intercept!(aws, θ̄ca)
            _series_table_reset!(sr.intercept, JAx, nothing, naff)
        end
        _series_chain_load!(sc.series, record.icpt_chain, naff)
        _series_table_extend!(sr.intercept, sc.series.s > 0 ? sc.series.m + 1 : sc.series.m)
        Phi = _series_block(sc.series.out, naff)
        @inbounds for i in 1:naff
            acc = zero(T)
            for q in 1:naff
                acc += Phi[q, i] * dINT_bar[q]
            end
            affine_bar[i] = acc
        end
        @inbounds for j in 1:naff, i in 1:naff
            aa1[i, j] = dINT_bar[i] * record.affine[j]
        end
        _series_intercept_table_pullback!(sr.intercept_bar, sc.series, sr.intercept,
            naff, aa1)
    else
        # --- dINT[1:naff] = JAx[1:naff,1:naff] \ s   (s̄ = JAx⁻ᵀ dINT_bar, M̄ = -s̄ dINT')
        # The affine block is a leading one, so everything below indexes directly
        # with no gather, and its two outer products go straight into the full
        # cotangents rather than through a sub-block that then has to be scattered.
        @inbounds for i in 1:naff; s̄[i] = dINT_bar[i]; end
        # `transpose(JAx) \ s̄` would be a LAPACK `getrf!`/`getrs!` pair, once per
        # prediction substep, on a matrix of the dynamic-state size. At that size
        # the call is mostly OpenBLAS's process-global buffer lock -- see
        # `small_linalg.jl` -- so it goes through the engine's own LU instead.
        @inbounds for j in 1:naff, i in 1:naff; aa1[i, j] = JAx[j, i]; end
        _solve_square_system_generic!(aa1, s̄, sc.piv, naff)
        @inbounds for j in 1:naff, i in 1:naff
            JAx_bar[i, j] -= s̄[i] * record.dINT_dynamic[j]
        end

        # --- s = -affine + A affine, on the same leading block
        @inbounds for i in 1:naff
            acc = zero(T)
            for q in 1:naff
                acc += A[q, i] * s̄[q]
            end
            affine_bar[i] = acc - s̄[i]
        end
        @inbounds for j in 1:naff, i in 1:naff
            Ā[i, j] += s̄[i] * record.affine[j]
        end
    end

    # --- affine[i] = CINT[i] + Σⱼ (DRIFT[i,j] - JAx[i,j]) x[j]
    # The cotangent's views bound once, outside the loops that index them.
    θ̄CINT = θ̄ca.CINT
    θ̄DRIFT = θ̄ca.DRIFT
    θ̄JAx = θ̄ca.JAx
    @inbounds for i in 1:naff
        θ̄CINT[i] += affine_bar[i]
    end
    @inbounds for j in 1:n, i in 1:naff
        contribution = affine_bar[i] * x[j]
        θ̄DRIFT[i, j] += contribution
        JAx_bar[i, j] -= contribution
        x̄_new[j] += (record.DRIFT[i, j] - JAx[i, j]) * affine_bar[i]
    end
    # The affine offset also carried a moment-matched substep's mean
    # correction, whose cotangent is the offset's (`_reverse_moment!`).
    record.moment && _moment_cbar!(aws, affine_bar, naff)

    # Scatter the dynamic-block contributions back into the full matrices.
    @inbounds for j in 1:k, i in 1:k
        Ā[dyn[i], dyn[j]] += Ād[i, j]
        JAx_bar[dyn[i], dyn[j]] += JAxd_bar[i, j]
    end

    @inbounds for j in 1:n, i in 1:n
        θ̄JAx[i, j] += JAx_bar[i, j]
    end

    # --- A = exp(JAx * dt). The direction is *queued* rather than pushed
    # through immediately, so substeps sharing a `JAx * dt` cost one block
    # exponential between them instead of one each -- see `_flush_frechet!`.
    # The exponential itself still goes through `Base.exp` rather than the
    # package's `my_exp!`, deliberately: buffering it and routing through
    # `my_exp!` was tried and made the adjoint *slower*, up to 2x on small
    # models, because `my_exp!` always runs the degree-13 Padé approximant
    # while `Base.exp` picks a lower degree adaptively from the matrix norm,
    # and `JAx * dt` is typically small-norm.
    scaled .= record.dt .* JAx
    # `==` on two matrices is an elementwise comparison and allocates nothing;
    # a `Val(n)`-dispatched helper would, because `n` is a runtime value here
    # and constructing the `Val` costs a dynamic dispatch per row.
    _frechet_enqueue!(aws, record.dt, scaled, Ā)

    # --- Qc = sdcovsqrt2cov(DIFFUSION); only the dynamic block was consumed.
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Qcd_bar[i, j]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, θ̄ca, n)

    copyto!(x̄, x̄_new)
    copyto!(P̄, P̄_new)
    return nothing
end

"""
    _reverse_stationary!(x̄, P̄, θ̄ca, record, dyn, n, aws)

Undo `_ctsem_stationary!`. The forward pass overwrote the mean over the leading
`naff` states and every entry of those rows and columns of the covariance, so
their cotangents end here; the `:init` entry before it sees only what the
overwrite left, which in those rows is nothing.

  * mean, `m = -A \\ c` with `A = DRIFT[1:naff, 1:naff]`: `w = A' \\ m̄`, then
    `Ā -= w m'` and `c̄ -= w`.
  * covariance, `A_D X + X A_D' + Q_D = 0` over the diffusion states: the
    prediction's own Lyapunov pullback, and `Q̄` deferred through DIFFUSION's
    construction exactly as a prediction's is, keyed on the same raw matrix.
"""
function _reverse_stationary!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMStationaryRecord{T}, dyn::AbstractVector{Int}, n::Int,
    aws) where {T}
    sc = aws.reverse_scratch
    naff = aws.affine_dim
    k = length(dyn)

    w = _rs(sc.av1, naff)
    At = _rs(sc.aa1, naff, naff)
    @inbounds for i in 1:naff
        w[i] = x̄[i]
    end
    @inbounds for j in 1:naff, i in 1:naff
        At[i, j] = record.DRIFT[j, i]
    end
    _solve_square_system_generic!(At, w, sc.piv, naff)
    θ̄CINT = θ̄ca.CINT
    θ̄DRIFT = θ̄ca.DRIFT
    @inbounds for i in 1:naff
        θ̄CINT[i] -= w[i]
    end
    @inbounds for j in 1:naff, i in 1:naff
        θ̄DRIFT[i, j] -= w[i] * record.mean[j]
    end

    Ad = _rs(sc.kk1, k, k)
    X̄ = _rs(sc.kk4, k, k)
    @inbounds for j in 1:k, i in 1:k
        Ad[i, j] = record.DRIFT[dyn[i], dyn[j]]
        X̄[i, j] = P̄[dyn[i], dyn[j]]
    end
    Ād = sc.kk7
    Q̄d = sc.kk8
    _ctsem_lyap_pullback!(Ād, Q̄d, sc.kk9, sc.kk10, Ad, record.X, X̄,
        aws.lyap_buffer)
    @inbounds for j in 1:k, i in 1:k
        θ̄DRIFT[dyn[i], dyn[j]] += Ād[i, j]
    end
    Qc_bar = _rs(sc.nn5, n, n)
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Q̄d[i, j]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, θ̄ca, n)

    @inbounds for i in 1:naff
        x̄[i] = zero(T)
    end
    @inbounds for j in 1:n, i in 1:naff
        P̄[i, j] = zero(T)
        P̄[j, i] = zero(T)
    end
    return nothing
end

"""
    _reverse_stationary_elements!(x̄, P̄, θ̄ca, record, dyn, n, aws, se)

Undo `_ctsem_stationary_elements_mean!` and `_ctsem_stationary_elements_cov!`.

  * the mean, `m = -A \\ c` over the closure system (`_ctsem_closure_drift`):
    `w = A' \\ m̄`, `m̄` the marked latents' cotangents, which end here; then
    `c̄ -= w` and `Ā -= w m'` within the closure. The decoupled entries are
    constants.
  * the covariance, `P = M P₀ M'`, `M` block diagonal with `M_B = Lx Lp⁻¹`:
    `P̄₀ = M' P̄ M`, `M̄ = (P̄ + P̄') M P₀` on each block, `L̄x = M̄ Lp⁻ᵀ` and
    `L̄p = -M' L̄x`, then the two Cholesky pullbacks (`_ctsem_chol_pullback!`):
    `Lp`'s into `P̄₀`, `Lx`'s into the stationary covariance's cotangent,
    which goes through the Lyapunov pullback as in `_reverse_stationary!`,
    within the closure.
"""
function _reverse_stationary_elements!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMStationaryRecord{T}, dyn::AbstractVector{Int}, n::Int, aws,
    se::CTSEMStationaryElements) where {T}
    sc = aws.reverse_scratch
    naff = aws.affine_dim
    k = length(dyn)
    inclosure = se.inclosure
    θ̄CINT = θ̄ca.CINT
    θ̄DRIFT = θ̄ca.DRIFT

    if !isempty(se.means)
        w = _rs(sc.av1, naff)
        At = _rs(sc.aa1, naff, naff)
        fill!(w, zero(T))
        @inbounds for i in se.means
            w[i] = x̄[i]
            x̄[i] = zero(T)
        end
        @inbounds for j in 1:naff, i in 1:naff
            At[i, j] = _ctsem_closure_drift(record.DRIFT, inclosure, j, i)
        end
        _solve_square_system_generic!(At, w, sc.piv, naff)
        @inbounds for i in 1:naff
            _ctsem_in_closure(inclosure, i) || continue
            θ̄CINT[i] -= w[i]
            for j in 1:naff
                _ctsem_in_closure(inclosure, j) || continue
                θ̄DRIFT[i, j] -= w[i] * record.mean[j]
            end
        end
    end
    isempty(se.blocks) && return nothing

    Pold = record.Pold
    M = record.M
    G = _rs(sc.nn1, n, n)
    @inbounds for j in 1:n, i in 1:n
        G[i, j] = P̄[i, j] + P̄[j, i]
    end
    # Y = M P₀: the blocks' rows transformed, the rest as they were.
    Y = _rs(sc.nn2, n, n)
    @inbounds for j in 1:n, i in 1:n
        Y[i, j] = Pold[i, j]
    end
    o = 0
    for B in se.blocks
        b = length(B)
        r = (o + 1):(o + b)
        mb = view(M, r, r)
        @inbounds for j in 1:n, p in 1:b
            s = zero(T)
            for q in 1:p
                s += mb[p, q] * Pold[B[q], j]
            end
            Y[B[p], j] = s
        end
        o += b
    end
    # M̄ on each block, packed as M is: the block's rows of G times its
    # columns of Y.
    Mbar = _rs(sc.nn3, n, n)
    o = 0
    for B in se.blocks
        b = length(B)
        r = (o + 1):(o + b)
        mbar = view(Mbar, r, r)
        @inbounds for q in 1:b, p in 1:b
            s = zero(T)
            for l in 1:n
                s += G[B[p], l] * Y[l, B[q]]
            end
            mbar[p, q] = s
        end
        o += b
    end
    # P̄₀ = M' P̄ M: the blocks' rows, then their columns.
    tmp = _rs(sc.nn4, n, n)
    o = 0
    for B in se.blocks
        b = length(B)
        r = (o + 1):(o + b)
        mb = view(M, r, r)
        @inbounds for j in 1:n, p in 1:b
            s = zero(T)
            for q in p:b
                s += mb[q, p] * P̄[B[q], j]
            end
            tmp[p, j] = s
        end
        @inbounds for j in 1:n, p in 1:b
            P̄[B[p], j] = tmp[p, j]
        end
        o += b
    end
    o = 0
    for B in se.blocks
        b = length(B)
        r = (o + 1):(o + b)
        mb = view(M, r, r)
        @inbounds for i in 1:n, p in 1:b
            s = zero(T)
            for q in p:b
                s += P̄[i, B[q]] * mb[q, p]
            end
            tmp[p, i] = s
        end
        @inbounds for i in 1:n, p in 1:b
            P̄[i, B[p]] = tmp[p, i]
        end
        o += b
    end
    # The factors' pullbacks: Lp's into P̄₀, Lx's into the stationary
    # covariance's cotangent over the diffusion states.
    X̄ = _rs(sc.kk4, k, k)
    fill!(X̄, zero(T))
    o = 0
    for B in se.blocks
        b = length(B)
        r = (o + 1):(o + b)
        lx = view(record.LX, r, r)
        lp = view(record.LP, r, r)
        mb = view(M, r, r)
        mbar = view(Mbar, r, r)
        Lxbar = _rs(sc.Abar, b, b)
        Lpbar = _rs(sc.Pbar_new, b, b)
        W = _rs(sc.Ps, b, b)
        Ab = _rs(sc.nn5, b, b)
        @inbounds for q in 1:b, p in 1:b
            Lxbar[p, q] = mbar[p, q]
        end
        _ctsem_right_solve_lower_t!(Lxbar, lp, b)
        @inbounds for q in 1:b, p in 1:b
            s = zero(T)
            for l in p:b
                s += mb[l, p] * Lxbar[l, q]
            end
            Lpbar[p, q] = -s
        end
        _ctsem_chol_pullback!(Ab, lp, Lpbar, b, W)
        @inbounds for q in 1:b, p in 1:b
            P̄[B[p], B[q]] += Ab[p, q]
        end
        _ctsem_chol_pullback!(Ab, lx, Lxbar, b, W)
        @inbounds for q in 1:b, p in 1:b
            pp = findfirst(==(B[p]), dyn)
            qq = findfirst(==(B[q]), dyn)
            (pp === nothing || qq === nothing) && continue
            X̄[pp, qq] += Ab[p, q]
        end
        o += b
    end

    Ad = _rs(sc.kk1, k, k)
    @inbounds for j in 1:k, i in 1:k
        Ad[i, j] = _ctsem_closure_drift(record.DRIFT, inclosure, dyn[i], dyn[j])
    end
    Ād = sc.kk7
    Q̄d = sc.kk8
    _ctsem_lyap_pullback!(Ād, Q̄d, sc.kk9, sc.kk10, Ad, record.X, X̄,
        aws.lyap_buffer)
    Qc_bar = _rs(sc.nn5, n, n)
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        (_ctsem_in_closure(inclosure, dyn[i]) &&
            _ctsem_in_closure(inclosure, dyn[j])) || continue
        θ̄DRIFT[dyn[i], dyn[j]] += Ād[i, j]
        Qc_bar[dyn[i], dyn[j]] = Q̄d[i, j]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, θ̄ca, n)
    return nothing
end

"""
Undo one TD-predictor impulse: `x⁺ = x + TDPREDEFFECT td`, `P⁺ = Jtd P Jtd'`.

Temporaries from `sc` and products by the engine's kernels, as in
`_reverse_predict!`: this runs once per row with a TD predictor, and the BLAS
products it used to take lock (`small_linalg.jl`).
"""
function _reverse_td!(x̄::Vector{T}, P̄::Matrix{T}, θ̄ca,
    record::CTSEMTDRecord{T}, n::Int, sc::CTSEMReverseScratch{T}) where {T}
    td = record.tdpreds
    θ̄TD = θ̄ca.TDPREDEFFECT
    @inbounds for j in eachindex(td), i in 1:n
        θ̄TD[i, j] += x̄[i] * td[j]
    end
    Ps = _symmetrize_into!(_rs(sc.Ps, n, n), P̄)
    Jtd = record.Jtd
    # Jtd_bar = Ps Jtd P' + Ps' Jtd P; `Ps` is symmetric, so both terms start
    # from `Ps Jtd`.
    PJ = _ctsem_mul!(_rs(sc.nn1, n, n), Ps, Jtd)
    Jtd_bar = _ctsem_mulNT!(_rs(sc.nn2, n, n), PJ, record.P_in)
    _ctsem_mul!(Jtd_bar, PJ, record.P_in, one(T), one(T))
    θ̄Jtd = θ̄ca.Jtd
    @inbounds for j in 1:n, i in 1:n
        θ̄Jtd[i, j] += Jtd_bar[i, j]
    end
    JP = _ctsem_mulTN!(_rs(sc.nn1, n, n), Jtd, Ps)
    _ctsem_mul!(P̄, JP, Jtd)
    # The mean update is a pure translation, so x̄ passes through unchanged.
    return nothing
end

"""
    _reverse_update!(x̄, P̄, Θ̄, θ̄ca, record, n, sc, hyb, xp̄)

Undo one measurement update, including its log-likelihood contribution. A
relinearised update (`record.relin > 0`, `_ekf_iterated_rows!`) scatters its
loading and intercept back to the cells they were built from, and its prior
mean's cotangent waits in `xp̄` (`_relin_route!`).

Forward (on the observed subset, with `ε` the Stan-matching ridge):

    Pr  = P + εI;   R̃ = R + εI
    S   = sym(H Pr H' + R) + εI  = H Pr H' + R̃
    ỹ   = y - (Λ x + μ)
    K   = Pr H' S⁻¹;   M = I - K H
    x⁺  = x + K ỹ
    P⁺  = M P M' + K R K'                (Joseph form, on the unridged P and R)
    ll  = -½ (m log 2π + logdet S + ỹ'S⁻¹ỹ)

The seed for `ll` is 1: the reverse pass differentiates the summed
log-likelihood, so every row contributes with unit weight.

# Stable in an explosive direction

After a long interval with a positive drift eigenvalue the prediction is huge
in one direction (1e10 on the SNSF pilot), so `S` is badly conditioned there:
an explicit `S⁻¹` is wrong there by about 1e-6 absolute, where the cotangent
`H'S̄H` it feeds is about 1e-10, and the reverse prediction then multiplies that
direction by up to 1e10. Derivatives of the Joseph form taken through `M P`
with `P` huge cancel the same way. So nothing here forms `S⁻¹` or multiplies by
`Pr`: with `S = U'U`, `e = U⁻ᵀỹ`, `a = S⁻¹ỹ`, `Kb = Pr H'U⁻¹`, `Wh = U⁻ᵀH`,

    P̄  = -½(Wh'Wh - (Wh'e)(Wh'e)') + sym(x̄⁺ (Wh'e)') - sym(H'K'x̄⁺ (Wh'e)') + M'P̄⁺M
    H̄  = -U⁻¹(I - ee')Kb' + a (P⁺x̄⁺)' - (K'x̄⁺)(Kỹ)' - 2K'P̄⁺P⁺
    R̄  = -½U⁻¹(I - ee')U⁻ᵀ - sym((K'x̄⁺) a') + K'P̄⁺K
    ỹ̄  = K'x̄⁺ - a

using the Joseph form's stationarity in `K` at the optimal gain (so its `K`
path contributes nothing) and `M Pr = P⁺` (so products with `Pr` become
products with `P⁺`). Both are exact for the posterior on the ridged `Pr` and
`R̃`, which is what is differentiated; the forward forms it on the unridged
`P` and `R` (Stan's choice), 1e-10 away, so the gradient is that of a
posterior 1e-10 from the one the forward carries. Putting the forward on the
ridged pair instead was tried: it moved a unit sitting exactly at curvature one
across the Laplace conditioning counts. Checked on the regression fixture's
explosive row against BigFloat: the amplified `P̄` to 3e-6 relative where the
previous derivation was 24% off, `H̄` to 2e-6 where it was 28% off.
"""
function _reverse_update!(x̄::Vector{T}, P̄::Matrix{T}, Θ̄::Matrix{T}, θ̄ca,
    record::CTSEMUpdateRecord{T}, n::Int, sc::CTSEMReverseScratch{T},
    hyb::Union{Nothing,CTSEMHybridTape{T}}, xp̄::Vector{T}) where {T}
    o = record.observed
    m = length(o)
    H = record.H
    Λ = record.Lambda
    R = record.R
    x = record.state_in

    # Every temporary is a view of a scratch buffer; see `CTSEMReverseScratch`.
    # The names and the order are the derivation's, unchanged -- what changed is
    # only where the memory comes from.
    Pr     = _rs(sc.Pr, n, n)
    PHt    = _rs(sc.PHt, n, m)
    S      = _rs(sc.S, m, m)
    Sinv   = _rs(sc.Sinv, m, m)
    G      = _rs(sc.G, n, m)
    M      = _rs(sc.M, n, n)
    S̄      = _rs(sc.Sbar, m, m)
    S̄0     = _rs(sc.Sbar0, m, m)
    R̄      = _rs(sc.Rbar, m, m)
    M̄      = _rs(sc.Mbar, n, n)
    P̄_new  = _rs(sc.Pnew, n, n)
    Ps     = _rs(sc.Ps, n, n)
    Ḡ      = _rs(sc.Gbar, n, m)
    H̄      = _rs(sc.Hbar, m, n)
    Λ̄      = _rs(sc.Lbar, m, n)
    PHt_bar = _rs(sc.PHtbar, n, m)
    ỹ      = _rs(sc.ytilde, m)
    ỹ̄      = _rs(sc.ybar, m)
    α      = _rs(sc.alpha, m)
    ᾱ      = _rs(sc.alphabar, m)
    β      = _rs(sc.beta, m)
    x̄_new  = _rs(sc.xnew, n)
    nn1    = _rs(sc.nn1, n, n)
    nn2    = _rs(sc.nn2, n, n)
    nm1    = _rs(sc.nm1, n, m)
    mm1    = _rs(sc.mm1, m, m)
    mm2    = _rs(sc.mm2, m, m)

    # Recompute the forward intermediates, in the same order and with the same
    # ridges as `_ekf_masked_update_step!`, so the derivative is taken at the
    # values the primal actually used.
    copyto!(Pr, record.P_in)
    _ridge_diagonal!(Pr, n, _CTSEM_RIDGE)
    _ctsem_mulNT!(PHt, Pr, H)                       # PHt = Pr H'
    _ctsem_mul!(mm1, H, PHt)                                 # mm1 = H PHt
    mm1 .+= R
    _symmetrize_into!(S, mm1)                         # S = sym(H PHt + R)
    _ridge_diagonal!(S, m, _CTSEM_RIDGE)
    # `S` is the innovation covariance: symmetric, positive definite, and the
    # size of the manifest vector. `inv` on it is a LAPACK `getrf!`/`getri!`
    # pair, which at this size is mostly OpenBLAS's process-global buffer lock
    # -- it was 6.4% of the reverse pass's samples and it does not thread. The
    # engine's own Cholesky plus `m` triangular solves is the same inverse
    # without the lock. See `small_linalg.jl`.
    copyto!(mm2, S)
    # The same size ternary the forward uses (`kalman_filters.jl`): above
    # `_CTSEM_SMALL_CHOLESKY` LAPACK's blocking earns its lock back, and using a
    # different factorization from the forward's is itself a way for the reverse
    # to decline an `S` the forward accepted.
    # One return type from both branches -- see the forward's copy of this in
    # `kalman_filters.jl` for why the `Union` cost what it did.
    F = if m <= _CTSEM_SMALL_CHOLESKY[]
        _ctsem_cholesky(mm2, m)
    else
        CTSEMCholesky(mm2, m, issuccess(cholesky!(mm2, check=false)))
    end
    # `Sinv` holds `U⁻¹` here, not `S⁻¹`: see the docstring for why nothing
    # below forms `S⁻¹`.
    Ui = Sinv
    if issuccess(F)
        _ctsem_cholesky_uinv!(Ui, F)
    else
        # Not `inv(S)`. The reverse recomputes `S` with its own kernels rather
        # than reusing the forward's factorization, so it can fail here on a
        # marginally definite `S` the forward accepted -- and `inv` of a
        # symmetric matrix that just failed to factorize returns something
        # large but finite, so the whole reverse pass would run on to a finite,
        # wrong gradient that `fg!` accepts and feeds to the L-BFGS curvature
        # update. Poison instead, which is how the forward already treats this
        # state: an invalid trial point.
        #
        # `θ̄ca` and not only `x̄`/`P̄`/`Θ̄`, because `_ctsem_regular_pullback!`
        # drops cotangents at non-mutable positions -- a poison confined to the
        # state cotangents could be filtered out of a model whose measurement
        # matrices are all fixed.
        fill!(x̄, _ctsem_nan(T))
        fill!(P̄, _ctsem_nan(T))
        fill!(Θ̄, _ctsem_nan(T))
        fill!(θ̄ca, _ctsem_nan(T))
        return nothing
    end
    copyto!(ỹ, record.manifestmeans)
    _ctsem_mulvec!(ỹ, Λ, x, one(T), one(T))                     # ỹ = Λ x + μ
    ỹ .= record.y .- ỹ

    # Names for the buffers this derivation uses; see the docstring.
    e   = α                                   # U⁻ᵀ ỹ
    a   = β                                   # S⁻¹ ỹ
    Kx  = ᾱ                                   # K' x̄⁺
    Kb  = PHt_bar                             # Pr H' U⁻¹
    K   = G                                   # Pr H' S⁻¹
    Wh  = Λ̄                                   # U⁻ᵀ H, until Λ̄ is formed
    Pp  = M̄                                   # P⁺
    We  = view(_rs(sc.nn3, n, n), :, 1)       # Wh' e = H' a
    HKx = view(_rs(sc.nn4, n, n), :, 1)       # H' K' x̄⁺
    Ky  = view(_rs(sc.nn5, n, n), :, 1)       # K ỹ
    Ppx = _rs(sc.nv1, n)                      # P⁺ x̄⁺
    oneT = one(T)

    _ctsem_mulTvec!(e, Ui, ỹ)
    _ctsem_mulvec!(a, Ui, e)
    _ctsem_mul!(Kb, PHt, Ui)
    _ctsem_mulNT!(K, Kb, Ui)
    # The integrated deviations' cotangents of the gain, the innovation
    # covariance, the innovation and LAMBDA's rows, folded in below where the
    # base reverse handles each (`_hybrid_reverse_update!`).
    hyb === nothing || _hybrid_reverse_update!(hyb, record, Ui, K, a, n)
    _ctsem_mul!(M, K, H)                                     # M = I - K H
    M .= .-M
    @inbounds for i in 1:n; M[i, i] += oneT; end
    _ctsem_mulTN!(Wh, Ui, H)
    _ctsem_mul!(nn1, M, Pr)                                  # P⁺ = M Pr M' + K R̃ K'
    _ctsem_mulNT!(Pp, nn1, M)
    _ctsem_mul!(nm1, K, R)
    _ctsem_mulNT!(Pp, nm1, K, oneT, oneT)
    _ctsem_mulNT!(Pp, K, K, T(_CTSEM_RIDGE), oneT)
    _symmetrize_into!(Ps, P̄)
    _ctsem_mulTvec!(Kx, K, x̄)
    _ctsem_mulTvec!(We, Wh, e)
    _ctsem_mulTvec!(HKx, H, Kx)
    _ctsem_mulvec!(Ky, Kb, e)
    _ctsem_mulvec!(Ppx, Pp, x̄)

    # --- P̄ of the prediction
    _ctsem_mulTN!(P̄_new, Wh, Wh, -T(0.5), false)
    @inbounds for j in 1:n, i in 1:n
        P̄_new[i, j] += T(0.5) * (We[i] * We[j] + x̄[i] * We[j] + We[i] * x̄[j] -
            HKx[i] * We[j] - We[i] * HKx[j])
    end
    _ctsem_mul!(nn1, Ps, M)
    _ctsem_mulTN!(P̄_new, M, nn1, oneT, oneT)                # + M' P̄⁺ M

    # The ridge. The forward's posterior is on the unridged P and R, which is
    # P⁺ - ε(MM' + KK'), so the correction's own derivative is added here,
    # exactly: its cotangents are ε-scaled, so their rounding cannot matter
    # even in an explosive direction. Without it the gradient was that of a
    # posterior 1e-10 away, which `test_continuation.jl`'s 1e-10 checks of
    # the exact Hessian against forward mode saw.
    ε = T(_CTSEM_RIDGE)
    Kc = Ḡ                                                   # K̄ of the correction
    Mc = nn2                                                 # M̄ of the correction
    _ctsem_mul!(Kc, Ps, K, -2ε, false)                       # -2ε P̄⁺ K
    _ctsem_mul!(Mc, Ps, M, -2ε, false)                       # -2ε P̄⁺ M
    _ctsem_mulNT!(Kc, Mc, H, -oneT, oneT)                    # through M = I - K H
    # The gain's cotangent from the integrated deviations rides the same
    # K = Pr H' S^-1 path as the ridge correction's.
    hyb === nothing || _hybrid_accumulate!(Kc, hyb.Kbar)
    _ctsem_mul!(nm1, Kc, Ui)
    _ctsem_mul!(nn1, nm1, Wh)                                # K̄ S⁻¹ H
    @inbounds for j in 1:n, i in 1:n
        P̄_new[i, j] += T(0.5) * (nn1[i, j] + nn1[j, i])
    end

    # --- H̄ (the covariance's loading, Jy)
    Z = Λ̄                                                    # Wh is done with
    @inbounds for j in 1:n, i in 1:m
        Z[i, j] = Kb[j, i] - e[i] * Ky[j]                    # (I - ee') Kb'
    end
    _ctsem_mul!(H̄, Ui, Z, -oneT, false)
    @inbounds for j in 1:n, i in 1:m
        H̄[i, j] += a[i] * Ppx[j] - Kx[i] * Ky[j]
    end
    _ctsem_mul!(nn1, Ps, Pp)
    _ctsem_mulTN!(H̄, K, nn1, -T(2), oneT)                   # - 2 K' P̄⁺ P⁺

    # --- R̄
    Y = S̄
    @inbounds for j in 1:m, i in 1:m
        Y[i, j] = Ui[j, i] - e[i] * a[j]                     # (I - ee') U⁻ᵀ
    end
    # -½ U⁻¹ Y with U⁻¹ upper triangular: only its k >= i entries, half the
    # products of a dense multiply (which summed the structural zeros too).
    @inbounds for j in 1:m, i in 1:m
        acc = zero(T)
        for k in i:m
            acc += Ui[i, k] * Y[k, j]
        end
        R̄[i, j] = -T(0.5) * acc
    end
    @inbounds for j in 1:m, i in 1:m
        R̄[i, j] -= T(0.5) * (Kx[i] * a[j] + a[i] * Kx[j])
    end
    _ctsem_mul!(nm1, Ps, K)
    _ctsem_mulTN!(R̄, K, nm1, oneT, oneT)                    # + K' P̄⁺ K

    # The ridge correction's H̄ and R̄, and its P̄ through S (see above).
    # `K̄ U⁻¹ U⁻ᵀ` once, n by m, serves both of its terms: `Pr` on the left for
    # H̄, `-K'` for S̄. Associated the other way -- `K'K̄` first, then the two
    # triangular factors -- S̄ was two m-by-m-by-m products, the cost of the
    # whole update on a row with many indicators, for a term scaled by ε.
    _ctsem_mulTN!(H̄, K, Mc, -oneT, oneT)                    # M = I - K H
    _ctsem_mul!(nm1, Kc, Ui)                                 # K = Pr H' S⁻¹:
    _ctsem_mulNT!(PHt, nm1, Ui)                              #   K̄ S⁻¹
    _ctsem_mul!(nm1, Pr, PHt)                                #   Pr K̄ S⁻¹
    @inbounds for j in 1:n, i in 1:m
        H̄[i, j] += nm1[j, i]
    end
    _ctsem_mulTN!(mm1, K, PHt, -oneT, false)                 # S̄ = -K'K̄ S⁻¹
    _symmetrize_into!(S̄0, mm1)
    hyb === nothing || _hybrid_accumulate!(S̄0, hyb.Rbar)
    _ctsem_mul!(Λ̄, S̄0, H)                                   # S = H Pr H' + R̃
    _ctsem_mulTN!(P̄_new, H, Λ̄, oneT, oneT)
    _ctsem_mul!(H̄, Λ̄, Pr, T(2), oneT)
    R̄ .+= S̄0

    # --- ỹ = y - (Λ x + μ), through x⁺ = x + K ỹ and the log-likelihood
    @inbounds for i in 1:m
        ỹ̄[i] = Kx[i] - a[i]
    end
    hyb === nothing || _hybrid_accumulate!(ỹ̄, hyb.ybar)
    copyto!(x̄_new, x̄)
    _ctsem_outer!(Λ̄, ỹ̄, x)
    Λ̄ .= .-Λ̄
    hyb === nothing || _hybrid_accumulate!(Λ̄, hyb.Lbar)
    _ctsem_mulTvec!(x̄_new, Λ, ỹ̄, -one(T), one(T))

    # --- scatter back into the full-size parameter and Θ cotangents
    θ̄MM = θ̄ca.MANIFESTMEANS
    θ̄LAMBDA = θ̄ca.LAMBDA
    θ̄Jy = θ̄ca.Jy
    if record.relin == 0
        @inbounds for i in 1:m
            θ̄MM[o[i]] -= ỹ̄[i]
            for j in 1:n
                θ̄LAMBDA[o[i], j] += Λ̄[i, j]
                θ̄Jy[o[i], j] += H̄[i, j]
            end
            for j in 1:m
                Θ̄[o[i], o[j]] += R̄[i, j]
            end
        end
    else
        # Relinearised (`_ekf_iterated_rows!`): the loading was Jy and the
        # intercept `c = MANIFESTMEANS + (LAMBDA - Jy) x_lin`, whose cotangent
        # is `-ỹ̄`, all three at the linearisation point.
        xl = record.x_lin
        @inbounds for i in 1:m
            c̄ = -ỹ̄[i]
            θ̄MM[o[i]] += c̄
            for j in 1:n
                θ̄LAMBDA[o[i], j] += c̄ * xl[j]
                θ̄Jy[o[i], j] += Λ̄[i, j] + H̄[i, j] - c̄ * xl[j]
            end
            for j in 1:m
                Θ̄[o[i], o[j]] += R̄[i, j]
            end
        end
    end

    copyto!(x̄, x̄_new)
    copyto!(P̄, P̄_new)
    record.relin > 0 && _relin_route!(x̄, xp̄, record.Lambda_lin, H, ỹ̄, n, m,
        record.relin == 1)
    return nothing
end

"""
    _relin_route!(x̄, xp̄, Λ, H, ỹ̄, n, m, first)

The end of a relinearised step's reverse. `x̄` holds the cotangent of the row's
prior mean, which joins `xp̄`, and leaves as that of the linearisation point,
through the intercept `c = MANIFESTMEANS + (Λ - H) x_lin` whose cotangent is
`-ỹ̄`. At the first step the point is the prior mean itself, so the two meet
there and `xp̄` is spent.
"""
function _relin_route!(x̄::Vector{T}, xp̄::Vector{T}, Λ, H, ỹ̄, n::Int, m::Int,
    first::Bool) where {T}
    @inbounds for j in 1:n
        xp̄[j] += x̄[j]
        acc = zero(T)
        for i in 1:m
            acc -= (Λ[i, j] - H[i, j]) * ỹ̄[i]
        end
        x̄[j] = acc
    end
    if first
        @inbounds for j in 1:n
            x̄[j] += xp̄[j]
            xp̄[j] = zero(T)
        end
    end
    return nothing
end

"""
    _reverse_relin!(x̄, xp̄, P̄, Θ̄, θ̄ca, record, n, sc)

Undo one step of an iterated update (`_ekf_relin_mean!`),

    x⁺ = x_p + Pr H' a,   a = S⁻¹ ỹ,   S = sym(H Pr H' + R) + εI,
    ỹ = y - H x_p - c,    c = MANIFESTMEANS + (LAMBDA - H) x_lin,

with `Pr = P + εI`, as the forward ridges, and the step `x_lin + α (x⁺ - x_lin)`
taken towards it. `x̄` is the cotangent of that step: the
prior mean's share goes to `xp̄`, the covariance's to `P̄` (the prior's, the
covariance at this point of the row), the innovation covariance's to `Θ̄`, the
measurement's to `θ̄ca`, and `x̄` leaves as the linearisation point's
(`_relin_route!`).
"""
function _reverse_relin!(x̄::Vector{T}, xp̄::Vector{T}, P̄::Matrix{T}, Θ̄::Matrix{T},
    θ̄ca, record::CTSEMRelinRecord{T}, n::Int, sc::CTSEMReverseScratch{T}) where {T}
    o = record.observed
    m = length(o)
    H = record.H
    xp = record.state_in
    xl = record.x_lin
    Pr  = _rs(sc.Pr, n, n)
    PHt = _rs(sc.PHt, n, m)
    mm1 = _rs(sc.mm1, m, m)
    S   = _rs(sc.S, m, m)
    mm2 = _rs(sc.mm2, m, m)
    S̄   = _rs(sc.Sbar, m, m)
    SH  = _rs(sc.Lbar, m, n)
    H̄   = _rs(sc.Hbar, m, n)
    ỹ   = _rs(sc.ytilde, m)
    ỹ̄   = _rs(sc.ybar, m)
    a   = _rs(sc.alpha, m)
    w   = _rs(sc.alphabar, m)
    v   = _rs(sc.xnew, n)
    Ha  = _rs(sc.nv1, n)
    keep = view(_rs(sc.nn3, n, n), :, 1)
    oneT = one(T)

    # The forward's intermediates, with its ridges and its factorization.
    copyto!(Pr, record.P_in)
    _ridge_diagonal!(Pr, n, _CTSEM_RIDGE)
    _ctsem_mulNT!(PHt, Pr, H)
    _ctsem_mul!(mm1, H, PHt)
    mm1 .+= record.R
    _symmetrize_into!(S, mm1)
    _ridge_diagonal!(S, m, _CTSEM_RIDGE)
    copyto!(mm2, S)
    F = if m <= _CTSEM_SMALL_CHOLESKY[]
        _ctsem_cholesky(mm2, m)
    else
        CTSEMCholesky(mm2, m, issuccess(cholesky!(mm2, check=false)))
    end
    if !issuccess(F)
        # Poisoned as `_reverse_update!` poisons a declined `S`.
        fill!(x̄, _ctsem_nan(T))
        fill!(xp̄, _ctsem_nan(T))
        fill!(P̄, _ctsem_nan(T))
        fill!(Θ̄, _ctsem_nan(T))
        fill!(θ̄ca, _ctsem_nan(T))
        return nothing
    end
    copyto!(ỹ, record.offset)
    _ctsem_mulvec!(ỹ, H, xp, oneT, oneT)
    ỹ .= record.y .- ỹ
    ldiv!(a, F, ỹ)

    # The step: `α x̄` reaches the mean `x⁺`, the rest the point it left.
    α = T(record.alpha)
    @inbounds for j in 1:n
        keep[j] = (oneT - α) * x̄[j]
        x̄[j] *= α
    end
    _ctsem_mulvec!(v, Pr, x̄)                         # Pr x̄
    _ctsem_mulvec!(w, H, v)                          # ā = H Pr x̄
    ldiv!(ỹ̄, F, w)                                   # ỹ̄ = S⁻¹ ā
    _ctsem_mulTvec!(Ha, H, a)                        # H' a
    @inbounds for j in 1:m, i in 1:m                 # S̄ = -sym(ỹ̄ a')
        S̄[i, j] = -T(0.5) * (ỹ̄[i] * a[j] + a[i] * ỹ̄[j])
    end
    @inbounds for j in 1:n, i in 1:n                 # P̄ += sym(x̄ (H'a)')
        P̄[i, j] += T(0.5) * (x̄[i] * Ha[j] + Ha[i] * x̄[j])
    end
    _ctsem_mul!(SH, S̄, H)
    _ctsem_mulTN!(P̄, H, SH, oneT, oneT)              # P̄ += H' S̄ H
    _ctsem_mul!(H̄, SH, Pr, T(2), false)             # H̄ = 2 S̄ H Pr
    @inbounds for j in 1:n, i in 1:m                 # + a v' + ỹ̄ (x_lin - x_p)'
        H̄[i, j] += a[i] * v[j] + ỹ̄[i] * (xl[j] - xp[j])
    end

    θ̄MM = θ̄ca.MANIFESTMEANS
    θ̄LAMBDA = θ̄ca.LAMBDA
    θ̄Jy = θ̄ca.Jy
    @inbounds for i in 1:m
        θ̄MM[o[i]] -= ỹ̄[i]
        for j in 1:n
            θ̄LAMBDA[o[i], j] -= ỹ̄[i] * xl[j]
            θ̄Jy[o[i], j] += H̄[i, j]
        end
        for j in 1:m
            Θ̄[o[i], o[j]] += S̄[i, j]
        end
    end
    _ctsem_mulTvec!(x̄, H, ỹ̄, -oneT, oneT)           # x_p's: x̄ - H' ỹ̄
    _relin_route!(x̄, xp̄, record.Lambda, H, ỹ̄, n, m, record.first)
    x̄ .+= keep
    return nothing
end

"""
    _ctsem_reverse_tape!(tape, sp, aws, n, m)

Walk one subject's tape backwards, *accumulating* into `aws.theta_bar`, the
cotangent of the materialized parameter matrices.

Deliberately stops there rather than going on to the free parameters: whether
`theta_bar` can be shared across subjects (and the parameter layer therefore
unwound once instead of once per subject) is a property of the model, not of
one tape. See `_ctsem_parameter_layer_shareable`.

`theta_bar` is *not* zeroed here. The caller owns that, because sharing is the
whole point.
"""
function _ctsem_reverse_tape!(tape::CTSEMAdjointTape{T},
    sp::EKFParameters, aws, n::Int, m::Int) where {T}
    x̄ = aws.x_bar
    P̄ = aws.P_bar
    Θ̄ = aws.manifest_cov_bar
    θ̄ = aws.theta_bar
    fill!(x̄, zero(T)); fill!(P̄, zero(T)); fill!(Θ̄, zero(T))
    resize!(tape.xp_bar, n)
    fill!(tape.xp_bar, zero(T))
    # A walk starts with nothing deferred. The end of a walk flushes both, but
    # a walk an exception left would hand its pending cotangent to the next
    # subject's, which with `_defer_diffusion!`'s comparison skipped would
    # merge into it unseen.
    aws.mvar_pending = false
    aws.dvar_pending = false
    aws.series_reverse.noise_bar.pending = false
    aws.series_reverse.intercept_bar.pending = false
    # Constructed here rather than read from the workspace on purpose: a
    # `ComponentVector` view is a free wrapper around `θ̄`, but a workspace
    # field holding one has to be `Any`-typed (its axis type is model
    # specific), and reading it would make every `θ̄ca.DRIFT[i,j] +=` in the
    # loop below a dynamic dispatch. Measured at 45% of the whole reverse
    # pass when this was done the other way.
    θ̄ca = CTSEMPars(θ̄, sp.parameter_axis)
    dyn = aws.diffusion_state_indices

    for entry_index in length(tape.program):-1:1
        kind, index = tape.program[entry_index]
        if kind === :update
            _reverse_update!(x̄, P̄, Θ̄, θ̄ca, tape.updates[index], n, aws.reverse_scratch,
                tape.hyb.active && tape.hyb.q > 0 ? tape.hyb : nothing, tape.xp_bar)
        elseif kind === :relin
            # A step of an iterated update, between the row's final update and
            # its update-group evaluations (`_ekf_iterated_rows!`).
            _reverse_relin!(x̄, tape.xp_bar, P̄, Θ̄, θ̄ca, tape.relins[index], n,
                aws.reverse_scratch)
        elseif kind === :hybrid
            # The integrated deviations' closing term; last in the forward, so
            # first here (`CTSEMHybrid`).
            _reverse_hybrid_marginal!(tape.hyb, n, _hybrid_items(sp.hybrid), sp.hybrid)
        elseif kind === :binary
            # After the Gaussian block of the same row, because the forward
            # applied binary observations before it. By the workspace's type,
            # so a model with no categorical indicator, whose tape never holds
            # one, does not compile the categorical reverse; see
            # `_ekf_categorical_call`.
            _ekf_categorical_call(aws.ekf_ws, _reverse_binary!, x̄, P̄, θ̄ca,
                tape.binaries[index], n)
        elseif kind === :td
            _reverse_td!(x̄, P̄, θ̄ca, tape.tds[index], n, aws.reverse_scratch)
        elseif kind === :predict
            _reverse_predict!(x̄, P̄, θ̄ca, tape.predicts[index], dyn, n, aws.lyap_buffer, aws,
                tape.hyb.active && tape.hyb.q > 0 ? tape.hyb : nothing)
        elseif kind === :moment
            # Between its substep's prediction and the group at the mean.
            _reverse_moment!(x̄, P̄, θ̄, θ̄ca, tape.moments[index], tape, sp, aws, n, m)
        elseif kind === :group
            # Only a group that writes a JAx cell can consume a queued Fréchet
            # contribution (its reverse zeroes that cell's cotangent).
            # Flushing before *every* group instead would break the batch on
            # ctsem's default model, whose individually varying MANIFESTMEANS
            # puts a group between every pair of prediction substeps in an
            # otherwise entirely linear model.
            (aws.groups_write_jax || aws.groups_write_diffusion) &&
                _flush_series!(aws, θ̄ca, n)
            aws.groups_write_jax && _flush_frechet!(aws)
            aws.groups_write_manifestvar && _flush_manifestvar!(aws, θ̄ca, m)
            aws.groups_write_diffusion && _flush_diffusion!(aws, θ̄ca, n)
            _reverse_group!(θ̄, x̄, tape.groups[index], sp, aws)
        elseif kind === :theta
            # Deferred rather than pushed back here; see `_defer_manifestvar!`.
            _defer_manifestvar!(aws, tape.thetas[index].MANIFESTVAR, Θ̄, θ̄ca, m)
            fill!(Θ̄, zero(T))
        elseif kind === :stationary
            _reverse_stationary!(x̄, P̄, θ̄ca, tape.stationaries[index], dyn, n, aws)
        elseif kind === :stationary_elements
            _reverse_stationary_elements!(x̄, P̄, θ̄ca, tape.stationaries[index], dyn,
                n, aws, sp.stationary_elements)
        elseif kind === :init
            θ̄T0M = θ̄ca.T0MEANS
            @inbounds for i in 1:n
                θ̄T0M[i] += x̄[i]
            end
            fill!(x̄, zero(T))
            # The population block first, while P̄ still holds its
            # cotangent, then mask those rows and columns out so the T0VAR
            # pullback sees only what survived the overwrite. See the note on
            # `_apply_population_block!` in kalman_filters.jl for the forward
            # order this mirrors.
            popidx = aws.sp.population_indices
            if !isempty(popidx)
                k = length(popidx)
                # The forward pass wrote `k_a k_b block[a,b]` into P, so the
                # cotangent picks up the same factors on its way back. Missing
                # them leaves the likelihood right and the gradient wrong.
                popscale = aws.sp.population_scale
                # Scratch rather than three new matrices per subject: nothing
                # else in the walk is live at its start, and `k <= n`.
                sc = aws.reverse_scratch
                popbar = _rs(sc.nn1, k, k)
                @inbounds for b in 1:k, a in 1:k
                    popbar[a, b] = popscale[a] * popscale[b] *
                        P̄[popidx[a], popidx[b]]
                end
                popraw_bar = _rs(sc.nn2, k, k)
                fill!(popraw_bar, zero(T))
                # The POPULATION code, not the model's: the forward pass
                # built this block with it, and a reverse pass differentiating
                # the other construction gives a gradient for a different
                # model than the likelihood.
                _sdcovsqrt2cov_pullback!(popraw_bar, tape.inits[index].POPRAW,
                    _symmetrize_into!(_rs(sc.nn3, k, k), popbar), k;
                    scratch=aws.covsqrt_scratch,
                    covmatcode=aws.sp.population_covmatcode)
                # Written through the flat cotangent by range rather than by
                # name: `θ̄ca.RAWPOPVAR` would have to compile for models
                # whose axis has no such block, and cannot.
                rng = aws.sp.population_range
                popslot = reshape(view(θ̄, rng), k, k)
                @inbounds for b in 1:k, a in 1:k
                    popslot[a, b] += popraw_bar[a, b]
                end
                @inbounds for sidx in popidx
                    for j in 1:n
                        P̄[sidx, j] = zero(T)
                    end
                    for i in 1:n
                        P̄[i, sidx] = zero(T)
                    end
                end
            end
            # Its own workspace slot rather than `aws.covsqrt_scratch`, which
            # the pullback on the next line consumes.
            t0var_bar = _rs(aws.reverse_scratch.t0var_bar, n, n)
            fill!(t0var_bar, zero(T))
            sym = _symmetrize_into!(_rs(aws.reverse_scratch.sym, n, n), P̄)
            _sdcovsqrt2cov_pullback!(t0var_bar, tape.inits[index].T0VAR, sym, n;
                scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
            θ̄T0V = θ̄ca.T0VAR
            @inbounds for j in 1:n, i in 1:n
                θ̄T0V[i, j] += t0var_bar[i, j]
            end
            fill!(P̄, zero(T))
        else
            throw(ArgumentError("adjoint: unknown tape entry $(kind)"))
        end
    end
    # Before the diffusion's release, which the noise table's feeds.
    _flush_series!(aws, θ̄ca, n)
    _flush_manifestvar!(aws, θ̄ca, m)
    _flush_diffusion!(aws, θ̄ca, n)

    return θ̄
end

"""
    _defer_manifestvar!(aws, mat, cov_bar, theta_bar_ca, m)

Hold this row's measurement-covariance cotangent instead of pushing it back.

`sdcovsqrt2cov`'s reverse is linear in the cotangent, so rows sharing a
MANIFESTVAR can sum their cotangents and pay for one pullback rather than one
each. Whether they share it is asked of the matrix itself rather than assumed
from the model: a state-dependent MANIFESTVAR changes between rows, and the
tape stores each row's own copy, so the comparison is exact and costs O(m^2)
against a pullback's O(m^3) plus a correlation square root.

The batch is released by `_flush_manifestvar!`, which the tape walk calls when
the matrix changes, before a transform group that writes MANIFESTVAR cells
(that group's reverse zeroes those cotangents, so a contribution arriving after
it would be attributed to the raw parameter instead of through the transform),
and at the end of the subject.
"""
function _defer_manifestvar!(aws, mat, cov_bar, θ̄ca, m::Int)
    if aws.mvar_pending && _ctsem_same_matrix(aws.mvar_mat, mat, m)
        @inbounds for j in 1:m, i in 1:m
            aws.mvar_bar[i, j] += cov_bar[i, j]
        end
        return nothing
    end
    _flush_manifestvar!(aws, θ̄ca, m)
    @inbounds for j in 1:m, i in 1:m
        aws.mvar_mat[i, j] = mat[i, j]
        aws.mvar_bar[i, j] = cov_bar[i, j]
    end
    aws.mvar_pending = true
    return nothing
end

"""
    _defer_diffusion!(aws, mat, cov_bar, theta_bar_ca, n)

`_defer_manifestvar!` for DIFFUSION, and the same argument applies with more
force: the prediction reverse runs once per *substep*, so a model with a
substep mesh pushes the same matrix back several times per row.
"""
function _defer_diffusion!(aws, mat, cov_bar, θ̄ca, n::Int)
    # No group writing DIFFUSION means one matrix for the whole walk, so the
    # comparison, per substep, can only succeed: it is skipped then.
    if aws.dvar_pending && (!aws.groups_write_diffusion ||
            _ctsem_same_matrix(aws.dvar_mat, mat, n))
        @inbounds for j in 1:n, i in 1:n
            aws.dvar_bar[i, j] += cov_bar[i, j]
        end
        return nothing
    end
    _flush_diffusion!(aws, θ̄ca, n)
    @inbounds for j in 1:n, i in 1:n
        aws.dvar_mat[i, j] = mat[i, j]
        aws.dvar_bar[i, j] = cov_bar[i, j]
    end
    aws.dvar_pending = true
    return nothing
end

"""
    _flush_series!(aws, θ̄ca, n)

Release the cotangents the series route's intervals left on their tables'
powers (`SeriesTable`): pull them through the power chains once, into JAx's
diffusing and affine blocks and, for the noise, into the diffusion covariance,
which goes on to `_defer_diffusion!` against the raw DIFFUSION the intervals
read. Linear in what was accumulated, so releasing once for a run of intervals
on one table is exact. The tape walk calls it when a table's drift changes,
before a transform group that writes JAx or DIFFUSION cells (whose reverse
zeroes those cotangents), and at the end of the subject.
"""
function _flush_series!(aws, θ̄ca, n::Int)
    _flush_series_noise!(aws, θ̄ca, n)
    _flush_series_intercept!(aws, θ̄ca)
    return nothing
end

function _flush_series_noise!(aws, θ̄ca, n::Int)
    sr = aws.series_reverse
    sr.noise_bar.pending || return nothing
    T = eltype(sr.Jbar)
    k = sr.noise.p
    dyn = aws.diffusion_state_indices
    Jb = view(sr.Jbar, 1:k, 1:k)
    Qb = view(sr.Qbar, 1:k, 1:k)
    fill!(Jb, zero(T))
    fill!(Qb, zero(T))
    _series_table_pullback!(sr.Jbar, sr.Qbar, sr.noise, sr.noise_bar)
    θ̄JAx = θ̄ca.JAx
    @inbounds for j in 1:k, i in 1:k
        θ̄JAx[dyn[i], dyn[j]] += Jb[i, j]
    end
    Qc = sr.Qfull
    fill!(Qc, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc[dyn[i], dyn[j]] = Qb[i, j]
    end
    _defer_diffusion!(aws, sr.diffusion, Qc, θ̄ca, n)
    return nothing
end

function _flush_series_intercept!(aws, θ̄ca)
    sr = aws.series_reverse
    sr.intercept_bar.pending || return nothing
    T = eltype(sr.Jbar)
    p = sr.intercept.p
    Jb = view(sr.Jbar, 1:p, 1:p)
    fill!(Jb, zero(T))
    _series_table_pullback!(sr.Jbar, nothing, sr.intercept, sr.intercept_bar)
    θ̄JAx = θ̄ca.JAx
    @inbounds for j in 1:p, i in 1:p
        θ̄JAx[i, j] += Jb[i, j]
    end
    return nothing
end

"""Push back whatever `_defer_diffusion!` is holding, if anything."""
function _flush_diffusion!(aws, θ̄ca, n::Int)
    aws.dvar_pending || return nothing
    aws.dvar_pending = false
    _sdcovsqrt2cov_pullback!(θ̄ca.DIFFUSION, aws.dvar_mat, aws.dvar_bar, n;
        scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
    return nothing
end

"""Push back whatever `_defer_manifestvar!` is holding, if anything."""
function _flush_manifestvar!(aws, θ̄ca, m::Int)
    aws.mvar_pending || return nothing
    aws.mvar_pending = false
    # Straight into the cotangent: the pullback accumulates rather than
    # overwrites, and writes only the lower triangle and the diagonal, which is
    # exactly what the intermediate buffer used to be copied for.
    _sdcovsqrt2cov_pullback!(θ̄ca.MANIFESTVAR, aws.mvar_mat, aws.mvar_bar, m;
        scratch=aws.covsqrt_scratch, covmatcode=aws.sp.covmatcode)
    return nothing
end

"""Elementwise equality over the leading `d` by `d` block."""
@inline function _ctsem_same_matrix(a::AbstractMatrix, b::AbstractMatrix, d::Int)
    @inbounds for j in 1:d, i in 1:d
        a[i, j] == b[i, j] || return false
    end
    return true
end

"""
    _ctsem_parameter_layer!(values_bar, θ̄, subject_values, sp, aws, tipreds, values)

Unwind the parameter layer: whatever cotangent is left on `all_params` was put
there by the regular (non-state-dependent) transforms, so push it through them
and then through the TI-predictor effects into `values_bar`.

`tipreds` is untyped rather than `::AbstractVector` because a subject with a
sampled (missing) TI predictor cell passes its `TIMissingRecipe` here instead
-- see `_ctsem_ti_pullback!`, which dispatches on it. `values` is the raw
trial point the forward pass ran at, needed only by that recipe method (the
plain-vector one ignores it).
"""
function _ctsem_parameter_layer!(values_bar::AbstractVector{T}, θ̄::AbstractVector{T},
    subject_values::AbstractVector, sp::EKFParameters, aws,
    tipreds, values::AbstractVector{T}) where {T}
    aws.defer_frechet || _flush_frechet!(aws)
    subject_values_bar = aws.subject_values_bar
    resize!(subject_values_bar, length(subject_values))
    fill!(subject_values_bar, zero(T))
    _ctsem_regular_pullback!(subject_values_bar, θ̄, subject_values, sp,
        aws.regular_supports, aws.regular_dual_scratch)
    _ctsem_ti_pullback!(values_bar, subject_values_bar, sp, tipreds, values)
    return values_bar
end

"""Keep the affine offset's cotangent for `_reverse_moment!` (`CTSEMMomentScratch`)."""
@inline function _moment_cbar!(aws, bar, k::Int)
    cbar = aws.moment.cbar
    fill!(cbar, zero(eltype(cbar)))
    @inbounds for i in 1:k
        cbar[i] = bar[i]
    end
    return nothing
end

"""
    _reverse_moment!(x̄, P̄, θ̄, θ̄ca, record, tape, sp, aws, n, m)

Undo a moment-matched substep's cubature (`_ekf_moment_cells!`). Forward, with
`P + eps I = U'U`, the offsets `d_i = +-sqrt(n) U[j, :]'` and points
`x_i = m + d_i`, weights `w = 1 / 2n`:

    cell_k = sum_i w cell_k(x_i)        the group's JAx, DRIFT and CINT cells
    c      = sum_i w DRIFT(x_i) d_i     the mean's correction

and the group's other cells hold their values at the mean. The prediction's
reverse has left the averaged cells' cotangent in `θ̄` and `c`'s in
`aws.moment.cbar`. Each point gets its share in its cells, `c`'s in its DRIFT
cells and its offset, and its group record is reversed with a state cotangent
of its own, which reaches the mean and, through the offset, the factor and so
`P`. The other cells' cotangent is set aside while the points run -- what they
send those cells flows on to the group at the mean, which wrote the values
they read -- and handed back for that group's own reverse, next on the walk.
"""
function _reverse_moment!(x̄::Vector{T}, P̄::Matrix{T}, θ̄::Vector{T}, θ̄ca,
    record::CTSEMMomentRecord{T}, tape::CTSEMAdjointTape{T}, sp::EKFParameters, aws,
    n::Int, m::Int) where {T}
    # As before any group that writes these cells; see `_ctsem_reverse_tape!`.
    (aws.groups_write_jax || aws.groups_write_diffusion) && _flush_series!(aws, θ̄ca, n)
    aws.groups_write_jax && _flush_frechet!(aws)
    aws.groups_write_manifestvar && _flush_manifestvar!(aws, θ̄ca, m)
    aws.groups_write_diffusion && _flush_diffusion!(aws, θ̄ca, n)
    ews = aws.ekf_ws
    avg = ews.moment_avg
    oth = ews.moment_oth
    ms = aws.moment
    cellbar = ms.cellbar
    savedoth = ms.savedoth
    @inbounds for k in eachindex(avg)
        cellbar[k] = θ̄[avg[k]]
        θ̄[avg[k]] = zero(T)
    end
    @inbounds for k in eachindex(oth)
        savedoth[k] = θ̄[oth[k]]
        θ̄[oth[k]] = zero(T)
    end
    cbar = ms.cbar
    U = record.U
    D = record.D
    DRIFTs = record.DRIFTs
    np = record.npoints
    w = inv(T(np))
    r = sqrt(T(n))
    Ubar = ms.Ubar
    fill!(Ubar, zero(T))
    mbar = ms.mbar
    fill!(mbar, zero(T))
    dbar = ms.dbar
    xs = ms.xs
    @inbounds for i in np:-1:1
        for k in eachindex(avg)
            θ̄[avg[k]] += w * cellbar[k]
        end
        for (k, a, b) in ms.drift_rc
            θ̄[avg[k]] += w * cbar[a] * D[b, i]
        end
        offset = (i - 1) * n
        for b in 1:n
            acc = zero(T)
            for a in 1:n
                acc += DRIFTs[a, offset + b] * cbar[a]
            end
            dbar[b] = w * acc
        end
        fill!(xs, zero(T))
        _reverse_group!(θ̄, xs, tape.groups[record.first + i - 1], sp, aws)
        for q in 1:n
            mbar[q] += xs[q]
            dbar[q] += xs[q]
        end
        j = (i + 1) >> 1
        scale = isodd(i) ? r : -r
        for q in j:n
            Ubar[j, q] += scale * dbar[q]
        end
    end
    @inbounds for k in eachindex(oth)
        θ̄[oth[k]] += savedoth[k]
    end
    @inbounds for q in 1:n
        x̄[q] += mbar[q]
    end
    _moment_cholesky_pullback!(P̄, U, Ubar, ms.B, ms.Z, n)
    return nothing
end

"""
    _moment_cholesky_pullback!(P̄, U, Ū, B, Z, n)

Add to `P̄` the cotangent of a symmetric `P` through its upper Cholesky factor
`P = U'U`, given the factor's cotangent `Ū` (upper triangles of both read).
With `L = U'`, `P̄ += sym(L^-T Phi(L' L̄) L^-1)` (Murray 2016), `Phi` the lower
triangle with its diagonal halved; here `L' L̄ = U Ū'` and `L^-T = U^-1`.
`B` and `Z` are scratch.
"""
function _moment_cholesky_pullback!(P̄::AbstractMatrix{T}, U::AbstractMatrix{T},
    Ū::AbstractMatrix{T}, B::AbstractMatrix{T}, Z::AbstractMatrix{T}, n::Int) where {T}
    @inbounds for j in 1:n
        for i in 1:(j - 1)
            B[i, j] = zero(T)
        end
        for i in j:n
            acc = zero(T)
            for q in i:n
                acc += U[i, q] * Ū[j, q]
            end
            B[i, j] = i == j ? acc / 2 : acc
        end
    end
    # Z = U^-1 B, by back substitution a column at a time.
    @inbounds for col in 1:n
        for i in n:-1:1
            acc = B[i, col]
            for q in (i + 1):n
                acc -= U[i, q] * Z[q, col]
            end
            Z[i, col] = acc / U[i, i]
        end
    end
    # S' = U^-1 Z', into B, whose Phi is no longer read.
    @inbounds for col in 1:n
        for i in n:-1:1
            acc = Z[col, i]
            for q in (i + 1):n
                acc -= U[i, q] * B[q, col]
            end
            B[i, col] = acc / U[i, i]
        end
    end
    @inbounds for j in 1:n, i in 1:n
        P̄[i, j] += (B[i, j] + B[j, i]) / 2
    end
    return P̄
end

"""Reverse one recorded state-dependent transform group."""
function _reverse_group!(θ̄::Vector{T}, x̄::Vector{T}, record::CTSEMGroupRecord{T},
    sp::EKFParameters, aws) where {T}
    # One call per group rather than a call on whichever triple an `if` picked:
    # the three transform tuples are different types, so the picked one is a
    # union and the call behind it was not static.
    if record.group == 1
        _reverse_group!(θ̄, x̄, record, sp, aws, sp.predict_transforms,
            aws.predict_indices, aws.predict_supports, aws.predict_tipred_supports)
    elseif record.group == 2
        _reverse_group!(θ̄, x̄, record, sp, aws, sp.td_transforms,
            aws.td_indices, aws.td_supports, aws.td_tipred_supports)
    else
        _reverse_group!(θ̄, x̄, record, sp, aws, sp.update_transforms,
            aws.update_indices, aws.update_supports, aws.update_tipred_supports)
    end
    return nothing
end

function _reverse_group!(θ̄::Vector{T}, x̄::Vector{T}, record::CTSEMGroupRecord{T},
    sp::EKFParameters, aws, transforms, indices, supports, tipred_supports) where {T}
    isempty(indices) && return nothing
    # `group_scratch` is a persistent full-length buffer whose entries outside
    # this group's relevant set are never written and stay `NaN` for the life
    # of the workspace (see `CTSEMAdjointWorkspace`). Scatter the compact
    # record into the relevant positions; anything a transform reads that we
    # failed to record therefore reads `NaN` and poisons the gradient loudly
    # rather than silently using a stale number.
    relevant = aws.tape.group_relevant[record.group]
    working = aws.group_scratch
    @inbounds for j in eachindex(relevant)
        working[relevant[j]] = record.params_before[j]
    end
    pars = CTSEMPars(working, sp.parameter_axis)
    ctx = CTSEMRowContext(record.state, pars, record.tdpreds, aws.tipreds,
        record.time, record.dt, 1, record.row)
    _ctsem_complex_group_pullback!(θ̄, x̄, transforms, indices, supports, ctx,
        relevant, aws.dual_context, (supports=tipred_supports, bar=aws.tipreds_bar,
            sampled=aws.tipred_sampled))
    return nothing
end

"""
    _reverse_predict_discrete!(XBAR, PBAR, THETA, record, dyn, n, aws)

Undo one prediction step of a discrete-time model.

Forward:

    A          = JAx
    dINT[i]    = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j]
    dDIFF[D,D] = Qc[D,D]
    x_next     = A x + dINT
    P_next     = A (P + eps I) A' + dDIFF

which is the continuous form with the exponential, the Lyapunov solve and the
intercept solve all collapsed -- a discrete model's DRIFT, CINT and DIFFUSION
are already the one-step quantities. The affine offset runs over every state,
matching the forward pass; `_compute_one_step_form!`'s docstring says why that
requires DRIFT's augmented diagonal to be 1 rather than 0.

Temporaries from the workspace's reverse scratch and products by the engine's
kernels, as `_reverse_predict!` takes them, and the DIFFUSION cotangent
deferred as its is (`_defer_diffusion!`): written the obvious way this was a
dozen arrays and four BLAS products per row.
"""
function _reverse_predict_discrete!(XBAR::Vector{T}, PBAR::Matrix{T}, THETA,
    record::CTSEMPredictRecord{T}, dyn::AbstractVector{Int}, n::Int, aws,
    hyb::Union{Nothing,CTSEMHybridTape{T}}=nothing) where {T}
    A = record.A
    JAx = record.JAx
    x = record.state_in
    k = length(dyn)
    sc = aws.reverse_scratch
    Abar = _rs(sc.Abar, n, n)
    JAx_bar = _rs(sc.nn3, n, n)
    Pbar_new = _rs(sc.Pbar_new, n, n)
    Ps = _rs(sc.Ps, n, n)
    Ptilde = _rs(sc.Pr, n, n)
    nn1 = _rs(sc.nn1, n, n)
    nn2 = _rs(sc.nn2, n, n)
    Qc_bar = _rs(sc.nn5, n, n)
    xbar_new = _rs(sc.xbar_new, n)
    dINT_bar = _rs(sc.nv1, n)
    fill!(JAx_bar, zero(T))
    fill!(Abar, zero(T))

    # --- mean: x_next = A x + dINT, and the integrated deviations'
    # X_next = A X + Xinc (`CTSEMHybrid`), whose intercept increment reads no
    # parameter here.
    _ctsem_outer!(Abar, XBAR, x, one(T), one(T))
    hyb === nothing || _hybrid_reverse_predict!(Abar, JAx_bar, A, JAx, hyb,
        aws.sp.hybrid, record, n, aws.affine_dim, sc, record.dt;
        continuous_time=false)
    copyto!(dINT_bar, XBAR)
    _ctsem_mulTvec!(xbar_new, A, XBAR)

    # --- covariance: P_next = A (P + eps I) A' + dDIFF
    _symmetrize_into!(Ps, PBAR)
    copyto!(Ptilde, record.P_in)
    _ridge_diagonal!(Ptilde, n, _CTSEM_RIDGE)
    _ctsem_mul!(nn1, Ps, A)
    _ctsem_mul!(Abar, nn1, Ptilde, T(2), one(T))
    _ctsem_mulTN!(nn2, A, Ps)
    _ctsem_mul!(Pbar_new, nn2, A)

    # --- dINT[i] = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j]
    θ̄CINT = THETA.CINT
    θ̄DRIFT = THETA.DRIFT
    θ̄JAx = THETA.JAx
    @inbounds for i in 1:n
        θ̄CINT[i] += dINT_bar[i]
    end
    @inbounds for j in 1:n, i in 1:n
        contribution = dINT_bar[i] * x[j]
        θ̄DRIFT[i, j] += contribution
        JAx_bar[i, j] -= contribution
        xbar_new[j] += (record.DRIFT[i, j] - JAx[i, j]) * dINT_bar[i]
    end
    record.moment && _moment_cbar!(aws, dINT_bar, n)

    # --- A is JAx itself, so its cotangent simply adds.
    @inbounds for j in 1:n, i in 1:n
        θ̄JAx[i, j] += JAx_bar[i, j] + Abar[i, j]
    end

    # --- dDIFF[D,D] = Qc[D,D]: the cotangent passes straight through, and
    # Qc = sdcovsqrt2cov(DIFFUSION) consumed only the dynamic block. `Ps` is
    # symmetric, so its dynamic block needs no symmetrising of its own.
    fill!(Qc_bar, zero(T))
    @inbounds for j in 1:k, i in 1:k
        Qc_bar[dyn[i], dyn[j]] = Ps[dyn[i], dyn[j]]
    end
    _defer_diffusion!(aws, record.DIFFUSION, Qc_bar, THETA, n)

    copyto!(XBAR, xbar_new)
    copyto!(PBAR, Pbar_new)
    return nothing
end
