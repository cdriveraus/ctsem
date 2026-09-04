"""
substep_mesh.jl

Choosing how many prediction substeps each observation interval needs.

# What substepping is for

Between two observations the filter forms the discrete-time system from the
continuous one -- `exp(JAx h)`, the intercept and the process covariance -- with
`JAx`, `DRIFT` and `CINT` frozen at the state where the interval starts. For a
linear model that is exact at any `h`. For a state-dependent model it is
exponential Rosenbrock-Euler on the mean and the linearised recursion on the
covariance, and its error grows with how much the drift's linearisation changes
across the step. Splitting the interval into substeps re-linearises at each one.

`maxtimestep` bounds every step by one global number, which is either wasted
work (linear elements, slow dynamics) or not enough (fast nonlinear ones), and
nothing tells the user which. This file measures the need instead.

# The estimator

After a step from `x0` to `x1` with the Jacobian `J = JAx(x0)` the step used,
re-evaluate the drift at `x1` and form the defect

    r = f(x1) - f(x0) - J (x1 - x0),      f(x) = DRIFT(x) x + CINT(x).

`r` is exactly the part of the drift's change over the step that the local
affine model did not see. It is identically zero for a linear model, so a
criterion built on it never substeps one, and it is per state, so it says which
elements are nonlinear without any per-model configuration. It costs one extra
evaluation of the predict-group transforms and no exponential.

The mean error the step commits is about `h |r| / 2`. It is compared with the
state's own predicted standard deviation `sqrt(P_predict[i,i])`: a linearisation
error far inside the uncertainty the filter already carries for that state
cannot move the likelihood. The indicator for a row is the largest such ratio
over its substeps and states, and a row is refined while its indicator exceeds
`tol`, by the factor the `h^2` scaling of the error predicts.

# Frozen, not adaptive

The count per row is decided by `ctsem_auto_substeps` and then held fixed: it is
data on the objective (`_ctsem_substeps`), so the likelihood stays smooth in the
parameters, the tape has a stable shape, and every consumer -- gradient, Hessian,
state-explicit dimension, generation -- integrates the same way. The R side
decides the mesh from the starting values, fits, decides again at the optimum
and refits once if the mesh moved. `maxtimestep` remains a floor on the count.
"""

export ctsem_auto_substeps

"""Default tolerance for the substep indicator; see `ctsem_auto_substeps`."""
# 0.01 from review/INTEGRATOR-substeps: on a simulated process with a
# state-dependent drift, one step per interval cost 5 to 64 log-likelihood
# units and biased estimates by 2 to 5 standard errors; 0.01 brought that
# to about 0.1 to 0.25 units and 0.2 standard errors at two to three substeps
# per interval on average, and 0.003 halved it again at roughly double the cost.
const _CTSEM_SUBSTEP_TOL = Ref(0.01)

"""
    CTSEMSubstepRecorder(T, nrows, n, npredict)

Per-subject storage for the defect estimator: the row indicators, the state the
step started from, the two drift evaluations, and the predict-group parameter
cells saved while the group is re-evaluated at the new state.
"""
mutable struct CTSEMSubstepRecorder{T}
    indicator::Vector{T}
    x0::Vector{T}
    f0::Vector{T}
    f1::Vector{T}
    saved::Vector{T}
end

CTSEMSubstepRecorder(::Type{T}, nrows::Int, n::Int, npredict::Int) where {T} =
    CTSEMSubstepRecorder{T}(zeros(T, nrows), zeros(T, n), zeros(T, n), zeros(T, n),
        zeros(T, npredict))

@inline _begin_substep!(::Nothing, ws) = nothing
@inline _record_substep_defect!(::Nothing, ws, pars, sp, all_params, h, t, ctx) = nothing

@inline function _begin_substep!(rec::CTSEMSubstepRecorder, ws)
    copyto!(rec.x0, ws.state)
    return nothing
end

"""`f .= DRIFT x + CINT` over the first `n` states, with the matrices as they stand."""
@inline function _ctsem_affine_field!(f, pars, x, n::Int)
    @inbounds for i in 1:n
        acc = pars.CINT[i]
        for j in 1:n
            acc += pars.DRIFT[i, j] * x[j]
        end
        f[i] = acc
    end
    return f
end

function _record_substep_defect!(rec::CTSEMSubstepRecorder{T}, ws, pars, sp, all_params,
    h, t::Int, ctx) where {T}
    x0 = rec.x0
    x1 = ws.state
    n = length(x0)
    idx = ws.predict_param_indices
    # The drift at the start state, with the matrices the step actually used.
    _ctsem_affine_field!(rec.f0, pars, x0, n)
    # Re-materialise the predict group at the new state. `ctx.state` is
    # `ws.state` itself, which the step has just moved, so the same context
    # evaluates the transforms at `x1`. The cells are restored afterwards, so
    # nothing downstream -- the TD and measurement groups, or a trace -- can
    # tell this ran.
    @inbounds for (k, i) in enumerate(idx)
        rec.saved[k] = all_params[i]
    end
    apply_complex_transforms_at_indices!(all_params, idx, sp.predict_transforms, ctx)
    _ctsem_affine_field!(rec.f1, pars, x1, n)
    @inbounds for (k, i) in enumerate(idx)
        all_params[i] = rec.saved[k]
    end
    # Defect against the linearisation the step used (JAx is restored above).
    P = ws.P_predict.data
    worst = zero(T)
    @inbounds for i in 1:n
        r = rec.f1[i] - rec.f0[i]
        for j in 1:n
            r -= pars.JAx[i, j] * (x1[j] - x0[j])
        end
        err = abs(r) * h / 2
        scale = sqrt(max(P[i, i], zero(T)) + T(1e-12))
        worst = max(worst, err / scale)
    end
    rec.indicator[t] = max(rec.indicator[t], worst)
    return nothing
end

"""
    ctsem_auto_substeps(objective, values; tol, max_substeps, passes, floor_rule)

Choose a substep mesh for every row of `objective`'s data at the parameter
values `values`.

Returns a NamedTuple: `mesh`, one `Int` per row in the order of the data (entry
1 of each subject is unused); `intervals`, the number of observation intervals;
`refined`, how many of them got more than the floor; `max_substeps`; `total`,
the sum over intervals; `passes`, how many filter passes it took; and `finite`,
false if the likelihood was not finite at `values`, in which case `mesh` is the
floor and nothing was measured.

`tol` is the largest acceptable indicator (see the file header); `max_substeps`
caps a row; `floor_rule` is the `maxtimestep` rule whose counts the mesh never
goes below, `Inf` for one step. A model with no state-dependent predict-group
transform has a defect of exactly zero everywhere and returns the floor without
running the filter.

`fallback`, a full mesh (one entry per row), is where a subject restarts when
the filter is not finite at the floor. Every pass begins from the floor so the
mesh can coarsen between calls, but at the optimum of a stiff model the
one-step filter may not be finite at all -- seen on a cubic-damping process at
dt = 2.5, where the re-mesh came back empty and the fit was lost to one step --
and the mesh the fit was found with is then the right place to start.

Each pass runs the primal filter once per subject with a recorder attached, so
the whole thing costs a few likelihood evaluations. The mesh is not installed
into `objective`; the caller passes it back as the `max_timestep` argument of
`ctsem_objective`, which is how it also reaches everything built from the same
specification later.
"""
function ctsem_auto_substeps(objective::CTSEMObjective, values::AbstractVector;
    tol::Real=_CTSEM_SUBSTEP_TOL[], max_substeps::Integer=64, passes::Integer=4,
    floor_rule::Real=Inf, fallback=nothing)
    tol > 0 || throw(ArgumentError("tol must be positive"))
    max_substeps >= 1 || throw(ArgumentError("max_substeps must be at least 1"))
    floor_rule > 0 || throw(ArgumentError("floor_rule must be positive"))
    subjects = objective.subject_objectives
    x = Vector{Float64}(values)
    meshes = Vector{Vector{Int}}(undef, length(subjects))
    finite = true
    used = 0
    offset = 0
    for (i, sub) in enumerate(subjects)
        nrows = length(sub.timesteps)
        fb = fallback === nothing ? nothing :
            Int[Int(fallback[offset + t]) for t in 1:nrows]
        offset += nrows
        mesh, ok, np = _ctsem_subject_mesh(sub, x, Float64(tol), Int(max_substeps),
            Int(passes), Float64(floor_rule), fb)
        meshes[i] = mesh
        finite &= ok
        used = max(used, np)
    end
    mesh = isempty(meshes) ? Int[] : reduce(vcat, meshes)
    counts = Int[]
    for m in meshes
        append!(counts, view(m, 2:length(m)))
    end
    floors = Int[]
    for sub in subjects
        ts = sub.timesteps
        for t in 2:length(ts)
            push!(floors, _ctsem_substeps(ts[t] - ts[t - 1], Float64(floor_rule), t))
        end
    end
    refined = count(k -> counts[k] > floors[k], eachindex(counts))
    return (mesh=mesh, intervals=length(counts), refined=refined,
        max_substeps=isempty(counts) ? 0 : maximum(counts), total=sum(counts),
        passes=used, finite=finite)
end

function _ctsem_subject_mesh(sub, x::Vector{Float64}, tol::Float64, cap::Int, passes::Int,
    floor_rule::Float64, fallback::Union{Nothing,Vector{Int}}=nothing)
    ts = sub.timesteps
    nrows = length(ts)
    mesh = ones(Int, nrows)
    @inbounds for t in 2:nrows
        mesh[t] = min(cap, _ctsem_substeps(ts[t] - ts[t - 1], floor_rule, t))
    end
    ws = _get_or_init_objective_workspace!(sub, Float64)
    # No predict-group transform: the drift is the same matrix at every state,
    # the defect is identically zero, and the floor is the answer.
    isempty(ws.predict_param_indices) && return (mesh, true, 0)
    tipred_vec = _ctsem_tipred_vector(sub.tipreds, x)
    rec = CTSEMSubstepRecorder(Float64, nrows, _val(ws.state_dim),
        length(ws.predict_param_indices))
    used = 0
    restarted = false
    pass = 0
    while pass < passes
        pass += 1
        fill!(rec.indicator, 0.0)
        ll = _extended_kalman_filter_continuous!(ws, x, sub.data, ts, sub.params,
            sub.tdpreds, tipred_vec, sub.subject, mesh, nothing, nothing, rec)
        used = pass
        if !isfinite(ll)
            # Not finite at this mesh. Once, restart from the mesh the caller
            # was already using; a second failure is reported as such.
            (restarted || fallback === nothing) && return (mesh, false, used)
            restarted = true
            copyto!(mesh, fallback)
            @inbounds for t in 2:nrows
                mesh[t] = min(cap, max(mesh[t], _ctsem_substeps(ts[t] - ts[t - 1], floor_rule, t)))
            end
            pass = 0
            continue
        end
        changed = false
        @inbounds for t in 2:nrows
            indicator = rec.indicator[t]
            if !(indicator <= tol) && mesh[t] < cap
                # The error scales as h^2, so the count wanted is the current one
                # times sqrt(indicator / tol); at least a doubling, so a row that
                # is only just over the line does not creep up one step per pass.
                # Bounded in floating point before it becomes an integer: a row
                # whose state has run away reports an indicator of 1e40 or NaN,
                # and the answer to that is the cap, not an InexactError.
                wanted = mesh[t] * sqrt(indicator / tol)
                mesh[t] = (isfinite(wanted) && wanted < cap) ? max(2 * mesh[t], ceil(Int, wanted)) : cap
                mesh[t] = min(cap, mesh[t])
                changed = true
            end
        end
        changed || break
    end
    return (mesh, true, used)
end
