using DataFrames

# This file is the acceptance-gate harness docs/src/adjoint-roadmap.md commits
# to: "compare it [the adjoint] with current ForwardDiff and central finite
# differences for linear, state/PARS-expression, partial-missing, fully
# missing, and unequal-length multi-subject models. Require ... gradient
# relative error within 1e-4." It now runs the full three-way comparison:
# every scenario checks ForwardDiff against central finite differences (which
# catches primal/dual inconsistencies) *and* the reverse-mode adjoint against
# ForwardDiff (which catches reverse-pass errors).
#
# The two comparisons get very different tolerances on purpose. Finite
# differences are only good to ~1e-4 relative, so they stay at the roadmap's
# stated gate. The adjoint, by contrast, computes the same quantity as
# ForwardDiff by a different route through the same arithmetic, so it should
# agree to near machine precision -- and it does. Holding it to 1e-4 would
# let a real reverse-pass error (a dropped term that happens to be small at
# the test point, say) slip through unnoticed, so it is held to 1e-9 instead.
# If a future model shape genuinely cannot meet that, loosen it deliberately
# for that scenario with a comment saying why, rather than globally.
#
# Each scenario below is built directly through `ekf_from_data_frame`, the
# same R-facing parser test_r_interface.jl exercises, rather than by hand-
# rolling EKFParameters -- so these models are constructed exactly the way
# real R-prepared models are, not a parallel, potentially-divergent test-only
# construction path.
#
# Note on the state-dependent scenario: `ctsem_validate_forward_gradient`
# only checks that ForwardDiff and finite differences agree on *whatever*
# function the code computes -- it does not independently verify that a
# hand-written JAx expression is the true mathematical Jacobian of a chosen
# DRIFT expression. That correctness is Stan's job (see
# `test-stan-julia-parity.R`'s nonlinear-predictor test on the R side); this
# harness's job is catching AD-consistency bugs, so the state-dependent model
# below picks a JAx expression matching its DRIFT expression textually
# rather than one independently re-derived as a Jacobian.

const _ADJOINT_GRADIENT_TOLERANCE = 1e-4   # ForwardDiff vs finite differences
const _ADJOINT_REVERSE_TOLERANCE = 1e-9    # adjoint vs ForwardDiff

"""
Run the three-way gradient comparison for one scenario and assert both gates.

Also asserts the adjoint reproduces the primal log-likelihood: the reverse
pass runs its own traced forward sweep, so a disagreement there would mean the
tape was recorded from a different computation than the one being checked.
"""
function _check_adjoint(objective, values)
    @test isfinite(objective(values))
    result = ContinuousTimeSEM.ctsem_validate_forward_gradient(objective, values)
    @test result.relative_error < _ADJOINT_GRADIENT_TOLERANCE
    @test result.adjoint !== nothing
    @test result.adjoint_relative_error < _ADJOINT_REVERSE_TOLERANCE
    @test isapprox(ContinuousTimeSEM.ctsem_adjoint_gradient(objective, values).value,
        objective(values); rtol=1e-12)
    return result
end

function _adjoint_test_dataframe(;
    drift::AbstractMatrix, jax::AbstractMatrix, cint::AbstractMatrix,
    diffusion::AbstractMatrix, lambda::AbstractMatrix, jy::AbstractMatrix,
    manifestmeans::AbstractMatrix, manifestvar::AbstractMatrix,
    t0var::AbstractMatrix, t0means::AbstractMatrix,
    pars::AbstractMatrix=zeros(1, 1),
    tdpredeffect::Union{Nothing,AbstractMatrix}=nothing,
    jtd::Union{Nothing,AbstractMatrix}=nothing,
    free::Dict=Dict(), predict::Dict=Dict(), update::Dict=Dict(),
)
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    parnumber = Union{Missing,Int}[]
    value = Union{Missing,Float64}[]
    transform = Union{Missing,String}[]
    predicttransform = Union{Missing,String}[]
    updatetransform = Union{Missing,String}[]

    function addmat!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            key = (name, i, j)
            push!(matrices, name); push!(rows, i); push!(cols, j)
            if haskey(free, key)
                pn, tf = free[key]
                push!(parnumber, pn); push!(value, missing); push!(transform, tf)
            else
                push!(parnumber, missing); push!(value, Float64(mat[i, j])); push!(transform, missing)
            end
            push!(predicttransform, get(predict, key, missing))
            push!(updatetransform, get(update, key, missing))
        end
    end

    addmat!(:DRIFT, drift); addmat!(:JAx, jax); addmat!(:CINT, cint)
    addmat!(:DIFFUSION, diffusion); addmat!(:LAMBDA, lambda); addmat!(:Jy, jy)
    addmat!(:MANIFESTMEANS, manifestmeans); addmat!(:MANIFESTVAR, manifestvar)
    addmat!(:T0VAR, t0var); addmat!(:T0MEANS, t0means); addmat!(:PARS, pars)
    if tdpredeffect !== nothing
        addmat!(:TDPREDEFFECT, tdpredeffect); addmat!(:Jtd, jtd)
    end

    DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber, value=value,
        transform=transform, predicttransform=predicttransform, updatetransform=updatetransform)
end

# Ties DRIFT[1,1] and JAx[1,1] to one negative-constrained free parameter and
# DIFFUSION[1,1] to one positive-constrained free parameter -- the smallest
# genuinely free (not all-fixed) linear model.
function _adjoint_linear_1d_parameters()
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DIFFUSION, 1, 1) => (2, "log1p_exp(param[2])"),
        ),
    )
    ekf_from_data_frame(df)
end

# A 2-latent, 2-manifest model with a free cross-effect, so the gradient
# check exercises coupled dynamics, not just a diagonal system.
#
# `continuous` is threaded through to `ekf_from_data_frame` so
# `test_discrete_time.jl` can reuse this exact cell layout for its discrete
# variant (F5) instead of writing a new fixture.
function _adjoint_cross_effect_2d_parameters(; continuous::Bool=true)
    df = _adjoint_test_dataframe(
        drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3],
        cint=[0.0; 0.0;;], diffusion=[0.2 0.0; 0.0 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
        t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DRIFT, 2, 2) => (2, "-log1p_exp(param[2])"),
            (:JAx, 2, 2) => (2, "-log1p_exp(param[2])"),
            (:DRIFT, 1, 2) => (3, "param[3]"),
            (:JAx, 1, 2) => (3, "param[3]"),
        ),
    )
    ekf_from_data_frame(df; continuous_time=continuous)
end

# PARS[1,1] is a free parameter feeding a state-dependent DRIFT/JAx
# expression -- exercises `apply_complex_transforms_at_indices!` under
# ForwardDiff, which is exactly the code path a hand-written or Enzyme-based
# adjoint has to reproduce correctly for nonlinear ctsem models.
function _adjoint_state_dependent_1d_parameters()
    expr = "PARS[1,1] * (1 + 0.1 * state[1])"
    df = _adjoint_test_dataframe(
        drift=[0.0;;], jax=[0.0;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;], pars=[0.0;;],
        free=Dict((:PARS, 1, 1) => (1, "param[1]")),
        predict=Dict((:DRIFT, 1, 1) => expr, (:JAx, 1, 1) => expr),
    )
    ekf_from_data_frame(df)
end

# Each of these three model *types* is built exactly once and reused across
# every @testset that needs it. `ekf_from_data_frame` parses transform
# strings via `eval`, so every call produces uniquely-typed closures and
# forces a full fresh compilation of the EKF pipeline for that type -- fine
# for real model fitting (once per R session) but wasteful if paid once per
# @testset here for what is really the same 3 model types varying only in
# data. Building each sp once cuts this file's compilation cost roughly in
# half versus one fresh `_adjoint_*_parameters()` call per testset.
# Everything the six original scenarios leave fixed -- T0VAR, MANIFESTVAR,
# MANIFESTMEANS, CINT, and the off-diagonal correlation-sqrt entries -- is free
# here, so the reverse paths through `sdcovsqrt2cov` (for all three covariance
# matrices), the discrete-intercept solve, and the measurement-mean term are
# actually exercised rather than multiplied by a structurally zero cotangent.
#
# `continuous` is threaded through for the same reason as in
# `_adjoint_cross_effect_2d_parameters` above (F5).
function _adjoint_free_covariance_2d_parameters(; continuous::Bool=true)
    df = _adjoint_test_dataframe(
        drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3],
        cint=[0.0; 0.0;;], diffusion=[0.2 0.0; 0.05 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.02 0.1],
        t0var=[1.0 0.0; 0.03 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DRIFT, 1, 2) => (2, "param[2]"),
            (:JAx, 1, 2) => (2, "param[2]"),
            (:CINT, 1, 1) => (3, "param[3]"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:DIFFUSION, 2, 1) => (5, "param[5]"),
            (:MANIFESTVAR, 2, 2) => (6, "log1p_exp(param[6])"),
            (:MANIFESTMEANS, 1, 1) => (7, "param[7]"),
            (:T0VAR, 1, 1) => (8, "log1p_exp(param[8])"),
            (:T0VAR, 2, 1) => (9, "param[9]"),
            (:T0MEANS, 2, 1) => (10, "param[10]"),
        ),
    )
    ekf_from_data_frame(df; continuous_time=continuous)
end

# A free TD-predictor effect, so `_apply_td_impulse!` and its reverse are
# covered, together with a TI predictor effect (parameter 2 shifted by the
# subject's TI predictor times the coefficient in parameter 3) so the
# `_materialize_subject_values!` layer contributes a real cross-term.
function _adjoint_td_ti_parameters()
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;],
        tdpredeffect=[0.0;;], jtd=[1.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:TDPREDEFFECT, 1, 1) => (2, "param[2]"),
        ),
    )
    ti_effects = DataFrame(parameter=[2], predictor=[1], coefficient=[3])
    ekf_from_data_frame(df, ti_effects)
end

_ADJOINT_SP_LINEAR_1D = _adjoint_linear_1d_parameters()
_ADJOINT_SP_CROSS_2D = _adjoint_cross_effect_2d_parameters()
_ADJOINT_SP_STATE_DEPENDENT = _adjoint_state_dependent_1d_parameters()
_ADJOINT_SP_FREE_COVARIANCE = _adjoint_free_covariance_2d_parameters()
_ADJOINT_SP_TD_TI = _adjoint_td_ti_parameters()

@testset "Forward-gradient validation: linear 1D" begin
    sp = _ADJOINT_SP_LINEAR_1D
    data = reshape([0.1, -0.2, 0.15], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    values = [0.3, -0.4]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: cross-effect 2D" begin
    sp = _ADJOINT_SP_CROSS_2D
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    values = [0.2, 0.1, -0.3]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: state/PARS-dependent expression" begin
    sp = _ADJOINT_SP_STATE_DEPENDENT
    data = reshape([0.1, -0.2, 0.15], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    values = [0.2]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: partial missingness" begin
    sp = _ADJOINT_SP_CROSS_2D
    # Middle row observes only the first manifest variable.
    data = reshape([0.1, -0.2, 0.15, NaN, -0.1, 0.2], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    values = [0.2, 0.1, -0.3]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: fully missing row" begin
    sp = _ADJOINT_SP_LINEAR_1D
    # Middle row fully missing.
    data = reshape([0.1, NaN, 0.15], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    values = [0.3, -0.4]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: unequal-length multi-subject" begin
    sp = _ADJOINT_SP_LINEAR_1D
    # Subject 1 has 3 observations, subject 2 has 2.
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1, 4], [0.0, 0.5, 1.2, 0.0, 0.7], data)
    values = [0.3, -0.4]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: free covariance and mean parameters" begin
    sp = _ADJOINT_SP_FREE_COVARIANCE
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.05], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2, 2.0], data)
    values = [0.3, 0.2, -0.1, 0.4, 0.15, -0.2, 0.25, 0.1, -0.15, 0.35]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: TD and TI predictors" begin
    sp = _ADJOINT_SP_TD_TI
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2], 1, :)
    tdpreds = reshape([0.0, 1.0, 0.0, 0.0, 0.5, 0.0], 1, :)
    tipreds = reshape([0.4, -0.7], 2, 1)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1, 4],
        [0.0, 0.5, 1.2, 0.0, 0.6, 1.4], data, tdpreds, tipreds)
    values = [0.3, -0.4, 0.25]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: bounded prediction substeps" begin
    # `max_timestep` below the observation spacing forces several substeps per
    # row, so the tape has to record (and the reverse pass has to unwind) more
    # than one prediction per transition.
    sp = _ADJOINT_SP_CROSS_2D
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data,
        zeros(0, 3), zeros(1, 0), 0.2)
    values = [0.2, 0.1, -0.3]
    _check_adjoint(objective, values)
end

@testset "gradient_method plumbing" begin
    # The selector has to survive the trip from R, which marshals character
    # vectors to `String` rather than `Symbol`, so both spellings are accepted
    # and must give the identical answer.
    sp = _ADJOINT_SP_FREE_COVARIANCE
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.05], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2, 2.0], data)
    values = [0.3, 0.2, -0.1, 0.4, 0.15, -0.2, 0.25, 0.1, -0.15, 0.35]

    # Explicitly `:forward` -- the default is now `:adjoint`, so relying on
    # the default here would silently turn this into adjoint-vs-adjoint.
    forward = ContinuousTimeSEM.ctsem_evaluate(objective, values; gradient_method=:forward)
    adjoint_symbol = ContinuousTimeSEM.ctsem_evaluate(objective, values; gradient_method=:adjoint)
    adjoint_string = ContinuousTimeSEM.ctsem_evaluate(objective, values; gradient_method="adjoint")

    @test adjoint_symbol.value ≈ forward.value
    @test isapprox(adjoint_symbol.gradient, forward.gradient; rtol=_ADJOINT_REVERSE_TOLERANCE)
    @test adjoint_string.gradient == adjoint_symbol.gradient
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_evaluate(objective, values; gradient_method=:reverse)

    # The default must be the adjoint (changed 2026-08-21); assert it rather
    # than leaving it implicit, since several checks above and in the R-side
    # parity suite exercise whichever path the default selects.
    @test ContinuousTimeSEM.ctsem_evaluate(objective, values).gradient ==
        adjoint_symbol.gradient

    # `gradient=false` must not build an adjoint workspace or a tape at all.
    @test ContinuousTimeSEM.ctsem_evaluate(objective, values;
        gradient=false, gradient_method=:adjoint).gradient === nothing

    # `ctsem_optimize` must actually drive L-BFGS with the adjoint gradient.
    #
    # The check is that the adjoint-driven run's *reported* value and gradient
    # agree with what forward mode computes at the point it actually reached.
    # It deliberately does not compare optimizer trajectories or final
    # minimizers between the two backends. These small models have a flat
    # plateau -- the drift transform `-log1p_exp(param)` saturates, so the
    # log-likelihood stops changing while the raw parameter runs off, and the
    # gradient underflows to ~1e-97 and then to exactly zero. Both backends
    # produce the *same* (zero) gradient there; L-BFGS's line search then
    # breaks ties on floating-point noise and the two runs stop at different
    # raw values with identical log-likelihoods. Comparing minimizers would
    # therefore test the line search, not the gradient feeding it.
    fit_sp = _ADJOINT_SP_LINEAR_1D
    fit_data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.25, -0.15], 1, :)
    fit_objective = ContinuousTimeSEM.ctsem_objective(fit_sp, [1, 5],
        [0.0, 0.5, 1.2, 2.0, 0.0, 0.6, 1.1, 1.9], fit_data)
    fit = ContinuousTimeSEM.ctsem_optimize(fit_objective, [0.1, 0.1]; maxiter=25,
        gradient_method=:adjoint)
    at_minimizer = ContinuousTimeSEM.ctsem_evaluate(fit_objective, fit.minimizer;
        gradient_method=:forward)
    @test isapprox(fit.maximum_loglik, at_minimizer.value; rtol=1e-12)
    @test isapprox(fit.gradient, at_minimizer.gradient; atol=1e-10)
    # ...and it must have improved on the starting point, i.e. the gradient was
    # actually used rather than silently ignored.
    @test fit.maximum_loglik > fit_objective([0.1, 0.1])
end

@testset "adjoint reports invalid trial points rather than hiding them" begin
    # An invalid point (Cholesky failure inside the update) must poison the
    # value *and* the gradient, the same way the ForwardDiff path does, so
    # `ctsem_optimize`'s isfinite guards reject the trial rather than accepting
    # a silently truncated gradient. This is the roadmap's "must never silently
    # fall back" requirement, checked rather than assumed.
    sp = _ADJOINT_SP_LINEAR_1D
    data = reshape([0.1, -0.2, 0.15], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data)
    # Drives DIFFUSION and the innovation covariance to extreme values.
    extreme = [1e5, -1e5]
    primal = objective(extreme)
    result = ContinuousTimeSEM.ctsem_adjoint_gradient(objective, extreme)
    if isfinite(primal)
        @test isfinite(result.value)
        @test all(isfinite, result.gradient)
    else
        @test !isfinite(result.value)
        @test all(!isfinite, result.gradient)
    end
end

# A 3-state model in which only states 1-2 are dynamic and state 3 is a static
# "carrier" -- the shape `.ctJuliaAugmentRandomEffects` produces for non-T0MEANS
# random effects (Stan's `extralatents`). This is the only scenario where
# `diffusion_state_indices != 1:n`, so it is the only one that exercises the
# reverse predict step's gather/scatter between the full state space and the
# dynamic sub-block (the Lyapunov solve, the discrete-intercept solve, and the
# `dDIFFUSION` term all live on that sub-block while `eJAx` spans everything).
function _adjoint_augmented_carrier_parameters()
    drift = [-0.5 0.2 0.3; 0.1 -0.4 0.0; 0.0 0.0 0.0]
    df = _adjoint_test_dataframe(
        drift=drift, jax=drift, cint=[0.0; 0.0; 0.0;;],
        diffusion=[0.2 0.0 0.0; 0.05 0.15 0.0; 0.0 0.0 0.0],
        lambda=[1.0 0.0 0.0; 0.0 1.0 0.0], jy=[1.0 0.0 0.0; 0.0 1.0 0.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
        t0var=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.02 0.0 1.0], t0means=[0.0; 0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            # The carrier state feeds the dynamic block: this is the cross
            # term that makes the static state matter at all.
            (:DRIFT, 1, 3) => (2, "param[2]"),
            (:JAx, 1, 3) => (2, "param[2]"),
            (:DIFFUSION, 2, 1) => (3, "param[3]"),
            (:T0VAR, 3, 3) => (4, "log1p_exp(param[4])"),
            (:T0VAR, 3, 1) => (5, "param[5]"),
        ),
    )
    # States 1:2 only -- state 3 has structurally zero diffusion.
    ekf_from_data_frame(df,
        DataFrame(parameter=Int[], predictor=Int[], coefficient=Int[]), [1, 2])
end

# A state-dependent MANIFESTVAR, which is the one construct whose reverse
# ordering differs between the first row and later rows: the forward pass
# builds Θ from MANIFESTVAR at initialisation on row 1 (i.e. *before* the
# update-transform group runs) but after that group on every later row. The
# tape records the order actually taken, and this asserts that is right.
function _adjoint_state_dependent_manifestvar_parameters()
    expr = "0.1 + 0.05 * PARS[1,1] * state[1] * state[1]"
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;], pars=[0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:PARS, 1, 1) => (2, "param[2]"),
        ),
        update=Dict((:MANIFESTVAR, 1, 1) => expr),
    )
    ekf_from_data_frame(df)
end

_ADJOINT_SP_AUGMENTED_CARRIER = _adjoint_augmented_carrier_parameters()
_ADJOINT_SP_STATEDEP_MANIFESTVAR = _adjoint_state_dependent_manifestvar_parameters()

@testset "Forward-gradient validation: augmented static carrier state" begin
    sp = _ADJOINT_SP_AUGMENTED_CARRIER
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.05], 2, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2, 2.0], data)
    values = [0.3, 0.2, 0.1, -0.15, 0.25]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: state-dependent MANIFESTVAR" begin
    sp = _ADJOINT_SP_STATEDEP_MANIFESTVAR
    data = reshape([0.1, -0.2, 0.15, 0.3], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2, 2.1], data)
    values = [0.3, 0.4]
    _check_adjoint(objective, values)
end

# Within a transform group the cells are written in ascending flattened-index
# order into the same buffer they read from, and `PARS` sorts *last* in the
# parameter axis. So a state-dependent `DRIFT` expression that references a
# state-dependent `PARS` cell reads that cell's value from the **previous
# row** -- the current row's value has not been written yet when DRIFT is
# evaluated. The reverse pass has to evaluate DRIFT's derivative at that older
# value, which is why `CTSEMGroupRecord` stores the pre-group snapshot and
# `_ctsem_complex_group_pullback!` replays the group forward.
#
# This scenario exists to make that concrete: recording the *post*-group
# parameter values instead makes it fail while every other scenario in this
# file still passes.
function _adjoint_intra_group_ordering_parameters()
    # Complex-transform expressions receive a `CTSEMRowContext`, not the raw
    # parameter vector, so they cannot reference `param[...]`; free parameters
    # reach them through a PARS cell instead. PARS[1,1] is the free drift
    # parameter, PARS[1,2] is the state-dependent cell written *after* DRIFT.
    drift_expr = "-log1p_exp(PARS[1,1]) - 0.1 * PARS[1,2] * PARS[1,2]"
    df = _adjoint_test_dataframe(
        drift=[-0.4;;], jax=[-0.4;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;], pars=zeros(1, 2),
        free=Dict(
            (:DIFFUSION, 1, 1) => (1, "log1p_exp(param[1])"),
            (:PARS, 1, 1) => (2, "param[2]"),
        ),
        predict=Dict(
            # Written first (DRIFT precedes PARS in the parameter axis), and
            # reads PARS[1,2] -- i.e. the value from the *previous* row.
            (:DRIFT, 1, 1) => drift_expr,
            (:JAx, 1, 1) => drift_expr,
            # Written second, from the state and the free PARS[1,1].
            (:PARS, 1, 2) => "PARS[1,1] * (1 + 0.3 * state[1])",
        ),
    )
    ekf_from_data_frame(df)
end

_ADJOINT_SP_INTRA_GROUP = _adjoint_intra_group_ordering_parameters()

@testset "Forward-gradient validation: intra-group transform ordering" begin
    sp = _ADJOINT_SP_INTRA_GROUP
    data = reshape([0.1, -0.2, 0.15, 0.3, -0.05], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1],
        [0.0, 0.4, 0.9, 1.5, 2.2], data)
    values = [0.2, 0.45]
    _check_adjoint(objective, values)
end

@testset "discretization cache cannot go stale" begin
    # `DiscretizationCache` reuses the matrix exponential and the asymptotic
    # diffusion solve across rows, guarded by a comparison against the inputs
    # that produced them. A stale-cache bug would be invisible at a single
    # parameter value, so this drives the same objective through *changing*
    # parameters and irregular time gaps and requires that every evaluation
    # match a pristine objective built for that call alone.
    sp = _ADJOINT_SP_CROSS_2D
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.05, 0.22, -0.11], 2, :)
    # Deliberately irregular: two gaps repeat (so the exp cache hits) and the
    # rest do not (so it misses), inside one subject.
    times = [0.0, 0.5, 1.0, 1.9, 2.4]
    reused = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)

    settings = ([0.2, 0.1, -0.3], [0.7, -0.4, 0.15], [0.2, 0.1, -0.3], [-0.5, 0.9, 0.0])
    for values in settings
        pristine = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)
        @test reused(values) == pristine(values)
        @test isapprox(ContinuousTimeSEM.ctsem_adjoint_gradient(reused, values).gradient,
            ContinuousTimeSEM.ctsem_adjoint_gradient(pristine, values).gradient;
            rtol=1e-12)
    end

    # Interleaving two parameter vectors is the sharpest version: each call
    # must invalidate what the previous one cached.
    a, b = [0.2, 0.1, -0.3], [0.7, -0.4, 0.15]
    reference_a = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)(a)
    reference_b = ContinuousTimeSEM.ctsem_objective(sp, [1], times, data)(b)
    for _ in 1:3
        @test reused(a) == reference_a
        @test reused(b) == reference_b
    end

    # A state-dependent model changes JAx every row, so the cache must miss
    # every time rather than reuse the previous row's exponential.
    nl = ContinuousTimeSEM.ctsem_objective(_ADJOINT_SP_STATE_DEPENDENT, [1],
        [0.0, 0.5, 1.2], reshape([0.1, -0.2, 0.15], 1, :))
    nl_pristine = ContinuousTimeSEM.ctsem_objective(_ADJOINT_SP_STATE_DEPENDENT, [1],
        [0.0, 0.5, 1.2], reshape([0.1, -0.2, 0.15], 1, :))
    @test nl([0.2]) == nl_pristine([0.2])
    _check_adjoint(nl, [0.35])
end
