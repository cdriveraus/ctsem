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

# Shared with `test_discrete_time.jl`; see `adjoint_fixtures.jl`.
isdefined(@__MODULE__, :_adjoint_test_dataframe) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))


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

# F6: a 2-latent TD-predictor scenario with `Jtd[2,1]` free. The only other TD
# scenario (`_adjoint_td_ti_parameters`, above) is 1-state with `Jtd` fixed, so
# `_reverse_td!`'s `Jtd_bar = Ps Jtd P' + Ps' Jtd P` line
# (`adjoint_ekf.jl:675`), where an operand-order or transpose error matters
# most, was never checked at a dimension where `P` is not scalar.
function _adjoint_td_cross_2d_parameters()
    drift = [-0.5 0.3; 0.1 -0.3]
    df = _adjoint_test_dataframe(
        drift=drift, jax=drift, cint=[0.0; 0.0;;],
        diffusion=[0.2 0.0; 0.0 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
        t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
        tdpredeffect=[0.0; 0.0;;], jtd=[1.0 0.0; 0.0 1.0],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:TDPREDEFFECT, 1, 1) => (2, "param[2]"),
            (:Jtd, 2, 1) => (3, "param[3]"),
        ),
    )
    ekf_from_data_frame(df)
end

# A T0MEANS and a T0VAR cell each materialised from *two* raw parameters, which
# is what `.ctJuliaResolveStaticRefs` composes for `T0VAR = 'log1p(par1 +
# par2)'`: the referenced PARS cells' own transforms are substituted into the
# user's expression at model-build time, leaving one regular transform that
# reads both parameters.
#
# The parameter-layer pullback used to require exactly one read per transform
# and threw at workspace construction otherwise, so this shape was refused in
# R before it could reach the engine. It now discovers the support and pushes a
# term back through each index, and this is the scenario that checks the terms
# are both there and both right -- a pullback that kept only the representative
# parameter would still produce a plausible gradient, just one missing
# `param[4]`'s contribution entirely.
#
# `parnumber` for each composed cell is its representative (3), the first
# parameter the expression references, which is what R sends; the support
# `[3, 4]` is discovered rather than declared.
function _adjoint_composed_static_1d_parameters()
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;], pars=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DIFFUSION, 1, 1) => (2, "log1p_exp(param[2])"),
            (:PARS, 1, 1) => (3, "param[3]"),
            (:PARS, 2, 1) => (4, "2 * param[4]"),
            (:T0VAR, 1, 1) => (3, "log1p_exp((param[3]) + (2 * param[4]))"),
            (:T0MEANS, 1, 1) => (3, "(param[3]) - 0.5 * (2 * param[4])"),
        ),
    )
    ekf_from_data_frame(df)
end

_ADJOINT_SP_LINEAR_1D = _adjoint_linear_1d_parameters()
_ADJOINT_SP_CROSS_2D = _adjoint_cross_effect_2d_parameters()
_ADJOINT_SP_STATE_DEPENDENT = _adjoint_state_dependent_1d_parameters()
_ADJOINT_SP_FREE_COVARIANCE = _adjoint_free_covariance_2d_parameters()
_ADJOINT_SP_TD_TI = _adjoint_td_ti_parameters()
_ADJOINT_SP_TD_CROSS_2D = _adjoint_td_cross_2d_parameters()
_ADJOINT_SP_COMPOSED_STATIC = _adjoint_composed_static_1d_parameters()

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
    values = [0.3, 0.2, -0.1, 0.4, 0.15, -0.2, 0.25, 0.1, -0.15, 0.35, 0.2]
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

@testset "Forward-gradient validation: TD predictor with off-diagonal Jtd (2-state)" begin
    # F6: `Jtd[2,1]` free at a 2-latent dimension, so `_reverse_td!`'s
    # `Jtd_bar` line is checked where `P` is not scalar and the operand order
    # actually matters.
    sp = _ADJOINT_SP_TD_CROSS_2D
    data = reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2], 2, :)
    tdpreds = reshape([0.0, 1.0, 0.0], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2], data, tdpreds)
    values = [0.2, -0.3, 0.15]
    _check_adjoint(objective, values)
end

@testset "Forward-gradient validation: T0 cells composed from two parameters" begin
    sp = _ADJOINT_SP_COMPOSED_STATIC
    data = reshape([0.1, -0.2, 0.15, 0.05], 1, :)
    objective = ContinuousTimeSEM.ctsem_objective(sp, [1], [0.0, 0.5, 1.2, 2.1], data)
    values = [0.3, -0.4, 0.25, -0.15]
    _check_adjoint(objective, values)

    # The support is discovered, not declared. Asserted directly as well as
    # through the gradient, because a support that silently lost `4` would
    # still pass every gate above for a model where `param[4]` happens to be
    # weakly identified -- and this fixture exists precisely to stop that.
    supports = ContinuousTimeSEM._ctsem_regular_transform_supports(sp, length(values))
    composed = [s for s in supports if length(s) > 1]
    @test length(composed) == 2
    @test all(s -> s == [3, 4], composed)

    # A transform that does not read its own parameter number is still a
    # fault, and still throws -- that is the half of the old one-parameter
    # assertion worth keeping, and the pullback would otherwise push this
    # cell's cotangent onto a parameter that never materialised it. Built as
    # its own tiny spec because the transforms are grouped by type into the
    # spec's own type parameter, so the closures of an existing spec cannot be
    # swapped out.
    misrendered = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            # Declares parameter 2 and reads parameter 1.
            (:DIFFUSION, 1, 1) => (2, "log1p_exp(param[1])"),
        ),
    )
    @test_throws ArgumentError ContinuousTimeSEM._ctsem_regular_transform_supports(
        ekf_from_data_frame(misrendered), 2)
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
    values = [0.3, 0.2, -0.1, 0.4, 0.15, -0.2, 0.25, 0.1, -0.15, 0.35, 0.2]

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
