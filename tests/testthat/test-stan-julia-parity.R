# Reuse ctsem's own precompiled generic Stan binary (`stanmodels$ctsm`) when
# `standata$recompile` says it's safe to, exactly the way `ctFit`'s own Stan
# backend decides this (see `stanoptimis.R`: `if(standata$recompile==0)
# smf <- stan_reinitsf(stanmodels$ctsm, standata)`). `recompile` is set in
# `ctFit.R` based on genuine structural requirements -- e.g. any state- or
# TD/TI-dependent calc not already expressible via `JAx[...]` (`ncalcsNoJ`,
# which includes the PARS-substituted-into-DRIFT and nonlinear-DRIFT models
# below) forces it to 1, because those need an analytic Jacobian the generic
# binary's finite-difference fallback doesn't have. An earlier version of
# this file called `stan_reinitsf(stanmodels$ctsm, ...)` unconditionally,
# which happened to work for the simpler models here but silently produced
# `multi_normal_cholesky_lpdf: Location parameter is nan`/`quad_form_sym: A
# is not symmetric, Inf` for the `recompile==1` ones -- checking the flag
# `ctFit` itself computes, instead of reimplementing (badly) a guess at when
# reuse is safe, is what actually fixes that.
#
# When `recompile==1`, several tests/calls here happen to reuse the exact
# same model definition (only `standata` differs -- e.g. test 4 builds two
# specs from one `model` object, and tests 2/3 share a model entirely), so
# the actual compile is additionally cached by generated-code hash to avoid
# paying for it more than once per distinct model.
.parity_stan_cache <- new.env(parent = emptyenv())

# Tiered, because the two tiers differ by three orders of magnitude in cost.
#
# A model shape with `recompile == 0` reuses ctsem's own precompiled binary and
# costs seconds. A shape with `recompile == 1` -- the nonlinear and
# state-dependent ones -- calls rstan::stan_model() and compiles C++ at test
# time: 15 to 25 minutes EACH on Windows, per distinct generated model. That is
# not a slow machine, it is a test file that compiles compilers' worth of code,
# and no amount of hardware fixes it.
#
# So the compiling ones skip by default and run on request:
#
#   CTSEM_PARITY_COMPILE=true Rscript -e '...test_file(...)'
#
# Run them before a release, and whenever anything touches the generated stan
# program, since that is what they cover. The rest run in any ordinary suite,
# which is the point: parity that never runs protects nothing, and that was
# this file's previous state for a different reason (it skipped on an
# environment variable nothing ever set).

.compiled_stan_fit <- function(stan_spec) {
  if (identical(stan_spec$standata$recompile, 0L)) {
    return(ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, stan_spec$standata))
  }
  testthat::skip_if(
    !identical(Sys.getenv("CTSEM_PARITY_COMPILE"), "true"),
    paste("this model shape compiles a fresh stan program (minutes);",
      "set CTSEM_PARITY_COMPILE=true to run it"))
  key <- digest::digest(stan_spec$stanmodeltext)
  if (!exists(key, envir = .parity_stan_cache, inherits = FALSE)) {
    assign(key, rstan::stan_model(model_code = stan_spec$stanmodeltext), envir = .parity_stan_cache)
  }
  ctsem:::stan_reinitsf(get(key, envir = .parity_stan_cache, inherits = FALSE), stan_spec$standata)
}

test_that("Stan and Julia agree for a linear likelihood", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  model <- ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  raw <- -.7

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for a linear augmented random effect", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, .1, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for a discrete-time model with augmented random effects", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # Every other test in this file is `type = "ct"`, and that is how a
  # discrete-time defect survived: the engine's one-step form
  # (`_compute_one_step_form!`) formed the local affine offset over every
  # state instead of only the diffusing ones. In continuous time that is
  # harmless, because DRIFT and JAx agree (both zero) on the static
  # coordinates a random effect augments the state with. In discrete time
  # `ctJacobian()` puts 1 on JAx's augmented diagonal -- a static state
  # carries forward one whole step -- while the DRIFT the engine holds stays
  # padded with zeros there, so the offset came out as `-x[i]` and cancelled
  # the state: every random-effect coordinate was zeroed at every step, and a
  # subject's intercept deviation applied to the first transition only.
  #
  # On the tutorial handbook's two-variable ESM model that moved the
  # autoregression from .71 to .99 and the intercept from -.01 to -.55, with
  # nothing failing: a near-unit root absorbing the individual differences in
  # level the augmented intercept could no longer carry.
  #
  # T0MEANS and CINT are `indvarying` by default, so the plain specification
  # below is the one a user writes -- 13 population means, 4 population SDs
  # and 6 population correlations.
  model <- suppressMessages(ctModel(
    type = "dt", n.latent = 2, LAMBDA = diag(2),
    DRIFT = matrix(c("d11", "d12", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c("diff11", 0, "diff21", "diff22"), 2, 2, byrow = TRUE),
    MANIFESTVAR = matrix(c("merr1", 0, 0, "merr2"), 2, 2, byrow = TRUE),
    MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(c("cint1", "cint2"), 2, 1),
    T0MEANS = matrix(c("t0m1", "t0m2"), 2, 1)))
  set.seed(11)
  data <- data.frame(id = rep(1:3, each = 5),
    time = rep(c(0, .125, .25, .5, 1), 3),
    Y1 = rnorm(15), Y2 = rnorm(15))
  data$Y2[4] <- NA

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))

  # The augmentation is what this test is about, so assert it happened rather
  # than silently checking parity on a two-state model.
  expect_equal(julia_spec$nlatent_augmented, 4L)
  npar <- rstan::get_num_upars(stan_fit)
  expect_equal(npar, 23L)

  set.seed(5)
  raw <- rnorm(npar, 0, .3)
  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  # Both gradient methods: the reverse pass has its own copy of the one-step
  # form (`_reverse_predict_discrete!`) and had the same defect.
  forward <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "forward")
  adjoint <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "adjoint")

  expect_equal(as.numeric(forward$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(forward$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
  expect_equal(as.numeric(adjoint$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and Julia agree for a latent with no diffusion of its own", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # Two sub-blocks of the state appear in the continuous form and conflating
  # them was a bug:
  #
  #   derrind    the states with their own diffusion, grown through JAx
  #              coupling (`.ctJuliaDerrind`). The Lyapunov solve lives here.
  #   1:nlatent  everything that is not a static random-effect carrier. The
  #              intercept solve and the affine offset live here.
  #
  # They coincide for every other model in this file, because a full DIFFUSION
  # makes every latent diffusion-reachable -- which is why using the first for
  # the second passed the whole suite. Latent 2 below has a zero DIFFUSION row
  # and column and no coupling to latent 1, so derrind excludes it, and the
  # affine offset is NOT zero on its row: `cint2` left the likelihood
  # altogether, with a gradient of exactly zero, so the optimiser held it at
  # its starting value and the summary reported that as an estimate. 1745 log
  # posterior units on this model.
  #
  # `drift22` came out wrong too, not just `cint2`: the predicted mean for the
  # indicator loading on latent 2 was wrong, so everything touching that state
  # was. Both are checked.
  #
  # Run in both time bases. The discrete form has no solve and so needs no
  # index set at all -- it runs over every row -- but it was briefly given the
  # same wrong one, so it is pinned here as well.
  for (type in c("ct", "dt")) {
    model <- suppressWarnings(suppressMessages(ctModel(
      type = type, n.latent = 2, LAMBDA = diag(2),
      DRIFT = matrix(c("drift11", 0, 0, "drift22"), 2, 2, byrow = TRUE),
      DIFFUSION = matrix(c("diff11", 0, 0, 0), 2, 2, byrow = TRUE),
      MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
      T0VAR = diag(2), T0MEANS = matrix(c("t0m1", "t0m2"), 2, 1),
      CINT = matrix(c("cint1", "cint2"), 2, 1))))
    model$pars$indvarying <- FALSE   # keep this about derrind, not augmentation
    set.seed(4)
    data <- data.frame(id = rep(1:3, each = 4), time = rep(c(0, 1, 2, 3), 3),
      Y1 = stats::rnorm(12, 0, .5), Y2 = stats::rnorm(12, 0, .5))

    stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
      priors = FALSE))
    stan_fit <- .compiled_stan_fit(stan_spec)
    julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
      priors = FALSE))

    # The premise of the test: derrind really is the smaller set here. Without
    # this the model could quietly become an ordinary one and prove nothing.
    expect_equal(julia_spec$dynamic_state_indices, 1L, info = type)
    expect_equal(julia_spec$nlatent, 2L, info = type)

    npar <- rstan::get_num_upars(stan_fit)
    set.seed(9)
    raw <- stats::rnorm(npar, 0, .3)
    stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE,
      gradient = TRUE)
    stan_grad <- as.numeric(attributes(stan_value)$gradient)
    forward <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "forward")
    adjoint <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "adjoint")

    # Relative, because this raw vector puts the model somewhere stiff: the
    # gradient entries run to 5e4, so 1e-7 relative is ~5e-3 absolute and is
    # accumulated floating point, not a difference in the model. Pre-fix the gap
    # was 1745 on the value and the whole of one gradient entry.
    expect_equal(as.numeric(forward$value), as.numeric(stan_value),
      tolerance = 1e-7, info = type)
    expect_equal(as.numeric(forward$gradient), stan_grad, tolerance = 1e-7, info = type)
    expect_equal(as.numeric(adjoint$gradient), stan_grad, tolerance = 1e-7, info = type)

    # The two entries the defect actually moved, named so a failure says which.
    which_cint2 <- which(julia_spec$parameter_table$param %in% "cint2")[1L]
    par_cint2 <- julia_spec$parameter_table$parnumber[which_cint2]
    expect_false(isTRUE(all.equal(as.numeric(adjoint$gradient)[par_cint2], 0)),
      info = type)
  }
})

test_that("Stan and Julia agree for a row with partial (not total) missingness", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # This exercises _ekf_update_observed!()'s partial-observation branch (some,
  # but not all, manifest variables observed in a row) -- the branch that
  # previously double-applied S^{-1} in the Kalman gain state update
  # (`gain * (factor \ innovation)` instead of `gain * innovation`), which the
  # other parity tests here never reached because their rows are always fully
  # observed or fully missing.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, NA, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for nonlinear predictors and augmented states", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2),
    PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c(
      "d11", "PARS[1,1] * (1 + .05 * eta1 + .02 * dose)",
      "d21", "d22"
    ), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1),
    n.TDpred = 1, TDpredNames = "dose", TDPREDEFFECT = matrix(c("impulse", 0), 2, 1),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "d11"] <- TRUE
  data <- data.frame(
    id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, .1, 0),
    dose = c(0, 1, 0, 1, 0, 1), group = rep(c(-.5, .75), each = 3)
  )

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  stan_fit <- .compiled_stan_fit(stan_spec)

  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  raw <- c(-1, .1, -1, .1, .2, -.4, .05)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))
  expect_equal(max(c(julia_spec$parameter_table$parnumber,
    julia_spec$ti_effects$coefficient), na.rm = TRUE), length(raw))

  zero_predictors <- data
  zero_predictors$dose <- 0
  zero_predictors$group <- 0
  zero_stan_spec <- suppressMessages(ctFit(zero_predictors, model, backend = "stan", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  zero_stan_fit <- .compiled_stan_fit(zero_stan_spec)
  zero_julia_spec <- suppressMessages(ctFit(zero_predictors, model, backend = "julia", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  zero_stan_value <- rstan::log_prob(zero_stan_fit, upars = raw,
    adjust_transform = FALSE, gradient = TRUE)
  zero_julia_value <- ctJuliaEvaluate(zero_julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(zero_julia_value$value), as.numeric(zero_stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(zero_julia_value$gradient),
    as.numeric(attributes(zero_stan_value)$gradient), tolerance = 1e-7)

  stan_value <- rstan::log_prob(stan_fit, upars = raw,
    adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)

  # Repeated nonlinear EKF linearisation accumulates a sub-micro likelihood
  # difference across Stan and Julia's otherwise equivalent matrix kernels.
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 2e-6)
  expect_equal(as.numeric(julia_value$gradient),
    as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for a T0MEANS-indvarying population SD with non-unit meanscale", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # Stan builds the population covariance in raw-parameter units, then
  # explicitly rescales each indvarying T0MEANS row/column by that parameter's
  # `multiplier*meanscale` to convert to state-space units before it's used as
  # a covariance (ctModelWriter.R: `T0cov[matsetup[ri,1], ] *= matvalues[ri,2]
  # * matvalues[ri,3]` and the matching column update). T0MEANS/CINT-type
  # custom pars default to meanscale=10 (`.ctModelDefaultFreePar`), so any
  # T0MEANS-indvarying parameter -- like `t0m||TRUE` below -- exercises this.
  # The Julia port originally omitted this rescaling entirely, so its fitted
  # population SD came out ~10x too large relative to Stan's (same maximum
  # likelihood, different raw parameter, since the two backends' unconstrained
  # spaces disagreed) -- this is the model shape that surfaced it, distilled
  # from ctsemTutorial.qmd's individual-differences example.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 1, LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix("t0m||TRUE", 1, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  raw <- c(0.4112875, -0.1694095, 0.1089385)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for a moderate-dimensional model mixing both kinds of random effect", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # Distilled from ctsemTutorial.qmd's individual-differences example, which
  # originally crashed Julia outright with a LAPACKException from the
  # Schur-based Lyapunov solver (`_lyap_solve_factorized!`/`trsyl!`) before
  # `dynamic_state_indices` was fixed in `.ctJuliaAugmentRandomEffects`
  # (ctsem/R/ctJuliaBackend.R), and separately mismatched Stan's T0VAR
  # population SD before the `multiplier*meanscale` rescaling above was
  # added. This model combines *both* kinds of population-varying state that
  # `augmented_indices`/`intoverpopindvaryingindex` conflates: two original,
  # still-dynamic states with directly indvarying T0MEANS (`t0a`, `t0b`) and
  # two newly-created static carrier states for non-T0MEANS random effects
  # (`cross21` on DRIFT, `B1`/`B2` on CINT) -- five augmented states total,
  # one more than the `n > 4` threshold that selects the Schur solver over
  # the small-system `ksolve!` path. If either fix regresses, this either
  # throws (dynamic_state_indices) or the value/gradient checks below fail
  # (meanscale rescaling).
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(1, 2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0VAR = diag(2),
    T0MEANS = c("t0a||TRUE", "t0b||TRUE"),
    CINT = c("B1||TRUE", "B2||TRUE"),
    DRIFT = matrix(c("auto1", "cross21||TRUE", "cross21||TRUE", "auto2"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15))
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, .1, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  expect_equal(julia_spec$nlatent_augmented, 5L)
  expect_equal(julia_spec$dynamic_state_indices, 1:2)

  raw <- c(0.6861741, -0.3590315, -0.2082878, -0.1236879, -0.291202, -0.284184,
    0.2244418, -0.03508657, 0.04579729, 0.6569934, 0.1070959, 0.8150255,
    0.6844356, 0.09720616, 0.5688201, 0.1403042, -0.2681402, -0.09219849,
    -0.001446727, 0.2964492, 0.2519251, 0.2116025)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-6)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-6)
})

test_that("Stan and Julia agree for a state/TD-dependent measurement equation with partial missingness", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # `_ekf_masked_update_step!`'s partial-observation branch must use LAMBDA
  # for the innovation/mean but its Jacobian (pars.Jy) for covariance
  # propagation -- these coincide for a fixed/linear LAMBDA (every other
  # test here, and the model that motivated this fix, use one), so a bug
  # specific to *nonlinear* LAMBDA combined with missingness would not have
  # been caught anywhere else. This was flagged as an untested combination
  # in HANDOFF.md rather than a known bug; this test confirms it's actually
  # correct (Stan and Julia agree here) rather than leaving it unverified.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = matrix(c(
      "1", "0",
      "1 + .1 * eta1 + .05 * dose", "1"
    ), 2, 2, byrow = TRUE),
    PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1),
    n.TDpred = 1, TDpredNames = "dose", TDPREDEFFECT = matrix(c("impulse", 0), 2, 1)
  ))
  data <- data.frame(
    id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, NA, 0),
    dose = c(0, 1, 0, 1, 0, 1)
  )
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  raw <- c(-0.1773093, 0.007978311, -0.4549659, -0.408796, 0.3535467, -0.2802454)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  # Same sub-micro nonlinear-EKF-linearisation tolerance as the other
  # nonlinear-model test above.
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 2e-5)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-5)
})

test_that("Stan and Julia agree for 3 original (not just augmented) latents", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # The "moderate-dimensional" test above has 5 *augmented* states but only
  # 2 *original* latents -- this checks the core EKF math (matrix-exponential
  # discretization, Lyapunov solve, Kalman recursion) directly at 3+ original
  # latents, fully cross-coupled DRIFT/DIFFUSION/T0VAR, with no indvarying
  # parameters or state-dependent calcs so this isolates dimensionality from
  # the random-effects machinery already covered elsewhere.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 3,
    LAMBDA = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, "cross32", 1
    ), 3, 3, byrow = TRUE),
    DRIFT = matrix(c(
      "d11", "d12", "d13",
      "d21", "d22", "d23",
      "d31", "d32", "d33"
    ), 3, 3, byrow = TRUE),
    DIFFUSION = matrix(c(
      "diff11", 0, 0,
      "diff21", "diff22", 0,
      "diff31", "diff32", "diff33"
    ), 3, 3, byrow = TRUE),
    MANIFESTVAR = diag(c(.1, .1, .1)),
    MANIFESTMEANS = matrix(0, 3, 1),
    T0VAR = matrix(c(
      "t0v11", 0, 0,
      "t0v21", "t0v22", 0,
      "t0v31", "t0v32", "t0v33"
    ), 3, 3, byrow = TRUE),
    T0MEANS = matrix(0, 3, 1)
  ))
  data <- data.frame(
    id = rep(1:2, each = 4), time = rep(c(0, .5, 1, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1, -.05, .05),
    Y2 = c(0, -.1, .1, .2, .1, 0, .1, -.05),
    Y3 = c(0, .05, -.05, .1, -.1, .05, 0, .1)
  )
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  raw <- c(-0.28858, -0.08775772, 0.07763646, -0.3456396, 0.05873485, 0.009037183,
    0.02562532, 0.3349831, -0.3656572, 0.3802106, -0.2234345, -0.3393656,
    -0.2149075, 0.07579571, 0.04561371, -0.09229693, -0.2859052, -0.1944728,
    0.3672941, 0.05994348, -0.1735451, -0.2826902)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia's actual optimizers converge to the same fit for TD/TI + individual differences", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # Every other test here checks log_prob/gradient agreement at one fixed
  # raw-parameter point -- necessary but not sufficient, since that's exactly
  # what the T0-SD meanscale bug could still pass (both backends reached the
  # *same maximum likelihood* from *different* raw parameters; a fixed-point
  # check at either backend's own optimum wouldn't by itself reveal that
  # unless you specifically evaluated at the *other* backend's point, as the
  # investigation that found it had to do). This test instead runs each
  # backend's own real optimizer (ctsem_optimize's Optim.LBFGS for Julia,
  # Stan's L-BFGS) on the combined TD/TI-predictor + both-kinds-of-
  # individual-differences model from the "moderate-dimensional" test above,
  # and checks that they land on matching loglik *and* matching raw
  # parameters -- the actual end-to-end guarantee a fixed point can't give.
  # Data is deliberately minimal (6 subjects, 4 waves) to keep this fast;
  # this model shape has standata$recompile==0, so no C++ compile is needed.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(1, 2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0VAR = diag(2),
    T0MEANS = c("t0a||TRUE", "t0b||TRUE"),
    CINT = c("B1||TRUE", "B2||TRUE"),
    DRIFT = matrix(c("auto1", "cross21||TRUE", "cross21||TRUE", "auto2"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)),
    n.TDpred = 1, TDpredNames = "dose", TDPREDEFFECT = matrix(c("impulse", 0), 2, 1),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "B1"] <- TRUE

  set.seed(21)
  NSubjects <- 6
  times <- c(0, .5, 1, 1.5)
  data <- data.frame()
  for (i in 1:NSubjects) {
    data <- rbind(data, data.frame(
      id = i, time = times,
      Y1 = rnorm(length(times), 0, .5), Y2 = rnorm(length(times), 0, .5),
      dose = c(0, 1, 0, 1),
      group = rep(rnorm(1), length(times))
    ))
  }

  jf <- suppressMessages(ctFit(data, model = model, backend = "julia",
    verbose = 0))
  sf <- suppressMessages(ctFit(data, model = model, backend = "stan",
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE),
    optimize = TRUE, verbose = 0, savescores = FALSE, cores = 1))

  expect_equal(jf$estimate$loglik, -sf$stanfit$optimfit$f, tolerance = 1e-3)
  expect_equal(jf$estimate$raw, sf$stanfit$rawest, tolerance = 1e-2)
})

test_that("Julia's adjoint gradient matches its forward gradient and Stan", {
  skip_without_julia()
  skip_if_not_installed("rstan")

  # The Julia-side suite already checks the adjoint against ForwardDiff and
  # finite differences across every model shape
  # (`test_adjoint_gradient_validation.jl`). What this adds is the part only
  # the R side can check: that `gradient_method` actually survives the
  # JuliaConnectoR boundary (R marshals character vectors to Julia `String`,
  # not `Symbol`), and that the adjoint agrees with *Stan* -- the independent
  # ground truth -- and not merely with the other Julia gradient, which shares
  # the same forward filter and so could in principle share a misreading of it.
  model <- ctModel(type = "ct", LAMBDA = diag(2),
    DRIFT = matrix(c("drift11", "drift21", "drift12", "drift22"), 2, 2),
    DIFFUSION = matrix(c("diff11", "diff21", 0, "diff22"), 2, 2),
    MANIFESTVAR = matrix(c(.1, 0, 0, .1), 2, 2),
    MANIFESTMEANS = matrix(c("mm1", "mm2"), 2, 1),
    T0VAR = matrix(c(1, 0, 0, 1), 2, 2), T0MEANS = matrix(c("t0m1", "t0m2"), 2, 1),
    CINT = matrix(c("cint1", 0), 2, 1))
  set.seed(3)
  data <- data.frame(id = rep(1:4, each = 4), time = rep(c(0, .4, 1.1, 2.0), 4),
    Y1 = rnorm(16, 0, .5), Y2 = rnorm(16, 0, .5))

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .compiled_stan_fit(stan_spec)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))

  npar <- max(c(julia_spec$parameter_table$parnumber,
    julia_spec$ti_effects$coefficient), na.rm = TRUE)
  set.seed(7)
  raw <- rnorm(npar, 0, .3)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  forward <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "forward")
  adjoint <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE, gradient_method = "adjoint")

  expect_equal(as.numeric(adjoint$value), as.numeric(forward$value), tolerance = 1e-12)
  expect_equal(as.numeric(adjoint$gradient), as.numeric(forward$gradient), tolerance = 1e-9)
  expect_equal(as.numeric(adjoint$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(adjoint$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})
