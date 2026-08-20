test_that("Stan and Julia agree for a linear likelihood", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("JuliaConnectoR")
  project <- Sys.getenv("CTSEM_JULIA_PROJECT", unset = "")
  skip_if(!nzchar(project) || !dir.exists(project),
    "Set CTSEM_JULIA_PROJECT to the local ContinuousTimeSEM project to run backend parity tests.")

  model <- ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_model <- rstan::stan_model(model_code = stan_spec$stanmodeltext)
  stan_fit <- ctsem:::stan_reinitsf(stan_model, stan_spec$standata)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE, backendcontrol = list(julia_project = project)))
  raw <- -.7

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for a linear augmented random effect", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("JuliaConnectoR")
  project <- Sys.getenv("CTSEM_JULIA_PROJECT", unset = "")
  skip_if(!nzchar(project) || !dir.exists(project),
    "Set CTSEM_JULIA_PROJECT to the local ContinuousTimeSEM project to run backend parity tests.")

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, .1, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_model <- rstan::stan_model(model_code = stan_spec$stanmodeltext)
  stan_fit <- ctsem:::stan_reinitsf(stan_model, stan_spec$standata)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE, backendcontrol = list(julia_project = project)))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient), tolerance = 1e-7)
})

test_that("Stan and Julia agree for nonlinear predictors and augmented states", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("JuliaConnectoR")
  project <- Sys.getenv("CTSEM_JULIA_PROJECT", unset = "")
  skip_if(!nzchar(project) || !dir.exists(project),
    "Set CTSEM_JULIA_PROJECT to the local ContinuousTimeSEM project to run backend parity tests.")

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
  stan_model <- rstan::stan_model(model_code = stan_spec$stanmodeltext)
  stan_fit <- ctsem:::stan_reinitsf(stan_model, stan_spec$standata)

  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25),
    backendcontrol = list(julia_project = project)))
  raw <- c(-1, .1, -1, .1, .2, -.4, .05)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))
  expect_equal(max(c(julia_spec$parameter_table$parnumber,
    julia_spec$ti_effects$coefficient), na.rm = TRUE), length(raw))

  zero_predictors <- data
  zero_predictors$dose <- 0
  zero_predictors$group <- 0
  zero_stan_spec <- suppressMessages(ctFit(zero_predictors, model, backend = "stan", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  zero_stan_fit <- ctsem:::stan_reinitsf(stan_model, zero_stan_spec$standata)
  zero_julia_spec <- suppressMessages(ctFit(zero_predictors, model, backend = "julia", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25),
    backendcontrol = list(julia_project = project)))
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
