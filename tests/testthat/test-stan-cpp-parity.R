# Stan vs. the hand-written C++ backend.
#
# Deliberately the same model shapes, data and fixed raw-parameter points as
# `test-stan-julia-parity.R`, because those shapes were chosen by the bugs they
# caught: partial (not total) row missingness, a nonlinear LAMBDA combined with
# missingness, the two structurally different kinds of random effect that
# `intoverpopindvaryingindex` conflates, and a T0MEANS-indvarying population SD
# with a non-unit meanscale. A third numerical backend has exactly the same
# opportunities to drift away from Stan as the second one did, so it is checked
# against exactly the same cases rather than against fresh, cleaner ones.
#
# Unlike the Julia suite this needs no external toolchain and no environment
# variable: the C++ engine is compiled into ctsem's own shared library, so these
# tests run wherever ctsem itself is installed. Only rstan is required, for the
# ground truth.

# Reuse ctsem's own precompiled generic Stan binary when `standata$recompile`
# says it is safe to, exactly the way ctFit's Stan backend decides this. When a
# real compile is unavoidable it is cached by generated-code hash, since several
# tests here share a model definition.
.cpp_parity_stan_cache <- new.env(parent = emptyenv())
.cpp_compiled_stan_fit <- function(stan_spec) {
  if (identical(stan_spec$standata$recompile, 0L)) {
    return(ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, stan_spec$standata))
  }
  key <- digest::digest(stan_spec$stanmodeltext)
  if (!exists(key, envir = .cpp_parity_stan_cache, inherits = FALSE)) {
    assign(key, rstan::stan_model(model_code = stan_spec$stanmodeltext),
      envir = .cpp_parity_stan_cache)
  }
  ctsem:::stan_reinitsf(get(key, envir = .cpp_parity_stan_cache, inherits = FALSE),
    stan_spec$standata)
}

.cpp_parity_skip <- function() {
  skip_if_not_installed("rstan")
  skip_if_not_installed("digest")
}

test_that("Stan and C++ agree for a linear likelihood", {
  .cpp_parity_skip()

  model <- ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- -.7

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for a linear augmented random effect", {
  .cpp_parity_skip()

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, .1, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for a row with partial (not total) missingness", {
  .cpp_parity_skip()

  # The branch where some, but not all, manifest variables are observed in one
  # row. In the Julia port this branch carried a double-applied S^-1 in the
  # Kalman-gain state update for months, because every other test's rows were
  # either fully observed or fully missing.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1), Y2 = c(0, -.1, .1, .2, NA, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for a fully missing row", {
  .cpp_parity_skip()

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2), PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, NA, -.1), Y2 = c(0, -.1, .1, .2, NA, 0))
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- c(-1, .1, -1, .2, -.4)

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for nonlinear predictors and augmented states", {
  .cpp_parity_skip()

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
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  raw <- c(-1, .1, -1, .1, .2, -.4, .05)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  # Repeated nonlinear EKF linearisation accumulates a sub-micro likelihood
  # difference across Stan's and this engine's otherwise equivalent kernels.
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 2e-6)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for a T0MEANS-indvarying population SD with non-unit meanscale", {
  .cpp_parity_skip()

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
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- c(0.4112875, -0.1694095, 0.1089385)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("Stan and C++ agree for a model mixing both kinds of random effect", {
  .cpp_parity_skip()

  # Five augmented states: two original, still-dynamic states with directly
  # indvarying T0MEANS, and three static carrier states for the non-T0MEANS
  # random effects. Five is one more than the threshold at which the Julia
  # backend switches Lyapunov solvers, which is how a wrong
  # `dynamic_state_indices` first surfaced there as a LAPACK exception.
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
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  expect_equal(cpp_spec$nlatent_augmented, 5L)
  expect_equal(cpp_spec$dynamic_state_indices, 1:2)

  raw <- c(0.6861741, -0.3590315, -0.2082878, -0.1236879, -0.291202, -0.284184,
    0.2244418, -0.03508657, 0.04579729, 0.6569934, 0.1070959, 0.8150255,
    0.6844356, 0.09720616, 0.5688201, 0.1403042, -0.2681402, -0.09219849,
    -0.001446727, 0.2964492, 0.2519251, 0.2116025)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-6)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-6)
})

test_that("Stan and C++ agree for a state/TD-dependent measurement equation with partial missingness", {
  .cpp_parity_skip()

  # A nonlinear LAMBDA combined with missingness: the update must use LAMBDA
  # for the innovation but its Jacobian (Jy) for covariance propagation. The
  # two coincide for a fixed or linear LAMBDA, so nothing else here can tell
  # them apart.
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
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  raw <- c(-0.1773093, 0.007978311, -0.4549659, -0.408796, 0.3535467, -0.2802454)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 2e-5)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-5)
})

test_that("Stan and C++ agree for 3 original (not just augmented) latents", {
  .cpp_parity_skip()

  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 3,
    LAMBDA = matrix(c(1, 0, 0, 0, 1, 0, 0, "cross32", 1), 3, 3, byrow = TRUE),
    DRIFT = matrix(c("d11", "d12", "d13", "d21", "d22", "d23", "d31", "d32", "d33"),
      3, 3, byrow = TRUE),
    DIFFUSION = matrix(c("diff11", 0, 0, "diff21", "diff22", 0, "diff31", "diff32", "diff33"),
      3, 3, byrow = TRUE),
    MANIFESTVAR = diag(c(.1, .1, .1)), MANIFESTMEANS = matrix(0, 3, 1),
    T0VAR = matrix(c("t0v11", 0, 0, "t0v21", "t0v22", 0, "t0v31", "t0v32", "t0v33"),
      3, 3, byrow = TRUE),
    T0MEANS = matrix(0, 3, 1)
  ))
  data <- data.frame(
    id = rep(1:2, each = 4), time = rep(c(0, .5, 1, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1, -.05, .05),
    Y2 = c(0, -.1, .1, .2, .1, 0, .1, -.05),
    Y3 = c(0, .05, -.05, .1, -.1, .05, 0, .1)
  )
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE))
  stan_fit <- .cpp_compiled_stan_fit(stan_spec)
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE))
  raw <- c(-0.28858, -0.08775772, 0.07763646, -0.3456396, 0.05873485, 0.009037183,
    0.02562532, 0.3349831, -0.3656572, 0.3802106, -0.2234345, -0.3393656,
    -0.2149075, 0.07579571, 0.04561371, -0.09229693, -0.2859052, -0.1944728,
    0.3672941, 0.05994348, -0.1735451, -0.2826902)
  expect_equal(rstan::get_num_upars(stan_fit), length(raw))

  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE, gradient = TRUE)
  cpp_value <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(cpp_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(cpp_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

test_that("The C++ adjoint gradient matches central finite differences of its own likelihood", {
  # The independent referee. Stan agreement checks the whole engine against
  # another implementation; this checks the reverse pass against the forward
  # pass it is supposed to mirror, which is the failure mode a hand-written
  # adjoint actually has -- a dropped term that happens to be small wherever
  # Stan agreement is loose. No rstan required.
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2),
    PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c(
      "d11", "PARS[1,1] * (1 + .05 * eta1 + .02 * dose)",
      "d21", "d22"
    ), 2, 2, byrow = TRUE),
    DIFFUSION = matrix(c("diff11", 0, "diff21", "diff22"), 2, 2, byrow = TRUE),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(c("mm1", "mm2"), 2, 1),
    T0VAR = matrix(c("t0v11", 0, "t0v21", "t0v22"), 2, 2, byrow = TRUE),
    T0MEANS = matrix(c("t0m1", "t0m2"), 2, 1),
    n.TDpred = 1, TDpredNames = "dose", TDPREDEFFECT = matrix(c("impulse", 0), 2, 1),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "d11"] <- TRUE
  set.seed(11)
  data <- data.frame(
    id = rep(1:3, each = 4), time = rep(c(0, .4, 1.1, 2.0), 3),
    Y1 = rnorm(12, 0, .5), Y2 = rnorm(12, 0, .5),
    dose = rep(c(0, 1, 0, 1), 3), group = rep(rnorm(3), each = 4)
  )
  data$Y2[5] <- NA   # partial row
  data$Y1[10] <- NA; data$Y2[10] <- NA  # fully missing row

  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE,
    priors = FALSE, nlcontrol = list(maxtimestep = .25)))
  npar <- max(c(cpp_spec$parameter_table$parnumber, cpp_spec$ti_effects$coefficient),
    na.rm = TRUE)
  set.seed(5)
  raw <- rnorm(npar, 0, .25)

  analytic <- ctCppEvaluate(cpp_spec, raw, gradient = TRUE)
  step <- 1e-6
  finite <- vapply(seq_len(npar), function(i) {
    up <- raw; down <- raw
    up[i] <- up[i] + step; down[i] <- down[i] - step
    (ctCppEvaluate(cpp_spec, up, gradient = FALSE)$value -
       ctCppEvaluate(cpp_spec, down, gradient = FALSE)$value) / (2 * step)
  }, numeric(1))

  expect_true(all(is.finite(analytic$gradient)))
  expect_lt(max(abs(analytic$gradient - finite)) / max(1, max(abs(finite))), 1e-6)
})

test_that("The C++ backend's own optimizer recovers the same optimum as Stan's", {
  .cpp_parity_skip()

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

  cf <- suppressMessages(ctFit(data, model = model, backend = "cpp", verbose = 0))
  sf <- suppressMessages(ctFit(data, model = model, backend = "stan",
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE),
    optimize = TRUE, verbose = 0, savescores = FALSE, cores = 1))

  expect_equal(cf$estimate$loglik, -sf$stanfit$optimfit$f, tolerance = 1e-3)
  expect_equal(cf$estimate$raw, sf$stanfit$rawest, tolerance = 1e-2)
})
