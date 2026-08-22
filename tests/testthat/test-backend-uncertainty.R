# ctOptimUncertainty() for backend='julia' and backend='cpp' fits.
#
# The point of these tests is that the backends do not get their *own*
# uncertainty machinery: they build a log-probability/gradient function and hand
# it to the same `ctOptimComputeUncertainty()` the Stan path uses. So the thing
# worth testing is that the numbers come out the same as Stan's on a model where
# they are well defined -- which also transitively checks that the engines'
# gradients are right in a way the fixed-point parity tests do not, since a
# finite-difference Hessian probes a whole neighbourhood of the optimum.
#
# The model below is deliberately small and *fully identified*: one latent OU
# process with a non-individually-varying manifest mean. That last part matters.
# ctsem's default makes MANIFESTMEANS individually varying, which adds a
# population-SD parameter; with data simulated without individual differences
# that parameter sits on its boundary, the likelihood is flat in it, and its
# standard error is arbitrary in every backend. Comparing it would be comparing
# noise, so the model avoids creating it rather than the test tolerating it.

.backend_uncertainty_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8; diffusion <- 0.5
  do.call(rbind, lapply(seq_len(60), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
}

.backend_uncertainty_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1),
    MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1),
    CINT = matrix(0, 1, 1)))
}

test_that("C++ fits get Hessian uncertainty matching Stan's", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()

  cpp_fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE), cores = 1, verbose = 0))

  # The two optimizers must have landed in the same place, or the Hessians are
  # not comparable and a failure below would say nothing about uncertainty.
  expect_equal(cpp_fit$estimate$raw, stan_fit$stanfit$rawest, tolerance = 1e-3)

  cpp_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(cpp_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))
  stan_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(stan_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))

  cpp_se <- sqrt(diag(cpp_unc$estimate$cov))
  stan_se <- sqrt(diag(stan_unc$stanfit$cov))
  expect_equal(cpp_se, stan_se, tolerance = 1e-3)

  # The fit carries usable uncertainty afterwards, not just a covariance.
  expect_equal(dim(cpp_unc$estimate$rawposterior), c(200L, length(cpp_se)))
  expect_equal(cpp_unc$estimate$se, cpp_se)
  expect_identical(cpp_unc$uncertainty$settings$method, "hessian")
  summarised <- summary(cpp_unc)
  expect_equal(dim(summarised$ci), c(length(cpp_se), 2L))
  expect_false(is.null(summarised$se))
})

test_that("Julia fits get Hessian uncertainty matching Stan's", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("JuliaConnectoR")
  skip_if(!isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
    "Julia is not available.")
  skip_on_cran()

  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()

  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE), cores = 1, verbose = 0))
  expect_equal(julia_fit$estimate$raw, stan_fit$stanfit$rawest, tolerance = 1e-3)

  julia_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(julia_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))
  stan_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(stan_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))

  expect_equal(sqrt(diag(julia_unc$estimate$cov)), sqrt(diag(stan_unc$stanfit$cov)),
    tolerance = 1e-3)
})

test_that("unsupported uncertainty methods are refused by name, not silently", {
  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()[1:24, ]
  cpp_fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))

  # The score-based and full-bootstrap methods need per-subject scores or
  # refits, which these engines do not expose. Refusing by name beats producing
  # a plausible-looking covariance from a method that did not actually run.
  for (method in c("opg", "sandwich", "bootstrap", "fullbootstrap")) {
    expect_error(ctOptimUncertainty(cpp_fit, uncertainty = method),
      "not available for backend", fixed = FALSE)
  }
})
