# ctOptimUncertainty() for backend='julia' fits.
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

test_that("Julia fits get Hessian uncertainty matching Stan's", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  skip_without_julia()

  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()

  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE), cores = 1, verbose = 0))

  # The two optimizers must have landed in the same place, or the Hessians are
  # not comparable and a failure below would say nothing about uncertainty.
  expect_equal(julia_fit$estimate$raw, stan_fit$stanfit$rawest, tolerance = 1e-3)

  julia_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(julia_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))
  stan_unc <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(stan_fit, uncertainty = "hessian", finishsamples = 200, verbose = 0)))

  julia_se <- sqrt(diag(julia_unc$estimate$cov))
  stan_se <- sqrt(diag(stan_unc$stanfit$cov))
  expect_equal(julia_se, stan_se, tolerance = 1e-3)

  # The fit carries usable uncertainty afterwards, not just a covariance.
  expect_equal(dim(julia_unc$estimate$rawposterior), c(200L, length(julia_se)))
  expect_equal(julia_unc$estimate$se, julia_se)
  expect_identical(julia_unc$uncertainty$settings$method, "hessian")
  # The summary reports on the transformed scale (see test-backend-summary.R),
  # so the raw-scale standard errors this test compares against Stan live on the
  # fit rather than in the printed summary; what the summary must show is an
  # interval per free parameter, earned from the draws.
  summarised <- summary(julia_unc)
  expect_equal(nrow(summarised$popmeans), length(julia_se))
  expect_true(all(c("2.5%", "97.5%") %in% colnames(summarised$popmeans)))
})

test_that("the Hessian is exact, and agrees with the finite difference it replaces", {
  skip_on_cran()
  skip_without_julia()

  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  exact <- suppressWarnings(suppressMessages(ctOptimUncertainty(fit,
    uncertainty = "hessian", finishsamples = 100, verbose = 0)))
  finite <- suppressWarnings(suppressMessages(ctOptimUncertainty(fit,
    uncertainty = "hessian", finishsamples = 100, verbose = 0,
    control = list(analyticHessian = FALSE))))

  # Agreement is the check that the exact one is right; the finite difference
  # is the independent referee, since it knows only the gradient's outputs.
  expect_equal(exact$uncertainty$hessian, finite$uncertainty$hessian,
    tolerance = 1e-4)
  expect_equal(sqrt(diag(exact$estimate$cov)), sqrt(diag(finite$estimate$cov)),
    tolerance = 1e-4)

  expect_identical(exact$uncertainty$hessian, t(exact$uncertainty$hessian))

  # And it is recorded as exact, so a fit says which one produced its intervals.
  expect_match(exact$uncertainty$details$hessian$source, "forward-mode")
  expect_null(finite$uncertainty$details$hessian)
})

test_that("unsupported uncertainty methods are refused by name, not silently", {
  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()[1:24, ]
  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

  # The score-based methods are supported now that the engines produce
  # per-subject gradients directly; `fullbootstrap` is not, because it
  # re-optimises each resample and so needs the model rebuilt rather than
  # re-evaluated. Refusing by name beats producing a plausible-looking
  # covariance from a method that did not actually run.
  expect_error(ctOptimUncertainty(julia_fit, uncertainty = "fullbootstrap"),
    "not available for backend")
  expect_true(all(c("opg", "sandwich", "bootstrap") %in%
      ctsem:::.ctBackendUncertaintySupported))
})
