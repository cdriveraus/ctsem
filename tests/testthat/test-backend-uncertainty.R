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

  # The draws, the covariance and the standard errors are labelled by raw
  # parameter, on both backends. All of them or none of them: the comparison
  # just above is julia's `estimate$cov` against stan's `stanfit$cov`, so a name
  # on one side only makes the same quantity a different object depending on the
  # backend; and `sd(rawposterior)` is asserted equal to `sqrt(diag(cov))`
  # further down this file, which naming one of that pair alone would break.
  # Spelled out rather than compared to each other, because two NULLs are equal.
  parnames <- c("drift", "diff", "mvar", "mmean", "t0v")
  expect_identical(colnames(julia_unc$estimate$cov), parnames)
  expect_identical(names(julia_unc$estimate$se), parnames)
  expect_identical(colnames(julia_unc$estimate$rawposterior), parnames)
  expect_identical(colnames(stan_unc$stanfit$cov), parnames)
  expect_identical(colnames(stan_unc$stanfit$rawposterior), parnames)

  # And already named on the fits ctFit() returned, not only after a manual
  # ctOptimUncertainty(): stanoptimis() computes uncertainty on a stub that
  # carries no model to read names from, so ctFit() names them once the object
  # is assembled.
  expect_identical(colnames(julia_fit$estimate$cov), parnames)
  expect_identical(colnames(stan_fit$stanfit$cov), parnames)

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

test_that("one entry point, one contract: the arguments mean the same on both backends", {
  skip_without_julia()

  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

  # draws='empirical' -- which uncertainty='bootstrap' resolves to -- must
  # actually return the bootstrap draws, not a normal cloud with the
  # bootstrap's covariance. The signature is exact: empirical draws are the
  # sample the covariance was computed from, so their SD *is* sqrt(diag(cov)).
  # Normal draws are a fresh sample and miss it by sampling error.
  boot <- suppressWarnings(suppressMessages(ctOptimUncertainty(fit,
    uncertainty = "bootstrap", finishsamples = 200, cores = 1, verbose = 0)))
  expect_identical(boot$uncertainty$settings$draws, "empirical")
  expect_equal(apply(boot$estimate$rawposterior, 2, stats::sd),
    sqrt(diag(boot$estimate$cov)), tolerance = 1e-10)

  # finishsamples=NULL reuses the fit's existing draw count on both backends.
  # This hardcoded 1000 on julia, so a second call silently resampled up.
  small <- suppressWarnings(suppressMessages(ctOptimUncertainty(fit,
    uncertainty = "hessian", finishsamples = 40, cores = 1, verbose = 0)))
  again <- suppressWarnings(suppressMessages(ctOptimUncertainty(small,
    uncertainty = "hessian", finishsamples = NULL, cores = 1, verbose = 0)))
  expect_equal(nrow(again$estimate$rawposterior), 40L)

  # control$parsteps is a stanoptimis() concept. It reached
  # ctOptimComputeUncertainty(), which never reads it, and was recorded in
  # $settings$control as though honoured -- same standard errors, a request
  # apparently granted. Refused by name instead.
  expect_error(ctOptimUncertainty(fit, uncertainty = "hessian",
    finishsamples = 20, control = list(parsteps = 1L)),
    "only available for backend='stan'", fixed = TRUE)

  # And the class guard names both backends' fits rather than one class.
  expect_error(ctOptimUncertainty(list(a = 1)), "ctJuliaFit")
})

test_that("unsupported uncertainty methods are refused by name, not silently", {
  skip_without_julia()
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

# uncertainty='stored' is backend-neutral, which is the point of it: the cheap
# redraw used to exist only as ctFitAddSamples(), which writes into
# fit$stanfit and so could never work here. Nothing about drawing from a
# covariance is backend-specific, and this asserts that the julia route reaches
# the same code and leaves the fit in the same shape.
test_that("uncertainty='stored' redraws a julia fit without touching the engine", {
  skip_without_julia()
  model <- .backend_uncertainty_model()
  data <- .backend_uncertainty_data()[1:96, ]
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

  cov0 <- fit$estimate$cov
  se0 <- fit$estimate$se
  hess0 <- fit$uncertainty$hessian

  set.seed(29)
  redrawn <- suppressMessages(ctOptimUncertainty(fit, uncertainty = "stored",
    finishsamples = 220, cores = 1))

  expect_equal(nrow(redrawn$estimate$rawposterior), 220L)
  expect_equal(unname(redrawn$estimate$cov), unname(cov0), tolerance = 0)
  expect_equal(unname(redrawn$estimate$se), unname(se0), tolerance = 0)
  expect_equal(redrawn$uncertainty$hessian, hess0, tolerance = 0)
  # `method` keeps naming what produced the covariance; `redrawn` says the
  # draws were regenerated from it.
  expect_identical(redrawn$uncertainty$settings$method, "hessian")
  expect_true(isTRUE(redrawn$uncertainty$settings$redrawn))
  expect_identical(colnames(redrawn$estimate$rawposterior),
    colnames(fit$estimate$rawposterior))

  # The same normals the shared helper would produce, so the two backends draw
  # identically given the same covariance and seed.
  set.seed(29)
  direct <- ctsem:::ctOptimNormalDraws(as.numeric(fit$estimate$raw), cov0, 220)
  expect_equal(unname(redrawn$estimate$rawposterior), unname(direct), tolerance = 0)

  # And no engine work: the exact Hessian is the expensive part of every other
  # method here, and this path must not ask for it.
  asked <- 0L
  trace(ctsem:::.ctBackendHessian, tracer = function() asked <<- asked + 1L,
    where = asNamespace("ctsem"), print = FALSE)
  on.exit(untrace(ctsem:::.ctBackendHessian, where = asNamespace("ctsem")),
    add = TRUE)
  suppressMessages(ctOptimUncertainty(fit, uncertainty = "stored",
    finishsamples = 10, cores = 1))
  expect_identical(asked, 0L)

  # Downstream reads it the way it reads any other draws.
  expect_silent(invisible(nrow(ctExtract(redrawn)$pop_DRIFT)))

  nocov <- fit
  nocov$estimate$cov <- NULL
  expect_error(ctOptimUncertainty(nocov, uncertainty = "stored"), "no usable one")
})

test_that("the value-only log probability is the value the gradient route returns", {
  skip_without_julia()

  # `imis_is` reads the log probability and nothing else, so it was paying for
  # a reverse pass per proposal draw and discarding it. Dropping that is only
  # free if the engine's value-only route returns the same number the adjoint's
  # own forward pass does -- they are separate code paths, and the adjoint sums
  # its chunk totals before adding the prior. Asserted bitwise, because
  # anything less would move a reported interval by an amount nobody could
  # then account for.
  fit <- suppressMessages(ctFit(.backend_uncertainty_data(),
    .backend_uncertainty_model(), backend = "julia", verbose = 0))
  withgrad <- ctsem:::.ctBackendLpgFunc(fit, gradient = TRUE)
  valueonly <- ctsem:::.ctBackendLpgFunc(fit, gradient = FALSE)
  est <- as.numeric(fit$estimate$raw)

  set.seed(3)
  npar <- length(est)
  points <- rbind(est, matrix(est, nrow = 15, ncol = npar, byrow = TRUE) +
      matrix(stats::rnorm(15 * npar, 0, .4), nrow = 15))
  a <- vapply(seq_len(nrow(points)),
    function(i) as.numeric(withgrad(points[i, ])), numeric(1))
  b <- vapply(seq_len(nrow(points)),
    function(i) as.numeric(valueonly(points[i, ])), numeric(1))
  expect_identical(a, b)

  # The gradient is there on one route and absent on the other, which is the
  # whole of the difference.
  expect_length(attr(withgrad(est), "gradient"), npar)
  expect_null(attr(valueonly(est), "gradient"))

  # And the invalid-point guard survives on both: `imis_is` needs a finite,
  # negligible weight at a draw the model cannot evaluate, not a failed batch.
  bad <- est; bad[1L] <- NaN
  expect_identical(as.numeric(withgrad(bad)), -1e100)
  expect_identical(as.numeric(valueonly(bad)), -1e100)
  expect_equal(attr(withgrad(bad), "gradient"), rep(0, npar))
})
