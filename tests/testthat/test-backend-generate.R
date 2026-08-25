# Posterior-predictive data generation for backend='julia'.
#
# Generation is not a recorder: it changes what the filter consumes, so the
# state carried forward is conditioned on the drawn data rather than the real
# data. That is what makes the result a draw from the model, and it is also what
# makes it easy to get subtly wrong -- a simulator that predicted each row from
# the *real* history would produce data that looks fine and is not a draw from
# anything.
#
# So the load-bearing test is an identity: re-running the ordinary likelihood on
# the generated dataset must reproduce the likelihood reported while generating
# it. Nothing that drew from the wrong covariance, or conditioned on the wrong
# history, can satisfy that.

.generate_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.generate_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8
  diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(20), function(i) {
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
  data$Y1[c(4, 20)] <- NA
  data
}

test_that("generated data is a draw from the model the filter conditions on", {
  skip_on_cran()
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  spec <- fit$model_spec
  nrows <- length(spec$times)

  set.seed(1)
  base <- matrix(stats::rnorm(nrows), 1, nrows)
  drawn <- ctsem:::.ctBackendGenerate(fit, fit$estimate$raw, base)

  expect_equal(dim(drawn$Y), c(1L, nrows))
  # The identity that defines it: the filter's own likelihood for the generated
  # data is the likelihood it reported while generating it.
  refit <- data
  refit$Y1 <- as.numeric(drawn$Y)
  respec <- suppressMessages(ctFit(refit, model, backend = "julia", fit = FALSE))
  regenerated <- ctJuliaEvaluate(ctsem:::.ctBackendAsModel(respec),
    fit$estimate$raw, gradient = FALSE)$value
  expect_equal(sum(drawn$subject_loglik), regenerated, tolerance = 1e-10)
  expect_equal(sum(drawn$llrow), regenerated, tolerance = 1e-10)

  # Missingness is preserved: an entry that was not observed is not invented.
  expect_equal(which(is.na(as.numeric(drawn$Y))), which(is.na(data$Y1)))
  expect_equal(drawn$llrow[which(is.na(data$Y1))], rep(0, sum(is.na(data$Y1))))

  # The draws drive it -- different normals, different data.
  set.seed(2)
  other <- ctsem:::.ctBackendGenerate(fit, fit$estimate$raw,
    matrix(stats::rnorm(nrows), 1, nrows))
  expect_false(isTRUE(all.equal(as.numeric(other$Y), as.numeric(drawn$Y))))
})

test_that("the observed data is an unremarkable draw from the fitted model", {
  skip_on_cran()
  skip_without_julia()
  # The calibration check the identity above cannot make: at the maximum
  # likelihood estimate, the observed data's log likelihood should sit somewhere
  # ordinary in the distribution of generated ones. A generator that drew from
  # too small a covariance, or that conditioned each row on the real history
  # instead of the drawn one, would put the observed value far into a tail.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  nrows <- length(fit$model_spec$times)

  set.seed(2)
  generated <- replicate(200, sum(ctsem:::.ctBackendGenerate(fit, fit$estimate$raw,
    matrix(stats::rnorm(nrows), 1, nrows))$subject_loglik))
  quantile <- mean(generated < fit$estimate$loglik)
  expect_gt(quantile, .05)
  expect_lt(quantile, .95)
  # And the spread is real, not a degenerate point mass.
  expect_gt(stats::sd(generated), 1)
})

test_that("ctGenerateFromFit returns what the posterior predictive tools expect", {
  skip_on_cran()
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  nrows <- length(fit$model_spec$times)

  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 20, cores = 1)
  expect_equal(dim(generated$generated$Y), c(20L, nrows, 1L))
  expect_identical(dimnames(generated$generated$Y)[[3]], model$manifestNames)
  expect_equal(dim(generated$generated$llrow), c(20L, nrows))
  # A row with nothing observed contributes nothing, and is reported as NA
  # rather than as a zero that would read as a likelihood.
  expect_true(all(is.na(generated$generated$llrow[, which(is.na(data$Y1))])))
  expect_true(all(is.na(generated$generated$Y[, which(is.na(data$Y1)), 1])))

  # fullposterior needs draws, and says so rather than silently using the mode.
  # A default fit has them, since ctFit() finishes with ctOptimUncertainty(), so
  # the refusal is only reachable for a fit that deliberately skipped it.
  estonly <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctGenerateFromFit(estonly, nsamples = 5, fullposterior = TRUE, cores = 1),
    "ctOptimUncertainty")

  fromposterior <- ctGenerateFromFit(fit, nsamples = 10, fullposterior = TRUE,
    cores = 1)
  expect_equal(dim(fromposterior$generated$Y), c(10L, nrows, 1L))
})

test_that("the posterior predictive tools run on a backend fit", {
  skip_on_cran()
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 20, cores = 1)

  predictive <- suppressMessages(ctsem:::ctPostPredData(generated))
  expect_true(all(c("row", "variable", "sample", "value", "id", "Time",
    "TimeInterval", "obsValue") %in% names(predictive)))
  # Every row of every sample, for each manifest plus the row likelihood.
  expect_equal(nrow(predictive), 20 * length(fit$model_spec$times) * 2)
  # The observed column is the data, not a copy of the generated one.
  observed <- predictive[predictive$variable == "Y1" & predictive$sample == 1, ]
  expect_equal(observed$obsValue[order(observed$row)], data$Y1)

  plots <- suppressWarnings(suppressMessages(ctPostPredPlots(generated)))
  expect_true(length(plots) > 0)
  expect_true(all(vapply(plots, function(p) inherits(p, "ggplot"), logical(1))))

  covcheck <- suppressWarnings(suppressMessages(
    ctFitCovCheck(generated, plot = FALSE, lags = 0:2, cores = 1, nsamples = 10)))
  expect_true(nrow(covcheck) > 0)
  expect_true("Sig" %in% names(covcheck))
})

test_that("ctPostPredData(residuals=TRUE) works, for stan too", {
  skip_on_cran()
  skip_without_julia()
  # This branch could never run: the residual rows lacked the id/time columns
  # the rbind below them needs, on every backend including stan. Fixed here, so
  # guarded here.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 3, cores = 1)

  predictive <- suppressMessages(ctsem:::ctPostPredData(generated, residuals = TRUE))
  expect_true("Y1 std. res." %in% predictive$variable)
  expect_false(anyNA(predictive$Time[predictive$variable == "Y1 std. res."]))

  skip_if_not_installed("rstan")
  stangenerated <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 2,
    cores = 1))
  stanpredictive <- suppressMessages(ctsem:::ctPostPredData(stangenerated,
    residuals = TRUE))
  expect_true(any(grepl("std. res.", stanpredictive$variable)))
})
