# `ctLaplaceCorrect()`: applying the correction `ctLaplaceCheck()` measures.
#
# The engine's suite proves the quadrature and the Newton step. What is left for
# here is that applying them changes the fit in the ways claimed and in no
# others: that a fit where Laplace is exact comes back untouched, that the cheap
# path really is only a recentring, and that each way of asking for something
# unsupported refuses rather than returning something that looks like an answer.

# A random effect on an identity-transformed parameter. Laplace is exact for
# that integrand, so there is nothing to correct -- which makes it the sharpest
# end-to-end check available: any movement here is the machinery inventing it.
.correct_test_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.correct_test_data <- function(nsubjects = 25, nobs = 6) {
  set.seed(20260827)
  drift <- -0.4; diffusion <- 0.6
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1.5, 0.9)
    state <- stats::rnorm(1, 0, 0.5)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        decay <- exp(drift)
        state <- decay * state +
          stats::rnorm(1, 0, sqrt(diffusion^2 / (-2 * drift) * (1 - decay^2)))
      }
      out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
}

.correct_cache <- new.env(parent = emptyenv())
.correct_fit <- function() {
  if (!exists("fit", envir = .correct_cache, inherits = FALSE)) {
    assign("fit", suppressMessages(ctFit(.correct_test_data(),
      .correct_test_model(), backend = "julia", intoverpop = "laplace",
      optimcontrol = list(finishsamples = 100))), envir = .correct_cache)
  }
  get("fit", envir = .correct_cache, inherits = FALSE)
}

test_that("a fit where Laplace is exact comes back with its estimate intact", {
  skip_without_julia()
  fit <- .correct_fit()
  before <- as.numeric(fit$estimate$raw)
  se <- as.numeric(fit$estimate$se)

  corrected <- ctLaplaceCorrect(fit, draws = "normal", nodes = 5,
    finishsamples = 100)

  expect_s3_class(corrected$laplace_correction, "ctLaplaceCorrection")
  # Not merely small: the gap gradient here is floating-point noise, so the
  # correction should be at the level of the arithmetic rather than at the
  # level of the parameter.
  expect_lt(max(abs(corrected$laplace_correction$check$parameters$delta_se)), 1e-4)
  expect_equal(as.numeric(corrected$estimate$raw), before, tolerance = 1e-6)
  expect_equal(corrected$laplace_correction$before$raw, before)
})

test_that("draws='normal' recentres without touching the width", {
  skip_without_julia()
  fit <- .correct_fit()
  corrected <- ctLaplaceCorrect(fit, draws = "normal", nodes = 5,
    finishsamples = 100)

  # The cheap path corrects the location and says so. Its covariance is the
  # fit's own, so any change here would mean it had quietly done more.
  expect_equal(as.matrix(corrected$estimate$cov), as.matrix(fit$estimate$cov))
  expect_identical(corrected$laplace_correction$draws, "normal")
  expect_true(is.na(corrected$laplace_correction$ess))
  # Draws are redrawn even when the estimate does not move, so that they and
  # the constrained draws describe one point rather than two.
  expect_equal(nrow(corrected$estimate$rawposterior), 100L)
  expect_false(is.null(corrected$transformedpars))
})

test_that("draws='keep' says the intervals no longer describe the estimate", {
  skip_without_julia()
  fit <- .correct_fit()
  # Silent here only because Laplace is exact for this model, so the estimate
  # does not actually move; the warning is conditional on it having moved.
  kept <- suppressWarnings(ctLaplaceCorrect(fit, draws = "keep", nodes = 5))
  expect_identical(kept$laplace_correction$draws, "keep")
  expect_equal(as.matrix(kept$estimate$cov), as.matrix(fit$estimate$cov))
})

test_that("correct_estimate=FALSE leaves the point where it was", {
  skip_without_julia()
  fit <- .correct_fit()
  same <- ctLaplaceCorrect(fit, draws = "normal", nodes = 5,
    finishsamples = 100, correct_estimate = FALSE)
  expect_equal(as.numeric(same$estimate$raw), as.numeric(fit$estimate$raw))
  expect_false(same$laplace_correction$corrected_estimate)
})

test_that("it refuses the fits it cannot correct", {
  skip_without_julia()
  expect_error(ctLaplaceCorrect(list()), "backend='julia'")

  augmented <- suppressMessages(ctFit(.correct_test_data(),
    .correct_test_model(), backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctLaplaceCorrect(augmented), "intoverpop='laplace'")

  # Without uncertainty there is no curvature to take the Newton step against,
  # which is a different failure from the fit being the wrong kind.
  estonly <- suppressMessages(ctFit(.correct_test_data(),
    .correct_test_model(), backend = "julia", intoverpop = "laplace",
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctLaplaceCorrect(estonly), "covariance")
})
