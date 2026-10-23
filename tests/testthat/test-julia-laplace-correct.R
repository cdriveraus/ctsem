# `ctFitUncertainty(fit, 'is' | 'checked', control = list(likelihood =
# 'quadrature'))`: the draws behind the intervals taken on the quadrature
# hybrid the continuation climbs (`.ctLaplaceContinue()`), by the method asked
# for.
#
# The engine's suite proves the rules and the climb. What is left for here is
# that the request changes the fit in the ways claimed and in no others: that
# a fit where Laplace is exact keeps its estimate and draws on it, that the
# record says what was done, and that each way of asking for something
# unsupported refuses rather than returning something that looks like an
# answer.

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

# A short run: one round of 40 draws is enough to show the machinery working
# end to end, and not to pin a posterior.
.quadrature_is <- function(fit, method = "is", ...) {
  suppressWarnings(ctFitUncertainty(fit, uncertainty = method, finishsamples = 60,
    control = list(likelihood = "quadrature", nodes = 5, isitersize = 40,
      imisMaxIter = 1L, isESS = 20, ...)))
}

test_that("a fit where Laplace is exact keeps its estimate, and records what was done", {
  skip_without_julia()
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .correct_fit()
  before <- as.numeric(fit$estimate$raw)

  set.seed(9)
  corrected <- .quadrature_is(fit)
  # The screen finds Laplace exact here (the file comment), so the hybrid is
  # the fit's own objective, nothing is continued, and the draws are taken on
  # it.
  record <- corrected$laplace$correction
  expect_identical(record$status, "exact")
  expect_equal(as.numeric(corrected$estimate$raw), before, tolerance = 1e-6)
  expect_equal(nrow(corrected$estimate$rawposterior), 60L)
  expect_equal(ncol(corrected$estimate$rawposterior), length(fit$estimate$raw))
  expect_false(is.null(corrected$uncertainty$imis))
  expect_equal(unname(corrected$estimate$se), unname(sqrt(diag(corrected$estimate$cov))))

  # The record every method writes, so summaries describe these draws.
  settings <- corrected$uncertainty$settings
  expect_identical(settings$method, "is")
  expect_identical(settings$draws, "imis")
  expect_identical(settings$control$likelihood, "quadrature")
  expect_identical(corrected$uncertainty$draws, "imis")
  expect_match(suppressMessages(summary(corrected))$uncertaintyNote,
    "likelihood='quadrature'", fixed = TRUE)

  # And by 'checked', on the same density.
  set.seed(9)
  checked <- .quadrature_is(fit, method = "checked")
  expect_identical(checked$uncertainty$settings$method, "checked")
  expect_identical(checked$uncertainty$settings$control$likelihood, "quadrature")
  expect_equal(as.numeric(checked$estimate$raw), before, tolerance = 1e-6)
  expect_equal(nrow(checked$estimate$rawposterior), 60L)
})

test_that("correctEstimate = FALSE leaves the point where it was", {
  skip_without_julia()
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .correct_fit()
  set.seed(10)
  same <- .quadrature_is(fit, correctEstimate = FALSE)
  expect_equal(as.numeric(same$estimate$raw), as.numeric(fit$estimate$raw))
  expect_false(isTRUE(same$laplace$correction$applied))
  expect_identical(same$uncertainty$settings$control$likelihood, "quadrature")
})

test_that("requests it cannot serve are refused by name", {
  skip_without_julia()
  fit <- .correct_fit()
  # Only the methods that draw take a likelihood.
  expect_error(ctFitUncertainty(fit, "hessian",
    control = list(likelihood = "quadrature")), "applies to uncertainty = 'is'")
  # Each likelihood takes its own settings and the shared ones, nothing else.
  expect_error(ctFitUncertainty(fit, "is",
    control = list(likelihood = "quadrature", particles = 100)),
    "control$particles is not used with likelihood = 'quadrature'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, "is", control = list(likelihood = "exact")),
    "should be one of")
  # The quadrature and particle filter are the julia engine's.
  expect_error(ctFitUncertainty(ctstantestfit, "is",
    control = list(likelihood = "quadrature")), "only available for backend='julia'")
  # The particle filter is sampled by 'is' alone.
  expect_error(ctFitUncertainty(fit, "checked",
    control = list(likelihood = "particle")), "sampled by uncertainty = 'is'")

  augmented <- suppressMessages(ctFit(.correct_test_data(),
    .correct_test_model(), backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctFitUncertainty(augmented, "is",
    control = list(likelihood = "quadrature")), "intoverpop = 'laplace'",
    fixed = TRUE)

  # Without uncertainty there is no curvature for the continuation's metric
  # or its certification, which is a different failure from the fit being
  # the wrong kind.
  estonly <- suppressMessages(ctFit(.correct_test_data(),
    .correct_test_model(), backend = "julia", intoverpop = "laplace",
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctFitUncertainty(estonly, "is",
    control = list(likelihood = "quadrature")), "covariance")
})
