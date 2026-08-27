# `ctLaplaceCheck()`: measuring how much of a Laplace fit is the approximation.
#
# The engine's own suite proves the quadrature -- that it is a Gauss-Hermite
# rule, that one node reproduces the Laplace value exactly, and that the
# correction is the Newton step it claims to be. What is left for here is the
# R side: that the wiring reaches the right engine entry points, that the
# parameter labels line up with the raw vector the correction is expressed in,
# and that every unsupported way of asking refuses rather than returning
# something that looks like an answer.

# A model with a random effect on an *identity-transformed* parameter. Laplace
# is exact for that integrand, so the quadrature must find essentially nothing
# to correct -- which is the sharpest available end-to-end check that the two
# are integrating the same thing through the whole R-to-engine path.
.check_test_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.check_test_data <- function(nsubjects = 25, nobs = 6) {
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

.check_cache <- new.env(parent = emptyenv())
.check_fit <- function() {
  if (!exists("fit", envir = .check_cache, inherits = FALSE)) {
    assign("fit", suppressMessages(ctFit(.check_test_data(), .check_test_model(),
      backend = "julia", intoverpop = "laplace",
      optimcontrol = list(finishsamples = 100))), envir = .check_cache)
  }
  get("fit", envir = .check_cache, inherits = FALSE)
}

test_that("the check finds nothing to correct where Laplace is exact", {
  skip_without_julia()
  fit <- .check_fit()
  check <- ctLaplaceCheck(fit, nodes = 5)

  expect_s3_class(check, "ctLaplaceCheck")
  expect_equal(check$nsubjects, 25L)
  # An identity-transformed random effect entering the state mean linearly makes
  # the integrand exactly Gaussian in z. The rule integrates the same Gaussian
  # whatever the node count, so the gap is quadrature error and nothing else.
  expect_lt(abs(check$gap), 1e-6)
  expect_equal(check$gap, check$quadrature - check$laplace)
  expect_equal(check$gap_per_subject, check$gap / 25)
  expect_lt(max(abs(check$parameters$delta_se)), 0.02)
})

test_that("the correction table describes the raw vector it corrects", {
  skip_without_julia()
  fit <- .check_fit()
  check <- ctLaplaceCheck(fit, nodes = 3)
  pars <- check$parameters

  expect_equal(nrow(pars), length(fit$estimate$raw))
  expect_equal(pars$estimate, as.numeric(fit$estimate$raw))
  expect_equal(pars$corrected, pars$estimate + pars$delta)
  # The population scale is the last raw position and must be named for the
  # parameter it is the scale of, not by its index.
  expect_true("popsd_mmean" %in% pars$parameter)
  expect_equal(which(pars$parameter == "popsd_mmean"),
    fit$model_spec$laplace$sd_index)
  expect_true("mmean" %in% pars$parameter)
  expect_false(any(duplicated(pars$parameter)))
  expect_true(all(is.finite(pars$se)))
})

test_that("more nodes cost more and change nothing here", {
  skip_without_julia()
  fit <- .check_fit()
  one <- ctLaplaceCheck(fit, nodes = 1, correction = FALSE)
  nine <- ctLaplaceCheck(fit, nodes = 9, correction = FALSE)
  # One node *is* Laplace, by construction rather than by approximation.
  expect_equal(one$quadrature, one$laplace, tolerance = 1e-9)
  expect_equal(nine$quadrature, one$quadrature, tolerance = 1e-6)
})

test_that("the check refuses what it cannot do", {
  skip_without_julia()
  fit <- .check_fit()
  expect_error(ctLaplaceCheck(list()), "backend='julia'")

  augmented <- suppressMessages(ctFit(.check_test_data(), .check_test_model(),
    backend = "julia", intoverpop = TRUE, optimcontrol = list(estonly = TRUE)))
  expect_error(ctLaplaceCheck(augmented), "intoverpop='laplace'")

  # Without a Hessian the correction cannot be formed, and saying so beats
  # returning a table of zeros.
  bare <- fit
  # Backend fits keep it at `fit$uncertainty`; the Stan path at `fit$stanfit`.
  bare$uncertainty$hessian <- NULL
  bare$stanfit$uncertainty$hessian <- NULL
  expect_warning(result <- ctLaplaceCheck(bare, nodes = 3), "Hessian")
  expect_null(result$parameters)
})

# Subjects nested in studies, both levels carrying an identity-transformed
# MANIFESTMEANS. Laplace is exact for that integrand at *both* depths, so the
# recursive rule has to agree with it -- which is the sharpest end-to-end test
# available for the nested case, because every constant in the recursion has to
# cancel against the joint block factorization rather than against itself.
.check_nested_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"),
    id = c("subject", "study"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model$pars$indvarying_study <- FALSE
  model$pars$indvarying_study[model$pars$param %in% "mmean"] <- TRUE
  model
}

.check_nested_data <- function(nstudy = 6, npersub = 4, nobs = 5) {
  set.seed(20260827)
  rows <- list(); sid <- 0
  for (g in seq_len(nstudy)) {
    studyeffect <- stats::rnorm(1, 0, 0.4)
    for (j in seq_len(npersub)) {
      sid <- sid + 1
      intercept <- 1.5 + studyeffect + stats::rnorm(1, 0, 0.5)
      state <- stats::rnorm(1, 0, 0.5); out <- numeric(nobs)
      for (t in seq_len(nobs)) {
        if (t > 1) {
          decay <- exp(-0.4)
          state <- decay * state + stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - decay^2)))
        }
        out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
      }
      rows[[length(rows) + 1L]] <- data.frame(subject = sid, study = g,
        time = seq_len(nobs) - 1, Y1 = out)
    }
  }
  do.call(rbind, rows)
}

test_that("the check works when subjects are nested in studies", {
  skip_without_julia()
  fit <- suppressMessages(ctFit(.check_nested_data(), .check_nested_model(),
    backend = "julia", intoverpop = "laplace",
    optimcontrol = list(finishsamples = 50)))
  expect_equal(as.integer(fit$model_spec$laplace$nlevels), 2L)

  check <- ctLaplaceCheck(fit, nodes = 3)
  expect_lt(abs(check$gap), 1e-5)
  # Per *unit*, not per subject: a study's integral does not decompose over its
  # members, so there are six terms for twenty-four subjects.
  expect_equal(check$nsubjects, 24L)
  expect_equal(nrow(check$parameters), length(fit$estimate$raw))
  # Two levels means two population scales, each named for its own level.
  scales <- grep("^popsd_", check$parameters$parameter, value = TRUE)
  expect_length(scales, 2L)
  expect_true(any(grepl("study", scales)))
  # Nothing to correct where Laplace is exact.
  expect_lt(max(abs(check$parameters$delta_se)), 0.05)
})

test_that("the fit records how hard the optimizer worked", {
  skip_without_julia()
  fit <- .check_fit()
  expect_false(isTRUE(fit$estimate$stalled))
  expect_true(is.finite(fit$estimate$f_calls))
  expect_gte(fit$estimate$f_calls, fit$estimate$iterations)
  expect_true(is.finite(fit$estimate$g_calls))
})
