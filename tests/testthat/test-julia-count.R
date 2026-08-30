# Count manifest variables: `manifesttype = 3`, Poisson with a log link.
#
# The engine's own suite proves the kernel -- that the log likelihood, score and
# information match the closed forms exactly, and that the 21-node rule matches
# dense numerical integration. What is left for here is the R side: that the
# model accepts and describes the type, that the data checks fire, that the
# backend that cannot fit one refuses instead of returning a number, and that
# the reverse-mode gradient is right on a real model.

.count_data <- function(nsubjects = 40, nobs = 8, nind = 2, seed = 5) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
      manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
      DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
      MANIFESTVAR = matrix(1e-6), T0VAR = matrix(1), T0MEANS = matrix(0),
      CINT = matrix(0), MANIFESTMEANS = matrix(0), Tpoints = nobs))
    one <- data.frame(ctGenerate(gen, n.subjects = 1, Tpoints = nobs,
      backend = "r"))
    one$id <- i
    one
  }))
  for (j in seq_len(nind)) {
    d[[paste0("c", j)]] <- stats::rpois(nrow(d), exp(1.2 + d$eta))
  }
  d$eta <- NULL
  d
}

.count_model <- function(nind = 2) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = nind, manifestNames = paste0("c", seq_len(nind)),
    latentNames = "eta1", manifesttype = rep(3L, nind),
    LAMBDA = matrix(1, nind, 1),
    MANIFESTMEANS = matrix("mm", nind, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, nind))))
  m$pars$indvarying <- FALSE
  m
}

# The same model with something for the Laplace route to integrate over.
# `intoverpop='laplace'` refuses a model with no varying parameter, correctly,
# so the two-route comparisons need one.
.count_model_varying <- function(nind = 2) {
  m <- .count_model(nind)
  m$pars$indvarying[m$pars$param %in% "mm"] <- TRUE
  m
}

test_that("ctModel takes manifesttype 3 and needs no categories for it", {
  m <- .count_model()
  expect_equal(unname(m$manifesttype), c(3L, 3L))
  # `ncategories` describes ordinal variables only, and a count must not be
  # made to invent one.
  expect_true(all(m$ncategories == 0L))
  expect_null(m$THRESHOLDS)
  expect_error(suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 1, n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 5L, LAMBDA = matrix(1)))), "manifesttype must be")
})

test_that("print names counts and the link they use", {
  described <- paste(utils::capture.output(print(.count_model())),
    collapse = " ")
  expect_match(described, "c1, c2 \\(count, Poisson log link\\)")
})

test_that("count data is checked against the model", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  m <- .count_model()

  negative <- d
  negative$c1[3] <- -1
  expect_error(suppressMessages(ctFit(negative, m, backend = "julia",
    fit = FALSE)), "negative values")

  fractional <- d
  fractional$c1[3] <- 1.5
  expect_error(suppressMessages(ctFit(fractional, m, backend = "julia",
    fit = FALSE)), "whole number")

  # Legal, and worth saying: a constant count identifies no rate.
  constant <- d
  constant$c1 <- 0L
  expect_warning(suppressMessages(ctFit(constant, m, backend = "julia",
    fit = FALSE)), "not identified")
})

test_that("stan refuses counts rather than treating them as Gaussian", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  expect_error(suppressMessages(ctFit(d, .count_model(), backend = "stan",
    optimcontrol = list(estonly = TRUE))), "need backend=\"julia\"")
})

test_that("the adjoint matches forward mode on a count model", {
  skip_without_julia()
  d <- .count_data()
  fit <- suppressWarnings(suppressMessages(ctFit(d, .count_model(),
    backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE))))
  handle <- structure(fit$model_spec,
    class = c("ctJuliaModel", "ctFitModel"))
  est <- as.numeric(fit$estimate$raw)

  # Anchored at the estimate rather than drawn from zero, which is what the
  # binary and ordinal suites do. Their likelihoods are bounded, so any raw
  # draw gives a representable objective; a count's linear predictor is
  # *exponentiated*, and the same draw from zero puts the Poisson rate near
  # 1e29. Comparing two gradients there compares two numbers that have lost
  # every significant digit -- measured, it reported a 50% disagreement where
  # the same code agrees to 1e-16 anywhere the objective is representable.
  set.seed(11)
  for (trial in 1:3) {
    at <- est + stats::rnorm(length(est), 0, 0.3)
    adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

test_that("a count model recovers what generated it, both routes", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()
  for (method in c("augmented", "laplace")) {
    fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = method, optimcontrol = list(estonly = TRUE))))
    expect_true(isTRUE(fit$estimate$converged))
    means <- summary(fit)$popmeans
    # Wide, because this is one dataset of forty subjects and the test is that
    # the estimator is pointed at the right place, not that it is precise.
    expect_equal(unname(means["drift_eta1", "mean"]), -0.4, tolerance = 0.35)
    expect_equal(unname(means["mm", "mean"]), 1.2, tolerance = 0.35)
  }
})

test_that("the two routes agree on a count model", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()
  augmented <- suppressWarnings(suppressMessages(ctFit(d, m,
    backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE))))
  laplace <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "laplace", optimcontrol = list(estonly = TRUE))))
  expect_equal(as.numeric(laplace$estimate$loglik),
    as.numeric(augmented$estimate$loglik), tolerance = 1)
})

test_that("ctGenerate draws counts rather than continuous values", {
  skip_without_julia()
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("c1", "c2"), latentNames = "eta1",
    manifesttype = c(3L, 3L), LAMBDA = matrix(1, 2, 1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
    MANIFESTVAR = diag(0, 2), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(1.0, 2, 1), Tpoints = 6)))
  set.seed(4)
  d <- data.frame(ctGenerate(gen, n.subjects = 15, Tpoints = 6,
    backend = "julia"))
  values <- c(d$c1, d$c2)
  values <- values[!is.na(values)]
  expect_true(length(values) > 0)
  # Whole numbers, none negative: the R generator has no measurement link and
  # would return continuous values here, which is the failure this guards.
  expect_true(all(values >= 0))
  expect_true(all(abs(values - round(values)) < 1e-8))
  expect_gt(length(unique(values)), 1)
})
