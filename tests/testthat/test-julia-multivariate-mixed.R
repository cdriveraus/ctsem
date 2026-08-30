# Several latents, several kinds of indicator, several random effects, at once.
#
# Each of those was covered separately and none of them together. The ordinal
# suite is entirely one-latent; the two-state adjoint regressions get their
# second state by augmenting a one-latent model whose indicators are all of one
# kind; and the mixed-measurement studies varied the measurement model while
# holding the latent process at one dimension.
#
# That matters because of what the gap has already hidden once. With a single
# latent state the covariance cotangent is a scalar and trivially symmetric, and
# the binary reverse pass assumed that symmetry in general -- so a wrong adjoint
# passed every gradient test the suite had until a two-state case was written.
# The combination here is where a repeat of that would live: a cross-lagged
# process, a measurement model mixing ordinal, binary and Gaussian indicators
# across different latents, and random effects on both.

.mvmix_data <- function(nsubjects = 30, nobs = 8, seed = 404) {
  set.seed(seed)
  tau <- c(-0.8, 0.3, 1.4)
  invlog <- function(x) 1 / (1 + exp(-x))
  drawcat <- function(eta) {
    cum <- sapply(seq_along(tau), function(k) invlog(tau[k] - eta))
    p <- cbind(cum, 1)
    p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
    apply(p, 1, function(pr) sample.int(length(pr), 1, prob = pmax(pr, 0)))
  }
  cint1 <- stats::rnorm(nsubjects, 0, 0.4)
  cint2 <- stats::rnorm(nsubjects, 0, 0.4)
  d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
      manifestNames = c("e1", "e2"), latentNames = c("eta1", "eta2"),
      LAMBDA = diag(2),
      DRIFT = matrix(c(-0.4, 0.15, -0.10, -0.6), 2, 2, byrow = TRUE),
      DIFFUSION = matrix(c(0.7, 0, 0.2, 0.6), 2, 2),
      MANIFESTVAR = diag(1e-6, 2), T0VAR = diag(1, 2),
      T0MEANS = matrix(0, 2, 1),
      CINT = matrix(c(cint1[i], cint2[i]), 2, 1),
      MANIFESTMEANS = matrix(0, 2, 1), Tpoints = nobs))
    one <- data.frame(ctGenerate(gen, n.subjects = 1, Tpoints = nobs,
      backend = "r"))
    one$id <- i
    one
  }))
  d$o1 <- drawcat(d$e1)
  d$o2 <- drawcat(d$e1)
  d$b1 <- stats::rbinom(nrow(d), 1, invlog(d$e1))
  d$y1 <- d$e2 + stats::rnorm(nrow(d), 0, 0.4)
  d$y2 <- d$e2 + stats::rnorm(nrow(d), 0, 0.4)
  d$e1 <- NULL
  d$e2 <- NULL
  d
}

# Ordinal and binary indicators load on the first latent, Gaussian ones on the
# second, so the measurement types are spread across the process rather than
# stacked on one state that a bug could average over.
.mvmix_model <- function() {
  lambda <- matrix(0, 5, 2)
  lambda[1:3, 1] <- 1
  lambda[4:5, 2] <- 1
  mvar <- diag(0, 5)
  mvar[4, 4] <- "mv1"
  mvar[5, 5] <- "mv2"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 5, manifestNames = c("o1", "o2", "b1", "y1", "y2"),
    latentNames = c("eta1", "eta2"), manifesttype = c(2L, 2L, 1L, 0L, 0L),
    ncategories = c(4L, 4L, 0L, 0L, 0L), LAMBDA = lambda,
    MANIFESTMEANS = matrix(0, 5, 1),
    CINT = matrix(c("cint1", "cint2"), 2, 1), T0MEANS = matrix(0, 2, 1),
    MANIFESTVAR = mvar)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% c("cint1", "cint2")] <- TRUE
  m
}

test_that("the adjoint matches forward mode with mixed indicators on two latents", {
  skip_without_julia()
  handle <- suppressWarnings(suppressMessages(ctFit(.mvmix_data(),
    .mvmix_model(), backend = "julia", intoverpop = "augmented", fit = FALSE)))
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  expect_gt(npar, 15)

  # Forward mode rather than a finite difference: it is exact, so a
  # disagreement is unambiguously the reverse pass rather than a step size.
  set.seed(7)
  for (trial in 1:2) {
    at <- stats::rnorm(npar, 0, 0.3)
    adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

# Thirty subjects and eight occasions, not fewer. This model has 23 free
# parameters, and below that the state-augmented route drives one of them to
# the flat region of its transform and is correctly reported as not converged:
# measured at 20 subjects, augmented saturates at both six and eight occasions
# while Laplace converges, with the two log likelihoods agreeing to 0.02 either
# way. That is a boundary case worth knowing about but not what this test is
# for, which is whether the two routes integrate the same thing.
test_that("laplace and augmented agree on a multivariate mixed model", {
  skip_without_julia()
  d <- .mvmix_data()
  m <- .mvmix_model()
  augmented <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "augmented", optimcontrol = list(estonly = TRUE))))
  laplace <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "laplace", optimcontrol = list(estonly = TRUE))))

  expect_true(isTRUE(augmented$estimate$converged))
  expect_true(isTRUE(laplace$estimate$converged))
  # The random effects are identity-transformed intercepts, so the Laplace
  # approximation is close to exact here and the two routes are integrating
  # essentially the same thing. A wide gap would mean one of them is not.
  expect_equal(as.numeric(laplace$estimate$loglik),
    as.numeric(augmented$estimate$loglik), tolerance = 0.5)
})

test_that("the Laplace correction handles more than one random effect", {
  skip_without_julia()
  d <- .mvmix_data()
  m <- .mvmix_model()
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "laplace", optimcontrol = list(finishsamples = 100))))
  skip_if_not(isTRUE(fit$estimate$converged), "fit did not converge")

  # Two random effects per unit means the quadrature is a `nodes^2` product
  # rule over the block rather than a line of nodes, which is the part of the
  # recursion a one-random-effect model never exercises.
  corrected <- suppressWarnings(ctLaplaceCorrect(fit, draws = "normal",
    nodes = 5, finishsamples = 100))
  expect_s3_class(corrected$laplace_correction, "ctLaplaceCorrection")
  expect_equal(corrected$laplace_correction$check$dropped_directions, 0L)
  expect_true(all(is.finite(as.numeric(corrected$estimate$raw))))
  # Identity-transformed intercepts again: there should be very little here to
  # correct, and a large correction would mean the recursion is not integrating
  # what the fit maximised.
  expect_lt(max(abs(corrected$laplace_correction$check$parameters$delta_se)), 0.5)
})
