# Binary observations on the julia backend.
#
# The measurement update integrates the Bernoulli likelihood against the
# predicted state rather than moment-matching it to a Gaussian, so the tests
# that matter are the ones with an answer known independently of the code:
# a likelihood that can be written down, and a gradient that finite differences
# can confirm.

.jbin_data <- function(nsubjects = 40, nobs = 12, nindicators = 3, seed = 11) {
  invlog <- function(x) exp(x) / (1 + exp(x))
  set.seed(seed)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.001),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = nobs))
  latent <- ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r")
  d <- data.frame(latent)
  for (i in seq_len(nindicators)) {
    d[[paste0("b", i)]] <- stats::rbinom(nrow(d), 1, invlog(d$eta))
  }
  d$eta <- NULL
  d
}

.jbin_model <- function(nindicators = 3, free = TRUE) {
  args <- list(type = "ct", n.latent = 1, n.manifest = nindicators,
    manifestNames = paste0("b", 1:nindicators), latentNames = "eta1",
    LAMBDA = matrix(1, nindicators, 1),
    MANIFESTMEANS = matrix(0, nindicators, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, nindicators))
  if (!free) {
    args$DRIFT <- matrix(-0.3); args$DIFFUSION <- matrix(0.8)
    args$T0VAR <- matrix(1)
  }
  m <- suppressMessages(do.call(ctModel, args))
  m$manifesttype[] <- 1L
  m$pars$indvarying <- FALSE
  m
}

.jbin_spec <- function(d, m) {
  spec <- ctsem:::.ctJuliaPrepare(d, m, priors = FALSE,
    intoverpop = "augmented")
  structure(spec, class = c("ctJuliaModel", "ctFitModel"))
}

test_that("one binary observation has the likelihood it can be shown to have", {
  skip_on_cran()
  skip_without_julia()
  # T0MEANS 0, T0VAR 1, a single observation of 1: the likelihood is
  # int inv_logit(eta) N(eta;0,1) deta, which is 0.5 by symmetry.
  d <- data.frame(id = 1L, time = 0, b1 = 1L, x = NA_real_)
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
    manifestNames = c("b1", "x"), latentNames = "eta1",
    LAMBDA = matrix(c(1, 0), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(0), T0MEANS = matrix(0), DRIFT = matrix(-0.3),
    DIFFUSION = matrix("dif"), T0VAR = matrix(1),
    MANIFESTVAR = diag(c(0, 1), 2)))
  m$manifesttype <- c(1L, 0L)
  m$pars$indvarying <- FALSE
  handle <- .jbin_spec(d, m)
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  value <- ctJuliaEvaluate(handle, rep(0, npar), gradient = FALSE)$value
  expect_equal(as.numeric(value), log(0.5), tolerance = 1e-8)
})

test_that("the gradient matches finite differences, both ways of computing it", {
  skip_on_cran()
  skip_without_julia()
  # The reverse pass is hand written, and a reverse pass that does not know a
  # branch returns a confident wrong number rather than failing -- before it
  # learned this one it was out by 4e6 while forward mode was exact.
  handle <- .jbin_spec(.jbin_data(nsubjects = 20), .jbin_model())
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  at <- c(0.9, 0.2, -0.1)[seq_len(npar)]
  value <- function(x) as.numeric(ctJuliaEvaluate(handle, x,
    gradient = FALSE)$value)
  step <- 1e-5
  fd <- vapply(seq_len(npar), function(k) {
    up <- at; down <- at; up[k] <- up[k] + step; down[k] <- down[k] - step
    (value(up) - value(down)) / (2 * step)
  }, numeric(1))
  for (method in c("forward", "adjoint")) {
    got <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = method)$gradient)
    expect_equal(got, fd, tolerance = 1e-6,
      info = paste("gradient_method =", method))
  }
})

test_that("the two gradient methods agree with each other", {
  skip_on_cran()
  skip_without_julia()
  handle <- .jbin_spec(.jbin_data(nsubjects = 15, nindicators = 4),
    .jbin_model(nindicators = 4))
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(3)
  at <- stats::rnorm(npar, 0, 0.4)
  a <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  f <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(a, f, tolerance = 1e-8)
})

test_that("a binary model recovers what generated it", {
  skip_on_cran()
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(
    ctFit(.jbin_data(nsubjects = 50), .jbin_model(), backend = "julia",
      cores = 1, optimcontrol = list(estonly = TRUE))))
  est <- summary(fit)$popmeans
  expect_true(isTRUE(fit$estimate$converged))
  # Three binary indicators carry much less information than one continuous
  # one, so these are loose -- the point is that it lands near the truth rather
  # than at a transform boundary, which is what it used to do.
  expect_equal(unname(est["drift_eta1", "mean"]), -0.3, tolerance = 0.15)
  expect_equal(unname(est["diff_eta1", "mean"]), 0.8, tolerance = 0.3)
})

test_that("a saturated optimum is not reported as converged", {
  skip_on_cran()
  skip_without_julia()
  # Every ctsem transform is flat to machine precision by |raw| ~ 20: the
  # exponential underflows and the gradient is exactly zero, indistinguishable
  # from an optimum. An optimiser that lands there used to report success, and
  # the `stalled` check could not catch it because the optimiser had moved --
  # it moved too far.
  #
  # Checked on a Gaussian model because the guard is model-independent and a
  # binary model at raw 22 is degenerate rather than merely saturated: its
  # likelihood cannot be evaluated at all, which is a different failure.
  set.seed(5)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 8))
  d <- ctGenerate(gen, n.subjects = 20, Tpoints = 8, backend = "r")
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0)))
  m$pars$indvarying <- FALSE
  npar <- sum(is.na(m$pars$value) & !is.na(m$pars$param))
  fit <- suppressWarnings(suppressMessages(
    ctFit(d, m, backend = "julia", cores = 1, inits = rep(24, npar),
      optimcontrol = list(estonly = TRUE),
      backendcontrol = list(maxiter = 1))))
  expect_gt(max(abs(fit$estimate$raw)), 20)
  expect_false(isTRUE(fit$estimate$converged))
})
