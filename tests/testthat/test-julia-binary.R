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

test_that("the binary adjoint is exact with more than one latent state", {
  skip_without_julia()
  # Every earlier gradient test used one latent state, and with one state the
  # covariance cotangent is a scalar and trivially symmetric. The binary
  # reverse assumed that symmetry in general, and `c = P lambda` breaks it --
  # so with two states the adjoint was about 1% wrong on DRIFT and DIFFUSION
  # and 12% wrong on the second state's variance, small enough to pass for
  # quadrature error. Forward mode is the comparison rather than a finite
  # difference because it is exact: a disagreement is then unambiguously the
  # reverse pass.
  d <- .jbin_data(nsubjects = 15, nobs = 12, nindicators = 2)
  d$v <- rep(stats::rnorm(15), each = 12)
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
    manifestNames = c("b1", "b2"), latentNames = "eta1",
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 2)))
  m$manifesttype[] <- 1L
  m$pars$indvarying <- FALSE
  # A random CINT expands the state by one under 'augmented', which is the
  # cheapest way to get a second state without changing the measurement model.
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  handle <- .jbin_spec(d, m)
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(3)
  at <- stats::rnorm(npar, 0, 0.3)
  adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
})

test_that("the two gradient methods agree with each other", {
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

# --- integration, found by exercising the package rather than the filter -----

.jbin_mixed <- function() {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 4,
    manifestNames = c("b1", "b2", "b3", "y1"), latentNames = "eta1",
    LAMBDA = matrix(1, 4, 1), MANIFESTMEANS = matrix(0, 4, 1),
    CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = diag(c(0, 0, 0, 1), 4)))
  m$manifesttype <- c(1L, 1L, 1L, 0L)
  m$pars$indvarying <- FALSE
  mm <- m$matrices; mm$MANIFESTVAR[4, 4] <- "mvar"; m$matrices <- mm
  m
}

# The same shape, as a *generating* model, with every matrix stated.
#
# `.jbin_mixed()` is a model to fit, so DRIFT, DIFFUSION and T0VAR are free in
# it. Handed to `ctGenerate()` those cells are filled from
# `.ctGenerateDefaults()`, and the data then moves whenever those defaults do --
# silently, under a `set.seed()` that reads as though it pinned everything. It
# already has: written on 2026-08-29 this generated a real process (DRIFT -0.5,
# DIFFUSION 1, T0VAR 1); under the defaults that followed, the latent is pinned
# at zero to within about 1e-3, which makes the "not degenerate" check below
# pass for exactly the wrong reason. Values match `.jbin_mixed_data()`,
# including a measurement sd of 0.5 on the gaussian indicator.
.jbin_mixed_genmodel <- function() {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 4,
    manifestNames = c("b1", "b2", "b3", "y1"), latentNames = "eta1",
    LAMBDA = matrix(1, 4, 1), MANIFESTMEANS = matrix(0, 4, 1),
    CINT = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8),
    MANIFESTVAR = diag(c(0, 0, 0, 0.5), 4)))
  m$manifesttype <- c(1L, 1L, 1L, 0L)
  m$pars$indvarying <- FALSE
  m
}

.jbin_mixed_data <- function(nsubjects = 30) {
  invlog <- function(x) exp(x) / (1 + exp(x))
  set.seed(11)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.001),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 12))
  d <- data.frame(ctGenerate(gen, n.subjects = nsubjects, Tpoints = 12,
    backend = "r"))
  for (i in 1:3) d[[paste0("b", i)]] <- stats::rbinom(nrow(d), 1, invlog(d$eta))
  d$y1 <- d$eta + stats::rnorm(nrow(d), 0, 0.5)
  d$eta <- NULL
  d
}

test_that("mixed binary and gaussian indicators fit together", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(
    ctFit(.jbin_mixed_data(), .jbin_mixed(), backend = "julia", cores = 2,
      optimcontrol = list(estonly = TRUE))))
  expect_true(isTRUE(fit$estimate$converged))
  est <- summary(fit)$popmeans
  expect_equal(unname(est["drift_eta1", "mean"]), -0.3, tolerance = 0.2)
  expect_equal(unname(est["diff_eta1", "mean"]), 0.8, tolerance = 0.35)
})

test_that("the functions that re-run the filter work on a binary fit", {
  skip_without_julia()
  # These all pass a `CTSEMKalmanTrace` rather than the adjoint tape, and the
  # binary recorder accepted it and then reached for a field it does not have.
  # One missing guard broke five functions at once, far from the cause.
  fit <- suppressWarnings(suppressMessages(
    ctFit(.jbin_mixed_data(), .jbin_mixed(), backend = "julia", cores = 2,
      optimcontrol = list(estonly = TRUE))))
  works <- function(expr) !inherits(
    try(suppressWarnings(suppressMessages(expr)), silent = TRUE), "try-error")
  expect_true(works(ctKalman(fit, subjects = 1)))
  expect_true(works(ctPredict(fit, subjects = 1)))
  expect_true(works(ctACFresiduals(fit)))
  expect_true(works(ctPostPredPlots(fit)))
  expect_true(works(ctLOO(fit, folds = 2, cores = 1)))
})

test_that("generation draws binary indicators as zeros and ones", {
  skip_without_julia()
  # The generate hook was passed only to the gaussian block, so binary columns
  # came back all NaN -- data that looks like a missing-data problem rather
  # than a bug in the generator.
  set.seed(4)
  d <- suppressWarnings(suppressMessages(
    ctGenerate(.jbin_mixed_genmodel(), n.subjects = 30, Tpoints = 10,
      backend = "julia")))
  for (nm in c("b1", "b2", "b3")) {
    expect_true(all(d[, nm] %in% c(0, 1)), info = nm)
    # Not degenerate: a mean-zero latent should give a mix.
    expect_gt(mean(d[, nm]), 0.2)
    expect_lt(mean(d[, nm]), 0.8)
  }
  expect_false(all(d[, "y1"] %in% c(0, 1)))
})

test_that("a binary model is routed away from the r generator", {
  # The r generator integrates a linear gaussian system and has no link, so it
  # produced continuous values for a manifest declared binary -- silently.
  m <- .jbin_mixed_genmodel()
  expect_warning(
    suppressMessages(ctGenerate(m, n.subjects = 3, Tpoints = 4,
      backend = "r")),
    "no measurement link")
})

test_that("state dependent measurement works alongside binary indicators", {
  skip_without_julia()
  d <- .jbin_mixed_data(nsubjects = 20)
  build <- function(lambda, means, mvar) {
    m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
      manifestNames = c("b1", "b2"), latentNames = "eta1",
      LAMBDA = lambda, MANIFESTMEANS = means, CINT = matrix(0),
      T0MEANS = matrix(0), MANIFESTVAR = mvar))
    m$manifesttype <- c(1L, 1L)
    m$pars$indvarying <- FALSE
    m
  }
  cases <- list(
    `state dependent LAMBDA` = build(matrix(c("1 + 0.2 * eta1", "1"), 2, 1),
      matrix(0, 2, 1), diag(0, 2)),
    `state dependent MANIFESTMEANS` = build(matrix(1, 2, 1),
      matrix(c("0.1 * eta1", "0"), 2, 1), diag(0, 2)),
    `state dependent MANIFESTVAR` = build(matrix(1, 2, 1), matrix(0, 2, 1),
      matrix(c("0.1 + 0.05 * eta1", "0", "0", "0"), 2, 2)))
  for (nm in names(cases)) {
    fit <- try(suppressWarnings(suppressMessages(
      ctFit(d, cases[[nm]], backend = "julia", cores = 2,
        optimcontrol = list(estonly = TRUE)))), silent = TRUE)
    expect_false(inherits(fit, "try-error"), info = nm)
    if (!inherits(fit, "try-error")) {
      expect_true(is.finite(fit$estimate$loglik), info = nm)
    }
  }
})

test_that("binary works with random effects and with laplace", {
  skip_without_julia()
  d <- .jbin_mixed_data(nsubjects = 25)
  d$age <- rep(stats::rnorm(25), each = 12)
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 3,
    manifestNames = c("b1", "b2", "b3"), latentNames = "eta1",
    LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 3),
    n.TIpred = 1, TIpredNames = "age"))
  m$manifesttype[] <- 1L
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  for (approach in c("augmented", "laplace")) {
    fit <- try(suppressWarnings(suppressMessages(
      ctFit(d, m, backend = "julia", cores = 2, intoverpop = approach,
        optimcontrol = list(estonly = TRUE)))), silent = TRUE)
    expect_false(inherits(fit, "try-error"), info = approach)
    if (!inherits(fit, "try-error")) {
      expect_true(is.finite(fit$estimate$loglik), info = approach)
    }
  }
})

test_that("ctModelLatex renders the link for binary indicators only", {
  m <- .jbin_mixed()
  txt <- paste(as.character(ctModelLatex(m, compile = FALSE)), collapse = "\n")
  # Three binary indicators get an inverse logit; the gaussian one does not.
  hits <- gregexpr("operatorname{logit}", txt, fixed = TRUE)[[1]]
  expect_equal(sum(hits > 0), 3L)
  expect_true(grepl("nu_", txt, fixed = TRUE))

  gaussian <- m
  gaussian$manifesttype <- rep(0L, 4)
  plain <- paste(as.character(ctModelLatex(gaussian, compile = FALSE)),
    collapse = "\n")
  expect_false(grepl("operatorname{logit}", plain, fixed = TRUE))
})
