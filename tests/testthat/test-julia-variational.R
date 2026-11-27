# ctFit(intoverpop = 'variational'): each unit's random effects integrated by
# a lower bound fitted at every evaluation (variational.jl in the engine,
# `.ctBackendVariational()`).
#
# What is pinned here: on a model where every unit's integrand is Gaussian the
# bound is the integral, so the fit is Laplace's (exact there); on a
# nonlinear one the fit reports the random-effects check's integral with the
# bound below it and the check's step from the estimate, its estimates are
# those of the refined fit, which integrates every unit numerically, and
# without the check it reports the bound; and the route is refused where it
# cannot apply.

.vf_linear_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.vf_linear_data <- function(nsubjects = 25, nobs = 6) {
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

# One random effect, on a `-log1p_exp` DRIFT, simulated here rather than by
# ctGenerate, whose draw stream moves.
.vf_nonlinear_model <- function() {
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  model$pars$indvarying <- model$pars$param %in% "drift"
  model
}

.vf_nonlinear_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
  set.seed(seed)
  drift <- -log1p(exp(-stats::rnorm(nsubjects, 1, 1)))
  rows <- lapply(seq_len(nsubjects), function(i) {
    a <- drift[i]; decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

.vf_cache <- new.env(parent = emptyenv())
.vf_fit <- function(which, route) {
  key <- paste(which, route)
  if (!exists(key, envir = .vf_cache, inherits = FALSE)) {
    data <- if (identical(which, "linear")) .vf_linear_data() else .vf_nonlinear_data()
    model <- if (identical(which, "linear")) .vf_linear_model() else .vf_nonlinear_model()
    set.seed(1)
    assign(key, suppressWarnings(suppressMessages(ctFit(data, model,
      backend = "julia", intoverpop = route, cores = 2))), envir = .vf_cache)
  }
  get(key, envir = .vf_cache)
}

test_that("on Gaussian units the variational fit is Laplace's", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("linear", "variational")
  lf <- .vf_fit("linear", "laplace")
  expect_identical(vf$args$resolved$intoverpop, "variational")
  expect_identical(vf$estimate$loglik_method, "variational")
  expect_equal(vf$estimate$raw, lf$estimate$raw, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik, lf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik_lowerbound, vf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$se, lf$estimate$se, tolerance = 0.02)
})

test_that("on a nonlinear model the bound lies below the reported log likelihood", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  qf <- .vf_fit("nonlinear", "quadrature")
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational")
  expect_true(is.finite(e$loglik) && is.finite(e$loglik_lowerbound))
  expect_lte(e$loglik_lowerbound, e$loglik + 4 * e$loglik_se)
  expect_equal(ctsem:::.ctFitOptimisedLoglik(vf), e$loglik_lowerbound)
  # The refined fit integrates every flagged unit numerically; the variational
  # estimate lies within a fifth of a standard error of it.
  expect_lt(max(abs(vf$estimate$raw - qf$estimate$raw) / qf$estimate$se), 0.2)
  expect_equal(e$loglik, qf$estimate$loglik, tolerance = 0.1)
  # The check's step from the estimate towards the integral's maximum, in the
  # bound's standard errors: small where the estimate is the refined fit's.
  expect_true(vf$variational$check$status %in% c("passed", "unresolved", "flagged"))
  expect_lt(vf$variational$check$max_shift, 0.2)
})

test_that("a variational fit's interval check and 'checked' run", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  check <- vf$uncertainty$details$normalcheck
  expect_null(check$failed)
  expect_gt(check$draws, 0)
  # 'checked' takes no coreset of an objective the engine cannot take one of.
  u <- suppressWarnings(suppressMessages(ctFitUncertainty(vf, "checked", cores = 2)))
  expect_true(all(is.finite(u$estimate$se)))
  # The joint sample's target is the Laplace objective of the same
  # specification; the bound's marginal is no sampler's.
  lf <- vf
  lf$model_spec$laplace$variational <- NULL
  joint <- ctsem:::.ctBackendSampleObjective(vf, list(marginal = FALSE))
  x <- as.numeric(vf$estimate$raw)
  module <- ctsem:::.ctJuliaModule(vf$model_spec$project)
  value <- function(o) as.numeric(ctsem:::.ctJuliaGet(module$ctsem_evaluate(o,
    ctsem:::.ctJuliaNumericVector(x), gradient = FALSE))$value)
  expect_equal(value(joint), value(ctsem:::.ctJuliaObjective(lf)))
  expect_error(ctsem:::.ctBackendSampleObjective(vf, list(marginal = TRUE)),
    "marginal cannot be sampled")
})

test_that("a variational fit is refined by quadrature from its own estimate", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  qf <- .vf_fit("nonlinear", "quadrature")
  rf <- suppressWarnings(suppressMessages(ctFitUncertainty(vf,
    intoverpop = "quadrature", cores = 2)))
  corr <- rf$laplace$correction
  expect_identical(corr$status, "continued")
  expect_true(corr$applied)
  # The climb started at the variational estimate, on the Laplace objective
  # with the units the check names against Laplace there flagged.
  expect_equal(corr$laplace_estimate, vf$estimate$raw)
  expect_null(rf$model_spec$laplace$variational)
  expect_gt(length(corr$check$carriers), 0)
  expect_true(all(corr$check$carriers %in% corr$flagged_units))
  expect_identical(rf$estimate$loglik_method, "quadrature")
  # It ends where ctFit(intoverpop = 'quadrature') ends from Laplace's
  # optimum, up to the units the screen flags at the two starts.
  expect_lt(max(abs(rf$estimate$raw - qf$estimate$raw) / qf$estimate$se), 0.1)
  expect_lt(abs(rf$estimate$loglik - qf$estimate$loglik), 0.1)
})

test_that("without the random-effects check a variational fit reports its bound", {
  skip_on_cran()
  skip_without_julia()
  set.seed(1)
  vf <- suppressWarnings(suppressMessages(ctFit(.vf_nonlinear_data(),
    .vf_nonlinear_model(), backend = "julia", intoverpop = "variational", cores = 2,
    optimcontrol = list(checks = FALSE))))
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational_bound")
  expect_identical(e$loglik, e$loglik_lowerbound)
  expect_null(vf$variational)
  expect_equal(e$loglik, .vf_fit("nonlinear", "variational")$estimate$loglik_lowerbound,
    tolerance = 1e-4)
})

test_that("leave-one-out and the score tests read a variational fit's bound", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  est <- as.numeric(vf$estimate$raw)
  # The units' terms are their bounds, and a subject with every row withheld
  # scores nothing.
  ut <- ctsem:::.ctBackendUnitTerms(vf$model_spec, est)
  expect_equal(sum(ut$unit_loglik), vf$estimate$loglik_lowerbound, tolerance = 1e-6)
  training <- vf$model_spec$data
  training[training$id == unique(training$id)[3], "Y1"] <- NA
  tf <- ctsem:::.ctFitReplaceData(vf, training)
  emptied <- ctsem:::.ctBackendUnitTerms(tf$model_spec, est, subjects = 3L)
  expect_equal(emptied$unit_loglik[3, 1], 0, tolerance = 1e-8)
  l <- suppressWarnings(suppressMessages(ctLOO(vf, folds = 4, cores = 2,
    casewiseApproximation = TRUE)))
  expect_true(all(is.finite(l$outsampleLogLikSubject)))
  expect_match(l$scoring, "variational bound")
  # A subject's row scores sum to its score on the bound, up to the inner
  # solves' tolerance.
  R <- ctsem:::.ctBackendRowScores(vf, est, diag(length(est)), 2)
  S <- ctsem:::.ctBackendScoreMatrix(vf, est, prior = FALSE)
  expect_lt(max(abs(rowsum(R, ctsem:::.ctFitRowSubject(vf)) - S)), 1e-3)
})

# Persons in studies: a manifest mean that varies at both levels (a Gaussian
# integrand, where the bound is the integral), and with `drift` a person-level
# drift through `-log1p_exp` beside it (not Gaussian).
.vf_nested_model <- function(drift = FALSE) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", silent = TRUE,
    id = c("id", "study"), manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), DRIFT = if (drift) "drift|-log1p_exp(-param)|TRUE" else
      matrix(-0.5), DIFFUSION = matrix(0.5), MANIFESTMEANS = matrix("mm"),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "mm" | (drift & m$pars$param %in% "drift")
  if (drift) m$pars$indvarying[m$pars$param %in% "drift"] <- "id"
  m
}

.vf_nested_data <- function(npersons = 40L, nstudies = 8L, ntimes = 8L, seed = 7L) {
  set.seed(seed)
  study <- rep(seq_len(nstudies), length.out = npersons)
  smean <- stats::rnorm(nstudies, 1, 0.7)
  rows <- lapply(seq_len(npersons), function(i) {
    a <- -log1p(exp(-stats::rnorm(1, 0.7, 0.6)))
    decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes)
    latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, study = study[i], time = seq_len(ntimes) - 1L,
      Y1 = smean[study[i]] + stats::rnorm(1, 0, 0.5) + latent +
        stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

.vf_nested_fit <- function(drift, route) {
  key <- paste("nested", drift, route)
  if (!exists(key, envir = .vf_cache, inherits = FALSE)) {
    set.seed(1)
    assign(key, suppressWarnings(suppressMessages(ctFit(.vf_nested_data(),
      .vf_nested_model(drift), backend = "julia", intoverpop = route, cores = 2))),
      envir = .vf_cache)
  }
  get(key, envir = .vf_cache)
}

test_that("persons in studies: the tree family is Laplace's on a Gaussian integrand", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_nested_fit(FALSE, "variational")
  lf <- .vf_nested_fit(FALSE, "laplace")
  expect_identical(vf$estimate$loglik_method, "variational")
  expect_equal(vf$estimate$raw, lf$estimate$raw, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik_lowerbound, vf$estimate$loglik, tolerance = 1e-3)
})

test_that("persons in studies with a random drift: the bound lies below the integral", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_nested_fit(TRUE, "variational")
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational")
  expect_true(is.finite(e$loglik) && is.finite(e$loglik_lowerbound))
  expect_lte(e$loglik_lowerbound, e$loglik + 4 * e$loglik_se + e$loglik_bound)
})

test_that("intoverpop = 'variational' is refused where it cannot apply", {
  skip_on_cran()
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "stan", intoverpop = "variational"))),
    "needs backend = 'julia'")
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "julia", intoverpop = "variational",
    optimize = FALSE))), "fits by optimizing")
})
