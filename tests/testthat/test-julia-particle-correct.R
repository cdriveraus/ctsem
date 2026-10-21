# `ctFitUncertainty(fit, 'is', control = list(likelihood = 'particle'))`:
# importance sampling a fit's posterior against the particle-filter likelihood.
#
# The engine's suite proves the particle filter; test-julia-particle.R proves
# ctParticleLik(). What is left for here is that sampling against it changes
# the fit in the ways claimed and no others: on a linear model, where the
# filter is exact, the particle and filter likelihoods differ by noise and the
# posterior is the filter's; every slot the correction claims to write is
# written consistently; and each unsupported request refuses rather than
# answering.

skip_on_cran()

.pc_cache <- new.env(parent = emptyenv())
.pc_fit <- function() {
  if (!exists("fit", envir = .pc_cache, inherits = FALSE)) {
    model <- suppressWarnings(ctModel(
      type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
      DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
      MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)))
    # A stationary process, so both parameters are identified. A random walk
    # here sent the drift to its saturation plateau (raw 21, se 1e4), where a
    # posterior has no centre to correct toward. And priors, because even the
    # stationary data leave that plateau only 4.6 log units below the mode
    # 3000 raw units out: without a prior the posterior is improper there,
    # which uncertainty = 'is' reaches.
    set.seed(4)
    dat <- data.frame(id = rep(1:8, each = 6), time = rep(0:5, 8),
      Y1 = as.vector(replicate(8, {
        x <- numeric(6)
        x[1] <- rnorm(1)
        for (t in 2:6) x[t] <- 0.6 * x[t - 1] + rnorm(1, 0, 0.5)
        x + rnorm(6, 0, 0.3)
      })))
    fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
      cores = 1, inits = c(0.1, -0.2), verbose = 0, priors = TRUE,
      optimcontrol = list(finishsamples = 200))))
    assign("fit", fit, envir = .pc_cache)
  }
  get("fit", envir = .pc_cache, inherits = FALSE)
}

.particle_is <- function(fit, ..., finishsamples = 60) {
  suppressWarnings(ctFitUncertainty(fit, uncertainty = "is",
    finishsamples = finishsamples, cores = 1,
    control = list(likelihood = "particle", substeps = 1, seed = 3, ...)))
}

test_that("sampling against the particle filter records everything it did", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .pc_fit()
  set.seed(3)
  out <- .particle_is(fit, particles = 300, isitersize = 40, imisMaxIter = 1L,
    isESS = 20)
  pc <- out$particle_correction
  expect_s3_class(pc, "ctParticleCorrection")
  expect_equal(pc$draws, "imis")
  expect_match(pc$proposal, "imis")
  expect_equal(pc$scale, 1.5)
  expect_true(is.finite(pc$ess) && pc$ess > 0)
  expect_true(nrow(pc$evaluations) >= 40L)
  expect_equal(pc$ndraws, nrow(pc$evaluations))
  expect_equal(sum(pc$weights), 1)
  # The exponential particle step is exact for a linear model, so the particle
  # and filter likelihoods differ by Monte Carlo noise only, of the size the
  # filter reports. The median rather than the mean: draws in the proposal's
  # tails can sit where few particles degenerate, and carry no weight anyway.
  expect_lt(stats::median(abs(pc$evaluations$difference)),
    1.5 * stats::median(pc$evaluations$se))
  expect_equal(pc$evaluations$difference, pc$evaluations$particle - pc$evaluations$filter)
  # The check is ctParticleLik() at the original estimate, same settings.
  expect_equal(pc$check$fit_loglik, fit$estimate$loglik)
  expect_equal(pc$check$particles, 300L)
  expect_equal(pc$check$substeps, 1L)

  # Every slot the correction claims to write, written consistently.
  expect_equal(nrow(out$estimate$rawposterior), 60L)
  expect_equal(colnames(out$estimate$rawposterior), colnames(fit$estimate$rawposterior))
  expect_equal(unname(out$estimate$se), unname(sqrt(diag(out$estimate$cov))))
  expect_equal(names(out$estimate$se), names(fit$estimate$se))
  expect_equal(names(pc$shift_se), colnames(fit$estimate$rawposterior))
  expect_equal(pc$before$raw, as.numeric(fit$estimate$raw))
  expect_false(identical(out$transformedpars, fit$transformedpars))
  expect_false(is.null(out$uncertainty$imis))
  expect_equal(out$uncertainty$draws, "imis")
  settings <- out$uncertainty$settings
  expect_identical(settings$method, "is")
  expect_identical(settings$control$likelihood, "particle")
  expect_equal(settings$finishsamples, 60L)
  expect_equal(out$uncertainty$details$particle_correction$ess, pc$ess)
  # The filter's likelihood is left as the filter's.
  expect_equal(out$estimate$loglik, fit$estimate$loglik)

  expect_output(print(pc), "Particle-filter correction applied")
  expect_output(print(pc), "effective sample size")
  expect_match(suppressMessages(summary(out))$uncertaintyNote,
    "likelihood='particle'", fixed = TRUE)
  expect_match(suppressMessages(summary(fit))$uncertaintyNote, "ctFitUncertainty")
})

test_that("correctEstimate = FALSE leaves the point where it was", {
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .pc_fit()
  set.seed(2)
  out <- .particle_is(fit, particles = 200, isitersize = 30, imisMaxIter = 1L,
    isESS = 10, correctEstimate = FALSE, finishsamples = 40)
  expect_equal(out$estimate$raw, fit$estimate$raw)
  expect_false(out$particle_correction$corrected_estimate)
  expect_equal(nrow(out$estimate$rawposterior), 40L)
  expect_output(print(out$particle_correction), "estimate left at the original point")
})

test_that("unsupported requests are refused by name", {
  fit <- .pc_fit()
  expect_error(.particle_is(fit, particles = 1), "at least 2")
  expect_error(ctFitUncertainty(fit, "is",
    control = list(likelihood = "particle", substeps = 0)), "at least 1")
  expect_error(ctFitUncertainty(fit, "is",
    control = list(likelihood = "particle", nodes = 5)),
    "control$nodes is not used with likelihood = 'particle'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, "hessian",
    control = list(likelihood = "particle")), "applies to uncertainty = 'is'")
  bare <- fit
  bare$estimate$cov <- NULL
  expect_error(ctFitUncertainty(bare, "is", control = list(likelihood = "particle")),
    "covariance")
  # ctParticleLik's cores argument does not change the answer.
  a <- ctParticleLik(fit, particles = 200, substeps = 1, seed = 5)
  b <- ctParticleLik(fit, particles = 200, substeps = 1, seed = 5, cores = 2)
  expect_equal(a$loglik, b$loglik)
})

test_that("where the filter is exact, sampling against the particle filter targets the filter posterior", {
  # The particle likelihood equals the filter's on a linear model, so the
  # correction reduces to importance sampling against the filter posterior --
  # which is exactly what uncertainty = 'is' computes against the fit's own
  # objective. The two must agree on the posterior's centre and width.
  skip_if_not_installed('mvtnorm'); skip_if_not_installed('diagis')
  fit <- .pc_fit()
  set.seed(6)
  reference <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "is", finishsamples = 400, verbose = 0,
    control = list(isitersize = 300, isESS = 60))))
  # uncertainty = 'is' on the julia route records the tail of its weights; on
  # a linear model the posterior has the proposal's shape and the tail is light.
  k <- reference$uncertainty$details$importance_sampling$pareto_k
  expect_true(is.finite(k))
  expect_lt(k, 0.7)
  set.seed(7)
  corrected <- .particle_is(fit, particles = 500, isitersize = 200, isESS = 60,
    finishsamples = 400)
  se <- as.numeric(fit$estimate$se)
  ref_mean <- colMeans(reference$estimate$rawposterior)
  ref_sd <- apply(reference$estimate$rawposterior, 2, stats::sd)
  expect_lt(max(abs(as.numeric(corrected$estimate$raw) - ref_mean) / se), 0.5)
  ratio <- as.numeric(corrected$estimate$se) / as.numeric(ref_sd)
  expect_true(all(ratio > 0.7 & ratio < 1.4))
})
