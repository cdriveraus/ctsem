# ctParticleCorrect(): importance-sampling correction of a fit's posterior
# draws against the particle-filter likelihood.
#
# The engine's suite proves the particle filter; test-julia-particle.R proves
# ctParticleLik(). What is left for here is the reweighting arithmetic and that
# applying it changes the fit in the ways claimed and no others: a linear model,
# where the filter is exact, comes back with nearly uniform weights and a
# barely moved estimate; every slot the correction claims to write is written
# consistently; and each unsupported request refuses rather than answering.

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
    # posterior has no centre to correct toward.
    set.seed(4)
    dat <- data.frame(id = rep(1:8, each = 6), time = rep(0:5, 8),
      Y1 = as.vector(replicate(8, {
        x <- numeric(6)
        x[1] <- rnorm(1)
        for (t in 2:6) x[t] <- 0.6 * x[t - 1] + rnorm(1, 0, 0.5)
        x + rnorm(6, 0, 0.3)
      })))
    fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
      cores = 1, inits = c(0.1, -0.2), verbose = 0,
      optimcontrol = list(finishsamples = 200))))
    assign("fit", fit, envir = .pc_cache)
  }
  get("fit", envir = .pc_cache, inherits = FALSE)
}

test_that("reweighting a linear fit changes little and records everything it did", {
  fit <- .pc_fit()
  expect_equal(fit$uncertainty$draws, "normal")
  set.seed(1)
  out <- ctParticleCorrect(fit, particles = 500, substeps = 2, seed = 3)
  pc <- out$particle_correction
  expect_s3_class(pc, "ctParticleCorrection")
  expect_equal(pc$draws, "reweight")
  expect_match(pc$proposal, "fresh normal draws")
  expect_equal(pc$scale, 1.5)
  expect_equal(pc$ndraws, 200L)
  expect_equal(dim(pc$draws_evaluated), dim(fit$estimate$rawposterior))
  expect_equal(colnames(pc$draws_evaluated), colnames(fit$estimate$rawposterior))
  # Fresh draws, wider than the fit's: not the fit's own rows.
  expect_false(isTRUE(all.equal(pc$draws_evaluated, fit$estimate$rawposterior)))
  expect_equal(length(pc$weights), 200L)
  expect_equal(sum(pc$weights), 1)
  expect_true(pc$ess > 20 && pc$ess <= 200)
  expect_equal(nrow(pc$evaluations), 200L)
  # The exponential particle step is exact for a linear model, so the particle
  # and filter likelihoods differ by Monte Carlo noise only, of the size the
  # filter reports. The median rather than the mean: draws in the proposal's
  # tails can sit where 500 particles degenerate, and carry no weight anyway.
  expect_lt(stats::median(abs(pc$evaluations$difference)),
    1.5 * stats::median(pc$evaluations$se))
  expect_equal(pc$evaluations$difference, pc$evaluations$particle - pc$evaluations$filter)
  # The check is ctParticleLik() at the original estimate, same settings.
  expect_equal(pc$check$fit_loglik, fit$estimate$loglik)
  expect_equal(pc$check$particles, 500L)
  expect_equal(pc$check$substeps, 2L)

  # Every slot the correction claims to write, written consistently.
  expect_equal(dim(out$estimate$rawposterior), dim(fit$estimate$rawposterior))
  expect_equal(colnames(out$estimate$rawposterior), colnames(fit$estimate$rawposterior))
  expect_equal(unname(out$estimate$se), unname(sqrt(diag(out$estimate$cov))))
  expect_equal(names(out$estimate$se), names(fit$estimate$se))
  expect_equal(unname(out$estimate$raw),
    unname(as.numeric(colSums(pc$weights * pc$draws_evaluated))))
  expect_equal(names(pc$shift_se), colnames(fit$estimate$rawposterior))
  # The shift here is the normal approximation being corrected toward a skewed
  # six-subject posterior, not the filter; the last test pins it against the
  # package's own importance sampler. This only says it is not absurd.
  expect_lt(pc$largest_shift_se, 3)
  expect_true(all(pc$width_ratio > 0.5 & pc$width_ratio < 2))
  expect_equal(pc$before$raw, as.numeric(fit$estimate$raw))
  expect_false(identical(out$transformedpars, fit$transformedpars))
  expect_equal(out$uncertainty$draws, "particle")
  expect_equal(out$uncertainty$settings$draws, "particle")
  expect_equal(out$uncertainty$settings$finishsamples, 200L)
  expect_equal(out$uncertainty$details$particle_correction$ess, pc$ess)
  # The filter's likelihood is left as the filter's.
  expect_equal(out$estimate$loglik, fit$estimate$loglik)

  expect_output(print(pc), "Particle-filter correction applied")
  expect_output(print(pc), "effective sample size")
  s <- suppressMessages(summary(out))
  expect_match(s$uncertaintyNote, "ctParticleCorrect")
  expect_match(suppressMessages(summary(fit))$uncertaintyNote, "ctOptimUncertainty")
})

test_that("a corrected fit reweights again as a posterior-distributed one", {
  fit <- .pc_fit()
  set.seed(2)
  once <- suppressWarnings(ctParticleCorrect(fit, particles = 300, substeps = 1, seed = 3,
    nsamples = 60, finishsamples = 80, correct_estimate = FALSE))
  pc <- once$particle_correction
  expect_equal(once$estimate$raw, fit$estimate$raw)
  expect_false(pc$corrected_estimate)
  expect_equal(pc$ndraws, 60L)
  expect_equal(nrow(pc$draws_evaluated), 60L)
  expect_equal(nrow(once$estimate$rawposterior), 80L)
  expect_equal(pc$finishsamples, 80L)
  expect_output(print(pc), "estimate left at the original point")
  # Its draws are now a resample of the particle posterior, so a second pass
  # reuses them and weights by the likelihood ratio alone, not against a normal.
  twice <- suppressWarnings(ctParticleCorrect(once, particles = 300, substeps = 1, seed = 3,
    nsamples = 50, finishsamples = 40))
  pc2 <- twice$particle_correction
  expect_match(pc2$proposal, "filter posterior")
  expect_equal(nrow(pc2$draws_evaluated), 50L)
  expect_true(all(pc2$draws_evaluated %in% once$estimate$rawposterior))
  expect_equal(nrow(twice$estimate$rawposterior), 40L)
})

test_that("draws = 'imis' samples against the particle posterior", {
  fit <- .pc_fit()
  set.seed(3)
  out <- suppressWarnings(ctParticleCorrect(fit, draws = "imis", particles = 300,
    substeps = 1, seed = 3, nbatch = 40, maxiter = 1L, target_ess = 20,
    finishsamples = 60))
  pc <- out$particle_correction
  expect_equal(pc$draws, "imis")
  expect_match(pc$proposal, "imis")
  expect_true(is.finite(pc$ess) && pc$ess > 0)
  expect_equal(nrow(out$estimate$rawposterior), 60L)
  expect_equal(ncol(out$estimate$rawposterior), length(fit$estimate$raw))
  expect_true(nrow(pc$evaluations) >= 40L)
  expect_equal(pc$ndraws, nrow(pc$evaluations))
  expect_equal(out$uncertainty$draws, "particle")
  expect_false(is.null(out$uncertainty$imis))
  expect_equal(unname(out$estimate$se), unname(sqrt(diag(out$estimate$cov))))
})

test_that("unsupported requests are refused by name", {
  fit <- .pc_fit()
  expect_error(ctParticleCorrect(list()), "backend = 'julia'")
  expect_error(ctParticleCorrect(fit, particles = 1), "at least 2")
  expect_error(ctParticleCorrect(fit, substeps = 0), "at least 1")
  bare <- fit
  bare$estimate$rawposterior <- NULL
  expect_error(ctParticleCorrect(bare), "no posterior draws")
  unknown <- fit
  unknown$uncertainty$draws <- "mystery"
  unknown$uncertainty$settings$draws <- "mystery"
  expect_error(ctParticleCorrect(unknown), "not recognised")
  # ctParticleLik gained the same cores argument; the count does not change the answer.
  a <- ctParticleLik(fit, particles = 200, substeps = 1, seed = 5)
  b <- ctParticleLik(fit, particles = 200, substeps = 1, seed = 5, cores = 2)
  expect_equal(a$loglik, b$loglik)
})

test_that("where the filter is exact, reweighting normal draws targets the filter posterior", {
  # The particle likelihood equals the filter's on a linear model, so the
  # correction reduces to importance sampling the normal draws against the
  # filter posterior -- which is exactly what ctOptimUncertainty('is') computes
  # by another route. The two must agree on the posterior's centre and width.
  fit <- .pc_fit()
  set.seed(6)
  reference <- suppressWarnings(suppressMessages(ctOptimUncertainty(fit,
    uncertainty = "is", finishsamples = 400, verbose = 0,
    control = list(isitersize = 300, isESS = 60))))
  set.seed(7)
  corrected <- suppressWarnings(ctParticleCorrect(fit, particles = 500, substeps = 1, seed = 3))
  se <- as.numeric(fit$estimate$se)
  ref_mean <- colMeans(reference$estimate$rawposterior)
  ref_sd <- apply(reference$estimate$rawposterior, 2, stats::sd)
  expect_lt(max(abs(as.numeric(corrected$estimate$raw) - ref_mean) / se), 0.5)
  ratio <- as.numeric(corrected$estimate$se) / as.numeric(ref_sd)
  expect_true(all(ratio > 0.7 & ratio < 1.4))
})
