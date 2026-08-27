# ctSample(): Hamiltonian sampling of a julia Laplace fit.
#
# The statistics of the sampler are tested in the engine suite, against
# posteriors known in closed form. What is tested here is the R side of it: that
# a sampled fit is a `ctJuliaFit` the existing summary machinery can read, that
# the diagnostics arrive intact, and that the failure modes are refused rather
# than half-performed.
#
# Deliberately small and short. A sample that mixes well enough to assert
# R-hat on takes minutes, which does not belong in a unit suite; the engine
# tests carry the correctness argument.

.sample_fixture <- function(nsub = 12L, tp = 6L, seed = 31L) {
  set.seed(seed)
  data <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    intercept <- stats::rnorm(1, 0, 0.8)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(tp)
    for (t in seq_len(tp)) {
      state <- 0.75 * state + stats::rnorm(1, 0, 0.4)
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(tp) - 1, Y1 = y)
  }))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE
  suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    cores = 1, intoverpop = "laplace", priors = TRUE,
    optimcontrol = list(finishsamples = 20))))
}

test_that("a sampled fit carries draws the summary machinery can read", {
  skip_on_cran()
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 80, draws = 80, cores = 1)))

  expect_s3_class(sampled, "ctJuliaFit")
  # The draws are where an optimised fit's normal-approximation draws live, so
  # everything downstream reads them without knowing which produced them.
  expect_equal(dim(sampled$estimate$rawposterior), c(160L, npar))
  expect_true(all(is.finite(sampled$estimate$rawposterior)))
  expect_equal(colnames(sampled$estimate$rawposterior),
    ctsem:::.ctBackendRawParameterNames(fit, npar))

  # The point estimate becomes the posterior mean, and the Laplace estimate it
  # started from is kept rather than overwritten.
  expect_equal(sampled$estimate$raw,
    as.numeric(colMeans(sampled$estimate$rawposterior)))
  expect_equal(sampled$estimate$laplace_raw, fit$estimate$raw)
  expect_equal(sampled$estimate$se,
    sqrt(diag(stats::cov(sampled$estimate$rawposterior))))

  # And the cached constrained draws describe the new draws, not the old ones.
  expect_false(is.null(sampled$transformedpars))
  expect_equal(dim(sampled$transformedpars$samples),
    dim(sampled$estimate$rawposterior))
  expect_no_error(suppressWarnings(summary(sampled)))
})

test_that("the diagnostics come back per parameter and per chain", {
  skip_on_cran()
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 80, draws = 80, cores = 1)))
  diagnostics <- sampled$sample

  expect_s3_class(diagnostics, "ctSampleDiagnostics")
  expect_length(diagnostics$rhat, npar)
  expect_length(diagnostics$ess, npar)
  expect_equal(names(diagnostics$rhat), colnames(sampled$estimate$rawposterior))
  expect_true(all(diagnostics$rhat > 0.99, na.rm = TRUE))
  expect_true(all(diagnostics$ess > 0, na.rm = TRUE))
  expect_length(diagnostics$stepsize, 2L)
  expect_true(all(diagnostics$stepsize > 0))
  expect_length(diagnostics$ebfmi, 2L)
  expect_true(diagnostics$divergent >= 0L)
  expect_length(diagnostics$accept, 160L)
  # The Laplace estimate the chains started from.
  expect_equal(diagnostics$start, fit$estimate$raw)
  expect_output(print(diagnostics), "ctsem Hamiltonian sample")
})

test_that("the effects come back summarised, or in full when asked for", {
  skip_on_cran()
  skip_without_julia()
  fit <- .sample_fixture()
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)

  # By default only the summary crosses the bridge: the draws themselves are
  # nsubjects x neffects x chains x draws numbers, and the bridge moves about
  # 1 MB/s.
  lean <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 1, warmup = 60, draws = 60, cores = 1)))
  expect_null(lean$sample$effects)
  expect_length(lean$sample$effect_mean, neffects)
  expect_length(lean$sample$effect_sd, neffects)
  expect_true(all(lean$sample$effect_sd > 0))

  full <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 1, warmup = 60, draws = 60, cores = 1,
      saveEffects = TRUE)))
  expect_equal(dim(full$sample$effects), c(60L, neffects))
  expect_true(all(is.finite(full$sample$effects)))
})

test_that("ctSample refuses what it cannot sample", {
  skip_on_cran()
  skip_without_julia()
  expect_error(ctSample(list()), "ctFit\\(backend='julia'\\)")

  # The augmented route carries the effects in the state, so there is no
  # separate posterior over them and no Laplace curvature to build a metric
  # from. Saying so beats sampling something else.
  set.seed(4)
  data <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * 0.5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  augmented <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_error(ctSample(augmented), "intoverpop='laplace'")
})
