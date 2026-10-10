# ctScoreTest(): score tests for relaxations of a fitted model.
#
# What is asserted is what the tests rest on, as properties rather than as
# plausible output: the per-subject pieces reproduce the engine's own
# derivatives (a covariate's score is the gradient in the TI effect it would
# add; the summed per-subject Hessians are the fit's Hessian); the model
# extended with state-dependent terms reproduces the fit at zero; and, on data
# built to contain each relaxation, the test finds it. Size and power were
# measured by simulation (CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md);
# repeating that here would cost hundreds of fits.

skip_without_julia()
skip_on_32bit()

# One latent, irregular intervals, a random manifest mean, and the three
# departures switched on one at a time. Simulated here, not by ctGenerate,
# whose draw stream moves under unrelated commits.
.score_data <- function(scenario = 'null', nsub = 80L, nobs = 25L, seed = 3L) {
  set.seed(seed)
  sigma <- 0.8; r <- 0.5; dtsub <- 0.02
  z <- stats::rnorm(nsub); w <- stats::rnorm(nsub)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    a <- -0.5; b <- 0
    if (scenario == 'het') a <- -0.5 * exp(0.6 * stats::rnorm(1))
    if (scenario == 'cov') a <- -0.5 * exp(0.3 * z[i])
    if (scenario == 'sdep') b <- 0.4
    mm <- stats::rnorm(1) + if (scenario == 'cov') 0.8 * z[i] else 0
    times <- cumsum(c(0, stats::rexp(nobs - 1L, 1)))
    x <- stats::rnorm(1, 0, sigma / sqrt(-2 * a))
    xs <- numeric(nobs); xs[1] <- x
    for (t in 2:nobs) {
      n <- max(1L, ceiling((times[t] - times[t - 1]) / dtsub))
      h <- (times[t] - times[t - 1]) / n
      for (s in seq_len(n)) x <- x + (a + b * tanh(x)) * x * h + sigma * sqrt(h) * stats::rnorm(1)
      xs[t] <- x
    }
    data.frame(id = i, time = times, Y1 = xs + mm + stats::rnorm(nobs, 0, r), z = z[i], w = w[i])
  }))
}

.score_model <- function() suppressMessages(ctModel(type = 'ct',
  manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1)))

.score_fit <- function(scenario, seed = 3L) fit_cached(paste('scoretest', scenario, seed), {
  d <- .score_data(scenario, seed = seed)
  list(data = d, fit = suppressMessages(ctFit(d[, c('id', 'time', 'Y1')], .score_model(),
    backend = 'julia', cores = 1, verbose = 0)))
})

test_that("the per-subject pieces reproduce the engine's derivatives", {
  x <- .score_fit('cov')
  fit <- x$fit
  est <- fit$estimate$raw
  ptab <- ctsem:::.ctScoreParameterTable(fit)
  model <- which(ptab$kind == 'model')
  pieces <- ctsem:::.ctScorePieces(fit, est, rows = model, second = integer(), cores = 1)
  # Each subject's Hessian row is differenced along the tested parameter alone
  # (its symmetry gives the row), so summed over subjects the rows must be the
  # fit's exact Hessian's.
  exact <- ctsem:::.ctBackendHessian(fit, est)
  for (a in seq_along(model)) expect_equal(colSums(pieces$H[, a, ]),
    as.numeric(exact[model[a], ]), tolerance = 1e-5)
  # A covariate's score is the gradient, at zero, in the TI effect it would add:
  # the same model with z as a TI predictor, effects on every parameter.
  withz <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1',
    latentNames = 'eta1', LAMBDA = matrix(1), TIpredNames = 'z'))
  args <- ctsem:::.ctFitReplayArgs(fit)
  args$model <- withz; args$fit <- FALSE; args$inits <- NULL
  args$datalong <- x$data[, c('id', 'time', 'Y1', 'z')]
  args$intoverpop <- ctsem:::.ctBackendIntOverPop(fit$model_spec)
  spec <- suppressMessages(do.call(ctFit, args))
  names2 <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), ctsem:::.ctBackendNpar(spec))
  at <- setNames(numeric(length(names2)), names2)
  at[ptab$name] <- est
  gradient <- setNames(ctJuliaEvaluate(spec, at, cores = 1)$gradient, names2)
  z <- unique(x$data[, c('id', 'z')])$z
  score <- setNames(as.numeric(crossprod(z, pieces$S[, model])), ptab$name[model])
  expect_equal(unname(gradient[paste0('rawtipredeffect_', names(score), '_z')]),
    unname(score), tolerance = 1e-6)
})

test_that("a covariate the model lacks is found and attributed; an unrelated one is not", {
  x <- .score_fit('cov')
  st <- suppressMessages(ctScoreTest(x$fit, covariates = unique(x$data[, c('id', 'z', 'w')]),
    tests = 'covariates', cores = 1, verbose = FALSE))
  expect_s3_class(st, 'ctScoreTest')
  expect_lt(st$covariates$z$joint$p, 1e-3)
  expect_gt(st$covariates$w$joint$p, 1e-3)
  par <- st$covariates$z$parameters
  mm <- par[par$parameter == 'mm_Y1', ]
  expect_lt(mm$p, 1e-3)
  # The held one-step is near what generated the data, 0.8 on the mean.
  expect_equal(mm$estimate_natural, 0.8, tolerance = 0.3)
  expect_output(print(st), 'Covariate z')
})

test_that("a parameter that varies over subjects is found", {
  x <- .score_fit('het')
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'heterogeneity', cores = 1,
    verbose = FALSE))
  h <- st$heterogeneity
  expect_setequal(h$parameter, c('drift_eta1', 'diff_eta1', 'mvarY1'))
  expect_lt(h$p[h$parameter == 'drift_eta1'], 0.01)
  expect_gt(h$sd_natural[h$parameter == 'drift_eta1'], 0)
})

test_that("a drift that depends on the state is found by the extended model", {
  # A strong dependence (the drift moves by 0.4 across the state's range). On
  # such data the fitted model's curvature at the null can be indefinite in the
  # term's direction: the one-step estimate is then NA and the model-based
  # variance unusable, which is what the robust variance is for -- over seeds
  # 3 to 5 it found the term at p = 0.03, 0.0004, 0.001.
  x <- .score_fit('sdep', seed = 4L)
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'statedependence', cores = 1,
    verbose = FALSE))
  terms <- st$statedependence$terms
  expect_equal(terms$term, 'eta1*eta1')
  expect_lt(terms$p, 0.01)
  expect_true(is.na(terms$estimate) || terms$estimate > 0)
  # Degree 3 adds the cubic term and keeps the quadratic one.
  st3 <- suppressMessages(ctScoreTest(x$fit, tests = 'statedependence', order = 3,
    cores = 1, verbose = FALSE))
  expect_setequal(st3$statedependence$terms$term, c('eta1*eta1', 'eta1*eta1*eta1'))
})

test_that("results are reproducible and leave the caller's random numbers alone", {
  x <- .score_fit('cov')
  set.seed(11)
  before <- .Random.seed
  a <- suppressMessages(ctScoreTest(x$fit, covariates = unique(x$data[, c('id', 'w')]),
    tests = 'covariates', cores = 1, verbose = FALSE, nboot = 199))
  expect_identical(.Random.seed, before)
  b <- suppressMessages(ctScoreTest(x$fit, covariates = unique(x$data[, c('id', 'w')]),
    tests = 'covariates', cores = 1, verbose = FALSE, nboot = 199))
  expect_identical(a$covariates$w$joint$p, b$covariates$w$joint$p)
})

test_that("what cannot be tested is refused by name", {
  x <- .score_fit('cov')
  expect_error(ctScoreTest(structure(list(), class = 'ctStanFit')), "backend = 'julia'")
  expect_error(ctScoreTest(x$fit, tests = 'covariates'), 'no covariates')
  expect_error(ctScoreTest(x$fit, covariates = data.frame(z = 1:3),
    tests = 'covariates', verbose = FALSE), 'one row per subject')
  expect_error(ctScoreTest(x$fit, parameters = 'nosuch', tests = 'heterogeneity',
    verbose = FALSE), 'Unknown parameters')
})

test_that("joint heterogeneity: a dimension test beside the per-parameter ones", {
  x <- .score_fit('het')
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'heterogeneity', cores = 1,
    verbose = FALSE, nboot = 499))
  d <- st$dimension
  expect_true(is.finite(d$lambda))
  expect_true(d$p > 0 && d$p <= 1)
  expect_setequal(names(d$direction), st$heterogeneity$parameter)
  # The largest eigenvalue of the standardised score matrix is at least its
  # largest diagonal, the per-parameter z.
  expect_gte(d$lambda, max(st$heterogeneity$z) - 1e-8)
  expect_true('class' %in% names(st$heterogeneity))
  expect_output(print(st), 'new dimension')
})

test_that("a threshold covariate finds the step and reports where", {
  x <- .score_fit('cov')
  covs <- unique(x$data[, c('id', 'z')])
  covs$step <- as.numeric(covs$z > 0.3)
  st <- suppressMessages(ctScoreTest(x$fit, covariates = covs[, c('id', 'z')],
    tests = 'covariates', form = 'threshold', cores = 1, verbose = FALSE, nboot = 499))
  j <- st$covariates$z$joint
  expect_identical(j$form, 'threshold')
  expect_true(is.finite(j$cut))
  expect_lt(j$p, 0.01)
  expect_error(ctScoreTest(x$fit, covariates = data.frame(id = covs$id,
    f = letters[seq_len(nrow(covs)) %% 3 + 1]), tests = 'covariates', form = 'threshold',
    verbose = FALSE), 'not numeric')
})

test_that("a family that fails is reported and the others stand", {
  x <- .score_fit('het')
  local_mocked_bindings(.ctScoreExtended = function(...) stop('cannot extend this model'))
  st <- suppressWarnings(suppressMessages(ctScoreTest(x$fit,
    tests = c('heterogeneity', 'statedependence'), cores = 1, verbose = FALSE, nboot = 199)))
  expect_false(is.null(st$heterogeneity))
  expect_null(st$statedependence)
  expect_match(st$info$errors[['statedependence']], 'cannot extend')
  expect_output(print(st), 'Not tested')
})

test_that("with priors on every parameter, the extended terms are tested on the likelihood", {
  x <- .score_fit('sdep', seed = 4L)
  args <- ctsem:::.ctFitReplayArgs(x$fit)
  args$model <- .score_model(); args$datalong <- x$data[, c('id', 'time', 'Y1')]
  args$priors <- TRUE; args$verbose <- 0
  fitp <- suppressMessages(do.call(ctFit, args))
  a <- suppressMessages(ctScoreTest(x$fit, tests = 'statedependence', cores = 1, verbose = FALSE))
  b <- suppressMessages(ctScoreTest(fitp, tests = 'statedependence', cores = 1, verbose = FALSE))
  expect_null(b$info$errors)
  # The prior moves the estimates a little; the statistic follows, no more.
  expect_equal(b$statedependence$terms$LM, a$statedependence$terms$LM, tolerance = 0.1)
})
