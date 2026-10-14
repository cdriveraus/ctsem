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
      # 'vol': the innovation sd rises with the state, exp(0.5 x) within two
      # sd, bounded beyond.
      g <- if (scenario == 'vol') exp(0.5 * 2 * tanh(x / 2)) else 1
      for (s in seq_len(n)) x <- x + (a + b * tanh(x)) * x * h + g * sigma * sqrt(h) * stats::rnorm(1)
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

test_that("innovations whose size depends on the state are found", {
  # The extended model reproduces the fit at zero (.ctScoreExtended stops
  # otherwise), so reaching the table is part of the assertion.
  x <- .score_fit('vol', seed = 4L)
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'volatility', cores = 1, verbose = FALSE))
  terms <- st$volatility$terms
  expect_equal(nrow(terms), 1L)
  expect_equal(terms$by, 'eta1')
  expect_lt(terms$p, 0.01)
  expect_true(is.na(terms$estimate) || terms$estimate > 0)
  expect_equal(st$volatility$joint$equation, c('eta1', 'all'))
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
  expect_gte(d$lambda, max(st$heterogeneity$z, na.rm = TRUE) - 1e-8)
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

test_that("row scores sum to the subjects' scores, on the filter and the Laplace route", {
  # Identities, not plausibility: a subject's rows' log predictive densities
  # sum to its log likelihood, so their derivatives must sum to its score.
  x <- .score_fit('het')
  fit <- x$fit
  est <- ctsem:::.ctFitRawEstimate(fit)
  ptab <- ctsem:::.ctScoreParameterTable(fit)
  model <- which(ptab$kind == 'model')
  D <- diag(length(est))[, model, drop = FALSE]
  R <- ctsem:::.ctBackendRowScores(fit, est, D, cores = 1)
  unit <- ctsem:::.ctAdequacyErrors(fit)$unit
  S <- ctsem:::.ctBackendScoreMatrix(fit, est)
  expect_equal(unname(rowsum(R, unit)), unname(S[, model, drop = FALSE]), tolerance = 1e-8)

  # The Laplace route: rows from the sequential predictions, each the
  # difference of the subject's terms on its record up to the row and before.
  lfit <- fit_cached('scoretest het laplace', suppressMessages(ctFit(
    x$data[, c('id', 'time', 'Y1')], .score_model(), backend = 'julia',
    intoverpop = 'laplace', cores = 1, verbose = 0)))
  lest <- ctsem:::.ctFitRawEstimate(lfit)
  lptab <- ctsem:::.ctScoreParameterTable(lfit)
  lmodel <- which(lptab$kind %in% c('model', 'popsd'))
  LD <- diag(length(lest))[, lmodel, drop = FALSE]
  LR <- ctsem:::.ctBackendRowScores(lfit, lest, LD, cores = 1)
  lunit <- ctsem:::.ctAdequacyErrors(lfit)$unit
  LS <- ctsem:::.ctBackendScoreMatrix(lfit, lest)
  expect_equal(unname(rowsum(LR, lunit)), unname(LS[, lmodel, drop = FALSE]), tolerance = 1e-6)

  # Each unit's posterior, as the quadrature's own rule. The effects here are
  # the intercept and the initial mean, which shift the means affinely, so
  # every posterior is exactly normal: three nodes a dimension and five give
  # the same marginal, mean and covariance, and the marginal is the
  # quadrature's unit term.
  d <- ctsem:::.ctBackendUnitDims(lfit)[1]
  expect_equal(d, 2L)
  p3 <- ctsem:::.ctBackendPosteriorNodes(lfit, lest, nodes = 3L, cores = 1)
  p5 <- ctsem:::.ctBackendPosteriorNodes(lfit, lest, nodes = 5L, cores = 1)
  expect_equal(dim(p3$points), c(nrow(p3$weights), 9L, 2L))
  expect_equal(rowSums(p5$weights), rep(1, nrow(p5$weights)))
  postmoments <- function(p) {
    w <- p$weights; a <- p$points[, , 1]; b <- p$points[, , 2]
    m1 <- rowSums(w * a); m2 <- rowSums(w * b)
    cbind(m1, m2, rowSums(w * a^2) - m1^2, rowSums(w * b^2) - m2^2,
      rowSums(w * a * b) - m1 * m2)
  }
  expect_equal(postmoments(p3), postmoments(p5), tolerance = 1e-6)
  expect_equal(p3$logmarginal, p5$logmarginal, tolerance = 1e-8)
  module <- ctsem:::.ctJuliaModule(lfit$model_spec$project)
  quad <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_quadrature(
    ctsem:::.ctJuliaObjective(lfit, hybrid = TRUE), ctsem:::.ctJuliaNumericVector(lest),
    nodes = 3L, contributions = TRUE))
  expect_equal(p3$logmarginal, as.numeric(quad$subject_loglik), tolerance = 1e-10)

  # The effect-distribution family runs there, and is quiet on normal effects.
  ad <- suppressMessages(ctFitAdequacy(lfit, families = 'effects', cores = 1, nboot = 199))
  expect_identical(ad$families$family, 'effects')
  expect_gt(ad$families$p, 0.01)

  # On this route the estimation correction can take its differences from
  # plug-in passes, and says that it did.
  ap <- suppressMessages(ctFitAdequacy(lfit, families = c('residuals', 'distribution'),
    projection = 'plugin', cores = 1, nboot = 199))
  expect_identical(ap$projection, 'plugin')
  expect_true(all(ap$families$p > 0 & ap$families$p <= 1))
  # The effects here shift the means, so each row's predictive is exactly
  # Gaussian and the row distances are tested; where effects reach the
  # dynamics they would not be (.ctAdequacyApproximate).
  expect_false(ctsem:::.ctAdequacyApproximate(lfit))
  expect_true('distance' %in% ap$statistics$statistic)

  # The model's exact sequential predictive by importance sampling over the
  # effects (what the error families use where the delta method is
  # approximate). Here the delta method is exact, so the two agree to the
  # sampler's error, and each subject's rows sum to its marginal likelihood.
  obj <- ctsem:::.ctJuliaObjective(lfit, hybrid = TRUE)
  old <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential(obj,
    ctsem:::.ctJuliaNumericVector(lest)))
  mix <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential_mixture(obj,
    ctsem:::.ctJuliaNumericVector(lest), particles = 512L))
  e <- (mix$y[1, , 1] - old$y[1, , 1]) / sqrt(old$ycov[1, , 1, 1])
  expect_lt(sqrt(mean(e^2)), 0.02)
  expect_lt(sqrt(mean((mix$ycov[1, , 1, 1] / old$ycov[1, , 1, 1] - 1)^2)), 0.04)
  expect_gt(min(mix$ess), 20)
  q7 <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_quadrature(obj,
    ctsem:::.ctJuliaNumericVector(lest), nodes = 7L, contributions = TRUE))
  expect_lt(abs(mean(tapply(mix$llrow, mix$subject, sum) - as.numeric(q7$subject_loglik))), 0.02)
})

test_that("stability along a predicted state finds a drift that depends on it", {
  x <- .score_fit('sdep', seed = 4L)
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'instability', cores = 1,
    verbose = FALSE, nboot = 499))
  j <- st$instability$joint
  expect_setequal(j$ordering, c('time', 'state eta1'))
  expect_lt(j$p[j$ordering == 'state eta1'], 0.05)
  # A drift that depends on the state also changes how widely the process
  # moves at the extremes, so the diffusion's score answers too; the drift is
  # among the two it points at.
  par <- st$instability$parameters
  expect_true('drift_eta1' %in% utils::head(par$parameter[par$ordering == 'state eta1'], 2))
  expect_output(print(st), 'Stability')
  # Not run unless asked for.
  st0 <- suppressMessages(ctScoreTest(x$fit, tests = 'heterogeneity', cores = 1,
    verbose = FALSE, nboot = 99))
  expect_null(st0$instability)
})

test_that("an innovation sd that varies over people is still tested for volatility", {
  # The volatility terms move each innovation sd into PARS. While random effects
  # took matsetup's order, which puts PARS first, that moved a varying sd ahead
  # of the drift's effect, the extended model's correlations were not the
  # fit's, and state dependence, volatility and fixed cells were all refused:
  # on the Fisher search, from round 6 to the end.
  # The structure is what is tested, not power: twenty people.
  d <- .score_data('het', nsub = 20L)
  m <- .score_model()
  m$pars$indvarying[m$pars$param %in% c('drift_eta1', 'diff_eta1')] <- TRUE
  # These data have no spread in the innovations, so its sd sits at the floor
  # and the fit says so; the model's shape is all that matters here.
  fit <- fit_cached('scoretest het random sd', suppressWarnings(suppressMessages(ctFit(
    d[, c('id', 'time', 'Y1')], m, backend = 'julia', intoverpop = 'laplace',
    cores = 1, verbose = 0))))
  st <- suppressMessages(ctScoreTest(fit, tests = c('statedependence', 'volatility'),
    cores = 1, verbose = FALSE))
  expect_null(st$info$errors)
  expect_false(is.null(st$statedependence))
  expect_false(is.null(st$volatility))
})
