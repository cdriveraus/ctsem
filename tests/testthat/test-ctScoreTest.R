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
    # 'level': persistence that differs with where the subject sits, its
    # manifest mean. Drawn only here, so the other scenarios' draws are as
    # they were.
    lev <- if (scenario == 'level') stats::rnorm(1) else NA_real_
    if (scenario == 'level') a <- -0.5 * exp(0.6 * lev)
    if (scenario == 'het') a <- -0.5 * exp(0.6 * stats::rnorm(1))
    if (scenario == 'cov') a <- -0.5 * exp(0.3 * z[i])
    if (scenario == 'sdep') b <- 0.4
    mm <- if (scenario == 'level') lev else
      stats::rnorm(1) + if (scenario == 'cov') 0.8 * z[i] else 0
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

test_that("a robust joint test with as many terms as the units span is not reported", {
  set.seed(4)
  # 12 units: their contributions span 11 dimensions at most.
  wide <- matrix(stats::rnorm(12 * 14), 12, 14)
  joint <- ctsem:::.ctScoreJoint(colSums(wide), ctsem:::.ctScoreSpread(wide), wide, 199)
  expect_true(is.na(joint$p))
  expect_true(is.na(joint$LM))
  narrow <- wide[, 1:3]
  joint <- ctsem:::.ctScoreJoint(colSums(narrow), ctsem:::.ctScoreSpread(narrow), narrow, 199)
  expect_equal(joint$df, 3L)
  expect_true(joint$p > 0 && joint$p <= 1)
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
  # exact Hessian's at the estimate. Not the fit's own: the finish may keep one
  # formed within `hessianreuse` se of the estimate, and where this fit's did,
  # 0.006 raw away, an element differed by 1.4 per cent.
  atest <- fit
  atest$uncertainty$hessian <- NULL
  exact <- suppressMessages(ctsem:::.ctBackendHessian(atest, est))
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
  p <- plot(st, top = 3)
  expect_s3_class(p, 'ggplot')
  expect_silent(ggplot2::ggplot_build(p))
  expect_true(all(startsWith(levels(p$data$key), 'covariates|')))
  expect_error(plot(st, family = 'fixed'), 'Nothing to plot')
})

test_that("a parameter that varies over subjects is found", {
  x <- .score_fit('het')
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'heterogeneity', cores = 1,
    verbose = FALSE))
  h <- st$heterogeneity
  expect_setequal(h$parameter, c('drift_eta1', 'diff_eta1', 'mvarY1'))
  # The drift varies independently of the initial state and the intercept, so
  # the variance alone carries it; p, the stronger of that and the test with
  # the covariances, pays for taking the two (.0115 here).
  expect_gt(h$z[h$parameter == 'drift_eta1'], stats::qnorm(0.99))
  expect_lt(h$p[h$parameter == 'drift_eta1'], 0.05)
  expect_gt(h$sd_natural[h$parameter == 'drift_eta1'], 0)
})

test_that("a varying intercept is tested with its covariance with the initial state", {
  # 'null' data have a random manifest mean (sd 1). Fitted with the level in
  # CINT and no random effects, CINT varying alone is subjects drifting apart
  # from a shared start; with its covariance with the initial mean it takes in
  # a random intercept, which is what the data hold. The fixed MANIFESTMEANS is
  # tested as a random intercept on the indicator, from the model with the
  # cell free.
  x <- .score_fit('null')
  model <- function(...) suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1',
    latentNames = 'eta1', LAMBDA = matrix(1), ...))
  fit <- function(m) suppressMessages(ctFit(x$data[, c('id', 'time', 'Y1')], m,
    backend = 'julia', cores = 1, verbose = 0))
  fc <- fit_cached('scoretest null cint level', fit(model(MANIFESTMEANS = matrix(0),
    CINT = matrix('cint||FALSE'), T0MEANS = matrix('t0m||FALSE'))))
  sc <- suppressMessages(ctScoreTest(fc, tests = 'heterogeneity', cores = 1,
    verbose = FALSE, nboot = 199))
  h <- sc$heterogeneity
  hc <- h[h$parameter == 'cint', ]
  expect_identical(hc$with, 't0m')
  expect_equal(hc$df, 2L)
  expect_lt(hc$p, 0.01)
  hm <- h[h$parameter == 'MANIFESTMEANS[Y1]', ]
  expect_equal(nrow(hm), 1L)
  expect_lt(hm$p, 0.01)
  expect_null(sc$info$errors)
  expect_output(print(sc), 'Covariances with: t0m')
  # ctModelImprove's edits: the intercept varies with the initial mean, and
  # the fixed indicator mean becomes a parameter varying about its value.
  base <- ctsem:::.ctFitBaseModel(fc)
  props <- ctsem:::.ctImproveProposals(sc, base, alpha = 0.05)
  keys <- vapply(props, `[[`, '', 'key')
  ec <- props[[match('het:cint', keys)]]$edit(base, x$data, NULL)
  expect_true(all(as.logical(ec$model$pars$indvarying[ec$model$pars$param %in% c('cint', 't0m')])))
  em <- props[[match('het:MANIFESTMEANS[Y1]', keys)]]$edit(base, x$data, NULL)
  expect_equal(em$model$fixedpopmeans, c(mm_Y1 = 0))
  expect_true(all(as.logical(em$model$pars$indvarying[em$model$pars$param %in% c('mm_Y1', 't0m')])))
  expect_no_error(suppressWarnings(suppressMessages(ctFit(x$data[, c('id', 'time', 'Y1')],
    em$model, backend = 'julia', fit = FALSE, cores = 1, verbose = 0))))
  # The covariance's per-subject contribution, s_j s_r + H_jr, and the
  # derivative its projection needs, against differences of the subjects' log
  # likelihoods.
  est <- fc$estimate$raw
  ptab <- ctsem:::.ctScoreParameterTable(fc)
  j <- which(ptab$name == 'cint')
  r <- which(ptab$name == 't0m')
  pieces <- ctsem:::.ctScorePieces(fc, est, rows = c(j, r), second = c(j, r),
    pairs = cbind(j, r), cores = 1)
  it <- ctsem:::.ctScoreHetItem(pieces, ptab, j, r)
  ll <- function(theta) as.numeric(ctsem:::.ctBackendKalmanRaw(fc, theta,
    subjectmatrices = FALSE, fields = 'subject_loglik')$subject_loglik)
  covariance <- function(theta, h = 1e-3) {
    ej <- replace(numeric(length(theta)), j, h)
    er <- replace(numeric(length(theta)), r, h)
    l <- function(d) ll(theta + d)
    (l(ej) - l(-ej)) / (2 * h) * (l(er) - l(-er)) / (2 * h) +
      (l(ej + er) - l(ej - er) - l(er - ej) + l(-ej - er)) / (4 * h^2)
  }
  expect_equal(it$s * pieces$S[, r] + it$H[, r], covariance(est), tolerance = 1e-4)
  k <- which(ptab$name == 'drift_eta1')
  g <- 1e-3
  Dk <- (sum(covariance(replace(est, k, est[k] + g))) -
    sum(covariance(replace(est, k, est[k] - g)))) / (2 * g)
  expect_equal(sum(it$H[, k] * pieces$S[, r]) +
    sum(it$s * ctsem:::.ctScoreHrow(pieces, r)[, k]) + it$Tmix[k, 1], Dk, tolerance = 1e-3)
  # With the initial mean fixed there is nothing to covary with: the variance
  # alone, as for any parameter of a model with no random effects and no free
  # initial mean.
  fs <- fit_cached('scoretest null cint fixed t0', fit(model(MANIFESTMEANS = matrix(0),
    CINT = matrix('cint||FALSE'), T0MEANS = matrix(0))))
  ss <- suppressMessages(ctScoreTest(fs, tests = 'heterogeneity', cores = 1,
    verbose = FALSE, nboot = 199))
  expect_equal(ss$heterogeneity$df[ss$heterogeneity$parameter == 'cint'], 1L)
})

test_that("a measurement error covariance is offered only between errors with a variance", {
  # With the error variances fixed at zero (the book's FI-CLPM) a covariance
  # between the errors has no direction, and printed as an empty row.
  m <- function(MANIFESTVAR) suppressMessages(ctModel(type = 'ct',
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), MANIFESTVAR = MANIFESTVAR))
  cells <- function(model) ctsem:::.ctScoreFixedModel(model)$terms$cell
  expect_false(any(grepl('^MANIFESTVAR', cells(m(matrix(0, 2, 2))))))
  expect_true('MANIFESTVAR[Y2,Y1]' %in%
    cells(m(matrix(c('mvar1', 0, 0, 'mvar2'), 2, 2, byrow = TRUE))))
  expect_false('MANIFESTVAR[Y2,Y1]' %in%
    cells(m(matrix(c('mvar1', 0, 0, 0), 2, 2, byrow = TRUE))))
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

test_that("a parameter that depends on the subject's level is found by the level terms", {
  # The state is the deviation from the subject's mean, which the model puts
  # in the random intercept, so the state dependence terms cannot see a
  # persistence that depends on the mean itself; the level terms ask.
  x <- .score_fit('level', seed = 4L)
  st <- suppressMessages(ctScoreTest(x$fit, tests = 'level', cores = 1, verbose = FALSE))
  expect_null(st$info$errors)
  lv <- st$level$terms
  expect_setequal(lv$parameter, c('drift_eta1', 'diff_eta1', 'mvarY1'))
  expect_equal(length(unique(lv$level)), 1L)
  expect_lt(lv$p[lv$parameter == 'drift_eta1'], 0.01)
  # Persistence shorter for the high: the drift more negative with the level.
  e <- lv$estimate[lv$parameter == 'drift_eta1']
  expect_true(is.na(e) || e < 0)
  expect_true(all(c(lv$level[1], 'all') %in% st$level$joint$level))
  # Nothing where nothing depends on the level.
  n <- .score_fit('null')
  st0 <- suppressMessages(ctScoreTest(n$fit, tests = 'level', cores = 1, verbose = FALSE))
  expect_gt(min(st0$level$terms$p_adj), 0.05)
  # A model with no random intercept has nothing to test.
  m <- .score_model()
  m$pars$indvarying <- FALSE
  fixedmean <- fit_cached('scoretest null no random mean', suppressMessages(ctFit(
    n$data[, c('id', 'time', 'Y1')], m, backend = 'julia', cores = 1, verbose = 0)))
  st1 <- suppressMessages(ctScoreTest(fixedmean, tests = 'level', cores = 1, verbose = FALSE))
  expect_null(st1$level)
  expect_null(st1$info$errors)
})

test_that("a PARS parameter the dynamics read through an expression has level terms", {
  # The same model with its drift written as a PARS parameter that the drift
  # cell reads in an expression: the term shifts the parameter where it is
  # read, so its test is the plain parameter's.
  x <- .score_fit('level', seed = 4L)
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1), DRIFT = matrix('a_eta * 1'),
    PARS = matrix('a_eta|-(1e-06 + 2 * log1p_exp(-(2 * param)))')))
  f <- fit_cached('scoretest level drift in PARS', suppressMessages(ctFit(
    x$data[, c('id', 'time', 'Y1')], m, backend = 'julia', cores = 1, verbose = 0)))
  st <- suppressMessages(ctScoreTest(f, tests = 'level', cores = 1, verbose = FALSE))
  plain <- suppressMessages(ctScoreTest(x$fit, tests = 'level', cores = 1, verbose = FALSE))
  expect_null(st$info$errors)
  a <- st$level$terms[st$level$terms$parameter == 'a_eta', ]
  d <- plain$level$terms[plain$level$terms$parameter == 'drift_eta1', ]
  expect_equal(nrow(a), 1L)
  expect_equal(a$LM, d$LM, tolerance = 1e-4)
  expect_equal(a$estimate, d$estimate, tolerance = 1e-4)
  # A second process held at a gap from the first reads the shifted
  # parameter too, and the gap is tested.
  two <- suppressMessages(ctModel(type = 'ct', manifestNames = c('y1', 'y2'),
    latentNames = c('eta', 'med'), LAMBDA = matrix(c(1, 0, 1, 1), 2, 2, byrow = TRUE),
    DRIFT = matrix(c('a_eta', 0, 0, 'a_eta - d_med'), 2, 2, byrow = TRUE),
    PARS = matrix(c('a_eta|-log1p_exp(-param)', 'd_med|log1p_exp(param)'), ncol = 1),
    CINT = matrix(c(ctPar('c_eta', indvarying = TRUE), 0)),
    MANIFESTMEANS = matrix(0, 2, 1)))
  lm <- ctsem:::.ctScoreLevelModel(two)
  expect_true(all(c('a_eta', 'd_med') %in% lm$terms$parameter))
  p <- lm$model$pars
  ta <- lm$terms$name[lm$terms$parameter == 'a_eta']
  td <- lm$terms$name[lm$terms$parameter == 'd_med']
  expect_equal(p$param[p$matrix == 'DRIFT' & p$row == 2 & p$col == 2],
    paste0('(a_eta + ', ta, ' * c_eta) - (d_med + ', td, ' * c_eta)'))
  expect_equal(p$param[p$matrix == 'DRIFT' & p$row == 1 & p$col == 1],
    paste0('a_eta + ', ta, ' * c_eta'))
  expect_identical(ctsem:::.ctScoreReads(c('a_eta - d', 'a_eta2', 'ba_eta', 'x.a_eta',
    'a_eta(1)', NA), 'a_eta'), c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
})

test_that("the level in the latent is the same linear model as the level in the indicator", {
  # With linear dynamics and the same initial distribution -- the deviation
  # from 0, the level from its stationary mean, each with its stationary
  # variance -- the two forms are one model, and fit to one likelihood.
  n <- .score_fit('null')
  m <- .score_model()
  m$pars$param[m$pars$matrix == 'T0MEANS'] <- NA
  m$pars$value[m$pars$matrix == 'T0MEANS'] <- 0
  m$pars$param[m$pars$matrix == 'T0VAR'] <- 'stationary'
  m$pars$value[m$pars$matrix == 'T0VAR'] <- NA
  m$pars$indvarying[m$pars$matrix %in% c('T0MEANS', 'T0VAR')] <- FALSE
  lm <- ctModelLatentLevel(m)
  expect_identical(attr(lm, 'moved'), 'eta1')
  fit <- function(model) suppressMessages(ctFit(n$data[, c('id', 'time', 'Y1')], model,
    backend = 'julia', intoverpop = 'laplace', cores = 1, verbose = 0))
  a <- fit_cached('scoretest null deviation form', fit(m))
  b <- fit_cached('scoretest null level form', fit(lm))
  expect_equal(b$estimate$loglik, a$estimate$loglik, tolerance = 1e-3)
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
  # largest diagonal, each variance's own z.
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

test_that("time varying covariates are tested beside the verdict, not in it", {
  # A noisy copy of the outcome at each row: the correct model's errors move
  # with it, so the family rejects, and says nothing against the model.
  x <- .score_fit('null')
  set.seed(1)
  tv <- data.frame(id = x$data$id, time = x$data$time,
    copy = x$data$Y1 + stats::rnorm(nrow(x$data), 0, 0.3))
  ad <- suppressMessages(ctFitAdequacy(x$fit, families = c('residuals', 'latent',
    'tvcovariates'), tvcovariates = tv, cores = 1, nboot = 199))
  f <- ad$families
  own <- f$family != 'tvcovariates'
  expect_equal(sum(!own), 1L)
  expect_lt(f$p[!own], 1e-6)
  expect_true(ad$informative)
  expect_equal(f$p_adj[!own], f$p[!own])
  expect_equal(f$p_adj[own], stats::p.adjust(f$p[own], 'holm'))
  expect_identical(ad$adequate, !any(f$p_adj[own] < ad$alpha))
  expect_output(print(ad), 'not counted above')
  without <- suppressMessages(ctFitAdequacy(x$fit, families = 'residuals', cores = 1,
    nboot = 199))
  expect_true(is.na(without$informative))
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
  # The nodes are the layer's, where both effects are coordinates.
  quad <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_quadrature(
    ctsem:::.ctJuliaObjective(lfit), ctsem:::.ctJuliaNumericVector(lest),
    nodes = 3L, contributions = TRUE))
  expect_equal(p3$logmarginal, as.numeric(quad$subject_loglik), tolerance = 1e-10)
  # Every unit's term is a concave quadratic here, so the nodes are placed from
  # the mode and curvature without evaluating any. The same objective without
  # that flag evaluates every node as the quadrature does.
  expect_true(isTRUE(lfit$model_spec$laplace$quadratic))
  evaluated <- lfit$model_spec
  evaluated$laplace$quadratic <- FALSE
  pe <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_posterior_nodes(
    ctsem:::.ctJuliaObjective(structure(evaluated, class = c("ctJuliaModel", "ctFitModel"))),
    ctsem:::.ctJuliaNumericVector(lest), nodes = 3L))
  expect_equal(p3$points, array(as.numeric(pe$points), dim(p3$points)), tolerance = 1e-8)
  expect_equal(p3$weights, as.matrix(pe$weights), tolerance = 1e-8)
  expect_equal(p3$logmarginal, as.numeric(pe$logmarginal), tolerance = 1e-8)
  # The effects family's statistics in closed form there, against the
  # quadrature's on the evaluated nodes, which takes out the rule's own error
  # for the Gaussian it is placed on: here every posterior is that Gaussian,
  # so the two agree to rounding (with the rule's error left in they differed
  # by up to 3e-3 a unit).
  lfe <- lfit
  lfe$model_spec <- evaluated
  closed <- ctsem:::.ctAdequacyEffectMoments(lfit, lest, 1L)$M
  quadrature <- ctsem:::.ctAdequacyEffectMoments(lfe, lest, 1L)$M
  expect_lt(max(abs(closed - quadrature)), 1e-8)

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
  # On the layer both effects are sampled; on the fit's own objective the
  # filter integrates them, so there is nothing to sample (ESS Inf) and the
  # mixture is exact.
  for (hybrid in c(FALSE, TRUE)) {
    obj <- ctsem:::.ctJuliaObjective(lfit, hybrid = hybrid)
    old <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential(obj,
      ctsem:::.ctJuliaNumericVector(lest)))
    mix <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_sequential_mixture(obj,
      ctsem:::.ctJuliaNumericVector(lest), particles = 512L))
    e <- (mix$y[1, , 1] - old$y[1, , 1]) / sqrt(old$ycov[1, , 1, 1])
    expect_lt(sqrt(mean(e^2)), 0.02)
    expect_lt(sqrt(mean((mix$ycov[1, , 1, 1] / old$ycov[1, , 1, 1] - 1)^2)), 0.04)
    if (hybrid) expect_true(all(is.infinite(mix$ess))) else expect_gt(min(mix$ess), 20)
    q7 <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_quadrature(obj,
      ctsem:::.ctJuliaNumericVector(lest), nodes = 7L, contributions = TRUE))
    expect_lt(abs(mean(tapply(mix$llrow, mix$subject, sum) -
      as.numeric(q7$subject_loglik))), 0.02)
  }
})

test_that("an effect the data leave at its prior is named and not tested", {
  # Six units, two effects; along the second every unit's posterior is the
  # prior, so the rule's nodes carry its own weights there. Uncorrected, the
  # rule's error (0.0065 at four nodes) sat in every unit with a spread of
  # rounding, and z reached 1e6 on a reduced-rank fit whose HC effect the data
  # did not inform. The posteriors are built on the product rule the family
  # asks for, placed at N(m, 0.16) along the first effect and the prior along
  # the second, and tilted along the first unless `tilt` is 0.
  units <- 6L
  m <- (seq_len(units) - 3.5) / 2
  tilt <- 0.3
  nodes_for <- function(nodes) {
    gh <- ctsem:::.ctBackendGaussHermite(nodes)
    grid <- as.matrix(expand.grid(gh$node, gh$node))
    G <- as.numeric(gh$weight %o% gh$weight)
    points <- array(0, c(units, nrow(grid), 2L))
    weights <- matrix(0, units, nrow(grid))
    for (u in seq_len(units)) {
      points[u, , 1] <- m[u] + 0.4 * grid[, 1]
      points[u, , 2] <- grid[, 2]
      w <- G * exp(tilt * grid[, 1]^2 / 2)
      weights[u, ] <- w / sum(w)
    }
    list(points = points, weights = weights, logmarginal = numeric(units))
  }
  testthat::local_mocked_bindings(
    .ctBackendUnitDims = function(fit) rep(2L, units),
    .ctBackendSpec = function(fit) list(laplace = list(quadratic = FALSE)),
    .ctBackendPosteriorNodes = function(fit, est, nodes, cores) nodes_for(nodes))
  e <- ctsem:::.ctAdequacyEffectMoments(list(), 0, 1L)
  expect_identical(e$degenerate, 'u2')
  u2 <- e$meta$variable == 'u2'
  expect_true(all(e$M[, u2] == 0))
  expect_true(all(abs(e$M[, !u2]) > 0))
  # Untilted, the posterior along the first is the rule's Gaussian, whose own
  # error is taken out: each unit's value is the closed form.
  tilt <- 0
  e2 <- ctsem:::.ctAdequacyEffectMoments(list(), 0, 1L)
  cuts <- c(-1.5, -0.75, 0, 0.75, 1.5)
  closed <- sapply(cuts, function(cc) pnorm((cc - m) / sqrt(0.75^2 + 0.16)) -
    pnorm(cc / sqrt(1 + 0.75^2)))
  expect_equal(unname(e2$M[, !u2]), unname(closed), tolerance = 1e-12)
})

test_that("covariances between blocks are tested only where the blocks leave them out", {
  # A dynamics block of two effects (one dimension) and a style block of four,
  # two of them the dynamics' (two dimensions): the style block's loadings
  # already move every covariance between the two, which Laplace's error
  # then scored as LM 66 on the hc-wellbeing fits.
  dirs <- function(rows, r) { m <- matrix(0, 5, r); m[rows, ] <- 1; m }
  plan <- list(list(name = 'dynamics', index = 1:2, dirs = dirs(1:2, 1)),
    list(name = 'style', index = 1:4, dirs = dirs(1:4, 2)))
  V <- cbind(plan[[1]]$dirs, plan[[2]]$dirs)
  vb <- c(1L, 2L, 2L)
  pr <- ctsem:::.ctScoreBlockPairs(plan, V, vb)
  expect_equal(NROW(pr), 0L)
  expect_equal(attr(pr, 'nested'), matrix(1:2, 1))
  # Disjoint blocks: each pair of dimensions is tested.
  plan[[2]] <- list(name = 'style', index = 3:4, dirs = dirs(3:4, 2))
  V <- cbind(plan[[1]]$dirs, plan[[2]]$dirs)
  pr <- ctsem:::.ctScoreBlockPairs(plan, V, vb)
  expect_equal(unclass(pr)[, 1:2], matrix(c(1L, 1L, 2L, 3L), 2))
  expect_null(attr(pr, 'nested'))
})

# A reduced-rank level in blocks, as the hc-wellbeing reanalysis fits them: a
# subject's level in one block, item intercepts in another named by matrix,
# one of them about a stated mean. The intercepts are generated from two
# dimensions, which their block's one does not span, independently of the
# level.
.score_blocks_fit <- function() fit_cached('scoretest blocks', {
  set.seed(21)
  lam <- c(1, 0.8, 1.2, 1)
  d <- do.call(rbind, lapply(seq_len(150L), function(i) {
    times <- cumsum(c(0, stats::runif(9L, 0.5, 1.5)))
    lev <- stats::rnorm(1, 0, 0.6)
    f <- stats::rnorm(2)
    mm <- c(0, 0.5 * f[1], 0.5 * f[1] + 0.5 * f[2], -0.5 * f[2])
    eta <- numeric(10L); eta[1] <- stats::rnorm(1)
    for (t in 2:10) {
      e <- exp(-0.6 * (times[t] - times[t - 1]))
      eta[t] <- e * eta[t - 1] + (1 - e) * lev / 0.6 +
        stats::rnorm(1, 0, sqrt((1 - e^2) / 1.2))
    }
    y <- sapply(1:4, function(j) lam[j] * eta + mm[j] + stats::rnorm(10L, 0, 0.4))
    data.frame(id = i, time = times, Y1 = y[, 1], Y2 = y[, 2], Y3 = y[, 3], Y4 = y[, 4])
  }))
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = paste0('Y', 1:4),
    latentNames = 'eta', LAMBDA = matrix(lam), DRIFT = matrix('dr'),
    CINT = matrix(ctPar('lev', indvarying = TRUE)), DIFFUSION = matrix('q'),
    MANIFESTMEANS = matrix(c(ctPar('mm1', indvarying = TRUE, mean = 0),
      ctPar('mm2', indvarying = TRUE), ctPar('mm3', indvarying = TRUE),
      ctPar('mm4', indvarying = TRUE))), T0MEANS = matrix(0), T0VAR = matrix(1)))
  blocks <- list(dynamics = list(rank = 1, effects = 'lev'),
    items = list(rank = 1, effects = 'MANIFESTMEANS'))
  list(data = d, blocks = blocks, fit = suppressWarnings(suppressMessages(ctFit(d, m,
    backend = 'julia', intoverpop = 'laplace', popblocks = blocks, cores = 1,
    verbose = 0))))
})

test_that("a reduced-rank level's blocks are tested for a further dimension and for covariances", {
  x <- .score_blocks_fit()
  fit <- x$fit
  est <- ctsem:::.ctFitRawEstimate(fit)
  ptab <- ctsem:::.ctScoreParameterTable(fit)
  plan <- ctsem:::.ctScoreBlockPlan(fit, est)
  expect_identical(vapply(plan, `[[`, '', 'name'), c('dynamics', 'items'))
  # The stated mean has no coordinate, so the items block tests three.
  expect_setequal(ptab$name[plan[[2]]$index], c('mm2', 'mm3', 'mm4'))
  # The directional pieces are the engine's: summed over subjects, H_i v is
  # the exact Hessian along v, and the mixed second derivative along a pair
  # is the exact Hessian's change along one direction applied to the other.
  V <- cbind(plan[[1]]$dirs, plan[[2]]$dirs)
  pieces <- ctsem:::.ctScorePieces(fit, est, rows = integer(), second = integer(),
    vdirs = V, vpairs = matrix(1:2, 1), cores = 1)
  atest <- fit
  atest$uncertainty$hessian <- NULL
  H <- suppressMessages(ctsem:::.ctBackendHessian(atest, est))
  expect_equal(colSums(pieces$Hv[, 2, ]), as.numeric(H %*% V[, 2]), tolerance = 1e-5)
  h <- 1e-3
  Hp <- suppressMessages(ctsem:::.ctBackendHessian(atest, est + h * V[, 2]))
  Hm <- suppressMessages(ctsem:::.ctBackendHessian(atest, est - h * V[, 2]))
  expect_equal(pieces$Vmix[, 1], as.numeric((Hp - Hm) %*% V[, 1]) / (2 * h),
    tolerance = 1e-3)

  st <- suppressMessages(ctScoreTest(fit, tests = c('heterogeneity', 'level'), cores = 1,
    nboot = 199, verbose = FALSE))
  expect_null(st$info$errors)
  # The items need a second dimension, which their own variances could not
  # give; the level is independent of them.
  dm <- st$blocks$dimension
  expect_identical(dm$block, 'items')
  expect_equal(dm$effects, 3L)
  expect_lt(dm$p, 0.05)
  cv <- st$blocks$covariance
  expect_identical(cv$blocks, 'dynamics x items')
  expect_equal(cv$df, 1L)
  expect_gt(cv$p, 0.05)
  expect_output(print(st), 'further shared dimension')
  # A block naming a matrix reaches the extended models, which move the
  # varying MANIFESTMEANS into PARS: the level family ran.
  expect_true(all(c('mm2', 'mm3', 'mm4', 'lev') %in% st$level$terms$level))
  named <- ctsem:::.ctScoreNamedLevels(fit, ctsem:::.ctFitReplayArgs(fit))
  expect_identical(named$popblocks$id$items$effects, paste0('mm', 1:4))
  # An extended model takes the fit's known modes with it: without them a
  # unit's search can miss a maximum the fit's reached, and the extension did
  # not reproduce the fit at zero (0.95 nats on a random-drift model).
  withknown <- fit
  withknown$model_spec$laplace$known <- list(units = 1L, values = matrix(0.5, 1, 1))
  fx <- ctsem:::.ctScoreFixedModel(ctsem:::.ctFitBaseModel(fit))
  ext <- ctsem:::.ctScoreExtendedSpec(withknown, est, fx$model, fx$terms$name, 'fixed',
    cores = 1)
  expect_identical(ext$spec$laplace$known, withknown$model_spec$laplace$known)
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
  # The fit's filter integrates its initial mean and item intercept beside the
  # random drift and noise; the extension's terms read the state, so its
  # filter cannot, and both are taken one-step, which the result says.
  expect_setequal(st$info$onestep, c('statedependence', 'volatility'))
  expect_output(print(st), 'On one-step Laplace')
  # By quadrature every subject's effects are integrated on its rule, so no
  # family is one-step; the units' scores on the placed objective sum to its
  # gradient, the prior's shares included.
  sq <- suppressMessages(ctScoreTest(fit, tests = 'volatility', cores = 1,
    verbose = FALSE, nboot = 99L, intoverpop = 'quadrature'))
  expect_null(sq$info$errors)
  expect_identical(sq$info$quadrature, 'volatility')
  expect_null(sq$info$onestep)
  expect_true(all(is.finite(sq$volatility$joint$LM)))
  expect_output(print(sq), 'By quadrature')
  vol <- ctsem:::.ctScoreVolatilityModel(ctsem:::.ctFitBaseModel(fit))
  ext <- ctsem:::.ctScoreExtendedSpec(fit, fit$estimate$raw, vol$model, vol$terms$name,
    'volatility', 1L)
  pl <- ctsem:::.ctScorePlaced(ext, 1L)
  expect_equal(colSums(pl$units), unname(pl$gradient), tolerance = 1e-10)
})

test_that("where effects reach the dynamics, only the unit-level families are projected", {
  # There the errors come from the mixture predictive and the scores from the
  # Laplace objective, and the row-level statistics' projected variance was
  # unreliable while unprojected they held their level (.ctAdequacyProject).
  # The same fit as above, cached.
  d <- .score_data('het', nsub = 20L)
  m <- .score_model()
  m$pars$indvarying[m$pars$param %in% c('drift_eta1', 'diff_eta1')] <- TRUE
  fit <- fit_cached('scoretest het random sd', suppressWarnings(suppressMessages(ctFit(
    d[, c('id', 'time', 'Y1')], m, backend = 'julia', intoverpop = 'laplace',
    cores = 1, verbose = 0))))
  expect_true(ctsem:::.ctAdequacyApproximate(fit, mixture = FALSE))
  check <- function(projection) suppressMessages(ctFitAdequacy(fit,
    families = c('residuals', 'distribution'), projection = projection, cores = 1, nboot = 199))
  exact <- check('derivative')
  plain <- check('none')
  expect_identical(exact$projection, 'derivative (unit-level families)')
  rowlevel <- exact$statistics$family == 'residuals'
  expect_true(any(rowlevel) && any(exact$statistics$family == 'distribution'))
  # The mixture's draws are seeded, so unprojected statistics match exactly.
  expect_equal(exact$statistics$z[rowlevel], plain$statistics$z[rowlevel])
  expect_false(isTRUE(all.equal(exact$statistics$z[!rowlevel], plain$statistics$z[!rowlevel])))
})

test_that("errors from the mixture predictive are its own distribution function's", {
  # The exact predictive integrates the random effects by importance sampling:
  # a mixture of normals, whose errors standardised by its mean and covariance
  # are heavier tailed than normal. Through the mixture's own distribution
  # function they are standard normal under the model.
  mix <- function(fit) {
    spec <- ctsem:::.ctBackendSpec(fit)
    k <- suppressMessages(ctsem:::.ctBackendKalmanRaw(fit, ctsem:::.ctFitRawEstimate(fit),
      subjectmatrices = FALSE, fields = 'subject', sequential = 'mixture', compact = TRUE))
    nrows <- length(spec$times)
    y <- as.numeric(spec$manifest_data)
    y[!is.finite(y) | y == 99999] <- NA
    list(standardised = (y - as.numeric(k$y1)) / sqrt(as.numeric(k$ycov1)),
      transformed = stats::qnorm(as.numeric(k$pitgauss)),
      sequential = stats::qnorm(as.numeric(k$pitseq)), ok = !is.na(y))
  }
  # Where the effects enter linearly (a random intercept, on the Laplace route,
  # where the exact predictive is the mixture) the predictive is normal, and
  # the transforms give the standardised errors back, to the sampling's
  # accuracy; with one item the two transforms are one.
  n <- .score_fit('null')
  lfit <- fit_cached('scoretest null laplace', suppressMessages(ctFit(
    n$data[, c('id', 'time', 'Y1')], .score_model(), backend = 'julia',
    intoverpop = 'laplace', cores = 1, verbose = 0)))
  lin <- mix(lfit)
  expect_true(all(is.finite(lin$transformed[lin$ok])))
  expect_equal(lin$transformed[lin$ok], lin$sequential[lin$ok])
  # The sampling's error is largest at a person's first row, where the
  # posterior is the prior and 512 draws cover its tails coarsely (0.2 per
  # cent of rows differ by more than 0.05 there, none at 4096 draws).
  dev <- abs(lin$transformed - lin$standardised)[lin$ok]
  expect_lt(stats::median(dev), 1e-3)
  expect_lt(stats::quantile(dev, 0.99, names = FALSE), 0.05)
  # Where they reach the dynamics the two differ, and the transforms' tails
  # are the lighter.
  d <- .score_data('het', nsub = 20L)
  m <- .score_model()
  m$pars$indvarying[m$pars$param %in% c('drift_eta1', 'diff_eta1')] <- TRUE
  fit <- fit_cached('scoretest het random sd', suppressWarnings(suppressMessages(ctFit(
    d[, c('id', 'time', 'Y1')], m, backend = 'julia', intoverpop = 'laplace',
    cores = 1, verbose = 0))))
  dyn <- mix(fit)
  expect_true(all(is.finite(dyn$transformed[dyn$ok])))
  expect_gt(max(abs(dyn$transformed - dyn$standardised)[dyn$ok]), 0.01)
  # The adequacy check takes them: its errors are the transforms.
  err <- ctsem:::.ctAdequacyErrors(fit)
  expect_equal(as.numeric(err$E[, 1])[dyn$ok], dyn$transformed[dyn$ok])
})

test_that("a skewed statistic and a family's sum are read against references of their own shape", {
  # A statistic of known skewness is read through a chi-square of that
  # skewness, so a chi-square's quantiles come back as the normal's.
  nu <- c(20, 200)
  for (pr in c(0.999, 0.01)) {
    q <- stats::qchisq(pr, nu)
    expect_equal(ctsem:::.ctAdequacySkewZ((q - nu) / sqrt(2 * nu), sqrt(8 / nu)),
      rep(stats::qnorm(pr), 2), tolerance = 1e-8)
  }
  # A family's sum of squared z, from its terms between units: k statistics
  # not few beside n units (here more), each unit's contributions independent
  # and of mean zero under the model. The p-values are uniform, and a small
  # shift in every statistic's mean, too small for any one z, is found.
  sump <- function(shift) withr::with_seed(1, vapply(seq_len(200), function(r) {
    X <- matrix(stats::rnorm(60 * 150), 60) + shift
    ctsem:::.ctAdequacySumP(X, rep(60, 150))
  }, 0))
  null <- sump(0)
  expect_gt(mean(null < 0.05), 0.02)
  expect_lt(mean(null < 0.05), 0.09)
  expect_gt(stats::ks.test(null, 'punif')$p.value, 0.01)
  # The shift puts each z's mean at 0.6; the largest of 150 such z is often
  # within the null's spread, the sum of their squares is not.
  expect_gt(mean(sump(0.6 / sqrt(60)) < 0.05), 0.8)
})

test_that("a unit with two maxima does not stop the extended families", {
  # Short records with a random drift (nullrate.R's hetA data set 51, in
  # scripts/burstict-2026-11). At this fit one person has two maxima a shallow
  # valley apart where every effect is a coordinate, and one where the
  # filter's intercepts are profiled out, so the one-step term over the
  # hybrid's coordinates is not the layer's, where the extension's terms put
  # every effect: the extension at zero must still reproduce the fit, so both
  # are read on the layer. Before, about one data set in twenty here refused.
  set.seed(1e6 + 51)
  d <- do.call(rbind, lapply(seq_len(100L), function(i) {
    a <- -0.5 * exp(0.6 * stats::rnorm(1))
    mm <- stats::rnorm(1)
    times <- cumsum(c(0, stats::rexp(7L, 1)))
    x <- numeric(8L)
    x[1] <- stats::rnorm(1, 0, 0.8 / sqrt(-2 * a))
    for (t in 2:8) {
      e <- exp(a * (times[t] - times[t - 1]))
      x[t] <- e * x[t - 1] + stats::rnorm(1, 0, sqrt(0.64 * (1 - e^2) / (-2 * a)))
    }
    data.frame(id = i, time = times, Y1 = x + mm + stats::rnorm(8L, 0, 0.5))
  }))
  m <- suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1), silent = TRUE))
  m$pars$indvarying[m$pars$matrix == 'DRIFT'] <- TRUE
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = 'julia', cores = 1,
    verbose = 0, optimcontrol = list(checks = FALSE, laplace_onestep = TRUE))))
  est <- fit$estimate$raw
  apart <- ctsem:::.ctJuliaEvaluate(fit, est, gradient = FALSE, cores = 1,
    hybrid = FALSE)$value - ctJuliaEvaluate(fit, est, gradient = FALSE, cores = 1)$value
  skip_if(abs(apart) < 1e-3, 'this fit stopped where no unit has two maxima')
  st <- suppressMessages(ctScoreTest(fit, tests = 'statedependence', cores = 1,
    nboot = 99L, verbose = FALSE))
  expect_length(st$info$errors, 0)
  j <- st$statedependence$joint
  expect_true(is.finite(j$LM[j$equation == 'all']))
})
