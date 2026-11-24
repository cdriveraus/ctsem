# What ctModelImprove(reduceOnly = TRUE) reads off a fit's identifiability
# report, the changes it proposes, and how it judges them. Fit-free: a model,
# a hand-built report, and the engine calls and refits mocked.

.reduce_model <- function() suppressMessages(ctModel(type = 'ct',
  manifestNames = c('y1', 'y2'), latentNames = c('trend', 'x'), LAMBDA = diag(2),
  DRIFT = matrix(c('trend_drift', 0, 0, 'drift_x'), 2),
  CINT = matrix(c('cint_trend', 0), 2),
  DIFFUSION = matrix(c(0, 0, 0, 'diff_x'), 2),
  MANIFESTVAR = matrix(c('mvar1', 'mvar21', 0, 'mvar2'), 2),
  MANIFESTMEANS = matrix(0, 2)))

.reduce_fit <- function(directions, effects = NULL, loglik = -100, npar = 10L,
    model = .reduce_model()) {
  structure(list(modelbase = model,
    estimate = list(loglik = loglik, raw = numeric(npar)),
    identifiability = list(nweak = length(directions), directions = directions,
      parameters = unique(unlist(lapply(directions, `[[`, 'parameters'))),
      effects = effects)), class = 'ctJuliaFit')
}

# A flat direction: its coordinates, their loadings (normalised), and its
# curvature relative to the largest.
.reduce_dir <- function(p, loadings = rep(1, length(p)), relative = 1e-6)
  list(parameters = p, loadings = loadings / sqrt(sum(loadings^2)), relative = relative)

# The ranges `.ctReduceValues()` would read off the transforms: an auto-effect
# ends at 0, an intercept or loading spans 0, a variance ends at 0.
.reduce_values <- function(fit, params) {
  out <- list(trend_drift = list(ends = 0, zero = FALSE),
    drift_x = list(ends = 0, zero = FALSE),
    cint_trend = list(ends = numeric(), zero = TRUE),
    mvar21 = list(ends = numeric(), zero = TRUE),
    mvar2 = list(ends = 0, zero = FALSE))
  out[intersect(params, names(out))]
}

# The criterion, as the verdict asks it (the reduced fit first): the fuller
# model is better where it fits more than 3 log likelihood units higher.
.reduce_compare <- function(current, refit, criterion, alpha, cores = 2L) {
  d <- refit$estimate$loglik - current$estimate$loglik
  list(accept = d > 3, dloglik = d, dpar = 1L, value = d, score = d, note = 'by mock')
}

test_that("a coordinate is undetermined by its share of the flat subspace", {
  # The report lists a direction twice where two eigenvalues are equal; the
  # shares are taken on an orthonormal basis, so each counts once.
  flat <- ctsem:::.ctReduceFlat(list(directions = list(
    .reduce_dir(c('a', 'b'), c(0.8, 0.6)),
    .reduce_dir(c('c', 'd'), c(1, 1), relative = 1e-22),
    .reduce_dir(c('d', 'c'), c(1, 1), relative = -1e-24),
    .reduce_dir('e', relative = 1e-15))))
  expect_equal(unname(flat$share[c('a', 'b', 'c', 'd', 'e')]), c(0.64, 0.36, 0.5, 0.5, 1))
  expect_identical(flat$unread, 'e')
  # An augmented route's own ridge is not the model's.
  partial <- .reduce_dir(c('popsd_a', 'rawcor_a__b'))
  partial$partial <- list(parameters = 'a')
  expect_length(ctsem:::.ctReduceFlat(list(directions = list(partial)))$share, 0L)
})

test_that("findings follow the shares, whatever the parameters are", {
  testthat::local_mocked_bindings(.ctReduceValues = .reduce_values)
  fit <- .reduce_fit(list(
    .reduce_dir('mvar21', relative = 1e-15),
    .reduce_dir(c('trend_drift', 'cint_trend'), c(0.8, 0.6)),
    .reduce_dir(c('rawcor_a__b', 'rawcor_a__c', 'rawcor_b__c', 'popsd_a'))))
  f <- ctsem:::.ctReduceFindings(fit)
  expect_identical(vapply(f, `[[`, '', 'kind'), c('parameters', 'covariance'))
  # cint_trend has 0.36 of its variation in the flat subspace: determined
  # mostly, so not offered. The covariance has a dimension in its coordinates
  # with no one coordinate undetermined.
  expect_identical(f[[1]]$parameters, c('mvar21', 'trend_drift'))
  expect_identical(f[[1]]$text, paste('mvar21 is not read by the likelihood;',
    'trend_drift is not determined by the data'))
  labels <- vapply(f[[1]]$candidates, `[[`, '', 'label')
  # Simple values only: fixed at the estimate would state no simpler model.
  expect_identical(labels, c('mvar21 fixed at 0', 'trend_drift fixed at 0'))
  expect_length(f[[2]]$candidates, 0L)
  expect_true(f[[2]]$reducible)
  # A warning reads the kinds without evaluating any range.
  testthat::local_mocked_bindings(.ctReduceValues = function(...) stop('evaluated'))
  expect_length(ctsem:::.ctReduceFindings(fit, candidates = FALSE)[[1]]$candidates, 0L)
})

test_that("a parameter whose transform saturated is at the end of its range", {
  # Zero curvature alone, as an unread parameter has, but the optimiser
  # reported its transform saturated: a variance estimated at 0.
  testthat::local_mocked_bindings(.ctReduceValues = .reduce_values)
  dir <- .reduce_dir('mvar2', relative = -1e-18)
  dir$saturated <- 'mvar2'
  fit <- .reduce_fit(list(dir, .reduce_dir('mvar21', relative = 1e-15)))
  flat <- ctsem:::.ctReduceFlat(fit$identifiability)
  expect_identical(flat$saturated, 'mvar2')
  expect_identical(flat$unread, 'mvar21')
  f <- ctsem:::.ctReduceFindings(fit)
  expect_identical(f[[1]]$text, paste('mvar21 is not read by the likelihood;',
    'mvar2 is at the end of its range'))
  expect_true('mvar2 fixed at 0' %in% vapply(f[[1]]$candidates, `[[`, '', 'label'))
})

test_that("the changes edit the model as they say", {
  testthat::local_mocked_bindings(.ctReduceValues = .reduce_values)
  m <- .reduce_model()
  m$pars$indvarying[m$pars$param %in% 'cint_trend'] <- TRUE
  f <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(c('trend_drift',
    'cint_trend'), c(0.6, 0.8))), model = m))
  cand <- f[[1]]$candidates
  expect_identical(vapply(cand, `[[`, '', 'label'), "cint_trend's mean fixed at 0")
  # A varying parameter's mean is stated, its random effect kept.
  edited <- cand[[1]]$edit(m, NULL, NULL)$model
  expect_equal(edited$fixedpopmeans, c(cint_trend = 0))
  expect_true(any(as.logical(edited$pars$indvarying[edited$pars$param %in% 'cint_trend'])))
  # At an end of its range no mean can be stated: fixed, and not varying.
  m$pars$indvarying[m$pars$param %in% 'trend_drift'] <- TRUE
  g <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(c('trend_drift',
    'cint_trend'), c(0.8, 0.6))), model = m))
  expect_identical(g[[1]]$candidates[[1]]$label, 'trend_drift fixed at 0, not varying')
  pars <- g[[1]]$candidates[[1]]$edit(m, NULL, NULL)$model$pars
  cell <- pars[pars$matrix == 'DRIFT' & pars$row == 1 & pars$col == 1, ]
  expect_true(is.na(cell$param))
  expect_equal(as.numeric(cell$value), 0)
})

test_that("a parameter's simple values are read off its transform", {
  drift <- function(x) -(1e-6 + 2 * log1p(exp(-2 * x)))
  fit <- .reduce_fit(list(), npar = 3L)
  fit$estimate$raw <- c(0.5, 1, -2)
  testthat::local_mocked_bindings(
    .ctFitNameRaw = function(fit, raw) stats::setNames(as.numeric(raw), c('d', 'c', 'v')),
    .ctBackendPopMeanSamples = function(fit, samples) list(parnumber = 1:3,
      values = cbind(drift(samples[, 1]), samples[, 2], exp(samples[, 3]))))
  v <- ctsem:::.ctReduceValues(fit, c('d', 'c', 'v', 'missing'))
  expect_identical(names(v), c('d', 'c', 'v'))
  # An auto-effect ends at 0 above and nowhere below.
  expect_identical(v$d$ends, 0)
  expect_false(v$d$zero)
  # An intercept has no end and spans zero.
  expect_length(v$c$ends, 0L)
  expect_true(v$c$zero)
  # A variance ends at 0 below.
  expect_identical(v$v$ends, 0)
})

test_that("undetermined sds and own variances are removed, together first", {
  m <- .reduce_model()
  m$pars$indvarying[m$pars$param %in% c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')] <- TRUE
  testthat::local_mocked_bindings(.ctReduceOwnVariances = function(fit, model)
    list(rank = 2L, own = c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')))
  # Two at their floor span a flat plane, whose directions mix them.
  f <- ctsem:::.ctReduceFindings(.reduce_fit(list(
    .reduce_dir(c('popunique_drift_x', 'popunique_diff_x'), c(0.6, 0.8), 1e-22),
    .reduce_dir(c('popunique_diff_x', 'popunique_drift_x'), c(0.6, -0.8), -1e-24),
    .reduce_dir('popsd_cint_trend')), model = m))
  expect_identical(vapply(f, `[[`, '', 'kind'), c('own variances', 'sds'))
  expect_match(f[[1]]$text, 'the data do not determine the own variances of')
  cand <- f[[1]]$candidates
  expect_length(cand, 2L)
  expect_identical(cand[[1]]$family, 'rank')
  expect_false(isTRUE(cand[[1]]$fallback))
  expect_true(cand[[2]]$fallback)
  both <- cand[[1]]$edit(m, NULL, NULL)$args
  expect_identical(both$poprank, 2L)
  expect_setequal(both$popunique, c('cint_trend', 'trend_drift'))
  expect_identical(f[[2]]$candidates[[1]]$label, 'cint_trend not varying over subjects')
  # Their coordinates are not a covariance with a dimension too many.
  expect_false(ctsem:::.ctImproveWeakCovariance(.reduce_fit(list(
    .reduce_dir(c('popunique_drift_x', 'popunique_diff_x'))), model = m)))
})

test_that("random effects the data hold at zero or barely determine are made fixed", {
  m <- .reduce_model()
  m$pars$indvarying[m$pars$param %in% c('cint_trend', 'drift_x')] <- TRUE
  effects <- list(table = data.frame(effect = c('cint_trend', 'drift_x'),
    level = 'id', weak = c(TRUE, FALSE), referenceinformation = c(NA, 50)))
  f <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir('popsd_drift_x')),
    effects = effects, model = m))
  expect_identical(vapply(f, `[[`, '', 'kind'), c('no spread', 'weak'))
  pars <- f[[2]]$candidates[[1]]$edit(m, NULL, NULL)$model$pars
  expect_false(any(as.logical(pars$indvarying[pars$param %in% 'cint_trend'])))
  # Their coordinates are not a covariance with a dimension too many.
  expect_false(ctsem:::.ctImproveWeakCovariance(.reduce_fit(
    list(.reduce_dir('popsd_drift_x')), effects = effects, model = m)))
})

test_that("a parameter whose simple values are refused stays free", {
  calls <- character()
  first <- function(model, m) model$pars$param[model$pars$matrix == m &
    model$pars$row == 1 & model$pars$col == 1]
  testthat::local_mocked_bindings(
    .ctReduceValues = .reduce_values,
    .ctImproveRefit = function(current, model, data, cores, start = NULL, args = NULL) {
      if (!is.null(model$fixedpopmeans)) stop('no mean stated here')
      fixed <- is.na(first(model, 'DRIFT'))
      calls <<- c(calls, if (fixed) 'zero' else 'none')
      # Fixed at 0 the trend ends far below and is refused. The current fit's
      # report still shows the direction, which the stage must not refit again.
      .reduce_fit(list(), loglik = -110, npar = 9L)
    },
    .ctImproveCompare = .reduce_compare,
    .ctImproveRankStep = function(base, current, kept, ...) list(current = current,
      kept = kept, data = NULL, args = NULL, history = list(), changed = FALSE))
  start <- .reduce_fit(list(.reduce_dir('trend_drift')))
  out <- ctsem:::.ctImproveReduce(.reduce_model(), start, list(), NULL, NULL, 'LR',
    0.05, 1L)
  expect_identical(calls, 'zero')
  expect_false(out$changed)
  expect_length(out$kept, 0L)
  h <- do.call(rbind, out$history)
  expect_identical(h$kept, FALSE)

  # Where a simple value is kept, the change is the model's.
  calls <- character()
  testthat::local_mocked_bindings(.ctImproveCompare = function(current, refit, ...)
    list(accept = FALSE, dloglik = 0, dpar = 1L, value = 0, score = 0, note = 'by mock'))
  out <- ctsem:::.ctImproveReduce(.reduce_model(), start, list(), NULL, NULL, 'LR',
    0.05, 1L)
  expect_identical(calls, 'zero')
  expect_identical(vapply(out$kept, `[[`, '', 'key'), 'reduce:fix:trend_drift=0')
})

test_that("a full covariance goes to shared dimensions with own variances", {
  # No more shared dimensions than keep it below the full covariance's count.
  expect_identical(ctsem:::.ctReduceFactorMax(10L), 5L)
  expect_identical(ctsem:::.ctReduceFactorMax(4L), 1L)
  # Three effects: one shared dimension with own variances is as many
  # parameters as the full covariance, so a lower rank goes without them.
  expect_identical(ctsem:::.ctReduceFactorMax(3L), 0L)
  set.seed(2)
  L <- matrix(stats::rnorm(16), 8)
  S <- tcrossprod(L) + diag(stats::runif(8, .2, .5))
  expect_identical(ctsem:::.ctReduceFactorRank(S), 2L)
  fa <- ctsem:::.ctReduceFactor(S, 2L)
  expect_equal(fa$shared + diag(fa$own^2), S, tolerance = 0.05)
  ch <- ctsem:::.ctImproveRankChange(3L, TRUE)
  expect_match(ch$label, 'rank 3 with own variances')
  m <- .reduce_model()
  m$pars$indvarying[m$pars$param %in% c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')] <- TRUE
  expect_identical(ch$edit(m, NULL, NULL)$args, list(poprank = 3L, popunique = TRUE))
})

test_that("a rank change drops own variances of parameters no longer varying", {
  m <- .reduce_model()
  m$pars$indvarying <- m$pars$param %in% c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')
  ch <- ctsem:::.ctImproveRankChange(2L, c('cint_trend', 'drift_x', 'diff_x'))
  m$pars$indvarying[m$pars$param %in% 'drift_x'] <- FALSE
  expect_identical(ch$edit(m, NULL, NULL)$args,
    list(poprank = 2L, popunique = c('cint_trend', 'diff_x')))
  m$pars$indvarying[m$pars$param %in% c('cint_trend', 'diff_x')] <- FALSE
  # One varying parameter left: the rank is no longer below it.
  expect_identical(ch$edit(m, NULL, NULL)$args, list(poprank = NA))
})

test_that("the fit's warning says what is flat and how to reduce it", {
  dirs <- list(.reduce_dir('mvar21', relative = 1e-15), .reduce_dir('trend_drift'),
    .reduce_dir(c('rawcor_a__b', 'rawcor_a__c', 'rawcor_b__c', 'popsd_a')))
  fit <- .reduce_fit(dirs)
  f <- ctsem:::.ctReduceFindings(fit, candidates = FALSE)
  expect_warning(ctsem:::.ctBackendIdentifyWarn(fit$identifiability, NULL, NULL,
    findings = f), paste0('mvar21 is not read by the likelihood; trend_drift is not ',
    'determined by the data; the population covariance has more dimensions than ',
    'the subjects identify.*ctModelImprove\\(fit, reduceOnly = TRUE\\)'))
  # Without findings, the bare list as before.
  expect_warning(ctsem:::.ctBackendIdentifyWarn(fit$identifiability, NULL, NULL),
    'Not estimable as the model stands')
  expect_identical(ctsem:::.ctReduceWarnText(list()), '')
  # A description that starts with a word starts the sentence.
  cov <- Filter(function(x) identical(x$kind, 'covariance'), f)
  expect_match(ctsem:::.ctReduceWarnText(cov), '^The population covariance has more')
})
