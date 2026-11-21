# What ctModelImprove(reduceOnly = TRUE) recognises in a fit's identifiability
# report, the changes it proposes, and how it judges them. Fit-free: a model,
# a hand-built report, and refits mocked where the stage decides.

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

.reduce_dir <- function(p, relative = 1e-6) list(parameters = p, relative = relative)

# The criterion, as the verdict asks it (the reduced fit first): the fuller
# model is better where it fits more than 3 log likelihood units higher.
.reduce_compare <- function(current, refit, criterion, alpha, cores = 2L) {
  d <- refit$estimate$loglik - current$estimate$loglik
  list(accept = d > 3, dloglik = d, dpar = 1L, value = d, score = d, note = 'by mock')
}

test_that("flat directions are sorted into the kinds the reduction has changes for", {
  fit <- .reduce_fit(list(
    .reduce_dir('mvar21', 1e-15),
    .reduce_dir(c('trend_drift', 'cint_trend', 'popsd_cint_trend')),
    .reduce_dir(c('drift_x', 'diff_x')),
    .reduce_dir(c('rawcor_a__b', 'rawcor_a__c', 'popsd_a')),
    .reduce_dir('mvar2', 1e-6)))
  f <- ctsem:::.ctReduceFindings(fit)
  kinds <- vapply(f, `[[`, '', 'kind')
  expect_identical(kinds, c('unread', 'drift-intercept', 'drift-diffusion',
    'covariance', 'other'))
  labels <- lapply(f, function(x) vapply(x$candidates, `[[`, '', 'label'))
  expect_identical(labels[[1]], 'mvar21 fixed at 0')
  # The trend has no noise of its own, so it can be a linear trend or a
  # constant level. The intercept's sd goes with the finding.
  expect_identical(labels[[2]], c('trend_drift fixed at 0 (a linear trend)',
    'trend a constant level (trend_drift and cint_trend fixed at 0)',
    "intercepts scaled by the drift (driftScaled = 'CINT')"))
  expect_identical(labels[[3]], "diffusions scaled by the drift (driftScaled = 'DIFFUSION')")
  # The rank step handles the covariance; nothing is known for the rest.
  expect_length(labels[[4]], 0L)
  expect_length(labels[[5]], 0L)
  expect_match(f[[2]]$text, 'trend_drift and cint_trend trade off: trend settles faster')

  # A zero-curvature direction on a variance is named, not changed; on an
  # auto-effect (at its bound near zero) it is fixed there.
  diag <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir('mvar2', 1e-15),
    .reduce_dir('trend_drift', -1e-15), .reduce_dir('drift_x', 1e-15))))
  expect_identical(vapply(diag, `[[`, '', 'kind'), rep('unread', 3))
  expect_length(diag[[1]]$candidates, 0L)
  expect_identical(diag[[2]]$candidates[[1]]$label, 'trend_drift fixed at 0 (a linear trend)')
  expect_identical(diag[[3]]$candidates[[1]]$label, 'drift_x fixed at 0 (a random walk)')
  # The same change under the trade-off and alone has one key.
  expect_identical(diag[[2]]$candidates[[1]]$key, f[[2]]$candidates[[1]]$key)

  # A trade-off carrying its intercept's spread is the trend's, not a
  # covariance's (the dis model at rank 3); one beside a covariance's surplus
  # dimensions is both.
  spread <- c('trend_drift', 'poploading_cint_trend_dim1', 'cint_trend',
    'popunique_cint_trend', 'poploading_cint_trend_dim2', 'poploading_cint_trend_dim3')
  expect_identical(vapply(ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(spread)))),
    `[[`, '', 'kind'), 'drift-intercept')
  mixed <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(c('rawcor_a__b',
    'rawcor_a__c', 'rawcor_b__c', 'trend_drift', 'popsd_cint_trend')))))
  expect_identical(vapply(mixed, `[[`, '', 'kind'), c('drift-intercept', 'covariance'))
  expect_identical(mixed[[2]]$parameters, c('rawcor_a__b', 'rawcor_a__c', 'rawcor_b__c'))
  # The intercept's spread in a direction of its own is the trade-off's too,
  # and so is the spread of the latent's start.
  apart <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(c('trend_drift',
    'cint_trend')), .reduce_dir(c('poploading_cint_trend_dim2', 'popunique_cint_trend')),
    .reduce_dir(c('popunique_T0m_trend', 'trend_drift', 'popunique_cint_trend')))))
  expect_identical(vapply(apart, `[[`, '', 'kind'), 'drift-intercept')

  # The augmented route's own ridge is not the model's.
  partial <- .reduce_dir(c('popsd_a', 'rawcor_a__b'))
  partial$partial <- list(parameters = 'a')
  expect_length(ctsem:::.ctReduceFindings(.reduce_fit(list(partial))), 0L)
})

test_that("the changes edit the model as they say", {
  m <- .reduce_model()
  f <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir('mvar21', 1e-15),
    .reduce_dir(c('trend_drift', 'cint_trend')), .reduce_dir(c('drift_x', 'diff_x')))))
  edit <- function(i, j) f[[i]]$candidates[[j]]$edit(m, NULL, NULL)$model$pars
  cell <- function(p, mat, r, c) p[p$matrix == mat & p$row == r & p$col == c, ]
  fixed <- cell(edit(1, 1), 'MANIFESTVAR', 2, 1)
  expect_true(is.na(fixed$param))
  expect_equal(as.numeric(fixed$value), 0)
  trend <- cell(edit(2, 1), 'DRIFT', 1, 1)
  expect_true(is.na(trend$param))
  expect_equal(as.numeric(trend$value), 0)
  scaled <- f[[2]]$candidates[[3]]$edit(m, NULL, NULL)$model
  expect_identical(scaled$driftScaled, 'CINT')
  expect_match(cell(scaled$pars, 'CINT', 1, 1)$param, 'DRIFT\\[1,1\\] / expm1')
  expect_match(cell(f[[3]]$candidates[[1]]$edit(m, NULL, NULL)$model$pars, 'DIFFUSION',
    2, 2)$param, 'DRIFT\\[2,2\\]')
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

test_that("the stage keeps the best change the criterion allows and stops", {
  start <- .reduce_fit(list(.reduce_dir(c('trend_drift', 'cint_trend'))))
  testthat::local_mocked_bindings(
    .ctImproveRefit = function(current, model, data, cores, start = NULL, args = NULL) {
      drift <- model$pars[model$pars$matrix == 'DRIFT' & model$pars$row == 1 &
        model$pars$col == 1, ]
      linear <- is.na(drift$param)
      level <- linear && is.na(model$pars$param[model$pars$matrix == 'CINT' &
        model$pars$row == 1])
      # The linear trend fits better with a parameter fewer; the constant level
      # ends far below; the restated model ends a little below the start.
      .reduce_fit(list(), loglik = if (level) -110 else if (linear) -95 else -100.5,
        npar = if (level) 8L else if (linear) 9L else 10L, model = model)
    },
    .ctImproveCompare = .reduce_compare,
    .ctImproveRankStep = function(base, current, kept, ...) list(current = current,
      kept = kept, data = NULL, args = NULL, history = list(), changed = FALSE))
  out <- ctsem:::.ctImproveReduce(.reduce_model(), start, list(), NULL, NULL, 'LR',
    0.05, 1L)
  expect_true(out$changed)
  expect_equal(out$current$estimate$loglik, -95)
  h <- do.call(rbind, out$history)
  # The constant level ends far below and is refused by the criterion.
  expect_identical(h$kept, c(TRUE, FALSE, FALSE))
  expect_match(h$note[1], 'better with fewer parameters')
  expect_match(h$note[3], 'restated')
  expect_identical(vapply(out$kept, `[[`, '', 'key'), 'reduce:drift0:trend_drift')
})

test_that("of the changes kept, the smallest wins, and none is refitted twice", {
  start <- .reduce_fit(list(.reduce_dir(c('trend_drift', 'cint_trend'))))
  calls <- character()
  testthat::local_mocked_bindings(
    .ctImproveRefit = function(current, model, data, cores, start = NULL, args = NULL) {
      first <- function(m) model$pars$param[model$pars$matrix == m &
        model$pars$row == 1 & model$pars$col == 1]
      linear <- is.na(first('DRIFT')); level <- linear && is.na(first('CINT'))
      calls <<- c(calls, if (level) 'level' else if (linear) 'linear' else 'scaled')
      # Both reductions cost little; the restatement is as good as the start.
      # Once reduced, the fit still shows the trade-off, so the stage would
      # propose the same changes again.
      .reduce_fit(list(.reduce_dir(c('trend_drift', 'cint_trend'))),
        loglik = if (level) -100.5 else if (linear) -100.2 else -100,
        npar = if (level) 8L else if (linear) 9L else 10L, model = .reduce_model())
    },
    .ctImproveCompare = .reduce_compare,
    .ctImproveRankStep = function(base, current, kept, ...) list(current = current,
      kept = kept, data = NULL, args = NULL, history = list(), changed = FALSE))
  out <- ctsem:::.ctImproveReduce(.reduce_model(), start, list(), NULL, NULL, 'LR',
    0.05, 1L)
  expect_identical(calls, c('linear', 'level', 'scaled'))
  expect_identical(vapply(out$kept, `[[`, '', 'key'), 'reduce:level:trend_drift')
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

test_that("the fit's warning says what each direction is and how to reduce it", {
  dirs <- list(.reduce_dir('mvar21', 1e-15), .reduce_dir(c('trend_drift', 'cint_trend')),
    .reduce_dir(c('rawcor_a__b', 'popsd_a')))
  fit <- .reduce_fit(dirs)
  f <- ctsem:::.ctReduceFindings(fit)
  expect_warning(ctsem:::.ctBackendIdentifyWarn(fit$identifiability, NULL, NULL,
    findings = f), paste0('mvar21 is not read by the likelihood; trend_drift and ',
    'cint_trend trade off.*more dimensions than the subjects identify.*',
    'ctModelImprove\\(fit, reduceOnly = TRUE\\)'))
  # Without findings, the bare list as before.
  expect_warning(ctsem:::.ctBackendIdentifyWarn(fit$identifiability, NULL, NULL),
    'Not estimable as the model stands')
  expect_identical(ctsem:::.ctReduceWarnText(list()), '')
  # A description that starts with a word starts the sentence.
  cov <- Filter(function(x) identical(x$kind, 'covariance'), f)
  expect_match(ctsem:::.ctReduceWarnText(cov), '^The population covariance has more')
})

test_that("own variances at zero leave the level's own variances, the rank kept", {
  m <- .reduce_model()
  m$pars$indvarying[m$pars$param %in% c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')] <- TRUE
  testthat::local_mocked_bindings(.ctReduceOwnVariances = function(fit, model)
    list(rank = 2L, own = c('cint_trend', 'drift_x', 'diff_x', 'trend_drift')))
  # Two at zero span a flat plane, whose directions mix them.
  f <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir(c('popunique_drift_x',
    'popunique_diff_x'), 1e-22), .reduce_dir(c('popunique_diff_x', 'popunique_drift_x'),
    -1e-24)), model = m))
  expect_identical(vapply(f, `[[`, '', 'kind'), 'own variance at zero')
  expect_match(f[[1]]$text, 'the own variances of drift_x, diff_x are at zero')
  ch <- f[[1]]$candidates[[1]]
  expect_identical(ch$family, 'rank')
  expect_identical(ch$label, 'no own variance for drift_x, diff_x')
  expect_identical(ch$edit(m, NULL, NULL)$args,
    list(poprank = 2L, popunique = c('cint_trend', 'trend_drift')))
})
