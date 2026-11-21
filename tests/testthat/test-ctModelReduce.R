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
  # The trend has no noise of its own, so it can be a linear trend.
  expect_identical(labels[[2]], c('trend_drift fixed at 0 (a linear trend)',
    "intercepts scaled by the drift (driftScaled = 'CINT')"))
  expect_identical(labels[[3]], "diffusions scaled by the drift (driftScaled = 'DIFFUSION')")
  # The rank step handles the covariance; nothing is known for the rest.
  expect_length(labels[[4]], 0L)
  expect_length(labels[[5]], 0L)
  expect_match(f[[2]]$text, 'trend_drift and cint_trend trade off: trend settles faster')

  # A zero-curvature direction on a diagonal is named, not changed.
  diag <- ctsem:::.ctReduceFindings(.reduce_fit(list(.reduce_dir('mvar2', 1e-15))))
  expect_identical(diag[[1]]$kind, 'unread')
  expect_length(diag[[1]]$candidates, 0L)

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
  scaled <- f[[2]]$candidates[[2]]$edit(m, NULL, NULL)$model
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
      # The linear trend fits better with a parameter fewer; the restated
      # model ends a little below the start.
      .reduce_fit(list(), loglik = if (linear) -95 else -100.5,
        npar = if (linear) 9L else 10L, model = model)
    },
    .ctImproveRankStep = function(base, current, kept, ...) list(current = current,
      kept = kept, data = NULL, args = NULL, history = list(), changed = FALSE))
  out <- ctsem:::.ctImproveReduce(.reduce_model(), start, list(), NULL, NULL, 'LR',
    0.05, 1L)
  expect_true(out$changed)
  expect_equal(out$current$estimate$loglik, -95)
  h <- do.call(rbind, out$history)
  expect_identical(h$kept, c(TRUE, FALSE))
  expect_match(h$note[1], 'better with fewer parameters')
  expect_match(h$note[2], 'restated')
  expect_identical(vapply(out$kept, `[[`, '', 'key'), 'reduce:lineartrend:trend_drift')
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
})
