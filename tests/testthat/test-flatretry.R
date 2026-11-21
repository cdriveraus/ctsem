# Hold-then-free for a model parameter a fit left run off where the likelihood
# is flat (R/ctFitFlatRetry.R). Fit-free: a model, a hand-built report, and the
# refits mocked.

.retry_model <- function() suppressMessages(ctModel(type = 'ct',
  manifestNames = c('y1', 'y2'), latentNames = c('trend', 'x'), LAMBDA = diag(2),
  DRIFT = matrix(c('trend_drift', 0, 0, 'drift_x'), 2),
  CINT = matrix(c('cint_trend', 0), 2),
  DIFFUSION = matrix(c(0, 0, 0, 'diff_x'), 2),
  MANIFESTVAR = matrix(c('mvar1', 0, 0, 'mvar2'), 2),
  MANIFESTMEANS = matrix(0, 2)))

.retry_raw <- c(trend_drift = -16, drift_x = 0.4, cint_trend = 9, diff_x = 0.2,
  mvar1 = -0.5, mvar2 = 0.1)

.retry_fit <- function(raw = .retry_raw, directions = list(), saturated = character(),
    loglik = -100) {
  structure(list(modelbase = .retry_model(),
    estimate = list(loglik = loglik, raw = unname(raw)),
    optim = list(saturated_parameters = saturated),
    identifiability = list(nweak = length(directions), directions = directions),
    args = list(input = list(tag = 'the call'), resolved = list())),
    class = 'ctJuliaFit', rawnames = names(raw))
}

.retry_dir <- function(p) list(parameters = p, loadings = rep(1, length(p)) /
  sqrt(length(p)), relative = 1e-9)

.retry_names <- function(fit, raw) stats::setNames(as.numeric(raw), attr(fit, 'rawnames'))

test_that("a parameter qualifies where it is flat and has run past the start's bound", {
  testthat::local_mocked_bindings(.ctFitNameRaw = .retry_names)
  # trend_drift undetermined and far out; drift_x undetermined but near the
  # start, not run off; cint_trend far out but determined.
  fit <- .retry_fit(directions = list(.retry_dir('trend_drift'), .retry_dir('drift_x')))
  expect_identical(ctsem:::.ctFitFlatRetryParameters(fit), 'trend_drift')
  # A saturated transform qualifies too.
  fit <- .retry_fit(saturated = 'cint_trend')
  expect_identical(ctsem:::.ctFitFlatRetryParameters(fit), 'cint_trend')
  expect_length(ctsem:::.ctFitFlatRetryParameters(.retry_fit()), 0L)
})

test_that("the parameter is held at the other side, then freed, and the better fit kept", {
  calls <- list()
  refit <- function(current, model, data, cores, start = NULL, optimcontrol = list()) {
    pars <- model$pars
    cell <- pars[pars$matrix == 'DRIFT' & pars$row == 1 & pars$col == 1, ]
    held <- is.na(cell$param)
    calls[[length(calls) + 1L]] <<- list(held = held, value = as.numeric(cell$value),
      start = start, optimcontrol = optimcontrol)
    # The held fit has no trend_drift coordinate; the freed one lands higher,
    # with the drift at its slow end and nothing left to retry.
    raw <- if (held) .retry_raw[names(.retry_raw) != 'trend_drift'] else
      replace(.retry_raw, 'trend_drift', 9)
    f <- .retry_fit(raw = raw, loglik = if (held) -95 else -80)
    list(fit = f, warnings = list(simpleWarning(if (held) 'held' else 'freed')))
  }
  testthat::local_mocked_bindings(
    .ctFitNameRaw = .retry_names,
    .ctFitRawEstimate = function(fit) fit$estimate$raw,
    .ctFitRawParNames = function(fit) attr(fit, 'rawnames'),
    .ctFitLongData = function(fit) data.frame(id = 1, time = 0, y1 = 0, y2 = 0),
    .ctFitOptimisedLoglik = function(fit) fit$estimate$loglik,
    .ctFitReplayArgs = function(fit) list(optimcontrol = list()),
    .ctFitFlatRetryValue = function(fit, param, rawvalue) -0.005 * sign(rawvalue),
    .ctFitFlatRetryRefit = refit)
  fit <- .retry_fit(directions = list(.retry_dir('trend_drift')))
  said <- testthat::capture_messages(out <- ctsem:::.ctFitFlatRetry(fit,
    list(simpleWarning('original')), cores = 1L))
  expect_true(any(grepl('held at', said)))
  expect_true(any(grepl('Kept the refit', said)))
  # Held at the other side of its range, without a Hessian.
  expect_true(calls[[1]]$held)
  expect_equal(calls[[1]]$value, -0.005)
  expect_true(isTRUE(calls[[1]]$optimcontrol$estonly))
  # Freed from there: the parameter's raw start on the held side, a full fit.
  expect_false(calls[[2]]$held)
  expect_equal(calls[[2]]$start[['trend_drift']], 3)
  expect_false(isTRUE(calls[[2]]$optimcontrol$estonly))
  # The freed fit is kept, with its own warnings and the caller's call.
  expect_equal(out$fit$estimate$loglik, -80)
  expect_identical(vapply(out$warnings, conditionMessage, ''), 'freed')
  expect_identical(out$fit$args$input$tag, 'the call')
  expect_identical(out$fit$flatretry$parameter, 'trend_drift')
  expect_true(out$fit$flatretry$kept)
  # Tried once: the freed fit leaves the drift saturated at its other end,
  # and it is not sent back.
  expect_length(calls, 2L)
})

test_that("a refit that ends lower leaves the fit and its warnings", {
  testthat::local_mocked_bindings(
    .ctFitNameRaw = .retry_names,
    .ctFitRawEstimate = function(fit) fit$estimate$raw,
    .ctFitRawParNames = function(fit) attr(fit, 'rawnames'),
    .ctFitLongData = function(fit) data.frame(id = 1, time = 0, y1 = 0, y2 = 0),
    .ctFitOptimisedLoglik = function(fit) fit$estimate$loglik,
    .ctFitReplayArgs = function(fit) list(optimcontrol = list()),
    .ctFitFlatRetryValue = function(fit, param, rawvalue) -0.005,
    .ctFitFlatRetryRefit = function(current, model, ...) list(
      fit = .retry_fit(loglik = -120), warnings = list(simpleWarning('refit'))))
  fit <- .retry_fit(directions = list(.retry_dir('trend_drift')))
  out <- suppressMessages(ctsem:::.ctFitFlatRetry(fit, list(simpleWarning('original')),
    cores = 1L))
  expect_equal(out$fit$estimate$loglik, -100)
  expect_identical(vapply(out$warnings, conditionMessage, ''), 'original')
  expect_false(out$fit$flatretry$kept)
})

test_that("the retry is a julia control", {
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = FALSE), 'julia'))
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = FALSE), 'stan'))
  expect_error(ctsem:::.ctFitCheckControls(list(flatretry = TRUE), 'stan'), 'flatretry')
})
