# A model parameter a stage left run off where the likelihood is flat is held
# at the other side of its range while every other coordinate starts again
# from the fit's start, and is then freed (`.ctBackendRunOff()`,
# `.ctBackendStallEscape()`, `.ctBackendCorrectResult()`). Fit-free: engine
# results built by hand, the optimiser and the certification mocked.

.runoff_names <- c('trend_drift', 'drift_x', 'cint_trend', 'popsd_drift_x')

# An engine result whose finish walked the flat directions `flat` (each a
# vector over the coordinates), at `x`.
.runoff_result <- function(x, flat = list(), saturated = 0L, loglik = -100) {
  n <- length(x)
  vectors <- diag(n)
  isflat <- rep(FALSE, n)
  for (j in seq_along(flat)) {
    vectors[, j] <- flat[[j]] / sqrt(sum(flat[[j]]^2))
    isflat[j] <- TRUE
  }
  list(minimizer = x, maximum_loglik = loglik, saturated_parameters = saturated,
    newton_flat_ran = TRUE, newton_flat_at = x, newton_flat_eigenvalues = numeric(n),
    newton_flat_eigenvectors = vectors, newton_flat_candidates = seq_len(n),
    newton_flat_change = numeric(n), newton_flat_flat = isflat,
    newton_flat_gain = 0, newton_flat_point = 0, newton_flat_direction = 0L,
    newton_flat_evaluations = 0L, newton_flat_complete = TRUE)
}

# The auto-effects end at zero, the sd at zero; the intercept has no end.
.runoff_mocks <- function(usable = TRUE, env = parent.frame()) testthat::local_mocked_bindings(
  .ctBackendRawParameterNames = function(fit, npar) .runoff_names[seq_len(npar)],
  .ctJuliaSaturatedNames = function(result, model_spec, npar, field)
    .runoff_names[as.integer(result[[field]])],
  .ctBackendRawEnds = function(model_spec, k, name = NA) if (k == 3L) c(FALSE, FALSE) else c(FALSE, TRUE),
  .ctDataStartUsable = function(...) usable,
  .env = env)

test_that("a flat, run-off parameter with an end qualifies", {
  .runoff_mocks()
  x <- c(-16, 0.4, 9, -20)
  # trend_drift flat and far out; drift_x flat but near the start; cint_trend
  # far out but not flat; the random effect's sd flat and at zero.
  r <- .runoff_result(x, flat = list(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1)))
  expect_identical(ctsem:::.ctBackendRunOff(r, list(), 4L), 1L)
  expect_identical(ctsem:::.ctBackendRunOff(r, list(), 4L, tried = 1L), 4L)
  expect_length(ctsem:::.ctBackendRunOff(r, list(), 4L, tried = c(1L, 4L)), 0L)
  # A saturated transform leads without a flat direction.
  expect_identical(ctsem:::.ctBackendRunOff(.runoff_result(c(-16, 5, 9, -20),
    saturated = 2L), list(), 4L), 2L)
  # Two directions spanning the drift and the intercept: both have all their
  # variation there, but the intercept has no end to be held at the other
  # side of.
  r <- .runoff_result(x, flat = list(c(0.8, 0, 0.6, 0), c(-0.6, 0, 0.8, 0)))
  expect_identical(ctsem:::.ctBackendRunOff(r, list(), 4L), 1L)
  expect_length(ctsem:::.ctBackendRunOff(r, list(), 4L, tried = 1L), 0L)
})

test_that("a transform's ends are read off its text", {
  spec <- list(parameter_table = data.frame(parnumber = c(1L, 2L, 3L),
    transform = c('-(1e-6 + 2 * log1p_exp(-2 * param[1]))', 'param[2] * 10',
      'exp(param[3])'), stringsAsFactors = FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 1L), c(FALSE, TRUE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 2L), c(FALSE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 3L), c(TRUE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 9L), c(FALSE, FALSE))
  # A random effect's sd or own variance is the engine's, zero below.
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 9L, 'popsd_drift'), c(TRUE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 9L, 'popunique_cint'), c(TRUE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 9L, 'poploading_cint_dim1'), c(FALSE, FALSE))
})

test_that("the escape point holds the parameter at the other side and restarts the rest", {
  .runoff_mocks()
  r <- .runoff_result(c(-16, 0.4, 9, -20), flat = list(c(0.8, 0, 0.6, 0)))
  on <- list(flatretry = TRUE)
  p <- ctsem:::.ctBackendStallEscape(r, on, list(), restart = c(0.1, 0.2, 0.3, 0.4))
  expect_identical(attr(p, 'kind'), 'runoff')
  expect_identical(attr(p, 'coordinates'), 1L)
  expect_equal(as.numeric(p), c(ctsem:::.ctDataStartRawBound, 0.2, 0.3, 0.4))
  expect_null(attr(p, 'unusable'))
  # Without a start to go back to, from where the stage ended.
  expect_equal(as.numeric(ctsem:::.ctBackendStallEscape(r, on, list())),
    c(ctsem:::.ctDataStartRawBound, 0.4, 9, -20))
  # Off unless asked for, and a measured better point goes first.
  expect_null(ctsem:::.ctBackendStallEscape(r, list(), list()))
  expect_null(ctsem:::.ctBackendStallEscape(r, list(flatretry = FALSE), list()))
  o <- r
  o$overshot <- TRUE
  o$overshoot_point <- c(-15, 0.4, 9, -20)
  o$overshoot_parameters <- 1L
  o$overshoot_gain <- 1
  expect_null(attr(ctsem:::.ctBackendStallEscape(o, on, list()), 'kind'))
})

test_that("a point the optimiser cannot start from is marked, for the caller to skip", {
  .runoff_mocks(usable = FALSE)
  r <- .runoff_result(c(-16, 0.4, 9, -20), flat = list(c(1, 0, 0, 0)))
  p <- ctsem:::.ctBackendStallEscape(r, list(flatretry = TRUE), list())
  expect_true(attr(p, 'unusable'))
  expect_identical(attr(p, 'coordinates'), 1L)
})

# The correction loop with its optimiser and certification mocked: `stages`
# are the results the optimiser returns in turn.
.runoff_loop <- function(first, stages, escape) {
  calls <- list()
  i <- 0L
  optimise <- function(from, endgame) {
    calls[[length(calls) + 1L]] <<- from
    i <<- i + 1L
    stages[[i]]
  }
  testthat::local_mocked_bindings(
    .ctJuliaRunCounts = function(r) c(newton_hessians = 0),
    .ctBackendEndgameOf = function(r, npar) NULL,
    .ctJuliaHoldModes = function(...) NULL,
    .ctBackendEndgameAt = function(...) NULL)
  said <- testthat::capture_messages(out <- ctsem:::.ctBackendCorrectResult(first,
    list(), npar = 4L, optimise = optimise, escape = escape, maxtries = 0L))
  list(out = out, calls = calls, said = said)
}

test_that("a run-off escape that wins is kept, and one that loses leaves the stage", {
  .runoff_mocks()
  first <- .runoff_result(c(-16, 0.4, 9, -20), flat = list(c(1, 0, 0, 0)))
  tried <- integer()
  escape <- function(r, runoffonly = FALSE) {
    p <- ctsem:::.ctBackendStallEscape(r, list(flatretry = TRUE), list(), tried = tried,
      runoffonly = runoffonly)
    if (identical(attr(p, 'kind'), 'runoff')) tried <<- c(tried, attr(p, 'coordinates'))
    p
  }
  # Freed, the drift runs off to the other end, flat again, and is not sent back.
  better <- .runoff_result(c(12, 0.4, 0.1, -20), flat = list(c(1, 0, 0, 0)), loglik = -80)
  run <- .runoff_loop(first, list(better), escape)
  # The escape runs with no budget for stall escapes or resumes (maxtries 0).
  expect_length(run$calls, 1L)
  expect_identical(attr(run$calls[[1]], 'kind'), 'runoff')
  expect_equal(as.numeric(run$out$result$maximum_loglik), -80)
  expect_identical(run$out$escapes, 1L)
  expect_true(any(grepl('trend_drift ran off', run$said)))
  expect_identical(run$out$corrections[[1]]$status, 'runoff')

  tried <- integer()
  worse <- .runoff_result(c(12, 0.4, 0.1, -20), loglik = -120)
  run <- .runoff_loop(first, list(worse), escape)
  expect_length(run$calls, 1L)
  expect_equal(as.numeric(run$out$result$maximum_loglik), -100)
  expect_false(run$out$corrections[[1]]$resumed)
  expect_equal(run$out$corrections[[1]]$loglik, -120)
  expect_length(run$said, 0L)
})

test_that("the retry is a julia control", {
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = FALSE), 'julia'))
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = TRUE), 'julia'))
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = FALSE), 'stan'))
  expect_error(ctsem:::.ctFitCheckControls(list(flatretry = TRUE), 'stan'), 'flatretry')
})

test_that("a fit that did not run the escape names every parameter it would hold", {
  .runoff_mocks()
  x <- c(-16, 0.4, 9, -20)
  r <- .runoff_result(x, flat = list(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1)))
  expect_identical(ctsem:::.ctBackendFlatRetryCandidates(r, list(), 4L), c(1L, 4L))
  expect_length(ctsem:::.ctBackendFlatRetryCandidates(.runoff_result(x), list(), 4L), 0L)
  one <- ctsem:::.ctBackendFlatRetryMessage('popsd_cint')
  expect_match(one, '^popsd_cint ran off to where the data say little about it[.]')
  expect_match(one, 'ctFitFlatRetry(fit)', fixed = TRUE)
  expect_match(one, 'its range')
  several <- ctsem:::.ctBackendFlatRetryMessage(c('a', 'b', 'c'))
  expect_match(several, '^a, b and c ran off to where the data say little about them[.]')
  expect_match(several, "each one's range", fixed = TRUE)
})

# A fit as ctFitFlatRetry() reads it, and ctFitUpdate() mocked to return
# `refit`, recording what it was asked.
.flatretry_fit <- function(names) structure(list(optim = list(
  flatretry_parameters = names), args = list(input = list(optimcontrol =
  list(stallretries = 1L)))), class = c('ctJuliaFit', 'ctFit'))

test_that("ctFitFlatRetry refits with the escape on and keeps only a win", {
  asked <- NULL
  refit <- NULL
  testthat::local_mocked_bindings(
    .ctFitIsJulia = function(fit) TRUE,
    .ctFitReplayArgs = function(fit) fit$args$input,
    ctFitUpdate = function(oldfit, ...) {
      asked <<- list(...)
      message('from the refit')
      refit
    })
  fit <- .flatretry_fit('popsd_cint')
  # Kept: the new fit, and what its refit said.
  refit <- list(optim = list(corrections = list(list(status = 'runoff', resumed = TRUE))))
  said <- testthat::capture_messages(out <- ctFitFlatRetry(fit, cores = 3L))
  expect_identical(out, refit)
  expect_true(asked$refit)
  expect_true(asked$optimcontrol$flatretry)
  expect_identical(asked$optimcontrol$stallretries, 1L)
  expect_identical(asked$cores, 3L)
  expect_true(any(grepl('from the refit', said)))
  # Not kept: the fit as it was, and the refit's messages left out.
  refit <- list(optim = list(corrections = list(list(status = 'runoff', resumed = FALSE))))
  said <- testthat::capture_messages(out <- ctFitFlatRetry(fit))
  expect_identical(out, fit)
  expect_false(any(grepl('from the refit', said)))
  expect_true(any(grepl('popsd_cint did not raise the log likelihood', said)))
  expect_null(asked$cores)
  # Nothing run off: no refit.
  asked <- NULL
  said <- testthat::capture_messages(out <- ctFitFlatRetry(.flatretry_fit(character())))
  expect_null(asked)
  expect_true(any(grepl('returned unchanged', said)))
})
