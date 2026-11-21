# Coordinates a stage left run off where the likelihood is flat: the one
# leading them is held at the other side of its range, those that ran off with
# it start again at zero, and the lead is then freed (`.ctBackendRunOff()`,
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
  .ctBackendRawEnds = function(model_spec, k) if (k == 3L) c(FALSE, FALSE) else c(FALSE, TRUE),
  .ctDataStartUsable = function(...) usable,
  .env = env)

test_that("the lead is a flat, run-off parameter with an end; what ran off with it follows", {
  .runoff_mocks()
  x <- c(-16, 0.4, 9, -20)
  # trend_drift flat and far out; drift_x flat but near the start; cint_trend
  # far out but not flat; the population sd far out and flat, but no lead.
  r <- .runoff_result(x, flat = list(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1)))
  expect_identical(ctsem:::.ctBackendRunOff(r, list(), 4L), list(lead = 1L, with = 4L))
  expect_null(ctsem:::.ctBackendRunOff(r, list(), 4L, tried = 1L))
  # A saturated transform leads without a flat direction.
  expect_identical(ctsem:::.ctBackendRunOff(.runoff_result(c(-16, 5, 9, -20),
    saturated = 2L), list(), 4L)$lead, 2L)
  # A direction shared by the drift and the intercept: the drift, at 0.64,
  # leads, and the intercept, at 0.36, ran off with it. The intercept has no
  # end to be held at the other side of, so it never leads.
  r <- .runoff_result(x, flat = list(c(0.8, 0, 0.6, 0)))
  expect_identical(ctsem:::.ctBackendRunOff(r, list(), 4L), list(lead = 1L, with = 3L))
  expect_null(ctsem:::.ctBackendRunOff(r, list(), 4L, tried = 1L))
})

test_that("a transform's ends are read off its text", {
  spec <- list(parameter_table = data.frame(parnumber = c(1L, 2L, 3L),
    transform = c('-(1e-6 + 2 * log1p_exp(-2 * param[1]))', 'param[2] * 10',
      'exp(param[3])'), stringsAsFactors = FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 1L), c(FALSE, TRUE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 2L), c(FALSE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 3L), c(TRUE, FALSE))
  expect_identical(ctsem:::.ctBackendRawEnds(spec, 9L), c(FALSE, FALSE))
})

test_that("the escape point holds the lead at the other side and restarts the rest", {
  .runoff_mocks()
  r <- .runoff_result(c(-16, 0.4, 9, -20), flat = list(c(0.8, 0, 0.6, 0)))
  p <- ctsem:::.ctBackendStallEscape(r, list(), list())
  expect_identical(attr(p, 'kind'), 'runoff')
  expect_identical(attr(p, 'coordinates'), 1L)
  expect_equal(as.numeric(p), c(ctsem:::.ctDataStartRawBound, 0.4, 0, -20))
  expect_null(attr(p, 'unusable'))
  # Off by its control, and a measured better point goes first.
  expect_null(ctsem:::.ctBackendStallEscape(r, list(flatretry = FALSE), list()))
  o <- r
  o$overshot <- TRUE
  o$overshoot_point <- c(-15, 0.4, 9, -20)
  o$overshoot_parameters <- 1L
  o$overshoot_gain <- 1
  expect_null(attr(ctsem:::.ctBackendStallEscape(o, list(), list()), 'kind'))
})

test_that("a point the optimiser cannot start from is marked, for the caller to skip", {
  .runoff_mocks(usable = FALSE)
  r <- .runoff_result(c(-16, 0.4, 9, -20), flat = list(c(1, 0, 0, 0)))
  p <- ctsem:::.ctBackendStallEscape(r, list(), list())
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
    p <- ctsem:::.ctBackendStallEscape(r, list(), list(), tried = tried,
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
  expect_true(any(grepl('trend_drift had run off', run$said)))
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
  expect_silent(ctsem:::.ctFitCheckControls(list(flatretry = FALSE), 'stan'))
  expect_error(ctsem:::.ctFitCheckControls(list(flatretry = TRUE), 'stan'), 'flatretry')
})
