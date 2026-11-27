# ctsemtestfit is a julia fit shipped as data. The examples summarise it where
# Julia is missing (CRAN), which works because a fit keeps its draws, and its
# estimate, pushed through the model's transforms (`fit$transformedpars`).
# When the engine's transforms or layout change, the stored values stop
# matching what the engine gives, and the fit must be made again:
# ctsem:::ctdataupdate("juliafit").

# Every way into Julia fails. `startable` says whether Julia counts as
# installed: FALSE is CRAN, TRUE a machine whose session has not started.
local_no_julia <- function(startable = FALSE, env = parent.frame()) {
  nojulia <- function(...) stop("Julia was asked for.", call. = FALSE)
  testthat::local_mocked_bindings(.ctJuliaStartable = function(...) startable,
    .ctJuliaSessionRunning = function(...) FALSE, .ctJuliaModule = nojulia,
    .ctJuliaCall = nojulia, .ctJuliaEval = nojulia, .package = "ctsem",
    .env = env)
}

test_that("the example fit's reports need no Julia", {
  local_no_julia()
  s <- summary(ctsemtestfit)
  expect_false("mode" %in% colnames(s$popmeans))
  # A single parameter's median under the fit's normal draws is its estimate.
  expect_equal(s$popmeans["drift_eta1", "50%"],
    unname(ctSummaryMatrices(ctsemtestfit)$estimate$DRIFT[1, 1]), tolerance = 1e-3)
  r2 <- ctDerived(ctsemtestfit, function(m) m$DRIFT[1, 1])
  expect_equal(summary(r2)[["50%"]], s$popmeans["drift_eta1", "50%"],
    tolerance = 1e-3)
  dp <- ctDiscretePars(ctsemtestfit, times = c(0, 1, 2), nsamples = 20)
  expect_equal(dim(dp), c(20L, 1L, 3L, 2L, 2L))
  expect_equal(nrow(summary(dp)), 16L)
  e <- ctExtract(ctsemtestfit, nsamples = 10)
  expect_equal(dim(e$pop_DRIFT), c(10L, 2L, 2L))
  expect_no_error(ctModelLatex(ctsemtestfit, compile = FALSE, open = FALSE))
  expect_no_error(ctNetwork(ctsemtestfit))
})

test_that("with Julia installed but not started, the reports do not start it", {
  local_no_julia(startable = TRUE)
  expect_no_error(ctDiscretePars(ctsemtestfit, times = c(0, 1), nsamples = 5))
  expect_no_error(ctExtract(ctsemtestfit, nsamples = 5))
  expect_no_error(summary(ctsemtestfit))
})

test_that("the example fit's stored values are the engine's", {
  skip_without_julia()
  stale <- paste("ctsemtestfit no longer matches the engine; make it again with",
    "ctsem:::ctdataupdate('juliafit').")
  tp <- ctsemtestfit$transformedpars
  # The estimate's entry, and the draws read from the cache by row, against
  # the engine's own constrain at those points.
  fresh <- .ctBackendConstrain(ctsemtestfit, tp$atestimate$samples)
  expect_equal(tp$atestimate[names(fresh)], fresh, tolerance = 1e-6, info = stale)
  rows <- c(1L, 500L)
  fresh <- .ctBackendConstrain(ctsemtestfit, tp$samples[rows, , drop = FALSE])
  expect_equal(.ctBackendConstrainRows(tp, rows)[names(fresh)], fresh,
    tolerance = 1e-6, info = stale)
  # And the engine still takes the stored fit.
  k <- ctKalman(ctsemtestfit, subjects = 1)
  expect_gt(nrow(k), 0)
})
