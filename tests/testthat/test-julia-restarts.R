# Random restarts for a julia fit that did not converge, and what a
# not-a-maximum verdict says. See R/ctBackendRestarts.R.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("restarts run only for a fit that is not converged and not asked for", {
  notmax <- list(certified = FALSE, status = "notmaximum")
  wanted <- function(...) ctsem:::.ctBackendRestartsWanted(...)
  expect_equal(wanted(list(converged = FALSE), notmax, list(), NULL, TRUE), 5L)
  expect_equal(wanted(list(converged = FALSE), notmax, list(restarts = 2), NULL, TRUE), 2L)
  # A certified fit, and a maximum some coordinate of which is undetermined,
  # are left alone; so is a state-explicit fit.
  expect_equal(wanted(list(), list(certified = TRUE, status = "certified"),
    list(), NULL, TRUE), 0L)
  expect_equal(wanted(list(), list(certified = FALSE, status = "unidentified"),
    list(), NULL, TRUE), 0L)
  expect_equal(wanted(list(converged = FALSE), notmax, list(), NULL, FALSE), 0L)
  # The user chose the start or capped the iterations: a non-converged fit is
  # then the fit that was asked for, and moving it would change the answer.
  expect_equal(wanted(list(converged = FALSE), notmax, list(), c(0, 0), TRUE), 0L)
  expect_equal(wanted(list(converged = FALSE), notmax, list(maxiter = 5), NULL, TRUE), 0L)
  expect_equal(wanted(list(converged = FALSE), notmax, list(restarts = 0), NULL, TRUE), 0L)
  # Without a certification the optimiser's own verdict decides.
  expect_equal(wanted(list(converged = TRUE), NULL, list(), NULL, TRUE), 0L)
  expect_equal(wanted(list(converged = FALSE), NULL, list(), NULL, TRUE), 5L)
})

test_that("a not-a-maximum message names the direction and the remedy", {
  fit <- list(model_spec = list(
    parameter_table = data.frame(param = c("drift_eta1", "cint1", "cint2"),
      parnumber = 1:3, stringsAsFactors = FALSE),
    random_effects = data.frame(parameter = c(4L, 5L), param = c("cint1",
      "cint2__cint1"), type = c("sd", "correlation"), stringsAsFactors = FALSE)),
    optim = list())
  # The rising direction is mostly the correlation and the scale.
  covariance <- list(status = "notmaximum",
    negative_vector = c(0.05, 0, 0, 0.6, 0.8))
  message <- ctsem:::.ctBackendNotMaximumMessage(fit, covariance)
  expect_match(message, "rawcor_cint2__cint1", fixed = TRUE)
  expect_match(message, "popsd_cint1", fixed = TRUE)
  expect_false(grepl("drift_eta1", message, fixed = TRUE))
  expect_match(message, "poprank = 1", fixed = TRUE)
  # A direction outside the covariance gets no poprank advice.
  fixed <- list(status = "notmaximum", negative_vector = c(1, 0.9, 0, 0, 0))
  message <- ctsem:::.ctBackendNotMaximumMessage(fit, fixed)
  expect_false(grepl("poprank", message, fixed = TRUE))
  # And it says when restarts were tried.
  fit$optim$restarts <- data.frame(start = 1:5)
  expect_match(ctsem:::.ctBackendNotMaximumMessage(fit, covariance),
    "5 random restarts found nothing better", fixed = TRUE)
})

test_that("restarts record every start and keep only a clear improvement", {
  skip_without_julia()
  set.seed(3)
  dat <- data.frame(id = rep(1:10, each = 6), time = rep(0:5, 10))
  dat$Y1 <- stats::rnorm(nrow(dat))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  spec <- suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "julia", fit = FALSE)))
  npar <- ctsem:::.ctBackendNpar(spec)
  set.seed(1)
  got <- ctsem:::.ctBackendRestarts(spec, rep(0, npar), npar, n = 2L,
    optimcontrol = list(), gradient = "adjoint", current = -Inf)
  expect_equal(nrow(got$table), 2L)
  expect_true(all(is.finite(got$table$logposterior)))
  expect_equal(sum(got$table$chosen), 1L)
  expect_false(is.null(got$best))
  expect_false(got$cancelled)
  # Against a current point better than any restart, none is kept.
  set.seed(1)
  none <- ctsem:::.ctBackendRestarts(spec, rep(0, npar), npar, n = 2L,
    optimcontrol = list(), gradient = "adjoint",
    current = max(got$table$logposterior) + 1)
  expect_null(none$best)
  expect_equal(sum(none$table$chosen), 0L)
})
