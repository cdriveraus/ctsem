
  library(ctsem)
  library(testthat)


test_that("ctCoverageCheck generation cores default to fitCores", {
  expect_identical(
    formals(ctCoverageCheck)$generateCores,
    quote(fitCores)
  )
})

# ctCoverageCheck() used to read fit$stanfit$rawest and
# fit$stanfit$rawposterior directly (truepars, and the quantiles it compares
# them against), and named parameters via ctFitgetparnamesfromraw(), which
# reads fit$setup$matsetup -- none of which a julia fit carries, so a
# backend='julia' entry in fitArgs would have errored partway through the
# first iteration. It now reads through .ctFitRawEstimate()/
# .ctFitRawPosterior()/.ctFitRawParNames() (R/ctBackendSummary.R), which
# branch on class(fit). This had no fitting test before this file (the one
# above only inspects formals()).

.coverage_julia_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.coverage_julia_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8; diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(20), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
  data
}

test_that("ctCoverageCheck runs end to end for a julia fit", {
  skip_without_julia()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("gridExtra")
  withr::local_pdf(NULL)

  dat <- .coverage_julia_data()
  model <- .coverage_julia_model()

  res <- suppressWarnings(suppressMessages(ctCoverageCheck(
    initialData = dat, fittingModel = model, niter = 1,
    fitArgs = list(hess = list(backend = "julia", verbose = 0)),
    cores = 1, fitCores = 1, plotEvery = 1)))

  expect_true(is.list(res))
  expect_true(nrow(res$results) > 0)
  expect_identical(res$successful_iterations, 1L)
  expect_true(all(c("truepars", "X2.5.", "X50.", "X97.5.", "par", "coverage") %in%
    names(res$results)))
  # par names come from .ctFitRawParNames(), in the same order as truepars
  # (.ctFitRawEstimate()) -- five free parameters in this model.
  expect_setequal(unique(res$results$par), c("drift", "diff", "mvar", "mmean", "t0v"))
  expect_s3_class(res$bias_plot, "ggplot")
  expect_s3_class(res$coverage_plot, "ggplot")
})
