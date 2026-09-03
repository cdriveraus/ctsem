# ctCheckFit(marginalcovcheck=TRUE) and ctCheckFit(trajectoryplot=TRUE).
#
# marginalcovcheck calls the existing ctFitCovCheck() machinery at lag 0, split
# by a per-subject early/late median split of the time variable (each
# subject's own first half against its own second half -- not a split on
# absolute time across the sample, which would put a whole subject's
# trajectory on one side under unbalanced start times). trajectoryplot reuses
# ctPredict()/ctKalman() for the implied (Kalman smoother) trajectory and its
# uncertainty; the computation is factored into .ctCheckFitTrajectory() so it
# can be asserted on directly rather than only checking the plotting call
# doesn't error.
#
# Both switches are built on machinery that already works for stan and julia
# fits (ctFitCovCheck(), ctPredict()/ctKalman()), so -- unlike the rest of
# ctCheckFit(), which still refuses julia fits -- they work on either backend
# as long as every other plot switch is off.

test_that("marginalcovcheck runs ctFitCovCheck at lag 0 with a per-subject early/late split (stan)", {
  skip_on_cran()

  ctmb <- ctstantestfit$ctstanmodelbase
  idname <- ctmb$subjectIDname
  timename <- ctmb$timeName
  covdat <- data.table::as.data.table(ctsem:::.ctFitLongData(ctstantestfit))
  covdat[, .period := ifelse(get(timename) <= stats::median(get(timename), na.rm = TRUE),
    "Early", "Late"), by = idname]

  # Every subject actually got split into both halves -- otherwise the split
  # is not doing what "early vs late" is supposed to mean.
  expect_true(all(c("Early", "Late") %in% covdat$.period))

  res <- suppressWarnings(suppressMessages(ctFitCovCheck(ctstantestfit, cor = TRUE,
    plot = FALSE, data = covdat, splitby = ".period", splitdata = covdat,
    split = "factor", lags = 0, nsamples = 5)))

  expect_true(nrow(res) > 0)
  expect_true(all(res$lag == 0))
  expect_true(all(res$split %in% c("Early", "Late")))
  expect_true("Sig" %in% names(res))
  expect_true(all(c("Y1", "Y2") %in% res$rowvar))

  # And the same thing, reached through ctCheckFit()'s new switch.
  withr::local_pdf(NULL)
  expect_no_error(ctCheckFit(ctstantestfit, data = FALSE, postpred = FALSE,
    marginalcovcheck = TRUE, nsamples = 5))
})

test_that(".ctCheckFitTrajectory returns observed and implied means with a positive-width band (stan)", {
  skip_on_cran()

  traj <- suppressMessages(ctsem:::.ctCheckFitTrajectory(ctstantestfit, by = "time", breaks = 4))

  expect_true(all(c("Row", ".TimeBin", "TimeMid", "Mean") %in% names(traj$observed)))
  expect_true(all(c("Row", ".TimeBin", "TimeMid", "Mean", "se", "q025", "q975") %in% names(traj$implied)))
  expect_setequal(unique(traj$observed$Row), c("Y1", "Y2"))
  expect_setequal(unique(traj$implied$Row), c("Y1", "Y2"))

  # se comes straight from the Kalman smoother covariance (see the function's
  # comment) -- it should be strictly positive, and the band it derives should
  # bracket the point estimate.
  expect_true(all(traj$implied$se > 0))
  expect_true(all(traj$implied$q975 > traj$implied$q025))
  expect_equal(traj$implied$q025, traj$implied$Mean - 1.96 * traj$implied$se)
  expect_equal(traj$implied$q975, traj$implied$Mean + 1.96 * traj$implied$se)

  withr::local_pdf(NULL)
  expect_no_error(ctCheckFit(ctstantestfit, data = FALSE, postpred = FALSE,
    trajectoryplot = TRUE))
})

test_that("ctCheckFit() defaults are unchanged and still refuse julia fits for the classic plots", {
  skip_on_cran()
  withr::local_pdf(NULL)
  expect_no_error(ctCheckFit(ctstantestfit))
  expect_identical(formals(ctCheckFit)$marginalcovcheck, FALSE)
  expect_identical(formals(ctCheckFit)$trajectoryplot, FALSE)
})

.checkfit_julia_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.checkfit_julia_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5)
  drift <- -0.8; diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(15), function(i) {
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

test_that("marginalcovcheck and trajectoryplot both work on a julia fit, while the classic gate still blocks it", {
  skip_on_cran()
  skip_without_julia()
  withr::local_pdf(NULL)

  fit <- suppressWarnings(suppressMessages(ctFit(.checkfit_julia_data(), .checkfit_julia_model(),
    backend = "julia", verbose = 0)))
  expect_s3_class(fit, "ctJuliaFit")

  # The classic ctFitMelt-based plots still refuse a julia fit ...
  err <- expect_error(ctCheckFit(fit))
  expect_match(conditionMessage(err), "julia backend fits")

  # ... but the two new switches, with every classic switch off, do not.
  expect_no_error(ctCheckFit(fit, data = FALSE, postpred = FALSE,
    marginalcovcheck = TRUE, nsamples = 5))
  expect_no_error(ctCheckFit(fit, data = FALSE, postpred = FALSE,
    trajectoryplot = TRUE, breaks = 3))

  # And the underlying computations return sensible, backend-agnostic values.
  mcov <- suppressWarnings(suppressMessages(ctFitCovCheck(fit, cor = TRUE,
    plot = FALSE, lags = 0, nsamples = 5)))
  expect_true(nrow(mcov) > 0)
  expect_true(all(mcov$lag == 0))

  traj <- suppressMessages(ctsem:::.ctCheckFitTrajectory(fit, by = "time", breaks = 3))
  expect_true(nrow(traj$observed) > 0)
  expect_true(nrow(traj$implied) > 0)
  expect_true(all(traj$implied$se > 0))
  expect_true(all(traj$implied$q975 > traj$implied$q025))
})
