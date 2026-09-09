# ctFitCheck() -- the whole dashboard -- on both backends, plus the two
# delegating switches marginalcovcheck and trajectoryplot.
#
# marginalcovcheck calls the existing ctFitCheckCov() machinery at lag 0, split
# by a per-subject early/late median split of the time variable (each
# subject's own first half against its own second half -- not a split on
# absolute time across the sample, which would put a whole subject's
# trajectory on one side under unbalanced start times). trajectoryplot reuses
# ctPredict()/ctKalman() for the implied (Kalman smoother) trajectory and its
# uncertainty; the computation is factored into .ctFitCheckTrajectory() so it
# can be asserted on directly rather than only checking the plotting call
# doesn't error.
#
# Every panel now reads the fit through the backend-agnostic accessors, so the
# whole function works for a julia fit. The one exception is the prior
# predictive: ctGenerateFromPriors(), which supplies the `$priorpred` that
# panel reads, is stan-only, so on a julia fit that panel is skipped with a
# message rather than failing the call.
#
# ctCheckFit / ctFitCovCheck are the pre-3.12 names, kept as aliases.

test_that("the old names are still exported and identical to the new ones", {
  expect_identical(ctCheckFit, ctFitCheck)
  expect_identical(ctFitCovCheck, ctFitCheckCov)
})

test_that("marginalcovcheck runs ctFitCheckCov at lag 0 with a per-subject early/late split (stan)", {
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

  res <- suppressWarnings(suppressMessages(ctFitCheckCov(ctstantestfit, cor = TRUE,
    plot = FALSE, data = covdat, splitby = ".period", splitdata = covdat,
    split = "factor", lags = 0, nsamples = 5)))

  expect_true(nrow(res) > 0)
  expect_true(all(res$lag == 0))
  expect_true(all(res$split %in% c("Early", "Late")))
  expect_true("Sig" %in% names(res))
  expect_true(all(c("Y1", "Y2") %in% res$rowvar))

  # And the same thing, reached through ctFitCheck()'s switch, which also hands
  # back the data it plotted rather than only printing.
  withr::local_pdf(NULL)
  out <- suppressWarnings(suppressMessages(ctFitCheck(ctstantestfit, data = FALSE,
    postpred = FALSE, marginalcovcheck = TRUE, nsamples = 5)))
  expect_true(all(out$marginalcov$lag == 0))
  expect_setequal(unique(out$marginalcov$split), c("Early", "Late"))
  expect_true(length(out$plots) > 0)
  expect_true(all(vapply(out$plots, inherits, logical(1), "ggplot")))
})

test_that(".ctFitCheckTrajectory returns observed and implied means with a positive-width band (stan)", {
  skip_on_cran()

  traj <- suppressMessages(ctsem:::.ctFitCheckTrajectory(ctstantestfit, by = "time", breaks = 4))

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
  out <- suppressWarnings(suppressMessages(ctFitCheck(ctstantestfit, data = FALSE,
    postpred = FALSE, trajectoryplot = TRUE)))
  expect_named(out$trajectory, c("observed", "implied"))
})

test_that("ctFitCheck() defaults are unchanged and it returns its panels invisibly (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  expect_identical(formals(ctFitCheck)$marginalcovcheck, FALSE)
  expect_identical(formals(ctFitCheck)$trajectoryplot, FALSE)

  # Printing is still the default behaviour, and the value is invisible -- so
  # `ctFitCheck(fit)` at the console shows plots and prints no list.
  expect_false(withVisible(suppressWarnings(suppressMessages(
    ctFitCheck(ctstantestfit))))$visible)

  out <- suppressWarnings(suppressMessages(ctFitCheck(ctstantestfit, covplot = TRUE,
    breaks = 2, nsamples = 5, statepred = TRUE, residuals = TRUE, fastcov = TRUE)))
  expect_named(out, c("plots", "data", "covariances", "marginalcov", "trajectory"))
  expect_true(all(vapply(out$plots, inherits, logical(1), "ggplot")))
  expect_true(all(c("Data", "StatePred", "Residuals") %in% names(out$covariances)))
  expect_true(nrow(out$data) > 0)
  expect_true(all(c("DataSource", "Sample", "WhichObs") %in% names(out$data)))
})

# The covariance heatmaps have two routes to a covariance matrix, and the
# default one is covml(). Its return element was renamed cp -> estimate in
# 2024 and this caller was not updated, so the default path read NULL: with
# corr=TRUE that surfaced as cov2cor()'s "'V' is not a square numeric matrix",
# and with corr=FALSE it was silent -- assigning NULL to a list element
# *removes* it, so corlist stayed empty, every downstream seq_along() loop ran
# zero times, and the function returned no covariances and no plots without
# complaint. Every existing covplot test passed fastcov=TRUE, which is why two
# years went by. So assert on the default route, and on the list being
# populated rather than only on the call not erroring.
test_that("ctFitCheck() covariance panels work on the default covml route, corr either way (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  expect_identical(formals(ctFitCheck)$fastcov, FALSE)
  expect_identical(formals(ctFitCheck)$corr, TRUE)

  for (cr in c(TRUE, FALSE)) {
    out <- suppressWarnings(suppressMessages(ctFitCheck(ctstantestfit, covplot = TRUE,
      breaks = 2, nsamples = 5, data = TRUE, postpred = TRUE, corr = cr)))
    expect_true(all(c("Data", "PostPred") %in% names(out$covariances)))
    expect_length(out$plots, 3L) # one heatmap per source, plus their difference
    for (cm in out$covariances) {
      expect_true(is.matrix(cm))
      expect_identical(nrow(cm), ncol(cm))
      if (cr) expect_equal(diag(cm), setNames(rep(1, nrow(cm)), rownames(cm)))
    }
  }
})

test_that("ctFitMelt reads the same one-step-ahead predictions the backend's own filter gives (stan)", {
  skip_on_cran()

  # The Residuals data source is Y - yprior. ctKalmanArray()'s errprior is the
  # same quantity by another route, and .ctFitMeltStates() is allowed to take
  # either route depending on what the fit carries, so they have to agree.
  melted <- suppressWarnings(suppressMessages(ctFitMelt(ctstantestfit, maxsamples = "all")))
  resid <- melted[melted$DataSource == "Residuals", ]
  k <- suppressWarnings(suppressMessages(ctKalmanArray(ctstantestfit, pointest = TRUE,
    subjectpars = FALSE)))
  errprior <- k$errprior[1, , ]
  for (vi in seq_along(ctstantestfit$ctstanmodelbase$manifestNames)) {
    v <- ctstantestfit$ctstanmodelbase$manifestNames[vi]
    expect_equal(unname(resid[[v]]), unname(errprior[, vi]), tolerance = 1e-8)
  }
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

test_that("the whole dashboard runs on a julia fit, and its panels match the backend's own filter", {
  skip_without_julia()
  withr::local_pdf(NULL)

  fit <- suppressWarnings(suppressMessages(ctFit(.checkfit_julia_data(), .checkfit_julia_model(),
    backend = "julia", verbose = 0)))
  expect_s3_class(fit, "ctJuliaFit")

  # This is the call that used to stop() with "not available for julia backend
  # fits". Nothing about it is stan-shaped any more.
  out <- suppressWarnings(suppressMessages(ctFitCheck(fit, nsamples = 5)))
  expect_true(length(out$plots) > 0)
  expect_true(all(vapply(out$plots, inherits, logical(1), "ggplot")))
  # data and postpred are the defaults; statepred and residuals are not, and
  # ctFitCheck() filters the melted data by those switches.
  expect_setequal(unique(as.character(out$data$DataSource)), c("Data", "PostPred"))

  # Covariance heatmaps, and the two delegating switches.
  outcov <- suppressWarnings(suppressMessages(ctFitCheck(fit, covplot = TRUE,
    breaks = 2, nsamples = 5, statepred = TRUE, residuals = TRUE, fastcov = TRUE)))
  expect_true(all(c("Data", "StatePred", "Residuals", "PostPred") %in%
      names(outcov$covariances)))
  outsw <- suppressWarnings(suppressMessages(ctFitCheck(fit, data = FALSE,
    postpred = FALSE, marginalcovcheck = TRUE, trajectoryplot = TRUE,
    breaks = 3, nsamples = 5)))
  expect_true(all(outsw$marginalcov$lag == 0))
  expect_named(outsw$trajectory, c("observed", "implied"))
  expect_true(all(outsw$trajectory$implied$se > 0))

  # The Residuals panel is Y - yprior from the julia engine's own forward pass,
  # so it must equal ctBackendKalman()'s errprior exactly rather than
  # approximately -- both come from .ctFitMeltStates().
  melted <- suppressWarnings(suppressMessages(ctFitMelt(fit, maxsamples = "all")))
  resid <- melted[melted$DataSource == "Residuals", ]
  k <- suppressWarnings(suppressMessages(ctBackendKalman(fit, pointest = TRUE)))
  expect_equal(unname(resid$Y1), as.numeric(k$errprior[1, , 1]), tolerance = 1e-10)

  # The prior predictive is the one panel with no julia route: it needs
  # ctGenerateFromPriors(), which is stan-only. Named, not silently dropped.
  expect_message(suppressWarnings(ctFitCheck(fit, priorpred = TRUE, data = FALSE,
    postpred = FALSE, statepred = TRUE, nsamples = 5)),
    "ctGenerateFromPriors")
})
