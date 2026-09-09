# ctFitCheckCov() -- the lagged covariance/correlation diagnostic.
#
# It is reachable two ways: on its own, which is the way that takes a `splitby`
# group, and as the LaggedCovariance panel of ctPostPredPlots(), which is
# tested in test-postpred-panels.R. The ctFitCheck() dashboard that used to
# wrap it was removed in 3.12.0.
#
# ctFitCovCheck is the pre-3.12 name, kept as an alias.

test_that("the old name is still exported and identical to the new one", {
  expect_identical(ctFitCovCheck, ctFitCheckCov)
})

test_that("ctFitCheckCov() splits at lag 0 by a per-subject early/late split (stan)", {
  skip_on_cran()

  ctmb <- ctstantestfit$ctstanmodelbase
  idname <- ctmb$subjectIDname
  timename <- ctmb$timeName
  covdat <- data.table::as.data.table(ctsem:::.ctFitLongData(ctstantestfit))
  covdat[, .period := ifelse(get(timename) <= stats::median(get(timename), na.rm = TRUE),
    "Early", "Late"), by = idname]

  # Every subject actually got split into both halves -- otherwise the split
  # is not doing what "early vs late" is supposed to mean. This is a
  # per-subject median split, not a split on absolute time across the sample:
  # with unbalanced start times the latter puts a whole subject's trajectory
  # on one side.
  expect_true(all(c("Early", "Late") %in% covdat$.period))

  res <- suppressWarnings(suppressMessages(ctFitCheckCov(ctstantestfit, cor = TRUE,
    plot = FALSE, data = covdat, splitby = ".period", splitdata = covdat,
    split = "factor", lags = 0, nsamples = 5)))

  expect_true(nrow(res) > 0)
  expect_true(all(res$lag == 0))
  expect_true(all(res$split %in% c("Early", "Late")))
  expect_true("Sig" %in% names(res))
  expect_true(all(c("Y1", "Y2") %in% res$rowvar))
})

test_that("ctFitCheckCov(plot=TRUE) gives one ggplot per row variable (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)

  gg <- suppressWarnings(suppressMessages(ctFitCheckCov(ctstantestfit, cor = TRUE,
    lags = 0:2, nsamples = 5, plot = TRUE)))
  expect_setequal(names(gg), c("Y1", "Y2"))
  expect_true(all(vapply(gg, inherits, logical(1), "ggplot")))
})
