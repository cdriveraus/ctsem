# ctFitCheckCov() -- the lagged covariance/correlation diagnostic.
#
# It is reachable two ways: on its own, which is the way that takes a `splitby`
# group, and as the LaggedCovariance panel of ctPostPredPlots(), which is
# tested in test-postpred-panels.R. The ctFitCheck() dashboard that used to
# wrap it was removed in 4.0.0.
#
# ctFitCovCheck is the pre-4.0 name, kept as an alias.

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

# The default view averages over observation indices, which hides misfit that
# differs by wave; x = 'obs' or a data column shows each wave or bin. Checked
# against cov() over each wave's subjects, computed here independently.
test_that("ctFitCheckCov(x = 'obs') is each wave's covariance, and 'lag' their mean (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)

  d <- data.table::as.data.table(ctsem:::.ctFitLongData(ctstantestfit))
  data.table::setorderv(d, c("id", "time"))
  d[, k := seq_len(.N), by = id]
  lag0 <- d[, .(direct = stats::cov(Y1, Y2, use = "complete.obs")), by = k]
  lag1 <- merge(d[, .(id, k, Y1)], d[, .(id, k = k + 1L, Y2)], by = c("id", "k"))[,
    .(direct = stats::cov(Y1, Y2, use = "complete.obs")), by = k]

  o <- suppressMessages(ctFitCheckCov(ctstantestfit, cor = FALSE, x = "obs",
    minpairn = 2, nsamples = 5, plot = FALSE))
  expect_setequal(unique(o$lag), 0:1)
  expect_identical(attr(o, "ctFitCovCheck_x"), "Observation")
  for (l in 0:1) {
    cmp <- merge(o[lag == l & rowvar == "Y1" & colvar == "Y2", .(k = x, empirical)],
      if (l == 0) lag0 else lag1, by = "k")
    expect_equal(nrow(cmp), if (l == 0) 10L else 9L)
    expect_equal(cmp$empirical, cmp$direct)
  }

  l <- suppressMessages(ctFitCheckCov(ctstantestfit, cor = FALSE, lags = 0:1,
    minpairn = 2, nsamples = 5, plot = FALSE))
  expect_false("x" %in% names(l))
  avg <- o[, .(mean = mean(empirical)), by = .(rowvar, colvar, lag)]
  m <- merge(l[, .(rowvar, colvar, lag, empirical)], avg,
    by = c("rowvar", "colvar", "lag"))
  expect_equal(m$empirical, m$mean)

  gg <- ctFitCovCheckPlot(o, cor = FALSE)
  expect_s3_class(gg$Y1$facet, "FacetGrid")
  expect_identical(gg$Y1$labels$x, "Observation")
})

test_that("ctFitCheckCov(x = <column>) bins a numeric column and keeps a factor's order (stan)", {
  skip_on_cran()

  tb <- suppressMessages(ctFitCheckCov(ctstantestfit, x = "time", xbins = 4,
    nsamples = 5, plot = FALSE))
  expect_lte(length(unique(tb$x)), 4L)
  expect_identical(attr(tb, "ctFitCovCheck_x"), "time (bin mean)")

  d <- data.table::as.data.table(ctsem:::.ctFitLongData(ctstantestfit))
  d[, phase := factor(ifelse(time < 3, "late", "early"), levels = c("late", "early"))]
  ph <- suppressMessages(ctFitCheckCov(ctstantestfit, data = d, x = "phase",
    nsamples = 5, plot = FALSE))
  expect_identical(levels(ph$x), c("late", "early"))
  expect_identical(attr(ph, "ctFitCovCheck_x"), "phase")

  expect_error(suppressMessages(ctFitCheckCov(ctstantestfit, x = "nosuch",
    nsamples = 5, plot = FALSE)), "column of the data")
})

test_that("ctFitCheckCov(lagrows = 'observed') lags a variable by its own observations (stan)", {
  skip_on_cran()

  # Y1 at every other row only, in the data and the generated sets alike.
  fit <- ctstantestfit
  d <- data.table::as.data.table(ctsem:::.ctFitLongData(fit))
  d[, k := seq_len(.N), by = id]
  even <- which(d$k %% 2 == 0)
  d[even, Y1 := NA]
  fit$generated$Y[, even, "Y1"] <- NA

  rowlag <- suppressMessages(ctFitCheckCov(fit, data = d, cor = FALSE,
    x = "obs", minpairn = 2, nsamples = 5, plot = FALSE))
  o <- suppressMessages(ctFitCheckCov(fit, data = d, cor = FALSE, x = "obs",
    minpairn = 2, nsamples = 5, plot = FALSE, lagrows = "observed"))
  expect_true(all(is.na(rowlag[lag == 1 & rowvar == "Y1" & colvar == "Y1", empirical])))

  y1 <- d[!is.na(Y1), .(id, k, Y1)]
  y1[, prev := data.table::shift(Y1), by = id]
  direct <- y1[!is.na(prev), .(direct = stats::cov(Y1, prev), n = .N), by = k]
  cmp <- merge(o[lag == 1 & rowvar == "Y1" & colvar == "Y1",
    .(k = x, empirical, n_empirical, q50)], direct, by = "k")
  expect_equal(nrow(cmp), 4L)
  expect_equal(cmp$empirical, cmp$direct)
  expect_equal(cmp$n_empirical, cmp$n)
  expect_true(all(is.finite(cmp$q50)))
  # Y2 is complete, so its lags are the rows' own.
  expect_equal(o[colvar == "Y2", empirical], rowlag[colvar == "Y2", empirical])
})
