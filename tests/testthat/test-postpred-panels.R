# The two structure panels added in 3.12.0 when the ctFitCheck() dashboard was
# removed and what was worth keeping in it moved into ctPostPredPlots().
#
# MeanTrajectory restores a panel lost when the two posterior predictive
# functions were merged in 201e8305. Its predecessor grouped by exact Time,
# which gives one observation per group -- and so a zero-width band and a
# "mean" that is the raw data -- for anything but a balanced panel, so the
# restored version bins time.
#
# LaggedCovariance is ctFitCheckCov() rendered as dashboard panels. The point
# of asserting on it here is that the panel and a direct ctFitCheckCov() call
# give the SAME numbers, rather than the panel quietly recomputing something
# adjacent.

test_that("both new panels are in the structure family and in 'all'", {
  expect_true(all(c("MeanTrajectory", "LaggedCovariance") %in%
      ctsem:::.ctPostPredGroups$structure))
  expect_true(all(c("MeanTrajectory", "LaggedCovariance") %in%
      ctsem:::.ctPostPredResolvePanels("all")))
  # And they are addressable individually, which is how a caller skips the
  # lagged covariance computation.
  expect_identical(ctsem:::.ctPostPredResolvePanels("MeanTrajectory"),
    "MeanTrajectory")
})

test_that("MeanTrajectory bands the mean and the data separately, on observed rows only (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  set.seed(3)
  f <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 20, cores = 1))

  p <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "MeanTrajectory", timebins = 10)))
  expect_named(p, "MeanTrajectory")
  d <- data.table::as.data.table(p$MeanTrajectory$data)
  expect_true(all(c("variable", "Bin", "TimeMid", "n", "Observed",
    "MeanLo", "MeanMid", "MeanHi", "ObsLo", "ObsHi") %in% names(d)))

  # Ten bins per variable, or fewer where time has too few distinct values to
  # cut that finely.
  expect_true(all(d[, .N, by = variable]$N <= 10))
  expect_true(all(d[, .N, by = variable]$N >= 2))

  # The two bands are different quantities and the panel is pointless if they
  # are not: the sampling distribution of a bin mean is narrower than the
  # spread of the observations it averages.
  expect_true(all(d$MeanHi - d$MeanLo > 0))
  expect_true(all(d$MeanHi - d$MeanLo <= d$ObsHi - d$ObsLo))
  expect_lt(median(d$MeanHi - d$MeanLo), median(d$ObsHi - d$ObsLo) / 2)

  # Both sides average exactly the rows that were observed. Generated values
  # exist for every row including the missing ones, so a model mean taken over
  # more rows than the observed mean would not be a comparison -- and the
  # per-observation dedup matters too, since the long table holds each observed
  # value once per draw.
  dd <- suppressMessages(ctPostPredData(f))
  nobs <- nrow(unique(dd[is.finite(dd$obsValue), c("variable", "row")]))
  expect_identical(sum(d$n), nobs)

  # The time axis is a property of the observation schedule, not of how many
  # draws the fit carries. dat holds each row once per draw, so binning it
  # directly put the bin edges -- and with them the plotted x positions and the
  # bin counts -- at the mercy of nsamples.
  p6 <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "MeanTrajectory", timebins = 10, nsamples = 6)))
  d6 <- data.table::as.data.table(p6$MeanTrajectory$data)
  data.table::setorder(d6, variable, Bin)
  dsorted <- data.table::copy(d); data.table::setorder(dsorted, variable, Bin)
  expect_identical(d6$n, dsorted$n)
  expect_equal(d6$TimeMid, dsorted$TimeMid)
  expect_equal(d6$Observed, dsorted$Observed)

  # The band is on the right scale: for a model fitted to this data most bin
  # means should land inside their own 95% interval. Loose, because it is a
  # sanity check on the reference rather than a calibration test -- that is
  # what the IntervalCoverage and PIT panels are for.
  inside <- d$Observed >= d$MeanLo & d$Observed <= d$MeanHi
  expect_gt(mean(inside), 0.6)
})

test_that("MeanTrajectory's observed means are the data's own, per time bin (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  set.seed(4)
  f <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 10, cores = 1))
  p <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "MeanTrajectory", timebins = 5)))
  d <- data.table::as.data.table(p$MeanTrajectory$data)

  # Recompute the observed side independently, from the fit's own long data
  # rather than from ctPostPredData(), binning the same way. This is the check
  # that would catch the panel averaging generated values into the observed
  # line, or binning one side on a different axis from the other.
  long <- data.table::as.data.table(ctsem:::.ctFitLongData(ctstantestfit))
  tname <- ctstantestfit$ctstanmodelbase$timeName
  for (v in c("Y1", "Y2")) {
    x <- long[is.finite(get(v)), .(t = get(tname), y = get(v))]
    x[, Bin := ctsem:::.ctPostPredBin(t, 5)]
    ref <- x[, .(Observed = mean(y), TimeMid = mean(t), n = .N), by = Bin]
    got <- d[variable == v]
    data.table::setorder(ref, Bin)
    data.table::setorder(got, Bin)
    expect_equal(got$Observed, ref$Observed)
    expect_equal(got$TimeMid, ref$TimeMid)
    expect_identical(got$n, ref$n)
  }
})

test_that("LaggedCovariance panels are ctFitCheckCov's own numbers (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  set.seed(5)
  # Generate once, so the panel and the direct call read the same draws and
  # the comparison is exact rather than approximate.
  f <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 8, cores = 1))

  p <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "LaggedCovariance", lags = 0:3, lagcor = TRUE)))
  expect_setequal(names(p), c("LaggedCovariance_Y1", "LaggedCovariance_Y2"))
  expect_true(all(vapply(p, inherits, logical(1), "ggplot")))

  direct <- suppressWarnings(suppressMessages(ctFitCheckCov(f, cor = TRUE,
    plot = FALSE, lags = 0:3, variables = c("Y1", "Y2"))))
  panel <- data.table::as.data.table(p$LaggedCovariance_Y1$data)
  ref <- data.table::as.data.table(direct)[rowvar == "Y1" & lag <= 3]
  data.table::setorder(panel, colvar, lag)
  data.table::setorder(ref, colvar, lag)
  expect_equal(panel$empirical, ref$empirical)
  expect_equal(panel$q025, ref$q025)
  expect_equal(panel$q975, ref$q975)

  # lagcor=FALSE is covariances, and the axis label says so -- an argument
  # accepted and then ignored is the trap this package has had before.
  pc <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "LaggedCovariance", lags = 0:1, lagcor = FALSE)))
  expect_match(pc$LaggedCovariance_Y1$labels$y, "^Covariance with")
  expect_match(p$LaggedCovariance_Y1$labels$y, "^Correlation with")

  # Same colour convention as every other panel in this dashboard -- model
  # blue, observed black. ctFitCovCheckPlot() drew the empirical value in red
  # and the model in black, which is the reverse of it, and the reverse of the
  # panels this one now sits beside.
  built <- ggplot2::ggplot_build(p$LaggedCovariance_Y1)
  drawn <- unique(unlist(lapply(built$data, function(z) as.character(z$colour))))
  expect_true(all(c("#111111", "#2166AC") %in% drawn))

  # A correlation with itself at lag 0 is 1 by construction; a covariance is
  # not. Cheap confirmation that the two calls really differ in the numbers
  # and not only in the label.
  pcd <- data.table::as.data.table(pc$LaggedCovariance_Y1$data)
  expect_equal(panel[colvar == "Y1" & lag == 0]$empirical, 1)
  expect_false(isTRUE(all.equal(pcd[colvar == "Y1" & lag == 0]$empirical, 1)))
})

test_that("the default panel set includes both, and each is skippable (stan)", {
  skip_on_cran()
  withr::local_pdf(NULL)
  set.seed(6)
  f <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 6, cores = 1))

  all_p <- suppressWarnings(suppressMessages(ctPostPredPlots(f, lags = 0:2)))
  expect_true("MeanTrajectory" %in% names(all_p))
  expect_true(any(grepl("^LaggedCovariance_", names(all_p))))
  expect_true(all(vapply(all_p, inherits, logical(1), "ggplot")))

  # LogLik rides along in MeanTrajectory as it does in Density, since a mean
  # row log-likelihood over time is a quantity like any other.
  expect_true("LogLik" %in% as.character(all_p$MeanTrajectory$data$variable))

  cal <- suppressWarnings(suppressMessages(
    ctPostPredPlots(f, panels = "calibration")))
  expect_false(any(grepl("MeanTrajectory|LaggedCovariance", names(cal))))
})
