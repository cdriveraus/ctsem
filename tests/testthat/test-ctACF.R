# The quantile band behind plotctACF().
#
# ctACFquantiles() tunes a qgam learning rate on the first variable that will
# take one and reuses it for the rest. Short panels have too few distinct time
# intervals for the spline basis and the tuning fails, so the three cases that
# matter are: tuning succeeds for every variable, for some, and for none. The
# last is the one with no natural exercise -- it happens on the residual panels
# of a small fit, through ctReport() -- so it is built here explicitly.

acfsamples <- function(ntime, nid = 10, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(nid), each = ntime),
    time = rep(seq_len(ntime), times = nid),
    Y1 = stats::rnorm(nid * ntime))
  suppressWarnings(suppressMessages(
    ctACF(d, varnames = "Y1", idcol = "id", timecol = "time",
      timestep = 1, time.max = ntime - 1, nboot = 5, plot = FALSE)))
}

qcolnames <- function(quantiles = c(.025, .5, .975)) paste0("Q", quantiles * 100, "%")

test_that("quantile bands are estimated when every variable can be tuned", {
  skip_on_cran()
  skip_if_not_installed("qgam")
  skip_if_not_installed("collapse")

  q <- suppressWarnings(suppressMessages(ctsem:::ctACFquantiles(acfsamples(8))))

  expect_true(all(qcolnames() %in% names(q)))
  expect_false(anyNA(q[[qcolnames()[1]]]))
  # the band has to bracket the median, or it is not a band
  expect_true(all(q[[qcolnames()[1]]] <= q[[qcolnames()[3]]]))
})

test_that("a panel too short to tune skips the band rather than erroring", {
  skip_on_cran()
  skip_if_not_installed("qgam")
  skip_if_not_installed("collapse")

  ac <- acfsamples(4)
  expect_lt(length(unique(ac$TimeInterval)), 4) # short enough that tuning fails

  # Before the fallback was fixed this errored with "$ operator is invalid for
  # atomic vectors", because the un-tuned learning rate stayed an atomic NA.
  leaked <- utils::capture.output(
    q <- suppressWarnings(suppressMessages(ctsem:::ctACFquantiles(ac))),
    type = "message")

  expect_s3_class(q, "data.table")
  expect_false(any(qcolnames() %in% names(q)))
  expect_equal(nrow(q), nrow(ac))
  # a handled try() must not print its error: this is what users saw as a bare
  # "Error in place.knots(x, nk)" on stderr from ctReport()'s residual panels.
  expect_equal(leaked, character(0))
  expect_message(suppressWarnings(ctsem:::ctACFquantiles(ac)),
    "could not be fitted")

  # and the plot falls back to the raw samples instead of failing
  leakedplot <- utils::capture.output(
    gg <- suppressWarnings(suppressMessages(plotctACF(ac))), type = "message")
  expect_s3_class(gg, "ggplot")
  expect_equal(leakedplot, character(0))
  expect_message(suppressWarnings(plotctACF(ac)), "plotting ACF samples")
})

test_that("a variable that cannot be tuned leaves its own band NA and keeps the rest", {
  skip_on_cran()
  skip_if_not_installed("qgam")
  skip_if_not_installed("collapse")

  long <- acfsamples(8)
  # 'short' carries only two distinct time intervals, so its spline cannot be
  # fitted; 'long' can. Put the failing one first, so the shared learning rate
  # also has to survive a failed first attempt.
  ac <- rbind(
    data.table::copy(long)[TimeInterval %in% c(1, 2), ][, Variable := "short"],
    data.table::copy(long)[, Variable := "long"])

  leaked <- utils::capture.output(
    q <- suppressWarnings(suppressMessages(ctsem:::ctACFquantiles(ac))),
    type = "message")

  expect_true(all(qcolnames() %in% names(q)))
  expect_false(anyNA(q[Variable == "long"][[qcolnames()[1]]]))
  expect_true(all(is.na(q[Variable == "short"][[qcolnames()[1]]])))
  expect_equal(leaked, character(0))
  expect_message(suppressWarnings(ctsem:::ctACFquantiles(ac)), "short")
})

test_that("cross correlations keep lag 0, autocorrelations drop it", {
  skip_if_not_installed("collapse")

  # Y2 shares Y1's value at the same occasion only, so the cross correlation is
  # large at lag 0 and near zero elsewhere. Lag 0 used to be dropped for cross
  # correlations as well as autocorrelations, hiding the largest value.
  set.seed(1)
  d <- data.frame(id = rep(1:20, each = 20), time = rep(1:20, times = 20),
    Y1 = stats::rnorm(400))
  d$Y2 <- d$Y1 + stats::rnorm(400)
  ac <- suppressMessages(ctACF(d, varnames = c("Y1", "Y2"), idcol = "id",
    timecol = "time", timestep = 1, time.max = 3, nboot = 0, plot = FALSE))

  expect_false(0 %in% ac[Variable == "Y1"]$TimeInterval)
  cc <- ac[Variable == "Y1_Y2"]
  expect_true(0 %in% cc$TimeInterval)
  expect_gt(cc[TimeInterval == 0]$ACF, 0.5)
  expect_true(all(abs(cc[TimeInterval != 0]$ACF) < 0.2))
})

# A ctACF-shaped table with a known truth: Sample 0 is the full data estimate,
# Samples 1..nboot scatter around the truth with standard deviation sd.
acftable <- function(truth, lags, sd, nboot = 100, seed = 1, variable = "Y1_Y2") {
  set.seed(seed)
  data.table::rbindlist(lapply(0:nboot, function(s) data.table::data.table(
    Sample = s, TimeInterval = lags, Variable = variable,
    ACF = truth + stats::rnorm(length(lags), 0, sd))))
}

test_that("the weighted spline follows a precisely estimated peak", {
  skip_if_not_installed("mgcv")

  # The quantile spline flattened a lag 0 cross correlation of about .4 to about
  # .05: a peak one lag wide, smoothed with the same stiffness as everything
  # else. Weighted by its precision, the peak is followed.
  lags <- -10:10
  truth <- ifelse(lags == 0, .4, 0)
  q <- suppressMessages(ctsem:::ctACFweightedSpline(acftable(truth, lags, sd = .03)))

  expect_true(all(qcolnames() %in% names(q)))
  med <- unique(q[, c("TimeInterval", "Q50%")])
  # thresholds hold across seeds (lowest peak about .25 over 30), while the
  # quantile spline gives about .05 here
  expect_gt(med[TimeInterval == 0][["Q50%"]], 0.2)
  expect_true(all(abs(med[abs(TimeInterval) >= 3][["Q50%"]]) < 0.15))
  expect_true(all(q[["Q2.5%"]] <= q[["Q50%"]] & q[["Q50%"]] <= q[["Q97.5%"]]))
})

test_that("the weighted spline smooths away noise consistent with its standard errors", {
  skip_if_not_installed("mgcv")

  # No signal, and one noisy lag far from zero: the spline should stay flat
  # rather than chase it, since the departure is within that lag's error.
  lags <- 1:12
  ac <- acftable(rep(0, 12), lags, sd = .05, seed = 2, variable = "Y1")
  sd7 <- .2 # one imprecise lag, as when few pairs fall at that interval
  ac[TimeInterval == 7, ACF := stats::rnorm(.N, 0, sd7)]
  ac[Sample == 0 & TimeInterval == 7, ACF := .3]
  q <- suppressMessages(ctsem:::ctACFweightedSpline(ac))

  med <- unique(q[, c("TimeInterval", "Q50%")])
  expect_true(all(abs(med[["Q50%"]]) < 0.15)) #half the stray estimate; about .1 at most over 30 seeds
})

test_that("without bootstrap samples the weighted spline falls back to the samples", {
  skip_if_not_installed("mgcv")

  ac <- acftable(rep(0, 10), 1:10, sd = .05, nboot = 0, variable = "Y1")
  expect_message(q <- ctsem:::ctACFweightedSpline(ac), "could not be fitted")
  expect_false(any(qcolnames() %in% names(q)))
  expect_equal(nrow(q), nrow(ac))
  expect_message(gg <- plotctACF(ac), "plotting ACF samples")
  expect_s3_class(gg, "ggplot")
})

test_that("plotctACF draws the weighted spline by default, the quantile spline on request", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  skip_if_not_installed("qgam")

  ac <- acftable(exp(-(1:12) / 4), 1:12, sd = .05, nboot = 20, variable = "Y1")
  gg <- suppressWarnings(suppressMessages(plotctACF(ac, reducedXlim = 0)))
  expect_s3_class(gg, "ggplot")
  # the full data estimates are drawn as points over the band
  expect_true(any(vapply(gg$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))))
  wq <- suppressMessages(ctsem:::ctACFweightedSpline(ac))
  built <- ggplot2::ggplot_build(gg)$data[[2]] # the spline line
  expect_equal(sort(unique(built$y)), sort(unique(wq[["Q50%"]])), tolerance = 1e-8)

  gq <- suppressWarnings(suppressMessages(plotctACF(ac, reducedXlim = 0, method = "quantile")))
  expect_s3_class(gq, "ggplot")
  expect_error(plotctACF(ac, method = "nonsense"))
})
