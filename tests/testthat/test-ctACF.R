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
