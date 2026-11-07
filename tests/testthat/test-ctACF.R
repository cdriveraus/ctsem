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

  q <- suppressWarnings(suppressMessages(ctsem:::ctACFquantiles(acfsamples(8))))

  expect_true(all(qcolnames() %in% names(q)))
  expect_false(anyNA(q[[qcolnames()[1]]]))
  # the band has to bracket the median, or it is not a band
  expect_true(all(q[[qcolnames()[1]]] <= q[[qcolnames()[3]]]))
})

test_that("a panel too short to tune skips the band rather than erroring", {
  skip_on_cran()
  skip_if_not_installed("qgam")

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

# The lag correlations, computed directly: every ordered pair of rows (a, b) of a
# subject, a's time minus b's rounded to the nearest whole timestep, and the
# correlation of x at a with y at b over the pairs at each lag, pooled over
# subjects. ctACF() once took the lags from a grid of doubles that drifted off the
# observed times, so most pairs fell at the wrong lag: on whole days in weeks, the
# lag of one day rested on 25 of 464 pairs and read 0.44 for white noise.
directlagcor <- function(d, x, y, timestep, nlag) {
  out <- list()
  for (s in unique(d$id)) {
    ds <- d[d$id == s, ]
    ab <- expand.grid(a = seq_len(nrow(ds)), b = seq_len(nrow(ds)))
    u <- ds$time[ab$a] - ds$time[ab$b]
    out[[length(out) + 1]] <- data.frame(lag = sign(u) * floor(abs(u) / timestep + 0.5),
      x = ds[[x]][ab$a], y = ds[[y]][ab$b])
  }
  p <- do.call(rbind, out)
  p <- p[!is.na(p$x) & !is.na(p$y) & abs(p$lag) <= nlag, ]
  if (x == y) p <- p[p$lag >= 1, ]
  res <- do.call(rbind, lapply(split(p, p$lag), function(q)
    data.frame(lag = q$lag[1], Pairs = nrow(q), ACF = stats::cor(q$x, q$y))))
  res[res$Pairs >= 3, ]
}

# Two gappy variables, x autocorrelated and y driven by it, at the given times.
gappypanel <- function(times, seed) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_along(times), function(i) {
    t <- times[[i]]
    x <- numeric(length(t))
    x[1] <- stats::rnorm(1)
    for (r in seq_along(t)[-1]) x[r] <- exp(-(t[r] - t[r - 1]) / 2) * x[r - 1] + stats::rnorm(1, 0, .6)
    data.frame(id = i, time = t, x = x + stats::rnorm(1), y = .6 * x + stats::rnorm(length(t)))
  }))
  d$x[stats::runif(nrow(d)) < .2] <- NA
  d$y[stats::runif(nrow(d)) < .2] <- NA
  d[sample(nrow(d)), ] # unsorted, as data may come
}

expect_matches_direct <- function(d, timestep, time.max) {
  ac <- suppressMessages(ctACF(d, varnames = c("x", "y"), timestep = timestep,
    time.max = time.max, nboot = 0, plot = FALSE))
  if (identical(timestep, "auto")) {
    d <- d[order(d$id, d$time), ]
    timestep <- 0.5 * stats::median(unlist(lapply(split(d$time, d$id), diff)))
  }
  nlag <- ceiling(time.max / timestep)
  for (v in list(c("x", "x"), c("y", "y"), c("x", "y"), c("y", "x"))) {
    name <- if (v[1] == v[2]) v[1] else paste0(v[1], "_", v[2])
    direct <- directlagcor(d, v[1], v[2], timestep, nlag)
    got <- ac[Variable == name]
    expect_equal(got$TimeInterval, direct$lag * timestep, tolerance = 1e-12)
    expect_equal(got$Pairs, direct$Pairs)
    expect_equal(got$ACF, direct$ACF, tolerance = 1e-10)
  }
  invisible(ac)
}

test_that("ctACF gives the pooled correlation of the pairs at each lag, on irregular gappy data", {
  # continuous, irregular times
  set.seed(3)
  times <- lapply(1:25, function(i) cumsum(stats::rexp(sample(8:20, 1), 1)))
  ac <- expect_matches_direct(gappypanel(times, 4), timestep = .5, time.max = 4)
  expect_gt(ac[Variable == "x" & TimeInterval == .5]$ACF, .4) # x is autocorrelated...
  expect_gt(ac[Variable == "x_y" & TimeInterval == 0]$ACF, .3) # ...and y moves with it

  # whole days in weeks, as the book's data: the timestep from the data is a
  # rounding error away from 1/7, and the lag of one day takes every pair a day apart
  set.seed(5)
  times <- lapply(1:25, function(i) sort(sample(0:60, sample(10:25, 1))) / 7)
  ac <- expect_matches_direct(gappypanel(times, 6), timestep = "auto", time.max = 2)
})

test_that("on white noise at whole days in weeks the lag one day correlation is near zero", {
  set.seed(7)
  # intervals of 1-3 days with a median of 2, so the timestep from the data is a day
  d <- do.call(rbind, lapply(1:150, function(i) {
    t <- cumsum(c(0, sample(c(1, 2, 2, 3), sample(5:30, 1), replace = TRUE))) / 7
    data.frame(id = i, time = t, X = stats::rnorm(length(t)))
  }))
  ac <- suppressMessages(ctACF(d, varnames = "X", nboot = 0, plot = FALSE))
  oneday <- ac[abs(TimeInterval - 1/7) < 1e-9]
  dd <- data.table::as.data.table(d)[order(id, time)]
  dd[, `:=`(dt = time - data.table::shift(time), Xl = data.table::shift(X)), by = id]
  expect_equal(oneday$Pairs, dd[abs(dt * 7 - 1) < 1e-9, .N])
  expect_lt(abs(oneday$ACF), 3 / sqrt(oneday$Pairs))
})

test_that("bootstrap samples resample subjects around the data's estimate", {
  set.seed(8)
  d <- data.frame(id = rep(1:60, each = 20), time = rep(1:20, 60), X = stats::rnorm(1200))
  ac <- suppressMessages(ctACF(d, varnames = "X", timestep = 1, time.max = 5, nboot = 200, plot = FALSE))
  expect_setequal(unique(ac$Sample), 0:200)
  boot <- ac[Sample > 0, list(m = mean(ACF), s = stats::sd(ACF)), by = TimeInterval]
  full <- ac[Sample == 0]
  # white noise: the spread over resampled subjects is the sampling error of a
  # correlation from that many pairs, and centres on the data's estimate
  expect_true(all(boot$s / (1 / sqrt(full$Pairs)) > .7 & boot$s / (1 / sqrt(full$Pairs)) < 1.4))
  expect_true(all(abs(boot$m - full$ACF) < 3 * boot$s / sqrt(200)))
})

test_that("center and scale act within subjects, whatever the id column is called", {
  set.seed(9)
  d <- data.frame(Subject = rep(1:30, each = 15), Time = rep(1:15, 30))
  d$X <- rep(stats::rnorm(30, 0, 2), each = 15) + stats::rnorm(450) * rep(stats::runif(30, .5, 2), each = 15)
  manual <- d
  manual$X <- stats::ave(d$X, d$Subject, FUN = function(x) (x - mean(x)) / stats::sd(x))
  args <- list(varnames = "X", idcol = "Subject", timecol = "Time", timestep = 1, time.max = 4,
    nboot = 0, plot = FALSE)
  both <- suppressMessages(do.call(ctACF, c(list(d, center = TRUE, scale = TRUE), args)))
  expect_equal(both, suppressMessages(do.call(ctACF, c(list(manual), args))))
  # person means, removed, take the stable differences out of every lag
  raw <- suppressMessages(do.call(ctACF, c(list(d), args)))
  centred <- suppressMessages(do.call(ctACF, c(list(d, center = TRUE), args)))
  expect_true(all(raw$ACF > .3))
  expect_true(all(abs(centred$ACF) < .2))
})

test_that("ctDiscretiseData puts every time on a grid of equal steps", {
  # whole days in weeks, with a timestep a rounding error away from a day, used
  # to leave observations beside grid rows at nearly the same time
  d <- data.frame(id = rep(1:2, c(4, 3)), time = c(1, 2, 4, 7, 3, 4, 9) / 7, X = 1:7)
  step <- 1/7 + 4e-16 # as half the median interval of the book's data came out
  g <- ctDiscretiseData(d, timestep = step)
  expect_equal(nrow(g), 7 + 7) # days 1-7 and 3-9
  for (i in 1:2) {
    gi <- g[g$id == i, ]
    expect_equal(diff(gi$time), rep(step, nrow(gi) - 1), tolerance = 1e-12)
    expect_equal(gi$time[1], min(d$time[d$id == i])) # each subject's grid starts at its first time
  }
  expect_equal(g$X[!is.na(g$X)], d$X)
  expect_equal(g$time[!is.na(g$X)], d$time, tolerance = 1e-12)

  # values that land on one grid point are averaged, and said so
  d2 <- data.frame(id = 1, time = c(0, .9, 1.1, 2.6), X = c(1, 2, 4, 8))
  expect_warning(g2 <- ctDiscretiseData(d2, timestep = 1), "1 cells of data removed")
  expect_equal(g2$time, 0:3)
  expect_equal(g2$X, c(1, 3, NA, 8))

  # a variable absent at an occupied time is missing, not the mean of nothing
  d3 <- data.frame(id = 1, time = 0:1, X = c(1, NA), Y = c(1, 2))
  g3 <- ctDiscretiseData(d3, timestep = 1)
  expect_identical(g3$X, c(1, NA))
})
