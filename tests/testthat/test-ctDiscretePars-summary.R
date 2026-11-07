# summary() of ctDiscretePars() output against curves with closed forms: an
# autoregression exp(-a t) halves at log(2) / a, and a cross effect t exp(-t)
# peaks at t = 1 with value exp(-1), halves after it where t exp(-t) =
# exp(-1) / 2, and has area 1 - (T + 1) exp(-T) up to T.

test_that("each curve's peak, its lag, its half-time and its area", {
  times <- seq(0, 12, by = 0.001)
  a <- c(0.5, 0.4)
  x <- array(NA_real_, c(2L, 1L, length(times), 2L, 2L), dimnames = list(
    Sample = 1:2, Subject = 'popmean', `Time interval` = times,
    row = c('eta1', 'eta2'), col = c('eta1', 'eta2')))
  for (s in 1:2) {
    x[s, 1, , 1, 1] <- exp(-a[s] * times)
    x[s, 1, , 2, 1] <- -times * exp(-times)
    x[s, 1, , 1, 2] <- 0
    x[s, 1, , 2, 2] <- exp(-0.01 * times)
  }
  class(x) <- c('ctDiscretePars', 'array')
  out <- summary(x)
  get <- function(r, c, q) out[out$row == r & out$col == c & out$quantity == q, ]
  ar <- get('eta1', 'eta1', 'halftime')
  lo <- log(2) / a[1]
  hi <- log(2) / a[2]
  expect_equal(c(ar$`2.5%`, ar$`97.5%`), c(lo + .025 * (hi - lo),
    hi - .025 * (hi - lo)), tolerance = 1e-5)
  expect_equal(ar$mean, mean(log(2) / a), tolerance = 1e-5)
  expect_equal(ar$beyond, 0)
  expect_equal(get('eta1', 'eta1', 'peaktime')$mean, 0)
  cross <- function(q) get('eta2', 'eta1', q)$mean
  expect_equal(cross('peak'), -exp(-1), tolerance = 1e-6)
  expect_equal(cross('peaktime'), 1, tolerance = 1e-6)
  late <- stats::uniroot(function(t) t * exp(-t) - exp(-1) / 2, c(1, 10))$root
  expect_equal(cross('halftime'), late, tolerance = 1e-5)
  expect_equal(cross('area'), -(1 - 13 * exp(-12)), tolerance = 1e-6)
  # A curve that does not halve within the lags computed says so.
  slow <- get('eta2', 'eta2', 'halftime')
  expect_equal(slow$beyond, 1)
  expect_true(is.na(slow$mean))
  # A flat curve has no peak to halve.
  expect_true(is.na(get('eta1', 'eta2', 'halftime')$mean))
})

test_that("summary() reads a fit's curves, and the array still prints and plots", {
  skip_on_cran()
  dp <- suppressMessages(ctDiscretePars(ctstantestfit, times = seq(0, 5, .1),
    nsamples = 20))
  expect_s3_class(dp, 'ctDiscretePars')
  expect_equal(length(dim(dp)), 5L)
  out <- summary(dp)
  expect_equal(nrow(out), 4L * 4L)
  expect_setequal(out$quantity, c('peak', 'peaktime', 'halftime', 'area'))
  expect_output(print(dp))
  expect_s3_class(suppressMessages(ctDiscreteParsPlot(dp)), 'ggplot')
})

test_that("the regressions of every draw at once agree with expm at each draw and lag", {
  # ctDiscreteParsDrift works on all draws together, each lag the last times
  # the step's expm; the reference is expm() of each draw at each lag. The
  # drifts include complex, defective, near-defective, zero and explosive ones,
  # and the lags are unsorted, uneven and repeated.
  drifts <- list(
    matrix(c(-.4, .2, .1, -.7), 2), matrix(c(-.1, -.8, .8, -.1), 2),
    matrix(c(-.5, 1, 0, -.5), 2), matrix(c(-.5, 1, 0, -.5000001), 2),
    matrix(0, 2, 2), matrix(c(.3, .1, 0, -2), 2))
  D <- array(NA_real_, c(length(drifts), 2, 2, 2))
  for (s in seq_along(drifts)) for (j in 1:2) D[s, j, , ] <- drifts[[s]] * j
  times <- c(3.7, 0, .1, 20, 1, .1)
  out <- ctsem:::ctDiscreteParsDrift(list(DRIFT = D), times, 'unit', FALSE, quiet = TRUE)
  for (s in seq_along(drifts)) for (j in 1:2) for (k in seq_along(times))
    expect_equal(out[s, j, k, , ], expm::expm(D[s, j, , ] * times[k]),
      tolerance = 1e-10, label = paste('draw', s, 'subject', j, 'lag', times[k]))
  # Discrete time: powers of the drift matrix.
  P <- D
  for (s in seq_along(drifts)) for (j in 1:2) P[s, j, , ] <- diag(2) + .4 * D[s, j, , ]
  out <- ctsem:::ctDiscreteParsDrift(list(DRIFT = P), c(0, 3, 1, 6), 'unit', FALSE,
    discreteInput = TRUE, quiet = TRUE)
  for (s in seq_along(drifts)) for (k in seq_along(c(0, 3, 1, 6))) {
    m <- diag(2)
    for (n in seq_len(c(0, 3, 1, 6)[k])) m <- m %*% P[s, 1, , ]
    expect_equal(out[s, 1, k, , ], m, tolerance = 1e-10)
  }
})

test_that("the plot's data are the quantiles of the draws at each lag", {
  set.seed(2)
  times <- seq(0, 4, .5)
  x <- array(stats::rnorm(200 * 2 * length(times) * 4), c(200, 2, length(times), 2, 2),
    dimnames = list(Sample = 1:200, Subject = c('p2', 'p1'), `Time interval` = times,
      row = c('X', 'Y'), col = c('X', 'Y')))
  x[1, 1, 2, 2, 1] <- NaN
  class(x) <- c('ctDiscretePars', 'array')
  attributes(x)$impulseType <- 'unit'
  attributes(x)$cov <- FALSE
  g <- ctDiscreteParsPlot(x, indices = rbind(c(2, 1)), quantiles = c(.1, .5, .9))
  d <- as.data.frame(g$data)
  expect_equal(nrow(d), 2 * length(times))
  for (s in c('p1', 'p2')) for (k in seq_along(times)) {
    z <- x[, s, k, 2, 1]
    z <- z[is.finite(z)]
    row <- d[d$Subject == s & d$`Time interval` == times[k], ]
    expect_equal(c(row$ymin, row$value, row$ymax),
      unname(stats::quantile(z, c(.1, .5, .9))))
  }
  # The interval first, the median line second, as code that adjusts the
  # layers by position expects.
  expect_s3_class(g$layers[[1]]$geom, 'GeomRibbon')
  expect_s3_class(g$layers[[2]]$geom, 'GeomLine')
  pooled <- as.data.frame(ctDiscreteParsPlot(x, indices = rbind(c(2, 1)),
    splitSubjects = FALSE)$data)
  z <- c(x[, , 3, 2, 1])
  expect_equal(pooled$value[pooled$`Time interval` == times[3]], stats::median(z))
  # Two limits zoom the plot, and the band keeps every interval: `ylim()`
  # dropped those reaching past the limits, and the scalar test of `ylim`
  # stopped every call that gave two.
  zoomed <- ctDiscreteParsPlot(x, indices = rbind(c(2, 1)), ylim = c(-.5, .5))
  expect_equal(zoomed$coordinates$limits$y, c(-.5, .5))
  built <- ggplot2::ggplot_build(zoomed)
  expect_equal(nrow(built$data[[1]]), nrow(d))
  expect_true(all(is.finite(built$data[[1]]$ymin)))
  expect_s3_class(ctDiscreteParsPlot(x, indices = rbind(c(2, 1)),
    facets = c('Subject', 'Effect')), 'ggplot')
})
