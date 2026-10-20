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
