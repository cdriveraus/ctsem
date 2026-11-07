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

test_that("every lag of a draw comes from one decomposition, and agrees with expm", {
  # .ctExpmTimes takes expm(D t) at all t from D's eigendecomposition, and
  # falls back to expm() at each t where the eigenvectors are ill-conditioned.
  times <- c(0, .1, 1, 3.7, 20)
  each <- function(D, discrete = FALSE) aperm(sapply(times, function(t)
    if (discrete) { m <- diag(nrow(D)); for (n in seq_len(t)) m <- m %*% D; m }
    else expm::expm(D * t), simplify = 'array'), c(3, 1, 2))
  drifts <- list(
    real = matrix(c(-.4, .2, .1, -.7), 2),
    complex = matrix(c(-.1, -.8, .8, -.1), 2),
    defective = matrix(c(-.5, 1, 0, -.5), 2),
    neardefective = matrix(c(-.5, 1, 0, -.5000001), 2),
    zero = matrix(0, 2, 2))
  for (nm in names(drifts)) expect_equal(ctsem:::.ctExpmTimes(drifts[[nm]], times),
    each(drifts[[nm]]), tolerance = 1e-10, label = nm)
  D <- diag(2) + .5 * drifts$complex
  expect_equal(ctsem:::.ctExpmTimes(D, 0:6, discrete = TRUE),
    { times <- 0:6; each(D, discrete = TRUE) }, tolerance = 1e-10)
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
})
