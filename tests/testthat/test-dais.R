# Doubly adaptive importance sampling as a stage of 'checked' (R/ctDais.R):
# the normal moved to the posterior's moments before anything is sampled.

# A skewed target with exact moments: the first k coordinates of a rotated
# frame are log-gamma(a) (density exp(a y - exp(y)); mean digamma(a),
# variance trigamma(a)), the rest standard normal. Its Laplace normal at the
# mode has the mean of a log-gamma(0.5) coordinate 1.3 of its sds off and
# its sd 0.64 of the truth.
.skewTarget <- function(d = 20L, a = c(0.5, 0.5, 1, 2), seed = 11L) {
  set.seed(seed)
  k <- length(a)
  Q <- qr.Q(qr(matrix(stats::rnorm(d * d), d)))
  lpg <- function(x) {
    y <- as.numeric(crossprod(Q, x))
    e <- exp(y[seq_len(k)])
    out <- sum(a * y[seq_len(k)] - e) - 0.5 * sum(y[-seq_len(k)]^2)
    attr(out, 'gradient') <- as.numeric(Q %*% c(a - e, -y[-seq_len(k)]))
    out
  }
  list(lpg = ctsem:::.ctImisPointGradbatch(lpg),
    mode = as.numeric(Q %*% c(log(a), rep(0, d - k))),
    lapcov = Q %*% diag(c(1 / a, rep(1, d - k))) %*% t(Q),
    mean = as.numeric(Q %*% c(digamma(a), rep(0, d - k))),
    sd = sqrt(diag(Q %*% diag(c(trigamma(a), rep(1, d - k))) %*% t(Q))))
}

test_that("DAIS moves a normal from the mode to a skewed posterior's moments", {
  t <- .skewTarget()
  # The Laplace normal's mean is off by over half a posterior sd somewhere
  # (1.3 of its own sds along a log-gamma(0.5) coordinate, diluted over the
  # rotated frame's parameters).
  expect_gt(max(abs(t$mode - t$mean) / t$sd), 0.5)
  set.seed(12)
  r <- ctsem:::.ctDaisRun(t$lpg, t$mode, t$lapcov, list(daisDraws = 64L))
  expect_true(r$used)
  expect_true(r$settled)
  expect_gt(r$averaged, 0L)
  expect_lt(max(abs(r$mean - t$mean) / t$sd), 0.15)
  ratio <- sqrt(diag(r$cov)) / t$sd
  expect_true(all(ratio > 0.85 & ratio < 1.1))
  # Even draws an iteration, counted; the trace carries one row per iteration.
  expect_equal(r$draws, 64L)
  expect_equal(r$gradients, 64L * r$iterations)
  expect_equal(nrow(r$trace), r$iterations)
  # A Gaussian target from its own normal: nothing to move, and it settles
  # within the window's two lengths plus the averaging.
  P <- solve(t$lapcov)
  g <- ctsem:::.ctImisPointGradbatch(function(x) {
    out <- -0.5 * sum(x * (P %*% x))
    attr(out, 'gradient') <- -as.numeric(P %*% x)
    out
  })
  set.seed(13)
  r0 <- ctsem:::.ctDaisRun(g, rep(0, 20), t$lapcov, list(daisDraws = 64L))
  expect_true(r0$settled)
  expect_lt(max(abs(r0$mean) / sqrt(diag(t$lapcov))), 0.1)
  expect_equal(r0$cov, t$lapcov, tolerance = 0.15)
})

test_that("DAIS runs inside 'checked' where the check names a direction, and not where it does not", {
  t <- .skewTarget(d = 6L, a = c(0.5, 1))
  set.seed(14)
  out <- suppressMessages(ctsem:::ctOptimComputeUncertainty(t$mode, list(nsubjects = 1),
    sm = NULL, lpgFunc = t$lpg, uncertainty = 'checked', finishsamples = 20, verbose = 0))
  expect_gt(ncol(out$details$normalcheck$directions), 0L)
  set.seed(15)
  drawn <- suppressWarnings(suppressMessages(ctsem:::.ctOptimDrawSamples(out, 'imis',
    list(), t$mode, 1000, lpg = t$lpg)))
  dais <- drawn$uncertaintyfit$details$dais
  expect_true(isTRUE(dais$used))
  expect_gt(dais$iterations, 0L)
  # The recheck's record, and the normal the stage leaves: its covariance is
  # the fit's proposal, and the draws' means sit on the posterior's.
  expect_true(is.numeric(dais$check$directions))
  expect_equal(dim(dais$cov), c(6L, 6L))
  expect_lt(max(abs(colMeans(drawn$samples) - t$mean) / t$sd), 0.15)
  expect_equal(sqrt(diag(drawn$uncertaintyfit$cov)) / t$sd, rep(1, 6), tolerance = 0.2)
  # Off by name, the Hessian's normal is where sampling starts, as before.
  set.seed(15)
  plain <- suppressWarnings(suppressMessages(ctsem:::.ctOptimDrawSamples(out, 'imis',
    list(dais = FALSE), t$mode, 1000, lpg = t$lpg)))
  expect_null(plain$uncertaintyfit$details$dais)
  # And a Gaussian posterior names nothing, so nothing runs.
  P <- solve(t$lapcov)
  g <- ctsem:::.ctImisPointGradbatch(function(x) {
    out <- -0.5 * sum(x * (P %*% x))
    attr(out, 'gradient') <- -as.numeric(P %*% x)
    out
  })
  set.seed(16)
  outg <- suppressMessages(ctsem:::ctOptimComputeUncertainty(rep(0, 6), list(nsubjects = 1),
    sm = NULL, lpgFunc = g, uncertainty = 'checked', finishsamples = 20, verbose = 0))
  expect_equal(ncol(outg$details$normalcheck$directions), 0L)
  drawng <- ctsem:::.ctOptimDrawSamples(outg, 'imis', list(), rep(0, 6), 100, lpg = g)
  expect_null(drawng$uncertaintyfit$details$dais)
  expect_identical(drawng$draws, 'normal')
})

test_that("the sampler's metric is the DAIS normal's precision where there is one", {
  set.seed(17)
  C <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  fit <- list(uncertainty = list(hessian = -diag(5), details = list(dais = list(cov = C),
    normalcheck = list(directions = matrix(1, 5, 1), widths = 4))))
  expect_equal(-ctsem:::.ctBackendSampleMetric(fit), solve(ctsem:::ctOptimSafeCov(C)),
    tolerance = 1e-10)
  fit$uncertainty$details$dais <- NULL
  # Without one, the Hessian rescaled to the check's widths, as before.
  expect_false(isTRUE(all.equal(-ctsem:::.ctBackendSampleMetric(fit), solve(C))))
})
