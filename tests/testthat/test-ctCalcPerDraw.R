# ctCalcPerDraw() and the medians the summaries take from the estimate. All on
# ctsemtestfit, an optimised fit with normal draws, which needs no Julia.

test_that("a median is the estimate's value only where one coordinate moves it", {
  set.seed(3)
  raw <- cbind(stats::rnorm(400, 1, .5), stats::rnorm(400, -2, .3))
  at <- c(exp(1), 1^2 + (-2)^2, -2)
  values <- cbind(exp(raw[, 1]), raw[, 1]^2 + raw[, 2]^2, raw[, 2])
  known <- ctsem:::.ctExactMedians(values, raw, at)
  expect_equal(known, c(exp(1), NA, -2))
  # Without draws symmetric about the estimate, nothing is known.
  expect_true(all(is.na(ctsem:::.ctExactMedians(values, NULL, at))))
})

test_that("ctCalcPerDraw computes a quantity draw by draw", {
  fit <- ctsemtestfit
  e <- ctExtract(fit)
  one <- ctCalcPerDraw(fit, function(m) m$DRIFT[1, 1])
  expect_s3_class(one, "ctCalcPerDraw")
  expect_equal(as.numeric(one), as.numeric(e$pop_DRIFT[, 1, 1]))

  both <- ctCalcPerDraw(fit, function(m) c(drift = m$DRIFT[1, 1],
    asym = m$asymDIFFUSIONcov[1, 1]))
  expect_equal(dim(both), c(dim(e$pop_DRIFT)[1L], 2L))
  s <- summary(both)
  expect_identical(rownames(s), c("drift", "asym"))
  estimate <- suppressMessages(ctSummaryMatrices(fit))$estimate
  expect_equal(unname(attr(both, "estimate")),
    c(estimate$DRIFT[1, 1], estimate$asymDIFFUSIONcov[1, 1]))
  # The drift's median is its estimate; the asymptotic variance, which bends
  # with every parameter of the drift and diffusion, keeps its draws' median.
  expect_equal(s["drift", "50%"], estimate$DRIFT[1, 1])
  expect_equal(s["asym", "50%"], stats::median(e$pop_asymDIFFUSIONcov[, 1, 1],
    na.rm = TRUE))

  # A matrix comes back by row and column; the Lyapunov equation holds in
  # every draw, where a matrix of medians does not.
  lyapunov <- ctCalcPerDraw(fit, function(m) {
    A <- m$DRIFT
    S <- m$asymDIFFUSIONcov
    A %*% S + S %*% t(A) + m$DIFFUSIONcov
  })
  expect_equal(dim(lyapunov)[-1L], c(2L, 2L))
  expect_lt(max(abs(lyapunov), na.rm = TRUE), 1e-6)
  expect_identical(colnames(summary(lyapunov))[1:2], c("row", "col"))

  few <- ctCalcPerDraw(fit, function(m) m$DRIFT[1, 1], nsamples = 10)
  expect_equal(dim(few), c(10L, 1L))
})

test_that("summary's 50% holds the estimate where it is the median", {
  fit <- ctsemtestfit
  s <- summary(fit, digits = 8)
  estimate <- suppressMessages(ctSummaryMatrices(fit))
  expect_false("mode" %in% colnames(s$parmatrices))
  pm <- s$parmatrices
  drift <- pm[pm$matrix == "DRIFT", ]
  expect_equal(drift[["50%"]], estimate$estimate$DRIFT[cbind(drift$row, drift$col)],
    tolerance = 1e-7)
  # ctSummaryMatrices()'s default, the median, gives the same.
  expect_equal(estimate$DRIFT, estimate$estimate$DRIFT)
})
