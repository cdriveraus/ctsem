# ctChisqTest() on stan and julia fits, and across backends.
#
# ctChisqTest() used to read fit$stanfit$optimfit$value/fit$stanfit$rawest
# directly, so it errored (opaquely -- "$ operator is invalid for atomic
# vectors") on a julia fit, which keeps the equivalent quantities under
# fit$estimate$logposterior (falling back to $loglik, the same fallback
# .ctBackendSummary() uses) and fit$estimate$raw. It now reads both through
# .ctFitOptimValue()/.ctFitRawEstimate() (R/ctBackendSummary.R), which branch
# on class(fit) once instead of the caller doing it, and had no test coverage
# of any kind before this file.

.chisq_stan_data <- function() {
  data.frame(id = 1, time = seq_along(sunspot.year[1:60]), Y1 = sunspot.year[1:60])
}

test_that("ctChisqTest compares two nested stan fits", {
  skip_on_cran()

  dat <- .chisq_stan_data()
  m1 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0)
  m2 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0, DRIFT = .9)
  f1 <- suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0))
  f2 <- suppressMessages(ctFit(dat, m2, cores = 1, verbose = 0))

  p <- ctChisqTest(f1, f2)
  expect_true(is.numeric(p))
  expect_true(!is.na(p))
  expect_true(p >= 0 && p <= 1)

  # The comparison is between models, not between the argument order they
  # were supplied in -- ctChisqTest() sorts by npars internally.
  expect_equal(p, ctChisqTest(f2, f1))
})

.chisq_julia_model <- function(fixedDrift = NA) {
  drift <- if (is.na(fixedDrift)) "drift" else fixedDrift
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix(drift, 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.chisq_julia_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8; diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(20), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
  data
}

test_that("ctChisqTest compares two nested julia fits, and a julia fit against a stan fit", {
  skip_on_cran()
  skip_without_julia()

  dat <- .chisq_julia_data()
  f_fixed <- suppressMessages(ctFit(dat, .chisq_julia_model(-0.5), backend = "julia", verbose = 0))
  f_free <- suppressMessages(ctFit(dat, .chisq_julia_model(), backend = "julia", verbose = 0))
  expect_s3_class(f_fixed, "ctJuliaFit")
  expect_s3_class(f_free, "ctJuliaFit")

  # f_free has exactly one more free parameter (drift) than f_fixed.
  expect_equal(length(ctsem:::.ctFitRawEstimate(f_free)),
    length(ctsem:::.ctFitRawEstimate(f_fixed)) + 1L)

  p <- ctChisqTest(f_fixed, f_free)
  expect_true(is.numeric(p))
  expect_true(!is.na(p))
  expect_true(p >= 0 && p <= 1)
  expect_equal(p, ctChisqTest(f_free, f_fixed))

  # Cross-backend: same model fit with backend='stan' should be numerically
  # close, since both are maximum-likelihood fits of the same model to the
  # same data, and ctChisqTest() should not error just because the two fits
  # came from different backends.
  f_free_stan <- suppressMessages(ctFit(dat, .chisq_julia_model(), cores = 1, verbose = 0))
  p_cross <- ctChisqTest(f_fixed, f_free_stan)
  expect_true(is.numeric(p_cross))
  expect_true(!is.na(p_cross))
})

test_that("ctChisqTest() refuses a sampled stan fit instead of an opaque row-count error, naming ctLOO()", {
  skip_on_cran()
  # A sampled fit (optimize=FALSE) never sets stanfit$optimfit -- that field
  # belongs to stanoptimis() alone -- so .ctFitOptimValue() returned NULL for
  # it and c(NULL, <value>) silently shrank `ll` to length 1 against `npars`
  # at length 2; data.frame() then failed with "arguments imply differing
  # number of rows", which names neither the sampled fit nor why. See
  # review/J7-sampled-fit-support.md.
  dat <- .chisq_stan_data()
  m1 <- ctModel(type = "dt", LAMBDA = diag(1), MANIFESTVAR = 0)
  f_opt <- suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0))
  f_samp <- suppressWarnings(suppressMessages(ctFit(dat, m1, cores = 1, verbose = 0,
    optimize = FALSE, chains = 1, iter = 60)))
  expect_gt(length(f_samp$stanfit$stanfit@sim), 0)  # genuinely sampled

  err <- tryCatch({
    ctChisqTest(f_opt, f_samp)
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_false(is.na(err))
  expect_match(err, "sampled", fixed = TRUE)
  expect_match(err, "ctLOO", fixed = TRUE)

  # Order shouldn't matter for the refusal either.
  err2 <- tryCatch({
    ctChisqTest(f_samp, f_opt)
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_false(is.na(err2))
  expect_match(err2, "sampled", fixed = TRUE)

  # Two optimized fits: unaffected, still works (regression guard for the
  # early-return check itself).
  p <- ctChisqTest(f_opt, f_opt)
  expect_true(is.numeric(p) && !is.na(p))
})
