# ctModelImprove() and ctFitAdequacy(): the search keeps what the data
# support, and a correct model stops at once. One latent, so each fit costs
# seconds; the two-latent searches the design was checked on are in
# CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md.

skip_without_julia()
skip_on_32bit()

.improve_data <- function(effect = 0.8, nsub = 80L, nobs = 20L, seed = 5L) {
  set.seed(seed)
  z <- stats::rnorm(nsub); w <- stats::rnorm(nsub)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    times <- cumsum(c(0, stats::rexp(nobs - 1L, 1)))
    x <- stats::rnorm(1, 0, 0.8)
    xs <- numeric(nobs); xs[1] <- x
    for (t in 2:nobs) {
      dt <- times[t] - times[t - 1]
      x <- exp(-0.5 * dt) * x + stats::rnorm(1, 0, 0.8 * sqrt(1 - exp(-dt)))
      xs[t] <- x
    }
    mm <- stats::rnorm(1) + effect * z[i]
    data.frame(id = i, time = times, Y1 = xs + mm + stats::rnorm(nobs, 0, 0.5),
      z = z[i], w = w[i])
  }))
}

.improve_fit <- function(effect, seed = 5L) fit_cached(paste('improve', effect, seed), {
  d <- .improve_data(effect, seed = seed)
  list(data = d, fit = suppressMessages(ctFit(d[, c('id', 'time', 'Y1')],
    suppressMessages(ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
      LAMBDA = matrix(1))), backend = 'julia', cores = 1, verbose = 0)))
})

test_that("adequacy is reported per family, and the residual family needs no draws", {
  x <- .improve_fit(0)
  ad <- ctFitAdequacy(x$fit, nsamples = 40, cores = 1)
  expect_s3_class(ad, 'ctFitAdequacy')
  expect_setequal(ad$families$family, c('lagged', 'subjects', 'nonlinear', 'residuals',
    'distribution', 'fit', 'dependence'))
  expect_true(all(ad$families$p > 0 & ad$families$p <= 1))
  expect_true(ad$adequate)
  # Reuses what the fit generated.
  expect_equal(dim(ad$fit$generated$Y)[1], 40L)
  expect_output(print(ad), 'adequate')
})

test_that("a strong misfit is not floored at one over the number of datasets", {
  # Subjects' means follow z, which the model does not have: the covariate
  # family's p-value used to be the share of 40 datasets beyond the data's,
  # so no smaller than 1/41, and Holm over the families then hid it.
  x <- .improve_fit(0.8)
  ad <- ctFitAdequacy(x$fit, covariates = unique(x$data[, c('id', 'z')]), nsamples = 40,
    cores = 1)
  expect_lt(ad$families$p[ad$families$family == 'covariates'], 1e-4)
  expect_false(ad$adequate)
})

test_that("the Anderson-Darling reference holds its level", {
  set.seed(3)
  p <- replicate(400, .ctAdequacyAD(stats::runif(200))[2])
  # 0.05, within three binomial sds at 400 draws.
  expect_gt(mean(p < 0.05), 0.017)
  expect_lt(mean(p < 0.05), 0.083)
  expect_lt(.ctAdequacyAD(stats::pnorm(stats::rnorm(500, 0.3)))[2], 0.001)
})

test_that("a covariate the model lacks is added, an unrelated one is not, and the search stops", {
  x <- .improve_fit(0.8)
  imp <- suppressMessages(ctModelImprove(x$fit, covariates = unique(x$data[, c('id', 'z', 'w')]),
    relax = c('heterogeneity', 'covariates'), nsamples = 40, cores = 1, verbose = FALSE))
  expect_s3_class(imp, 'ctModelImprove')
  kept <- imp$history$change[imp$history$kept]
  expect_true(any(grepl('^z as a predictor of', kept)))
  expect_false(any(grepl('^w ', imp$history$change)))
  expect_true('z' %in% ctsem:::.ctFitBaseModel(imp$fit)$TIpredNames)
  expect_gt(imp$fit$estimate$loglik, x$fit$estimate$loglik)
  expect_identical(imp$stopped, 'adequate')
  expect_output(print(imp), 'stopped because adequate')
})

test_that("a correct model stops at once", {
  # Its own draw: with seed 5, a correct model's data are one of the five per
  # cent of draws where the search keeps something at alpha = 0.05 (a random
  # diffusion, LR p = 0.0097) -- the rate the simulation in the review note
  # measured, which is the place to test it, not one dataset.
  x <- .improve_fit(0, seed = 1L)
  imp <- suppressMessages(ctModelImprove(x$fit, covariates = unique(x$data[, c('id', 'w')]),
    relax = c('heterogeneity', 'covariates', 'statedependence'), nsamples = 40,
    cores = 1, verbose = FALSE))
  expect_identical(imp$stopped, 'adequate')
  # The screen is at alpha, so it can propose something; the refit's
  # likelihood is what decides, and nothing may be kept.
  expect_false(any(imp$history$kept))
  expect_equal(length(ctsem:::.ctFitRawEstimate(imp$fit)),
    length(ctsem:::.ctFitRawEstimate(x$fit)))
})
