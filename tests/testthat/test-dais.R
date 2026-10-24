# Doubly adaptive importance sampling as a stage of 'checked': the normal
# moved to the posterior's moments before anything is sampled. The method
# runs in the engine (`ctsem_dais`, dais.jl), whose tests hold it to closed
# forms (test_dais.jl: a skewed target's moments, a held direction, unusable
# draws); here, that 'checked' runs it on a julia fit and reads what it
# returns, and the sampler's metric.

test_that("'checked' runs DAIS in the engine on a julia fit, and not where asked not to", {
  skip_without_julia()
  # The shared Laplace fixture at 100 subjects, whose check names directions.
  fit <- laplace_fixture(nsub = 100L, tp = 6L, seed = 31L)
  set.seed(14)
  u <- suppressWarnings(ctFitUncertainty(fit, "checked", cores = 1))
  d <- u$uncertainty$details
  expect_gt(NCOL(d$normalcheck$directions), 0L)
  dais <- d$dais
  expect_true(isTRUE(dais$used))
  expect_gt(dais$iterations, 0L)
  expect_equal(dais$gradients, dais$draws * dais$iterations)
  expect_equal(dais$cost, 2 * dais$gradients)
  expect_equal(nrow(dais$trace), dais$iterations)
  expect_true(all(c("gamma", "natural", "noise", "usable") %in% colnames(dais$trace)))
  npar <- length(fit$estimate$raw)
  expect_length(dais$centre, npar)
  expect_equal(dim(dais$cov), c(npar, npar))
  expect_true(is.numeric(dais$check$directions))
  # The same seed, the same stage.
  set.seed(14)
  again <- suppressWarnings(ctFitUncertainty(fit, "checked", cores = 1))
  expect_equal(again$uncertainty$details$dais$cov, dais$cov)
  # Off by name: no stage, the check's normal as before.
  set.seed(14)
  off <- suppressWarnings(ctFitUncertainty(fit, "checked", cores = 1,
    control = list(dais = FALSE)))
  expect_null(off$uncertainty$details$dais)
})

test_that("a density without an engine has no DAIS stage", {
  # A closed form built in R, as a stan fit's density is: nothing to run it on.
  lpg <- ctsem:::.ctImisPointGradbatch(function(x) {
    out <- -0.5 * sum(x^2)
    attr(out, 'gradient') <- -x
    out
  })
  expect_null(ctsem:::.ctDaisEngine(lpg, rep(0, 3), diag(3)))
})

test_that("the sampler's metric is the DAIS normal's precision where there is one", {
  set.seed(17)
  C <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  fit <- list(uncertainty = list(hessian = -diag(5), details = list(dais = list(cov = C),
    normalcheck = list(directions = matrix(1, 5, 1), widths = 4))))
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "dais")
  expect_equal(-ctsem:::.ctBackendSampleMetric(fit), solve(C), tolerance = 1e-8)
  # The importance draws' covariance comes first where they were used,
  # weighted, and number at least twice the dimension in effective size; a
  # smaller effective size, or draws not used, leave DAIS's.
  Cis <- crossprod(matrix(rnorm(25), 5)) + 2 * diag(5)
  fit$estimate <- list(raw = rep(0, 5), cov = Cis)
  fit$uncertainty$details$importance_sampling <- list(used = TRUE, weighted = TRUE, ess = 50)
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "is")
  expect_equal(-ctsem:::.ctBackendSampleMetric(fit), solve(Cis), tolerance = 1e-8)
  fit$uncertainty$details$importance_sampling$ess <- 9
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "dais")
  fit$uncertainty$details$importance_sampling <- list(used = FALSE, weighted = TRUE, ess = 50)
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "dais")
  fit$estimate <- NULL
  # A direction the covariance does not hold gets no precision, so the engine
  # floors it wide, as it does for the Hessian's unidentified directions.
  v <- c(1, 0, 0, 0, 0)
  Cnull <- C - C %*% v %*% t(v) %*% C / as.numeric(t(v) %*% C %*% v)
  fit$uncertainty$details$dais$cov <- Cnull
  Pn <- -ctsem:::.ctBackendSampleMetric(fit)
  expect_equal(as.numeric(Pn %*% v), rep(0, 5), tolerance = 1e-6)
  expect_equal(Pn %*% Cnull %*% Pn, Pn, tolerance = 1e-6)
  fit$uncertainty$details$dais <- NULL
  # Without one, the Hessian rescaled to the check's widths, as before.
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "check")
  expect_false(isTRUE(all.equal(-ctsem:::.ctBackendSampleMetric(fit), solve(C))))
  fit$uncertainty$details$normalcheck <- NULL
  expect_identical(ctsem:::.ctBackendMetricSource(fit), "hessian")
  expect_null(ctsem:::.ctBackendSampleMetric(fit))
})
