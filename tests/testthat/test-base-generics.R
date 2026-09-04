# coef(), logLik() and print() have to mean the same thing whichever backend
# produced the fit. They were registered for ctJuliaFit only, so the same call
# worked on one and fell through to the default on the other -- the asymmetry
# was the bug, not the values. These tests check the shape agrees; the julia
# stan-vs-julia value comparisons live in test-stan-julia-parity.R.

.generics_julia_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    set.seed(1)
    n <- 8L
    t <- 6L
    dat <- data.frame(
      id = rep(seq_len(n), each = t),
      time = rep(seq_len(t), times = n),
      Y1 = as.numeric(stats::rnorm(n * t)))
    model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
      manifestNames = 'Y1', latentNames = 'eta1', LAMBDA = matrix(1))
    cached <<- suppressWarnings(suppressMessages(ctFit(dat, model,
      backend = 'julia', cores = 1,
      optimcontrol = list(estonly = TRUE))))
    cached
  }
})

test_that("coef returns a plain numeric raw parameter vector on both backends", {
  stanfit <- ctstantestfit
  stanraw <- coef(stanfit)
  expect_true(is.numeric(stanraw))
  expect_null(dim(stanraw))
  expect_equal(length(stanraw), length(ctsem:::.ctFitRawEstimate(stanfit)))

  skip_on_cran()
  juliafit <- .generics_julia_fit()
  juliaraw <- coef(juliafit)
  expect_true(is.numeric(juliaraw))
  expect_null(dim(juliaraw))
  expect_equal(length(juliaraw), length(ctsem:::.ctFitRawEstimate(juliafit)))
})

test_that("logLik returns a logLik with df and nobs on both backends", {
  stanll <- logLik(ctstantestfit)
  expect_s3_class(stanll, 'logLik')
  expect_equal(length(as.numeric(stanll)), 1L)
  expect_equal(attr(stanll, 'df'), length(coef(ctstantestfit)))
  expect_equal(attr(stanll, 'nobs'), as.integer(ctstantestfit$standata$ndatapoints))
  # AIC() is the reason the attributes have to be there.
  expect_equal(AIC(stanll),
    -2 * as.numeric(stanll) + 2 * attr(stanll, 'df'))

  skip_on_cran()
  juliafit <- .generics_julia_fit()
  juliall <- logLik(juliafit)
  expect_s3_class(juliall, 'logLik')
  expect_equal(length(as.numeric(juliall)), 1L)
  expect_equal(attr(juliall, 'df'), length(coef(juliafit)))
  expect_true(is.finite(attr(juliall, 'nobs')))
  expect_equal(AIC(juliall),
    -2 * as.numeric(juliall) + 2 * attr(juliall, 'df'))
})

test_that("print gives a short summary rather than the whole object", {
  stanlines <- capture.output(print(ctstantestfit))
  expect_lt(length(stanlines), 10L)
  expect_match(stanlines[1], 'ctsem Stan fit')

  skip_on_cran()
  julialines <- capture.output(print(.generics_julia_fit()))
  expect_lt(length(julialines), 10L)
  expect_match(julialines[1], 'ctsem Julia fit')
})
