# Sampling weights (`ctFit(weights = )`, R/ctWeights.R).
#
# The independent answer: with integer weights, the weighted objective is the
# unweighted objective of the data with each subject repeated as many times as
# its weight, up to the normalising factor. So at any parameter vector the
# value, the gradient (by the adjoint, whose chunk rescales its accumulators
# around a weighted subject, and by forward mode), the exact Hessian and each
# subject's score must equal the duplicated data's, the score summed over a
# subject's copies. One comparison per route, at a fixed point, fitting
# nothing; then one weighted fit for what a fit records, and the refusals.

.wt_data <- function() {
  set.seed(7)
  nsub <- 12; nobs <- 6
  d <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    t <- cumsum(stats::runif(nobs, .5, 1.5)); mu <- stats::rnorm(1, 0, 1)
    x <- numeric(nobs); x[1] <- stats::rnorm(1)
    for (j in 2:nobs) x[j] <- x[j - 1] * exp(-0.5 * (t[j] - t[j - 1])) +
      stats::rnorm(1, 0, .6)
    data.frame(id = i, time = t, y1 = mu + x + stats::rnorm(nobs, 0, .5))
  }))
  w <- rep(c(1, 2, 3), length.out = nsub)
  d$w <- w[d$id]
  dup <- do.call(rbind, lapply(seq_len(nsub), function(i) do.call(rbind,
    lapply(seq_len(w[i]), function(k) { s <- d[d$id == i, ]; s$id <- i * 100 + k; s }))))
  list(d = d, dup = dup, w = w, factor = nsub / sum(w))
}

.wt_model <- function() suppressMessages(ctModel(LAMBDA = matrix(1), type = "ct",
  manifestNames = "y1", latentNames = "eta1", MANIFESTMEANS = matrix("mm|||TRUE"),
  CINT = matrix(0), T0MEANS = matrix(0), silent = TRUE))

test_that("a weighted objective is the duplicated data's, on both routes", {
  skip_without_julia()
  data <- .wt_data()
  for (route in c("augmented", "laplace")) {
    spec <- function(dat, ...) structure(suppressWarnings(suppressMessages(ctFit(dat,
      .wt_model(), backend = "julia", intoverpop = route, priors = FALSE,
      fit = FALSE, ...))), class = c("ctJuliaModel", "ctFitModel"))
    weighted <- spec(data$d, weights = "w")
    duplicated <- spec(data$dup)
    # Normalised to a mean of one.
    expect_equal(weighted$subject_weights, data$w * data$factor)
    module <- ctsem:::.ctJuliaModule(weighted$project)
    set.seed(3)
    x <- ctsem:::.ctJuliaNumericVector(stats::rnorm(ctsem:::.ctBackendNpar(weighted), 0, .3))
    get <- function(value) ctsem:::.ctBackendJuliaValue(value)
    evaluate <- function(s, method = "adjoint") get(module$ctsem_evaluate(
      ctsem:::.ctJuliaObjective(s), x, gradient = TRUE, gradient_method = method))
    a <- evaluate(weighted); b <- evaluate(duplicated)
    expect_equal(as.numeric(a$value) / data$factor, as.numeric(b$value), tolerance = 1e-10)
    expect_equal(as.numeric(a$gradient) / data$factor, as.numeric(b$gradient),
      tolerance = 1e-10)
    expect_equal(as.numeric(evaluate(weighted, "forward")$gradient),
      as.numeric(a$gradient), tolerance = 1e-8)
    hessian <- function(s) as.matrix(get(if (route == "laplace")
      module$ctsem_laplace_hessian(ctsem:::.ctJuliaObjective(s), x) else
      module$ctsem_hessian(ctsem:::.ctJuliaObjective(s), x)))
    expect_equal(hessian(weighted) / data$factor, hessian(duplicated), tolerance = 1e-8)
    scores <- function(s) as.matrix(get(module$ctsem_subject_gradients(
      ctsem:::.ctJuliaObjective(s), x, prior = FALSE))$scores)
    expect_equal(unname(scores(weighted) / data$factor),
      unname(rowsum(scores(duplicated), rep(seq_along(data$w), data$w))), tolerance = 1e-10)
  }
})

test_that("a weighted fit records its weighted log likelihood and takes the sandwich", {
  skip_without_julia()
  data <- .wt_data()
  fit <- suppressWarnings(suppressMessages(ctFit(data$d, .wt_model(), backend = "julia",
    intoverpop = "laplace", priors = FALSE, weights = "w", cores = 1, verbose = 0)))
  expect_identical(fit$uncertainty$method, "sandwich")
  expect_equal(fit$estimate$loglik_weighted,
    sum(fit$model_spec$subject_weights * fit$estimate$subject_loglik))
  expect_equal(fit$estimate$loglik, sum(fit$estimate$subject_loglik))
  expect_error(ctLOO(fit), "sampling weights")
  expect_error(ctFitUncertainty(fit, "sample"), "sampling weights")
  # Reported as weighted: no posterior that is only the weighted objective,
  # and no AIC, which a pseudo-likelihood does not have.
  s <- summary(fit)
  # The summary rounds to three decimals.
  expect_equal(s$loglik_weighted, fit$estimate$loglik_weighted, tolerance = 1e-5)
  expect_match(s$weights, "12 subjects")
  expect_null(s$aic)
  expect_null(s$logposterior)
  expect_error(logLik(fit), "sampling weights")
  expect_error(AIC(fit), "sampling weights")
  expect_output(print(fit), "weighted log likelihood")
})

test_that("weights are refused where they do not apply, and must be one per subject", {
  skip_without_julia()
  data <- .wt_data()
  m <- .wt_model()
  expect_error(ctFit(data$d, m, backend = "julia", optimize = FALSE, weights = "w"),
    "optimised fits")
  expect_error(ctFit(data$d, m, backend = "stan", weights = "w"), "backend='julia'")
  varying <- data$d
  varying$w <- ifelse(varying$time > 3, 1, 2)
  expect_error(suppressMessages(ctFit(varying, m, backend = "julia", weights = "w",
    fit = FALSE)), "one value per subject")
  expect_error(suppressMessages(ctFit(data$d, m, backend = "julia", weights = c(1, 2),
    fit = FALSE)), "one value per subject")
  named <- stats::setNames(data$w, seq_along(data$w))
  spec <- suppressMessages(ctFit(data$d, m, backend = "julia", weights = rev(named),
    fit = FALSE))
  expect_equal(spec$subject_weights, data$w * data$factor)
})
