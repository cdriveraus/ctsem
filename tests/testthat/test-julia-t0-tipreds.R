# A T0MEANS or T0VAR cell that reads a time independent predictor (julia).
#
# Those cells are evaluated from the subject's parameter vector before any row,
# so a predictor they read is carried at the end of that vector
# (`EKFParameters.ti_tail`, `.ctJuliaResolveStaticRefs()`); the row groups read
# theirs from the row context. Checked here against answers derived in R: with
# one prompt per subject the likelihood is a normal density in closed form, and
# the gradient -- the adjoint's, through the tail -- against central
# differences of the objective, on the plain filter, with a predictor sampled,
# and on the Laplace route.

.t0ti_model <- function(indvarying = FALSE) {
  model <- suppressMessages(ctModel(type = "ct", manifestNames = c("y1", "y2"),
    latentNames = c("e1", "e2"), LAMBDA = diag(2), DRIFT = diag(-1, 2),
    DIFFUSION = diag(0.5, 2), CINT = matrix(0, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(0.6, 2), T0VAR = diag(1, 2),
    # The first reads only the predictor; the second reads it with two PARS.
    T0MEANS = matrix(c("x0 * 1", "a * x0 + b"), 2, 1),
    PARS = c("a|param", "b|param"), TIpredNames = "x0", tipredDefault = FALSE,
    silent = TRUE))
  model$pars$indvarying <- indvarying & model$pars$param %in% "b"
  model
}

.t0ti_data <- function(n = 12L, missing = integer()) {
  set.seed(4)
  x0 <- round(seq(-1.5, 1.5, length.out = n), 3)
  d <- data.frame(id = seq_len(n), time = 0, x0 = x0,
    y1 = x0 + stats::rnorm(n, 0, 1.2), y2 = 0.7 * x0 + 0.3 + stats::rnorm(n, 0, 1.2))
  d$x0[missing] <- NA
  d
}

.t0ti_spec <- function(data, ...) structure(suppressMessages(ctFit(data,
  .t0ti_model(...), backend = "julia", priors = FALSE, fit = FALSE, cores = 1)),
  class = c("ctJuliaModel", "ctFitModel"))

.t0ti_central <- function(spec, x, h = 1e-5) vapply(seq_along(x), function(k) {
  e <- replace(numeric(length(x)), k, h)
  (ctJuliaEvaluate(spec, x + e, gradient = FALSE, cores = 1)$value -
    ctJuliaEvaluate(spec, x - e, gradient = FALSE, cores = 1)$value) / (2 * h)
}, numeric(1))

test_that("a T0MEANS cell reads the subject's TI predictor, alone or with PARS", {
  skip_without_julia()
  data <- .t0ti_data()
  spec <- .t0ti_spec(data)
  # The two cells resolve to transforms reading the predictor, the first with
  # no parameter behind it.
  pt <- as.data.frame(spec$parameter_table)
  t0 <- pt[pt$matrix == "T0MEANS", ]
  expect_true(all(grepl("tipreds[1]", t0$transform, fixed = TRUE)))
  expect_true(is.na(t0$parnumber[t0$row == 1]))
  expect_equal(.ctBackendNpar(spec), 2L)
  x <- c(0.4, -0.2)
  out <- ctJuliaEvaluate(spec, x, gradient = TRUE, cores = 1)
  # One prompt each: y ~ N(T0MEANS, T0VAR + MANIFESTVAR), the two independent.
  sd <- sqrt(1 + 0.6^2)
  closed <- sum(stats::dnorm(data$y1, data$x0, sd, log = TRUE) +
    stats::dnorm(data$y2, x[1] * data$x0 + x[2], sd, log = TRUE))
  expect_equal(out$value, closed, tolerance = 1e-10)
  expect_equal(as.numeric(out$gradient), c(
    sum((data$y2 - x[1] * data$x0 - x[2]) * data$x0) / sd^2,
    sum(data$y2 - x[1] * data$x0 - x[2]) / sd^2), tolerance = 1e-8)
})

test_that("a sampled TI predictor read by a T0MEANS cell gets its derivative", {
  skip_without_julia()
  data <- .t0ti_data(missing = c(3L, 8L))
  spec <- structure(suppressMessages(ctFit(data, .t0ti_model(), backend = "julia",
    priors = FALSE, fit = FALSE, cores = 1,
    optimcontrol = list(tipredMissing = "sample"))), class = c("ctJuliaModel", "ctFitModel"))
  npar <- .ctBackendNpar(spec)
  # Two PARS and the two missing values.
  expect_equal(npar, 4L)
  x <- c(0.4, -0.2, 0.5, -0.9)
  out <- ctJuliaEvaluate(spec, x, gradient = TRUE, cores = 1)
  expect_true(all(is.finite(out$gradient)))
  expect_equal(as.numeric(out$gradient), .t0ti_central(spec, x), tolerance = 1e-6)
  # The sampled values do read as the predictor: both cells move with them.
  expect_gt(min(abs(out$gradient[3:4])), 1e-3)
})

test_that("on the Laplace route the T0MEANS predictor reaches value and gradient", {
  skip_without_julia()
  data <- .t0ti_data(n = 16L)
  # Three prompts each, so a person's effect on b is identified.
  data <- do.call(rbind, lapply(0:2, function(k) transform(data, time = k * 0.5,
    y1 = y1 + stats::rnorm(nrow(data), 0, 0.3), y2 = y2 + stats::rnorm(nrow(data), 0, 0.3))))
  data <- data[order(data$id, data$time), ]
  spec <- structure(suppressMessages(ctFit(data, .t0ti_model(indvarying = TRUE),
    backend = "julia", intoverpop = "laplace", priors = FALSE, fit = FALSE, cores = 1)),
    class = c("ctJuliaModel", "ctFitModel"))
  npar <- .ctBackendNpar(spec)
  x <- c(0.4, -0.2, rep(-0.5, npar - 2L))
  out <- ctJuliaEvaluate(spec, x, gradient = TRUE, cores = 1)
  expect_true(is.finite(out$value))
  expect_equal(as.numeric(out$gradient), .t0ti_central(spec, x), tolerance = 1e-5)
  # Not the predictor's absence: moving it moves the objective.
  shifted <- transform(data, x0 = x0 + 1)
  spec2 <- structure(suppressMessages(ctFit(shifted, .t0ti_model(indvarying = TRUE),
    backend = "julia", intoverpop = "laplace", priors = FALSE, fit = FALSE, cores = 1)),
    class = c("ctJuliaModel", "ctFitModel"))
  expect_gt(abs(ctJuliaEvaluate(spec2, x, gradient = FALSE, cores = 1)$value - out$value), 1)
})

test_that("a model whose T0MEANS reads a TI predictor fits, not only evaluates", {
  skip_without_julia()
  # The fitted specification carries no `standata`; the predictor count the
  # engine needs came from there and stopped every such fit before its first
  # evaluation, while the tests above, on `fit = FALSE`, passed.
  data <- .t0ti_data(n = 40L)
  fit <- suppressMessages(ctFit(data, .t0ti_model(), backend = "julia", priors = FALSE,
    cores = 1, verbose = 0))
  expect_true(is.finite(fit$estimate$loglik))
  spec <- .t0ti_spec(data)
  expect_equal(fit$estimate$loglik, ctJuliaEvaluate(spec, ctsem:::.ctFitRawEstimate(fit),
    gradient = FALSE, cores = 1)$value, tolerance = 1e-6)
})
