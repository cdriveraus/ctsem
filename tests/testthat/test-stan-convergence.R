# Does a stan fit that reports itself finished actually sit at an optimum?
#
# Nine julia test files assert `fit$estimate$converged`; nothing asserted any
# stan-side equivalent, and the suite once passed for months while ordinary
# fits were not converging, because nothing looked (review J11, J15/R8).
#
# The thing not to trust is the termination reason. Both stan optimisers report
# "no step found" as a termination, and the measurement below shows the same
# `terminate$what == 'abs_tol'` on a fit stopped after two iterations, 28 nats
# short, as on the converged one. What discriminates is the gradient at the
# reported optimum, and the optimiser already records it:
# `fit$stanfit$optimfit$ginfn` is its infinity norm. Measured here, converged
# against not: 1.3e-3 against 66.4.
#
# Two small fits, ~2 s total. `ginfn` is free wherever a stan fit already
# exists, so the one-line form of this check also sits in test-knownFits.R.

skip_on_cran()
skip_on_32bit()

.conv_data <- function() {
  gm <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), DRIFT = matrix(-.4), DIFFUSION = matrix(.6),
    MANIFESTVAR = matrix(.2), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(0), Tpoints = 6))
  set.seed(3)
  data.frame(ctGenerate(gm, n.subjects = 15, Tpoints = 6, burnin = 0))
}

.conv_model <- function() {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), MANIFESTMEANS = matrix(0), CINT = matrix("cint")))
  m$pars$indvarying <- FALSE
  m
}

.conv_fit <- function(...) {
  set.seed(1)
  suppressWarnings(suppressMessages(ctFit(.conv_data(), .conv_model(),
    cores = 1, verbose = 0, savescores = FALSE, ...)))
}

test_that("a stan fit reports a gradient that is actually at an optimum", {
  fit <- .conv_fit()
  lpg <- ctsem:::ctOptimDataLpgFunc(fit$stanmodel, fit$standata, cores = 1)$lpg

  atopt <- lpg(fit$stanfit$rawest)
  gradopt <- attr(atopt, "gradient")

  # The recorded number is the thing it claims to be: recomputing the gradient
  # at the reported estimate reproduces `ginfn` exactly. Without this the
  # assertion below would only be testing that the optimiser agrees with
  # itself.
  expect_equal(max(abs(gradopt)), fit$stanfit$optimfit$ginfn,
    tolerance = 1e-8)

  # Measured max|gradient| here is 1.07e-3, identical over three repeats. The
  # bound has ~10x headroom and is four orders of magnitude below the 59.5
  # seen at the starting values below.
  expect_lt(fit$stanfit$optimfit$ginfn, .01)

  # The log posterior at the optimum also has to beat the starting values by a
  # real margin -- measured 163 nats.
  set.seed(1)
  start <- rnorm(length(fit$stanfit$rawest), 0, .01)
  atstart <- lpg(start)
  expect_gt(max(abs(attr(atstart, "gradient"))), 1)
  expect_gt(atopt[1], atstart[1] + 10)
})

test_that("a fit that stopped early terminates the same way and fails the check", {
  # A loose `tol` makes the optimiser stop after two iterations. This is the
  # shape the suite was blind to: a fit object that looks finished. The
  # termination reason is identical to the converged fit's, so a boolean
  # derived from it would report success; the gradient does not.
  stalled <- .conv_fit(optimcontrol = list(tol = 1e6, stochastic = FALSE,
    carefulfit = FALSE))
  converged <- .conv_fit(optimcontrol = list(tol = 1e-8, stochastic = FALSE,
    carefulfit = FALSE))

  expect_identical(stalled$stanfit$optimfit$terminate$what,
    converged$stanfit$optimfit$terminate$what)
  expect_lt(stalled$stanfit$optimfit$iter, converged$stanfit$optimfit$iter)

  # Measured: 66.4 against 1.26e-3, and a log posterior 28 nats worse.
  expect_gt(stalled$stanfit$optimfit$ginfn, 1)
  expect_lt(converged$stanfit$optimfit$ginfn, .01)
  expect_gt(converged$stanfit$optimfit$value,
    stalled$stanfit$optimfit$value + 10)
})
