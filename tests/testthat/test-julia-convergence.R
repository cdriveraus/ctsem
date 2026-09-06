# Does `fit$estimate$converged` mean what it says?
#
# The stan side has `test-stan-convergence.R`, which checks the same thing the
# same way: the recorded gradient is the gradient at the reported estimate, it
# is small on a fit that arrived, and it is large on one that stopped early.
# This is the julia equivalent, and it exists because the flag was wrong in the
# other direction. Over 64 optimisation replications of a benchmark whose log
# likelihoods matched stan's to the digit, 45 reported `converged = FALSE`.
#
# All of them for the same reason. `converged` was keyed on `saturated` -- a
# raw coordinate whose materialising transform has gone flat -- and saturation
# covers two outcomes that share a zero gradient:
#
#   * the optimizer overstepped into the flat region and stopped somewhere
#     that is not a maximum (a real failure, measured once at 16 log units
#     below the profile peak), and
#   * the data do not identify that coordinate, so it correctly ran to the
#     edge while everything else converged. A population standard deviation
#     with no individual differences behind it is the common case, and it is a
#     finding rather than a fault.
#
# The engine now tells them apart by pulling the flagged coordinate back and
# asking whether the objective improves -- at a maximum, nothing improves it.
# See `_ctsem_overshot` in inst/julia/ContinuousTimeSEM/src/ctsem_backend.jl.
#
# Two fits of one small model, plus one gradient evaluation.

skip_on_cran()
skip_on_32bit()

.jconv_data <- function() {
  # Two manifest means, generated with no individual differences in either.
  # ctModel makes MANIFESTMEANS individually varying by default, so the fitted
  # model asks for two population SDs and a correlation that the data cannot
  # supply -- which is the ordinary shape this used to fail on, not a contrived
  # one.
  gm <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
    LAMBDA = matrix(c(1, 1), 2, 1), DRIFT = matrix(-.4), DIFFUSION = matrix(.6),
    MANIFESTVAR = diag(.2, 2), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(c(1, -.5), 2, 1), Tpoints = 8))
  set.seed(3)
  data.frame(ctGenerate(gm, n.subjects = 25, Tpoints = 8, burnin = 0))
}

.jconv_model <- function() suppressMessages(ctModel(type = "ct", n.latent = 1,
  n.manifest = 2, LAMBDA = matrix(c(1, 1), 2, 1), CINT = matrix(0),
  T0MEANS = matrix(0)))

test_that("a collapsed population scale is not reported as a failure to converge", {
  fit <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0)))
  e <- fit$estimate

  # The precondition. Without this the test below would pass for the boring
  # reason that nothing saturated, and would stop testing anything the day the
  # default changed.
  expect_true(e$saturated)
  expect_true(any(grepl("^popsd_", e$saturated_parameters)))

  # The recorded number is the thing it claims to be, so the assertion after it
  # is not just the optimizer agreeing with itself. Measured: 9.7e-07 against a
  # tolerance of 2.0e-04.
  expect_equal(max(abs(e$gradient)), e$gradient_norm, tolerance = 1e-10)
  expect_lt(e$gradient_norm, e$gradient_tolerance)

  # Pulling the flat coordinate back does not improve the objective, so the
  # point is a maximum and the fit converged. This is the assertion that fails
  # before the fix, where `saturated` alone made `converged` FALSE.
  expect_false(e$overshot)
  expect_equal(e$overshoot_gain, 0)
  expect_true(e$converged)

  # And the flat coordinates are still reported -- as a statement about
  # identification, which is what they are, rather than about convergence.
  expect_gt(fit$identifiability$nweak, 0)
  expect_true(any(grepl("^popsd_", fit$identifiability$parameters)))
  # Rounding-scale negative eigenvalues are not saddles. A symmetric
  # eigendecomposition returns them for a direction whose true curvature is
  # zero, and counting them fired "this is not a maximum" on fits sitting at
  # one.
  expect_equal(fit$identifiability$negative, 0L)
})

test_that("a julia fit stopped early reports a large gradient and does not converge", {
  # The contrast that makes the assertion above mean something: the flag is not
  # simply TRUE everywhere now. Measured: |g| 204.9 against a tolerance of
  # 2.8e-04, and a log likelihood 74 nats worse.
  capped <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(maxiter = 2, carefulfit = FALSE))))
  expect_false(capped$estimate$converged)
  expect_gt(capped$estimate$gradient_norm, 1)
  expect_equal(max(abs(capped$estimate$gradient)),
    capped$estimate$gradient_norm, tolerance = 1e-10)
})
