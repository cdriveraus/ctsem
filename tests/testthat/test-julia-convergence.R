# backend='r' pinned, not left at 'auto'. This test characterises what the model recovers from one particular dataset, and 'auto' prefers the julia engine now, which generates different data for the same seed. The generator is not what is under test here.
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
  data.frame(ctGenerate(gm, n.subjects = 25, Tpoints = 8, burnin = 0,
    backend = 'r'))
}

.jconv_model <- function() suppressMessages(ctModel(type = "ct", n.latent = 1,
  n.manifest = 2, LAMBDA = matrix(c(1, 1), 2, 1), CINT = matrix(0),
  T0MEANS = matrix(0)))

test_that("the optimiser reaches a maximum instead of collapsing a population scale", {
  fit <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0)))
  e <- fit$estimate

  # This fixture used to end with its population standard deviation in the flat
  # region of its own transform, and the file was built around explaining why
  # that still counted as converged. It was not the data: it was the optimiser
  # stepping in raw coordinates whose scales differ by a factor of ten, with
  # only a scalar metric to go on. Measured over eight fits of this model,
  # with the diagonal metric and without:
  #
  #                 best ll      saturated   certification
  #   with metric   -202.0846        0/8     certified x8
  #   without       -202.2331        8/8     unidentified x8
  #
  # So the collapse is now the contrast rather than the subject. Asserted
  # directly, because a regression in the metric would bring it straight back
  # and nothing else in this file would notice.
  expect_false(e$saturated)
  expect_equal(fit$uncertainty$certification$status, "certified")
  expect_true(fit$uncertainty$certification$certified)
  expect_true(e$converged)

  # The recorded number is the thing it claims to be, so the assertions around
  # it are not just the optimizer agreeing with itself.
  expect_equal(max(abs(e$gradient)), e$gradient_norm, tolerance = 1e-10)

  # With no flat direction left, the optimiser's own estimate of what remains
  # is accurate rather than meaningless -- `1/2 g'Bg` is 8e-12 here against an
  # exact gap of 3e-16, where on the saturated fit it read 1.002 against 2e-15.
  # That is the case `_ctsem_optimise_verdict`'s second condition exists for,
  # and `test-backend-optimgap.R` keeps it under unit test now that a real fit
  # no longer produces it.
  expect_true(is.finite(e$predicted_gain))
  expect_lt(e$predicted_gain, e$convergence_tolerance)

  # Nothing saturated, so nothing to overstep.
  expect_false(e$overshot)
  expect_equal(e$overshoot_gain, 0)

  # Rounding-scale negative eigenvalues are not saddles. A symmetric
  # eigendecomposition returns them for a direction whose true curvature is
  # zero, and counting them fired "this is not a maximum" on fits sitting at
  # one.
  expect_equal(fit$identifiability$negative, 0L)
})

test_that("without the metric the same fit collapses its population scale", {
  # The other half of the measurement above, and the reason the assertions in
  # the previous test are worth making: the difference is the metric and not
  # the model, the data or the seed.
  fit <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(precondition = FALSE))))
  e <- fit$estimate
  expect_true(e$saturated)
  expect_true(any(grepl("^popsd_", e$saturated_parameters)))
  # And it is still a maximum with an unidentified coordinate rather than a
  # failed fit -- the distinction this file has always been about.
  expect_false(e$overshot)
  expect_true(e$converged)
  expect_equal(fit$uncertainty$certification$status, "unidentified")
})

test_that("a julia fit stopped early reports a large gradient and does not converge", {
  # The contrast that makes the assertion above mean something: the flag is not
  # simply TRUE everywhere.
  #
  # One iteration, not two. `maxiter` caps a *stage*, and a fit the curvature
  # says is short is continued from a damped Newton step, so a cap of two now
  # reaches the optimum on many starts -- measured, it certified on some draws
  # and not others, which is a flaky test rather than a weak optimiser. At one
  # iteration the fit is decisively short: |g| around 7, an exact gap of 5e-3,
  # and `notmaximum`.
  capped <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(maxiter = 1, carefulfit = FALSE))))
  expect_false(capped$estimate$converged)
  # Stated as objective still available rather than as a large gradient. With
  # the metric even one iteration brings the gradient to around 0.1, so "short"
  # and "large gradient" have stopped being the same thing -- which is the
  # reason the criterion moved to nats in the first place. The exact gap is
  # 5e-3 here against a tolerance of 1e-6.
  cert <- capped$uncertainty$certification
  expect_true(cert$status %in% c("notmaximum", "suboptimal"))
  expect_gt(cert$gap, cert$tolerance)
  expect_equal(max(abs(capped$estimate$gradient)),
    capped$estimate$gradient_norm, tolerance = 1e-10)
})


test_that("what the fit reports is what the curvature measured", {
  # The complaint: a fit could pass certification -- the optimum bounded within
  # `gaptol` of the estimate, no warning raised, the summary saying so -- and
  # still report `converged = FALSE`, because that field was the optimiser's
  # own gradient verdict and nothing ever revisited it. Two verdicts on one
  # fit, and the weaker one was the one a user reads.
  fit <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0)))
  certification <- fit$uncertainty$certification

  # The precondition: something was measured. Without this the assertion below
  # would hold vacuously on a fit that certified nothing.
  expect_true(is.list(certification) && length(certification$status) == 1L)
  expect_true(is.finite(certification$gap))

  # The one rule, on a real fit: converged iff the curvature says this is a
  # maximum. `unidentified` is a maximum with a coordinate the data does not
  # determine; the rest of the not-certified statuses are not maxima.
  expect_equal(isTRUE(fit$estimate$converged),
    certification$status %in% c("certified", "unidentified"))

  # The held complaint is gone rather than left behind to be warned about by a
  # later call on the same fit.
  expect_null(fit$estimate$convergence_pending)
})
