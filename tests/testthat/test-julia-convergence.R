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
# The engine now tells them apart by pulling coordinates back and asking
# whether the objective improves -- at a maximum, nothing improves it. See
# `_ctsem_overshot` in inst/julia/ContinuousTimeSEM/src/ctsem_backend.jl.
#
# ## Every fit here starts from zeros
#
# It did not, and that made two of these tests assert a coin flip. `inits =
# NULL` draws `rnorm(npar, 0, .01)` from whatever RNG state the preceding tests
# left, and on this fixture the starting point decides which of two optima the
# fit reaches. Over four seeded starts, with the diagonal metric and without:
#
#   start   with metric                  without
#     1     -202.0846  certified         -202.2331  unidentified
#     2     -202.2331  unidentified      -202.2331  certified
#     3     -202.0846  certified         -202.2331  certified
#     4     -202.0846  certified         -202.2331  certified
#
# So `expect_false(e$saturated)` failed on start 2, and the file's second test,
# which asserted the opposite of it, failed on three starts of four. Both
# passed in practice only because of where in the file they ran. From zeros --
# the neutral start, not the one that passes -- every fit below is reproducible
# to the digit.
#
# Three fits of one small model: the full fit is computed once and shared, and
# the other two are the contrasts that make its assertions mean something.

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

# Computed once each: three tests want the same start and two want the same
# fit, and nothing here mutates either.
.jconv_cache <- new.env(parent = emptyenv())

# Zeros of the right length. `fit = FALSE` builds the parameter table without
# optimising anything, which is where the count comes from.
.jconv_zeros <- function() {
  if (is.null(.jconv_cache$zeros)) {
    spec <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
      backend = "julia", fit = FALSE)))
    .jconv_cache$zeros <- rep(0, ctsem:::.ctBackendNpar(spec))
  }
  .jconv_cache$zeros
}

# The full fit, computed once: fitting it twice was a minute spent reproducing
# a number we already had.
.jconv_fit <- function() {
  if (is.null(.jconv_cache$fit)) {
    .jconv_cache$fit <- suppressWarnings(suppressMessages(
      ctFit(.jconv_data(), .jconv_model(), backend = "julia", cores = 1,
        verbose = 0, inits = .jconv_zeros())))
  }
  .jconv_cache$fit
}

test_that("the optimiser reaches a maximum and says so", {
  fit <- .jconv_fit()
  e <- fit$estimate

  # The optimum, in nats. This is the assertion that guards the optimiser: the
  # other basin on this fixture is -202.2331, so a change that sent the fit
  # there fails here by 0.15 rather than by a flag that may or may not have
  # tripped. -202.0846 from zeros, 265 iterations, reproducible.
  expect_equal(e$loglik, -202.0846, tolerance = 1e-4)
  expect_true(e$converged)
  expect_equal(fit$uncertainty$certification$status, "certified")
  expect_true(fit$uncertainty$certification$certified)

  # Nothing flat at this optimum, so nothing to overstep. Asserted from a
  # pinned start, where it is a property of the estimate rather than of the
  # draw that reached it.
  expect_false(e$saturated)
  expect_false(e$overshot)
  expect_equal(e$overshoot_gain, 0)

  # The recorded number is the thing it claims to be, so the assertions around
  # it are not just the optimizer agreeing with itself.
  expect_equal(max(abs(e$gradient)), e$gradient_norm, tolerance = 1e-10)

  # With no flat direction left, the optimiser's own estimate of what remains
  # is accurate rather than meaningless -- `1/2 g'Bg` against an exact gap of
  # 3.6e-13 here, where on a saturated fit it has read 1.002 against 2e-15.
  # That is the case `_ctsem_optimise_verdict`'s second condition exists for,
  # and `test-backend-optimgap.R` keeps it under unit test.
  expect_true(is.finite(e$predicted_gain))
  expect_lt(e$predicted_gain, e$convergence_tolerance)

  # Rounding-scale negative eigenvalues are not saddles. A symmetric
  # eigendecomposition returns them for a direction whose true curvature is
  # zero, and counting them fired "this is not a maximum" on fits sitting at
  # one.
  expect_equal(fit$identifiability$negative, 0L)
})

test_that("the fit is the same one without the diagonal metric", {
  # What `optimcontrol$precondition = FALSE` is for, and the only claim about
  # the metric this fixture actually supports.
  #
  # It used to claim more: that without the metric this fit collapses a
  # population scale and comes back `unidentified`. From a pinned start it does
  # not -- both routes reach -202.0846 and both certify, and the difference is
  # 265 iterations against 344. The collapse the old test asserted was the
  # random start landing in the other basin, which happens with the metric too
  # (one start in four, above). Asserting it was asserting the draw.
  #
  # The iteration count is not asserted either. It is a real measurement and it
  # is in the commit message, but pinning 265 <= 344 would fail on any
  # legitimate change to the optimiser and teach whoever hit it nothing. What
  # is worth holding is that switching the metric off still fits, and fits the
  # same model rather than a differently conditioned approximation to it.
  plain <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0, inits = .jconv_zeros(),
    optimcontrol = list(precondition = FALSE))))
  expect_equal(plain$estimate$loglik, .jconv_fit()$estimate$loglik,
    tolerance = 1e-4)
  expect_true(plain$estimate$converged)
})

test_that("a julia fit stopped early does not converge, and says what is left", {
  # The contrast that makes the assertions above mean something: the flag is
  # not simply TRUE everywhere.
  #
  # One iteration, not two. `maxiter` caps a *stage*, and a fit the curvature
  # says is short is continued from a damped Newton step, so a cap of two now
  # reaches the optimum on many starts -- measured, it certified on some draws
  # and not others, which is a flaky test rather than a weak optimiser.
  capped <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0, inits = .jconv_zeros(),
    optimcontrol = list(maxiter = 1, carefulfit = FALSE))))
  expect_false(capped$estimate$converged)
  # Stated as objective still available rather than as a large gradient. With
  # the metric even one iteration brings the gradient down a long way, so
  # "short" and "large gradient" have stopped being the same thing -- which is
  # the reason the criterion moved to nats in the first place.
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
  fit <- .jconv_fit()
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
