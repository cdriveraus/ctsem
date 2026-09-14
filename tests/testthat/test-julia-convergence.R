# backend='r' pinned, not left at 'auto'. This test characterises what the model recovers from one particular dataset, and 'auto' prefers the julia engine now, which generates different data for the same seed. The generator is not what is under test here.
# Does `fit$optim$converged` mean what it says?
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
# ## The fixture this file used to use was not identified
#
# It generated data with the manifest means fixed at c(1, -0.5) for every
# subject -- both population standard deviations exactly zero -- and fitted a
# model estimating two of them and their correlation. Two of the tests here
# asserted which optimum that fit reached and which flags came with it, and
# both passed only because of where in the file they ran. Ten starts drawn the
# way `inits = NULL` draws them, at sd 0.01:
#
#   -202.0846   8/10   popsd 0.382, 0.374   correlation ~1.000000
#   -202.2331   2/10   popsd 0.000, 0.000   correlation ~0.999
#
# The reason is in the third column, and it is a property of the model rather
# than of the optimiser or the data. Two *perfectly correlated* manifest-mean
# random effects, on two manifests that both load 1 on a single latent, are a
# latent-level random intercept -- which that model already has, through T0VAR
# and DRIFT. So the between-person spread has two homes and the likelihood has
# a maximum for each: the first solution puts it in the manifest means, the
# second puts it in a slower drift (-0.396 against -0.341) and a larger T0VAR.
# 0.15 nats apart, with 25 subjects and 8 timepoints.
#
# Nothing about convergence can be tested on that. A fit landing in either
# basin is a correct answer to an under-identified question, and an assertion
# about which one is an assertion about the starting draw. The fixture is gone
# rather than pinned to a start that happens to pass: a fixture that only
# behaves from a start we chose is a fixture that is hiding something, and the
# hidden thing here was the model.
#
# What went with it is the only place in the suite where `unidentified` came
# from a real fit. That status is still covered, deterministically and without
# a fit: `test-ctidentify.R` and `test-ctOptimUncertainty.R` drive
# `.ctBackendIntervalCheck()` from information matrices directly, which is
# where the behaviour actually lives.
#
# The fixture below is the same model over data that has the individual
# differences it asks for, drawn independently for the two manifests so the
# correlation is estimable and away from its boundary. One maximum: ten starts,
# ten times -251.8548, spread 5.7e-14 nats, correlation -0.011, population sds
# 0.30 and 0.29. Nothing here pins a starting value.

skip_on_cran()
skip_on_32bit()

# Each subject gets its own manifest means, drawn around c(1, -0.5) and
# independently of each other, so the fitted model's two population SDs and
# their correlation are all things the data can speak to. Generated a subject
# at a time because the individual differences are the point.
.jconv_data <- function() {
  set.seed(3)
  rows <- lapply(seq_len(25), function(i) {
    mm <- c(1, -0.5) + stats::rnorm(2, 0, 0.5)
    gm <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
      LAMBDA = matrix(c(1, 1), 2, 1), DRIFT = matrix(-.4), DIFFUSION = matrix(.6),
      MANIFESTVAR = diag(.2, 2), T0VAR = matrix(1), T0MEANS = matrix(0),
      CINT = matrix(0), MANIFESTMEANS = matrix(mm, 2, 1), Tpoints = 8))
    d <- suppressMessages(data.frame(ctGenerate(gm, n.subjects = 1, Tpoints = 8,
      burnin = 0, backend = 'r')))
    d$id <- i
    d
  })
  do.call(rbind, rows)
}

.jconv_model <- function() suppressMessages(ctModel(type = "ct", n.latent = 1,
  n.manifest = 2, LAMBDA = matrix(c(1, 1), 2, 1), CINT = matrix(0),
  T0MEANS = matrix(0)))

# The fit, computed once. Two tests ask different questions of it and nothing
# here mutates it, so fitting it twice was a minute spent reproducing a number
# we already had.
.jconv_cache <- new.env(parent = emptyenv())
.jconv_fit <- function() {
  if (is.null(.jconv_cache$fit)) {
    .jconv_cache$fit <- suppressWarnings(suppressMessages(
      ctFit(.jconv_data(), .jconv_model(), backend = "julia", cores = 1,
        verbose = 0)))
  }
  .jconv_cache$fit
}

test_that("the optimiser reaches the maximum and says so", {
  fit <- .jconv_fit()
  # `e` is what was estimated, `o` the run that found it.
  e <- fit$estimate
  o <- fit$optim

  # The optimum, in nats, from whatever start the RNG supplies. Ten starts
  # reached this to within 6e-14, which is what makes it assertable without
  # pinning anything -- and what the old fixture could not offer.
  expect_equal(e$loglik, -251.8548, tolerance = 1e-4)
  expect_true(o$converged)
  expect_equal(fit$uncertainty$certification$status, "certified")
  expect_true(fit$uncertainty$certification$certified)

  # Nothing is sitting on a boundary here -- the population SDs are 0.30 and
  # 0.29 and the correlation is -0.011 -- so nothing is flat and there is
  # nothing to overstep. That is the property the old fixture lacked.
  expect_false(o$saturated)
  expect_false(o$overshot)
  expect_equal(o$overshoot_gain, 0)

  # The recorded number is the thing it claims to be, so the assertions around
  # it are not just the optimizer agreeing with itself.
  expect_equal(max(abs(o$gradient)), o$gradient_norm, tolerance = 1e-10)

  # With no flat direction left, the optimiser's own estimate of what remains
  # is accurate rather than meaningless -- `1/2 g'Bg` against an exact gap of
  # 1e-13 here, where on a saturated fit it has read 1.002 against 2e-15. That
  # is the case `_ctsem_optimise_verdict`'s second condition exists for, and
  # `test-backend-optimgap.R` keeps it under unit test.
  expect_true(is.finite(o$predicted_gain))
  expect_lt(o$predicted_gain, o$convergence_tolerance)

  # Rounding-scale negative eigenvalues are not saddles. A symmetric
  # eigendecomposition returns them for a direction whose true curvature is
  # zero, and counting them fired "this is not a maximum" on fits sitting at
  # one.
  expect_equal(fit$identifiability$negative, 0L)
})

test_that("the diagonal metric is not what finds it", {
  # What `optimcontrol$precondition = FALSE` is for, and all this fixture
  # supports saying about it.
  #
  # The claim used to be much stronger -- that without the metric the fit
  # collapses a population scale and comes back `unidentified` -- and it was
  # made on the unidentified fixture, where that describes the starting value
  # rather than the metric: the collapse happened with the metric too, on two
  # starts in ten. On data that identifies the scales, both routes land on the
  # same optimum from any start and the difference is iterations. The counts
  # are in the commit message; pinning them here would fail on any legitimate
  # change to the optimiser and teach whoever hit it nothing.
  plain <- suppressWarnings(suppressMessages(ctFit(.jconv_data(),
    .jconv_model(), backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(precondition = FALSE))))
  expect_equal(plain$estimate$loglik, .jconv_fit()$estimate$loglik,
    tolerance = 1e-4)
  expect_true(plain$optim$converged)
})

test_that("a julia fit stopped early does not converge, and says what is left", {
  # The contrast that makes the assertions above mean something: the flag is
  # not simply TRUE everywhere.
  #
  # One iteration, not two. `maxiter` caps a *stage*, and a fit the curvature
  # says is short is continued from a damped Newton step, so a cap of two now
  # reaches the optimum on many starts -- measured, it certified on some draws
  # and not others, which is a flaky test rather than a weak optimiser.
  capped <- suppressWarnings(suppressMessages(ctFit(.jconv_data(),
    .jconv_model(), backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(maxiter = 1, carefulfit = FALSE))))
  expect_false(capped$optim$converged)
  # Stated as objective still available rather than as a large gradient. With
  # the metric even one iteration brings the gradient down a long way, so
  # "short" and "large gradient" have stopped being the same thing -- which is
  # the reason the criterion moved to nats in the first place.
  cert <- capped$uncertainty$certification
  expect_true(cert$status %in% c("notmaximum", "suboptimal"))
  expect_gt(cert$gap, cert$tolerance)
  expect_equal(max(abs(capped$optim$gradient)),
    capped$optim$gradient_norm, tolerance = 1e-10)
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
  expect_equal(isTRUE(fit$optim$converged),
    certification$status %in% c("certified", "unidentified"))

  # The held complaint is gone rather than left behind to be warned about by a
  # later call on the same fit.
  expect_null(fit$optim$convergence_pending)
})
