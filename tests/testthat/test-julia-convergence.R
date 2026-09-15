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

# Both of the optimiser's own stopping rules were switched off for every fit in
# the package by a single line, and no test saw it. `.ctJuliaOptimise()` fills
# in a default objective and then asked `is.null(objective)` to decide whether
# the caller had brought one -- which is FALSE from that line onwards, so both
# rules took the state-explicit branch on every route. `gap_tol` was 0, meaning
# nothing stopped a fit whose next step was predicted to gain 1e-8, and
# `stall_window` was 0, which is the whole stall check.
#
# Nothing failed. A rule that is switched off is indistinguishable, from the
# outside, from a rule that never had cause to fire -- which is exactly how the
# post-merge check on the stall work read "never fires on any test fit" as good
# news. So the settings are now reported on the fit and asserted here: this is
# the fit-free half, and it costs one read of a fit two other tests already
# made.
test_that("the stopping rules a marginal fit runs are the ones it asked for", {
  skip_without_julia()
  o <- .jconv_fit()$optim
  # The marginal route certifies, so both rules are live at their defaults.
  expect_equal(o$stall_window, .ctBackendStallWindow(list()))
  expect_gt(o$stall_window, 0L)
  expect_equal(o$gap_tol, .ctBackendInnerGapTol(list(), intoverstates = TRUE))
  expect_gt(o$gap_tol, 0)
  # And a caller who turns one off gets it off, which is the other half of
  # showing the argument reaches the engine rather than a default doing so.
  off <- suppressWarnings(suppressMessages(ctFit(.jconv_data(), .jconv_model(),
    backend = "julia", cores = 1, verbose = 0,
    optimcontrol = list(stallwindow = 0L, innergaptol = 0, maxiter = 5L))))
  expect_equal(off$optim$stall_window, 0L)
  expect_equal(off$optim$gap_tol, 0)
  # And a fit that arrives under its own power is not escaped from. The escape
  # runs only on a stage that stalled or that its own probe says is not a
  # maximum, so a healthy fit pays nothing for it -- the other half of the
  # claim the flat-transform test below makes.
  expect_equal(o$stall_escapes, 0L)
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

# Does any of this recover anything, or does it only report?
#
# Everything above asks what the optimiser says about where it stopped. This
# asks whether it stops somewhere better, and it needs an extreme starting
# value to ask it: from an ordinary start this model fits without incident, so
# a guard that only ever runs on healthy fits can be switched off for a year
# without anything going red. That is not hypothetical. Both of the
# optimiser's stopping rules were off for every fit in the package, for months,
# and the suite was green the whole time -- which is what the test above this
# one exists to stop happening again, and what this one exists to make
# meaningful.
#
# The start puts `drift` at raw 8, inside its own transform's flat region:
# `-log1p_exp(-param)` has a derivative of 3.4e-4 there against 0.5 at zero, so
# the gradient the optimiser sees in that coordinate is three orders down on
# what the rest of the model gives it. The fit climbs, stops, and reports that
# it converged -- at a point its own pullback probe can beat by 215 nats.
#
# Measured on this fixture, 25 subjects and 25 timepoints, local machine:
#
#   escape unavailable   -1021.1986   22 s
#   escape available      -805.8592   39 s, one escape
#
# The control arm turns the *probe* off rather than the escape, because the
# probe is what finds the point: with it off there is nothing to escape to.
# It is also why the worse arm still reports `converged = TRUE` -- `overshot`
# is one of the three things that can falsify that flag, and turning the probe
# off removes it. A converged fit 215 nats low is the failure this is about.
#
# `estonly` so the comparison is of the optimisation alone: the correction and
# uncertainty phases can move an estimate too, and they are not what is under
# test here.
.jconv_flat_data <- function(nsubjects = 25L, ntimes = 25L) {
  set.seed(13)
  baseline <- stats::rnorm(nsubjects, 2, 2)
  t0m <- stats::rnorm(nsubjects, baseline / 2, 1)
  effect <- -log1p(exp(-stats::rnorm(nsubjects, baseline / 2, 0.5)))
  rows <- lapply(seq_len(nsubjects), function(i) {
    gm <- suppressMessages(ctModel(silent = TRUE, Tpoints = ntimes,
      LAMBDA = matrix(1), DRIFT = c(effect[i]), T0MEANS = c(t0m[i]),
      DIFFUSION = c(0.5), MANIFESTVAR = 0.5, T0VAR = c(0),
      CINT = c(baseline[i]), MANIFESTMEANS = 0))
    d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,
      n.subjects = 1, burnin = 0, dtmean = 1, logdtsd = 0)))
    d$id <- i
    d
  })
  do.call(rbind, rows)
}

# Individually varying, nonlinear drift: the transform is the point, and
# `intoverpop = 'laplace'` is the route whose population scales and
# correlations the relative-flatness detector cannot see on its own.
.jconv_flat_model <- function() suppressMessages(ctModel(silent = TRUE,
  type = "ct", CINT = "cint", MANIFESTMEANS = 0, LAMBDA = matrix(1),
  DRIFT = "drift|-log1p_exp(-param)|TRUE"))

test_that("a fit started inside a flat transform gets back out of it", {
  skip_without_julia()
  data <- .jconv_flat_data()
  model <- .jconv_flat_model()
  # Which raw coordinate `drift` is, without paying for a fit to find out.
  spec <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", intoverpop = "laplace", fit = FALSE)))
  npar <- .ctBackendNpar(spec)
  names <- .ctBackendRawParameterNames(list(model_spec = spec), npar)
  drift <- which(names == "drift")
  expect_length(drift, 1L)

  inits <- rep(0, npar)
  inits[drift] <- 8
  # The transform is genuinely flat there, which is the premise of the whole
  # test rather than an incidental detail.
  expect_lt((1 / (1 + exp(8))) / 0.5, 1e-3)

  fit_from <- function(...) suppressWarnings(suppressMessages(
    ctFit(datalong = data, model = model, backend = "julia", cores = 1,
      verbose = 0, intoverpop = "laplace", inits = inits,
      optimcontrol = list(estonly = TRUE, ...))))
  stuck <- fit_from(overshoot = "off", stallwindow = 0L)
  free <- fit_from()

  expect_true(is.finite(stuck$estimate$loglik))
  expect_true(is.finite(free$estimate$loglik))
  # The whole claim, in one line: the same model, the same data and the same
  # starting values land somewhere materially better. 215 nats when measured;
  # a hundred is the bar, so a change that costs most of the effect fails here
  # rather than silently halving it.
  expect_gt(free$estimate$loglik - stuck$estimate$loglik, 100)
  # And by the route this is supposed to take, not by luck.
  expect_gte(free$optim$stall_escapes, 1L)
  expect_true(isTRUE(free$optim$converged))
  # The control really did have nothing to escape with, so the difference is
  # the escape rather than a second optimisation from anywhere at all.
  expect_equal(stuck$optim$stall_escapes, 0L)
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
