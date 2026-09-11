# Convergence certification: the arithmetic and the three verdicts.
#
# All of it without a fit. The quantity is a function of a Hessian and a
# gradient, so a constructed pair says whether it is right, and a constructed
# pair is the only way to put a genuinely flat or genuinely negative direction
# in front of it on demand.

context("backend-optimgap")

test_that("the gap is the predicted improvement, exactly, for a quadratic", {
  # A quadratic in two coordinates with known curvature. Its optimum is at
  # `H^-1 g` from here and the log likelihood there is higher by `0.5 g'H^-1g`,
  # so the gap is not an estimate of anything -- it is that number.
  H <- diag(c(4, 25))
  g <- c(2, 5)
  out <- ctsem:::.ctBackendOptimGap(-H, g)
  expect_true(out$ok)
  expect_equal(out$gap, 0.5 * (2^2 / 4 + 5^2 / 25))
  expect_equal(out$lambda, sqrt(2 * out$gap))
  # And the step is the displacement that realises it.
  expect_equal(out$step, c(2 / 4, 5 / 25))
  expect_equal(out$ntrusted, 2L)
  expect_equal(out$nflat, 0L)
  expect_equal(out$nnegative, 0L)
})

test_that("the gap is invariant to rescaling a parameter, and the gradient is not", {
  # The whole reason for the quantity. Rescale a coordinate and the gradient
  # and the curvature move in opposite directions; the gap does not move at
  # all, while `max|g|` -- the current criterion -- changes by the same factor.
  H <- matrix(c(4, 1, 1, 25), 2, 2)
  g <- c(2, 5)
  A <- diag(c(1, 1000))            # theta -> A theta
  Hs <- solve(A) %*% H %*% solve(A)
  gs <- as.numeric(solve(A) %*% g)
  bare <- ctsem:::.ctBackendOptimGap(-H, g)
  scaled <- ctsem:::.ctBackendOptimGap(-Hs, gs)
  expect_equal(scaled$gap, bare$gap)
  expect_equal(scaled$lambda, bare$lambda)
  expect_false(isTRUE(all.equal(max(abs(gs)), max(abs(g)))))
})

test_that("a flat direction is excluded from the gap and its gradient kept", {
  # The unsafe version of this computes the gap on the positive-curvature
  # subspace and stops. Here the second coordinate has no curvature and a live
  # gradient: the gap must not swallow it, and the residual must carry it out
  # whole so the caller can measure what it is worth.
  H <- diag(c(4, 0))
  g <- c(2, 7)
  out <- ctsem:::.ctBackendOptimGap(-H, g)
  expect_equal(out$ntrusted, 1L)
  expect_equal(out$nflat, 1L)
  expect_equal(out$gap, 0.5 * 2^2 / 4)
  expect_equal(out$residual, c(0, 7))
  expect_equal(out$residual_norm, 7)
  # And the Newton step stays inside the trusted subspace rather than dividing
  # by a curvature of zero.
  expect_true(all(is.finite(out$step)))
  expect_equal(out$step[2L], 0)
})

test_that("negative curvature is counted and reported as not a maximum", {
  H <- diag(c(4, -9))
  out <- ctsem:::.ctBackendOptimGap(-H, c(1, 1))
  expect_equal(out$nnegative, 1L)
  verdict <- ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01)
  expect_equal(verdict$status, "notmaximum")
  expect_false(verdict$certified)
  expect_match(verdict$reason, "negative curvature")
})

test_that("a small gap with a live flat direction is not certified", {
  # The case the projection alone would certify: nothing to gain where the
  # curvature is trusted, and a whole log likelihood unit available where it
  # is not.
  out <- ctsem:::.ctBackendOptimGap(-diag(c(4, 0)), c(1e-8, 7))
  expect_lt(out$gap, 1e-10)
  verdict <- ctsem:::.ctBackendCertify(out, probe = list(gain = 1.4),
    tolerance = 0.01)
  expect_equal(verdict$status, "uncertified")
  expect_match(verdict$reason, "no curvature")
  # The same point with nothing to be had along that direction is certified:
  # a flat direction is not by itself a failure, it is a flat direction.
  quiet <- ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01)
  expect_equal(quiet$status, "certified")
})

test_that("a saturated transform is never certified by a gradient that underflowed", {
  # A saturated coordinate reports no gradient and no curvature, so it passes
  # every test here for the wrong reason. The fit's own verdict is what stops
  # that, which is why it is consulted rather than recomputed.
  out <- ctsem:::.ctBackendOptimGap(-diag(c(4, 0)), c(1e-9, 0))
  expect_equal(ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01, saturated = TRUE)$status, "uncertified")
  expect_equal(ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01, overshot = TRUE)$status, "notmaximum")
  expect_equal(ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01)$status, "certified")
})

test_that("a gap above tolerance says how far, in log likelihood", {
  out <- ctsem:::.ctBackendOptimGap(-diag(c(4, 25)), c(2, 5))
  verdict <- ctsem:::.ctBackendCertify(out, probe = list(gain = 0),
    tolerance = 0.01)
  expect_equal(verdict$status, "suboptimal")
  expect_match(verdict$reason, "log likelihood above this estimate")
  expect_false(verdict$certified)
})

test_that("the probe reports the best actual gain, and zero when there is none", {
  # A likelihood rising along the direction: the probe has to find it and say
  # how far it went, because that number is what a norm of the gradient cannot
  # give.
  rising <- function(x) 3 * x[2L]
  found <- ctsem:::.ctBackendOptimGapProbe(rising, at = c(0, 0),
    direction = c(0, 1), value = 0)
  expect_equal(found$gain, 12)      # the longest step, 4, times 3
  expect_equal(found$length, 4)
  falling <- function(x) -3 * x[2L]
  none <- ctsem:::.ctBackendOptimGapProbe(falling, at = c(0, 0),
    direction = c(0, 1), value = 0)
  expect_equal(none$gain, 0)
  # A direction of length zero is nothing to probe, not an error.
  expect_equal(ctsem:::.ctBackendOptimGapProbe(rising, c(0, 0), c(0, 0), 0)$gain, 0)
  # An objective that cannot be evaluated there is a step not taken, not a
  # failed certification.
  broken <- function(x) stop("no")
  expect_equal(ctsem:::.ctBackendOptimGapProbe(broken, c(0, 0), c(0, 1), 0)$gain, 0)
})

test_that("a hessian that cannot be used says so rather than certifying", {
  expect_false(ctsem:::.ctBackendOptimGap(NULL, 1)$ok)
  expect_false(ctsem:::.ctBackendOptimGap(matrix(NA_real_, 2, 2), c(1, 1))$ok)
  # Wrong length gradient is a caller error, not a verdict.
  expect_false(ctsem:::.ctBackendOptimGap(-diag(2), c(1, 1, 1))$ok)
  verdict <- ctsem:::.ctBackendCertify(ctsem:::.ctBackendOptimGap(NULL, 1))
  expect_equal(verdict$status, "unknown")
  expect_false(verdict$certified)
})

test_that("the trusted subspace uses the same rule the interval check does", {
  # Two rules for "this direction has no curvature" that can disagree is how a
  # fit comes to be described one way by its intervals and another by its
  # convergence. Same relative tolerance, same answer on the same matrix.
  information <- diag(c(1, 1e-14))
  split <- ctsem:::.ctBackendInformationSplit(-information)
  mass <- ctsem:::.ctBackendNullMass(information)
  expect_equal(sum(!split$trusted), 1L)
  expect_gt(mass[2L], 0.9)
  expect_lt(mass[1L], 1e-6)
})

test_that("a gap is never reported through fixed-point rounding", {
  # Found end to end: `summary()` rounds every numeric element to `digits`, so
  # a certified fit whose gap was 1.9e-04 reported `optimgap = 0` -- a
  # measurement turned into a claim of exactness, on the one quantity this
  # whole mechanism exists to state honestly. A gap worth reporting spans many
  # orders of magnitude, so it goes out as a sentence and the exact value stays
  # on the fit.
  expect_equal(round(1.911e-04, 3), 0)
  note <- paste0("Predicted objective still available at this estimate: ",
    signif(1.911e-04, 3))
  expect_equal(roundSummaryCtStanFitValue(note, digits = 3), note)
  expect_match(note, "0.000191")
})

test_that("the damped step backtracks as far as the arithmetic allows", {
  # The case that refused a correction in a real fit: an ascent direction whose
  # step is far too long, because the trusted curvature spans nine orders and
  # the Newton step is 50 units in coordinates where the parameters are order
  # one. A four-rung ladder stopped at 1/8, which still overshot; 1/16 improved.
  # Nothing here is tuned to that fit -- the floor is the objective's own
  # resolution, so a differently conditioned model simply takes a different
  # number of rungs.
  #
  # A quadratic with its optimum a sixteenth of the way along the step.
  peak <- 1 / 16
  objective <- function(x) -100 * (x[1L] - peak)^2
  stepped <- ctsem:::.ctBackendDampedStep(objective, at = 0, step = 1,
    value = objective(0), directional = 2 * 100 * peak^2)
  expect_false(is.null(stepped$accepted))
  expect_lte(stepped$accepted$alpha, 0.125)
  expect_gt(stepped$accepted$value, objective(0))
  expect_gt(stepped$achievable, 0)
})

test_that("a direction that offers nothing is refused, and says how much", {
  falling <- function(x) -abs(x[1L])
  stepped <- ctsem:::.ctBackendDampedStep(falling, at = 0, step = 1,
    value = 0, directional = 1)
  expect_null(stepped$accepted)
  expect_equal(stepped$achievable, 0)
  # A direction that is not an ascent direction is not stepped along at all:
  # the curvature said this cannot help, so the objective is never called.
  called <- 0
  counted <- function(x) { called <<- called + 1; 0 }
  expect_null(ctsem:::.ctBackendDampedStep(counted, 0, 1, 0, -1)$accepted)
  expect_equal(called, 0)
})

test_that("an increase the objective cannot represent is not a correction", {
  # Armijo alone does not rule this out, and a first draft of this test assumed
  # it did: sufficient increase scales with the step, so an increase of 1e-14
  # satisfies it once alpha is around 1e-10. Accepting that spends a resumed
  # optimisation on a point indistinguishable from the one it started at, so
  # acceptance also requires the increase to exceed the objective's resolution.
  crumbs <- function(x) if (x[1L] > 0) -2705 + 1e-20 else -2705
  stepped <- ctsem:::.ctBackendDampedStep(crumbs, at = 0, step = 1,
    value = -2705, directional = 1)
  expect_null(stepped$accepted)
  expect_equal(stepped$achievable, 1e-20)
  # An increase that is representable and Armijo-sufficient is taken, at
  # whatever length it first holds.
  real <- function(x) if (x[1L] > 0) -2705 + 1e-3 else -2705
  taken <- ctsem:::.ctBackendDampedStep(real, at = 0, step = 1, value = -2705,
    directional = 1)
  expect_false(is.null(taken$accepted))
  expect_equal(taken$accepted$value, -2705 + 1e-3)
})

test_that("the ladder terminates on an objective that never improves", {
  # No rung count bounds this loop, so the floor has to. An objective that is
  # flat everywhere must still return, and in a bounded number of calls.
  # Flat *at the value*: a first draft of this returned 0 against a baseline of
  # -2705, which is an improvement of 2705 and was accepted on the first rung --
  # the test failing rather than passing for the wrong reason.
  called <- 0
  flat <- function(x) { called <<- called + 1; -2705 }
  stepped <- ctsem:::.ctBackendDampedStep(flat, at = 0, step = 1, value = -2705,
    directional = 0.0436)
  expect_null(stepped$accepted)
  # log2(alpha * directional / (|value| * eps)) rungs, which is about 50 here
  # and can never be unbounded: each halving doubles the distance to the floor.
  expect_lt(called, 80)
  expect_gt(called, 10)
})

test_that("the resumed gradient tolerance is the one that would close the gap", {
  # `gap <= n |g|_inf^2 / (2 lambda_min)` over the trusted subspace, so a
  # gradient at the derived tolerance cannot leave a gap above the target. Both
  # directions asserted: the derivation, and the bound it claims.
  tol <- 0.01; npar <- 24L; lambda_min <- 6.573e-06
  derived <- ctsem:::.ctBackendGapGradientTolerance(lambda_min, npar, tol)
  expect_equal(derived, sqrt(2 * tol * lambda_min / npar))
  expect_lte(npar * derived^2 / (2 * lambda_min), tol * (1 + 1e-12))

  # A worse-conditioned model needs a tighter gradient for the same gap, which
  # is the whole reason the tolerance cannot be a constant: this is the measured
  # curvature of a real fit, nine orders below its own largest.
  expect_lt(derived, ctsem:::.ctBackendGapGradientTolerance(1, npar, tol))

  # Nothing to derive it from is NA, not a number: a model with no trusted
  # curvature needs a different report, not a tighter tolerance.
  expect_true(is.na(ctsem:::.ctBackendGapGradientTolerance(0, npar, tol)))
  expect_true(is.na(ctsem:::.ctBackendGapGradientTolerance(NA_real_, npar, tol)))
  expect_true(is.na(ctsem:::.ctBackendGapGradientTolerance(1, 0L, tol)))
})

# --- the reasoning behind the two tolerances, and what verifies it ----------
#
# The choices here are not measurements. Which tolerance to stop an optimiser
# on is a question about every model anyone might fit, and a number read off
# one of them would be fitted to that one. So the reasoning is written down in
# R/ctBackendOptimGap.R and these are the tests that it holds.

test_that("the optimiser's rule and the certification estimate the same thing", {
  # Why there is one tolerance and not two. The line search hands back
  # `dphi0 = g'p` with `p = -Bg`, so `-dphi0/2` is `g'Bg/2` -- L-BFGS's
  # estimate of the objective still available. With `B` exact that is the
  # certification's own quantity, to the last bit. A separate tolerance for the
  # optimiser would be a second answer to one question.
  H <- matrix(c(4, 1, 1, 25), 2, 2)
  g <- c(2, 5)
  exact <- ctsem:::.ctBackendOptimGap(-H, g)
  proxy <- 0.5 * as.numeric(t(g) %*% solve(H) %*% g)   # B = H^-1
  expect_equal(proxy, exact$gap)
  # And `dphi0` is its negative for the minimised objective, which is the sign
  # the engine reads: predicted gain is `abs(dphi0)/2`.
  dphi0 <- -as.numeric(t(g) %*% solve(H) %*% g)
  expect_equal(abs(dphi0) / 2, exact$gap)
})

test_that("only the exact quantity is invariant, which is why only it certifies", {
  # The architectural split, as a property rather than as a claim. Rescale a
  # coordinate: the exact gap is unchanged, and an L-BFGS approximation carried
  # over from before the rescaling is not. So the proxy can stop an optimiser
  # and cannot certify a fit.
  H <- matrix(c(4, 1, 1, 25), 2, 2)
  g <- c(2, 5)
  A <- diag(c(1, 1000))
  Hs <- solve(A) %*% H %*% solve(A)
  gs <- as.numeric(solve(A) %*% g)
  expect_equal(ctsem:::.ctBackendOptimGap(-Hs, gs)$gap,
    ctsem:::.ctBackendOptimGap(-H, g)$gap)
  # `B` from the old coordinates, applied in the new ones, as a limited-memory
  # approximation built before a rescaling would be.
  stale <- 0.5 * as.numeric(t(gs) %*% solve(H) %*% gs)
  expect_false(isTRUE(all.equal(stale,
    ctsem:::.ctBackendOptimGap(-H, g)$gap)))
})

test_that("lambda is a displacement in standard errors, which is what sets the scale", {
  # One dimension, where `lambda` literally is the displacement in standard
  # errors: this pins the arithmetic the reasoning rests on, whatever the bar is
  # later set to.
  information <- matrix(4, 1, 1)          # se = 1/2
  g <- 0.1 * 4 * 0.5                      # a gradient 0.1 se from the optimum
  out <- ctsem:::.ctBackendOptimGap(-information, g)
  expect_equal(out$lambda, 0.1)
  expect_equal(out$gap, 0.005)
})

test_that("the bar is set by the diagnostics, not by what an inference needs", {
  # An inference needs `lambda` around 0.1, which would be a gap of 0.005. The
  # curvature-based reports need far more: identifiability classifies
  # eigenvalues against a relative threshold, and one near it moves with the
  # estimate -- measured on a degenerate fixture, the flat direction is missed
  # at a gap of 3e-04 and found at 4e-07. So the bar is the tighter of the two
  # requirements, and this asserts which one won.
  expect_equal(ctsem:::.ctBackendGapTolerance(list()), 1e-6)
  expect_lt(ctsem:::.ctBackendGapTolerance(list()), 0.005)
})

test_that("the optimiser aims inside the bar, so a correction stays exceptional", {
  # If the optimiser targeted the bar itself, every fit whose limited-memory
  # proxy is slightly optimistic would fail the check and take a correction --
  # and a correction is a Hessian plus a resumed optimisation, which is the
  # expensive way to gain the last of the objective: 5197 objective calls
  # against 2729 for simply optimising further, on the same fixture. So the
  # inner target sits two orders inside.
  expect_lt(ctsem:::.ctBackendInnerGapTol(list()),
    ctsem:::.ctBackendGapTolerance(list()))
  expect_equal(ctsem:::.ctBackendInnerGapTol(list()), 1e-8)
  # And it follows the bar when the bar moves, rather than being a second
  # number to keep in step by hand.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(gaptol = 1e-4)), 1e-6)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(gaptol = 0.02)), 2e-4)
})

test_that("a loose stop is only used where something will check it", {
  # The coupling that makes the proxy safe: it can stop but not certify, so
  # with no certification running there is nothing to catch a stop that was too
  # early. Off in each of the three cases where the check does not run, and the
  # certification's own tolerance otherwise -- not a second number.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list()), 1e-8)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(estonly = TRUE)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(certify = FALSE)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(), intoverstates = FALSE), 0)
  # An explicit value wins, including zero and including where the default
  # would be off: the licence is a rule about what to do unasked, and a caller
  # who names the argument has asked. Accepting it and ignoring it would be the
  # worst of the three.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = 0.5)), 0.5)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = 0)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(
    list(innergaptol = 0.5, certify = FALSE)), 0.5)
  expect_equal(ctsem:::.ctBackendInnerGapTol(
    list(innergaptol = 0.5, estonly = TRUE)), 0.5)
  # And a nonsense value is off rather than an error: this decides how hard an
  # optimiser works, and refusing to fit over it would be the wrong trade.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(gaptol = -1)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = NA)), 0)
})

test_that("stopping is not certifying: the optimiser's verdict is not consulted", {
  # `.ctBackendCertify()` takes no convergence flag, by construction -- the
  # question it answers is about the curvature and not about how the optimiser
  # felt when it stopped. Asserted because the previous rule conflated the two,
  # and a future edit passing `converged` in here would undo the whole split.
  expect_false("converged" %in% names(formals(ctsem:::.ctBackendCertify)))
  short <- ctsem:::.ctBackendOptimGap(-diag(c(4, 25)), c(2, 5))
  expect_equal(ctsem:::.ctBackendCertify(short, probe = list(gain = 0),
    tolerance = 1e-6)$status, "suboptimal")
})
