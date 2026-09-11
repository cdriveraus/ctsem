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
