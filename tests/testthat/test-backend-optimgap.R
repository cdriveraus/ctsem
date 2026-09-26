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
  expect_equal(verdict$status, "notstationary")
  expect_match(verdict$reason, "a flat direction still gains 1.4", fixed = TRUE)
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
    tolerance = 0.01, saturated = TRUE)$status, "saturated")
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

test_that("a flat direction that still gains says what the probe measured", {
  # AnomAuth S1 from default starts, in miniature: a flat direction along one
  # population sd whose probe gains 5.5e-05 at a quarter of a unit and less
  # further out. The verdict stays `notstationary` -- a small gain that turned
  # within the probe did not mean a maximum on S2, which had 1.6 exact nats
  # further out -- and what changes is that the reader is told what was seen:
  # which parameters, how much, at what length, and how far the probe looked.
  # The probe itself is the engine's now (`_ctsem_flat_probe`, tested in
  # test_ctsem_backend.jl); these are its numbers for that bump.
  probe <- list(gain = 5.5e-5, length = 0.25, longest = 4,
    direction = c(0, 1), ok = TRUE)
  gap <- ctsem:::.ctBackendOptimGap(-diag(c(4, 0)), c(1e-9, 3))
  verdict <- ctsem:::.ctBackendCertify(gap, probe, tolerance = 1e-6,
    parnames = c("drift", "popsd_cint"))
  expect_equal(verdict$status, "notstationary")
  expect_identical(verdict$reason, paste0("a flat direction (popsd_cint) ",
    "still gains 5.5e-05 within 0.25 raw units of the estimate; the probe ",
    "looked no further than 4"))

  # And the warning a fit raises from it, which names what would settle it.
  record <- ctsem:::.ctBackendProbeRecord(probe, c("drift", "popsd_cint"))
  expect_identical(record$residual_parameters, "popsd_cint")
  expect_equal(record$residual_longest, 4)
  fit <- list(optim = list(converged = FALSE),
    uncertainty = list(certification = c(verdict, record)))
  expect_warning(ctsem:::.ctBackendCertifyWarn(fit), paste0("Not converged: ",
    "a flat direction (popsd_cint) still gains 5.5e-05 within 0.25 raw units ",
    "of the estimate; the probe looked no further than 4. More iterations, ",
    "other starts, or ctFitProfile() on popsd_cint would say whether it keeps ",
    "rising."), fixed = TRUE)
})

test_that("a direction is described by its largest loadings, in one rule", {
  # A share of the largest loading, not an absolute bar: spread over ten
  # coordinates a unit vector has loadings near 0.32, over twenty near 0.22,
  # and an absolute 0.25 named all of the first and none of the second.
  even <- rep(1 / sqrt(20), 20)
  expect_length(ctsem:::.ctBackendLoadedCoordinates(even), 20L)
  # Largest first, and a coordinate well under a third of the largest is left
  # out.
  expect_identical(ctsem:::.ctBackendLoadedCoordinates(c(0.1, -0.9, 0.4, 0.01)),
    c(2L, 3L))
  expect_identical(ctsem:::.ctBackendLoadedCoordinates(c(0.5, 0.5, 0.5, 0.5),
    most = 2L), c(1L, 2L))
  expect_length(ctsem:::.ctBackendLoadedCoordinates(numeric(3)), 0L)
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

# --- where a Hessian may be reused -------------------------------------------
#
# The engine's finish keeps the Hessian it took at the hand-over when its steps
# moved the estimate less than a hundredth of a standard error, and says where
# it was taken. Everything that reuses a stored Hessian asks the same question
# with the same arithmetic, which is what these pin.

test_that("distance is measured in the standard errors the Hessian implies", {
  # One dimension: information 4, standard error a half, so a move of 0.005 is
  # a hundredth of one.
  expect_equal(ctsem:::.ctBackendHessianDistance(matrix(-4, 1, 1), 0, 0.005),
    0.01)
  expect_equal(ctsem:::.ctBackendHessianDistance(matrix(-4, 1, 1), 0.3, 0.3), 0)
  # The largest over the coordinates, each against its own standard error: a
  # correlated information, whose inverse is the covariance a fit reports.
  info <- matrix(c(4, 1, 1, 25), 2, 2)
  se <- sqrt(diag(solve(info)))
  d <- c(0.01, -0.02)
  expect_equal(ctsem:::.ctBackendHessianDistance(-info, c(0, 0), d),
    max(abs(d) / se))
  # A coordinate the trusted curvature says nothing about has no standard
  # error there, so moving it at all is infinitely far -- and not moving it
  # costs nothing.
  flat <- -diag(c(4, 0))
  expect_equal(ctsem:::.ctBackendHessianDistance(flat, c(0, 0), c(0, 1e-9)), Inf)
  expect_equal(ctsem:::.ctBackendHessianDistance(flat, c(0, 0), c(0.005, 0)),
    0.01)
  # A matrix that cannot be decomposed describes nothing.
  expect_equal(ctsem:::.ctBackendHessianDistance(matrix(NA_real_, 1, 1), 0, 1),
    Inf)
  expect_equal(ctsem:::.ctBackendHessianReuse(), 0.01)
})

test_that("a stored Hessian is reused within a hundredth of a standard error, and no further", {
  fit <- list(uncertainty = list(hessian = matrix(-4, 1, 1), evaluated_at = 1))
  stored <- fit$uncertainty$hessian
  expect_identical(ctsem:::.ctBackendStoredHessian(fit, 1), stored)
  expect_identical(ctsem:::.ctBackendStoredHessian(fit, 1.004), stored)
  expect_null(ctsem:::.ctBackendStoredHessian(fit, 1.006))
  # A posterior mean a tenth of a standard error from the Laplace point is
  # somewhere else, which is the case `evaluated_at` exists for.
  expect_null(ctsem:::.ctBackendStoredHessian(fit, 1.05))
  # A Hessian that does not say where it was taken predates the field and is
  # not reused; nor is one of the wrong size.
  expect_null(ctsem:::.ctBackendStoredHessian(
    list(uncertainty = list(hessian = stored)), 1))
  expect_null(ctsem:::.ctBackendStoredHessian(fit, c(1, 1)))
})

# --- the resume rule ----------------------------------------------------------
#
# `.ctBackendCorrectResult()` certifies what the engine's finish handed back and
# decides one thing: whether to resume. Fake results and a fake optimiser, so
# the rule is tested without an engine.

.fake_result <- function(hessian, gradient, value, x = c(0, 0), ...) {
  utils::modifyList(list(minimizer = x, gradient = gradient,
    maximum_loglik = value, hessian = hessian, hessian_evaluated_at = x,
    hessian_distance = 0, probe_ran = FALSE, probe_gain = 0, probe_length = 0,
    probe_longest = 0, probe_direction = 0, newton_hessians = 1L,
    iterations = 10L, f_calls = 12L, g_calls = 11L,
    trace = list(objective = c(value - 5, value))), list(...))
}

test_that("a fit short of its optimum is resumed under its own rules, with its progress carried", {
  # A gap of one nat: the finish did not close it (a line search that ran out,
  # say), so the certification says suboptimal.
  short <- .fake_result(-diag(c(4, 25)), c(2, 5), value = -100)
  done <- .fake_result(-diag(c(4, 25)), c(0, 0), value = -99, x = c(0.5, 0.2))
  seen <- list()
  optimise <- function(from, carried) {
    seen[[length(seen) + 1L]] <<- list(from = from, carried = carried)
    done
  }
  out <- ctsem:::.ctBackendCorrectResult(short, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 2L, optimise = optimise)
  # Once, from where it stopped, carrying the five nats the fit had made -- and
  # nothing else: no tolerance tightened, no cap raised, no rule switched off.
  expect_length(seen, 1L)
  expect_equal(seen[[1L]]$from, c(0, 0))
  expect_equal(seen[[1L]]$carried, 5)
  expect_identical(out$result, done)
  expect_equal(out$certification$status, "certified")
  expect_length(out$corrections, 1L)
  expect_true(out$corrections[[1L]]$resumed)
  expect_equal(out$corrections[[1L]]$status, "suboptimal")
  # Work over both stages, and the Hessians both finishes formed.
  expect_equal(out$hessians, 2L)
  expect_equal(unname(out$totals[["iterations"]]), 20)
  # The matrix the fit keeps is the one the last certification decided on,
  # with where it was evaluated.
  expect_identical(out$hessian, done$hessian)
  expect_equal(out$evaluated_at, c(0.5, 0.2))
})

test_that("a run the fit does not keep is still counted in its work", {
  # The stage a stall escape replaced, an escape that did not come out ahead,
  # the run before a substep refit: each did its work, and the fit's counts are
  # the whole of what it ran, as the corrections' totals are. They were the kept
  # run's alone; the optimiser bench's tally of engine runs is what showed it.
  kept <- list(iterations = 40L, f_calls = 90L, g_calls = 60L,
    newton_steps = 2L, newton_hessians = 1L)
  dropped <- list(iterations = 75L, f_calls = 200L, g_calls = 150L,
    newton_steps = 3L, newton_hessians = 2L, newton_subset_hessians = 0L)
  out <- ctsem:::.ctJuliaAddRunCounts(kept, dropped)
  expect_equal(out$iterations, 115)
  expect_equal(out$f_calls, 290)
  expect_equal(out$g_calls, 210)
  expect_equal(out$newton_steps, 5)
  expect_equal(out$newton_hessians, 3)
  # The kept run's own count stays readable, since its trace is the one the
  # fit reports; a second dropped run adds to the totals and leaves it alone.
  expect_equal(out$stage_iterations, 40L)
  again <- ctsem:::.ctJuliaAddRunCounts(out,
    ctsem:::.ctJuliaRunCounts(dropped))
  expect_equal(again$iterations, 190)
  expect_equal(again$stage_iterations, 40L)
  # Nothing dropped, nothing changed.
  expect_identical(ctsem:::.ctJuliaAddRunCounts(kept, list()), kept)
})

test_that("a resume that does not improve is not taken", {
  short <- .fake_result(-diag(c(4, 25)), c(2, 5), value = -100)
  worse <- .fake_result(-diag(c(4, 25)), c(0, 0), value = -101)
  out <- ctsem:::.ctBackendCorrectResult(short, spec = list(), npar = 2L,
    tolerance = 1e-6, optimise = function(from, carried) worse)
  expect_identical(out$result, short)
  expect_equal(out$certification$status, "suboptimal")
  expect_false(out$corrections[[1L]]$resumed)
})

test_that("rounds are capped, and a certified fit takes none", {
  short <- .fake_result(-diag(c(4, 25)), c(2, 5), value = -100)
  calls <- 0L
  again <- function(from, carried) {
    calls <<- calls + 1L
    .fake_result(-diag(c(4, 25)), c(2, 5), value = -100 + calls)
  }
  out <- ctsem:::.ctBackendCorrectResult(short, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 2L, optimise = again)
  expect_equal(calls, 2L)
  expect_length(out$corrections, 2L)
  calls <- 0L
  none <- ctsem:::.ctBackendCorrectResult(short, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 0L, optimise = again)
  expect_equal(calls, 0L)
  expect_equal(none$certification$status, "suboptimal")
  fine <- .fake_result(-diag(c(4, 25)), c(0, 0), value = -100)
  kept <- ctsem:::.ctBackendCorrectResult(fine, spec = list(), npar = 2L,
    tolerance = 1e-6, optimise = again)
  expect_equal(calls, 0L)
  expect_equal(kept$certification$status, "certified")
  expect_length(kept$corrections, 0L)
})

test_that("a flat direction that still gains stops, and so does a saddle the finish could not leave", {
  never <- function(from, carried) stop("resumed")
  # `notstationary` stops: continuing along the probe's direction walked
  # AnomAuth S1 from its good optimum into the spurious basin.
  flatgain <- .fake_result(-diag(c(4, 0)), c(0, 3), value = -100,
    probe_ran = TRUE, probe_gain = 0.5, probe_length = 1, probe_longest = 4,
    probe_direction = c(0, 1))
  out <- ctsem:::.ctBackendCorrectResult(flatgain, spec = list(), npar = 2L,
    tolerance = 1e-6, optimise = never)
  expect_equal(out$certification$status, "notstationary")
  expect_equal(out$certification$residual_gain, 0.5)
  expect_length(out$corrections, 0L)
  # A saddle whose negative curvature the finish already tried, and found
  # nothing along: the optimiser resumed from there would hand the same point
  # back, and on a Laplace fit each such round is a Hessian.
  saddle <- .fake_result(diag(c(1, -1)), c(0, 0), value = -100,
    newton_saddle = TRUE, newton_ladder_tried = TRUE)
  stuck <- ctsem:::.ctBackendCorrectResult(saddle, spec = list(), npar = 2L,
    tolerance = 1e-6, optimise = never)
  expect_equal(stuck$certification$status, "notmaximum")
  expect_length(stuck$corrections, 0L)
  # Without that record the saddle is resumed, as any point short of a maximum.
  resumed <- 0L
  once <- function(from, carried) {
    resumed <<- resumed + 1L
    .fake_result(-diag(c(1, 1)), c(0, 0), value = -99)
  }
  saddle$newton_ladder_tried <- FALSE
  ctsem:::.ctBackendCorrectResult(saddle, spec = list(), npar = 2L,
    tolerance = 1e-6, optimise = once)
  expect_equal(resumed, 1L)
  # And a tried saddle with the gap in its trusted directions still open is
  # resumed too: it is short of its optimum, whatever the negative curvature
  # says. Gated-gaps config B8 was stopped 10.8 nats short without this.
  open <- .fake_result(diag(c(1, -1)), c(0, 2), value = -100,
    newton_saddle = TRUE, newton_ladder_tried = TRUE)
  resumed <- 0L
  out <- ctsem:::.ctBackendCorrectResult(open, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 1L, optimise = once)
  expect_equal(resumed, 1L)
  expect_equal(out$corrections[[1L]]$status, "notmaximum")
})

test_that("the probe the finish ran is what the certification reads", {
  # The engine's numbers go straight into the verdict: no probe of R's own.
  gaining <- .fake_result(-diag(c(4, 0)), c(1e-9, 3), value = -100,
    probe_ran = TRUE, probe_gain = 5.5e-5, probe_length = 0.25,
    probe_longest = 4, probe_direction = c(0, 1))
  out <- ctsem:::.ctBackendCorrectResult(gaining, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 0L)
  expect_equal(out$certification$status, "notstationary")
  expect_equal(out$certification$residual_length, 0.25)
  expect_equal(out$certification$residual_longest, 4)
  # A probe that did not run is no gain measured.
  quiet <- .fake_result(-diag(c(4, 0)), c(1e-9, 3), value = -100)
  expect_equal(ctsem:::.ctBackendCorrectResult(quiet, spec = list(), npar = 2L,
    tolerance = 1e-6, maxtries = 0L)$certification$status, "certified")
  # And the engine's "none" sentinels are not read as a direction.
  expect_null(ctsem:::.ctBackendProbeFields(list(probe_ran = TRUE,
    probe_gain = 1, probe_direction = 0), 2L))
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

test_that("the bar is read from wherever the fit keeps its controls", {
  # `ctFit()` stores them under `$args$resolved` and `$args$input`. The
  # backend's own `$args`, which a fit carries while it is being built, and a
  # stored fit from before that split keep them at the top -- the only place
  # this used to look, so a finished fit was re-certified at the default.
  asked <- function(tol) list(optimcontrol = list(gaptol = tol))
  tolerance <- function(args) ctsem:::.ctBackendGapTolerance(list(args = args))
  expect_equal(tolerance(list(input = asked(1e-4), resolved = asked(1e-3))), 1e-3)
  expect_equal(tolerance(list(input = asked(1e-4))), 1e-4)
  expect_equal(tolerance(asked(1e-5)), 1e-5)
  expect_equal(tolerance(list(input = list(optimcontrol = list()),
    resolved = list(optimcontrol = list()))), 1e-6)
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
  # early. Off in each of the two cases where the check does not run, and the
  # certification's own tolerance otherwise -- not a second number.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list()), 1e-8)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(certify = FALSE)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(), intoverstates = FALSE), 0)
  # `estonly` is not one of them, and this asserts that rather than assuming
  # it. It asks for the estimate without the uncertainty and correction phases,
  # which is a statement about what happens *after* the optimisation -- and
  # while it was in the list the same model fitted with and without it ran
  # different stopping rules and could stop in different places. Every other
  # reading of `estonly` in the package skips a post-fit step.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(estonly = TRUE)), 1e-8)
  # An explicit value wins, including zero and including where the default
  # would be off: the licence is a rule about what to do unasked, and a caller
  # who names the argument has asked. Accepting it and ignoring it would be the
  # worst of the three.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = 0.5)), 0.5)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = 0)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(
    list(innergaptol = 0.5, certify = FALSE)), 0.5)
  # And a nonsense value is off rather than an error: this decides how hard an
  # optimiser works, and refusing to fit over it would be the wrong trade.
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(gaptol = -1)), 0)
  expect_equal(ctsem:::.ctBackendInnerGapTol(list(innergaptol = NA)), 0)
})

# Where a stage that stopped short is resumed from, decided without a fit.
#
# Two routes to a point and they are not interchangeable. The in-flight stall
# check finds one while the optimiser is still running, and can only fire on a
# run that presents a stalled window. A fit that climbs steadily and then stops
# dead -- a line search that finds nothing, which is how most short fits end --
# never presents one, and for that case the post-fit probe has already found an
# improving point by the time the result is assembled. Both were measured on
# the same model from two starting values inside a flat transform: the first
# recovered 387 nats, the second 169, and neither covered the other's case.
test_that("a stage that stopped short is resumed from the better point", {
  base <- list(minimizer = c(1, 2, 3), stall_point = c(0, 0, 0),
    overshoot_point = c(0, 0, 0), overshoot_gain = 5,
    overshoot_parameters = 1L, stall_parameters = 1L)

  # The in-flight route: the engine stopped here and brought the point with it.
  stalled <- modifyList(base, list(stopped_by_stall = TRUE,
    stall_point = c(0.5, 2, 3)))
  escaped <- ctsem:::.ctBackendStallEscape(stalled, list(), NULL)
  expect_equal(as.numeric(escaped), c(0.5, 2, 3))
  # Which coordinates were moved rides along on the point, because the caller
  # that pins them needs to know and every other caller does not.
  expect_equal(attr(escaped, "coordinates"), 1L)

  # The post-fit route: it stopped for its own reasons, and its own probe says
  # the point it stopped at is not a maximum.
  overshot <- modifyList(base, list(stopped_by_stall = FALSE, overshot = TRUE,
    overshoot_point = c(0.25, 2, 3)))
  escaped <- ctsem:::.ctBackendStallEscape(overshot, list(), NULL)
  expect_equal(as.numeric(escaped), c(0.25, 2, 3))
  expect_equal(attr(escaped, "coordinates"), 1L)
  # The engine's `0` sentinel for "none" is dropped rather than reaching a pin
  # as a coordinate index.
  none <- ctsem:::.ctBackendEscapeCoordinates(c(1, 2, 3), 0L, 3L)
  expect_length(attr(none, "coordinates"), 0L)

  # A fit with nothing wrong with it is left where it is. This is the case that
  # runs on almost every fit in the package, so it is the one that has to be
  # free.
  fine <- modifyList(base, list(stopped_by_stall = FALSE, overshot = FALSE))
  expect_null(ctsem:::.ctBackendStallEscape(fine, list(), NULL))

  # `[0]` and `[0.0]` are the engine's "none" sentinels, because a zero-length
  # vector deadlocks the R bridge. A one-element point is not a point for a
  # three-parameter model and must not be read as one.
  sentinel <- modifyList(base, list(stopped_by_stall = TRUE,
    stall_point = 0, overshot = TRUE, overshoot_point = 0))
  expect_null(ctsem:::.ctBackendStallEscape(sentinel, list(), NULL))

  # Nor is a point with a non-finite entry, which is what a probe that ran off
  # the edge of the objective hands back.
  broken <- modifyList(base, list(stopped_by_stall = FALSE, overshot = TRUE,
    overshoot_point = c(0.25, NaN, 3)))
  expect_null(ctsem:::.ctBackendStallEscape(broken, list(), NULL))

  # And nothing escapes on the state-explicit route, whatever the probe found.
  # The joint mode is degenerate -- the innovations re-optimise to absorb
  # almost any parameter change -- so a pullback can nearly always find
  # something and "not a maximum" stops carrying information.
  expect_null(ctsem:::.ctBackendStallEscape(overshot, list(), NULL,
    escapes = FALSE))
  expect_null(ctsem:::.ctBackendStallEscape(stalled, list(), NULL,
    escapes = FALSE))
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


# --- the fit's own verdict ---------------------------------------------------

.verdict_fit <- function(status, certified = FALSE, pending = TRUE) {
  list(optim = list(converged = FALSE, convergence_pending = pending),
    uncertainty = list(certification = list(status = status,
      certified = certified, reason = "because")))
}

test_that("the curvature's verdict replaces the optimiser's, in both directions", {
  # The complaint this answers: a fit could be certified -- the optimum bounded
  # within `gaptol` of the estimate -- and still report `converged = FALSE`,
  # because `converged` was keyed on a gradient bar computed before anything
  # knew the curvature. Whichever way they disagree, the measurement wins.
  certified <- ctsem:::.ctBackendCertifiedVerdict(
    .verdict_fit("certified", certified = TRUE))
  expect_true(certified$optim$converged)
  # And the held complaint is dropped rather than left to be warned about.
  expect_null(certified$optim$convergence_pending)

  # The other direction: the optimiser was happy, the curvature is not.
  happy <- .verdict_fit("suboptimal", pending = FALSE)
  happy$optim$converged <- TRUE
  expect_false(ctsem:::.ctBackendCertifiedVerdict(happy)$optim$converged)
})

test_that("converged says maximum, and the two findings that are not failures", {
  # `saturated` is a maximum with a coordinate the data does not determine,
  # which is a result and not a failure -- reporting it as a failure to
  # converge is what put `converged = FALSE` on 45 of 64 fits whose log
  # likelihoods matched stan's to the digit. `notstationary` is the other half
  # of what used to share that name, and it *is* a failure: stepping along the
  # direction gains likelihood, so the point is not a maximum at all.
  # `unidentified` is `saturated`'s name before 2026-09-25, and a stored fit
  # still carries it.
  expected <- c(certified = TRUE, saturated = TRUE, unidentified = TRUE,
    suboptimal = FALSE, notstationary = FALSE, notmaximum = FALSE,
    unknown = FALSE)
  got <- vapply(names(expected), function(status)
    isTRUE(ctsem:::.ctBackendCertifiedVerdict(
      .verdict_fit(status, certified = identical(status, "certified"))
    )$optim$converged), logical(1))
  expect_equal(got, expected)
})

test_that("a stored fit's 'unidentified' reads as 'saturated' everywhere", {
  status <- ctsem:::.ctBackendCertificationStatus
  expect_identical(status(list(status = "unidentified")), "saturated")
  expect_identical(status(list(status = "saturated")), "saturated")
  expect_identical(status(list(status = "certified")), "certified")
  expect_identical(status(NULL), character())
  expect_identical(status(list(status = character())), character())
  # And the warning a stored fit raises is the finding's, not a failure's.
  old <- .verdict_fit("unidentified")
  old$uncertainty$certification$reason <- "a parameter transform has saturated"
  expect_warning(ctsem:::.ctBackendCertifyWarn(old),
    "This fit is a maximum, but a parameter transform has saturated.",
    fixed = TRUE)
})

test_that("a Laplace spec is differentiated by ctsem_hessian whatever gradient says", {
  # `ctsem_hessian_forward` has no method for the Laplace objective, and on that
  # route `gradient` chooses nothing about differentiation. See
  # test-julia-laplace.R for what asking for it used to cost a fit.
  pick <- ctsem:::.ctBackendHessianFunction
  expect_identical(pick(list(laplace = list(levels = list())), "forward"),
    "ctsem_hessian")
  expect_identical(pick(list(laplace = list(levels = list())), "adjoint"),
    "ctsem_hessian")
  expect_identical(pick(list(), "forward"), "ctsem_hessian_forward")
  expect_identical(pick(list(), "adjoint"), "ctsem_hessian")
})

test_that("a fit that certified nothing keeps the optimiser's verdict", {
  # `estonly`, or `certify = FALSE`: there is no better measurement, so there is
  # nothing to replace it with, and inventing one would be worse than the bit
  # the optimiser can honestly supply. `convergence_pending` is what says so.
  bare <- list(optim = list(converged = TRUE, convergence_pending = TRUE))
  kept <- ctsem:::.ctBackendCertifiedVerdict(bare)
  expect_true(kept$optim$converged)
  expect_true(kept$optim$convergence_pending)
  # An empty certification is the same case, not a verdict of FALSE.
  empty <- list(optim = list(converged = TRUE),
    uncertainty = list(certification = list(status = character(0))))
  expect_true(ctsem:::.ctBackendCertifiedVerdict(empty)$optim$converged)
})
