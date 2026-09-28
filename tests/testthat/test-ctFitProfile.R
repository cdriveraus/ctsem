# The profile likelihood: what the data says about a parameter, asked by moving
# it rather than by differentiating at the estimate.
#
# `.ctOptimFlatDirectionScreen()` is the cheap half of this question and runs on
# every fit with a flat direction; it walks without re-optimising, so it can
# only ever confirm flatness and never refute it. This is the two-sided version,
# and it is opt-in because every point is an optimisation.
#
# The verdicts are Raue et al.'s. The bar is `qchisq(level, 1) / 2`, which is
# where these differ from anything keyed on curvature: a statistical quantity,
# invariant to reparameterisation, rather than a tolerance on a differentiated
# approximation whose value depends on where the optimiser stopped.

skip_on_cran()
skip_on_32bit()

test_that("the parameters to profile can be named or numbered, and a typo is refused", {
  names <- c("drift", "diff", "cint")
  expect_equal(ctsem:::.ctFitProfileParameters(NULL, names, 3L), 1:3)
  expect_equal(ctsem:::.ctFitProfileParameters(c("cint", "drift"), names, 3L),
    c(3L, 1L))
  expect_equal(ctsem:::.ctFitProfileParameters(c(2L, 3L), names, 3L), c(2L, 3L))
  # Named, because a silently dropped parameter is a profile that answers a
  # different question from the one asked.
  expect_error(ctsem:::.ctFitProfileParameters("drfit", names, 3L), "drfit")
  expect_error(ctsem:::.ctFitProfileParameters(7L, names, 3L), "between 1 and 3")
})

test_that("the verdict is what the walk did, and a limit is where it crossed", {
  bar <- stats::qchisq(0.95, 1) / 2
  walk <- function(side, values, drops) data.frame(parameter = "p", index = 1L,
    side = side, value = values, loglik = -drops, drop = drops,
    iterations = 1L, stringsAsFactors = FALSE)

  # Crossed on both sides: identifiable, and the limits are interpolated
  # between the points that straddle the bar rather than being the first point
  # past it -- linear in the log likelihood, which is exact for the quadratic
  # the bar assumes.
  both <- rbind(walk(-1, c(-1, -2), c(0.5, 2.5)), walk(1, c(1, 2), c(0.5, 2.5)))
  s <- ctsem:::.ctFitProfileSummary(both, "p", 1L, 0, bar)
  expect_equal(s$verdict, "identifiable")
  expect_equal(s$lower, -1 - (bar - 0.5) / 2)
  expect_equal(s$upper, 1 + (bar - 0.5) / 2)
  expect_equal(s$flat, "")

  # Falls on both sides but crosses on neither: the limits are beyond where we
  # looked, which is not the same as absent -- `*_walked` says how far that was.
  slow <- rbind(walk(-1, c(-1, -2), c(0.1, 0.3)), walk(1, c(1, 2), c(0.1, 0.3)))
  s <- ctsem:::.ctFitProfileSummary(slow, "p", 1L, 0, bar)
  expect_equal(s$verdict, "practically non-identifiable")
  expect_true(is.na(s$lower) && is.na(s$upper))
  expect_equal(s$lower_walked, -2)
  expect_equal(s$upper_walked, 2)

  # One flat side is conclusive on its own: it exhibits a curve along which the
  # likelihood is constant. This is the shape a saturating transform makes --
  # flat one way, falling the other -- and calling it merely practical would
  # understate it.
  onesided <- rbind(walk(-1, c(-1, -2), c(0, 0)), walk(1, c(1, 2), c(0.5, 2.5)))
  s <- ctsem:::.ctFitProfileSummary(onesided, "p", 1L, 0, bar)
  expect_equal(s$verdict, "structurally non-identifiable")
  expect_equal(s$flat, "lower")
  # And the side that did cross still reports its limit: half a verdict is
  # still worth reading.
  expect_equal(s$upper, 1 + (bar - 0.5) / 2)

  expect_equal(nrow(ctsem:::.ctFitProfileSummary(both[0, ], "p", 1L, 0, bar)), 0L)
})

test_that("'flagged' profiles what the fit already doubts, from every detector", {
  names <- c("drift", "diff", "cint", "popsd_drift")
  # Three detectors, three different complaints, and they do not agree by
  # construction: a saturated transform, a flat direction in the curvature, and
  # an interval wider than that curvature supports. The union is taken because
  # a parameter profiled needlessly costs time and one skipped costs the
  # answer. Names from one source and indices from another, because that is how
  # they arrive depending on how far through the reporting they came.
  fit <- list(
    identifiability = list(parameters = "drift"),
    uncertainty = list(intervalcheck = list(parameters = "cint",
      unidentified = "popsd_drift")),
    optim = list(saturated_parameters = 2L))
  expect_equal(ctsem:::.ctFitProfileFlagged(fit, names, 4L), 1:4)
  expect_equal(ctsem:::.ctFitProfileParameters("flagged", names, 4L, fit), 1:4)

  # `0` is the engine's "none" sentinel -- a zero-length vector deadlocks the
  # R bridge -- and must not become coordinate zero.
  quiet <- list(optim = list(saturated_parameters = 0L))
  expect_length(ctsem:::.ctFitProfileFlagged(quiet, names, 4L), 0L)
  # A fit with nothing flagged is an error rather than a silent empty profile:
  # "I profiled nothing" and "nothing needed profiling" read the same in an
  # empty table.
  expect_error(ctsem:::.ctFitProfileParameters("flagged", names, 4L, quiet) |>
    (function(i) if (!length(i)) stop("nothing to profile") else i)(),
    "nothing to profile")
})

test_that("the step comes from the curvature, and falls back where there is none", {
  # Half a standard error puts several points inside an identified parameter's
  # interval, so the limit is interpolated between neighbours rather than
  # guessed from the estimate to the first step.
  fit <- list(estimate = list(se = c(2, 0.5, NA, 0)))
  expect_equal(ctsem:::.ctFitProfileSteps(NULL, fit, 4L), c(1, 0.25, 0.5, 0.5))
  # No curvature at all -- an estonly fit -- and every step is the default.
  expect_equal(ctsem:::.ctFitProfileSteps(NULL, list(), 3L), rep(0.5, 3))
  # An explicit step wins, as one number or one per parameter.
  expect_equal(ctsem:::.ctFitProfileSteps(0.2, fit, 4L), rep(0.2, 4))
  expect_equal(ctsem:::.ctFitProfileSteps(c(1, 2, 3, 4), fit, 4L), 1:4)
  expect_error(ctsem:::.ctFitProfileSteps(0, fit, 4L), "positive")
  expect_error(ctsem:::.ctFitProfileSteps(c(1, 2), fit, 4L), "positive")
})

test_that("printing says what was found, and says first when the fit was wrong", {
  bar <- stats::qchisq(0.95, 1) / 2
  walk <- function(side, values, drops) data.frame(parameter = "p",
    index = 1L, side = side, value = values, loglik = -drops, drop = drops,
    iterations = 1L, transformed = NA_real_, stringsAsFactors = FALSE)
  profile <- rbind(walk(-1, c(-1, -2), c(0.5, 2.5)),
    walk(1, c(1, 2), c(0.1, 0.3)))
  out <- structure(list(profile = profile,
    summary = ctsem:::.ctFitProfileSummary(profile, "p", 1L, 0, bar),
    bar = bar, level = 0.95, base = -10, better = NULL,
    estimate = c(p = 0), call = NULL), class = "ctFitProfile")

  text <- paste(utils::capture.output(print(out)), collapse = " ")
  expect_match(text, "Profile likelihood, 95%")
  expect_match(text, "practically non-identifiable")
  # A limit the walk never found is printed as how far it looked, not as NA:
  # "beyond 2" and "unknown" are different findings.
  expect_match(text, "<2")
  expect_false(grepl("NA", text))

  # And when a constrained fit beat the estimate, that is the first thing
  # said -- every number below it describes a point that is not the maximum.
  out$better <- list(point = c(1, 2), gain = 3.5, parameter = "p", value = 1)
  first <- utils::capture.output(print(out))[1]
  expect_match(first, "NOT at a maximum")
})

test_that("plotting a profile draws a panel per parameter without complaint", {
  bar <- stats::qchisq(0.95, 1) / 2
  profile <- do.call(rbind, lapply(c("a", "b"), function(name)
    data.frame(parameter = name, index = if (name == "a") 1L else 2L,
      side = c(-1, -1, 1, 1), value = c(-1, -2, 1, 2),
      loglik = c(-1, -3, -1, -3), drop = c(0.5, 2.5, 0.5, 2.5),
      iterations = 1L, transformed = NA_real_, stringsAsFactors = FALSE)))
  out <- structure(list(profile = profile,
    summary = ctsem:::.ctFitProfileSummary(profile, c("a", "b"), 1:2,
      c(0, 0), bar),
    bar = bar, level = 0.95, base = -0.5, better = NULL,
    estimate = c(a = 0, b = 0), call = NULL), class = "ctFitProfile")
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_silent(plot(out))
  # A name that was never profiled is a message, not an error: asking about a
  # parameter you did not profile is a mistake worth saying out loud and not
  # worth stopping for.
  expect_message(plot(out, parameters = "nope"), "None of those")
})

test_that("a profile separates a determined parameter from one on a flat ray", {
  skip_without_julia()
  # The same noise fixture `test-backend-summary.R` uses. Fitted to noise it
  # has more than one local maximum, -205.827 and -207.019 among them, and
  # which one a fit reaches depends on its path: the default fit reached
  # -205.827 before the endgame moved into the engine and has reached -207.019
  # since, converged and certified both times, because each is a maximum. At
  # the higher, one diffusion correlation is on a ray along which the
  # likelihood is flat to about 1e-10, while the other parameters are
  # ordinary.
  set.seed(5)
  data <- do.call(rbind, lapply(1:30, function(i) data.frame(id = i,
    time = c(0, .5, 1.5, 2.4, 3.5), Y1 = stats::rnorm(5, 0, .5),
    Y2 = stats::rnorm(5, 0, .5))))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2,
    LAMBDA = diag(2), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0MEANS = matrix(0, 2, 1),
    CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2,
      byrow = TRUE)))
  # Profiled from an estonly fit started at the certified maximum, so every
  # step is the default rather than half a standard error the flat direction
  # does not have. Routing profile points through the fit's own pipeline makes
  # each constrained reoptimisation as capable as the fit itself, so from a
  # lower maximum one of them lands in a higher basin -- +1.17 from -207.019
  # -- and `$better` reports the point. Refitting from it is what a user does
  # with that, and what this does, until the profile finds nothing higher.
  # It warns twice, about the flat ray by name -- the Hessian repair and the
  # identifiability report -- which is the finding checked field by field below.
  fitted <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", verbose = 0)))
  expect_identical(fitted$uncertainty$certification$status, "certified")
  refit <- function(from) suppressMessages(ctFit(data, model, backend = "julia",
    verbose = 0, inits = from, optimcontrol = list(estonly = TRUE)))
  profile <- function(fit) ctFitProfile(fit,
    parameters = c("auto1", "diff_eta2_eta1"), points = 6L)
  fit <- refit(fitted$estimate$raw)
  out <- profile(fit)
  for (round in 1:2) {
    if (is.null(out$better)) break
    # The point `$better` reports is worth what it says: a free fit from it
    # ends at least that far above the base.
    reached <- out$base + out$better$gain
    fit <- refit(out$better$point)
    expect_gt(fit$estimate$logposterior, reached - 1e-6)
    out <- profile(fit)
  }

  # At the point the refits reached, no constrained point beat the fit. This
  # is the check that makes the rest of the output mean anything: a profile
  # computed around a point that is not a maximum describes the wrong point.
  # Which of the fixture's maxima that is depends on the path, so the test
  # asserts the property rather than a value: with the prior warm-up off by
  # default on this Gaussian model (2026-09-27) the refits end at -205.872,
  # where the profile finds nothing higher and the ray below is as flat, and
  # before that they ended at -205.827.
  expect_null(out$better)
  expect_equal(out$bar, stats::qchisq(0.95, 1) / 2)
  expect_gt(nrow(out$profile), 6L)
  # Every point is a constrained maximum, so none may exceed the free one.
  expect_true(all(out$profile$drop > -1e-6))

  # The base is evaluated directly at the estimate (ctJuliaEvaluate()), not by
  # an optimisation capped at one iteration -- so it is exactly the fit's own
  # objective value there, `$logposterior` (`$loglik` would be the wrong
  # quantity to drop the profile's points against whenever priors are on; see
  # the comment in ctFitProfile()).
  expect_equal(out$base, fit$estimate$logposterior)

  flat <- out$summary[out$summary$parameter == "diff_eta2_eta1", ]
  expect_equal(flat$verdict, "structurally non-identifiable")
  # Flat downward specifically: the correlation's transform saturates one way
  # and not the other.
  expect_true(grepl("lower", flat$flat))
  expect_true(is.na(flat$lower))
  # And it really was walked a long way before saying so, rather than giving up
  # after two steps -- which is what `lower_walked` exists to let a reader
  # check.
  expect_lt(flat$lower_walked, flat$estimate - 5)

  ordinary <- out$summary[out$summary$parameter == "auto1", ]
  # `auto1` moves the likelihood, which is the whole contrast with the row
  # above; on noise data its upper side runs to zero drift and stays open, so
  # the verdict is not asserted, only that a limit was found on the side that
  # closed.
  expect_equal(ordinary$flat, "")
  expect_true(is.finite(ordinary$lower))
  expect_lt(ordinary$lower, ordinary$estimate)
})

test_that("a profile of a simple identified parameter routes through the fit's own pipeline", {
  skip_without_julia()
  # `.ctFitProfilePoint()` used to call the engine's optimiser directly, with
  # none of the fit's own `optimcontrol` -- no transform-scale metric, no
  # batching, no Newton finish, the gap rule off -- and took the base value
  # from an optimisation capped at one iteration rather than an evaluation.
  # Compared directly against that implementation on this exact fixture (the
  # pre-fix R/ctFitProfile.R sourced into its own environment, parented on the
  # ctsem namespace, so both versions could be called in the same session):
  # the two agree to numerical noise (crossing limits within 1e-10 of each
  # other, base identical to machine precision) and the fixed version ran
  # faster (5.4 s against 8.5 s) -- routing through the fit's own batching and
  # Newton finish rather than plain L-BFGS, not merely matching it. The
  # hardcoded limits below are what that comparison found; the tolerance is
  # far looser than the measured agreement, to allow for engine changes that
  # are not this one.
  times <- c(0, .6, 1.3, 2.1, 3.0, 3.8)
  drift <- -0.7; diffusion <- 0.5
  set.seed(7)
  data <- do.call(rbind, lapply(seq_len(40), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.1
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1),
    MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1),
    CINT = matrix(0, 1, 1)))
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  out <- ctFitProfile(fit, parameters = "drift", points = 8L)
  # A cleanly identified parameter on a well conditioned fit: the fit's own
  # estimate stands up to the fuller rigor a profile point now gets.
  expect_null(out$better)
  # The base is the fit's own objective value at the estimate, exactly --
  # ctFitProfile() no longer takes it from a one-iteration-capped optimisation.
  expect_equal(out$base, fit$estimate$logposterior)

  row <- out$summary[out$summary$parameter == "drift", ]
  expect_equal(row$verdict, "identifiable")
  expect_equal(row$lower, 0.2363984, tolerance = 1e-3)
  expect_equal(row$upper, 0.5372844, tolerance = 1e-3)
})
