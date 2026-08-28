# Generating individual differences.
#
# The engine already draws them -- an individually varying parameter becomes a
# state under the augmented layout and its population standard deviation an
# entry of T0VAR -- so what these tests check is the mapping: that a standard
# deviation stated on the parameter's natural scale comes back at that size,
# that an unstated one still produces something, and that the machinery leaves
# fixed-effects models exactly as they were.

.re_model <- function(sd = NA, tpoints = 20) {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.2), MANIFESTVAR = matrix(0.05),
    T0VAR = matrix(0.2), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix("mm"), Tpoints = tpoints))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "mm"] <- TRUE
  if (!is.na(sd)) m$pars$indvaryingsd[m$pars$param %in% "mm"] <- sd
  m
}

# The observed subject means carry within-subject noise as well, so the
# between-subject component has to be separated out before comparing.
.re_between <- function(d) {
  mu <- tapply(d[, "Y1"], d[, "id"], mean)
  within <- mean(tapply(d[, "Y1"], d[, "id"], stats::var))
  n <- nrow(d) / length(mu)
  sqrt(max(0, stats::var(mu) - within / n))
}

test_that("the specification carries a population spread, defaulting to unstated", {
  m <- .re_model()
  expect_true("indvaryingsd" %in% names(m$pars))
  # NA rather than a number: an unstated spread must leave every existing model
  # fitting and generating exactly as it did.
  expect_true(all(is.na(m$pars$indvaryingsd)))
})

test_that("TI predictor effect sizes get a column beside the effect flags", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    n.TIpred = 2, TIpredNames = c("TI1", "TI2"), Tpoints = 5))
  expect_true(all(c("TI1_effect", "TI2_effect") %in% names(m$pars)))
  expect_true(all(c("TI1_effectsize", "TI2_effectsize") %in% names(m$pars)))
  expect_true(all(is.na(m$pars$TI1_effectsize)))
})

test_that("a varying model generates rather than being refused", {
  skip_on_cran()
  skip_without_julia()
  set.seed(3)
  d <- suppressMessages(ctGenerate(.re_model(), n.subjects = 100,
    Tpoints = 20, backend = "julia"))
  expect_equal(nrow(d), 2000L)
  expect_true(all(is.finite(d[, "Y1"])))
  # Unstated spread uses the prior's centre, which is a real spread rather than
  # zero -- a model asking for individual differences must not generate data
  # without them.
  expect_gt(.re_between(d), 0)
})

test_that("a stated population sd comes back at that size", {
  skip_on_cran()
  skip_without_julia()
  for (target in c(0.5, 2)) {
    set.seed(3)
    d <- suppressMessages(ctGenerate(.re_model(target), n.subjects = 400,
      Tpoints = 20, backend = "julia"))
    # 400 subjects puts the standard error of an estimated sd near 3.5%, so the
    # tolerance is sampling noise rather than slack in the mapping.
    expect_equal(.re_between(d), target, tolerance = 0.12)
  }
})

test_that("a varying parameter stays free through preparation", {
  # The bug this guards: resolving free parameters to values made the varying
  # ones fixed, a fixed parameter is not augmented, and the population structure
  # then did not exist at all -- npar came out -Inf and every stated spread was
  # silently ignored.
  m <- .re_model(0.5)
  resolved <- ctsem:::.ctGenerateResolveFree(m, quiet = TRUE)
  expect_true(is.na(resolved$pars$value[resolved$pars$param %in% "mm"]))
  # The value it would have assigned is recorded instead, for the mean.
  expect_true("mm" %in% names(attr(resolved, "ctGenerateMeans")))
})

test_that("a negative population sd is refused", {
  skip_on_cran()
  skip_without_julia()
  expect_error(
    suppressMessages(ctGenerate(.re_model(-1), n.subjects = 5, Tpoints = 5,
      backend = "julia")),
    "cannot be")
})

test_that("transform inversion round trips, and reports linearity", {
  linear <- "10 * param"
  expect_true(ctsem:::.ctGenerateTransformIsLinear(linear))
  expect_equal(ctsem:::.ctGenerateTransformSlope(linear, 0), 10,
    tolerance = 1e-6)
  expect_equal(ctsem:::.ctGenerateTransformInvert(linear, 5), 0.5,
    tolerance = 1e-8)

  positive <- "1e-10 + 1 * log1p_exp(2 * param - 1)"
  expect_false(ctsem:::.ctGenerateTransformIsLinear(positive))
  root <- ctsem:::.ctGenerateTransformInvert(positive, 0.75)
  expect_equal(ctsem:::.ctGenerateTransformAt(positive, root), 0.75,
    tolerance = 1e-8)

  # Unreachable targets are NA rather than a silently clamped endpoint: a
  # negative value for a positive-constrained transform is a specification
  # error, and returning the nearest attainable number would hide it.
  expect_true(is.na(ctsem:::.ctGenerateTransformInvert(positive, -1)))
})

test_that("a nonlinear transform is flagged as first order rather than silent", {
  skip_on_cran()
  skip_without_julia()
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix("dr"), DIFFUSION = matrix(0.2), MANIFESTVAR = matrix(0.05),
    T0VAR = matrix(0.2), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 20))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "dr"] <- TRUE
  m$pars$indvaryingsd[m$pars$param %in% "dr"] <- 0.1
  expect_message(
    suppressWarnings(ctGenerate(m, n.subjects = 20, Tpoints = 20,
      backend = "julia")),
    "first order")
})

test_that("fixed effects models generate exactly as before", {
  skip_on_cran()
  skip_without_julia()
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 8))
  set.seed(1)
  d <- suppressMessages(ctGenerate(m, n.subjects = 60, Tpoints = 8,
    backend = "julia"))
  expect_equal(nrow(d), 480L)
  expect_true(all(is.finite(d[, "Y1"])))
  # No varying parameters, so nothing in this branch should touch it.
  expect_equal(sd(d[, "Y1"]), 0.7755, tolerance = 1e-3)
})
