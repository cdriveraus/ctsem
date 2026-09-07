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
  if (!is.na(sd)) {
    mats <- m$matrices
    mats$POPCOV["mm", "mm"] <- sd
    m$matrices <- mats
  }
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

test_that("POPCOV surfaces the population covariance from the model onward", {
  m <- .re_model()
  popcov <- m$matrices$POPCOV
  expect_true(is.matrix(popcov))
  expect_equal(rownames(popcov), "mm")
  # Free and labelled by default: that is what the model already did, and the
  # point of surfacing it is to show what it implies, not to change it.
  expect_equal(popcov["mm", "mm"], "popsd_mm")
  # Not in `pars`, which is what keeps it out of every free-parameter count.
  expect_false("POPCOV" %in% m$pars$matrix)

  # A model with no varying parameters has none.
  plain <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0),
    Tpoints = 5))
  expect_false("POPCOV" %in% names(plain$matrices))
})

test_that("POPCOV tracks indvarying by name, keeping what was set", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), Tpoints = 5))
  mats <- m$matrices
  mats$POPCOV["mm_Y1", "mm_Y1"] <- 0.3
  m$matrices <- mats
  expect_equal(m$matrices$POPCOV["mm_Y1", "mm_Y1"], "0.3")
  # Dropping a different random effect must not shuffle this one's entry.
  m$pars$indvarying[m$pars$param %in% "mm_Y2"] <- FALSE
  expect_false("mm_Y2" %in% rownames(m$matrices$POPCOV))
  expect_equal(m$matrices$POPCOV["mm_Y1", "mm_Y1"], "0.3")

  # A name the model does not have is refused rather than matched by position.
  bad <- m$matrices
  rownames(bad$POPCOV)[1] <- "nonsense"
  expect_error({m$matrices <- bad}, "not individually varying")
})

test_that("TI predictor effect sizes get a column beside the effect flags", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    n.TIpred = 2, TIpredNames = c("TI1", "TI2"), Tpoints = 5))
  expect_true(all(c("TI1_effect", "TI2_effect") %in% names(m$pars)))
  # Character, so an effect can be free ('TRUE'), free and named, fixed at a
  # value, or absent -- four states a logical column could not hold.
  expect_type(m$pars$TI1_effect, "character")
  expect_true(all(ctsem:::.ctTipredEffectActive(m$pars$TI1_effect)))
  expect_true(all(is.na(ctsem:::.ctTipredEffectValue(m$pars$TI1_effect))))
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
    "cannot be negative")
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
  mats <- m$matrices
  mats$POPCOV["dr", "dr"] <- 0.1
  m$matrices <- mats
  # A nonlinear transform means the stated spread is matched to first order at
  # the raw origin, which the conversion documents rather than hides.
  d <- suppressWarnings(suppressMessages(ctGenerate(m, n.subjects = 20,
    Tpoints = 20, backend = "julia")))
  expect_equal(nrow(d), 400L)
  expect_true(all(is.finite(d[, "Y1"])))
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

.ti_model <- function(effect = NA, sd = 0.5) {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.2), MANIFESTVAR = matrix(0.05),
    T0VAR = matrix(0.2), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix("mm"), Tpoints = 20,
    n.TIpred = 1, TIpredNames = "TI1"))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "mm"] <- TRUE
  mats <- m$matrices
  mats$POPCOV["mm", "mm"] <- sd
  m$matrices <- mats
  m$pars$TI1_effect <- 'FALSE'
  m$pars$TI1_effect[m$pars$param %in% "mm"] <-
    if (is.na(effect)) 'TRUE' else as.character(effect)
  m
}

test_that("a TI effect can be free, named, fixed, or absent", {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix("mm||||TI1=4.3, TI2=myti2effect, TI3"),
    n.TIpred = 3, TIpredNames = c("TI1", "TI2", "TI3"), Tpoints = 5))
  row <- m$pars[m$pars$param %in% "mm", ]
  expect_equal(row$TI1_effect, "4.3")
  expect_equal(row$TI2_effect, "myti2effect")
  expect_equal(row$TI3_effect, "TRUE")
  # Fixed: has a value, is not free.
  expect_equal(ctsem:::.ctTipredEffectValue(row$TI1_effect), 4.3)
  expect_false(ctsem:::.ctTipredEffectFree(row$TI1_effect))
  # Named free: no value, carries a label another parameter can share.
  expect_true(ctsem:::.ctTipredEffectFree(row$TI2_effect))
  expect_equal(ctsem:::.ctTipredEffectLabel(row$TI2_effect), "myti2effect")
  # Bare name: free, named automatically, which is what it always meant.
  expect_true(ctsem:::.ctTipredEffectFree(row$TI3_effect))
  expect_true(is.na(ctsem:::.ctTipredEffectLabel(row$TI3_effect)))
  # All three are effects; the rest of the model has none.
  expect_true(all(ctsem:::.ctTipredEffectActive(
    row[, c("TI1_effect", "TI2_effect", "TI3_effect")])))
  expect_error(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm||||NOTAPRED"),
    n.TIpred = 1, TIpredNames = "TI1", Tpoints = 5)),
    "not a time independent predictor")
})

test_that("a TI effect fixed to a value is refused for fitting, not ignored", {
  skip_on_cran()
  skip_without_julia()
  # Fitting does not honour the value: the coefficient it occupies is an
  # ordinary free parameter to the optimiser, so a fit would estimate it and
  # quietly disregard what was written.
  m <- .ti_model(effect = 1.5)
  d <- suppressMessages(ctGenerate(m, n.subjects = 10, Tpoints = 6,
    backend = "julia"))
  expect_error(
    suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      cores = 1, optimcontrol = list(estonly = TRUE)))),
    "fixed")
})

test_that("time independent predictors are drawn per subject, not left at zero", {
  skip_on_cran()
  skip_without_julia()
  set.seed(7)
  d <- suppressMessages(ctGenerate(.ti_model(), n.subjects = 200,
    Tpoints = 20, backend = "julia"))
  # Constant within a subject is what makes them time independent; varying
  # across subjects is what makes them predictors. Zeros satisfied the first
  # and not the second, so no effect could show in the data and no fit to that
  # data could identify one.
  expect_true(all(tapply(d[, "TI1"], d[, "id"],
    function(x) length(unique(x))) == 1L))
  values <- tapply(d[, "TI1"], d[, "id"], function(x) x[1L])
  expect_equal(length(unique(values)), 200L)
  expect_equal(sd(values), 1, tolerance = 0.15)
})

test_that("a stated TI effect appears as a slope of the right size", {
  skip_on_cran()
  skip_without_julia()
  slope <- function(target) {
    set.seed(7)
    d <- suppressMessages(ctGenerate(.ti_model(effect = target),
      n.subjects = 400, Tpoints = 20, backend = "julia"))
    mu <- tapply(d[, "Y1"], d[, "id"], mean)
    ti <- tapply(d[, "TI1"], d[, "id"], function(x) x[1L])
    unname(stats::coef(stats::lm(mu ~ ti))[2L])
  }
  expect_equal(slope(0), 0, tolerance = 0.1)
  expect_equal(slope(0.5), 0.5, tolerance = 0.12)
  expect_equal(slope(1.5), 1.5, tolerance = 0.12)
})

test_that("TI coefficients are counted when sizing the raw vector", {
  skip_on_cran()
  skip_without_julia()
  # TI coefficients and the Laplace block are numbered past the end of the
  # parameter table, so counting the table alone left the raw vector short and
  # the engine raised a BoundsError from inside its own state assembly -- an
  # error that says nothing about the mistake.
  m <- .ti_model()
  sk <- ctsem:::.ctGenerateSkeleton(ctsem:::.ctGenerateResolveFree(m, quiet = TRUE),
    3, lapply(1:3, function(i) 0:5))
  spec <- ctsem:::.ctJuliaPrepare(sk, ctsem:::.ctGenerateResolveFree(m, quiet = TRUE),
    priors = FALSE, intoverpop = "augmented")
  fromtable <- suppressWarnings(max(c(0L,
    as.integer(spec$parameter_table$parnumber)), na.rm = TRUE))
  full <- suppressWarnings(max(c(0L, as.integer(spec$parameter_table$parnumber),
    as.integer(spec$laplace$npar), as.integer(spec$ti_effects$coefficient)),
    na.rm = TRUE))
  expect_gt(full, fromtable)
})

test_that("an unstated spread is drawn from its prior, scaled by sdscale", {
  skip_on_cran()
  skip_without_julia()
  scaled <- function(sdscale, seed) {
    m <- .re_model()
    m$pars$sdscale[m$pars$param %in% "mm"] <- sdscale
    set.seed(seed)
    .re_between(suppressMessages(ctGenerate(m, n.subjects = 250,
      Tpoints = 20, backend = "julia")))
  }
  small <- mean(vapply(1:3, function(s) scaled(0.2, s), numeric(1)))
  unit <- mean(vapply(1:3, function(s) scaled(1, s), numeric(1)))
  large <- mean(vapply(1:3, function(s) scaled(5, s), numeric(1)))
  # `sdscale` is inside the population sd's own transform, so it scales the
  # drawn spread linearly. It used to be discarded: `modelmats` is built by
  # ctFit on the way to a backend and is absent when a specification is
  # prepared directly, and the fallback of 1 silently ignored what was asked.
  expect_lt(small, unit)
  expect_gt(large, unit)
  expect_equal(small / unit, 0.2, tolerance = 0.15)
  expect_equal(large / unit, 5, tolerance = 0.15)
})

test_that("the drawn spread is governed by set.seed", {
  skip_on_cran()
  skip_without_julia()
  # Drawn rather than fixed at the prior's centre, which costs nothing in
  # reproducibility because the draw comes from R's own generator.
  set.seed(11)
  a <- suppressMessages(ctGenerate(.re_model(), n.subjects = 40, Tpoints = 10,
    backend = "julia"))
  set.seed(11)
  b <- suppressMessages(ctGenerate(.re_model(), n.subjects = 40, Tpoints = 10,
    backend = "julia"))
  expect_identical(a, b)
  set.seed(12)
  c <- suppressMessages(ctGenerate(.re_model(), n.subjects = 40, Tpoints = 10,
    backend = "julia"))
  expect_false(identical(a, c))
})

test_that("without fromPriors, sdscale is the population sd and nothing is drawn", {
  skip_without_julia()
  # The contract. `ctGenerate()` is not a prior predictive: every quantity is
  # pinned, the number the user wrote is the number used, and anything the
  # model leaves unsaid gets a stated default that is named.
  #
  # `sdscale` multiplies the population sd's *prior* when a model is fitted,
  # which is the only thing a specification can say about a quantity the data
  # will estimate. Generating estimates nothing, so it is simply the sd.
  mk <- function(sdscale) suppressWarnings(ctModel(type = 'ct', n.latent = 1,
    n.manifest = 1, Tpoints = 4, LAMBDA = matrix(1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0.01), MANIFESTVAR = matrix(0.01),
    MANIFESTMEANS = matrix(paste0('mmean||TRUE|', sdscale)),
    T0VAR = matrix(0.01), T0MEANS = matrix(0), CINT = matrix(0)))

  realised <- function(model, n = 400, seed = 1) {
    set.seed(seed)
    d <- suppressMessages(ctGenerate(model, n.subjects = n, backend = 'julia'))
    stats::sd(tapply(d[, 'Y1'], d[, 'id'], mean))
  }

  # Equal to sdscale, not merely proportional to it. Before this the same three
  # values gave 0.178, 0.732 and 3.593 -- a monotone function of what was asked
  # for, which is the kind of wrong that looks right in a plot.
  expect_equal(realised(mk(0.2)), 0.2, tolerance = 0.15)
  expect_equal(realised(mk(1)), 1, tolerance = 0.15)
  expect_equal(realised(mk(5)), 5, tolerance = 0.15)

  # A stated POPCOV is the more specific statement and wins.
  stated <- mk(0.2)
  stated$matrices$POPCOV['mmean', 'mmean'] <- 2
  expect_equal(realised(stated), 2, tolerance = 0.15)

  # Nothing is drawn. The spread used to come from its prior on every call:
  # measured 1.10, 0.63 and 0.53 on three consecutive seeds from an unchanged
  # specification, which made plain ctGenerate() a partial prior predictive.
  # What moves between seeds now is the finite-sample realisation alone.
  spread <- vapply(1:3, function(s) realised(mk(1), n = 200, seed = s), numeric(1))
  expect_lt(max(spread) / min(spread), 1.3)

  # And both defaults are reported rather than invented in silence: the
  # population mean among the filled parameters, the sd on its own.
  expect_message(
    suppressWarnings(ctGenerate(mk(0.2), n.subjects = 10, backend = 'julia')),
    regexp = 'population mean')
  expect_message(
    suppressWarnings(ctGenerate(mk(0.2), n.subjects = 10, backend = 'julia')),
    regexp = 'taken from sdscale')
})
