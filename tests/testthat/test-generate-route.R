# User side generation: which route `ctGenerate()` takes, and what it returns.
#
# Distinct from test-backend-generate.R, which is the posterior predictive from
# a *fit* and has to match that fit's own conventions. This file is the
# specification-to-data direction, where the conventions are `ctGenerate`'s own.

.route_model <- function(...) {
  suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 5,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(1), MANIFESTVAR = matrix(0.1), T0VAR = matrix(0.5),
    MANIFESTMEANS = matrix(0), ...)))
}

# The id and time columns are named as the model names them.
#
# Both generators hardcoded 'id' and 'time'. On the julia route that made every
# model with a custom name ungeneratable rather than merely oddly named: the
# skeleton's columns are looked up by `model$subjectIDname` and
# `model$timeName`, neither was found, and preparation died inside `order()`
# with "argument 1 is not a vector". On the R route the data came back with
# columns the model could not be fitted to without renaming them first.
test_that("generated column names follow the model, backend r", {
  d <- ctGenerate(.route_model(id = "subject", time = "age"), n.subjects = 3,
    Tpoints = 4, backend = "r")
  expect_true(all(c("subject", "age") %in% colnames(d)))
  expect_false(any(c("id", "time") %in% colnames(d)))
})

test_that("generated column names follow the model, backend julia", {
  skip_without_julia()
  d <- suppressMessages(ctGenerate(.route_model(id = "subject", time = "age"),
    n.subjects = 3, Tpoints = 4, backend = "julia"))
  expect_true(all(c("subject", "age") %in% colnames(d)))
  expect_false(any(c("id", "time") %in% colnames(d)))
  expect_equal(nrow(d), 12L)
  expect_true(all(is.finite(d[, "Y1"])))
})

# A default-named model is the case every existing caller has, and it keeps the
# names it always had -- the point of the fix is the lookup, not a rename.
test_that("default names are unchanged on both routes", {
  d <- ctGenerate(.route_model(), n.subjects = 3, Tpoints = 4, backend = "r")
  expect_true(all(c("id", "time") %in% colnames(d)))
  skip_without_julia()
  dj <- suppressMessages(ctGenerate(.route_model(), n.subjects = 3,
    Tpoints = 4, backend = "julia"))
  expect_true(all(c("id", "time") %in% colnames(dj)))
})

# Burnin trimming and the wide conversion both indexed the time column by name,
# so they are the two places a custom name would have failed after the skeleton
# was fixed.
test_that("burnin and wide output honour a custom time name", {
  skip_without_julia()
  m <- .route_model(id = "subject", time = "age")
  d <- suppressMessages(ctGenerate(m, n.subjects = 2, Tpoints = 4,
    burnin = 3, backend = "julia"))
  expect_equal(nrow(d), 8L)
  # Time restarts at zero for each subject once the burnin is dropped.
  expect_equal(unname(d[1, "age"]), 0)
  expect_equal(unname(d[5, "age"]), 0)
  w <- suppressMessages(ctGenerate(m, n.subjects = 2, Tpoints = 4,
    wide = TRUE, backend = "julia"))
  expect_equal(nrow(w), 2L)
})

# Fixed values only: no random effects at any level.
#
# The between-subject spread used to be drawn from the augmented layout's own
# T0 draw -- measured 2.18 against a stated 2 -- through a RAWPOPVAR entry read on
# the parameter's natural scale, which is not the scale a fit reports. Rather
# than carry a convention that disagrees with fitting, user side generation
# draws nothing and says so; ctGenerateFromFit() is where random effects in
# generated data come from.
.route_varying <- function(sd = NA) {
  m <- suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 8,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), DRIFT = matrix(-0.4),
    DIFFUSION = matrix(0.2), MANIFESTVAR = matrix(0.05), T0VAR = matrix(0.2),
    MANIFESTMEANS = matrix("mm"))))
  m$pars$indvarying <- m$pars$param %in% "mm"
  m <- ctsem:::.ctModelRawPopVarSync(m)
  if (!is.na(sd)) {
    mats <- m$matrices
    mats$RAWPOPVAR["mm", "mm"] <- sd
    m$matrices <- mats
  }
  m
}

test_that("a declared random effect is ignored, and named", {
  skip_without_julia()
  expect_message(ctGenerate(.route_varying(sd = 2), n.subjects = 3,
    Tpoints = 5, backend = "julia"),
    "Individual differences are ignored for mm")
})

test_that("no between-subject spread is generated whatever the stated sd", {
  skip_without_julia()
  spread <- function(sd) {
    set.seed(9)
    d <- suppressMessages(ctGenerate(.route_varying(sd = sd), n.subjects = 40,
      Tpoints = 8, backend = "julia"))
    stats::sd(tapply(d[, "Y1"], d[, "id"], mean))
  }
  # A stated sd of 4 is large against the within-subject scale here, so if any
  # of it reached the data the two would differ far beyond sampling noise.
  # Identical, because the same seed generates the same fixed-effects dataset.
  expect_equal(spread(0.5), spread(4))
})

# A multilevel model generates rather than erroring, from its fixed values.
test_that("a grouping level generates from fixed values and says so", {
  skip_without_julia()
  m <- suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 5,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(1), MANIFESTVAR = matrix(0.1), T0VAR = matrix(0.5),
    MANIFESTMEANS = matrix("mm"), id = c("subject", "study"))))
  m$pars$indvarying <- FALSE
  m$pars$indvarying_study <- m$pars$param %in% "mm"
  expect_message(d <- ctGenerate(m, n.subjects = 4, Tpoints = 5,
    backend = "julia"), "Individual differences are ignored for mm")
  expect_true(all(c("subject", "study", "Y1") %in% colnames(d)))
  expect_true(all(is.finite(d[, "Y1"])))
})

# Which route 'auto' takes, and how the two compare on one specification,
# lives in test-julia-intoverstates.R -- that file owns the comparison.

# A time independent predictor effect goes with the random effects, and the
# silence would be worse than the loss: the predictor column is still drawn
# and still varies between subjects, so the data looks like data with a
# predictor in it. Measured on a model stating TI1=4.3, the correlation
# between the subject means and the predictor came out -0.19 over 60
# subjects, which is noise.
test_that("a TI predictor effect is ignored, and named", {
  skip_without_julia()
  m <- suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 6,
    manifestNames = "Y1", latentNames = "eta1", n.TIpred = 1,
    TIpredNames = "TI1", LAMBDA = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), DRIFT = matrix(-0.4), DIFFUSION = matrix(0.2),
    MANIFESTVAR = matrix(0.05), T0VAR = matrix(0.2),
    MANIFESTMEANS = matrix("mm||TRUE|1|TI1=4.3"))))
  expect_message(ctGenerate(m, n.subjects = 4, Tpoints = 6,
    backend = "julia"), "Effects of TI1 are ignored")

  set.seed(2)
  d <- suppressMessages(ctGenerate(m, n.subjects = 60, Tpoints = 6,
    backend = "julia"))
  # The column is there and varies, and carries nothing.
  expect_gt(stats::sd(tapply(d[, "TI1"], d[, "id"], mean)), 0.5)
  expect_lt(abs(stats::cor(tapply(d[, "Y1"], d[, "id"], mean),
    tapply(d[, "TI1"], d[, "id"], mean))), 0.4)
})
