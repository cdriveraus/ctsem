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
  d <- ctGenerate(.route_model(id = "subject", time = "age"), n = 3,
    Tpoints = 4, backend = "r")
  expect_true(all(c("subject", "age") %in% colnames(d)))
  expect_false(any(c("id", "time") %in% colnames(d)))
})

test_that("generated column names follow the model, backend julia", {
  skip_without_julia()
  d <- suppressMessages(ctGenerate(.route_model(id = "subject", time = "age"),
    n = 3, Tpoints = 4, backend = "julia"))
  expect_true(all(c("subject", "age") %in% colnames(d)))
  expect_false(any(c("id", "time") %in% colnames(d)))
  expect_equal(nrow(d), 12L)
  expect_true(all(is.finite(d[, "Y1"])))
})

# A default-named model is the case every existing caller has, and it keeps the
# names it always had -- the point of the fix is the lookup, not a rename.
test_that("default names are unchanged on both routes", {
  d <- ctGenerate(.route_model(), n = 3, Tpoints = 4, backend = "r")
  expect_true(all(c("id", "time") %in% colnames(d)))
  skip_without_julia()
  dj <- suppressMessages(ctGenerate(.route_model(), n = 3,
    Tpoints = 4, backend = "julia"))
  expect_true(all(c("id", "time") %in% colnames(dj)))
})

# Burnin trimming and the wide conversion both indexed the time column by name,
# so they are the two places a custom name would have failed after the skeleton
# was fixed.
test_that("burnin and wide output honour a custom time name", {
  skip_without_julia()
  m <- .route_model(id = "subject", time = "age")
  d <- suppressMessages(ctGenerate(m, n = 2, Tpoints = 4,
    burnin = 3, backend = "julia"))
  expect_equal(nrow(d), 8L)
  # Time restarts at zero for each subject once the burnin is dropped.
  expect_equal(unname(d[1, "age"]), 0)
  expect_equal(unname(d[5, "age"]), 0)
  w <- suppressMessages(ctGenerate(m, n = 2, Tpoints = 4,
    wide = TRUE, backend = "julia"))
  expect_equal(nrow(w), 2L)
})

# Random effects are drawn from the population the model states: a
# RAWPOPVAR number is the raw-scale sd, exactly as a fit reports it.
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

test_that("a stated population sd reaches the data, and is named", {
  skip_without_julia()
  expect_message(ctGenerate(.route_varying(sd = 2), n = 3,
    Tpoints = 5, backend = "julia"), "population sd \\(raw\\): mm 2 \\[id\\]")
  spread <- function(sd) {
    set.seed(9)
    d <- suppressMessages(ctGenerate(.route_varying(sd = sd), n = 40,
      Tpoints = 8, backend = "julia"))
    stats::sd(tapply(d[, "Y1"], d[, "id"], mean))
  }
  # One seed gives the same standard normals, so the subject means scale with
  # the stated sd; the within-subject noise is small beside either.
  expect_equal(spread(4) / spread(0.5), 8, tolerance = 0.05)
})

test_that("a grouping level draws its own effects", {
  skip_without_julia()
  m <- suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 5,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(1), MANIFESTVAR = matrix(0.1), T0VAR = matrix(0.5),
    MANIFESTMEANS = matrix("mm"), id = c("subject", "study"))))
  m$pars$indvarying <- FALSE
  m$pars$indvarying_study <- m$pars$param %in% "mm"
  set.seed(1)
  expect_message(d <- ctGenerate(m, n = c(40, 4), Tpoints = 5,
    backend = "julia"), "mm 1 \\[study\\]")
  expect_true(all(c("subject", "study", "Y1") %in% colnames(d)))
  # Subjects within a study share its effect: the study means spread far more
  # than the subjects within one.
  means <- tapply(d[, "Y1"], d[, "subject"], mean)
  study <- tapply(d[, "study"], d[, "subject"], `[`, 1)
  expect_gt(stats::sd(tapply(means, study, mean)),
    3 * mean(tapply(means, study, stats::sd)))
})

# A TI predictor effect fixed in the model is applied as a fit applies it: a
# shift of the raw parameter per unit of the predictor.
test_that("a fixed TI predictor effect reaches the data", {
  skip_without_julia()
  m <- suppressMessages(suppressWarnings(ctModel(type = "ct", Tpoints = 6,
    manifestNames = "Y1", latentNames = "eta1", n.TIpred = 1,
    TIpredNames = "TI1", LAMBDA = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), DRIFT = matrix(-0.4), DIFFUSION = matrix(0.2),
    MANIFESTVAR = matrix(0.05), T0VAR = matrix(0.2),
    MANIFESTMEANS = matrix("mm||TRUE|1|TI1=0.3"))))
  m$RAWPOPVAR["mm", "mm"] <- 0.01
  set.seed(2)
  d <- suppressMessages(ctGenerate(m, n = 60, Tpoints = 6,
    backend = "julia"))
  ti <- tapply(d[, "TI1"], d[, "id"], mean)
  expect_gt(stats::sd(ti), 0.5)
  # mm = 10 * param: 0.3 raw is a slope of 3.
  expect_equal(unname(stats::coef(stats::lm(tapply(d[, "Y1"], d[, "id"], mean) ~
    ti))[2]), 3, tolerance = 0.1)
})
