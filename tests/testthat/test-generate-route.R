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
