# `ctPar()` writes the `|`-separated cell specification from named arguments.
#
# The property that matters is that it is not a second specification language:
# whatever it emits must be read by `ctModel()` exactly as the hand-written
# string would be. So the tests compare the two rather than checking the string
# in isolation.

.ctpar_model <- function(spec) {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(spec), n.TIpred = 2,
    TIpredNames = c("age", "sex"), Tpoints = 5))
}

.ctpar_row <- function(model) {
  columns <- c("param", "transform", "indvarying", "sdscale",
    "age_effect", "sex_effect")
  model$pars[model$pars$param %in% "mm", intersect(columns, names(model$pars))]
}

test_that("it writes the fields in the positions the parser reads", {
  expect_equal(ctPar("mm", indvarying = TRUE), "mm||TRUE")
  expect_equal(ctPar("d11", transform = "-exp(param)"), "d11|-exp(param)")
  expect_equal(ctPar("mm", indvarying = TRUE, sdscale = 0.5), "mm||TRUE|0.5")
  expect_equal(ctPar("mm", indvarying = TRUE, tipreds = c("age", "sex")),
    "mm||TRUE||age,sex")
  # Trailing empties carry nothing and only add pipes to read past.
  expect_equal(ctPar("plain"), "plain")
})

test_that("what it writes is read exactly as the hand-written string", {
  built <- .ctpar_model(ctPar("mm", indvarying = TRUE, sdscale = 0.5,
    tipreds = "age"))
  written <- .ctpar_model("mm||TRUE|0.5|age")
  expect_identical(.ctpar_row(built), .ctpar_row(written))
})

test_that("a transform survives the round trip", {
  built <- .ctpar_model(ctPar("mm", transform = "exp(param)"))
  written <- .ctpar_model("mm|exp(param)")
  expect_identical(.ctpar_row(built), .ctpar_row(written))
  expect_equal(.ctpar_row(built)$transform, "exp(param)")
})

test_that("a missing or empty name is refused", {
  expect_error(ctPar(), "needs a parameter name")
  expect_error(ctPar(""), "needs a parameter name")
  expect_error(ctPar(NA), "needs a parameter name")
})

test_that("tipreds takes names, and says so when given something else", {
  # A value there would mean a fixed effect, which the cell syntax expresses
  # and this argument deliberately does not: it would read as an effect to
  # estimate under a different name.
  expect_error(ctPar("mm", tipreds = 4.3), "predictor names")
})
