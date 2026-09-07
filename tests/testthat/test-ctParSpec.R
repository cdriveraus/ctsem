# `ctParSpec()` writes the `|`-separated cell specification from named arguments.
#
# The property that matters is that it is not a second specification language:
# whatever it emits must be read by `ctModel()` exactly as the hand-written
# string would be. So the tests compare the two rather than checking the string
# in isolation.

.ctParSpec_model <- function(spec) {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(spec), n.TIpred = 2,
    TIpredNames = c("age", "sex"), Tpoints = 5))
}

.ctParSpec_row <- function(model) {
  columns <- c("param", "transform", "indvarying", "sdscale",
    "age_effect", "sex_effect")
  model$pars[model$pars$param %in% "mm", intersect(columns, names(model$pars))]
}

test_that("it writes the fields in the positions the parser reads", {
  expect_equal(ctParSpec("mm", indvarying = TRUE), "mm||TRUE")
  expect_equal(ctParSpec("d11", transform = "-exp(param)"), "d11|-exp(param)")
  expect_equal(ctParSpec("mm", indvarying = TRUE, sdscale = 0.5), "mm||TRUE|0.5")
  expect_equal(ctParSpec("mm", indvarying = TRUE, tipreds = c("age", "sex")),
    "mm||TRUE||age,sex")
  # Trailing empties carry nothing and only add pipes to read past.
  expect_equal(ctParSpec("plain"), "plain")
})

test_that("what it writes is read exactly as the hand-written string", {
  built <- .ctParSpec_model(ctParSpec("mm", indvarying = TRUE, sdscale = 0.5,
    tipreds = "age"))
  written <- .ctParSpec_model("mm||TRUE|0.5|age")
  expect_identical(.ctParSpec_row(built), .ctParSpec_row(written))
})

test_that("a transform survives the round trip", {
  built <- .ctParSpec_model(ctParSpec("mm", transform = "exp(param)"))
  written <- .ctParSpec_model("mm|exp(param)")
  expect_identical(.ctParSpec_row(built), .ctParSpec_row(written))
  expect_equal(.ctParSpec_row(built)$transform, "exp(param)")
})

test_that("a missing or empty name is refused", {
  expect_error(ctParSpec(), "needs a parameter name")
  expect_error(ctParSpec(""), "needs a parameter name")
  expect_error(ctParSpec(NA), "needs a parameter name")
})

test_that("tipreds takes names, and says so when given something else", {
  # A value there would mean a fixed effect, which the cell syntax expresses
  # and this argument deliberately does not: it would read as an effect to
  # estimate under a different name.
  expect_error(ctParSpec("mm", tipreds = 4.3), "predictor names")
})

# The named fields can also be written straight into a cell. The property that
# matters is the same one: whatever a user writes, the model must come out the
# same. So every test below compares the two forms rather than inspecting one.

.ctParSpec_blankmodel <- function() {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    n.TIpred = 2, TIpredNames = c("age", "sex"), Tpoints = 5))
}

test_that("a named cell and a pipe cell build the same model", {
  expect_identical(
    .ctParSpec_row(.ctParSpec_model("mm, indvarying=TRUE, sdscale=0.5, tipreds=age")),
    .ctParSpec_row(.ctParSpec_model("mm||TRUE|0.5|age")))
  expect_identical(
    .ctParSpec_row(.ctParSpec_model("mm, transform=exp(param)")),
    .ctParSpec_row(.ctParSpec_model("mm|exp(param)")))
  expect_identical(
    .ctParSpec_row(.ctParSpec_model("mm, indvarying=TRUE, tipreds=c(age, sex)")),
    .ctParSpec_row(.ctParSpec_model("mm||TRUE||age,sex")))
})

test_that("a bare name is still a bare name", {
  expect_equal(.ctCellSpecToPipe("mm"), "mm")
  expect_equal(.ctParSpec_row(.ctParSpec_model("mm"))$param, "mm")
})

test_that("a comparison operator in a cell is not a named specification", {
  # Detection keys on the field names, not on a bare `=`, so a state-dependent
  # cell keeps meaning what it says.
  expect_equal(.ctCellSpecToPipe("state[1]>=2"), "state[1]>=2")
  expect_equal(.ctCellSpecToPipe("PARS[1,1]*(state[1]==0)"),
    "PARS[1,1]*(state[1]==0)")
  expect_false(.ctCellSpecIsNamed("dr11"))
})

test_that("a comma inside brackets belongs to the value", {
  expect_equal(.ctCellSpecToPipe("mm, transform=pnorm(param, 0, 1)"),
    "mm|pnorm(param, 0, 1)")
})

test_that("the same cell works through model$matrices <- ", {
  piped <- .ctParSpec_blankmodel()
  piped$matrices$MANIFESTMEANS[1, 1] <- "mm||TRUE|0.5|age"
  named <- .ctParSpec_blankmodel()
  named$matrices$MANIFESTMEANS[1, 1] <- "mm, indvarying=TRUE, sdscale=0.5, tipreds=age"
  expect_identical(.ctParSpec_row(piped), .ctParSpec_row(named))
  expect_true(.ctParSpec_row(named)$indvarying)
  # Read through the accessor, because the TI effect columns are character
  # rather than logical: an effect is TRUE, FALSE, or a value for generation to
  # use ('age=0.4'), and a logical column cannot hold the third. The columns
  # still read 'TRUE'/'FALSE' here, so the change is in the type and not in
  # what this cell says.
  expect_true(ctsem:::.ctTipredEffectActive(.ctParSpec_row(named)$age_effect))
  expect_false(ctsem:::.ctTipredEffectActive(.ctParSpec_row(named)$sex_effect))
})

test_that("an unusable field is refused rather than becoming a parameter name", {
  expect_error(.ctCellSpecToPipe("mm, indvarying=maybe"), "TRUE or FALSE")
  expect_error(.ctCellSpecToPipe("mm, sdscale=big"), "must be a number")
  expect_error(.ctCellSpecToPipe("mm, transform=exp(param), transform=param"),
    "given twice")
  model <- .ctParSpec_blankmodel()
  expect_error(model$matrices$MANIFESTMEANS[1, 1] <- "mm, tipreds=notapredictor",
    "not a time independent predictor")
})
