# `ctPar()` writes the `|`-separated cell specification from named arguments.
#
# The property that matters is that it is not a second specification language:
# whatever it emits must be read by `ctModel()` exactly as the hand-written
# string would be. So the tests compare the two rather than checking the string
# in isolation.

.ctPar_model <- function(spec) {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(spec), n.TIpred = 2,
    TIpredNames = c("age", "sex"), Tpoints = 5))
}

.ctPar_row <- function(model) {
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
  built <- .ctPar_model(ctPar("mm", indvarying = TRUE, sdscale = 0.5,
    tipreds = "age"))
  written <- .ctPar_model("mm||TRUE|0.5|age")
  expect_identical(.ctPar_row(built), .ctPar_row(written))
})

test_that("a transform survives the round trip", {
  built <- .ctPar_model(ctPar("mm", transform = "exp(param)"))
  written <- .ctPar_model("mm|exp(param)")
  expect_identical(.ctPar_row(built), .ctPar_row(written))
  expect_equal(.ctPar_row(built)$transform, "exp(param)")
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

# The named fields can also be written straight into a cell. The property that
# matters is the same one: whatever a user writes, the model must come out the
# same. So every test below compares the two forms rather than inspecting one.

.ctPar_blankmodel <- function() {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    n.TIpred = 2, TIpredNames = c("age", "sex"), Tpoints = 5))
}

test_that("a named cell and a pipe cell build the same model", {
  expect_identical(
    .ctPar_row(.ctPar_model("mm, indvarying=TRUE, sdscale=0.5, tipreds=age")),
    .ctPar_row(.ctPar_model("mm||TRUE|0.5|age")))
  expect_identical(
    .ctPar_row(.ctPar_model("mm, transform=exp(param)")),
    .ctPar_row(.ctPar_model("mm|exp(param)")))
  expect_identical(
    .ctPar_row(.ctPar_model("mm, indvarying=TRUE, tipreds=c(age, sex)")),
    .ctPar_row(.ctPar_model("mm||TRUE||age,sex")))
})

test_that("a bare name is still a bare name", {
  expect_equal(.ctCellSpecToPipe("mm"), "mm")
  expect_equal(.ctPar_row(.ctPar_model("mm"))$param, "mm")
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
  piped <- .ctPar_blankmodel()
  piped$matrices$MANIFESTMEANS[1, 1] <- "mm||TRUE|0.5|age"
  named <- .ctPar_blankmodel()
  named$matrices$MANIFESTMEANS[1, 1] <- "mm, indvarying=TRUE, sdscale=0.5, tipreds=age"
  expect_identical(.ctPar_row(piped), .ctPar_row(named))
  expect_true(.ctPar_row(named)$indvarying)
  # Read through the accessor, because the TI effect columns are character
  # rather than logical: an effect is TRUE, FALSE, or a value for generation to
  # use ('age=0.4'), and a logical column cannot hold the third. The columns
  # still read 'TRUE'/'FALSE' here, so the change is in the type and not in
  # what this cell says.
  expect_true(ctsem:::.ctTipredEffectActive(.ctPar_row(named)$age_effect))
  expect_false(ctsem:::.ctTipredEffectActive(.ctPar_row(named)$sex_effect))
})

test_that("an unusable field is refused rather than becoming a parameter name", {
  # Not TRUE or FALSE, so a level name -- which only the model can check.
  expect_equal(.ctCellSpecToPipe("mm, indvarying=maybe"), "mm||maybe")
  expect_error(.ctPar_model("mm, indvarying=maybe"), "not an id level")
  expect_error(.ctCellSpecToPipe("mm, sdscale=big"), "must be a number")
  expect_error(.ctCellSpecToPipe("mm, transform=exp(param), transform=param"),
    "given twice")
  model <- .ctPar_blankmodel()
  expect_error(model$matrices$MANIFESTMEANS[1, 1] <- "mm, tipreds=notapredictor",
    "not a time independent predictor")
})

test_that("tipreds names the predictors acting on a cell, written either way", {
  expect_equal(.ctCellSpecToPipe("mm, tipreds=age"), "mm||||age")
  # A bare name, a quoted name and a c() of either all mean the same list.
  expect_equal(.ctCellSpecToPipe("mm, tipreds=c(age,sex)"), "mm||||age,sex")
  expect_equal(.ctCellSpecToPipe("mm, tipreds=c('age','sex')"), "mm||||age,sex")
  expect_equal(.ctCellSpecToPipe('mm, tipreds="age"'), "mm||||age")
  # And they reach the model as the ordered form does.
  expect_identical(
    .ctPar_row(.ctPar_model("mm, tipreds=c(age,sex)")),
    .ctPar_row(.ctPar_model("mm||||age,sex")))
  one <- .ctPar_row(.ctPar_model("mm, tipreds=age"))
  expect_true(ctsem:::.ctTipredEffectActive(one$age_effect))
  # Naming one predictor switches the others off: the list is the whole set of
  # effects on the cell, not an addition to the default.
  expect_false(ctsem:::.ctTipredEffectActive(one$sex_effect))
})

test_that("a misspelled field is refused rather than becoming a parameter name", {
  # This was silent. 'mm, topreds=age' is not the ordered form and named no
  # known field, so it fell through as a parameter literally called
  # "mm, topreds=age", took the model's default tipred effects, and failed much
  # later inside the symbolic differentiator with `Could not retrieve body of
  # '=()'` -- which names neither the cell nor the typo.
  expect_error(.ctCellSpecToPipe("mm, topreds=age"), "'tipreds'")
  expect_error(.ctCellSpecToPipe("mm, tipred=age"), "'tipreds'")
  expect_error(.ctCellSpecToPipe("mm, TIpreds=age"), "'tipreds'")
  expect_error(.ctCellSpecToPipe("mm, indvaring=TRUE"), "'indvarying'")
  expect_error(.ctCellSpecToPipe("mm, sdscal=0.5"), "'sdscale'")
  expect_error(.ctCellSpecToPipe("mm, tranform=exp(param)"), "'transform'")
  # Nothing near enough to suggest still says what the valid fields are.
  expect_error(.ctCellSpecToPipe("mm, notafield=3"), "is not a field")
  # Through ctModel(), not only through the helper.
  expect_error(.ctPar_model("mm, topreds=age"), "'tipreds'")
})

test_that("only a field boundary makes a named specification", {
  # A named argument inside a call belongs to the value, so the bracket-aware
  # split is what decides -- not the presence of an `=` anywhere.
  expect_false(.ctCellSpecIsNamed("mm|pnorm(param, mean=0)"))
  expect_equal(.ctCellSpecToPipe("mm, transform=pnorm(param, mean=0)"),
    "mm|pnorm(param, mean=0)")
  expect_equal(.ctCellSpecToPipe("state[1]<=2"), "state[1]<=2")
  expect_equal(.ctCellSpecToPipe("a!=b"), "a!=b")
  # The ordered form's tipreds field may itself carry `name=value`, so a `|`
  # anywhere settles which form the cell is and the named parser never sees it.
  expect_false(.ctCellSpecIsNamed("mm||||age=4.3,sex=2"))
  expect_equal(.ctCellSpecToPipe("mm||||age=4.3,sex=2"), "mm||||age=4.3,sex=2")
})

test_that("a predictor the model does not have is named in the error", {
  expect_error(.ctPar_model("mm, tipreds=aeg"),
    "not a time independent predictor")
  # And asking for effects in a model with no predictors says that, rather
  # than counting separators the user never wrote.
  nopreds <- function(spec) suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), MANIFESTMEANS = matrix(spec), Tpoints = 5))
  expect_error(nopreds("mm, tipreds=age"),
    "no time independent predictors")
  expect_error(nopreds("mm||||age"), "no time independent predictors")
})

# --- indvarying names levels ------------------------------------------------
#
# With id = c('person', 'couple') a cell says where its parameter varies:
# TRUE at both levels, FALSE at neither, or the levels by name. Every way of
# writing a cell has to land on the same two columns.

.ctPar_nested <- function(...) suppressMessages(ctModel(type = "ct",
  manifestNames = c("Y1", "Y2"), latentNames = c("a", "b"),
  id = c("person", "couple"), LAMBDA = diag(2), MANIFESTVAR = diag(.3, 2),
  DIFFUSION = diag(2), T0VAR = diag(2), CINT = 0, T0MEANS = 0, ...))

.ctPar_levels <- function(m, p) {
  row <- match(p, m$pars$param)
  c(person = m$pars$indvarying[row], couple = m$pars$indvarying_couple[row])
}

test_that("indvarying resolves to every level, none, or the levels named", {
  m <- .ctPar_nested(MANIFESTMEANS = c("m1||TRUE", "m2||couple"),
    DRIFT = c("d1||FALSE", "c12||person,couple", "c21||person", "d2"))
  expect_equal(.ctPar_levels(m, "m1"), c(person = TRUE, couple = TRUE))
  expect_equal(.ctPar_levels(m, "m2"), c(person = FALSE, couple = TRUE))
  expect_equal(.ctPar_levels(m, "d1"), c(person = FALSE, couple = FALSE))
  expect_equal(.ctPar_levels(m, "c12"), c(person = TRUE, couple = TRUE))
  expect_equal(.ctPar_levels(m, "c21"), c(person = TRUE, couple = FALSE))
  # Unspecified, the defaults: intercepts at the subject level, nothing above.
  expect_equal(.ctPar_levels(m, "d2"), c(person = FALSE, couple = FALSE))
})

test_that("the named form, ctPar() and a matrix assignment say the same", {
  piped <- .ctPar_nested(MANIFESTMEANS = c("m1||person,couple", "m2||couple|0.5"))
  named <- .ctPar_nested(MANIFESTMEANS = c("m1, indvarying=c(person, couple)",
    "m2, indvarying=couple, sdscale=0.5"))
  helper <- .ctPar_nested(MANIFESTMEANS = c(
    ctPar("m1", indvarying = c("person", "couple")),
    ctPar("m2", indvarying = "couple", sdscale = 0.5)))
  expect_equal(ctPar("m2", indvarying = "couple", sdscale = 0.5), "m2||couple|0.5")
  cols <- c("param", "indvarying", "indvarying_couple", "sdscale", "sdscale_couple")
  rows <- function(m) m$pars[m$pars$matrix %in% "MANIFESTMEANS", cols]
  expect_equal(rows(named), rows(piped))
  expect_equal(rows(helper), rows(piped))
  # sdscale applies wherever the parameter varies.
  expect_equal(rows(piped)$sdscale_couple[2], 0.5)

  assigned <- .ctPar_nested(MANIFESTMEANS = c("m1", "m2"))
  assigned$matrices$MANIFESTMEANS[1, 1] <- "m1||TRUE"
  assigned$matrices$MANIFESTMEANS[2, 1] <- "m2, indvarying=couple, sdscale=0.5"
  expect_equal(rows(assigned), rows(piped))
  # A cell fixed afterwards varies at no level.
  assigned$matrices$MANIFESTMEANS[2, 1] <- 0
  expect_false(assigned$pars$indvarying_couple[assigned$pars$matrix %in%
    "MANIFESTMEANS" & assigned$pars$row == 2])
})

test_that("a level the model does not have is named, with the ones it has", {
  expect_error(.ctPar_nested(MANIFESTMEANS = c("m1||house", "m2")),
    "names house, which is not an id level of this model. Its levels are: person, couple")
  expect_error(ctPar("m1", indvarying = "person,couple"), "TRUE, FALSE or the names")
})

test_that("with a single id, TRUE and the id's name are the same", {
  one <- function(spec) suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "a", LAMBDA = matrix(1),
    DRIFT = matrix(spec)))
  expect_equal(one("dd||TRUE")$pars, one("dd||id")$pars)
  expect_true(one("dd||id")$pars$indvarying[one("dd||id")$pars$param %in% "dd"])
})
