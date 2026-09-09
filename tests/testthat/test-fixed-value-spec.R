# A number in a matrix cell's first element fixes that cell, and nothing after
# it can apply.
#
# A fixed value cannot individually differ, cannot carry a covariate and cannot
# have a transform: there is no parameter left for any of those to act on. What
# the parser did instead depended on which number was written. A non-negative
# integer slipped through as a *free parameter named after its digits* --
# `'2|param|TRUE|3'` gave param "2", value NA, indvarying TRUE -- while a
# negative or non-integer value hit the check for non-word characters (a minus
# sign and a decimal point are both non-word) and errored. One specification,
# two wrong answers, neither of them a fixed value.

.fv_model <- function(drift, ...) suppressMessages(ctModel(Tpoints = 5,
  LAMBDA = matrix(1), DRIFT = matrix(drift), DIFFUSION = 0.5,
  MANIFESTVAR = 0.5, T0VAR = 1, CINT = 0, MANIFESTMEANS = 0,
  T0MEANS = 0, ...))

.fv_drift <- function(m) m$pars[m$pars$matrix == "DRIFT", ]

test_that("a numeric first element fixes the cell whatever its sign or form", {
  for (spec in list("2|param|TRUE|3", "-0.5|param|TRUE|3",
      "-2|-log1p_exp(-param)|TRUE|.5", "1e-3|param|TRUE|1")) {
    m <- suppressWarnings(.fv_model(spec))
    row <- .fv_drift(m)
    expected <- as.numeric(strsplit(spec, "|", fixed = TRUE)[[1]][1])
    expect_equal(row$value, expected)
    expect_true(is.na(row$param))
    expect_true(is.na(row$transform))
    expect_false(row$indvarying %in% TRUE)
    expect_true(is.na(row$sdscale))
    # And no population spread is offered for it.
    expect_null(m$matrices$RAWPOPVAR)
  }
})

test_that("what was ignored is named in a warning", {
  expect_warning(.fv_model("2|param|TRUE|3"),
    "fixed to 2", fixed = TRUE)
  expect_warning(.fv_model("2|param|TRUE|3"),
    "transform param, indvarying TRUE, sdscale 3", fixed = TRUE)
  # A predictor effect is dropped too, and named.
  expect_warning(m <- suppressMessages(ctModel(Tpoints = 5,
    LAMBDA = matrix(1), DRIFT = matrix("2||||TI1=4.3"), DIFFUSION = 0.5,
    MANIFESTVAR = 0.5, T0VAR = 1, CINT = 0, MANIFESTMEANS = 0, T0MEANS = 0,
    n.TIpred = 1, TIpredNames = "TI1")), "tipred effects TI1=4.3", fixed = TRUE)
  expect_equal(.fv_drift(m)$TI1_effect, "FALSE")
})

test_that("a bare separator with nothing after it warns about nothing", {
  expect_no_warning(m <- .fv_model("-0.5|"))
  expect_equal(.fv_drift(m)$value, -0.5)
})

test_that("a named parameter with separators is untouched", {
  m <- suppressWarnings(.fv_model("dr|param|TRUE|3"))
  row <- .fv_drift(m)
  expect_equal(row$param, "dr")
  expect_true(is.na(row$value))
  expect_equal(row$transform, "param")
  expect_true(row$indvarying)
  expect_equal(as.numeric(row$sdscale), 3)
  expect_equal(rownames(m$matrices$RAWPOPVAR), "dr")
})

test_that("a plain fixed numeric is unchanged", {
  row <- .fv_drift(.fv_model(-0.5))
  expect_equal(row$value, -0.5)
  expect_true(is.na(row$param))
  expect_true(is.na(row$transform))
})

# RAWPOPVAR describes the random effects the model has, so a cell fixed after
# the fact loses its row. The surface and .ctVaryingRows() disagreed before:
# fixing a varying parameter's value removed its random effect everywhere that
# acts on one, while RAWPOPVAR went on offering a spread for it.
test_that("fixing a varying parameter's value withdraws its RAWPOPVAR row", {
  m <- suppressWarnings(.fv_model("dr|param|TRUE|3"))
  expect_equal(rownames(m$matrices$RAWPOPVAR), "dr")
  m$pars$value[m$pars$param %in% "dr"] <- -0.5
  expect_length(ctsem:::.ctModelRawPopVarNames(m$pars), 0L)
  expect_length(ctsem:::.ctVaryingParams(m), 0L)
  expect_null(ctsem:::.ctModelRawPopVarSync(m)[["RAWPOPVAR"]])
})

# The same contradiction reached by assignment rather than by spec string.
#
# Setting `indvarying` directly on `pars` is how nearly every multilevel model
# here is written, so a value assigned to a row that is also flagged varying is
# reachable. `.ctVaryingRows()` treats such a cell as fixed, so neither
# preparation route augments it and the request used to evaporate in silence.
test_that("a fixed named parameter asked to vary says so at fit time", {
  skip_without_julia()
  m <- suppressWarnings(.fv_model("dr|param|TRUE|3"))
  m$pars$value[m$pars$param %in% "dr"] <- -0.5
  dat <- data.frame(id = rep(1:6, each = 5), time = rep(0:4, 6),
    Y1 = stats::rnorm(30))
  expect_warning(suppressMessages(ctFit(dat, m, backend = "julia",
    fit = FALSE)), "cannot vary between individuals", fixed = TRUE)
  expect_warning(suppressMessages(ctFit(dat, m, backend = "julia",
    fit = FALSE)), "dr", fixed = TRUE)
})

test_that("the blanket indvarying idiom stays quiet", {
  skip_without_julia()
  # `model$pars$indvarying <- TRUE` means "vary whatever can vary", and every
  # fixed LAMBDA and T0VAR cell picks up the flag too. Naming those back at
  # the caller would be noise, so only a *named* fixed parameter is reported:
  # a number written in a matrix has no parameter name, whereas a value
  # assigned to a named parameter got there deliberately.
  m <- suppressWarnings(suppressMessages(ctModel(Tpoints = 5,
    LAMBDA = matrix(1), DRIFT = matrix("dr2"), DIFFUSION = 0.5,
    MANIFESTVAR = 0.5, T0VAR = 1, CINT = 0, MANIFESTMEANS = "mm",
    T0MEANS = 0)))
  m$pars$indvarying <- TRUE
  dat <- data.frame(id = rep(1:6, each = 5), time = rep(0:4, 6),
    Y1 = stats::rnorm(30))
  expect_no_warning(suppressMessages(ctFit(dat, m, backend = "julia",
    fit = FALSE)))
})
