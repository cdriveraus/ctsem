suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

# `T0VARredundancies()` (R/ctFit.R) fixes the free T0VAR cells that indvarying
# T0MEANS makes redundant. It has to clear the `<TIpred>_effect` columns along
# with `param`, `transform` and `indvarying`, because `ctStanData()` reads those
# columns for every row of a matrix with no free-parameter filter: a stale
# effect puts T0VAR into `standata$subindices`, and stan then builds a
# per-subject T0VAR whose value cannot vary by subject -- the same answer,
# computed once per subject. The julia path reused the same cells for its
# population parameters and the stale flag was worse there; see
# test-backend-priors-scores.R.
#
# The column holds a character spec rather than a logical (R/ctTipredEffect.R),
# so ask `.ctTipredEffectActive()` instead of testing it directly: a cell
# written `TI1=4.3` is an active effect, and `any()` on that string is NA.
#
# The T0VAR slot of `subindices` is indexed by the matrix code from
# `.ctMatricesList()$base` (R/ctModelWriter.R), which is 8. T0MEANS is 1.
.t0varred_T0VAR_slot <- 8L
.t0varred_T0MEANS_slot <- 1L

.t0varred_model <- function(...) {
  suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), MANIFESTVAR = diag(.2, 2), MANIFESTMEANS = matrix(0, 2, 1),
    n.TIpred = 1, TIpredNames = "TI1", ...))
}

.t0varred_data <- function() {
  set.seed(3)
  do.call(rbind, lapply(1:10, function(i) data.frame(id = i, time = 0:4,
    Y1 = rnorm(5), Y2 = rnorm(5), TI1 = rep(rnorm(1), 5))))
}

.t0varred_standata <- function(model) {
  suppressMessages(ctFit(.t0varred_data(), model, backend = "stan",
    fit = FALSE, priors = FALSE, verbose = 0))$standata
}

# The free T0VAR cells, and whether each still carries an active TI effect.
.t0varred_free <- function(model) {
  which(model$pars$matrix == "T0VAR" & is.na(model$pars$value))
}
.t0varred_active <- function(pars, rows) {
  ctsem:::.ctTipredEffectActive(pars[rows, "TI1_effect"])
}

test_that("T0VARredundancies clears the TI effect flags on the cells it fixes", {
  m <- .t0varred_model()
  free <- .t0varred_free(m)
  expect_gt(length(free), 0)
  expect_true(all(.t0varred_active(m$pars, free)))  # candidates before fixing

  red <- ctsem:::T0VARredundancies(m)
  expect_true(all(is.na(red$pars$param[free])))
  expect_true(all(is.na(red$pars$transform[free])))
  expect_false(any(red$pars$indvarying[free]))
  expect_false(any(.t0varred_active(red$pars, free)))
})

# A cell can carry a fixed effect size or an effect name, not just TRUE
# (R/ctTipredEffect.R). Those are the specs a user writes to generate data or to
# constrain two parameters to one effect, and they have to be cleared too -- a
# T0VAR cell that is no longer a parameter cannot carry either.
test_that("T0VARredundancies clears a fixed-value and a named TI effect spec", {
  m <- .t0varred_model(T0VAR = matrix(c("t0var11||||TI1=4.3", 0,
    "t0var21||||TI1=sharedeffect", "t0var22"), 2, 2, byrow = TRUE))
  free <- .t0varred_free(m)
  expect_setequal(m$pars$TI1_effect[free], c("4.3", "sharedeffect", "TRUE"))

  red <- ctsem:::T0VARredundancies(m)
  expect_false(any(.t0varred_active(red$pars, free)))
})

test_that("fixed T0VAR cells do not ask stan for a per-subject T0VAR", {
  skip_on_cran()
  sd <- .t0varred_standata(.t0varred_model())
  expect_equal(sd$subindices[.t0varred_T0VAR_slot], 0L)
  # not a blanket zero -- T0MEANS does still vary by subject in this model
  expect_equal(sd$subindices[.t0varred_T0MEANS_slot], 1L)
})

test_that("a T0VAR that stayed free still gets a per-subject T0VAR", {
  skip_on_cran()
  m <- .t0varred_model()
  m$pars$indvarying <- FALSE  # nothing left for T0VARredundancies to fix
  sd <- .t0varred_standata(m)
  expect_equal(sd$subindices[.t0varred_T0VAR_slot], 1L)
})
