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

# T0cov is two covariances in one matrix, and the indices interleave.
#
# The population block is indexed by `intoverpopindvaryingindex`, and for an
# individually varying T0MEANS that index is the *main latent state* -- no
# carrier state is appended for it. So with eta1 varying and eta2 not, row 1 of
# T0cov comes from RAWPOPVAR and row 2 from T0VAR, and a random effect on a
# non-T0MEANS cell adds a third row from RAWPOPVAR after the latents. Getting
# that mapping wrong reads a plausible covariance out of the wrong matrix, so
# it is checked directly, at a fixed raw vector, with no optimiser involved.
test_that("T0cov takes the main latents from T0VAR and the population block from RAWPOPVAR", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  model <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), MANIFESTVAR = diag(.2, 2),
    MANIFESTMEANS = matrix(0, 2, 1),
    # eta1 varies, eta2 does not; and a non-T0MEANS effect that does get its
    # own appended state.
    T0MEANS = matrix(c("t0a||TRUE", 0), 2, 1),
    CINT = matrix(c("b1||TRUE", 0), 2, 1),
    T0VAR = matrix(c("t0v11", 0, "t0v21", "t0v22"), 2, 2, byrow = TRUE),
    DRIFT = matrix(c("dr1", 0, 0, "dr2"), 2, 2)))
  set.seed(4)
  data <- data.frame(id = rep(1:8, each = 4), time = rep(0:3, 8),
    Y1 = stats::rnorm(32), Y2 = stats::rnorm(32))

  prep <- suppressMessages(suppressWarnings(
    ctFit(data, model, fit = FALSE, cores = 1L, verbose = 0L)))
  sdat <- prep$standata
  popidx <- as.integer(sdat$intoverpopindvaryingindex)
  # eta1 is a varying T0MEANS, so the population block starts at state 1
  # rather than after the latents. That is the overlap this test exists for.
  expect_true(1L %in% popidx)
  expect_false(2L %in% popidx)
  expect_equal(length(popidx), 2L)

  sf <- ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, sdat)
  npar <- rstan::get_num_upars(sf)
  set.seed(9)
  cp <- rstan::constrain_pars(sf, stats::rnorm(npar, 0, .4))
  grab <- function(x) if (length(dim(x)) == 3) x[1, , ] else drop(x)
  T0cov <- grab(cp$pop_T0cov)
  popcov <- grab(cp$rawpopcov)

  # The population block is the constructed RAWPOPVAR, and no conversion
  # happens here.
  #
  # Each state has a scale: 10 for a varying T0MEANS, whose state carries
  # natural units, and 1 for an appended carrier, which keeps raw units
  # because the augmentation writes its T0MEANS with the identity transform
  # and the cell reading the state does the scaling. Those factors are applied
  # inside the population standard deviation -- in the diagonal element's own
  # transform -- so rawpopcov is already in state units and T0cov's block
  # equals it outright. They used to be applied here instead, to T0cov's rows
  # and columns after construction, which is why the two backends disagreed
  # about pop_T0VAR and why this block used to differ by 100, 10 and 1.
  ms <- sdat$matsetup
  mv <- sdat$matvalues
  t0meansrows <- which(ms[, 7] == .t0varred_T0MEANS_slot & ms[, 2] == 1L)
  scale <- vapply(popidx, function(state) {
    ri <- t0meansrows[ms[t0meansrows, 1] == state][1L]
    if (is.na(ri)) 1 else mv[ri, 2] * mv[ri, 3]
  }, numeric(1))
  expect_equal(scale[1], 10, tolerance = 1e-8)   # eta1's T0MEANS
  expect_equal(scale[2], 1, tolerance = 1e-8)    # b1's carrier
  block <- popcov[seq_along(popidx), seq_along(popidx), drop = FALSE]
  expect_equal(T0cov[popidx, popidx], block, tolerance = 1e-12)
  # The correlation is what a user reads and is scale free, so it must match
  # regardless of the units above.
  expect_equal(cov2cor(T0cov[popidx, popidx]), cov2cor(block),
    tolerance = 1e-10)
  # eta2 keeps its own T0VAR variance, which RAWPOPVAR knows nothing about
  expect_gt(T0cov[2, 2], 0)
  # and the pair spanning both matrices is the disabled off-diagonal: zero,
  # because RAWPOPVAR does not span it and T0VAR was told not to state it
  expect_equal(T0cov[1, 2], 0, tolerance = 1e-12)
  expect_equal(T0cov[2, 1], 0, tolerance = 1e-12)
  # T0cov stays a covariance matrix across the join
  expect_true(all(is.finite(T0cov)))
  expect_gt(min(eigen(T0cov, symmetric = TRUE, only.values = TRUE)$values), 0)
})

# The rule is level-specific and that is load-bearing.
#
# A T0MEANS that varies over studies but not over subjects must KEEP its
# T0VAR: within a study every subject shares that T0MEANS, so T0VAR is still
# what disperses their initial states, and the study effect moves the whole
# study rather than widening any individual. Disabling on "varies at any
# level" would delete a parameter the data can identify, and nothing else in
# the suite would notice.
#
# Note that T0MEANS is indvarying by default, so the no-effect case has to
# turn it off explicitly -- leaving it alone tests the subject-level case
# twice, which is how this test first came out green for the wrong reason.
test_that("only a subject level T0MEANS random effect disables T0VAR", {
  skip_on_cran()
  mk <- function() suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 2, manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2),
    MANIFESTVAR = diag(.2, 2), MANIFESTMEANS = matrix(0, 2, 1),
    DRIFT = matrix(c("dr1", 0, 0, "dr2"), 2, 2),
    T0MEANS = matrix(c("t0a", 0), 2, 1),
    T0VAR = matrix(c("t0v11", 0, "t0v21", "t0v22"), 2, 2, byrow = TRUE),
    id = c("id", "study")))

  m <- mk()
  # the higher level is expressible at all, which the rest of this depends on
  expect_true("indvarying_study" %in% names(m$pars))
  expect_identical(m$groupIDnames, "study")

  nfree <- function(mm) sum(mm$pars$matrix %in% "T0VAR" & is.na(mm$pars$value))

  none <- mk(); none$pars$indvarying[none$pars$param %in% "t0a"] <- FALSE
  subj <- mk(); subj$pars$indvarying[subj$pars$param %in% "t0a"] <- TRUE
  study <- mk(); study$pars$indvarying[study$pars$param %in% "t0a"] <- FALSE
  study$pars$indvarying_study[study$pars$param %in% "t0a"] <- TRUE

  expect_equal(nfree(suppressMessages(ctsem:::T0VARredundancies(none))),
    nfree(none))
  expect_lt(nfree(suppressMessages(ctsem:::T0VARredundancies(subj))),
    nfree(subj))
  expect_equal(nfree(suppressMessages(ctsem:::T0VARredundancies(study))),
    nfree(study))

  # and precisely: the subject level case drops eta1 row and column, keeping
  # eta2 own variance, which RAWPOPVAR knows nothing about
  out <- suppressMessages(ctsem:::T0VARredundancies(subj))
  kept <- out$pars[out$pars$matrix %in% "T0VAR" & is.na(out$pars$value), ]
  expect_identical(as.character(kept$param), "t0v22")
})
