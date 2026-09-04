# Names for the raw (unconstrained) parameter vector, julia backend.
#
# These are what a user is shown when something has to point at one raw
# coordinate: the identifiability warning, `fit$identifiability`,
# `ctIdentify()`, `ctReport()`'s profile table, the Laplace correction table,
# the column names of `estimate$rawposterior`. A positional `raw[17]` in that
# list is not a name -- it tells the reader nothing about which part of their
# model is involved -- and three of the four parameters one real identifiability
# warning named were exactly that.
#
# The invariant, which is what most of this file asserts, is cheap and covers
# classes nobody has thought of yet: for a model of any shape the backend can
# build, no reported name is positional. Every model here is prepared with
# `fit = FALSE`, so none of it costs a fit or needs the engine.

.rawnames_data <- function(nsub = 6L, tp = 4L, seed = 3L) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(nsub), function(i)
    data.frame(id = i, time = seq_len(tp) - 1, Y1 = stats::rnorm(tp))))
  d$TI1 <- rep(stats::rnorm(nsub), each = tp)
  d$TI2 <- rep(stats::rnorm(nsub), each = tp)
  d
}

.rawnames_model <- function(tipreds = c("TI1", "TI2")) {
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    TIpredNames = tipreds, T0MEANS = matrix(0))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$matrix == "MANIFESTMEANS"] <- TRUE
  model
}

.rawnames_of <- function(spec) {
  ctsem:::.ctBackendRawParameterNames(list(model_spec = spec),
    ctsem:::.ctBackendNpar(spec))
}

.rawnames_prepare <- function(data, model, ...) {
  suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    fit = FALSE, cores = 1, ...)))
}

test_that("no raw parameter is reported by position, on either random-effect route", {
  data <- .rawnames_data()
  model <- .rawnames_model()
  for (route in c("augmented", "laplace")) {
    spec <- .rawnames_prepare(data, model, intoverpop = route)
    names <- .rawnames_of(spec)
    expect_length(names, ctsem:::.ctBackendNpar(spec))
    # The whole point: not "most of them".
    expect_false(any(grepl("^raw\\[", names)), info = route)
    expect_false(any(is.na(names) | !nzchar(names)), info = route)
    expect_equal(anyDuplicated(names), 0L, info = route)
  }
})

test_that("TI-predictor coefficients are named after the parameter and the predictor", {
  spec <- .rawnames_prepare(.rawnames_data(), .rawnames_model(),
    intoverpop = "augmented")
  names <- .rawnames_of(spec)
  effects <- spec$ti_effects
  expect_gt(nrow(effects), 0L)
  # `rawtipredeffect_<parameter>_<predictor>`, the spelling the stan path's
  # ctFitgetparnamesfromraw() already uses for the same block.
  expect_equal(names[effects$coefficient],
    paste0("rawtipredeffect_", names[effects$parameter], "_",
      spec$model$TIpredNames[effects$predictor]))
  expect_true("rawtipredeffect_drift_eta1_TI1" %in% names)
  expect_true("rawtipredeffect_drift_eta1_TI2" %in% names)
})

test_that("the augmented route's population parameters get their model names, not cell names", {
  data <- .rawnames_data()
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 2, n.manifest = 1, manifestNames = "Y1",
    latentNames = c("eta1", "eta2"), LAMBDA = matrix(c(1, 0), 1, 2),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTMEANS = matrix(0, 1, 1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$matrix == "DRIFT" &
    model$pars$row == model$pars$col] <- TRUE
  spec <- .rawnames_prepare(data, model, intoverpop = "augmented")
  names <- .rawnames_of(spec)
  # The augmented route stores its population scales and correlations as
  # rewritten T0VAR cells, whose own `param` is an internal label
  # (julia_popcov_3_3). That label must not reach a report.
  expect_false(any(grepl("julia_popcov", names)))
  effects <- spec$random_effects
  expect_gt(nrow(effects), 0L)
  sds <- effects[effects$type %in% "sd", , drop = FALSE]
  expect_equal(names[sds$parameter], paste0("popsd_", sds$param))
  correlations <- effects[effects$type %in% "correlation", , drop = FALSE]
  expect_equal(names[correlations$parameter],
    paste0("rawcor_", correlations$param))
  # And the same conceptual parameter is spelled the same way on the Laplace
  # route, which reaches it through its own index vectors instead.
  laplace <- .rawnames_of(.rawnames_prepare(data, model, intoverpop = "laplace"))
  expect_true(all(paste0("popsd_", sds$param) %in% laplace))
  expect_true(all(paste0("rawcor_", correlations$param) %in% laplace))
})

test_that("a sampled missing TI-predictor value names its predictor and subject", {
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix("t0m", 1, 1), n.TIpred = 1, TIpredNames = "grp",
    tipredDefault = FALSE)))
  model$pars$grp_effect[model$pars$param == "t0m"] <- TRUE
  data <- data.frame(id = rep(1:4, each = 3), time = rep(0:2, 4), Y1 = 0,
    grp = rep(c(-1, 2, NA, 0.5), each = 3))
  spec <- .rawnames_prepare(data, model, optimize = FALSE,
    intoverpop = "augmented")
  expect_equal(nrow(spec$ti_missing), 1L)
  names <- .rawnames_of(spec)
  expect_false(any(grepl("^raw\\[", names)))
  expect_match(names[spec$ti_missing$parameter], "^tipredvalue_grp_")
})

test_that("a name that cannot be derived is reported rather than shipped as a name", {
  spec <- .rawnames_prepare(.rawnames_data(), .rawnames_model(),
    intoverpop = "augmented")
  npar <- ctsem:::.ctBackendNpar(spec)
  # Two indices past everything the layout accounts for stand in for a block
  # nobody taught the namer about, which is the only way the fallback can fire.
  expect_warning(
    names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec),
      npar + 2L),
    "could not be named")
  expect_equal(names[npar + 1L], paste0("raw[", npar + 1L, "]"))
  # A fit restored without its parameter table has no naming gap -- it has no
  # spec -- so it must stay quiet.
  expect_silent(bare <- ctsem:::.ctBackendRawParameterNames(list(), 3L))
  expect_equal(bare, c("raw[1]", "raw[2]", "raw[3]"))
})
