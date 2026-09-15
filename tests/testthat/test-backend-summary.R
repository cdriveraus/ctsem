# Transformed-parameter summaries for backend='julia'.
#
# The load-bearing claim of this architecture is that the engines produce the
# *same* pop_* arrays Stan does, so that ctSummaryMatrices(), summary() and
# ctDiscretePars() can be one implementation rather than three. The first test
# checks exactly that, at a fixed raw parameter vector, against Stan's own
# constrain step -- which is a much sharper check than comparing two fits, since
# it removes the optimizer from the comparison entirely.
#
# pop_T0VAR is in that comparison, and its being there is the point of the
# population covariance having become a matrix of its own. It used to be
# excluded: the engines folded the population scale into the augmented T0VAR
# cells, so their T0VAR was the matrix whose sdcovsqrt2cov equalled the reported
# T0cov, and stan's was the model's own -- two different quantities under one
# name, comparable only through T0cov. Both sides now report the model's own
# T0VAR and agree to machine precision.

.summary_model <- function() {
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2),
    T0MEANS = c("t0a||TRUE", "t0b||TRUE"), CINT = c("B1||TRUE", "B2||TRUE"),
    DRIFT = matrix(c("auto1", "cross21||TRUE", "cross21||TRUE", "auto2"), 2, 2,
      byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE))
  model$pars$group_effect[model$pars$param == "B1"] <- TRUE
  model
}

.summary_data <- function() {
  set.seed(5)
  do.call(rbind, lapply(1:12, function(i) data.frame(id = i, time = c(0, .5, 1.5, 2.4),
    Y1 = stats::rnorm(4, 0, .5), Y2 = stats::rnorm(4, 0, .5),
    group = rep(stats::rnorm(1), 4))))
}

.summary_pointfit <- function(spec, model, raw, backend) {
  structure(list(model_spec = spec, model = model, backend = backend,
    estimate = list(raw = raw, loglik = NA_real_)),
    class = c("ctJuliaFit", "ctFit"))
}

test_that("Julia pop_* arrays match Stan's constrained parameters", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .summary_model()
  data <- .summary_data()

  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient), na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE))
  stan_pop <- suppressMessages(ctsem:::stan_constrainsamples(sm = ctsem:::stanmodels$ctsm,
    standata = stan_spec$standata, samples = matrix(raw, nrow = 1), cores = 1,
    pcovn = 500, dokalman = FALSE, savesubjectmatrices = FALSE))

  fit <- .summary_pointfit(spec, model, raw, "julia")
  backend_pop <- ctsem:::.ctBackendPopArrays(fit)

  compared <- intersect(grep("^pop_", names(stan_pop), value = TRUE),
    names(backend_pop))
  # A model with an intoverpop augmentation, TI predictors and a state-dependent
  # DRIFT: if this list ever shrinks, the comparison below has stopped covering
  # the interesting matrices and the test has quietly weakened.
  expect_true(all(c("pop_DRIFT", "pop_DIFFUSIONcov", "pop_T0cov", "pop_T0VAR",
    "pop_asymCINT", "pop_asymDIFFUSIONcov", "pop_CINT", "pop_LAMBDA") %in%
    compared))
  for (name in compared) {
    expect_equal(dim(backend_pop[[name]]), dim(stan_pop[[name]]), info = name)
    expect_equal(as.numeric(backend_pop[[name]]), as.numeric(stan_pop[[name]]),
      tolerance = 1e-8, info = name)
  }

  # A state whose T0MEANS is a random effect has no initial covariance of its
  # own: T0VAR's row and column for it are disabled, and its entry in T0cov
  # comes from the population block. Both T0MEANS are individually varying
  # here, so both latent rows are disabled, and so is every carrier state --
  # which leaves T0VAR entirely zero and T0cov entirely population. Asserting
  # it this way rather than by index keeps the test honest if the fixture
  # changes: `random_effects` is where the population rows are named.
  population_rows <- sort(unique(spec$random_effects$row[
    spec$random_effects$type %in% "sd"]))
  t0var <- drop(backend_pop$pop_T0VAR)
  t0cov <- drop(backend_pop$pop_T0cov)
  expect_gt(length(population_rows), 0)
  expect_equal(as.numeric(t0var[population_rows, ]),
    rep(0, length(population_rows) * ncol(t0var)))
  expect_equal(as.numeric(t0var[, population_rows]),
    rep(0, nrow(t0var) * length(population_rows)))
  # And the block T0VAR does not state is a covariance nonetheless, so the
  # zeros above are a reparameterisation and not a lost variance.
  expect_true(all(diag(t0cov)[population_rows] > 0))

  # The state-unit conversion, at the one point it happens. A carrier state
  # holds its effect in the units the consuming cell reads, so T0cov's
  # diagonal for a carrier is that effect's raw population sd times the
  # carrier factor: the cell multiplier*meanscale for an individually varying
  # T0MEANS, whose carrier is the model latent itself, and one for an appended
  # carrier, whose T0MEANS uses the identity transform. Exact arithmetic, so no
  # tolerance, and fit-free -- which is the point, because the failure mode
  # here leaves the fit correct and only what a user reads is wrong.
  sds <- spec$random_effects$type %in% "sd"
  scales <- as.numeric(spec$random_effects$scale[sds])
  carrier_rows <- as.integer(spec$random_effects$row[sds])
  rawsd <- as.numeric(drop(stan_pop$rawpopsd))
  # Without a factor other than one the two assertions below say nothing.
  expect_true(any(scales != 1))
  expect_equal(sqrt(diag(drop(stan_pop$pop_T0cov))[carrier_rows]),
    rawsd * scales)

  # And it is applied once. `popsd` is the spread of the *transformed*
  # parameter, drawn from rawpopmeans and the factorisation of the population
  # covariance, and a varying T0MEANS has a linear cell transform whose slope
  # is exactly that factor -- so its reported spread is the raw sd times the
  # factor and not the factor twice. Loose because it is a standard deviation
  # over pcovn draws; the quantity being ruled out is an order of magnitude,
  # which is what this reported when the factorisation was taken from the
  # scaled block rather than the raw one.
  varying_t0 <- which(scales != 1)
  expect_equal(as.numeric(drop(stan_pop$popsd))[varying_t0],
    (rawsd * scales)[varying_t0], tolerance = 0.1)
})

test_that("every covmattransform means the same thing on both backends", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  data <- .summary_data()

  # The integer the two backends have to agree on, so a mismatch is reported as
  # the code rather than as a pile of differing arrays.
  wanted <- c(rawcorr = 0L, cholesky = 1L, z = 2L)

  for (transform in names(wanted)) {
    model <- .summary_model()
    model$covmattransform <- transform
    spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
    stan_spec <- suppressMessages(ctFit(data, model, backend = "stan",
      fit = FALSE))
    expect_equal(as.integer(spec$covmatcode), wanted[[transform]],
      info = transform)
    expect_equal(as.integer(stan_spec$standata$choleskymats),
      wanted[[transform]], info = transform)

    npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient),
      na.rm = TRUE)
    set.seed(8)
    raw <- stats::rnorm(npar, 0, .3)

    stan_pop <- suppressMessages(ctsem:::stan_constrainsamples(
      sm = ctsem:::stanmodels$ctsm, standata = stan_spec$standata,
      samples = matrix(raw, nrow = 1), cores = 1, pcovn = 10,
      dokalman = FALSE, savesubjectmatrices = FALSE))
    backend_pop <- ctsem:::.ctBackendPopArrays(
      .summary_pointfit(spec, model, raw, "julia"))

    compared <- intersect(grep("^pop_", names(stan_pop), value = TRUE),
      names(backend_pop))
    # The covariances are the arrays a construction can differ on, so require
    # them by name: a layout change that dropped one would otherwise leave the
    # loop comparing only the matrices no construction touches.
    expect_true(all(c("pop_T0cov", "pop_DIFFUSIONcov", "pop_MANIFESTcov") %in%
      compared), info = transform)
    for (name in compared) {
      expect_equal(as.numeric(backend_pop[[name]]),
        as.numeric(stan_pop[[name]]), tolerance = 1e-10,
        info = paste(transform, name))
    }
  }
})

# summary() printed `+Inf` and `NaN` in the rawpopcorr mean column on a
# degenerate fit. `popsd` and `rawpopcorr` report the spread of the
# *transformed* parameter by quadrature, so a population sd estimated large
# enough to push its parameter onto the flat part of its own transform gives a
# transformed spread of numerically zero -- and then the correlation is 0/0.
# Observed on a 6-subject, 5-random-effect fit whose optimiser walked one drift
# sd to about 13 on the raw scale.
#
# An estimate of positive infinity is worse than a blank, because a blank with
# a sentence beside it sends the reader to fit$identifiability and an Inf sends
# them nowhere. `.ctBackendMarkNoWidth` already makes the same judgement one
# column over, for a parameter whose width the curvature cannot supply.
#
# Unit and julia-free: arranging a fit that lands somewhere degenerate is
# neither cheap nor reliable, and the decision is what matters.
test_that("a reported value that is not a number is blanked and explained", {
  ordinary <- data.frame(mean = c(0.4, -0.2), sd = c(0.1, 0.2),
    `2.5%` = c(0.2, -0.6), `97.5%` = c(0.6, 0.2), check.names = FALSE,
    row.names = c("rawcor_a__b", "rawcor_c__b"))

  # Untouched, and no note: every entry is a number.
  clean <- ctsem:::.ctBackendMarkNotFinite(ordinary)
  expect_equal(clean, ordinary, ignore_attr = TRUE)
  expect_null(attr(clean, "nonfinite"))
  expect_null(ctsem:::.ctBackendNotFiniteNote(clean, "correlation"))

  # Inf and NaN both go, in every numeric column of the affected row, and NA
  # that was already there is left as it was -- a blanked width is not a
  # non-finite value and must not be counted as one.
  degenerate <- ordinary
  degenerate$mean <- c(Inf, -0.2)
  degenerate$sd <- c(NA_real_, NaN)
  marked <- ctsem:::.ctBackendMarkNotFinite(degenerate)
  expect_true(is.na(marked$mean[1]))
  expect_false(is.nan(marked$sd[2]))
  expect_true(is.na(marked$sd[2]))
  # Both rows carried a non-finite entry, so both are counted.
  expect_equal(attr(marked, "nonfinite"), 2L)
  # And the untouched row's own numbers survive.
  expect_equal(marked$mean[2], -0.2)

  note <- ctsem:::.ctBackendNotFiniteNote(marked, "correlation")
  expect_true(is.character(note))
  expect_match(note, "^2 correlations are not reported")
  expect_match(note, "fit$identifiability", fixed = TRUE)
  # Singular agreement, since a one-row note reading "1 correlations are" is
  # the kind of thing nobody fixes later.
  one <- ctsem:::.ctBackendMarkNotFinite(
    data.frame(mean = c(Inf, 0.3), row.names = c("rawcor_a__b", "rawcor_c__b")))
  expect_equal(attr(one, "nonfinite"), 1L)
  expect_match(ctsem:::.ctBackendNotFiniteNote(one, "correlation"),
    "^1 correlation is not reported")

  # A table with no numeric column at all, and an empty one: both return
  # unchanged rather than erroring, since this runs on every summary.
  expect_silent(ctsem:::.ctBackendMarkNotFinite(
    data.frame(label = "a", stringsAsFactors = FALSE)))
  expect_silent(ctsem:::.ctBackendMarkNotFinite(ordinary[0, , drop = FALSE]))
})

test_that("ctTIpredEffects reports a julia fit's own TI predictor effect (parmatrices=TRUE)", {
  skip_without_julia()
  # .summary_model() gives 'group' (TIpred 1) an effect on B1 only, and B1 is
  # the CINT[1,1] cell with the model's default meanscale=10 transform
  # (10*param), so the effect on CINT[1,1] is analytically 10*coefficient*tipred
  # -- checked exactly below, not just "runs without error".
  model <- .summary_model()
  data <- .summary_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient), na.rm = TRUE)
  raw <- rep(0.1, npar)
  fit <- .summary_pointfit(spec, model, raw, "julia")

  b1parnum <- spec$parameter_table$parnumber[spec$parameter_table$param %in% "B1"][1]
  coefrow <- spec$ti_effects[spec$ti_effects$parameter == b1parnum, ]
  expect_equal(nrow(coefrow), 1L)
  coefval <- raw[coefrow$coefficient]

  res <- suppressMessages(ctTIpredEffects(fit, parmatrices = TRUE,
    whichpars = "CINT", nsubjects = 5, whichTIpreds = 1))

  expect_named(res, c("y", "x"))
  expect_equal(dim(res$y), c(5L, 2L, 3L))
  expect_equal(dimnames(res$y)[[2]], c("CINT[1,1]", "CINT[2,1]"))
  expected <- 10 * (raw[b1parnum] + coefval * res$x[, 1])
  expect_equal(unname(res$y[, "CINT[1,1]", "Quantile0.5"]), unname(expected))
  # B2 has no TI-predictor effect in this model, so it does not move with the
  # covariate at all.
  expect_equal(unname(res$y[, "CINT[2,1]", "Quantile0.5"]),
    rep(10 * raw[spec$parameter_table$parnumber[spec$parameter_table$param %in% "B2"][1]], 5))
})

test_that("ctTIpredEffects reports a julia fit's own TI predictor effect (parmatrices=FALSE)", {
  skip_without_julia()
  model <- .summary_model()
  data <- .summary_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient), na.rm = TRUE)
  raw <- rep(0.1, npar)
  fit <- .summary_pointfit(spec, model, raw, "julia")

  b1parnum <- spec$parameter_table$parnumber[spec$parameter_table$param %in% "B1"][1]
  coefrow <- spec$ti_effects[spec$ti_effects$parameter == b1parnum, ]
  coefval <- raw[coefrow$coefficient]

  abs <- suppressMessages(ctTIpredEffects(fit, parmatrices = FALSE,
    whichpars = b1parnum, nsubjects = 4, whichTIpreds = 1))
  expected <- 10 * (raw[b1parnum] + coefval * abs$x[, 1])
  expect_equal(unname(abs$y[, 1, "Quantile0.5"]), unname(expected))

  # returndifference=TRUE subtracts the no-effect (population) value, leaving
  # exactly the covariate's own contribution.
  diff <- suppressMessages(ctTIpredEffects(fit, parmatrices = FALSE,
    whichpars = b1parnum, nsubjects = 4, whichTIpreds = 1, returndifference = TRUE))
  expect_equal(unname(diff$y[, 1, "Quantile0.5"]), unname(10 * coefval * diff$x[, 1]))
})

test_that("ctTIpredEffects on a julia fit still explains a model with no TI predictors", {
  skip_without_julia()
  nopred <- suppressWarnings(ctModel(type = "ct", n.latent = 1,
    LAMBDA = matrix(1), MANIFESTVAR = matrix(.5)))
  data <- data.frame(id = rep(1:3, each = 2), time = rep(c(0, 1), 3),
    Y1 = stats::rnorm(6))
  spec <- suppressMessages(ctFit(data, nopred, backend = "julia", fit = FALSE))
  fit <- .summary_pointfit(spec, nopred, rep(0.1, 3), "julia")

  expect_error(ctTIpredEffects(fit), "no time independent predictors")
})

test_that("a Julia model with no state-dependent cells summarises", {
  skip_without_julia()
  skip_without_julia()
  # Regression: the engine used to return the state-dependent cells as vectors
  # from the layout call, and JuliaConnectoR *hangs* -- not errors -- marshalling
  # a zero-length one. A linear model with no random effects has no such cells,
  # which makes the simplest possible model the one that deadlocked.
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
  data <- data.frame(id = rep(1:6, each = 4), time = rep(c(0, .5, 1.2, 2), 6),
    Y1 = stats::rnorm(24, 0, .5))
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  fit <- .summary_pointfit(spec, model, rep(0.1, 5), "julia")

  matrices <- ctBackendParMatrices(fit)
  expect_equal(nrow(attr(matrices, "stateDependent")), 0L)
  expect_equal(dim(matrices$DRIFT), c(1L, 1L))
  expect_true(matrices$DRIFT[1, 1] < 0)
})

test_that("state-dependent cells are named and follow the state they are given", {
  skip_without_julia()
  model <- .summary_model()
  data <- .summary_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient), na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  fit <- .summary_pointfit(spec, model, raw, "julia")

  at_default <- ctBackendParMatrices(fit)
  statedep <- attr(at_default, "stateDependent")
  # This model's DRIFT and CINT are individually varying, which is implemented
  # as a state dependence on the augmented carrier states -- so the cells must
  # be reported as conditional rather than as constants.
  expect_true(nrow(statedep) > 0)
  expect_true(all(c("DRIFT", "CINT") %in% statedep$matrix))

  augmented <- ctBackendParMatrices(fit, trim = FALSE)
  expect_equal(nrow(augmented$DRIFT), spec$nlatent_augmented)
  expect_equal(nrow(at_default$DRIFT), spec$nlatent)

  moved <- ctBackendParMatrices(fit, state = rep(1, spec$nlatent_augmented))
  expect_false(isTRUE(all.equal(at_default$DRIFT, moved$DRIFT)))
  # A cell with no state dependence must not move.
  expect_equal(at_default$LAMBDA, moved$LAMBDA)
})

# J9/F1: `ctsem_parameter_matrices` (summary_matrices.jl) is an independent,
# hand-written copy of "apply predict, then td, then update, in that order"
# -- the same block that drifted from the filter's own row 1 three times
# (kalman_filters.jl, state_sampling.jl). It is currently correct, but nothing
# ties its group order to the filter's, so a future change to one would not
# fail anything. This is the model the review report itself probed: PARS is
# in the predict group and MANIFESTVAR reads it, so this can only agree if
# predict runs before update.
#
# This test is deliberately built to fail if the summary path's group order
# regresses: reordering `apply_complex_transforms_at_indices!` in
# `ctsem_parameter_matrices` (e.g. update before predict) makes MANIFESTVAR
# read PARS before it is materialized from `state`, which is the same
# UNSET_PARAMETER-sentinel failure mode the filter's row 1 bug had -- not a
# silent near-miss.
test_that("ctBackendParMatrices runs predict before update, so an update-group cell sees a predict-group value", {
  skip_without_julia()

  t0 <- 1.5
  pars_val <- 0.4

  .m <- function(manifestvar) suppressWarnings(ctModel(
    type = "ct",
    LAMBDA = diag(1),
    PARS = matrix("mvp||TRUE", 1, 1),
    DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix(manifestvar, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(t0, 1, 1)))

  model <- .m("PARS[1,1]")
  set.seed(11)
  dat <- data.frame(id = 1:8, time = 0, Y1 = stats::rnorm(8, t0, 1))
  spec <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  fit <- .summary_pointfit(spec, model, rep(-0.5, npar), "julia")

  # `state` is given directly here (bypassing the raw-to-state materialization
  # ctFit's optimizer would do), so MANIFESTVAR's value depends only on
  # whether the predict group (which writes PARS from state[2]) ran before
  # the update group (which reads PARS[1,1] into MANIFESTVAR) -- exactly the
  # ordering this function must get right.
  matrices <- ctBackendParMatrices(fit, state = c(t0, pars_val))
  expect_equal(unname(matrices$MANIFESTVAR[1, 1]), pars_val, tolerance = 1e-10)

  # And it has to actually move with PARS, or a bug that just returned a
  # constant would pass the check above too.
  other <- ctBackendParMatrices(fit, state = c(t0, pars_val * 3))
  expect_false(isTRUE(all.equal(unname(matrices$MANIFESTVAR[1, 1]),
    unname(other$MANIFESTVAR[1, 1]))))
})

test_that("summary reports fixed effects and system matrices, with intervals only when earned", {
  skip_without_julia()
  set.seed(5)
  data <- do.call(rbind, lapply(1:30, function(i) data.frame(id = i,
    time = c(0, .5, 1.5, 2.4, 3.5), Y1 = stats::rnorm(5, 0, .5),
    Y2 = stats::rnorm(5, 0, .5))))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2, byrow = TRUE)))
  # estonly: ctFit() now finishes with ctOptimUncertainty() as the Stan path
  # does, and these assertions are about the point-estimate-only fit -- the
  # one whose summary must not print an interval it has not earned.
  #
  # `innergaptol = 0` holds the optimiser's cheap stopping rule off, and it is
  # load-bearing rather than tidying. What this test is about is what the
  # *reporting* does with a fit that ran into a direction the data does not
  # determine, and that is decided downstream by an eigenvalue of the
  # information at whatever point the fit stopped -- `.ctBackendNullMass()`,
  # against `.ctFlatDirectionRtol()`. The diffusion correlation here is on a
  # flat ray: the fit walks out along it gaining nothing, -207.01897 to five
  # decimals wherever it stops, and only how far it walked decides which side
  # of that threshold the curvature lands. With the rule on it stops at raw
  # -9.27 after 220 iterations and the direction reads as determined; with it
  # off it reaches -14.46 after 278 and reads as undetermined. Same estimate,
  # same likelihood, opposite diagnosis.
  #
  # So the rule is pinned here to hold that variable still, not because either
  # answer is wrong. The sensitivity itself is worth knowing about: two
  # detectors describe this coordinate and they can disagree --
  # `identifiability$parameters` flags it from the transform in both cases,
  # `intervalcheck$unidentified` from the curvature in only one -- which is
  # exactly the split `.ctFlatDirectionRtol()`'s comment warns about.
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE, innergaptol = 0)))

  point <- summary(fit)
  expect_s3_class(point, "summary.ctStanFit")
  expect_equal(nrow(point$popmeans), length(fit$estimate$raw))
  expect_identical(colnames(point$popmeans), "mean")
  # DIFFUSION and T0VAR are parameterisation, not result, and are dropped.
  expect_false(any(c("DIFFUSION", "T0VAR") %in% point$parmatrices$matrix))
  expect_true(all(c("DRIFT", "DIFFUSIONcov", "T0cov", "dtDRIFT") %in%
      point$parmatrices$matrix))
  expect_true(grepl("point estimates only", point$uncertaintyNote))
  expect_null(point$nsamples)

  # A summary must not print a zero-width interval as if it were an interval.
  expect_false("2.5%" %in% colnames(point$parmatrices))

  uncertain <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(fit, uncertainty = "hessian", finishsamples = 100, verbose = 0)))
  interval <- summary(uncertain)
  expect_identical(colnames(interval$popmeans), c("mean", "sd", "2.5%", "50%", "97.5%"))
  expect_true(all(c("Mean", "sd", "2.5%", "50%", "97.5%") %in% colnames(interval$parmatrices)))
  # `ndraws`, not `nsamples`: an optimised fit never sampled, and reporting a
  # sample count for uncertainty draws read as MCMC output.
  expect_equal(interval$ndraws, 100)
  expect_null(interval$nsamples)
  # The intervals come from the draws, so they must have width -- except along
  # a direction the data does not determine, which is not inverted at all and
  # so contributes none. This model, fitted to noise, drives its diffusion
  # correlation onto the boundary and gives exactly one such direction. The
  # covariance used to floor that eigenvalue at `ridge` instead, which put a
  # standard error of 1e4 on the coordinate and, through the draws, an interval
  # of [-1, 1] on a correlation estimated at -1 -- width invented by the ridge.
  #
  # NA rather than zero, which is what this asserted before. A zero-width
  # interval is still an interval to read, and on a coordinate that only
  # partly lies in the flat direction the same projection gives a *small*
  # width instead of none -- measured elsewhere as sd 0.009 and z 65.3 on a
  # correlation whose profile likelihood is bit-identical from 0.60 to 0.998.
  # There is no width that reports that honestly, so none is reported; the
  # note and `fit$uncertainty$intervalcheck` say which coordinates and why.
  flat <- uncertain$identifiability$parameters
  expect_length(flat, 1L)
  widths <- stats::setNames(interval$popmeans[, "97.5%"] -
      interval$popmeans[, "2.5%"], rownames(interval$popmeans))
  expect_true(all(widths[setdiff(names(widths), flat)] > 0))
  expect_true(is.na(unname(widths[flat])))
  expect_true(is.na(interval$popmeans[flat, "sd"]))
  expect_equal(uncertain$uncertainty$intervalcheck$unidentified, flat)
  # And the reader is told, in the note that is always there rather than in a
  # section that comes and goes.
  expect_match(interval$uncertaintyNote, "No curvature at the estimate along")
  expect_match(interval$uncertaintyNote, flat, fixed = TRUE)

  expect_output(print(interval), "System Matrices")
  expect_output(print(interval), "Fixed-effects")
})

test_that("a default julia fit carries uncertainty, as an optimized Stan fit does", {
  skip_without_julia()
  set.seed(5)
  data <- do.call(rbind, lapply(1:30, function(i) data.frame(id = i,
    time = c(0, .5, 1.5, 2.4, 3.5), Y1 = stats::rnorm(5, 0, .5),
    Y2 = stats::rnorm(5, 0, .5))))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2, byrow = TRUE)))
  # No optimcontrol: the point of this test is what a user gets by default.
  fit <- suppressWarnings(suppressMessages(
    ctFit(data, model, backend = "julia", verbose = 0,
      optimcontrol = list(finishsamples = 50))))

  expect_equal(nrow(fit$estimate$rawposterior), 50L)
  expect_equal(fit$uncertainty$settings$method, "hessian")
  expect_true(all(is.finite(fit$estimate$se)))

  out <- summary(fit, parmatrices = FALSE)
  expect_identical(colnames(out$popmeans), c("mean", "sd", "2.5%", "50%", "97.5%"))
  expect_false(is.null(out$residCovStd))
  expect_false(is.null(out$logposterior))
  # The filter output summary() reads for that is cached at fit time, as the
  # Only the prior prediction errors are cached now, not a whole filter pass:
  # every summary that reads one reads exactly `errprior`, and on the julia
  # backend the rest is a bridge transfer rather than a memory cost.
  expect_false(is.null(fit$priorerrors))
  expect_null(fit$kalman)
  expect_equal(dim(fit$priorerrors), dim(ctsem:::.ctFitObservedY(fit)))
  # The narrow path must agree with the filter it replaces, exactly. An
  # off-by-one in the prior/updated/smoothed stacking would otherwise show up
  # only as a quietly wrong residual covariance.
  full <- suppressMessages(ctKalmanArray(fit, pointest = TRUE))
  expect_equal(as.numeric(fit$priorerrors), as.numeric(full$errprior))
  expect_true(nzchar(attr(fit$priorerrors, "conditioning")))
})

# An actual OU process, so the variance parameters sit in the interior. Fitting
# a continuous-time process model to white noise drives DIFFUSION to its zero
# boundary, which is the right answer for that data but a poor place to check
# that a summary reports transformed values.
.summary_ou_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8
  diffusion <- 0.5
  do.call(rbind, lapply(seq_len(40), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
}

test_that("summary reports transformed values, not the raw parameters", {
  skip_without_julia()
  data <- .summary_ou_data()
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  popmeans <- summary(fit)$popmeans
  expect_identical(rownames(popmeans), c("drift", "diff", "mvar", "mmean", "t0v"))
  # The transform is what makes these summaries worth having: drift is
  # -log1p_exp(raw), so it is negative whatever the raw value is, and the
  # variance parameters are positive. A summary that reported the raw vector
  # would fail both.
  expect_true(popmeans["drift", "mean"] < 0)
  expect_true(all(popmeans[c("diff", "mvar", "t0v"), "mean"] > 0))
  expect_false(isTRUE(all.equal(as.numeric(popmeans[, "mean"]), fit$estimate$raw)))
  # And it is the value the engine actually put in the matrix.
  expect_equal(popmeans["drift", "mean"],
    round(ctBackendParMatrices(fit)$DRIFT[1, 1], 3))
})

test_that("ctSummaryMatrices and ctDiscretePars work on backend fits", {
  skip_without_julia()
  set.seed(5)
  data <- do.call(rbind, lapply(1:30, function(i) data.frame(id = i,
    time = c(0, .5, 1.5, 2.4, 3.5), Y1 = stats::rnorm(5, 0, .5),
    Y2 = stats::rnorm(5, 0, .5))))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2, byrow = TRUE)))
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  matrices <- ctSummaryMatrices(fit)
  expect_true(all(c("DRIFT", "DIFFUSIONcov", "T0cov", "asymDIFFUSIONcov", "dtDRIFT") %in%
      names(matrices)))
  expect_identical(dimnames(matrices$DRIFT), list(model$latentNames, model$latentNames))
  # MANIFESTVAR is dropped in favour of MANIFESTcov, as for Stan fits.
  expect_null(matrices$MANIFESTVAR)
  # The default calcfunc is the median over samples; with one sample it is that
  # sample, so this must equal the engine's own matrix rather than merely being
  # the right shape.
  expect_equal(unname(matrices$DRIFT), unname(ctBackendParMatrices(fit)$DRIFT),
    tolerance = 1e-10)

  discrete <- ctDiscretePars(fit, times = c(0, 1, 2))
  expect_equal(dim(discrete), c(1L, 1L, 3L, 2L, 2L))
  # A zero time interval is the identity: the regression of a process on itself
  # at no elapsed time.
  expect_equal(unname(drop(discrete[1, 1, 1, , ])), diag(2), tolerance = 1e-10)

  uncertain <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(fit, uncertainty = "hessian", finishsamples = 60, verbose = 0)))
  sampled <- ctDiscretePars(uncertain, times = c(0, 1), nsamples = 20)
  expect_equal(dim(sampled)[1], 20L)
})

test_that("ctExtract returns pop_* arrays sized by the posterior", {
  skip_without_julia()
  data <- .summary_ou_data()
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

  point <- ctExtract(fit)
  expect_equal(dim(point$pop_DRIFT), c(1L, 1L, 1L))

  uncertain <- suppressWarnings(suppressMessages(
    ctOptimUncertainty(fit, uncertainty = "hessian", finishsamples = 100, verbose = 0)))
  posterior <- ctExtract(uncertain)
  expect_equal(dim(posterior$pop_DRIFT), c(100L, 1L, 1L))
  expect_equal(dim(posterior$popmeans), c(100L, 5L))
  expect_equal(dim(posterior$rawpars), c(100L, 5L))
  # Subsampling is honoured rather than ignored.
  expect_equal(dim(ctExtract(uncertain, nsamples = 25)$pop_DRIFT), c(25L, 1L, 1L))
})

test_that("the Stan summary path is unchanged by the shared refactor", {
  skip_without_julia()
  # ctSummaryMatrices.ctStanFit and ctDiscretePars now route through shared
  # helpers; this is the regression guard that they still work for Stan fits.
  matrices <- ctSummaryMatrices(ctstantestfit)
  expect_true(all(c("DRIFT", "DIFFUSIONcov", "T0cov", "dtDRIFT") %in% names(matrices)))
  expect_equal(dim(matrices$DRIFT),
    rep(length(ctstantestfit$ctstanmodel$latentNames), 2))

  discrete <- ctDiscretePars(ctstantestfit, times = c(0, 1), nsamples = 10)
  expect_equal(dim(discrete)[3], 2L)
  expect_equal(unname(drop(discrete[1, 1, 1, , ])),
    diag(length(ctstantestfit$ctstanmodel$latentNames)), tolerance = 1e-10)
})

# The constrain step asks the engine for the value of a handful of parameter
# cells, several times over -- five nodes of a Gauss-Hermite quadrature per
# random-effect level, two more per TI predictor. It used to ask for *every*
# cell each time and keep the two it wanted, which was ~97% waste on a bridge
# that moves about 1 MB/s: 3.7 s of a 16.7 s fit, transferring 3.3 MB to use
# 66 KB of it. `ctsem_parameter_matrices(rows = ...)` selects engine-side now.
#
# The risk in that change is an off-by-one in the flat position arithmetic,
# which would not error -- it would return a neighbouring cell's value and be
# visible only as a wrong number in a summary. So the test is that selecting
# engine-side gives bit-for-bit what selecting in R off the full array gives.
test_that("engine-side cell selection returns exactly what full transfer did", {
  skip_without_julia()
  set.seed(11)
  data <- do.call(rbind, lapply(1:15, function(i)
    data.frame(id = i, time = 0:5, Y1 = cumsum(stats::rnorm(6)) * .5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  fit <- suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    cores = 1, optimcontrol = list(finishsamples = 20))))

  cells <- ctsem:::.ctBackendFreeParameterCells(fit)
  layout <- ctsem:::.ctBackendSummaryLayout(fit)
  samples <- ctsem:::.ctBackendRawSamples(fit)
  full <- ctsem:::.ctBackendParMatricesFlat(fit, t(samples))

  # Every cell, and then a scattered subset -- a contiguous one would pass even
  # if the offsets were wrong by a constant.
  for (selection in list(cells, cells[c(2L, 5L), , drop = FALSE],
    cells[rev(seq_len(nrow(cells))), , drop = FALSE])) {
    expect_equal(
      ctsem:::.ctBackendPopCellValues(fit, samples, selection, layout),
      ctsem:::.ctBackendPopCellsFromFlat(full, selection, layout))
  }

  # And the narrowing is real, not just harmless.
  narrow <- ctsem:::.ctBackendParMatricesFlat(fit, t(samples),
    rows = ctsem:::.ctBackendCellPositions(cells[2L, , drop = FALSE], layout))
  expect_equal(nrow(narrow), 1L)
  expect_lt(nrow(narrow), nrow(full))
})

# A coordinate with no curvature reaches the reader through the summary tables,
# not through the diagnostics object, so the mapping from raw coordinate names
# to table row names is the part that has to be right. The two vocabularies
# differ by a prefix and, on a multilevel Laplace fit, by the level: the raw
# vector says `popsd_a.study` where the section is `popsd.study` and the row is
# `a`.
test_that("no-width coordinates are matched to the rows that report them", {
  expect_equal(
    ctsem:::.ctBackendNoWidthRows("popsd_diff_eta1",
      c("drift_eta1", "diff_eta1"), "popsd_"),
    c(FALSE, TRUE))
  expect_equal(
    ctsem:::.ctBackendNoWidthRows("rawcor_diff_eta1__drift_eta1",
      c("diff_eta1__drift_eta1"), "rawcor_"),
    TRUE)
  # The level travels on the raw name and in the section heading, so it has to
  # be put back before the comparison or a multilevel fit would match nothing.
  expect_equal(
    ctsem:::.ctBackendNoWidthRows("popsd_a.study", "a", "popsd_", "study"),
    TRUE)
  expect_equal(
    ctsem:::.ctBackendNoWidthRows("popsd_a.subject", "a", "popsd_", "study"),
    FALSE)
  # A model parameter is spelled the same way in both, so no prefix.
  expect_equal(ctsem:::.ctBackendNoWidthRows("lambda", c("lambda", "drift"), ""),
    c(TRUE, FALSE))
  expect_equal(ctsem:::.ctBackendNoWidthRows(character(), c("a", "b"), "popsd_"),
    c(FALSE, FALSE))
})

test_that("marking a row blanks its width and keeps its estimate", {
  table <- data.frame(mean = c(0.6, -0.4), sd = c(0.009, 0.02),
    `2.5%` = c(0.579, -0.44), `50%` = c(0.597, -0.4),
    `97.5%` = c(0.614, -0.36), z = c(65.3, -20),
    row.names = c("diff_eta1__drift_eta1", "drift_eta1"), check.names = FALSE)
  marked <- ctsem:::.ctBackendMarkNoWidth(table, c(TRUE, FALSE))
  # The estimate stays: it is an arbitrary point on the ridge, but it is what
  # the optimiser returned and there is nothing else to print. The sd, the
  # interval and the z are the ones claiming something the data does not say.
  expect_equal(marked$mean, c(0.6, -0.4))
  expect_true(all(is.na(unlist(marked[1, c("sd", "2.5%", "50%", "97.5%", "z")]))))
  # And the identified row is untouched, column for column.
  expect_equal(marked[2, ], table[2, ])
  # Nothing flagged, nothing changed.
  expect_equal(ctsem:::.ctBackendMarkNoWidth(table, c(FALSE, FALSE)), table)
})
