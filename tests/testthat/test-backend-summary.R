# Transformed-parameter summaries for backend='julia'.
#
# The load-bearing claim of this architecture is that the engines produce the
# *same* pop_* arrays Stan does, so that ctSummaryMatrices(), summary() and
# ctDiscretePars() can be one implementation rather than three. The first test
# checks exactly that, at a fixed raw parameter vector, against Stan's own
# constrain step -- which is a much sharper check than comparing two fits, since
# it removes the optimizer from the comparison entirely.
#
# pop_T0VAR is deliberately excluded from that comparison. Stan computes
# T0cov = sdcovsqrt2cov(T0VAR) and *then* rescales T0cov's indvarying T0MEANS
# rows and columns by the parameter's multiplier and meanscale, leaving T0VAR
# itself unscaled; the engines fold that scale into T0VAR, so their T0VAR is the
# one whose sdcovsqrt2cov actually equals the reported T0cov. Both routes give
# an identical T0cov, which is the quantity summaries report -- summary() drops
# T0VAR from the system matrices table for precisely this parameterisation
# reason. The test asserts the agreement on T0cov rather than papering over the
# difference on T0VAR.

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
  skip_on_cran()
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
    pcovn = 10, dokalman = FALSE, savesubjectmatrices = FALSE))

  fit <- .summary_pointfit(spec, model, raw, "julia")
  backend_pop <- ctsem:::.ctBackendPopArrays(fit)

  compared <- setdiff(intersect(grep("^pop_", names(stan_pop), value = TRUE),
    names(backend_pop)), "pop_T0VAR")
  # A model with an intoverpop augmentation, TI predictors and a state-dependent
  # DRIFT: if this list ever shrinks, the comparison below has stopped covering
  # the interesting matrices and the test has quietly weakened.
  expect_true(all(c("pop_DRIFT", "pop_DIFFUSIONcov", "pop_T0cov", "pop_asymCINT",
    "pop_asymDIFFUSIONcov", "pop_CINT", "pop_LAMBDA") %in% compared))
  for (name in compared) {
    expect_equal(dim(backend_pop[[name]]), dim(stan_pop[[name]]), info = name)
    expect_equal(as.numeric(backend_pop[[name]]), as.numeric(stan_pop[[name]]),
      tolerance = 1e-8, info = name)
  }

  # The reported T0VAR differs by parameterisation, but the covariance it stands
  # for does not, and that identity is what makes the difference harmless.
  # 1e-4 rather than machine precision: sdcovsqrt2cov's correlation constraint
  # carries a 1e-5 ridge, so the implied SD is the parameter plus that ridge.
  expect_equal(sqrt(diag(drop(backend_pop$pop_T0cov))),
    diag(drop(backend_pop$pop_T0VAR)), tolerance = 1e-4)
})

test_that("a Julia model with no state-dependent cells summarises", {
  skip_without_julia()
  skip_on_cran()
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

test_that("summary reports fixed effects and system matrices, with intervals only when earned", {
  skip_on_cran()
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
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))

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
  expect_equal(interval$nsamples, 100)
  # The intervals come from the draws, so they must have width.
  expect_true(all(interval$popmeans[, "97.5%"] > interval$popmeans[, "2.5%"]))

  expect_output(print(interval), "System Matrices")
  expect_output(print(interval), "Fixed-effects")
})

test_that("a default julia fit carries uncertainty, as an optimized Stan fit does", {
  skip_on_cran()
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
  # Stan path caches stanfit$kalman.
  expect_false(is.null(fit$kalman$errprior))
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
