# Prediction and Kalman output for backend='julia'.
#
# The design claim is that prediction rides on the pass the likelihood already
# makes, and that everything downstream of the engine is the code Stan fits
# already use. So there are two things worth testing and they are different:
#
#   * the engine produces what Stan's `dosmoother` branch produces -- checked at
#     a *fixed* raw vector against `stan_constrainsamples()`, which takes the
#     optimizer out of the comparison entirely;
#   * tracing does not perturb the filter -- checked by the traced likelihood
#     against the untraced one, which is the failure a parity check against Stan
#     would not necessarily catch, since both sides would move together.
#
# The intoverpop model below is the one that matters for subject parameters: an
# individually varying parameter is an augmented latent state, so a subject's
# value for it is only known after that subject's data has been smoothed.

.kalman_linear_model <- function() {
  suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2, byrow = TRUE)))
}

.kalman_linear_data <- function() {
  set.seed(3)
  data <- do.call(rbind, lapply(1:6, function(i) data.frame(id = i,
    time = c(0, .4, 1.1, 2.0), Y1 = stats::rnorm(4, 0, .5),
    Y2 = stats::rnorm(4, 0, .5))))
  # Partial missingness, and one row with nothing at all: the fully missing row
  # is the case the measurement update returns early from, so it is the one
  # where a recording hook placed inside the update would silently skip.
  data$Y2[3] <- NA
  data$Y1[7] <- NA
  data[10, c("Y1", "Y2")] <- NA
  data
}

.kalman_indvar_model <- function() {
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

.kalman_indvar_data <- function() {
  set.seed(5)
  data <- do.call(rbind, lapply(1:10, function(i) data.frame(id = i,
    time = c(0, .5, 1.5, 2.4), Y1 = stats::rnorm(4, 0, .5),
    Y2 = stats::rnorm(4, 0, .5), group = rep(stats::rnorm(1), 4))))
  data$Y1[7] <- NA
  data
}

.kalman_pointfit <- function(spec, model, raw, backend = "julia") {
  structure(list(model_spec = spec, model = model, backend = backend,
    estimate = list(raw = raw, rawposterior = matrix(raw, nrow = 1),
      loglik = NA_real_)),
    class = c("ctJuliaFit", "ctFit"))
}

.kalman_stan_scores <- function(model, data, raw) {
  spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE))
  suppressMessages(ctsem:::stan_constrainsamples(sm = ctsem:::stanmodels$ctsm,
    standata = spec$standata, samples = matrix(raw, nrow = 1), cores = 1,
    pcovn = 5, savescores = TRUE, savesubjectmatrices = TRUE))
}

test_that("Julia per-row Kalman output matches Stan's", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)

  stan <- .kalman_stan_scores(model, data, raw)
  julia <- ctsem:::.ctBackendKalmanRaw(
    ctsem:::.ctBackendAsModel(spec), raw)

  # Prior, filtered and smoothed, all three, at every row -- checked separately
  # rather than as one array so a failure names which pass is wrong.
  for (kind in 1:3) {
    expect_equal(julia$eta[kind, , ], stan$etaa[1, kind, , ], tolerance = 1e-8,
      info = kind)
    expect_equal(julia$etacov[kind, , , ], stan$etacova[1, kind, , , ],
      tolerance = 1e-8, info = kind)
    expect_equal(julia$y[kind, , ], stan$ya[1, kind, , ], tolerance = 1e-8, info = kind)
    expect_equal(julia$ycov[kind, , , ], stan$ycova[1, kind, , , ], tolerance = 1e-8,
      info = kind)
  }
  expect_equal(as.numeric(julia$llrow), as.numeric(stan$llrow[1, ]), tolerance = 1e-8)
})

test_that("tracing does not change the filter", {
  skip_without_julia()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  asmodel <- ctsem:::.ctBackendAsModel(spec)

  plain <- ctJuliaEvaluate(asmodel, raw, gradient = FALSE)$value
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  # Same likelihood two ways: summed per subject, and summed per row.
  expect_equal(sum(traced$subject_loglik), plain, tolerance = 1e-12)
  expect_equal(sum(traced$llrow), plain, tolerance = 1e-10)

  # And the gradient is untouched by the recording code existing at all.
  expect_equal(ctJuliaEvaluate(asmodel, raw, gradient = TRUE)$value, plain,
    tolerance = 1e-14)

  # A row with nothing observed cannot update anything, and must contribute
  # nothing -- the measurement update returns early there, which is exactly
  # where a recorder placed inside it would fail to fire.
  missingrow <- which(apply(is.na(data[, c("Y1", "Y2")]), 1, all))
  expect_length(missingrow, 1L)
  expect_equal(traced$eta[2, missingrow, ], traced$eta[1, missingrow, ])
  expect_equal(traced$etacov[2, missingrow, , ], traced$etacov[1, missingrow, , ])
  expect_equal(traced$llrow[missingrow], 0)
})

# The stan and julia fits the two whole-R-path comparisons below share. Both
# compare the filters with the julia fit placed at stan's estimate, so where
# either optimiser stopped does not enter; they were fitted twice, identically,
# and are fitted once now (`fit_cached()`, helper-julia.R).
.kalman_indvar_fits <- function() fit_cached("kalman_indvar_fits", {
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  # 'augmented' by name: this compares the filters at one raw vector, and stan
  # has only the augmented route, while 'auto' takes laplace for this model's
  # DRIFT random effect.
  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    intoverpop = "augmented"))
  # Compare the filters, not the optimizers: run both at Stan's estimate.
  julia_fit$estimate$raw <- stan_fit$stanfit$rawest
  list(stan = stan_fit, julia = julia_fit)
})

# The linear model's julia fit, which removeObs and ctPredictTIP's refusal both
# read. Fitted once, as above.
.kalman_linear_fit <- function() fit_cached("kalman_linear_fit",
  suppressMessages(ctFit(.kalman_linear_data(), .kalman_linear_model(),
    backend = "julia", verbose = 0)))

test_that("ctKalmanArray matches Stan through the whole R path", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  fits <- .kalman_indvar_fits()
  stan_fit <- fits$stan
  julia_fit <- fits$julia

  stan <- suppressMessages(ctKalmanArray(stan_fit, subjects = seq_len(10),
    standardisederrors = TRUE))
  julia <- suppressMessages(ctKalmanArray(julia_fit, subjects = "all",
    standardisederrors = TRUE))

  expect_identical(names(stan)[names(stan) %in% names(julia)],
    names(julia)[names(julia) %in% names(stan)])
  # This model's augmented filter differs from Stan's at the 1e-7 level in the
  # likelihood itself (see test-stan-julia-parity.R), so the states it implies
  # cannot be closer than that.
  for (name in c("yprior", "yupd", "ysmooth", "etaprior", "etaupd", "etasmooth",
    "errprior", "errupd", "errsmooth", "errstdprior", "errstdupd", "errstdsmooth")) {
    expect_equal(as.numeric(julia[[name]]), as.numeric(stan[[name]]), tolerance = 1e-5,
      info = name)
  }
  for (name in c("ypriorcov", "yupdcov", "ysmoothcov", "etapriorcov", "etaupdcov",
    "etasmoothcov")) {
    expect_equal(as.numeric(julia[[name]]), as.numeric(stan[[name]]), tolerance = 1e-4,
      info = name)
  }
  expect_equal(as.numeric(julia$time), as.numeric(stan$time))
  expect_equal(as.numeric(julia$id), as.numeric(stan$id))
  expect_equal(as.numeric(julia$y), as.numeric(stan$y))
})

test_that("ctPredict interpolates a time grid the same way Stan does", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  data <- .kalman_indvar_data()
  fits <- .kalman_indvar_fits()
  stan_fit <- fits$stan
  julia_fit <- fits$julia

  stan <- suppressMessages(ctPredict(stan_fit, subjects = 4, timestep = .3))
  julia <- suppressMessages(ctPredict(julia_fit, subjects = 4, timestep = .3))
  expect_equal(nrow(julia), nrow(stan))

  key <- intersect(c("Element", "Time", "Row", "Col", "Subject"), names(stan))
  stan <- stan[do.call(order, as.list(stan[key])), ]
  julia <- julia[do.call(order, as.list(julia[key])), ]
  expect_identical(lapply(julia[key], as.character), lapply(stan[key], as.character))
  expect_equal(julia$value, stan$value, tolerance = 1e-5)

  # The grid is what makes this a prediction rather than a refit: times that are
  # not in the data must be there, with the model's expectation and no
  # observation.
  expect_true(any(!round(unique(julia$Time), 8) %in% round(data$time, 8)))
  interpolated <- julia[julia$Element == "y" & !round(julia$Time, 8) %in% round(data$time, 8), ]
  expect_true(nrow(interpolated) > 0)
  expect_true(all(is.na(interpolated$value)))
  predicted <- julia[julia$Element == "ysmooth" &
      !round(julia$Time, 8) %in% round(data$time, 8), ]
  expect_true(all(is.finite(predicted$value)))
})

test_that("subject matrices match Stan's, and only the varying ones vary", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  # Augmented by name: compared with stan's subject matrices, and 'auto'
  # takes laplace for this model's random DRIFT.
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    intoverpop = "augmented"))
  npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient),
    na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  fit <- .kalman_pointfit(spec, spec$model, raw)

  extracted <- ctExtract(fit, subjectMatrices = TRUE)
  stan <- .kalman_stan_scores(model, data, raw)
  compared <- intersect(grep("^subj_", names(extracted), value = TRUE), names(stan))
  # Stan only stores the matrices that can vary, so this list is short by
  # construction -- but it must contain the ones this model varies.
  expect_true(all(c("subj_T0MEANS", "subj_DRIFT", "subj_CINT") %in% compared))

  # The matrices a carrier state is *computed from* no longer agree, and that
  # is the point rather than a tolerance to widen.
  #
  # Both backends record a subject's matrices at the end of its pass, and both
  # then overwrite T0MEANS with the smoothed state while leaving every cell
  # derived from that state where the forward pass left it -- one observation
  # short, since the transforms run at the *start* of a row from the state
  # before that row's update. Stan's writer says so in as many words:
  # "t0means updated, other pars as per final time point" (R/ctModelWriter.R).
  # So the two agreed because julia mirrored stan, which is what makes a parity
  # check blind to an inherited mistake.
  #
  # julia re-materialises at the smoothed state now. The evidence that this is
  # the right target rather than merely a different one: ctSubjectPars() on an
  # augmented fit and on the Laplace fit of the same data -- a filtered carrier
  # against a per-subject Newton solve, sharing no code -- went from a
  # regression slope of 0.926 to exactly 1. Stan is the one behind; it is
  # deprecated, and fixing it means changing the generated program.
  #
  # The shortfall is one observation, so it is worst where there are fewest:
  # this fixture has four occasions per subject.
  behind <- c("subj_DRIFT", "subj_CINT", "subj_asymDIFFUSIONcov", "subj_asymCINT")
  for (name in compared) {
    expect_equal(dim(extracted[[name]]), dim(stan[[name]]), info = name)
    if (name %in% behind) next
    expect_equal(as.numeric(extracted[[name]]), as.numeric(stan[[name]]),
      tolerance = 1e-6, info = name)
  }

  # What replaces the comparison for those: the returned object has to be
  # consistent with itself. A carrier state has no drift and no diffusion, so
  # the CINT reported beside it must be the CINT that state implies -- which is
  # exactly what was false while the two came from different rows.
  nlatent <- fit$model_spec$nlatent
  carrier <- (nlatent + 1L):fit$model_spec$nlatent_augmented
  population <- suppressMessages(ctsem:::ctBackendParMatrices(fit, trim = FALSE))
  implied <- vapply(seq_len(dim(extracted$subj_CINT)[2]), function(si) {
    state <- as.numeric(population$T0MEANS)
    state[carrier] <- extracted$subj_T0MEANS[1, si, carrier, 1]
    as.numeric(suppressMessages(ctsem:::ctBackendParMatrices(fit, filterstate = state,
      trim = FALSE))$CINT)[seq_len(nlatent)]
  }, numeric(nlatent))
  expect_equal(as.numeric(t(implied)),
    as.numeric(extracted$subj_CINT[1, , , 1]), tolerance = 1e-5)

  # LAMBDA is fixed in this model, so it cannot differ between subjects; DRIFT
  # has an individually varying cross effect, so it must.
  for (i in seq_len(2)) for (j in seq_len(2)) {
    expect_equal(stats::sd(extracted$subj_LAMBDA[1, , i, j]), 0, tolerance = 1e-12)
  }
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 2, 1]) > 1e-6)
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 1, 1]) < 1e-12)
})

test_that("on the default route only the varying subject matrices vary", {
  skip_without_julia()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  # The default route for this model's random DRIFT is laplace, which has no
  # carrier states and no stan counterpart, so what carries over from the test
  # above is its last claim: a cell with a random effect differs between
  # subjects, and a fixed one does not.
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  expect_false(is.null(spec$laplace))
  set.seed(8)
  raw <- stats::rnorm(ctsem:::.ctBackendNpar(spec), 0, .3)
  fit <- .kalman_pointfit(spec, spec$model, raw)
  extracted <- ctExtract(fit, subjectMatrices = TRUE)
  for (i in seq_len(2)) for (j in seq_len(2)) {
    expect_equal(stats::sd(extracted$subj_LAMBDA[1, , i, j]), 0, tolerance = 1e-12)
  }
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 2, 1]) > 1e-6)
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 1, 1]) < 1e-12)
})

test_that("removeObs withholds observations without withholding covariates", {
  skip_without_julia()
  fit <- .kalman_linear_fit()

  kept <- suppressMessages(ctKalmanArray(fit, subjects = "all"))
  withheld <- suppressMessages(ctKalmanArray(fit, subjects = "all", removeObs = TRUE))

  # With nothing observed the filter can only propagate: the update must be the
  # identity everywhere, and every row's likelihood contribution zero.
  expect_equal(withheld$etaupd, withheld$etaprior)
  expect_equal(sum(withheld$llrow), 0)
  expect_false(isTRUE(all.equal(kept$etaupd, kept$etaprior)))
  # The data itself is still reported, so a caller can compare prediction to
  # observation -- which is the entire point of withholding.
  expect_equal(withheld$y, kept$y)
})

test_that("standardised residuals feed ctResiduals and ctACFresiduals", {
  skip_without_julia()
  # ctResiduals() and everything on top of it (ctACFresiduals, and the residual
  # diagnostics in the tutorial) go through ctKalmanArray(standardisederrors=
  # TRUE), so they work for these backends without their own code path. The
  # check that they are *right* rather than merely present is the definition of
  # a standardised residual: at the maximum likelihood estimate of a correctly
  # specified model they have unit variance.
  set.seed(11)
  times <- seq(0, 7, 1)
  drift <- -0.8
  diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(25), function(i) {
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
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

  residuals <- suppressMessages(ctsem:::ctResiduals(fit))
  expect_equal(nrow(residuals), nrow(data))
  expect_true(all(c("Subject", "Time", "Y1") %in% names(residuals)))
  expect_equal(stats::sd(as.numeric(residuals[["Y1"]])), 1, tolerance = .1)
  expect_equal(mean(as.numeric(residuals[["Y1"]])), 0, tolerance = .1)
})

# Three things the engines now do exactly where Stan is knowingly approximate.
# Each is checked by an identity the exact version satisfies rather than by a
# comparison against the code that produced the numbers, and each is checked to
# leave the likelihood alone -- none of them is allowed to reach it.

.kalman_indvarmeans_model <- function() {
  # MANIFESTMEANS is individually varying, which is ctsem's default. That makes
  # the measurement intercept an augmented latent state, so the measurement
  # equation is linear in the augmented state: y = Jy x, exactly.
  suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1),
    CINT = matrix(0, 2, 1), DIFFUSION = diag(c(.2, .15)),
    DRIFT = matrix(c("auto1", "cross12", "cross21", "auto2"), 2, 2, byrow = TRUE)))
}

.kalman_indvarmeans_data <- function() {
  set.seed(5)
  do.call(rbind, lapply(1:8, function(i) data.frame(id = i, time = c(0, .5, 1.5, 2.4),
    Y1 = stats::rnorm(4, 0, .5), Y2 = stats::rnorm(4, 0, .5))))
}

test_that("the measurement model is re-evaluated at the updated state", {
  skip_without_julia()
  model <- .kalman_indvarmeans_model()
  data <- .kalman_indvarmeans_data()
  expect_true(all(model$pars$indvarying[model$pars$matrix == "MANIFESTMEANS"]))

  # No priors, so ctJuliaEvaluate() below is the likelihood the trace reports.
  # The default 'randomCorr' puts a prior on the two intercepts' correlation,
  # and the evaluated objective would then be a posterior.
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  asmodel <- ctsem:::.ctBackendAsModel(spec)
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  Jy <- ctBackendParMatrices(asmodel, raw, trim = FALSE)$Jy

  # y = Jy x has to hold at prior, filtered *and* smoothed. Reporting the
  # filtered observation with a pre-update measurement intercept -- which is
  # what Stan does -- breaks it at the last two.
  for (kind in 1:3) {
    implied <- t(apply(traced$eta[kind, , , drop = TRUE], 1, function(x) Jy %*% x))
    expect_equal(as.numeric(traced$y[kind, , ]), as.numeric(implied), tolerance = 1e-10,
      info = kind)
  }

  # And it is not a rounding-level difference: the intercept state moves at the
  # update, so the stale version is wrong by LAMBDA-free carrier terms.
  nrows <- dim(traced$eta)[2]
  stale <- vapply(seq_len(nrows), function(r) {
    max(abs(traced$y[2, r, ] - (Jy %*% traced$eta[1, r, ] +
        Jy[, 1:2] %*% (traced$eta[2, r, 1:2] - traced$eta[1, r, 1:2]))))
  }, numeric(1))
  expect_true(max(stale) > 1e-6)

  # None of which may touch the likelihood.
  expect_equal(sum(traced$subject_loglik),
    ctJuliaEvaluate(asmodel, raw, gradient = FALSE)$value, tolerance = 1e-12)
})

test_that("the interval transition is the Jacobian of the interval", {
  skip_without_julia()
  skip_if_not_installed("Matrix")
  model <- .kalman_indvarmeans_model()
  data <- .kalman_indvarmeans_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  asmodel <- ctsem:::.ctBackendAsModel(spec)
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  matrices <- ctBackendParMatrices(asmodel, raw, trim = FALSE)

  times <- spec$times
  starts <- spec$subject_starts
  naug <- nrow(matrices$JAx)
  # This model has no TD predictors, so Jtd is the identity and the factor is
  # invisible here; the non-identity case, which is where dropping it changes
  # the answer, is covered in the engine suite (test_kalman_trace.jl).
  Jtd <- if (is.null(matrices$Jtd)) diag(naug) else unname(matrices$Jtd)
  for (row in seq_along(times)) {
    if (row %in% starts) {
      # Nothing precedes a subject's first row.
      expect_equal(traced$transition[row, , ], diag(naug), info = row)
      next
    }
    expected <- Jtd %*%
      as.matrix(Matrix::expm(unname(matrices$JAx) * (times[row] - times[row - 1])))
    expect_equal(traced$transition[row, , ], unname(expected), tolerance = 1e-10,
      info = row)
  }
})

test_that("the improved reports differ from Stan only where Stan is approximate", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  # The divergence is deliberate, so it is asserted rather than tolerated: the
  # prior estimates and the likelihood still match Stan exactly, and only the
  # filtered and smoothed observation estimates move -- by the amount the stale
  # measurement intercept accounts for.
  model <- .kalman_indvarmeans_model()
  data <- .kalman_indvarmeans_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  asmodel <- ctsem:::.ctBackendAsModel(spec)
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  stan <- .kalman_stan_scores(model, data, raw)

  # Latent states and the likelihood are untouched -- the filter itself did not
  # change.
  for (kind in 1:3) {
    expect_equal(traced$eta[kind, , ], stan$etaa[1, kind, , ], tolerance = 1e-6,
      info = kind)
  }
  expect_equal(as.numeric(traced$llrow), as.numeric(stan$llrow[1, ]), tolerance = 1e-6)
  expect_equal(traced$y[1, , ], stan$ya[1, 1, , ], tolerance = 1e-6)

  # The reported observations do differ, and by enough to matter.
  expect_false(isTRUE(all.equal(traced$y[2, , ], stan$ya[1, 2, , ], tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(traced$y[3, , ], stan$ya[1, 3, , ], tolerance = 1e-3)))

  # Stan's filtered estimate is the prior intercept plus the updated process,
  # which is exactly the stale quantity this replaced.
  Jy <- ctBackendParMatrices(asmodel, raw, trim = FALSE)$Jy
  nlatent <- spec$nlatent
  stalecheck <- t(vapply(seq_len(dim(traced$eta)[2]), function(r) {
    as.numeric(Jy[, -(1:nlatent), drop = FALSE] %*% traced$eta[1, r, -(1:nlatent)] +
        Jy[, 1:nlatent, drop = FALSE] %*% traced$eta[2, r, 1:nlatent])
  }, numeric(dim(traced$y)[3])))
  expect_equal(stalecheck, stan$ya[1, 2, , ], tolerance = 1e-6)
})

test_that("ctPredictTIP builds its covariate grid on a backend fit", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  # ctPredictTIP predicts at chosen covariate values by constructing a dataset
  # of pseudo-subjects, one per value, and asking the fitted model for its
  # expectation. That is a *data* operation, so once the fit can be re-prepared
  # against a new data frame it works for any backend -- and its dynamics panels
  # go through ctDiscretePars with per-subject matrices, which exercises the
  # subject-parameter path as well.
  set.seed(5)
  data <- do.call(rbind, lapply(1:20, function(i) {
    group <- stats::rnorm(1)
    data.frame(id = i, time = c(0, .6, 1.4, 2.5, 3.3),
      Y1 = stats::rnorm(5, group * .5, .5), group = group)
  }))
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 1, LAMBDA = matrix(1, 1, 1),
    MANIFESTVAR = matrix(.1, 1, 1), MANIFESTMEANS = matrix("mm||FALSE", 1, 1),
    T0MEANS = matrix("t0||TRUE", 1, 1), CINT = matrix("b||TRUE", 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE))
  model$pars$group_effect[model$pars$param == "b"] <- TRUE

  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  # Compare the predictions, not the optimizers.
  julia_fit$estimate$raw <- stan_fit$stanfit$rawest

  stan <- suppressWarnings(suppressMessages(ctPredictTIP(stan_fit, tipreds = "group",
    doDynamics = FALSE, plot = FALSE, timestep = .5)))
  julia <- suppressWarnings(suppressMessages(ctPredictTIP(julia_fit, tipreds = "group",
    doDynamics = FALSE, plot = FALSE, timestep = .5)))
  expect_equal(nrow(julia), nrow(stan))

  key <- intersect(c("Element", "Time", "Row", "Col", "Subject"), names(stan))
  stan <- stan[do.call(order, as.list(stan[key])), ]
  julia <- julia[do.call(order, as.list(julia[key])), ]
  expect_identical(lapply(julia[key], as.character), lapply(stan[key], as.character))
  expect_equal(as.numeric(julia[["value"]]), as.numeric(stan[["value"]]),
    tolerance = 1e-6)
  expect_equal(as.numeric(julia[["sd"]]), as.numeric(stan[["sd"]]), tolerance = 1e-6)

  # One pseudo-subject per requested covariate value, named by it.
  expect_equal(length(unique(julia$Subject)), 3L)
  expect_true(all(grepl("^group = ", as.character(unique(julia$Subject)))))

  # The dynamics panels need per-subject matrices, so this covers that path too.
  plots <- suppressWarnings(suppressMessages(ctPredictTIP(julia_fit, tipreds = "group",
    doDynamics = TRUE, plot = TRUE, timestep = .5)))
  expect_true(all(c("Process", "Dynamics") %in% names(plots)))
  expect_true(inherits(plots$Process$Observed[[1]], "ggplot"))
  expect_true(length(plots$Dynamics$Independent) > 0)

  # A model without TI predictors is refused by name rather than failing
  # somewhere inside the grid construction.
  plain <- .kalman_linear_fit()
  expect_error(ctPredictTIP(plain, plot = FALSE), "no time independent predictors")
})

# That prediction warns for an intoverstates=FALSE julia fit, as Stan's does,
# is asserted in test-julia-intoverstates.R, on the joint-density fit that file
# already makes; it was fitted a second time here, identically, to ask it.

# How the fit represents and integrates its random effects -------------------
#
# Three methods, two representations. 'augmented' carries every varying
# parameter as a latent state; 'laplace' and 'none' both keep them as separate
# coordinates and prepare identically, differing only in whether the effects
# are then integrated out or sampled -- which is why the description lives in
# `spec$laplace` for both and why the presence of that field cannot say which
# method was used. Deriving the method from it reported every sampled fit as
# 'laplace'.

test_that('the integration method is read from the fit rather than guessed', {
  skip_without_julia()
  set.seed(83)
  generating <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    MANIFESTVAR = matrix(0.4), CINT = matrix(0), T0MEANS = matrix(0),
    T0VAR = matrix(1), MANIFESTMEANS = matrix(0))
  dat <- do.call(rbind, lapply(1:8, function(i) {
    m <- generating
    m$matrices$CINT <- matrix(stats::rnorm(1, 0, 0.4))
    d <- as.data.frame(suppressMessages(ctGenerate(m, n.subjects = 1,
      burnin = 10, Tpoints = 8, backend = 'r')))
    d$id <- i
    d
  }))
  model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    MANIFESTVAR = matrix(0.4), CINT = matrix('cint1'), T0MEANS = matrix(0),
    T0VAR = matrix(1), MANIFESTMEANS = matrix(0))
  model$pars$indvarying <- model$pars$matrix == 'CINT'

  cases <- list(
    augmented = list(intoverpop = "augmented", optimize = TRUE),
    laplace = list(intoverpop = 'laplace', optimize = TRUE),
    none = list(intoverpop = FALSE, optimize = FALSE))

  for (nm in names(cases)) {
    args <- c(list(datalong = dat, model = model, backend = 'julia', cores = 1,
      verbose = 0, fit = FALSE), cases[[nm]])
    prepared <- suppressWarnings(suppressMessages(do.call(ctFit, args)))
    spec <- ctsem:::.ctBackendSpec(prepared)
    expect_equal(ctsem:::.ctBackendIntOverPop(spec), nm, label = nm)
    # And the representation, which two of the three share.
    expect_equal(ctsem:::.ctSpecEffectsAreCoordinates(spec), nm != 'augmented',
      label = nm)
    # The structure reads the same whichever way the effects are handled.
    levels <- ctsem:::.ctFitRandomEffectLevels(prepared)
    expect_length(levels, 1L)
    expect_equal(levels[[1L]]$params, 'cint1', label = nm)
    expect_equal(levels[[1L]]$nunits, 8L, label = nm)
  }
})

# TI-predictor effects on the Laplace route -----------------------------------
#
# A Laplace subject is filtered at the population vector shifted by its random
# effects, and the filter adds the subject's TI-predictor effects on top. A
# filter handed a vector that already carries them applies them twice, and
# nothing fails: the likelihood and the estimates do not come through these
# consumers. Comparing the consumers with one another cannot see that when
# they share the helper that builds the vector, so each is checked here
# against values derived in R from the raw vector alone -- the TI effects added
# by hand, the random-effect mode in closed form (with the effect on a manifest
# mean the model is linear Gaussian in it), and the likelihood as a
# multivariate normal.

.kalman_ti_model <- function() {
  model <- suppressMessages(ctModel(type = "ct", manifestNames = "Y1",
    latentNames = "eta1", LAMBDA = matrix(1), DRIFT = matrix("drift"),
    DIFFUSION = matrix("diffusion"), MANIFESTVAR = matrix(.3),
    MANIFESTMEANS = matrix("mm"), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), TIpredNames = "TI1"))
  model$pars$indvarying <- model$pars$param %in% "mm"
  model
}

.kalman_ti_data <- function() {
  set.seed(11)
  ti <- seq(-1, 1, length.out = 10)
  do.call(rbind, lapply(seq_along(ti), function(i) data.frame(id = i,
    time = cumsum(c(0, stats::runif(5, .5, 1.5))),
    Y1 = stats::rnorm(6, 2 * ti[i], .7), TI1 = ti[i])))
}

# Subject `i`'s manifest mean and log likelihood at standardised random effect
# `u`, or at its mode when `u` is NULL. The engine supplies only population
# quantities: the transforms, at a raw vector shifted here, and the effect's sd.
.kalman_ti_reference <- function(fit, data, i, u = NULL) {
  spec <- ctsem:::.ctBackendSpec(fit)
  raw <- fit$estimate$raw
  rows <- data$id == i
  y <- data$Y1[rows]
  time <- data$time[rows]
  ti <- spec$ti_effects
  own <- raw
  own[ti$parameter] <- own[ti$parameter] + raw[ti$coefficient] * data$TI1[rows][1]
  matrices <- function(v) suppressMessages(ctBackendParMatrices(fit, raw = v))
  at <- matrices(own)
  module <- ctsem:::.ctJuliaModule(spec$project)
  sd <- sqrt(as.numeric(ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_popcov(
    ctsem:::.ctJuliaObjective(fit), ctsem:::.ctJuliaNumericVector(raw), 1L))))
  shifted <- own
  index <- spec$laplace$levels[[1L]]$re_index
  shifted[index] <- shifted[index] + sd
  slope <- as.numeric(matrices(shifted)$MANIFESTMEANS - at$MANIFESTMEANS)

  # The latent covariance at the observation times, for a one-state process
  # starting at zero with no intercept.
  a <- at$DRIFT[1, 1]
  n <- length(time)
  v <- at$T0cov[1, 1]
  for (k in seq_len(n)[-1L]) {
    growth <- exp(2 * a * (time[k] - time[k - 1L]))
    v[k] <- growth * v[k - 1L] + at$DIFFUSIONcov[1, 1] * (growth - 1) / (2 * a)
  }
  latent <- outer(seq_len(n), seq_len(n), function(j, k)
    exp(a * abs(time[k] - time[j])) * v[pmin(j, k)])
  V <- at$LAMBDA[1, 1]^2 * latent + diag(at$MANIFESTcov[1, 1], n)
  Vinv <- solve(V)
  base <- as.numeric(at$MANIFESTMEANS)
  if (is.null(u)) u <- slope * sum(Vinv %*% (y - base)) / (1 + slope^2 * sum(Vinv))
  mean <- base + slope * u
  residual <- y - mean
  list(mean = mean, loglik = -0.5 * (n * log(2 * pi) +
    as.numeric(determinant(V)$modulus) + sum(residual * (Vinv %*% residual))))
}

test_that("each Laplace filter consumer applies a subject's TI effects once", {
  skip_without_julia()
  data <- .kalman_ti_data()
  spec <- suppressMessages(ctFit(data, .kalman_ti_model(), backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  # Drift, diffusion and the mean all carry a TI1 effect, so a subject filtered
  # at twice its effects is wrong in the dynamics as well as in the level.
  expect_length(spec$ti_effects$parameter, 3L)
  raw <- c(.3, -1, .1, .2, .2, -.2, .15)
  expect_equal(ctsem:::.ctBackendNpar(spec), length(raw))
  fit <- .kalman_pointfit(spec, spec$model, raw)
  ids <- sort(unique(data$id))
  atmode <- lapply(ids, function(i) .kalman_ti_reference(fit, data, i))
  population <- lapply(ids, function(i) .kalman_ti_reference(fit, data, i, u = 0))
  firstmean <- function(ref) vapply(ref, function(x) x$mean, numeric(1))
  loglik <- function(ref) vapply(ref, function(x) x$loglik, numeric(1))

  # The first-row prior mean and each subject's summed log likelihood. On the
  # fit's own rows ctKalmanArray solves the modes in the engine. ctKalman
  # rebuilds the rows it filters, even at timestep = 'asdata', and takes the
  # modes from the fitted specification instead, which is the path every
  # prediction takes.
  fromarray <- function(...) {
    k <- suppressMessages(ctKalmanArray(fit, ...))
    list(first = as.numeric(k$yprior[1L, !duplicated(k$id), 1L]),
      loglik = as.numeric(tapply(k$llrow[1L, ], k$id, sum)))
  }
  kalman <- function(...) {
    k <- suppressMessages(ctKalman(fit, subjects = ids, ...))
    prior <- k[k$Element == "yprior", ]
    prior <- prior[order(prior$Subject, prior$Time), ]
    rows <- k$Element == "llrow"
    list(first = prior$value[!duplicated(prior$Subject)],
      loglik = as.numeric(tapply(k$value[rows], k$Subject[rows], sum,
        na.rm = TRUE)))
  }
  for (result in list(fromarray(), kalman(timestep = "asdata"),
    kalman(timestep = 0.25))) {
    expect_equal(result$first, firstmean(atmode), tolerance = 1e-6)
    expect_equal(result$loglik, loglik(atmode), tolerance = 1e-6)
  }
  for (result in list(fromarray(randomEffects = "population"),
    kalman(timestep = "asdata", randomEffects = "population"))) {
    expect_equal(result$first, firstmean(population), tolerance = 1e-8)
    expect_equal(result$loglik, loglik(population), tolerance = 1e-8)
  }

  # Generation with every deviate zero: the state starts at zero and stays
  # there, so every generated row is the subject's manifest mean. At the modes,
  # and at a given draw of the effects, on both generators.
  rowsubject <- ctsem:::.ctFitRowSubject(fit)
  zeros <- matrix(0, 1L, nrow(data))
  nz <- ctsem:::.ctBackendStateDimension(fit)
  bysubject <- function(Y) as.numeric(tapply(as.numeric(Y), rowsubject, mean))
  spread <- function(Y) max(tapply(as.numeric(Y), rowsubject,
    function(x) diff(range(x))))
  filtered <- ctsem:::.ctBackendGenerate(fit, raw, zeros)$Y
  expect_equal(bysubject(filtered), firstmean(atmode), tolerance = 1e-6)
  expect_lt(spread(filtered), 1e-8)
  states <- ctsem:::.ctBackendGenerateStates(fit, raw, numeric(nz), zeros)$Y
  expect_equal(bysubject(states), firstmean(atmode), tolerance = 1e-6)
  expect_lt(spread(states), 1e-8)
  # One effect per unit, in unit order, which is the sampler's layout.
  units <- ctsem:::.ctFitRandomEffectLevels(fit)[[1L]]$units
  effects <- seq(-1.5, 1.2, length.out = max(units))
  given <- vapply(ids, function(i) .kalman_ti_reference(fit, data, i,
    u = effects[units[i]])$mean, numeric(1))
  expect_equal(bysubject(ctsem:::.ctBackendGenerate(fit, raw, zeros,
    effects = effects)$Y), given, tolerance = 1e-8)
  expect_equal(bysubject(ctsem:::.ctBackendGenerateStates(fit, raw, numeric(nz),
    zeros, effects = effects)$Y), given, tolerance = 1e-8)

  # The effect draws behind leave-one-row-out, concentrated at the modes.
  module <- ctsem:::.ctJuliaModule(spec$project)
  draws <- ctsem:::.ctBackendJuliaValue(module$ctsem_laplace_effect_draws(
    ctsem:::.ctJuliaObjective(fit), ctsem:::.ctJuliaNumericVector(raw), 1L,
    seed = 1L, scale = 1e-9))
  expect_equal(as.numeric(tapply(as.numeric(draws$llrow), rowsubject, sum)),
    loglik(atmode), tolerance = 1e-6)
})
