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

test_that("ctKalmanArray matches Stan through the whole R path", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()

  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  # Compare the filters, not the optimizers: run both at Stan's estimate.
  julia_fit$estimate$raw <- stan_fit$stanfit$rawest

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
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  julia_fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  julia_fit$estimate$raw <- stan_fit$stanfit$rawest

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
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
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
  for (name in compared) {
    expect_equal(dim(extracted[[name]]), dim(stan[[name]]), info = name)
    expect_equal(as.numeric(extracted[[name]]), as.numeric(stan[[name]]),
      tolerance = 1e-6, info = name)
  }

  # LAMBDA is fixed in this model, so it cannot differ between subjects; DRIFT
  # has an individually varying cross effect, so it must.
  for (i in seq_len(2)) for (j in seq_len(2)) {
    expect_equal(stats::sd(extracted$subj_LAMBDA[1, , i, j]), 0, tolerance = 1e-12)
  }
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 2, 1]) > 1e-6)
  expect_true(stats::sd(extracted$subj_DRIFT[1, , 1, 1]) < 1e-12)
})

test_that("removeObs withholds observations without withholding covariates", {
  skip_without_julia()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

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

  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
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
  plain <- suppressMessages(ctFit(.kalman_linear_data(), .kalman_linear_model(),
    backend = "julia", verbose = 0))
  expect_error(ctPredictTIP(plain, plot = FALSE), "no time independent predictors")
})

test_that("prediction warns for an intoverstates=FALSE julia fit, as Stan's does", {
  skip_without_julia()
  # The smallest intoverstates=FALSE julia fit already used in the suite
  # (tests/testthat/test-julia-intoverstates.R, "a fit over the joint density
  # runs and carries its trajectory"), reused here rather than built fresh: a
  # state-explicit fit's point estimate is the mode of the joint density of
  # states and data, not of the marginal the Kalman filter below computes, so
  # filtering through it looks like an ordinary result and is not one.
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 3, manifestNames = c("o1", "b1", "c1"),
    latentNames = c("eta1", "eta2"),
    manifesttype = c(2L, 1L, 3L), ncategories = c(4L, 0L, 0L),
    LAMBDA = matrix(c(1, 0, 1, 0, 0, 1), 3, 2, byrow = TRUE),
    DRIFT = matrix(c(-0.4, 0.25, -0.15, -0.7), 2, 2, byrow = TRUE),
    DIFFUSION = diag(1, 2), T0VAR = diag(1, 2),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTVAR = diag(0, 3),
    MANIFESTMEANS = matrix(c(0, 0, 0.8), 3, 1), Tpoints = 10)))
  rows <- model$pars$matrix %in% "THRESHOLDS" & model$pars$row %in% 1L
  model$pars$value[rows] <- c(-1, 1, 1)[model$pars$col[rows]]
  model$pars$indvarying <- FALSE
  set.seed(2)
  data <- data.frame(suppressMessages(ctGenerate(model, n.subjects = 12,
    Tpoints = 8, backend = "julia", intoverstates = FALSE)))

  fitmodel <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 2, n.manifest = 3, manifestNames = c("o1", "b1", "c1"),
    latentNames = c("eta1", "eta2"), manifesttype = c(2L, 1L, 3L),
    ncategories = c(4L, 0L, 0L),
    LAMBDA = matrix(c(1, 0, 1, 0, 0, 1), 3, 2, byrow = TRUE),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTVAR = diag(0, 3),
    MANIFESTMEANS = matrix(c(0, 0, "cmean"), 3, 1), Tpoints = 8)))
  fitmodel$pars$indvarying <- FALSE

  fit <- suppressWarnings(suppressMessages(ctFit(data, fitmodel,
    backend = "julia", intoverstates = FALSE, verbose = 0,
    optimcontrol = list(estonly = TRUE))))
  expect_identical(fit$args$resolved$intoverstates, FALSE)

  expect_warning(suppressMessages(ctPredict(fit, subjects = 1)),
    "system noise represents prior")
  expect_warning(suppressMessages(ctKalman(fit, subjects = 1)),
    "system noise represents prior")
})
