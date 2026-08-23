# Prediction and Kalman output for backend='julia' and backend='cpp'.
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

.kalman_pointfit <- function(spec, model, raw, backend = "cpp") {
  structure(list(model_spec = spec, model = model, backend = backend,
    estimate = list(raw = raw, rawposterior = matrix(raw, nrow = 1),
      loglik = NA_real_)),
    class = c(if (identical(backend, "cpp")) "ctCppFit" else "ctJuliaFit", "ctFit"))
}

.kalman_stan_scores <- function(model, data, raw) {
  spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE))
  suppressMessages(ctsem:::stan_constrainsamples(sm = ctsem:::stanmodels$ctsm,
    standata = spec$standata, samples = matrix(raw, nrow = 1), cores = 1,
    pcovn = 5, savescores = TRUE, savesubjectmatrices = TRUE))
}

test_that("C++ per-row Kalman output matches Stan's", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)

  stan <- .kalman_stan_scores(model, data, raw)
  cpp <- ctsem:::.ctBackendKalmanRaw(
    ctsem:::.ctBackendAsModel(spec, "cpp"), raw)

  # Prior, filtered and smoothed, all three, at every row -- checked separately
  # rather than as one array so a failure names which pass is wrong.
  for (kind in 1:3) {
    expect_equal(cpp$eta[kind, , ], stan$etaa[1, kind, , ], tolerance = 1e-8,
      info = kind)
    expect_equal(cpp$etacov[kind, , , ], stan$etacova[1, kind, , , ],
      tolerance = 1e-8, info = kind)
    expect_equal(cpp$y[kind, , ], stan$ya[1, kind, , ], tolerance = 1e-8, info = kind)
    expect_equal(cpp$ycov[kind, , , ], stan$ycova[1, kind, , , ], tolerance = 1e-8,
      info = kind)
  }
  expect_equal(as.numeric(cpp$llrow), as.numeric(stan$llrow[1, ]), tolerance = 1e-8)
})

test_that("tracing does not change the filter", {
  skip_on_cran()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  asmodel <- ctsem:::.ctBackendAsModel(spec, "cpp")

  plain <- ctCppEvaluate(asmodel, raw, gradient = FALSE)$value
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  # Same likelihood two ways: summed per subject, and summed per row.
  expect_equal(sum(traced$subject_loglik), plain, tolerance = 1e-12)
  expect_equal(sum(traced$llrow), plain, tolerance = 1e-10)

  # And the gradient is untouched by the recording code existing at all.
  expect_equal(ctCppEvaluate(asmodel, raw, gradient = TRUE)$value, plain,
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

test_that("Julia and C++ produce identical Kalman output", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if(!isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
    "Julia is not available.")
  skip_on_cran()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  cpp_spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE))
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  npar <- max(c(cpp_spec$parameter_table$parnumber, cpp_spec$ti_effects$coefficient),
    na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)

  cpp <- ctsem:::.ctBackendKalmanRaw(ctsem:::.ctBackendAsModel(cpp_spec, "cpp"), raw)
  julia <- ctsem:::.ctBackendKalmanRaw(ctsem:::.ctBackendAsModel(julia_spec, "julia"), raw)
  for (name in c("eta", "etacov", "y", "ycov", "llrow", "subject_loglik",
    "subject_matrices")) {
    expect_equal(as.numeric(cpp[[name]]), as.numeric(julia[[name]]),
      tolerance = 1e-11, info = name)
  }
  expect_identical(as.integer(cpp$subject), as.integer(julia$subject))
})

test_that("ctKalmanArray matches Stan through the whole R path", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()

  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  cpp_fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))
  # Compare the filters, not the optimizers: run both at Stan's estimate.
  cpp_fit$estimate$raw <- stan_fit$stanfit$rawest

  stan <- suppressMessages(ctKalmanArray(stan_fit, subjects = seq_len(10),
    standardisederrors = TRUE))
  cpp <- suppressMessages(ctKalmanArray(cpp_fit, subjects = "all",
    standardisederrors = TRUE))

  expect_identical(names(stan)[names(stan) %in% names(cpp)],
    names(cpp)[names(cpp) %in% names(stan)])
  # This model's augmented filter differs from Stan's at the 1e-7 level in the
  # likelihood itself (see test-stan-cpp-parity.R), so the states it implies
  # cannot be closer than that.
  for (name in c("yprior", "yupd", "ysmooth", "etaprior", "etaupd", "etasmooth",
    "errprior", "errupd", "errsmooth", "errstdprior", "errstdupd", "errstdsmooth")) {
    expect_equal(as.numeric(cpp[[name]]), as.numeric(stan[[name]]), tolerance = 1e-5,
      info = name)
  }
  for (name in c("ypriorcov", "yupdcov", "ysmoothcov", "etapriorcov", "etaupdcov",
    "etasmoothcov")) {
    expect_equal(as.numeric(cpp[[name]]), as.numeric(stan[[name]]), tolerance = 1e-4,
      info = name)
  }
  expect_equal(as.numeric(cpp$time), as.numeric(stan$time))
  expect_equal(as.numeric(cpp$id), as.numeric(stan$id))
  expect_equal(as.numeric(cpp$y), as.numeric(stan$y))
})

test_that("ctPredict interpolates a time grid the same way Stan does", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  stan_fit <- suppressMessages(ctFit(data, model, backend = "stan", optimize = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE, finishsamples = 10),
    cores = 1, verbose = 0))
  cpp_fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))
  cpp_fit$estimate$raw <- stan_fit$stanfit$rawest

  stan <- suppressMessages(ctPredict(stan_fit, subjects = 4, timestep = .3))
  cpp <- suppressMessages(ctPredict(cpp_fit, subjects = 4, timestep = .3))
  expect_equal(nrow(cpp), nrow(stan))

  key <- intersect(c("Element", "Time", "Row", "Col", "Subject"), names(stan))
  stan <- stan[do.call(order, as.list(stan[key])), ]
  cpp <- cpp[do.call(order, as.list(cpp[key])), ]
  expect_identical(lapply(cpp[key], as.character), lapply(stan[key], as.character))
  expect_equal(cpp$value, stan$value, tolerance = 1e-5)

  # The grid is what makes this a prediction rather than a refit: times that are
  # not in the data must be there, with the model's expectation and no
  # observation.
  expect_true(any(!round(unique(cpp$Time), 8) %in% round(data$time, 8)))
  interpolated <- cpp[cpp$Element == "y" & !round(cpp$Time, 8) %in% round(data$time, 8), ]
  expect_true(nrow(interpolated) > 0)
  expect_true(all(is.na(interpolated$value)))
  predicted <- cpp[cpp$Element == "ysmooth" &
      !round(cpp$Time, 8) %in% round(data$time, 8), ]
  expect_true(all(is.finite(predicted$value)))
})

test_that("subject matrices match Stan's, and only the varying ones vary", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  model <- .kalman_indvar_model()
  data <- .kalman_indvar_data()
  spec <- suppressMessages(ctFit(data, model, backend = "cpp", fit = FALSE))
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
  skip_on_cran()
  model <- .kalman_linear_model()
  data <- .kalman_linear_data()
  fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))

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
  skip_on_cran()
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
  fit <- suppressMessages(ctFit(data, model, backend = "cpp", verbose = 0))

  residuals <- suppressMessages(ctsem:::ctResiduals(fit))
  expect_equal(nrow(residuals), nrow(data))
  expect_true(all(c("Subject", "Time", "Y1") %in% names(residuals)))
  expect_equal(stats::sd(as.numeric(residuals[["Y1"]])), 1, tolerance = .1)
  expect_equal(mean(as.numeric(residuals[["Y1"]])), 0, tolerance = .1)
})
