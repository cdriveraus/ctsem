# Discrete-time models for backend='julia'.
#
# Discrete time is not a second filter in these engines, only a different
# discretization: a discrete model's DRIFT, CINT and DIFFUSION are already the
# one-step quantities, so the matrix exponential, the Lyapunov solve and the
# discrete-intercept solve all collapse. Everything downstream -- the
# measurement update, the smoother, subject parameters, generation -- is the
# same code.
#
# The comparison against Stan is made at a *fixed* raw parameter vector rather
# than between two fits, so the optimizer plays no part in it.

.discrete_model <- function() {
  suppressWarnings(ctModel(type = "dt", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    DRIFT = matrix(c("a11", "a12", "a21", "a22"), 2, 2, byrow = TRUE)))
}

.discrete_data <- function(times = 0:4) {
  set.seed(3)
  data <- do.call(rbind, lapply(1:10, function(i) data.frame(id = i, time = times,
    Y1 = stats::rnorm(length(times), 0, .5), Y2 = stats::rnorm(length(times), 0, .5))))
  data$Y2[3] <- NA
  data[10, c("Y1", "Y2")] <- NA
  data
}

.discrete_raw <- function(spec) {
  set.seed(8)
  stats::rnorm(max(spec$parameter_table$parnumber, na.rm = TRUE), 0, .3)
}

test_that("discrete-time likelihood and gradient match Stan", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .discrete_model()
  data <- .discrete_data()
  expect_false(model$continuoustime)

  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  expect_false(spec$continuoustime)
  raw <- .discrete_raw(spec)

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE))
  stan <- rstan::log_prob(ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm,
    stan_spec$standata), upars = raw, adjust_transform = FALSE, gradient = TRUE)
  julia <- ctJuliaEvaluate(spec, raw, gradient = TRUE)

  expect_equal(as.numeric(julia$value), as.numeric(stan), tolerance = 1e-8)
  expect_equal(as.numeric(julia$gradient), as.numeric(attributes(stan)$gradient),
    tolerance = 1e-8)

  # The adjoint's discrete branch against central differences of its own
  # likelihood: an independent check that does not share code with it.
  finite <- vapply(seq_along(raw), function(i) {
    step <- rep(0, length(raw))
    step[i] <- 1e-5
    (ctJuliaEvaluate(spec, raw + step, gradient = FALSE)$value -
        ctJuliaEvaluate(spec, raw - step, gradient = FALSE)$value) / 2e-5
  }, numeric(1))
  expect_equal(as.numeric(julia$gradient), finite, tolerance = 1e-6)
})

test_that("in discrete time the recorded intervals do not matter", {
  skip_without_julia()
  # The property that distinguishes a discrete model from a continuous one: each
  # row advances exactly one step whatever the interval says. This is also the
  # check that `dt` really has been taken out of the discretization rather than
  # left in somewhere.
  model <- .discrete_model()
  even <- suppressMessages(ctFit(.discrete_data(0:4), model, backend = "julia",
    fit = FALSE))
  uneven <- suppressMessages(ctFit(.discrete_data(c(0, 2.5, 7, 9, 20)), model,
    backend = "julia", fit = FALSE))
  raw <- .discrete_raw(even)

  expect_equal(ctJuliaEvaluate(even, raw, gradient = TRUE)$value,
    ctJuliaEvaluate(uneven, raw, gradient = TRUE)$value, tolerance = 1e-12)
  expect_equal(ctJuliaEvaluate(even, raw, gradient = TRUE)$gradient,
    ctJuliaEvaluate(uneven, raw, gradient = TRUE)$gradient, tolerance = 1e-10)

  # And a continuous reading of the same data does depend on them, so the
  # comparison above is not vacuous.
  continuous <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1)))
  even_ct <- suppressMessages(ctFit(.discrete_data(0:4), continuous, backend = "julia",
    fit = FALSE))
  uneven_ct <- suppressMessages(ctFit(.discrete_data(c(0, 2.5, 7, 9, 20)), continuous,
    backend = "julia", fit = FALSE))
  rawct <- .discrete_raw(even_ct)
  expect_false(isTRUE(all.equal(ctJuliaEvaluate(even_ct, rawct, gradient = FALSE)$value,
    ctJuliaEvaluate(uneven_ct, rawct, gradient = FALSE)$value)))
})

test_that("discrete-time prediction matches Stan and the asymptotics are discrete", {
  skip_if_not_installed("rstan")
  skip_without_julia()
  model <- .discrete_model()
  data <- .discrete_data()
  spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  raw <- .discrete_raw(spec)
  asmodel <- ctsem:::.ctBackendAsModel(spec)

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE))
  stan <- suppressMessages(ctsem:::stan_constrainsamples(sm = ctsem:::stanmodels$ctsm,
    standata = stan_spec$standata, samples = matrix(raw, nrow = 1), cores = 1,
    pcovn = 5, savescores = TRUE, savesubjectmatrices = TRUE))
  traced <- ctsem:::.ctBackendKalmanRaw(asmodel, raw)
  for (kind in 1:3) {
    expect_equal(traced$eta[kind, , ], stan$etaa[1, kind, , ], tolerance = 1e-7,
      info = kind)
    expect_equal(traced$etacov[kind, , , ], stan$etacova[1, kind, , , ],
      tolerance = 1e-7, info = kind)
  }
  expect_equal(as.numeric(traced$llrow), as.numeric(stan$llrow[1, ]), tolerance = 1e-7)

  # The interval transition is JAx itself, not an exponential of it.
  matrices <- ctBackendParMatrices(asmodel, raw, trim = FALSE)
  starts <- spec$subject_starts
  for (row in setdiff(seq_along(spec$times), starts)) {
    expect_equal(traced$transition[row, , ], unname(matrices$JAx), tolerance = 1e-12,
      info = row)
  }

  # And the asymptotic forms are the discrete ones: (I - A) x = c, and
  # X = A X A' + Q rather than the continuous Lyapunov equation.
  A <- unname(matrices$DRIFT)
  expect_equal((diag(2) - A) %*% unname(matrices$asymCINT), unname(matrices$CINT),
    tolerance = 1e-9)
  X <- unname(matrices$asymDIFFUSIONcov)
  expect_equal(X, A %*% X %*% t(A) + unname(matrices$DIFFUSIONcov), tolerance = 1e-9)
})

test_that("summary, ctKalmanArray and generation work in discrete time", {
  skip_without_julia()
  model <- .discrete_model()
  data <- .discrete_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))

  summarised <- summary(fit)
  expect_true("DRIFT" %in% summarised$parmatrices$matrix)
  # dtDRIFT is a continuous-time construction -- there is nothing to
  # discretise -- so it must be absent, as it is for a Stan discrete fit.
  expect_false("dtDRIFT" %in% summarised$parmatrices$matrix)
  expect_true(all(c("asymCINT", "asymDIFFUSIONcov") %in% summarised$parmatrices$matrix))

  kalman <- suppressMessages(ctKalmanArray(fit, subjects = "all",
    standardisederrors = TRUE))
  expect_equal(dim(kalman$etaprior)[2], length(fit$model_spec$times))
  expect_true(all(is.finite(kalman$etasmooth)))

  # Interpolation is meaningless for a discrete model and is refused rather
  # than silently producing steps at fractional times.
  expect_error(ctKalmanArray(fit, subjects = "all", timestep = .5),
    "Discrete time")

  set.seed(5)
  generated <- ctGenerateFromFit(fit, nsamples = 5, cores = 1)
  expect_equal(dim(generated$generated$Y), c(5L, length(fit$model_spec$times), 2L))
  # Missingness is preserved here too.
  observed <- t(fit$model_spec$manifest_data)
  expect_equal(which(is.na(generated$generated$Y[1, , 1])),
    which(!is.finite(observed[, 1])))
})
