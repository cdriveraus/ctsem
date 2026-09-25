# nlcontrol$nsubsteps = 'auto': the julia engine chooses the number of
# prediction substeps per observation interval from how nonlinear each interval
# turns out to be, at the starting values and again at the optimum.

skip_without_julia()

.substep_linear_model <- function() {
  suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(2, 1, 1)))
}

# dx = drift (1 + 0.5 x) x dt + diffusion dW: the drift's slope moves with the
# state, so one exponential step per interval is not exact.
.substep_nonlinear_model <- function() {
  suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), PARS = matrix("drift", 1, 1),
    DRIFT = matrix("PARS[1,1] * (1 + 0.5 * eta1)", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(2, 1, 1)))
}

# Simulate the nonlinear process with a fine Euler scheme; observed with noise.
.substep_data <- function(nsub = 30, nrow = 10, dt = 1, nonlinear = TRUE, seed = 3) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(nsub)) {
    x <- 2 + rnorm(1)
    y <- numeric(nrow)
    for (t in seq_len(nrow)) {
      y[t] <- x + 0.3 * rnorm(1)
      if (t == nrow) break
      for (k in 1:200) {
        slope <- if (nonlinear) -0.6 * (1 + 0.5 * x) else -0.6
        x <- x + slope * x * dt / 200 + 0.5 * sqrt(dt / 200) * rnorm(1)
      }
    }
    rows[[s]] <- data.frame(id = s, time = (seq_len(nrow) - 1) * dt, Y1 = y)
  }
  do.call(rbind, rows)
}

test_that("nsubsteps = 'auto' is refused outside the julia backend and outside 'auto'", {
  dat <- .substep_data(nsub = 2, nrow = 3, nonlinear = FALSE)
  expect_error(ctFit(dat, .substep_linear_model(), backend = "stan",
    nlcontrol = list(nsubsteps = "auto")), "requires backend = 'julia'")
  expect_error(ctFit(dat, .substep_linear_model(), backend = "julia",
    nlcontrol = list(nsubsteps = 3)), "NULL or 'auto'")
})

test_that("a linear model keeps one step per interval and the same fit", {
  dat <- .substep_data(nonlinear = FALSE)
  inits <- c(0.1, -0.2)
  plain <- suppressMessages(ctFit(dat, .substep_linear_model(), backend = "julia",
    cores = 1, inits = inits, verbose = 0))
  auto <- suppressMessages(ctFit(dat, .substep_linear_model(), backend = "julia",
    cores = 1, inits = inits, verbose = 0, nlcontrol = list(nsubsteps = "auto")))
  expect_null(plain$substeps)
  s <- auto$substeps
  expect_true(s$finite)
  expect_equal(s$refined, 0L)
  expect_equal(s$intervals, 30L * 9L)
  expect_equal(s$max_substeps, 1L)
  expect_false(s$refit)
  # The mesh is the spec's substep slot, one entry per row.
  expect_equal(length(auto$model_spec$max_timestep), nrow(dat))
  expect_true(all(auto$model_spec$max_timestep == 1L))
  expect_equal(auto$estimate$loglik, plain$estimate$loglik, tolerance = 1e-6)
})

test_that("a nonlinear model gets refined intervals and a fit that runs", {
  dat <- .substep_data(nonlinear = TRUE)
  auto <- suppressMessages(ctFit(dat, .substep_nonlinear_model(), backend = "julia",
    cores = 1, inits = c(0.1, -0.2), verbose = 0,
    nlcontrol = list(nsubsteps = "auto", substeptol = 0.02)))
  s <- auto$substeps
  expect_true(s$finite)
  expect_gt(s$refined, 0L)
  expect_lte(s$max_substeps, 64L)
  expect_gt(s$total, s$intervals)
  expect_true(is.logical(s$refit))
  expect_equal(length(auto$model_spec$max_timestep), nrow(dat))
  expect_true(is.finite(auto$estimate$loglik))
  # A fit with the mesh already chosen reproduces itself: the mesh is data on
  # the spec, so refitting from the stored spec integrates the same way.
  expect_true(any(auto$model_spec$max_timestep > 1L))
})

# A fitted mesh goes with the rows it was chosen for --------------------------
#
# The mesh is one count per row of the fitted data. Prediction and
# cross-validation re-prepare the fit over other rows -- ctPredictTIP()'s
# pseudo-subjects, a subset of subjects, a time grid, withheld observations.
# Copying the whole mesh there failed in the engine, or handed one subject's
# counts to another; dropping it filtered at the maxtimestep rule. The
# references below are not the code under test: the fit's own mesh and filter
# where the rows are the fit's, and otherwise the engine's choice over the
# fitted rows at the estimate. It chooses one subject at a time, so a subject's
# counts there are what any rows holding that subject's data should get.

# `TI1` is for ctPredictTIP(); the simulated process ignores it.
.substep_tip_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    dat <- .substep_data(nsub = 12, nrow = 8)
    dat$TI1 <- rep(seq(-1, 1, length.out = 12), each = 8)
    model <- suppressWarnings(suppressMessages(ctModel(
      type = "ct", LAMBDA = diag(1), PARS = matrix("drift", 1, 1),
      DRIFT = matrix("PARS[1,1] * (1 + 0.5 * eta1)", 1, 1),
      DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
      MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
      T0MEANS = matrix(2, 1, 1), TIpredNames = "TI1")))
    cached <<- withr::with_seed(1, suppressWarnings(suppressMessages(ctFit(dat,
      model, backend = "julia", cores = 1, verbose = 0,
      nlcontrol = list(nsubsteps = "auto", substeptol = 0.02)))))
    cached
  }
})

.substep_estimate <- function(fit) {
  fit$estimate$raw[seq_len(ctsem:::.ctBackendNpar(fit$model_spec))]
}

.substep_llrow <- function(spec, est) {
  as.numeric(ctsem:::.ctBackendKalmanRaw(ctsem:::.ctBackendAsModel(spec), est,
    subjectmatrices = FALSE, fields = "llrow")$llrow)
}

test_that("replaced data on the fitted rows keeps the fitted mesh", {
  fit <- .substep_tip_fit()
  spec <- fit$model_spec
  expect_true(is.integer(spec$max_timestep))
  # Observations withheld, as a cross-validation fold builds its training
  # data: the rows are the fit's, so the mesh is too.
  fold <- spec$data
  fold$Y1[c(3, 20, 41)] <- NA
  expect_identical(ctsem:::.ctFitReplaceData(fit, fold)$model_spec$max_timestep,
    spec$max_timestep)
})

test_that("replaced data on fewer subjects keeps their counts and their filter", {
  fit <- .substep_tip_fit()
  spec <- fit$model_spec
  est <- .substep_estimate(fit)
  # The whole mesh was copied, one count for every fitted row, and the engine
  # refused it.
  kept <- spec$data$id %in% c(2, 5, 9)
  fewer <- ctsem:::.ctFitReplaceData(fit, spec$data[kept, ])$model_spec
  expect_identical(fewer$max_timestep, spec$max_timestep[kept])
  expect_equal(.substep_llrow(fewer, est), .substep_llrow(spec, est)[kept],
    tolerance = 1e-10)
})

test_that("replaced data on new subjects gets the mesh chosen for their rows", {
  fit <- .substep_tip_fit()
  spec <- fit$model_spec
  est <- .substep_estimate(fit)
  # As many rows as the fit had, which is the case that did not fail: each new
  # subject was given the counts of the fitted subject in its position.
  # Relabelled in reverse, subject 113 - s holds subject s's data, so it needs
  # subject s's mesh at the estimate.
  relabelled <- spec$data
  relabelled$id <- 113 - relabelled$id
  moved <- ctsem:::.ctFitReplaceData(fit, relabelled)$model_spec
  # Chosen from the maxtimestep rule, as the relabelled rows are; the fitted
  # mesh would be the engine's fallback here and not theirs.
  rule <- spec
  rule$max_timestep <- spec$substeps$floor
  chosen <- ctsem:::.ctJuliaAutoSubsteps(rule, est)$spec$max_timestep
  expect_identical(moved$max_timestep,
    unlist(lapply(12:1, function(s) chosen[spec$data$id == s])))
})

test_that("prediction filters with the fit's mesh, or one chosen for its rows", {
  fit <- .substep_tip_fit()
  spec <- fit$model_spec
  est <- .substep_estimate(fit)
  five <- spec$data$id == 5

  # One subject at its observed times: the fit's rows, and so the fit's
  # filter. Re-preparing put the mesh back to the maxtimestep rule.
  one <- ctsem:::.ctBackendKalmanSpec(fit, subjects = 5, timestep = "asdata")
  expect_identical(one$max_timestep, spec$max_timestep[five])
  expect_equal(.substep_llrow(one, est), .substep_llrow(spec, est)[five],
    tolerance = 1e-10)

  # Withholding observations keeps every row.
  withheld <- ctsem:::.ctBackendKalmanSpec(fit, removeObs = TRUE)
  expect_identical(withheld$max_timestep, spec$max_timestep)

  # An interpolation grid is rows the fit never had, and gets a mesh of its
  # own: ctPredict()'s default.
  grid <- ctsem:::.ctBackendKalmanSpec(fit, subjects = 5, timestep = "auto")
  expect_gt(length(grid$times), sum(five))
  expect_true(is.integer(grid$max_timestep))
  expect_length(grid$max_timestep, length(grid$times))
})

test_that("ctPredictTIP filters its pseudo-subjects with a mesh chosen for them", {
  fit <- .substep_tip_fit()
  # What reaches the filter. Its pseudo-subjects are rows the fit never had,
  # and they were filtered at the maxtimestep rule.
  filtered <- list()
  original <- ctsem:::.ctBackendKalmanRaw
  testthat::local_mocked_bindings(.ctBackendKalmanRaw = function(fit, raw, ...) {
    filtered[[length(filtered) + 1L]] <<- ctsem:::.ctBackendSpec(fit)
    original(fit, raw, ...)
  }, .package = "ctsem")
  for (timestep in c("asdata", "auto")) {
    filtered <- list()
    predicted <- suppressWarnings(suppressMessages(ctPredictTIP(fit,
      plot = FALSE, doDynamics = FALSE, timestep = timestep)))
    expect_gt(nrow(predicted), 0L)
    expect_true(all(is.finite(predicted$value)))
    expect_gt(length(filtered), 0L)
    for (s in filtered) {
      expect_true(is.integer(s$max_timestep))
      expect_length(s$max_timestep, length(s$times))
    }
  }
})
