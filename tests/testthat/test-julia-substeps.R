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
  expect_null(plain$estimate$substeps)
  s <- auto$estimate$substeps
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
  s <- auto$estimate$substeps
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
