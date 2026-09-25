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

# Under intoverpop = 'laplace' a subject is filtered at the population vector
# shifted by its random effects, so the path a mesh is measured along depends
# on them. The engine had no method for a Laplace objective, and every such fit
# failed at its start -- including the default route of any model with a
# random effect above the subject, where intoverpop = 'auto' is 'laplace'.
#
# The data: the nonlinear process observed around a mean that varies between
# subjects. The means are spread evenly rather than drawn, so the offsets that
# separate a path at the modes from one at zero effects are there whatever the
# seed.
.substep_laplace_data <- function(nsub = 12, nrow = 8, dt = 1, seed = 3) {
  set.seed(seed)
  means <- seq(-1.5, 1.5, length.out = nsub)
  rows <- list()
  for (s in seq_len(nsub)) {
    x <- 2 + rnorm(1)
    y <- numeric(nrow)
    for (t in seq_len(nrow)) {
      y[t] <- x + means[s] + 0.3 * rnorm(1)
      if (t == nrow) break
      for (k in 1:200) {
        x <- x - 0.6 * (1 + 0.5 * x) * x * dt / 200 + 0.5 * sqrt(dt / 200) * rnorm(1)
      }
    }
    rows[[s]] <- data.frame(id = s, time = (seq_len(nrow) - 1) * dt, Y1 = y)
  }
  do.call(rbind, rows)
}

.substep_laplace_model <- function(id = "id") {
  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = diag(1),
    PARS = matrix("drift", 1, 1), DRIFT = matrix("PARS[1,1] * (1 + 0.5 * eta1)", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix("mm", 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(2, 1, 1), id = id))
  model$pars$indvarying <- model$pars$param %in% "mm"
  model
}

# `estonly`: nothing here reads a standard error, and the Laplace Hessian is
# most of what a fit this size costs.
.substep_laplace_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- withr::with_seed(1, suppressWarnings(
      suppressMessages(ctFit(.substep_laplace_data(), .substep_laplace_model(),
        backend = "julia", intoverpop = "laplace", cores = 1, verbose = 0,
        optimcontrol = list(estonly = TRUE),
        nlcontrol = list(nsubsteps = "auto", substeptol = 0.02)))))
    cached
  }
})

test_that("intoverpop = 'laplace' chooses a mesh and fits with it", {
  fit <- .substep_laplace_fit()
  expect_false(is.null(fit$model_spec$laplace))
  s <- fit$substeps
  expect_true(s$finite)
  expect_gt(s$refined, 0L)
  expect_equal(s$intervals, 12L * 7L)
  expect_true(is.integer(fit$model_spec$max_timestep))
  expect_length(fit$model_spec$max_timestep, 12L * 8L)
  expect_true(is.finite(fit$estimate$loglik))
})

test_that("the Laplace mesh is measured with each subject at its random-effect modes", {
  fit <- .substep_laplace_fit()
  spec <- fit$model_spec
  est <- fit$estimate$raw[seq_len(ctsem:::.ctBackendNpar(spec))]
  # Each subject's raw vector at its modes: the population vector shifted by
  # its effects, as the Laplace term filters it. With no TI predictors that is
  # the same vector before their effects and after.
  module <- ctsem:::.ctJuliaModule(spec$project)
  persubject <- as.matrix(ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_subject_values(ctsem:::.ctJuliaObjective(fit),
      ctsem:::.ctJuliaNumericVector(est))))
  # Without its Laplace layer the specification is the plain filter, which
  # takes one vector for every subject: measure each subject at its own.
  plain <- spec
  plain$laplace <- NULL
  ids <- unique(spec$data$id)
  expect_equal(nrow(persubject), length(ids))
  atmodes <- integer(length(spec$max_timestep))
  for (s in seq_along(ids)) {
    rows <- spec$data$id == ids[s]
    atmodes[rows] <- ctsem:::.ctJuliaAutoSubsteps(plain,
      persubject[s, ])$spec$max_timestep[rows]
  }
  chosen <- ctsem:::.ctJuliaAutoSubsteps(spec, est)$spec$max_timestep
  expect_identical(chosen, atmodes)
  # And the modes decided it. At zero effects the filter moves each latent
  # state to absorb its subject's offset, and measures that path instead.
  expect_false(identical(chosen,
    ctsem:::.ctJuliaAutoSubsteps(plain, est)$spec$max_timestep))
})

test_that("a study-level effect takes the Laplace route by default, and meshes", {
  # intoverpop = 'auto' is 'laplace' whenever a parameter varies above the
  # subject, which is how the failure was reached without asking for it. Units
  # of three subjects each, where the fit above has units of one.
  dat <- .substep_laplace_data()
  dat$study <- (dat$id - 1) %/% 3 + 1
  model <- .substep_laplace_model(id = c("id", "study"))
  model$pars$indvarying_study <- model$pars$param %in% "mm"
  fit <- withr::with_seed(1, suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "julia", cores = 1, verbose = 0, optimcontrol = list(estonly = TRUE),
    nlcontrol = list(nsubsteps = "auto", substeptol = 0.02)))))
  expect_equal(ctsem:::.ctBackendIntOverPop(fit$model_spec), "laplace")
  s <- fit$substeps
  expect_true(s$finite)
  expect_gt(s$refined, 0L)
  expect_length(fit$model_spec$max_timestep, nrow(dat))
  expect_true(is.finite(fit$estimate$loglik))
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
