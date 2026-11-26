# ctFit(intoverpop = 'variational'): each unit's random effects integrated by
# a lower bound fitted at every evaluation (variational.jl in the engine,
# `.ctBackendVariational()`).
#
# What is pinned here: on a model where every unit's integrand is Gaussian the
# bound is the integral, so the fit is Laplace's (exact there); on a
# nonlinear one the fit reports the importance-sampled log likelihood with
# the bound below it, and its estimates are those of the refined fit, which
# integrates every unit numerically; and the route is refused where it cannot
# apply.

.vf_linear_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.vf_linear_data <- function(nsubjects = 25, nobs = 6) {
  set.seed(20260827)
  drift <- -0.4; diffusion <- 0.6
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1.5, 0.9)
    state <- stats::rnorm(1, 0, 0.5)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        decay <- exp(drift)
        state <- decay * state +
          stats::rnorm(1, 0, sqrt(diffusion^2 / (-2 * drift) * (1 - decay^2)))
      }
      out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
}

# One random effect, on a `-log1p_exp` DRIFT, simulated here rather than by
# ctGenerate, whose draw stream moves.
.vf_nonlinear_model <- function() {
  model <- suppressMessages(ctModel(silent = TRUE, type = "ct", CINT = 0,
    MANIFESTMEANS = 0, LAMBDA = matrix(1), T0MEANS = matrix(0),
    DRIFT = "drift|-log1p_exp(-param)|TRUE"))
  model$pars$indvarying <- model$pars$param %in% "drift"
  model
}

.vf_nonlinear_data <- function(seed = 3L, nsubjects = 40L, ntimes = 8L) {
  set.seed(seed)
  drift <- -log1p(exp(-stats::rnorm(nsubjects, 1, 1)))
  rows <- lapply(seq_len(nsubjects), function(i) {
    a <- drift[i]; decay <- exp(a)
    innovation <- sqrt(0.25 * (exp(2 * a) - 1) / (2 * a))
    latent <- numeric(ntimes); latent[1] <- stats::rnorm(1, 0, 1)
    for (t in seq_len(ntimes - 1L)) latent[t + 1L] <- decay * latent[t] +
      stats::rnorm(1, 0, innovation)
    data.frame(id = i, time = seq_len(ntimes) - 1L,
      Y1 = latent + stats::rnorm(ntimes, 0, 0.3))
  })
  do.call(rbind, rows)
}

.vf_cache <- new.env(parent = emptyenv())
.vf_fit <- function(which, route) {
  key <- paste(which, route)
  if (!exists(key, envir = .vf_cache, inherits = FALSE)) {
    data <- if (identical(which, "linear")) .vf_linear_data() else .vf_nonlinear_data()
    model <- if (identical(which, "linear")) .vf_linear_model() else .vf_nonlinear_model()
    set.seed(1)
    assign(key, suppressWarnings(suppressMessages(ctFit(data, model,
      backend = "julia", intoverpop = route, cores = 2))), envir = .vf_cache)
  }
  get(key, envir = .vf_cache)
}

test_that("on Gaussian units the variational fit is Laplace's", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("linear", "variational")
  lf <- .vf_fit("linear", "laplace")
  expect_identical(vf$args$resolved$intoverpop, "variational")
  expect_identical(vf$estimate$loglik_method, "variational")
  expect_equal(vf$estimate$raw, lf$estimate$raw, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik, lf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$loglik_lowerbound, vf$estimate$loglik, tolerance = 1e-3)
  expect_equal(vf$estimate$se, lf$estimate$se, tolerance = 0.02)
})

test_that("on a nonlinear model the bound lies below the reported log likelihood", {
  skip_on_cran()
  skip_without_julia()
  vf <- .vf_fit("nonlinear", "variational")
  qf <- .vf_fit("nonlinear", "quadrature")
  e <- vf$estimate
  expect_identical(e$loglik_method, "variational")
  expect_true(is.finite(e$loglik) && is.finite(e$loglik_lowerbound))
  expect_lte(e$loglik_lowerbound, e$loglik + 4 * e$loglik_se)
  expect_equal(ctsem:::.ctFitOptimisedLoglik(vf), e$loglik_lowerbound)
  # The refined fit integrates every flagged unit numerically; the variational
  # estimate lies within a fifth of a standard error of it.
  expect_lt(max(abs(vf$estimate$raw - qf$estimate$raw) / qf$estimate$se), 0.2)
  expect_equal(e$loglik, qf$estimate$loglik, tolerance = 0.1)
})

test_that("intoverpop = 'variational' is refused where it cannot apply", {
  skip_on_cran()
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "stan", intoverpop = "variational"))),
    "needs backend = 'julia'")
  expect_error(suppressWarnings(suppressMessages(ctFit(.vf_linear_data(5, 3),
    .vf_linear_model(), backend = "julia", intoverpop = "variational",
    optimize = FALSE))), "fits by optimizing")
})
