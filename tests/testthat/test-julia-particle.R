# ctParticleLik: the particle-filter reference likelihood against a julia fit,
# and nlcontrol$transition for the state-explicit path. Kept small: one fit.

skip_without_julia()

test_that("ctParticleLik agrees with the filter on a linear fit and reports per row", {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)))
  set.seed(4)
  dat <- data.frame(id = rep(1:6, each = 5), time = rep(0:4, 6),
    Y1 = as.vector(replicate(6, cumsum(rnorm(5, 0, .5)))))
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia", cores = 1,
    inits = c(0.1, -0.2), verbose = 0)))
  out <- ctParticleLik(fit, particles = 2000, substeps = 2, seed = 2)
  expect_true(is.finite(out$loglik) && is.finite(out$se))
  expect_equal(out$fit_loglik, fit$estimate$loglik)
  # Exact transition for a linear model: the difference is Monte Carlo only.
  expect_lt(abs(out$difference), 4 * out$se + 0.1)
  expect_equal(nrow(out$rows), nrow(dat))
  expect_equal(sum(out$rows$particle), out$loglik, tolerance = 1e-8)
  expect_equal(sum(out$rows$filter), out$fit_loglik, tolerance = 1e-6)
  expect_equal(out$transition, "exponential")
  euler <- ctParticleLik(fit, particles = 1000, substeps = 10, transition = "euler", seed = 2)
  expect_lt(abs(euler$difference), 4 * euler$se + 0.3)
  expect_error(ctParticleLik(list()), "backend = 'julia'")
})

test_that("nlcontrol$transition reaches the spec and is validated", {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix(.3, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)))
  dat <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2), Y1 = 0)
  spec <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE,
    nlcontrol = list(transition = "euler", maxtimestep = 0.5)))
  expect_equal(spec$transition, "euler")
  plain <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  expect_equal(plain$transition, "exponential")
  expect_error(ctFit(dat, model, backend = "julia", fit = FALSE,
    nlcontrol = list(transition = "midpoint")), "'exponential' or 'euler'")
  expect_error(ctFit(dat, model, backend = "stan", nlcontrol = list(transition = "euler")),
    "requires backend = 'julia'")
})

test_that("ctParticleLik runs a two-latent state-dependent model with a binary indicator", {
  # Triangular DRIFT, so the process is stable at every state, with the cross
  # effect of eta1 on eta2 depending on eta2; Y1 binary, Y2 Gaussian.
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 2, manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    manifesttype = c(1L, 0L), LAMBDA = diag(2), PARS = c("cross"),
    DRIFT = matrix(c("drift11", "cross * (1 + 0.3 * eta2)", 0, -0.4), 2, 2),
    CINT = matrix(0, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(c(0, .3), 2), DIFFUSION = diag(.5, 2),
    T0VAR = diag(1, 2), T0MEANS = matrix(0, 2, 1))))
  model$pars$indvarying <- FALSE
  set.seed(5)
  dat <- data.frame(id = rep(1:6, each = 5), time = rep(0:4, 6),
    Y1 = rbinom(30, 1, 0.5),
    Y2 = as.vector(replicate(6, cumsum(rnorm(5, 0, .5)))))
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia", cores = 1,
    inits = c(-0.3, 0.2), verbose = 0, optimcontrol = list(estonly = TRUE))))
  out <- ctParticleLik(fit, particles = 1000, substeps = 4, seed = 3)
  expect_true(is.finite(out$loglik) && is.finite(out$se) && out$ess_min > 0)
  expect_equal(out$fit_loglik, fit$estimate$loglik)
  expect_equal(nrow(out$rows), nrow(dat))
  expect_equal(out$rows$subject, dat$id)
  expect_equal(sum(out$rows$particle), out$loglik, tolerance = 1e-8)
  expect_equal(sum(out$rows$filter), out$fit_loglik, tolerance = 1e-6)
  # Nothing exact to compare with here; the two transitions agree.
  euler <- ctParticleLik(fit, particles = 1000, substeps = 16, transition = "euler", seed = 3)
  expect_true(is.finite(euler$loglik))
  expect_lt(abs(out$loglik - euler$loglik), 4 * (out$se + euler$se) + 0.5)
})
