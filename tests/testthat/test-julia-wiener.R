# Choice and response time indicators (manifesttype 5): the drift diffusion
# model's Wiener first-passage density as one observation of the drift
# (inst/julia/ContinuousTimeSEM/src/wiener_measurement.jl), and ctModelDDM().
#
# The engine's density is checked against an independent reference, the
# large-time series to 200 terms, written here; its gradient against central
# differences on the Laplace route, which reaches the hand-written partials,
# their second order at the sweeps' duals and the third order in the jets.

# The large-time series of Navarro & Fuss (2009), many terms: a reference
# where the decision time is not small against the boundary squared.
.wiener_reference <- function(rt, up, v, a, tau, w, K = 200) {
  t <- rt - tau
  if (up) {
    v <- -v
    w <- 1 - w
  }
  k <- seq_len(K)
  g <- pi * sum(k * exp(-k^2 * pi^2 * t / (2 * a^2)) * sin(k * pi * w))
  log(g) - 2 * log(a) - v * a * w - v^2 * t / 2
}

# Trials simulated by small Euler steps, enough for data to fit.
.wiener_data <- function(nsub = 8, ntrial = 60, a = 1.2, tau = 0.25,
    seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    v <- stats::rnorm(1, 1, 0.4) + 0.5 * (x <- sample(c(-0.5, 0.5), ntrial, TRUE))
    pos <- rep(a / 2, ntrial)
    time <- rep(0, ntrial)
    open <- rep(TRUE, ntrial)
    dt <- 1e-3
    while (any(open)) {
      pos[open] <- pos[open] + v[open] * dt + stats::rnorm(sum(open), 0, sqrt(dt))
      time[open] <- time[open] + dt
      open <- open & pos > 0 & pos < a
    }
    data.frame(id = i, time = seq_len(ntrial), x = x,
      rt = ifelse(pos >= a, 1, -1) * (time + tau))
  }))
}

.wiener_model <- function(...) suppressMessages(ctModel(type = "ct",
  manifestNames = "rt", manifesttype = 5, latentNames = "v", LAMBDA = matrix(1),
  MANIFESTVAR = matrix(0), DRIFT = matrix(0), DIFFUSION = matrix(0),
  CINT = matrix(0), T0MEANS = matrix("v0||TRUE"), T0VAR = matrix(1e-3), ...))

test_that("a choice and response time indicator gets WIENER defaults and checks", {
  w <- .wiener_model()$pars
  w <- w[w$matrix %in% "WIENER", ]
  w <- w[order(w$col), ]
  expect_identical(w$param, c("boundary_rt", "ndt_rt", NA))
  expect_equal(w$value[3], 0.5)
  expect_error(suppressMessages(ctModel(type = "ct", manifestNames = "y",
    latentNames = "v", LAMBDA = matrix(1), WIENER = c(1, 0.2, 0.5))),
    "WIENER applies")
  expect_error(.wiener_model(WIENER = c("a", "tau", 1.5)), "between 0 and 1")
  expect_error(.wiener_model(WIENER = c(-1, "tau", 0.5)), "positive")
  ctm <- list(manifestNames = "rt", manifesttype = 5L)
  expect_warning(ctsem:::.ctDataWiener(data.frame(rt = c(0.5, 0.7)), ctm),
    "upper boundary")
  expect_error(ctsem:::.ctDataWiener(data.frame(rt = c(0.5, 0)), ctm), "zero")
  expect_error(suppressMessages(ctFit(data.frame(id = 1, time = 1:3,
    rt = c(0.5, -0.6, 0.7)), .wiener_model(), backend = "stan", fit = FALSE)),
    "julia")
})

test_that("ctModelDDM writes the drift, boundary, ndt and bias cells", {
  m <- ctModelDDM(name = "rt", drift = ~ x, boundary = ~ x, bias = ~ x,
    dynamic = TRUE, sv = TRUE, random = c("drift", "ndt"))
  p <- m$pars
  cell <- function(mat, col = 1) p$param[p$matrix %in% mat & p$col == col]
  expect_identical(m$manifesttype, 5)
  expect_identical(cell("WIENER", 1), "boundary * exp(boundary_x * x)")
  expect_identical(cell("WIENER", 3), "inv_logit(bias + bias_x * x)")
  expect_identical(cell("MANIFESTMEANS"), "drift_x * x")
  expect_identical(cell("MANIFESTVAR"), "sv")
  expect_identical(cell("CINT"), "drift_cint")
  expect_true(p$indvarying[p$matrix %in% "CINT"])
  expect_true(p$indvarying[p$matrix %in% "WIENER" & p$col == 2])
  expect_false(any(p$indvarying[p$param %in% c("boundary", "bias")]))
  s <- ctModelDDM(name = "rt")$pars
  expect_identical(s$param[s$matrix %in% "T0MEANS"], "drift_mean")
  expect_equal(s$value[s$matrix %in% "WIENER" & s$col == 3], 0.5)
  expect_error(ctModelDDM(drift = ~ x:z), "interactions")
  expect_error(ctModelDDM(random = "bias"), "fixed bias")
  expect_error(ctModelDDM(LAMBDA = matrix(2)), "set by ctModelDDM")
})

test_that("the engine's Wiener density matches the large-time series", {
  skip_on_cran()
  skip_without_julia()
  ctsem:::.ctJuliaModule()
  cases <- expand.grid(rt = c(0.6, 1.1, 2.5), up = c(0, 1), v = c(-1, 0.4, 2),
    a = c(0.9, 1.6), tau = 0.2, w = c(0.4, 0.6))
  ours <- vapply(seq_len(nrow(cases)), function(i) with(cases[i, ],
    as.numeric(ctsem:::.ctJuliaEval(sprintf(
      "ContinuousTimeSEM._wiener_loglikelihood(%.17g, %.17g, (0.0, %.17g, %.17g, %.17g))",
      v, if (up == 1) rt else -rt, a, tau, w)))), numeric(1))
  ref <- vapply(seq_len(nrow(cases)), function(i) with(cases[i, ],
    .wiener_reference(rt, up == 1, v, a, tau, w)), numeric(1))
  expect_equal(ours, ref, tolerance = 1e-6)
  # Drift variability: the density integrated over a normal drift by
  # quadrature in R.
  sv <- 0.7
  z <- seq(-8, 8, length.out = 801)
  dens <- vapply(z, function(e) exp(.wiener_reference(0.9, TRUE, 0.5 + sv * e,
    1.3, 0.2, 0.45)), numeric(1))
  ref_sv <- log(sum(dens * stats::dnorm(z)) * (z[2] - z[1]))
  expect_equal(as.numeric(ctsem:::.ctJuliaEval(
    "ContinuousTimeSEM._wiener_loglikelihood(0.5, 0.9, (0.7, 1.3, 0.2, 0.45))")),
    ref_sv, tolerance = 1e-7)
})

test_that("the Wiener objective's gradient matches differences on the Laplace route", {
  skip_on_cran()
  skip_without_julia()
  dat <- .wiener_data(nsub = 4, ntrial = 30)
  m <- ctModelDDM(name = "rt", drift = ~ x, boundary = ~ x, bias = ~ 1,
    sv = TRUE, dynamic = TRUE)
  prep <- suppressMessages(ctFit(dat, m, fit = FALSE))
  np <- ctsem:::.ctBackendNpar(prep)
  set.seed(5)
  x <- stats::rnorm(np, 0, 0.3)
  g <- as.numeric(ctJuliaEvaluate(prep, pars = x, cores = 2)$gradient)
  fd <- vapply(seq_len(np), function(k) {
    h <- 1e-5
    xp <- xm <- x
    xp[k] <- xp[k] + h
    xm[k] <- xm[k] - h
    (ctJuliaEvaluate(prep, pars = xp, gradient = FALSE, cores = 2)$value -
      ctJuliaEvaluate(prep, pars = xm, gradient = FALSE, cores = 2)$value) / (2 * h)
  }, numeric(1))
  expect_equal(g, fd, tolerance = 1e-6)
})

test_that("a ctModelDDM fit recovers its parameters, predicts choices and generates", {
  skip_on_cran()
  skip_without_julia()
  dat <- .wiener_data()
  fit <- fit_cached("wiener-static", suppressWarnings(suppressMessages(ctFit(dat,
    ctModelDDM(name = "rt", drift = ~ x), cores = 2))))
  pm <- summary(fit)$popmeans
  expect_true(abs(pm["boundary", "50%"] - 1.2) < 0.2)
  expect_true(abs(pm["ndt", "50%"] - 0.25) < 0.03)
  expect_true(abs(pm["drift_x", "50%"] - 0.5) < 0.3)
  k <- suppressMessages(ctKalmanArray(fit, timestep = "asdata", subjects = 1))
  expect_setequal(unique(as.numeric(k$y)), c(0, 1))
  expect_true(all(k$yprior > 0 & k$yprior < 1))
  g <- suppressMessages(ctGenerateFromFit(fit, nsamples = 1))
  y <- as.numeric(g$generated$Y)
  expect_true(all(is.finite(y) & y != 0))
  expect_true(mean(y > 0) > 0.55 && mean(y > 0) < 0.95)
  expect_true(min(abs(y)) > 0.2)
})
