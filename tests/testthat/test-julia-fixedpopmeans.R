# Stated population means: ctPar(mean=), 'name=value|...', model$fixedpopmeans
# (R/ctFixedPopMeans.R, `CTSEMAnchors` in the engine).
#
# The independent answer is the same model with the mean free: evaluated with
# its mean at the stated value, its objective is the stated model's, and its
# gradient in the other coordinates is the stated model's gradient. That holds
# exactly, so it is checked to rounding on every route that builds the
# objective differently -- Laplace, augmented, the loading form -- before
# anything is fitted. Then one fit and what reports it, one short sampled fit,
# and generation.

.stated_data <- function(nsub = 40L, nt = 8L, seed = 11L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    times <- cumsum(c(0, stats::runif(nt - 1L, 0.5, 1.5)))
    age <- stats::rnorm(1)
    b1 <- 0.7 + 0.4 * age + stats::rnorm(1, 0, 0.8)
    b2 <- 1 + stats::rnorm(1, 0, 0.5)
    x <- matrix(0, nt, 2); x[1, ] <- stats::rnorm(2)
    dr <- c(0.5, 0.8)
    for (t in 2:nt) {
      a <- exp(-dr * (times[t] - times[t - 1]))
      x[t, ] <- a * x[t - 1, ] + stats::rnorm(2, 0, sqrt((1 - a^2) / (2 * dr)))
    }
    data.frame(id = i, time = times, Y1 = b1 + x[, 1] + stats::rnorm(nt, 0, 0.4),
      Y2 = b2 + x[, 2] + stats::rnorm(nt, 0, 0.4), age = age)
  }))
}

.stated_model <- function(mm1 = "mm1=0.7||TRUE||age", ti = TRUE) {
  suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), TIpredNames = if (ti) "age",
    LAMBDA = diag(2), DRIFT = matrix(c("d11", 0, 0, "d22"), 2),
    DIFFUSION = matrix(c("q11", 0, 0, "q22"), 2),
    MANIFESTVAR = matrix(c("e1", 0, 0, "e2"), 2),
    MANIFESTMEANS = matrix(c(mm1, if (ti) "mm2||TRUE|1|" else "mm2||TRUE"), 2),
    T0MEANS = matrix(0, 2), T0VAR = diag(2), CINT = matrix(0, 2)))
}

test_that("a stated mean is written in a cell, by name, by ctPar or on the model", {
  expect_equal(ctPar("lev", indvarying = TRUE, mean = 0), "lev=0||TRUE")
  expect_equal(ctPar("lev", mean = -0.25), "lev=-0.25|")
  expect_error(ctPar("lev", mean = "a"), "one finite number")
  for (cell in c("lev=0.5||TRUE", "lev, mean=0.5, indvarying=TRUE",
      ctPar("lev", indvarying = TRUE, mean = 0.5), "lev=0.5")) {
    m <- suppressMessages(ctModel(type = "ct", manifestNames = "Y1",
      latentNames = "eta1", LAMBDA = matrix(1), MANIFESTMEANS = matrix(cell)))
    expect_equal(m$fixedpopmeans, c(lev = 0.5), info = cell)
    expect_equal(m$pars$param[m$pars$matrix %in% "MANIFESTMEANS"], "lev",
      info = cell)
  }
  expect_error(suppressMessages(ctModel(type = "ct", manifestNames = c("Y1", "Y2"),
    latentNames = c("eta1", "eta2"), LAMBDA = diag(2),
    MANIFESTMEANS = matrix(c("lev=0||TRUE", "lev=1||TRUE"), 2))),
    "Two population means")
  # The value is on summary()'s scale: through the cell's own transform.
  m <- .stated_model()
  expect_equal(ctsem:::.ctFixedPopMeans(m)$raw, 0.07)
})

test_that("a stated mean is refused where it cannot apply", {
  dat <- .stated_data(nsub = 6L, nt = 3L)
  m <- .stated_model()
  expect_error(suppressMessages(ctFit(dat, m, backend = "stan", fit = FALSE)),
    "backend='julia'")
  still <- .stated_model("mm1||FALSE")
  still$fixedpopmeans <- c(mm1 = 0)
  expect_error(ctsem:::.ctFixedPopMeans(still), "does not vary")
  unknown <- m
  unknown$fixedpopmeans <- c(nope = 0)
  expect_error(ctsem:::.ctFixedPopMeans(unknown), "nope, which is not a parameter")
  expect_error(suppressMessages(ctGenerate(m, n = 3, Tpoints = 3, backend = "r")),
    "backend='julia'")
})

test_that("with the stated value in place, the objective is the free-mean model's", {
  skip_without_julia()
  dat <- .stated_data()
  routes <- list(laplace = list(intoverpop = "laplace"),
    augmented = list(intoverpop = "augmented"),
    loading = list(intoverpop = "laplace", poprank = 1))
  for (route in names(routes)) {
    prep <- function(m) do.call(ctFit, c(list(dat, m, backend = "julia",
      fit = FALSE, cores = 2), routes[[route]]))
    p1 <- suppressMessages(prep(.stated_model()))
    p0 <- suppressMessages(prep(.stated_model("mm1||TRUE||age")))
    n1 <- ctsem:::.ctBackendNpar(p1)
    expect_equal(n1, ctsem:::.ctBackendNpar(p0) - 1L, info = route)
    expect_equal(p1$anchors$index, n1 + 1L, info = route)
    expect_equal(p1$anchors$raw, 0.07, info = route)
    table <- as.data.frame(p0$parameter_table)
    at <- table$parnumber[match("mm1", table$param)]
    set.seed(5)
    for (k in 1:2) {
      theta <- stats::rnorm(n1, 0, 0.3)
      e1 <- ctJuliaEvaluate(p1, theta, cores = 2)
      e0 <- ctJuliaEvaluate(p0, append(theta, 0.07, after = at - 1L), cores = 2)
      expect_equal(e1$value, e0$value, tolerance = 1e-10, info = route)
      expect_equal(e1$gradient, e0$gradient[-at], tolerance = 1e-9, info = route)
    }
  }
})

test_that("a fit holds the stated mean, and reports it as stated", {
  skip_without_julia()
  dat <- .stated_data()
  fit <- suppressMessages(ctFit(dat, .stated_model(), backend = "julia",
    intoverpop = "laplace", cores = 2))
  npar <- length(fit$estimate$raw)
  # The free-mean model at the fit, mean at the stated value: stationary in
  # every other coordinate.
  p0 <- suppressMessages(ctFit(dat, .stated_model("mm1||TRUE||age"),
    backend = "julia", intoverpop = "laplace", fit = FALSE, cores = 2))
  table <- as.data.frame(p0$parameter_table)
  at <- table$parnumber[match("mm1", table$param)]
  e0 <- ctJuliaEvaluate(p0, append(fit$estimate$raw, 0.07, after = at - 1L),
    cores = 2)
  expect_equal(e0$value, fit$estimate$logposterior, tolerance = 1e-6)
  expect_lt(max(abs(e0$gradient[-at])), 0.05)
  expect_equal(dim(fit$uncertainty$hessian), c(npar, npar))

  s <- summary(fit)
  expect_false("mm1" %in% rownames(s$popmeans))
  expect_true(all(c("mm1", "mm2") %in% rownames(s$popsd)))
  expect_match(s$popNote, "not estimated: mm1 = 0.7", fixed = TRUE)
  expect_true("tip_age_mm1" %in% rownames(s$tipreds))
  sp <- ctSubjectPars(fit)
  expect_equal(dim(sp)[2L], length(unique(dat$id)))
  expect_true(all(is.finite(sp)))
  expect_silent(invisible(suppressMessages(ctKalman(fit))))
  generated <- suppressMessages(ctGenerateFromFit(fit, nsamples = 2, cores = 1,
    effects = "fitted"))$generated$Y
  expect_equal(dim(generated)[1L], 2L)
  expect_true(all(is.finite(generated[!is.na(generated)])))

  # The augmented route fits the same function: a linear Gaussian model, so
  # both integrals are exact.
  aug <- suppressMessages(ctFit(dat, .stated_model(), backend = "julia",
    intoverpop = "augmented", cores = 2))
  expect_equal(aug$estimate$loglik, fit$estimate$loglik, tolerance = 1e-3)
})

test_that("a sampled fit draws the free parameters and holds the stated mean", {
  skip_without_julia()
  dat <- .stated_data(nsub = 30L, nt = 6L)
  fit <- suppressWarnings(suppressMessages(ctFit(dat,
    .stated_model("mm1=0.7||TRUE", ti = FALSE),
    backend = "julia", optimize = FALSE, priors = TRUE, cores = 1,
    sampleControl = list(chains = 1L, warmup = 30L, draws = 30L, minESS = 0,
      saveEffects = TRUE, processes = FALSE))))
  expect_equal(ncol(fit$estimate$rawposterior), ctsem:::.ctBackendNpar(fit))
  expect_true(all(is.finite(fit$estimate$rawposterior)))
  s <- summary(fit)
  expect_false("mm1" %in% rownames(s$popmeans))
  expect_match(s$popNote, "mm1 = 0.7", fixed = TRUE)
})

test_that("generation draws persons about the stated mean", {
  skip_without_julia()
  m <- .stated_model("mm1=2||TRUE", ti = FALSE)
  set.seed(4)
  d <- suppressMessages(ctGenerate(m, n = 300, Tpoints = 2, values = list(
    popmeans = c(d11 = -1, d22 = -1, q11 = 0.1, q22 = 0.1, e1 = 0.01, e2 = 0.01,
      mm2 = 0), popsd = c(mm1 = 0.3, mm2 = 0.3))))
  persons <- attr(d, "subjectpars")
  # The spread is the one stated about the stated mean, read from its own
  # position past the free parameters.
  expect_equal(stats::sd(persons$mm1), 0.3, tolerance = 0.15)
  within <- function(x, centre) abs(mean(x) - centre) < 4 * stats::sd(x) / sqrt(length(x))
  expect_true(within(persons$mm1, 2))
  expect_true(within(persons$mm2, 0))
  expect_error(suppressMessages(ctGenerate(m, n = 3, Tpoints = 2,
    values = c(mm1 = 1))), "the model states")
})
