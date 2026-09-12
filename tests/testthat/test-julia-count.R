# Count manifest variables: `manifesttype = 3`, Poisson with a log link.
#
# Unlike test-julia-binary.R and test-julia-ordinal.R, there is no closed-form
# single-observation check here, at this level or in the engine's own suite --
# no file under inst/julia/ContinuousTimeSEM/test mentions count, censored, or
# the quadrature functions that carry them (`_binary_moments`, `_binary_mode`,
# `_ekf_binary_update!`) at all. What this file does check: that the model
# accepts and describes the type, that the data checks fire, that the backend
# that cannot fit one refuses instead of returning a number, that the
# reverse-mode gradient is right on a real model, and that a fit recovers
# generating parameters.

.count_data <- function(nsubjects = 40, nobs = 8, nind = 2, seed = 5) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
      manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
      DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
      MANIFESTVAR = matrix(1e-6), T0VAR = matrix(1), T0MEANS = matrix(0),
      CINT = matrix(0), MANIFESTMEANS = matrix(0), Tpoints = nobs))
    one <- data.frame(ctGenerate(gen, n.subjects = 1, Tpoints = nobs,
      backend = "r"))
    one$id <- i
    one
  }))
  for (j in seq_len(nind)) {
    d[[paste0("c", j)]] <- stats::rpois(nrow(d), exp(1.2 + d$eta))
  }
  d$eta <- NULL
  d
}

.count_model <- function(nind = 2) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = nind, manifestNames = paste0("c", seq_len(nind)),
    latentNames = "eta1", manifesttype = rep(3L, nind),
    LAMBDA = matrix(1, nind, 1),
    MANIFESTMEANS = matrix("mm", nind, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, nind))))
  m$pars$indvarying <- FALSE
  m
}

# The same model with something for the Laplace route to integrate over.
# `intoverpop='laplace'` refuses a model with no varying parameter, correctly,
# so the two-route comparisons need one.
.count_model_varying <- function(nind = 2) {
  m <- .count_model(nind)
  m$pars$indvarying[m$pars$param %in% "mm"] <- TRUE
  m
}

test_that("ctModel takes manifesttype 3 and needs no categories for it", {
  m <- .count_model()
  expect_equal(unname(m$manifesttype), c(3L, 3L))
  # `ncategories` describes ordinal variables only, and a count must not be
  # made to invent one.
  expect_true(all(m$ncategories == 0L))
  expect_null(m$THRESHOLDS)
  expect_error(suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 1, n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 5L, LAMBDA = matrix(1)))), "manifesttype must be")
})

test_that("print names counts and the link they use", {
  described <- paste(utils::capture.output(print(.count_model())),
    collapse = " ")
  expect_match(described, "c1, c2 \\(count, Poisson log link\\)")
})

test_that("count data is checked against the model", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  m <- .count_model()

  negative <- d
  negative$c1[3] <- -1
  expect_error(suppressMessages(ctFit(negative, m, backend = "julia",
    fit = FALSE)), "negative values")

  fractional <- d
  fractional$c1[3] <- 1.5
  expect_error(suppressMessages(ctFit(fractional, m, backend = "julia",
    fit = FALSE)), "whole number")

  # Legal, and worth saying: a constant count identifies no rate.
  constant <- d
  constant$c1 <- 0L
  expect_warning(suppressMessages(ctFit(constant, m, backend = "julia",
    fit = FALSE)), "not identified")
})

test_that("stan refuses counts rather than treating them as Gaussian", {
  d <- .count_data(nsubjects = 5, nobs = 4)
  expect_error(suppressMessages(ctFit(d, .count_model(), backend = "stan",
    optimcontrol = list(estonly = TRUE))), "need backend=\"julia\"")
})

test_that("the adjoint matches forward mode on a count model", {
  skip_without_julia()
  d <- .count_data()
  fit <- suppressWarnings(suppressMessages(ctFit(d, .count_model(),
    backend = "julia", intoverpop = "augmented",
    optimcontrol = list(estonly = TRUE))))
  handle <- structure(fit$model_spec,
    class = c("ctJuliaModel", "ctFitModel"))
  est <- as.numeric(fit$estimate$raw)

  # Anchored at the estimate rather than drawn from zero, which is what the
  # binary and ordinal suites do. Their likelihoods are bounded, so any raw
  # draw gives a representable objective; a count's linear predictor is
  # *exponentiated*, and the same draw from zero puts the Poisson rate near
  # 1e29. Comparing two gradients there compares two numbers that have lost
  # every significant digit -- measured, it reported a 50% disagreement where
  # the same code agrees to 1e-16 anywhere the objective is representable.
  set.seed(11)
  for (trial in 1:3) {
    at <- est + stats::rnorm(length(est), 0, 0.3)
    adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

test_that("a count model recovers what generated it, and laplace does it better", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()

  # A random start, deliberately: `inits = NULL` draws `rnorm(npar, 0, .01)`
  # from wherever the session's RNG has reached, so this asks whether the fit
  # depends on where it began. It used to. Over twenty such starts the laplace
  # route reached the optimum 14 times and the augmented route 9, and eleven of
  # the augmented failures reported `converged = TRUE` at a largest gradient of
  # 2.7e4; both are 20/20 now. The two defects behind that are a count
  # intercept carrying a location parameter's `meanscale` (see
  # `.ctModelDefaultFreePar`) and an inner mode solve that reported failure for
  # reaching the limit of double precision (see `_laplace_newton_unit_mode`).
  #
  # So this test is not pinned, and should not be: a pinned start would pass
  # over either defect returning.
  fits <- lapply(c(augmented = "augmented", laplace = "laplace"), function(route) {
    suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, optimcontrol = list(estonly = TRUE))))
  })
  for (route in names(fits)) expect_true(isTRUE(fits[[route]]$estimate$converged))

  got <- vapply(fits, function(f) {
    means <- summary(f)$popmeans
    c(mm = unname(means["mm", "mean"]),
      drift = unname(means["drift_eta1", "mean"]))
  }, numeric(2))

  # One tolerance for both routes is what this used to have, and it hid the
  # result rather than testing it: the two do not have the same accuracy here
  # and the difference is the point of having both.
  #
  # The laplace route integrates the random intercept; the augmented route
  # carries it as a state through a linearised filter, which biases a nonlinear
  # model -- the same shrinkage `test-julia-multivariate-mixed.R` records for a
  # nonlinear sd, here on a Poisson intercept. Measured, both converged:
  #
  #   route      mm (truth 1.2)      drift (truth -0.4)
  #   laplace    1.242   3.5% off    -0.402   0.5% off
  #   augmented  1.716    43% off    -0.219    45% off
  #
  # So each route is asked for what it delivers, with room over the
  # measurement, and the ordering is asserted separately because that is the
  # claim worth defending: a tolerance that covered both would pass whichever
  # way round they came out.
  expect_equal(unname(got["mm", "laplace"]), 1.2, tolerance = 0.1)
  expect_equal(unname(got["drift", "laplace"]), -0.4, tolerance = 0.1)
  expect_equal(unname(got["mm", "augmented"]), 1.2, tolerance = 0.6)
  expect_equal(unname(got["drift", "augmented"]), -0.4, tolerance = 0.6)

  # The direction, which a pair of tolerances cannot state.
  expect_lt(abs(got["mm", "laplace"] - 1.2), abs(got["mm", "augmented"] - 1.2))
  expect_lt(abs(got["drift", "laplace"] + 0.4),
    abs(got["drift", "augmented"] + 0.4))
})

test_that("the two routes agree on a count model", {
  skip_without_julia()
  d <- .count_data()
  m <- .count_model_varying()
  # A COMMON, FIXED starting point, because the difference below is read as a
  # methodological gap and that reading needs both routes to have started from
  # the same place. `inits = NULL` starts each at rnorm(npar, 0, .01) from
  # whatever RNG state it inherits, so the laplace fit's start depended on how
  # much RNG the augmented fit above it had consumed -- and on this model that
  # matters far more than 0.01 suggests. From fixed starts:
  #
  #   route      zeros          rnorm sd .01     rnorm sd .3
  #   augmented  -1483.480158   -1483.480158     -1483.480158
  #   laplace    -1471.584356   -1471.584356     -1.47e14
  #
  # The augmented route is flat over all three; the laplace route reaches the
  # intended optimum from a small start, runs away from a wider one, and from
  # some sd-0.01 starts lands near -2175. That is what made this fail
  # intermittently, and it is a property of the model rather than of either
  # construction -- both recorded values below are reproduced exactly from
  # zeros, so zeros is the neutral choice and not the one that passes.
  fit <- function(route) {
    spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, fit = FALSE)))
    suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
      intoverpop = route, inits = rep(0, ctsem:::.ctBackendNpar(spec)),
      optimcontrol = list(estonly = TRUE))))
  }
  augmented <- fit("augmented")
  laplace <- fit("laplace")
  # The two routes integrate the same random effect differently, so the gap is
  # a real methodological difference and not noise: -1471.584 (laplace) against
  # -1483.480 (augmented), 11.896 apart, 0.0080 relative. `expect_equal`
  # tolerances are relative, so the old value of 1 permitted a difference of
  # 1483 -- the two could have had nothing to do with each other and passed.
  # 0.02 keeps 2.5x headroom over the measurement.
  expect_equal(as.numeric(laplace$estimate$loglik),
    as.numeric(augmented$estimate$loglik), tolerance = 0.02)
  # And the laplace route is the better fit here, which is the direction the
  # methodology predicts: it integrates the random effect rather than carrying
  # it as a state through a linearised filter. A gap inside tolerance but the
  # wrong way round would be worth knowing about.
  expect_gt(as.numeric(laplace$estimate$loglik),
    as.numeric(augmented$estimate$loglik))
})

test_that("ctGenerate draws counts rather than continuous values", {
  skip_without_julia()
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("c1", "c2"), latentNames = "eta1",
    manifesttype = c(3L, 3L), LAMBDA = matrix(1, 2, 1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6),
    MANIFESTVAR = diag(0, 2), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(1.0, 2, 1), Tpoints = 6)))
  set.seed(4)
  d <- data.frame(ctGenerate(gen, n.subjects = 15, Tpoints = 6,
    backend = "julia"))
  values <- c(d$c1, d$c2)
  values <- values[!is.na(values)]
  expect_true(length(values) > 0)
  # Whole numbers, none negative: the R generator has no measurement link and
  # would return continuous values here, which is the failure this guards.
  expect_true(all(values >= 0))
  expect_true(all(abs(values - round(values)) < 1e-8))
  expect_gt(length(unique(values)), 1)
})
