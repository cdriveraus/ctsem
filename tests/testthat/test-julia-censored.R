# Censored manifest variables: `manifesttype = 4`, Gaussian within known limits.
#
# Censored is the one non-Gaussian type with a measurement error of its own, so
# it is the one whose standard deviation stays free and has to be differentiated.
# That is what most of this file is about: the engine's suite proves the kernel
# against closed forms, and what is left for here is that the standard deviation
# reaches the engine, comes back correctly in the gradient, and that everything
# the R side has to know about limits is checked rather than assumed.

.censored_data <- function(nsubjects = 40, nobs = 8, lower = 0, upper = 5,
  seed = 7) {
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
  latent <- 2.5 + d$eta + stats::rnorm(nrow(d), 0, 0.5)
  d$y <- pmin(pmax(latent, lower), upper)
  d$eta <- NULL
  d
}

.censored_model <- function(lower = 0, upper = 5, sd = NULL) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 4L, censormin = lower, censormax = upper,
    LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm"), CINT = matrix(0),
    T0MEANS = matrix(0),
    MANIFESTVAR = matrix(if (is.null(sd)) "mvar" else sd))))
  m$pars$indvarying <- FALSE
  m
}

test_that("ctModel takes manifesttype 4 and its limits", {
  m <- .censored_model()
  expect_equal(unname(m$manifesttype), 4L)
  expect_equal(unname(m$censormin), 0)
  expect_equal(unname(m$censormax), 5)

  build <- function(...) suppressWarnings(suppressMessages(ctModel(
    type = "ct", n.latent = 1, n.manifest = 1, manifestNames = "y",
    latentNames = "eta1", LAMBDA = matrix(1), ...)))
  # Censored nowhere is a Gaussian variable, and saying so is more useful than
  # fitting an identical model under a different name.
  expect_error(build(manifesttype = 4L), "at least one")
  expect_error(build(manifesttype = 4L, censormin = 5, censormax = 0),
    "must be below")
  # A limit left on a variable that is no longer censored must not apply.
  plain <- build(manifesttype = 0L, censormin = 1, censormax = 2)
  expect_equal(unname(plain$censormin), -Inf)
  expect_equal(unname(plain$censormax), Inf)
})

test_that("a censored variable keeps its measurement standard deviation", {
  skip_without_julia()
  # Every other non-Gaussian type has its randomness supplied by the link, so
  # ctFit fixes MANIFESTVAR for them. Censored is Gaussian within its limits and
  # would have no scale at all if that happened here.
  fit <- suppressWarnings(suppressMessages(ctFit(.censored_data(),
    .censored_model(), backend = "julia", intoverpop = "augmented",
    fit = FALSE)))
  free <- fit$parameter_table[fit$parameter_table$matrix %in% "MANIFESTVAR", ]
  expect_true(any(!is.na(free$parnumber)))
})

test_that("print names censored variables with their limits", {
  described <- paste(utils::capture.output(print(.censored_model())),
    collapse = " ")
  expect_match(described, "y \\(censored, min 0, max 5\\)")
  onesided <- paste(utils::capture.output(print(
    .censored_model(lower = 0, upper = Inf))), collapse = " ")
  expect_match(onesided, "no max")
})

test_that("censored data is checked against the declared limits", {
  d <- .censored_data(nsubjects = 6, nobs = 4)
  m <- .censored_model()

  beyond <- d
  beyond$y[3] <- -1
  expect_error(suppressMessages(ctFit(beyond, m, backend = "julia",
    fit = FALSE)), "below its censormin")

  above <- d
  above$y[3] <- 6
  expect_error(suppressMessages(ctFit(above, m, backend = "julia",
    fit = FALSE)), "above its censormax")

  # Legal but worth saying: censoring that never applies, and censoring that
  # applies to everything, are both models the user probably did not mean.
  never <- d
  never$y <- pmin(pmax(never$y, 0.5), 4.5)
  expect_warning(suppressMessages(ctFit(never, m, backend = "julia",
    fit = FALSE)), "no observations at either limit")
})

test_that("stan refuses censored rather than ignoring the censoring", {
  expect_error(suppressMessages(ctFit(.censored_data(nsubjects = 5, nobs = 4),
    .censored_model(), backend = "stan",
    optimcontrol = list(estonly = TRUE))), "need backend=\"julia\"")
})

test_that("the adjoint matches forward mode with censoring present", {
  skip_without_julia()
  # The standard deviation's cotangent is the new path here: it is read from
  # MANIFESTVAR's own diagonal and handed straight back there, rather than going
  # through the covariance construction the Gaussian rows use.
  handle <- suppressWarnings(suppressMessages(ctFit(.censored_data(),
    .censored_model(), backend = "julia", intoverpop = "augmented",
    fit = FALSE)))
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(2)
  for (trial in 1:3) {
    at <- stats::rnorm(npar, 0, 0.25)
    adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "adjoint")$gradient)
    forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
      gradient_method = "forward")$gradient)
    expect_equal(adjoint, forward, tolerance = 1e-8)
  }
})

test_that("the gradient is right when most observations are at a limit", {
  skip_without_julia()
  # The regression this pins. ForwardDiff breaks comparison ties on the
  # partials, so an observation exactly at a limit compared against a limit
  # carrying a derivative seed took the interior branch during differentiation
  # while the forward pass took the censored one. The error grew with the number
  # of censored observations: 1% at two of them, 263% at two hundred and
  # eighty-nine.
  d <- .censored_data(lower = 2.4, upper = 2.6)
  m <- .censored_model(lower = 2.4, upper = 2.6)
  atlimit <- sum(d$y <= 2.4 + 1e-8) + sum(d$y >= 2.6 - 1e-8)
  expect_gt(atlimit, nrow(d) / 2)

  handle <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "augmented", fit = FALSE)))
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(2)
  at <- stats::rnorm(npar, 0, 0.25)
  adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE)$gradient)
  eps <- 1e-5
  differenced <- vapply(seq_len(npar), function(i) {
    up <- at; up[i] <- up[i] + eps
    down <- at; down[i] <- down[i] - eps
    (as.numeric(ctJuliaEvaluate(handle, up)$value) -
        as.numeric(ctJuliaEvaluate(handle, down)$value)) / (2 * eps)
  }, numeric(1))
  expect_equal(adjoint, differenced, tolerance = 1e-5)
})

test_that("a censored model recovers what generated it", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctFit(.censored_data(nsubjects = 50,
    nobs = 10), .censored_model(), backend = "julia",
    intoverpop = "augmented", optimcontrol = list(estonly = TRUE))))
  expect_true(isTRUE(fit$estimate$converged))
  means <- summary(fit)$popmeans
  expect_equal(unname(means["drift_eta1", "mean"]), -0.4, tolerance = 0.35)
  expect_equal(unname(means["mm", "mean"]), 2.5, tolerance = 0.3)
  expect_equal(unname(means["mvar", "mean"]), 0.5, tolerance = 0.25)
})

test_that("ctGenerate respects the censoring limits", {
  skip_without_julia()
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y", latentNames = "eta1",
    manifesttype = 4L, censormin = 0, censormax = 3, LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.5),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(1.5), Tpoints = 6)))
  set.seed(4)
  d <- data.frame(ctGenerate(gen, n.subjects = 20, Tpoints = 6,
    backend = "julia"))
  values <- d$y[!is.na(d$y)]
  expect_true(length(values) > 0)
  expect_true(all(values >= 0 - 1e-8))
  expect_true(all(values <= 3 + 1e-8))
  # Some should land on the limits and some inside, or the limits are doing
  # nothing and the test would pass on a generator that ignored them.
  expect_gt(sum(values <= 1e-8 | values >= 3 - 1e-8), 0)
  expect_gt(sum(values > 1e-8 & values < 3 - 1e-8), 0)
})
