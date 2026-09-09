# Binary indicators, end to end.
#
# `manifesttype = 1` changes the measurement model from Gaussian to a Bernoulli
# link, which touches the filter, the measurement Jacobian, and what a residual
# means. The two older files in this directory exercise it on larger models;
# this one is small enough to run routinely and checks the things most likely
# to rot: that the path still recovers what generated the data, that the julia
# backend refuses rather than quietly doing something else, and that the
# functions people reach for after a fit still work on one.

.binary_data <- function(nsubjects = 60, nobs = 12, nindicators = 3) {
  invlog <- function(x) exp(x) / (1 + exp(x))
  set.seed(11)
  gm <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8),
    MANIFESTVAR = matrix(0.001), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(0), Tpoints = nobs))
  latent <- ctGenerate(gm, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r")
  d <- data.frame(latent)
  for (i in seq_len(nindicators)) {
    d[[paste0("b", i)]] <- stats::rbinom(nrow(d), 1, invlog(d$eta))
  }
  d$eta <- NULL
  d
}

.binary_model <- function(nindicators = 3, manifestvar = 0) {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = nindicators, manifestNames = paste0("b", 1:nindicators),
    latentNames = "eta1", LAMBDA = matrix(1, nindicators, 1),
    MANIFESTMEANS = matrix(0, nindicators, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(manifestvar, nindicators)))
  m$manifesttype[] <- 1L
  m$pars$indvarying <- FALSE
  m
}

test_that("a binary model recovers what generated it", {
  skip_on_cran()
  fit <- suppressWarnings(suppressMessages(
    ctFit(.binary_data(), .binary_model(), cores = 1, verbose = 0)))
  est <- summary(fit)$popmeans
  # Intervals rather than point estimates: three binary indicators carry much
  # less information than one continuous one, so the point estimate is allowed
  # to be off as long as the interval covers.
  expect_lt(est["drift_eta1", "2.5%"], -0.3)
  expect_gt(est["drift_eta1", "97.5%"], -0.3)
  expect_lt(est["diff_eta1", "2.5%"], 0.8)
  expect_gt(est["diff_eta1", "97.5%"], 0.8)
})

test_that("the julia backend fits binary indicators rather than refusing them", {
  skip_without_julia()
  # This asserted a refusal, from back when the julia filter had no measurement
  # link and treating a binary indicator as Gaussian would have produced a fit
  # that looked fine and answered a different question. The filter integrates
  # the observation now, so the contract is the opposite one: it fits, and what
  # is refused is a manifest type beyond ordinal (test-julia-backend.R).
  #
  # Small and degenerate on purpose -- five subjects and four occasions is
  # where the marshalling of integer manifest columns used to fail.
  fit <- suppressWarnings(suppressMessages(
    ctFit(.binary_data(nsubjects = 5, nobs = 4), .binary_model(),
      backend = "julia", cores = 1)))
  expect_s3_class(fit, "ctJuliaFit")
  expect_true(is.finite(fit$estimate$loglik))
})

test_that("a fixed non-zero measurement variance on a binary indicator warns", {
  skip_on_cran()
  # Free variances are fixed silently because there is a right answer. A value
  # the user stated is theirs, so it is questioned rather than overwritten: the
  # Bernoulli link already supplies the randomness.
  expect_warning(
    suppressMessages(ctFit(.binary_data(nsubjects = 5, nobs = 4),
      .binary_model(manifestvar = 0.5), cores = 1, fit = FALSE)),
    "measurement link")
})

test_that("the deprecated binomial argument no longer disables the filter", {
  skip_on_cran()
  # It used to set `intoverstates = FALSE`, which the very next check in ctFit
  # warns is unusable for optimization -- the documented shortcut put a user
  # into the state the code calls unreliable.
  d <- .binary_data(nsubjects = 5, nobs = 4)
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 3,
    manifestNames = paste0("b", 1:3), latentNames = "eta1",
    LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
    CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 3)))
  m$pars$indvarying <- FALSE
  warnings <- character()
  prepared <- withCallingHandlers(
    suppressMessages(ctFit(d, m, binomial = TRUE, cores = 1, fit = FALSE)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("binomial argument is deprecated", warnings)))
  expect_false(any(grepl("intoverstates=TRUE required", warnings)))
  # It still does the thing it is for, checked on what reaches the engine
  # rather than on the model object -- `ctstanmodelbase` is the user's original
  # and deliberately keeps its zeros.
  expect_true(all(prepared$standata$manifesttype == 1L))
  expect_true(all(prepared$ctstanmodel$manifesttype == 1L))
})

test_that("post-fit functions work on a binary fit", {
  skip_on_cran()
  fit <- suppressWarnings(suppressMessages(
    ctFit(.binary_data(), .binary_model(), cores = 1, verbose = 0)))
  works <- function(expr) {
    result <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
    !inherits(result, "try-error")
  }
  expect_true(works(summary(fit)))
  expect_true(works(ctSummaryMatrices(fit)))
  expect_true(works(ctKalman(fit, subjects = 1)))
  expect_true(works(ctPredict(fit, subjects = 1)))
  expect_true(works(ctExtract(fit)))
  expect_true(works(ctPostPredPlots(fit)))
})

test_that("ctTIpredEffects explains itself when there are no predictors", {
  skip_on_cran()
  # Not binary-specific -- it did this on any model without predictors -- but
  # found here, and "subscript out of bounds" names an internal rather than the
  # problem.
  fit <- suppressWarnings(suppressMessages(
    ctFit(.binary_data(), .binary_model(), cores = 1, verbose = 0)))
  expect_error(suppressMessages(ctTIpredEffects(fit)),
    "no time independent predictors")
})
