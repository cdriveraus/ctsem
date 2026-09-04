# priors=TRUE and per-subject scores for backend='julia'.
#
# Priors are checked against Stan rather than against a hand-computed density,
# because the thing that can actually go wrong is not the normal log-density --
# it is the *mapping* from ctsem's prior semantics onto the raw parameter
# vector. The engines are handed a list of (index, scale) pairs built by
# `.ctBackendPriorSpec()` from `standata`, so an off-by-one in the population-SD
# or TI-effect block would produce a perfectly well-formed but wrong posterior.
# Comparing the whole log probability against Stan's is what catches that, and
# the second model below deliberately has both blocks plus correlations.
#
# The score matrix is checked by the identity that defines it: the rows must sum
# to the full gradient. The summed gradient takes two per-subject shortcuts (a
# shared parameter layer, a batched matrix-exponential Frechet term) that the
# score path has to switch off, so a leak across a subject boundary shows up
# here and nowhere else.

.prior_stan_cache <- new.env(parent = emptyenv())
.prior_stan_fit <- function(stan_spec) {
  if (identical(stan_spec$standata$recompile, 0L)) {
    return(ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, stan_spec$standata))
  }
  key <- digest::digest(stan_spec$stanmodeltext)
  if (!exists(key, envir = .prior_stan_cache, inherits = FALSE)) {
    assign(key, rstan::stan_model(model_code = stan_spec$stanmodeltext),
      envir = .prior_stan_cache)
  }
  ctsem:::stan_reinitsf(get(key, envir = .prior_stan_cache, inherits = FALSE),
    stan_spec$standata)
}

.prior_simple_model <- function() {
  ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diff", 1, 1), MANIFESTVAR = matrix("mvar", 1, 1),
    MANIFESTMEANS = matrix("mm||FALSE", 1, 1), T0VAR = matrix("t0v", 1, 1),
    T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1))
}

.prior_simple_data <- function() {
  set.seed(2)
  data.frame(id = rep(1:6, each = 4), time = rep(c(0, .5, 1.2, 2), 6),
    Y1 = stats::rnorm(24, 0, .5))
}

# Random effects (population SDs and a correlation) plus a TI predictor effect,
# so every block of the prior mapping is exercised, not just rawpopmeans.
.prior_full_model <- function() {
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 2, LAMBDA = diag(2),
    MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2),
    T0MEANS = c("t0a||TRUE", "t0b||TRUE"), CINT = c("B1||TRUE", "B2||TRUE"),
    DRIFT = matrix(c("auto1", "cross21||TRUE", "cross21||TRUE", "auto2"), 2, 2,
      byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)),
    n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE))
  model$pars$group_effect[model$pars$param == "B1"] <- TRUE
  model
}

.prior_full_data <- function() {
  set.seed(5)
  do.call(rbind, lapply(1:5, function(i) data.frame(id = i, time = c(0, .5, 1.5),
    Y1 = stats::rnorm(3, 0, .5), Y2 = stats::rnorm(3, 0, .5),
    group = rep(stats::rnorm(1), 3))))
}

test_that("Stan and Julia agree with priors=TRUE, without random effects", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("digest")
  model <- .prior_simple_model()
  data <- .prior_simple_data()
  raw <- c(-.3, .2, -.5, .1, .4)

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = TRUE))
  stan_value <- rstan::log_prob(.prior_stan_fit(stan_spec), upars = raw,
    adjust_transform = FALSE, gradient = TRUE)
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = TRUE))
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)

  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)

  # The prior has to actually change the answer, or the test above would pass
  # just as well with priors ignored entirely.
  nopriors <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))
  expect_false(isTRUE(all.equal(ctJuliaEvaluate(nopriors, raw, gradient = FALSE)$value,
    julia_value$value)))
})

test_that("Stan and Julia agree with priors=TRUE, with random effects and a TI predictor", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("digest")
  model <- .prior_full_model()
  data <- .prior_full_data()

  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = TRUE))
  npar <- max(c(julia_spec$parameter_table$parnumber, julia_spec$ti_effects$coefficient),
    na.rm = TRUE)
  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)

  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = TRUE))
  stan_fit <- .prior_stan_fit(stan_spec)
  expect_equal(rstan::get_num_upars(stan_fit), npar)
  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE,
    gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)

  # Same tolerance as this model shape achieves without priors: the residual is
  # the augmented/nonlinear filter difference, not the prior term.
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-7)
  expect_equal(as.numeric(julia_value$gradient), as.numeric(attributes(stan_value)$gradient),
    tolerance = 1e-7)
})

# A plain `ctModel()` with a TI predictor, which nothing above covers: the two
# models there fix T0VAR (`T0VAR = diag(2)`) and one sets `tipredDefault =
# FALSE`, and it takes free T0VAR *plus* default TI-predictor effects to reach
# this. T0MEANS is individually varying by default, so `T0VARredundancies()`
# fixes the now-redundant free T0VAR cells -- clearing `param`, `transform` and
# `indvarying`, but not the `<TIpred>_effect` columns -- and the julia
# augmentation then reuses those same cells for the population SD and
# correlation parameters. The stale flag gave each of those a TI-predictor
# coefficient that the generated Stan model has no counterpart for, with a
# different symptom per setting: `priors = TRUE` refused the fit because its
# Stan-derived layout was three parameters short, while `priors = FALSE`
# silently estimated the three extras and reported them in `summary()$tipreds`
# as `tip_TI1_julia_popcov_1_1` and friends.
#
# MANIFESTMEANS is fixed here on purpose. Leaving it free makes it individually
# varying too, which adds carrier states, and `ctStanModelIntOverPop()` then
# rebuilds T0VAR from scratch with the effect columns already FALSE -- so the
# bug disappears. The narrow case is the one where every random effect is a
# T0MEANS row and the user's own T0VAR rows survive into the augmentation.
.prior_tipred_default_model <- function() {
  suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), MANIFESTVAR = diag(.2, 2), MANIFESTMEANS = matrix(0, 2, 1),
    n.TIpred = 1, TIpredNames = "TI1"))
}

.prior_tipred_default_data <- function() {
  set.seed(3)
  do.call(rbind, lapply(1:10, function(i) data.frame(id = i, time = 0:4,
    Y1 = stats::rnorm(5), Y2 = stats::rnorm(5), TI1 = rep(stats::rnorm(1), 5))))
}

test_that("redundant free T0VAR gives the population parameters no TI-predictor effects", {
  skip_on_cran()
  skip_without_julia()
  model <- .prior_tipred_default_model()
  data <- .prior_tipred_default_data()

  # `priors = TRUE` is what refused outright; that this returns at all is half
  # the assertion.
  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = TRUE))
  npar <- max(c(julia_spec$parameter_table$parnumber, julia_spec$ti_effects$coefficient),
    na.rm = TRUE)

  # Not vacuous: this model shape does have population parameters and does have
  # TI-predictor effects. They must simply not overlap -- Stan's population
  # block is `rawpopsdbase`/`sqrtpcov` and takes no TI effects at all.
  expect_gt(nrow(julia_spec$random_effects), 0L)
  expect_gt(nrow(julia_spec$ti_effects), 0L)
  expect_length(intersect(julia_spec$ti_effects$parameter,
    julia_spec$random_effects$parameter), 0L)

  skip_if_not_installed("rstan")
  skip_if_not_installed("digest")
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = TRUE))
  stan_fit <- .prior_stan_fit(stan_spec)
  # Stated against Stan rather than against a literal count, because the
  # property is that the two backends fit the same model. This is what catches
  # the `priors = FALSE` case, where nothing errors and the only evidence is
  # three parameters that Stan does not have.
  expect_equal(npar, rstan::get_num_upars(stan_fit))

  set.seed(8)
  raw <- stats::rnorm(npar, 0, .3)
  stan_value <- rstan::log_prob(stan_fit, upars = raw, adjust_transform = FALSE,
    gradient = TRUE)
  julia_value <- ctJuliaEvaluate(julia_spec, raw, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-7)
  expect_equal(as.numeric(julia_value$gradient),
    as.numeric(attributes(stan_value)$gradient), tolerance = 1e-6)
})

test_that("Laplace priors are refused rather than silently treated as normal", {
  # Tested against `.ctBackendPriorSpec()` directly rather than through a
  # fitted model: setting `model$laplaceprior` stores the field but does not
  # reach `standata$laplaceprior`, which stays all-zero, so there is no way from
  # here to drive the guard end to end. Asserting on the guard itself is honest
  # about what is being checked; a `ctFit()` call that silently produced normal
  # priors would pass an end-to-end test just as well.
  standata <- list(nparams = 3L, nindvarying = 0L, nindvaryingoffdiagonals = 0L,
    ntipredeffects = 0L, priormod = 1, nsubsets = 1,
    laplaceprior = c(0L, 1L, 0L), laplacetipreds = 0L, laplaceprioronly = 0L)
  expect_error(ctsem:::.ctBackendPriorSpec(standata, 3L),
    "Laplace priors are not implemented")

  standata$laplaceprior <- c(0L, 0L, 0L)
  standata$laplacetipreds <- 1L
  expect_error(ctsem:::.ctBackendPriorSpec(standata, 3L),
    "Laplace priors on TI predictor effects")

  standata$laplacetipreds <- 0L
  expect_equal(ctsem:::.ctBackendPriorSpec(standata, 3L)$index, 1:3)

  # A layout that does not account for every free parameter must refuse rather
  # than silently apply priors to the wrong ones.
  expect_error(ctsem:::.ctBackendPriorSpec(standata, 5L),
    "Cannot map ctsem's priors")
})

test_that("per-subject scores sum to the gradient, with and without priors", {
  model <- .prior_full_model()
  data <- .prior_full_data()

  for (priors in c(FALSE, TRUE)) {
    spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
      priors = priors))
    npar <- max(c(spec$parameter_table$parnumber, spec$ti_effects$coefficient),
      na.rm = TRUE)
    set.seed(8)
    raw <- stats::rnorm(npar, 0, .3)
    fit <- structure(list(model_spec = spec, backend = "julia"),
      class = c("ctJuliaFit", "ctFit"))

    scores <- ctsem:::.ctBackendScoreMatrix(fit, raw)
    gradient <- ctJuliaEvaluate(spec, raw, gradient = TRUE)$gradient
    expect_equal(dim(scores), c(length(spec$subject_starts), npar))
    expect_equal(colSums(scores), gradient, tolerance = 1e-9)
    # Every subject contributes something; an all-zero row would mean a subject
    # was filtered but never unwound.
    expect_true(all(apply(scores, 1, function(row) any(row != 0))))
  }
})

test_that("score-based uncertainty methods work for backend fits", {
  skip_on_cran()
  skip_without_julia()
  model <- .prior_full_model()
  # More subjects than parameters, so the score covariance is not rank limited.
  set.seed(9)
  data <- do.call(rbind, lapply(1:40, function(i) data.frame(id = i,
    time = c(0, .5, 1.2, 2), Y1 = stats::rnorm(4, 0, .5), Y2 = stats::rnorm(4, 0, .5),
    group = rep(stats::rnorm(1), 4))))

  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  for (method in c("opg", "sandwich", "bootstrap")) {
    updated <- suppressWarnings(suppressMessages(
      ctOptimUncertainty(fit, uncertainty = method, finishsamples = 50, verbose = 0)))
    covariance <- updated$estimate$cov
    expect_equal(dim(covariance), c(length(fit$estimate$raw), length(fit$estimate$raw)))
    expect_true(all(is.finite(covariance)))
    expect_true(all(diag(covariance) > 0))
    expect_identical(updated$uncertainty$settings$method, method)
  }
  # fullbootstrap still needs refits per resample, and says so.
  expect_error(ctOptimUncertainty(fit, uncertainty = "fullbootstrap"),
    "not available for backend")
})
