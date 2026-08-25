# intoverpop='laplace': subject-level parameter random effects integrated out
# per subject rather than by augmenting the latent state.
#
# The engine's own suite (inst/julia/.../test/test_laplace.jl) proves the
# approximation itself -- exactness for a Gaussian integrand against a
# closed-form reference, and that the outer gradient is the gradient of the
# value returned. This file proves the R side: that the metadata describing the
# random effects to the engine is right, that the two routes agree where theory
# says they must, and that the ways of asking for something unsupported all
# fail rather than quietly doing something else.

# A linear model whose only random effect is an identity-transformed
# MANIFESTMEANS. The transform matters: a random effect entering the state mean
# linearly makes the integrand exactly Gaussian in z, so Laplace is exact and
# the two routes are describing the same model rather than two similar ones.
#
# T0VAR is *fixed*, and that is not incidental. A free T0VAR and a random
# MANIFESTMEANS both describe between-subject spread, and they separate only
# through how fast the initial state decays -- so whenever the estimated drift
# comes out weak the two trade off along a ridge, the population scale wanders
# up it, and the tests below end up measuring L-BFGS's path rather than
# anything about the Laplace approximation. Fixing T0VAR at the value the data
# are generated with removes the confound and leaves the population scale
# cleanly identified, which is what these tests are actually about.
.laplace_test_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.laplace_test_data <- function(nsubjects = 30, nobs = 6) {
  set.seed(20260825)
  drift <- -0.4; diffusion <- 0.6
  rows <- lapply(seq_len(nsubjects), function(i) {
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
  })
  do.call(rbind, rows)
}


# Fitting is the expensive part of this file -- the exact Laplace gradient costs
# a forward sweep over the parameters on top of each subject's reverse pass --
# and four tests below want the same two fits. They are computed once, on first
# use, and shared. Nothing here mutates a fit, so sharing is safe.
.laplace_fit_cache <- new.env(parent = emptyenv())
.laplace_cached <- function(key, expr) {
  if (!exists(key, envir = .laplace_fit_cache, inherits = FALSE)) {
    assign(key, force(expr), envir = .laplace_fit_cache)
  }
  get(key, envir = .laplace_fit_cache, inherits = FALSE)
}
.laplace_exact_fit <- function() .laplace_cached("exact", suppressMessages(
  ctFit(.laplace_test_data(), .laplace_test_model(), backend = "julia",
    intoverpop = "laplace", optimcontrol = list(finishsamples = 100))))
.laplace_augmented_fit <- function() .laplace_cached("augmented", suppressMessages(
  ctFit(.laplace_test_data(), .laplace_test_model(), backend = "julia",
    intoverpop = TRUE, optimcontrol = list(estonly = TRUE))))

test_that("the Laplace route describes random effects without augmenting the state", {
  model <- .laplace_test_model()
  dat <- .laplace_test_data()

  laplace <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  augmented <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = TRUE, fit = FALSE))

  # The point of the route: the filtered system stays at its single-subject
  # size, where augmenting grows it by one state per random effect.
  expect_equal(laplace$nlatent_augmented, laplace$nlatent)
  expect_equal(augmented$nlatent_augmented, augmented$nlatent + 1L)

  spec <- laplace$laplace
  expect_equal(spec$nrandom, 1L)
  expect_equal(spec$param, "mmean")
  # One scale parameter, no correlations, appended after the model's own.
  expect_equal(spec$sd_index, spec$base_npar + 1L)
  expect_length(spec$cor_index, 0L)
  expect_equal(spec$npar, spec$base_npar + 1L)
  # The varying parameter is an ordinary model parameter, not a new one.
  expect_true(spec$re_index <= spec$base_npar)
  expect_true(spec$re_index %in% laplace$parameter_table$parnumber)
})

test_that("the population parameter block is laid out the way Stan lays it out", {
  # `.ctBackendPriorSpec` maps ctsem's priors onto raw positions by assuming
  # Stan's ordering -- means, then scales, then lower-triangular correlations,
  # then TI-predictor effects. If the Laplace block ever drifts from that, every
  # prior lands on the wrong parameter, so the assumption is asserted here
  # rather than left implicit.
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("Y1", "Y2"), latentNames = c("e1", "e2"),
    LAMBDA = diag(2))))
  model$pars$indvarying <- FALSE
  free <- unique(model$pars$param[!is.na(model$pars$param) &
    is.na(model$pars$value) & !grepl("[", model$pars$param, fixed = TRUE)])
  model$pars$indvarying[model$pars$param %in% free[1:3]] <- TRUE

  dat <- do.call(rbind, lapply(1:6, function(i)
    data.frame(id = i, time = 0:4, Y1 = stats::rnorm(5), Y2 = stats::rnorm(5))))
  spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))$laplace

  expect_equal(spec$nrandom, 3L)
  expect_equal(spec$sd_index, spec$base_npar + 1:3)
  expect_equal(spec$cor_index, spec$base_npar + 3L + 1:3)   # 3 choose 2
  expect_equal(spec$npar, spec$base_npar + 6L)
})

test_that("Laplace and augmented agree where the integrand is exactly Gaussian", {
  skip_without_julia()
  model <- .laplace_test_model()
  dat <- .laplace_test_data()

  laplace <- .laplace_exact_fit()
  augmented <- .laplace_augmented_fit()

  expect_true(laplace$estimate$converged)
  expect_true(augmented$estimate$converged)
  # A random effect on an identity-transformed MANIFESTMEANS is the same model
  # either way, so this is an equality, not a comparison. The tolerance is set
  # by the two routes' different numerical offsets in building the population
  # covariance (1e-8 versus none), not by the approximation.
  expect_equal(laplace$estimate$loglik, augmented$estimate$loglik, tolerance = 1e-5)

  # The estimates themselves are not compared coordinate by coordinate. The
  # augmented model carries an extra state, so its raw vector is a different
  # length in a different order, and matching by sorted value is a check on the
  # optimiser's path rather than on the model: on a weakly identified ridge two
  # routes can agree on the likelihood to eight digits while sitting at
  # different points on it. The log likelihood is the claim being made.

  expect_true(laplace$laplace$inner_converged)
  expect_false(any(laplace$laplace$hessian_repaired))
})

test_that("the random-effect population sd is reported and recovers its value", {
  skip_without_julia()
  fit <- .laplace_exact_fit()
  expect_true(fit$laplace$inner_converged)
  summarised <- summary(fit)

  expect_true("popsd" %in% names(summarised))
  expect_true("mmean" %in% rownames(summarised$popsd))
  # Simulated with a between-subject sd of 0.9. This is a recovery check, not a
  # precision claim: it fails if the population scale is being read off the
  # wrong parameter or on the wrong scale, which is what it is for.
  expect_gt(summarised$popsd["mmean", "mean"], 0.5)
  expect_lt(summarised$popsd["mmean", "mean"], 1.5)
  expect_lt(summarised$popsd["mmean", "97.5%"], 3.0)
})

test_that("unsupported ways of asking for Laplace fail rather than doing something else", {
  model <- .laplace_test_model()
  dat <- .laplace_test_data(nsubjects = 4, nobs = 4)

  expect_error(suppressMessages(ctFit(dat, model, backend = "stan",
    intoverpop = "laplace", fit = FALSE)), "requires backend='julia'")

  nonvarying <- model
  nonvarying$pars$indvarying <- FALSE
  expect_error(suppressMessages(ctFit(dat, nonvarying, backend = "julia",
    intoverpop = "laplace", fit = FALSE)), "nothing to integrate over")

})

test_that("existing intoverpop values keep their existing meanings", {
  model <- .laplace_test_model()
  dat <- .laplace_test_data(nsubjects = 4, nobs = 4)

  # TRUE and 'augmented' name the same thing, and 'auto' still resolves to it
  # when optimizing a model with declared random effects.
  logical_spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = TRUE, fit = FALSE))
  named_spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "augmented", fit = FALSE))
  auto_spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "auto", fit = FALSE))

  expect_equal(logical_spec$nlatent_augmented, named_spec$nlatent_augmented)
  expect_equal(logical_spec$nlatent_augmented, auto_spec$nlatent_augmented)
  expect_null(logical_spec$laplace)
  expect_null(named_spec$laplace)
  expect_null(auto_spec$laplace)
  expect_equal(logical_spec$intoverpop, "augmented")
})
