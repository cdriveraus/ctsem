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

test_that("cores splits the Laplace subject loop without changing the answer", {
  skip_without_julia()
  model <- .laplace_test_model()
  dat <- .laplace_test_data()

  # `cores` reaches the Laplace path the same way it reaches the ordinary one,
  # as a chunk count for the subject loop. Splitting only reorders a sum, so
  # the two runs differ by floating-point association and nothing else -- if
  # they differed by more, chunks would be sharing scratch they should own.
  serial <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", cores = 1, optimcontrol = list(estonly = TRUE)))
  split <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", cores = 2, optimcontrol = list(estonly = TRUE)))

  expect_equal(split$estimate$loglik, serial$estimate$loglik, tolerance = 1e-8)
  expect_equal(split$estimate$raw, serial$estimate$raw, tolerance = 1e-6)
  expect_true(split$laplace$inner_converged)
})

test_that("verbose reports the optimiser trace and the inner solve", {
  skip_without_julia()
  model <- .laplace_test_model()
  dat <- .laplace_test_data(nsubjects = 8, nobs = 5)

  # Neither half of the trace travels on stdout, and they do not travel
  # together, so this captures both channels separately and asserts on the one
  # each line is actually on. The engine's `verbose` lines are raw writes to
  # Julia's `stderr`, which JuliaConnectoR relays to R's `stderr()`, so a sink
  # catches them and `suppressMessages()` never could. The progress line goes
  # through `.ctProgressSink()`, which emits it with `message()` -- a condition,
  # which testthat muffles as it records it, so under testthat it does not reach
  # `stderr()` for a sink to catch either. Capturing stdout and suppressing
  # messages, as this did, saw neither: three assertions passed vacuously for as
  # long as the greps happened to be wrong.
  relayed <- character(); emitted <- character()
  invisible(capture.output(
    relayed <- capture.output(
      emitted <- capture_messages(ctFit(dat, model, backend = "julia",
        intoverpop = "laplace", verbose = 1,
        optimcontrol = list(estonly = TRUE))),
      type = "message")))

  # The inner solve is part of the objective, so its status belongs in the
  # trace rather than only on the fit object. Asserted nowhere else: the engine
  # suite proves the solve, and this is the only check that its status is
  # reported.
  expect_true(any(grepl("Laplace: inner modes", relayed, fixed = TRUE)))
  # The progress reporter replaced Optim's own trace, so the assertion is on
  # what it emits: a labelled line carrying the objective. `Iter` was the old
  # format's column heading and says nothing about whether the fit is going
  # anywhere.
  #
  # The line's content is `test_progress.jl`'s subject and its delivery as a
  # message is `test-julia-session.R`'s, so what is left for a real fit to prove
  # is that the reporter is *wired in* -- which is a failure neither of those can
  # see, and one that has happened: `progress = isTRUE(progress)` named an
  # element of the list being built rather than anything in scope, so the sink
  # was never called.
  expect_true(any(grepl("optimise", emitted, fixed = TRUE)))
  expect_true(any(grepl("logpost", emitted, fixed = TRUE)))
})

# --- subjects nested in studies ----------------------------------------------

.laplace_nested_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"),
    id = c("subject", "study"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  # A grouping level above the subject declares its varying parameters in
  # `indvarying_<idname>`, the same way the subject level uses `indvarying`.
  model$pars$indvarying_study <- FALSE
  model$pars$indvarying_study[model$pars$param %in% "mmean"] <- TRUE
  model
}

.laplace_nested_data <- function(nstudy = 8, npersub = 5, nobs = 5) {
  set.seed(20260825)
  rows <- list(); sid <- 0
  for (g in seq_len(nstudy)) {
    studyeffect <- stats::rnorm(1, 0, 0.10)
    for (j in seq_len(npersub)) {
      sid <- sid + 1
      intercept <- 10 * (0.15 + studyeffect + stats::rnorm(1, 0, 0.06))
      state <- stats::rnorm(1, 0, 0.5); out <- numeric(nobs)
      for (t in seq_len(nobs)) {
        if (t > 1) {
          decay <- exp(-0.4)
          state <- decay * state + stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - decay^2)))
        }
        out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
      }
      rows[[length(rows) + 1L]] <- data.frame(subject = sid, study = g,
        time = seq_len(nobs) - 1, Y1 = out)
    }
  }
  do.call(rbind, rows)
}

test_that("a vector of id columns builds one level per id", {
  model <- .laplace_nested_model()
  dat <- .laplace_nested_data()
  spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  laplace <- spec$laplace

  expect_equal(laplace$nlevels, 2L)
  expect_equal(vapply(laplace$levels, function(x) x$name, character(1)),
    c("subject", "study"))
  expect_equal(laplace$levels[[1]]$ngroups, 40L)   # subjects
  expect_equal(laplace$levels[[2]]$ngroups, 8L)    # studies
  expect_equal(laplace$levels[[1]]$param, "mmean")
  expect_equal(laplace$levels[[2]]$param, "mmean")

  # Every level's scale gets its own slot in the raw vector, and `npar` counts
  # them all. Sizing the raw vector from the subject level alone left the study
  # scale past its end, where it silently never moved.
  expect_equal(laplace$levels[[1]]$sd_index, laplace$base_npar + 1L)
  expect_equal(laplace$levels[[2]]$sd_index, laplace$base_npar + 2L)
  expect_equal(laplace$npar, laplace$base_npar + 2L)
  expect_length(ctsem:::.ctJuliaInitialValues(laplace$npar), laplace$npar)
})

test_that("strict nesting is required, and said so when it is not", {
  model <- .laplace_nested_model()
  dat <- .laplace_nested_data()
  dat$study[dat$subject == 1][1] <- 99   # one subject in two studies
  expect_error(suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE)), "strictly nested")

  missingcol <- .laplace_nested_data()
  missingcol$study <- NULL
  expect_error(suppressMessages(ctFit(missingcol, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE)), "not found in the data")
})

test_that("both levels of a nested fit are estimated and reported separately", {
  skip_without_julia()
  model <- .laplace_nested_model()
  dat <- .laplace_nested_data(nstudy = 12, npersub = 5)

  fit <- suppressMessages(suppressWarnings(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", optimcontrol = list(finishsamples = 60))))
  expect_true(fit$laplace$inner_converged)
  expect_equal(fit$laplace$nlevels, 2L)

  summarised <- summary(fit, residualcov = FALSE)
  # With more than one level nothing goes in the unlabelled `popsd` slot: a bare
  # table would leave the reader guessing whether it described spread between
  # subjects or between studies.
  expect_null(summarised$popsd)
  expect_true("popsd.subject" %in% names(summarised))
  expect_true("popsd.study" %in% names(summarised))
  expect_true(!is.null(summarised$randomEffects))
  expect_equal(length(summarised$randomEffects), 2L)

  # Simulated with a subject sd of 0.6 and a study sd of 1.0 on the transformed
  # scale. A recovery check, not a precision claim: it fails if a level is read
  # off the wrong parameter, on the wrong scale, or -- as it did -- pinned at
  # its starting value because the raw vector was too short to hold it.
  subject <- summarised$popsd.subject["mmean", "mean"]
  study <- summarised$popsd.study["mmean", "mean"]
  expect_gt(subject, 0.3); expect_lt(subject, 1.1)
  expect_gt(study, 0.4);   expect_lt(study, 2.2)

  # The two levels must not be reported as the *same* number, which is what a
  # mis-wired level index would produce. Requiring them to differ by some
  # margin instead would be wrong: two genuinely distinct levels can land close
  # together by chance, and this assertion is about wiring, not about spacing.
  expect_false(isTRUE(all.equal(subject, study)))
})

test_that("the Laplace prior layout reproduces the Stan one at a single level", {
  # `.ctBackendPriorSpec` encodes the generated Stan model's raw ordering and is
  # what every single-level fit has used. The Laplace construction builds the
  # same thing from its own layout so it can carry more than one level. Where
  # both apply they must agree exactly, or a fit would silently change its
  # priors depending on which route built them.
  standata <- list(nparams = 4L, nindvarying = 2L, nindvaryingoffdiagonals = 1L,
    ntipredeffects = 2L, tipredeffectscale = 0.5, priormod = 1, nsubsets = 1,
    laplaceprior = 0L, laplacetipreds = 0L, laplaceprioronly = 0L)
  laplace <- list(npar = 7L, levels = list(list(sd_index = 5:6, cor_index = 7L)))

  stanform <- ctsem:::.ctBackendPriorSpec(standata, 9L)
  laplaceform <- ctsem:::.ctBackendLaplacePriorSpec(standata, laplace, 9L)
  expect_equal(laplaceform, stanform)
})

test_that("every level gets the same prior shape, and its own sdscale", {
  standata <- list(nparams = 4L, nindvarying = 1L, nindvaryingoffdiagonals = 0L,
    ntipredeffects = 0L, priormod = 1, nsubsets = 1,
    laplaceprior = 0L, laplacetipreds = 0L, laplaceprioronly = 0L)
  # Two levels, one effect each: a scale slot per level, no correlations.
  laplace <- list(npar = 6L, levels = list(list(sd_index = 5L, cor_index = integer()),
    list(sd_index = 6L, cor_index = integer())))
  spec <- ctsem:::.ctBackendLaplacePriorSpec(standata, laplace, 6L)
  expect_equal(spec$index, 1:6)
  # The same shape at both levels -- the parameterisation is identical, so the
  # prior is. That it does more work over few studies than over many subjects
  # is a property of the design, not something the prior should try to fix.
  expect_equal(spec$scale, rep(1, 6))

  # And a level scales its own population sd rather than borrowing level one's.
  model <- .laplace_nested_model()
  model$pars$sdscale_study[model$pars$param %in% "mmean"] <- 3
  dat <- .laplace_nested_data()
  levels <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))$laplace$levels
  expect_equal(levels[[1]]$sd_scale, 1)
  expect_equal(levels[[2]]$sd_scale, 3)
})

test_that("a grouping id creates its level columns, defaulting sensibly", {
  model <- .laplace_nested_model()
  expect_true("sdscale_study" %in% names(model$pars))
  expect_true("indvarying_study" %in% names(model$pars))
  # sdscale defaults to the subject level's, which is 1 for free parameters.
  free <- !is.na(model$pars$param) & is.na(model$pars$value)
  expect_true(all(model$pars$sdscale_study[free] == 1))

  # An outer level varies only where asked. Defaulting it the way `indvarying`
  # defaults -- TRUE for T0MEANS, MANIFESTMEANS and CINT -- would make every
  # model with a grouping id enormously parameterised without the user asking.
  plain <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    id = c("subject", "study"))))
  expect_false(any(plain$pars$indvarying_study %in% TRUE))
  expect_true(any(plain$pars$indvarying %in% TRUE))
})

test_that("a TI predictor operates at a grouping level when it varies there", {
  skip_without_julia()
  # A covariate constant within study and varying between studies is a
  # study-level predictor: ctsem applies TI effects additively to each
  # subject's raw vector, so a study-constant shift *is* a study-level fixed
  # effect. Nothing special is needed for that -- but the raw vector has to be
  # long enough to hold the coefficient, which sits after every level's
  # population block, and it was not.
  set.seed(31)
  nstudy <- 10; npersub <- 4; nobs <- 5
  rows <- list(); sid <- 0
  for (g in seq_len(nstudy)) {
    Z <- stats::rnorm(1)
    studyeffect <- 1.2 * Z + stats::rnorm(1, 0, 0.04)
    for (j in seq_len(npersub)) {
      sid <- sid + 1
      intercept <- 10 * (0.15 + studyeffect + stats::rnorm(1, 0, 0.06))
      x <- stats::rnorm(1, 0, 0.5); out <- numeric(nobs)
      for (t in seq_len(nobs)) {
        if (t > 1) {
          decay <- exp(-0.4)
          x <- decay * x + stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - decay^2)))
        }
        out[t] <- x + intercept + stats::rnorm(1, 0, 0.3)
      }
      rows[[length(rows) + 1L]] <- data.frame(subject = sid, study = g,
        time = seq_len(nobs) - 1, Y1 = out, Z = Z)
    }
  }
  dat <- do.call(rbind, rows)

  build <- function(withZ) {
    model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
      manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
      T0MEANS = matrix(0), CINT = matrix(0), T0VAR = matrix(0.5),
      MANIFESTMEANS = matrix("mmean"), id = c("subject", "study"),
      TIpredNames = if (withZ) "Z" else NULL)))
    model$pars$indvarying <- FALSE
    model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
    model$pars$indvarying_study[model$pars$param %in% "mmean"] <- TRUE
    if (withZ) {
      model$pars$Z_effect <- FALSE
      model$pars$Z_effect[model$pars$param %in% "mmean"] <- TRUE
    }
    model
  }

  spec <- suppressMessages(ctFit(dat, build(TRUE), backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  # The coefficient sits past every level's population block, and the raw
  # vector must reach it.
  expect_gt(spec$ti_effects$coefficient, spec$laplace$npar)
  expect_length(ctsem:::.ctJuliaInitialValues(
    max(spec$laplace$npar, spec$ti_effects$coefficient)),
    spec$ti_effects$coefficient)

  fits <- lapply(c(FALSE, TRUE), function(withZ)
    suppressMessages(suppressWarnings(ctFit(dat, build(withZ), backend = "julia",
      intoverpop = "laplace", optimcontrol = list(estonly = TRUE)))))
  sdof <- function(fit, level) {
    lv <- fit$model_spec$laplace$levels[[level]]
    10 * (log1p(exp(2 * fit$estimate$raw[lv$sd_index] - 1)) + 1e-10) * lv$sd_scale
  }

  # The substantive check: a real study-level predictor soaks up between-study
  # variance, so the *study* random-effect sd collapses while the subject one
  # does not. Asserting only that the fit ran would not distinguish a predictor
  # that works from one silently applied at the wrong level.
  expect_lt(sdof(fits[[2]], 2), sdof(fits[[1]], 2) / 4)
  expect_equal(sdof(fits[[2]], 1), sdof(fits[[1]], 1), tolerance = 0.3)
  expect_gt(fits[[2]]$estimate$loglik, fits[[1]]$estimate$loglik)
})

# --- a grouping level with no subject level effects --------------------------

# One intercept per study with subjects exchangeable inside it is an ordinary
# model, and nothing about the method needs a level below the varying one. It
# was refused anyway: every "does this model have random effects" test read the
# subject level's `indvarying` column and no other.

.laplace_studyonly_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"),
    id = c("subject", "study"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying_study <- model$pars$param %in% "mmean"
  model
}

test_that("a study level effect needs no subject level effect", {
  skip_without_julia()
  dat <- .laplace_nested_data()
  spec <- suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))

  expect_equal(spec$laplace$nlevels, 2L)
  expect_equal(vapply(spec$laplace$levels, function(x) x$nrandom, integer(1)),
    c(0L, 1L))
  expect_equal(spec$laplace$nrandom, 1L)

  fit <- suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "julia", intoverpop = "laplace",
    optimcontrol = list(estonly = TRUE)))
  expect_true(is.finite(fit$estimate$loglik))
})

test_that("one subject per study reproduces the subject level fit exactly", {
  skip_without_julia()
  # The check that is not parity between two of our own paths: with one subject
  # in every study the two specifications are the same model, so the study
  # level route has to return what the subject level route -- tested against
  # closed forms elsewhere in this file -- already returns. Anything the empty
  # inner level perturbs shows up here as a difference in the likelihood.
  dat <- .laplace_nested_data(nstudy = 20, npersub = 1, nobs = 6)

  subjectlevel <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"), id = "subject")))
  subjectlevel$pars$indvarying <- subjectlevel$pars$param %in% "mmean"

  set.seed(99)
  one <- suppressMessages(ctFit(dat, subjectlevel, backend = "julia",
    intoverpop = "laplace", cores = 1))
  set.seed(99)
  two <- suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "julia", intoverpop = "laplace", cores = 1))

  expect_equal(two$estimate$loglik, one$estimate$loglik)
  expect_equal(two$estimate$raw, one$estimate$raw)
  expect_equal(summary(two)$popsd, summary(one)$popsd)
})

test_that("a varying grouping level chooses its own route", {
  skip_without_julia()
  dat <- .laplace_nested_data(nstudy = 4, npersub = 3, nobs = 4)

  # 'auto' resolved to the augmented route, which reads `indvarying` alone and
  # would have dropped the study effect without saying so.
  auto <- suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "julia", fit = FALSE))
  expect_identical(auto$args$resolved$intoverpop, "laplace")
  expect_equal(auto$laplace$nlevels, 2L)

  # Asked for explicitly, the same route is refused rather than silently
  # narrowed, and the message names the parameter it could not carry.
  expect_error(suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "julia", intoverpop = "augmented")), "mmean")
  # Same for the backend that has one grouping level, whatever the route.
  expect_error(suppressMessages(ctFit(dat, .laplace_studyonly_model(),
    backend = "stan", optimize = FALSE, fit = FALSE)), "mmean")
})


test_that("subject parameters are available at draws, with spread", {
  skip_without_julia()
  fit <- .laplace_exact_fit()

  point <- ctsem:::.ctBackendSubjectPars(fit, pointest = TRUE)
  draws <- ctsem:::.ctBackendSubjectPars(fit, pointest = FALSE)
  expect_equal(dim(point)[1], 1L)
  expect_gt(dim(draws)[1], 1L)
  expect_equal(dim(draws)[2:3], dim(point)[2:3])
  expect_equal(dimnames(draws)$param, dimnames(point)$param)

  # The draw distribution should sit around the point estimate...
  centres <- apply(draws[, , 1, drop = FALSE], 2, mean)
  expect_lt(max(abs(centres - point[1, , 1])), 0.5)
  # ...and actually have spread. Reporting only the point estimate was the gap:
  # a subject parameter carries uncertainty from the population parameters as
  # well as from its own conditional distribution.
  expect_true(all(apply(draws[, , 1, drop = FALSE], 2, sd) > 0))
})

test_that("the raw vector reaches every index the layout references", {
  # Three separate bugs of this shape have shipped into this branch: a level's
  # population scales, then the TI-predictor coefficients, each indexed past
  # the end of a vector sized from only part of the layout. Neither failed
  # usefully. The invariant is asserted directly so the next one cannot.
  model <- .laplace_nested_model()
  dat <- .laplace_nested_data()
  spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  laplace <- spec$laplace

  reach <- c(spec$parameter_table$parnumber, spec$ti_effects$coefficient,
    unlist(lapply(laplace$levels, function(x)
      c(x$re_index, x$sd_index, x$cor_index))))
  expect_lte(max(reach, na.rm = TRUE),
    max(laplace$npar, spec$ti_effects$coefficient, na.rm = TRUE))

  # And the check itself fires rather than sitting inert.
  broken <- laplace
  broken$levels[[2]]$sd_index <- 999L
  expect_error(ctsem:::.ctJuliaCheckLayout(spec$parameter_table, broken,
    spec$ti_effects, laplace$npar), "reach raw parameter 999")
})

test_that("ctKalman on a Laplace fit selects which levels it conditions on", {
  skip_without_julia()
  model <- .laplace_nested_model()
  dat <- .laplace_nested_data(nstudy = 5, npersub = 4)
  fit <- suppressMessages(suppressWarnings(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", optimcontrol = list(estonly = TRUE))))

  # A Laplace fit can be filtered, so its prior residuals are available -- and
  # they are conditioned on each subject's mode from their whole record, which
  # is a different quantity from the augmented route's filtered residuals and
  # is labelled as such.
  expect_false(is.null(fit$priorerrors))
  expect_match(attr(fit$priorerrors, "conditioning"), "whole record")

  study <- dat$study[match(unique(dat$subject), dat$subject)]
  first_of <- function(k) { first <- !duplicated(k$id); k$yprior[1, first, 1] }
  # The *prior* prediction at each subject's first row, because that is driven
  # by the parameters alone. The smoothed one uses the subject's own data
  # whatever the parameters, so it would look nearly identical at every level
  # and would not test anything.
  levels <- lapply(c("subject", "study", "population"), function(lv)
    first_of(suppressMessages(ctKalmanArray(fit, pointest = TRUE,
      randomEffects = lv))))
  names(levels) <- c("subject", "study", "population")

  # Naming a level includes it and everything outside it, and zeroes what is
  # inside: subjects differ; study-level trajectories are shared within a study
  # but differ between them; population is one trajectory for everyone.
  expect_gt(stats::sd(levels$subject), 0)
  expect_equal(max(tapply(levels$study, study, function(x) diff(range(x)))), 0)
  expect_gt(stats::sd(levels$study), 0)
  expect_equal(stats::sd(levels$population), 0)
  # Including the subject level must actually add something over the study one.
  expect_gt(stats::sd(levels$subject - levels$study[match(study, study)]), 0)

  expect_error(suppressMessages(ctKalmanArray(fit, pointest = TRUE,
    randomEffects = "nope")), "must be 'population' or one of the model's id")
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

# ---------------------------------------------------------------------------
# The user-facing surface
# ---------------------------------------------------------------------------
#
# Everything above this point asserts on numbers the engine returns. That left
# a gap: four defects shipped in the path between a fit and what a user
# actually sees, and none of them was a wrong number. Each produced a
# well-formed object that was wrong in a way only visible on use -- a
# specification quietly rebuilt as a different model, a data frame whose
# `value` column became a list, a print method that emitted thousands of rows
# of draws, an argument that was documented but never wired up.
#
# The tests below are the invariants that would have caught them, and they are
# deliberately structural rather than numeric: the failures were all "right
# shape, wrong content" or "wrong shape, no error", which value assertions on
# the engine cannot see.

.laplace_nested_fit <- function() .laplace_cached("nested", suppressMessages(
  ctFit(.laplace_nested_data(), .laplace_nested_model(), backend = "julia",
    intoverpop = "laplace", optimcontrol = list(finishsamples = 100))))

test_that("re-preparing a specification preserves what the fit is", {
  skip_without_julia()
  fit <- .laplace_nested_fit()
  fitted <- ctsem:::.ctBackendSpec(fit)
  # Any non-default `subjects` or `timestep` makes `.ctBackendKalmanSpec`
  # rebuild the specification from the data rather than reuse the fit's. That
  # rebuild is the one step in prediction that can disagree with the fit, and
  # it silently did: `.ctJuliaPrepare` defaults to augmenting the state, so a
  # Laplace fit came back as a different model.
  rebuilt <- ctsem:::.ctBackendKalmanSpec(fit, subjects = 1:5)

  expect_false(is.null(rebuilt$laplace))
  expect_equal(ctsem:::.ctBackendIntOverPop(rebuilt),
    ctsem:::.ctBackendIntOverPop(fitted))
  expect_equal(rebuilt$nlatent_augmented, rebuilt$nlatent)

  # Group membership is allowed to change -- five subjects are not thirty --
  # but nothing describing the *model* may.
  for (l in seq_along(fitted$laplace$levels)) {
    for (field in c("name", "re_index", "sd_index", "cor_index", "sd_scale",
        "param", "nrandom")) {
      expect_equal(rebuilt$laplace$levels[[l]][[field]],
        fitted$laplace$levels[[l]][[field]],
        info = paste("level", l, field))
    }
  }
  expect_equal(rebuilt$laplace$npar, fitted$laplace$npar)
})

test_that("a prediction does not depend on how it was asked for", {
  skip_without_julia()
  fit <- .laplace_nested_fit()
  # The whole-fit path returns the fit's own specification; naming subjects
  # rebuilds one. Asking the same question two ways has to give one answer, and
  # this is the cheapest assertion that covers every property re-preparation
  # could drop -- present and future, without naming them.
  all <- suppressMessages(ctPredict(fit, subjects = 1:40, timestep = "asdata"))
  some <- suppressMessages(ctPredict(fit, subjects = c(2, 9, 31), timestep = "asdata"))

  keep <- all[all$Subject %in% c(2, 9, 31) & all$Element %in% "yprior", ]
  got <- some[some$Element %in% "yprior", ]
  keep <- keep[order(keep$Subject, keep$Time), ]
  got <- got[order(got$Subject, got$Time), ]
  expect_equal(nrow(got), nrow(keep))
  expect_equal(got$value, keep$value, tolerance = 1e-8)
})

test_that("the random effects used are the fitted ones, whatever is filtered", {
  skip_without_julia()
  fit <- .laplace_nested_fit()
  pars <- ctSubjectPars(fit)
  subs <- c(1, 12, 27)

  # At a subject's first row the expectation is driven by its parameters and
  # not yet by its data, so it reports the random effects the trajectory was
  # built from. Withholding observations must not change them: a mode is
  # defined relative to the data it was *estimated* from, and re-solving it
  # against withheld data collapses it to zero -- which silently turned every
  # subject into the average one, exactly when a per-subject prediction was
  # what was being asked for.
  firstrow <- function(...) {
    k <- suppressMessages(ctPredict(fit, subjects = subs, ...))
    v <- k[k$Element %in% "yprior", ]
    v <- v[order(v$Subject, v$Time), ]
    unlist(lapply(split(v, v$Subject), function(d) d$value[1L]))
  }
  baseline <- firstrow(timestep = "asdata")
  expect_equal(firstrow(timestep = 0.5, timerange = c(0, 4)), baseline,
    tolerance = 1e-8)
  expect_equal(firstrow(timestep = "asdata", removeObs = TRUE), baseline,
    tolerance = 1e-8)
  expect_equal(firstrow(timestep = 0.5, timerange = c(0, 4), removeObs = TRUE),
    baseline, tolerance = 1e-8)

  # And they are the same effects `ctSubjectPars` reports, rather than a second
  # copy computed some other way.
  mm <- pars[1, subs, "mmean"]
  expect_equal(unname(baseline), unname(mm), tolerance = 1e-6)
})

test_that("randomEffects reaches the engine and selects a level", {
  skip_without_julia()
  fit <- .laplace_nested_fit()
  # It was documented on `ctBackendKalman` and unreachable from `ctPredict`,
  # which sends `...` to plot(). A test that only called the inner function
  # would have passed.
  subs <- c(1, 6, 11, 16)
  first <- function(lev) {
    k <- suppressMessages(ctPredict(fit, subjects = subs, timestep = "asdata",
      randomEffects = lev))
    v <- k[k$Element %in% "yprior", ]
    v <- v[order(v$Subject, v$Time), ]
    unname(unlist(lapply(split(v, v$Subject), function(d) d$value[1L])))
  }
  subject <- first("subject"); study <- first("study"); population <- first("population")

  # Subject effects give every subject its own; population gives them all one.
  expect_gt(stats::sd(subject), 0)
  expect_equal(population, rep(population[1L], length(population)),
    tolerance = 1e-8)
  # The study level sits between: not the subject answer, and these four
  # subjects are in four different studies, so not the population one either.
  expect_false(isTRUE(all.equal(study, subject)))
  expect_false(isTRUE(all.equal(study, population)))

  expect_error(suppressMessages(ctPredict(fit, subjects = 1,
    randomEffects = "nosuchlevel")), "not 'nosuchlevel'")
})

test_that("a prediction data frame is numeric and its plot prints", {
  skip_without_julia()
  skip_if_not_installed("ggplot2")
  fit <- .laplace_nested_fit()
  # `meltkalman` melts everything it is handed, and `ctKalmanArray` also
  # returns a list of per-subject matrices. Melting that made `value` a list
  # column, which the final rbind then spread over the whole frame. It only
  # happened with more than one subject, and `plot()` still *returned* a
  # ggplot -- only printing it failed -- so nothing noticed.
  for (subs in list(1L, c(1L, 2L, 3L))) {
    k <- suppressMessages(ctPredict(fit, subjects = subs, timestep = "asdata"))
    expect_true(is.numeric(k$value),
      info = paste(length(subs), "subject(s)"))
    expect_false(is.list(k$value))
    p <- plot(k)
    expect_s3_class(p, "ggplot")
    # Rendering, not merely constructing: a list `value` column builds a
    # ggplot happily and only fails when something tries to draw it. Warnings
    # are expected -- interpolated rows carry no observation.
    tmp <- tempfile(fileext = ".png")
    expect_error(suppressWarnings(suppressMessages(
      ggplot2::ggsave(tmp, p, width = 5, height = 3))), NA)
    expect_true(file.exists(tmp))
    unlink(tmp)
  }
})

test_that("a multilevel summary prints its tables and not its payload", {
  skip_without_julia()
  fit <- .laplace_nested_fit()
  s <- summary(fit)
  # `$randomEffects` is the machine-readable per-level payload. It stayed
  # reachable, and the printed summary stayed a summary: printing the draws put
  # thousands of lines ahead of the tables that summarise them.
  expect_false(is.null(s$randomEffects))
  expect_equal(length(s$randomEffects), 2L)

  out <- utils::capture.output(print(s))
  expect_lt(length(out), 200L)
  expect_true(any(grepl("Random-effects standard deviations (subject)", out,
    fixed = TRUE)))
  expect_true(any(grepl("Random-effects standard deviations (study)", out,
    fixed = TRUE)))
  # Asked for by name it is still printable, which is what keeps the exclusion
  # a display choice rather than a removal.
  expect_gt(length(utils::capture.output(
    print(s, sections = "randomEffects"))), 10L)
})

test_that("optimising without integrating the random effects is refused", {
  skip_without_julia()
  # `intoverpop=FALSE` leaves each subject's random effects as free parameters
  # of the objective, so maximising it maximises over those effects too and the
  # population variance it lands on is whatever makes them most likely --
  # zero, or as near as the data allow. Sampling is what handles that model,
  # which is why `intoverpop='auto'` resolves to FALSE exactly when `optimize`
  # is. Refused rather than silently corrected, because both readings of the
  # request are plausible and guessing answers a different question.
  set.seed(3)
  d <- do.call(rbind, lapply(1:8, function(i) data.frame(id = i,
    time = 0:4, Y1 = stats::rnorm(5))))
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix("cint"), T0MEANS = matrix(0)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  expect_error(suppressWarnings(suppressMessages(
    ctFit(d, m, backend = "julia", cores = 1, intoverpop = FALSE,
      optimcontrol = list(estonly = TRUE)))),
    "intoverpop=FALSE with optimize=TRUE")
  # The same model with the effects integrated is fine, which is what makes
  # the refusal a guard rather than a limitation.
  expect_error(suppressWarnings(suppressMessages(
    ctFit(d, m, backend = "julia", cores = 1, intoverpop = "augmented",
      optimcontrol = list(estonly = TRUE)))), NA)
})
