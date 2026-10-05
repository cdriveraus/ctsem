# Posterior-predictive data generation for backend='julia'.
#
# Generation is not a recorder: it changes what the filter consumes, so the
# state carried forward is conditioned on the drawn data rather than the real
# data. That is what makes the result a draw from the model, and it is also what
# makes it easy to get subtly wrong -- a simulator that predicted each row from
# the *real* history would produce data that looks fine and is not a draw from
# anything.
#
# So the load-bearing test is an identity: re-running the ordinary likelihood on
# the generated dataset must reproduce the likelihood reported while generating
# it. Nothing that drew from the wrong covariance, or conditioned on the wrong
# history, can satisfy that.

.generate_model <- function() {
  suppressWarnings(ctModel(type = "ct", LAMBDA = matrix(1, 1, 1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diff", 1, 1),
    MANIFESTVAR = matrix("mvar", 1, 1), MANIFESTMEANS = matrix("mmean||FALSE", 1, 1),
    T0VAR = matrix("t0v", 1, 1), T0MEANS = matrix(0, 1, 1), CINT = matrix(0, 1, 1)))
}

.generate_data <- function() {
  set.seed(11)
  times <- c(0, .5, 1, 1.7, 2.5, 3.4)
  drift <- -0.8
  diffusion <- 0.5
  data <- do.call(rbind, lapply(seq_len(20), function(i) {
    state <- stats::rnorm(1, 0, .6)
    y <- numeric(length(times))
    for (t in seq_along(times)) {
      if (t > 1) {
        dt <- times[t] - times[t - 1]
        state <- exp(drift * dt) * state +
          stats::rnorm(1, 0, diffusion * sqrt((1 - exp(2 * drift * dt)) / (-2 * drift)))
      }
      y[t] <- state + stats::rnorm(1, 0, .3) + 1.2
    }
    data.frame(id = i, time = times, Y1 = y)
  }))
  data$Y1[c(4, 20)] <- NA
  data
}

test_that("generated data is a draw from the model the filter conditions on", {
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  spec <- fit$model_spec
  nrows <- length(spec$times)

  set.seed(1)
  base <- matrix(stats::rnorm(nrows), 1, nrows)
  drawn <- ctsem:::.ctBackendGenerate(fit, fit$estimate$raw, base)

  expect_equal(dim(drawn$Y), c(1L, nrows))
  # The identity that defines it: the filter's own likelihood for the generated
  # data is the likelihood it reported while generating it.
  refit <- data
  refit$Y1 <- as.numeric(drawn$Y)
  respec <- suppressMessages(ctFit(refit, model, backend = "julia", fit = FALSE))
  regenerated <- ctJuliaEvaluate(ctsem:::.ctBackendAsModel(respec),
    fit$estimate$raw, gradient = FALSE)$value
  expect_equal(sum(drawn$subject_loglik), regenerated, tolerance = 1e-10)
  expect_equal(sum(drawn$llrow), regenerated, tolerance = 1e-10)

  # Missingness is preserved: an entry that was not observed is not invented.
  expect_equal(which(is.na(as.numeric(drawn$Y))), which(is.na(data$Y1)))
  expect_equal(drawn$llrow[which(is.na(data$Y1))], rep(0, sum(is.na(data$Y1))))

  # The draws drive it -- different normals, different data.
  set.seed(2)
  other <- ctsem:::.ctBackendGenerate(fit, fit$estimate$raw,
    matrix(stats::rnorm(nrows), 1, nrows))
  expect_false(isTRUE(all.equal(as.numeric(other$Y), as.numeric(drawn$Y))))
})

test_that("the observed data is an unremarkable draw from the fitted model", {
  skip_without_julia()
  # The calibration check the identity above cannot make: at the maximum
  # likelihood estimate, the observed data's log likelihood should sit somewhere
  # ordinary in the distribution of generated ones. A generator that drew from
  # too small a covariance, or that conditioned each row on the real history
  # instead of the drawn one, would put the observed value far into a tail.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  nrows <- length(fit$model_spec$times)

  set.seed(2)
  generated <- replicate(200, sum(ctsem:::.ctBackendGenerate(fit, fit$estimate$raw,
    matrix(stats::rnorm(nrows), 1, nrows))$subject_loglik))
  quantile <- mean(generated < fit$estimate$loglik)
  expect_gt(quantile, .05)
  expect_lt(quantile, .95)
  # And the spread is real, not a degenerate point mass.
  expect_gt(stats::sd(generated), 1)
})

test_that("ctGenerateFromFit returns what the posterior predictive tools expect", {
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  nrows <- length(fit$model_spec$times)

  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 20, cores = 1)
  expect_equal(dim(generated$generated$Y), c(20L, nrows, 1L))
  expect_identical(dimnames(generated$generated$Y)[[3]], model$manifestNames)
  expect_equal(dim(generated$generated$llrow), c(20L, nrows))
  # A row with nothing observed contributes nothing, and is reported as NA
  # rather than as a zero that would read as a likelihood.
  expect_true(all(is.na(generated$generated$llrow[, which(is.na(data$Y1))])))
  expect_true(all(is.na(generated$generated$Y[, which(is.na(data$Y1)), 1])))

  # fullposterior needs draws, and says so rather than silently using the mode.
  # A default fit has them, since ctFit() finishes with ctOptimUncertainty(), so
  # the refusal is only reachable for a fit that deliberately skipped it.
  estonly <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0,
    optimcontrol = list(estonly = TRUE)))
  expect_error(ctGenerateFromFit(estonly, nsamples = 5, fullposterior = TRUE, cores = 1),
    "ctFitUncertainty")

  fromposterior <- ctGenerateFromFit(fit, nsamples = 10, fullposterior = TRUE,
    cores = 1)
  expect_equal(dim(fromposterior$generated$Y), c(10L, nrows, 1L))
})

test_that("the posterior predictive tools run on a backend fit", {
  skip_without_julia()
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 20, cores = 1)

  predictive <- suppressMessages(ctsem:::ctPostPredData(generated))
  expect_true(all(c("row", "variable", "sample", "value", "id", "Time",
    "TimeInterval", "obsValue") %in% names(predictive)))
  # Every row of every sample, for each manifest plus the row likelihood.
  expect_equal(nrow(predictive), 20 * length(fit$model_spec$times) * 2)
  # The observed column is the data, not a copy of the generated one.
  observed <- predictive[predictive$variable == "Y1" & predictive$sample == 1, ]
  expect_equal(observed$obsValue[order(observed$row)], data$Y1)

  plots <- suppressWarnings(suppressMessages(ctPostPredPlots(generated)))
  expect_true(length(plots) > 0)
  expect_true(all(vapply(plots, function(p) inherits(p, "ggplot"), logical(1))))
  # Both panel families, from the merged function.
  expect_true(all(c("Density", "PIT", "CalibrationBySubject") %in% names(plots)))
  # Subject-level log likelihood: totals per subject, compared against the same
  # curve from each generated dataset. Both sides must total the same rows.
  expect_true("SubjectLogLik" %in% names(plots))
  ll <- predictive[as.character(predictive$variable) == "LogLik" &
    is.finite(predictive$obsValue) & is.finite(predictive$value), ]
  gensub <- ll[, .(tot = sum(value)), by = .(sample, id)]
  obssub <- unique(ll[, .(row, id, obsValue)])[, .(tot = sum(obsValue)), by = id]
  expect_equal(nrow(obssub), length(unique(ll$id)))
  expect_equal(nrow(gensub), nrow(obssub) * length(unique(ll$sample)))
  expect_true(all(is.finite(gensub$tot)) && all(is.finite(obssub$tot)))
  expect_true(any(grepl("^ChangeByValue_", names(plots))))
  # Time intervals run forwards. They did not: the generated table was ordered
  # by row label rather than row number, so diff(Time) crossed subject
  # boundaries for every subject whose rows span a digit-count boundary.
  expect_true(all(predictive$TimeInterval >= 0, na.rm = TRUE))
  # Notes are captions, and can be turned off.
  bare <- suppressWarnings(suppressMessages(
    ctPostPredPlots(generated, panels = "PIT", notes = FALSE)))
  expect_null(bare$PIT$labels$caption)
  expect_false(is.null(plots$PIT$labels$caption))
  # datarows selects rows rather than erroring, which it used to.
  sub <- suppressWarnings(suppressMessages(
    ctPostPredPlots(generated, panels = "PIT", datarows = 5:40)))
  expect_true(inherits(sub$PIT, "ggplot"))

  covcheck <- suppressWarnings(suppressMessages(
    ctFitCheckCov(generated, plot = FALSE, lags = 0:2, cores = 1, nsamples = 10)))
  expect_true(nrow(covcheck) > 0)
  expect_true("Sig" %in% names(covcheck))
})

test_that("ctPostPredData(residuals=TRUE) works, for stan too", {
  skip_without_julia()
  # This branch could never run: the residual rows lacked the id/time columns
  # the rbind below them needs, on every backend including stan. Fixed here, so
  # guarded here.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 3, cores = 1)

  predictive <- suppressMessages(ctsem:::ctPostPredData(generated, residuals = TRUE))
  expect_true("Y1 std. res." %in% predictive$variable)
  expect_false(anyNA(predictive$Time[predictive$variable == "Y1 std. res."]))

  skip_if_not_installed("rstan")
  stangenerated <- suppressMessages(ctGenerateFromFit(ctstantestfit, nsamples = 2,
    cores = 1))
  stanpredictive <- suppressMessages(ctsem:::ctPostPredData(stangenerated,
    residuals = TRUE))
  expect_true(any(grepl("std. res.", stanpredictive$variable)))
})

# User side generation draws each subject from the population the model
# states, at every level, and simulates it. The values it has to choose are
# the priors' scales: raw zero for a free mean, sdscale for a population sd.
.gen_model <- function(id = "id", ti = 0) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = id,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix("drift"), DIFFUSION = matrix(.5), MANIFESTVAR = matrix(.3),
    MANIFESTMEANS = matrix("mm"), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), n.TIpred = ti,
    TIpredNames = if (ti) paste0("TI", seq_len(ti)))))
  m$pars$indvarying <- m$pars$param %in% "mm"
  m
}
.gen_subject_means <- function(d) tapply(d[, "Y1"], d[, "id"], mean)

test_that("n counts each id level, innermost first", {
  m <- .gen_model(id = c("id", "study"))
  expect_equal(ctsem:::.ctGenerateLevelCounts(m, c(12, 3)), c(id = 12L, study = 3L))
  expect_equal(ctsem:::.ctGenerateLevelCounts(m, 12), c(id = 12L, study = 1L))
  expect_error(ctsem:::.ctGenerateLevelCounts(m, c(3, 12)), "must not grow")
  expect_error(ctsem:::.ctGenerateLevelCounts(m, c(12, 3, 1)), "2 id levels")
  skeleton <- ctsem:::.ctGenerateSkeleton(m, c(12, 3), lapply(1:12, function(i) 0:1))
  expect_equal(as.numeric(table(skeleton$study[!duplicated(skeleton$id)])), c(4, 4, 4))
  expect_error(ctGenerate(m, n = c(12, 3), Tpoints = 2, backend = "r"),
    "n must be one count")
  # The former name still works, and says so.
  flat <- .gen_model()
  set.seed(1)
  expect_warning(old <- suppressMessages(ctGenerate(flat, n.subjects = 3,
    Tpoints = 2, backend = "r")), "is deprecated; use n", fixed = TRUE)
  set.seed(1)
  expect_equal(old, suppressMessages(ctGenerate(flat, n = 3, Tpoints = 2,
    backend = "r")))
})

test_that("ctGenerate draws subjects from the population at every level", {
  skip_without_julia()
  m <- .gen_model(id = c("id", "study"))
  m$pars$indvarying_study <- m$pars$param %in% "mm"
  m$pars$sdscale_study[m$pars$param %in% "mm"] <- 0.5
  set.seed(3)
  expect_message(d <- ctGenerate(m, n = c(600, 30), Tpoints = 6,
    backend = "julia"), "mm 1 \\[id\\], mm 0.5 \\[study\\]")
  means <- .gen_subject_means(d)
  study <- tapply(d[, "study"], d[, "id"], `[`, 1)
  # mm = 10 * param, so raw sds of 1 and 0.5 are 10 and 5 in the data. A study
  # mean also carries its 20 subjects' spread, 10 / sqrt(20).
  expect_equal(mean(tapply(means, study, sd)), 10, tolerance = 0.1)
  expect_equal(sd(tapply(means, study, mean)), sqrt(25 + 100 / 20), tolerance = 0.25)
})

test_that("popmeans, RAWPOPVAR and fixed TI effects set what they state", {
  skip_without_julia()
  m <- .gen_model(ti = 1)
  m$RAWPOPVAR["mm", "mm"] <- 0.2
  m$pars$TI1_effect <- ifelse(m$pars$param %in% "mm", "0.3", "FALSE")
  set.seed(4)
  expect_message(d <- ctGenerate(m, n = 400, Tpoints = 8,
    backend = "julia", popmeans = c(mm = 3, drift = -0.5)),
    "from popmeans: mm, drift; population sd \\(raw\\): mm 0.2 \\[id\\]")
  means <- .gen_subject_means(d)
  ti <- tapply(d[, "TI1"], d[, "id"], `[`, 1)
  regression <- stats::lm(means ~ ti)
  # A fixed effect of 0.3 shifts the raw parameter, 3 in the data; the
  # residual spread is RAWPOPVAR's 0.2 raw, 2 in the data.
  expect_equal(unname(stats::coef(regression)), c(3, 3), tolerance = 0.1)
  expect_equal(stats::sigma(regression), 2, tolerance = 0.15)

  free <- .gen_model(ti = 1)
  free$pars$TI1_effect <- ifelse(free$pars$param %in% "mm", "TRUE", "FALSE")
  set.seed(4)
  expect_message(ctGenerate(free, n = 4, Tpoints = 3, backend = "julia"),
    "at prior centres \\(raw 0\\): drift=-1.386, mm=0.*TI effects at zero: TI1 on mm")
  expect_error(ctGenerate(free, n = 4, Tpoints = 3, backend = "julia",
    popmeans = c(nothere = 1)), "not a free parameter")
})

test_that("state dependent generation carries the dependence into the data", {
  skip_without_julia()
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("y1", "y2"), latentNames = "eta1",
    LAMBDA = matrix(c(1, "0.9 + 0.35 * eta1"), 2, 1),
    MANIFESTVAR = matrix(c(0.3, 0, 0, 0.3), 2, 2),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(1.5), T0VAR = matrix(1.7),
    T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(c(0, 0), 2, 1), Tpoints = 8)))
  set.seed(9)
  d <- data.frame(suppressMessages(ctGenerate(m, n = 40, Tpoints = 8,
    backend = "julia")))

  # The loading on y2 rises with the state, so regressing y2 on y1 where y1 is
  # high must give a steeper slope than where it is low. A filled-in constant
  # loading gives the same slope in both halves.
  high <- d$y1 > stats::median(d$y1)
  slope <- function(rows) unname(stats::coef(stats::lm(y2 ~ y1, d[rows, ]))[2])
  expect_gt(slope(high), slope(!high) + 0.2)
})

test_that("ctPostPredict() runs on a julia backend fit", {
  skip_without_julia()
  # It used to refuse one, because it read `standata` and `data$Y` directly. It
  # is now an alias for ctPostPredPlots(), which goes through the
  # backend-agnostic accessors, so the refusal is gone.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  generated <- suppressMessages(ctGenerateFromFit(fit, nsamples = 10, cores = 1))
  plots <- suppressWarnings(suppressMessages(
    ctPostPredict(generated, plot = FALSE, panels = "calibration")))
  expect_true(length(plots) > 0)
  expect_true(all(vapply(plots, function(p) inherits(p, "ggplot"), logical(1))))
})

test_that("ctGenerateFromPriors() generates over a julia backend fit's own design", {
  skip_without_julia()
  # A julia fit used to be refused, because it carried only the form
  # .ctModelIntOverPop() produced, and re-preparing from that would augment it
  # twice. It now keeps the model it was given, so it is unwrapped as a stan
  # fit is: that model, over the fit's own subjects and times.
  model <- .generate_model()
  data <- .generate_data()
  fit <- suppressMessages(ctFit(data, model, backend = "julia", verbose = 0))
  # The julia default is a prior on the random-effect correlations only, and
  # the draws use every parameter's prior, so that is said.
  expect_warning(out <- suppressMessages(ctsem:::.ctGenerateFromPriors(fit,
    nsamples = 3, cores = 1, backend = "julia")), regexp = "randomCorr")
  expect_equal(dim(out$Y), c(3L, nrow(data), 1L))
  expect_true(all(is.finite(out$Y)))

  # A julia fit made before fits kept their model has only the augmented form,
  # and is still refused, with the way past it.
  old <- fit
  old[c("modelbase", "ctstanmodelbase")] <- NULL
  expect_null(ctsem:::.ctFitBaseModel(old))
  expect_error(ctsem:::.ctGenerateFromPriors(old), regexp = "does not carry")
  expect_error(ctsem:::.ctGenerateFromPriors(old), regexp = "Pass the model instead")
})

# Generation on the Laplace random-effect route (intoverpop='laplace') ------
#
# `ctGenerateFromFit()` used to error here: the engine's `ctsem_generate`
# refused any `CTSEMLaplaceObjective` with "not implemented ... yet". The
# augmented route (intoverpop='augmented') carries individual differences as
# augmented latent states the filter integrates over, so one shared raw
# parameter vector is enough to generate from. The Laplace route has no such
# vector -- each subject's random effect is a conditional mode estimated from
# that subject's own data -- so generation needs the per-subject parameter
# matrix `ctsem_kalman(laplace, ...)` already computes via
# `ctsem_laplace_subject_values`, and `ctsem_generate` did not accept one.
#
# A generator that silently fell back to the population vector for every
# subject would produce data that looks entirely reasonable -- smooth,
# correctly scaled, correctly missing -- and would simply have no individual
# differences in it. That is the specific wrong answer these tests are aimed
# at: not "does it error", but "does each subject's own random effect
# actually reach the generated data".

# The fixture is `laplace_fixture()` (helper-julia.R), the fit
# test-julia-sample.R samples from: one latent, a random manifest mean, the
# Laplace route. It was written out here and fitted twice, identically.
test_that("ctGenerateFromFit works on a Laplace fit and matches the augmented route", {
  skip_without_julia()
  data <- laplace_fixture_data()
  model <- laplace_fixture_model()
  nsub <- length(unique(data$id))
  tp <- sum(data$id == data$id[1])

  fit_laplace <- laplace_fixture()
  fit_augmented <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, intoverpop = "augmented", priors = TRUE)))

  # New subjects from the fitted population on both routes, the default.
  set.seed(123)
  gen_laplace <- ctGenerateFromFit(fit_laplace, nsamples = 20)
  set.seed(123)
  gen_augmented <- ctGenerateFromFit(fit_augmented, nsamples = 20)
  expect_error(ctGenerateFromFit(fit_augmented, nsamples = 2, effects = "fitted"),
    "intoverpop = 'augmented'")

  # Same shape as the other three routes: same names, same dimensions, same
  # class, so downstream tools cannot tell which route produced the fit.
  expect_s3_class(gen_laplace, "ctJuliaFit")
  expect_s3_class(gen_laplace, "ctFit")
  expect_identical(class(gen_laplace), class(gen_augmented))
  expect_identical(names(gen_laplace$generated), names(gen_augmented$generated))
  expect_identical(dim(gen_laplace$generated$Y), dim(gen_augmented$generated$Y))
  expect_identical(dimnames(gen_laplace$generated$Y)[[3]], model$manifestNames)
  expect_identical(dim(gen_laplace$generated$llrow), dim(gen_augmented$generated$llrow))

  # The two routes fit the same model to the same data, so new subjects from
  # either should look alike: the same grand mean and the same spread of
  # subject means. Compared per dataset, against the standard error the 20
  # datasets themselves give -- rows within a subject are dependent, so a test
  # that treats every value as independent would be far too strict.
  row_subject <- rep(seq_len(nsub), each = tp)
  per_dataset <- function(gen, f) apply(gen$generated$Y[, , 1], 1, f)
  closeness <- function(f) {
    l <- per_dataset(gen_laplace, f)
    a <- per_dataset(gen_augmented, f)
    abs(mean(l) - mean(a)) / sqrt(stats::var(l) / length(l) + stats::var(a) / length(a))
  }
  expect_lt(closeness(mean), 4)
  expect_lt(closeness(function(y) stats::sd(tapply(y, row_subject, mean))), 4)

  # effects = 'fitted' keeps the fitted subjects, each at its own conditional
  # mode, so a subject's generated mean tracks its observed one; the new
  # subjects of the default do not. A generator that put every subject at the
  # population vector would fail the first.
  obs_subject_mean <- tapply(data$Y1, data$id, mean)
  subject_means <- function(gen) {
    gen_y <- gen$generated$Y[, , 1]
    vapply(seq_len(nsub), function(s) mean(gen_y[, row_subject == s]), numeric(1))
  }
  set.seed(123)
  gen_fitted <- ctGenerateFromFit(fit_laplace, nsamples = 20, effects = "fitted")
  expect_gt(stats::cor(obs_subject_mean, subject_means(gen_fitted)), 0.8)
  # Averaged over 20 datasets of new subjects, each subject's mean is near the
  # population's, whatever that subject's data said.
  expect_lt(abs(stats::cor(obs_subject_mean, subject_means(gen_laplace))), 0.6)
})

test_that("the posterior predictive tools run on a Laplace backend fit", {
  skip_without_julia()
  fit <- laplace_fixture()

  set.seed(7)
  generated <- ctGenerateFromFit(fit, nsamples = 10, cores = 1)

  predictive <- suppressMessages(ctsem:::ctPostPredData(generated))
  expect_true(all(c("row", "variable", "sample", "value", "id", "Time",
    "TimeInterval", "obsValue") %in% names(predictive)))

  plots <- suppressWarnings(suppressMessages(ctPostPredPlots(generated)))
  expect_true(length(plots) > 0)

  # The residual branch re-filters the generated data conditional on modes
  # re-estimated from it -- the same design ctPostPredData() already uses on
  # the other three routes -- so it should run, not merely the default path.
  withresiduals <- suppressMessages(ctsem:::ctPostPredData(generated, residuals = TRUE))
  expect_true("Y1 std. res." %in% withresiduals$variable)

  # The two structure panels that came from the removed ctFitCheck() dashboard.
  # LaggedCovariance goes through ctFitCheckCov(), which is documented to work
  # on a julia fit; confirm both do once $generated is populated from the
  # Laplace route.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  panels <- suppressWarnings(suppressMessages(ctPostPredPlots(generated,
    panels = c("MeanTrajectory", "LaggedCovariance"), lags = 0:2)))
  expect_true("MeanTrajectory" %in% names(panels))
  expect_true(any(grepl("^LaggedCovariance_", names(panels))))
})

test_that("missingness survives Laplace-route generation unchanged", {
  skip_without_julia()
  data <- laplace_fixture_data()
  data$Y1[c(3, 40)] <- NA
  model <- laplace_fixture_model()
  fit <- suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    cores = 1, intoverpop = "laplace", priors = TRUE,
    optimcontrol = list(finishsamples = 20))))

  generated <- ctGenerateFromFit(fit, nsamples = 5, cores = 1)
  missing_rows <- which(is.na(data$Y1))
  expect_true(all(is.na(generated$generated$Y[, missing_rows, 1])))
  expect_true(all(is.na(generated$generated$llrow[, missing_rows])))
  expect_false(any(is.na(generated$generated$Y[, -missing_rows, 1])))
})
