# How much of each random effect a subject's own data determine
# (R/ctBackendEffectInformation.R, `ctsem_effect_information` in the engine).
#
# What is pinned here: that the two routes' computations -- the Laplace unit
# curvature on one side, the filter's smoothed covariance on the other --
# agree where both are exact, which is the check neither can supply for
# itself; that a population sd come out near zero is measured again at the
# starting spread, and that this is what decides; that a fit takes it at its
# estimate and keeps it through every later rebuild of the
# identifiability report; that a weakly informed effect is said once, in the
# words `summary()` and `ctReport()` repeat; and that an informed one is not
# said at all.

# A random CINT and a random initial level: one static carrier and one
# dynamic state on the augmented route, so its backward pass is exercised for
# both kinds of effect.
.ei_model <- function(random = c("cint", "t0m")) {
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix("cint"), T0MEANS = matrix("t0m"),
    silent = TRUE)))
  model$pars$indvarying <- model$pars$param %in% random
  model
}

# Simulated here rather than by ctGenerate, whose draw stream moves.
.ei_data <- function(nsubjects = 20L, nobs = 6L, noise = 0.3) {
  set.seed(20260927)
  drift <- -0.5; diffusion <- 0.6
  decay <- exp(drift)
  do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    cint <- stats::rnorm(1, 0.5, 0.4)
    state <- stats::rnorm(1, 1, 0.8)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) state <- decay * state + cint * (decay - 1) / drift +
        stats::rnorm(1, 0, sqrt(diffusion^2 * (decay^2 - 1) / (2 * drift)))
      out[t] <- state + stats::rnorm(1, 0, noise)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
}

# The data the fits use: enough measurement error that MANIFESTVAR stays
# identified. At the default 0.3 with six occasions it collapses towards
# zero, and the fit then spends a minute on a Hessian it cannot finish --
# with or without the check, so not what these tests are about.
.ei_fitdata <- function() .ei_data(nsubjects = 30L, nobs = 8L, noise = 0.5)

.ei_specs <- function(dat = .ei_data()) {
  lapply(c(augmented = "augmented", laplace = "laplace"), function(route)
    suppressWarnings(suppressMessages(ctFit(dat, .ei_model(), backend = "julia",
      intoverpop = route, priors = FALSE, fit = FALSE))))
}

# The same point on both routes, by name: they lay the raw vector out
# differently but name it the same way, and the population scales and the
# correlation have the same transforms. Diffusion and measurement error near
# the simulated ones (raw 0.1 would put them at 8 and 4, which drowns every
# subject's data and makes every effect weakly informed).
.ei_point <- function(spec, popsd_cint = -0.2, popsd_t0m = 0.4, cor = 0.3) {
  npar <- ctsem:::.ctBackendNpar(spec)
  names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
  values <- stats::setNames(rep(0.1, npar), names)
  values[grepl("^diff", names)] <- -2
  values[grepl("^mvar", names)] <- -2
  values[names == "popsd_cint"] <- popsd_cint
  values[names == "popsd_t0m"] <- popsd_t0m
  values[grepl("^rawcor_", names)] <- cor
  values
}

test_that("the two routes agree on how much each subject's data determine, where both are exact", {
  skip_without_julia()
  # Linear and Gaussian, so the Laplace term is the exact marginal and the
  # filter's smoothed covariance the exact posterior: the two are independent
  # derivations of the same number, at any parameter vector.
  specs <- .ei_specs()
  points <- lapply(specs, .ei_point)
  expect_setequal(names(points$augmented), names(points$laplace))
  # The same model at the same point, or the comparison below means nothing.
  # To 1e-6 relative, not machine precision: at this point, with small
  # diffusion and measurement error, the two objectives differ by 1.5e-4 of
  # 724 (2e-7), where at raw 0.1 for both they agreed to 1e-8. A different
  # model would differ by far more; the shares below agree to 1e-6.
  value <- function(spec, values) {
    module <- ctsem:::.ctJuliaModule(spec$project)
    as.numeric(JuliaConnectoR::juliaGet(module$ctsem_evaluate(
      ctsem:::.ctJuliaObjective(spec), ctsem:::.ctJuliaNumericVector(values),
      gradient = FALSE))$value)
  }
  expect_equal(value(specs$laplace, points$laplace),
    value(specs$augmented, points$augmented), tolerance = 1e-6)

  info <- Map(function(spec, values) ctsem:::.ctEffectInformation(spec,
    values, point = "a test point"), specs, points)
  laplace <- info$laplace$groups[[1L]]
  augmented <- info$augmented$groups[[1L]]
  expect_setequal(colnames(laplace), c("cint", "t0m"))
  expect_setequal(colnames(augmented), c("cint", "t0m"))
  expect_equal(dim(laplace), c(20L, 2L))
  # Elementwise, per subject and effect.
  expect_equal(augmented[, colnames(laplace)], laplace, tolerance = 1e-6)
  # And a share, not something else of the right shape: between zero and one
  # for a likelihood concave in the effects.
  expect_true(all(laplace > 0 & laplace < 1))
  expect_identical(info$laplace$point, "a test point")
  expect_equal(info$laplace$values, unname(points$laplace))
  expect_equal(info$laplace$units$n, 20L)
})

test_that("an sd near zero is measured again at the starting spread, and that decides", {
  skip_without_julia()
  # The AnomAuth shape in miniature: a population sd at the bottom of its
  # transform, where the share at the point is zero whatever the data say. Six
  # occasions of a Gaussian indicator pin each subject's intercept, so at a
  # spread large enough to see, the share is high: no individual differences,
  # seen clearly, which is a finding and not a weakly informed effect.
  # Uncorrelated, so that nothing of the intercept is learned through the
  # initial level: with a correlation, the share of a collapsed effect is
  # the share of it the other effect explains, not zero.
  specs <- .ei_specs()
  info <- lapply(specs, function(spec) {
    values <- .ei_point(spec, popsd_cint = -6, cor = 0)
    list(values = values, record = ctsem:::.ctEffectInformation(spec, values,
      point = "a test point"))
  })
  for (route in names(info)) {
    table <- info[[route]]$record$table
    cint <- table[table$effect == "cint", ]
    t0m <- table[table$effect == "t0m", ]
    expect_lt(cint$determined, 1e-6)
    expect_gt(cint$reference, 0.5)
    expect_false(cint$weak)
    # At raw zero, the spread every fit starts from: log1p_exp(-1).
    expect_equal(cint$referencesd, log1p(exp(-1)), tolerance = 1e-6)
    # Only an sd below the starting spread is measured again.
    expect_true(is.na(t0m$reference))
    # And the second point is the first with that one coordinate at zero.
    values <- info[[route]]$values
    values[names(values) == "popsd_cint"] <- 0
    again <- ctsem:::.ctEffectInformation(specs[[route]], values,
      point = "a test point")$table
    expect_equal(cint$reference, again$determined[again$effect == "cint"])
  }
  # Both routes find the same starting-spread share, where both are exact.
  shares <- vapply(info, function(x) x$record$table$reference[
    x$record$table$effect == "cint"], numeric(1))
  expect_equal(shares[["augmented"]], shares[["laplace"]], tolerance = 1e-6)
  # An sd below the starting spread whose share already clears the bar is
  # not measured again: raising it could not change the verdict.
  modest <- ctsem:::.ctEffectInformation(specs$laplace,
    .ei_point(specs$laplace, popsd_cint = -0.2), point = "a test point")$table
  cint <- modest[modest$effect == "cint", ]
  expect_gt(cint$determined, ctsem:::.ctEffectThresholds()$determined)
  expect_true(is.na(cint$reference))
  expect_false(cint$weak)
})

test_that("AnomAuth's random drift is weakly informed and its random intercept is not", {
  skip_without_julia()
  # The case the check exists for, on the data it was calibrated on: anomia,
  # subjects with three to five observed waves -- the bench's AnomAuth cells
  # (dev/optimbench/cells.R, gg_genS), on the first 100 of those subjects.
  e <- new.env()
  utils::data("AnomAuth", package = "ctsem", envir = e)
  long <- suppressMessages(ctWideToLong(e$AnomAuth, Tpoints = 5,
    n.manifest = 2, manifestNames = c("Y1", "Y2")))
  long <- as.data.frame(suppressMessages(ctDeintervalise(long)))
  long <- long[!is.na(long$Y1), c("id", "time", "Y1")]
  waves <- table(long$id)
  keep <- utils::head(as.numeric(names(waves)[waves >= 3]), 100)
  long <- long[long$id %in% keep, ]
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "anom", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix("cint||TRUE"),
    DRIFT = matrix("drift|-log1p_exp(-param)|TRUE"), T0MEANS = matrix("t0m"),
    silent = TRUE)))
  model$pars$indvarying <- model$pars$param %in% c("drift", "cint")
  spec <- suppressWarnings(suppressMessages(ctFit(long, model,
    backend = "julia", intoverpop = "laplace", fit = FALSE)))
  # The best-known point of the bench's 800-subject cell S1
  # (dev/optimbench/starts.R, hist_anomS1), by name: both population sds at
  # the bottom of their transform, which is where fits of these data end.
  npar <- ctsem:::.ctBackendNpar(spec)
  names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
  best <- c(t0m = 0.2584, drift = 3.9345, diff_anom = -2.0670,
    mvarY1 = -1.1698, cint = 0.0071, T0var_anom = -0.9870,
    popsd_drift = -4.4930, popsd_cint = -6.6331, rawcor_cint__drift = 0)
  expect_setequal(names, names(best))
  record <- ctsem:::.ctEffectInformation(spec, best[names], point = "x")
  table <- record$table
  drift <- table[table$effect == "drift", ]
  cint <- table[table$effect == "cint", ]
  # Each subject's three to five waves pin its intercept, and say almost
  # nothing about its rate of change.
  expect_lt(drift$reference, 0.05)
  expect_true(drift$weak)
  expect_gt(cint$reference, 0.9)
  expect_false(cint$weak)
  advice <- ctsem:::.ctEffectAdvice(record)
  expect_length(advice, 1L)
  expect_match(advice, "Individual differences in drift are barely informed",
    fixed = TRUE)
  expect_match(advice, "indvarying = FALSE for drift", fixed = TRUE)
})

test_that("a fit takes it at its estimate, keeps it, and says nothing when every effect is informed", {
  skip_without_julia()
  messages <- character()
  # The random CINT alone on this route: with a random initial level as well,
  # the Laplace Hessian's inner solves stop short of their tolerance on these
  # data (with or without the check), and the fit spends a minute on it.
  fit <- withCallingHandlers(suppressWarnings(ctFit(.ei_fitdata(),
    .ei_model("cint"), backend = "julia", intoverpop = "laplace", cores = 1,
    verbose = 0)),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  effects <- fit$identifiability$effects
  expect_identical(effects$point, "at the estimate")
  expect_equal(effects$values, fit$estimate$raw)
  expect_identical(effects$route, "laplace")
  expect_identical(effects$table$effect, "cint")
  expect_false(any(effects$table$weak))
  expect_false(any(grepl("barely informed", messages)))
  # The uncertainty stage and the fit's last step both rebuild the report from
  # a new curvature; neither drops what the check found.
  again <- suppressWarnings(suppressMessages(ctFitUncertainty(fit, "hessian")))
  expect_identical(again$identifiability$effects, effects)
  # Nothing to say, so the summary says nothing about it.
  printed <- utils::capture.output(print(summary(fit)))
  expect_false(any(grepl("barely informed", printed)))
})

test_that("a weakly informed effect is said once, and summary() and ctReport() repeat it", {
  skip_without_julia()
  # Every effect called weak, so the plumbing is tested on a fit that is cheap
  # and well behaved rather than on one that is weak in fact; which effects
  # the real rule calls weak is calibrated on the bench, not here.
  testthat::local_mocked_bindings(
    .ctEffectThresholds = function() list(determined = 1.01), .package = "ctsem")
  messages <- character()
  fit <- withCallingHandlers(suppressWarnings(ctFit(.ei_fitdata(), .ei_model(),
    backend = "julia", intoverpop = "augmented", cores = 1, verbose = 0,
    optimcontrol = list(estonly = TRUE))),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  said <- grep("barely informed", messages, value = TRUE)
  expect_length(said, 1L)
  expect_match(said, "Individual differences in cint are barely informed",
    fixed = TRUE)
  expect_match(said, "Consider indvarying = FALSE for t0m", fixed = TRUE)
  expect_true(all(fit$identifiability$effects$table$weak))
  printed <- paste(utils::capture.output(print(summary(fit))), collapse = " ")
  expect_match(printed, "barely informed", fixed = TRUE)
  expect_match(printed, "fit$identifiability$effects", fixed = TRUE)
  # The identification component of ctReport(), which is where it is printed.
  report <- ctsem:::.ctReportIdentificationLines(fit, summary(fit))
  expect_true(any(grepl("Random-effect information, at the estimate",
    report, fixed = TRUE)))
  expect_true(any(grepl("barely informed", report, fixed = TRUE)))
})

test_that("ctIdentify() takes it at supplied inits only", {
  skip_without_julia()
  dat <- .ei_data()
  model <- .ei_model()
  plain <- suppressWarnings(suppressMessages(ctIdentify(dat, model,
    intoverpop = "laplace", nstart = 1L)))
  expect_null(plain$effects)
  spec <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", priors = FALSE, fit = FALSE)))
  at <- .ei_point(spec)
  given <- suppressWarnings(suppressMessages(ctIdentify(dat, model,
    intoverpop = "laplace", inits = unname(at))))
  expect_identical(given$effects$point, "at the supplied inits")
  expect_equal(given$effects$values, unname(at))
  expect_setequal(given$effects$table$effect, c("cint", "t0m"))
  expect_true(all(is.finite(given$effects$table$determined)))
  printed <- utils::capture.output(print(given))
  expect_true(any(grepl("informed by its own group's data", printed,
    fixed = TRUE)))
})

test_that("the wording names the level, its switch, and which spread it measured", {
  effects <- list(table = data.frame(level = c("id", "study", "id"),
    effect = c("drift", "cint", "t0m"), popsd = c(0.1, 0.2, 1e-8),
    groups = c(40L, 8L, 40L), determined = c(0.6, 0.004, 1e-9),
    reference = c(NA, NA, 0.02), referencesd = c(NA, NA, 0.313),
    widened = c(0, 0.5, 0), unit = c("subject", "study", "subject"),
    switch = c("indvarying", "indvarying_study", "indvarying"),
    weak = c(FALSE, TRUE, TRUE), stringsAsFactors = FALSE))
  advice <- ctsem:::.ctEffectAdvice(effects)
  expect_length(advice, 2L)
  expect_match(advice[1L], "Individual differences in cint (study level)",
    fixed = TRUE)
  expect_match(advice[1L], paste0("at its estimated raw-scale population sd ",
    "of 0.2 a typical study's own data determine almost none"), fixed = TRUE)
  expect_match(advice[1L], "indvarying_study = FALSE for cint", fixed = TRUE)
  expect_match(advice[1L], "more observations per study", fixed = TRUE)
  # Measured at the starting spread, and saying so -- both sds by value.
  expect_match(advice[2L], "t0m (id level)", fixed = TRUE)
  expect_match(advice[2L], paste0("its raw-scale population sd was estimated ",
    "at 1e-08, and even at 0.31 a typical subject's own data would determine ",
    "2% of its value"), fixed = TRUE)
  effects$table$weak <- FALSE
  expect_length(ctsem:::.ctEffectAdvice(effects), 0L)
  expect_length(ctsem:::.ctEffectAdvice(NULL), 0L)
})
