# ctFitUpdate() rebuilds a fit from its own call: the arguments it was made
# with, the model it was given and the data it was fitted to, or new data. Two
# ways that goes wrong without anything erroring, which is why most of what
# follows compares the rebuilt object against the original rather than only
# checking that a call returns:
#
#   * the call is replayed with a value that means something different from
#     what the fit ran with -- `priors = TRUE` for a julia fit made under
#     `'randomCorr'` is a prior on every coordinate, a different estimator, and
#     the refit would have moved to it without a word;
#   * a field of the updated fit is left describing the old data while the one
#     beside it describes the new.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

.quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# The warning for the deprecated `iter`/`chains`/`control` spellings is said
# once per session, so a test that asserts it is not drawn has to start from a
# session where it has not been.
.local_fresh_sample_deprecation <- function(env = parent.frame()) {
  flag <- ctsem:::.ct_sample_deprecation
  said <- flag$said
  flag$said <- NULL
  withr::defer(flag$said <- said, envir = env)
}

# Stan ------------------------------------------------------------------------
#
# `ctstantestfit` is optimised, so `refit = FALSE` rebuilds its data without
# fitting -- a second or two, and no compilation.
#
# What an update should rebuild is what ctFit() builds now from the call the
# fit was made with (R/ctdataupdate.R), not the stored fit's own `$standata`,
# whose layout is whatever the format was when the fixture was last saved.
.stan_rebuilt <- function() .quietly(ctFit(ctstantestdat,
  ctstantestfit$ctstanmodelbase, priors = TRUE, fit = FALSE, backend = 'stan'))

test_that("a stan fit updates without refitting, unchanged and without warnings", {
  .local_fresh_sample_deprecation()
  updated <- NULL
  # The replayed `iter`, `chains` and `control` hold what the fit resolved
  # them to, and passed back as arguments they read as the deprecated
  # spellings -- a warning at a caller who never used them.
  suppressMessages(expect_no_warning(
    updated <- ctFitUpdate(ctstantestfit, refit = FALSE)))
  expect_s3_class(updated, "ctStanFit")
  rebuilt <- .stan_rebuilt()
  expect_equal(updated$standata, rebuilt$standata)
  expect_equal(updated$data, rebuilt$data)
})

test_that("replacement data is used, and the old data is not rebuilt first", {
  long <- ctsem:::.ctFitLongData(ctstantestfit)
  kept <- unique(long$id)[1:3]
  # `length(data==1)` was TRUE for any data set, so the old data was rebuilt on
  # every call and thrown away whenever new data came with it.
  testthat::local_mocked_bindings(standatatolong = function(...)
    stop("the old data was rebuilt"), .package = "ctsem")
  updated <- suppressMessages(ctFitUpdate(ctstantestfit,
    data = long[long$id %in% kept, ], refit = FALSE))
  expect_equal(updated$standata$nsubjects, 3L)
  # `$data` is `$standata` as a user reads it, and has to follow it.
  expect_equal(updated$data$nsubjects, 3L)
  observed <- updated$standata$Y
  observed[observed == 99999] <- NA
  expect_equal(updated$data$Y, observed)
})

test_that("an update keeps the subject ids the fit was given", {
  # As if the fit had been made with ids 101 to 130: they are recorded in the
  # id map and nowhere else. An update that rebuilt its data with the internal
  # 1:N ids replaced the map, and every subject was renamed.
  fit <- ctstantestfit
  fit$standata$idmap$original <- fit$standata$idmap$original + 100
  updated <- suppressMessages(ctFitUpdate(fit, refit = FALSE))
  expect_equal(updated$standata$idmap, fit$standata$idmap)
})

test_that("a fit saved by ctsem 3.11.1 is updated with its own arguments", {
  # 3.11.1 kept the call in `$args` itself rather than in `$args$input`, and
  # had no `backend`, `poprank` or `sampleControl`. It did have `vb`, which
  # every fit it made carries and ctFit() now refuses by name. With no
  # `backend` it is replayed on stan, which made it, not on today's default.
  v3111 <- c("stanmodeltext", "iter", "intoverstates", "binomial", "fit",
    "intoverpop", "sameInitialTimes", "stationary", "plot", "derrind",
    "optimize", "optimcontrol", "nlcontrol", "nopriors", "priors", "chains",
    "cores", "inits", "compileArgs", "forcerecompile", "saveCompile",
    "savescores", "savesubjectmatrices", "saveComplexPars", "gendata",
    "control", "verbose", "datavars")
  old <- ctstantestfit
  old$args <- c(ctstantestfit$args$input[v3111], list(vb = FALSE))
  expect_true(isTRUE(old$args$priors))

  updated <- suppressMessages(ctFitUpdate(old, refit = FALSE))
  # Replayed at their defaults, the priors this fit was estimated with would
  # be gone from the data it is evaluated against.
  expect_equal(updated$standata$priors, ctstantestfit$standata$priors)
  expect_equal(updated$standata, .stan_rebuilt()$standata)
})

# ctFit() sets a model up for one subject without random effects and with
# T0VAR fixed. An update with one subject was set up so, and read the fit's
# estimates as that smaller model's parameters: predictions from neither
# model, and no error. The subject's predictions are those it gets beside
# another.
.one_and_two <- function(fit, long, ids, ...) {
  one <- .quietly(ctFitUpdate(fit, data = long[long$id %in% ids[1], ], refit = FALSE))
  two <- .quietly(ctFitUpdate(fit, data = long[long$id %in% ids[1:2], ], refit = FALSE))
  predict <- function(f) .quietly(ctPredict(f, subjects = ids[1], ...))$value
  list(one = one, two = two, predicted = list(predict(one), predict(two)))
}

test_that("one subject updates under the fitted model, as it predicts beside another", {
  long <- ctsem:::.ctFitLongData(ctstantestfit)
  got <- .one_and_two(ctstantestfit, long, unique(long$id), timestep = 'asdata')
  expect_equal(got$one$standata$nsubjects, 1L)
  expect_equal(got$one$standata$nindvarying, ctstantestfit$standata$nindvarying)
  expect_equal(got$predicted[[1]], got$predicted[[2]])
})

test_that("an update that keeps the estimates keeps their backend", {
  expect_error(suppressMessages(ctFitUpdate(ctstantestfit, backend = "julia")),
    "refit=TRUE")
})

# Refusals that need no julia session --------------------------------------

test_that("a julia fit refuses recompile by name", {
  # `model_spec`, not the class, is what marks a julia fit (.ctFitIsJulia).
  fit <- structure(list(modelbase = ctstantestfit$ctstanmodelbase,
    model_spec = list()), class = c("ctJuliaFit", "ctFit"))
  expect_error(ctFitUpdate(fit, recompile = TRUE), "recompile")
})

test_that("a julia fit without the model it was built from says how to refit", {
  # What a julia fit made before fits carried `$modelbase` looks like
  # to this function. Its `$model` is the prepared form, which ctFit() cannot
  # prepare again, so the only way on is the model the user still has.
  fit <- structure(list(estimate = list(raw = 0), model_spec = list()),
    class = c("ctJuliaFit", "ctFit"))
  expect_error(ctFitUpdate(fit), "ctFit\\(datalong, model, inits = fit\\$estimate\\$raw")
})

# Julia -----------------------------------------------------------------------

# Ids that are not 1:N, so an update that lost them would show.
.update_data <- function(nsub = 10L, tp = 6L) {
  withr::with_seed(11, do.call(rbind, lapply(seq_len(nsub), function(i) {
    intercept <- stats::rnorm(1, 0, 0.8)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(tp)
    for (t in seq_len(tp)) {
      state <- 0.75 * state + stats::rnorm(1, 0, 0.4)
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = 100L + i, time = seq_len(tp) - 1, Y1 = y)
  })))
}

.update_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE
  model
}

# Built once and reused. `optimised` takes every default but the prior, which
# is `'randomCorr'`: the scope a replay can mistake for `priors = TRUE`, since
# the logical it is reduced to reads TRUE. `meshed` chooses its substeps per
# row.
.update_fits <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    data <- .update_data()
    model <- .update_model()
    optimised <- .quietly(ctFit(data, model, backend = "julia", cores = 1,
      priors = "randomCorr",
      optimcontrol = list(finishsamples = 20)))
    sampled <- .quietly(ctFit(data, model, backend = "julia", cores = 1,
      intoverpop = "laplace", optimize = FALSE,
      sampleControl = list(chains = 1, warmup = 25, draws = 25)))
    meshed <- .quietly(ctFit(data, model, backend = "julia", cores = 1,
      nlcontrol = list(nsubsteps = "auto"),
      optimcontrol = list(finishsamples = 20)))
    cached <<- list(optimised = optimised, sampled = sampled, meshed = meshed)
    cached
  }
})

test_that("a julia fit carries the model it was given", {
  skip_without_julia()
  fits <- .update_fits()
  for (fit in fits) {
    expect_identical(ctsem:::.ctFitBaseModel(fit), .update_model())
    # And every reader of the model the fit runs still gets that one. Stored
    # as `$ctstanmodelbase`, the stan fit's name, it was what `fit$ctstanmodel`
    # returned on a julia fit -- `$` matches a unique prefix -- and
    # prediction was re-prepared from the unprepared model.
    expect_null(fit$ctstanmodel)
    expect_identical(ctsem:::.ctFitModelObject(fit), fit$model_spec$model)
  }
  # Which is not the model it was given: handing that back to ctFit() would
  # prepare it a second time.
  expect_gt(nrow(fits$optimised$model$pars), nrow(.update_model()$pars))
})

test_that("a prepared julia model keeps the data frame its accessors read", {
  skip_without_julia()
  data <- .update_data()
  spec <- .quietly(ctFit(data, .update_model(), backend = "julia", fit = FALSE))
  # `ctFit(fit = FALSE)` returns the specification itself, whose `$data` is
  # the long frame; it used to be overwritten with the prepared list.
  expect_s3_class(spec$data, "data.frame")
  expect_equal(unique(spec$data$id), unique(data$id))
  expect_equal(nrow(spec$data), nrow(data))
  expect_identical(ctsem:::.ctFitBaseModel(spec), .update_model())
})

test_that("a julia fit updates without refitting, unchanged", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  updated <- .quietly(ctFitUpdate(fit, refit = FALSE))
  expect_s3_class(updated, "ctJuliaFit")
  expect_identical(names(updated), names(fit))
  # The specification is rebuilt from the call and compared whole. Its prior
  # is where a replayed `priors = TRUE` would show.
  expect_equal(updated$model_spec, fit$model_spec)
  expect_equal(updated$standata, fit$standata)
  expect_equal(updated$data, fit$data)
  expect_identical(updated$model, updated$model_spec$model)
  expect_identical(updated$estimate, fit$estimate)
})

test_that("new data replaces a julia fit's data and keeps its estimates", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  data <- .update_data()
  kept <- c(101L, 104L, 107L)
  updated <- .quietly(ctFitUpdate(fit, data = data[data$id %in% kept, ],
    refit = FALSE))
  # Every description of the data moves together.
  expect_equal(length(updated$model_spec$subject_starts), 3L)
  expect_equal(updated$standata$nsubjects, 3L)
  expect_equal(updated$data$nsubjects, 3L)
  expect_equal(unique(ctsem:::.ctFitLongData(updated)$id), kept)
  expect_equal(ctsem:::.ctFitIdMap(updated)$original, kept)
  expect_identical(updated$estimate, fit$estimate)
  # And the result is a fit the rest of the package can use, by the ids the
  # user knows its subjects by.
  predicted <- .quietly(ctKalman(updated, subjects = kept))
  expect_setequal(unique(as.character(predicted$Subject)), as.character(kept))
})

test_that("one subject updates a julia fit under its model, as beside another", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  data <- .update_data()
  got <- .one_and_two(fit, data, c(103L, 106L), timerange = c(0, 7), timestep = 1)
  expect_equal(ctsem:::.ctBackendNpar(got$one$model_spec), length(fit$estimate$raw))
  expect_equal(got$predicted[[1]], got$predicted[[2]])
})

test_that("ctPredict(newdata) predicts new and cut-short records under the fit", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  data <- .update_data()
  # A fitted subject's own rows as new data are that subject, predicted.
  own <- data[data$id == 104L, ]
  at <- function(k) k[k$Element %in% c('yprior', 'yupd', 'ysmooth'), 'value']
  expect_equal(at(.quietly(ctPredict(fit, newdata = own, timestep = 0.5))),
    at(.quietly(ctPredict(fit, subjects = 104L, timestep = 0.5))))
  # The same record cut short at time 2 and continued to 7, beside a new
  # subject with nothing observed: every subject in newdata by default, the
  # same predictions to time 2, and the population's for the new one.
  cut <- own[own$time <= 2, ]
  cut <- rbind(cut, transform(cut[1, ], time = 7, Y1 = NA))
  blank <- transform(cut, id = 999L, Y1 = NA)
  k <- .quietly(ctPredict(fit, newdata = rbind(cut, blank), timerange = c(0, 7),
    timestep = 0.5))
  expect_setequal(as.character(unique(k$Subject)), c('104', '999'))
  full <- .quietly(ctPredict(fit, newdata = own, timerange = c(0, 7), timestep = 0.5))
  early <- function(k) k[k$Subject == '104' & k$Element == 'yprior' & k$Time <= 2, 'value']
  expect_equal(early(k), early(full))
  population <- .quietly(ctPredict(fit, newdata = blank, timerange = c(0, 7),
    timestep = 0.5, randomEffects = 'population'))
  expect_equal(k[k$Subject == '999' & k$Element == 'ysmooth', 'value'],
    population[population$Element == 'ysmooth', 'value'])
  # And from the fit's draws, through ctKalmanArray.
  draws <- .quietly(ctKalmanArray(fit, newdata = blank, pointest = FALSE, nsamples = 5))
  expect_equal(dim(draws$yprior), c(5L, nrow(blank), 1L))
  expect_true(all(is.finite(draws$yprior)))
})

test_that("a julia refit starts from the estimate and keeps the estimator", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  refitted <- .quietly(ctFitUpdate(fit, refit = TRUE))
  expect_s3_class(refitted, "ctJuliaFit")
  expect_equal(refitted$args$input$inits, as.numeric(fit$estimate$raw))
  # Starting values turn the prior warm-up off, as they do on stan.
  expect_false(isTRUE(refitted$optim$carefulfit))
  # The prior the fit was estimated under, not one on every coordinate.
  expect_identical(refitted$args$input$priorscope, "randomCorr")
  expect_equal(refitted$model_spec$priors, fit$model_spec$priors)
  # Same data, same objective, started at its optimum: each certified within
  # the bar's standard errors of the maximum (`.ctBackendGapSe`), so the two
  # within twice that of each other. A raw tolerance of 1e-6 held only while
  # the bar was a gap of 1e-6; at 0.01 se the refit moved 2e-5.
  expect_lt(ctsem:::.ctBackendHessianDistance(fit$uncertainty$hessian,
    fit$estimate$raw, refitted$estimate$raw), 2 * ctsem:::.ctBackendGapSe)
})

test_that("a sampled julia fit is updated rather than refitted", {
  skip_without_julia()
  fit <- .update_fits()$sampled
  updated <- NULL
  suppressWarnings(suppressMessages(expect_message(
    updated <- ctFitUpdate(fit, refit = TRUE), "not refitted")))
  expect_false(is.null(updated$sample))
  expect_equal(updated$model_spec, fit$model_spec)
  expect_identical(updated$estimate, fit$estimate)
})

test_that("an automatic substep mesh goes with the rows it was chosen for", {
  skip_without_julia()
  fit <- .update_fits()$meshed
  # One substep count per row, which a specification built fresh from the
  # call does not have: the fit chose it at its estimate.
  expect_true(is.integer(fit$model_spec$max_timestep))
  expect_length(fit$model_spec$max_timestep, length(fit$model_spec$times))

  same <- .quietly(ctFitUpdate(fit, refit = FALSE))
  expect_equal(same$model_spec, fit$model_spec)

  data <- .update_data()
  fewer <- .quietly(ctFitUpdate(fit, data = data[data$id %in% c(102L, 105L), ],
    refit = FALSE))
  expect_true(is.integer(fewer$model_spec$max_timestep))
  expect_length(fewer$model_spec$max_timestep, length(fewer$model_spec$times))
})

test_that("a julia fit keeps its estimates on its own backend", {
  skip_without_julia()
  fit <- .update_fits()$optimised
  expect_error(.quietly(ctFitUpdate(fit, refit = FALSE, backend = "stan")),
    "refit=TRUE")
})
