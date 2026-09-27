# A julia fit is built in two places, and they have to agree on what a fit is.
#
# Until review/OPTIM-consolidation-plan-2026-09-25.md P5, `.ctFitJuliaBackendImpl()`
# built the object for `ctFit(optimize = TRUE)` and `.ctJuliaSampleFit()` built a
# second one for `ctFit(optimize = FALSE)`, from scratch, independently, so a
# field added to one was missing from the other and nothing said so -- the
# reader got NULL and took whatever its own default was. `$collapsedScales` was
# that: the optimised path computed it and warned, and the sampled path did
# neither.
#
# `.ctJuliaSampleFit()` now places the sampler by running
# `.ctJuliaOptimiseFit()` -- the same pipeline `ctFit(optimize = TRUE)` runs --
# and hands *that* fit to the assembler, so most of what this file used to
# catch is now structurally impossible: the two routes share one constructor
# up to `$sample` and the few `$estimate` fields a posterior mean cannot carry.
# This file is kept because "structurally impossible" is a claim worth a test,
# not an assumption, and because `ctFitUncertainty(fit, 'sample')`'s own field
# additions (`$estimate$rawposterior`, `$uncertainty`, `$sample`) are a second
# assembler `.ctBackendSampleAssemble()` that both routes now share and that
# this file does not exercise directly.
#
# (There is a third route, `ctFitUncertainty(fit, 'sample')` on an *already
# optimised* fit, which mutates that fit rather than building one and so
# cannot show a constructor divergence -- which is why the comparison below
# is between the two `ctFit()` routes and not against that one.)
#
# The allowlist below is the point -- a field that is legitimately on one
# route only has to be named here, with a reason, rather than being absent
# because someone forgot.

suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

# Deliberately tiny. This tests the *shape* of the two objects, not the
# statistics of either, so the fits only have to complete.
.shape_data <- function(nsub = 8L, tp = 5L, seed = 11L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    intercept <- stats::rnorm(1, 0, 0.8)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(tp)
    for (t in seq_len(tp)) {
      state <- 0.75 * state + stats::rnorm(1, 0, 0.4)
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(tp) - 1, Y1 = y)
  }))
}

.shape_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE
  model
}

# Both fits, built once and reused: each costs an optimisation and the sampled
# one costs a short chain on top.
.shape_fits <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    data <- .shape_data()
    model <- .shape_model()
    optimised <- suppressWarnings(suppressMessages(ctFit(data, model,
      backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
      optimcontrol = list(finishsamples = 20))))
    sampled <- suppressWarnings(suppressMessages(ctFit(data, model,
      backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
      optimize = FALSE,
      sampleControl = list(chains = 1, warmup = 25, draws = 25))))
    cached <<- list(optimised = optimised, sampled = sampled)
    cached
  }
})

# Top-level fields that are legitimately on one route only. Each needs a reason,
# and the reason is what a reader checks when this list grows.
#
# Empty now, on both sides. The sampled route places the sampler by running
# `.ctJuliaOptimiseFit()` -- the same pipeline the optimising route runs -- and
# hands *that* fit (not a second, bespoke one) to the same assembler
# `ctFitUncertainty(fit, 'sample')` uses
# (review/OPTIM-consolidation-plan-2026-09-25.md P5), so every top-level field
# the placement computed, including `$substeps`, is
# already on the sampled fit before the assembler adds `$sample`. What used to
# be here was not a reason either time it shrank: `$optim$trace`, `$laplace`
# and `$collapsedScales` were on an optimised fit and missing from a sampled
# one merely because the two constructors had drifted, and `$substeps` was the
# last of that shape -- the sampled route never re-chose a mesh, but nothing
# stopped it reporting the one the placement chose.
.SHAPE_OPTIMISED_ONLY <- character(0)

.SHAPE_SAMPLED_ONLY <- c(
  # Chain diagnostics: R-hat, ESS, divergences, step sizes, tree depths. Its
  # presence is what marks a fit as genuinely sampled, which
  # `ctOptimUncertainty()` tests in order to refuse one.
  "sample"
)

# `$estimate` sub-fields on the optimised route only.
#
# Empty now: `subject_loglik`, `states`, `innovations` and `loglik_type` are
# set by `.ctJuliaOptimiseFit()` regardless of route, so the placement fit
# already carries them (NULL-valued where they do not apply) before
# `.ctBackendSampleAssemble()` runs, the same way `$estimate$loglik` has
# always described the placement point on a sampled fit rather than the
# posterior mean that replaces it in `$raw`. `.SHAPE_ESTIMATE_SAMPLED_ONLY`
# below is the only remaining asymmetry, and it is a real one: an optimised
# fit has no comparable "point the chain started from" to report.
.SHAPE_ESTIMATE_OPTIMISED_ONLY <- character(0)

# `$optim` sub-fields on the optimised route only.
#
# Empty now, for the same reason as `.SHAPE_OPTIMISED_ONLY`: `$optim` is the
# placement run's own `optim_run` list either way, not a summary rebuilt for
# the sampled route, so every field the optimising route reports (gradient,
# Newton and stall counts, the correction stage's Hessian count, and so on)
# describes the run that placed the sampler on a sampled fit. That is a
# reportable fact about the placement, not a claim about the posterior --
# `$optim$converged` in particular is the placement's own verdict on both
# routes, per decision 5 of review/OPTIM-consolidation-plan-2026-09-25.md,
# and `$sample$converged` is the one that describes whether the chains agreed.
.SHAPE_OPTIM_OPTIMISED_ONLY <- character(0)

.SHAPE_OPTIM_SAMPLED_ONLY <- character(0)

.SHAPE_ESTIMATE_SAMPLED_ONLY <- c(
  # The Laplace point the chain started from, kept so it can be told from the
  # posterior mean that replaced it in `$raw`.
  "laplace_raw")

test_that("both ctFit routes produce a fit of the same class", {
  skip_without_julia()
  fits <- .shape_fits()
  for (f in fits) {
    expect_s3_class(f, "ctJuliaFit")
    expect_s3_class(f, "ctFit")
  }
})

test_that("the two ctFit routes agree on what a julia fit carries", {
  skip_without_julia()
  fits <- .shape_fits()
  opt <- names(fits$optimised)
  smp <- names(fits$sampled)

  missing_from_sampled <- setdiff(setdiff(opt, smp), .SHAPE_OPTIMISED_ONLY)
  missing_from_optimised <- setdiff(setdiff(smp, opt), .SHAPE_SAMPLED_ONLY)

  # Named rather than counted, so a failure says which field drifted and the
  # reader can decide whether it belongs on both or belongs in the allowlist
  # above with a reason beside it.
  expect_equal(missing_from_sampled, character(0),
    info = paste("on the optimised fit but not the sampled one:",
      paste(missing_from_sampled, collapse = ", ")))
  expect_equal(missing_from_optimised, character(0),
    info = paste("on the sampled fit but not the optimised one:",
      paste(missing_from_optimised, collapse = ", ")))
})

test_that("the fields both routes carry hold the same kind of thing", {
  skip_without_julia()
  fits <- .shape_fits()
  common <- intersect(names(fits$optimised), names(fits$sampled))
  # A field present on both but holding a list on one and a numeric on the
  # other is the same defect one level down, and reads as a data problem rather
  # than a construction one.
  for (nm in common) {
    expect_equal(class(fits$optimised[[nm]])[1], class(fits$sampled[[nm]])[1],
      info = nm)
  }
})

test_that("a sampled fit reports its identifiability findings, not just stores them", {
  skip_without_julia()
  fits <- .shape_fits()
  # Both routes compute `$identifiability`; only the optimised one calls
  # `.ctBackendIdentifyWarn()` on it (ctJuliaBackend.R). So a sampled fit knows
  # about a flat direction and never says so.
  #
  # Not asserting that a warning *fires* -- this fixture need not be
  # unidentified -- only that the finding is present in the same shape on both,
  # which is what a caller or a GUI reads.
  for (f in fits) {
    expect_true(!is.null(f$identifiability))
    expect_true(is.numeric(f$identifiability$nweak) ||
      is.integer(f$identifiability$nweak))
  }
  # The random-effect check, on both, saying where it was taken: on the
  # sampled fit the estimate is now the posterior mean, and the check
  # describes the optimum that placed the sampler.
  opt <- fits$optimised$identifiability$effects
  smp <- fits$sampled$identifiability$effects
  expect_identical(opt$point, "at the estimate")
  expect_identical(smp$point, "at the optimum that placed the sampler")
  expect_equal(smp$values, fits$sampled$estimate$laplace_raw)
  expect_identical(smp$table$effect, opt$table$effect)
})

test_that("the two routes agree on what $estimate carries", {
  skip_without_julia()
  fits <- .shape_fits()
  opt <- names(fits$optimised$estimate)
  smp <- names(fits$sampled$estimate)

  missing_from_sampled <- setdiff(setdiff(opt, smp), .SHAPE_ESTIMATE_OPTIMISED_ONLY)
  missing_from_optimised <- setdiff(setdiff(smp, opt), .SHAPE_ESTIMATE_SAMPLED_ONLY)

  expect_equal(missing_from_sampled, character(0),
    info = paste("in optimised$estimate but not sampled$estimate:",
      paste(missing_from_sampled, collapse = ", ")))
  expect_equal(missing_from_optimised, character(0),
    info = paste("in sampled$estimate but not optimised$estimate:",
      paste(missing_from_optimised, collapse = ", ")))
})

test_that("the two routes agree on what $optim carries", {
  skip_without_julia()
  fits <- .shape_fits()
  # Both routes optimise -- sampling begins by placing the sampler -- so both
  # carry a `$optim`, and a field added to one constructor and not the other is
  # the drift this file exists to catch. It caught `overshoot_parameters` once
  # already, on `$estimate`, before the two objects were separated.
  for (route in names(fits)) {
    expect_false(is.null(fits[[route]]$optim), info = route)
    expect_true(is.logical(fits[[route]]$optim$converged), info = route)
  }
  opt <- names(fits$optimised$optim)
  smp <- names(fits$sampled$optim)
  missing_from_sampled <- setdiff(setdiff(opt, smp), .SHAPE_OPTIM_OPTIMISED_ONLY)
  missing_from_optimised <- setdiff(setdiff(smp, opt), .SHAPE_OPTIM_SAMPLED_ONLY)
  expect_equal(missing_from_sampled, character(0),
    info = paste("in optimised$optim but not sampled$optim:",
      paste(missing_from_sampled, collapse = ", ")))
  expect_equal(missing_from_optimised, character(0),
    info = paste("in sampled$optim but not optimised$optim:",
      paste(missing_from_optimised, collapse = ", ")))
})

test_that("$estimate holds the estimate and nothing about the run", {
  skip_without_julia()
  fits <- .shape_fits()
  # The point of the split. Any name here that describes the search rather than
  # what it found has landed in the wrong object -- which is how `$estimate`
  # grew to thirty-five fields in the first place.
  run_shaped <- c("converged", "convergence_pending", "convergence_tolerance",
    "predicted_gain", "last_gain", "gradient", "gradient_norm", "iterations",
    "stage_iterations", "f_calls", "g_calls", "chunks", "linesearch",
    "batch_sizes", "batch_iterations", "newton_steps", "newton_hessians",
    "newton_subset_hessians", "newton_escapes", "restarts", "restarts_cancelled", "stalled", "stopped_by_gap", "overshot", "overshoot_gain",
    "overshoot_parameters", "saturated", "saturated_parameters", "carefulfit",
    "carefulfit_iterations", "carefulfit_skipped", "corrections", "hessians",
    "hessian_evaluations", "hessian_distance", "hessian", "hessian_profile",
    "trace", "substeps")
  for (route in names(fits)) {
    expect_equal(intersect(names(fits[[route]]$estimate), run_shaped),
      character(0), info = route)
  }
})

test_that("a sampled fit reports the convergence of the run that placed it", {
  skip_without_julia()
  fits <- .shape_fits()
  # Was hardcoded TRUE, which said "converged" about a run whose result was
  # sitting unread two lines away. A sampler placed from a point the optimiser
  # did not reach is worth knowing about, because the metric is built there too.
  expect_true(is.logical(fits$sampled$optim$converged))
  expect_length(fits$sampled$optim$converged, 1L)
  expect_false(is.na(fits$sampled$optim$converged))
})

test_that("the laplace block is the same shape whichever route built it", {
  skip_without_julia()
  fits <- .shape_fits()
  # Both fixtures are intoverpop='laplace', so both must carry it.
  for (route in names(fits)) {
    expect_false(is.null(fits[[route]]$laplace), info = route)
  }
  # `boundary` (a correlation capped at the optimum) and `correction` (the
  # quadrature correction stage, or its "off" status -- see
  # `.ctLaplaceAutoCorrect()`) both used to be optimised-only, because the
  # sampled route called `.ctJuliaFitLaplaceBlock()` directly and never added
  # what came after it. The placement now runs the same pipeline that does add
  # them -- `laplace_correct` forced off, so `correction$status` is `"off"`
  # rather than absent -- so the two are identical in shape.
  expect_equal(sort(names(fits$optimised$laplace)),
    sort(names(fits$sampled$laplace)))
})

test_that("estimate slots a reader depends on are present on both", {
  skip_without_julia()
  fits <- .shape_fits()
  # The fields every downstream reader treats as the fit's answer, whichever
  # route produced it.
  for (nm in c("raw", "rawposterior", "loglik")) {
    for (route in names(fits)) {
      expect_false(is.null(fits[[route]]$estimate[[nm]]),
        info = paste(route, nm))
    }
  }
  npar <- length(fits$optimised$estimate$raw)
  expect_equal(length(fits$sampled$estimate$raw), npar)
  expect_equal(ncol(fits$sampled$estimate$rawposterior), npar)
  expect_equal(ncol(fits$optimised$estimate$rawposterior), npar)
})

test_that("nothing in $uncertainty is a prefix of a sibling there", {
  skip_without_julia()
  fits <- .shape_fits()
  # R's `$` partial-matches on lists. An exact match always wins, so a prefix
  # pair is harmless while both names are present -- `$estimate$raw` beside
  # `rawposterior` has never been a problem. It becomes a hazard only where the
  # SHORTER name is optional: with `hessian` removed and `hessian_at` left,
  # `x$hessian` returned the evaluation vector instead of NULL, so every
  # `is.null(fit$uncertainty$hessian)` guard took the wrong branch and handed a
  # length-npar vector to code expecting a matrix. That is why the field is
  # `evaluated_at`.
  #
  # Checked on `$uncertainty` alone, and that narrowness is the point. A
  # blanket check over the whole fit fails on names that are fine: `model` and
  # `model_spec` at top level, and `gradient`/`gradient_norm`,
  # `saturated`/`saturated_parameters`, `carefulfit`/`carefulfit_iterations`
  # in `$optim`, none of which is ever present without its partner.
  # `$uncertainty$hessian` is the one field the package tests for absence and
  # then indexes as a matrix, so it is the one slot where the rule must hold.
  offenders <- function(nms) {
    nms <- nms[nzchar(nms)]
    out <- character(0)
    for (a in nms) {
      hit <- setdiff(nms[startsWith(nms, a)], a)
      if (length(hit)) out <- c(out, paste0(a, " <- ", paste(hit, collapse = "/")))
    }
    out
  }
  for (route in names(fits)) {
    nms <- names(fits[[route]]$uncertainty)
    skip_if(is.null(nms), "no $uncertainty on this route")
    expect_equal(offenders(nms), character(0), info = route)
  }
  # And the guard the naming protects still fires: remove the Hessian and the
  # slot reads as absent rather than as its neighbour.
  bare <- fits$optimised
  bare$uncertainty$hessian <- NULL
  expect_null(bare$uncertainty$hessian)
})
