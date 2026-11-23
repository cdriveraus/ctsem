# The maximum check: is there a higher maximum than the estimate?
#
# Everything else that judges an estimate is local. The certification, the
# flat-direction walk and the interval check read the curvature at the
# estimate and the objective a few standard errors out, and a higher maximum
# behind a dip shows in none of them. On the affect pilot's tied fit (819
# parameters, 2026-11-21, review/OPEN.md section 3) the estimate was certified
# with every direction trusted; a maximum 41.6 nats higher lay where its
# Hessian put the objective 1,899 nats lower, spread over its fifty softest
# directions, and on the straight line between the two the objective fell 69
# nats before it rose.
#
# What did show it was the joint sampler. Its chains start at the estimate and
# go where the posterior's mass is, on a density whose gradient costs far less
# than the fit's (no inner solve: 25 times less there). So the check reads the
# sampler's first learning window, which every map-sampled fit runs anyway
# (`ctsem_sample_map`'s `first_only`), in three steps:
#
#  1. Where each chain's draws are centred, as the squared distance of their
#     mean from the estimate in the wider, direction by direction, of the
#     Hessian's normal and the window's own spread within chains
#     (`ctsem_map_centre`, `_jmap_window_distances`). The mean, not the draws:
#     a posterior wider than the Hessian's normal puts every draw far out and
#     leaves the mean where it was, and in many dimensions a width a few per
#     cent off moves a draw's distance by several of its standard deviations.
#     The window's spread, because along a direction the posterior is several
#     times wider than the normal and the chains cross slowly, a window's mean
#     is near one draw: on test-julia-sample.R's random-drift fixture, centred
#     at the estimate with sds up to 4.6 times the Hessian's, the Hessian's
#     metric alone flagged 23 runs of 80, each climb finding nothing, and this
#     one none. A mean past `.ctMaximumBar()` is where neither has a
#     thousandth of its mass. No evaluation is spent on this, and a run that
#     passes goes on from the window's draws as if it had not been asked.
#  2. Only then a climb of the fit's own objective from those means
#     (`ctsem_climb_above`): L-BFGS with no Hessian, stopped once it is
#     `.ctMaximumMargin` above the estimate's value. A candidate's own value
#     decides nothing -- the pilot's chain means were 77 and 114 nats below the
#     estimate.
#  3. What a higher point means is the caller's: `ctFit(optimize = FALSE)` fits
#     again from it before sampling (`.ctJuliaSampleFit()`, `stop`), and
#     `ctFitUncertainty(fit, 'sample')`, which was handed its estimate, warns,
#     keeps the point (`$sample$maximumcheck$raw`) and samples on.
#
# A point's value rests on which of each unit's maxima its run found. A fit
# keeps its run's as its specification's known modes; a climb's stay on the
# engine objective it ran on. So an end is read by a search on that objective
# (`end_values`), `higher` is judged on that, and the maxima the search ended
# on come back as `modes`, which the fit made again takes as its first
# search's further starts (`.ctJuliaOptimiseFit()`'s `resume`, the route a
# checkpoint's modes take) and so starts on the value the check read. Without
# them the pilot's fit from a climb's end started 6.7 nats below the climb.
#
# It finds a maximum the sampler's chains reach in their first window. One
# behind a region of low mass is not found, and neither is anything where the
# first draws stay centred on the estimate.

# The bar on a chain mean's squared distance from the estimate: the 0.999
# quantile of what one draw from the Hessian's normal has.
.ctMaximumBar <- function(npar) stats::qchisq(0.999, max(1L, as.integer(npar)))

# How far above the estimate's value a climb must get: half the 95% point of a
# chi-square on one degree of freedom, the bar the flat-direction walks use
# for "the same value" (`.ctOptimFlatDirectionScreen()`). The two values are
# read with different sets of each unit's maxima, and a climb that ends this
# close to the estimate's value has shown nothing.
.ctMaximumMargin <- 1.92

# How often `ctFit(optimize = FALSE)` fits again from a higher point before it
# samples from where it is and says so.
.ctMaximumRounds <- 2L

# The fit's objective at its estimate, read with the modes its run ended with
# held (R/ctJuliaFitModes.R), which is the value the estimate maximises; the
# climb also reads the objective searched there and takes the larger.
.ctMaximumHeldValue <- function(fit, cores) {
  .ctJuliaHoldModes(fit$model_spec, .ctJuliaFitModes(fit))
  value <- tryCatch(as.numeric(.ctJuliaEvaluate(fit, gradient = FALSE,
    cores = cores)[["value"]])[1L], error = function(e) NA_real_)
  if (isTRUE(is.finite(value))) value else -Inf
}

# The condition `ctFit(optimize = FALSE)` catches to fit again from a higher
# point (`.ctJuliaSampleFit()`).
.ctMaximumCondition <- function(record) {
  structure(class = c("ctHigherMaximum", "condition"),
    list(message = "a higher maximum was found", call = NULL, record = record))
}

# The check, from the sampler's first window (`centre`, `ctsem_map_centre`'s
# record). Returns its record: the bar, each chain's distance, and, where a
# chain was past the bar, the climb -- `higher`, the estimate's value
# (`reference`), the highest value reached (`value`), `gain`, the point
# (`raw`) and the maxima a search there ended on (`modes`: `units`, one per
# column of `values`; absent where the route has none), and per start its
# value, its end's (`end_values`, searched; `climb_values` the climb's own),
# iterations and seconds. With `stop`, a higher point is signalled as a
# `ctHigherMaximum` condition carrying the record instead of warned about.
#' @keywords internal
.ctBackendMaximumCheck <- function(fit, centre, cores = 1L, verbose = 0,
  progress = .ctVerboseOn(verbose), stop = FALSE, maxiter = NULL,
  finish = isTRUE(stop)) {
  flagged <- as.logical(centre[["flagged"]])
  flagged[is.na(flagged)] <- FALSE
  record <- list(bar = as.numeric(centre[["bar"]])[1L],
    distance = as.numeric(centre[["distance"]]), flagged = flagged,
    draws = as.integer(centre[["draws"]])[1L],
    secs = as.numeric(centre[["secs"]])[1L],
    gradients = as.numeric(centre[["gradients"]])[1L],
    climbed = FALSE, higher = FALSE)
  if (!any(flagged)) return(record)

  means <- matrix(as.numeric(centre[["means"]]), ncol = length(flagged))
  npar <- nrow(means)
  estimate <- as.numeric(fit$estimate$raw)
  if (length(estimate) != npar) return(record)
  # Said to someone watching, since a climb on a large model takes a while;
  # a climb that finds nothing is otherwise only on record.
  if (isTRUE(progress)) message("The sampler's first draws are centred away ",
    "from the estimate (squared distance ",
    paste(format(round(record$distance[flagged])), collapse = ", "),
    " where under ", format(round(record$bar)),
    " is expected). Looking for a higher maximum from there.")

  objective <- .ctJuliaObjective(fit)
  module <- .ctJuliaModule(fit$model_spec$project)
  # The climb's curvature is its own, and what it leaves is no later run's.
  forget <- function() try(.ctJuliaCall("ContinuousTimeSEM.ctsem_lbfgs_forget!",
    objective, .defer = TRUE), silent = TRUE)
  forget()
  on.exit({
    try(.ctJuliaCall("ContinuousTimeSEM.ctsem_laplace_clear_anchors!",
      objective, .defer = TRUE), silent = TRUE)
    forget()
  }, add = TRUE)
  if (is.null(maxiter)) maxiter <- .ctJuliaOr(fit$args$input$optimcontrol$maxiter, 1000L)
  held <- .ctMaximumHeldValue(fit, cores)
  # The estimate's covariance by its Hessian is each climb's first metric, as
  # it is a certification's resume's: the starts are in another basin, but
  # most of the curvature is the same there (the pilot's standard errors at
  # its two maxima differed by 2-3 per cent at the quartiles).
  arguments <- list(objective, .ctJuliaPut(means[, flagged, drop = FALSE]),
    .ctJuliaNumericVector(estimate), value = held, margin = .ctMaximumMargin,
    maxiter = as.integer(maxiter)[1L], verbose = isTRUE(progress),
    progress_overwrite = .ctProgressOverwrite(verbose),
    progress_sink = if (isTRUE(progress)) .ctProgressSink(
      .ctProgressOverwrite(verbose)) else NULL, finish = isTRUE(finish))
  inverse <- .ctBackendResumeInverse(list(hessian = fit$uncertainty$hessian))
  if (!is.null(inverse) && nrow(inverse) == npar) {
    storage.mode(inverse) <- "double"
    arguments$inverse <- .ctJuliaPut(inverse)
  }
  climb <- .ctBackendWithMaxChunks(cores,
    .ctJuliaGet(do.call(module$ctsem_climb_above, arguments)))

  record$climbed <- TRUE
  record$higher <- isTRUE(climb[["higher"]])
  record$reference <- as.numeric(climb[["reference"]])[1L]
  record$value <- as.numeric(climb[["value"]])[1L]
  record$gain <- record$value - record$reference
  record$raw <- as.numeric(climb[["point"]])
  record$starts <- which(flagged)
  record$start_values <- as.numeric(climb[["start_values"]])
  record$end_values <- as.numeric(climb[["end_values"]])
  record$climb_values <- as.numeric(climb[["climb_values"]])
  units <- as.integer(climb[["mode_units"]])
  values <- as.matrix(climb[["mode_values"]])
  if (length(units) && all(units > 0L) && ncol(values) == length(units)) {
    record$modes <- list(units = units, values = values)
  }
  record$iterations <- as.integer(climb[["iterations"]])
  record$climb_secs <- as.numeric(climb[["secs"]])
  record$ends <- matrix(as.numeric(climb[["ends"]]), nrow = npar)
  record$passed <- as.integer(climb[["passed"]])
  record$passed_secs <- as.numeric(climb[["passed_secs"]])
  record$metric <- isTRUE(climb[["metric"]])
  if (!record$higher) {
    if (isTRUE(progress)) message("No higher maximum was found from there ",
      "(the highest log posterior reached is ", format(round(-record$gain, 2)),
      " below the estimate's). Sampling goes on.")
    return(record)
  }
  if (isTRUE(stop)) base::stop(.ctMaximumCondition(record))
  warning("The estimate is not the highest maximum. Optimising from where the ",
    "sampler's first draws are centred reached a log posterior ",
    format(round(record$gain, 2)), " above the estimate's after ",
    sum(record$iterations), " iterations",
    if (isTRUE(finish)) "" else ", and was stopped there", ". The fit's ",
    "estimate and standard errors describe the lower maximum. The sampled ",
    "posterior does not depend on it. ctFitUpdate(fit, inits = ",
    "fit$sample$maximumcheck$raw) fits from the higher point.", call. = FALSE)
  record
}
