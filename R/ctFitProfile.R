# Profile likelihood for a julia fit.
#
# The curvature at an estimate is a local quadratic approximation, and on a
# direction the data does not determine it is approximating something that is
# not there: the eigenvalue it reports is a residue of the transform's own
# derivative at whatever raw value the optimiser stopped at, not a property of
# the model and the data. `.ctOptimFlatDirectionScreen()` answers the cheap
# version of the question -- is the likelihood flat along this direction, right
# here -- by walking it without re-optimising. That is a slice, and it is
# deliberately one-sided: flat exhibits a curve along which the likelihood is
# constant and settles the matter, while a rise settles nothing, because the
# flat manifold may be curved and a profile would have followed it.
#
# This is the two-sided version. Each point fixes one coordinate and maximises
# over all the others, which is the profile likelihood (Raue et al. 2009), and
# the reported verdict is theirs:
#
#   * identifiable -- the profile crosses the likelihood-ratio bound on both
#     sides, and the crossings are the confidence limits;
#   * practically non-identifiable -- it falls but does not cross within the
#     range walked, so the limit is beyond it, possibly infinite;
#   * structurally non-identifiable -- one side of it does not fall at all,
#     which exhibits a curve along which the likelihood is constant and is
#     conclusive whatever the other side does.
#
# The bar is `qchisq(level, 1) / 2`, which is a statistical quantity rather
# than a tolerance on an approximation: invariant to reparameterisation and
# comparable across models, which no threshold on an eigenvalue is. The general
# reason not to threshold eigenvalues is sloppiness (Gutenkunst et al. 2007) --
# spectra are typically spread over many orders with no gap to cut at.
#
# ## A profile can find a better optimum, and then it is not a profile
#
# Every point here is a constrained optimisation, so every point is also an
# escape attempt: fix a coordinate somewhere other than where the optimiser
# left it, and the rest of the model has to accommodate it rather than absorb
# the change. If a constrained maximum comes out *above* the unconstrained one,
# the fit was not at a maximum, and everything computed from it -- this profile
# included -- describes the wrong point.
#
# So that is checked at every point and aborts the run rather than being
# reported as a curiosity. `$better` carries the point it found, for a refit.
# This is not a rare accident: it is the same mechanism the escape loop in
# `.ctJuliaOptimise()` uses, which recovered 387 and 169 log likelihood units
# on two fits started inside a flat transform.
#
# ## Cost
#
# One constrained optimisation per point per side per parameter, and each is
# roughly a fit. What makes it affordable is continuation: every point starts
# from the previous point's estimate rather than from the fit's, so a step is a
# short move rather than a fresh optimisation, and a direction stops as soon as
# it crosses the bar. It is opt-in for that reason and will never run as part
# of a fit.

#' Profile likelihood for a ctsem julia fit
#'
#' Fixes one raw parameter at a sequence of values, re-optimising all the
#' others at each, and reports what the likelihood does. This answers the
#' identification question directly rather than through the curvature at the
#' estimate, which on a flat direction measures the transform rather than the
#' data -- see the note in \code{ctOptimUncertainty} on flat directions.
#'
#' Each point is a constrained optimisation, so a profile can discover that the
#' fit was not at a maximum. When that happens the run stops and the better
#' point is returned in \code{$better}; refit from there and profile again.
#'
#' @param fit a \code{ctJuliaFit}.
#' @param parameters which parameters to profile: raw parameter names, indices,
#'   \code{NULL} for all of them, or \code{"flagged"} for only those the fit's
#'   own diagnostics already doubt -- the saturated coordinates, the
#'   unidentified directions, and the intervals wider than the curvature
#'   supports. \code{"flagged"} is the cheap internal use; naming parameters is
#'   the usual one. Profiling every parameter of a large model is a great many
#'   optimisations.
#' @param points maximum points per side.
#' @param step first displacement, in raw units. \code{NULL} takes half of each
#'   parameter's own standard error where the fit has one, so an identified
#'   parameter gets several points inside its interval instead of crossing the
#'   bar on the first step; it falls back to 0.5 raw, which is a small move in
#'   ctsem's standardised raw coordinates, when there is no curvature to ask.
#' @param growth factor the step is multiplied by after a point that moved the
#'   likelihood by less than the target, so a flat direction is walked out
#'   quickly rather than in a hundred equal steps.
#' @param level confidence level setting the likelihood-ratio bar.
#' @param maxiter iteration cap for each constrained optimisation.
#' @param verbose 0 for silence, 1 to report each point as it is computed.
#'
#' @return An object of class \code{ctFitProfile}, which
#'   \code{print} and \code{plot} both understand. \code{$profile} is one row
#'   per point (parameter, raw value, transformed value where the coordinate
#'   has one, log likelihood, drop from the estimate); \code{$summary} is one
#'   row per parameter (verdict, limits, how far each side was walked); and
#'   \code{$better} is \code{NULL} unless a constrained fit beat the estimate,
#'   in which case it holds that point and the run stopped there.
#'
#' @examples
#' \donttest{
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   # Profile only what the fit itself is unsure about.
#'   # result <- ctFitProfile(fit, parameters = "flagged")
#'   # print(result); plot(result)
#' }
#' }
#' @export
ctFitProfile <- function(fit, parameters = NULL, points = 8L, step = NULL,
  growth = 1.6, level = 0.95, maxiter = 500L, verbose = 0) {
  # `.ctFitIsJulia()` is the named predicate for this question; the class
  # literal is the same question spelled out, and `test-duplication-ratchet.R`
  # counts the spellings that are not.
  if (!.ctFitIsJulia(fit)) {
    stop("ctFitProfile() needs a fit from backend = 'julia'.", call. = FALSE)
  }
  spec <- fit$model_spec
  # Priming rather than unused: the same cached engine objective every point
  # below re-fetches (`.ctJuliaOptimise()` calls `.ctJuliaObjective()` itself,
  # and the base value below evaluates it through the fit directly), fetched
  # once here so a model shape Julia has not compiled for announces itself
  # before the loop rather than partway through it.
  invisible(.ctJuliaObjective(fit))
  # The fit's own resolved controls, the way `.ctBackendGapTolerance()` reads
  # them: a profile point used to run on the engine's bare defaults -- no
  # transform-scale metric, no batching, no Newton finish, the gap rule off --
  # so a point could land somewhere the fit itself never would, and could take
  # far more iterations getting there. `certify = FALSE` is added on top: a
  # profile point is a constrained optimum and nothing certifies it afterwards
  # (see `.ctFitProfilePoint`), so the optimiser's cheap stopping rule should
  # not assume a Hessian will close the rest of the gap the way it may for the
  # fit itself.
  optimcontrol <- utils::modifyList(.ctFitProfileOptimcontrol(fit),
    list(certify = FALSE))
  gradient <- .ctJuliaOr(optimcontrol$gradient, "adjoint")
  cores <- suppressWarnings(as.integer(fit$args$resolved$cores)[1L])
  if (!length(cores) || is.na(cores) || cores < 1L) cores <- 1L
  estimate <- as.numeric(fit$estimate$raw)
  npar <- length(estimate)
  names <- .ctBackendRawParameterNames(fit, npar)
  bar <- stats::qchisq(level, 1) / 2

  index <- .ctFitProfileParameters(parameters, names, npar, fit)
  if (!length(index)) {
    stop("ctFitProfile(): nothing to profile. With parameters = 'flagged' ",
      "that means the fit's own diagnostics found nothing to doubt.",
      call. = FALSE)
  }
  # One step size per parameter, because the parameters are not on one scale
  # even in raw coordinates. Half a standard error puts a handful of points
  # inside the interval for an identified parameter, which is what makes the
  # reported limit an interpolation between neighbours rather than a single
  # linear guess from the estimate to the first step -- measured on a fitted
  # correlation, a flat 0.5 crossed the bar on point one and reported the limit
  # off one interpolation.
  steps <- .ctFitProfileSteps(step, fit, npar)
  # The base value is evaluated directly at the estimate (`ctJuliaEvaluate()`),
  # not by an optimisation capped at one iteration as this used to do: nothing
  # is pinned, so nothing needs optimising, and an evaluate cannot take a step
  # a real fit would not have. `$loglik` is the likelihood and the objective
  # may be the posterior (`$logposterior`); a profile whose drops are taken
  # against a different quantity from the one being maximised would report
  # crossings that are not crossings.
  base <- .ctFitProfilePoint(fit, spec, optimcontrol, gradient, cores,
    estimate, integer(), numeric(), maxiter)
  if (!is.finite(base$value)) {
    stop("ctFitProfile(): the objective is not finite at the estimate.",
      call. = FALSE)
  }

  rows <- list()
  points_raw <- list()
  better <- NULL
  for (k in index) {
    if (!is.null(better)) break
    for (side in c(-1, 1)) {
      # Continuation: each point starts from the last one's estimate, so a
      # step is a short move rather than a fresh optimisation. The first step
      # of each side starts from the fit.
      from <- estimate
      at <- estimate[k]
      size <- steps[k]
      for (point in seq_len(as.integer(points))) {
        at <- at + side * size
        start <- from
        start[k] <- at
        got <- .ctFitProfilePoint(fit, spec, optimcontrol, gradient, cores,
          start, k, at, maxiter)
        drop <- base$value - got$value
        rows[[length(rows) + 1L]] <- data.frame(
          parameter = names[k], index = k, side = side, value = at,
          loglik = got$value, drop = drop, iterations = got$iterations,
          stringsAsFactors = FALSE)
        points_raw[[length(points_raw) + 1L]] <- got$par
        if (verbose > 0) {
          message(sprintf("%s = %+.4f  loglik %.5f  drop %.5f", names[k], at,
            got$value, drop))
        }
        if (!is.finite(got$value)) break
        # A constrained maximum above the unconstrained one. The fit was not at
        # a maximum, so this profile describes the wrong point and the only
        # honest thing to do is stop and say where the better one is.
        if (drop < -.ctFitProfileTolerance()) {
          better <- list(point = got$par, gain = -drop,
            parameter = names[k], value = at)
          break
        }
        from <- got$par
        # Grow the step while the likelihood is barely moving, so a flat
        # direction is walked out rather than crawled along. The target is a
        # fraction of the bar per point, which is what makes the ladder adapt
        # to the model rather than to whoever chose `step`.
        #
        # A step that changed *nothing* is the case the multiplier cannot
        # rescue on its own. The starting step is half a standard error, and
        # the coordinates whose profile matters most are exactly those whose
        # standard error is near zero -- so the ladder would start at 1e-6 and
        # need forty doublings to reach anywhere. Measured: a diffusion
        # correlation on a flat ray walked 0.003 raw units over six points,
        # against 35 with a flat starting step, and reported a verdict on
        # almost no evidence. A step that produced no measurable change says
        # nothing about what the step should be, so fall back to the model-free
        # default and grow from there.
        if (drop < bar / as.integer(points)) {
          size <- if (drop < bar * 1e-6) max(size * growth, .ctFitProfileStep())
            else size * growth
        }
        if (drop >= bar) break
      }
      if (!is.null(better)) break
    }
  }

  profile <- if (length(rows)) do.call(rbind, rows) else
    data.frame(parameter = character(), index = integer(), side = numeric(),
      value = numeric(), loglik = numeric(), drop = numeric(),
      iterations = integer(), stringsAsFactors = FALSE)
  points_raw <- if (length(points_raw)) do.call(rbind, points_raw) else NULL
  profile$transformed <- .ctFitProfileTransformed(fit, profile, points_raw)
  out <- list(profile = profile,
    summary = .ctFitProfileSummary(profile, names, index, estimate, bar),
    bar = bar, level = level, base = base$value, better = better,
    estimate = stats::setNames(estimate, names), call = match.call())
  class(out) <- "ctFitProfile"
  if (!is.null(better)) {
    warning("ctFitProfile(): a constrained fit beat the estimate by ",
      signif(better$gain, 3), " while profiling ", better$parameter,
      ", so this fit was not at a maximum and the profile stopped. Refit from ",
      "$better$point.", call. = FALSE)
  }
  out
}

# Each profiled point on the scale a reader reports, where the coordinate has
# one.
#
# A raw value is what was profiled and is unambiguous, but it is not what
# anybody writes down. Every point here is a whole parameter vector, so the
# model quantities at it can be *evaluated* rather than back-solved -- which
# also sidesteps the trap that a raw coordinate can feed more than one cell:
# nothing is inverted and nothing is reported as "the" value of anything.
#
# `NA` for a coordinate that is not a fixed-effect cell. Population standard
# deviations and correlations are built from the random-effect block rather
# than read off it, so `.ctBackendPopMeanSamples()` has no column for them and
# inventing one here would mean duplicating that construction. Their raw scale
# is what this reports, and `$summary` says so.
#' @keywords internal
.ctFitProfileTransformed <- function(fit, profile, points_raw) {
  out <- rep(NA_real_, nrow(profile))
  if (!nrow(profile) || is.null(points_raw)) return(out)
  mapped <- try(.ctBackendPopMeanSamples(fit, samples = points_raw),
    silent = TRUE)
  if (inherits(mapped, "try-error") || is.null(mapped$values)) return(out)
  values <- as.matrix(mapped$values)
  if (nrow(values) != nrow(profile)) return(out)
  for (k in unique(profile$index)) {
    column <- match(k, as.integer(mapped$parnumber))
    if (is.na(column) || column > ncol(values)) next
    rows <- which(profile$index == k)
    out[rows] <- as.numeric(values[rows, column])
  }
  out
}

# The first displacement for each parameter.
#
# Half a standard error where the fit has curvature to ask, and 0.5 raw where
# it does not -- an `estonly` fit has no `se`, and neither has a coordinate the
# uncertainty stage projected out, which is exactly a coordinate whose profile
# is the interesting one. A non-finite or zero se falls back the same way.
#' @keywords internal
.ctFitProfileSteps <- function(step, fit, npar,
  default = .ctFitProfileStep()) {
  if (!is.null(step)) {
    value <- suppressWarnings(as.numeric(step))
    if (length(value) == 1L) value <- rep(value, npar)
    if (length(value) != npar || any(!is.finite(value)) || any(value <= 0)) {
      stop("ctFitProfile(): step must be one positive number, or one per ",
        "parameter.", call. = FALSE)
    }
    return(value)
  }
  se <- suppressWarnings(as.numeric(fit$estimate$se))
  if (length(se) != npar) return(rep(default, npar))
  ifelse(is.finite(se) & se > 0, se / 2, default)
}

# Which coordinates to profile.
#
# `"flagged"` is the internal use and the cheap one: profile what this fit
# already doubts rather than all of it. Three separate detectors contribute and
# they do not agree by construction -- a saturated transform, a flat direction
# in the curvature, and an interval wider than that curvature supports are
# three different complaints -- so the union is taken. That is the right
# direction to err in here: a parameter profiled needlessly costs time, one
# skipped costs the answer.
#' @keywords internal
.ctFitProfileParameters <- function(parameters, names, npar, fit = NULL) {
  # A dispatch rather than a ladder of early returns, so the four ways of
  # naming parameters sit beside each other and none of them reads as the
  # special case.
  if (is.null(parameters)) return(seq_len(npar))
  if (identical(parameters, "flagged")) {
    flagged <- .ctFitProfileFlagged(fit, names, npar)
    return(flagged)
  }
  if (is.character(parameters)) {
    named <- match(parameters, names)
    if (anyNA(named)) {
      stop("ctFitProfile(): no parameter named ",
        paste(sQuote(parameters[is.na(named)]), collapse = ", "), ".",
        call. = FALSE)
    }
    return(named)
  }
  index <- suppressWarnings(as.integer(parameters))
  if (anyNA(index) || any(index < 1L) || any(index > npar)) {
    stop("ctFitProfile(): parameter indices must be between 1 and ", npar, ".",
      call. = FALSE)
  }
  index
}

# What the fit's own diagnostics doubt, as raw coordinate indices.
#
# Each source may hold names or indices depending on how far through the
# reporting it came, so both are accepted rather than one being assumed. `0` is
# the engine's "none" sentinel -- a zero-length vector deadlocks the R bridge --
# and is dropped here.
#' @keywords internal
.ctFitProfileFlagged <- function(fit, names, npar) {
  resolve <- function(x) {
    if (is.null(x) || !length(x)) return(integer())
    if (is.character(x)) return(match(x, names))
    index <- suppressWarnings(as.integer(x))
    index
  }
  found <- c(
    resolve(fit$identifiability$parameters),
    resolve(fit$uncertainty$intervalcheck$parameters),
    resolve(fit$uncertainty$intervalcheck$unidentified),
    resolve(fit$optim$saturated_parameters))
  found <- found[!is.na(found) & found >= 1L & found <= npar]
  sort(unique(found))
}

# A rise smaller than this is arithmetic, not a better optimum.
#' @keywords internal
.ctFitProfileTolerance <- function() 1e-6

# The step to use when the curvature cannot suggest one. Raw units, where
# ctsem's coordinates are standardised by construction -- priors are
# normal(0, 1) and each transform carries its own scale -- so half a unit is a
# small move for any model and needs to know nothing about this one.
#' @keywords internal
.ctFitProfileStep <- function() 0.5

# The fit's own optimcontrol, read the way `.ctBackendGapTolerance()` reads it
# for the same reason: `ctFit()` stores it under `$args$resolved` and
# `$args$input`, and the backend's own `$args` -- what a fit carries while it
# is being built, and what a stored fit from before that split has -- keeps it
# at the top. Reading only the top found nothing once `ctFit()` had returned,
# which is `ctFitProfile()`'s only caller.
#' @keywords internal
.ctFitProfileOptimcontrol <- function(fit) {
  args <- fit$args
  optimcontrol <- args$resolved$optimcontrol
  if (is.null(optimcontrol)) optimcontrol <- args$input$optimcontrol
  if (is.null(optimcontrol)) optimcontrol <- args$optimcontrol
  if (is.null(optimcontrol)) list() else as.list(optimcontrol)
}

# One constrained optimum: `index` pinned at `value`, everything else free --
# or, when `index` is empty, the objective evaluated at `start` with nothing
# pinned, which is how the base value is taken (see `ctFitProfile()`); the
# same evaluation `.ctJuliaOptimise()`'s pinned points are optimised against,
# so the base and the points cannot come from two different quantities.
#
# Routed through `.ctJuliaOptimise()` with the fit's own resolved
# `optimcontrol`, rather than a bare `module$ctsem_optimize()` call on engine
# defaults: a profile point used to be optimised with none of the fit's
# controls -- no transform-scale metric, no batching, no Newton finish, the
# gap rule off -- so a point could land somewhere, and take however long
# getting there, that the fit itself never would. The pin is
# `.ctJuliaOptimise()`'s own pin path (`optimise_once(from, damp, pin)`),
# exposed as its `pin` argument, which keeps the metric, the batching and the
# stall escapes in play around the fixed coordinate exactly as they are for
# the fit. Profile points are never certified afterwards (`ctFitProfile()`
# sets `optimcontrol$certify = FALSE`), and `$better` is what stands in for
# certification here: a constrained point that beats the unconstrained
# estimate is the one failure a profile can still catch.
#' @keywords internal
.ctFitProfilePoint <- function(fit, spec, optimcontrol, gradient, cores,
  start, index, value, maxiter) {
  if (!length(index)) {
    out <- try(ctJuliaEvaluate(fit, pars = as.numeric(start), gradient = FALSE,
      gradient_method = gradient), silent = TRUE)
    if (inherits(out, "try-error") || !length(out$value)) {
      return(list(value = NA_real_, par = start, iterations = NA_integer_))
    }
    return(list(value = as.numeric(out$value)[1L], par = start,
      iterations = 0L))
  }
  out <- try(.ctJuliaOptimise(spec, as.numeric(start),
    optimcontrol = optimcontrol, gradient = gradient, cores = cores,
    maxiter = as.integer(maxiter), verbose = 0L,
    pin = list(index = as.integer(index), value = as.numeric(value))),
    silent = TRUE)
  if (inherits(out, "try-error")) {
    return(list(value = NA_real_, par = start, iterations = NA_integer_))
  }
  par <- as.numeric(out$minimizer)
  # Restored rather than trusted, as `.ctJuliaOptimise()`'s own pin path
  # already does internally -- defensive here too, since a caller must be able
  # to rely on the pinned coordinate reading back exactly.
  par[index] <- as.numeric(value)
  list(value = as.numeric(out$maximum_loglik)[1L], par = par,
    iterations = as.integer(.ctJuliaOr(out$iterations, NA_integer_)))
}

# The verdict per parameter, and the limits where there are any.
#
# Raue et al.'s three outcomes, decided by what the walk actually did rather
# than by a threshold on curvature. A side that crossed the bar has a limit,
# found by interpolating between the last two points -- linear in the log
# likelihood, which is exact for the quadratic the bar assumes and is the usual
# reading of a profile.
#
# `structurally non-identifiable` needs only *one* flat side, not both, and
# that is the same one-sided logic `.ctOptimFlatDirectionScreen()` uses: a side
# along which the likelihood does not move has exhibited a curve along which it
# is constant, and that settles the matter whatever the other side does. The
# case is not hypothetical -- a diffusion correlation on noise data walks flat
# to raw -44 and falls on the other side, because the transform saturates in
# one direction and not the other. Calling that merely practical would
# understate it.
#
# A limit of `NA` means the walk did not find one, which is not the same as
# there being none: `lower_walked` and `upper_walked` say how far it looked, so
# a reader can tell "unbounded as far as raw -44" from "we took two steps".
#' @keywords internal
.ctFitProfileSummary <- function(profile, names, index, estimate, bar) {
  empty <- data.frame(parameter = character(), index = integer(),
    estimate = numeric(), lower = numeric(), upper = numeric(),
    flat = character(), lower_walked = numeric(), upper_walked = numeric(),
    verdict = character(), stringsAsFactors = FALSE)
  if (!nrow(profile)) return(empty)
  rows <- lapply(index, function(k) {
    mine <- profile[profile$index == k, , drop = FALSE]
    if (!nrow(mine)) return(NULL)
    limit <- function(side) {
      part <- mine[mine$side == side, , drop = FALSE]
      part <- part[is.finite(part$drop), , drop = FALSE]
      if (!nrow(part)) return(list(at = NA_real_, crossed = FALSE, moved = FALSE))
      crossed <- which(part$drop >= bar)
      moved <- max(part$drop) > .ctFitProfileTolerance()
      if (!length(crossed)) {
        return(list(at = NA_real_, crossed = FALSE, moved = moved))
      }
      first <- crossed[1L]
      # The point before the crossing, or the estimate itself when the very
      # first step already crossed.
      previous <- if (first > 1L) part[first - 1L, ] else
        data.frame(value = estimate[k], drop = 0)
      span <- part$drop[first] - previous$drop
      at <- if (!is.finite(span) || span <= 0) part$value[first] else
        previous$value + (part$value[first] - previous$value) *
          (bar - previous$drop) / span
      list(at = at, crossed = TRUE, moved = TRUE)
    }
    low <- limit(-1)
    high <- limit(1)
    verdict <- if (low$crossed && high$crossed) "identifiable" else
      if (!low$moved || !high$moved) "structurally non-identifiable" else
        "practically non-identifiable"
    walked <- function(side) {
      part <- mine[mine$side == side & is.finite(mine$drop), , drop = FALSE]
      if (!nrow(part)) return(NA_real_)
      part$value[which.max(abs(part$value - estimate[k]))]
    }
    data.frame(parameter = names[k], index = k, estimate = estimate[k],
      lower = low$at, upper = high$at,
      flat = paste(c(if (!low$moved) "lower", if (!high$moved) "upper"),
        collapse = "+"),
      lower_walked = walked(-1), upper_walked = walked(1),
      verdict = verdict, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

#' Print a profile likelihood
#'
#' @param x a \code{ctFitProfile}.
#' @param digits significant digits for the reported values.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.ctFitProfile <- function(x, digits = 3, ...) {
  # The better optimum first, and before anything else, because every other
  # number below it describes a point that is not the maximum. A reader who
  # stops after one line should stop after the right one.
  if (!is.null(x$better)) {
    cat("This fit was NOT at a maximum.\n")
    cat(sprintf("  Profiling %s at %s found a point %s log units better.\n",
      x$better$parameter, signif(x$better$value, digits),
      signif(x$better$gain, digits)))
    cat("  Refit from $better$point; the profile below stopped there and\n")
    cat("  describes the wrong point.\n\n")
  }
  cat(sprintf("Profile likelihood, %g%% (%s log units), %d point%s\n",
    100 * x$level, signif(x$bar, digits), nrow(x$profile),
    if (nrow(x$profile) == 1L) "" else "s"))
  if (!nrow(x$summary)) {
    cat("Nothing profiled.\n")
    return(invisible(x))
  }
  table <- x$summary[, c("parameter", "estimate", "lower", "upper",
    "verdict"), drop = FALSE]
  for (column in c("estimate", "lower", "upper")) {
    table[[column]] <- signif(table[[column]], digits)
  }
  # An absent limit is not a missing number, it is a statement -- the walk went
  # this far and did not find one -- so it is printed as where the walk got to
  # rather than as NA.
  edge <- function(limit, walked) ifelse(is.na(limit),
    paste0("<", signif(walked, digits)), format(limit))
  table$lower <- edge(x$summary$lower, x$summary$lower_walked)
  table$upper <- edge(x$summary$upper, x$summary$upper_walked)
  print(table, row.names = FALSE)
  cat("\nLimits are on the raw scale, which is what was profiled.")
  if (any(is.na(x$summary$lower) | is.na(x$summary$upper))) {
    cat("\n`<value` means the profile was walked to there without crossing",
      "\nthe bar, so the limit is beyond it and may be unbounded.")
  }
  cat("\n")
  invisible(x)
}

#' Plot a profile likelihood
#'
#' One panel per parameter: the drop in log likelihood against the value the
#' parameter was held at, with the likelihood-ratio bar drawn across. A curve
#' that crosses the bar on both sides has a confidence interval; one that runs
#' along the bottom is a direction the data does not determine.
#'
#' @param x a \code{ctFitProfile}.
#' @param parameters which to draw; defaults to all that were profiled.
#' @param ... passed to \code{plot}.
#' @return \code{x}, invisibly.
#' @export
plot.ctFitProfile <- function(x, parameters = NULL, ...) {
  profile <- x$profile
  if (!nrow(profile)) {
    message("Nothing to plot: no profile points.")
    return(invisible(x))
  }
  which <- if (is.null(parameters)) unique(profile$parameter) else
    intersect(parameters, unique(profile$parameter))
  if (!length(which)) {
    message("None of those parameters were profiled.")
    return(invisible(x))
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  rows <- ceiling(sqrt(length(which)))
  graphics::par(mfrow = c(rows, ceiling(length(which) / rows)),
    mar = c(4, 4, 2, 1))
  for (name in which) {
    part <- profile[profile$parameter == name, , drop = FALSE]
    part <- part[is.finite(part$value) & is.finite(part$drop), , drop = FALSE]
    estimate <- unname(x$estimate[name])
    # The estimate is a profile point too -- a drop of zero by definition --
    # and including it is what makes the curve meet the axis rather than
    # starting a step away from it.
    part <- rbind(part[0, ], data.frame(parameter = name, index = NA_integer_,
      side = 0, value = estimate, loglik = x$base, drop = 0,
      iterations = NA_integer_, transformed = NA_real_,
      stringsAsFactors = FALSE), part)
    part <- part[order(part$value), , drop = FALSE]
    graphics::plot(part$value, part$drop, type = "b", pch = 16,
      xlab = paste(name, "(raw)"), ylab = "drop in log likelihood",
      main = name, ylim = range(c(0, part$drop, x$bar * 1.1), finite = TRUE),
      ...)
    graphics::abline(h = x$bar, lty = 2)
    graphics::abline(v = estimate, lty = 3)
    limits <- x$summary[x$summary$parameter == name, , drop = FALSE]
    if (nrow(limits)) {
      for (at in c(limits$lower, limits$upper)) {
        if (is.finite(at)) graphics::abline(v = at, col = "grey40")
      }
    }
  }
  invisible(x)
}
