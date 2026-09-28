# Certifying convergence, rather than asserting it from a gradient ------------
#
# A gradient is not a convergence criterion. `max|g|` has units of log
# likelihood per unit parameter, so what counts as small depends on how each
# parameter happens to be scaled, and the engine's own rule compared it against
# `1e-6 * max(1, |loglik|)` -- a threshold in log likelihood units, which is
# dimensionally a different thing, carries the likelihood's arbitrary additive
# constants, and therefore *loosens* when the manifest variables are rescaled.
#
# What can be certified, once the curvature is known, is the quantity the
# optimizer is actually trying to make zero:
#
#     gap = 1/2 g' H^-1 g
#
# with `H` the observed information. Under the local quadratic approximation
# that is the log likelihood still available -- the predicted improvement from
# taking the Newton step -- so it is in log likelihood units, it is invariant
# to any nonsingular linear reparameterisation, and a tolerance set on it means
# the same thing for every model. `lambda = sqrt(2 * gap)` is the Mahalanobis
# length of that Newton displacement in the information metric: in one
# dimension it is literally the displacement in standard errors, and in several
# it is the joint distance, which is not the same claim.
#
# ## Flat directions are checked, not discarded
#
# `H^-1` does not exist when a direction has no curvature, and the obvious
# repair -- compute the gap on the positive-curvature subspace and ignore the
# rest -- is unsafe. The gradient can have a large component in a near-null
# eigenvector, and dropping it would report a tiny gap at a point that is
# plainly not stationary: the flatter the direction, the more of the remaining
# likelihood it can hold, and the more certain the projected answer looks.
#
# So the excluded directions are *measured* instead. `gap` covers the subspace
# whose curvature the data supports; the component of the gradient outside it
# is stepped along and the actual change in log likelihood recorded. Measuring
# rather than thresholding matters because a norm in those directions has the
# same units problem this whole file exists to remove, and because in a flat
# direction the quadratic model is exactly what cannot be trusted -- the
# engine's own `_ctsem_overshot` answers its question the same way, by taking
# the step and looking.
#
# ## Six outcomes, not two
#
# In the order `.ctBackendCertify()` checks them, and the first that applies
# wins -- except that a gap which is not a number is also `unknown`, found
# after `saturated` because it is only read there.
#
#   unknown        the curvature could not be decomposed, or the gap is not a
#                  number, so nothing is claimed either way
#   notmaximum     the overshoot probe found better, or a direction has genuine
#                  negative curvature: the point is a saddle whatever the gap
#                  says
#   notstationary  a flat direction with a live gradient: stepping along it
#                  gains likelihood, so this is not a maximum, and the gap
#                  cannot see it because that direction was excluded
#   saturated      a transform has saturated, where the gradient underflows to
#                  zero and the curvature with it, so no tolerance means
#                  anything for that coordinate. The point is still a maximum
#   suboptimal     the gap says the optimum is measurably above this estimate
#   certified      the gap is below tolerance and the excluded directions hold
#                  no material likelihood
#
# `saturated` was called `unidentified` until 2026-09-25. It never meant more
# than a saturated transform -- a flat direction whose probe gains nothing is
# `certified`, and `fit$identifiability` is what names it -- so the old word
# claimed a verdict this status does not make. A stored fit may still carry it;
# `.ctBackendCertificationStatus()` reads it as the new one.
#
# `notstationary` and `saturated` were one status, and they are not one
# finding: the first is a fit nobody should use and the second is a fit with a
# result in it -- a population scale with no individual differences behind it is
# the usual cause, and reporting that as a failure to converge is the mistake
# `test-julia-convergence.R` exists to prevent.
#
# The saturation verdict is consulted rather than recomputed. A saturated
# parameter has both a zero gradient and no curvature, so it lands in the
# excluded subspace with nothing to report, and a certification that did not
# ask would contradict the warning the fit already carries.
#
# ## Arithmetic in the engine, verdicts here
#
# The Hessian, the Newton steps, the escape from a saddle and the probe along
# the flat directions are the engine's: its finish (`_ctsem_newton_finish`)
# ends the optimiser's run on them and hands them back, so a fit is certified
# without another engine call in the common case. What stays here is what is
# cheap on a matrix and tested without an engine -- the gap, the split into
# trusted and flat directions, the verdict and its words -- and the one
# decision about what to do next (`.ctBackendCorrectResult()`).

# Eigen-decomposition of the observed information, with the trusted subspace
# marked.
#
# `rtol` is `.ctBackendNullMass()`'s, deliberately: two rules for "this
# direction has no curvature" that can disagree is how a fit comes to be
# described one way by its intervals and another by its convergence.
#' @keywords internal
.ctBackendInformationSplit <- function(hessian, rtol = .ctFlatDirectionRtol(),
  negative = 1e-8) {
  if (is.null(hessian)) return(NULL)
  hessian <- as.matrix(hessian)
  if (nrow(hessian) != ncol(hessian) || !nrow(hessian)) return(NULL)
  if (!all(is.finite(hessian))) return(NULL)
  information <- -(hessian + t(hessian)) / 2
  decomposition <- try(eigen(information, symmetric = TRUE), silent = TRUE)
  if (inherits(decomposition, "try-error")) return(NULL)
  values <- decomposition$values
  scale <- max(values)
  if (!is.finite(scale) || scale <= 0) {
    # Not a maximum in any direction, so there is no trusted subspace at all
    # and the caller needs to hear that rather than an empty gap.
    scale <- max(abs(values))
    if (!is.finite(scale) || scale <= 0) return(NULL)
    return(list(values = values, vectors = decomposition$vectors,
      scale = scale, trusted = rep(FALSE, length(values)),
      negative = values < -negative * scale))
  }
  list(values = values, vectors = decomposition$vectors, scale = scale,
    trusted = values > rtol * scale,
    negative = values < -negative * scale)
}

# The predicted remaining log likelihood, and what is outside it.
#
# Returns the gap over the trusted subspace, the Newton displacement that would
# realise it, and the gradient component the projection left behind -- as a
# vector, because the caller measures it rather than comparing its norm to
# anything.
#' @keywords internal
.ctBackendOptimGap <- function(hessian, gradient, rtol = .ctFlatDirectionRtol(),
  negative = 1e-8) {
  gradient <- as.numeric(gradient)
  split <- .ctBackendInformationSplit(hessian, rtol = rtol,
    negative = negative)
  empty <- list(gap = NA_real_, lambda = NA_real_, step = NULL,
    residual = NULL, residual_norm = NA_real_, ntrusted = NA_integer_,
    nflat = NA_integer_, nnegative = NA_integer_, ok = FALSE)
  if (is.null(split) || length(gradient) != length(split$values)) return(empty)
  if (!all(is.finite(gradient))) return(empty)

  keep <- split$trusted
  vectors <- split$vectors[, keep, drop = FALSE]
  values <- split$values[keep]
  # Coefficients of the gradient in the trusted eigenbasis. The gap is a sum of
  # squares over curvature, which is where the invariance comes from: rescale a
  # parameter and both the coefficient and the eigenvalue move to compensate.
  coefficients <- as.numeric(crossprod(vectors, gradient))
  gap <- if (!length(values)) 0 else 0.5 * sum(coefficients^2 / values)
  step <- if (!length(values)) rep(0, length(gradient)) else
    as.numeric(vectors %*% (coefficients / values))
  residual <- gradient - as.numeric(vectors %*% coefficients)
  list(gap = gap, lambda = sqrt(2 * max(gap, 0)), step = step,
    residual = residual, residual_norm = sqrt(sum(residual^2)),
    ntrusted = sum(keep), nflat = sum(!keep & !split$negative),
    nnegative = sum(split$negative),
    # The direction the point is not a maximum along -- the most negative
    # curvature -- so a report can say which parameters it runs through.
    negative_vector = if (any(split$negative))
      split$vectors[, which.min(split$values)] else NULL,
    # The smallest curvature still trusted, which is what sets how tight a
    # gradient has to be before the gap can be under a given tolerance.
    lambda_min = if (length(values)) min(values) else NA_real_, ok = TRUE)
}

# What the excluded directions are actually worth, in log likelihood.
#
# A norm of the leftover gradient cannot answer this: it has the units the rest
# of this file exists to avoid, and in a direction with no curvature the
# quadratic model that would convert it into a likelihood is precisely the one
# that does not hold. So step along it and look. The stepping is the engine's
# (`_ctsem_flat_probe`, at 0.25, 1 and 4 raw units), through the route's own
# predicate, so a Laplace point whose inner solve did not converge is refused
# rather than compared -- which the R-side probe this replaced could not see.
#
# A fit's certification reads the probe the optimiser's finish already ran
# (`.ctBackendEndgameOf()`); this is for a certification assembled from a
# Hessian R holds (`.ctBackendCertification()`). NULL when there was nothing to
# probe or the engine could not, which `.ctBackendCertify()` reads as no gain
# measured.
#
# Besides the gain it returns what a reader needs to weigh it: the length that
# gave it, the longest length the objective could be evaluated at, and the unit
# direction, so a message can name the parameters it runs through. A gain of
# 5.5e-05 at a quarter of a unit, found by a probe that looked no further than
# four, says something different from a ridge, and the verdict alone could not
# tell them apart.
#' @keywords internal
.ctBackendFlatProbe <- function(x, at, direction) {
  direction <- as.numeric(direction)
  size <- sqrt(sum(direction^2))
  if (!is.finite(size) || size <= 0) return(NULL)
  module <- .ctJuliaModule(.ctBackendSpec(x)$project)
  out <- try(.ctJuliaGet(module$ctsem_flat_probe(
    .ctJuliaObjective(x), .ctJuliaNumericVector(as.numeric(at)),
    .ctJuliaNumericVector(direction))), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  .ctBackendProbeFields(out, length(direction))
}

# The probe's numbers off an engine result, in the shape `.ctBackendCertify()`
# reads, or NULL when it did not run. The engine never sends an empty vector --
# it would hang the bridge -- so `probe_ran` is what says there is a direction.
#' @keywords internal
.ctBackendProbeFields <- function(r, npar) {
  if (!isTRUE(r$probe_ran)) return(NULL)
  direction <- suppressWarnings(as.numeric(r$probe_direction))
  if (length(direction) != npar || !all(is.finite(direction))) return(NULL)
  list(gain = as.numeric(r$probe_gain)[1L],
    length = as.numeric(r$probe_length)[1L],
    longest = as.numeric(r$probe_longest)[1L], direction = direction,
    ok = TRUE)
}

# The parameters a direction runs through: those carrying at least `share` of
# its largest component, largest first, at most `most` of them.
#
# One rule for the three places that describe a direction by its parameters --
# the not-a-maximum and flat-gain messages below, and the identifiability
# report -- so a direction is described the same way wherever it is named. A
# share of the largest rather than an absolute loading: a unit vector spread
# over ten coordinates has loadings near 0.32 and one over twenty near 0.22, so
# an absolute bar of 0.25 names every coordinate of the first and none of the
# second, and on a ridge whose loadings shift as the optimiser walks it the
# same bar named ten parameters at one stopping point and nine at another.
#' @keywords internal
.ctBackendLoadedCoordinates <- function(vector, share = 1 / 3, most = Inf) {
  size <- abs(as.numeric(vector))
  if (!length(size) || !any(is.finite(size)) || max(size, na.rm = TRUE) <= 0) {
    return(integer())
  }
  keep <- which(is.finite(size) & size >= share * max(size, na.rm = TRUE))
  keep <- keep[order(size[keep], decreasing = TRUE)]
  keep[seq_len(min(length(keep), most))]
}

# What a flat direction that still gains is worth, in the words a reader needs:
# which parameters, how much, how far along, and how far the probe looked.
#
# The verdict itself stays `notstationary` -- a small gain that turns within the
# probe did not mean a maximum on AnomAuth S2, which had 1.6 exact nats further
# out -- so what changes is that the message says what was measured rather than
# asserting more than that. `parnames` is optional; without it the direction is
# described by its size alone.
#' @keywords internal
.ctBackendFlatGainReason <- function(probe, parnames = NULL) {
  gain <- suppressWarnings(as.numeric(probe$gain)[1L])
  at <- suppressWarnings(as.numeric(probe$length)[1L])
  longest <- suppressWarnings(as.numeric(probe$longest)[1L])
  involved <- character()
  if (!is.null(parnames) && length(probe$direction) == length(parnames)) {
    involved <- parnames[.ctBackendLoadedCoordinates(probe$direction, most = 4L)]
  }
  paste0("a flat direction",
    if (length(involved)) paste0(" (", paste(involved, collapse = ", "), ")"),
    " still gains ", signif(gain, 2),
    if (isTRUE(at > 0)) paste0(" within ", signif(at, 3),
      " raw units of the estimate"),
    if (isTRUE(longest > 0)) paste0("; the probe looked no further than ",
      signif(longest, 3)))
}

# The verdict, in the terms a reader needs.
#
# `saturated` is the fit's own, not recomputed here: a saturated transform
# reports a zero gradient and no curvature, so it passes every test in this
# file for the wrong reason, and the fit already says so.
#
# `parnames`, when given, lets a `notstationary` reason name the parameters the
# flat direction runs through. Everything else here is arithmetic on the gap
# and the probe, which is what keeps it testable without a fit.
#' @keywords internal
.ctBackendCertify <- function(gap, probe = NULL, tolerance = 0.01,
  saturated = FALSE, overshot = FALSE, parnames = NULL) {
  if (is.null(gap) || !isTRUE(gap$ok)) {
    return(list(status = "unknown", certified = FALSE,
      reason = "the curvature at the estimate could not be decomposed"))
  }
  residual_gain <- if (is.null(probe)) NA_real_ else as.numeric(probe$gain)[1L]
  if (isTRUE(overshot)) {
    return(list(status = "notmaximum", certified = FALSE,
      reason = "the optimizer overstepped into a flat region of a transform"))
  }
  if (isTRUE(gap$nnegative > 0L)) {
    return(list(status = "notmaximum", certified = FALSE,
      reason = paste0(gap$nnegative, " direction",
        if (gap$nnegative > 1L) "s have" else " has",
        " negative curvature, so this point is not a maximum")))
  }
  material <- is.finite(residual_gain) && residual_gain > tolerance
  if (material) {
    return(list(status = "notstationary", certified = FALSE,
      reason = .ctBackendFlatGainReason(probe, parnames)))
  }
  if (isTRUE(saturated) && isTRUE(gap$nflat > 0L)) {
    return(list(status = "saturated", certified = FALSE,
      reason = paste0("a parameter transform has saturated, where the ",
        "gradient underflows to zero and the curvature with it, so no ",
        "tolerance here means anything for that coordinate. Nothing here says ",
        "the estimate is not a maximum -- see fit$identifiability for which ",
        "coordinate the data does not determine")))
  }
  if (!is.finite(gap$gap)) {
    return(list(status = "unknown", certified = FALSE,
      reason = "the predicted gap is not a number"))
  }
  if (gap$gap > tolerance) {
    return(list(status = "suboptimal", certified = FALSE,
      reason = paste0("the optimum is about ", signif(gap$gap, 3),
        " log likelihood above this estimate")))
  }
  list(status = "certified", certified = TRUE,
    reason = paste0("the optimum is within ", signif(tolerance, 3),
      " log likelihood of this estimate"))
}

# The certification's number, in the shape `ctLaplaceCheck()` and
# `ctParticleLik()` also report -- see R/ctFitGap.R for why one shape.
#' @keywords internal
.ctBackendCertifyGap <- function(gap, tolerance = 0.01) {
  .ctFitGap("curvature", gap = if (is.null(gap)) NA_real_ else gap$gap,
    tolerance = tolerance,
    remedy = paste0("A gap above the bar means the optimiser stopped short; ",
      "ctFit() continues from the curvature automatically unless ",
      "optimcontrol$certify = FALSE."))
}

# The bar a fit has to clear, in objective units.
#
# Two candidate principles, and the binding one is not the obvious one.
#
# For the *estimate*, `lambda = sqrt(2 * gap)` is the displacement in the
# information metric -- standard errors -- so numerical error is negligible
# against statistical error once `lambda` is around 0.1: a tenth of a standard
# error is about one percent of the variance and moves no inference. That would
# put the bar at 0.005.
#
# For the *curvature-based reports* it is nowhere near enough. Identifiability
# and the interval check classify eigenvalues against a relative threshold, and
# an eigenvalue near that threshold is exquisitely sensitive to where the
# estimate sits: on a deliberately degenerate fixture -- six subjects, four
# waves, ten population correlations -- the flat direction and the five
# correlations it names go unreported at a gap of 3e-04 and are found at 4e-07.
# A ridge that goes unreported is the expensive kind of wrong answer here,
# because the fit looks fine and the numbers look plausible.
#
# So the diagnostics set the bar, at 1e-06, and the estimate gets far better
# precision than it needs as a side effect. `lambda` goes out beside the gap so
# a reader can have either reading.
#
# The controls are read from wherever the fit keeps them. `ctFit()` stores them
# under `$args$resolved` and `$args$input`; the backend's own `$args` -- what a
# fit carries while it is being built, and what a stored fit from before that
# split has -- keeps them at the top. Reading only the top found nothing once
# `ctFit()` had returned, so `ctOptimUncertainty(fit)` re-certified every fit
# at the default, whatever `gaptol` it was fitted with.
#' @keywords internal
.ctBackendGapTolerance <- function(fit, default = 1e-6) {
  args <- fit$args
  optimcontrol <- args$resolved$optimcontrol
  if (is.null(optimcontrol)) optimcontrol <- args$input$optimcontrol
  if (is.null(optimcontrol)) optimcontrol <- args$optimcontrol
  .ctBackendConvergeTol(optimcontrol, default = default)
}

# The bar itself, from the controls rather than from a fit.
#
# `.ctBackendGapTolerance()` needs a fit and the optimiser needs the same number
# before there is one, so the reading lives here and both ask it. One default in
# one place: a second copy of `1e-6` is how a tolerance gets changed in one
# reader and not the other.
#' @keywords internal
.ctBackendConvergeTol <- function(optimcontrol = list(), default = 1e-6) {
  value <- if (is.null(optimcontrol)) NULL else optimcontrol$gaptol
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(value) || value <= 0) return(default)
  value
}

# Which directions the post-fit overshoot probe pulls back, from the controls.
#
# `"magnitude"` (the default) pulls back prefixes of the coordinates ordered by
# |raw|, jointly; `"saturation"` pulls back each coordinate the saturation
# detector flagged, one at a time, which is what this did before; `"off"` skips
# it. See `_ctsem_overshot` in the engine for what each costs and catches.
#
# Off is a real option and not a footgun to be hidden: the probe costs up to
# `8 * npar` value-only evaluations per optimisation stage (one per pullback
# fraction per prefix), and on a large model that is a visible fraction of the
# fit -- 3% to 38% on the optimiser bench. What it buys is that a fit stopped in
# a degenerate corner says so instead of reporting convergence. The engine
# skips it where every coordinate is small and nothing is flagged
# (`_CTSEM_OVERSHOOT_MIN_RAW`), which is where it never found anything.
#' @keywords internal
.ctBackendOvershootProbe <- function(optimcontrol = list(),
    default = "magnitude") {
  value <- if (is.null(optimcontrol)) NULL else optimcontrol$overshoot
  if (is.null(value)) return(default)
  value <- as.character(value)[1L]
  allowed <- c("magnitude", "saturation", "off")
  if (!value %in% allowed) {
    stop("optimcontrol$overshoot must be one of ",
      paste(allowed, collapse = ", "), ", not ", sQuote(value), call. = FALSE)
  }
  value
}

# When the engine looks for a stalled fit, and how readily.
#
# The window is a number of iterations and the fraction is a share of the
# progress the fit has already made -- not a number of nats, which are not
# comparable across models, and not a share of what is predicted to remain,
# which is the one quantity that goes wrong exactly here. See `_ctsem_stalled`
# in the engine.
#
# This is only half of the test. Firing it costs a derivative pass over the
# transforms, and the fit is stopped only if that finds one flat -- so the
# fraction is deliberately loose, and it is the conjunction rather than this
# bar that decides anything. See `_ctsem_stall_verdict!`.
#
# `stallwindow = 0` switches the whole check off.
#
# Provenance: the window of 80 iterations and the fraction of 1e-2 were set by
# hand on the flat-transform fixtures the watch was built for (fcae6459,
# 2026-09-14) -- a drift started inside `-log1p_exp`'s flat region, now the
# slow tier of `test-julia-convergence.R`, and `test_state_sampling.jl`'s count
# model for the case the conjunction must leave alone. No sweep over either is
# on record, nor which other values were tried, and neither has been checked on
# categorical or multilevel models (plan of 2026-09-25, Appendix B and 2i). The
# watch's other four constants -- cooldown 30, tightening 0.1 at most twice,
# flat ratio 1e-3 -- have the same provenance and are read inline in
# `.ctJuliaOptimise()`.
#' @keywords internal
.ctBackendStallWindow <- function(optimcontrol = list(), default = 80L) {
  value <- if (is.null(optimcontrol)) NULL else optimcontrol$stallwindow
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.integer(value)[1L])
  if (is.na(value) || value < 0L) default else value
}

# The fraction's provenance is the window's, above: set by hand on the
# flat-transform fixtures, not swept.
#' @keywords internal
.ctBackendStallFraction <- function(optimcontrol = list(), default = 1e-2) {
  value <- if (is.null(optimcontrol)) NULL else optimcontrol$stallfraction
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(value) || value < 0) default else value
}

# Where to restart a fit that stopped on a boundary, or NULL if there is
# nowhere better to go.
#
# The engine only stops a stage when its own probe has already found a better
# point -- stalled, flat, *and* somewhere to go, all three -- so that point
# comes back on the result and is used as it stands. Recomputing the ladder
# here would cost `8 * npar` objective evaluations for an answer already in
# hand, which on a laplace fit near a degenerate corner is a couple of minutes.
#
# The fallback is for a fit that finished normally and is nonetheless sitting
# somewhere it should not be: zero the coordinates whose transforms have gone
# flat and refit. That cannot be decided cheaply -- the zeroed point is
# deliberately *worse* on the spot, 9.5 nats down on the fit it was measured
# on, and a refit from it landed 7.2 nats up -- so it costs a fit and is off by
# default. `optimcontrol$escapesaturated = TRUE` turns it on.
#
# Zeroing is not neutral and is not claimed to be: raw zero is a correlation of
# 0 but a drift of -0.693 and an sd of 0.693. It is a starting value for a
# refit, not an answer, and the refit is kept only if it wins.
#' @keywords internal
.ctBackendStallEscape <- function(result, optimcontrol, model_spec,
    verbose = 0, escapes = TRUE) {
  # `escapes = FALSE` on the state-explicit route, where the joint mode is
  # degenerate: the innovations re-optimise to absorb almost any parameter
  # change, so a pullback can nearly always find something and "not a maximum"
  # stops being informative. The same reason the two stopping rules are off
  # there -- see `.ctJuliaOptimise()`.
  if (!isTRUE(escapes)) return(NULL)
  npar <- length(as.numeric(result$minimizer))
  point <- suppressWarnings(as.numeric(result$stall_point))
  if (isTRUE(result$stopped_by_stall) && length(point) == npar &&
      all(is.finite(point))) {
    if (verbose > 0) {
      message("Stopped on a boundary; pulling raw parameter(s) ",
        paste(.ctJuliaSaturatedNames(list(p = result$stall_parameters),
          model_spec, npar, "p"), collapse = ", "), " back gains ",
        signif(as.numeric(result$stall_gain)[1L], 3), ", refitting from there")
    }
    return(.ctBackendEscapeCoordinates(point, result$stall_parameters, npar))
  }
  # A stage that ended *without* stalling and without converging, whose own
  # post-fit probe found a better point. That probe runs at the end of every
  # fit and has already paid for the point, so resuming from it costs nothing
  # extra -- and `overshot` is one of the three things that make `converged`
  # false, so this only ever fires on a fit that has already failed.
  #
  # The in-flight check cannot cover this case and is not meant to. It asks
  # whether the last `stallwindow` iterations got anywhere, so a fit that
  # climbs steadily and then stops dead -- a line search that finds nothing,
  # which is how most of these end -- never presents a stalled window at all.
  # Measured on a laplace fit started at drift raw 12, inside its own flat
  # transform: 188 iterations, not converged, 169 nats short, and the progress
  # test never fired once.
  point <- suppressWarnings(as.numeric(result$overshoot_point))
  if (isTRUE(result$overshot) && length(point) == npar &&
      all(is.finite(point))) {
    if (verbose > 0) {
      message("Stopped somewhere that is not a maximum; pulling raw ",
        "parameter(s) ", paste(.ctJuliaSaturatedNames(
          list(p = result$overshoot_parameters), model_spec, npar, "p"),
          collapse = ", "), " back gains ",
        signif(as.numeric(result$overshoot_gain)[1L], 3),
        ", refitting from there")
    }
    return(.ctBackendEscapeCoordinates(point,
      result$overshoot_parameters, npar))
  }
  if (!isTRUE(optimcontrol$escapesaturated)) return(NULL)
  flat <- suppressWarnings(as.integer(result$saturated_parameters))
  flat <- flat[!is.na(flat) & flat >= 1L & flat <= npar]
  if (!length(flat)) return(NULL)
  from <- as.numeric(result$minimizer)
  from[flat] <- 0
  if (verbose > 0) {
    message("Saturated with no pullback available; zeroing raw parameter(s) ",
      paste(.ctJuliaSaturatedNames(list(p = result$saturated_parameters),
        model_spec, npar, "p"), collapse = ", "),
      " and refitting from there to see whether it beats this")
  }
  .ctBackendEscapeCoordinates(from, flat, npar)
}

# Which coordinates the escape moved, carried on the point it returns.
#
# An attribute rather than a second return value because every caller of
# `.ctBackendStallEscape()` wants the point and only one wants this, and
# because `NULL` for "no point" is the contract that decides whether to escape
# at all. `0` is the engine's "none" sentinel -- a zero-length vector deadlocks
# the R bridge -- so it is dropped here rather than reaching the pin as a
# coordinate index.
#' @keywords internal
.ctBackendEscapeCoordinates <- function(point, coordinates, npar) {
  index <- suppressWarnings(as.integer(coordinates))
  index <- index[!is.na(index) & index >= 1L & index <= npar]
  attr(point, "coordinates") <- unique(index)
  point
}

# The certification for one fit, at its own estimate, against a Hessian R holds.
#
# Used by the uncertainty stage when it could not reuse the certification the
# fit was made with -- a Hessian from somewhere else, or a fit certified with
# none (`estonly`). The probe is what makes the flat directions safe to exclude
# from the gap: their gradient is stepped along and the objective measured,
# rather than a norm being compared to a threshold that would carry the units
# this whole approach exists to remove.
#' @keywords internal
.ctBackendCertification <- function(fit, hessian, tolerance = 0.01) {
  gradient <- as.numeric(fit$optim$gradient)
  gap <- .ctBackendOptimGap(hessian, gradient)
  probe <- NULL
  if (isTRUE(gap$ok) && isTRUE(gap$residual_norm > 0)) {
    probe <- .ctBackendFlatProbe(fit, as.numeric(fit$estimate$raw),
      gap$residual)
  }
  parnames <- try(.ctBackendRawParameterNames(fit, length(gradient)),
    silent = TRUE)
  if (inherits(parnames, "try-error")) parnames <- NULL
  .ctBackendCertificationRecord(gap, probe, tolerance = tolerance,
    saturated = isTRUE(fit$optim$saturated),
    overshot = isTRUE(fit$optim$overshot), parnames = parnames)
}

# The certification list a fit carries, from the gap and the probe.
#
# One assembly for both places a certification is made -- the correction loop,
# and the uncertainty stage when it cannot reuse that one -- so a field added
# here reaches whichever made the certification a reader ends up with. There
# were two assemblies until 2026-09-25, and a field added to one never reached
# a reader of the other: the first attempt at `$verdict` got a NULL that way.
#' @keywords internal
.ctBackendCertificationRecord <- function(gap, probe, tolerance = 0.01,
  saturated = FALSE, overshot = FALSE, parnames = NULL) {
  verdict <- .ctBackendCertify(gap, probe, tolerance = tolerance,
    saturated = saturated, overshot = overshot, parnames = parnames)
  out <- list(status = verdict$status, certified = verdict$certified,
    reason = verdict$reason, tolerance = tolerance,
    # The same `$verdict` shape `ctLaplaceCheck()` and `ctParticleLik()` report;
    # see R/ctFitGap.R.
    verdict = .ctBackendCertifyGap(gap, tolerance),
    gap = gap$gap, lambda = gap$lambda)
  c(out, .ctBackendProbeRecord(probe, parnames),
    list(ntrusted = gap$ntrusted, nflat = gap$nflat, nnegative = gap$nnegative,
      lambda_min = gap$lambda_min, negative_vector = gap$negative_vector,
      # The displacement the gap predicts, kept so a caller can try it without
      # decomposing the Hessian again.
      step = gap$step))
}

# What the flat-direction probe measured, as certification fields.
#
# One assembly for both places a certification is built, so a field added here
# reaches whichever of them wrote the certification a reader ends up with --
# the note on `.ctBackendCertification()` above is about exactly that failure.
# `residual_parameters` are the coordinates the probed direction runs through,
# by `.ctBackendLoadedCoordinates()`'s rule, and `residual_direction` is the
# unit vector itself, so a caller can profile along it.
#' @keywords internal
.ctBackendProbeRecord <- function(probe, parnames = NULL) {
  if (is.null(probe)) {
    return(list(residual_gain = 0, residual_length = 0, residual_longest = 0,
      residual_direction = NULL, residual_parameters = character()))
  }
  direction <- probe$direction
  involved <- if (!is.null(parnames) && length(direction) == length(parnames))
    parnames[.ctBackendLoadedCoordinates(direction, most = 4L)] else character()
  list(residual_gain = probe$gain, residual_length = probe$length,
    residual_longest = .ctJuliaOr(probe$longest, 0),
    residual_direction = direction, residual_parameters = involved)
}

# The verdict on the fit, once the curvature has been measured.
#
# `converged` used to be the engine's own: a gradient against a bar, computed
# before anything knew the curvature. The certification is the same question
# answered properly -- how much objective is still available, in objective
# units, invariantly to reparameterisation -- and it was being computed, stored
# and reported while `converged` went on saying what the gradient bar thought.
# A fit could be certified and report `converged = FALSE`, which is the one
# reading a user takes at face value.
#
# So the better measurement wins, in both directions. A certified fit is
# converged. A fit whose curvature says the optimum is above it is not, however
# small its gradient -- which is the case a gradient bar cannot see at all,
# since a flat direction reports no gradient and no curvature.
#
# Applied wherever a certification is attached: after the correction loop in
# `.ctJuliaOptimise()`, and after `ctOptimUncertainty()` computes one for a fit
# that had none. Where nothing certified -- `estonly`, `certify = FALSE` -- the
# engine's verdict is all there is, and `convergence_pending` says so.
#' @keywords internal
.ctBackendCertifiedVerdict <- function(fit) {
  certification <- fit$uncertainty$certification
  if (is.null(certification) || !length(certification$status)) return(fit)
  # `converged` answers "is this a maximum", which is what a reader takes it
  # for. `certified` is the stronger claim that also bounds how far the optimum
  # can be, and a saturated coordinate defeats that bound without saying
  # anything against the maximum -- so the two statuses that are findings
  # rather than failures map to TRUE.
  fit$optim$converged <- .ctBackendCertificationStatus(certification) %in%
    c("certified", "saturated")
  # Superseded rather than answered: the optimizer's complaint was held for
  # this measurement, and the measurement has now been made. Leaving it would
  # make `.ctBackendCertifyWarn()` warn about a gradient on a fit whose
  # curvature already said better.
  fit$optim$convergence_pending <- NULL
  fit
}

# The certification status, with the spelling a stored fit may carry read as
# the current one.
#
# `unidentified` became `saturated` on 2026-09-25, because a saturated
# transform is all it ever meant (see the header of this file). A fit stored
# before then still says `unidentified`, so every reader of the status goes
# through here rather than comparing against a literal, and the two spellings
# cannot be read two ways.
#' @keywords internal
.ctBackendCertificationStatus <- function(certification) {
  status <- if (is.list(certification)) certification$status else NULL
  if (!length(status)) return(character())
  status <- as.character(status)[1L]
  if (identical(status, "unidentified")) "saturated" else status
}

# The one convergence statement a fit makes.
#
# Five cases, and they are genuinely different things to tell a user:
#
#   saturated        the coordinate the data does not determine, once,
#                    pointing at the identifiability report rather than
#                    repeating it here.
#   certified        nothing. The optimizer's own stopping rule was superseded
#                    by a criterion it does not know about, and repeating its
#                    complaint would send someone chasing a fit that is right.
#   notstationary    what the flat-direction probe measured -- the parameters,
#                    the gain, how far along it and how far the probe looked --
#                    and what would settle it. The verdict is right and stays:
#                    a small gain that turns within the probe did not mean a
#                    maximum on AnomAuth S2. What it cannot do is tell a
#                    5.5e-05 bump from a ridge, so it says what it saw.
#   not certified    what the curvature says, in objective units, because that
#                    is actionable where a gradient is not.
#   no certification the gradient, and the fact that nothing checked further.
#' @keywords internal
.ctBackendCertifyWarn <- function(fit) {
  certification <- fit$uncertainty$certification
  pending <- isTRUE(fit$optim$convergence_pending)
  if (is.null(certification) || !length(certification$status)) {
    if (!pending) return(invisible(NULL))
    warning("The optimizer stopped without meeting its convergence criterion: ",
      "largest gradient ", signif(as.numeric(fit$optim$gradient_norm), 3),
      ". No curvature was computed, so how far this is from the optimum is ",
      "unknown -- fit without optimcontrol$estonly, or call ",
      "ctOptimUncertainty(), to have it certified.", call. = FALSE)
    return(invisible(NULL))
  }
  if (isTRUE(certification$certified)) return(invisible(NULL))
  status <- .ctBackendCertificationStatus(certification)
  # A maximum with a coordinate the data does not determine is a finding, and
  # the finding is `fit$identifiability`'s to report. Warning "not converged"
  # here is what sent 45 of 64 good fits back to be re-run.
  if (identical(status, "saturated")) {
    warning("This fit is a maximum, but ", certification$reason, ".",
      call. = FALSE)
    return(invisible(NULL))
  }
  if (identical(status, "notstationary")) {
    involved <- as.character(certification$residual_parameters)
    warning("Not converged: ", certification$reason, ". More iterations, ",
      "other starts, or ctFitProfile() ",
      if (length(involved)) paste0("on ", involved[1L]) else "along it",
      " would say whether it keeps rising.", call. = FALSE)
    return(invisible(NULL))
  }
  if (identical(status, "notmaximum") &&
      length(certification$negative_vector)) {
    warning(.ctBackendNotMaximumMessage(fit, certification), call. = FALSE)
    return(invisible(NULL))
  }
  warning("This fit is not certified as converged: ", certification$reason,
    ". See fit$uncertainty$certification.", call. = FALSE)
  invisible(NULL)
}

# What a not-a-maximum verdict tells a user: which parameters the likelihood
# still rises along, and what they are doing there -- read off the ascent
# direction (the most negative curvature, signed by the gradient) and the
# current raw values, rather than a stock remedy. The parameters are those
# carrying at least a third of the direction's largest component, at most four.
# Two readings are specific enough to say out loud, because each points at the
# model rather than the optimiser: a random-effect sd whose raw value is
# falling (the data may not support that random effect), and a random-effect
# correlation moving away from zero (heading for +-1, where two effects act as
# one -- the one case a lower poprank describes).
#' @keywords internal
.ctBackendNotMaximumMessage <- function(fit, certification) {
  v <- as.numeric(certification$negative_vector)
  npar <- length(v)
  names <- .ctBackendRawParameterNames(fit, npar)
  gradient <- as.numeric(fit$optim$gradient)
  raw <- as.numeric(fit$estimate$raw)
  if (length(gradient) >= npar && sum(gradient[seq_len(npar)] * v) < 0) v <- -v
  keep <- .ctBackendLoadedCoordinates(v, most = 4L)
  involved <- names[keep]
  readings <- character()
  for (k in seq_along(keep)) {
    i <- keep[k]; name <- involved[k]
    if (startsWith(name, "popsd_") && v[i] < 0) {
      readings <- c(readings, paste0("the sd of ", sub("^popsd_", "", name),
        " is shrinking towards zero, so the data may not support that random ",
        "effect"))
    } else if (startsWith(name, "rawcor_") && length(raw) >= i &&
        is.finite(raw[i]) && sign(v[i]) == sign(raw[i]) && raw[i] != 0) {
      pair <- strsplit(sub("^rawcor_", "", name), "__", fixed = TRUE)[[1L]]
      readings <- c(readings, paste0("the correlation between ",
        paste(pair, collapse = " and "), " is heading towards ",
        if (raw[i] > 0) "+1" else "-1", ", where the two random effects act as ",
        "one (a lower poprank)"))
    }
  }
  restarts <- fit$optim$restarts
  tried <- if (is.data.frame(restarts) && nrow(restarts))
    paste0(" ", nrow(restarts), " random restart",
      if (nrow(restarts) > 1L) "s" else "", " found nothing better.") else ""
  paste0("This fit is not a maximum: the likelihood still rises along ",
    paste(involved, collapse = ", "), ".",
    if (length(readings)) paste0(" ", toupper(substr(readings[1L], 1L, 1L)),
      substring(paste(readings, collapse = "; "), 2L), ".") else
      " Those parameters may not be separately determined by the data.",
    tried, " See fit$uncertainty$certification.")
}

# Which engine function gives the curvature of this spec's objective.
#
# `ctsem_hessian_forward` exists for one case, a sampled TI predictor value on
# the augmented route, and has a method for the plain objective only. A Laplace
# spec builds a `CTSEMLaplaceObjective`, whose only Hessian is `ctsem_hessian`
# -- and on that route `gradient` does not choose how anything is
# differentiated, since the Laplace gradient is a forward sweep over the
# reverse pass either way (see `.ctJuliaOptimise()`). So `gradient = 'forward'`
# on a Laplace fit sent the certification to a method that does not exist:
# the call raised a MethodError inside a `try`, the Hessian came back NULL,
# and the correction loop stopped without a word -- on its first round
# whenever the optimiser's finish had handed back no Hessian, so with no
# verdict to continue on. Measured on `test-julia-laplace.R`'s 30-subject
# fixture with `maxiter = 3`: no Hessian at all, where `'adjoint'` computed
# one. (What the loop then does is the verdict's: in one session `'adjoint'`
# continued to the optimum 60 nats higher, in a fresh one it stopped on
# `notstationary` -- inferred, not verified, to be the inner modes the
# Laplace Hessian warm-starts from.) The Laplace route
# therefore takes `ctsem_hessian` whatever `gradient` says, and so does the
# sampler's metric on the same spec (`.ctBackendHessian()`).
#' @keywords internal
.ctBackendHessianFunction <- function(spec, gradient = "adjoint") {
  if (!is.null(spec$laplace)) return("ctsem_hessian")
  if (identical(gradient, "forward")) "ctsem_hessian_forward" else "ctsem_hessian"
}

# The curvature at one point, from the spec rather than from a fit.
#
# `.ctBackendHessian()` is the same call with a fit around it; this is for a
# caller holding a spec and no fit. The certification no longer is one: its
# Hessian comes from the engine's finish, or from `ctsem_endgame` where no
# finish ran (`.ctBackendEndgameAt()`).
#' @keywords internal
.ctBackendHessianAt <- function(spec, est, gradient = "adjoint") {
  # Classed here because the spec reaches this point unclassed: `ctFit()`
  # carries `model_spec` as a plain list and `.ctJuliaOptimise()` does the same
  # thing for the same reason.
  spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  name <- .ctBackendHessianFunction(spec, gradient)
  available <- isTRUE(tryCatch(is.function(module[[name]]),
    error = function(e) FALSE))
  if (!available) return(NULL)
  est <- as.numeric(est)
  out <- try(.ctBackendJuliaValue(module[[name]](.ctJuliaObjective(spec),
    .ctJuliaVector(est))), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  matrix(as.numeric(out), nrow = length(est), ncol = length(est))
}

# How far one point is from another, in the standard errors a Hessian implies.
#
# The largest `|to[i] - from[i]| / se[i]`, with `se` the square root of the
# diagonal of the inverse information over its trusted directions -- the
# covariance a fit reports. A coordinate the trusted curvature says nothing
# about has a standard error of zero there, so moving it at all is infinitely
# far, which is the safe side for the one question this answers: whether a
# Hessian evaluated at `from` still counts as the curvature at `to`. The engine
# asks the same question with the same arithmetic (`_ctsem_hessian_distance`)
# when its finish decides whether to keep the Hessian it took at the hand-over.
#' @keywords internal
.ctBackendHessianDistance <- function(hessian, from, to,
  rtol = .ctFlatDirectionRtol()) {
  d <- as.numeric(to) - as.numeric(from)
  if (!length(d) || all(d == 0)) return(0)
  split <- .ctBackendInformationSplit(hessian, rtol = rtol)
  if (is.null(split) || length(split$values) != length(d)) return(Inf)
  vectors <- split$vectors[, split$trusted, drop = FALSE]
  values <- split$values[split$trusted]
  variance <- if (length(values)) rowSums(sweep(vectors^2, 2L, values, "/")) else
    rep(0, length(d))
  ratio <- ifelse(d == 0, 0, ifelse(variance > 0, abs(d) / sqrt(variance), Inf))
  max(ratio)
}

# How far, in those standard errors, a Hessian may have been evaluated from a
# point and still be the curvature there: a hundredth. The engine's finish
# keeps the Hessian it took at the hand-over when its steps moved the estimate
# less than this, and everything that reuses a stored Hessian asks the same
# question (`.ctBackendStoredHessian()`), so one number serves both sides and
# the engine is passed this one. Chosen, not measured: decision 3 of
# review/OPTIM-consolidation-plan-2026-09-25.md, for the bench to test.
#' @keywords internal
.ctBackendHessianReuse <- function() 0.01

# The fit's stored Hessian when it describes the curvature at `est`, else NULL.
#
# It does when it was evaluated at `est` exactly, or within
# `.ctBackendHessianReuse()` standard errors of it by its own inverse: the
# engine's finish keeps a Hessian taken at the hand-over as the final one when
# its steps moved the estimate less than that, and says where it was taken in
# `$uncertainty$evaluated_at`. Reading `evaluated_at` rather than assuming
# `$estimate$raw` is what stops a sampled fit's Hessian, taken at the Laplace
# estimate, being reused at its posterior mean -- a different point, and far
# more than a hundredth of a standard error away.
#' @keywords internal
.ctBackendStoredHessian <- function(fit, est) {
  stored <- fit$uncertainty$hessian
  at <- fit$uncertainty$evaluated_at
  est <- as.numeric(est)
  if (is.null(stored) || !is.matrix(stored) || is.null(at) ||
      nrow(stored) != length(est) || ncol(stored) != length(est) ||
      length(at) != length(est)) {
    return(NULL)
  }
  if (isTRUE(all.equal(as.numeric(at), est, tolerance = 0))) return(stored)
  moved <- .ctBackendHessianDistance(stored, at, est)
  if (isTRUE(moved <= .ctBackendHessianReuse())) stored else NULL
}

# The curvature's numbers from an optimiser result: the Hessian the engine's
# finish ended on (`_ctsem_newton_finish`), where it was evaluated, how far the
# estimate is from there in the standard errors it implies, and what the probe
# along the untrusted directions found. NULL when the result carries no usable
# Hessian -- a stage that ended on its iteration cap or a stall, where no finish
# ran, or a finish that could not form one.
#' @keywords internal
.ctBackendEndgameOf <- function(result, npar) {
  hessian <- result$hessian
  if (!is.matrix(hessian) || nrow(hessian) != npar || ncol(hessian) != npar ||
      !all(is.finite(hessian))) {
    return(NULL)
  }
  at <- suppressWarnings(as.numeric(result$hessian_evaluated_at))
  if (length(at) != npar || !all(is.finite(at))) {
    at <- as.numeric(result$minimizer)[seq_len(npar)]
  }
  distance <- suppressWarnings(as.numeric(result$hessian_distance)[1L])
  list(hessian = hessian, evaluated_at = at,
    distance = if (length(distance) && is.finite(distance)) distance else 0,
    probe = .ctBackendProbeFields(result, npar), hessians = 0L)
}

# The same numbers at a point the optimiser left without a finish: the engine's
# `ctsem_endgame`, which forms the exact Hessian there and runs the probe, and
# takes no step. Through the chunk wrapper like every engine run, with the
# ceiling left as it is, so an instrument that counts engine runs counts this.
#' @keywords internal
.ctBackendEndgameAt <- function(spec, est, gradient = "adjoint", verbose = 0) {
  spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  available <- isTRUE(tryCatch(is.function(module$ctsem_endgame),
    error = function(e) FALSE))
  if (!available) return(NULL)
  est <- as.numeric(est)
  npar <- length(est)
  # A certification that forms its own Hessian is exactly the case a finish
  # never ran for -- the optimiser stopped on its iteration cap or a stall --
  # and on the Laplace route that Hessian is `2 * npar` gradients, the same
  # slow loop as the finish's own (see `ctsem_laplace_hessian`). Reported
  # through the same line and sink as everything else, at the same default
  # verbosity, so it does not go quiet exactly where the finish's fix does not
  # reach: this path bypasses the finish entirely.
  reporting <- .ctBackendReporting(verbose)
  out <- try(.ctBackendWithMaxChunks(NA_integer_, .ctJuliaGet(
    module$ctsem_endgame(.ctJuliaObjective(spec), .ctJuliaNumericVector(est),
      gradient_method = gradient, flat_rtol = .ctFlatDirectionRtol(),
      progress = reporting, progress_overwrite = .ctProgressOverwrite(verbose),
      progress_sink = .ctBackendProgressSink(verbose), progress_label = "certify"))),
    silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  hessian <- out$hessian
  if (!is.matrix(hessian) || nrow(hessian) != npar || ncol(hessian) != npar ||
      !all(is.finite(hessian))) {
    return(NULL)
  }
  list(hessian = hessian, evaluated_at = est, distance = 0,
    probe = .ctBackendProbeFields(out, npar),
    hessians = as.integer(.ctJuliaOr(out$newton_hessians, 1L)))
}

# Certify the optimiser's result, and resume the optimiser when the curvature
# says it fell short.
#
# The numbers come from the engine. The finish that ended the optimiser's run
# (`_ctsem_newton_finish`) hands back the Hessian it ended on, where that was
# evaluated, and what the probe along the untrusted directions found, so this
# certifies without another engine call in the common case; a run that ended
# without a finish -- its iteration cap, a stall -- gets the same numbers from
# `ctsem_endgame`, with no step taken. The arithmetic on them, the gap and the
# verdict, stays here, where it is tested on plain matrices.
#
# What happens next is one rule. `suboptimal` or `notmaximum`, with rounds
# left, resumes the optimiser from the point it reached, under the same rules as
# the first stage, with the stall watch carrying the fit's progress so a
# resumed stage that stops gaining is stopped (`CTSEMStallWatch` in the
# engine). Nothing is switched off and no budget enlarged: the resumed stage's
# own finish is what closes a gap. The rule this replaced resumed after its own
# damped Newton step and negative-curvature step, on a Hessian formed again each
# round, with the predicted-gain stop switched off, the gradient tolerance
# tightened and the iteration cap quadrupled -- which on AnomAuth turned a
# corner of the objective into hours, and from its stored spurious start did
# not return in 1 h 47 min.
#
# Two cases stop instead of resuming. `notstationary`, a flat direction that
# still gains: continuing along it walked AnomAuth S1 from its good optimum into
# the spurious basin (decision 1 of review/OPTIM-consolidation-plan-2026-09-25).
# And a saddle the finish already tried to leave along its negative curvature
# and could not, with nothing else left -- the gap over the trusted directions
# within tolerance: the resumed optimiser would start where the finish
# converged and hand the same point back, and on a Laplace fit each such round
# is a Hessian. A saddle with the gap still open is resumed like any point short
# of its optimum; gated-gaps config B8 (bench, seed 1) was stopped 10.8 nats
# short when this rule did not ask.
#
# `maxtries` is small on purpose: a point still short after two resumes needs a
# different start, which is a decision for whoever runs the fit.
#
# Returns the result to carry on with, the certification that decided, the
# Hessian it decided on with where it was evaluated and how far that is from the
# estimate (so the uncertainty stage need not form the same matrix again), a
# record of each resume, the counts over every stage, and how many Hessians
# were formed.
#' @keywords internal
.ctBackendCorrectResult <- function(result, spec, npar, tolerance = 0.01,
  maxtries = 2L, gradient = "adjoint", optimise = NULL, verbose = 0) {
  spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  history <- list()
  endgame <- NULL
  certification <- NULL
  # The raw coordinates by name, for a verdict that has to say which parameters
  # a flat direction runs through. Positional names if that fails: a message
  # naming `raw[7]` is worse than one naming the parameter and better than none.
  parnames <- tryCatch(.ctBackendRawParameterNames(list(model_spec = spec), npar),
    error = function(e) paste0("raw[", seq_len(npar), "]"))
  # Work done across every stage, not the last one's: a fit that optimised and
  # resumed twice would otherwise report the tail of that and hide the rest.
  stage_counts <- function(r) c(
    iterations = as.numeric(.ctJuliaOr(r$iterations, 0)),
    f_calls = as.numeric(.ctJuliaOr(r$f_calls, 0)),
    g_calls = as.numeric(.ctJuliaOr(r$g_calls, 0)))
  stage_hessians <- function(r) as.integer(.ctJuliaOr(r$newton_hessians, 0L))
  totals <- stage_counts(result)
  hessians <- stage_hessians(result)
  # Where the fit started, so a resumed stage's stall watch measures its
  # progress against the whole fit's rather than its own, which starts at zero.
  started_at <- suppressWarnings(as.numeric(result$trace$objective)[1L])
  if (!length(started_at) || !is.finite(started_at)) {
    started_at <- as.numeric(result$maximum_loglik)[1L]
  }
  for (attempt in seq_len(max(0L, as.integer(maxtries)) + 1L)) {
    est <- as.numeric(result$minimizer)[seq_len(npar)]
    endgame <- .ctBackendEndgameOf(result, npar)
    if (is.null(endgame)) {
      endgame <- .ctBackendEndgameAt(spec, est, gradient = gradient, verbose = verbose)
      if (!is.null(endgame)) hessians <- hessians + endgame$hessians
    }
    if (is.null(endgame)) break
    gap <- .ctBackendOptimGap(endgame$hessian,
      as.numeric(result$gradient)[seq_len(npar)])
    certification <- .ctBackendCertificationRecord(gap, endgame$probe,
      tolerance = tolerance, saturated = isTRUE(result$saturated),
      overshot = isTRUE(result$overshot), parnames = parnames)
    status <- certification$status
    if (!status %in% c("suboptimal", "notmaximum")) break
    if (attempt > as.integer(maxtries) || is.null(optimise)) break
    if (identical(status, "notmaximum") && !isTRUE(result$overshot) &&
        isTRUE(result$newton_ladder_tried) &&
        isTRUE(is.finite(gap$gap) && gap$gap <= tolerance)) {
      break
    }
    value <- as.numeric(result$maximum_loglik)[1L]
    if (verbose > 0) {
      message("Not yet a certified maximum (", status, ", predicted ",
        signif(gap$gap, 3), "); resuming the optimiser")
    }
    resumed <- try(optimise(est, max(0, value - started_at)), silent = TRUE)
    ok <- !inherits(resumed, "try-error") &&
      is.finite(as.numeric(resumed$maximum_loglik)[1L]) &&
      as.numeric(resumed$maximum_loglik)[1L] >= value
    history[[length(history) + 1L]] <- list(attempt = attempt,
      status = status, predicted = gap$gap,
      total_gain = if (ok) as.numeric(resumed$maximum_loglik)[1L] - value else
        NA_real_,
      resumed = ok,
      iterations = if (ok) as.integer(resumed$iterations) else NA_integer_,
      # Why the resumed stage stopped, and whether it was the carried progress
      # watch -- a resume that stopped gaining -- that stopped it.
      stop_reason = if (ok) as.character(.ctJuliaOr(resumed$stop_reason,
        NA_character_))[1L] else NA_character_,
      stopped_progress = if (ok) isTRUE(resumed$stopped_by_progress) else NA)
    # Never accept a resume that did not improve on what we had: the fit can
    # only move forward here.
    if (!ok) break
    totals <- totals + stage_counts(resumed)
    hessians <- hessians + stage_hessians(resumed)
    result <- resumed
  }
  list(result = result, certification = certification,
    hessian = if (is.null(endgame)) NULL else endgame$hessian,
    evaluated_at = if (is.null(endgame)) NULL else endgame$evaluated_at,
    distance = if (is.null(endgame)) NA_real_ else endgame$distance,
    corrections = history, totals = totals, hessians = hessians)
}

# What the optimiser's own stopping rule is set to, and when it is on at all.
#
# Inside the bar, not equal to it. The two do estimate the same quantity -- the
# line search's `g'Bg` is L-BFGS's estimate of the objective still available and
# `g'H^-1 g` is the exact one -- but `B` is a limited-memory approximation, so
# setting the optimiser's target *at* the bar means every fit whose proxy is
# even slightly optimistic fails the check and takes a correction. A correction
# is a Hessian and a resumed optimisation, and measured against simply
# optimising further it is the expensive way to gain the last of the objective:
# on that same fixture, reaching the point where the diagnostics work cost 2729
# objective calls through the optimiser and 5197 through corrections.
#
# So the optimiser aims two orders inside the bar and the correction is what it
# was meant to be -- the rare case where a fit genuinely stopped short, not the
# normal route to precision.
#
# Rechecked on 480 fits under the preconditioner and the short first step --
# `dev/simstudies/simstudy-gaptol.R`, four measurement types and both
# intoverpop routes, paired on the same data and the same starting values, with
# the constant swept over `bar/100`, `bar/10`, the bar itself and off. Scored on
# objective calls, which is the deterministic half of the trade; the wall clock
# there is not readable, because the four settings run back to back in one
# worker with the default first, so it absorbs Julia's per-shape specialisation.
#
#   setting     median f_calls   cells landing >0.01 below the reference
#   bar/100 (default)   32.0     0   (largest shortfall 1.3e-08)
#   bar/10              31.5     0   (largest shortfall 1.4e-07)
#   at the bar          30.5     1   (2.90 log units, 2 corrections, and it
#                                     still reported not converged)
#   off                 44.5     0
#
# Which is the same answer as the fixture, from the other side: moving the rule
# out to the bar buys about 5% of the calls and puts a fit in sixty into exactly
# the loop this constant exists to avoid -- proxy fails the exact check, two
# Hessians and two resumed optimisations, and it lands 2.9 log units low
# anyway. `bar/10` is safe and saves 1.5%, which is not a reason to move.
#
# The rule as a whole is worth having, incidentally: 32 calls against 44.5 with
# it off, so it removes about 28% of the optimiser's work for a shortfall five
# orders of magnitude inside the bar.
#
# The licence applies to the *default*, not to a request. Nothing will check a
# stop when certification is off -- `certify = FALSE`, or the state-explicit
# route, whose only curvature is the profile's -- and the proxy can stop but
# never certify, so defaulting it on there would be a fit that quit early with
# nothing to say so. A caller who names `innergaptol` has asked for exactly
# that and is given it: accepting the argument and ignoring it would be worse
# than either answer.
#
# `estonly` used to be in that list and is not. It asks for the estimate
# without the uncertainty and correction phases, which is a statement about
# what happens *after* the optimisation -- and it was silently changing the
# optimisation too, so the same model fitted with and without it ran different
# stopping rules and could stop in different places. Every other reading of
# `estonly` in the package skips a post-fit step; this one did not, which is
# the shape of the trap CLAUDE.md names. `certify = FALSE` stays, because that
# one does say the curvature will not be computed.
#' @keywords internal
.ctBackendInnerGapTol <- function(optimcontrol = list(), intoverstates = TRUE) {
  explicit <- optimcontrol$innergaptol
  if (!is.null(explicit)) {
    value <- suppressWarnings(as.numeric(explicit)[1L])
    return(if (is.finite(value) && value >= 0) value else 0)
  }
  if (!isTRUE(intoverstates)) return(0)
  if (identical(optimcontrol$certify, FALSE)) return(0)
  bar <- .ctBackendConvergeTol(optimcontrol)
  if (!is.finite(bar) || bar <= 0) 0 else bar / 100
}
