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
# ## Four outcomes, not two
#
#   certified      the gap is below tolerance and the excluded directions hold
#                  no material likelihood
#   suboptimal     the gap says the optimum is measurably above this estimate
#   notstationary  a flat direction with a live gradient: stepping along it
#                  gains likelihood, so this is not a maximum, and the gap
#                  cannot see it because that direction was excluded
#   unidentified   a transform has saturated, where the gradient underflows to
#                  zero and the curvature with it, so no tolerance means
#                  anything for that coordinate. The point is still a maximum
#   notmaximum     a direction of genuine negative curvature, where the point
#                  is a saddle whatever the gap says
#
# `notstationary` and `unidentified` were one status, and they are not one
# finding: the first is a fit nobody should use and the second is a fit with a
# result in it -- a population scale with no individual differences behind it is
# the usual cause, and reporting that as a failure to converge is the mistake
# `test-julia-convergence.R` exists to prevent.
#
# The saturation verdict is consulted rather than recomputed. A saturated
# parameter has both a zero gradient and no curvature, so it lands in the
# excluded subspace with nothing to report, and a certification that did not
# ask would contradict the warning the fit already carries.

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
    # The smallest curvature still trusted, which is what sets how tight a
    # gradient has to be before the gap can be under a given tolerance.
    lambda_min = if (length(values)) min(values) else NA_real_, ok = TRUE)
}

# What the excluded directions are actually worth, in log likelihood.
#
# A norm of the leftover gradient cannot answer this: it has the units the rest
# of this file exists to avoid, and in a direction with no curvature the
# quadratic model that would convert it into a likelihood is precisely the one
# that does not hold. So step along it and look.
#
# The step lengths are in raw parameter units, where ctsem's coordinates are
# standardised by construction -- priors are normal(0,1) and each transform
# carries its own scale -- so a quarter, one and four span "a small move" to "a
# large one" without needing to know anything about the model. The best actual
# improvement found is what is reported; a direction that gains nothing at any
# of them holds nothing worth continuing for.
#' @keywords internal
.ctBackendOptimGapProbe <- function(evaluate, at, direction, value,
  lengths = c(0.25, 1, 4)) {
  norm <- sqrt(sum(direction^2))
  if (!is.finite(norm) || norm <= 0) {
    return(list(gain = 0, length = 0, ok = TRUE))
  }
  unit <- direction / norm
  best <- 0
  bestlength <- 0
  for (len in lengths) {
    probe <- try(evaluate(at + len * unit), silent = TRUE)
    if (inherits(probe, "try-error")) next
    probe <- as.numeric(probe)[1L]
    if (!is.finite(probe)) next
    if (probe - value > best) {
      best <- probe - value
      bestlength <- len
    }
  }
  list(gain = best, length = bestlength, ok = TRUE)
}

# The verdict, in the terms a reader needs.
#
# `saturated` is the fit's own, not recomputed here: a saturated transform
# reports a zero gradient and no curvature, so it passes every test in this
# file for the wrong reason, and the fit already says so.
#' @keywords internal
.ctBackendCertify <- function(gap, probe = NULL, tolerance = 0.01,
  saturated = FALSE, overshot = FALSE) {
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
      reason = paste0("stepping along the directions with no curvature gains ",
        signif(residual_gain, 3), " log likelihood, so the estimate is not ",
        "stationary in a direction the data does not determine")))
  }
  if (isTRUE(saturated) && isTRUE(gap$nflat > 0L)) {
    return(list(status = "unidentified", certified = FALSE,
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
#' @keywords internal
.ctBackendGapTolerance <- function(fit, default = 1e-6) {
  .ctBackendConvergeTol(fit$args$optimcontrol, default = default)
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
# `4 * npar` value-only evaluations per optimisation stage, and on a large model
# that is a visible fraction of the fit. What it buys is that a fit stopped in a
# degenerate corner says so instead of reporting convergence.
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
#' @keywords internal
.ctBackendStallWindow <- function(optimcontrol = list(), default = 80L) {
  value <- if (is.null(optimcontrol)) NULL else optimcontrol$stallwindow
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.integer(value)[1L])
  if (is.na(value) || value < 0L) default else value
}

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
# here would cost `4 * npar` objective evaluations for an answer already in
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

# The certification for one fit, at its own estimate, against one Hessian.
#
# The probe is what makes the flat directions safe to exclude from the gap:
# their gradient is stepped along and the objective measured, rather than a
# norm being compared to a threshold that would carry the units this whole
# approach exists to remove.
#' @keywords internal
.ctBackendCertification <- function(fit, hessian, tolerance = 0.01) {
  gradient <- as.numeric(fit$optim$gradient)
  gap <- .ctBackendOptimGap(hessian, gradient)
  probe <- NULL
  if (isTRUE(gap$ok) && isTRUE(gap$residual_norm > 0)) {
    evaluate <- try(.ctBackendLpgFunc(fit, gradient = FALSE), silent = TRUE)
    if (!inherits(evaluate, "try-error")) {
      probe <- .ctBackendOptimGapProbe(evaluate, as.numeric(fit$estimate$raw),
        gap$residual, as.numeric(fit$estimate$logposterior)[1L])
    }
  }
  verdict <- .ctBackendCertify(gap, probe, tolerance = tolerance,
    saturated = isTRUE(fit$optim$saturated),
    overshot = isTRUE(fit$optim$overshot))
  list(status = verdict$status, certified = verdict$certified,
    reason = verdict$reason, tolerance = tolerance,
    # The same `$verdict` shape `ctLaplaceCheck()` and `ctParticleLik()` report;
    # see R/ctFitGap.R. Set here *and* in `.ctBackendCorrectResult()`, because
    # the certification list is assembled in both. This one runs later, from the
    # fit, and overwrites that one -- so a field added only there never reaches
    # a reader, which is how the first attempt at this got a NULL. The two
    # assemblies are not a drop-in merge: one takes an optimiser `result` and
    # the other a fit.
    verdict = .ctBackendCertifyGap(gap, tolerance),
    gap = gap$gap, lambda = gap$lambda,
    residual_gain = if (is.null(probe)) 0 else probe$gain,
    residual_length = if (is.null(probe)) 0 else probe$length,
    ntrusted = gap$ntrusted, nflat = gap$nflat, nnegative = gap$nnegative,
    # The displacement the gap predicts, kept so a caller can try it without
    # decomposing the Hessian again.
    step = gap$step)
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
  fit$optim$converged <- certification$status %in%
    c("certified", "unidentified")
  # Superseded rather than answered: the optimizer's complaint was held for
  # this measurement, and the measurement has now been made. Leaving it would
  # make `.ctBackendCertifyWarn()` warn about a gradient on a fit whose
  # curvature already said better.
  fit$optim$convergence_pending <- NULL
  fit
}

# The one convergence statement a fit makes.
#
# Three cases, and they are genuinely different things to tell a user:
#
#   unidentified     the coordinate the data does not determine, once,
#                    pointing at the identifiability report rather than
#                    repeating it here.
#   certified        nothing. The optimizer's own stopping rule was superseded
#                    by a criterion it does not know about, and repeating its
#                    complaint would send someone chasing a fit that is right.
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
  # A maximum with a coordinate the data does not determine is a finding, and
  # the finding is `fit$identifiability`'s to report. Warning "not converged"
  # here is what sent 45 of 64 good fits back to be re-run.
  if (identical(certification$status, "unidentified")) {
    warning("This fit is a maximum, but ", certification$reason, ".",
      call. = FALSE)
    return(invisible(NULL))
  }
  warning("This fit is not certified as converged: ", certification$reason,
    ". See fit$uncertainty$certification.", call. = FALSE)
  invisible(NULL)
}

# The curvature at one point, from the spec rather than from a fit.
#
# `.ctBackendHessian()` is the same call with a fit around it; this exists
# because the correction runs before there is a fit to pass, and both go
# through the one engine entry point rather than two.
#' @keywords internal
.ctBackendHessianAt <- function(spec, est, gradient = "adjoint") {
  # Classed here because the spec reaches this point unclassed: `ctFit()`
  # carries `model_spec` as a plain list and `.ctJuliaOptimise()` does the same
  # thing for the same reason.
  spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  name <- if (identical(gradient, "forward")) "ctsem_hessian_forward" else
    "ctsem_hessian"
  available <- isTRUE(tryCatch(is.function(module[[name]]),
    error = function(e) FALSE))
  if (!available) return(NULL)
  est <- as.numeric(est)
  out <- try(.ctBackendJuliaValue(module[[name]](.ctJuliaObjective(spec),
    .ctJuliaVector(est))), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  matrix(as.numeric(out), nrow = length(est), ncol = length(est))
}

# Certify the optimiser's result, and continue from a corrected point when it
# falls short.
#
# Returns the result to carry on with -- the original when nothing was wrong or
# nothing could be improved -- plus the certification that decided, the Hessian
# it decided on (so the uncertainty stage need not compute the same matrix
# again), and a record of every correction attempted.
#
# `maxtries` is small on purpose. Each round costs a Hessian and an
# optimisation, and a point that is still short after two corrections is not
# going to be rescued by a third: what it needs is a different starting value,
# which is a decision for whoever is running the fit.
#' @keywords internal
.ctBackendCorrectResult <- function(result, spec, npar, tolerance = 0.01,
  maxtries = 2L, gradient = "adjoint", optimise = NULL, maxiter = NA_integer_,
  gtol = 1e-8, verbose = 0) {
  # Carried across attempts: a stage that needed a bigger budget once needs it
  # again, and a tolerance tightened once must not slacken on the next round.
  overrides <- list()
  spec <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  # The laplace objective when the model has one: `.ctJuliaObjective()` wraps
  # the marginal objective for such a spec, so the curvature certified here is
  # the curvature of what was actually maximised.
  objective <- .ctJuliaObjective(spec)
  value_at <- function(pars) {
    # `ctJuliaEvaluate()`'s call, without needing a fit to hang it on.
    out <- try(JuliaConnectoR::juliaGet(module$ctsem_evaluate(objective,
      .ctJuliaNumericVector(as.numeric(pars)), gradient = FALSE,
      contributions = FALSE, gradient_method = gradient)), silent = TRUE)
    if (inherits(out, "try-error")) return(NA_real_)
    value <- suppressWarnings(as.numeric(out$value)[1L])
    if (!length(value)) NA_real_ else value
  }
  history <- list()
  hessian <- NULL
  certification <- NULL
  # Work done across every stage, not the last one's. `iterations`, `f_calls`
  # and `g_calls` come off whichever Optim run finished, so a fit that
  # optimised, corrected, resumed and corrected again reported the tail of that
  # and hid the cost of the rest. The totals are what a comparison between
  # stopping rules has to be made on.
  stage_counts <- function(r) c(
    iterations = as.numeric(.ctJuliaOr(r$iterations, 0)),
    f_calls = as.numeric(.ctJuliaOr(r$f_calls, 0)),
    g_calls = as.numeric(.ctJuliaOr(r$g_calls, 0)))
  totals <- stage_counts(result)
  hessians <- 0L
  for (attempt in seq_len(max(0L, as.integer(maxtries)) + 1L)) {
    est <- as.numeric(result$minimizer)[seq_len(npar)]
    hessian <- .ctBackendHessianAt(spec, est, gradient = gradient)
    if (is.null(hessian)) break
    # Counted, because a Hessian is `ceil(npar / chunksize)` sweeps and is the
    # price of certifying at all -- paid once here and reused by the
    # uncertainty stage, so a reader comparing arrangements needs to see it.
    hessians <- hessians + 1L
    gap <- .ctBackendOptimGap(hessian, as.numeric(result$gradient)[seq_len(npar)])
    probe <- NULL
    if (isTRUE(gap$ok) && isTRUE(gap$residual_norm > 0)) {
      probe <- .ctBackendOptimGapProbe(value_at, est, gap$residual,
        as.numeric(result$maximum_loglik)[1L])
    }
    certification <- .ctBackendCertify(gap, probe, tolerance = tolerance,
      saturated = isTRUE(result$saturated), overshot = isTRUE(result$overshot))
    certification$verdict <- .ctBackendCertifyGap(gap, tolerance)
    certification$gap <- gap$gap
    certification$lambda <- gap$lambda
    certification$ntrusted <- gap$ntrusted
    certification$nflat <- gap$nflat
    certification$nnegative <- gap$nnegative
    certification$residual_gain <- if (is.null(probe)) 0 else probe$gain
    certification$tolerance <- tolerance
    # Continue for either reason the curvature gives. `suboptimal` is objective
    # left on the table; `notmaximum` is a direction of negative curvature,
    # which is a stronger reason to carry on and not a verdict to stop at --
    # measured on a six-subject fixture, the cheap rule stopped at a point with
    # a curvature of -0.327 and a log likelihood 0.3 worse than the same fit
    # run on, and the loop declined to act because it only looked for
    # `suboptimal`. The trusted subspace still gives an ascent direction, and
    # where it gives little the resumed stage with a tightened rule is what
    # moves off the saddle.
    if (!certification$status %in% c("suboptimal", "notmaximum")) break
    if (attempt > as.integer(maxtries) || is.null(optimise)) break

    # Damped, and accepted only on an increase that is actually observed and
    # large enough to be one. `directional` is `g'd`, which is `2 * gap` for
    # the Newton step over the trusted subspace, so the first-order gain at
    # `alpha` is `alpha * directional`; halving stops when that falls below what
    # the objective can represent, because succeeding there is indistinguishable
    # from not trying. Derived from the arithmetic rather than from a rung
    # count, which would be fitted to whichever model it was measured on.
    value <- as.numeric(result$maximum_loglik)[1L]
    directional <- sum(as.numeric(result$gradient)[seq_len(npar)] * gap$step)
    stepped <- .ctBackendDampedStep(value_at, est, gap$step, value, directional)
    accepted <- stepped$accepted
    best <- stepped$achievable
    if (is.null(accepted)) {
      # The step predicted an improvement and no scaling of it delivered one,
      # which says the quadratic model is wrong here rather than that the fit
      # is finished. Recorded, because it is the interesting case.
      #
      # Then resume from where we are rather than stopping: at a saddle the
      # trusted subspace can have nothing to offer while the point is still not
      # a maximum, and the optimiser with a tightened rule is what leaves it.
      # Stopping here because one step failed would be giving up for the wrong
      # reason.
      history[[length(history) + 1L]] <- list(attempt = attempt,
        predicted = gap$gap, step_gain = NA_real_, ratio = NA_real_,
        achievable = best, total_gain = NA_real_, alpha = NA_real_,
        resumed = FALSE)
      accepted <- list(par = est, value = value, alpha = 0)
    }
    # Tighten whatever ended the last stage, or the resume stops there again.
    # Which one it was is not a guess: `iterations` is the engine's own count,
    # not Optim's, which stops being updated when a callback ends the run, and
    # `stopped_by_gap` says outright whether the cheap rule was what stopped it.
    hitcap <- is.finite(maxiter) && maxiter > 0 &&
      as.integer(result$iterations) >= as.integer(maxiter)
    overrides <- .ctBackendResumeOverrides(overrides, hitcap = hitcap,
      stoppedbygap = isTRUE(result$stopped_by_gap),
      lambda_min = gap$lambda_min, npar = npar, tolerance = tolerance,
      maxiter = maxiter, gtol = gtol)
    if (verbose > 0) {
      message("Continuing from a Newton correction: predicted ",
        signif(gap$gap, 3), ", step gained ", signif(accepted$value - value, 3),
        if (hitcap) paste0(", iterations raised to ", overrides$maxiter) else
          paste0(", gradient tolerance tightened to ",
            signif(overrides$g_tol, 3),
            if (identical(overrides$innergaptol, 0))
              ", predicted-gain stop switched off" else ""))
    }
    resumed <- try(optimise(accepted$par, overrides), silent = TRUE)
    ok <- !inherits(resumed, "try-error") &&
      is.finite(as.numeric(resumed$maximum_loglik)[1L]) &&
      as.numeric(resumed$maximum_loglik)[1L] >= value
    # `ratio` is the step's own gain over what the quadratic predicted, and
    # nothing else: it exists to say whether the local quadratic could be
    # trusted here, and measuring it after the resume answered a different
    # question -- it read above one, because L-BFGS carries on improving past
    # the quadratic's optimum. `total_gain` is what the pair were worth.
    history[[length(history) + 1L]] <- list(attempt = attempt,
      predicted = gap$gap,
      step_gain = accepted$value - value,
      ratio = (accepted$value - value) / gap$gap,
      achievable = best,
      total_gain = if (ok) as.numeric(resumed$maximum_loglik)[1L] - value else
        accepted$value - value,
      alpha = accepted$alpha, resumed = ok,
      # What the resumed stage was given, so a fit that still comes up short
      # says what was already tried.
      resume_maxiter = .ctJuliaOr(overrides$maxiter, maxiter),
      resume_gtol = .ctJuliaOr(overrides$g_tol, gtol),
      # NA when it was left alone, which is the case where the stage ended at
      # something other than the predicted-gain rule.
      resume_innergaptol = .ctJuliaOr(overrides$innergaptol, NA_real_),
      # What stopped the resumed stage. A resume that returns after a couple of
      # iterations having met the optimiser's own criterion is the futile case:
      # the rule that was just shown to be inadequate is what ended it, and
      # another correction would meet the same wall. That is what switching the
      # predicted-gain rule off above prevents; this is how a resume that still
      # comes up short for some other reason says so.
      iterations = if (ok) as.integer(resumed$iterations) else NA_integer_,
      stopped_converged = if (ok) isTRUE(resumed$converged) else NA)
    # Never accept a resume that did not improve on what we had: the corrected
    # point is already better than the estimate, so the fit can only move
    # forward here.
    if (ok) {
      totals <- totals + stage_counts(resumed)
      result <- resumed
    } else break
  }
  list(result = result, certification = certification, hessian = hessian,
    corrections = history, totals = totals, hessians = hessians)
}

# What to give the resumed stage, so that it does not stop where the last one
# did.
#
# Separated from the loop above for the same reason as `.ctBackendDampedStep`:
# the decision is branching on three numbers, and the cases worth pinning -- a
# stage that ran out of iterations, one the cheap rule ended, one whose derived
# bound is looser than the rule already in force -- are three lines each and
# take no engine at all.
#
# There are three rules that can end a stage and each is answered in its own
# terms.
#
# **The iteration cap.** Quadrupled, capped at 1e6. A stage that needed a
# bigger budget once needs it again, which is why `overrides` is carried across
# attempts rather than rebuilt.
#
# **The gradient bound.** Tightened to `.ctBackendGapGradientTolerance()`, the
# gradient below which the gap cannot exceed the target -- and only ever
# tightened, because a derivation that came out looser than the rule already in
# force would be licensing the stop that just failed.
#
# **The predicted-gain rule.** Switched off, not tightened by some factor. The
# proxy is `1/2 g'Bg` against the exact `1/2 g'H^-1 g`, and `B` is limited
# memory, so no factor provably makes a second failure impossible -- picking
# one would be fitting a constant to whichever model it was measured on. What
# can be said is that on *this* model the proxy has just been shown wrong, so
# it has forfeited its licence here; the resumed stage still has the gradient
# bound, the cap, and the exact check that follows it.
#
# That last one was missing, and it is what made a second correction possible.
# Measured on the one cell in 480 of `dev/simstudies/simstudy-gaptol.R` that
# needed more than one: with the rule left in force both resumes ran a single
# iteration and reported the optimiser's own criterion met -- the futile case
# named where the history is recorded -- and the fit ended 2.9 log units low
# and not converged. With it switched off the same fit takes one correction,
# its resume runs 108 iterations for 2.9 log units, and it lands on the
# log likelihood the default stopping rule reaches.
#
# The cap branch is exclusive: a stage that exhausted its iterations did not
# stop at either tolerance, so neither is what needs loosening.
#' @keywords internal
.ctBackendResumeOverrides <- function(overrides = list(), hitcap = FALSE,
  stoppedbygap = FALSE, lambda_min = NA_real_, npar = 1L, tolerance = 0.01,
  maxiter = NA_integer_, gtol = 1e-8) {
  if (isTRUE(hitcap)) {
    overrides$maxiter <- as.integer(min(4 * as.numeric(
      .ctJuliaOr(overrides$maxiter, maxiter)), 1e6))
    return(overrides)
  }
  if (isTRUE(stoppedbygap)) overrides$innergaptol <- 0
  needed <- .ctBackendGapGradientTolerance(lambda_min, npar, tolerance)
  if (is.finite(needed)) {
    overrides$g_tol <- min(needed, .ctJuliaOr(overrides$g_tol, gtol))
  }
  overrides
}

# Backtrack along an ascent direction until the objective actually increases.
#
# Separated from the loop above so that it can be tested on a function rather
# than on a fit: the ladder is arithmetic, and the cases worth pinning -- a step
# far too long, a direction that offers nothing, an increase too small to be one
# -- are constructed in two lines each and take no engine at all.
#
# `directional` is `g'd`, so the first-order gain at `alpha` is
# `alpha * directional`. Halving stops once that falls below what the objective
# can represent (`|f| * eps`), because succeeding below it is indistinguishable
# from not trying; there is no rung count, which would be fitted to whichever
# model it was measured on. Acceptance is Armijo at 1e-4 rather than any
# increase, so a rounding error is not recorded as a correction.
#' @keywords internal
.ctBackendDampedStep <- function(value_at, at, step, value, directional,
  c1 = 1e-4) {
  achievable <- 0
  if (!is.finite(directional) || directional <= 0) {
    return(list(accepted = NULL, achievable = achievable))
  }
  floor <- max(abs(value), 1) * .Machine$double.eps
  alpha <- 1
  while (alpha * directional > floor) {
    got <- value_at(at + alpha * step)
    if (is.finite(got)) {
      achievable <- max(achievable, got - value)
      # Armijo, and representable. Sufficient increase scales with the step,
      # so at a small enough alpha an increase of 1e-14 satisfies it -- and
      # accepting that spends a whole resumed optimisation on a point that is
      # not distinguishable from where it started. The floor is the objective's
      # own resolution, the same bound that ends the loop.
      if (got - value >= c1 * alpha * directional && got - value > floor) {
        return(list(accepted = list(par = at + alpha * step, value = got,
          alpha = alpha), achievable = achievable))
      }
    }
    alpha <- alpha / 2
  }
  list(accepted = NULL, achievable = achievable)
}

# The gradient tolerance that would have made this gap small enough.
#
# `gap = 1/2 g'H^-1 g` is at most `n |g|_inf^2 / (2 lambda_min)` over the
# trusted subspace, so a gradient below `sqrt(2 tol lambda_min / n)` cannot
# leave a gap above `tol`. Conservative -- it charges every coordinate the worst
# direction's curvature -- and that is the right side to be on for a rule whose
# job is to stop a resumed stage halting where the last one did.
#
# NA when there is no trusted curvature to derive it from, which is a model with
# nothing to certify rather than one needing a tighter tolerance.
#' @keywords internal
.ctBackendGapGradientTolerance <- function(lambda_min, npar, tolerance) {
  if (!isTRUE(is.finite(lambda_min)) || lambda_min <= 0 || npar < 1) {
    return(NA_real_)
  }
  sqrt(2 * tolerance * lambda_min / npar)
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
