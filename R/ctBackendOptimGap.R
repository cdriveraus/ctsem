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
# ## Three outcomes, not two
#
#   certified      the gap is below tolerance and the excluded directions hold
#                  no material likelihood
#   uncertified    the gap is small but something outside the trusted subspace
#                  is not: flat directions with a live gradient, or a
#                  transform that has saturated, where the gradient underflows
#                  to zero and every tolerance passes for the wrong reason
#   notmaximum     a direction of genuine negative curvature, where the point
#                  is a saddle whatever the gap says
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
.ctBackendInformationSplit <- function(hessian, rtol = 1e-12,
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
.ctBackendOptimGap <- function(hessian, gradient, rtol = 1e-12,
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
    nnegative = sum(split$negative), ok = TRUE)
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
    return(list(status = "uncertified", certified = FALSE,
      reason = paste0("stepping along the directions with no curvature gains ",
        signif(residual_gain, 3), " log likelihood, so the estimate is not ",
        "stationary in a direction the data does not determine")))
  }
  if (isTRUE(saturated) && isTRUE(gap$nflat > 0L)) {
    return(list(status = "uncertified", certified = FALSE,
      reason = paste0("a parameter transform has saturated, where the ",
        "gradient underflows to zero and the curvature with it, so no ",
        "tolerance here means anything for that coordinate")))
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

# The tolerance, in objective units. One number, settable per fit.
#
# 0.01 is a hundredth of a log likelihood unit, well inside any difference that
# would change a conclusion -- an AIC comparison turns on 2 -- and far enough
# above the arithmetic's own noise to be reachable. It is not scaled by the
# parameter count: the gap is the joint improvement still available, which is
# the quantity of interest whatever the dimension.
#' @keywords internal
.ctBackendGapTolerance <- function(fit, default = 0.01) {
  control <- fit$args$optimcontrol
  value <- if (is.null(control)) NULL else control$gaptol
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(value) || value <= 0) return(default)
  value
}

# The certification for one fit, at its own estimate, against one Hessian.
#
# The probe is what makes the flat directions safe to exclude from the gap:
# their gradient is stepped along and the objective measured, rather than a
# norm being compared to a threshold that would carry the units this whole
# approach exists to remove.
#' @keywords internal
.ctBackendCertification <- function(fit, hessian, tolerance = 0.01) {
  gradient <- as.numeric(fit$estimate$gradient)
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
    saturated = isTRUE(fit$estimate$saturated),
    overshot = isTRUE(fit$estimate$overshot))
  list(status = verdict$status, certified = verdict$certified,
    reason = verdict$reason, tolerance = tolerance,
    gap = gap$gap, lambda = gap$lambda,
    residual_gain = if (is.null(probe)) 0 else probe$gain,
    residual_length = if (is.null(probe)) 0 else probe$length,
    ntrusted = gap$ntrusted, nflat = gap$nflat, nnegative = gap$nnegative,
    # The displacement the gap predicts, kept so a caller can try it without
    # decomposing the Hessian again.
    step = gap$step)
}

# The one convergence statement a fit makes.
#
# Three cases, and they are genuinely different things to tell a user:
#
#   certified        nothing. The optimizer's own stopping rule was superseded
#                    by a criterion it does not know about, and repeating its
#                    complaint would send someone chasing a fit that is right.
#   not certified    what the curvature says, in objective units, because that
#                    is actionable where a gradient is not.
#   no certification the gradient, and the fact that nothing checked further.
#' @keywords internal
.ctBackendCertifyWarn <- function(fit) {
  certification <- fit$uncertainty$certification
  pending <- isTRUE(fit$estimate$convergence_pending)
  if (is.null(certification) || !length(certification$status)) {
    if (!pending) return(invisible(NULL))
    warning("The optimizer stopped without meeting its convergence criterion: ",
      "largest gradient ", signif(as.numeric(fit$estimate$gradient_norm), 3),
      ". No curvature was computed, so how far this is from the optimum is ",
      "unknown -- fit without optimcontrol$estonly, or call ",
      "ctOptimUncertainty(), to have it certified.", call. = FALSE)
    return(invisible(NULL))
  }
  if (isTRUE(certification$certified)) return(invisible(NULL))
  warning("This fit is not certified as converged: ", certification$reason,
    ". See fit$uncertainty$certification.", call. = FALSE)
  invisible(NULL)
}
