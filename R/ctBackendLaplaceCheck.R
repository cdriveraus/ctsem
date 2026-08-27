# How much of a Laplace fit is the approximation.
#
# `intoverpop='laplace'` replaces each subject's integral over its random
# effects by a single Gaussian fitted at the mode. That is exact when the
# varying parameters have identity transforms and enter the state mean
# linearly, and it is not exact otherwise -- most visibly for a random effect on
# DRIFT, whose transform is `-log1p_exp`.
#
# The error matters because it is not constant. It grows with the population
# scale, so it tilts the profile in that direction and the maximum moves down.
# On a simulated 40-subject, 24-wave model with a random DRIFT, the Laplace
# profile peaks at a population sd of 0.67 where a nine-node adaptive
# Gauss-Hermite rule peaks at 0.89, against a generating 1.00 -- two thirds of
# the shortfall belongs to the approximation rather than to the data.
#
# The same quadrature that measures the error also corrects it, and to first
# order it does so without refitting: at the Laplace estimate the Laplace
# gradient is zero, so the quadrature objective's gradient there is the
# derivative of the *gap* between them, and one Newton step against the
# curvature already computed for the standard errors gives the shift. That
# costs `2 * npar` quadrature evaluations; a full refit against the quadrature
# objective costs the same per gradient and many gradients.

#' Measure, and optionally correct, the Laplace approximation error in a fit
#'
#' Recomputes a \code{backend='julia'}, \code{intoverpop='laplace'} fit's log
#' marginal likelihood by adaptive Gauss-Hermite quadrature over the same
#' random-effect integral the Laplace approximation approximates, and reports
#' how far the estimate would move if that integral were taken exactly.
#'
#' The quadrature is adaptive in the usual sense: it is centred at each
#' subject's inner mode and scaled by the inverse of the inner curvature there,
#' both of which the fit has already computed. With \code{nodes = 1} it
#' reproduces the Laplace value exactly, which is what makes the comparison
#' meaningful -- the two are integrating the same thing by different rules.
#'
#' The correction is a single Newton step, \code{delta = solve(-hessian,
#' gradient of (quadrature - laplace))}, evaluated at the estimate. Read it
#' against the standard errors: \code{delta / se} well below one says the
#' approximation is not what limits the answer, and near or above one says the
#' point estimate is approximation-limited and the interval will not cover
#' whatever width it has.
#'
#' Nested groupings are handled. A group's integral does not factor over its
#' members, but it does factor \emph{conditionally} -- given the group effect
#' the members are independent -- so the rule recurses over the same block tree
#' the fit already builds, at a cost linear in the number of groups at each
#' level rather than exponential in the number of members.
#'
#' @param fit A \code{ctJuliaFit} fitted with \code{intoverpop='laplace'}.
#' @param nodes Quadrature nodes per random effect. Cost is
#'   \code{nodes^k} process log likelihoods per block, where \code{k} is that
#'   level's number of random effects; 5 is enough to locate the maximum in the
#'   cases tested, 9 to settle the value.
#' @param correction Compute the first-order correction to the estimate. Costs
#'   \code{2 * npar} quadrature evaluations, and needs the fit's Hessian, so it
#'   is skipped when uncertainty was not computed.
#' @param step Finite-difference step for the gap gradient, on the raw scale.
#' @param refine Refit against the quadrature objective rather than correcting
#'   linearly. Much slower -- each gradient is another \code{2 * npar}
#'   quadrature evaluations -- and worth it only when \code{delta / se} is
#'   large enough that a linear correction is not credible.
#' @param maxiter Iteration cap for \code{refine}.
#' @param cores Engine threads for the quadrature.
#' @param verbose Integer; 1 or more prints progress.
#'
#' @return A list with \code{gap} (quadrature minus Laplace log marginal at the
#'   estimate), \code{gap_per_subject}, \code{nodes}, and, when
#'   \code{correction} is \code{TRUE}, a \code{parameters} data frame with one
#'   row per raw parameter giving \code{estimate}, \code{delta},
#'   \code{corrected}, \code{se} and \code{delta_se}. With \code{refine} it also
#'   carries \code{refined}.
#' @export
ctLaplaceCheck <- function(fit, nodes = 5L, correction = TRUE, step = 1e-3,
  refine = FALSE, maxiter = 50L, cores = NULL, verbose = 0L) {

  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctLaplaceCheck applies to backend='julia' fits.", call. = FALSE)
  }
  if (is.null(fit$model_spec$laplace)) {
    stop("ctLaplaceCheck applies to fits made with intoverpop='laplace'. ",
      "The augmented route carries the random effects in the state and has no ",
      "separate integral to check.", call. = FALSE)
  }
  est <- as.numeric(fit$estimate$raw)
  module <- .ctJuliaModule(fit$model_spec$project)
  objective <- .ctJuliaObjective(fit)
  # Restored on exit: this is a diagnostic, and it has no business changing how
  # fast everything else in the session runs afterwards.
  if (!is.null(cores)) {
    previous <- tryCatch(as.integer(.ctBackendJuliaValue(JuliaConnectoR::juliaEval(
      "ContinuousTimeSEM.ctsem_max_chunks().max_chunks"))), error = function(e) NA_integer_)
    JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!",
      as.integer(max(1L, cores)))
    if (!is.na(previous)) {
      on.exit(try(JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_set_max_chunks!",
        previous), silent = TRUE), add = TRUE)
    }
  }

  if (verbose > 0) message("Quadrature at the estimate (", nodes, " nodes)")
  quadrature <- .ctBackendJuliaValue(module$ctsem_laplace_quadrature(
    objective, .ctJuliaNumericVector(est), nodes = as.integer(nodes))$value)
  laplacevalue <- .ctBackendJuliaValue(module$ctsem_laplace_evaluate(
    objective, .ctJuliaNumericVector(est), gradient = FALSE)$value)
  nsubjects <- length(fit$model_spec$subject_starts)
  out <- list(gap = quadrature - laplacevalue,
    quadrature = quadrature, laplace = laplacevalue,
    gap_per_subject = (quadrature - laplacevalue) / max(1L, nsubjects),
    nodes = as.integer(nodes), nsubjects = nsubjects)
  class(out) <- "ctLaplaceCheck"

  # Backend fits keep their uncertainty at `fit$uncertainty` and their raw
  # covariance at `fit$estimate$cov`; the Stan path uses `fit$stanfit$...`.
  # Both are checked so this reads either without the caller knowing which.
  hessian <- fit$uncertainty$hessian
  if (is.null(hessian)) hessian <- fit$stanfit$uncertainty$hessian
  if (!correction) return(out)
  if (is.null(hessian)) {
    warning("No Hessian on the fit, so the correction cannot be formed. ",
      "Run ctOptimUncertainty(fit) first, or use correction=FALSE.",
      call. = FALSE)
    return(out)
  }

  if (verbose > 0) message("Gap gradient (", 2 * length(est), " quadrature evaluations)")
  result <- JuliaConnectoR::juliaGet(module$ctsem_laplace_correction(
    objective, .ctJuliaNumericVector(est),
    JuliaConnectoR::juliaPut(as.matrix(hessian)),
    nodes = as.integer(nodes), step = as.numeric(step)))
  covariance <- fit$estimate$cov
  if (is.null(covariance)) covariance <- fit$stanfit$cov
  se <- sqrt(abs(diag(as.matrix(covariance))))
  delta <- as.numeric(result$delta)
  out$parameters <- data.frame(
    parameter = .ctBackendRawParameterNames(fit, length(est)),
    estimate = est, delta = delta, corrected = est + delta,
    se = se, delta_se = ifelse(se > 0, delta / se, NA_real_),
    stringsAsFactors = FALSE)
  out$corrected <- est + delta
  out$gap_gradient <- as.numeric(result$gap_gradient)
  # How many directions the information matrix did not identify well enough to
  # correct along. Reported rather than silently absorbed: a dropped direction
  # means "no estimable correction here", which is a different statement from
  # "the correction is zero", and the difference matters to anyone reading the
  # table to decide whether the estimate is approximation-limited.
  out$dropped_directions <- if (is.null(result$dropped_directions)) 0L else
    as.integer(result$dropped_directions)

  if (refine) {
    if (verbose > 0) message("Refitting against the quadrature objective")
    refined <- JuliaConnectoR::juliaGet(module$ctsem_laplace_refine(
      objective, .ctJuliaNumericVector(est), nodes = as.integer(nodes),
      maxiter = as.integer(maxiter), verbose = verbose > 1L))
    out$refined <- as.numeric(refined$minimizer)
    out$refined_converged <- isTRUE(refined$converged)
    out$parameters$refined <- out$refined
  }
  out
}

#' @export
print.ctLaplaceCheck <- function(x, ...) {
  cat("Laplace approximation check,", x$nodes, "quadrature nodes per effect\n")
  cat("  log marginal: laplace", format(x$laplace, digits = 8),
    " quadrature", format(x$quadrature, digits = 8), "\n")
  cat("  gap:", format(x$gap, digits = 4), "log units over", x$nsubjects,
    "subjects (", format(x$gap_per_subject, digits = 3), "each )\n")
  if (!is.null(x$parameters)) {
    worst <- x$parameters[order(-abs(x$parameters$delta_se)), ]
    cat("  largest first-order corrections, in standard errors:\n")
    print(utils::head(worst[, c("parameter", "estimate", "corrected", "se",
      "delta_se")], 6), row.names = FALSE, digits = 3)
    if (isTRUE(x$dropped_directions > 0L)) {
      cat("  ", x$dropped_directions, " direction(s) of the information matrix ",
        "were too weakly
  identified to correct along, and are reported as zero.
",
        sep = "")
    }
    if (max(abs(worst$delta_se), na.rm = TRUE) > 0.5) {
      cat("  A correction near or above one standard error means the point ",
        "estimate is\n  limited by the approximation, not by the data.\n", sep = "")
    }
  }
  invisible(x)
}

# Raw parameter labels for the correction table.
#
# The Laplace layout is: the model's own free parameters, then one block per
# level of population scales followed by that level's unconstrained
# correlations, then TI-predictor effects. `.ctJuliaLaplaceSpec` builds it and
# records the indices, so the names are read off those rather than recounted
# here.
.ctBackendRawParameterNames <- function(fit, npar) {
  names <- paste0("raw[", seq_len(npar), "]")
  table <- fit$model_spec$parameter_table
  if (!is.null(table) && nrow(table)) {
    free <- !is.na(table$parnumber) & table$parnumber > 0
    number <- as.integer(table$parnumber[free])
    label <- ifelse(is.na(table$param[free]), paste0("param", number),
      as.character(table$param[free]))
    keep <- !duplicated(number) & number <= npar
    names[number[keep]] <- label[keep]
  }
  laplace <- fit$model_spec$laplace
  if (!is.null(laplace)) {
    for (level in laplace$levels) {
      varying <- as.character(level$param)
      sd_index <- as.integer(level$sd_index)
      if (length(sd_index) && length(varying) == length(sd_index)) {
        names[sd_index] <- paste0("popsd_", varying,
          if (length(laplace$levels) > 1L) paste0(".", level$name) else "")
      }
      cor_index <- as.integer(level$cor_index)
      if (length(cor_index) && length(varying) > 1L) {
        pairs <- which(lower.tri(matrix(0, length(varying), length(varying))),
          arr.ind = TRUE)
        pairs <- pairs[order(pairs[, "col"], pairs[, "row"]), , drop = FALSE]
        n <- min(nrow(pairs), length(cor_index))
        names[cor_index[seq_len(n)]] <- paste0("rawcor_",
          varying[pairs[seq_len(n), "row"]], "_", varying[pairs[seq_len(n), "col"]],
          if (length(laplace$levels) > 1L) paste0(".", level$name) else "")
      }
    }
  }
  names
}
