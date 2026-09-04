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
#'
#'   \code{dropped_directions} counts the directions of the information matrix
#'   too weakly identified to correct along, which are reported as zero rather
#'   than as the arbitrarily large step an unguarded solve would produce. On a
#'   model where the Laplace approximation is already exact the gap gradient is
#'   floating-point noise, and dividing that by a near-singular curvature is how
#'   a meaningless correction gets a plausible-looking number.
#'
#' @seealso \code{\link{ctSample}} removes the approximation instead of
#'   measuring it, by sampling the joint posterior; much slower, and the right
#'   answer when \code{delta_se} says a linear correction is not credible.
#'
#' @examples
#' \dontrun{
#' data <- ctstantestdat
#' model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
#'   LAMBDA = matrix(1))
#' model$pars$indvarying <- model$pars$matrix %in% 'MANIFESTMEANS'
#' fit <- ctFit(data, model, backend = 'julia', intoverpop = 'laplace')
#'
#' check <- ctLaplaceCheck(fit)
#' check                     # gap, and the largest corrections in standard errors
#' check$parameters          # the full per-parameter table
#' }
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
    previous <- .ctBackendSetMaxChunks(max(1L, cores))
    on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  }

  nsubjects <- length(fit$model_spec$subject_starts)

  # Where `ctOptimUncertainty()` leaves it. A `fit$stanfit$uncertainty$hessian`
  # fallback used to sit here "so this reads either without the caller knowing
  # which" -- but the caller is known: this function refuses anything that is
  # not a ctJuliaFit thirty lines above, and nothing ever puts a `$stanfit` on
  # one, so the fallback could not fire.
  hessian <- fit$uncertainty$hessian
  if (isTRUE(correction) && is.null(hessian)) {
    warning("No Hessian on the fit, so the correction cannot be formed. ",
      "Run ctOptimUncertainty(fit) first, or use correction=FALSE.",
      call. = FALSE)
  }
  do_correction <- isTRUE(correction) && !is.null(hessian)

  if (do_correction) {
    # `ctsem_laplace_correction` makes this same quadrature and Laplace
    # evaluation itself and returns both, along with their difference -- so
    # take `quadrature`, `laplace` and `gap` from its result rather than
    # evaluating both a second time. The `correction=FALSE` path below still
    # needs to make these two calls itself.
    if (verbose > 0) message("Gap gradient (", 2 * length(est), " quadrature evaluations)")
    result <- JuliaConnectoR::juliaGet(module$ctsem_laplace_correction(
      objective, .ctJuliaNumericVector(est),
      JuliaConnectoR::juliaPut(as.matrix(hessian)),
      nodes = as.integer(nodes), step = as.numeric(step)))
    quadrature <- as.numeric(result$quadrature)
    laplacevalue <- as.numeric(result$laplace)
    gap <- as.numeric(result$gap)
  } else {
    if (verbose > 0) message("Quadrature at the estimate (", nodes, " nodes)")
    quadrature <- .ctBackendJuliaValue(module$ctsem_laplace_quadrature(
      objective, .ctJuliaNumericVector(est), nodes = as.integer(nodes))$value)
    laplacevalue <- .ctBackendJuliaValue(module$ctsem_laplace_evaluate(
      objective, .ctJuliaNumericVector(est), gradient = FALSE)$value)
    gap <- quadrature - laplacevalue
  }

  out <- list(gap = gap, quadrature = quadrature, laplace = laplacevalue,
    gap_per_subject = gap / max(1L, nsubjects),
    nodes = as.integer(nodes), nsubjects = nsubjects)
  class(out) <- "ctLaplaceCheck"

  if (!do_correction) return(out)
  covariance <- fit$estimate$cov
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

# A name for every element of the raw (unconstrained) parameter vector.
#
# These names are what a user is shown when something has to point at one raw
# coordinate: the identifiability warning and `fit$identifiability`,
# `ctIdentify()`, the Laplace correction table, `ctReport()`'s profile
# component, the column names of `fit$estimate$rawposterior`. So the whole
# vector has to be namable, not most of it: a positional `raw[17]` in that
# list says nothing about which part of the model is involved, which is the
# one thing the reader needs.
#
# `.ctJuliaCheckLayout()` (R/ctJuliaBackend.R) enumerates the blocks the raw
# vector is made of, and they tile `1:npar`. Every block therefore has a name
# here, spelled to match what the same quantity is called elsewhere:
#
#   model parameters              the cell's own name         drift_eta2_eta1
#   population scales             popsd_<parameter>           popsd_drift_eta1
#   population correlations       rawcor_<par>__<par>         rawcor_a__b
#   TI-predictor coefficients     rawtipredeffect_<par>_<ti>  as stan's
#                                                             ctFitgetparnamesfromraw()
#   sampled TI-predictor values   tipredvalue_<ti>_<subject>
#
# Two of those need the spec rather than the parameter table. The augmented
# route (`intoverpop='augmented'`) keeps its population scales and correlations
# *in* the table, under internal cell names (`julia_popcov_3_3`), so they are
# relabelled from `random_effects` to the same `popsd_`/`rawcor_` spelling the
# Laplace route uses -- one conceptual parameter, one name, whichever route
# produced it. TI-predictor coefficients and sampled TI-predictor values are
# not matrix cells at all and appear only on the spec.
#
# A leftover `raw[i]` is therefore a block nobody named rather than a name, and
# is reported as such: see `.ctBackendReportUnnamedParameters()` below.
.ctBackendRawParameterNames <- function(fit, npar) {
  names <- paste0("raw[", seq_len(npar), "]")
  spec <- fit$model_spec
  table <- spec$parameter_table
  if (!is.null(table) && nrow(table)) {
    free <- !is.na(table$parnumber) & table$parnumber > 0
    number <- as.integer(table$parnumber[free])
    label <- .ctBackendParamLabel(table$param[free], number)
    keep <- !duplicated(number) & number <= npar
    names[number[keep]] <- label[keep]
  }

  # The augmented route's population scales and correlations, which the loop
  # above just labelled with their internal cell names.
  effects <- spec$random_effects
  if (!is.null(effects) && length(effects) && nrow(effects)) {
    index <- as.integer(effects$parameter)
    label <- ifelse(effects$type %in% "correlation",
      paste0("rawcor_", as.character(effects$param)),
      paste0("popsd_", as.character(effects$param)))
    keep <- !is.na(index) & index >= 1L & index <= npar & !is.na(effects$param)
    names[index[keep]] <- label[keep]
  }

  laplace <- spec$laplace
  if (!is.null(laplace)) {
    for (level in laplace$levels) {
      varying <- as.character(level$param)
      suffix <- if (length(laplace$levels) > 1L) paste0(".", level$name) else ""
      sd_index <- as.integer(level$sd_index)
      if (length(sd_index) && length(varying) == length(sd_index)) {
        names[sd_index] <- paste0("popsd_", varying, suffix)
      }
      cor_index <- as.integer(level$cor_index)
      if (length(cor_index) && length(varying) > 1L) {
        pairs <- which(lower.tri(matrix(0, length(varying), length(varying))),
          arr.ind = TRUE)
        pairs <- pairs[order(pairs[, "col"], pairs[, "row"]), , drop = FALSE]
        n <- min(nrow(pairs), length(cor_index))
        # Two underscores between the pair, as `summary()`'s `rawpopcorr`
        # column names use: parameter names contain single underscores
        # themselves (`drift_eta2_eta1`), so a single one does not say where
        # the first name stops.
        names[cor_index[seq_len(n)]] <- paste0("rawcor_",
          varying[pairs[seq_len(n), "row"]], "__", varying[pairs[seq_len(n), "col"]],
          suffix)
      }
    }
  }

  # TI-predictor coefficients, named after the parameter they act on and the
  # predictor they carry -- the same pair, in the same order, that the stan
  # path's `ctFitgetparnamesfromraw()` spells `rawtipredeffect_<par>_<ti>`.
  # Named last of the three parameter blocks because the target parameter's own
  # name has to be resolved first.
  tipred_names <- spec$TIpredNames
  if (is.null(tipred_names)) tipred_names <- spec$model$TIpredNames
  ti_effects <- spec$ti_effects
  if (!is.null(ti_effects) && length(ti_effects) && nrow(ti_effects)) {
    index <- as.integer(ti_effects$coefficient)
    target <- as.integer(ti_effects$parameter)
    predictor <- .ctBackendPredictorLabel(tipred_names, ti_effects$predictor)
    keep <- !is.na(index) & index >= 1L & index <= npar &
      !is.na(target) & target >= 1L & target <= npar
    names[index[keep]] <- paste0("rawtipredeffect_", names[target[keep]], "_",
      predictor[keep])
  }

  # Sampled values for missing TI-predictor cells: one parameter per cell, so
  # the name says which subject's which predictor. The subject label is the
  # user's own id where the spec still carries the data (`ti_missing$subject`
  # indexes subjects in first-appearance order, as `.ctJuliaTIData()` builds
  # them), and the position otherwise.
  ti_missing <- spec$ti_missing
  if (!is.null(ti_missing) && length(ti_missing) && nrow(ti_missing)) {
    index <- as.integer(ti_missing$parameter)
    predictor <- .ctBackendPredictorLabel(tipred_names, ti_missing$predictor)
    subject <- .ctBackendSubjectLabel(fit, ti_missing$subject)
    keep <- !is.na(index) & index >= 1L & index <= npar
    names[index[keep]] <- paste0("tipredvalue_", predictor[keep], "_",
      subject[keep])
  }

  .ctBackendReportUnnamedParameters(names, spec)
  names
}

# A predictor's own name, or its position when the spec carries no name for it.
.ctBackendPredictorLabel <- function(tipred_names, predictor) {
  predictor <- as.integer(predictor)
  label <- if (is.null(tipred_names)) rep(NA_character_, length(predictor)) else
    as.character(tipred_names)[predictor]
  ifelse(is.na(label) | !nzchar(label), paste0("TI", predictor), label)
}

# The subject id a raw-vector index belongs to, as the user wrote it.
#
# `.ctFitIdMap()` (R/ctBackendKalman.R) is the position-to-original-id mapping
# every other user-facing report speaks through, so the name here says the same
# id the rest of the output does. It needs the long data frame, which a fit
# carries on its spec and a *prepared* spec does not -- `ctFit()` overwrites
# the top-level `$data` of what it returns with the prepared standata -- so the
# position is the fallback, spelled so it cannot be mistaken for an id.
.ctBackendSubjectLabel <- function(fit, subject) {
  subject <- as.integer(subject)
  map <- try(.ctFitIdMap(fit), silent = TRUE)
  label <- if (inherits(map, "try-error") || !is.data.frame(map) || !nrow(map)) {
    rep(NA_character_, length(subject))
  } else as.character(map$original)[match(subject, map$new)]
  ifelse(is.na(label), paste0("subject", subject), label)
}

# Say so when a raw index had no name to give.
#
# `raw[17]` is a defensible last resort, but it is not a name, and shipping it
# inside a list of names reads as one. The blocks above cover every block
# `.ctJuliaCheckLayout()` knows about, so a gap here means the raw vector grew
# a block that nobody taught this function to name -- a developer's problem,
# and one that otherwise reaches a user as an unreadable identifiability
# warning. Reported only when the spec is complete enough to have been
# namable: a fit restored without its parameter table has no gap, it has no
# spec.
.ctBackendReportUnnamedParameters <- function(names, spec) {
  table <- spec$parameter_table
  if (is.null(table) || !nrow(table)) return(invisible(NULL))
  unnamed <- which(startsWith(names, "raw["))
  if (!length(unnamed)) return(invisible(NULL))
  warning("Raw parameter ", paste(unnamed, collapse = ", "),
    " could not be named, and is reported by position. This is a gap in ",
    "ctsem's parameter naming rather than a problem with the model; please ",
    "report it.", call. = FALSE)
  invisible(NULL)
}
