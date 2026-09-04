# Apply the Laplace correction, and draw from the corrected posterior.
#
# `ctLaplaceCheck()` measures the approximation and reports what it would cost
# to ignore it. This applies that finding: the estimate moves to the corrected
# point, and the draws behind every interval the fit reports are drawn against
# the quadrature posterior rather than the Laplace one.
#
# The two halves matter for different reasons. Correcting the estimate alone
# moves the centre and leaves the width -- and the width is wrong in the same
# direction, because the Laplace profile is tilted rather than merely shifted.
# Re-drawing alone would sample a target the proposal is offset from by exactly
# the bias being corrected, which is the situation importance sampling handles
# worst. Doing both is what makes the interval mean what it says.
#
# The target is available because `ctsem_laplace_quadrature` evaluates the same
# adaptive Gauss-Hermite rule at an arbitrary parameter vector, and includes the
# prior term, so it is a log posterior on the same scale as the Laplace
# objective the fit maximised. Importance sampling needs a density and not a
# gradient, which is why this costs one quadrature evaluation per draw rather
# than the `2 * npar` a gradient would.

#' Correct a Laplace fit's estimate and posterior draws by quadrature
#'
#' Applies the first-order correction \code{\link{ctLaplaceCheck}} measures, so
#' that the returned fit reports the corrected estimate and every downstream
#' summary follows from it. Optionally also redraws the posterior sample by
#' importance sampling against the adaptive Gauss-Hermite posterior rather than
#' the Laplace approximation to it.
#'
#' The two are separate because they cost very different amounts. Correcting the
#' estimate is one Newton step and costs \code{2 * npar} quadrature evaluations
#' whatever the model size. Importance sampling costs one quadrature evaluation
#' \emph{per draw}, which for a normal request is three orders of magnitude
#' more. Nobody should have to pay the second price to get the first, so
#' \code{draws='normal'} is the default: the estimate is corrected, and the
#' draws are recentred on it carrying the curvature the fit already had.
#'
#' \code{intoverpop='laplace'} replaces each unit's integral over its random
#' effects by a single Gaussian fitted at the mode. Where that is inexact the
#' error is not constant -- it grows with the population scale, so it tilts the
#' profile rather than shifting it, and both the estimate and the width of the
#' interval around it are affected. On a simulated model with a random DRIFT the
#' population standard deviation's correction is about 0.9 standard errors,
#' which is the whole of that parameter's under-coverage.
#'
#' The estimate moves by one Newton step against the curvature already computed
#' for the standard errors. What that fixes is the \emph{location}: with
#' \code{draws='normal'} the interval is the fit's own width, recentred. With
#' \code{draws='imis'} the draws come instead from incremental mixture
#' importance sampling against the quadrature posterior, which corrects the shape
#' as well, with the proposal centred at the \emph{corrected} point -- a proposal
#' centred at the uncorrected estimate is displaced from its target by the bias
#' being corrected, and a displaced proposal is what importance sampling handles
#' least well.
#'
#' Under \code{draws='imis'}, read \code{ess} on the result. Importance sampling
#' reports its own reliability, and a low effective sample size means the
#' corrected draws are carried by a few points and the interval should not be
#' trusted. When that happens the approximation is too far from the target to
#' repair by reweighting, and \code{\link{ctSample}} is the answer rather than
#' this.
#'
#' @param fit A \code{ctJuliaFit} fitted with \code{intoverpop='laplace'}, with
#'   uncertainty already computed -- the correction needs the Hessian and the
#'   proposal needs the covariance.
#' @param draws How to obtain the posterior draws behind the fit's intervals.
#'   \code{'normal'} (the default) recentres normal draws on the corrected
#'   estimate, keeping the fit's covariance: cheap, and it corrects the location
#'   but not the shape. \code{'imis'} importance samples against the quadrature
#'   posterior, correcting both, at one quadrature evaluation per draw.
#'   \code{'keep'} corrects the estimate and leaves the existing draws alone,
#'   which is only sensible when the draws are about to be replaced anyway.
#' @param nodes Quadrature nodes per random effect, as in
#'   \code{\link{ctLaplaceCheck}}. Cost is \code{nodes^k} process log
#'   likelihoods per block for every evaluation, so this is the main cost
#'   control.
#' @param finishsamples Draws to return. Defaults to as many as the fit already
#'   carries.
#' @param scale \code{draws='imis'} only. Proposal scale multiplier on the
#'   fit's covariance. Above one by default: the Hessian covariance is typically
#'   \emph{narrower} than the posterior, and a proposal narrower than its target
#'   cannot correct it.
#' @param target_ess \code{draws='imis'} only. Effective sample size at which
#'   sampling stops.
#' @param nbatch \code{draws='imis'} only. Draws per importance-sampling
#'   iteration.
#' @param maxiter \code{draws='imis'} only. Iteration cap.
#' @param correct_estimate Move the point estimate to the corrected point. When
#'   \code{FALSE} only the draws are redrawn, which is occasionally what you
#'   want when comparing.
#' @param cores Engine threads for the quadrature.
#' @param verbose Integer; 1 or more prints progress.
#'
#' @return The fit, with \code{estimate$raw} at the corrected point,
#'   \code{estimate$rawposterior} and \code{transformedpars} redrawn around it
#'   (and \code{estimate$cov} and \code{estimate$se} replaced too under
#'   \code{draws='imis'}), and a \code{laplace_correction} entry recording what
#'   was done: the \code{ctLaplaceCheck} it was based on, which draws were used,
#'   the effective sample size where one applies, and the estimate before
#'   correction.
#'
#' @seealso \code{\link{ctLaplaceCheck}} measures without applying.
#'   \code{\link{ctSample}} removes the approximation instead of correcting it.
#'
#' @examples
#' \dontrun{
#' fit <- ctFit(data, model, backend = 'julia', intoverpop = 'laplace')
#'
#' corrected <- ctLaplaceCorrect(fit)              # estimate corrected, cheap
#' corrected$laplace_correction                    # gap, corrections in se
#' summary(corrected)                              # intervals about the corrected point
#'
#' # Correct the shape of the posterior too, at one quadrature evaluation per draw
#' resampled <- ctLaplaceCorrect(fit, draws = 'imis')
#' resampled$laplace_correction$ess                # check before trusting it
#' }
#' @export
ctLaplaceCorrect <- function(fit, draws = c("normal", "imis", "keep"),
  nodes = 5L, finishsamples = NULL,
  scale = 1.5, target_ess = 100, nbatch = NULL, maxiter = 10L,
  correct_estimate = TRUE, cores = NULL, verbose = 0L) {

  draws <- match.arg(draws)

  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctLaplaceCorrect applies to backend='julia' fits.", call. = FALSE)
  }
  if (is.null(fit$model_spec$laplace)) {
    stop("ctLaplaceCorrect applies to fits made with intoverpop='laplace'. ",
      "The augmented route carries the random effects in the state and has no ",
      "separate integral to correct.", call. = FALSE)
  }
  covariance <- fit$estimate$cov
  if (is.null(covariance) || !all(is.finite(as.matrix(covariance)))) {
    stop("ctLaplaceCorrect needs the fit's parameter covariance: the ",
      "correction is a Newton step against the curvature computed for the ",
      "standard errors. Run ctOptimUncertainty(fit) first, or refit without ",
      "optimcontrol$estonly.", call. = FALSE)
  }
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)

  # The check does the correction arithmetic and the reporting, and it already
  # handles the chunk-count restore, nested groupings, and the near-singular
  # directions it refuses to correct along. Repeating any of that here would be
  # a second implementation of it.
  check <- ctLaplaceCheck(fit, nodes = nodes, correction = TRUE,
    cores = cores, verbose = verbose)
  if (is.null(check$corrected)) {
    stop("The correction could not be formed; see the warning from ",
      "ctLaplaceCheck.", call. = FALSE)
  }
  centre <- if (isTRUE(correct_estimate)) as.numeric(check$corrected) else est

  module <- .ctJuliaModule(fit$model_spec$project)
  objective <- .ctJuliaObjective(fit)
  # The quadrature posterior, as a scalar log density. `imis_is` asks for
  # nothing else -- importance weights need the target's density and not its
  # slope, which is what keeps this to one quadrature evaluation per draw.
  #
  # An invalid point returns a large finite penalty rather than aborting, as
  # the Laplace path's own `lpgFunc` does. Here it is also the right answer
  # statistically: a proposal draw that lands where the model cannot be
  # evaluated should carry no weight, and `exp(-1e100 - c)` is zero.
  quadlp <- function(parm) {
    value <- try(.ctBackendJuliaValue(module$ctsem_laplace_quadrature(
      objective, .ctJuliaNumericVector(as.numeric(parm)),
      nodes = as.integer(nodes))$value), silent = TRUE)
    value <- if (inherits(value, "try-error")) NaN else as.numeric(value)[1L]
    if (!is.finite(value)) value <- -1e100
    value
  }

  if (is.null(finishsamples)) {
    finishsamples <- if (!is.null(fit$estimate$rawposterior))
      nrow(fit$estimate$rawposterior) else 1000L
  }
  if (is.null(nbatch)) nbatch <- max(200L, as.integer(finishsamples))

  # Correcting the estimate and resampling the posterior are separate jobs with
  # very different prices, so they are separately chosen. The correction is one
  # Newton step and costs `2 * npar` quadrature evaluations however many draws
  # are involved. Importance sampling costs one quadrature evaluation *per
  # draw*, which is three orders of magnitude more on a typical request --
  # worth it when the shape of the posterior matters, and not something to make
  # anyone pay for a corrected point estimate.
  is_res <- NULL
  ess <- NA_real_
  newcov <- as.matrix(covariance)
  samples <- fit$estimate$rawposterior

  if (identical(draws, "imis")) {
    if (verbose > 0) message("Importance sampling against the quadrature ",
      "posterior (", nodes, " nodes, ", nbatch, " draws per iteration)")
    is_res <- imis_is(quadlp, mu_hat = centre,
      Sigma_hat = as.matrix(covariance) * scale^2,
      cl = NA, n_batch = as.integer(nbatch), target_ess = target_ess,
      max_iter = as.integer(maxiter),
      # The proposal is already scaled above, so the sampler's own initial
      # scaling is left at one rather than compounding with it.
      scale_init = 1, tail_scale = 1.2, df = Inf,
      finishsamples = as.integer(finishsamples), verbose = verbose > 0)
    samples <- is_res$theta
    if (is.null(samples) || !nrow(samples)) {
      stop("Importance sampling returned no usable draws against the ",
        "quadrature posterior. The Laplace approximation is likely too far ",
        "from the target to repair by reweighting; use ctSample() instead.",
        call. = FALSE)
    }
    newcov <- if (!is.null(is_res$covariance) && all(is.finite(is_res$covariance)))
      ctOptimSafeCov(is_res$covariance) else ctOptimSafeCov(stats::cov(samples))
    ess <- if (is.null(is_res$ess)) NA_real_ else as.numeric(is_res$ess)[1L]
    # Said plainly rather than left in a list nobody prints. A corrected
    # interval resting on a handful of effective draws is worse than the
    # uncorrected one, because it looks like it has been improved.
    if (is.finite(ess) && ess < target_ess / 2) {
      warning("Importance sampling reached an effective sample size of ",
        round(ess, 1), " against a target of ", target_ess,
        ". The corrected draws rest on few points, so treat the corrected ",
        "interval as indicative. ctSample() samples the joint posterior ",
        "directly and does not rely on the approximation being close.",
        call. = FALSE)
    }
  } else if (identical(draws, "normal")) {
    # The centre is corrected and the width is not: these draws carry the
    # Laplace curvature, moved to the corrected point. That is the honest
    # description of what one Newton step buys, and it is the right default --
    # the location error is what `ctLaplaceCheck` measures and what tilts an
    # estimate, while correcting the *shape* costs a quadrature evaluation per
    # draw and is asked for with `draws='imis'`.
    samples <- ctOptimNormalDraws(centre, newcov, finishsamples)
  }

  before <- list(raw = est, cov = fit$estimate$cov, se = fit$estimate$se)
  if (isTRUE(correct_estimate)) fit$estimate$raw <- centre
  if (!identical(draws, "keep")) {
    fit$estimate$cov <- newcov
    fit$estimate$se <- sqrt(diag(newcov))
    fit$estimate$rawposterior <- samples
    fit <- .ctFitNameRawUncertainty(fit)
    # The constrained draws describe whatever raw draws they were built from,
    # so they are refreshed here rather than left to disagree with the ones
    # above.
    fit$transformedpars <- .ctBackendConstrain(fit)
  } else if (isTRUE(correct_estimate)) {
    # The estimate moved and the draws did not, so they now describe a point
    # the fit no longer reports. Saying so is better than leaving it to be
    # discovered in a summary.
    warning("The estimate was corrected but draws='keep' left the posterior ",
      "draws as they were, so intervals still describe the uncorrected point. ",
      "Use draws='normal' to recentre them.", call. = FALSE)
  }
  fit$laplace_correction <- list(
    check = check, nodes = as.integer(nodes), ess = ess, draws = draws,
    effective_fraction = if (is.finite(ess)) ess / nrow(samples) else NA_real_,
    corrected_estimate = isTRUE(correct_estimate),
    proposal_scale = scale, before = before,
    largest_delta_se = if (is.null(check$parameters)) NA_real_ else
      max(abs(check$parameters$delta_se), na.rm = TRUE))
  class(fit$laplace_correction) <- "ctLaplaceCorrection"
  if (!is.null(fit$uncertainty) && !identical(draws, "keep")) {
    fit$uncertainty$cov <- newcov
    # The draws on the fit are now importance-sampled, so the record of how
    # they were made has to say so. Downstream code reads these to describe the
    # fit, and leaving them at whatever the original uncertainty pass wrote
    # would have every summary claim normal draws around a Hessian covariance.
    fit$uncertainty$draws <- draws
    if (!is.null(fit$uncertainty$settings)) {
      fit$uncertainty$settings$draws <- draws
      fit$uncertainty$settings$finishsamples <- nrow(samples)
    }
    if (!is.null(is_res)) {
      fit$uncertainty$proposal_cov <- as.matrix(covariance) * scale^2
      fit$uncertainty$imis <- is_res
      fit$uncertainty$details$importance_sampling <- list(ess = ess,
        df_used = is_res$df_used,
        covariance = "weighted importance-sampling covariance")
    }
    fit$uncertainty$details$laplace_correction <- list(
      nodes = as.integer(nodes), ess = ess,
      target = "adaptive Gauss-Hermite posterior")
  }
  fit
}

#' @export
print.ctLaplaceCorrection <- function(x, ...) {
  cat("Laplace correction applied\n")
  cat("  quadrature nodes      ", x$nodes, "\n", sep = "")
  cat("  log marginal gap      ", format(x$check$gap, digits = 4),
    " (", format(x$check$gap_per_subject, digits = 3), " per subject)\n", sep = "")
  cat("  largest correction    ", format(x$largest_delta_se, digits = 3),
    " standard errors\n", sep = "")
  cat("  draws                 ", x$draws,
    if (identical(x$draws, "normal"))
      "  (recentred; width is the fit's own curvature)" else "", "\n", sep = "")
  if (is.finite(x$ess)) {
    cat("  effective sample size ", format(x$ess, digits = 4),
      if (is.finite(x$effective_fraction))
        paste0(" (", format(100 * x$effective_fraction, digits = 3), "% of draws)")
      else "", "\n", sep = "")
  }
  if (!x$corrected_estimate) {
    cat("  estimate left at the uncorrected point; draws redrawn only\n")
  }
  invisible(x)
}
