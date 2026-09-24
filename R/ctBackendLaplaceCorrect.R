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
#' Laplace fits are now corrected when they are fitted, by default
#' (\code{optimcontrol$laplace_correct}, see \code{\link{ctFit}}): the estimate
#' moves and the draws are recentred exactly as \code{draws='normal'} does here.
#' On such a fit (\code{fit$laplace$correction$applied}) this function refuses
#' anything that would apply the step a second time. \code{draws='imis'} is
#' still available there, and redraws around the corrected estimate without
#' moving it. To correct by hand instead, refit with
#' \code{laplace_correct = FALSE}.
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

  # A fit the default correction (optimcontrol$laplace_correct) already moved
  # must not be moved again by the same step. What is still worth asking for
  # there is the shape, so draws='imis' redraws around the corrected estimate
  # and leaves it where it is; anything that would move the estimate or only
  # recentre the draws is refused, because that is what was already done.
  already <- .ctLaplaceIsCorrected(fit)
  if (already) {
    if (!identical(draws, "imis") ||
        (!missing(correct_estimate) && isTRUE(correct_estimate))) {
      stop("This fit was already corrected by quadrature when it was fitted ",
        "(fit$laplace$correction), so ctLaplaceCorrect() would apply the step ",
        "twice. Use draws='imis' to redraw around the corrected estimate, or ",
        "refit with optimcontrol$laplace_correct = FALSE to correct by hand.",
        call. = FALSE)
    }
    correct_estimate <- FALSE
  }

  # The check does the correction arithmetic and the reporting, and it already
  # handles the chunk-count restore, nested groupings, and the near-singular
  # directions it refuses to correct along. Repeating any of that here would be
  # a second implementation of it. On an already-corrected fit only its gap is
  # needed, which is one quadrature evaluation rather than `2 * npar`.
  check <- ctLaplaceCheck(fit, nodes = nodes, correction = !already,
    cores = cores, verbose = verbose)
  if (isTRUE(correct_estimate) && is.null(check$corrected)) {
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
    # `.ctOptimImisDraws()` is shared with `ctParticleCorrect()` and with the
    # uncertainty stage's `uncertainty='is'`: same `imis_is` call, same weighted
    # covariance with an unweighted fallback, same effective-size check. What is
    # this function's own is the density (`quadlp`) and the remedy named below.
    #
    # The proposal is already widened by `scale^2` here, so `scaleInit = 1`
    # rather than compounding two scalings.
    drawn <- .ctOptimImisDraws(quadlp, centre = centre,
      cov = as.matrix(covariance) * scale^2,
      finishsamples = finishsamples, nbatch = nbatch, target_ess = target_ess,
      maxiter = maxiter, scaleInit = 1, tailScale = 1.2, df = Inf,
      verbose = verbose,
      # Said plainly rather than left in a list nobody prints. A corrected
      # interval resting on a handful of effective draws is worse than the
      # uncorrected one, because it looks like it has been improved.
      remedy = paste0("Treat the corrected interval as indicative. ctSample() ",
        "samples the joint posterior directly and does not rely on the ",
        "approximation being close."))
    is_res <- drawn$is_res
    samples <- drawn$samples
    if (is.null(samples) || !nrow(samples)) {
      stop("Importance sampling returned no usable draws against the ",
        "quadrature posterior. The Laplace approximation is likely too far ",
        "from the target to repair by reweighting; use ctSample() instead.",
        call. = FALSE)
    }
    newcov <- drawn$cov
    ess <- drawn$ess
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

# The correction every Laplace fit gets by default ---------------------------
#
# `optimcontrol$laplace_correct` (default TRUE) runs `ctsem_laplace_autocorrect`
# at the end of an optimised `intoverpop='laplace'` fit: after the optimiser,
# after the curvature certification and its Newton continuation have converged
# the Laplace objective, and after the uncertainty stage has built the fit's
# Hessian and draws. It computes no Hessian of its own -- the step's metric is
# the one the standard errors came from, at the point it was evaluated -- and
# pays in tiers, each skipped when the one before finds nothing:
#
#   screen  one quadrature evaluation: the per-unit gap summed in absolute
#           value. Where every unit's likelihood is quadratic in its effects
#           Laplace is exact, the gap is rounding, and the fit is left
#           untouched to the bit.
#   step    `2 * npar` quadrature evaluations for the gradient, then a line
#           search on the quadrature objective. Accepted only on an increase.
#   repeat  up to `maxsteps` (3), same metric, until the predicted gain is
#           small or a step moved no parameter by `step_tol` (0.1) of its
#           standard error. On a 40-subject random-DRIFT fixture one step
#           took 80% of the quadrature objective's gain and three took 99%;
#           see review/LAPLACE-default-correction-2026-09-24.md.
#
# What moves when it applies, and what deliberately does not:
#
#   fit$estimate$raw             the corrected point
#   fit$estimate$loglik          the quadrature log likelihood there (and
#                                `logposterior`, `subject_loglik` with it);
#                                `loglik_laplace` keeps the Laplace value at
#                                the Laplace optimum, `loglik_method` says which
#   fit$estimate$rawposterior    shifted by the step: recentred, not reshaped
#   fit$estimate$cov, $se        unchanged -- the Laplace curvature, at
#                                `fit$uncertainty$evaluated_at`
#   fit$laplace$correction       what was done, and where
.ctLaplaceCorrectDefaults <- list(nodes = 5L, tolerance = 0.01, maxsteps = 3L,
  gain_tol = 1e-3, step_tol = 0.1, step = 1e-3, material = 0.1)

# Whether this fit is corrected, refusing by name a request that cannot apply.
# FALSE is accepted anywhere, since it describes what every other route does.
.ctLaplaceCorrectResolve <- function(optimcontrol, intoverpop, optimize,
  intoverstates) {
  value <- optimcontrol$laplace_correct
  explicit <- !is.null(value)
  if (explicit && !(is.logical(value) && length(value) == 1L && !is.na(value))) {
    stop("optimcontrol$laplace_correct must be TRUE or FALSE.", call. = FALSE)
  }
  if (explicit && isFALSE(value)) return(FALSE)
  why <- if (!identical(as.character(intoverpop)[1L], "laplace")) {
    paste0("applies to intoverpop='laplace' only; with intoverpop='",
      intoverpop, "' there is no Laplace term to correct")
  } else if (!isTRUE(optimize)) {
    paste0("corrects an optimised estimate, and a sampled fit ",
      "(optimize=FALSE) has none to move")
  } else if (!isTRUE(intoverstates)) {
    paste0("corrects the marginal Laplace estimate, which an ",
      "intoverstates=FALSE fit does not make")
  } else if (isTRUE(optimcontrol$estonly)) {
    "steps against the fit's Hessian, which optimcontrol$estonly skips"
  } else NULL
  if (is.null(why)) return(TRUE)
  if (explicit) {
    stop("optimcontrol$laplace_correct ", why, ". Drop it.", call. = FALSE)
  }
  FALSE
}

# The raw estimate the Laplace objective was maximised at: the fit's estimate,
# unless the default correction moved it. Anything that is about the Laplace
# objective itself -- cross-validating it, importance sampling its posterior --
# reads this rather than `fit$estimate$raw`.
.ctLaplaceOptimum <- function(fit) {
  corr <- fit$laplace$correction
  if (isTRUE(corr$applied) &&
      length(corr$laplace_estimate) == length(fit$estimate$raw)) {
    return(as.numeric(corr$laplace_estimate))
  }
  as.numeric(fit$estimate$raw)
}

# Whether the default correction moved this fit's estimate.
.ctLaplaceIsCorrected <- function(fit) isTRUE(fit$laplace$correction$applied)

.ctLaplaceAutoCorrect <- function(fit, cores = 1L, verbose = 0L,
  control = .ctLaplaceCorrectDefaults) {
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  hessian <- fit$uncertainty$hessian
  at <- fit$uncertainty$evaluated_at
  usable <- is.matrix(hessian) && nrow(hessian) == npar &&
    ncol(hessian) == npar && length(at) == npar &&
    isTRUE(all.equal(as.numeric(at), est, tolerance = 0))
  # A 1x1 NaN rather than an empty matrix: an empty one deadlocks the bridge,
  # and the engine reads any size mismatch as "no Hessian".
  if (!usable) hessian <- matrix(NaN, 1L, 1L)
  failed <- function(phrase) {
    warning("Laplace correction skipped: ", phrase, ". The uncorrected fit is ",
      "returned; see fit$laplace$correction.", call. = FALSE)
    fit$laplace$correction <- list(status = "failed", applied = FALSE,
      message = phrase, nodes = as.integer(control$nodes))
    fit
  }
  # The standard errors scale the stopping rule; never an empty vector across
  # the bridge, so an absent one is a vector of NaN, which the engine ignores.
  se <- as.numeric(fit$estimate$se)
  if (length(se) != npar) se <- rep(NaN, npar)
  se[!is.finite(se)] <- NaN
  chunks <- suppressWarnings(as.integer(fit$optim$chunks)[1L])
  if (is.na(chunks) || chunks < 1L) chunks <- max(1L, as.integer(cores)[1L])
  if (verbose > 0) message("Laplace correction: quadrature screen (",
    control$nodes, " nodes)")
  res <- try(.ctBackendWithMaxChunks(chunks, JuliaConnectoR::juliaGet(
    .ctJuliaModule(fit$model_spec$project)$ctsem_laplace_autocorrect(
      .ctJuliaObjective(fit), .ctJuliaNumericVector(est),
      JuliaConnectoR::juliaPut(as.matrix(hessian)),
      nodes = as.integer(control$nodes), step = as.numeric(control$step),
      tolerance = as.numeric(control$tolerance),
      maxsteps = as.integer(control$maxsteps),
      gain_tol = as.numeric(control$gain_tol),
      scale = .ctJuliaNumericVector(se),
      step_tol = as.numeric(control$step_tol)))), silent = TRUE)
  if (inherits(res, "try-error")) {
    return(failed("the quadrature could not be evaluated"))
  }
  status <- as.character(res$status)
  if (identical(status, "quadrature_failed")) {
    return(failed("the quadrature was not finite at the estimate"))
  }
  if (identical(status, "nonfinite_step")) {
    return(failed("the correction step was not finite"))
  }
  steps <- as.integer(res$steps)
  newest <- as.numeric(res$estimate)
  delta <- newest - est
  delta_se <- if (length(se) == npar) ifelse(se > 0, delta / se, NA_real_) else
    rep(NA_real_, npar)
  names(delta) <- names(delta_se) <- names(fit$estimate$se)
  record <- list(status = status, applied = FALSE,
    nodes = as.integer(res$nodes), tolerance = as.numeric(res$tolerance),
    # sum over units of |quadrature - Laplace| at the Laplace optimum, and the
    # signed total the fit's log likelihood was off by there.
    screen = as.numeric(res$screen), gap = as.numeric(res$gap_start),
    laplace_estimate = est,
    loglik_laplace = as.numeric(fit$estimate$loglik),
    steps = steps,
    predicted_gain = if (steps > 0L) as.numeric(res$predicted)[seq_len(steps)] else
      numeric(0),
    alpha = if (steps > 0L) as.numeric(res$alpha)[seq_len(steps)] else numeric(0),
    first_delta = as.numeric(res$first_delta),
    delta = delta, delta_se = delta_se,
    dropped_directions = as.integer(res$dropped_directions),
    # Where the step's metric was evaluated: the Laplace optimum, whose
    # Hessian the standard errors also come from.
    hessian_at = if (usable) "laplace_estimate" else "none",
    seconds = c(screen = as.numeric(res$screen_seconds),
      total = as.numeric(res$seconds)))
  moved <- identical(status, "corrected") && all(is.finite(newest))
  # The log likelihood is the quadrature one wherever the screen found a gap,
  # whether or not a step was accepted. On AnomAuth's spurious optimum no step
  # raised the quadrature objective, and the Laplace value left on the fit was
  # 17 nats above the true optimum's while its quadrature value was 10 below:
  # the reported number is what makes such a fit comparable, so it cannot wait
  # on the estimate moving.
  units <- as.numeric(res$quadrature_units)
  quadrature_ok <- status %in% c("corrected", "no_gain", "no_hessian") &&
    is.finite(as.numeric(res$quadrature)) && all(is.finite(units))
  if (quadrature_ok) {
    record$loglik_quadrature <- sum(units)
    record$logposterior_quadrature <- as.numeric(res$quadrature)
    record$logposterior_laplace <- as.numeric(res$laplace)
    # Quadrature minus Laplace at the point the fit now reports.
    record$gap_reported <- as.numeric(res$quadrature) - as.numeric(res$laplace)
    fit$estimate$loglik_laplace <- record$loglik_laplace
    fit$estimate$loglik <- record$loglik_quadrature
    fit$estimate$logposterior <- record$logposterior_quadrature
    subjects <- as.numeric(res$quadrature_subjects)
    if (length(subjects) == length(fit$estimate$subject_loglik) &&
        all(is.finite(subjects))) {
      fit$estimate$subject_loglik <- subjects
    }
    fit$estimate$loglik_method <- "quadrature"
  }
  if (!moved) {
    fit$laplace$correction <- record
    return(fit)
  }
  # Applied. The draws are moved by the step and not redrawn, so they keep the
  # Laplace curvature's shape -- recentred, not reshaped, which is what
  # `draws='imis'` in ctLaplaceCorrect() is for.
  record$applied <- TRUE
  record$material <- isTRUE(max(abs(delta_se), na.rm = TRUE) >= control$material)
  record$draws <- "recentred"
  fit$estimate$raw <- newest
  post <- fit$estimate$rawposterior
  if (!is.null(post) && ncol(post) == npar) {
    fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
  }
  if (!is.null(fit$uncertainty)) {
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("recentred on the corrected estimate; the covariance is ",
        "the Laplace curvature at fit$uncertainty$evaluated_at"),
      nodes = record$nodes)
  }
  fit$laplace$correction <- record
  fit
}
