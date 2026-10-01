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
#' repair by reweighting, and \code{\link{ctFitUncertainty}} with
#' \code{uncertainty = 'sample'} is the answer rather than this.
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
#'   sampling stops; 200 by default, the target every draw-producing route
#'   shares (see \code{\link{ctFitUncertainty}}). Ending short of it, at
#'   \code{maxiter}, is warned about.
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
#'   \code{\link{ctFitUncertainty}} (\code{uncertainty = 'sample'}) removes the
#'   approximation instead of correcting it.
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
  scale = 1.5, target_ess = 200, nbatch = NULL, maxiter = 10L,
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

  # A fit the fit-time correction (optimcontrol$laplace_correct, either
  # method) already moved must not be moved again: this would be a second
  # correction of the same error. What is still worth asking for
  # there is the shape, so draws='imis' redraws around the corrected estimate
  # and leaves it where it is; anything that would move the estimate or only
  # recentre the draws is refused, because that is what was already done.
  already <- .ctLaplaceIsCorrected(fit)
  if (already) {
    if (!identical(draws, "imis") ||
        (!missing(correct_estimate) && isTRUE(correct_estimate))) {
      stop("This fit was already corrected by quadrature when it was fitted ",
        "(fit$laplace$correction), so ctLaplaceCorrect() would correct it ",
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
  paretok <- NA_real_
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
      remedy = paste0("Treat the corrected interval as indicative. ",
        "ctFitUncertainty(fit, 'sample') samples the joint posterior directly ",
        "and does not rely on the approximation being close."))
    is_res <- drawn$is_res
    samples <- drawn$samples
    if (is.null(samples) || !nrow(samples)) {
      stop("Importance sampling returned no usable draws against the ",
        "quadrature posterior. The Laplace approximation is likely too far ",
        "from the target to repair by reweighting; use ",
        "ctFitUncertainty(fit, 'sample') instead.",
        call. = FALSE)
    }
    newcov <- drawn$cov
    ess <- drawn$ess
    paretok <- drawn$k
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
        pareto_k = paretok, df_used = is_res$df_used,
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
#
# `'quadrature'`, the default, runs a different stage in the same place: see
# `.ctLaplaceContinue()` below.
#
# Provenance of the constants (Appendix B of
# review/OPTIM-consolidation-plan-2026-09-25.md):
#
#   nodes      5   the quadrature's own default since ctLaplaceCheck: enough to
#                  locate the maximum in the cases tested, 9 to settle the
#                  value (quadrature.jl's docstring measured the gap on one
#                  40-subject dataset).
#   tolerance  0.01 nats, the screen's pass mark, summed absolutely over units.
#                  A tolerance in objective units rather than a measured value:
#                  the linear fixture screens at 3e-14, so any model with a
#                  nonlinear effect clears it by orders of magnitude.
#   maxsteps   3   measured, on the 40-subject nonlinear fixture only: three
#                  steps take 99% of the quadrature objective's gain and stop
#                  0.1 se short of the refine optimum
#                  (review/LAPLACE-default-correction-2026-09-24.md 3a).
#   gain_tol   1e-3 nats predicted, and
#   step_tol   0.1 se: the stopping rules, chosen so the next gradient is not
#                  paid for; set on the same fixture.
#   step       1e-3 raw, the finite-difference step of the gap gradient. Not
#                  tuned; differences of the gap, which is smooth and small.
#   material   0.1 se: when print() says the estimate moved. Presentation.
.ctLaplaceCorrectDefaults <- list(nodes = 5L, tolerance = 0.01, maxsteps = 3L,
  gain_tol = 1e-3, step_tol = 0.1, step = 1e-3, material = 0.1)

# Which correction this fit gets -- FALSE, "step" or "quadrature" -- refusing
# by name a request that cannot apply. FALSE is accepted anywhere, since it
# describes what every other route does. TRUE, and the default, are
# "quadrature", the continuation: on the known-broken cases of
# review/OPTIM-consolidation-plan-2026-09-25.md section 10 it ended higher in
# exact log likelihood than the step correction on 9 of 11 and level on the
# other 2, and never below the Laplace optimum, where the step correction
# ended 1.4 nats below it on the gated-gaps A14 config.
.ctLaplaceCorrectResolve <- function(optimcontrol, intoverpop, optimize,
  intoverstates) {
  value <- optimcontrol$laplace_correct
  explicit <- !is.null(value)
  valid <- (is.logical(value) && length(value) == 1L && !is.na(value)) ||
    (is.character(value) && length(value) == 1L &&
      value %in% c("step", "quadrature"))
  if (explicit && !valid) {
    stop("optimcontrol$laplace_correct must be TRUE or FALSE, or 'step' or ",
      "'quadrature'.", call. = FALSE)
  }
  if (explicit && isFALSE(value)) return(FALSE)
  method <- if (is.character(value)) value else "quadrature"
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
  if (is.null(why)) return(method)
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
  # The curvature at the estimate: evaluated there, or within a hundredth of a
  # standard error, where the optimiser's finish keeps the Hessian it took at
  # the hand-over (see `.ctBackendStoredHessian()`).
  usable <- !is.null(.ctBackendStoredHessian(fit, est))
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
  res <- try(.ctBackendWithMaxChunks(chunks, .ctJuliaGet(
    .ctJuliaModule(fit$model_spec$project)$ctsem_laplace_autocorrect(
      .ctJuliaObjective(fit), .ctJuliaNumericVector(est),
      .ctJuliaPut(as.matrix(hessian)),
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
  record <- list(method = "step", status = status, applied = FALSE,
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

# The quadrature continuation (`optimcontrol$laplace_correct = 'quadrature'`,
# the default) --
#
# Where `'step'` takes up to three Newton steps on the quadrature objective by
# finite differences, `'quadrature'` climbs it: from the Laplace optimum, with the
# engine's own optimiser, on an objective whose nodes are held fixed so that
# its gradient is exact (the Fisher identity; see laplace_continuation.jl), in
# rounds between which the nodes are re-placed. bigIRT's `laplaceRefine` is the
# model (../bigIRT/optimisation-review-2026-09.md, R3.4 and R3.5); what differs
# here, and why, is in review/OPTIM-consolidation-plan-2026-09-25.md, section
# 10.
#
#   screen   the rule at the Laplace optimum against Laplace, unit by unit, as
#            `'step'` screens. A fit that passes is left untouched to the bit.
#   flags    only the units the screen flags carry nodes: largest gap first,
#            until what is left carries at most the screen's tolerance. Every
#            other unit keeps its Laplace term and its exact gradient.
#   rounds   each a trust region in the Laplace fit's own standard errors:
#            the round's L-BFGS runs in the fit's identified directions,
#            whitened by its curvature, so it starts from the Newton step the
#            fit's Hessian implies and is confined to `radius` standard errors
#            of the round's centre; then the nodes are re-placed where it
#            ended.
#   target   the FIXED POINT of that iteration: the estimate at which the
#            gradient of the objective with its nodes placed there is zero.
#            With the nodes held, that gradient is the quadrature's own
#            estimate of the exact score -- a node-weighted average of the
#            per-node scores -- so the fixed point solves the quadrature
#            version of the exact likelihood equation. It is not the maximum
#            of the quadrature VALUE with the nodes re-placed at every point
#            (what the finite-difference refinement found), whose gradient
#            also carries how the rule's error moves with its centre. The two
#            differ by that term, and on the 40-subject fixture the fixed point
#            is the nearer to the exact optimum: 0.07 of a standard error
#            against 0.14 (the exact marginal by a dense grid; see
#            test-julia-laplace-continue.R). bigIRT's refinement targets the
#            same fixed point. It is not always the better point, which is
#            why the rounds only approach it while the re-placed value holds
#            (`accept`): where the two part, the estimate is the last point
#            that did not lose value, status `stalled`, and the fit's
#            certification says it is not stationary.
#   accept   a round is kept when the fixed-point residual falls -- half the
#            squared whitened gradient at the re-placed point, the Newton gain
#            it predicts -- AND the quadrature value with the nodes re-placed
#            there falls by no more than `value_tol`. Either failing, the
#            round's nodes go back and the radius shrinks to a quarter of the
#            step. Both conditions are measured, one at a time:
#              - the value alone stalls short of the fixed point, because near
#                it the re-placed value can fall along the score (by 2.6e-4
#                against a promised +5e-5 on the 40-subject fixture);
#              - the residual alone reaches a worse point. On the gated-gaps
#                A14 config (random drift, T0MEANS and CINT: three effects a
#                subject, so the soft rule) it took the residual from 1.3 to
#                0.007 while the re-placed value fell by 0.9 nats, and ended
#                at an exact log posterior of -398.81, where both conditions
#                stop at -397.97 (Laplace -399.23, the step correction -400.60).
#                Before `_continuation_stiff_rule` gave the soft rule's stiff
#                complement a full product rule it was worse: the complement
#                sat at one node, whose fixed-node gradient leaves out its log
#                determinant, and the estimate ended 4.9 exact nats below the
#                Laplace optimum.
#   no worse the estimate is only reported if the re-placed value there is not
#            below the value at the Laplace optimum by more than `value_tol`;
#            otherwise the fit keeps the Laplace optimum (status `no_gain`),
#            with the quadrature log likelihood there reported, as the step
#            correction does when no step raises its objective.
#   skip     before any round: one gradient gives the gain the fixed-node
#            model predicts from the Laplace optimum (the residual there);
#            below `skip_gain` nothing moves (status `skipped`).
#   stop     when that residual is below the fit's certification tolerance,
#            or after `rounds` kept rounds or `attempts` rounds in all.
#            `stop_gain` > 0 stops instead on predicted gain in nats: the
#            residual below it once the last kept round realised less (a
#            round raising the re-placed objective by more being kept even if
#            the residual rose), a region shrunk until it promises less, or a
#            kept round gaining less. It is off (0) by default: on the bench's
#            paired grid (review/bench, quadcost 2026-09-27) it cut the
#            correction's time to a quarter and gained up to 2.5 exact nats on
#            gB8 and 0.7 on gC8, but lost 0.46 to 0.52 on gD3, 0.026 on gC2,
#            0.031 on gN3 and 0.06 to 0.16 on the AnomAuth default starts
#            against the rounds run to the tolerance.
#   guard    a continuation that moves the quadrature objective by more than
#            max(50, N/2) nats is reverted to the Laplace optimum with a
#            warning, keeping the rejected point, as bigIRT does.
#   Hessian  of the hybrid at the final point, its nodes placed there: exact
#            for the flagged units (the Louis identity over their fixed
#            nodes, each node's Hessian by forward-over-reverse sweeps),
#            central differences for the units left to Laplace and the prior
#            (see `ctsem_laplace_continuation_hessian`), and only where the
#            rounds reached the fixed point; the
#            covariance and draws then come from it, not from the Laplace
#            curvature at another point. Where they stopped short, the
#            estimate is not stationary on that objective and its Hessian need
#            not be concave, so none is taken: the fit's own uncertainty stays
#            and its draws are recentred, as the step correction's are.
#
# What moves when it applies: `fit$estimate$raw` is the continuation's
# estimate; `loglik`, `logposterior` and `subject_loglik` the quadrature values
# there (every unit by the rule, nodes placed there), with `loglik_laplace`
# kept; `cov`, `se` and the draws, when the rounds converged, from the
# continuation's Hessian, which is `fit$uncertainty$hessian` with
# `evaluated_at` the new estimate, and whose certification is the fit's
# (`hessian_at` in the record says which, as the step correction's does). `fit$laplace$correction` says what
# ran.
#
# Provenance of the constants. `nodes`, `tolerance` and `material` are the
# step correction's (see `.ctLaplaceCorrectDefaults`). `product_maxdim` (2),
# `soft_tau` (3.5) and `soft_maxdirs` (2) are the engine's, explained in
# laplace_continuation.jl. `rounds` (10) is twice bigIRT's cap of 5: its units
# are flat, and a nested unit here holds its leaves' conditional modes fixed
# with its nodes, so the iteration contracts more slowly there (the nested
# fixture's residual fell by a factor of 0.3 to 0.7 a round, where the
# one-level fixture's fell by 20 to 300). `attempts` (15) bounds the rounds a
# trust region may reject. `radius` (1 se, doubling to at most `radius_max` =
# 16 when a kept round ends on the boundary having gained at least three
# quarters of what its model promised, cut to a quarter of the step when a
# round is rejected) is the textbook trust-region schedule. Growing it on the
# residual instead -- only when a round cut it to a quarter -- held gated-gaps
# D3 at 0.25 se for six rounds whose gains matched their model's to 2%, and
# the rounds ran out 0.7 exact nats short of the best-known point. `maxiter`
# (100) caps one round; bigIRT's whole continuation took 30 to 60
# evaluations. Five is enough when `stop_gain` is on (gated-gaps A1, N1, C8:
# gradients from 54, 52 and 90 to 34, 52 and 79, end points within 0.04
# exact nats), not with the rounds run to the certification tolerance.
# `guard` and `guard_per_subject` are bigIRT's max(50, N/2).
# `skip_gain` (5e-3 nats) was set on the optimiser bench's default Laplace
# cells (dev/lapcontinue/calibrate-cost.R, dev1, 2026-09-27, one start each,
# the AnomAuth cells from their spurious maxima). Below it were cf_mixed,
# cf_ordinal, mvmix, ordinal, cf_binary and gD1, whose whole correction gained
# 6e-6 to 4.9e-3 exact nats -- at most 1.7 times the prediction -- and moved
# the estimate at most 0.075 se; the smallest prediction above it was jflat's
# 0.07 (gain 0.078). It is also the gain of a whitened Newton step of 0.1 se,
# the move `material` calls worth a line in print(). On the paired grid the
# twelve cells it skips kept their exact log likelihood to 0.005 nats at 0.05
# to 0.24 of the correction's time.
# `value_tol` (1e-3 nats) is the step correction's `gain_tol`: a change in the
# objective below it is not one to act on either way.
# `rtol` (1e-8) is the identifiability report's: a direction the Laplace
# curvature does not identify is held where the fit left it. `maxdim` (5) is
# the widest unit the correction scores, explained with the engine's default.
.ctLaplaceContinueDefaults <- list(nodes = 5L, tolerance = 0.01,
  product_maxdim = 2L, soft_tau = 3.5, soft_maxdirs = 2L, rounds = 10L,
  attempts = 15L, radius = 1, radius_max = 16, maxiter = 100L,
  material = 0.1, guard = 50, guard_per_subject = 0.5, rtol = 1e-8,
  value_tol = 1e-3, maxdim = 5L, stop_gain = 0, skip_gain = 5e-3)

# The directions a round moves in: the Laplace curvature's identified ones,
# each scaled to one of its standard errors, so the round's L-BFGS starts from
# that curvature and its trust region is measured in standard errors.
#
# Left out, and held where the fit left them: the directions the likelihood was
# measured flat along (`flat`, the columns of
# `fit$uncertainty$details$flatdirections$vectors`), and those whose curvature
# is below `rtol` of the largest -- the identifiability report's two rules for
# a direction the data do not identify. A direction of negative curvature is
# kept, whitened by its magnitude: the quadrature objective may rise along it,
# which is how a continuation leaves a saddle of the Laplace objective.
.ctLaplaceContinueBasis <- function(hessian, flat = NULL,
  rtol = .ctLaplaceContinueDefaults$rtol) {
  information <- -(hessian + t(hessian)) / 2
  n <- nrow(information)
  projector <- diag(n)
  measured <- 0L
  if (is.matrix(flat) && nrow(flat) == n && ncol(flat) > 0L &&
      all(is.finite(flat))) {
    decomposition <- qr(flat)
    q <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
    projector <- projector - q %*% t(q)
    measured <- decomposition$rank
  }
  kept <- projector %*% information %*% projector
  e <- eigen((kept + t(kept)) / 2, symmetric = TRUE)
  scale <- max(abs(e$values))
  if (!is.finite(scale) || scale <= 0) return(NULL)
  keep <- abs(e$values) > rtol * scale
  basis <- e$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(abs(e$values[keep])), nrow = sum(keep), ncol = sum(keep))
  list(basis = basis, kept = sum(keep), frozen = n - sum(keep),
    measured_flat = measured, negative = sum(e$values[keep] < 0))
}

# The rounds, shared by the fit (`.ctLaplaceContinue()`) and by
# `ctLaplaceCheck(refine = TRUE)`. `cont` is the engine's continuation object,
# placed at `est`; on return it is placed at the point returned. See the notes
# above `.ctLaplaceContinueDefaults` for what a round is, what it aims at and
# when one is kept. `residual` is half the squared whitened gradient with the
# nodes placed at the point: the Newton gain it predicts, zero at the answer.
.ctLaplaceContinueRun <- function(module, cont, est, basis, tol,
  control = .ctLaplaceContinueDefaults, verbose = 0L, skip_gain = 0,
  sink = NULL) {
  get <- .ctJuliaGet
  # `.ctLaplaceContinue()` passes its own sink, so the rounds continue the
  # same in-place line the screen started. `ctLaplaceCheck(refine = TRUE)`
  # passes none, so one is built here -- both report at the same default
  # verbosity rather than only when `verbose > 0`.
  if (is.null(sink)) sink <- .ctBackendProgressSink(verbose)
  x <- as.numeric(est)
  # Every stopping decision is a predicted gain in nats against one bar:
  # `stop_gain`, or the fit's certification tolerance where `stop_gain` is 0,
  # which is how the rounds ran before 2026-09-27 (they then stopped only on
  # that tolerance, 1e-6 by default, or by running out of rounds or attempts).
  stop_gain <- as.numeric(.ctJuliaOr(control$stop_gain, 0))
  bar <- if (stop_gain > 0) stop_gain else as.numeric(tol)
  start <- get(module$ctsem_laplace_continuation_info(cont))
  value <- as.numeric(start$quadrature)
  radius <- as.numeric(control$radius)
  kept <- 0L
  attempts <- 0L
  status <- "rounds"
  rows <- list()
  B <- .ctJuliaPut(as.matrix(basis))
  optimise <- function(from, stationary = FALSE, tol = bar) get(
    module$ctsem_laplace_continuation_optimize(cont, .ctJuliaNumericVector(from),
      B, radius, maxiter = as.integer(control$maxiter), tol = tol,
      stationary_only = stationary))
  row <- function(round, residual_after, keep, flagged, gain = NA_real_) data.frame(
    round = attempts, radius = radius,
    residual = as.numeric(round$start_gain), residual_after = residual_after,
    fixed_gain = as.numeric(round$value) - as.numeric(round$start_value),
    gain = gain,
    moved = as.numeric(round$moved), iterations = as.integer(round$iterations),
    f_calls = as.integer(round$f_calls), g_calls = as.integer(round$g_calls),
    kept = keep, flagged = as.integer(flagged))
  # The gain the fixed-node model predicts from the start: half the squared
  # whitened gradient, the Newton gain in nats. One gradient, and it decides
  # whether any round runs at all (`skip_gain`).
  residual <- as.numeric(optimise(x, stationary = TRUE)$start_gain)
  predicted <- residual
  # What the last kept round raised the re-placed objective by; none yet.
  realised <- 0
  if (!is.finite(residual)) {
    status <- "failed"
  } else if (residual < skip_gain) {
    status <- "skipped"
  } else repeat {
    # The residual is a Newton gain in the Laplace fit's metric, and where the
    # quadrature objective is flatter than that it understates what is left:
    # gated-gaps C2 stopped at a residual of 0.003 with 0.1 exact nats still to
    # gain. So a residual under the bar is converged only when the last round
    # also realised less than it; while rounds keep realising more, they go
    # on, each taking its step on the curvature its own L-BFGS measures.
    rising <- stop_gain > 0 && realised >= stop_gain
    if (residual < bar && !rising) { status <- "converged"; break }
    if (kept >= control$rounds) { status <- "rounds"; break }
    if (attempts >= control$attempts) { status <- "attempts"; break }
    # What the next round can promise inside its region: the Newton gain when
    # the whitened Newton step fits, else the gain of the whitened gradient's
    # step to the boundary on the same model. A region that rejections have
    # shrunk until it promises less than the bar is as done as a residual
    # below it; before this, those rounds ran on until the attempts ran out.
    reach <- sqrt(2 * residual)
    promise <- if (reach <= radius) residual else radius * reach - radius^2 / 2
    if (stop_gain > 0 && promise < stop_gain && !rising) {
      status <- "stalled"; break
    }
    attempts <- attempts + 1L
    # Below the bar, the round's own tolerance goes below the residual, so
    # that it takes its step rather than reporting itself stationary.
    round <- optimise(x, tol = if (residual < bar) residual / 100 else bar)
    residual <- as.numeric(round$start_gain)
    if (isTRUE(round$stationary)) {
      status <- "converged"
      rows[[length(rows) + 1L]] <- row(round, NA_real_, FALSE, start$nflagged)
      break
    }
    if (!is.finite(residual)) { status <- "failed"; break }
    xn <- as.numeric(round$minimizer)
    fixed <- as.numeric(round$value) - as.numeric(round$start_value)
    if (!is.finite(fixed) || fixed <= 0) {
      # The fixed-node model had a gradient to offer and no step along it: a
      # line search that found no decrease. There is nothing to re-place for.
      status <- "no_progress"
      rows[[length(rows) + 1L]] <- row(round, NA_real_, FALSE, start$nflagged)
      break
    }
    placed <- get(module[["ctsem_laplace_continuation_recentre!"]](cont,
      .ctJuliaNumericVector(xn)))
    after <- as.numeric(optimise(xn, stationary = TRUE)$start_gain)
    gain <- as.numeric(placed$quadrature) - value
    # Kept when the residual fell without the re-placed value falling, or --
    # whatever the residual did -- when that value rose by more than
    # `value_tol` (or `stop_gain` where that is on). On gated-gaps D3 two
    # rounds that raised it by 0.087 and 0.020 were rejected for a residual
    # that rose, the region shrank, and the rounds stopped 0.6 exact nats
    # short; with this, the paired grid of 2026-09-27 gained 0.9 to 2.5 exact
    # nats on gB8, 0.6 to 0.7 on gC8 and 2.3 on AnomAuth S1 from its spurious
    # start. A kept rise is a rise in the objective the fit reports, so this
    # cannot walk downhill, and each such round gains at least the bar, so it
    # cannot cycle.
    rise <- max(stop_gain, as.numeric(control$value_tol))
    keep <- is.finite(after) && is.finite(gain) &&
      ((after < residual && gain >= -as.numeric(control$value_tol)) ||
        gain > rise)
    rows[[length(rows) + 1L]] <- row(round, after, keep, placed$nflagged, gain)
    # At the same default verbosity as the rest of the correction: a round can
    # run for minutes (Charles's ordinal fixture, 2026-09-28), and this used to
    # print only at verbose > 0, so a default fit's screen showed nothing for
    # the whole of it. Kept as one detailed line rather than split into a
    # separate terse default and a separate detailed verbose>0 one -- it is
    # already a phrase, not a paragraph.
    if (!is.null(sink)) {
      sink(sprintf(paste0("Laplace continuation round %d: radius %.3g, ",
        "residual %.3g -> %.3g, value %+.4g, %s"), attempts, radius, residual,
        after, gain, if (keep) "kept" else "rejected"), "update")
    }
    realised <- if (keep) gain else 0
    if (keep) {
      kept <- kept + 1L
      x <- xn
      value <- as.numeric(placed$quadrature)
      # Grown when the round ended on the boundary and the re-placed value
      # rose by three quarters or more of what the fixed-node model promised:
      # the model is right, and the region is what holds the rounds back.
      if (isTRUE(round$boundary) && gain >= 0.75 * fixed) {
        radius <- min(2 * radius, as.numeric(control$radius_max))
      }
      residual <- after
      # Kept, but the objective with its nodes re-placed rose by less than
      # the bar: the rounds are closing on the fixed point without gaining
      # anything a fit could report. On gated-gaps A1 eight such rounds in a
      # row each took 0.015 se and lost 2e-4 to 9e-4 nats, until the rounds
      # ran out.
      if (stop_gain > 0 && gain < stop_gain && residual >= bar) {
        status <- "stalled"; break
      }
    } else {
      get(module[["ctsem_laplace_continuation_revert!"]](cont))
      radius <- 0.25 * max(as.numeric(round$moved), 1e-12)
      if (radius < 1e-3) { status <- "stalled"; break }
    }
  }
  info <- get(module$ctsem_laplace_continuation_info(cont))
  list(x = x, status = status, rounds = kept, attempts = attempts,
    radius = radius, residual = residual, predicted = predicted, bar = bar,
    trace = if (length(rows)) do.call(rbind, rows) else NULL,
    start = start, info = info, value = value)
}

# A log-posterior function on the continuation's hybrid objective, in the
# `lpgFunc` contract `.ctBackendUncertainty()` takes: see `.ctBackendLpgFunc()`,
# whose guard it copies.
.ctLaplaceContinueLpg <- function(module, cont) {
  function(fit, gradient = TRUE) {
    wantgrad <- isTRUE(gradient)
    function(parm) {
      result <- try(.ctJuliaGet(
        module$ctsem_laplace_continuation_evaluate(cont,
          .ctJuliaNumericVector(as.numeric(parm)), gradient = wantgrad)),
        silent = TRUE)
      failed <- inherits(result, "try-error") || !isTRUE(result$converged)
      value <- if (failed) NaN else as.numeric(result$value)[1L]
      grad <- if (failed || !wantgrad) NULL else as.numeric(result$gradient)
      if (!is.finite(value) || (wantgrad && (is.null(grad) ||
          length(grad) != length(parm) || any(!is.finite(grad))))) {
        value <- -1e100
        grad <- if (wantgrad) rep(0, length(parm)) else NULL
      }
      if (wantgrad) attributes(value) <- list(gradient = grad)
      value
    }
  }
}

.ctLaplaceContinue <- function(fit, cores = 1L, verbose = 0L,
  control = .ctLaplaceContinueDefaults) {
  started <- proc.time()[["elapsed"]]
  seconds <- function() proc.time()[["elapsed"]] - started
  est <- as.numeric(fit$estimate$raw)
  npar <- length(est)
  nsubjects <- length(fit$model_spec$subject_starts)
  module <- .ctJuliaModule(fit$model_spec$project)
  get <- .ctJuliaGet
  # The whole correction's progress line: the screen, each round and the final
  # Hessian all report through this one sink, at the same default verbosity as
  # the optimiser's own line, rather than only when `verbose > 0` -- which is
  # how a 22-minute correction (screen plus one round 4.8 min, then the
  # Hessian 18.8 min) printed nothing at all by default. One sink for the
  # whole function, not one per phase, so a shorter later phrase does not
  # leave the tail of a longer earlier one behind it when a line is
  # overwritten in place.
  sink <- .ctBackendProgressSink(verbose)
  # Unconditional and registered before any early return below: whatever the
  # last thing printed through `sink` was, this closes it so a warning or the
  # ordinary R prompt does not land on the same line. Safe to call with
  # nothing open -- `.ctProgressSink()`'s "break" is then a no-op.
  if (!is.null(sink)) on.exit(sink("", "break"), add = TRUE)
  failed <- function(phrase) {
    if (!is.null(sink)) sink("", "break")
    warning("Laplace continuation skipped: ", phrase, ". The uncorrected fit is ",
      "returned; see fit$laplace$correction.", call. = FALSE)
    fit$laplace$correction <- list(method = "quadrature", status = "failed",
      applied = FALSE, message = phrase, nodes = as.integer(control$nodes))
    fit
  }
  chunks <- suppressWarnings(as.integer(fit$optim$chunks)[1L])
  if (is.na(chunks) || chunks < 1L) chunks <- max(1L, as.integer(cores)[1L])
  previous <- .ctBackendSetMaxChunks(chunks)
  on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  # `sink` is non-NULL exactly when `.ctBackendReporting(verbose)` is, which
  # `verbose > 0` alone already satisfies -- so there is no separate case left
  # for a bare `verbose > 0` to reach for a plainer `message()` here.
  if (!is.null(sink)) {
    sink(sprintf("Laplace continuation screen (%d nodes) | %8s",
      control$nodes, .ctDuration(seconds())), "update")
  }
  cont <- try(module$ctsem_laplace_continuation(.ctJuliaObjective(fit),
    .ctJuliaNumericVector(est), nodes = as.integer(control$nodes),
    tolerance = as.numeric(control$tolerance),
    product_maxdim = as.integer(control$product_maxdim),
    soft_tau = as.numeric(control$soft_tau),
    soft_maxdirs = as.integer(control$soft_maxdirs),
    maxdim = as.integer(control$maxdim)), silent = TRUE)
  if (inherits(cont, "try-error")) return(failed("the quadrature could not be evaluated"))
  start <- get(module$ctsem_laplace_continuation_info(cont))
  record <- list(method = "quadrature", status = "exact", applied = FALSE,
    nodes = as.integer(start$nodes), tolerance = as.numeric(start$tolerance),
    screen = as.numeric(start$screen), gap = as.numeric(start$gap),
    laplace_estimate = est, loglik_laplace = as.numeric(fit$estimate$loglik),
    flagged = as.integer(start$nflagged), units = as.integer(start$nunits),
    rule_failures = as.integer(start$rule_failures),
    wide = as.integer(start$nwide), maxdim = as.integer(start$maxdim))
  # Wall seconds by stage as each is passed -- the screen, the rounds, the
  # continuation's Hessian, its certification, the uncertainty redrawn from it
  # -- and `total` when the record is written.
  stages <- c(screen = seconds())
  timed <- function() c(stages, total = seconds())
  since <- function(from) proc.time()[["elapsed"]] - from
  # Units wider than `maxdim` are not scored and keep the Laplace term (see
  # laplace_continuation.jl for the measured cost). Said in one line, since
  # the correction then covers less than the fit, or none of it.
  if (record$wide > 0L) {
    if (!is.null(sink)) sink("", "break")
    message(sprintf(paste0("Quadrature correction: %d of %d %s kept the ",
      "Laplace term, having more than %d random effects."), record$wide,
      record$units, if (record$units == nsubjects) "subjects" else "groups",
      record$maxdim))
  }
  if (record$wide == record$units) {
    record$status <- "too_wide"
    record$seconds <- timed()
    fit$laplace$correction <- record
    return(fit)
  }
  if (!is.finite(record$screen)) {
    return(failed("the quadrature was not finite at the estimate"))
  }
  # Passed, as the step correction's screen passes: Laplace is exact here to
  # the tolerance, and the fit is left alone to the bit.
  if (record$screen <= record$tolerance) {
    record$seconds <- timed()
    fit$laplace$correction <- record
    return(fit)
  }
  # The quadrature log likelihood wherever the screen found a gap, as the step
  # correction reports it, whether or not anything moves.
  report_quadrature <- function(fit, info) {
    units <- as.numeric(info$quadrature_units)
    if (!all(is.finite(units))) return(fit)
    fit$estimate$loglik_laplace <- record$loglik_laplace
    fit$estimate$loglik <- sum(units)
    fit$estimate$logposterior <- as.numeric(info$quadrature)
    subjects <- as.numeric(info$quadrature_subjects)
    if (length(subjects) == length(fit$estimate$subject_loglik) &&
        all(is.finite(subjects))) fit$estimate$subject_loglik <- subjects
    fit$estimate$loglik_method <- "quadrature"
    fit
  }
  hessian <- .ctBackendStoredHessian(fit, est)
  if (is.null(hessian) || !all(is.finite(hessian))) {
    record$status <- "no_hessian"
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- record
    return(fit)
  }
  flat <- fit$uncertainty$details$flatdirections$vectors
  basis <- .ctLaplaceContinueBasis(hessian, flat, rtol = control$rtol)
  if (is.null(basis) || basis$kept < 1L) {
    record$status <- "no_directions"
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- record
    return(fit)
  }
  tol <- .ctBackendGapTolerance(fit)
  # What the continuation cost, in the engine's own counts -- whole-objective
  # values and gradients, placements of the rule, and the member likelihood
  # values and reverse sweeps they took -- and its seconds in each kind of
  # call (`engine_seconds`, the screen's placement among the placements),
  # beside the wall seconds by stage.
  account <- function(record) {
    after <- get(module$ctsem_laplace_continuation_info(cont))
    record$evaluations <- c(values = as.integer(after$value_calls),
      gradients = as.integer(after$gradient_calls),
      placements = as.integer(after$recentres),
      member_values = as.numeric(after$member_values),
      member_sweeps = as.numeric(after$member_sweeps),
      refused = as.integer(after$refused))
    record$engine_seconds <- c(values = as.numeric(after$seconds_values),
      gradients = as.numeric(after$seconds_gradients),
      placements = as.numeric(after$seconds_placements),
      hessian = as.numeric(after$seconds_hessian))
    record$seconds <- timed()
    record
  }
  rounds_started <- proc.time()[["elapsed"]]
  run <- try(.ctLaplaceContinueRun(module, cont, est, basis$basis, tol,
    control = control, verbose = verbose, sink = sink,
    skip_gain = as.numeric(.ctJuliaOr(control$skip_gain, 0))), silent = TRUE)
  stages["rounds"] <- since(rounds_started)
  if (inherits(run, "try-error")) {
    return(failed(paste0("a round could not be evaluated (",
      trimws(as.character(run)), ")")))
  }
  record$predicted_gain <- as.numeric(run$predicted)
  record$skip_gain <- as.numeric(.ctJuliaOr(control$skip_gain, 0))
  record$stop_gain <- as.numeric(run$bar)
  # Skipped: the fixed-node model promised less than `skip_gain` from the
  # Laplace optimum, so no round ran. The estimate stays, and the log
  # likelihood reported is the quadrature one there, as where a round ran and
  # gained nothing.
  if (identical(run$status, "skipped")) {
    record$status <- "skipped"
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    record$logposterior_quadrature <- as.numeric(start$quadrature)
    record$gap_reported <- as.numeric(start$quadrature) - as.numeric(start$laplace)
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- account(record)
    return(fit)
  }
  info <- run$info
  change <- as.numeric(info$quadrature) - as.numeric(start$quadrature)
  limit <- max(control$guard, control$guard_per_subject * nsubjects)
  record$rounds <- run$rounds
  record$attempts <- run$attempts
  record$radius <- run$radius
  record$stationarity <- run$residual
  record$trace <- run$trace
  record$continuation <- run$status
  record$frozen_directions <- as.integer(basis$frozen)
  record$flagged <- as.integer(info$nflagged)
  record$flagged_units <- as.integer(info$flagged)[seq_len(info$nflagged)]
  record$bar <- as.numeric(info$bar)
  record$residual <- as.numeric(info$residual)
  record$soft_blocks <- as.integer(info$soft_blocks)
  record$guard <- list(limit = limit, change = change,
    fired = !is.finite(change) || abs(change) > limit)
  # The run moved nothing: no round was kept, the rounds ended lower on the
  # objective than they began, or the guard reverts it. The quadrature log
  # likelihood at the Laplace optimum is reported, as the step correction
  # reports it when no step raised the objective.
  worse <- is.finite(change) && change < -as.numeric(control$value_tol)
  record$ended_lower <- worse
  moved <- run$rounds > 0L && any(run$x != est) && !worse
  if (isTRUE(record$guard$fired) || !moved) {
    if (isTRUE(record$guard$fired)) {
      if (!is.null(sink)) sink("", "break")
      warning(sprintf(paste0("The Laplace continuation moved the quadrature ",
        "objective by %.3g nats, which no approximation error of this size ",
        "explains; the fit keeps the Laplace optimum, and the rejected point is ",
        "fit$laplace$correction$rejected_estimate."), change), call. = FALSE)
      record$rejected_estimate <- run$x
      record$status <- "reverted"
    } else {
      record$status <- "no_gain"
    }
    record$loglik_quadrature <- sum(as.numeric(start$quadrature_units))
    record$logposterior_quadrature <- as.numeric(start$quadrature)
    record$gap_reported <- as.numeric(start$quadrature) - as.numeric(start$laplace)
    fit <- report_quadrature(fit, start)
    fit$laplace$correction <- account(record)
    return(fit)
  }
  x <- run$x
  # The continuation's own curvature only where it reached its fixed point.
  # Stopped short, the estimate is not stationary on the objective the Hessian
  # is of, and on gated-gaps A14 that Hessian had two directions of positive
  # curvature: ten parameters would have reported no interval, and the fit
  # `converged: FALSE`, where the Laplace curvature at the fit's optimum gives
  # them all one. So a continuation that stopped short keeps the fit's
  # uncertainty and recentres its draws, as the step correction does, and
  # takes no Hessian at all: on wide blocks that is a large share of a
  # correction's cost, and `stationarity` in the record already says how far
  # from stationary the estimate is.
  reached <- identical(run$status, "converged")
  hc <- NULL
  if (reached) {
    # By forward differences this Hessian was 18.8 minutes of a 22-minute
    # correction on the ordinal fixture that found the silence; it is exact
    # now (laplace_continuation.jl) and can still take minutes. A start and a
    # done line around the call, since the engine reports no count from
    # inside it; the sink reports at every verbosity that reports at all.
    if (!is.null(sink)) {
      sink(sprintf("Laplace continuation hessian | %8s",
        .ctDuration(seconds())), "update")
    }
    hessian_started <- proc.time()[["elapsed"]]
    hc <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_laplace_continuation_hessian(cont, .ctJuliaNumericVector(x)))),
      npar, npar), silent = TRUE)
    stages["hessian"] <- since(hessian_started)
    if (!is.null(sink)) {
      sink(sprintf("Laplace continuation hessian done in %s",
        .ctDuration(stages[["hessian"]])), "done")
    }
  }
  final <- get(module$ctsem_laplace_continuation_evaluate(cont,
    .ctJuliaNumericVector(x), gradient = TRUE))
  record$status <- "continued"
  record$applied <- TRUE
  se <- as.numeric(fit$estimate$se)
  delta <- x - est
  delta_se <- if (length(se) == npar) ifelse(se > 0, delta / se, NA_real_) else
    rep(NA_real_, npar)
  names(delta) <- names(delta_se) <- names(fit$estimate$se)
  record$delta <- delta
  record$delta_se <- delta_se
  record$material <- isTRUE(max(abs(delta_se), na.rm = TRUE) >= control$material)
  record$loglik_quadrature <- sum(as.numeric(info$quadrature_units))
  record$logposterior_quadrature <- as.numeric(info$quadrature)
  record$gap_reported <- as.numeric(info$quadrature) - as.numeric(info$laplace)
  # The objective the estimate maximises: the hybrid, whose unflagged units are
  # Laplace. It differs from the quadrature value above by at most `residual`.
  record$objective <- as.numeric(final$value)
  fit$estimate$raw <- x
  fit <- report_quadrature(fit, info)
  usable <- reached && !inherits(hc, "try-error") && all(is.finite(hc))
  method <- .ctJuliaOr(fit$uncertainty$settings$method, "hessian")
  # Present even when empty: `$` partial-matches, and without it
  # `correction$hessian` would answer with `hessian_at`.
  record["hessian"] <- list(NULL)
  if (usable) {
    # The certification of the estimate on the objective it maximises: its
    # gradient against its own Hessian, with the flat-direction probe walking
    # the hybrid.
    certification_started <- proc.time()[["elapsed"]]
    gradient <- as.numeric(final$gradient)
    gap <- .ctBackendOptimGap(hc, gradient)
    probe <- NULL
    # Only along a direction the curvature does not trust. With every one
    # trusted the residual is the rounding left by projecting onto a complete
    # eigenbasis, and the probe walked that: three values of the hybrid, a
    # fifteenth of ord4's correction, and a gain of 0 on every converged
    # continuation of the bench (controbust-c2a67356). A negative direction
    # makes the verdict `notmaximum` whatever a probe finds.
    if (isTRUE(gap$ok) && isTRUE(gap$nflat > 0L) &&
        isTRUE(gap$residual_norm > 0)) {
      out <- try(get(module$ctsem_flat_probe(cont, .ctJuliaNumericVector(x),
        .ctJuliaNumericVector(gap$residual))), silent = TRUE)
      if (!inherits(out, "try-error")) probe <- .ctBackendProbeFields(out, npar)
    }
    parnames <- try(.ctBackendRawParameterNames(fit, npar), silent = TRUE)
    if (inherits(parnames, "try-error")) parnames <- NULL
    # Certified to the bar the rounds stopped on, not the Laplace fit's: the
    # quadrature rule is not resolved to 1e-6 nats, and a continuation that
    # met its own stopping rule is not a failure to converge.
    certification <- .ctBackendCertificationRecord(gap, probe, tolerance = run$bar,
      saturated = isTRUE(fit$optim$saturated), overshot = FALSE,
      parnames = parnames)
    record$certification <- certification
    record$hessian <- hc
    stages["certification"] <- since(certification_started)
  }
  own <- usable && identical(method, "hessian")
  record$hessian_at <- if (own) "estimate" else "laplace_estimate"
  if (own) {
    uncertainty_started <- proc.time()[["elapsed"]]
    settings <- fit$uncertainty$settings
    finishsamples <- .ctJuliaOr(settings$finishsamples,
      if (!is.null(fit$estimate$rawposterior)) nrow(fit$estimate$rawposterior) else 1000L)
    fit$uncertainty <- list(hessian = hc, evaluated_at = x,
      certification = certification)
    fit <- .ctBackendUncertainty(fit, uncertainty = "hessian", draws = "normal",
      finishsamples = finishsamples, cores = chunks,
      control = .ctJuliaOr(settings$control, list()), verbose = verbose,
      lpg = .ctLaplaceContinueLpg(module, cont))
    record$draws <- "redrawn"
    fit$uncertainty$details$laplace_correction <- list(
      draws = paste0("normal draws about the continuation's estimate, from the ",
        "Hessian of its quadrature objective there"),
      nodes = record$nodes)
    stages["uncertainty"] <- since(uncertainty_started)
  } else {
    # Recentred, as the step correction does, when the continuation stopped
    # short of its fixed point, has no Hessian to offer, or the fit's
    # uncertainty did not come from one.
    post <- fit$estimate$rawposterior
    if (!is.null(post) && ncol(post) == npar) {
      fit$estimate$rawposterior <- sweep(post, 2L, delta, "+")
    }
    record$draws <- "recentred"
    if (!is.null(fit$uncertainty)) {
      fit$uncertainty$details$laplace_correction <- list(
        draws = paste0("recentred on the continuation's estimate; the ",
          "covariance is the fit's own at fit$uncertainty$evaluated_at"),
        nodes = record$nodes)
    }
  }
  fit$laplace$correction <- account(record)
  fit
}
