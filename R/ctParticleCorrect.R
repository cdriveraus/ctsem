# `ctFitUncertainty(fit, 'is', control = list(likelihood = 'particle'))`:
# applying the correction ctParticleLik() measures.
#
# The filter's posterior draws -- however they were produced -- approximate the
# posterior under the filter's likelihood, which for a nonlinear or non-Gaussian
# model is itself an approximation. Importance sampling closes both gaps in one
# step: weight each draw by the exact posterior density (the particle-filter
# likelihood times the prior) over the density the draw came from. That the
# particle likelihood is a noisy but unbiased *estimate* costs variance in the
# weights, not bias (Tran, Scharth, Pitt & Kohn 2014, "importance sampling
# squared"), so the corrected draws still target the exact posterior; the
# effective sample size says how much the noise and the approximation cost.
#
# The proposal is the adaptive mixture `uncertainty = 'is'` uses against the
# fit's own objective (`imis_is`), started from a normal `scale` times wider
# than the fit's covariance: a proposal no wider than the Hessian cannot correct
# a posterior wider than it, the same finding that set `imisScaleInit = 1.5`.
# Each proposal draw costs one particle-filter evaluation.
#
# Two things that were tried and are wrong, recorded so they stay untried. A
# shared particle-filter seed across draws (common random numbers) makes the
# filter's Monte Carlo error a fixed random function of the parameters, which
# tilts the weighted posterior by that seed's error surface instead of
# averaging out; more draws do not cure it, and it cost 0.45 standard errors on
# a two-parameter test fit. Independent seeds per draw are what the unbiasedness argument
# needs. And weights of `exp(particle - filter)` on normal draws would correct
# the filter's approximation while leaving the normal one in place, a mixture
# that is neither posterior.

# Returns the fit with the draws, covariance and standard errors replaced, the
# estimate at the importance-weighted posterior mean unless `correct_estimate =
# FALSE` (as `uncertainty = 'sample'` moves it to the posterior mean), and a
# printable `particle_correction` recording the ctParticleLik() check at the
# original estimate, the effective sample size, every evaluation's particle
# and filter log likelihoods, each parameter's shift in original standard
# errors and its ratio of new to old standard error, and the estimate before.
# `fit$estimate$loglik` stays the filter's value at the original estimate.
.ctParticleCorrect <- function(fit, particles = 1000, substeps = 10,
  transition = c("exponential", "euler"), seed = 1, finishsamples = NULL,
  scale = 1.5, target_ess = 200, nbatch = 200L, maxiter = 10L,
  correct_estimate = TRUE, cores = NULL, verbose = 0L) {

  transition <- match.arg(transition)
  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$laplace)) {
    stop("likelihood = 'particle' is not available for intoverpop = 'laplace' or ",
      "'none' fits: the particles carry random effects only as augmented states. ",
      "likelihood = 'quadrature' is the Laplace route's exact reference.",
      call. = FALSE)
  }
  particles <- as.integer(particles)[1L]
  substeps <- as.integer(substeps)[1L]
  if (is.na(particles) || particles < 2L) stop("particles must be at least 2.", call. = FALSE)
  if (is.na(substeps) || substeps < 1L) stop("substeps must be at least 1.", call. = FALSE)
  posterior <- fit$estimate$rawposterior
  hasdraws <- !is.null(posterior) && length(dim(posterior)) == 2L && nrow(posterior) >= 2L
  covariance <- fit$estimate$cov
  if (is.null(covariance)) {
    stop("likelihood = 'particle' needs the fit's covariance as its first ",
      "proposal. Run ctFitUncertainty(fit) first.", call. = FALSE)
  }
  est <- as.numeric(fit$estimate$raw)
  parnames <- if (hasdraws && !is.null(colnames(posterior))) colnames(posterior) else
    tryCatch(as.character(.ctFitRawParNames(fit)), error = function(e) NULL)
  if (length(parnames) != length(est)) parnames <- NULL

  if (!is.null(cores)) {
    previous <- .ctBackendSetMaxChunks(max(1L, cores))
    on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  }
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  settings <- list(particles = particles, substeps = substeps, transition = transition,
    seed = as.integer(seed)[1L])

  # The check at the original estimate, with the same settings, so the gap
  # reported and the numbers used to close it agree on what "the particle
  # likelihood" is.
  check <- ctParticleLik(fit, particles = particles, substeps = substeps,
    transition = transition, seed = seed)
  if (verbose > 0) {
    message(sprintf(paste("Particle likelihood at the estimate %.3f (se %.3f),",
      "filter %.3f, difference %.3f"), check$loglik, check$se, check$fit_loglik,
      check$difference))
  }
  se0 <- fit$estimate$se
  if (is.null(se0) || length(se0) != length(est)) {
    se0 <- if (hasdraws) apply(posterior, 2, stats::sd) else rep(NA_real_, length(est))
  }
  before <- list(raw = est, cov = covariance, se = as.numeric(se0))
  # Every evaluation is kept, so the record can show the particle and filter
  # likelihoods per proposal draw; `imis_is` itself only keeps the target
  # density.
  memo <- new.env(parent = emptyenv())
  memo$rows <- list()
  lp <- function(parm) {
    # Its own seed per proposal draw: the estimates have to be independent
    # across draws for the self-normalised weights to be consistent (above).
    r <- .ctParticleEvaluate(module, objective, matrix(as.numeric(parm), nrow = 1L),
      settings, seeds = settings$seed + length(memo$rows))
    memo$rows[[length(memo$rows) + 1L]] <- r
    value <- r$posterior - r$filter + r$particle
    # A proposal draw the model cannot be evaluated at should carry no
    # weight, and `exp(-1e100 - c)` is zero: the Laplace path's own guard.
    if (!is.finite(value)) -1e100 else value
  }
  if (is.null(finishsamples)) finishsamples <- if (hasdraws) nrow(posterior) else 1000L
  if (verbose > 0) {
    message("Importance sampling against the particle posterior (", nbatch,
      " draws per iteration, one particle filter each)")
  }
  # Shared with `uncertainty='is'` against the fit's own objective; see
  # .ctOptimImisDraws(). This function's own parts are the
  # density (`lp`, one particle filter per proposal draw) and the remedy.
  #
  # The proposal is already widened by `scale^2` here, so `scaleInit = 1`
  # rather than compounding two scalings.
  drawn <- .ctOptimImisDraws(lp, centre = est,
    cov = ctOptimSafeCov(as.matrix(covariance)) * scale^2,
    finishsamples = finishsamples, nbatch = nbatch, target_ess = target_ess,
    maxiter = maxiter, scaleInit = 1, tailScale = 1.2, df = Inf,
    verbose = verbose, diagPlots = FALSE,
    remedy = paste0("Treat the corrected intervals as indicative; a wider ",
      "scale or more iterations may help."))
  is_res <- drawn$is_res
  samples <- drawn$samples
  if (is.null(samples) || !nrow(samples)) {
    stop("Importance sampling returned no usable draws against the particle ",
      "posterior. Widen the proposal (scale) or raise nbatch.", call. = FALSE)
  }
  theta <- is_res$full_theta
  w <- as.numeric(is_res$full_weights)
  ess <- drawn$ess
  paretok <- drawn$k
  newmean <- as.numeric(is_res$mean)
  newcov <- drawn$cov
  ev <- do.call(rbind, memo$rows)
  evaluations <- data.frame(particle = ev$particle, filter = ev$filter,
    difference = ev$particle - ev$filter, se = ev$se,
    weight = if (nrow(ev) == length(w)) w else NA_real_)

  shift_se <- (newmean - est) / before$se
  width_ratio <- sqrt(diag(newcov)) / before$se
  if (!is.null(parnames)) {
    names(shift_se) <- parnames
    names(width_ratio) <- parnames
    names(newmean) <- parnames
  }

  if (isTRUE(correct_estimate)) fit$estimate$raw <- as.numeric(newmean)
  fit$estimate$cov <- newcov
  fit$estimate$se <- sqrt(diag(newcov))
  fit$estimate$rawposterior <- samples
  fit <- .ctFitNameRawUncertainty(fit)
  # The constrained draws describe whatever raw draws they were built from, so
  # they are refreshed here rather than left to disagree with the ones above.
  fit$transformedpars <- .ctBackendConstrain(fit)

  fit$particle_correction <- structure(list(
    check = check, draws = "imis", proposal = "adaptive mixture (imis)",
    scale = scale, particles = particles, substeps = substeps, transition = transition,
    seed = settings$seed, ndraws = length(w), ess = ess,
    effective_fraction = if (is.finite(ess)) ess / length(w) else NA_real_,
    draws_evaluated = theta, weights = w, evaluations = evaluations, shift_se = shift_se,
    width_ratio = width_ratio,
    largest_shift_se = if (all(is.na(shift_se))) NA_real_ else max(abs(shift_se), na.rm = TRUE),
    corrected_estimate = isTRUE(correct_estimate), before = before,
    finishsamples = nrow(samples)), class = "ctParticleCorrection")

  if (!is.null(fit$uncertainty)) {
    fit$uncertainty$cov <- newcov
    # The draws on the fit are now importance-sampled, so the record of how
    # they were made has to say so; summaries describe the fit from it.
    fit$uncertainty$draws <- "imis"
    fit$uncertainty$proposal_cov <- ctOptimSafeCov(as.matrix(covariance)) * scale^2
    fit$uncertainty$imis <- is_res
    fit$uncertainty$details$importance_sampling <- list(ess = ess,
      pareto_k = paretok, df_used = is_res$df_used,
      covariance = "weighted importance-sampling covariance")
    fit$uncertainty$details$particle_correction <- list(draws = "imis", ess = ess,
      pareto_k = paretok, particles = particles, substeps = substeps, transition = transition,
      target = "particle-filter posterior")
  }
  fit
}

#' @export
print.ctParticleCorrection <- function(x, ...) {
  cat("Particle-filter correction applied\n")
  cat("  particle filter       ", x$particles, " particles, ", x$substeps,
    " substeps, ", x$transition, " transition\n", sep = "")
  cat("  at the estimate       particle ", format(x$check$loglik, digits = 8),
    " (se ", format(x$check$se, digits = 3), ")  filter ",
    format(x$check$fit_loglik, digits = 8), "  difference ",
    format(x$check$difference, digits = 4), "\n", sep = "")
  cat("  draws                 ", x$draws, " (", x$proposal, "; ", x$ndraws,
    " evaluated, ", x$finishsamples, " returned)\n", sep = "")
  if (is.finite(x$ess)) {
    cat("  effective sample size ", format(x$ess, digits = 4),
      if (is.finite(x$effective_fraction))
        paste0(" (", format(100 * x$effective_fraction, digits = 3), "% of draws)")
      else "", "\n", sep = "")
  }
  if (length(x$shift_se) && !all(is.na(x$shift_se))) {
    worst <- order(-abs(x$shift_se))[seq_len(min(5L, length(x$shift_se)))]
    cat("  largest shifts, in original standard errors, with the ratio of new to old:\n")
    print(data.frame(parameter = if (is.null(names(x$shift_se))) worst else names(x$shift_se)[worst],
      shift_se = unname(x$shift_se[worst]), se_ratio = unname(x$width_ratio[worst])),
      row.names = FALSE, digits = 3)
  }
  if (!x$corrected_estimate) {
    cat("  estimate left at the original point; draws corrected only\n")
  }
  invisible(x)
}

# The particle and filter likelihoods at each row of `theta`, in batches of one
# bridge call each, so a long run reports progress and a short one does not pay
# a round trip per draw.
.ctParticleEvaluate <- function(module, objective, theta, settings, seeds,
  verbose = 0L, batch = 25L) {
  n <- nrow(theta)
  seeds <- as.integer(seeds)
  if (length(seeds) != n) stop("one particle-filter seed per draw is required", call. = FALSE)
  starts <- seq.int(1L, n, by = batch)
  out <- vector("list", length(starts))
  began <- Sys.time()
  for (b in seq_along(starts)) {
    rows <- starts[b]:min(n, starts[b] + batch - 1L)
    values <- t(theta[rows, , drop = FALSE])
    res <- .ctJuliaGet(module$ctsem_particle_batch(objective, values,
      particles = settings$particles, substeps = settings$substeps,
      transition = settings$transition, seed = seeds[rows]))
    out[[b]] <- data.frame(particle = as.numeric(res$particle), se = as.numeric(res$se),
      ess_min = as.numeric(res$ess_min), filter = as.numeric(res$filter),
      posterior = as.numeric(res$posterior))
    if (verbose > 0 && n > batch) {
      done <- max(rows)
      elapsed <- as.numeric(difftime(Sys.time(), began, units = "secs"))
      message(sprintf("  %d of %d draws evaluated, %.0f s elapsed, about %.0f s remaining",
        done, n, elapsed, elapsed / done * (n - done)))
    }
  }
  do.call(rbind, out)
}
