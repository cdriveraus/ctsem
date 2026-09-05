# ctParticleCorrect(): applying the correction ctParticleLik() measures.
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
# The density the draws came from is the one thing this has to get right, and
# it depends on how the fit was finished. Curvature methods draw from a normal
# around the estimate with the fit's covariance; `ctOptimUncertainty('is')`,
# the bootstrap and `ctSample()` leave draws distributed as the filter posterior
# itself. `.ctParticleProposal()` reads that off the fit's own record rather
# than assuming, because assuming the wrong one gives weights that look fine
# and describe nothing.
#
# Two draw strategies, as in `ctLaplaceCorrect()`: one batch of weighted draws,
# or the adaptive importance sampler the package already uses (`imis_is`,
# shared with `ctOptimUncertainty('is')`). One batch is the default because it
# costs one particle-filter evaluation per draw and nothing else; the sampler is
# the remedy when too few of those draws survive.
#
# For a curvature fit the batch is drawn afresh from a normal `scale` times
# wider than the fit's covariance rather than reusing the fit's own draws. The
# fit's draws are not data, they are a cheap normal sample, and a proposal no
# wider than the Hessian cannot correct a posterior wider than it -- the same
# finding that set `imisScaleInit = 1.5` for `ctOptimUncertainty('is')`.
# Measured here on a two-parameter linear fit: exact weights on the fit's own
# 200 Hessian draws left the mean 0.24 standard errors short of the true
# posterior mean and the spread at 0.87 of the true spread.
#
# Two things that were tried and are wrong, recorded so they stay untried. A
# shared particle-filter seed across draws (common random numbers) makes the
# filter's Monte Carlo error a fixed random function of the parameters, which
# tilts the weighted posterior by that seed's error surface instead of
# averaging out; more draws do not cure it, and it cost 0.45 standard errors on
# the same fit. Independent seeds per draw are what the unbiasedness argument
# needs. And weights of `exp(particle - filter)` on normal draws would correct
# the filter's approximation while leaving the normal one in place, a mixture
# that is neither posterior.

#' Correct a julia fit's posterior draws against the particle-filter likelihood
#'
#' \code{\link{ctParticleLik}} measures how far the filter's likelihood is from
#' the exact one at the estimate. This function applies the correction to the
#' posterior draws behind the fit's intervals, by importance sampling: each draw
#' is weighted by the exact posterior density (the particle-filter likelihood
#' times the prior) over the density the draw came from, and the weighted draws
#' describe the exact posterior. The particle likelihood is an unbiased
#' estimate, so its Monte Carlo noise widens the weights without biasing them.
#'
#' Two ways to obtain the draws. \code{draws = 'reweight'} (the default) weights
#' one batch of draws, and what the batch is depends on how the fit's own draws
#' were made, which is read from the fit's record. After a curvature method
#' (\code{\link{ctOptimUncertainty}} with \code{'hessian'}, \code{'opg'},
#' \code{'sandwich'} or \code{'surrogate'}) the batch is drawn afresh from a
#' normal around the estimate \code{scale} times wider than the fit's
#' covariance, and the weights correct both the normal approximation and the
#' filter's. The fit's own draws are not reused, because a proposal no wider
#' than the curvature cannot correct a posterior wider than it. After
#' \code{ctOptimUncertainty(uncertainty = 'is')}, a bootstrap, or
#' \code{\link{ctSample}}, the fit's draws already sit on the filter posterior,
#' so they are reused and the weight is the likelihood ratio alone. \code{draws
#' = 'imis'} draws afresh, adaptively, from a mixture that starts at the fit's
#' covariance and moves toward the particle posterior, at one particle-filter
#' evaluation per proposal draw; it is the same sampler
#' \code{ctOptimUncertainty(uncertainty = 'is')} and \code{\link{ctLaplaceCorrect}}
#' use, and the remedy when one batch leaves too few effective draws.
#'
#' The cost is one particle-filter evaluation per draw, so this takes minutes to
#' hours where \code{ctParticleLik} took seconds. Each draw gets its own
#' particle-filter seed (\code{seed}, \code{seed + 1}, ...): the weights are
#' consistent only if the estimates are independent across draws, and a shared
#' seed would tilt the whole posterior by that seed's Monte Carlo error. Fewer
#' \code{particles} per draw cost weight variance rather than bias, so a few
#' hundred often serves here where the check wanted thousands. Read \code{ess}
#' before trusting the result: the effective sample size is how many draws the
#' corrected posterior rests on.
#'
#' @param fit A \code{ctJuliaFit} with posterior draws (\code{\link{ctFit}}
#'   computes them by default) and the default \code{intoverpop}; random
#'   effects carried as augmented states are integrated by the particles.
#'   Laplace-route fits are refused, as in \code{ctParticleLik}.
#' @param draws \code{'reweight'} reuses the fit's draws; \code{'imis'} draws
#'   afresh by adaptive importance sampling. See Details.
#' @param particles Particles per particle-filter evaluation.
#' @param substeps Transition steps per observation interval in the particle
#'   filter, floored by the fit's own \code{maxtimestep} or mesh.
#' @param transition The particle step, as in \code{ctParticleLik}.
#' @param seed Seed of the particle filter at the first draw; draw \code{i} uses
#'   \code{seed + i - 1}.
#' @param nsamples \code{draws = 'reweight'} only. Draws to evaluate: fresh
#'   normal draws after a curvature fit, or a random subset of this many of the
#'   fit's draws otherwise. Defaults to as many draws as the fit carries. The
#'   draws, the subset and the resampling below use R's generator, so
#'   \code{set.seed()} fixes them.
#' @param finishsamples Draws to return. Defaults to as many as the fit carries.
#'   Under \code{'reweight'} they are a systematic resample of the weighted
#'   draws, so a small effective sample size shows as repeated rows.
#' @param scale Proposal scale multiplier on the fit's covariance: for the
#'   fresh normal draws of \code{'reweight'} after a curvature fit, and for the
#'   first proposal of \code{'imis'}. Above one by default, because a proposal
#'   narrower than its target cannot correct it.
#' @param target_ess \code{draws = 'imis'} only. Effective sample size at which
#'   sampling stops.
#' @param nbatch \code{draws = 'imis'} only. Proposal draws, and so
#'   particle-filter evaluations, per iteration.
#' @param maxiter \code{draws = 'imis'} only. Iteration cap.
#' @param correct_estimate Move the point estimate to the importance-weighted
#'   posterior mean, as \code{ctSample} moves it to the posterior mean. When
#'   \code{FALSE} only the draws, covariance and standard errors change.
#' @param cores Engine threads. The particle filter splits the subjects across
#'   them; its result does not depend on the count.
#' @param verbose Integer; 1 or more reports progress with the time remaining.
#'
#' @return The fit, with \code{estimate$rawposterior} replaced by the corrected
#'   draws, \code{estimate$cov} and \code{estimate$se} by their weighted
#'   covariance, \code{estimate$raw} at the weighted mean unless
#'   \code{correct_estimate = FALSE}, \code{transformedpars} refreshed, and a
#'   printable \code{particle_correction} entry recording: the
#'   \code{ctParticleLik} check at the original estimate; how the draws were
#'   obtained and which proposal density was used; the number of draws
#'   evaluated, the effective sample size and the fraction of draws it makes;
#'   the draws evaluated (\code{draws_evaluated}), their normalised weights and,
#'   per draw, the particle and filter log likelihoods and their difference
#'   (\code{evaluations}); each parameter's
#'   shift in original standard errors and the ratio of new to old standard
#'   error; and the estimate, covariance and standard errors before correction.
#'   \code{fit$estimate$loglik} stays the filter's value at the original
#'   estimate; \code{ctParticleLik} on the returned fit gives the particle
#'   likelihood at the corrected one.
#'
#' @seealso \code{\link{ctParticleLik}} measures without applying.
#'   \code{\link{ctLaplaceCorrect}} is the same shape of correction for the
#'   Laplace route's random effects.
#'
#' @examples
#' \dontrun{
#' fit <- ctFit(data, model, backend = 'julia')
#' ctParticleLik(fit)                       # the gap at the estimate
#' corrected <- ctParticleCorrect(fit, particles = 500, substeps = 5)
#' corrected$particle_correction            # effective sample size, shifts
#' summary(corrected)                       # intervals from the corrected draws
#'
#' # If too few draws survive the reweighting, sample afresh against the
#' # particle posterior
#' corrected <- ctParticleCorrect(fit, draws = 'imis', particles = 500, substeps = 5)
#' }
#' @export
ctParticleCorrect <- function(fit, draws = c("reweight", "imis"), particles = 1000,
  substeps = 10, transition = c("exponential", "euler"), seed = 1,
  nsamples = NULL, finishsamples = NULL, scale = 1.5, target_ess = 100,
  nbatch = 200L, maxiter = 10L, correct_estimate = TRUE, cores = NULL,
  verbose = 0L) {

  draws <- match.arg(draws)
  transition <- match.arg(transition)
  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctParticleCorrect needs a fit from ctFit(..., backend = 'julia').", call. = FALSE)
  }
  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$laplace)) {
    stop("ctParticleCorrect is not available for intoverpop = 'laplace' or 'none' fits: ",
      "the particles carry random effects only as augmented states. Refit with the ",
      "default intoverpop.", call. = FALSE)
  }
  particles <- as.integer(particles)[1L]
  substeps <- as.integer(substeps)[1L]
  if (is.na(particles) || particles < 2L) stop("particles must be at least 2.", call. = FALSE)
  if (is.na(substeps) || substeps < 1L) stop("substeps must be at least 1.", call. = FALSE)
  posterior <- fit$estimate$rawposterior
  hasdraws <- !is.null(posterior) && length(dim(posterior)) == 2L && nrow(posterior) >= 2L
  if (identical(draws, "reweight") && !hasdraws) {
    stop("The fit carries no posterior draws to reweight. Run ctOptimUncertainty() ",
      "first, or use draws = 'imis'.", call. = FALSE)
  }
  covariance <- fit$estimate$cov
  if (identical(draws, "imis") && is.null(covariance)) {
    stop("draws = 'imis' needs the fit's covariance as its first proposal. Run ",
      "ctOptimUncertainty() first.", call. = FALSE)
  }
  proposal <- if (identical(draws, "reweight")) .ctParticleProposal(fit) else "imis"
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
  is_res <- NULL

  if (identical(draws, "reweight")) {
    n <- if (is.null(nsamples)) nrow(posterior) else as.integer(nsamples)[1L]
    if (is.na(n) || n < 2L) stop("nsamples must be at least 2.", call. = FALSE)
    if (identical(proposal, "normal")) {
      # Fresh draws, wider than the curvature says. The fit's own normal draws
      # are a cheap sample around the mode, not data, and a proposal no wider
      # than the Hessian cannot correct a posterior wider than it -- exact
      # weights on the fit's own draws left a two-parameter test posterior's
      # mean 0.24 standard errors short and its spread at 0.87 of the truth.
      propcov <- ctOptimSafeCov(as.matrix(covariance)) * scale^2
      theta <- ctOptimNormalDraws(est, propcov, n)
      colnames(theta) <- colnames(posterior)
      logq <- mvtnorm::dmvnorm(theta, mean = est, sigma = propcov, log = TRUE)
    } else {
      theta <- posterior
      if (n < nrow(theta)) theta <- theta[sort(sample.int(nrow(theta), n)), , drop = FALSE]
      logq <- NULL
    }
    # One seed per draw. The estimates have to be independent across draws for
    # the self-normalised weights to be consistent; a shared seed makes the
    # filter's Monte Carlo error a fixed function of the parameters and tilts
    # the whole posterior by it (0.45 standard errors on the same test fit).
    ev <- .ctParticleEvaluate(module, objective, theta, settings,
      seeds = settings$seed + seq_len(nrow(theta)) - 1L, verbose = verbose)
    # The exact log posterior at each draw: the filter's log posterior with the
    # filter's likelihood swapped for the particle one (prior and any other
    # term ride along in `posterior - filter`).
    target <- ev$posterior - ev$filter + ev$particle
    # Draws that already sit on the filter posterior have that as their density.
    if (is.null(logq)) logq <- ev$posterior
    logw <- target - logq
    logw[!is.finite(logw)] <- -Inf
    if (!any(is.finite(logw))) {
      stop("The particle filter could not evaluate any of the draws; nothing to ",
        "reweight.", call. = FALSE)
    }
    w <- exp(logw - max(logw))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)
    if (is.null(finishsamples)) finishsamples <- nrow(posterior)
    samples <- theta[.ctParticleResample(w, as.integer(finishsamples)[1L]), , drop = FALSE]
    newmean <- as.numeric(colSums(w * theta))
    newcov <- ctOptimSafeCov(stats::cov.wt(theta, wt = w, method = "ML")$cov)
    evaluations <- data.frame(particle = ev$particle, filter = ev$filter,
      difference = ev$particle - ev$filter, se = ev$se, weight = w)
    # Said plainly rather than left in a list nobody prints: a corrected
    # interval resting on a handful of effective draws is worse than the
    # uncorrected one, because it looks like it has been improved.
    if (ess < max(50, 0.1 * nrow(theta))) {
      warning("Reweighting left an effective sample size of ", round(ess, 1),
        " from ", nrow(theta), " draws. The corrected draws rest on few points; ",
        "treat the corrected intervals as indicative, or use draws = 'imis' to ",
        "sample against the particle posterior directly.", call. = FALSE)
    }
  } else {
    # Every evaluation is kept, so the record can show the particle and filter
    # likelihoods per proposal draw as it does under reweighting; `imis_is`
    # itself only keeps the target density.
    memo <- new.env(parent = emptyenv())
    memo$rows <- list()
    lp <- function(parm) {
      # Its own seed per proposal draw, for the reason given under 'reweight'.
      r <- .ctParticleEvaluate(module, objective, matrix(as.numeric(parm), nrow = 1L),
        settings, seeds = settings$seed + length(memo$rows))
      memo$rows[[length(memo$rows) + 1L]] <- r
      value <- r$posterior - r$filter + r$particle
      # A proposal draw the model cannot be evaluated at should carry no
      # weight, and `exp(-1e100 - c)` is zero. Same guard as ctLaplaceCorrect.
      if (!is.finite(value)) -1e100 else value
    }
    if (is.null(finishsamples)) finishsamples <- if (hasdraws) nrow(posterior) else 1000L
    if (verbose > 0) {
      message("Importance sampling against the particle posterior (", nbatch,
        " draws per iteration, one particle filter each)")
    }
    is_res <- imis_is(lp, mu_hat = est,
      Sigma_hat = ctOptimSafeCov(as.matrix(covariance)) * scale^2,
      cl = NA, n_batch = as.integer(nbatch), target_ess = target_ess,
      max_iter = as.integer(maxiter),
      # The proposal is already scaled above, so the sampler's own initial
      # scaling is left at one rather than compounding with it.
      scale_init = 1, tail_scale = 1.2, df = Inf,
      finishsamples = as.integer(finishsamples), verbose = verbose > 0,
      diag_plots = FALSE)
    samples <- is_res$theta
    if (is.null(samples) || !nrow(samples)) {
      stop("Importance sampling returned no usable draws against the particle ",
        "posterior. Widen the proposal (scale) or raise nbatch.", call. = FALSE)
    }
    theta <- is_res$full_theta
    w <- as.numeric(is_res$full_weights)
    ess <- if (is.null(is_res$ess)) NA_real_ else as.numeric(is_res$ess)[1L]
    newmean <- as.numeric(is_res$mean)
    newcov <- if (!is.null(is_res$covariance) && all(is.finite(is_res$covariance))) {
      ctOptimSafeCov(is_res$covariance)
    } else {
      ctOptimSafeCov(stats::cov(samples))
    }
    ev <- do.call(rbind, memo$rows)
    evaluations <- data.frame(particle = ev$particle, filter = ev$filter,
      difference = ev$particle - ev$filter, se = ev$se,
      weight = if (nrow(ev) == length(w)) w else NA_real_)
    if (is.finite(ess) && ess < target_ess / 2) {
      warning("Importance sampling reached an effective sample size of ",
        round(ess, 1), " against a target of ", target_ess,
        ". The corrected draws rest on few points, so treat the corrected ",
        "intervals as indicative; a wider scale or more iterations may help.",
        call. = FALSE)
    }
  }

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
    check = check, draws = draws,
    proposal = switch(proposal,
      normal = sprintf("fresh normal draws, scale %s on the fit's covariance", format(scale)),
      posterior = "the fit's draws, on the filter posterior",
      imis = "adaptive mixture (imis)"),
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
    # The draws on the fit now describe the particle posterior, so the record
    # of how they were made has to say so; `.ctParticleProposal()` reads it back
    # if the fit is corrected again, and summaries describe the fit from it.
    fit$uncertainty$draws <- "particle"
    if (!is.null(fit$uncertainty$settings)) {
      fit$uncertainty$settings$draws <- "particle"
      fit$uncertainty$settings$finishsamples <- nrow(samples)
    }
    if (!is.null(is_res)) {
      fit$uncertainty$proposal_cov <- ctOptimSafeCov(as.matrix(covariance)) * scale^2
      fit$uncertainty$imis <- is_res
    }
    fit$uncertainty$details$particle_correction <- list(draws = draws, ess = ess,
      particles = particles, substeps = substeps, transition = transition,
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
  if (is.finite(x$effective_fraction) && x$effective_fraction < 0.1) {
    cat("  Few draws carry the weight; draws = 'imis' samples against the particle",
      "posterior directly.\n")
  }
  if (!x$corrected_estimate) {
    cat("  estimate left at the original point; draws corrected only\n")
  }
  invisible(x)
}

# Which density the fit's draws came from, read off the fit's own record.
# `ctOptimUncertainty()` writes `fit$uncertainty$draws` for every method,
# `ctLaplaceCorrect()` and `ctParticleCorrect()` rewrite it when they replace
# the draws, and `ctSample()` marks a sampled fit with `fit$sample`.
.ctParticleProposal <- function(fit) {
  if (!is.null(fit$sample) || identical(fit$uncertainty$settings$method, "sampling")) {
    return("posterior")
  }
  kind <- fit$uncertainty$draws
  if (is.null(kind)) kind <- fit$uncertainty$settings$draws
  if (is.null(kind)) {
    stop("The fit does not record how its posterior draws were produced ",
      "(fit$uncertainty$draws), so their density is unknown. Recompute them with ",
      "ctOptimUncertainty(), or use draws = 'imis'.", call. = FALSE)
  }
  switch(as.character(kind)[1L],
    normal = {
      if (is.null(fit$estimate$cov)) {
        stop("The fit's draws are normal but it carries no covariance to describe ",
          "them with. Recompute them with ctOptimUncertainty().", call. = FALSE)
      }
      "normal"
    },
    imis = "posterior", empirical = "posterior", particle = "posterior",
    stop("Draws of kind '", kind, "' are not recognised here. Recompute them with ",
      "ctOptimUncertainty(), or use draws = 'imis'.", call. = FALSE))
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
    res <- JuliaConnectoR::juliaGet(module$ctsem_particle_batch(objective, values,
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

# Systematic resampling by normalised weights, as `imis_is` does it.
.ctParticleResample <- function(w, n) {
  cs <- cumsum(w / sum(w))
  idx <- findInterval((stats::runif(1) + seq_len(n) - 1) / n, cs) + 1L
  pmin(idx, length(w))
}
