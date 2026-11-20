#' Check a julia-backend fit against a particle-filter likelihood
#'
#' The Kalman filter ctsem fits with approximates the distribution of the
#' latent states by a Gaussian. It linearises the drift over each prediction
#' step, and approximates non-Gaussian measurement updates by Gaussian ones.
#' For a linear model with Gaussian indicators it is exact, and otherwise its
#' error cannot be seen from the fit alone. This function computes the same
#' marginal log likelihood at the fit's estimates with a bootstrap particle
#' filter, which makes none of these approximations, and reports the
#' difference, in total and per row.
#'
#' Each particle moves through the process with drift, diffusion and every
#' state dependent parameter evaluated at its own state, is weighted by the
#' density of the row's observations given that state (Gaussian, binary,
#' ordinal, count and censored indicators are all handled), and is resampled
#' when the weights become uneven. The estimate is unbiased, and its Monte
#' Carlo error falls with the square root of \code{particles}.
#'
#' It is a check at fixed parameters, and does not estimate them. The
#' comparison can be trusted when \code{se} is well below the difference that
#' matters (raise \code{particles}, or set \code{replicates} for an empirical
#' standard error), and doubling \code{substeps} does not change the answer.
#'
#' @param fit A fit from \code{ctFit(..., backend = 'julia')} with the default
#'   \code{intoverpop}. Random effects carried as extra states are integrated
#'   by the particles. Fits using the Laplace approximation are refused.
#' @param particles Particles per subject.
#' @param substeps Transition steps per observation interval, independent of
#'   the fit's own settings, but never coarser than the fit's
#'   \code{maxtimestep} or automatic step size.
#' @param transition \code{'exponential'} uses the filter's own
#'   discretisation at each particle's state, which is exact for the linear
#'   parts of the model, so the comparison shows the error of the filter's
#'   Gaussian approximation. \code{'euler'} is plain Euler-Maruyama, which
#'   shares no approximation with the filter but needs several times more
#'   substeps.
#' @param seed Seed for every draw, so the result is reproducible.
#' @param replicates Number of independent runs (with seeds \code{seed},
#'   \code{seed + 1}, ...). With one run, an approximate standard error is
#'   computed from the effective sample sizes, and with more, the standard
#'   deviation across runs is reported.
#' @param cores Engine threads. Each subject has its own random stream, so
#'   the result does not depend on the number of cores.
#' @return A list: \code{loglik}, the particle estimate (mean over replicates);
#'   \code{se}; \code{fit_loglik}, the filter's log likelihood from the fit;
#'   \code{difference}, particle minus filter; \code{rows}, a data frame with
#'   the subject, time, particle and filter log likelihood of every row and
#'   their difference; \code{replicates}, the individual estimates;
#'   \code{ess_min}, the smallest effective sample size any row saw; and the
#'   settings.
#' @seealso \code{\link{ctFitUncertainty}} with \code{method = 'is'} and
#'   \code{filter = 'particle'} applies the correction this
#'   measures, to the fit's posterior draws.
#' @examples
#' \donttest{
#' # Needs a working julia backend, so the example is inert where Julia is
#' # absent -- including on CRAN, whose check machines have none.
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
#'     LAMBDA = matrix(1))
#'   model$pars$indvarying <- FALSE
#'   fit <- ctFit(ctstantestdat, model, backend = 'julia', cores = 1)
#'   # A linear Gaussian model, where the filter is exact: the difference
#'   # should be within a few standard errors of zero.
#'   check <- ctParticleLik(fit, particles = 1000, substeps = 5)
#'   check[c('loglik', 'se', 'fit_loglik', 'difference')]
#'   # Julia keeps running between calls, holding what it has compiled;
#'   # this ends it.
#'   ctJuliaKill("session")
#' }
#' }
#' @export
ctParticleLik <- function(fit, particles = 2000, substeps = 20,
  transition = c("exponential", "euler"), seed = 1, replicates = 1, cores = getOption("mc.cores", 2L)) {
  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctParticleLik needs a fit from ctFit(..., backend = 'julia').", call. = FALSE)
  }
  transition <- match.arg(transition)
  if (!is.null(cores)) {
    previous <- .ctBackendSetMaxChunks(max(1L, cores))
    on.exit(.ctBackendRestoreMaxChunks(previous), add = TRUE)
  }
  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$laplace)) {
    stop("ctParticleLik is not available for intoverpop = 'laplace' or 'none' fits: ",
      "the particles carry random effects only as augmented states. Refit with ",
      "intoverpop = 'augmented'.", call. = FALSE)
  }
  particles <- as.integer(particles)[1L]
  substeps <- as.integer(substeps)[1L]
  replicates <- max(1L, as.integer(replicates)[1L])
  if (is.na(particles) || particles < 2L) stop("particles must be at least 2.", call. = FALSE)
  if (is.na(substeps) || substeps < 1L) stop("substeps must be at least 1.", call. = FALSE)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  values <- .ctJuliaNumericVector(as.numeric(fit$estimate$raw))
  runs <- lapply(seq_len(replicates), function(r) {
    .ctJuliaGet(module$ctsem_particle_loglik(objective, values,
      particles = particles, substeps = substeps, transition = transition,
      seed = as.integer(seed)[1L] + r - 1L))
  })
  estimates <- vapply(runs, function(x) as.numeric(x$loglik), numeric(1))
  rowmat <- vapply(runs, function(x) as.numeric(x$row_loglik), numeric(length(spec$times)))
  rowmean <- if (is.matrix(rowmat)) rowMeans(rowmat) else as.numeric(rowmat)
  se <- if (replicates > 1L) stats::sd(estimates) / sqrt(replicates) else as.numeric(runs[[1L]]$se)
  # The filter's own per-row increments, at the same values, for the same rows.
  # `juliaCall` hands a plain numeric vector straight back; only the NamedTuple
  # above needs `juliaGet`.
  filter_rows <- as.numeric(.ctJuliaCall(
    "ContinuousTimeSEM._ctsem_row_loglikelihood", objective, values))
  starts <- as.integer(spec$subject_starts)
  counts <- diff(c(starts, length(spec$times) + 1L))
  rows <- data.frame(subject = rep(seq_along(starts), counts), time = as.numeric(spec$times),
    particle = rowmean, filter = filter_rows, difference = rowmean - filter_rows)
  fit_loglik <- as.numeric(fit$estimate$loglik)
  # The one place in this function that says what the number means rather than
  # what it is. `$difference` is the same quantity `ctLaplaceCheck()` calls
  # `$gap` and the curvature certification calls `$gap`; `$verdict` is the
  # shape all three now share. See R/ctFitGap.R.
  #
  # This is the only one of the three whose reference is stochastic, so it is
  # the only one whose verdict can say whether the gap is distinguishable from
  # its own noise -- which `$se` allowed and nothing previously asked.
  verdict <- .ctFitGap("particle",
    gap = mean(estimates) - fit_loglik, gap_se = se,
    nsubjects = .ctFitNsubjects(fit),
    remedy = paste0("A resolved positive gap means the filter's Gaussian ",
      "assumption is costing likelihood; ctFitUncertainty(fit, 'is', filter = ",
      "'particle') samples the posterior against the ",
      "particle filter."))
  list(verdict = verdict,
    loglik = mean(estimates), se = se, fit_loglik = fit_loglik,
    difference = mean(estimates) - fit_loglik, rows = rows, replicates = estimates,
    ess_min = min(vapply(runs, function(x) as.numeric(x$ess_min), numeric(1))),
    particles = particles, substeps = substeps, transition = transition)
}
