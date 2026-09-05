#' Check a julia-backend fit against a particle-filter likelihood
#'
#' The filter ctsem fits with is an assumed-density approximation: it carries a
#' Gaussian for the latent states, linearises the drift over each prediction
#' step, and projects non-Gaussian measurement updates back onto a Gaussian.
#' For a linear model with Gaussian indicators it is exact; for anything else
#' the size of its error is invisible from inside it. This function estimates
#' the same marginal log likelihood at the fit's estimates with a bootstrap
#' particle filter that makes none of those approximations, and reports the
#' difference, in total and per row.
#'
#' Each particle steps through the process with drift, diffusion and every
#' state-dependent parameter evaluated at its own state, is weighted by the
#' conditional density of the row's observations given that state (Gaussian,
#' binary, ordinal, count and censored indicators are all handled), and is
#' resampled when the weights degenerate. The estimate is unbiased with Monte
#' Carlo error falling as one over the square root of \code{particles}.
#'
#' It is a diagnostic at fixed parameters, not an estimator. Two things decide
#' whether to trust the comparison: \code{se} should be well below the
#' difference you care about (raise \code{particles}, or set \code{replicates}
#' for an empirical standard error), and doubling \code{substeps} should not
#' move the answer.
#'
#' @param fit A fit from \code{ctFit(..., backend = 'julia')} with the default
#'   \code{intoverpop}; random effects carried as augmented states are
#'   integrated by the particles. Laplace-route fits are refused.
#' @param particles Particles per subject.
#' @param substeps Transition steps per observation interval. This is the
#'   reference mesh, independent of anything the fit used, and floored by the
#'   fit's own \code{maxtimestep} or automatic mesh.
#' @param transition \code{'exponential'} applies the filter's own
#'   discretisation at each particle's state: exact for the linear parts of the
#'   model, second order otherwise, so the comparison isolates the filter's
#'   Gaussian assumption. \code{'euler'} is plain Euler-Maruyama, which shares no
#'   approximation with the filter at all but needs several times more substeps.
#' @param seed Fixes every draw; the result is reproducible given it.
#' @param replicates Independent runs (with seeds \code{seed}, \code{seed + 1},
#'   ...). One run reports an approximate standard error from the effective
#'   sample sizes; more than one reports the standard deviation across runs.
#' @param cores Engine threads. The particle filter splits the subjects across
#'   them, each subject with its own random stream, so the result does not
#'   depend on the count.
#' @return A list: \code{loglik}, the particle estimate (mean over replicates);
#'   \code{se}; \code{fit_loglik}, the filter's log likelihood from the fit;
#'   \code{difference}, particle minus filter; \code{rows}, a data frame with
#'   the subject, time, particle and filter log likelihood of every row and
#'   their difference; \code{replicates}, the individual estimates;
#'   \code{ess_min}, the smallest effective sample size any row saw; and the
#'   settings.
#' @seealso \code{\link{ctParticleCorrect}} applies the correction this
#'   measures, to the fit's posterior draws.
#' @export
ctParticleLik <- function(fit, particles = 2000, substeps = 20,
  transition = c("exponential", "euler"), seed = 1, replicates = 1, cores = NULL) {
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
      "the particles carry random effects only as augmented states. Refit with the ",
      "default intoverpop.", call. = FALSE)
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
    JuliaConnectoR::juliaGet(module$ctsem_particle_loglik(objective, values,
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
  filter_rows <- as.numeric(JuliaConnectoR::juliaCall(
    "ContinuousTimeSEM._ctsem_row_loglikelihood", objective, values))
  starts <- as.integer(spec$subject_starts)
  counts <- diff(c(starts, length(spec$times) + 1L))
  rows <- data.frame(subject = rep(seq_along(starts), counts), time = as.numeric(spec$times),
    particle = rowmean, filter = filter_rows, difference = rowmean - filter_rows)
  fit_loglik <- as.numeric(fit$estimate$loglik)
  list(loglik = mean(estimates), se = se, fit_loglik = fit_loglik,
    difference = mean(estimates) - fit_loglik, rows = rows, replicates = estimates,
    ess_min = min(vapply(runs, function(x) as.numeric(x$ess_min), numeric(1))),
    particles = particles, substeps = substeps, transition = transition)
}
