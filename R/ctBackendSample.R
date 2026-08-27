# Hamiltonian sampling of a julia backend fit.
#
# `intoverpop='laplace'` approximates each unit's integral by a Gaussian at its
# mode. That is exact when the integrand is Gaussian in the random effects and
# otherwise wrong by an amount that grows with the population scale, which tilts
# the profile and shrinks the scale estimate -- `ctLaplaceCheck()` measures that
# error and corrects it to first order. `ctSample()` removes it instead, by
# sampling the joint posterior over population parameters *and* random effects
# with no Gaussian assumption anywhere.
#
# It takes a fitted object rather than a model and data, and that is not merely
# convenience. The fit supplies the starting point *and* the metric: the engine
# reads the sampler's initial mass matrix off the Laplace curvature, block by
# block, so the chain begins as well conditioned as the approximation can make
# it and warmup refines rather than discovers. Sampling from scratch would work
# and would be substantially slower.

#' Sample the posterior of a julia backend fit
#'
#' Draws from the joint posterior over population parameters and random effects
#' by Hamiltonian Monte Carlo (the No-U-Turn sampler), starting from a fit made
#' with \code{intoverpop='laplace'} and using that fit's curvature as the
#' sampler's metric.
#'
#' This is the exact counterpart of the Laplace approximation rather than a
#' replacement for it: where \code{\link{ctLaplaceCheck}} measures how wrong the
#' Gaussian approximation is, this does not make it. The cost is time --
#' thousands of gradient evaluations rather than hundreds -- and the return is a
#' posterior rather than a point estimate with a normal approximation around it.
#'
#' The result is a \code{ctJuliaFit} carrying \code{estimate$rawposterior}, so
#' \code{\link{summary}}, \code{\link{ctExtract}}, \code{ctKalman} and the
#' system-matrix helpers all read it the way they read an optimised fit's
#' normal-approximation draws.
#'
#' @param fit A \code{ctJuliaFit} made with \code{intoverpop='laplace'}.
#' @param chains Number of chains. Run concurrently when the Julia session has
#'   at least that many threads; see \code{\link{ctJuliaSetup}}.
#' @param warmup Warmup iterations per chain, used to adapt the step size and
#'   refine the metric, and discarded.
#' @param draws Retained iterations per chain.
#' @param cores Ceiling on the engine's parallelism. With several chains they
#'   take a thread each; with one chain the subject loop is split instead.
#' @param saveEffects Return every draw of every random effect, not just their
#'   posterior mean and standard deviation. Off by default because the draws are
#'   \code{nsubjects * neffects * chains * draws} numbers and the R-to-Julia
#'   bridge moves about 1 MB/s -- for a hundred subjects that transfer takes
#'   longer than many fits do.
#' @param seed Random seed; each chain uses \code{seed + chain}.
#' @param control A list of sampler settings: \code{maxdepth} (default 10),
#'   \code{target_accept} (0.8), \code{adapt_metric} (TRUE), \code{init_scale}
#'   (1), \code{maxdelta} (1000).
#' @param verbose Print the sampler's configuration before it starts.
#'
#' @return The fit, with \code{estimate$rawposterior} holding the draws and
#'   \code{$sample} holding the diagnostics: split R-hat and effective sample
#'   size per parameter, divergences, tree depths, step sizes and E-BFMI.
#' @export
ctSample <- function(fit, chains = 4L, warmup = 500L, draws = 500L, cores = 1L,
  saveEffects = FALSE, seed = 20260828L, control = list(), verbose = FALSE) {

  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctSample applies to fits made with ctFit(backend='julia').", call. = FALSE)
  }
  if (is.null(fit$model_spec$laplace)) {
    stop("ctSample needs a fit made with intoverpop='laplace'. The augmented ",
      "route carries the random effects in the state, so there is no separate ",
      "posterior over them to sample.", call. = FALSE)
  }
  chains <- max(1L, as.integer(chains)[1L])
  warmup <- max(0L, as.integer(warmup)[1L])
  draws <- max(1L, as.integer(draws)[1L])
  cores <- max(1L, as.integer(cores)[1L])

  module <- .ctJuliaModule(fit$model_spec$project)
  objective <- .ctJuliaObjective(fit)
  estimate <- as.numeric(fit$estimate$raw)
  npar <- length(estimate)

  # Chains are the parallel axis, and they can only be concurrent if the session
  # was started with threads for them. Said once, here, because the alternative
  # is a user concluding the sampler is slow when it is running four chains on
  # one thread.
  threads <- tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
    error = function(e) NA_integer_)
  if (!is.na(threads) && chains > 1L && threads < chains) {
    message("The Julia session has ", threads, " thread(s) and ", chains,
      " chains were asked for, so they will run one after another. ",
      "ctJuliaSetup(threads = ", chains, ", force = TRUE) before fitting ",
      "runs them together.")
  }

  # The fit's Hessian, when it has one: the sampler would otherwise recompute
  # it to build the metric, at 2n gradient evaluations it need not spend.
  hessian <- fit$uncertainty$hessian
  arguments <- list(objective, .ctJuliaNumericVector(estimate),
    npar = as.integer(npar), nchains = chains, nwarmup = warmup,
    ndraws = draws, seed = as.integer(seed)[1L],
    save_effects = isTRUE(saveEffects), verbose = isTRUE(verbose),
    maxdepth = as.integer(.ctJuliaOr(control$maxdepth, 10L)),
    target_accept = as.numeric(.ctJuliaOr(control$target_accept, 0.8)),
    maxdelta = as.numeric(.ctJuliaOr(control$maxdelta, 1000)),
    init_scale = as.numeric(.ctJuliaOr(control$init_scale, 1)),
    adapt_metric = isTRUE(.ctJuliaOr(control$adapt_metric, TRUE)))
  if (!is.null(hessian)) {
    arguments$hessian <- JuliaConnectoR::juliaPut(as.matrix(hessian))
  }

  result <- .ctBackendWithMaxChunks(cores,
    JuliaConnectoR::juliaGet(do.call(module$ctsem_sample, arguments)))

  # The engine returns draws parameter-major; every R-side consumer wants them
  # draw-major, which is also the shape `ctOptimUncertainty` leaves behind.
  #
  # The row count is `ndim` when the effects were saved and `npar` when they
  # were not, so it is read off the result rather than assumed -- reshaping an
  # effects-carrying matrix to `npar` rows would silently interleave parameters
  # and effects into plausible-looking nonsense.
  kept <- if (isTRUE(saveEffects)) as.integer(result$ndim) else as.integer(result$npar)
  raw <- matrix(as.numeric(result$draws), nrow = kept)
  posterior <- t(raw[seq_len(npar), , drop = FALSE])
  colnames(posterior) <- .ctBackendRawParameterNames(fit, npar)

  out <- fit
  out$estimate$rawposterior <- posterior
  # The posterior mean, not the Laplace mode, is now the point estimate: it is
  # what the draws describe, and leaving `raw` at the mode would make ctKalman()
  # and the system matrices report a different fit from the one summarised.
  out$estimate$laplace_raw <- estimate
  out$estimate$raw <- as.numeric(colMeans(posterior))
  out$estimate$cov <- stats::cov(posterior)
  out$estimate$se <- sqrt(diag(out$estimate$cov))
  out$uncertainty <- list(method = "sampling", hessian = hessian,
    settings = list(chains = chains, warmup = warmup, draws = draws))

  out$sample <- list(
    chains = chains, warmup = warmup, draws = draws,
    rhat = stats::setNames(as.numeric(result$rhat)[seq_len(npar)], colnames(posterior)),
    ess = stats::setNames(as.numeric(result$ess)[seq_len(npar)], colnames(posterior)),
    divergent = as.integer(result$ndivergent),
    warmup_divergent = as.integer(result$warmup_divergent),
    saturated = as.integer(result$nsaturated),
    max_depth = as.integer(result$max_depth),
    stepsize = as.numeric(result$stepsize),
    ebfmi = as.numeric(result$ebfmi),
    accept = as.numeric(result$accept),
    depth = as.integer(result$depth),
    energy = as.numeric(result$energy),
    effect_mean = as.numeric(result$effect_mean),
    effect_sd = as.numeric(result$effect_sd),
    start = estimate)
  if (isTRUE(saveEffects)) {
    out$sample$effects <- t(raw[-seq_len(npar), , drop = FALSE])
  }
  class(out$sample) <- "ctSampleDiagnostics"

  # Constrained draws describe the *new* draws, so the cached ones are stale.
  out$transformedpars <- NULL
  out$transformedpars <- .ctBackendConstrained(out)
  out$kalman <- suppressMessages(ctKalmanArray(out, pointest = TRUE))

  .ctSampleWarn(out$sample)
  out
}

# The three failures worth interrupting for, in the order a user should read
# them. Deliberately not silent: a divergent transition means the sampler could
# not follow the geometry there, and a posterior summarised over draws it could
# not reach is wrong in a way no amount of averaging fixes.
.ctSampleWarn <- function(diagnostics) {
  total <- diagnostics$chains * diagnostics$draws
  if (diagnostics$divergent > 0L) {
    warning(diagnostics$divergent, " of ", total, " transitions diverged. The ",
      "sampler could not follow the posterior's geometry there, so these draws ",
      "under-represent whatever it could not reach -- most often a population ",
      "standard deviation near zero. Raising control$target_accept towards ",
      "0.95 shortens the steps and often clears it.", call. = FALSE)
  }
  worst <- suppressWarnings(max(diagnostics$rhat, na.rm = TRUE))
  if (is.finite(worst) && worst > 1.01) {
    warning("Largest R-hat is ", signif(worst, 4), ". The chains have not ",
      "agreed on the same distribution, so the draws are not yet a posterior. ",
      "Run longer, and see fit$sample$rhat.", call. = FALSE)
  }
  fewest <- suppressWarnings(min(diagnostics$ess, na.rm = TRUE))
  if (is.finite(fewest) && fewest < 100) {
    warning("Smallest effective sample size is ", round(fewest), ", from ",
      total, " draws. Interval estimates from this few are unreliable; see ",
      "fit$sample$ess.", call. = FALSE)
  }
  if (diagnostics$saturated > 0L) {
    message(diagnostics$saturated, " of ", total, " transitions hit the maximum ",
      "tree depth of ", diagnostics$max_depth, ". That costs efficiency rather ",
      "than correctness; control$maxdepth raises it.")
  }
  invisible(diagnostics)
}

#' @export
print.ctSampleDiagnostics <- function(x, ...) {
  total <- x$chains * x$draws
  cat("ctsem Hamiltonian sample\n")
  cat("  ", x$chains, " chains x ", x$draws, " draws (", x$warmup,
    " warmup discarded)\n", sep = "")
  cat("  divergent: ", x$divergent, " of ", total,
    "   max tree depth reached: ", x$saturated, "\n", sep = "")
  cat("  step size: ", paste(signif(x$stepsize, 3), collapse = ", "),
    "\n", sep = "")
  cat("  E-BFMI:    ", paste(signif(x$ebfmi, 3), collapse = ", "),
    if (any(x$ebfmi < 0.3, na.rm = TRUE)) "  (below 0.3 suggests a funnel)" else "",
    "\n", sep = "")
  worst <- order(-x$rhat)[seq_len(min(5L, length(x$rhat)))]
  cat("  worst R-hat and effective sample size:\n")
  print(data.frame(parameter = names(x$rhat)[worst],
    rhat = round(x$rhat[worst], 4), ess = round(x$ess[worst])),
    row.names = FALSE)
  invisible(x)
}
