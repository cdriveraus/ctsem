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

# Which (subject, parameter) each entry of the flat effect vector belongs to.
#
# The engine lays the effects out unit by unit, and a unit's own layout is its
# block tree -- one block per (level, group) it contains. With a single level
# that degenerates to one unit per subject and `k` effects each, in subject
# order, which is determinable from the R side alone.
#
# With more than one level it is not: the blocks interleave a group's own
# effects with its members', and reconstructing that here would mean
# reimplementing `_laplace_build_units` in R and keeping the two in step.
# Returning nothing is better than returning a plausible mislabelling, which
# would attach the wrong subject's name to a number and never announce itself.
#' @keywords internal
.ctBackendEffectIndex <- function(fit) {
  laplace <- fit$model_spec$laplace
  if (is.null(laplace) || is.null(laplace$levels)) return(NULL)
  if (length(laplace$levels) != 1L) return(NULL)
  level <- laplace$levels[[1L]]
  parameters <- as.character(level$param)
  if (!length(parameters)) return(NULL)
  subjects <- fit$model_spec$subject_starts
  nsubjects <- if (is.null(subjects)) 0L else length(subjects)
  if (nsubjects < 1L) return(NULL)
  ids <- .ctBackendSubjectIds(fit, nsubjects)
  index <- expand.grid(parameter = parameters, subject = seq_len(nsubjects),
    stringsAsFactors = FALSE)
  if (is.null(ids)) {
    index$label <- paste0(index$parameter, "_subject", index$subject)
    return(index[, c("subject", "parameter", "label")])
  }
  index$id <- ids[index$subject]
  index$label <- paste0(index$parameter, "_", index$id)
  index[, c("subject", "id", "parameter", "label")]
}

# The user's own identifiers, when the fit kept them, and the internal index
# otherwise. A label of "mmean_7" is only useful if 7 is the id the user knows.
#' @keywords internal
.ctBackendSubjectIds <- function(fit, nsubjects) {
  # In first-appearance order, which is the order the engine numbers subjects.
  if (!is.null(fit$data) && !is.null(fit$data$id)) {
    original <- unique(fit$data$id)
    if (length(original) == nsubjects) return(as.character(original))
  }
  ids <- fit$model_spec$subject_ids
  if (!is.null(ids) && length(ids) == nsubjects) return(as.character(ids))
  # No map found. The internal index is still a correct label, and calling it
  # `id` when it is not the user's id would be worse than not having one.
  NULL
}

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
#' The reported \code{estimate$se} of an optimised fit is the curvature of the
#' approximated marginal posterior at its mode, so it describes a normal
#' approximation rather than the posterior itself. That approximation is what a
#' sample replaces, and the two differ most where the posterior is skewed --
#' variance-like parameters at modest subject counts. Measured on a model whose
#' Laplace integral is exact, the sampled standard deviations were 1.1 to 2.3
#' times the reported standard errors at forty subjects and within 13\% of them
#' at two hundred.
#'
#' @section Diagnostics:
#' Divergent transitions, R-hat above 1.01 and effective sample sizes below 100
#' warn rather than pass quietly. A divergence means the sampler could not
#' follow the posterior's geometry somewhere, most often a population standard
#' deviation near zero, and draws that miss such a region are wrong in a way
#' averaging does not fix. \code{fit$sample} carries the per-parameter R-hat and
#' effective sample size, the per-draw acceptance statistic, tree depth and
#' energy, and the per-chain step size and E-BFMI. An E-BFMI below about 0.3
#' indicates a funnel the metric could not straighten.
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
#'   \code{target_accept} (0.8), \code{adapt_metric} (TRUE),
#'   \code{adapt_effects} (FALSE), \code{init_scale} (1), \code{maxdelta}
#'   (1000). \code{adapt_effects} controls whether warmup re-estimates the
#'   random-effect blocks of the metric as well as the population block; they
#'   start from a conditional covariance that is exact for a linear model, so
#'   replacing one with an estimate from a few hundred draws can add more noise
#'   than it removes.
#' @param verbose Print the sampler's configuration before it starts, and
#'   report progress while it runs. Progress overwrites a single line where the
#'   output is going to a console and prints occasional separate lines where it
#'   is not; set \code{options(ctsem.progress.overwrite = FALSE)} if that
#'   detection is wrong for your front end, or \code{TRUE} to force it on.
#'
#' @return The fit, with \code{estimate$rawposterior} holding the draws and
#'   \code{$sample} holding the diagnostics: split R-hat and effective sample
#'   size per parameter, divergences, tree depths, step sizes and E-BFMI.
#'
#' @seealso \code{\link{ctLaplaceCheck}} measures the Laplace approximation's
#'   error and corrects it to first order, at a small fraction of the cost;
#'   \code{\link{ctFit}} for the fit this starts from, and
#'   \code{\link{ctJuliaSetup}} for the thread count that decides whether
#'   chains run concurrently.
#'
#' @examples
#' \dontrun{
#' data <- ctstantestdat
#' model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
#'   LAMBDA = matrix(1))
#' model$pars$indvarying <- model$pars$matrix %in% 'MANIFESTMEANS'
#'
#' # Four threads so the four chains run together rather than in turn.
#' ctJuliaSetup(threads = 4, force = TRUE)
#' fit <- ctFit(data, model, backend = 'julia', intoverpop = 'laplace')
#'
#' sampled <- ctSample(fit, chains = 4, warmup = 1000, draws = 1000, cores = 4)
#' sampled$sample                  # convergence and geometry diagnostics
#' summary(sampled)                # reads the draws, not a normal approximation
#'
#' # How far the Laplace approximation itself is from exact, for comparison.
#' ctLaplaceCheck(fit)
#' }
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
    adapt_metric = isTRUE(.ctJuliaOr(control$adapt_metric, TRUE)),
    adapt_effects = isTRUE(.ctJuliaOr(control$adapt_effects, FALSE)))
  if (!is.null(hessian)) {
    arguments$hessian <- JuliaConnectoR::juliaPut(as.matrix(hessian))
  }

  arguments$progress_overwrite <- .ctProgressOverwrite(verbose)
  result <- .ctBackendWithMaxChunks(cores,
    JuliaConnectoR::juliaGet(do.call(module$ctsem_sample, arguments)))
  .ctBackendSampleAssemble(fit, result, npar, saveEffects, chains, warmup,
    as.integer(result$ndraws), hessian, estimate)
}

# Turn an engine sample result into a fit object.
#
# Shared by `ctSample()` and by `ctFit(optimize = FALSE)`, which differ only in
# which engine entry point produced the draws: the joint sampler returns
# population parameters and effects, the marginal ones return population
# parameters with the effects already integrated out. Everything after that --
# where the draws go, what becomes the point estimate, which diagnostics warn --
# is the same, and was worth having in one place rather than two that drift.
#' @keywords internal
.ctBackendSampleAssemble <- function(fit, result, npar, saveEffects, chains,
  warmup, draws, hessian, startvalues) {

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
  # The posterior mean, not the mode, is now the point estimate: it is what the
  # draws describe, and leaving `raw` at the mode would make ctKalman() and the
  # system matrices report a different fit from the one summarised.
  out$estimate$laplace_raw <- as.numeric(startvalues)
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
    marginal = identical(as.integer(result$ndim), as.integer(result$npar)),
    start = as.numeric(startvalues))
  if (length(out$sample$effect_mean)) {
    out$sample$effectIndex <- .ctBackendEffectIndex(fit)
    labels <- out$sample$effectIndex$label
    if (length(labels) == length(out$sample$effect_mean)) {
      names(out$sample$effect_mean) <- labels
      names(out$sample$effect_sd) <- labels
    }
  }
  if (isTRUE(saveEffects) && kept > npar) {
    out$sample$effects <- t(raw[-seq_len(npar), , drop = FALSE])
    if (!is.null(out$sample$effectIndex) &&
        ncol(out$sample$effects) == nrow(out$sample$effectIndex)) {
      colnames(out$sample$effects) <- out$sample$effectIndex$label
    }
  }
  class(out$sample) <- "ctSampleDiagnostics"

  # Constrained draws describe the *new* draws, so the cached ones are stale.
  out$transformedpars <- NULL
  out$transformedpars <- .ctBackendConstrained(out)
  out$priorerrors <- .ctBackendPriorErrors(out)

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

# `ctFit(backend='julia', optimize=FALSE)`: fit by sampling rather than by
# maximising.
#
# Which sampler depends on `intoverpop`, and the three are genuinely different
# targets rather than three settings of one:
#
#   'none'       the joint posterior over population parameters *and* every
#                subject's random effects. Exact whatever the model, and
#                `npar + sum_U dim(u_U)` dimensions, so the dimension grows
#                with the subject count.
#   'laplace'    the population parameters, with the effects integrated by the
#                Laplace approximation. `npar` dimensions.
#   'augmented'  the population parameters, with the effects integrated by the
#                filter itself. `npar` dimensions, and the cheapest gradient of
#                the three.
#
# Measured on a model where all three are exact, 4 chains of 1000 warmup and
# 1000 draws: at 40 subjects the joint sampler is the more efficient (0.28
# effective draws per second against the Laplace marginal's 0.23), and at 200
# subjects the marginal is (1.19 against 0.59). Dimension is why -- the joint
# target grows from 47 coordinates to 207 while the marginal stays at 7 -- so
# the crossover moves with the subject count and neither is right everywhere.
#
# An optimisation runs first regardless, and is not merely a convenience: the
# sampler reads its metric from the fit's curvature, which is the difference
# between a chain that starts well conditioned and one that spends its warmup
# discovering what the model could have told it.
#
# What is optimised is always an *integrated* objective, never the joint one,
# and that distinction matters more than it looks. Maximising over individual
# random effects is not a defensible thing to do: the joint density of
# parameters and effects has no interior maximum in the scale direction --
# drive a population standard deviation to zero with the effects at their
# centre and it diverges -- so a joint mode is an artefact of where the
# optimiser stopped rather than a location worth starting from.
#
# For `intoverpop='none'`, which samples the effects, the metric therefore
# comes from the *Laplace* fit of the same specification. That works because
# 'none' and 'laplace' prepare identically -- same parameters, same ordering,
# same number of them -- and differ only in whether the effects are integrated
# or sampled afterwards. So the population block of the metric is a Laplace
# outer curvature, the effect blocks are the per-unit conditional curvatures,
# and the parameter vector they index is the one being sampled. The integration
# approach used for the metric is deliberately not the one being sampled, and
# it does not have to be.
#' @keywords internal
.ctJuliaSampleFit <- function(model_spec, datalong, model, inits, cores,
  backendcontrol, optimcontrol, chains, iter, control, priors, intoverpop,
  gradient, verbose) {

  npar <- max(c(model_spec$parameter_table$parnumber, model_spec$laplace$npar,
    model_spec$ti_effects$coefficient), na.rm = TRUE)
  start <- .ctJuliaInitialValues(npar, inits)

  # Stan's vocabulary, because these are Stan's arguments: `iter` counts warmup
  # and sampling together and warmup is half of it unless said otherwise.
  warmup <- as.integer(.ctJuliaOr(control$warmup, max(1L, floor(iter / 2))))
  draws <- max(1L, as.integer(iter) - warmup)
  maxdepth <- as.integer(.ctJuliaOr(control$max_treedepth, 10L))
  target <- as.numeric(.ctJuliaOr(control$adapt_delta, 0.8))
  seed <- as.integer(.ctJuliaOr(control$seed, 20260828L))
  saveEffects <- isTRUE(optimcontrol$saveEffects)

  if (verbose > 0L) {
    message("Optimising the ",
      if (identical(intoverpop, "augmented")) "augmented" else "Laplace",
      " objective to initialise the sampler and its metric",
      if (identical(intoverpop, "none"))
        " (the effects are integrated for this step, and sampled after it)"
      else "", ".")
  }
  optimised <- .ctJuliaOptimise(model_spec, start, backendcontrol = backendcontrol,
    gradient = gradient, cores = cores, verbose = verbose)
  estimate <- as.numeric(optimised$minimizer)

  spec <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(model_spec$project)
  objective <- .ctJuliaObjective(spec)
  hessian <- try(.ctBackendHessian(list(model_spec = model_spec,
    estimate = list(raw = estimate), backend = "julia"), estimate,
    verbose = verbose), silent = TRUE)
  if (inherits(hessian, "try-error")) hessian <- NULL

  arguments <- list(objective, .ctJuliaNumericVector(estimate),
    nchains = as.integer(chains), nwarmup = warmup, ndraws = draws,
    maxdepth = maxdepth, target_accept = target, seed = seed,
    verbose = verbose > 0L,
    progress_overwrite = .ctProgressOverwrite(verbose))
  # Sampling targets, when asked for. Left at zero the sampler takes exactly the
  # draws it was told to; set, it keeps going until the effective sample size is
  # there or the budget runs out, which is usually what a user wanted from a
  # draw count they had to guess at.
  if (!is.null(control$minEss)) arguments$min_ess <- as.numeric(control$minEss)
  if (!is.null(control$meanEss)) arguments$mean_ess <- as.numeric(control$meanEss)
  if (!is.null(control$maxDraws)) arguments$max_draws <- as.integer(control$maxDraws)
  if (!is.null(control$rhatTarget)) arguments$rhat_target <- as.numeric(control$rhatTarget)
  # `settleTol` ends warmup early once the metric stops moving between windows.
  # It is off by default and should stay off: measured on the N=200 augmented
  # marginal route it cost 1795 s for min ESS 142.8 where the fixed schedule
  # spent 559 s for min ESS 246.3, a factor of 5.5 against. A settled metric is
  # not a good metric, and the sampling phase pays for the shortened warmup on
  # every draw.
  if (!is.null(control$settleTol)) arguments$settle_tol <- as.numeric(control$settleTol)
  if (!is.null(hessian)) {
    arguments$hessian <- JuliaConnectoR::juliaPut(as.matrix(hessian))
  }
  joint <- identical(intoverpop, "none")
  if (joint) {
    arguments$npar <- as.integer(npar)
    arguments$save_effects <- saveEffects
    arguments$adapt_effects <- isTRUE(control$adapt_effects)
  }
  entry <- if (joint) module$ctsem_sample else module$ctsem_sample_marginal

  result <- .ctBackendWithMaxChunks(cores,
    JuliaConnectoR::juliaGet(do.call(entry, arguments)))

  # The shell the assembler fills, matching what an optimised fit carries so
  # that everything downstream reads a sampled fit the same way.
  subject_loglik <- as.numeric(optimised$subject_loglik)
  out <- list(backend = "julia", model = model, model_spec = model_spec,
    data = datalong,
    estimate = list(raw = estimate,
      loglik = if (length(subject_loglik)) sum(subject_loglik) else
        as.numeric(optimised$maximum_loglik),
      logposterior = as.numeric(optimised$maximum_loglik),
      converged = TRUE, chunks = as.integer(optimised$chunks)),
    engine = model_spec$engine,
    args = list(backend = "julia", backendcontrol = backendcontrol,
      optimcontrol = optimcontrol, cores = cores, priors = priors,
      intoverpop = intoverpop, optimize = FALSE))
  class(out) <- c("ctJuliaFit", "ctFitModel")
  out <- .ctBackendSampleAssemble(out, result, npar, saveEffects && joint,
    as.integer(chains), warmup, draws, hessian, estimate)
  out$identifiability <- .ctBackendIdentifiability(hessian,
    .ctBackendRawParameterNames(out, npar))
  out
}
