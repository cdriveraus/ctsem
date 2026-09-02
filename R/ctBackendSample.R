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
# effects with its members'. Reconstructing that here would mean reimplementing
# `_laplace_build_units` in R and keeping the two in step, so the engine is
# asked instead -- `ctsem_laplace_effect_layout` reports, for each position,
# which unit and level it belongs to and which subject the block starts at.
# That is the same information without the second implementation, and a
# plausible mislabelling would attach the wrong subject's name to a number and
# never announce itself.
#' @keywords internal
.ctBackendEffectIndex <- function(fit) {
  laplace <- fit$model_spec$laplace
  if (is.null(laplace) || is.null(laplace$levels)) return(NULL)
  if (length(laplace$levels) != 1L) return(.ctBackendEffectIndexNested(fit))
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

# The nested case, from the engine's own account of the layout.
#
# A block covering one member is that subject's; a block covering several is the
# group's, and the group is named by the grouping id its members share. The
# level's parameter names come from the fit, and `within` says which of them a
# position is -- so a label is the parameter, the level, and whose it is.
#' @keywords internal
.ctBackendEffectIndexNested <- function(fit) {
  laplace <- fit$model_spec$laplace
  layout <- try(JuliaConnectoR::juliaGet(
    .ctJuliaModule(fit$model_spec$project)$ctsem_laplace_effect_layout(
      .ctJuliaObjective(fit))), silent = TRUE)
  if (inherits(layout, "try-error") || is.null(layout$position)) return(NULL)
  level_index <- as.integer(layout$level)
  # A level the fit does not describe means the two have gone out of step, and
  # an unlabelled effect vector is better than a confidently wrong one.
  if (any(level_index < 1L) || any(level_index > length(laplace$levels))) return(NULL)
  within <- as.integer(layout$within)
  subject <- as.integer(layout$first_member)
  nmembers <- as.integer(layout$nmembers)

  parameter <- vapply(seq_along(within), function(i) {
    pars <- as.character(laplace$levels[[level_index[i]]]$param)
    if (within[i] >= 1L && within[i] <= length(pars)) pars[within[i]] else
      paste0("effect", within[i])
  }, character(1))
  levelname <- vapply(level_index, function(l) {
    nm <- laplace$levels[[l]]$name
    if (is.null(nm) || !nzchar(nm)) paste0("level", l) else as.character(nm)
  }, character(1))

  nsubjects <- length(fit$model_spec$subject_starts)
  ids <- .ctBackendSubjectIds(fit, nsubjects)
  who <- ifelse(subject >= 1L & subject <= nsubjects,
    if (is.null(ids)) as.character(subject) else ids[pmax(subject, 1L)],
    NA_character_)
  # For a grouping block the label should name the group, not one of its
  # members, so the group's own identifier is looked up from the data where the
  # fit kept it.
  group <- .ctBackendGroupIds(fit, laplace, levelname, subject, nmembers)
  label <- ifelse(nmembers > 1L,
    paste0(parameter, "_", levelname, "_", group),
    paste0(parameter, "_", who))
  data.frame(position = as.integer(layout$position), level = level_index,
    levelname = levelname, subject = ifelse(nmembers > 1L, NA_integer_, subject),
    id = ifelse(nmembers > 1L, NA_character_, who), group = group,
    parameter = parameter, label = label, stringsAsFactors = FALSE)
}

# The grouping identifier a block's members share, when the fit kept the column.
#' @keywords internal
.ctBackendGroupIds <- function(fit, laplace, levelname, subject, nmembers) {
  out <- rep(NA_character_, length(subject))
  d <- fit$data
  if (is.null(d)) return(out)
  # The subject identifier is not always called `id`: a nested model names its
  # levels, and the first level's name *is* the subject column. Looking only
  # for `id` is what left every group labelled NA.
  idname <- if (length(laplace$levels)) laplace$levels[[1L]]$name else NULL
  idcol <- if (!is.null(d$id)) d$id else
    if (!is.null(idname) && idname %in% names(d)) d[[idname]] else NULL
  if (is.null(idcol)) return(out)
  for (l in unique(levelname)) {
    if (!l %in% names(d)) next
    # First appearance order, which is how the engine numbers subjects.
    firstrow <- match(unique(idcol), idcol)
    bysubject <- as.character(d[[l]][firstrow])
    take <- levelname == l & nmembers > 1L & subject >= 1L &
      subject <= length(bysubject)
    out[take] <- bysubject[subject[take]]
  }
  out
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
#'   (1000). Stan's spellings \code{max_treedepth} and \code{adapt_delta},
#'   which \code{\link{ctFit}} takes for the same two settings, are
#'   accepted here as well.
#'
#'   Sampling takes exactly the draws it was asked for unless it is given a
#'   target to reach: \code{minEss} and \code{meanEss} are effective sample
#'   sizes to keep drawing towards, \code{rhatTarget} (1.01) the R-hat to reach
#'   alongside them, and \code{maxDraws} the budget that stops it. All are off
#'   by default, so nothing runs longer than asked without being told to, and
#'   setting \code{minEss} without \code{maxDraws} does nothing -- the budget
#'   is what the loop checks it against. Worth setting when a draw count had to
#'   be guessed at; not worth setting when a warning says a parameter is
#'   unidentified, because no number of draws fixes an improper posterior.
#'
#'   \code{adapt_effects} controls whether warmup re-estimates the
#'   random-effect blocks of the metric as well as the population block; they
#'   start from a conditional covariance that is exact for a linear model, so
#'   replacing one with an estimate from a few hundred draws can add more noise
#'   than it removes.
#' @param processes Run each chain in its own R process rather than its own
#'   thread, so that chains share no allocator and no garbage collector.
#'   \code{TRUE} by default whenever there is more than one chain.
#'
#'   A worker must start Julia and compile the engine for this model's
#'   dimensions before it can draw anything -- 26-43 seconds, unavoidable per
#'   process and not shareable between them. That is paid once against a
#'   sampling run that is normally minutes to hours, so it is worth it for any
#'   real run; set \code{FALSE} for very short ones, where it is the larger
#'   cost. Note also that each worker holds its own copy of the data and the
#'   adjoint tape, so memory scales with the number of chains.
#'
#'   Draws match the in-process path to about 1e-10 on the first draw and
#'   diverge chaotically from there, which is inherent rather than a defect --
#'   see \code{.ctBackendSampleProcesses} for why. Results will therefore not be
#'   bit-identical to a run made before this became the default. Needs the
#'   \pkg{future} package; without it, or if a worker fails, sampling falls back
#'   to this session.
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
  saveEffects = FALSE, seed = 20260828L, control = list(), verbose = FALSE,
  processes = TRUE) {

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

  # Chains in separate processes, which is the default when there is more than
  # one chain to separate. Dispatched here rather than deeper because the
  # process path does not share the engine call below at all: each worker runs
  # this same function with `chains = 1`, and the parent pools what comes back.
  #
  # On by default because the arithmetic is not close. A worker costs 26-43 s of
  # Julia startup and engine compilation, against a sampling run that is
  # normally minutes to hours -- the startup is noise at any realistic draw
  # count, and only dominates on the short runs used for testing. What it buys
  # is chains that contend for neither the allocator nor the garbage collector.
  #
  # A `NULL` back means it could not run and sampling continues here rather than
  # failing -- a slower answer beats none. The missing-package case is checked
  # separately because it is the only one that would otherwise be silent: a
  # failing worker warns on its way out, but an absent `future` just returns
  # nothing, and someone who asked for processes should be told why they did not
  # get them.
  if (isTRUE(processes) && chains > 1L) {
    # Silent when `future` is simply absent and the default put us here: that
    # is not the user's doing and there is nothing for them to act on. Said out
    # loud only when they asked for processes explicitly, via a call that named
    # the argument.
    if (!.ctBackendCanWarm()) {
      if ("processes" %in% names(match.call())) {
        message("processes = TRUE needs the future package, which is not ",
          "installed. Sampling in this session instead.")
      }
    } else {
      out <- .ctBackendSampleProcesses(fit, chains = chains, warmup = warmup,
        draws = draws, cores = cores, control = control,
        saveEffects = saveEffects, seed = seed, verbose = verbose)
      if (!is.null(out)) return(out)
      message("Sampling in this session instead.")
    }
  }

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
  .ctBackendSampleRun(fit, objective, estimate, npar, chains = chains,
    warmup = warmup, draws = draws, cores = cores, saveEffects = saveEffects,
    seed = seed, control = control, verbose = verbose, hessian = hessian)
}

# The sampler settings, from either spelling of the control list.
#
# `ctFit(optimize = FALSE)` took Stan's names for two of these and `ctSample()`
# takes the engine's, so both are read here rather than each entry point
# quietly ignoring what the other documents. The rest are spelled the same on
# both, and the whole list is assembled in one place so that a knob added for
# one cannot go missing from the other -- which is how `minEss` and its three
# companions came to be documented on `ctSample()` and passed only by `ctFit()`.
#' @keywords internal
.ctBackendSampleControl <- function(control) {
  control <- .ctJuliaOr(control, list())
  settings <- list(
    maxdepth = as.integer(.ctJuliaOr(control$maxdepth,
      .ctJuliaOr(control$max_treedepth, 10L))),
    target_accept = as.numeric(.ctJuliaOr(control$target_accept,
      .ctJuliaOr(control$adapt_delta, 0.8))),
    maxdelta = as.numeric(.ctJuliaOr(control$maxdelta, 1000)),
    # 2, not 1. Chains are dispersed by drawing from the Laplace approximation,
    # which has the right shape and the wrong width: Laplace understates spread
    # wherever the posterior is skewed or heavy-tailed, which is the case
    # `ctLaplaceCorrect` exists to repair. Starting every chain from that
    # narrower distribution makes R-hat compare chains that began already
    # agreeing, so it reads low exactly where the approximation is worst and the
    # diagnostic is needed most. Over-dispersing costs a little warmup and
    # would restore the comparison Stan gets from its uniform [-2, 2] start.
    #
    # **Left at 1 until that is measured rather than argued.** The reasoning
    # above is the same shape as the reasoning that made `settle_tol` look
    # obviously right, and `settle_tol` lost by a factor of 5.5. Raising this
    # moves every user's diagnostics, so it needs a comparison on an identified
    # model -- the one attempted here diverged on every transition at 1, 1.5 and
    # 2 alike, so it discriminated nothing. `control$init_scale` takes any value
    # meanwhile.
    init_scale = as.numeric(.ctJuliaOr(control$init_scale, 1)),
    adapt_metric = isTRUE(.ctJuliaOr(control$adapt_metric, TRUE)),
    adapt_effects = isTRUE(.ctJuliaOr(control$adapt_effects, FALSE)))

  # Sampling targets, when asked for, and absent from the call when not: the
  # engine reads zero as "no target", so an unset element here and an omitted
  # argument there mean the same thing. Left unset the sampler takes exactly
  # the draws it was told to; set, it keeps going until the effective sample
  # size is there or the budget runs out, which is usually what a user wanted
  # from a draw count they had to guess at.
  #
  # `settleTol` ends warmup early once the metric stops moving between windows.
  # It is off by default and should stay off: measured on the N=200 augmented
  # marginal route it cost 1795 s for min ESS 142.8 where the fixed schedule
  # spent 559 s for min ESS 246.3, a factor of 5.5 against. A settled metric is
  # not a good metric, and the sampling phase pays for the shortened warmup on
  # every draw.
  if (!is.null(control$minEss)) settings$min_ess <- as.numeric(control$minEss)
  if (!is.null(control$meanEss)) settings$mean_ess <- as.numeric(control$meanEss)
  if (!is.null(control$maxDraws)) settings$max_draws <- as.integer(control$maxDraws)
  if (!is.null(control$rhatTarget)) settings$rhat_target <- as.numeric(control$rhatTarget)
  if (!is.null(control$settleTol)) settings$settle_tol <- as.numeric(control$settleTol)
  settings
}

# Call the engine's sampler and assemble the fit it produced.
#
# The one place either entry point reaches the sampler from. `ctSample()` comes
# here with a fit it was handed and `ctFit(optimize = FALSE)` with a shell it
# has just optimised, and that -- which object gets filled, and what placed the
# sampler -- is the whole difference between them. Everything after it was
# duplicated until it drifted: the two argument lists had diverged over five
# settings, and the ones only `ctFit()` passed were documented on `ctSample()`
# as though they worked.
#
# `marginal` selects `ctsem_sample_marginal` over `ctsem_sample`: the
# parameters alone, with whatever integration the objective already does,
# against the joint posterior over parameters and random effects. It is the
# only structural difference in the call, because `npar`, `save_effects` and
# `adapt_effects` all describe an effect block the marginal entry does not
# have.
#
# `estimate` may be longer than `npar` -- the state-explicit route samples the
# trajectory alongside the parameters -- so the parameter block is sliced out
# for the assembler rather than assumed to be the whole vector.
#
# `progress` is separate from `verbose` because the two paths decide it
# differently: a flag on `ctSample()`, and on the fitting path anyone watching
# a console, since sampling there follows an optimisation that has already been
# printing and silence after it reads as a finished run rather than a running
# one.
#' @keywords internal
.ctBackendSampleRun <- function(fit, objective, estimate, npar, chains, warmup,
  draws, cores, saveEffects, seed, control, verbose, hessian, marginal = FALSE,
  progress = isTRUE(verbose)) {

  settings <- .ctBackendSampleControl(control)
  module <- .ctJuliaModule(fit$model_spec$project)
  arguments <- list(objective, .ctJuliaNumericVector(estimate),
    nchains = as.integer(chains), nwarmup = as.integer(warmup),
    ndraws = as.integer(draws), seed = as.integer(seed)[1L],
    maxdepth = settings$maxdepth, target_accept = settings$target_accept,
    maxdelta = settings$maxdelta, init_scale = settings$init_scale,
    adapt_metric = settings$adapt_metric,
    verbose = isTRUE(progress),
    progress_overwrite = .ctProgressOverwrite(verbose))
  for (name in c("min_ess", "mean_ess", "max_draws", "rhat_target", "settle_tol")) {
    if (!is.null(settings[[name]])) arguments[[name]] <- settings[[name]]
  }
  if (!marginal) {
    arguments$npar <- as.integer(npar)
    arguments$save_effects <- isTRUE(saveEffects)
    arguments$adapt_effects <- settings$adapt_effects
  }
  if (!is.null(hessian)) {
    arguments$hessian <- JuliaConnectoR::juliaPut(as.matrix(hessian))
  }
  entry <- if (marginal) module$ctsem_sample_marginal else module$ctsem_sample

  result <- .ctBackendWithMaxChunks(cores,
    JuliaConnectoR::juliaGet(do.call(entry, arguments)))
  # `result$ndraws` rather than the count asked for: with an effective sample
  # size target the sampler decides when to stop, and reporting the request
  # would describe a run that did not happen.
  .ctBackendSampleAssemble(fit, result, npar, isTRUE(saveEffects) && !marginal,
    as.integer(chains), warmup, as.integer(result$ndraws), hessian,
    as.numeric(estimate)[seq_len(npar)])
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

  # Which sampled coordinates reached the flat region of their transform.
  # Computed here because the warner sees only the diagnostics, and the draws
  # and the parameter names both live at this level.
  out$sample$unidentified <- .ctBackendFlatParameters(out, npar)
  .ctSampleWarn(out$sample)
  out
}

# Past |raw| ~ 20 every ctsem transform is flat to machine precision. A sampled
# coordinate that gets there is not mixing badly, it has no posterior to mix
# over: the likelihood cannot separate those values and, without a prior, the
# density is improper in that direction.
#' @keywords internal
.ctBackendFlatParameters <- function(fit, npar) {
  draws <- fit$estimate$rawposterior
  if (is.null(draws) || !length(draws)) return(character(0))
  reach <- suppressWarnings(apply(abs(as.matrix(draws)), 2, max, na.rm = TRUE))
  flat <- which(is.finite(reach) & reach >= 20)
  flat <- flat[flat <= npar]
  if (!length(flat)) return(character(0))
  labels <- .ctBackendRawParameterNames(fit, npar)
  if (length(labels) < npar) labels <- paste0("par", seq_len(npar))
  labels[flat]
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
  # Checked before R-hat, because it changes what a bad R-hat means.
  #
  # Past |raw| ~ 20 every ctsem transform is flat to machine precision, so the
  # likelihood cannot distinguish one value from another and, with no prior to
  # hold it, the posterior is improper in that direction. A chain does the only
  # thing it can: it random-walks away. Measured on a model sampled with
  # priors=FALSE, two of six parameters had per-chain means of -1439, -672865
  # and -546340 while the other four returned R-hat 1.000 at an effective size
  # of 3000 -- so the fit was perfectly good apart from directions that had no
  # posterior at all.
  #
  # Reported separately because R-hat sends the reader to the wrong problem.
  # Running longer cannot fix an improper posterior, and neither can a smaller
  # step size; a prior can, and so can removing the parameter.
  flat <- if (is.null(diagnostics$unidentified)) character(0) else
    diagnostics$unidentified
  if (length(flat)) {
    warning(length(flat), " sampled parameter(s) reached the region where ",
      "their transform is flat to machine precision: ",
      paste(utils::head(flat, 5), collapse = ", "),
      if (length(flat) > 5) ", ..." else "",
      ". The likelihood cannot tell those values apart, so with no prior the ",
      "posterior is improper there and the chains wander rather than mix. ",
      "More draws will not help. Set priors=TRUE, or fix or remove the ",
      "parameter. See fit$estimate$rawposterior.", call. = FALSE)
  }
  worst <- suppressWarnings(max(diagnostics$rhat, na.rm = TRUE))
  if (is.finite(worst) && worst > 1.01) {
    # Naming the arguments, because "run longer" is advice the reader then has
    # to go and look up. Both entry points land here and they are controlled
    # differently: `ctFit` takes `iter`, which counts warmup and sampling
    # together, and `ctSample` takes `draws` directly.
    warning("Largest R-hat is ", signif(worst, 4), ". The chains have not ",
      "agreed on the same distribution, so the draws are not yet a posterior. ",
      if (length(flat))
        "That is expected for the unidentified parameter(s) named above, which more draws cannot fix. For the rest, "
      else "",
      if (length(flat)) "raise" else "Raise", " the draw count -- iter in ctFit (now ",
      diagnostics$warmup + diagnostics$draws, ", of which ",
      diagnostics$warmup, " is warmup, leaving ", diagnostics$draws,
      " per chain) or draws in ctSample -- or set control$minEss together with ",
      "control$maxDraws to keep sampling until an effective size is reached. ",
      "See fit$sample$rhat.", call. = FALSE)
  }
  fewest <- suppressWarnings(min(diagnostics$ess, na.rm = TRUE))
  if (is.finite(fewest) && fewest < 100) {
    warning("Smallest effective sample size is ", round(fewest), ", from ",
      total, " draws. Interval estimates from this few are unreliable. Raise ",
      "the draw count (iter in ctFit, draws in ctSample), or set ",
      "control$minEss with control$maxDraws to keep sampling until an ",
      "effective size is reached. See fit$sample$ess.", call. = FALSE)
  }
  if (diagnostics$saturated > 0L) {
    # Raising the cap is the mechanical answer and rarely the right first one.
    # Saturation means the sampler wanted longer trajectories than it was
    # allowed, and it wants them because the geometry is hard -- a badly scaled
    # metric, or a funnel. Raising `maxdepth` buys those trajectories at double
    # the cost per draw without touching the cause, so the cause is named first.
    message(diagnostics$saturated, " of ", total, " transitions hit the maximum ",
      "tree depth of ", diagnostics$max_depth, ", so the sampler was cut off ",
      "before it finished exploring. That costs efficiency rather than ",
      "correctness. It usually reflects difficult posterior geometry rather ",
      "than a cap set too low",
      if (diagnostics$divergent > 0L)
        " -- the divergences above point the same way" else "",
      "; check the divergence count and any near-zero population standard ",
      "deviation before raising control$maxdepth.")
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
  gradient, verbose, intoverstates = TRUE) {

  npar <- .ctBackendNpar(model_spec)
  # As in the optimising path: the zero keeps `max` from warning and returning
  # -Inf on a fully fixed model, and the refusal replaces the "invalid
  # arguments" that -Inf produced two lines later. A sampler with no
  # dimensions to move in is worse than a sentence saying so.
  if (npar < 1L) {
    stop("This model has no free parameters, so there is nothing to sample. ",
      "Free a parameter, or use fit = FALSE to prepare the model without ",
      "fitting it.", call. = FALSE)
  }
  start <- .ctJuliaInitialValues(npar, inits)

  # Stan's vocabulary, because these are Stan's arguments: `iter` counts warmup
  # and sampling together and warmup is half of it unless said otherwise.
  warmup <- as.integer(.ctJuliaOr(control$warmup, max(1L, floor(iter / 2))))
  draws <- max(1L, as.integer(iter) - warmup)
  seed <- as.integer(.ctJuliaOr(control$seed, 20260828L))
  saveEffects <- isTRUE(optimcontrol$saveEffects)

  # Whenever the progress line below it will be drawn, not only at `verbose`.
  #
  # Sampling begins with an optimisation, to place the sampler and build its
  # metric, and that optimisation prints a progress line labelled "optimise".
  # Without this sentence in front of it a user who asked for HMC watches an
  # optimiser run and concludes the sampler never started -- which is exactly
  # what was reported. The explanation costs one line and was previously hidden
  # behind `verbose > 0`, which is not the default.
  announce <- verbose > 0L || .ctProgressConsole()
  if (announce) {
    message("Sampling: optimising first, to place the sampler and build its ",
      "metric. Sampling follows.")
  }
  # The state-explicit target, when asked for. This is the estimator the path
  # is really for: NUTS over the parameters *and* the states is the exact
  # posterior of both, with no Gaussian assumption about the state anywhere,
  # where optimising the same density gives its joint mode and the downward
  # bias in the variances that comes with maximising over what should be
  # integrated.
  jointobjective <- NULL
  nstate <- 0L
  if (!isTRUE(intoverstates)) {
    jointobjective <- .ctJuliaJointObjective(model_spec, npar)
    nstate <- .ctJuliaStateDimension(model_spec)
    start <- c(start, numeric(nstate))
  }
  optimised <- .ctJuliaOptimise(model_spec, start, backendcontrol = backendcontrol,
    gradient = gradient, cores = cores, verbose = verbose,
    objective = jointobjective)
  estimate <- as.numeric(optimised$minimizer)

  spec <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(model_spec$project)
  objective <- .ctJuliaObjective(spec)
  # `spec`, not a bare list carrying `model_spec`: .ctBackendHessian() reaches
  # the objective through .ctJuliaObjective(), which requires a classed
  # ctJuliaModel/ctJuliaFit and errors on anything else. An unclassed list made
  # that error every time, and .ctBackendHessian() catches its own errors and
  # returns NULL -- so every sampled fit warned that the engine could not
  # differentiate its gradient, when nothing had been asked of the engine at
  # all. `spec` is the same object the objective is already cached under.
  # The metric's curvature. On the state-explicit target that is the whole
  # arrow-shaped joint Hessian rather than the profiled one the standard
  # errors use: the sampler moves in every coordinate, so it needs the
  # curvature of every coordinate.
  hessian <- if (is.null(jointobjective)) {
    try(.ctBackendHessian(spec, estimate, verbose = verbose), silent = TRUE)
  } else {
    try(matrix(as.numeric(.ctBackendJuliaValue(module$ctsem_joint_hessian(
      jointobjective, .ctJuliaNumericVector(estimate), profile = FALSE))),
      nrow = length(estimate), ncol = length(estimate)), silent = TRUE)
  }
  if (inherits(hessian, "try-error")) hessian <- NULL

  joint <- identical(intoverpop, "none")

  # The shell the assembler fills, matching what an optimised fit carries so
  # that everything downstream reads a sampled fit the same way.
  subject_loglik <- as.numeric(optimised$subject_loglik)
  # The population block alone, as everywhere else: `estimate$raw` means the
  # parameters on every fit, and the assembler overwrites it with the
  # posterior mean of exactly those.
  theta <- estimate[seq_len(npar)]
  out <- list(backend = "julia", model = model, model_spec = model_spec,
    data = datalong,
    estimate = list(raw = theta,
      loglik = if (length(subject_loglik)) sum(subject_loglik) else
        as.numeric(optimised$maximum_loglik),
      logposterior = as.numeric(optimised$maximum_loglik),
      converged = TRUE, chunks = as.integer(optimised$chunks)),
    engine = model_spec$engine,
    args = list(backend = "julia", backendcontrol = backendcontrol,
      optimcontrol = optimcontrol, cores = cores, priors = priors,
      intoverpop = intoverpop, optimize = FALSE,
      intoverstates = isTRUE(intoverstates)))
  class(out) <- c("ctJuliaFit", "ctFitModel")
  out <- .ctBackendSampleRun(out, .ctJuliaOr(jointobjective, objective),
    estimate, npar, chains = chains, warmup = warmup, draws = draws,
    cores = cores, saveEffects = saveEffects, seed = seed, control = control,
    verbose = verbose, hessian = hessian, marginal = !joint,
    # On when someone is watching, matching the optimiser rather than
    # differing from it. The two run one after the other in this same call,
    # and having the first print progress by default while the second stayed
    # silent is what made a running sampler look like a finished
    # optimisation: the visible output stopped at "Computing exact Hessian"
    # and nothing followed it for several minutes.
    progress = verbose > 0L || .ctProgressConsole())
  # The identifiability report is about the parameters, so it is given the
  # parameter block's own curvature -- the profiled one on the state route,
  # not the corner of the joint matrix, which describes the parameters at a
  # trajectory held fixed.
  identhessian <- hessian
  if (!is.null(jointobjective)) {
    identhessian <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_joint_hessian(jointobjective,
        .ctJuliaNumericVector(estimate), profile = TRUE))),
      nrow = npar, ncol = npar), silent = TRUE)
    if (inherits(identhessian, "try-error")) identhessian <- NULL
    out$estimate$innovations <- estimate[npar + seq_len(nstate)]
    out$estimate$states <- try(.ctJuliaJointStates(model_spec,
      jointobjective, estimate), silent = TRUE)
    if (inherits(out$estimate$states, "try-error")) out$estimate$states <- NULL
    out$estimate$loglik_type <- "joint"
  }
  out$identifiability <- .ctBackendIdentifiability(identhessian,
    .ctBackendRawParameterNames(out, npar))
  out
}
