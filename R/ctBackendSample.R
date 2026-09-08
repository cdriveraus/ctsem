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
  # `.ctBackendSpec(fit)$data`, not the top-level `fit$data`: the latter is now
  # the sentinel-cleaned `standata` structure (for `$data`/`$standata` parity
  # with a stan fit, see R/ctFit.R) and no longer a data.frame with the user's
  # original column names. `model_spec$data` is the long data.frame
  # `.ctJuliaPrepare()` kept for exactly this kind of lookup and is unaffected
  # by that change.
  d <- .ctBackendSpec(fit)$data
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
  # `.ctBackendSpec(fit)$data`, not the top-level `fit$data` -- see the same
  # note in `.ctBackendGroupIds()` above.
  d <- .ctBackendSpec(fit)$data
  if (!is.null(d) && !is.null(d$id)) {
    original <- unique(d$id)
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
#' normal-approximation draws. Those draws are the thing this is not: an
#' optimised fit's \code{rawposterior} holds pseudo-posterior draws from a
#' covariance fitted around the mode, written by
#' \code{\link{ctOptimUncertainty}}, and no number of them makes a posterior
#' sample. This writes genuine posterior draws into the same slot and marks the
#' fit with \code{$sample}, which is how the two are told apart afterwards.
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
#' Because a warning is only seen by whoever is at the console, and any batch
#' script wraps its fitting call in \code{suppressWarnings()}, the verdict is
#' also kept on the object: \code{fit$sample$converged} says whether the chains
#' agreed on one distribution and \code{fit$sample$diagnosis} lists what went
#' wrong, both of which \code{print()} shows. \code{summary()} reports
#' \code{n_eff} and \code{Rhat} beside every estimate, exactly as it does for
#' \code{backend='stan'}, and opens with a line naming the worst of each.
#'
#' The two R-hats a fit carries are not the same statistic and are not meant to
#' be. \code{fit$sample$rhat} is the engine's split R-hat over the raw,
#' unconstrained coordinates -- the one the sampler's own \code{rhatTarget}
#' stopping rule reads, and unbounded, so a badly failed run shows a number in
#' the thousands. The \code{Rhat} column in \code{summary()} is
#' \code{rstan::monitor}'s rank-normalised split R-hat over the transformed
#' quantities in the table, which is what \code{backend='stan'} reports and is
#' deliberately robust rather than dramatic. Both cross 1.01 on the same runs.
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
#' @param control \strong{Deprecated} -- use \code{sampleControl}. Still
#'   honoured, with a warning.
#' @param sampleControl A list of sampler settings: \code{maxdepth} (default 10),
#'   \code{target_accept} (0.8), \code{adapt_metric} (FALSE),
#'   \code{adapt_effects} (FALSE), \code{init_scale} (1), \code{stepsize},
#'   \code{maxdelta}
#'   (1000). Stan's spellings \code{max_treedepth} and \code{adapt_delta},
#'   which \code{\link{ctFit}} takes for the same two settings, are
#'   accepted here as well.
#'
#'   Sampling takes exactly the draws it was asked for unless it is given a
#'   target to reach: \code{minESS} and \code{meanESS} are effective sample
#'   sizes to draw towards and \code{rhatTarget} (1.01) the R-hat to reach
#'   alongside them. Given one, the count asked for becomes a \emph{budget}
#'   rather than an instruction -- the run stops as soon as the target is met,
#'   and never draws more than was asked unless \code{maxDraws} says so
#'   explicitly. It is checked in batches: the first is sized from the target
#'   rather than from the budget -- effective size cannot exceed the draws
#'   behind it, so \code{minESS} needs at least \code{minESS / chains} of them,
#'   with a floor of 50 because R-hat and effective size read off fewer are too
#'   noisy to stop on.
#'
#'   A small target does not buy a short run, and this is the half that
#'   surprises: the rule is min ESS \emph{and} mean ESS \emph{and}
#'   \code{rhatTarget}, so at a small size target R-hat is what binds. Asked
#'   for \code{minESS = 100} on five chains, a well behaved model met the size
#'   target at the first check of 50 draws (250 effective) and went on to 86
#'   because R-hat was still 1.027 there; it stopped with 430. Raise
#'   \code{rhatTarget} to let the size target decide alone.
#'
#'   That is a change: \code{minEss} used to do nothing at all unless
#'   \code{maxDraws} was also set, because the budget defaulted to exactly the
#'   draws asked for and the loop had nothing to extend into. Setting an
#'   effective size and watching the sampler run to the end regardless is what
#'   this fixes.
#'
#'   Worth setting when a draw count had to be guessed at; not worth setting
#'   when a warning says a parameter is unidentified, because no number of
#'   draws fixes an improper posterior.
#'
#'   \code{stepsize} fixes the step size every chain starts from, instead of
#'   each chain estimating its own from a single trial leapfrog step -- which
#'   answers differently in every chain, and is the whole of what a chain keeps
#'   when \code{warmup} is 0. Dual averaging moves it from there unless warmup
#'   is 0, so this is a starting point rather than a setting of the step size
#'   itself.
#'
#'   \code{adapt_metric} re-estimates the metric during warmup. It is off by
#'   default: the metric starts as the inverse of the exact Hessian at the
#'   mode, and a sample covariance from a few hundred warmup draws is measured
#'   to be worse -- a smaller step size, deeper trees, 15-40% more time for the
#'   same effective sample, on every model tried. Worth turning on where the
#'   starting curvature had to be repaired or floored, which is where the
#'   estimate has something to improve on.
#'
#'   \code{adapt_effects} controls whether warmup re-estimates the
#'   random-effect blocks of the metric as well as the population block; they
#'   start from a conditional covariance that is exact for a linear model, so
#'   replacing one with an estimate from a few hundred draws can add more noise
#'   than it removes.
#'
#'   \code{control$callback} is a function called while sampling runs, with
#'   \code{(phase, iteration, total, logp, divergent)}: \code{phase} is
#'   \code{"warmup"} or \code{"sampling"}, \code{iteration}/\code{total} count
#'   against the current phase, and \code{logp}/\code{divergent} are the log
#'   posterior and divergence count so far. It is for a front end that wants
#'   to draw progress live; the engine calls it on a time cadence rather than
#'   once per iteration (see \code{optimcontrol$callback} in
#'   \code{\link{ctFit}}), and always once more when a phase ends. An error
#'   inside it disables it and warns, leaving the sample unaffected. With
#'   several chains only the first calls back, matching the printed line --
#'   several chains calling into R at once is not just unreadable, it is
#'   unsafe. Under \code{processes = TRUE} it is not called at all, because a
#'   worker process cannot call back into this session's callback; the
#'   parent's own per-chain lines (see \code{verbose}) are what cover that
#'   case instead.
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
#'   diverge chaotically from there. That is inherent rather than a defect:
#'   each unit's mode comes from an inner Newton solve warm-started from
#'   whatever the objective last held, and this session carries an
#'   optimisation's worth of that history where a fresh worker carries one
#'   warm-up evaluation. Both reach the same mode to solver tolerance rather
#'   than to the last bit, and NUTS is chaotic, so 1e-10 becomes order 1
#'   within a few dozen transitions. Results will therefore not be
#'   bit-identical to a run made before this became the default. Needs the
#'   \pkg{future} package; without it, or if a worker fails, sampling falls back
#'   to this session.
#' @param verbose Report progress while sampling runs: warmup and sampling
#'   separately, iterations against the total, and an estimated time
#'   remaining, the same shape \code{\link{ctFit}}'s progress line has. A
#'   logical flag is accepted as well as a level, as elsewhere in ctsem --
#'   \code{FALSE}/\code{0} silent (the default), \code{TRUE}/\code{1} the
#'   progress just described, \code{2} the same reporting kept as scrolling
#'   history rather than overwritten in place, which is what \code{verbose =
#'   2} means throughout the julia backend and there is nothing further to add
#'   for the sampler specifically. With \code{chains > 1} and \code{processes
#'   = TRUE} (the default above one chain), each worker's own printed line
#'   never reaches this session, so the parent prints one line per chain
#'   instead, polling what the workers have done so far -- the single-process
#'   line's content, relayed rather than duplicated.
#'
#'   Progress overwrites a single line where the output is going to a console
#'   and prints occasional separate lines where it is not; set
#'   \code{options(ctsem.progress.overwrite = FALSE)} if that detection is
#'   wrong for your front end, or \code{TRUE} to force it on. The per-chain
#'   lines under \code{processes = TRUE} are never overwritten in place --
#'   several chains share the console, and one finishing should not erase an
#'   earlier line that is still current for another.
#'
#' @return The fit, with \code{estimate$rawposterior} holding the draws and
#'   \code{$sample} holding the diagnostics: split R-hat and effective sample
#'   size per parameter, divergences, tree depths, step sizes and E-BFMI, plus
#'   \code{converged} and \code{diagnosis} summarising them.
#'   \code{estimate$raw} is set to the per-parameter posterior mean of the
#'   draws -- unlike \code{backend='stan'}'s sampled point estimate
#'   (\code{ctFit(..., optimize=FALSE)}'s \code{stanfit$rawest}), which is the
#'   per-parameter median.
#'
#' @seealso \code{\link{ctLaplaceCheck}} measures the Laplace approximation's
#'   error and corrects it to first order, at a small fraction of the cost;
#'   \code{\link{ctOptimUncertainty}} for pseudo-posterior draws from a fitted
#'   covariance, which is the cheap approximation this replaces rather than
#'   extends; \code{\link{ctFit}} for the fit this starts from, and
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
  saveEffects = FALSE, seed = 20260828L, sampleControl = list(),
  verbose = FALSE, processes = TRUE, control = list()) {
  # `control` renamed to `sampleControl`, as on `ctFit()`, where the same list
  # had to be told apart from rstan's `control`. Still accepted, at the end of
  # the signature so that nobody's positional call quietly means something new,
  # and under `sampleControl` where both name a setting.
  if ("control" %in% names(match.call())) {
    warning("ctSample(control = ) is deprecated: it is sampleControl now. ",
      "What was passed still takes effect.", call. = FALSE)
    for (name in setdiff(names(control), names(sampleControl))) {
      sampleControl[[name]] <- control[[name]]
    }
  }

  # The five settings that are also arguments here can be written in either
  # place -- `ctFit()` has only the list, so a script moving between the two
  # should not have to move them. Both at once is refused rather than resolved
  # by a precedence rule nobody would remember.
  supplied <- names(match.call())
  for (name in c("chains", "warmup", "draws", "seed", "saveEffects",
      "processes")) {
    if (is.null(sampleControl[[name]])) next
    if (name %in% supplied) {
      stop("ctSample(", name, " = ) and sampleControl$", name,
        " were both given. Use one.", call. = FALSE)
    }
    assign(name, sampleControl[[name]])
    sampleControl[[name]] <- NULL
  }
  # Wrapped as `ctFit(backend='julia')` is wrapped, and it was not: an
  # interrupted `ctSample()` left the Julia session desynchronised and the
  # worker pool full of orphaned chains, with nothing to put either right. This
  # is the entry point most likely to be interrupted -- it is the one that runs
  # for an hour.
  # Whether the caller *named* `processes` is decided here and passed on: it is
  # a fact about this call, and `match.call()` one frame down would report the
  # arguments this line writes rather than the ones the user wrote -- so the
  # message below would have fired for everyone rather than for the caller who
  # asked for something they are not getting.
  .ctJuliaInterruptSafe(.ctSampleImpl(fit, chains = chains, warmup = warmup,
    draws = draws, cores = cores, saveEffects = saveEffects, seed = seed,
    control = sampleControl, verbose = verbose, processes = processes,
    processes_named = "processes" %in% supplied))
}

#' @keywords internal
.ctSampleImpl <- function(fit, chains, warmup, draws, cores, saveEffects, seed,
  control, verbose, processes, processes_named = FALSE) {

  if (!inherits(fit, "ctJuliaFit")) {
    stop("ctSample applies to fits made with ctFit(backend='julia').", call. = FALSE)
  }
  .ctBackendSampleCheckControl(control)
  if (is.null(fit$model_spec$laplace)) {
    stop("ctSample needs a fit made with intoverpop='laplace'. The augmented ",
      "route carries the random effects in the state, so there is no separate ",
      "posterior over them to sample.", call. = FALSE)
  }
  chains <- max(1L, as.integer(chains)[1L])
  warmup <- max(0L, as.integer(warmup)[1L])
  draws <- max(1L, as.integer(draws)[1L])
  cores <- max(1L, as.integer(cores)[1L])

  # Whether the chains get processes is decided in `.ctBackendSampleRun()`,
  # which both entry points share. The one part that cannot move is this
  # message: it depends on whether *this* call named the argument. Silent when
  # `future` is simply absent and the default put us here -- that is not the
  # user's doing and there is nothing for them to act on -- and said out loud
  # when they asked for processes and are not getting them.
  if (isTRUE(processes) && chains > 1L && !.ctBackendCanWarm() &&
      isTRUE(processes_named)) {
    message("processes = TRUE needs the future package, which is not ",
      "installed. Sampling in this session instead.")
  }

  # The joint posterior over parameters and effects, started from the fit's own
  # estimate and metered by its curvature -- the sampler would otherwise
  # recompute that Hessian, at 2n gradient evaluations it need not spend.
  # `ctSample()` samples a Laplace fit and this is the target it is the exact
  # counterpart of; `ctFit(optimize = FALSE)` reaches the same runner with
  # whichever target its `intoverpop` and `intoverstates` chose.
  target <- .ctBackendSampleTarget(estimate = as.numeric(fit$estimate$raw),
    npar = length(fit$estimate$raw), hessian = fit$uncertainty$hessian)
  .ctBackendSampleRun(fit, target, chains = chains, warmup = warmup,
    draws = draws, cores = cores, saveEffects = saveEffects, seed = seed,
    control = control, verbose = verbose, processes = processes)
}

# Every name the sampler reads out of `control`, in one place.
#
# Kept beside `.ctBackendSampleControl()` because that is what reads most of
# them: a knob added there and not added here is refused as a typo, which is a
# loud failure and the right way round. The four it does not read are read by
# `.ctJuliaSampleFit()` (`warmup`, `seed`, `processes`) and
# `.ctBackendSampleEngine()` (`callback`).
# Whether the deprecation has already been said this session.
.ct_sample_deprecation <- new.env(parent = emptyenv())

.CT_SAMPLE_CONTROL_NAMES <- c(
  "maxdepth", "max_treedepth", "target_accept", "adapt_delta", "maxdelta",
  "init_scale", "adapt_metric", "adapt_effects",
  "minESS", "meanESS", "maxDraws", "rhatTarget", "settleTol",
  # How much to draw and how, which `ctFit()` used to take as arguments of its
  # own. `iter` is kept because it is what `chains` and `warmup` were always
  # expressed against, and because a script that passed it should keep working
  # through `sampleControl` as well as through the deprecated argument.
  "iter", "chains", "warmup", "draws", "seed", "saveEffects", "processes",
  "stepsize",
  "callback")

#' Fold the deprecated sampling arguments into \code{sampleControl}
#'
#' \code{ctFit()} took \code{iter}, \code{chains} and \code{control} as
#' arguments of its own, which put three of a sampler's settings in the
#' signature and the rest in a list -- and made \code{control} mean rstan's
#' control list on one backend and the julia sampler's settings on the other.
#' They are all entries of \code{sampleControl} now.
#'
#' Resolved here rather than threaded through: everything downstream, on both
#' backends, still receives \code{iter}, \code{chains} and \code{control} as
#' it always did, so the rename cannot change what a fit does. In particular
#' \code{control} still reaches \code{\link[rstan]{stan}} unchanged on the
#' stan path.
#'
#' Precedence is \code{sampleControl} first, because it is the argument being
#' kept. A deprecated argument that was *explicitly supplied* is used only for
#' what \code{sampleControl} does not say, and is reported either way -- a
#' silent deprecation teaches nobody, and a silently ignored one is worse.
#'
#' @param sampleControl The new list.
#' @param given Names the caller actually used, from \code{names(match.call())}
#'   in the calling function -- a default cannot be told from a value that
#'   happens to equal it any other way.
#' @param iter,chains,control The deprecated arguments, as they arrived.
#' @return A list of \code{iter}, \code{chains} and \code{control} to carry on
#'   with.
#' @keywords internal
.ctSampleControlResolve <- function(sampleControl = list(), given = character(0),
  iter = 1000L, chains = 2L, control = list()) {
  sampleControl <- .ctJuliaOr(sampleControl, list())
  control <- .ctJuliaOr(control, list())
  if (!is.list(sampleControl)) {
    stop("sampleControl must be a list of named sampler settings.", call. = FALSE)
  }
  deprecated <- intersect(c("iter", "chains", "control"), given)
  if (length(deprecated) && !isTRUE(.ct_sample_deprecation$said)) {
    .ct_sample_deprecation$said <- TRUE
    warning("ctFit(", paste(paste0(deprecated, " = "), collapse = ", "),
      ") is deprecated: sampling settings are entries of sampleControl now, ",
      "as in sampleControl = list(iter = 2000, chains = 4, warmup = 500). ",
      "What was passed here still takes effect. Said once per session.",
      call. = FALSE)
  }

  # The deprecated `control` sits *under* `sampleControl`: an entry in both is
  # the new argument's.
  merged <- sampleControl
  for (name in setdiff(names(control), names(sampleControl))) {
    merged[[name]] <- control[[name]]
  }

  # `iter` and `chains` move into the list, so the list is where they are read
  # from; the arguments fill in only when the list is silent about them.
  resolved_iter <- if (!is.null(merged$iter)) merged$iter else iter
  resolved_chains <- if (!is.null(merged$chains)) merged$chains else chains
  merged$iter <- NULL
  merged$chains <- NULL

  list(iter = as.integer(resolved_iter)[1L],
    chains = as.integer(resolved_chains)[1L], control = merged)
}

# Refuse a control name nothing reads.
#
# `control$minESS` on `list(minEss = 100)` is NULL -- `$` on a list matches
# exactly or by unique prefix, and neither folds case. So the
# setting was accepted, ignored, and the run reported nothing about it: reported
# from a real session as a sampler that would not stop at an effective size the
# user had asked for. That is the shape this package refuses by name elsewhere
# (see the note on arguments honoured by one backend and ignored by the other),
# and there is no reason for the control list to be the exception.
#
# An error rather than a warning, and before any fitting rather than at the
# sampler: the whole cost of getting this wrong is a long run that answers a
# different question, and the cure is one character. A near match is named
# because case is what goes wrong most.
#' @keywords internal
.ctBackendSampleCheckControl <- function(control) {
  control <- .ctJuliaOr(control, list())
  if (!length(control)) return(invisible(NULL))
  if (is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("Every entry of control must be named. Accepted names: ",
      paste(sort(.CT_SAMPLE_CONTROL_NAMES), collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(names(control), .CT_SAMPLE_CONTROL_NAMES)
  if (!length(unknown)) return(invisible(NULL))
  hint <- vapply(unknown, function(name) {
    near <- .CT_SAMPLE_CONTROL_NAMES[tolower(.CT_SAMPLE_CONTROL_NAMES) ==
        tolower(name)]
    if (length(near)) paste0(" (did you mean ", near[1L], "?)") else ""
  }, character(1))
  stop("control has ", if (length(unknown) > 1L) "entries" else "an entry",
    " the julia sampler does not read: ",
    paste0(unknown, hint, collapse = ", "),
    ". Accepted names: ", paste(sort(.CT_SAMPLE_CONTROL_NAMES),
      collapse = ", "), ".", call. = FALSE)
}

# The sampler settings, from either spelling of the control list.
#
# `ctFit(optimize = FALSE)` took Stan's names for two of these and `ctSample()`
# takes the engine's, so both are read here rather than each entry point
# quietly ignoring what the other documents. The rest are spelled the same on
# both, and the whole list is assembled in one place so that a knob added for
# one cannot go missing from the other -- which is how `minESS` and its three
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
    # A step size for every chain to start from, instead of each estimating
    # its own. `_init_stepsize` doubles or halves from 1 until one leapfrog
    # step under one momentum draw crosses an acceptance of a half, at that
    # chain's own starting point, so its answer differs between chains for
    # reasons that carry no information. Warmup erases that; `warmup = 0` has
    # nothing to erase it with, which is why chains asked for no warmup came
    # back some fast and divergent and others fine.
    #
    # Supplied here rather than shared inside the engine, deliberately: one
    # estimate computed per run there answers differently in a worker process
    # than in this session, because the chunk tuner shifts the last bits of
    # the density and the ladder's threshold turns that into a factor of two.
    stepsize = as.numeric(.ctJuliaOr(control$stepsize, 0)),
    # `adapt_metric` defaults to FALSE, which is a change, and it is measured
    # rather than argued. Warmup re-estimates the metric from 150 iterations
    # up; the metric it replaces is `inv(-H)` for the *exact* Hessian at the
    # mode. Stan shrinks its estimate toward the identity because the identity
    # is all it starts with -- here the starting point is already the right
    # answer, and a few hundred draws cannot improve on it.
    #
    # Three models on dev1, 4 chains x 200 draws, adapt against Laplace-only:
    #
    #   model        warmup   eps (adapt)   eps (laplace)   s (adapt)   s (lap)
    #   informative     200   0.533-0.580   0.708-0.726         393       324
    #   informative     500   0.510-0.668   0.762-0.781         650       533
    #   sparse          200   0.198-0.525   0.613-0.629          99        75
    #   sparse          500   0.243-0.503   0.632-0.699         175       126
    #   wider           200   0.284-0.431   0.616-0.665         218       169
    #   wider           500   0.428-0.443   0.638-0.676         320       275
    #
    # Every cell: a smaller step size, deeper trees, 15-40% more time, and
    # never more effective draws (`sparse` at 200 gave 800 against 595, and
    # `wider` at 500 782 against 688). Below 150 the two are the same run,
    # which is the control -- and those cells agree.
    #
    # It is a setting rather than a removal because the reasoning for it is
    # sound wherever the starting metric is *not* exact, which is any route
    # whose curvature had to be repaired or floored.
    adapt_metric = isTRUE(.ctJuliaOr(control$adapt_metric, FALSE)),
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
  if (!is.null(control$minESS)) settings$min_ess <- as.numeric(control$minESS)
  if (!is.null(control$meanESS)) settings$mean_ess <- as.numeric(control$meanESS)
  if (!is.null(control$maxDraws)) settings$max_draws <- as.integer(control$maxDraws)
  if (!is.null(control$rhatTarget)) settings$rhat_target <- as.numeric(control$rhatTarget)
  if (!is.null(control$settleTol)) settings$settle_tol <- as.numeric(control$settleTol)
  settings
}

# What to sample, in a form that survives a process boundary.
#
# A Julia objective is a handle into one session and cannot be sent anywhere, so
# a worker has to build its own. A chain is therefore described by what that
# rebuild needs -- where to start, how many of those coordinates are parameters,
# which of the engine's two entry points, and whether the objective is the
# state-explicit one -- and `.ctBackendSampleObjective()` turns the description
# back into a handle wherever it lands. The Hessian travels as the plain matrix
# it is, so the workers metre their chains with the parent's rather than each
# spending 2n gradients recomputing it.
#
# `marginal` and `state_explicit` are separate because they are different
# questions. `marginal` says the random effects are integrated out rather than
# sampled, which selects the engine's marginal entry point; `state_explicit`
# says the latent trajectory is sampled alongside the parameters, which changes
# the objective rather than the entry. The second implies the first -- the state
# route cannot be combined with the Laplace effect route at all -- but not the
# reverse.
#' @keywords internal
.ctBackendSampleTarget <- function(estimate, npar, marginal = FALSE,
  state_explicit = FALSE, hessian = NULL, gradient = "adjoint") {
  list(estimate = as.numeric(estimate), npar = as.integer(npar)[1L],
    marginal = isTRUE(marginal), state_explicit = isTRUE(state_explicit),
    hessian = if (is.null(hessian)) NULL else as.matrix(hessian),
    gradient = gradient)
}

# The objective a target names, in whichever process asks for it. Both branches
# are cached per process and keyed on content, so a worker pays for the build
# once and the parent's own handle is reused rather than rebuilt.
#' @keywords internal
.ctBackendSampleObjective <- function(fit, target) {
  if (isTRUE(target$state_explicit)) {
    .ctJuliaJointObjective(fit, target$npar)
  } else {
    .ctJuliaObjective(fit)
  }
}

# Call the engine's sampler, and return what it returned.
#
# The one place either entry point reaches the sampler from, and the one place a
# worker reaches it from too. Everything above it was duplicated until it
# drifted: the two argument lists had diverged over five settings, and the ones
# only `ctFit()` passed were documented on `ctSample()` as though they worked.
#
# `progress` is separate from `verbose` because the two paths decide it
# differently -- a flag on `ctSample()`, and on the fitting path anyone watching
# a console, since sampling there follows an optimisation that has already been
# printing and silence after it reads as a finished run rather than a running
# one.
#' @keywords internal
.ctBackendSampleEngine <- function(fit, target, chains, warmup, draws, cores,
  saveEffects, seed, control, verbose, progress = .ctVerboseOn(verbose),
  callback = control$callback) {

  settings <- .ctBackendSampleControl(control)

  # An effective-size target turns the draw count into a budget.
  #
  # It used to be inert without `maxDraws`: the engine takes `ndraws` and then
  # extends towards `max_draws`, which defaulted to `ndraws` itself, so there
  # was nothing to extend into and the target could only ever be reported after
  # the fact. Someone who set `minESS` watched the sampler run to the end
  # regardless -- reported exactly that way.
  #
  # So the count asked for becomes the budget. `maxDraws` is still there to say
  # "keep going past what I asked for", and given explicitly it wins.
  #
  # The first batch is sized from the *target*, not from the budget, and that
  # distinction is the whole of whether the target does anything. A quarter of
  # the budget was the first rule here and it is useless for a small target: on
  # a 1900-draw budget with `minESS = 100` the first batch was 475 draws, which
  # on five chains is already about 2300 effective ones, so the target was met
  # before it could ever bind and the run returned 29 times what was asked for.
  #
  # Effective size cannot exceed the draws behind it, so `minESS` needs at least
  # `minESS / nchains` draws per chain and there is no point asking for fewer.
  # The floor of 50 is about the estimators rather than the target: split R-hat
  # and effective size read off a few dozen draws are too noisy to stop on, and
  # `rhatTarget` is ANDed with the size target so a batch that cannot support
  # an R-hat estimate just spends a round of scheduling.
  #
  # Which is worth saying plainly, because it is the other half of the surprise:
  # a small `minESS` does not buy a short run. The rule is min ESS *and* mean
  # ESS *and* R-hat, and at a small size target R-hat is what binds -- so a run
  # asked for `minESS = 100` stops when the chains agree, with whatever
  # effective size that took, which is usually far more than 100.
  if (is.null(settings$max_draws) &&
      (!is.null(settings$min_ess) || !is.null(settings$mean_ess))) {
    settings$max_draws <- as.integer(draws)
    # Not `target`: that is this function's own argument, and shadowing it
    # replaced the sampling target with a number, which surfaced as
    # "$ operator is invalid for atomic vectors" from a line nowhere near here.
    ess_target <- max(.ctJuliaOr(settings$min_ess, 0),
      .ctJuliaOr(settings$mean_ess, 0))
    first <- max(50L, ceiling(ess_target / max(1L, as.integer(chains))))
    draws <- max(1L, min(as.integer(draws), as.integer(first)))
  }
  module <- .ctJuliaModule(fit$model_spec$project)
  objective <- .ctBackendSampleObjective(fit, target)

  # Chains are the parallel axis, and in one session they can only be concurrent
  # if it was started with threads for them. Said once, here, because the
  # alternative is a user concluding the sampler is slow when it is running four
  # chains on one thread. Not said in a worker, which runs a single chain.
  if (chains > 1L) {
    threads <- tryCatch(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")),
      error = function(e) NA_integer_)
    if (!is.na(threads) && threads < chains) {
      message("The Julia session has ", threads, " thread(s) and ", chains,
        " chains were asked for, so they will run one after another. ",
        "ctJuliaSetup(threads = ", chains, ", force = TRUE) before fitting ",
        "runs them together.")
    }
  }

  arguments <- list(objective, .ctJuliaNumericVector(target$estimate),
    nchains = as.integer(chains), nwarmup = as.integer(warmup),
    ndraws = as.integer(draws), seed = as.integer(seed)[1L],
    maxdepth = settings$maxdepth, target_accept = settings$target_accept,
    maxdelta = settings$maxdelta, init_scale = settings$init_scale,
    stepsize = settings$stepsize,
    adapt_metric = settings$adapt_metric,
    verbose = isTRUE(progress),
    progress_overwrite = .ctProgressOverwrite(verbose),
    progress_sink = if (isTRUE(progress)) .ctProgressSink(
      .ctProgressOverwrite(verbose)) else NULL)
  for (name in c("min_ess", "mean_ess", "max_draws", "rhat_target", "settle_tol")) {
    if (!is.null(settings[[name]])) arguments[[name]] <- settings[[name]]
  }
  # A live callback into R while the chains run, mirroring
  # `optimcontrol$callback` on `.ctJuliaOptimise()`: a front end that wants to
  # draw sampling progress rather than read it afterwards. Only chain 1 of an
  # in-process multi-chain run ever calls it -- see `sample_run.jl` -- so this
  # is the same "one representative chain" contract the printed line already
  # has, not a second one invented for the callback.
  callback_failure <- NULL
  if (!is.null(callback)) {
    if (!is.function(callback)) {
      stop("control$callback must be a function of (phase, iteration, ",
        "total, logp, divergent).", call. = FALSE)
    }
    # Wrapped exactly as the optimiser's callback is: an error thrown out of
    # an R callback does not reach the engine, it desynchronises the
    # JuliaConnectoR bridge, and a reporting convenience must never be able to
    # take the sample down with it. The message is stored rather than warned
    # immediately, for the same reason as there: `options(warn = 2)` would
    # turn this warning into exactly the error it exists to prevent.
    alive <- TRUE
    arguments$progress_callback <- function(phase, iteration, total, logp,
      divergent) {
      if (alive) {
        tryCatch(callback(phase, iteration, total, logp, divergent),
          error = function(e) {
            alive <<- FALSE
            callback_failure <<- conditionMessage(e)
          })
      }
      NULL
    }
  }
  # `npar`, `save_effects` and `adapt_effects` all describe an effect block the
  # marginal entry does not have, which is the whole structural difference
  # between the two calls.
  if (!isTRUE(target$marginal)) {
    arguments$npar <- as.integer(target$npar)
    arguments$save_effects <- isTRUE(saveEffects)
    arguments$adapt_effects <- settings$adapt_effects
  } else if (!identical(target$gradient, "adjoint")) {
    # `ctsem_sample_marginal`'s own default is `:adjoint`; only said
    # explicitly when something asked for the other one -- currently only a
    # model with a sampled TI predictor value, which forces 'forward' because
    # the adjoint path has no cotangent for it (see
    # `.ctFitJuliaBackendImpl`). `ctsem_sample` (the joint entry, the `!
    # isTRUE(marginal)` branch above) has no such keyword at all, so this is
    # deliberately only reachable on the marginal path.
    arguments$gradient_method <- target$gradient
  }
  if (!is.null(target$hessian)) {
    arguments$hessian <- JuliaConnectoR::juliaPut(as.matrix(target$hessian))
  }
  # How many of the sampled coordinates are model parameters. Only differs
  # from all of them on the state-explicit route, where the vector is
  # `[theta; innovations]` -- and there the engine needs to know, or it builds
  # one dense block over every innovation in the data by inverting an arrow
  # Hessian that is indefinite at the point it is handed. See the note at the
  # metric in `ctsem_sample_marginal`.
  if (isTRUE(target$marginal) && isTRUE(target$state_explicit)) {
    arguments$nparameters <- as.integer(target$npar)
  }
  entry <- if (isTRUE(target$marginal)) module$ctsem_sample_marginal else
    module$ctsem_sample

  result <- .ctBackendWithMaxChunks(cores,
    JuliaConnectoR::juliaGet(do.call(entry, arguments)))
  if (!is.null(callback_failure)) {
    warning("The progress callback failed and was disabled after the first ",
      "error; sampling itself is unaffected. The error was: ",
      callback_failure, call. = FALSE)
  }
  result
}

# Sample, in this session or in one process per chain, and assemble the fit.
#
# The process branch returns a finished fit because the pooling has to happen
# before the assembly can: R-hat and effective sample size are properties of the
# whole run and cannot be averaged from per-chain values. A `NULL` back means
# the workers could not be used and sampling continues here rather than failing
# -- a slower answer beats none.
#
# `processes` is on by default above one chain because the arithmetic is not
# close. A worker costs 26-43 s of Julia startup and engine compilation, against
# a sampling run that is normally minutes to hours -- the startup is noise at
# any realistic draw count, and only dominates on the short runs used for
# testing. What it buys is chains that contend for neither the allocator nor the
# garbage collector.
#' @keywords internal
.ctBackendSampleRun <- function(fit, target, chains, warmup, draws, cores,
  saveEffects, seed, control, verbose, progress = .ctVerboseOn(verbose),
  processes = FALSE, handles = NULL) {

  # What `warmup = 0` costs, said where it is asked for.
  #
  # Warmup is not only a burn-in here: it is the whole of the adaptation. With
  # none, dual averaging never updates, so each chain samples for its whole run
  # at the step size `_init_stepsize` guessed -- doubling or halving from 1
  # until a *single* leapfrog step, under a *single* momentum draw, crosses an
  # acceptance of one half. One step rather than a trajectory, aimed at 0.5
  # rather than at `target_accept`, and evaluated only where the chain starts,
  # so it is a different guess in every chain. Chains then differ in speed and
  # in divergences, which is easy to read as one chain being broken when it is
  # the setting.
  #
  # The metric is not what is lost. It is the Laplace one to begin with, and
  # warmup only re-estimates it at `warmup >= 150` -- 75 of initial buffer, a
  # 25-iteration window, 50 of terminal buffer -- so every shorter warmup
  # already keeps the curvature the fit measured.
  #
  # `target_accept` is what dual averaging aims at, so with no warmup it is not
  # used at all. Said only when the caller set it, because that is the case
  # where something was asked for and is not happening.
  if (warmup < 1L) {
    message("warmup = 0, so nothing adapts: each chain keeps the step size ",
      "found by a single trial leapfrog step, so chains will differ in speed ",
      "and in divergences. The metric is the fit's own curvature either way.",
      if (!is.null(control$target_accept) || !is.null(control$adapt_delta))
        " target_accept only reaches the dual averaging that warmup runs, so it is unused here." else "")
  }

  if (isTRUE(processes) && chains > 1L && .ctBackendCanWarm()) {
    out <- .ctBackendSampleProcesses(fit, target, chains = chains,
      warmup = warmup, draws = draws, cores = cores, handles = handles,
      control = control, saveEffects = saveEffects, seed = seed,
      verbose = verbose, progress = progress)
    if (!is.null(out)) return(out)
    message("Sampling in this session instead.")
  }

  result <- .ctBackendSampleEngine(fit, target, chains = chains,
    warmup = warmup, draws = draws, cores = cores, saveEffects = saveEffects,
    seed = seed, control = control, verbose = verbose, progress = progress)
  # `result$ndraws` rather than the count asked for: with an effective sample
  # size target the sampler decides when to stop, and reporting the request
  # would describe a run that did not happen.
  .ctBackendSampleAssemble(fit, result, target$npar,
    isTRUE(saveEffects) && !isTRUE(target$marginal), as.integer(chains), warmup,
    as.integer(result$ndraws), target$hessian,
    target$estimate[seq_len(target$npar)])
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
  # Named here rather than through `.ctFitNameRawUncertainty()`, because the
  # covariance and the standard errors below are computed *from* this matrix and
  # inherit its names. Same contract, one step earlier.
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
    settings = list(chains = chains, warmup = warmup, draws = draws,
      processes = FALSE))

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
    processes = FALSE,
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
  # The verdict, kept on the fit rather than only shouted once. A warning is
  # seen by whoever is at the console at the time and by nobody afterwards --
  # and `suppressWarnings()` around a fitting call, which any batch script has,
  # removes it entirely. A failed sample has to be readable from the object.
  verdict <- .ctSampleDiagnosis(out$sample)
  out$sample$converged <- verdict$converged
  out$sample$diagnosis <- verdict$problems
  .ctSampleWarn(out$sample)
  out
}

# What went wrong with this sample, as short phrases, and whether the chains
# agreed at all.
#
# `converged` is about whether the draws are a posterior, not about how precise
# one they are: a bad R-hat or a divergent transition means the chains did not
# describe the same distribution, while a small effective sample size or a
# saturated tree depth means they did and did so inefficiently. Only the first
# two set it to FALSE.
#' @keywords internal
.ctSampleDiagnosis <- function(diagnostics) {
  total <- diagnostics$chains * diagnostics$draws
  problems <- character(0)
  divergent <- suppressWarnings(as.integer(diagnostics$divergent)[1L])
  if (!is.na(divergent) && divergent > 0L) {
    problems <- c(problems,
      paste0(divergent, " of ", total, " transitions diverged"))
  }
  flat <- if (is.null(diagnostics$unidentified)) character(0) else
    diagnostics$unidentified
  if (length(flat)) {
    problems <- c(problems, paste0(length(flat),
      " parameter(s) reached a flat region of their transform: ",
      paste(utils::head(flat, 5), collapse = ", ")))
  }
  worst <- suppressWarnings(max(diagnostics$rhat, na.rm = TRUE))
  if (is.finite(worst) && worst > 1.01) {
    problems <- c(problems, paste0("largest R-hat ", signif(worst, 4)))
  }
  fewest <- suppressWarnings(min(diagnostics$ess, na.rm = TRUE))
  if (is.finite(fewest) && fewest < 100) {
    problems <- c(problems,
      paste0("smallest effective sample size ", round(fewest)))
  }
  saturated <- suppressWarnings(as.integer(diagnostics$saturated)[1L])
  if (!is.na(saturated) && saturated > 0L) {
    problems <- c(problems, paste0(saturated, " of ", total,
      " transitions hit the maximum tree depth"))
  }
  # NA, not TRUE, when nothing could be computed: every R-hat is NaN only when
  # every chain sat still, and calling that convergence is the failure mode this
  # whole field exists to avoid.
  converged <- if (!is.finite(worst)) {
    problems <- c(problems, "R-hat could not be computed for any parameter")
    NA
  } else !(worst > 1.01 || (!is.na(divergent) && divergent > 0L))
  list(converged = converged, problems = problems)
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
      " per chain) or draws in ctSample -- or set sampleControl$minESS with ",
      "control$maxDraws to keep sampling until an effective size is reached. ",
      "See fit$sample$rhat.", call. = FALSE)
  }
  fewest <- suppressWarnings(min(diagnostics$ess, na.rm = TRUE))
  if (is.finite(fewest) && fewest < 100) {
    warning("Smallest effective sample size is ", round(fewest), ", from ",
      total, " draws. Interval estimates from this few are unreliable. Raise ",
      "the draw count (iter in ctFit, draws in ctSample), or set ",
      "sampleControl$minESS to keep sampling until an ",
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
  # Last, because it is the conclusion. A reader who stops at the table above
  # has to know what the numbers in it mean; this says it.
  if (!is.null(x$converged)) {
    cat("  chains converged: ", if (is.na(x$converged)) "unknown" else
      as.character(isTRUE(x$converged)), "\n", sep = "")
    if (length(x$diagnosis)) {
      cat("  ", paste(x$diagnosis, collapse = "; "), "\n", sep = "")
    }
  }
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
  optimcontrol, chains, iter, control, priors, intoverpop,
  gradient, verbose, intoverstates = TRUE) {

  # First, because everything after it takes minutes and this takes none.
  .ctBackendSampleCheckControl(control)
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
  start <- .ctJuliaInitialValues(npar, inits,
    initsd = .ctJuliaOr(optimcontrol$initsd, .01))

  # Stan's vocabulary, because these were Stan's arguments: `iter` counts warmup
  # and sampling together and warmup is half of it unless said otherwise.
  #
  # `draws` says the post-warmup count directly, which is what a user means
  # nine times in ten: `iter` and `warmup` together to express "500 draws" is
  # arithmetic nobody should have to do, and getting it wrong silently changes
  # how much of the run is kept. Given, it wins and `iter` is not consulted.
  warmup <- as.integer(.ctJuliaOr(control$warmup, max(1L, floor(iter / 2))))
  draws <- if (!is.null(control$draws)) max(1L, as.integer(control$draws)[1L]) else
    max(1L, as.integer(iter) - warmup)
  seed <- as.integer(.ctJuliaOr(control$seed, 20260828L))
  # `optimcontrol$saveEffects` is where this lived, which was always the wrong
  # list -- it is a sampling setting, and the optimiser has no effects to save.
  # Still read, because scripts pass it.
  saveEffects <- isTRUE(.ctJuliaOr(control$saveEffects,
    optimcontrol$saveEffects))

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
  # The state-explicit sampler is placed from the *integrated* fit, not from
  # the joint mode.
  #
  # Optimising the joint density over parameters and states is the one thing
  # this route must not do: that density has no interior maximum in the state
  # directions, so its "mode" is wherever the optimiser stopped, and the
  # curvature there is not a covariance. Measured on a 12-subject,
  # 10-occasion model, the arrow Hessian at that point had **16 non-positive
  # eigenvalues of 155**, which `_bounded_inverse` floors -- so the metric
  # claimed a standard deviation of about 8.6 along sixteen directions in
  # which the density increases. That is a worse start than no information.
  #
  # So the optimisation and the Hessian both come from the filter, with the
  # states integrated out, which is a proper Laplace approximation of the
  # parameter posterior: positive definite, npar x npar, and already what
  # `intoverstates = TRUE` computes. The innovations then start at zero --
  # their prior mode, and the trajectory the parameters alone imply -- and the
  # engine meters them at identity, which is exactly their prior scale since
  # they are standardised. `nparameters` is what tells it where the parameters
  # stop; see the metric note in `ctsem_sample_marginal`.
  jointobjective <- NULL
  nstate <- 0L
  if (!isTRUE(intoverstates)) {
    jointobjective <- .ctJuliaJointObjective(model_spec, npar)
    nstate <- .ctJuliaStateDimension(model_spec)
  }
  # The workers, started before the optimisation rather than after it. A worker
  # costs 26-43 s of Julia startup and engine compilation, and the optimisation
  # that has to run first anyway is where that cost belongs: measured, workers
  # warmed alongside a 39.8 s optimisation were ready with 0.0 s of waiting.
  # This is the path `.ctBackendWarmWorkers()` was written for. `ctSample()`
  # starts from a fit that is already optimised and so has nothing to overlap,
  # and pays the compile serially.
  #
  # Any finite point compiles the same code, so the pre-optimisation start is as
  # good as the estimate for this -- but only because `ctsem_sample_metric`
  # re-solves the inner modes at the vector being sampled from. Until it did,
  # this warm decided where every process-parallel chain of this path started:
  # the worker's objective retains the modes of whatever it last evaluated, so
  # the chains sampled with theta at the estimate and each subject's effects at
  # the values that were modal for `start`, which put a real 100-subject fit's
  # chains at a joint density of -1e5 and -4e14 against a mode of -5300.
  spec <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  processes <- isTRUE(.ctJuliaOr(control$processes, TRUE))
  handles <- if (processes && chains > 1L && .ctBackendCanWarm()) {
    .ctBackendWarmWorkers(spec, workers = chains, values = start)
  } else NULL

  # `objective = NULL` even on the state route: the integrated objective is
  # what is optimised, for the reason given where `jointobjective` is built.
  optimised <- .ctJuliaOptimise(model_spec, start, optimcontrol = optimcontrol,
    gradient = gradient, cores = cores, verbose = verbose)
  estimate <- as.numeric(optimised$minimizer)

  module <- .ctJuliaModule(model_spec$project)
  # `spec`, not a bare list carrying `model_spec`: .ctBackendHessian() reaches
  # the objective through .ctJuliaObjective(), which requires a classed
  # ctJuliaModel/ctJuliaFit and errors on anything else. An unclassed list made
  # that error every time, and .ctBackendHessian() catches its own errors and
  # returns NULL -- so every sampled fit warned that the engine could not
  # differentiate its gradient, when nothing had been asked of the engine at
  # all. `spec` is the same object the objective is already cached under.
  # The metric's curvature: the integrated parameter Hessian, on every route.
  # The state route used to take the arrow-shaped joint Hessian at the joint
  # mode instead, which is where its metric went wrong.
  hessian <- try(.ctBackendHessian(spec, estimate, verbose = verbose,
    gradient = gradient), silent = TRUE)
  if (inherits(hessian, "try-error")) hessian <- NULL

  # The innovations join the vector after the optimisation rather than before
  # it, at zero: their prior mode, and the trajectory these parameters alone
  # imply. The sampler moves them from there.
  if (!is.null(jointobjective)) estimate <- c(estimate, numeric(nstate))

  joint <- identical(intoverpop, "none")

  # The shell the assembler fills, matching what an optimised fit carries so
  # that everything downstream reads a sampled fit the same way.
  subject_loglik <- as.numeric(optimised$subject_loglik)
  # The population block alone, as everywhere else: `estimate$raw` means the
  # parameters on every fit, and the assembler overwrites it with the
  # posterior mean of exactly those.
  theta <- estimate[seq_len(npar)]
  # `$data` is not set here: see the identical note in .ctFitJuliaBackendImpl()
  # (R/ctJuliaBackend.R) -- ctFit.R attaches the sentinel-cleaned `standata`
  # after this call returns, and `model_spec$data` already covers what code in
  # this package needs from the verbatim long frame.
  out <- list(backend = "julia", model = model, model_spec = model_spec,
    estimate = list(raw = theta,
      loglik = if (length(subject_loglik)) sum(subject_loglik) else
        as.numeric(optimised$maximum_loglik),
      logposterior = as.numeric(optimised$maximum_loglik),
      converged = TRUE, chunks = as.integer(optimised$chunks)),
    engine = model_spec$engine,
    args = list(backend = "julia",
      optimcontrol = optimcontrol, cores = cores, priors = priors,
      intoverpop = intoverpop, optimize = FALSE,
      intoverstates = isTRUE(intoverstates)))
  # "ctFit", not "ctFitModel": the latter marks an unfitted model spec
  # (see .ctFitModelObject), and the optimised path gives its fits "ctFit".
  # Nothing dispatches on either today, so this is consistency rather than a
  # behaviour change.
  class(out) <- c("ctJuliaFit", "ctFit")
  # What the runner needs to know, and all a worker needs to rebuild it: the
  # state-explicit route samples the trajectory alongside the parameters, so
  # `estimate` is longer than `npar` there and the objective is the joint one.
  target <- .ctBackendSampleTarget(estimate = estimate, npar = npar,
    marginal = !joint, state_explicit = !isTRUE(intoverstates),
    hessian = hessian, gradient = gradient)
  out <- .ctBackendSampleRun(out, target, chains = chains, warmup = warmup,
    draws = draws, cores = cores, saveEffects = saveEffects, seed = seed,
    control = control, verbose = verbose,
    # On when someone is watching, matching the optimiser rather than
    # differing from it. The two run one after the other in this same call,
    # and having the first print progress by default while the second stayed
    # silent is what made a running sampler look like a finished
    # optimisation: the visible output stopped at "Computing exact Hessian"
    # and nothing followed it for several minutes.
    progress = verbose > 0L || .ctProgressConsole(),
    processes = processes, handles = handles)
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
    .ctBackendRawParameterNames(out, npar), fit = out, at = theta)
  out
}
