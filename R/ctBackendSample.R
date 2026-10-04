# Sampling a julia backend fit: NUTS, or SAEM's kernel on the joint posterior.
#
# `intoverpop='laplace'` approximates each unit's integral by a Gaussian at its
# mode. That is exact when the integrand is Gaussian in the random effects and
# otherwise wrong by an amount that grows with the population scale, which tilts
# the profile and shrinks the scale estimate -- `ctLaplaceCheck()` measures that
# error and corrects it to first order. Sampling removes it: the joint
# posterior over population parameters *and* random effects makes no Gaussian
# assumption anywhere, and the fit is then used only to place the chains --
# by default through SAEM's state, which starts from it.
#
# `control$target='auto'`, the default, samples the exact posterior the fit's
# route can reach: the joint posterior for an `intoverpop='laplace'` or
# `'none'` fit, and the marginal for an `intoverpop='augmented'` fit, whose
# filter integrates the effects itself. An approximate posterior is worth
# nothing as a sampling target when the exact one is in reach; its job is to
# place the sampler. `target='marginal'` on a Laplace fit samples the Laplace
# marginal all the same -- `npar` dimensions whatever the subject count -- for
# when the approximation is trusted and the joint is too large. Two entry
# points reach the same targets by the same names: `ctFit(optimize=FALSE)` and
# `ctFitUncertainty(fit, 'sample')`. Until 2026-10-03 'auto' meant the Laplace
# marginal on a Laplace fit (decision 5 of
# review/OPTIM-consolidation-plan-2026-09-25.md).
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
# asked instead -- `ctsem_laplace_deviation_layout` reports, for each position,
# which unit and level it belongs to, which subject the block starts at, and
# which of the level's parameters it is. The positions are natural deviations,
# one per parameter a block moves, not the sampler's coordinates: a
# reduced-rank block's coordinates are fewer than its parameters and are no
# parameter's effect.
# That is the same information without the second implementation, and a
# plausible mislabelling would attach the wrong subject's name to a number and
# never announce itself.
#' @keywords internal
.ctBackendEffectIndex <- function(fit) {
  laplace <- fit$model_spec$laplace
  if (is.null(laplace) || is.null(laplace$levels)) return(NULL)
  if (length(laplace$levels) != 1L) return(.ctBackendEffectIndexNested(fit))
  parameters <- .ctSpecRandomEffectLevels(fit$model_spec)[[1L]]$params
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
  layout <- try(.ctJuliaGet(
    .ctJuliaModule(fit$model_spec$project)$ctsem_laplace_deviation_layout(
      .ctJuliaObjective(fit))), silent = TRUE)
  if (inherits(layout, "try-error") || is.null(layout$position)) return(NULL)
  level_index <- as.integer(layout$level)
  # A level the fit does not describe means the two have gone out of step, and
  # an unlabelled effect vector is better than a confidently wrong one.
  if (any(level_index < 1L) || any(level_index > length(laplace$levels))) return(NULL)
  within <- as.integer(layout$within)
  subject <- as.integer(layout$first_member)
  nmembers <- as.integer(layout$nmembers)

  structure <- .ctSpecRandomEffectLevels(fit$model_spec)
  parameter <- vapply(seq_along(within), function(i) {
    pars <- structure[[level_index[i]]]$params
    if (within[i] >= 1L && within[i] <= length(pars)) pars[within[i]] else
      paste0("effect", within[i])
  }, character(1))
  levelname <- vapply(level_index, function(l) {
    nm <- structure[[l]]$name
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

# The grouping identifier a block's members share.
#
# Read off the level structure rather than looked up in the data. The
# specification already carries both halves -- `units`, which unit of a level
# each subject is in, and `labels`, the identifier the user wrote for each unit
# -- and reconstructing them here meant guessing which column of `spec$data`
# held the subject id, which is what once left every group labelled NA. One
# fewer place that has to know how a hierarchy is stored.
#' @keywords internal
.ctBackendGroupIds <- function(fit, laplace, levelname, subject, nmembers) {
  out <- rep(NA_character_, length(subject))
  structure <- .ctSpecRandomEffectLevels(fit$model_spec)
  if (!length(structure)) return(out)
  names(structure) <- vapply(structure, function(x) x$name, character(1))
  for (l in unique(levelname)) {
    level <- structure[[l]]
    if (is.null(level) || is.null(level$labels)) next
    # A block's members share a unit, so the first member's unit names it.
    take <- levelname == l & nmembers > 1L & subject >= 1L &
      subject <= length(level$units)
    if (!any(take)) next
    unit <- level$units[subject[take]]
    out[take] <- ifelse(unit >= 1L & unit <= length(level$labels),
      level$labels[unit], NA_character_)
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

# Sample a fit's posterior by MCMC: `uncertainty = 'sample'`
# on `ctFitUncertainty()`, and what `ctFit(backend = 'julia', optimize =
# FALSE)` calls once its placement optimisation
# (`.ctJuliaOptimiseFit()`/`.ctJuliaSampleFit()`) has produced a fit to
# sample. One runner (`.ctBackendSampleRun()`) for both entry points, so a
# field one of them adds and the other does not is the bug
# `test-julia-fit-shape.R` catches rather than something a caller discovers
# later.
#
# `control$target` names which posterior: `'auto'` (the default) is the exact
# one the fit's route (`.ctBackendIntOverPop()`) can reach -- the joint
# posterior over parameters and random effects for `intoverpop = 'laplace'` or
# `'none'`, the filter's marginal for `'augmented'`. `'marginal'`/`'joint'` ask
# for one explicitly;
# `'joint'` needs the Laplace structure (`intoverpop = 'laplace'` or
# `'none'`) to have somewhere to put the effects, and is refused by name on
# an augmented fit rather than silently sampling the marginal instead.
#
# `state_explicit` is not a `control` entry -- a caller of
# `ctFitUncertainty()` never sets it, because a fit whose states are not
# integrated out does not normally exist to be handed here
# (`optimize = TRUE` with `intoverstates = FALSE` is refused unless
# `optimcontrol$estonly` asked for it). `.ctJuliaSampleFit()` passes it
# explicitly for `ctFit(optimize = FALSE, intoverstates = FALSE)`, and
# `handles` when it has already warmed worker processes alongside the
# placement optimisation that produced `fit`.
#' @keywords internal
.ctBackendUncertaintySample <- function(fit, control = list(), cores = 1L,
  verbose = 0, state_explicit = FALSE, handles = NULL) {

  .ctBackendSampleCheckControl(control)
  target_arg <- .ctJuliaOr(control$target, "auto")
  if (!identical(target_arg, "auto") && !target_arg %in% c("marginal", "joint")) {
    stop("control$target must be 'auto', 'marginal' or 'joint', not '",
      target_arg, "'.", call. = FALSE)
  }
  route <- .ctBackendIntOverPop(fit$model_spec)
  marginal <- switch(target_arg,
    marginal = TRUE,
    joint = FALSE,
    identical(route, "augmented"))
  if (!marginal && is.null(fit$model_spec$laplace)) {
    stop("The joint posterior needs a fit made with intoverpop = 'laplace' ",
      "or 'none': the augmented route carries the random effects in the ",
      "state, so there is no separate posterior over them to sample. Ask ",
      "for control = list(target = 'marginal') to sample the population ",
      "parameters alone.", call. = FALSE)
  }
  placement <- .ctBackendPlacementName(control$placement, marginal || isTRUE(state_explicit))
  if (identical(control$sampler, "saem") && (marginal || isTRUE(state_explicit))) {
    stop("control$sampler = 'saem' samples the joint posterior over population ",
      "parameters and random effects, and this run asks for ",
      if (isTRUE(state_explicit)) "the joint posterior over the latent states. "
      else "the marginal posterior. ",
      "Use the default sampler, or a fit made with intoverpop = 'laplace' or ",
      "'none' and control$target = 'joint'.", call. = FALSE)
  }
  # The sampler, resolved against the target and written back into `control`,
  # which is what the engine call and every worker read it from: the default
  # differs by target, and a worker sees only `control`.
  control$sampler <- .ctBackendSamplerName(control$sampler,
    joint = !marginal && !isTRUE(state_explicit))
  # Resolved once here for its refusals, before any worker is started or the
  # engine reached: a setting refused deep in a worker process would surface
  # as a failed chain and a fallback rather than as the error it is.
  invisible(.ctBackendSampleControl(control))

  npar <- length(fit$estimate$raw)
  chains <- max(1L, as.integer(.ctJuliaOr(control$chains, 4L))[1L])
  warmup <- max(0L, as.integer(.ctJuliaOr(control$warmup, 500L))[1L])
  draws <- max(1L, as.integer(.ctJuliaOr(control$draws, 500L))[1L])
  seed <- as.integer(.ctJuliaOr(control$seed, 20260828L))
  saveEffects <- isTRUE(.ctJuliaOr(control$saveEffects, FALSE))
  processes <- isTRUE(.ctJuliaOr(control$processes, TRUE))

  # Whether the chains get processes is decided in `.ctBackendSampleRun()`,
  # which every caller of this function shares. The one part that cannot move
  # is this message: it depends on whether `control` *named* `processes`.
  # Silent when `future` is simply absent and the default put us here -- that
  # is not the caller's doing and there is nothing to act on -- and said out
  # loud when they asked for processes and are not getting them.
  if (isTRUE(processes) && chains > 1L && !.ctBackendCanWarm() &&
      "processes" %in% names(control)) {
    message("processes = TRUE needs the future package, which is not ",
      "installed. Sampling in this session instead.")
  }

  # The state-explicit target, when asked for: the joint density of the
  # parameters and the innovations that build the latent states, over
  # `[theta; z]`. The innovations join the vector at zero -- their prior
  # mode, and the trajectory the parameters alone imply -- and the engine
  # meters them at identity, exactly their prior scale since they are
  # standardised. `nparameters` (set inside `.ctBackendSampleEngine()`) is
  # what tells it where the parameters stop.
  estimate <- as.numeric(fit$estimate$raw)
  jointobjective <- NULL
  nstate <- 0L
  if (isTRUE(state_explicit)) {
    jointobjective <- .ctJuliaJointObjective(fit$model_spec, npar)
    nstate <- .ctJuliaStateDimension(fit$model_spec)
    estimate <- c(estimate, numeric(nstate))
  }

  # The passed Hessian is a fallback for the joint route, not the metric
  # itself: `_conditional_population_covariance` (sample_nuts.jl) always
  # recomputes the population block fresh, by central differences of the
  # joint gradient at `2 * npar` evaluations, because that is the conditional
  # curvature the joint metric needs and the marginal Hessian passed in is
  # the wrong matrix for it; the fit's Hessian is used only if that
  # difference is not finite.
  sampletarget <- .ctBackendSampleTarget(estimate = estimate, npar = npar,
    marginal = marginal, state_explicit = isTRUE(state_explicit),
    hessian = fit$uncertainty$hessian, placement = placement)

  out <- .ctBackendSampleRun(fit, sampletarget, chains = chains, warmup = warmup,
    draws = draws, cores = cores, saveEffects = saveEffects, seed = seed,
    control = control, verbose = verbose,
    progress = .ctBackendReporting(verbose),
    processes = processes, handles = handles)

  # State-explicit only: the identifiability report is about the parameters,
  # so it needs the curvature that lets the states respond to them -- the
  # *profiled* joint Hessian, by the implicit-function theorem at this
  # (theta, z=0) point -- not the plain marginal `fit` carries (the filter's
  # own linearisation, which the state-explicit route exists to avoid
  # depending on for the answer) and not the corner of the naive joint matrix
  # either, which describes the parameters at a trajectory held fixed. Every
  # other route keeps the identifiability the placement already computed and
  # warned about, which `.ctBackendSampleAssemble()` carried through
  # unchanged above.
  if (!is.null(jointobjective)) {
    theta <- estimate[seq_len(npar)]
    module <- .ctJuliaModule(fit$model_spec$project)
    identhessian <- try(matrix(as.numeric(.ctBackendJuliaValue(
      module$ctsem_joint_hessian(jointobjective,
        .ctJuliaNumericVector(estimate), profile = TRUE))),
      nrow = npar, ncol = npar), silent = TRUE)
    if (inherits(identhessian, "try-error")) identhessian <- NULL
    out$estimate$innovations <- estimate[npar + seq_len(nstate)]
    out$estimate$states <- try(.ctJuliaJointStates(fit$model_spec,
      jointobjective, estimate), silent = TRUE)
    if (inherits(out$estimate$states, "try-error")) out$estimate$states <- NULL
    out$estimate$loglik_type <- "joint"
    out$identifiability <- .ctBackendIdentifiability(identhessian,
      .ctBackendRawParameterNames(out, npar), fit = out, at = theta)
    out$collapsedScales <- .ctBackendCollapsedScales(out)
    .ctBackendIdentifyWarn(out$identifiability, out$collapsedScales, NULL)
  }
  out
}


# Every name the sampler reads out of `control`, in one place.
#
# Kept beside `.ctBackendSampleControl()` because that is what reads most of
# them: a knob added there and not added here is refused as a typo, which is a
# loud failure and the right way round. The five it does not read are read by
# `.ctBackendUncertaintySample()`/`.ctJuliaSampleFit()` (`warmup`, `seed`,
# `processes`, `target`) and `.ctBackendSampleEngine()` (`callback`).
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
  # `target` says which posterior -- 'auto'/'marginal'/'joint' -- read by
  # `.ctBackendUncertaintySample()`, the only caller that resolves it to
  # something other than the fit's own route.
  "iter", "chains", "warmup", "draws", "seed", "saveEffects", "processes",
  "target", "stepsize",
  # Which kernel draws the joint posterior: 'saem' (SAEM's sweeps for the
  # effects, NUTS for the parameters given them; the default there) or 'nuts'
  # (NUTS on the whole joint vector); resolved against the target by
  # `.ctBackendUncertaintySample()` and read by `.ctBackendSampleControl()`.
  "sampler",
  # Where the chains start on the joint target: 'saem' (SAEM's own run and
  # state, the default there) or 'fit' (the fit's estimate); read by
  # `.ctBackendUncertaintySample()`.
  "placement",
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

# The kernel `control$sampler` names. 'saem' runs SAEM's kernel
# (`ctsem_saem_sample`): the effects by SAEM's sweeps, the parameters given
# them by NUTS, placed and stopped exactly as NUTS is. 'nuts' runs NUTS on the
# whole vector. Both draw the same posterior.
#
# 'saem' is the default on the joint target because it is the more reliable
# of the two there, not only the faster: through `ctFitUncertainty()` on 13
# models (dev2, review/POSTERIOR-race-2026-10-02.md) its mean worst error
# against a long reference was 0.34 posterior sds to NUTS's 0.50, at 0.6x the
# time, and NUTS on the joint vector never reached ESS 200 on three variance
# and nested models where it reached it everywhere. A marginal or
# state-explicit target has no effects for its sweeps, so NUTS is the
# sampler there.
#' @keywords internal
.ctBackendSamplerName <- function(sampler, joint = TRUE) {
  if (is.null(sampler)) return(if (isTRUE(joint)) "saem" else "nuts")
  if (!is.character(sampler) || length(sampler) != 1L ||
      !sampler %in% c("nuts", "saem")) {
    stop("control$sampler must be 'nuts' or 'saem'.", call. = FALSE)
  }
  sampler
}

# Where the chains start, `control$placement`. 'saem', the default on the joint
# target, runs SAEM from the fit's estimate (or, under `ctFit(optimize =
# FALSE)`, from the start and prior warm-up alone) and starts the chains from
# its state: its estimate, which targets the exact marginal posterior where the
# fit's optimum targets the Laplace approximation to it, and its chains'
# effects. 'fit' starts them around the fit's own estimate, as before. A
# marginal or state-explicit target has no effects to start, and takes 'fit'.
#' @keywords internal
.ctBackendPlacementName <- function(placement, marginal) {
  if (is.null(placement)) return(if (isTRUE(marginal)) "fit" else "saem")
  if (!is.character(placement) || length(placement) != 1L ||
      !placement %in% c("saem", "fit")) {
    stop("control$placement must be 'saem' or 'fit'.", call. = FALSE)
  }
  if (isTRUE(marginal) && identical(placement, "saem")) {
    stop("control$placement = 'saem' starts the chains from SAEM's draws of the ",
      "random effects, and this run samples the marginal or the latent states, ",
      "which has none to start. Use placement = 'fit'.", call. = FALSE)
  }
  placement
}

# The sampler settings, from either spelling of the control list.
#
# `ctFit(optimize = FALSE)` took Stan's names for two of these and
# `ctFitUncertainty(fit, 'sample')` takes the engine's, so both are read here
# rather than each entry point quietly ignoring what the other documents. The
# rest are spelled the same on both, and the whole list is assembled in one
# place so that a knob added for one cannot go missing from the other -- which
# is how `minESS` and its three companions came to be documented on the
# sampler alone and passed only by `ctFit()`.
#' @keywords internal
.ctBackendSampleControl <- function(control) {
  control <- .ctJuliaOr(control, list())
  settings <- list(
    maxdepth = as.integer(.ctJuliaOr(control$maxdepth,
      .ctJuliaOr(control$max_treedepth, 10L))),
    target_accept = as.numeric(.ctJuliaOr(control$target_accept,
      .ctJuliaOr(control$adapt_delta, 0.8))),
    maxdelta = as.numeric(.ctJuliaOr(control$maxdelta, 1000)),
    # How widely the chains' starts are spread about the placement's centre.
    # Under SAEM's placement (the joint target's default) the spread is one
    # draw from SAEM's complete-data information, which is narrower than the
    # posterior by the information the random effects hide; under
    # `placement = 'fit'` it is one draw from the Laplace approximation.
    #
    # 2, not 1? The Laplace draw has the right shape and the wrong width: Laplace understates spread
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
    adapt_effects = isTRUE(.ctJuliaOr(control$adapt_effects, FALSE)),
    sampler = .ctBackendSamplerName(control$sampler))

  # The SAEM kernel draws the effects by SAEM's own sweeps and keeps the
  # parameters' metric at the conditional curvature it was placed with, so
  # the joint metric's adaptation and the warmup's early stop have nothing to
  # act on there. Refused by name rather than accepted and ignored.
  if (identical(settings$sampler, "saem")) {
    unused <- c(adapt_metric = isTRUE(control$adapt_metric),
      adapt_effects = isTRUE(control$adapt_effects),
      settleTol = !is.null(control$settleTol))
    if (any(unused)) {
      stop("control$", paste(names(unused)[unused], collapse = ", control$"),
        " is a setting of the NUTS sampler's joint metric, which the SAEM ",
        "kernel (control$sampler = 'saem', the default on the joint ",
        "posterior) does not use. Add sampler = 'nuts' to use it.",
        call. = FALSE)
    }
  }

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
  # `.ctEssTarget` (200) by default, the target every draw-producing route
  # shares: "stop once every parameter has 200 effective draws and the chains
  # agree", within a budget of `.ctSampleBudgetMultiple` times the draws asked
  # for (see `.ctBackendSampleBudget()`), and a run that ends short of it says
  # so (`.ctSampleWarn()`).
  #
  # It is the *worst* coordinate that has to reach it, not the mean, which is
  # what makes 200 a defensible floor rather than a loose one -- the 2.5% and
  # 97.5% quantiles `summary()` reports are the part that needs the draws, and
  # they need them for every parameter, not on average. `minESS = 0` turns it
  # off and takes exactly the draws asked for.
  settings$min_ess <- as.numeric(.ctJuliaOr(control$minESS, .ctEssTarget))
  if (!isTRUE(settings$min_ess > 0)) settings$min_ess <- NULL
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
  state_explicit = FALSE, hessian = NULL, gradient = "adjoint", placement = "fit") {
  list(estimate = as.numeric(estimate), npar = as.integer(npar)[1L],
    marginal = isTRUE(marginal), state_explicit = isTRUE(state_explicit),
    hessian = if (is.null(hessian)) NULL else as.matrix(hessian),
    placement = placement,
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

# How many draws per chain a run starts with and may extend to, given the count
# asked for and the effective-size target.
#
# A target turns the draw count into a size to aim at rather than a count to
# take. It used to be inert without `maxDraws`: the engine takes `ndraws` and
# then extends towards `max_draws`, which defaulted to `ndraws` itself, so
# there was nothing to extend into and the target could only ever be reported
# after the fact. Then the count asked for became the budget, which could only
# shorten a run -- and at the defaults it stopped short on most models that
# needed sampling at all: job M's runs (dev2, 2026-09-29) ended at min ESS 100
# (cf6), 70 (gA1p), 165 (gD1p), 107 (gN1p) and 32 (smallp) against the target
# of 200. So the budget is `.ctSampleBudgetMultiple` times the count asked
# for, which those numbers would have reached on all but smallp -- whose
# chains disagreed (R-hat 1.12), where more draws are not the answer and the
# shortfall is warned about instead (`.ctSampleWarn()`). `maxDraws`, given,
# is the budget instead; `minESS = 0` takes exactly the count asked for.
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
.ctSampleBudgetMultiple <- 4L

#' @keywords internal
.ctBackendSampleBudget <- function(draws, chains, settings) {
  draws <- max(1L, as.integer(draws)[1L])
  # Not `target`: that name is the sampling target elsewhere in this file, and
  # shadowing it once replaced the target with a number.
  ess_target <- max(.ctJuliaOr(settings$min_ess, 0), .ctJuliaOr(settings$mean_ess, 0))
  if (!isTRUE(ess_target > 0) || !is.null(settings$max_draws)) {
    return(list(first = draws, max_draws = settings$max_draws))
  }
  first <- max(50L, ceiling(ess_target / max(1L, as.integer(chains))))
  list(first = max(1L, min(draws, as.integer(first))),
    max_draws = .ctSampleBudgetMultiple * draws)
}

# The target a run aimed for, as recorded on the fit (`$sample$ess_target`), or
# NA when it had none.
#' @keywords internal
.ctBackendSampleTargetESS <- function(control) {
  .ctJuliaOr(.ctBackendSampleControl(control)$min_ess, NA_real_)
}

# Call the engine's sampler, and return what it returned.
#
# The one place either entry point reaches the sampler from, and the one place a
# worker reaches it from too. Everything above it was duplicated until it
# drifted: the two argument lists had diverged over five settings, and the ones
# only `ctFit()` passed were documented on the standalone sampler as though
# they worked there too.
#
# `progress` is separate from `verbose` because the two paths decide it
# differently -- a `control` entry on `ctFitUncertainty(fit, 'sample')`, and on
# the fitting path anyone watching
# a console, since sampling there follows an optimisation that has already been
# printing and silence after it reads as a finished run rather than a running
# one.
#' @keywords internal
.ctBackendSampleEngine <- function(fit, target, chains, warmup, draws, cores,
  saveEffects, seed, control, verbose, progress = .ctVerboseOn(verbose),
  callback = control$callback) {

  settings <- .ctBackendSampleControl(control)
  budget <- .ctBackendSampleBudget(draws, chains, settings)
  draws <- budget$first
  settings$max_draws <- budget$max_draws
  module <- .ctJuliaModule(fit$model_spec$project)
  objective <- .ctBackendSampleObjective(fit, target)

  # Chains are the parallel axis, and in one session they can only be concurrent
  # if it was started with threads for them. Said once, here, because the
  # alternative is a user concluding the sampler is slow when it is running four
  # chains on one thread. Not said in a worker, which runs a single chain.
  if (chains > 1L) {
    threads <- tryCatch(as.integer(.ctJuliaEval("Threads.nthreads()")),
      error = function(e) NA_integer_)
    if (!is.na(threads) && threads < chains) {
      message("The Julia session has ", threads, " thread(s) and ", chains,
        " chains were asked for, so they will run one after another. ",
        "ctJuliaSetup(threads = ", chains, ", force = TRUE) before fitting ",
        "runs them together.")
    }
  }

  saem <- identical(settings$sampler, "saem")
  arguments <- list(objective, .ctJuliaNumericVector(target$estimate),
    nchains = as.integer(chains), nwarmup = as.integer(warmup),
    ndraws = as.integer(draws), seed = as.integer(seed)[1L],
    maxdepth = settings$maxdepth, target_accept = settings$target_accept,
    maxdelta = settings$maxdelta, init_scale = settings$init_scale,
    stepsize = settings$stepsize,
    verbose = isTRUE(progress),
    progress_overwrite = .ctProgressOverwrite(verbose),
    progress_sink = if (isTRUE(progress)) .ctProgressSink(
      .ctProgressOverwrite(verbose)) else NULL)
  # The joint metric's adaptation belongs to NUTS on the joint vector; the
  # SAEM kernel has none (`.ctBackendSampleControl()` refuses it there).
  if (!saem) arguments$adapt_metric <- settings$adapt_metric
  for (name in c("min_ess", "mean_ess", "max_draws", "rhat_target",
    if (!saem) "settle_tol")) {
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
    if (!saem) arguments$adapt_effects <- settings$adapt_effects
    if (identical(target$placement, "saem")) arguments$saem <- TRUE
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
    arguments$hessian <- .ctJuliaPut(as.matrix(target$hessian))
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
    if (saem) module$ctsem_saem_sample else module$ctsem_sample

  result <- .ctBackendWithMaxChunks(cores,
    .ctJuliaGet(do.call(entry, arguments)))
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
      "and in divergences. The metric keeps its starting value either way.",
      if (!is.null(control$target_accept) || !is.null(control$adapt_delta))
        " target_accept only reaches the dual averaging that warmup runs, so it is unused here." else "")
  }

  # What this run will do, in one line, before it does it.
  #
  # The effective-size target belongs in it: it decides when the run stops, so
  # a reader who does not know it is there cannot tell a sample that finished
  # early from one that was cut short. Reported wherever progress is --
  # someone watching a console is exactly the someone who needs it.
  if (isTRUE(progress) || .ctProgressConsole()) {
    settings <- .ctBackendSampleControl(control)
    target_ess <- max(.ctJuliaOr(settings$min_ess, 0),
      .ctJuliaOr(settings$mean_ess, 0))
    # Which posterior this run draws from, said once here rather than left to
    # be inferred from `$sample$target` after the fact:
    # `ctFitUncertainty(fit, 'sample')`'s default target depends on the fit's
    # route (joint for 'laplace' and 'none', marginal for 'augmented'), and
    # it changed on 2026-10-03, so it is not always what a reader remembers
    # from an earlier call.
    targetlabel <- if (isTRUE(target$state_explicit))
        "the joint posterior over parameters and the latent states"
      else if (isTRUE(target$marginal))
        "the marginal posterior over population parameters"
      else "the joint posterior over population parameters and random effects"
    budget <- .ctBackendSampleBudget(draws, chains, settings)
    message("Sampling ", targetlabel,
      if (identical(settings$sampler, "saem")) " with the SAEM kernel" else "",
      ": ", chains, " chain",
      if (chains == 1L) "" else "s",
      ", ", warmup, " warmup + ",
      if (target_ess > 0) paste0("up to ", .ctJuliaOr(budget$max_draws, draws))
        else draws,
      " draws each",
      if (target_ess > 0) paste0(", stopping once min ESS ",
        .ctJuliaOr(settings$min_ess, target_ess), " and R-hat ",
        .ctJuliaOr(settings$rhat_target, 1.01), " are reached") else "", ".")
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
    target$estimate[seq_len(target$npar)], marginal = isTRUE(target$marginal),
    ess_target = .ctBackendSampleTargetESS(control))
}

# Turn an engine sample result into a fit object.
#
# Shared by `ctFitUncertainty(fit, 'sample')` and by `ctFit(optimize = FALSE)`,
# which differ only in which engine entry point produced the draws: the joint
# sampler returns
# population parameters and effects, the marginal ones return population
# parameters with the effects already integrated out. Everything after that --
# where the draws go, what becomes the point estimate, which diagnostics warn --
# is the same, and was worth having in one place rather than two that drift.
#' @keywords internal
.ctBackendSampleAssemble <- function(fit, result, npar, saveEffects, chains,
  warmup, draws, hessian, startvalues, marginal = FALSE, ess_target = NA_real_) {

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
  # Where the chains were placed: SAEM's estimate under `placement = 'saem'`,
  # the fit's otherwise. `startvalues` stays the fit's estimate, which is where
  # any Hessian on the fit was taken.
  placed <- if (!is.null(result$placement$theta))
    as.numeric(result$placement$theta)[seq_len(npar)] else as.numeric(startvalues)
  out$estimate$placed_raw <- placed
  out$estimate$raw <- as.numeric(colMeans(posterior))
  out$estimate$cov <- stats::cov(posterior)
  out$estimate$se <- sqrt(diag(out$estimate$cov))
  # `evaluated_at` is the Laplace estimate, not `$estimate$raw`: the exact
  # Hessian was taken at the point the sampler was placed from, and `$raw` is
  # now the posterior mean. Saying so is what stops it being reused as
  # curvature at the mean -- see `.ctBackendHessian()` -- and what lets a
  # reader of the conditional SEs know which point they belong to.
  out$uncertainty <- list(method = "sampling", hessian = hessian,
    evaluated_at = as.numeric(startvalues),
    settings = list(chains = chains, warmup = warmup, draws = draws,
      processes = FALSE))
  # The same for the random-effect check the placement fit made "at the
  # estimate": that estimate was the optimum, and the fit's estimate is now
  # the posterior mean. Its `values` say exactly where; this says it in words.
  if (!is.null(out$identifiability[["effects"]])) {
    out$identifiability$effects$point <- "at the optimum that placed the sampler"
  }

  out$sample <- list(
    chains = chains, warmup = warmup, draws = draws,
    rhat = stats::setNames(as.numeric(result$rhat)[seq_len(npar)], colnames(posterior)),
    ess = stats::setNames(as.numeric(result$ess)[seq_len(npar)], colnames(posterior)),
    # The tails' effective size (`ctsem_sample_diagnostics`): what the 5% and
    # 95% points rest on, and the one a chain that has missed a tail lowers.
    # The verdict, the warning and the run's stopping rule take the worse of
    # the two.
    ess_tail = if (is.null(result$ess_tail)) NULL else
      stats::setNames(as.numeric(result$ess_tail)[seq_len(npar)], colnames(posterior)),
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
    # Which posterior the caller asked for, in the vocabulary
    # `ctFitUncertainty(fit, 'sample', control=list(target=))` and
    # `ctFit(optimize=FALSE)` share -- distinct from
    # `marginal` above, which asks a narrower question (whether the sampled
    # space was exactly `npar`-dimensional) that reads FALSE for both the
    # `'none'` route and the state-explicit one alike.
    target = if (isTRUE(marginal)) "marginal" else "joint",
    # Which kernel drew it (`control$sampler`). Under 'saem' the step size,
    # tree depth, divergences and energy describe the parameters' NUTS
    # transitions given the effects, E-BFMI is not computed (an energy taken
    # at different effects every iteration does not measure it), and the two
    # acceptance rates are the level scale moves' (centred and non-centred).
    sampler = if (is.null(result$sampler)) "nuts" else as.character(result$sampler),
    scale_accept = as.numeric(result$scale_accept),
    ncp_accept = as.numeric(result$ncp_accept),
    processes = FALSE,
    placement = if (is.null(result$placement)) list(method = "fit") else
      list(method = "saem", saem_iterations = as.integer(result$placement$saem_iterations),
        saem_settled = isTRUE(result$placement$saem_settled),
        saem_seconds = as.numeric(result$placement$saem_secs)),
    start = placed)
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
  # The effective size the run aimed for, so the verdict, the warning and the
  # summary judge it against that rather than against a number of their own.
  out$sample$ess_target <- as.numeric(ess_target)[1L]
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
  fewest <- .ctSampleFewest(diagnostics)
  if (is.finite(fewest) && fewest < .ctSampleEssFloor(diagnostics)) {
    problems <- c(problems,
      paste0("smallest effective sample size ", round(fewest),
        if (is.finite(.ctJuliaOr(diagnostics$ess_target, NA_real_)))
          paste0(" (target ", diagnostics$ess_target, ")") else ""))
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

# The smallest effective sample size of a run, bulk or tail: the one the
# engine's stopping rule compared with the target. A fit sampled before tail
# sizes were recorded has only the bulk ones.
#' @keywords internal
.ctSampleFewest <- function(diagnostics) {
  suppressWarnings(min(c(diagnostics$ess, diagnostics$ess_tail), na.rm = TRUE))
}

# The smallest effective sample size a run should end with: the target it
# aimed for (`$ess_target`, `.ctEssTarget` by default), or 100 for a run with
# none -- `minESS = 0`, or a fit sampled before the target was recorded. One
# rule for the verdict, the warning and the summary's note.
#' @keywords internal
.ctSampleEssFloor <- function(diagnostics) {
  target <- suppressWarnings(as.numeric(diagnostics$ess_target)[1L])
  if (length(target) && is.finite(target) && target > 0) target else 100
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
    # together, and `ctFitUncertainty(fit, 'sample')` takes `control$draws`
    # directly.
    # With an effective-size target the run chose its own length within its
    # budget, so the budget is what to raise; without one, the count asked for.
    remedy <- if (is.finite(.ctJuliaOr(diagnostics$ess_target, NA_real_)))
      paste0("aise sampleControl$maxDraws in ctFit (control$maxDraws in ",
        "ctFitUncertainty) past the ", diagnostics$draws, " draws per chain ",
        "this run took. ")
    else paste0("aise the draw count -- iter in ctFit (now ",
      diagnostics$warmup + diagnostics$draws, ", of which ",
      diagnostics$warmup, " is warmup, leaving ", diagnostics$draws,
      " per chain) or control$draws -- or set sampleControl$minESS to keep ",
      "sampling until an effective size is reached. ")
    warning("Largest R-hat is ", signif(worst, 4), ". The chains have not ",
      "agreed on the same distribution, so the draws are not yet a posterior. ",
      if (length(flat))
        "That is expected for the unidentified parameter(s) named above, which more draws cannot fix. For the rest, r"
      else "R",
      remedy, "See fit$sample$rhat.", call. = FALSE)
  }
  fewest <- .ctSampleFewest(diagnostics)
  target <- .ctJuliaOr(diagnostics$ess_target, NA_real_)
  if (is.finite(fewest) && fewest < .ctSampleEssFloor(diagnostics)) {
    # With a target the run stopped at its budget, so the budget is what to
    # raise; without one the draw count asked for was all it was going to take.
    warning("Smallest effective sample size is ", round(fewest), ", from ",
      total, " draws",
      if (is.finite(target)) paste0(", short of the target of ", target,
        ": the sampler stopped at its budget") else "",
      ". Interval estimates rest on that many effective draws. ",
      if (is.finite(target)) paste0("Raise sampleControl$maxDraws in ctFit ",
        "(control$maxDraws in ctFitUncertainty) to let it run longer. ")
      else paste0("Raise the draw count (iter in ctFit, control$draws ",
        "otherwise), or set sampleControl$minESS to keep sampling until an ",
        "effective size is reached. "),
      "See fit$sample$ess and fit$sample$ess_tail.", call. = FALSE)
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
  saem <- identical(x$sampler, "saem")
  cat(if (saem) "ctsem sample, SAEM kernel (effects by SAEM's sweeps, parameters by NUTS given them)\n"
    else "ctsem Hamiltonian sample\n")
  if (!is.null(x$target)) cat("  target: ", x$target, " posterior\n", sep = "")
  cat("  ", x$chains, " chains x ", x$draws, " draws (", x$warmup,
    " warmup discarded)\n", sep = "")
  cat("  divergent: ", x$divergent, " of ", total,
    "   max tree depth reached: ", x$saturated, "\n", sep = "")
  cat("  step size: ", paste(signif(x$stepsize, 3), collapse = ", "),
    "\n", sep = "")
  if (saem) {
    cat("  scale moves accepted, centred / non-centred: ",
      paste(signif(x$scale_accept, 2), collapse = ", "), " / ",
      paste(signif(x$ncp_accept, 2), collapse = ", "), "\n", sep = "")
  } else {
    cat("  E-BFMI:    ", paste(signif(x$ebfmi, 3), collapse = ", "),
      if (any(x$ebfmi < 0.3, na.rm = TRUE)) "  (below 0.3 suggests a funnel)" else "",
      "\n", sep = "")
  }
  worst <- order(-x$rhat)[seq_len(min(5L, length(x$rhat)))]
  cat("  worst R-hat and effective sample size (bulk, tail):\n")
  table <- data.frame(parameter = names(x$rhat)[worst],
    rhat = round(x$rhat[worst], 4), ess = round(x$ess[worst]))
  if (length(x$ess_tail) == length(x$ess)) table$ess_tail <- round(x$ess_tail[worst])
  print(table, row.names = FALSE)
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
# That placement optimisation is `.ctJuliaOptimiseFit()` (R/ctJuliaBackend.R)
# -- the same pipeline `ctFit(optimize=TRUE)` runs: start, prior warm-up,
# substep mesh, the approach, the endgame's certification and its resume, and
# opt-in restarts -- not a lesser one, and not a second fit constructor: the
# placement fit it returns carries `$optim`, `$laplace`, `$uncertainty` and
# `$identifiability` exactly as an optimised fit does, and
# `.ctBackendSampleAssemble()` (shared with `ctFitUncertainty(fit, 'sample')`)
# turns it into the sampled fit by replacing only `$estimate$raw` and adding
# `$sample`. Its
# Hessian is the endgame's own -- reused via `.ctBackendStoredHessian()`
# inside `ctOptimUncertainty()` rather than recomputed here -- so the common
# case costs one Hessian for the whole call, not one for the fit and a second
# for the sampler. `laplace_correct` is forced off for the placement: the
# Laplace fit only places the chains, from its own optimum and curvature, and
# the quadrature-corrected point answers a question the placement does not
# ask.
# review/OPTIM-consolidation-plan-2026-09-25.md P5.
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
.ctJuliaSampleFit <- function(model_spec, datalong, model, prepared_data,
  inits, cores, optimcontrol, chains, iter, control, priors, priorscope,
  intoverpop, gradient, verbose, intoverstates = TRUE) {

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

  # Stan's vocabulary, because these were Stan's arguments: `iter` counts warmup
  # and sampling together and warmup is half of it unless said otherwise.
  # `draws` says the post-warmup count directly and wins when given; see the
  # measurements this cap and default were set from in `dev/` history. 200 is
  # the point at which step-size adaptation has converged and chains agree;
  # capped by half of `iter` so a small `iter` still splits sensibly.
  warmup <- as.integer(.ctJuliaOr(control$warmup,
    max(1L, min(200L, floor(iter / 2)))))
  draws <- if (!is.null(control$draws)) max(1L, as.integer(control$draws)[1L]) else
    max(1L, as.integer(iter) - warmup)
  seed <- as.integer(.ctJuliaOr(control$seed, 20260828L))
  # `optimcontrol$saveEffects` is where this lived, which was always the wrong
  # list -- it is a sampling setting, and the optimiser has no effects to save.
  # Still read, because scripts pass it.
  saveEffects <- isTRUE(.ctJuliaOr(control$saveEffects,
    optimcontrol$saveEffects))
  processes <- isTRUE(.ctJuliaOr(control$processes, TRUE))

  # Whenever the progress line below it will be drawn, not only at `verbose`.
  #
  # Sampling begins with an optimisation, to place the sampler and build its
  # metric, and that optimisation prints its own progress lines (warm-up,
  # approach, endgame) exactly as `optimize=TRUE` does. Without this sentence
  # in front of it a user who asked for HMC watches what looks like an
  # ordinary fit and concludes the sampler never started -- which is exactly
  # what was reported when the placement was a single bare call. The
  # explanation costs one line and was previously hidden behind `verbose > 0`,
  # which is not the default.
  announce <- .ctBackendReporting(verbose)
  if (announce) {
    message("Sampling: placing the sampler first, through the same pipeline ",
      "an optimised fit runs (warm-up, mesh, endgame), which is also where ",
      "its metric comes from. Sampling follows.")
  }

  # The workers, started before the placement rather than after it. A worker
  # costs 26-43 s of Julia startup and engine compilation, and the placement
  # that has to run first anyway is where that cost belongs: measured, workers
  # warmed alongside a 39.8 s optimisation were ready with 0.0 s of waiting.
  # `ctFitUncertainty(fit, 'sample')` starts from a fit that is already
  # optimised and so has nothing to overlap, and pays the compile serially.
  #
  # Any finite point compiles the same code, so a pre-placement start is as
  # good as the placement's own estimate for this -- only compilation is being
  # bought here, not a value the chains will actually start from.
  #
  # Any substep mesh does too, but the maxtimestep rule does not: the engine's
  # subject objective is typed on its substep policy, a `Float64` rule or a
  # per-row `Vector{Int}` mesh. Under `nsubsteps = 'auto'` the placement
  # replaces the rule with a mesh and the chains filter with that, so a worker
  # warmed at the rule compiled code no chain calls, then compiled the filter
  # and its gradient again when its chain began. Ones are a valid mesh on any
  # rows. Only a start where the likelihood is not finite keeps the rule, and
  # pays that.
  start0 <- .ctJuliaInitialValues(npar, inits,
    initsd = .ctJuliaOr(optimcontrol$initsd, .01))
  spec0 <- structure(model_spec, class = c("ctJuliaModel", "ctFitModel"))
  if (!is.null(spec0$substeps)) spec0$max_timestep <- rep(1L, length(spec0$times))
  handles <- if (processes && chains > 1L && .ctBackendCanWarm()) {
    .ctBackendWarmWorkers(spec0, workers = chains, values = start0)
  } else NULL

  # Placement: `.ctJuliaOptimiseFit()` (R/ctJuliaBackend.R) is the whole
  # pipeline `ctFit(optimize=TRUE)` runs -- start, prior warm-up, substep
  # mesh, the approach, the endgame's certification and its resume, restarts
  # only if asked -- not a bare `.ctJuliaOptimise()` call. `intoverstates` is
  # forced TRUE here regardless of what was asked: optimising the joint
  # density of the parameters and the latent states is a category error, not
  # a lesser estimate (`state-explicit-generation.md`) -- the joint mode is
  # degenerate and not a place to start a sampler from. `laplace_correct` is
  # forced off and `uncertainty` forced to `'hessian'`: under
  # `intoverpop='laplace'` the fit only places the chains, from its optimum
  # and curvature -- not the quadrature-corrected point, which answers a
  # question the placement is not asking -- so only the Hessian is needed
  # here, not importance draws. `finishsamples = 2` (the least
  # `ctFitUncertainty()` accepts) for the same reason: it would otherwise
  # draw a thousand Gaussian pseudo-posterior samples around the placement
  # point only for `.ctBackendSampleAssemble()` to overwrite them with the
  # real ones below.
  placementcontrol <- utils::modifyList(optimcontrol,
    list(laplace_correct = FALSE, uncertainty = "hessian", estonly = FALSE,
      finishsamples = 2L))
  # When SAEM places the chains (`control$placement = 'saem'`, the default on
  # the joint target), it runs from this fit's point, so this fit needs only
  # the start and the prior warm-up: the Laplace approximation's optimum is
  # not where the exact posterior is, and SAEM is what finds that. No
  # Hessian, so no identifiability report from a placement optimum either; the
  # sampler's own flat-region check stands in.
  jointtarget <- isTRUE(intoverstates) && switch(.ctJuliaOr(control$target, "auto"),
    joint = TRUE, marginal = FALSE, intoverpop %in% c("laplace", "none"))
  if (jointtarget && identical(.ctBackendPlacementName(control$placement, FALSE), "saem")) {
    placementcontrol <- utils::modifyList(placementcontrol,
      list(estonly = TRUE, maxiter = 0L))
  }
  placementfit <- .ctJuliaOptimiseFit(model_spec = model_spec, datalong = datalong,
    model = model, prepared_data = prepared_data, inits = inits, cores = cores,
    optimcontrol = placementcontrol, verbose = verbose, priors = priors,
    priorscope = priorscope, intoverpop = intoverpop, intoverstates = TRUE,
    gradient = gradient, correctlaplace = FALSE)

  # From here, sampling is exactly `ctFitUncertainty(placementfit, 'sample')`:
  # one function for both entry points (`.ctBackendUncertaintySample()`,
  # above), so a field it adds is on both routes and a bug in it is one bug,
  # not two to keep in step. `chains`/`warmup`/`draws`/`seed`/`saveEffects`/
  # `processes` were resolved above from `ctFit()`'s own `iter`/`optimcontrol`
  # vocabulary; folded into `control` here so the shared function reads them
  # exactly as it reads a caller's own `control` list. `target` is left at
  # `.ctBackendUncertaintySample()`'s default ('auto'): its own route
  # inference from `.ctBackendIntOverPop(placementfit$model_spec)` already
  # gives the joint posterior for `intoverpop='none'` and the marginal
  # otherwise, which is what `intoverpop` says here too.
  samplecontrol <- utils::modifyList(control,
    list(chains = chains, warmup = warmup, draws = draws, seed = seed,
      saveEffects = saveEffects, processes = processes))
  .ctBackendUncertaintySample(placementfit, control = samplecontrol,
    cores = cores, verbose = verbose, state_explicit = !isTRUE(intoverstates),
    handles = handles)
}
