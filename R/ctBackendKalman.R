# Prediction and Kalman output for backend='julia' ---------------------------
#
# `ctKalman()` / `ctPredict()` read four arrays: prior, filtered and smoothed
# states and observations for every data row. The engines now produce those from
# the same forward pass they use for the likelihood (see `kalman.hpp` and
# `kalman_trace.jl`), so everything downstream -- unpacking, standardised
# residuals, naming, melting, plotting -- is the code Stan fits already use,
# reached through `.ctKalmanArrayAssemble()`.
#
# What is genuinely backend-specific is only the two steps around that:
#
#   * getting the arrays out of the engine (`.ctBackendKalmanRaw`), and
#   * preparing the data to filter, when the caller wants a subset of subjects,
#     a finer time grid than the data has, or observations withheld
#     (`.ctBackendKalmanSpec`).
#
# The second is done by rebuilding the *data*, not by editing a prepared
# structure. Stan interpolates by inserting rows into `standata` and remerging;
# the engines take a `datalong`, so inserting the rows there and re-preparing is
# both simpler and impossible to get subtly out of step with what a fresh fit to
# the same frame would do.

# Whether `fit` is a julia backend fit. `model_spec` is set only by the julia
# path (see .ctFitJuliaBackendImpl / .ctJuliaSampleFit in R/ctJuliaBackend.R
# and R/ctBackendSample.R), so its presence is a reliable backend test -- and,
# since a julia fit now also carries `$standata` (for `$data`/`$standata`
# parity with a stan fit, see R/ctFit.R), no longer one that presence of
# `$standata` can make. The functions below use this to keep a julia fit
# routed through its own `model_spec`, which is what stays in the row order
# the engine actually used: ctPostPredData() zips several of these accessors
# together by row position, and mixing one that reads `$standata` with others
# that read `model_spec` would misalign rows the moment the two disagree on
# ordering, silently.
.ctFitIsJulia <- function(fit) !is.null(fit$model_spec)

# The original-id to internal-index mapping, for the backends that do not carry
# a `standata`. ctPredict() speaks in the user's own subject ids and the filter
# speaks in positions, so something has to hold the correspondence.
.ctFitIdMap <- function(fit) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata$idmap)) return(fit$standata$idmap)
  spec <- .ctBackendSpec(fit)
  ids <- unique(spec$data[[.ctFitModelObject(fit)$subjectIDname]])
  data.frame(original = ids, new = seq_along(ids), stringsAsFactors = FALSE)
}

# Which levels of random effect a Laplace trajectory is built from.
#
# Naming a level means "this level and every level outside it": the innermost
# id gives each subject its own trajectory, an outer id gives the group mean --
# every subject in a study sharing one -- and 'population' gives the trajectory
# implied by the population parameters alone. The rule is the same one
# `_laplace_restrict_levels` applies, stated once here in the user's terms.
.ctBackendLaplaceLevel <- function(spec, effects) {
  if (is.null(spec$laplace)) return(1L)
  names <- vapply(.ctSpecRandomEffectLevels(spec), function(x) x$name,
    character(1))
  if (identical(effects, "population")) return(length(names) + 1L)
  position <- match(effects, names)
  if (is.na(position)) {
    stop("randomEffects must be 'population' or one of the model's id names (",
      paste(names, collapse = ", "), "), not '", effects, "'.", call. = FALSE)
  }
  as.integer(position)
}

# What a Laplace fit's trajectories are conditional on, for the calls that
# return trajectories. Silent for every other kind of fit.
.ctBackendLaplaceTrajectoryNote <- function(spec, randomEffects) {
  if (is.null(spec$laplace)) return(invisible(NULL))
  if (is.null(randomEffects)) randomEffects <- .ctSpecRandomEffectLevels(spec)[[1L]]$name
  # Validates the name here rather than leaving a bad one to error further in,
  # where the message would be about a level index.
  .ctBackendLaplaceLevel(spec, randomEffects)
  # A Laplace fit's effects are conditional modes; a 'none' fit samples them,
  # and calling a draw a mode would be a different claim about the same number.
  sampled <- identical(.ctBackendIntOverPop(spec), "none")
  message(if (sampled) "Sampled random effects: " else "Laplace fit: ",
    "trajectories are conditional on random effects ",
    if (sampled) "drawn for each subject" else
      "estimated from each subject's whole record",
    ", so they are the smoothed equivalent rather than filtered. ",
    "randomEffects='", as.character(randomEffects), "'.")
  invisible(NULL)
}

.ctBackendKalmanRaw <- function(fit, raw, subjectmatrices = TRUE,
  randomEffects = NULL, fields = NULL) {
  raw <- as.numeric(raw)
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  arguments <- list(.ctJuliaObjective(fit), .ctJuliaNumericVector(raw),
    subject_matrices = isTRUE(subjectmatrices))
  if (length(fields)) arguments$fields <- .ctJuliaVector(as.character(fields))

  # The Laplace route is the same call with two more keywords: which levels of
  # random effect to build the trajectory from, and -- when this specification
  # filters different rows than the fit did -- the fitted modes to use. The
  # engine dispatches on the objective type, so `from_level` must not be sent
  # on the augmented route at all, which is what the branch is for; everything
  # else about the call is shared and used to be written out twice.
  if (!is.null(spec$laplace)) {
    if (is.null(randomEffects)) randomEffects <- .ctSpecRandomEffectLevels(spec)[[1L]]$name
    from <- .ctBackendLaplaceLevel(spec, randomEffects)
    arguments$from_level <- as.integer(from)
    # Modes from the fit, not from whatever rows this call happens to filter
    # over. `.ctBackendKalmanSpec` attaches the fitted specification when it
    # rebuilds; without one, this specification *is* the fit's and there is
    # nothing to carry.
    source <- attr(fit, "laplaceSource")
    if (!is.null(source)) {
      subjects <- attr(fit, "laplaceSubjects")
      fitted <- .ctBackendJuliaValue(module$ctsem_laplace_subject_values(
        .ctJuliaObjective(source), .ctJuliaNumericVector(raw),
        from_level = as.integer(from)))
      arguments$subject_values <- JuliaConnectoR::juliaPut(
        fitted[subjects, , drop = FALSE])
    }
  }

  result <- .ctBackendJuliaValue(do.call(module$ctsem_kalman, arguments))
  if (!is.null(result$subject)) result$subject <- as.integer(result$subject)
  result
}

# Whether this fit's random effects are separate coordinates rather than
# augmented latent states.
#
# `spec$laplace` is named after one of the two methods that use it and is built
# for both: `.ctJuliaPrepare()` prepares 'none' exactly as it prepares
# 'laplace', because what that structure holds is the *description* of the
# random effects -- which raw parameters vary, at which level, with which
# population scale -- and that description is needed whether they are then
# integrated out or sampled. So the presence of the field says how the effects
# are represented and not what the fit did with them, and the two questions
# have to be asked separately.
.ctSpecEffectsAreCoordinates <- function(spec) !is.null(spec$laplace)

# How a fit integrates its random effects, so that re-preparing a specification
# keeps doing what the fit did. `.ctJuliaPrepare` defaults to "augmented",
# which for a Laplace fit is not a slower path to the same answer -- it is a
# different model, with the random effects carried as latent states.
#
# Read from the specification rather than derived from it. Deriving gave
# "laplace" for every fit carrying the structure above, which includes every
# `intoverpop='none'` fit -- the NUTS-over-parameters-and-effects route a
# sampled fit with random effects takes (`ctFit(optimize=FALSE)`, see the table
# in R/ctFit.R). Nothing breaks today, because 'none' and 'laplace' prepare
# identically and differ only in what happens afterwards, so re-preparing one
# as the other rebuilds the same structure. It is correct by coincidence
# though, and this function's whole job is to say what the fit did.
.ctBackendIntOverPop <- function(spec) {
  recorded <- spec$intoverpop
  if (!is.null(recorded) && nzchar(as.character(recorded)[1L])) {
    return(as.character(recorded)[1L])
  }
  if (.ctSpecEffectsAreCoordinates(spec)) "laplace" else "augmented"
}

# Rebuild the prepared specification over the subjects, times and observations a
# prediction call asks for.
.ctBackendKalmanSpec <- function(fit, subjects = "all", timestep = "asdata",
  maxtime = "asdata", removeObs = FALSE) {
  spec <- .ctBackendSpec(fit)
  # Nothing to change: hand back the specification the fit already has, rather
  # than rebuilding an identical one. That is not only faster (it keeps the
  # cached engine objective) -- re-preparation is the one step in this file that
  # could differ from what the fit itself did, so not doing it when it buys
  # nothing is the safer default.
  if (identical(subjects, "all") && identical(timestep, "asdata") &&
      (identical(removeObs, FALSE) || identical(removeObs, 0))) {
    return(.ctBackendAsModel(spec))
  }
  model <- .ctFitModelObject(fit)
  idname <- model$subjectIDname
  timename <- model$timeName
  dat <- spec$data
  ids <- unique(dat[[idname]])

  if (!identical(subjects, "all")) {
    wanted <- if (is.numeric(subjects)) ids[subjects] else subjects
    wanted <- wanted[!is.na(wanted)]
    if (!length(wanted)) stop("No subjects selected.", call. = FALSE)
    dat <- dat[dat[[idname]] %in% wanted, , drop = FALSE]
  }

  if (!identical(timestep, "asdata")) {
    if (!isTRUE(model$continuoustime)) {
      stop("Discrete time model fits must use timestep = 'asdata'.", call. = FALSE)
    }
    step <- if (identical(timestep, "auto")) {
      gaps <- diff(dat[[timename]])
      gaps <- gaps[gaps > 0]
      if (!length(gaps)) stop("Cannot choose a timestep automatically.", call. = FALSE)
      stats::median(gaps) / 5
    } else as.numeric(timestep)[1L]
    if (!is.finite(step) || step <= 0) stop("timestep must be a positive number.", call. = FALSE)
    upper <- if (identical(maxtime, "asdata")) max(dat[[timename]]) else max(as.numeric(maxtime))
    grid <- seq(min(dat[[timename]]), upper, by = step)
    dat <- .ctBackendFillTime(dat, grid, idname, timename, model$manifestNames,
      model$TDpredNames)
  }

  dat <- dat[order(match(dat[[idname]], unique(dat[[idname]])), dat[[timename]]), ,
    drop = FALSE]

  # removeObs withholds observations from the *filter*, not from the report:
  # the point of it is to compare a prediction against the observations it was
  # not given, so what comes back still carries the data. TRUE withholds
  # everything, a positive integer N keeps every Nth. Covariates are untouched
  # either way. This runs after the time grid is built so that the reported
  # observations line up with the rows that were actually filtered.
  reported <- t(as.matrix(dat[, model$manifestNames, drop = FALSE]))
  withhold <- isTRUE(removeObs) || (is.numeric(removeObs) && removeObs > 0)
  if (withhold) {
    keep <- if (is.numeric(removeObs) && removeObs > 1) {
      unlist(lapply(split(seq_len(nrow(dat)), dat[[idname]]),
        function(rows) rows[seq(1, length(rows), by = as.integer(removeObs))]))
    } else integer()
    dat[setdiff(seq_len(nrow(dat)), keep), model$manifestNames] <- NA
  }

  # With the fit's integration step. Re-preparing put a mesh from
  # `nsubsteps = 'auto'` back to the maxtimestep rule, so every prediction that
  # selected subjects or interpolated a grid -- ctPredict()'s defaults do both
  # -- was integrated more coarsely than the fit.
  prepared <- .ctBackendAsModel(.ctJuliaCarrySubsteps(.ctJuliaPrepare(dat, model,
    project = spec$project, intoverpop = .ctBackendIntOverPop(spec),
    laplacecontrol = spec$laplace$inner), spec, .ctFitRawEstimate(fit)))
  if (withhold) attr(prepared, "reportManifest") <- reported

  # A Laplace fit's random effects are conditional modes, and a mode is only
  # defined relative to the data it was estimated from. This specification
  # filters over different rows than the fit did, so re-solving here would
  # answer a different question -- and with `removeObs` it would answer none at
  # all, because a subject with no observations left has nothing to condition
  # on and its mode collapses to zero. The fitted specification is carried
  # along so the modes can be taken from it instead, together with the map from
  # this specification's subjects back to its own.
  if (!is.null(spec$laplace)) {
    attr(prepared, "laplaceSource") <- .ctBackendAsModel(spec)
    attr(prepared, "laplaceSubjects") <- match(
      unique(dat[[idname]]), unique(spec$data[[idname]]))
  }
  prepared
}

.ctBackendAsModel <- function(spec) {
  structure(spec, class = c("ctJuliaModel", "ctFitModel"))
}

# Extra rows at the requested times, with every manifest missing and every TD
# predictor zero, so the filter produces an expectation there without being told
# anything about it. A subject's own first observation is kept as its t0 rather
# than being displaced by a grid point below it, matching `standataFillTime`.
.ctBackendFillTime <- function(dat, grid, idname, timename, manifestNames, TDpredNames) {
  added <- lapply(split(seq_len(nrow(dat)), dat[[idname]]), function(rows) {
    block <- dat[rows, , drop = FALSE]
    present <- round(block[[timename]], 10)
    wanted <- grid[!round(grid, 10) %in% present & grid > min(block[[timename]])]
    if (!length(wanted)) return(NULL)
    filler <- block[rep(1L, length(wanted)), , drop = FALSE]
    filler[[timename]] <- wanted
    filler[, manifestNames] <- NA
    if (length(TDpredNames)) filler[, TDpredNames] <- 0
    filler
  })
  added <- added[!vapply(added, is.null, logical(1))]
  if (!length(added)) return(dat)
  rbind(dat, do.call(rbind, added))
}

# Raw parameter draws for a prediction call, in ctKalmanArray's terms.
.ctBackendKalmanSamples <- function(fit, pointest, nsamples, collapsefunc, ...) {
  if (isTRUE(pointest)) return(matrix(as.numeric(fit$estimate$raw), nrow = 1L))
  samples <- .ctBackendRawSamples(fit)
  if (!is.na(nsamples) && nsamples < nrow(samples)) {
    samples <- samples[sample(seq_len(nrow(samples)), nsamples), , drop = FALSE]
  }
  if (is.function(collapsefunc)) {
    samples <- matrix(apply(samples, 2, collapsefunc, ...), nrow = 1L)
  }
  samples
}

#' Kalman filter and smoother estimates from a julia backend fit
#'
#' Prior, filtered and smoothed estimates of the latent states and the
#' observations, for every row of data, from the same forward pass the engine
#' uses for the likelihood.
#' \code{\link{ctKalmanArray}} (alias \code{ctStanKalman}) dispatches here for
#' a \code{ctJuliaFit}, so call that unless you want the julia route
#' specifically. \code{\link{ctPredict}} (alias \code{ctKalman}) returns the
#' same information as a long data frame.
#'
#' Internal. \code{ctBackend*} means internal throughout the package, so this is
#' not exported; users reach it through \code{\link{ctKalmanArray}} and
#' \code{\link{ctPredict}}, which dispatch here for a \code{ctJuliaFit}. The
#' \code{randomEffects} notes below describe behaviour those two inherit.
#'
#' @param fit A \code{ctJuliaFit}.
#' @param subjects \code{'all'}, a vector of subject ids, or integer positions
#'   into the fitted subjects.
#' @param timestep \code{'asdata'} to use the observed times, \code{'auto'} to
#'   choose an interpolation step, or a positive number.
#' @param maxtime Only used when interpolating: the largest time to compute for.
#' @param removeObs \code{TRUE} withholds every observation so that only
#'   expectations given parameters and covariates are returned; a positive
#'   integer N keeps every Nth.
#' @param pointest Use the point estimate as the single sample.
#' @param nsamples Number of posterior draws to use when \code{pointest=FALSE}.
#' @param collapsefunc Function applied over draws, e.g. \code{mean}.
#' @param standardisederrors Also return residuals standardised by the prior
#'   observation covariance.
#' @param subjectpars Also return each subject's own model matrices.
#' @param indvarstates Keep the augmented individual-difference states in the
#'   latent output rather than trimming to the real processes.
#' @param ... Passed to \code{collapsefunc}.
#' @return A list of arrays, the same shape \code{\link{ctKalmanArray}} returns
#'   for Stan fits, with \code{subjectMatrices} attached when
#'   \code{subjectpars=TRUE}.
#' @details Individual-difference parameters come out of the smoothed initial
#'   state: ctsem carries an individually varying parameter as an augmented
#'   latent state with no drift and no diffusion, so its smoothed t0 estimate is
#'   the subject's value for it.
#' @examples
#' \donttest{
#' # ctBackendKalman(fit, timestep = .1)
#' }
#' @param randomEffects For an `intoverpop='laplace'` fit, which levels of
#'   random effect the trajectories are built from. Naming one of the model's
#'   id columns includes that level and every level outside it, so the
#'   innermost id gives each subject its own trajectory and an outer id gives
#'   the group mean, shared by every subject in the group. `'population'`
#'   includes none, giving the trajectory implied by the population parameters
#'   alone. Defaults to the innermost id. Ignored for other fits.
#'
#'   These are the *smoothed* equivalent whichever level is chosen: the random
#'   effects are modes estimated from each subject's whole record, so unlike an
#'   augmented fit -- whose carrier states are updated observation by
#'   observation -- there is no version of them from before a subject's later
#'   data arrived. A message says so once per call.
#'
#'   The modes are the ones the fit arrived at, and they stay fixed however
#'   this call changes the rows being filtered. Selecting subjects,
#'   interpolating a time grid or withholding observations with
#'   \code{removeObs} therefore leaves each subject's parameters alone: a
#'   prediction with every observation withheld still uses that subject's own
#'   random effects, which is what makes it a prediction *for that subject*
#'   rather than for the average one.
#' @keywords internal
ctBackendKalman <- function(fit, subjects = "all", timestep = "asdata",
  maxtime = "asdata", removeObs = FALSE, pointest = TRUE, nsamples = NA,
  collapsefunc = NA, standardisederrors = FALSE, subjectpars = FALSE,
  indvarstates = FALSE, randomEffects = NULL, ...) {

  # Same hazard as .ctBackendGenerateFromFit's statepath branch further down
  # this file: a state-explicit fit's point estimate is the mode of the joint
  # density of states and data, not of the marginal the filter below computes.
  # Filtering through it at that estimate is a draw from a density this fit
  # did not maximise, and it looks entirely reasonable.
  if (isFALSE(fit$args$resolved$intoverstates)) {
    warning('Kalman filter operation unreliable when states were sampled -- system noise represents prior while point estimates represent posterior / smoothed')
  }

  spec <- .ctBackendKalmanSpec(fit, subjects = subjects, timestep = timestep,
    maxtime = maxtime, removeObs = removeObs)
  # Said here, at the one entry point that returns trajectories, because the
  # difference is easy to miss and changes what the picture means. An
  # augmented fit's random effects are carrier *states*, updated observation by
  # observation, so its filtered output shows an effect being learned. A
  # Laplace mode is estimated from all of a subject's data at once, so these
  # trajectories are the smoothed equivalent throughout -- there is no "before
  # this subject's later data arrived" version of them.
  #
  # Not in `.ctBackendKalmanRaw`, where it used to be. Every julia fit runs a
  # filter pass for its prior residuals, `ctLOO` runs one per fold and the
  # summary runs one per posterior draw -- so a note about how to read a
  # trajectory was printed by things that return no trajectory, and every
  # Laplace fit said it whether or not anyone asked for a prediction.
  .ctBackendLaplaceTrajectoryNote(spec, randomEffects)
  model <- .ctFitModelObject(fit)
  samples <- .ctBackendKalmanSamples(fit, pointest, nsamples, collapsefunc, ...)

  nrows <- length(spec$times)
  naug <- spec$nlatent_augmented
  nmanifest <- length(model$manifestNames)
  iterations <- nrow(samples)

  etaa <- array(0, dim = c(iterations, 3L, nrows, naug))
  etacova <- array(0, dim = c(iterations, 3L, nrows, naug, naug))
  ya <- array(0, dim = c(iterations, 3L, nrows, nmanifest))
  ycova <- array(0, dim = c(iterations, 3L, nrows, nmanifest, nmanifest))
  llrow <- matrix(0, iterations, nrows)
  subjectmatrices <- NULL
  id <- NULL

  for (iteration in seq_len(iterations)) {
    scores <- .ctBackendKalmanRaw(spec, samples[iteration, ],
      subjectmatrices = isTRUE(subjectpars),
      # Every field this loop reads, and nothing else -- in particular not
      # `transition`, which nothing here or downstream looks at.
      fields = c("eta", "etacov", "y", "ycov", "llrow", "subject", "subject_loglik"),
      randomEffects = randomEffects)
    etaa[iteration, , , ] <- scores$eta
    etacova[iteration, , , , ] <- scores$etacov
    ya[iteration, , , ] <- scores$y
    ycova[iteration, , , , ] <- scores$ycov
    llrow[iteration, ] <- scores$llrow
    if (is.null(id)) id <- scores$subject
    # A subject whose filter failed leaves its smoothed rows untouched, so this
    # has to be said rather than returned as zeros that look like estimates.
    failed <- which(!is.finite(scores$subject_loglik))
    if (length(failed)) {
      warning("The filter failed for subject(s) ", paste(failed, collapse = ", "),
        "; their smoothed estimates are not usable.", call. = FALSE)
    }
    if (isTRUE(subjectpars)) {
      if (is.null(subjectmatrices)) {
        subjectmatrices <- array(0, dim = c(iterations, dim(scores$subject_matrices)))
      }
      subjectmatrices[iteration, , ] <- scores$subject_matrices
    }
  }

  latentNames <- model$latentNames
  nlatent <- spec$nlatent
  if (isTRUE(indvarstates)) {
    nlatent <- naug
    if (naug > length(latentNames)) {
      latentNames <- c(latentNames, .ctBackendCarrierNames(spec, naug - length(latentNames)))
    }
  }

  # The manifests as reported, which is not the same as the manifests filtered
  # when removeObs withheld some (see .ctBackendKalmanSpec).
  reported <- attr(spec, "reportManifest")
  Y <- t(if (is.null(reported)) spec$manifest_data else reported)
  Y[!is.finite(Y)] <- NA
  colnames(Y) <- model$manifestNames

  out <- .ctKalmanArrayAssemble(
    list(etaa = etaa, etacova = etacova, ya = ya, ycova = ycova, llrow = llrow),
    time = spec$times, Y = Y, id = id, nlatent = nlatent, latentNames = latentNames,
    manifestNames = model$manifestNames, standardisederrors = standardisederrors)

  if (isTRUE(subjectpars)) {
    out$subjectMatrices <- .ctBackendSubjectMatrices(spec, subjectmatrices)
  }
  out
}

# Names for the augmented carrier states, which are the individually varying
# parameters themselves.
.ctBackendCarrierNames <- function(spec, count) {
  effects <- spec$random_effects
  # Only the states beyond the real latent processes are carriers: a random
  # effect on T0MEANS varies a genuine state, and has a `random_effects` row
  # too, so taking the first `count` rows named the wrong parameters whenever
  # both kinds were present.
  if (!is.null(effects) && length(effects) && !is.null(effects$param) &&
      !is.null(spec$nlatent)) {
    carriers <- effects[effects$type %in% "sd" & effects$row > spec$nlatent, , drop = FALSE]
    carriers <- carriers[order(carriers$row), , drop = FALSE]
    if (nrow(carriers) == count && !anyNA(carriers$param)) {
      return(as.character(carriers$param))
    }
  }
  paste0("indvar", seq_len(count))
}

# Subject matrices, reshaped from the engine's flat layout into the `subj_*`
# arrays ctExtract() reports for Stan fits: [iteration, subject, row, col].
.ctBackendSubjectMatrices <- function(spec, flat) {
  if (is.null(flat)) return(NULL)
  layout <- .ctBackendSummaryLayout(spec)
  out <- lapply(seq_along(layout$matrix), function(index) {
    count <- layout$nrow[index] * layout$ncol[index]
    block <- flat[, layout$offset[index] + seq_len(count), , drop = FALSE]
    value <- aperm(array(block, dim = c(dim(flat)[1], layout$nrow[index],
      layout$ncol[index], dim(flat)[3])), c(1L, 4L, 2L, 3L))
    .ctBackendTrimAugmented(value, layout$matrix[index], spec, margin = c(3L, 4L))
  })
  names(out) <- paste0("subj_", layout$matrix)
  out
}


# Posterior-predictive data generation ----------------------------------------
#
# The engines draw each row's observation from its own prior predictive as the
# filter reaches it, so the generated data is a draw from the model rather than
# a sequence of independent one-step predictions. The standard normals are drawn
# *here*, with R's RNG, so that `set.seed()` means what a user expects and the
# two engines produce identical data for the same seed.
#
# One deviate per cell covers every manifest type but the count. A count given
# its predictor is a Poisson, and given a dispersion a Poisson mixed over the
# predictor, so the draw needs a second, independent source of randomness rather
# than a second reading of the first -- reusing the one deviate for both stages
# ties them together and comes out over-dispersed. The julia engine draws that
# part from a stream of its own, seeded by the integer below, which is drawn
# here with R's RNG for exactly the reason the normals are: `set.seed()` has to
# govern the whole dataset, not most of it.
#
# `.Machine$integer.max` rather than anything wider because the engine splits
# the seed into two 32-bit halves and a subject index, matching the particle
# filter's streams.

.ctBackendGenerate <- function(fit, raw, base, effects = NULL) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  raw <- .ctJuliaNumericVector(as.numeric(raw))
  # A draw of the random effects, when the caller has one.
  #
  # Without it `ctsem_generate(laplace, ...)` puts each subject at its
  # *conditional mode*, solved from that subject's data. For a fit that
  # integrated the effects out that is the right quantity. For a fit that
  # sampled them it is not: the draws exist, and using modes conditions the
  # predictive on a point estimate of each effect rather than integrating over
  # its posterior, which understates between-subject uncertainty. Measured on a
  # 12-subject sampled fit, the generated subject means carried no relationship
  # to the effect draws they were supposed to come from -- mean |correlation|
  # 0.123 against a shuffled null of 0.087 -- while tracking the *mean* effect
  # across subjects at 0.965, which is the signature of modes.
  seed <- sample.int(.Machine$integer.max, 1L)
  if (!is.null(effects)) {
    persubject <- .ctBackendJuliaValue(module$ctsem_laplace_subject_values(
      objective, raw, .ctJuliaNumericVector(as.numeric(effects))))
    return(.ctBackendJuliaValue(module$ctsem_generate(objective, raw,
      JuliaConnectoR::juliaPut(base), seed = seed,
      subject_values = JuliaConnectoR::juliaPut(persubject))))
  }
  .ctBackendJuliaValue(module$ctsem_generate(objective, raw,
    JuliaConnectoR::juliaPut(base), seed = seed))
}


# The state-explicit route -----------------------------------------------------
#
# `intoverstates=FALSE`. Rather than drawing each row from the filter's
# one-step-ahead predictive and conditioning on the draw, this samples the
# latent trajectory from the process itself and then each observation from its
# conditional distribution given the state at its row.
#
# For a linear Gaussian model the two agree in distribution, because the
# filter's predictive is then exact. For a non-Gaussian one they do not: the
# filter's categorical update is an assumed-density projection that moves the
# state it conditions on, and with an unbounded indicator -- a count -- an
# improbable draw moves it far enough that the next row is drawn from a rate
# that has already run. Nothing of the kind can happen here, because no
# observation touches a state.
#
# Two blocks of standard normals, both drawn on the R side so `set.seed()`
# governs the result: one per latent innovation, and one per manifest cell.

.ctBackendStateDimension <- function(fit) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  as.integer(.ctBackendJuliaValue(
    module$ctsem_state_dimension(.ctJuliaObjective(fit))))
}

# Where each subject's innovations sit inside the vector
# `.ctBackendGenerateStates()` takes. Zero-based offsets, and the first
# `nlatent` entries of a subject's block are its initial state draw -- which is
# what lets a caller pin a person and redraw only its path.
.ctBackendStateLayout <- function(fit) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  layout <- .ctBackendJuliaValue(module$ctsem_state_layout(.ctJuliaObjective(fit)))
  lapply(layout, as.numeric)
}

.ctBackendGenerateStates <- function(fit, raw, z, base, effects = NULL) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  arguments <- list(.ctJuliaObjective(fit),
    .ctJuliaNumericVector(as.numeric(raw)),
    .ctJuliaNumericVector(as.numeric(z)), JuliaConnectoR::juliaPut(base),
    transition = .ctJuliaOr(spec$transition, "exponential"),
    seed = sample.int(.Machine$integer.max, 1L))
  # A draw of the random effects, so each subject's trajectory is drawn at that
  # subject's own parameters rather than at the population vector. Only the
  # random-effect objective takes it; the plain one has no effects to place.
  if (!is.null(effects)) {
    arguments$effects <- .ctJuliaNumericVector(as.numeric(effects))
  }
  .ctBackendJuliaValue(do.call(module$ctsem_generate_states, arguments))
}

# The joint density of states and data, and its gradient with respect to the
# two of them stacked. `z` is `.ctBackendStateDimension()` long.
.ctBackendJointDensity <- function(fit, raw, z, gradient = TRUE) {
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  .ctBackendJuliaValue(module$ctsem_joint_evaluate(.ctJuliaObjective(fit),
    .ctJuliaNumericVector(as.numeric(raw)),
    .ctJuliaNumericVector(as.numeric(z)), gradient = isTRUE(gradient)))
}

.ctBackendGenerateFromFit <- function(fit, nsamples = 200, fullposterior = FALSE,
  cores = 2, intoverstates = "fit") {
  spec <- .ctBackendSpec(fit)
  # `ctsem_generate` (and `ctsem_generate_states` below it) pass each
  # subject's `tipreds` straight to the extended Kalman filter without the
  # `TIMissingRecipe` substitution the adjoint/gradient path performs for a
  # sampled (missing) TI predictor value (see `_ctsem_subject_gradient_chunk!`
  # in adjoint.jl) -- so a subject carrying one dispatches to no matching
  # `_extended_kalman_filter_continuous!` method and fails with a raw Julia
  # `MethodError` naming an internal workspace type, not this model. Caught
  # here rather than left to surface that way.
  if (!is.null(spec$ti_missing) && nrow(spec$ti_missing)) {
    stop("ctGenerateFromFit()/ctPostPredict() are not available for a fit ",
      "with a sampled (missing) TI predictor value: the engine's generate ",
      "routine does not yet substitute the sampled value the way the ",
      "gradient does. See fit$model_spec$ti_missing.", call. = FALSE)
  }
  model <- .ctFitModelObject(fit)
  manifestNames <- model$manifestNames
  nmanifest <- length(manifestNames)
  nrows <- length(spec$times)

  # Random effect draws, paired with the parameter draws they were taken with.
  #
  # A fit that sampled its random effects carries one draw of every effect
  # alongside every draw of the parameters, and the two belong together: draw
  # `s` of the predictive is the model at draw `s` of everything. Pairing them
  # is what makes this integrate over the random effects. Without it the engine
  # falls back to conditional modes, which is a point estimate of each effect.
  #
  # `sample$effects` exists only when the fit was sampled with
  # `saveEffects=TRUE` -- the draws are nsubjects x neffects x chains x draws
  # numbers and the bridge moves about 1 MB/s, so they are summarised by
  # default. Where they are absent this says so and uses the modes, rather than
  # substituting them in silence.
  effectdraws <- fit$sample$effects
  sampledeffects <- !is.null(fit$sample) &&
    (!is.null(effectdraws) || length(fit$sample$effect_mean))

  if (isTRUE(fullposterior)) {
    posterior <- fit$estimate$rawposterior
    if (is.null(posterior)) {
      stop("fullposterior=TRUE needs posterior draws; run ctOptimUncertainty() first, ",
        "or use fullposterior=FALSE to generate from the point estimate.", call. = FALSE)
    }
    rows <- sample(seq_len(nrow(posterior)), nsamples, replace = nsamples > nrow(posterior))
    samples <- posterior[rows, , drop = FALSE]
    # The same rows, so parameters and effects come from one draw. They are
    # only alignable when the two have the same number of draws; a fit whose
    # `rawposterior` came from somewhere other than the sampler (a Hessian
    # draw, say) is not paired with these effects and must not be treated as
    # though it were.
    effects <- if (!is.null(effectdraws) &&
        nrow(effectdraws) == nrow(posterior)) effectdraws[rows, , drop = FALSE]
  } else {
    samples <- matrix(as.numeric(fit$estimate$raw), nrow = nsamples,
      ncol = length(fit$estimate$raw), byrow = TRUE)
    # A point estimate of the parameters takes the point estimate of the
    # effects, which is their posterior mean -- available whether or not the
    # draws themselves were kept.
    mean_effects <- if (length(fit$sample$effect_mean))
      as.numeric(fit$sample$effect_mean) else
        if (!is.null(effectdraws)) colMeans(effectdraws)
    effects <- if (!is.null(mean_effects))
      matrix(mean_effects, nrow = nsamples, ncol = length(mean_effects),
        byrow = TRUE)
  }

  if (sampledeffects && is.null(effects)) {
    message("This fit sampled its random effects, but the draws needed to pair ",
      "them with the parameter draws are not on it, so each subject is ",
      "generated at its conditional mode instead -- a point estimate of its ",
      "effect rather than a draw from its posterior. Re-sample with ",
      "saveEffects=TRUE to generate from the draws.")
  }

  generated <- array(NA_real_, dim = c(nsamples, nrows, nmanifest))
  llrow <- matrix(0, nsamples, nrows)
  # Which route generates the states, and it is a choice rather than a
  # consequence.
  #
  # `'fit'`, the default, mirrors the fit: a state-explicit fit gets a
  # state-explicit predictive, because generating through the filter at its
  # estimate would be a draw from a density that fit did not maximise, and it
  # would look entirely reasonable. It also keeps the identity that the
  # generated dataset's `llrow` is the likelihood reported while generating it.
  #
  # `FALSE` resamples the trajectory whatever the fit did -- drawn from the
  # process at that subject's own parameters, then each observation from its
  # conditional distribution given the state. That is a draw from the model
  # rather than from the filter's approximation of it, and the two differ for a
  # non-Gaussian indicator, where the measurement update is an assumed-density
  # projection, or for state-dependent dynamics. `TRUE` forces the filter
  # route.
  #
  # Either way the states are *regenerated*: nothing a fit sampled is carried
  # over. What comes from the fit is the individual differences.
  statepath <- if (identical(as.character(intoverstates)[1L], "fit"))
    isFALSE(fit$args$resolved$intoverstates) else
      !isTRUE(as.logical(intoverstates)[1L])
  nz <- if (statepath) .ctBackendStateDimension(fit) else 0L
  for (iteration in seq_len(nsamples)) {
    # Innovations first, then the observation deviates, matching the order
    # `.ctGenerateJulia` draws them in.
    z <- if (statepath) stats::rnorm(nz) else NULL
    base <- matrix(stats::rnorm(nmanifest * nrows), nmanifest, nrows)
    drawn <- if (statepath) {
      .ctBackendGenerateStates(fit, samples[iteration, ], z, base,
        effects = if (!is.null(effects)) effects[iteration, ])
    } else .ctBackendGenerate(fit, samples[iteration, ], base,
      effects = if (!is.null(effects)) effects[iteration, ])
    generated[iteration, , ] <- t(matrix(as.numeric(drawn$Y), nmanifest, nrows))
    llrow[iteration, ] <- as.numeric(drawn$llrow)
  }
  dimnames(generated) <- list(sample = seq_len(nsamples), row = seq_len(nrows),
    manifestNames)
  llrow[llrow == 0] <- NA
  fit$generated <- list(Y = generated, llrow = llrow)
  fit
}


# Accessors the posterior-predictive machinery needs, for whichever backend ----
#
# ctPostPredData() reaches into `standata` for the observed data, the row-to-
# subject map and the fitted row likelihoods. These give the same four things
# from either a ctStanFit or a backend fit, so that function has one body.

.ctFitObservedY <- function(fit) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata$Y)) {
    observed <- fit$standata$Y
    observed[observed == 99999] <- NA
    colnames(observed) <- .ctFitModelObject(fit)$manifestNames
    return(observed)
  }
  spec <- .ctBackendSpec(fit)
  observed <- t(spec$manifest_data)
  observed[!is.finite(observed)] <- NA
  colnames(observed) <- .ctFitModelObject(fit)$manifestNames
  observed
}

.ctFitRowSubject <- function(fit) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata$subject)) return(as.integer(fit$standata$subject))
  spec <- .ctBackendSpec(fit)
  starts <- spec$subject_starts
  rep(seq_along(starts), diff(c(starts, length(spec$times) + 1L)))
}

.ctFitRowTime <- function(fit) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata$time)) return(as.numeric(fit$standata$time))
  as.numeric(.ctBackendSpec(fit)$times)
}

# Each row's log likelihood at the fitted estimate -- the quantity a posterior
# predictive check compares the generated ones against.
.ctFitObservedRowLoglik <- function(fit) {
  if (!is.null(fit$stanfit$transformedparsfull$llrow)) {
    return(as.numeric(fit$stanfit$transformedparsfull$llrow[1, ]))
  }
  spec <- .ctBackendAsModel(.ctBackendSpec(fit))
  as.numeric(.ctBackendKalmanRaw(spec, fit$estimate$raw, subjectmatrices = FALSE,
    fields = "llrow")$llrow)
}

# A copy of the fit whose observed data has been replaced, so that the residual
# branch of ctPostPredData() can filter a generated dataset.
.ctFitReplaceY <- function(fit, Y) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata$Y)) {
    fit$standata$Y <- matrix(Y, ncol = ncol(fit$standata$Y))
    return(fit)
  }
  data <- .ctBackendSpec(fit)$data
  data[, .ctFitModelObject(fit)$manifestNames] <- Y
  .ctFitReplaceData(fit, data)
}

# The observed data as a long data frame in its original structure. A ctStanFit
# has to reconstruct it from `standata`; a julia fit now carries a `standata`
# of its own too (see R/ctFit.R), so this reconstructs it the same way for
# both -- unlike the accessors above, nothing here reads another accessor's
# output by row position, so there is no ordering hazard in preferring
# `standata` whichever backend produced it.
#
# `standatatolong()`'s id column is `standata$subject`, the ascending 1:N
# `makeNumericIDs()` (R/ctsemUtils.R) assigned -- never the user's own ids,
# which only `standata$idmap` (original/new, built right before that
# remapping in .ctPrepareData(), R/ctData.R) still knows. Left unmapped, a stan
# fit returned the internal index and a julia fit returned the user's real
# ids (from the verbatim frame this used to fall back to), so the two
# backends silently disagreed about what a subject is called -- a
# `merge(userdata, ..., by='id')` downstream would mismatch every row on one
# backend and not the other. Mapped back here, both return the ids the user
# supplied.
.ctFitLongData <- function(fit) {
  if (!is.null(fit$standata)) {
    long <- standatatolong(standata = fit$standata, ctm = .ctFitModelObject(fit),
      origstructure = TRUE)
    idname <- .ctFitModelObject(fit)$subjectIDname
    idmap <- fit$standata$idmap
    if (!is.null(idmap) && idname %in% names(long)) {
      long[[idname]] <- idmap$original[match(long[[idname]], idmap$new)]
    }
    return(long)
  }
  .ctBackendSpec(fit)$data
}

# TI predictor values per subject, for whichever backend. ctPredictTIP() reads
# these to pick the covariate values it predicts at.
.ctFitTIpredData <- function(fit) {
  model <- .ctFitModelObject(fit)
  values <- if (!.ctFitIsJulia(fit) && !is.null(fit$standata$tipredsdata)) {
    as.matrix(fit$standata$tipredsdata)
  } else as.matrix(.ctBackendSpec(fit)$tipred_data)
  if (ncol(values) == length(model$TIpredNames)) colnames(values) <- model$TIpredNames
  values
}

# What varies between units in this fit, however it is represented ------------
#
# ctsem has two representations of individual differences and they share no
# field. `intoverpop='augmented'` carries each varying parameter as a latent
# state and describes the set in `spec$random_effects`; `'laplace'` keeps them
# as separate coordinates, describes them in `spec$laplace$levels`, and leaves
# `spec$random_effects` empty. A stan fit is always the augmented kind and
# describes itself through `standata` instead.
#
# A consumer that knows one representation and not the other does not fail on
# the other. It reports a model with no individual differences in it -- which is
# a legitimate model, so the answer looks like an answer. That is how
# ctVarianceDecomposition() came to return a between person variance of exactly
# zero for every Laplace fit, which is every multilevel model, while the
# information sat unread in `spec$laplace$levels`.
#
# So the question "what varies, between which units" is answered once, here,
# and the three representations differ only inside this function. Returns a list
# of levels, innermost first, each with
#
#   name    the id column that indexes the level
#   units   which unit of this level each subject belongs to, one entry per
#           subject -- `seq_len(nsubjects)` at the subject level
#   nunits  how many units the level has
#   params  the parameters that vary at it
#
# and an empty list when the model declares no individual differences at all.
# That emptiness is a statement about the model rather than about what this
# function could find, which is the distinction the failure above turned on:
# `.ctFitHasRandomEffects()` below is what a caller should ask before trusting a
# between-unit number it computed some other way.
# The julia layer of it. `spec` is a prepared specification, not a fit, which is
# what the call sites inside the backend have to hand.
#
# `labels` is left NULL at the subject level here: naming a subject needs the
# id map, which belongs to the fit. `.ctFitRandomEffectLevels()` fills it.
.ctSpecRandomEffectLevels <- function(spec, idname = NULL) {
  if (!is.null(spec$laplace) && length(spec$laplace$levels)) {
    return(lapply(spec$laplace$levels, function(level) {
      # `group` is one entry per subject, saying which unit of this level that
      # subject is in; `labels` is one entry per unit, the identifier the user
      # wrote. Both are built where the hierarchy is read (R/ctJuliaBackend.R),
      # so nothing downstream has to reconstruct either from the data -- which
      # is what the two call sites that used to do it got wrong in different
      # ways.
      list(name = as.character(level$name)[1L],
        units = as.integer(level$group),
        nunits = as.integer(level$ngroups)[1L],
        params = as.character(level$param),
        labels = if (is.null(level$labels)) NULL else as.character(level$labels))
    }))
  }
  effects <- spec$random_effects
  if (is.null(effects) || !length(effects) || !NROW(effects)) return(list())
  params <- unique(as.character(effects$param[effects$type %in% "sd"]))
  params <- params[!is.na(params)]
  if (!length(params)) return(list())
  n <- length(spec$subject_starts)
  list(list(name = idname, units = seq_len(n), nunits = n, params = params,
    labels = NULL))
}

.ctFitRandomEffectLevels <- function(fit) {
  idname <- .ctFitModelObject(fit)$subjectIDname

  # Selected by what the object *is*, not by the absence of a julia field.
  # `.ctFitIsJulia()` asks for `$model_spec`, which a prepared model from
  # `ctFit(fit = FALSE)` does not have -- it *is* the specification -- so that
  # test sent every prepared model down the stan branch, where `getparnames()`
  # fails and the answer came back as "no random effects". Which is this
  # function's whole failure mode, reintroduced inside it.
  levels <- if (inherits(fit, 'ctStanFit')) {
    varying <- tryCatch(getparnames(fit, subjvariationonly = TRUE),
      error = function(e) character())
    if (!length(varying)) list() else {
      n <- as.integer(fit$standata$nsubjects)
      list(list(name = idname, units = seq_len(n), nunits = n,
        params = as.character(varying), labels = NULL))
    }
  } else {
    spec <- .ctBackendSpec(fit)
    # An empty answer has to mean "this model declares no individual
    # differences", so anything this cannot read says so instead of agreeing.
    if (is.null(spec$parameter_table)) {
      stop('Cannot read the random effect structure of this object: it is ',
        'neither a stan fit nor a prepared julia specification.', call. = FALSE)
    }
    .ctSpecRandomEffectLevels(spec, idname = idname)
  }

  # The innermost level's units are the subjects, so its labels are their ids.
  if (length(levels) && is.null(levels[[1L]]$labels)) {
    map <- try(.ctFitIdMap(fit), silent = TRUE)
    if (!inherits(map, "try-error") && nrow(map) == levels[[1L]]$nunits) {
      levels[[1L]]$labels <- as.character(map$original)
    }
  }
  levels
}

# Does this fit declare individual differences at all?
#
# The question worth asking before reporting a between-unit quantity computed
# some other way: an empty answer here is a statement about the model, where
# "I found no carrier states" is a statement about one representation.
.ctFitHasRandomEffects <- function(fit) length(.ctFitRandomEffectLevels(fit)) > 0L

# A copy of the fit re-prepared against a different long data frame. This is
# what lets ctPredictTIP() build its covariate grid: it constructs a dataset of
# pseudo-subjects, one per covariate value, and asks the fitted model to predict
# for them.
.ctFitReplaceData <- function(fit, datalong) {
  if (!.ctFitIsJulia(fit) && !is.null(fit$standata)) {
    fit$standata <- suppressMessages(.ctPrepareData(fit$ctstanmodel, datalong, optimize = TRUE))
    return(fit)
  }
  spec <- .ctBackendSpec(fit)
  prepared <- .ctJuliaPrepare(datalong, .ctFitModelObject(fit),
    project = spec$project, intoverpop = .ctBackendIntOverPop(spec),
    laplacecontrol = spec$laplace$inner)
  # Carry across the two things that are properties of the *fit* rather than of
  # the data, and that re-preparation would otherwise silently drop: the prior
  # specification (a function of the model, and `priors` defaults to FALSE
  # here) and the integration step. Cross-validation re-optimises against a
  # re-prepared specification, and a refit that quietly lost its priors or its
  # step size would not be the same model. A step chosen by
  # `nsubsteps = 'auto'` is a count per fitted row, and goes only with those
  # rows: ctPredictTIP()'s pseudo-subjects are none of them.
  prepared$priors <- spec$priors
  fit$model_spec <- .ctJuliaCarrySubsteps(prepared, spec, .ctFitRawEstimate(fit))
  fit
}
