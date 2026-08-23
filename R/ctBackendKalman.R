# Prediction and Kalman output for backend='julia' and backend='cpp' ---------
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

# The original-id to internal-index mapping, for the backends that do not carry
# a `standata`. ctPredict() speaks in the user's own subject ids and the filter
# speaks in positions, so something has to hold the correspondence.
.ctFitIdMap <- function(fit) {
  if (!is.null(fit$standata$idmap)) return(fit$standata$idmap)
  spec <- .ctBackendSpec(fit)
  ids <- unique(spec$data[[.ctFitModelObject(fit)$subjectIDname]])
  data.frame(original = ids, new = seq_along(ids), stringsAsFactors = FALSE)
}

.ctBackendKalmanRaw <- function(fit, raw, subjectmatrices = TRUE) {
  raw <- as.numeric(raw)
  if (identical(.ctBackendEngineKind(fit), "cpp")) {
    return(.ctsemCppKalman(.ctCppObjective(fit), raw, isTRUE(subjectmatrices)))
  }
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  result <- .ctBackendJuliaValue(module$ctsem_kalman(.ctJuliaObjective(fit),
    .ctJuliaNumericVector(raw), subject_matrices = isTRUE(subjectmatrices)))
  result$subject <- as.integer(result$subject)
  result
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
    return(.ctBackendAsModel(spec, .ctBackendEngineKind(fit)))
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

  kind <- .ctBackendEngineKind(fit)
  prepared <- if (identical(kind, "cpp")) {
    .ctCppPrepare(dat, model)
  } else {
    .ctJuliaPrepare(dat, model, project = spec$project)
  }
  prepared <- .ctBackendAsModel(prepared, kind)
  if (withhold) attr(prepared, "reportManifest") <- reported
  prepared
}

.ctBackendAsModel <- function(spec, kind) {
  structure(spec,
    class = c(if (identical(kind, "cpp")) "ctCppModel" else "ctJuliaModel", "ctFitModel"))
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

#' Kalman filter and smoother estimates from a julia or cpp backend fit
#'
#' Prior, filtered and smoothed estimates of the latent states and the
#' observations, for every row of data, from the same forward pass the engine
#' uses for the likelihood.
#'
#' @param fit A \code{ctJuliaFit} or \code{ctCppFit}.
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
#' @export
ctBackendKalman <- function(fit, subjects = "all", timestep = "asdata",
  maxtime = "asdata", removeObs = FALSE, pointest = TRUE, nsamples = NA,
  collapsefunc = NA, standardisederrors = FALSE, subjectpars = FALSE,
  indvarstates = FALSE, ...) {

  spec <- .ctBackendKalmanSpec(fit, subjects = subjects, timestep = timestep,
    maxtime = maxtime, removeObs = removeObs)
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
      subjectmatrices = isTRUE(subjectpars))
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
  if (!is.null(effects) && nrow(effects) >= count && !is.null(effects$param)) {
    return(as.character(effects$param)[seq_len(count)])
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
