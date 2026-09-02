# Transformed-parameter summaries for backend='julia' ------------------------
#
# ctsem's summary and plot functions all read the same thing: samples of the
# *model matrices* implied by the raw parameter vector, as `pop_DRIFT`,
# `pop_DIFFUSIONcov` and so on, each an [iteration, row, column] array. Given
# those, `ctSummaryMatrices()`, `summary()` and the plotting helpers work
# unchanged.
#
# So the architecture here is one primitive and a thin stack on top of it:
#
#   .ctBackendParMatricesFlat()  ask the engine to materialize every model
#                                matrix for a matrix of raw vectors
#   ctBackendParMatrices()       the same, reshaped and named, for one vector
#   ctExtract()                  the same over the posterior, as pop_* arrays
#   ctSummaryMatrices()          collapse those arrays -- the *identical* code
#                                path Stan fits use, factored out below
#   summary()                    fixed effects plus the system matrices
#
# The primitive runs inside the engine rather than re-deriving the transforms in
# R. The engine already materializes every model matrix from the raw vector on
# its way to a log likelihood, so asking it for that same materialization is the
# only way to guarantee a summary reports what the likelihood actually used. A
# second, R-side implementation of ctsem's transforms is exactly the kind of
# duplication that has let this package's backends drift apart before.
#
# Two properties of these fits shape what can honestly be reported:
#
#   * They carry `estimate$rawposterior`, because `ctFit()` finishes with
#     `ctOptimUncertainty()` as `stanoptimis()` does, and every interval below
#     comes from pushing those draws through the transforms -- which is what
#     Stan's own optimized-and-sampled path does too. A fit built with
#     `optimcontrol$estonly=TRUE` has one "sample" instead, and the interval
#     columns are omitted rather than filled with a zero-width interval that
#     would read as certainty.
#
#   * State-dependent cells are functions of the latent state, so no single
#     number describes them. They are evaluated at a state (T0MEANS by default)
#     and reported as conditional on it, not as if they were constants.

.ctBackendSpec <- function(fit) {
  if (!is.null(fit$model_spec)) fit$model_spec else fit
}

# Model metadata and sample count, for the summary and plot functions that are
# shared across backends. A ctStanFit keeps the model in `$ctstanmodel` and its
# samples in `$stanfit`; a julia fit keeps them in `$model` and `$estimate`.
# These two accessors are what let those functions stop reaching into either
# layout directly.
.ctFitModelObject <- function(fit) {
  if (!is.null(fit$ctstanmodel)) return(fit$ctstanmodel)
  # The specification's own model, ahead of the fit's, because that is the
  # object the parameter table was actually derived from. It matters when
  # something re-prepares the data (ctPredict interpolating a time grid, say):
  # ctFit hands the backends a model that has already been through
  # ctStanModelIntOverPop, and re-deriving the augmentation from a model that
  # has not been produces an algebraically equivalent but differently written
  # parameter table -- with a different mapping from the raw vector.
  if (!is.null(fit$model_spec$model)) return(fit$model_spec$model)
  if (!is.null(fit$model)) return(fit$model)
  stop("The fit does not carry the model it was built from.", call. = FALSE)
}

.ctFitNsubjects <- function(fit) {
  if (!is.null(fit$standata$subject)) return(length(unique(fit$standata$subject)))
  length(.ctBackendSpec(fit)$subject_starts)
}

.ctFitNsamples <- function(fit) {
  if (!is.null(fit$stanfit$transformedpars$pop_DRIFT)) {
    return(dim(fit$stanfit$transformedpars$pop_DRIFT)[1])
  }
  posterior <- fit$estimate$rawposterior
  if (!is.null(posterior) && length(dim(posterior)) == 2L) return(nrow(posterior))
  1L
}

.ctBackendModel <- .ctFitModelObject

# Names and dimensions of everything the engine can materialize, plus which
# cells are state dependent. A property of the model, not of the parameter
# values -- so it is computed once per objective and cached, rather than
# re-asked every time something wants to know where a matrix lives.
#
# It is cheap engine-side and was described as cheap, but it is one to two
# JuliaConnectoR round trips, and on this backend a round trip is 40-80 ms
# whatever it carries. `.ctBackendConstrain` alone asks for it once per call:
# measured on dev1 at 0.29 s on a 4-latent model and 0.33 s on a 2-latent one,
# against a whole constrain step of 2.6 s and 2.1 s.
#
# The key is the objective's, so a model whose data or parameter table
# changed gets a fresh layout for the same reason it gets a fresh objective.
.ctBackendSummaryLayout <- function(fit) {
  spec <- .ctBackendSpec(fit)
  key <- .ctJuliaObjectiveKey(spec)
  if (exists(key, envir = .ct_julia_cache$layouts, inherits = FALSE)) {
    return(get(key, envir = .ct_julia_cache$layouts, inherits = FALSE))
  }
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  raw <- .ctBackendJuliaValue(module$ctsem_parameter_layout(objective))
  # Fetched separately, and only when there are any: JuliaConnectoR hangs
  # marshalling a zero-length vector, and a model with no state-dependent cells
  # is the common case, not an edge one.
  statedep <- if (as.integer(raw$n_statedep)[1L] > 0L) {
    cells <- .ctBackendJuliaValue(module$ctsem_state_dependent_cells(objective))
    data.frame(matrix = as.character(cells$matrix), row = as.integer(cells$row),
      col = as.integer(cells$col), stringsAsFactors = FALSE)
  } else {
    data.frame(matrix = character(), row = integer(), col = integer(),
      stringsAsFactors = FALSE)
  }
  layout <- list(matrix = as.character(raw$matrix), nrow = as.integer(raw$nrow),
    ncol = as.integer(raw$ncol), offset = as.integer(raw$offset),
    size = as.integer(raw$size)[1L], nlatent = as.integer(raw$nlatent)[1L],
    nmanifest = as.integer(raw$nmanifest)[1L], statedep = statedep)
  assign(key, layout, envir = .ct_julia_cache$layouts)
  layout
}

# `raw` is npar x nsamples; the result is (flat layout) x nsamples. The whole
# posterior travels in one call: the Julia bridge marshals a numeric array in
# one transfer but a list element by element, and per-sample calls were the
# single largest avoidable cost measured in this backend.
.ctBackendParMatricesFlat <- function(fit, raw, tipreds = NULL, state = NULL,
  time = 0, dt = 0, rows = NULL) {
  raw <- if (is.matrix(raw)) raw else matrix(as.numeric(raw), ncol = 1L)
  storage.mode(raw) <- "double"
  spec <- .ctBackendSpec(fit)
  module <- .ctJuliaModule(spec$project)
  arguments <- list(.ctJuliaObjective(fit), JuliaConnectoR::juliaPut(raw),
    time = as.numeric(time)[1L], dt = as.numeric(dt)[1L])
  # Optional keywords are passed only when non-empty: JuliaConnectoR hangs
  # marshalling a zero-length vector.
  if (length(tipreds)) arguments$tipreds <- .ctJuliaVector(as.numeric(tipreds))
  if (length(state)) arguments$state <- .ctJuliaVector(as.numeric(state))
  # `rows` selects flat positions engine-side. The materialization is under
  # 0.02 s either way; what this saves is the bridge, which moves about 1 MB/s
  # and does not care that 98% of what it is carrying will be dropped on
  # arrival.
  if (length(rows)) arguments$rows <- .ctJuliaVector(as.integer(rows))
  result <- .ctBackendJuliaValue(do.call(module$ctsem_parameter_matrices, arguments))
  matrix(as.numeric(result), nrow = nrow(result), ncol = ncol(result))
}

# JuliaConnectoR materialises plain arrays and named tuples on the R side
# already and only hands back a proxy for what it cannot; `juliaGet` errors on
# the former, so ask for it only when there is something to fetch.
.ctBackendJuliaValue <- function(x) {
  if (inherits(x, "JuliaProxy")) JuliaConnectoR::juliaGet(x) else x
}

# Which matrices carry a latent dimension, and on which side.
#
# An `intoverpop` model augments the latent state with a carrier state per
# individually-varying parameter, and the engines work throughout in that
# augmented space. Stan does too, but what it *reports* is trimmed: only the
# T0 matrices keep the carrier rows, because that is where the random effects
# live; everything else is reported over the real latent processes. Matching
# that here rather than in the collapse is deliberate -- it makes `pop_DRIFT`
# mean the same thing whichever backend produced it, which is the entire point
# of sharing the summary code below.
.ctBackendLatentRows <- c("DRIFT", "DIFFUSION", "CINT", "TDPREDEFFECT", "JAx", "Jtd",
  "DIFFUSIONcov", "asymCINT", "asymDIFFUSIONcov")
.ctBackendLatentCols <- c("DRIFT", "DIFFUSION", "LAMBDA", "JAx", "Jtd", "Jy",
  "DIFFUSIONcov", "asymDIFFUSIONcov")

.ctBackendTrimAugmented <- function(value, name, spec, margin = c(1L, 2L)) {
  nlatent <- spec$nlatent
  augmented <- spec$nlatent_augmented
  if (is.null(nlatent) || is.null(augmented) || identical(nlatent, augmented)) return(value)
  index <- lapply(dim(value), seq_len)
  if (name %in% .ctBackendLatentRows && dim(value)[margin[1L]] == augmented) {
    index[[margin[1L]]] <- seq_len(nlatent)
  }
  if (name %in% .ctBackendLatentCols && dim(value)[margin[2L]] == augmented) {
    index[[margin[2L]]] <- seq_len(nlatent)
  }
  do.call(`[`, c(list(value), index, list(drop = FALSE)))
}

.ctBackendReshape <- function(flat, layout, index) {
  count <- layout$nrow[index] * layout$ncol[index]
  block <- flat[layout$offset[index] + seq_len(count), , drop = FALSE]
  # (rows*cols) x iterations, column-major within a sample, to [iter, row, col].
  aperm(array(block, dim = c(layout$nrow[index], layout$ncol[index], ncol(block))),
    c(3L, 1L, 2L))
}

#' Model-implied parameter matrices from a julia backend fit
#'
#' Materialise every model matrix -- and the covariance and asymptotic matrices
#' derived from them -- from a raw (unconstrained) parameter vector, using the
#' same engine code the likelihood uses.
#'
#' @param fit A \code{ctJuliaFit}, or a prepared model from
#'   \code{ctFit(..., fit=FALSE)}.
#' @param raw Raw parameter vector. Defaults to the fitted estimate.
#' @param tipreds Time-independent predictor values for the subject to
#'   materialise. Defaults to all zero, i.e. the population values.
#' @param state Latent state at which to evaluate state-dependent cells.
#'   Defaults to \code{T0MEANS}.
#' @param time,dt Time and time interval passed to state-dependent expressions
#'   that use them.
#' @param trim Report matrices over the real latent processes (the default), as
#'   ctsem's summaries do, rather than over the augmented state an
#'   \code{intoverpop} model filters in. \code{trim=FALSE} returns the matrices
#'   the filter actually uses, carrier states included.
#' @return A named list of matrices, with an attribute \code{stateDependent}
#'   giving the cells whose values are conditional on \code{state}.
#' @details For a linear model no cell depends on the state and the result is
#'   exact. For a nonlinear one, the state-dependent cells are evaluated at
#'   \code{state} and are conditional on it; \code{attr(x, 'stateDependent')}
#'   names them.
#' @examples
#' \donttest{
#' # ctBackendParMatrices(fit)$DRIFT
#' }
#' @export
ctBackendParMatrices <- function(fit, raw = NULL, tipreds = NULL, state = NULL,
  time = 0, dt = 0, trim = TRUE) {
  if (is.null(raw)) {
    raw <- fit$estimate$raw
    if (is.null(raw)) stop("raw must be supplied for a fit without an estimate.", call. = FALSE)
  }
  layout <- .ctBackendSummaryLayout(fit)
  spec <- .ctBackendSpec(fit)
  flat <- .ctBackendParMatricesFlat(fit, raw, tipreds = tipreds, state = state,
    time = time, dt = dt)
  out <- lapply(seq_along(layout$matrix), function(index) {
    value <- matrix(
      flat[layout$offset[index] + seq_len(layout$nrow[index] * layout$ncol[index]), 1L],
      nrow = layout$nrow[index], ncol = layout$ncol[index])
    if (isTRUE(trim)) {
      value <- .ctBackendTrimAugmented(value, layout$matrix[index], spec)
    }
    value
  })
  names(out) <- layout$matrix
  out <- .ctBackendNameMatrices(out, .ctBackendModel(fit))
  attr(out, "stateDependent") <- layout$statedep
  out <- .ctContextAttach(out, fit)
  # Only when the caller did not choose the point themselves.
  if (is.null(state)) .ctContextMessage(fit, .ctContextPopLabel)
  out
}

# Samples of the raw parameters: the posterior draws from ctOptimUncertainty()
# when present, otherwise the point estimate as a single row.
.ctBackendRawSamples <- function(fit) {
  posterior <- fit$estimate$rawposterior
  if (!is.null(posterior) && length(dim(posterior)) == 2L && nrow(posterior) > 0L) {
    return(matrix(as.numeric(posterior), nrow = nrow(posterior)))
  }
  matrix(as.numeric(fit$estimate$raw), nrow = 1L)
}

# The `pop_*` arrays every downstream summary reads, in Stan's shape.
#
# Split from the engine call because a summary needs the same arrays collapsed
# five different ways, and materializing them is the expensive half: ~2 s per
# thousand draws, against ~0.1 s for a collapse. `summary()` therefore calls the
# engine once and reshapes five times.
.ctBackendPopArraysFromFlat <- function(flat, layout, spec) {
  out <- lapply(seq_along(layout$matrix), function(index) {
    .ctBackendTrimAugmented(.ctBackendReshape(flat, layout, index), layout$matrix[index],
      spec, margin = c(2L, 3L))
  })
  names(out) <- paste0("pop_", layout$matrix)
  out
}

.ctBackendPopArrays <- function(fit, samples = NULL, tipreds = NULL, state = NULL,
  time = 0, dt = 0) {
  # The fit's cached constrain step covers the default evaluation point only:
  # asking for another `state`, `tipreds`, `time` or `dt` is asking for
  # different matrices, so those go to the engine.
  if (is.null(tipreds) && is.null(state) && identical(time, 0) && identical(dt, 0)) {
    constrained <- .ctBackendConstrained(fit, samples)
    return(.ctBackendPopArraysFromFlat(constrained$flat, constrained$layout,
      .ctBackendSpec(fit)))
  }
  layout <- .ctBackendSummaryLayout(fit)
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  flat <- .ctBackendParMatricesFlat(fit, t(samples), tipreds = tipreds, state = state,
    time = time, dt = dt)
  .ctBackendPopArraysFromFlat(flat, layout, .ctBackendSpec(fit))
}

.ctBackendExtract <- function(object, subjectMatrices = FALSE, nsamples = "all",
  subjects = "all", ...) {
  samples <- .ctBackendRawSamples(object)
  if (!identical(nsamples, "all")) {
    wanted <- min(nrow(samples), as.integer(nsamples)[1L])
    samples <- samples[round(seq(1, nrow(samples), length.out = wanted)), , drop = FALSE]
  }
  # Subject matrices need the filter, not just the transforms: an individually
  # varying parameter is an augmented latent state, so a subject's value for it
  # is only known once that subject's data has been filtered and smoothed.
  #
  # A draw from a normal approximation can land somewhere the filter cannot go
  # -- a prior covariance the smoother's solve finds singular, say -- and one
  # such draw must not take the whole extract with it. `stan_constrainsamples()`
  # has always dropped inadmissable samples and reported the proportion; this
  # does the same, and drops them from the rest of the extract too, so every
  # array it returns is over the same surviving draws.
  subject <- NULL
  if (isTRUE(subjectMatrices)) {
    spec <- .ctBackendKalmanSpec(object, subjects = subjects)
    computed <- lapply(seq_len(nrow(samples)), function(iteration) {
      # Only `subject_matrices` is read below, which `fields` does not gate --
      # it is assembled separately from the filter trace regardless of which
      # named arrays are requested. The bridge cannot marshal a zero-length
      # vector, so `fields` still needs one name; `subject_loglik` is the
      # smallest of the eight.
      scores <- try(.ctBackendKalmanRaw(spec, samples[iteration, ],
        subjectmatrices = TRUE, fields = "subject_loglik"), silent = TRUE)
      if (inherits(scores, "try-error")) NULL else scores$subject_matrices
    })
    admissable <- !vapply(computed, is.null, logical(1L))
    if (!any(admissable)) stop("No admissable samples!?", call. = FALSE)
    if (any(!admissable)) {
      message(round(mean(!admissable) * 100, 1), "% of samples inadmissable")
      samples <- samples[admissable, , drop = FALSE]
      computed <- computed[admissable]
    }
    flat <- array(0, dim = c(length(computed), dim(computed[[1L]])))
    for (iteration in seq_along(computed)) flat[iteration, , ] <- computed[[iteration]]
    subject <- .ctBackendSubjectMatrices(spec, flat)
  }

  arrays <- .ctBackendPopArrays(object, samples = samples, ...)
  popmeans <- .ctBackendPopMeanSamples(object, samples = samples)
  c(list(rawpars = samples, popmeans = popmeans$values,
    loglik = object$estimate$loglik, gradient = object$estimate$gradient,
    subject_loglik = object$estimate$subject_loglik), arrays, subject)
}

.ctBackendNameMatrices <- function(out, model) {
  latent <- model$latentNames
  manifest <- model$manifestNames
  tdpred <- model$TDpredNames
  setnames <- function(name, rows, cols) {
    if (is.null(out[[name]]) || !length(dim(out[[name]]))) return(invisible(NULL))
    if (length(rows) == nrow(out[[name]]) && length(cols) == ncol(out[[name]])) {
      dimnames(out[[name]]) <<- list(rows, cols)
    } else if (length(rows) == nrow(out[[name]])) {
      rownames(out[[name]]) <<- rows
    }
    invisible(NULL)
  }
  for (name in c("DRIFT", "DIFFUSION", "DIFFUSIONcov", "T0VAR", "T0cov",
    "asymDIFFUSIONcov", "JAx")) setnames(name, latent, latent)
  for (name in c("MANIFESTVAR", "MANIFESTcov", "Jy")) setnames(name, manifest, manifest)
  setnames("LAMBDA", manifest, latent)
  if (length(tdpred)) setnames("TDPREDEFFECT", latent, tdpred)
  for (name in c("T0MEANS", "CINT", "asymCINT")) setnames(name, latent, "")
  setnames("MANIFESTMEANS", manifest, "")
  out
}


# Shared summary-matrix collapse ---------------------------------------------
#
# Factored out of `ctSummaryMatrices.ctStanFit()` unchanged so that the Stan
# and Julia backends collapse identical inputs with identical code. Anything
# that changes about how ctsem summarises system matrices changes here, once.

.ctSummaryMatricesFromArrays <- function(e, continuoustime, latentNames, manifestNames,
  TDpredNames, calcfunc = quantile, calcfuncargs = list(probs = 0.5), timeinterval = 1) {

  mats <- ctStanMatricesList()
  mats <- c(names(mats$base), names(mats$asymptotic), names(mats$extra))
  if (isTRUE(continuoustime)) {
    d <- list(DRIFT = e$pop_DRIFT)
    dd <- ctDiscreteParsDrift(d, timeinterval, observational = FALSE, standardise = FALSE,
      cov = FALSE, quiet = TRUE)
    e$pop_dtDRIFT <- array(dd, dim = dim(dd)[-2:-3])
    mats <- c(mats, "dtDRIFT")
  }

  out <- list()
  for (matname in mats) {
    try({
      calcfuncargs$collapsemargin <- 1
      calcfuncargs$collapsefunc <- calcfunc
      calcfuncargs$na.rm <- TRUE
      calcfuncargs$inarray <- e[[paste0("pop_", matname)]]
      out[[matname]] <- array(do.call(ctCollapse, calcfuncargs),
        dim = dim(calcfuncargs$inarray)[-1])
    }, silent = TRUE)
  }

  if (nrow(out$T0MEANS) > nrow(out$CINT)) { # then intoverpop used...
    nlatent <- nrow(out$CINT)
    out$T0MEANS <- out$T0MEANS[1:nlatent, 1, drop = FALSE]
    out$DRIFT <- out$DRIFT[1:nlatent, 1:nlatent, drop = FALSE]
    out$T0VAR <- out$T0VAR[1:nlatent, 1:nlatent, drop = FALSE]
    out$T0cov <- out$T0cov[1:nlatent, 1:nlatent, drop = FALSE]
  }

  ln <- latentNames
  mn <- manifestNames
  tdn <- TDpredNames
  dimnames(out$DRIFT) <- list(ln, ln)
  dimnames(out$DIFFUSIONcov) <- list(ln, ln)
  dimnames(out$DIFFUSION) <- list(ln, ln)
  dimnames(out$T0cov) <- list(ln, ln)
  dimnames(out$asymDIFFUSIONcov) <- list(ln, ln)
  rownames(out$CINT) <- ln
  rownames(out$MANIFESTMEANS) <- mn
  rownames(out$T0MEANS) <- ln
  dimnames(out$T0VAR) <- list(ln, ln)
  dimnames(out$LAMBDA) <- list(mn, ln)

  if (!is.null(e$pop_MANIFESTVAR)) {
    dimnames(out$MANIFESTVAR) <- list(mn, mn)
    dimnames(out$MANIFESTcov) <- list(mn, mn)
  }
  if (!is.null(e$pop_TDPREDEFFECT)) dimnames(out$TDPREDEFFECT) <- list(ln, tdn)

  out$MANIFESTVAR <- NULL
  out
}

.ctBackendSummaryMatrices <- function(fit, calcfunc = quantile,
  calcfuncargs = list(probs = 0.5), timeinterval = 1, ...) {
  model <- .ctBackendModel(fit)
  e <- .ctBackendPopArrays(fit, ...)
  .ctSummaryMatricesFromArrays(e, continuoustime = model$continuoustime,
    latentNames = model$latentNames, manifestNames = model$manifestNames,
    TDpredNames = model$TDpredNames, calcfunc = calcfunc,
    calcfuncargs = calcfuncargs, timeinterval = timeinterval)
}


# Fixed effects on the transformed scale --------------------------------------
#
# Each free raw parameter is defined by a transform into one or more model
# matrix cells, so the transformed value of raw parameter j *is* the value the
# engine wrote into a cell that j occupies. Reading it back is therefore exact
# and needs no second implementation of the transform; the only choice is which
# cell to read when a parameter is shared across several, and any of them gives
# the same number by construction.

# Jacobian blocks are derivatives of the model matrices, not model matrices, so
# a cell in one of these never reports a parameter's own value.
.ctBackendJacobianMatrices <- c("JAx", "Jy", "Jtd")

# The model-matrix cell whose value *is* each free parameter's population value.
#
# For most parameters that is the first cell the parameter occupies. An
# `intoverpop` random effect is the exception, and not a rare one: ctsem
# represents such a parameter as a carrier latent state, so the cell the
# parameter itself occupies is `T0MEANS[k]` -- the *raw* value, carrying none of
# the parameter's transform -- while the transform lives in whichever model
# matrix reads `state[k]`. Reading the carrier cell reports a raw number where
# the transformed one belongs: a random-effects CINT parameter came back ten
# times too small, its `10*param` transform missing. Stan resolves the same
# ambiguity the same way, following its `pr2` reference from the parameter's own
# matsetup row to the state-dependent one before applying `tform`.
.ctBackendFreeParameterCells <- function(fit) {
  spec <- .ctBackendSpec(fit)
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  free <- table[!is.na(table$parnumber), , drop = FALSE]
  free <- free[!duplicated(free$parnumber), , drop = FALSE]
  free <- free[order(free$parnumber), , drop = FALSE]

  nlatent <- spec$nlatent
  carrier <- if (is.null(nlatent)) rep(FALSE, nrow(free)) else {
    free$matrix %in% "T0MEANS" & free$row > nlatent
  }
  if (any(carrier)) {
    blank <- function(x) if (is.null(x)) rep("", nrow(table)) else replace(x, is.na(x), "")
    expressions <- paste(blank(table$predicttransform), blank(table$updatetransform),
      blank(table$tdtransform))
    reportable <- !table$matrix %in% .ctBackendJacobianMatrices
    references <- paste0("state[", free$row, "]")
    for (index in which(carrier)) {
      candidates <- which(reportable & grepl(references[index], expressions, fixed = TRUE))
      if (!length(candidates)) next
      free[index, c("matrix", "row", "col")] <- table[candidates[1L], c("matrix", "row", "col")]
    }
  }
  # The population sd / correlation parameters `.ctJuliaAugmentRandomEffects`
  # appends are free parameters too, but they belong in the random-effects
  # sections rather than among the fixed effects -- as they do for Stan, whose
  # `popmeans` covers only `nparams`.
  free$randomeffect <- free$parnumber %in% .ctBackendRandomEffectParameters(spec)
  free
}

.ctBackendRandomEffectParameters <- function(spec) {
  effects <- spec$random_effects
  if (is.null(effects) || !length(effects) || !nrow(effects)) return(integer())
  as.integer(effects$parameter)
}

# A parameter's own name, or `paramN` when the model gave the cell none.
#
# Four places wanted this rule and had three spellings of it, differing only in
# where the number came from -- `parnumber` here and in the raw-label table,
# `re_index` at a Laplace level. Same rule, one place.
.ctBackendParamLabel <- function(param, number) {
  label <- as.character(param)
  unnamed <- is.na(label)
  label[unnamed] <- paste0("param", as.integer(number)[unnamed])
  label
}

.ctBackendParameterNames <- function(cells) {
  .ctBackendParamLabel(cells$param, cells$parnumber)
}

# The value of every parameter's population cell, for a whole matrix of raw
# vectors: nsamples x nrow(cells). One engine call, whatever the sample count.
.ctBackendPopCellValues <- function(fit, samples, cells, layout, ...) {
  # Only these cells' rows cross the bridge. Callers ask for a handful of cells
  # out of the eighty-odd a model has, and the whole array used to come back
  # every time: five nodes of the random-effect quadrature spent 3.7 s
  # transferring 3.3 MB to keep 66 KB of it.
  flat <- .ctBackendParMatricesFlat(fit, t(samples),
    rows = .ctBackendCellPositions(cells, layout), ...)
  values <- t(flat)
  colnames(values) <- .ctBackendParameterNames(cells)
  values
}

# Where each cell sits in the flat parameter-matrix column.
.ctBackendCellPositions <- function(cells, layout) {
  index <- match(cells$matrix, layout$matrix)
  layout$offset[index] + (cells$col - 1L) * layout$nrow[index] + cells$row
}

.ctBackendPopCellsFromFlat <- function(flat, cells, layout) {
  values <- t(flat[.ctBackendCellPositions(cells, layout), , drop = FALSE])
  colnames(values) <- .ctBackendParameterNames(cells)
  values
}

.ctBackendPopMeanSamples <- function(fit, samples = NULL, ...) {
  if (!length(list(...))) {
    constrained <- .ctBackendConstrained(fit, samples)
    cells <- constrained$cells[!constrained$cells$randomeffect, , drop = FALSE]
    return(list(
      values = .ctBackendPopCellsFromFlat(constrained$flat, cells, constrained$layout),
      parnumber = as.integer(cells$parnumber)))
  }
  cells <- .ctBackendFreeParameterCells(fit)
  cells <- cells[!cells$randomeffect, , drop = FALSE]
  layout <- .ctBackendSummaryLayout(fit)
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  values <- .ctBackendPopCellValues(fit, samples, cells, layout, ...)
  list(values = values, parnumber = as.integer(cells$parnumber))
}

# Mean / sd / quantiles of a sample matrix, in Stan's summary column order. Uses
# base R rather than rstan's `monitor()`: with one point-estimate "sample" there
# is nothing to monitor, and with a normal-approximation posterior the
# convergence diagnostics `monitor()` adds would be meaningless anyway.
.ctBackendSampleSummary <- function(values, digits = 3, z = FALSE) {
  if (nrow(values) < 2L) {
    out <- data.frame(mean = as.numeric(values[1L, ]), row.names = colnames(values))
    return(round(out, digits))
  }
  probs <- c(.025, .5, .975)
  quantiles <- t(apply(values, 2L, stats::quantile, probs = probs, na.rm = TRUE))
  out <- data.frame(mean = colMeans(values, na.rm = TRUE),
    sd = apply(values, 2L, stats::sd, na.rm = TRUE),
    quantiles, check.names = FALSE, row.names = colnames(values))
  colnames(out) <- c("mean", "sd", "2.5%", "50%", "97.5%")
  if (isTRUE(z)) out$z <- out$mean / out$sd
  round(out, digits)
}

# Gauss-Hermite nodes and weights for a standard normal, by Golub-Welsch on the
# probabilists' Hermite recurrence. Used to integrate a parameter's transform
# over its population distribution (see .ctBackendRandomEffectSummary), where a
# fixed quadrature is both cheaper and steadier than Stan's 5000 random draws.
#
# Five nodes integrates a degree-9 polynomial exactly, which for a second moment
# of a smooth transform is far inside Stan's own Monte Carlo error; each node
# costs one engine call over the whole posterior, so this is also the term that
# sets what a summary costs.
.ctBackendGaussHermite <- function(nodes = 5L) {
  index <- seq_len(nodes - 1L)
  jacobi <- matrix(0, nodes, nodes)
  jacobi[cbind(index, index + 1L)] <- sqrt(index)
  jacobi[cbind(index + 1L, index)] <- sqrt(index)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(node = decomposition$values, weight = decomposition$vectors[1L, ]^2)
}

# Random-effects standard deviations and raw correlations -- Stan's `popsd` and
# `rawpopcorr` sections.
#
# Both come from the population covariance of the carrier states, `T0cov`, but
# they report it on different scales, and that is the whole subtlety:
#
#   * The **correlations** are between the *raw* parameters, so the state
#     scaling `.ctJuliaAugmentRandomEffects` folds into the sd transform cancels
#     and `cov2cor(T0cov)` is already the reported quantity.
#
#   * The **standard deviations** are of the *transformed* parameter, which for
#     a nonlinear transform is not the transform of the standard deviation. Stan
#     draws 5000 subjects from the raw population distribution, pushes each
#     through the transform and takes the sd; this does the same integral by
#     Gauss-Hermite quadrature over the same distribution, so a linear transform
#     is exact and a nonlinear one is far steadier than a random cloud. The
#     quadrature displaces the *raw* parameter and reads the population cell
#     back through the engine, so whatever the transform is, it is applied once,
#     by the code that owns it.
#
# All the varying parameters are displaced together on a shared node, which is
# what keeps this to one engine call per node rather than one per parameter:
# only each parameter's own marginal spread is read, so the perfect correlation
# a shared node induces between them never enters an answer.
.ctBackendRandomEffectDraws <- function(fit, samples, cells, layout, flat) {
  spec <- .ctBackendSpec(fit)
  populations <- if (!is.null(spec$laplace)) {
    .ctBackendLaplacePopulations(fit, spec, samples)
  } else {
    p <- .ctBackendAugmentedPopulation(spec, samples, layout, flat)
    if (is.null(p)) NULL else list(p)
  }
  if (is.null(populations) || !length(populations)) return(NULL)

  # One set of results per level. With a single level this is the ordinary
  # random-effects summary and `levels` has one entry; with a study level above
  # the subjects it has two, and each is a population in its own right -- a
  # study sd is the spread between studies, not between subjects, and averaging
  # them together would describe neither.
  out <- list(levels = list())
  for (population in populations) {
    out$levels[[length(out$levels) + 1L]] <-
      .ctBackendRandomEffectLevel(fit, population, cells, layout, samples)
  }
  # The innermost level stays at the top level of the result, so every existing
  # single-level caller of `$popsd` and `$rawpopcorr` keeps working unchanged.
  out$popsd <- out$levels[[1]]$popsd
  out$rawpopcorr <- out$levels[[1]]$rawpopcorr
  return(out)
}

# The quadrature for one level. Split out because it is identical at every
# level: what differs is only which raw parameters vary and how spread out they
# are, both of which arrive in `population`.
.ctBackendRandomEffectLevel <- function(fit, population, cells, layout, samples) {
  parnumber <- population$parnumber
  parname <- population$param
  rawsd <- population$rawsd

  out <- list(level = population$level)
  column <- match(parnumber, cells$parnumber)
  quadrature <- .ctBackendGaussHermite()
  # Narrowed before the call, not after: asking for every cell and keeping
  # `column` was five full transfers of the whole parameter-matrix array per
  # level, ~97% of it discarded on arrival.
  wanted <- cells[column, , drop = FALSE]
  # One engine call for all five nodes, not one per node.
  #
  # Every node asks for the same cells of the same posterior, displaced by a
  # different multiple of the same sd, so the five calls differ only in the
  # numbers they send. `ctsem_parameter_matrices` already takes a matrix of
  # raw vectors and treats each column independently, so stacking the five
  # displaced posteriors and splitting the result afterwards computes exactly
  # the same values in exactly the same way.
  #
  # It is worth doing because this phase is not compute bound. Measured on
  # dev1: a `ctsem_parameter_matrices` call costs about 0.25 s of which about
  # 0.01 s is the engine -- the rest is JuliaConnectoR round trips, whose cost
  # is per call and nearly independent of how much is in each. Five narrow
  # calls cost 1.27 s; one call five times as wide costs about a fifth of
  # that. `rows` keeps the reply narrow either way.
  ndraws <- nrow(samples)
  stacked <- do.call(rbind, lapply(quadrature$node, function(node) {
    perturbed <- samples
    perturbed[, parnumber] <- perturbed[, parnumber, drop = FALSE] + rawsd * node
    perturbed
  }))
  together <- .ctBackendPopCellValues(fit, stacked, wanted, layout)
  displaced <- lapply(seq_along(quadrature$node), function(index) {
    together[seq_len(ndraws) + (index - 1L) * ndraws, , drop = FALSE]
  })
  centre <- Reduce(`+`, Map(function(value, weight) value * weight,
    displaced, quadrature$weight))
  spread <- Reduce(`+`, Map(function(value, weight) weight * (value - centre)^2,
    displaced, quadrature$weight))
  spread <- matrix(sqrt(pmax(spread, 0)), nrow = nrow(samples))
  colnames(spread) <- parname
  out$popsd <- spread

  if (!is.null(population$rawcorr) && ncol(population$rawcorr)) {
    lower <- which(lower.tri(diag(length(parnumber))), arr.ind = TRUE)
    correlation <- population$rawcorr
    colnames(correlation) <- paste0(parname[lower[, 1L]], "__", parname[lower[, 2L]])
    out$rawpopcorr <- correlation
  }
  out
}

# The augmented route's population sds, and the raw parameter each one belongs
# to. NULL when the model has none.
#
# `.ctJuliaAugmentRandomEffects` records an sd against its carrier state's
# T0VAR row, so the parameter that state varies is whatever its T0MEANS cell
# holds. The random-effects summary and ctSubjectPars both need that lookup and
# had a copy each, which is one lookup too many for a mapping this indirect.
.ctBackendAugmentedSds <- function(spec) {
  effects <- spec$random_effects
  if (is.null(effects) || !length(effects) || !nrow(effects)) return(NULL)
  sds <- effects[effects$type %in% "sd", , drop = FALSE]
  if (!nrow(sds)) return(NULL)
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  t0means <- table[table$matrix %in% "T0MEANS" & table$col == 1L, , drop = FALSE]
  position <- match(sds$row, t0means$row)
  list(sds = sds, parnumber = as.integer(t0means$parnumber[position]),
    param = t0means$param[position])
}

# Which parameters this fit can report per subject.
#
# The two routes reach the same set from different starting points -- the
# augmented one through its carrier states, the Laplace one through its own
# `re_index` -- and both then added the TI-predictor effects, filtered against
# the reportable cells, and raised the same message when nothing survived. Only
# the first line differs, so only the first line is written twice.
.ctBackendVaryingParameters <- function(spec, cells) {
  varying <- if (!is.null(spec$laplace)) as.integer(spec$laplace$re_index) else {
    augmented <- .ctBackendAugmentedSds(spec)
    if (is.null(augmented)) integer() else augmented$parnumber
  }
  if (!is.null(spec$ti_effects) && nrow(spec$ti_effects)) {
    varying <- c(varying, as.integer(spec$ti_effects$parameter))
  }
  varying <- sort(unique(varying[!is.na(varying)]))
  varying <- varying[varying %in% cells$parnumber]
  if (!length(varying)) stop("No individually varying parameters in model!", call. = FALSE)
  varying
}

# The augmented route's population scales and correlations, read out of the
# filtered T0 covariance the carrier states live in. `scale` divides out the
# state-unit factor `.ctJuliaAugmentRandomEffects` folded into the sd transform,
# because what is wanted here is the sd on the *raw parameter* scale.
.ctBackendAugmentedPopulation <- function(spec, samples, layout, flat) {
  augmented <- .ctBackendAugmentedSds(spec)
  if (is.null(augmented)) return(NULL)
  sds <- augmented$sds
  parnumber <- augmented$parnumber
  if (any(is.na(parnumber))) return(NULL)
  parname <- .ctBackendParamLabel(augmented$param, parnumber)

  t0cov <- .ctBackendReshape(flat, layout, match("T0cov", layout$matrix))
  variance <- matrix(vapply(sds$row, function(row) t0cov[, row, row],
    numeric(nrow(samples))), nrow = nrow(samples))
  scale <- if (is.null(sds$scale)) rep(1, nrow(sds)) else as.numeric(sds$scale)
  rawsd <- sweep(sqrt(pmax(variance, 0)), 2L, scale, "/")

  rawcorr <- NULL
  if (nrow(sds) > 1L) {
    lower <- which(lower.tri(diag(nrow(sds))), arr.ind = TRUE)
    rawcorr <- matrix(vapply(seq_len(nrow(lower)), function(entry) {
      i <- sds$row[lower[entry, 1L]]
      j <- sds$row[lower[entry, 2L]]
      t0cov[, i, j] / sqrt(t0cov[, i, i] * t0cov[, j, j])
    }, numeric(nrow(samples))), nrow = nrow(samples))
  }
  list(parnumber = parnumber, param = parname, rawsd = rawsd, rawcorr = rawcorr,
    level = spec$model$subjectIDname)
}

# Subject parameters on the Laplace route.
#
# The augmented route reads these off the filter, because there the random
# effects *are* states and the filter estimates them. Here they come from the
# inner modes instead, which is the same conditional-mode quantity by a
# different route: the engine assembles each subject's raw vector -- population
# values, plus its own random effects, plus its own TI-predictor effects -- and
# it is pushed through the model's transforms by the same code that produces the
# population values, so a subject's parameter and the population parameter are
# reported on the same scale by construction.
#
# Point estimate only. Posterior draws would need the inner mode re-solved at
# every draw, which is a real computation rather than a lookup, and returning
# the point-estimate modes against varying population draws would silently
# understate the spread it is being asked for.
.ctBackendLaplaceSubjectPars <- function(fit, spec, pointest = TRUE,
  nsamples = "all") {
  cells <- .ctBackendFreeParameterCells(fit)
  cells <- cells[!cells$randomeffect, , drop = FALSE]
  varying <- .ctBackendVaryingParameters(spec, cells)

  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  estimate <- .ctJuliaNumericVector(fit$estimate$raw)
  index <- match(varying, cells$parnumber)
  selected <- cells[index, , drop = FALSE]
  layout <- .ctBackendSummaryLayout(fit)
  parnames <- .ctBackendParameterNames(cells)[index]
  alphabetical <- order(parnames)

  draws <- if (isTRUE(pointest)) NULL else .ctBackendRawSamples(fit)
  if (!is.null(draws) && identical(nrow(draws), 1L)) draws <- NULL
  if (!is.null(draws) && !identical(nsamples, "all")) {
    keep <- unique(round(seq(1, nrow(draws), length.out = min(nrow(draws),
      as.integer(nsamples)))))
    draws <- draws[keep, , drop = FALSE]
  }

  if (is.null(draws)) {
    subject_raw <- .ctBackendJuliaValue(module$ctsem_laplace_subject_values(
      objective, estimate))
    subject_raw <- matrix(as.numeric(subject_raw), ncol = length(fit$estimate$raw))
    values <- .ctBackendPopCellValues(fit, subject_raw, selected, layout)
    out <- array(as.numeric(values[, alphabetical, drop = FALSE]),
      dim = c(1L, nrow(subject_raw), length(varying)))
    dimnames(out) <- list(iter = 1L, subject = seq_len(nrow(subject_raw)),
      param = parnames[alphabetical])
    return(out)
  }

  # Draws. The population factors are rebuilt exactly at every draw; only the
  # random-effect mode is linearised around the estimate, which is what makes
  # this a matrix-vector product per draw rather than a Newton solve. See
  # `ctsem_laplace_subject_values`: it is an approximation, and it is documented
  # as one wherever it surfaces.
  raw <- .ctBackendJuliaValue(module$ctsem_laplace_subject_values(objective,
    JuliaConnectoR::juliaPut(as.matrix(draws)), estimate))
  nsubjects <- length(spec$subject_starts)
  raw <- array(as.numeric(raw), dim = c(nrow(draws), nsubjects,
    length(fit$estimate$raw)))
  flat <- matrix(aperm(raw, c(2L, 1L, 3L)), nrow = nsubjects * nrow(draws))
  values <- .ctBackendPopCellValues(fit, flat, selected, layout)
  out <- array(NA_real_, dim = c(nrow(draws), nsubjects, length(varying)))
  for (position in seq_along(alphabetical)) {
    out[, , position] <- matrix(values[, alphabetical[position]],
      nrow = nrow(draws), ncol = nsubjects, byrow = TRUE)
  }
  dimnames(out) <- list(iter = seq_len(nrow(draws)), subject = seq_len(nsubjects),
    param = parnames[alphabetical])
  out
}

# The Laplace route's population scales and correlations. There is no carrier
# state to read them off, and no state-unit rescaling to undo: the population
# covariance is built on the raw parameter scale in the first place, so the
# engine is asked for it directly, once for the whole posterior sample.
.ctBackendLaplacePopulations <- function(fit, spec, samples) {
  laplace <- spec$laplace
  if (is.null(laplace) || !laplace$nrandom) return(NULL)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  draws <- JuliaConnectoR::juliaPut(as.matrix(samples))
  out <- list()
  for (l in seq_along(laplace$levels)) {
    level <- laplace$levels[[l]]
    if (!level$nrandom) next
    result <- JuliaConnectoR::juliaGet(module$ctsem_laplace_population(
      objective, draws, as.integer(l)))
    parname <- .ctBackendParamLabel(level$param, level$re_index)
    out[[length(out) + 1L]] <- list(
      parnumber = as.integer(level$re_index), param = parname,
      rawsd = matrix(as.numeric(result$sd), nrow = nrow(samples)),
      rawcorr = matrix(as.numeric(result$correlation), nrow = nrow(samples)),
      level = level$name)
  }
  if (!length(out)) NULL else out
}

# Time-independent predictor effects, on the transformed parameters -- Stan's
# `linearTIPREDEFFECT`.
#
# The raw coefficient is an effect on the *unconstrained* parameter, which is
# not a quantity anyone reads off a summary. Stan reports the effect the
# covariate has on the parameter itself, linearised at the population mean: the
# transform evaluated a hundredth of an effect either side of the raw mean,
# differenced, and scaled back up. This does exactly that, through the engine's
# own transforms rather than a second copy of them, and one predictor at a time
# so the cost is two engine calls per predictor rather than two per effect.
.ctBackendTipredDraws <- function(fit, samples, cells, layout) {
  spec <- .ctBackendSpec(fit)
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) return(NULL)
  predictor_names <- .ctBackendModel(fit)$TIpredNames
  parameter_names <- stats::setNames(.ctBackendParameterNames(cells), cells$parnumber)

  predictors <- sort(unique(effects$predictor))
  perpredictor <- lapply(predictors, function(predictor)
    effects[effects$predictor %in% predictor, , drop = FALSE])
  ndraws <- nrow(samples)

  # One engine call for every predictor and both directions, for the same
  # reason the random-effect quadrature takes one: the calls differ only in the
  # numbers they send, and on this backend the cost is per call rather than per
  # byte. This was two calls per predictor -- six on a three-predictor model,
  # measured at 1.72 s, which was 67% of the whole constrain step.
  #
  # The reply is still narrowed by `rows`, now to the union of the cells any
  # predictor asks for. In the ordinary case every predictor perturbs the same
  # parameters, so that union is one predictor's set and the batched call
  # returns exactly the bytes the separate calls did between them. A model
  # whose predictors touch disjoint parameters would send more, which is the
  # trade this makes knowingly: a round trip costs 40-80 ms before it carries
  # anything, and the bridge moves about 1.7 MB/s.
  parameters <- sort(unique(unlist(lapply(perpredictor, function(rows) rows$parameter))))
  wanted <- cells[match(parameters, cells$parnumber), , drop = FALSE]

  blocks <- vector("list", 2L * length(predictors))
  position <- 0L
  for (rows in perpredictor) {
    step <- samples[, rows$coefficient, drop = FALSE] * .01
    for (direction in c(1, -1)) {
      perturbed <- samples
      perturbed[, rows$parameter] <- perturbed[, rows$parameter, drop = FALSE] +
        direction * step
      position <- position + 1L
      blocks[[position]] <- perturbed
    }
  }
  together <- .ctBackendPopCellValues(fit, do.call(rbind, blocks), wanted, layout)

  # Block 2i-1 is predictor i displaced up and block 2i the same displaced
  # down, in the order they were stacked. `take` puts the union's columns back
  # into this predictor's own order, which is what the names below assume.
  linear <- lapply(seq_along(predictors), function(index) {
    rows <- perpredictor[[index]]
    take <- match(rows$parameter, parameters)
    up <- together[seq_len(ndraws) + (2L * index - 2L) * ndraws, take, drop = FALSE]
    down <- together[seq_len(ndraws) + (2L * index - 1L) * ndraws, take, drop = FALSE]
    values <- matrix((up - down) / .02, nrow = ndraws)
    colnames(values) <- paste0("tip_", predictor_names[predictors[index]], "_",
      parameter_names[as.character(rows$parameter)])
    values
  })
  do.call(cbind, linear)
}

# Per-subject parameter values -- ctSubjectPars() for a julia backend fit.
#
# Every parameter that varies over subjects, for every subject, read out of the
# `subj_*` matrices the filter already produced. Which parameters those are is a
# property of the specification (a random effect, a TI predictor effect, or
# both), and where each one lives is the same population cell the fixed-effects
# summary reads -- so a random-effects CINT parameter is reported from `subj_
# CINT`, transform applied, rather than from its raw carrier state.
.ctBackendSubjectPars <- function(fit, pointest = TRUE, nsamples = "all") {
  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$laplace)) {
    return(.ctBackendLaplaceSubjectPars(fit, spec, pointest, nsamples))
  }
  cells <- .ctBackendFreeParameterCells(fit)
  cells <- cells[!cells$randomeffect, , drop = FALSE]
  varying <- .ctBackendVaryingParameters(spec, cells)

  if (isTRUE(pointest)) fit$estimate$rawposterior <- NULL
  extracted <- .ctBackendExtract(fit, subjectMatrices = TRUE, nsamples = nsamples)

  index <- match(varying, cells$parnumber)
  parnames <- .ctBackendParameterNames(cells)[index]
  reference <- extracted[[paste0("subj_", cells$matrix[index[1L]])]]
  out <- array(NA_real_, dim = c(dim(reference)[1L], dim(reference)[2L], length(varying)))
  for (position in seq_along(index)) {
    cell <- index[position]
    values <- extracted[[paste0("subj_", cells$matrix[cell])]]
    out[, , position] <- values[, , cells$row[cell], cells$col[cell]]
  }
  alphabetical <- order(parnames)
  out <- out[, , alphabetical, drop = FALSE]
  dimnames(out) <- list(iter = seq_len(dim(out)[1L]), subject = seq_len(dim(out)[2L]),
    param = parnames[alphabetical])
  out
}


# Standardised residual covariance -- Stan's `residCovStd`, from the same
# quantity (the filter's prior errors at the estimate) reached through the
# backend-neutral accessors.
.ctBackendResidCovStd <- function(object, digits = 3) {
  errors <- .ctFitPriorErrors(object)
  if (is.null(errors)) return(NULL)
  observed <- .ctFitObservedY(object)
  obscov <- stats::cov(observed, use = "pairwise.complete.obs")
  standardise <- diag(1 / sqrt(diag(obscov)), ncol(obscov))
  rescov <- stats::cov(matrix(errors, ncol = ncol(obscov)),
    use = "pairwise.complete.obs")
  unavailable <- which(is.na(rescov))
  rescov[unavailable] <- 0
  out <- round(standardise %*% rescov %*% standardise, digits)
  out[unavailable] <- NA
  manifest <- .ctBackendModel(object)$manifestNames
  dimnames(out) <- list(manifest, manifest)
  # Carried on the result so `print.summary` can say it without recomputing it.
  attr(out, "conditioning") <- attr(errors, "conditioning")
  out
}

# The constrain step, done once ------------------------------------------------
#
# A Stan fit carries `stanfit$transformedpars`: its draws pushed through the
# model's transforms once, at fit time, so `summary()` is a collapse over
# something already computed. This is the same thing for the julia backend, and
# for the same reason -- without it every `summary()` call re-materialized every
# model matrix for every draw, which on a 28-parameter model was about sixteen
# seconds, every time.
#
# What it holds is everything downstream reads that costs an engine call:
#
#   flat        every model matrix for every draw, in the engine's flat layout
#   popsd       the population sd of each varying parameter, per draw
#   rawpopcorr  the raw population correlations, per draw
#   tipreds     the linearised TI predictor effects, per draw
#
# The last two are here rather than derived on demand because they are read at
# *displaced* parameter values -- the quadrature nodes and the linearisation
# steps -- so they cannot be recovered from `flat` afterwards. Stan computes its
# equivalents (`popsd`, `linearTIPREDEFFECT`) in generated quantities at the
# same moment, for the same reason.
#
# `samples` is kept alongside so the cache can be checked rather than trusted:
# anything that changes the posterior (a fresh `ctOptimUncertainty()`, say)
# leaves a cache that no longer matches, and is recomputed rather than silently
# describing the previous draws.
.ctBackendConstrain <- function(fit, samples = NULL) {
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  layout <- .ctBackendSummaryLayout(fit)
  cells <- .ctBackendFreeParameterCells(fit)
  flat <- .ctBackendParMatricesFlat(fit, t(samples))
  randomeffects <- .ctBackendRandomEffectDraws(fit, samples = samples, cells = cells,
    layout = layout, flat = flat)
  list(samples = samples, layout = layout, cells = cells, flat = flat,
    popsd = randomeffects$popsd, rawpopcorr = randomeffects$rawpopcorr,
    randomeffectlevels = randomeffects$levels,
    tipreds = .ctBackendTipredDraws(fit, samples = samples, cells = cells,
      layout = layout))
}

# The cached constrain step if it describes these draws, otherwise a fresh one.
.ctBackendConstrained <- function(fit, samples = NULL) {
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  cached <- fit$transformedpars
  if (!is.null(cached) && identical(dim(cached$samples), dim(samples)) &&
      isTRUE(all.equal(cached$samples, samples, check.attributes = FALSE))) {
    return(cached)
  }
  .ctBackendConstrain(fit, samples)
}

.ctBackendSummary <- function(object, timeinterval = 1, digits = 3, parmatrices = TRUE,
  residualcov = TRUE, ...) {
  has_posterior <- !is.null(object$estimate$rawposterior)

  out <- list()

  constrained <- .ctBackendConstrained(object)
  samples <- constrained$samples
  cells <- constrained$cells
  layout <- constrained$layout
  flat <- constrained$flat

  if (isTRUE(residualcov)) {
    residCovStd <- .ctBackendResidCovStd(object, digits = digits)
    if (!is.null(residCovStd)) {
      out$residCovStd <- residCovStd
      # As a field rather than only an attribute on the matrix: attributes do
      # not reliably survive the way a summary is assembled and printed, and
      # this is not decoration. The augmented, Laplace and sampled routes
      # condition the residuals on different things -- effects learned
      # observation by observation, effects from a subject's whole record, and
      # the posterior mean respectively -- so the same table means three
      # different quantities and a reader comparing fits needs to know which.
      out$residCovStdConditioning <- attr(residCovStd, "conditioning")
    }
  }

  if (length(constrained$randomeffectlevels) > 1L) {
    # More than one level, so nothing goes in the unlabelled slot: a bare
    # "Random-effects correlations" table would leave the reader guessing
    # whether it described spread between subjects or between studies. Each
    # level gets its own named section instead, and `$randomEffects` carries
    # them all for programmatic use.
    out$randomEffects <- constrained$randomeffectlevels
    for (lv in constrained$randomeffectlevels) {
      if (!is.null(lv$rawpopcorr)) {
        out[[paste0("rawpopcorr.", lv$level)]] <-
          .ctBackendSampleSummary(lv$rawpopcorr, digits = digits)
      }
    }
  } else if (!is.null(constrained$rawpopcorr)) {
    out$rawpopcorr <- .ctBackendSampleSummary(constrained$rawpopcorr,
      digits = digits, z = nrow(samples) > 1L)
    out$rawpopcorrNote <-
      "These reflect correlations between the raw / unconstrained parameters."
  }

  if (!is.null(constrained$tipreds)) {
    out$tipreds <- .ctBackendSampleSummary(constrained$tipreds, digits = digits,
      z = nrow(samples) > 1L)
    out$tipredsNote <- "Approximate (linearised) effects on the transformed parameters."
  }

  if (isTRUE(parmatrices)) {
    # One materialization, five collapses -- not five calls to
    # .ctBackendSummaryMatrices(), each of which would ask the engine for the
    # same arrays again.
    model <- .ctBackendModel(object)
    arrays <- .ctBackendPopArraysFromFlat(flat, layout, .ctBackendSpec(object))
    collapse <- function(calcfunc, calcfuncargs) {
      .ctSummaryMatricesFromArrays(arrays, continuoustime = model$continuoustime,
        latentNames = model$latentNames, manifestNames = model$manifestNames,
        TDpredNames = model$TDpredNames, calcfunc = calcfunc,
        calcfuncargs = calcfuncargs, timeinterval = timeinterval)
    }
    collapsed <- list(Mean = collapse(mean, list()))
    if (has_posterior) {
      collapsed$sd <- collapse(stats::sd, list(na.rm = TRUE))
      for (probability in c(.025, .5, .975)) {
        collapsed[[paste0(probability * 100, "%")]] <-
          collapse(stats::quantile, list(probs = probability))
      }
      names(collapsed) <- c("Mean", "sd", "2.5%", "50%", "97.5%")
    }
    d <- data.frame(ctModelUnlist(collapsed$Mean, matnames = names(collapsed$Mean)))
    colnames(d)[colnames(d) %in% "value"] <- "Mean"
    for (name in setdiff(names(collapsed), "Mean")) {
      d[[name]] <- round(ctModelUnlist(collapsed[[name]], names(collapsed$Mean))$value, digits)
    }
    d$param <- NULL
    d$Mean <- round(d$Mean, digits)
    d <- d[!d$matrix %in% c("DIFFUSION", "T0VAR"), ]
    # A model with no PARS still carries a 1x1 PARS matrix, which printed as a
    # row of zeros among the estimates -- a matrix the user never mentioned,
    # reported as though it were a result. Dropped when it says nothing.
    parsrows <- d$matrix %in% "PARS"
    if (any(parsrows)) {
      values <- as.matrix(d[parsrows, intersect(names(d),
        c("Mean", "sd", "2.5%", "50%", "97.5%")), drop = FALSE])
      if (all(!is.finite(values) | values == 0)) d <- d[!parsrows, , drop = FALSE]
    }
    out$parmatrices <- d

    # The note comes from the shared classifier (R/ctContextDependence.R)
    # rather than from the engine's `statedep` count, because that count does
    # not distinguish an individually varying parameter's carrier state -- which
    # is reported exactly and needs no caveat -- from a genuinely dynamic
    # reference, which does.
    out$parmatNote <- .ctContextNote(.ctFitConditionalCells(object),
      .ctContextPopLabel, .ctContextRemedy(object))
  }

  if (length(constrained$randomeffectlevels) > 1L) {
    for (lv in constrained$randomeffectlevels) {
      if (!is.null(lv$popsd)) {
        out[[paste0("popsd.", lv$level)]] <-
          .ctBackendSampleSummary(lv$popsd, digits = digits)
      }
    }
  } else if (!is.null(constrained$popsd)) {
    out$popsd <- .ctBackendSampleSummary(constrained$popsd, digits = digits)
  }

  fixed <- cells[!cells$randomeffect, , drop = FALSE]
  out$popmeans <- .ctBackendSampleSummary(
    .ctBackendPopCellsFromFlat(flat, fixed, layout), digits = digits)
  out$popNote <- paste0("Population values on the transformed scale. ",
    "Covariance parameters appear in sd / unconstrained correlation form; ",
    "see System Matrices (or ctSummaryMatrices()) for cor/cov.")

  logposterior <- object$estimate$logposterior
  if (is.null(logposterior)) logposterior <- object$estimate$loglik
  out$logposterior <- logposterior
  # Without priors these are the same number; the shared print method (not
  # this builder) skips the duplicate row -- see
  # summaryCtStanFitLoglikDuplicatesPosterior() in R/summary.ctStanFit.R.
  # The likelihood is always returned, priors or not.
  out$loglik <- object$estimate$loglik
  out$npars <- length(object$estimate$raw)
  out$aic <- 2 * out$npars - 2 * object$estimate$loglik
  # Named for what they are. "Number of samples" on an optimised fit meant the
  # uncertainty draws and read as MCMC samples -- a fit that never sampled
  # reporting a sample count.
  if (has_posterior) {
    out[[if (isTRUE(object$args$optimize %in% FALSE)) "nsamples" else "ndraws"]] <-
      nrow(object$estimate$rawposterior)
  }
  out$uncertaintyNote <- if (has_posterior) {
    paste0("Julia backend; intervals from ctOptimUncertainty(uncertainty='",
      object$uncertainty$settings$method, "') draws pushed through the transforms.")
  } else {
    paste0("Julia backend; point estimates only. ",
      "Run ctOptimUncertainty() for standard errors and intervals.")
  }

  # Matrices become data frames exactly as summary.ctStanFit does, so the shared
  # print method treats a one-manifest residual covariance as a table rather
  # than as a bare scalar.
  out <- lapply(out, function(x) {
    if ("matrix" %in% class(x)) x <- data.frame(x, check.names = FALSE)
    roundSummaryCtStanFitValue(x, digits = digits)
  })
  attr(out, "digits") <- digits
  # Two classes, one print method: the sections are named as
  # print.summary.ctStanFit expects, so it prints these unchanged and there is
  # no second printer to keep in step. The extra class is there for code that
  # wants to tell a backend summary apart -- notably, one whose intervals come
  # from a normal approximation rather than from HMC.
  class(out) <- c("summary.ctBackendFit", "summary.ctStanFit")
  out
}

# Prior prediction errors, and what they were conditioned on.
#
# `errprior = y - yprior` is the only part of a filter pass any summary reads --
# `.ctBackendResidCovStd()` here and `summary.ctStanFit()` are the two callers,
# and both want exactly this. Caching the whole filter output for it was 18
# arrays where one was wanted, and on the julia backend that is a bridge
# transfer rather than a memory cost: `etacov` alone is `nrows * nlatent^2`.
#' @keywords internal
.ctBackendPriorErrors <- function(fit) {
  observed <- .ctFitObservedY(fit)
  # `y` from the engine is the *prior prediction*; the observations are already
  # on this side. So only that one array need cross, not the covariances.
  raw <- try(.ctBackendKalmanRaw(fit, as.numeric(fit$estimate$raw),
    subjectmatrices = FALSE, fields = "y"), silent = TRUE)
  if (inherits(raw, "try-error") || is.null(raw$y)) return(NULL)
  # The engine hands this back manifest-major, and sometimes with a leading
  # sample dimension, so the shape is rebuilt from the length rather than
  # trusted. `.ctBackendPriorErrors` is checked against a full filter pass in
  # the test suite, which is what pins the ordering.
  # The engine stacks three predictions -- prior, updated, smoothed -- as
  # `3 x nrows x nmanifest`. The first is the one a prior residual is defined
  # against.
  predicted <- raw$y
  if (length(dim(predicted)) == 3L) {
    predicted <- predicted[1L, , , drop = TRUE]
  }
  predicted <- matrix(as.numeric(predicted), nrow = nrow(observed))
  if (!identical(dim(predicted), dim(observed))) return(NULL)
  errors <- observed - predicted
  attr(errors, "conditioning") <- .ctFitResidualConditioning(fit)
  errors
}

# Whichever backend and route produced the fit, say what the residuals are
# conditional on. The three answers are genuinely different quantities and a
# reader comparing them across fits needs to know which they have.
#' @keywords internal
.ctFitResidualConditioning <- function(fit) {
  sampled <- !is.null(fit$sample) ||
    identical(fit$uncertainty$settings$method, "sampling")
  laplace <- !is.null(fit$model_spec$laplace)
  if (sampled && laplace) {
    return(paste("at the posterior mean of the population parameters, with each",
      "subject's random effects re-estimated from their whole record"))
  }
  if (laplace) {
    return(paste("on each subject's random effects estimated from their whole",
      "record, so these are smoothed rather than filtered residuals"))
  }
  paste("on random effects carried as states and learned observation by",
    "observation, so early residuals for a subject are larger than late ones")
}

# The prior errors, from wherever this fit keeps them.
#' @keywords internal
.ctFitPriorErrors <- function(fit) {
  if (!is.null(fit$priorerrors)) return(fit$priorerrors)
  # Older objects, and the Stan path, keep a whole filter pass.
  cached <- if (!is.null(fit$kalman)) fit$kalman else fit$stanfit$kalman
  if (!is.null(cached$errprior)) {
    errors <- matrix(cached$errprior, ncol = ncol(.ctFitObservedY(fit)))
    attr(errors, "conditioning") <- .ctFitResidualConditioning(fit)
    return(errors)
  }
  if (inherits(fit, "ctJuliaFit")) return(.ctBackendPriorErrors(fit))
  NULL
}
