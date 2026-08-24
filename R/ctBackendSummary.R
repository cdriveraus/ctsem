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
#   * They are point estimates unless `ctOptimUncertainty()` has been run. With
#     uncertainty they carry `estimate$rawposterior`, and every interval below
#     comes from pushing those draws through the transforms -- which is what
#     Stan's own optimized-and-sampled path does too. Without it there is one
#     "sample", and the interval columns are omitted rather than filled with a
#     zero-width interval that would read as certainty.
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
# cells are state dependent. Cheap, but queried once per call rather than per
# sample; it is a property of the model, not of the parameter values.
.ctBackendSummaryLayout <- function(fit) {
  spec <- .ctBackendSpec(fit)
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
  list(matrix = as.character(raw$matrix), nrow = as.integer(raw$nrow),
    ncol = as.integer(raw$ncol), offset = as.integer(raw$offset),
    size = as.integer(raw$size)[1L], nlatent = as.integer(raw$nlatent)[1L],
    nmanifest = as.integer(raw$nmanifest)[1L], statedep = statedep)
}

# `raw` is npar x nsamples; the result is (flat layout) x nsamples. The whole
# posterior travels in one call: the Julia bridge marshals a numeric array in
# one transfer but a list element by element, and per-sample calls were the
# single largest avoidable cost measured in this backend.
.ctBackendParMatricesFlat <- function(fit, raw, tipreds = NULL, state = NULL,
  time = 0, dt = 0) {
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
.ctBackendPopArrays <- function(fit, samples = NULL, tipreds = NULL, state = NULL,
  time = 0, dt = 0) {
  layout <- .ctBackendSummaryLayout(fit)
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  flat <- .ctBackendParMatricesFlat(fit, t(samples), tipreds = tipreds, state = state,
    time = time, dt = dt)
  spec <- .ctBackendSpec(fit)
  out <- lapply(seq_along(layout$matrix), function(index) {
    .ctBackendTrimAugmented(.ctBackendReshape(flat, layout, index), layout$matrix[index],
      spec, margin = c(2L, 3L))
  })
  names(out) <- paste0("pop_", layout$matrix)
  out
}

.ctBackendExtract <- function(object, subjectMatrices = FALSE, nsamples = "all",
  subjects = "all", ...) {
  samples <- .ctBackendRawSamples(object)
  if (!identical(nsamples, "all")) {
    wanted <- min(nrow(samples), as.integer(nsamples)[1L])
    samples <- samples[round(seq(1, nrow(samples), length.out = wanted)), , drop = FALSE]
  }
  arrays <- .ctBackendPopArrays(object, samples = samples, ...)
  popmeans <- .ctBackendPopMeanSamples(object, samples = samples)
  out <- c(list(rawpars = samples, popmeans = popmeans$values,
    loglik = object$estimate$loglik, gradient = object$estimate$gradient,
    subject_loglik = object$estimate$subject_loglik), arrays)

  # Subject matrices need the filter, not just the transforms: an individually
  # varying parameter is an augmented latent state, so a subject's value for it
  # is only known once that subject's data has been filtered and smoothed.
  if (isTRUE(subjectMatrices)) {
    spec <- .ctBackendKalmanSpec(object, subjects = subjects)
    flat <- NULL
    for (iteration in seq_len(nrow(samples))) {
      scores <- .ctBackendKalmanRaw(spec, samples[iteration, ], subjectmatrices = TRUE)
      if (is.null(flat)) {
        flat <- array(0, dim = c(nrow(samples), dim(scores$subject_matrices)))
      }
      flat[iteration, , ] <- scores$subject_matrices
    }
    out <- c(out, .ctBackendSubjectMatrices(spec, flat))
  }
  out
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

.ctBackendFreeParameterCells <- function(fit) {
  spec <- .ctBackendSpec(fit)
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  free <- table[!is.na(table$parnumber), , drop = FALSE]
  free <- free[!duplicated(free$parnumber), , drop = FALSE]
  free[order(free$parnumber), , drop = FALSE]
}

.ctBackendPopMeanSamples <- function(fit, samples = NULL, ...) {
  cells <- .ctBackendFreeParameterCells(fit)
  layout <- .ctBackendSummaryLayout(fit)
  if (is.null(samples)) samples <- .ctBackendRawSamples(fit)
  flat <- .ctBackendParMatricesFlat(fit, t(samples), ...)

  index <- match(cells$matrix, layout$matrix)
  position <- layout$offset[index] + (cells$col - 1L) * layout$nrow[index] + cells$row
  values <- t(flat[position, , drop = FALSE])
  colnames(values) <- ifelse(is.na(cells$param), paste0("param", cells$parnumber),
    as.character(cells$param))
  list(values = values, parnumber = as.integer(cells$parnumber))
}

# Mean / sd / quantiles of a sample matrix, in Stan's summary column order. Uses
# base R rather than rstan's `monitor()`: with one point-estimate "sample" there
# is nothing to monitor, and with a normal-approximation posterior the
# convergence diagnostics `monitor()` adds would be meaningless anyway.
.ctBackendSampleSummary <- function(values, digits = 3) {
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
  round(out, digits)
}

# Raw-scale TI predictor effects, the analogue of Stan's `tipreds` section.
.ctBackendTipredSummary <- function(fit, digits = 3) {
  spec <- .ctBackendSpec(fit)
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) return(NULL)
  samples <- .ctBackendRawSamples(fit)
  cells <- .ctBackendFreeParameterCells(fit)
  parameter_names <- stats::setNames(as.character(cells$param), cells$parnumber)
  predictor_names <- .ctBackendModel(fit)$TIpredNames

  values <- samples[, effects$coefficient, drop = FALSE]
  colnames(values) <- paste0("tip_", predictor_names[effects$predictor], "_",
    ifelse(is.na(parameter_names[as.character(effects$parameter)]),
      paste0("param", effects$parameter),
      parameter_names[as.character(effects$parameter)]))
  .ctBackendSampleSummary(values, digits = digits)
}

.ctBackendSummary <- function(object, timeinterval = 1, digits = 3, parmatrices = TRUE,
  ...) {
  has_posterior <- !is.null(object$estimate$rawposterior)

  out <- list()

  popmeans <- .ctBackendPopMeanSamples(object)
  out$popmeans <- .ctBackendSampleSummary(popmeans$values, digits = digits)
  out$popNote <- paste0("Population values on the transformed scale. ",
    "Covariance parameters appear in sd / unconstrained correlation form; ",
    "see System Matrices (or ctSummaryMatrices()) for cor/cov.")

  tipreds <- .ctBackendTipredSummary(object, digits = digits)
  if (!is.null(tipreds)) {
    out$tipreds <- tipreds
    out$tipredsNote <- "Raw-scale effects on the unconstrained parameters."
  }

  if (isTRUE(parmatrices)) {
    collapsed <- list(
      Mean = .ctBackendSummaryMatrices(object, calcfunc = mean, calcfuncargs = list(),
        timeinterval = timeinterval))
    if (has_posterior) {
      collapsed$sd <- .ctBackendSummaryMatrices(object, calcfunc = stats::sd,
        calcfuncargs = list(na.rm = TRUE), timeinterval = timeinterval)
      for (probability in c(.025, .5, .975)) {
        collapsed[[paste0(probability * 100, "%")]] <- .ctBackendSummaryMatrices(object,
          calcfunc = stats::quantile, calcfuncargs = list(probs = probability),
          timeinterval = timeinterval)
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
    out$parmatrices <- d

    statedep <- attr(ctBackendParMatrices(object), "stateDependent")
    if (!is.null(statedep) && nrow(statedep)) {
      out$parmatNote <- paste0("State-dependent cells (",
        paste0(unique(statedep$matrix), collapse = ", "),
        ") were evaluated at the T0MEANS state and are conditional on it; ",
        "use ctBackendParMatrices(fit, state=) for another.")
    }
  }

  out$loglik <- object$estimate$loglik
  out$npars <- length(object$estimate$raw)
  out$aic <- 2 * out$npars - 2 * out$loglik
  if (has_posterior) out$nsamples <- nrow(object$estimate$rawposterior)
  out$uncertaintyNote <- if (has_posterior) {
    paste0("Julia backend; intervals from ctOptimUncertainty(uncertainty='",
      object$uncertainty$settings$method, "') draws pushed through the transforms.")
  } else {
    paste0("Julia backend; point estimates only. ",
      "Run ctOptimUncertainty() for standard errors and intervals.")
  }

  out <- lapply(out, function(x) roundSummaryCtStanFitValue(x, digits = digits))
  attr(out, "digits") <- digits
  # Two classes, one print method: the sections are named as
  # print.summary.ctStanFit expects, so it prints these unchanged and there is
  # no second printer to keep in step. The extra class is there for code that
  # wants to tell a backend summary apart -- notably, one whose intervals come
  # from a normal approximation rather than from HMC.
  class(out) <- c("summary.ctBackendFit", "summary.ctStanFit")
  out
}
