# Time independent predictor effects for backend='julia' fits -----------------
#
# ctTIpredEffects() (R/ctStanTIpredeffects.R) reads stan-only structures
# (fit$stanfit$rawest, fit$setup$popsetup, stan_constrainsamples()) to build a
# grid of TI-predictor values, add each predictor's raw effect to the affected
# parameters, push the result through the model's transforms, and collapse to
# quantiles. Every one of those steps has a julia counterpart already built for
# ctBackendParMatrices()/ctSummaryMatrices()/summary() (R/ctBackendSummary.R),
# so this reuses them rather than re-deriving the transform:
#
#   fit$estimate$raw / $rawposterior   stan's stanfit$rawest / $rawposterior
#   spec$ti_effects                    which raw parameter each TI coefficient
#                                       displaces, and where the coefficient
#                                       itself lives in the raw vector
#   spec$tipred_data                   stan's fit$data$tipredsdata
#   .ctBackendParMatricesFlat(tipreds=) stan_constrainsamples(), for one grid
#                                       point across the whole sample at once
#   .ctSummaryMatricesFromArrays()     the same matrix-family collapse
#                                       summary() and ctSummaryMatrices() use
#
# One engine call per grid point (not per grid point per sample, and not per
# grid point per quantile) carries the whole sample matrix in one transfer, so
# the julia round-trip cost here is nsubjects calls, not
# nsubjects*nsamples*length(probs).

# Raw parameter cells affected by the requested predictor(s), and the default
# whichpars set when none was named: every parameter any selected predictor has
# a coefficient for. Matches the stan path's
# `which(apply(tieffect,2,function(x) any(x!=0)))`, read here from the model's
# own ti_effects table rather than from a posterior draw.
.ctBackendTIpredWantedPars <- function(effects, whichTIpreds, whichpars) {
  if (identical(whichpars, "all")) {
    wanted <- sort(unique(effects$parameter[effects$predictor %in% whichTIpreds]))
    if (!length(wanted)) {
      stop("None of the selected TI predictor(s) (whichTIpreds) have an ",
        "effect on any parameter.", call. = FALSE)
    }
    return(wanted)
  }
  as.integer(whichpars)
}

# Every model matrix, quantile-collapsed, at one TI-predictor grid point.
#
# `.ctBackendPopArrays(tipreds=)` materialises the pop_* arrays for the whole
# sample at this one covariate vector -- the engine call -- and
# `.ctSummaryMatricesFromArrays()` is then called once per requested quantile,
# which is cheap (an R-side collapse of an already-materialised array).
.ctBackendTIpredEffectsMatrices <- function(fit, samples, fullpreds, probs,
  whichpars, timeinterval) {

  model <- .ctBackendModel(fit)
  nsub <- nrow(fullpreds)
  nprobs <- length(probs)
  selection <- NULL
  parnames <- NULL
  out <- NULL

  # Off-diagonal cells of a symmetric covariance are the same value twice; the
  # stan path drops the redundant (lower-triangle) half unless a caller asked
  # for a specific bracketed cell, and the corresponding julia-side matrix
  # family (.ctSummaryMatricesFromArrays) already omits MANIFESTVAR entirely,
  # so that half of the stan rule has nothing left to apply to.
  symmetric <- c("MANIFESTcov", "T0cov", "DIFFUSIONcov", "dtDIFFUSIONcov",
    "asymDIFFUSIONcov")
  bracketed <- !identical(whichpars, "all") && any(grepl("[", whichpars, fixed = TRUE))

  for (i in seq_len(nsub)) {
    arrays_i <- .ctBackendPopArrays(fit, samples = samples, tipreds = fullpreds[i, ])
    for (j in seq_len(nprobs)) {
      matlist <- .ctSummaryMatricesFromArrays(arrays_i,
        continuoustime = model$continuoustime, latentNames = model$latentNames,
        manifestNames = model$manifestNames, TDpredNames = model$TDpredNames,
        calcfunc = stats::quantile, calcfuncargs = list(probs = probs[j]),
        timeinterval = timeinterval)
      flat <- ctModelUnlist(matlist, matnames = names(matlist))

      if (is.null(selection)) {
        keep <- rep(TRUE, nrow(flat))
        if (!bracketed) keep <- keep & !(flat$matrix %in% symmetric & flat$row < flat$col)
        if (!identical(whichpars, "all")) {
          cellname <- paste0(flat$matrix, "[", flat$row, ",", flat$col, "]")
          matched <- lapply(whichpars, function(x) which(flat$matrix %in% x | cellname %in% x))
          keep <- keep & seq_len(nrow(flat)) %in% unlist(matched)
        }
        selection <- which(keep)
        if (!length(selection)) {
          stop("whichpars did not match any reported parameter matrix or cell.",
            call. = FALSE)
        }
        parnames <- paste0(flat$matrix[selection], "[", flat$row[selection], ",",
          flat$col[selection], "]")
        out <- array(NA_real_, dim = c(nprobs, length(parnames), nsub))
      }
      out[j, , i] <- flat$value[selection]
    }
  }
  dimnames(out) <- list(Quantile = paste0("Quantile", probs), param = parnames,
    subject = NULL)
  out
}

# Individual raw parameters' own transformed value, at one TI-predictor grid
# point -- the parmatrices=FALSE path. `.ctBackendFreeParameterCells()` maps a
# raw parameter number onto the single model-matrix cell it fills, so reading
# that cell through the same engine assembly parmatrices=TRUE uses gives
# exactly the raw parameter's own transform: there is no separate reduction
# (Cholesky, quadrature) between a free parameter and its cell.
#
# Not filtered to non-random-effect cells: an individually varying parameter
# is exactly the usual target of a TI-predictor effect (a covariate shifting
# a random effect's population mean), and `.ctBackendFreeParameterCells()`
# already resolves such a parameter's carrier row to the matrix cell that
# actually reads it, so its population-mean transform is reported correctly
# either way.
.ctBackendTIpredEffectsPars <- function(fit, samples, fullpreds, effects,
  whichTIpreds, whichpars, probs, returndifference) {

  cells <- .ctBackendFreeParameterCells(fit)
  wantedpar <- .ctBackendTIpredWantedPars(effects, whichTIpreds, whichpars)
  wantedcells <- cells[match(wantedpar, cells$parnumber), , drop = FALSE]
  if (anyNA(wantedcells$parnumber)) {
    stop("whichpars references a parameter number this model does not have; ",
      "see fit$model_spec$parameter_table$parnumber.", call. = FALSE)
  }
  layout <- .ctBackendSummaryLayout(fit)
  parnames <- .ctBackendParameterNames(wantedcells)

  nsub <- nrow(fullpreds)
  nprobs <- length(probs)
  out <- array(NA_real_, dim = c(nprobs, length(parnames), nsub))

  # The stan path's "noeffect" baseline: the same parameter, transformed with
  # no TI-predictor contribution added, so the difference isolates the
  # covariate's effect from the population value it is added to.
  baseline <- if (isTRUE(returndifference)) {
    .ctBackendPopCellValues(fit, samples, wantedcells, layout)
  } else NULL

  for (i in seq_len(nsub)) {
    values <- .ctBackendPopCellValues(fit, samples, wantedcells, layout,
      tipreds = fullpreds[i, ])
    if (!is.null(baseline)) values <- values - baseline
    for (j in seq_len(nprobs)) {
      out[j, , i] <- apply(values, 2, stats::quantile, probs = probs[j], na.rm = TRUE)
    }
  }
  dimnames(out) <- list(Quantile = paste0("Quantile", probs), param = parnames,
    subject = NULL)
  out
}

# The julia counterpart of ctTIpredEffects() (R/ctStanTIpredeffects.R). Same
# arguments, same return shape (list(y=, x=)); see that function's roxygen for
# the user-facing documentation, and the header comment above for what each
# piece is read from.
.ctBackendTIpredEffects <- function(fit, returndifference = FALSE,
  probs = c(.025, .5, .975), includeMeanUncertainty = FALSE, whichTIpreds = 1,
  parmatrices = TRUE, whichpars = "all", nsamples = 100, timeinterval = 1,
  nsubjects = 20, filter = NA, plot = FALSE) {

  model <- .ctBackendModel(fit)
  spec <- .ctBackendSpec(fit)

  # Named rather than caught downstream as a bounds error: the model having no
  # TI predictors at all is the ordinary case for someone exploring, and
  # deserves the same message the stan path gives.
  if (is.null(model$n.TIpred) || model$n.TIpred < 1) {
    stop("This model has no time independent predictors, so there are no ",
      "effects to report. Add them with n.TIpred and TIpredNames in ctModel().",
      call. = FALSE)
  }
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) {
    stop("This model has no time independent predictors, so there are no ",
      "effects to report. Add them with n.TIpred and TIpredNames in ctModel().",
      call. = FALSE)
  }
  ntipred <- ncol(spec$tipred_data)
  # A smooth term is one predictor here: named by its covariate, or reached
  # through any of its columns, it is swept over the covariate and its columns
  # evaluated from the stored basis at each point (R/ctTIpredSmooth.R).
  smooth <- .ctTIsmoothSelected(model, whichTIpreds)
  if (is.character(whichTIpreds)) {
    whichTIpreds <- if (!is.null(smooth)) match(smooth$columns, model$TIpredNames) else
      match(whichTIpreds, model$TIpredNames)
    if (anyNA(whichTIpreds)) stop("whichTIpreds names a predictor or smooth ",
      "covariate this model does not have.", call. = FALSE)
  } else if (!is.null(smooth)) {
    whichTIpreds <- match(smooth$columns, model$TIpredNames)
  }
  if (max(whichTIpreds) > ntipred || min(whichTIpreds) < 1) {
    stop("whichTIpreds selects a predictor outside 1:", ntipred,
      ", the predictors this model has.", call. = FALSE)
  }

  # Samples of the raw parameter vector: the posterior draws from
  # ctOptimUncertainty()/ctFitUncertainty(fit, 'sample') when present,
  # otherwise the point estimate repeated -- the same fallback
  # ctBackendParMatrices() and friends use.
  samples <- .ctBackendRawSamples(fit)
  niter <- nrow(samples)
  if (identical(nsamples, "all") || nsamples > niter) nsamples <- niter
  samplerows <- sample.int(niter, nsamples)

  samples <- samples[samplerows, , drop = FALSE]
  # Without the mean parameters' uncertainty, they are held at their median
  # and the TI-predictor effects keep theirs, as documented. Holding every
  # coordinate at the median, coefficients included, left no uncertainty at all.
  if (!isTRUE(includeMeanUncertainty)) {
    held <- setdiff(seq_len(ncol(samples)), as.integer(effects$coefficient))
    samples[, held] <- matrix(apply(samples[, held, drop = FALSE], 2, stats::median),
      nrow = nsamples, ncol = length(held), byrow = TRUE)
  }

  tipreds <- spec$tipred_data
  if (any(!is.na(filter))) {
    tipreds <- eval(parse(text = paste0("tipreds[tipreds[,", filter[1], "]",
      filter[2], ",,drop=FALSE]")))
  }
  tipreds <- tipreds[, whichTIpreds, drop = FALSE]
  if (identical(nsubjects, "all")) nsubjects <- nrow(tipreds)
  # If real subject-level tipred data will be used (more than one predictor,
  # so an interaction can't be swept as a smooth grid), there's no point using
  # more subjects than the data has.
  if (nsubjects > nrow(tipreds) && length(whichTIpreds) > 1) nsubjects <- nrow(tipreds)

  covariate <- NULL
  if (!is.null(smooth)) {
    # The covariate's own values, one per person, are the x axis.
    persons <- spec$data[[smooth$covariate]][!duplicated(spec$data[[model$subjectIDname]])]
    covariate <- seq(min(persons, na.rm = TRUE), max(persons, na.rm = TRUE),
      length.out = nsubjects)
    tipreds <- .ctTIsmoothEval(smooth$basis, covariate)
  } else if (length(whichTIpreds) > 1) {
    tipreds <- tipreds[sample.int(nrow(tipreds), nsubjects), , drop = FALSE]
  }
  if (is.null(smooth) && length(whichTIpreds) == 1) {
    tipreds <- cbind(seq(from = min(tipreds), to = max(tipreds, na.rm = TRUE),
      length.out = nsubjects))
  }
  tiorder <- if (is.null(smooth)) order(tipreds[, 1]) else seq_len(nrow(tipreds))
  tipreds <- tipreds[tiorder, , drop = FALSE]

  # Every predictor not selected stays at zero -- the same as the stan path
  # never adding its effect in.
  fullpreds <- matrix(0, nrow = nrow(tipreds), ncol = ntipred)
  fullpreds[, whichTIpreds] <- tipreds

  message(sprintf("Getting %s samples by %s subjects for %s total samples",
    nsamples, nrow(tipreds), nsamples * nrow(tipreds)))
  message("Calculating time independent predictor effects...")

  if (isTRUE(parmatrices)) {
    out <- .ctBackendTIpredEffectsMatrices(fit, samples = samples,
      fullpreds = fullpreds, probs = probs, whichpars = whichpars,
      timeinterval = timeinterval)
    # Population matrices are evaluated at the T0MEANS state with TD
    # predictors at zero, as every other population-matrix summary is; said
    # once here, as the stan path does for this same function.
    .ctContextMessage(fit, .ctContextPopLabel)
  } else {
    out <- .ctBackendTIpredEffectsPars(fit, samples = samples,
      fullpreds = fullpreds, effects = effects, whichTIpreds = whichTIpreds,
      whichpars = whichpars, probs = probs, returndifference = returndifference)
  }
  dimnames(out)$subject <- tiorder

  colnames(tipreds) <- colnames(spec$tipred_data)[whichTIpreds]
  xname <- if (is.null(smooth)) paste(colnames(tipreds), collapse = "") else
    smooth$covariate
  names(attributes(out)$dimnames) <- c("Parameter", "param", xname)
  out <- list(y = aperm(out, c(3, 2, 1)), x = if (is.null(smooth))
    tipreds[, 1, drop = FALSE] else matrix(covariate, ncol = 1L,
      dimnames = list(NULL, smooth$covariate)))
  # For a smooth, each person's own value of what is plotted, at their own
  # covariate value: the curve plus their random effects, so the spread about
  # the curve is the person differences it leaves (`$persons`).
  if (!is.null(smooth)) {
    out$persons <- .ctTIsmoothPersons(fit, smooth, dimnames(out$y)[[2L]],
      parmatrices = parmatrices)
    attr(out, "smooth") <- .ctTIsmoothSubtitle(fit, smooth)
  }

  if (!isTRUE(plot)) return(out)
  g <- ctPlotArrayGG(out)
  if (!is.null(out$persons) && nrow(out$persons)) {
    x <- value <- variable <- NULL
    g <- g + ggplot2::geom_point(data = out$persons,
        ggplot2::aes(x = x, y = value, colour = variable), alpha = .35,
        size = .9, inherit.aes = FALSE) +
      ggplot2::geom_rug(data = out$persons[!duplicated(out$persons$x), ],
        ggplot2::aes(x = x), sides = "b", alpha = .3, inherit.aes = FALSE) +
      ggplot2::labs(x = smooth$covariate, subtitle = attr(out, "smooth"))
  }
  g
}
