# Sampling weights (`ctFit(weights = )`).
#
# A weight multiplies its subject's contribution to the objective: on the
# filter's route the subject's log marginal likelihood, on the Laplace route
# its unit's term. That is the weighted pseudo-likelihood of survey analysis,
# with the weight at the level the data are sampled at -- the person -- and
# never on the prior. The engine applies it where it sums subjects
# (`subject_weights` on the objective), so the value, the gradient, the
# Hessian and the per-subject scores all carry it; the per-subject log
# likelihoods it reports stay the data's own.
#
# The weights are normalised to a mean of one over the subjects that have one,
# so the data keep the weight the prior is balanced against, and a fit's
# estimates do not depend on the weights' scale. Standard errors come from the
# sandwich estimator unless another uncertainty is asked for: the Hessian alone
# describes the pseudo-likelihood's curvature, not the estimate's variability
# across samples.
#
# Not yet: sampling (a weighted likelihood is not a posterior), the state-
# explicit route, the stan backend, and weights that differ within a group of
# an outer level of nested ids (those need weights inside the group's
# integral).

# One weight per subject in the specification's order (`.ctFitIdMap()`),
# normalised; NULL for none.
.ctResolveSubjectWeights <- function(spec, weights, datalong) {
  if (is.null(weights)) return(NULL)
  idname <- spec$model$subjectIDname
  ids <- unique(spec$data[[idname]])
  key <- as.character(ids)
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(datalong)) stop("weights = '", weights, "' is not a ",
      "column of the data.", call. = FALSE)
    values <- datalong[[weights]]
    subject <- as.character(datalong[[idname]])
    per <- tapply(values, factor(subject, levels = key), function(v) {
      v <- unique(v[!is.na(v)])
      if (length(v) == 1L) v else NA_real_
    })
    w <- as.numeric(per[key])
    bad <- key[is.na(w)]
    if (length(bad)) stop("weights must take one value per subject, and the ",
      "column '", weights, "' is missing or varies within subject(s) ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) ", ...", ".", call. = FALSE)
  } else if (is.numeric(weights)) {
    if (!is.null(names(weights))) {
      w <- as.numeric(weights[key])
      if (anyNA(w)) stop("weights is named, but not for every subject id.",
        call. = FALSE)
    } else if (length(weights) == length(key)) {
      w <- as.numeric(weights)
    } else stop("weights needs one value per subject (", length(key),
      "), named by id or in the order subjects first appear in the data, or ",
      "the name of a column of the data.", call. = FALSE)
  } else stop("weights must be the name of a column of the data, or one ",
    "number per subject.", call. = FALSE)
  if (any(!is.finite(w)) || any(w < 0) || !any(w > 0)) stop("weights must be ",
    "finite and non-negative, and not all zero.", call. = FALSE)
  w <- w / mean(w[w > 0])
  .ctWeightsNestedCheck(spec, w)
  w
}

# On the Laplace route a unit is a group of the outermost level, integrated as
# a whole, so its subjects' weights have to agree.
.ctWeightsNestedCheck <- function(spec, w) {
  levels <- spec$laplace$levels
  if (is.null(levels) || length(levels) < 2L) return(invisible(NULL))
  outer <- levels[[length(levels)]]$name
  idname <- spec$model$subjectIDname
  data <- spec$data
  ids <- as.character(unique(data[[idname]]))
  group <- as.character(data[[outer]][match(ids, as.character(data[[idname]]))])
  varies <- tapply(w, group, function(x) diff(range(x)) > 1e-12 * max(abs(x)))
  if (any(varies)) stop("With nested ids, weights must be one value within each ",
    "'", outer, "' group: that group's random effects are integrated as a whole, ",
    "and weights within it would have to enter that integral. Weight the '",
    outer, "' groups, or fit without nesting.", call. = FALSE)
  invisible(NULL)
}

# Whether a fit was weighted.
.ctFitWeighted <- function(fit) {
  spec <- tryCatch(.ctBackendSpec(fit), error = function(e) NULL)
  length(spec$subject_weights) > 0L
}

# Functions not yet defined for a weighted fit refuse by name rather than
# treating it as unweighted.
.ctWeightsRefuse <- function(fit, what) {
  if (.ctFitWeighted(fit)) stop(what, " is not available for a fit with sampling ",
    "weights yet: it would treat the subjects as unweighted.", call. = FALSE)
  invisible(NULL)
}
