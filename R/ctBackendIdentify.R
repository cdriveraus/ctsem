# What a fit can honestly say about its own credibility.
#
# An optimiser can stop somewhere the data does not determine, and the fit then
# reports estimates and standard errors that look exactly like any other fit's.
# Observed on a six-subject model: population standard deviations driven to
# 1e-8, process parameters at plus or minus a hundred, `converged = TRUE`.
#
# The temptation is to call such a fit absurd and refuse it. That is the wrong
# move, because "absurd" is a judgement about the model and the data that this
# code is in no position to make -- a population standard deviation of zero is a
# legitimate finding (no detectable individual differences), an extreme raw
# value can be a perfectly ordinary value on a transformed scale, and a
# likelihood that is flat in some direction may still be exactly what the
# researcher wants to report.
#
# What can be said without assuming anything is narrower and more useful: which
# directions the data does not identify, and therefore which reported intervals
# do not mean what a reader will take them to mean. That is a property of the
# information matrix, checkable, and it does not require deciding whether the
# answer is sensible.

#' @keywords internal
.ctBackendIdentifiability <- function(hessian, parnames = NULL, rtol = 1e-8,
  loading = 0.25) {
  empty <- list(nweak = 0L, condition = NA_real_, directions = list(),
    parameters = character())
  if (is.null(hessian)) return(empty)
  hessian <- as.matrix(hessian)
  if (!all(is.finite(hessian)) || nrow(hessian) != ncol(hessian)) return(empty)
  n <- nrow(hessian)
  if (is.null(parnames) || length(parnames) != n) parnames <- paste0("par", seq_len(n))

  # The information matrix, symmetrised. Eigenvalues rather than a determinant
  # or a condition number alone: which direction is flat is the actionable part,
  # and only the eigenvectors say that.
  information <- -(hessian + t(hessian)) / 2
  decomposition <- try(eigen(information, symmetric = TRUE), silent = TRUE)
  if (inherits(decomposition, "try-error")) return(empty)
  values <- decomposition$values
  scale <- max(abs(values))
  if (!is.finite(scale) || scale <= 0) return(empty)

  # A direction counts as unidentified when its curvature is negligible against
  # the sharpest direction, or when it is negative -- a negative eigenvalue of
  # the information matrix means the optimiser stopped somewhere that is not a
  # maximum in that direction at all.
  weak <- which(values <= rtol * scale)
  directions <- lapply(weak, function(k) {
    loadings <- decomposition$vectors[, k]
    involved <- order(-abs(loadings))
    involved <- involved[abs(loadings[involved]) >= loading]
    if (!length(involved)) involved <- which.max(abs(loadings))
    list(eigenvalue = values[k], relative = values[k] / scale,
      parameters = parnames[involved], loadings = loadings[involved])
  })
  list(
    nweak = length(weak),
    condition = scale / max(min(values[values > 0], na.rm = TRUE), .Machine$double.xmin),
    negative = sum(values < 0),
    directions = directions,
    parameters = unique(unlist(lapply(directions, `[[`, "parameters"))))
}

# Population standard deviations that have collapsed to the floor of their
# transform.
#
# Reported as a finding rather than a fault. A zero population standard
# deviation says the data show no individual differences in that parameter,
# which is a result; what it also says is that the parameter sits at the edge of
# its own transform, where the curvature is zero and the reported interval is
# therefore not a confidence statement about anything.
#' @keywords internal
.ctBackendCollapsedScales <- function(fit, tolerance = 1e-6) {
  spec <- fit$model_spec
  if (is.null(spec$laplace)) return(data.frame())
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  estimate <- .ctJuliaNumericVector(as.numeric(fit$estimate$raw))
  levels <- max(1L, length(spec$laplace$levels))
  rows <- list()
  for (l in seq_len(levels)) {
    covariance <- try(.ctBackendJuliaValue(module$ctsem_laplace_popcov(
      objective, estimate, as.integer(l))), silent = TRUE)
    if (inherits(covariance, "try-error")) next
    covariance <- as.matrix(covariance)
    sds <- sqrt(abs(diag(covariance)))
    name <- if (!is.null(spec$laplace$levels) &&
        !is.null(spec$laplace$levels[[l]]$name)) {
      as.character(spec$laplace$levels[[l]]$name)
    } else as.character(l)
    collapsed <- which(sds <= tolerance)
    if (length(collapsed)) {
      rows[[length(rows) + 1L]] <- data.frame(level = name,
        effect = collapsed, sd = sds[collapsed], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Say it once, at the end of a fit, in the terms a reader needs.
#' @keywords internal
.ctBackendIdentifyWarn <- function(identify, collapsed) {
  if (!is.null(identify) && identify$nweak > 0L) {
    involved <- paste(utils::head(identify$parameters, 6), collapse = ", ")
    if (length(identify$parameters) > 6) involved <- paste0(involved, ", ...")
    warning("The data do not identify ", identify$nweak, " direction",
      if (identify$nweak > 1L) "s" else "", " of this model. The estimates ",
      "are still whatever the optimiser found, but the standard errors along ",
      "those directions are arbitrary rather than small or large, and any ",
      "interval built from them will be too. Parameters involved: ", involved,
      ". See fit$identifiability, and ctIdentify(data, model) to check this ",
      "before spending a fit next time.", call. = FALSE)
  }
  if (!is.null(identify) && isTRUE(identify$negative > 0L)) {
    warning(identify$negative, " direction",
      if (identify$negative > 1L) "s have" else " has",
      " negative curvature at the estimate, so it is not a maximum there. ",
      "Treat the estimate as a stopping point rather than a solution.",
      call. = FALSE)
  }
  if (is.data.frame(collapsed) && nrow(collapsed)) {
    message(nrow(collapsed), " population standard deviation",
      if (nrow(collapsed) > 1L) "s were" else " was",
      " estimated at zero (level ",
      paste(unique(collapsed$level), collapse = ", "),
      "). That is a finding -- no detectable individual differences -- but the ",
      "parameter sits at the edge of its transform, where the curvature is ",
      "zero, so its reported interval is not a confidence statement. See ",
      "fit$collapsedScales.")
  }
  invisible(NULL)
}
