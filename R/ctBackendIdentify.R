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
    # Negative *against the scale of the matrix*, not against zero. An
    # eigenvalue of -3e-16 where the largest is 9e5 is what a symmetric
    # eigendecomposition does to a direction whose true curvature is zero; it
    # is the flat direction `nweak` already counts, not a saddle. Counting it
    # as negative fired "this is not a maximum" on ordinary fits sitting at
    # their optimum with a gradient of 5e-10, which is how a warning worth
    # reading gets ignored.
    negative = sum(values < -rtol * scale),
    directions = directions,
    parameters = unique(unlist(lapply(directions, `[[`, "parameters"))))
}

# Is each reported interval as wide as the curvature at the estimate supports?
#
# Two standard errors can be computed for a parameter from the same Hessian.
# The *conditional* one, `1 / sqrt(information[i, i])`, is what the curvature in
# that one coordinate supports with every other parameter held. The *marginal*
# one, `sqrt(cov[i, i])`, is what gets reported, and it is the conditional one
# divided by `sqrt(1 - R^2)`, where `R^2` is how well the other parameters
# reproduce this one in the information metric. So their ratio is a pure number
# saying how much of the reported width comes from the data and how much from
# the parameter not being separable from the rest: a ratio of 10 is `R^2` of
# .99, and 100 is .9999.
#
# Why it earns its place. A benchmark fit reached the right optimum -- the same
# log likelihood to eight decimal places as its twin, and Stan's -- and reported
# a drift interval a hundred times too wide, with a point estimate that had
# wandered with it. Nothing else about the fit looked wrong; it was caught only
# because two runs on identical data could be compared, and a user gets one run.
# This ratio separates the two cleanly: 1.0 to 1.5 across every parameter of the
# healthy fits measured here, against 1e4 on the parameter that had gone.
#
# Cheap: one diagonal and one square root, no extra evaluation of anything.
# Computed for julia fits, where the exact Hessian is already on the fit;
# nothing prevents the stan path from using it, and `ctReport()` is where that
# would show.
#' @keywords internal
.ctBackendIntervalCheck <- function(hessian, se, parnames = NULL,
  threshold = 100) {
  empty <- list(threshold = threshold, nflagged = 0L, parameters = character(),
    table = data.frame(param = character(), se = numeric(),
      curvature_se = numeric(), ratio = numeric(), stringsAsFactors = FALSE))
  if (is.null(hessian) || is.null(se)) return(empty)
  hessian <- as.matrix(hessian)
  se <- as.numeric(se)
  if (nrow(hessian) != ncol(hessian) || nrow(hessian) != length(se)) return(empty)
  if (!all(is.finite(hessian))) return(empty)
  n <- length(se)
  if (is.null(parnames) || length(parnames) != n) parnames <- paste0("par", seq_len(n))
  information <- diag(-(hessian + t(hessian)) / 2)
  # A non-positive diagonal is not a wider interval, it is no curvature at all;
  # `.ctBackendIdentifiability()` is what reports that, so it is left NA here
  # rather than counted as a ratio of infinity and reported twice.
  curvature <- ifelse(information > 0, 1 / sqrt(information), NA_real_)
  ratio <- se / curvature
  table <- data.frame(param = as.character(parnames), se = se,
    curvature_se = curvature, ratio = ratio, stringsAsFactors = FALSE)
  flagged <- which(is.finite(ratio) & ratio > threshold)
  list(threshold = threshold, nflagged = length(flagged),
    parameters = as.character(parnames[flagged]),
    table = table[order(-ifelse(is.finite(ratio), ratio, -Inf)), , drop = FALSE])
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
.ctBackendIdentifyWarn <- function(identify, collapsed, intervals = NULL) {
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
  # Said even when no direction is flat enough to count as unidentified, which
  # is the case this exists for: a parameter can be separable in principle and
  # still have almost all of its reported width come from its entanglement with
  # the others, and then the interval moves by orders of magnitude between two
  # runs that reached the same optimum.
  if (!is.null(intervals) && isTRUE(intervals$nflagged > 0L)) {
    involved <- paste(utils::head(intervals$parameters, 6), collapse = ", ")
    if (length(intervals$parameters) > 6) involved <- paste0(involved, ", ...")
    widest <- max(intervals$table$ratio[is.finite(intervals$table$ratio)])
    warning(intervals$nflagged, " reported interval",
      if (intervals$nflagged > 1L) "s are" else " is",
      " far wider than the curvature at the estimate supports -- up to ",
      signif(widest, 3), " times the width that parameter's own curvature ",
      "gives. That width comes from the parameter not being separable from ",
      "the others rather than from the data, and it is not stable: it can ",
      "move by orders of magnitude between two fits that reach the same ",
      "optimum. Parameters involved: ", involved,
      ". See fit$uncertainty$intervalcheck.", call. = FALSE)
  }
  invisible(NULL)
}
