# How far the reported answer is from a better one.
#
# ctsem asks this three times, against three different notions of "better", and
# until now answered it in three shapes under three names:
#
#   reference                                    measure              the number
#   -------------------------------------------  -------------------  ----------
#   the exact Hessian, against the optimiser's    .ctBackendOptimGap   $gap
#     implicit curvature
#   adaptive quadrature, against the Laplace      ctLaplaceCheck       $gap
#     approximation
#   a particle filter, against the assumed-       ctParticleLik        $difference
#     density filter
#
# All three are the same quantity: a log likelihood difference between what the
# fit reports and what a more accurate computation says. Two call it `gap` and
# one calls it `difference`; one normalises per subject and two do not; one
# states a tolerance and two leave the reader to decide. So a user asking "how
# good is this approximation" gets three answers of three shapes, and a GUI
# asking it has to know which function it called.
#
# `.ctFitGap()` is the one shape. It is attached *alongside* each function's
# existing return rather than replacing it, because two of the three are
# exported with a documented structure that scripts read.
#
# Deliberately does not invent a verdict where the package never had one. Only
# the curvature route states a tolerance (0.01 log likelihood, see
# `.ctBackendCertify`); for the other two `tolerance` is NA and
# `exceeds_tolerance` is NA with it. Reporting "immaterial" against a bar nobody
# chose would be worse than reporting the number and saying no bar is set.

#' A comparison between a fit's own answer and a more accurate one
#'
#' The common shape behind [ctLaplaceCheck()], [ctParticleLik()] and the
#' curvature certification on `fit$uncertainty$certification`.
#'
#' @param reference short name for what the fit was compared against:
#'   `"curvature"`, `"quadrature"` or `"particle"`.
#' @param gap log likelihood difference, better minus reported. Positive means
#'   the better computation finds more likelihood than the fit reports.
#'
#'   One asymmetry between the three, which is a property of the references and
#'   not of this shape: the curvature gap is `1/2 g' H^-1 g` and so is
#'   non-negative by construction, while the quadrature and particle gaps are
#'   signed. A negative particle gap says the assumed-density filter reports
#'   *more* likelihood than the particle filter finds, which is a finding rather
#'   than an error -- hence `resolved` tests `abs(gap)`, not `gap`.
#' @param tolerance the bar, in log likelihood units, or `NA` where the route
#'   states none.
#' @param gap_se standard error of `gap`, where the reference is stochastic.
#' @param nsubjects used to report the gap per subject, which is what says
#'   whether a total is large.
#' @param remedy one sentence on what to do if the gap matters.
#' @param detail optional per-parameter or per-row table the route already
#'   produces.
#'
#' @return an object of class `ctFitGap`.
#' @keywords internal
.ctFitGap <- function(reference, gap, tolerance = NA_real_, gap_se = NA_real_,
  nsubjects = NA_integer_, remedy = NULL, detail = NULL) {

  gap <- suppressWarnings(as.numeric(gap)[1L])
  gap_se <- suppressWarnings(as.numeric(gap_se)[1L])
  tolerance <- suppressWarnings(as.numeric(tolerance)[1L])
  nsubjects <- suppressWarnings(as.integer(nsubjects)[1L])

  # Two separate questions, kept separate. "Is it bigger than the bar" and "is
  # it bigger than its own noise" are different, and collapsing them into one
  # `material` would hide which of the two failed -- a gap of 3 with a standard
  # error of 4 and a gap of 0.001 with no error at all are both "not material"
  # and want different responses.
  exceeds <- if (is.finite(gap) && is.finite(tolerance)) gap > tolerance else NA
  resolved <- if (is.finite(gap) && is.finite(gap_se) && gap_se > 0)
    abs(gap) > 2 * gap_se else NA

  structure(list(
    reference = as.character(reference)[1L],
    gap = gap,
    gap_se = gap_se,
    gap_per_subject = if (is.finite(gap) && !is.na(nsubjects) && nsubjects > 0L)
      gap / nsubjects else NA_real_,
    nsubjects = nsubjects,
    tolerance = tolerance,
    exceeds_tolerance = exceeds,
    resolved = resolved,
    remedy = remedy,
    detail = detail), class = "ctFitGap")
}

.CT_GAP_REFERENCE <- c(
  curvature = "the exact curvature at the estimate",
  quadrature = "adaptive quadrature over the random effects",
  particle = "a bootstrap particle filter")

#' @export
print.ctFitGap <- function(x, ...) {
  what <- .CT_GAP_REFERENCE[[x$reference]]
  if (is.null(what) || !length(what)) what <- x$reference
  cat("Comparison against ", what, "\n", sep = "")
  if (!is.finite(x$gap)) {
    cat("  gap: not available\n")
    return(invisible(x))
  }
  cat("  gap: ", format(x$gap, digits = 4), " log likelihood",
    if (is.finite(x$gap_se)) paste0(" (se ", format(x$gap_se, digits = 3), ")") else "",
    "\n", sep = "")
  if (is.finite(x$gap_per_subject)) {
    cat("       ", format(x$gap_per_subject, digits = 3), " per subject over ",
      x$nsubjects, "\n", sep = "")
  }
  if (is.finite(x$tolerance)) {
    cat("  bar: ", signif(x$tolerance, 3), " -- ",
      if (isTRUE(x$exceeds_tolerance)) "exceeded" else "not exceeded", "\n", sep = "")
  } else {
    # Said rather than left blank. A reader who sees no verdict should know it
    # is because none is defined here, not because the gap passed something.
    cat("  bar: none stated for this comparison; the number above is the finding\n")
  }
  if (isFALSE(x$resolved)) {
    cat("  the gap is not distinguishable from zero at twice its standard error\n")
  }
  if (!is.null(x$remedy)) cat("  ", x$remedy, "\n", sep = "")
  invisible(x)
}
