# What a model says about itself.
#
# `ctModel()` returned a plain list, so printing one dumped 22 fields starting
# with a 26-row data frame of every cell and its transform expression. That is
# the first thing anyone does after building a model, and what it showed was
# `-(1e-06 + 2 * log1p_exp(-(2 * param)))` repeated down a column.
#
# The information that matters is not in any one field. A model's most
# consequential property is usually its *implications*: `indvarying` defaults to
# TRUE for T0MEANS, MANIFESTMEANS and CINT, so most models estimate a population
# covariance their author never asked for and could not see. That was visible in
# the dump as a logical column, and invisible as a fact.

#' @keywords internal
.ctModelFreeNames <- function(pars) {
  free <- is.na(pars$value) & !is.na(pars$param)
  unique(as.character(pars$param[free]))
}

# Names on the axes, because `[,1]` and `[2,]` make the reader count rows to
# find out which latent a cell belongs to -- the same counting problem the
# pipe-separated cell syntax has.
#' @keywords internal
.ctModelPrintMatrix <- function(pars, name, model,
  width = getOption("width", 80L)) {
  rows <- pars$matrix %in% name
  if (!any(rows)) return(invisible(NULL))
  mat <- listOfMatrices(pars[rows, , drop = FALSE])[[name]]
  latent <- model$latentNames
  manifest <- model$manifestNames
  axes <- switch(name,
    T0MEANS = list(latent, ""),
    CINT = list(latent, ""),
    T0VAR = list(latent, latent),
    DRIFT = list(latent, latent),
    DIFFUSION = list(latent, latent),
    LAMBDA = list(manifest, latent),
    MANIFESTMEANS = list(manifest, ""),
    MANIFESTVAR = list(manifest, manifest),
    TDPREDEFFECT = list(latent, model$TDpredNames),
    NULL)
  if (!is.null(axes) && length(axes[[1L]]) == nrow(mat) &&
      (identical(axes[[2L]], "") || length(axes[[2L]]) == ncol(mat))) {
    dimnames(mat) <- list(axes[[1L]],
      if (identical(axes[[2L]], "")) rep("", ncol(mat)) else axes[[2L]])
  }
  cat("\n", name, "\n", sep = "")
  print(mat, quote = FALSE, right = TRUE)
  invisible(NULL)
}

#' Print a ctsem model
#'
#' Summarises what the model is, what it will estimate, and what its
#' individually varying parameters imply, rather than printing the underlying
#' list.
#'
#' @param x A model from \code{\link{ctModel}}.
#' @param matrices Whether to print the model matrices. Defaults to TRUE for
#'   models small enough to read.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @details \code{x$pars} remains the canonical specification and can be
#'   printed directly; \code{\link{ctModelMatrices}} gives the matrix view.
#'
#' @export
print.ctStanModel <- function(x, matrices = NULL, ...) {
  pars <- x$pars
  kind <- if (isTRUE(x$continuoustime)) "continuous" else "discrete"
  cat("ctsem model (", kind, " time)\n", sep = "")
  cat("  ", x$n.latent, " latent (", paste(x$latentNames, collapse = ", "),
    "), ", x$n.manifest, " manifest (", paste(x$manifestNames, collapse = ", "),
    ")\n", sep = "")

  if (isTRUE(x$n.TIpred > 0)) {
    cat("  ", x$n.TIpred, " time independent predictor",
      if (x$n.TIpred > 1L) "s" else "", ": ",
      paste(x$TIpredNames, collapse = ", "), "\n", sep = "")
  }
  if (isTRUE(x$n.TDpred > 0)) {
    cat("  ", x$n.TDpred, " time dependent predictor",
      if (x$n.TDpred > 1L) "s" else "", ": ",
      paste(x$TDpredNames, collapse = ", "), "\n", sep = "")
  }
  if (length(x$groupIDnames)) {
    cat("  grouping above subject: ", paste(x$groupIDnames, collapse = ", "),
      "\n", sep = "")
  }
  if (any(x$manifesttype > 0)) {
    described <- character()
    asym <- .ctAsymptoteRows(x$pars, x$manifesttype) > 0
    binary <- x$manifestNames[x$manifesttype == 1 & !asym]
    if (length(binary)) described <- c(described,
      paste0(paste(binary, collapse = ", "), " (binary)"))
    binary <- x$manifestNames[x$manifesttype == 1 & asym]
    if (length(binary)) described <- c(described,
      paste0(paste(binary, collapse = ", "), " (binary, with ASYMPTOTES)"))
    # Named one at a time, because the category count differs per variable and
    # is the thing a reader most wants to check against their data.
    ordinal <- which(x$manifesttype == 2)
    if (length(ordinal)) described <- c(described,
      paste0(x$manifestNames[ordinal], " (ordinal, ",
        x$ncategories[ordinal], " categories)"))
    counts <- x$manifestNames[x$manifesttype == 3]
    if (length(counts)) described <- c(described,
      paste0(paste(counts, collapse = ", "), " (count, Poisson log link)"))
    # Named one at a time with their limits, which are the thing a reader most
    # wants to check against the range their data actually takes.
    censored <- which(x$manifesttype == 4)
    if (length(censored)) described <- c(described,
      paste0(x$manifestNames[censored], " (censored, ",
        ifelse(is.finite(x$censormin[censored]),
          paste0("min ", signif(x$censormin[censored], 4)), "no min"), ", ",
        ifelse(is.finite(x$censormax[censored]),
          paste0("max ", signif(x$censormax[censored], 4)), "no max"),
        ifelse(censored %in% .ctResponseRows(x$pars), ", with RESPONSE", ""),
        ")"))
    cat("  non-Gaussian indicators: ", paste(described, collapse = ", "),
      "\n", sep = "")
  }

  # Counted as ctFit() will see them: T0VAR cells made redundant by a varying
  # T0MEANS are fixed there, and carry neither a parameter nor an effect. The
  # matrices below still show the cells as written.
  shownpars <- pars
  pars <- suppressMessages(T0VARredundancies(x))$pars
  free <- .ctModelFreeNames(pars)
  cat("\n", length(free), " free parameter", if (length(free) != 1L) "s" else "",
    "\n", sep = "")

  # The implication, stated rather than left to be inferred from a column --
  # for every level of the id hierarchy, each of which estimates its own
  # population covariance.
  levels <- c(x$subjectIDname, x$groupIDnames)
  columns <- .ctVaryingColumns(x)
  varying <- lapply(columns, function(cl) unique(as.character(pars$param[
    .ctVaryingRows(x, cl) & !is.na(pars$param)])))
  names(varying) <- if (length(levels) == length(columns)) levels else columns
  counts <- function(n) paste0(n, " standard deviation", if (n != 1L) "s" else "",
    if (n > 1L) paste0(" and ", n * (n - 1L) / 2L, " correlation",
      if (n * (n - 1L) / 2L != 1L) "s" else "") else "")
  nonempty <- varying[lengths(varying) > 0L]
  if (length(columns) == 1L && length(nonempty)) {
    n <- length(nonempty[[1L]])
    cat("\n", n, " individually varying: ", paste(nonempty[[1L]], collapse = ", "),
      "\n", sep = "")
    cat("  This is a multilevel model. It estimates a ", n, " x ", n,
      " population covariance\n  over those parameters -- ", counts(n),
      " on top of the ", length(free), " above.\n", sep = "")
    cat("  Set indvarying = FALSE in a cell (see ?ctPar) for any that should not vary.\n")
  } else if (length(nonempty)) {
    cat("\nVarying parameters, by level:\n")
    width <- max(nchar(names(varying)))
    for (l in seq_along(varying)) {
      cat("  ", formatC(names(varying)[l], width = -width), "  ",
        if (length(varying[[l]])) paste(varying[[l]], collapse = ", ") else "(none)",
        "\n", sep = "")
    }
    cat("  Each level estimates a population covariance over its own parameters:\n  ",
      paste0(vapply(nonempty, function(v) counts(length(v)), character(1)),
        " (", names(nonempty), ")", collapse = ", "),
      ", on top of the ", length(free), " above.\n", sep = "")
    cat("  Set the levels with indvarying, in a cell or in pars$indvarying: TRUE,",
      " FALSE or level names (see ?ctPar).\n", sep = "")
  } else {
    cat("\nNo individually varying parameters: a fixed effects model.\n")
  }

  # Which parameters the time independent predictors act on, which otherwise
  # takes reading a column per predictor.
  effects <- paste0(x$TIpredNames, "_effect")
  effects <- effects[effects %in% names(pars)]
  if (length(effects)) {
    freerows <- which(is.na(pars$value) & !is.na(pars$param) &
      !duplicated(pars$param))
    acting <- vapply(freerows, function(r) {
      on <- .ctTipredEffectActive(pars[r, effects, drop = FALSE])
      if (any(on)) paste0(pars$param[r], " (",
        paste(sub("_effect$", "", effects[on]), collapse = ", "), ")") else NA_character_
    }, character(1))
    acting <- acting[!is.na(acting)]
    if (length(acting)) cat("\nPredictor effects: ", paste(acting, collapse = ", "),
      "\n", sep = "")
  }

  if (is.null(matrices)) matrices <- x$n.latent <= 4L && x$n.manifest <= 6L
  if (isTRUE(matrices)) {
    for (name in c("T0MEANS", "T0VAR", "DRIFT", "DIFFUSION", "CINT",
      "LAMBDA", "MANIFESTMEANS", "MANIFESTVAR", "TDPREDEFFECT", "PARS")) {
      .ctModelPrintMatrix(shownpars, name, x)
    }
  } else {
    cat("\nMatrices not shown; see model$matrices or ctModelMatrices(model).\n")
  }
  invisible(x)
}
