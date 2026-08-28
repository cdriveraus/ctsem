# The population covariance, as a matrix of its own.
#
# ## Why a separate matrix rather than a bigger T0VAR
#
# Under the augmented layout a varying parameter becomes a state, so its
# population covariance really does end up inside T0VAR, and extending T0VAR to
# expose it is the obvious move. It does not work, and the reason is worth
# keeping: **the population covariance is not an ordinary free parameter.** Stan
# parameterises it through `rawpopcovbase`/`rawpopsd`, the julia backend builds
# its own `julia_popcov_*` entries, and adding labelled T0VAR cells counts it a
# second time -- one varying parameter took a one-latent model from 4 free
# parameters to 7 and the Stan round trip failed with "Number of unconstrained
# parameters does not match that of the model (8 vs 23)".
#
# Worse, the cells could not be told apart from the placeholders the augmentation
# already writes there, because `ctStanModelIntOverPop()` creates them with
# `value = 0` and a user fixing a covariance to zero writes exactly the same
# thing.
#
# A matrix of its own solves both at once: its *name* is the marker. Nothing that
# walks `pars` looking for free parameters can see it, because it is not in
# `pars` at all -- it is a field on the model, surfaced through the matrix view.
# That is what makes this safe where the T0VAR route was not, and it is why the
# cells here are the specification rather than placeholders standing in for one.
#
# ## What it holds
#
# One row and column per individually varying parameter, named for that
# parameter. Written like any other ctsem matrix: a number fixes a cell, a
# character labels a free one.
#
#   m$matrices$POPCOV['mm', 'mm'] <- 0.3     # this spread, exactly
#   m$matrices$POPCOV['mm', 'T0m_eta1'] <- 0 # uncorrelated with that effect
#
# Diagonal entries are standard deviations and off-diagonal entries below it are
# correlations, which is the form both backends already work in (`popsd` and
# `rawpopcorr` in every summary) rather than a covariance needing decomposition.
# Above the diagonal is a fixed zero, as T0VAR's own upper triangle is.

#' @keywords internal
.ctModelPopCovNames <- function(pars) {
  if (is.null(pars$indvarying)) return(character())
  varying <- !is.na(pars$indvarying) & pars$indvarying & !is.na(pars$param)
  if (!any(varying)) return(character())
  unique(as.character(pars$param[varying]))
}

# The default specification: everything free, labelled after the parameter it
# belongs to.
#
# Free rather than fixed because that is what the model already does -- a varying
# parameter's population spread has always been estimated -- and the point of
# surfacing it is to show what the model implies, not to change it.
#' @keywords internal
.ctModelPopCov <- function(pars) {
  names <- .ctModelPopCovNames(pars)
  if (!length(names)) return(NULL)
  n <- length(names)
  out <- matrix("0", n, n, dimnames = list(names, names))
  for (i in seq_len(n)) {
    out[i, i] <- paste0("popsd_", names[i])
    if (i > 1L) for (j in seq_len(i - 1L)) {
      out[i, j] <- paste0("popcorr_", names[i], "__", names[j])
    }
  }
  out
}

# Keep a POPCOV in step with the parameters that are varying now.
#
# `indvarying` can change after `ctModel()` returns -- setting it directly is
# how nearly every multilevel model in the tests is written -- so the matrix has
# to be rebuilt on demand rather than only at construction. Cells the user set
# are carried across by *name*, not position, so adding a random effect does not
# shuffle the specification of the ones already there.
#' @keywords internal
.ctModelPopCovSync <- function(model) {
  fresh <- .ctModelPopCov(model$pars)
  previous <- model[["POPCOV"]]
  if (is.null(fresh)) {
    model[["POPCOV"]] <- NULL
    return(model)
  }
  if (!is.null(previous) && !is.null(dimnames(previous))) {
    shared_row <- intersect(rownames(fresh), rownames(previous))
    shared_col <- intersect(colnames(fresh), colnames(previous))
    if (length(shared_row) && length(shared_col)) {
      fresh[shared_row, shared_col] <- previous[shared_row, shared_col]
    }
  }
  model[["POPCOV"]] <- fresh
  model
}

# What a cell says: a number, or a label for something to estimate.
#' @keywords internal
.ctModelPopCovValue <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

#' @keywords internal
.ctModelPopCovFree <- function(x) {
  is.na(.ctModelPopCovValue(x)) & !is.na(x) & nzchar(as.character(x))
}

# The entry for one varying parameter pair, by name, or NA when the model has
# nothing to say about it.
#' @keywords internal
.ctModelPopCovEntry <- function(model, rowname, colname = rowname) {
  popcov <- model[["POPCOV"]]
  if (is.null(popcov)) return(NA_character_)
  if (!rowname %in% rownames(popcov) || !colname %in% colnames(popcov)) {
    return(NA_character_)
  }
  as.character(popcov[rowname, colname])
}

# Take a user's POPCOV, checked against what the model actually has.
#
# By name rather than position: a matrix written for one set of random effects
# and assigned to a model with another would otherwise silently attach each
# value to the wrong parameter.
#' @keywords internal
.ctModelPopCovAssign <- function(model, value) {
  model <- .ctModelPopCovSync(model)
  current <- model[["POPCOV"]]
  if (is.null(current)) {
    stop("This model has no individually varying parameters, so there is no ",
      "population covariance to set.", call. = FALSE)
  }
  value <- as.matrix(value)
  if (is.null(dimnames(value)) || is.null(rownames(value))) {
    if (!identical(dim(value), dim(current))) {
      stop("POPCOV must be ", nrow(current), " x ", ncol(current),
        " for this model's ", nrow(current), " varying parameter",
        if (nrow(current) > 1L) "s" else "", ".", call. = FALSE)
    }
    dimnames(value) <- dimnames(current)
  }
  unknown <- setdiff(c(rownames(value), colnames(value)), rownames(current))
  if (length(unknown)) {
    stop("POPCOV names ", paste(unknown, collapse = ", "),
      ", which ", if (length(unknown) > 1L) "are" else "is",
      " not individually varying in this model.", call. = FALSE)
  }
  storage.mode(value) <- "character"
  current[rownames(value), colnames(value)] <- value
  model[["POPCOV"]] <- current
  model
}
