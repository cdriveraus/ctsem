#' Thresholds for ordinal manifest variables.
#'
#' An ordinal variable with K categories is modelled by the cumulative logit
#'
#'     P(y <= k | eta) = inv_logit(tau_k - eta)
#'
#' so it needs K-1 thresholds, which have to increase. The THRESHOLDS matrix
#' holds them one variable per row, but it does not hold the thresholds
#' themselves: column 1 is tau_1 and every later column is the *gap* to the
#' previous threshold, constrained positive by its transform. The engine
#' accumulates them.
#'
#' Storing gaps rather than thresholds is what makes the ordering constraint
#' free. Handed K-1 unconstrained cells an optimiser will cross them, and a
#' crossed pair gives the category between them probability zero -- the
#' likelihood is -Inf and there is no gradient pointing back out. Writing
#' `tau_1 + exp(delta_2) + ...` into each cell's transform instead would keep
#' real thresholds in the matrix, but a transform reading several free
#' parameters is not something the julia adjoint's parameter layer supports,
#' and supporting it would cost every model a gradient per cell.
#'
#' The consequence for the user is only in how the summary reads: row m column
#' 1 is that variable's first threshold, and columns 2+ are gaps.
#'
#' @noRd
NULL

#' Check and normalise the ncategories argument.
#' @noRd
.ctCheckNcategories <- function(ncategories, manifesttype, manifestNames) {
  n <- length(manifesttype)
  ordinal <- manifesttype %in% 2
  if (is.null(ncategories)) {
    if (!any(ordinal)) return(rep(0L, n))
    stop('manifesttype 2 (ordinal) needs ncategories -- the number of ',
      'categories of each ordinal variable. Ordinal variable(s): ',
      paste(manifestNames[ordinal], collapse = ', '), call. = FALSE)
  }
  if (length(ncategories) == 1L) ncategories <- rep(ncategories, n)
  if (length(ncategories) != n) stop('ncategories must have one entry per ',
    'manifest variable (', n, '), or be a single value', call. = FALSE)
  ncategories <- as.integer(ncategories)
  ncategories[is.na(ncategories)] <- 0L
  # Only the ordinal entries mean anything; zero the rest so nothing downstream
  # has to remember to ignore them.
  ncategories[!ordinal] <- 0L
  bad <- ordinal & ncategories < 3L
  if (any(bad)) stop('an ordinal variable needs at least 3 categories -- use ',
    'manifesttype 1 for a binary variable. Check: ',
    paste(manifestNames[bad], collapse = ', '), call. = FALSE)
  ncategories
}

#' Build the THRESHOLDS matrix for a model with ordinal variables.
#'
#' One row per manifest variable and as many columns as the widest ordinal
#' variable needs. Variables with fewer categories leave their trailing cells
#' fixed at zero and the engine never reads them, which is what lets ordinal
#' variables with different category counts share one rectangular matrix.
#' @noRd
.ctThresholdMatrix <- function(ncategories, manifesttype, manifestNames) {
  n <- length(manifesttype)
  ncol <- max(ncategories) - 1L
  out <- matrix(0, n, ncol,
    dimnames = list(manifestNames, paste0('threshold', seq_len(ncol))))
  for (i in seq_len(n)) {
    if (!manifesttype[i] %in% 2) next
    for (j in seq_len(ncategories[i] - 1L)) {
      out[i, j] <- paste0('threshold_', manifestNames[i], '_', j)
    }
  }
  out
}

#' Category counts implied by the data, for checking a model against it.
#' @noRd
.ctDataCategories <- function(datalong, ctm) {
  ordinal <- which(ctm$manifesttype %in% 2)
  if (!length(ordinal)) return(invisible(NULL))
  for (i in ordinal) {
    name <- ctm$manifestNames[i]
    v <- datalong[[name]]
    v <- v[!is.na(v)]
    if (!length(v)) next
    if (any(v != round(v)) || any(v < 1)) stop('ordinal variable ', name,
      ' must be coded as consecutive integers from 1 -- found values outside ',
      'that (min ', min(v), '). Shift 0-based codes up by one.', call. = FALSE)
    k <- ctm$ncategories[i]
    if (max(v) > k) stop('ordinal variable ', name, ' has values up to ',
      max(v), ' but the model was built with ncategories = ', k, call. = FALSE)
    # Fewer observed categories than declared is legal but rarely intended: the
    # thresholds bounding an empty category are not identified by the data.
    missingcats <- setdiff(seq_len(k), unique(v))
    if (length(missingcats)) warning('ordinal variable ', name,
      ' has no observations in categor', if (length(missingcats) > 1)
        'ies ' else 'y ', paste(missingcats, collapse = ', '),
      '; the thresholds bounding those are not identified by the data')
  }
  invisible(NULL)
}
