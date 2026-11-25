# Choice and response time indicators (manifesttype 5): the drift diffusion
# model's Wiener first-passage density, integrated exactly by the julia filter
# (inst/julia/ContinuousTimeSEM/src/wiener_measurement.jl).
#
# The indicator's value is the signed response time in seconds: positive for a
# response at the upper boundary, negative at the lower. The drift is the
# indicator's linear predictor, LAMBDA eta + MANIFESTMEANS; trial to trial
# drift variability is the MANIFESTVAR diagonal (a standard deviation, as a
# count's dispersion is); boundary separation, nondecision time and relative
# start are the indicator's row of WIENER.

.ctWienerColumns <- c('boundary', 'ndt', 'bias')
.ctWienerDefaults <- c('1', '0', '0.5')

#' Check the WIENER matrix argument.
#'
#' One row per manifest variable, columns boundary, ndt and bias, so the engine
#' reads a row by the manifest's own index. The rows of other indicators are
#' unread and sit at fixed defaults. NA, or no matrix, gives a Wiener row a free
#' boundary and nondecision time and an unbiased start. A length 3 vector fills
#' every Wiener row.
#' @noRd
.ctCheckWienerMatrix <- function(WIENER, manifesttype, manifestNames) {
  wiener <- manifesttype %in% 5L
  given <- !(is.null(WIENER) || all(is.na(WIENER)))
  if (!any(wiener)) {
    if (given) stop('WIENER applies to choice and response time indicators ',
      '(manifesttype 5), and this model has none', call. = FALSE)
    return(NULL)
  }
  n <- length(manifestNames)
  out <- matrix(.ctWienerDefaults, n, 3, byrow = TRUE,
    dimnames = list(manifestNames, .ctWienerColumns))
  out[wiener, 1] <- paste0('boundary_', manifestNames[wiener])
  out[wiener, 2] <- paste0('ndt_', manifestNames[wiener])
  if (given) {
    if (!is.matrix(WIENER)) {
      if (length(WIENER) != 3L) stop('WIENER must be a matrix with one row per ',
        'manifest variable and three columns (boundary, ndt, bias), or a ',
        'length 3 vector applied to every choice and response time indicator',
        call. = FALSE)
      fill <- matrix(as.character(WIENER), sum(wiener), 3, byrow = TRUE)
      blank <- is.na(fill) | trimws(fill) %in% 'NA'
      out[wiener, ][!blank] <- fill[!blank]
    } else {
      if (!identical(dim(WIENER), c(n, 3L))) stop('WIENER needs one row per ',
        'manifest variable (', n, ') and three columns (boundary, ndt, bias); ',
        'got ', paste(dim(WIENER), collapse = ' * '), call. = FALSE)
      given <- matrix(as.character(WIENER), n, 3)
      blank <- is.na(given) | trimws(given) %in% 'NA'
      if (any(!blank[!wiener, ])) stop('WIENER applies to choice and response ',
        'time indicators only; leave the rows of ',
        paste(manifestNames[!wiener & rowSums(!blank) > 0], collapse = ', '),
        ' as NA', call. = FALSE)
      out[!blank] <- given[!blank]
    }
  }
  .ctCheckWienerValues(out[wiener, , drop = FALSE])
  out
}

# Fixed values a Wiener row cannot take.
.ctCheckWienerValues <- function(cells) {
  num <- suppressWarnings(matrix(as.numeric(cells), nrow(cells), 3))
  if (any(num[, 1] <= 0, na.rm = TRUE)) stop('a fixed WIENER boundary ',
    'separation must be positive', call. = FALSE)
  if (any(num[, 2] < 0, na.rm = TRUE)) stop('a fixed WIENER nondecision time ',
    'cannot be negative', call. = FALSE)
  if (any(num[, 3] <= 0 | num[, 3] >= 1, na.rm = TRUE)) stop('a fixed WIENER ',
    'bias is the start as a share of the boundary separation, so lies ',
    'strictly between 0 and 1', call. = FALSE)
  invisible(NULL)
}

#' The fit-time half of `.ctCheckWienerMatrix`: `pars` and `manifesttype` can
#' both be edited after ctModel.
#' @noRd
.ctCheckWienerPars <- function(pars, manifesttype) {
  wiener <- which(manifesttype %in% 5L)
  r <- pars[pars$matrix %in% 'WIENER', , drop = FALSE]
  if (!length(wiener)) {
    if (nrow(r) && any(!is.na(r$param))) stop('WIENER has free cells but the ',
      'model has no choice and response time indicators (manifesttype 5)',
      call. = FALSE)
    return(invisible(NULL))
  }
  if (!nrow(r) || max(r$row) != length(manifesttype)) stop('choice and response ',
    'time indicators need a WIENER matrix with one row per manifest variable; ',
    'rebuild the model with ctModel()', call. = FALSE)
  if (any(!is.na(r$param) & !r$row %in% wiener)) stop('WIENER applies to choice ',
    'and response time indicators only (manifesttype 5); rows ',
    paste(unique(r$row[!is.na(r$param) & !r$row %in% wiener]), collapse = ', '),
    ' are not', call. = FALSE)
  fixed <- r[r$row %in% wiener & is.na(r$param), , drop = FALSE]
  if (nrow(fixed)) {
    cells <- matrix(NA_character_, length(manifesttype), 3)
    cells[cbind(fixed$row, fixed$col)] <- as.character(fixed$value)
    .ctCheckWienerValues(cells[wiener, , drop = FALSE])
  }
  invisible(NULL)
}

#' Choice and response time data: each value the signed response time in
#' seconds, positive at the upper boundary.
#' @noRd
.ctDataWiener <- function(datalong, ctm) {
  for (name in ctm$manifestNames[ctm$manifesttype %in% 5L]) {
    y <- suppressWarnings(as.numeric(datalong[[name]]))
    y <- y[!is.na(y)]
    if (!length(y)) next
    if (any(!is.finite(y) | y == 0)) stop('choice and response time indicator ',
      name, ' has values that are zero or not finite; each value is the ',
      'response time in seconds, negative for a response at the lower ',
      'boundary', call. = FALSE)
    if (all(y > 0) || all(y < 0)) warning('every response of ', name, ' is at ',
      if (all(y > 0)) 'the upper' else 'the lower', ' boundary. Values are ',
      'signed response times: negative for a response at the lower boundary.',
      call. = FALSE)
    if (stats::median(abs(y)) > 20) warning('the median response time of ',
      name, ' is ', signif(stats::median(abs(y)), 3), '; response times are in ',
      'seconds', call. = FALSE)
  }
  invisible(NULL)
}
