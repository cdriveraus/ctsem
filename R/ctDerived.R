#' Quantities computed from a fit's matrices, draw by draw
#'
#' Computes a quantity from the model's matrices for each of a fit's draws, so
#' that it has its own median and interval. Use it for quantities ctsem does
#' not report, such as an explained variance or an indirect effect.
#'
#' @param fit A fit from \code{\link{ctFit}}, on either backend.
#' @param FUN A function of a list of the model's matrices, named as
#'   \code{\link{ctSummaryMatrices}} names them (\code{DRIFT},
#'   \code{DIFFUSIONcov}, \code{asymDIFFUSIONcov}, \code{dtDRIFT},
#'   \code{popcov} and so on). It returns a number, a vector or a matrix, of
#'   the same size for every draw.
#' @param ... Further arguments passed to \code{FUN}.
#' @param timeinterval Time interval for \code{dtDRIFT}.
#' @param state The latent state at which state dependent matrix cells are
#'   evaluated, as in \code{\link{ctSummaryMatrices}}. Needs
#'   \code{backend='julia'}.
#' @param nsamples \code{'all'}, or the number of draws to use.
#' @return An array of class \code{ctDerived}, with one row per draw and the
#'   quantity's own dimensions after it. For an optimised fit,
#'   \code{attr(x, 'estimate')} holds the quantity at the fit's estimate.
#'   \code{summary()} gives the mean, sd, 2.5\%, 50\% and 97.5\% of each
#'   element, and printing shows that summary.
#' @details The 50\% column is the point value. A quantity built from
#'   summaries of the matrices, such as the drift matrix of medians, is
#'   neither the median of that quantity nor a quantity of any one draw, so
#'   build it here instead. For an optimised fit whose draws come from the
#'   normal approximation at its estimate, an element that depends on a single
#'   parameter has its median at the estimate, and the 50\% column holds that
#'   value.
#' @examples
#' # The share of each process's stationary variance accounted for by the
#' # state one time unit earlier.
#' r2 <- ctDerived(ctsemtestfit, function(m) {
#'   A <- m$dtDRIFT
#'   S <- m$asymDIFFUSIONcov
#'   diag(A %*% S %*% t(A)) / diag(S)
#' })
#' r2
#'
#' # The draws themselves, one row per draw.
#' dim(r2)
#' @export
ctDerived <- function(fit, FUN, ..., timeinterval = 1, state = NULL,
  nsamples = 'all') {
  if (!inherits(fit, 'ctFit')) stop('fit must be a fit from ctFit().', call. = FALSE)
  FUN <- match.fun(FUN)
  model <- .ctFitModelObject(fit)
  build <- function(e, population) c(.ctSummaryArrays(e,
    continuoustime = model$continuoustime, latentNames = model$latentNames,
    manifestNames = model$manifestNames, TDpredNames = model$TDpredNames,
    timeinterval = timeinterval,
    binaryNames = model$manifestNames[model$manifesttype %in% 1],
    jacobians = TRUE), .ctPopulationArrays(population))

  if (.ctFitIsJulia(fit)) {
    resolved <- .ctResolveState(fit, state)
    samples <- .ctBackendRawSamples(fit)
    if (!identical(nsamples, 'all')) samples <- samples[seq_len(min(nrow(samples),
      as.integer(nsamples))), , drop = FALSE]
    drawn <- .ctBackendDrawArrays(fit, samples, state = resolved$state)
    arrays <- build(drawn$e, drawn$population)
    raw <- .ctFitEstimateRaw(fit)
    point <- if (!is.null(raw)) .ctBackendDrawArrays(fit, raw, state = resolved$state)
    atarrays <- if (!is.null(point)) build(point$e, point$population)
    centred <- if (!is.null(point)) .ctFitCentredDraws(fit, drawn$raw)
    .ctContextMessage(fit, resolved$label)
  } else {
    .ctContextRequireStateSupport(fit, state)
    e <- ctExtract(fit, cores = 1)
    arrays <- build(e, e)
    if (!identical(nsamples, 'all')) {
      keep <- seq_len(min(dim(arrays[[1L]])[1L], as.integer(nsamples)))
      arrays <- lapply(arrays, function(a)
        do.call(`[`, c(list(a, keep), rep(list(TRUE), length(dim(a)) - 1L),
          list(drop = FALSE))))
    }
    estimate <- .ctFitEstimateArrays(fit)
    atarrays <- if (!is.null(estimate)) build(estimate, estimate)
    centred <- NULL
  }

  ndraws <- dim(arrays[[1L]])[1L]
  value <- function(arrs, i) {
    n <- dim(arrs[[1L]])[1L]
    FUN(lapply(arrs, function(a) {
      d <- dim(a)
      array(a[i + n * (seq_len(prod(d[-1L])) - 1L)], dim = d[-1L],
        dimnames = dimnames(a)[-1L])
    }), ...)
  }
  first <- value(arrays, 1L)
  if (!is.numeric(first) && !is.logical(first)) stop('FUN must return numbers.',
    call. = FALSE)
  shape <- if (length(dim(first)) == 2L) dim(first) else length(first)
  labels <- if (length(dim(first)) == 2L) dimnames(first) else
    list(names(first))
  failed <- 0L
  draws <- vapply(seq_len(ndraws), function(i) {
    v <- tryCatch(value(arrays, i), error = function(err) NULL)
    if (is.null(v) || length(v) != prod(shape)) {
      failed <<- failed + 1L
      return(rep(NA_real_, prod(shape)))
    }
    as.numeric(v)
  }, numeric(prod(shape)))
  if (failed) message('FUN failed or changed size on ', failed, ' of ', ndraws,
    ' draws, which are NA.')
  draws <- matrix(draws, nrow = prod(shape))
  out <- array(t(draws), dim = c(ndraws, shape),
    dimnames = c(list(NULL), if (is.null(labels)) rep(list(NULL), length(shape))
      else labels))
  estimate <- if (!is.null(atarrays)) tryCatch(value(atarrays, 1L),
    error = function(err) NULL)
  if (!is.null(estimate) && length(estimate) == prod(shape)) {
    estimate <- if (length(shape) == 2L) matrix(as.numeric(estimate), shape[1L],
      shape[2L], dimnames = labels) else stats::setNames(as.numeric(estimate),
      labels[[1L]])
    attr(out, 'estimate') <- estimate
    attr(out, 'exactmedian') <- .ctExactMedians(t(draws), centred,
      as.numeric(estimate))
  }
  class(out) <- c('ctDerived', 'array')
  out
}

#' Summarise a quantity computed by ctDerived
#'
#' @param object Output of \code{\link{ctDerived}}.
#' @param probs Quantiles of each element over draws.
#' @param ... Not used.
#' @return A data frame with one row per element of the quantity, giving the
#'   mean, sd and quantiles of its draws. For a matrix, the \code{row} and
#'   \code{col} columns name the element.
#' @export
summary.ctDerived <- function(object, probs = c(.025, .5, .975), ...) {
  d <- dim(object)
  n <- d[1L]
  draws <- matrix(unclass(object), nrow = n)
  quantiles <- t(apply(draws, 2L, stats::quantile, probs = probs, na.rm = TRUE,
    names = FALSE))
  if (length(probs) == 1L) quantiles <- t(quantiles)
  colnames(quantiles) <- paste0(100 * probs, '%')
  known <- attr(object, 'exactmedian')
  middle <- which(abs(probs - .5) < 1e-12)
  if (length(middle) && length(known) == ncol(draws)) {
    quantiles[!is.na(known), middle] <- known[!is.na(known)]
  }
  stats <- data.frame(mean = colMeans(draws, na.rm = TRUE),
    sd = apply(draws, 2L, stats::sd, na.rm = TRUE), quantiles,
    check.names = FALSE)
  labels <- dimnames(object)[-1L]
  if (length(d) == 3L) {
    rows <- if (is.null(labels[[1L]])) seq_len(d[2L]) else labels[[1L]]
    cols <- if (is.null(labels[[2L]])) seq_len(d[3L]) else labels[[2L]]
    return(data.frame(row = rep(rows, d[3L]), col = rep(cols, each = d[2L]),
      stats, check.names = FALSE))
  }
  if (!is.null(labels[[1L]])) rownames(stats) <- labels[[1L]]
  stats
}

#' @export
print.ctDerived <- function(x, digits = 4L, ...) {
  s <- summary(x)
  values <- vapply(s, is.numeric, logical(1L)) & !names(s) %in% c('row', 'col')
  s[values] <- lapply(s[values], round, digits)
  cat('Medians and 95% intervals over', dim(x)[1L], 'draws.\n')
  print(s)
  invisible(x)
}
