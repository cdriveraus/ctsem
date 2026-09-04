# Base generics for the stan backend.
#
# `coef`, `logLik` and `print` were registered for `ctJuliaFit` only, so the
# same call worked on a julia fit and fell through to the default -- an error,
# or a screenful of list -- on a stan fit. Nothing here computes a new
# statistic: each method reads what `summary.ctStanFit()` and the existing stan
# accessors already produce, so the two backends answer the same question the
# same way.
#
# This whole file goes when the stan backend goes. The julia counterparts live
# at the end of R/ctJuliaBackend.R.

#' @export
coef.ctStanFit <- function(object, ...) {
  raw <- .ctFitRawEstimate(object)
  # A sampled fit has no point estimate in `$rawest`; the posterior mean of the
  # raw draws is what `summary()` reports for it, so it is what is returned
  # here too.
  if (is.null(raw)) raw <- colMeans(ctStanRawSamples(object))
  as.numeric(raw)
}

#' @export
logLik.ctStanFit <- function(object, ...) {
  ll <- object$stanfit$transformedparsfull$ll
  if (is.null(ll) || !length(ll)) {
    stop('logLik needs a single log likelihood, which a sampled fit does not ',
      'have -- summary(fit)$loglik gives its distribution instead.',
      call. = FALSE)
  }
  nobs <- object$standata$ndatapoints
  structure(as.numeric(ll)[1L],
    df = length(.ctFitRawEstimate(object)),
    nobs = if (is.null(nobs)) NA_integer_ else as.integer(nobs),
    class = 'logLik')
}

#' @export
print.ctStanFit <- function(x, ...) {
  cat('ctsem Stan fit\n')
  sampled <- !is.null(x$stanfit$stanfit) && length(x$stanfit$stanfit@sim) > 0
  ll <- x$stanfit$transformedparsfull$ll
  if (!is.null(ll) && length(ll)) cat('  log likelihood:', format(ll[1L]), '\n')
  if (sampled) {
    cat('  sampled\n')
  } else {
    # The reason optimisation stopped, verbatim: stan's optimisers report
    # 'no step found' as a termination, so there is no honest boolean here to
    # match print.ctJuliaFit's `converged`.
    optimfit <- x$stanfit$optimfit
    cat('  terminated:',
      if (is.null(optimfit$terminate$what)) 'unknown' else optimfit$terminate$what,
      ' iterations:', if (is.null(optimfit$iter)) NA else optimfit$iter, '\n')
  }
  invisible(x)
}
