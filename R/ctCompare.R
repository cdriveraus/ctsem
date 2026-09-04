#' Chi Square test wrapper for ctsem fit objects.
#'
#' Likelihood ratio test between two nested model fits. Works for \code{ctStanFit}
#' and \code{ctJuliaFit} objects, including comparisons between one of each, via
#' \code{.ctFitOptimValue()}/\code{.ctFitRawEstimate()} (see \code{R/ctBackendSummary.R}),
#' which read the objective value actually optimised (stan's \code{optimfit$value},
#' julia's \code{estimate$logposterior} falling back to \code{estimate$loglik} --
#' both are the log posterior including any priors, not the pure log likelihood)
#' and the raw parameter vector for each backend.
#'
#' @param fit1 One of the fits to be compared (better fit is assumed as base for comparison)
#' @param fit2 Second fit to be compared
#'
#' @return Numeric probability
#' @export
#'
#' @examples
#' \donttest{
#'     df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)
#'
#'     m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
#'     m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)
#'
#'     f1 <- ctFit(df,m1,cores=1)
#'     f2 <- ctFit(df,m2,cores=1)
#'
#'     ctChisqTest(f1,f2)
#' }

ctChisqTest<-function(fit1,fit2){
  # A likelihood-ratio test needs the objective value *at a point estimate*,
  # which a sampled fit does not have -- it has a posterior. Stan's sampling
  # path never sets `stanfit$optimfit` (that field belongs to `stanoptimis()`
  # alone; see R/ctFit.R), so `.ctFitOptimValue()` silently returns NULL for a
  # sampled ctStanFit, and the `c(NULL, <value>)` below used to shrink `ll` to
  # length 1 against `npars` at length 2 -- `data.frame()` then failed with
  # "arguments imply differing number of rows", which names neither fit nor
  # says why. Checked and named here instead, for both backends: julia's
  # equivalent marker is `fit$sample` (set only by the routine
  # `ctFit(optimize=FALSE)` and `ctSample()` share).
  sampled <- function(fit) {
    if (inherits(fit, 'ctJuliaFit')) !is.null(fit$sample) else
      length(fit$stanfit$stanfit@sim) > 0
  }
  if (sampled(fit1) || sampled(fit2)) {
    stop("ctChisqTest() compares point estimates and needs an optimized fit ",
      "for both models; a sampled fit (ctFit(optimize=FALSE) or ctSample()) ",
      "has a posterior instead of a single objective value. Refit with ",
      "optimize=TRUE, or compare sampled fits with ctLOO() instead.",
      call.=FALSE)
  }
  d <- data.frame(
    ll=c(.ctFitOptimValue(fit1),
      .ctFitOptimValue(fit2)),
    npars=c(length(.ctFitRawEstimate(fit1)),
      length(.ctFitRawEstimate(fit2)))
  )
  d <- d[order(d$npars,decreasing = FALSE),]
  stats::pchisq(q =  diff(2*d$ll),df = diff(d$npars),lower.tail = FALSE)
}
