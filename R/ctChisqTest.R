#' Chi Square test wrapper for ctStanFit objects.
#'
#' @param fit1 One of the fits to be compared (better fit is assumed as base for comparison)
#' @param fit2 Second fit to be compared
#'
#' @return Numeric probability
#' @export
#'
#' @examples
#' \donttest{
#'   if(w32chk()){ #skips on 32 bit systems
#'     df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)
#'     
#'     m1 <- ctModel(type='standt', LAMBDA=diag(1),MANIFESTVAR=0)
#'     m2 <- ctModel(type='standt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)
#'     
#'     f1 <- ctStanFit(df,m1,cores=1)
#'     f2 <- ctStanFit(df,m2,cores=1)
#'     
#'     ctChisqTest(f1,f2)
#'   }
#' }

ctChisqTest<-function(fit1,fit2){
  d <- data.frame(
    ll=c(fit1$stanfit$optimfit$value,
      fit2$stanfit$optimfit$value),
    npars=c(length(fit1$stanfit$rawest),
      length(fit2$stanfit$rawest))
  )
  d <- d[order(d$npars,decreasing = FALSE),]
  stats::pchisq(q =  diff(2*d$ll),df = diff(d$npars),lower.tail = FALSE)
}
