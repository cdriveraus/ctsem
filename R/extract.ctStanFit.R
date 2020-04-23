#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#' @param ... additional arguments to pass to \code{rstan::extract}.
#' @return Array of posterior samples.
#' @aliases extract
#' @examples
#' \donttest{
#' if (!exists("ctstantestfit")) example(ctstantestfit)
#' e = ctExtract(ctstantestfit)
#' }
#' @export
ctExtract <- function(object,...){
  if(!class(object) %in% c('ctStanFit', 'stanfit')) stop('Not a ctStanFit or stanfit object')
  if('stanfit' %in% class(object)) out <- rstan::extract(object,...) else{
  if('stanfit' %in% class(object$stanfit)) out <- rstan::extract(object$stanfit,...)
  if(!'stanfit' %in% class(object$stanfit)) out <- object$stanfit$transformedpars
  out$Ygen[out$Ygen==99999] <- NA
  }
  return(out)
}
