#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#'
#' @return Array of posterior samples.
#' @export
#' @examples
#' e = extract.ctStanFit(ctstantestfit)
#' head(e)
extract.ctStanFit <- function(object){
  if(class(object)!='ctStanFit') stop('Not a ctStanFit object')
  if(class(object$stanfit)=='stanfit') out <- rstan::extract(object$stanfit)
  if(class(object$stanfit)!='stanfit') out <- object$stanfit$transformedpars
  out$Ygen[out$Ygen==99999] <- NA
  return(out)
}
