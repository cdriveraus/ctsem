#' Extract samples from a ctStanFit object
#'
#' @param fit ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#'
#' @return Array of posterior samples.
#' @export
#'
#' @examples
extract.ctStanFit <- function(fit){
  if(class(fit)!='ctStanFit') stop('Not a ctStanFit object')
  if(class(fit$stanfit)=='stanfit') out <- rstan::extract(fit$stanfit)
  if(class(fit$stanfit)!='stanfit') out <- fit$stanfit$transformedpars
  out$Ygen[out$Ygen==99999] <- NA
  return(out)
}
