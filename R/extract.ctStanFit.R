#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#' @param subjectMatrices Calculate subject specific system matrices?
#' @param cores Only used if subjectMatrices = TRUE . For faster computation use more cores.
#' @return Array of posterior samples.
#' @aliases extract
#' @examples
#' \donttest{
#' e = ctExtract(ctstantestfit())
#' }
#' @export
ctExtract <- function(object,subjectMatrices=FALSE,cores=2){
  if(!class(object) %in% c('ctStanFit', 'stanfit')) stop('Not a ctStanFit or stanfit object')
  
  if(subjectMatrices & object$standata$savesubjectmatrices ==0){
    if(!'stanfit' %in% class(object$stanfit)) samps <- object$stanfit$rawposterior
    if('stanfit' %in% class(object$stanfit)) samps <- stan_unconstrainsamples(object$stanfit)
    out = stan_constrainsamples(sm = object$stanmodel,standata = object$standata,
      samples = samps,
      cores = cores,savescores = FALSE,savesubjectmatrices = TRUE,
      dokalman = TRUE,onlyfirstrow = TRUE)
    
  } else{
  if('stanfit' %in% class(object)) out <- rstan::extract(object) else{
  if('stanfit' %in% class(object$stanfit)) out <- rstan::extract(object$stanfit)
  if(!'stanfit' %in% class(object$stanfit)) out <- object$stanfit$transformedpars
  }
  out$Ygen[out$Ygen==99999] <- NA
  }
  return(out)
}
