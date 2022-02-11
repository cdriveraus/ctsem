#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#' @param subjectMatrices Calculate subject specific system matrices?
#' @param cores Only used if subjectMatrices = TRUE . For faster computation use more cores.
#' @param nsamples either 'all' or an integer denoting number of random samples to extract.
#' @return Array of posterior samples.
#' @aliases extract
#' @examples
#' \donttest{
#' e = ctExtract(ctstantestfit)
#' }
#' @export
ctExtract <- function(object,subjectMatrices=FALSE,cores=2,nsamples='all'){
  if(!class(object) %in% c('ctStanFit', 'stanfit')) stop('Not a ctStanFit or stanfit object')
  
  
  
  if(length(object$stanfit$stanfit@sim)==0){
    samps <- object$stanfit$rawposterior
    if(!nsamples %in% 'all') samps <- samps[sample(1:nrow(samps),nsamples),,drop=FALSE]
    if(subjectMatrices && object$standata$savesubjectmatrices==0){
      out = stan_constrainsamples(sm = object$stanmodel,standata = object$standata,
        samples = samps,
        cores = cores,savescores = FALSE,savesubjectmatrices = subjectMatrices,
        dokalman = TRUE,onlyfirstrow = FALSE)
    } else out <- object$stanfit$transformedpars
  }
  
  if(length(object$stanfit$stanfit@sim)>0){
    if(subjectMatrices & object$standata$savesubjectmatrices!=1){
      samps <- t(stan_unconstrainsamples(object$stanfit$stanfit,standata=object$standata))
      if(!nsamples %in% 'all') samps <- samps[sample(1:nrow(samps),nsamples),,drop=FALSE]
      out = stan_constrainsamples(sm = object$stanmodel,standata = object$standata,
        samples = samps,
        cores = cores,savescores = FALSE,savesubjectmatrices = subjectMatrices)
    } else  out <- rstan::extract(object$stanfit$stanfit)
  } 
  
  out$Ygen[out$Ygen==99999] <- NA
  
  # if(!is.null(out$rawpopc)){
  #   out$rawpopcov <- array(out$rawpopc[,4,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcorr <-  array(out$rawpopc[,3,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcovchol <-  array(out$rawpopc[,2,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcovbase <-  array(out$rawpopc[,1,,],dim=dim(out$rawpopc)[-2])
  # }
  
  return(out)
}
