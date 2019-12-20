#' Generate data from a ctstanmodel object
#'
#' @param ctm \code{\link{ctStanModel}} object.
#' @param datastruct long format data structure as used by ctsem.
#' @param optimize Whether to optimize or use Stan's HMC sampler
#' @param is If optimizing, follow up with importance sampling? 
#' @param fullposterior Generate from the full posterior or just the mean?
#' @param nsamples How many samples to generate?
#' @param parsonly If TRUE, only return samples of raw parameters, don't generate data.
#' @param ... arguments to pass to stanoptimis 
#'
#' @return Array of nsamples x time points x manifest variables.
#' @export
#'
#' @examples
#' \dontrun{
#' m1 <- ctModel(type = 'stanct',
#' manifestNames = c('Exercise1'), 
#' latentNames=c('Exercise'),
#' DRIFT= 0,
#' DIFFUSION=0,
#' CINT='cint1',
#' T0MEANS='t0m1',
#' T0VAR=0, #only need to set this when t0means not individually varying
#' LAMBDA = 1,
#' MANIFESTMEANS=0,
#' MANIFESTVAR='merror')
#'
#' 
#' #generate and plot samples from prior predictive
#' priorpred <- ctStanGenerate(ctm = m1,datastruct = exfitdat,cores=6,nsamples = 50)
#'}
ctStanGenerate <- function(ctm,datastruct, optimize=TRUE, is=FALSE, fullposterior=TRUE, nsamples=200, parsonly=FALSE,...){
  datastruct[,ctm$manifestNames] <- NA
  dots <- list(...)
  dots$carefulfit=FALSE
  dots$is <- is
  dots$tol=1e-18
  if(is.null(dots$finishsamples) && parsonly) dots$finishsamples=nsamples
  #problem with multiple cores inside function?
  pp<-ctStanFit(datalong = datastruct[c(1,nrow(datastruct)),,drop=FALSE], 
    ctstanmodel = ctm,optimize=optimize, optimcontrol=dots,cores=1,verbose=0)

  if(parsonly) dat <- pp else{
  
  filled <- datastruct
  filled[,ctm$manifestNames] <- -99
  # browser()
  ppf <- ctStanFit(datalong = filled, ctstanmodel = ctm,optimize=optimize, optimcontrol=dots,fit=FALSE)
  # pp$standata$Y <- ppf$standata$Y
  ppf$stanfit <- pp$stanfit
  class(ppf) <- c('ctStanFit',class(ppf))
  ppf <- ctStanGenerateFromFit(fit = ppf,nsamples = nsamples,fullposterior = fullposterior)
  
  dat <- ppf$generated$Y
  dimnames(dat) <- list(datapoints=1:dim(dat)[1], samples=1:dim(dat)[2], manifests = ctm$manifestNames)
  }
  
  
  return(dat)
}
  
