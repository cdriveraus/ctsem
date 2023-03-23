

#'ctStanContinuousPars
#'
#'Returns the continuous time parameter matrices of a ctStanFit fit object
#'
#'@param fit fit object from \code{\link{ctStanFit}}
#'@param calcfunc Function to apply over samples, must return a single value. 
#'By default the median over all samples is returned using the \code{\link[stats]{quantile}} function, 
#'but one might also be interested in the \code{\link[base]{mean}} or \code{\link[stats]{sd}}, for instance.
#'@param calcfuncargs A list of additional parameters to pass to calcfunc. 
#'For instance, with the default of calcfunc = quantile, 
#'the probs argument is needed to ensure only a single value is returned.
#'@param timeinterval time interval for discrete time parameter matrix computation.
#'@examples
#'\donttest{
#'#posterior median over all subjects (also reflects mean of unconstrained pars)
#'ctStanContinuousPars(ctstantestfit)
#'}
#'@export
ctStanContinuousPars <- function(fit,
  calcfunc=quantile,calcfuncargs=list(probs=0.5),timeinterval=1){

  if(!'ctStanFit' %in% class(fit)) stop(paste0('Not an object of class ctStanFit! Instead is ',paste0(class(fit),collapse=', ')))
  
  e<-ctExtract(fit,cores=1) #Qfit$stanfit$transformedpars #first dim of subobjects is iter, 2nd subjects
  niter=dim(e$pop_DRIFT)[1]

  


  mats <- ctStanMatricesList()
  mats <- c(names(mats$base), names(mats$asymptotic),names(mats$extra))
  if(fit$ctstanmodel$continuoustime){
    d=list(DRIFT=e$pop_DRIFT)
    dd=ctStanDiscreteParsDrift(d,timeinterval,observational = FALSE,standardise = FALSE,cov = FALSE,quiet=TRUE)
    e$pop_dtDRIFT <- array(dd,dim=dim(dd)[-2:-3])
    mats <- c(mats, 'dtDRIFT')
  }

  out <- list()
  for(matname in (mats)){
    try({
      calcfuncargs$collapsemargin = 1
    calcfuncargs$collapsefunc=calcfunc

      calcfuncargs$inarray = e[[paste0('pop_',matname)]]
      out[[matname]] <- array(do.call(ctCollapse,calcfuncargs),
        dim=dim(calcfuncargs$inarray)[-1])
    },silent=TRUE)
  }

  if(nrow(out$T0MEANS) > nrow(out$CINT)){ #then intoverpop used...
    nlatent <- nrow(out$CINT)
    out$T0MEANS <- out$T0MEANS[1:nlatent,1,drop=FALSE]
    out$DRIFT <- out$DRIFT[1:nlatent,1:nlatent,drop=FALSE]
    out$T0VAR <- out$T0VAR[1:nlatent,1:nlatent,drop=FALSE]
    out$T0cov <- out$T0cov[1:nlatent,1:nlatent,drop=FALSE]
  }
  
  ln=fit$ctstanmodel$latentNames
  mn=fit$ctstanmodel$manifestNames
  tdn=fit$ctstanmodel$TDpredNames
  dimnames(out$DRIFT)=list(ln,ln)
  dimnames(out$DIFFUSIONcov)=list(ln,ln)
  dimnames(out$DIFFUSION)=list(ln,ln)
  dimnames(out$T0cov)=list(ln,ln)
  dimnames(out$asymDIFFUSIONcov)=list(ln,ln)
  rownames(out$CINT)=ln
  rownames(out$MANIFESTMEANS)=mn
  rownames(out$T0MEANS)=ln
  
  dimnames(out$T0VAR)=list(ln,ln)
  dimnames(out$asymDIFFUSIONcov)=list(ln,ln)
  dimnames(out$LAMBDA)=list(mn,ln)

  
  if(!is.null(e$pop_MANIFESTVAR)) {
    dimnames(out$MANIFESTVAR)=list(mn,mn)
    dimnames(out$MANIFESTcov)=list(mn,mn)
    # out$MANIFESTVAR=out$MANIFESTVAR %*% t(out$MANIFESTVAR) #cholesky factor inside stanfit...
    
  }
  
  if(!is.null(e$pop_TDPREDEFFECT)) {
    dimnames(out$TDPREDEFFECT)=list(ln,tdn)
  }

  out$MANIFESTVAR <- NULL ; 
  
  
  return(out)
}

