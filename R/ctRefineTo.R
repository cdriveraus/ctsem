#' ctRefineTo
#' 
#' Fits a ctsem model in a stepwise fashion to help with difficult optimization.
#' 
#' This function fits a sequence of ctsem models increasing in complexity, 
#' starting with a model involving fixed and relatively strong auto effects, no cross effects, no predictors, and no off-diagonal covariances.
#' For many models this can improve the speed and robustness of fitting
#' @param datawide Data in ctsem wide format
#' @param ctmodelobj A continuous time model specified via the \code{\link{ctModel}} function.
#' @param ... additional parameters to pass to \code{\link{ctFit}}.
#' @param modfunc function to run prior to each optimization step, that takes ctsem fit object, modifies it as desired, and returns the fit object.
#' @return Returns a fitted ctsem object in the same manner as \code{\link{ctFit}}.
#' @export

 

ctRefineTo<-function(datawide,ctmodelobj,modfunc=NULL,...){
  
  message('Fitting simplified model with any free DRIFT diagonals fixed to -.3, no free cross effects, no trait / diffusion correlations, no predictors, LAMBDA loadings fixed to 1')

  model<-ctmodelobj
  if(!is.null(model$TRAITVAR) & model$n.latent > 1) model$TRAITVAR[!diag(model$n.latent)]<-0
  if( model$n.latent > 1) model$DIFFUSION[!diag(model$n.latent)]<-0
  if( model$n.manifest > 1) model$MANIFESTVAR[!diag(model$n.manifest)]<-0
  #if( model$n.latent > 1) model$T0VAR[!diag(model$n.latent)]<-0
  if(!is.null(model$MANIFESTTRAITVAR) & model$n.manifest > 1) model$MANIFESTTRAITVAR[!diag(model$n.manifest)]<-0
  model$DRIFT[is.na(suppressWarnings(as.numeric(model$DRIFT)))]<-diag(-.3,model$n.latent)[is.na(suppressWarnings(as.numeric(model$DRIFT)))]
  model$DRIFT[row(model$DRIFT) != col(model$DRIFT)][is.na(suppressWarnings(as.numeric(model$DRIFT[row(model$DRIFT) != col(model$DRIFT)])))]<- 0
  model$LAMBDA[suppressWarnings(is.na(as.numeric(model$LAMBDA)))]<-1
  model$TDPREDEFFECT<-matrix(0,nrow=model$n.latent,ncol=model$n.TDpred)
  model$T0TDPREDCOV<-matrix(0,nrow=model$n.latent,ncol=(model$n.TDpred*(model$Tpoints-1)))
  model$TIPREDEFFECT<-matrix(0,nrow=model$n.latent,ncol=model$n.TIpred)
  model$T0TIPREDEFFECT<-matrix(0,nrow=model$n.latent,ncol=model$n.TIpred)
  model$TDTIPREDCOV<-matrix(0,nrow=(model$n.TDpred*(model$Tpoints-1)),ncol=model$n.TIpred)
  model$TRAITTDPREDCOV<-matrix(0,nrow=model$n.latent,ncol=(model$n.TDpred*(model$Tpoints-1)))
  if(model$n.TDpred > 0) model$TDPREDVAR[lower.tri(model$TDPREDVAR)]<-0
  if(model$n.TIpred > 0) model$TIPREDVAR[lower.tri(model$TIPREDVAR)]<-0
  fit<-ctFit(datawide,model,nofit=TRUE,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=F, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=F, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
  startValues<-OpenMx::omxGetParameters(fit$mxobj)
  
  
  
    message('Adding correlations, autoeffects, and any predictors, but no free cross effects')
  model<-ctmodelobj
  model$DRIFT[row(model$DRIFT) != col(model$DRIFT)][is.na(suppressWarnings(as.numeric(model$DRIFT[row(model$DRIFT) != col(model$DRIFT)])))]<- -.001
  fit<-ctFit(datawide,model,nofit=TRUE,omxStartValues=startValues,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=F, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=F, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
    startValues<-OpenMx::omxGetParameters(fit$mxobj)
    oldfit<-fit
  
    if(any(is.na(suppressWarnings(as.numeric(ctmodelobj$DRIFT[!diag(ctmodelobj$n.latent)]))))){
  message('Fitting complete user specified model')
  model<-ctmodelobj
  fit<-ctFit(datawide,model,nofit=TRUE,omxStartValues=startValues,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=F, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=F, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
    }
  
return(fit)
}
  
  
