#' Fits a ctsem model in a stepwise fashion to potentially improve optimization.
#' This function fits a sequence of ctsem models increasing in complexity, 
#' starting with a model involving fixed and high auto effects, no cross effects, and no trait parameters,
#' then introducing traits, then introducing any cross effects and freeing the parameters as appropriate to the specified model.
#' @param datawide Data in ctsem wide format
#' @param ctmodelobj A continuous time model specified via the \code{\link{ctModel}} function.
#' @return Returns a fitted ctsem object in the same manner as \code{\link{ctFit}}

 

ctFitStepwise<-function(datawide,ctmodelobj,...){
  
  message('Fitting simplified model with fixed high autoregression of .95, no cross effects, no traits, no predictors')

  model<-ctmodelobj
  model$TRAITVAR<-matrix(0,model$n.latent,model$n.latent)
#   model$MANIFESTVAR<-diag(.1,model$n.manifest)
  model$MANIFESTTRAITVAR<-matrix(0,model$n.latent,model$n.latent)
  model$DRIFT<-diag(-.05,model$n.latent)
  model$TDPRED<-matrix(0,nrow=model$n.latent,ncol=model$n.TDpred)
  model$T0DPREDEFFECT<-matrix(0,nrow=model$n.latent,ncol=model$n.TDpred)
  model$TIPRED<-matrix(0,nrow=model$n.latent,ncol=model$n.TIpred)
  model$T0TIPREDEFFECT<-matrix(0,nrow=model$n.latent,ncol=model$n.TIpred)
  fit<-ctFit(datawide,model,...)
  inits<-ctGetInits(fit)
  print(summary(fit$mxobj))
  
  
  
  if(any(ctmodelobj$TRAITVAR != 0) | any(ctmodelobj$MANIFESTTRAITVAR != 0)){
    message('Fitting simplified model with fixed autoeffects, no cross effects')
  model<-ctmodelobj
  model$DRIFT<-diag(-.05,model$n.latent)
    model$inits<-inits    
  fit<-ctFit(datawide,model,carefulFit=F...)
    inits<-ctGetInits(fit)
    print(summary(fit$mxobj))
  }
  
  message('Fitting complete user specified model')
  model<-ctmodelobj
  model$inits<-inits
  fit<-ctFit(datawide,model,carefulFit=F,...)
  
return(fit)
}
  
  