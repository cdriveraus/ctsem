#'ctDiscreteToContinuous
#'Converts values in ctModel matrices from discrete time to continuous time parameters.
#'@param ctmodelobj model specified via ctModel function. Must contain only fixed numeric values, 
#'parameter values should represent those from a discrete time model.
#'@param timeInterval time interval the provided discrete parameters represent. 1 may be appropriate.
#'@details Does not convert T0TRAITEFFECT or T0TIPREDEFFECT matrices yet.
#'@export
ctDiscreteToContinuous <- function(ctmodelobj,timeInterval){
  
  
  for(i in 1:length(ctmodelobj)){ #convert matrices to numeric
    if(is.matrix(ctmodelobj[[i]])){
      ctmodelobj[[i]] <- matrix(as.numeric(ctmodelobj[[i]]),nrow=nrow(ctmodelobj[[i]]))
    }
  }
  
  n.latent<-ctmodelobj$n.latent
  II<-diag(n.latent)
  
  ctmodelobj$DRIFT <- OpenMx::logm(ctmodelobj$DRIFT)
  
  DRIFTHATCH<-(ctmodelobj$DRIFT %x% diag(n.latent) + diag(n.latent) %x% ctmodelobj$DRIFT)
  
  ctmodelobj$DIFFUSION <- matrix(DRIFTHATCH %*% solve((OpenMx::expm(DRIFTHATCH %x% timeInterval)) - (II%x%II)) %*% 
    OpenMx::rvectorize(ctmodelobj$DIFFUSION),nrow(II))
  
  ctmodelobj$CINT <- ctmodelobj$DRIFT %*% solve(OpenMx::expm(ctmodelobj$DRIFT %x% timeInterval) - II) %*% 
      (ctmodelobj$CINT)
  
  smalltraitloadings <- solve(ctmodelobj$DRIFT) %*% (-OpenMx::expm(ctmodelobj$DRIFT %x% timeInterval)) - solve(ctmodelobj$DRIFT)
  
  if(!is.null(ctmodelobj$TRAITVAR)) ctmodelobj$TRAITVAR<-matrix(
    solve(smalltraitloadings %x% smalltraitloadings) %*% c(ctmodelobj$TRAITVAR),nrow=n.latent)
  
  if(ctmodelobj$n.TIpred > 0) ctmodelobj$TIPREDEFFECT <- -ctmodelobj$DRIFT %*% 
    solve(diag(n.latent) - OpenMx::expm(ctmodelobj$DRIFT %x% timeInterval) )    %*%  (ctmodelobj$TIPREDEFFECT)
  
  if(ctmodelobj$n.TDpred > 0) ctmodelobj$TDPREDEFFECT <- 
    solve(OpenMx::expm(ctmodelobj$DRIFT %x% timeInterval)) %*% ctmodelobj$TDPREDEFFECT

  return(ctmodelobj)
}