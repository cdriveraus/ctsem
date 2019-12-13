#' Raise the order of a ctsem model object of type 'omx'.
#'
#' @param ctm ctModel
#' @param indices Vector of integers, which latents to raise the order of.
#' @param diffusion Shift the diffusion parameters / values to the higher order?
#' @param crosseffects Shift cross coupling parameters of the DRIFT matrix to the higher order?
#' @param cint shift continuous intercepts to higher order?
#' @param explosive Allow explosive (non equilibrium returning) processes?
#'
#' @return extended ctModel
#' @export
#'
#' @examples
#' om <- ctModel(LAMBDA=diag(1,2),DRIFT=0, 
#'   MANIFESTMEANS=0,type='omx',Tpoints=4)
#'   
#' om <- ctModelHigherOrder(om,1:2)
#' print(om$DRIFT)
#' 
#' m <- ctStanModel(om)
#' print(m$pars)
ctModelHigherOrder <- function(ctm, indices,diffusion=TRUE, crosseffects=FALSE,cint=FALSE, explosive=FALSE){
  ctm$latentNames <- c(ctm$latentNames,paste0('d',ctm$latentNames[indices]))
  nl <- ctm$n.latent
  
  for(i in 1:length(indices)){
    for(m in c('DRIFT','DIFFUSION','T0VAR')){
    ctm[[m]] <- rbind(cbind(ctm[[m]],0),0)
    }
    for(m in c('CINT','T0MEANS')){
      ctm[[m]] <- rbind(ctm[[m]],0)
    }
    ctm$LAMBDA <- cbind(ctm$LAMBDA,0)
    if(ctm$n.TDpred > 0) ctm$TDPREDEFFECT <- rbind(ctm$TDPREDEFFECT,0)
  }
  
  for(i in 1:length(indices)){
    # browser()
    if(!crosseffects) ctm$DRIFT[i+nl,i+nl] <- ctm$DRIFT[indices[i],indices[i]] #move ar param / value to higher order
    if(crosseffects){
      m <- 'DRIFT'
      # browser()
      ctm[[m]][i+nl,] <- ctm[[m]][indices[i],] #move crosseffect params / value to higher order
      ctm[[m]][,i+nl] <- ctm[[m]][,indices[i]] #move crosseffect params / value to higher order
      ctm[[m]][indices[i],] <- 0 #set order 1 crosseffect to 0
      ctm[[m]][,indices[i]] <- 0 #set order 1 crosseffect to 0
    }
    ctm$DRIFT[indices[i],indices[i]] <- 0 #set ar param / value to 0
    ctm$DRIFT[indices[i],i+nl] <- 1 #set effect of higher order to 1
    ctm$DRIFT[i+nl,indices[i]] <-paste0('drift_',ctm$latentNames[i+nl],'_',
      ctm$latentNames[indices[i]],
      ifelse(!explosive,'|-5*log1p(exp(-param*2))','|param*5-2')) #estimate effect of 1st order on 2nd
    if(indices[i]==tail(indices,1)) rownames(ctm$DRIFT) <- ctm$latentNames
    if(indices[i]==tail(indices,1)) colnames(ctm$DRIFT) <- ctm$latentNames
    
    for(m in c('T0VAR','DIFFUSION')){
      if(diffusion && m=='DIFFUSION'){
        ctm[[m]][i+nl,] <- ctm[[m]][indices[i],] #move diffusion params / value to higher order
        ctm[[m]][,i+nl] <- ctm[[m]][,indices[i]] #move diffusion params / value to higher order
        ctm[[m]][indices[i],] <- 0 #set order 1 diffusion to 0
        ctm[[m]][,indices[i]] <- 0 #set order 1 diffusion to 0
      }
      if(indices[i]==tail(indices,1)) colnames(ctm[[m]]) <- ctm$latentNames
      if(indices[i]==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
    
    for(m in c('CINT','T0MEANS')){
      # browser()
      if(m %in% 'T0MEANS') ctm[[m]][i+nl,1] <- paste0('T0mean_d_',ctm$latentNames[indices[i]])
      if(!m %in% 'T0MEANS' && cint){
        ctm[[m]][indices[i]+nl,] <- ctm[[m]][indices[i],]
        ctm[[m]][indices[i],] <- 0 #set order 1 param to 0
        }
      if(indices[i]==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
  }
  # browser()
  for(i in c(nl+seq_along(indices))){
    for(j in c(1:nl,nl+seq_along(indices))){
      if(i >= j) ctm$T0VAR[i,j] <- paste0('T0var_',ctm$latentNames[i],'_',
        ctm$latentNames[j])
    }
  }
  ctm$n.latent <- ctm$n.latent + length(indices)
  
  return(ctm)
}

