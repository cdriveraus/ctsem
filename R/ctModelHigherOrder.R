#' Raise the order of a ctsem model object of type 'omx'.
#'
#' @param ctm ctModel
#' @param indices Vector of integers, which latents to raise the order of.
#' @param diffusion Shift the diffusion parameters / values to the higher order?
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
ctModelHigherOrder <- function(ctm, indices,diffusion=TRUE){
  ctm$latentNames <- c(ctm$latentNames,paste0('d',ctm$latentNames[indices]))
  nl <- ctm$n.latent
  
  for(i in indices){
    ctm$DRIFT <- rbind(cbind(ctm$DRIFT,0),0)
    ctm$DRIFT[i+nl,i+nl] <- ctm$DRIFT[i,i] #move ar param / value to higher order
    ctm$DRIFT[i,i] <- 0 #set ar param / value to 0
    ctm$DRIFT[i,i+nl] <- 1 #set effect of higher order to 1
    ctm$DRIFT[i+nl,i] <-paste0('drift_',ctm$latentNames[i+nl],'_',ctm$latentNames[i],'|-log1p(exp(-param*2))') #estimate effect of 1st order on 2nd
    if(i==tail(indices,1)) rownames(ctm$DRIFT) <- ctm$latentNames
    if(i==tail(indices,1)) colnames(ctm$DRIFT) <- ctm$latentNames
    
    for(m in c('T0VAR','DIFFUSION')){
      ctm[[m]] <- rbind(cbind(ctm[[m]],0),0)
      if(diffusion && m=='DIFFUSION'){
      ctm[[m]][i+nl,] <- ctm[[m]][i,] #move diffusion params / value to higher order
      ctm[[m]][,i+nl] <- ctm[[m]][,i] #move diffusion params / value to higher order
      ctm[[m]][i,] <- 0 #set order 1 diffusion to 0
      ctm[[m]][,i] <- 0 #set order 1 diffusion to 0
      }
      if(i==tail(indices,1)) colnames(ctm[[m]]) <- ctm$latentNames
      if(i==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
    
    for(m in c('CINT','T0MEANS')){
      ctm[[m]] <- rbind(ctm[[m]],0)
      if(m %in% 'T0MEANS') ctm[[m]][i+nl,1] <- paste0('T0mean_d_',ctm$latentNames[i])
      if(!m %in% 'T0MEANS') ctm[[m]][i,] <- 0 #set order 1 param to 0
      if(i==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
 
    ctm$LAMBDA <- cbind(ctm$LAMBDA,0)
  }
  
  for(i in c(indices,indices+nl)){
    for(j in c(indices,indices+nl)){
      if(i >= j) ctm$T0VAR[i,j] <- paste0('T0var_',ctm$latentNames[i],'_',
        ctm$latentNames[j])
    }
  }
  ctm$n.latent <- ctm$n.latent + length(indices)
  
  return(ctm)
}

