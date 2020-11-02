ctModeltoNumeric <- function(ctmodelobj){
  ###read in model
  #set any matrices to numeric elements
  domessage<-FALSE
  sapply(names(ctmodelobj), function(x){
    if(is.matrix(ctmodelobj[[x]])){
      m <- ctmodelobj[[x]]
      if(any(suppressWarnings(is.na(as.numeric(m))))){
        domessage <<-TRUE
      }
      m[suppressWarnings(is.na(as.numeric(m)))] <- 0
      ctmodelobj[[x]] <<- matrix(as.numeric(m),nrow=nrow(m), ncol=ncol(m))
    }
  })
  if(domessage) message('Some free parameters in matrices set to zero!')
  
  return(ctmodelobj)
}

#' ctGenerate
#' 
#' This function generates data according to the specified ctsem model object. 
#' 
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dtmean Positive numeric. Average time interval (delta T) to use.
#' @param logdtsd Numeric. Standard deviation for variability of the time interval.
#' @param dtmat Either NA, or numeric matrix of n.subjects rows and burnin+Tpoints-1 columns, 
#' containing positive numeric values for all time intervals between measurements. 
#' If not NA, dtmean and logdtsd are ignored.
#' @param wide Logical. Output in wide format?
#' @details TRAITVAR and MANIFESTRAITVAR are treated as Cholesky factor covariances 
#' of CINT and MANIFESTMEANS, respectively. 
#' TRAITTDPREDCOV and TIPREDCOV matrices are not accounted for, at present. 
#' The first 1:n.TDpred rows and columns of TDPREDVAR are used for generating
#' tdpreds at each time point. 
#' @examples 
#' #generate data for 2 process model, each process measured by noisy indicator, 
#' #stable individual differences in process levels.
#' 
#' generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'  MANIFESTVAR=diag(.1,2),
#'  LAMBDA=diag(1,2),
#'  DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
#'  TRAITVAR=matrix(c(.5,.2,0,.8),nrow=2),
#'  DIFFUSION=matrix(c(1,.2,0,4),2),
#'  CINT=matrix(c(1,0),nrow=2),
#'  T0MEANS=matrix(0,ncol=1,nrow=2),
#'  T0VAR=diag(1,2))
#'
#' data<-ctGenerate(generatingModel,n.subjects=15,burnin=10)
#' @export

ctGenerate<-function(ctmodelobj,n.subjects=100,burnin=0,dtmean=1,logdtsd=0,dtmat=NA,
  wide=FALSE){
  
  ctmodelobj <- ctModeltoNumeric(ctmodelobj)
  
  m <- ctmodelobj
  
  fullTpoints<-burnin+m$Tpoints

  for(si in 1:n.subjects){
    
    if(is.na(dtmat[1])) dtvec<- exp(rnorm(fullTpoints,log(dtmean),logdtsd))
    if(!is.na(dtmat[1])) dtvec <- dtmat[si,,drop=FALSE]
    time=rep(0,fullTpoints)
    for(t in 2:fullTpoints) time[t] = round(time[t-1] + dtvec[t-1],6)
    
    if(m$n.TDpred > 0) {
      tdpreds <- rbind(matrix(0,nrow=1+(burnin),ncol=m$n.TDpred)[-1,,drop=FALSE], #additional row added then removed in case no burnin
        matrix(m$TDPREDMEANS + m$TDPREDVAR %*% rnorm(nrow(m$TDPREDVAR),0,1),ncol=m$n.TDpred))
    }
    
    sm=m
    if(any(m$TRAITVAR != 0)) {
      traits = m$TRAITVAR %*% rnorm(m$n.latent,0,1)
      sm$CINT = sm$CINT +  traits
      sm$T0MEANS = sm$T0MEANS + m$T0TRAITEFFECT %*% traits
    }
    
    if(any(m$MANIFESTTRAITVAR != 0)) {
      sm$MANIFESTMEANS = sm$MANIFESTMEANS + m$MANIFESTTRAITVAR %*% rnorm(m$n.manifest,0,1)
    }
    
    if(m$n.TIpred > 0) {
      tipreds <- m$TIPREDMEANS + m$TIPREDVAR %*% rnorm(m$n.TIpred,0,1)
      sm$CINT = sm$CINT + m$TIPREDEFFECT %*% tipreds
    }
    
    manifests<-matrix(NA,fullTpoints,m$n.manifest)
    latents<-matrix(NA,fullTpoints,m$n.latent)
    
    sdat <- cbind(si,time,manifests,
      if(m$n.TDpred > 0) tdpreds,
      if(m$n.TIpred > 0) matrix(tipreds,byrow=TRUE,nrow=fullTpoints,ncol=m$n.TIpred))
    
    colnames(sdat) <- c('id','time',m$manifestNames,m$TDpredNames,m$TIpredNames)
    sdat <- data.frame(sdat)
    
    # sdat[,manifestNames] <- Kalman(kpars=skpars,datalong=sdat,manifestNames=manifestNames,TDpredNames=TDpredNames,
    #   latentNames=latentNames,
    #   imputeMissings=TRUE)$y
    
    latents[1,] <- sm$T0MEANS+sm$T0VAR %*% rnorm(m$n.latent)
    Qinf <- fQinf(sm$DRIFT,sm$DIFFUSION)
    for(i in 2:nrow(latents)){
      dtA=expm::expm(sm$DRIFT * (sdat$time[i]-sdat$time[i-1]))
      latents[i,] <- dtA %*% latents[i-1,] +
        solve(sm$DRIFT,(dtA - diag(m$n.latent))) %*% sm$CINT + 
        t(chol(fdtQ(Qinf,dtA))) %*% rnorm(m$n.latent)
      # browser()
      if(m$n.TDpred > 0) latents[i,] <- latents[i,] + sm$TDPREDEFFECT %*% 
        t(as.matrix(sdat[i,m$TDpredNames,drop=FALSE]))
    }
    
    for(i in 1:nrow(sdat)){
      sdat[i,m$manifestNames] <- sm$LAMBDA %*% latents[i,] + sm$MANIFESTMEANS + 
        sm$MANIFESTVAR %*% rnorm(m$n.manifest)
    }
        
    
    
    
    sdat=sdat[(burnin+1):fullTpoints,]
    
    sdat[,'time'] = sdat[,'time'] - sdat[1,'time'] 
    
    if(si==1) datalong <- sdat else datalong <- as.matrix(rbind(datalong,sdat))
  }
  
  
  if(wide==FALSE) return(datalong) else {
    datawide <- ctLongToWide(datalong = datalong,id = 'id',time = 'time',
      manifestNames = m$manifestNames, TDpredNames = m$TDpredNames,TIpredNames = m$TIpredNames)
    datawide <- ctIntervalise(datawide = datawide,Tpoints = m$Tpoints,n.manifest = m$n.manifest,n.TDpred = m$n.TDpred,n.TIpred = m$n.TIpred,
      manifestNames=m$manifestNames,TDpredNames=m$TDpredNames,TIpredNames=m$TIpredNames)
    return(datawide)
  }
}
