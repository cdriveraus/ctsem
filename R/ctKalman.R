#' ctKalman
#'
#' Takes list containing ctsem subject matrices, as well as long form data object, and calculates 
#' predicted and updated latent states, likelihoods, and predicted observations using the Kalman filter.

#' @param kpars list object containing DRIFT,T0VAR,DIFFUSION,CINT,T0MEANS,TDPREDEFFECT,
#' MANIFESTMEANS, LAMBDA, and MANIFESTVAR matrices, with list elements named accordingly. 
#' Such a list is returned by \code{\link{ctStanContinuousPars}}.
#' @param datalong long format data object as used by \code{\link{ctStanFit}}, 
#' but must contain only a single subjects' data and does not need an id column.
#' @param subject integer denoting which subjects data to use.  
#' @param manifestNames String vector of names of manifest variables to use from datalong.
#' @param latentNames String vector of names of latent variables.
#' @param TDpredNames If model contains time dependent predictors, 
#' string vector of their names in the data.
#' @param imputeMissings Logical. If TRUE, randomly generate any missing observations of 
#' manifest variables according to model.
#' @param continuoustime Logical, whether to use a continuous time Kalman filter or discrete time. 
#' Refers only to latent states, observations are always at discrete time points.
#' @param timecol name of time column in datalong. Note that time column must be an ascending sequence
#' of numeric values from row 1 to row n. Ignored if continuoustime=FALSE.
#' @return Returns a list containing matrix objects etaprior, etaupd, etasmooth, y, yprior, 
#' yupd, ysmooth, prederror, time, loglik,  with values for each time point in each row. 
#' eta refers to latent states and y to manifest indicators - y itself is thus just 
#' the input data. 
#' Covariance matrices etapriorcov, etaupdcov, etasmoothcov, ypriorcov, yupdcov, ysmoothcov,  
#' are returned in a row * column * time array. 
#' @examples
#ctstantestfit is a dummy ctStanFit object with 2 manifest indicators,
# 4 latents, and 1 time dependent predictor.
#' 
#get parameter matrices
#' kpars <- ctStanContinuousPars(ctstantestfit)
#' 
#' #construct dummy data
#' datalong <- cbind(0:9, matrix(rnorm(20,2,1),ncol=2))
#' datalong[c(1:3,9:10),2:3]<-NA #missing data to pre/fore cast
#' colnames(datalong) <- c('time', paste0('Y',1:2))
#' print(datalong)
#' 
#' #obtain Kalman filtered estimates
#' kout <- ctKalman(kpars=kpars, datalong=datalong,
#'   subject=1, manifestNames=paste0('Y',1:2),
#'   latentNames=paste0('eta',1:4))
#' 
#' #print and plot smoothed estimates (conditional on all states) of indicators.
#' print(kout$ysmooth)
#' matplot(kout$time,kout$ysmooth,type='l')
#' matplot(kout$time,datalong[,2:3],type='p',add=TRUE,pch=1)
#' @export

ctKalman<-function(kpars,datalong,
  manifestNames,latentNames,imputeMissings=FALSE,
  TDpredNames=NULL,
  continuoustime=TRUE,
  timecol='time'){
  
  nmanifest=length(manifestNames)
  nlatent=length(latentNames)
  ntdpred=length(TDpredNames)
  
  if(any(c(nmanifest,nlatent) < 1)) stop('Length of manifestNames and latentNames must be greater than 0!')

  
  Y<-datalong[,manifestNames,drop=FALSE]
  if(ntdpred > 0) {
    tdpreds<-datalong[,TDpredNames,drop=FALSE]
    if(any(is.na(tdpreds))) stop('missingness in time dependent predictors! ctKalman cannot run.')
  }
  
  if(continuoustime) {
    DRIFTHATCH <- (kpars$DRIFT %x% diag(nlatent) + diag(nlatent) %x% kpars$DRIFT) 
    asymDIFFUSION <- matrix(-solve(DRIFTHATCH) %*% c(kpars$DIFFUSION), nrow=nrow(kpars$DRIFT))
  }
  
  etaprior<-list()
  etapriorcov<-list()
  etaupd<-list()
  etaupdcov<-list()
  err<-list()
  yprior<-list()
  ypriorcov<-list()
  yupd<-list()
  yupdcov<-list()
  
  etaprior[[1]]<-kpars$T0MEANS
  if(ntdpred > 0) etaprior[[1]]<-etaprior[[1]]+kpars$TDPREDEFFECT %*% t(datalong[1,TDpredNames,drop=FALSE])
  etapriorcov[[1]]<-kpars$T0VAR
  loglik<-rep(0,nrow(datalong))
  observed<-list()
  
  discreteDRIFT<- list()
  discreteCINT<- list()
  discreteDIFFUSION <- list()
  
  for(rowi in 1:(nrow(datalong))){
    
    if(continuoustime){
      dt<-datalong[rowi,timecol]-datalong[max(c(1,rowi-1)),timecol]
      discreteDRIFT[[rowi]] <- OpenMx::expm(kpars$DRIFT * dt)
      discreteCINT[[rowi]] <- solve(kpars$DRIFT) %*% (discreteDRIFT[[rowi]] - diag(nlatent)) %*% kpars$CINT
      discreteDIFFUSION[[rowi]] <- asymDIFFUSION - (discreteDRIFT[[rowi]] %*% asymDIFFUSION %*% t(discreteDRIFT[[rowi]]))
    }
    if(!continuoustime){
      discreteDRIFT[[rowi]] <-kpars$DRIFT
      discreteCINT[[rowi]]<- kpars$CINT
      discreteDIFFUSION[[rowi]] <- kpars$DIFFUSION
    }
    # }
    if(rowi>1){
      etaprior[[rowi]] <- discreteCINT[[rowi]]  + discreteDRIFT[[rowi]] %*% etaupd[[rowi-1]]
      if(ntdpred > 0) etaprior[[rowi]] <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
      etapriorcov[[rowi]] <-  discreteDRIFT[[rowi]] %*% etaupdcov[[rowi-1]] %*% t(discreteDRIFT[[rowi]])  + discreteDIFFUSION[[rowi]] #check transpose
    }
    
    if(imputeMissings) Y[rowi,] <- 0 #etaprior[[rowi]] + t(chol(etapriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
    
    nafilter<-!is.na(Y[rowi,])
    observed[[rowi]]<-nafilter
    
    # // one step ahead predictive distribution of y
    yprior[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaprior[[rowi]]
    ypriorcov[[rowi]] <- kpars$LAMBDA %*% etapriorcov[[rowi]] %*% t(kpars$LAMBDA) + kpars$MANIFESTVAR
    
    if(imputeMissings) Y[rowi,] <- yprior[[rowi]] + t(chol(ypriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
    
    y <- Y[rowi,,drop=FALSE][,nafilter,drop=FALSE]
    
    # // forecast error
    err[[rowi]]<-matrix(NA,nrow=nmanifest)
    err[[rowi]][nafilter] <- as.numeric(y - yprior[[rowi]][nafilter,])
    
    #if all missing...
    if(all(!nafilter)){
      etaupd[[rowi]] <- etaprior[[rowi]]
      etaupdcov[[rowi]] <- etapriorcov[[rowi]]
      yupd[[rowi]] <- yprior[[rowi]]
      yupdcov[[rowi]] <- ypriorcov[[rowi]]
    }
    
    #if any not missing
    if(any(nafilter)){
      
      # // Kalman gain
      invypriorcov <- solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE])
      K<-matrix(0,nrow=nlatent,ncol=nmanifest)
      K[,nafilter] <-  etapriorcov[[rowi]] %*%   t(kpars$LAMBDA[nafilter,,drop=FALSE]) %*% invypriorcov
      
      # // updated distribution 
      # etaupd[[rowi + 1]] <- etaprior[[rowi]]
      
      etaupd[[rowi]] <- etaprior[[rowi]] + K[,nafilter,drop=FALSE] %*% (err[[rowi]][nafilter,,drop=FALSE])
      
      # etaupdcov[[rowi + 1]] <- etapriorcov[[rowi]]
      etaupdcov[[rowi]] <- (diag(ncol(kpars$LAMBDA)) - K[,nafilter,drop=FALSE] %*% kpars$LAMBDA[nafilter,,drop=FALSE]) %*% etapriorcov[[rowi]]
      
      # // log likelihood
      loglik[rowi] <- - 0.5 * (nrow(kpars$LAMBDA[nafilter,,drop=FALSE]) * log(2 * pi)  + 
          log(det(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]))    + 
          t(err[[rowi]][nafilter,,drop=FALSE]) %*% invypriorcov %*% (err[[rowi]][nafilter,,drop=FALSE]))
      
      yupd[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaupd[[rowi]]
      yupdcov[[rowi]] <- kpars$LAMBDA %*% etaupdcov[[rowi]] %*% t(kpars$LAMBDA) + kpars$MANIFESTVAR
      
    }
    # }
  }
  
  
  #smoother
  etasmooth<-list()
  etasmoothcov<-list()
  ysmooth<-list()
  ysmoothcov<-list()
  
  for(rowi in nrow(datalong):1){
    if(rowi==nrow(datalong)) {
      etasmooth[[rowi]]<-etaupd[[rowi]]
      etasmoothcov[[rowi]]<-etaupdcov[[rowi]]
    } else{
      smoother<- etaupdcov[[rowi]] %*% t(discreteDRIFT[[rowi+1]]) %*% #is the rowi+1 correct?
        solve(etapriorcov[[rowi+1]])
      etasmooth[[rowi]]<-etaupd[[rowi]]+smoother %*% (etasmooth[[rowi+1]] - etaprior[[rowi+1]])
      etasmoothcov[[rowi]]<-etaupdcov[[rowi]] + smoother %*% ( etasmoothcov[[rowi+1]] - etapriorcov[[rowi+1]])
    }
    ysmooth[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etasmooth[[rowi]]
    ysmoothcov[[rowi]] <- kpars$LAMBDA %*% etasmoothcov[[rowi]] %*% t(kpars$LAMBDA) + kpars$MANIFESTVAR
  }
  
  timedims=paste0('t',datalong[,timecol])
  
  etaprior<-matrix(unlist(etaprior),byrow=T,ncol=nlatent,
    dimnames=list(timedims,latentNames))

  etaupd<-matrix(unlist(etaupd),byrow=T,ncol=nlatent,
    dimnames=list(timedims,latentNames))
  
  etasmooth<-matrix(unlist(etasmooth),byrow=T,ncol=nlatent,
    dimnames=list(timedims,latentNames))
  
  yprior<-matrix(unlist(yprior),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))
  
  yupd<-matrix(unlist(yupd),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))
  
  y<-Y #matrix(datalong[,manifestNames,drop=FALSE],ncol=nmanifest)
  colnames(y)<-manifestNames
  
  if(ntdpred>0) {
    tdpreds<-datalong[,TDpredNames,drop=FALSE]
    colnames(tdpreds)<-TDpredNames
  }
  
  err<-matrix(unlist(err),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))

  ysmooth<-matrix(unlist(ysmooth),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))
  
  etapriorcov <- plyr::alply(etapriorcov,1,function(x) x,.dims=TRUE)
  dimnames(etapriorcov)
  
  
  ypriorcov <- array(unlist(ypriorcov),
    dim=c(nmanifest,nmanifest,length(timedims)),
    dimnames=list(manifestNames,manifestNames,timedims))
  
  yupdcov <- array(unlist(yupdcov),
    dim=c(nmanifest,nmanifest,length(timedims)),
    dimnames=list(manifestNames,manifestNames,timedims))
  
  ysmoothcov <- array(unlist(ysmoothcov),
    dim=c(nmanifest,nmanifest,length(timedims)),
    dimnames=list(manifestNames,manifestNames,timedims))
  
  etapriorcov <- array(unlist(etapriorcov),
    dim=c(nlatent,nlatent,length(timedims)),
    dimnames=list(latentNames,latentNames,timedims))
  
  etaupdcov <- array(unlist(etaupdcov),
    dim=c(nlatent,nlatent,length(timedims)),
    dimnames=list(latentNames,latentNames,timedims))
  
  etasmoothcov <- array(unlist(etasmoothcov),
    dim=c(nlatent,nlatent,length(timedims)),
    dimnames=list(latentNames,latentNames,timedims))
  
  names(loglik) = timedims
  
  out<-list(observed,etaprior=etaprior,etapriorcov=etapriorcov,
    etaupd=etaupd,etaupdcov=etaupdcov,loglik=loglik, 
    prederror=err,y=y,yprior=yprior,ypriorcov=ypriorcov,
    yupd=yupd,yupdcov=yupdcov,
    etasmooth=etasmooth, etasmoothcov=etasmoothcov, ysmooth=ysmooth, ysmoothcov=ysmoothcov)
  
  if(ntdpred > 0) out$tdpreds=tdpreds
  if(continuoustime) out$time=datalong[,timecol]
  
  return(out)
}


