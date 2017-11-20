#' Kalman
#'
#' Takes list containing ctsem subject matrices, as well as long form data object, and calculates 
#' predicted and updated latent states, likelihoods, and predicted observations using the Kalman filter.
#'
#' @param kpars list object containing DRIFT,T0VAR,DIFFUSION,CINT,T0MEANS,TDPREDEFFECT,
#' MANIFESTMEANS, LAMBDA, and MANIFESTVAR matrices, with list elements named accordingly. 
#' Such a list is returned by \code{\link{ctStanContinuousPars}}.
#' @param datalong long format data object as used by \code{\link{ctStanFit}}, 
#' but must contain only a single subjects' data and does not need an id column.
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
#' @param diffusionindices vector of integers denoting which latent variables are involved in covariance calcs.
#' @param optimize Set to TRUE when using for optimization.
#' @param binary Set to TRUE when using binary data.
#' @return When optimize=TRUE, returns log likelihood. Else, 
#' returns a list containing matrix objects etaprior, etaupd, etasmooth, y, yprior, 
#' yupd, ysmooth, prederror, time, loglik,  with values for each time point in each row. 
#' eta refers to latent states and y to manifest indicators - y itself is thus just 
#' the input data. 
#' Covariance matrices etapriorcov, etaupdcov, etasmoothcov, ypriorcov, yupdcov, ysmoothcov,  
#' are returned in a row * column * time array. 
#' @examples
#' ### ctstantestfit is a dummy ctStanFit object with 2 manifest indicators,
#' ###  4 latents, and 1 time dependent predictor.
#' 
#' ### get parameter matrices
#' kpars <- ctStanContinuousPars(ctstantestfit)
#' 
#' #construct dummy data
#' datalong <- cbind(0:9, matrix(rnorm(20,2,1),ncol=2))
#' datalong[c(1:3,9:10),2:3]<-NA #missing data to pre/fore cast
#' colnames(datalong) <- c('time', paste0('Y',1:2))
#' print(datalong)
#' 
#' #obtain Kalman filtered estimates
#' kout <- Kalman(kpars=kpars, datalong=datalong,
#'   manifestNames=paste0('Y',1:nrow(kpars$MANIFESTMEANS)),
#'   latentNames=paste0('eta',1:nrow(kpars$DRIFT)))
#' 
#' #print and plot smoothed estimates (conditional on all states) of indicators.
#' print(kout$ysmooth)
#' matplot(kout$time,kout$ysmooth,type='l')
#' matplot(kout$time,datalong[,2:3],type='p',add=TRUE,pch=1)
#' @export

Kalman<-function(kpars,datalong,
  manifestNames,latentNames,imputeMissings=FALSE,
  TDpredNames=NULL,
  continuoustime=TRUE,
  timecol='time', diffusionindices='all',optimize=FALSE,binary=FALSE){
  
  datalong=as.matrix(datalong)
  
  nmanifest=length(manifestNames)
  nlatent=length(latentNames)
  ntdpred=length(TDpredNames)
  
  if(any(c(nmanifest,nlatent) < 1)) stop('Length of manifestNames and latentNames must be greater than 0!')
  if(all(diffusionindices=='all')) diffusionindices=1:nlatent
  ndiffusion=length(diffusionindices)
  
  Y<-datalong[,manifestNames,drop=FALSE]
  if(ntdpred > 0) {
    tdpreds<-datalong[,TDpredNames,drop=FALSE]
    tdpreds[is.na(tdpreds)] <- 0
    # if(any(is.na(tdpreds))) stop('missingness in time dependent predictors! Kalman cannot run.')
  }
  
  if(continuoustime) {
    DRIFTHATCH <- (kpars$DRIFT[diffusionindices,diffusionindices] %x% diag(ndiffusion) + 
        diag(ndiffusion) %x% kpars$DRIFT[diffusionindices,diffusionindices,drop=FALSE]) 
    asymDIFFUSION <- matrix(-solve(DRIFTHATCH, c(kpars$DIFFUSION[diffusionindices,diffusionindices,drop=FALSE])), nrow=ndiffusion)
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
  
  etaprior[[1]]<-kpars$T0MEANS #tdpreds added below
  etapriorcov[[1]]<-kpars$T0VAR[diffusionindices,diffusionindices,drop=FALSE]
  loglik<-rep(0,nrow(datalong))
  observed<-list()
  
  discreteDRIFT<- list()
  discreteCINT<- list()
  discreteDIFFUSION <- list()
  
  for(rowi in 1:(nrow(datalong))){
    
    if(!binary){
    
    if(continuoustime){
      dt<-datalong[rowi,timecol]-datalong[max(c(1,rowi-1)),timecol]
      discreteDRIFT[[rowi]] <- expm(kpars$DRIFT * dt)
      discreteCINT[[rowi]] <- solve(kpars$DRIFT, (discreteDRIFT[[rowi]] - diag(nlatent))) %*% kpars$CINT
      discreteDIFFUSION[[rowi]] <- asymDIFFUSION - (discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*% 
          asymDIFFUSION %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE]))
    }
    if(!continuoustime){
      discreteDRIFT[[rowi]] <-kpars$DRIFT
      discreteCINT[[rowi]]<- kpars$CINT
      discreteDIFFUSION[[rowi]] <- kpars$DIFFUSION[diffusionindices,diffusionindices,drop=FALSE]
    }
    
    
    if(rowi>1){
      etaprior[[rowi]] <- discreteCINT[[rowi]]  + discreteDRIFT[[rowi]] %*% etaupd[[rowi-1]]
      etapriorcov[[rowi]] <-  discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*% 
        etaupdcov[[rowi-1]] %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE])  + discreteDIFFUSION[[rowi]] #check transpose
    }
    
    if(ntdpred > 0) etaprior[[rowi]] <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
    
    if(imputeMissings) Y[rowi,] <- 0 #etaprior[[rowi]] + t(chol(etapriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
    
    nafilter<-!is.na(Y[rowi,])
    observed[[rowi]]<-nafilter
    err[[rowi]]<-matrix(NA,nrow=nmanifest) #init prediction errors for this row
    
    # // one step ahead predictive distribution of y
    yprior[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaprior[[rowi]]
    ypriorcov[[rowi]] <- kpars$LAMBDA[,diffusionindices,drop=FALSE] %*% etapriorcov[[rowi]] %*% 
      t(kpars$LAMBDA[,diffusionindices,drop=FALSE]) + kpars$MANIFESTVAR

    if(imputeMissings) Y[rowi,] <- yprior[[rowi]] + t(chol(ypriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
    
    y <- Y[rowi,,drop=FALSE][,nafilter,drop=FALSE]
    
    #if all missing...
    if(all(!nafilter)){
      etaupd[[rowi]] <- etaprior[[rowi]]
      etaupdcov[[rowi]] <- etapriorcov[[rowi]]
      yupd[[rowi]] <- yprior[[rowi]]
      yupdcov[[rowi]] <- ypriorcov[[rowi]]
    }
    
    #if any not missing
    if(any(nafilter)){
      
      # // forecast error
      err[[rowi]][nafilter] <- as.numeric(y - yprior[[rowi]][nafilter,])
      
      # // Kalman gain
      # invypriorcov <- solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE])
      K<-matrix(0,nrow=nlatent,ncol=nmanifest)
      
      K[diffusionindices,nafilter] <-  t(solve(t(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]),
        t(etapriorcov[[rowi]] %*%  t(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) ) ))
      
      # updated distribution 
      etaupd[[rowi]] <- etaprior[[rowi]] + K[,nafilter,drop=FALSE] %*% (err[[rowi]][nafilter,,drop=FALSE])
      
      # etaupdcov[[rowi + 1]] <- etapriorcov[[rowi]]
      etaupdcov[[rowi]] <- (diag(ndiffusion) - K[diffusionindices,nafilter,drop=FALSE] %*% 
          kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) %*% etapriorcov[[rowi]]
      
      # // log likelihood
      loglik[rowi] <- - 0.5 * (nrow(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) * log(2 * pi)  + 
          log(det(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]))    + 
          t(err[[rowi]][nafilter,,drop=FALSE]) %*% 
          solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE],err[[rowi]][nafilter,,drop=FALSE]))
      
      yupd[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaupd[[rowi]]
      yupdcov[[rowi]] <- kpars$LAMBDA[,diffusionindices,drop=FALSE] %*% etaupdcov[[rowi]] %*% 
        t(kpars$LAMBDA[,diffusionindices,drop=FALSE]) + kpars$MANIFESTVAR
      }
      
    }
      if(binary){
        # browser()
        if(rowi==1) {
          etaprior[[rowi]] <- list(sig1=etaprior[[rowi]],sig2=etaprior[[rowi]] )
        yprior[[rowi]] <- etaprior[[rowi]]
        } 
        
        if(rowi > 1) {
          etaprior[[rowi]] <- list(sig1=etaprior[[rowi-1]],sig2=etaprior[[rowi-1]] )
          yprior[[rowi]] <- etaprior[[rowi]]
        } 
        
        if(continuoustime){
          dt<-datalong[rowi,timecol]-datalong[max(c(1,rowi-1)),timecol]
          discreteDRIFT[[rowi]] <- expm(kpars$DRIFT * dt)
          discreteCINT[[rowi]] <- solve(kpars$DRIFT, (discreteDRIFT[[rowi]] - diag(nlatent))) %*% kpars$CINT
          discreteDIFFUSION[[rowi]] <- asymDIFFUSION - (discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*% 
              asymDIFFUSION %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE]))
        }
        if(!continuoustime){
          discreteDRIFT[[rowi]] <-kpars$DRIFT
          discreteCINT[[rowi]]<- kpars$CINT
          discreteDIFFUSION[[rowi]] <- kpars$DIFFUSION[diffusionindices,diffusionindices,drop=FALSE]
        }
        
        
        if(rowi>1){
          for(sigi in 1:2){
            sigm <- ifelse(sigi==1,-1,1)
          etaprior[[rowi]][[sigi]] <- discreteCINT[[rowi]]  + discreteDRIFT[[rowi]] %*% (etaupd[[rowi-1]] + sigm*t(chol(discreteDIFFUSION[[rowi]])))
          }
          etapriorcov[[rowi]] <-  discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*% 
            etaupdcov[[rowi-1]] %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE])  + discreteDIFFUSION[[rowi]] #check transpose
        }
        
        if(ntdpred > 0) stop('binary!') #etaprior[[rowi]] <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
        
        if(imputeMissings) Y[rowi,] <- 0 #etaprior[[rowi]] + t(chol(etapriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
        
        nafilter<-!is.na(Y[rowi,])
        observed[[rowi]]<-nafilter
        err[[rowi]]<-matrix(NA,nrow=nmanifest) #init prediction errors for this row
        
        # // one step ahead predictive distribution of y
        for(sigi in 1:2){
        yprior[[rowi]][[sigi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaprior[[rowi]][[sigi]]
        }
        ypriorcov[[rowi]] <- kpars$LAMBDA[,diffusionindices,drop=FALSE] %*% etapriorcov[[rowi]] %*% 
          t(kpars$LAMBDA[,diffusionindices,drop=FALSE]) + kpars$MANIFESTVAR
        
        if(imputeMissings) Y[rowi,] <- yprior[[rowi]] + t(chol(ypriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
        
        y <- Y[rowi,,drop=FALSE][,nafilter,drop=FALSE]
        
        #if all missing...
        if(all(!nafilter)){
          etaupd[[rowi]] <- etaprior[[rowi]]
          etaupdcov[[rowi]] <- etapriorcov[[rowi]]
          yupd[[rowi]] <- mean(yprior[[rowi]][[1]],yprior[[rowi]][[2]])
          yupdcov[[rowi]] <- ypriorcov[[rowi]]
        }
        
        #if any not missing
        if(any(nafilter)){
          ypriorlinked <- list()
          err[[rowi]] <- list(sig1=err[[rowi]],sig2=err[[rowi]])
          for(sigi in 1:2){
        ypriorlinked[[sigi]] <- 1 / (1+exp(-yprior[[rowi]][[sigi]][nafilter,]))
          }
        
        # loglik[rowi] <- log( (1-(y*(1-ypriorlinked))) * ( ( 1- ypriorlinked) * y) )
        loglik[rowi] <- mean(sum(log( y*(ypriorlinked[[1]]) + (1-y)*(1-ypriorlinked[[1]]))),
          sum(log( y*(ypriorlinked[[2]]) + (1-y)*(1-ypriorlinked[[2]]))) )
        # browser()
        
        for(sigi in 1:2){
        err[[rowi]][[sigi]][nafilter] <- as.numeric(y - ypriorlinked[[sigi]])
        }
        # err[[rowi]][nafilter] <- log(err[[rowi]][nafilter] / (1-err[[rowi]][nafilter]))
        
        # // Kalman gain
        # invypriorcov <- solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE])
        K<-matrix(0,nrow=nlatent,ncol=nmanifest)
        
        K[diffusionindices,nafilter] <-  t(solve(t(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]),
          t(etapriorcov[[rowi]] %*%  t(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) ) ))
        
        # updated distribution 
        # browser()

        etaupd[[rowi]] <- ( (etaprior[[rowi]][[1]] +  (K[,nafilter,drop=FALSE] *1) %*% (err[[rowi]][[1]][nafilter,,drop=FALSE])) +
          (etaprior[[rowi]][[2]] + (K[,nafilter,drop=FALSE] *1) %*% (err[[rowi]][[2]][nafilter,,drop=FALSE])) ) /2 
        
        # etaupdcov[[rowi + 1]] <- etapriorcov[[rowi]]
        etaupdcov[[rowi]] <- (diag(ndiffusion) - K[diffusionindices,nafilter,drop=FALSE] %*% 
            kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) %*% etapriorcov[[rowi]]
        
        # // log likelihood
        
        
        yupd[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaupd[[rowi]]
        yupdcov[[rowi]] <- kpars$LAMBDA[,diffusionindices,drop=FALSE] %*% etaupdcov[[rowi]] %*% 
          t(kpars$LAMBDA[,diffusionindices,drop=FALSE]) + kpars$MANIFESTVAR
        
 
      }
      
    }
    
    # }
  }
  if(binary){
    # browser()
  yprior2 <- array(unlist(yprior),dim=c(length(yprior[[1]][[1]]),2,nrow(Y)))
  plot(1 / (1+exp(-yprior2[1,1,])),type='b',col='red',ylim=c(0,1))
  points(1 / (1+exp(-yprior2[1,2,])),type='b',col='blue')
  points(apply(Y,1,mean))
  }
  
  #extra output
  if(!optimize){
  
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
      smoother<-diag(0,nlatent)
      smoother[diffusionindices,diffusionindices]<- t(solve(t(etapriorcov[[rowi+1]]), t( etaupdcov[[rowi]] %*% 
        t(discreteDRIFT[[rowi+1]][diffusionindices,diffusionindices,drop=FALSE]) ) )) #is the rowi+1 correct?
        
      etasmooth[[rowi]]<-etaupd[[rowi]]+smoother %*% (etasmooth[[rowi+1]] - etaprior[[rowi+1]])
      etasmoothcov[[rowi]]<-etaupdcov[[rowi]] + smoother[diffusionindices,diffusionindices,drop=FALSE] %*% 
        ( etasmoothcov[[rowi+1]] - etapriorcov[[rowi+1]])
    }
    ysmooth[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etasmooth[[rowi]]
    ysmoothcov[[rowi]] <- kpars$LAMBDA[,diffusionindices,drop=FALSE] %*% etasmoothcov[[rowi]] %*% 
      t(kpars$LAMBDA[,diffusionindices,drop=FALSE]) + kpars$MANIFESTVAR
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
  
  y<-Y #matrix(datalong[,manifestNames,drop=FALSE],ncol=nmanifest) #
  colnames(y)<-manifestNames
  
  if(ntdpred>0) {
    tdpreds<-datalong[,TDpredNames,drop=FALSE]
    colnames(tdpreds)<-TDpredNames
  }
  
  err<-matrix(unlist(err),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))
  
  ysmooth<-matrix(unlist(ysmooth),byrow=T,ncol=nmanifest,
    dimnames=list(timedims,manifestNames))
  
  
  insertdiffusionzeros <- function(x){
    y <- array(0,c(nlatent,nlatent,length(timedims)))
    y[diffusionindices,diffusionindices,] <- x
    dimnames(y) <- list(latentNames,latentNames,timedims)
    return(y)
  }
  
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

  etapriorcov <- insertdiffusionzeros(array(unlist(etapriorcov),
    dim=c(length(diffusionindices),length(diffusionindices),length(timedims))))
  
  etaupdcov <- insertdiffusionzeros(array(unlist(etaupdcov),
    dim=c(length(diffusionindices),length(diffusionindices),length(timedims))))
  
  etasmoothcov <- insertdiffusionzeros(array(unlist(etasmoothcov),
    dim=c(length(diffusionindices),length(diffusionindices),length(timedims))))
  
  names(loglik) = timedims
  
  out<-list(observed,
    etaprior=etaprior,etapriorcov=etapriorcov,
    etaupd=etaupd,etaupdcov=etaupdcov,
    loglik=loglik, 
    prederror=err,y=y,
    yprior=yprior,ypriorcov=ypriorcov,
    yupd=yupd,yupdcov=yupdcov,
    etasmooth=etasmooth, etasmoothcov=etasmoothcov, 
    ysmooth=ysmooth, ysmoothcov=ysmoothcov)
  
  if(ntdpred > 0) out$tdpreds=tdpreds
  if(continuoustime) out$time=datalong[,timecol]
  }
  
  if(optimize) out <- sum(loglik,na.rm=TRUE)
  
  return(out)
}


