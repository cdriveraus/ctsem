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
#' @param ekf set to TRUE to use the extended Kalman filter, for non-linear measurement models.
#' @param binary Set to TRUE when using binary data.
#' @param plotoptim set to TRUE to plot / print optimization steps.
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
  timecol='time', diffusionindices='all',optimize=FALSE,ekf=FALSE, binary=FALSE,plotoptim=FALSE){
  
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
  
  
  loglik<-rep(0,nrow(datalong))
  observed<-list()
  
  discreteDRIFT<- list()
  discreteCINT<- list()
  discreteDIFFUSION <- list()
  
  t0check <- c(1,as.numeric(datalong[-1,'id']!=datalong[-nrow(datalong),'id']))
  
  if(continuoustime) dt<-c(-9,datalong[-1,timecol]-datalong[-nrow(datalong),timecol])
  
  if(!ekf){
    for(rowi in 1:(nrow(datalong))){
      
      if(continuoustime){
        discreteDRIFT[[rowi]] <- expm(kpars$DRIFT * dt[rowi])
        discreteCINT[[rowi]] <- solve(kpars$DRIFT, (discreteDRIFT[[rowi]] - diag(nlatent))) %*% kpars$CINT
        discreteDIFFUSION[[rowi]] <- asymDIFFUSION - (discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*% 
            asymDIFFUSION %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE]))
      }
      if(!continuoustime){
        discreteDRIFT[[rowi]] <-kpars$DRIFT
        discreteCINT[[rowi]]<- kpars$CINT
        discreteDIFFUSION[[rowi]] <- kpars$DIFFUSION[diffusionindices,diffusionindices,drop=FALSE]
      }
      
      if(t0check[rowi]==1){
        etaprior[[rowi]]<-kpars$T0MEANS #tdpreds added below
        etapriorcov[[rowi]]<-kpars$T0VAR[diffusionindices,diffusionindices,drop=FALSE]
      }
      
      
      if(t0check[rowi]==0){
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
  }
  
  
  if(ekf){
    
    if(binary) measurementmodel <- function(x){
      1/(1+exp(-kpars$LAMBDA %*% x))
    }
    
    if(!binary) measurementmodel <- function(x){
      kpars$LAMBDA %*% x
    }
    
    dynamicmodel <- function(x,jacobian=FALSE){
      # browser()
      if(!is.null(PARMEANS)) {
        PARMEANS <- x[(nlatent+1):length(x)]
        x<-x[1:nlatent]
        # kpars[kpars$calcindices] <- unlist(lapply(kpars$calcs,function(x) eval(parse(text=x))))
      }
      
      if(!is.null(kpars$calcs)) {
        for(calci in kpars$calcs){
        eval(parse(text=calci))
      }}
      
        discreteDRIFT <- expm(kpars$DRIFT * dt[rowi])
        discreteCINT <- solve(kpars$DRIFT, (discreteDRIFT - diag(nlatent))) %*% kpars$CINT
      out <- discreteCINT  + discreteDRIFT %*% x
      if(ntdpred > 0) out <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
      if(jacobian) out <- c(out,PARMEANS)
      return(out)
    }
    
    PARMEANS <- as.numeric(kpars$PARMEANS)
    nvarpar <- length(PARMEANS)
    
    
    
      for(rowi in 1:(nrow(datalong))){
        if(continuoustime){
    
          if(t0check[rowi]==1 || dt[rowi]!=dt[rowi-1]){
            discreteDRIFT[[rowi]] <- expm(kpars$DRIFT * dt[rowi])
            discreteCINT[[rowi]] <- solve(kpars$DRIFT, (discreteDRIFT[[rowi]] - diag(nlatent))) %*% kpars$CINT
            discreteDIFFUSION[[rowi]] <- asymDIFFUSION - (discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE] %*%
                asymDIFFUSION %*% t(discreteDRIFT[[rowi]][diffusionindices,diffusionindices,drop=FALSE]))
          } else {
            discreteDRIFT[[rowi]] <-  discreteDRIFT[[rowi-1]]
            discreteCINT[[rowi]] <- discreteCINT[[rowi-1]]
            discreteDIFFUSION[[rowi]] <-  discreteDIFFUSION[[rowi-1]]
          }
        }
        
        if(!continuoustime){
          discreteDRIFT[[rowi]] <-kpars$DRIFT
          discreteCINT[[rowi]]<- kpars$CINT
          discreteDIFFUSION[[rowi]] <- kpars$DIFFUSION[diffusionindices,diffusionindices,drop=FALSE]
        }
        
        if(t0check[rowi]==1){
          etaprior[[rowi]]<-kpars$T0MEANS 
          if(ntdpred > 0) etaprior[[rowi]] <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
          etapriorcov[[rowi]]<-kpars$T0VAR[diffusionindices,diffusionindices,drop=FALSE]
        }
        
        if(t0check[rowi]==0){

          etaprior[[rowi]] <- dynamicmodel(c(etaupd[[rowi-1]],PARMEANS))


            if(is.null(kpars$PARVAR)) {
              dF <- numDeriv::jacobian(func = dynamicmodel, x= c(etaprior[[rowi]]),method='simple',
                method.args=list(eps=sqrt(c(diag(discreteDIFFUSION[[rowi]])))),jacobian=FALSE)
              etapriorcov[[rowi]] <-  dF %*% etaupdcov[[rowi-1]] %*% t(dF)  + discreteDIFFUSION[[rowi]]
            }
              
            if(!is.null(kpars$PARVAR)) {
              dF <- numDeriv::jacobian(func = dynamicmodel, x= c(etaprior[[rowi]],PARMEANS),method='simple',
                method.args=list(eps=sqrt(c(diag(discreteDIFFUSION[[rowi]]),diag(kpars$PARVAR)))),jacobian=TRUE)
             etaupdcov <- cbind( rbind(etaupdcov[[rowi-1]],rep(0,nvarpar)), rbind(rep(0,nlatent),kpars$PARVAR))
             discreteDIFF <- cbind( rbind(discreteDIFFUSION[[rowi]],rep(0,nvarpar)), matrix(0,nrow=nlatent+nvarpar))
              etapriorcov[[rowi]] <- (dF %*% etaupdcov %*% t(dF)  + discreteDIFF )[1:nlatent,1:nlatent]
            }
        }

        if(imputeMissings) Y[rowi,] <- 0 #etaprior[[rowi]] + t(chol(etapriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
        
        nafilter<-!is.na(Y[rowi,])
        observed[[rowi]]<-nafilter
        err[[rowi]]<-matrix(NA,nrow=nmanifest) #init prediction errors for this row
        
        # linked LAMBDA
        # H0 <- 1/(1+exp(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]])))
        # H1 <- 1/(1+exp(-kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]] + etapriorcov[[rowi]])))
        # H2 <- 1/(1+exp(-kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]] - etapriorcov[[rowi]])))
        # H <- (H1 + H2) /2
        # H <- t(solve(t(etaprior[[rowi]]), t(H)))
        
        # H0 <- exp(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]]))
        # H1 <- exp(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]] + etapriorcov[[rowi]]))
        # H2 <- exp(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE] %*% (etaprior[[rowi]] - etapriorcov[[rowi]]))
        # H <- (H1 + H2) /2
        # H <- t(solve(t(etaprior[[rowi]]), t(H)))

        H <- numDeriv::jacobian(func = measurementmodel, etaprior[[rowi]] ,
          method='simple',method.args=list(eps= sqrt(diag(kpars$DIFFUSION))))
        
        # LAMBDA <- kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]
        
        # // one step ahead predictive distribution of y
        yprior[[rowi]] <- kpars$MANIFESTMEANS + measurementmodel(etaprior[[rowi]])
        # yprior[[rowi]] <- ( kpars$MANIFESTMEANS + measurementmodel(etaprior[[rowi]]+sqrt(diag(kpars$DIFFUSION))) +
        #   kpars$MANIFESTMEANS + measurementmodel(etaprior[[rowi]]+sqrt(diag(kpars$DIFFUSION))) ) /2
        
        if(binary) kpars$MANIFESTVAR[row(kpars$MANIFESTVAR)==col(kpars$MANIFESTVAR)] <- .5^2 - (yprior[[rowi]]-.5)^2
        ypriorcov[[rowi]] <- H %*% etapriorcov[[rowi]] %*% t(H) + kpars$MANIFESTVAR
  
        # if(imputeMissings) Y[rowi,] <- yprior[[rowi]] + t(chol(ypriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
        
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
          
          
          if(binary) loglik[rowi] <-  sum(log( y*(yprior[[rowi]][nafilter,]) + (1-y)*(1-yprior[[rowi]][nafilter,]))) #sum over log is correct, yes?
          
          
          
          err[[rowi]][nafilter] <- as.numeric(y - yprior[[rowi]][nafilter,])
          
          K<-matrix(0,nrow=nlatent,ncol=nmanifest)
          
          K[diffusionindices,nafilter] <-  t(solve(t(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]),
            t(etapriorcov[[rowi]] %*%  t(H) ) ))
          
          # updated distribution 
          etaupd[[rowi]] <- etaprior[[rowi]] + K[,nafilter,drop=FALSE] %*% (err[[rowi]][nafilter,,drop=FALSE])
          
          # etaupdcov[[rowi + 1]] <- etapriorcov[[rowi]]
          etaupdcov[[rowi]] <- (diag(ndiffusion) - K[diffusionindices,nafilter,drop=FALSE] %*% H) %*% etapriorcov[[rowi]]
          
          # // log likelihood
          if(!binary) loglik[rowi] <- - 0.5 * (nrow(kpars$LAMBDA[nafilter,diffusionindices,drop=FALSE]) * log(2 * pi)  +
              log(det(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]))    +
              t(err[[rowi]][nafilter,,drop=FALSE]) %*%
              solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE],err[[rowi]][nafilter,,drop=FALSE]))
          
          if(!optimize){
            yupd[[rowi]] <- kpars$MANIFESTMEANS + H %*% etaupd[[rowi]]
            yupdcov[[rowi]] <- H %*% etaupdcov[[rowi]] %*% t(H) + kpars$MANIFESTVAR
          }
          
        }
        
      }
      
      # }
    }
    if(plotoptim && runif(1) > .95){
      yprior2 <- array(unlist(yprior),dim=c(length(yprior[[1]]),1,nrow(Y)))
      plot(apply(Y,1,mean),type='l')
      # points(1 / (1+exp(-yprior2[1,2,])),type='b',col='blue')
      points(yprior2[1,1,],type='l',col='red')
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
  
  
  
