#' ctGenerate
#' 
#' This function generates data according to the specified ctsem model object. 
#' 
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dtmean Positive numeric. Average time interval (delta T) to use.
#' @param logdtsd Numeric. Standard deviation for variability of the time interval.
#' @param wide Logical. Output in wide format?
#' @param simultdpredeffect logical - whether time dependent predictors impact 
#' instantaneously, or an instant *after* instantaneously. 
#' Switch reflects difference between ctStanFit and ctFit.
#' @details TRAITTDPREDCOV and TIPREDCOV matrices are not accounted for, at present. 
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

ctGenerate<-function(ctmodelobj,n.subjects=1000,burnin=0,dtmean=1,logdtsd=0,wide=TRUE,simultdpredeffect=FALSE){
  
  
  ###read in model
  for(i in 1:length(ctmodelobj)){ #this loop reads in the specified continuous time model
    assign(names(ctmodelobj[i]),ctmodelobj[[i]])

    if(is.matrix(ctmodelobj[[i]])){ #if this element is a matrix, continue on...
     
    if(any(is.na(suppressWarnings(as.numeric(get(names(ctmodelobj[i]))))))){ #if it contains character labels
      assign(names(ctmodelobj[i]),matrix(0,nrow=nrow(get(names(ctmodelobj[i]))), #set the values to 0 instead
                                      ncol=ncol(get(names(ctmodelobj[i])))))
      message(paste0(names(ctmodelobj[i])," contained character labels - setting matrix to 0"))
    }
    
     #set any matrices to numeric elements
      assign(names(ctmodelobj[i]), 
        matrix(as.numeric(get(names(ctmodelobj[i]))),nrow=nrow(get(names(ctmodelobj[i]))),
                                      ncol=ncol(get(names(ctmodelobj[i])))))
    }
  }

  
  #lower triangular transform to full and keep cholesky (named appropriately) 
  for(tempmatname in c('T0VAR','MANIFESTVAR', 'DIFFUSION', 'TRAITVAR','MANIFESTTRAITVAR','TDPREDVAR','TIPREDVAR')){

    tryCatch(assign(paste0(tempmatname,'chol'),get(tempmatname) ), error=function(e) {
      assign(tempmatname,NULL)})
    
    tryCatch(assign(tempmatname,get(tempmatname) %*% t(get(tempmatname))), error=function(e) {
      assign(tempmatname,NULL)})
  }

  
  kpars=list(DRIFT=DRIFT,T0VAR=T0VAR,DIFFUSION=DIFFUSION,CINT=CINT,T0MEANS=T0MEANS,
    MANIFESTMEANS=MANIFESTMEANS,MANIFESTVAR=MANIFESTVAR, LAMBDA=LAMBDA)
  
  if(!is.null(TDPREDEFFECT)) kpars$TDPREDEFFECT<-TDPREDEFFECT

  
  datalong <- matrix(NA,nrow=n.subjects*Tpoints,ncol= 1 + 1 + n.manifest + n.TDpred + n.TIpred)
  colnames(datalong) <- c('id','time',manifestNames,TDpredNames,TIpredNames)
  
  fullTpoints<-burnin+Tpoints
  
  if(n.TDpred > 0) {
    TDPREDMEANS <- rbind(matrix(0,nrow=(burnin+ifelse(simultdpredeffect,0,1))*n.TDpred),
    TDPREDMEANS)
    if(simultdpredeffect) TDPREDMEANS=rbind(TDPREDMEANS,0)
  }
  
  # TRAITVARchol = t(chol(-solve(DRIFT) %*% TRAITVAR %*% -t(solve(DRIFT))))
  
  for(si in 1:n.subjects){
    dt<- exp(rnorm(fullTpoints,log(dtmean),logdtsd))
    time=rep(0,fullTpoints)
    for(t in 2:fullTpoints) time[t] = round(time[t-1] + dt[t],3)
    if(n.TDpred > 0) tdpreds <- matrix(sapply(1:(fullTpoints),function(x) 
      TDPREDMEANS[(1+(x-1)*n.TDpred):(x*n.TDpred)] + TDPREDVARchol[1:n.TDpred,1:n.TDpred] %*% rnorm(n.TDpred,0,1)),byrow=TRUE, ncol=n.TDpred)
    
    skpars=kpars
    if(any(TRAITVARchol != 0)) {
      traits = rnorm(n.latent,0,1)
      skpars$CINT = skpars$CINT + TRAITVARchol %*% traits
      skpars$T0MEANS = skpars$T0MEANS + T0TRAITEFFECT %*% traits
    }
    
    if(n.TIpred > 0) {
      tipreds <- TIPREDMEANS + TIPREDVARchol %*% rnorm(n.TIpred,0,1)
      skpars$CINT = skpars$CINT + TIPREDEFFECT %*% tipreds
    }
    
    manifests<-matrix(NA,fullTpoints,n.manifest)
    
    sdat <- cbind(si,time,manifests,
      if(n.TDpred > 0) tdpreds,
      if(n.TIpred > 0) matrix(tipreds,byrow=TRUE,nrow=fullTpoints,ncol=n.TIpred))
    
    colnames(sdat) <- colnames(datalong)
    
    sdat[,manifestNames] <- ctKalman(kpars=skpars,datalong=sdat,manifestNames=manifestNames,TDpredNames=TDpredNames,
      latentNames=latentNames,
      imputeMissings=TRUE)$y
    
    sdat=sdat[(burnin+1):fullTpoints,]
    
    sdat[,'time'] = sdat[,'time'] - sdat[1,'time']
    
    datalong[(1+(si-1)*Tpoints):(si*Tpoints),]<-sdat
  }
  

  if(wide==FALSE) return(datalong) else {
    datawide <- ctLongToWide(datalong = datalong,id = 'id',time = 'time',
      manifestNames = manifestNames, TDpredNames = TDpredNames,TIpredNames = TIpredNames)
    
    datawide <- ctIntervalise(datawide = datawide,Tpoints = Tpoints,n.manifest = n.manifest,n.TDpred = n.TDpred,n.TIpred = n.TIpred,
      manifestNames=manifestNames,TDpredNames=TDpredNames,TIpredNames=TIpredNames)
    return(datawide)
  }
}
