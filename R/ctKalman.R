#' Get Kalman filter estimates from a ctStanFit object
#'
#' @param fit fit object from \code{\link{ctFit}}.
#' @param nsamples either NA (to extract all) or a positive integer from 1 to maximum samples in the fit.
#' @param cores Integer number of cpu cores to use. Only needed if savescores was set to FALSE when fitting.
#' @param collapsefunc function to apply over samples, such as \code{mean}
#' @param pointest If TRUE, uses the posterior mode as the single sample.
#' @param standardisederrors If TRUE, computes standardised errors for prior, upd, smooth conditions.
#' @param subjectpars if TRUE, state estimates are not returned, instead, predictions of each subjects parameters
#' are returned, for parameters that had random effects specified.
#' @param tformsubjectpars if FALSE, subject level parameters are returned in raw, pre transformation form.
#' @param indvarstates if TRUE, do not remove indvarying states from output
#' @param removeObs Logical or integer. If TRUE, observations (but not covariates)
#' are set to NA, so only expectations based on parameters and covariates are returned. If a positive integer N, 
#' every N observations are retained while others are set NA for computing model expectations -- useful for observing prediction performance
#' forward further in time than one observation.
#' @param subjects integer vector of subjects to compute for.
#' @param timestep Either a positive numeric value, 'asdata' to use the times in the dataset, or 'auto' to select 
#' a timestep automatically (resulting in some interpolation but not excessive computation).
#' @param maxtime only relevant if timestep is not 'asdata'. Positive numeric denoting max time for computations.
#' @param ... additional arguments to collpsefunc.
#'
#' @return list containing Kalman filter elements, each element in array of
#' iterations, data row, variables. llrow is the log likelihood for each row of data.
#' @export
#'
#' @examples 
#' k=ctKalmanArray(ctstantestfit,subjectpars=TRUE,collapsefunc=mean)
ctKalmanArray <- function(fit,nsamples=NA,pointest=TRUE, collapsefunc=NA,cores=1,
  subjects=1:max(fit$standata$subject), timestep='asdata',maxtime='asdata',
  standardisederrors=FALSE, subjectpars=TRUE, tformsubjectpars=TRUE, indvarstates=FALSE,removeObs=F,...){
  
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
  if(fit$standata$intoverstates==0){
    warning('Kalman filter operation unreliable when states were sampled -- system noise represents prior while point estimates represent posterior / smoothed')
  }
  message('Computing state estimates..')
  # standata <- fit$standata
  if(pointest)  samples <- matrix(fit$stanfit$rawest,nrow=1)
  
  if(!pointest){
    samples<-ctStanRawSamples(fit)
    if(!is.na(nsamples)) samples <- samples[sample(1:nrow(samples),nsamples),,drop=FALSE] else nsamples <- nrow(samples)
    if(is.function(collapsefunc)) samples = matrix(apply(samples,2,collapsefunc,...),ncol=ncol(samples))
  }
  
  # use only selected subjects data -----------------------------------------
  if(length(subjects)!=fit$standata$nsubjects){
    idstore <- fit$standata$subject
    if(length(fit$stanfit$stanfit@sim)==0) { #only select subjects if optimized fit
      fit$standata <- standatact_specificsubjects(fit$standata, subjects = subjects)
      subjects <- sort(unique(fit$standata$subject))
    }
  }
  
  if(removeObs || removeObs > 0){
    if(is.numeric(removeObs)) skipn <- removeObs else skipn <- 1
    for(x in c('nobs_y','nbinary_y','ncont_y','whichobs_y','whichbinary_y','whichcont_y')){
      if(skipn > 1) {
        if(length(dim(fit$standata[[x]]))==2) fit$standata[[x]][(1:skipn)!=1,] <- 0L else fit$standata[[x]][(1:skipn)!=1] <-0L 
      }
      if(skipn==1) {
        if(length(dim(fit$standata[[x]]))==2) fit$standata[[x]][,] <- 0L else fit$standata[[x]][] <-0L 
      }
    }
  }
  
  # timerange ---------------------------------------------------------------
  
  if(timestep=='auto'){
    timediff <- diff(fit$standata$time)
    timediff <- timediff[timediff > 0]
    if(fit$standata$intoverstates==1) timestep=median(timediff)/5 else timestep ='asdata'
  }
  if(all(maxtime == 'asdata')) maxtime <- max(fit$standata$time[fit$standata$subject %in% subjects]) 
  if(is.na(timestep)) timestep=0
  
  if(!'asdata' %in% timestep && fit$ctstanmodel$continuoustime) {
    if(fit$ctstanmodel$continuoustime != TRUE) stop('Discrete time model fits must use timestep = "asdata"')
    times <- seq(min(fit$standata$time[fit$standata$subject %in% subjects]) ,maxtime,timestep)
    fit$standata <- standataFillTime(fit$standata,times,subject=subjects,maintainT0=TRUE)
  }
  
  #only do computations for rows of subjects requested, only really necessary for sampling approach because otherwise dropped from data
  fit$standata$dokalmanrows <- as.integer(fit$standata$subject %in% subjects)
  
  
  e=stan_constrainsamples(sm = fit$stanmodel,standata = fit$standata,
    savesubjectmatrices = subjectpars,
    samples = samples,cores=cores,savescores=TRUE,pcovn=5)
  
  
  e$ya <- e$ya[,,fit$standata$dokalmanrows==1,,drop=FALSE]
  e$ycova <- e$ycova[,,fit$standata$dokalmanrows==1,,,drop=FALSE]
  e$etaa <- e$etaa[,,fit$standata$dokalmanrows==1,,drop=FALSE]
  e$etacova <- e$etacova[,,fit$standata$dokalmanrows==1,,,drop=FALSE]
  
  nsamples <-nrow(samples) #in case it was set NA, compute nsamples
  e$yprior <- array(e$ya[,1,,,drop=FALSE],dim=dim(e$ya)[-2])
  e$yupd <-  array(e$ya[,2,,,drop=FALSE],dim=dim(e$ya)[-2])
  e$ysmooth<-  array(e$ya[,3,,,drop=FALSE],dim=dim(e$ya)[-2])
  e$etaprior <-  array(e$etaa[,1,,,drop=FALSE],dim=dim(e$etaa)[-2])
  e$etaupd <-  array(e$etaa[,2,,,drop=FALSE],dim=dim(e$etaa)[-2])
  e$etasmooth <-  array(e$etaa[,3,,,drop=FALSE],dim=dim(e$etaa)[-2])
  e$ypriorcov <-  array(e$ycova[,1,,,,drop=FALSE],dim=dim(e$ycova)[-2])
  e$yupdcov <-  array(e$ycova[,2,,,,drop=FALSE],dim=dim(e$ycova)[-2])
  e$ysmoothcov <-  array(e$ycova[,3,,,,drop=FALSE],dim=dim(e$ycova)[-2])
  e$etapriorcov <-  array(e$etacova[,1,,,,drop=FALSE],dim=dim(e$etacova)[-2])
  e$etaupdcov <-  array(e$etacova[,2,,,,drop=FALSE],dim=dim(e$etacova)[-2])
  e$etasmoothcov <-  array(e$etacova[,3,,,,drop=FALSE],dim=dim(e$etacova)[-2])
  
  
  
  nlatent <- ifelse(!indvarstates, fit$standata$nlatent,fit$standata$nlatentpop)
  latentNames <- fit$ctstanmodel$latentNames
  if(indvarstates) latentNames <- c(latentNames,
    # paste0('indvar',1:(fit$standata$nlatentpop-fit$standata$nlatent))
    getparnames(fit,popstatesonly=TRUE)
  )
  nmanifest <- fit$standata$nmanifest
  
  
  
  out=list(time=cbind(fit$standata$time[fit$standata$dokalmanrows==1]), 
    y=matrix(fit$standata$Y[fit$standata$dokalmanrows==1,,drop=FALSE],ncol=ncol(fit$standata$Y),dimnames = list(NULL,fit$ctstanmodel$manifestNames)), 
    llrow=e$llrow[,fit$standata$dokalmanrows==1,drop=FALSE])
  out$y[out$y==99999] <- NA
  for(basei in c('y','eta')){
    for(typei in c('prior','upd','smooth')){
      for(typex in c('','cov')){
        ref=paste0(basei,typei,typex)
        out[[ref]] <- e[[ref]]
        out[[ref]][out[[ref]] == 99999] <- NA
        if(basei=='y') {
          dimnames(out[[ref]]) <- list(NULL, NULL, fit$ctstanmodel$manifestNames) 
        } 
        if(basei=='eta'){
          if(typex=='') {
            out[[ref]] <- out[[ref]][,,1:nlatent,drop=FALSE] 
            dimnames(out[[ref]]) <- list(NULL, NULL, latentNames)
          } else { #for cov
            out[[ref]] <- out[[ref]][,,1:nlatent,1:nlatent,drop=FALSE]
          }
        }
      }
    }
  }
  
  for(typei in c('prior','upd','smooth')){
    out[[paste0('err',typei)]] <- aaply(out[[paste0('y',typei)]],1, function(yp) array(out$y-yp,dim=dim(out$y)),.drop=FALSE,.inform=TRUE)
  } 
  # 
  if(standardisederrors){
    for(typei in c('prior','upd','smooth')){
      arr <- array(sapply(1:dim(out$yprior)[1], function(i){
        array(sapply(1:nrow(out$y), function(r){
          tmp <- matrix(NA,nmanifest)
          if(sum(!is.na(out$y[r,])) > 0) tmp[which(!is.na(out$y[r,]))] <- 
              matrix(solve(
                t(chol(matrix(out[[paste0('ypriorcov')]][i,r,,],ncol=nmanifest) + diag(1e-10,nmanifest)))[
                  !is.na(out$y[r,]),!is.na(out$y[r,])], 
                out[[paste0('err',typei)]][i,r,!is.na(out$y[r,])]), nrow=sum(!is.na(out$y[r,])))
          return(tmp)
        },simplify = 'array'), dim=c(nmanifest,1,nrow(out$y)))
      },simplify = 'array'), dim=c(nmanifest,1,nrow(out$y),nsamples))
      
      out[[paste0('errstd',typei)]] <- array(aperm(arr, c(4,3,1,2)),dim=dim(arr)[c(4,3,1)])
    }
  }
  
  mindex <- grep('(^y)|(^err)|(^ll)',names(out))
  lindex <- grep('^eta',names(out))
  nosampindex <- which(names(out) %in% c('time','y'))
  out$llrow <- matrix(out$llrow,dim(out$llrow)[1],dim(out$llrow)[2])
  
  for(i in 1:length(out)){
    d<-list()
    if(!i %in% nosampindex){
      ds <- 1:dim(out[[i]])[1]
      d <- c(d,Sample=list(ds))
    }
    do <- 1:dim(out[[i]])[ifelse(i %in% nosampindex,1,2)]#obs
    d <- c(d,Obs=list(do))
    
    if(names(out)[i] %in% 'time') d <- c(d,Row=list('Time'))
    if(names(out)[i] %in% 'y') d <- c(d,Row = list(fit$ctstanmodelbase$manifestNames))
    
    
    if(length(dim(out[[i]])) > 2){
      if(i %in% mindex) dr <- fit$ctstanmodelbase$manifestNames
      if(i %in% lindex) dr <- latentNames
      d <- c(d,Row=list(dr))
      if(length(dim(out[[i]])) > 3) d <- c(d,Col=list(dr))
    }
    
    dimnames(out[[i]]) <- d
  }
  out$id <- fit$standata$subject[fit$standata$dokalmanrows==1]
  
  return(out)
}

#' Backward-compatible alias for \code{ctKalmanArray}.
#' @rdname ctKalmanArray
#' @export
ctStanKalman <- ctKalmanArray




#' ctPredictTIP
#' 
#' Outputs the estimated effect of time independent predictors (covariate moderators) on the expected observations.
#' 
#' @param sf A fitted ctStanFit object from the ctsem package.
#' @param tipreds A character vector specifying which time independent predictors to use. Default is 'all', which uses all time independent predictors in the model.
#' @param subject An integer value specifying the internal ctsem subject ID (mapping visible under myfit$setup$idmap) for which predictions are made. 
#' This is relevant only when time dependent predictors are also included in the model. 
#' @param doDynamics A logical value indicating whether to plot the effects of time independent predictors on the dynamics of the system. Default is TRUE. 
#' Can be problematic for systems with many dimensions.
#' @param timestep A numeric value specifying the time step for predictions. Default is 'auto', which tries to automatically determine an appropriate time step.
#' @param plot A logical value indicating whether to ggplot the results instead of returning a data.frame of predictions. Default is TRUE.
#' @param quantiles A numeric vector specifying the quantiles of the time independent predictors to plot. Default is 1SD either side and the median, c(.32,.5,.68).
#' @param discreteTimeQuantiles a numeric vector of length 3 specifying the quantiles of the discrete time points to plot, when 
#' showUncertainty is TRUE.
#' @param dynamicsControl A named list of additional arguments to pass to
#' \code{\link{ctDiscretePars}} and \code{\link{ctDiscreteParsPlot}}
#' when \code{doDynamics=TRUE}. Arguments matching
#' \code{ctDiscreteParsPlot} are automatically routed to the plot call.
#' Internally controlled arguments \code{fit}, \code{ctstanfitobj}, \code{plot},
#' \code{subjects}, \code{observational}, \code{x}, \code{quantiles}, and
#' \code{splitSubjects} are ignored if supplied.
#' @param showUncertainty A logical value indicating whether to plot the uncertainty of the predictions. Default is TRUE.
#' @param TIPvalues An nvalue * nTIpred numeric matrix specifying the fixed values for each time independent predictor effect to plot. 
#' Default is NA, which instead relies on the quantiles specified in the quantiles argument.
#' 
#' @return If plot is TRUE, a list of ggplot objects showing the estimated effects of covariate moderators. Otherwise, a data frame with the predictions.
#' 
#' @details This function estimates the effects of covariate moderators on the expected process 
#' and observations for a specified subject in a dynamic system. The covariate moderators are defined at the specified quantiles, 
#' and their effects on the trajectory are plotted or returned as a data frame.
#' 
#' @examples
#' # Example usage:
#' ctPredictTIP(ctstantestfit, tipreds='all', doDynamics=FALSE, plot=TRUE)
#' @export
ctPredictTIP <- function(sf,tipreds='all',subject=1,timestep='auto',doDynamics=TRUE, plot=TRUE,
  quantiles=c(.16,.5,.84), discreteTimeQuantiles=c(.025, .5, .975),
  dynamicsControl=list(),
  showUncertainty=TRUE, 
  TIPvalues=NA){
  if(!is.list(dynamicsControl)) stop('dynamicsControl must be a list')
  if(length(dynamicsControl) > 0 && (is.null(names(dynamicsControl)) || any(names(dynamicsControl) == ''))){
    stop('dynamicsControl must be a named list')
  }
  dynamicsControl[names(dynamicsControl) %in% c('fit','ctstanfitobj','plot','subjects','observational','x',
    'quantiles','splitSubjects')] <- NULL
  dynamicsPlotControl <- dynamicsControl[names(dynamicsControl) %in% names(formals(ctDiscreteParsPlot))]
  dynamicsControl <- dynamicsControl[!names(dynamicsControl) %in% names(dynamicsPlotControl)]
  if(tipreds[1] %in% 'all') tipreds <- sf$ctstanmodel$TIpredNames
  if(length(subject) > 1) stop('>1 subject!')
  if(length(unique(sf$standata$subject)) < 3) stop('With fewer than 3 subjects in the data, these predictions are not possible')
  
  if(all(is.na(TIPvalues))){
    TIPvalues = apply(sf$standata$tipredsdata[,tipreds,drop=FALSE],2,quantile,probs=quantiles)
  }
  #check for duplicates in columns of TIPvalues and stop if found
  if(!all(apply(TIPvalues,2,function(x) length(unique(x))==length(x)))){
    stop('Duplicate values for TI predictors -- if using categorical / dummy predictors, specify values using TIPvalues arg')
  }
  
  sdat <- standatact_specificsubjects(standata = sf$standata,subjects = subject)
  
  dat <- standatatolong(sdat,origstructure = TRUE,ctm=sf$ctstanmodelbase)
  dat[,sf$ctstanmodelbase$manifestNames] <- NA #set all manifest obs to missing
  
  TIPvalues <- matrix(apply(TIPvalues,2,function(x) sort(x)),ncol=ncol(TIPvalues)) #sort ascending to get plot colours correct
  
  fulldat <- dat[c(),]
  for(tipi in 1:length(tipreds)){ #for each tipred
    for(vali in 1:nrow(TIPvalues)){ #for each tipred value
      tdat <- dat #copy the data
      tdat[,sf$ctstanmodelbase$TIpredNames] <- 0 #set all covariates to zero
      tipvali <- TIPvalues[vali,tipi] #get the tipred value
      tdat[[sf$ctstanmodelbase$subjectIDname]] <- paste0(tipreds[tipi],' = ',round(tipvali,2)) #modify the subject ID
      tdat[,tipreds[tipi]] <- tipvali #set the tipred value to the quantile or value specified
      if(nrow(fulldat)==0){
        fulldat <- tdat #output to full dataset
      } else fulldat <- rbind(fulldat,tdat)
    }
  }
  
  sf$standata <- suppressMessages(ctStanData(sf$ctstanmodel,fulldat,optimize=TRUE))
  k=ctKalman(fit = sf,subjects=unique(fulldat[[sf$ctstanmodelbase$subjectIDname]]),realid=TRUE,timestep=timestep)
  k = k[k$Element %in% c('etaprior','yprior','ypriorcov','etapriorcov'),]
  k$V1 <- k$variable <- NULL

  
  if(!plot) return(k) else{
    Sample <- Col <- CovariateValue <- Effect <- Subject<- NULL #local variables for ggplot
    gglist <- list(Process=list(Observed=list(),Latent=list()),Dynamics=list(Independent=list(),Correlated=list())) #empty list
    for(tipi in 1:length(tipreds)){ #for each tipred
      ks <- k[grepl(paste0('^',tipreds[tipi],' = '),k$Subject),] #subset to tipred
      ks$Subject <- factor(as.numeric(gsub('^.* = ','',ks$Subject))) #extract the subject value

      
      
      for(elementi in c('yprior','etaprior')){
        ksp=plot.ctKalmanDF(ks,plot=FALSE,kalmanvec=elementi,
          polygonsteps = ifelse(showUncertainty,10,0))
        # if(tail(ksp$data$Subject,1)==99999){ #if we added a dummy subject, remove from plot now
        #   ksp$data <- ksp$data[-nrow(ksp$data),]
        #   # ksp$data$Subject <- factor(ksp$data$Subject)
        # }
        
        ksp <- ksp +
          theme(legend.position = 'bottom') +
          guides(colour=guide_legend(title=tipi)) +
          aes(linetype=NULL,colour=Subject,fill=Subject)+
          scale_colour_manual(values=colorRampPalette(c("red",'black', "blue"))(length(unique(ksp$data$Subject))))+
          scale_fill_manual(values=colorRampPalette(c("red",'black', "blue"))(length(unique(ksp$data$Subject))))+
          theme(legend.position = 'bottom')+
          guides(colour=guide_legend(title=tipreds[tipi]))+
          ggtitle(paste0('Effect of ',tipreds[tipi],' on ',ifelse(elementi=='yprior','observed','latent'),' trajectory'))+
          facet_wrap(vars(Variable),scales ='free_y')
        
        gglist[['Process']][[ifelse(elementi=='yprior','Observed','Latent')]][[tipreds[tipi]]] <- ksp
      }
      
      #include ctstandiscretepars plots
      if(doDynamics){
        for(typei in c('Independent','Correlated')){
          discreteParsArgs <- c(dynamicsControl, list(
            fit=sf,
            plot=FALSE,
            subjects=(tipi-1)*nrow(TIPvalues) + 1:nrow(TIPvalues),
            observational = !typei %in% 'Independent'))
          ctd=do.call(ctDiscretePars, discreteParsArgs)
          if(showUncertainty) ctdQuantiles <- discreteTimeQuantiles else ctdQuantiles <- c(.5,.5,.5)
          ctdp=do.call(ctDiscreteParsPlot, c(dynamicsPlotControl,
            list(x=ctd,quantiles=ctdQuantiles,splitSubjects = TRUE)))
          ctdp$data$CovariateValue <- factor(round(TIPvalues[ctdp$data$Subject,tipi],3))
          
          ctdp=ctdp+aes(colour=CovariateValue,fill=CovariateValue)+
            facet_wrap(vars(Effect))+theme(legend.position = 'bottom')+
            scale_colour_manual(values=colorRampPalette(c("red",'black', "blue"))(nrow(TIPvalues)))+
            scale_fill_manual(values=colorRampPalette(c("red",'black', "blue"))(nrow(TIPvalues)))+
            theme(legend.position = 'bottom') +
            guides(colour=guide_legend(title=tipreds[tipi]),
              fill=guide_legend(title=tipreds[tipi]))
          
          gglist[['Dynamics']][[typei]][[tipreds[tipi]]] <- ctdp
        }
      }
    }
    return(gglist)
  }
}

#' ctKalman 
#'
#' Outputs predicted, updated, and smoothed estimates of manifest indicators and latent states, 
#' with covariances, for specific subjects from data fit with \code{\link{ctFit}},
#' based on either the mode (if optimized) or mean (if sampled) of parameter distribution.
#' 
#' @param fit fit object as generated by \code{\link{ctFit}}.
#' @param timerange Either 'asdata' to just use the observed data range, or a numeric vector of length 2 denoting start and end of time range, 
#' allowing for estimates outside the range of observed data. Ranges smaller than the observed data are ignored.
#' @param timestep Either 'asdata' to just use the observed data 
#' (which also requires 'asdata' for timerange) or a positive numeric value
#' indicating the time step to use for interpolating values. Lower values give a more accurate / smooth representation,
#' but take a little more time to calculate. 
#' @param subjects vector of integers denoting which subjects (from 1 to N) to plot predictions for. 
#' @param removeObs Logical or integer. If TRUE, observations (but not covariates)
#' are set to NA, so only expectations based on parameters and covariates are returned. If a positive integer N, 
#' every N observations are retained while others are set NA for computing model expectations -- useful for observing prediction performance
#' forward further in time than one observation.
#' @param standardisederrors if TRUE, also include standardised error output (based on covariance
#' per time point).
#' @param plot Logical. If TRUE, plots output instead of returning it. 
#' See \code{\link{plot.ctKalmanDF}} 
#' (Stan based fit) for the possible arguments.
#' @param realid use original (not necessarily integer sequence) subject id's? Otherwise use integers 1:N.
#' @param ... additional arguments to pass to \code{\link{plot.ctKalmanDF}}.
#' @return Returns a list containing matrix objects etaprior, etaupd, etasmooth, y, yprior, 
#' yupd, ysmooth, prederror, time, loglik,  with values for each time point in each row. 
#' eta refers to latent states and y to manifest indicators - y itself is thus just 
#' the input data. 
#' Covariance matrices etapriorcov, etaupdcov, etasmoothcov, ypriorcov, yupdcov, ysmoothcov,  
#' are returned in a row * column * time array. 
#' Some outputs are unavailable for ctStan fits at present.
#' If plot=TRUE, nothing is returned but a plot is generated.
#' @examples
#' \donttest{
#' 
#' #Basic
#' ctKalman(ctstantestfit, timerange=c(0,60), plot=TRUE)
#' 
#' #Multiple subjects, y and yprior, showing plot arguments
#' plot1<-ctKalman(ctstantestfit, timerange=c(0,60), timestep=.1, plot=TRUE,
#'   subjects=2:3, 
#'   kalmanvec=c('y','yprior'),
#'   errorvec=c(NA,'ypriorcov')) #'auto' would also have achieved this
#'   
#'  #modify plot as per normal with ggplot
#'  print(plot1+ggplot2::coord_cartesian(xlim=c(0,10)))
#'  
#'  #or generate custom plot from scratch:#'  
#'  k=ctKalman(ctstantestfit, timerange=c(0,60), timestep=.1, subjects=2:3)
#'  library(ggplot2)
#'  ggplot(k[k$Element %in% 'yprior',],
#'    aes(x=Time, y=value,colour=Subject,linetype=Row)) +
#'    geom_line() +
#'    theme_bw()
#'
#'  }
#' @export

ctKalman<-function(fit, timerange='asdata', timestep='auto',
  subjects=fit$standata$idmap[1,1], removeObs = FALSE, plot=FALSE, 
  standardisederrors=FALSE,realid=TRUE,...){
  
  
  if('ctsemFit' %in% class(fit)) stop('This function is no longer supported with ctsemOMX, try ctsem')
  if(!'ctStanFit' %in% class(fit)) stop('fit object is not from ctStanFit!')
  
  # get subjects ------------------------------------------------------------
  idmap <- fit$standata$idmap #store now because we may reduce it
  if('factor' %in% class(idmap$original)) idmap$original <- as.character(idmap$original)
  if('factor' %in% class(subjects)) subjects <- as.character(subjects)
  subjectsarg <- subjects
  if(realid) subjects <- idmap[which(idmap[,1] %in% subjects),2]
  
  if(length(subjects) == 0){
    if(all(!is.na(as.integer(subjectsarg)))){ #if all subjects specified as integers
      subjects <- as.integer(subjectsarg)
      warning('Specified subjects not found in original id set -- assuming integers correspond to internal integer mapping. Consider setting realid=FALSE')
      realid=FALSE
    } else stop('Specified subjects not found in original id set, and (some) are not integers...')
  }
  subjects <- sort(subjects) #in case not entered in ascending order
  
  out <- ctKalmanArray(fit,pointest=TRUE,
    removeObs=removeObs, subjects=subjects,timestep = timestep,maxtime=max(timerange),
    collapsefunc=mean, indvarstates = FALSE,standardisederrors = standardisederrors) #extract state predictions
  
  out <- meltkalman(out)
  out[['Subject']] <- factor(subjects[out[['Subject']]]) #correct for subjects being set 1:Nsub by ctKalmanArray
  if(!all(timerange %in% 'asdata')) out <- out[out$Time >= min(timerange) & out$Time <= max(timerange),]
  if(realid){
    out$Subject <- factor(idmap[
      match(out$Subject,idmap[,2]),1])
  }
  
  # out=out[!(out$Subject %in% subjects) %in% FALSE,]
  
  
  if(plot) {
    plot(x=out,subjects=subjectsarg,...)
  } else return(out)
}




#' Plots Kalman filter output from ctKalman.
#'
#' @param x Output from \code{\link{ctKalman}}. In general it is easier to call 
#' \code{\link{ctKalman}} directly with the \code{plot=TRUE} argument, which calls this function.
#' @param subjects vector of integers denoting which subjects (from 1 to N) to plot predictions for. 
#' @param kalmanvec string vector of names of any elements of the output you wish to plot, 
#' the defaults of 'y' and 'ysmooth' plot the original data, 'y', 
#' and the estimates of the 'true' value of y given all data. Replacing 'y' by 'eta' will 
#' plot latent states instead (though 'eta' alone does not exist) and replacing 'smooth' 
#' with 'upd' or 'prior' respectively plots updated (conditional on all data up to current time point)
#' or prior (conditional on all previous data) estimates.
#' @param errorvec vector of names indicating which kalmanvec elements to plot uncertainty bands for. 
#' 'auto' plots all possible.
#' @param elementNames if NA, all relevant object elements are included -- e.g. if yprior is in the kalmanvec
#' argument, all manifest variables are plotted, and likewise for latent states if etasmooth was specified.
#' Alternatively, a character vector specifying the manifest and latent names to plot explicitly can be specified.
#' @param plot if FALSE, plots are not generated and the ggplot object is simply returned invisibly.
#' @param errormultiply Numeric denoting the multiplication factor of the std deviation of errorvec objects. 
#' Defaults to 1.96, for 95\% intervals.
#' @param polygonsteps Number of steps to use for uncertainty band shading. 
#' @param polygonalpha Numeric for the opacity of the uncertainty region.
#' @param facets when multiple subjects are included in multivariate plots, the default is to facet plots 
#' by variable type. This can be set to NA for no facets, or "Subject" for facetting by subject.
#' @param ... not used.
#' @return A ggplot2 object. Side effect -- Generates plots.
#' @method plot ctKalmanDF
#' @export plot.ctKalmanDF
#' @export
#' @examples
#' ### Get output from ctKalman
#' x<-ctKalman(ctstantestfit,subjects=2,timestep=.01)
#' 
#' ### Plot with plot.ctKalmanDF
#' plot(x, subjects=2)
#' 
#' ###Single step procedure:
#' ctKalman(ctstantestfit,subjects=2,
#'   kalmanvec=c('y','yprior'),
#'   elementNames=c('Y1','Y2'), 
#'   plot=TRUE,timestep=.01)
plot.ctKalmanDF<-function(x, subjects=unique(x$Subject), kalmanvec=c('y','yprior'),
  errorvec='auto', errormultiply=1.96,plot=TRUE,elementNames=NA,
  polygonsteps=10,polygonalpha=.1,
  facets='Variable',
  ...){
  if(!'ctKalmanDF' %in% class(x)) stop('not a ctKalmanDF object')
  
  if(FALSE) Time <- Value <- Subject <- Row <- Variable <- Element <- NULL
  colnames(x)[colnames(x) %in% 'Row'] <- "Variable"
  colnames(x)[colnames(x) %in% 'value'] <- "Value"
  
  if(any(!is.na(elementNames))) x <- subset(x,Variable %in% elementNames)
  x <- subset(x,Subject %in% subjects)
  x<-subset(x,Element %in% kalmanvec)
  
  klines <- kalmanvec[grep('(prior)|(upd)|(smooth)',kalmanvec)]
  if(all(errorvec %in% 'auto')) errorvec <- klines
  errorvec <- errorvec[grep('(prior)|(upd)|(smooth)',errorvec)]
  colvec='Variable'
  # if(length(subjects) > 1){
  colvec=ifelse(facets %in% 'Subject','Variable','Subject')
  # }
  ltyvec <- setNames( rep(0,length(kalmanvec)),kalmanvec)
  ltyvec[kalmanvec %in% klines] = setNames(1:length(klines),klines)
  
  shapevec<-ltyvec
  shapevec[shapevec %in% 0] <- 1
  shapevec[shapevec>0] <- NA
  shapevec[ltyvec %in% 0] <- 19
  
  g <- ggplot(x,
    aes(x=Time,y=Value,colour=.data[[colvec]],linetype=Element,shape=Element)) +
    scale_linetype_manual(breaks=names(ltyvec),values=ltyvec)+
    guides(linetype='none',shape='none')+
    scale_shape_manual(breaks=names(shapevec),values=shapevec)
  
  if(length(unique(ltyvec[!is.na(ltyvec)]))<1) g<-g+guides(linetype='none')
  
  if(!is.na(facets[1]) && length(unique(subset(x,Element %in% kalmanvec)$Variable)) > 1){
    g <- g+ facet_wrap(facets,scales = 'free') 
  } 
  g <- g + ylab(ifelse(length(unique(x$Variable))==1, unique(x$Variable),'Variable'))
  
  polygonsteps <- polygonsteps + 1
  polysteps <- seq(errormultiply,0,-errormultiply/(polygonsteps+1))[c(-polygonsteps+1)]
  if(any(!is.na(errorvec))){
    for(si in polysteps){
      
      d2 <- subset(x,Element %in% errorvec)
      d2$sd <- d2$sd *si
      
      if(facets %in% 'Subject'){
        g<- g+ 
          geom_ribbon(data=d2,aes(ymin=(Value-sd),x=Time,
            ymax=(Value+sd),fill=(Variable)),
            alpha=ifelse(si== polysteps[1],.05,polygonalpha/polygonsteps),
            inherit.aes = FALSE,linetype=0)
        if(si== polysteps[1]) g <- g + 
            geom_line(data=d2,aes(y=(Value-sd),colour=Variable),linetype='dotted',alpha=.4) + 
            geom_line(data=d2,aes(y=(Value+sd),colour=Variable),linetype='dotted',alpha=.4)
      } else {
        g <- g+ 
          geom_ribbon(data=d2,aes(ymin=(Value-sd),x=Time,
            ymax=(Value+sd),fill=(Subject)),inherit.aes = FALSE,
            alpha=polygonalpha/polygonsteps,linetype=0)
        if(si== polysteps[1]) g <- g + 
            geom_line(data=d2,aes(y=(Value-sd),colour=Subject),linetype='dotted',alpha=.7) + 
            geom_line(data=d2,aes(y=(Value+sd),colour=Subject),linetype='dotted',alpha=.7)
      }
    } 
  }
  g <- g + 
    geom_line()+
    geom_point()+
    theme_minimal()+
    guides(fill='none')+
    theme(legend.position = 'bottom')
  
  return(g)
  
}


