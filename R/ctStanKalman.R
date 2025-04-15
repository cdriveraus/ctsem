#' Get Kalman filter estimates from a ctStanFit object
#'
#' @param fit fit object from \code{\link{ctStanFit}}.
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
#' k=ctStanKalman(ctstantestfit,subjectpars=TRUE,collapsefunc=mean)
ctStanKalman <- function(fit,nsamples=NA,pointest=TRUE, collapsefunc=NA,cores=1, 
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



