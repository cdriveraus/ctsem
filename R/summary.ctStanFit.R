ctStanRawSamples<-function(fit){
  if(length(fit$stanfit$stanfit@sim)==0) {
    samples = fit$stanfit$rawposterior
  } else {
    samples = t(stan_unconstrainsamples(fit$stanfit$stanfit,fit$standata))
  }
  return(samples)
}


#' Extract an array of subject specific parameters from a ctStanFit object.
#'
#' @param fit fit object
#' @param pointest if TRUE, returns only the set of individual difference parameters
#' based on the max a posteriori estimate (or the median if sampling approaches were used).
#' @param cores Number of cores to use.
#' @param nsamples Number of samples to calculate parameters for. Not used if pointest=TRUE.
#' 
#' @details This function returns the estimates of individual parameters, taking into account any
#' covariates and random effects. 
#'
#' @return an nsamples by nsubjects by nparams array.
#' @export
#'
#' @examples
#' indpars <- ctStanSubjectPars(ctstantestfit)
#' dimnames(indpars)
#' plot(indpars[1,,'cint1'],indpars[1,,'cint2'])
ctStanSubjectPars <- function(fit,pointest=TRUE,cores=2,nsamples='all'){
  
  if(!nsamples[1] %in% 'all') fit$stanfit$rawposterior <- 
      fit$stanfit$rawposterior[sample(1:nrow(fit$stanfit$rawposterior),nsamples),,drop=FALSE]
  pnames <- getparnames(fit,subjvariationonly = TRUE)
  if(length(pnames)==0) stop('No individually varying parameters in model!')
  m <- fit$ctstanmodelbase$pars
  if(pointest) tfp <- fit$stanfit$transformedparsfull else {
    gc()
    tfp <- ctExtract(fit,subjectMatrices = TRUE,cores=cores,nsamples=nsamples)
  }
  p <- array(NA, dim=c(dim(tfp$pop_DRIFT)[1],fit$standata$nsubjects,length(pnames)))
  dimnames(p) <- list(iter=1:dim(p)[1],subject=1:dim(p)[2],param=1:dim(p)[3])
  co<-0
  for(i in 1:nrow(m)){
    if(!is.na(m$param[i]) & is.na(m$value[i]) & #if a free param, not complex and not already collected...
        !grepl('[',m$param[i],fixed=TRUE) & !m$param[i] %in% dimnames(p)[[3]]){
      if(m$param[i] %in% pnames){ #if ind differences
        co = co+1
        p[,,co] <- tfp[[paste0('subj_',m$matrix[i])]] [
          , ,m$row[i],m$col[i]]
        dimnames(p)[[3]][co] <- m$param[i]
      }
    }
  }
  p=p[,,order(dimnames(p)[[3]]),drop=FALSE]
  return(p)
}


getparnames <- function(fit,reonly=FALSE, subjvariationonly=FALSE, popstatesonly=FALSE){
  ms <- fit$setup$matsetup
  
  if(popstatesonly)  indices=ms$param > 0 & ms$copyrow <1 & ms$matrix==1 & ms$indvarying > 0 & ms$row > fit$standata$nlatent
  if(!popstatesonly)  indices=ms$when %in% c(0,-1) & ms$param > 0 & ms$copyrow < 1
  if(subjvariationonly) indices = ms$when %in% c(0,-1) & ms$param > 0 & ms$copyrow < 1 & (ms$tipred >0 | ms$indvarying > 0)
  pars <- data.frame(parnames = ms$parname[indices],  parindices = ms$param[indices])
  
  pars<-pars[!duplicated(pars$parnames),]
  pars<-pars[order(pars$parindices),]
  parnames <- pars[pars$parindices >0, 1]
  if(reonly)  parnames <- parnames[fit$standata$indvaryingindex]
  
  return(parnames)
}

#' summary.ctStanFit
#'
#' Summarise a ctStanFit object that was fit using \code{\link{ctStanFit}}. 
#' 
#' @param object fit object from \code{\link{ctStanFit}}, of class ctStanFit.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter calculations
#' reported in summary. 
#' @param digits integer denoting number of digits to report.
#' @param parmatrices if TRUE, also return additional parameter matrices -- can be slow to compute
#' for large models with many samples.
#' @param priorcheck Whether or not to use \code{ctsem:::priorchecking} to compare posterior mean and sd to prior mean and sd.
#' @param residualcov Whether or not to show standardised residual covariance. Takes a little longer to compute.
#' @param ... Additional arguments to pass to \code{ctsem:::priorcheckreport}, such as \code{meanlim}, or \code{sdlim}.
#' @return List containing summary items.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=4,parmatrices=TRUE,priorcheck=TRUE,residualcov = TRUE,...){
  
  if(!'ctStanFit' %in% class(object)) stop('Not a ctStanFit object!')
  
  monitor <- function(dat,...){
    out <- rstan::monitor(dat,...)
    class(out) <- 'data.frame'
    return(out)
  }
  
  
  out=list()
  monvars <- c('mean','sd','2.5%','50%','97.5%')
  
  optimize=length(object$stanfit$stanfit@sim)==0
  
  if(!optimize){ 
    smr<-suppressWarnings(getMethod('summary','stanfit')(object$stanfit$stanfit))
    if('98%' %in% colnames(smr$summary)) colnames(smr$summary)[colnames(smr$summary)=='98%'] <- '97.5%'
  }
  
  
  e <- ctExtract(object) 
  
  if(residualcov){ #cov of residuals
    obscov <- cov(object$data$Y,use='pairwise.complete.obs')
    idobscov <- diag(1/sqrt(diag(obscov)),ncol(obscov))
    rescov <- cov(matrix(object$stanfit$kalman$errprior,ncol=ncol(obscov)),use='pairwise.complete.obs')
    narescov <- which(is.na(rescov))
    rescov[narescov] <- 0
    
    out$residCovStd <- round(idobscov %*% rescov %*% idobscov ,3)
    out$residCovStd[narescov] <- NA
    dimnames(out$residCovStd) <- list(object$ctstanmodel$manifestNames,object$ctstanmodel$manifestNames)
    out$resiCovStdNote <- 'Standardised covariance of residuals'
  }
  
  ms=object$setup$matsetup
  parnames = getparnames(object)
  # parnames <- unique(parnames)
  parnamesiv <- getparnames(object,reonly = TRUE)
  
  #### generate covcor matrices of raw and transformed subject level params
  
  iter=dim(e$rawpopcov)[1]
  if(!is.null(iter)){ #then there is some individual variation so continue
    nindvarying=dim(e$rawpopcov)[2]
    
    if(nindvarying>1){
      
      #raw pop distribution params
      dimrawpopcorr <- dim(e$rawpopcorr)
      rawpopcorr= array(e$rawpopcorr,dim=c(dimrawpopcorr[1],1,dimrawpopcorr[2] * dimrawpopcorr[3]))
      rawpopcorrout <- suppressWarnings(monitor(rawpopcorr, digits_summary=digits,warmup=0,print = FALSE)[lower.tri(diag(nindvarying)),c(monvars,'n_eff','Rhat'),drop=FALSE])
      if(optimize) rawpopcorrout <- rawpopcorrout[,-which(colnames(rawpopcorrout) %in% c('n_eff','Rhat')),drop=FALSE]
      
      rownames(rawpopcorrout) <- matrix(paste0('',parnamesiv,'__',rep(parnamesiv,each=length(parnamesiv))),
        length(parnamesiv),length(parnamesiv))[lower.tri(diag(nindvarying)),drop=FALSE]
      
      rawpopcorrout <- cbind(rawpopcorrout,rawpopcorrout[,'mean'] / rawpopcorrout[,'sd'])
      colnames(rawpopcorrout)[ncol(rawpopcorrout)] <- 'z'
      
      out$rawpopcorr = round(rawpopcorrout,digits)
    }
  }
  
  if(priorcheck & !object$standata$priors) {
    priorcheckres <- priorcheckreport(object,...)
    if(nrow(priorcheckres$priorcheck) > 0) out = c(out,priorcheckres)
  }
  
  if(object$ctstanmodel$n.TIpred > 0) {
    
    if(!optimize){
      rawtieffect <- rstan::extract(object$stanfit$stanfit,permuted=FALSE,pars='TIPREDEFFECT')
      tidiags <- suppressWarnings(monitor(rawtieffect,warmup=0,digits_summary = digits,print = FALSE))
    }
    
    tieffect <- array(e$linearTIPREDEFFECT,dim=c(dim(e$linearTIPREDEFFECT)[1], 1, length(parnames) * dim(e$linearTIPREDEFFECT)[3]))
    tieffectnames <- paste0('tip_',rep(object$ctstanmodel$TIpredNames,each=length(parnames)),'_',parnames)
    dimnames(tieffect)<-list(c(),c(),tieffectnames)
    tipreds = suppressWarnings(monitor(tieffect,warmup = 0,print = FALSE)[,monvars,drop=FALSE])
    if(!optimize) tipreds <- cbind(tipreds,tidiags[,c('n_eff','Rhat'),drop=FALSE])
    tipreds <- tipreds[c(object$data$TIPREDEFFECTsetup)>0,,drop=FALSE]
    z = tipreds[,'mean'] / tipreds[,'sd'] 
    out$tipreds= round(cbind(tipreds,z),digits) #[order(abs(z)),]
  }
  
  
  
  if(parmatrices){
    Mean=ctStanContinuousPars(object,calcfunc = mean,calcfuncargs=list(),timeinterval=timeinterval)
    sd=ctStanContinuousPars(object,calcfunc = sd,calcfuncargs = list(na.rm=TRUE),timeinterval=timeinterval)
    `2.5%` = ctStanContinuousPars(object,calcfunc = quantile,calcfuncargs = list(probs=.025),timeinterval=timeinterval)
    `50%` = ctStanContinuousPars(object,calcfunc = quantile,calcfuncargs = list(probs=.5),timeinterval=timeinterval)
    `97.5%` = ctStanContinuousPars(object,calcfunc = quantile,calcfuncargs = list(probs=.975),timeinterval=timeinterval)
    
    d <- data.frame(ctModelUnlist(Mean,matnames = names(Mean)))
    colnames(d)[colnames(d) %in% 'value'] <- 'Mean'
    lapply(c('sd','2.5%','50%','97.5%'),function(x){
      d[[x]] <<-round(ctModelUnlist(get(x),names(Mean))$value,digits)
      })
    d$param <- NULL
    d$Mean <- round(d$Mean,digits)
    rm(sd)
    d <- d[!d$matrix %in% c('DIFFUSION','T0VAR'),]
    
    out$parmatrices=d
  }
  
  
  
  if(!optimize){
    # browser()
    popsd=smr$summary[c(grep('^popsd',rownames(smr$summary),fixed=FALSE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE] #[ object$data$indvaryingindex,,drop=FALSE]
    rownames(popsd)=parnames[ object$data$indvaryingindex]
    popmeans=smr$summary[c(grep('popmeans[', rownames(smr$summary),fixed=TRUE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=popmeans[(nrow(popmeans)/2+1):nrow(popmeans),,drop=FALSE]
    rownames(popmeans) <- parnames
    
    loglik = data.frame(mean=mean(e$ll),sd=sd(e$ll), max=max(e$ll))
    logposterior=smr$summary[c(grep('lp',rownames(smr$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
  }
  
  if(optimize){ #if optimized / importance sampled
    
    if(!is.null(iter)){ popsd <- suppressWarnings(monitor(array(
      e$popsd,dim=c(dim(e$popsd)[1],1,dim(e$popsd)[2])),warmup=0,print=FALSE))
    popsd=popsd[, monvars,drop=FALSE]
    rownames(popsd)=parnamesiv
    }
    popmeans=suppressWarnings(monitor(array(e$popmeans,dim=c(dim(e$popmeans)[1],1,dim(e$popmeans)[2])),warmup=0,print=FALSE))
    rownames(popmeans) = parnames #names(e)[grep('hmean_',names(e))]
    popmeans = popmeans[,monvars,drop=FALSE]
    
    loglik = object$stanfit$transformedparsfull$ll
    logposterior = object$stanfit$optimfit$value
    npars = length(object$stanfit$rawest)
    aic = 2* npars - 2*loglik
  }
  
  if(!is.null(iter)) out$popsd=round(popsd,digits=digits)
  
  out$popmeans=round(popmeans,digits=digits)
  
  out$popNote=paste0('popmeans are reported as specified in ctModel -- covariance related matrices are in sd / unconstrained correlation form -- see $parmatrices for simpler interpretations!')
  
  if(optimize) {
    
    out$loglik=loglik
    out$npars = npars
    out$aic = aic
  }
  out$logposterior=logposterior
  if(optimize) out$nsamples <- nrow(object$stanfit$samples)
  
  if(!parmatrices) out$parmatNote <- 'For additional summary matrices, use argument: parmatrices = TRUE'
  
  out <- lapply(out,function(x){
    if('matrix' %in% class(x)){
      x <- data.frame(x,check.names=FALSE)
    }
    x})
  
  return(out)
}
