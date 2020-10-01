ctStanRawSamples<-function(fit){
  if(!class(fit$stanfit) %in% 'stanfit') {
    samples = fit$stanfit$rawposterior
  } else {
    samples = t(stan_unconstrainsamples(fit$stanfit,fit$standata))
  }
  return(samples)
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
#' if(w32chk()){
#'
#' summary(ctstantestfit)
#' }
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=4,parmatrices=TRUE,priorcheck=TRUE,residualcov = TRUE,...){
  
  if(!'ctStanFit' %in% class(object)) stop('Not a ctStanFit object!')
  
  out=list()
  monvars <- c('mean','sd','2.5%','50%','97.5%')
  
  if('stanfit' %in% class(object$stanfit)){ 
    smr<-suppressWarnings(getMethod('summary','stanfit')(object$stanfit))
    if('98%' %in% colnames(smr$summary)) colnames(smr$summary)[colnames(smr$summary)=='98%'] <- '97.5%'
  }

  e <- ctExtract(object) 
 
  if(residualcov){ #cov of residuals
  obscov <- cov(object$data$Y,use='pairwise.complete.obs')
  idobscov <- diag(1/sqrt(diag(obscov)),ncol(obscov))
  rescov <- cov(matrix(object$kalman$errprior,ncol=ncol(obscov)),use='pairwise.complete.obs')
  narescov <- which(is.na(rescov))
  rescov[narescov] <- 0
  
  out$residCovStd <- round(idobscov %*% rescov %*% idobscov ,3)
  out$residCovStd[narescov] <- NA
  dimnames(out$residCovStd) <- list(object$ctstanmodel$manifestNames,object$ctstanmodel$manifestNames)
  out$resiCovStdNote <- 'Standardised covariance of residuals'
  }
  
  parnames <- object$setup$matsetup$parname[object$setup$matsetup$when==0 & object$setup$matsetup$param > 0]
  parindices <- object$setup$matsetup$param[object$setup$matsetup$when==0 & object$setup$matsetup$param > 0]
  pars <- cbind(parnames,parindices)
  pars<-pars[!duplicated(pars[,1,drop=FALSE]),,drop=FALSE]
  parnames <- pars[as.numeric(pars[,2,drop=FALSE]) >0, 1]
  # parnames <- unique(parnames)
  parnamesiv <- parnames[object$data$indvaryingindex]
  
  #### generate covcor matrices of raw and transformed subject level params
  
  iter=dim(e$rawpopcov)[1]
  if(!is.null(iter)){ #then there is some individual variation so continue
    nindvarying=dim(e$rawpopcov)[2]
    
    if(nindvarying>1){
      
      #raw pop distribution params
      # browser()
      dimrawpopcorr <- dim(e$rawpopcorr)
      # if(!'stanfit' %in% class(object$stanfit)) 
        rawpopcorr= array(e$rawpopcorr,dim=c(dimrawpopcorr[1],1,dimrawpopcorr[2] * dimrawpopcorr[3]))
      # if('stanfit' %in% class(object$stanfit)) rawpopcorr= rstan::extract(object$stanfit,pars='rawpopcorr',permuted=FALSE)

      rawpopcorrout <- suppressWarnings(monitor(rawpopcorr, digits_summary=digits,warmup=0,print = FALSE)[lower.tri(diag(nindvarying)),c(monvars,'n_eff','Rhat'),drop=FALSE])
      if(!'stanfit' %in% class(object$stanfit)) rawpopcorrout <- rawpopcorrout[,-which(colnames(rawpopcorrout) %in% c('n_eff','Rhat')),drop=FALSE]

      rownames(rawpopcorrout) <- matrix(paste0('',parnamesiv,'__',rep(parnamesiv,each=length(parnamesiv))),
        length(parnamesiv),length(parnamesiv))[lower.tri(diag(nindvarying)),drop=FALSE]
      
      rawpopcorrout <- cbind(rawpopcorrout,rawpopcorrout[,'mean'] / rawpopcorrout[,'sd'])
      colnames(rawpopcorrout)[ncol(rawpopcorrout)] <- 'z'
      
      out$rawpopcorr = round(rawpopcorrout,digits)
    }
  }
  
  if(priorcheck & object$standata$nopriors==0) {
    priorcheckres <- priorcheckreport(object,...)
    if(nrow(priorcheckres$priorcheck) > 0) out = c(out,priorcheckres)
  }
  
  if(object$ctstanmodel$n.TIpred > 0) {
    
    if('stanfit' %in% class(object$stanfit)){
      rawtieffect <- rstan::extract(object$stanfit,permuted=FALSE,pars='TIPREDEFFECT')
      tidiags <- suppressWarnings(monitor(rawtieffect,warmup=0,digits_summary = digits,print = FALSE))
    }
    
    tieffect <- array(e$linearTIPREDEFFECT,dim=c(dim(e$linearTIPREDEFFECT)[1], 1, length(parnames) * dim(e$linearTIPREDEFFECT)[3]))
    tieffectnames <- paste0('tip_',rep(object$ctstanmodel$TIpredNames,each=length(parnames)),'_',parnames)
    dimnames(tieffect)<-list(c(),c(),tieffectnames)
    tipreds = suppressWarnings(monitor(tieffect,warmup = 0,print = FALSE)[,monvars,drop=FALSE])
    if('stanfit' %in% class(object$stanfit)) tipreds <- cbind(tipreds,tidiags[,c('n_eff','Rhat'),drop=FALSE])
    tipreds <- tipreds[c(object$data$TIPREDEFFECTsetup)>0,,drop=FALSE]
    z = tipreds[,'mean'] / tipreds[,'sd'] 
    out$tipreds= round(cbind(tipreds,z),digits) #[order(abs(z)),]
  }
  

  
  if(parmatrices){
    
    # #check if stanfit object can be used
    # sf <- object$stanfit
    # npars <- try(get_num_upars(sf),silent=TRUE) #$stanmodel)
    # 
    # if(class(npars)=='try-error'){ #in case R has been restarted or similar
    #   standataout <- object$standata
    #   smf <- stan_reinitsf(object$stanmodel,standataout)
    # }
    object$standata$savescores <- 0L
    object$standata$gendata <- 0L
    object$standata$dokalman <- 0L
    object$standata$popcovn <- 5L
    sf <- stan_reinitsf(object$stanmodel,data=object$standata)
    parmatlists <- try(apply(ctStanRawSamples(object),
      # sample(x = 1:dim(e$rawpopmeans)[1],
      #   size = dim(e$rawpopmeans)[1],  #min(parmatsamples,
      #   replace = FALSE),],
      1,ctStanParMatrices,fit=object,timeinterval=timeinterval,sf=sf))
    
    if(!'try-error' %in% class(parmatlists)[1]){
      parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
      parmats <- matrix(NA,nrow=length(unlist(parmatlists[[1]])),ncol=7)
      rownames(parmats) <- paste0('r',1:nrow(parmats))
      counter=0
      for(mati in 1:length(parmatlists[[1]])){
        if(all(dim(parmatlists[[1]][[mati]]) > 0)){
          for(coli in 1:ncol(parmatlists[[1]][[mati]])){
            for(rowi in 1:nrow(parmatlists[[1]][[mati]])){
              counter=counter+1
              new <- matrix(c(
                rowi,
                coli,
                mean(parmatarray[counter,],na.rm=TRUE),
                sd(parmatarray[counter,],na.rm=TRUE),
                quantile(parmatarray[counter,],probs=c(.025,.5,.975),na.rm=TRUE)),
                nrow=1)
              try(rownames(parmats)[counter] <- names(parmatlists[[1]])[mati])
              try(parmats[counter,]<-new)
            }}}}
      
      colnames(parmats) <- c('Row','Col', 'Mean','Sd','2.5%','50%','97.5%')
      
      #remove certain parmatrices lines
      removeindices <- which(rownames(parmats) == 'MANIFESTVAR' & parmats[,'Row'] != parmats[,'Col'])
      
      removeindices <- c(removeindices,which((rownames(parmats) %in% c('MANIFESTVAR','T0VAR','DIFFUSION','dtDIFFUSION','asymDIFFUSION',
        'T0VARcor','DIFFUSIONcor','DIFFUSIONcov','dtDIFFUSIONcor','asymDIFFUSIONcor') &  parmats[,'Row'] < parmats[,'Col'])))
      
      removeindices <- c(removeindices,which((rownames(parmats) %in% c('T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') & 
          parmats[,'Row'] == parmats[,'Col'])))
      
      parmats <- parmats[-removeindices,]
      
      out$parmatrices=round(parmats,digits=digits)
      
      out$parmatNote=paste0('Population mean parameter matrices calculated with time interval of ', timeinterval,' for discrete time (dt) matrices. ',
        'Covariance related matrices shown as covariance matrices, correlations have (cor) suffix. Asymptotic (asym) matrices based on infinitely large time interval.')
    }
    if('try-error' %in% class(parmatlists)[1]) out$parmatNote = 'Could not calculate parameter matrices'
  }
    
  
  
  if('stanfit' %in% class(object$stanfit)){
    # browser()
    popsd=smr$summary[c(grep('^popsd',rownames(smr$summary),fixed=FALSE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE] #[ object$data$indvaryingindex,,drop=FALSE]
    rownames(popsd)=parnames[ object$data$indvaryingindex]
    # popmeans=smr$summary[c(grep('hmean_',rownames(smr$summary))),
    #   c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=smr$summary[c(grep('popmeans[', rownames(smr$summary),fixed=TRUE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=popmeans[(nrow(popmeans)/2+1):nrow(popmeans),,drop=FALSE]
    rownames(popmeans) <- parnames
    

    logposterior=smr$summary[c(grep('lp',rownames(smr$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
  }
  
  if(!'stanfit' %in% class(object$stanfit)){ #if optimized / importance sampled
    
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
  
  if(!'stanfit' %in% class(object$stanfit)) {
    
    out$loglik=loglik
    out$npars = npars
    out$aic = aic
  }
  out$logposterior=logposterior
  if(!'stanfit' %in% class(object$stanfit)) out$nsamples <- nrow(object$stanfit$samples)
  
  if(!parmatrices) out$parmatNote <- 'For additional summary matrices, use argument: parmatrices = TRUE'
  
  
  
  # out$posteriorpredictive=round(smr$summary[c(grep('stateppll',rownames(smr$summary))),
  #     c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)
  # }
  
  
  # if(!'stanfit' %in% class(object$stanfit)){ #optimization summary
  #   out=list()
  #   out$popmeans=object$stanfit$transformedpars[grep('hmean_',rownames(object$stanfit$transformedpars)),]
  #   out$popsd=object$stanfit$transformedpars[grep('hsd_',rownames(object$stanfit$transformedpars)),]
  #   out$logprob=round(-object$stanfit$optimfit$value)
  # }
  
  return(out)
}
