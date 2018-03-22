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
#' @param ... Unused at present.
#' @return List containing summary items.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=3,parmatrices=FALSE,...){
  if(class(object) != 'ctStanFit') stop('Not a ctStanFit object!')
  
  if(class(object$stanfit)=='stanfit'){ #summary of samples

    
    s<-getMethod('summary','stanfit')(object$stanfit)
    
    if('98%' %in% colnames(s$summary)) colnames(s$summary)[colnames(s$summary)=='98%'] <- '97.5%'
    
    parnames=gsub('hsd_','',object$stanfit@model_pars[grep('hsd',object$stanfit@model_pars)])
    
    #### generate covcor matrices of raw and transformed subject level params
    out=list()
    e=extract(object$stanfit)
    
    iter=dim(e$rawpopcorrsqrt)[1]
    if(!is.null(iter)){ #then there is some individual variation so continue
      nindvarying=dim(e$rawpopcorrsqrt)[2]
      
      if(nindvarying>1){
        
        getMean=function(myarray){
          out=matrix(NA,nrow=nindvarying,ncol=nindvarying)
          for(i in 1:nrow(out)){
            for(j in 1:ncol(out)){
              out[i,j]<-round(mean(myarray[i,j,]),digits=digits)
            }}
          return(out)
        }
        
        getSd=function(myarray){
          out=matrix(NA,nrow=nindvarying,ncol=nindvarying)
          for(i in 1:nrow(out)){
            for(j in 1:ncol(out)){
              out[i,j]<-round(sd(myarray[i,j,]),digits=digits)
            }}
          return(out)
        }
        

        #transformed subject level params
        rawpopcorr_transformed= array(sapply(1:iter, function(x) cor(e$indparams[x,,])),dim=c(nindvarying,nindvarying,iter))
        rawpopcov_transformed= array(sapply(1:iter, function(x) cov(e$indparams[x,,])),dim=c(nindvarying,nindvarying,iter))
        
        rawpopcorr_transformedmean=getMean(rawpopcorr_transformed)
        rawpopcorr_transformedsd=getSd(rawpopcorr_transformed)
        
        rawpopcov_transformedmean=getMean(rawpopcov_transformed)
        rawpopcov_transformedsd=getSd(rawpopcov_transformed)
        
        rawpopcovcor_transformedmean=rawpopcov_transformedmean
        rawpopcovcor_transformedmean[lower.tri(diag(nindvarying))]=rawpopcorr_transformedmean[lower.tri(diag(nindvarying))]
        
        rawpopcovcor_transformedsd=rawpopcov_transformedsd
        rawpopcovcor_transformedsd[lower.tri(diag(nindvarying))]=rawpopcorr_transformedsd[lower.tri(diag(nindvarying))]
        
        dimnames(rawpopcovcor_transformedsd)<-list(parnames,parnames)
        dimnames(rawpopcovcor_transformedmean)<-list(parnames,parnames)
        
        #raw subject level params
        rawpopcorr= e$rawpopcorr 
    
        rawpopcorrout <- ctCollapse(rawpopcorr,1,mean)[lower.tri(diag(nindvarying))]
        rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,sd)[lower.tri(diag(nindvarying))])
        rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.025))[lower.tri(diag(nindvarying))])
        rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.5))[lower.tri(diag(nindvarying))])
        rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.975))[lower.tri(diag(nindvarying))])
        colnames(rawpopcorrout) <- c('mean','sd','2.5%','50%','97.5%')
        rownames(rawpopcorrout) <- matrix(paste0('corr_',parnames,'__',rep(parnames,each=length(parnames))),
          length(parnames),length(parnames))[lower.tri(diag(nindvarying))]
        rawpopcorrout <- round(rawpopcorrout,digits=digits)
        
        rawpopcorrout <- cbind(rawpopcorrout,rawpopcorrout[,'mean'] / rawpopcorrout[,'sd'])
        colnames(rawpopcorrout)[ncol(rawpopcorrout)] <- 'z'
        
        rawpopcorrout <- rawpopcorrout[order(abs(rawpopcorrout[,'z'])),,drop=FALSE]
        
        rawpopcorrmean= ctCollapse(e$rawpopcorr,1,mean)
        rawpopcorrsd= ctCollapse(e$rawpopcorr,1,sd)
        
        dimnames(rawpopcorrmean)<-list(parnames,parnames)
        dimnames(rawpopcorrsd)<-list(parnames,parnames)
        
        out=list(note1='The following matrix contains the posterior means of the raw parameter population distribution correlation matrix',
          rawpopcorr_mean=rawpopcorrmean,
          note2='The following matrix contains the posterior std dev. of the raw parameter population distribution correlation matrix',
          rawpopcorr_sd=rawpopcorrsd,
          note3=paste0('The following matrix is the posterior mean of the correlation and covariance matrix of subject level parameters,', 
            ' with correlations on the lower triangle'),
          popcovcor_mean=rawpopcovcor_transformedmean,
          note4=paste('The following matrix is the posterior std dev. of the correlation and covariance matrix of subject level parameters,', 
            'with correlations on the lower triangle'),
          popcovcor_sd=rawpopcovcor_transformedsd)
        
        out$rawpopcorr = rawpopcorrout
      }
    }
    
    
    
    if(object$ctstanmodel$n.TIpred > 0) {
      out$tipreds=round(s$summary[c(grep('tipred_',rownames(s$summary))),
        c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
      z = out$tipreds[,'mean'] / out$tipreds[,'sd'] 
      out$tipreds= cbind(out$tipreds,z)[order(abs(z)),]
    }
    
  if(parmatrices){
    parmatlists <- try(apply(e$rawpopmeans,1,ctStanParMatrices,model=object,timeinterval=timeinterval))
    if(class(parmatlists)!='try-error'){
    parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
    parmats <- matrix(NA,nrow=length(unlist(parmatlists[[1]])),ncol=7)
    rownames(parmats) <- paste0('r',1:nrow(parmats))
    counter=0
    for(mati in 1:length(parmatlists[[1]])){
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
        }}}
    colnames(parmats) <- c('Row','Col', 'Mean','Sd','2.5%','50%','97.5%')
    
    #remove certain parmatrices lines
    removeindices <- which(rownames(parmats) == 'MANIFESTVAR' & parmats[,'Row'] != parmats[,'Col'])
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('MANIFESTVAR','T0VAR','DIFFUSION','dtDIFFUSION','asymDIFFUSION',
      'T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') &  parmats[,'Row'] < parmats[,'Col'])))
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') & 
        parmats[,'Row'] == parmats[,'Col'])))
    
    parmats <- parmats[-removeindices,]
    
    
    out$parmatrices=round(parmats,digits=digits)
    
    out$parmatNote=paste0('Population mean parameter matrices calculated with time interval of ', timeinterval,' for discrete time (dt) matrices. ',
      'Covariance related matrices shown as covariance matrices, correlations have (cor) suffix. Asymptotic (asym) matrices based on infinitely large time interval.')
    }
    if(class(parmatlists)=='try-error') out$parmatNote = 'Could not calculate parameter matrices'
  }
    
    out$popsd=round(s$summary[c(grep('hsd_',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    out$popmeans=round(s$summary[c(grep('hmean_',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    out$popNote=paste0('popmeans and popsd are reported as specified in ctModel -- diagonals of covariance related matrices are std. deviations, ',
      'off-diagonals are unconstrained correlation square roots.')
    
    out$logprob=round(s$summary[c(grep('lp',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    if(!parmatrices) out$parmatNote <- 'For additional summary matrices, use argument: parmatrices = TRUE'
    
    
    
    # out$posteriorpredictive=round(s$summary[c(grep('stateppll',rownames(s$summary))),
    #     c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)
  }
  
  
  if(class(object$stanfit)!='stanfit'){ #optimization summary
    out=list()
    out$popmeans=object$stanfit$transformedpars[grep('hmean_',rownames(object$stanfit$transformedpars)),]
    out$popsd=object$stanfit$transformedpars[grep('hsd_',rownames(object$stanfit$transformedpars)),]
    out$logprob=round(-object$stanfit$optimfit$value)
  }

  return(out)
}
