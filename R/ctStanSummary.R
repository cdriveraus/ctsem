#' summary.ctStanFit
#'
#' Summarise a ctStanFit object that was fit using \code{\link{ctStanFit}}. 
#' 
#' @param object fit object from \code{\link{ctStanFit}}, of class ctStanFit.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter calculations
#' reported in summary. 
#' @param digits integer denoting number of digits to report.
#' @param ... Unused at present.
#' @return List containing summary items.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=3,...){
  if(class(object) != 'ctStanFit') stop('Not a ctStanFit object!')
  
  if(class(object$stanfit)=='stanfit'){ #summary of samples
    
    
    s<-getMethod('summary','stanfit')(object$stanfit)
    
    if('98%' %in% colnames(s$summary)) colnames(s$summary)[colnames(s$summary)=='98%'] <- '97.5%'
    
    parnames=gsub('hsd_','',object$stanfit@model_pars[grep('hsd',object$stanfit@model_pars)])
    
    #### generate covcor matrices of raw and transformed subject level params
    out=list()
    e=extract(object$stanfit)
    
    iter=dim(e$hypercorrsqrt)[1]
    if(!is.null(iter)){ #then there is some individual variation so continue
      npars=dim(e$hypercorrsqrt)[2]
      
      if(npars>1){
        
        getMean=function(myarray){
          out=matrix(NA,nrow=npars,ncol=npars)
          for(i in 1:nrow(out)){
            for(j in 1:ncol(out)){
              out[i,j]<-round(mean(myarray[i,j,]),digits=digits)
            }}
          return(out)
        }
        
        getSd=function(myarray){
          out=matrix(NA,nrow=npars,ncol=npars)
          for(i in 1:nrow(out)){
            for(j in 1:ncol(out)){
              out[i,j]<-round(sd(myarray[i,j,]),digits=digits)
            }}
          return(out)
        }
        
        
        #transformed subject level params
        hypercorr_transformed= array(sapply(1:iter, function(x) cor(e$indparams[x,,])),dim=c(npars,npars,iter))
        hypercov_transformed= array(sapply(1:iter, function(x) cov(e$indparams[x,,])),dim=c(npars,npars,iter))
        
        hypercorr_transformedmean=getMean(hypercorr_transformed)
        hypercorr_transformedsd=getSd(hypercorr_transformed)
        
        hypercov_transformedmean=getMean(hypercov_transformed)
        hypercov_transformedsd=getSd(hypercov_transformed)
        
        hypercovcor_transformedmean=hypercov_transformedmean
        hypercovcor_transformedmean[lower.tri(diag(npars))]=hypercorr_transformedmean[lower.tri(diag(npars))]
        
        hypercovcor_transformedsd=hypercov_transformedsd
        hypercovcor_transformedsd[lower.tri(diag(npars))]=hypercorr_transformedsd[lower.tri(diag(npars))]
        
        dimnames(hypercovcor_transformedsd)<-list(parnames,parnames)
        dimnames(hypercovcor_transformedmean)<-list(parnames,parnames)
        
        #raw subject level params
        hypercorr=array(unlist(lapply(1:iter,function(x){ #get array of hypercorr samples
          hypercorrsqrt=matrix(e$hypercorrsqrt[x,,],nrow=npars)
          hypercorrsqrt%*% t(hypercorrsqrt)
        })),dim=c(npars,npars,iter))
        
        
        popcorr <- ctCollapse(hypercorr_transformed,3,mean)[lower.tri(diag(dim(hypercorr_transformed)[1]))]
        popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,sd)[lower.tri(diag(dim(hypercorr_transformed)[1]))])
        popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.025))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
        popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.5))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
        popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.975))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
        colnames(popcorr) <- c('mean','sd','2.5%','50%','97.5%')
        rownames(popcorr) <- matrix(paste0('corr_',parnames,'__',rep(parnames,each=length(parnames))),
          length(parnames),length(parnames))[lower.tri(diag(dim(hypercorr)[1]))]
        popcorr <- round(popcorr,digits=digits)
        
        popcorr <- cbind(popcorr,popcorr[,'mean'] / popcorr[,'sd'])
        colnames(popcorr)[ncol(popcorr)] <- 'z'
        
        popcorr <- popcorr[order(abs(popcorr[,'z'])),,drop=FALSE]
        
        hypercorrmean= ctCollapse(array(apply(e$hypercorrsqrt,1,function(x) x%*% t(x)),dim = dim(e$hypercorrsqrt)[c(2,3,1)]),3,mean)
        hypercorrsd= ctCollapse(array(apply(e$hypercorrsqrt,1,function(x) x%*% t(x)),dim = dim(e$hypercorrsqrt)[c(2,3,1)]),3,sd)
        
        dimnames(hypercorrmean)<-list(parnames,parnames)
        dimnames(hypercorrsd)<-list(parnames,parnames)
        
        out=list(note1='The following matrix is the posterior means of the raw parameter population distribution correlation matrix',
          hypercorr_mean=hypercorrmean,
          note2='The following matrix is the posterior std dev. of the raw parameter population distribution correlation matrix',
          hypercorr_sd=hypercorrsd,
          note3=paste0('The following matrix is the posterior mean of the correlation and covariance matrix of subject level parameters,', 
            ' with correlations on the lower triangle'),
          hypercovcor_transformedmean=hypercovcor_transformedmean,
          note4=paste('The following matrix is the posterior std dev. of the correlation and covariance matrix of subject level parameters,', 
            'with correlations on the lower triangle'),
          hypercovcor_transformedsd=hypercovcor_transformedsd)
        
        out$popcorr = popcorr
      }
    }
    
    
    
    if(object$ctstanmodel$n.TIpred > 0) {
      out$tipreds=round(s$summary[c(grep('tipred_',rownames(s$summary))),
        c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
      z = out$tipreds[,'mean'] / out$tipreds[,'sd'] 
      out$tipreds= cbind(out$tipreds,z)[order(abs(z)),]
    }
    
  
    parmatlists <- apply(e$hypermeans,1,ctStanParMatrices,model=object,timeinterval=timeinterval)
    parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
    parmats <- matrix(0,nrow=0,ncol=7)
    counter=0
    for(mati in 1:length(parmatlists[[1]])){
      for(coli in 1:nrow(parmatlists[[1]][[mati]])){
        for(rowi in 1:ncol(parmatlists[[1]][[mati]])){
          counter=counter+1
          new <- matrix(c(
            rowi,
            coli,
            mean(parmatarray[counter,]),
            sd(parmatarray[counter,]),
            quantile(parmatarray[counter,],probs=c(.025,.5,.975))),
            nrow=1)
          rownames(new) = names(parmatlists[[1]])[mati]
          parmats<-rbind(parmats, new)
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
    
    
    out$popsd=round(s$summary[c(grep('hsd_',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    out$popmeans=round(s$summary[c(grep('hmean_',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    out$popNote=paste0('popmeans and popsd are reported as specified in ctModel -- diagonals of covariance related matrices are std. deviations, ',
      'off-diagonals are unconstrained correlation square roots.')
    
    out$logprob=round(s$summary[c(grep('lp',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],digits=digits)
    
    
    
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
