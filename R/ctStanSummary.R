#' summary.ctStanFit
#'
#' Summarise a ctStanFit object that was fit using \code{\link{ctStanFit}}. 
#' 
#' @param object fit object from \code{\link{ctStanFit}}, of class ctStanFit.
#' @param ... Unused at present.
#' @return List containing summary items.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,...){

  if(class(object$stanfit)=='stanfit'){ #summary of samples
  
  
s<-getMethod('summary','stanfit')(object$stanfit)

if('98%' %in% colnames(s$summary)) colnames(s$summary)[colnames(s$summary)=='98%'] <- '97.5%'

parnames=gsub('hsd_','',object$stanfit@model_pars[grep('hsd',object$stanfit@model_pars)])

#### generate covcor matrices of raw and transformed subject level params
out=list()
e=extract(object$stanfit)

iter=dim(e$hypercorrchol)[1]
if(!is.null(iter)){ #then there is some individual variation so continue
npars=dim(e$hypercorrchol)[2]

if(npars>1){

getMean=function(myarray){
  out=matrix(NA,nrow=npars,ncol=npars)
  for(i in 1:nrow(out)){
    for(j in 1:ncol(out)){
      out[i,j]<-round(mean(myarray[i,j,]),3)
    }}
  return(out)
  }

getSd=function(myarray){
  out=matrix(NA,nrow=npars,ncol=npars)
  for(i in 1:nrow(out)){
    for(j in 1:ncol(out)){
      out[i,j]<-round(sd(myarray[i,j,]),3)
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
  hypercorrchol=matrix(e$hypercorrchol[x,,],nrow=npars)
  hypercorrchol%*% t(hypercorrchol)
})),dim=c(npars,npars,iter))


popcorr <- ctCollapse(hypercorr_transformed,3,mean)[lower.tri(diag(dim(hypercorr_transformed)[1]))]
popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,sd)[lower.tri(diag(dim(hypercorr_transformed)[1]))])
popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.025))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.5))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
popcorr <- cbind(popcorr,ctCollapse(hypercorr_transformed,3,quantile,probs=c(.975))[lower.tri(diag(dim(hypercorr_transformed)[1]))])
colnames(popcorr) <- c('mean','sd','2.5%','50%','97.5%')
rownames(popcorr) <- matrix(paste0('corr_',parnames,'__',rep(parnames,each=length(parnames))),
  length(parnames),length(parnames))[lower.tri(diag(dim(hypercorr)[1]))]
popcorr <- round(popcorr,3)

popcorr <- cbind(popcorr,popcorr[,'mean'] / popcorr[,'sd'])
colnames(popcorr)[ncol(popcorr)] <- 'z'

popcorr <- popcorr[order(popcorr[,'z']),,drop=FALSE]

hypercorrmean= ctCollapse(array(apply(e$hypercorrchol,1,function(x) x%*% t(x)),dim = dim(e$hypercorrchol)[c(2,3,1)]),3,mean)
hypercorrsd= ctCollapse(array(apply(e$hypercorrchol,1,function(x) x%*% t(x)),dim = dim(e$hypercorrchol)[c(2,3,1)]),3,sd)

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
  c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)
  z = out$tipreds[,'mean'] / out$tipreds[,'sd'] 
  out$tipreds= cbind(out$tipreds,z)
}

out$popsd=round(s$summary[c(grep('hsd_',rownames(s$summary))),
  c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)

out$popmeans=round(s$summary[c(grep('hmean_',rownames(s$summary))),
  c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)

out$note1=paste0('Parameters are reported as specified in ctModel -- diagonals of covariance related matrices are std. deviations, ',
'off-diagonals are sqrt(logit(x/2+.5)) of partial correlations. Full covariance matrices can be attained using ctStanContinuousPars')

out$logprob=round(s$summary[c(grep('lp',rownames(s$summary))),
  c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)

  
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
