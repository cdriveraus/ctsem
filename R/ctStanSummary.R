#' ctStanSummary
#'
#' Summarise a Stan fit object fit using \link{\code{ctStanFit}}. 
#' 
#' @param ctstanfitobject Stan fit object from \link{\code{ctStanFit}}. 
#' @import rstan
#' @export

ctStanSummary<-function(ctstanfitobject){
  
  trace=FALSE
  density=FALSE
  
  fit<-ctstanfitobject

if(trace==TRUE) {

  getMethod('traceplot','stanfit')(fit,fit@model_pars[grep('output_hmean',fit@model_pars,fixed=TRUE)],inc_warmup=F)
  getMethod('traceplot','stanfit')(fit,fit@model_pars[grep('output_hsd',fit@model_pars,fixed=TRUE)],inc_warmup=F)
}
if( density == TRUE) {
  densfunc<-rstan:::stan_dens
  densfunc(fit,c('lp__',fit@model_pars[grep('output_hmean',fit@model_pars,fixed=TRUE)]),inc_warmup=F)
  getMethod('stan_dens','stanfit')(fit,c('lp__',fit@model_pars[grep('output_hsd',fit@model_pars,fixed=TRUE)]),inc_warmup=F)
}
  
s<-getMethod('summary','stanfit')(fit)

parnames=gsub('output_hsd_','',fit@model_pars[grep('hsd',fit@model_pars)])

#### generate covcor matrices of raw and transformed subject level params
out=list()
e=extract(fit)
# browser()
iter=dim(e$hypercorrchol)[1]
if(!is.null(iter)){ #then there is some individual variation so continue
npars=dim(e$hypercorrchol)[2]

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

# browser()


hypercorrmean=getMean(hypercorr)
hypercorrsd=getSd(hypercorr)

hypercorrcholmeanpars=s$summary[grep('hypercorrchol',rownames(s$summary)),'mean']
hypercorrcholmean=matrix(hypercorrcholmeanpars,byrow=T,nrow=sqrt(length(hypercorrcholmeanpars)))

hypercorrcholsdpars=s$summary[grep('hypercorrchol',rownames(s$summary)),'sd']
hypercorrcholsd=matrix(hypercorrcholsdpars,byrow=T,nrow=sqrt(length(hypercorrcholsdpars)))

dimnames(hypercorrcholmean)<-list(parnames,parnames)
dimnames(hypercorrcholsd)<-list(parnames,parnames)

out=list(note1='The following matrix is the posterior means for the raw subject level parameters correlation matrix',
  hypercorr_means=hypercorrcholmean %*% t(hypercorrcholmean),
  note2='The following matrix is the posterior std dev. for the raw subject level parameters correlation matrix',
  hypercorr_sd=hypercorrcholsd %*% t(hypercorrcholsd),
  note3=paste0('The following matrix is the posterior mean for the post transformation subject level parameters correlation and covariance matrix,', 
    ' with correlations on the lower triangle'),
  hypercovcor_transformedmean=hypercovcor_transformedmean,
  note4=paste('The following matrix is the posterior std dev. for the post transformation subject level parameters correlation and covariance matrix,', 
    'with correlations on the lower triangle'),
  hypercovcor_transformedsd=hypercovcor_transformedsd)

}

out$outputpars=s$summary[c(grep('output',rownames(s$summary)),grep('lp',rownames(s$summary))),
  c('mean','sd','n_eff','Rhat')]

return(out)
}
