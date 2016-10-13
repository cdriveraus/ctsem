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
return(s$summary[c(grep('output',rownames(s$summary)),grep('lp',rownames(s$summary))),
  c('mean','sd','n_eff','Rhat')])
}
