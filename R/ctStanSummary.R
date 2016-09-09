#' ctStanSummary
#'
#' Summarise a Stan fit object fit using \link{\code{ctStanFit}}. 
#' 
#' @param ctstanfitobject Stan fit object from \link{\code{ctStanFit}}. 
#' @import rstan
#' @export

ctStanSummary<-function(ctstanfitobject,trace=TRUE,density=TRUE){
  
  fit<-ctstanfitobject

if(trace==TRUE) {
  rstan:::traceplot(fit,fit@model_pars[grep('output_hmean',fit@model_pars,fixed=TRUE)],inc_warmup=F)
  rstan:::traceplot(fit,fit@model_pars[grep('output_hsd',fit@model_pars,fixed=TRUE)],inc_warmup=F)
}
if( density == TRUE) {
  rstan:::stan_dens(fit,c('lp__',fit@model_pars[grep('output_hmean',fit@model_pars,fixed=TRUE)]),inc_warmup=F)
  rstan:::stan_dens(fit,c('lp__',fit@model_pars[grep('output_hsd',fit@model_pars,fixed=TRUE)]),inc_warmup=F)
}
s<-getMethod('summary','stanfit')(fit)
return(s$summary[c(grep('output',rownames(s$summary)),grep('lp',rownames(s$summary))),
  c('mean','sd','n_eff','Rhat')])
}
