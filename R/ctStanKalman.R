#' Get Kalman filter estimates from a ctStanFit object
#'
#' @param fit fit object from \code{\link{ctStanFit}}.
#' @param nsamples either NA (to extract all) or a positive integer from 1 to maximum samples in the fit.
#' @param cores Integer number of cpu cores to use. Only needed if savescores was set to FALSE when fitting.
#' @param collapsefunc function to apply over samples, such as \code{mean}
#' @param ... additional arguments to collpsefunc.
#'
#' @return list containing Kalman filter elements, each element in array of
#' iterations, data row, variables. llrow is the log likelihood for each row of data.
#' @export
#'
#' @examples 
#' \donttest{
#' k=ctStanKalman(ctstantestfit)
#' }
ctStanKalman <- function(fit,nsamples=NA,collapsefunc=NA,cores=2,...){
  if(class(fit)!='ctStanFit') stop('Not a ctStanFit object')
 if(class(collapsefunc) %in% 'function' ) e=extract(fit)

  if(!class(collapsefunc) %in% 'function' || length(dim(e$k))==0){
    message('State estimates not saved, computing...')
    standata <- fit$standata
    standata$savescores <- 1L
    # smf <- stan_reinitsf(fit$stanmodel, standata)
    if(!class(fit$stanfit) %in% 'stanfit') {
      samples = fit$stanfit$rawposterior
    } else {
      samples = t(stan_unconstrainsamples(fit$stanfit,fit$standata))
    }
    if(!is.na(nsamples)) samples <- samples[sample(1:nrow(samples),nsamples),,drop=FALSE]
    if(class(collapsefunc) %in% 'function') samples = matrix(apply(samples,2,collapsefunc,...),ncol=ncol(samples))
    e=stan_constrainsamples(sm = fit$stanmodel,standata = standata,samples = samples,cores=cores)
  }
  
  k=e$kalman
  
  k[k==99999] <- NA #for missingness
  nlatent <- fit$standata$nlatent
  nmanifest <- fit$standata$nmanifest
  dimnames(k) = list(iter=1:dim(k)[1],drow=1:dim(k)[2],
    kalman=paste0(c(rep('lln',nmanifest),
      rep('llscale',nmanifest),rep('err',nmanifest),rep('yprior',nmanifest),rep('etaprior',nlatent),rep('etaupd',nlatent)),
      c(1:nmanifest,1:nmanifest,1:nmanifest,1:nmanifest,1:nlatent,1:nlatent)))
  
  
  
  obj <- c('lln','llscale','err','yprior','etaprior','etaupd')
  
  
  lln=k[,,1:nmanifest,drop=FALSE]
  llscale=k[,,(nmanifest*1+1):(nmanifest*1+nmanifest),drop=FALSE]
  err=k[,,(nmanifest*2+1):(nmanifest*2+nmanifest),drop=FALSE]
  yprior=k[,,(nmanifest*3+1):(nmanifest*3+nmanifest),drop=FALSE]
  etaprior=k[,,(nmanifest*4+1):(nmanifest*4+nlatent),drop=FALSE]
  etaupd=k[,,(nmanifest*4+nlatent+1):(nmanifest*4+nlatent*2),drop=FALSE]
  
  llvec = apply(lln,1:2,function(x) {
    sum(dnorm(x[!is.na(x)],log = TRUE))
  })
  llrow = llvec - apply(llscale, 1:2, function(x) sum(x,na.rm=TRUE))
  
  #covariance
  etapriorcov <- e$etapriorcov
  etaupdcov <- e$etaupdcov
  etapriorcov[etapriorcov==99999] <- NA
  etaupdcov[etaupdcov==99999] <- NA
  ypriorcov <- e$ypriorcov
  ypriorcov[ypriorcov==99999] <- NA
  etasmoothcov <- e$etasmoothcov[,,1:nlatent,1:nlatent,drop=FALSE]
  etasmoothcov[etasmoothcov==99999] <- NA
  etasmooth <- e$etasmooth[,,1:nlatent,drop=FALSE]
  etasmooth[etasmooth==99999] <- NA
  
    return(list(time=cbind(fit$standata$time), lln=lln,llscale=llscale,err=err,
      y=matrix(fit$data$Y,ncol=ncol(fit$data$Y)), yprior=yprior,ypriorcov=ypriorcov,
      etaprior=etaprior, etapriorcov=etapriorcov,
      etaupd=etaupd,etaupdcov=etaupdcov,
      etasmooth=etasmooth,etasmoothcov=etasmoothcov,
      llrow=llrow))
}



