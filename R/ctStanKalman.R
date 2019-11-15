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
  # if(class(collapsefunc) %in% 'function' ) e=extract(fit)
  
  # if(!class(collapsefunc) %in% 'function' || length(dim(e$k))==0){
  message('Computing state estimates..')
  standata <- fit$standata
  standata$savescores <- 1L
  # smf <- stan_reinitsf(fit$stanmodel, standata)
  samples<-ctStanRawSamples(fit)
  if(!is.na(nsamples)) samples <- samples[sample(1:nrow(samples),nsamples),,drop=FALSE] else nsamples <- nrow(samples)
  if(is.function(collapsefunc)) samples = matrix(apply(samples,2,collapsefunc,...),ncol=ncol(samples))
  e=stan_constrainsamples(sm = fit$stanmodel,standata = standata,samples = samples,cores=cores)
  # }
  
  k=e$kalman
  
  k[k==99999] <- NA #for missingness
  nlatent <- fit$standata$nlatent
  nmanifest <- fit$standata$nmanifest
  dimnames(k) = list(iter=1:dim(k)[1],drow=1:dim(k)[2],
    kalman=paste0(c(rep('lln',nmanifest),
      rep('llscale',nmanifest),rep('stderr',nmanifest),rep('yprior',nmanifest),rep('etaprior',nlatent),rep('etaupd',nlatent)),
      c(1:nmanifest,1:nmanifest,1:nmanifest,1:nmanifest,1:nlatent,1:nlatent)))
  
  
  lln=k[,,1:nmanifest,drop=FALSE]
  llscale=k[,,(nmanifest*1+1):(nmanifest*1+nmanifest),drop=FALSE]
  stderr=k[,,(nmanifest*2+1):(nmanifest*2+nmanifest),drop=FALSE]
  e$yprior=k[,,(nmanifest*3+1):(nmanifest*3+nmanifest),drop=FALSE]
  # etaprior=k[,,(nmanifest*4+1):(nmanifest*4+nlatent),drop=FALSE]
  # etaupd=k[,,(nmanifest*4+nlatent+1):(nmanifest*4+nlatent*2),drop=FALSE]
  
  llvec = apply(lln,1:2,function(x) {
    sum(dnorm(x[!is.na(x)],log = TRUE))
  })
  llrow = llvec - apply(llscale, 1:2, function(x) sum(x,na.rm=TRUE))
  
  
  
  y=matrix(fit$standata$Y,ncol=ncol(fit$standata$Y),dimnames = list(NULL,fit$ctstanmodel$manifestNames))
  y[y==99999] <- NA
  
  
  out=list(time=cbind(fit$standata$time), lln=lln,llscale=llscale,stderr=stderr,
    y=y, 
    llrow=llrow)
  
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
            dimnames(out[[ref]]) <- list(NULL, NULL, fit$ctstanmodel$latentNames[1:nlatent])
          } else { #for cov
            out[[ref]] <- out[[ref]][,,1:nlatent,1:nlatent,drop=FALSE]
          }
        }
      }
    }
  }
  
  for(typei in c('prior','upd','smooth')){
    out[[paste0('err',typei)]] <- aaply(out[[paste0('y',typei)]],1, function(yp) out$y-yp,.drop=FALSE)
  } 
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
  
  
  mindex <- grep('(^y)|(^err)|(^ll)',names(out))
  lindex <- grep('^eta',names(out))
  nosampindex <- which(names(out) %in% c('time','y'))
  
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
      if(i %in% lindex) dr <- fit$ctstanmodelbase$latentNames
      d <- c(d,Row=list(dr))
      if(length(dim(out[[i]])) > 3) d <- c(d,Col=list(dr))
    }
    
    dimnames(out[[i]]) <- d
  }
  out$id <- fit$standata$subject

  return(out)
}



