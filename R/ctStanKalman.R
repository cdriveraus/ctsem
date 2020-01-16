#' Get Kalman filter estimates from a ctStanFit object
#'
#' @param fit fit object from \code{\link{ctStanFit}}.
#' @param nsamples either NA (to extract all) or a positive integer from 1 to maximum samples in the fit.
#' @param cores Integer number of cpu cores to use. Only needed if savescores was set to FALSE when fitting.
#' @param collapsefunc function to apply over samples, such as \code{mean}
#' @param standardisederrors If TRUE, computes standardised errors for prior, upd, smooth conditions.
#' @param subjectpars if TRUE, state estimates are not returned, instead, predictions of each subjects parameters
#' are returned, for parameters that had random effects specified.
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
ctStanKalman <- function(fit,nsamples=NA,collapsefunc=NA,cores=2,standardisederrors=FALSE, subjectpars=FALSE,...){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
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
  
  if(subjectpars){
    p=standata$indvaryingindex
    p=p[p > standata$nlatent]
    if(length(p)==0) stop('No random effects have been integrated over! Did you specify individual variation, and did you optimize
    or set intoverpop = TRUE ?')
    lastsub <- rep(FALSE,standata$ndatapoints)
    for(i in 2:standata$ndatapoints){
      if(standata$subject[i] != standata$subject[i-1]) lastsub[i-1] <- TRUE
      if(i == standata$ndatapoints) lastsub[i] <- TRUE
    }
    pars <- e$etaprior[,lastsub,p,drop=FALSE]
    dimnames(pars) <- list(iter=1:dim(pars)[1],id = unique(standata$subject), par = 1:dim(pars)[3])
    ms <- fit$setup$matsetup
    mv <- fit$setup$matvalues
    mats <- ctStanMatricesList()$all
    for(i in 1:nrow(ms)){
      if(ms$when[i] > 0 && ms$param[i] %in% p){
        if(ms$transform[i] >=0){ #havent done custom transforms here
        pars[,,ms$param[i]-standata$nlatent] <- tform(
        param = pars[,,ms$param[i]-standata$nlatent], transform = ms$transform[i],
          multiplier = mv$multiplier[i],meanscale = mv$meanscale[i],offset = mv$offset[i],
          inneroffset = mv$inneroffset[i])
        dimnames(pars)$par[ms$param[i]-standata$nlatent] <- paste0(names(mats[ms$matrix[i]]),'[',ms$row[i],',',ms$col[i],']')
        }
        if(ms$transform[i] <0) message('custom transforms not implemented yet...sorry')
      }
    }
    return(pars)
  } else{

  nlatent <- fit$standata$nlatent
  nmanifest <- fit$standata$nmanifest

  
  y=matrix(fit$standata$Y,ncol=ncol(fit$standata$Y),dimnames = list(NULL,fit$ctstanmodel$manifestNames))
  y[y==99999] <- NA
  
  out=list(time=cbind(fit$standata$time), 
    y=y, 
    llrow=e$llrow)
  
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
  # 
  if(standardisederrors){
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
  }
  
  
  mindex <- grep('(^y)|(^err)|(^ll)',names(out))
  lindex <- grep('^eta',names(out))
  nosampindex <- which(names(out) %in% c('time','y'))
  out$llrow <- matrix(out$llrow,dim(out$llrow)[1],dim(out$llrow)[2])
  
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
}



