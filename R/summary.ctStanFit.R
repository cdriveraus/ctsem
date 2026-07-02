priorchecker <- function(sf,pars=c('rawpopmeans','rawpopsdbase','tipredeffectparams'),digits=2){
  e=ctExtract(sf)
  funcs <- c(base::mean,stats::sd)
  pars=unlist(lapply(pars,function(x) if(!is.null(dim(e[[x]]))) x))
  out=round(do.call(cbind,lapply(funcs, function(fn) do.call(c,
    lapply(pars, function(obji) apply(e[[obji]],2,fn,na.rm=TRUE)) ))),digits)
  rownames(out)=do.call(c,c(lapply(pars, function(obji) paste0(obji,'_',1:ncol(e[[obji]])))))
  out=data.frame(out,do.call(c,c(lapply(pars, function(obji) 1:ncol(e[[obji]])))))
  out=data.frame(out,do.call(c,c(lapply(pars, function(obji) rep(obji, ncol(e[[obji]]))))),stringsAsFactors = FALSE)
  colnames(out)=c('mean','sd', 'param', 'object')
  
  
  getparnamesfromraw <- function(priorcheck, sf){
    newnames=rownames(priorcheck)
    for(ni in 1:nrow(priorcheck)){
      if(priorcheck$object[ni] %in% 'rawpopmeans'){
        newnames[ni]=paste0('rawpop_',sf$setup$popsetup$parname[sf$setup$popsetup$param %in% priorcheck$param[ni]][1])
      }
      if(priorcheck$object[ni] %in% 'tipredeffectparams'){
        newnames[ni]=paste0('rawtipredeffect_',paste0(
          which(sf$standata$TIPREDEFFECTsetup == priorcheck$param[ni],arr.ind = TRUE),collapse='_'))
      }
    }
    return(newnames)
  }
  
  rownames(out) = getparnamesfromraw(priorcheck=out,sf=sf)
  return(out)
}

ctFitgetparnamesfromraw <- function(sf){
  names=list()
  ms <- sf$setup$matsetup
  names$popmeans <- unique(ms$parname[ms$when==0 & ms$param > 0 & ms$copyrow ==0])
  if(sum(sf$ctstanmodelbase$pars$indvarying) > 0){
    names$popsd <- paste0(ms$parname[ms$when==0 & ms$param > 0 & ms$copyrow == 0 & ms$indvarying > 0],'_SD')
  } 
  if(sum(sf$ctstanmodelbase$pars$indvarying) > 1){
    popcorrmat <- matrix(paste0(gsub('_SD','',rep(names$popsd,times=length(names$popsd))), '_', gsub('_SD','',rep(names$popsd,each=length(names$popsd))),'_corr'),ncol=length(names$popsd),nrow=length(names$popsd))
    names$popcorr <- as.vector(popcorrmat[lower.tri(popcorrmat)])
  } 
  
  if(sf$ctstanmodelbase$n.TIpred > 0){
    names$tipreds <- sapply(sf$standata$TIPREDEFFECTsetup[sf$standata$TIPREDEFFECTsetup!=0],function(x){
      index <- which(sf$standata$TIPREDEFFECTsetup == x,arr.ind = TRUE)
      paste0('rawtipredeffect_',paste0(
        names$popmeans[index[1]],'_',sf$ctstanmodelbase$TIpredNames[index[2]]))
    })
  }
  return(names)
  
}


priorcheckreport <- function(sf, meanlim = 2, sdlim= .2,digits=2){
  p=priorchecker(sf)
  ps=sf$setup$popsetup
  p=p[abs(p$mean) > meanlim | p$sd > sdlim,]
  out<-p[,c('mean','sd')]
  return(out)
}



ctStanRawSamples<-function(fit){
  if(length(fit$stanfit$stanfit@sim)==0) {
    samples = fit$stanfit$rawposterior
  } else {
    samples = t(stan_unconstrainsamples(fit$stanfit$stanfit,fit$standata))
  }
  return(samples)
}



#' ctSummaryMatrices
#'
#' Summarise model-implied parameter matrices from a ctsem fit object
#'
#'@param fit fit object from \code{\link{ctFit}}
#'@param calcfunc Function to apply over samples, must return a single value. 
#'By default the median over all samples is returned using the \code{\link[stats]{quantile}} function, 
#'but one might also be interested in the \code{\link[base]{mean}} or \code{\link[stats]{sd}}, for instance.
#'@param calcfuncargs A list of additional parameters to pass to calcfunc. 
#'For instance, with the default of calcfunc = quantile, 
#'the probs argument is needed to ensure only a single value is returned.
#'@param timeinterval time interval for discrete time parameter matrix computation.
#'@examples
#'\donttest{
#'#posterior median over all subjects (also reflects mean of unconstrained pars)
#'ctSummaryMatrices(ctstantestfit)
#'}
#'@export
ctSummaryMatrices <- function(fit,
  calcfunc=quantile,calcfuncargs=list(probs=0.5),timeinterval=1){
  
  if(!'ctStanFit' %in% class(fit)) stop(paste0('Not an object of class ctStanFit! Instead is ',paste0(class(fit),collapse=', ')))
  
  e<-ctExtract(fit,cores=1) #Qfit$stanfit$transformedpars #first dim of subobjects is iter, 2nd subjects
  niter=dim(e$pop_DRIFT)[1]
  
  
  
  
  mats <- ctStanMatricesList()
  mats <- c(names(mats$base), names(mats$asymptotic),names(mats$extra))
  if(fit$ctstanmodel$continuoustime){
    d=list(DRIFT=e$pop_DRIFT)
    dd=ctDiscreteParsDrift(d,timeinterval,observational = FALSE,standardise = FALSE,cov = FALSE,quiet=TRUE)
    e$pop_dtDRIFT <- array(dd,dim=dim(dd)[-2:-3])
    mats <- c(mats, 'dtDRIFT')
  }
  
  out <- list()
  for(matname in (mats)){
    try({
      calcfuncargs$collapsemargin = 1
      calcfuncargs$collapsefunc=calcfunc
      calcfuncargs$na.rm=TRUE
      
      calcfuncargs$inarray = e[[paste0('pop_',matname)]]
      out[[matname]] <- array(do.call(ctCollapse,calcfuncargs),
        dim=dim(calcfuncargs$inarray)[-1])
    },silent=TRUE)
  }
  
  if(nrow(out$T0MEANS) > nrow(out$CINT)){ #then intoverpop used...
    nlatent <- nrow(out$CINT)
    out$T0MEANS <- out$T0MEANS[1:nlatent,1,drop=FALSE]
    out$DRIFT <- out$DRIFT[1:nlatent,1:nlatent,drop=FALSE]
    out$T0VAR <- out$T0VAR[1:nlatent,1:nlatent,drop=FALSE]
    out$T0cov <- out$T0cov[1:nlatent,1:nlatent,drop=FALSE]
  }
  
  ln=fit$ctstanmodel$latentNames
  mn=fit$ctstanmodel$manifestNames
  tdn=fit$ctstanmodel$TDpredNames
  dimnames(out$DRIFT)=list(ln,ln)
  dimnames(out$DIFFUSIONcov)=list(ln,ln)
  dimnames(out$DIFFUSION)=list(ln,ln)
  dimnames(out$T0cov)=list(ln,ln)
  dimnames(out$asymDIFFUSIONcov)=list(ln,ln)
  rownames(out$CINT)=ln
  rownames(out$MANIFESTMEANS)=mn
  rownames(out$T0MEANS)=ln
  
  dimnames(out$T0VAR)=list(ln,ln)
  dimnames(out$LAMBDA)=list(mn,ln)
  
  
  if(!is.null(e$pop_MANIFESTVAR)) {
    dimnames(out$MANIFESTVAR)=list(mn,mn)
    dimnames(out$MANIFESTcov)=list(mn,mn)
    # out$MANIFESTVAR=out$MANIFESTVAR %*% t(out$MANIFESTVAR) #cholesky factor inside stanfit...
    
  }
  
  if(!is.null(e$pop_TDPREDEFFECT)) {
    dimnames(out$TDPREDEFFECT)=list(ln,tdn)
  }
  
  out$MANIFESTVAR <- NULL ; 
  
  
  return(out)
}

#' Backward-compatible alias for \code{ctSummaryMatrices}.
#' @rdname ctSummaryMatrices
#' @export
ctStanContinuousPars <- ctSummaryMatrices




#' Extract an array of subject specific parameters from a ctStanFit object.
#'
#' @param fit fit object
#' @param pointest if TRUE, returns only the set of individual difference parameters
#' based on the max a posteriori estimate (or the median if sampling approaches were used).
#' @param cores Number of cores to use.
#' @param nsamples Number of samples to calculate parameters for. Not used if pointest=TRUE.
#' 
#' @details This function returns the estimates of individual parameters, taking into account any
#' covariates and random effects. 
#'
#' @return an nsamples by nsubjects by nparams array.
#' @export
#'
#' @examples
#' indpars <- ctSubjectPars(ctstantestfit)
#' dimnames(indpars)
#' plot(indpars[1,,'cint1'],indpars[1,,'cint2'])
ctSubjectPars <- function(fit,pointest=TRUE,cores=2,nsamples='all'){
  
  if(!nsamples[1] %in% 'all') fit$stanfit$rawposterior <- 
      fit$stanfit$rawposterior[sample(1:nrow(fit$stanfit$rawposterior),nsamples),,drop=FALSE]
  pnames <- getparnames(fit,subjvariationonly = TRUE)
  if(length(pnames)==0) stop('No individually varying parameters in model!')
  m <- fit$ctstanmodelbase$pars
  if(pointest) tfp <- fit$stanfit$transformedparsfull else {
    gc()
    tfp <- ctExtract(fit,subjectMatrices = TRUE,cores=cores,nsamples=nsamples)
  }
  p <- array(NA, dim=c(dim(tfp$pop_DRIFT)[1],fit$standata$nsubjects,length(pnames)))
  dimnames(p) <- list(iter=1:dim(p)[1],subject=1:dim(p)[2],param=1:dim(p)[3])
  co<-0
  for(i in 1:nrow(m)){
    if(!is.na(m$param[i]) & is.na(m$value[i]) & #if a free param, not complex and not already collected...
        !grepl('[',m$param[i],fixed=TRUE) & !m$param[i] %in% dimnames(p)[[3]]){
      if(m$param[i] %in% pnames){ #if ind differences
        co = co+1
        p[,,co] <- tfp[[paste0('subj_',m$matrix[i])]] [
          , ,m$row[i],m$col[i]]
        dimnames(p)[[3]][co] <- m$param[i]
      }
    }
  }
  p=p[,,order(dimnames(p)[[3]]),drop=FALSE]
  return(p)
}

#' Backward-compatible alias for \code{ctSubjectPars}.
#' @rdname ctSubjectPars
#' @export
ctStanSubjectPars <- ctSubjectPars


getparnames <- function(fit,reonly=FALSE, subjvariationonly=FALSE, popstatesonly=FALSE){
  ms <- fit$setup$matsetup
  
  if(popstatesonly)  indices=ms$param > 0 & ms$copyrow <1 & ms$matrix==1 & ms$indvarying > 0 & ms$row > fit$standata$nlatent
  if(!popstatesonly)  indices=ms$when %in% c(0,-1) & ms$param > 0 & ms$copyrow < 1 & !grepl('[',ms$parname,fixed=TRUE)
  if(subjvariationonly) indices = ms$when %in% c(0,-1) & ms$param > 0 & ms$copyrow < 1 & (ms$tipred >0 | ms$indvarying > 0)
  pars <- data.frame(parnames = ms$parname[indices],  parindices = ms$param[indices])
  
  pars<-pars[!duplicated(pars$parnames),]
  pars<-pars[order(pars$parindices),]
  parnames <- pars[pars$parindices >0, 1]
  if(reonly)  parnames <- parnames[fit$standata$indvaryingindex]
  
  return(parnames)
}

#' summary.ctStanFit
#'
#' Summarise a ctStanFit object that was fit using \code{\link{ctFit}}.
#' 
#' @param object fit object from \code{\link{ctFit}}, of class ctStanFit.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter calculations
#' reported in summary. 
#' @param digits integer denoting number of digits to report.
#' @param parmatrices if TRUE, also return additional parameter matrices -- can be slow to compute
#' for large models with many samples.
#' @param priorcheck Whether or not to use \code{ctsem:::priorchecking} to compare posterior mean and sd to prior mean and sd.
#' @param residualcov Whether or not to show standardised residual covariance. Takes a little longer to compute.
#' @param ... Additional arguments to pass to \code{ctsem:::priorcheckreport}, such as \code{meanlim}, or \code{sdlim}.
#' @return List containing summary items, with a \code{print} method for readable console and knitr output.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=3,parmatrices=TRUE,priorcheck=TRUE,residualcov = TRUE,...){
  
  if(!'ctStanFit' %in% class(object)) stop('Not a ctStanFit object!')
  
  monitor <- function(dat,...){
    out <- rstan::monitor(dat,...)
    class(out) <- 'data.frame'
    return(out)
  }
  
  out=list()
  monvars <- c('mean','sd','2.5%','50%','97.5%')
  
  optimize=length(object$stanfit$stanfit@sim)==0
  
  if(!optimize){ 
    smr<-suppressWarnings(getMethod('summary','stanfit')(object$stanfit$stanfit))
    if('98%' %in% colnames(smr$summary)) colnames(smr$summary)[colnames(smr$summary)=='98%'] <- '97.5%'
  }
  
  
  e <- ctExtract(object) 
  
  if(residualcov){ #cov of residuals
    obscov <- cov(object$data$Y,use='pairwise.complete.obs')
    idobscov <- diag(1/sqrt(diag(obscov)),ncol(obscov))
    rescov <- cov(matrix(object$stanfit$kalman$errprior,ncol=ncol(obscov)),use='pairwise.complete.obs')
    narescov <- which(is.na(rescov))
    rescov[narescov] <- 0
    
    out$residCovStd <- round(idobscov %*% rescov %*% idobscov ,digits)
    out$residCovStd[narescov] <- NA
    dimnames(out$residCovStd) <- list(object$ctstanmodel$manifestNames,object$ctstanmodel$manifestNames)
  }
  
  ms=object$setup$matsetup
  parnames = getparnames(object)
  # parnames <- unique(parnames)
  parnamesiv <- getparnames(object,reonly = TRUE)
  
  #### generate covcor matrices of raw and transformed subject level params
  
  iter=dim(e$rawpopcov)[1]
  if(!is.null(iter)){ #then there is some individual variation so continue
    nindvarying=dim(e$rawpopcov)[2]
    
    if(nindvarying>1){
      
      #raw pop distribution params
      dimrawpopcorr <- dim(e$rawpopcorr)
      rawpopcorr= array(e$rawpopcorr,dim=c(dimrawpopcorr[1],1,dimrawpopcorr[2] * dimrawpopcorr[3]))
      rawpopcorrout <- suppressWarnings(monitor(rawpopcorr, digits_summary=digits,warmup=0,print = FALSE)[lower.tri(diag(nindvarying)),c(monvars,'n_eff','Rhat'),drop=FALSE])
      if(optimize) rawpopcorrout <- rawpopcorrout[,-which(colnames(rawpopcorrout) %in% c('n_eff','Rhat')),drop=FALSE]
      
      rownames(rawpopcorrout) <- matrix(paste0('',parnamesiv,'__',rep(parnamesiv,each=length(parnamesiv))),
        length(parnamesiv),length(parnamesiv))[lower.tri(diag(nindvarying)),drop=FALSE]
      
      rawpopcorrout <- cbind(rawpopcorrout,rawpopcorrout[,'mean'] / rawpopcorrout[,'sd'])
      colnames(rawpopcorrout)[ncol(rawpopcorrout)] <- 'z'
      
      out$rawpopcorr = round(rawpopcorrout,digits)
      out$rawpopcorrNote = 'These reflect correlations between the raw / unconstrained parameters.'
    }
  }
  
  if(priorcheck && object$standata$priors) {
    priorcheckres <- priorcheckreport(object,...)
    if(nrow(priorcheckres) > 0){
      out$priorcheck=priorcheckres
      out$priorcheckNote='These posteriors exceeded arbitrary limits re normal(0,1) -- priors / transforms are likely somewhat informative (Not necessarily a problem).'
    }
  }
  
  if(object$ctstanmodel$n.TIpred > 0) {
    
    if(!optimize){
      rawtieffect <- rstan::extract(object$stanfit$stanfit,permuted=FALSE,pars='TIPREDEFFECT')
      tidiags <- suppressWarnings(monitor(rawtieffect,warmup=0,digits_summary = digits,print = FALSE))
    }
    
    tieffect <- array(e$linearTIPREDEFFECT,dim=c(dim(e$linearTIPREDEFFECT)[1], 1, length(parnames) * dim(e$linearTIPREDEFFECT)[3]))
    tieffectnames <- paste0('tip_',rep(object$ctstanmodel$TIpredNames,each=length(parnames)),'_',parnames)
    dimnames(tieffect)<-list(c(),c(),tieffectnames)
    tipreds = suppressWarnings(monitor(tieffect,warmup = 0,print = FALSE)[,monvars,drop=FALSE])
    if(!optimize) tipreds <- cbind(tipreds,tidiags[,c('n_eff','Rhat'),drop=FALSE])
    tipreds <- tipreds[c(object$data$TIPREDEFFECTsetup)>0,,drop=FALSE]
    z = tipreds[,'mean'] / tipreds[,'sd'] 
    out$tipreds= round(cbind(tipreds,z),digits) #[order(abs(z)),]
    out$tipredsNote = 'Approximate (linearised) effects on the transformed parameters.'
  }
  
  
  
  if(parmatrices){
    Mean=ctSummaryMatrices(object,calcfunc = mean,calcfuncargs=list(),timeinterval=timeinterval)
    sd=ctSummaryMatrices(object,calcfunc = sd,calcfuncargs = list(na.rm=TRUE),timeinterval=timeinterval)
    `2.5%` = ctSummaryMatrices(object,calcfunc = quantile,calcfuncargs = list(probs=.025),timeinterval=timeinterval)
    `50%` = ctSummaryMatrices(object,calcfunc = quantile,calcfuncargs = list(probs=.5),timeinterval=timeinterval)
    `97.5%` = ctSummaryMatrices(object,calcfunc = quantile,calcfuncargs = list(probs=.975),timeinterval=timeinterval)
    
    d <- data.frame(ctModelUnlist(Mean,matnames = names(Mean)))
    colnames(d)[colnames(d) %in% 'value'] <- 'Mean'
    lapply(c('sd','2.5%','50%','97.5%'),function(x){
      d[[x]] <<-round(ctModelUnlist(get(x),names(Mean))$value,digits)
      })
    d$param <- NULL
    d$Mean <- round(d$Mean,digits)
    rm(sd)
    d <- d[!d$matrix %in% c('DIFFUSION','T0VAR'),]
    
    out$parmatrices=d
  }
  
  
  
  if(!optimize){
    # browser()
    popsd=smr$summary[c(grep('^popsd',rownames(smr$summary),fixed=FALSE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE] #[ object$data$indvaryingindex,,drop=FALSE]
    rownames(popsd)=parnames[ object$data$indvaryingindex]
    popmeans=smr$summary[c(grep('popmeans[', rownames(smr$summary),fixed=TRUE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=popmeans[(nrow(popmeans)/2+1):nrow(popmeans),,drop=FALSE]
    rownames(popmeans) <- parnames
    
    loglik = data.frame(mean=mean(e$ll),sd=sd(e$ll), max=max(e$ll))
    logposterior=smr$summary[c(grep('lp',rownames(smr$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
  }
  
  if(optimize){ #if optimized / importance sampled
    
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
  
  out$popNote=paste0('covariance pars in sd / unconstrained cor form, see System Matrices (or ctSummaryMatrices function) for cor/cov.')
  
  out$logposterior=logposterior
  if(optimize) {
    out$loglik=loglik
    out$npars = npars
    out$aic = aic
  }

  if(optimize) out$nsamples <- nrow(object$stanfit$samples)
  
  if(!parmatrices) out$parmatNote <- 'For additional summary matrices, use argument: parmatrices = TRUE'
  
  out <- lapply(out,function(x){
    if('matrix' %in% class(x)){
      x <- data.frame(x,check.names=FALSE)
    }
    x <- roundSummaryCtStanFitValue(x,digits=digits)
    x})
  
  attr(out,'digits') <- digits
  class(out) <- 'summary.ctStanFit'
  return(out)
}

roundSummaryCtStanFitValue <- function(x,digits){
  if(is.data.frame(x)){
    for(col in seq_along(x)){
      if(is.numeric(x[[col]])) x[[col]] <- round(x[[col]],digits)
    }
    return(x)
  }
  
  if(is.numeric(x)) return(round(x,digits))
  x
}

summaryCtStanFitLabel <- function(x){
  labels <- c(
    residCovStd = 'Standardised residual covariance',
    rawpopcorr = 'Random-effects correlations',
    rawpopcorrNote = 'Note',
    priorcheck = 'Prior check',
    priorcheckNote = 'Note',
    tipreds = 'Time-independent predictor effects',
    tipredsNote = 'Note',
    parmatrices = 'System Matrices',
    popsd = 'Random-effects standard deviations',
    popmeans = 'Fixed-effects / Population means',
    popNote = 'Note',
    loglik = 'Log likelihood',
    npars = 'Number of parameters',
    aic = 'AIC',
    logposterior = 'Log posterior',
    nsamples = 'Number of samples',
    parmatNote = 'Summary matrices note')
  if(x %in% names(labels)) return(unname(labels[x]))
  gsub('([a-z])([A-Z])', '\\1 \\2', x)
}

printSummaryCtStanFitValue <- function(x,width,row.names,...){
  if(is.character(x) && is.null(dim(x))){
    wrapped <- unlist(strwrap(x,width=width),use.names=FALSE)
    if(length(wrapped)) cat(paste(wrapped,collapse='\n'),'\n',sep='')
    return(invisible(x))
  }
  
  if(is.data.frame(x)){
    print(x,row.names=row.names,...)
    return(invisible(x))
  }
  
  print(x,...)
  invisible(x)
}

summaryCtStanFitIsScalar <- function(x){
  is.atomic(x) && length(x)==1
}

summaryCtStanFitIsNote <- function(section){
  grepl('Note$',section)
}

printSummaryCtStanFitScalar <- function(label,x,width){
  prefix <- paste0(label,': ')
  if(is.character(x)){
    wrapped <- unlist(strwrap(x,width=width,initial=prefix,exdent=nchar(prefix)),use.names=FALSE)
    if(length(wrapped)) cat(paste(wrapped,collapse='\n'),'\n',sep='')
    return(invisible(x))
  }
  
  cat(prefix,format(as.vector(x),trim=TRUE),'\n',sep='')
  invisible(x)
}

#' Print ctStanFit summaries
#'
#' @param x Object returned by \code{\link{summary.ctStanFit}}.
#' @param width Console width to use for wrapping note text and formatting tables.
#' @param sections Optional character vector of summary sections to print.
#' @param row.names Logical. Print row names for data frame sections?
#' @param ... Additional arguments passed to \code{\link{print}} for individual sections.
#' @return Invisibly returns \code{x}.
#' @method print summary.ctStanFit
#' @export
print.summary.ctStanFit <- function(x,width=getOption('width'),sections=names(x),row.names=TRUE,...){
  width <- suppressWarnings(as.integer(width)[1])
  if(is.na(width)) width <- getOption('width')
  width <- max(20,width)
  oldwidth <- getOption('width')
  options(width=width)
  on.exit(options(width=oldwidth),add=TRUE)
  
  sections <- intersect(sections,names(x))
  for(section in sections){
    label <- summaryCtStanFitLabel(section)
    if(summaryCtStanFitIsScalar(x[[section]])){
      if(summaryCtStanFitIsNote(section)) cat('\n')
      printSummaryCtStanFitScalar(label,x[[section]],width=width)
      if(summaryCtStanFitIsNote(section)) {
        cat(paste(rep('-',min(width,72)),collapse=''),'\n',sep='')
      }
      next
    }
    cat('\n',label,'\n',paste(rep('-',nchar(label)),collapse=''),'\n',sep='')
    printSummaryCtStanFitValue(x[[section]],width=width,row.names=row.names,...)
  }
  invisible(x)
}
