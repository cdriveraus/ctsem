#' Append pseudo-posterior samples to an optimized ctStanFit object (deprecated)
#'
#' Appends normal draws from the fit's already-computed covariance to its
#' raw posterior. Deprecated in favour of \code{\link{ctOptimUncertainty}} with
#' \code{uncertainty = 'stored'}, which does the same on both backends.
#'
#' @param fit fit object
#' @param nsamples number of samples desired
#' @param cores number of cores to use
#'
#' @return fit object with extra samples
#' @aliases ctAddSamples
#' @details These are pseudo-posterior draws from the fitted covariance, not
#'   posterior draws from a sampler. \code{\link{ctSample}} is the latter --
#'   Hamiltonian Monte Carlo from an optimized \code{ctJuliaFit} -- and is a
#'   different object, not another route to this one.
#'
#'   Deprecated in favour of \code{ctOptimUncertainty(fit, uncertainty =
#'   'stored', finishsamples = n)}, which draws from the same covariance for
#'   the same zero model evaluations, works on a \code{ctJuliaFit} as well as a
#'   \code{ctStanFit}, and records the new sample count in
#'   \code{$uncertainty$settings}. It also replaces the draws rather than
#'   appending to them: appending normal draws to a posterior that came from
#'   \code{uncertainty='is'} or \code{'bootstrap'} leaves part of each in one
#'   matrix with nothing recording the mixture, which is what this function
#'   does.
#'
#'   \code{ctAddSamples} is the same function under its pre-3.11 name; both are
#'   deprecated together.
#' @seealso \code{\link{ctOptimUncertainty}}, \code{\link{ctSample}}
#' @export
#'
#' @examples
#' \dontrun{
#' newfit <- ctOptimUncertainty(ctstantestfit, uncertainty = 'stored',
#'   finishsamples = 30, cores = 1)
#' }
ctFitAddSamples <- function(fit,nsamples,cores=2){

  # `.Deprecated(msg=)` rather than a bare `warning()`: it carries the same
  # text but as a `deprecatedWarning`, so a caller can silence or catch it by
  # class. The body below is deliberately unchanged -- it draws in a different
  # RNG order from `ctOptimNormalDraws()` (one `rnorm(npar)` per sample rather
  # than one `rnorm(n*npar)` filled by column), so routing it through the
  # replacement would move every number this has ever produced for a given
  # seed.
  .Deprecated(msg = paste0(
    "ctFitAddSamples() is deprecated. Use ctOptimUncertainty(fit, ",
    "uncertainty = 'stored', finishsamples = n), which redraws from the same ",
    "covariance and works on both backends."))

  if(length(fit$stanfit$stanfit@sim) > 0) stop('ctStanFit object was sampled and not optimized, cannot add samples!')

  mchol <- t(chol(fit$stanfit$cov))
  resamples <- matrix(unlist(lapply(1:nsamples,function(x){
    fit$stanfit$rawest + (mchol) %*% t(matrix(rnorm(length(fit$stanfit$rawest)),nrow=1))
  } )),byrow=TRUE,ncol=length(fit$stanfit$rawest))
  
  fit$stanfit$rawposterior <- rbind(fit$stanfit$rawposterior,resamples)
  
  fit$stanfit$transformedpars=stan_constrainsamples(sm = fit$stanmodel,
    standata = fit$standata,samples=fit$stanfit$rawposterior,
    savescores = fit$standata$savescores,
    savesubjectmatrices=as.logical(fit$standata$savesubjectmatrices),
    # Either flag needs the filter pass; savesubjectmatrices already forces
    # savescores on at ctFit.R:1176, so savescores is the one to read.
    dokalman=as.logical(fit$standata$savescores),
    cores=cores)
  return(fit)
}

# Pre-3.11 name, documented on ctFitAddSamples' page via @aliases. Kept
# because it is on CRAN; both are deprecated together and both go out with the
# stan backend.
#' @export
ctAddSamples <- ctFitAddSamples

#' Quickly initialise stanfit object from model and data
#'
#' @param model stanmodel
#' @param data standata
#' @param fast Use cut down form for speed
#'
#' @return stanfit object
#' @export
#'
#' @examples
#' sf <- stan_reinitsf(ctstantestfit$stanmodel,ctstantestfit$standata)
stan_reinitsf <- function(model, data,fast=FALSE){
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))
  
  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<- 
      rstan::sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE,
        control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
  
  return(sf)
}

# Function to compute numeric Hessian using finite differences
numericHessianFunc <- function(pars, step=1e-3, whichpars='all',
  lpdifmin=1e-8, lpdifmax=.1, cl=NA, verbose=1, directions=c(-1,1), parsteps=c(), 
  lpgFunc, base_value, base_gradient=NULL) {
  
  if('all' %in% whichpars) whichpars <- 1:length(pars)
  if(is.null(base_gradient)){
    base_eval <- suppressMessages(suppressWarnings(lpgFunc(pars)))
    base_gradient <- attributes(base_eval)$gradient
    if(missing(base_value) || is.null(base_value)) base_value <- base_eval[1]
  }
  if(is.null(base_gradient) || length(base_gradient) != length(pars) ||
      any(!is.finite(base_gradient))) {
    base_gradient <- rep(0, length(pars))
  }
  
  hessout <- sapply(whichpars, function(i){
    
    message(paste0("\rEstimating Hessian, par ",i,',', 
      as.integer(i/length(pars)*50+ifelse(directions[1]==1,0,50)),
      '%'),appendLF = FALSE)
    if(verbose) message('### Par ',i,'###')
    stepsize = step
    uppars<-rep(0,length(pars))
    uppars[i]<-1
    accepted <- FALSE
    count <- 0
    lp <- list()
    steplist <- list()
    for(di in 1:length(directions)){
      count <- 0
      accepted <- FALSE
      stepchange = 0
      stepchangemultiplier = 1
      while(!accepted && (count==0 || 
          ( 
            (count < 30 && any(is.na(attributes(lp[[di]])$gradient))) || #if NA gradient, try for 30 attempts
              (count < 15 && all(!is.na(attributes(lp[[di]])$gradient))))#if gradient ok, stop after 15
      )){ 
        stepchangemultiplier <- max(stepchangemultiplier,.11)
        count <- count + 1
        lp[[di]] <-  suppressMessages(suppressWarnings(lpgFunc(pars+uppars*stepsize*directions[di])))
        accepted <- !'try-error' %in% class(lp[[di]]) && all(!is.na(attributes(lp[[di]])$gradient))
        if(accepted){
          lpdiff <- base_value[1] - lp[[di]][1]
          # if(lpdiff > 1e100) 
          if(lpdiff < lpdifmin) {
            if(verbose) message('Increasing step')
            if(stepchange == -1) stepchangemultiplier = stepchangemultiplier*.5
            stepchange <- 1
            stepsize <- stepsize*(1-stepchangemultiplier)+ (stepsize*10)*stepchangemultiplier
          }
          if(lpdiff > lpdifmax){
            if(verbose) message('Decreasing step')
            
            
            if(stepchange == 1) stepchangemultiplier = stepchangemultiplier * .5
            stepchange <- -1
            stepsize <- stepsize*(1-stepchangemultiplier)+ (stepsize*.1)*stepchangemultiplier
          }
          if(lpdiff > lpdifmin && lpdiff < lpdifmax && lpdiff > 0) accepted <- TRUE else accepted <- FALSE
          if(lpdiff < 0){
            base_value <<- lp[[di]]
            base_gradient <<- attributes(lp[[di]])$gradient
            if(verbose) message('Better log probability found during Hessian estimation...')
            accepted <- FALSE
            stepchangemultiplier <- 1
            stepchange=0
            count <- 0
            di <- 1
          }
        } else stepsize <- stepsize * 1e-3
      }
      if(stepsize < step) step <<- step *.1
      if(stepsize > step) step <<- step *10
      steplist[[di]] <- stepsize
    }
    
    grad<- (attributes(lp[[1]])$gradient - base_gradient) /
      steplist[[1]] * directions[1]
    if(any(is.na(grad))){
      warning('NA gradient encountered at param ',i,immediate. =TRUE)
    }
    if(length(directions) > 1) grad <- (grad +
        (attributes(lp[[2]])$gradient - base_gradient) /
          (steplist[[2]] * directions[2]))/2
    return(grad)
  }) #end sapply
  
  out=(hessout+t(hessout))/2
  return(out)
}

# Function to process and clean Hessian matrices
processHessianMatrices <- function(hess1, hess2, verbose, matsetup) {
  hess <- hess1
  hess[is.na(hess)] <- 0 #set hess1 NA's to 0
  hess[!is.na(hess2)] <- hess[!is.na(hess2)] + hess2[!is.na(hess2)] #add hess2 non NA's
  hess[!is.na(hess1) & !is.na(hess2)] <- hess[!is.na(hess1) & !is.na(hess2)] /2 #divide items where both hess1 and 2 used by 2
  hess[is.na(hess1) & is.na(hess2)] <- NA #set NA when both hess1 and 2 NA
  
  onesided <- which(xor(is.na(diag(hess1)), is.na(diag(hess2))))
  if(any(is.na(c(diag(hess1),diag(hess2))))){
    if(any(is.na(hess))) message ('Problems computing Hessian...')
  }
  
  #make symmetric
  hess <- (t(hess)+hess)/2
  hess[upper.tri(hess)] <- t(hess)[upper.tri(hess)] 
  probpars <- which(is.finite(diag(hess)) & diag(hess) > -1e-6)
  
  if(length(c(probpars,onesided)) > 0){
    parlabels <- function(ii){
      if(length(ii) < 1) return(character())
      labels <- as.character(ii)
      if(('data.frame' %in% class(matsetup) || !all(is.na(matsetup[1])))) {
        ms=matsetup
        ms=ms[ms$param > 0 & ms$when == 0,]
        ms=ms[!duplicated(ms$param),]
        matched <- ms$parname[match(ii, ms$param)]
        labels[!is.na(matched)] <- matched[!is.na(matched)]
      }
      labels
    }
    if('data.frame' %in% class(matsetup) || !all(is.na(matsetup[1]))){
      ms=matsetup
      ms=ms[ms$param > 0 & ms$when == 0,]
      ms=ms[!duplicated(ms$param),]
    }
    if(length(onesided) > 0){
      message ('One sided Hessian used for params: ',
        paste0(parlabels(onesided),collapse=', '))
    }
    if(length(probpars) > 0){
      message('***These params "may" be not identified: ',
        paste0(parlabels(probpars),collapse=', '))
    }
  }
  
  return(list(hess = hess, probpars = probpars, onesided = onesided))
}


# =============================================================================
# ADDITIONAL HELPER FUNCTIONS
# =============================================================================

# PARALLEL AND UTILITY FUNCTIONS
# =============================================================================

parallelStanSetup <- function(cl, standata,split=TRUE,nsubsets=1,smfile=NA){
  cores <- length(cl)
  if(split) stanindices <- split(unique(standata$subject),(unique(standata$subject) %% min(standata$nsubjects,cores))) #disabled sorting so subset works parallel
  if(!split) stanindices <- lapply(1:cores,function(x) unique(standata$subject))
  if(length(stanindices) < cores){
    for(i in (length(stanindices)+1):cores){
      stanindices[[i]] <- NA
    }
  }
  
  standata$nsubsets <- as.integer(nsubsets)
  if(!split) cores <- 1 #for prior mod
  
  parallel::clusterExport(cl,c('standata','stanindices','cores'),envir=environment())
  
  #if smfile not exported and recompile needed, export it
  if(standata$recompile > 0 && !all(unlist(parallel::clusterEvalQ(cl,{exists('smfile')})))){ 
    parallel::clusterExport(cl,'smfile',envir=environment())
  }
  
  parallel::clusterEvalQ(cl,{
    # g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    # environment(parlptext) <- g
    if(standata$recompile > 0) load(file=smfile) else sm <- utils::getFromNamespace("stanmodels", "ctsem")$ctsm
    # eval(parse(text=parlptext))
    # assign("parlp",parlp,pos=g)
    if(FALSE){ sm=99;smfile=NULL;nodeid=NULL} #global variables
    parlp <- function(parm){
      a=Sys.time()
      out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
      if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
        outerr <- out
        out <- -1e100
        attributes(out)$gradient <- rep(NaN, length(parm))
        attributes(out)$err <- outerr
      }
      if(any(is.infinite(attributes(out)$gradient))) { #clip infinite gradients
        attributes(out)$gradient[is.infinite(attributes(out)$gradient)] <- 1000 * sign(attributes(out)$gradient[is.infinite(attributes(out)$gradient)])
      }
      attributes(out)$time <- Sys.time()-a
      if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
      return(out)
    }
    if(length(stanindices[[nodeid]]) < length(unique(standata$subject))) standata <-  utils::getFromNamespace("standatact_specificsubjects", "ctsem")(standata,stanindices[[nodeid]])
    standata$priormod <- 1/cores
    
    smf=utils::getFromNamespace("stan_reinitsf", "ctsem")(sm,standata)
  })
  NULL
}

singlecoreStanSetup <-function(standata, nsubsets,sm){
  cores <- 1
  standata$nsubsets <- as.integer(nsubsets)
  # if(!is.null(standata$recompile)) standata$recompile <- 0 #no recompile on single core
  if(standata$recompile == 0) smf <- stan_reinitsf(stanmodels$ctsm,standata)
  if(standata$recompile > 0) smf <- stan_reinitsf(sm,standata)
  return(eval(parse(text=parlptext)))#create parlp function
}

#create as text because of parallel communication weirdness
parlptext <- 
  'parlp <- function(parm){
     a=Sys.time()
          out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
          outerr <- out
          out <- -1e100
          attributes(out)$gradient <- rep(NaN, length(parm))
          attributes(out)$err <- outerr
        }
        attributes(out)$time <- Sys.time()-a
        if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
        return(out)
        }'

#based on rstan function, very cut down, may fail in some cases...
#' @importFrom Rcpp cpp_object_initializer
getcxxfun <- function(object) {
  if (length(object@dso_saved) == 0){
    return(function(...) stop("this function should not be called"))
  }  else  return(object@.CXXDSOMISC$cxxfun)
}

flexlapply <- function(cl, X, fn,cores=1,...){
  if(cores > 1) parallel::parLapply(cl,X,fn,...) else lapply(X, fn,...)
}

# Results come back in input order, and a NULL stays a NULL.
#
# Neither held at `cores > 1`, and together they ended whole calls. Each worker
# evaluates `nodeindices[[nodeid]]`, and `nodeid` is a rank drawn inside
# `makeClusterID()` rather than the node's position in `cl`, so the
# concatenated results arrive in an order the caller cannot reconstruct.
# `unlist()` then *drops* the NULL a failed element returns, so the result was
# both permuted and short, with nothing left to say which elements were
# missing.
#
# `stan_constrainsamples()` is the only caller, and it reads both properties:
# it takes `which(!nulls)[1]` as a row of `samples` to build the skeleton from.
# With the NULLs dropped that index is always 1, so an inadmissible first draw
# was handed straight to an unprotected `rstan::constrain_pars()` and its
# exception -- `quad_form_sym: A is not symmetric. A[1,2] = inf`, or
# `mdivide_left_spd: Matrix A is not positive definite` -- ended the call.
# Intermittent, because it needed draw 1 in particular to be inadmissible, and
# invisible at `cores = 1`, where `lapply()` keeps the NULLs and the two
# indices agree. That is the \donttest example on ctGenerateFromPriors failing
# under `R CMD check --as-cran`, and no seed suppresses it.
#
# Carrying the index alongside each result fixes both, and leaves the two
# branches returning the same shape. `ordered[i] <- list(v)`, not
# `ordered[[i]] <- v`: the latter *deletes* element i when v is NULL, which is
# exactly the case being preserved.
flexlapplytext <- function(cl, X, fn,cores=1,...){
  if(cores > 1) {
    nodeindices <- split(1:length(X), sort((1:length(X))%%cores))
    nodeindices<-nodeindices[1:cores]
    clusterIDexport(cl,c('nodeindices'))

    out <-unlist(clusterIDeval(cl,paste0(
      'lapply(nodeindices[[nodeid]], function(.i) list(.idx = .i, .val = ',fn,'(.i)))')),
      recursive = FALSE)
    ordered <- vector('list', length(X))
    for(el in out) ordered[el$.idx] <- list(el$.val)
    out <- ordered
    # out2<-parallel::parLapply(cl,X,tparfunc,...)
  } else out <- lapply(X, eval(parse(text=fn),envir =parent.frame()),...)
  return(out)
}

#' Adjust standata from ctsem to only use specific subjects
#'
#' @param standata standata
#' @param subjects vector of subjects
#' @param timestep ignored at present
#'
#' @return list of updated structure
#' @export
#'
#' @examples
#' d <- standatact_specificsubjects(ctstantestfit$standata, 1:2)
standatact_specificsubjects <- function(standata, subjects,timestep=NA){
  long <- standatatolong(standata)
  long <- long[long$subject %in% subjects,]
  standatamerged <- standatalongremerge(long=long, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(long))
  if(standata$ntipred > 0) standatamerged$tipredsdata <- standatamerged$tipredsdata[unique(standatamerged$subject),,drop=FALSE]
  standatamerged$nsubjects <- as.integer(length(unique(standatamerged$subject)))
  standatamerged$subject <- array(as.integer(factor(standatamerged$subject)))
  standatamerged$idmap <- standata$idmap[standata$idmap$new %in% subjects,]
  return(standatamerged)
}  

standatalongobjects <- function() {
  longobjects <- c('subject','time','dokalmanrows','nobs_y','ncont_y','nbinary_y',#'nordinal_y','whichordinal_y',
    'Y','tdpreds', 'whichobs_y','whichbinary_y','whichcont_y')
  return(longobjects)
}

standatatolong <- function(standata, origstructure=FALSE,ctm=NA){
  long <- lapply(standatalongobjects(),function(x) as.matrix(standata[[x]]))
  names(long) <- standatalongobjects()
  
  if(origstructure){
    if(is.na(ctm[1])) stop('Missing ctm arg in standatatolong()')
    colnames(long[['Y']]) <- ctm$manifestNames#colnames(standata$Y)
    long[['Y']][long[['Y']] %in% 99999] <- NA
    colnames(long[['subject']]) <- ctm$subjectIDname
    colnames(long[['time']]) <- ctm$timeName
    longout <- data.frame(long[['subject']],long[['time']],long[['Y']])
    if(standata$ntdpred > 0){
      colnames(long[['tdpreds']]) <- colnames(standata$tdpreds)
      if(!is.na(ctm[1])) colnames(long[['tdpreds']]) <- ctm$TDpredNames
      longout <- cbind(longout,long[['tdpreds']])
    }
    if(standata$ntipred > 0){
      tipreds <- standata$tipredsdata[longout[[ctm$subjectIDname]],,drop=FALSE]
      tipreds[tipreds %in% 99999] <- NA
      if(!is.na(ctm[1])) colnames(tipreds) <- ctm$TIpredNames
      longout <- cbind(longout,tipreds)
    }
  } else longout <- data.frame(long) 
  
  #,simplify=data.frame(subject=standata$subject, time=standata$time
  # colnames(long)[colnames(long) %in% 'Y'] <- paste0('Y.1'
  # colnames(long)[colnames(long) %in% 'tdpreds'] <- paste0('tdpreds.1'
  return(longout)
}

standatalongremerge <- function(long, standata){ #merge an updated long portion of standata into original standata
  n = names(standata)
  standatamerged <- lapply(names(standata), function(x) {
    if(x %in% standatalongobjects()){
      objdims <- dim(standata[[x]])
      if(is.null(objdims)) objdims <- c()
      objdims[1] <- nrow(long)
      xdat <- unlist(long[,grep(paste0('^',x),colnames(long)),drop=FALSE])
      if(is.null(xdat)) xdat <- NA
      return(array(xdat, dim = objdims))
    } else return(standata[[x]])
  })
  names(standatamerged) <- n
  return(standatamerged)
}

standataFillTime <- function(standata, times, subject, maintainT0=FALSE){
  long <- standatatolong(standata)
  
  if(any(!times %in% long$time)){ #if missing any times, add empty rows
    nlong <- do.call(rbind,
      lapply(subject, function(si){
        mintime <- min(long$time[long$subject==si])
        originaltimes <- round(long$time[long$subject==si],10)
        stimes <- times[(!times %in% originaltimes)]
        if(maintainT0) stimes <-stimes[stimes > mintime]
        data.frame(subject=si,time=stimes)
      })
    )
    
    nlong <- suppressWarnings(data.frame(nlong,long[1,!colnames(long) %in% c('subject','time')]))
    nlong[,grep('(^nobs)|(^which)|(^ncont)|(^nbin)',colnames(nlong))] <- 0L
    nlong[,grep('^dokalman',colnames(nlong))] <- 1L
    nlong[,grep('^Y',colnames(nlong))] <- 99999
    nlong[,grep('^tdpreds',colnames(nlong))] <- 0
    
    long <- rbind(long,nlong)
  } #end empty rows addition
  
  long <- long[order(long$subject,long$time),]
  standatamerged <- standatalongremerge(long=long, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(long))
  return(standatamerged)
}

stan_constrainsamples<-function(sm,standata, samples,cores=2, cl=NA,
  savescores=FALSE,
  savesubjectmatrices=TRUE,
  dokalman=TRUE,
  onlyfirstrow=FALSE, #ifelse(any(savesubjectmatrices,savescores),FALSE,TRUE),
  pcovn=2000,
  quiet=FALSE){
  if((savesubjectmatrices || savescores) && !dokalman){
    dokalman <- TRUE
    warning('savescores or savesubjectmatrices = TRUE requires dokalman=TRUE also!')
  }
  standata$savescores <- as.integer(savescores)
  standata$dokalman <- as.integer(dokalman)
  standata$savesubjectmatrices<-as.integer(savesubjectmatrices)
  if(onlyfirstrow) standata$dokalmanrows <- as.integer(c(1,diff(standata$subject)))
  
  if(!quiet) message('Computing quantities for ', nrow(samples),' samples...')
  if(nrow(samples)==1) cores <- 1
  
  if(cores > 1) {
    if(all(is.na(cl))){
      cl <- makeClusterID(cores)
      on.exit(try(parallel::stopCluster(cl),silent=TRUE),add = TRUE)
    }
    clusterIDexport(cl, c('sm','standata','samples'))
    clusterIDeval(cl,list(
      'require(data.table)',
      'smf <- ctsem::stan_reinitsf(sm,standata)',
      'tparfunc <- function(x){ 
         out <- try(data.table::as.data.table(lapply(1:length(x),function(li){
        unlist(rstan::constrain_pars(smf, upars=samples[x[li],]))
      })))
      if(!"try-error" %in% class(out)) return(out)
  }'))
    
  }
  
  if(cores ==1){
    smf <- stan_reinitsf(sm,standata) 
    tparfunc <- function(x){ 
      out <- try(data.table::as.data.table(lapply(1:length(x),function(li){
        unlist(rstan::constrain_pars(smf, upars=samples[x[li],]))
      })))
      if(!'try-error' %in% class(out)) return(out)
    }
  }
  transformedpars <- try(flexlapplytext(cl, 
    1:nrow(samples),
    'tparfunc',cores=cores))
  nulls <- unlist(lapply(transformedpars,is.null))
  if(any(nulls==FALSE)) transformedpars <- transformedpars[!nulls] else stop('No admissable samples!?')
  if(sum(nulls)>0) message(paste0(sum(nulls)/length(nulls)*100,'% of samples inadmissable'))
  
  if(cores >1) smf <- stan_reinitsf(sm,standata) #needs to be after, weird parallel stuff...
  skel= rstan::constrain_pars(smf, upars=samples[which(!nulls)[1],,drop=FALSE]) 
  transformedpars <- t(data.table::as.data.table(transformedpars))
  
  nasampscount <- nrow(transformedpars)-nrow(samples) 
  
  
  if(nasampscount > 0) {
    message(paste0(nasampscount,' NAs generated during final sampling of ', nrow(samples), '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
  }
  if(nasampscount < nrow(samples)){ 
    nresamples <- nrow(samples) - nasampscount
  } else{
    message('All samples contain NAs -- returning anyway')
    nresamples <- nrow(samples) 
  }
  transformedpars=tostanarray(flesh=transformedpars, skeleton = skel)
  
  return(transformedpars)
}

tostanarray <- function(flesh, skeleton){
  skelnames <- names(skeleton)
  skelstruc <- lapply(skeleton,dim)
  count=1
  npars <- ncol(flesh)
  niter=nrow(flesh)
  out <- list()
  for(ni in skelnames){
    if(prod(skelstruc[[ni]])>0){
      if(!is.null(skelstruc[[ni]])){
        out[[ni]] <- array(flesh[,count:(count+prod(skelstruc[[ni]])-1)],dim = c(niter,skelstruc[[ni]]))
        count <- count + prod(skelstruc[[ni]])
      } else {
        out[[ni]] <- array(flesh[,count],dim = c(niter))
        count <- count + 1
      }
    }
  }
  return(out)
}

# Workers, quietly.
#
# Two sources of console noise, both fixed here rather than downstream.
#
# `outfile = ""` forwards every worker's stdout and stderr to the console, so a
# two-core fit announced `starting worker pid=...` twice per optimisation pass
# and nothing else ever used the channel. `ctsem` in `default_packages` loaded
# it before any code could run, so its dependencies' startup warnings -- most
# visibly `package 'Rcpp' was built under R version ...` -- arrived through that
# same channel, once per worker per pass. Measured on a two-core fit: six
# copies of one warning.
#
# Loading ctsem through `clusterEvalQ` instead puts it inside something that can
# be silenced, and the workers have it before any work is dispatched either way.
# `options(ctsem.cluster.outfile = "")` restores the old behaviour for anyone
# debugging a worker, which is the only thing it was useful for.
makeClusterID <- function(cores = parallel::detectCores()) {
  outfile <- getOption("ctsem.cluster.outfile", NULL)
  arguments <- list(cores, useXDR = FALSE,
    # Workers otherwise search their own default .libPaths(), not the
    # caller's -- so `library(ctsem)` two lines down can silently load a
    # different install than the one running this code (e.g. a stale
    # globally-installed package while a development tree is under test).
    # Passing the caller's own search path down closes that gap; see
    # parallelly::makeClusterPSOCK's rscript_libs documentation.
    rscript_libs = .libPaths(),
    default_packages = c("datasets", "utils", "grDevices", "graphics",
      "stats", "methods"))
  if (!is.null(outfile)) arguments$outfile <- outfile
  cl <- do.call(parallelly::makeClusterPSOCK, arguments)
  invisible(parallel::clusterEvalQ(cl,
    suppressWarnings(suppressPackageStartupMessages(library(ctsem)))))
  duplicateNodeIDs <- TRUE
  while(duplicateNodeIDs){ 
    nodeids=unlist(parallel::clusterEvalQ(cl,{
      assign('nodeid',runif(1,0,99999999))#,envir = globalenv())
    }))
    duplicateNodeIDs <- any(duplicated(nodeids))
  }
  nodeids <- cbind(nodeids,order(nodeids))
  parallel::clusterExport(cl, varlist = "nodeids",envir   = environment())
  nodeids=unlist(parallel::clusterEvalQ(cl,{
    assign('nodeid',nodeids[nodeids[,1] %in% nodeid,2])#,envir = globalenv())
  }))
  return(invisible(cl))
}

clusterIDexport <- function(cl, vars){
  parallel::clusterExport(cl,vars,envir = parent.frame())
}

clusterIDeval <- function(cl,commands){
  clusterIDexport(cl,'commands')
  unlist(parallel::clusterEvalQ(cl = cl, 
    lapply(commands,function(x){
      eval(parse(text=x),envir = globalenv())
    })),
    recursive = FALSE)
}

# `mizemaxiter`, `grad_tol`, `step_tol` and `lbfgs_memory` are formals rather
# than part of `...` deliberately: `...` goes to sgd(), which takes none of them,
# and `maxiter` there is sgd's own and a different number from mize's.
ctOptim <- function(init, lpgFunc, tol, nsubsets, stochastic, stochasticTolAdjust,
  bfgsType='mize', mizemaxiter=99999L, grad_tol=0, step_tol=0, lbfgs_memory=100L, ...){
  if(nsubsets > 1) stochastic <- TRUE #if nsubsets > 1, use stochastic
  if(stochastic){
    args <- list(...)
    args$itertol=tol* stochasticTolAdjust
    args$lpgFunc <- lpgFunc
    args$init <- init
    args$nsubsets = nsubsets
    f <- try(do.call(sgd,args))
  }
  if(!stochastic || 'try-error' %in% class(f)){
    if(bfgsType == 'optim'){ #consider using this when importance sampling, or maybe use it to attain hessian?
      opt <- optim(par=init, fn=function(x) -lpgFunc(x)[1], gr=function(x) -attributes(lpgFunc(x))$gradient,
        method='BFGS', hessian=TRUE, control=list(maxit=99999, reltol=tol, trace=0))
      Hinv <- tryCatch(solve(opt$hessian), error=function(e) MASS::ginv(opt$hessian))
      f <- list(par=opt$par, value=-opt$value, hessian=opt$hessian, hessian_inv=Hinv,
        convergence=opt$convergence, message=opt$message, counts=opt$counts)
    } 
    if(bfgsType == 'mize') {
      mizelpg=list(  # create log prob and gradient list of functions needed for mize optim
        fg=function(pars){
          r=-lpgFunc(pars)
          r=list(fn=r[1],gr= -attributes(r)$gradient)
          return(r)
        },
        fn=function(x) -lpgFunc(x),
        gr=function(pars) -attributes(lpgFunc(pars))$gradient
      )
      f=mize::mize(init, fg=mizelpg, max_iter=mizemaxiter,
        method="L-BFGS",memory=lbfgs_memory,
        line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
        abs_tol=tol,grad_tol=grad_tol,rel_tol=0,step_tol=step_tol,ginf_tol=0)
      f$value = -f$f #reverse because mize minimizes
    }
  } #end bfgs section
  return(f)
}

carefulfitFunc <- function(cl, standata, sm, optimcores, subsamplesize, nsubsets,optimArgs,notipredsfirstpass,smfile){
  
  message('1st pass optimization (carefulfit)...')
  if(subsamplesize < 1){
    smlnsub <- min(standata$nsubjects,max(min(30,optimcores*2),ceiling(standata$nsubjects * subsamplesize)))
    standatasml <- standatact_specificsubjects(standata,
      sample(unique(standata$subject),smlnsub))
  } else standatasml <- standata
  standatasml$priors <- 1L
  standatasml$nsubsets <- as.integer(nsubsets)
  
  if(standatasml$ntipredeffects > 0 && notipredsfirstpass){
    TIPREDEFFECTsetup <- standatasml$TIPREDEFFECTsetup
    standatasml$TIPREDEFFECTsetup[,] <- 0L
    ntipredeffects <- standatasml$ntipredeffects
    standatasml$ntipredeffects <- 0L
    ninit <- length(optimArgs$init)-max(TIPREDEFFECTsetup)
    optimArgs$init <- optimArgs$init[1:ninit] #remove tipred inits
  }
  
  if(optimcores > 1) parallelStanSetup(cl = cl,standata = standatasml,split=TRUE,nsubsets = nsubsets,smfile=smfile)
  if(optimcores==1) optimArgs$lpgFunc <- singlecoreStanSetup(standata = standatasml, nsubsets = nsubsets,sm=sm)
  optimArgs$tol <- optimArgs$tol * 1000 
  optimArgs$maxiter <- 50
  optimArgs$nsubsets= nsubsets
  optimArgs$worsecountconverge <- 20
  fit = do.call(ctOptim,optimArgs)
  
  if(standata$ntipredeffects > 0 && notipredsfirstpass && !standata$TIpredAuto){
    message('Including tipred effects...')
    standata$TIPREDEFFECTsetup <- TIPREDEFFECTsetup
    standata$ntipredeffects <- ntipredeffects
    optimArgs$init <- c(fit$par,rep(0,max(TIPREDEFFECTsetup)))
    if(optimcores > 1) parallelStanSetup(cl = cl,standata = standata,split=TRUE,nsubsets = nsubsets,smfile=smfile)
    if(optimcores==1) optimArgs$lpgFunc <- singlecoreStanSetup(standata = standata, nsubsets = nsubsets,sm=sm)
    fit = do.call(ctOptim,optimArgs)
  }
  
  return(fit)
}

autoTIpredsFunc <- function(cl, standata, sm, optimArgs, parsteps, optimcores, cores) {
  initbase <- optimArgs$init
  optimArgs$maxiter=500
  optimArgs$worsecountconverge=20
  tifinished <- FALSE
  found <- 0
  nbasepars <-length(optimArgs$init)-standata$ntipredeffects
  optimArgsReduced <- optimArgs
  optimArgs$init <- optimArgs$init#[1:(npars-standata$ntipredeffects)] #remove tipred inits
  standatabase <- standata
  while (!tifinished) {
    message('Looking for tipred effects...')
    oldtia <- standata$TIPREDEFFECTsetup
    fit <- list(stanfit = list(rawest = optimArgs$init), #use full size init vec with updated inits as new tipreds included
      standata = standatabase,
      stanmodel = sm)
    tia <- ctTIauto(fit,cores=cores) #problem here on second time around, sm seems to be updated
    tia[tia > .05] <- 0
    tia[tia > 0] <- seq_along(tia[tia > 0]) #assign sequential numbers to new predictors
    
    if (max(tia) > found) { #if new predictors found
      standata$TIPREDEFFECTsetup <- array(as.integer(tia), dim = dim(tia)) #coerce tia into TIPREDEFFECTsetup
      found <- max(tia) #update number found
      message('Found ', found, ' viable TIpred effects')
      standata$ntipredeffects <- as.integer(found)
      optimArgsReduced$init <- head(optimArgs$init, nbasepars+found) #remove unused ti inits from init
      # optimArgs$init <- head(optimArgs$init, length(initbase))
      if((found + nbasepars) == length(optimArgs$init)) tifinished <- TRUE
    } else {
      tifinished <- TRUE
      message('No further predictors found, finishing optimization...')
    }
    if (optimcores > 1) parallelStanSetup(cl = cl, standata = standata, split = TRUE, nsubsets = 1)
    if (optimcores == 1) optimArgsReduced$lpgFunc <- singlecoreStanSetup(standata = standata, nsubsets = 1,sm=sm)
    iter <- 0L
    optimfit <- do.call(ctOptim, optimArgsReduced)
    optimArgs$init[1:(nbasepars+found)] <- optimfit$par #update full init vec
  }
  
  return(list(standata   = standata,optimfit   = optimfit))
} # end ti pred auto function

imis_is <- function(parlp,
  mu_hat,
  Sigma_hat,
  cl,
  n_batch       = 1000,
  target_ess    = 100,
  max_iter      = 10,
  scale_init    = 1.5,
  tail_scale    = 1.2,
  df            = Inf,
  ridge         = 1e-8,
  finishsamples = 1000,
  verbose       = TRUE,
  diag_plots    = TRUE) {
  
  for (pkg in c("mvtnorm", "diagis", "gridExtra", "ggplot2", "grid"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Install '%s' first.", pkg))
  
  if (verbose)
    message(sprintf(
      "Importance sampling: target ESS = %d, max_iter = %d, batch = %d",
      target_ess, max_iter, n_batch))
  
  ## ── helpers ──────────────────────────────────────────────────────────
  safe_pd <- function(S, eps = ridge) {
    S2 <- S + diag(eps, nrow(S))
    while (any(eigen(S2, TRUE, TRUE)$values <= 0))
      S2 <- S2 + diag(eps, nrow(S))
    S2
  }
  Sigma_hat <- safe_pd(Sigma_hat)
  
  logplus <- function(a, b) {
    idx <- a > b
    r <- numeric(length(a))
    r[idx]  <- a[idx] + log1p(exp(b[idx] - a[idx]))
    r[!idx] <- b[!idx] + log1p(exp(a[!idx] - b[!idx]))
    r
  }
  ess   <- function(w) diagis::ess(w)
  rsamp <- function(w, N) {
    cs <- cumsum(w / sum(w))
    findInterval((runif(1) + 0:(N - 1)) / N, cs) + 1L
  }
  
  ## ── containers ───────────────────────────────────────────────────────
  comp_mu  <- list(mu_hat)
  # `scale_init^2 * Sigma`, the whole matrix. This was written
  # `Sigma_hat * (diag(scale_init^2-1, n) + 1)`, an elementwise product with a
  # matrix carrying `scale_init^2` on the diagonal and *1* off it -- so it
  # inflated the variances, left the covariances untouched, and thereby divided
  # every proposal correlation by `scale_init^2`. That is not a wider proposal
  # but a differently shaped one, and along the correlated directions it is
  # narrower than `Sigma` itself, which is the opposite of what a scale above
  # one is for. Computed exactly for a Gaussian target in the nine identified
  # dimensions of a 400-subject fit, ESS/n at `scale_init = 1.5` was 0.058 the
  # old way against 0.190 this way; run end to end on that fit at the julia
  # defaults and a fixed seed, the old form spent all 51,000 evaluations to
  # reach an effective sample of 7.7 and this one reached 144 in 4,000.
  comp_cov <- list(Sigma_hat * scale_init^2)
  T_comp   <- 1L
  
  samples   <- matrix(0, 0, length(mu_hat))
  log_p_all <- numeric(0)
  log_qsum  <- numeric(0)
  w_raw     <- numeric(0)
  ess_now   <- 0
  interrupted <- FALSE
  
  # A multivariate t rather than a normal, when `df` is finite.
  #
  # Available, and not the default, because it was measured and did not help.
  # The theory says an importance-sampling proposal wants heavier tails than its
  # target; on a 40-subject ctsem model the t was consistently *worse* than the
  # normal at the same scale -- mean standard error 0.70 of a reference sample's
  # against the normal's 0.75 -- and matched it only at a scale wide enough to
  # drop the effective sample size from 411 to 76, which is not a trade worth
  # making. Kept as a knob because that is one model.
  #
  # `Sigma` is the t's *scale* matrix here, not rescaled so its covariance
  # equals `Sigma`. That rescaling by `(df-2)/df` is the obvious-looking move
  # and it is wrong: it shrinks the bulk to pay for the heavy tails, which is
  # the opposite of the point. The first version did it and came out narrower
  # than the normal it was meant to widen on.
  #
  # Importance sampling needs a proposal with *heavier* tails than the target,
  # because a region the proposal never visits cannot be upweighted however
  # large its weight would have been. A normal proposal fitted to the Laplace
  # curvature has lighter tails than the posterior it is approximating, which is
  # exactly backwards, and the failure is quiet: the effective sample size looks
  # healthy because the draws that exist agree with each other, while the
  # answer stays close to the proposal. Measured on a 40-subject model, a normal
  # proposal at scale 1.1 returned standard errors within 10% of the Hessian's
  # where the true posterior was up to twice as wide.
  rprop <- function(n, mu, Sigma) {
    if (!is.finite(df)) return(mvtnorm::rmvnorm(n, mu, Sigma))
    mvtnorm::rmvt(n, sigma = Sigma, df = df, delta = mu, type = "shifted")
  }
  dprop <- function(x, mu, Sigma) {
    if (!is.finite(df)) return(mvtnorm::dmvnorm(x, mu, Sigma, log = TRUE))
    mvtnorm::dmvt(x, delta = mu, sigma = Sigma, df = df, log = TRUE,
      type = "shifted")
  }
  draw_mix <- function(n) {
    if (T_comp == 1L)
      rprop(n, comp_mu[[1]], comp_cov[[1]])
    else {
      sel <- sample.int(T_comp, n, TRUE)
      do.call(rbind, lapply(seq_len(T_comp), function(k) {
        m <- sum(sel == k)
        if (m) rprop(m, comp_mu[[k]], comp_cov[[k]])
      }))
    }
  }
  
  
  
  
  ## ── main loop ─────────────────────────────────────────────────────────
  for (it in 0:max_iter) {
    
    x_new <- draw_mix(n_batch)
    
    ## ---------- log-p with interrupt guard -----------------------------
    log_p_new <- {
      if (!is.null(cl) && length(cl) > 1) {
        parallel::clusterExport(cl, "x_new", envir = environment())
        unlist(parallel::parLapply(
          cl, seq_len(nrow(x_new)), \(i) parlp(x_new[i, ])), FALSE)
      } else {
        vapply(seq_len(nrow(x_new)), \(i) parlp(x_new[i, ]), numeric(1))
      }
    }
    
    
    ## ---------- mixture log-q -----------------------------------------
    lq_new <- rep.int(-Inf, n_batch)
    for (k in seq_len(T_comp))
      lq_new <- logplus(lq_new, dprop(x_new, comp_mu[[k]], comp_cov[[k]]))
    
    samples   <- rbind(samples, x_new)
    log_p_all <- c(log_p_all, log_p_new)
    log_qsum  <- c(log_qsum,  lq_new)
    
    log_mix <- log_qsum - log(T_comp)
    log_w   <- log_p_all - log_mix
    log_w   <- log_w - max(log_w)
    w_raw   <- exp(log_w)
    ess_now <- ess(w_raw)
    
    if (verbose)
      message(sprintf("\r iter %2d | n %6d | ESS %8.1f",
        it + 1, nrow(samples), ess_now),
        appendLF = FALSE)
    
    if (interactive() && diag_plots) {
      g <- diagis::weight_plot(w_raw)
      gridExtra::grid.arrange(
        g, top = grid::textGrob(
          sprintf("IMIS diagnostics - ESS %.1f / %d",
            ess_now, nrow(samples)),
          gp = grid::gpar(fontface = "bold", fontsize = 14)))
    }
    
    if (ess_now >= target_ess || it == max_iter) break
    
    ## ---------- add new component -------------------------------------
    top_idx <- w_raw > quantile(w_raw, 0.9)
    comp_mu[[T_comp + 1L]] <- diagis::weighted_mean(
      samples[top_idx,,drop=FALSE], w_raw[top_idx])
    # `tail_scale^2 *` the weighted covariance, for the same reason the initial
    # component is scaled that way above: the elementwise form this replaces
    # left the covariances at their unscaled values and so shrank the
    # correlations of every added component.
    comp_cov[[T_comp + 1L]] <- safe_pd(
      diagis::weighted_var(
        samples[top_idx,,drop=FALSE], w_raw[top_idx]) * tail_scale^2)
    T_comp <- T_comp + 1L
    
    lq_newcomp <- mvtnorm::dmvnorm(samples, comp_mu[[T_comp]],
      comp_cov[[T_comp]], log = TRUE)
    log_qsum <- logplus(log_qsum, lq_newcomp)
  } #end main loop
  
  if (verbose) message("")   # newline after progress line
  
  w_norm <- if (length(w_raw)) w_raw / sum(w_raw) else numeric(0)
  idx_eq <- if (length(w_norm)) rsamp(w_norm, finishsamples) else integer(0)
  
  list(theta        = if (length(idx_eq)) samples[idx_eq,,drop=FALSE] else samples,
    lpsamples    = if (length(idx_eq)) log_p_all[idx_eq] else log_p_all,
    weights      = if (length(idx_eq)) rep(1/finishsamples, length(idx_eq)) else numeric(0),
    full_theta   = samples,
    full_weights = w_norm,
    ess          = ess_now,
    mean         = if (length(w_norm))
      as.numeric(diagis::weighted_mean(samples, w_norm))
    else rep(NA_real_, length(mu_hat)),
    covariance   = if (length(w_norm))
      diagis::weighted_var(samples, w_norm)
    else matrix(NA_real_, length(mu_hat), length(mu_hat)),
    df_used      = df)
  
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Optimize / importance sample a stan or ctStan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init vector of unconstrained parameter values, or character string 'random' to initialise with
#' random values very close to zero.
#' @param initsd positive numeric specifying sd of normal distribution governing random sample of init parameters,
#' if init='random' .
#' @param stochastic Logical. Use stochastic gradient descent as main optimizer. Always finishes (double checks) with mize (bfgs) optimizer.
#' @param plot Logical. If TRUE, plot iteration details. Probably slower.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param tol objective tolerance -- the optimizer stops when the objective
#' changes by less than this between evaluations.
#' @param g_tol gradient tolerance: stop when the l2 norm of the gradient falls
#' below this. \code{NULL} (the default) leaves it off, which is what this
#' optimizer has always done -- \code{tol} is its criterion. The julia backend
#' reads the same name and defaults it to 1e-8 instead.
#' @param x_tol parameter-step tolerance: stop when the update to the parameter
#' vector is smaller than this. \code{NULL} leaves it off.
#' @param maxiter maximum optimizer iterations. \code{NULL} leaves the existing
#' caps in place (99999 for the bfgs optimizer, 5000 for the stochastic one).
#' @param lbfgs_memory number of curvature pairs L-BFGS keeps. \code{NULL} uses
#' 100.
#' @param priors logical. If TRUE, a priors integer is set to 1 (TRUE) in the standata object -- only has an effect if 
#' the stan model uses this value. 
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when priors=FALSE
#' @param stallretries Integer. Number of times to restart the optimizer from fresh random
#' values when it finishes somewhere that is not a maximum -- see \code{stalltol}. Only applies
#' when \code{init='random'}; with supplied inits a failed fit is reported but not retried.
#' Set to 0 to disable retrying.
#' @param stalltol Gradient per data point above which the optimizer is taken to have stopped
#' short rather than converged. Rough likelihoods -- binary indicators are the case this was
#' built for -- can leave both optimizers unable to find an improving step, which they report as
#' convergence, returning the starting values as the estimate. Converged fits sit near 1e-5 per
#' data point, stalled ones at 1 or more.
#' @param subsamplesize value between 0 and 1 representing proportion of subjects to include in first pass fit. 
#' @param cores Number of cpu cores to use, should be at least 2.
#' @param uncertainty Character string selecting the optimized-fit uncertainty
#' approximation. Options are \code{'hessian'}, \code{'surrogate'},
#' \code{'is'}, \code{'bootstrap'}, \code{'fullbootstrap'}, \code{'sandwich'},
#' and \code{'opg'}.
#' @param uncertaintyDraws Character string controlling approximate
#' raw-parameter draws
#' from the approximate uncertainty. \code{'auto'} uses empirical draws for
#' \code{uncertainty='bootstrap'} or \code{uncertainty='fullbootstrap'} and
#' normal draws otherwise. \code{'normal'} draws from a multivariate normal
#' using the selected covariance, \code{'empirical'} uses empirical draws when
#' available, and \code{'imis'} runs importance sampling.
#' @param uncertaintyControl List of method-specific options passed to
#' \code{\link{ctOptimUncertainty}} internals. Score-based methods use
#' subject-level score contributions when there are at least two subjects;
#' single-subject models warn and use case-level contributions. Score-based
#' methods warn when there are fewer than ten independent subjects or no more
#' score rows than raw parameters. Full bootstrap requires at least two
#' subjects and warns below ten independent subjects. Bootstrap-style methods
#' require at least two returned samples / refits.
#' @param finishsamples Number of samples to draw (either from hessian
#' based covariance or posterior distribution) for final results computation.
#' @param parsteps ordered list of vectors of integers denoting which parameters should begin fixed
#' at zero, and freed sequentially (by list order). Useful for complex models, e.g. keep all cross couplings fixed to zero 
#' as a first step, free them in second step.
#' @param matsetup subobject of ctStanFit output. If provided, parameter names instead of numbers are output for any problem indications.
#' @param nsubsets number of subsets for stochastic optimizer. Subsets are further split across cores, 
#' but each subjects data remains whole -- processed by one core in one subset.
#' @param stochasticTolAdjust Multiplier for stochastic optimizer tolerance. 
#' @param lproughnesstarget target log posterior roughness for stochastic optimizer (suggest between .05 and .4).
#' @return list containing fit elements
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @importFrom parallelly makeClusterPSOCK

stanoptimis <- function(standata, sm, init='random',initsd=.01,
  estonly=FALSE,tol=1e-8,
  stochastic = TRUE,
  priors=TRUE,carefulfit=TRUE,
  stallretries=2,stalltol=1e-2,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg'),
  uncertaintyDraws='auto',
  uncertaintyControl=list(),
  subsamplesize=1,
  parsteps=c(),
  plot=FALSE,
  finishsamples=1000,
  lproughnesstarget=.2,
  verbose=0,
  cores=2,
  matsetup=NA,
  nsubsets=1,
  stochasticTolAdjust=1000,
  # Appended rather than placed beside `tol`, so that a positional call to this
  # exported function keeps meaning what it did. NULL is "leave the optimizer as
  # it was"; see the vocabulary table at the top of R/ctFit.R.
  g_tol=NULL, x_tol=NULL, maxiter=NULL, lbfgs_memory=NULL){
  
  
  
  # initial checks ----------------------------------------------------------
  
  if(!interactive()) plot <- FALSE #if not interactive, don't plot
  
  if(!is.null(standata$verbose)) {
    if(verbose > 1) standata$verbose=as.integer(verbose) else standata$verbose=0L
  }
  standata$priors=as.integer(priors)
  
  if(nsubsets > (standata$nsubjects/10)) nsubsets <- ceiling(standata$nsubjects/10) #restrict to 10 subsets per 100 subjects
  if(nsubsets > (standata$nsubjects/cores)) nsubsets <- max(1,ceiling(standata$nsubjects/cores)) #minimum 1 subset
  
  if(is.null(init)) init <- 'random' # if no inits are given, use random initialisation
  if(init[1] !='random') carefulfit <- FALSE # if inits are given, carefulfit is not needed
  
  savesubjectmatrices <- standata$savesubjectmatrices
  standata$savesubjectmatrices <- 0L #reinsert when saving samples
  
  optimArgs <- list(init=init,
    lpgFunc=NA,
    tol=tol,
    stochastic=stochastic,
    plot=plot,
    nsubsets=nsubsets,
    maxiter=5000,
    stochasticTolAdjust=stochasticTolAdjust,
    lproughnesstarget=lproughnesstarget,
    parrangetol=1e-6,
    whichignore=integer())

  # The three tolerances mize also offers, and its iteration cap and memory.
  # ctsem pinned grad_tol and step_tol at zero and max_iter at 99999, which is
  # why the gradient criterion looked like something only the julia engine had;
  # mize has had it all along. NULL means "leave it as it was", so an existing
  # call optimises exactly as before. `optimArgs$maxiter` is left alone because
  # it is sgd's, not mize's, and the two have never been the same number.
  optimArgs$mizemaxiter <- if(is.null(maxiter)) 99999L else as.integer(maxiter)[1]
  optimArgs$grad_tol <- if(is.null(g_tol)) 0 else as.numeric(g_tol)[1]
  optimArgs$step_tol <- if(is.null(x_tol)) 0 else as.numeric(x_tol)[1]
  optimArgs$lbfgs_memory <- if(is.null(lbfgs_memory)) 100L else as.integer(lbfgs_memory)[1]
  if(!is.null(maxiter)) optimArgs$maxiter <- as.integer(maxiter)[1]
  
  smf <- stan_reinitsf(sm,standata)
  npars=rstan::get_num_upars(smf)
  
  if(stochastic=='auto' && npars > 50){
    message('> 50 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
    stochastic <- TRUE
  } else if(stochastic=='auto') stochastic <- FALSE
  if(length(parsteps)>0 && !stochastic){
    stochastic=TRUE
    message('Stochastic optimizer used for data driven parameter inclusion') 
  }
  
  optimcores <- ifelse(length(unique(standata$subject)) < cores, length(unique(standata$subject)),cores)
  if(optimcores > 1) rm(smf)
  
  if(plot > 0 && .Platform$OS.type=="windows" && interactive()) {
    dev.new(noRStudioGD = TRUE)
    on.exit(expr = {try({dev.off()})},add = TRUE)
  }
  
  message('Using ',cores,'/', parallel::detectCores(),' CPU cores')
  if(cores > parallel::detectCores()) warning('More cores requested than available, if errors occur, try reducing cores')
  
  storedPars <- as.numeric(c())
  storedLp <- c()
  
  optimfinished <- FALSE
  on.exit({
    if(!optimfinished){
      message('Optimization cancelled -- restart from current point by including this argument:')
      message((paste0(c('inits = c(',   paste0(round(storedPars,5),collapse=', '), ')'    ))))
    }},add=TRUE)
  
  
  # initial values ----------------------------------------------------------
  
  randominit <- all(optimArgs$init %in% 'random')
  if(randominit){
    optimArgs$init <- rnorm(npars, 0, initsd)
    if(length(parsteps)>0) optimArgs$init[unlist(parsteps)] <- 0 
  }
  
  if(all(optimArgs$init %in% 0)) optimArgs$init <- rep(0,npars)
  
  if(length(optimArgs$init) != npars){
    warning('Initialisation vector length does not match number of parameters in model, extending with zeros')
    optimArgs$init=c(optimArgs$init[1:min(length(optimArgs$init),npars)],rep(0,abs(npars-length(optimArgs$init))))
  }
  
  if(any(is.na(optimArgs$init))){
    warning('Initialisation vector contains NAs, replacing with zeros')
    optimArgs$init[is.na(optimArgs$init)] <- 0
  }
  
  # initialise cluster ------------------------------------------------------
  if(cores > 1){
    if(standata$recompile > 0){
      smfile <- file.path(tempdir(),paste0('ctsem_sm_',ceiling(runif(1,0,100000)),'.rda'))
      save(sm,file=smfile,eval.promises = FALSE,precheck = FALSE)
      on.exit(add = TRUE,expr = {file.remove(smfile)})
    } else smfile <- ''
  } #end smfile setup
  
  if(optimcores > 1){ #for parallelised computation
    clctsem=makeClusterID(optimcores)
    on.exit(try({parallel::stopCluster(clctsem)},silent=TRUE),add=TRUE)
  }
  
  ######log prob function setup#######
  
  lpg_single<-function(parm) { #single core log prob function, used for importance sampling and single core optimization
    a=Sys.time()
    out<- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
    b=Sys.time()
    if('try-error' %in% class(out) || is.nan(out)) {
      out=-1e100
      attributes(out) <- list(gradient=rep(0,length(parm)))
    }
    storedPars <<- parm
    evaltime <- b-a
    if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(evaltime,2)),digits=14)
    return(out)
  }
  
  lpg_parallel<-function(parm) {
    a=Sys.time()
    clusterIDexport(clctsem,'parm')
    out2<-  parallel::clusterEvalQ(cl = clctsem,parlp(parm))
    error <- FALSE
    tmp<-sapply(1:length(out2),function(x) {
      if(!is.null(attributes(out2[[x]])$err)){
        if(!error & length(out2) > 1 && as.logical(verbose)){
          message('Error on core ', x,' but continuing:')
          error <<- TRUE
          message(attributes(out2[[x]])$err)
        }
      }
    })
    out <- try(sum(unlist(out2)),silent=TRUE)
    coretimes <- sapply(out2,function(x) round(attributes(x)$time,3))
    for(i in seq_along(out2)){
      if(i==1) attributes(out)$gradient <- attributes(out2[[1]])$gradient
      if(i>1) attributes(out)$gradient <- attributes(out)$gradient+attributes(out2[[i]])$gradient
    }
    b=Sys.time()
    if('try-error' %in% class(out) || is.nan(out)) {
      out=-1e100
      attributes(out) <- list(gradient=rep(0,length(parm)))
    }
    if(plot > 0 && interactive() && ( (!stochastic &&!carefulfit && nsubsets ==1))){
      if(out[1] > (-1e99)) storedLp <<- c(storedLp,out[1])
      iter <<- iter+1
      g=log(abs(attributes(out)$gradient))*sign(attributes(out)$gradient)
      if(iter %% plot == 0){
        par(mfrow=c(1,3))
        plot(parm,xlab='param',ylab='par value',col=1:length(parm))
        plot(log(1+tail(-storedLp,500)-min(tail(-storedLp,500))),ylab='target',type='l')
        plot(g,type='p',col=1:length(parm),ylab='gradient',xlab='param')
      }
      if(verbose==0) message(paste('\rlp= ',out,' ,    iter time = ',round(b-a,3), '; core times = ',
        paste0(coretimes,collapse=', ')),appendLF = FALSE) #if not verbose, print lp when plotting
    }
    storedPars <<- parm
    if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(b-a,3), '; core times = ',
      paste0(coretimes,collapse=', '))) #if not verbose, print lp when plotting
    return(out)
  }
  
  if(optimcores > 1) optimArgs$lpgFunc <- lpg_parallel else optimArgs$lpgFunc <- lpg_single 
  iter <-0
  
  if(carefulfit) {
    iter <-0
    storedLp <- c()
    optimfit <- carefulfitFunc(cl=clctsem,standata=standata, sm=sm, optimcores=optimcores, 
      nsubsets=nsubsets,subsamplesize=subsamplesize,optimArgs=optimArgs,notipredsfirstpass=TRUE,smfile=smfile)
    optimArgs$init[1:length(optimfit$par)] <- optimfit$par #update non ti pred inits
  } #end carefulfit
  
  # end subsetting / carefulfit -----------------------------------------
  
  standata$nsubsets <- 1L
  optimArgs$nsubsets <- 1L
  if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=TRUE,smfile=smfile)
  if(optimcores==1) smf<-stan_reinitsf(sm,standata)
  
  
  # tipredauto --------------------------------------------------------------
  
  
  if(standata$ntipred > 0 && !is.null(standata$TIpredAuto) && standata$TIpredAuto){
    if((length(parsteps) > 0)) stop('parsteps not supported with TIpredAuto')
    
    # insert TI-pred logic via our new function
    ti_res <- autoTIpredsFunc(
      cl           = clctsem,
      standata     = standata,
      sm           = sm,
      optimArgs    = optimArgs,
      parsteps     = parsteps,
      optimcores   = optimcores,
      cores=cores
    )
    # unpack results of tipred auto
    standata  <- ti_res$standata
    optimfit  <- ti_res$optimfit
    optimArgs$init <- optimfit$par #update inits with ti pred inits
    npars <- length(optimfit$par) #update npars after ti pred auto
    if (optimcores == 1) optimArgs$lpgFunc <- singlecoreStanSetup(standata = standata, nsubsets = 1,sm=sm) #update single core lpg (parallel already updated)
  } #end ti pred auto total loop
  
  
  ##parameter stepwise / selection
  if(length(parsteps) > 0){
    message('Freeing parameters...')
    parstepsfinished <- FALSE
    while(!parstepsfinished && length(parsteps)>0){
      if(length(parsteps)>1) parsteps <- parsteps[-1] else parsteps <- c()

      optimArgs$tol <- tol * 1000 #increase tolerance for parameter freeing
      iter <-0
      optimfit <- do.call(ctOptim,optimArgs)

      if(length(parsteps)>0){
        optimArgs$init[-unlist(parsteps)] = optimfit$par
      }else{
        parstepsfinished <- TRUE
        optimArgs$init = optimfit$par
      }
    }
  }

  message('Optimizing...')

  optimArgs$nsubsets <- 1
  optimArgs$parrangetol <- tol*100
  optimArgs$whichignore <- unlist(parsteps)
  iter <-0
  optimfit <- do.call(ctOptim,optimArgs)


  if(!'try-error' %in% class(optimfit) & !'NULL' %in% class(optimfit)){
    if(length(parsteps)>0) optimArgs$init[-unlist(parsteps)] = optimfit$par else optimArgs$init=optimfit$par
  }

  #use bfgs to double check stochastic fit (or just use bfgs if requested)...
  if(stochastic){
    message('Finishing optimization...')
    optimArgs$stochastic <- FALSE
    iter <-0
    optimfit <- do.call(ctOptim,optimArgs)
  }
  optimArgs$init = optimfit$par

  # A fit can finish exactly where it started. Where the likelihood is rough --
  # binary indicators found this: the linearised measurement update lets the
  # filtered states run out to |eta| ~ 100, the logit saturates, and the log
  # likelihood then swings by hundreds over parameter changes of 1e-4 -- the
  # step acceptance in sgd rejects every proposal and permanently collapses its
  # step size, and the L-BFGS pass that follows stops on an unchanged log
  # probability. Both report success, so the starting values come back as the
  # estimate with a Hessian computed about them. The gradient tells the two
  # apart with orders of magnitude to spare: measured over converged fits it is
  # ~1e-5 per data point, and 1 to 100 per data point when the fit has stalled.
  if(length(parsteps)==0){
    stallthreshold <- max(1, stalltol * standata$ndatapoints)
    stallstate <- function(pars){
      lpg <- suppressWarnings(try(optimArgs$lpgFunc(pars),silent=TRUE))
      g <- attributes(lpg)$gradient
      if('try-error' %in% class(lpg) || is.null(g) || any(!is.finite(g)))
        return(list(lp=-Inf, maxg=Inf))
      list(lp=lpg[1], maxg=max(abs(g)))
    }
    best <- stallstate(optimArgs$init)
    bestfit <- optimfit
    attempt <- 0
    while(best$maxg > stallthreshold && randominit && attempt < stallretries){
      attempt <- attempt + 1
      message('Optimization stopped with a gradient of ',signif(best$maxg,3),
        ' -- that is not a maximum. Restarting from new values (',attempt,' of ',stallretries,')...')
      optimArgs$init <- rnorm(npars, 0, initsd)
      optimArgs$stochastic <- stochastic
      iter <- 0
      newfit <- try(do.call(ctOptim,optimArgs))
      if(!'try-error' %in% class(newfit) && !'NULL' %in% class(newfit) && stochastic){
        optimArgs$init <- newfit$par
        optimArgs$stochastic <- FALSE
        iter <- 0
        newfit <- try(do.call(ctOptim,optimArgs))
      }
      if('try-error' %in% class(newfit) || 'NULL' %in% class(newfit)) next
      new <- stallstate(newfit$par)
      if(new$maxg <= stallthreshold || new$lp > best$lp){
        best <- new
        bestfit <- newfit
      }
    }
    optimfit <- bestfit
    optimArgs$init <- optimfit$par
    if(best$maxg > stallthreshold) warning(paste0(
      'Optimization finished with a gradient of ',signif(best$maxg,3),
      ', which is not a maximum -- the estimates are wherever the optimizer stopped, ',
      'and the uncertainty is computed about that point. ',
      ifelse(randominit,
        'Restarting from new values did not help. ',
        'Supplied inits are not retried automatically -- try init="random". '),
      'Rough likelihoods do this; with binary indicators, backend="julia" integrates ',
      'the observation instead of linearising it and does not.'),immediate. = TRUE)
  }

  
  est2=optimArgs$init #because init contains the fixed values #unconstrain_pars(smf, est1)
  if(length(parsteps)>0) est2[-parsteps] = optimfit$par else est2=optimfit$par
  
  npars = length(est2)
  
  if(estonly) {
    smf <- stan_reinitsf(sm,standata)
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2,parsteps=parsteps)
    optimfinished <- TRUE #disable exit message re pars
    return(stanfit)
  }
  
  uncertainty <- match.arg(uncertainty,
    c('hessian','surrogate','is','bootstrap','fullbootstrap','sandwich','opg'))
  if(length(parsteps) > 0 && uncertainty == 'fullbootstrap') {
    stop('fullbootstrap uncertainty is not currently supported with parsteps')
  }
  uncertaintyDraws <- match.arg(uncertaintyDraws,
    c('auto','normal','empirical','imis'))
  if(uncertaintyDraws == 'auto') {
    if(uncertainty == 'is') {
      uncertaintyDraws <- 'imis'
    } else if(uncertainty %in% c('bootstrap','fullbootstrap')) {
      uncertaintyDraws <- 'empirical'
    } else {
      uncertaintyDraws <- 'normal'
    }
  }
  standata$savesubjectmatrices=savesubjectmatrices #if we save subject matrices, we need to use the full standata
  uncertaintyControl$parsteps <- parsteps
  stanfit=list(optimfit=optimfit,stanfit=stan_reinitsf(sm,standata),
    rawest=est2, rawposterior=NULL, cov=NULL,
    standata=list(TIPREDEFFECTsetup=standata$TIPREDEFFECTsetup,ntipredeffects = standata$ntipredeffects))
  fit <- list(standata=standata, stanmodel=sm,
    setup=list(matsetup=matsetup), stanfit=stanfit)
  class(fit) <- 'ctStanFit'
  fit <- ctOptimUncertainty(fit=fit, uncertainty=uncertainty,
    draws=uncertaintyDraws, finishsamples=finishsamples, cores=cores,
    control=uncertaintyControl, verbose=verbose)
  stanfit <- fit$stanfit
  
  
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}





