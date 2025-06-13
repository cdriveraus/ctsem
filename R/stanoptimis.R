
#' Sample more values from an optimized ctstanfit object
#'
#' @param fit fit object
#' @param nsamples number of samples desired
#' @param cores number of cores to use
#'
#' @return fit object with extra samples
#' @export
#'
#' @examples
#' \dontrun{
#' newfit <- ctAddSamples(ctstantestfit, 10, 1)
#' }
ctAddSamples <- function(fit,nsamples,cores=2){
  
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
    dokalman=as.logical(fit$standata$savesubjectmatrices),
    cores=cores)
  return(fit)
}


parallelStanSetup <- function(cl, standata,split=TRUE,nsubsets=1){
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
  
  parallel::clusterEvalQ(cl,{
    # g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    # environment(parlptext) <- g
    if(standata$recompile > 0) load(file=smfile) else sm <- ctsem:::stanmodels$ctsm
    # eval(parse(text=parlptext))
    # assign("parlp",parlp,pos=g)
    parlp <- function(parm){
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
    }
    if(length(stanindices[[nodeid]]) < length(unique(standata$subject))) standata <- ctsem:::standatact_specificsubjects(standata,stanindices[[nodeid]])
    standata$priormod <- 1/cores
    if(FALSE) sm=99
    smf=ctsem:::stan_reinitsf(sm,standata)
  })
  NULL
}

singlecoreStanSetup <-function(standata, nsubsets){
  cores <- 1
  standata$nsubsets <- as.integer(nsubsets)
  if(!is.null(standata$recompile)) standata$recompile <- 0 #no recompile on single core
  smf <- ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm,standata)
  return(eval(parse(text=ctsem:::parlptext)))#create parlp function)
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

flexsapply <- function(cl, X, fn,cores=1){
  if(cores > 1) parallel::parSapply(cl,X,fn) else sapply(X, fn)
}

flexlapply <- function(cl, X, fn,cores=1,...){
  if(cores > 1) parallel::parLapply(cl,X,fn,...) else lapply(X, fn,...)
}

flexlapplytext <- function(cl, X, fn,cores=1,...){
  if(cores > 1) {
    nodeindices <- split(1:length(X), sort((1:length(X))%%cores))
    nodeindices<-nodeindices[1:cores]
    clusterIDexport(cl,c('nodeindices'))
    
    out <-unlist(clusterIDeval(cl,paste0('lapply(nodeindices[[nodeid]],',fn,')')),recursive = FALSE)
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
  # colnames(long)[colnames(long) %in% 'tdpreds'] <- 'tdpreds.1'
  return(longout)
}

# stanlongtostandatasml <- function(long){
#   standatasml <- sapply( standatalongobjects(), function(x) long[,grep(paste0('^',x),colnames(long)),drop=FALSE])
#   names(standatasml) <- standatalongobjects()
#   return(standatasml)
# }

standatalongremerge <- function(long, standata){ #merge an updated long portion of standata into original standata
  n = names(standata)
  standatamerged <- lapply(names(standata), function(x) {
    if(x %in% standatalongobjects()){
      objdims <- dim(standata[[x]])
      if(is.null(objdims)) objdims <- c()
      objdims[1] <- nrow(long)
      xdat <- unlist(long[,grep(paste0('^',x),colnames(long)),drop=FALSE])
      if(is.null(xdat)) xdat <- NA
      return(array(xdat, dim=objdims))
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
  if(savesubjectmatrices && !dokalman){
    dokalman <- TRUE
    warning('savesubjectmatrices = TRUE requires dokalman=TRUE also!')
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


makeClusterID <- function(cores = parallel::detectCores()) {
  cl <- parallelly::makeClusterPSOCK(cores,
    useXDR      = FALSE,
    outfile     = "") 
  duplicateNodeIDs <- TRUE
  while(duplicateNodeIDs){ 
    nodeids=unlist(parallel::clusterEvalQ(cl,{
      assign('nodeid',runif(1,0,99999999),envir = globalenv())
    }))
    duplicateNodeIDs <- any(duplicated(nodeids))
  }
  nodeids <- cbind(nodeids,order(nodeids))
  parallel::clusterExport(cl, varlist = "nodeids",envir   = environment())
  nodeids=unlist(parallel::clusterEvalQ(cl,{
    assign('nodeid',nodeids[nodeids[,1] %in% nodeid,2],envir = globalenv())
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




ctOptim <- function(init, lpgFunc, tol, nsubsets, stochastic,stochasticTolAdjust,...){
  # browser()
  if(stochastic || nsubsets > 1){
    args <- list(...)
    args$itertol=tol* stochasticTolAdjust
    args$lpgFunc <- lpgFunc
    args$init <- init
    args$nsubsets = nsubsets
    f <- try(do.call(sgd,args))
  }
  if(!stochastic || 'try-error' %in% class(f)){
    mizelpg=list(  # create log prob and gradient list of functions needed for mize optim
      fg=function(pars){
        r=-lpgFunc(pars)
        r=list(fn=r[1],gr= -attributes(r)$gradient)
        return(r)
      },
      fn=function(x) -lpgFunc(x),
      gr=function(pars) -attributes(lpgFunc(pars))$gradient
    )
    f=mize(init, fg=mizelpg, max_iter=99999,
      method="L-BFGS",memory=100,
      line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
      abs_tol=tol,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
    f$value = -f$f #reverse because mize minimizes
  }
  return(f)
}

carefulfitFunc <- function(cl, standata, sm, optimcores, subsamplesize, nsubsets,optimArgs,notipredsfirstpass){
  
  if(standata$ntipredeffects > 0 && notipredsfirstpass){
    TIPREDEFFECTsetup <- standata$TIPREDEFFECTsetup
    standata$TIPREDEFFECTsetup[,] <- 0L
    ntipredeffects <- standata$ntipredeffects
    standata$ntipredeffects <- 0L
    ninit <- length(optimArgs$init)-max(TIPREDEFFECTsetup)
    optimArgs$init <- optimArgs$init[1:ninit] #remove tipred inits
  }
  
  message('1st pass optimization (carefulfit)...')
  if(subsamplesize < 1){
    smlnsub <- min(standata$nsubjects,max(min(30,optimcores*2),ceiling(standata$nsubjects * subsamplesize)))
    standatasml <- standatact_specificsubjects(standata,
      sample(unique(standata$subject),smlnsub))
  } else standatasml <- standata
  standatasml$priors <- 1L
  standatasml$nsubsets <- as.integer(nsubsets)
  if(optimcores > 1) parallelStanSetup(cl = cl,standata = standatasml,split=TRUE,nsubsets = nsubsets)
  if(optimcores==1) optimArgs$lpgFunc <- singlecoreStanSetup(standata = standatasml, nsubsets = nsubsets)
  optimArgs$tol <- optimArgs$tol * 100 
  optimArgs$maxiter <- 500
  optimArgs$nsubsets= nsubsets
  optimArgs$worsecountconverge <- 20
  fit = do.call(ctOptim,optimArgs)
  
  if(standata$ntipredeffects > 0 && notipredsfirstpass && !standata$TIpredAuto){
    message('Including tipred effects...')
    standata$TIPREDEFFECTsetup <- TIPREDEFFECTsetup
    standata$TIPREDEFFECTsetup[,] 
    standata$ntipredeffects <- ntipredeffects
    optimArgs$init <- c(optimArgs$init,rep(0,max(TIPREDEFFECTsetup)))
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
    # browser()
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
    if (optimcores == 1) optimArgsReduced$lpgFunc <- singlecoreStanSetup(standata = standata, nsubsets = 1)
    iter <- 0L
    # browser()
    optimfit <- do.call(ctOptim, optimArgsReduced)
    optimArgs$init[1:(nbasepars+found)] <- optimfit$par #update full init vec
  }
  
  return(list(standata   = standata,optimfit   = optimfit))
} # end ti pred auto function


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
#' @param tol objective tolerance.
#' @param priors logical. If TRUE, a priors integer is set to 1 (TRUE) in the standata object -- only has an effect if 
#' the stan model uses this value. 
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when priors=FALSE
#' @param subsamplesize value between 0 and 1 representing proportion of subjects to include in first pass fit. 
#' @param cores Number of cpu cores to use, should be at least 2.
#' @param bootstrapUncertainty Logical. If TRUE, subject wise gradient contributions are resampled to estimate the hessian, 
#' for computing standard errors or initializing importance sampling.
#' @param is Logical. Use importance sampling, or just return map estimates?
#' @param isloopsize Number of samples of approximating distribution per iteration of importance sampling.
#' @param finishsamples Number of samples to draw (either from hessian
#' based covariance or posterior distribution) for final results computation.
#' @param finishmultiply Importance sampling stops once available samples reach \code{finishsamples * finishmultiply} , then the final samples are drawn
#' without replacement from this set.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent
#' importance sampling, at risk of truncating tails.
#' @param parsteps ordered list of vectors of integers denoting which parameters should begin fixed
#' at zero, and freed sequentially (by list order). Useful for complex models, e.g. keep all cross couplings fixed to zero 
#' as a first step, free them in second step. 
#' @param parstepsAutoModel if TRUE, determines model structure for the parameters specified in parsteps automatically. If 'group', determines this on a group level first and then a subject level. Primarily for internal ctsem use, see \code{?ctFitAuto}.
#' @param groupFreeThreshold threshold for determining whether a parameter is free in a group level model. If the proportion of subjects with a non-zero parameter is above this threshold, the parameter is considered free. Only used with parstepsAutoModel = 'group'.
#' @param chancethreshold drop iterations of importance sampling where any samples are chancethreshold times more likely to be drawn than expected.
#' @param matsetup subobject of ctStanFit output. If provided, parameter names instead of numbers are output for any problem indications.
#' @param nsubsets number of subsets for stochastic optimizer. Subsets are further split across cores, 
#' but each subjects data remains whole -- processed by one core in one subset.
#' @param stochasticTolAdjust Multiplier for stochastic optimizer tolerance. 
#' @param lproughnesstarget target log posterior roughness for stochastic optimizer (suggest between .05 and .4).
#' @return list containing fit elements
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp

stanoptimis <- function(standata, sm, init='random',initsd=.01,
  estonly=FALSE,tol=1e-8,
  stochastic = TRUE,
  priors=TRUE,carefulfit=TRUE,
  bootstrapUncertainty=FALSE,
  subsamplesize=1,
  parsteps=c(),
  parstepsAutoModel=FALSE,
  groupFreeThreshold=.5,
  plot=FALSE,
  is=FALSE, 
  isloopsize=1000, 
  finishsamples=1000, 
  tdf=10,
  chancethreshold=100,
  finishmultiply=5,
  lproughnesstarget=.2,
  verbose=0,
  cores=2,
  matsetup=NA,
  nsubsets=10, 
  stochasticTolAdjust=1000){
  
  
  
  # initial checks ----------------------------------------------------------
  
  
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
  
  if(cores<2 && parstepsAutoModel %in% 'group') stop('parstepsAutoModel = "group" requires cores > 1')
  
  optimcores <- ifelse(length(unique(standata$subject)) < cores, length(unique(standata$subject)),cores)
  if(optimcores > 1) rm(smf)
  
  if(plot > 0 && .Platform$OS.type=="windows") {
    dev.new(noRStudioGD = TRUE)
    on.exit(expr = {try({dev.off()})},add = TRUE)
  }
  
  message('Using ',cores,'/', parallel::detectCores(),' CPU cores')
  
  storedPars <- as.numeric(c())
  storedLp <- c()
  
  optimfinished <- FALSE
  on.exit({
    if(!optimfinished){
      message('Optimization cancelled -- restart from current point by including this argument:')
      message((paste0(c('inits = c(',   paste0(round(storedPars,5),collapse=', '), ')'    ))))
    }},add=TRUE)
  
  
  # initial values ----------------------------------------------------------
  
  if(all(optimArgs$init %in% 'random')){
    optimArgs$init <- rnorm(npars, 0, initsd)
    if(length(parsteps)>0) optimArgs$init[unlist(parsteps)] <- 0 
  }
  
  if(all(optimArgs$init == 0)) optimArgs$init <- rep(0,npars)
  
  if(length(optimArgs$init) != npars){
    warning('Initialisation vector length does not match number of parameters in model, extending with zeros')
    optimArgs$init=c(optimArgs$init[1:min(length(optimArgs$init),npars)],rep(0,abs(npars-length(optimArgs$init))))
  }
  
  if(any(is.na(optimArgs$init))){
    warning('Initialisation vector contains NAs, replacing with zeros')
    optimArgs$init[is.na(optimArgs$init)] <- 0
  }
  
  # if(notipredsfirstpass && standata$ntipredeffects > 0){
  #   optimArgs$init[length(optimArgs$init):(length(optimArgs$init)+1-standata$ntipredeffects)] <- 0
  # }
  
  
  # initialise cluster ------------------------------------------------------
  
  if(optimcores > 1){ #for parallelised computation
    clctsem=makeClusterID(optimcores)
    on.exit(try({parallel::stopCluster(clctsem)},silent=TRUE),add=TRUE)
    
    if(standata$recompile > 0){
      smfile <- file.path(tempdir(),paste0('ctsem_sm_',ceiling(runif(1,0,100000)),'.rda'))
      save(sm,file=smfile,eval.promises = FALSE,precheck = FALSE)
      on.exit(add = TRUE,expr = {file.remove(smfile)})
    } else smfile <- ''
    
    parallel::clusterExport(clctsem,c('cores','smfile'),envir=environment())
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
    if(plot > 0 && ( (!stochastic &&!carefulfit && nsubsets ==1))){
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
      nsubsets=nsubsets,subsamplesize=subsamplesize,optimArgs=optimArgs,notipredsfirstpass=TRUE)
    optimArgs$init[1:length(optimfit$par)] <- optimfit$par #update non ti pred inits
  } #end carefulfit
  
  
  
  
  
  
  
  
  
  # end subsetting / carefulfit -----------------------------------------
  
  if(nsubsets > 1){ # need to reinit without subsets
    standata$nsubsets <- 1L
    optimArgs$nsubsets <- 1L
    if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=TRUE)
    if(optimcores==1) smf<-stan_reinitsf(sm,standata)
  }
  
  
  # tipredauto --------------------------------------------------------------
  
  
  if(standata$ntipred > 0 && !is.null(standata$TIpredAuto) && standata$TIpredAuto){
    if((length(parsteps) > 0)) stop('parsteps not supported with TIpredAuto')
    
    # insert TI‐pred logic via our new function
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
    if (optimcores == 1) optimArgs$lpgFunc <- singlecoreStanSetup(standata = standata, nsubsets = 1) #update single core lpg (parallel already updated)
  } #end ti pred auto total loop
  
  
  ##parameter stepwise / selection
  if(length(parsteps) > 0){
    if(parstepsAutoModel %in% FALSE){
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
    if(parstepsAutoModel %in% TRUE){
      # -----------------------------
      # Assume the following are available:
      #   - parFreeList: a list of length 2
      #         parFreeList[[1]]: vector of indices for basic parameters (always estimated)
      #         parFreeList[[2]]: vector of candidate parameter indices that are initially fixed
      #   - init: the current parameter vector (from the previous optimization stage)
      #   - lpgFunc: a function that returns the log probability, with an attribute "gradient"
      #   - jac: a function to compute a finite-difference approximation of a parameter's Hessian (diagonal element)
      #   - sgd: your stochastic optimizer that accepts an argument 'whichignore' (the parameters to keep fixed)
      #   - tol, stochasticTolAdjust, nsubsets, plot: parameters as in your code.
      #
      # Set a Wald threshold:
      wald_threshold <- 1.96  
      parstepsAutoModelOptimArgs <- optimArgs
      
      jacPars <- function(pars, step = 1e-3, whichpars) {
        # Initialize a vector to store the Hessian diagonal estimates for the specified parameters.
        hess_diag <- numeric(length(whichpars))
        # Loop over each requested parameter index.
        for (i in seq_along(whichpars)) {
          idx <- whichpars[i]
          # Create perturbed parameter vectors: one for the forward difference and one for the backward difference.
          pars_forward <- pars
          pars_backward <- pars
          pars_forward[idx] <- pars_forward[idx] + step
          pars_backward[idx] <- pars_backward[idx] - step
          # Evaluate the lpgFunc function at both perturbed vectors.
          # It is assumed the 'lpgFunc' function returns an object with the gradient as an attribute "gradient".
          forward_val <- optimArgs$lpgFunc(pars_forward)
          backward_val <- optimArgs$lpgFunc(pars_backward)
          # Extract the gradient for the current parameter.
          grad_forward <- attributes(forward_val)$gradient[idx]
          grad_backward <- attributes(backward_val)$gradient[idx]
          # Estimate the second derivative via the central difference formula.
          hess_diag[i] <- (grad_forward - grad_backward) / (2 * step)
        }
        return(hess_diag)
      }
      
      # Start with all candidate parameters fixed:
      currentFixed <- parsteps[[1]]
      # The permanently free parameter indices are:
      freePars <- (1:length(parstepsAutoModelOptimArgs$init))[!parsteps[[1]] %in% (1:length(parstepsAutoModelOptimArgs$init))]
      # Track which candidate parameters have been freed (initially, none)
      freed_candidates <- c()
      
      continueFreeing <- TRUE
      improvement_threshold <- 1.96
      while (continueFreeing && length(currentFixed) > 0) {
        
        # Evaluate the current log probability and obtain the gradient.
        parstepsAutoModelOptimArgs$whichignore <- currentFixed
        iter <-0
        optimfit <- do.call(ctOptim, parstepsAutoModelOptimArgs)
        
        # Create the list of indices to ignore in the next optimization step.
        # If no candidates remain, use an empty vector.
        ignore_indices <- if (length(currentFixed) > 0) currentFixed else integer(0)
        
        parstepsAutoModelOptimArgs$init[-ignore_indices] = optimfit$par
        obj_val <- optimArgs$lpgFunc(parstepsAutoModelOptimArgs$init)
        grad_vec <- attributes(obj_val)$gradient
        
        # For each candidate parameter currently fixed, compute expected improvement.
        improvement_est <- numeric(length(currentFixed))
        for (i in seq_along(currentFixed)) {
          idx <- currentFixed[i]
          
          # Use our simple jacPars to compute the approximate second derivative for this parameter.
          hess_est <- jacPars(parstepsAutoModelOptimArgs$init, step = 1e-6, whichpars = idx)
          
          # For stability, if the second derivative is NA or non-negative (which is
          # unexpected at a maximum), force a small negative value.
          if (is.na(hess_est) || hess_est >= 0) {
            hess_val <- -1e-6
          } else {
            hess_val <- hess_est
          }
          
          # Estimate expected improvement: 0.5 * (grad^2 / |H_ii|)
          improvement_est[i] <- 0.5 * (grad_vec[idx]^2 / abs(hess_val))
        }
        
        # Find the candidate with the highest expected improvement.
        best_candidate_index <- which.max(improvement_est)
        best_improvement <- improvement_est[best_candidate_index]
        if (best_improvement >= improvement_threshold) {
          best_param <- currentFixed[best_candidate_index]
          ms=data.frame(standata$matsetup)
          bestpar_ms <- ms[ms$param == best_param,,drop=FALSE][1,]
          message(paste0("Freeing parameter ",names(sort(ctStanMatricesList()$all))[bestpar_ms$matrix],'[',
            bestpar_ms$row,',',bestpar_ms$col,'] with expected improvement ',round(best_improvement,3)))
          
          # Mark the best candidate as freed.
          freed_candidates <- c(freed_candidates, best_param)
          # Remove it from the list of currently fixed candidate parameters.
          currentFixed <- currentFixed[-best_candidate_index]
          
        } else {
          message("No candidate parameter meets the improvement threshold. Terminating freeing sequence.")
          continueFreeing <- FALSE
          parsteps <- currentFixed
        }
      }
      
      # At this point, all parameters in `freePars` (the basic ones) and those in `freed_candidates`
      # are free (and have been re–optimized), while those remaining in `currentFixed` are kept fixed.
      # You can now proceed with the rest of your model estimation using 'init' as the final parameter estimate.
      parsteps <- currentFixed
    }
  }
  
  if (parstepsAutoModel %in% 'group') {
    # --------------------------------------------------
    # Group‐level stepwise freeing, with per‐subject re‐fit and init averaging
    groupParStepsOptimArgs <- optimArgs
    # thresholds
    improvement_threshold <- 1.96   # per‐subject ΔLL must exceed this
    
    # initial fixed candidates and subjects
    currentFixed <- parsteps[[1]]
    subject_ids  <- unique(standata$subject)
    
    
    parallel::clusterExport(clctsem, c(
      "tol", "stochasticTolAdjust",  "standatact_specificsubjects","subject_ids"
    ), envir = environment())
    
    continueFreeing <- TRUE
    subj_init <- matrix(rep(groupParStepsOptimArgs$init, length(subject_ids)), nrow = length(subject_ids), byrow = TRUE)
    while (continueFreeing && length(currentFixed) > 0) {
      # export updated init and fixed set to workers
      parallel::clusterExport(clctsem, c("currentFixed",'subj_init'), envir = environment())
      # 1) fit each subject with currentFixed held fixed
      subj_pars <- parallel::parLapply(clctsem, subject_ids, function(sid) {
        sd    <- standatact_specificsubjects(standata, sid)
        smf   <- stan_reinitsf(sm, sd)
        parlp <- function(parm){
          out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
          if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
            outerr <- out
            out <- -1e100
            attributes(out)$gradient <- rep(NaN, length(parm))
            attributes(out)$err <- outerr
          }
          if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
          return(out)
        }
        groupParStepsOptimArgs$init <- subj_init[sid,]
        groupParStepsOptimArgs$lpgFunc <- parlp
        iter <-0
        fit <- do.call(ctOptim, optimArgs)
        as.numeric(fit$par)
      })
      subj_mat <- do.call(rbind, subj_pars)
      subj_init[,-currentFixed] <- subj_mat # update init for next iteration
      
      # 2) average into init
      groupParStepsOptimArgs$init[-currentFixed] <- colMeans(subj_mat)
      
      # 3) compute per‐subject expected ΔLL for each candidate
      subj_imp <- parallel::parLapply(clctsem,seq_along(subject_ids), function(i) {
        sid  <- subject_ids[i]
        pvec <- groupParStepsOptimArgs$init
        pvec[-currentFixed] <- subj_mat[i, ]
        sd   <- standatact_specificsubjects(standata, sid)
        smf  <- stan_reinitsf(sm, sd)
        tgt  <- function(p) log_prob(smf, upars = p, adjust_transform = TRUE, gradient = TRUE)
        grad <- attributes(tgt(pvec))$gradient
        # helper: finite‐difference Hessian diagonal
        jacPars <- function(pars, step = 1e-3, whichpars) {
          sapply(whichpars, function(idx) {
            pf <- pars; pb <- pars
            pf[idx] <- pf[idx] + step
            pb[idx] <- pb[idx] - step
            gf <- attributes(tgt(pf))$gradient[idx]
            gb <- attributes(tgt(pb))$gradient[idx]
            (gf - gb) / (2 * step)
          })
        }
        sapply(currentFixed, function(idx) {
          h <- jacPars(pvec, step = 1e-6, whichpars = idx)
          if (is.na(h) || h >= 0) h <- -1e-6
          0.5 * (grad[idx]^2 / abs(h))
        })
      })
      # coerce to matrix if needed
      imp_vecs <- lapply(subj_imp, as.numeric)
      imp_mat  <- do.call(rbind, imp_vecs)
      if (is.null(dim(imp_mat))) {
        imp_mat <- matrix(imp_mat, nrow = length(imp_vecs), byrow = TRUE)
      }
      
      # 4) identify group‐level candidates
      prop_above  <- colMeans(imp_mat >= improvement_threshold)
      group_cands <- currentFixed[prop_above > groupFreeThreshold]
      
      if (length(group_cands) == 0) {
        message("No group‐level parameters exceed threshold; stopping.")
        break
      }
      
      ## pick group candidate with highest mean ΔLL
      # means      <- colMeans(imp_mat[, currentFixed %in% group_cands, drop = FALSE])
      # best_param <- group_cands[which.max(means)]
      
      # pick group candidate with highest proportion significant
      best_param <- currentFixed[which.max(prop_above)]
      
      message(sprintf(
        "Freeing parameter %d (%.0f%% subjects ΔLL ≥ %.2f)",
        best_param,
        100 * prop_above[currentFixed == best_param],
        improvement_threshold
      ))
      
      # update fixed set only
      currentFixed <- setdiff(currentFixed, best_param)
    }
    groupFixed <- currentFixed
    parallel::clusterExport(clctsem, c("groupFixed", "subj_init"), envir = environment())
    
    subj_res <- parallel::parLapply(clctsem, seq_along(subject_ids), function(i) {
      # for(i in 1:length(subject_ids)){
      # lapply(seq_along(subject_ids), function(i) {
      sid    <- subject_ids[i]
      sd     <- standatact_specificsubjects(standata, sid)
      smf    <- stan_reinitsf(sm, sd)
      p_i    <- subj_init[sid,]
      free_i <- setdiff(seq_along(groupParStepsOptimArgs$init), groupFixed)
      subjFixed <- groupFixed
      freed_i   <- integer(0)
      parlp <- function(parm){
        out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
          outerr <- out
          out <- -1e100
          attributes(out)$gradient <- rep(NaN, length(parm))
          attributes(out)$err <- outerr
        }
        if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
        return(out)
      }
      repeat {
        lpinit <- parlp(p_i)
        grad <- attributes(lpinit)$gradient
        impr <- sapply(subjFixed, function(idx) {
          pf <- p_i; pb <- p_i
          pf[idx] <- pf[idx] + 1e-6
          pb[idx] <- pb[idx] - 1e-6
          gf <- attributes(parlp(pf))$gradient[idx]
          gb <- attributes(parlp(pb))$gradient[idx]
          h  <- (gf - gb) / (2 * 1e-6)
          if (is.na(h) || h >= 0) h <- -1e-6
          0.5 * (grad[idx]^2 / abs(h))
        })
        if (all(impr < improvement_threshold)) break
        best_idx  <- subjFixed[which.max(impr)]
        freed_i   <- c(freed_i, best_idx)
        subjFixed <- setdiff(subjFixed, freed_i)
        if(length(subjFixed) == 0) break
        fit <- sgd(
          p_i+rnorm(length(p_i),0,.01), #init away from old max
          fitfunc = parlp,
          itertol     = tol * stochasticTolAdjust,
          maxiter     = 5000,
          whichignore = subjFixed,
          worsecountconverge = 20
        )
        p_i[sort(c(free_i,freed_i))] <- fit$par
      }
      list(par = p_i, freed = freed_i)
    })
    
    # build output matrices
    subjPars  <- t(sapply(subj_res, `[[`, "par"))
    subjFreed <- t(sapply(subj_res, function(x) groupFixed %in% x$freed))
    rownames(subjFreed) <- subject_ids
    colnames(subjFreed) <- as.character(groupFixed)
    
    # finalize
    parsteps  <- groupFixed
    optimfit  <- list(
      par       = groupParStepsOptimArgs$init[-parsteps],
      subjPars  = subjPars,
      subjFreed = subjFreed
    )
  }
  
  
  
  
  if(parstepsAutoModel %in% FALSE){
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
  } #end if not auto model parsteps
  
  
  bestfit <-optimfit$value
  est2=optimArgs$init #because init contains the fixed values #unconstrain_pars(smf, est1)
  if(length(parsteps)>0) est2[-parsteps] = optimfit$par else est2=optimfit$par
  
  npars = length(est2)
  
  if(!estonly){
    
    message('Estimating Hessian',appendLF = bootstrapUncertainty)
    plot=FALSE
    
    
    if(length(parsteps)>0) grinit= est2[-parsteps] else grinit = est2
    
    if(bootstrapUncertainty %in% 'auto'){
      if(standata$nsubjects < 5){
        bootstrapUncertainty <- FALSE
        message('Too few subjects (<5) for bootstrap uncertainty, classical hessian used')
      } else bootstrapUncertainty=TRUE
    }
    
    if(bootstrapUncertainty){
      scores <- scorecalc(standata = standata,est = grinit,stanmodel = sm,
        subjectsonly = standata$nsubjects > 5,returnsubjectlist = F,cores=cores)
      num_bootstrap_samples <- max(c(finishsamples,1000))
      alpha_max = 100 # Maximum bootstrap sample size factor
      alpha_min = 1 # Minimum bootstrap sample size factor
      n_threshold=1000 # Threshold n for alpha correction
      alpha <- alpha_max - (alpha_max - alpha_min) * (min(1000,nrow(scores))  / n_threshold)  # Bootstrap sample size factor
      num_bootstrap_samples  # Total number of bootstrap samples
      n <- nrow(scores)  # Number of observations
      p <- ncol(scores)  # Number of parameters
      
      # Create a bootstrap resampling matrix
      resample_matrix <- matrix(sample(1:n, size = round(alpha * n) * num_bootstrap_samples, replace = TRUE),
        nrow = num_bootstrap_samples, ncol = round(alpha * n))
      
      # Generate random weights for smoothing
      weights <- matrix(runif(length(resample_matrix), min = 0.1, max = 2),
        nrow = num_bootstrap_samples)
      
      # Aggregate gradients using matrix multiplication
      gradsamples <- matrix(0, nrow = num_bootstrap_samples, ncol = p)  # Initialize gradsamples
      
      # Compute the weighted sum of gradients for each bootstrap sample
      for (i in 1:num_bootstrap_samples) {
        gradsamples[i, ] <- colSums(scores[resample_matrix[i, ], , drop = FALSE] * weights[i, ])
      }
      
      # Trim outliers
      trim_percent <- 0.05  # Proportion of outliers to trim
      if (trim_percent > 0) {
        lower <- apply(gradsamples, 2, quantile, probs = trim_percent, na.rm = TRUE)
        upper <- apply(gradsamples, 2, quantile, probs = 1 - trim_percent, na.rm = TRUE)
        gradsamples <- pmax(pmin(gradsamples, upper), lower)        }
      
      #  Compute the Hessian with alpha correction
      hess <- -cov(gradsamples) / alpha
      
    }
    
    
    
    
    
    
    if(!bootstrapUncertainty){
      jac<-function(pars,step=1e-3,whichpars='all',
        lpdifmin=1e-5,lpdifmax=5, cl=NA,verbose=1,directions=c(-1,1),parsteps=c()){
        if('all' %in% whichpars) whichpars <- 1:length(pars)
        base <- optimfit$value
        
        
        hessout <- sapply( whichpars, function(i){
          
          
          
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
              lp[[di]] <-  suppressMessages(suppressWarnings(optimArgs$lpgFunc(pars+uppars*stepsize*directions[di])))
              accepted <- !'try-error' %in% class(lp[[di]]) && all(!is.na(attributes(lp[[di]])$gradient))
              if(accepted){
                lpdiff <- base[1] - lp[[di]][1]
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
                  base <- lp[[di]]
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
          
          grad<- attributes(lp[[1]])$gradient / steplist[[di]] * directions[di]
          if(any(is.na(grad))){
            warning('NA gradient encountered at param ',i,immediate. =TRUE)
          }
          if(length(directions) > 1) grad <- (grad + attributes(lp[[2]])$gradient / (steplist[[di]]*-1))/2
          return(grad)
        }
        ) #end sapply
        
        out=(hessout+t(hessout))/2
        return(out)
      }
      
      hess1 <- jac(pars = grinit,parsteps=parsteps,
        step = 1e-3,cl=NA,verbose=verbose,directions=1)
      hess2 <- jac(pars = grinit,parsteps=parsteps,#fgfunc = fgfunc,
        step = 1e-3,cl=NA,verbose=verbose,directions=-1)
      message('') #to create new line due to overwriting progress bar
      
      
      
      probpars <- c()
      onesided <- c()
      
      hess <- hess1
      hess[is.na(hess)] <- 0 #set hess1 NA's to 0
      hess[!is.na(hess2)] <- hess[!is.na(hess2)] + hess2[!is.na(hess2)] #add hess2 non NA's
      hess[!is.na(hess1) & !is.na(hess2)] <- hess[!is.na(hess1) & !is.na(hess2)] /2 #divide items where both hess1 and 2 used by 2
      hess[is.na(hess1) & is.na(hess2)] <- NA #set NA when both hess1 and 2 NA
      
      if(any(is.na(c(diag(hess1),diag(hess2))))){
        if(any(is.na(hess))) message ('Problems computing Hessian...')
        onesided <- which(sum(is.na(diag(hess1) & is.na(diag(hess2)))) %in% 1)
      }
      
      hess <- (t(hess)+hess)/2
      
      
      if(length(c(probpars,onesided)) > 0){
        if('data.frame' %in% class(matsetup) || !all(is.na(matsetup[1]))){
          ms=matsetup
          ms=ms[ms$param > 0 & ms$when == 0,]
          ms=ms[!duplicated(ms$param),]
        }
        if(length(onesided) > 0){
          onesided=paste0(ms$parname[ms$param %in% onesided],collapse=', ')
          message ('One sided Hessian used for params: ', onesided)
        }
        if(length(probpars) > 0){
          probpars=paste0(ms$parname[ms$param %in% c(probpars)],collapse=', ')
          message('***These params "may" be not identified: ', probpars)
        }
      }
    } #end classical hessian
    
    # cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
    
    # if('try-error' %in% class(cholcov) && !is) message('Approximate hessian used for std error estimation.')
    
    mcov=try(solve(-hess),silent=TRUE)
    if('try-error' %in% class(mcov)){
      mcov=MASS::ginv(-hess) 
      warning('***Generalized inverse required for Hessian inversion -- interpret standard errors with caution. Consider simplification, priors, or alternative uncertainty estimators',call. = FALSE,immediate. = TRUE)
      probpars=which(diag(hess) > -1e-6)
    }
    
    
    mcovtmp=try({as.matrix(Matrix::nearPD(mcov,conv.norm.type = 'F')$mat)})
    if(any(class(mcovtmp) %in% 'try-error')) stop('Hessian could not be computed')
    mcov <- diag(1e-10,npars)
    if(length(parsteps)>0) mcov[-parsteps,-parsteps] <- mcovtmp else mcov <- mcovtmp
    mchol = t(chol(mcov))
    
    mcovl<-list()
    mcovl[[1]]=mcov
    delta=list()
    delta[[1]]=est2
    samples <-matrix(NA)
    resamples <- c()
    prop_dens <-c()
    target_dens<-c()
    sample_prob<-c()
    ess <- 0
    qdiag<-0
    
    if(!is) {
      nresamples = finishsamples
      resamples <- matrix(unlist(lapply(1:nresamples,function(x){
        delta[[1]] + (mchol) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
      } )),byrow=TRUE,ncol=length(delta[[1]]))
    }
    message('')
    
    
    
    # finish split and sum across core approach, now importance sampling / par computations ----------------------------------
    
    #configure each node with full dataset for adaptive sampling
    if(cores > optimcores) clctsem=makeClusterID(cores) #reinit cluster if more cores available than used for optimising
    if(cores > 1)   parallelStanSetup(cl=clctsem,standata,split=FALSE)    
    
    
    if(is){
      
      message('Importance sampling...')
      
      log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
      }
      
      
      targetsamples <- finishsamples * finishmultiply
      j <- 0
      while(nrow(samples) < targetsamples){
        j<- j+1
        if(j==1){
          samples <- newsamples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
        } else {
          delta[[j]]=colMeans(resamples)
          mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
          
          for(iteri in 1:j){
            if(iteri==1) mcov=mcovl[[j]] else mcov = mcov*.5+mcovl[[j]]*.5 #smoothed covariance estimator
            if(iteri==1) mu=delta[[j]] else mu = mu*.5+delta[[j]]*.5 #smoothed means estimator
          }
          newsamples <- mvtnorm::rmvt(isloopsize, delta = mu, sigma = mcov,   df = tdf)
          samples <- rbind(samples, newsamples)
        }
        
        prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
        
        if(cores > 1) parallel::clusterExport(clctsem,c('samples'),envir = environment())
        
        
        target_dens[[j]] <- unlist(flexlapplytext(cl = clctsem, 
          X = 1:isloopsize, 
          fn = "function(x){parlp(samples[x,])}",cores=cores))
        
        
        
        target_dens[[j]][is.na(target_dens[[j]])] <- -1e200
        if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
        
        
        
        targetvec <- unlist(target_dens)
        
        target_dens2 <- target_dens[[j]] 
        target_dens2[!is.finite(target_dens[[j]])] <- -1e30
        weighted_dens <- target_dens2 - prop_dens
        
        
        newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
        counter <- 1
        isfinished <- FALSE
        
        sample_prob <- c(sample_prob,newsampleprob) #sum to 1 for each iteration, normalise later
        sample_prob[!is.finite(sample_prob)] <- 0
        sample_prob[is.na(sample_prob)] <- 0
        
        if(nrow(samples) >= targetsamples && (max(sample_prob)/ sum(sample_prob)) < chancethreshold) {
          message ('Finishing importance sampling...')
          nresamples <- finishsamples 
        }else nresamples = max(5000,nrow(samples)/5)
        
        
        resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(nrow(samples) > targetsamples,FALSE,TRUE),
          prob = sample_prob / sum(sample_prob))
        
        message(paste0('Importance sample loop ',j,', ',length(unique(resample_i)), ' unique samples, from ', nresamples,' resamples of ', nrow(samples),' actual, prob sd = ', round(sd(sample_prob),4),
          ', max chance = ',max(sample_prob) * isloopsize))
        if(length(unique(resample_i)) < 100) {
          message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
        }
        
        resamples <- samples[resample_i, , drop = FALSE]
        
        ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
        qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
          (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) )
        })))
        
      }
    }
  }
  
  if(!estonly){
    if(!is) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]
    
    message('Computing posterior with ',nrow(resamples),' samples')
    standata$savesubjectmatrices=savesubjectmatrices
    
    if(!savesubjectmatrices) sdat=standatact_specificsubjects(standata,1) #only use 1 subject
    if(savesubjectmatrices) sdat=standata
    
    transformedpars=stan_constrainsamples(sm = sm,standata = sdat,
      savesubjectmatrices = savesubjectmatrices, savescores = standata$savescores,
      dokalman=as.logical(standata$savesubjectmatrices),
      samples=resamples,cores=cores, cl=clctsem, quiet=TRUE)
    
    if(cores > 1) {
      parallel::stopCluster(clctsem)
      smf <- stan_reinitsf(sm,standata)
    }
    
    
    
    # transformedparsfull=stan_constrainsamples(sm = sm,standata = standata,
    #   savesubjectmatrices = TRUE, dokalman=TRUE,savescores = TRUE,
    #   samples=matrix(est2,nrow=1),cores=1, quiet = TRUE)
    
    
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if('try-error' %in% class(sds)[1]) sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, upars=lest)),
      unlist(constrain_pars(smf, upars= est2)),
      unlist(constrain_pars(smf, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    stanfit=list(optimfit=optimfit,stanfit=stan_reinitsf(sm,standata), rawest=est2, rawposterior = resamples, cov=mcov,
      transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      # transformedparsfull=transformedparsfull,
      standata=list(TIPREDEFFECTsetup=standata$TIPREDEFFECTsetup,ntipredeffects = standata$ntipredeffects),
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
    if(bootstrapUncertainty) stanfit$subjectscores <- scores #subjectwise gradient contributions
  }
  
  if(estonly) {
    smf <- stan_reinitsf(sm,standata)
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2,parsteps=parsteps)
  }
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}



