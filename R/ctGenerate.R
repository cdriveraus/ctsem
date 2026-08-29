ctModeltoNumeric <- function(ctmodelobj){
  ###read in model
  #set any matrices to numeric elements
  #
  # Free parameters need *a* value before anything can be simulated, and zero is
  # the wrong one nearly everywhere. A zero DRIFT is singular, so `fQinf()`
  # cannot solve for the asymptotic covariance and generation fails outright; a
  # zero DIFFUSION is a process with no innovation; a zero MANIFESTVAR is
  # noiseless measurement; a zero LAMBDA disconnects a manifest from its latent.
  # Only the location matrices -- means and intercepts -- are naturally zero.
  #
  # `.ctGenerateDefaults()` holds the per-matrix choices, shared with the julia
  # generation path so both produce the same kind of data from an underspecified
  # model. They are defaults for *simulation*, not estimates: the point is that
  # generated data looks like data rather than like an artefact, and anything a
  # user cares about they should set.
  defaults <- .ctGenerateDefaults()
  filled <- character()
  sapply(names(ctmodelobj), function(x){
    if(is.matrix(ctmodelobj[[x]])){
      m <- ctmodelobj[[x]]
      free <- suppressWarnings(is.na(as.numeric(m)))
      dim(free) <- dim(m)
      if(any(free)){
        spec <- defaults[[x]]
        idx <- which(free, arr.ind=TRUE)
        for(k in seq_len(nrow(idx))){
          i <- idx[k,1]; j <- idx[k,2]
          value <- if(is.null(spec)) 0 else
            if(i == j) spec$diagonal else spec$offdiagonal
          m[i,j] <- value
          if(value != 0) filled <<- c(filled, sprintf('%s[%d,%d]=%s', x, i, j,
            format(value)))
        }
      }
      ctmodelobj[[x]] <<- matrix(as.numeric(m),nrow=nrow(m), ncol=ncol(m))
    }
  })
  if(length(filled)) message('Free parameters were given generating values: ',
    paste(utils::head(filled, 8), collapse=', '),
    if(length(filled) > 8) ', ...' else '',
    '. Others were set to zero. Set them in the model if they matter.')
  
  return(ctmodelobj)
}

#' Generate data from a ctstanmodel object
#'
#' @param cts \code{\link{ctModelConvertOMX}}, \code{\link{ctModel}}, or
#' \code{\link{ctStanFit}} object.
#' @param datastruct long format data structure as used by ctsem. 
#' Not used if cts is a ctStanFit object.
#' @param is If optimizing, follow up with importance sampling? 
#' @param fullposterior Generate from the full posterior or just the (unconstrained) mean?
#' @param nsamples How many samples to generate?
#' @param parsonly If TRUE, only return samples of raw parameters, don't generate data.
#' @param cores Number of cpu cores to use.
#'
#' @return List contining Y, and array of nsamples by data rows by manifest variables, 
#' and llrow, an array of nsamples by data rows log likelihoods.
#' @aliases ctStanGenerate
#' @export
#'
#' @examples
#' \donttest{
#' #generate and plot samples from prior predictive
#' priorpred <- ctGenerateFromPriors(cts = ctstantestfit,cores=2,nsamples = 50)
#'}
ctGenerateFromPriors <- function(cts,datastruct=NA, is=FALSE,
  fullposterior=TRUE, nsamples=200, parsonly=FALSE,cores=2){
  
  # includePreds <- FALSE #old argument, could reinstate some day...
  #update this function to also generate posterior predictive
  
  # nopriors <- FALSE # update this when creating posterior predictive, go to TRUE if fullposterior=F and fit object had no priors
  
  if('ctStanFit' %in% class(cts)){
    # if(!fullposterior && cts$standata$nopriors==1) nopriors <- TRUE #generate from point estimate
    priors <- cts$args$priors
    datastruct <- standatatolong(cts$standata, origstructure=TRUE, ctm=cts$ctstanmodelbase)
    
    cts <- cts$ctstanmodelbase
    # browser() 
    #  if(cts$setup$recompile){ #then temporarily attach compiled stanmodels to search path to avoid recompiling
    #    ctsem.compiledmodel <- new.env()
    #    ctsem.compiledmodel$fitmodel <- cts$stanmodel
    #    if(!is.null(cts$generated)) ctsem.compiledmodel$genmodel <- cts$generated$stanmodel
    #    attach(ctsem.compiledmodel)
    #    on.exit(add = TRUE, {detach(name = 'ctsem.compiledmodel')})
    #    }
    
  } else priors<-TRUE
  
  datastruct[,cts$manifestNames] <- NA #remove manifest variables
  optimcontrol<- list()
  optimcontrol$carefulfit=FALSE
  optimcontrol$is <- is
  optimcontrol$stochastic=FALSE
  optimcontrol$finishsamples=nsamples
  
  
  cts$TIpredAuto <- 0L
  
  ds <- data.table(datastruct)
  ds[,WhichObs:=(1:.N),by=eval(cts$subjectIDname)]
  datadummy= data.frame(datastruct)[ds$WhichObs==1,]
  datadummy[,cts$TIpredNames] <- 0
  
  
  args <- cts$args
  args$optimcontrol=optimcontrol
  args$optimize=TRUE
  args$cores=cores
  args$model <- cts
  args$ctstanmodel <- NULL
  args$intoverstates <- TRUE
  args$intoverpop <- TRUE
  args$inits=0
  args$datalong=datadummy
  args$priors <- priors
  args$optimcontrol=list(stochastic=FALSE,carefulfit=FALSE)
  if(!is.null(args$priors) && !as.logical(args$priors)) stop('Priors disabled, cannot sample from prior!')
  
  #fit to empty data 
  message('Fitting model to empty dataset...')
  
  pp<-do.call(ctFit,args)
  
  if(parsonly) dat <- pp else{
    
    datastruct[,cts$manifestNames] <- -99
    
    #get filled standata object
    pp$standata<-ctStanData(ctm=pp$ctstanmodel, datalong=datastruct,optimize=TRUE)
    
    ppf <- ctGenerateFromFit(fit = pp,nsamples = nsamples,fullposterior = fullposterior,cores=cores)
    
    #collect generated stuff
    dat <-list()
    dat$Y <- ppf$generated$Y
    dimnames(dat$Y) <- list(datapoints=1:dim(dat$Y)[1], samples=1:dim(dat$Y)[2], manifests = cts$manifestNames)
    dat$llrow <- ppf$generated$llrow
  }
  
  
  return(dat)
}

#' @export
ctStanGenerate <- ctGenerateFromPriors



#' ctGenerate
#' 
#' This function generates data according to the specified ctsem model object. 
#' 
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dtmean Positive numeric. Average time interval (delta T) to use.
#' @param logdtsd Numeric. Standard deviation for variability of the time interval.
#' @param dtmat Either NA, or numeric matrix of n.subjects rows and Tpoints-1 columns, 
#' containing positive numeric values for all time intervals between measurements. 
#' If not NA, dtmean and logdtsd are ignored.
#' @param Tpoints Optional number of time points to generate. If supplied, this overrides
#' any \code{Tpoints} stored in \code{ctmodelobj}. If not supplied, \code{ctGenerate}
#' uses \code{ctmodelobj$Tpoints} when available.
#' @param wide Logical. Output in wide format?
#' @details Covariance related matrices are treated as Cholesky factors. 
#' TRAITTDPREDCOV and TIPREDCOV matrices are not accounted for, at present. 
#' The first 1:n.TDpred rows and columns of TDPREDVAR are used for generating
#' tdpreds at each time point. 
#' @examples 
#' #generate data for 2 process model, each process measured by noisy indicator, 
#' #stable individual differences in process levels.
#' 
#' generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'  MANIFESTVAR=diag(.1,2),
#'  LAMBDA=diag(1,2),
#'  DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
#'  DIFFUSION=matrix(c(1,.2,0,4),2),
#'  CINT=matrix(c(1,0),nrow=2),
#'  T0MEANS=matrix(0,ncol=1,nrow=2),
#'  T0VAR=diag(1,2))
#'
#' nsubjects <- 15
#' traitChol <- matrix(c(.5,.2,0,.8),nrow=2)
#' subjectCint <- t(replicate(nsubjects, as.numeric(traitChol %*% rnorm(2))))
#' datalist <- vector("list", nsubjects)
#' for(i in seq_len(nsubjects)){
#'   subjectModel <- generatingModel
#'   subjectModel$CINT <- matrix(subjectCint[i,], ncol = 1)
#'   d <- ctGenerate(subjectModel,n.subjects=1,burnin=10)
#'   d[,'id'] <- i
#'   datalist[[i]] <- d
#' }
#' data <- do.call(rbind, datalist)
#' @export

ctGenerate<-function(ctmodelobj,n.subjects=100,burnin=0,dtmean=1,logdtsd=0,dtmat=NA,
  Tpoints=NULL, wide=FALSE, backend=c('auto','r','julia')){
  backend <- match.arg(backend)
  # `auto` routes to julia only what the generator below cannot do. That
  # generator integrates the linear system with a matrix exponential, which is
  # exact for a linear model and simply inapplicable to a state-dependent one;
  # the engine filters and generates both. Defaulting to julia for everything
  # would change the numbers under every existing caller for no gain on the
  # models they use, so the split is by capability rather than by preference.
  nonlinear <- isTRUE(try(ctModelIsNonlinear(ctmodelobj), silent=TRUE))
  # A binary indicator is the same situation as a nonlinear one: the generator
  # below integrates a linear Gaussian system and has no notion of a link, so
  # it produces continuous values for a manifest the model declares binary --
  # silently, which is the worst of both. Measured before this line existed: a
  # model with `manifesttype = 1` generated values with a mean of 0.063 and no
  # zeros or ones among them.
  binary <- !is.null(ctmodelobj$manifesttype) && any(ctmodelobj$manifesttype > 0)
  if(backend == 'auto') backend <- if(nonlinear || binary) 'julia' else 'r'
  if(backend == 'r' && binary){
    warning('This model declares binary indicators (manifesttype = 1) and ',
      "backend='r' generates continuous values for them: the R generator has ",
      'no measurement link. Use backend="julia" for binary data.',
      call.=FALSE)
  }
  if(backend == 'r' && nonlinear) {
    stop("This model is nonlinear, and ctGenerate's own generator integrates a ",
      "linear system: it has no way to apply a state-dependent specification. ",
      "Use backend='julia' (the default for such models), which generates ",
      "through the same filter that fits them.", call.=FALSE)
  }
  if(backend == 'julia'){
    if(!'ctStanModel' %in% class(ctmodelobj)) {
      stop("backend='julia' generation needs a ctModel(type='ct'/'dt') object, ",
        "which carries the parameter specification the engine reads. The ",
        "matrix-list form does not.", call.=FALSE)
    }
    modelTpoints <- if(!is.null(Tpoints) && !is.na(Tpoints[1])) Tpoints[1] else
      if(!is.null(ctmodelobj$Tpoints) && !is.na(ctmodelobj$Tpoints[1]))
        ctmodelobj$Tpoints[1] else
          stop('Tpoints not found in ctmodelobj and no Tpoints argument supplied. Provide Tpoints explicitly.')
    fullTpoints <- burnin + as.integer(modelTpoints)
    times <- lapply(seq_len(n.subjects), function(si){
      dtvec <- if(is.na(dtmat[1])) exp(rnorm(fullTpoints,log(dtmean),logdtsd)) else
        c(rep(1,burnin), dtmat[si,,drop=TRUE])
      tv <- numeric(fullTpoints)
      for(t in 2:fullTpoints) tv[t] <- round(tv[t-1] + dtvec[t-1], 6)
      tv
    })
    out <- .ctGenerateJulia(ctmodelobj, n.subjects, times)
    if(burnin > 0){
      keep <- unlist(lapply(seq_len(n.subjects), function(si)
        (si-1)*fullTpoints + (burnin+1):fullTpoints))
      out <- out[keep, , drop=FALSE]
      # Time restarts at zero for each subject once the burnin is dropped, as
      # the generator below does. Leaving it running from the burnin would make
      # the first observed interval look like the whole burnin period.
      for(si in seq_len(n.subjects)){
        rows <- (si-1)*(fullTpoints-burnin) + seq_len(fullTpoints-burnin)
        out[rows,'time'] <- out[rows,'time'] - out[rows[1],'time']
      }
    }
    if(wide) return(ctLongToWide(out, id='id', time='time',
      manifestNames=ctmodelobj$manifestNames,
      TDpredNames=ctmodelobj$TDpredNames, TIpredNames=ctmodelobj$TIpredNames))
    return(out)
  }
  if('ctStanModel' %in% class(ctmodelobj)){
    # Reconstruct matrix-style slots when a ctStanModel is supplied.
    mlist <- listOfMatrices(ctmodelobj$pars)
    for(nm in names(mlist)){
      ctmodelobj[[nm]] <- mlist[[nm]]
    }
  }
  
  ctmodelobj <- ctModeltoNumeric(ctmodelobj)
  
  m <- ctmodelobj
  
  modelTpoints <- NULL
  if(!is.null(m$Tpoints) && !is.na(m$Tpoints[1])) modelTpoints <- m$Tpoints[1]
  if(!is.null(Tpoints) && !is.na(Tpoints[1])) modelTpoints <- Tpoints[1]
  
  if(is.null(modelTpoints)) stop('Tpoints not found in ctmodelobj and no Tpoints argument supplied. Provide Tpoints explicitly.')
  if(length(modelTpoints) != 1 || !is.finite(modelTpoints) || modelTpoints < 1) stop('Tpoints must be a single finite value >= 1.')
  
  modelTpoints <- as.integer(modelTpoints)
  m$Tpoints <- modelTpoints
  fullTpoints<-burnin+m$Tpoints

  for(si in 1:n.subjects){
    
    if(is.na(dtmat[1])) dtvec<- exp(rnorm(fullTpoints,log(dtmean),logdtsd))
    if(!is.na(dtmat[1])) dtvec <- c(rep(1,burnin),dtmat[si,,drop=FALSE])
    time=rep(0,fullTpoints)
    for(t in 2:fullTpoints) time[t] = round(time[t-1] + dtvec[t-1],6)
    
    if(m$n.TDpred > 0) {
      tdpreds <- rbind(matrix(0,nrow=1+(burnin),ncol=m$n.TDpred)[-1,,drop=FALSE], #additional row added then removed in case no burnin
        matrix(m$TDPREDMEANS,ncol=m$n.TDpred))
    }
    
    # #convert to triangular...
    # m$T0VAR <- t(chol(m$T0VAR))
    # m$MANIFESTVAR <- t(chol(m$MANIFESTVAR))
    
    sm=m
    if(any(m$TRAITVAR != 0)) {
      traits = m$TRAITVAR %*% rnorm(m$n.latent,0,1)
      sm$CINT = m$CINT +  traits
      sm$T0MEANS = m$T0MEANS + m$T0TRAITEFFECT %*% traits
    }
    
    if(any(m$MANIFESTTRAITVAR != 0)) {
      sm$MANIFESTMEANS = m$MANIFESTMEANS + m$MANIFESTTRAITVAR %*% rnorm(m$n.manifest,0,1)
    }
    
    if(m$n.TIpred > 0) {
      tipreds <- m$TIPREDMEANS + m$TIPREDVAR %*% rnorm(m$n.TIpred,0,1)
      sm$CINT = m$CINT + m$TIPREDEFFECT %*% tipreds
    }
    
    manifests<-matrix(NA,fullTpoints,m$n.manifest)
    latents<-matrix(NA,fullTpoints,m$n.latent)
    
    sdat <- cbind(si,time,manifests,
      if(m$n.TDpred > 0) tdpreds,
      if(m$n.TIpred > 0) matrix(tipreds,byrow=TRUE,nrow=fullTpoints,ncol=m$n.TIpred))
    
    colnames(sdat) <- c('id','time',m$manifestNames,m$TDpredNames,m$TIpredNames)
    sdat <- data.frame(sdat)

    
    latents[1,] <- sm$T0MEANS+m$T0VAR %*% rnorm(m$n.latent)
    Qinf <- fQinf(sm$DRIFT,sm$DIFFUSION)
    
    for(i in 2:nrow(latents)){
      dtA=expm::expm(sm$DRIFT * (sdat$time[i]-sdat$time[i-1]))
      # message('dtA')
      # print(dtA)
      # message('dtCINT')
      # print(solve(sm$DRIFT,(dtA - diag(m$n.latent))) %*% sm$CINT)
      # message('dtG')
      # print(t(chol(fdtQ(Qinf,dtA))))
      latents[i,] <- dtA %*% latents[i-1,] +
        solve(sm$DRIFT,(dtA - diag(m$n.latent))) %*% sm$CINT + 
        t(chol(fdtQ(Qinf,dtA))) %*% rnorm(m$n.latent)
      # browser()
      if(m$n.TDpred > 0) latents[i,] <- latents[i,] + sm$TDPREDEFFECT %*% 
        t(as.matrix(sdat[i,m$TDpredNames,drop=FALSE]))
    }
    
    for(i in seq_len(nrow(sdat))){
      sdat[i,m$manifestNames] <- sm$LAMBDA %*% latents[i,] + sm$MANIFESTMEANS + 
        sm$MANIFESTVAR %*% rnorm(m$n.manifest)
    }
        
    
    
    
    sdat=sdat[(burnin+1):fullTpoints,]
    
    sdat[,'time'] = sdat[,'time'] - sdat[1,'time'] 
    
    if(si==1) datalong <- sdat else datalong <- rbind(datalong,sdat)
  }
  
  datalong<-as.matrix(datalong)
  
  
  if(wide==FALSE) return(datalong) else {
    datawide <- ctLongToWide(datalong = datalong,id = 'id',time = 'time',
      manifestNames = m$manifestNames, TDpredNames = m$TDpredNames,TIpredNames = m$TIpredNames)
    datawide <- ctIntervalise(datawide = datawide,Tpoints = m$Tpoints,n.manifest = m$n.manifest,n.TDpred = m$n.TDpred,n.TIpred = m$n.TIpred,
      manifestNames=m$manifestNames,TDpredNames=m$TDpredNames,TIpredNames=m$TIpredNames)
    return(datawide)
  }
}
