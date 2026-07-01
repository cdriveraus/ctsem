ctModeltoNumeric <- function(ctmodelobj){
  ###read in model
  #set any matrices to numeric elements
  domessage<-FALSE
  sapply(names(ctmodelobj), function(x){
    if(is.matrix(ctmodelobj[[x]])){
      m <- ctmodelobj[[x]]
      if(any(suppressWarnings(is.na(as.numeric(m))))){
        domessage <<-TRUE
      }
      m[suppressWarnings(is.na(as.numeric(m)))] <- 0
      ctmodelobj[[x]] <<- matrix(as.numeric(m),nrow=nrow(m), ncol=ncol(m))
    }
  })
  if(domessage) message('Some free parameters in matrices set to zero!')
  
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

#' Backward-compatible alias for \code{ctGenerateFromPriors}.
#' @rdname ctGenerateFromPriors
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
  Tpoints=NULL, wide=FALSE){
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
