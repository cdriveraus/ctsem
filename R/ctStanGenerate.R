generator2 <- function(gm,nsubjects,
  maxtime, genstep=1e-3, obsstep=100,
  interventiontimes=c(), burnin=0){
  
  times <- seq(0,maxtime,genstep) #create the vector of observation times
  
  intervention <- rep(0,length(times)) #at nearly all times, there is no intervention happening
  if(length(interventiontimes)>1) intervention[interventiontimes / genstep+1] <- 1
  
  for(si in 1:nsubjects){
  #latent process loop
  latents <- matrix(gm$T0MEANS,nrow=1) #initialise our latent processes with the start value
  for(i in 2:length(times)){ #then calculate forward in time the subjects process values
    latents <- rbind(latents, #put the previous values of the latent processes on top with row bind
      c(latents[i-1,] + #then the new latents are based off the previous values
          (gm$DRIFT %*% latents[i-1,] + gm$CINT) * genstep + #deterministic change
          sysnoise %*% rnorm(length(starts),0,sqrt(genstep))))
    if(length(interventiontimes)>0) latents[i,]=latents[i,]+gm$TDPREDEFFECT %*% intervention[i]
  }#end latent loop
  
  #observations
  for(i in 1:length(times)){
    newobs <- gm$LAMBDA %*% latents[i,] +#our factor loadings transfer the latent process into observed measures
      gm$MANIFESTMEANS + gm$MANIFESTVAR %*% rnorm(length(obsNames),0,1) #which also can have some intercept and error variation
    if(i==1) obs <- matrix(newobs,nrow=1) else obs <- rbind(obs, c(newobs)) #if first row, we need to create the object, else we append below again
  }
  
  colnames(latents)<-gm$latentNames
  colnames(obs) <- gm$manifestNames
  
  
  dat <- data.frame(id=si,time=times, obs) #put all our data together to output it
  if(length(interventiontimes)>0) dat <- cbind(dat,intervention)
  dat <- dat[dat$time %in% seq(0,maxtime,obsstep),]
  }
  if(si==1) fulldat <- dat else fulldat <- rbind(fulldat,dat)
  fulldat <- fulldat[fulldat$time >= burnin,]
  return(fulldat)
}



#' Generate data from a ctstanmodel object
#'
#' @param cts \code{\link{ctStanModel}} , or \code{\link{ctStanFit}},object.
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
#' priorpred <- ctStanGenerate(cts = ctstantestfit,cores=2,nsamples = 50)
#'}
ctStanGenerate <- function(cts,datastruct=NA, is=FALSE, 
  fullposterior=TRUE, nsamples=200, parsonly=FALSE,cores=2){
  
  # includePreds <- FALSE #old argument, could reinstate some day...
  #update this function to also generate posterior predictive
  
  # nopriors <- FALSE # update this when creating posterior predictive, go to TRUE if fullposterior=F and fit object had no priors
  derrind <- 'all' #possibly update below
  if('ctStanFit' %in% class(cts)){
    # if(!fullposterior && cts$standata$nopriors==1) nopriors <- TRUE #generate from point estimate
    derrind <- cts$standata$derrind
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
  args$cores=cores #problem with multiple cores inside function?
  args$ctstanmodel <- cts
  args$intoverstates <- TRUE
  args$intoverpop <- TRUE
  args$inits=1e-10
  args$datalong=datadummy
  args$priors <- priors
  if(!is.null(args$priors) && !as.logical(args$priors)) stop('Priors disabled, cannot sample from prior!')

  #fit to empty data 
  message('Fitting model to empty dataset...')
  pp<-do.call(ctStanFit,args)
  
  if(parsonly) dat <- pp else{

    datastruct[,cts$manifestNames] <- -99

    #get filled standata object
    pp$standata<-ctStanData(ctm=pp$ctstanmodel, datalong=datastruct,optimize=TRUE,derrind= derrind)

    ppf <- ctStanGenerateFromFit(fit = pp,nsamples = nsamples,fullposterior = fullposterior,cores=cores)
    
    #collect generated stuff
    dat <-list()
    dat$Y <- ppf$generated$Y
    dimnames(dat$Y) <- list(datapoints=1:dim(dat$Y)[1], samples=1:dim(dat$Y)[2], manifests = cts$manifestNames)
    dat$llrow <- ppf$generated$llrow
  }
  
  
  return(dat)
}

