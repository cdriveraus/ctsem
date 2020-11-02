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
#' @param ctm \code{\link{ctStanModel}} object.
#' @param datastruct long format data structure as used by ctsem.
#' @param optimize Whether to optimize or use Stan's HMC sampler
#' @param is If optimizing, follow up with importance sampling? 
#' @param fullposterior Generate from the full posterior or just the mean?
#' @param nsamples How many samples to generate?
#' @param parsonly If TRUE, only return samples of raw parameters, don't generate data.
#' @param includePreds if TRUE, the prior for covariate effects (TD and TI predictors)
#' is included, as well as the TD and TI pred data. Else the effects are set to zero.
#' @param ... arguments to pass to stanoptimis 
#'
#' @return Array of nsamples x time points x manifest variables.
#' @export
#'
#' @examples
#' \donttest{
#' #generate and plot samples from prior predictive
#' priorpred <- ctStanGenerate(ctm = ctstantestfit$ctstanmodelbase,
#'   datastruct = ctstantestdat,cores=2,nsamples = 50)
#'}
ctStanGenerate <- function(ctm,datastruct, optimize=TRUE, is=FALSE, 
  fullposterior=TRUE, nsamples=200, parsonly=FALSE,includePreds = FALSE,...){
  
  datastruct[,ctm$manifestNames] <- NA
  dots <- list(...)
  dots$carefulfit=FALSE
  dots$is <- is
  dots$tol=1e-18
  dots$stochastic=FALSE
  if(is.null(dots$finishsamples) && parsonly) dots$finishsamples=nsamples
  if(!includePreds){
    ctm$n.TDpred <- 0
    ctm$TDpredNames <- NULL
    ctm$n.TIpred <- 0
    ctm$TIpredNames <- NULL
  }
  ctm$TIpredAuto <- 0L
  
  #problem with multiple cores inside function?
  dots1=dots
  dots1$cores=1
  
  datadummy= datastruct[c(1,nrow(datastruct)),,drop=FALSE]
  datadummy[,ctm$TIpredNames] <- 0
  pp<-ctStanFit(datalong =datadummy, #reenable multi core -- check parallel craziness in stanoptimis
    ctstanmodel = ctm,optimize=optimize, optimcontrol=dots1,verbose=0, nopriors=FALSE)
  
  if(parsonly) dat <- pp else{
    
    filled <- datastruct
    filled[,ctm$manifestNames] <- -99
    # browser()
    ppf <- ctStanFit(datalong = filled, ctstanmodel = ctm,optimize=optimize, 
      optimcontrol=dots1,fit=FALSE,nopriors=FALSE)
    # pp$standata$Y <- ppf$standata$Y
    ppf$stanfit <- pp$stanfit
    class(ppf) <- c('ctStanFit',class(ppf))
    ppf <- ctStanGenerateFromFit(fit = ppf,nsamples = nsamples,fullposterior = fullposterior)
    
    dat <- ppf$generated$Y
    dimnames(dat) <- list(datapoints=1:dim(dat)[1], samples=1:dim(dat)[2], manifests = ctm$manifestNames)
  }
  
  
  return(dat)
}

