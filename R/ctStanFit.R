T0VARredundancies <- function(ctm) {
  whichT0VAR_T0MEANSindvarying <- ctm$pars$matrix %in% 'T0VAR'  &  
    is.na(ctm$pars$value) &
    (ctm$pars$row %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying] |
        ctm$pars$col %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying])
  if(any(whichT0VAR_T0MEANSindvarying)){
    message('Free T0VAR parameters as well as indvarying T0MEANS -- fixing T0VAR pars to diag matrix 1e-3')
    ctm$pars$value[whichT0VAR_T0MEANSindvarying & ctm$pars$col == ctm$pars$row ] <- 1e-3
    ctm$pars$value[whichT0VAR_T0MEANSindvarying & ctm$pars$col != ctm$pars$row ] <- 0
    ctm$pars$param[whichT0VAR_T0MEANSindvarying] <- NA
    ctm$pars$transform[whichT0VAR_T0MEANSindvarying] <- NA
    ctm$pars$indvarying[whichT0VAR_T0MEANSindvarying] <- FALSE
  }
  return(ctm)
}

stansubjectdata <- function(ctsmodel, datalong,maxtimestep){
  #t0 index
  T0check<-rep(1,nrow(datalong))
  for(i in 2:nrow(datalong)){
    T0check[i]<- ifelse(datalong[i,ctsmodel$subjectIDname] != datalong[i-1,ctsmodel$subjectIDname], 1, 0)
  }
  if (!(ctsmodel$timeName %in% colnames(datalong))) stop(paste('time column', omxQuotes(ctsmodel$timeName), "not found in data"))
  if(any(is.na(datalong[,ctsmodel$timeName]))) stop('Missings in time column!')
  #check id and calculate intervals, discrete matrix indices
  # driftindex<-rep(0,nrow(datalong))
  # diffusionindex<-driftindex
  # cintindex<-driftindex
  oldsubi<-datalong[1,ctsmodel$subjectIDname]-1
  dT<-rep(-1,length(datalong[,ctsmodel$timeName]))
  
  for(rowi in 1:length(datalong[,ctsmodel$timeName])) {
    subi<-datalong[rowi,ctsmodel$subjectIDname]
    # if(rowi==1 && subi!=1) stop('subject id column must ascend from 1 to total subjects without gaps')
    # if(oldsubi!=subi && subi-oldsubi!=1) stop('subject id column must ascend from 1 to total subjects without gaps')
    if(subi - oldsubi == 1) {
      dT[rowi]<-0
      subistartrow<-rowi
    }
    if(subi - oldsubi == 0) {
      dT[rowi]<-datalong[rowi,ctsmodel$timeName] - datalong[rowi-1,ctsmodel$timeName]
      if(dT[rowi] < 0) stop(paste0('A time interval of ', dT[rowi],' was found at row ',rowi))
      if(dT[rowi] == 0) warning(paste0('A time interval of ', dT[rowi],' was found at row ',rowi))
      
    }
    oldsubi<-subi
  }
  
  if(mean(dT) > 3) message('Average time interval greater than 3 -- if using default priors, consider rescaling time data...')
  
  
  if(ctsmodel$n.TDpred > 0) {
    tdpreds <- datalong[,ctsmodel$TDpredNames,drop=FALSE]
    if(any(is.na(tdpreds))) {
      # if(NAtdpreds == 'error'
      message('Missingness in TDpreds! Replaced by zeroes...')
      tdpreds[is.na(tdpreds)] <-0 ## rough fix for missingness
    }
  }
  
  
  subdata <- list(
    Y=cbind(as.matrix(datalong[,ctsmodel$manifestNames])),
    subject=as.integer(datalong[,ctsmodel$subjectIDname]),
    time=datalong[,ctsmodel$timeName], #not used in model but used elsewhere
    ndatapoints=as.integer(nrow(datalong)),
    nobs_y=array(as.integer(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) length(x[x!=99999]))),dim=nrow(datalong)),
    whichobs_y=matrix(as.integer(t(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which(x!=99999))
      if(length(out)==0) out<-rep(0,ctsmodel$n.manifest)
      if(length(out)<ctsmodel$n.manifest) out<-c(out,rep(0,ctsmodel$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctsmodel$n.manifest)),
    nbinary_y=array(as.integer(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) 
      length(x[ctsmodel$manifesttype==1 & x!=99999]))),dim=nrow(datalong)),
    whichbinary_y=matrix(as.integer(t(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which(ctsmodel$manifesttype==1 & x!=99999)) #conditional on whichobs_y
      if(length(out)==0) out<-rep(0,ctsmodel$n.manifest)
      if(length(out)<ctsmodel$n.manifest) out<-c(out,rep(0,ctsmodel$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctsmodel$n.manifest)),
    ncont_y=array(as.integer(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) 
      length(x[(ctsmodel$manifesttype==0 | ctsmodel$manifesttype==2) & x!=99999]))),dim=nrow(datalong)),
    whichcont_y=matrix(as.integer(t(apply(datalong[,ctsmodel$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which( (ctsmodel$manifesttype==0 | ctsmodel$manifesttype==2) & x!=99999)) #conditional on whichobs_y
      if(length(out)==0) out<-rep(0,ctsmodel$n.manifest)
      if(length(out)<ctsmodel$n.manifest) out<-c(out,rep(0,ctsmodel$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctsmodel$n.manifest))
  )
  
  if(ctsmodel$n.TDpred ==0) tdpreds <- matrix(0,subdata$ndatapoints,0) #subdata$ndatapoints,
  subdata$tdpreds=array(as.matrix(tdpreds),dim=c(nrow(tdpreds),ncol(tdpreds)))
  
  #subset selection
  if(is.null(ctsmodel$dokalmanrows)) subdata$dokalmanrows <- 
    rep(1L, subdata$ndatapoints) else subdata$dokalmanrows <- as.integer(ctsmodel$dokalmanrows)
  subdata$dokalmanpriormodifier = sum(subdata$dokalmanrows)/subdata$ndatapoints
  
  return(subdata)
}

verbosify<-function(sf,verbose=2){
  sm <- sf$stanmodel
  sd <- sf$standata
  sd$verbose=as.integer(verbose)
  
  sfr <- stan_reinitsf(sm,sd)
  log_prob(sfr,sf$stanfit$rawest)
}

#' ctStanFit
#'
#' Fits a ctsem model specified via \code{\link{ctModel}} with type either 'stanct' or 'standt', using Bayseian inference software
#' Stan. 
#' 
#' @param datalong long format data containing columns for subject id (numeric values, 1 to max subjects), manifest variables, 
#' any time dependent (i.e. varying within subject) predictors, 
#' and any time independent (not varying within subject) predictors.
#' @param ctstanmodel model object as generated by \code{\link{ctModel}} with type='stanct' or 'standt', for continuous or discrete time
#' models respectively.
#' @param stanmodeltext already specified Stan model character string, generally leave NA unless modifying Stan model directly.
#' (Possible after modification of output from fit=FALSE)
#' @param intoverstates logical indicating whether or not to integrate over latent states using a Kalman filter. 
#' Generally recommended to set TRUE unless using non-gaussian measurement model. 
#' @param binomial Deprecated. Logical indicating the use of binary rather than Gaussian data, as with IRT analyses.
#' This now sets \code{intoverstates = FALSE} and the \code{manifesttype} of every indicator to 1, for binary.
#' @param fit If TRUE, fit specified model using Stan, if FALSE, return stan model object without fitting.
#' @param intoverpop if TRUE, integrates over population distribution of parameters rather than full sampling.
#' Allows for optimization of non-linearities and random effects.
#' @param plot if TRUE, for sampling, a Shiny program is launched upon fitting to interactively plot samples. 
#' May struggle with many (e.g., > 5000) parameters. For optimizing, various optimization details are plotted -- in development.
#' @param derrind vector of integers denoting which latent variables are involved in dynamic error calculations.
#' latents involved only in deterministic trends or input effects can be removed from matrices (ie, that
#' obtain no additional stochastic inputs after first observation), speeding up calculations. 
#' If unsure, leave default of 'all' ! Ignored if intoverstates=FALSE.
#' @param optimize if TRUE, use \code{\link{stanoptimis}} function for maximum a posteriori / importance sampling estimates, 
#' otherwise use the HMC sampler from Stan, which is (much) slower, but generally more robust, accurate, and informative.
#' @param optimcontrol list of parameters sent to \code{\link{stanoptimis}} governing optimization / importance sampling.
#' @param nopriors logical. If TRUE, any priors are disabled -- sometimes desirable for optimization. 
#' @param iter number of iterations, half of which will be devoted to warmup by default when sampling.
#' When optimizing, this is the maximum number of iterations to allow -- convergence hopefully occurs before this!
#' @param inits vector of parameter start values, as returned by the rstan function \code{rstan::unconstrain_pars} for instance. 
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling. Unless the cores
#' argument is also set, the number of chains determines the number of cpu cores used, up to 
#' the maximum available minus one. Irrelevant when \code{optimize=TRUE}.
#' @param cores number of cpu cores to use. Either 'maxneeded' to use as many as available minus one,
#' up to the number of chains, or a positive integer. If \code{optimize=TRUE}, more cores are generally faster.
#' @param control List of arguments sent to \code{\link[rstan]{stan}} control argument, 
#' regarding warmup / sampling behaviour. Unless specified, values used are:
#' list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001)
#' @param nlcontrol List of non-linear control parameters. 
#' \code{nldynamics} defaults to "auto", but may also be a logical. Set to FALSE to use estimator that assumes linear dynamics, 
#' TRUE to use non-linear estimator. "auto" selects linear when the model is obviously linear, 
#' otherwise nonlinear -- nonlinear is slower.
#' \code{nlmeasurement} defaults to "auto", but may also be a logical. Set to TRUE to use non linear measurement model estimator, 
#' FALSE to use linear model. "auto" selects linear if appropriate, otherwise nonlinear. Non-linear methods are slower but applicable to both linear
#' and non linear cases.
#' \code{maxtimestep} must be a positive numeric,  specifying the largest time
#' span covered by the numerical integration. The large default ensures that for each observation time interval, 
#' only a single step of exponential integration is used. When \code{maxtimestep} is smaller than the observation time interval, 
#' the integration is nested within an Euler like loop. 
#' Smaller values may offer greater accuracy, but are slower and not always necessary. Given the exponential integration,
#' linear model elements are fit exactly with only a single step. 
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param stationary Logical. If TRUE, T0VAR and T0MEANS input matrices are ignored, 
#' the parameters are instead fixed to long run expectations. More control over this can be achieved
#' by instead setting parameter names of T0MEANS and T0VAR matrices in the input model to 'stationary', for
#' elements that should be fixed to stationarity.
#' @param forcerecompile logical. For development purposes. 
#' If TRUE, stan model is recompiled, regardless of apparent need for compilation.
#' @param savescores Logical. If TRUE, output from the Kalman filter is saved in output. For datasets with many variables
#' or time points, will increase file size substantially.
#' @param savesubjectmatrices Logical. If TRUE, subject specific matrices are saved -- only relevant when either time dependent predictors
#' are used, or individual differences are obtained via sampling (not via optimization, where they are integrated over).
#' @param gendata Logical -- If TRUE, uses provided data for only covariates and a time and missingness structure, and 
#' generates random data according to the specified model / priors. 
#' Generated data is in the $Ygen subobject after running \code{extract} on the fit object.
#' For datasets with many manifest variables or time points, file size may be large.
#' To generate data based on the posterior of a fitted model, see \code{\link{ctStanGenerateFromFit}}.
#' @param ... additional arguments to pass to \code{\link[rstan]{stan}} function.
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' \donttest{
#' #test data with 2 manifest indicators measuring 1 latent process each, 
#' # 1 time dependent predictor, 3 time independent predictors
#' head(ctstantestdat) 
#' 
#' #generate a ctStanModel
#' model<-ctModel(type='stanct',
#' n.latent=2, latentNames=c('eta1','eta2'),
#' n.manifest=2, manifestNames=c('Y1','Y2'),
#' n.TDpred=1, TDpredNames='TD1', 
#' n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
#' LAMBDA=diag(2))
#' 
#' #set all parameters except manifest means to be fixed across subjects
#' model$pars$indvarying[-c(19,20)] <- FALSE
#' 
#' #fit model to data (takes a few minutes - but insufficient 
#' # iterations and max_treedepth for inference!)
#' fit<-ctStanFit(ctstantestdat, model, iter=200, chains=2, 
#' control=list(max_treedepth=6))
#' 
#' #output functions
#' summary(fit) 
#' 
#' plot(fit,wait=FALSE)
#' 
#' }
#' \dontrun{
#' library(ctsem)
#' set.seed(3)
#' 
#' #  Data generation (run this, but no need to understand!) -----------------
#' 
#' Tpoints <- 20
#' nmanifest <- 4
#' nlatent <- 2
#' nsubjects<-20
#' 
#' #random effects
#' age <- rnorm(nsubjects) #standardised
#' cint1<-rnorm(nsubjects,2,.3)+age*.5
#' cint2 <- cint1*.5+rnorm(nsubjects,1,.2)+age*.5
#' tdpredeffect <- rnorm(nsubjects,5,.3)+age*.5
#' 
#' for(i in 1:nsubjects){
#'   #generating model
#'   gm<-ctModel(Tpoints=Tpoints,n.manifest = nmanifest,n.latent = nlatent,n.TDpred = 1,
#'     LAMBDA = matrix(c(1,0,0,0, 0,1,.8,1.3),nrow=nmanifest,ncol=nlatent),
#'     DRIFT=matrix(c(-.3, .2, 0, -.5),nlatent,nlatent),
#'     TDPREDEFFECT=matrix(c(tdpredeffect[i],0),nrow=nlatent),
#'     TDPREDMEANS=matrix(c(rep(0,Tpoints-10),1,rep(0,9)),ncol=1),
#'     DIFFUSION = matrix(c(1, 0, 0, .5),2,2),
#'     CINT = matrix(c(cint1[i],cint2[i]),ncol=1),
#'     T0VAR=diag(2,nlatent,nlatent),
#'     MANIFESTVAR = diag(.5, nmanifest))
#' 
#'   #generate data
#'   newdat <- ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 2,
#'     dtmat<-rbind(c(rep(.5,8),3,rep(.5,Tpoints-9))),
#'     wide = FALSE)
#'   newdat[,'id'] <- i #set id for each subject
#'   newdat <- cbind(newdat,age[i]) #include time independent predictor
#'   if(i ==1) {
#'     dat <- newdat[1:(Tpoints-10),] #pre intervention data
#'     dat2 <- newdat #including post intervention data
#'   }
#'   if(i > 1) {
#'     dat <- rbind(dat, newdat[1:(Tpoints-10),])
#'     dat2 <- rbind(dat2,newdat)
#'   }
#' }
#' colnames(dat)[ncol(dat)] <- 'age'
#' colnames(dat2)[ncol(dat)] <- 'age'
#' 
#' 
#' #plot generated data for sanity
#' plot(age)
#' matplot(dat[,gm$manifestNames],type='l',pch=1)
#' plotvar <- 'Y1'
#' plot(dat[dat[,'id']==1,'time'],dat[dat[,'id']==1,plotvar],type='l',
#'   ylim=range(dat[,plotvar],na.rm=TRUE))
#' for(i in 2:nsubjects){
#'   points(dat[dat[,'id']==i,'time'],dat[dat[,'id']==i,plotvar],type='l',col=i)
#' }
#' 
#' 
#' dat2[,gm$manifestNames][sample(1:length(dat2[,gm$manifestNames]),size = 100)] <- NA
#' 
#' 
#' #data structure
#' head(dat2)
#' 
#' 
#' # Model fitting -----------------------------------------------------------
#' 
#' ##simple univariate default model
#' 
#' m <- ctModel(type = 'stanct', manifestNames = c('Y1'), LAMBDA = diag(1))
#' ctModelLatex(m)
#' 
#' #Specify univariate linear growth curve
#' 
#' m1 <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   DRIFT=matrix(-.0001,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#' 
#' ctModelLatex(m1)
#' 
#' #fit
#' f1 <- ctStanFit(datalong = dat2, ctstanmodel = m1, optimize=TRUE, nopriors=TRUE)
#' 
#' summary(f1)
#' 
#' #plots of individual subject models v data
#' ctKalman(f1,plot=TRUE,subjects=1,kalmanvec=c('y','yprior'),timestep=.01)
#' ctKalman(f1,plot=TRUE,subjects=1:3,kalmanvec=c('y','ysmooth'),timestep=.01,errorvec=NA)
#' 
#' ctStanPostPredict(f1, wait=FALSE) #compare randomly generated data from posterior to observed data
#' 
#' cf<-ctCheckFit(f1) #compare mean and covariance of randomly generated data to observed cov
#' plot(cf,wait=FALSE)
#' 
#' 
#' 
#' #Include intervention
#' m2 <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#' 
#' f2 <- ctStanFit(datalong = dat2, ctstanmodel = m2, optimize=TRUE)
#' 
#' summary(f2)
#' 
#' ctKalman(f2,plot=TRUE,subjects=1,kalmanvec=c('y','ysmooth'))
#' ctKalman(f2,plot=TRUE,subjects=1:3,kalmanvec=c('y','ysmooth'),errorvec=NA,legend=FALSE)
#' 
#' ctStanPostPredict(f2, datarows=1:100, wait=FALSE)
#' 
#' 
#' 
#' #Individual differences in intervention, Bayesian estimation, covariates
#' m2i <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   TIpredNames = 'age',
#'   TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect||TRUE'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#'   
#' plot(m2i)
#' 
#' f2i <- ctStanFit(datalong = dat2, ctstanmodel = m2i,intoverpop=TRUE,
#'   iter=200,chains=2,control=list(max_treedepth=7))
#' summary(f2i)
#' ctStanPlotPost(f2i)
#' ctKalman(f2i,kalmanvec=c('y','ysmooth'),subjects=2:4,plot=TRUE,errorvec=NA)
#' 
#' 
#' #Including covariate effects
#' m2ic <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TIpred = 1, TIpredNames = 'age',
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#' 
#' m2ic$pars$indvarying[m2ic$pars$matrix %in% 'TDPREDEFFECT'] <- TRUE
#' 
#' plot(m2ic)
#' 
#' f2ic <- ctStanFit(datalong = dat2, ctstanmodel = m2ic,optimize=TRUE)
#' summary(f2ic)
#' 
#' ctStanTIpredeffects(fit = f2ic,includeMeanUncertainty = TRUE,whichpars = 'TDPREDEFFECT',
#'   plot=TRUE,probs = c(.025,.5,.975))
#' 
#' #Include deterministic dynamics
#' m3 <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
#' 
#' ctModelLatex(m3)
#' 
#' f3 <- ctStanFit(datalong = dat2, ctstanmodel = m3, optimize=TRUE)
#' 
#' summary(f3)
#' 
#' ctKalman(f3,plot=TRUE,subjects=1,kalmanvec=c('y','ysmooth'))
#' ctKalman(f3,plot=TRUE,subjects=1:3,kalmanvec=c('y','ysmooth'),errorvec=NA)
#' 
#' 
#' 
#' 
#' 
#' #Add system noise to allow for fluctuations that persist in time
#' m3n <- ctModel(type = 'stanct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix('diffusion',nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c(0),nrow=1,ncol=1))
#' 
#' ctModelLatex(m3n)
#' 
#' f3n <- ctStanFit(datalong = dat2, ctstanmodel = m3n, optimize=TRUE)
#' 
#' summary(f3n)
#' 
#' k=ctKalman(f3n,plot=T,subjects=1,kalmanvec=c('y','etasmooth'),timestep=.01)
#' ctKalman(f3n,plot=TRUE,subjects=1:3,kalmanvec=c('y','etasmooth'),errorvec=NA)
#' 
#' 
#' 
#' 
#' 
#' #include 2nd latent process
#' 
#' m4 <- ctModel(n.manifest = 2,n.latent = 2, type = 'stanct',
#'   manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),
#'   n.TDpred=1,TDpredNames = 'TD1',
#'   TDPREDEFFECT=matrix(c('tdpredeffect1','tdpredeffect2'),nrow=2,ncol=1),
#'   DRIFT=matrix(c('drift11','drift21','drift12','drift22'),nrow=2,ncol=2),
#'   DIFFUSION=matrix(c('diffusion11','diffusion21',0,'diffusion22'),nrow=2,ncol=2),
#'   CINT=matrix(c('cint1','cint2'),nrow=2,ncol=1),
#'   T0MEANS=matrix(c('t0m1','t0m2'),nrow=2,ncol=1),
#'   T0VAR=matrix(c('t0var11','t0var21',0,'t0var22'),nrow=2,ncol=2),
#'   LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2),
#'   MANIFESTMEANS=matrix(c(0,0),nrow=2,ncol=1),
#'   MANIFESTVAR=matrix(c('merror1',0,0,'merror2'),nrow=2,ncol=2))
#' 
#' f4 <- ctStanFit(datalong = dat2, ctstanmodel = m4,optimize=TRUE)
#' 
#' summary(f4)
#' 
#' ctStanDiscretePars(f4,plot=TRUE) #auto and cross regressive plots over time
#' 
#' ctKalman(f4,plot=TRUE,subjects=1,kalmanvec=c('y','ysmooth'))
#' ctKalman(f4,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'),errorvec=NA)
#' 
#' 
#' 
#' #dynamic factor model -- fixing CINT to 0 and freeing indicator level intercepts
#' 
#' m3df <- ctModel(type = 'stanct',
#'   manifestNames = c('Y2','Y3'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix('diffusion',nrow=1,ncol=1),
#'   CINT=matrix(c(0),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = matrix(c(1,'Y3loading'),nrow=2,ncol=1),
#'   MANIFESTMEANS=matrix(c('Y2_int','Y3_int'),nrow=2,ncol=1),
#'   MANIFESTVAR=matrix(c('Y2residual',0,0,'Y3residual'),nrow=2,ncol=2))
#' 
#' ctModelLatex(m3df)
#' 
#' f3df <- ctStanFit(datalong = dat2, ctstanmodel = m3df, optimize=TRUE)
#' 
#' summary(f3df)
#' 
#' ctKalman(f3df,plot=TRUE,subjects=1,kalmanvec=c('y','ysmooth'),errorvec=NA)
#' ctKalman(f3df,plot=TRUE,subjects=1:3,kalmanvec=c('y','ysmooth'),errorvec=NA)
#' 
#' 
#' 
#' }

ctStanFit<-function(datalong, ctstanmodel, stanmodeltext=NA, iter=1000, intoverstates=TRUE, binomial=FALSE,
  fit=TRUE, intoverpop=FALSE, stationary=FALSE,plot=FALSE,  derrind='all',
  optimize=FALSE,  optimcontrol=list(),
  nlcontrol = list(), nopriors=FALSE, chains=2,
  cores=ifelse(optimize,getOption("mc.cores", 2L),'maxneeded'),
  inits=NULL,
  forcerecompile=FALSE,savescores=FALSE,savesubjectmatrices=TRUE,
  gendata=FALSE,
  control=list(),verbose=0,...){
  if(.Machine$sizeof.pointer == 4) message('Bayesian functions not available on 32 bit systems') else{
    if(!'ctStanModel' %in% class(ctstanmodel)) stop('not a ctStanModel object')
    
    if(is.null(nlcontrol$ukfspread)) nlcontrol$ukfspread = 1e-1
    if(is.null(nlcontrol$maxtimestep)) nlcontrol$maxtimestep = 999999
    if(is.null(nlcontrol$nldynamics)) nlcontrol$nldynamics = 'auto'
    if(is.null(nlcontrol$nlmeasurement)) nlcontrol$nlmeasurement = 'auto'
    if(is.null(nlcontrol$Jstep)) nlcontrol$Jstep = 1e-6
    nldynamics <- nlcontrol$nldynamics
    
    args=match.call()
    
    ctm <- ctstanmodel
    
    idName<-ctm$subjectIDname
    timeName<-ctm$timeName
    if(!ctm$timeName %in% colnames(datalong) && !ctm$continuoustime) {
      datalong[[ctm$timeName]] <- 1:nrow(datalong)
    }
    continuoustime<-ctm$continuoustime
    ctstanmodelbase <- ctm
    ctm <- ctModel0DRIFT(ctm, continuoustime)
    ctm <- ctModelStatesAndPARS(ctm)
    # 
    ctm <- ctModelTransformsToNum(ctm)
    ctm$pars <- ctStanModelCleanctspec(ctm$pars)
    # ctm <- ctStanModelIntOverPop(ctm)
    
    if('data.table' %in% class(datalong)) datalong <- data.frame(datalong)
    
    ###stationarity
    # if(intoverpop && (any(ctm$pars$param %in% 'stationary') || stationary)) stop('stationarity incompatible with random effects integration')
    if(stationary) {
      ctm$pars$param[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- 'stationary'
      ctm$pars$value[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- NA
      ctm$pars$indvarying[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- FALSE
    }
    
    #collect individual stationary elements and update ctm$pars
    t0varstationary <- as.matrix(rbind(ctm$pars[which(ctm$pars$param %in% 'stationary' & ctm$pars$matrix %in% 'T0VAR'),c('row','col')]))
    if(nrow(t0varstationary) > 0){ #ensure upper tri is consistent with lower
      for(i in 1:nrow(t0varstationary)){
        if(t0varstationary[i,1] != t0varstationary[i,2]) t0varstationary <- rbind(t0varstationary,t0varstationary[i,c(2,1)])
      }}
    t0varstationary = unique(t0varstationary) #remove any duplicated rows
    
    t0meansstationary <- as.matrix(rbind(ctm$pars[which(ctm$pars$param[ctm$pars$matrix %in% 'T0MEANS'] %in% 'stationary'),c('row','col')]))
    ctm$pars$value[ctm$pars$param %in% 'stationary'] <- -99 #does this get inserted?
    ctm$pars$indvarying[ctm$pars$param %in% 'stationary'] <- FALSE
    ctm$pars$transform[ctm$pars$param %in% 'stationary'] <- NA
    ctm$pars$param[ctm$pars$param %in% 'stationary'] <- NA
    
    nt0varstationary <- nrow(t0varstationary)
    nt0meansstationary <- nrow(t0meansstationary)
    
    ####
    
    
    if(length(unique(datalong[,idName]))==1 && any(ctm$pars$indvarying[is.na(ctm$pars$value)]==TRUE) && 
        is.null(ctm$fixedrawpopmeans) && is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
      ctm$pars$indvarying <- FALSE
      message('Individual variation not possible as only 1 subject! indvarying set to FALSE on all parameters')
    }
    
    
    if(length(unique(datalong[,idName]))==1 & any(!is.na(ctm$pars$value[ctm$pars$matrix %in% 'T0VAR'])) & 
        is.null(ctm$fixedrawpopmeans) & is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
      for(ri in 1:nrow(ctm$pars)){
        if(is.na(ctm$pars$value[ri]) && ctm$pars$matrix[ri] %in% 'T0VAR'){
          ctm$pars$value[ri] <- ifelse(ctm$pars$row[ri] == ctm$pars$col[ri], 1, 0)
        }
      }
      message('Free T0VAR parameters fixed to diagonal matrix of 1\'s as only 1 subject - consider appropriateness!')
    }
    
    if(binomial){
      message('Binomial argument deprecated -- in future set manifesttype in the model object to 1 for binary indicators')
      intoverstates <- FALSE
      ctm$manifesttype[] <- 1
    }
    
    recompile <- FALSE
    if(optimize && !intoverstates) stop('intoverstates=TRUE required for optimization!')
    
    if(optimize && !intoverpop && any(ctm$pars$indvarying[is.na(ctm$pars$value)]) && 
        is.null(ctm$fixedrawpopchol) && is.null(ctm$fixedsubpars)){
      intoverpop <- TRUE
      message('Setting intoverpop=TRUE to enable optimization of random effects...')
    }
    
    # if(intoverpop & any(ctm$pars$indvarying[is.na(ctm$pars$value)]) & nldynamics=='auto') {
    #   nldynamics=TRUE 
    # } else 
    if(intoverpop==TRUE && !any(ctm$pars$indvarying[is.na(ctm$pars$value)])) {
      message('No individual variation -- disabling intoverpop switch'); intoverpop <- FALSE
    } 
    # else 
    # if(intoverpop & any(ctm$pars$indvarying[is.na(ctm$pars$value)]) & nldynamics==FALSE) { stop('nldynamics cannot be set FALSE if intoverpop is TRUE')
    # }
    
    if(intoverpop)   ctm <- ctStanModelIntOverPop(ctm)
    
    ctm$jacobian <- ctJacobian(ctm)
    
    jl <- ctModelUnlist(ctm$jacobian,names(ctm$jacobian))
    jl2 <- as.data.frame(rbind(data.table(ctm$pars[1,]),data.table(jl),fill=TRUE))[-1,]
    
    jl3=ctModelTransformsToNum(list(pars=data.frame(jl2)))
    jl3$pars$indvarying<-FALSE
    ctm$pars <- rbind(ctm$pars, jl3$pars)
    
    
    if(naf(!is.na(ctm$rawpopsdbaselowerbound))) recompile <- TRUE
    if(ctm$rawpopsdbase != 'normal(0,1)') recompile <- TRUE
    if(ctm$rawpopsdtransform != 'log1p(exp(2*rawpopsdbase-1)) .* sdscale') recompile <- TRUE
    
    
    if(cores=='maxneeded') cores=max(1,min(c(chains,parallel::detectCores()-1))) else cores <-max(1, min(cores,parallel::detectCores()-1))
    
    
    ctm <- T0VARredundancies(ctm)
    
    
    mats <- ctStanMatricesList()
    
    
    if(!all(ctm$pars$transform[!is.na(suppressWarnings(as.integer(ctm$pars$transform)))] %in% c(0,1,2,3,4))) stop('Unknown transform specified -- integers should be 0 to 4')
    
    # if(binomial) {
    #   ctm$pars<-ctm$pars[ctm$pars$matrix != 'MANIFESTVAR',]
    #   message(paste0('MANIFESTVAR matrix is ignored when binomial=TRUE'))
    # }
    
    
    
    manifesttype=ctm$manifesttype
    
    #fix binary manifestvariance
    if(any(manifesttype==1)){ #if any non continuous variables, (with free parameters)...
      if(any(is.na(as.numeric(c(ctm$pars$value[ctm$pars$matrix=='MANIFESTVAR'][ctm$pars$row[ctm$pars$matrix=='MANIFESTVAR'] %in% which(manifesttype==1)],
        ctm$pars$value[ctm$pars$matrix=='MANIFESTVAR'][ctm$pars$col[ctm$pars$matrix=='MANIFESTVAR'] %in% which(manifesttype!=0)]))))){
        message('Fixing any free MANIFESTVAR parameters for binary indicators to deterministic calculation')
        ctm$pars$value[ctm$pars$matrix=='MANIFESTVAR'][ctm$pars$row[ctm$pars$matrix=='MANIFESTVAR'] %in% which(manifesttype==1)] <- 0
        ctm$pars$value[ctm$pars$matrix=='MANIFESTVAR'][ctm$pars$col[ctm$pars$matrix=='MANIFESTVAR'] %in% which(manifesttype==1)] <- 0
        ctm$pars$value[ctm$pars$matrix=='MANIFESTVAR' & ctm$pars$row %in% which(manifesttype==1) & ctm$pars$row == ctm$pars$col] <- 1e-5
      }}
    ctm$pars <- ctStanModelCleanctspec(ctm$pars)
    
    #adjust transforms for optimization
    if(1==99 && optimize) {
      message('Adapting standard deviation transforms for optimization')
      # ctm$rawpopsdtransform <- 'log(1+exp(rawpopsdbase))*10'
      indices <- ctm$pars$matrix %in% c('DIFFUSION','MANIFESTVAR','T0VAR') & is.na(ctm$pars$value) & ctm$pars$row == ctm$pars$col
      ctm$pars$transform[indices] <- 1
      ctm$pars$multiplier[indices] <- 1
      ctm$pars$meanscale[indices] <- 10
      ctm$pars$offset[indices] <- 0
    }
    
    
    
    # if(nt0meansstationary ==0) t0meansstationary <- matrix(-99,ncol=2)
    
    
    nsubjects <- length(unique(datalong[, idName])) 
    
    #create random effects indices for each matrix
    for(mati in names(mats$base)){
      if( (!intoverpop && any(ctm$pars$indvarying[ctm$pars$matrix==mati])) || 
          (ctm$n.TIpred >0 && any(unlist(ctm$pars[ctm$pars$matrix==mati,paste0(ctm$TIpredNames,'_effect')])))) subindex <- 1 else subindex <- 0
          assign(paste0(mati,'subindex'), subindex)
    }
    if(stationary || nt0varstationary > 0) T0VARsubindex <- max(c(T0VARsubindex,DRIFTsubindex,DIFFUSIONsubindex))
    if(stationary || nt0meansstationary > 0) T0MEANSsubindex <- max(c(T0MEANSsubindex,DRIFTsubindex,CINTsubindex))
    asymCINTsubindex <- max(c(CINTsubindex,DRIFTsubindex))
    asymDIFFUSIONsubindex <- max(c(DIFFUSIONsubindex,DRIFTsubindex))
    DIFFUSIONcovsubindex <- max(c(DIFFUSIONsubindex))
    
    #simply exponential?
    driftdiagonly <- ifelse(all(!is.na(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col]) &
        all(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col] == 0) ), 1, 0)
    
    n.latent<-ctm$n.latent
    n.manifest<-ctm$n.manifest
    n.TDpred<-ctm$n.TDpred
    n.TIpred<-ctm$n.TIpred
    
    ###data checks
    
    if(any(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))) stop(paste0('
      variables: ', paste0(c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames)[
        which(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))], ', '),' not in data'))
    
    if (!(idName %in% colnames(datalong))) stop(paste('id column', omxQuotes(idName), "not found in data"))
    
    
    if(nopriors==FALSE){
      if(ctm$n.TIpred > 1 && any(abs(colMeans(datalong[,c(ctm$TIpredNames),drop=FALSE],na.rm=TRUE)) > .3)){
        message('Uncentered TI predictors noted -- interpretability may be hindered and default priors may not be appropriate')
      }
      
      #scale check
      if(naf(any(abs(colMeans(datalong[,c(ctm$manifestNames,ctm$TDpredNames),drop=FALSE],na.rm=TRUE)) > 5))){
        message('Uncentered data noted -- default priors *may* not be appropriate')
      }
      
      
      if(naf(any(abs(apply(datalong[,c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames),drop=FALSE],2,sd,na.rm=TRUE)) > 3))){
        message('Unscaled data noted -- default priors may not be appropriate')
      }
    }
    
    #fit spec checks
    # if(binomial & any(intoverstates)) stop('Binomial only possible with intoverstates=FALSE')
    
    #id mapping
    
    original <- unique(datalong[,idName])
    datalong <- makeNumericIDs(datalong,idName,timeName)
    new <- unique(datalong[,idName])
    idmap <- cbind(original, new)
    
    
    
    
    
    #generate model matrix lists for stan
    for(i in 1:nrow(ctm$pars)){ #simplify any calcs
      if(grepl('[',ctm$pars$param[i],fixed=TRUE)) ctm$pars$param[i] <- Deriv::Simplify(ctm$pars$param[i])
    }
    ctsmodelmats <- ctStanModelMatrices(ctm)
    matsetup <- ctsmodelmats$matsetup
    matvalues <- ctsmodelmats$matvalues
    extratforms <- ctsmodelmats$extratforms
    TIPREDEFFECTsetup=ctsmodelmats$TIPREDEFFECTsetup
    matrixdims <- ctsmodelmats$matrixdims
    ctm$calcs <- ctsmodelmats$calcs
    
    
    #get extra calculations and adjust model spec as needed
    ctm <- ctStanCalcsList(ctm)
    
    if(sum(sapply(ctm$calcs[!names(ctm$calcs) %in% c('jacobian','measurement')],length)) > 0 || 
        any(matsetup$when %in% c(1,2,3)) ||
        length(ctm$calcs$jacobian) - sum(grepl('sJy[',unlist(ctm$calcs$jacobian),fixed=TRUE)) > 0 #non measurement jacobians
        ){
      if(nldynamics == FALSE) warning('Linear model requested but nonlinear model specified! May be a poor approximation') else nldynamics <- TRUE 
    }
    
    ncalcs <- length(unlist(ctm$calcs)) 
    ncalcsNoJ<- length(unlist(ctm$calcs)[!grepl('JAx[',unlist(ctm$calcs),fixed=TRUE)])
    if(ncalcsNoJ > 0) recompile <- TRUE
    if(!recompile && ncalcsNoJ >0)  message('Finite difference jacobian used to avoid recompiling -- use forcerecompile=TRUE for analytic jacobians')
    if( (nt0varstationary + nt0meansstationary) >0 && 
        length(c(ctm$calcs$driftcint, ctm$calcs$diffusion)) > 0) message('Stationarity assumptions based on initial states when using non-linear dynamics')
    
    nlmeasurement <- nlcontrol$nlmeasurement

    if(length(ctm$calcs$measurement) > 0 || 
        any(matsetup$when %in% c(4)) ||
        (intoverpop && any(ctstanmodelbase$pars$indvarying[ctstanmodelbase$pars$matrix %in% names(mats$measurement)]))) { 
      if(nlmeasurement == FALSE) warning('Linear measurement model requested but nonlinear measurement specified!') else nlmeasurement <- TRUE
      
    }
    
    if(nldynamics==TRUE && !intoverstates) stop('intoverstates must be TRUE for nonlinear dynamics')
    
    if(nlmeasurement=='auto') nlmeasurement <- FALSE
    if(nlmeasurement) message('Using nonlinear Kalman filter for measurement update');
    if(!nlmeasurement) message('Using linear Kalman filter for measurement update');
    
    
    if(nldynamics == 'auto') nldynamics <- FALSE
    if(nldynamics) message('Using nonlinear Kalman filter for dynamics')
    if(!nldynamics) message('Using linear Kalman filter for dynamics')
    
    if(intoverstates==FALSE || all(derrind=='all') ) derrind = 1:n.latent
    # if(all(derrind=='all')) derrind = sort(unique(ctm$pars$col[
    #   ctm$pars$matrix=='DIFFUSION' & (!is.na(ctm$pars$param) | ctm$pars$value!=0)]))
    derrind = as.integer(derrind)
    if(any(derrind > n.latent)) stop('derrind > than n.latent found!')
    if(length(derrind) > n.latent) stop('derrind vector cannot be longer than n.latent!')
    if(length(unique(derrind)) < length(derrind)) stop('derrind vector cannot contain duplicates or!')
    ndiffusion=length(derrind)
    # message(paste(ndiffusion ,'/',n.latent,'latent variables needed for covariance calculations'))
    
    
    
    
    nindvarying <- max(matsetup$indvarying)
    nparams <- max(matsetup$param[matsetup$when==0])
    nmatrices <- length(mats$base)
    matsetup[which(matsetup$indvarying > 0),]
    indvaryingindex <- matsetup$param[which(matsetup$indvarying > 0)]
    indvaryingindex <- array(indvaryingindex[!duplicated(indvaryingindex)])
    
    sdscale <- array(matvalues$sdscale[match(indvaryingindex,matsetup$param)])
    
    
    if(any(matsetup[,'transform'] < -10)) recompile <- TRUE #if custom transforms needed
    
    
    
    
    
    if(is.na(stanmodeltext)) {
      stanmodeltext<- ctStanModelWriter(ctm, gendata, extratforms,matsetup)
    }  else recompile<-TRUE
    
    # message('using ols!!!!!!!!!')
    #  stanmodeltext<- ctolsStanModelWriter(ctm, gendata, extratforms)
    # recompile=TRUE
    
    
    # out<-list(stanmodeltext=stanmodeltext)
    
    #tipred data
    if(ctm$n.TIpred > 0) {
      tipreds <- datalong[match(unique(datalong[,ctm$subjectIDname]),datalong[,ctm$subjectIDname]),ctm$TIpredNames,drop=FALSE]
      if(any(is.na(tipreds))) {
        if(!optimize){
          message(paste0('Missingness in TIpreds - sampling ', sum(is.na(tipreds)),' values'))
          tipreds[is.na(tipreds)] = 99999
        }
        if(optimize){
          message(paste0('Missingness in TIpreds - single imputing ', sum(is.na(tipreds)),'  NA\'s to allow optimization -- TI predictor effect estimates will be overly confident.'))
          tipreds[is.na(tipreds)] = 0
          
          meandat <- data.table((datalong))[ , lapply(.SD, function(x) 
            mean(x,na.rm=TRUE)) , 
            by=c("id")]
          sddat <- data.table((datalong))[ , lapply(.SD, function(x) 
            sd(x,na.rm=TRUE)) , 
            by=c("id")]
          sddat<-sddat[,!colnames(sddat) %in% ctm$subjectIDname,with=FALSE]
          meandat <- meandat[,apply(meandat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]
          sddat <- sddat[,apply(sddat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]
          colnames(sddat) <- paste0('sd_',colnames(sddat))
          meandat <- cbind(meandat,sddat)
          
          for(i in 1:n.TIpred){
            lmform = formula(paste0(ctm$TIpredNames[i],' ~ 1 + ',
              paste0(colnames(meandat)[-which(colnames(meandat) %in% ctm$TIpredNames[i])],
                collapse=' + ')))
            # lmr <- MASS::lm.ridge(formula = lmform, data = meandat,lambda=0,na.action=na.exclude)
            lmf <- lm(formula = lmform,data = meandat,na.action=na.exclude)
            # lmf$coefficients[-1] <- lmr$coef
            # lmf$coefficients[1] <- lmr$ym
            # lmr$coef * 
            #   apply(meandat[,-which(colnames(meandat) %in% ctm$TIpredNames[i]),with=FALSE],
            #     2,sd,na.rm=TRUE) - 
            #   lmf$coefficients[-1]
            # 
            plot(c(meandat[,ctm$TIpredNames[i],with=FALSE])[[1]],predict(lmf),main=ctm$TIpredNames[i])
            tipreds[is.na(tipreds[,1]),1] <- predict(lmf)[is.na(tipreds[,1])]
          }
        }
      }
    }
    
    datalong[,c(ctm$manifestNames,ctm$TIpredNames)][is.na(datalong[,c(ctm$manifestNames,ctm$TIpredNames)])]<-99999 #missing data
    
    
    
    standata<-c(stansubjectdata(ctsmodel = ctm,datalong = datalong, maxtimestep = nlcontrol$maxtimestep), 
      list(
        nsubjects=as.integer(nsubjects),
        nmanifest=as.integer(n.manifest),
        nlatentpop = as.integer(ifelse(intoverpop ==1, max(ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS']),  n.latent)),
        nldynamics=as.integer(nldynamics),
        Jstep = nlcontrol$Jstep,
        maxtimestep = nlcontrol$maxtimestep,
        dokalman=as.integer(is.null(ctm$dokalman)),
        intoverstates=as.integer(intoverstates),
        verbose=as.integer(verbose),
        manifesttype=array(as.integer(manifesttype),dim=length(manifesttype)),
        indvaryingindex=array(as.integer(indvaryingindex)),
        intoverpopindvaryingindex=array(as.integer(ctm$intoverpopindvaryingindex)),
        notindvaryingindex=array(as.integer(which(!(1:nparams) %in% indvaryingindex))),
        continuoustime=as.integer(sum(continuoustime)),
        nlatent=as.integer(n.latent),
        ntipred=as.integer(n.TIpred),
        ntdpred=as.integer(n.TDpred),
        binomial=as.integer(binomial),
        nparams=as.integer(nparams),
        gendata=as.integer(gendata),
        nindvarying=as.integer(nindvarying),
        nindvaryingoffdiagonals=as.integer((nindvarying^2-nindvarying)/2),
        nt0varstationary=as.integer(nt0varstationary),
        nt0meansstationary=as.integer(nt0meansstationary),
        t0varstationary=matrix(as.integer(t0varstationary),ncol=2),
        t0meansstationary=matrix(as.integer(t0meansstationary),ncol=2),
        derrind=array(as.integer(derrind),dim=ndiffusion),
        ndiffusion=as.integer(ndiffusion),
        driftdiagonly = as.integer(driftdiagonly),
        intoverpop=as.integer(intoverpop),
        nlmeasurement=as.integer(nlmeasurement),
        nopriors=as.integer(nopriors),
        savescores=as.integer(savescores),
        savesubjectmatrices=as.integer(savesubjectmatrices)
      ))
    
    if(ctm$n.TIpred == 0) tipreds <- array(0,c(0,0))
    standata$tipredsdata <- as.matrix(tipreds)
    standata$nmissingtipreds <- as.integer(length(tipreds[tipreds== 99999]))
    
    standata$ntipredeffects <- as.integer(ifelse(n.TIpred > 0, as.integer(max(TIPREDEFFECTsetup)), 0))
    standata$TIPREDEFFECTsetup <- apply(TIPREDEFFECTsetup,c(1,2),as.integer,.drop=FALSE)
    standata$tipredsimputedscale <- ctm$tipredsimputedscale
    standata$tipredeffectscale <- ctm$tipredeffectscale
    
    
    #drift, jacobian off diagonal check
    mx=listOfMatrices(ctm$pars)
    driftcint <- rbind(cbind(mx$DRIFT[1:n.latent,1:n.latent,drop=FALSE],mx$CINT),0)
    z=c()
    for(ri in 1:nrow(driftcint)){
      if(all(driftcint[ri,-ri] %in% 0) && all(driftcint[-ri,ri] %in% 0)) z[ri]=0L else z[ri]=1L
    }
    standata$drcintoffdiag <- array(as.integer(c(z,1)),dim=nrow(driftcint))

    jac <- mx$JAx
        z=c()
    for(ri in 1:nrow(jac)){
      if(all(jac[ri,-ri] %in% 0) && all(jac[-ri,ri] %in% 0)) z[ri]=0L else z[ri]=1L
    }
    standata$jacoffdiag <- array(as.integer(c(z,1)),dim=nrow(jac))
    standata$jacoffdiagindex <- array(as.integer(sort(unique(c(1:n.latent,which(standata$jacoffdiag ==1))))))
    standata$njacoffdiagindex <- as.integer(length(standata$jacoffdiagindex))
    # browser()
    standata$sJycolindex <- array(unique(as.integer(apply(mx$Jy,1, function(x) which(!x%in%0)))))
    standata$sJycolindex<-array(standata$sJycolindex[!is.na(standata$sJycolindex)])
    standata$sJycolindexsize <- as.integer(length(standata$sJycolindex))
    
    # #jacobians
    # standata$sJAxdrift <- array(jacobianelements(ctm$jacobian$JAx,mats=mx,
    #   remove='drift',
    #   ntdpred=ctm$n.TDpred,when=2,matsetup=matsetup,returndriftonly=TRUE),
    #   dim=c(standata$nlatentpop,standata$nlatentpop))
    # 
    # standata$sJylambda <- array(jacobianelements(ctm$jacobian$Jy,mats=mx,
    #   remove='lambda',
    #   ntdpred=ctm$n.TDpred,when=4,matsetup=matsetup,returnlambdaonly=TRUE),
    #   dim=c(standata$nmanifest,standata$nlatentpop))
    
    standata$difftype <- 2L;
    standata$dotipred <- 1L;
    
    if(!recompile){ #then use finite diffs for some elements
      # standata$sJAxfinite <- array(as.integer(unique(c(which(matrix(ctm$jacobian$JAx %in% #which rows of jacobian are not simply drift / fixed / state refs
      #   jacobianelements(ctm$jacobian$JAx,mats=mx,remove=c('drift','fixed'),
      #     ntdpred=ctm$n.TDpred,when=2,matsetup=matsetup),standata$nlatentpop,standata$nlatentpop), 
      # arr.ind = TRUE))))) #[,'col'] maybe split up into row / column?
      
      #collect row and column of complicated jacobian elements into vector
      standata$sJAxfinite <- array(as.integer(unique(
        unlist(matsetup[matsetup$matrix %in% 52 & matsetup$when == -999 & matsetup$copyrow < 1,c('row','col')]))))# ])))
      
      #if any needed, set all as temp workaround for bad subsetting in finite diff jacobian
      if(length(standata$sJAxfinite) > 0) standata$sJAxfinite <- 1:standata$nlatentpop
      
      standata$nsJAxfinite <- length(standata$sJAxfinite)
      whichfinite <- matsetup$row %in% standata$sJAxfinite & matsetup$matrix %in% 52
      matsetup <- matsetup[!whichfinite,]
      matvalues <- matvalues[!whichfinite,]
    }
    if(recompile){
      standata$sJAxfinite <- array(as.integer(c()))
      standata$nsJAxfinite <- 0L
    }
    
    
    
    #add subject variability indices to data
    for(mati in c(names(mats$base),'asymCINT','asymDIFFUSION','DIFFUSIONcov')){
      sname <- paste0(mati,'subindex')
      standata[[sname]] <- as.integer((get(sname)))
    }
    
    #state dependence
    statedependence=rep(0L,4)
    multiplicativenoise = 0L
    if(any(matsetup$when == 2) ||
        any(grepl('state[',ctm$jacobian$JAx,fixed=TRUE)) ||
        any(grepl('state[',ctm$calcs$driftcint,fixed=TRUE)) ||
        any(grepl('state[',ctm$calcs$diffusion,fixed=TRUE)) 
    ) statedependence[2] = 1L
    
    if(any(matsetup$when == 2 & matsetup$matrix == 4) ||
        any(grepl('state[',ctm$calcs$diffusion,fixed=TRUE)) 
    ) multiplicativenoise = 1L
    
    standata$statedependence <- statedependence
    standata$multiplicativenoise <- multiplicativenoise
    standata$choleskymats<- ifelse(ctm$covmattransform=='unconstrainedcorr',0L,1L)
    if(!ctm$covmattransform %in% c('unconstrainedcorr','cholesky')) stop('covtransform must be either "unconstrainedcorr" or "cholesky"')
    
    standata$matsetup <- apply(matsetup[,-1],c(1,2),as.integer,.drop=FALSE) #remove parname and convert to int
    standata$matvalues <- apply(matvalues,c(1,2),as.numeric)
    standata$nmatrices <- as.integer(nmatrices)
    standata$matrixdims <- matrixdims
    standata$nrowmatsetup <- as.integer(nrow(matsetup))
    
    standata$sdscale <- array(as.numeric(sdscale),dim=length(sdscale))
    
    standata$approxct <- 0L
    if(!is.null(ctm$approxct)) standata$approxct <- as.integer(ctm$approxct)
    
    standata$taylorheun <- 0L
    if(!is.null(ctm$taylorheun)) standata$taylorheun <- as.integer(ctm$taylorheun)
    
    #fixed hyper pars
    if(!is.null(ctm$fixedrawpopchol)) {
      standata$fixedrawpopmeans = array(ctm$fixedrawpopmeans)
      standata$fixedrawpopchol= ctm$fixedrawpopchol
    }
    standata$fixedhyper <- as.integer(ifelse(is.null(ctm$fixedrawpopchol),0,1))
    if(is.null(ctm$fixedrawpopchol)) {
      standata$fixedrawpopmeans = array(0,dim = 0)
      standata$fixedrawpopchol = matrix(0,0,0)
    }
    standata$fixedsubpars <- as.integer(!is.null(ctm$fixedsubpars))
    if(!is.null(ctm$fixedsubpars)) standata$fixedindparams <- 
      ctm$fixedsubpars else standata$fixedindparams <-array(0,dim=c(0,0))
    
    if(fit){
      # if(gendata && stanmodels$ctsmgen@model_code != stanmodeltext) recompile <- TRUE
      # if(!gendata && paste0(stanmodels$ctsm@model_code) != paste0(stanmodeltext)) recompile <- TRUE
      
      STAN_NUM_THREADS <- Sys.getenv('STAN_NUM_THREADS',unset=NA)
      Sys.setenv(STAN_NUM_THREADS=cores)
      
      if(recompile || forcerecompile) {
        message('Compiling model...') 
        sm <- stan_model(model_name='ctsem', model_code = c(stanmodeltext),auto_write=TRUE
          # ,allow_undefined = TRUE,verbose=TRUE,
          # includes = paste0(
          #   '\n#include "', file.path(getwd(), 'syl2.hpp'),'"',
          #   '\n')
          )
      }
      if(!recompile && !forcerecompile) {
        if(!gendata) sm <- stanmodels$ctsm else sm <- stanmodels$ctsmgen
      }
      if(!is.null(inits) & any(inits!=0)){
        sf <- stan_reinitsf(sm,standata)
        staninits <- list(stan1=constrain_pars(sf,inits))
        # staninits=list()
        if(chains > 1){
          for(i in 2:chains){
            staninits[[i]]<-staninits[[1]]
          }
        }
      }
      
      if(!is.null(inits) & class(inits) !='list') staninits=inits #if 0 init
      if(is.null(inits)){
        staninits=list()
        # if(chains > 0){
        # init=0
        # message('Finding good start values...')
        # if(is.null(ctm$fixedsubpars)) freesubpars <- TRUE else freesubpars <- FALSE
        for(i in 1:(chains)){
          #   if(nindvarying > 0 && !intoverpop & !optimize) {
          #     ctm$fixedsubpars <- matrix( rnorm(nindvarying*nsubjects,0,.1),ncol=nindvarying)
          #     if(i==1){
          #     sf <- stan_reinitsf(sm,standata)
          #     fitb=suppressMessages(ctStanFit(datalong = datalong,ctm = ctm,optimize=TRUE,fit=TRUE,inits=init,
          #       savescores=FALSE,gendata = FALSE,
          #       optimcontrol = list(estonly=TRUE,deoptim=FALSE,isloops=0,finishsamples=2,tol=1e-4),verbose=0,...))
          #     init <- c(fitb$stanfit$rawest)+ rnorm(length(init),0,.05)
          #     } else init = init + rnorm(length(init),0,.05)
          #     staninits[[i]] <- constrain_pars(sf,c(init,ctm$fixedsubpars))
          #     if(i==chains & freesubpars) ctm$fixedsubpars <- NULL
          #   } else {
          staninits[[i]]=list(
            # baseindparams=array(rnorm(ifelse(intoverpop,0,nsubjects*nindvarying),0,.1),dim = c(ifelse(intoverpop,0,nsubjects),ifelse(intoverpop,0,nindvarying))),
            # eta=array(stats::rnorm(nrow(datalong)*n.latent,0,.1),dim=c(nrow(datalong),n.latent)),
            # tipredeffectparams=array(rnorm(standata$ntipredeffects,0,.1)) 
          )
          
          if(standata$fixedhyper==0){
            staninits[[i]]$rawpopmeans=array(rnorm(nparams,0,.1))
            staninits[[i]]$rawpopsdbase=array(rnorm(nindvarying,0,.1))
            staninits[[i]]$sqrtpcov=array(rnorm((nindvarying^2-nindvarying)/2,0,.1))
          }
          if(!is.na(ctm$rawpopsdbaselowerbound) & standata$fixedhyper==0) staninits[[i]]$rawpopsdbase=exp(staninits[[i]]$rawpopsdbase)
          # }
          # }
        }
      }
      
      if(!optimize){
        
        #control arguments for rstan
        # if(is.null(control$adapt_term_buffer)) control$adapt_term_buffer <- min(c(iter/10,max(iter-20,75)))
        if(is.null(control$adapt_delta)) control$adapt_delta <- .8
        if(is.null(control$adapt_window)) control$adapt_window <- 20
        if(is.null(control$max_treedepth)) control$max_treedepth <- 10
        if(is.null(control$adapt_init_buffer)) control$adapt_init_buffer=20
        if(is.null(control$stepsize)) control$stepsize=.001
        if(is.null(control$metric)) control$metric='dense_e'
        
        
        message('Sampling...')
        
        stanargs <- list(object = sm, 
          # enable_random_init=TRUE,
          init_r=.1,
          init=staninits,
          refresh=20,
          iter=iter,
          data = standata, chains = chains, control=control,
          cores=cores,
          ...) 
        if(plot==TRUE) stanfit <- do.call(stanWplot,stanargs) else stanfit <- do.call(sampling,stanargs)
      }
      
      if(optimize==TRUE) {
        optimcontrol$cores <- cores
        optimcontrol$verbose <- verbose
        optimcontrol$nopriors <- as.logical(nopriors)
        opcall <- paste0('stanoptimis(standata = standata,sm = sm,init = inits,plot=plot,',
          paste0(gsub('list(','',paste0(deparse(optimcontrol),collapse=''),fixed=TRUE)))
        stanfit <- eval(parse(text=opcall))
        # stanfit <- rlang::exec(stanoptimis,!!!optimcontrol,standata = standata,sm = sm,init = inits, cores=cores, verbose=verbose,nopriors=as.logical(nopriors))
        # stanfit <- stanoptimis(standata = standata,sm = sm,init = inits, cores=cores, verbose=verbose,nopriors=as.logical(nopriors))
      }
      
      if(is.na(STAN_NUM_THREADS)) Sys.unsetenv('STAN_NUM_THREADS') else Sys.setenv(STAN_NUM_THREADS = STAN_NUM_THREADS) #reset sys env
    } # end if fit==TRUE
    #convert missings back to NA's for data output
    standataout<-unlist(standata)
    standataout[standataout==99999] <- NA
    standataout <- utils::relist(standataout,skeleton=standata)
    
    setup=list(recompile=recompile,idmap=idmap,matsetup=matsetup,matvalues=matvalues,
      popsetup=matsetup[matsetup$when==0 & matsetup$param > 0,],
      popvalues=matvalues[matvalues$when==0 & matvalues$param > 0,],
      extratforms=extratforms)
    if(fit) {
      out <- list(args=args,
        setup=setup, 
        stanmodeltext=stanmodeltext, data=standataout, ctdatastruct=datalong[c(1,nrow(datalong)),],standata=standata, ctstanmodelbase=ctstanmodelbase, ctstanmodel=ctm,stanmodel=sm, stanfit=stanfit)
      class(out) <- 'ctStanFit'
    }
    
    if(!fit) out=list(args=args,setup=setup,
      stanmodeltext=stanmodeltext,data=standataout,  ctdatastruct=datalong[c(1,nrow(datalong)),],standata=standata, ctstanmodelbase=ctstanmodelbase,  ctstanmodel=ctm)
    
    return(out)
  }
}


