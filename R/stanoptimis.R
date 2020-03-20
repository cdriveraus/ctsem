# parlp <- function()  out <- try(rstan::log_prob(sf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)

parallelStanSetup <- function(cl, standata,split=TRUE){
  cores <- length(cl)
  if(split) stansubjectindices <- split(unique(standata$subject),sort(unique(standata$subject) %% min(standata$nsubjects,cores)))
  if(!split) stansubjectindices <- lapply(1:cores,function(x) unique(standata$subject))
  if(!split && length(stansubjectindices) < cores){
    for(i in (length(stansubjectindices)+1):cores){
      stansubjectindices[[i]] <- NA
    }
  }
  
  parallel::clusterExport(cl,c('standata'),envir = environment())
  parallel::clusterApply(cl,stansubjectindices,function(subindices) {
    # require(Rcpp)
    library(ctsem)
    if(length(subindices) < length(unique(standata$subject))) standata <- standatact_specificsubjects(standata,subindices)
    if(!1 %in% subindices) standata$nopriors <- 1L
    # future(globals = c('sm','standata'),
    #   packages=c('ctsem','rstan'),
    if(1==99) sm=1
    g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    assign('smf',stan_reinitsf(sm,standata),pos = g)
    
    # env <- new.env(parent=globalenv())
    # environment(parlp) <- env
    # assign('parlp',parlp,pos = g)
    
    NULL
  })
  NULL
}



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
#' \donttest{
#' sf <- stan_reinitsf(ctstantestfit$stanmodel,ctstantestfit$standata)
#' }
stan_reinitsf <- function(model, data,fast=FALSE){
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))
  
  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<-sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE,
    control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
  
  return(sf)
}

flexsapply <- function(cl, X, fn,cores=1){
  if(cores > 1) parallel::parSapply(cl,X,fn) else sapply(X, fn)
}

flexlapply <- function(cl, X, fn,cores=1,...){
  if(cores > 1) parallel::parLapply(cl,X,fn,...) else lapply(X, fn,...)
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
  return(standatamerged)
}  


standatalongobjects <- function() {
  longobjects <- c('subject','time','dokalmanrows','nobs_y','ncont_y','nbinary_y','Y','tdpreds', 'whichobs_y','whichbinary_y','whichcont_y')
  return(longobjects)
}

standatatolong <- function(standata){
  long <- lapply(standatalongobjects(),function(x) standata[[x]])
  names(long) <- standatalongobjects()
  long <- data.frame(long) #,simplify=data.frame(subject=standata$subject, time=standata$time
  # colnames(long)[colnames(long) %in% 'Y'] <- paste0('Y.1'
  # colnames(long)[colnames(long) %in% 'tdpreds'] <- 'tdpreds.1'
  return(long)
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

standataFillTime <- function(standata, times){
  long <- standatatolong(standata)
  nlong <- do.call(rbind,
    lapply(1:max(long$subject), function(si){
      stimes <- times[!times %in% long$time[long$subject==si]]
      # nadata <- matrix(99999, nrow=length(stimes),ncol=standata$nmanifest)
      # nadata <- cbind(nadata,matrix(0, nrow=length(stimes),ncol=standata$ntdpred))
      # colnames(nadata) <- colnames(long)[c(grep('^Y.',colnames(long)),grep('^tdpreds.',colnames(long)))] #c(paste0('Y.',1:standata$nmanifest),paste0(rep('tdpreds.',standata$ntdpred),seq_len(standata$ntdpred)))
      # subjdata <- data.frame(subject=si,time=stimes, 
      #   dokalmanrows=1L,nobs_y=0L,ncont_y=0L,nbinary_y=0L,whichobs_y=0L,whichbinary_y=0L,whichcont_y=0L,
      #   nadata)
      data.frame(subject=si,time=stimes)
    })
  )
  
  nlong <- suppressWarnings(data.frame(nlong,long[1,!colnames(long) %in% c('subject','time')]))
  nlong[,grep('(^nobs)|(^which)|(^ncont)|(^nbin)',colnames(nlong))] <- 0L
  nlong[,grep('^dokalman',colnames(nlong))] <- 1L
  nlong[,grep('^Y',colnames(nlong))] <- 99999
  nlong[,grep('^tdpreds',colnames(nlong))] <- 0
  
  
  mlong <- rbind(long,nlong)
  mlong <- mlong[order(mlong$subject,mlong$time),]
  standatamerged <- standatalongremerge(long=mlong, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(mlong))
  return(standatamerged)
}




stan_constrainsamples<-function(sm,standata, samples,cores=2, cl=NA){
  # smf <- stan_reinitsf(model = sm,data = standata)
  message('Computing quantities for ', nrow(samples),' samples...')
  # est1=NA
  # class(est1)<-'try-error'
  # i=0
  if(nrow(samples)==1) cores <- 1
  # while(i < nrow(samples) && class(est1)=='try-error'){
  #   i=i+1
  #   est1=try(constrain_pars(smf, upars=samples[i,]))
  # }
  # if(class(est1)[1]=='try-error') stop('All samples generated errors! Respecify, try stochastic optimizer, try again?')
  # 
  
  #create via eval due to parallel communication rubbish
  eval(parse(text="tparfunc <- function(x, parallel = TRUE){ 
    if(parallel) library(ctsem)
    if(!is.null(standata$savescores) && !standata$savescores){
      standata$dokalmanrows <- 
        as.integer(c(1,standata$subject[-1] - standata$subject[-standata$ndatapoints]))
    }
    smf <- stan_reinitsf(sm,standata)
    out <- list()
    for(li in 1:length(x)){
      out[[li]] <- try(rstan::constrain_pars(smf, upars=samples[x[li],]))
      if(any(
        sapply(out[[li]], function(x){
          test <- length(x) > 0 && 
            (any(c(is.nan(x),is.infinite(x),is.na(x))))
          # if(any(test)) print(x)
          return(test)
        })
      )){
        class(out[[li]]) <- c(class(out[[li]]),'try-error')
      }
    }
    return(out)
  }"))
  
  env <- new.env(parent = globalenv(),hash = TRUE)
  environment(tparfunc) <- env
  env$standata <- standata
  env$sm <- sm
  env$samples <- samples
  if(cores > 1 && all(is.na(cl))){
    cl <- parallel::makeCluster(cores, type = "PSOCK",useXDR=TRUE)
    on.exit(parallel::stopCluster(cl),add = TRUE)
    # parallel::clusterExport(cl2, c('sm','standata','samples'),environment())
    # parallel::clusterApply(cl2,1:cores, function(x) library(ctsem))
  }
  transformedpars <- try(flexlapply(cl, 
    split(1:nrow(samples), sort((1:nrow(samples))%%cores)),tparfunc,cores=cores,parallel=cores > 1))
  
  #fix this hack
  # browser()
  if(!is.null(transformedpars[[1]][[1]]$popmeans)) transformedpars=unlist(transformedpars,recursive = FALSE)
  est1=transformedpars[[1]]
  missingsamps <-sapply(transformedpars, function(x) 'try-error' %in% class(x))
  nasampscount <- sum(missingsamps) 
  
  
  if(nasampscount > 0) {
    browser()
    # a=lapply(transformedpars[[1]],function(x) any(is.na(x)));a[unlist(a)]
    message(paste0(nasampscount,' NAs generated during final sampling of ', nrow(samples), '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
  }
  if(nasampscount < nrow(samples)){
    transformedpars <- transformedpars[!missingsamps] 
    nresamples <- nrow(samples) - nasampscount
  } else{
    message('All samples contain NAs -- returning anyway')
    nresamples <- nrow(samples) 
  }
  
  
  #this seems inefficient and messy, should be a better way...  
  transformedpars=try(tostanarray(flesh=matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = est1))
  
  # parallel::stopCluster(cl)
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





#' Optimize / importance sample a stan or ctStan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init vector of unconstrained parameter values, or character string 'random' to initialise with
#' random values very close to zero.
#' @param initsd positive numeric specifying sd of normal distribution governing random sample of init parameters,
#' if init='random' .
#' @param sampleinit either NA, or an niterations * nparams matrix of samples to initialise importance sampling.
#' @param deoptim Do first pass optimization using differential evolution? Slower, but better for cases with multiple
#' minima / difficult optimization.
#' @param stochastic Logical. Use stochastic gradient descent instead of mize (bfgs) optimizer.
#' Still experimental, worth trying for either robustness checks or problematic, high dimensional, nonlinear, problems.
#' @param plot Logical. If TRUE, plot iteration details. Probably slower.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{DEoptim.control}.
#' @param tol objective tolerance.
#' @param nopriors logical. If TRUE, a nopriors integer is set to 1 (TRUE) in the standata object -- only has an effect if 
#' the stan model uses this value. 
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when nopriors=TRUE.
#' @param subsamplesize value between 0 and 1 representing proportion of subjects to include in first pass fit. 
#' @param cores Number of cpu cores to use, should be at least 2.
#' @param is Logical. Use importance sampling, or just return map estimates?
#' @param isloopsize Number of samples of approximating distribution per iteration of importance sampling.
#' @param finishsamples Number of samples to draw for final results of importance sampling.
#' @param finishmultiply Importance sampling stops once available samples reach \code{finishsamples * finishmultiply} , then the final samples are drawn
#' without replacement from this set.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent
#' importance sampling, at risk of truncating tails.
#' @param chancethreshold drop iterations of importance sampling where any samples are chancethreshold times more likely to be drawn than expected.
#'
#' @return list containing fit elements
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @import rstan
#' @export
#' @examples
#' \donttest{
#'
#' library(rstan)
#' scode <- "
#' parameters {
#'   real y[2];
#' }
#' model {
#'   y[1] ~ normal(0, 1);
#'   y[2] ~ double_exponential(0, 2);
#' }
#' "
#'
#' sm <- stan_model(model_code=scode)
#' fit <- sampling(sm, iter = 10000)
#' summary(fit)$summary
#'
#' ## extract samples as a list of arrays
#' e <- ctExtract(fit, permuted = TRUE)
#'
#' #for ml or map estimates
#' optimis <- stanoptimis(standata = list(),sm = sm,finishsamples = 3000,cores=2)
#' optimis$optimfit
#'
#' #for posterior distributions
#' optimis <- stanoptimis(standata = list(),sm = sm,finishsamples = 3000,cores=2,tdf=5)
#'
#' apply(optimis$rawposterior,2,mean)
#' apply(optimis$rawposterior,2,sd)
#' isdiag(optimis)
#'
#' plot(density(optimis$rawposterior[,2],bw=.05))
#' points(density(e$y[,2],bw=.05),type='l',col=2)
#' }
stanoptimis <- function(standata, sm, init='random',initsd=.01,sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-12,
  decontrol=list(),
  stochastic = "auto",
  nopriors=FALSE,carefulfit=TRUE,
  subsamplesize=.5,
  plot=FALSE,
  is=FALSE, isloopsize=1000, finishsamples=500, tdf=10,chancethreshold=100,finishmultiply=5,
  verbose=0,cores=2){
  if(!is.null(standata$verbose)) {
    if(standata$verbose > 1) standata$verbose=as.integer(verbose) else standata$verbose=0L
  }
  standata$nopriors=as.integer(nopriors)
  
  if(is.null(decontrol$steptol)) decontrol$steptol=3
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-2
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  if(is.null(init)) init <- 'random'
  if(init[1] !='random') carefulfit <- FALSE
  
  clctsem <- NA #placeholder for flexsapply usage
  
  if(length(unique(standata$subject)) < cores) optimcores <- length(unique(standata$subject))  else optimcores <- cores
  
  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  #remove tipreds for first pass
  TIPREDEFFECTsetup <- standata$TIPREDEFFECTsetup
  standata$TIPREDEFFECTsetup[,] <- 0L
  standata$ntipredeffects <- 0L
  # npars <- npars-sum(TIPREDEFFECTsetup)
  
  if(standata$nindvarying > 0 && standata$intoverpop==0){ #detect subject level pars
    stochastic <- TRUE
    standata$doonesubject <- 1L
    message('Using combined stochastic gradient descent + mcmc')
    a1=standata$nparams+standata$nindvarying+
      (standata$nindvarying^2-standata$nindvarying)/2
    whichmcmcpars <- (a1+1):(a1+standata$nindvarying) #*standata$nsubjects)
  } else whichmcmcpars <- NA
  
  
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    smf <- stan_reinitsf(sm,standata)
    npars=rstan::get_num_upars(smf)
    # browser()
    
    if(all(init %in% 'random')) init <- rnorm(npars, 0, initsd)
    if(all(init == 0)) init <- rep(0,npars)
    if(length(init) != npars) init=init[1:npars]
    init[is.na(init)] <- 0
    if(standata$ntipredeffects > 0) init[length(init):(length(init)+1-standata$ntipredeffects)] <- 0
    

    
    
    if(is.na(sampleinit[1])){
      
      
      
      storedPars <- c()#matrix(0,nrow=npars,ncol=0)
      gradstore <- rep(0,npars)
      gradmem <- .9
      storedLp <- c()
      
      optimfinished <- FALSE
      on.exit({
        if(!optimfinished){
          message('Optimization cancelled -- restart from current point by including this argument:')
          message((paste0(c('init = c(',   paste0(round(storedPars,5),collapse=', '), ')'    ))))
          # message('Return inits? Y/N')
          # if(readline() %in% c('Y','y')) returnValue(storedPars)
        }},add=TRUE)
      
      #create parlp via eval because of parallel communication weirdness
      eval(parse(text=
          'parlp <- function(parm){
          out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)

        
        if("try-error" %in% class(out)) {
          out <- -1e100
          attributes(out)$gradient <- rep(NaN, length(parm))
        }
        if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
        attributes(out)$gradient[is.nan(attributes(out)$gradient)] <-
          rnorm(length(attributes(out)$gradient[is.nan(attributes(out)$gradient)]),0,100)
          
        return(out)
        }'))
      
      singletarget<-function(parm,gradnoise=TRUE) {
        a=Sys.time()
        out<- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        
        if('try-error' %in% class(out) || is.nan(out)) {
          out=-1e100
          attributes(out) <- list(gradient=rep(0,length(parm)))
        }
        storedPars <<- parm
        # storedLp <<- c(storedLp,out[1])
        b=Sys.time()
        evaltime <- b-a
        if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(evaltime,2)))
        if(gradnoise) attributes(out)$gradient <-attributes(out)$gradient *
          exp(rnorm(length(parm),0,1e-3))
        return(out)
      }
      
      
      if(optimcores==1) target = singletarget #we use this for importance sampling
      if(optimcores > 1) {
        
        clctsem=parallel::makeCluster(optimcores,useXDR=TRUE)
        on.exit({parallel::stopCluster(clctsem)},add=TRUE)
        
        #crazy trickery to avoid parallel communication pauses
        env <- new.env(parent = globalenv(),hash = TRUE)
        environment(parlp) <- env
        
        parallel::clusterExport(clctsem,varlist = 'sm',envir = environment())
        
        target<-function(parm,gradnoise=TRUE) {
          whichframe <- which(sapply(lapply(sys.frames(),ls),function(x){ 'clctsem' %in% x}))
          a=Sys.time()
          out2<- parallel::clusterCall(clctsem,parlp,parm)
          b=Sys.time()
          b-a
          
          out <- try(sum(unlist(out2)),silent=TRUE)
          attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
          
          if('try-error' %in% class(out) || is.nan(out)) {
            out=-1e100
            attributes(out) <- list(gradient=rep(0,length(parm)))
          } 
          
          # if(plot){
          #   storedLp <<- c(storedLp,out[1])
          #   # attributes(out)$gradient <- (1-gradmem)*attributes(out)$gradient + gradmem*gradstore
          #   gradstore <<- cbind(gradstore,attributes(out)$gradient)
          #   gstore <- log(abs(gradstore))*sign(gradstore)
          #  if(runif(1,0,1) > .9) {
          #    par(mfrow=c(1,1))
          #    matplot(tail(t(gstore),50),type='l')
          #  }
          # }
          storedPars <<- parm
          # storedLp <<- c(storedLp,out[1])
          if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',b-a))
          if(gradnoise) attributes(out)$gradient <-attributes(out)$gradient * 
            exp( rnorm(length(attributes(out)$gradient),0,1e-3))
          return(out)
        }
        
      } #end multicore setup
      

      
      if(deoptim){ #init with DE
        message('Using differential evolution for initialization')
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
          
          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]
          
          lp2 = function(parm) {
            out<- -target(parm)
            attributes(out)$gradient <- -attributes(out)$gradient
            return(out)
          }
          
          deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
          deinit[2,] <- rnorm(npars,0,.0002)
          if(length(init)>1 & try2) { #if better fit found during IS
            deinit[1,] <- unconstrain_pars(smf,init)
            if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
          }
          decontrollist$initialpop=deinit
          decontrollist$NP = NP
          
          if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata)
          if(optimcores==1) smf<-stan_reinitsf(sm,standata)
          
          optimfitde <- suppressWarnings(DEoptim::DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
            control = decontrollist))
          # init=constrain_pars(object = smf,optimfitde$optim$bestmem)
          init=optimfitde$optim$bestmem
          bestfit = -optimfitde$optim$bestval
        } else stop(paste0('use install.packages(\"DEoptim\") to use deoptim')) #end require deoptim
      }
      
      
      
      
      mizelpg=list( #single core mize functions
        fg=function(pars){
          r=-target(pars)
          r=list(fn=r[1],gr= -attributes(r)$gradient)
          return(r)
        },
        fn=function(x) -target(x),
        gr=function(pars) -attributes(target(pars))$gradient
      )
      
      if(stochastic=='auto' && npars > 50){
        message('> 50 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE
      
      if(carefulfit && !deoptim){ #init using priors
        message('Doing 1st pass with priors on reduced data set')
        nopriorsbak <- standata$nopriors
        taylorheun <- standata$taylorheun
        standata$nopriors <- as.integer(0)
        standata$taylorheun <- 1L
        sdscale <- standata$sdscale
        # standata$sdscale <- sdscale * 1e-6
        # tipredeffectscale <- standata$tipredeffectscale
        # standata$tipredeffectscale <- tipredeffectscale* 1e-5
        smlnsub <- min(standata$nsubjects,max(min(30,cores*2),ceiling(standata$nsubjects * subsamplesize)))
        standatasml <- standatact_specificsubjects(standata,
          sample(unique(standata$subject),smlnsub))
        # smlndat <- min(standatasml$ndatapoints,ceiling(max(standatasml$nsubjects * 10, standatasml$ndatapoints*.5)))
        # standatasml$dokalmanrows[sample(1:standatasml$ndatapoints,smlndat)] <- 0L
        # standatasml$dokalmanrows[match(unique(standatasml$subject),standatasml$subject)] <- 1L #ensure first obs is included for t0var consistency
        if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standatasml)
        if(optimcores==1) smf<-stan_reinitsf(sm,standatasml)
        
        if(!stochastic) {
          optimfit <- mize(init, fg=mizelpg, max_iter=99999,
            method="L-BFGS",memory=100,
            line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
            abs_tol=1e-4,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
          optimfit$value = optimfit$f
        }

        if(stochastic) {
          optimfit <- sgd(init, fitfunc = function(x) target(x),
            whichmcmcpars=whichmcmcpars,nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
          plot=plot,itertol=1e-1,deltatol=1e-2,maxiter=500)
        }
        
        standata$nopriors <- as.integer(nopriorsbak)
        standata$taylorheun <- as.integer(taylorheun)
        standata$sdscale <- sdscale
        
        init = optimfit$par
      } #end carefulfit init
      
      
      message('Optimizing...')
      if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata)
      if(optimcores==1) smf<-stan_reinitsf(sm,standata)
      
      if(!stochastic){
        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          method="L-BFGS",memory=100,
          # check_conv_every = 5,
          # try_newton_step=TRUE,
          line_search='Schmidt',c1=1e-4,c2=.9,step0='schmidt',ls_max_fn=999,
          abs_tol=NULL,grad_tol=NULL,rel_tol=tol,step_tol=NULL,ginf_tol=NULL)
        optimfit$value = -optimfit$f
        init = optimfit$par
        if(is.infinite(bestfit)) stochastic<-TRUE
      }
      
      if(stochastic){#  || #carefulfit
        if(is.infinite(bestfit)) message('Switching to stochastic optimizer -- failed initialisation with bfgs')
        if(!stochastic && carefulfit) message('carefulfit = TRUE , so checking for improvements with stochastic optimizer')
        optimfit <- sgd(init, fitfunc = target,
          whichmcmcpars=whichmcmcpars,nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
          ndatapoints=standata$ndatapoints,plot=plot)
      }
      
      if(standata$ntipred > 0){ 
        # 
        # standata$tipredeffectscale <- tipredeffectscale
        standata$dokalmanrows[] <- 1L
        # init = optimfit$par
        optimbase = optimfit$par#[-(length(optimfit$par):(length(optimfit$par)-max(TIPREDEFFECTsetup)+1))]
        finished <- FALSE
        found <- 0
        tia=standata$TIPREDEFFECTsetup
        while(!finished){
          if(!is.null(standata$TIpredAuto) && standata$TIpredAuto){
            message ('Looking for tipreds...')
            oldtia=tia
            fit <- list(stanfit=list(rawest=optimfit$par),standata=standata,stanmodel=sm)
            tia=ctTIauto(fit)  
            tia[tia > .05] = 0L
            a <- 0
            for(ci in 1:ncol(tia)){
              for(ri in 1:nrow(tia)){
                if(tia[ri,ci] +oldtia[ri,ci] > 0){
                  a <- a+1
                  tia[ri,ci] <- as.integer(a)
                }
              }
            }
            
            
            if(max(tia) > found){
              tiinit=rep(0,max(tia))
              oldtiinit = tail(optimfit$par,found)[oldtia[oldtia>0]]
              tiinit[tia[tia>0&oldtia>0]] <- oldtiinit
              standata$TIPREDEFFECTsetup <- array(as.integer(tia),dim=dim(tia))
              found <- max(tia)
              message('Found ',max(tia),' viable TIpreds')
              standata$ntipredeffects <- as.integer(max(tia))
              init = head(optimfit$par,length(optimbase))
              init = c(init,tiinit)
            } else finished <- TRUE
          }
          # if(!finished){
            
            if(!standata$TIpredAuto){
              finished <- TRUE
              init = c(init,rep(0,max(TIPREDEFFECTsetup)))
              standata$TIPREDEFFECTsetup[,] <- TIPREDEFFECTsetup
              standata$ntipredeffects <- as.integer(max(TIPREDEFFECTsetup))
            }
            
            if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata)
            if(optimcores==1) smf<-stan_reinitsf(sm,standata)
            if(!stochastic) {
              optimfit <- mize(init, fg=mizelpg, max_iter=99999,
                method="L-BFGS",memory=100,
                line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
                abs_tol=tol*ifelse(finished,1,10000),grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
              optimfit$value = -optimfit$f
            }
            if(stochastic && finished) optimfit <- sgd(init, fitfunc = target, 
              whichmcmcpars=whichmcmcpars,nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
              plot=plot)
            if(stochastic && !finished) optimfit <- sgd(init, fitfunc = target, 
              whichmcmcpars=whichmcmcpars,nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
              plot=plot,itertol=1e-1,deltatol=1e-2)
          # }
        }
        smf <- stan_reinitsf(sm,standata)
      }
      
      
      
      bestfit <-optimfit$value
      est2=optimfit$par #unconstrain_pars(smf, est1)
      
      if(standata$nindvarying > 0 && standata$intoverpop==0){ #recompose into single model
        standata$doonesubject <- 0L
        a1=standata$nparams+standata$nindvarying+
          (standata$nindvarying^2-standata$nindvarying)/2
        whichmcmcpars <- (a1+1):(a1+standata$nindvarying*standata$nsubjects)
        est3=est2[-whichmcmcpars]
        est3=c(est3,(optimfit$mcmcpars))
        est2=est3
        if(standata$ntipredeffects > 0) est3 <- c(est3,est2[(a1+1):length(a1)])
        if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata)
        if(optimcores==1) smf<-stan_reinitsf(sm,standata)
      } 
      npars = length(est2)
    }
    
    if(!estonly){
      
      # if(cores > 1) parallelStanSetup(cl=clctsem,sm,standata,split=FALSE) #parallel::clusterExport(clctsem,varlist = c('
      base <- target(est2,gradnoise = FALSE)
      basegrad <- attributes(base)$gradient
      
      grmat<-function(pars,step=1e-5,lpdifmin=1e-8, 
        lpdifmax=1e-4, direction=1,whichpars='all',gradmod=FALSE){
        if('all' %in% whichpars) whichpars <- 1:length(pars)
        hessout <- flexsapply(cl = clctsem, cores = 1, whichpars, function(i) {
          stepsize <- step * direction
          colout <- NA
          dolpchecks <- TRUE #set to true to try the log prob checks again..
          # 
          while(any(is.na(colout)) && abs(stepsize) > 1e-20){
            stepsize <- stepsize * .1
            lpdifok<-FALSE
            lpdifcount <- 0
            lpdifdirection <- 0
            lpdifmultiplier <- 1
            # message('par',i)
            while(!lpdifok & lpdifcount < 15){
              # message(paste(i,'  col=',colout,'  lpdifmultiplier=',lpdifmultiplier, '  stepsize=',stepsize))
              lpdifok <- TRUE
              lpdifcount <- lpdifcount + 1
              uppars<-pars
              uppars[i]<-pars[i]+stepsize / ifelse(gradmod,(abs(basegrad[i])),1)
              
              suppressMessages(suppressWarnings(try({uplp<- target(uppars,gradnoise=FALSE)},silent=TRUE))) #try(smf$log_prob(upars=uppars,adjust_transform=TRUE,gradient=TRUE)) #lpg(uppars)
              # storedPars <<- cbind(storedPars,matrix(uppars))
              # storedLp <<- c(storedLp,uplp[1])
              
              if('try-error' %in% class(uplp)){
                lpdifok <- TRUE
                upgrad <- rep(NA,length(pars))
                dolpchecks <- FALSE
              } else{
                upgrad= attributes(uplp)$gradient -basegrad
                
                # print((  (upgrad[i] * (-abs(uppars[i]-pars[i]))) / (uplp[1]-bestfit[1]) ))
                # upgrad = upgrad / (  (upgrad[i] * (-abs(uppars[i]-pars[i]))) / (uplp[1]-bestfit[1]) ) #linearised upgrad
                
                if(dolpchecks){
                  if(abs(base-uplp) > lpdifmax) {
                    # message(paste0('decreasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== 1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * (1e-2 * lpdifmultiplier)
                    lpdifdirection <- -1
                  }
                  if(abs(uplp-base) < lpdifmin ) { #include sufficient gradient[i] checks #|| upgrad[i] < 1e-2
                    # 
                    # message(paste0('increasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== -1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * (100 * lpdifmultiplier)
                    lpdifdirection <- 1
                  }
                  if(any(is.na(c(uplp)))) stepsize = stepsize * .1
                }
              }
            }
            colout<- (upgrad) /(stepsize) * ifelse(gradmod,(abs(basegrad[i])),1)
          }
          rbind(colout)
        })
        return((hessout))# + t(hessout))/2)
      }
      
      # A more numerically stable way of calculating log( sum( exp( x ))) Source:
      # http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
      log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
      }
      
      
      if(is.na(sampleinit[1])){
        
        message('Estimating Hessian')
        
        # hesslist <- list()
        # remainingpars <- c()
        # for(direction in c(1,-1)){
        #   whichpars <- 1:npars
        #   hess = matrix(0,npars,npars)
        #   hessbad <- hess
        #   steps=c(1e-5)
        #   for(stepi in 1:length(steps) ){
        #     # message(steps[stepi])
        #     if(length(whichpars) > 0){
        #       hesstry <- hess
        #       hesstry[,whichpars]= at(pars=est2,step=steps[stepi], direction=direction,gradmod = TRUE,whichpars=whichpars)
        #       h=(hesstry+t(hesstry))/2
        #       rankifremoved <- sapply(1:ncol(h), function (x) qr(h[,-x,drop=FALSE])$rank)
        #       hessgood<- which(rankifremoved < max(rankifremoved))
        #       if(length(whichpars) < length(which(rankifremoved == max(rankifremoved)))){
        #         # message('Worse!')
        #         next #if worse, don't use this stepsize
        #       }
        #       whichpars <- which(rankifremoved == max(rankifremoved)) #update whichpars
        #       
        #       # hn=jacobian(function(x) attributes(target(x))$gradient,est2)
        #       # eig <- eigen(h)$values
        #       # hessgood <- which(abs(eig) > 1e-4) #which columns / pars are non singular
        #       # bothgood <- which(apply(hess,2,mean,na.rm=TRUE) != 0) #which columns of old good hess are non zero
        #       # bothgood <- bothgood[bothgood %in% hessgood]
        #       hess[,hessgood] <- hesstry[,hessgood] 
        #       # newgood <- hessgood[!hessgood%in% bothgood]
        #       # hess[,newgood] = hesstry[,newgood]
        #       
        #       # print(whichpars)
        #       hessbad <- hessbad/stepi * (stepi-1) + hesstry/stepi #cumulative update
        #     }
        #   }
        #   if(direction==1){
        #     remainingpars <- whichpars
        #     hesslist[[1]] <- hess
        #     hesslist[[2]] <- hessbad
        #   } else {
        #     hesslist[[3]] <- hess
        #     hesslist[[4]] <- hessbad
        #   }
        # }
        # badpars <- unique(c(whichpars,remainingpars) )
        # onesidedpars <- c(whichpars[!whichpars %in% remainingpars],
        #   remainingpars[!remainingpars %in% whichpars])
        # badpars <- badpars[!badpars %in% onesidedpars]
        # 
        # mats <- ctStanMatricesList()$all
        # nameref <- function(p){
        #   r <- which(standata$matsetup[,'param'] %in% p)
        #   
        #   out <- c()
        #   for(ri in r){
        #     out <- c(out,paste0(names(mats)[standata$matsetup[ri,'matrix']],'[',
        #       standata$matsetup[ri,'row'],',',
        #       standata$matsetup[ri,'col'],'] ')
        #     )
        #   }
        #   if(length(r) < length(p)) out <-c(out, p[!p%in% standata$matsetup[,'param']])
        #   return(out)
        # }        
        # 
        # hess <- (hesslist[[1]] + hesslist[[3]] ) /2
        # hess[,onesidedpars] <- hess[,onesidedpars] * 2
        # hess[,badpars] <- (( hesslist[[2]] + hesslist[[4]]) /2)[,badpars]
        # hess <- (t(hess)+hess )/2
        # rankifremoved <- sapply(1:ncol(hess), function (x) qr(hess[,-x,drop=FALSE])$rank)
        # # newbadpars <- which(rankifremoved == max(rankifremoved))
        
        
        whichhesspars = (1:npars)
        if(standata$intoverpop == 0 && standata$nindvarying > 0) whichhesspars <- whichhesspars[-whichmcmcpars]
        hessup= grmat(pars=est2,step=1e-8, direction=1,whichpars = whichhesspars,gradmod = FALSE)[whichhesspars,,drop=FALSE]
        hessdown= grmat(pars=est2,step=1e-8, direction=1,whichpars = whichhesspars,gradmod = FALSE)[whichhesspars,,drop=FALSE]
        hess = (hessup+hessdown)/2
        cholcovup = try(suppressWarnings(t(chol(solve(-hessup)))),silent = TRUE)
        cholcovdown = try(suppressWarnings(t(chol(solve(-hessdown)))),silent = TRUE)
        cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
        # 
        if('try-error' %in% class(cholcov)){
          if(!'try-error' %in% class(cholcovup)){
            hess = hessup
            message('One sided hessian used for std error estimation')
          } else if(!'try-error' %in% class(cholcovdown)){
            message('One sided hessian used for std error estimation')
            hess = hessdown
          } else if(!is){
            warning('Approximate hessian used for std error estimation --  treat SE\'s with caution, consider respecification / priors / sampling.')
          }
        }
      
        mcov=MASS::ginv(-hess) #-optimfit$hessian)
        mcov=as.matrix(Matrix::nearPD(mcov,conv.norm.type = 'F')$mat)
        mchol = t(chol(mcov))
        if(standata$intoverpop == 0 && standata$nindvarying > 0){
          mcmcsamps <- matrix(rnorm(finishsamples*length(whichmcmcpars)),nrow=finishsamples,ncol=length(whichmcmcpars))
          mcmccov = cov(mcmcsamps)
          mcmcchol = diag(1,length(whichmcmcpars))
          fullchol <- matrix(0,nrow(mchol)+nrow(mcmcchol),nrow(mchol)+nrow(mcmcchol))
          fullchol[whichmcmcpars,whichmcmcpars] <- mcmcchol
          fullchol[-whichmcmcpars,-whichmcmcpars] <- mchol
          mcov <- fullchol %*% t(fullchol)
          mchol <- fullchol
        }
      }
      
      
      if(!is.na(sampleinit[1])){
        mcov = cov(sampleinit)*1.5+diag(1e-6,ncol(sampleinit))
        est2 = apply(sampleinit,2,mean)
        bestfit = 9e100
        optimfit <- suppressWarnings(list(par=sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
      }
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
      
      # domix=F
      
      if(!is) {
        nresamples = finishsamples
        # if(!domix){  
        message('Getting ',finishsamples,' samples from Hessian for interval estimates...')
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + (mchol) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        # }
        
        # if(domix){
        #   message('Getting ',finishsamples,' samples from estimated density for interval estimates...')
        #   lplim <- 1
        #   dmix=densitymixt::tmix(dfvec = 50,fixeddf = FALSE,
        #     samples = t(storedPars[,storedLp > max(storedLp)-lplim,drop=FALSE]),
        #     lpsamps = storedLp[storedLp > max(storedLp)-lplim],
        #     ngen = finishsamples,priors = FALSE,stochastic=FALSE,tol=1e-8)
        #   resamples <- dmix$fit$gensamples
        # }
        
        
        # for(i in 1:ncol(resamples)){
        # plot(density(resamples[,i]),main=i)
        # points(density(dmix$fit$gensamples[,i]),type='l',col='red')
        # }
        # resamples <- dmix$fit$gensamples
      }
      if(is){
        message('Importance sampling...')
        if(optimcores > 1)  parallelStanSetup(cl=clctsem,standata,split=FALSE)
        targetsamples <- finishsamples * finishmultiply
        # message('Adaptive importance sampling, loop:')
        j <- 0
        while(nrow(samples) < targetsamples){
          j<- j+1
          if(j==1){
            # if(!npd)
            # if(!domix) 
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            
            # if(domix){
            #   lplim <- 10
            #   dmix=densitymixt::tmix(dfvec = 50,fixeddf = F,
            #     samples = t(storedPars[,storedLp > max(storedLp)-lplim,drop=FALSE]),
            #     lpsamps = storedLp[storedLp > max(storedLp)-lplim],
            #     targettype = 0,
            #     ngen = isloopsize,priors = FALSE,stochastic=F,tol=1e-12,plot=T)
            #   samples <- dmix$fit$gensamples
            #   prop_dens <- dmix$fit$genlp
            #   
            # }
            
          } else {
            # if(!domix){
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf))
            # }
            # if(domix){
            #   lplim=3
            #   # 
            # dmix=densitymixt::tmix(dfvec = 10,fixeddf = F,
            #   samples = samples[targetvec > max(targetvec)-lplim,],
            #   lpsamps = targetvec[targetvec > max(targetvec)-lplim],
            #   ngen = isloopsize,priors = FALSE,stochastic=F,tol=1e-8,plot=T)
            # 
            # # delta[[j]]=colMeans(resamples)
            # # mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            # samples <- rbind(samples,dmix$fit$gensamples)
            # prop_dens <- dmix$fit$genlp
            # }
            
            
            # samples <- rbind(samples, MASS::mvrnorm(isloopsize, delta[[j]],mcovl[[j]]))
          }
          # if(j > 1 || !npd)
          # if(!domix || j==1) 
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          # prop_dens <- mvtnorm::dmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]],log = TRUE)
          # prop_dens <- ctdmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]])
          
          
          standata$verbose <- 0L
          if(optimcores > 1) parallel::clusterExport(clctsem,'samples',envir = environment())
          
          # if(cores==1) parlp <- function(parm){ #remove this duplication somehow
          #   out <- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
          #   if('try-error' %in% class(out)) {
          #     out[1] <- -1e100
          #     attributes(out)$gradient <- rep(NaN, length(parm))
          #   }
          #   if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
          #   attributes(out)$gradient[is.nan(attributes(out)$gradient)] <- 
          #     rnorm(length(attributes(out)$gradient[is.nan(attributes(out)$gradient)]),0,100)
          #   return(out)
          # }
          
          target_dens[[j]] <- unlist(flexlapply(cl = clctsem, split(1:isloopsize, sort((1:isloopsize)%%cores)), function(x){
            # future(globals=c('x','samples','isloopsize'),expr={
            eval(parse(text='apply(tail(samples,isloopsize)[x,],1,parlp)'))
          },cores=cores))
          
          
          
          
          target_dens[[j]][is.na(target_dens[[j]])] <- -1e200
          if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
          if(any(target_dens[[j]] > bestfit)){
            oldfit <- bestfit
            try2 <- TRUE
            bestfit<-max(target_dens[[j]],na.rm=TRUE)
            betterfit<-TRUE
            init = samples[which(unlist(target_dens) == bestfit),]
            message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
            break
          }
          targetvec <- unlist(target_dens)
          
          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)
          
          newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          counter <- 1
          isfinished <- FALSE
          
          # while(counter < 30 && length(unique(sample(x=1:length(newsampleprob),size=length(newsampleprob),prob=newsampleprob,replace=TRUE))
          # ) < 20) {
          #   if(counter==1) message ('Sampling problematic -- trying to recover... ')
          #   counter = counter + 1
          #   if(counter == 30) stop('Importance sampling failed -- either posterior mode not found, or mode is inadequate starting point for sampling')
          #   weighted_dens <- weighted_dens /2
          #   newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          #   # plot(newsampleprob)
          # }
          sample_prob <- c(sample_prob,newsampleprob) #sum to 1 for each iteration, normalise later
          
          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
          
          
          #check max resample probability and drop earlier samples if too high
          dropuntil <- ceiling(max(c(0,which(sample_prob > (chancethreshold / isloopsize)) / isloopsize),na.rm=TRUE))*isloopsize
          if((isloopsize - dropuntil) > isloopsize) dropuntil <- dropuntil -isloopsize
          if(nrow(samples)-dropuntil < isloopsize*2) dropuntil <- nrow(samples)-isloopsize*2
          # if(length(unique(resample_i)) < 200) dropuntil <- 0 
          if(nrow(samples) <= isloopsize *2) dropuntil <- 0
          
          if(dropuntil > 0){
            # 
            # resamples <- resamples[-(0:dropuntil),,drop=FALSE]
            targetvec <- targetvec[-(0:dropuntil)]
            sample_prob <- sample_prob[-(0:dropuntil)]
            samples <- samples[-(0:dropuntil),,drop=FALSE]
          }
          
          
          # if(plot&runif(1,0,1) > .98){
          #   par(mfrow=c(1,1))
          #   lprob <- Var2 <-value <- NULL
          #   msamps = cbind(data.table(lprob=rep(log(sample_prob),ncol(samples))),
          #     data.table::melt(data.table(cbind(samples))))
          #   ggplot2::ggplot(dat=msamps,aes(y=lprob,x=value,alpha=exp(msamps$lprob))) + 
          #     geom_point(colour='red') +
          #     scale_alpha_continuous(range=c(0,1))+
          #     theme_minimal()+
          #     facet_wrap(vars(Var2),scales='free') 
          # }
          
          
          if(nrow(samples) >= targetsamples && (max(sample_prob)/ sum(sample_prob)) < chancethreshold) {
            message ('Finishing importance sampling...')
            nresamples <- finishsamples 
          }else nresamples = max(5000,nrow(samples)/5)
          
          
          # points(target_dens2[sample_prob> (1/isloopsize * 10)], prop_dens[sample_prob> (1/isloopsize * 10)],col='red')
          resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(nrow(samples) > targetsamples,FALSE,TRUE),
            prob = sample_prob / sum(sample_prob))
          # resample_i <- sample(tail(1:nrow(samples),isloopsize), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE),
          #   prob = tail(sample_prob,isloopsize) / sum(tail(sample_prob,isloopsize) ))
          
          message(paste0('Importance sample loop ',j,', ',length(unique(resample_i)), ' unique samples, from ', nresamples,' resamples of ', nrow(samples),' actual, prob sd = ', round(sd(sample_prob),4),
            ', max chance = ',max(sample_prob) * isloopsize))
          if(length(unique(resample_i)) < 100) {
            message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
            # return(est)
          }
          
          resamples <- samples[resample_i, , drop = FALSE]
          
          ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
          qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
            (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) )
          })))
          
        }
      }
    }
  }#end while no better fit
  
  if(!estonly){
    if(!is) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]
    
    transformedpars=stan_constrainsamples(sm = sm,standata = standata,samples=resamples,cores=cores, cl=clctsem)
    
    
    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if('try-error' %in% class(sds)[1]) sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, upars=lest)),
      unlist(constrain_pars(smf, upars= est2)),
      unlist(constrain_pars(smf, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2, rawposterior = resamples, cov=mcov,
      transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      standata=list(TIPREDEFFECTsetup=standata$TIPREDEFFECTsetup,ntipredeffects = standata$ntipredeffects),
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}


