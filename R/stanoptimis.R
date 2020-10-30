


# hesscalc <- function(sf,step){ #input stanfit, get hessian
#   smf <- stan_reinitsf(sf$stanmodel,sf$standata)
#   fgfunc <- function(x) rstan::log_prob(smf,x,gradient=TRUE)
#   # step=findstepsize(est = sf$stanfit$rawest,func = fgfunc,eps = 1e-2,
#   #   maxlpd = 2*length(sf$stanfit$rawest),minlpd=.1*length(sf$stanfit$rawest))
#   H=jac(pars = sf$stanfit$rawest,fgfunc = fgfunc,step = rep(step,length(sf$stanfit$rawest)))
# }

findstepsize <- function(est,func,eps=1e-3,maxlpd=3,minlpd=1e-3){
  
  fmax=func(est)
  epsbase=.01
  i=0
  n<-100
  devs <- matrix(NA,nrow=n,ncol=length(est))
  lp<-c()
  epscount <- 0
  goodcount <- 0
  
  while(i < n && goodcount < 5){
    i=i+1
    accepted <- FALSE
    devs[i,] <- rnorm(length(est))*eps
    count <- 0
    while(!accepted){
      count = count + 1
      if(count==20) stop('Couldnt determine step size!')
      fg<-try(func(devs[i,]+est))
      if('try-error' %in% class(fg)) eps <- eps * .5+1e-6 else accepted <- TRUE
    }
    
    bad=(lp<quantile(lp,.2,na.rm = TRUE) & (fmax-lp) > maxlpd) | 
      (lp>quantile(lp, .95,na.rm = TRUE) & (fmax-lp) < minlpd*.001)
    
    good <- c(1:i)[!bad]
    lp <- c(lp,fg[1])
    fl=lp-min(lp)
    fl <- sqrt(fl)
    
    frange <- fmax[1]-lp[i]
    print(frange)
    if(length(good) <=5){
      if(frange < minlpd) eps <- eps * 2
      if(frange > maxlpd) eps <- eps /2
    }
    # message(i)
    # print(eps)
    if(frange > minlpd && frange < maxlpd) goodcount <- goodcount + 1 else goodcount <- 0
    # 
    if(length(good)>5){
      if(frange < minlpd) epsbase <- epsbase * 2
      if(frange > maxlpd) epsbase <- epsbase / 2
      eps=epsbase/(sqrt(abs(apply(devs,2,function(x) coefficients(lm(fl[good]~abs(x[good])))[-1])))+1e-8)
      epscount <- epscount + 1
      print(epsbase)
    }
  }
  return(eps)
}

gradcheck <- function(sf,min=1e-5,seqlength=10,whichpars=NA,
  offdiags=FALSE,scale=TRUE,finite=FALSE){
  p <- sf$stanfit$rawest
  smf <- stan_reinitsf(sf$stanmodel,sf$standata)
  if(is.na(whichpars[1])) whichpars <- 1:length(p)
  
  func<-function(x) log_prob(smf,x)
  b=func(p)
  
  for(pi in whichpars){
    g <- list()
    gf=c()
    devs <- min*2^(1:seqlength)
    devs=sort(c(devs,-devs,0))
    for(i in seq_along(devs)){
      px <- p
      px[pi] <- p[pi] + devs[i]
      g[[i]] <- log_prob(smf,px,gradient=TRUE)
      if(finite) gf[i]<-(func(px)-b)/devs[i]
    }
    lp=sapply(g,function(x) x[1])
    plot(devs,lp,type='l',lwd=2,main=pi)
    abline(h=lp[devs==0],lty=2)
    abline(v=0,lty=2)
    gmat<- sapply(g,function(x) attributes(x)$gradient)
    if(finite) gfmat=sapply(gf,function(x) x*2)
    if(scale) gmat <- t(t(gmat)/devs)
    if(finite && scale) gfmat=t(t(gfmat)/devs)
    plot(devs,gmat[pi,],type='l',lwd=2,main=pi,ylab='grad')
    if(finite)  points(devs,gfmat,type='l',lwd=2,lty=2,col=2,main=pi,ylab='grad')
    abline(v=0,lty=2)
    abline(h=0,lty=2)
    if(offdiags){
      for(i in 1:length(p)){
        if(i != pi) plot(devs,gmat[i,],type='l',lwd=2,main=paste0(pi,'  ',i),ylab='grad')
        abline(v=0,lty=2)
        abline(h=0,lty=2)
      }
    }
  }
}


jacrandom <- function(grfunc, est, eps=1e-4,
  maxlpd=.5, minlpd=1e-5){
  
  fmax=grfunc(est)
  n <- max(length(est)*2,200)
  mcholtest=c()
  epsbase=.0001
  mcholtest <- FALSE
  i=0
  epsadapted <- FALSE
  pd<-rep(FALSE,n)
  H <- NA
  
  while(i < n && suppressWarnings(min(which(pd))) > (i-20)){
    # print(i)
    if(i==0){
      devs <- matrix(NA,nrow=n,ncol=length(est))
      covstore <- devs
      fg<-list()
      f<-c()
      grm <- devs
    }
    i=i+1
    # 
    accepted <- FALSE
    trycount <-0
    while(!accepted){
      trycount <- trycount + 1
      if(trycount > 20) stop('Errors encountered estimating Hessian')
      devs[i,] <- rnorm(length(est),0,eps)
      fg[[i]] <- try(grfunc((devs[i,])+est))
      if('try-error' %in% class(fg[[i]])) eps <- eps * .5 else accepted <- TRUE
    }
    
    
    grm[i,] <- attributes(fg[[i]])$gradient
    f[i] <- fg[[i]][1]
    
    # plot(f)
    bad=(f<quantile(f,.2,na.rm = TRUE) & (fmax-f) > maxlpd) | 
      (f>quantile(f, .95,na.rm = TRUE) & (fmax-f) < minlpd*.001)
    # if(i > length(est)+5) bad <- unique(c(1:(i-length(est)-5)),(bad))
    
    good <- c(1:i)[!bad]
    fl=f-min(f)
    fl <- sqrt(fl)
    
    if(length(good)>5){
      frange <- fmax-f[i]
      if(frange < minlpd) epsbase <- epsbase * 1.1
      if(frange > maxlpd) epsbase <- epsbase / 1.1
      eps=epsbase/abs(apply(devs,2,function(x) coefficients(lm(fl[good]~abs(x[good])))[-1]))
      # if(length(good)>length(est)){
      #   eps=epsbase/abs(coefficients(lm(fl[good]~abs(devs[good,])))[-1])
      #   }
      # print(epsbase)
      
      if(i==10 && !epsadapted){
        epsadapted <- TRUE
        i <- 0
      }
    }
    # print(round(epsbase,4))
    
    if(length(good) > length(est) && i > length(est)) {
      lf=t(apply(grm[good,],2,function(y) coefficients(lm((y) ~ (devs[good,])))[-1])) #without intercept
      #     lf=t(apply(grm[good,],2,function(y){
      #       apply(devs[good,],2,function(x) coefficients(lm((y) ~ x-1)))
      #       }))
      # print(i)
      
      lf=(lf+t(lf))/2
      # mc=MASS::ginv(-lf)
      # mc=as.matrix(Matrix::nearPD(mc)$mat)
      mchol = try(t(chol(solve(-lf))),silent=TRUE)
      pd[i] <- !'try-error' %in% class(mchol)
      if(pd[i]) H <- lf
      # covstore[i,] <- diag(mc)
    }
  }
  if(is.na(H[1])) H <- lf #return last non pos def if no good ones found
  return(H)
}


#' Sample more values from an optimized ctstanfit object
#'
#' @param fit fit object
#' @param nsamples number of extra samples desired
#' @param cores number of cores to use
#'
#' @return fit object with extra samples
#' @export
#'
#' @examples
#' if(w32chk()) newfit <- ctAddSamples(ctstantestfit, 10, 1)
ctAddSamples <- function(fit,nsamples,cores=2){
  mchol <- t(chol(fit$stanfit$cov))
  
  resamples <- matrix(unlist(lapply(1:nsamples,function(x){
    fit$stanfit$rawest + (mchol) %*% t(matrix(rnorm(length(fit$stanfit$rawest)),nrow=1))
  } )),byrow=TRUE,ncol=length(fit$stanfit$rawest))
  
  fit$stanfit$rawposterior <- rbind(fit$stanfit$rawposterior,resamples)
  
  fit$stanfit$transformedpars=stan_constrainsamples(sm = fit$stanmodel,
    standata = fit$standata,samples=fit$stanfit$rawposterior,
    savesubjectmatrices=as.logical(fit$standata$savesubjectmatrices),
    dokalman=as.logical(fit$standata$savesubjectmatrices),
    cores=cores)
  return(fit)
}

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
    library(ctsem)
    if(length(subindices) < length(unique(standata$subject))) standata <- standatact_specificsubjects(standata,subindices)
    if(!1 %in% subindices) standata$nopriors <- 1L
    if(1==99) sm=1
    g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    assign('smf',stan_reinitsf(sm,standata),pos = g)
    
    rm(standata,env=g)
    rm(sm,env=g)
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
#' if(w32chk()){
#'
#' sf <- stan_reinitsf(ctstantestfit$stanmodel,ctstantestfit$standata)
#' }
stan_reinitsf <- function(model, data,fast=FALSE){
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))
  
  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<-rstan::sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE,
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
#' if(w32chk()){
#'
#' d <- standatact_specificsubjects(ctstantestfit$standata, 1:2)
#' }
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
  longobjects <- c('subject','time','dokalmanrowsdata','nobs_y','ncont_y','nbinary_y','Y','tdpreds', 'whichobs_y','whichbinary_y','whichcont_y')
  return(longobjects)
}

standatatolong <- function(standata, origstructure=FALSE,ctm=NA){
  long <- lapply(standatalongobjects(),function(x) as.matrix(standata[[x]]))
  names(long) <- standatalongobjects()
  
  if(origstructure){
    colnames(long[['Y']]) <- ctm$manifestNames#colnames(standata$Y)
    long[['Y']][long[['Y']] %in% 99999] <- NA
    if(!is.na(ctm[1])){
      colnames(long[['subject']]) <- ctm$subjectIDname
      colnames(long[['time']]) <- ctm$timeName
    }
    longout <- data.frame(long[['subject']],long[['time']],long[['Y']])
    if(standata$ntdpred > 0){
      colnames(long[['tdpreds']]) <- colnames(standata$tdpreds)
      if(!is.na(ctm[1])) colnames(long[['tdpreds']]) <- ctm$TDpredNames
      longout <- cbind(longout,long[['tdpreds']])
    }
    if(standata$ntipred > 0){
      tipreds <- standata$tipredsdata[longout$id,,drop=FALSE]
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




stan_constrainsamples<-function(sm,standata, samples,cores=2, cl=NA,
  savescores=FALSE,
  savesubjectmatrices=TRUE,
  dokalman=TRUE,
  onlyfirstrow=ifelse(any(savesubjectmatrices,savescores),FALSE,TRUE),
  pcovn=500,
  quiet=FALSE){
  
  if(savesubjectmatrices && !dokalman){
    dokalman <- TRUE
    warning('savesubjectmatrices = TRUE requires dokalman=TRUE also!')
  }
  
  standata$savescores <- as.integer(savescores)
  standata$dokalman <- as.integer(dokalman)
  standata$savesubjectmatrices<-as.integer(savesubjectmatrices)
  if(onlyfirstrow) standata$dokalmanrowsdata <- as.integer(c(1,diff(standata$subject)))
  
  if(!quiet) message('Computing quantities for ', nrow(samples),' samples...')
  if(nrow(samples)==1) cores <- 1
  
  #create via eval due to parallel communication rubbish
  eval(parse(text="tparfunc <- function(x, parallel = TRUE){ 
    if(parallel) library(ctsem)
    smf <- stan_reinitsf(sm,standata)
    out <- list()
    for(li in 1:length(x)){
      out[[li]] <- try(rstan::constrain_pars(smf, upars=samples[x[li],]))
      if(any( #check for invalid values in transformed pars
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
    on.exit(try(parallel::stopCluster(cl),silent=TRUE),add = TRUE)
    # parallel::clusterExport(cl2, c('sm','standata','samples'),environment())
    # parallel::clusterApply(cl2,1:cores, function(x) library(ctsem))
  }
  if(cores > 1) parallel::clusterExport(cl, 'standata',environment())
  transformedpars <- try(flexlapply(cl, 
    split(1:nrow(samples), sort((1:nrow(samples))%%cores)),
    tparfunc,cores=cores,parallel=cores > 1))
  
  
  #fix this hack
  # 
  if(!is.null(transformedpars[[1]][[1]]$popmeans)) transformedpars=unlist(transformedpars,recursive = FALSE)
  est1=transformedpars[[1]]
  missingsamps <-sapply(transformedpars, function(x) 'try-error' %in% class(x))
  nasampscount <- sum(missingsamps) 
  
  
  if(nasampscount > 0) {
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
  
  transformedpars=try(tostanarray(flesh=matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = est1))
  
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
#' @param finishsamples Number of samples to draw (either from hessian
#' based covariance or posterior distribution) for final results computation.
#' @param finishmultiply Importance sampling stops once available samples reach \code{finishsamples * finishmultiply} , then the final samples are drawn
#' without replacement from this set.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent
#' importance sampling, at risk of truncating tails.
#' @param parsteps ordered list of vectors of integers denoting which parameters should begin fixed
#' at zero, and freed sequentially (by list order). Useful for complex models, e.g. keep all cross couplings fixed to zero 
#' as a first step, free them in second step. 
#' @param finitediff Either 'ask', TRUE, or FALSE. Whether to use the slow finite difference calculations
#' for the Hessian (used for confidence intervals) if other approaches do not give a positive definite result. 
#' @param chancethreshold drop iterations of importance sampling where any samples are chancethreshold times more likely to be drawn than expected.
#'
#' @return list containing fit elementsF
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp

stanoptimis <- function(standata, sm, init='random',initsd=.01,sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-12,
  decontrol=list(),
  stochastic = TRUE,
  nopriors=FALSE,carefulfit=TRUE,
  subsamplesize=.5,finitediff=FALSE,
  parsteps=c(),
  plot=FALSE,
  is=FALSE, isloopsize=1000, finishsamples=1000, tdf=10,chancethreshold=100,finishmultiply=5,
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
  
  notipredsfirstpass <- FALSE
  if(init[1] =='random'){# #remove tipreds for first pass
    notipredsfirstpass <- TRUE
    TIPREDEFFECTsetup <- standata$TIPREDEFFECTsetup
    standata$TIPREDEFFECTsetup[,] <- 0L
    standata$ntipredeffects <- 0L
  }
  
  if(standata$nindvarying > 0 && standata$intoverpop==0){ #detect subject level pars
    stochastic <- TRUE
    standata$doonesubject <- 1L
    message('Using combined stochastic gradient descent + mcmc')
    a1=standata$nparams+standata$nindvarying+
      (standata$nindvarying^2-standata$nindvarying)/2
    whichmcmcpars <- (a1+1):(a1+standata$nindvarying) #*standata$nsubjects)
  } else whichmcmcpars <- NA
  
  smf <- stan_reinitsf(sm,standata)
  npars=rstan::get_num_upars(smf)
  if(cores > 1 & !deoptim) rm(smf)
  
  if(stochastic=='auto' && npars > 50){
    message('> 50 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
    stochastic <- TRUE
  } else if(stochastic=='auto') stochastic <- FALSE
  if(length(parsteps)>0 && !stochastic){
    stochastic=TRUE
    message('Stochastic optimizer used for data driven parameter inclusion') 
  }
  
  parsets <- 1
  if(length(unique(standata$subject)) < cores){
    if(length(unique(standata$subject))==1){
      optimcores <- cores
      parsets <- cores
      if(cores >1) stochastic=TRUE
    } else  optimcores <- length(unique(standata$subject))
  } else optimcores <- cores
  
  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  
  message('Using ',cores,'/', parallel::detectCores(),' logical CPU cores')
  
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    
    # if(standata$TIpredAuto >0) parsteps <- c(parsteps,(npars:(npars-max(standata$TIPREDEFFECTsetup)+1)))
    
    if(all(init %in% 'random')){
      init <- rnorm(npars, 0, initsd)
      if(length(parsteps)>0) init[unlist(parsteps)] <- 0 
    }
    if(all(init == 0)) init <- rep(0,npars)
    if(length(init) != npars) init=init[1:npars]
    init[is.na(init)] <- 0
    if(notipredsfirstpass && standata$ntipredeffects > 0) init[length(init):(length(init)+1-standata$ntipredeffects)] <- 0
    
    
    
    
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
        if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(evaltime,2)),digits=14)
        if(gradnoise) attributes(out)$gradient <-attributes(out)$gradient *
          exp(rnorm(length(parm),0,1e-3))
        return(out)
      }
      
      if(optimcores==1) target = singletarget #we use this for importance sampling
      if(cores > 1){ #for parallelised computation after fitting, if only single subject
        clctsem=parallel::makeCluster(cores,useXDR=TRUE,type='PSOCK')
        on.exit(try({parallel::stopCluster(clctsem)},silent=TRUE),add=TRUE)
        parallel::clusterExport(clctsem,varlist = 'sm',envir = environment())
      }
      if(optimcores > 1) {
        
        #crazy trickery to avoid parallel communication pauses
        env <- new.env(parent = globalenv(),hash = TRUE)
        environment(parlp) <- env
        
        target<-function(parm,gradnoise=TRUE) {
          # whichframe <- which(sapply(lapply(sys.frames(),ls),function(x){ 'clctsem' %in% x}))
          a=Sys.time()
          if('list' %in% class(parm)){
            out2 <- parallel::parLapply(cl = clctsem,X = parm,parlp)
            bestset <- which(unlist(out2)==max(unlist(out2)))#sapply(out,max))
            # print(bestset)
            bestset <- bestset[1]
            parm <- parm[[bestset]]
            out <- out2[[bestset]]
            attributes(out)$bestset <- bestset
          } else {
            out2<- parallel::clusterCall(clctsem,parlp,parm)
            out <- try(sum(unlist(out2)),silent=TRUE)
            attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
          }
          b=Sys.time()
          b-a
          
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
          
          if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=parsets<2)
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
      
      nopriorsbak <- standata$nopriors
      # taylorheun <- standata$taylorheun
      finished <- FALSE
      
      if(carefulfit && !deoptim){ #init using priors
        message('Doing 1st pass with priors on reduced data set')
        
        standata$nopriors <- as.integer(0)
        # standata$taylorheun <- 1L
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
        if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standatasml,split=parsets<2)
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
            parsets=parsets,
            whichignore = unlist(parsteps),
            whichmcmcpars=whichmcmcpars,nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
            plot=plot,itertol=1e-1,deltatol=1e-2,maxiter=500)
        }
        
        standata$nopriors <- as.integer(nopriorsbak)
        # if(standata$ntipred ==0)
        # standata$taylorheun <- as.integer(taylorheun)
        standata$sdscale <- sdscale
        
        if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else init=optimfit$par
      } #end carefulfit init
      
      
      message('Optimizing...')
      if(length(parsteps) <1 && (standata$ntipred ==0 || notipredsfirstpass ==FALSE)) finished <- TRUE
      if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=parsets<2)
      if(optimcores==1) smf<-stan_reinitsf(sm,standata)
      
      if(!stochastic){
        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          method="L-BFGS",memory=100,
          # check_conv_every = 5,
          # try_newton_step=TRUE,
          line_search='Schmidt',c1=1e-4,c2=.9,step0='schmidt',ls_max_fn=999,
          abs_tol=NULL,grad_tol=NULL,
          rel_tol=tol*ifelse(!finished,100,1),
          step_tol=NULL,ginf_tol=NULL)
        optimfit$value = -optimfit$f
        init = optimfit$par
        if(is.infinite(bestfit)) stochastic<-TRUE
      }
      
      if(stochastic){#  || #carefulfit
        if(is.infinite(bestfit)) message('Switching to stochastic optimizer -- failed initialisation with bfgs')
        if(!stochastic && carefulfit) message('carefulfit = TRUE , so checking for improvements with stochastic optimizer')
        optimfit <- sgd(init, fitfunc = target,
          parsets=parsets,
          itertol = ifelse(!finished,1e-1,1e-3),
          deltatol=ifelse(!finished,1e-2,1e-5),
          # itertol = ifelse(standata$ntipred >0,1e-1,1e-3),
          # deltatol=ifelse(standata$ntipred >0,1e-2,1e-5),
          whichignore = unlist(parsteps),whichmcmcpars=whichmcmcpars,
          nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
          ndatapoints=standata$ndatapoints,plot=plot)
        if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else init=optimfit$par
      }
      
      # browser()
      if(length(parsteps) > 0){
        message('Freeing parameters...')
        finished <- FALSE
        while(!finished && length(parsteps)>0){
          if(length(parsteps)>1) parsteps <- parsteps[-1] else parsteps <- c()#parsteps[pstat> pcheck]
          
          optimfit <- sgd(init, fitfunc = target,
            parsets=parsets,
            itertol = 1e-2, deltatol= 1e-2,
            whichignore = unlist(parsteps),whichmcmcpars=whichmcmcpars,
            nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
            ndatapoints=standata$ndatapoints,plot=plot)
          
          if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else{
            finished <- TRUE
            init = optimfit$par
          }
        }
      }
      
      
      if(standata$ntipred > 0 && notipredsfirstpass){
        message('Including TI predictors...')
        # standata$tipredeffectscale <- tipredeffectscale
        standata$dokalmanrows[] <- 1L
        # init = optimfit$par
        optimbase = init#[-(length(optimfit$par):(length(optimfit$par)-max(TIPREDEFFECTsetup)+1))]
        finished <- FALSE
        found <- 0
        tia=standata$TIPREDEFFECTsetup
        while(!finished){
          if(!is.null(standata$TIpredAuto) && standata$TIpredAuto){
            message ('Looking for tipred effects...')
            oldtia=tia
            fit <- list(stanfit=list(rawest=init),standata=standata,stanmodel=sm)
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
              oldtiinit = tail(init,found)[oldtia[oldtia>0]]
              tiinit[tia[tia>0&oldtia>0]] <- oldtiinit
              standata$TIPREDEFFECTsetup <- array(as.integer(tia),dim=dim(tia))
              found <- max(tia)
              message('Found ',max(tia),' viable TIpred effects')
              standata$ntipredeffects <- as.integer(max(tia))
              init = head(init,length(optimbase))
              init = c(init,tiinit)
            } else {
              finished <- TRUE
              message('No further predictors found, finishing optimization...')
              # standata$taylorheun <- as.integer(taylorheun)
            }
          }
          # if(!finished){
          
          if(!standata$TIpredAuto){
            finished <- TRUE
            # standata$taylorheun <- as.integer(taylorheun)
            init = c(optimbase,rep(0,max(TIPREDEFFECTsetup)))
            standata$TIPREDEFFECTsetup[,] <- TIPREDEFFECTsetup
            standata$ntipredeffects <- as.integer(max(TIPREDEFFECTsetup))
          }
          
          if(optimcores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=parsets<2)
          if(optimcores==1) smf<-stan_reinitsf(sm,standata)
          if(!stochastic) {
            optimfit <- mize(init, fg=mizelpg, max_iter=99999,
              method="L-BFGS",memory=100,
              line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
              abs_tol=tol*ifelse(finished,1,10000),grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
            optimfit$value = -optimfit$f
          }
          if(stochastic && finished) optimfit <- sgd(init, fitfunc = target,
            parsets=parsets,
            whichignore = parsteps,whichmcmcpars=whichmcmcpars,
            nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
            plot=plot)
          if(stochastic && !finished) optimfit <- sgd(init, fitfunc = target,
            parsets=parsets,
            whichignore = parsteps,whichmcmcpars=whichmcmcpars,
            nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
            plot=plot,itertol=1e-1,deltatol=1e-2)
          if(length(parsteps)>0) init[-parsteps] = optimfit$par else init=optimfit$par
          # }
        }
        # smf <- stan_reinitsf(sm,standata)
        # 
        #   optimfit <- sgd(init, fitfunc = target,
        #     itertol = 1e-3, deltatol= 1e-3,
        #     whichignore = parsteps,whichmcmcpars=whichmcmcpars,
        #     nsubjects=ifelse(is.na(whichmcmcpars[1]),NA,standata$nsubjects),
        #     ndatapoints=standata$ndatapoints,plot=plot)
      } #end ti pred auto total loop
      
      
      
      bestfit <-optimfit$value
      est2=init #because init contains the fixed values #unconstrain_pars(smf, est1)
      if(length(parsteps)>0) est2[-parsteps] = optimfit$par else est2=optimfit$par
      
      if(standata$nindvarying > 0 && standata$intoverpop==0){ #recompose into single model
        standata$doonesubject <- 0L
        a1=standata$nparams+standata$nindvarying+
          (standata$nindvarying^2-standata$nindvarying)/2
        whichmcmcpars <- (a1+1):(a1+standata$nindvarying*standata$nsubjects)
        est3=est2[-whichmcmcpars]
        est3=c(est3,(optimfit$mcmcpars))
        est2=est3
        if(standata$ntipredeffects > 0) est3 <- c(est3,est2[(a1+1):length(a1)])
      }
      npars = length(est2)
    }
    
    if(!estonly){
      
      if(standata$nsubjects > cores){
      hesscl <- NA
      } else {
        hesscl <- clctsem
        if(cores > 1) parallelStanSetup(cl = clctsem,standata = standata,split=FALSE)#,split=parsets<2)
        if(cores==1) smf<-stan_reinitsf(sm,standata)
      }
      
      # base <- lapply(1:10,function(x) target(est2 + rnorm(length(est2),0,1e-12),
      #   gradnoise = FALSE))
      # basegrad <- mean(unlist(lapply(base,function(x) attributes(x)$gradient)))
      # base <- mean(unlist(base))
      
      # grmat<-function(pars,step=1e-1,lpdifmin=1e-1, 
      #   lpdifmax=3, direction=1,whichpars='all',gradmod=FALSE,parmod=TRUE){
      #   if('all' %in% whichpars) whichpars <- 1:length(pars)
      #   hessout <- flexsapply(cl = clctsem, cores = 1, whichpars, function(i) {
      #     stepsize <- step * direction
      #     if(parmod) stepsize <- min(step, abs(pars[i]*1e-4)) * direction
      #     colout <- NA
      #     dolpchecks <- TRUE #set to true to try the log prob checks again..
      #     # 
      #     while(any(is.na(colout)) && abs(stepsize) > 1e-20){
      #       stepsize <- stepsize * .1
      #       lpdifok<-FALSE
      #       lpdifcount <- 0
      #       lpdifdirection <- 0
      #       lpdifmultiplier <- 1
      #       # message('par',i)
      #       while(!lpdifok & lpdifcount < 15){
      #         # message(paste(i,'  col=',colout,'  lpdifmultiplier=',lpdifmultiplier, '  stepsize=',stepsize))
      #         lpdifok <- TRUE
      #         lpdifcount <- lpdifcount + 1
      #         uppars<-pars
      #         uppars[i]<-pars[i]+stepsize / ifelse(gradmod,(abs(basegrad[i])),1)
      #         
      #         suppressMessages(suppressWarnings(try({uplp<- target(uppars,gradnoise=FALSE)},silent=TRUE))) #try(smf$log_prob(upars=uppars,adjust_transform=TRUE,gradient=TRUE)) #lpg(uppars)
      #         # storedPars <<- cbind(storedPars,matrix(uppars))
      #         # storedLp <<- c(storedLp,uplp[1])
      #         
      #         if('try-error' %in% class(uplp)){
      #           lpdifok <- TRUE
      #           upgrad <- rep(NA,length(pars))
      #           dolpchecks <- FALSE
      #         } else{
      #           upgradnew= ((attributes(uplp)$gradient -basegrad) / stepsize)
      #           upgrad = upgradnew * ifelse(abs(base-uplp) < 2,.5,0) + 
      #             upgradnew*ifelse(lpdifcount==1 ||abs(base-uplp) < 2,.5,0)
      #           
      #           
      #           # print((  (upgrad[i] * (-abs(uppars[i]-pars[i]))) / (uplp[1]-bestfit[1]) ))
      #           # upgrad = upgrad / (  (upgrad[i] * (-abs(uppars[i]-pars[i]))) / (uplp[1]-bestfit[1]) ) #linearised upgrad
      #           
      #           if(dolpchecks){
      #             if(abs(base-uplp) > lpdifmax) {
      #               # print(abs(base-uplp))
      #               message(paste0('decreasing step for ', i))
      #               lpdifok <- FALSE
      #               if(lpdifdirection== 1) {
      #                 lpdifmultiplier = lpdifmultiplier * .5
      #               }
      #               stepsize = stepsize * (1e-2 * lpdifmultiplier)
      #               lpdifdirection <- -1
      #             }
      #             if(abs(base-uplp) < lpdifmin ) { #include sufficient gradient[i] checks #|| upgrad[i] < 1e-2
      #               # print(abs(base-uplp))
      #               message(paste0('increasing step for ', i))
      #               lpdifok <- FALSE
      #               if(lpdifdirection== -1) {
      #                 lpdifmultiplier = lpdifmultiplier * .5
      #               }
      #               stepsize = stepsize * (100 * lpdifmultiplier)
      #               lpdifdirection <- 1
      #             }
      #             if(any(is.na(c(uplp)))) stepsize = stepsize * .1
      #           }
      #         }
      #       }
      #       # print(abs(base-uplp))
      #       colout<- (upgrad)  * ifelse(gradmod,(abs(basegrad[i])),1)
      #     }
      #     rbind(colout)
      #   })
      #   return((hessout))# + t(hessout))/2)
      # }
      
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
        
        
        # whichhesspars = (1:npars)
        # if(standata$intoverpop == 0 && standata$nindvarying > 0) whichhesspars <- whichhesspars[-whichmcmcpars]
        # hessup= grmat(pars=est2,step=1e-3, direction=1,whichpars = whichhesspars,gradmod = FALSE)[whichhesspars,,drop=FALSE]
        # hessdown= grmat(pars=est2,step=1e-3, direction=1,whichpars = whichhesspars,gradmod = FALSE)[whichhesspars,,drop=FALSE]
        # hess = (hessup+hessdown)/2
        # cholcovup = try(suppressWarnings(t(chol(solve(-hessup)))),silent = TRUE)
        # cholcovdown = try(suppressWarnings(t(chol(solve(-hessdown)))),silent = TRUE)
        
        
        # grfunc <- function(x){
        #   if(length(parsteps)>0){
        #     pars <- est2
        #     pars[-parsteps] <- x
        #   } else pars <- x
        #   fg=target(pars)
        #   if(length(parsteps)>0){ 
        #     attributes(fg)$gradient <- attributes(fg)$gradient[-parsteps]
        #   }
        #   return(attributes(fg)$gradient)
        # }
        
        # fgfunc <- function(x){
        #   if(length(parsteps)>0){
        #     pars <- est2
        #     pars[-parsteps] <- x
        #   } else pars <- x
        #   fg=target(pars)
        #   if(length(parsteps)>0){ 
        #     attributes(fg)$gradient <- attributes(fg)$gradient[-parsteps]
        #   }
        #   # print(fg[1])
        #   if(fg[1] < -1e99) class(fg) <- c('try-error',class(fg))
        #   return(fg)
        # }
        
        # gfunc <- function(x) attributes(fgfunc(x))$gradient
        
        if(length(parsteps)>0) grinit= est2[-parsteps] else grinit = est2
        # hess <- numDeriv::jacobian(grfunc, grinit,method.args=list(eps=1e-3,r=3,v=10))
        # 
        # 
        # lpstore <- optimfit$lpstore
        # gstore <- t(optimfit$gstore[,(bestfit-lpstore) < .5])
        # parstore <- t(optimfit$parstore[,(bestfit-lpstore) < .5] - est2)
        # 
        # lmf <- apply(gstore,2,function(y) coefficients(lm(y ~ parstore-1)))
        
        # hess <- jacrandom(fgfunc, grinit)
        # 
        # eps <- findstepsize(grinit,fgfunc)
        
        jac<-function(pars,step=1e-3,whichpars='all',
          lpdifmin=1e-2,lpdifmax=5, cl=NA,verbose=1,directions=c(-1,1),parsteps=c()){
          if('all' %in% whichpars) whichpars <- 1:length(pars)
          base <- optimfit$value
          
          if(!is.na(hesscl[1]))  target <- function(x){ #if cluster is passed, create target function, otherwise use old target
            if(length(parsteps)>0){
              pars <- est2
              pars[-parsteps] <- x
            } else pars <- x
            fg=parlp(pars)
            if(length(parsteps)>0){ 
              attributes(fg)$gradient <- attributes(fg)$gradient[-parsteps]
            }
            # print(fg[1])
            if(fg[1] < -1e99) class(fg) <- c('try-error',class(fg))
            return(fg)
          }
          
          hessout <- flexsapply(cl = cl, cores = length(cl), whichpars, function(i){
          
          # if(is.na(cl[1])) fgfunc <- target
            
            # for(i in whichpars){
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
              while(!accepted && count < 15){
                if(count>8) stepsize=stepsize*-1 #is this good?
                stepchangemultiplier <- max(stepchangemultiplier,.11)
                count <- count + 1
                lp[[di]] <- target(pars+uppars*stepsize*directions[di])
                accepted <- !'try-error' %in% class(lp[[di]])
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
                # if(count > 1) print(paste0(count,'___',stepsize,'___',lpdiff))
              }
              step <<- stepsize
              steplist[[di]] <- stepsize
            }
            # 
            grad<- attributes(lp[[1]])$gradient / steplist[[di]] * directions[di]
            if(length(directions) > 1) grad <- (grad + attributes(lp[[2]])$gradient / (steplist[[di]]*-1))/2
            
            # hessout[,i] <-rbind(grad)
            return(grad)
          }
          ) #end flexapply
          # 
          out=(hessout+t(hessout))/2
          # attributes(out)$steplist=steplist
          # })
          return(out)
        }

        hess1 <- jac(pars = grinit,parsteps=parsteps,
          step = rep(1e-3,length(grinit)),cl=hesscl,verbose=verbose,directions=1)
        hess2 <- jac(pars = grinit,parsteps=parsteps,#fgfunc = fgfunc,
          step = rep(1e-3,length(grinit)),cl=hesscl,verbose=verbose,directions=-1)
        # 
        
        hess1good <- diag(hess1) < -1e-8
        hess2good <- diag(hess2) < -1e-8
        
        hess <- hess1 
        hess[,!hess1good & !hess2good] <- (hess1[,!hess2good & !hess1good] + 
            hess2[,!hess2good & !hess1good]) /2
        hess[,hess2good & hess1good] <- (hess1[,hess2good & hess1good] + 
            hess2[,hess2good & hess1good]) /2
        hess[,hess1good & !hess2good] <- hess1[,hess1good & !hess2good]
        hess[,hess2good & !hess1good] <- hess2[,hess2good & !hess1good]
        
        if(sum(hess1good)+sum(hess2good) < nrow(hess)*2){
          if(any(is.na(hess))) message ('Problems computing Hessian...')
          onesided <- c(which(!hess1good[hess2good]), which(!hess2good[hess1good]))
          if(length(onesided))  message ('One sided Hessian used for pars: ', onesided)
        }
        
        hess <- (t(hess)+hess)/2
        cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
        # 
        # if('try-error' %in% class(cholcov)){
        #   message('Trying harder...')
        #   hess <- numDeriv::jacobian(grfunc, grinit,method.args=list(r=5))
        #   cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
        # }
        # if('try-error' %in% class(cholcov)){
        # if(!'try-error' %in% class(cholcovup)){
        #   hess = hessup
        #   message('One sided hessian used for std error estimation')
        # } else if(!'try-error' %in% class(cholcovdown)){
        #   message('One sided hessian used for std error estimation')
        #   hess = hessdown
        # } else{
        
        # time2 <- Sys.time()-time1
        # esttime <- round(time2*npars/60,0)
        # continue <- ifelse(finitediff %in% FALSE,'N','Y')
        # if(time2 * npars > 5*60 && interactive() && finitediff=='ask' && length(parsteps)==0) continue <- readline(
        #   paste0('Hessian not positive definite, std errors may be invalid, try slower approach? Est ',esttime,' mins.  Y/N?'))
        # # if(!continue %in% c('n','N','no','No') || finitediff==TRUE){
        #   message('Trying finite differences Richardson Hessian')
        #   hess = numDeriv::hessian(target,est2,method.args=list(r=2))
        #   cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
        # }
        if('try-error' %in% class(cholcov) && !is) message('Approximate hessian used for std error estimation.')
        # }
        # }
        
        mcov=try(solve(-hess),silent=TRUE)
        if('try-error' %in% class(mcov)){
          mcov=MASS::ginv(-hess) 
          warning('Generalized inverse required for Hessian inversion -- standard errors not trustworthy')
        }
        
        # if(robust){
        #   message('Getting scores for robust std errors...')
        #   mcovnonrobust <-mcov
        #   fit <- list(standata=standata,stanmodel=sm) #this is hacky...
        #   fit$stanfit$rawest=est2
        #   
        #   score=as.matrix(ctsem:::scorecalc(fit,subjectsonly = F,returnsubjectlist = F,cores=cores))
        #   mcov=tcrossprod(mcov,MASS::ginv(crossprod(score)))
        # }
        # 
        mcovtmp=as.matrix(Matrix::nearPD(mcov,conv.norm.type = 'F')$mat)
        mcov <- diag(1e-10,npars)
        if(length(parsteps)>0) mcov[-parsteps,-parsteps] <- mcovtmp else mcov <- mcovtmp
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
        optimfit <- suppressWarnings(list(par=rstan::sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
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
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + (mchol) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        
      }
      if(is){
        message('Importance sampling...')
        if(cores > 1)  parallelStanSetup(cl=clctsem,standata,split=FALSE)
        targetsamples <- finishsamples * finishmultiply
        # message('Adaptive importance sampling, loop:')
        j <- 0
        while(nrow(samples) < targetsamples){
          j<- j+1
          if(j==1){
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
          } else {
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf))
          }
          
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          
          if(cores > 1) parallel::clusterExport(clctsem,c('samples'),envir = environment())
          
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
          
          # 
          target_dens[[j]] <- unlist(flexlapply(cl = clctsem, 
            split(1:isloopsize, sort((1:isloopsize)%%cores)), 
            function(x){
              # future(globals=c('x','samples','isloopsize'),expr={
              eval(parse(text='apply(tail(samples,isloopsize)[x,],1,parlp)'))
            },cores=cores))
          # 
          
          
          
          
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
    
    transformedpars=stan_constrainsamples(sm = sm,standata = standata,
      savesubjectmatrices = standata$savesubjectmatrices, 
      dokalman=standata$savesubjectmatrices,
      samples=resamples,cores=cores, cl=clctsem)
    
    if(cores > 1) {
      parallel::stopCluster(clctsem)
      smf <- stan_reinitsf(sm,standata)
    }
    
    transformedparsfull=stan_constrainsamples(sm = sm,standata = standata,
      savesubjectmatrices = TRUE, dokalman=TRUE,savescores = TRUE,
      samples=matrix(est2,nrow=1),cores=1, quiet = TRUE)
    
    
    
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
      transformedparsfull=transformedparsfull,
      standata=list(TIPREDEFFECTsetup=standata$TIPREDEFFECTsetup,ntipredeffects = standata$ntipredeffects),
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  
  if(estonly) {
    smf <- stan_reinitsf(sm,standata)
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  }
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}


