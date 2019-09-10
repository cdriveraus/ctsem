


parallelStanSetup <- function(cl, sm, standata,split=TRUE){
  cores <- length(cl)
  if(split) stansubjectindices <- split(unique(standata$subject),sort(unique(standata$subject) %% min(standata$nsubjects,cores)))
  if(!split) stansubjectindices <- lapply(1:cores,function(x) unique(standata$subject))
  if(!split && length(stansubjectindices) < cores){
    for(i in (length(stansubjectindices)+1):cores){
      stansubjectindices[[i]] <- NA
    }
  }
  parallel::clusterExport(cl,c('standata','sm'),envir = environment())
  parallel::clusterApply(cl,stansubjectindices,function(subindices) {
    # require(Rcpp)
    if(length(subindices) < length(unique(standata$subject))) standata <- standatact_specificsubjects(standata,subindices)
    # future(globals = c('sm','standata'),
    #   packages=c('ctsem','rstan'),
    # gc=FALSE,expr = {
    g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    assign('smfnode',stan_reinitsf(sm,standata,fast=TRUE),pos = g)
    
    parlp <- function(parm){
      out <- try(smfnode$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
      if(class(out)=='try-error') {
        out[1] <- -1e100
        attributes(out)$gradient <- rep(NaN, length(parm))
      }
      if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
      attributes(out)$gradient[is.nan(attributes(out)$gradient)] <- 
        rnorm(length(attributes(out)$gradient[is.nan(attributes(out)$gradient)]),0,100)
      return(out)
    }
    assign('parlp',parlp,pos = globalenv())
    
    NULL
  })
  NULL
}

parallelStanFunctionCreator <- function(cl, verbose){
  cores=length(cl)
  neglpgf<-function(parm) {
    a=Sys.time()
    out2<- parallel::clusterApply(cl,
      lapply(1:cores,function(x) parm),
      function(p){
        parlp(p)
      })
    b=Sys.time()
    
    out <- try(sum(unlist(out2)),silent=TRUE)
    attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
    
    if(class(out)=='try-error' || is.nan(out)) {
      out=-1e100
      attributes(out) <- list(gradient=rep(0,length(parm)))
    }
    assign('storedPars', parm,pos= sys.frame(4))
    if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',b-a))
    return(-out)
  }
  
  
  mizelpg=list(
    fg=function(pars){
      r=neglpgf(pars)
      r=list(fn=r[1],gr= -attributes(r)$gradient)
      return(r)
    },
    fn=neglpgf,
    gr=function(pars) -attributes(neglpgf(pars))$gradient
  )
  return(list(neglpgf=neglpgf,mizelpg=mizelpg))
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

flexlapply <- function(cl, X, fn,cores=1){
  if(cores > 1) parallel::parLapply(cl,X,fn) else lapply(X, fn)
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
#' standatact_specificsubjects(ctstantestfit$standata, 1:2)
standatact_specificsubjects <- function(standata, subjects,timestep=NA){
  long <- standatatolong(standata)
  long <- long[long$subject %in% subjects,]
  standatamerged <- standatalongremerge(long=long, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(long))
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
  nlong[,grep('^tdpreds.',colnames(nlong))] <- 0
  
  
  mlong <- rbind(long,nlong)
  mlong <- mlong[order(mlong$subject,mlong$time),]
  standatamerged <- standatalongremerge(long=mlong, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(mlong))
  return(standatamerged)
}




stan_constrainsamples<-function(sm,standata, samples,cores=2){
  smf <- stan_reinitsf(model = sm,data = standata)
  message('Computing quantities for ', nrow(samples),' samples...')
  est1=NA
  class(est1)<-'try-error'
  i=0
  while(i < nrow(samples) && class(est1)=='try-error'){
    i=i+1
    est1=try(constrain_pars(smf, upars=samples[i,]))
  }
  if(class(est1)=='try-error') stop('All samples generated errors! Respecify, try stochastic optimizer, try again?')
  
  cl2 <- parallel::makeCluster(cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl2),add = TRUE)
  parallel::clusterExport(cl2, c('sm','standata','samples','est1'),environment())
  
  transformedpars <- try(flexsapply(cl2, parallel::clusterSplit(cl2,1:nrow(samples)), function(x){ #could pass smaller samples
    # Sys.sleep(.1)
    if(!is.null(standata$savescores) && !standata$savescores) standata$dokalmanrows <- as.integer(c(1,standata$subject[-1] - standata$subject[-standata$ndatapoints]))
    smf <- stan_reinitsf(sm,standata)
    # Sys.sleep(.1)
    out <- list()
    skeleton=est1
    for(li in 1:length(x)){
      out[[li]] <- try(constrain_pars(smf, upars=samples[x[li],]))
      if(any(sapply(out[[li]], function(x) any(c(is.nan(x),is.infinite(x),is.na(x)))))) class(out[[li]]) <- c(class(out[[li]]),'try-error')
    }
    return(out)
  },cores=cores))
  transformedpars=unlist(transformedpars,recursive = FALSE)
  missingsamps <-sapply(transformedpars, function(x) 'try-error' %in% class(x))
  nasampscount <- sum(missingsamps) 
  
  transformedpars <- transformedpars[!missingsamps]
  nresamples <- nrow(samples) - nasampscount
  if(nasampscount > 0) {
    message(paste0(nasampscount,' NAs generated during final sampling of ', nrow(samples), '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
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


sgd <- function(init,fitfunc,ndatapoints,plotsgd=FALSE,stepbase=1e-4,gmeminit=ifelse(is.na(startnrows),.9,.9),gmemmax=.95,maxparchange = .5,
  startnrows=NA,gsmoothness = 1,roughnessmemory=.95,groughnesstarget=.5,lproughnesstarget=.1,
  gsmoothroughnesstarget=.2,
  warmuplength=50,
  minparchange=1e-20,maxiter=50000,nconvergeiter=20, itertol=1e-3, deltatol=1e-5){
  pars=init
  bestpars = pars
  maxpars=pars
  minpars=pars
  changepars=pars
  step=rep(stepbase,length(init))
  # parallelStanSetup(cores,sm,standata)
  
  g= -fitfunc(init) # try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
  if(class(g)=='try-error') {
    i = 0
    message('Problems initialising, trying random values...')
    while(i < 50 && class(g)=='try-error'){
      if(i %%5 == 0) init = rep(0,length(init))
      init=init+rnorm(length(init),0,abs(init)+ .1)
      g= -fitfunc(init) #try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE)
      i = i + 1
    }
  }
  g=sign(g)#*sqrt(abs(g))
  gsmooth=g
  oldg=g
  groughness = rep(groughnesstarget,length(g))
  gsmoothroughness = rep(gsmoothroughnesstarget,length(g))
  deltasmoothsq=.01
  lproughness=lproughnesstarget
  gmemory <- gmeminit
  oldgmemory <- gmemory
  oldlpdif <- 0
  lpdif <- 0
  maxlp <- -Inf
  i=0
  lp<-c()
  oldlp <- -Inf
  converged <- FALSE
  nrows <- ifelse(is.na(startnrows),ndatapoints, min(ndatapoints, startnrows))
  
  while(!converged && i < maxiter){
    print
    i = i + 1
    accepted <- FALSE
    lproughnesstarget2 = ifelse(nrows==ndatapoints,lproughnesstarget,.49)
    while(!accepted){
      newpars = bestpars
      delta =   step  * (gsmooth*(gsmoothness) + g*(1-(gsmoothness))) * exp((rnorm(length(g),0,.01)))
      delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
      newpars = newpars + delta
      
      #sub sampling
      # if(!is.na(startnrows) || (nrows!=ndatapoints)){
      #   subjects <- sample(1:ndatapoints,nrows,replace = FALSE)
      #   standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,ndatapoints) #
      #   parallelStanSetup(cores = cores,sm = sm,standata = standata) #could be improved for subsampling
      # }
      
      # lpg = try(smf$log_prob(upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
      lpg= -fitfunc(newpars)
      
      if(lpg > -1e99 && class(lpg) !='try-error' && !is.nan(lpg[1]) && all(!is.nan(attributes(lpg)$gradient))){
        accepted <- TRUE
      } 
      else {
        step <- step * .1
      }
      
      if(is.na(startnrows) && i < warmuplength && i > 1 && lpg[1] < lp[1]) {
        accepted <- FALSE
        if(plotsgd) print('not accepted!')
        # step = step * .5
        gsmooth=gsmooth*.5
      }
    }
    lp[i]=lpg[1]
    pars <- newpars
    
    oldg=g
    g=attributes(lpg)$gradient
    g=sign(g)*sqrt(abs(g))
    oldgsmooth = gsmooth
    gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
    gsmooth= gsmooth*gmemory2 + (1-gmemory2)*g#^2 #should it really be squared? sgd algorithms do so
    roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)
    
    # stdgdifold = (g-oldg) * step
    # stdgdifsmooth = (g-gsmooth) * step
    groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(g)!=sign(oldg))
    gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
    if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] >= (lp[i] + sd(tail(lp,min(i,3)))))
    
    # print(stdgdif)
    # step=exp(mean(log(step))+(.99*(log(step)-mean(log(step)))))
    # step[oldsigng == signg] = step[which(oldsigng == signg)] * sqrt(2-gmemory) #exp((1-gmemory)/2)
    # step[oldsigng != signg] = step[which(oldsigng != signg)] / sqrt(2-gmemory) #ifelse(nrows == ndatapoints, (2-gmemory),1.1) #1.2 #exp((1-gmemory)/2)
    
    signdifmod = step
    signdifmod[sign(oldg) == sign(g)] =  .1 #/ (1.5-inv_logit(abs(stdgdif[oldsigng == signg])))^4 #* (1/ ( ( (roughness*.05+.95)^2) ))
    signdifmod[sign(oldg) != sign(g)]  = -.1 #10* ((1.5-inv_logit(abs(stdgdifold[sign(oldg) != sign(g)])))-1) #* ( ( (roughness*.05+.95)^2) )
    signdifmod[is.nan(signdifmod)] <- .5 #oldstep[is.nan(step)] #because of overflow in some cases
    
    deltasmoothsq = deltasmoothsq * gmemory + (1-gmemory)*delta^2
    lproughnessmod= 2 * ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    gmemoryupd = min(gmemmax,max(.1,gmemory /  ( (1/(-mean(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ))
    gsmoothroughnessmod =  .5 *(( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    groughnessmod = .5 *( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    rmsstepmod = sqrt(abs(gsmooth+1e-7))/step -1 #like adagrad but with decaying gradient
    
    step = (step
      # + sqrt(deltasmoothsq)/abs(gsmooth) /2
      # + step
      + step*signdifmod #* min(sqrt(deltasmoothsq),1)
      + step*lproughnessmod
      + step* gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      + step* groughnessmod# * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    )
    
    # gmemory = gmemory * roughnessmemory2 + (1-roughnessmemory2)*gmemoryupd
    
    if(lp[i] >= max(lp)) {
      # step = step * sqrt(2-gmemory) #exp((1-gmemory)/8)
      if(i > warmuplength/2) {
        gsmooth[pars>maxpars | pars < minpars] <- gsmooth[pars>maxpars | pars < minpars]  + .2*delta[pars>maxpars | pars < minpars] /step[pars>maxpars | pars < minpars]
        step[pars>maxpars | pars < minpars] <- step[pars>maxpars | pars < minpars] * 2  #+ pars[pars>maxpars | pars < minpars]
        changepars=pars
        changepars[!(pars>maxpars | pars < minpars)] <- NA
        lproughness = lproughness * .9
      }
      # pars <- newpars
      bestpars <- pars
      
      maxpars[pars>maxpars] <-pars[pars>maxpars]
      minpars[pars<minpars] <-pars[pars<minpars]
      
      
    }
    if(i > 1 && lp[i] < lp[i-1]) {
      #slowly forget old max and mins, allow fast re exploration of space
      maxpars <- maxpars - (maxpars-minpars)*.01
      minpars <- minpars + (maxpars-minpars)*.01
      
      # step[sign(oldgsmooth) != signg] = .5 * step[sign(oldgsmooth) != signg]
      # gsmooth[sign(oldgsmooth) != signg] = .05 * g[sign(oldgsmooth) != signg]
      # gmemory=min(.8,gmemory)
      # signg <- oldsigng
      # gsmooth = oldgsmooth
      # pars <- bestpars
      # if(nrows == ndatapoints) {
      # gsmooth[sign(gsmooth) != signg] = gsmooth[sign(gsmooth) != signg] #* .5 + g[sign(gsmooth) != signg] * .5
      # step = step/ (2-gmemory)#step  / max( 1.5, (-10*(lp[i] - lp[i-1]) / sd(head(tail(lp,20),10)))) #exp((1-gmemory)/4)
      # } else {
      #   step = step / 1.06
      # }
    }
    # if(i %%10 ==0) gmemory = min(gmemory+.1,gmemmax)# * exp(mean(sign(diff(tail(lp,20)))))
    # if(i %%20 ==0) gmemory =  max(gmeminit, min(gmemmax, 1.6*(1-(log(sd(tail(lp,20)) ) -log(itertol)) / (log(sd(head(lp,20)))-log(itertol)))* (1-gmeminit) + gmeminit))
    
    if(i > 30 && i %% 20 == 0) {
      lpdif <- sum(diff(tail(lp,10)))
      oldlpdif <- sum(diff(head(tail(lp,10),20)))
      if(oldlpdif >= lpdif) gmemory <- oldgmemory
      proposal = gmemory*2-oldgmemory
      gmemory <- min(gmemmax, max(0, proposal + runif(1,-.05,.1)))
      oldgmemory <- gmemory
    }
    
    step[step > maxparchange] <- maxparchange
    step[step < minparchange] <- minparchange
    
    if(plotsgd){
      parbase=par(no.readonly=TRUE)
      on.exit(do.call(par,parbase),add=TRUE)
      par(mfrow=c(2,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
      plot(pars)
      points(changepars,pch=17,col='red')
      
      plot(log(step))
      
      plot(groughness,col='red',ylim=c(0,1))
      abline(h=mean(gsmoothroughness),col='blue',lty=2)
      abline(h=(gsmoothroughnesstarget),col='blue',lty=1,lwd=2)
      points(gsmoothroughness,ylim=c(0,1),col='blue')
      abline(h=mean(groughness),col='red',lty=2)
      # abline(h=(groughnesstarget),col='red',lty=1)
      
      abline(h=lproughnesstarget,lty=1,col='green')
      abline(h=lproughness, col='green',lty=2)
      
      plot(tail(log(-(lp-max(lp)-1)),500),type='l')
      plot(gsmooth,ylim= c(-max(abs(gsmooth)),max(abs(gsmooth))))
      
      matplot(cbind(signdifmod,gsmoothroughnessmod),col=c('black','blue'),pch=1,ylim=c(-1,1))
      points(groughnessmod,col='red')
      abline(h=lproughnessmod,col='green')
      
      # message(paste0('Iter = ',i, '   Best LP = ', max(lp),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
    }
    
    #check convergence
    if(i > 30){
      if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
      # print(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)))
      if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
      if(nrows < ndatapoints && (length(lp) - match(min(lp),lp)) > nconvergeiter) converged <- TRUE
    }
    if(converged & nrows != ndatapoints){
      converged <- FALSE
      nrows <- min(ndatapoints, nrows * 4)
      if(nrows > ndatapoints/2){
        nrows <- ndatapoints
        i=0
        lp=c()
      }
      message(paste0('nrows now ',nrows, ' out of ',ndatapoints),' total')
      
    }
  }
  return(list(itervalues = lp, value = max(lp),par=bestpars) )
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
#' @param plotsgd Logical. If TRUE, plot iteration details when using stochastic optimizer.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{DEoptim.control}.
#' @param tol objective tolerance.
#' @param nopriors logical. If TRUE, a nopriors integer is set to 1 (TRUE) in the standata object -- only has an effect if 
#' the stan model uses this value. 
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when nopriors=TRUE.
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
#' e <- extract(fit, permuted = TRUE)
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
  deoptim=FALSE, estonly=FALSE,tol=1e-14,
  decontrol=list(),
  stochastic = FALSE, #'auto',
  nopriors=FALSE,carefulfit=TRUE,
  plotsgd=FALSE,
  is=FALSE, isloopsize=1000, finishsamples=500, tdf=10,chancethreshold=100,finishmultiply=5,
  verbose=0,cores=2){
  
  if(!is.null(standata$verbose)) {
    if(standata$verbose > 1) standata$verbose=as.integer(verbose) else standata$verbose=0L
  }
  standata$nopriors=as.integer(nopriors)
  
  if(is.null(decontrol$steptol)) decontrol$steptol=5
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  if(is.null(init)) init <- 'random'
  if(init[1] !='random') carefulfit <- FALSE

  # Preserve execution plan
  # oplan <- future::plan()
  # on.exit(future::plan(oplan), add = TRUE)
  
  # if(cores > 1) plan(multisession,persistent=TRUE,workers=cores,globals=FALSE,gc=FALSE) else {
  #   plan(transparent,local=FALSE,globals=FALSE)
  # }
  if(length(unique(standata$subject)) < cores) optimcores <- length(unique(standata$subject))  else optimcores <- cores
  
  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    
    smf <- stan_reinitsf(sm,standata)
    npars=rstan::get_num_upars(smf)
    
    if(all(init %in% 'random')) init <- rnorm(npars, 0, initsd)
    if(all(init == 0)) init <- rep(0,npars)
    
    if(is.na(sampleinit[1])){
      
      if(deoptim){ #init with DE
        message('Using differential evolution for initialization')
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
          
          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]
          
          lp2 = function(parm) {
            out<-try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
            if(class(out)=='try-error') {
              out=-1e200
            }
            return(-out)
          }
          
          deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
          deinit[2,] <- rnorm(npars,0,.0002)
          if(length(init)>1 & try2) {
            deinit[1,] <- unconstrain_pars(smf,init)
            if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
          }
          decontrollist$initialpop=deinit
          decontrollist$NP = NP
          optimfitde <- suppressWarnings(DEoptim::DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
            control = decontrollist))
          # init=constrain_pars(object = smf,optimfitde$optim$bestmem)
          init=optimfitde$optim$bestmem
          bestfit = -optimfitde$optim$bestval
        } else stop(paste0('use install.packages(\"DEoptim\") to use deoptim')) #end require deoptim
      }
      
      # if(optimcores > 1){ #check if parallel helps optimizing
      #   a <- Sys.time()
      #   # a=log_prob(smf,init)
      #   if( (Sys.time() - a) < .2) optimcores <- 1
      # }
      
      
      
      storedPars <- c()
      smfnode <- NULL #global variables issue...
      parlp <- NULL 
      optimfinished <- FALSE
      on.exit({
        if(!optimfinished){
          message('Optimization cancelled -- restart from current point by including this argument:')
          message((paste0(c('init = c(',   paste0(round(storedPars,5),collapse=', '), ')'    ))))
        }},add=TRUE)
      if(cores > 1) {
        cl<-parallel::makeCluster(cores,type='PSOCK',useXDR=FALSE)
        on.exit(parallel::stopCluster(cl),add = TRUE)
      }
      if(optimcores > 1) {
        fitfuncs <- parallelStanFunctionCreator(cl=cl,verbose = verbose)  
        mizelpg <- fitfuncs$mizelpg
        neglpgf=fitfuncs$neglpgf
      }
      if(optimcores == 1) {
        neglpgf<-function(parm) {
          a=Sys.time()
          out<- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
          
          if(class(out)=='try-error' || is.nan(out)) {
            out=-1e100
            attributes(out) <- list(gradient=rep(0,length(parm)))
          }
          storedPars <<- parm
          b=Sys.time()
          evaltime <- b-a
          if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(evaltime,2)))
          return(-out)
        }
        
        mizelpg=list( #list also created in parallel function creator
          fg=function(pars){
            r=neglpgf(pars)
            r=list(fn=r[1],gr= -attributes(r)$gradient)
            return(r)
          },
          fn=neglpgf,
          gr=function(pars) -attributes(neglpgf(pars))$gradient
        )
        
      }
      
      
      
      
      if(stochastic=='auto' && npars > 100){
        message('> 100 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE
      
      if(carefulfit && !deoptim & standata$nopriors == 1 ){ #init using priors
        message('carefulfit=TRUE , doing 1st pass with priors')
        standata$nopriors <- as.integer(0)
        tipredscale <- standata$tipredscale
        standata$tipredscale <- .0001
        if(optimcores > 1) parallelStanSetup(cl = cl,sm = sm,standata = standata)
        if(optimcores==1) smf<-stan_reinitsf(sm,standata)
        if(!stochastic) {
          optimfit <- mize(init, fg=mizelpg, max_iter=99999,
            method="L-BFGS",memory=100,
            line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
            abs_tol=tol*1e8,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
          optimfit$value = optimfit$f
          
        }
        if(stochastic) optimfit <- sgd(init, fitfunc = neglpgf,ndatapoints=standata$ndatapoints,plotsgd=plotsgd,itertol=1,maxiter=300)
        standata$nopriors <- as.integer(1)
        # smf<-stan_reinitsf(sm,standata)
        init = optimfit$par #+ rnorm(length(optimfit$par),0,abs(init/8)+1e-3)#rstan::constrain_pars(object = smf, optimfit$par)
        standata$tipredscale <- tipredscale
      } #end carefulfit
      
      
      message('Optimizing...')
      if(optimcores > 1) parallelStanSetup(cl = cl,sm = sm,standata = standata)
      if(optimcores==1) smf<-stan_reinitsf(sm,standata)
      
      if(!stochastic) {
        
        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          method="L-BFGS",memory=100,
          line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
          abs_tol=tol,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
        optimfit$value = optimfit$f
        init = optimfit$par
        bestfit <- -optimfit$value
        optimfit$value <- -optimfit$value
      }
      
      if(stochastic || is.infinite(bestfit) || carefulfit){
        if(is.infinite(bestfit)) {
          message('Switching to stochastic optimizer -- failed initialisation with bfgs')
        }
        if(!stochastic && carefulfit) message('carefulfit = TRUE , so checking for improvements with stochastic optimizer')
        optimfit <- sgd(init, fitfunc = neglpgf,ndatapoints=standata$ndatapoints,plotsgd=plotsgd)
        bestfit <-optimfit$value
      }
      
      est2=optimfit$par #unconstrain_pars(smf, est1)
    }
    
    
    if(!estonly){
      
      grmat<-function(pars,step=1e-5,lpdifmin=1e-8, lpdifmax=1e-3){
        hessout <- sapply(1:length(pars), function(i) {
          # smf <- stan_reinitsf(sm,standata,fast=TRUE)
          stepsize <- step *10
          colout <- NA
          dolpchecks <- FALSE #set to true to try the log prob checks again...
          while(any(is.na(colout)) && stepsize > 1e-16){
            stepsize <- stepsize * .1
            lpdifok<-FALSE
            lpdifcount <- 0
            lpdifdirection <- 0
            lpdifmultiplier <- 1
            # message('par',i)
            while(!lpdifok & lpdifcount < 15){
              # message(stepsize)
              # message(paste(i,'  col=',colout,'  lpdifmultiplier=',lpdifmultiplier, '  stepsize=',stepsize))
              lpdifok <- TRUE
              lpdifcount <- lpdifcount + 1
              uppars<-pars
              downpars<-pars
              uppars[i]<-pars[i]+stepsize
              downpars[i]<-pars[i]-stepsize
              uplp= -neglpgf(uppars) #try(smf$log_prob(upars=uppars,adjust_transform=TRUE,gradient=TRUE)) #lpg(uppars)
              downlp = -neglpgf(downpars) #try(smf$log_prob(upars=downpars,adjust_transform=TRUE,gradient=TRUE)) #lpg(downpars)
              if(class(uplp)=='try-error' || class(downlp)=='try-error'){
                lpdifok <- TRUE
                upgrad <- rep(NA,length(pars))
                downgrad <- rep(NA,length(pars))
                dolpchecks <- FALSE
                # if(stepsize < 1e-12) 
              } else{
                upgrad= attributes(uplp)$gradient
                downgrad = attributes(downlp)$gradient
                
                if(dolpchecks){
                  if(abs(uplp-downlp) > lpdifmax) {
                    # message(paste0('decreasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== 1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * 1e-2 * lpdifmultiplier
                    lpdifdirection <- -1
                  }
                  if(abs(uplp-downlp) < lpdifmin) {
                    # 
                    # message(paste0('increasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== -1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * 100 * lpdifmultiplier
                    lpdifdirection <- 1
                  }
                  if(any(is.na(c(uplp,downlp)))) stepsize = stepsize * .1
                }
              }
            }
            # hessout[i,]<- (upgrad-downgrad) /stepsize/2
            colout<- (upgrad-downgrad) /stepsize/2
            # print(colout)
          }
          rbind(colout)
        })
        return(t(hessout))
      }
      
      # 
      # hess1s<-function(pars,direction=1,step=1e-5,lpdifmin=1e-6, lpdifmax=1e-1){
      #   hessout<-matrix(NA,nrow=length(pars),ncol=length(pars))
      #   bestlp=lpg(pars)
      #   basegrad=gradout
      #   for(i in 1:length(pars)){
      #     stepsize <- step
      #     lpdifok<-FALSE
      #     lpdifcount <- 0
      #     lpdifdirection <- 0
      #     lpdifmultiplier <- 1
      #     while(!lpdifok & lpdifcount < 15){
      #       lpdifok <- TRUE
      #       lpdifcount <- lpdifcount + 1
      #       uppars<-pars
      #       uppars[i]<-pars[i]+stepsize*direction
      #       uplp=lpg(uppars)
      #       upgrad=gradout
      #       if(abs(uplp-bestlp) > lpdifmax) {
      #         # message(paste0('decreasing step for ', i))
      #         lpdifok <- FALSE
      #         if(lpdifdirection== 1) {
      #           lpdifmultiplier = lpdifmultiplier * .5
      #         }
      #         stepsize = stepsize * 1e-2 * lpdifmultiplier
      #         lpdifdirection <- -1
      #       }
      #       if(abs(uplp-bestlp) < lpdifmin) {
      #         # message(paste0('increasing step for ', i))
      #         lpdifok <- FALSE
      #         if(lpdifdirection== -1) {
      #           lpdifmultiplier = lpdifmultiplier * .5
      #         }
      #         stepsize = stepsize * 100 * lpdifmultiplier
      #         lpdifdirection <- 1
      #       }
      #       hessout[i,]<- (upgrad-basegrad) /stepsize*direction
      #     }
      #     
      #     # }
      #   }
      #   return(t(hessout))
      # }
      
      
      # A more numerically stable way of calculating log( sum( exp( x ))) Source:
      # http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
      log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
      }
      
      
      if(is.na(sampleinit[1])){
        # hessup=hess1s(pars = est2,direction = 1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hessdown=hess1s(pars = est2,direction = -1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hess=(hessup+hessdown)/2
        
        hess=grmat(pars=est2,step=1e-12)
        for(ri in 1:nrow(hess)){
          for(ci in 1:ncol(hess)){
            if(is.na(hess[ri,ci])) hess[ri,ci] <- hess[ci,ri]
          }}
        if(any(is.na(hess))) warning(paste0('Hessian could not be computed for pars: ', paste0(which(apply(hess,1,function(x) any(is.na(x)))),collapse=', '), ' -- standard errors will be nonsense, model adjustment may be needed.',collapse=''))
        diag(hess)[is.na(diag(hess))]<- -1
        hess[is.na(hess)] <- 0
        hess = ((hess) + t(hess))/2
        # neghesschol = try(chol(-hess),silent=TRUE)
        
        mchol=try(t(chol(solve(-hess))),silent=TRUE)
        if(class(mchol)=='try-error') {
          message('Hessian not positive-definite so approximating, treat SE\'s with caution, consider respecification / priors.')
          npd <- TRUE
        } else npd <- FALSE
        # if(class(mchol)=='try-error') {
        mcov=MASS::ginv(-hess) #-optimfit$hessian)
        mcov=as.matrix(Matrix::nearPD(mcov,conv.norm.type = 'F')$mat)
      }
      
      if(!is.na(sampleinit[1])){
        mcov = cov(sampleinit)*1.5+diag(1e-6,ncol(sampleinit))
        est2 = apply(sampleinit,2,mean)
        bestfit = 9e100
        optimfit <- suppressWarnings(list(par=sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
      }
      
      mcovl <- list()
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
          delta[[1]] + t(chol(mcovl[[1]])) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        message('Getting ',finishsamples,' samples from Hessian for interval estimates...')
      }
      
      if(is){
        targetsamples <- finishsamples * finishmultiply
        # message('Adaptive importance sampling, loop:')
        j <- 0
        while(nrow(samples) < targetsamples){
          j<- j+1
          if(j==1){
            # if(!npd)
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            
          } else {
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf))
            # samples <- rbind(samples, MASS::mvrnorm(isloopsize, delta[[j]],mcovl[[j]]))
          }
          # if(j > 1 || !npd)
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          # prop_dens <- mvtnorm::dmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]],log = TRUE)
          # prop_dens <- ctdmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]])
          
          
          
          standata$verbose <- 0L
          if(cores > 1){
          parallelStanSetup(cl=cl,sm,standata,split=FALSE)
          neglpgf <- parallelStanFunctionCreator(cl=cl,verbose = 0)$neglpgf
          parallel::clusterExport(cl,'samples',envir = environment())
          }
          target_dens[[j]] <- flexsapply(cl = cl, split(1:isloopsize, sort((1:isloopsize)%%cores)), function(x){
            # future(globals=c('x','samples','isloopsize'),expr={
              eval(parse(text='apply(tail(samples,isloopsize)[x,],1,parlp)'))
            },cores=cores)
          
          
          
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
          
          if(nrow(samples) >= targetsamples) nresamples <- finishsamples else nresamples = min(5000,length(samples)/5)
          
          
          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)
          
          newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          counter <- 1
          while(counter < 30 &&
              length(unique(sample(x=1:length(newsampleprob),size=length(newsampleprob),prob=newsampleprob,replace=TRUE))
              ) < 20) {
            if(counter==1) message ('Sampling problematic -- trying to recover... ')
            counter = counter + 1
            if(counter == 30) stop('Importance sampling failed -- either posterior mode not found, or mode is inadequate starting point for sampling')
            weighted_dens <- weighted_dens /2
            newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
            # plot(newsampleprob)
          }
          sample_prob <- c(sample_prob,newsampleprob) #sum to 1 for each iteration, normalise later
          
          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
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
          
          #check max resample probability and drop earlier samples if too high
          dropuntil <- ceiling(max(c(0,which(sample_prob > (chancethreshold / isloopsize)) / isloopsize),na.rm=TRUE))*isloopsize
          if((isloopsize - dropuntil) > isloopsize) dropuntil <- dropuntil -isloopsize
          # 
          if(length(unique(resample_i)) < 200) dropuntil <- 0 
          if(nrow(samples) < isloopsize *2) dropuntil <- 0
          
          if(dropuntil > 0){
            resamples <- resamples[-(0:dropuntil),,drop=FALSE]
            sample_prob <- sample_prob[-(0:dropuntil)]
            samples <- samples[-(0:dropuntil),,drop=FALSE]
          }
          
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
    
    transformedpars=stan_constrainsamples(sm = sm,standata = standata,samples=resamples,cores=cores)
    
    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if(class(sds)=='try-error') sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, upars=lest)),
      unlist(constrain_pars(smf, upars= est2)),
      unlist(constrain_pars(smf, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2, rawposterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}

