standatact_specificsubjects <- function(standata, subjects){
  standata$dokalmanrows <- as.integer(standata$dokalmanrows * (standata$subject %in% subjects))
  standata2=standata
  standata2$dokalmanpriormodifier <- sum(standata$dokalmanrows)/standata$ndatapoints
  for(pi in c('time','subject','dokalmanrows','nobs_y','nbinary_y','ncont_y')){
    standata2[[pi]] <- standata2[[pi]][standata$dokalmanrows == 1]
  }
  for(pi in c('whichcont_y','Y','whichbinary_y','whichobs_y','tdpreds')){
    standata2[[pi]] <- standata2[[pi]][standata$dokalmanrows == 1,,drop=FALSE]
  }
  standata2$ndatapoints=as.integer(nrow(standata2$Y))
  return(standata2)
}  

parlp <- function(parm,subjects=NA){
  if(!is.na(subjects[1])) standata <- standatact_specificsubjects(standata,subjects) 
  smf<-stan_reinitsf(sm,standata,fast=TRUE) #list[[corei]]
  try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
}


stan_constrainsamples<-function(sm,standata, samples,cores=2){
  smfull <- stan_reinitsf(model = sm,data = standata)
  message('Computing quantities for ', samples,' samples...')
  est1=NA
  class(est1)<-'try-error'
  i=0
  while(class(est1)=='try-error'){
    i=i+1
    est1=try(constrain_pars(smfull, upars=samples[i,]),silent=TRUE)
  }
  if(class(est1)=='try-error') stop('All samples generated errors! Respecify, try stochastic optimizer, try again?')

  cl2 <- parallel::makeCluster(cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl2))
  parallel::clusterExport(cl2, c('sm','standata','samples','est1'),environment())
  parallel::clusterApply(cl2,1:cores,function(x) require(ctsem))
  # target_dens <- c(target_dens,
  
  transformedpars <- try(parallel::parLapply(cl2, parallel::clusterSplit(cl2,1:nrow(samples)), function(x){ #could pass smaller samples
    # Sys.sleep(.1)
    if(!is.null(standata$savescores) && !standata$savescores) standata$dokalmanrows <- as.integer(c(1,standata$subject[-1] - standata$subject[-standata$ndatapoints]))
    smfull <- stan_reinitsf(sm,standata)
    # Sys.sleep(.1)
    out <- list()
    skeleton=est1
    for(li in 1:length(x)){
      out[[li]] <- try(constrain_pars(smfull, upars=samples[x[li],]))
      if(any(sapply(out[[li]], function(x) any(c(is.nan(x),is.infinite(x),is.na(x)))))) class(out[[li]]) <- c(class(out[[li]]),'try-error')
    }
    return(out)
  }))
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
#' @param nopriors logical.f If TRUE, any priors are disabled -- sometimes desirable for optimization. 
#' @param startnrows Either NA to ignore, or an integer specifying the number of subjects to start the
#' stochastic optimization with. This will be increased to the full amount by convergence time.
#' @param tol objective tolerance. 
#' @param cores Number of cpu cores to use.
#' @param isloops Number of iterations of adaptive importance sampling to perform after optimization.
#' @param isloopsize Number of samples per iteration of importance sampling.
#' @param finishsamples Number of samples to use for final results of importance sampling.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent 
#' importance sampling, at risk of truncating tails.
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when nopriors=TRUE. 
#'
#' @return ctStanFit object
#' @importFrom mize mize
#' @importFrom Matrix bdiag
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @export
#'
#' @examples
#' \donttest{
#'  sunspots<-sunspot.year
#'  sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#'  id <- 1
#'  time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#' 
#' #setup model
#'  ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1,
#'   manifestNames='sunspots',
#'   latentNames=c('ss_level', 'ss_velocity'),
#'    LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#'    DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
#'    MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#'    MANIFESTVAR=diag(0,1),
#'    CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'    T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
#'    DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#' 
#'  ssmodel$pars$indvarying<-FALSE #Because single subject
#'  ssmodel$pars$offset[14]<- 44 #Because not mean centered
#'  ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality
#' 
#' #fit using optimization without importance sampling
#' ssfit <- ctStanFit(datalong[1:50,], #limited data for example
#'   ssmodel, optimize=TRUE,optimcontrol=list(isloops=0,finishsamples=50))
#' 
#' #output
#' summary(ssfit)
#' }
optimstan <- function(standata, sm, init='random',initsd=.01,sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-14,
  decontrol=list(),
  stochastic = FALSE, #'auto',
  startnrows=NA,
  plotsgd=FALSE,
  isloops=0, isloopsize=1000, finishsamples=500, tdf=50,carefulfit=TRUE,
  verbose=0,nopriors=FALSE,cores=2){
  
  standata$verbose=as.integer(verbose)
  standata$nopriors=as.integer(nopriors)
  
  if(is.null(decontrol$steptol)) decontrol$steptol=5 
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  
  
  
  sgd <- function(init,stepbase=1e-4,gmeminit=ifelse(is.na(startnrows),.9,.9),gmemmax=.95,maxparchange = .5,
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
    g=try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
    if(class(g)=='try-error') {
      i = 0
      message('Problems initialising, trying random values...')
      while(i < 50 && class(g)=='try-error'){
        if(i %%5 == 0) init = rep(0,length(init))
        init=init+rnorm(length(init),0,abs(init)+ .1)
        g=try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE)
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
    nrows <- ifelse(is.na(startnrows),standata$ndatapoints, min(standata$ndatapoints, startnrows))
    
    while(!converged && i < maxiter){
      print
      i = i + 1
      accepted <- FALSE
      lproughnesstarget2 = ifelse(nrows==standata$ndatapoints,lproughnesstarget,.49)
      while(!accepted){
        newpars = bestpars
        delta =   step  * (gsmooth*(gsmoothness) + g*(1-(gsmoothness))) * exp((rnorm(length(g),0,.01)))
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        newpars = newpars + delta
        
        #sub sampling
        if(!is.na(startnrows) || (nrows!=standata$ndatapoints)){
          subjects <- sample(1:standata$ndatapoints,nrows,replace = FALSE)
          standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,standata$ndatapoints) #
          # smf<-stan_reinitsf(sm,standata,fast=TRUE)
        }
        
        # lpg = try(smf$log_prob(upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        lpg= -neglpgf(newpars)
        if(class(lpg) !='try-error' && !is.nan(lpg[1]) && all(!is.nan(attributes(lpg)$gradient))) accepted <- TRUE else step <- step * .1
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
      
      stdgdifold = (g-oldg) * step
      stdgdifsmooth = (g-gsmooth) * step
      groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(g)!=sign(oldg))
      gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
      if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] >= (lp[i] + sd(tail(lp,min(i,3)))))
      
      # print(stdgdif)
      # step=exp(mean(log(step))+(.99*(log(step)-mean(log(step)))))
      # step[oldsigng == signg] = step[which(oldsigng == signg)] * sqrt(2-gmemory) #exp((1-gmemory)/2)
      # step[oldsigng != signg] = step[which(oldsigng != signg)] / sqrt(2-gmemory) #ifelse(nrows == standata$ndatapoints, (2-gmemory),1.1) #1.2 #exp((1-gmemory)/2)
      
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
          step[pars>maxpars | pars < minpars] <- step[pars>maxpars | pars < minpars] * 1.5  #+ pars[pars>maxpars | pars < minpars]
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
        # step[sign(oldgsmooth) != signg] = .5 * step[sign(oldgsmooth) != signg]
        # gsmooth[sign(oldgsmooth) != signg] = .05 * g[sign(oldgsmooth) != signg]
        # gmemory=min(.8,gmemory)
        # signg <- oldsigng
        # gsmooth = oldgsmooth
        # pars <- bestpars
        # if(nrows == standata$ndatapoints) {
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
        if(nrows < standata$ndatapoints && (length(lp) - match(min(lp),lp)) > nconvergeiter) converged <- TRUE
      }
      if(converged & nrows != standata$ndatapoints){
        converged <- FALSE
        nrows <- min(standata$ndatapoints, nrows * 4)
        if(nrows > standata$ndatapoints/2){
          nrows <- standata$ndatapoints
          i=0
          lp=c()
        }
        message(paste0('nrows now ',nrows, ' out of ',standata$ndatapoints),' total')
        
      }
    }
    return(list(itervalues = lp, value = max(lp),par=bestpars) )
  }
  
  
  
  
  
  
  
  
  message('Optimizing...')

  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    
    smfull <- stan_reinitsf(sm,standata)
    smf <- stan_reinitsf(sm,standata,fast=TRUE)
    npars=smf$num_pars_unconstrained()
    
    if(all(init %in% 'random')) init <- rnorm(npars, 0, initsd)
    if(all(init == 0)) init <- rep(0,npars)
    
    if(is.na(sampleinit[1])){
      
      if(deoptim){ #init with DE
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
          
          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]
          
          lp2 = function(parm) {
            out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
            if(class(out)=='try-error') {
              out=-1e200
            }
            return(-out)
          }
          
          deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
          deinit[2,] <- rnorm(npars,0,.0002)
          if(length(init)>1 & try2) {
            deinit[1,] <- unconstrain_pars(smfull,init)
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
      
      gradout <- c()
      bestlp <- -Inf
      
      cl <- parallel::makeCluster(cores, type = "PSOCK") #should vary this depending on os for memory management
      parallel::clusterExport(cl = cl, varlist = c('sm','standata','parlp','standatact_specificsubjects'),envir = environment())
      parallel::clusterApply(cl,1:cores,function(x) library(ctsem))
      on.exit(parallel::stopCluster(cl))
      
      if(cores > 1){
        a=Sys.time()
        evaltime <- smf$log_prob(upars=init, adjust_transform=TRUE, gradient=TRUE)
        b=Sys.time()
        evaltime <- b-a
        # print(evaltime)
      }
      
      neglpgf<-function(parm) {

        if(cores > 1 && evaltime > .1){
          
          out2 <- parallel::clusterApply(cl, 
            split(1:standata$nsubjects,sort(1:standata$nsubjects %% min(standata$nsubjects,cores))), function(subjects) parlp(parm,subjects))
          out <- try(sum(unlist(out2)),silent=TRUE)
          if(standata$verbose > 0) print(out)
          attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
        } else {
          out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        }
        
        if(class(out)=='try-error' || is.nan(out)) {
          out[1]=-Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          if(out[1] > bestlp) {
            bestlp <<- out[1]
            gradout <<- attributes(out)$gradient
          }
        }
        return(-out)
      }
      
      grffromlp<-function(parm) {
        return(-gradout)
      }
      
      parbase=par()
      
      
      # require(mize)
      mizelpg=list(
        fg=function(pars){
          r=neglpgf(pars)
          r=list(fn=r[1],gr= -attributes(r)$gradient)
          return(r)
        },
        fn=neglpgf,
        gr=function(pars) -attributes(neglpgf(pars))$gradient
      )
      
      
      if(stochastic=='auto' && npars > 100){
        message('> 100 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE
      
      if(carefulfit && !deoptim & standata$nopriors == 1 ){ #init using priors
        standata$nopriors <- as.integer(0)
        smf<-stan_reinitsf(sm,standata,fast=TRUE)
        parallel::clusterExport(cl = cl, varlist = c('smf','standata'),envir = environment())
        if(!stochastic) {
          
          # optimfit <- ucminf(init,fn = neglpgf,gr = grffromlp,control=list(xtol=tol*1e8,maxeval=300))
          optimfit <- mize(init, fg=mizelpg, max_iter=99999,
            method="L-BFGS",memory=100,
            line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',
            abs_tol=tol*1e8,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
        optimfit$value = optimfit$f
        }
        if(stochastic) optimfit <- sgd(init, itertol=1,startnrows=startnrows,maxiter=300)
        standata$nopriors <- as.integer(1)
        smf<-stan_reinitsf(sm,standata,fast=TRUE)
        parallel::clusterExport(cl = cl, varlist = c('smf','standata'),envir = environment())
        init = optimfit$par #+ rnorm(length(optimfit$par),0,abs(init/8)+1e-3)#rstan::constrain_pars(object = smf, optimfit$par)
      }

      if(!stochastic) {
        # optimfit <- ucminf(init,fn = neglpgf,gr = grffromlp,control=list(grtol=1e-99,xtol=tol,maxeval=10000),hessian=2)
        
        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          # method = 'NAG', nest_q = .01, #nest_convex_approx=TRUE,
          method="L-BFGS",memory=100,
          # method='SR1',
          line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',
          abs_tol=tol,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
        optimfit$value = optimfit$f
        init = optimfit$par
        bestfit <- -optimfit$value
        optimfit$value <- -optimfit$value
       
      }
      
      if(stochastic || is.infinite(bestfit)){
        if(is.infinite(bestfit)) {
          message('Switching to stochastic optimizer -- failed initialisation with bfgs')
        }
        optimfit <- sgd(init)
        bestfit <-optimfit$value
      }
      
      
      est2=optimfit$par #unconstrain_pars(smf, est1)
      # parallel::stopCluster(cl)
    }
    
    if(!estonly){
      lpg<-function(parm) {
        out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          gradout <<- attributes(out)$gradient
        }
        return(out)
      }
      
      grmat<-function(pars,step=1e-5,lpdifmin=1e-8, lpdifmax=1e-3){
        
        # hessout <- parallel::parSapply(cl=cl, X = 1:length(pars), function(i) {
          hessout <- sapply(X = 1:length(pars), function(i) {
          smf <- stan_reinitsf(sm,standata,fast=TRUE)
          # for(i in 1:length(pars)){
          stepsize <- step #*10
          colout <- NA
          while(any(is.na(colout)) && stepsize > 1e-14 && stepsize < 1e8){
          lpdifok<-FALSE
          lpdifcount <- 0
          lpdifdirection <- 0
          lpdifmultiplier <- 1
          while(!lpdifok & lpdifcount < 15){
            # message(paste(i,'  col=',colout,'  lpdifmultiplier=',lpdifmultiplier, '  stepsize=',stepsize))
            lpdifok <- TRUE
            lpdifcount <- lpdifcount + 1
            uppars<-pars
            downpars<-pars
            uppars[i]<-pars[i]+stepsize
            downpars[i]<-pars[i]-stepsize
            uplp= smf$log_prob(upars=uppars,adjust_transform=TRUE,gradient=TRUE) #lpg(uppars)
            upgrad= attributes(uplp)$gradient #gradout
            downlp = smf$log_prob(upars=downpars,adjust_transform=TRUE,gradient=TRUE) #lpg(downpars)
            downgrad = attributes(downlp)$gradient #gradout
            
            
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
              # browser()
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
            # hessout[i,]<- (upgrad-downgrad) /stepsize/2
            colout<- (upgrad-downgrad) /stepsize/2
          }
          return(colout)
        })
        return(t(hessout))
      }
      
      
      hess1s<-function(pars,direction=1,step=1e-5,lpdifmin=1e-6, lpdifmax=1e-1){
        hessout<-matrix(NA,nrow=length(pars),ncol=length(pars))
        bestlp=lpg(pars)
        basegrad=gradout
        for(i in 1:length(pars)){
          stepsize <- step 
          lpdifok<-FALSE
          lpdifcount <- 0
          lpdifdirection <- 0
          lpdifmultiplier <- 1
          while(!lpdifok & lpdifcount < 15){
            lpdifok <- TRUE
            lpdifcount <- lpdifcount + 1
            uppars<-pars
            uppars[i]<-pars[i]+stepsize*direction
            uplp=lpg(uppars)
            upgrad=gradout
            if(abs(uplp-bestlp) > lpdifmax) {
              # message(paste0('decreasing step for ', i))
              lpdifok <- FALSE
              if(lpdifdirection== 1) {
                lpdifmultiplier = lpdifmultiplier * .5
              }
              stepsize = stepsize * 1e-2 * lpdifmultiplier
              lpdifdirection <- -1
            }
            if(abs(uplp-bestlp) < lpdifmin) {
              # message(paste0('increasing step for ', i))
              lpdifok <- FALSE
              if(lpdifdirection== -1) {
                lpdifmultiplier = lpdifmultiplier * .5
              }
              stepsize = stepsize * 100 * lpdifmultiplier
              lpdifdirection <- 1
            }
            hessout[i,]<- (upgrad-basegrad) /stepsize*direction
          }
          
          # }
        }
        return(t(hessout))
      }
      
      
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
        hess=grmat(pars=est2,step=1e-4)
        if(any(is.na(hess))) stop(paste0('Hessian could not be computed for pars ', paste0(which(apply(hess,1,function(x) any(is.na(x))))), ' -- consider reparameterising.',collapse=''))
        hess = (hess/2) + t(hess/2)
        # neghesschol = try(chol(-hess),silent=TRUE)
        
        mchol=try(t(chol(solve(-hess))),silent=TRUE)
        if(class(mchol)=='try-error') {
          message('Hessian not positive-definite so approximating, treat SE\'s with caution, consider respecification / priors.')
          npd <- TRUE
        } else npd <- FALSE
        # if(class(mchol)=='try-error') {
        mcov=MASS::ginv(-hess) #-optimfit$hessian)
        mcov=as.matrix(Matrix::nearPD(mcov)$mat)
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
      samples <-c()
      resamples <- c()
      prop_dens <-c()
      target_dens<-c()
      sample_prob<-c()
      counter <- 0
      ess <- 0
      qdiag<-0
      
      if(isloops == 0) {
        nresamples = finishsamples
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + t(chol(mcovl[[1]])) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        message('Importance sampling not done -- interval estimates via hessian based sampling only')
      }
      
      if(isloops > 0){
        message('Adaptive importance sampling, loop:')
        j <- 0
        while(j < isloops){
          j<- j+1
          message(paste0('  ', j, ' / ', isloops, '...'))
          if(j==1){
            # if(!npd) 
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            # samples <- mvtnorm::rmvnorm(isloopsize, delta[[j]], sigma = mcovl[[j]])
            # 
            # gensigstates <- function(t0means, t0chol,steps){
            #   out <- NA
            #   nlatent=nrow(t0chol)
            #   t0states <- matrix(t0means,byrow=TRUE,nlatent*2,nlatent)
            #   t0base <- matrix(t0means,byrow=TRUE,nlatent,nlatent)
            #   for(sqrtukfadjust in steps){
            #     sigpoints <- t(chol(as.matrix(Matrix::bdiag((t0chol%*%t(t0chol))))))*sqrtukfadjust
            #     t0states[1:(nlatent),] =  t0base + t(sigpoints)
            #     t0states[(1+nlatent):(nlatent*2),] = t0base - t(sigpoints)
            #     if(is.na(out[1])) out <- t0states else out <- rbind(out,t0states)
            #   }
            #   return(out)
            # }
            # samples=gensigstates(delta[[j]],mcovl[[j]],5000*(.1^(exp(seq(0,1,length.out=ceiling(isloopsize/npars/2))))))
            # samples=samples[sample(1:nrow(samples),isloopsize),]
            # 
            # if(npd){
            #   samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
            #   prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize*10), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
            #   samples <- samples[prop_dens > quantile(prop_dens,probs = .9),]
            #   prop_dens <- prop_dens[prop_dens > quantile(prop_dens,probs = .9)]
            # }
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
          
          parallel::clusterExport(cl, c('samples'),environment())
          
          target_dens[[j]] <- unlist(parallel::parLapply(cl, parallel::clusterSplit(cl,1:isloopsize), function(x){
            # eval(parse(text=paste0('library(rstan)')))
            
            smf <- stan_reinitsf(sm,standata, fast=TRUE)
            
            lp<-function(parm) {
              out<-try(smf$log_prob(upars=parm, adjust_transform = TRUE, gradient=FALSE),silent = TRUE)
              if(class(out)=='try-error') {
                out=-Inf
              }
              return(out)
            }
            out <- apply(tail(samples,isloopsize)[x,],1,lp)
            
            try(dyn.unload(file.path(tempdir(), paste0(smf@stanmodel@dso@dso_filename, .Platform$dynlib.ext))),silent = TRUE)
            return(out)
            
          }))
          
          if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
          if(any(target_dens[[j]] > bestfit && (j < isloops && !try2))){
            oldfit <- bestfit
            try2 <- TRUE
            bestfit<-max(target_dens[[j]],na.rm=TRUE)
            betterfit<-TRUE
            init = samples[which(unlist(target_dens) == bestfit),]
            message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
            break
          }
          nresamples = ifelse(j==isloops,finishsamples,5000)
          
          
          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)
          
          sample_prob <- c(sample_prob,exp((weighted_dens - log_sum_exp(weighted_dens)))) #sum to 1 for each iteration, normalise later
          # if(j==isloops) isloopsize = length(sample_prob) #on last loop use all samples for resampling
          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
          # points(target_dens2[sample_prob> (1/isloopsize * 10)], prop_dens[sample_prob> (1/isloopsize * 10)],col='red')
          resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(j == isloops,FALSE,TRUE),
            prob = sample_prob / sum(sample_prob))
          # resample_i <- sample(tail(1:nrow(samples),isloopsize), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE), 
          #   prob = tail(sample_prob,isloopsize) / sum(tail(sample_prob,isloopsize) ))
          if(j < isloops){
            message(paste0(length(unique(resample_i)), ' unique samples drawn, from ', nresamples,' resamples of ', nrow(samples),' actual, probability sd = ', sd(sample_prob)))
            if(length(unique(resample_i)) < 100) {
              message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
              # return(est)
            }
          }
          resamples <- samples[resample_i, , drop = FALSE]
          # points(target_dens2[resample_i],prop_dens[resample_i],col='blue')
          # resamples=mcmc(resamples)
          
          ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
          qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
            (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) ) 
          })))
          
        }
      }
    }
  }#end while no better fit
  if(!estonly){
    if(isloops==0) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]
    
    # parallel::stopCluster(cl)
    transformedpars=stan_constrainsamples(sm = sm,standata = standata,samples=resamples,cores=cores)
    
    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if(class(sds)=='try-error') sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smfull, upars=lest)),
      unlist(constrain_pars(smfull, upars= est2)),
      unlist(constrain_pars(smfull, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    
    # parallel::stopCluster(cl)
    #to return proper stanfit object
    
    stanfit=list(optimfit=optimfit,stanfit=smfull, rawest=est2, rawposterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  suppressWarnings(do.call(par,parbase)) #reset par in case plots done
  return(stanfit)
}

