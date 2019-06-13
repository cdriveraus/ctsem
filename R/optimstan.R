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
#' @param stochastic Logical. Use stochastic gradient descent instead of ucminf (bfgs with trust region) optimizer.
#' Generally more robust, a little slower with few parameters, faster with many.
#' @param plotsgd Logical. If TRUE, plot iteration details when using stochastic optimizer.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{DEoptim.control}.
#' @param nopriors logical.f If TRUE, any priors are disabled -- sometimes desirable for optimization. 
#' @param startnsubjects Either NA to ignore, or an integer specifying the number of subjects to start the
#' stochastic optimization with. This will be increased to the full amount by convergence time.
#' @param tol absolute object tolerance
#' @param cores Number of cpu cores to use.
#' @param isloops Number of iterations of adaptive importance sampling to perform after optimization.
#' @param isloopsize Number of samples per iteration of importance sampling.
#' @param finishsamples Number of samples to use for final results of importance sampling.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent 
#' importance sampling, at risk of truncating tails.
#'
#' @return ctStanFit object
#' @importFrom ucminf ucminf
#' @importFrom Matrix bdiag
#' @importFrom utils head tail
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
  deoptim=FALSE, estonly=FALSE,tol=1e-12,
  decontrol=list(),
  stochastic = 'auto',
  startnsubjects=NA,
  plotsgd=FALSE,
  isloops=0, isloopsize=1000, finishsamples=500, tdf=50,
  verbose=0,nopriors=FALSE,cores=1){
  
  standata$verbose=as.integer(verbose)
  standata$nopriors=as.integer(nopriors)
  
  if(is.null(decontrol$steptol)) decontrol$steptol=5 
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  
  message('Optimizing...')
  
  betterfit<-TRUE
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    # if(nopriors){
    #   standata$nopriors <- 0
    #   suppressWarnings(suppressOutput(optimfit <- optimizing(sm,standata, hessian=FALSE, iter=400, init=0,as_vector=FALSE,
    #     tol_obj=1e-8, tol_rel_obj=0,init_alpha=.000001, tol_grad=0,tol_rel_grad=1e7,tol_param=1e-5,history_size=100),verbose=verbose))
    #   init=optimfit$par
    #   standata$nopriors <- 1
    # }
    
    # suppressMessages(suppressWarnings(suppressOutput(smf<-sampling(sm,iter=0,chains=0,init=0,data=standata,check_data=FALSE, 
    # control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
    smf <- stan_reinitsf(sm,standata)
    npars=get_num_upars(smf)
    if(all(init %in% 'random')) init <- rnorm(npars, 0, initsd)
    if(all(init == 0)) init <- rep(0,npars)
    
    if(is.na(sampleinit[1])){
      
      if(deoptim){ #init with DE
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
          
          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]
          
          lp2 = function(parm) {
            out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
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
        } else stop(paste0('use install.packages(\"DEoptim\") to use deoptim')) #end require deoptim
      }
      
      
      # suppressWarnings(suppressOutput(optimfit <- optimizing(sm,standata, hessian=FALSE, iter=1e6, init=init,as_vector=FALSE,draws=0,constrained=FALSE,
      #   tol_obj=tol, tol_rel_obj=0,init_alpha=.001, tol_grad=0,tol_rel_grad=0,tol_param=0,history_size=50,verbose=verbose),verbose=verbose))
      
      gradout <- c()
      bestlp <- -Inf
      
      lp<-function(parm) {
        # print((parm)
        out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=-Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          if(out[1] > bestlp) {
            bestlp <<- out[1]
            gradout <<- attributes(out)$gradient
          }
        }
        return(-out[1])
      }
      
      grffromlp<-function(parm) {
        return(-gradout)
      }
      
      parbase=par()
      
      sgd <- function(init,stepbase=1e-4,nsubjects=1,gmeminit=ifelse(is.na(startnsubjects),.7,.9),gmemmax=.9,maxparchange = .1,
        startnsubjects=NA,
        minparchange=1e-16,maxiter=50000,perturbpercent=0,nconvergeiter=20, itertol=1e-3, deltatol=1e-5){
        
        pars=init
        bestpars = pars
        step=rep(stepbase,length(init))
        g=try(grad_log_prob(smf, upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
        gsmooth=g
        gpersist=g
        signg=sign(g)
        oldsigng = g
        oldg=g
        gmemory <- gmeminit
        oldgmemory <- gmemory
        oldlpdif <- 0
        lpdif <- 0
        maxlp <- -Inf
        i=0
        lp<-c()
        oldlp <- -Inf
        converged <- FALSE
        nsubjects <- ifelse(is.na(startnsubjects),standata$nsubjects, min(standata$nsubjects, startnsubjects))
        
        while(!converged && i < maxiter){
          i = i + 1
          accepted <- FALSE
          while(!accepted){
            whichpars = sample(1:length(pars), ceiling(length(pars)*perturbpercent))
            newpars = bestpars #* .5 + bestpars * .5
            parupdate =  exp(rnorm(length(step),0,step)) * step  * gsmooth #(sign(g) +  sign(g)*(abs(gsmooth) / mean(abs(gsmooth)))) #
            if(i %% 28 == 0) step[whichpars] = step[whichpars] * 10
            parupdate[abs(parupdate) > maxparchange] <- maxparchange*sign(parupdate[abs(parupdate) > maxparchange])
            newpars = newpars + parupdate
            
            #sub sampling
            if(!is.na(startnsubjects) || (nsubjects==standata$nsubjects && i > 2)){
              subjects <- sample(unique(standata$subject),nsubjects,replace = FALSE)
              standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,standata$ndatapoints) #
              smf<-stan_reinitsf(sm,standata)
            }
            
            lpg = try(log_prob(smf, upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
            if(class(lpg) !='try-error' && !is.nan(lpg[1]) && all(!is.nan(attributes(lpg)$gradient))) accepted <- TRUE else step <- step * .1
            if(is.na(startnsubjects) && i < 20 && i > 1 && lpg[1] < lp[1]) {
              accepted <- FALSE
              if(plotsgd) print('not accepted!')
              step = step * .5
              gsmooth=gsmooth*.5
            }
          }
          lp[i]=lpg[1]
          pars <- newpars
          oldsigng=signg
          oldg=g
          g=attributes(lpg)$gradient
          oldgsmooth = gsmooth
          gmemory2 = gmemory * min(i/20,1)
          gsmooth= gsmooth*gmemory2 + (1-gmemory2)*g
          signg=sign(g) #or gsmooth?
          stdgdif = (g-oldg) * step
          # print(stdgdif)
          # step=exp(mean(log(step))+(.99*(log(step)-mean(log(step)))))
          # step[oldsigng == signg] = step[which(oldsigng == signg)] * sqrt(2-gmemory) #exp((1-gmemory)/2)
          # step[oldsigng != signg] = step[which(oldsigng != signg)] / sqrt(2-gmemory) #ifelse(nsubjects == standata$nsubjects, (2-gmemory),1.1) #1.2 #exp((1-gmemory)/2)
          
          step[oldsigng == signg] = step[oldsigng == signg] * 1.1 #(.5+inv_logit(abs(stdgdif[oldsigng == signg])))
          step[oldsigng != signg]  = step[oldsigng != signg] * (1.5-inv_logit(abs(stdgdif[oldsigng != signg])))^2
          
          if(lp[i] >= max(lp)) {
            step = step * sqrt(2-gmemory) #exp((1-gmemory)/8)
            bestpars <- pars
          } 
          if(i > 1 && lp[i] < lp[i-1]) {
            # signg <- oldsigng
            # gsmooth = oldgsmooth
            # pars <- bestpars
            # if(nsubjects == standata$nsubjects) {
              # gsmooth=gsmooth*.9
              step = step/ (2-gmemory)#step  / max( 1.5, (-10*(lp[i] - lp[i-1]) / sd(head(tail(lp,20),10)))) #exp((1-gmemory)/4)
            # } else {
            #   step = step / 1.06
            # }
          }
          # if(i %%10 ==0) gmemory = min(gmemory+.1,gmemmax)# * exp(mean(sign(diff(tail(lp,20)))))
          # if(i %%20 ==0) gmemory =  max(gmeminit, min(gmemmax, 1.6*(1-(log(sd(tail(lp,20)) ) -log(itertol)) / (log(sd(head(lp,20)))-log(itertol)))* (1-gmeminit) + gmeminit))
          
          if(i > 30 && i %% 10 == 0) {
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
            par(mfrow=c(3,1))
            plot(pars)
            plot(log(step))
            plot(tail(log(-(lp-max(lp)-1)),500),type='l')
            message(paste0('Iter = ',i, '   Best LP = ', max(lp),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
          }
          
          #check convergence
          if(i > 30){
            if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
            # print(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)))
            if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
            if(nsubjects < standata$nsubjects && (length(lp) - match(min(lp),lp)) > nconvergeiter) converged <- TRUE
          }
          if(converged & nsubjects != standata$nsubjects){
            converged <- FALSE
            nsubjects <- min(standata$nsubjects, nsubjects * 2)
            if(nsubjects > standata$nsubjects/2) nsubjects <- standata$nsubjects
            message(paste0('nsubjects now ',nsubjects, ' out of ',standata$nsubjects),' total')
            i=0
            lp=c()
          }
        }
        return(list(itervalues = lp, value = max(lp),par=bestpars) )
      }
      # browser()
      # if(!deoptim & standata$nopriors == 0 ) init='random'
      
      if(stochastic=='auto' && npars > 50){
        message('> 20 parameters and stochastic="auto" so stochastic gradient descent used')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE
      
      if(!deoptim & standata$nopriors == 1 ){ #init using priors
        standata$nopriors <- as.integer(0)
        smf <- stan_reinitsf(sm,standata)
        if(!stochastic) optimfit <- ucminf(init,fn = lp,gr = grffromlp,control=list(xtol=tol*1e6,maxeval=10000))
        if(stochastic) optimfit <- sgd(init, itertol=1,startnsubjects=startnsubjects)
        standata$nopriors <- as.integer(1)
        smf <- stan_reinitsf(sm,standata)
        init = optimfit$par #rstan::constrain_pars(object = smf, optimfit$par)
      }
      
      if(!stochastic) {
        optimfit <- ucminf(init,fn = lp,gr = grffromlp,control=list(grtol=1e-99,xtol=tol,maxeval=10000),hessian=2)
        init = optimfit$par
        optimfit$value <- -optimfit$value
        ucminfcov <- optimfit$invhessian
      }

      if(stochastic){
      optimfit <- sgd(init, nsubjects=standata$nsubjects,startnsubjects=startnsubjects, perturbpercent = .0) 
      }
   
      est1=constrain_pars(smf,optimfit$par)
      bestfit <-optimfit$value
      
      est2=optimfit$par #unconstrain_pars(smf, est1)
    }
    
    if(!estonly){
      lpg<-function(parm) {
        out<-try(log_prob(smf, upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          gradout <<- attributes(out)$gradient
        }
        return(out[1])
      }
      
      grf<-function(parm,...) {
        out=try(grad_log_prob(smf, upars=parm, adjust_transform = TRUE))
        if(class(out)=='try-error') {
          out=rep(NA,length(parm))
        }
        return(out)
      }

      grmat<-function(pars,step=1e-5,lpdifmin=1e-8, lpdifmax=1e-3){
        hessout<-matrix(NA,nrow=length(pars),ncol=length(pars))
        for(i in 1:length(pars)){
          stepsize <- step #*10
          # while((any(is.na(hessout[i,])) || hessout[i,i] >=0)  && stepsize > 1e-12){
            # stepsize <- step * .1
            lpdifok<-FALSE
            lpdifcount <- 0
            lpdifdirection <- 0
            lpdifmultiplier <- 1
            while(!lpdifok & lpdifcount < 15){
              lpdifok <- TRUE
              lpdifcount <- lpdifcount + 1
              uppars<-pars
              downpars<-pars
              uppars[i]<-pars[i]+stepsize
              downpars[i]<-pars[i]-stepsize
              uplp=lpg(uppars)
              upgrad=gradout
              downlp = lpg(downpars)
              downgrad = gradout
              # print(abs(uplp-downlp))
              # print(upgrad)
              # print(downgrad)
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
                # message(paste0('increasing step for ', i))
                lpdifok <- FALSE
                if(lpdifdirection== -1) {
                  lpdifmultiplier = lpdifmultiplier * .5
                }
                stepsize = stepsize * 100 * lpdifmultiplier
                lpdifdirection <- 1
              }
              hessout[i,]<- (upgrad-downgrad) /stepsize/2
            }
            
          # }
        }
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
        # browser()
        # hessup=hess1s(pars = est2,direction = 1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hessdown=hess1s(pars = est2,direction = -1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hess=(hessup+hessdown)/2
        hess=grmat(pars=est2,step=1e-4)
        if(any(is.na(hess))) stop(paste0('Hessian could not be computed for pars ', paste0(which(apply(hess,1,function(x) any(is.na(x))))), ' -- consider reparameterising.',collapse=''))
        hess = (hess/2) + t(hess/2)
        # neghesschol = try(chol(-hess),silent=TRUE)
        
        mchol=try(t(chol(solve(-hess))),silent=TRUE)
        if(class(mchol)=='try-error') {
          message('Hessian not positive-definite -- check importance sampling convergence with isdiag')
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
      
      ctdmvnorm <- function(samples, mu, sigma){
        sigi <- solve(sigma)
        sigdet <- det(sigma)
        d <- 
          log(1/(sqrt((2*pi)^length(mu)*sigdet))) + apply(samples,1, function(x) {
            (-1/2*( c(x-mu) %*% sigi %*% c(x-mu)))
          })
        return(d)
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
      
      cl <- parallel::makeCluster(cores, type = "PSOCK")
      parallel::clusterExport(cl, c('sm','standata'),environment())
      
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
            # browser()
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
            eval(parse(text=paste0('library(rstan)')))
            
            smf <- stan_reinitsf(sm,standata)
            
            lp<-function(parm) {
              out<-try(log_prob(smf, upars=parm, adjust_transform = TRUE, gradient=FALSE),silent = TRUE)
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
            # init = rstan::constrain_pars(object = smf, samples[which(unlist(target_dens) == bestfit),])
            init = samples[which(unlist(target_dens) == bestfit),]
            message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
            break
          }
          nresamples = ifelse(j==isloops,finishsamples,5000)
          
          
          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # browser()
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
    message('Computing quantities...')
    
    # cl <- parallel::makeCluster(min(cores,chains), type = "PSOCK")
    parallel::clusterExport(cl, c('relistarrays','resamples','sm','standata','optimfit'),environment())
    
    # target_dens <- c(target_dens,
    
    transformedpars <- try(parallel::parLapply(cl, parallel::clusterSplit(cl,1:nresamples), function(x){
      require(ctsem)
      Sys.sleep(.1)
      smf <- stan_reinitsf(sm,standata)
      Sys.sleep(.1)
      # smf<-new(sm@mk_cppmodule(sm),standata,0L,rstan::grab_cxxfun(sm@dso))
      out <- list()
      skeleton=est1
      for(li in 1:length(x)){
        out[[li]] <- try(rstan::constrain_pars(smf, resamples[x[li],]))
        if(any(sapply(out[[li]], function(x) any(c(is.nan(x),is.infinite(x),is.na(x)))))) class(out[[li]]) <- c(class(out[[li]]),'try-error')
      }
      return(out)
    }))
    transformedpars=unlist(transformedpars,recursive = FALSE)
    missingsamps <-sapply(transformedpars, function(x) 'try-error' %in% class(x))
    nasampscount <- sum(missingsamps) 
    
    transformedpars <- transformedpars[!missingsamps]
    nresamples <- nresamples - nasampscount
    if(nasampscount > 0) {
      message(paste0(nasampscount,' NAs generated during final sampling of ', finishsamples, '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
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
    
    
    transformedpars=try(tostanarray(flesh=matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = est1))
    
    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if(class(sds)=='try-error') sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, lest)),
      unlist(constrain_pars(smf, est2)),
      unlist(constrain_pars(smf, uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    
    parallel::stopCluster(cl)
    
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2, rawposterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  suppressWarnings(do.call(par,parbase)) #reset par in case plots done
  return(stanfit)
}

