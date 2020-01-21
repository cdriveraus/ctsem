sgd <- function(init,fitfunc,ndatapoints=NA,plot=FALSE,
  stepbase=1e-1,gmeminit=ifelse(is.na(startnrows),.8,.8),gmemmax=.95,maxparchange = 50,
  startnrows=NA,gsmoothness = 1,roughnessmemory=.95,groughnesstarget=.5,lproughnesstarget=.2,
  gsmoothroughnesstarget=.1,
  warmuplength=20,
  minparchange=1e-40,maxiter=50000,nconvergeiter=30, itertol=1e-3, deltatol=1e-5){
  pars=init
  ngstore <- length(pars)*2
  pstore <- gstore <- matrix(NA,ngstore, length(pars))
  bestpars = pars
  maxpars=pars
  minpars=pars
  changepars=pars
  step=rep(stepbase,length(init))
  # parallelStanSetup(cores,sm,standata)

  g= fitfunc(init) # try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
  if('try-error' %in% class(g)) {
    i = 0
    message('Problems initialising, trying random values...')
    while(i < 50 && 'try-error' %in% class(g)){
      if(i %%5 == 0) init = rep(0,length(init))
      init=init+rnorm(length(init),0,abs(init)+ .1)
      g= fitfunc(init) #try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE)
      i = i + 1
    }
  }
  g=sign(attributes(g)$gradient)#*.1*sqrt(abs(g))
  gsmooth=g
  oldg=g
  groughness = rep(groughnesstarget,length(g))
  gsmoothroughness = rep(gsmoothroughnesstarget,length(g))
  # deltasmoothsq=.01
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
  # nrows <- ifelse(is.na(startnrows),ndatapoints, min(ndatapoints, startnrows))
  while(!converged && i < maxiter){
    i = i + 1
    accepted <- FALSE
    lproughnesstarget2 = lproughnesstarget #ifelse(nrows==ndatapoints,lproughnesstarget,.49)
    notacceptedcount <- 0
    while(!accepted){
      notacceptedcount <- notacceptedcount+1
      if(notacceptedcount > 50) {
        stop('Cannot optimize! Problematic model, or bug?')
        print(lpg)
      }
      newpars = bestpars
      if(i > 1){
      delta =   step  * gsmooth * exp((rnorm(length(g),0,.01)))
      delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
      newpars = newpars + delta
      }
      
      #sub sampling
      # if(!is.na(startnrows) || (nrows!=ndatapoints)){
      #   subjects <- sample(1:ndatapoints,nrows,replace = FALSE)
      #   standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,ndatapoints) #
      #   parallelStanSetup(cores = cores,sm = sm,standata = standata) #could be improved for subsampling
      # }
      
      # lpg = try(smf$log_prob(upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
      lpg= fitfunc(newpars)
      #regular check
      if(lpg > -1e99 && 
          class(lpg) !='try-error' && 
          !is.nan(lpg[1]) && 
          all(!is.nan(attributes(lpg)$gradient)) &&
           (i < warmuplength || lp[i-1]- lpg[1] < sd(tail(lp,5))*4+1e-6)){
        accepted <- TRUE
      } 
      else {
        step <- step * .5
      }
      #warmup check
      if(is.na(startnrows) && i < warmuplength && i > 1 && lpg[1] < lp[1]) {
        accepted <- FALSE
        step = step * .5
        gsmooth=gsmooth*.5
      }
      if(plot && !accepted) print(paste0('iter ', i,' not accepted!'))
    }
    lp[i]=lpg[1]
    pars <- newpars

    gstore[(i-1) %% ngstore +1,] <- g
    pstore[(i-1) %% ngstore +1,] <- pars
    oldg=g
    g=attributes(lpg)$gradient
    g=sign(g)*.1*sqrt(abs(g))
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
    signdifmod[sign(oldg) == sign(g)] =  1 #/ (1.5-inv_logit(abs(stdgdif[oldsigng == signg])))^4 #* (1/ ( ( (roughness*.05+.95)^2) ))
    signdifmod[sign(oldg) != sign(g)]  = -1 #10* ((1.5-inv_logit(abs(stdgdifold[sign(oldg) != sign(g)])))-1) #* ( ( (roughness*.05+.95)^2) )
    signdifmod[is.nan(signdifmod)] <- 0 #oldstep[is.nan(step)] #because of overflow in some cases
    
    # deltasmoothsq = deltasmoothsq * gmemory + (1-gmemory)*delta^2
    lproughnessmod=  ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    # gmemoryupd = min(gmemmax,max(.1,gmemory /  ( (1/(-mean(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ))
    gsmoothroughnessmod =  .1 *(( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    groughnessmod = .5 *( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    # rmsstepmod = sqrt(abs(gsmooth+1e-7))/step -1 #like adagrad but with decaying gradient
    
    step = (step
      + step*signdifmod *.1#* min(sqrt(deltasmoothsq),1)
      + step*lproughnessmod
      + step* 0*gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      + step* groughnessmod# * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    )

    if(lp[i] >= max(lp)) {
      # step = step * sqrt(2-gmemory) #exp((1-gmemory)/8)
      if(i > warmuplength/2) {
        ##max/min par update extra
        gsmooth[pars>maxpars | pars < minpars] <- gsmooth[pars>maxpars | pars < minpars]  * 1.1#*delta[pars>maxpars | pars < minpars] /step[pars>maxpars | pars < minpars]
        step[pars>maxpars | pars < minpars] <- step[pars>maxpars | pars < minpars] * 1.2  #+ pars[pars>maxpars | pars < minpars]
        changepars=pars
        changepars[!(pars>maxpars | pars < minpars)] <- NA
        lproughness = lproughness * .9
      }
      # pars <- newpars
      bestpars <- pars
      bestiter <- i
      
      maxpars[pars>maxpars] <-pars[pars>maxpars]
      minpars[pars<minpars] <-pars[pars<minpars]
      
      
    }
    if(i > 1 && runif(1,0,1) > .95) {
      # #slowly forget old max and mins, allow fast re exploration of space
      rndchange <- runif(length(maxpars),0,1) > .95
      maxpars[rndchange] <- pars[rndchange]
      minpars[rndchange] <- pars[rndchange]
    }

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
    
    if(plot){
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
      
      message(paste0('Iter = ',i, '   Best LP = ', max(lp),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
    }
    
    #check convergence
    if(i > 30){
      if( (i - bestiter) > nconvergeiter*5) converged <- TRUE #time since best
      if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
      if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
    }
  }
  return(list(itervalues = lp, value = max(lp),par=bestpars,gstore=gstore,pstore=pstore) )
}
