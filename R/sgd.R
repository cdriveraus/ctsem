logit = function(x) log(x)-log((1-x))

sgd <- function(init,fitfunc,ndatapoints=NA,plot=FALSE,
  stepbase=1e-3,gmeminit=ifelse(is.na(startnrows),.8,.8),gmemmax=.9,maxparchange = .50,
  startnrows=NA,roughnessmemory=.95,groughnesstarget=.4,lproughnesstarget=.3,
  gsmoothroughnesstarget=.05,
  warmuplength=20,
  minparchange=1e-800,maxiter=50000,nconvergeiter=30, itertol=1e-3, deltatol=1e-5){
  
  errsum = function(x) sqrt(sum(abs(x)))
  pars=init
  delta=rep(0,length(pars))
  bestpars = newpars=maxpars=minpars=changepars=pars

  step=rep(stepbase,length(init))
  lpg= fitfunc(init) # try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
  if('try-error' %in% class(lpg)) {
    i = 0
    message('Problems initialising, trying random values...')
    while(i < 50 && 'try-error' %in% class(lpg)){
      if(i %%5 == 0) init = rep(0,length(init))
      init=init+rnorm(length(init),0,abs(init)+ .1)
      lpg= fitfunc(init) 
      i = i + 1
    }
  }
  g=sign(attributes(lpg)$gradient)*(abs(lpg))^(1/8)
  gsmooth=oldgsmooth=oldg=gmid=g
  ghatmix=oldghatmix=g
  
  dghatsmooth=dghat=dghatmix=rep(1,length(g))
  dghatsmoothpar = rep(1,length(pars))
  dghatweight = rep(.2,length(pars)) #ideally one step length...
  dghatdownerr = dghatmixerr = dghatuperr = rep(0,length(g))
  
  lprdif = lpdif = 0

  groughness = rep(groughnesstarget,length(g))
  gsmoothroughness = rep(gsmoothroughnesstarget,length(g))
  lproughness=oldlproughnesstarget=lproughnesstarget
  gmemory <- gmeminit
  oldgmemory  <- gmemory
  ghatsmoothpar <- rep(.8,length(g))
  oldlpdif <- 0
  lpdif <- 0
  maxlp <- -Inf
  i=0
  lp<-c()
  oldlp <- -Inf
  converged <- FALSE
  while(!converged && i < maxiter){
    i = i + 1
    accepted <- FALSE
    lproughnesstarget2 = lproughnesstarget 
    notacceptedcount <- 0
    
    while(!accepted){
      notacceptedcount <- notacceptedcount+1
      if(notacceptedcount > 50) {
        stop('Cannot optimize! Problematic model, or bug?')
        print(lpg)
      }
      
      gdelta =  (ghatmix +dghatmix*step*dghatweight/2)
      if(i > 1){
        delta =   step  *sign(gdelta)*(abs(gdelta))  * exp((rnorm(length(g),0,.02)))
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        newpars = pars + delta
        newpars = newpars  + delta/2 - deltaold/2 #+ delta - deltaold #
      }
      
      lpg= fitfunc(newpars)
      
      if(lpg > -1e99 &&       #regular check
          class(lpg) !='try-error' && 
          !is.nan(lpg[1]) && 
          all(!is.nan(attributes(lpg)$gradient)) 
        && (i < warmuplength || (lp[i-1]- lpg[1]) < sd(tail(lp,20))*8+1e-6)
      ){
        accepted <- TRUE
      } 
      else {
        ghatmix=g*.5
        # gsmooth= gsmooth*gmemory2 + (1-gmemory2) * g #increase influence of last gradient at inflections
        step <- step * .5
        # pars=bestpars
      }
      #warmup check
      if(is.na(startnrows) && i < warmuplength && i > 1 && lpg[1] < lp[1]-5) {
        # 
        accepted <- FALSE
        step = step * .1
        pars=bestpars
      }
      if(plot && !accepted) {
        print(paste0('iter ', i,' not accepted!'))
        # 
      }
    } #end acceptance loop
    
    #once accepted 
    lp[i]=lpg[1]
    pars=newpars
    
    
    
    deltaold=delta
    oldg=g
    g=attributes(lpg)$gradient
    g=sign(g)*(abs(g))^(1/4)#sqrt
    gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
    roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)

   
    oldgmid=gmid
    gmid = (oldg+g)/2
    dg=(gmid-oldgmid)/step
    
    #compute pred errors before new predictions

    # dghatdown = (dghatsmooth * (dghatsmoothpar-.01) + dghat*(1-(dghatsmoothpar-.01)))
    # dghatup = (dghatsmooth * (dghatsmoothpar+.01) + dghat*(1-(dghatsmoothpar+.01)))
    dghatdownerr <- dghatdownerr * gmemory2 + (1-gmemory2)*(abs(dg - (dghatsmooth * (dghatsmoothpar-.01) + dghat*(1-(dghatsmoothpar-.01)))))
    dghatuperr <- dghatuperr * gmemory2 + (1-gmemory2)* (abs(dg - (dghatsmooth * (dghatsmoothpar+.01) + dghat*(1-(dghatsmoothpar+.01)))))
    dghatmixerr <- dghatmixerr * gmemory2 + (1-gmemory2)* (abs(dg-dghatmix))
    down = dghatmixerr > dghatdownerr
    up= dghatmixerr > dghatuperr
    # dghatsmoothpar[down] = dghatsmoothpar[down]-.01
    # dghatsmoothpar[up] = dghatsmoothpar[up]+.01
    
    # if(i==30) browser()
    ghatmixuperr <- abs(g- ( (inv_logit(logit(ghatsmoothpar)+.1))*gsmooth + (1-(inv_logit(logit(ghatsmoothpar)+.1)))*oldg ))
    ghatmixdownerr <- abs(g- ( (inv_logit(logit(ghatsmoothpar)-.1))*gsmooth + (1-(inv_logit(logit(ghatsmoothpar)-.1)))*oldg ))
    ghatmixerr <- abs(g-ghatmix)
    ghatsmoothpar[ghatmixerr > ghatmixdownerr]=(inv_logit(logit(ghatsmoothpar[ghatmixerr > ghatmixdownerr])-.1))
    ghatsmoothpar[ghatmixerr > ghatmixuperr]=(inv_logit(logit(ghatsmoothpar[ghatmixerr > ghatmixuperr])+.1))
    # print(mean(as.numeric(ghatmixerr > ghatmixuperr)))
    # print(mean(as.numeric(ghatmixerr > ghatmixdownerr)))
    
    #  if(i %% 10 ==0){
    # if(ghatsmoothpar > .95 && gmemory < .98) gmemory = inv_logit(logit(gmemory)+.2)
    # if(ghatsmoothpar < .2&& gmemory > .2) gmemory = inv_logit(logit(gmemory)-.2)
    #  }
    
# test changes in dghatweight

    up = which(abs( g-(ghatmix + step * dghatmix * dghatweight)) > abs(g-(ghatmix + step * dghatmix * (dghatweight+.01))) )
    down = which(abs( g-(ghatmix + step * dghatmix * dghatweight)) > abs(g-(ghatmix + step * dghatmix * (dghatweight-.01))) )
    # dghatweight[up]=dghatweight[up]+.01
    # dghatweight[down]=dghatweight[down]-.01
    
    
    
    #predictions
    oldgsmooth = gsmooth
    gsmooth= gsmooth*gmemory2 + (1-gmemory2) * g 
    dghat = dg
    dghatsmooth = gmemory2*dghatsmooth +(1-gmemory2)*dghat
    dghatmix = dghatsmooth * dghatsmoothpar + dghat*(1-dghatsmoothpar)
    
    
    
    
    
    groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gmid)!=sign(oldgmid))
    gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
    if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] > (lp[i]))
    
    lproughnessmod=  ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    gsmoothroughnessmod =  (( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    groughnessmod = ( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    
    step = (step + .5*(
      step* .8*lproughnessmod
      # + step* .1*gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      + step* .3*groughnessmod# * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    ))

    # step[gsmoothroughness < gsmoothroughnesstarget] <- step[gsmoothroughness < gsmoothroughnesstarget] + 
      # step[gsmoothroughness < gsmoothroughnesstarget] * .1*gsmoothroughnessmod[gsmoothroughness < gsmoothroughnesstarget]
    signdif= sign(gmid)!=sign(gdelta)
    if(lp[i] >= max(lp)) {
      step = step * 1.1 #sqrt(2-gmemory) #exp((1-gmemory)/8)
      if(i > warmuplength/2) {
        ##max/min par update extra
        # gsmooth[pars>maxpars | pars < minpars] <- gsmooth[pars>maxpars | pars < minpars]  * 1.1#*delta[pars>maxpars | pars < minpars] /step[pars>maxpars | pars < minpars]
        # step[pars>maxpars | pars < minpars] <- step[pars>maxpars | pars < minpars] * 1.2  #+ pars[pars>maxpars | pars < minpars]
        changepars=pars
        changepars[!(pars>maxpars | pars < minpars)] <- NA
        # lproughness = lproughness * .9
      }
      bestpars <- pars <- newpars
      bestg <- g
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
    
    # gmemory <- gmemory * gsmoothroughnessmod
    if(i > 30 && i %% 20 == 0) {
      oldlpdif <- lpdif# sum(diff(head(tail(lp,10),20)))
      lpdif <- sum(diff(tail(lp,10)))
      if(oldlpdif > lpdif) gmemory <- oldgmemory
      proposal = gmemory*2-oldgmemory
      oldgmemory <- gmemory
      gmemory <- min(gmemmax, max(0, proposal + runif(1,-.05,.1)))
    }
    
    if(i > 30 && i %% 30 == 0) {
      oldlprdif <- lprdif
      lprdif <- sum(diff(tail(lp,30)))
      if(oldlprdif > lprdif) lproughnesstarget <- oldlproughnesstarget
      lprproposal = lproughnesstarget*2-oldlproughnesstarget
      oldlproughnesstarget <- lproughnesstarget
      lproughnesstarget <- min(.5, max(.05, lprproposal + runif(1,-.05,.05)))
    }
    
    step[step > maxparchange] <- maxparchange
    step[step < minparchange] <- minparchange
    
    if(i > warmuplength && lp[i] < lp[i-1]) { #if worsening, update gradient faster
      step[signdif]=step[signdif]*.5
      # step=step*.5
      gsmooth[signdif]= gsmooth[signdif]*gmemory2 + (1-gmemory2) * g[signdif] #increase influence of gradient at inflections
    }
    
    oldghatmix=ghatmix
    ghatmix = ghatsmoothpar*gsmooth + (1-ghatsmoothpar)*g
    # message('ghatsmoothpar = ',ghatsmoothpar)
    
    
    if(plot){
      parbase=par(no.readonly=TRUE)
      on.exit(do.call(par,parbase),add=TRUE)
      par(mfrow=c(3,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
      plot(bestpars)
      points(changepars,pch=17,col='red')
      
      plot(log(step))
      
      plot(groughness,col='red',ylim=c(0,1))
      abline(h=mean(gsmoothroughness),col='blue',lty=2)
      abline(h=(gsmoothroughnesstarget),col='blue',lty=1,lwd=2)
      points(gsmoothroughness,ylim=c(0,1),col='blue')
      abline(h=mean(groughness),col='red',lty=2)
      abline(h=(groughnesstarget),col='red',lty=1)
      
      abline(h=lproughnesstarget,lty=1,col='green')
      abline(h=lproughness, col='green',lty=2)
      
      plot(tail(log(-(lp-max(lp)-1)),500),type='l')
      gsmoothsqrt=sign(gsmooth) * sqrt(abs(gsmooth))
      # plot(gsmoothsqrt,ylim=c(-max(abs(gsmoothsqrt)),max(abs(gsmoothsqrt))))
      plot(dghatmix*step,ylim=c(-max(abs(dghatmix*step)),max(abs(dghatmix*step))))
      plot(ghatsmoothpar,ylim=c(0,1))
      plot(dghatsmoothpar)
      abline(h=mean(dghatsmoothpar))
      plot(dghatweight)
      abline(h=mean(dghatweight))
      
      # matplot(cbind(signdifmod,gsmoothroughnessmod),col=c('black','blue'),pch=1,ylim=c(-1,1))
      # points(groughnessmod,col='red')
      # abline(h=lproughnessmod,col='green')
      
      message(paste0('Iter = ',i, '   Best LP = ', max(lp),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
    }
    
    #check convergence
    if(i > 30){
      if( (i - bestiter) > nconvergeiter*5) converged <- TRUE #time since best
      if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
      if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
    }
  }
  return(list(itervalues = lp, value = max(lp),par=bestpars))#,gstore=gstore,pstore=pstore) )
}
