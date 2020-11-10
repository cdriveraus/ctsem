logit = function(x) log(x)-log((1-x))

sgd <- function(init,fitfunc,whichignore=c(),whichmcmcpars=NA,mcmcstep=.01,nsubjects=NA,ndatapoints=NA,plot=FALSE,
  stepbase=1e-3,gmeminit=ifelse(is.na(startnrows),.8,.8),gmemmax=.93, maxparchange = .50,
  startnrows=NA,roughnessmemory=.9,groughnesstarget=.4,roughnesschangemulti = 2,
  lproughnesstarget=ifelse(is.na(whichmcmcpars[1]),ifelse(parsets==1,.2,.1),.2),parsets=1,
  # gamiter=50000,
  gsmoothroughnesstarget=.05,
  warmuplength=20,nstore=max(100,length(init)),
  minparchange=1e-800,maxiter=50000,
  nconvergeiter=ifelse(is.na(whichmcmcpars[1]),30,60), 
  itertol=1e-3, deltatol=1e-5, parsdtol=1e-3){
  
  initfull=init #including ignored params start values
  if(length(whichignore)>0) init=init[-whichignore]
  
  errsum = function(x) sqrt(sum(abs(x)))
  
  combinepars <- function(pars,mcmcpars,whichmcmcpars,subject=NA){
    out <- rep(NA,length(c(pars,mcmcpars)))
    out[whichmcmcpars] <- mcmcpars
    out[-whichmcmcpars] <- pars
    if(!is.na(subject)) out <- c(out,subject)
    return(out)
  }
  
  if(plot){
    parbase=par(no.readonly=TRUE)
    on.exit(do.call(par,parbase),add=TRUE)
  }
  
  fitfunc2 <- function(x,whichmcmcpars=NA){
    out <- fitfunc(x)
    if(!is.na(whichmcmcpars[1])) attributes(out)$gradient <- attributes(out)$gradient[-whichmcmcpars]
    return(out)
  }
  
  mcmcfitfunc <- function(newpars,newmcmcpars,mcmcconverged){
    nsubsinsample <- ceiling(ifelse(mcmcconverged,1,1)*nsubjects)
    mcmclpg <- list()
    subjects <- sample(1:nsubjects,nsubsinsample)
    for(subi in subjects){
      mcmcaccepted <- mcmcconverged
      # 
      mcmccount <- 0
      mcmclpgb= fitfunc2(combinepars(newpars,newmcmcpars[,subi],whichmcmcpars,subi),
        c(whichmcmcpars,length(init))) #fit with updated pars and old mcmcpars
      while(!mcmcaccepted && mcmccount < 50){ #sample new mcmcpars until accepted iteration
        mcmccount <- mcmccount + 1
        newmcmcpars[,subi] <- mcmcpars[,subi] + 
          # propchol %*% 
          rnorm(length(mcmcpars[,subi]),0,mcmcstep)
        mcmclpg[[subi]]= fitfunc2(combinepars(newpars,newmcmcpars[,subi],whichmcmcpars,subi),
          c(whichmcmcpars,length(init)))
        if(runif(1,0,1) < exp(mcmclpg[[subi]]-mcmclpgb)){
          mcmcaccepted <- TRUE
          ar<<-.95*ar+.05
        } else ar<<-.95*ar
      }
      if(mcmcconverged) mcmclpg[[subi]] <- mcmclpgb
    }
    # 
    # mcmcdiff <- newmcmcpars[,subjects,drop=FALSE] - mcmcpars[,subjects,drop=FALSE]
    # mcmcdiv <- (1-apply(newmcmcpars[,subjects,drop=FALSE],1,sum) / 
    #   apply(mcmcpars[,subjects,drop=FALSE],1,sum)) * nsubsinsample/nsubjects * .1
    # newmcmcpars[,-subjects] <- newmcmcpars[,-subjects] * (1+mcmcdiv) + 
    #   apply(mcmcdiff,1,mean) *nsubsinsample/nsubjects
    
    mcmclpg <- mcmclpg[!sapply(mcmclpg,is.null)]
    
    mcmclpg[[length(mcmclpg)+1]] = fitfunc2(combinepars(newpars,newmcmcpars[,1],whichmcmcpars,0),
      c(whichmcmcpars,length(init)))
    lpg <- sum(unlist(mcmclpg)) * nsubjects/nsubsinsample
    attributes(lpg)$gradient <- apply(sapply(mcmclpg,function(x) attributes(x)$gradient,simplify='matrix'),1,sum)
    return(list(lpg=lpg,newmcmcpars=newmcmcpars))
  }
  pars=init
  mcmcconverged <- TRUE
  if(!is.na(whichmcmcpars[1])){
    artarget=.23
    ar=artarget
    mcmcconverged<-FALSE
    pars=init[-whichmcmcpars]
    pars=pars[-length(pars)] #remove subject specifier
    newmcmcpars=matrix(rnorm(length(whichmcmcpars)*nsubjects),length(whichmcmcpars),nsubjects)#init[whichmcmcpars]
    mcmcstore <- array(rnorm(nrow(newmcmcpars)*nsubjects*nstore),
      dim=c(nrow(newmcmcpars),nsubjects,nstore))
    mcmclpg <- list()
  } 
  delta=deltaold=rep(0,length(pars))
  bestpars = newpars=maxpars=minpars=changepars=pars
  gstore=parstore = matrix(rnorm(length(bestpars)*nstore),length(bestpars),nstore)
  # gamweights = rep(.5,length(pars))
  
  step=rep(stepbase,length(pars))
  bestiter=1
  
  # lpg= fitfunc2(init,whichmcmcpars) # try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
  # if('try-error' %in% class(lpg)) {
  #   i = 0
  #   message('Problems initialising, trying random values...')
  #   while(i < 50 && 'try-error' %in% class(lpg)){
  #     if(i %%5 == 0) init = rep(0,length(init))
  #     init=init+rnorm(length(init),0,abs(init)+ .1)
  #     lpg= fitfunc2(init,whichmcmcpars) 
  #     i = i + 1
  #   }
  # }
  lpg=-999999
  attributes(lpg)$gradient <- rep(0,length(pars))
  
  
  g=sign(attributes(lpg)$gradient)*(abs(lpg))^(1/8)
  gsmooth=oldgsmooth=oldg=gmid=g
  ghatmix=oldghatmix=g
  mcmccount <- 0
  
  dghatsmooth=dghat=dghatmix=rep(0,length(g))
  dghatsmoothpar = rep(1,length(pars))
  dghatweight = rep(.2,length(pars)) #ideally one step length...
  dghatdownerr = dghatmixerr = dghatuperr = rep(0,length(g))
  
  lprdif = lpdif = 0
  parscore=rep(0,length(pars))
  
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
      
      # if(i > nstore && i %% gamiter == 0 && notacceptedcount<2){
      #   print(gamiter)
      #   if(i-bestiter > gamiter) gamweights<-gamweights*.1
      #   accepted <- TRUE #force acceptance
      #   lpgood <- tail(lp,nstore)
      #   # lpgood <- c(TRUE,sapply(2:nstore, function(x) lpgood[x]>=max(lpgood[1:(x-1)])))
      #   lpgood <- order(lpgood)
      #   gampars <- gamparsold <- pars
      #   gamg=gamgold=gsmooth
      #   gamindices <- sample(1:length(pars),ceiling(.1*length(pars)))
      #   if(sum(lpgood) > 10){
      #     for(pari in gamindices){
      #       dat <- data.frame(x=
      #           # log1p(-lp[1:nstore]+max(lp[1:nstore])),
      #           c(1:max(lpgood)),
      #         # g=gstore[pari,lpgood],
      #         y=c((parstore[pari,lpgood]-mean(parstore[pari,lpgood]))/sd(parstore[pari,lpgood])))
      #       weighting <- rep(1,nrow(dat))
      #       # dat <- dat[lpgood,]
      #       if(sd(parstore[pari,lpgood]) < .0001) next
      #       
      #       dat <- rbind(dat,data.frame(x=max(dat$x)*1.5,y=9999)) #for gam only
      #       weighting=c(seq(.1,1,length.out=nrow(dat)-1),0)
      #       
      #       # 
      #       # gamf <- suppressWarnings(try(gamboostLSS::gamboostLSS(formula = formula(y ~ x),
      #       #   control=boost_control(nu=.01),dfbase=4,
      #       #   # control=gamlss::gamlss.control(trace=FALSE)
      #       #   # ,sigma.formula=formula(~s(x,df=3))
      #       #   data = dat,weights = weighting),silent=TRUE))
      #       # 
      #       # if('try-error' %in% class(gamf)){
      #       #   
      #       if(1==1){
      #         gamf <- try(gam::gam(formula = formula(y ~ s(x)),
      #         control=gam::gam.control(trace=FALSE),
      #         data = dat,weights = weighting),silent=TRUE)
      #         
      #         if(!'try-error' %in% class(gamf)){
      #         gampars[pari] = 
      #           tail(predict(gamf),1)
      #         
      #         # plot(dat$x,predict(gamf),ylim=range(c(dat$y[-length(dat$y)],predict(gamf))))
      #         # points(dat$x,dat$y,col='red')
      #         }
      #         
      #       } else{
      #         gampars[pari] = 
      #           tail(predict(gamf)$mu,1)
      #         
      #         # plot(dat$x,predict(gamf)$mu,ylim=range(c(dat$y[-length(dat$y)],predict(gamf)$mu)))
      #         # points(dat$x,dat$y,col='red')
      #       }
      #       if(!'try-error' %in% class(gamf)) gampars[pari]=gampars[pari]*sd(parstore[pari,lpgood])+mean(parstore[pari,lpgood])
      #       
      #       # gamfg <- try(gam::gam(formula = formula(g ~ s(x)),
      #       #   control=gam::gam.control(trace=FALSE),
      #       #   data = dat,weights = weighting))
      #       # gamf <- lm(y ~ poly(x,5),data = dat) 
      # 
      #       # gamg[pari] = 
      #       #   tail(predict(gamfg),1)
      # 
      #     }
      #     # 
      #   }}
      
      gdelta =  (ghatmix +dghatmix*dghatweight/2) #removed step multiply because divide by zero elsewhere
      if(i > 1){
        delta =   step  *sign(gdelta)*(abs(gdelta))  #* exp((rnorm(length(g),0,.02)))
        # if(runif(1) > .95) {
        #   parextra=sample(1:length(pars),floor(.05*length(pars)))
        #   delta[parextra] <- step*sqrt(abs(gsmooth[parextra]))*10
        # }
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        newpars = pars + delta
        newpars = newpars  + delta/2 - deltaold/2 #+ delta - deltaold #
      }
      
      # if(i > nstore && (i%%gamiter)==0){
      #   # gampars = newpars*(1-min(1,gamweights))+gampars*gamweights
      #   newpars[gamindices]=newpars[gamindices] + 
      #     (gampars[gamindices]-gamparsold[gamindices])*gamweights[gamindices]
      #   # gsmooth=gsmooth*.8 + (gampars-gamparsold)/step*.2#(gamg-gamgold)*gamweights
      # }
      
      if(any(is.na(newpars))) browser() 
        if(i==1) itertime <- Sys.time()
      
      
      if(!mcmcconverged){
        mcmcpars <- newmcmcpars #store current mcmcpars
        # 
        propchol <- t(chol(matrix(
          apply(
            apply(aperm(mcmcstore,c(3,2,1)),2,cov),
            1,mean),
          dim(mcmcstore)[1])))
        fitout <- mcmcfitfunc(newpars,newmcmcpars,mcmcconverged)
        lpg <- fitout$lpg
        newmcmcpars <- fitout$newmcmcpars
        mcmcstore[,,(i-1)%%nstore+1] <- newmcmcpars
        # ar = ar*.95+.05*i/mcmccount
        mcmcstep = mcmcstep * (1-2*( ( (1/(-ar-artarget)) / (1/-artarget) + .5) -1))
        # message(paste0('ar = ',ar,'   mcmcstep = ',mcmcstep))
      }
      if(!is.na(whichmcmcpars[1]) && mcmcconverged){
        lpg= mcmcfitfunc(newpars,newmcmcpars,mcmcconverged)$lpg #fitfunc2(combinepars(newpars,newmcmcpars,whichmcmcpars),whichmcmcpars)
      }
      
      fullnewpars <- initfull
      if(length(whichignore)>0) fullnewpars[-whichignore] <- newpars else fullnewpars <- newpars
      
      
      if(parsets > 1){
        parmod <- (exp(1*((1:parsets)-(floor(parsets/2)))))
        fullnewpars <- lapply(parmod,function(x) {
          newpars <-  pars+x*
            # sample(c(rep(1,parsets),parmod),length(newpars),replace = TRUE) *
            (newpars-pars)
          fullnewpars <- initfull
          if(length(whichignore)>0) fullnewpars[-whichignore] <- newpars else fullnewpars <- newpars
          return(fullnewpars)
        })
      }
      
      if(is.na(whichmcmcpars[1])) lpg= fitfunc2(fullnewpars,whichmcmcpars)
      if(length(whichignore)>0) attributes(lpg)$gradient <- attributes(lpg)$gradient[-whichignore]
      if(!is.null(attributes(lpg)$bestset)){
        bestset <- attributes(lpg)$bestset
        if(bestset > ceiling(parsets/2)) lproughnesstarget <- lproughnesstarget +.01
        if(bestset < floor(parsets/2)) lproughnesstarget <- lproughnesstarget -.01
        newpars <- fullnewpars[[bestset]]
        if(length(whichignore)>0) newpars <- newpars[-whichignore]
      }
      
      
      # if(i==1) {
      #   itertime <- as.numeric(Sys.time()-itertime)
      #   gamiter <- max(ceiling(10/itertime),gamiter)#
      #   #/30.01*length(pars)
      # }
      
      if(lpg > -1e99 &&       #regular check
          class(lpg) !='try-error' && 
          !is.nan(lpg[1]) && 
          all(!is.nan(attributes(lpg)$gradient)) &&
          (i < warmuplength || (!mcmcconverged || (lp[i-1]- lpg[1]) < sd(tail(lp,100))*8+1e-3))
      ){
        accepted <- TRUE
      } 
      else {
        ghatmix=g*.5
        dghatmix=dghatmix*.5
        # gsmooth= gsmooth*gmemory2 + (1-gmemory2) * g #increase influence of last gradient at inflections
        step <- step * .5
        deltaold <- deltaold * .5
        # pars=bestpars
      }
      #warmup check
      if(mcmcconverged && is.na(startnrows) && 
          i < warmuplength && i > 1 && lpg[1] < lp[1]-5) {
        accepted <- FALSE
        step = step * .1
        deltaold <- deltaold * .1
        pars=bestpars
        ghatmix=.1*ghatmix
        dghatmix=.1*dghatmix
        
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
    g=sign(g)*(abs(g))^(1/2)#sqrt
    gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
    roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)
    
    
    oldgmid=gmid
    gmid = (oldg+g)/2
    dg=(gmid-oldgmid)#/step #removed step divide due to divide by zero problems
    
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
    
    
    parstore[,1+(i-1) %% nstore] = pars
    gstore[,1+(i-1) %% nstore] = g
    
    
    
    # if(i > (nstore+10) && (i %% gamiter)==10){# == (gamiter-2)) {
    #   # 
    #   print((lp[i]-lp[i-12]) > (lp[i-12]-lp[i-24]))
    #   gamup <- sign(gampars[gamindices]-gamparsold[gamindices]) == 
    #     sign(bestpars[gamindices] - gampars[gamindices])
    #   if((lp[i]-lp[i-12]) > (lp[i-12]-lp[i-24])) {
    #     gamweights[gamindices] <- gamweights[gamindices] * 1.1 
    #   } else{
    #     gamweights[gamindices] = gamweights[gamindices] * .9
    #     # gamiter=gamiter*1.5 #increase time between gams
    #   }
    
    # gamweights[gamindices][gamup] <- gamweights[gamindices][gamup] * 1.2
    # gamweights[gamindices][!gamup] <- gamweights[gamindices][!gamup] * .8
    # }
    
    # if(!i %% gamiter==0){
    groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gmid)!=sign(oldgmid))
    gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
    if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] > (lp[i]))#exp(-1/(i-bestiter+.1))
    
    lproughnessmod=  ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    gsmoothroughnessmod =  (( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    groughnessmod = ( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    
    step = (step + roughnesschangemulti*(
      step* .8*lproughnessmod
      # + step* .1*gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      + step* .3*groughnessmod# * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    ))
    
    step[gsmoothroughness < gsmoothroughnesstarget] <- step[gsmoothroughness < gsmoothroughnesstarget] * 1.2
    gsmooth[gsmoothroughness < gsmoothroughnesstarget] <- gsmooth[gsmoothroughness < gsmoothroughnesstarget] * 1.1
    # step[gsmoothroughness < gsmoothroughnesstarget] * .1*gsmoothroughnessmod[gsmoothroughness < gsmoothroughnesstarget]
    signdif= sign(gmid)!=sign(gdelta)
    if(i > 1 && lp[i] >= max(head(lp,length(lp)-1))) {
      step = step * 1.1 #sqrt(2-gmemory) #exp((1-gmemory)/8)
      if(i > warmuplength) {
        ##max/min par update extra
        parscore <- parscore * .98
        whichmax <- which(pars > maxpars | pars < minpars)
        if(length(whichmax) > 0){
          parscore[whichmax] <- parscore[whichmax]+.1*(as.numeric(pars[whichmax]>maxpars[whichmax])*2-1)
          # gsmooth[whichmax] <- gsmooth[whichmax]  * 1.2*(1+abs(parscore[whichmax]))#*delta[whichmax] /step[whichmax]
          # step[whichmax] <- step[whichmax] * 2*(1+abs(parscore[whichmax]))  #+ pars[whichmax]
          pars[pars>maxpars] <- pars[pars>maxpars]+10*(1+abs(parscore[pars>maxpars]))*(pars[pars>maxpars]-maxpars[pars>maxpars] )
          pars[pars< minpars] <- pars[pars< minpars]+10*(1+abs(parscore[pars<minpars]))*(pars[pars< minpars]-minpars[pars< minpars] )
          
          maxpars[pars>maxpars] <-pars[pars>maxpars]
          minpars[pars<minpars] <-pars[pars<minpars]
        }
        changepars=pars
        if(length(whichmax)) changepars[-whichmax] <- NA else changepars[]<-NA
        # lproughness = lproughness * .9
      }
      # gmemory <- gmemory +(1-gmemory)*1.001
      bestpars <- pars <- newpars
      bestg <- g
      bestiter <- i
      
    }
    
    # }#end if not gam
    
    if(i > 1 && runif(1,0,1) > .95) {
      # #slowly forget old max and mins, allow fast re exploration of space
      rndchange <- runif(length(maxpars),0,1) > .95
      # step[rndchange] <- stepbase
      if(any(rndchange)){
        maxpars[rndchange] <- max(parstore[rndchange,]+1e-6)
        minpars[rndchange] <- min(parstore[rndchange,]-1e-6)
      }
    }
    
    # gmemory <- gmemory * gsmoothroughnessmod
    if(i > 45 && i %% 40 == 0) {
      oldlpdif <- lpdif# sum(diff(head(tail(lp,10),20)))
      sublp <- tail(lp,45)
      lpdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      if(oldlpdif > lpdif) gmemory <- oldgmemory
      proposal = gmemory*2-oldgmemory
      oldgmemory <- gmemory
      gmemory <- min(gmemmax, max(0, proposal + runif(1,-.025,.05)))
      if(gmemory < .9) gmemory <- gmemory + .02
    }
    
    if(i > 31 && i %% 30 == 0 && is.na(whichmcmcpars[1])) {
      oldlprdif <- lprdif
      sublp <- tail(lp,35)
      lprdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      if(oldlprdif > lprdif) lproughnesstarget <- oldlproughnesstarget
      lprproposal = lproughnesstarget*2-oldlproughnesstarget
      oldlproughnesstarget <- lproughnesstarget
      lproughnesstarget <- min(.7, max(.05, lprproposal + .025 * (-1+2*rbinom(n = 1,size = 1,prob = .5))))
    }
    
    step[step > maxparchange] <- maxparchange
    step[step < minparchange] <- minparchange
    
    if(i > warmuplength && lp[i] < lp[i-1] && mcmcconverged) { #if worsening, update gradient faster
      step[signdif]=step[signdif]*lproughnesstarget
      if(lp[i] < lp[i-10]) gmemory <- gmemory * .995
      # step=step*.5
      gsmooth[signdif]= gsmooth[signdif]*gmemory2 + (1-gmemory2) * g[signdif] #increase influence of gradient at inflections
    }
    
    oldghatmix=ghatmix
    ghatmix = gsmooth #ghatsmoothpar*gsmooth + (1-ghatsmoothpar)*g
    # message('ghatsmoothpar = ',ghatsmoothpar)
    # print(as.numeric(plot))
    if(plot && i %% as.numeric(plot) ==0){
      par(mfrow=c(2,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
      plot(pars,col=1:length(pars))
      points(changepars,pch=17,col='red')
      
      plot(log(abs(step*gsmooth)),col=1:length(pars))
      plot(tail(log(-(lp-max(lp)-1)),500),type='l')
      # plot(gamweights,col=1:length(pars))
      parsd=(apply(parstore,1,sd,na.rm=T))
      plot(pars,col=1:length(pars))
      abline(h=(parsdtol))
      matplot(t(parstore[
        which(parsd > sort(parsd,decreasing = TRUE)[min(c(length(pars),5))]),,drop=FALSE]),
        type='l')
      if(1==2){
        plot(groughness,col='red',ylim=c(0,1))
        abline(h=mean(gsmoothroughness),col='blue',lty=2)
        abline(h=(gsmoothroughnesstarget),col='blue',lty=1,lwd=2)
        points(gsmoothroughness,ylim=c(0,1),col='blue')
        abline(h=mean(groughness),col='red',lty=2)
        abline(h=(groughnesstarget),col='red',lty=1)
        
        abline(h=lproughnesstarget,lty=1,col='green')
        abline(h=lproughness, col='green',lty=2)
        
        
        # gsmoothsqrt=sign(gsmooth) * sqrt(abs(gsmooth))
        # plot(gsmoothsqrt,ylim=c(-max(abs(gsmoothsqrt)),max(abs(gsmoothsqrt))))
        
        # plot(dghatmix*step,ylim=c(-max(abs(dghatmix*step)),max(abs(dghatmix*step))))
        
        plot(ghatsmoothpar,ylim=c(0,1))
        if(!is.na(whichmcmcpars[1])) plot(c(newmcmcpars))
      }
      Sys.sleep(.03)
      # plot(dghatsmoothpar)
      # abline(h=mean(dghatsmoothpar))
      # plot(dghatweight)
      # abline(h=mean(dghatweight))
      
      # matplot(cbind(signdifmod,gsmoothroughnessmod),col=c('black','blue'),pch=1,ylim=c(-1,1))
      # points(groughnessmod,col='red')
      # abline(h=lproughnessmod,col='green')
      
      message(paste0('Iter = ',i, '   LP = ', (lp[i]),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory,'  lprt = ',lproughnesstarget))
    }
    
    #check convergence
    if(i > 30){
      if(is.na(whichmcmcpars[1])){
        if( (i - bestiter) > nconvergeiter*5 && 
            mean(sign(diff(tail(lp,nconvergeiter)))) < .3) converged <- TRUE #time since best
        if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
        if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
        if(max(apply(parstore,1,sd)) < parsdtol) converged <- TRUE
      } else {
        mu1=apply(parstore[,1:ceiling(nstore/2),drop=FALSE],1,mean,na.rm=TRUE)
        mu2=apply(parstore[,(1+ceiling(nstore/2)):nstore,drop=FALSE],1,mean,na.rm=TRUE)
        sd1=apply(parstore[,1:ceiling(nstore/2),drop=FALSE],1,sd,na.rm=TRUE)
        sd2=apply(parstore[,(1+ceiling(nstore/2)):nstore,drop=FALSE],1,sd,na.rm=TRUE)
        plot(abs((mu1-mu2))/sd2,ylim=c(0,max(c(abs((mu1-mu2))/sd2,abs(sd1-sd2),sd2))))
        points(abs(sd1-sd2),col='red')
        points(sd2,col='green')
        if( all(abs((mu1-mu2))/sd2 < itertol) && abs(sd1-sd2) < itertol) converged <- TRUE
      }
    }
    if(converged && !mcmcconverged){
      message('mcmc converged')
      # bestpars=apply(parstore,1,mean,na.rm=TRUE)
      
      # bestiter<-iter<-i
      # newmcmcpars <- apply(mcmcstore,c(1,2),mean,na.rm=TRUE)
      # converged <- FALSE
      # mcmcconverged <- TRUE
    }
  }
  # if(!is.na(whichmcmcpars[1])) bestpars = combinepars(bestpars,newmcmcpars,whichmcmcpars)
  out=list(itervalues = lp, value = max(lp),
    par=bestpars,parstore=parstore,gstore=gstore,lpstore=tail(lp,nstore))
  if(!is.na(whichmcmcpars[1])){
    out$mcmcpars=newmcmcpars
    out$bestpars = mu2
    out$sd = sd2#apply(parstore,1,mean,na.rm=TRUE)
  }
  return(out)#,gstore=gstore,pstore=pstore) )
}

