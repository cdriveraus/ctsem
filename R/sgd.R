logit = function(x) log(x)-log((1-x))

sgd <- function(init,fitfunc,whichignore=c(),nsubjects=NA,ndatapoints=NA,plot=FALSE,
  stepbase=1e-3,gmeminit=ifelse(is.na(startnrows),.9,.8),gmemmax=.95, maxparchange = .50,
  startnrows=NA,roughnessmemory=.9,groughnesstarget=.4,roughnesschangemulti = 2,
  lproughnesstarget=ifelse(parsets==1,.4,.2),parsets=1,
  gsmoothroughnesstarget=.05,
  warmuplength=20,nstore=max(100,length(init)),
  minparchange=1e-800,maxiter=50000,
  nconvergeiter=30, 
  itertol=1e-3, deltatol=1e-5, parsdtol=1e-3){
  
  initfull=init #including ignored params start values
  if(length(whichignore)>0) init=init[-whichignore]
  
  errsum = function(x) sqrt(sum(abs(x)))
  
  if(plot){
    parbase=par(no.readonly=TRUE)
    on.exit(do.call(par,parbase),add=TRUE)
  }
  
  
  pars=init

  delta=deltaold=rep(0,length(pars))
  bestpars = newpars=maxpars=minpars=changepars=pars
  gstore=parstore = deltastore=matrix(rnorm(length(bestpars)*nstore),length(bestpars),nstore)
  
  step=rep(stepbase,length(pars))
  bestiter=1
  
  lpg=fitfunc(init)#-999999
  # attributes(lpg)$gradient <- rep(0,length(pars))
  
  
  g=sign(attributes(lpg)$gradient)*(abs(lpg))^(1/8)
  gsmooth=oldgsmooth=oldg=gmid=dgsmooth=gdelta=g
  
  lprdif = lpdif = 0
  parscore=rep(0,length(pars))
  
  groughness = rep(groughnesstarget,length(g))
  gsmoothroughness = rep(gsmoothroughnesstarget,length(g))
  lproughness=oldlproughnesstarget=lproughnesstarget
  gmemory <- gmeminit
  oldgmemory  <- gmemory
  oldlpdif <- 0
  lpdif <- 0
  maxlp <- -Inf
  i=0
  lp<-c()
  oldlp <- -Inf
  converged <- FALSE
  S<-1
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
      
      
     
      if(i > 1){
        # browser()
        delta =   step * (gsmooth+dgsmooth/2) 
        # deltadgdp= step * dgsmooth * sum(abs(gsmooth))/sum(abs(dgsmooth))# +dgsmooth *sum(abs(gsmooth))/sum(abs(dgsmooth)))   #* exp((rnorm(length(g),0,.02)))
        # delta=delta+deltadgdp
        
        if(rbinom(1,1,.05)==1) delta = delta * exp((rnorm(length(g),0,.5)))
        s <- which(rbinom(n = length(pars),size = 1,prob= .1)==1)
        # delta[s] <- delta[s] + sign(newpars[s])*abs(rnorm(length(s),0,abs(gdelta[s])*step[s]))
        # s <- which(rbinom(n = length(pars),size = 1,prob= .1)==1)
        # delta[s] <- delta[s] + rnorm(length(s),0,abs(gdelta[s])*step[s]*2)
        
        # if(runif(1) > .95) {
        #   parextra=sample(1:length(pars),floor(.05*length(pars)))
        #   delta[parextra] <- step*sqrt(abs(gsmooth[parextra]))*10
        # }
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        delta = delta +  delta/2 - deltaold/2
        
        newpars = pars + delta 

        
        # #random jumping
        # if(i %% 10 ==0){
        #   s <- which(rbinom(n = length(pars),size = 1,prob= .8)==1)
        #   newpars[s] <- newpars[s] + rnorm(length(s),0,abs(delta[s])*1)
        # }

      }
      
      
      if(any(is.na(newpars))) browser() 
      if(i==1) itertime <- Sys.time()
      
      
      
      
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
      
      lpg= fitfunc(fullnewpars)
      if(length(whichignore)>0) attributes(lpg)$gradient <- attributes(lpg)$gradient[-whichignore]
      if(!is.null(attributes(lpg)$bestset)){
        bestset <- attributes(lpg)$bestset
        if(bestset > ceiling(parsets/2)) lproughnesstarget <- lproughnesstarget +.01
        if(bestset < floor(parsets/2)) lproughnesstarget <- lproughnesstarget -.01
        newpars <- fullnewpars[[bestset]]
        if(length(whichignore)>0) newpars <- newpars[-whichignore]
      }
      
      
      if(lpg > -1e99 &&       #regular check
          class(lpg) !='try-error' && 
          !is.nan(lpg[1]) && 
          all(!is.nan(attributes(lpg)$gradient)) &&
          (i < warmuplength || ( exp(lpg[1] - lp[i-1]) > runif(1,0,1))) #sd(tail(lp,100))*8+
      ){
        accepted <- TRUE
      } 
      else {
        gsmooth= gsmooth*gmemory2^2 + (1-gmemory2^2) * g #increase influence of last gradient at inflections
        step <- step * .5
        deltaold <- deltaold * .5
        # if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) ##exp(-1/(i-bestiter+.1))
        # pars=bestpars
      }
      #warmup check
      if(is.na(startnrows) && 
          i < warmuplength && i > 1 && lpg[1] < lp[1]-5) {
        accepted <- FALSE
        step = step * .1
        deltaold <- deltaold * .1
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
    g=sign(g)*(abs(g))^(1/2)#sqrt
    gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
    roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)
    
    
    oldgmid=gmid
    gmid = g#(oldg+g)/2
    # dg=(g-oldg)#/step #removed step divide due to divide by zero problems
    
    
    #predictions
    oldgsmooth = gsmooth
    gsmooth= gsmooth*gmemory2 + (1-gmemory2) * g 
    dgsmooth = gmemory2*dgsmooth +(1-gmemory2)*(gsmooth-oldgsmooth)
    
    
    parstore[,1+(i-1) %% nstore] = pars
    gstore[,1+(i-1) %% nstore] = g
    

    if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] > (lp[i]))#because accepted here, also see non accepted version
    groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gmid)!=sign(oldgmid))
    gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
    
    lproughnessmod=  ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    gsmoothroughnessmod =  (( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    groughnessmod = ( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    
    step = (step + roughnesschangemulti*(
      step* .6*lproughnessmod
      # + step* .05*gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      + step* .6*groughnessmod# * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    ))
    
    signdif= sign(gsmooth)!=sign(gmid)
    if(i > warmuplength) step[gsmoothroughness < gsmoothroughnesstarget & !signdif] <- step[gsmoothroughness < gsmoothroughnesstarget& !signdif] *2
    
    # gsmooth[gsmoothroughness < gsmoothroughnesstarget] <- gsmooth[gsmoothroughness < gsmoothroughnesstarget] * 1.2
    # step[gsmoothroughness < gsmoothroughnesstarget] * .1*gsmoothroughnessmod[gsmoothroughness < gsmoothroughnesstarget]
    
    if(i > 1 && lp[i] >= max(head(lp,length(lp)-1))) {
      # step[!signdif] = step[!signdif] * 1.5 #sqrt(2-gmemory) #exp((1-gmemory)/8)
      # step = step * 1.2
      if(i > warmuplength) {
        ##max/min par update extra
        parscore <- parscore * .98
        whichmax <- which(pars > maxpars | pars < minpars)
        if(length(whichmax) > 0){
          parscore[whichmax] <- parscore[whichmax]+.1*(as.numeric(pars[whichmax]>maxpars[whichmax])*2-1)
          # gsmooth[whichmax] <- gsmooth[whichmax]  * 1.2*(1+abs(parscore[whichmax]))#*delta[whichmax] /step[whichmax]
          # step[whichmax] <- step[whichmax] * 2*(1+abs(parscore[whichmax]))  #+ pars[whichmax]
          # pars[pars>maxpars] <- pars[pars>maxpars]+10*(1+abs(parscore[pars>maxpars]))*(pars[pars>maxpars]-maxpars[pars>maxpars] )
          # pars[pars< minpars] <- pars[pars< minpars]+10*(1+abs(parscore[pars<minpars]))*(pars[pars< minpars]-minpars[pars< minpars] )

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
    if(i > 25 && i %% 20 == 0) {
      oldlpdif <- lpdif# sum(diff(head(tail(lp,10),20)))
      sublp <- tail(lp,20)
      lpdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      if(oldlpdif > lpdif) gmemory <- oldgmemory
      proposal = gmemory*2-oldgmemory
      oldgmemory <- gmemory
      gmemory <- min(gmemmax, max(0, proposal + runif(1,-.025,.05)))
      if(gmemory < .95) gmemory <- gmemory + .02
    }
    
    if(i > 31 && i %% 30 == 0) {
      oldlprdif <- lprdif
      sublp <- tail(lp,15)
      lprdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      if(oldlprdif > lprdif) lproughnesstarget <- oldlproughnesstarget
      lprproposal = lproughnesstarget*2-oldlproughnesstarget
      oldlproughnesstarget <- lproughnesstarget
      if(max(lp) > max(tail(lp,30))) lprproposal <- min(.2,.5 * lprproposal)
      lproughnesstarget <- min(.5, max(.2, lprproposal + .025 * (-1+2*rbinom(n = 1,size = 1,prob = .5))))
      
    }
    
    step[step > maxparchange] <- maxparchange
    step[step < minparchange] <- minparchange
    
    if(i > warmuplength && lp[i] < lp[i-1]) { #if worsening, update gradient faster
      step[signdif]=step[signdif]*.5#lproughnesstarget
      # step = step * lproughnesstarget
      if(lp[i] < lp[i-10]) gmemory <- gmemory * .995
      # step=step*.5
      gsmooth[signdif]= gsmooth[signdif]*gmemory2^2 + (1-gmemory2^2) * g[signdif] #increase influence of gradient at inflections
      # gsmooth[signdif]= gsmooth[signdif]*.5 + .5 * g[signdif] #increase influence of gradient at inflections
    }
    
    if(plot && i %% as.numeric(plot) ==0){
      par(mfrow=c(2,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
      plot(pars,col=1:length(pars))
      points(changepars,pch=17,col='red')
      plot(log(abs(step*gsmooth)+1e-50),col=1:length(pars))
      plot(tail(log(-(lp-max(lp)-1)),500),type='l')
      # plot(gamweights,col=1:length(pars))
      parsd=(apply(parstore,1,sd,na.rm=T))
      abline(h=(parsdtol))
      matplot(t(parstore[
        which(parsd > sort(parsd,decreasing = TRUE)[min(c(length(pars),5))]),,drop=FALSE]),
        type='l')
      if(1==1){
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
        # plot(gsmooth,ylim=c(-max(abs(gsmooth)),max(abs(gsmooth))))
        
      }
      Sys.sleep(.03)
      
      
      # matplot(cbind(signdifmod,gsmoothroughnessmod),col=c('black','blue'),pch=1,ylim=c(-1,1))
      # points(groughnessmod,col='red')
      # abline(h=lproughnessmod,col='green')
      
      message(paste0('Iter = ',i, '   LP = ', (lp[i]),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory,'  lprt = ',lproughnesstarget))
    }
    
    #check convergence
    if(i > 30 && max(tail(lp,nconvergeiter)) ==max(lp)){
      # if( (i - bestiter) > nconvergeiter*5 && 
      #     mean(sign(diff(tail(lp,nconvergeiter)))) < .3) converged <- TRUE #time since best
      lpdiff=max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter))
      if(lpdiff < itertol & lpdiff > 0) converged <- TRUE
      if(abs(max(diff(tail(lp,nconvergeiter)))) < deltatol) converged <- TRUE
      if(max(apply(parstore,1,sd)) < parsdtol) converged <- TRUE
      # if(converged) browser()
    }
  }
  out=list(itervalues = lp, value = max(lp),
    par=bestpars,parstore=parstore,gstore=gstore,lpstore=tail(lp,nstore))
  
  return(out)#,gstore=gstore,pstore=pstore) )
}

