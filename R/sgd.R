logit = function(x) log(x)-log((1-x))

sgd <- function(init,fitfunc,whichignore=c(),nsubsets=1,nsubjects=NA,plot=FALSE,
  stepbase=1e-3,gmeminit=ifelse(is.na(startnrows),.8,.8),gmemmax=.95, maxparchange = .50,
  startnrows=NA,roughnessmemory=.9,groughnesstarget=.5,roughnesschangemulti = .2,
  parsets=1,
  lproughnesstarget=ifelse(parsets==1,.2,.5),
  gsmoothroughnesstarget=.1,
  warmuplength=20,nstore=100,
  minparchange=1e-800,maxiter=50000,
  nconvergeiter=30, 
  worsecountconverge=1000,
  lpnoisethresh= .1,#length(init)*.01,
  itertol=1e-3, parrangetol=1e-3){
  
  
  if(nsubsets>1){
    lpchange=0
    progressEst=0
    oldsubsetilp <- -Inf
    bestsubsetlp <- -Inf
    plot=plot*nsubsets
    worsecount = 0
    warmuplength = warmuplength * nsubsets
    nstore=nsubsets
    oldsubsetlp <- -1e100
    gmeminit=.8
    stepbase=1e-3
    maxiter=maxiter*nsubsets
    lproughnesstarget=.45#1/nsubsets + .05
    gsmoothroughnesstarget=.05
    groughnesstarget=.5
    roughnesschangemulti= .1#0r .5?
    roughnessmemory=.9
    subsetorder=1:nsubsets
    subsetlps <- rep(-Inf, nsubsets)
    subsetlpstore <- c()
  }
  
  initfull=init #including ignored params start values
  if(length(whichignore)>0) init=init[-whichignore]
  
  fitfunc2 <- function(p,subset=NA){
    ptemp=initfull
    if(length(whichignore)>0) ptemp[-whichignore] <- p else ptemp <- p
    if(!is.na(subset)) ptemp <- c(ptemp,subset)
    lpg=fitfunc(ptemp)
    if(length(whichignore)>0) attributes(lpg)$gradient <- attributes(lpg)$gradient[-whichignore]
    if(nsubsets > 1) attributes(lpg)$gradient <- 
      attributes(lpg)$gradient[-length(attributes(lpg)$gradient)]
    return(lpg)
  }
  
  errsum = function(x) sqrt(sum(abs(x)))
  
  if(plot){
    parbase=par(no.readonly=TRUE)
    on.exit(do.call(par,parbase),add=TRUE)
  }
  
  pars=init
  zeroG=0
  delta=deltaold=rep(0,length(pars))
  bestpars = newpars=maxpars=minpars=changepars=pars
  gstore=parstore = deltastore=matrix(rnorm(length(bestpars)*nstore),length(bestpars),nstore)
  
  step=rep(stepbase,length(pars))
  bestiter=1
  
  
  lpg=fitfunc2(init,subset=ifelse(nsubsets>1,1,NA))#-999999
  # if(any(is.na(attributes(lpg)$gradient))) browser()
  
  g=sign(attributes(lpg)$gradient)*(abs(lpg))^(1/2)
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
  converged <- 0
  S<-1
  while(converged < 1 && i < maxiter){
    i = i + 1
    accepted <- FALSE
    lproughnesstarget2 = lproughnesstarget 
    notacceptedcount <- 0
    
    while(!accepted){
      notacceptedcount <- notacceptedcount+1
      if(notacceptedcount > 20) {
        warning('Cannot optimize stochastically! Problematic model, or bug?')
        print(lpg)
        accepted<-TRUE
      }
      
      if(i > 1){
        delta =   step * sign(gsmooth+dgsmooth/4) * (abs(gsmooth+dgsmooth/4))^(1/2)
        # deltadgdp= step * dgsmooth * sum(abs(gsmooth))/sum(abs(dgsmooth))# +dgsmooth *sum(abs(gsmooth))/sum(abs(dgsmooth)))   #* exp((rnorm(length(g),0,.02)))
        # delta=delta+deltadgdp
        
        
        # s <- which(rbinom(n = length(pars),size = 1,prob= .1)==1)
        # if(nsubsets > 1 && rbinom(1,1,.05)==1) delta[s] = delta[s] * exp((rnorm(length(g[s]),0,.5)))
        # delta[s] <- delta[s] + sign(newpars[s])*abs(rnorm(length(s),0,abs(gdelta[s])*step[s]))
        # s <- which(rbinom(n = length(pars),size = 1,prob= .1)==1)
        # delta[s] <- delta[s] + rnorm(length(s),0,abs(gdelta[s])*step[s]*2)
        
        # if(runif(1) > .95) {
        #   parextra=sample(1:length(pars),floor(.05*length(pars)))
        #   delta[parextra] <- step*sqrt(abs(gsmooth[parextra]))*10
        # }
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        # delta = delta +  delta/2 - deltaold/2
        
        newpars = pars + delta
        
        
        #random jumping
        # if(nsubsets==1 && i %% ceiling(runif(1,1,30)) ==0){
        #   s <- which(rbinom(n = length(pars),size = 1,prob= .8)==1)
        #   newpars[s] <- newpars[s] + rnorm(length(s),0,abs(delta[s])*1)
        # }
        
      }
      
      
      if(any(is.na(newpars))) stop('NA in parameter proposal!') 
      if(i==1) itertime <- Sys.time()
      
      fullnewpars <- newpars
      
      if(parsets > 1){
        parmod <- (exp(1*((1:parsets)-(floor(parsets/2)))))
        fullnewpars <- lapply(parmod,function(x) {
          newpars <-  pars+x*(newpars-pars)
          return(newpars)
        })
      }
      
      if(nsubsets > 1){
        if(i > 1) subsetiold <- subseti
        subseti <- subsetorder[(i-1) %% (nsubsets) + 1]
      }

      lpg= fitfunc2(fullnewpars,subset=ifelse(nsubsets > 1, subseti, NA)) #fit function!
      
      if(!is.null(attributes(lpg)$bestset)){
        bestset <- attributes(lpg)$bestset
        if(bestset > ceiling(parsets/2)) lproughnesstarget <- lproughnesstarget +.01
        if(bestset < floor(parsets/2)) lproughnesstarget <- lproughnesstarget -.01
        newpars <- fullnewpars[[bestset]]
      }
      
      
      if(lpg > -1e99 &&       #regular check
          !'try-error' %in% class(lpg) && 
          !is.nan(lpg[1]) && 
          all(!is.nan(attributes(lpg)$gradient)) &&
          (nsubsets > 1 || i ==1 || lpg[1] > (min(tail(lp,20))-notacceptedcount-1))  #no subset lp check
        # (i < warmuplength || ( exp(lpg[1] - lp[i-1]) > runif(1,0,1))) #sd(tail(lp,100))*8+
      ){
        accepted <- TRUE
      } 
      if(!accepted){
        if(nsubsets==1 && i > warmuplength) gsmooth= gsmooth*gmemory2^2 + (1-gmemory2^2) * g #increase influence of last gradient at inflections
        step <- step / (exp(notacceptedcount))
        deltaold <- 0#deltaold * 0
        if(nsubsets > 1) pars =bestpars* .8 + apply(parstore,1,mean)*.2
        # if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) ##exp(-1/(i-bestiter+.1))
        # pars=bestpars
      }
      #warmup check
      if(is.na(startnrows) && nsubsets==1 && 
          i < warmuplength && i > 1 && lpg[1] < lp[1]-5) {
        accepted <- FALSE
        step = step * .1
        deltaold <- deltaold * 0
        pars=bestpars
        
      }
      # if(plot && !accepted) {
      #   message(paste0('\rIter ', i,' not accepted! lp = ', lpg[1],', continuing...                     '))
      #   # 
      # }
    } #end acceptance loop
    
    #once accepted 
    lp[i]=lpg[1]
    pars=newpars
    
    
    
    
    deltaold=delta
    oldg=g
    g=attributes(lpg)$gradient
    if(nsubsets == 1 && zeroG < 11){ #check zero gradient count
      if(any(g %in% 0)){
        message(paste0('Gradient of parameter ',paste0(which(g %in% 0),collapse=', '), ' is exactly zero -- if this repeats or optimization fails, maybe model problem?'))
        zeroG <- zeroG + 1
        if(zeroG >  10) message ('Gradient warning limit exceeded, silencing...')
      }
    }
    if(any(is.na(g))) warning(paste0('Gradient of parameter ',paste0(which(is.na(g)),collapse=', '), ' is NA, if this repeats or optimization fails, maybe model problem?'))
    # g=sign(g)*(abs(g))#^(1/2)#sqrt
    gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
    roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)
    
    
    oldgmid=gmid
    gmid = (oldg+g)/2
    # dg=(g-oldg)#/step #removed step divide due to divide by zero problems
    
    
    #predictions
    oldgsmooth = gsmooth
    gsmooth= gsmooth*gmemory2 + (1-gmemory2) * g 
    dgsmooth = gmemory2*dgsmooth +(1-gmemory2)*(gsmooth-oldgsmooth)
    
    
    parstore[,1+(i-1) %% nstore] = pars
    gstore[,1+(i-1) %% nstore] = g
    
    
    if(i > 1) {
      if(nsubsets==1) lproughness = lproughness * (roughnessmemory2) + 
          (1-(roughnessmemory2)) * as.numeric(lp[i-1] > (lp[i]))#because accepted here, also see non accepted version
      
      if(nsubsets > 1) lproughness = lproughness * (roughnessmemory2) + 
          (1-(roughnessmemory2)) * as.numeric( (lp[i-1]-oldsubsetilp) > (lp[i]-subsetlps[subseti]))#because accepted here, also see non accepted version
    }
    # if(i==101) 
    if(nsubsets > 1){
      oldsubsetilp <- subsetlps[subseti]
      subsetlps[subseti] <- lpg[1]
    }
    
    
    
    groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gmid)!=sign(oldgmid))
    gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
    
    lproughnessmod=  ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
    
    gsmoothroughnessmodmax=(( ( (1/(-(0)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)^4
    gsmoothroughnessmod =  (( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
    gsmoothroughnessmod=(gsmoothroughnessmod^4*sign(gsmoothroughnessmod))/gsmoothroughnessmodmax
    # gsmoothroughnessmod =  (max(c(0,gsmoothroughnesstarget-gsmoothroughness))^4)/gsmoothroughnesstarget^4
    groughnessmod = ( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
    
    if(i > warmuplength) {
      if(nsubsets==1) step = step + roughnesschangemulti*(step* .5*lproughnessmod) #(ifelse(nsubsets > 0,0,1))
      step = step + roughnesschangemulti*step* 2*gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
      if(nsubsets==1) step = step + roughnesschangemulti*step* .6*groughnessmod # * min(sqrt(deltasmoothsq),1)
      # + step * rmsstepmod
    }
    
    signdif= sign(gsmooth)!=sign(gmid)
    # if(i > warmuplength*5) 
    # step[gsmoothroughness < gsmoothroughnesstarget & !signdif] <-
    # step[gsmoothroughness < gsmoothroughnesstarget& !signdif] *(1+roughnesschangemulti*.01)
    
    # gsmooth[gsmoothroughness < gsmoothroughnesstarget] <- gsmooth[gsmoothroughness < gsmoothroughnesstarget] * 1.2
    # step[gsmoothroughness < gsmoothroughnesstarget] * .1*gsmoothroughnessmod[gsmoothroughness < gsmoothroughnesstarget]
    
    if(i > 1 && lp[i] >= max(head(lp,length(lp)-1))) {
      # step[!signdif] = step[!signdif] * 1.5 #sqrt(2-gmemory) #exp((1-gmemory)/8)
      # step = step * 1.05
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
    if(nsubsets==1 &&i > 25 && i %% 20 == 0) {
      # oldlpdif <- lpdif# sum(diff(head(tail(lp,10),20)))
      # sublp <- tail(lp,20)
      # lpdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      # if(oldlpdif > lpdif) gmemory <- oldgmemory
      # proposal = gmemory*2-oldgmemory
      # oldgmemory <- gmemory
      # gmemory <- min(gmemmax, max(0, proposal + runif(1,-.025,.1)))
      gmemory=min(c(gmemmax,max(c(.8,gmemory+runif(1,-.01,.01)))))
      # if(gmemory > .95) gmemory = .95
    }
    
    if(nsubsets==1 && i %% 10 == 0) {
      # oldlprdif <- lprdif
      # sublp <- tail(lp,15)
      # lprdif <- diff(c(max(head(sublp,5)),max(tail(sublp,5))))
      # if(oldlprdif > lprdif) lproughnesstarget <- oldlproughnesstarget
      # lprproposal = lproughnesstarget*2-oldlproughnesstarget
      # oldlproughnesstarget <- lproughnesstarget
      # if(max(lp) > max(tail(lp,30))) lprproposal <- min(.2,.5 * lprproposal)
      lproughnesstarget <- min(.4, max(.2, .01 + lproughnesstarget + .02 * (-1+2*rbinom(n = 1,size = 1,prob = .5))))
      # if(sd(tail(lp,30)) < lpnoisethresh) lproughnesstarget <- lproughnesstarget + .01
      
    }
    
    # if(i %% 200== 0){ #big changes in lp roughness
    #   lproughnesstarget <- sample(x = c(.1,.3,.6),1)
    #   oldlproughnesstarget <- lproughnesstarget
    # }
    
    # 
    # step[step > maxparchange] <- maxparchange
    # step[step < minparchange] <- minparchange
    
    if(nsubsets==1 && i > warmuplength && lp[i] < lp[i-1]) { #if worsening, update gradient faster
      # step[signdif]=step[signdif]*.5#lproughnesstarget
      # step = step * lproughnesstarget
      if(lp[i] < lp[i-10]) gmemory <- gmemory * .995
      # step=step*.5
      gsmooth[signdif]= gsmooth[signdif]*gmemory2^2 + (1-gmemory2^2) * g[signdif] #increase influence of gradient at inflections
      # gsmooth[signdif]= gsmooth[signdif]*.5 + .5 * g[signdif] #increase influence of gradient at inflections
    }
    
    if(plot && i %% as.numeric(plot) ==0 && (nsubsets==1 || i > nsubsets)){
      par(mfrow=c(2,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
      plot(pars,col=1:length(pars))
      points(changepars,pch=17,col='red')
      plot(log(abs(step*gsmooth)+1e-50),col=1:length(pars))
      if(nsubsets==1) plot(tail(log(-(lp-max(lp)-1)),100),type='l')
      if(nsubsets>1)  plot(tail(log(-(subsetlpstore-max(subsetlpstore)-1)),500),type='l',
        ylim=c(0,max(tail(log(-(subsetlpstore-max(subsetlpstore)-1)),500))))
      # plot(gamweights,col=1:length(pars))
      parrange=(apply(parstore,1,function(x) abs(diff(range(x)))))
      abline(h=(parrangetol))
      matplot(t(parstore[
        which(parrange > sort(parrange,decreasing = TRUE)[min(c(length(pars),5))]),,drop=FALSE]),
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
      }
    }
    
    #check convergence
    if(i==1) lpdiff1 <- NA
    if(nsubsets==1 && 
        i > nconvergeiter && 
        max(tail(lp,nconvergeiter)) ==max(lp)){
      
      if( (i - bestiter) > nconvergeiter*3 &&
          mean(sign(diff(tail(lp,nconvergeiter)))) < .3) converged <- TRUE #time since best
      
      lpdiff=(max(tail(lp,nconvergeiter)) - max(lp[1:(i-nconvergeiter)]))/nconvergeiter #avg progress over convergerange
      
      #convergence progress update and itertol check
      
      if(is.na(lpdiff1) || lpdiff1 < lpdiff) lpdiff1 <- lpdiff
      progressEst <- round(min(100,100*(1 - log(lpdiff/itertol) / log(lpdiff1/itertol))),2)
      
      message(paste0('\rProgress est. = ',progressEst,'%, LPchange = ',signif(lpdiff,digits = 2), 
        ', Iter = ',i, ', LP = ',lp[bestiter]),appendLF = FALSE)
      # ', grad = ', mean(sqrt((g^2))), ', gmem = ', gmemory,'  lprt = ',lproughnesstarget))
      
      if(lpdiff < itertol & lpdiff > 0) converged <- 1
      
      prevbest=max(head(lp,length(lp)-nconvergeiter))
      
      if(i > (nconvergeiter*2)) if(max(lp)-prevbest > 0 && max(lp)-prevbest < itertol) converged <- 5
      
      if(!is.na(parrangetol)){
        if(max(apply(parstore,1,function(x) diff(range(x)))) < parrangetol) converged <- 3
      }
      
    }
    
    if(nsubsets > 1 && (i %% (nsubsets))==0){ #store subset lps
      if(i >=(nsubsets*2)){ #then oldsubset exists
        oldsubsetlp <- subsetlp
        oldsubsetpars <- pars
      }
      subsetlp <- sum(tail(lp,nsubsets))
      subsetlpstore <- c(subsetlpstore,subsetlp)
      subsetpars <- apply(parstore,1,mean)
      
      if(subsetlp > bestsubsetlp) {
        bestsubsetlp <- subsetlp
        bestsubsetpars <- subsetpars
      }
    }
    
    
    
    if(nsubsets > 1 && #subsampling convergence check
        i > (nconvergeiter*nsubsets) && 
        (i %% (nsubsets))==0){ 
      
      lpdiff <- (subsetlp-max(subsetlpstore[1:((i/nsubsets)-nconvergeiter)]))/nconvergeiter
      
      if(i >=(nsubsets*2)){ #then oldsubset exists
        
        if(lpdiff >0){ #epochwise lproughness update!
          step = step*2 
        } else {
          step = step*.1 
          gsmooth=gsmooth*.5+g*.5
        }
        
        if(lpdiff >= 0){
        #convergence progress update and itertol check
        if(is.na(lpdiff1) || lpdiff1 < lpdiff) lpdiff1 <- lpdiff
        progressEst <- round(min(100,100*(1 - log(lpdiff/itertol) / log(lpdiff1/itertol))),2)
        lpchange=signif(lpdiff,digits = 2)
        }
        
        message(paste0('\rProgress est. = ',progressEst,'%, LPchange = ',lpchange, 
          ', Iter = ',i, ', LP = ',bestsubsetlp,', worseCount = ',worsecount,'           '),appendLF = FALSE)
        # ', grad = ', mean(sqrt((g^2))), ', gmem = ', gmemory,'  lprt = ',lproughnesstarget))
        
        if(lpdiff < itertol & lpdiff > 0 && bestsubsetlp==subsetlp) converged <- 1
        
        
        if(subsetlp < bestsubsetlp) worsecount <- worsecount + 1 else worsecount = 0
        
        if(worsecount > worsecountconverge) converged = 3
        # if(plot > 0) message(paste0(subsetlp, ' , bestlp = ',bestsubsetlp,' , worsecount = ',worsecount))
        # message(bestsubsetlp)
        # message(paste0(worsecount))
        
        # if(subsetlp <= bestsubsetlp &&
        #     all(abs(subsetpars-oldsubsetpars) < (itertol/1000))) converged <- 4
        # if(plot > 0 && i %%plot==0) plot(abs(subsetpars-oldsubsetpars))
        # if(subsetlp < (oldsubsetlp-lpnoisethresh)) step = step * .8 else step=step*1.1
        
        subsetorder = sample(1:nsubsets)
        if(converged) bestpars=bestsubsetpars
        
      }
    } #end subsampling check
    
  }
  message('')
  convergemessages <- c(
    paste0('Converged -- lp change within tol(',itertol,')'),
    'Converged -- parameter changes within parrangetol',
    'Converged -- count of non-improving iterations exceeded',
    'Converged -- lp and par change within tolerances',
    'Converged -- lp change within itertol')
  if(converged > 0) message(convergemessages[converged]) else message('Max iterations reached')
  out=list(itervalues = lp, value = max(lp),
    par=bestpars,parstore=parstore,gstore=gstore,lpstore=tail(lp,nstore))
  
  return(out)#,gstore=gstore,pstore=pstore) )
}

