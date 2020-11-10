ctLongToWideSF <- function(fit){
  dat <- data.frame(standatatolong(standata = fit$standata,
    ctm = fit$ctstanmodel,origstructure = TRUE))
  dat <- dat[,c(fit$ctstanmodelbase$subjectIDname,
    fit$ctstanmodelbase$timeName,
    fit$ctstanmodelbase$manifestNames,
    if(fit$ctstanmodelbase$n.TDpred > 0) fit$ctstanmodelbase$TDpredNames,
    if(fit$ctstanmodelbase$n.TIpred > 0) fit$ctstanmodelbase$TIpredNames)]
  dat=data.table(dat)
  dat[ ,WhichObs:=0:(.N-1),by=eval(fit$ctstanmodelbase$subjectIDname)]
  dat=melt(dat,id.vars = c(fit$ctstanmodelbase$subjectIDname,'WhichObs'))

  dat$WhichObs <- paste0('T',dat$WhichObs);
  dat=dcast(dat,'id~variable+WhichObs')
  lapply(fit$ctstanmodelbase$TIpredNames, function(x){
    colnames(dat)[grep(paste0('\\b',x,'\\_T0'),colnames(dat))] <<- x
  })
  dat=data.frame(dat)
  return(dat)
}

ctSaturatedFitConditional<-function(dat,ucols,reg,hmc=FALSE,covf=NA,verbose=0){
  o1 = ucols
  o2=(1:ncol(dat))[-o1]
  
  if(all(is.na(covf))) covf=covml(dat,reg = reg,hmc=hmc,verbose=verbose)
  if(length(o2) > 0){
    m=(covf$cp$mu)
    m1=m[o1]
    m2=m[o2]
    s=covf$cp$covm
    s12=s[o1,o2,drop=FALSE]
    isig2=MASS::ginv(s[o2,o2,drop=FALSE])
    sigma=s[o1,o1,drop=FALSE]-s[o1,o2,drop=FALSE] %*% isig2 %*% s[o2,o1,drop=FALSE]
    sigma=solve(solve(s)[o1,o1,drop=FALSE])
    llrow=sapply(1:nrow(dat),function(i){
      llr=NA
      a=as.numeric(c(dat[i,o2]))
      d2=!is.na(a)
      d1 = !is.na(dat[i,o1])
      if(any(is.na(a))) isig2=MASS::ginv(s[o2,o2,drop=FALSE][d2,d2,drop=FALSE])
      if(sum(d1)>0){
        mu=c(m1[d1])
        if(sum(d2)>0) mu=mu+s12[d1,d2,drop=FALSE] %*% isig2 %*% c(a[d2]-m2[d2])
        
        llr=mvtnorm::dmvnorm(x = dat[i,o1][d1],mean = mu,sigma = sigma[d1,d1,drop=FALSE],log=TRUE)
      }
      return(llr)
    })
    # covf$ll_unconditional <- covf$ll
    covf$ll <- sum(llrow,na.rm=TRUE)
  } else { #if not conditional
    
    covdat <- covdata(ndat = dat[,o1,drop=FALSE],reg=reg)
    smf <- rstan::sampling(stanmodels$cov,data=covdat,chains=0)
    cp <- rstan::constrain_pars(smf,covf$fit$par)
    
    #   covf$llrow <- sapply(1:nrow(dat),function(i){
    #   llr=NA
    #   d1 = !is.na(dat[i,o1])
    #   if(sum(d1)>0){
    #     mu=c(covf$cp$mu)[d1]
    #     llr=try(mvtnorm::dmvnorm(x = as.numeric(dat[i,o1[d1]]),mean = mu,sigma = covf$cp$covm[d1,d1,drop=FALSE],log=TRUE))
    #     if('try-error' %in% class(llr)) 
    #   }
    #   return(llr)
    # })
    covf$cp <- cp
    covf$ll <- sum(cp$llrow,na.rm=TRUE)
    # print( covf$ll)
    # 
  }
  
  return(covf)
}

ctSaturatedFit <- function(fit,conditional=FALSE,reg=0, hmc=FALSE,
  time=FALSE, oos=TRUE, folds=10, cores=2,verbose=0){
  dat <-ctLongToWideSF(fit)
  dat=dat[,-1]
  if(!time) dat <- dat[,-grep(paste0('\\b',fit$ctstanmodelbase$timeName,'_T'),colnames(dat))]
  
  if(!conditional){
    dat <- data.frame(dat)[,unique(c(sapply(fit$ctstanmodelbase$manifestNames,function(x){
      grep(paste0('\\b',x,'\\_T'),colnames(dat))
    }))),drop=FALSE]
    dat <- dat[,apply(dat,2,function(x) any(!is.na(x))),drop=FALSE] #drop na columns
    ucols=1:ncol(dat)
  } else { #if conditional
    #get unconditional columns
    dat <- dat[,apply(dat,2,function(x) any(!is.na(x))),drop=FALSE] #drop na columns
    ucols=unique(c(unlist(sapply(fit$ctstanmodelbase$manifestNames,function(x){
      grep(paste0('\\b',x,'\\_T'),colnames(dat))
    }))))
  }
  
  # message('Min obs per columns: ',paste(apply(dat,2,function(x) nrow(dat)-sum(is.na(x))),collapse=', '))
  
  if(cores > 1){
    cl=parallel::makeCluster(cores,'PSOCK')
    parallel::clusterExport(cl,c('ucols','dat','reg'),environment())
    on.exit({parallel::stopCluster(cl)},add = TRUE)
  }
  
  bootstrap=FALSE
  
  counter=0
  minobs=0
  while(counter <10 && minobs < 3){
    srows <- sample(1:nrow(dat),
      nrow(dat),# *ifelse(bootstrap,folds,1),
      replace=bootstrap)
    
    srows <- split(
      srows,
      sort(1:length(srows) %% folds))
    
    minobspercol <- unlist(apply(dat,2,function(d) min(unlist(lapply(srows,function(sr) ( sum(!is.na(d[-sr]))))))))
    minobs<-min(minobspercol)
    if(minobs <5)  message('Minimum number of observations in a column = ', minobs, ', sampling again...')
    counter = counter + 1
  }
  if(counter==10) stop('Too few observations per fold -- try more folds?')
  message('Min obs per col: ', paste0(minobspercol,collapse=', '))
  
  covf = ctSaturatedFitConditional(dat = dat,ucols = ucols,reg = reg,verbose=0)
  covfindep = covml(ndat = dat[,ucols,drop=FALSE],reg = reg,
    hmc=hmc,independent=TRUE)
  covf$llindependent <- covfindep$ll
  
  #independence fit of residuals
  err=ctStanKalman(fit = fit,collapsefunc = mean)$errprior
  err=data.table(data.frame(id=fit$standata$subject,matrix(err,ncol=dim(err)[3])))
  err[ ,WhichObs:=(1:.N),by=id]
  err = data.table::dcast(data = data.table(melt(err,id.vars = c('id','WhichObs'))),
    formula='id~WhichObs+variable')
  err=data.frame(err)
  err <- err[,apply(err,2,function(x) any(!is.na(x))),drop=FALSE][,-1,drop=FALSE] #remove na and id cols
  covfindepresid = covml(ndat = err,reg = reg,
    hmc=hmc,independent=TRUE)
  covf$llresid <- covfindepresid$ll
  
  
  
  if(oos){
    sf = flexlapply(cl,srows,function(x){
      library(ctsem)
      datheldout <- dat
      datheldout[x,ucols] <- NA
      #saturated model oos
      fin=ctSaturatedFitConditional(dat = datheldout,ucols = ucols,reg = reg,hmc=hmc)
      fitoos=ctSaturatedFitConditional(dat = dat[x,,drop=FALSE],ucols = ucols,reg = reg,covf=fin)
      fin$lloos <- fitoos$cp$llrow
      plot(fitoos$cp$llrow)
      
      #independence model oos
      findep=covml(ndat = datheldout[,ucols,drop=FALSE],reg = reg,
        hmc=hmc,independent=TRUE)
      covdat <- covdata(ndat = dat[x,ucols,drop=FALSE],reg=reg,independent = TRUE)
      smf <- rstan::sampling(stanmodels$cov,data=covdat,chains=0)
      cp <- rstan::constrain_pars(smf,findep$fit$par)
      fin$llindependentoos <- cp$llrow
      
      return(fin)
    },cores = cores)
    
    #out of sample saturated ll
    llfold=unlist(lapply(sf,function(x) x$lloos))
    llfold=llfold[match(1:nrow(dat),(unlist(srows)))]
    # llfold=llfold[unlist(srows)]#sapply(unlist(srows),function(x) which(unlist(srows) %in% x))]
    
    
    plot(llfold,covf$cp$llrow, main='OOS vs IS LogLik',ylab='In sample LL',
      xlab='Out of sample LL')
    abline(a=0,b=1)
    
    covf$lloos = sum(llfold)
    covf$llrowoos = llfold
    # 
    
    #out of sample independent ll
    covf$llrowindependentoos <-  unlist(lapply(sf,function(x) x$llindependentoos))[match(1:nrow(dat),(unlist(srows)))]
    plot(covf$llrowindependentoos,llfold, main='LL OOS indep. vs OOS sat.',ylab='Out of sample saturated LL',
      xlab='Out of sample independence LL')
    abline(a=0,b=1)
    covf$llindependentoos <- sum(unlist(sapply(sf,function(x) x$llindependentoos)),na.rm=TRUE)
  }
  
  d=nrow(covf$cp$covm)
  covf$npars=length(covf$cp$mu)+(d^2-d)/2-length(fit$stanfit$rawest)
  fit$stanfit$saturated <- covf
  
  covf$test <- 1-pchisq(
    covf$ll-fit$stanfit$transformedparsfull$ll,
    df=covf$npars)
  covf$aicdiff = 
    (2*length(fit$stanfit$rawest) - 2*fit$stanfit$transformedparsfull$ll) -
    (2 * covf$npars - 2* covf$ll)
  message('Model ll: ',fit$stanfit$transformedparsfull$ll,'; Saturated ll: ',covf$ll,
    '; OOS ll = ',covf$lloos,'; P = ',
    covf$test,'; AIC diff = ',covf$aicdiff)
  return(fit)
}

ctDataMelt <- function(dat,id='id',by='time', combinevars=NULL){
  # setnames(dat,id,'Subject')
  # by <- by[-which(by %in% 'WhichObs')]
  if(!is.null(combinevars)) {
    # dat[ ,WhichObstmp:=(1:.N),by=eval(id)]
    dat <- ctDataCombineSplit(dat = dat,idvars = c(id,by),vars = combinevars)
    # dat <- dcast(dat,formula = 'id +WhichObstmp~ variable',value.var = 'value')
    # dat$WhichObstmp <- NULL
  }
  if(is.null(combinevars)) dat <- melt(dat,id.vars = c(id,by))
  return(dat)
}

#takes list of vars, combines according to structure
ctDataCombineSplit <- function(dat, idvars, vars){ #no splits yet
  ldat=as.data.table(dat)
  
  idvarcombine <- idvars[which(idvars %in% names(vars) & !idvars %in% unlist(vars))]
  if(length(idvarcombine) > 0){
    for(vi in 1:length(idvarcombine)){
      # ldat[,eval(idvarcombine[vi]):= mean(eval(vars[[idvarcombine[vi]]]))]
      ldat[[eval(idvarcombine[vi])]] <- apply(data.frame(ldat)[,vars[[idvarcombine[vi]]]],1,mean,na.rm=TRUE)
    }
  }
  
  ldat <- melt(data = ldat,id.vars = idvars)
  
  #remove unneeded variables
  # vars <- fit$ctstanmodel$manifestNames
  ldat <- ldat[variable %in% eval(unlist(vars))]
  
  #combine variables
  # combinevars=list(ithr=c('T1ithr','T2ithr','T3ithr','T4ithr'),
  #   itfa=c('T1itfa','T2itfa','T3itfa','T4itfa'),
  #   ashr=c('T1ashr','T2ashr','T3ashr','T4ashr'),
  #   asfa=c('T1asfa','T2asfa','T3asfa','T4asfa') )
  
  ldat$variable <- as.character(ldat$variable)
  ldat <- as.data.frame(ldat)
  for(vi in 1:length(vars)){
    ldat[ldat$variable %in% vars[[vi]],'variable'] <- names(vars)[vi]
  }
  ldat <- as.data.table(ldat)
  
}


meltcov <- function(covm){
  if(is.null(dimnames(covm))) dimnames(covm) = list(paste0('v',1:nrow(covm)),paste0('v',1:nrow(covm)))
  covm=data.frame(row=factor(rownames(covm),levels=rownames(covm)),covm)
  
  # # browser()
  # div = 4:6
  # divisor = which(nrow(covm)%%div == min(nrow(covm)%%div))+3 #find minimum remainder for splitting factor
  # 
  # if(is.na(factors)) factors <- suppressWarnings(factor(c(matrix(paste0('g',1:(ceiling(nrow(covm)/divisor))),nrow=nrow(covm)))))
  # 
  # covm=cbind(covm,Group=factors)
  o=data.table::melt(data.table(covm),id.vars=c('row'))
  colnames(o)[colnames(o) %in% c('row','variable')] <- c('Var1','Var2')
  return(o)
}

ctStanFitMelt <- function(fit, by,combinevars=NULL,maxsamples='all'){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
  
  
  datasources <- c('Data','StatePred','Residuals')
  if(!is.null(fit$generated)) datasources <- c(datasources,'PostPred')
  if(!is.null(fit$priorpred)) datasources <- c(datasources,'PriorPred')
  
  datbase <- data.table(standatatolong(standata = fit$standata,origstructure = TRUE,ctm = fit$ctstanmodel))
  datbase[ ,WhichObs:=(1:.N),by=eval(fit$ctstanmodelbase$subjectIDname)]
  datbase <- data.frame(datbase)
  
  dat <- matrix(NA,nrow=0,ncol=0) #blank placeholder
  
  for(dsi in datasources){
    dexists <- FALSE
    if(dsi == 'Data'){
      dexists<-TRUE
      d<-list()
      d$Y<- array(t(t(datbase[,(fit$ctstanmodelbase$manifestNames)])), dim=c(1,dim(fit$standata$Y)))
      d$llrow <- fit$stanfit$transformedparsfull$llrow
    }
    
    if(dsi=='PriorPred'){
      dexists<-TRUE
      d <- fit$priorpred
    }
    if(dsi=='PostPred'){
      dexists<-TRUE
      d <- fit$generate
    }
    
    if(dsi== 'StatePred'){ #use kalman predictions
      dexists<-TRUE
      d<-list()
      d$Y<- array(fit$stanfit$transformedparsfull$ya[1,1,,],dim=dim(fit$stanfit$transformedparsfull$ya[1,1,,,drop=FALSE])[-1])
      d$Y[rep(fit$standata$Y,each=dim(d$Y)[1])==99999] <- NA
      d$llrow <- d$llrow <- fit$stanfit$transformedparsfull$llrow # array(NA,dim=c(1,dim(d$Y)[2])) 
    }
    
    if(dsi=='Residuals'){
      dexists<-TRUE
      d<-list()
      d$Y<- d$Y<- array(t(t(datbase[,(fit$ctstanmodelbase$manifestNames)])), dim=c(1,dim(fit$standata$Y))) -
        array(fit$stanfit$transformedparsfull$ya[1,1,,],dim=dim(fit$stanfit$transformedparsfull$ya[1,1,,,drop=FALSE])[-1])
      d$llrow <- d$llrow <- fit$stanfit$transformedparsfull$llrow #array(NA,dim=c(1,dim(d$Y)[2]))
    }
    
    if(dexists){
      if(maxsamples=='all') samples=1:dim(d$Y)[1] else samples = sample(1:dim(d$Y)[1],min(maxsamples,dim(d$Y)[1]))
      gdat <- data.table(datbase[rep(seq_len(nrow(datbase)), length(samples)),])
      
      gdat[, (fit$ctstanmodelbase$manifestNames)]   <- #confusing aperm needed here, wish I understood why...
        data.table(matrix(aperm(d$Y[samples,,,drop=FALSE],c(2,1,3)), ncol=dim(d$Y)[3]))
      gdat$Sample=rep(samples,each=dim(d$Y)[2])
      gdat$LogLik<-c(d$llrow[samples,])
      gdat$DataSource <- dsi
      
      if(nrow(dat)==0) dat <- gdat else dat <- rbind(dat, gdat)
    }
  }
  
  by=c(by,'Sample','WhichObs','DataSource')
  if(is.na(combinevars[1])){
    combinevars = setNames(
    colnames(datbase)[!colnames(datbase) %in% by],colnames(datbase)[!colnames(datbase) %in% by])
  combinevars<-c(combinevars,LogLik='LogLik')
  }
  
  dat <- ctDataMelt(dat=dat,id=fit$ctstanmodelbase$subjectIDname, by=by,combinevars = combinevars)
  dat$Sample <- factor(dat$Sample)
  dat$DataSource <- factor(dat$DataSource)
  
  return(dat)
}


ctCheckFit2 <- function(fit, 
  data=TRUE, postpred=TRUE, priorpred=FALSE, statepred=FALSE, residuals=FALSE,
  by=fit$ctstanmodelbase$timeName,
  TIpredNames=fit$ctstanmodelbase$TIpredNames,
  nsamples=10, covplot=FALSE, corr=TRUE, combinevars=NA, fastcov=FALSE,
  aggfunc=mean,aggregate=TRUE,
  groupby='split', byNA=TRUE,lag=0,
  smooth=TRUE, k=10,breaks=4,entropy=FALSE,reg=FALSE,verbose=0){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
  
  covORcor <- function(m){
    if(corr) return(cov2cor(m)) else return(m)
  }
  
  # combinevars<-c(combinevars,LogLik='LogLik')
  
  dat <- ctStanFitMelt(fit = fit,combinevars = combinevars, by=by,maxsamples = nsamples)
  
  if(is.na(combinevars[1])) combinevars <- setNames(unique(dat$variable), unique(dat$variable))
  
  if(!data) dat<-dat[!DataSource %in% 'Data']
  if(!priorpred) dat<-dat[!DataSource %in% 'PriorPred']
  if(!postpred) dat<-dat[!DataSource %in% 'PostPred']
  if(!statepred) dat<-dat[!DataSource %in% 'StatePred']
  if(!residuals) dat<-dat[!DataSource %in% 'Residuals']
  
  # if(!'WhichObs' %in% by) mdat$WhichObs <- NULL
  
  dat = dat[!is.na(value),]
  
  
  if(lag!=0){ 
    
    wdat=dcast(dat,paste0(
      fit$ctstanmodelbase$subjectIDname,
      '+Sample+DataSource+WhichObs+time','~variable'))
    
    lagcols=colnames(wdat)[!colnames(wdat) %in% c('id','Sample','DataSource','WhichObs')]
    
    wdat[, (paste0("lag",  rep(lag, times = length(lagcols)),'_',rep(lagcols, each = length(lag)))) :=  
        shift(.SD, lag), by=c('id','Sample','DataSource')]
    
    dat <- melt(wdat,id.vars=c('id','Sample','DataSource','WhichObs','time'))
    
  }
  
  breaks = min(breaks,length(unique(dat[[by]][!is.na(dat[[by]])])))
  # k = min(k,length(unique(dat[[by]][!is.na(dat[[by]])]))-1)
  
  if(covplot ||entropy){
    discdat=dat #make copy for use in covs
    discdat[variable %in% TIpredNames & WhichObs > 1] <- NA #remove later tipreds to avoid duplicate cov columns

    nontivars <- unique(discdat$variable)[!unique(discdat$variable) %in% TIpredNames]
    
    if(is.double(dat[[by]])){
      if(requireNamespace('arules')){
        discdat[[by]] <- arules::discretize(dat[[by]], #discretize
          method='cluster',breaks = breaks,labels=FALSE)
      } else stop('arules package needed for discretization!')
    }
    if(covplot){
      corlist <- list()
      datasources <- as.character(unique(dat$DataSource))
      for(dsi in unique(dat$DataSource)){ #make wide discretized data and get cov
        
        wdat <- dcast(discdat[DataSource==dsi],
          paste0(
            fit$ctstanmodelbase$subjectIDname,
            '+Sample',if(!aggregate) '+WhichObs','~variable + ',by),fun.aggregate=aggfunc,na.rm=TRUE)
  

        
        #put loglik on one side
        if('LogLik_1' %in% colnames(wdat)) wdat <- cbind(wdat[,colnames(wdat) %in% paste0('LogLik_',1:breaks),with=FALSE],
          wdat[,!colnames(wdat) %in% paste0('LogLik_',1:breaks),with=FALSE])
        
        if(byNA) breakset <-c(1:breaks,NA) else breakset=1:breaks #some values can't be placed in a by column if by = NA
          
        if(groupby == 'split'){
          colorder <- unlist(lapply(breakset, function(b) grep(paste0('\\_',b,'$'),x=colnames(wdat))))
          wdat <- wdat[,colorder,with=FALSE]
        } else if(!byNA) wdat <- wdat[,-grep('\\_NA%',colnames(wdat)),with=FALSE] #if not seperating by group, still remove NA by relations if needed
        
        sapply(TIpredNames,function(x){ #fix tipred column names
          if(paste0(x,'_1') %in% colnames(wdat)){
            colnames(wdat)[ colnames(wdat) %in% paste0(x,'_1')] <<- x
          }
        })
        
        wdat <- cbind(wdat[,c(which(colnames(wdat) %in% TIpredNames), which(!colnames(wdat) %in% TIpredNames)),with=FALSE])
        
        if(!fastcov) corlist[[dsi]] <- covml(
          data.frame(wdat)[,!colnames(wdat) %in% c('Sample',fit$ctstanmodelbase$subjectIDname),drop=FALSE],
          reg=reg,
          verbose=verbose)$cp$covm
        
        if(fastcov) corlist[[dsi]] <- cov(
          data.frame(wdat)[,!colnames(wdat) %in% c('Sample',fit$ctstanmodelbase$subjectIDname),drop=FALSE],
          use='pairwise.complete.obs')
        
        # if(covsplitdiffs){
        #   browser()
        #   splitcovs <- plyr::laply(breakset[!is.na(breakset)], function(b){
        #     corlist[[dsi]][,colnames(corlist[[dsi]]) %in% c(TIpredNames, paste0(nontivars,'_',b)), drop=FALSE]
        #   })
        #   commonnames <- x
        #   splitcovs <- x
        #   cplot=corplotmelt(meltcov(corlist[[datasources[i]]]),
        #     label=paste0(' corr.'),limits=c(-1,1))
        #   cplot = cplot + ggplot2::labs(title=paste0(datasources[i],' correlations, split by ', by))
        #   print(cplot)
        # }
          
      } #finish datasource loop
      # browser()
      if(corr) covplotlims <- c(-1,1) else covplotlims <- NA
      # browser()
      for(i in seq_along(corlist)){ #regular corplots
        cplot=corplotmelt(meltcov(corlist[[datasources[i]]]),
          label=paste0(ifelse(corr,'Corr.','Cov.')),limits=covplotlims)
        cplot = cplot + ggplot2::labs(title=paste0(datasources[i],
          ifelse(corr,' correlations',' covariances'),', split by ', by))
        print(cplot)
      }
      for(i in seq_along(corlist)){
        for(j in seq_along(corlist)){ #difference corplots
          if(j > i){ #only plot unique differences
            cplot=corplotmelt(meltcov(corlist[[datasources[i]]] - corlist[[datasources[j]]]),
              label=paste0(ifelse(corr,'Corr.','Cov.'),' diff.'),limits=covplotlims)
            cplot = cplot + ggplot2::labs(title=paste0(datasources[i],' - ',datasources[j] ,' correlations, split by ', by))
            print(cplot)
          }
        }
      }
    }
    if(entropy){
      # 
      dentropy <- -mean(covdat$cp$llrow)#lp/nrow(wdat)
      
      gdentropy <- sapply(unique(wgdat$Sample), function(x){
        print(x)
        -mean(covml(data.frame(wgdat)[wgdat$Sample==x,!colnames(wgdat) %in% 
            c('Sample',fit$ctstanmodelbase$subjectIDname)],verbose=verbose)$cp$llrow)
      }) 
      print(sum(covdat$cp$llrow))
      message('Entropy: Fit = ', -sum(fit$stanfit$transformedparsfull$ll)/nrow(wdat),'; Data = ',dentropy,'; Mean gen.data = ',mean(gdentropy),'  SD gen. data = ',sd(gdentropy))
    }
  }
  
  if(!covplot){ #then plot bivariate relations
    vars <- names(combinevars)
    dat$manifest <- dat$variable %in% vars
    
    # dat$DataSource <- factor(as.character(dat$DataSource), levels = c( "Data", "StatePred","PostPred","PriorPred"))
    # browser()
    
    g=ggplot(data = dat[manifest==TRUE],
      mapping = aes_string(x=by,y='value',
        colour='DataSource'
        ,fill='DataSource'
        # , group='Sample'
        # ,alpha='DataSource'
      )
    ) +theme_bw() +
      # scale_color_brewer(palette = 'Set1') +
      # scale_fill_brewer(palette = 'Set1')
      # +
      scale_colour_manual(values = c("blue", "red","green",  "Orange","Black"))+
      scale_fill_manual(values = c("blue", "red","green",  "Orange",'Black'))#+
    # scale_alpha_manual(values = c(1, max(.05,(.3/nsamples)),max(.05,(.3/nsamples)),1))
    
    if(smooth) {
      g = g + geom_smooth(data=dat[manifest==TRUE & !DataSource %in% c('Data','StatePred','Residuals')],
        aes(group=Sample),
        alpha= max(.05,(.3/nsamples)),
        linetype=0,#3,#ifelse(nsamples==1,1,0),
        stat="smooth",#se = FALSE,#nsamples==1,
        size=.1
        ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
      
      g = g + geom_smooth(data=dat[manifest==TRUE & DataSource %in% c('Data','StatePred','Residuals')],
        # aes(colour=DataSource,alpha=NULL),
        stat="smooth",se = TRUE,size=1,alpha=.3
        ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
    }
    if(!smooth) {
      # if(nsamples > 1) g = g + stat_summary(fun=mean,geom = "line",size=1) #geom_line(stat=mean)
      
      # if(nsamples ==1)
      g = g +  stat_summary(data=dat[manifest==TRUE & !DataSource %in% c('Data','StatePred','Residuals')],
        aes(group=Sample),
        fun.data = function(x) list(y=mean(x),
          ymin=mean(x,na.rm=TRUE)-sd(x)/sqrt(length(x)), 
          ymax=mean(x)+sd(x)/sqrt(length(x))),
        geom = "ribbon",
        alpha= max(.05,max(.05,sqrt(.2/nsamples))),
        linetype=ifelse(nsamples==1,1,0))
      
      g = g +  stat_summary(data=dat[manifest==TRUE & DataSource %in% c('Data','StatePred','Residuals')],
        fun.data = function(x) list(y=mean(x),
          ymin=mean(x,na.rm=TRUE)-sd(x)/sqrt(length(x)), 
          ymax=mean(x)+sd(x)/sqrt(length(x))),
        geom = "ribbon",alpha=.3)
    }
    g=g+facet_wrap(facets = vars(variable),scales = 'free')
    print(g)
  }
}


#' Check absolute fit of ctFit or ctStanFit object.
#'
#' @param fit ctsem fit object.
#' @param niter number of data generation iterations to use to calculate quantiles.
#' @param probs 3 digit vector of quantiles to return and to test significance.
#'
#' @return List containing a means and cov object, computed by sorting data into discrete time points.
#' cov is a numeric matrix containing measures of the covariance matrices for observed and simulated data. 
#' The MisspecRatio column shows Z score difference for each lower triangular index of the covariance matrix of data --
#' observed covariance minus mean of generated, weighted by sd of generated covariance.
#' means contains the empirical and generated data means.
#' @export
#' @importFrom data.table dcast
#' 
#' @details for plotting help see \code{\link{plot.ctsemFitMeasure}}
#'
#' @examples
#' \donttest{
#' scheck <- ctCheckFit(ctstantestfit,niter=50)
#' }
ctCheckFit <- function(fit, niter=500,probs=c(.025,.5,.975)){
  id=NULL #global warnings
  if(!class(fit) %in% c('ctStanFit','ctsemFit')) stop('not a ctsemFit or ctStanFit object!')
  
  
  if(class(fit)=='ctsemFit'){
    manifestNames = fit$ctmodelobj$manifestNames
    nmanifest=fit$ctmodelobj$n.manifest
    maxtp=fit$ctmodelobj$Tpoints
    if('Kalman' %in% fit$ctfitargs$objective) {
      suppressMessages(wdat <- ctLongToWide(fit$mxobj@data$observed,id='id',time='time',
        manifestNames = manifestNames)[,paste0(manifestNames,'_T',1),drop=FALSE])
    } else  wdat <- fit$mxobj@data$observed[,paste0(rep(manifestNames,each=fit$ctmodelobj$Tpoints),'_T',
      0:(maxtp-1)),drop=FALSE]
  }
  
  if(class(fit)=='ctStanFit') {
    if(fit$data$nsubjects==1) stop('Only for nsubjects > 1!')
    # 
    manifestNames=fit$ctstanmodel$manifestNames
    nmanifest=fit$ctstanmodel$n.manifest
    ldat <- cbind(fit$data$subject,fit$data$time,fit$data$Y)
    tpoints <- max(unlist(lapply(unique(fit$data$subject),function(x) length(fit$data$subject[fit$data$subject==x]))))
    colnames(ldat) <- c('id','time', manifestNames)
    dt = cbind(data.table(id=fit$data$subject),data.table(fit$data$Y))[ ,.(discrete.time.point=1:.N),by=id]
    discrete.time.point=NULL #global variable complaint
    maxtp=max(dt[,discrete.time.point])
    suppressMessages(wdat <- ctLongToWide(ldat,id='id',time='time',
      manifestNames = manifestNames)[,paste0(manifestNames,'_T',
        rep(0:(tpoints-1),each=nmanifest)),drop=FALSE][,paste0(rep(manifestNames,each=maxtp),'_T',
          0:(maxtp-1))])
  }
  
  ecov <- cov(wdat,use = "pairwise.complete.obs")
  emeans <- (matrix(apply(wdat,2,mean,na.rm=TRUE),ncol=nmanifest))
  colnames(emeans) = manifestNames
  rownames(emeans) = paste0('T',0:(maxtp-1))
  
  covarray<-array(NA,dim = c(dim(ecov),niter))
  means <- array(NA,dim=c(maxtp,nmanifest,niter))
  
  if(class(fit)=='ctStanFit'){
    if(is.null(fit$generated) || dim(fit$generated$Y)[2] < niter){
      ygen <- ctStanGenerateFromFit(fit,fullposterior=FALSE,nsamples=niter)$generated$Y #array(e$Ygen,dim=c(ygendim[1] * ygendim[2],ygendim[-1:-2]))
    } else ygen <- fit$generated$Y
    wide <- matrix(NA, nrow=length(unique(fit$data$subject)),ncol=dim(ygen)[3]*maxtp)
    itervec <- sample(1:dim(ygen)[1],niter)
    
    for(i in 1:niter){
      idat <- data.table(ygen[i,,,drop=TRUE])
      colnames(idat) = manifestNames
      w <- dcast(data = cbind(dt,idat),
        formula= id ~ discrete.time.point,value.var=manifestNames)
      w=w[,-1]
      covarray[,,i] <- cov(w, use='pairwise.complete.obs')
      means[,,i] <- t(matrix(apply(w,2,mean,na.rm=TRUE),byrow=TRUE,ncol=nmanifest))
    }
    
  }
  
  if(class(fit)=='ctsemFit'){
    stop('OpenMx based fit objects not supported -- try ctModel types standt or stanct!')
  }
  # 
  covql <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[1],na.rm=TRUE)
  covqm <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[2],na.rm=TRUE)
  covqh <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[3],na.rm=TRUE)
  covmean <- ctCollapse(covarray,collapsemargin = 3,mean,na.rm=TRUE)
  covsd <- ctCollapse(covarray,collapsemargin = 3,sd,na.rm=TRUE)
  
  test<-matrix(NA,ncol=8,nrow=(nrow(covql)^2+nrow(covql))/2)
  counter=0
  rowname <- c()
  colname <- c()
  
  for(i in 1:nrow(covql)){
    for(j in 1:nrow(covql)){
      if(j <=i){
        counter=counter+1
        rowname <- c(rowname, rownames(ecov)[i])
        colname <- c(colname, colnames(ecov)[j])
        test[counter,] <- c(i,j,covmean[i,j],covql[i,j],covqm[i,j],covqh[i,j],ecov[i,j],
          ifelse((ecov[i,j] > covqh[i,j] || ecov[i,j] < covql[i,j]), TRUE,FALSE))
      }}}
  
  colnames(test) <- c('row','col','mean',paste0(probs*100,'%'), 'observed', 'significant')
  MisspecRatio <- (test[,'observed'] - test[,'mean']) / covsd[lower.tri(diag(nrow(covql)),diag = TRUE)] #((test[,'97.5%'] - test[,'2.5%']))^2
  sd <- covsd[lower.tri(diag(nrow(covql)),diag = TRUE)]
  test<- cbind(rowname,colname,as.data.frame(cbind(test,sd)),MisspecRatio)
  check <- list(cov=test,means=list(empirical=emeans, simulated=means))
  class(check) <- c('ctsemFitMeasure',class(check))
  return(check)
}

#' Misspecification plot using ctCheckFit output
#'
#' @param x Object output from ctCheckFit function.
#' @param indices Either 'all' or a vector of integers denoting which observations to 
#' include (from 1 to n.manifest * maximum number of obs for a subject, blocked by manifest).
#' @param covtype Column name of \code{$cov} sub object
#' @param cov Logical -- plot simulated cov vs observed?
#' @param means Logical -- plot simulated means vs observed?
#' @param cov2cor Logical -- convert covariances to correlations?
#' @param separatemeans Logical -- means from different variables on same or different plots?
#' @param ggcorrArgs List of arguments to GGally::ggcorr .
#' @param wait Logical -- wait for input before new plot?
#' @param ... not used.
#'
#' @return Nothing, just plots.
#' @export
#' @method plot ctsemFitMeasure
#'
#' @examples
#' if(w32chk()){
#' 
#'
#' scheck <- ctCheckFit(ctstantestfit,niter=50)
#' plot(scheck,wait=FALSE)
#' 
#' }
plot.ctsemFitMeasure <- function(x,indices='all', means=TRUE,separatemeans=TRUE, 
  cov=TRUE,covtype='MisspecRatio',cov2cor=FALSE,wait=TRUE,
  ggcorrArgs=list(data=NULL, cor_matrix =  get(covtype),
    limits=limits, geom = 'circle',max_size = 10,name=covtype),...){
  
  if(!covtype %in% colnames(x$cov)) stop('covtype not in column names of x!')
  
  if(means){
    if(separatemeans) manifests <- 1:dim(x$means$empirical)[2] else manifests <- 'all'
    for(mani in manifests){
      if(mani[1] > 1 && wait) readline('Press enter to continue')
      if(mani == 'all') mani <- 1:dim(x$means$empirical)[2]
      simd <- matrix(x$means$simulated[,mani,], nrow=dim(x$means$simulated)[1])
      if(length(mani)>1) vpar=mani else vpar=1
      matplot(simd,type='l',xlab='Time point',ylab='Mean',
        col=adjustcolor(rep(vpar, dim(x$means$simulated)[3]),alpha.f = .1),
        lty = vpar)
      matplot(x$means$empirical[,mani,drop=FALSE],type='l',col=vpar,lwd=2,add=TRUE,lty=vpar)
      legend('topright',colnames(x$means$empirical)[mani], col=vpar,text.col=vpar,bty='n',lty=vpar)
    }
  }
  
  if(cov){
    if(requireNamespace('GGally')){ 
      if(means && wait) readline('Press enter to continue')
      n=x$cov$rowname[match(x = unique(1:max(x$cov$row)),x$cov$row)]
      for(xi in covtype){
        mat <- matrix(NA,max(x$cov[,'row']),max(x$cov[,'row']))
        mat[upper.tri(mat,diag = TRUE)] = x$cov[,'MisspecRatio']
        mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
        dimnames(mat) <- dimnames(mat) <- list(n,n)
        if(indices[1]=='all') indices <-1:nrow(mat)
        if(cov2cor) mat <- cov2cor(mat)
        assign(x = xi, mat[indices,indices,drop=FALSE])
      }
      
      # main <- '(Observed - implied) / sd(implied)'
      if(cov2cor) limits <-c(-1,1) else limits <- range(get(covtype))
      # 
      
      corm <- meltcov(ggcorrArgs$cor_matrix)
      
      corplotmelt(corm,'MissSpec. Ratio',limits)
      
      
      # do.call(GGally::ggcorr,ggcorrArgs) #(data=NULL,cor_matrix =  get(covtype),limits=limits, geom = 'circle',max_size = 13,name=covtype,...)
    }
  }
  
}

corplotmelt <- function(corm, label,limits=NA){
  # browser()
  if(is.na(limits[1])) limits <- range(corm[,'value'],na.rm=TRUE)
  ggplot(data=(corm),aes_string(x='Var1',y='Var2',fill=('value')))+ #ifelse(groups,NULL,'Group')))+
    geom_tile( width=1,height=1,colour='black')+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
      midpoint = 0, limits = limits, space = "Lab", 
      name=label)  +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90))
}

