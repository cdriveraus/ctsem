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
  if(!is.null(combinevars)) {
    dat <- ctDataCombineSplit(dat = dat,idvars = c(id,by),vars = combinevars)
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
  if(is.null(dimnames(covm))) dimnames(covm) = list(paste0('v',1:nrow(covm)),
    paste0('v',1:ncol(covm)))
  
  colnames(covm)[colnames(covm) %in% ''] <- 
    paste0('v',(1:ncol(covm))[colnames(covm) %in% ''])
  rownames(covm)[rownames(covm) %in% ''] <- 
    paste0('v',(1:nrow(covm))[rownames(covm) %in% ''])
  covm=data.frame(row=factor(rownames(covm),levels=rownames(covm)),covm)
  
  o=data.table::melt(data.table(covm),id.vars=c('row'))
  colnames(o)[colnames(o) %in% c('row','variable')] <- c('Var1','Var2')
  return(o)
}

ctStanFitMelt <- function(fit, maxsamples='all'){
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
      d$Y<- array(t(t(datbase[,(fit$ctstanmodelbase$manifestNames)])), dim=c(1,dim(fit$standata$Y))) -
        array(fit$stanfit$transformedparsfull$ya[1,1,,],dim=dim(fit$stanfit$transformedparsfull$ya[1,1,,,drop=FALSE])[-1])
      d$llrow <- fit$stanfit$transformedparsfull$llrow #array(NA,dim=c(1,dim(d$Y)[2]))
    }
    
    if(dexists){
      if(maxsamples=='all') samples=1:dim(d$Y)[1] else samples = sample(1:dim(d$Y)[1],min(maxsamples,dim(d$Y)[1]))
      gdat <- data.table(datbase[rep(seq_len(nrow(datbase)), length(samples)),])
      
      gdat[, (fit$ctstanmodelbase$manifestNames)]   <- #confusing aperm needed here, wish I understood why...
        data.table(matrix(aperm(d$Y[samples,,,drop=FALSE],c(2,1,3)), ncol=dim(d$Y)[3]))
      gdat$Sample=rep(samples,each=dim(d$Y)[2])
      gdat$LogLik<-c(d$llrow[samples,])
      gdat[[fit$ctstanmodelbase$timeName]]<-rep(datbase[[fit$ctstanmodelbase$timeName]],times=length(samples))
      gdat$DataSource <- dsi
      
      if(nrow(dat)==0) dat <- gdat else dat <- rbind(dat, gdat)
    }
  }
  
  return(dat)
}


#' Visual model fit diagnostics for ctsem fit objects.
#'
#' @param fit ctStanFit object.
#' @param data Include empirical data in plots?
#' @param postpred Include post predictive (conditional on estimated parameters and covariates) distribution data in plots? 
#' @param priorpred Include prior predictive (conditional on priors) distribution data in plots? 
#' @param statepred Include one step ahead (conditional on estimated parameters, covariates, and earlier data points) distribution data in plots? 
#' @param residuals Include one step ahead error (conditional on estimated parameters, covariates, and earlier data points) in plots? 
#' @param by Variable name to split or plot by. 'time', 'LogLik', and 'WhichObs' are also possibilities. 
#' @param TIpredNames Since time independent predictors do not change with time, by default observations after the first are ignored. 
#' For observing attrition it can be helpful to set this to NULL, or when the combinevars argument is used, specifying different names may be useful.
#' @param nsamples Number of samples (when applicable) to include in plots.
#' @param covplot Splits variables in the model by the 'by' argument, according to the number of breaks (breaks argument), 
#' and shows the covariance (or correlation) for the different data sources selected, as well as the differences between each pair. 
#' @param corr Turns the covplot into a correlation plot. Usually easier to make sense of visually. 
#' @param combinevars Can be a list of (possibly new) variable names, where each named element of the list contains a character vector 
#' of one or more variable names in the fit object, to combine into the one variable. By default, the mean is used, but see the aggfunc argument.
#'  The combinevars argument can also be used to ensure that only certain variables are plotted.  
#' @param fastcov Uses base R cov function for computing covariances. Not recommended with missing data.
#' @param aggfunc Function to use for aggregation, if needed.
#' @param aggregate If TRUE, duplicate observation types are aggregated over using aggfunc. For example,
#' if by = 'time' and there are 8 time points per subject, but breaks = 2, 
#' there will be 4 duplicate observation types per 'row' that will be collapsed. In most cases it is helpful to not collapse. 
#' @param groupbysplit Logical. Affects variable ordering in covariance plots. 
#' Defaults to FALSE, grouping by variable, and within variable by split.
#' @param byNA Logical. Create an extra break for when the split variable is missing?
#' @param lag Integer vector. lag = 1 creates additional variables for plotting, prefixed by 'lag1_', containing the prior row of observations
#' for that subject.
#' @param smooth For bivariate plots, use a smoother for estimation?
#' @param k Integer denoting number of knots to use in the smoothing spline.
#' @param breaks Integer denoting number of discrete breaks to split variables by (when covariance plotting).
#' @param entropy Still in development. 
#' @param reg Logical. Use regularisation when estimating covariance matrices? Can be necessary / faster for some problems.
#' @param verbose Logical. If TRUE, shows optimization output when estimating covariances.
#' @param indlines Integer number of individual subject lines to draw per data type. 
#'
#' @return Nothing. Just plots. 
#' @export
#'
#' @examples
#' if(w32chk()){
#' ctCheckFit(ctstantestfit)
#' }
ctCheckFit <- function(fit, 
  data=TRUE, postpred=TRUE, priorpred=FALSE, statepred=FALSE, residuals=FALSE,
  by=fit$ctstanmodelbase$timeName,
  TIpredNames=fit$ctstanmodelbase$TIpredNames,
  nsamples=10, covplot=FALSE, corr=TRUE, combinevars=NA, fastcov=FALSE,
  lagcovplot=FALSE,
  aggfunc=mean,aggregate=TRUE,
  groupbysplit=FALSE, byNA=TRUE,lag=0,
  smooth=TRUE, k=4,breaks=4,entropy=FALSE,reg=FALSE,verbose=0, indlines=30){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
  covORcor <- function(m){
    if(corr) return(cov2cor(m)) else return(m)
  }
  
  DataSource <-Sample<-NULL
  
  dat <- ctStanFitMelt(fit = fit,maxsamples = nsamples)
  
  
  
  byc=unique(c('Sample','WhichObs','DataSource',fit$ctstanmodel$timeName))
  bycid=c(fit$ctstanmodelbase$subjectIDname,byc)

  if(is.na(combinevars[1])){
    combinevars = setNames(
      colnames(dat)[!colnames(dat) %in% bycid],colnames(dat)[!colnames(dat) %in% bycid])
    # combinevars<-c(combinevars,LogLik='LogLik')
  }
  dat <- ctDataMelt(dat=dat,id=fit$ctstanmodelbase$subjectIDname, by=byc,combinevars = combinevars)
  dat$Sample <- factor(dat$Sample)
  dat$DataSource <- factor(dat$DataSource)

  # ?acf
  
  
  wdat=dcast(dat,paste0(
   paste0(bycid,collapse='+'),
    '~variable'),value.var = 'value',fun.aggregate = mean,na.rm=TRUE)
  
  if(any(lag!=0)){ 
    if(aggregate){
      if(covplot){
        message('disabling aggregation -- incoherent with lags!')
        aggregate=FALSE
      }
    }
    prelagcols=colnames(wdat)[!colnames(wdat) %in% c(bycid,TIpredNames)]
    lagcols=(paste0("lag",  rep(lag, times = length(prelagcols)),'_',rep(prelagcols, each = length(lag))))
    
    wdat[, eval(lagcols) :=  
        shift(.SD, n=lag,type='lead'), by=c(fit$ctstanmodel$subjectIDname,'Sample','DataSource'),.SDcols=prelagcols]
  }
  
  dat <- melt(wdat,id.vars=unique(c(bycid,by)))
  
  if(!data) dat<-dat[!DataSource %in% 'Data']
  if(!priorpred) dat<-dat[!DataSource %in% 'PriorPred']
  if(!postpred) dat<-dat[!DataSource %in% 'PostPred']
  if(!statepred) dat<-dat[!DataSource %in% 'StatePred']
  if(!residuals) dat<-dat[!DataSource %in% 'Residuals']
  
  # if(!'WhichObs' %in% by) mdat$WhichObs <- NULL
  
  dat = dat[!is.na(value),]
  

  
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
        
        if(groupbysplit){
          colorder <- unlist(lapply(breakset, function(b) grep(paste0('\\_',b,'$'),x=colnames(wdat))))
          wdat <- wdat[,colorder,with=FALSE]
        } else if(!byNA) wdat <- wdat[,-grep('\\_NA%',colnames(wdat)),with=FALSE] #if not seperating by group, still remove NA by relations if needed
        
        sapply(TIpredNames,function(x){ #fix tipred column names
          if(paste0(x,'_1') %in% colnames(wdat)){
            colnames(wdat)[ colnames(wdat) %in% paste0(x,'_1')] <<- x
          }
        })
        if(breaks==1) colnames(wdat) <- gsub("\\_\\d+$","",colnames(wdat))
        wdat <- cbind(wdat[,c(which(colnames(wdat) %in% TIpredNames), which(!colnames(wdat) %in% TIpredNames)),with=FALSE])
        
        if(!fastcov) corlist[[dsi]] <- covml(
          data.frame(wdat)[,!colnames(wdat) %in% c('Sample',fit$ctstanmodelbase$subjectIDname),drop=FALSE],
          reg=reg,
          verbose=verbose)$cp$covm
        
        if(fastcov) corlist[[dsi]] <- cov(
          data.frame(wdat)[,!colnames(wdat) %in% c('Sample',fit$ctstanmodelbase$subjectIDname),drop=FALSE],
          use='pairwise.complete.obs')
        
        if(corr) corlist[[dsi]] <- cov2cor(corlist[[dsi]])
        
        # if(covsplitdiffs){
        #   
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
      
      if(corr) covplotlims <- c(-1,1) else covplotlims <- NA
      
      if(!lagcovplot){
      for(i in seq_along(corlist)){ #regular corplots
        cplot=corplotmelt(meltcov(covORcor(corlist[[datasources[i]]])),
          label=paste0(ifelse(corr,'Corr.','Cov.')),limits=covplotlims)
        cplot = cplot + ggplot2::labs(title=paste0(datasources[i],
          ifelse(corr,' correlations',' covariances'),if(breaks > 1) paste0(', split by ', by)))
        print(cplot)
      }
      for(i in seq_along(corlist)){
        for(j in seq_along(corlist)){ #difference corplots
          if(j > i){ #only plot unique differences
            coli <- colnames(corlist[[datasources[i]]])[
              colnames(corlist[[datasources[i]]]) %in% colnames(corlist[[datasources[j]]])
            ]
            cplot=corplotmelt(meltcov(corlist[[datasources[i]]][coli,coli,drop=FALSE] - 
                corlist[[datasources[j]]][coli,coli,drop=FALSE]),
              label=paste0(ifelse(corr,'Corr.','Cov.'),' diff.'),limits=covplotlims)
            cplot = cplot + ggplot2::labs(title=paste0(datasources[i],' - ',datasources[j],
              ifelse(corr,' correlations',' covariances'),if(breaks > 1) paste0(', split by ', by)))
            print(cplot)
          }
        }
      }
      }
     
      if(lagcovplot){
        for(ci in seq_along(corlist)){
          if(ci==1) lagdat <-data.frame(corlist[[ci]],DataSource = names(corlist)[ci], 
            Row=rownames(corlist[[ci]])) else{
          lagdat <- rbind(lagdat,data.frame(corlist[[ci]],
            DataSource = names(corlist)[ci],Row=rownames(corlist[[ci]])))
          }
        }
        lagdat$Lag <- 0
        l=regmatches(x = lagdat$Row,m = regexpr(pattern = '^lag\\d+_',lagdat$Row))
        l=gsub('lag','',l)
        l=gsub('\\_','',l)
        lagdat$Lag[grepl('^lag\\d+\\_',lagdat$Row)] <- as.numeric(l)
        lagdat <- lagdat[,-grep('^lag\\d+\\_',colnames(lagdat))]
        lagdat$Row <- gsub('^lag\\d+\\_','',lagdat$Row)
        lagdat <- melt(data.table(lagdat),id.vars = c('Row','Lag','DataSource'))
        lagdat$variable <- gsub('\\_1$','',lagdat$variable)
        lagdat <- lagdat[!lagdat$Row %in% TIpredNames,]
        
        combinevars <- combinevars[!names(combinevars) %in% c(TIpredNames,by)] #because discretized / not in cov

        for(vi in names(combinevars)){    
g=ggplot(data = lagdat[lagdat$variable %in% vi,],
  mapping = aes(x=Lag,y=value,colour=DataSource, linetype=DataSource))+geom_line(size=1)+
  facet_wrap(facets = vars(Row))+theme_minimal()+ylab(label = ifelse(corr,'Correlation','Covariance'))+
  scale_y_continuous(minor_breaks = c(-.5,.5),breaks=seq(-1,1,1))+
  scale_x_continuous(minor_breaks=NULL,breaks=seq(0,max(lag), ceiling(max(lag)/5)))+
    ggtitle(label = paste0(vi,' ',ifelse(corr,'Correlations','Covariances')))
if(corr) g <- g+coord_cartesian(ylim=c(-1,1))
print(g)
# if('try-error' %in% class(g) ) browser()
        }

        
      }
      
    }
    # if(entropy){
    #   # 
    #   dentropy <- -mean(covdat$cp$llrow)#lp/nrow(wdat)
    #   
    #   gdentropy <- sapply(unique(wgdat$Sample), function(x){
    #     print(x)
    #     -mean(covml(data.frame(wgdat)[wgdat$Sample==x,!colnames(wgdat) %in% 
    #         c('Sample',fit$ctstanmodelbase$subjectIDname)],verbose=verbose)$cp$llrow)
    #   }) 
    #   print(sum(covdat$cp$llrow))
    #   message('Entropy: Fit = ', -sum(fit$stanfit$transformedparsfull$ll)/nrow(wdat),'; Data = ',dentropy,'; Mean gen.data = ',mean(gdentropy),'  SD gen. data = ',sd(gdentropy))
    # }
  }
  
  if(!covplot){ #then plot bivariate relations
    vars <- names(combinevars)
    if(any(lag>0)) vars <- c(vars,lagcols)
    dat$manifest <- dat$variable %in% vars
    
    
    g=ggplot(data = dat[manifest==TRUE],
      mapping = aes_string(x=by,y='value',
        colour='DataSource'
        ,fill='DataSource'
      )
    ) +theme_bw() +
      # scale_color_brewer(palette = 'Set1') +
      # scale_fill_brewer(palette = 'Set1')
      # +
      scale_colour_manual(values = c("blue", "red","green",  "Orange","Black"))+
      scale_fill_manual(values = c("blue", "red","green",  "Orange",'Black'))#+
    # scale_alpha_manual(values = c(1, max(.05,(.3/nsamples)),max(.05,(.3/nsamples)),1))
    
    if(smooth) {
      # g = g + geom_smooth(data=dat[manifest==TRUE & !DataSource %in% c('Data','StatePred','Residuals')],
      #   aes(group=Sample),
      #   alpha= max(.05,(.3/nsamples)),
      #   linetype=0,#3,#ifelse(nsamples==1,1,0),
      #   stat="smooth",#se = FALSE,#nsamples==1,
      #   size=.1
      #   ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
      
      g = g + geom_smooth(data=dat[manifest==TRUE 
        # & DataSource %in% c('Data','StatePred','Residuals')
        ],
        # aes(colour=DataSource,alpha=NULL),
        stat="smooth",se = TRUE,size=1,alpha=.3
        ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
      
      
  
    }
    if(!smooth) {
      # if(nsamples > 1) g = g + stat_summary(fun=mean,geom = "line",size=1) #geom_line(stat=mean)
      
      # if(nsamples ==1)
      # g = g +  stat_summary(data=dat[manifest==TRUE 
      #   & !DataSource %in% c('Data','StatePred','Residuals')],
      #   aes(group=Sample),
      #   fun.data = function(x) list(y=mean(x),
      #     ymin=mean(x,na.rm=TRUE)-sd(x)/sqrt(length(x)), 
      #     ymax=mean(x)+sd(x)/sqrt(length(x))),
      #   geom = "ribbon",
      #   alpha= max(.05,max(.05,sqrt(.2/nsamples))),
      #   linetype=ifelse(nsamples==1,1,0))
      
      g = g +  stat_summary(data=dat[manifest==TRUE 
        # & DataSource %in% c('Data','StatePred','Residuals')
        ],
        fun.data = function(x) list(y=mean(x),
          ymin=mean(x,na.rm=TRUE)-sd(x)/sqrt(length(x)), 
          ymax=mean(x)+sd(x)/sqrt(length(x))),
        geom = "ribbon",alpha=.3)
    }
    
    if(indlines > 0){
      individualInt <- eval(parse(text=paste0('interaction(dat$',fit$ctstanmodelbase$subjectIDname,',dat$Sample, dat$DataSource)')))
      dat$individualInt <- individualInt
      ids <- unique(dat[[fit$ctstanmodelbase$subjectIDname]])
      ids=sample(ids,  min(indlines,length(ids)),replace = FALSE)
      ids <- unlist(lapply(unique(dat$DataSource),function(ds){
        dat$individualInt[dat$DataSource %in% ds][ 
          match(ids, dat[[fit$ctstanmodelbase$subjectIDname]][dat$DataSource %in% ds])]
      }))
      g = g + geom_line(data=dat[manifest==TRUE &dat$individualInt %in% ids],
        mapping = aes(group=individualInt),
        alpha= max(.2,(.8/sqrt(indlines))),
        linetype=1,#3,#ifelse(nsamples==1,1,0),
        # stat="smooth",#se = FALSE,#nsamples==1,
        size=.1
        # ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')'))
      )
      
    }
    g=g+facet_wrap(facets = vars(variable),scales = 'free')
    print(g)
  }
}


corplotmelt <- function(corm, label='Coef.',limits=NA,title=''){
  
  if(is.na(limits[1])) limits <- range(corm[,'value'],na.rm=TRUE)
  ggplot(data=(corm),aes_string(x='Var1',y='Var2',fill=('value')))+ #ifelse(groups,NULL,'Group')))+
    geom_tile( width=1,height=1,colour='black')+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
      midpoint = 0, limits = limits, space = "Lab", 
      name=label)  +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90)) + ggtitle(title)
}


