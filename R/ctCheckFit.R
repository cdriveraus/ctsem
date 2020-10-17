ctLongToWideSF <- function(fit){
  dat <- data.frame(standatatolong(standata = fit$standata,
    ctm = fit$ctstanmodel,origstructure = TRUE))
  # dat[[fit$ctstanmodelbase$subjectIDname]] <- NULL
  # dat[[fit$ctstanmodelbase$timeName]] <- NULL
  dat <- dat[,c(fit$ctstanmodelbase$subjectIDname,
    fit$ctstanmodelbase$timeName,
    fit$ctstanmodelbase$manifestNames,
    if(fit$ctstanmodelbase$n.TDpred > 0) fit$ctstanmodelbase$TDpredNames,
    if(fit$ctstanmodelbase$n.TIpred > 0) fit$ctstanmodelbase$TIpredNames)]
  dat=data.table(dat)
  dat[ ,WhichObs:=0:(.N-1),by=eval(fit$ctstanmodelbase$subjectIDname)]
  dat=melt(dat,id.vars = c(fit$ctstanmodelbase$subjectIDname,'WhichObs'))
  if(fit$ctstanmodelbase$n.TIpred > 0){
    dat <- dat[WhichObs <1 | !variable%in% eval(fit$ctstanmodelbase$TIpredNames)]
  }
  dat$WhichObs <- paste0('T',dat$WhichObs);
  dat=dcast(dat,'id~variable+WhichObs')
  lapply(fit$ctstanmodelbase$TIpredNames, function(x){
    colnames(dat)[grep(paste0('\\b',x,'\\_T0'),colnames(dat))] <<- x
  })
  return(data.frame(dat))
}

ctSaturatedFitConditional<-function(dat,ucols,reg,covf=NA,verbose=0){
  o1 = ucols
  o2=(1:ncol(dat))[-o1]
  
  if(all(is.na(covf))) covf=covml(dat,reg = reg,verbose=verbose)
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
    covf$ll <- sum(llrow)
  } else { #if not conditional
    covf$ll <- sum(sapply(1:nrow(dat),function(i){
    llr=NA
    d1 = !is.na(dat[i,o1])
    if(sum(d1)>0){
      mu=c(covf$cp$mu[d1])
      llr=mvtnorm::dmvnorm(x = dat[i,o1][d1],mean = mu,sigma = covf$cp$covm[d1,d1,drop=FALSE],log=TRUE)
    }
    return(llr)
  }))
  }
  
  return(covf)
}

ctSaturatedFit <- function(fit,conditional=FALSE,reg=FALSE, 
  time=FALSE, oos=TRUE, folds=10, cores=2,verbose=0){
  dat <-ctLongToWideSF(fit)
  dat=dat[,-1]
  if(!time) dat <- dat[,-grep(paste0('\\b',fit$ctstanmodelbase$timeName,'_T'),colnames(dat))]
  
  #get unconditional columns
  ucols=unique(c(unlist(sapply(fit$ctstanmodelbase$manifestNames,function(x){
    grep(paste0('\\b',x,'\\_T'),colnames(dat))
  }))))
  
  if(cores > 1){
    cl=parallel::makeCluster(cores,'PSOCK')
    parallel::clusterExport(cl,c('ucols','dat','reg'),environment())
    on.exit({parallel::stopCluster(cl)},add = TRUE)
  }
  
  bootstrap=FALSE
  srows <- sample(1:nrow(dat),
    nrow(dat),# *ifelse(bootstrap,folds,1),
    replace=bootstrap)
  
  srows <- split(
    srows,
    sort(1:length(srows) %% folds))
  
  if(!conditional){
    dat <- data.frame(dat)[,unique(c(sapply(fit$ctstanmodelbase$manifestNames,function(x){
      grep(paste0('\\b',x,'\\_T'),colnames(dat))
    }))),drop=FALSE]
    o1=1:ncol(dat)
  }
  
  sf = flexlapply(cl,srows,function(x){
    library(ctsem)
    datheldout <- dat
    datheldout[x,ucols] <- NA
    fit=ctSaturatedFitConditional(dat = datheldout,ucols = ucols,reg = reg)
    # covdata <- covdata(dat[x,,drop=FALSE],reg=reg)
    fitoos=ctSaturatedFitConditional(dat = dat[x,,drop=FALSE],ucols = ucols,reg = reg,covf=fit)
      # suppressMessages(sampling(object = stanmodels$cov,iter=1,chains=0,check_data=FALSE,data=covdata))
    # f$lloos=rstan::constrain_pars(sfoos, f$fit$par)$llrow
    fit$lloos <- fitoos$ll
    return(fit)
  },cores = cores)
  # browser()
  
  covf = ctSaturatedFitConditional(dat = dat,ucols = ucols,reg = reg)
  covf$lloos = sum(unlist(lapply(sf,function(x) x$lloos)),na.rm=TRUE)
  
  
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
  # ldat[ ,WhichObs:=(1:.N),by=Subject]
  
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




ctCheckFit2 <- function(fit, by=fit$ctstanmodelbase$timeName,
  TIpredNames=fit$ctstanmodelbase$TIpredNames,
  nsamples=10, covplot=TRUE, combinevars=NULL,
  smooth=TRUE, k=10,breaks=4,entropy=covplot,reg=FALSE,verbose=0){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object')
  if(is.null(fit$generated)) stop('First use ctStanGenerateFromFit() on fit object!')
  if(nsamples > dim(fit$generated$Y)[2]) stop('Not enough samples in generated values! Try fewer samples?')
  
  dat <- data.table(standatatolong(standata = fit$standata,ctm = fit$ctstanmodel,origstructure = TRUE))
  dat$Sample <- 1
  dat[ ,WhichObs:=(1:.N),by=eval(fit$ctstanmodelbase$subjectIDname)]
  # dat$WhichObs <- as.integ(dat$WhichObs)
  
  gdat <- data.frame(dat[rep(seq_len(nrow(dat)), nsamples)])
  samples <- sample(1:dim(fit$generated$Y)[2], nsamples,replace = FALSE)
  gdat[,(fit$ctstanmodelbase$manifestNames)] <- data.frame(matrix(
    fit$generated$Y[, samples, , drop=FALSE], 
    # nrow=dim(fit$generated$Y)[1],
    ncol=dim(fit$generated$Y)[3]))
  gdat$Sample=rep(samples,each=dim(fit$generated$Y)[1])
  gdat <- data.table(gdat)
  
  if(is.null(combinevars)) combinevars = setNames(
    lapply(fit$ctstanmodelbase$manifestNames,function(x) x),fit$ctstanmodelbase$manifestNames)
  
  dat <- ctDataMelt(dat=dat,id=fit$ctstanmodelbase$subjectIDname, by=c(by,'Sample','WhichObs'),combinevars = combinevars)
  gdat <- ctDataMelt(dat=gdat,id=fit$ctstanmodelbase$subjectIDname, by=c(by,'Sample','WhichObs'),combinevars = combinevars)
  # browser()
  
  dat <- dat[WhichObs <2 | !variable%in% eval(TIpredNames)]
  gdat <- gdat[WhichObs <2 | !variable%in% eval(TIpredNames)]
  
  # dat<- data.frame(dat)
  # dat <- dat[!(dat$WhichObs > 1 & dat$variable %in% TIpredNames),]
  # dat=data.table(dat)
  # 
  # gdat<- data.frame(gdat)
  # gdat <- gdat[!(gdat$WhichObs > 1 & gdat$variable %in% TIpredNames),]
  # gdat=data.table(gdat)
  if(!'WhichObs' %in% by){
    dat$WhichObs <- NULL
    gdat$WhichObs <- NULL
  }
  
  dat = dat[!is.na(value),]
  gdat = gdat[!is.na(value),]
  
  if(covplot ||entropy){
    if(is.double(dat[[by]])){
      datb=dat
      gdatb=gdat
      dat[[by]] <- arules::discretize(dat[[by]],method='cluster',breaks = breaks,labels=FALSE)
      gdat[[by]] <- arules::discretize(gdat[[by]],method='cluster',breaks = breaks,labels=FALSE)
    }
    if(covplot){
      # browser()
      wdat <- dcast(dat,paste0(fit$ctstanmodelbase$subjectIDname, '~variable + ',by),fun.aggregate=mean,na.rm=TRUE)
      # wdat <- wdat[,apply(wdat,2,function(x) any(!is.na(x))),drop=FALSE]
      covdat <- covml(
        data.frame(wdat)[,!colnames(wdat) %in% fit$ctstanmodelbase$subjectIDname,drop=FALSE],
        reg=reg,
        verbose=verbose)#,
      # use='pairwise.complete.obs')
      
      wgdat <- dcast(gdat,paste0(fit$ctstanmodelbase$subjectIDname, '+Sample~variable + ',by),fun.aggregate=mean,na.rm=TRUE)
      # wdat <- wdat[,apply(wdat,2,function(x) any(!is.na(x))),drop=FALSE]
      covgdat <- covml(
        data.frame(wgdat)[,!colnames(wgdat) %in% c('Sample',fit$ctstanmodelbase$subjectIDname)],
        reg=reg,
        verbose=verbose)#,
      # use='pairwise.complete.obs')
      # covdatmx <- OpenMx::mxRefModels(data.frame(wdat)[,!colnames(wdat) %in% fit$ctstanmodelbase$subjectIDname,drop=FALSE],run=TRUE)
      # covgdat / covdat
      # browser()
      print(corplotmelt(reshape2::melt(cov2cor(covdat$cp$covm)),label='Orig. corr.',limits=c(-1,1)))
      print(corplotmelt(reshape2::melt(cov2cor(covgdat$cp$covm)),label='Gen. corr.',limits=c(-1,1)))
      
      print(corplotmelt(reshape2::melt(cov2cor(covgdat$cp$covm)-cov2cor(covdat$cp$covm)),label='Gen. corr. diff.',limits=c(-1,1)))
    }
    if(entropy){
      # browser()
      # wdat <- data.frame(wdat)[,!(colnames(wdat) %in%  #remove covariates not used in fit lp
      #     c(fit$ctstanmodelbase$TIpredNames,fit$ctstanmodelbase$TDpredNames))]
      # wgdat <- data.frame(wgdat)[,!(colnames(wdat) %in%
      #     c(fit$ctstanmodelbase$TIpredNames,fit$ctstanmodelbase$TDpredNames))]
      
      # mxRefModels(wdat,)
      dentropy <- -mean(covdat$cp$llrow)#lp/nrow(wdat)
      
      gdentropy <- sapply(unique(wgdat$Sample), function(x){
        print(x)
        -mean(covml(data.frame(wgdat)[wgdat$Sample==x,!colnames(wgdat) %in% 
            c('Sample',fit$ctstanmodelbase$subjectIDname)],verbose=verbose)$cp$llrow)
      }) 
      print(sum(covdat$cp$llrow))
      message('Entropy: Fit = ', -sum(fit$stanfit$transformedparsfull$ll)/nrow(wdat),'; Data = ',dentropy,'; Mean gen.data = ',mean(gdentropy),'  SD gen. data = ',sd(gdentropy))
    }
    if(is.double(dat[[by]])){
      dat=datb
      gdat=gdatb
    }
  } 
  
  gdat$Original <- FALSE
  dat$Original <- TRUE
  dat <- rbind(dat,gdat)
  dat$Sample <- factor(dat$Sample)
  
  if(!is.null(combinevars)) vars <- names(combinevars) else vars <- fit$ctstanmodelbase$manifestNames
  
  dat$manifest <- dat$variable %in% vars
  # if(varapprox) dat[manifest==TRUE,value:=value^2]# <- dat[manifest==TRUE,value]^2
  
  # gam
  
  g=ggplot(data = dat[manifest==TRUE],
    mapping = aes_string(x=by,y='value',
      colour='Original',
      group='Sample',
      alpha='Original')) +theme_bw()+
    # geom_point(alpha=.1)+
    scale_alpha_ordinal(range = c(0.3, 1))
  
  if(smooth) {
    g = g + geom_line(data=dat[manifest==TRUE & Original==FALSE],aes(alpha=NULL),alpha=.3,
      stat="smooth",se = FALSE,size=1
      ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
    g = g + geom_smooth(data=dat[manifest==TRUE & Original==TRUE],aes(alpha=NULL),
      stat="smooth",se = TRUE,size=1,alpha=.3
      ,method='gam', formula= as.formula(paste0('y ~ s(x,bs="cr",k=',k,')')))
  }
  if(!smooth) {
    g = g + stat_summary(fun=mean,geom = "line",size=1) #geom_line(stat=mean)
    g = g +  stat_summary(data=dat[manifest==TRUE & Original==TRUE],
      fun.data = function(x) list(y=mean(x),
        ymin=mean(x,na.rm=TRUE)-sd(x)/sqrt(length(x)), 
        ymax=mean(x)+sd(x)/sqrt(length(x))),
      # fun.min = function(x) mean(x) - sd(x), 
      # fun.max = function(x) mean(x) + sd(x), 
      geom = "ribbon",alpha=.2)
  }
  g=g+facet_wrap(facets = vars(variable),scales = 'free')
  print(g)
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
    # browser()
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
      ygen <- aperm(ctStanGenerateFromFit(fit,fullposterior=FALSE,nsamples=niter)$generated$Y,c(2,1,3)) #array(e$Ygen,dim=c(ygendim[1] * ygendim[2],ygendim[-1:-2]))
    } else ygen <- aperm(fit$generated$Y,c(2,1,3))
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
  # browser()
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
      # browser()
      
      corm <- reshape2::melt(ggcorrArgs$cor_matrix)
      
      corplotmelt(corm,'MissSpec. Ratio',limits)
      
      
      # do.call(GGally::ggcorr,ggcorrArgs) #(data=NULL,cor_matrix =  get(covtype),limits=limits, geom = 'circle',max_size = 13,name=covtype,...)
    }
  }
  
}

corplotmelt <- function(corm, label,limits=NULL){
  ggplot(data=(corm),aes(x=Var1,y=Var2,fill=(value)))+geom_tile()+
    geom_tile(color = "black")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
      midpoint = 0, limits = limits, space = "Lab", 
      name=label)  +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90))
}
