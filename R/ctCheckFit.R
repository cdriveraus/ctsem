#' Check absolute fit of ctFit or ctStanFit object.
#'
#' @param fit ctFit or ctStanFit object.
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
#' data(ctExample1)
#' traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
#'   manifestNames=c('LeisureTime', 'Happiness'), 
#'   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' traitfit <- ctFit(dat=ctExample1, ctmodelobj=traitmodel)
#' 
#' check <- ctCheckFit(traitfit,niter=5)
#' plot(check, wait=FALSE)
#' }
ctCheckFit <- function(fit, niter=500,probs=c(.025,.5,.975)){
  
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
      ygen <- aperm(ctStanGenerateFromFit(fit,fullposterior=TRUE,nsamples=niter)$generated$Y,c(2,1,3)) #array(e$Ygen,dim=c(ygendim[1] * ygendim[2],ygendim[-1:-2]))
    } else ygen <- aperm(fit$generated$Y,c(2,1,3))
    wide <- matrix(NA, nrow=length(unique(fit$data$subject)),ncol=dim(ygen)[3]*maxtp)
    itervec <- sample(1:dim(ygen)[1],niter)
    dimnames(ygen)<-list(iter=1:dim(ygen)[1],row=1:dim(ygen)[2],manifestNames)
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
    for(i in 1:niter){
      ndat <- ctGenerateFromFit(fit = fit,n.subjects = nrow(wdat))
      ndat[is.na(wdat)] <- NA #match missingness
      covarray[,,i] <- cov(ndat[,paste0(rep(manifestNames,each=fit$ctmodelobj$Tpoints),'_T',
        0:(fit$ctmodelobj$Tpoints-1)),drop=FALSE], use='pairwise.complete.obs')
      means[,,i] <- t(matrix(apply(ndat[,paste0(rep(manifestNames,each=fit$ctmodelobj$Tpoints),'_T',
        0:(fit$ctmodelobj$Tpoints-1)),drop=FALSE],2,mean,na.rm=TRUE),byrow=TRUE,ncol=nmanifest))
    }
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
#' \donttest{
#' 
#' if (!exists("ctstantestfit")) ctstantestfit <- ctstantestfitgen()
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
      do.call(GGally::ggcorr,ggcorrArgs) #(data=NULL,cor_matrix =  get(covtype),limits=limits, geom = 'circle',max_size = 13,name=covtype,...)
    }
  }
  
}
