#' Check absolute fit of ctFit or ctStanFit object.
#'
#' @param fit ctFit or ctStanFit object.
#' @param niter number of data generation iterations to use to calculate quantiles.
#' @param probs 3 digit vector of quantiles to return and to test significance.
#'
#' @return numeric matrix showing fit for each lower triangular index of the covariance matrix of data.
#' @export
#'
#' @examples
#' \dontrun{
#' data(ctExample1)
#' traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
#'   manifestNames=c('LeisureTime', 'Happiness'), 
#'   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
#' 
#' check <- ctCheckFit(traitfit,niter=5)
#' plot(check)
#' }
ctCheckFit <- function(fit, niter=50,probs=c(.025,.5,.975)){
  
  if(class(fit)=='ctsemFit'){
    if('Kalman' %in% fit$ctfitargs$objective) {
      wdat <- ctLongToWide(fit$mxobj@data$observed,id='id',time='time',
        manifestNames = fit$ctmodelobj$manifestNames)[,paste0(fit$ctmodelobj$manifestNames,'_T',1),drop=FALSE]
    } else  wdat <- fit$mxobj@data$observed[,paste0(fit$ctmodelobj$manifestNames,'_T',
      rep(0:(fit$ctmodelobj$Tpoints-1),each=fit$ctmodelobj$n.manifest)),drop=FALSE]
  }
  
  if(class(fit)=='ctStanFit') {
    ldat <- cbind(fit$data$subject,fit$data$time,fit$data$Y)
    tpoints <- max(unlist(lapply(unique(fit$data$subject),function(x) length(fit$data$subject[fit$data$subject==x]))))
    colnames(ldat)[1:2] <- c('subject','time')
    wdat <- ctLongToWide(ldat,id='subject',time='time',
      manifestNames = fit$ctstanmodel$manifestNames)[,paste0(fit$ctstanmodel$manifestNames,'_T',
        rep(0:(tpoints-1),each=fit$ctstanmodel$n.manifest)),drop=FALSE]
  }
  
  ecov <- cov(wdat,use = "pairwise.complete.obs")
  
  if(class(fit)=='ctStanFit'){
    e <- extract.ctStanFit(fit)
    ygendim <- dim(e$Ygen)
    ygen <- array(e$Ygen,dim=c(ygendim[1] * ygendim[2],ygendim[-1:-2]))
    covarray<-array(NA,dim = c(dim(ecov),dim(ygen)[1]))

    for(i in 1:dim(ygen)[1]){
    ndat <- cbind(fit$data$subject,fit$data$time,ygen[i,,])
    colnames(ndat) <- c('subject','time',fit$ctstanmodel$manifestNames)
    ndat <- suppressMessages(ctLongToWide(ndat,id='subject',time='time',
      manifestNames = fit$ctstanmodel$manifestNames)) #[,paste0( fit$ctstanmodel$manifestNames,'_T',1),drop=FALSE]
     covarray[,,i] <- cov(ndat[,paste0(fit$ctstanmodel$manifestNames,'_T',
        rep(0:(tpoints-1),each=fit$ctstanmodel$n.manifest)),drop=FALSE], use='pairwise.complete.obs')
    }
  }
  
  if(class(fit)=='ctsemFit'){
    covarray<-array(NA,dim = c(dim(ecov),niter))
    for(i in 1:niter){
      ndat <- ctGenerateFromFit(fit = fit,n.subjects = nrow(wdat))
      ndat[is.na(wdat)] <- NA #match missingness
      covarray[,,i] <- cov(ndat[,paste0(fit$ctmodelobj$manifestNames,'_T',
        rep(0:(fit$ctmodelobj$Tpoints-1),each=fit$ctmodelobj$n.manifest)),drop=FALSE], use='pairwise.complete.obs')
    }
  }
  
  covql <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[1],na.rm=TRUE)
  covqm <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[2],na.rm=TRUE)
  covqh <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[3],na.rm=TRUE)
  covmean <- ctCollapse(covarray,collapsemargin = 3,mean,na.rm=TRUE)
  
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
  MisspecRatio <- (test[,'observed'] - test[,'50%']) / (test[,'97.5%'] - test[,'2.5%'])
  test<- cbind(rowname,colname,as.data.frame(test),MisspecRatio)
  class(test) <- c('ctsemFitMeasure',class(test))
  return(test)
}

#' Misspecification plot using ctCheckFit output
#'
#' @param x Object output from ctsemFitMeasure function.
#' @param corrplotargs Extra arguments to pass to corrplot function.
#' @param ... not used.
#'
#' @return Nothing, just plots.
#' @export
#' @importFrom corrplot corrplot
#' @method plot ctsemFitMeasure
#'
#' @examples
#' \dontrun{
#' data(ctExample1)
#' traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
#'   manifestNames=c('LeisureTime', 'Happiness'), 
#'   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
#' 
#' check <- ctCheckFit(traitfit,niter=5)
#' plot(check)
#' }
plot.ctsemFitMeasure <- function(x,corrplotargs = list(method='square',is.corr=FALSE),...){
  ratiomat <- matrix(NA,max(x[,'row']),max(x[,'row']))
  ratiomat[upper.tri(ratiomat,diag = TRUE)] = x[,'MisspecRatio']
  ratiomat[lower.tri(ratiomat)] = t(ratiomat)[lower.tri(ratiomat)]
  
  colnames(ratiomat) <- unique(x[,'colname'])
  rownames(ratiomat) <- unique(x[,'rowname'])
  
  corrplotargs$corr <- ratiomat
  do.call(corrplot,corrplotargs)
}
