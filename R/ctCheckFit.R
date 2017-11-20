#' Check absolute fit of ctFit object.
#'
#' @param fit ctFit object.
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
#' print(ctCheckFit(traitfit))
#' }
ctCheckFit <- function(fit, niter=50,probs=c(.025,.5,.975)){
  if('Kalman' %in% fit$ctfitargs$objective) {
    wdat <- ctLongToWide(fit$mxobj@data$observed,
      manifestNames = fit$ctmodelobj$manifestNames)[,paste0(fit$ctmodelobj$manifestNames,'_T',1),drop=FALSE]
  } else  wdat <- fit$mxobj@data$observed[,paste0(fit$ctmodelobj$manifestNames,'_T',
    rep(0:(fit$ctmodelobj$Tpoints-1),times=fit$ctmodelobj$n.manifest)),drop=FALSE]
ecov <- cov(wdat)
covarray<-array(NA,dim = c(dim(ecov),niter))

for(i in 1:niter){
  ndat <- ctGenerateFromFit(fit = fit,n.subjects = nrow(wdat))
  covarray[,,i] <- cov(ndat[,paste0(fit$ctmodelobj$manifestNames,'_T',
    rep(0:(fit$ctmodelobj$Tpoints-1),each=fit$ctmodelobj$n.manifest)),drop=FALSE])
}
# browser()
covql <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[1])
covqm <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[2])
covqh <- ctCollapse(covarray,collapsemargin = 3,quantile,probs=probs[3])
covmean <- ctCollapse(covarray,collapsemargin = 3,mean)

test<-matrix(NA,ncol=8,nrow=(nrow(covql)^2+nrow(covql))/2)
counter=0
for(i in 1:nrow(covql)){
  for(j in 1:nrow(covql)){
    if(j <=i){
      counter=counter+1
      test[counter,] <- c(i,j,covmean[i,j],covql[i,j],covqm[i,j],covqh[i,j],ecov[i,j],
        ifelse((ecov[i,j] > covqh[i,j] || ecov[i,j] < covql[i,j]), TRUE,FALSE))
    }}}
colnames(test) <- c('row','col','mean',paste0(probs*100,'%'), 'observed', 'significant')

return(test)
}
