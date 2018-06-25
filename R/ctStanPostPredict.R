
#' Compares model implied density and values to observed, for a ctStanFit object.
#'
#' @param fit ctStanFit object.
#' @param legend Logical, whether to plot a legend.
#' @param wait Logical, if TRUE, waits for input before plotting next plot.
#' @return Nothing, only plots.
#' @export
#' @details This function relies on the data generated during each iteration of fitting to approximate the
#' model implied distributions -- thus, when limited iterations are available, the approximation will be worse.
#'
#' @examples
#' ctStanPostPredict(ctstantestfit)
ctStanPostPredict <- function(fit,legend=TRUE,wait=FALSE){
  e<-extract.ctStanFit(fit)
  probs=c(.025,.5,.975)
  Ygen <- e$Ygen
  # Ygen[Ygen==99999] <- NA
  ctDensityList(x=list(fit$data$Y,Ygen),plot=T,main='All variables',lwd=2,legend = c('Observed','Model implied'),xlab='Value')
  
  y<-aaply(Ygen,c(2,3,4),quantile,na.rm=TRUE,probs=probs,.drop=FALSE)
  # y<-array(unlist(lapply(c(.025,.5,.975),function(p) ctCollapse(Ygen,1,quantile,probs=p,na.rm=TRUE))),dim = c(3,dim(Ygen)[3],dim(Ygen)[4]))
  y<-array(y,dim=dim(y)[-1])  
  # y<-aperm(y,c(2,3,1))
  dimnames(y) <- list(NULL,fit$ctstanmodel$manifestNames,paste0(probs*100,'%'))
  x=1:fit$data$ndatapoints
  if(wait) readline("Press [return] for next plot.")
  
  for(i in 1:dim(Ygen)[4]){
    notmissing <- !is.na(c(y[,i,1]))
    ctPlotArray(list(y=y[notmissing,i,,drop=FALSE],x=x[notmissing]),legend=FALSE,plotcontrol=list(xlab='Observation',main=dimnames(y)[[2]][i]))
    ocol <- rgb(0,0,.7,.7)
    points(x[notmissing],fit$data$Y[,i][notmissing],type='l',lwd=1,col=ocol)
    if(legend) legend('topright',c('Model implied','Observed'),text.col=c('red',ocol))
    if(i < dim(Ygen)[4])  if(wait) readline("Press [return] for next plot.")
  }
}
