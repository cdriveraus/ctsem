
#' Compares model implied density and values to observed, for a ctStanFit object.
#'
#' @param fit ctStanFit object.
#' @param legend Logical, whether to plot a legend.
#' @param diffsize Integer > 0. Number of discrete time lags to use for data viz.
#' @param samples Logical -- plot all generated data points? Otherwise, plot shaded polygon based on quantile estimate. 
#' @param probs Vector of length 3 containing quantiles to plot -- should be rising numeric values between 0 and 1. 
#' @param wait Logical, if TRUE, waits for input before plotting next plot.
#' @param jitter Positive numeric between 0 and 1, if TRUE, jitters empirical data by specified proportion of std dev.
#' @param datarows integer vector specifying rows of data to plot. Otherwise 'all' uses all data.
#' @param ... extra arguments to pass to plot function.
#' @return If plot=FALSE, an array containing quantiles of generated data. If plot=TRUE, nothing, only plots.
#' @export
#' @details This function relies on the data generated during each iteration of fitting to approximate the
#' model implied distributions -- thus, when limited iterations are available, the approximation will be worse.
#'
#' @examples
#' ctStanPostPredict(ctstantestfit,wait=FALSE, samples=FALSE, datarows=1:10)
ctStanPostPredict <- function(fit,legend=TRUE,diffsize=1,jitter=.02, wait=TRUE,probs=c(.025,.5,.975),samples=TRUE, datarows='all',...){
  e<-extract.ctStanFit(fit)
  
  if(datarows[1]=='all') datarows <- 1:nrow(fit$data$Y)
  
  
  
  xmeasure=1:fit$data$ndatapoints
  xtime=fit$data$time
  
  Ygen <- array(e$Ygen,dim=c(dim(e$Ygen)[1],1,dim(e$Ygen)[-1]))
  
  ctDensityList(x=list(fit$data$Y[datarows,,drop=FALSE],Ygen[,,datarows,,drop=FALSE]),plot=TRUE,
    main='All variables',lwd=2,legend = c('Observed','Model implied'),xlab='Value',...)
  
  y<-aaply(Ygen,c(2,3,4),quantile,na.rm=TRUE,probs=probs,.drop=FALSE)
  y<-array(y,dim=dim(y)[-1])  
  dimnames(y) <- list(NULL,fit$ctstanmodel$manifestNames,paste0(probs*100,'%'))
  
  for(typei in c('obs','change')){
    for(i in 1:dim(Ygen)[4]){
      if(wait) readline("Press [return] for next plot.")
      xsmeasure=rep(xmeasure,each=dim(Ygen)[1])
      xstime=rep(xtime,each=dim(Ygen)[1])
      ycut=quantile(Ygen,c(.005,.995),na.rm=TRUE)
      ys=Ygen
      xsmeasure=xsmeasure[ys>ycut[1] & ys<ycut[2]]
      xstime=xstime[ys>ycut[1] & ys<ycut[2]]
      ys[ys<ycut[1] | ys>ycut[2]] <- NA
      
      
      
      if(typei=='obs'){
        
        ctDensityList(x=list(fit$data$Y[datarows,i,drop=FALSE],Ygen[,,datarows,i,drop=FALSE]),plot=TRUE,
          main=fit$ctstanmodel$manifestNames[i],lwd=2,legend = c('Observed','Model implied'),xlab='Value',...)
        
        if(wait) readline("Press [return] for next plot.")
        
        for(subtypei in c('Time','Observation')){
          if(subtypei=='Observation') x <- xmeasure
          if(subtypei=='Time') x <- xtime
          
          notmissing <- which(!is.na(c(y[datarows,i,1])))
          
          if(samples) {
            
            xs=rep(x,each=dim(Ygen)[1])
            ycut=quantile(Ygen[,,,i],c(.005,.995),na.rm=TRUE)
            ysamps=Ygen[,,,i]
            xs=xs[ysamps>ycut[1] & ysamps<ycut[2]]
            ysamps=ysamps[ysamps>ycut[1] & ysamps<ycut[2]]
            graphics::smoothScatter(xs,ysamps,nbin=256,colramp=grDevices::colorRampPalette(colors=c(rgb(1,1,1,0),rgb(1,.4,.4,.3))),nrpoints=0,
              transformation=function(x) x,ylim=range(c(y[,i,],quantile(ysamps,probs = c(.01,.99),na.rm=TRUE)),na.rm=TRUE),
              xlab=subtypei,ylab=dimnames(y)[[2]][i])
          }
          
          if(subtypei=='Observation') ctPlotArray(list(y=y[notmissing,i,,drop=FALSE],x=x[notmissing]),legend=FALSE,add=samples,polygon=!samples,
            plotcontrol=list(xlab=subtypei,main=dimnames(y)[[2]][i],...))
          
          # if(subtypei=='Time')
          
          ocol <- rgb(0,0,.7,.7)
          points(x[notmissing],
            fit$data$Y[,i][notmissing] +  rnorm(length(fit$data$Y[,i][notmissing]),0, jitter * sd(fit$data$Y[,i][notmissing],na.rm=TRUE)),
            type=ifelse(subtypei=='Time','p','l'),lwd=2,lty=1,pch=17, col=ocol)
          if(legend) legend('topright',c('Model implied','Observed'),text.col=c('red',ocol))
          if(i < dim(Ygen)[4])  if(wait) readline("Press [return] for next plot.")
        }
      }
      
      
      if(typei=='change'){
        yp<-aperm(ys,c(3,4,1,2))
        
        for(cdiffsize in diffsize){
          diffindex <- c() 
          for(diffi in 1:cdiffsize){
            diffindex <- c(diffindex,which(as.logical(fit$data$T0check[-1]))-(diffi-1),
              fit$data$ndatapoints-(diffi-1))
          }
          dygen<-diff(yp[,i,,,drop=TRUE],lag = cdiffsize) #drop true set here if looking for problems!
          # yp[-1,i,,,drop=FALSE] - yp[-fit$data$ndatapoints,i,,,drop=FALSE]
          dygendt <- dygen / diff(fit$data$time,lag = cdiffsize)
          dygendt<-dygendt[-diffindex,,drop=FALSE]
          
          dydt<-diff(fit$data$Y[,i], lag = cdiffsize)/diff(fit$data$time,lag = cdiffsize)
          dydt <- dydt[-diffindex]
          dydt <- dydt +  rnorm(length(dydt),0, jitter * sd(fit$data$Y[,i],na.rm=TRUE)) #add jitter
          time <- fit$data$time[-diffindex]
          # smoothScatter(matrix(yp[-fit$data$ndatapoints,i,,,drop=FALSE],ncol=1),
          #   matrix(dygendt[,,,,drop=FALSE],ncol=1),
          #   nbin=512,colramp=colorRampPalette(colors=c('white',rgb(1,0,0,1))),nrpoints=0,
          #   # transformation=function(x) x,
          #   ylim=range(c(quantile(c(dygendt),probs = c(.01,.99),na.rm=TRUE)),na.rm=TRUE),
          #   xlab='Observation',ylab=dimnames(y)[[2]][i])
          samps<-sample(1:length(dygendt),size=50000,replace=TRUE)
          plot(matrix(yp[-diffindex,i,,,drop=FALSE][samps],ncol=1),
            matrix(dygendt[,,drop=FALSE][samps],ncol=1),
            ylab=paste0('dy/dt, diff=',cdiffsize),xlab='y', main=dimnames(y)[[2]][i],
            pch=16,cex=.2,col=rgb(1,0,0,.1),...)
          points( fit$data$Y[-diffindex,i],
            dydt,
            col=rgb(0,0,1,.5),pch=17,...)
          
          if(wait) readline("Press [return] for next plot.")  
          
          plot(
            rep(time,(dim(dygendt)[2]))[samps],
            matrix(dygendt[,,drop=FALSE],ncol=1)[samps],
            ylab=paste0('dy/dt, diff=',cdiffsize),xlab='time', main=dimnames(y)[[2]][i],
            pch=16,cex=.1,col=rgb(1,0,0,.3),...)
          points(time,
            dydt,
            col=rgb(0,0,1,.5),pch=17,...)
          
        }
      }
    }
  }
  
}
