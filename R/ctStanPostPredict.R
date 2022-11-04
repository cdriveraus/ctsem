
#' Compares model implied density and values to observed, for a ctStanFit object.
#'
#' @param fit ctStanFit object. 
#' @param diffsize Integer > 0. Number of discrete time lags to use for data viz.
#' @param probs Vector of length 3 containing quantiles to plot -- should be rising numeric values between 0 and 1. 
#' @param wait Logical, if TRUE and \code{plot=TRUE}, waits for input before plotting next plot.
#' @param jitter Positive numeric between 0 and 1, if TRUE, jitters empirical data by specified proportion of std dev.
#' @param datarows integer vector specifying rows of data to plot. Otherwise 'all' uses all data.
#' @param nsamples Number of datasets to generate for comparisons, if fit object does not contain generated
#' data already.
#' @param resolution Positive integer, the number of rows and columns to split plots into for shading.
#' @param plot logical. If FALSE, a list of ggplot objects is returned.
#' @return If plot=FALSE, an array containing quantiles of generated data. If plot=TRUE, nothing, only plots.
#' @export
#' @details This function relies on the data generated during each iteration of fitting to approximate the
#' model implied distributions -- thus, when limited iterations are available, the approximation will be worse.
#' @return if plot=TRUE, nothing is returned and plots are created. Otherwise, a list containing ggplot objects is returned 
#' and may be customized as desired.
#' @examples
#' if(w32chk()){
#'
#' ctStanPostPredict(ctstantestfit,wait=FALSE, diffsize=2,resolution=100)
#' }
ctStanPostPredict <- function(fit,diffsize=1,jitter=.02, wait=TRUE,probs=c(.025,.5,.975),
  datarows='all',nsamples=500,resolution=100,plot=TRUE){
  
  plots <-list()
  if(datarows[1]=='all') datarows <- 1:nrow(fit$data$Y)
  xmeasure=data.table(id=fit$standata$subject[datarows])
  if(1==99) id <- count <- Density <- param <- NULL
  xmeasure= xmeasure[, count := seq(.N), by = .(id)]$count
  
  if(is.null(fit$generated$Y) || dim(fit$generated$Y)[1] < nsamples){
    Ygen <- ctStanGenerateFromFit(fit,fullposterior=TRUE,nsamples=nsamples)$generated$Y
    } else Ygen <- fit$generated$Y
  
  # Ygen<-aperm(Ygen,c(2,1,3)) #previously necessary but fixed generator
  Ygen <- Ygen[,datarows,,drop=FALSE]
  time <- fit$standata$time[datarows]
  Ydat <- fit$data$Y[datarows,,drop=FALSE]
  
  
  dens <- ctDensityList(x=list(Observed=Ydat[datarows,,drop=FALSE],Model=Ygen[,,,drop=FALSE]),
    main='',xlab='All variables',colvec=c('blue','red'),plot=FALSE,probs=c(.01,.99))
  densdt <- data.table(x= dens$density[[1]]$x, Density=dens$density[[1]]$y,source='Observed',param='All Variables')
  densdt <- rbind(densdt, data.table(x= dens$density[[2]]$x, Density=dens$density[[2]]$y,source='Model',param='All Variables'))
  
  
  y<-aaply(Ygen,c(2,3),quantile,na.rm=TRUE,probs=probs,.drop=FALSE)
  # y<-array(y,dim=dim(y)[-1])  
  dimnames(y) <- list(NULL,fit$ctstanmodel$manifestNames,paste0(probs*100,'%'))
  
  ycut=quantile(Ygen,c(.005,.995),na.rm=TRUE)
  ys=Ygen
  ys[ys<ycut[1] | ys>ycut[2]] <- NA 
  
  for(i in 1:dim(Ygen)[3]){ #dim 3 is indicator, dim 2 is datarows, dim 1 iterations
    plots[[i]] <- list()
    for(typei in c('obs','change')){
      
      # xsmeasure=rep(xmeasure,each=dim(Ygen)[1])
      # xstime=rep(xtime,each=dim(Ygen)[1])
      # xsmeasure=xsmeasure[ys>ycut[1] & ys<ycut[2]]
      # xstime=xstime[ys>ycut[1] & ys<ycut[2]]
      # 
      
      if(typei=='obs'){
        dens <- ctDensityList(x=list(Observed=Ydat[,i,drop=FALSE],Model=Ygen[,datarows,i,drop=FALSE]),
          main=fit$ctstanmodel$manifestNames[i],
          probs=c(.01,.99),
          xlab=dimnames(y)[[2]][i],
          colvec=c('blue','red'),plot=FALSE)
        densdt <- rbind(densdt, data.table(x= dens$density[[1]]$x, 
          Density=dens$density[[1]]$y,source='Observed', param=dimnames(y)[[2]][i]))
        densdt <- rbind(densdt, data.table(x= dens$density[[2]]$x, 
          Density=dens$density[[2]]$y,source='Model', param=dimnames(y)[[2]][i]))
        
        for(subtypei in c('Time','Observation')){
          if(subtypei=='Observation') x <- xmeasure
          if(subtypei=='Time') x <- time
          
        
          
          plots[[i]] <- c(plots[[i]],
            list(
              plotdensity2dby2(x1 = rep(x,each=dim(Ygen)[1]),
                y1 = c(Ygen[,,i]), #(dygendt[,,drop=FALSE]),
                x2=x,
                y2=Ydat[,i],
                xlab = subtypei,ylab = dimnames(y)[[2]][i],title='',
                grouplab = c('Observed','Model'),colours=c('blue','red'),
                resolution=resolution)
            ))
          
        }
      }
      
      
      
      if(typei=='change'){
        
        yp<-aperm(matrix(ys[,,i,drop=TRUE],nrow=dim(ys)[1],ncol=dim(ys)[2]),c(2,1))#drop true set here if looking for problems!
        
        for(cdiffsize in diffsize){
          # diffindex <- c() 
          # for(diffi in 1:cdiffsize){
          #   diffindex <- c(diffindex,which(as.logical(fit$data$T0check[datarows][-1]))-(diffi-1))
          # }
          
          subdiff <- c( which(diff(fit$data$subject[datarows],cdiffsize) != 0), datarows[length(datarows):1][1:cdiffsize])
          
          dygen<-diff(yp,lag = cdiffsize) 
          # yp[-1,i,,,drop=FALSE] - yp[-fit$data$ndatapoints,i,,,drop=FALSE]
          dygendt <- dygen / matrix(diff(time,lag = cdiffsize),nrow=length(time)-cdiffsize,ncol=nsamples)
          dygendt<-dygendt[-subdiff,,drop=FALSE]
          
          # dydt<-diff(Ydat[,i], lag = cdiffsize)/diff(time,lag = cdiffsize)
          dydt <- diff(Ydat[,i],lag=cdiffsize)
          dydt <- (dydt/ diff(time,lag = cdiffsize))[-subdiff]
          dydt <- dydt +  rnorm(length(dydt),0, jitter * sd(Ydat[,i],na.rm=TRUE)) #add jitter
          xtime <- time[-subdiff]
          # smoothScatter(matrix(yp[-fit$data$ndatapoints,i,,,drop=FALSE],ncol=1),
          #   matrix(dygendt[,,,,drop=FALSE],ncol=1),
          #   nbin=512,colramp=colorRampPalette(colors=c('white',rgb(1,0,0,1))),nrpoints=0,
          #   # transformation=function(x) x,
          #   ylim=range(c(quantile(c(dygendt),probs = c(.01,.99),na.rm=TRUE)),na.rm=TRUE),
          #   xlab='Observation',ylab=dimnames(y)[[2]][i])
          samps<-sample(1:length(dygendt),size=min(50000,length(dygendt)),replace=FALSE)
          # plot(matrix(yp[-subdiff,,drop=FALSE][samps],ncol=1),
          #   matrix(dygendt[,,drop=FALSE][samps],ncol=1),
          #   ylab=paste0('dy/dt, diff=',cdiffsize),xlab='y', main=dimnames(y)[[2]][i],
          #   pch=16,cex=.2,col=rgb(1,0,0,.1),...)
          # points( Ydat[-subdiff,i],
          #   dydt,
          #   col=rgb(0,0,1,.5),pch=17,...)
          # if(length(Ydat[-subdiff,i]) > 500) datasample <- sample(1:length(Ydat[-subdiff,i]),500) else datasample <- 1:length(Ydat[-subdiff,i])
          
          plots[[i]] <- c(plots[[i]],
            list(
              plotdensity2dby2(x1 = c(yp[-subdiff,,drop=FALSE]),
                y1 = c(dygendt[,,drop=FALSE]),
                x2=c(Ydat[-subdiff,i]),
                y2=dydt,
                xlab = dimnames(y)[[2]][i],ylab = paste0('d (',dimnames(y)[[2]][i],') / dt'), title='',
                grouplab = c('Observed','Model'),colours=c('blue','red'),
                resolution=resolution)
            ))
          
          
          # pd <- data.table(Source='Model',
          #   y=c(matrix(yp[-subdiff,,drop=FALSE][samps],ncol=1)),
          #   dydt=c(matrix(dygendt[,,drop=FALSE][samps],ncol=1)))
          # pd <- rbind(pd,data.table(Source='Observed', 
          #   y=Ydat[-subdiff,i][datasample],
          #   dydt=dydt))
          # 
          # lims <- lapply(c('dydt','y'),function(b) sapply(c('Observed','Model'),function (a) quantile(pd[Source==a,..b],c(.01,.99),na.rm=TRUE)))
          # names(lims) <- c('dydt','y')
          # lims <- lapply(lims, function(x) c(min(x[1,]),max(x[2,])))
          # 
          # suppressWarnings(print(ggplot(data=pd,aes(y=dydt,x=y,shape=Source,colour=Source)) + 
          #     stat_density_2d(data = subset(pd, Source=='Model'),geom="raster", 
          #       aes(alpha=..density..,fill = ..density..),show.legend = FALSE, contour = FALSE) +
          #     stat_density_2d(data = subset(pd, Source=='Observed'),
          #       aes(alpha=..level..),linetype='dotted',show.legend = FALSE, contour = TRUE) +
          #     coord_cartesian(xlim = lims$y, ylim = lims$dydt) +
          #     geom_point(data = subset(pd, Source=='Observed'),show.legend = TRUE) +
          #     scale_fill_gradient (low = "#FFFFFF", high = "#FF0000",guide='none')+
          #     labs(x='y',y='dy/dt',title = dimnames(y)[[2]][i], color  = "Source", shape = "Source")+ 
          #     theme_minimal() + 
          #     scale_alpha(guide = 'none') +
          #     scale_colour_manual(values=c("Observed"="blue", "Model"="red"))+
          #     scale_shape_manual(values=c("Observed" = 20, "Model" = 19)) +
          #     theme(legend.title = element_blank())))
          
          plots[[i]] <- c(plots[[i]],
            list(
              plotdensity2dby2(
                x1 = rep(xtime,(dim(dygendt)[2]))[samps],
                y1 = matrix(dygendt[,,drop=FALSE],ncol=1)[samps],
                x2=xtime,
                y2=dydt,
                xlab = 'time',ylab = paste0('d (',dimnames(y)[[2]][i],') / dt'),title='',
                grouplab = c('Observed','Model'),colours=c('blue','red'),
                resolution=resolution)
            ))
          
          # plot(
          #   rep(xtime,(dim(dygendt)[2]))[samps],
          #   matrix(dygendt[,,drop=FALSE],ncol=1)[samps],
          #   ylab=paste0('dy/dt, diff=',cdiffsize),xlab='time', main=dimnames(y)[[2]][i],
          #   pch=16,cex=.1,col=rgb(1,0,0,.3),...)
          # points(xtime,
          #   dydt,
          #   col=rgb(0,0,1,.5),pch=17,...)
          
        }
      }
    }
  }
  
  plots[[length(plots)+1]]<-c(list(
    ggplot(densdt,aes(x=x,fill=source,ymax=Density,y=Density) )+
      geom_line(alpha=.3) +
      geom_ribbon(alpha=.4,ymin=0) +
      scale_fill_manual(values=c('red','blue')) +
      theme_minimal()+
      theme(legend.title = element_blank(),
        panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = .2),
        strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm"))) +
      facet_wrap(vars(param),scales = 'free')
  )) #add density plots to front
  
  
  if(plot) {
    if(!requireNamespace('gridExtra')) plots <- unlist(plots,recursive = FALSE)
    firstplot=TRUE
    lapply(plots,function(x){
      if(wait && !firstplot) readline("Press [return] for next plot.")
      firstplot <<- FALSE
      if(requireNamespace('gridExtra')){
        if(length(x) > 1){
          for(li in seq_along(x)[-2]){
            x[[li]] <- x[[li]]  + theme(legend.position = "none")
          }
        }
        suppressWarnings(gridExtra::grid.arrange(grobs=x))
      } else suppressWarnings(print(x))
    })
    return(NULL)
  } else return(invisible(plots))
  
}
