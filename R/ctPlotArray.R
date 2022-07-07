ctPlotArrayGG <- function(input,
  ggeval = paste0("ggplot(dt[as.numeric(dt$Quantile) < mean(qu),],aes(x=x,y=value,colour=variable)) +
    ylab(names(dimnames(input$y)[3])) +xlab(names(dimnames(input$y)[1])) +
    geom_ribbon(aes(ymin=value,ymax=max,fill=variable),alpha=.1,linetype=0) +
    geom_line(lwd=.2,linetype='dotted')+
    geom_line(data=dt[as.numeric(dt$Quantile) > mean(qu),],lwd=.2,linetype='dotted')+
    geom_line(data=dt[as.numeric(dt$Quantile) == mean(qu),],lwd=1)+
    theme_minimal() + 
    theme(legend.title = element_blank())")){

  names <- names(attributes(input$y)$dimnames)
  
  dt=unlist(apply(input$y[,,1:dim(input$y)[3],drop=FALSE],3,function(x) list(x)),recursive = FALSE)
  dt=data.table(plyr::ldply(dt,.id=names(dimnames(input$y)[3])))

  dt$x <- rep(input$x,dim(input$y)[3])
  dt$Quantile <- factor(match(dt[[names[3] ]],unique(dt[[names[3] ]])))
  dt <- melt(dt[,!names(dt) %in% names[3],with=FALSE],id.vars = c('Quantile','x'))
  
  dt$max <- -9999
  qu <- as.numeric(unique(dt$Quantile))
  qalpha <- abs(mean(qu) - qu) 
  dt<-data.frame(dt)
  for(qi in 1:(floor(dim(input$y)[3]/2))){
    dt[dt$Quantile==qi,'max'] <- c(dt[dt$Quantile==qi+(2*qi),'value'])
  }
  
  if(1==99) variable <- value <- x <- NULL
  
  g<-eval(parse(text=ggeval))
  
  return(invisible(g))
  
}



#' Plots three dimensional y values for quantile plots
#' 
#' 1st margin of $Y sets line values, 2nd sets variables, 3rd quantiles.
#'
#' @param input list containing 3 dimensional array to use for Y values, \code{$y}
#' and vector of corresponding x values \code{$x}.
#' @param grid Logical. Plot with a grid?
#' @param add Logical. If TRUE, plotting is overlayed on current plot, without creating new plot.
#' @param colvec color vector of same length as 2nd margin.
#' @param lwdvec lwd vector of same length as 2nd margin.
#' @param ltyvec lty vector of same length as 2nd margin.
#' @param typevec type vector of same length as 2nd margin.
#' @param plotcontrol list of arguments to pass to plot.
#' @param legend Logical. Draw a legend?
#' @param legendcontrol list of arguments to pass to \code{\link[graphics]{legend}}.
#' @param polygon Logical. Draw the uncertainty polygon?
#' @param polygonalpha Numeric, multiplier for alpha (transparency) of the 
#' uncertainty polygon.
#' @param polygoncontrol list of arguments to pass to \code{\link{ctPoly}}
#'
#' @return Nothing. Generates plots.
#' @export
#'
#' @examples
#' \donttest{#'
#' input<-ctStanTIpredeffects(ctstantestfit, plot=FALSE, whichpars='CINT', 
#'  nsamples=10,nsubjects=10)
#'     
#' ctPlotArray(input=input)
#' }
ctPlotArray <- function(input,
  grid=FALSE,add=FALSE,
  colvec='auto',lwdvec='auto',ltyvec='auto',typevec='auto',
  plotcontrol=list(ylab='Array values',xaxs='i'),
  legend=TRUE,legendcontrol=list(),
  polygon=TRUE, polygonalpha=.1,polygoncontrol=list(steps=25)){
  if(!w32chk()) message('Bayesian functions not available on 32 bit systems') else{  
    if(!'list' %in% class(input)) stop('Input must be a list containing y and x subobjects!')
    
    x <- input$x
    if(length(dim(input$y))!=3) stop('$y subobject must have 3 dimensions!')
    if(length(x) != length(input$y[,1,1])) stop ('Length of 1st margin of $y and $x must be the same!')
    
    separate=FALSE
    nvars<-dim(input$y)[2]
    
    plotcontrolpars <- c('ylab','xlab')
    plotcontroldefaults <- c('Array values',colnames(input$x))
    for(i in 1:length(plotcontrolpars)) {
      if(is.null(plotcontrol[[plotcontrolpars[i]]])) plotcontrol[[plotcontrolpars[i]]] <- plotcontroldefaults[i]
    }
    
    legendcontrolpars <- c('x')
    legendcontroldefaults <- c('topright')
    for(i in 1:length(legendcontrolpars)) {
      if(is.null(legendcontrol[[legendcontrolpars[i]]])) legendcontrol[[legendcontrolpars[i]]] <- legendcontroldefaults[i]
    }
    
    if(is.null(legendcontrol$legend)){
      if(!is.null(dimnames(input$y)[2])) {
        legendcontrol$legend <- dimnames(input$y)[2] 
      } else  legend <- FALSE
    }
    
    # polygoncontrolpars <- c('border')
    # polygoncontroldefaults <- c(NA)
    # for(i in 1:length(polygoncontrolpars)) {
    #   if(is.null(polygoncontrol[[polygoncontrolpars[i]]])) polygoncontrol[[polygoncontrolpars[i]]] <- polygoncontroldefaults[i]
    # }
    
    if(all(ltyvec=='auto')) ltyvec = 1:nvars
    if(all(lwdvec=='auto')) lwdvec = rep(3,nvars)
    if(all(colvec=='auto')) colvec = rainbow(nvars,v=.9)
    if(all(typevec=='auto')) typevec = rep('l',nvars)
    # if(all(mainvec=='auto')){
    #   if(separate) mainvec=dimnames(input$y)[[2]] else mainvec =rep(ifelse(is.null(plotcontrol$main),'',plotcontrol$main),nvars)
    # }
    
    plotargs<-plotcontrol
    plotargs$x <- x
    if(!separate && is.null(plotcontrol$ylim)) plotargs$ylim = range(input$y,na.rm=TRUE)+ c(0,
      (max(input$y,na.rm=TRUE) - min(input$y,na.rm=TRUE)) /3)
    plotargs$xlim = range(x,na.rm=TRUE)
    
    ctpolyargs<-polygoncontrol
    
    legargs<-legendcontrol
    
    
    #blank plot
    if(add==FALSE){
      blankargs=plotargs
      blankargs$y=NA
      blankargs$x=NA
      do.call(plot,blankargs)
      if(grid) {
        grid()
        par(new=TRUE)
        do.call(plot,blankargs)
        par(new=FALSE)
      }
    }
    
    #confidence
    if(polygon) {
      for(pari in c(1:dim(input$y)[2],dim(input$y)[2]:1)){
        ctpolyargs$col=adjustcolor(colvec[pari],alpha.f=max(c(.004,polygonalpha/(4*sqrt(ctpolyargs$steps)))))
        ctpolyargs$x=plotargs$x[!is.na(plotargs$x) & !is.na(input$y[,pari,2])]
        ctpolyargs$y=input$y[,pari,2][!is.na(plotargs$x) & !is.na(input$y[,pari,2])] #predict(loess(input$y[,pari,2]~ctpolyargs$x))
        ctpolyargs$yhigh = input$y[,pari,3][!is.na(plotargs$x) & !is.na(input$y[,pari,2])] #predict(loess(input$y[,pari,3]~ctpolyargs$x))
        ctpolyargs$ylow = input$y[,pari,1][!is.na(plotargs$x) & !is.na(input$y[,pari,2])]#predict(loess(input$y[,pari,1]~ctpolyargs$x))
        
        do.call(ctPoly,ctpolyargs) 
      }
    }
    
    for(pari in c(1:dim(input$y)[2],dim(input$y)[2]:1)){
      for(qi in 1:3){
        plotargs$y = input$y[,pari,qi]#predict(loess(input$y[,pari,qi]~plotargs$x))
        plotargs$col=colvec[pari]
        plotargs$lwd=ifelse(qi==2,lwdvec[pari],1)
        plotargs$lty=ltyvec[pari]
        if(qi!=2) plotargs$col =grDevices::adjustcolor(plotargs$col,alpha.f=.5)
        plotargs$type=typevec[pari]
        # plotargs$main=mainvec[pari]
        do.call(points,plotargs)
      }
      
      if(separate && legend) {
        legargs$legend=dimnames(input$y)[[2]][pari]
        legargs$col = plotargs$col
        legargs$lty = plotargs$lty
        legargs$text.col = plotargs$col
        legargs$lwd = plotargs$lwd
        do.call(graphics::legend,legargs) 
      }
    }
    
    if(!separate & legend) {
      legargs$legend=dimnames(input$y)[[2]]
      legargs$col = colvec
      legargs$lty = ltyvec
      legargs$text.col = colvec
      legargs$lwd = lwdvec
      do.call(graphics::legend,legargs) 
    }
    
  } 
}
