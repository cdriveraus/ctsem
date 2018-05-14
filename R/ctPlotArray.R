
#' Plots three dimensional y values for quantile plots
#' 
#' 1st margin of $Y sets line values, 2nd sets variables, 3rd quantiles.
#'
#' @param input list containing 3 dimensional array to use for Y values, \code{$y}
#' and vector of corresponding x values \code{$x}.
#' @param grid Logical. Plot with a grid?
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
#' y <- ctStanTIpredeffects(ctstantestfit,plot=FALSE,whichpars='dtDRIFT',nsamples=10)
#' x<-ctstantestfit$data$tipreds[order(ctstantestfit$data$tipreds[,1]),1]
#' ctPlotArray(y,x)
ctPlotArray <- function(input,
  grid=FALSE,
  colvec='auto',lwdvec='auto',ltyvec='auto',typevec='auto',
  plotcontrol=list(ylab='Array values', xlab='X values',xaxs='i'),
  legend=TRUE,legendcontrol=list(x='topright'),
  polygon=TRUE, polygonalpha=.1,polygoncontrol=list(steps=25)){
  
  if(class(input)!='list') stop('Input must be a list containing y and x subobjects!')
  yarray <- input$y
  x <- input$x
  if(length(dim(yarray))!=3) stop('$y subobject must have 3 dimensions!')
  if(length(x) != length(yarray[,1,1])) stop ('Length of 1st margin of $y and $x must be the same!')
  
  separate=FALSE
  nvars<-dim(yarray)[2]
  
  plotcontrolpars <- c('ylab','xlab')
  plotcontroldefaults <- c('Array values','X values')
  for(i in 1:length(plotcontrolpars)) {
    if(is.null(plotcontrol[[plotcontrolpars[i]]])) plotcontrol[[plotcontrolpars[i]]] <- plotcontroldefaults[i]
  }
  
  legendcontrolpars <- c('x')
  legendcontroldefaults <- c('topright')
  for(i in 1:length(legendcontrolpars)) {
    if(is.null(legendcontrol[[legendcontrolpars[i]]])) legendcontrol[[legendcontrolpars[i]]] <- legendcontroldefaults[i]
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
  #   if(separate) mainvec=dimnames(yarray)[[2]] else mainvec =rep(ifelse(is.null(plotcontrol$main),'',plotcontrol$main),nvars)
  # }
  
  plotargs<-plotcontrol
  plotargs$x <- x
  if(!separate && is.null(plotcontrol$ylim)) plotargs$ylim = range(yarray,na.rm=TRUE)+ c(0,
    (max(yarray,na.rm=TRUE) - min(yarray,na.rm=TRUE)) /3)
  plotargs$xlim = range(x,na.rm=TRUE)
  
  ctpolyargs<-polygoncontrol
  
  legargs<-legendcontrol
  
  
  #blank plot
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
  
  #confidence
  if(polygon) {
    for(pari in c(1:dim(yarray)[2],dim(yarray)[2]:1)){
      ctpolyargs$col=adjustcolor(colvec[pari],alpha.f=max(c(.004,polygonalpha/(4*sqrt(ctpolyargs$steps)))))
      ctpolyargs$x=plotargs$x[!is.na(plotargs$x) & !is.na(yarray[,pari,2])]
      ctpolyargs$y=yarray[,pari,2][!is.na(plotargs$x) & !is.na(yarray[,pari,2])] #predict(loess(yarray[,pari,2]~ctpolyargs$x))
      ctpolyargs$yhigh = yarray[,pari,3][!is.na(plotargs$x) & !is.na(yarray[,pari,2])] #predict(loess(yarray[,pari,3]~ctpolyargs$x))
      ctpolyargs$ylow = yarray[,pari,1][!is.na(plotargs$x) & !is.na(yarray[,pari,2])]#predict(loess(yarray[,pari,1]~ctpolyargs$x))
      
      do.call(ctPoly,ctpolyargs) 
    }
  }
  
  for(pari in c(1:dim(yarray)[2],dim(yarray)[2]:1)){
    for(qi in 1:3){
      plotargs$y = yarray[,pari,qi]#predict(loess(yarray[,pari,qi]~plotargs$x))
      plotargs$col=colvec[pari]
      plotargs$lwd=ifelse(qi==2,lwdvec[pari],1)
      plotargs$lty=ltyvec[pari]
      if(qi!=2) plotargs$col =grDevices::adjustcolor(plotargs$col,alpha.f=.5)
      plotargs$type=typevec[pari]
      # plotargs$main=mainvec[pari]
      do.call(points,plotargs)
    }
    
    if(separate & legend) {
      legargs$legend=dimnames(yarray)[[2]][pari]
      legargs$col = plotargs$col
      legargs$lty = plotargs$lty
      legargs$text.col = plotargs$col
      legargs$lwd = plotargs$lwd
      do.call(graphics::legend,legargs) 
    }
  }
  
  if(!separate & legend) {
    legargs$legend=dimnames(yarray)[[2]]
    legargs$col = colvec
    legargs$lty = ltyvec
    legargs$text.col = colvec
    legargs$lwd = lwdvec
    do.call(graphics::legend,legargs) 
  }
  
  
}
