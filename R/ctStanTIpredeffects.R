#' Get time independent predictor effect estimates
#' 
#' COmputes and plots combined effects and quantiles for effects of time independent predictors
#' on subject level parameters of a ctStanFit object.
#'
#' @param fit fit object from \code{\link{ctStanFit}}
#' @param probs numeric vector of quantile probabilities from 0 to 1. Specify 3
#' values if plotting, the 2nd will be drawn as a line with uncertainty polygon
#' based on 1st and 3rd.
#' @param whichTIpreds integer vector specifying which of the tipreds in the fit object you want to
#' use to calculate effects. Unless quadratic / higher order versions of predictors have been 
#' included, selecting more than one probably doesn't make sense. If for instance a squared
#' predictor has been included, then you can specify both the linear and squared version. 
#' The x axis of the plot (if generated) will be based off the first indexed predictor. To 
#' check what predictors are in the model, run \code{fit$ctstanmodel$TIpredNames}.
#' @param whichpars Integer vector specifying Which of the individually varying subject
#' level parameters to compute effects on. 
#' 'auto' uses all available, which is equivalent to 
#' \code{1:sum(fit$ctstanmodel$pars$indvarying)}.
#' The integer corresponding to specific parameters can be found as follows, replacing \code{fit} as appropriate:
#' \code{fit$ctstanmodel$pars[sf$ctstanmodel$pars$indvarying,'param']}.
#' @param plot Logical. If TRUE, nothing is returned but instead \code{\link{ctPlotArray}}
#' is used to plot the output instead.
#' @param ... arguments to pass to \code{\link{ctPlotArray}} for plotting.
#' @return Either a three dimensional array of predictor effects, or nothing with a plot
#' generated.
#' @export
#'
#' @examples
#' ctStanTIpredeffects(ctstantestfit,plot=TRUE)
#' ctStanTIpredeffects(ctstantestfit,plot=TRUE,separate=TRUE)
ctStanTIpredeffects<-function(fit,probs=c(.25,.5,.75),
  whichTIpreds=1,whichpars='all',
  plot=FALSE,...){
  
  if(all(whichpars=='all')) whichpars=1:sum(fit$ctstanmodel$pars$indvarying)
  
e<-extract(fit$stanfit)
hypermeans <- e$hypermeans[,
  fit$ctstanmodel$pars$indvarying[is.na(fit$ctstanmodel$pars$value)],
  drop=FALSE] #first get indvarying hypermeans
hypermeans <- hypermeans[,whichpars,drop=FALSE]  #then just the ones in whichpars
tieffect<-e$tipredeffect[,whichpars,whichTIpreds,drop=FALSE]
tipreds<-fit$data$tipreds[,whichTIpreds,drop=FALSE]
tiorder<-order(tipreds[,1])
tipreds<-tipreds[tiorder,,drop=FALSE] #order tipreds according to first one
npars<-length(whichpars)
niter<-dim(e$hypermeans)[1]
nsubjects <- nrow(tipreds)
spec <- fit$ctstanmodel$pars[fit$ctstanmodel$pars$indvarying,,drop=FALSE]
spec <- spec[whichpars,,drop=FALSE]

message('Calculating time independent predictor effects...')
  
  raweffect <- aaply(1:niter,1,function(iterx) { #for every iter
    aaply(tipreds,1,function(tix){ #and every distinct tipred vector
      hypermeans[iterx,] + matrix(tieffect[iterx,,],nrow=dim(tieffect)[2]) %*% tix
    },.drop=FALSE)
  })

  
effect<-aaply(1:npars, 1,function(pari){ #for each param
param=raweffect[,,pari]
  out=eval(parse(text=spec$transform[pari]))
  return(out)
})

noeffect<-aaply(1:npars, 1,function(pari){ #for each param
  param <- hypermeans[,pari]
  out=eval(parse(text=spec$transform[pari]))
  return(out)
})

deffect<-effect-array(noeffect,dim=dim(effect))

out<-aaply(probs,1,function(x) ctCollapse(deffect,2,quantile,probs=x))
dimnames(out)=list(Quantile=paste0('Quantile',probs),
  popmean=spec$param,
  subject=tiorder #subjects reordered because tipreds were at top
)

out <- aperm(out, c(3,2,1))

if(!plot) return(out) else {
  ctPlotArray(yarray=out,x=tipreds[,1],
  plotcontrol=list(ylab='Effect',xlab=colnames(tipreds)[1]))
}
}





#' Plots a three dimensional array
#' 
#' 1st margin sets line values, 2nd sets variables, 3rd quantiles.
#'
#' @param yarray 3 dimensional array to use for Y values
#' @param x numeric vector specifying x axis
#' @param grid Logical. Plot with a grid?
#' @param separate Logical. Generate a plot per variable (2nd margin of array)
#' or plot all on one?
#' @param colvec color vector of same length as 2nd margin.
#' @param lwdvec lwd vector of same length as 2nd margin.
#' @param ltyvec lty vector of same length as 2nd margin.
#' @param typevec type vector of same length as 2nd margin.
#' @param mainvec main vector of same length as 2nd margin, only used if 
#' separate=TRUE.
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
#' y <- ctStanTIpredeffects(ctstantestfit,plot=FALSE)
#' x<-ctstantestfit$data$tipreds[order(ctstantestfit$data$tipreds[,1]),1]
#' ctPlotArray(y,x,separate=TRUE)
ctPlotArray <- function(yarray,x,
grid=TRUE,separate=FALSE,
  colvec='auto',lwdvec='auto',ltyvec='auto',typevec='auto',
  mainvec='auto',
  plotcontrol=list(ylab='Array values', xlab='X values',xaxs='i'),
legend=TRUE,legendcontrol=list(x='topright'),
polygon=TRUE, polygonalpha=.1,polygoncontrol=list(border=NA,steps=50)){

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
  
  polygoncontrolpars <- c('border')
  polygoncontroldefaults <- c(NA)
  for(i in 1:length(polygoncontrolpars)) {
    if(is.null(polygoncontrol[[polygoncontrolpars[i]]])) polygoncontrol[[polygoncontrolpars[i]]] <- polygoncontroldefaults[i]
  }
  
  if(all(ltyvec=='auto')) ltyvec = 1:nvars
  if(all(lwdvec=='auto')) lwdvec = rep(3,nvars)
  if(all(colvec=='auto')) colvec = rainbow(nvars)
  if(all(typevec=='auto')) typevec = rep('l',nvars)
  if(all(mainvec=='auto')){
    if(separate) mainvec=dimnames(yarray)[[2]] else mainvec =rep(ifelse(is.null(plotcontrol$main),'',plotcontrol$main),nvars)
  }

  plotargs<-plotcontrol
  plotargs$x <- x
  if(!separate && is.null(plotcontrol$ylim)) plotargs$ylim = range(yarray)
   plotargs$xlim = range(x,na.rm=TRUE)
   
   ctpolyargs<-polygoncontrol
   
   legargs<-legendcontrol
  
   
   #blank plot
   blankargs=plotargs
   blankargs$y=NA
   blankargs$x=NA
   do.call(plot,blankargs)
   if(grid) grid()
   
   #confidence
   if(polygon) {
     for(pari in c(1:dim(yarray)[2],dim(yarray)[2]:1)){
     ctpolyargs$col=adjustcolor(colvec[pari],alpha.f=max(c(.004,polygonalpha/ctpolyargs$steps)))
     ctpolyargs$x=plotargs$x
     ctpolyargs$y=yarray[,pari,2]
     ctpolyargs$yhigh = yarray[,pari,3]
     ctpolyargs$ylow = yarray[,pari,1]
     
     do.call(ctPoly,ctpolyargs) 
     }
   }
   
for(pari in c(1:dim(yarray)[2],dim(yarray)[2]:1)){
  for(qi in 1:3){
  plotargs$y = yarray[,pari,qi]
  plotargs$col=colvec[pari]
  plotargs$lwd=ifelse(qi==2,lwdvec[pari],1)
  plotargs$lty=ltyvec[pari]
  if(qi!=2) plotargs$col =grDevices::adjustcolor(plotargs$col,alpha.f=.5)
  plotargs$type=typevec[pari]
  plotargs$main=mainvec[pari]
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

