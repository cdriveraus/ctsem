#' ctIndplot
#' 
#' Convenience function to simply plot individual subject trajectories from ctsem wide format data
#' @param datawide ctsem wide format data
#' @param n.subjects Number of subjects to randomly select for plotting.
#' @param Tpoints Number of discrete time points per case in data structure
#' @param n.manifest Number of manifest variables in data structure
#' @param colourby set  plot colours by "subject" or "variable"
#' @param vars either 'all' or a numeric vector specifying which manifest variables to plot.
#' @param opacity Opacity of plot lines
#' @param varnames vector of variable names for legend (defaults to NULL)
#' @param xlab X axis label.
#' @param ylab Y axis label.
#' @param start Measurement occasion to start plotting from - defaults to T0.
#' @param legendposition Where to position the legend.
#' @param ... additional plotting parameters.
#' @examples 
#' 
#' data(ctExample1)
#' ctIndplot(ctExample1,n.subjects=1, n.manifest=2,Tpoints=6, colourby='variable')
#' 
#' @export
ctIndplot<-function(datawide,n.subjects,n.manifest,Tpoints,colourby="subject",
  vars='all',opacity=1,varnames=NULL,xlab='Time',ylab='Value',start=0,
  legendposition='topright',...){

  
  
  if(length(vars)==1 && vars=='all') vars<-1:n.manifest
  
  if(colourby=="variable") colourvector <- grDevices::rainbow(length(vars),alpha=opacity)
  if(colourby=="subject") colourvector <- grDevices::rainbow(n.subjects,alpha=opacity)
 
  
  ymin<-min(datawide[1:nrow(datawide),cseq(vars,n.manifest*Tpoints,n.manifest)],na.rm=T)
  ymax<-max(datawide[1:nrow(datawide),cseq(vars,n.manifest*Tpoints,n.manifest)],na.rm=T)
  
  # browser()
  individuals<-sample(1:nrow(datawide),n.subjects)
  times<-matrix(unlist(lapply(1:(Tpoints-1),function(x){
    apply(datawide[individuals,,drop=FALSE][,paste0('dT',1:x),drop=FALSE],1,sum,na.rm=T)
  })),ncol=(Tpoints-1))
  
  graphics::plot(0:0,ylim=c(ymin,ymax),type="l",xlim=c(start,max(times)),
    ylab=ylab,xlab=xlab,...)
  
 
  message(c('Plotting rows ',paste0(individuals,", ")))
  for(i in 1:n.subjects){

    for(j in 1:length(vars)){
      graphics::points(c(0,times[i,]),
        datawide[individuals[i],seq(vars[j],n.manifest*Tpoints,n.manifest)],
        col=ifelse(colourby=="variable",colourvector[j],colourvector[i]),type="b",pch=j,lty=1,...) 
    }}
  
  if(is.null(varnames)) varnames <- paste0("Y",1:n.manifest) #set varnames for legend
  
  if(colourby=="variable") {
    graphics::legend(legendposition,varnames,pch=vars,col=colourvector,text.col=colourvector,bty="n")
  }
  if(colourby=="subject") {
    graphics::legend(legendposition,varnames,pch=vars,bty="n")
  }
}
