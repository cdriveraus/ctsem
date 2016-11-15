#' ctDensity
#'
#' Wrapper for base R density function that removes outliers and computes 'reasonable' bandwidth and x and y limits.
#' Used for ctsem density plots.
#' 
#' @param x numeric vector to compute density for.
#' @examples
#' y <- ctDensity(exp(rnorm(80)))
#' plot(y$density,xlim=y$xlim,ylim=y$ylim)
#' 
#' #### Compare to base defaults:
#' par(mfrow=c(1,2))
#' y=exp(rnorm(10000))
#' ctdens<-ctDensity(y)
#' plot(ctdens$density, ylim=ctdens$ylim,xlim=ctdens$xlim)
#' plot(density(y))
#' @export

ctDensity<-function(x){
  xlims=quantile(x,probs=c(.02,.98))
  mid=mean(c(xlims[2],xlims[1]))
  xlims[1] = xlims[1] - (mid-xlims[1])
  xlims[2] = xlims[2] + (xlims[2]-mid)
  x=x[x>xlims[1] & x<xlims[2]]
  bw=(max(x)-min(x)) / length(x)^.4 *.3
  
  xlims=quantile(x,probs=c(.01,.99))
  mid=mean(c(xlims[2],xlims[1]))
  xlims[1] = xlims[1] - (mid-xlims[1])/8
  xlims[2] = xlims[2] + (xlims[2]-mid)/8
  
  out1<-density(x,bw=bw,n=5000)
  out3=c(0,max(out1$y)*1.1)
  return(list(density=out1,xlim=xlims,ylim=out3))
}

