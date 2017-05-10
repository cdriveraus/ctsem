#' Get time independent predictor effect estimates
#' 
#' Computes and plots combined effects and quantiles for effects of time independent predictors
#' on subject level parameters of a ctStanFit object.
#'
#' @param fit fit object from \code{\link{ctStanFit}}
#' @param returndifference logical. If FALSE, absolute parameter values are returned. 
#' If TRUE, only the effect of the covariate (i.e. without the average value of the parameter)
#' are returned. The former can be easier to interpret, but the latter are more likely to fit on 
#' a single plot.
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
ctStanTIpredeffects<-function(fit,returndifference=FALSE, probs=c(.025,.5,.975),
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

if(returndifference) effect<-effect-array(noeffect,dim=dim(effect))

out<-aaply(probs,1,function(x) ctCollapse(effect,2,quantile,probs=x),.drop=FALSE)
dimnames(out)=list(Quantile=paste0('Quantile',probs),
  popmean=spec$param,
  subject=tiorder #subjects reordered because tipreds were at top
)

out <- aperm(out, c(3,2,1))

if(!plot) return(out) else {
  ctPlotArray(yarray=out,x=tipreds[,1],
  plotcontrol=list(ylab='Effect',xlab=colnames(tipreds)[1]),...)
}
}

