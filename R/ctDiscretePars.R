#' ctStanParnames
#' 
#' Gets internal stan parameter names of a ctStanFit object sampled via stan based on specified substrings.
#'
#' @param x ctStanFit object
#' @param substrings vector of character strings, parameter names of the stan model
#' containing any of these strings will be returned. Useful strings may be 'pop_' for 
#' population means, 'popsd' for population standard deviations,
#'  or specific combinations such as 'pop_DRIFT' for the population
#' means of temporal dynamics parameters
#' @return vector of character strings.
#' @examples
#' \donttest{
#' sunspots<-sunspot.year
#' sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#' id <- 1
#' time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#'
#' #setup model
#' ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
#'  manifestNames='sunspots', 
#'  latentNames=c('ss_level', 'ss_velocity'),
#'  LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
#'  DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
#'  MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
#'  MANIFESTVAR=diag(0,1), #As per original spec
#'  CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'  DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#'
#' #fit
#' ssfit <- ctStanFit(datalong, ssmodel, iter=2, 
#'   optimize=FALSE, chains=1)
#' ctStanParnames(ssfit,substrings=c('pop_','popsd'))
#' }
#' 
#' @export
ctStanParnames <- function(x,substrings=c('pop_','popsd')){
  if(!'stanfit' %in% class(x$stanfit)) stop('Doesnt contain sampled stanfit object')
  out<-c()
  for(subsi in substrings){
    out<- c(out, x$stanfit@model_pars[grep(paste0('^',subsi),x$stanfit@model_pars)])
  }
  return(out)
}




#'ctDiscretePars
#'
#'Generate discrete time parameters for a sequence of times based on a list containing coninuous time
#'parameter matrices as used in ctsem.
#'
#'@param ctpars List of continuous time parameter matrices.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param type String. 'all' returns all possible outputs in a list. 
#''discreteDRIFT' returns only discrete time auto and cross regression effects.
#''latentMeans' returns only the expected latent means, given initial (T0MEANS) level, 
#'latent intercept (CINT) and temporal effects (DRIFT).
#'@examples
#'\donttest{
#'pars <- ctStanContinuousPars(ctstantestfit)
#'ctDiscretePars(pars,times=c(.5,1))
#'}
#'
#'@export
ctDiscretePars<-function(ctpars,times=seq(0,10,.1),type='all'){
  
  if(type=='all') type=c('discreteDRIFT','latentMeans') #this needs to match with ctStanDiscretePars
  nlatent=nrow(ctpars$DRIFT)
  latentNames=paste0('eta',1:nlatent)
  
  out<-list()
  discreteDRIFT = array(unlist(lapply(times, function(x) 
    expm::expm(ctpars$DRIFT*x))),
    dim=c(nlatent,nlatent,length(times)),
    dimnames=list(latentNames,latentNames,paste0('t',times)))
  
  if('discreteDRIFT' %in% type) out$discreteDRIFT = discreteDRIFT
  
  if('latentMeans' %in% type) out$latentMeans = array(unlist(lapply(1:length(times), function(x)
    discreteDRIFT[,,x] %*% ctpars$T0MEANS)),
    dim=c(nlatent,length(times)),
    dimnames=list(latentNames,paste0('t',times)))
  
  # if('impulse' %in% type) out$impulse <- 
  
  return(out)
}


#'ctStanDiscretePars
#'
#'Calculate model implied regressions for a sequence of time intervals based on a continuous time model fit
#'from ctStanFit, for specified subjects.
#'
#'@param ctstanfitobj Continuous time model fit from \code{\link{ctStanFit}}
#'@param subjects Either 'all', to take the average over all subjects, or a vector of integers denoting which
#'subjects.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param quantiles Which quantiles to return. If plotting, specify 3 quantiles, 
#'the 2nd will be plotted as a line with 1 and 3 as uncertainty bounds.
#'@param nsamples Number of samples from the stanfit to use for plotting. Higher values will
#'increase smoothness / accuracy, at cost of plotting speed. Values greater than the total
#'number of samples will be set to total samples.
#'@param observational Logical. If TRUE, outputs expected change in processes *conditional on observing* a 1 unit change in each -- 
#'this change is correlated according to the DIFFUSION matrix. If FALSE, outputs expected regression values -- also interpretable as
#'an independent 1 unit change on each process, giving the expected response under a 1 unit experimental impulse.
#'@param standardise Logical. If TRUE, output is standardised according to expected total within subject variance, given by the 
#'asymDIFFUSION matrix.
#'@param cov Logical. If TRUE, covariances are returned instead of regression coefficients.
#'@param plot Logical. If TRUE, plots output using \code{\link{ctStanDiscreteParsPlot}}
#'instead of returning output. 
#'@param ... additional plotting arguments to control \code{\link{ctStanDiscreteParsPlot}}
#'@examples
#'if(w32chk()){
#'
#' ctStanDiscretePars(ctstantestfit,times=seq(.5,4,.1), 
#'  plot=TRUE,indices='all')
#'  
#'#modify plot
#'require(ggplot2)
#'g=ctStanDiscretePars(ctstantestfit,times=seq(.5,4,.1), 
#'  plot=TRUE,indices='CR')
#'g= g+ labs(title='Cross effects')
#'print(g)
#'
#'}
#'@export
ctStanDiscretePars<-function(ctstanfitobj, subjects='all', times=seq(from=0,to=10,by=.1), 
  quantiles = c(.025, .5, .975),nsamples=500,observational=FALSE,standardise=FALSE, 
  cov=FALSE, plot=FALSE,...){

  type='discreteDRIFT'
  collapseSubjects=TRUE #consider this for a switch
  e<-ctExtract(ctstanfitobj,subjectMatrices = subjects[1]!='all')
  
  # if(type=='all') type=c('discreteDRIFT','latentMeans') #must match with ctDiscretePars
  
  if(subjects[1] != 'all' && any(!is.integer(as.integer(subjects)))) stop('
  subjects argument must be either "all" or an integer denoting specific subjects')
  
  nsubjects <- dim(e$subj_DRIFT)[2]
  if(is.null(nsubjects)) nsubjects=1
  if('all' %in% subjects) subjects='all' 
  
  if(is.null(e$subj_DRIFT) && any(!subjects %in% 'all')) stop('No individual variation in DRIFT matrix found?? Try subjects="all"')
  
  niter=dim(e$pop_DRIFT)[1]
  nlatent=ctstanfitobj$standata$nlatent#outdims[3]
  latentNames=ctstanfitobj$ctstanmodel$latentNames
  
  if(nsamples > niter) nsamples <- niter
  
  out<-list()
  
  #get all ctparameter matrices at once and remove unneeded subjects
  ctpars <- list()
  # browser()
  for(matname in c('DRIFT','DIFFUSIONcov','asymDIFFUSION')){ #,'CINT','T0MEANS', 'T0VAR','MANIFESTMEANS',if(!is.null(e$MANIFESTVAR)) 'MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    if('all' %in% subjects || is.null(e[[paste0('subj_',matname)]])){
      ctpars[[matname]] <- e[[paste0('pop_',matname)]]
    } else {
      # browser()
      ctpars[[matname]] <- e[[paste0('subj_',matname)]][,subjects,,,drop=FALSE]
      ctpars[[matname]]<-  array(ctpars[[matname]],dim=c(prod(dim(ctpars[[matname]])[1:2]),dim(ctpars[[matname]])[-1:-2]))
    }
    ctpars[[matname]] <- ctpars[[matname]][sample(1:dim(ctpars[[matname]])[1],nsamples),,,drop=FALSE]
  }
  
  # 
  nlatent <- dim(ctpars$asymDIFFUSION)[2]
  ctpars$DRIFT <- ctpars$DRIFT[,1:nlatent,1:nlatent,drop=FALSE] #intoverpop

  
  out <- ctStanDiscreteParsDrift(ctpars,times, observational, standardise, cov=cov)
  out <- apply(out,c(1,2,3),quantile,probs=quantiles)
    
    dimnames(out)<- list(quantiles=paste0('quantile_',quantiles),
      row=latentNames,
      col=latentNames,
      times=paste0('t',times)
    )
    
    out=aperm(out,c(2,3,4,1))
  
  if(plot) {
    
    out <- ctStanDiscreteParsPlot(out,times=times,latentNames=ctstanfitobj$ctstanmodel$latentNames,...)
  } 
  return(out)
}



ctStanDiscreteParsDrift<-function(ctpars,times, observational,  standardise,cov=FALSE){
  nl=dim(ctpars$DRIFT)[3]
  message('Computing temporal regression coefficients for ', dim(ctpars$DRIFT)[1],' samples, may take a moment...')
  discreteDRIFT <- array(sapply(1:(dim(ctpars$DRIFT)[1]),function(d){
    if(observational|standardise){
      asymDIFFUSIONdiag <- diag(matrix(ctpars$asymDIFFUSION[d,,],nl,nl))
    asymDIFFUSIONdiag[rl(asymDIFFUSIONdiag <= 0) ] <- 1
    }
    DRIFT <- matrix(ctpars$DRIFT[d,,],nl,nl)
    if(observational) {
      g <- matrix(ctpars$DIFFUSIONcov[d,,],nl,nl)
      g <- cov2cor(g)^2 * sign(g)
      g[is.nan(g)] <- 0
    }
    sapply(times, function(ti){ 
      out <-expm::expm(DRIFT *ti)
      if(standardise) out <- out * matrix(rep(sqrt((asymDIFFUSIONdiag)),each=nl) / 
          rep((sqrt(asymDIFFUSIONdiag)),times=nl),nl)
      if(observational) out <- out %*% g
      # browser()
      if(cov) out <- (out %*% t(out))
      return(matrix(out,ncol(out),ncol(out)))
    },simplify = 'array')
  },simplify = 'array'),dim=c(nl,nl,length(times),dim(ctpars$DRIFT)[1]))
}

#'ctStanDiscreteParsPlot
#'
#'Plots model implied regression strengths at specified times for 
#'continuous time models fit with ctStanFit.
#'
#'@param x list object returned from \code{\link{ctStanDiscretePars}}.
#'@param indices Either a string specifying type of plot to create, or an n by 2
#'matrix specifying which indices of the output matrix to plot.
#''AR' specifies all diagonals, for discrete time autoregression parameters.
#''CR' specifies all off-diagonals,for discrete time cross regression parameters.
#''all' plots all AR and CR effects at once.
#'@param add Logical. If FALSE, a new plot is generated, if TRUE, specified plot/s are
#'overlayed on existing plot.
#'@param gg Logical -- use GGplot2 or not? if TRUE, other graphical parameters are ignored, and the
#'ggplot object is returned and may be modified further. 
#'@param plot Logical. Only relevant with gg=TRUE. 
#'@param legend Logical. If TRUE, generates a legend.
#'@param grid Logical. Plot with a grid?
#'@param polygon Logical. If TRUE, fills a polygon between the first and last specified quantiles.
#'@param quantiles numeric vector of length 3, with values between 0 and 1, specifying which quantiles to plot.
#'The default of c(.05,.5,.95) plots 95\% credible intervals and the posterior median at 50\%. 
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param latentNames Vector of character strings denoting names for the latent variables. 
#''auto' just uses eta1 eta2 etc.
#'@param lwdvec Either 'auto', or a vector of positive integers denoting line widths for each quantile.
#''auto' specifies c(1,3,1) if there are 3 quantiles to be plotted (default), otherwise simply 3.
#'@param ltyvec Either 'auto', or a vector of line type integers (as for the lty parameter normally)
#' denoting line types for each quantile.
#' 'auto' specifies c(3, 1, 3) if there are 3 quantiles to be plotted (default), otherwise simply 1.
#'@param colvec Either 'auto', or a vector of color values denoting colors for each index to be plotted.
#''auto' generates colors using the \code{grDevices::rainbow} function.
#'@param plotcontrol list of arguments to pass to plot function. 
#'The following arguments are ignored: ylim,lwd,lty,col,x,y.
#'@param legendcontrol list of arguments to pass to legend function. 'legend=' and 'text.col=' arguments
#'will be ignored.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha (transparency) of colvec by for 
#'the fill polygon.
#'@param polygoncontrol list of arguments to pass to ctPoly function (if polygon=TRUE).
#'x,y, and col arguments will be ignored. Steps specifies the number of polygons to overlay to 
#'create a graduated transparency. Set to 1 for a flat looking plot.
#'@param ... for plot adjustments a ggeval argument can be added, which should be based on the default code
#'found in the ctsem:::ctPlotArrayGG function.
#'@examples
#'if(w32chk()){
#'x <- ctStanDiscretePars(ctstantestfit)
#'ctStanDiscreteParsPlot(x, indices='CR')
#'
#'#to modify plot:
#'g <- ctStanDiscreteParsPlot(x, indices='CR',plot=FALSE) + 
#'  ggplot2::labs(title='My ggplot modification')
#'print(g)
#'}
#'
#'@export

ctStanDiscreteParsPlot<- function(x,indices='all',add=FALSE,legend=TRUE, polygon=TRUE, 
  gg=TRUE,plot=TRUE,
  quantiles=c(.025,.5,.975), times=seq(0,10,.1),latentNames='auto',
  lwdvec='auto',colvec='auto',ltyvec='auto',
  plotcontrol=list(ylab='Value',xlab='Time interval',
    main='Regression coefficients',type='l', xaxs='i'),grid=FALSE,
  legendcontrol=list(x='topright',bg='white'),
  polygonalpha=.1,
  polygoncontrol=list(steps=20),...){

  if(is.data.frame(indices)) indices <- as.matrix(indices)
  
  nlatent=dim(x)[1]
  
  if(latentNames[1]=='auto') latentNames=dimnames(x)$row
  
  if(is.null(plotcontrol$ylab)) plotcontrol$ylab  <- 'Value'
  if(is.null(plotcontrol$xlab)) plotcontrol$xlab  <- 'Time interval'
  if(is.null(plotcontrol$main)) plotcontrol$main  <- 'Regression coefficients'
  if(is.null(plotcontrol$type)) plotcontrol$type  <- 'l'
  
  
  if(is.null(legendcontrol$x)) legendcontrol$x = 'topright'
  if(is.null(legendcontrol$bg)) legendcontrol$bg = 'white'

  if(all(indices=='AR')) indices <- matrix(1:nlatent,nrow=nlatent,ncol=2)
  
  if(all(indices=='CR')) indices <- cbind(
    rep(1:nlatent,nlatent)[-seq(1,nlatent^2,nlatent+1)],
    rep(1:nlatent,each=nlatent-1))
  
  if(indices[1]=='all') indices <- cbind(
    rep(1:nlatent,nlatent),
    rep(1:nlatent,each=nlatent))
  
  if(!'array' %in% class(indices) && !'matrix' %in% class(indices)){#interpret as individual columns
    indices <- cbind(
      rep(1:nlatent,length(unique(indices))),
      rep(unique(indices),each=length(unique(indices))))
  }
  
  if(ltyvec[1]=='auto') ltyvec=1:nrow(indices)
  if(lwdvec[1]=='auto') lwdvec= rep(3,nrow(indices))
  
  if(colvec[1]=='auto') colvec=grDevices::rainbow(nrow(indices),alpha=.8,v=.9)
  
  if(is.null(plotcontrol$ylim)) {
    plotcontrol$ylim=range(plyr::aaply(x,c(3,4),function(x) 
      x[indices]),na.rm=TRUE) #range of diagonals
    if(legend) plotcontrol$ylim[2] <- plotcontrol$ylim[2] + sd(plotcontrol$ylim)/3
  }
  
  if(gg){
    parnames<-paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
    y = x
    y = array(y,dim = c(dim(y)[1]^2,dim(y)[c(3,4)]))
    y <- y[matrix(1:dim(y)[1],sqrt(dim(y)[1]),sqrt(dim(y)[1]))[indices],,,drop=FALSE]
    dimn <- list(Index=parnames)
    dimn <- c(dimn,dimnames(x)[3:4])
    names(dimn) <- c('Index','Time interval','Effect')
    dimnames(y) <- dimn
    g <- ctPlotArrayGG(list(x=times,y=aperm(y,c(2,1,3))),...)
    # if(plot) print(g)
    # if(!plot) return(invisible(g))
    return(g)
  } else {
    
    
    
    #blank plot
    blankargs=plotcontrol
    blankargs$xlim=range(times)
    blankargs$y=NA
    blankargs$x=NA
    do.call(plot,blankargs)
    if(grid) {
      grid()
      par(new=TRUE)
      do.call(plot,blankargs)
      par(new=FALSE)
    }
    
    ####plotting confidence region
    if(polygon) {
      cc=0
      ccup=TRUE
      for(indexi in c(1:nrow(indices))){
        cc=ifelse(ccup,cc+1,cc-1)
        if(indexi==nrow(indices)) ccup=FALSE
        ri=indices[indexi,1]
        ci=indices[indexi,2]
        polygonargs<-polygoncontrol
        polygonargs$x=times
        polygonargs$y=x[ri,ci,,2]
        polygonargs$ylow=x[ri,ci,,1]
        polygonargs$yhigh=x[ri,ci,,length(quantiles)]
        polygonargs$col=grDevices::adjustcolor(colvec[cc],alpha.f=max(c(.004,polygonalpha/(2*sqrt(polygonargs$steps)))))
        do.call(ctPoly,polygonargs)
      }
    }
    
    ####plotting quantile lines
    for(qi in 1:3){
      cc=0
      for(indexi in 1:nrow(indices)){
        cc=cc+1
        ri=indices[indexi,1]
        ci=indices[indexi,2]
        
        plotargs<-plotcontrol
        plotargs$x=times
        plotargs$y=x[ri,ci,,qi]
        plotargs$lty=ltyvec[cc]
        plotargs$col=ifelse(qi!= 2,grDevices::adjustcolor(colvec[cc],alpha.f=.5) ,colvec[cc])
        plotargs$lwd=ifelse(qi!= 2,1, lwdvec[cc])
        do.call(points,plotargs)
      }}
    
    
    
    legendcontrol$legend=paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
    legendcontrol$text.col=colvec
    legendcontrol$col=colvec
    legendcontrol$lty = ltyvec
    legendcontrol$lwd=lwdvec
    
    if(legend) do.call(graphics::legend,legendcontrol)
  }#end if not gg
}
