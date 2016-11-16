#' ctStanParnames
#' 
#' Gets internal stan parameter names of a ctStanFit object based on specified substrings.
#'
#' @param x ctStanFit object
#' @param substrings vector of character strings, parameter names of the stan model
#' containing any of these strings will be returned. Useful strings may be 'hmean_' for 
#' hyper (population) means, 'hsd_' for hyper standard deviations, 'tipred_' for time
#' independent predictors, or specific combinations such as 'hmean_drift' for the population
#' means of temporal dynamics parameters
#' @return vector of character strings.
#' @examples
#' ctStanParnames(ctstantestfit,substrings=c('hmean_','hsd_'))
#' @export
ctStanParnames <- function(x,substrings=c('hmean_','hsd_')){
  out<-c()
  for(subsi in substrings){
    out<- c(out, x$stanfit@model_pars[grep(subsi,x$stanfit@model_pars)])
  }
  return(out)
}


#'ctStanContinuousPars
#'
#'Returns the continuous time parameter matrices for specified subjects of a ctStanFit fit object
#'
#'@param ctstanfitobj fit object from \code{\link{ctStanFit}}
#'@param subjects Either 'all', or integers denoting which subjects to perform the calculation over. 
#'When multiple subjects are specified, the returned matrices will be a mean over subjects.
#'@param iter Either character string 'all' which will then use all post-warmup iterations, 
#'or an integer specifying which iteration/s to use.
#'@param calcfunc Function to apply over samples, must return a single value. 
#'By default the mean over all samples is returned, but one might also be interested in
#'the \code{\link[stats]{sd}} or \code{\link[stats]{quantile}} functions.
#'@param ... additional parameters to pass to calcfunc. For instance, with calcfunc = quantile, 
#'the probs argument is needed to ensure only a single value is returned.
#'@examples
#'#posterior mean over all subjects
#'ctStanContinuousPars(ctstantestfit)
#'
#'#posterior 95% quantiles for subject 2
#'ctStanContinuousPars(ctstantestfit, subjects=2, calcfunc=quantile, probs = .95)
#'@export
ctStanContinuousPars <- function(ctstanfitobj,subjects='all',iter='all',
  calcfunc=mean,...){
  
  if(subjects[1] != 'all' && !is.integer(as.integer(subjects))) stop('
    subjects argument must be either "all" or an integer denoting specific subjects')
  
  if(class(ctstanfitobj)!='ctStanFit') stop('Not an object of class ctStanFit')
  
  e<-rstan::extract(ctstanfitobj$stanfit) #first dim of subobjects is iter, 2nd subjects
  niter=dim(e$DRIFT)[1]
  
  if(iter!='all') {
    if(!any(iter %in% 1:niter)) stop('Invalid iteration specified!')
    e=lapply(e, function(x) {
      xdims=dim(x)
      out=array(eval(parse(text=
          paste0('x[iter',if(length(xdims)>1) paste0(rep(',',length(xdims)-1),collapse=''),']')
      )),dim=c(length(iter),xdims[-1]))
      return(out)
    }
    )
  }
  
  nsubjects <- dim(e$indparams)[2]
  if(is.null(nsubjects)) nsubjects=1
  
  if(subjects[1]=='all') subjects=1:nsubjects
  
  collapsemargin<-c(1,2)
  # if(collapseIterations) collapsemargin=1
  # if(collapseSubjects) collapsemargin=c(collapsemargin,2)
  
  for(matname in c('DRIFT','DIFFUSION','CINT','T0MEANS', 
    'T0VAR','MANIFESTMEANS','MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    
    if(dim(e[[matname]])[2] > 1) subselection <- subjects else subselection <- 1
    
    vector <- FALSE
    
    if(matname %in% c('T0MEANS','CINT', 'MANIFESTMEANS')) vector <- TRUE
    
    if(!vector) assign(matname, 
      array(ctCollapse(inarray = e[[matname]][,subselection,,,drop=FALSE],
        collapsemargin = collapsemargin, 
        collapsefunc = calcfunc,...),dim=dim(e[[matname]])[-c(1,2)])
    )
    
    if(vector) assign(matname,
      array(ctCollapse(inarray = e[[matname]][,subselection,,drop=FALSE],
        collapsemargin = collapsemargin, 
        collapsefunc = calcfunc, 
        ...),dim=c(dim(e[[matname]])[-c(1,2)],1))
    )
    
  }
  
  model<-list(DRIFT=DRIFT,T0VAR=T0VAR,DIFFUSION=DIFFUSION,CINT=CINT,T0MEANS=T0MEANS,
    MANIFESTMEANS=MANIFESTMEANS,MANIFESTVAR=MANIFESTVAR, LAMBDA=LAMBDA)
  
  if(!is.null(e$TDPREDEFFECT)) model$TDPREDEFFECT<-TDPREDEFFECT
  
  return(model)
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
#'pars <- ctStanContinuousPars(ctstantestfit)
#'ctDiscretePars(pars,times=c(.5,1))
#'
#'@export
ctDiscretePars<-function(ctpars,times=seq(0,10,.1),type='all'){
  
  if(type=='all') type=c('discreteDRIFT','latentMeans') #this needs to match with ctStanDiscretePars
  nlatent=nrow(ctpars$DRIFT)
  latentNames=paste0('eta',1:nlatent)
  
  out<-list()
  
  discreteDRIFT = array(unlist(lapply(times, function(x) 
    OpenMx::expm(ctpars$DRIFT*x))),
    dim=c(nlatent,nlatent,length(times)),
    dimnames=list(latentNames,latentNames,paste0('t',times)))
  
  if('discreteDRIFT' %in% type) out$discreteDRIFT = discreteDRIFT
  
  if('latentMeans' %in% type) out$latentMeans = array(unlist(lapply(1:length(times), function(x)
    discreteDRIFT[,,x] %*% ctpars$T0MEANS)),
    dim=c(nlatent,length(times)),
    dimnames=list(latentNames,paste0('t',times)))
  
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
#'@param quantiles Which quantiles to return
#'@param plot Logical. If TRUE, plots output using \code{\link{ctStanDiscreteParsPlot}}
#'instead of returning output. 
#'@param ... additional plotting arguments to control \code{\link{ctStanDiscreteParsPlot}}
#'@examples
#'ctStanDiscretePars(ctstantestfit,times=seq(.5,4,.1), 
#'plot=TRUE,indices='all')
#'@export
ctStanDiscretePars<-function(ctstanfitobj, subjects='all', times=seq(from=0,to=10,by=.1), 
  quantiles = c(.05, .5, .95),plot=FALSE,...){
  
  type='discreteDRIFT'
  collapseSubjects=TRUE #consider this for a switch
  
  e<-rstan::extract(ctstanfitobj$stanfit)
  if(type=='all') type=c('discreteDRIFT','latentMeans') #must match with ctDiscretePars
  
  if(subjects[1] != 'all' && !is.integer(as.integer(subjects))) stop('
  subjects argument must be either "all" or an integer denoting specific subjects')
  
  nsubjects <- dim(e$indparams)[2]
  if(is.null(nsubjects)) nsubjects=1
  if(subjects[1]=='all') subjects=1:nsubjects
  
  outdims=dim(e$DRIFT)
  niter=outdims[1]
  nlatent=outdims[3]
  latentNames=paste0('eta',1:nlatent)
  
  out<-list()
  
  #get all ctparameter matrices at once and remove unneeded subjects
  ctpars <- list()
  for(matname in c('DRIFT','DIFFUSION','CINT','T0MEANS', 
    'T0VAR','MANIFESTMEANS','MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    
    vector <- FALSE
    if(matname %in% c('T0MEANS','CINT', 'MANIFESTMEANS')) vector <- TRUE
    
    xdims=dim(e[[matname]])
    dimout=xdims
    dimout[2]=min(length(subjects),xdims[2])
    
    ctpars[[matname]] <- array(eval(parse(text=
        paste0('e[[matname]][, ',
          ifelse(xdims[2] > 1, 'subjects',1),if(length(xdims)>1) paste0(rep(', ',length(xdims)-2),collapse=''),']')
    )),dim=dimout)
    
    ctpars[[matname]] <-aperm(ctpars[[matname]],c(3,if(length(dim(ctpars[[matname]]))>3) 4, 1, 2))
  }
  
  nsubjects <- length(subjects)
  
  
  for(typei in 1:length(type)){ #for all types of discrete parameter requested, compute pars
    message('Getting ',typei,' / ', length(type), ', ',type[typei])
    matrixtype=ifelse(type[typei] %in% c('discreteDRIFT'),TRUE, FALSE) #include any matrix type outputs
    
    out[[typei]] <- plyr::aaply(1:niter,1,function(iterx){
      if(collapseSubjects) subjectvec=1 else subjectvec=1:nsubjects #average over subjects before computing? much faster but answers dif question.
      plyr::aaply(subjectvec,1,function(subjecty){
        ctparsxy <- plyr::llply(ctpars, function(obji) {
          ismatrix = length(dim(obji)) > 3 #check if obji is a matrix or vector
          out=eval(parse(text=paste0('obji[,',
            if(ismatrix) ',',
            iterx,',',
            ifelse(dim(obji)[ifelse(ismatrix,4,3)] > 1, ifelse(collapseSubjects, '1:nsubjects', 'subjecty'),1),
            ']')))
          
          if(collapseSubjects & dim(obji)[ifelse(ismatrix,4,3)] > 1) out <- ctCollapse(out,
            collapsemargin = ifelse(ismatrix,3,2),
            collapsefunc = mean)
          
          return(out)
        })
        
        
        
        discreteparsxy <- ctDiscretePars(ctparsxy,
          times=times,
          type=type[typei])[[1]]
        
        return(discreteparsxy)
      },.drop=FALSE)
    },.drop=FALSE,.progress='text')
    
    out[[typei]] = plyr::aaply(quantiles, 1, function(quantx) 
      ctCollapse(inarray=out[[typei]],collapsemargin=c(1,if(collapseSubjects) 2),collapsefunc=quantile,probs=quantx)
    )
    
    
    dimlist<- list(quantiles=paste0('quantile_',quantiles),
      row=dimnames(out[[typei]])[[2]],
      col=dimnames(out[[typei]])[[3]],
      times=paste0('t',times)
    )
    if(!matrixtype) dimlist[[3]]=NULL
    
    dimnames(out[[typei]])=dimlist
    
    out[[typei]]=aperm(out[[typei]],c(2,3,if(matrixtype) 4,1))
  }
  
  names(out) <- type
  
  if(plot) {
    ctStanDiscreteParsPlot(out,times=times,latentNames=ctstanfitobj$ctstanmodel$latentNames,...)
  } else return(out)
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
#'@param legend Logical. If TRUE, generates a legend.
#'@param polygon Logical. If TRUE, fills a polygon between the first and last specified quantiles.
#'@param quantiles numeric vector of length 3, with values between 0 and 1, specifying which quantiles to plot.
#'The default of c(.05,.5,.95) plots 95\% credible intervals and the posterior median at 50\%. 
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param latentNames Vector of character strings denoting names for the latent variables. 
#''auto' just uses eta1 eta2 etc.
#'@param ylim Either 'auto' to determine automatically, or vector of length 2 specifying
#'upper and lower limits of y axis.
#'@param lwdvec Either 'auto', or a vector of positive integers denoting line widths for each quantile.
#''auto' specifies c(1,3,1) if there are 3 quantiles to be plotted (default), otherwise simply 3.
#'@param ltyvec Either 'auto', or a vector of line type integers (as for the lty parameter normally)
#' denoting line types for each quantile.
#' 'auto' specifies c(3, 1, 3) if there are 3 quantiles to be plotted (default), otherwise simply 1.
#'@param colvec Either 'auto', or a vector of color values denoting colors for each index to be plotted.
#''auto' generates colors using the \code{\link[grDevices]{rainbow}} function.
#'@param plotcontrol list of arguments to pass to plot function. 
#'The following arguments are ignored: ylim,lwd,lty,col,x,y.
#'@param legendcontrol list of arguments to pass to legend function. 'legend=' and 'text.col=' arguments
#'will be ignored.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha (transparency) of colvec by for 
#'the fill polygon.
#'@param polygoncontrol list of arguments to pass to polgyon function (if polygon=TRUE).
#'x,y, and col arguments will be ignored.
#'@examples
#'x <- ctStanDiscretePars(ctstantestfit)
#'
#'ctStanDiscreteParsPlot(x, 'CR')
#'@export

ctStanDiscreteParsPlot<- function(x,indices='all',add=FALSE,legend=TRUE, polygon=TRUE, 
  quantiles=c(.05,.5,.95), times=seq(0,10,.1),latentNames='auto',
  ylim='auto',lwdvec='auto',colvec='auto',ltyvec='auto',
  plotcontrol=list(ylab='Value',xlab='Time interval',
    main='Regression coefficients',type='l'),
  legendcontrol=list(x='topright',bg='white'),
  polygonalpha=.1,
  polygoncontrol=list(border=NA)){
  
  input <- x[[1]] #ctStanDiscretePars(x,type='discreteDRIFT',times=times,quantiles=quantiles,...)[[1]]
  
  nlatent=dim(input)[1]
  
  if(latentNames[1]=='auto') latentNames=paste0('eta',1:nlatent)
  
  if(is.null(plotcontrol$ylab)) plotcontrol$ylab  <- 'Value'
  if(is.null(plotcontrol$xlab)) plotcontrol$xlab  <- 'Time interval'
  if(is.null(plotcontrol$main)) plotcontrol$main  <- 'Regression coefficients'
  if(is.null(plotcontrol$type)) plotcontrol$type  <- 'l'
  
  if(is.null(legendcontrol$x)) legendcontrol$x = 'topright'
  if(is.null(legendcontrol$bg)) legendcontrol$bg = 'white'
  
  
  if(indices[1]=='AR') indices <- matrix(1:nlatent,nrow=nlatent,ncol=2)
  
  if(indices[1]=='CR') indices <- cbind(
    rep(1:nlatent,nlatent)[-seq(1,nlatent^2,nlatent+1)],
    rep(1:nlatent,each=nlatent-1))
  
  if(indices[1]=='all') indices <- cbind(
    rep(1:nlatent,nlatent),
    rep(1:nlatent,each=nlatent))
  
  if(ltyvec[1]=='auto') ltyvec=1:nrow(indices)
  if(lwdvec[1]=='auto') lwdvec= rep(3,nrow(indices))
  
  if(colvec[1]=='auto') colvec=grDevices::rainbow(nrow(indices),alpha=.8)
  
  if(ylim=='auto') ylim=range(plyr::aaply(input,c(3,4),function(x)
    x[indices]),na.rm=TRUE) #range of diagonals
  
  ####plotting quantiles - now just using polygon so only plotting middle quantile
  for(qi in 1:3){
    cc=0
    for(indexi in 1:nrow(indices)){
      cc=cc+1
      ri=indices[indexi,1]
      ci=indices[indexi,2]
      
      plotargs<-plotcontrol
      plotargs$x=times
      plotargs$y=input[ri,ci,,qi]
      plotargs$lty=ltyvec[cc]
      plotargs$col=ifelse(qi!= 2,grDevices::adjustcolor(colvec[cc],alpha.f=sqrt(polygonalpha)) ,colvec[cc])
      plotargs$lwd=ifelse(qi!= 2,1, lwdvec[cc])
      plotargs$ylim=ylim
      
      if(qi==1 && indexi==1 & !add) do.call(plot,plotargs) else do.call(points,plotargs)
    }}
  
  if(polygon) {
    cc=0
    backwardstimesindex=order(times,decreasing=TRUE)
    for(indexi in 1:nrow(indices)){
      cc=cc+1
      ri=indices[indexi,1]
      ci=indices[indexi,2]
      polygonargs<-polygoncontrol
      polygonargs$x=c(times,times[backwardstimesindex])
      polygonargs$y=c(input[ri,ci,,1], input[ri,ci,,length(quantiles)][backwardstimesindex])
      polygonargs$col=grDevices::adjustcolor(colvec[cc],alpha.f=polygonalpha)
      do.call(graphics::polygon,polygonargs)
    }
  }
  
  legendcontrol$legend=paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
  legendcontrol$text.col=colvec
  legendcontrol$col=colvec
  legendcontrol$lty = ltyvec
  legendcontrol$lwd=lwdvec
  
  if(legend) do.call(graphics::legend,legendcontrol)
  
}
