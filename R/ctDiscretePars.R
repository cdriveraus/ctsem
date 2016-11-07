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
#'pars <- ctStanGetPars(ctstantestfit)
#'ctDiscretePars(pars,times=c(.5,1))
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
  
  return(out)
}


#'ctStanDiscretePars
#'
#'Generate discrete time parameters for a sequence of times based on a continuous time model fit
#'from ctStanFit, for specified subjects.
#'
#'@param ctstanfitobj Continuous time model fit from \link{\code{ctStanFit}}
#'@param subjects Either 'all', to take the average over all subjects, or a vector of integers denoting which
#'subjects.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param type String. 'all' returns all possible outputs in a list. 
#''discreteDRIFT' returns only discrete time auto and cross regression effects.
#''latentMeans' returns only the expected latent means, given initial (T0MEANS) level, 
#'latent intercept (CINT) and temporal effects (DRIFT).
#'@examples
#'ctStanDiscretePars(ctstantestfit,times=c(.5,1))
#'@export
ctStanDiscretePars<-function(ctstanfitobj, subjects='all', times=seq(from=0,to=10,by=.1), 
  quantiles = c(.05, .5, .95),type='all'){
  
  collapseSubjects=TRUE
  e<-extract(ctstanfitobj)
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
    
    out[[typei]] <- aaply(1:niter,1,function(iterx){
      if(collapseSubjects) subjectvec=1 else subjectvec=1:nsubjects #average over subjects before computing? much faster but answers dif question.
      aaply(subjectvec,1,function(subjecty){
        ctparsxy <- llply(ctpars, function(obji) {
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
  
  return(out)
}

#'ctStanPlot
#'
#'Plotting functions for continuous time models fit with ctStanFit
#'@param x fit object from \link{\code{ctStanFit}}
#'@param indices Either a string specifying type of plot to create, or an n by 2
#'matrix specifying which indices of the output matrix to plot.
#''AR' specifies all diagonals, for discrete time autoregression parameters.
#''CR' specifies all off-diagonals,for discrete time cross regression parameters.
#''all' plots all AR and CR effects at once.
#''latentMeans' for expected latent means, given initial (T0MEANS) level, 
#'latent intercept (CINT) and temporal effects (DRIFT).
#'@param add Logical. If FALSE, a new plot is generated, if TRUE, specified plot/s are
#'overlayed on existing plot.
#'@param legend Logical. If TRUE, generates a legend.
#'@param polygon Logical. If TRUE, fills a polygon between the first and last specified quantiles.
#'@param quantiles numeric vector of values between 0 and 1, specifying which quantiles to plot.
#'The default of c(.05,.5,.95) plots 95\% credible intervals and the posterior mean at 50\%. 
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each.
#'@param ylim Either 'auto' to determine automatically, or vector of length 2 specifying
#'upper and lower limits of y axis.
#'@param lwdvec Either 'auto', or a vector of positive integers denoting line widths for each quantile.
#''auto' specifies c(1,3,1) if there are 3 quantiles to be plotted (default), otherwise simply 3.
#'@param ltyvec Either 'auto', or a vector of line type integers (as for the lty parameter normally)
#' denoting line types for each quantile.
#' 'auto' specifies c(3, 1, 3) if there are 3 quantiles to be plotted (default), otherwise simply 1.
#'@param colorvec Either 'auto', or a vector of color values denoting colors for each index to be plotted.
#''auto' generates colors using the \link{\code{rainbow}} function.
#'@param plotcontrol list of arguments to pass to plot function. 
#'The following arguments are ignored: ylim,lwd,lty,col,x,y.
#'@param legendcontrol list of arguments to pass to legend function. 'legend=' and 'text.col=' arguments
#'will be ignored.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha (transparency) of colorvec by for 
#'the fill polygon.
#'@param polygoncontrol list of arguments to pass to polgyon function (if polygon=TRUE).
#'x,y, and col arguments will be ignored.
#'@param ... additional arguments to pass to \link{\code{ctStanDiscretePars}},
#'such as subjects, or times.
#'@examples
#'ctStanPlot(ctstantestfit, 'AR', 2)
#'@export

ctStanPlot<- function(x,indices='AR',add=FALSE,legend=TRUE, polygon=TRUE, 
  quantiles=c(.05,.5,.95), times=seq(0,10,.1),
  ylim='auto',lwdvec='auto',colorvec='auto',ltyvec='auto',
  plotcontrol=list(ylab='Value',xlab='Time',type='l'),
  legendcontrol=list(x='topright',bg='white'),
  polygonalpha=.05,
  polygoncontrol=list(border=NA),
  ...){
  
  input <- ctStanDiscretePars(x,type='discreteDRIFT',times=times,quantiles=quantiles,...)[[1]]
  
  nlatent=dim(input)[1]
  
  latentNames=paste0('eta',1:nlatent)
  
  if(ltyvec=='auto' & length(quantiles)==3) ltyvec=c(3,1,3) else ltyvec=1
  if(lwdvec=='auto' & length(quantiles)==3) lwdvec=c(1,3,1) else lwdvec=3
  
  if(indices[1]=='AR') indices <- matrix(1:nlatent,nrow=nlatent,ncol=2)
  
  if(indices[1]=='CR') indices <- cbind(
    rep(1:nlatent,nlatent)[-seq(1,nlatent^2,nlatent+1)],
    rep(1:nlatent,each=nlatent-1))
  
  if(indices[1]=='all') indices <- cbind(
    rep(1:nlatent,nlatent),
    rep(1:nlatent,each=nlatent))
  
  if(colorvec[1]=='auto') colorvec=rainbow(nrow(indices),alpha=.8)
  
  if(ylim=='auto') ylim=range(plyr::aaply(input,c(3,4),function(x)
    x[indices]),na.rm=TRUE) #range of diagonals
  
  for(qi in 1:length(quantiles)){
    cc=0
    for(indexi in 1:nrow(indices)){
      cc=cc+1
      ri=indices[indexi,1]
      ci=indices[indexi,2]
      
      plotargs<-plotcontrol
      plotargs$x=times
      plotargs$y=input[ri,ci,,qi]
      plotargs$lty=ltyvec[qi]
      plotargs$col=colorvec[cc]
      plotargs$lwd=lwdvec[qi]
      plotargs$lty=ltyvec[qi]
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
      polygonargs$col=grDevices::adjustcolor(colorvec[cc],alpha.f=polygonalpha)
      do.call(graphics::polygon,polygonargs)
    }
  }
  
  legendcontrol$legend=paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
  legendcontrol$text.col=colorvec
  
  if(legend) do.call(graphics::legend,legendcontrol)
  
}
