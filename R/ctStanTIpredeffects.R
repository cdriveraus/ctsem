#' Get time independent predictor effect estimates
#' 
#' Computes and plots combined effects and quantiles for effects of time independent predictors
#' on subject level parameters of a ctStanFit object.
#'
#' @param fit fit object from \code{\link{ctStanFit}}
#' @param returndifference logical. If FALSE, absolute parameter values are returned. 
#' If TRUE, only the effect of the covariate (i.e. without the average value of the parameter)
#' are returned. The former can be easier to interpret, but the latter are more likely to fit multiple plots together. 
#' Not used if \code{parmatrices=TRUE}.
#' @param probs numeric vector of quantile probabilities from 0 to 1. Specify 3
#' values if plotting, the 2nd will be drawn as a line with uncertainty polygon
#' based on 1st and 3rd.
#' @param includeMeanUncertainty if TRUE, output includes sampling variation in the mean parameters. If FALSE,
#' mean parameters are fixed at their median, only uncertainty in time independent predictor effects is included.  
#' @param whichTIpreds integer vector specifying which of the tipreds in the fit object you want to
#' use to calculate effects. Unless quadratic / higher order versions of predictors have been 
#' included, selecting more than one probably doesn't make sense. If for instance a squared
#' predictor has been included, then you can specify both the linear and squared version. 
#' The x axis of the plot (if generated) will be based off the first indexed predictor. To 
#' check what predictors are in the model, run \code{fit$ctstanmodel$TIpredNames}.
#' @param parmatrices Logical. If TRUE (default), the \code{\link{ctStanParMatrices}} function
#' is used to return an expanded range of possible matrices of interest.
#' @param whichpars if parmatrices==TRUE, character vector specifying which matrices, and potentially which 
#' indices of the matrices, to plot. c('dtDRIFT[2,1]', 'DRIFT') would output for row 2 and column 1 of 
#' the discrete time drift matrix, as well as all indices of the continuous time drift matrix. 
#' If parmatrices==FALSE, integer vector specifying which of the individually varying subject
#' level parameters to compute effects on. In either case 'all' uses all available.
#' The integer corresponding to specific parameters can be found as follows, replacing \code{fit} as appropriate:
#' \code{fit$ctstanmodel$pars[fit$ctstanmodel$pars$indvarying,'param']}.
#' @param niterations Positive integer specifying the maximum number of saved iterations to use. 
#' Character string 'all' can also be used.
#' @param nsubjects Positive integer specifying the number of subjects to compute values for.
#' Character string 'all' can also be used. Time taken is a function of nsubjects*niterations.
#' @param plot Logical. If TRUE, nothing is returned but instead \code{\link{ctPlotArray}}
#' is used to plot the output instead.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter matrices,
#' if \code{parmatrices=TRUE}.
#' @param ... arguments to pass to \code{\link{ctPlotArray}} for plotting.
#' @return Either a three dimensional array of predictor effects, or nothing with a plot
#' generated.
#' @export
#'
#' @examples
#' #samples reduced here for speed
#' ctStanTIpredeffects(ctstantestfit,plot=TRUE,whichpars='dtDRIFT',nsamples=10)
ctStanTIpredeffects<-function(fit,returndifference=FALSE, probs=c(.025,.5,.975),
  includeMeanUncertainty=FALSE,
  whichTIpreds=1,parmatrices=TRUE, whichpars='all', nsamples=50, timeinterval=1,
  nsubjects=50,
  plot=FALSE,...){

  #drop fixed and duplicated params
  spec_nofixed <- fit$ctstanmodel$pars[is.na(fit$ctstanmodel$pars$value),,drop=FALSE]
  spec_nofixed_noduplicates <- spec_nofixed[!duplicated(spec_nofixed$param),]
  
  #get indvarying rawpopmeans
  e<-extract(fit$stanfit)
  rawpopmeans <- e$rawpopmeans
  
  niter<-dim(e$rawpopmeans)[1]
  if(nsamples > niter) nsamples <- niter
  rawpopmeans <- rawpopmeans[sample(x = 1:niter, nsamples,replace = FALSE),]
  
  if(!includeMeanUncertainty) rawpopmeans <- matrix(apply(rawpopmeans,2,median),byrow=TRUE,nrow=nrow(rawpopmeans),ncol=ncol(rawpopmeans))
  rawpopmeansindvarying <- rawpopmeans[,spec_nofixed_noduplicates$indvarying,drop=FALSE] 
  
  tipreds<-fit$data$tipreds[,whichTIpreds,drop=FALSE]
  if(nsubjects=='all') nsubjects = nrow(tipreds)
  if(nsubjects > nrow(tipreds)) nsubjects <- nrow(tipreds)
  tipreds <- tipreds[sample(x = 1:nsubjects, nsubjects,replace = FALSE),]

  tieffect<-e$tipredeffect[,,whichTIpreds,drop=FALSE]
  
  
  if(!parmatrices){ #update ctspec to only indvarying and those in whichpars
    spec_nofixed_noduplicates_indvarying <- spec_nofixed_noduplicates[spec_nofixed_noduplicates$indvarying,]
   if(all(whichpars=='all')) whichpars=1:sum(spec_nofixed_noduplicates_indvarying$indvarying)
    spec_nofixed_noduplicates_indvarying <- spec_nofixed_noduplicates_indvarying[whichpars,,drop=FALSE]
    rawpopmeansindvarying <- rawpopmeansindvarying[,whichpars,drop=FALSE]  #then just the ones in whichpars
   tieffect<-tieffect[,whichpars,,drop=FALSE] #updating...
    npars<-length(whichpars)
  }
  
  tiorder<-order(tipreds[,1])
  tipreds<-tipreds[tiorder,,drop=FALSE] #order tipreds according to first one
  
  message('Calculating time independent predictor effects...')

  raweffect <- aaply(1:nrow(rawpopmeans),1,function(iterx) { #for every iter
    aaply(tipreds,1,function(tix){ #and every distinct tipred vector
      rawpopmeansindvarying[iterx,] + matrix(tieffect[iterx,,],nrow=dim(tieffect)[2]) %*% tix
    },.drop=FALSE)
  })
  
  
  if(!parmatrices) {
    effect<-aaply(1:npars, 1,function(pari){ #for each param
      param=raweffect[,,pari]
      out=eval(parse(text=spec_nofixed_noduplicates_indvarying$transform[pari]))
      return(out)
    })
    if(returndifference){ #if only returning differences from zero
      noeffect<-aaply(1:npars, 1,function(pari){ #for each param
        param <- rawpopmeans[,pari]
        out=eval(parse(text=spec_nofixed_noduplicates_indvarying$transform[pari]))
        return(out)
      })
      effect<-effect-array(noeffect,dim=dim(effect))
    }
  }
  

  if(parmatrices)  {
    rawpopmeans <- rawpopmeans[rep(1:nrow(rawpopmeans),each=nsubjects),] #match rows of rawpopmeans and raweffect
    raweffect <- matrix(raweffect,ncol=dim(raweffect)[3])
    
    parmatlists<-lapply(1:nrow(rawpopmeans), function(x) { #for each param vector
      parvec = rawpopmeans[x,]
      parvec[spec_nofixed_noduplicates$indvarying] <- raweffect[x,]
      out = ctStanParMatrices(fit,parvec,timeinterval=timeinterval)
      return(out)
    })
    
    
    # parmatlists <- apply(e$rawpopmeans,1,ctStanParMatrices,model=object,timeinterval=timeinterval)
    parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
    parmats <- matrix(0,nrow=0,ncol=2)
    counter=0
    for(mati in 1:length(parmatlists[[1]])){
        for(coli in 1:ncol(parmatlists[[1]][[mati]])){
          for(rowi in 1:nrow(parmatlists[[1]][[mati]])){
          counter=counter+1
          new <- matrix(c(
            rowi,
            coli),
            nrow=1)
          rownames(new) = paste0(names(parmatlists[[1]])[[mati]])
          parmats<-rbind(parmats, new)
        }}}
    colnames(parmats) <- c('Row','Col') 

    rownames(parmatarray) <- paste0(rownames(parmats),'[',parmats[,'Row'],',',parmats[,'Col'],']')
    
    #remove certain parmatrices lines
    removeindices <- which(rownames(parmats) == 'MANIFESTVAR' & parmats[,'Row'] != parmats[,'Col'])
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('MANIFESTVAR','T0VAR','DIFFUSION','dtDIFFUSION','asymDIFFUSION',
      'T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') &  parmats[,'Row'] < parmats[,'Col'])))
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') & 
        parmats[,'Row'] == parmats[,'Col'])))
    
    parmatarray <- parmatarray[-removeindices,]
    
    effect <- array(parmatarray,dim=c(nrow(parmatarray),nsamples,nsubjects))
    rownames(effect) <- rownames(parmatarray)
    if(any(whichpars !='all')) {
      selection <- unlist(lapply(whichpars,function(x) grep(paste0('^\\Q',x,'\\E'),dimnames(effect)[[1]])))
      effect <- effect[selection,,,drop=FALSE]
    }
    
  }    

  
  out<-aaply(probs,1,function(x) ctCollapse(effect,2,quantile,probs=x,na.rm=TRUE),.drop=FALSE)
  
  if(!parmatrices) dimnames(out)=list(Quantile=paste0('Quantile',probs),
    popmean=spec_nofixed_noduplicates_indvarying$param,
    subject=tiorder #subjects reordered because tipreds were at top
  )
  if(parmatrices) dimnames(out)=list(Quantile=paste0('Quantile',probs),
    param=rownames(effect),
    subject=tiorder #subjects reordered because tipreds were at top
  )
  
  out <- list(y=aperm(out, c(3,2,1)), x=tipreds[,1])

  if(!plot) return(out) else {
    dots <- list(...)
    dots$yarray=out
    dots$x=tipreds[,1]
    if(is.null(dots$plotcontrol)) dots$plotcontrol=list(
      ylab=ifelse(!returndifference,'Par. Value','Effect'),
      xlab=colnames(tipreds)[1],
      xaxs='i')
    
    do.call(ctPlotArray,dots)
  }
}

