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
  if(length(x$stanfit$stanfit@sim)==0) stop('Doesnt contain sampled stanfit object')
  out<-c()
  for(subsi in substrings){
    out<- c(out, x$stanfit$stanfit@model_pars[grep(paste0('^',subsi),x$stanfit$stanfit@model_pars)])
  }
  return(out)
}



#'ctStanDiscretePars
#'
#'Calculate model implied regressions for a sequence of time intervals (if ct) or steps (if dt) based on
#'a ctStanFit object, for specified subjects.
#'
#'@param ctstanfitobj model fit from \code{\link{ctStanFit}}
#'@param subjects Either 'popmean', to use the population mean parameter, or a vector of integers denoting which
#'subjects.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each. If the fit 
#'object is a discrete time model, these should be positive integers.
#'@param nsamples Number of samples from the stanfit to use for plotting. Higher values will
#'increase smoothness / accuracy, at cost of plotting speed. Values greater than the total
#'number of samples will be set to total samples.
#'@param observational Logical. If TRUE, outputs expected change in processes *conditional on observing* a 1 unit change in each -- 
#'this change is correlated according to the DIFFUSION matrix. If FALSE, outputs expected regression values -- also interpretable as
#'an independent 1 unit change on each process, giving the expected response under a 1 unit experimental impulse.
#'@param standardise Logical. If TRUE, output is standardised according to expected total within subject variance, given by the 
#'asymDIFFUSIONcov matrix.
#'@param cov Logical. If TRUE, covariances are returned instead of regression coefficients.
#'@param plot Logical. If TRUE, plots output using \code{\link{ctStanDiscreteParsPlot}}
#'instead of returning output. 
#'@param cores Number of cpu cores to use for computing subject matrices. 
#'If subject matrices were saved during fiting, not used. 
#'@param ... additional plotting arguments to control \code{\link{ctStanDiscreteParsPlot}}
#'@examples
#'if(w32chk()){
#'
#' ctStanDiscretePars(ctstantestfit,times=seq(.5,4,.1), 
#'  plot=TRUE,indices='CR')
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
ctStanDiscretePars<-function(ctstanfitobj, subjects='popmean', 
  times=seq(from=0,to=10,by=.1), 
  nsamples=100,observational=FALSE,standardise=FALSE, 
  cov=FALSE, plot=FALSE,cores=2,...){
  
  if(!ctstanfitobj$ctstanmodel$continuoustime) times <- unique(round(times))
  type='discreteDRIFT'
  collapseSubjects=TRUE #consider this for a switch
  e<-ctExtract(ctstanfitobj,subjectMatrices = subjects[1]!='popmean',cores=cores)
  
  if(subjects[1] != 'popmean' && any(!is.integer(as.integer(subjects)))) stop('
  subjects argument must be either "all" or an integer denoting specific subjects')
  
  nsubjects <- dim(e$subj_DRIFT)[2]
  if(is.null(nsubjects)) nsubjects=1
  if('popmean' %in% subjects) subjects='popmean' 
  
  if(is.null(e$subj_DRIFT) && any(!subjects %in% 'popmean')) stop('No individual variation in DRIFT matrix found?? Try subjects="all"')
  
  niter=dim(e$pop_DRIFT)[1]
  nlatent=ctstanfitobj$standata$nlatent#outdims[3]
  latentNames=ctstanfitobj$ctstanmodel$latentNames
  
  if(nsamples > niter) nsamples <- niter
  
  out<-list()
  
  #get all ctparameter matrices at once and remove unneeded subjects
  ctpars <- list()
  
  samples <- sample(1:dim(e[['pop_DRIFT']])[1],nsamples)
  
  # browser()
  for(matname in c('DRIFT','DIFFUSIONcov','asymDIFFUSIONcov')){ #,'CINT','T0MEANS', 'T0VAR','MANIFESTMEANS',if(!is.null(e$MANIFESTVAR)) 'MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    if('popmean' %in% subjects || is.null(e[[paste0('subj_',matname)]])){
      ctpars[[matname]] <- e[[paste0('pop_',matname)]][samples,,,drop=FALSE]
    } else {
      # browser()
      ctpars[[matname]] <- e[[paste0('subj_',matname)]][samples,subjects,,,drop=FALSE]
      # ctpars[[matname]]<-  array(ctpars[[matname]],dim=c(prod(dim(ctpars[[matname]])[1:2]),dim(ctpars[[matname]])[-1:-2]))
    }
    # ctpars[[matname]] <- ctpars[[matname]][sample(1:dim(ctpars[[matname]])[1],nsamples),,,drop=FALSE]
  }
  
  
  out <- ctStanDiscreteParsDrift(ctpars,times, observational, standardise, cov=cov,discreteInput = ctstanfitobj$ctstanmodel$continuoustime==FALSE)
  
  dimnames(out)<- list(Sample=samples, Subject=subjects,
    `Time interval`=times, row=latentNames, col=latentNames)
  
  attributes(out)$observational <- observational
  attributes(out)$cov <- cov
  
  if(plot) {
    
    out <- ctStanDiscreteParsPlot(out,
      # latentNames=ctstanfitobj$ctstanmodel$latentNames,
      ...)
  } 
  return(out)
}



ctStanDiscreteParsDrift<-function(ctpars,times, observational,  standardise,cov=FALSE,
  types='dtDRIFT',discreteInput=FALSE, quiet=FALSE){
  
  nl=dim(ctpars$DRIFT)[3]
  
  if(!quiet) message('Computing temporal regression coefficients for ', dim(ctpars$DRIFT)[1],' samples, may take a moment...')
  
  lapply(names(ctpars),function(x){ #add in extra dim if only 3 dims (e.g. when not individually varying)
    dm=dim(ctpars[[x]])
    if(length(dm)==3){
      ctpars[[x]] <<- array(ctpars[[x]],dim=c(dm[1],1,dm[2:3]))
    }
  })
  
  nsubs <- lapply(ctpars,function(x) dim(x)[2])
  
  
  if('dtDRIFT' %in% types){ 
    ctpars$dtDRIFT <- array(NA, dim=c(dim(ctpars$DRIFT)[1],max(unlist(nsubs)),length(times),dim(ctpars$DRIFT)[3:4]))
    
    mpow <- function(m,n){
      if(n==0) return(diag(1,nrow(m))) else{
        if(n>1){
          mo <-m
        for(i in 2:n){
          m <- m %*% mo
        }
        }
        return(m)
      }}
    
    for(i in 1:dim(ctpars$DRIFT)[1]){
      for(j in 1:dim(ctpars$DRIFT)[2]){
        for(ti in 1:length(times)){
          if(!discreteInput) ctpars$dtDRIFT[i,j,ti,,] <- expm::expm(as.matrix(ctpars$DRIFT[i,min(j,nsubs$DRIFT),,] * times[ti]))
          if(discreteInput) ctpars$dtDRIFT[i,j,ti,,] <- mpow(as.matrix(ctpars$DRIFT[i,min(j,nsubs$DRIFT),,]),times[ti])
          if(standardise) {
            if(any(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,]) < 0)) stop(
              "Asymptotic diffusion matrix has negative diagonals -- I don't know what non stationary standardization looks like")
            ctpars$dtDRIFT[i,j,ti,,] <- ctpars$dtDRIFT[i,j,ti,,] * 
              matrix(rep(sqrt(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,])+1e-10),each=nl) / 
                  rep((sqrt(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,]))),times=nl),nl)
          }
          if(observational){
            Qcor<-cov2cor(matrix(ctpars$DIFFUSIONcov[i,min(j,nsubs$DIFFUSIONcov),,],nl,nl)+diag(1e-8,nl)) 
            Qcor <- Qcor^2 * sign(Qcor) #why is this squared?
            ctpars$dtDRIFT[i,j,ti,,]  <- ctpars$dtDRIFT[i,j,ti,,]  %*% Qcor
          }
          if(cov) ctpars$dtDRIFT[i,j,ti,,]  <- tcrossprod(ctpars$dtDRIFT[i,j,ti,,] )
        }
      }
    }
  } #end dtdrift
  
  
  return(ctpars$dtDRIFT)
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
#'@param quantiles numeric vector of length 3, with values between 0 and 1, specifying which quantiles to plot.
#'The default of c(.05,.5,.95) plots 95\% credible intervals and the posterior median at 50\%. 
#'@param latentNames Vector of character strings denoting names for the latent variables. 
#''auto' just uses eta1 eta2 etc.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha of 
#'the fill.
#'@param ylab y label.
#'@param xlab x label.
#'@param ylim Custom ylim.
#'@param facets May be 'Subject' or 'Effect'.
#'@param colour Character string denoting how colour varies. 'Effect' or 'Subject'.
#'@param title Character string.
#'@param splitSubjects if TRUE, subjects are plotted separately, if FALSE they are combined.
#'@param ggcode if TRUE, returns a list containing the data.table to plot, and a character string that can be
#'evaluated (with the necessary arguments such as ylab etc filled in). For modifying plots.
#'@return A ggplot2 object. This can be modified by the various ggplot2 functions, or displayed using print(x).
#'@examples
#'if(w32chk()){
#'x <- ctStanDiscretePars(ctstantestfit)
#'ctStanDiscreteParsPlot(x, indices='CR')
#'
#'#to modify plot:
#'g <- ctStanDiscreteParsPlot(x, indices='CR') + 
#'  ggplot2::labs(title='My ggplot modification')
#'print(g)
#'}
#'
#'@export

ctStanDiscreteParsPlot<- function(x,indices='all',
  quantiles=c(.025,.5,.975), latentNames='auto',
  ylab='Coefficient',xlab='Time interval',ylim=NA,facets=NA,splitSubjects=TRUE,
  colour='Effect',title='Temporal regressions | independent shock of 1.0',
  polygonalpha=.1,ggcode=NA){
  
  if(is.data.frame(indices)) indices <- as.matrix(indices)
  
  title= paste0('Temporal ',ifelse(attributes(x)$cov,'covariance','regressions'),
    ' | ',ifelse(attributes(x)$observational,'correlated','independent'), ' shock of 1.0')
  
  nlatent=dim(x)[5]
  
  if(latentNames[1]=='auto') latentNames=dimnames(x)$row
  
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
  
  
  ym <- as.data.table(x)
  ym$row <- factor(ym$row)
  ym$col <- factor(ym$col)
  ym$`Time interval` <- as.numeric(ym$`Time interval`)
  
  ym$Effect <- interaction(ym$row,ym$col)
  ym$Subject <- factor(ym$Subject)
  if(!splitSubjects) ym$Subject <- factor(1)
  ym$Sample <- factor(ym$Sample)
  
  #remove rows not in indices
  ym <- ym[paste0(row,'_',col) %in% apply(indices,1,function(x) paste0(latentNames[x],collapse='_'))]
  
  # title <- paste0('Temporal ', ifelse(cov,'covariance','regressions'),' | ',
  #   ifelse(cov,'correlated','uncorrelated'), 'shock of 1.0')
  
  g<-'ggplot2::ggplot(data = ym,mapping=aes(y=value,x=`Time interval`,
    colour=Effect,
    fill=Effect,
    type=Subject))+
    theme_bw()+ylab(ylab)+
    ggplot2::labs(title = title)+  
    stat_summary( #ribbon
      fun.data = function(x) list(
        y=quantile(x,quantiles[2]),
        ymin=quantile(x,quantiles[1]), 
        ymax=quantile(x,quantiles[3])
      ),
      geom = "ribbon",
      alpha= polygonalpha,
      linetype=3)+
    stat_summary( #center line
      fun.data = function(x) list(
        y=quantile(x,quantiles[2])
      ),
      geom = "line")'
  
  if(!is.na(facets)) g <- paste0(g,'+ facet_wrap(facets)')
  
  if(!is.na(ylim)) g <- paste0(g,' + ylim(ylim)')
  
  if(is.na(ggcode)) g <- eval(parse(text=g)) else g <- list(dt=ym,ggcode=g)
  
  return(g)
}
# } else {

#   
#   #blank plot
#   rm(plot) #otherwise can't find function
#   blankargs=plotcontrol
#   blankargs$xlim=range(times)
#   blankargs$y=NA
#   blankargs$x=NA
#   do.call(plot,blankargs)
#   if(grid) {
#     grid()
#     par(new=TRUE)
#     do.call(plot,blankargs)
#     par(new=FALSE)
#   }
#   
#   ####plotting confidence region
#   if(polygon) {
#     cc=0
#     ccup=TRUE
#     for(indexi in c(1:nrow(indices))){
#       cc=ifelse(ccup,cc+1,cc-1)
#       if(indexi==nrow(indices)) ccup=FALSE
#       ri=indices[indexi,1]
#       ci=indices[indexi,2]
#       polygonargs<-polygoncontrol
#       polygonargs$x=times
#       polygonargs$y=x[ri,ci,,2]
#       polygonargs$ylow=x[ri,ci,,1]
#       polygonargs$yhigh=x[ri,ci,,length(quantiles)]
#       polygonargs$col=grDevices::adjustcolor(colvec[cc],alpha.f=max(c(.004,polygonalpha/(2*sqrt(polygonargs$steps)))))
#       do.call(ctPoly,polygonargs)
#     }
#   }
#   
#   ####plotting quantile lines
#   for(qi in 1:3){
#     cc=0
#     for(indexi in 1:nrow(indices)){
#       cc=cc+1
#       ri=indices[indexi,1]
#       ci=indices[indexi,2]
#       
#       plotargs<-plotcontrol
#       plotargs$x=times
#       plotargs$y=x[ri,ci,,qi]
#       plotargs$lty=ltyvec[cc]
#       plotargs$col=ifelse(qi!= 2,grDevices::adjustcolor(colvec[cc],alpha.f=.5) ,colvec[cc])
#       plotargs$lwd=ifelse(qi!= 2,1, lwdvec[cc])
#       do.call(points,plotargs)
#     }}
#   
#   
#   
#   legendcontrol$legend=paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
#   legendcontrol$text.col=colvec
#   legendcontrol$col=colvec
#   legendcontrol$lty = ltyvec
#   legendcontrol$lwd=lwdvec
#   
#   if(legend) do.call(graphics::legend,legendcontrol)
# }#end if not gg
# }
