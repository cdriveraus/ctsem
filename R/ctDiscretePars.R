#' ctRawParnames
#'
#' Gets internal stan parameter names of a ctStanFit object sampled via stan based on specified substrings.
#' \code{ctStanParnames} is maintained as a backward-compatible alias.
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
#' ssmodel <- ctModel(type='ct', n.latent=2, n.manifest=1,
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
#' ssfit <- ctFit(datalong, ssmodel, iter=2,
#'   optimize=FALSE, chains=1)
#' ctRawParnames(ssfit,substrings=c('pop_','popsd'))
#' }
#'
#' @aliases ctStanParnames
#' @export
ctRawParnames <- function(x,substrings=c('pop_','popsd')){
  if(length(x$stanfit$stanfit@sim)==0) stop('Doesnt contain sampled stanfit object')
  out<-c()
  for(subsi in substrings){
    out<- c(out, x$stanfit$stanfit@model_pars[grep(paste0('^',subsi),x$stanfit$stanfit@model_pars)])
  }
  return(out)
}

#' @export
ctStanParnames <- ctRawParnames



#'ctDiscretePars
#'
#'Calculate model implied regressions for a sequence of time intervals (if ct) or steps (if dt) based on
#'a ctStanFit object, for specified subjects. Wrap with print() when used inside for loops!
#'\code{ctStanDiscretePars} is maintained as a backward-compatible alias.
#'
#'@param fit model fit from \code{\link{ctFit}}
#'@param ctstanfitobj Deprecated. Use \code{fit}.
#'@param subjects Either 'popmean', to use the population mean parameter, or a vector of integers denoting which
#'subjects.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each. If the fit
#'object is a discrete time model, these should be positive integers.
#'@param nsamples Number of samples from the stanfit to use for plotting. Higher values will
#'increase smoothness / accuracy, at cost of plotting speed. Values greater than the total
#'number of samples will be set to total samples.
#'@param observational What a "one unit change in process c" is taken to bring
#'with it. Every variant is \code{dtDRIFT(t) \%*\% C} for a different companion
#'matrix \code{C}, whose column c says what else moves.
#'\describe{
#'  \item{\code{FALSE}, \code{'experimental'}}{\code{C = I}. Nothing else
#'    moves, because the change was imposed. This is the \emph{partial}
#'    regression \code{E[x(t+u)|x(t)]}, and the only variant that is a property
#'    of the dynamics alone.}
#'  \item{\code{TRUE}, \code{'observational'}}{\code{C = Sigma diag(Sigma)^-1}
#'    with \code{Sigma = asymDIFFUSIONcov}. Process c is \emph{observed} one
#'    unit above expectation and the others are not held fixed, so they move by
#'    \code{Sigma_rc/Sigma_cc}. This is the \emph{simple} regression -- what
#'    regressing the observed data on x_c alone recovers.}
#'  \item{\code{'shock'}}{\code{C = Q diag(Q)^-1} with
#'    \code{Q = DIFFUSIONcov}. One hypothetical system noise innovation arrives
#'    in c, with companion innovations \code{Q_rc/Q_cc}. An impulse response
#'    rather than a regression.}
#'  \item{\code{'orthogonal'}}{From the Cholesky factor of DIFFUSION, so
#'    shocks are uncorrelated. Depends on the order of the latent processes,
#'    which is a modelling assumption rather than a property of the fit.}
#'}
#'The companion matrices are asymmetric, and must be: a one unit move in a
#'wide-spread process implies more movement in a narrow one than the reverse.
#'\code{'observational'} and \code{'shock'} use different covariance matrices
#'because they ask different questions, and coincide only under isotropic decay
#'with no cross effects.
#'@param standardise Logical. If TRUE, output is standardised according to expected total within subject variance, given by the
#'asymDIFFUSIONcov matrix.
#'@param cov Logical. If TRUE, covariances are returned instead of regression coefficients.
#'@param plot Logical. If TRUE, ggplots output using \code{\link{ctDiscreteParsPlot}}
#'instead of returning output.
#'@param cores Number of cpu cores to use for computing subject matrices.
#'If subject matrices were saved during fiting, not used.
#'@param state Evaluation point for any matrix cell that depends on the latent
#'state or a time dependent predictor -- see \code{\link{ctContextDependence}}.
#'\code{NULL} (the default) or \code{'T0MEANS'} uses the population initial
#'state, \code{'mean'} the mean smoothed latent state, \code{'asymptotic'} the
#'system's own fixed point; or supply a numeric state vector. Requires
#'\code{backend='julia'} and \code{subjects='popmean'}. For a linear model the
#'argument makes no difference, since no cell depends on the state.
#'@param method How the interval regressions are obtained. \code{'linearise'}
#'(the default) exponentiates the DRIFT matrix, \code{expm(DRIFT*t)} -- exact
#'for a linear model, and for a model whose DRIFT depends on the state, the
#'response of the model linearised at \code{state}. \code{'simulate'}
#'integrates the nonlinear system from \code{state} and from \code{state} plus
#'one unit on each process and differences the trajectories, which is the actual
#'model implied regression; it starts at 1 and decays to 0 for a stable process,
#'and reduces exactly to \code{'linearise'} when no cell depends on the state.
#'\code{'simulate'} requires \code{backend='julia'} and is far slower, so it
#'uses \code{nsamples} posterior draws with a small default.
#'@param ... additional plotting arguments to control \code{\link{ctDiscreteParsPlot}}
#'@examples
#' data.table::setDTthreads(1) #ignore this line
#' ctDiscretePars(ctstantestfit,times=seq(.5,4,.1),
#'  plot=TRUE,indices='CR')
#'
#'#modify plot
#'require(ggplot2)
#'g=ctDiscretePars(ctstantestfit,times=seq(.5,4,.1),
#'  plot=TRUE,indices='CR')
#'g= g+ labs(title='Cross effects')
#'print(g)
#'@details
#'Two of these combinations are quantities you may already know under other
#'names. \code{observational=FALSE} is the discrete time autoregression and
#'cross-lagged matrix, \code{expm(DRIFT*t)}. \code{observational=TRUE} with
#'\code{standardise=TRUE} is the model implied \emph{cross-correlation
#'function}, \code{Cor(x_r(t+u), x_c(t))} -- the model's counterpart to the
#'empirical autocorrelations \code{\link{ctACF}} computes from the data, and
#'equal to the latent correlation matrix at \code{t=0}.
#'
#'If plot=TRUE, the function will return a ggplot2 object
#'(and hence needs to be printed if intended to display within a loop).
#'This can be modified by the various ggplot2 functions, or displayed using print(x).
#'@aliases ctStanDiscretePars
#'@export
ctDiscretePars<-function(fit, subjects='popmean',
  times=seq(from=0,to=10,by=.1),
  nsamples=200,observational=FALSE,standardise=FALSE,
  cov=FALSE, plot=FALSE,cores=1,state=NULL,method='linearise',..., ctstanfitobj){

  if(missing(fit)){
    if(missing(ctstanfitobj)) stop('fit must be supplied')
    warning('ctstanfitobj argument is deprecated, use fit instead')
    fit <- ctstanfitobj
  } else if(!missing(ctstanfitobj)) {
    stop('Use only one of fit or deprecated ctstanfitobj')
  }

  method <- match.arg(method, c('linearise','simulate'))

  # method='simulate' answers a different question and takes a different route
  # to it: integrate the nonlinear system from a state and from that state plus
  # one unit, and difference. It reduces exactly to the linearised answer when
  # no cell depends on the state, so it is safe to ask for either way.
  if(method %in% 'simulate'){
    if(!'popmean' %in% subjects) stop(call.=FALSE,
      "method='simulate' applies to subjects='popmean' only.")
    ctmS <- .ctFitModelObject(fit)
    if(!ctmS$continuoustime) stop(call.=FALSE,
      "method='simulate' is for continuous time models.")
    out <- .ctDiscreteParsSimulate(fit, times=times,
      state=if(is.null(state)) 'asymptotic' else state, nsamples=nsamples,
      observational=observational, standardise=standardise)
    times <- attr(out,'times')
    dimnames(out) <- list(Sample=seq_len(dim(out)[1]), Subject='popmean',
      `Time interval`=times, row=ctmS$latentNames, col=ctmS$latentNames)
    attributes(out)$observational <- observational
    attributes(out)$cov <- FALSE
    attributes(out)$method <- 'simulate'
    out <- .ctContextAttach(out, fit)
    if(plot) out <- ctDiscreteParsPlot(out, ...)
    return(out)
  }

  # `fit` may be a ctStanFit or a ctJuliaFit. Everything this
  # function needs is either the pop_* arrays that ctExtract() now returns for
  # all three, or model metadata; .ctFitModelObject() supplies the latter so the
  # body no longer reaches into $ctstanmodel and $standata directly.
  ctm <- .ctFitModelObject(fit)
  if(!ctm$continuoustime) times <- unique(round(times))
  type='discreteDRIFT'
  collapseSubjects=TRUE #consider this for a switch

  if(subjects[1] != 'popmean' && any(!is.integer(as.integer(subjects)))) stop('
  subjects argument must be either "popmean" or an integer denoting specific subjects')

  extractSubjects <- subjects
  if('popmean' %in% extractSubjects) extractSubjects <- 'all'
  # An explicit evaluation point only means something for the population
  # matrices: a subject's are the ones its own filter pass ended with.
  stateArgs <- list()
  if(!is.null(state)){
    if(!inherits(fit,'ctJuliaFit')) stop(call.=FALSE,
      paste0("state= requires backend='julia'. ",
        "The stan model materialises its matrices during the fit, at T0MEANS, ",
        "and cannot be asked for another point afterwards."))
    if(!'popmean' %in% subjects) stop(call.=FALSE,
      paste0("state= applies to subjects='popmean' only; a subject's matrices ",
        "come from its own filter pass."))
    stateArgs$state <- .ctResolveState(fit,state)$state
  }
  e<-do.call(ctExtract,c(list(fit,subjectMatrices = subjects[1]!='popmean',cores=cores,
    nsamples = min(nsamples,.ctFitNsamples(fit)),
    subjects=extractSubjects),stateArgs))

  nsubjects <- dim(e$subj_DRIFT)[2]
  if(is.null(nsubjects)) nsubjects=1
  if('popmean' %in% subjects) subjects='popmean'

  if(is.null(e$subj_DRIFT) && any(!subjects %in% 'popmean')) stop('No individual variation in DRIFT matrix found?? Try subjects="popmean"')

  niter=dim(e$pop_DRIFT)[1]
  latentNames=ctm$latentNames
  nlatent=length(latentNames)

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
      if(dim(e[[paste0('subj_',matname)]])[2] != .ctFitNsubjects(fit)){ #if we computed subject parameters for only the specified subjects
        parsubjects <- 1:length(subjects)
      } else parsubjects <- subjects
      ctpars[[matname]] <- e[[paste0('subj_',matname)]][samples,parsubjects,,,drop=FALSE]
      # ctpars[[matname]]<-  array(ctpars[[matname]],dim=c(prod(dim(ctpars[[matname]])[1:2]),dim(ctpars[[matname]])[-1:-2]))
    }
    # ctpars[[matname]] <- ctpars[[matname]][sample(1:dim(ctpars[[matname]])[1],nsamples),,,drop=FALSE]
  }


  # Which evaluation point the DRIFT being exponentiated below actually came
  # from -- the population pass, or each subject's own filter pass. They are
  # different points and the reader is told which one they got.
  .ctContextMessage(fit,
    if(!'popmean' %in% subjects) .ctContextSubjectLabel else
      if(is.null(state)) .ctContextPopLabel else .ctResolveState(fit,state)$label,
    paste0("expm(DRIFT*t) is therefore the transition of the model linearised ",
      "there, not the nonlinear system's own interval regression."))

  out <- ctDiscreteParsDrift(ctpars,times, observational, standardise, cov=cov,discreteInput = ctm$continuoustime==FALSE)

  dimnames(out)<- list(Sample=samples, Subject=subjects,
    `Time interval`=times, row=latentNames, col=latentNames)

  attributes(out)$observational <- observational
  attributes(out)$cov <- cov
  attributes(out)$method <- 'linearise'
  out <- .ctContextAttach(out, fit)

  if(plot) {

    out <- ctDiscreteParsPlot(out,
      # latentNames=ctstanfitobj$ctstanmodel$latentNames,
      ...)
  }
  return(out)
}

#' @export
ctStanDiscretePars <- ctDiscretePars



ctDiscreteParsDrift<-function(ctpars,times, observational,  standardise,cov=FALSE,
  types='dtDRIFT',discreteInput=FALSE, quiet=FALSE){

  nl=dim(ctpars$DRIFT)[3]

  if(!quiet) message('Computing temporal regression coefficients for ', dim(ctpars$DRIFT)[1],' samples')

  lapply(names(ctpars),function(x){ #add in extra dim if only 3 dims (e.g. when not individually varying)
    dm=dim(ctpars[[x]])
    if(length(dm)==3){
      ctpars[[x]] <<- array(ctpars[[x]],dim=c(dm[1],1,dm[2:3]))
    }
  })

  nsubs <- lapply(ctpars,function(x) dim(x)[2])


  nonstationary <- 0L
  companionType <- .ctCompanionType(observational)

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
          if(!identical(companionType,'experimental')){
            # Every variant is dtDRIFT %*% C; only the companion matrix C
            # differs. See R/ctCompanionShock.R for what each one means.
            C <- .ctCompanionMatrix(companionType,
              matrix(ctpars$DIFFUSIONcov[i,min(j,nsubs$DIFFUSIONcov),,],nl,nl),
              matrix(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,],nl,nl),
              nl)
            if(is.null(C)){
              nonstationary <- nonstationary + 1L
              ctpars$dtDRIFT[i,j,ti,,] <- NaN
              next
            }
            ctpars$dtDRIFT[i,j,ti,,] <- ctpars$dtDRIFT[i,j,ti,,] %*% C
          }

          if(standardise) {
            # Response of r per one sd_c change in c, in sd_r units: S^-1 M S,
            # with the sds of the processes themselves (asymDIFFUSION). Applied
            # after the companion step so the two compose.
            #
            # A frozen DRIFT can be non-stationary at the point it was frozen at
            # while the nonlinear system it came from is perfectly well behaved
            # -- often the point of the nonlinearity -- so this reports NaN for
            # the affected draw and names the likely cause rather than aborting
            # and blaming the model, as julia's _ctsem_asymptotics already does.
            asym <- matrix(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,],nl,nl)
            if(any(diag(asym) < 0)){
              nonstationary <- nonstationary + 1L
              ctpars$dtDRIFT[i,j,ti,,] <- NaN
              next
            }
            sdv <- sqrt(diag(asym) + 1e-10)
            ctpars$dtDRIFT[i,j,ti,,] <- ctpars$dtDRIFT[i,j,ti,,] *
              matrix(rep(sdv,each=nl) / rep(sdv,times=nl),nl)
          }
          if(cov) ctpars$dtDRIFT[i,j,ti,,]  <- tcrossprod(ctpars$dtDRIFT[i,j,ti,,] )
        }
      }
    }
  } #end dtdrift

  if(nonstationary > 0 && !quiet) message(
    nonstationary, ' of ', length(ctpars$dtDRIFT)/prod(dim(ctpars$DRIFT)[3:4]),
    ' standardisations returned NaN: the asymptotic diffusion has negative ',
    'diagonals, so there is no stationary variance to standardise by. For a ',
    'model whose DRIFT depends on the state this usually reflects the point it ',
    'was linearised at rather than the model -- try another state, or ',
    'standardise=FALSE.')

  return(ctpars$dtDRIFT)
}

ctStanDiscreteParsDrift <- ctDiscreteParsDrift

#'ctDiscreteParsPlot
#'
#'Plots the output from \code{\link{ctDiscretePars}}, for model implied regression strengths at specified times for continuous time models fit with ctStanFit.
#'\code{ctStanDiscreteParsPlot} is maintained as a backward-compatible alias.
#'
#'@param x list object returned from \code{\link{ctDiscretePars}}.
#'@param indices Either a string specifying type of plot to create, or an n by 2
#'matrix specifying which indices of the output matrix to plot.
#''AR' specifies all diagonals, for discrete time autoregression parameters.
#''CR' specifies all off-diagonals,for discrete time cross regression parameters.
#''all' plots all AR and CR effects at once.
#'@param quantiles numeric vector of length 3, with values between 0 and 1, specifying which quantiles to plot.
#'The default plots 95\% credible intervals and the posterior median at 50\%.
#'@param latentNames Vector of character strings denoting names for the latent variables.
#''auto' just uses eta1 eta2 etc.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha of
#'the fill.
#'@param ylab y label.
#'@param xlab x label.
#'@param ylim Custom ylim.
#'@param facets May be 'Subject' or 'Effect'.
#'@param colour Character string denoting how colour varies. 'Effect' or 'Subject'.
#'@param title Character string. 'auto' generates automatically, NULL can be used to disable title.
#'@param splitSubjects if TRUE, subjects are plotted separately, if FALSE they are combined.
#'@param ggcode if TRUE, returns a list containing the data.table to plot, and a character string that can be
#'evaluated (with the necessary arguments such as ylab etc filled in). For modifying plots.
#'@return A ggplot2 object. This can be modified by the various ggplot2 functions, or displayed using print(x).
#'@examples
#' data.table::setDTthreads(1) #ignore this line
#'x <- ctDiscretePars(ctstantestfit)
#'ctDiscreteParsPlot(x, indices='CR')
#'
#'#to modify plot:
#'g <- ctDiscreteParsPlot(x, indices='CR') +
#'  ggplot2::labs(title='My ggplot modification')
#'print(g)
#'
#'@aliases ctStanDiscreteParsPlot
#'@export

ctDiscreteParsPlot<- function(x,indices='all',
  quantiles=c(.025,.5,.975), latentNames='auto',
  ylab='Coefficient',xlab='Time interval',ylim=NA,facets=NA,splitSubjects=TRUE,
  colour='Effect',title='auto',
  polygonalpha=.1,ggcode=FALSE){

  if(is.data.frame(indices)) indices <- as.matrix(indices)

  if(!is.null(title)){
    if(title %in% 'auto'){
      title= paste0('Temporal ',ifelse(attributes(x)$cov,'covariance','regressions'),
        ' | ',ifelse(attributes(x)$observational,'correlated','independent'), ' shock of 1.0')
    }
  }else title=title

  nlatent=dim(x)[5]

  if(latentNames[1]=='auto') latentNames=dimnames(x)$row
  dimnames(x)$row <- dimnames(x)$col <- latentNames

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
      rep(unique(indices),each=nlatent))
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

  g<-paste0('ggplot2::ggplot(data = ym,mapping=aes(y=value,x=`Time interval`,
    colour=Effect,group=interaction(Effect, Subject),
    fill=Effect))+
    theme_bw()+
    ylab("',ylab,'")+
    xlab("',xlab,'")+
    ggplot2::labs(title = "',title,'")+
    stat_summary( #ribbon
      fun.data = function(x) list(
        y=quantile(x,',quantiles[2],'),
        ymin=quantile(x,',quantiles[1],'),
        ymax=quantile(x,',quantiles[3],')
      ),
      geom = "ribbon",
      alpha= ',polygonalpha,',
      linetype=3)+
    stat_summary( #center line
      fun.data = function(x) list(y=quantile(x,',quantiles[2],')),
      geom = "line")')

  #needs data input as well...
  # if(rug) g <- paste0(g,'+ geom_rug(data= as.data.table(hist(data.table(ym)[order(Subject,`Time interval`),][,.(timeDiff=c(diff(`Time interval`))),by=Subject]$timeDiff,plot=F)[2:4])[counts > 0,],
  #  aes(x = mids,alpha=density,size=density/max(density)),inherit.aes=F,linetype=1,show.legend = FALSE)')

  if(!is.na(facets)) g <- paste0(g,'+ facet_wrap(facets)')

  if(!is.na(ylim)) g <- paste0(g,' + ylim(ylim)')

  if(!ggcode) g <- eval(parse(text=g)) else g <- list(dt=ym,ggcode=g)

  return(g)
}

#' @export
ctStanDiscreteParsPlot <- ctDiscreteParsPlot
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
