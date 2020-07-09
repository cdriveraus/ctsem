#' Returns population system matrices from a ctStanFit object, and vector of values for free parameters.
#'
#' @param fit ctStanFit object.
#' @param parvalues vector of parameter values to assign to free parameters in the model
#' @param timeinterval time interval to use for discrete time (dt) matrix calculations.
#' @param sf stanfit object. Generally not necessary, but for repeated calls to this function, can speed things up.
#'
#' @return A list containing various matrices related to a continuous time dynamic model. 
#' Matrices with "dt" in front refers to discrete time, "asym" refers to asymptotic (time interval = infinity), 
#' and "cor" refers to correlations. 
#' @export
#'
#' @examples
#' if(w32chk()){
#'
#' ctStanParMatrices(ctstantestfit,
#'   rnorm(length(ctstantestfit$stanfit$rawest),0,.1))
#' }
ctStanParMatrices <- function(fit, parvalues, timeinterval=1, sf=NA){
  
  if(!'ctStanFit' %in% class(fit)) stop('not a ctStanFit object')

  model <- fit$ctstanmodel
  fit$standata$savescores <- 0L
  fit$standata$gendata <- 0L
  fit$standata$dokalman <- 0L
  nlatent = fit$standata$nlatent
  # browser()
  # if(length(parvalues)!=fit$data$nparams) stop('length of parvalues != number of free params (',fit$data$nparams,') in model!')
  if(suppressWarnings(is.na(sf))) sf <- stan_reinitsf(fit$stanmodel,data=fit$standata) #suppressOutput(sf <- suppressWarnings(sampling(,iter=1,control=list(max_treedepth=1),chains=1)))
  # npars <- get_num_upars(sf)
  # pars <- c(parvalues,rep(0,npars - fit$data$nparams))
  sfc <- constrain_pars(sf, parvalues)
  
  
   whichmatrices='all'
  
  if(whichmatrices[1] == 'all') {
    whichmatrices <- sort(c(unique(fit$ctstanmodel$pars$matrix),'DIFFUSIONcov','DIFFUSIONcor','asymDIFFUSION','asymDIFFUSIONcor','T0VARcor','asymCINT'))
    if(fit$ctstanmodel$continuoustime) whichmatrices <- c(whichmatrices, 
      'dtDRIFT','dtDIFFUSION','dtCINT')
  } else whichmatrices <- unique(c(whichmatrices, fit$setup$matrices$base)) #need base matrices for computations
  
  # stanmats <- c(fit$setup$matrices$base,'asymDIFFUSION')[c(fit$setup$matrices$base,'asymDIFFUSION') %in% whichmatrices]
  
  out <- list()
  for(m in whichmatrices){
    # assign(m,ctCollapse(sfc[[paste0('pop_',m)]],1,mean)) 
    if(!is.null(sfc[[paste0('pop_',m)]])) out[[m]] <- sfc[[paste0('pop_',m)]]
  }

  #because of intoverpop
  out$DRIFT <- out$DRIFT[1:nlatent,1:nlatent,drop=FALSE]
  out$T0VAR <- out$T0VAR[1:nlatent,1:nlatent,drop=FALSE]
  out$T0MEANS <- out$T0MEANS[1:nlatent,,drop=FALSE]
  
  
  #cholesky factor fix
  out$MANIFESTVAR=out$MANIFESTVAR %*% t(out$MANIFESTVAR) #cholesky factor inside stanfit...
  
  #dimension naming (latent row object, manifest column object, etc
  for(lro in c('DRIFT','DIFFUSION','CINT','T0VAR','T0MEANS','asymDIFFUSION',if('TDPREDEFFECT' %in% model$pars$matrix) 'TDPREDEFFECT')){
    if(lro %in% whichmatrices) rownames(out[[lro]]) <- model$latentNames
  }
  for(lco in c('DRIFT','DIFFUSION','T0VAR','asymDIFFUSION','LAMBDA')){
    if(lco %in% whichmatrices) colnames(out[[lco]]) <- model$latentNames
  }
  for(mro in c('LAMBDA','MANIFESTVAR','MANIFESTMEANS')){
    if(mro %in% whichmatrices) rownames(out[[mro]]) <- model$manifestNames
  }
  for(mco in c('MANIFESTVAR')){
    if(mco %in% whichmatrices) colnames(out[[mco]]) <- model$manifestNames
  }

  if('TDPREDEFFECT' %in% model$pars$matrix)  colnames(out$TDPREDEFFECT) <- model$TDpredNames
  
  
  choltrue <- FALSE #!as.logical(fit$data$lineardynamics)
  
  # if(choltrue) DIFFUSION = msquare(DIFFUSION) #sdcovchol2cov(DIFFUSION,0)
  if('DIFFUSIONcor' %in% whichmatrices){
    out$DIFFUSIONcor = suppressWarnings(stats::cov2cor(out$DIFFUSIONcov))
    out$DIFFUSIONcor[is.na(out$DIFFUSIONcor)] <- 0
  }

  if('asymDIFFUSIONcor' %in% whichmatrices){
    out$asymDIFFUSIONcor = suppressWarnings(stats::cov2cor(out$asymDIFFUSION))
    out$asymDIFFUSIONcor[is.na(out$asymDIFFUSIONcor)] <- 0
  }

  if(fit$ctstanmodel$continuoustime) out$dtDRIFT=as.matrix(Matrix::expm(out$DRIFT * timeinterval))
  if('dtDIFFUSION' %in% whichmatrices) out$dtDIFFUSION = out$asymDIFFUSION - (out$dtDRIFT %*% out$asymDIFFUSION %*% t(out$dtDRIFT ))
  if('dtDIFFUSIONcor' %in% whichmatrices) out$dtDIFFUSIONcor = cov2cor(out$dtDIFFUSION)
  if('dtCINT' %in% whichmatrices) out$dtCINT = (solve(out$DRIFT, out$dtDRIFT - diag(nrow(out$DRIFT))) %*% (out$CINT))
  if('asymCINT' %in% whichmatrices) out$asymCINT = matrix(out$asymCINT,ncol=1)#-solve(out$DRIFT) %*% out$CINT
  
  if('T0VARcor' %in% whichmatrices) {
    out$T0VARcor = suppressWarnings(stats::cov2cor(out$T0VAR))
    out$T0VARcor[is.na(out$T0VARcor)] <- 0
  }

  return(out)
}
