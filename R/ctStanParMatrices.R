
ctStanParMatrices <- function(fit, parvalues, parindex, timeinterval=1){
  if(!'ctStanFit' %in% class(fit)) stop('not a ctStanFit object')
  
  model <- fit$ctstanmodel
  fit$standata$savescores <- 0L
  fit$standata$gendata <- 0L
  fit$standata$dokalman <- 0L
  fit$standata$dokalmanrows[-1] <- 0L
  nlatent = fit$standata$nlatent

  whichmatrices <- c('PARS','T0MEANS','LAMBDA','DRIFT','MANIFESTMEANS','CINT',
    'DIFFUSIONcov','DIFFUSIONcor','asymDIFFUSIONcov','asymDIFFUSIONcor','T0cov','T0cor','asymCINT','MANIFESTcov')
  if(fit$ctstanmodel$continuoustime) whichmatrices <- c(whichmatrices, 
    'dtDRIFT','dtDIFFUSION','dtCINT')
  
  out <- list()
  for(m in whichmatrices){
    if(!is.null(parvalues[[paste0('pop_',m)]])){
      arrcommas <- paste0(rep(',',times=length(dim(parvalues[[paste0('pop_',m)]]))-1),collapse='')
      out[[m]] <- 
        eval(parse(text=paste0(
          "array(parvalues[[paste0('pop_',m)]][parindex",arrcommas,"],dim=dim(parvalues[[paste0('pop_',m)]])[-1])")))
    }
  }
  
  out$T0cov <- out$T0cov[1:nlatent,1:nlatent,drop=FALSE]
  out$T0MEANS <- out$T0MEANS[1:nlatent,,drop=FALSE]
  
  
  # #cholesky factor fix
  # out$MANIFESTVAR=out$MANIFESTcov #cholesky factor inside stanfit...
  # out$DIFFUSION=out$DIFFUSIONcov
  
  # #dimension naming (latent row object, manifest column object, etc
  # for(lro in c('DRIFT','DIFFUSION','CINT','T0VAR','T0MEANS','asymDIFFUSION',if('TDPREDEFFECT' %in% model$pars$matrix) 'TDPREDEFFECT')){
  #   if(lro %in% whichmatrices) rownames(out[[lro]]) <- model$latentNames
  # }
  # for(lco in c('DRIFT','DIFFUSIONcov','T0VAR','asymDIFFUSION','LAMBDA')){
  #   if(lco %in% whichmatrices) colnames(out[[lco]]) <- model$latentNames
  # }
  # for(mro in c('LAMBDA','MANIFESTVAR','MANIFESTMEANS')){
  #   if(mro %in% whichmatrices) rownames(out[[mro]]) <- model$manifestNames
  # }
  # for(mco in c('MANIFESTVAR')){
  #   if(mco %in% whichmatrices) colnames(out[[mco]]) <- model$manifestNames
  # }
  # 
  # if('TDPREDEFFECT' %in% model$pars$matrix)  colnames(out$TDPREDEFFECT) <- model$TDpredNames
  # 
  
  # choltrue <- FALSE #!as.logical(fit$data$lineardynamics)
  
  # if(choltrue) DIFFUSION = msquare(DIFFUSION) #sdcovchol2cov(DIFFUSION,0)
  if('DIFFUSIONcor' %in% whichmatrices){
    out$DIFFUSIONcor = suppressWarnings(stats::cov2cor(out$DIFFUSIONcov))
    out$DIFFUSIONcor[is.na(out$DIFFUSIONcor)] <- 0
  }
  
  if('asymDIFFUSIONcor' %in% whichmatrices){
    out$asymDIFFUSIONcor = suppressWarnings(stats::cov2cor(out$asymDIFFUSIONcov))
    out$asymDIFFUSIONcor[is.na(out$asymDIFFUSIONcor)] <- 0
  }
  
  if(fit$ctstanmodel$continuoustime) out$dtDRIFT=as.matrix(Matrix::expm(out$DRIFT * timeinterval))
  if('dtDIFFUSION' %in% whichmatrices) out$dtDIFFUSION = out$asymDIFFUSIONcov - (out$dtDRIFT %*% out$asymDIFFUSIONcov %*% t(out$dtDRIFT ))
  if('dtDIFFUSIONcor' %in% whichmatrices) out$dtDIFFUSIONcor = cov2cor(out$dtDIFFUSION)
  if('dtCINT' %in% whichmatrices) out$dtCINT = (solve(out$DRIFT, out$dtDRIFT - diag(nrow(out$DRIFT))) %*% (out$CINT))
  if('asymCINT' %in% whichmatrices) out$asymCINT = matrix(out$asymCINT,ncol=1)#-solve(out$DRIFT) %*% out$CINT
  
  if('T0cor' %in% whichmatrices) {
    out$T0cor = suppressWarnings(stats::cov2cor(out$T0cov))
    out$T0cor[is.na(out$T0cor)] <- 0
  }
  
  return(out)
}
