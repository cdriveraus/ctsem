ctExtractComplexMatrix <- function(matname, fit, time,id){
  ks=ctStanKalman(fit)
  if(!id %in% ks$id) stop('Specified id not found in dataset')
  if(!time %in% ks$time[id %in% ks$id]) stop(paste0('Specified time not found in data from subject', id))
  s=summary(fit,parmatrices=FALSE,residualcov=FALSE,priorcheck=FALSE)
  list2env(data.frame(t(s$popmeans[,'50%',drop=FALSE])),envir = environment())
  list2env(data.frame(t(ks$etasmooth[1,which(ks$time %in% time & ks$id %in% id)[1], ])),envir = environment())
  mat <- listOfMatrices(fit$ctstanmodelbase$pars)[[matname]]
  fit$ctstanmodel$LAMBDA[,1,drop=FALSE]
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j] <- eval(parse(text=paste0(mat[i,j])))
    }
  }
  mat <- matrix(as.numeric(mat),nrow(mat))
}
