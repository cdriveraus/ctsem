#'hyperparfit
#'@export



hyperparFit<-function(hyperpars,hyperparnames,mxobj,OOSiterations,submodels,
  missingsList, OOSthreads=1,...){
  
  manifestNames<-rownames(mxobj[[
    paste0('penalised_',submodels[1])]][[
      submodels[1]]][[
        mxobj[[
          paste0('penalised_',submodels[1])]][[
            submodels[1]]]$expectation$C]]$labels)
  
  names(hyperpars)<-hyperparnames

  mxobj<-omxSetParameters(mxobj, labels=mxobj$hyperpars$labels,free=FALSE, values=hyperpars)
  
  mxobj<-ctmxTryHard(multif,initialTolerance=1e-10,verbose=0,
    showInits=T,
    bestInitsOutput=FALSE,
    extraTries=0,loc=1,scale=.2,paste=FALSE,iterationSummary=TRUE)

   CL = makeCluster(OOSthreads)
  clusterExport(cl = CL, c('mxobj', 
    'submodels', 'missingsList'), envir=environment()) 
  parLapply(CL,1:OOSthreads,function(x) {
    library(ctsem)
    library(OpenMx)
    NULL
  })

  fits <- parLapply(CL,1:length(submodels), function(x) {

      submodel<-submodels[x]
    fit <- kalmanOOS(mxobj=mxobj[[paste0('penalised_',submodel)]][[submodel]],nmissing=nmissing,
      missingsList=missingsList[[x]])
  })
  
  stopCluster(CL)
    fits<-unlist(fits)
    
    print(hyperpars)
    # message(paste0('fits = ', paste(fits,sep=', ')))
    # message('fit variance = ', paste(sd(c(fits))^2))
    message('overall fit = ', mean(fits))
    return(mean(fits,na.rm=T))
}





