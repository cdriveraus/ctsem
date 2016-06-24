#' Fits an individually varying continuous time model.
#' 
#' 
#' @param datawide Wide format data, as used in \code{\link{ctFit}}.  See \code{\link{ctLongToWide}} to
#' easily convert long format data.
#' @param ctmodelobj Continuous time model to fit, specified via \code{\link{ctModel}} function.
#' @param fixedmodel Modified version of ctmodelobj, wherein any parameters you wish to keep 
#' fixed over individuals should be given the value 'groupfixed'.  
#' If specified, all other parameters will be free across groups.
#' @param freemodel Modified version of ctmodelobj, wherein any parameters you wish to free across individuals
#' should be given the label 'groupfree'.  
#' If specified, all other parameters will be fixed across individuals  
#' If left NULL, the default, all parameters are free across individuals
#' @param confidenceintervals Character vector of parameter labels to estimate confidence intervals for.  
#' Unlike with \code{\link{ctFit}}, entire matrices cannot be specified.
#' @param showInits Displays start values prior to optimization
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to encourage parameters to remain closer to 0, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help with optimization in some cases, though results in user specified inits being ignored for the final fit.
#' @param retryattempts Number of fit retries to make.
#' @param plotOptimization plots graphs of optimization progress after fitting, uses OpenMx checkpointing.
#' @param ... additional arguments to pass to \code{\link{ctFit}}.
#' @return Returns either an OpenMx fit object (if fitasgroup=TRUE), or a list containing 
#' individual ctsem fit objects.
#' @details Additional \code{\link{ctFit}} parameters may be specified as required.
#' @seealso \code{\link{ctFit}} and \code{\link{ctModel}}
#' @export


ctIndividualFit<-function(datawide,ctmodelobj,hyperpars='estimate',
  meanParams='estimate',optimizeMeans=TRUE, optimizePenalties=TRUE,roughPass=TRUE,
  OOSiterations=4,percentmissing=0.1,threads=4,
  retryattempts=5,...){
  
  flexLapply<-function(cluster, object, fn){
    
    if(threads > 1) out<-parLapply(cluster, object, fn)
    if(threads == 1) out<-lapply(object, fn)
    return(out)
  }
  
  library(optimx)
  library(parallel)
  library(OpenMx)
  
  if(threads > 1) CL=makeCluster(threads)
  
  n.subjects<-nrow(datawide)
  manifestNames<-ctmodelobj$manifestNames
  
 
  ########base models and data
  model<-ctFit(datawide[1,,drop=F],ctmodelobj,discreteTime=F,
    asymptotes=T,stationary=c('T0VAR','T0MEANS','T0TRAITEFFECT','T0TIPREDEFFECT'))
  if(threads > 1) clusterExport(cl = CL, c("datawide", "ctmodelobj", 'n.subjects'), envir=environment()) 
  base<- flexLapply(CL, 1:n.subjects, function(x) {
    library(ctsem)
    model<-ctFit(datawide[x,,drop=F],ctmodelobj,discreteTime=F,
      asymptotes=T,stationary=c('T0VAR','T0MEANS','T0TRAITEFFECT','T0TIPREDEFFECT'))
    params<-omxGetParameters(model$mxobj)
    data<-model$mxobj$data$observed
    out<-list(params,data)
  }
  )
  
  data<-lapply(1:nrow(datawide), function(x) base[[x]][[2]])
  indParams<-lapply(1:nrow(datawide), function(x) base[[x]][[1]])

  indParams<-matrix(unlist(indParams),byrow=T,ncol=length(indParams[[1]]), dimnames=list(c(),names(indParams[[1]])))
  if(meanParams == 'estimate') {
    meanParams <- apply(indParams,2,median,na.rm=T)
  }
  names(meanParams) <- colnames(indParams)
  
  
  ##### plot base densities
  # originalmfrow<-par()$mfrow
#   par(mfrow=c(3,3))
#   lapply(1:length(meanParams), function(x) plot(density(indParams[,x]),main=names(meanParams)[x]))
#   
#   
  
  
  ######penalised model 
#   if(threads > 1) clusterExport(cl = CL, c("meanParams",'indParams'), envir=environment()) 
#   fullmodels<-flexLapply(CL, 1:nrow(datawide), function(x) {
#     library(OpenMx)

    fullmodel<-mxModel(paste0('i',1,'_'), 
      
      model$mxobj,
      
      mxMatrix(name='meanParams', 
        values=meanParams, 
        labels=paste0('mean_',names(meanParams)),
        free=F,nrow=1,ncol=length(meanParams)),
      
      mxMatrix(name='indParams', 
        labels=paste0(names(meanParams)),
        values=meanParams, #indParams[x,], #
        free=T,nrow=1,ncol=length(meanParams)),
      
      mxMatrix(name='hyperpars',
        values=rep(1,length(meanParams)),
        labels=paste0('penalty_',names(meanParams)),
        free=F, nrow=1, ncol=length(meanParams)),
      
      mxAlgebra(name='algFit', ctsem.objective + sum(
        (-indParams + meanParams) * (-indParams + meanParams) * exp(hyperpars))),
      
      mxFitFunctionAlgebra('algFit')
    )
    fullmodel<-omxAssignFirstParameters(fullmodel)
  
  #####create list of missings to use
  
  missingsList<-lapply(1:nrow(datawide), function(i) {
    lapply(1:OOSiterations,function(x){
      observedIndices<-which(!is.na(data[[i]][,manifestNames]))
      sample(observedIndices, max(1,floor(percentmissing * length(observedIndices))))
    })
  })
  
  
  
  if(hyperpars=='estimate') hyperpars<- rep(-1,length(meanParams))
  if(threads > 1) clusterExport(cl = CL, c('fullmodel','data', 'retryattempts',"meanParams", 'hyperpars',
    'missingsList', 'OOSiterations', 'manifestNames','n.subjects'), envir=environment())
  

  
  
  ####### out of sample fit function
  OOSfit<-function(params,fullmodel, data,missingsList, OOSiterations,optimizeMeans,optimizePenalties,rough,...){

    if(optimizePenalties==TRUE && optimizeMeans==TRUE) {
      hyperpars<-params[1:(length(params)/2)]
      meanParams<-params[((length(params)/2)+1):length(params)]
    }
    if(optimizePenalties==FALSE && optimizeMeans==TRUE) meanParams<-params
    
    if(optimizePenalties==TRUE && optimizeMeans==FALSE) hyperpars <- params
    
    manifestNames<-rownames(fullmodel$ctsem[[fullmodel$ctsem$expectation$C]]$labels)
    
    if(threads > 1) clusterExport(cl = CL, c('retryattempts',"meanParams", 'hyperpars',
      'OOSiterations', 'manifestNames','n.subjects'), envir=environment()) 
    
    #   if(meanParams!='estimate'){
    #   ###calculate mean params for new hyperparameter
    #   params<-flexLapply(CL,1:n.subjects, function(x) {
    #     library(OpenMx)
    #     library(ctsem)
    #     fullmodel<-fullmodels[[x]]
    #     params<-list()
    #   
    #       fullmodel$hyperpars$values<-hyperpars
    #       fullmodel<-mxModel(fullmodel,
    #         mxMatrix(name='meanParams', 
    #           values=meanParams, 
    #           labels=paste0('mean_',names(meanParams)),
    #           free=F,nrow=1,ncol=length(meanParams))
    #       )
    #       
    #       fullmodel<-suppressMessages(ctmxTryHard(fullmodel, initialTolerance=1e-18, checkHess=TRUE, greenOK=FALSE, 
    #         iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=0,
    #         extraTries=retryattempts, loc=1, scale=2, paste=FALSE))
    #       params<-omxGetParameters(fullmodel)
    #     })
    #     meanParams<-colMeans(matrix(unlist(params),byrow=T,ncol=length(params[[1]]),dimnames=list(c(),names(params[[1]]))),na.rm=T)
    #     ##### plot new densities over base
    #     #   lapply(1:length(meanParams), function(x) {
    #     #     plot(density(indParams[,x]),main=names(meanParams)[x])
    #     #     points(density(matrix(unlist(params),byrow=T,ncol=length(params[[1]]))[,x]),main=names(meanParams)[x],type='l',col='red')
    #     #   })
    #   }
    
    
    print(hyperpars)
    print(meanParams)
    
    
    
    
    ### calculate OOS error for new hyperparameter
    OOSerror<-flexLapply(CL, 1:n.subjects, function(x) {

      
      m2ll<-c()
      for( OOSit in 1:OOSiterations){
        retry<-TRUE
        while(retry==TRUE){
          retry<-FALSE
        model<-fullmodel
        model$ctsem$data$observed<-data[[x]]
        model$ctsem$data$observed[1,'dT1']<-.05
        missings<-rep(FALSE,length(model$ctsem$data$observed[,manifestNames]))
        missings[missingsList[[x]][[OOSit]]] <- TRUE
        fulldata<-model$ctsem$data$observed
        model$ctsem$data$observed[,manifestNames][missings] <- NA
        model$hyperpars$values<-hyperpars
        model$meanParams$values <- meanParams
        model$indParams$values <- meanParams
        model<-omxAssignFirstParameters(model)
        
#         model<-mxModel(model,
#           mxMatrix(name='meanParams', 
#             values=meanParams, 
#             labels=paste0('meanParam_',1:length(meanParams)),
#             free=F,nrow=1,ncol=length(meanParams))
        # )
        
        model<-suppressMessages(ctmxTryHard(model, initialTolerance=1e-18, checkHess=TRUE, greenOK=FALSE, 
          iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=0,
          extraTries=30, loc=1, scale=2, paste=FALSE))
        
        if(class(model) =='try-error' || is.na(model$output$fit)) {
          m2ll <- NA
          return(m2ll)
        }
        else{
        

        print(x)
        # OOSpredictions<-model$output$fit # try(kalmanPredict(model$ctsem,mxstyle=T))
        # OOS<-model
        # model<-OOS
        model$ctsem$data$observed <-fulldata
        model$ctsem$data$observed[,manifestNames][!missings] <- NA
        # model<-mxModel(model, mxComputeOnce('fitfunction', 'fit'))
#         model$ctsem<-mxModel(model$ctsem, mxComputeOnce('fitfunction', 'fit'))
        # model$ctsem<-mxRun(model$ctsem)

        model<-omxSetParameters(model,names(omxGetParameters(model)),free=FALSE)
        model$ctsem<-mxModel(model$ctsem,mxMatrix(name='tempdummy',free=T,nrow=1,ncol=1,values=999))
        model<-suppressMessages(ctmxTryHard(model$ctsem, initialTolerance=1e-18, checkHess=FALSE, greenOK=TRUE, 
          iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=0,
          extraTries=0, loc=1, scale=2, paste=FALSE))
        
        # ISpredictions<-model$output$fit #try(kalmanPredict(model$ctsem,mxstyle=T))
        # mse<- c(mse, mean( (predictions[missings] - models[[x]]$ctsem$data$observed[,manifestNames][missings])^2,na.rm=T) )
        # m2ll<-ISpredictions - OOSpredictions
        
        if(is.na(model$output$fit) && retry == FALSE) retry<-TRUE
        }
      }
      }
  
      if(!exists('m2ll')) m2ll<-NA
      # m2ll<-sum(m2ll)
      m2ll <- c(m2ll,model$output$fit)
    })
   
    if(any(is.na(OOSerror))) message(length(is.na(unlist(OOSerror))==TRUE), ' missing likelihoods!')
    OOSerror <- mean(unlist(OOSerror),na.rm=T)
    print(OOSerror)
    return(OOSerror)
  }
  
  
#   mxOOSfit<-function(model,state,...){
#    OOSfit(...) 
#   }
#   
#   model<-mxModel('hypermodel',
#     
#     mxFitFunctionR(mxOOSfit
  
  ####optimize commands
  
  if(optimizePenalties==TRUE || optimizeMeans==TRUE) { #then optimize
    
    #set bounds for hyper and mean params
    lbound<-c()
    params<-c()
    ndeps<-c()
    fitdifference <- Inf

    # while(fitdifference > .000001){
      

      
    
    if(optimizePenalties==TRUE){  
      params <- hyperpars
      lbound<-rep(0,length(meanParams))
      ndeps<-rep(0.1,length(meanParams))
      
#       fit<-optimx(params,OOSfit,fullmodels=fullmodels,
#         lower=lbound,optimizeMeans=F,optimizePenalties=T,
#         missingsList=missingsList,
#         # method='bobyqa',
#         # rough=TRUE,
#         OOSiterations=OOSiterations,
#         control=list(reltol=1e-3, ndeps=ndeps,
# #           rhobeg=10,
# #           rhoend=.01,
# #           kkt=F, 
# #           # fnscale=ndeps*10, 
# #           # factr=1e-1, 
#           lmm=50,
#           abstol=.001))
      
    }
      
      if(optimizeMeans==TRUE) {
        params <- c(params,meanParams)
        lbound<-c(lbound,rep(-9999999,length(meanParams)))
        ndeps<-c(ndeps,rep(0.01,length(meanParams)))
        
        #         fit<-optimx(meanParams,OOSfit,fullmodels=fullmodels,
        #           # lower=lbound,
        #           optimizeMeans=T,optimizePenalties=F,
        #           missingsList=missingsList,
        #           method='L-BFGS-B',
        #           # rough=TRUE,
        #           OOSiterations=OOSiterations,
        #           control=list(reltol=1e-3, ndeps=rep(.001,length(meanParams)),
        #             #           rhobeg=10,
        #             #           rhoend=.01,
        #             #           kkt=F, 
        #             #           # fnscale=ndeps*10, 
        #             #           # factr=1e-1, 
        #                       lmm=50,
        #             abstol=.0001))
        
        
      }
      
   
    
    
    ##rough optimize hyper and mean params
    if(roughPass==TRUE){
      
      fit<-optimx(params,OOSfit,fullmodels=fullmodels,
        # lower=lbound,
        missingsList=missingsList,
        # method='bobyqa',
        rough=TRUE,
        OOSiterations=1,
        control=list(reltol=1e-2, ndeps=ndeps,
          rhobeg=10,
          rhoend=.01,
          kkt=F, 
          # fnscale=ndeps*10, 
          # factr=1e-1, 
          lmm=50,
          abstol=.1))
      
      

      
      params<-unlist(lapply(1: length(params),  function(x) fit[[x]][1]))
    }
    
    
#     negOOSfit<-function(params=params,...){
#       out<- -.5*OOSfit(params,...)
#     }
#     OOSMCMC<-MCMC
#     environment(OOSMCMC)<-environment()
#     
#     fit<-OOSMCMC(params,negOOSfit,
#       fullmodels=fullmodels,OOSiterations=OOSiterations,rough=F,missingsList=missingsList,
#       optimizeMeans=optimizeMeans, optimizePenalties=optimizePenalties,
#       mcmcthreads=1,plotfreq=10)
#     
    
    #fine optimize
    fit<-optimx(params,OOSfit,fullmodel=fullmodel,data=data,optimizeMeans=optimizeMeans, optimizePenalties=optimizePenalties,
      # lower=lbound,
      missingsList=missingsList,
      method='BFGS',
      # method='Nelder-Mead',
      # method='BFGS',
      # method='newuoa',
      rough=FALSE,
      OOSiterations=OOSiterations,
    
      control=list(reltol=1e-6,
        pgtol=1e-6,
        # ndeps=ndeps,
        rhobeg=.1,
        rhoend=.0001,
        npt=length(params*2+1),
        
#         kkt=F, 
        parscale=ndeps*100, 
        factr=1e10, 
        lmm=10,abstol=.00001))
    
    #extract hyper and mean params from fit
    if(optimizePenalties==TRUE){  
      hyperpars<- unlist(lapply(1:ifelse(optimizeMeans==TRUE, length(params)/2, length(params)), 
        function(x) fit[[x]][1]))
      if(any(is.na(hyperpars))) browser()
    }
    if(optimizeMeans==TRUE) {
      meanParams<- unlist(lapply(ifelse(optimizePenalties==TRUE, length(params)/2+1, 1):length(params), 
        function(x) fit[[x]][1]))
      if(any(is.na(meanParams))) browser()
    }
  }
  
  
  
  
  #######final model fit
  if(threads > 1) clusterExport(cl = CL, c("hyperpars", 'meanParams'), envir=environment()) 
  fits<-flexLapply(CL,1:n.subjects, function(x) {

    fullmodel$ctsem$data$observed<-data[[x]]
    params<-list()
    
    fullmodel$hyperpars$values<-hyperpars
    fullmodel$meanParams$values <- meanParams
    fullmodel$indParams$values <- meanParams
    fullmodel<-omxAssignFirstParameters(fullmodel)
#     fullmodel<-mxModel(fullmodel,
#       mxMatrix(name='meanParams', 
#         values=meanParams, 
#         labels=paste0('mean_',1:length(meanParams)),
#         free=F,nrow=1,ncol=length(meanParams))
    # )
    
    
    fullmodel<-suppressMessages(ctmxTryHard(fullmodel, initialTolerance=1e-18, checkHess=TRUE, greenOK=FALSE, 
      iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=0,
      extraTries=retryattempts, loc=1, scale=2, paste=FALSE))
  })
  
  freeParams<-matrix(unlist(lapply(1:n.subjects, function(x) omxGetParameters(fits[[x]]))),nrow=n.subjects)
  colnames(freeParams) <- names(omxGetParameters(fits[[1]]))
  
out<-list(freeParams=freeParams,fits=fits, data=data, hyperParams=hyperpars, meanParams=meanParams)

  if(threads > 1)  stopCluster(CL)
  return(out)
}


