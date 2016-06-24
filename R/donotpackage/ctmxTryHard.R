ctmxTryHard<-function (model, extraTries = 10, greenOK = FALSE, loc = 1, 
  scale = 0.25,  initialGradientStepSize = .00001, initialGradientIterations = 1,
  initialTolerance=1e-12, checkHess = TRUE, fit2beat = Inf, paste = TRUE,
  iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals=FALSE){
  
    defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
  
  lastNoError<-TRUE
  
  generalTolerance<- 1e-5 #used for lowest min check
  
  gradientStepSize<- initialGradientStepSize
  tolerance <- initialTolerance
  lastBestFitCount<-0 #number of consecutive improvements in fit
  gradientIterations <- initialGradientIterations  
  stopflag <- FALSE #should the iterative optimization process stop
  numdone <- 0
  lowestminsofar<-Inf 
  inits<-omxGetParameters(model) 
  
  
  while (!stopflag) {
    
    # if(iterationSummary==TRUE) 
      message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
    
    if(lastNoError==TRUE) params <- omxGetParameters(model)
    
    if(lastBestFitCount == 0 & numdone > 0){ #if the last fit was not the best
      if(exists('bestfit')) params<-bestfit.params #if bestfit exists use this instead
      if(numdone %% 4 == 0) params<-inits #sometimes, use initial start values instead
      
      model <- omxSetParameters(model, labels = names(params), 
        values = params * rnorm(length(params),loc,scale) + rnorm(length(params),0,scale)
        )
      
      gradientStepSize <- initialGradientStepSize
                  tolerance <- initialTolerance
      gradientIterations<-initialGradientIterations
    }#end if last fit not best section
    
   
    if(lastBestFitCount > 0){ #if the last fit was the best so far
      if(exists('bestfit')) {
        params<-bestfit.params      
        model<-bestfit
      }
      
      if(defaultComputePlan==TRUE){
        
        if(lastBestFitCount == 2) gradientStepSize <- gradientStepSize *.1
        if(lastBestFitCount == 3) gradientStepSize <- gradientStepSize *10
        if(lastBestFitCount == 5) gradientStepSize <- gradientStepSize *.1
        if(lastBestFitCount  > 0) tolerance<-tolerance * .001 
        if(lastBestFitCount  > 0) gradientIterations<-gradientIterations+2
        # if(lastBestFitCount  %in% seq(4,100,4)) gradientIterations<-gradientIterations-1
        if(lastBestFitCount > 4) model <- omxSetParameters(model, labels = names(bestfit.params), 
          values = bestfit.params * rnorm(length(params),loc,scale/10) + 
            rnorm(length(params),0,scale / 10)
          )
      }
      
      if(defaultComputePlan==FALSE){
        model <- omxSetParameters(model, labels = names(bestfit.params), 
          values = bestfit.params * rnorm(length(bestfit.params),loc,scale/10) + 
            rnorm(length(params),0,scale / 10)
        )
      }
    }#end if last fit was best section
    
    
    
    
    if(defaultComputePlan==TRUE) model<-mxModel(model, mxComputeSequence(list(
      mxComputeGradientDescent(verbose=verbose, gradientStepSize = gradientStepSize, 
        nudgeZeroStarts=FALSE,   gradientIterations = gradientIterations, tolerance=tolerance, 
        maxMajorIter=3000),
      mxComputeNumericDeriv(), mxComputeStandardError(),  
      mxComputeReportDeriv())))
    
    if(showInits==TRUE) {
      message('Starting values:  ')
      message(paste0(names(omxGetParameters(model)),' : ', omxGetParameters(model),'\n'))
    }
    
    if(any(is.na(omxGetParameters(model)))) browser()
    fit <- suppressWarnings(try(mxRun(model, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
    
    numdone <- numdone + 1
    
    
    if( class(fit) == "try-error" || is.na(fit$output$minimum) || fit$output$status$status== -1) {
      lastBestFitCount <- 0
      lastNoError<-FALSE
      message('\n Fit attempt generated errors') 
    }
    
    if(class(fit) != "try-error" & !is.na(fit$output$minimum) & fit$output$status$status != -1) { #if fit was not an error
      
      if (fit$output$minimum >= lowestminsofar + generalTolerance) {
        lastBestFitCount <- 0
        lastNoError<-TRUE
        message(paste0('\n Fit attempt worse than current best:  ',fit$output$minimum ,' vs ', lowestminsofar )) 
      }
      
      if(fit$output$minimum >= lowestminsofar) lastBestFitCount<-0
      
      if (fit$output$minimum < lowestminsofar) { #if this is the best fit so far
        message(paste0('\n Lowest minimum so far:  ',fit$output$minimum) )
        lastBestFitCount<-lastBestFitCount+1 
        lowestminsofar <- fit$output$minimum
        lastNoError<-TRUE
        bestfit <- fit
        bestfit.params <- omxGetParameters(bestfit)
      }

      if (fit$output$minimum <= lowestminsofar + generalTolerance) { #if this is the best fit or equal best so far, check the following
        ###########stopflag checks
        stopflag<-TRUE
        
          if(fit$output$status[[1]] > greenOK) {
            stopflag<-FALSE
            message(paste0('\n OpenMx status code ', fit$output$status[[1]], ' greater than ', as.numeric(greenOK)))
            
          }
          
   
          if(fit$output$minimum > fit2beat) {
            message(paste0('\n Fit value of ', fit$output$minimum, ' greater than fit2beat of ', fit2beat))
            stopflag<-FALSE
          }
          
          if(fit$output$minimum > lowestminsofar + generalTolerance){
            message(paste0('\n Fit value of ', fit$output$minimum, ' greater than best so far of ', lowestminsofar))
            stopflag<-FALSE
          }
        
        if(checkHess==TRUE) {
          hessEigenval <- try(eigen(fit$output$calculatedHessian, symmetric = T, only.values = T)$values)
          
          if(class(hessEigenval)=='try-error') {
            message(paste0('\n Eigenvalues of hessian could not be calculated'))
            stopflag<-FALSE
          }
          
          if(class(hessEigenval)!='try-error' && any(hessEigenval < 0)) {
            message(paste0('\n Not all eigenvalues of hessian are greater than ', 0,': ', paste(hessEigenval,collapse=', ')))
            stopflag<-FALSE
          }
          
#           if(stopflag==TRUE && any(hessEigenval < 0)) {
# message('Some hessian eigenvalues were negative but within tolerance of ', generalTolerance, ', hessian adjusted to provide standard errors.')
#             newhess<- eigen(fit$output$calculatedHessian)$vectors %*% diag(abs(hessEigenval)) %*% solve(eigen(fit$output$calculatedHessian)$vectors)
#             fit@output$calculatedHessian[,]<-newhess
#           }
          
          if(stopflag ==TRUE) bestfit <- fit #message(paste(hessEigenval,collapse=', '))
        }
        
      }#end stopflag checks and if lowest min section
      
      if (!stopflag) {        
        if(iterationSummary==TRUE){
          message(paste0("\n Attempt ",numdone," fit:  "))
          message(paste(names(params),": ", fit$output$estimate,"\n"))
          message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
        }
      }
      
      if(stopflag){
        message('\nSolution found\n')
        if(length(bestfit$intervals)>0){ #only calculate confidence intervals once the best fit is established
          
          message("Estimating confidence intervals\n") 
          
          if(defaultComputePlan==TRUE) bestfit <- OpenMx::mxModel(bestfit, 
            mxComputeSequence(list(
              mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
                gradientIterations=gradientIterations, tolerance=tolerance, 
                maxMajorIter=3000),
                constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')),
              mxComputeNumericDeriv(), mxComputeStandardError(), 
              mxComputeReportDeriv())))
          
          cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
          
          
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors\n')
          } else {
            if (length(OpenMx::summary(cifit)$npsolMessage) > 0) message('Warning messages generated from confidence interval refit\n')
            bestfit<-cifit
          }
          
        }
        if (length(OpenMx::summary(bestfit)$npsolMessage) > 0) {
          warning(OpenMx::summary(bestfit)$npsolMessage)
        }
        
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
        
      }
    } #end 'if fit not an error' section
    
    
    
    
    
    
    if (numdone > extraTries & stopflag==FALSE) { #added stopflag==FALSE
      message('\nRetry limit reached')
      stopflag <- TRUE
      if (exists("bestfit")) {
        
        if(length(bestfit$intervals)>0){ #calculate intervals for best fit, even though imperfect
          message("Estimate confidence intervals for imperfect solution\n") 
          
          if(defaultComputePlan==TRUE) bestfit <- OpenMx::mxModel(bestfit, 
            mxComputeSequence(list(
              mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
                gradientIterations=gradientIterations, tolerance=tolerance, 
                maxMajorIter=3000),
                constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')),
              mxComputeNumericDeriv(), mxComputeStandardError(), 
              mxComputeReportDeriv())))
          
          cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
          
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors, returning fit without confidence intervals\n')
          } else {
            bestfit<-cifit
          }
        }
        if (length(bestfit$output$status$statusMsg) > 0) { 
          warning(bestfit$output$status$statusMsg)
        }
        if(bestfit$output$status$code==6) message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
      }
    }
  } #end while loop
  
  
  
  if(bestInitsOutput & exists("bestfit")){
    bestfit.params <- omxGetParameters(bestfit)
    message("\nStart values from best fit:")
    if(paste) message(paste(bestfit.params, sep=",", collapse = ",")) 
    if(!paste)  message(paste(names(bestfit.params),": ", bestfit.params,"\n"))
  }
  
  if (!exists("bestfit")) {
    if(class(fit) == 'try-error') warning(fit[[length(fit)]])
    message('All fit attempts resulted in errors - check starting values or model specification')
    bestfit<-fit
  }
  
  return(bestfit)
}