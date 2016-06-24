#' internal ctsem function to fit an independence model to data

indFit  <- function(fittedmodel){
  if(class(fittedmodel)!="ctsemFit"){
    datawide<-fittedmodel@data@observed[,-grep("dT",
      colnames(fittedmodel@data@observed)),drop=FALSE]
  }
  
  if(class(fittedmodel)=="ctsemFit"){ 
    datawide<-fittedmodel$fit@data@observed[,-grep("dT",
      colnames(fittedmodel$fit@data@observed)),drop=FALSE]
  }  
  
  datawide<-data.frame(datawide)
  if(nrow(datawide)==1) { #if only 1 subject, use saturated model function as it is equivalent
    indfit <- satFit(fittedmodel)
    return(indfit)
  }
  
  if(nrow(datawide)>1) { #otherwise continue with independence model
  
  indModel  <- OpenMx::mxModel("indModel",
    
    mxMatrix(type = "Full",nrow = 1,ncol = ncol(datawide),free = TRUE,
      values = sapply(datawide,mean,na.rm = TRUE),name = "M",
      dimnames=list("row",colnames(datawide))),
    
    mxMatrix(type = "Symm",nrow = ncol(datawide),ncol = ncol(datawide),
      free = diag(ncol(datawide))==1,values = diag(ncol(datawide)),name = "S"),
    
    mxMatrix(type = "Full",nrow = ncol(datawide),ncol = ncol(datawide),
      free = FALSE,values = 0,name = "A"),
    
    mxMatrix(type = "Full",nrow = ncol(datawide),ncol = ncol(datawide),
      free = FALSE,values = diag(1,ncol(datawide)),name = "F",
      
      dimnames=list(colnames(datawide),colnames(datawide))),
    
    mxData(observed = datawide,type = "raw"),
    type="RAM",
    mxExpectationRAM(M = "M"),
    mxFitFunctionML()
  )
  indfit <- OpenMx::mxRun(indModel)
  summaryindfit<-summary(indfit)
  llInd <- summaryindfit$Minus2LogLikelihood
  dfInd <- summaryindfit$degreesOfFreedom
  indfit@output$IndependenceLikelihood <- llInd
  indfit@output$IndependenceDoF <- dfInd
  
  return(indfit)
}
#   summary_SatFit <- summary(SatModelFit)
#   return(c(summary_SatFit$Minus2LogLikelihood,summary_SatFit$degreesOfFreedom))
}