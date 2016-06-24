#' internal function for ctsem to calculate saturated fit

satFit  <- function(fittedmodel){  

    datawide<-fittedmodel@data@observed[,-grep("dT",
      colnames(fittedmodel@data@observed)),drop=FALSE]

  
  datawide<-data.frame(datawide)
  
  mvalues<-sapply(datawide,mean,na.rm = TRUE) #set mean values
  mvalues[is.na(mvalues)]<-0 #set nonexisting means of data to 0
  

  SatModel  <- OpenMx::mxModel("SatModel",
    
    mxMatrix(type = "Full",nrow = 1,ncol = ncol(datawide),free = TRUE,
      values = mvalues,name = "M",
      dimnames=list("row",colnames(datawide))),
    
    mxMatrix(type = "Full",nrow = ncol(datawide),ncol = ncol(datawide),
      free = FALSE,values = 0,name = "A"),
    
    mxMatrix(type = "Full",nrow = ncol(datawide),ncol = ncol(datawide),
      free = FALSE,values = diag(1,ncol(datawide)),name = "F",
      dimnames=list(colnames(datawide),colnames(datawide))),
    
    mxData(observed = datawide,type = "raw"),
    mxMatrix(type = "Iden",nrow = ncol(datawide),ncol = ncol(datawide),name = "bigI"), #identity
    
    mxExpectationRAM(M="M"),
    type="RAM"
  )
  
  if(nrow(datawide)==1){ #use custom r fiml row likelihood function
    SatModel <- OpenMx::mxModel(SatModel,
      
      mxMatrix(type = "Symm",nrow = ncol(datawide), ncol = ncol(datawide),
        free = FALSE,  #if only single subject then no covariance possible
        values = 0,name = "S"),
    
    mxFitFunctionR(CondLogLik),
    mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F),name = "expCov"),
    mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)),name = "expMean")
    )
  }

if(nrow(datawide)>1){ #use openmx fiml likelihood
   
  sdiag<-diag(ncol(datawide))
  svalues<-var(datawide,use="pairwise.complete.obs")
  svalues[is.na(svalues)]<-sdiag[is.na(svalues)]
  SatModel <- OpenMx::mxModel(SatModel,
    
    mxMatrix(name = "S", type = "Full",nrow = ncol(datawide),ncol = ncol(datawide),
      free = TRUE,  #if multiple subjects then covariance possible
      values = svalues,
      labels=indexMatrix(dimension=nrow(var(datawide,use="pairwise.complete.obs")),
        symmetrical=TRUE,sep="_",starttext="cov"))
#     
#     mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F),name = "expCov"),
#     mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)),name = "expMean"),
#     mxExpectationNormal(covariance = "expCov",means = "expMean",dimnames = colnames(datawide)),
#     mxFitFunctionML()
  )
}
satfit <- OpenMx::mxRun(SatModel,useOptimizer=ifelse(nrow(datawide)>1,TRUE,FALSE)) #if only single subject then we don't optimize, just take means model
return(satfit)
}
