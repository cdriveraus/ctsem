#'kalmanPredict
#'@export
kalmanPredict<-function(mxobj,submodel=NULL,mxstyle=TRUE,timeVaryingParams=FALSE,...){
 
  KalmanFilter <- function(A, B, C, D, Q, R, x, y, u, P){ #Mike Hunter's Kalman Filter from OpenMx, copied here for convenience.
browser()
    x <- A %*% x + B %*% u
    P <- A %*% P %*% t(A) + Q
    # if(any(diag(P) < 0)) P[diag(P)<0]<-.0001
    x.pred <- x
    P.pred <- P
    r <- y - (C %*% x + D %*% u)
    notMiss <- !is.na(r)
    r[!notMiss] <- 0
    if(length(r)==sum(!notMiss)){#all missing row
      m2ll <- log(det(C %*% P %*% t(C) + R))
      return(list(x.pred=x.pred, P.pred=P.pred, x.upda=x.pred, P.upda=P.pred, m2ll=m2ll, L=exp(m2ll/-2) ))
    } else {
      Cf <- C[notMiss,,drop=F]
      Rf <- R[notMiss, notMiss,drop=F]
      S <- C %*% P %*% t(C) + R # was Cf %*% P %*% t(Cf) + Rf 
      Sf <- S[notMiss, notMiss,drop=F] #added S filtered - is this correct?

      Sinv <- try(solve(S) )
      if(class(Sinv) == 'try-error') Sinv <- solve(chol(S)) %*% t(solve(chol(S)))
      Sfinv <-try(solve(Sf),silent=T) #added S filtered - correct?
      if(class(Sfinv) == 'try-error') Sfinv <- solve(chol(Sf)) %*% t(solve(chol(Sf)))
      rf <- matrix(r[notMiss], ncol=1)
      K <- P %*% t(C) %*% Sinv
      Kf <- K[,notMiss,drop=F] #added K filtered - correct?
      xu <- x + Kf %*% rf #included k filtered 
      Pu <- P - K %*% C %*% P #included K filtered
      
      const <- length(rf)*log(2*pi)
      m2ll <- log(det(Sf)) + t(rf) %*% Sfinv %*% rf + const
      print(m2ll)
      return(list(x.pred=x.pred, P.pred=P.pred, x.upda=xu, P.upda=Pu, m2ll=m2ll, L=exp(m2ll/-2) ))
    }
  }


    data<-mxobj$data$observed
  manifestNames<-rownames(mxobj[[mxobj$expectation$C]]$labels)
  latentNames<-colnames(mxobj[[mxobj$expectation$C]]$labels)
  expectation <- mxobj$expectation

  
manifestForecast<-matrix(NA,ncol=length(manifestNames),nrow=nrow(data)+1)
manifestUpdated<-matrix(NA,ncol=length(manifestNames),nrow=nrow(data)+1)
latentForecast<-matrix(NA,ncol=length(latentNames),nrow=nrow(data)+1)
latentUpdated<-matrix(NA,ncol=length(latentNames),nrow=nrow(data)+1)

  C<-mxEvalByName(expectation$C, mxobj,compute=T,defvar.row=2)
x0<-mxEvalByName(expectation$x0, mxobj, compute=T,defvar.row=2)
P0<-mxEvalByName(expectation$P0, mxobj, compute=TRUE,defvar.row=2)
PPredicted<-list()
PUpdated<-list()
PPredicted[[1]]<-P0
PUpdated[[1]]<-P0
m2ll<-numeric(nrow(data)+1)
m2ll[1]<-0

  manifestUpdated[1,] <- data[1,manifestNames]
  manifestUpdated[1,][is.na(manifestUpdated[1,])] <- (C %*% x0)[is.na(manifestUpdated[1,])]
  manifestForecast[1,] <- C %*% x0
  
latentUpdated[1,]<- x0
  

  A<-mxEvalByName(expectation$A,mxobj,compute=T,defvar.row=2)
  B<-mxEvalByName(expectation$B,mxobj,compute=T,defvar.row=2)
  C<-mxEvalByName(expectation$C,mxobj,compute=T,defvar.row=2)
  D<-mxEvalByName(expectation$D,mxobj,compute=T,defvar.row=2)
  Q<-mxEvalByName(expectation$Q,mxobj,compute=T,defvar.row=2)
  R<-mxEvalByName(expectation$R,mxobj,compute=T,defvar.row=2)
  u<-mxEvalByName(expectation$u,mxobj,compute=T,defvar.row=2)

  for(i in 2:(nrow(data)+1)){

    if(timeVaryingParams==TRUE){
      A<-mxEvalByName(expectation$A,mxobj,compute=T,defvar.row=i)
    B<-mxEvalByName(expectation$B,mxobj,compute=T,defvar.row=i)
    C<-mxEvalByName(expectation$C,mxobj,compute=T,defvar.row=i)
    D<-mxEvalByName(expectation$D,mxobj,compute=T,defvar.row=i)
    Q<-mxEvalByName(expectation$Q,mxobj,compute=T,defvar.row=i)
    R<-mxEvalByName(expectation$R,mxobj,compute=T,defvar.row=i)
    u<-mxEvalByName(expectation$u,mxobj,compute=T,defvar.row=i)
    P0<-mxEvalByName(expectation@P0, mxobj, compute=TRUE,defvar.row=i)
    }

    if(mxstyle==TRUE){
      
      kalman<-try(KalmanFilter(A=A, B=B, C=C, D=D, Q=Q, R=R, x=matrix(latentUpdated[i-1,]), 
      y=matrix(unlist(data[i-1,rownames(C)])), u=u, P=PUpdated[[i-1]]))
      if(class(kalman)=='try-error') browser()
      
      latentForecast[i,]<-kalman$x.pred
      m2ll[i]<-kalman$m2ll
      latentForecast[i,]<- A %*% latentUpdated[i-1,] + B %*% u
      
      latentUpdated[i,] <- kalman$x.upda
      PPredicted[[i]] <- kalman$P.pred
      PUpdated[[i]] <- kalman$P.upda
      
      manifestForecast[i,]<-C %*% latentForecast[i,] + D %*% u
      
      manifestUpdated[i,] <- C %*% latentForecast[i,] + D %*% u
      # manifestUpdated[i,][is.na(manifestUpdated[i,])] <- manifestForecast[i,][is.na(manifestUpdated[i,])]
    }
    
    if(mxstyle==FALSE){
    latentForecast[i,]<- A %*% latentUpdated[i-1,] + B %*% u
   
    manifestForecast[i,]<-C %*% latentForecast[i,] + D %*% u
    
    manifestUpdated[i,] <- data[i,manifestNames]
    manifestUpdated[i,][is.na(manifestUpdated[i,])] <- manifestForecast[i,][is.na(manifestUpdated[i,])]
    
    latentUpdated[i,] <- solve(C) %*% manifestUpdated[i,]
    }
    

    
    
    
    

  }
  
  return(list(data=data,manifestUpdated=manifestUpdated,manifestForecast=manifestForecast,
    latentUpdated=latentUpdated,latentForecast=latentForecast, m2ll=m2ll))
}


# predictor<-models[[1]]
# predictor<-mxModel(predictor$mxobj,
#   # mxAlgebra(name='latentForecast', DRIFT %*% data.Y1 + intd1 %*% u),
#   mxAlgebra(name='latentForecast',  data.Y1),
#   # mxAlgebra(name='rowAlgebra', LAMBDA %*% latentForecast + D %*% u),
#   mxAlgebra(name='rowAlgebra', data.Y1 ),
#   mxAlgebra(name='reduceAlgebra', sum(rowResults)),
#   mxFitFunctionRow('rowAlgebra', 'reduceAlgebra', 'Y1'),
#   mxComputeOnce('fitfunction', 'fit')
# )
# predictor$expectation<-NULL
# predictions<-mxRun(predictor)
# 
# latentForecast[i,]<- A %*% latentUpdated[i-1,] + B %*% u
# 
# manifestForecast[i,]<-C %*% latentForecast[i,] + D %*% u
    
    