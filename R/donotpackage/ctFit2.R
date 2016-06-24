# fd<-fdHess(start.values, modell.2)
# g <- 1/fd$gradient
# out<-optim(start.values, modell.2, method="SANN", hessian=TRUE, control=list(trace=2, parscale=g)) 


ctOptim<-function(data,ctmodelobjvalues,ctmodelobjlabels,...){

  pars<-c()
  for(i in 1:length(ctmodelobjlabels)){
    if(is.matrix(ctmodelobjlabels[[i]])){
#       
newparnames<-unique(ctmodelobjlabels[[i]][is.na(suppressWarnings(as.numeric(ctmodelobjlabels[[i]])))])
      newpars<-ctmodelobjvalues[[i]][match(newparnames, ctmodelobjlabels[[i]])]
      names(newpars)<-newparnames
      pars<-c(pars,newpars)
}
  }
  message(paste0(names(pars),': ',pars,'\n'))
  
#   
# browser()
#   grad<-grad(func=ctOptimFit,x=pars,method='simple',data=data,ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels)
#   
#   gradf<-  function(x){#,data=data,ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels){
#     grad<-grad(func=ctOptimFit,x=pars,method='simple',data=data,ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels)
#     return(grad)
#   }
  
  fit<-optim(par=pars,fn=ctOptimFit,data=data,ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels)
#     control=list(parscale=1/grad),...)
#   browser()
#   fit<-nloptr(x0=pars,eval_f=ctOptimFit,eval_grad_f=gradf,data=data,ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels,...)
#   fit<-nloptr(x0=pars,eval_f=ctOptimFit,data=data,
#     ctmodelobjvalues=ctmodelobjvalues,ctmodelobjlabels=ctmodelobjlabels,
#     opts=list('algorithm'='NLOPT_LN_NELDERMEAD'))
}
  
  
  ctOptimFit<-function(ctOptimFitpars,data,ctmodelobjvalues,ctmodelobjlabels){
#     
    for(i in 1:length(ctmodelobjlabels)){
      if(is.matrix(ctmodelobjlabels[[i]]) & any(is.na(suppressWarnings(as.numeric(ctmodelobjlabels[[i]]))))){
#         
         ctmodelobjvalues[[i]][is.na(suppressWarnings(as.numeric(ctmodelobjlabels[[i]])))] <-
#           ctOptimFitpars[names(ctOptimFitpars) %in% ctmodelobjlabels[[i]][is.na(suppressWarnings(as.numeric(ctmodelobjlabels[[i]])))]]
          ctOptimFitpars[ctmodelobjlabels[[i]][is.na(suppressWarnings(as.numeric(ctmodelobjlabels[[i]])))]]
      }
    }
    fit<-try(ctFit2(data=data,ctmodelobj=ctmodelobjvalues))
#     
    return(fit)
  }
    
    
    

ctFit2<-function(data,ctmodelobj){
  for(i in 1:length(ctmodelobj)){ #this loop reads in the specified continuous time model so the objects are available
    assign(names(ctmodelobj[i]),ctmodelobj[[i]])
    
    if(is.matrix(ctmodelobj[[i]])) { #if it's a matrix then ensure that it is numeric
      suppressWarnings(assign(names(ctmodelobj[i]), matrix(as.numeric(ctmodelobj[[i]]),nrow(ctmodelobj[[i]]))))
    }
  }
  
  manifests<-data[,1:(Tpoints*n.manifest)]
  intervals<-data[,(Tpoints*n.manifest + (Tpoints-1)*n.TDpred + 1) : (Tpoints*n.manifest + (Tpoints-1)*n.TDpred + Tpoints-1)]
  
  S<-matrix(0,nrow = (n.latent + n.manifest) * Tpoints+n.latent, ncol = (n.latent + n.manifest) * Tpoints+n.latent)
  A<-matrix(0,nrow = (n.latent + n.manifest) * Tpoints+n.latent, ncol = (n.latent + n.manifest) * Tpoints+n.latent)
  M<-matrix(0,nrow = 1, ncol = (n.latent + n.manifest) * Tpoints+n.latent)
  filter<-cbind(matrix(0,nrow=n.manifest*Tpoints,ncol=n.latent*Tpoints+n.latent), diag(n.manifest*Tpoints))
  bigI<-diag(nrow(A))
  II<-diag(n.latent)
  DRIFTHATCH <- DRIFT %x% II + II %x% DRIFT 
  log2pi <- log(2*pi)
  ll<-0
 
  
  for(i in 1:nrow(data)){
    
    if(i==1 | (i>1 & identical(intervals[i,], intervals[i-1,]))){ #if intervals are different to last subject then recalculate discrete matrices
        discreteDRIFT<- matrix(apply( intervals[i,,drop=F],2,function(x) expm(DRIFT %x% x)),nrow=n.latent)
        
        
        discreteDIFFUSION<- matrix(apply( intervals[i,,drop=F],2,function(x) {
          solve(DRIFTHATCH) %*% ((expm(DRIFTHATCH %x% x)) - (II %x% II)) %*% rvectorize(DIFFUSION) 
        }
        ),nrow=n.latent)
        
        discreteCINT<- matrix(apply( intervals[i,,drop=F],2,function(x) {
          solve(DRIFT) %*% (expm(DRIFT) - II) %*% CINT  
        }
        ),nrow=1)
        
        discreteTRAITtemp<- matrix(apply( intervals[i,,drop=F],2,function(x) {
          solve(DRIFT) %*% (expm(DRIFT %x% x) - II)
        }
        ),nrow=n.latent)
      
      discreteTRAIT<-t(discreteTRAITtemp)
      for(j in 1:(Tpoints-1)){
        discreteTRAIT[((j-1)*n.latent+1) : (j*n.latent), ] <-discreteTRAITtemp[,((j-1)*n.latent+1) : (j*n.latent)]
      }
    
      
      
    discreteDRIFT-> A[(row(A) - 1 - n.latent)%/%n.latent ==  # when the rows and columns,grouped by n.latent,
        (col(A) - 1)%/%n.latent & row(A) <= n.latent*Tpoints]#are equal,and not greater than the total number of latent variables
    
    LAMBDA -> A[(row(A) - 1 - (Tpoints*n.latent+n.latent))%/%n.manifest == #insert LAMBDA loadings into Avalues
        (col(A) - 1)%/%n.latent & 
        col(A) <= n.latent*Tpoints+n.latent]
    
    MANIFESTVAR -> S[(row(S) - n.latent*Tpoints - 1)%/%n.manifest ==  
        (col(S) - n.latent*Tpoints - 1)%/%n.manifest &
        col(S) >  n.latent*Tpoints &
        row(S) >  n.latent*Tpoints]
    
    discreteDIFFUSION -> S[(row(S) - 1)%/%n.latent == (col(S) - 1)%/%n.latent & #initial DIFFUSION values with fixed coding      
        col(S) <= n.latent*Tpoints &
        col(S) >n.latent]
    
    T0VAR -> S[1:n.latent,1:n.latent]
      
    T0TRAITEFFECT -> A[1:n.latent,(n.latent*Tpoints+1):(n.latent*Tpoints+n.latent)]
      
      discreteTRAIT -> A[(n.latent+1) : (n.latent*Tpoints), (n.latent*Tpoints+1):(n.latent*Tpoints+n.latent)]
      
      TRAITVAR -> S[(n.latent * Tpoints + 1) : (n.latent*Tpoints+n.latent),  (n.latent * Tpoints + 1) : (n.latent*Tpoints+n.latent)]
    
    T0MEANS -> M[1:n.latent]
    
    discreteCINT -> M[(n.latent+1):(n.latent*Tpoints)]
    
    MANIFESTMEANS -> M[(n.latent*Tpoints+n.latent+1) : (n.latent*Tpoints++n.latent+n.manifest*Tpoints)]
    }
    
    browser()
    
    existenceVec<-!is.na(manifests[i,])
    if(i==1 | (i > 1 && (existenceVec != oldExistenceVec))) { #if this subject has different missingness to previous
    expCov<- filter %*% solve(bigI - A) %*% S%*% t( solve(bigI - A)) %*% t(filter)
      
#       F * (I-A)-1 * S * ((I-A) -1)' * F'
      
    expMean<- t(filter%*%(solve(bigI - A)) %*%t(M))
      
      expCovi<-expCov[existenceVec,existenceVec,drop=FALSE]
      expMeani<-expMean[,existenceVec,drop=FALSE]
    }
      


    
    expCovCheck<-det(expCovi) > 0 & all(eigen(expCovi)$values > 0)

    if(expCovCheck){
    rowll <- length(existenceVec[existenceVec==TRUE]) * log2pi +
      log(det(expCov)) + 
      (manifests[i, ] - expMeani) %*% ((solve(chol(expCovi))) %*% t(solve(chol(expCovi)))) %*% t(manifests[i, ] - expMeani)
    }
    if(expCovCheck==FALSE) {
      ll<-9999999
      message(ll)
      return(ll)
    }
    
    oldExistenceVec<-existenceVec
    ll<-ll+rowll
  }
  message(ll)
 
  return(ll)
}


