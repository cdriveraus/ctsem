#' ctPSMfit
#' 
#' Fits a specified ctsem model (without predictors) using the PSM package. 
#'
#' @param datawide a dataset in the wide format used by ctsem.
#' @param ctmodelobj A ctsem model specified using \code{\link{ctModel}}. As predictors (covariates) are 
#' not implemented currently, ensure the model contains none of these or the fit may fail.
#' @param omxStartValues Named vector of start values, equivalent to this argument for \code{\link{ctFit}}.
#' @param ... Additional parameters to pass to PSM.#'
#' @return PSM fit data
#' @examples
#' ## Examples set to 'dontrun' because they take longer than 5s. 
#' 
#' \dontrun{
#' generatingModel <- ctModel(n.latent = 1, n.manifest = 1, Tpoints = 10,
#' LAMBDA = diag(1), DRIFT = matrix(-.3, nrow = 1),
#' MANIFESTVAR = diag(1),
#' CINT = matrix(3, 1, 1),
#' DIFFUSION = t(chol(diag(5, 1))))
#' 
#' dat <- ctGenerate(generatingModel, n.subjects=10, burnin=300)
#' 
#' ### ctsem model and fit
#' ctsemModel <- ctModel(n.latent=1, n.manifest = 1, Tpoints = 10,
#'   LAMBDA = diag(1))
#' ctsemFit <- ctFit(dat, ctsemModel, stationary = c('T0VAR'))
#' 
#' ### fit with PSM
#' psmFit <- ctPSMfit(dat, omxStartValues = 
#'     omxGetParameters(ctsemFit$mxobj), ctsemModel)
#' }
#' @export

ctPSMfit<-function(datawide,ctmodelobj,omxStartValues=NULL, ...){
  
  if (!requireNamespace("PSM", quietly = TRUE)) {
    stop("PSM package needed for this function to work. Please install it.",
      call. = FALSE)
  }
  
  psmdata<-lapply(1:nrow(datawide), function(x) {
    
    long<-ctWideToLong(datawide[x,,drop=F], Tpoints=ctmodelobj$Tpoints, 
      n.manifest=ctmodelobj$n.manifest, manifestNames=ctmodelobj$manifestNames)
    
    long<-suppressMessages(ctDeintervalise(long))
    colnames(long)[which(colnames(long)=='AbsTime')] <-'Time'
    out<-list(Time=long[,'Time'], Y=t(long[,ctmodelobj$manifestNames]))
    out$U=matrix(1,nrow=1,ncol=ncol(out$Y))
    return(out)
  })
  
  #fixed value transformations
  if(any(!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$MANIFESTVAR)))))){
    diag(ctmodelobj$MANIFESTVAR)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$MANIFESTVAR))))] <- 
      log(as.numeric(diag(ctmodelobj$MANIFESTVAR)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$MANIFESTVAR))))]))
    diag(ctmodelobj$MANIFESTVAR)[diag(ctmodelobj$MANIFESTVAR)== -Inf][!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$MANIFESTVAR))))] <- -999
  }
  
  if(any(!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DRIFT)))))){
    if(any(diag(ctmodelobj$DRIFT)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DRIFT))))] >=0)) message(
      'non negative DRIFT diagonal specified, setting to -.00001.')
    diag(ctmodelobj$DRIFT)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DRIFT))))] <- 
      suppressWarnings(log(-as.numeric(diag(ctmodelobj$DRIFT)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DRIFT))))])) + .00001)
    
    diag(ctmodelobj$DRIFT)[is.nan(diag(ctmodelobj$DRIFT)) | 
        diag(ctmodelobj$DRIFT) == -Inf ] <- -999
  }
  
  if(any(!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DIFFUSION)))))){
    if(any(diag(ctmodelobj$DIFFUSION)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DIFFUSION))))] <=0)) message(
      'non positive DIFFUSION diagonal specified, setting to .00001.')
    diag(ctmodelobj$DIFFUSION)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DIFFUSION))))] <- 
      suppressWarnings(log(as.numeric(diag(ctmodelobj$DIFFUSION)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$DIFFUSION))))]) - .00001) )
    diag(ctmodelobj$DIFFUSION)[is.nan(diag(ctmodelobj$DIFFUSION)) | diag(ctmodelobj$DIFFUSION)== -Inf] <- -999
  }
  
  if(any(!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$TRAITVAR)))))){
    if(any(diag(ctmodelobj$TRAITVAR)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$TRAITVAR))))] <=0)) message(
      'non positive TRAITVAR diagonal specified, setting to .00001.')
    diag(ctmodelobj$TRAITVAR)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$TRAITVAR))))] <- 
      suppressWarnings(log(as.numeric(diag(ctmodelobj$TRAITVAR)[!is.na(suppressWarnings(as.numeric(diag(ctmodelobj$TRAITVAR))))]) - .00001) )
    diag(ctmodelobj$TRAITVAR)[is.nan(diag(ctmodelobj$TRAITVAR)) | diag(ctmodelobj$TRAITVAR)== -Inf] <- -999
  }
  
  
  
  flatcm<-unlist(ctmodelobj[-which(names(ctmodelobj) %in% c('latentNames', 'manifestNames', 'T0VAR'))])
  flatcm<-flatcm[-suppressWarnings(which(!is.na(as.numeric(flatcm))))]
  names(flatcm)<-flatcm
  inits<-unlist(lapply(flatcm, function(x) stats::rnorm(1,.3,.01)))
  if(!is.null(omxStartValues)) {
    nameMatch<-names(inits)[which(names(inits) %in% names(omxStartValues))] 
    inits<- omxStartValues[nameMatch]
  }
  
  MyModel <- vector(mode="list")
  MyModel$Matrices=function(phi) {
    
    matA<-ctmodelobj$DRIFT
    
    matA[which(matA %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matA[which(matA %in% names(phi))],names(phi))]))
    mode(matA)<-'numeric'
    diag(matA)<- -exp(diag(matA))
    
    matB<-ctmodelobj$CINT
    matB[which(matB %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matB[which(matB%in% names(phi))],names(phi))]))
    mode(matB)<-'numeric'
    
    matC<-ctmodelobj$LAMBDA
    matC[which(matC %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matC[which(matC %in% names(phi))],names(phi))]))
    mode(matC)<-'numeric'
    
    matD<-ctmodelobj$MANIFESTMEANS
    matD[which(matD %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matD[which(matD%in% names(phi))],names(phi))]))
    mode(matD)<-'numeric'
    
    out<-list(matA=matA,
      matB=matB,
      matC=matC,
      matD=matD
    )
    
    out
  }
  
  MyModel$h = function(eta,theta,covar) {
    
    
    if(!is.null(ctmodelobj$TRAITVAR)){
      
      phi <- theta
      
      matA<-ctmodelobj$DRIFT
      matA[which(matA %in% names(phi))] <- 
        as.numeric(unlist(phi[match(matA[which(matA %in% names(phi))],names(phi))]))
      mode(matA)<-'numeric'
      diag(matA)<- -exp(diag(matA))
      
      
      # message(paste0('eta = ', eta))
      names(theta)<-flatcm[-which(flatcm %in% ctmodelobj$TRAITVAR )]
      
      matB<-ctmodelobj$CINT
      matB[which(matB %in% names(theta))] <- 
        as.numeric(unlist(theta[match(matB[which(matB %in% names(theta))],names(theta))]))
      mode(matB)<-'numeric'
      
      matB<- matB -matA  %*% eta
      
      
      phi[which(names(phi) %in% ctmodelobj$CINT)] <- matB #this assignment may not be reliable - relies on ordering of phi matching matB
      
      
      X0<-ctmodelobj$T0MEANS
      X0[which(X0 %in% names(theta))] <- 
        as.numeric(unlist(theta[match(X0[which(X0 %in% names(theta))],names(theta))]))
      mode(X0)<-'numeric'
      
      X0<- X0 + eta
      
      phi[which(names(phi) %in% ctmodelobj$T0MEANS)] <- X0 
      
      return(phi)
    }
    
    
    if(is.null(ctmodelobj$TRAITVAR)){
      phi <- theta
      return(phi)
    }
    
  }
  
  MyModel$S = function(phi) {
    S<-ctmodelobj$MANIFESTVAR
    
    S[which(S %in% names(phi))] <- 
      as.numeric(unlist(phi[match(S[which(S %in% names(phi))],names(phi))]))
    mode(S)<-'numeric'
    diag(S)<-exp(diag(S))
    S<-t(S) %*% S
    S
    
  }
  MyModel$SIG = function(phi) {
    SIG<-ctmodelobj$DIFFUSION
    
    SIG[which(SIG %in% names(phi))] <- 
      as.numeric(unlist(phi[match(SIG[which(SIG %in% names(phi))],names(phi))]))
    mode(SIG)<-'numeric'
    
    diag(SIG)<-exp(diag(SIG))
    
    # SIG<-t(SIG) %*% SIG
    # SIG[,]<-0.000
    SIG
  }
  MyModel$X0 = function(Time,phi,U) {
    X0<-ctmodelobj$T0MEANS
    
    X0[which(X0 %in% names(phi))] <- 
      as.numeric(unlist(phi[match(X0[which(X0 %in% names(phi))],names(phi))]))
    mode(X0)<-'numeric'
    X0
  }
  MyModel$ModelPar = function(THETA) {
    
    names(THETA)<-flatcm
    
    if( !is.null(ctmodelobj$TRAITVAR)){
      
      OMEGA<-ctmodelobj$TRAITVAR
      out<-list(theta=lapply(THETA[-which(names(THETA) %in% OMEGA )], function(x) x))
      names(out$theta)<-flatcm[-which(flatcm %in% OMEGA )]
      
      OMEGA[which(OMEGA %in% names(THETA))] <- 
        as.numeric(unlist(THETA[match(OMEGA[which(OMEGA %in% names(THETA))],names(THETA))]))
      mode(OMEGA)<-'numeric'
      diag(OMEGA)<-exp(diag(OMEGA))
      OMEGA<-t(OMEGA) %*% OMEGA
      
      #     matA<-ctmodelobj$DRIFT
      #     matA[which(matA %in% names(THETA))] <- 
      #       as.numeric(unlist(THETA[match(matA[which(matA %in% names(THETA))],names(THETA))]))
      #     mode(matA)<-'numeric'
      #     diag(matA)<- -exp(diag(matA))
      #     
      #     OMEGA <- matA %*% OMEGA %*% t(matA)
      out$OMEGA<-OMEGA
      
    }
    
    
    if( is.null(ctmodelobj$TRAITVAR)){
      out<-list(theta=lapply(THETA, function(x) x))
      names(out$theta)<-flatcm
    }
    
    out
  }
  
  
  fit<-PSM::PSM.estimate(MyModel, psmdata, Par=list(Init=inits), 
    control=list(optimizer='ucminf', maxeval=3000),fast=F,...)
  
  return(list(PSMdata=psmdata, PSMfit=fit, PSMmodel=MyModel, PSMinitpars=list(Init=inits)))
}
