#' Fit ctsem models with Kalman filter R code -- slower!
#'
#' @param datalong long data
#' @param ctmodel ctsem model of type 'omx'
#' @param ... arguments to pass to ctsem::Kalman
#'
#' @return matrix of estimates and standard errors
#' @export
#'
#' @examples
#' \dontrun{
#' Tpoints<-250
#' n.manifest=5
#' gm<-ctModel(type='omx',n.latent=1,n.manifest=n.manifest,Tpoints=Tpoints,LAMBDA=matrix(rep(1,n.manifest),ncol=1),
#'   DRIFT=diag(-.3,1),
#'   T0VAR=diag(1),
#'   MANIFESTMEANS=matrix(0,nrow=n.manifest),
#'   MANIFESTVAR=t(chol(diag(.001,n.manifest))),
#'   DIFFUSION=t(chol(diag(1,1))))
#' 
#' cd<-ctGenerate(gm,n.subjects=1,burnin=300,wide=FALSE)
#' 
#' gm2=gm
#' gm2$DRIFT[1]='drift11'
#' 
#' cfit=ctFit(dat=cd,ctmodelobj=gm2,dataform='long')
#' cfit$mxobj$DRIFT$values
#' 
#' rfit=ctFitR(datalong=cd,ctmodel=gm2)
#' rfit
#' }
ctFitR<-function(datalong, ctmodel,binary=FALSE,carefulfit=TRUE,regulariser=1){
  
  # ctmodel <- list(ctmodel)
  kparsskeleton <- list(DRIFT=ctmodel$DRIFT,DIFFUSION=ctmodel$DIFFUSION,MANIFESTVAR=ctmodel$MANIFESTVAR,CINT=ctmodel$CINT,
    T0VAR=ctmodel$T0VAR,T0MEANS=ctmodel$T0MEANS,MANIFESTMEANS=ctmodel$MANIFESTMEANS,LAMBDA=ctmodel$LAMBDA)
  if(ctmodel$n.TDpred > 0) kparsskeleton$TDPREDS=ctmodel$TDPREDS
  kpars <- unlist(as.relistable(kparsskeleton))
  
  optimfunc <- function(pars){
    kpars[suppressWarnings(is.na(as.numeric(kpars)))] <- pars
    kpars <- relist(kpars)
    kpars <- lapply(kpars,function(x) matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x)))
    kpars$DIFFUSION <- kpars$DIFFUSION %*% t(kpars$DIFFUSION)
    kpars$T0VAR <- kpars$T0VAR %*% t(kpars$T0VAR)
    kpars$MANIFESTVAR <- kpars$MANIFESTVAR %*% t(kpars$MANIFESTVAR)
    kpars$DRIFT[row(kpars$DRIFT) == col(kpars$DRIFT)] <- -exp(diag(kpars$DRIFT))
  lp=-try(Kalman(kpars=kpars,datalong=datalong,manifestNames=ctmodel$manifestNames,
    latentNames=ctmodel$latentNames,TDpredNames=ctmodel$TDpredNames,timecol='time',optimize=TRUE,binary=binary))
  if(carefulfit) lp <- lp +abs(-1 -(sum(diag(kpars$DRIFT))))*regulariser
  }
  
  fit <- optim(par= c(rnorm(sum(suppressWarnings(is.na(as.numeric(kpars)))))),fn=optimfunc,hessian=TRUE,
    method = c("L-BFGS-B"),control=list(trace=3,REPORT=1,
      ndeps=rep(1e-6,sum(suppressWarnings(is.na(as.numeric(kpars)))))))
  
  parnames <- kpars[suppressWarnings(is.na(as.numeric(kpars)))]
  
  kpars[suppressWarnings(is.na(as.numeric(kpars)))] <- fit$par
  kpars <- relist(kpars)
  kpars <- lapply(kpars,function(x) matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x)))
  kpars$DIFFUSION <- kpars$DIFFUSION %*% t(kpars$DIFFUSION)
  kpars$T0VAR <- kpars$T0VAR %*% t(kpars$T0VAR)
  kpars$MANIFESTVAR <- kpars$MANIFESTVAR %*% t(kpars$MANIFESTVAR)
  kpars$DRIFT[row(kpars$DRIFT) == col(kpars$DRIFT)] <- -exp(diag(kpars$DRIFT))
  
  est <- unlist(kpars)[suppressWarnings(is.na(as.numeric(unlist(kparsskeleton))))]
  stderrors <- sqrt(diag(-fit$hessian))
  
  out <- matrix(c(est,stderrors),ncol=2)
  rownames(out) <- parnames
  colnames(out) <- c('Estimate', 'SE')
  return(out)
}
