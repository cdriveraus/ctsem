#' Generates data according to the model estimated in a ctsemFit object.#' 
#'
#' @param fit object of class ctsemFit as returned from \code{\link{ctFit}}
#' @param n.subjects integer. Number of subjects worth of data to generate
#' @param predictorSubjects vector of integers, or string 'all', defining which 
#' subjects to sample time dependent and independent predictors from.
#' @param ... parameters to pass to ctGenerate function, such as wide=FALSE.
#'
#' @return matrix of generated data
#' @export
#'
#' @examples
#' 
#' data(AnomAuth) 
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#'   Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#' 
#' dwide <- ctGenerateFromFit(AnomAuthfit,n.subjects=100)
#' 
#' ctIndplot(datawide = dwide,n.subjects = 100,n.manifest = 2,Tpoints = 4)
#' ctIndplot(datawide = AnomAuth+rnorm(length(AnomAuth),n.subjects = 100,n.manifest = 2,Tpoints = 4)
#' 
ctGenerateFromFit<-function(fit,n.subjects=100,timestep=.1,timerange='asdata',
  predictorSubjects='all',...){
s=summary(fit,verbose=TRUE)

gm=fit$ctmodelobj



#fix ctmodel matrices to fitted matrices
gm$LAMBDA=s$LAMBDA
gm$DRIFT=s$DRIFT
gm$DIFFUSION=t(chol(Matrix::nearPD(s$DIFFUSION+diag(1e-8,gm$n.latent))$mat))
gm$CINT=s$CINT
gm$T0MEANS=s$T0MEANS
gm$MANIFESTMEANS=s$MANIFESTMEANS
gm$T0VAR=t(chol(Matrix::nearPD(s$T0VAR)$mat))
gm$MANIFESTVAR=t(chol(Matrix::nearPD(s$MANIFESTVAR+diag(1e-8,gm$n.manifest))$mat))

if(!is.null(gm$TRAITVAR)) { #adjust traitvar from asymptotic form to cint variance form
  gm$TRAITVAR<- (gm$DRIFT) %*% s$TRAITVAR %*% t(gm$DRIFT)
  gm$TRAITVAR=t(chol(Matrix::nearPD(gm$TRAITVAR+diag(1e-8,gm$n.latent))$mat))
}

if(!is.null(gm$MANIFESTTRAITVAR)) gm$MANIFESTTRAITVAR=t(chol(Matrix::nearPD(s$MANIFESTTRAITVAR+diag(1e-8,gm$n.manifest))$mat))

if(gm$n.TDpred > 0) gm$TDPREDEFFECT=s$TDPREDEFFECT
if(gm$n.TIpred > 0) gm$TIPREDEFFECT=s$TIPREDEFFECT



if(!is.null(fit$mxobj$expectation$P0)) { #if fit with kalman filter then data needs rearranging
  dat=suppressMessages(ctLongToWide(datalong = fit$mxobj$data$observed,
  id = 'id',time = 'dT1',manifestNames = gm$manifestNames,TDpredNames = gm$TDpredNames,
  TIpredNames = gm$TIpredNames))
  
  dat=dat[,-which(colnames(dat)=='T0')]
  colnames(dat)[colnames(dat) %in% c(paste0('T',1:(gm$Tpoints-1)))]=paste0('dT',1:(gm$Tpoints-1))
} else dat=fit$mxobj$data$observed

if(predictorSubjects=='all') predictorSubjects=1:(nrow(dat))

if(timerange=='asdata') timerange=c(0,max(apply(dat[,paste0('dT',1:(gm$Tpoints-1))],1,sum,na.rm=TRUE)))

gm$Tpoints=length(seq(timerange[1],timerange[2],timestep))



out=c()
for(i in 1:n.subjects){
  predSub=sample(x = predictorSubjects,size = 1)
  if(gm$n.TDpred > 0) gm$TDPREDMEANS=matrix(dat[predSub,paste0(gm$TDpredNames,'_T',rep(0:(gm$Tpoints-1),each=gm$n.TDpred))],ncol=1)
  if(gm$n.TIpred > 0) gm$TIPREDMEANS=matrix(dat[predSub,paste0(gm$TIpredNames),drop=FALSE],ncol=1)
  
  new=suppressMessages(ctGenerate(ctmodelobj = gm,n.subjects = 1,dtmean=timestep,...))
  new[,'id']=i
   out=rbind(out,new)
   # if(i==1 & n.subjects > 1) out=rbind(out,matrix(NA,nrow=nrow(out)*(n.subjects-1),ncol=ncol(out))) #preallocate
}
return(out)
}


