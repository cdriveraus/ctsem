#' Simulate continuous time data
#' 
#' This function generates data according to the specified ctsem model object. 
#' Not all T0 matrices are included at present,
#' safest to use a high burnin (where 'high' is sufficient for the process to forget starting values)
#' 
#' 
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dT Time interval (delta T) to use, defaults to 1.
#' @param asymptotes Are the parameters provided asymptotic paramters, or the regular continuous time parameters?
#' @details TRAITTDPREDCOV and TIPREDCOV matrices are not accurately accounted for, at present.
#' @examples 
#' #generate data for 2 process model, each process measured by noisy indicator, 
#' #stable individual differences in process levels.
#' 
#' generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'  MANIFESTVAR=diag(.1,2),
#'  LAMBDA=diag(1,2),
#'  DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
#'  TRAITVAR=matrix(c(.5,.2,0,.8),nrow=2),
#'  DIFFUSION=matrix(c(1,.2,0,4),2),
#'  CINT=matrix(c(1,0),nrow=2),
#'  T0MEANS=matrix(0,ncol=1,nrow=2),
#'  T0VAR=diag(1,2))
#'
#' data<-ctGenerate(generatingModel,n.subjects=150,burnin=500)
#'
#' ## Further examples set to 'dontrun' because they take longer than 5s. 
#' 
#' \dontrun{
#' ctIndplot(data,n.manifest=2,Tpoints=4,n.subjects=10)
#'
#' model<-ctModel(Tpoints=8, TRAITVAR='auto', n.latent=2,
#'  n.manifest=2, LAMBDA=diag(2))
#'
#' checkf<-ctFit(data,model,stationary=c('T0VAR','T0MEANS'))
#' summary(checkf,verbose=TRUE)
#'  
#'  
#' #### with 2 process from 4 indicators, latent trait, TDpred and TIpred 
#' Tpoints=8
#' n.latent=2
#' n.manifest=4
#' n.TDpred=1
#' n.TIpred=1
#' 
#' generatingModel<-ctModel(Tpoints=Tpoints,n.latent=n.latent,
#'  n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#'  LAMBDA=matrix(c(1,.4,.8,0,0,0,0,1),nrow=n.manifest,ncol=n.latent),
#'  MANIFESTVAR=diag(c(.2),n.manifest),
#'  MANIFESTTRAITVAR=matrix(c(.3,.1,0,.2,  0,.2,0,.15, 
#'  0,0,0,0, 0,0,0,.4) ,n.manifest,n.manifest),
#'  MANIFESTMEANS=matrix(c(0,0,0,0),n.manifest,1),
#'  DRIFT=matrix(c(-.23,.1,.0,-.4),n.latent),
#'  DIFFUSION=matrix(c(8.3,5.1,0,8.4),n.latent,n.latent),
#'  CINT=matrix(c(0,.4),n.latent,1),
#'  TDPREDEFFECT=matrix(c(1.2,-.4),nrow=n.latent,ncol=n.TDpred),
#'  TIPREDEFFECT=matrix(c(.32,-.08),nrow=n.latent,ncol=n.TIpred),
#'  TDPREDMEANS=matrix(c(0,0,1,rep(0, (Tpoints-1-3)*n.TDpred)),nrow=n.TDpred*(Tpoints-1)),
#'  TIPREDMEANS=matrix(0,nrow=n.TIpred),
#'  TDPREDVAR=diag(0.1,n.TDpred*(Tpoints-1)),
#'  TIPREDVAR=diag(.4,n.TIpred),
#'  T0MEANS=matrix(0,ncol=1,nrow=n.latent))
#'  
#' data<-ctGenerate(generatingModel,n.subjects=100,burnin=500)
#' 
#' model<-ctModel(n.latent=n.latent, n.TDpred=n.TDpred, n.TIpred=n.TIpred, 
#' n.manifest=n.manifest, 
#' LAMBDA=matrix(c(1,"l2","l3",0,0,0,0,1),n.manifest,n.latent),  
#' TRAITVAR='auto',Tpoints=Tpoints)
#' 
#' 
#' fit<-ctFit(data,model, stationary=c('T0VAR','T0MEANS', 'T0TIPREDEFFECT'))
#' summary(checkf)  
#' }
#' @export

ctGenerate<-function(ctmodelobj,n.subjects=1000,burnin=0,dT=1,asymptotes=FALSE,simulTDpredeffect= FALSE){
  
  
  
  
  ###read in model
  for(i in 1:length(ctmodelobj)){ #this loop reads in the specified continuous time model
    assign(names(ctmodelobj[i]),ctmodelobj[[i]])

    if(is.matrix(ctmodelobj[[i]])){ #if this element is a matrix, continue on...
     
    if(any(is.na(suppressWarnings(as.numeric(get(names(ctmodelobj[i]))))))){ #if it contains character labels
      assign(names(ctmodelobj[i]),matrix(0,nrow=nrow(get(names(ctmodelobj[i]))), #set the values to 0 instead
                                      ncol=ncol(get(names(ctmodelobj[i])))))
      message(paste0(names(ctmodelobj[i])," contained character labels - setting matrix to 0"))
    }
    
     #set any matrices to numeric elements
      assign(names(ctmodelobj[i]), 
        matrix(as.numeric(get(names(ctmodelobj[i]))),nrow=nrow(get(names(ctmodelobj[i]))),
                                      ncol=ncol(get(names(ctmodelobj[i])))))
    }
  }
  
  # #lower triangular transform
  # for(tempmatname in c('T0VAR','MANIFESTVAR', 'DIFFUSION', 'TRAITVAR','MANIFESTTRAITVAR','TDPREDVAR','TIPREDVAR')){
  # 
  #   tryCatch(assign(tempmatname,get(tempmatname) %*% t(get(tempmatname))), error=function(e) {
  #     assign(tempmatname,NULL)})
  # }
  # 
  
  #set up extra matrices
  DRIFTHATCH <- DRIFT %x% diag(n.latent) + diag(n.latent) %x% DRIFT #generate drifthatch
  if(asymptotes==FALSE) dynresidualcov <- matrix(solve(DRIFTHATCH) %*% ((OpenMx::expm(DRIFTHATCH %x% dT)) - #generate dynamic error cov from continuous value
                                                  diag(1,n.latent^2)) %*% OpenMx::rvectorize(DIFFUSION %*% t(DIFFUSION)),nrow=n.latent)
  if(asymptotes==TRUE) dynresidualcov <- matrix((diag(n.latent^2) - OpenMx::expm(DRIFT %x% dT) %x% OpenMx::expm(DRIFT %x% dT)) %*% c(DIFFUSION),nrow=n.latent)
  
  
  if(!all(is.numeric(T0VAR))) { #if T0VAR does not have all values fixed
    print("No T0VAR specified - generated process will be at equilibrium")
    T0VAR<-diag(1,n.latent) #arbitrarily set it
    if(burnin < 200){ #and if burnin has not been set
      burnin <- 500 #ensure burnin is high enough that arbitrary phit1 doesn't matter
    }
  }
  
  if(!all(is.numeric(T0MEANS))) T0MEANS<-matrix(0,ncol=n.latent) #if T0MEANS has not been fixed arbitrarily set it

  
  
  
  
  ####traits
  traiteffect<-matrix(0,nrow=n.subjects,ncol=n.latent) #create traiteffect matrix with 0 effect
  trait<-matrix(0,nrow=n.subjects,ncol=n.latent)
  if(!is.null(TRAITVAR[1])) { #if traits are specified
 trait <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=TRAITVAR %*% t(TRAITVAR),tol=1) #generate trait effects    
    traiteffect<- trait %*% t(LAMBDA)
 }
  
  #####predictors
  
  TDpreds<-matrix(NA,nrow=n.subjects,ncol=n.TDpred*(Tpoints-1))
  TDpredeffects <- matrix(0,nrow=n.subjects,ncol=n.latent*(Tpoints-1)) #create 0 TDpredeffects in case no TDpreds
  if (n.TDpred>0) { #but if TDpreds exist
    if(simulTDpredeffect==FALSE) TDpredparam <- OpenMx::expm(DRIFT %x% dT) %*%  TDPREDEFFECT #calculate effect size
    if(simulTDpredeffect==TRUE) TDpredparam <-TDPREDEFFECT #calculate effect size
    TDpreds <- MASS::mvrnorm(n=n.subjects,mu=TDPREDMEANS, #generate TDpred variables from TDPREDMEANS and TDPREDVAR
                       Sigma=TDPREDVAR %*% t(TDPREDVAR) ,tol=1)
    
    TDpreds <- TDpreds + trait %*% TRAITTDPREDCOV
    
    for(l in 1:n.latent){ #for each latent process
      for(i in 1:(Tpoints-1)){#for each latent state variable except the final ones
        #     TDpredeffects[,(i-1)*n.latent+l] <-  TDpreds[,seq(i,n.TDpred*(Tpoints-1),(Tpoints-1))] %*% t(TDpredparam[l, ]) #create TDpred effects
        TDpredeffects[,(i-1)*n.latent+l] <-  TDpreds[,seq(i,n.TDpred*(Tpoints-1),(Tpoints-1))] %*% t(TDpredparam[l, ,drop=F]) #create TDpred effects
      }
    }
  }    
  
  TIpreds<-matrix(NA,nrow=n.subjects,ncol=n.TIpred)
  TIpredeffects<-matrix(0,nrow=n.subjects,ncol=n.latent) #set 0 effects in case no TIpreds
  if (n.TIpred>0) { #if TIpreds exist
    if(asymptotes==FALSE) TIPREDEFFECTdiscrete <-  solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1, n.latent)) %*%  TIPREDEFFECT #calculte their effect
    if(asymptotes==TRUE) TIPREDEFFECTdiscrete <-  (diag(1, n.latent) - OpenMx::expm(DRIFT %x% dT) ) %*%  TIPREDEFFECT #calculate their effect
    
    TIpreds <- matrix(
      MASS::mvrnorm(n=n.subjects,
              mu=TIPREDMEANS, #generate TIpreds
              Sigma=TIPREDVAR %*% t(TIPREDVAR), tol=1)
      ,nrow=n.subjects)
      #     TIpredeffects<-matrix(TIpreds*rep(TIpredparam,n.subjects*n.latent),nrow=n.subjects,ncol=n.TIpred) #generate effects on the latents
      TIpredeffects<-matrix(TIpreds %*% t(TIPREDEFFECTdiscrete),nrow=n.subjects,ncol=n.latent) #generate effects on the latents
  }
  
  ####cint
  
  cinteffect <- matrix(solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1, n.latent)) %*% 
                         CINT,byrow=T,nrow=n.subjects,ncol=n.latent) #continuous intercept effect
  
  
  
  Tpoints<-Tpoints+burnin #add burnin to Tpoints (after we checked if extra burnin was needed for T0 cov, and after predictor generation)
  T0VAReffect<-MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=(T0VAR %*% T0VAR),tol=1) #create effect of non-trait variance at T0
  
  
  
  
  
  
  ########create latents  
  latents <- matrix(,nrow=n.subjects,ncol=n.latent*Tpoints) #create latent matrix
  latents[,1:n.latent] <- matrix(T0MEANS,nrow=n.subjects,ncol=n.latent,byrow=T)+ T0VAReffect #create T0 latents including all possible effects
  
  #burnin
  if(burnin>0){
    for(i in seq(n.latent+1,n.latent*(burnin+1),n.latent)){ #for every time block of latents after the first
      drifteffect <- latents[,(i-n.latent):(i-1)] %*% t(OpenMx::expm(DRIFT %x% dT))#effect of past time points    
      qeffect <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=dynresidualcov,tol=1) #effect of q noise
      
      latents[,i:(i+n.latent-1)]<-drifteffect+cinteffect+qeffect+TIpredeffects#sum of all constant effects at t to create latents
    }

    
    latents<-as.matrix(latents[,-1:-(burnin*n.latent),drop=FALSE]) #remove burnin from latents
    Tpoints<-Tpoints-burnin #remove burnin from Tpoints
  }
  
  
  #post burnin latents (here TDpredictors are added)
  for(i in seq(n.latent+1,n.latent*Tpoints,n.latent)){ #for every time block of latents after the first
    drifteffect <- latents[,(i-n.latent):(i-1)] %*% t(OpenMx::expm(DRIFT %x% dT))#effect of past time points    
    qeffect <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=dynresidualcov,tol=1) #effect of q noise
    latents[,i:(i+n.latent-1)]<-drifteffect+cinteffect+qeffect+TIpredeffects+
      TDpredeffects[,(i-n.latent) : (i-1)]  #sum of all effects at t to create latents
  }
  
  
  #generate indicators from latents

  manifests<-matrix(0,nrow=n.subjects,ncol=n.manifest*Tpoints) #create manifests matrix 
  
  #   
  if(is.null(MANIFESTTRAITVAR[1])) MANIFESTTRAITVAR <- diag(0,n.manifest) #generate 0 matrix if needed
  manifesttraiteffects<-MASS::mvrnorm(n.subjects,mu=rep(0,n.manifest),Sigma= MANIFESTTRAITVAR %*% t(MANIFESTTRAITVAR))

  for(i in 1:Tpoints){
    manifests[,((i-1)*n.manifest+1):(i*n.manifest)] <- (latents[,((i-1)*n.latent+1):(i*n.latent)] + trait) %*% t(LAMBDA) + 
      MASS::mvrnorm(n.subjects,mu=MANIFESTMEANS,Sigma=MANIFESTVAR %*% t(MANIFESTVAR)) + # manifest means and error variance
      manifesttraiteffects # manifest traits
  }  

  
  intervals<-matrix(dT,nrow=n.subjects,ncol=(Tpoints-1)) #add intervals to latent output
  out<-as.matrix(cbind(manifests,TDpreds,intervals,TIpreds)) #output variables
  colnames(out)<-ctWideNames(Tpoints=Tpoints,n.manifest=n.manifest,n.TIpred=n.TIpred,n.TDpred=n.TDpred)
  rownames(out)<-1:nrow(out)
  return(out)
}
