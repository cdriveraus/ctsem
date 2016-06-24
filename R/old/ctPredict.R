#'Plots original data against predicted - only valid for single factors without 
#'measurement error
#'@param DRIFT matrix of continuous (or discrete if discrete=T) time drift coefficients.  
ctPredict<-function(individualdatawide,ctmodelobj,pause=TRUE,n.subjects=1000){
  
  ###read in model
  for(i in 1:length(ctmodelobj)){ #this loop reads in the specified continuous time model
    assign(names(ctmodelobj[i]),eval(parse(text = paste0("ctmodelobj","$",names(ctmodelobj[i])))))
    if(any(is.na(as.numeric(eval(parse(text=names(ctmodelobj[i]))))))){
      assign(names(ctmodelobj[i]),matrix(0,nrow=nrow(eval(parse(text=names(ctmodelobj[i])))),
        ncol=ncol(eval(parse(text=names(ctmodelobj[i]))))))
      message(paste0(names(ctmodelobj[i])," contained character labels - setting matrix to 0"))
    }
  }
  

  DRIFTHATCH <- DRIFT %x% diag(n.latent) + diag(n.latent) %x% DRIFT #generate drifthatch

  
  
original<-matrix(NA,ncol=n.manifest,nrow=Tpoints)
estimated<-matrix(NA,ncol=n.manifest,nrow=Tpoints)

for(i in 1:n.manifest){
original[,i]<-individualdatawide[,seq(i,Tpoints*n.manifest,n.manifest)]
}

timing<-c(0,individualdatawide[1,(n.manifest*Tpoints+1):(n.manifest*Tpoints+Tpoints-1)])#columnise intervals

  estimatedlong<-c()
  for(subject in 1:n.subjects){  
estimated[1,]<-T0MEANS #set first wave to T0MEANS
 
for(j in 2:Tpoints){
  
  dynresidualcov <- matrix(solve(DRIFTHATCH)%*%((expm(DRIFTHATCH * timing[j])) - #generate dynamic error cov
      diag(1,n.latent^2))%*%rvectorize(DIFFUSION),nrow=n.latent)
  qeffect <- mvrnorm(1,mu=rep(0,n.latent),Sigma=dynresidualcov,tol=1) #effect of noise
  if(subject==1) qeffect<-matrix(0,ncol=n.latent) #for first estimation, estimate without dynamic residual - best guess
  estimated[j,]<-original[j-1,]%*%t(expm(DRIFT*timing[j]))+ #drift
    CINT+ #intercept
    +qeffect
}
    
    singlelong<-cbind(subject,estimated,timing)
    estimatedlong<-rbind(estimatedlong,singlelong)
}

for(i in 2:length(timing)){ #adjust timing to absolute for plot representation
  timing[i]<-sum(timing[(i-1):i])
}

  for(i in (1:n.latent)+1){ #plus one because of id column

plot(timing,original[,(i-1)],type="b",main=paste0("Y",i-1))
    colindex<-2:(ncol(estimatedlong)-1)
    
points(timing,estimatedlong[1:Tpoints,i],col="red",type="b")
    for(subject in 2:n.subjects){
      points(timing,estimatedlong[(1:Tpoints)+Tpoints*(subject-1),i],col=rgb(.9,.6,.2,.05),pch=16)
    }
mae<-mean(sqrt((estimatedlong[1:Tpoints,i]-original[,(i-1)])^2),na.rm=T)
rmse<-sqrt(mean((estimatedlong[1:Tpoints,i]-original[,(i-1)])^2,na.rm=T))
legend("topright",c("Original","Estimated",paste0("MAE = ",round(mae,3)),
                    paste0("RMSE = ",round(rmse,3))),text.col=1:4,bty="n")
if(pause==T){
  cat ("Press [enter] to show next variable")
  readline()
}
plot(timing,original[,(i-1)]-estimated[,(i-1)],type="b",main=paste0("mean residuals V",i-1))
if(pause==T && i<n.latent) {
  cat ("Press [enter] to show next variable")
  readline()
}
}
}
