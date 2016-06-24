#' simulate panel data
#' @examples 
#' #Bivariate model, no measurement error, no trait effect
#' bivd<-generate(lag1=matrix(c(.6,.3,.1,.8),nrow=2,ncol=2),lag0=matrix(c(1,0,0,1),nrow=2),
#' n.variables=2,n.manifest=2,n.TDpred=0,subjects=1,Tpoints=5,interpoints=1,intercepts=c(0,0),
#' slopes=c(0,0),burnin=0,linearvar=c(0,0))
#' 
#' #Univariate model, N=1, measurement error, no trait effect, continuous intercept
#' unid<-generate(lag1=matrix(.8,nrow=1,ncol=1),lag0=matrix(c(1),nrow=1),
#' n.variables=1,n.manifest=1,subjects=1,Tpoints=50,intercepts=c(0),
#' slopes=c(2),burnin=0,linearvar=c(.2))
#' 
generate<-function(n.variables,n.manifest,n.TDpred=0,n.TIpred=0,subjects=3000,Tpoints=5,intercepts=rep(0,each=n.variables),
  variance=matrix(0,nrow=n.variables,ncol=n.variables),burnin=0,
  lag1=matrix(0,nrow=n.variables,ncol=n.variables),lag0=matrix(0,nrow=n.variables,ncol=n.variables),
  lagx=matrix(0,nrow=n.variables,ncol=n.variables),lagxorder=1,
  MA1=matrix(0,nrow=n.variables,ncol=n.variables),Qlag1=matrix(0,nrow=n.variables,ncol=n.variables),
  TIpredeffect=matrix(0,nrow=n.TIpred,ncol=n.variables),TIpredvar=c(0,length=n.TIpred),
  slopes=rep(0,each=n.variables),linearvar=matrix(0,nrow=n.variables,ncol=n.variables),intervals=T){
  `[` <- function(..., drop=FALSE) base::`[`(...,drop=drop) #to set DROP=FALSE on all bracket subset operations
  
  Tpointsoriginal<-Tpoints
  Tpoints<-Tpoints+burnin
  
subjects<-subjects+1 #to avoid N=1 problems, remove a row later

  
  if(n.TIpred>0){
    TIpredvar<-matrix(rnorm(subjects*n.TIpred,0,rep(TIpredvar,each=subjects)),nrow=subjects, ncol=n.TIpred  )
#     TIpredeffect<-matrix(rep(TIpredeffect,each=subjects),nrow=subjects,ncol=n.TIpred*n.variables)
  }
  
  intercepts<-matrix(1,nrow=subjects) %*% intercepts 
  
  slopes<-matrix(slopes,nrow=subjects,ncol=n.variables,byrow=T)
  if(sum(diag(linearvar))>0) {
    slopes<-matrix(rnorm(subjects*n.variables,0,1),nrow=subjects,ncol=n.variables) %*%  chol(linearvar)+slopes    
  }
    
  y<-matrix(0,nrow=subjects,ncol=Tpoints*n.variables)
  errormatrix<-matrix(rnorm((subjects*Tpoints*n.variables),0,1),nrow=subjects,ncol=Tpoints*n.variables)
  innovmatrix<-matrix(rnorm((subjects*Tpoints*n.variables),0,1),nrow=subjects,ncol=Tpoints*n.variables)
  y[,1:n.variables]<-intercepts+(innovmatrix[,1:n.variables]) %*% chol(lag0) #T1
  yslopesonly<-y
  yslopesonly[,1:n.variables]<-intercepts
  
  #browser()
  for(j in 1:(Tpoints-1)){ 
    
    #set innovation matrix
    innovmatrix[,(j*n.variables+1):(j*n.variables+n.variables)]<-innovmatrix[,(j*n.variables+1):(j*n.variables+n.variables)]%*%chol(lag0)+
      innovmatrix[,((j-1)*n.variables+1):((j-1)*n.variables+n.variables)]%*%Qlag1
    
    y[,(j*n.variables+1):(j*n.variables+n.variables)]<-slopes+ 
      y[,(j*n.variables+1-n.variables):(j*n.variables)]%*%t(lag1)+ #lag1
      innovmatrix[,(j*n.variables+1):(j*n.variables+n.variables)]+ #instant effects
      innovmatrix[,((j-1)*n.variables+1):((j-1)*n.variables+n.variables)]%*%sqrt(lag0)%*%(MA1)
    
    if(all((  (j*n.variables+1-(n.variables*(lagxorder-1))) : (j*n.variables-(n.variables*(lagxorder-2)))  )>0)){
      y[,(j*n.variables+1):(j*n.variables+n.variables)]<-y[,(j*n.variables+1):(j*n.variables+n.variables)]+
        y[,(j*n.variables+1-(n.variables*(lagxorder-1))):(j*n.variables-(n.variables*(lagxorder-2)))]%*%lagx #lagx
    }
    
  }
  
  if(n.TIpred>0) {
    y[,1:(n.variables*Tpoints)]<- #TIpredictors
      y[,1:(n.variables*Tpoints)] + rep(TIpredvar%*%(TIpredeffect),times=Tpoints) #needs fixing
    
  }
  
  #errors
  for(j in 0:(Tpoints-1)){ 
    for (i in 1:n.variables){
      y[,(j*n.variables+i)]<-y[,(j*n.variables+i)]+
        errormatrix[,(j*n.variables+1):((j+1)*n.variables),drop=F] %*% sqrt(variance[,i]) #error variance 
    }
  }
    
  if(burnin>0) y<-y[,-0:-(burnin*n.variables)]
  
  Tpoints<-Tpointsoriginal
  ymanifest<-matrix(y[  ,cumsum(c(rep(1,times=n.manifest),
    (rep(c((n.TDpred+1),rep(1,times=n.manifest-1)),times=Tpoints-1))))  ],nrow=nrow(y))
  colnames(ymanifest)<-c(paste0("V",1:(n.manifest*Tpoints)))
  
 if(intervals==TRUE){
   ymanifest<-cbind(ymanifest,matrix(1,ncol=Tpoints-1,nrow=nrow(ymanifest)))
  colnames(ymanifest)<-c(paste0("V",1:(n.manifest*Tpoints)),paste0("I",1:(Tpoints-1)))
 }
  yout<-ymanifest
  
  if(n.TDpred>0){  
    yexogenous<-NULL
    for(i in 1:n.TDpred){
      temp<-y[,cumsum(c(n.manifest+n.TDpred+n.manifest+i,(rep(n.manifest+n.TDpred,times=Tpoints-2))))]
      colnames(temp)<-c(paste0("P",i,"_",1:(Tpoints-1)))
#       if(intervals==TRUE) { #removed as predictors now use regular intervals
#         temp<-cbind(temp,matrix(1,nrow=nrow(temp),ncol=ncol(temp)))
#       colnames(temp)<-c(paste0("P",i,"_",1:(Tpoints-1)),paste0("PI",i,"_",1:(Tpoints-1)))    #change 1 to i when interval names sorted
#       }
      yexogenous<-cbind(yexogenous,temp)
    }
    yout<-cbind(yout,yexogenous)
  }
  
  if(n.TIpred>0) {
    colnames(TIpredvar)<-paste0("Z",1:n.TIpred)
    yout<-cbind(yout,TIpredvar)
  }
  
  y<-as.data.frame(yout[-1,,drop=FALSE])
  
  return(y)
  
}
