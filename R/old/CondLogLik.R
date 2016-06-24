#' internal ctsem function to calculate log likelihood row by row when requested
CondLogLik=function(rammodel,state){

  data<-matrix(as.matrix(rammodel$data$observed), #[,-grep("dT",colnames(rammodel$data$observed)),drop=F]
    nrow=dim(rammodel$data$observed)[1])

  Amat<-mxEval(A,model=rammodel,compute=TRUE)
  Amat[is.na(Amat)]<-0
  Smat<-mxEval(S,model=rammodel,compute=TRUE)
  Smat[is.na(Smat)]<-0
  Fmat<-rammodel$F@values
  Mmat<-matrix(mxEval(M,model=rammodel,compute=TRUE),nrow=1)
  Mmat[is.na(Mmat)]<-0
  
  bigI<-diag(nrow(Amat))
  expcov<-Fmat%*%solve(bigI-Amat)%*%Smat%*%t(solve(bigI-Amat))%*%t(Fmat)
  expmeans<-t(Fmat%*%(solve(bigI-Amat))%*%t(Mmat))
#   expcov<-rammodel$expcov$result 
#   expmeans<-rammodel$expmeans@result 
  
  totalvalue=0
#   if(!exists("Ltotal")) Ltotal<<-0 #delete this later, for incremental plotting
  
for(i in 1:nrow(data)){
  filter<-!is.na(data[i,,drop=FALSE])
k<-length(data[i,filter,drop=FALSE])
  x<-data[i,filter,drop=FALSE]
  expcovfilt<-expcov[filter,filter,drop=FALSE]
  expmeansfilt<-expmeans[,filter,drop=FALSE]
  
  
  detexpcovfilt<-det(expcovfilt)
  if(detexpcovfilt<=0) {

  expcovfilt<-diag(k)
  detexpcovfilt<-30
  }
  
  if(any(eigen(expcovfilt)$values==0)){
    expcovfilt<-diag(k)
  }
  
  meandif <- x-expmeansfilt
  rowvalue=k*log(2*pi)+log(detexpcovfilt)+(meandif)%*%solve(expcovfilt)%*%t(meandif)
  if(is.nan(rowvalue)) browser()
  totalvalue<-totalvalue+rowvalue
}
#   Ltotal<<-c(Ltotal,LL)
#   plot(Ltotal[-1])
#   print(Ltotal)
  
   return(totalvalue)}

