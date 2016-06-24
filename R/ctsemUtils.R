# helper function to generate an index matrix, or return unique elements of a matrix
indexMatrix<-function(dimension,symmetrical=FALSE,upper=FALSE,lowerTriangular=FALSE, sep=NULL,starttext=NULL,endtext=NULL,
  unique=FALSE,rowoffset=0,coloffset=0,indices=FALSE,diagonal=TRUE,namesvector=NULL){
  if(is.null(namesvector)) namesvector=1:9999
  if(indices==T) sep<-c(",")
  tempmatrix<-matrix(paste0(starttext,namesvector[1:dimension+rowoffset],sep,rep(namesvector[1:dimension+coloffset],each=dimension),endtext),nrow=dimension,ncol=dimension)
  if(upper==TRUE) tempmatrix<-t(tempmatrix)
  if(symmetrical==TRUE) tempmatrix[col(tempmatrix)>row(tempmatrix)] <-t(tempmatrix)[col(tempmatrix)>row(tempmatrix)]
  if(unique==TRUE && symmetrical==TRUE) tempmatrix<-tempmatrix[lower.tri(tempmatrix,diag=diagonal)]
  if(lowerTriangular==TRUE) tempmatrix[col(tempmatrix) > row(tempmatrix)] <- 0
  if(indices==T){
    tempmatrix<-matrix(c(unlist(strsplit(tempmatrix,","))[seq(1,length(tempmatrix)*2,2)],
      unlist(strsplit(tempmatrix,","))[seq(2,length(tempmatrix)*2,2)]),ncol=2)
  }
  return(tempmatrix)
}




#' ctWideNames
#' sets default column names for wide ctsem datasets. Primarily intended for internal ctsem usage.
#' @param n.manifest number of manifest variables per time point in the data.
#' @param Tpoints Maximum number of discrete time points (waves of data, or measurement occasions) 
#' for an individual in the input data structure.
#' @param n.TDpred number of time dependent predictors in the data structure.
#' @param n.TIpred number of time independent predictors in the data structure.
#' @param manifestNames vector of character strings giving column names of manifest indicator variables
#' @param TDpredNames vector of character strings giving column names of time dependent predictor variables
#' @param TIpredNames vector of character strings giving column names of time independent predictor variables
#' @export

ctWideNames<-function(n.manifest,Tpoints,n.TDpred=0,n.TIpred=0,manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){
  
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  manifestnames<-paste0(manifestNames,"_T",rep(0:(Tpoints-1),each=n.manifest))
  if(n.TDpred > 0 && Tpoints > 1) {
      TDprednames<-paste0(TDpredNames,"_T",rep(0:(Tpoints-2),each=n.TDpred))
  } else {
      TDprednames<-NULL
  }
  if (Tpoints > 1) {
      intervalnames<-paste0("dT",1:(Tpoints-1))
  } else {
      intervalnames <- NULL
  }
  if(n.TIpred>0) TIprednames <- paste0(TIpredNames) else TIprednames <- NULL
  return(c(manifestnames,TDprednames,intervalnames,TIprednames))
}

# generates more complex sequences than seq
cseq <- function(from, to, by){
  temp<-c()
  for(i in from){
    temp<-c(temp,seq(i,to,by))
  }
  temp<-sort(temp)
  return(temp)
}

