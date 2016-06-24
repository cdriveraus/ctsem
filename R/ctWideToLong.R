#' ctWideToLong
#' Convert ctsem wide to long format
#' @param datawide ctsem wide format data
#' @param Tpoints number of measurement occasions in data
#' @param n.manifest number of manifest variables
#' @param n.TDpred number of time dependent predictors
#' @param n.TIpred number of time independent predictors
#' @param manifestNames Character vector of manifest variable names.
#' @param TDpredNames Character vector of time dependent predictor names.
#' @param TIpredNames Character vector of time independent predictor names.
#' @examples 
#'  #First load the example ctsem wide format data with absolute times
#'  data('datastructure')
#'  datastructure #contains two time intervals (dTx), therefore 3 time points.
#'  #Then convert to long format
#'  longexample <- ctWideToLong(datawide = datastructure, Tpoints=3, 
#'  n.manifest=3, manifestNames = c("Y1", "Y2", "Y3"),
#'  n.TDpred=1, TDpredNames = "TD1", 
#'  n.TIpred=2, TIpredNames = c("TI1", "TI2"))
#'
#'  #Then convert the time intervals to absolute time
#'  long <- ctDeintervalise(datalong = longexample, id='id', dT='dT')
#'  long
#'
#' 
#' @export

ctWideToLong<-function(datawide,Tpoints,n.manifest,n.TDpred=0, n.TIpred=0, 
  manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){ 
  
  #names
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
  
  datawide<-as.matrix(datawide,nrow=nrow(datawide),ncol=ncol(datawide)) #set to matrix for manipulations
  n.subjects<-nrow(datawide) #calculate number of subjects in datawide
  
  
#   datalong<-matrix(NA,nrow=n.subjects*Tpoints,ncol=1+n.manifest+1+n.TDpred+n.TIpred) #create blank datalong matrix
#   colnames(datalong)<-paste0("name",1:ncol(datalong)) #create default colnames
#   colnames(datalong)[1:(1+n.manifest)]<-c("id",manifestNames) #set column names for manifests
#   colnames(datalong)[ncol(datalong)-n.TIpred]<-"dT" #set colnames for time intervals
  
  # if(n.TDpred>0)  colnames(datalong)[(1+n.manifest+1):(1+n.manifest+n.TDpred)]<-TDpredNames #set TDpredictor colnames    
  # if(n.TIpred>0)  colnames(datalong)[(ncol(datalong)-(n.TIpred-1)):ncol(datalong)]<-TIpredNames #set TIpredictor colnames 
 
  
  manifests<-matrix(t(datawide[,1:(n.manifest*Tpoints),drop=FALSE]),byrow=T,nrow=n.subjects*Tpoints)
  colnames(manifests)<-manifestNames
  
  times<-matrix(t(cbind(0,datawide[,(Tpoints*n.manifest+(Tpoints-1)*n.TDpred+1) : 
      (Tpoints*n.manifest+(Tpoints-1)*n.TDpred+ (Tpoints-1)),drop=FALSE])),byrow=T,nrow=n.subjects*Tpoints)
  colnames(times)<-'dT'
  
  id<-matrix(rep(1:n.subjects,each=Tpoints),ncol=1)
  colnames(id)<-'id'
  
  datalong<-cbind(id,times,manifests)
  
  if(n.TDpred>0) {
  tdpreds<-matrix(NA,nrow=n.subjects*(Tpoints-1),ncol=n.TDpred)
  for(tdpredi in 1:n.TDpred){
  tdpreds[,tdpredi]<-t(datawide[,(Tpoints*n.manifest+(Tpoints-1)*(tdpredi-1)+1) : (Tpoints*n.manifest+(Tpoints-1)*(tdpredi)),drop=FALSE])
  }
  tdpredsfull<-matrix(NA,nrow=n.subjects*Tpoints,ncol=n.TDpred)
  tdpredsfull[-seq(Tpoints,n.subjects*Tpoints,Tpoints),]<-tdpreds
    colnames(tdpredsfull)<-TDpredNames
    datalong<-cbind(datalong,tdpredsfull)
  }
  
  if(n.TIpred > 0){
   
  tipreds<- matrix(t(datawide[,(Tpoints*n.manifest+(Tpoints-1)*n.TDpred+(Tpoints)) : 
      (Tpoints*n.manifest+(Tpoints-1)*n.TDpred+ Tpoints -1 + n.TIpred),drop=FALSE]),byrow=T,nrow=n.subjects)
  tipredsfull<-matrix(NA,nrow=n.subjects*Tpoints,ncol=n.TIpred)
  tipredsfull<-tipreds[rep(1:n.subjects,each=Tpoints),,drop=F]
  colnames(tipredsfull)<-TIpredNames
  datalong<-cbind(datalong,tipredsfull)
  }
    
  
  
#   
#   pb<-txtProgressBar(1,n.subjects)
#   for(j in 1:n.subjects){ #for each subject
#     setTxtProgressBar(pb,j)
#     for(i in 1:Tpoints){ #at each time point
#       
#       datalong[(j-1)*Tpoints+i, 1:(1+n.manifest)]<-
#         c(j, #add a subject indicator /id to datalong
#           datawide[j,(i-1)*n.manifest+(1:n.manifest)]) #add manifests to datalong
#       
#       if(i==1) datalong[(j-1)*Tpoints+i,(ncol(datalong)-n.TIpred)]<-0 #if first obs, add interval of 0
#       
#       if(i!=1) datalong[(j-1)*Tpoints+i,(ncol(datalong)-n.TIpred)]<- #if not first obs, set time interval to
#         datawide[j,Tpoints*n.manifest+(Tpoints-1)*n.TDpred+(i-1)] #observed interval
#       
#       if(n.TDpred>0) for(p in 1:n.TDpred){ #if time dependent predictors exist, for each:
#         datalong[(j-1)*Tpoints+i,(1+n.manifest+p)]<-c(
#           ifelse(i==Tpoints,NA,datawide[j,(n.manifest*Tpoints+(p-1)*(Tpoints-1)+(i))])) #add time dependent predictor to datalong if not last obs
#       }
#       
#       if(n.TIpred>0) { #if time independent predictors exist
#         datalong[(j-1)*Tpoints+i,(1+n.manifest+n.TDpred+1+1:n.TIpred)]<-
#           datawide[j,(n.manifest*Tpoints+(Tpoints-1)+n.TDpred*(Tpoints-1)+1:n.TIpred),drop=FALSE] #add time independent predictor to all obs
#       }        
#     }
#   }
  
  
  return(datalong)
}