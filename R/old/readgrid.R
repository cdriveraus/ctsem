
#' read grid output files into single output matrix
#' @param fileprefix Name of file excluding extension
#' @param rprefix Name of R object in file, excluding numeric identifier
#' @param filelocation Defaults to c:/gridout/
#' @param keepcodereds Defaults to T, keeping model fits regardles of mxstatus of 0 or 1, or 6 etc.
readgrid<-function(fileprefix,rprefix,n.subjects,filelocation="c:/gridout/",keepcodereds=F){
  
  load(file=paste0(filelocation,fileprefix[1],".RData"))
  output<-matrix(nrow=length(fileprefix),ncol=length(get(paste0(rprefix[1]))))
  colnames(output)<-colnames(get(paste0(rprefix[1])))
  for(i in 1:length(fileprefix)){
    errorToNA(load(file=paste0(filelocation,fileprefix[i],".RData")))
    
    temp<-tryCatch(c(get(paste0(rprefix[i]))),
      error=function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        NA}
    )
    if(!is.na(temp[1])) output[i,]<-temp
  }
  rm(list=paste0(rprefix))
  # output[,1]<-output[,1]+1 #for wave correction
  print(head(output))
  print(paste0(length(which(!is.na(output[,1])))," records loaded"))
  
  # output<-cbind(rangelist[,2]-rangelist[,1],rangelist[,1],output)
  # colnames(output)[1:2]<-c("waves","startyear")
  
  #remove code reds
  if(keepcodereds==F) output<-output[which(output[,3]<2),]
  
  return(output)
}
