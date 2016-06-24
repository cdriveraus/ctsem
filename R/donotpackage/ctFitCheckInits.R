ctModelSelectInits<-function(ctmodelobj,datawide,paramValues,paramMatrix,...){
#   require(iterpc)
  
  n.params<-length(ctmodelobj[[paramMatrix]] #length of paramMatrix...
    [is.na(suppressWarnings(as.numeric(ctmodelobj[[paramMatrix]])))])
#   
#   iterates<-iterpc(n=length(paramValues), r=n.freeparams,  replace=TRUE,ordered=TRUE)
#   iterates<-getall(iterates)
#   
#   paramCombos<-matrix(apply(iterates,1,function(x){
#     paramValues[x]
#   }),nrow=n.freeparams)
#   message(paste0('Fitting over ',ncol(paramCombos),' fixed parameter combinations'))
  
#   paramCombos<-combn(paramValues, #generate all possible combinations of parameter values
#     length(ctmodelobj[[paramMatrix]] #length of paramMatrix...
#       [is.na(suppressWarnings(as.numeric(ctmodelobj[[paramMatrix]])))] #where the elements are free...
# #       [row(ctmodelobj[[paramMatrix]])==col(ctmodelobj[[paramMatrix]])] #and on the diagonal
#   ))
    
    
#     offdiagParamCombos<-combn(offdiagParamValues, #generate all possible combinations of parameter values
#       length(ctmodelobj[[paramMatrix]] #length of paramMatrix...
#         [is.na(suppressWarnings(as.numeric(ctmodelobj[[paramMatrix]])))] #where the elements are free...
#         [row(ctmodelobj[[paramMatrix]])!=col(ctmodelobj[[paramMatrix]])] #and not on the diagonal
#       )
    
    
  ll<-c()

  i<-0
  for(param in 1:n.params){ #apply values to model and estimate with parameters fixed
    for(pvalue in 1:length(paramValues)){
      ctmodelobjnew<-ctmodelobj
      i<-i+1

    message(paste0('Fixing ', 
      ctmodelobjnew[[paramMatrix]] [is.na(suppressWarnings(as.numeric(ctmodelobjnew[[paramMatrix]])))] [param] ,
      ' to ', paramValues[pvalue]))
      
#     fixedparamnames<-ctmodelobjnew[[paramMatrix]] [is.na(suppressWarnings(as.numeric(ctmodelobjnew[[paramMatrix]])))]
#     message(paste0(fixedparamnames, ': ', paramCombos[,i], '\n'))
    
    ctmodelobjnew[[paramMatrix]] [is.na(suppressWarnings(as.numeric(ctmodelobjnew[[paramMatrix]])))] [param] <-paramValues[pvalue]
    
#      if(length(ll[!is.na(ll)]) > 0)  ctmodelobjnew$inits<-bestinits #if we have model fits already, use the best init values

    
    fit<-ctFit(datawide, ctmodelobjnew,...)

      ll[i]<-NA
      if(!is.null(fit$mxobj$output$minimum)) {
        ll[i]<-fit$mxobj$output$minimum
        if(ll[i] == min(ll,na.rm=T)) {
          bestinits<-ctGetInits(fit)
        bestinits<-rbind(bestinits, matrix(c(
          ctmodelobj[[paramMatrix]] [is.na(suppressWarnings(as.numeric(ctmodelobj[[paramMatrix]])))] [param],
          paramValues[param]),ncol=2))
        }
      }
  }
  }
  
  #output model object with inits from lowest ll parameter combinations
  if(!exists('bestinits')) stop('All combinations failed to fit!')
  message(paste0('Best fit value: ', min(ll,na.rm=T)))
#   inits<-cbind(bestinits,
#     matrix(c(ctmodelobj[[paramMatrix]][is.na(suppressWarnings(as.numeric(ctmodelobj[[paramMatrix]])))],
#     paramCombos[,which(ll==min(ll,na.rm=T))]),ncol=2) )
  ctmodelobj$inits<-bestinits
    return(ctmodelobj)
}
  