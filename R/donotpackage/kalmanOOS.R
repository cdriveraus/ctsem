#'kalmanOOS
#'
#'@export
kalmanOOS<-function(mxobj, nmissing,
  missingsList='random', plot=FALSE,plotlength=50,...){
  
  manifestNames<-rownames(mxobj[[mxobj$expectation$C]]$labels)
  
  data<-mxobj$data$observed[,manifestNames,drop=F]

  observedIndices<- which(!is.na(data)) 
    
    
    require(ctsem)
    require(OpenMx)
  if(any(missingsList != 'random') & any(is.na(data[missingsList]))) stop('missingsList indexes already missing data!')
    if(all(missingsList == 'random')) missingsList<-sample(observedIndices,nmissing)
    tempdata<-data
    tempdata[missingsList] <- NA

      mxobj$data$observed[,manifestNames] <- tempdata
    mxobj<-ctmxTryHard(mxobj,initialTolerance=1e-18,
      bestInitsOutput=FALSE,
      extraTries=1,loc=1,scale=.2,paste=FALSE,...)
    
    if(class(mxobj) == 'try-error' || mxobj$output$status[[1]] > 1 | mxobj$output$status$status == -1) fit <-NA
    
    if(mxobj$output$status[[1]] <= 1){
    predictions<-kalmanPredict(mxobj,mxstyle=T)$manifestForecast
    
   browser()
    fit<- mean((predictions[missingsList] - data[missingsList])^2)
   

   
    if(plot==TRUE){
      OOSpredictions<-predictions
      OOSpredictions[-missingsList]<-NA
    plot(data[1:plotlength],type='b')
    points(predictions[1:plotlength],col='red')
    points( OOSpredictions[1:plotlength],col='blue')
    }
    
    }
  
  return(fit)
}