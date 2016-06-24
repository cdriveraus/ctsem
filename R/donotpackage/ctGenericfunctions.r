#plots means of ctsem wide panel data
meanplot<-function(data,n.manifest,Tpoints,...){
  plot(1:Tpoints,rep(0,times=Tpoints),ylim=c(min(colMeans(data[,1:(n.manifest*Tpoints)],na.rm=T)),
    max(colMeans(data[,1:(n.manifest*Tpoints)],na.rm=T))),...)
  for(i in 1:n.manifest){
    points(1:Tpoints,colMeans(data[,seq(i,n.manifest*Tpoints,n.manifest)],na.rm=T),col=(i+1),type="b")
  }}

#Allows easy creation of matrices with a diagonal of character vectors.

chardiag<-function(x,dim=x){
  out<-diag(dim)
  diag(out)<-x
  return(out)
}

#' ctGetInits
#' 
#' Extracts estimates from a fitted ctsem model and returns in ctsem init matrix layout
#' @param ctfitobject ctsem fit object to extract new starting values from

ctGetInits<-function(ctfitobject,...){
  if(class(ctfitobject)!='ctsemFit') stop('Specified fitobject is not of class ctsemFit')
  inits<-matrix(OpenMx::omxGetParameters(ctfitobject$mxobj,...),ncol=1)
  inits<-cbind(names(OpenMx::omxGetParameters(ctfitobject$mxobj,...)),inits)
  return(inits)
}



#Remove observations at random from a wide ctsem format dataset
# @param retainpercent Percent of each row to retain. Must be between 1 and 100, retainpercent/100*Tpoints must generate a whole number.
ctRemoveObservations<-function(datawide,Tpoints,n.manifest,n.TDpred=0,n.TIpred=0,retainpercent=100,
  manifestNames="auto",TDpredNames="auto",TIpredNames='auto'){
  
  datalong<-ctWideToLong(datawide,
    Tpoints=Tpoints,
    n.manifest=n.manifest,
    n.TIpred=n.TIpred,
    n.TDpred=n.TDpred) #convert to long
  datalong<-ctDeintervalise(datalong) #convert to absolute time
  
  samplelist<-c()
  for(i in 1:nrow(datawide)){
    samplelist<-c(samplelist,
      sample(  ((i-1)*Tpoints+1) : 
          (i*Tpoints), 
        ceiling(retainpercent/100 * Tpoints))) #keep retainpercent rows of each individual at random
  }
  
  datalong<-datalong[samplelist,,drop=F] #retain only samplelist rows
  
  datalong<-datalong[order(datalong[,"id"],datalong[,"AbsTime"]),] #fix ordering of rows by time
  
  if(manifestNames=='auto') manifestNames <- paste0("Y",1:n.manifest)
  if(n.TDpred >0 && TDpredNames=='auto') TDpredNames <- paste0("TD",1:n.TDpred) 
  if(n.TDpred ==0) TDpredNames <- NULL
  if(n.TIpred >0 && TIpredNames=='auto') TIpredNames <- paste0("TI",1:n.TIpred) 
  if(n.TIpred == 0) TIpredNames <- NULL
  
  datawide<-ctLongToWide(datalong,
    id="id",
    manifestNames=manifestNames,
    TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,
    time="AbsTime") 
  
  datawide<-ctIntervalise(datawide,
    manifestNames=manifestNames,
    TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,
    Tpoints=Tpoints*retainpercent/100,
    n.manifest=n.manifest,
    n.TDpred=n.TDpred,
    n.TIpred=n.TIpred,
    mininterval=.001,
    individualRelativeTime=TRUE)
  rownames(datawide)<-1:nrow(datawide)
  return(datawide)
}

# Prints errors but does not stop, outputs NA instead
errorToNA<-function(call){
  tryCatch(call,error=function(e) {
    message(paste0("ERROR :",conditionMessage(e)))
    NA
  }
  )
}


# Removes average between person differences and group trends in time from ctsem wide data
ctDemean<-function(datawide,Tpoints,n.manifest,n.TDpred=0,persons=TRUE,trends=TRUE,grand=TRUE){
#   message('Caution: Predictor mean removal not yet implemented')
  namesvec<-colnames(datawide)
  
  if(persons==TRUE){
    message('Removing individual specific means from variables')
    for(i in 1:n.manifest){
      datawide<-
        matrix(apply(datawide,1,function(x){
          
          x[seq(i,Tpoints*n.manifest,n.manifest)] <- 
            x[seq(i,Tpoints*n.manifest,n.manifest)] - 
            mean(x[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T)
          
          return(x)}),byrow=T,nrow=nrow(datawide))
    }
    
    if(n.TDpred >0){
      
    tdpreds <- datawide[,(Tpoints*n.manifest+1) : (Tpoints*n.manifest + n.TDpred*(Tpoints-1)),drop=F] #extract tdpreds
    
    for(i in 1:n.TDpred){
      tdpreds<-matrix(apply(tdpreds,1,function(x){
        x[((i-1)*(Tpoints-1)+1) : ((i)*(Tpoints-1))] <- #for all variables of tdpred i
          x[((i-1)*(Tpoints-1)+1) : ((i)*(Tpoints-1))] - mean( #subtract the mean of tdpred i from them
            x[((i-1)*(Tpoints-1)+1) : ((i)*(Tpoints-1))], na.rm=T)
      }#end apply function
      ),nrow=nrow(tdpreds),byrow=T)
    } #end tdpred loop
    
    datawide[,(Tpoints*n.manifest+1) : (Tpoints*n.manifest + n.TDpred*(Tpoints-1))] <- tdpreds #reinsert tdpreds
    
  } #end tdpred if
  } #end persons if
  
  if(trends==TRUE){ #then remove time trends
    message('Removing group level trends over time in variables')
    
    datawide[,1:(Tpoints*n.manifest)]<-datawide[,1:(Tpoints*n.manifest)] -
      rep(apply(datawide[,1:(Tpoints*n.manifest),drop=F],2,mean,na.rm=T),each=nrow(datawide))
  }
  
  if(grand==TRUE){
    message('Removing grand mean over time in variables')
    for(i in 1:n.manifest){
      datawide[,seq(i,n.manifest*Tpoints,n.manifest)]<- datawide[,seq(i,n.manifest*Tpoints,n.manifest)] - 
        mean(datawide[,seq(i,n.manifest*Tpoints,n.manifest)],na.rm=T)
    }
  }
  
  colnames(datawide)<-namesvec
  return(datawide)
}

# Removes ceiling and floor effects from ctsem wide format data.
# Sets the observations for individual on a variable to NA if more than threshold percent are at max or min for that variable in the panel
ctDelimit<-function(datawide,Tpoints,n.manifest,threshold=.8){
  nameslist<-colnames(datawide)
  message(paste0('Setting any variables where more than ',threshold*100,'% of observations are at ceiling or floor to NA'))
  
  for(i in 1:n.manifest){
    maxv<-max(datawide[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T) #ceiling given data
    minv<-min(datawide[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T) #floor given data
    datawide<-
      matrix(apply(datawide,1,function(x){
        obs<-x[seq(i,Tpoints*n.manifest,n.manifest)] #observations of this manifest
        limit<- ceiling (length(!is.na(obs))*threshold) #if this many observations are at a limit
        if(length(obs) >1) { #if there is only one observation for an individual it may not be a ceiling or floor effect
          if(length(obs[obs==maxv]) >= limit | length(obs[obs==minv]) >= limit ){
            x[seq(i,Tpoints*n.manifest,n.manifest)] <- NA #set all observations on this variable to NA
          }
        }
        return(x)}),byrow=T,nrow=nrow(datawide))
  }
  colnames(datawide)<-nameslist
  
  return(datawide)
}

#Plot the distribution of variance on variables for individuals in a ctsem wide data set.
ctVarplot<-function(datawide,Tpoints,n.manifest){
  
  for(i in 1:n.manifest){
    indvar<-apply(datawide,1,function(x){ #variance on variable
      return(var(x[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T))
    }
    )
    #     browser()
    plot(density(indvar,na.rm=T),main=paste0('Variance on variable ',i))
    
    
    
    ranges<-apply(datawide,1,function(x){ #range on variable
      return(max(x[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T) - 
          min(x[seq(i,Tpoints*n.manifest,n.manifest)],na.rm=T))
    }
    )
    
    plot(density(ranges,na.rm=T),main=paste0('Range on variable ',i))
    
  }  
  return(datawide)
}

# Converts a string to mxAlgebra - From T Brick
stringToMxAlgebra <- function(algString, name=NA, dimnames=NA) { 
  eval(substitute(mxAlgebra(tExp, name=name, dimnames=dimnames), list(tExp = parse(text=algString)[[1]]))) 
}

# Plots panel ctsem wide data for either difference or squared difference
panelvis<-function(data,variables,Tpoints,type="difference",minrecords=0,spread=.1,sortby="variable"){
  
  data<-as.matrix(data[,-grep("I",colnames(data))])
  included<-apply(data,1,function(x){length(x[!is.na(x)])>minrecords})
  data<-data[included,]
  
  if(type=="difference" | type=="difsquared"){
    dif<-matrix(NA,nrow=nrow(data),ncol=variables*Tpoints)
    dif[,(variables+1):(Tpoints*variables)]<-
      data[,(variables+1):(Tpoints*variables)]-data[,1:((Tpoints-1)*variables)]
    
    
    if(type=="difference"){
      coldif<-as.matrix(round(
        pbeta(  (dif - min(dif,na.rm=T))  /  diff(range(dif,na.rm=T))  ,200,200),
        digits=3),ncol=1,nrow=nrow(data))
    }
    if(type=="difsquared"){
      dif<-scale(dif^2)
      coldif<-as.matrix(round(
        pbeta(  (dif - min(dif,na.rm=T))  /  diff(range(dif,na.rm=T))  ,2,200),
        digits=3),ncol=1,nrow=nrow(data))
    }
    
    
    #sort data
    data<-data[,cseq(1:variables,variables*Tpoints,variables)]
    coldif<-coldif[,cseq(1:variables,variables*Tpoints,variables)]
    
    #flatten coldif
    coldif<-matrix(coldif,ncol=1)
    
    coldif<-cbind(coldif,.2)
    coldif[which(is.na(coldif[,1])),2]<-0
    coldif[which(is.na(coldif[,1])),1]<-0.5
    coldif<-rgb(1-coldif[,1],0,coldif[,1],alpha=coldif[,2])
    
    plot(rep(1:(Tpoints*variables),each=nrow(data))+rnorm(nrow(data),0,spread),
      data+rnorm(nrow(data),0,spread),
      col=coldif,xlab="Time",ylab="Observed")
    legend("topright",legend=list("red = low, blue = high"))
    #   plot(rep(1:(Tpoints*variables),each=nrow(data))+rnorm(nrow(data),0,.1),
    #     data+rnorm(nrow(data),0,.1),
    #     col=coldif,xlab="Time",ylab="Observed")
    
  } #end difference types
  
  if(type=="indvar"){
    varind<-apply(data,1,var,na.rm=T)
    includedcount<-apply(data,1,function(x){length(x[!is.na(x)])})
    colind<-rgb(0,0,0,alpha=includedcount,maxColorValue=max(includedcount))
    plot(varind,col=colind)
  }
  if(type=="wavevar"){
    varwave<-apply(data,2,var,na.rm=T)
    includedcount<-apply(data,2,function(x){length(x[!is.na(x)])})
    colwave<-rgb(0,0,0,alpha=includedcount,maxColorValue=max(includedcount))
    plot(varwave,col=colwave)
  }
  
  if(type=="colbyvar"){
    varind<-apply(data,1,var,na.rm=T)
    varind[varind<0]<-0
    varind<-pbeta(varind,3,3)
    #     colvarind<-as.matrix(round(
    #       pbeta(  varind - min(varind,na.rm=T),200,200),
    #       digits=3),ncol=1,nrow=nrow(data))
    
    colvarind<-rep(varind,each=Tpoints)
    colvarind2<-rgb(max(varind,na.rm=T)-colvarind,0,colvarind,alpha=.2)
    plot(rep(1:(Tpoints*variables),each=nrow(data))+rnorm(nrow(data),0,spread),
      data+rnorm(nrow(data),0,spread),
      col=colvarind2,xlab="Time",ylab="Observed")
    legend("topright",legend=list("red = low, blue = high"))
  }
  
}










