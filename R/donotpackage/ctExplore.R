ctExplore<-function(datalong,id, time, manifestNames, TIpredNames, maxProcesses){
  

dlong<-datalong[order(datalong[,id]),]
sortedid<-dlong[,id,drop=F]
colnames(sortedid)<-'id'
uniqueid<-matrix(unique(sortedid),ncol=1)
colnames(uniqueid)<-'id'
uniqueid<-cbind(1:length(uniqueid),uniqueid)
#set simple id's
newid<-merge(uniqueid,sortedid)[,2]
  
  browser()
  n.subjects<-length(unique(sortedid))
# max(unlist(lapply(unique(id),function(x) length((id==x)==T))))

  
  
  #determine max Tpoints
Tpoints<-0
newMax<-0
oldMax<-0
for(i in 2:length(newid)){
if(newid[i] != newid[i-1]) {
  if(newMax > oldMax) oldMax <- newMax
  newMax <- 1
}
  if(newid[i] == newid[i-1]) newMax <- newMax + 1
}
if(newMax > oldMax) oldMax <- newMax
Tpoints<-oldMax

  
  #calculate means of manifests
  manifestmeans<- cbind(unique(dlong[,id]))
  
  for(manifest in manifestNames){
    print(manifest)
    
    obs<-rep(list(list()),n.subjects)
    for(i in 1:nrow(datalong)){
      obs[[newid[i]]]<-c(unlist(obs[[newid[i]]]),c(dlong[i,manifest]))
    }
    
    manifestmeans<-cbind(manifestmeans,unlist(lapply(obs,mean,na.rm=T)))
  
  }
}