ctStanModel<-function(ctmodelobj, type, indvarying='all'){
if(type=='stanct') continuoustime<-TRUE
if(type=='standt') continuoustime<-FALSE
  
  if(!is.null(ctmodelobj$timeVarying)) stop('Time varying parameters not allowed for ctsem Stan model at present! Correct ctModel spec')
  
  unlistCtModel<-function(ctmodelobj){
    out<-matrix(NA,nrow=0,ncol=5)
    out<-as.data.frame(out,stringsAsFactors =FALSE)
    objectlist<-c('T0MEANS','LAMBDA','DRIFT','DIFFUSION','MANIFESTVAR','MANIFESTMEANS', 'CINT', if(n.TDpred > 0) 'TDPREDEFFECT', 'T0VAR')
    for(obji in objectlist){
      if(!is.na(ctmodelobj[[obji]][1])){
      for(rowi in 1:nrow(ctmodelobj[[obji]])){
        for(coli in 1:ncol(ctmodelobj[[obji]])){
          out<-rbind(out,data.frame(obji,rowi,coli,
            ifelse(is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is character string
              ctmodelobj[[obji]][rowi,coli],
              NA),
            ifelse(!is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is numeric
              as.numeric(ctmodelobj[[obji]][rowi,coli]),
              NA),
            stringsAsFactors =FALSE
          ))
        }
      }
      }
    }
    colnames(out)<-c('matrix','row','col','param','value')
    return(out)
  }
  
  
  #read in ctmodel values
  n.latent<-ctmodelobj$n.latent
  n.manifest<-ctmodelobj$n.manifest
  Tpoints<-ctmodelobj$Tpoints
  n.TDpred<-ctmodelobj$n.TDpred
  n.TIpred<-ctmodelobj$n.TIpred
  
  manifestNames<-ctmodelobj$manifestNames
  latentNames<-ctmodelobj$latentNames
  TDpredNames<-ctmodelobj$TDpredNames
  TIpredNames<-ctmodelobj$TIpredNames
  
  ctspec<-unlistCtModel(ctmodelobj)
  
  freeparams<-is.na(ctspec[,'value'])
  
  ctspec$transform<-NA
  
  ######### STAN parameter transforms
  ctspec$transform[ctspec$matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT','CINT') & 
      freeparams] <- '(param) * 10'
  
  ctspec$transform[ctspec$matrix %in% c('LAMBDA') &  freeparams] <- '(param+.5) * 10'
  
  ctspec$transform[ctspec$matrix %in% c('DIFFUSION','MANIFESTVAR', 'T0VAR') & 
      freeparams & ctspec$row != ctspec$col] <- 'inv_logit(param)*2-1'
  
  ctspec$transform[ctspec$matrix %in% c('DIFFUSION','MANIFESTVAR', 'T0VAR') & 
      freeparams & ctspec$row == ctspec$col] <- 'exp(param*2) +.00001' #'1/(.1+exp(param*1.8))*10+.001'
  
  if(continuoustime==TRUE){
    ctspec$transform[ctspec$matrix %in% c('DRIFT') & 
        freeparams & ctspec$row == ctspec$col] <- '-log(exp(-param*1.5)+1)-.00001' #'log(1/(1+(exp(param*-1.5))))'
    
    ctspec$transform[ctspec$matrix %in% c('DRIFT') & freeparams & 
        ctspec$row != ctspec$col] <- '(param)*.5'
  }
  if(continuoustime==FALSE){
    ctspec$transform[ctspec$matrix %in% c('DRIFT') & freeparams & 
        ctspec$row == ctspec$col] <- '1/(1+exp((param)*-1.5))'
    
    ctspec$transform[ctspec$matrix %in% c('DRIFT') & freeparams & 
        ctspec$row != ctspec$col] <- '(param)'
  }
  
  nparams<-sum(freeparams)
  
  if(all(indvarying=='all'))  indvarying<-rep(TRUE,nparams)
  if(length(indvarying) != sum(freeparams)) stop('indvarying must be ', nparams,' long!')
  nindvarying <- sum(indvarying)
  
  
  ctspec$indvarying<-NA
  ctspec$indvarying[is.na(ctspec$value)]<-indvarying
  ctspec$indvarying[!is.na(ctspec$value)]<-FALSE
  
  ctspec$sdscale<-1
  
  
 
  if(n.TIpred > 0) {
    tipredspec<-matrix(TRUE,ncol=n.TIpred,nrow=1)
    colnames(tipredspec)<-paste0(TIpredNames,'_effect')
    ctspec<-cbind(ctspec,tipredspec,stringsAsFactors=FALSE)
    ctspec[!ctspec$indvarying,paste0(TIpredNames,'_effect')]<-FALSE
    for(predi in TIpredNames){
      class(ctspec[,paste0(predi,'_effect')])<-'logical'
    }
    if(sum(unlist(ctspec[,paste0(TIpredNames,'_effect')]))==0) stop('TI predictors included but no effects specified!')
  }
  out<-list(parameters=ctspec,n.latent=n.latent,n.manifest=n.manifest,n.TIpred=n.TIpred,n.TDpred=n.TDpred,
    latentNames=latentNames,manifestNames=manifestNames,TIpredNames=TIpredNames,TDpredNames=TDpredNames,subjectIDname='id',
    timeName='time',
    continuoustime=continuoustime)
  class(out)<-'ctsemStanModel'
  
  if(n.TIpred > 0) out$tipredeffectprior <- 'normal(0,1)'
  # out$hypersdpriorscale <- 1
  # out$hypersdtransform <- 'hypersd .* sdscale'
  out$stationarymeanprior <- NA
  out$stationaryvarprior <- NA
  
 return(out)
}
