ctdataupdate<-function(forcerecompile=FALSE){
  
  message(paste0('Updating from ',(getwd()),', continue T / F?'))
  continue <- readline()
  if(continue){
    set.seed(1)
    
  #two process, one time dependent predictor example
  Tpoints=20
  manifestNames<-c('LeisureTime','Happiness')
  TDpredNames<-'MoneyInt'
  testm<-ctModel(Tpoints=Tpoints,n.latent=3,n.TDpred=1,n.TIpred=0,n.manifest=2,
    LAMBDA=cbind(diag(1,2),0),
    MANIFESTVAR=diag(.1,2),
    DRIFT=matrix(c(-.3,.12,0,  -.02,-.3,0, 1,-.3,-.0001  ),nrow=3,ncol=3),
    TRAITVAR=t(chol(matrix(c(.2,-.1,0,  -.1,.21,0,  0,0,0.00001),ncol=3,nrow=3))),
    DIFFUSION=t(chol(diag(c(1.2,.6,0.0001),3))),
    CINT=matrix(c(1,.3,0),nrow=3),
    T0MEANS=matrix(0,ncol=1,nrow=3),
    T0VAR=diag(c(1,1,0),3),
    TDPREDEFFECT=matrix(c(.6,.4,1),nrow=3),
    TDPREDVAR=diag(c(rep(0,Tpoints)),Tpoints),
    TDPREDMEANS=matrix(c(0,0,0,0,0,1,rep(0,Tpoints-6)),ncol=1,nrow=(Tpoints)))
  testd<-ctGenerate(testm,n.subjects=10,burnin=10) #generate data

  ctIndplot(testd,Tpoints=Tpoints,n.manifest=2,n.subjects=10,colourby="variable")

  timestokeep=c(0,1,4,5,7,8,16,19)
  deltaT<-timestokeep[-1] - timestokeep[-8]
  testd<-testd[,c(paste0('Y',1:2,'_T',rep(timestokeep,each=2)),paste0('TD1_T',timestokeep))]
  testd<-cbind(testd,matrix(deltaT,nrow=nrow(testd),ncol=length(deltaT),byrow=TRUE))

  colnames(testd)<-ctWideNames(n.manifest=2,Tpoints=8,n.TDpred=1,
  manifestNames=manifestNames,TDpredNames=TDpredNames)
  ctExample2<-testd
  save(ctExample2,file=".\\data\\ctExample2.rda")
  
  
  
  Tpoints=20
  subjects=20
  full<-c()
  for(i in 1:20){
    LAMBDA<-matrix(c(1,.7,ifelse(i >(subjects/2),.2,1.4)))
    # print(LAMBDA)
    testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=0,n.TIpred=0,n.manifest=3,
      MANIFESTVAR=diag(.2,3),
      # TRAITVAR=diag(.2,1),
      LAMBDA=LAMBDA,
      DRIFT=matrix(c(-.1),nrow=1,ncol=1),
      DIFFUSION=diag(c(.12),1),
      MANIFESTMEANS=matrix(c(0,.42,1.3),ncol=1),
      CINT=matrix(c(.2),ncol=1))

    testd<-ctGenerate(testm,n.subjects=1,burnin=300)
    full<-rbind(full,testd)
  }

  ctExample4<-full
  save(ctExample4,file=".\\data\\ctExample4.rda") #save wide format example
  
  
  
  
  
  
  Tpoints=30
  testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,
    LAMBDA=matrix(1,ncol=1,nrow=3),
    DRIFT=diag(-.3,1),
    DIFFUSION=diag(.1,1),
    CINT=diag(2,1),
    MANIFESTVAR=diag(1,3),
    TDPREDEFFECT=diag(.2,1),
    TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
    TDPREDVAR=diag(1,1*(Tpoints)),
    TIPREDVAR=diag(1,2)
  )
  longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:34),],2)
  longexample[2,c(2,7)]<-NA
  longexample[4,c(3)]<-NA
  datastructure <- ctLongToWide(datalong = longexample,id='id',time='time',
    manifestNames = testm$manifestNames,TDpredNames = testm$TDpredNames,
    TIpredNames=testm$TIpredNames)
  datastructure<-ctIntervalise(datawide = datastructure,
    Tpoints = 3,n.manifest = testm$n.manifest,n.TDpred = testm$n.TDpred,
    n.TIpred=testm$n.TIpred)
  save(datastructure,file='.\\data\\datastructure.rda')

  
  
  
  
  

  #long example (using datastructure base)
  Tpoints=30
  testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,
    LAMBDA=matrix(1,ncol=1,nrow=3),
    DRIFT=diag(-.3,1),
    DIFFUSION=diag(.1,1),
    CINT=diag(2,1),
    MANIFESTVAR=diag(1,3),
    TDPREDEFFECT=diag(.2,1),
    TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
    TDPREDVAR=diag(1,1*(Tpoints)),
    TIPREDVAR=diag(1,2)
  )
  longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:35),],2)
  longexample[2,c(2,7)]<-NA
  longexample[4,c(3)]<-NA
  longexample
  save(longexample,file='.\\data\\longexample.rda')




  set.seed(1)

  Tpoints=10
  n.manifest=2
  n.TDpred=1
  n.TIpred=3
  n.latent=2
  n.subjects=30
  gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
  n.TDpred=n.TDpred,
    n.TIpred=n.TIpred,
    n.manifest=n.manifest,
    MANIFESTVAR=diag(0.5,2),
    TIPREDEFFECT=matrix(c(.5,0,0,-.7,0,2),nrow=2),
    TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
    TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints),ncol=n.TDpred*(Tpoints)),
    TDPREDMEANS=matrix(round(exp(rnorm(n.TDpred*(Tpoints),-1.9,1)),0),
     nrow=n.TDpred*(Tpoints)),
     TDPREDEFFECT = matrix(c(1,-1),ncol=1),
    LAMBDA=diag(1,2),
    DRIFT=matrix(c(-.3,.2,0,-.2),nrow=2),
    DIFFUSION=matrix(c(2,1,0,2),2),
    CINT=matrix(c(0,0),nrow=2),
    T0MEANS=matrix(10,ncol=1,nrow=2),
    T0VAR=diag(1,2))

  ctstantestdat<-ctGenerate(gm,n.subjects=n.subjects,burnin=3,
  wide=FALSE,logdtsd=.4,dtmean = .3)

  ctstantestdat[2,'Y1'] <- NA
  ctstantestdat[ctstantestdat[,'id']==2,'TI1'] <- NA
  ctstantestdat[2,'TD1'] <- NA

  save(ctstantestdat,file='.\\data\\ctstantestdat.rda')

  checkm<-ctModel(type='stanct',
    n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=matrix(0,nrow=n.manifest),
    CINT=matrix(c('cint1','cint2'),ncol=1),
    n.manifest=n.manifest,LAMBDA=diag(2))
  
  # checkm$pars$indvarying[-c(7,13)] <- FALSE
  # checkm$pars$sdscale <- .1
  # 
  # checkm$pars[c(-1,-2, -21,-22) ,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
  
  ctstantestfit<-ctStanFit(ctstantestdat,checkm,iter=300, warmup=260,thin=2,chains=2,
    # plot=TRUE,
    forcerecompile=forcerecompile,save_warmup=FALSE,savescores=FALSE,
    control=list(max_treedepth=8,adapt_delta=.8))
  ctstantestfit <- ctStanGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE)
  print( summary(ctstantestfit))

  save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
  
  paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
  tools::resaveRdaFiles(paths)
}
}
