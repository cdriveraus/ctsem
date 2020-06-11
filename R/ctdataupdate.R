ctdataupdate<-function(forcerecompile=FALSE){
  
  message(paste0('Updating from ',(getwd()),', continue T / F?'))
  continue <- readline()
  if(continue){
    set.seed(1)
    
  
  


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

  ctstantestdat<-ctGenerate(gm,n.subjects=n.subjects,burnin=3,logdtsd=.4,dtmean = .3)

  ctstantestdat[2,'Y1'] <- NA
  ctstantestdat[ctstantestdat[,'id']==2,'TI1'] <- NA
  ctstantestdat[2,'TD1'] <- NA

  save(ctstantestdat,file='.\\data\\ctstantestdat.rda')

  # checkm<-ctModel(type='stanct',
  #   n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
  #   MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
  #   MANIFESTMEANS=matrix(0,nrow=n.manifest),
  #   CINT=matrix(c('cint1','cint2'),ncol=1),
  #   n.manifest=n.manifest,LAMBDA=diag(2))
  # 
  # # checkm$pars$indvarying[-c(7,13)] <- FALSE
  # # checkm$pars$sdscale <- .1
  # # 
  # # checkm$pars[c(-1,-2, -21,-22) ,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
  # 
  # ctstantestfit<-ctStanFit(ctstantestdat,checkm,
  #   optimize = TRUE,optimcontrol=list(finishsamples=20),
  #   iter=300, warmup=260,thin=2,chains=2,verbose=0,
  #   # plot=TRUE,
  #   # forcerecompile=forcerecompile,
  #   save_warmup=TRUE,savescores=FALSE,
  #   control=list(max_treedepth=8,adapt_delta=.8))
  # ctstantestfit <- ctStanGenerateFromFit(ctstantestfit(),nsamples = 20,fullposterior = TRUE)
  # print( summary(ctstantestfit()))
  # 
  # save(ctstantestfit(),file='.\\data\\ctstantestfit.rda')
  
  paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
  tools::resaveRdaFiles(paths)
}
}
