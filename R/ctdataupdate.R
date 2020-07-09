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

  checkm<-ctModel(
    type='stanct',
    n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=0,
    DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
    CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
    LAMBDA=diag(2),tipredDefault=FALSE)  
  
  ctstantestfit<-ctStanFit(ctsem::ctstantestdat,checkm,cores=1,
    inits = c(0.748310681869536,0.945659953796114,0.0964592332562144,
      0.029153487981562,0.651471066485501,0.0314778013950629,
      0.217818608752396,1.10441297459423,-0.801320300354595,
      0.647010811111734,-0.7344068376597,-1.04150782976995,
      0.0558819480347101,-0.108435212373754,-0.225029736388403,
      -0.203457959897841,-0.736264486213394,-0.687369939087293,
      0.576641002392084,0.248625561427667,-0.0683189297539777,
      -0.230342395895042,0.205299380670756,-0.34522281922735,
      0.0829819407118698,0.0137367678089216,-0.0611697475527028),
    optimize = TRUE,optimcontrol=list(finishsamples=20),nopriors=FALSE)
  
  ctstantestfit <- ctStanGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE)
  
  save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
  
  paths <- sort(Sys.glob(c("data/ctstantestfit.rda","data/ctstantestdat.rda")))
  tools::resaveRdaFiles(paths)
}
}
