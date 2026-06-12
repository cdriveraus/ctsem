ctdataupdate<-function(forcerecompile=FALSE){
  
  message(paste0('Updating from ',(getwd()),', continue T / F?'))
  continue <- readline()
  if(continue){
    set.seed(1)

  Tpoints=10
  n.manifest=2
  n.TDpred=1
  n.TIpred=3
  n.latent=2
  n.subjects=30
  tipredEffect <- matrix(c(.5,0,0,-.7,0,2),nrow=2)
  tipredVar <- matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3)
  gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
  n.TDpred=n.TDpred,
    n.TIpred=0,
    n.manifest=n.manifest,
    MANIFESTVAR=diag(0.5,2),
    TDPREDMEANS=matrix(round(exp(rnorm(n.TDpred*(Tpoints),-1.9,1)),0),
     nrow=n.TDpred*(Tpoints)),
     TDPREDEFFECT = matrix(c(1,-1),ncol=1),
    LAMBDA=diag(1,2),
    DRIFT=matrix(c(-.3,.2,0,-.2),nrow=2),
    DIFFUSION=matrix(c(2,1,0,2),2),
    CINT=matrix(c(0,0),nrow=2),
    T0MEANS=matrix(10,ncol=1,nrow=2),
    T0VAR=diag(1,2))

  for(i in seq_len(n.subjects)){
    gm_i <- gm
    tipreds_i <- tipredVar %*% rnorm(n.TIpred)
    gm_i$CINT <- gm$CINT + tipredEffect %*% tipreds_i
    dat_i <- ctGenerate(gm_i,n.subjects=1,burnin=3,logdtsd=.4,dtmean = .3)
    dat_i[,'id'] <- i
    dat_i <- cbind(dat_i, matrix(tipreds_i, nrow=nrow(dat_i), ncol=n.TIpred, byrow=TRUE))
    if(i == 1) ctstantestdat <- dat_i else ctstantestdat <- rbind(ctstantestdat, dat_i)
  }
  colnames(ctstantestdat)[(ncol(ctstantestdat)-n.TIpred+1):ncol(ctstantestdat)] <- paste0('TI',1:n.TIpred)

  ctstantestdat[2,'Y1'] <- NA
  ctstantestdat[ctstantestdat[,'id']==2,'TI1'] <- NA
  ctstantestdat[2,'TD1'] <- NA

  save(ctstantestdat,file='.\\data\\ctstantestdat.rda')

  ## now in zzz.R
  checkm<-ctModel(
    type='ct',
    n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=0,
    DRIFT=c('dr1','dr12','dr21||||TI1','dr22'),
    DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
    CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
    LAMBDA=diag(2),tipredDefault=FALSE)

  ctstantestfit<-ctFit(ctstantestdat,checkm,cores=1,inits=0,
    optimize = TRUE,optimcontrol=list(finishsamples=20,stochastic=T,tol=1e-5),priors=TRUE)

  ctstantestfit <- ctGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE,cores=1)

  save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
  
  paths <- sort(Sys.glob(c("data/ctstantestfit.rda","data/ctstantestdat.rda"))) #
  tools::resaveRdaFiles(paths)
}
}
