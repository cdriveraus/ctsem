#'ctStanFit example fit
#'
#'@name ctstantestfit
#'@export
#'@examples
#'\donttest{
#'testfit <- ctstantestfitgen()
#'}
#'

ctstantestfitgen<-function(){
  checkm<-ctModel(
    type='stanct',
    n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=0,
    DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
    CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
    LAMBDA=diag(2),tipredDefault=FALSE)  
  
  
  ctstantestfit<-ctStanFit(ctstantestdat,checkm,cores=1,
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
  
  return(ctstantestfit)
}
