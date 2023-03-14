ctsmupdate<-function(usecurrentwd=FALSE,scat=FALSE){
  
  sunspots<-datasets::sunspot.year
  sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
  id <- 1
  time <- 1749:1924
  datalong <- cbind(id, time, sunspots)
  
  #setup model
  model <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
    manifestNames='sunspots', 
    latentNames=c('ss_level', 'ss_velocity'),
    LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
    DRIFT=matrix(c(-.0001, 'a21', 1, 'a22'), nrow=2, ncol=2),
    MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
    CINT=matrix(c(0, 0), nrow=2, ncol=1),
    MANIFESTVAR=diag(.001,1),
    T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
    DIFFUSION=matrix(c(.0001, 0, 0, "diffusion"), ncol=2, nrow=2))
  
  #fit
  
  sm <- ctStanFit(datalong, model,fit=FALSE,gendata=FALSE,forcerecompile=TRUE)$stanmodeltext
  
  if(scat) scat(sm)
  
  stanc(model_code = sm,verbose = TRUE)
  
  smgen <- ctStanFit(datalong, model,fit=FALSE,gendata=TRUE,forcerecompile=TRUE)$stanmodeltext
  stanc(model_code = smgen,verbose = TRUE)
  
  
  
  # model$w32 <- TRUE
  # smgen32 <- ctStanFit(datalong, model,fit=FALSE,gendata=TRUE)$stanmodeltext
  # stanc(model_code = smgen32,verbose = TRUE)
  # 
  # sm32 <- ctStanFit(datalong, model,fit=FALSE,gendata=FALSE)$stanmodeltext
  # stanc(model_code = sm32,verbose = TRUE)
  
  message(paste0('Update files? T / F?'))
  continue <- readline()
  if(continue){
    pkgdir=ifelse(usecurrentwd, paste0(getwd()),'~/../sync/CT-SEM/ctsem')
    pathbase <- ifelse(usecurrentwd, paste0(getwd(),'/inst/stan'),'~/../sync/CT-SEM/ctsem/inst/stan')
    
    file.rename(file.path(pathbase,'ctsm.stan'), file.path(pathbase,'ctsm.bak'))
    file.rename(file.path(pathbase,'ctsmgen.stan'), file.path(pathbase,'ctsmgen.bak'))
    sink(file=file.path(pathbase,'ctsm.stan'))
    cat(sm)
    sink()
    sink(file=file.path(pathbase,'ctsmgen.stan'))
    cat(smgen)
    sink()
    rstantools::rstan_config(pkgdir = pkgdir)
  }
}
