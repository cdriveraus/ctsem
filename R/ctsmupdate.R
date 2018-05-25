ctsmupdate<-function(){
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
sm <- ctStanFit(datalong, model,fit=FALSE)$stanmodeltext
stanc(model_code = sm,verbose = TRUE)
message(paste0('Updating from ',(getwd()),', continue T / F?'))
continue <- readline()
if(continue){
  file.rename('./src/stan_files/ctsm.stan', './src/stan_files/ctsm.bak')
sink(file='./src/stan_files/ctsm.stan')
cat(sm)
sink()
# 
# message('All ok? finish this...')
# compile <- readline()
# if(compile) eval(parse(text=paste0('eval(devtools::install(local=FALSE),envir = globalenv())')))
}
}
