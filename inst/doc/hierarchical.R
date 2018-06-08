## ----setup, include = FALSE, cache = FALSE, echo = FALSE------------------------------------------
library('ctsem')
#library(knitr)
set.seed(22)
knit_hooks$set(crop = hook_pdfcrop)
opts_chunk$set(warning = FALSE, fig.align = 'center', width.cutoff = 100, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, comment = NA, tidy = FALSE, autodep=TRUE, out.truncate = 100, size='small', crop=TRUE, fig.pos="htbp")
options(width = 100, scipen = 12, digits = 3)

set.seed(1)
# Tpoints=50
# n.manifest=2
# n.TDpred=1
# n.TIpred=3
# n.latent=2
# n.subjects=20
# gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#   MANIFESTVAR=diag(0.5,2),
#   TIPREDEFFECT=matrix(c(.5,0,0,-.5,0,0),nrow=2),
#   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
#   TDPREDEFFECT=matrix(c(.1,-.2),nrow=2),
#   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints-1),ncol=n.TDpred*(Tpoints-1)),
#   TDPREDMEANS=matrix(rnorm(n.TDpred*(Tpoints-1),0,1),nrow=n.TDpred*(Tpoints-1)),
#   LAMBDA=diag(1,2), 
#   # DRIFT=matrix(c(-.6+rnorm(1,0,.15),-.2+rnorm(1,0,.1),.12+rnorm(1,0,.1),-.3+rnorm(1,0,.05)),nrow=2),
#   DRIFT=matrix(c(-.3,.2,-.1,-.2),nrow=2),
#   TRAITVAR=t(chol(matrix(c(4,3,3,4),nrow=2))),
#   # T0TRAITEFFECT=diag(3,n.latent),
#   DIFFUSION=matrix(c(.3,.1,0,.2),2),CINT=matrix(c(0,0),nrow=2),T0MEANS=matrix(0,ncol=1,nrow=2),
#   T0VAR=diag(100,2))
# 
# cd<-ctGenerate(gm,n.subjects=n.subjects,burnin=300, dT=1,asymptotes=F,simulTDpredeffect = T)
# model<-ctModel(type='stanct',n.latent=n.latent,n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,LAMBDA=diag(n.latent))
# long<-ctWideToLong(cd,Tpoints,n.manifest=model$n.manifest,manifestNames = model$manifestNames, 
#   n.TDpred=n.TDpred,n.TIpred=n.TIpred,TDpredNames = model$TDpredNames,TIpredNames = model$TIpredNames)
# long<-ctDeintervalise(long)
# long[is.na(long)]<-0
# ctstantestdat <- long


## ----install,eval=FALSE---------------------------------------------------------------------------
#  install.packages("ctsem")
#  library("ctsem")

## ----data,echo=FALSE,size='footnotesize'----------------------------------------------------------
ctstantestdat[c(3:5,17:22),]

## ----scaling,echo=TRUE----------------------------------------------------------------------------
ctstantestdat[,c('Y1','Y2','TI1','TI2','TI3')] <- 
 scale(ctstantestdat[,c('Y1','Y2','TI1','TI2','TI3')])

## ----model----------------------------------------------------------------------------------------
model<-ctModel(type='stanct',
  n.latent=2, latentNames=c('eta1','eta2'),
  n.manifest=2, manifestNames=c('Y1','Y2'),
  n.TDpred=1, TDpredNames='TD1', 
  n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
  LAMBDA=diag(2))

## ----modelpars,size='footnotesize'----------------------------------------------------------------
head(model$pars,7)

## ----transform, fig.width=8, fig.height=5,fig.cap="Prior distribution density plots."-------------
par(mfrow=c(1,2))
plot(model,rows=7,rawpopsd=1)
model$pars$transform[7]<- 0
model$pars$meanscale[7] <- 2
model$pars$multiplier[7] <- 1
model$pars$offset[7] <- -1
plot(model, rows=7,rawpopsd=1)

## ----transform2, fig.width=8, fig.height=5,fig.cap="Prior distribution density plots of auto-effects, with default (left) and adjusted (right) scale parameter for population standard deviation."----
par(mfrow=c(1,2))
plot(model,rows=7)
model$pars$sdscale<- 0.1
plot(model, rows=7)

## ----restrictbetween------------------------------------------------------------------------------
model$pars$indvarying[!(model$pars$matrix %in% c('DRIFT','MANIFESTMEANS'))] <- FALSE

## ----restricttipred-------------------------------------------------------------------------------
model$pars[,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
model$pars[model$pars$matrix == 'DRIFT', 
  c('TI1_effect','TI2_effect','TI3_effect')] <- TRUE

## ----fitting,include=FALSE,cache=TRUE-----------------------------------------------------------------------
suppressWarnings(fit<-ctStanFit(datalong = ctstantestdat, ctstanmodel = model, iter=200,
  # ukfpop=TRUE,optimize=TRUE,cores=2,isloops=10,verbose=1,
  chains=2, plot=FALSE))
# fit <- ctstantestfit

## ----fittingshow,include=TRUE,eval=FALSE--------------------------------------------------------------------
#  fit<-ctStanFit(datalong = ctstantestdat, ctstanmodel = model, iter=300,
#    control=list(max_treedepth=6), chains=2, plot=FALSE)
#  # fit <- ctstantestfit

## ----output,eval=FALSE--------------------------------------------------------------------------------------
#  summary(fit,timeinterval = 1)

## ----ctStanContinuousPars,eval=FALSE------------------------------------------------------------------------
#  ctStanContinuousPars(fit,subjects = 3, calcfunc = quantile, calcfuncargs = list(probs=.975))

## ----plots1,echo=TRUE,fig.width=8, fig.height=5,fig.cap='Discrete-time cross-effect dynamics of the estimated system for a range of time intervals, with 95\\% credible intervals.'----
ctStanDiscretePars(fit, plot=TRUE, indices = 'CR', subjects = 'all')

## ----outputposterior, fig.width=8, fig.height=6, fig.cap='Prior and posterior densities relevant to the second process auto effect.'----
ctStanPlotPost(obj = fit, rows=3) 

## ----kalmanplot,echo=TRUE, fig.width=8, fig.height=5, fig.cap='Predicted and smoothed estimates for one subject with two processes. Uncertainty shown is a 95\\% credible interval comprising both process and measurement error.'----
par(mfrow=c(1,2))
 
ctKalman(fit, subjects=2, timerange=c(0,30), kalmanvec=c('y', 'yprior'), timestep=.01, 
plot=TRUE, plotcontrol=list(xaxs='i', main = 'Predicted'))
 
ctKalman(fit, subjects=2, timerange=c(0,30), kalmanvec=c('y', 'ysmooth'), timestep=.01, 
plot=TRUE, plotcontrol=list(xaxs='i',main = 'Smoothed'))

## ----tipredeffects,echo=TRUE, fig.width=7, fig.height=5, fig.cap='Expectations for individuals parameter values change depending on their score on time independent predictors.'----
ctStanTIpredeffects(fit, plot = TRUE, whichpars=c('dtDRIFT','MANIFESTVAR[2,2]'), 
 timeinterval = .5, whichTIpreds = 3, includeMeanUncertainty = FALSE, nsubjects=10,
  nsamples = 50)

## ----sunspots,include=FALSE, cache=TRUE---------------------------------------------------------------------
#get data
 sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
   CINT=matrix(c(0, 0), nrow=2, ncol=1),
   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
   DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
 
 ssmodel$pars$indvarying<-FALSE #Because single subject
 ssmodel$pars$offset[14]<- 44 #Because not mean centered
 ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality 

#fit
ssfit <- ctStanFit(datalong, ssmodel, iter=300, chains=2)

#output
summary(ssfit)$popmeans

## ----sunspotsshow,include=TRUE,eval=FALSE-------------------------------------------------------------------
#  #get data
#   sunspots<-sunspot.year
#   sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#   id <- 1
#   time <- 1749:1924
#  datalong <- cbind(id, time, sunspots)
#  
#  #setup model
#   ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1,
#    manifestNames='sunspots',
#    latentNames=c('ss_level', 'ss_velocity'),
#     LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#     DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
#     MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#     CINT=matrix(c(0, 0), nrow=2, ncol=1),
#     T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
#     DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#  
#   ssmodel$pars$indvarying<-FALSE #Because single subject
#   ssmodel$pars$offset[14]<- 44 #Because not mean centered
#   ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality
#  
#  #fit
#  ssfit <- ctStanFit(datalong, ssmodel, iter=300, chains=2)
#  
#  #output
#  summary(ssfit)$popmeans

## ----transforms,  fig.width=8, fig.height=6, fig.cap='Depiction of the prior distributions and sampling process through which individual specific parameters are determined.'----
#set plotting parameters
par(mfrow=c(2,2), lwd=3, yaxs='i', mgp=c(1.8,.5,0), 
  mar=c(3,3,3,1)+.1)
bw=.03 

n <- 999999 #number of samples to draw to from prior for plotting purposes
nsubjects <- 4 #number of subjects

#parameter specific transform
tform <- function(x) -log(exp(-1.5 * x) + 1) #default drift auto effect transform

#raw pop sd transform
sdscale <- 1 #default
rawsdtform <- function(x) exp(x * 2 -2) * sdscale #default

#sd approximation function
sdapprox <- function(means,sds,tform) {
  for(i in 1:length(means)){
    sds[i] <- ((tform(means[i]+sds[i]*3) - tform(means[i]-sds[i]*3))/6 + 
      (tform(means[i]+sds[i]) - tform(means[i]-sds[i]))/2) /2
  }
  return(sds)
}

#raw population mean parameters
rawpopmeans_prior <- rnorm(n, 0, 1) #prior distribution for rawpopmeans
rawpopmeans_sample <- -.3 #hypothetical sample
sdscale <- 1 #default

#population mean parameters after parameter specific transform
popmeans_prior <- tform(rawpopmeans_prior)
popmeans_sample <- tform(rawpopmeans_sample)

#plot pop means
plot(density(rawpopmeans_prior), ylim=c(0,1), xlim=c(-5,2),
  xlab='Parameter value', main='Population means')
points(density(popmeans_prior, bw=bw),col=2,type='l')
segments(y0=0,y1=.5,x0=c(rawpopmeans_sample,popmeans_sample),lty=3,col=1:2) 
legend('topleft',c('Raw pop. mean prior', 'Pop. mean prior', 
  'Raw pop. mean sample', 'Pop. mean sample'),lty=c(1,1,3,3), col=1:2, bty='n')

#population standard deviation parameters
rawpopsd_prior <- rawsdtform(rnorm(n, 0, 1)) #raw population sd prior

popsd_prior <- sdapprox(rawpopmeans_prior,rawpopsd_prior,tform)

#sample population standard deviation posterior
rawpopsd_sample <- rawsdtform(.9) #hypothetical sample
popsd_sample <- sdapprox(means=rawpopmeans_sample,  #transform sample to actual pop sd
  sds=rawpopsd_sample,tform=tform)

#plot pop sd
plot(density(rawpopsd_prior,from=-.2,to=10,na.rm=TRUE, bw=bw), xlab='Parameter value', 
  xlim=c(-.1,3), ylim=c(0,2), main='Population sd')
points(density(popsd_prior,from=-.2,to=10,na.rm=TRUE, bw=bw),type='l', col=2)
segments(y0=0,y1=1,x0=c(rawpopsd_sample, popsd_sample), col=1:2,lty=3)
legend('topright',c('Raw pop. sd prior','Pop. sd prior', 
  'Raw pop. sd sample','Pop. sd sample'), col=1:2, lty=c(1,1,3,3),bty='n')

#individual level parameters

#marginal individual level parameters (given all possible values for mean and sd)
rawindparams_margprior <- rawpopmeans_prior + rawpopsd_prior * rnorm(n, 0, 1) 
indparams_margprior <- tform(rawindparams_margprior)

plot(density(rawindparams_margprior,from=-10,to=10,bw=bw), xlab='Parameter value', 
  xlim=c(-5,2), ylim=c(0,1), main='Marginal dist. individual parameters')
points(density(indparams_margprior,from=-10,to=.2,bw=bw),type='l',col=2)
legend('topleft',c('Raw individual parameters prior','Individual parameters prior'), 
  col=1:2,lty=1,bty='n')


#conditional individual level parameters (given sampled values for mean and sd)
rawindparams_condprior<- rawpopmeans_sample + rawpopsd_sample * rnorm(n,0,1) 
rawindparams_condsample<- rawpopmeans_sample + rawpopsd_sample * rnorm(nsubjects,0,1) 
indparams_condprior<- tform(rawindparams_condprior)
indparams_condsample<- tform(rawindparams_condsample)


plot(density(rawindparams_condprior), xlab='Parameter value', xlim=c(-5,2),
  ylim=c(0,1), main='Conditional dist. individual parameters')
points(density(indparams_condprior),type='l',col=2)
segments(y0=0,y1=.5,x0=c(rawindparams_condsample, indparams_condsample),
  col=rep(1:2,each=nsubjects),lty=3, lwd=2)
legend('topleft',c('Raw ind. pars. prior','Ind. pars. prior', 
  'Raw ind. pars. samples','Ind. pars. samples'), col=1:2, lty=c(1,1,3,3),bty='n')


