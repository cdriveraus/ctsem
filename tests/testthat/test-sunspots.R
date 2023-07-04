library(ctsem)
library(testthat)
set.seed(1)

context("sunspots") 

test_that("sunspots", { 
sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
    # n.TDpred = 1,
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 0), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 'a21|-log1p(exp(param))', 1, 'a22'), nrow=2, ncol=2),
   TDPREDEFFECT=matrix(c('tdeffect',0),2),
   MANIFESTMEANS=0, CINT=c(0,'cint'),
   # MANIFESTVAR=diag(0,1),
   T0VAR=matrix(c(1e-3,0,0,1e-3), nrow=2, ncol=2), #Because single subject
   DIFFUSION=c(0, 'ma1|param-1', 
     0,'diff'))

 #td preds for testing only -- no real effect
 TD1 <- 0
 datalong <- cbind(datalong,TD1)
 datalong[seq(10,150,10),'TD1'] = 1
 
 ssmodel$covmattransform='cholesky'

ssfit1 <- ctStanFit(datalong, ssmodel,cores=1,verbose=0)
# ssfit <- ctStanFit(datalong, ssmodel,cores=1,verbose=0,optimcontrol = list(hessianType='stochastic', stochasticHessianEpsilon=1e-1))
ssfit2 <- ctStanFit(datalong, ssmodel,cores=2,verbose=0)
ssfit3 <- ctStanFit(datalong, ssmodel,cores=1,nlcontrol=list(maxtimestep=.3))
ssfit4 <- ctStanFit(datalong, ssmodel,chains=2,cores=2,iter=300,optimize=F,verbose=0,priors = F,inits='optimize',intoverpop = T)
# ssfit5 <- ctStanFit(datalong, ssmodel,chains=3,iter=300,intoverstates = FALSE,optimize=F,verbose=0)

for(i in 2:4){
testthat::expect_equivalent(get(paste0('ssfit',i))$stanfit$transformedparsfull$ll,
   get(paste0('ssfit',i-1))$stanfit$transformedparsfull$ll,tol=1e-2)
}

for(i in 2:4){
   testthat::expect_equivalent(
      ctStanContinuousPars(get(paste0('ssfit',i)))$DRIFT,
      ctStanContinuousPars(get(paste0('ssfit',i-1)))$DRIFT,tol=1e-1)
}

})
