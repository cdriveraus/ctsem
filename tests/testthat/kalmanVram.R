require(ctsem)
require(testthat)

context("kalmanVram")

test_that("time calc", {
  set.seed(4)
nsubjects=20
  #2 latent 2 td preds
  
gm=ctModel(Tpoints=10,n.latent=2,n.manifest=3,LAMBDA=matrix(c(1,0,0,0,1,.5),ncol=2),
  DRIFT=matrix(c(-.3,0,0,-.05),2,2),DIFFUSION=diag(2,2),
  # TRAITVAR=diag(2),
  CINT=matrix(c(5,3),ncol=1),
  TDPREDEFFECT=matrix(c(-4,-3,-7,-10),2,2),
  MANIFESTVAR=diag(.3,3),
  n.TDpred=2,
  TDPREDMEANS=matrix(rep(c(1,rep(0,4),1,rep(0,4)),2),ncol=1))

gd=ctGenerate(gm,nsubjects,burnin=5)
ctIndplot(datawide = gd,n.manifest = 3,Tpoints = 10)

m=ctModel(Tpoints=10,n.latent=2,n.manifest=3,LAMBDA=matrix(c(1,0,0,0,1,.5),ncol=2),
  # DRIFT=matrix(c(-.3,0,0,-.05),2,2),
  # MANIFESTVAR=diag(.3,3),
  # DIFFUSION=diag(2,2),
  # TDPREDEFFECT=matrix(c(-4,-3,-7,-10),2,2),
  # TRAITVAR='auto',
  CINT=matrix(c(5,3),ncol=1),
  MANIFESTMEANS=matrix(c(0,0,0),ncol=1),
  n.TDpred=2)

f1=ctFit(datawide = gd,m,retryattempts = 0,stationary=c('T0VAR'),objective='Kalman',carefulFit = F,carefulFitWeight = 1)
# f1$mxobj=mxRun(f1$mxobj)
f2=ctFit(datawide = gd,m,retryattempts = 0,stationary=c('T0VAR'),objective='mxRAM',carefulFit = F,carefulFitWeight = 1,
  plotOptimization = F)

expect_equal(f1$mxobj$output$estimate,f2$mxobj$output$estimate)
}
