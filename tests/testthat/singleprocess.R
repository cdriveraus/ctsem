require(ctsem)
require(testthat)

context("singleprocess")

test_that("time calc", {
set.seed(4)
Tpoints<-60
n.manifest=1
nsubjects=60
n.latent=1

DRIFT=matrix(c(-.3), nrow=n.latent, ncol=n.latent)

genm=ctModel(Tpoints=Tpoints,
  n.latent=n.latent, n.manifest=n.manifest,
  LAMBDA=matrix(c(1), nrow=n.manifest, ncol=n.latent),
  DRIFT=DRIFT,
  T0MEANS=matrix(3),
  MANIFESTTRAITVAR=diag(2,1),
  T0VAR=matrix(6),
  DIFFUSION=matrix(c(2), byrow=TRUE, nrow=n.latent, ncol=n.latent),
  MANIFESTVAR=matrix(c(2), nrow=n.manifest, ncol=n.manifest))

cd=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=0, dtmean=1, 
  logdtsd=0,wide=TRUE)

long=ctWideToLong(datawide = cd,Tpoints = Tpoints,n.manifest = n.manifest)
long=ctDeintervalise(datalong = long)
long=long[-seq(3,length(long),3),]
wide=ctLongToWide(datalong = long,id='id',time='time',manifestNames= genm$manifestNames)
Tpoints=Tpoints*2/3
wide=ctIntervalise(datawide = wide,Tpoints = Tpoints,n.manifest = n.manifest)

sm<-ctModel(Tpoints=Tpoints,type='stanct',
  n.latent=n.latent,n.manifest=n.manifest,
  CINT=matrix('cint'),
  MANIFESTVAR=diag(2,1),
  MANIFESTMEANS=matrix(0),
  LAMBDA=diag(1,n.manifest))

sm$pars$indvarying[!(sm$pars$matrix %in% c('T0MEANS','CINT'))]=FALSE
# sm$pars$transform[(sm$pars$matrix %in% c('DRIFT'))]='-exp(param)-.0001'
# sm$pars$transform[(sm$pars$matrix %in% c('DIFFUSION','T0VAR','MANIFESTVAR'))]='log(exp(param)+1)*500'
# sm$pars$sdscale=1
# sm$pars$indvarying=FALSE


# sm$pars$transform[5]='exp(param*10) +.00001'

sf=ctStanFit(datalong=long,ctstanmodel=sm,
  # stanmodeltext=sm2,
  iter=300,chains=3,plot=FALSE,fit=TRUE,optimize=FALSE,kalman=TRUE)
sfsum <- summary(sf)


m<-ctModel(Tpoints=Tpoints,type='omx',
  n.latent=n.latent,n.manifest=n.manifest,
  # CINT=matrix('cint'),
  MANIFESTVAR=diag(2,1),
  TRAITVAR='auto',
  # MANIFESTMEANS=matrix(0),
  LAMBDA=diag(1,n.manifest))

f=ctFit(dat=wide,ctmodelobj=m)
summary(f)$ctparameters

stantraitvar <- (-1/(sfsum$popmeans[2,1])) %*% sfsum$popsd[2,1]^2 %*% (-1/(sfsum$popmeans[2,1]))

#check drift using different fit approaches
expect_equal(sfsum$popmeans[2,1],summary(f)$ctparameters[2,1],tolerance=1e-2)

})
