require(ctsem)
require(testthat)

context("differenttraits")

test_that("time calc", {
set.seed(4)
Tpoints<-30
n.manifest=2
nsubjects=1000
n.latent=2

DRIFT=matrix(c(-.3, .2, 0, -0.5), byrow=TRUE, nrow=n.latent, ncol=n.latent)

genm=ctModel(Tpoints=Tpoints,
  n.latent=n.latent, n.manifest=n.manifest,
  LAMBDA=matrix(c(1, 0,0,1), nrow=n.manifest, ncol=n.latent),
  DRIFT=DRIFT,
  DIFFUSION=matrix(c(2, 0, 0, 1), byrow=TRUE, nrow=n.latent, ncol=n.latent),
  MANIFESTVAR=matrix(c(1, 0,0,.5), nrow=n.manifest, ncol=n.manifest),
  TRAITVAR=matrix(c(1,.5,0,.8),n.latent,n.latent))

cd=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=50, dtmean=.9, 
  logdtsd=0,simultdpredeffect=TRUE,wide=TRUE)

long=ctWideToLong(datawide = cd,Tpoints = Tpoints,n.manifest = n.manifest)
long=ctDeintervalise(datalong = long)
long=long[-seq(3,length(long),3),]
wide=ctLongToWide(datalong = long,id='id',time='time',manifestNames= genm$manifestNames)

Tpoints=20
wide=ctIntervalise(datawide = wide,Tpoints = Tpoints,n.manifest = n.manifest)

mltrait<-ctModel(Tpoints=Tpoints,n.latent=n.latent,n.manifest=n.manifest,
  LAMBDA=diag(1,n.manifest),
  TRAITVAR='auto')
mmtrait<-ctModel(Tpoints=Tpoints,n.latent=n.latent,n.manifest=n.manifest,
  LAMBDA=diag(1,n.manifest),
  MANIFESTTRAITVAR='auto')
mptrait<-ctModel(Tpoints=Tpoints,n.latent=4,n.manifest=n.manifest,
  LAMBDA=matrix(c(1,0, 0,1, 0,0, 0,0),2,4),
  DRIFT=matrix(c(
    'dr11','dr12',1,0,
    'dr21','dr22',0,1,
    0,0,.0001,0,
    0,0,0,.0001),byrow=TRUE,4,4),
  DIFFUSION=matrix(c(
    'df11',0,0,0,
    'df21','df22',0,0,
    0,0,.0001,0,
    0,0,0,.0001),byrow=TRUE,4,4),
  T0MEANS=matrix(c('t1','t2',0,0),ncol=1))

fmlstrait=ctFit(datawide = wide,ctmodelobj = mltrait,retryattempts = 5,stationary='T0TRAITEFFECT')
fmlstrait=ctFit(datawide = wide,ctmodelobj = mltrait,retryattempts = 5,stationary='')
fmmtrait=ctFit(datawide = wide,ctmodelobj = mmtrait,retryattempts = 5)
fmptrait=ctFit(datawide = wide,ctmodelobj = mptrait,retryattempts = 5)
summary(fmltrait,verbose=TRUE)
summary(fmmtrait)
summary(fmptrait)

#check traits using different fit approaches
expect_equal(rep(0,4),c(fmltrait$mxobj$DRIFT$values-fmmtrait$mxobj$DRIFT$values),tolerance=1e-2)
expect_equal(rep(0,4),c(fmlstrait$mxobj$DRIFT$values-fmptrait$mxobj$DRIFT$values[1:2,1:2]),tolerance=1e-2)

#check DRIFT is reasonably estimated
expect_equal(rep(0,4),c(fmltrait$mxobj$DRIFT$values-DRIFT),tolerance=.1)

})
