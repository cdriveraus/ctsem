# backend='r' pinned, not left at 'auto'. This test characterises what the model recovers from one particular dataset, and 'auto' prefers the julia engine now, which generates different data for the same seed. The generator is not what is under test here.
skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged
  
  library(ctsem)
  library(testthat)
  cores=2
  
  context("ukfpopcheck")

  # Part of the family whose design is written down at the top of
  # test-tdeffectvariation_covtest.R. Both blocks here are C1 (equivalence):
  # ctsem's own random effect against a hand-written state augmentation, on
  # CINT here and on a nonlinear DRIFT below.
  #
  # Both sizes are left as they are, and that is also measured. For
  # ukfpopcheck1 the loglik, correlation and DRIFT gaps have room at every n
  # tried (n = 160/80/40: loglik 1.0e-4/1.3e-4/2.0e-4 against 1e-3, correlation
  # 7.2e-4/2.6e-4/8.1e-5 against 1e-2, DRIFT 9.8e-5/3.4e-5/2.6e-5 against
  # 1e-2), but the pop_T0cov-to-popcov RATIO check does not: max |ratio - 1| is
  # 0.861 at n = 160, 0.069 at 80 and 0.185 at 40, against a tolerance of 1.
  # It is erratic because the ratio has near-zero elements in its denominator,
  # so it can neither be tightened nor relied on to survive a smaller n. The
  # 12 s a shrink would save is not worth a test that fails on a seed.
  test_that("ukfpopcheck1", {
    set.seed(2)
    Tpoints<-10
    n.latent=2
    n.manifest=3
    nsubjects=160
    burnin=3
    dtmat = matrix(exp(rnorm(burnin+Tpoints-1,-.3,.5)),1)
    par1<-rnorm(nsubjects,-.3,.8)
    par2 <- par1*.4 + rnorm(nsubjects,0,.3)
    for(i in 1:nsubjects){
      gm<-suppressMessages(ctModel(type='omx',n.latent=2,n.manifest=n.manifest,Tpoints=Tpoints,LAMBDA=matrix(c(1,.4,0,0,0,1),3,ncol=2),
        DRIFT=diag(-.4,2),
        CINT=matrix(c(par1[i],par2[i]),2),
        T0VAR=diag(1,2),
        T0MEANS=matrix(c(30,50),n.latent),
        MANIFESTVAR=t(chol(diag(.5,n.manifest))),
        DIFFUSION=t(chol(diag(3,2)))))
      if(i==1) cd<-suppressMessages(ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat,backend='r')) else {
        newdat <- suppressMessages(ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat,backend='r'))
        newdat[,'id'] <- i
        cd<-rbind(cd,newdat)
      }
    }
    

    m1<-ctModel(type='ct',n.latent=4,n.manifest=n.manifest,Tpoints=Tpoints,
      LAMBDA=cbind(gm$LAMBDA,0,0),
      DRIFT=matrix(c(
        'dr11','dr12',1,0,
        'dr21','dr22',0,1,
        0,0,-1e-5,0,
        0,0,0,-1e-5),4,4,byrow=TRUE),
      DIFFUSION=matrix(c(
        'diff11',0,0,0,
        'diff21','diff22',0,0,
        0,0,0,0,
        0,0,0,0),4,4,byrow=TRUE),
      MANIFESTMEANS=matrix(0,nrow=n.manifest)
    )
    
    #model for ukf ctsem
    m2<-ctModel(type='ct',n.latent=2,n.manifest=n.manifest,Tpoints=Tpoints,
      MANIFESTMEANS=matrix(0,n.manifest),
      CINT=matrix(paste0('cint',1:gm$n.latent)),
      LAMBDA=gm$LAMBDA
    )
    
    
    
    #fit 1
    sm1 <- m1
    sm1$pars$indvarying <- FALSE
    sf1 <- ctFit(cd,sm1,cores=cores,verbose=0)
    sf1d=sf1$stanfit$transformedparsfull$pop_DRIFT[1,1:2,1:2]
    sf1ll=sf1$stanfit$optimfit$value
    
    #fit 2
    sm2 <- m2
    sm2$pars$sdscale <- .2
    sf2 <- ctFit(cd,sm2,cores=cores)
    sf2d=sf2$stanfit$transformedparsfull$pop_DRIFT[1,1:2,1:2]
    sf2ll=sf2$stanfit$optimfit$value
    
    test_isclose(sf1ll,sf2ll,tol=1e-3)
    
    test_isclose(sf1$stanfit$transformedparsfull$pop_T0cov[1,,]/
      sf2$stanfit$transformedparsfull$popcov[1,,],
      matrix(1,4,4),tol=1)
    
    test_isclose(cov2cor(sf1$stanfit$transformedparsfull$pop_T0cov[1,,]),
      cov2cor(sf2$stanfit$transformedparsfull$rawpopcov[1,,]),tol=1e-2)

    test_isclose(sf1d,sf2d,tol=1e-2)
    
  })
  
  
  
  test_that("ukfpopcheck2_randomdrift", {
    set.seed(1)
    nsubjects=60
    Tpoints=20
    driftsd=.2
    driftmu= log(.5)
    dt=1
    drift= -2*log1p(exp(2*(rnorm(nsubjects,driftmu,driftsd))))
    
    for(si in 1:nsubjects){
      m=suppressMessages(ctModel(LAMBDA=diag(1), Tpoints=Tpoints, DRIFT=matrix(drift[si]),T0MEANS = matrix(3), 
        T0VAR=matrix(sqrt(.1)), DIFFUSION=diag(1,1), CINT=matrix(-2),MANIFESTVAR=matrix(sqrt(.1))))
      d=suppressMessages(ctGenerate(m,n.subjects = 1,burnin = 0,wide = FALSE,dtmean = dt,
        backend='r'))[,-1]
      if(si==1) dat=cbind(si,d) else dat=rbind(dat,cbind(si,d))
    }
    colnames(dat)[1]='id'
    
    dtm <- ctModel(LAMBDA=diag(1), type='ct',
      DRIFT=matrix('dr11| -2*log1p(exp(-2*param))'),
      CINT=matrix('cint'),
      MANIFESTMEANS = matrix(0)
    )
    
    dtm$pars$indvarying <- FALSE
    dtm$pars$indvarying[dtm$pars$matrix %in% 'DRIFT'] <- TRUE
    
    dtf = ctFit(datalong = dat,model= dtm,optimize=TRUE,
      verbose=0,optimcontrol=list(estonly=F),savescores = F)
    s1=summary(dtf,parmatrices = F,priorcheck = F,residualcov = F)
    
    dtm2 <- ctModel(LAMBDA=matrix(c(1,0),1,2), type='ct',
      DIFFUSION=matrix(c('diff',0,0,0),2,2),DRIFT=matrix(c('-2*log1p(exp(-2*state[2]))',0,0,-.00001),2,2),
      CINT=matrix(c('cint',0),2,1),T0MEANS=matrix(c('t0m1','t0m2|param'),2,1),
      T0VAR=matrix(c('t0v11',0,0,paste0('t0var22 | ',gsub(".* sdscale","",gsub('rawpopsdbase','param',dtm$rawpopsdtransform),fixed=TRUE))),2,2),
      MANIFESTMEANS = matrix(0))
    
    dtm2$pars$indvarying <- FALSE
    dtf2=ctFit(datalong = dat,model= dtm2,optimize = TRUE)
    s2=summary(dtf2,parmatrices = F,priorcheck = F,residualcov = F)
    
    test_isclose(s1$loglik,s2$loglik,tol=1e-3)
    test_isclose(sort(dtf2$stanfit$rawest),sort(dtf$stanfit$rawest),tol=1e-3) #sorting is an ugly hack! could improve...
    
  })
  
}
