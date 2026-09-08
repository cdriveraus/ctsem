skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged
  
  library(ctsem)
  library(testthat)
  set.seed(2)
  
  context("tipredcheck")
  
  test_that("simpleTIpredcheck", {
    Tpoints=10
    n.manifest=1
    n.TDpred=0
    n.TIpred=1
    n.latent=1
    n.subjects=50
    TI1 <- rnorm(n.subjects)
    gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
      n.TDpred=n.TDpred,n.manifest=n.manifest,
      MANIFESTVAR=diag(0.5,1),
      LAMBDA=diag(1,1),T0MEANS=100,
      DRIFT=matrix(c(-.3),nrow=1),
      DIFFUSION=matrix(c(2),1),
      T0VAR=diag(10,1))
    
    for(i in 1:n.subjects){
      gm$CINT[1,1] <- TI1[i]*5+rnorm(1,0,.6)
      ndat<-suppressMessages(ctGenerate(gm,n.subjects=1,burnin=10,logdtsd=.4))
      ndat <- cbind(ndat,TI1[i])
      ndat[,1] <- i
      if(i>1) tdat <- rbind(tdat,ndat) else tdat <- ndat
    }
    colnames(tdat)[4] <- 'TI1'
    
    tdat[2,'Y1'] <- NA
    tdat[tdat[,'id']==2,'TI1'] <- NA
    
    checkm<-suppressMessages(ctModel(type='ct',Tpoints=Tpoints,
      MANIFESTVAR=diag(0.5,1),
      DRIFT=matrix(c(-.3),nrow=1),
      DIFFUSION=matrix(c(2),1),
      n.latent=n.latent,n.TDpred=n.TDpred,
      n.TIpred=n.TIpred,
      MANIFESTMEANS=matrix(0,nrow=n.manifest),
      CINT=matrix(c('cint1'),ncol=1),
      n.manifest=n.manifest,LAMBDA=diag(1)))
    
    # checkm$pars$indvarying <- FALSE
    
    checkm$pars[c(-1,-7) ,c('TI1_effect')] <- FALSE
    
    tfit1<-ctFit(tdat,checkm,optimize=TRUE,
      optimcontrol=list(uncertainty='is',carefulfit=F),
      priors=TRUE,verbose=0)
    s1=summary(tfit1)
    
    test_isclose(s1$tipreds[2,'mean'],5,tol=.2)
    test_isclose(s1$popsd[2,'50%'],.6,tol=.2)
    
    tfit2<-ctFit(tdat,checkm,optimize=TRUE,cores=2,verbose=0,priors=TRUE)
    s2=summary(tfit2)
    
    test_isclose(s2$tipreds[2,'mean'],5,tol=.2)
    test_isclose(s2$popsd[2,'50%'],.6,tol=.2)
    
    tfit3<-suppressWarnings(ctFit(tdat,checkm,iter=300,chains=2,optimize=FALSE,
      sampleControl=list(adapt_delta=.8,max_treedepth=6),plot=FALSE))
    s3=summary(tfit3)
    
    test_isclose(s3$tipreds[2,'mean'],5,tol=.5)
    test_isclose(s3$popsd[2,'50%'],.6,tol=.5)
  })

  #stanoptimis() under estonly=TRUE used to return no $standata, so ctFit's
  #copy-back was a no-op and the TIpred block hit apply(NULL,1,sum).
  test_that("estonly stan fit reports tipred effects", {
    m <- suppressMessages(ctModel(type='ct',n.latent=2,n.manifest=2,
      LAMBDA=diag(2),manifestNames=c('Y1','Y2'),
      TIpredNames=c('TI1','TI2')))
    m$pars$indvarying <- FALSE

    f <- ctFit(ctstantestdat, m, backend='stan', cores=1, verbose=0,
      optimcontrol=list(estonly=TRUE, carefulfit=FALSE))

    expect_s3_class(f, 'ctStanFit')
    expect_true(sum(f$setup$matsetup$tipred) > 0)
    expect_true(f$standata$ntipredeffects > 0)

    #one flagged population parameter per nonzero row of TIPREDEFFECTsetup
    ms <- f$setup$matsetup
    expect_equal(
      length(unique(ms$param[ms$tipred == 1L])),
      sum(apply(f$standata$TIPREDEFFECTsetup,1,sum) > 0))

    #and the estonly path agrees with the full one
    ffull <- ctFit(ctstantestdat, m, backend='stan', cores=1, verbose=0,
      optimcontrol=list(carefulfit=FALSE))
    expect_equal(ms$tipred, ffull$setup$matsetup$tipred)
    expect_equal(f$standata$ntipredeffects, ffull$standata$ntipredeffects)
  })
}
