if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  cores=2
  
  context("randomEffects")
  

  test_that("randomEffectsTDPREDEFFECT", {
  nsubjects <- 1000
  ntimes <- 20
  
  baseline <- rnorm(nsubjects,2, 2)
  t0m <- rnorm(nsubjects,baseline/2,1)
  effect <- rnorm(nsubjects, 5-baseline/3, 0.5)
  
  for(i in 1:nsubjects){
    gm <- suppressMessages(ctModel(silent=TRUE,Tpoints=ntimes,
      LAMBDA=matrix(c(1,0),1,2), 
      DRIFT= c(-1,1,
        0,-.5),
      T0MEANS = c(t0m[i],0),
      DIFFUSION=c(0.5,0,0,1e-6),
      MANIFESTVAR = 0.5,
      T0VAR = c(0,0,0,0),
      TDPREDMEANS = matrix(c(rep(0,9),1,rep(0,ntimes-10))),
      TDPREDEFFECT = matrix(c(0,effect[i]),2),
      MANIFESTMEANS = baseline[i]))
    
    d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 0,dtmean = 1,logdtsd = 0)))
    d$id <- i
    if(i==1) dat <- d else dat <- rbind(dat,d)
  }
  
  #regular bw effect approach
  m <- ctModel(silent=TRUE,type='stanct',
    LAMBDA=matrix(c(1,0),1,2), 
    DRIFT= c('drift',1,
      0,-0.5),
    T0MEANS = c('t0m',0),
    T0VAR = diag(1e-3,2),
    DIFFUSION=c('diffusion',0,0,0),
    TDPREDEFFECT = matrix(c(0,'tdpredeffect|param|TRUE')))
  
  #manual bw effects
  m2 <- ctModel(silent=TRUE,type='omx',Tpoints=3,
    LAMBDA=matrix(c(1,0,0,0),1,4), 
    DRIFT= c('drift',1,0,0,
      0,-0.5,0,0,
      0,0,-1e-6,0,
      0,0,0,-1e-6),
    DIFFUSION=c('diffusion',0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0),
    T0MEANS = c('t0m',0,'mm','tdpredeffect'),
    TDPREDEFFECT = matrix(c(0,'state[4]',0,0)),
    MANIFESTMEANS='state[3]')
  m2$T0VAR[2,]=m2$T0VAR[,2]=0
  m2=ctStanModel(m2,type='stanct')
  m2$pars$indvarying=F
  
  f <- ctStanFit(datalong = dat,ctstanmodel = m,cores=cores)
  s=summary(f)
  subjpars=ctStanSubjectPars(f)[1,,] #calculate subject specific parameter estimates
  
  f2 <- ctStanFit(datalong = dat,ctstanmodel = m2,cores=cores)
  s2=summary(f2)
  cp2=ctStanContinuousPars(f2)
  
  
  
  # checks ------------------------------------------------------------------
  
  # 
  # plot(subjpars[,1],baseline)
  # abline(0,1)
  # plot(subjpars[,2],t0m)
  # abline(0,1)
  # plot(subjpars[,3],effect)
  # abline(0,1)
  # 
  # 
  # f$stanfit$transformedparsfull$popsd
  # log1p_exp(2*f$stanfit$transformedparsfull$rawpopsdbase-1)
  
  #loglik checks
  testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)
  
  #sd checks
  dfsd=data.frame(trueSample=c(sd(baseline),sd(t0m),sd(effect)),  #sample sd
    subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
    f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,]))[c(3,1,4)],
    s$popsd[c('mm_Y1','t0m','tdpredeffect'),]) #population estimate
  
  #test sd of ctsem between subjects setup vs manual specification
  testthat::expect_true(all(abs(dfsd$f2est - dfsd$X50.) < .01))
  
  #test sd of ctsem between subjects setup vs subject specific pars
  testthat::expect_true(all(abs(dfsd$mean - dfsd$subjPars) < .3))
  
  #test sd of ctsem between subjects setup vs true sample sd -- why is t0means sd overestimated?
  testthat::expect_true(all(abs(dfsd[,'trueSample'] - dfsd[,'X50.']) < .1))
  
  # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,4,4]))) #distribution of pop sd estimates
  # points(density(f$stanfit$transformedpars$popsd[,3]),col=2,type='l') #distribution of pop sd estimates
  
  # #cov checks
  # f$stanfit$transformedparsfull$popcov[1,,]
  # cov(subjpars)
  # cov(cbind(baseline,t0m,effect)) #true sample cov
  # 
  #corr checks
  dfcorr <- data.frame(trueSample=cor(cbind(t0m,baseline,effect))[lower.tri(diag(3))], #true sample cor,
    subjPars=cor(subjpars[,c('t0m','mm_Y1','tdpredeffect')])[lower.tri(diag(3))],
    f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
    f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,-2,-2])[lower.tri(diag(3))],
    s$rawpopcorr  )
  
  #test corr of ctsem between subjects setup vs manual specification
  testthat::expect_true(all(abs(dfcorr$f2est - dfcorr$X50.) < .01))
  
  #test corr of ctsem between subjects setup vs subject specific pars
  # testthat::expect_equivalent(dfcorr$mean,dfcorr$subjPars,tol=.2)
  
  #test corr of ctsem between subjects setup vs true sample sd
  testthat::expect_true(all(abs(dfcorr[,'trueSample'] -dfcorr[,'X50.']) < .1))
  
  
})
  
  
  test_that("randomEffectsLambda", {
    set.seed(1)
    nsubjects <- 1000
    ntimes <- 20
    
    baseline <- rnorm(nsubjects,2, 2)
    t0m <- rnorm(nsubjects,baseline/2,1)
    effect <- rnorm(nsubjects, 5-baseline/3, 0.5)
    
    for(i in 1:nsubjects){
      gm <- suppressMessages(ctModel(silent=TRUE,Tpoints=ntimes,
        LAMBDA=matrix(c(1,effect[i]),1,2), 
        DRIFT= c(-1,0,
          0,-.5),
        T0MEANS = c(t0m[i],0),
        DIFFUSION=c(0.5,0,0,1e-6),
        MANIFESTVAR = 0.5,
        T0VAR = c(0,0,0,0),
        TDPREDMEANS = matrix(c(rep(0,9),1,rep(0,ntimes-10))),
        TDPREDEFFECT = matrix(c(0,1),2),
        MANIFESTMEANS = baseline[i]))
      
      d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 0,dtmean = 1,logdtsd = 0)))
      d$id <- i
      if(i==1) dat <- d else dat <- rbind(dat,d)
    }
    
    #regular bw effect approach
    m <- ctModel(silent=TRUE,type='stanct',
      LAMBDA=matrix(c(1,'tdpredeffect|param|TRUE'),1,2), 
      DRIFT= c('drift',0,
        0,-0.5),
      T0MEANS = c('t0m',0),
      T0VAR = diag(1e-3,2),
      DIFFUSION=c('diffusion',0,0,0),
      TDPREDEFFECT = matrix(c(0,1)))
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='omx',Tpoints=3,
      LAMBDA=matrix(c(1,'state[4]',0,0),1,4), 
      DRIFT= c('drift',0,0,0,
        0,-0.5,0,0,
        0,0,-1e-6,0,
        0,0,0,-1e-6),
      DIFFUSION=c('diffusion',0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0),
      T0MEANS = c('t0m',0,'mm','tdpredeffect'),
      TDPREDEFFECT = matrix(c(0,1,0,0)),
      MANIFESTMEANS='state[3]')
    m2$T0VAR[2,]=m2$T0VAR[,2]=0
    m2=ctStanModel(m2,type='stanct')
    m2$pars$indvarying=F
    
    f <- ctStanFit(datalong = dat,ctstanmodel = m,cores=cores)
    s=summary(f)
    subjpars=ctStanSubjectPars(f)[1,,] #calculate subject specific parameter estimates
    
    f2 <- ctStanFit(datalong = dat,ctstanmodel = m2,cores=cores)
    s2=summary(f2)
    cp2=ctStanContinuousPars(f2)
    
    
    
    # checks ------------------------------------------------------------------
    
    
    # f$stanfit$transformedparsfull$popsd
    # log1p_exp(2*f$stanfit$transformedparsfull$rawpopsdbase-1)
    # 
    
    #loglik checks
    testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)
    
    #sd checks
    dfsd=data.frame(trueSample=c(sd(baseline),sd(t0m),sd(effect)),  #sample sd
      subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,]))[c(3,1,4)],
      s$popsd[c('mm_Y1','t0m','tdpredeffect'),]) #population estimate
    
    #test sd of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(dfsd$f2est - dfsd$X50.) < .01))
    
    #test sd of ctsem between subjects setup vs subject specific pars
    testthat::expect_true(all(abs(dfsd$mean -dfsd$subjPars) < .2))
    
    #test sd of ctsem between subjects setup vs true sample sd -- why is t0means sd overestimated?
    testthat::expect_true(all(abs(dfsd[,'trueSample'] - dfsd[,'X50.']) <.1))
    
    # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,4,4]))) #distribution of pop sd estimates
    # points(density(f$stanfit$transformedpars$popsd[,2]),col=2,type='l') #distribution of pop sd estimates
    
    #cov checks
    
    # f$stanfit$transformedparsfull$popcov[1,,]
    # cov(subjpars)
    # cov(cbind(baseline,t0m,effect)) #true sample cov
    
    #corr checks
    dfcorr <- data.frame(trueSample=cor(cbind(baseline,t0m,effect))[lower.tri(diag(3))][c(3,1,2)], #true sample cor,
      subjPars=cor(subjpars[,c('t0m','tdpredeffect','mm_Y1')])[lower.tri(diag(3))],
      f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
      f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,-2,-2])[lower.tri(diag(3))][c(2,1,3)],
      s$rawpopcorr  )
    
    #test corr of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(dfcorr$f2est - dfcorr$X50.) < .01))
    
    #test corr of ctsem between subjects setup vs subject specific pars
    testthat::expect_true(all(abs(dfcorr$mean - dfcorr$subjPars) <.2))
    
    #test corr of ctsem between subjects setup vs true sample sd
    testthat::expect_true(all(abs(dfcorr[,'trueSample'] - dfcorr[,'X50.']) < .1))
    
    
  })
  
  
  
  
  
  
  test_that("randomEffectsDRIFT", {
    set.seed(1)
    nsubjects <- 400
    ntimes <- 50
    
    baseline <- rnorm(nsubjects,2, 2)
    t0m <- rnorm(nsubjects,baseline/2,1)
    raweffect <- rnorm(nsubjects, baseline/2, 0.5)
    effect <- -log1p(exp(-raweffect))
    
    for(i in 1:nsubjects){
      gm <- suppressMessages(ctModel(silent=TRUE,Tpoints=ntimes,
        LAMBDA=matrix(1), 
        DRIFT= c(effect[i]),
        T0MEANS = c(t0m[i]),
        DIFFUSION=c(0.5),
        MANIFESTVAR = 0.5,
        T0VAR = c(0),
        CINT = c(baseline[i]),MANIFESTMEANS=0))
      
      d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 0,dtmean = 1,logdtsd = 0)))
      d$id <- i
      if(i==1) dat <- d else dat <- rbind(dat,d)
    }
    
    #regular bw effect approach
    m <- ctModel(silent=TRUE,type='stanct',
      CINT='cint',MANIFESTMEANS=0,
      LAMBDA=matrix(1),DRIFT='drift|-log1p_exp(-param)|TRUE')
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='omx',Tpoints=3,
      LAMBDA=matrix(c(1,0,0),ncol=3), 
      DRIFT= c('-log1p_exp(-state[2])',0,0,
        0,-1e-6,0,
        0,0,-1e-6),
      DIFFUSION=c('diffusion',0,0,
        0,0,0,
        0,0,0),
      T0MEANS = c('t0m','drift','cint'),
      CINT=c('state[3]',0,0),MANIFESTMEANS=0)
    m2=ctStanModel(m2,type='stanct')
    m2$pars$indvarying=F
    
    f <- ctStanFit(datalong = dat,ctstanmodel = m,cores=cores)
    s=summary(f)
    s
    subjpars=ctStanSubjectPars(f)[1,,c('T0m_eta1','drift','cint')] #calculate subject specific parameter estimates
    
    f2 <- ctStanFit(datalong = dat,ctstanmodel = m2,cores=cores)
    s2=summary(f2)
    s2
    cp2=ctStanContinuousPars(f2)
    
    
    
    # checks ------------------------------------------------------------------
    
    # mean(effect)
    # 
    # plot(subjpars[,'cint'],baseline)
    # abline(0,1)
    # plot(subjpars[,'drift'],effect)
    # abline(0,1)
    # plot(subjpars[,'T0m_eta1'],t0m)
    # abline(0,1)
    # 
    # plot(subjpars[,c('cint','drift')])
    # 
    # f$stanfit$transformedparsfull$popsd
    # log1p_exp(2*f$stanfit$transformedparsfull$rawpopsdbase-1)
    
    #loglik checks
    testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)
    
    #sd checks (f2 differences expected)
    dfsd=data.frame(trueSample=c(sd(t0m),sd(effect),sd(baseline)),  #sample sd
      subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])),
      s$popsd[c('T0m_eta1','drift','cint'),]) #population estimate
    
    #test sd of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])) -
      f$stanfit$transformedparsfull$rawpopsd *c(10,1,10)) < 1e-2))
    
    #test sd of ctsem between subjects setup vs subject specific pars
    testthat::expect_true(all(abs(dfsd$mean - dfsd$subjPars) < .2*dfsd$mean))
    
    #test sd of ctsem between subjects setup vs true sample sd 
    testthat::expect_true(all(abs(dfsd[,'trueSample'] - dfsd[,'X50.']) < .2*dfsd[,'X50.']))
    # 
    # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2]))) #distribution of pop sd estimates
    # points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
    # 
    # #cov checks
    # f$stanfit$transformedparsfull$popcov[1,,]
    # cov(subjpars)
    # cov(cbind(baseline,t0m,effect)) #true sample cov
    
    #corr checks
    dfcorr <-
      data.frame(trueSample=cor(cbind(t0m,raweffect,baseline))[lower.tri(diag(3))], #true (raw) sample cor,
      subjPars=cor(subjpars)[lower.tri(diag(3))],
      f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
      f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
      s$rawpopcorr  )
    
    #test corr of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(dfcorr$f2est -dfcorr$X50) <1e-2))
    
    #test corr of ctsem between subjects setup vs subject specific pars
    testthat::expect_true(all(abs(dfcorr$mean - dfcorr$subjPars) < .1))
    
    #test corr of ctsem between subjects setup vs true sample sd -- why is t0means sd overestimated?
    testthat::expect_true(all(abs(dfcorr[,'trueSample'] - dfcorr[,'X50.']) <1e-1))
    
    
  })
  
  test_that("randomEffectsDIFFUSION", {
    if(F){ #skip for now
    set.seed(1)
    nsubjects <- 400
    ntimes <- 50
    
    baseline <- rnorm(nsubjects,2, 2)
    t0m <- rnorm(nsubjects,baseline/2,1)
    raweffect <- rnorm(nsubjects,-baseline/3, .1)
    effect <- log1p(exp(raweffect))
    
    for(i in 1:nsubjects){
      gm <- suppressMessages(ctModel(silent=TRUE,Tpoints=ntimes,
        LAMBDA=matrix(1), 
        DRIFT= -1,
        T0MEANS = c(t0m[i]),
        DIFFUSION=effect[i],
        MANIFESTVAR = 0.5,
        T0VAR = c(0),
        CINT = baseline[i],
        MANIFESTMEANS=0))
      
      d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 0,dtmean = .1,logdtsd = 0)))
      d$id <- i
      if(i==1) dat <- d else dat <- rbind(dat,d)
    }
    
    #regular bw effect approach
    m <- ctModel(silent=TRUE,type='stanct',
      T0MEANS='t0m|param',
      MANIFESTVAR=.5,
      MANIFESTMEANS=0,CINT='cint|param',
      LAMBDA=matrix(1),DIFFUSION='diffusion|log1p_exp(param)|TRUE')
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='omx',Tpoints=3,
      LAMBDA=matrix(c(1,0,0),ncol=3), 
      DRIFT= c('drift',0,0,
        0,-1e-12,0,
        0,0,-1e-12),
      DIFFUSION=c('log1p_exp(state[2])',0,0,
        0,0,0,
        0,0,0),
    MANIFESTVAR=.5,
      T0MEANS = c('t0m|param','diffusion|param','cint|param'),
      CINT=c('state[3]',0,0),
      MANIFESTMEANS=0)
    m2$T0VAR[diag(3)==1] <- paste0(m2$T0VAR[diag(3)==1] ,'|log1p_exp(2*param-1)')
    m2=ctStanModel(m2,type='stanct')
    m2$pars$indvarying=F
    
    f <- ctStanFit(datalong = dat,ctstanmodel = m,cores=cores)
    s=summary(f)
    s
    subjpars=ctStanSubjectPars(f)[1,,c('t0m','diffusion','cint')] #calculate subject specific parameter estimates
    
    f2 <- ctStanFit(datalong = dat,ctstanmodel = m2,cores=cores)
    s2=summary(f2)
    s2
    cp2=ctStanContinuousPars(f2)
    
    
    
    # checks ------------------------------------------------------------------
    # mean(effect)
    # 
    # plot(subjpars[,'cint'],baseline)
    # abline(0,1)
    # plot(subjpars[,'diffusion'],effect)
    # abline(0,1)
    # plot(subjpars[,'t0m'],t0m)
    # abline(0,1)
    # 
    # plot(subjpars[,c('diffusion','cint')])
    
    # f$stanfit$transformedparsfull$popsd
    # log1p_exp(2*f$stanfit$transformedparsfull$rawpopsdbase-1)
    
    #loglik checks
    testthat::expect_true(all(abs(s$loglik -s2$loglik) < 1e-2))
    
    #sd checks
    dfsd=data.frame(trueSample=c(sd(t0m),sd(raweffect),sd(baseline)),  #sample sd
      # subjPars=sqrt(diag(cov(subjpars)))[c(2,3,1)], #sd of individual effect point estimates
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])),
      f1est=c(f$stanfit$transformedparsfull$rawpopsd)) #population estimate
    
    dfsdtf=data.frame(trueSample=c(sd(t0m),sd(effect),sd(baseline)),  #sample sd
      subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
      s$popsd) #population estimate
    
    #test sd of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])) -
        sqrt(diag(f$stanfit$transformedparsfull$popcov[1,,])) < .05))
    
    #test sd of ctsem between subjects setup vs subject specific pars
    testthat::expect_equivalent(dfsdtf$X50.,dfsdtf$subjPars,tol=.1)
    
    #test sd of ctsem between subjects setup vs true sample sd 
    testthat::expect_equivalent(dfsd[,'trueSample'],dfsd[,'f1est'],tol=.1)
    testthat::expect_equivalent(dfsdtf$X50.,dfsdtf[,'trueSample'],tol=.1) 
    
    plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2])),type='l') #distribution of pop sd estimates
    points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
    
    #cov checks
    f$stanfit$transformedparsfull$popcov[1,,]
    cov(subjpars)
    cov(cbind(baseline,t0m,effect)) #true sample cov
    
    #rawcorr checks
    dfcorr <- data.frame(trueSample=cor(cbind(t0m,raweffect,baseline))[lower.tri(diag(3))], #true sample cor,
        # subjPars=cor(subjpars[,c(3,1,2)])[lower.tri(diag(3))],
        # f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
        f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
        s$rawpopcorr  )

    
    #test corr of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(dfcorr$f2est -dfcorr$X50) <.05))
    
    #test corr of ctsem between subjects setup vs subject specific pars
    # testthat::expect_true(all(abs(dfcorr$mean - dfcorr$subjPars) < .1))
    
    #test corr of ctsem between subjects setup vs true sample sd -- why is t0means sd overestimated?
    testthat::expect_true(all(abs(dfcorr[,'trueSample'] - dfcorr[,'X50.']) <1e-1))
    
    } #end skip
  })
  
  test_that("randomEffectsMANIFESTVAR", {
    set.seed(1)
    nsubjects <- 400
    ntimes <- 50
    
    baseline <- rnorm(nsubjects,2, 2)
    t0m <- rnorm(nsubjects,baseline/2,1)
    raweffect <- rnorm(nsubjects,-baseline/5, .3)
    effect <- log1p(exp(raweffect))
    
    for(i in 1:nsubjects){
      gm <- suppressMessages(ctModel(silent=TRUE,Tpoints=ntimes,
        LAMBDA=matrix(1), 
        DRIFT= -1,
        T0MEANS = c(t0m[i]),
        MANIFESTVAR=effect[i],
        DIFFUSION = 0.5,
        T0VAR = c(0),
        CINT = baseline[i],
        MANIFESTMEANS=0))
      
      d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 0,dtmean = .1,logdtsd = 0)))
      d$id <- i
      if(i==1) dat <- d else dat <- rbind(dat,d)
    }
    
    #regular bw effect approach
    m <- ctModel(silent=TRUE,type='stanct',
      T0MEANS='t0m|param',
      DIFFUSION=.5,
      MANIFESTMEANS=0,CINT='cint|param',
      LAMBDA=matrix(1),MANIFESTVAR='errsd|log1p_exp(param)|TRUE')
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='omx',Tpoints=3,
      LAMBDA=matrix(c(1,0,0),ncol=3), 
      DRIFT= c('drift',0,0,
        0,-1e-12,0,
        0,0,-1e-12),
      DIFFUSION=c(.5,0,0,
        0,0,0,
        0,0,0),
      MANIFESTVAR='log1p_exp(state[2])',
      T0MEANS = c('t0m|param','errsd|param','cint|param'),
      CINT=c('state[3]',0,0),
      MANIFESTMEANS=0)
    m2$T0VAR[diag(3)==1] <- paste0(m2$T0VAR[diag(3)==1] ,'|log1p_exp(2*param-1)')
    m2=ctStanModel(m2,type='stanct')
    m2$pars$indvarying=F
    
    f <- ctStanFit(datalong = dat,ctstanmodel = m,cores=cores)
    s=summary(f)
    # s
    subjpars=ctStanSubjectPars(f)[1,,c('t0m','errsd','cint')] #calculate subject specific parameter estimates
    
    f2 <- ctStanFit(datalong = dat,ctstanmodel = m2,cores=cores)
    s2=summary(f2)
    # s2
    cp2=ctStanContinuousPars(f2)
    
    
    
    # checks ------------------------------------------------------------------
    
    # mean(effect)
    # 
    # plot(subjpars[,'cint'],baseline)
    # abline(0,1)
    # plot(subjpars[,'errsd'],effect)
    # abline(0,1)
    # plot(subjpars[,'t0m'],t0m)
    # abline(0,1)
    # 
    # plot(subjpars[,c('errsd','cint')])
    # 
    # f$stanfit$transformedparsfull$popsd
    # log1p_exp(2*f$stanfit$transformedparsfull$rawpopsdbase-1)
    
    #loglik checks
    testthat::expect_true(abs(s$loglik -s2$loglik) < .01)
    
    if(F){ #skip for now
    #sd checks
    dfsd=data.frame(trueSample=c(sd(t0m),sd(raweffect),sd(baseline)),  #sample sd
      # subjPars=sqrt(diag(cov(subjpars)))[c(2,3,1)], #sd of individual effect point estimates
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])),
      f1est=c(f$stanfit$transformedparsfull$rawpopsd)) #population estimate
    
    dfsdtf=data.frame(trueSample=c(sd(t0m),sd(effect),sd(baseline)),  #sample sd
      subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
      s$popsd) #population estimate
    
    #test sd of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])) -
      f$stanfit$transformedparsfull$rawpopsd) < .1))
    
    #test sd of ctsem between subjects setup vs subject specific pars
    testthat::expect_equivalent(dfsdtf$X50.,dfsdtf$subjPars,tol=.1)
    
    #test sd of ctsem between subjects setup vs true sample sd 
    testthat::expect_equivalent(dfsd[,'trueSample'],dfsd[,'f1est'],tol=.1)
    testthat::expect_equivalent(dfsdtf$X50.,dfsdtf[,'trueSample'],tol=.1)
    
    plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2])),type='l') #distribution of pop sd estimates
    points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
    
    #cov checks
    f$stanfit$transformedparsfull$popcov[1,,]
    cov(subjpars)
    cov(cbind(baseline,t0m,effect)) #true sample cov
    
    #rawcorr checks
    dfcorr <- data.frame(trueSample=cor(cbind(t0m,raweffect,baseline))[lower.tri(diag(3))], #true sample cor,
      # subjPars=cor(subjpars[,c(3,1,2)])[lower.tri(diag(3))],
      # f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
      f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
      s$rawpopcorr  )
    
    
    
    #test corr of ctsem between subjects setup vs manual specification
    testthat::expect_equivalent(dfcorr$f2est,dfcorr$X50.,tol=.1)
    
    # #test corr of ctsem between subjects setup vs subject specific pars
    # testthat::expect_equivalent(dfcorr$mean,dfcorr$subjPars,tol=.2)
    
    #test corr of ctsem between subjects setup vs true sample sd 
    testthat::expect_equivalent(dfcorr[,'trueSample'],
      dfcorr[,'X50.'],tol=.2)
    
    }
  })
}
