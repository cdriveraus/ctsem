if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  cores=2
  
  context("randomEffects")
  
  
  test_that("randomEffectsTDPREDEFFECT", {
    set.seed(1)
    nsubjects <- 1000
    ntimes <- 20
    
    baseline <- rnorm(nsubjects,2, 2)
    t0m <- rnorm(nsubjects,baseline/2,1)
    effect <- rnorm(nsubjects, 5-baseline/3, 0.5)
    
    for(i in 1:nsubjects){
      gm <- suppressMessages(ctModel(type='omx',silent=TRUE,Tpoints=ntimes,
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
    m <- ctModel(silent=TRUE,type='ct',
      LAMBDA=matrix(c(1,0),1,2), 
      DRIFT= c('drift',1,
        0,-0.5),
      T0MEANS = c('t0m',0),
      T0VAR = diag(1e-3,2),
      DIFFUSION=c('diffusion',0,0,0),
      TDPREDEFFECT = matrix(c(0,'tdpredeffect|param|TRUE')))
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='ct',Tpoints=3,
      LAMBDA=matrix(c(1,0,0,0),1,4), 
      DRIFT= c('drift',1,0,0,
        0,-0.5,0,0,
        0,0,-1e-6,0,
        0,0,0,-1e-6),
      DIFFUSION=c('diffusion',0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0),
      T0VAR=matrix(c(
        't0v11',0,0,0,
        0,0,0,0,
        't0v31',0,'t0v33',0,
        't0v41',0,'t0v43','t0v44'),4,4,byrow=TRUE),
      T0MEANS = c('t0m',0,'mm','tdpredeffect'),
      TDPREDEFFECT = matrix(c(0,'state[4]',0,0)),
      MANIFESTMEANS='state[3]')
    m2$pars$indvarying=F
    
    f <- ctFit(datalong = dat,model= m,cores=cores)
    s=summary(f)
    subjpars=ctSubjectPars(f)[1,,] #calculate subject specific parameter estimates
    
    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)
    cp2=ctSummaryMatrices(f2)
    
    
    
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
    # test_isclose(dfcorr$mean,dfcorr$subjPars,tol=.2)
    
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
      gm <- suppressMessages(ctModel(type='omx',silent=TRUE,Tpoints=ntimes,
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
    m <- ctModel(silent=TRUE,type='ct',
      LAMBDA=matrix(c(1,'tdpredeffect|param|TRUE'),1,2), 
      DRIFT= c('drift',0,
        0,-0.5),
      T0MEANS = c('t0m',0),
      T0VAR = diag(1e-3,2),
      DIFFUSION=c('diffusion',0,0,0),
      TDPREDEFFECT = matrix(c(0,1)))
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='ct',Tpoints=3,
      LAMBDA=matrix(c(1,'state[4]',0,0),1,4), 
      DRIFT= c('drift',0,0,0,
        0,-0.5,0,0,
        0,0,-1e-6,0,
        0,0,0,-1e-6),
      DIFFUSION=c('diffusion',0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0),
      T0VAR=matrix(c(
        't0v11',0,0,0,
        0,0,0,0,
        't0v31',0,'t0v33',0,
        't0v41',0,'t0v43','t0v44'),4,4,byrow=TRUE),
      T0MEANS = c('t0m',0,'mm','tdpredeffect'),
      TDPREDEFFECT = matrix(c(0,1,0,0)),
      MANIFESTMEANS='state[3]')
    m2$pars$indvarying=F
    
    f <- ctFit(datalong = dat,model= m,cores=cores)
    s=summary(f)
    subjpars=ctSubjectPars(f)[1,,] #calculate subject specific parameter estimates
    
    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)
    cp2=ctSummaryMatrices(f2)
    
    
    
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
  
  
  
  
  
  
  # The DRIFT population, shared by the two DRIFT blocks below so that they
  # describe the same subjects.  The drift is drawn on the RAW scale and pushed
  # through the model's own transform, so `raweffect` and `effect` are different
  # quantities -- sd 1.102 against sd 0.382 for seed 1 -- and every comparison
  # below has to say which of the two it is using.
  driftREdata <- function(nsubjects=400, ntimes=50, seed=1){
    set.seed(seed)
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
    list(dat=dat, baseline=baseline, t0m=t0m, raweffect=raweffect, effect=effect)
  }

  #regular bw effect approach
  driftREmodel <- function() ctModel(silent=TRUE,type='ct',
    CINT='cint',MANIFESTMEANS=0,
    LAMBDA=matrix(1),DRIFT='drift|-log1p_exp(-param)|TRUE')

  test_that("randomEffectsDRIFT", {
    g <- driftREdata()
    dat <- g$dat; baseline <- g$baseline; t0m <- g$t0m
    raweffect <- g$raweffect; effect <- g$effect

    m <- driftREmodel()

    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='ct',Tpoints=3,
      LAMBDA=matrix(c(1,0,0),ncol=3), 
      DRIFT= c('-log1p_exp(-state[2])',0,0,
        0,-1e-6,0,
        0,0,-1e-6),
      DIFFUSION=c('diffusion',0,0,
        0,0,0,
        0,0,0),
      T0MEANS = c('t0m','drift','cint'),
      CINT=c('state[3]',0,0),MANIFESTMEANS=0)
    m2$pars$indvarying=F
    
    f <- ctFit(datalong = dat,model= m,cores=cores)
    s=summary(f)
    #s
    subjpars=ctSubjectPars(f)[1,,c('T0m_eta1','drift','cint')] #calculate subject specific parameter estimates
    
    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)
    #s2
    cp2=ctSummaryMatrices(f2)
    
    
    
    # checks ------------------------------------------------------------------

    # Scales.  `drift` is the one parameter in this file whose transform is
    # nonlinear, so raw and transformed are different quantities for it:
    #
    #   raw          the parameter the population distribution is normal on.
    #                f2 carries it as state 2 of pop_T0cov; f reports it as
    #                `rawpopsd`, before the multiplier that T0MEANS and CINT
    #                carry by default (10) and that DRIFT's custom transform
    #                does not (1) -- hence the c(10,1,10) below.  sd(raweffect)
    #                is the truth on this scale.
    #   transformed  the drift itself.  `summary()$popsd` is on this scale on
    #                both backends -- stan by Monte Carlo over the raw
    #                population, julia by quadrature over the same -- and so is
    #                sd(effect).  `ctSubjectPars` is on this scale too.
    #
    # sd(raweffect) is 1.102 and sd(effect) is 0.382 here: mixing the two up is
    # a factor of three, not a rounding difference.  `$rawpopcorr` is a raw
    # quantity despite sitting beside popsd, so the corr checks below pair the
    # drift row with raweffect and the sd checks pair it with effect.

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

    #sd checks -- TRANSFORMED scale in every column (f2's raw sds are checked
    #against f's raw sds separately, below, and are not comparable to these)
    dfsd=data.frame(trueSample=c(sd(t0m),sd(effect),sd(baseline)),  #sample sd
      subjPars=sqrt(diag(cov(subjpars))), #sd of individual effect point estimates
      s$popsd[c('T0m_eta1','drift','cint'),]) #population estimate

    #test sd of ctsem between subjects setup vs manual specification, RAW scale.
    #Both sides are deterministic functions of the fitted parameters and agree
    #to ~1e-5, so 1e-2 is loose.  expect_length first because a fit without a
    #`$stanfit` (any non-stan backend) makes the right hand side length zero,
    #and all(logical(0)) is TRUE -- the comparison would pass having tested
    #nothing.
    rawpopsd_f <- c(f$stanfit$transformedparsfull$rawpopsd) * c(10,1,10)
    testthat::expect_length(rawpopsd_f, 3L)
    testthat::expect_true(all(abs(sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])) -
        rawpopsd_f) < 1e-2))

    #test sd of ctsem between subjects setup vs subject specific pars.  These
    #are shrunken point estimates, so their sd is not the population sd and the
    #band is a shrinkage allowance rather than an identity.  What they do have
    #to do is track the truth, which the band alone does not check:
    testthat::expect_true(all(diag(cor(subjpars,cbind(t0m,effect,baseline))) > .8))
    testthat::expect_true(all(abs(dfsd$mean - dfsd$subjPars) < .2*dfsd$mean))

    #test sd of ctsem between subjects setup vs true sample sd, per row.
    #T0MEANS and CINT transform linearly and come back within a few percent.
    #DRIFT does not: stan integrates the population over the nonlinearity by
    #turning the random effect into a latent state and linearising the DRIFT
    #cell around it, which shrinks the reported sd.  Measured here, against a
    #true population sd of 0.369: stan 0.311 (-16%) on this seed and 0.328
    #against 0.370 (-11%) on another, while julia's intoverpop='laplace', which
    #does not linearise, returns 0.368 and 0.367.  The .3 is what the stan
    #route can deliver on this row, not the accuracy of the quantity -- the
    #randomEffectsDRIFT_julia block holds the same row to .1.
    testthat::expect_true(all(abs(dfsd[,'trueSample'] - dfsd[,'X50.']) <
        c(.1,.3,.1)*dfsd[,'X50.']))
    #
    # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2]))) #distribution of pop sd estimates
    # points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
    # 
    # #cov checks
    # f$stanfit$transformedparsfull$popcov[1,,]
    # cov(subjpars)
    # cov(cbind(baseline,t0m,effect)) #true sample cov
    
    #corr checks -- RAW scale on every side, so the drift row pairs with
    #raweffect rather than effect
    dfcorr <-
      data.frame(trueSample=cor(cbind(t0m,raweffect,baseline))[lower.tri(diag(3))], #true (raw) sample cor,
        subjPars=cor(subjpars)[lower.tri(diag(3))],
        f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
        f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
        s$rawpopcorr  )

    #test corr of ctsem between subjects setup vs manual specification
    testthat::expect_true(all(abs(dfcorr$f2est -dfcorr[,'X50.']) <1e-2))

    #test corr of ctsem between subjects setup vs subject specific pars
    testthat::expect_true(all(abs(dfcorr$mean - dfcorr$subjPars) < .1))

    #test corr of ctsem between subjects setup vs true sample sd -- why is t0means sd overestimated?
    testthat::expect_true(all(abs(dfcorr[,'trueSample'] - dfcorr[,'X50.']) <1e-1))


  })

  # The nonlinear transform is where the two ways of integrating the population
  # over a random effect stop agreeing, and this is the block that says which of
  # them is right.  stan makes the random effect a latent state and linearises
  # the DRIFT cell around it; julia's intoverpop='laplace' does the integral
  # without linearising, and reaches a marginal likelihood ~34 nats higher on
  # the same data.  On the population sd of the transformed drift it recovers
  # 0.368 against a true 0.369 (0.367 against 0.370 on a second seed), where the
  # stan route in the block above gives 0.311 and 0.328.  Hence the tolerance
  # here is .1 rather than the .3 that row needs on stan: if this block starts
  # failing, it is the estimate that moved, not the tolerance that was optimistic.
  test_that("randomEffectsDRIFT_julia", {
    skip_without_julia()
    g <- driftREdata()

    f <- ctFit(datalong = g$dat, model = driftREmodel(), cores = cores,
      backend = 'julia', intoverpop = 'laplace')
    s <- summary(f)

    #sd checks -- TRANSFORMED scale on both sides: sd(effect) is the spread of
    #the subjects' own drifts, s$popsd the population sd of that same quantity
    dfsd <- data.frame(trueSample = c(sd(g$t0m), sd(g$effect), sd(g$baseline)),
      s$popsd[c('T0m_eta1','drift','cint'),])
    testthat::expect_true(all(abs(dfsd[,'trueSample'] - dfsd[,'X50.']) <
        .1*dfsd[,'X50.']))

    #corr checks -- RAW scale, so the drift row pairs with raweffect
    dfcorr <- data.frame(
      trueSample = cor(cbind(g$t0m, g$raweffect, g$baseline))[lower.tri(diag(3))],
      s$rawpopcorr)
    testthat::expect_true(all(abs(dfcorr[,'trueSample'] - dfcorr[,'X50.']) < .1))
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
      m <- ctModel(silent=TRUE,type='ct',
        T0MEANS='t0m|param',
        MANIFESTVAR=.5,
        MANIFESTMEANS=0,CINT='cint|param',
        LAMBDA=matrix(1),DIFFUSION='diffusion|log1p_exp(param)|TRUE')
      
      #manual bw effects
      m2 <- ctModel(silent=TRUE,type='ct',Tpoints=3,
        LAMBDA=matrix(c(1,0,0),ncol=3), 
        DRIFT= c('drift',0,0,
          0,-1e-12,0,
          0,0,-1e-12),
        DIFFUSION=c('log1p_exp(state[2])',0,0,
          0,0,0,
          0,0,0),
        MANIFESTVAR=.5,
        T0MEANS = c('t0m|param','diffusion|param','cint|param'),
        T0VAR=matrix(c(
          't0var11 | log1p_exp(2*param-1)',0,0,
          't0var21','t0var22 | log1p_exp(2*param-1)',0,
          't0var31','t0var32','t0var33 | log1p_exp(2*param-1)'),3,3,byrow=TRUE),
        CINT=c('state[3]',0,0),
        MANIFESTMEANS=0)
      m2$pars$indvarying=F
      
      f <- ctFit(datalong = dat,model= m,cores=cores)
      s=summary(f)
      s
      subjpars=ctSubjectPars(f)[1,,c('t0m','diffusion','cint')] #calculate subject specific parameter estimates
      
      f2 <- ctFit(datalong = dat,model= m2,cores=cores)
      s2=summary(f2)
      s2
      cp2=ctSummaryMatrices(f2)
      
      
      
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
          sqrt(diag(f$stanfit$transformedparsfull$popcov[1,,])) < .05)))
      
      #test sd of ctsem between subjects setup vs subject specific pars
      test_isclose(dfsdtf$X50.,dfsdtf$subjPars,tol=.1)
      
      #test sd of ctsem between subjects setup vs true sample sd 
      test_isclose(dfsd[,'trueSample'],dfsd[,'f1est'],tol=.1)
      test_isclose(dfsdtf$X50.,dfsdtf[,'trueSample'],tol=.1) 
      # 
      # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2])),type='l') #distribution of pop sd estimates
      # points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
      # 
      # #cov checks
      # f$stanfit$transformedparsfull$popcov[1,,]
      # cov(subjpars)
      # cov(cbind(baseline,t0m,effect)) #true sample cov
      # 
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
    m <- ctModel(silent=TRUE,type='ct',
      T0MEANS='t0m|param',
      DIFFUSION=.5,
      MANIFESTMEANS=0,CINT='cint|param',
      LAMBDA=matrix(1),MANIFESTVAR='errsd|log1p_exp(param)|TRUE')
    
    #manual bw effects
    m2 <- ctModel(silent=TRUE,type='ct',Tpoints=3,
      LAMBDA=matrix(c(1,0,0),ncol=3), 
      DRIFT= c('drift',0,0,
        0,-1e-12,0,
        0,0,-1e-12),
      DIFFUSION=c(.5,0,0,
        0,0,0,
        0,0,0),
      MANIFESTVAR='log1p_exp(state[2])',
      T0MEANS = c('t0m|param','errsd|param','cint|param'),
      T0VAR=matrix(c(
        't0var11 | log1p_exp(2*param-1)',0,0,
        't0var21','t0var22 | log1p_exp(2*param-1)',0,
        't0var31','t0var32','t0var33 | log1p_exp(2*param-1)'),3,3,byrow=TRUE),
      CINT=c('state[3]',0,0),
      MANIFESTMEANS=0)
    m2$pars$indvarying=F
    
    f <- ctFit(datalong = dat,model= m,cores=cores)
    s=summary(f)
    # s
    subjpars=ctSubjectPars(f)[1,,c('t0m','errsd','cint')] #calculate subject specific parameter estimates
    
    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)
    # s2
    cp2=ctSummaryMatrices(f2)
    
    
    
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
      test_isclose(dfsdtf$X50.,dfsdtf$subjPars,tol=.1)
      
      #test sd of ctsem between subjects setup vs true sample sd 
      test_isclose(dfsd[,'trueSample'],dfsd[,'f1est'],tol=.1)
      test_isclose(dfsdtf$X50.,dfsdtf[,'trueSample'],tol=.1)
      # 
      # plot(density(sqrt(f2$stanfit$transformedpars$pop_T0cov[,2,2])),type='l') #distribution of pop sd estimates
      # points(density(f$stanfit$transformedpars$rawpopsd[,2]),col=2,type='l') #distribution of pop sd estimates
      
      # #cov checks
      # f$stanfit$transformedparsfull$popcov[1,,]
      # cov(subjpars)
      # cov(cbind(baseline,t0m,effect)) #true sample cov
      
      #rawcorr checks
      dfcorr <- data.frame(trueSample=cor(cbind(t0m,raweffect,baseline))[lower.tri(diag(3))], #true sample cor,
        # subjPars=cor(subjpars[,c(3,1,2)])[lower.tri(diag(3))],
        # f1popCovbased=cov2cor(f$stanfit$transformedparsfull$popcov[1,,])[lower.tri(diag(3))],
        f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
        s$rawpopcorr  )
      
      
      
      #test corr of ctsem between subjects setup vs manual specification
      test_isclose(dfcorr$f2est,dfcorr$X50.,tol=.1)
      
      # #test corr of ctsem between subjects setup vs subject specific pars
      # test_isclose(dfcorr$mean,dfcorr$subjPars,tol=.2)
      
      #test corr of ctsem between subjects setup vs true sample sd 
      test_isclose(dfcorr[,'trueSample'], dfcorr[,'X50.'],tol=.2)
      
    }
  })
}
