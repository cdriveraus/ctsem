# STAYS ON STAN, for now, and this is the most expensive file that does: 425 s
# and eight stan fits. Its C1 blocks compare ctsem's own random-effect
# parameterisation against a hand-written one by reading
# `f2$stanfit$transformedparsfull$pop_T0cov` -- stan's own population
# covariance -- and checking it against `popsd` and `rawpopcorr`. Choosing the
# julia equivalent of that object is not mechanical, and choosing it wrongly
# gives a test that compares the wrong quantity and passes. It converts once
# the population-covariance work settles; see `helper-julia.R` for the
# harness the other expensive files use.
#
# The random-effect equivalence family: what is claimed, and where.
#
# Four files check that ctsem's own parameterisation of a random effect agrees
# with a hand-written model that carries the same effect as an extra latent
# state with ~zero drift: this one (TDPREDEFFECT, LAMBDA, DRIFT, MANIFESTVAR),
# test-ukfpoptest.R (CINT, DRIFT), test-dtVct.R (CINT) and test-timevarying.R
# (state-dependent LAMBDA). Three different claims wear that one costume, and
# they do not need the same amount of data:
#
#   C1 EQUIVALENCE. The two specifications are the same model, so at their
#      optima they give the same maximised loglik, the same population
#      covariance and the same population correlation. This is algebra, not
#      statistics: it does not become more true with more subjects. n only has
#      to be enough to identify the model, so that both optimisers land on the
#      same optimum rather than at two points of one flat ridge. Every block in
#      this file asserts C1, at NEQ subjects.
#
#   C2 RECOVERY. The estimated population sd is the sd of the subjects' actual
#      parameters. This one does need subjects, and one instance of it is worth
#      about as much as six. It is asserted once, at 400 subjects, in
#      randomEffectsDRIFT_julia -- on the julia/Laplace route, because that is
#      the route that gets the nonlinear case right; see the comment there.
#
#   C3 CT VERSUS DT, in test-dtVct.R and test-timevarying.R, is a separate
#      claim and carries its own n.
#
# Why NEQ = 200. Measured, seed 1, this machine: for each design the two fits
# were run at n = 400, 200, 100, 50 and the three C1 gaps recorded. The loglik
# gap is 0 to 3e-3 everywhere except drift at n = 50, where it is 0.149 and
# fails the 1e-1 assertion. What breaks first is the population correlation,
# because at small n those parameters sit on a flat ridge and the two
# parameterisations settle at different points of it:
#
#   design        n=400            n=200            n=100            n=50
#                 dsd    dcorr     dsd    dcorr     dsd    dcorr     dsd    dcorr
#   TDPREDEFFECT  .0032  .0057     .0038  .0041     .0119  .0109 X   .016   .029 X
#   LAMBDA        .0015  .0022     .0039  .0075     .0059  .0054     .0111  .0124 X
#   DRIFT (raw)   .0001  .0019     .0001  .0017     .0001  .0028     .0038  .0148 X
#
# X marks a value past the 1e-2 the assertions use. The first failure is at
# n = 100, so NEQ is twice that. Do not read the small non-monotonicities as
# signal: they are one seed each.
#
# MANIFESTVAR is the exception and its block says why: there only the loglik
# equivalence is assertable.
skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged
  library(ctsem)
  library(testthat)
  cores=2

  # Subjects for every C1 (equivalence) block. See the header for the
  # measurement behind the number.
  NEQ <- 200

  context("randomEffects")

  
  # C1 only. The recovery assertions this block used to carry (population sd
  # and correlation against the generating draws) are C2 and live once, in
  # randomEffectsDRIFT_julia.
  test_that("randomEffectsTDPREDEFFECT", {
    set.seed(1)
    nsubjects <- NEQ
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

    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)

    # checks: C1 equivalence -------------------------------------------------

    #loglik. Measured gap at NEQ: 0 (below print precision) on seed 1.
    testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)

    #sd of ctsem between subjects setup vs manual specification.
    #Measured gap at NEQ: .0038.
    dfsd=data.frame(
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,]))[c(3,1,4)],
      s$popsd[c('mm_Y1','t0m','tdpredeffect'),]) #population estimate
    testthat::expect_true(all(abs(dfsd$f2est - dfsd$X50.) < .01))

    #corr of ctsem between subjects setup vs manual specification.
    #Measured gap at NEQ: .0041.
    dfcorr <- data.frame(
      f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,-2,-2])[lower.tri(diag(3))],
      s$rawpopcorr  )
    testthat::expect_true(all(abs(dfcorr$f2est - dfcorr$X50.) < .01))

  })
  
  
  # C1 only, as above.
  test_that("randomEffectsLambda", {
    set.seed(1)
    nsubjects <- NEQ
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

    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)

    # checks: C1 equivalence -------------------------------------------------

    #loglik. Measured gap at NEQ: 0 (below print precision) on seed 1.
    testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)

    #sd of ctsem between subjects setup vs manual specification.
    #Measured gap at NEQ: .0039.
    dfsd=data.frame(
      f2est=sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,]))[c(3,1,4)],
      s$popsd[c('mm_Y1','t0m','tdpredeffect'),]) #population estimate
    testthat::expect_true(all(abs(dfsd$f2est - dfsd$X50.) < .01))

    #corr of ctsem between subjects setup vs manual specification.
    #Measured gap at NEQ: .0075.
    dfcorr <- data.frame(
      f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,-2,-2])[lower.tri(diag(3))][c(2,1,3)],
      s$rawpopcorr  )
    testthat::expect_true(all(abs(dfcorr$f2est - dfcorr$X50.) < .01))

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

  # C1 only. The recovery half of this block -- the population sd against the
  # sd of the generating draws -- is C2 and is asserted in
  # randomEffectsDRIFT_julia below, which is the arm that gets it right.
  test_that("randomEffectsDRIFT", {
    g <- driftREdata(nsubjects = NEQ)
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
    # At the 400 subjects randomEffectsDRIFT_julia uses, sd(raweffect) is 1.102
    # and sd(effect) is 0.382: mixing the two up is a factor of three, not a
    # rounding difference.  `$rawpopcorr` is a raw quantity despite sitting
    # beside popsd, so the corr checks pair the drift row with raweffect and
    # the sd checks pair it with effect.

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

    #loglik checks. Measured gap at NEQ: .003 on seed 1; at n = 50 it is .149
    #and would fail this, which is the measurement NEQ comes from.
    testthat::expect_true(abs(s$loglik-s2$loglik) < 1e-1)

    #test sd of ctsem between subjects setup vs manual specification, RAW scale.
    #Both sides are deterministic functions of the fitted parameters and agree
    #to ~1e-5, so 1e-2 is loose.  expect_length first because a fit without a
    #`$stanfit` (any non-stan backend) makes the right hand side length zero,
    #and all(logical(0)) is TRUE -- the comparison would pass having tested
    #nothing.
    #Measured gap at NEQ: 6.9e-5.
    rawpopsd_f <- c(f$stanfit$transformedparsfull$rawpopsd) * c(10,1,10)
    testthat::expect_length(rawpopsd_f, 3L)
    testthat::expect_true(all(abs(sqrt(diag(f2$stanfit$transformedparsfull$pop_T0cov[1,,])) -
        rawpopsd_f) < 1e-2))

    #test sd of ctsem between subjects setup vs subject specific pars.  These
    #are shrunken point estimates, so their sd is not the population sd and a
    #band on it would be a shrinkage allowance rather than an identity.  What
    #they do have to do is track the truth:
    testthat::expect_true(all(diag(cor(subjpars,cbind(t0m,effect,baseline))) > .8))

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
      data.frame(
        f2est=cov2cor(f2$stanfit$transformedparsfull$pop_T0cov[1,,])[lower.tri(diag(3))],
        s$rawpopcorr  )

    #test corr of ctsem between subjects setup vs manual specification.
    #Measured gap at NEQ: .0017.
    testthat::expect_true(all(abs(dfcorr$f2est -dfcorr[,'X50.']) <1e-2))

    #The same correlations against the generating draws are C2 and are asserted
    #in randomEffectsDRIFT_julia.

  })

  # The nonlinear transform is where the two ways of integrating the population
  # over a random effect stop agreeing, and this is the block that says which of
  # them is right.  stan makes the random effect a latent state and linearises
  # the DRIFT cell around it; julia's intoverpop='laplace' does the integral
  # without linearising, and reaches a marginal likelihood ~34 nats higher on
  # the same data.  On the population sd of the transformed drift it recovers
  # 0.368 against a true 0.369 (0.367 against 0.370 on a second seed), where the
  # stan route measured 0.311 (-16%) and 0.328 (-11%) on the same data.  Those
  # two stan numbers are a recorded measurement, not a live assertion: the
  # block above is now equivalence-only.
  #
  # THIS IS THE FAMILY'S ONLY RECOVERY TEST (C2) and the only reason it is at
  # 400 subjects rather than NEQ.  It is what stops anyone "fixing" the
  # accurate number to match the biased one, so do not shrink it, and if it
  # starts failing it is the estimate that moved, not a tolerance that was
  # optimistic.
  test_that("randomEffectsDRIFT_julia", {
    skip_without_julia()
    g <- driftREdata()   # nsubjects = 400, deliberately

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

  # randomEffectsDIFFUSION lived here.  It was 130 lines inside `if(F)`: it ran
  # never, asserted nothing, and was read by everyone.  It is not revived,
  # because the equivalence it checked beyond the loglik cannot hold reliably
  # -- see the note on randomEffectsMANIFESTVAR, which has the same shape and
  # the same flat direction.  Measured at n = 400/200/100/50 the loglik gap is
  # 0 at every n while the population correlation gap is 0.3 to 1.9.  The study
  # itself is kept, runnable, at dev/simstudies/simstudy-randomeffects-
  # variance-matrices.R.

  # The exception in this file: here only the LOGLIK half of C1 is assertable,
  # and that is why this block's other checks were sitting inside `if(F)`.
  #
  # A random effect on a variance matrix through log1p_exp leaves the
  # population sd of that effect close to unidentified -- the fit says so
  # itself, warning that "some direction of this model is close to
  # unidentified", and summary() reports the errsd population sd with a 95%
  # interval of 0 to 33.  The two parameterisations therefore agree on the
  # likelihood exactly and can still settle at different points of that flat
  # direction.  Measured at n = 100, two runs of the same data and models
  # differing only in the random starting values: the loglik gap was 0 in both,
  # while the population covariance agreed to 3e-6 in one run and differed by
  # 0.04 in the errsd variance in the other, everything else agreeing to 2e-5.
  # So an assertion on the population covariance here would pass or fail on the
  # starting values.  It is the loglik that carries the claim, and it does so
  # at every n measured (400, 200, 100, 50: gap 0 at each).
  test_that("randomEffectsMANIFESTVAR", {
    set.seed(1)
    nsubjects <- NEQ
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

    f2 <- ctFit(datalong = dat,model= m2,cores=cores)
    s2=summary(f2)

    # checks: C1, loglik half only (see the comment above the block) ---------

    #Measured gap at NEQ: 0, to the precision summary() prints.
    testthat::expect_true(abs(s$loglik -s2$loglik) < .01)

    #A loglik that agrees because both fits failed the same way would pass the
    #line above, so say the fits are real ones.
    testthat::expect_true(is.finite(s$loglik))
    testthat::expect_true(all(is.finite(f$stanfit$rawest)))
    testthat::expect_true(all(is.finite(f2$stanfit$rawest)))

  })
}
