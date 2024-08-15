if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  context("ctBinaryGaussianMix") #develop some expectations here!
  
  test_that("ctBinaryGaussianMix1", {
    set.seed( 1234 )
    cores=2

    invlog=function (x) exp(x)/(1 + exp(x))
    
    #gen data
    gm <- ctModel(DRIFT= c(-.2, .2, 
      0,-.1),
      DIFFUSION=c(.3,0,
        0,.4), 
      CINT=c(.1,.1),
      # TRAITVAR=diag(.3,2), #old approach to allow individual variation 
      LAMBDA= diag(1,2),
      n.latent=2,n.manifest=2,Tpoints=50)
    
    d=ctGenerate(gm,n.subjects = 50,logdtsd=.2,dtmean = .2,burnin = 20)
    d[,gm$manifestNames[1]] <- d[,gm$manifestNames[1]] + rnorm(nrow(d),0,.2)
    d=data.frame(d)
    for(i in 1:10){
      d[[paste0('b',i)]] <- rbinom(nrow(d),size=1,prob=invlog(d[,gm$manifestNames[2]]))
    }
    
   
    
    # plot(invlog(d[,gm$manifestNames[2]])[1:100],type='l',col=2)
    # points( apply(d[,paste0('b',1:10)],1,function(x) mean(x))[1:100],type='l')

    MANIFESTVAR = diag(c(1,rep(0,10)),11)
    MANIFESTVAR[1]='mvar1'
    m <- ctModel(type='stanct',
      manifestNames = c('Y1',paste0('b',1:10)),
      LAMBDA=rbind(diag(1,2),cbind(rep(0,9),rep(1,9))),
      MANIFESTMEANS = 0,
      MANIFESTVAR = MANIFESTVAR,
      CINT=c('CINT1','cint2'))
    m$manifesttype[2:11]=1 #set type to binary
    m$pars$indvarying=F
    
    #fit with integration (linearised approximation)
    f <- ctStanFit( datalong = d, ctstanmodel = m,cores=cores,plot=10)

  
    #test if the estimated model pars 95% confidence intervals contain true pars
    lowmats <- ctStanContinuousPars(f,calcfuncargs = list(probs=.025))
    upmats <- ctStanContinuousPars(f,calcfuncargs = list(probs=.975))
    
    gmn <- ctsem:::ctModeltoNumeric(gm)
    gmn$DIFFUSIONcov <- tcrossprod(gmn$DIFFUSION)
    
    mats <- c('DRIFT','DIFFUSIONcov','CINT')
    for(i in 1:length(mats)){
      expect_true(all(lowmats[[mats[i]]] < gmn[[mats[i]]] & upmats[[mats[i]]] > gmn[[mats[i]]]))
    }
    
  })
  
}
