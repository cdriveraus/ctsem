if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  
  library(ctsem)
  library(testthat)
  
  context("ukfpopcheck") 
  
  test_that("ukfpopcheck1", {
    set.seed(1)
    Tpoints<-10
    n.latent=2
    n.manifest=3
    nsubjects=60
    burnin=30
    dtmat = matrix(exp(rnorm(burnin+Tpoints-1,-.3,.5)),1)
    par1<-rnorm(nsubjects,-.3,.8)
    par2 <- par1*.4 + rnorm(nsubjects,0,.3)
    for(i in 1:nsubjects){
      gm<-ctModel(type='omx',n.latent=2,n.manifest=n.manifest,Tpoints=Tpoints,LAMBDA=matrix(c(1,.4,0,0,0,1),3,ncol=2),
        DRIFT=diag(-.4,2),
        CINT=matrix(c(par1[i],par2[i]),2),
        T0VAR=diag(.1,2),
        T0MEANS=matrix(c(3,5),n.latent),
        # MANIFESTMEANS=diag(par1[i],1), #matrix(mmeans,nrow=n.manifest),
        MANIFESTVAR=t(chol(diag(.0001,n.manifest))),
        DIFFUSION=t(chol(diag(3,2))))
      if(i==1) cd<-ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat) else {
        newdat <- ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat)
        newdat[,'id'] <- i
        cd<-rbind(cd,newdat)
      }
    }
    
    cd[,gm$manifestNames]<-(cd[,gm$manifestNames]) + rnorm(length(cd[,gm$manifestNames]),0, .9^2)
    
    
    m1<-ctModel(type='omx',n.latent=4,n.manifest=3,Tpoints=Tpoints,
      LAMBDA=cbind(gm$LAMBDA,0,0),
      MANIFESTMEANS=matrix(0,nrow=n.manifest),
      # CINT=matrix(c('cint1',0),2,1),
      # T0MEANS=matrix(c('t0m1',0),2),
      # T0VAR=matrix(c('t0var11',0, 0,'t0var22'),2,2),
      # DIFFUSION=matrix(c('diff11',0,0,1e-5),2,2),
      # DRIFT=matrix(c('dr11',0,1,-1e-5),2,2)
    )
    m1$DRIFT[3:4,] <- 0
    m1$DRIFT[,3:4] <- 0
    diag(m1$DRIFT)[3:4] <- -1e-5
    m1$DRIFT[1,3] = 1
    m1$DRIFT[2,4] = 1
    m1$DIFFUSION[3:4,] <- 0
    m1$DIFFUSION[,3:4] <- 0
    # m1$T0VAR[3:4,1:2] <- 0
    
    #model for ukf ctsem
    m2<-ctModel(type='omx',n.latent=2,n.manifest=n.manifest,Tpoints=Tpoints,
      # T0MEANS=matrix(c(0),1),
      # T0VAR=matrix(c(1.5),1),
      MANIFESTMEANS=matrix(0,n.manifest),
      CINT=matrix(paste0('cint',1:gm$n.latent)),
      # DIFFUSION=matrix(c(1.4),1),
      # DRIFT=matrix(c(-.4),1),
      TRAITVAR='auto',
      LAMBDA=gm$LAMBDA
    )
    
    
    
    # #original ctsem
    cfit1 <- ctRefineTo(dat = cd,dataform = 'long',ctmodelobj = m1,retryattempts = 1)
    ct1d=cfit1$mxobj$DRIFT$values[1:2,1:2]
    # cfit1$mxobj$DIFFUSION$result
    # cfit1$mxobj$T0VAR$result
    ctll1=cfit1$mxobj$output$fit *-.5
    
    cfit2 <- ctRefineTo(dat = cd,dataform = 'long',ctmodelobj = m2,retryattempts = 0,carefulFit=T,stationary='')
    ct2d=cfit2$mxobj$DRIFT$values
    # cfit2$mxobj$DIFFUSION$result
    # cfit2$mxobj$T0VAR$result
    ctll2=cfit2$mxobj$output$fit *-.5
    
    # summary(cfit1)$ctparameters
    # summary(cfit2)$ctparameters
    
    #bayesian / ukf ctsem
    sm1 <- ctStanModel(m1)
    sm1$pars$indvarying <- FALSE
    # sm1$pars$indvarying[!sm1$pars$matrix %in% c('MANIFESTMEANS')] <- FALSE
    
    # sink(file='../sf1.txt')
    sf1 <- ctStanFit(cd,sm1,iter=200,
      optimize=TRUE,verbose=0,nopriors = TRUE,
      optimcontrol=list(deoptim = F,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:4,
      nlcontrol=list(nldynamics=FALSE))
    # sink()
    # summary(sf1)
    
    sf1d=matrix(sf1$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf1ll=sf1$stanfit$optimfit$value
    
    sf1_derrind<- ctStanFit(cd,sm1,iter=200,
      optimize=TRUE,verbose=0,nopriors = TRUE,
      optimcontrol=list(deoptim = FALSE,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:2,
      nlcontrol=list(nldynamics=FALSE))
    sf1_derrindd=matrix(sf1_derrind$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1_derrind$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf1_derrindll=sf1_derrind$stanfit$optimfit$value
    
    sf1nld<- ctStanFit(cd,sm1,iter=200,
      optimize=TRUE,verbose=0,nopriors = TRUE,
      optimcontrol=list(deoptim = FALSE,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:4,
      nlcontrol=list(nldynamics=TRUE))
    sf1nldd=matrix(sf1nld$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1nld$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf1nldll=sf1nld$stanfit$optimfit$value
    
    # sink(file='../sinkoutg.txt')
    sf1nld_derrind<- ctStanFit(cd,sm1,iter=200,
      optimize=TRUE,verbose=0,nopriors = TRUE,
      init=0,
      optimcontrol=list(deoptim = FALSE,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:2,
      nlcontrol=list(nldynamics=TRUE))
    # sink()
    sf1nld_derrindd=matrix(sf1nld_derrind$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1nld_derrind$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf1nld_derrindll=sf1nld_derrind$stanfit$optimfit$value
    
    
    sf1nl_derrind<- ctStanFit(cd,sm1,iter=200,
      optimize=TRUE,verbose=0,nopriors = TRUE,
      optimcontrol=list(deoptim = FALSE,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:2,
      nlcontrol=list(nldynamics=TRUE,nlmeasurement=TRUE))
    sf1nl_derrindd=matrix(sf1nl_derrind$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1nl_derrind$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf1nl_derrindll=sf1nl_derrind$stanfit$optimfit$value
    
    
    # m2$T0VAR[,]=0
    # m2$T0MEANS[,]=2.588
    sm2 <- ctStanModel(m2)
    sm2$pars$indvarying[!(sm2$pars$matrix %in% c('CINT','T0MEANS'))] <- FALSE
    # sm2$pars$multiplier[(sm2$pars$matrix %in% c('CINT','T0MEANS'))] =1
    # sm2$pars$meanscale[(sm2$pars$matrix %in% c('CINT','T0MEANS'))] =1
    
    # sink(file='../sinkout.txt')
    sf2 <- ctStanFit(cd,sm2,iter=10,fit=T,
      optimize=T,verbose=0,nopriors = TRUE,intoverpop = TRUE,
      # init=0,
      optimcontrol=list(deoptim = FALSE,is=FALSE,isloopsize=50,finishsamples=50),
      derrind=1:2)
    # sink()
    sf2d=matrix(sf2$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf2$stanfit$transformedpars_old)),'mean'],4,4)[1:2,1:2]
    sf2ll=sf2$stanfit$optimfit$value
    # summary(sf2)$popmeans
    
    dvec=c('ct1d','ct2d','sf1d','sf2d','sf1nl_derrindd','sf1nld_derrindd')
    llvec=c('ctll1','ctll2','sf1ll','sf1nldll','sf1nld_derrindll','sf1nl_derrindll')
    
    sapply(dvec,get,envir=sys.frame(sys.parent(0)))
    sapply(llvec,get,envir=sys.frame(sys.parent(0)))
    
    for(di in 2:length(dvec)){
      expect_equivalent(get(dvec[di]),get(dvec[di-1]),tol=1e-1)
      expect_equivalent(get(llvec[di]),get(llvec[di-1]),tol=1e-3)
    }
    
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
      m=ctModel(LAMBDA=diag(1), Tpoints=Tpoints, DRIFT=matrix(drift[si]),T0MEANS = matrix(3), 
        T0VAR=matrix(sqrt(.1)), DIFFUSION=diag(1,1), CINT=matrix(-2),MANIFESTVAR=matrix(sqrt(.1)))
      d=suppressMessages(ctGenerate(m,n.subjects = 1,burnin = 0,wide = FALSE,dtmean = dt))[,-1]
      if(si==1) dat=cbind(si,d) else dat=rbind(dat,cbind(si,d))
    }
    
    colnames(dat)[1]='id'
    
      dtm <- ctModel(LAMBDA=diag(1), type='stanct',
      DRIFT=matrix('dr11| -2*log1p(exp(-2*param))'),
      CINT=matrix('cint'),
      MANIFESTMEANS = matrix(0)
    )
    
    dtm$pars$indvarying <- FALSE
    dtm$pars$indvarying[dtm$pars$matrix %in% 'DRIFT'] <- TRUE
    
    dtf = ctStanFit(datalong = dat,ctstanmodel = dtm,optimize=TRUE,
      verbose=0,optimcontrol=list(estonly=F),savescores = F,nopriors=T)
    
       dtm2 <- ctModel(LAMBDA=matrix(c(1,0),1,2), type='stanct',
      DIFFUSION=matrix(c('diff',0,0,0),2,2),DRIFT=matrix(c('-2*log1p(exp(-2*state[2]))',0,0,-.00001),2,2),
      CINT=matrix(c('cint',0),2,1),T0MEANS=matrix(c('t0m1','t0m2|param'),2,1),
      T0VAR=matrix(c('t0v11',0,0,paste0('t0var22 | ',gsub(".* sdscale","",gsub('rawpopsdbase','param',dtm$rawpopsdtransform),fixed=TRUE))),2,2),
      MANIFESTMEANS = matrix(0))
    
    dtm2$pars$indvarying <- FALSE
    dtf2=ctStanFit(datalong = dat,ctstanmodel = dtm2,optimize = TRUE,nopriors = T)
    
    expect_equivalent(dtf2$stanfit$rawest,dtf$stanfit$rawest,tol=1e-2)

  })
  
}
