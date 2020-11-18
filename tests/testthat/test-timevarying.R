if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  # Sys.setenv(NOT_CRAN = 'true')
  
  library(ctsem)
  library(testthat)
  set.seed(1)
  
  context("timevarying")
  
  test_that("varyingLAMBDA", {
    set.seed(1)
    s=list()
    nsubjects=50
    Tpoints=50
    lambdafactor = .3
    dt=1
    
    for(subi in 1:nsubjects){
      gm=ctModel(LAMBDA=diag(2), Tpoints=Tpoints, DRIFT=diag(-.2,2),T0MEANS = matrix(c(3,2)), 
        DIFFUSION=diag(.5,2),
        T0VAR=diag(2))
      d=suppressMessages(ctGenerate(gm,n.subjects = 1,burnin = 3,dtmean = dt))
      d[,'id'] <- subi
      if(subi==1) dat=d else dat=rbind(dat,d)
    }
    dat <- as.matrix(dat)
    dat[,'Y1'] <-  dat[,'Y1'] * (1+ lambdafactor * dat[,'Y2']) #state dependent lambda
    dat[,c('Y1')] <- dat[,c('Y1')] + rnorm(nrow(dat),0,.5) #measurement error
    dat[,c('Y2')] <- dat[,c('Y2')] + rnorm(nrow(dat),0,.5) #measurement error
    
    colnames(dat)[1]='id'
    
    cm <- ctModel(LAMBDA=matrix(c('lbystate * eta2 + 1',0,0,1),2,2),  T0MEANS=c('t0m1','t0m2'),
      T0VAR=matrix(c('t0v11',0,0,'t0v22'),2,2),
      PARS=c('lbystate','lbystate * eta2 + 1'),type='stanct')
    
    cm$pars$indvarying <- FALSE
    # cm$pars$indvarying[cm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    dm<- ctModel(LAMBDA=matrix(c('PARS[2,1]',0,0,1),2,2), T0MEANS=c('t0m1','t0m2|log1p_exp(param)'),
      T0VAR=matrix(c('t0v11',0,0,'t0v22'),2,2),
      PARS=c('lbystate'),type='standt')
    
    dm$pars$indvarying <- FALSE
    # dm$pars$indvarying[dm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    
    
    if(1==99){ 
      f=ctStanFit(datalong = dat,ctstanmodel = cm,optimize=TRUE,
        verbose=0,optimcontrol=list(estonly=F,stochastic=T),savescores = F,nopriors=F)
    }
    
    
    
    for(m in c('cm','dm')){
      argslist <- list(
        ml=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(),
          verbose=0,optimcontrol=list(estonly=F,stochastic=F),savescores = F,nopriors=F)

        #, mlis=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(Jstep=1e-6), 
        #   verbose=0,optimcontrol=list(plot=F,estonly=F,isloops=1,stochastic=F),savescores = F,nopriors=T)
        # ,mapis=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(Jstep=1e-6),
        # verbose=0,optimcontrol=list(plot=F,estonly=F,isloops=1,stochastic=F),savescores = F,nopriors=F)
        # ,hmcintoverpop=list(datalong = dat,ctstanmodel = get(m),optimize=F,iter=500,chains=3, nlcontrol=list(Jstep=1e-6), 
        #   verbose=0,optimcontrol=list(plot=F,estonly=F,stochastic=F),savescores = F,nopriors=F,intoverpop=T),
        # hmc=list(datalong = dat,ctstanmodel = get(m),optimize=F,iter=500,chains=3, nlcontrol=list(),
        # verbose=0,optimcontrol=list(plot=F,estonly=F,stochastic=F),savescores = F,nopriors=F,control=list(max_treedepth=7))
      )
      
      
      for(argi in names(argslist)){
        f = do.call(ctStanFit,argslist[[argi]])
        if(is.null(s[[argi]])) s[[argi]] = list()
        s[[argi]][[m]] <- summary(f,parmatrices=TRUE)
      }
    }
    
    
    dtpars=lapply(s, function(argi) {
      ct=argi$cm$parmatrices[c(grep('(^dt|^asym)',rownames(argi$cm$parmatrices)),
        which(rownames(argi$cm$parmatrices) %in% c('MANIFESTVAR','LAMBDA','T0MEANS','T0VAR'))),c('Mean','Sd')]
      dt=argi$dm$parmatrices[c(grep('(^asym)',rownames(argi$dm$parmatrices)),
        which(rownames(argi$dm$parmatrices) %in% c('DRIFT','CINT','DIFFUSIONcov','MANIFESTVAR','LAMBDA','T0MEANS','T0VAR'))),c('Mean','Sd')]
      rownames(dt)[rownames(dt) %in% c('DRIFT','CINT','DIFFUSIONcov')] <- paste0('dt',rownames(dt)[rownames(dt) %in% c('DRIFT','CINT','DIFFUSIONcov')])
      rownames(dt)[rownames(dt) %in% c('dtDIFFUSIONcov')] <- 'dtDIFFUSION'
      list(ct=ct,dt=dt)
    }
    )
    dtpars=unlist(dtpars,recursive = FALSE)
    dtpars=lapply(dtpars,function(x) x[order(rownames(x)),])
    
    for(ri in 1:nrow(dtpars$ml.ct)){
      for(ci in 1:ncol(dtpars$ml.ct)){
        par=sapply(dtpars, function(y) y[ri,ci])
        for(dimi in 2:length(par)){
          expect_equivalent(par[dimi],par[dimi-1],tol=ifelse(ci==1,1e-1,5e-1))
        }
      }}
    
    
    ll=unlist(lapply(s, function(argi) lapply(argi, function(m) m$loglik)))
    
    for(dimi in 2:length(ll)){
      expect_equivalent(ll[dimi],ll[dimi-1],tol=1e-3)
    }
    
    #do.call(cbind,dtpars)
    
    
    #check time varying lambda estimation
    sapply(s$ml,function(x) expect_equivalent(lambdafactor,x$popmeans[rownames(x$popmeans) %in% 'lbystate','mean'],tol=1e-1))
    
  }) 
}
