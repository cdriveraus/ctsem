if(identical(Sys.getenv("NOT_CRAN"), "true") & .Machine$sizeof.pointer != 4){
  # Sys.setenv(NOT_CRAN='true')

  set.seed(1)
  library(ctsem)
  library(testthat)
  
  context("dtVct_lVnl")
  
  test_that("dtVct_CINTheterogeneity", {
    set.seed(1)
    s=list()
    nsubjects=200
    Tpoints=15
    parsd=1.4
    parmu= -10.4
    dt=1
    par= (rnorm(nsubjects,parmu,parsd))
    mean(par)
    sd(par)
    
    for(subi in 1:nsubjects){
      gm=ctModel(LAMBDA=diag(1), Tpoints=Tpoints, DRIFT=matrix(-1),T0MEANS = matrix(4), 
        CINT=matrix(par[subi]),DIFFUSION=matrix(2),
        T0VAR=matrix(2), MANIFESTVAR=matrix(.8))
      d=suppressMessages(ctGenerate(gm,n.subjects = 1,burnin = 3,dtmean = dt))
      if(subi==1) dat=cbind(subi,d) else dat=rbind(dat,cbind(subi,d))
    }
    
    colnames(dat)[1]='id'
    
    cm <- ctModel(LAMBDA=diag(1), type='stanct',
      CINT=matrix('cint'),
      MANIFESTMEANS = matrix(0))
    
    
    cm$pars$indvarying <- FALSE
    cm$pars$indvarying[cm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    dm <- ctModel(LAMBDA=diag(1), type='standt',
      CINT=matrix('cint'),
      MANIFESTMEANS = matrix(0))
    
    dm$pars$indvarying <- FALSE
    dm$pars$indvarying[dm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
 
    for(m in c('cm','dm')){
      argslist <- list(
        ml=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE,
          verbose=0,optimcontrol=list(estonly=FALSE,stochastic=F),savescores = FALSE,nopriors=TRUE)
        # ,mlnl=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(nldynamics=TRUE),
        #   verbose=0,optimcontrol=list(plot=F,estonly=F,stochastic=F),savescores = FALSE,nopriors=TRUE)
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

        if(!interactive()){
          ctKalman(f,plot=TRUE)
        plot(f,wait=FALSE)
        ctModelLatex(f)
        p=ctStanKalman(f,collapsefunc = mean,subjectpars = TRUE)
        }
      }
    }
    
    dtpars=lapply(s, function(argi) {
      print(argi)
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

  }) #end cint heterogeneity
    
    
    test_that("dtVct_noheterogeneity", {
      set.seed(1)
      s=list()
      nsubjects=100
      Tpoints=10
      parsd=0
      parmu= -1.4
      dt=1
      par= (rnorm(nsubjects,parmu,parsd))
      mean(par)
      sd(par)
      
      for(subi in 1:nsubjects){
        gm=ctModel(LAMBDA=diag(1), Tpoints=Tpoints, DRIFT=matrix(-.5),T0MEANS = matrix(4), 
          CINT=matrix(par[subi]),DIFFUSION=matrix(2),
          T0VAR=matrix(2), MANIFESTVAR=matrix(2))
        d=suppressMessages(ctGenerate(gm,n.subjects = 1,burnin = 10,dtmean = dt))
        if(subi==1) dat=cbind(subi,d) else dat=rbind(dat,cbind(subi,d))
      }
      
      colnames(dat)[1]='id'
      
      cm <- ctModel(LAMBDA=diag(1), type='stanct',
        CINT=matrix('cint'),
        MANIFESTMEANS = matrix(0))
      
      cm$pars$indvarying <- FALSE
      
      dm <- ctModel(LAMBDA=diag(1), type='standt',
        CINT=matrix('cint'),
        MANIFESTMEANS = matrix(0))
      
      dm$pars$indvarying <- FALSE
      
      for(m in c('cm','dm')){
        argslist <- list(ml=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(),
          verbose=0,optimcontrol=list(estonly=FALSE,stochastic=F),savescores = F,nopriors=T)
          ,mlnl=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(nldynamics=TRUE,nlmeasurement=TRUE),
            verbose=0,optimcontrol=list(estonly=F,stochastic=F),savescores = F,nopriors=T)
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
      
      
    } #end no heterogeneity
      
      
    )
}
