if(identical(Sys.getenv("NOT_CRAN"), "true") & .Machine$sizeof.pointer != 4){
  # Sys.setenv(NOT_CRAN='true')

  set.seed(1)
  library(ctsem)
  library(testthat)
  
  context("dtVct_lVnl")
  
  test_that("dtVct_CINTheterogeneity", {
    set.seed(1)
    s=list()
    nsubjects=500
    Tpoints=15
    parsd=0#1.4
    parmu= -3.4
    dt=1
    par= (rnorm(nsubjects,parmu,parsd))
    mean(par)
    sd(par)
    
    for(subi in 1:nsubjects){
      gm=ctModel(LAMBDA=diag(1), Tpoints=Tpoints, DRIFT=matrix(-.3),T0MEANS = matrix(4), 
        CINT=matrix(par[subi]),DIFFUSION=matrix(2),
        T0VAR=matrix(2), MANIFESTVAR=matrix(.3))
      d=suppressMessages(ctGenerate(gm,n.subjects = 1,burnin = 0,dtmean = dt))
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
          verbose=0,optimcontrol=list(estonly=FALSE,stochastic=F),savescores = FALSE,nopriors=TRUE))

      for(argi in names(argslist)){
        f = do.call(ctStanFit,argslist[[argi]])
        if(is.null(s[[argi]])) s[[argi]] = list()
        s[[argi]][[m]] <- summary(f,parmatrices=TRUE)

        # if(!interactive()){
        #   ctKalman(f,plot=TRUE)
        # plot(f,wait=FALSE)
        # ctModelLatex(f)
        # p=ctStanKalman(f,collapsefunc = mean,subjectpars = TRUE)
        # }
      }
    }
    
    ctpars=s[[1]]$cm$parmatrices
    ctpars <- ctpars[!ctpars$matrix %in% c('DRIFT','CINT','DIFFUSIONcov'),]
    dtpars=s[[1]]$dm$parmatrices
    dtpars$matrix[dtpars$matrix %in% 'DRIFT'] <- 'dtDRIFT'
    
    for(ri in 1:nrow(dtpars)){
      i <- which(apply(ctpars,1,function(x) all(x[1:3] == dtpars[ri,1:3])))
      if(length(i)>0){
        for(ti in 4:5){
          print(c(ctpars[i,ti],dtpars[ri,ti]))
        testthat::expect_equivalent(ctpars[i,ti],dtpars[ri,ti],tol=ifelse(ti==4,1e-2,1e-1))
        }
      }
    }
    
   
    
    ll=unlist(lapply(s, function(argi) lapply(argi, function(m) m$loglik)))
    
    for(dimi in 2:length(ll)){
      expect_equivalent(ll[dimi],ll[dimi-1],tol=1e-3)
    }

  }) #end cint heterogeneity
    
    
    test_that("dtVct_noheterogeneity", {
      set.seed(1)
      s=list()
      nsubjects=200
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
          ,mlnl=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(nldynamics=TRUE),
            verbose=0,optimcontrol=list(estonly=F,stochastic=F),savescores = F,nopriors=T)
        )
        
        
        for(argi in names(argslist)){
          f = do.call(ctStanFit,argslist[[argi]])
          if(is.null(s[[argi]])) s[[argi]] = list()
          s[[argi]][[m]] <- summary(f,parmatrices=TRUE)
        }
      }
      ctpars=s[[1]]$cm$parmatrices
      ctpars <- ctpars[!ctpars$matrix %in% c('DRIFT','CINT','DIFFUSIONcov'),]
      dtpars=s[[1]]$dm$parmatrices
      dtpars$matrix[dtpars$matrix %in% 'DRIFT'] <- 'dtDRIFT'
      
      for(ri in 1:nrow(dtpars)){
        i <- which(apply(ctpars,1,function(x) all(x[1:3] == dtpars[ri,1:3])))
        if(length(i)>0){
          for(ti in 4:5){
            # print(c(ctpars[i,ti],dtpars[ri,ti]))
            testthat::expect_equivalent(ctpars[i,ti],dtpars[ri,ti],tol=ifelse(ti==4,1e-1,1e-1))
          }
        }
      }
      
      
      ll=unlist(lapply(s, function(argi) lapply(argi, function(m) m$loglik)))
      
      for(dimi in 2:length(ll)){
        expect_equivalent(ll[dimi],ll[dimi-1],tol=1e-3)
      }
      
      
    } #end no heterogeneity
      
      
    )
}
