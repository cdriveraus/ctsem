if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  
  library(ctsem)
  library(testthat)
  set.seed(1)
  
  context("timevarying")
  
  test_that("varyingLAMBDA", {
    set.seed(1)
    s=list()
    nsubjects=10
    Tpoints=150
    lambdafactor = .3
    dt=1

    for(subi in 1:nsubjects){
      gm=ctModel(LAMBDA=diag(2), Tpoints=Tpoints, DRIFT=diag(-.2,2),T0MEANS = matrix(c(3,2)), 
        DIFFUSION=diag(.5,2),
        T0VAR=diag(2))
      d=suppressMessages(ctGenerate(gm,n.subjects = 1,burnin = 3,wide = FALSE,dtmean = dt))
      d[,'id'] <- subi
      if(subi==1) dat=d else dat=rbind(dat,d)
    }
    dat[,'Y1'] <-  dat[,'Y1'] * (1+ lambdafactor * dat[,'Y2']) #state dependent lambda
    dat[,c('Y1','Y2')] <- dat[,c('Y1','Y2')] + rnorm(nrow(dat)*2,0,.5) #measurement error
    
    colnames(dat)[1]='id'
    
    cm <- ctModel(LAMBDA=matrix(c('1 + PARS[1,1] * state[2]',0,0,1),2,2), PARS=matrix('lbystate'),type='stanct')
    
    cm$pars$indvarying <- FALSE
    # cm$pars$indvarying[cm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    dm<- ctModel(LAMBDA=matrix(c('1 + PARS[1,1] * state[2]',0,0,1),2,2), PARS=matrix('lbystate'),, type='standt')
    
    dm$pars$indvarying <- FALSE
    # dm$pars$indvarying[dm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    
    for(m in c('cm','dm')){
      argslist <- list(
        ml=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(),
        verbose=0,optimcontrol=list(plotsgd=F,estonly=F,stochastic=F),savescores = F,nopriors=F)
        #, mlis=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(Jstep=1e-6), 
        #   verbose=0,optimcontrol=list(plotsgd=F,estonly=F,isloops=1,stochastic=F),savescores = F,nopriors=T)
        # ,mapis=list(datalong = dat,ctstanmodel = get(m),optimize=TRUE, nlcontrol=list(Jstep=1e-6),
        # verbose=0,optimcontrol=list(plotsgd=F,estonly=F,isloops=1,stochastic=F),savescores = F,nopriors=F)
        # ,hmcintoverpop=list(datalong = dat,ctstanmodel = get(m),optimize=F,iter=500,chains=3, nlcontrol=list(Jstep=1e-6), 
        #   verbose=0,optimcontrol=list(plotsgd=F,estonly=F,stochastic=F),savescores = F,nopriors=F,intoverpop=T),
        # hmc=list(datalong = dat,ctstanmodel = get(m),optimize=F,iter=500,chains=3, nlcontrol=list(),
          # verbose=0,optimcontrol=list(plotsgd=F,estonly=F,stochastic=F),savescores = F,nopriors=F,control=list(max_treedepth=7))
      )
      
      
      for(argi in names(argslist)){
        f = do.call(ctStanFit,argslist[[argi]])
        if(is.null(s[[argi]])) s[[argi]] = list()
        s[[argi]][[m]] <- summary(f,parmatrices=TRUE)
      }
    }
    
    
     dtpars=lapply(s, function(argi) {
        ct=argi$cm$parmatrices[grep('dt',rownames(argi$cm$parmatrices)),c('Mean','Sd')]
        dt=argi$dm$parmatrices[rownames(argi$dm$parmatrices) %in% c('DRIFT','CINT','DIFFUSION'),c('Mean','Sd')]
        rownames(dt) <- paste0('dt',rownames(dt))
        list(ct=ct,dt=dt)
      }
      )
      
      for(ri in 1:nrow(dtpars[[1]]$ct)){
        for(ci in 1:ncol(dtpars[[1]]$ct)){
          par=sapply(dtpars, function(x) sapply(x, function(y) y[ri,ci]))
          for(dimi in 2:length(par)){
            expect_equivalent(par[dimi],par[dimi-1],tol=1e-1)
          }
        }}
      
      
      ll=unlist(lapply(s, function(argi) lapply(argi, function(m) m$logprob)))
      
      for(dimi in 2:length(ll)){
        expect_equivalent(ll[dimi],ll[dimi-1],tol=1e-3)
      }
      
      #check time varying lambda estimation
      sapply(s$ml,function(x) expect_equivalent(lambdafactor,x$popmeans[rownames(x$popmeans) %in% 'lbystate','mean'],tol=1e-1))
    
  }) 
}
