if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4 & 
    !(.Platform$OS.type=="windows" && R.version$major %in% 4 && as.numeric(R.version$minor) >= 2 &&
    unlist(utils::packageVersion('rstan'))[2] < 25)){
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
    
    cm <- ctModel(LAMBDA=matrix(c('lbystate * eta2 + 1',0,0,1),2,2),  T0MEANS=c('t0m1','t0m2|log1p_exp(param)'),
      T0VAR=matrix(c('t0v11',0,0,'t0v22'),2,2),
      PARS=c('lbystate|log1p_exp(param)','lbystate * eta2 + 1'),type='stanct')
    
    cm$pars$indvarying <- FALSE
    # cm$pars$indvarying[cm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    dm <- ctModel(LAMBDA=matrix(c('lbystate * eta2 + 1',0,0,1),2,2),  T0MEANS=c('t0m1','t0m2|log1p_exp(param)'),
      T0VAR=matrix(c('t0v11',0,0,'t0v22'),2,2),
      PARS=c('lbystate|log1p_exp(param)','lbystate * eta2 + 1'),type='standt')
    
    dm$pars$indvarying <- FALSE
    # dm$pars$indvarying[dm$pars$matrix %in% c('CINT','T0MEANS')] <- TRUE
    
    
    
    
    for(m in c('cm','dm')){
      argslist <- list(
        ml=list(datalong = dat,ctstanmodel = get(m))
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
          # print(paste0(ctpars[i,'matrix'],' ', ctpars[i,'row'],',', ctpars[i,'col'],' ',
          # colnames(ctpars)[ti],' = ', ctpars[i,ti],', ',dtpars[ri,ti]))
          testthat::expect_equivalent(ctpars[i,ti],dtpars[ri,ti],tol=ifelse(ti==4,1e-1,1e-1))
        }
      }
    }
    
    
    ll=unlist(lapply(s, function(argi) lapply(argi, function(m) m$loglik)))
    
    for(dimi in 2:length(ll)){
      expect_equivalent(ll[dimi],ll[dimi-1],tol=1e-3)
    }
    
    #do.call(cbind,dtpars)
    
    
    #check time varying lambda estimation
    sapply(s$ml,function(x) expect_equivalent(lambdafactor,x$popmeans[rownames(x$popmeans) %in% 'lbystate','mean'],tol=1e-1))
    
  }) 
  
  test_that("higherDimNonLinearCompileCheck", {
    test_ <- ctModel(type='stanct',
      n.latent=3, n.manifest=3,
      manifestNames=c("X", "Y", "Z"),
      latentNames = c("X_", "Y_", "Z_"),
      Tpoints = 4,
      time = "time",
      LAMBDA=diag(3),
      TRAITVAR = "auto",
      DRIFT=matrix(c('a11', 'a12','a13',
        '(a + b * Z_)', 'a22', 'a23',
        'a31', 'a32', 'a33'), nrow=3, ncol=3, byrow=TRUE),
      DIFFUSION='auto',
      T0VAR='auto',
      CINT=0,
      T0MEANS = 'auto',
      MANIFESTMEANS = 'auto',
      MANIFESTVAR = 0,
      PARS=c('a', 'b'))
    
    
    gm <- ctModel(LAMBDA=diag(2), #diagonal factor loading, 2 latents 2 observables
      Tpoints = 5,
      DRIFT=matrix(c(-1,.5,0,-1),2,2), #temporal dynamics
      TRAITVAR = diag(.5,2), #stable latent intercept variance (cholesky factor)
      DIFFUSION=diag(2)) #within person covariance 
    
    d <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 100,
      burnin = 20,dtmean = 1))
    
    d<-data.frame(d)
    d$Z <- d$Y1 + rnorm(nrow(d))
    d$X <- d$Y1
    d$Y <- d$Y2
    
    f <- ctStanFit(datalong = d,ctstanmodel = test_)
    testthat::expect_equivalent(class(f),'ctStanFit')
  })
  
}
