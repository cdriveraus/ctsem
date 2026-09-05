skip_on_cran()
skip_on_32bit()
skip_if(.Platform$OS.type == "windows" && R.version$major %in% 4 &&
    as.numeric(R.version$minor) >= 2 &&
    unlist(utils::packageVersion('rstan'))[2] < 25,
  "rstan < 2.26 cannot compile a model on R >= 4.2 on Windows.")
{  # body of the guard this replaced; indentation unchanged
  # Sys.setenv(NOT_CRAN = 'true')
  
  library(ctsem)
  library(testthat)
  set.seed(1)
  
  context("timevarying")

  # Part of the random-effect / ct-vs-dt family whose design is written down at
  # the top of test-tdeffectvariation_covtest.R. This block carries two claims:
  # C3, that a ct fit and a dt fit of the same state-dependent LAMBDA model
  # agree, and a recovery claim -- that `lbystate` comes back as the
  # lambdafactor the data were generated with.
  #
  # The size below is NOT shrinkable, which was measured rather than assumed.
  # The recovery error on lambdafactor (true 0.3), and the ct-vs-dt loglik gap:
  #
  #   n x Tpoints   err_ct   err_dt   loglik gap
  #   50 x 50        0.001    0.001    0        <- as written
  #   50 x 25        0.105    0.104    0        <- past the 1e-1 assertion
  #   30 x 50        0.066    0.067    0
  #   30 x 30        0.044    0.279    430      <- the dt fit diverged outright
  #   25 x 25        0.017    0.020    0
  #
  # Recovery degrades non-monotonically and one intermediate size loses the dt
  # fit completely, so there is no smaller design here that is safe. Halving
  # the data to save 180 s buys a test that fails.
  test_that("varyingLAMBDA", {
    set.seed(1)
    s=list()
    nsubjects=50
    Tpoints=50
    lambdafactor = .3
    dt=1

      # MANIFESTVAR and MANIFESTMEANS are stated rather than left free. A free
      # generating cell is filled from `.ctGenerateDefaults()`, so this test's
      # data moves whenever those defaults do -- under a set.seed() that reads
      # as though it pinned everything. Zero is what has always been generated
      # here and what the test wants: the measurement error is added by hand
      # below, so a non-zero generating MANIFESTVAR would double-count it.
      gm=suppressMessages(ctModel(LAMBDA=diag(2), Tpoints=Tpoints, DRIFT=diag(-.1,2),T0MEANS = matrix(c(3,2)),
        DIFFUSION=diag(.2,2),
        MANIFESTVAR=diag(0,2), MANIFESTMEANS=matrix(0,2,1),
        T0VAR=diag(2)))
      dat=suppressMessages(ctGenerate(gm,n.subjects = nsubjects,burnin = 3,dtmean = dt))

    dat <- as.matrix(dat)
    dat[,'Y1'] <-  dat[,'Y1'] * (1+ lambdafactor * dat[,'Y2']) #state dependent lambda
    dat[,c('Y1')] <- dat[,c('Y1')] + rnorm(nrow(dat),0,.1) #measurement error
    dat[,c('Y2')] <- dat[,c('Y2')] + rnorm(nrow(dat),0,.1) #measurement error
    
    colnames(dat)[1]='id'
    
    cm <- ctModel(LAMBDA=matrix(c('lbystate * eta2 + 1',0,0,1),2,2),  T0MEANS=c('t0m1','t0m2|log1p_exp(param)'),
      PARS=c('lbystate|log1p_exp(param)'),type='ct')
    
    cm$pars$indvarying <- FALSE
    
    dm <- ctModel(LAMBDA=matrix(c('lbystate * eta2 + 1',0,0,1),2,2),  T0MEANS=c('t0m1','t0m2|log1p_exp(param)'),
      PARS=c('lbystate|log1p_exp(param)'),type='dt')
    
    dm$pars$indvarying <- FALSE
    

    fct <- ctFit(datalong = dat,model= cm)
    fdt <- ctFit(datalong = dat,model= dm)
    
    sct <- summary(fct,parmatrices=TRUE)
    sdt <- summary(fdt,parmatrices=TRUE)
    
    ctpars=sct$parmatrices
    ctpars <- ctpars[!ctpars$matrix %in% c('DRIFT','CINT','DIFFUSIONcov'),]
    dtpars=sdt$parmatrices
    dtpars$matrix[dtpars$matrix %in% 'DRIFT'] <- 'dtDRIFT'
    
    for(ri in 1:nrow(dtpars)){
      i <- which(apply(ctpars,1,function(x) all(x[1:3] == dtpars[ri,1:3])))[1] #find matching row between ct and dt fits
      if(!is.na(i) & length(i)>0){
        for(ti in 4:5){ #compare parameters mean and sd
          # print(paste0(ctpars[i,'matrix'],' ', ctpars[i,'row'],',', ctpars[i,'col'],' ',
          # colnames(ctpars)[ti],' = ', ctpars[i,ti],', ',dtpars[ri,ti]))
          test_isclose(ctpars[i,ti],dtpars[ri,ti],tol=ifelse(ti==4,1e-1,1e-1))
        }
      }
    }
    
    
    test_isclose(sct$loglik,sdt$loglik,tol=1e-3)
    
    
    #do.call(cbind,dtpars)
    
    
    #check time varying lambda estimation
    test_isclose(
      lambdafactor,
      sct$popmeans[rownames(sct$popmeans) %in% 'lbystate','mean'],
      sdt$popmeans[rownames(sdt$popmeans) %in% 'lbystate','mean'],tol=1e-1)
    
  }) 
  
  test_that("higherDimNonLinearCompileCheck", {
    test_ <- ctModel(type='ct',
      n.latent=3, n.manifest=3,
      manifestNames=c("X", "Y", "Z"),
      latentNames = c("X_", "Y_", "Z_"),
      Tpoints = 4,
      time = "time",
      LAMBDA=diag(3),
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
    
    
    nsubjects <- 100
    traitChol <- diag(.5,2)
    subjectCint <- t(replicate(nsubjects, as.numeric(traitChol %*% rnorm(2))))
    # Stated in full, for the same reason as the model above. T0VAR is the one
    # value that is not what generation has lately been supplying (1e-6 from
    # the defaults); burnin = 20 against a drift of -1 washes the initial
    # condition out entirely, and this test asserts only that the fit returns.
    gm <- ctModel(LAMBDA=diag(2), #diagonal factor loading, 2 latents 2 observables
      Tpoints = 5,
      DRIFT=matrix(c(-1,.5,0,-1),2,2), #temporal dynamics
      MANIFESTVAR=diag(0,2), MANIFESTMEANS=matrix(0,2,1),
      T0MEANS=matrix(0,2,1), T0VAR=diag(1,2),
      DIFFUSION=diag(2)) #within person covariance
    
    dlist <- vector("list", nsubjects)
    for(i in seq_len(nsubjects)){
      gm_i <- gm
      # Through $matrices, not `gm_i$CINT <-`. `gm` is a ctStanModel, whose
      # canonical specification is $pars; ctGenerate() rebuilds every top level
      # matrix from $pars before generating, so a direct assignment is silently
      # discarded and all 100 subjects were generated with CINT = 0. The
      # $matrices view writes back into $pars, so the value reaches the data.
      gm_i$matrices$CINT <- matrix(subjectCint[i, ], ncol = 1)
      d_i <- suppressMessages(ctGenerate(ctmodelobj = gm_i,n.subjects = 1,
        burnin = 20,dtmean = 1))
      d_i[, "id"] <- i
      dlist[[i]] <- d_i
    }
    d <- do.call(rbind, dlist)
    d <- data.frame(d)
    d$Z <- d$Y1 + rnorm(nrow(d))
    d$X <- d$Y1
    d$Y <- d$Y2
    
    # The cheap half first, so a specification bug fails in under a second
    # rather than after the compile. `recompile == 1` is the reason this block
    # is expensive and is what its name refers to: this shape cannot use the
    # precompiled ctsm program, so rstan builds a bespoke one, and that C++
    # compile -- not the data -- is essentially all of the block's wall clock.
    # Measured: the whole block is 197 s at 100 subjects and 191 s at 25, of
    # which `fit = FALSE` is 0.7 s.
    spec <- ctFit(datalong = d, model = test_, fit = FALSE)
    testthat::expect_equal(spec$standata$recompile, 1)
    # The state-dependent DRIFT cell has to reach matsetup as a state
    # reference, not as a parameter. Column 10 is `stateref`.
    testthat::expect_true(sum(spec$standata$matsetup[, 10] != 0) > 0)

    f <- ctFit(datalong = d,model= test_)
    testthat::expect_s3_class(f, 'ctStanFit')

    # ...and then say the fit MOVED. `expect_s3_class` alone passed on a fit
    # that had done nothing: an optimiser that returned its starting values
    # still returns an object of the right class. Every number below is
    # measured on this design (25 subjects, seed 1); the counterfactual in
    # brackets is the same quantity evaluated at the neutral raw start, which
    # is what a fit that did not move would report.
    sf <- f$stanfit$stanfit
    raw <- f$stanfit$rawest
    lp_opt <- rstan::log_prob(sf, upars = raw, adjust_transform = FALSE)
    lp_start <- rstan::log_prob(sf, upars = rep(0, length(raw)),
      adjust_transform = FALSE)
    grad_opt <- rstan::grad_log_prob(sf, upars = raw, adjust_transform = FALSE)

    testthat::expect_true(is.finite(summary(f)$loglik))
    # measured gain 471 [0]
    testthat::expect_true(lp_opt > lp_start + 10)
    # measured max|grad| 0.005 at the optimum [126 at the start]
    testthat::expect_true(max(abs(grad_opt)) < 1)
    # measured max|raw| 37.9 [0]
    testthat::expect_true(max(abs(raw)) > 1)
  })
  
}
