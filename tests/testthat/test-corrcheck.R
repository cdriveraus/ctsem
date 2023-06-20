if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  cores=2
  
  # context("corrCheck")
  
  #anomauth
  test_that("corrCheck", {
    set.seed(1)
    gm <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),
      T0VAR=matrix(c(5,-5,-5,0,1,-1,0,0,2),3,3),
      MANIFESTTRAITVAR=matrix(c(2,-1,-1, 0,1,1,0,0,2),3,3),
      DIFFUSION=matrix(c(2,1,1,0,4,-2,0,0,2),3,3),Tpoints=30)
    d <- ctGenerate(ctmodelobj = gm,n.subjects = 600,burnin = 0)
    
    m <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),type='stanct',
      # MANIFESTMEANS = 0,
      MANIFESTVAR = 0)
    
    f <- ctStanFit(datalong = d,ctstanmodel = m,priors=T,verbose=0,cores=cores)
    
    p <- ctStanContinuousPars(f)
    
    diffcov <- p$DIFFUSIONcov
    diffcor <- cov2cor(p$DIFFUSIONcov)
    ediffcov <- tcrossprod(gm$DIFFUSION)
    ediffcor <- cov2cor(tcrossprod(gm$DIFFUSION))
    
    s=summary(f)
    s$rawpopcorr
    f$stanfit$transformedparsfull$rawpopcorr[1,,]
    tcrossprod(gm$MANIFESTTRAITVAR)
    cov2cor(tcrossprod(gm$MANIFESTTRAITVAR))
    
    #check diagonal of 1's for corr
    testthat::expect_equivalent(diag(f$stanfit$transformedparsfull$rawpopcorr[1,,]),
      rep(1,nrow(f$stanfit$transformedparsfull$rawpopcorr[1,,])),tol=1e-5)
    
    #cov check
    testthat::expect_equivalent(f$stanfit$transformedparsfull$popcov[1,4:6,4:6],
      tcrossprod(gm$MANIFESTTRAITVAR),tol=.5)
    
    #cor check
    testthat::expect_equivalent(f$stanfit$transformedparsfull$rawpopcorr[1,4:6,4:6],
      cov2cor(tcrossprod(gm$MANIFESTTRAITVAR)),tol=1e-1)
    
    
  })
  
  if(FALSE) test_that("corrCheckHighDim", {
    set.seed(1)
    
    cmat <- diag(.5,10) + 1 
    cmat=t(chol(cmat %*% t(cmat)))
    
    cmat2 <- diag(.5,10) + 1
    cmat2[5:10,] <- cmat2[5:10,] * -1
    cmat2=t(chol(cmat2 %*% t(cmat2)))
    cov2cor(tcrossprod(cmat2))
    
    gm <- ctModel(LAMBDA = diag(1,10),DRIFT=diag(-1,10),
      T0VAR=cmat,
      DIFFUSION=diag(1,10),Tpoints=2)
    d1 <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1000,burnin = 0))
    
    gm <- ctModel(LAMBDA = diag(1,10),DRIFT=diag(-1,10),
      T0VAR=cmat2,
      DIFFUSION=diag(1,10),Tpoints=2)
    d2 <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1000,burnin = 0))
    
    d2$id <- d2$id + 2000
    d <- rbind(d1,d2)
    d$TI1 <- 0
    d$TI1[d$id > 2000] <- 1
  
    
    
    
    m <- ctModel(LAMBDA = diag(1,10),DRIFT=diag(-1,10),type='standt',
      # MANIFESTMEANS = 0,
      DIFFUSION=diag(1,10),T0MEANS=0,
      TIpredNames = 'TI1',
      MANIFESTMEANS=0,
      MANIFESTVAR = 0)
    
    f <- ctStanFit(datalong = d,ctstanmodel = m,priors=TRUE,cores=cores,verbose=0)
    
    p <- ctStanContinuousPars(f)
    
    t0cov <- p$T0cov
    t0cor <- cov2cor(t0cov)
    et0cov <- tcrossprod(gm$T0VAR)
    et0cor <- cov2cor(tcrossprod(gm$T0VAR))
    
    f$stanfit$rawposterior=f$stanfit$rawposterior[1:5,]
    
    e=ctExtract(f,subjectMatrices = T,cores=1)
    cov2cor(e$subj_T0cov[1,1,,])
    cov2cor(e$subj_T0cov[1,301,,])
    
    
    #cov check
    testthat::expect_equivalent(cov2cor(e$subj_T0cov[1,1001,,]),
      cov2cor(tcrossprod(cmat2)),tol=.005)
    
    #cor check
    testthat::expect_equivalent(cov2cor(e$subj_T0cov[1,1,,]),
      cov2cor(tcrossprod(cmat)),tol=.005)
    
    
  })
}

