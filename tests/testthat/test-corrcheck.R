skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged
  library(ctsem)
  library(testthat)
  cores=2
  
  # context("corrCheck")
  
  #anomauth
  test_that("corrCheck", {
    set.seed(1)
    nsubjects <- 600
    manifestTraitChol <- matrix(c(2,-1,-1, 0,1,1,0,0,2),3,3)
    # MANIFESTVAR and T0MEANS are stated rather than left free. A free
    # generating cell is filled from `.ctGenerateDefaults()`, so this test's
    # data moves whenever those defaults do -- under a set.seed() that reads as
    # though it pinned everything. Zero is what has always been generated here
    # and what the fitted model below assumes: it fixes MANIFESTVAR to 0, so a
    # non-zero generating value would be a misspecification.
    gm <- ctModel(type='omx',LAMBDA = diag(1,3),DRIFT=diag(-1,3),
      T0VAR=matrix(c(5,-5,-5,0,1,-1,0,0,2),3,3),
      MANIFESTVAR=diag(0,3), T0MEANS=matrix(0,3,1),
      DIFFUSION=matrix(c(2,1,1,0,4,-2,0,0,2),3,3),Tpoints=30)
    
    subjectManifestMeans <- t(replicate(nsubjects, as.numeric(manifestTraitChol %*% rnorm(3))))
    targetManifestCov <- stats::cov(subjectManifestMeans)
    dlist <- vector("list", nsubjects)
    for(i in seq_len(nsubjects)){
      gm_i <- gm
      gm_i$MANIFESTMEANS <- matrix(subjectManifestMeans[i, ], ncol = 1)
      d_i <- ctGenerate(ctmodelobj = gm_i, n.subjects = 1, burnin = 0)
      d_i[, "id"] <- i
      dlist[[i]] <- d_i
    }
    d <- do.call(rbind, dlist)
    
    m <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),type='ct',
      # MANIFESTMEANS = 0,
      MANIFESTVAR = 0)
    
    f <- ctFit(datalong = d,model= m,priors=T,verbose=0,cores=cores)
    
    p <- ctSummaryMatrices(f)
    
    diffcov <- p$DIFFUSIONcov
    diffcor <- cov2cor(p$DIFFUSIONcov)
    ediffcov <- tcrossprod(gm$DIFFUSION)
    ediffcor <- cov2cor(tcrossprod(gm$DIFFUSION))
    
    s=summary(f)
    s$rawpopcorr
    f$stanfit$transformedparsfull$rawpopcorr[1,,]
    targetManifestCov
    cov2cor(targetManifestCov)
    
    #check diagonal of 1's for corr
    test_isclose(diag(f$stanfit$transformedparsfull$rawpopcorr[1,,]),
      rep(1,nrow(f$stanfit$transformedparsfull$rawpopcorr[1,,])),tol=1e-3)
    
    #cov check
    test_isclose(f$stanfit$transformedparsfull$popcov[1,4:6,4:6],
      targetManifestCov,tol=1)
    
    #cor check
    test_isclose(f$stanfit$transformedparsfull$rawpopcorr[1,4:6,4:6],
      cov2cor(targetManifestCov),tol=.1)
    
    
  })
  
  if(FALSE) test_that("corrCheckHighDim", {
    set.seed(1)
    
    cmat <- diag(.5,10) + 1 
    cmat=t(chol(cmat %*% t(cmat)))
    
    cmat2 <- diag(.5,10) + 1
    cmat2[5:10,] <- cmat2[5:10,] * -1
    cmat2=t(chol(cmat2 %*% t(cmat2)))
    cov2cor(tcrossprod(cmat2))
    
    gm <- ctModel(type='omx',LAMBDA = diag(1,10),DRIFT=diag(-1,10),
      T0VAR=cmat,
      DIFFUSION=diag(1,10),Tpoints=2)
    d1 <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1000,burnin = 0))
    
    gm <- ctModel(type='omx',LAMBDA = diag(1,10),DRIFT=diag(-1,10),
      T0VAR=cmat2,
      DIFFUSION=diag(1,10),Tpoints=2)
    d2 <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1000,burnin = 0))
    
    d2$id <- d2$id + 2000
    d <- rbind(d1,d2)
    d$TI1 <- 0
    d$TI1[d$id > 2000] <- 1
  
    
    
    
    m <- ctModel(LAMBDA = diag(1,10),DRIFT=diag(-1,10),type='dt',
      # MANIFESTMEANS = 0,
      DIFFUSION=diag(1,10),T0MEANS=0,
      TIpredNames = 'TI1',
      MANIFESTMEANS=0,
      MANIFESTVAR = 0)
    
    f <- ctFit(datalong = d,model= m,priors=TRUE,cores=cores,verbose=0)
    
    p <- ctSummaryMatrices(f)
    
    t0cov <- p$T0cov
    t0cor <- cov2cor(t0cov)
    et0cov <- tcrossprod(gm$T0VAR)
    et0cor <- cov2cor(tcrossprod(gm$T0VAR))
    
    f$stanfit$rawposterior=f$stanfit$rawposterior[1:5,]
    
    e=ctExtract(f,subjectMatrices = T,cores=1)
    cov2cor(e$subj_T0cov[1,1,,])
    cov2cor(e$subj_T0cov[1,301,,])
    
    
    #cov check
    test_isclose(cov2cor(e$subj_T0cov[1,1001,,]),
      cov2cor(tcrossprod(cmat2)),tol=.005)
    
    #cor check
    test_isclose(cov2cor(e$subj_T0cov[1,1,,]),
      cov2cor(tcrossprod(cmat)),tol=.005)
    
    
  })
}

