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
    skip_without_julia()
    nsubjects <- 600
    # The generating covariances, stated as covariances. ctCov() writes the
    # cells that give them under the model's construction; the manifest traits
    # are individually varying MANIFESTMEANS, drawn by ctGenerate from the
    # population RAWPOPVAR states -- on the raw scale, and MANIFESTMEANS is
    # 10 * raw, so the raw covariance is a hundredth of the trait covariance.
    # MANIFESTVAR is zero, as the fitted model below assumes.
    targetManifestCov <- tcrossprod(matrix(c(2,-1,-1, 0,1,1,0,0,2),3,3))
    diffusionCov <- tcrossprod(matrix(c(2,1,1,0,4,-2,0,0,2),3,3))
    gm <- ctModel(LAMBDA = diag(1,3), DRIFT = diag(-1,3),
      T0VAR = ctCov(tcrossprod(matrix(c(5,-5,-5,0,1,-1,0,0,2),3,3))),
      DIFFUSION = ctCov(diffusionCov), MANIFESTVAR = diag(0,3),
      T0MEANS = matrix(0,3,1), CINT = matrix(0,3,1),
      MANIFESTMEANS = matrix(c('mm1','mm2','mm3')), Tpoints = 30)
    gm$pars$indvarying <- gm$pars$param %in% c('mm1','mm2','mm3')
    gm$RAWPOPVAR <- ctCov(targetManifestCov / 100)
    d <- suppressMessages(ctGenerate(gm, n = nsubjects, burnin = 0,
      backend = 'julia'))
    
    m <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),type='ct',
      # MANIFESTMEANS = 0,
      MANIFESTVAR = 0)
    
    f <- ctFit(datalong = d,model= m,priors=T,verbose=0,cores=cores, backend = 'stan')
    
    p <- ctSummaryMatrices(f)
    
    diffcov <- p$DIFFUSIONcov
    diffcor <- cov2cor(p$DIFFUSIONcov)
    ediffcov <- diffusionCov
    ediffcor <- cov2cor(diffusionCov)
    
    s=summary(f)
    s$rawpopcorr
    f$stanfit$transformedparsfull$rawpopcorr[1,,]
    targetManifestCov
    cov2cor(targetManifestCov)
    
    #check diagonal of 1's for corr
    test_isclose(diag(f$stanfit$transformedparsfull$rawpopcorr[1,,]),
      rep(1,nrow(f$stanfit$transformedparsfull$rawpopcorr[1,,])),tol=1e-3)
    
    # The traits the data were generated with. 600 draws leave their own
    # covariance up to about 0.8 from the target (a variance of 6.75 for 6 at
    # this seed), which is not the fit's error, so the fit is compared with
    # the traits' covariance and the traits with the target.
    traits <- stats::cov(as.matrix(attr(d, 'subjectpars')[, c('mm1','mm2','mm3')]))
    expect_lt(max(abs(traits - targetManifestCov)), 1)

    #cov check
    expect_lt(max(abs(f$stanfit$transformedparsfull$popcov[1,4:6,4:6] - traits)), 0.5)

    #cor check
    expect_lt(max(abs(f$stanfit$transformedparsfull$rawpopcorr[1,4:6,4:6] -
      cov2cor(traits))), 0.1)

    # The initial correlations and the diffusion covariance the data were
    # generated with. Until 2026-11-22 the generator read T0VAR's cells through
    # another construction than the model's, and the initial correlations of
    # -0.98, -0.91 and 0.86 were not the data's.
    expect_lt(max(abs(cov2cor(p$T0cov) - cov2cor(tcrossprod(
      matrix(c(5,-5,-5,0,1,-1,0,0,2),3,3))))), 0.05)
    expect_lt(max(abs(diffcov - ediffcov)), 0.5)

    
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
    d1 <- data.frame(ctGenerate(ctmodelobj = gm,n = 1000,burnin = 0))
    
    gm <- ctModel(type='omx',LAMBDA = diag(1,10),DRIFT=diag(-1,10),
      T0VAR=cmat2,
      DIFFUSION=diag(1,10),Tpoints=2)
    d2 <- data.frame(ctGenerate(ctmodelobj = gm,n = 1000,burnin = 0))
    
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
    
    f <- ctFit(datalong = d,model= m,priors=TRUE,cores=cores,verbose=0, backend = 'stan')
    
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

