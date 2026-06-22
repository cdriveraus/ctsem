suppressWarnings(suppressPackageStartupMessages(library(ctsem)))

test_that("ctEmpiricalBayesFit summary adjusts transforms from raw point estimates", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    DRIFT=matrix('drift|param|FALSE',1,1),
    DIFFUSION=matrix(0,1,1),
    CINT=matrix(0,1,1),
    T0MEANS=matrix(0,1,1),
    T0VAR=matrix(0,1,1),
    LAMBDA=matrix(1,1,1),
    MANIFESTMEANS=matrix(0,1,1),
    MANIFESTVAR=matrix('merror',1,1),
    silent=TRUE)
  model$pars$indvarying <- FALSE
  subjectmodel <- model
  subjectmodel$pars$indvarying <- FALSE
  
  fakefit <- function(rawest){
    fit <- list(stanfit=list(rawest=rawest))
    class(fit) <- 'ctStanFit'
    fit
  }
  
  eb <- list(
    subjects=c(1,2,3),
    initialfits=list(
      '1'=fakefit(c(-1, 0)),
      '2'=fakefit(c(0, 1)),
      '3'=fakefit(c(1, 2))),
    fits=list(
      '1'=fakefit(c(9, 9)),
      '2'=fakefit(c(10, 10)),
      '3'=fakefit(c(11, 11))),
    parnames=c('drift','merror'),
    ebUse='rawest',
    ctstanmodel=model,
    subjectmodel=subjectmodel)
  class(eb) <- 'ctEmpiricalBayesFit'
  
  s <- summary(eb, use='rawest', sdscale='unit', digits=6)
  
  expect_s3_class(s, 'summary.ctEmpiricalBayesFit')
  expect_equal(s$initialrawstats$mean, c(0, 1))
  expect_equal(s$rawstats$mean, c(10, 10))
  expect_equal(s$rawstats$sd, c(1, 1))
  expect_equal(s$adjustedmodel$pars$sdscale[
    s$adjustedmodel$pars$param %in% c('drift','merror')], c(1, 1))
  expect_false(any(s$adjustedmodel$pars$indvarying))
  expect_false(is.na(suppressWarnings(as.numeric(
    s$adjustedmodel$pars$transform[
      s$adjustedmodel$pars$param %in% 'merror']))))
})

test_that("ctEmpiricalBayesFit summary can use raw empirical SDs for sdscale", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    DRIFT=matrix('drift',1,1),
    DIFFUSION=matrix(0,1,1),
    CINT=matrix(0,1,1),
    T0MEANS=matrix(0,1,1),
    T0VAR=matrix(0,1,1),
    LAMBDA=matrix(1,1,1),
    MANIFESTMEANS=matrix(0,1,1),
    MANIFESTVAR=matrix('merror',1,1),
    silent=TRUE)
  model$pars$indvarying <- TRUE
  subjectmodel <- model
  subjectmodel$pars$indvarying <- FALSE
  
  fakefit <- function(rawest){
    fit <- list(stanfit=list(rawest=rawest))
    class(fit) <- 'ctStanFit'
    fit
  }
  
  eb <- list(
    subjects=c(1,2,3),
    initialfits=list(
      '1'=fakefit(c(0, 0)),
      '2'=fakefit(c(0, 2)),
      '3'=fakefit(c(0, 4))),
    fits=list(
      '1'=fakefit(c(0, 0)),
      '2'=fakefit(c(0, 2)),
      '3'=fakefit(c(0, 4))),
    parnames=c('drift','merror'),
    ebUse='rawest',
    ctstanmodel=model,
    subjectmodel=subjectmodel)
  class(eb) <- 'ctEmpiricalBayesFit'
  
  s <- summary(eb, use='rawest', sdscale='rawsd', minsd=.001, digits=6)
  
  expect_equal(s$adjustedmodel$pars$sdscale[
    s$adjustedmodel$pars$param %in% 'drift'], .001)
  expect_equal(s$adjustedmodel$pars$sdscale[
    s$adjustedmodel$pars$param %in% 'merror'], 2)
})

test_that("ctEmpiricalBayesFit EB adjustment keeps known transforms numeric", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    DRIFT=matrix('drift|-log1p_exp(-param)|FALSE',1,1),
    DIFFUSION=matrix('diffusion|log1p_exp(param)|FALSE',1,1),
    CINT=matrix(0,1,1),
    T0MEANS=matrix('t0m|param|FALSE',1,1),
    T0VAR=matrix(0,1,1),
    LAMBDA=matrix(1,1,1),
    MANIFESTMEANS=matrix(0,1,1),
    MANIFESTVAR=matrix('merror|log1p_exp(param)|FALSE',1,1),
    silent=TRUE)
  rawstats <- data.frame(
    param=c('t0m','drift','diffusion','merror'),
    mean=c(.2, -.4, .1, .3),
    sd=c(.7, .5, .8, .6))
  
  adjusted <- ctsem:::ctEBadjustModel(model, rawstats)
  fitsetup <- ctsem:::ctModelTransformsToNum(adjusted)
  rows <- fitsetup$pars$param %in% rawstats$param
  
  expect_true(all(!is.na(suppressWarnings(as.numeric(
    fitsetup$pars$transform[rows])))))
  expect_false(any(suppressWarnings(as.numeric(
    fitsetup$pars$transform[rows])) < -10))
  expect_equal(fitsetup$pars$meanscale[
    fitsetup$pars$param %in% 'diffusion'], .8)
  expect_equal(fitsetup$pars$inneroffset[
    fitsetup$pars$param %in% 'diffusion'], .1)
})

test_that("ctEmpiricalBayesFit robust raw handling winsorizes or removes extremes", {
  raw <- matrix(c(-100, 0, 100, 0, 2, 4), ncol=2)
  colnames(raw) <- c('drift','merror')
  
  winsorized <- ctsem:::ctEBrobustRaw(raw,
    outlierMAD=NULL,
    outlierQuantiles=c(.25,.75),
    winsorize=TRUE)
  expect_equal(range(winsorized$raw[, 'drift']), c(-50, 50))
  expect_equal(winsorized$report$nchanged[1], 2)
  
  removed <- ctsem:::ctEBrobustRaw(raw,
    outlierMAD=NULL,
    outlierQuantiles=c(.25,.75),
    winsorize=FALSE)
  expect_equal(sum(is.na(removed$raw[, 'drift'])), 2)
})

test_that("ctEmpiricalBayesFit maps later EB raw estimates back to original raw scale", {
  raw <- matrix(c(-1, 0, 1, 0, 1, 2), ncol=2)
  colnames(raw) <- c('drift','merror')
  rawmap <- data.frame(param=c('drift','merror'), mean=c(10, 20), sd=c(2, 3))
  
  mapped <- ctsem:::ctEBmapRaw(raw, rawmap)
  
  expect_equal(mapped[, 'drift'], c(8, 10, 12))
  expect_equal(mapped[, 'merror'], c(20, 23, 26))
})

test_that("ctEmpiricalBayesFit summary uses final pass map and prior stats", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    DRIFT=matrix('drift|param|FALSE',1,1),
    DIFFUSION=matrix(0,1,1),
    CINT=matrix(0,1,1),
    T0MEANS=matrix(0,1,1),
    T0VAR=matrix(0,1,1),
    LAMBDA=matrix(1,1,1),
    MANIFESTMEANS=matrix(0,1,1),
    MANIFESTVAR=matrix('merror|param|FALSE',1,1),
    silent=TRUE)
  subjectmodel <- model
  subjectmodel$pars$indvarying <- FALSE
  
  fakefit <- function(rawest){
    fit <- list(stanfit=list(rawest=rawest))
    class(fit) <- 'ctStanFit'
    fit
  }
  
  eb <- list(
    subjects=c(1,2,3),
    initialfits=list(
      '1'=fakefit(c(0, 0)),
      '2'=fakefit(c(0, 0)),
      '3'=fakefit(c(0, 0))),
    fits=list(
      '1'=fakefit(c(-1, 0)),
      '2'=fakefit(c(0, 1)),
      '3'=fakefit(c(1, 2))),
    parnames=c('drift','merror'),
    ebUse='rawest',
    ctstanmodel=model,
    subjectmodel=subjectmodel,
    passrawstats=list(
      data.frame(param=c('drift','merror'), mean=c(0, 0), sd=c(1, 1)),
      data.frame(param=c('drift','merror'), mean=c(10, 20), sd=c(2, 3))),
    passrawmaps=list(
      data.frame(param=c('drift','merror'), mean=c(0, 0), sd=c(1, 1)),
      data.frame(param=c('drift','merror'), mean=c(0, 0), sd=c(1, 1)),
      data.frame(param=c('drift','merror'), mean=c(10, 20), sd=c(2, 3))))
  class(eb) <- 'ctEmpiricalBayesFit'
  
  s <- summary(eb, use='rawest', digits=6)
  
  expect_equal(unname(s$originalraw[, 'drift']), c(8, 10, 12))
  expect_equal(unname(s$originalraw[, 'merror']), c(20, 23, 26))
  expect_equal(s$rawstats$mean, c(10, 23))
  expect_equal(s$adjustedmodel$empiricalbayes$rawstats$mean, c(10, 20))
  expect_equal(s$adjustedmodel$empiricalbayes$rawstats$sd, c(2, 3))
})

test_that("ctEmpiricalBayesFit optimization defaults avoid stochastic first pass hessian", {
  args <- ctsem:::ctEBfitArgsOptimDefaults(list(), firstpass=FALSE)
  expect_false(args$optimcontrol$stochastic)
  expect_null(args$optimcontrol$estonly)
  
  firstargs <- ctsem:::ctEBfitArgsOptimDefaults(list(), firstpass=TRUE)
  expect_false(firstargs$optimcontrol$stochastic)
  expect_true(firstargs$optimcontrol$estonly)
  
  userargs <- ctsem:::ctEBfitArgsOptimDefaults(
    list(optimcontrol=list(stochastic=TRUE)),
    firstpass=TRUE)
  expect_true(userargs$optimcontrol$stochastic)
  expect_true(userargs$optimcontrol$estonly)
})

test_that("ctEmpiricalBayesFit rejects TI predictor models", {
  model <- ctModel(type='ct',
    n.latent=1, latentNames='eta1',
    n.manifest=1, manifestNames='Y1',
    TIpredNames='x',
    LAMBDA=matrix(1,1,1),
    silent=TRUE)
  dat <- data.frame(id=c(1,1,2,2), time=c(0,1,0,1), Y1=rnorm(4), x=0)
  
  expect_error(ctEmpiricalBayesFit(dat, model), 'Time independent predictors')
})
