library(ctsem)
library(testthat)

test_that("ctOptimComputeUncertainty recovers quadratic covariance", {
  A <- diag(c(2, 4))
  lpg <- function(p){
    out <- -0.5 * drop(t(p) %*% A %*% p)
    attr(out, 'gradient') <- -drop(A %*% p)
    out
  }
  
  hess <- ctsem:::ctOptimComputeUncertainty(c(0, 0), list(nsubjects=1),
    sm=NULL, lpgFunc=lpg, uncertainty='hessian', finishsamples=20,
    verbose=0)
  expect_equal(diag(hess$cov), c(.5, .25), tolerance=.01)
  
  surrogate <- ctsem:::ctOptimComputeUncertainty(c(0, 0), list(nsubjects=1),
    sm=NULL, lpgFunc=lpg, uncertainty='surrogate', finishsamples=20,
    control=list(initialCov=diag(2), surrogateNpoints=20), verbose=0)
  expect_equal(diag(surrogate$cov), c(.5, .25), tolerance=.01)
  
  surrogate_from_hessian <- ctsem:::ctOptimComputeUncertainty(c(0, 0),
    list(nsubjects=1), sm=NULL, lpgFunc=lpg, uncertainty='surrogate',
    finishsamples=20, control=list(surrogateNpoints=20), verbose=0)
  expect_equal(diag(surrogate_from_hessian$cov), c(.5, .25),
    tolerance=.01)
  
  is_proposal <- ctsem:::ctOptimComputeUncertainty(c(0, 0),
    list(nsubjects=1), sm=NULL, lpgFunc=lpg, uncertainty='is',
    finishsamples=20, verbose=0)
  expect_equal(diag(is_proposal$cov), c(.5, .25), tolerance=.01)
  expect_equal(is_proposal$method, 'is')
})

test_that("ctOptimNormalDraws returns requested dimensions", {
  draws <- ctsem:::ctOptimNormalDraws(c(1, 2), diag(c(.2, .3)), n=12)
  expect_equal(dim(draws), c(12, 2))
  expect_true(all(is.finite(draws)))
})

test_that("Hessian covariance reports numerical repairs", {
  hess <- -diag(c(1, 0))
  expect_warning(
    cov <- ctsem:::ctOptimCovFromHessian(hess, ridge=1e-6,
      context='test Hessian'),
    'required numerical repair')
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_true(diagnostics$infoRidgeApplied)
  expect_equal(diagnostics$ridge, 1e-6)
})

test_that("Hessian covariance tries raw inversion before repair", {
  expect_silent(
    cov <- ctsem:::ctOptimCovFromHessian(-diag(2), context='clean Hessian'))
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'solve')
  expect_true(diagnostics$rawSolveSucceeded)
  expect_true(diagnostics$rawCholSucceeded)
  expect_false(diagnostics$infoRidgeApplied)
  expect_false(diagnostics$usedGinv)
  
  hess <- diag(c(-1, 1))
  expect_warning(
    cov <- ctsem:::ctOptimCovFromHessian(hess, context='indefinite Hessian'),
    'nearPD')
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'nearPD_cov')
  expect_true(diagnostics$rawSolveSucceeded)
  expect_false(diagnostics$rawCholSucceeded)
  expect_true(diagnostics$usedNearPD)
  expect_false(diagnostics$infoRidgeApplied)
})

test_that("Hessian processing reports one-sided and weak curvature parameters", {
  h1 <- -diag(c(1, 2, .Machine$double.eps))
  h2 <- h1
  h2[2, 2] <- NA_real_
  matsetup <- data.frame(param=1:3, when=0, parname=paste0('p', 1:3))
  expect_message(
    out <- ctsem:::processHessianMatrices(h1, h2, verbose=0,
      matsetup=matsetup),
    'One sided Hessian used for params: p2')
  expect_message(
    ctsem:::processHessianMatrices(h1, h2, verbose=0,
      matsetup=matsetup),
    'may.*not identified: p3')
  expect_equal(out$onesided, 2)
  expect_equal(out$probpars, 3)
})

test_that("surrogate design scales with dimension and filters lp outliers", {
  set.seed(1)
  p <- 12
  lpg <- function(x){
    out <- -0.5 * sum(x^2)
    attr(out, 'gradient') <- -x
    out
  }
  
  surrogate <- ctsem:::ctOptimSurrogateHessian(rep(0, p), lpgFunc=lpg,
    cov=diag(p), npoints=NULL, scale=.5, verbose=0)
  
  expect_equal(nrow(surrogate$design), max(4 * p, 50))
  expect_true(all(surrogate$drops >= surrogate$dropRange[1]))
  expect_true(all(surrogate$drops <= surrogate$dropRange[2]))
  expect_equal(diag(surrogate$hessian), rep(-1, p), tolerance=.01)
})

test_that("surrogate whitening handles correlated proposal covariance", {
  set.seed(2)
  A <- matrix(c(2, .4, .4, 1), 2, 2)
  lpg <- function(x){
    out <- -0.5 * drop(t(x) %*% A %*% x)
    attr(out, 'gradient') <- -drop(A %*% x)
    out
  }
  propcov <- matrix(c(.8, .3, .3, .6), 2, 2)
  surrogate <- ctsem:::ctOptimSurrogateHessian(c(0, 0), lpgFunc=lpg,
    cov=propcov, npoints=60, scale=.5, verbose=0)
  expect_equal(surrogate$hessian, -A, tolerance=.03)
})

test_that("surrogate profiles all fitted curvature directions", {
  set.seed(4)
  lpg <- function(x){
    theta <- log1p(exp(x[1] - 8))
    out <- -0.5 * (theta / .5)^2 - 0.5 * x[2]^2
    grad <- c(-(theta / .25) * plogis(x[1] - 8), -x[2])
    attr(out, 'gradient') <- grad
    out
  }
  surrogate <- ctsem:::ctOptimSurrogateHessian(c(0, 0), lpgFunc=lpg,
    cov=diag(c(100, 1)), npoints=40, scale=.5,
    profile=TRUE, verbose=0)
  cov <- ctsem:::ctOptimCovFromHessian(surrogate$hessian)
  expect_equal(surrogate$profile$nProfiled, 2)
  expect_gt(surrogate$profile$nAdjusted, 0)
  expect_true(any(surrogate$profile$profiles$reached))
  expect_true(all(is.finite(cov)))
})

test_that("surrogate profiling uses drop magnitude for flat directions", {
  lpg <- function(x){
    out <- -0.5 * .001 * x[1]^2
    attr(out, 'gradient') <- -.001 * x
    out
  }
  prof <- ctsem:::ctOptimSurrogateProfileDirections(est=0, lpgFunc=lpg,
    cholcov=matrix(1), directions=matrix(1), targetDrop=2,
    maxStep=64, verbose=0)
  expect_true(all(prof$reached))
  expect_equal(prof$step, rep(sqrt(2 * 2 / .001), 2), tolerance=1)
})

test_that("surrogate profiling expands to surrogate-implied flat target", {
  lpg <- function(x){
    out <- -0.5 * .00025 * x[1]^2
    attr(out, 'gradient') <- -.00025 * x
    out
  }
  profiled <- ctsem:::ctOptimSurrogateProfileCurvature(
    hessWhite=matrix(-.00025), est=0, lpgFunc=lpg, cholcov=matrix(1),
    targetDrop=2, maxStep=64, verbose=0)
  expect_true(all(profiled$profiles$reached))
  expect_gt(max(profiled$profiles$step), 64)
})

test_that("surrogate profiling keeps expanding when observed profile is flatter", {
  lpg <- function(x){
    out <- -0.5 * .00001 * x[1]^2
    attr(out, 'gradient') <- -.00001 * x
    out
  }
  profiled <- ctsem:::ctOptimSurrogateProfileCurvature(
    hessWhite=matrix(-.00025), est=0, lpgFunc=lpg, cholcov=matrix(1),
    targetDrop=2, maxStep=64, verbose=0)
  expect_true(all(profiled$profiles$reached))
  expect_gt(max(profiled$profiles$expansions), 0)
  expect_gt(max(profiled$profiles$step), 250)
})

test_that("optimized uncertainty API uses explicit method and draw names", {
  fit_args <- names(formals(ctsem:::stanoptimis))
  update_args <- names(formals(ctOptimUncertainty))
  
  expect_true(all(c('uncertainty', 'uncertaintyDraws') %in% fit_args))
  expect_false(any(c('is', 'isESS', 'isitersize') %in% fit_args))
  expect_false('bootstrapUncertainty' %in% fit_args)
  expect_true(all(c('uncertainty', 'draws') %in% update_args))
  expect_false('sampleMethod' %in% update_args)
  expect_true('is' %in% eval(formals(ctsem:::stanoptimis)$uncertainty))
  expect_true('is' %in% eval(formals(ctOptimUncertainty)$uncertainty))
  expect_true('opg' %in% eval(formals(ctOptimUncertainty)$uncertainty))
  expect_true('fullbootstrap' %in% eval(formals(ctOptimUncertainty)$uncertainty))
  expect_false('score' %in% eval(formals(ctOptimUncertainty)$uncertainty))
})

test_that("full bootstrap standata resamples and reindexes subjects", {
  standata <- list(
    subject = as.integer(c(1, 1, 2, 2, 3)),
    time = c(0, 1, 0, 1, 0),
    dokalmanrows = as.integer(c(1, 1, 1, 1, 1)),
    nobs_y = as.integer(c(1, 1, 1, 1, 1)),
    ncont_y = as.integer(c(1, 1, 1, 1, 1)),
    nbinary_y = as.integer(c(0, 0, 0, 0, 0)),
    Y = matrix(seq_len(5), ncol=1),
    tdpreds = matrix(seq_len(5), ncol=1),
    whichobs_y = matrix(1L, nrow=5, ncol=1),
    whichbinary_y = matrix(0L, nrow=5, ncol=1),
    whichcont_y = matrix(1L, nrow=5, ncol=1),
    ntipred = 1L,
    tipredsdata = matrix(c(10, 20, 30), ncol=1),
    idmap = data.frame(original=letters[1:3], new=1:3),
    ndatapoints = 5L,
    nsubjects = 3L
  )
  
  boot <- ctsem:::ctOptimBootstrapStandata(standata, c(2, 2, 1))
  
  expect_equal(boot$nsubjects, 3L)
  expect_equal(boot$ndatapoints, 6L)
  expect_equal(as.integer(boot$subject), c(1, 1, 2, 2, 3, 3))
  expect_equal(as.numeric(boot$Y[,1]), c(3, 4, 3, 4, 1, 2))
  expect_equal(as.numeric(boot$tipredsdata[,1]), c(20, 20, 10))
  expect_equal(boot$idmap$new, 1:3)
})

test_that("uncertainty data checks report small sample limitations", {
  expect_error(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=1L, ndatapoints=1L),
      uncertainty='bootstrap', finishsamples=20),
    'at least two'
  )
  
  expect_warning(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=3L, ndatapoints=12L),
      uncertainty='sandwich', finishsamples=20, npars=2L),
    'fewer than ten'
  )
  
  expect_warning(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=1L, ndatapoints=12L),
      uncertainty='sandwich', finishsamples=20, npars=2L),
    'case-level score contributions|single-subject'
  )
  
  expect_warning(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=12L, ndatapoints=60L),
      uncertainty='opg', finishsamples=20, npars=12L),
    'rank limited'
  )
  
  expect_error(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=1L, ndatapoints=12L),
      uncertainty='fullbootstrap', finishsamples=20),
    'at least two subjects'
  )
  
  expect_error(
    ctsem:::ctOptimCheckUncertaintyData(
      list(nsubjects=8L, ndatapoints=40L),
      uncertainty='fullbootstrap', finishsamples=1),
    'at least two samples'
  )
})

# uncertainty='stored' is the cheap redraw: same covariance, new draws, no model
# evaluations. It replaces `ctFitAddSamples()`, which did this on stan alone and
# appended rather than replaced.

test_that("uncertainty='stored' redraws from the fit's covariance and computes nothing", {
  skip_if_not_installed('rstan')
  fit <- ctstantestfit
  before <- fit$stanfit$cov

  set.seed(41)
  redrawn <- suppressMessages(ctOptimUncertainty(fit, uncertainty='stored',
    finishsamples=37, cores=1))

  expect_equal(nrow(redrawn$stanfit$rawposterior), 37L)
  # The covariance, the Hessian and the recorded method survive: this path is
  # about the draws and nothing else.
  expect_equal(unname(redrawn$stanfit$cov), unname(before), tolerance=0)
  expect_equal(redrawn$stanfit$uncertainty$hessian, fit$stanfit$uncertainty$hessian,
    tolerance=0)
  expect_identical(redrawn$stanfit$uncertainty$settings$method, 'hessian')
  expect_true(isTRUE(redrawn$stanfit$uncertainty$settings$redrawn))
  expect_identical(redrawn$stanfit$uncertainty$settings$finishsamples, 37)

  # The draws are exactly `ctOptimNormalDraws()` on the stored covariance --
  # which is the whole claim, and it also pins what a given seed produces.
  set.seed(41)
  direct <- ctsem:::ctOptimNormalDraws(fit$stanfit$rawest, before, 37)
  expect_equal(unname(redrawn$stanfit$rawposterior), unname(direct), tolerance=0)

  # No log-probability evaluations. Counted rather than timed: a wall clock on a
  # shared machine says nothing, and the count is the property that matters.
  calls <- 0L
  trace(ctsem:::ctOptimFitLpgFunc, tracer=function() calls <<- calls + 1L,
    where=asNamespace('ctsem'), print=FALSE)
  on.exit(untrace(ctsem:::ctOptimFitLpgFunc, where=asNamespace('ctsem')), add=TRUE)
  suppressMessages(ctOptimUncertainty(fit, uncertainty='stored', finishsamples=5,
    cores=1))
  expect_identical(calls, 0L)
})

test_that("uncertainty='stored' refuses a fit with no covariance and warns on non-normal draws", {
  skip_if_not_installed('rstan')
  nocov <- ctstantestfit
  nocov$stanfit$cov <- NULL
  expect_error(ctOptimUncertainty(nocov, uncertainty='stored'),
    'no usable one')

  # Redrawing an importance-sampled or bootstrapped posterior gives normal
  # draws from that covariance, which is not the same distribution. Said out
  # loud, because `ctFitAddSamples()` used to mix the two silently.
  wasimis <- ctstantestfit
  wasimis$stanfit$uncertainty$settings$draws <- 'imis'
  wasimis$stanfit$uncertainty$settings$method <- 'is'
  expect_warning(suppressMessages(ctOptimUncertainty(wasimis, uncertainty='stored',
    finishsamples=5, cores=1)), "came from 'imis'")
})

test_that("ctFitAddSamples is deprecated and its draws have not moved", {
  skip_if_not_installed('rstan')
  fit <- ctstantestfit

  # Verbatim reimplementation of the pre-deprecation body. If the function is
  # ever routed through ctOptimUncertainty() this fails, which is the point:
  # `ctOptimNormalDraws()` consumes the same normals in a different order (one
  # `rnorm(n*npar)` filled by column against one `rnorm(npar)` per sample), so
  # every number a given seed used to produce would change.
  set.seed(707)
  mchol <- t(chol(fit$stanfit$cov))
  reference <- matrix(unlist(lapply(1:6, function(x){
    fit$stanfit$rawest + mchol %*% t(matrix(rnorm(length(fit$stanfit$rawest)), nrow=1))
  })), byrow=TRUE, ncol=length(fit$stanfit$rawest))

  set.seed(707)
  added <- suppressWarnings(suppressMessages(ctFitAddSamples(fit, nsamples=6, cores=1)))
  appended <- added$stanfit$rawposterior[-seq_len(nrow(fit$stanfit$rawposterior)), ,
    drop=FALSE]
  expect_equal(unname(appended), unname(reference), tolerance=0)
  # And it still appends rather than replaces.
  expect_equal(unname(added$stanfit$rawposterior[seq_len(nrow(fit$stanfit$rawposterior)), ]),
    unname(fit$stanfit$rawposterior), tolerance=0)

  expect_warning(suppressMessages(ctFitAddSamples(fit, nsamples=2, cores=1)),
    'deprecated')
  expect_warning(suppressMessages(ctAddSamples(fit, nsamples=2, cores=1)),
    'deprecated')
})
