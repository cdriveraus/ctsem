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
})

test_that("ctOptimNormalDraws returns requested dimensions", {
  draws <- ctsem:::ctOptimNormalDraws(c(1, 2), diag(c(.2, .3)), n=12)
  expect_equal(dim(draws), c(12, 2))
  expect_true(all(is.finite(draws)))
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

test_that("optimized uncertainty API uses explicit method and draw names", {
  fit_args <- names(formals(ctsem:::stanoptimis))
  update_args <- names(formals(ctOptimUncertainty))
  
  expect_true(all(c('uncertainty', 'uncertaintyDraws') %in% fit_args))
  expect_false('bootstrapUncertainty' %in% fit_args)
  expect_true(all(c('uncertainty', 'draws') %in% update_args))
  expect_false('sampleMethod' %in% update_args)
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
