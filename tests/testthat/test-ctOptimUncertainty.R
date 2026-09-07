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
  expect_true(diagnostics$usedNullProjection)
  expect_equal(diagnostics$nullDirections, 1L)
  expect_equal(diagnostics$nullParameters, 2L)
  expect_equal(diagnostics$ridge, 1e-6)
  # The second coordinate has no curvature, so it gets no spread -- not
  # `1 / ridge`, which is a property of the ridge and not of the data. The
  # first is untouched by the repair.
  expect_equal(unname(diag(cov)), c(1, 0))
})

test_that("Hessian covariance tries raw inversion before repair", {
  expect_silent(
    cov <- ctsem:::ctOptimCovFromHessian(-diag(2), context='clean Hessian'))
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'solve')
  expect_true(diagnostics$rawSolveSucceeded)
  expect_true(diagnostics$rawCholSucceeded)
  expect_false(diagnostics$usedNullProjection)
  expect_false(diagnostics$usedGinv)
  
  # An information matrix with a direction of negative curvature is not
  # inverted either. `solve()` returns a finite answer for it and only the
  # subsequent Cholesky notices, which is one rounding step away from not
  # noticing -- so the eigenvalues decide, before anything is inverted. The
  # answer is the same one nearPD reached by a longer route (unit variance on
  # the well-curved coordinate, none on the other); what is new is that the
  # dropped direction is named rather than smoothed away.
  hess <- diag(c(-1, 1))
  expect_warning(
    cov <- ctsem:::ctOptimCovFromHessian(hess, context='indefinite Hessian'),
    'no curvature')
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'nullprojection')
  expect_false(diagnostics$rawSolveSucceeded)
  expect_false(diagnostics$rawCholSucceeded)
  expect_false(diagnostics$usedNearPD)
  expect_true(diagnostics$usedNullProjection)
  expect_equal(diagnostics$nullParameters, 2L)
  expect_equal(unname(diag(cov)), c(1, 0))
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
  # maxStep is 200 rather than the default 64 so the cap cannot be what
  # produces the answer: the expansion grows by factors of 8 and lands on
  # exactly 64 of its own accord, which under maxStep=64 was indistinguishable
  # from running out of room. `expect_lt` is the assertion that fails if the
  # expansion ever does reach the cap.
  #
  # Measured step is exactly 64.0 -- identical over three repeats and at
  # maxStep 64, 200 and 1000 -- against a target of sqrt(2*2/.001) = 63.2456,
  # a relative difference of 0.0119. The tolerance below has 1.7x headroom
  # over that; the old tolerance of 1 had 84x and could discriminate nothing.
  prof <- ctsem:::ctOptimSurrogateProfileDirections(est=0, lpgFunc=lpg,
    cholcov=matrix(1), directions=matrix(1), targetDrop=2,
    maxStep=200, verbose=0)
  expect_true(all(prof$reached))
  expect_lt(max(prof$step), 200)
  expect_equal(prof$step, rep(sqrt(2 * 2 / .001), 2), tolerance=.02)
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

# --- a flat direction must not manufacture a standard error ------------------
#
# The defect these two cover, in its measured form. A benchmark fit reached the
# same optimum as its twin to eight decimal places -- same data, same starting
# values, same log likelihood -- and reported a drift interval a hundred times
# wider, with a point estimate that had moved with it. Its information matrix
# had three eigenvalues at 1e-16 against a largest of 8.8e5: a population SD
# with no individual differences behind it, and the correlation that goes with
# it.
#
# Flooring those eigenvalues at `ridge` and inverting put 1/ridge = 1e8 along
# each of them. The orientation of a null eigenvector is set by rounding error,
# because the block it spans is numerically zero, so it carries an arbitrary
# small component of the identified parameters and 1e8 multiplies that
# component. That is the leak, and it is not reproducible: the same Hessian
# perturbed by 1e-12 of its own scale gave standard errors of 0.019, 0.24 and
# 0.087 for the same well-determined parameter.

.leaky_information <- function(leak = 1e-3, curvature = c(1e4, 1e3)) {
  # Three parameters. The third has almost no curvature of its own, and the
  # direction it dominates carries a `leak` component of the first, which is
  # well determined. Built from its own eigendecomposition so the null
  # direction is exact rather than merely small.
  v3 <- c(leak, 0, 1); v3 <- v3 / sqrt(sum(v3^2))
  v1 <- c(1, 0, -leak); v1 <- v1 / sqrt(sum(v1^2))
  v2 <- c(0, 1, 0)
  V <- cbind(v1, v2, v3)
  V %*% diag(c(curvature, 0)) %*% t(V)
}

test_that("a direction with no curvature is projected out, not floored", {
  info <- .leaky_information()

  # What the old repair did, reconstructed here so the contrast is in the test
  # rather than only in the commit message: eigenvalues floored at the ridge,
  # then inverted.
  floored <- solve(ctsem:::ctOptimSafeCov(info, ridge = 1e-8))
  expect_gt(sqrt(diag(floored))[1], 5)      # measured 10, from 1e-3 * 1e4

  cov <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE)))
  diagnostics <- attr(cov, 'ctOptimCovFromHessian')
  expect_equal(diagnostics$method, 'nullprojection')
  expect_equal(diagnostics$nullDirections, 1L)
  expect_equal(diagnostics$nullParameters, 3L)

  # The identified parameter keeps the width its own curvature supports, and
  # the flat one gets none rather than a fabricated 1e4.
  expect_equal(sqrt(diag(cov))[1], 1 / sqrt(1e4), tolerance = 1e-6)
  expect_lt(sqrt(diag(cov))[3], 1e-3)

  # And it is stable. Rounding-scale noise in the information matrix used to
  # move the first parameter's standard error by an order of magnitude, because
  # the floored eigenvalue's 1e8 multiplied whatever component of it the
  # rotated null vector happened to pick up: measured on the real Hessian this
  # came from, 0.019 became 0.24 and then 0.087 under perturbations of 1e-12 of
  # its scale. The perturbation here is 1e-14 -- two orders below the tolerance
  # that decides which directions are dropped, so the classification is not
  # what is being tested; the arithmetic after it is.
  set.seed(11)
  scale <- max(abs(info))
  perturbed <- replicate(5, {
    E <- matrix(stats::rnorm(9), 3, 3); E <- (E + t(E)) / 2
    cv <- suppressWarnings(suppressMessages(
      ctsem:::ctOptimCovFromHessian(-(info + 1e-14 * scale * E), warn = FALSE)))
    sqrt(diag(cv))[1]
  })
  expect_equal(max(perturbed) / min(perturbed), 1, tolerance = 1e-5)
})

test_that("an interval wider than the curvature supports is detected and named", {
  info <- .leaky_information()
  parnames <- c('drift', 'diffusion', 'popsd')

  # The reported width under the old repair, against the width the curvature at
  # the estimate supports. This is the check a user with a single fit now has:
  # nothing else about such a fit looks wrong.
  leaked <- ctsem:::.ctBackendIntervalCheck(-info,
    sqrt(diag(solve(ctsem:::ctOptimSafeCov(info, ridge = 1e-8)))), parnames)
  expect_gt(leaked$nflagged, 0)
  expect_true('drift' %in% leaked$parameters)
  expect_gt(max(leaked$table$ratio), 100)

  # And it is quiet on the covariance that does not leak. Measured across the
  # healthy benchmark fits, every ratio sat between 1.0 and 2.5.
  cov <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE)))
  clean <- ctsem:::.ctBackendIntervalCheck(-info, sqrt(diag(cov)), parnames)
  expect_equal(clean$nflagged, 0L)
  expect_equal(clean$parameters, character())
  # A ratio is only defined where there is curvature to compare against; a
  # non-positive diagonal is `.ctBackendIdentifiability()`'s business.
  expect_equal(nrow(clean$table), 3L)
  expect_true(all(clean$table$ratio[clean$table$param != 'popsd'] < 2))
})

# The projection's own failure mode, which is the opposite of the one above and
# reads as a result rather than as a problem: a coordinate lying *along* a
# dropped direction inherits almost none of the variance that is left, so it is
# reported with a tight interval and a large z on a quantity the likelihood does
# not distinguish at all.
#
# The geometry is the measured one. On a one-latent model with `indvarying` on
# DRIFT and DIFFUSION under `intoverpop='augmented'` the flat direction mixes
# the population sd of the diffusion effect with its correlation to the drift
# effect -- the profile likelihood is bit-identical from r = 0.597 to r = 0.998
# while the sd compensates to hold their product fixed -- and `summary()`
# reported that correlation as 0.597 with sd 0.009 and z 65.3. Here the same
# direction is written down exactly rather than fitted, so the test is
# deterministic and needs no engine.
.mixed_null_information <- function(curvature = c(1e4, 1e3)) {
  # Two coordinates share the flat direction, 0.36 of one and 0.64 of the
  # other. Neither is flat on its own axis, which is the whole difficulty:
  # each has real curvature of its own and each is reported with a plausible
  # standard error.
  flat <- c(0.6, 0.8, 0)
  V <- cbind(c(-0.8, 0.6, 0), c(0, 0, 1), flat)
  V %*% diag(c(curvature, 0)) %*% t(V)
}

test_that("a parameter along a projected-out direction is reported as unidentified", {
  info <- .mixed_null_information()
  parnames <- c('popsd', 'rawcor', 'drift')
  cov <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info, warn = FALSE)))
  se <- sqrt(diag(cov))
  check <- ctsem:::.ctBackendIntervalCheck(-info, se, parnames)

  # The fault, stated as the test's premise: the reported spread is small and
  # entirely fictitious. 0.006 on a coordinate whose asymptotic variance is
  # infinite is a z of 100 at an estimate of 0.6.
  expect_equal(unname(se[2]), 0.006, tolerance = 1e-8)

  # And the check that now says so. The share of the dropped subspace is the
  # statement, and it is basis-invariant: any rotation within the null space
  # leaves these two numbers alone, which no single eigenvector's loading does.
  expect_equal(check$table$nullmass[match(parnames, check$table$param)],
    c(0.36, 0.64, 0), tolerance = 1e-8)
  expect_equal(check$nunidentified, 2L)
  expect_setequal(check$unidentified, c('popsd', 'rawcor'))

  # The existing width check cannot see it, which is why this is a second
  # branch rather than a lower threshold on the first: the ratio only grows.
  expect_equal(check$nflagged, 0L)
  # Below one, in fact -- the reported width is *narrower* than the curvature
  # in that one coordinate supports, which no genuine marginal standard error
  # can be.
  expect_true(all(check$table$ratio[check$table$param != 'drift'] < 1))
  # The identified coordinate keeps everything it had.
  expect_equal(check$table$ratio[check$table$param == 'drift'], 1,
    tolerance = 1e-8)

  # Undetermined coordinates come first, or a reader ordering by ratio meets
  # them last: their ratio is small precisely because their interval collapsed.
  expect_setequal(check$table$param[1:2], c('popsd', 'rawcor'))

  # `ctOptimCovFromHessian()` carries the same number, so the stan path -- which
  # has no interval check of its own -- can still say how many parameters lost
  # their spread rather than only how many directions were dropped.
  expect_equal(attr(cov, 'ctOptimCovFromHessian')$nullMass, c(0.36, 0.64, 0),
    tolerance = 1e-8)
})

test_that("a covariance with no flat direction flags nothing as unidentified", {
  # The margin matters as much as the verdict. A well conditioned information
  # matrix must give a null mass of exactly zero, not merely a small one, or
  # the threshold would be reading rounding.
  info <- diag(c(1e4, 1e3, 1e2))
  check <- ctsem:::.ctBackendIntervalCheck(-info, 1 / sqrt(diag(info)),
    c('a', 'b', 'c'))
  expect_equal(check$nunidentified, 0L)
  expect_equal(check$unidentified, character())
  expect_equal(check$table$nullmass, rep(0, 3))
})

test_that("imisScaleInit scales the whole proposal covariance, not its diagonal", {
  skip_if_not_installed('mvtnorm')
  skip_if_not_installed('diagis')
  skip_if_not_installed('gridExtra')
  skip_if_not_installed('ggplot2')

  # A strongly correlated Gaussian target, which is what makes the two
  # spellings of "scale the proposal" differ: the elementwise
  # `Sigma * (diag(s^2-1, n) + 1)` this replaced inflated the variances and
  # left the covariances alone, so it divided every proposal correlation by
  # s^2 -- narrower than Sigma along the correlated directions.
  d <- 5
  Sigma <- matrix(0.9, d, d); diag(Sigma) <- 1
  target <- function(x) mvtnorm::dmvnorm(x, rep(0, d), Sigma, log = TRUE)
  # max_iter = 0 is one batch drawn from the initial component alone, so the
  # effective sample size measures that component and nothing else.
  run <- function(S, s) {
    set.seed(7)
    ctsem:::imis_is(target, mu_hat = rep(0, d), Sigma_hat = S, max_iter = 0,
      scale_init = s, tail_scale = 1.2, df = Inf, target_ess = 1e9,
      n_batch = 4000, cl = NA, finishsamples = 100, verbose = FALSE,
      diag_plots = FALSE)
  }

  # `scale_init = s` and a proposal pre-scaled by s^2 are the same proposal.
  # They agree only to the ridge `safe_pd` adds, which is why this is a
  # tolerance rather than `expect_identical`.
  scaled <- run(Sigma, 1.5)
  prescaled <- run(Sigma * 1.5^2, 1)
  expect_equal(scaled$ess, prescaled$ess, tolerance = 1e-6)
  expect_equal(scaled$covariance, prescaled$covariance, tolerance = 1e-6)

  # And it matters: for the same nominal scale the elementwise form is a
  # markedly worse proposal on a correlated target.
  elementwise <- run(Sigma * (diag(1.5^2 - 1, d) + 1), 1)
  expect_gt(scaled$ess, 2 * elementwise$ess)

  # `scale_init = 1` still leaves the proposal alone, which is what
  # `ctParticleCorrect()` and `ctLaplaceCorrect()` rely on when they pre-scale
  # their own proposal and pass 1: the weighted covariance recovers the target.
  unscaled <- run(Sigma, 1)
  expect_equal(unscaled$covariance, Sigma, tolerance = 0.1)
})
